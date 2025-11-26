#ifndef DMC_ROUTE_ROUTER_HPP
#define DMC_ROUTE_ROUTER_HPP

#include "types.hpp"
#include "network.hpp"
#include <vector>
#include <functional>

namespace dmc {

/**
 * Configuration for the router.
 */
struct RouterConfig {
    double dt = 3600.0;                    // Timestep [s]
    bool enable_gradients = true;          // Record tape for AD
    double min_flow = 1e-6;                // Minimum flow [m³/s]
    double x_lower_bound = 0.0;            // Muskingum X bounds
    double x_upper_bound = 0.5;
    bool adaptive_substepping = false;     // Substep if Courant > 1
    int max_substeps = 10;
};

/**
 * Core Muskingum-Cunge routing engine.
 * 
 * Supports automatic differentiation for parameter learning.
 */
class MuskingumCungeRouter {
public:
    explicit MuskingumCungeRouter(Network& network, RouterConfig config = {});
    
    // =========== Core Routing ===========
    
    /**
     * Route a single reach for one timestep.
     * Returns the outflow Q(t+Δt).
     */
    Real route_reach(Reach& reach, double dt);
    
    /**
     * Route entire network for one timestep.
     * Processes reaches in topological order.
     */
    void route_timestep();
    
    /**
     * Route for multiple timesteps.
     */
    void route(int num_timesteps);
    
    // =========== Gradient Computation ===========
    
    /**
     * Enable/disable gradient recording.
     */
    void enable_gradients(bool enable);
    
    /**
     * Start recording operations to tape.
     * Call before simulation loop.
     */
    void start_recording();
    
    /**
     * Stop recording and prepare for backward pass.
     * Call after simulation loop.
     */
    void stop_recording();
    
    /**
     * Compute gradients via reverse AD.
     * @param gauge_reaches Reach IDs where we have observations
     * @param dL_dQ Gradient of loss w.r.t. discharge at each gauge
     */
    void compute_gradients(const std::vector<int>& gauge_reaches,
                           const std::vector<double>& dL_dQ);
    
    /**
     * Get accumulated gradients for all parameters.
     */
    std::unordered_map<std::string, double> get_gradients() const;
    
    /**
     * Reset tape and gradients.
     */
    void reset_gradients();
    
    // =========== State Management ===========
    
    /**
     * Set lateral inflow for a reach (from rainfall-runoff model).
     */
    void set_lateral_inflow(int reach_id, double inflow);
    
    /**
     * Set lateral inflow for all reaches (same order as topological order).
     */
    void set_lateral_inflows(const std::vector<double>& inflows);
    
    /**
     * Get current discharge at a reach.
     */
    double get_discharge(int reach_id) const;
    
    /**
     * Get discharge at all reaches (same order as topological order).
     */
    std::vector<double> get_all_discharges() const;
    
    /**
     * Reset state to initial conditions.
     */
    void reset_state();
    
    // =========== Time Management ===========
    
    double current_time() const { return current_time_; }
    void set_time(double t) { current_time_ = t; }
    
    // =========== Access ===========
    
    Network& network() { return network_; }
    const Network& network() const { return network_; }
    const RouterConfig& config() const { return config_; }
    
private:
    Network& network_;
    RouterConfig config_;
    double current_time_ = 0.0;
    bool recording_ = false;
    
    // Gauge output storage for gradient computation
    std::vector<int> gauge_reach_ids_;
    std::vector<Real> gauge_outputs_;
    
    /**
     * Compute inflow to a reach from upstream junction.
     */
    Real compute_reach_inflow(const Reach& reach);
    
    /**
     * Compute Muskingum K parameter.
     */
    Real compute_K(const Reach& reach, const Real& Q_ref);
    
    /**
     * Compute Muskingum X parameter.
     */
    Real compute_X(const Reach& reach, const Real& Q_ref, const Real& K);
    
    /**
     * Compute routing coefficients C1, C2, C3, C4.
     */
    void compute_routing_coefficients(const Real& K, const Real& X, double dt,
                                       Real& C1, Real& C2, Real& C3, Real& C4);
};

// ==================== Implementation ====================

inline MuskingumCungeRouter::MuskingumCungeRouter(Network& network, RouterConfig config)
    : network_(network), config_(std::move(config)) {
    network_.build_topology();
}

inline Real MuskingumCungeRouter::compute_K(const Reach& reach, const Real& Q_ref) {
    // Hydraulic radius
    Real R_h = reach.geometry.hydraulic_radius(Q_ref);
    
    // Velocity from Manning's equation: v = (1/n) * R^(2/3) * S^(1/2)
    Real velocity = (1.0 / reach.manning_n) * 
                    safe_pow(R_h, 2.0/3.0) * 
                    safe_sqrt(Real(reach.slope));
    
    // Wave celerity: c = (5/3) * v for wide rectangular channel
    Real celerity = (5.0 / 3.0) * velocity;
    celerity = safe_max(celerity, Real(0.1));  // Prevent division by zero
    
    // K = Δx / c
    return Real(reach.length) / celerity;
}

inline Real MuskingumCungeRouter::compute_X(const Reach& reach, const Real& Q_ref, 
                                             const Real& K) {
    Real width = reach.geometry.width(Q_ref);
    Real celerity = Real(reach.length) / K;
    
    // X = 0.5 - Q / (2 * c * B * S₀ * Δx)
    Real X = 0.5 - Q_ref / (2.0 * celerity * width * 
                            Real(reach.slope) * Real(reach.length));
    
    // Clamp to valid range [0, 0.5]
    return clamp(X, Real(config_.x_lower_bound), Real(config_.x_upper_bound));
}

inline void MuskingumCungeRouter::compute_routing_coefficients(
    const Real& K, const Real& X, double dt,
    Real& C1, Real& C2, Real& C3, Real& C4) {
    
    Real denom = 2.0 * K * (1.0 - X) + dt;
    
    C1 = (dt - 2.0 * K * X) / denom;
    C2 = (dt + 2.0 * K * X) / denom;
    C3 = (2.0 * K * (1.0 - X) - dt) / denom;
    C4 = 2.0 * dt / denom;
}

inline Real MuskingumCungeRouter::compute_reach_inflow(const Reach& reach) {
    if (reach.upstream_junction_id < 0) {
        // Headwater - no upstream inflow
        return Real(0.0);
    }
    
    const Junction& junc = network_.get_junction(reach.upstream_junction_id);
    Real total_inflow = junc.external_inflow;
    
    for (int up_reach_id : junc.upstream_reach_ids) {
        const Reach& up_reach = network_.get_reach(up_reach_id);
        total_inflow = total_inflow + up_reach.outflow_curr;
    }
    
    return total_inflow;
}

inline Real MuskingumCungeRouter::route_reach(Reach& reach, double dt) {
    // Reference discharge for hydraulic calculations
    Real Q_ref = safe_max(reach.outflow_prev, Real(config_.min_flow));
    
    // Compute Muskingum parameters
    Real K = compute_K(reach, Q_ref);
    Real X = compute_X(reach, Q_ref, K);
    
    // Store for diagnostics
    reach.K = K;
    reach.X = X;
    reach.celerity = Real(reach.length) / K;
    reach.velocity = reach.celerity * (3.0 / 5.0);
    
    // Compute routing coefficients
    Real C1, C2, C3, C4;
    compute_routing_coefficients(K, X, dt, C1, C2, C3, C4);
    
    // Apply Muskingum-Cunge equation
    Real outflow = C1 * reach.inflow_curr + 
                   C2 * reach.inflow_prev + 
                   C3 * reach.outflow_prev +
                   C4 * reach.lateral_inflow;
    
    // Ensure non-negative
    outflow = safe_max(outflow, Real(config_.min_flow));
    
    return outflow;
}

inline void MuskingumCungeRouter::route_timestep() {
    // Process reaches in topological order
    for (int reach_id : network_.topological_order()) {
        Reach& reach = network_.get_reach(reach_id);
        
        // Update inflow from upstream
        reach.inflow_curr = compute_reach_inflow(reach);
        
        // Route through reach
        reach.outflow_curr = route_reach(reach, config_.dt);
    }
    
    // Advance state: current becomes previous
    for (int reach_id : network_.topological_order()) {
        Reach& reach = network_.get_reach(reach_id);
        reach.inflow_prev = reach.inflow_curr;
        reach.outflow_prev = reach.outflow_curr;
    }
    
    current_time_ += config_.dt;
}

inline void MuskingumCungeRouter::route(int num_timesteps) {
    for (int t = 0; t < num_timesteps; ++t) {
        route_timestep();
    }
}

inline void MuskingumCungeRouter::enable_gradients(bool enable) {
    config_.enable_gradients = enable;
}

inline void MuskingumCungeRouter::start_recording() {
    if (!config_.enable_gradients || !AD_ENABLED) return;
    
    activate_tape();
    recording_ = true;
    
    // Register all parameters as inputs
    for (Real* param : network_.get_all_parameters()) {
        register_input(*param);
    }
}

inline void MuskingumCungeRouter::stop_recording() {
    if (!recording_) return;
    
    // Don't deactivate tape yet - we need it active for gradient computation
    // Just mark that the forward pass is done
    recording_ = false;
}

inline void MuskingumCungeRouter::compute_gradients(
    const std::vector<int>& gauge_reaches,
    const std::vector<double>& dL_dQ) {
    
    if (!AD_ENABLED) return;
    
    // Register outputs and seed with loss gradients
    for (size_t i = 0; i < gauge_reaches.size(); ++i) {
        Reach& reach = network_.get_reach(gauge_reaches[i]);
        register_output(reach.outflow_curr);
        set_gradient(reach.outflow_curr, dL_dQ[i]);
    }
    
    // Reverse pass
    evaluate_tape();
    
    // Collect gradients from tape
    network_.collect_gradients();
    
    // NOW deactivate the tape
    deactivate_tape();
}

inline std::unordered_map<std::string, double> MuskingumCungeRouter::get_gradients() const {
    std::unordered_map<std::string, double> grads;
    
    for (int reach_id : network_.topological_order()) {
        const Reach& reach = network_.get_reach(reach_id);
        std::string prefix = "reach_" + std::to_string(reach_id) + "_";
        
        grads[prefix + "manning_n"] = reach.grad_manning_n;
        grads[prefix + "width_coef"] = reach.grad_width_coef;
        grads[prefix + "width_exp"] = reach.grad_width_exp;
        grads[prefix + "depth_coef"] = reach.grad_depth_coef;
        grads[prefix + "depth_exp"] = reach.grad_depth_exp;
    }
    
    return grads;
}

inline void MuskingumCungeRouter::reset_gradients() {
    reset_tape();
    network_.zero_gradients();
}

inline void MuskingumCungeRouter::set_lateral_inflow(int reach_id, double inflow) {
    network_.get_reach(reach_id).lateral_inflow = Real(inflow);
}

inline void MuskingumCungeRouter::set_lateral_inflows(const std::vector<double>& inflows) {
    const auto& order = network_.topological_order();
    for (size_t i = 0; i < order.size() && i < inflows.size(); ++i) {
        network_.get_reach(order[i]).lateral_inflow = Real(inflows[i]);
    }
}

inline double MuskingumCungeRouter::get_discharge(int reach_id) const {
    return to_double(network_.get_reach(reach_id).outflow_curr);
}

inline std::vector<double> MuskingumCungeRouter::get_all_discharges() const {
    std::vector<double> discharges;
    for (int reach_id : network_.topological_order()) {
        discharges.push_back(to_double(network_.get_reach(reach_id).outflow_curr));
    }
    return discharges;
}

inline void MuskingumCungeRouter::reset_state() {
    for (int reach_id : network_.topological_order()) {
        Reach& reach = network_.get_reach(reach_id);
        reach.inflow_prev = Real(0.0);
        reach.inflow_curr = Real(0.0);
        reach.outflow_prev = Real(0.0);
        reach.outflow_curr = Real(0.0);
        reach.lateral_inflow = Real(0.0);
    }
    current_time_ = 0.0;
    reset_gradients();
}

} // namespace dmc

#endif // DMC_ROUTE_ROUTER_HPP
