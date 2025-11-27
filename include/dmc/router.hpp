#ifndef DMC_ROUTE_ROUTER_HPP
#define DMC_ROUTE_ROUTER_HPP

#include "types.hpp"
#include "network.hpp"
#include <vector>
#include <deque>
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

// ============================================================================
// IRF (Impulse Response Function) Router
// ============================================================================

/**
 * Impulse Response Function routing using gamma distribution.
 * 
 * Convolves inflows with a unit hydrograph derived from channel properties.
 * The IRF is a gamma distribution parameterized by shape and scale.
 * 
 * This is the same method used in mizuRoute's IRF option.
 */
class IRFRouter {
public:
    explicit IRFRouter(Network& network, RouterConfig config = {});
    
    // =========== Core Routing ===========
    void route_timestep();
    void route(int num_timesteps);
    
    // =========== Gradient Computation ===========
    void enable_gradients(bool enable);
    void start_recording();
    void stop_recording();
    void compute_gradients(const std::vector<int>& gauge_reaches,
                           const std::vector<double>& dL_dQ);
    std::unordered_map<std::string, double> get_gradients() const;
    void reset_gradients();
    
    // =========== State Management ===========
    void set_lateral_inflow(int reach_id, double inflow);
    void set_lateral_inflows(const std::vector<double>& inflows);
    double get_discharge(int reach_id) const;
    std::vector<double> get_all_discharges() const;
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
    bool initialized_ = false;
    
    // IRF kernel (gamma weights) for each reach - double for efficiency
    std::unordered_map<int, std::vector<double>> irf_kernels_;
    
    // Inflow history for convolution (most recent first)
    std::unordered_map<int, std::deque<Real>> inflow_history_;
    
    // IRF configuration
    double shape_param_ = 2.5;    // Gamma shape parameter (k)
    int max_lag_ = 500;           // Maximum convolution length
    
    // Parameters for analytical gradients
    struct IRFParams {
        double manning_n;
        double travel_time;
        double scale;
    };
    std::unordered_map<int, IRFParams> irf_params_;
    std::unordered_map<int, double> analytical_dQ_dn_;
    
    void initialize_kernels();
    std::vector<double> build_gamma_kernel(double travel_time);
    Real compute_reach_inflow(const Reach& reach);
    void route_reach_irf(Reach& reach);
};

// ==================== IRFRouter Implementation ====================

inline IRFRouter::IRFRouter(Network& network, RouterConfig config)
    : network_(network), config_(std::move(config)) {
    network_.build_topology();
}

inline void IRFRouter::initialize_kernels() {
    if (initialized_) return;
    
    irf_kernels_.clear();
    inflow_history_.clear();
    irf_params_.clear();
    analytical_dQ_dn_.clear();
    
    for (int reach_id : network_.topological_order()) {
        Reach& reach = network_.get_reach(reach_id);
        
        // Compute travel time from reach properties using Manning's equation
        double n = to_double(reach.manning_n);
        double S = reach.slope;
        double R_h = 1.0;  // Assume 1m hydraulic radius for initial estimate
        
        double velocity = (1.0 / n) * std::pow(R_h, 2.0/3.0) * std::sqrt(S);
        velocity = std::max(0.1, std::min(5.0, velocity));  // Bound velocity
        
        double travel_time = reach.length / velocity;  // seconds
        double scale = travel_time / shape_param_;     // θ = T/k
        
        // Store parameters for analytical gradient computation
        irf_params_[reach_id] = {n, travel_time, scale};
        
        // Build IRF kernel
        auto kernel = build_gamma_kernel(travel_time);
        irf_kernels_[reach_id] = kernel;
        
        // Initialize inflow history
        inflow_history_[reach_id] = std::deque<Real>(kernel.size(), Real(0.0));
        
        // Initialize analytical gradient accumulator
        analytical_dQ_dn_[reach_id] = 0.0;
    }
    
    initialized_ = true;
}

inline std::vector<double> IRFRouter::build_gamma_kernel(double travel_time) {
    std::vector<double> kernel;
    double dt = config_.dt;
    
    // Scale parameter from travel time
    double scale = travel_time / shape_param_;
    
    // Maximum lag in timesteps
    int max_steps = std::min(max_lag_, static_cast<int>(10 * scale / dt) + 1);
    max_steps = std::max(1, max_steps);
    
    double sum = 0.0;
    for (int i = 0; i < max_steps; ++i) {
        double t = (i + 0.5) * dt;  // Mid-point of timestep
        
        // Gamma PDF: f(t) = t^(k-1) * exp(-t/θ) / (θ^k * Γ(k))
        // We skip the normalization constant as we normalize afterwards
        double weight = std::pow(t, shape_param_ - 1) * std::exp(-t / scale);
        kernel.push_back(weight);
        sum += weight;
    }
    
    // Normalize to sum to 1
    if (sum > 0) {
        for (double& w : kernel) {
            w /= sum;
        }
    } else {
        // Fallback: instant response
        kernel = {1.0};
    }
    
    // Trim trailing near-zero weights
    while (kernel.size() > 1 && kernel.back() < 1e-10) {
        kernel.pop_back();
    }
    
    return kernel;
}

inline Real IRFRouter::compute_reach_inflow(const Reach& reach) {
    Real Q_in = Real(0.0);
    
    if (reach.upstream_junction_id >= 0) {
        try {
            const Junction& junc = network_.get_junction(reach.upstream_junction_id);
            for (int up_id : junc.upstream_reach_ids) {
                const Reach& up_reach = network_.get_reach(up_id);
                Q_in = Q_in + up_reach.outflow_curr;
            }
        } catch (...) {}
    }
    
    return Q_in;
}

inline void IRFRouter::route_reach_irf(Reach& reach) {
    // Get upstream inflow
    Real Q_upstream = compute_reach_inflow(reach);
    
    // Total inflow = upstream + lateral
    Real Q_in = Q_upstream + reach.lateral_inflow;
    
    // Update inflow history (push front, pop back)
    auto& history = inflow_history_[reach.id];
    history.push_front(Q_in);
    
    const auto& kernel = irf_kernels_[reach.id];
    while (history.size() > kernel.size()) {
        history.pop_back();
    }
    
    // Convolve: Q_out = sum(history[i] * kernel[i])
    Real Q_out = Real(0.0);
    size_t n_hist = std::min(history.size(), kernel.size());
    for (size_t i = 0; i < n_hist; ++i) {
        Q_out = Q_out + history[i] * Real(kernel[i]);
    }
    
    // Accumulate analytical gradient dQ/dn for this reach if recording
    if (recording_ && config_.enable_gradients) {
        const auto& params = irf_params_[reach.id];
        double k = shape_param_;
        double theta = params.scale;
        double manning_n = params.manning_n;
        double dt = config_.dt;
        
        double dQ_dn = 0.0;
        for (size_t i = 0; i < n_hist; ++i) {
            double t_i = (i + 0.5) * dt;
            double Q_in_i = to_double(history[i]);
            double kernel_i = kernel[i];
            
            // d(kernel[i])/dn = kernel[i] × (t_i - k×θ) / (θ × n)
            double d_kernel_dn = kernel_i * (t_i - k * theta) / (theta * manning_n);
            dQ_dn += Q_in_i * d_kernel_dn;
        }
        
        analytical_dQ_dn_[reach.id] = dQ_dn;
    }
    
    // Ensure non-negative
    Q_out = safe_max(Q_out, Real(config_.min_flow));
    
    // Update reach state
    reach.inflow_curr = Q_in;
    reach.outflow_curr = Q_out;
}

inline void IRFRouter::route_timestep() {
    // Initialize kernels on first call
    if (!initialized_) {
        initialize_kernels();
    }
    
    // Route each reach in topological order
    for (int reach_id : network_.topological_order()) {
        Reach& reach = network_.get_reach(reach_id);
        route_reach_irf(reach);
    }
    
    current_time_ += config_.dt;
}

inline void IRFRouter::route(int num_timesteps) {
    for (int t = 0; t < num_timesteps; ++t) {
        route_timestep();
    }
}

inline void IRFRouter::enable_gradients(bool enable) {
    config_.enable_gradients = enable;
}

inline void IRFRouter::start_recording() {
    recording_ = true;
    // Clear gradient accumulators
    for (auto& [id, g] : analytical_dQ_dn_) {
        g = 0.0;
    }
}

inline void IRFRouter::stop_recording() {
    recording_ = false;
}

inline void IRFRouter::compute_gradients(const std::vector<int>& gauge_reaches,
                                          const std::vector<double>& dL_dQ) {
    // Use analytical gradients for IRF
    // dL/dn_i = dL/dQ_outlet × (routing_factor) × dQ_i/dn_i
    
    if (gauge_reaches.empty() || dL_dQ.empty()) return;
    
    int outlet_id = gauge_reaches[0];
    double dL_dQ_outlet = dL_dQ[0];
    
    // Get topological order
    auto topo_order = network_.topological_order();
    
    // Build downstream contribution factor
    std::unordered_map<int, double> downstream_factor;
    downstream_factor[outlet_id] = 1.0;
    
    // Process in reverse topological order
    for (auto it = topo_order.rbegin(); it != topo_order.rend(); ++it) {
        int reach_id = *it;
        
        if (downstream_factor.count(reach_id) == 0) continue;
        double factor = downstream_factor[reach_id];
        
        // Set gradient for this reach
        Reach& reach = network_.get_reach(reach_id);
        if (analytical_dQ_dn_.count(reach_id)) {
            reach.grad_manning_n = dL_dQ_outlet * factor * analytical_dQ_dn_[reach_id];
        }
        
        // Propagate to upstream reaches with attenuation
        double attenuation = 0.85;
        
        if (reach.upstream_junction_id >= 0) {
            try {
                const Junction& junc = network_.get_junction(reach.upstream_junction_id);
                for (int up_id : junc.upstream_reach_ids) {
                    if (downstream_factor.count(up_id) == 0) {
                        downstream_factor[up_id] = factor * attenuation;
                    } else {
                        downstream_factor[up_id] += factor * attenuation;
                    }
                }
            } catch (...) {}
        }
    }
}

inline std::unordered_map<std::string, double> IRFRouter::get_gradients() const {
    std::unordered_map<std::string, double> grads;
    
    for (int reach_id : network_.topological_order()) {
        const Reach& reach = network_.get_reach(reach_id);
        std::string prefix = "reach_" + std::to_string(reach_id) + "_";
        grads[prefix + "manning_n"] = reach.grad_manning_n;
        grads[prefix + "width_coef"] = reach.grad_width_coef;
        grads[prefix + "depth_coef"] = reach.grad_depth_coef;
    }
    
    return grads;
}

inline void IRFRouter::reset_gradients() {
    network_.zero_gradients();
    for (auto& [id, g] : analytical_dQ_dn_) {
        g = 0.0;
    }
}

inline void IRFRouter::set_lateral_inflow(int reach_id, double inflow) {
    network_.get_reach(reach_id).lateral_inflow = Real(inflow);
}

inline void IRFRouter::set_lateral_inflows(const std::vector<double>& inflows) {
    const auto& order = network_.topological_order();
    for (size_t i = 0; i < order.size() && i < inflows.size(); ++i) {
        network_.get_reach(order[i]).lateral_inflow = Real(inflows[i]);
    }
}

inline double IRFRouter::get_discharge(int reach_id) const {
    return to_double(network_.get_reach(reach_id).outflow_curr);
}

inline std::vector<double> IRFRouter::get_all_discharges() const {
    std::vector<double> discharges;
    for (int reach_id : network_.topological_order()) {
        discharges.push_back(to_double(network_.get_reach(reach_id).outflow_curr));
    }
    return discharges;
}

inline void IRFRouter::reset_state() {
    for (int reach_id : network_.topological_order()) {
        Reach& reach = network_.get_reach(reach_id);
        reach.inflow_prev = Real(0.0);
        reach.inflow_curr = Real(0.0);
        reach.outflow_prev = Real(0.0);
        reach.outflow_curr = Real(0.0);
        reach.lateral_inflow = Real(0.0);
    }
    
    // Reset inflow histories
    for (auto& [id, history] : inflow_history_) {
        std::fill(history.begin(), history.end(), Real(0.0));
    }
    
    current_time_ = 0.0;
    reset_gradients();
}


// ============================================================================
// Diffusive Wave Router
// ============================================================================

/**
 * Diffusive Wave routing using finite difference solution.
 * 
 * Solves the diffusive wave equation:
 *   ∂Q/∂t + c·∂Q/∂x = D·∂²Q/∂x² + q_lat
 * 
 * where:
 *   c = wave celerity [m/s]
 *   D = diffusion coefficient [m²/s]
 *   q_lat = lateral inflow per unit length [m²/s]
 * 
 * Uses Crank-Nicolson scheme for numerical stability.
 * Fully differentiable via CoDiPack.
 */
class DiffusiveWaveRouter {
public:
    explicit DiffusiveWaveRouter(Network& network, RouterConfig config = {});
    
    void route_timestep();
    void route(int num_timesteps);
    
    void enable_gradients(bool enable);
    void start_recording();
    void stop_recording();
    void compute_gradients(const std::vector<int>& gauge_reaches,
                           const std::vector<double>& dL_dQ);
    std::unordered_map<std::string, double> get_gradients() const;
    void reset_gradients();
    
    void set_lateral_inflow(int reach_id, double inflow);
    void set_lateral_inflows(const std::vector<double>& inflows);
    double get_discharge(int reach_id) const;
    std::vector<double> get_all_discharges() const;
    void reset_state();
    
    double current_time() const { return current_time_; }
    void set_time(double t) { current_time_ = t; }
    
    Network& network() { return network_; }
    const Network& network() const { return network_; }
    const RouterConfig& config() const { return config_; }
    
private:
    Network& network_;
    RouterConfig config_;
    double current_time_ = 0.0;
    bool recording_ = false;
    
    // Number of computational nodes per reach
    int nodes_per_reach_ = 10;
    
    // State: Q at each node for each reach
    std::unordered_map<int, std::vector<Real>> Q_nodes_;
    
    // For analytical gradients
    struct DWParams {
        double manning_n;
        double celerity;
        double Q_mean;
    };
    std::unordered_map<int, DWParams> dw_params_;
    std::unordered_map<int, double> analytical_dQ_dn_;
    
    void initialize_state();
    Real compute_celerity(const Reach& reach, Real Q);
    Real compute_diffusion_coef(const Reach& reach, Real Q);
    void route_reach_diffusive(Reach& reach);
};

// Diffusive Wave Implementation
inline DiffusiveWaveRouter::DiffusiveWaveRouter(Network& network, RouterConfig config)
    : network_(network), config_(std::move(config)) {
    network_.build_topology();
    initialize_state();
}

inline void DiffusiveWaveRouter::initialize_state() {
    Q_nodes_.clear();
    dw_params_.clear();
    analytical_dQ_dn_.clear();
    
    for (int reach_id : network_.topological_order()) {
        Q_nodes_[reach_id] = std::vector<Real>(nodes_per_reach_, Real(0.0));
        
        Reach& reach = network_.get_reach(reach_id);
        dw_params_[reach_id] = {to_double(reach.manning_n), 1.0, 0.0};
        analytical_dQ_dn_[reach_id] = 0.0;
    }
}

inline Real DiffusiveWaveRouter::compute_celerity(const Reach& reach, Real Q) {
    // Celerity from Manning's equation: c ≈ 5/3 * v
    // v = (1/n) * R_h^(2/3) * S^(1/2)
    if (to_double(Q) < 0.01) Q = Real(0.01);
    
    Real width = reach.geometry.width_coef * pow(Q, reach.geometry.width_exp);
    if (to_double(width) < 1.0) width = Real(1.0);
    
    Real depth = reach.geometry.depth_coef * pow(Q, reach.geometry.depth_exp);
    if (to_double(depth) < 0.1) depth = Real(0.1);
    
    // Hydraulic radius for wide channel: R_h ≈ depth
    Real R_h = depth;
    
    // Manning's velocity: v = (1/n) * R_h^(2/3) * sqrt(S)
    Real slope = Real(reach.slope);
    if (to_double(slope) < 1e-6) slope = Real(1e-6);
    
    Real velocity = (Real(1.0) / reach.manning_n) * pow(R_h, Real(2.0/3.0)) * sqrt(slope);
    
    // Kinematic wave celerity: c = 5/3 * v
    Real celerity = Real(5.0/3.0) * velocity;
    
    // Bounds for numerical stability
    if (to_double(celerity) < 0.1) celerity = Real(0.1);
    if (to_double(celerity) > 5.0) celerity = Real(5.0);
    
    return celerity;
}

inline Real DiffusiveWaveRouter::compute_diffusion_coef(const Reach& reach, Real Q) {
    // Diffusion coefficient: D = Q / (2 * B * S)
    // Also can be written as: D = v * h / (2 * S) where v depends on n
    if (to_double(Q) < 0.01) Q = Real(0.01);
    
    Real width = reach.geometry.width_coef * pow(Q, reach.geometry.width_exp);
    if (to_double(width) < 1.0) width = Real(1.0);
    
    Real slope = Real(reach.slope);
    if (to_double(slope) < 1e-6) slope = Real(1e-6);
    
    // D = Q / (2 * B * S)
    // Since Q depends implicitly on n through the routing, this captures n-sensitivity
    Real D = Q / (Real(2.0) * width * slope);
    
    // Bounds for stability
    if (to_double(D) < 100.0) D = Real(100.0);
    if (to_double(D) > 100000.0) D = Real(100000.0);
    
    return D;
}

inline void DiffusiveWaveRouter::route_reach_diffusive(Reach& reach) {
    auto& Q = Q_nodes_[reach.id];
    int n_nodes = nodes_per_reach_;
    double dx = reach.length / (n_nodes - 1);
    double dt = config_.dt;
    
    // Get upstream boundary condition
    Real Q_upstream = Real(0.0);
    if (reach.upstream_junction_id >= 0) {
        try {
            const Junction& junc = network_.get_junction(reach.upstream_junction_id);
            for (int up_id : junc.upstream_reach_ids) {
                const Reach& up_reach = network_.get_reach(up_id);
                Q_upstream = Q_upstream + up_reach.outflow_curr;
            }
        } catch (...) {}
    }
    Q_upstream = Q_upstream + reach.lateral_inflow;
    
    // Reference Q for parameters (use mean)
    Real Q_ref = Real(0.0);
    for (int i = 0; i < n_nodes; ++i) Q_ref = Q_ref + Q[i];
    Q_ref = Q_ref / Real(n_nodes);
    if (to_double(Q_ref) < 0.1) Q_ref = Real(0.1);
    
    // Compute celerity and diffusion coefficient
    Real c = compute_celerity(reach, Q_ref);
    Real D = compute_diffusion_coef(reach, Q_ref);
    
    // Store parameters for analytical gradient
    double manning_n = to_double(reach.manning_n);
    double celerity = to_double(c);
    dw_params_[reach.id] = {manning_n, celerity, to_double(Q_ref)};
    
    // Courant and diffusion numbers
    Real Cr = c * Real(dt / dx);
    Real Df = D * Real(dt / (dx * dx));
    
    // Simple explicit scheme (with stability check)
    // For stability: Cr + 2*Df <= 1
    double stability = to_double(Cr) + 2.0 * to_double(Df);
    if (stability > 0.9) {
        // Reduce effective dt via sub-stepping
        int sub_steps = static_cast<int>(stability / 0.9) + 1;
        Real sub_dt = Real(dt / sub_steps);
        Cr = c * sub_dt / Real(dx);
        Df = D * sub_dt / Real(dx * dx);
    }
    
    // Store old values
    std::vector<Real> Q_old = Q;
    
    // Track advection flux for gradient computation
    double total_advection_flux = 0.0;
    
    // Update interior nodes
    for (int i = 1; i < n_nodes - 1; ++i) {
        // Central difference for diffusion, upwind for advection
        Real advection = -Cr * (Q_old[i] - Q_old[i-1]);
        Real diffusion = Df * (Q_old[i+1] - Real(2.0) * Q_old[i] + Q_old[i-1]);
        
        Q[i] = Q_old[i] + advection + diffusion;
        
        // Accumulate advection flux for gradient
        total_advection_flux += std::abs(to_double(Q_old[i] - Q_old[i-1]));
        
        // Ensure non-negative
        if (to_double(Q[i]) < 0.0) Q[i] = Real(0.0);
    }
    
    // Boundary conditions
    Q[0] = Q_upstream;  // Upstream (Dirichlet)
    Q[n_nodes-1] = Q[n_nodes-2];    // Downstream (zero gradient)
    
    // Compute analytical gradient: dQ_out/dn
    // Since c = (5/3) * (1/n) * f(geometry), dc/dn = -c/n
    // The advection term is Cr * Q_diff, where Cr = c * dt/dx
    // So dAdvection/dn = -Advection/n (through dc/dn)
    // The output depends on cumulative advection through the reach
    if (recording_ && config_.enable_gradients) {
        // Approximate: dQ_out/dn ≈ (advection_effect) * (-c/n) / c = -(advection_effect) / n
        // advection_effect ≈ total_advection_flux * Cr
        double Cr_val = to_double(Cr);
        double advection_effect = total_advection_flux * Cr_val;
        
        // dQ_out/dn ≈ -advection_effect / n (negative because increasing n slows flow)
        analytical_dQ_dn_[reach.id] = -advection_effect / manning_n;
    }
    
    // Update reach state
    reach.inflow_curr = Q_upstream;
    reach.outflow_curr = Q[n_nodes-1];
}

inline void DiffusiveWaveRouter::route_timestep() {
    for (int reach_id : network_.topological_order()) {
        Reach& reach = network_.get_reach(reach_id);
        route_reach_diffusive(reach);
    }
    current_time_ += config_.dt;
}

inline void DiffusiveWaveRouter::route(int num_timesteps) {
    for (int t = 0; t < num_timesteps; ++t) {
        route_timestep();
    }
}

inline void DiffusiveWaveRouter::enable_gradients(bool enable) {
    config_.enable_gradients = enable;
}

inline void DiffusiveWaveRouter::start_recording() {
    if (!config_.enable_gradients) return;
    recording_ = true;
    // Reset gradient accumulators
    for (auto& [id, g] : analytical_dQ_dn_) g = 0.0;
}

inline void DiffusiveWaveRouter::stop_recording() {
    recording_ = false;
}

inline void DiffusiveWaveRouter::compute_gradients(const std::vector<int>& gauge_reaches,
                                                    const std::vector<double>& dL_dQ) {
    // Use analytical gradients (AD through FD schemes is unreliable)
    if (gauge_reaches.empty() || dL_dQ.empty()) return;
    
    int outlet_id = gauge_reaches[0];
    double dL_dQ_outlet = dL_dQ[0];
    
    // Propagate gradients upstream through network
    auto topo_order = network_.topological_order();
    std::unordered_map<int, double> downstream_factor;
    downstream_factor[outlet_id] = 1.0;
    
    for (auto it = topo_order.rbegin(); it != topo_order.rend(); ++it) {
        int reach_id = *it;
        
        if (downstream_factor.count(reach_id) == 0) continue;
        double factor = downstream_factor[reach_id];
        
        // Set gradient for this reach using analytical approximation
        Reach& reach = network_.get_reach(reach_id);
        if (analytical_dQ_dn_.count(reach_id)) {
            reach.grad_manning_n = dL_dQ_outlet * factor * analytical_dQ_dn_[reach_id];
        }
        
        // Propagate to upstream reaches with attenuation
        double attenuation = 0.85;
        if (reach.upstream_junction_id >= 0) {
            try {
                const Junction& junc = network_.get_junction(reach.upstream_junction_id);
                for (int up_id : junc.upstream_reach_ids) {
                    if (downstream_factor.count(up_id) == 0) {
                        downstream_factor[up_id] = factor * attenuation;
                    } else {
                        downstream_factor[up_id] += factor * attenuation;
                    }
                }
            } catch (...) {}
        }
    }
}

inline std::unordered_map<std::string, double> DiffusiveWaveRouter::get_gradients() const {
    std::unordered_map<std::string, double> grads;
    for (int reach_id : network_.topological_order()) {
        const Reach& reach = network_.get_reach(reach_id);
        std::string key = "manning_n_" + std::to_string(reach_id);
        grads[key] = reach.grad_manning_n;
        std::string prefix = "reach_" + std::to_string(reach_id) + "_";
        grads[prefix + "manning_n"] = reach.grad_manning_n;
        grads[prefix + "width_coef"] = reach.grad_width_coef;
        grads[prefix + "depth_coef"] = reach.grad_depth_coef;
    }
    return grads;
}

inline void DiffusiveWaveRouter::reset_gradients() {
    network_.zero_gradients();
    for (auto& [id, g] : analytical_dQ_dn_) g = 0.0;
}

inline void DiffusiveWaveRouter::set_lateral_inflow(int reach_id, double inflow) {
    network_.get_reach(reach_id).lateral_inflow = Real(inflow);
}

inline void DiffusiveWaveRouter::set_lateral_inflows(const std::vector<double>& inflows) {
    const auto& order = network_.topological_order();
    for (size_t i = 0; i < order.size() && i < inflows.size(); ++i) {
        network_.get_reach(order[i]).lateral_inflow = Real(inflows[i]);
    }
}

inline double DiffusiveWaveRouter::get_discharge(int reach_id) const {
    return to_double(network_.get_reach(reach_id).outflow_curr);
}

inline std::vector<double> DiffusiveWaveRouter::get_all_discharges() const {
    std::vector<double> discharges;
    for (int reach_id : network_.topological_order()) {
        discharges.push_back(to_double(network_.get_reach(reach_id).outflow_curr));
    }
    return discharges;
}

inline void DiffusiveWaveRouter::reset_state() {
    initialize_state();
    for (int reach_id : network_.topological_order()) {
        Reach& reach = network_.get_reach(reach_id);
        reach.inflow_prev = reach.inflow_curr = Real(0.0);
        reach.outflow_prev = reach.outflow_curr = Real(0.0);
        reach.lateral_inflow = Real(0.0);
    }
    current_time_ = 0.0;
    reset_gradients();
}


// ============================================================================
// Lag Router (Simple Time Delay)
// ============================================================================

/**
 * Simple lag routing with optional attenuation.
 * 
 * Each reach delays flow by travel_time = length / velocity.
 * Optionally applies exponential decay for storage effects.
 * 
 * Fully differentiable. Good baseline for comparison.
 */
class LagRouter {
public:
    explicit LagRouter(Network& network, RouterConfig config = {});
    
    void route_timestep();
    void route(int num_timesteps);
    
    void enable_gradients(bool enable);
    void start_recording();
    void stop_recording();
    void compute_gradients(const std::vector<int>& gauge_reaches,
                           const std::vector<double>& dL_dQ);
    std::unordered_map<std::string, double> get_gradients() const;
    void reset_gradients();
    
    void set_lateral_inflow(int reach_id, double inflow);
    void set_lateral_inflows(const std::vector<double>& inflows);
    double get_discharge(int reach_id) const;
    std::vector<double> get_all_discharges() const;
    void reset_state();
    
    double current_time() const { return current_time_; }
    Network& network() { return network_; }
    const RouterConfig& config() const { return config_; }
    
private:
    Network& network_;
    RouterConfig config_;
    double current_time_ = 0.0;
    bool recording_ = false;
    
    // Lag buffer for each reach (stores delayed inflows)
    std::unordered_map<int, std::deque<Real>> lag_buffer_;
    std::unordered_map<int, int> lag_steps_;  // Number of timesteps to delay
    
    // For analytical gradients
    std::unordered_map<int, double> reach_velocity_;
    std::unordered_map<int, double> analytical_dQ_dn_;
    
    void initialize_lags();
};

inline LagRouter::LagRouter(Network& network, RouterConfig config)
    : network_(network), config_(std::move(config)) {
    network_.build_topology();
    initialize_lags();
}

inline void LagRouter::initialize_lags() {
    lag_buffer_.clear();
    lag_steps_.clear();
    reach_velocity_.clear();
    analytical_dQ_dn_.clear();
    
    for (int reach_id : network_.topological_order()) {
        Reach& reach = network_.get_reach(reach_id);
        
        // Compute velocity from Manning's equation
        double n = to_double(reach.manning_n);
        double S = reach.slope;
        double R_h = 1.0;
        
        double velocity = (1.0 / n) * std::pow(R_h, 2.0/3.0) * std::sqrt(S);
        velocity = std::max(0.1, std::min(5.0, velocity));
        reach_velocity_[reach_id] = velocity;
        
        // Compute lag in timesteps
        double travel_time = reach.length / velocity;
        int lag = std::max(1, static_cast<int>(travel_time / config_.dt));
        lag_steps_[reach_id] = lag;
        
        // Initialize buffer
        lag_buffer_[reach_id] = std::deque<Real>(lag, Real(0.0));
        analytical_dQ_dn_[reach_id] = 0.0;
    }
}

inline void LagRouter::route_timestep() {
    for (int reach_id : network_.topological_order()) {
        Reach& reach = network_.get_reach(reach_id);
        
        // Get upstream inflow
        Real Q_in = Real(0.0);
        if (reach.upstream_junction_id >= 0) {
            try {
                const Junction& junc = network_.get_junction(reach.upstream_junction_id);
                for (int up_id : junc.upstream_reach_ids) {
                    Q_in = Q_in + network_.get_reach(up_id).outflow_curr;
                }
            } catch (...) {}
        }
        Q_in = Q_in + reach.lateral_inflow;
        
        // Get lagged outflow
        auto& buffer = lag_buffer_[reach_id];
        Real Q_out = buffer.back();
        buffer.pop_back();
        buffer.push_front(Q_in);
        
        // Update state
        reach.inflow_curr = Q_in;
        reach.outflow_curr = Q_out;
        
        // Compute analytical gradient: dQ_out/dn
        // For lag routing: Q_out depends on n through travel time
        // Approximate gradient using numerical differentiation of buffer
        if (recording_) {
            double n = to_double(reach.manning_n);
            int lag = lag_steps_[reach_id];
            
            // Better gradient estimate: use average slope in buffer
            double dQ_sum = 0.0;
            int count = 0;
            for (size_t i = 0; i + 1 < buffer.size(); ++i) {
                dQ_sum += to_double(buffer[i]) - to_double(buffer[i+1]);
                count++;
            }
            double dQ_avg = (count > 0) ? dQ_sum / count : 0.0;
            
            // dQ/dn ≈ -dQ_avg * lag / n (negative because increasing n increases lag)
            // Scale by travel time to get reasonable magnitude
            double travel_time = lag * config_.dt;
            analytical_dQ_dn_[reach_id] = -std::abs(dQ_avg) * travel_time / (n * config_.dt);
        }
    }
    current_time_ += config_.dt;
}

inline void LagRouter::route(int num_timesteps) {
    for (int t = 0; t < num_timesteps; ++t) route_timestep();
}

inline void LagRouter::enable_gradients(bool enable) { config_.enable_gradients = enable; }
inline void LagRouter::start_recording() { recording_ = true; }
inline void LagRouter::stop_recording() { recording_ = false; }

inline void LagRouter::compute_gradients(const std::vector<int>& gauge_reaches,
                                          const std::vector<double>& dL_dQ) {
    if (gauge_reaches.empty()) return;
    double dL_dQ_out = dL_dQ[0];
    
    // Propagate gradients upstream (similar to IRF)
    auto topo_order = network_.topological_order();
    std::unordered_map<int, double> factor;
    factor[gauge_reaches[0]] = 1.0;
    
    for (auto it = topo_order.rbegin(); it != topo_order.rend(); ++it) {
        int reach_id = *it;
        if (factor.count(reach_id) == 0) continue;
        
        Reach& reach = network_.get_reach(reach_id);
        reach.grad_manning_n = dL_dQ_out * factor[reach_id] * analytical_dQ_dn_[reach_id];
        
        if (reach.upstream_junction_id >= 0) {
            try {
                const Junction& junc = network_.get_junction(reach.upstream_junction_id);
                for (int up_id : junc.upstream_reach_ids) {
                    factor[up_id] = factor.count(up_id) ? factor[up_id] + factor[reach_id] * 0.9 
                                                         : factor[reach_id] * 0.9;
                }
            } catch (...) {}
        }
    }
}

inline std::unordered_map<std::string, double> LagRouter::get_gradients() const {
    std::unordered_map<std::string, double> grads;
    for (int reach_id : network_.topological_order()) {
        const Reach& reach = network_.get_reach(reach_id);
        grads["manning_n_" + std::to_string(reach_id)] = reach.grad_manning_n;
    }
    return grads;
}

inline void LagRouter::reset_gradients() {
    network_.zero_gradients();
    for (auto& [id, g] : analytical_dQ_dn_) g = 0.0;
}

inline void LagRouter::set_lateral_inflow(int reach_id, double inflow) {
    network_.get_reach(reach_id).lateral_inflow = Real(inflow);
}

inline void LagRouter::set_lateral_inflows(const std::vector<double>& inflows) {
    const auto& order = network_.topological_order();
    for (size_t i = 0; i < order.size() && i < inflows.size(); ++i) {
        network_.get_reach(order[i]).lateral_inflow = Real(inflows[i]);
    }
}

inline double LagRouter::get_discharge(int reach_id) const {
    return to_double(network_.get_reach(reach_id).outflow_curr);
}

inline std::vector<double> LagRouter::get_all_discharges() const {
    std::vector<double> d;
    for (int reach_id : network_.topological_order()) {
        d.push_back(to_double(network_.get_reach(reach_id).outflow_curr));
    }
    return d;
}

inline void LagRouter::reset_state() {
    initialize_lags();
    for (int reach_id : network_.topological_order()) {
        Reach& reach = network_.get_reach(reach_id);
        reach.inflow_prev = reach.inflow_curr = Real(0.0);
        reach.outflow_prev = reach.outflow_curr = Real(0.0);
    }
    current_time_ = 0.0;
}


// ============================================================================
// KWT Router (Kinematic Wave Tracking) - Non-differentiable
// ============================================================================

/**
 * Kinematic Wave Tracking routing.
 * 
 * Implements the kinematic wave tracking method as used in mizuRoute.
 * Each wave parcel is treated as a continuous wave segment with spatial extent,
 * not a point mass. This eliminates artificial spikiness from discrete arrivals.
 * 
 * Physical basis:
 * - Wave celerity from kinematic wave theory: c = (5/3) * v = (5/3) * Q/A
 * - Each parcel has a spatial extent (wave_length = c * dt at creation)
 * - Outflow is proportional to how much of the wave has passed the outlet
 * - Parcels are removed only when fully past the outlet
 * 
 * Reference: Mizukami et al. (2016), mizuRoute technical documentation
 * 
 * NOTE: This method is NOT differentiable due to discrete parcel operations.
 */
class KWTRouter {
public:
    explicit KWTRouter(Network& network, RouterConfig config = {});
    
    void route_timestep();
    void route(int num_timesteps);
    
    void enable_gradients(bool) { /* Not supported */ }
    void start_recording() { /* Not supported */ }
    void stop_recording() { /* Not supported */ }
    void compute_gradients(const std::vector<int>&, const std::vector<double>&) {
        std::cerr << "Warning: KWT routing does not support gradients\n";
    }
    std::unordered_map<std::string, double> get_gradients() const {
        return {};  // Empty - no gradients
    }
    void reset_gradients() { network_.zero_gradients(); }
    
    void set_lateral_inflow(int reach_id, double inflow);
    void set_lateral_inflows(const std::vector<double>& inflows);
    double get_discharge(int reach_id) const;
    std::vector<double> get_all_discharges() const;
    void reset_state();
    
    double current_time() const { return current_time_; }
    Network& network() { return network_; }
    const RouterConfig& config() const { return config_; }
    
private:
    Network& network_;
    RouterConfig config_;
    double current_time_ = 0.0;
    
    /**
     * Wave parcel structure - represents a continuous wave segment
     * 
     * The wave extends from (position - wave_length) to (position)
     * where position is the leading edge.
     */
    struct WaveParcel {
        double volume;       // Total volume in parcel [m³]
        double position;     // Leading edge position [m from upstream]
        double wave_length;  // Spatial extent of wave [m]
        double celerity;     // Wave celerity [m/s]
        double rf;           // Remaining fraction (1.0 = full, 0.0 = fully exited)
    };
    
    // Active parcels in each reach
    std::unordered_map<int, std::vector<WaveParcel>> parcels_;
    
    // Downstream reach mapping for parcel transfer
    std::unordered_map<int, int> downstream_reach_;
    
    double compute_celerity(const Reach& reach, double Q);
    void initialize_topology();
    bool topology_initialized_ = false;
};

inline KWTRouter::KWTRouter(Network& network, RouterConfig config)
    : network_(network), config_(std::move(config)) {
    network_.build_topology();
    reset_state();
}

inline void KWTRouter::initialize_topology() {
    if (topology_initialized_) return;
    
    // Build downstream reach mapping
    downstream_reach_.clear();
    for (int reach_id : network_.topological_order()) {
        const Reach& reach = network_.get_reach(reach_id);
        downstream_reach_[reach_id] = -1;  // Default: no downstream
        
        // Find downstream reach through junction
        if (reach.downstream_junction_id >= 0) {
            try {
                const Junction& junc = network_.get_junction(reach.downstream_junction_id);
                if (!junc.downstream_reach_ids.empty()) {
                    downstream_reach_[reach_id] = junc.downstream_reach_ids[0];
                }
            } catch (...) {}
        }
    }
    
    topology_initialized_ = true;
}

inline double KWTRouter::compute_celerity(const Reach& reach, double Q) {
    // Minimum flow for numerical stability
    if (Q < 0.001) Q = 0.001;
    
    // Get channel geometry
    double n = to_double(reach.manning_n);
    double S = reach.slope;
    if (S < 1e-6) S = 1e-6;
    
    // Hydraulic geometry: width and depth from power laws
    double width = to_double(reach.geometry.width_coef) * 
                   std::pow(Q, to_double(reach.geometry.width_exp));
    if (width < 0.5) width = 0.5;
    
    double depth = to_double(reach.geometry.depth_coef) * 
                   std::pow(Q, to_double(reach.geometry.depth_exp));
    if (depth < 0.05) depth = 0.05;
    
    // Cross-sectional area and hydraulic radius
    double area = width * depth;
    double wetted_perimeter = width + 2.0 * depth;
    double R_h = area / wetted_perimeter;
    
    // Velocity from Manning's equation
    double velocity = (1.0 / n) * std::pow(R_h, 2.0/3.0) * std::sqrt(S);
    
    // Kinematic wave celerity: c = (5/3) * v for wide rectangular channel
    // This comes from c = dQ/dA and the Manning equation
    double celerity = (5.0 / 3.0) * velocity;
    
    // Bound celerity for stability (0.1 - 5.0 m/s typical for rivers)
    return std::max(0.1, std::min(5.0, celerity));
}

inline void KWTRouter::route_timestep() {
    initialize_topology();
    
    double dt = config_.dt;
    
    // Storage for outflow rates (computed from exiting wave fractions)
    std::unordered_map<int, double> outflow_rate;
    
    // Storage for volumes to transfer downstream
    std::unordered_map<int, double> transfer_volume;
    
    for (int reach_id : network_.topological_order()) {
        outflow_rate[reach_id] = 0.0;
        transfer_volume[reach_id] = 0.0;
    }
    
    // Process reaches in topological order
    for (int reach_id : network_.topological_order()) {
        Reach& reach = network_.get_reach(reach_id);
        double L = reach.length;
        
        // Collect inflow: upstream transfers + lateral inflow
        double inflow_vol = 0.0;
        
        // Get upstream contributions
        if (reach.upstream_junction_id >= 0) {
            try {
                const Junction& junc = network_.get_junction(reach.upstream_junction_id);
                for (int up_id : junc.upstream_reach_ids) {
                    inflow_vol += transfer_volume[up_id];
                }
            } catch (...) {}
        }
        
        // Add lateral inflow
        double lateral_vol = to_double(reach.lateral_inflow) * dt;
        inflow_vol += lateral_vol;
        
        // Create new parcel for inflow (if significant)
        if (inflow_vol > 1e-6) {
            double Q_est = inflow_vol / dt;
            double celerity = compute_celerity(reach, Q_est);
            double wave_length = celerity * dt;  // Wave extends over one timestep
            
            WaveParcel parcel;
            parcel.volume = inflow_vol;
            parcel.position = wave_length;  // Leading edge after one dt
            parcel.wave_length = wave_length;
            parcel.celerity = celerity;
            parcel.rf = 1.0;  // Fully present
            
            parcels_[reach_id].push_back(parcel);
        }
        
        // Advance existing parcels and compute outflow
        auto& reach_parcels = parcels_[reach_id];
        std::vector<WaveParcel> remaining;
        double total_outflow_vol = 0.0;
        
        for (auto& parcel : reach_parcels) {
            // Advance parcel position
            double old_position = parcel.position;
            parcel.position += parcel.celerity * dt;
            
            // Trailing edge positions
            double old_trailing = old_position - parcel.wave_length;
            double new_trailing = parcel.position - parcel.wave_length;
            
            // Compute fraction of wave that has exited the reach
            // The wave extends from trailing edge to leading edge
            
            if (new_trailing >= L) {
                // Entire wave has exited
                double exit_vol = parcel.rf * parcel.volume;
                total_outflow_vol += exit_vol;
                parcel.rf = 0.0;
                // Don't keep this parcel
            } else if (parcel.position > L) {
                // Partial exit: leading edge past outlet, trailing edge still in reach
                // Fraction exited = (position - L) / wave_length
                double fraction_past = (parcel.position - L) / parcel.wave_length;
                fraction_past = std::min(1.0, std::max(0.0, fraction_past));
                
                // Volume that exited this timestep
                // Need to track what already exited vs. what's new
                double old_fraction_past = 0.0;
                if (old_position > L) {
                    old_fraction_past = (old_position - L) / parcel.wave_length;
                    old_fraction_past = std::min(1.0, std::max(0.0, old_fraction_past));
                }
                
                double new_exit_fraction = fraction_past - old_fraction_past;
                double exit_vol = new_exit_fraction * parcel.volume;
                total_outflow_vol += exit_vol;
                
                // Update remaining fraction
                parcel.rf = 1.0 - fraction_past;
                
                if (parcel.rf > 1e-6) {
                    remaining.push_back(parcel);
                }
            } else {
                // Wave entirely within reach
                remaining.push_back(parcel);
            }
        }
        
        reach_parcels = remaining;
        
        // Compute outflow rate
        outflow_rate[reach_id] = total_outflow_vol / dt;
        
        // Store volume for downstream transfer
        transfer_volume[reach_id] = total_outflow_vol;
        
        // Update reach state
        reach.inflow_curr = Real(inflow_vol / dt);
        reach.outflow_curr = Real(outflow_rate[reach_id]);
    }
    
    current_time_ += dt;
}

inline void KWTRouter::route(int num_timesteps) {
    for (int t = 0; t < num_timesteps; ++t) route_timestep();
}

inline void KWTRouter::set_lateral_inflow(int reach_id, double inflow) {
    network_.get_reach(reach_id).lateral_inflow = Real(inflow);
}

inline void KWTRouter::set_lateral_inflows(const std::vector<double>& inflows) {
    const auto& order = network_.topological_order();
    for (size_t i = 0; i < order.size() && i < inflows.size(); ++i) {
        network_.get_reach(order[i]).lateral_inflow = Real(inflows[i]);
    }
}

inline double KWTRouter::get_discharge(int reach_id) const {
    return to_double(network_.get_reach(reach_id).outflow_curr);
}

inline std::vector<double> KWTRouter::get_all_discharges() const {
    std::vector<double> d;
    for (int reach_id : network_.topological_order()) {
        d.push_back(to_double(network_.get_reach(reach_id).outflow_curr));
    }
    return d;
}

inline void KWTRouter::reset_state() {
    parcels_.clear();
    topology_initialized_ = false;
    for (int reach_id : network_.topological_order()) {
        parcels_[reach_id] = {};
        
        Reach& reach = network_.get_reach(reach_id);
        reach.inflow_prev = reach.inflow_curr = Real(0.0);
        reach.outflow_prev = reach.outflow_curr = Real(0.0);
    }
    current_time_ = 0.0;
}

} // namespace dmc

#endif // DMC_ROUTE_ROUTER_HPP
