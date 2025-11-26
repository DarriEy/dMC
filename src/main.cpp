/**
 * dMC-Route: Differentiable Muskingum-Cunge Routing
 * 
 * Standalone executable for testing and demonstration.
 */

#include <dmc/types.hpp>
#include <dmc/network.hpp>
#include <dmc/network_io.hpp>
#include <dmc/topology_nc.hpp>
#include <dmc/runoff_forcing.hpp>
#include <dmc/router.hpp>
#include <dmc/bmi.hpp>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <memory>
#include <cmath>
#include <getopt.h>

using namespace dmc;

void print_banner() {
    std::cout << R"(
╔═══════════════════════════════════════════════════════════════╗
║  dMC-Route: Differentiable Muskingum-Cunge Routing v0.1.0     ║
║  Automatic differentiation for hydrological parameter learning║
╚═══════════════════════════════════════════════════════════════╝
)" << std::endl;
}

void print_usage(const char* prog) {
    std::cout << "Usage: " << prog << " [OPTIONS]\n\n"
              << "Options:\n"
              << "  -n, --network FILE     Network file (topology.nc, GeoJSON, or CSV)\n"
              << "  -f, --forcing FILE     Runoff forcing file (NetCDF or CSV)\n"
<<<<<<< HEAD
              << "  -c, --config FILE      Forcing configuration file (mizuRoute-style)\n"
=======
              << "  -c, --config FILE      Control file (can specify all options)\n"
>>>>>>> b5bae52 (Initial commit of dmc_design codebase)
              << "  -p, --preset NAME      Use preset config: summa, fuse, gr4j\n"
              << "  -o, --output FILE      Output file for results (CSV)\n"
              << "  -t, --timestep SECS    Timestep in seconds (default: 3600)\n"
              << "  -g, --gradients        Enable gradient computation\n"
<<<<<<< HEAD
              << "  -j, --jacobian FILE    Output Jacobian matrix to FILE (requires -g)\n"
=======
              << "  -j, --jacobian FILE    Output Jacobian matrix to FILE (enables -g)\n"
>>>>>>> b5bae52 (Initial commit of dmc_design codebase)
              << "  -d, --demo             Run built-in demo\n"
              << "  -h, --help             Show this help\n"
              << "\n"
              << "Examples:\n"
              << "  " << prog << " --demo\n"
              << "  " << prog << " -n topology.nc -f summa_output.nc -p summa -o discharge.csv\n"
              << "  " << prog << " -n topology.nc -f fuse_output.nc -c forcing_config.txt\n"
<<<<<<< HEAD
              << "  " << prog << " -n topology.nc -f summa_output.nc -p summa -g -j jacobian.csv\n"
=======
              << "  " << prog << " -c full_config.txt                 (all options in config)\n"
              << "  " << prog << " -c config.txt -j jacobian.csv      (config + CLI override)\n"
>>>>>>> b5bae52 (Initial commit of dmc_design codebase)
              << std::endl;
}

/**
 * Run simulation from files.
 */
void run_from_files(const std::string& network_file,
                    const std::string& forcing_file,
                    const std::string& config_file,
                    const std::string& preset,
                    const std::string& output_file,
                    const std::string& jacobian_file,
                    double dt_override,
                    bool enable_gradients) {
    
    std::cout << "Loading network from: " << network_file << "\n";
    
    Network network;
    HRUMapping hru_mapping;
    bool has_mapping = false;
    
    // Check if this is a topology.nc file (mizuRoute format)
    bool is_topology_nc = (network_file.find("topology") != std::string::npos &&
                           network_file.find(".nc") != std::string::npos);
    
    if (is_topology_nc) {
#ifdef DMC_USE_NETCDF
        // Load from mizuRoute topology.nc
        TopologyNCReader topo_reader;
        network = topo_reader.load_network(network_file);
        
        // Extract HRU mapping from topology
        for (const auto& info : topo_reader.get_hru_info()) {
            hru_mapping.hru_to_reach[info.hru_id] = info.segment_id;
            hru_mapping.hru_areas[info.hru_id] = info.area_m2;
        }
        has_mapping = true;
        
        std::cout << "  Loaded topology.nc: " << network.num_reaches() << " reaches, "
                  << hru_mapping.hru_to_reach.size() << " HRUs\n";
#else
        throw std::runtime_error("NetCDF support required for topology.nc files");
#endif
    } else {
        // Load from GeoJSON or CSV
        NetworkIO io;
        if (network_file.find(".geojson") != std::string::npos ||
            network_file.find(".json") != std::string::npos) {
            network = io.load_geojson(network_file);
        } else {
            network = io.load_csv(network_file);
        }
        std::cout << "  Loaded " << network.num_reaches() << " reaches, "
                  << network.num_junctions() << " junctions\n";
    }
    
    // Load forcing configuration
    RunoffForcingConfig forcing_config;
    
    if (!config_file.empty()) {
        std::cout << "Loading forcing config from: " << config_file << "\n";
        forcing_config = RunoffForcingConfig::load_from_control_file(config_file);
    } else if (!preset.empty()) {
        std::cout << "Using preset: " << preset << "\n";
        if (preset == "summa") {
            forcing_config = RunoffForcingConfig::summa_preset();
        } else if (preset == "fuse") {
            forcing_config = RunoffForcingConfig::fuse_preset();
        } else if (preset == "gr4j") {
            forcing_config = RunoffForcingConfig::gr4j_preset();
        } else {
            std::cerr << "Unknown preset: " << preset << ". Using defaults.\n";
        }
    }
    
    // Override timestep if specified
    double dt = (dt_override > 0) ? dt_override : forcing_config.dt_qsim;
    
    // Load forcing
    std::cout << "Loading forcing from: " << forcing_file << "\n";
    
    LateralInflowData forcing_data;
    
    if (forcing_file.find(".nc") != std::string::npos) {
#ifdef DMC_USE_NETCDF
        RunoffForcingReader reader(forcing_config);
        
        if (has_mapping) {
            reader.set_mapping(hru_mapping);
        } else {
            // Create identity mapping
            std::vector<int> reach_ids(network.topological_order().begin(),
                                       network.topological_order().end());
            reader.set_mapping(HRUMapping::create_identity(reach_ids));
        }
        
        forcing_data = reader.read(forcing_file);
        std::cout << "  Loaded " << forcing_data.num_timesteps() << " timesteps\n";
#else
        throw std::runtime_error("NetCDF support required for .nc forcing files");
#endif
    } else {
        // CSV forcing - simple format: time,reach1,reach2,...
        std::ifstream file(forcing_file);
        if (!file.is_open()) {
            throw std::runtime_error("Cannot open forcing file: " + forcing_file);
        }
        
        std::string line;
        bool header = true;
        std::vector<int> reach_ids;
        
        while (std::getline(file, line)) {
            if (line.empty()) continue;
            
            std::stringstream ss(line);
            std::string token;
            std::vector<std::string> tokens;
            
            while (std::getline(ss, token, ',')) {
                tokens.push_back(token);
            }
            
            if (header) {
                header = false;
                // Parse reach IDs from header
                for (size_t i = 1; i < tokens.size(); ++i) {
                    // Extract numeric part from column name like "Q_123" or just "123"
                    std::string col = tokens[i];
                    size_t pos = col.find_last_of("_");
                    if (pos != std::string::npos) {
                        col = col.substr(pos + 1);
                    }
                    reach_ids.push_back(std::stoi(col));
                }
                continue;
            }
            
            // Parse data row
            double time_val = std::stod(tokens[0]);
            forcing_data.times.push_back(time_val);
            
            for (size_t i = 1; i < tokens.size() && i - 1 < reach_ids.size(); ++i) {
                int reach_id = reach_ids[i - 1];
                double inflow = std::stod(tokens[i]);
                forcing_data.reach_inflows[reach_id].push_back(inflow);
            }
        }
        
        std::cout << "  Loaded CSV: " << forcing_data.num_timesteps() << " timesteps\n";
    }
    
    if (forcing_data.num_timesteps() == 0) {
        std::cerr << "Error: No forcing data loaded (0 timesteps)\n";
        return;
    }
    
    // Setup router
    RouterConfig router_config;
    router_config.dt = dt;
    router_config.enable_gradients = enable_gradients && AD_ENABLED;
    
    MuskingumCungeRouter router(network, router_config);
    
    std::cout << "\nRunning " << forcing_data.num_timesteps() 
              << " timesteps (dt=" << dt << "s)...\n";
    
    if (router_config.enable_gradients) {
        router.start_recording();
    }
    
    // Open output file
    std::ofstream outfile;
    if (!output_file.empty()) {
        outfile.open(output_file);
        outfile << "time";
        for (int reach_id : network.topological_order()) {
            outfile << ",Q_" << reach_id;
        }
        outfile << "\n";
    }
    
    // Run simulation
    size_t num_timesteps = forcing_data.num_timesteps();
    size_t print_interval = std::max(size_t(1), num_timesteps / 10);
    
    for (size_t t = 0; t < num_timesteps; ++t) {
        // Set lateral inflows
        for (int reach_id : network.topological_order()) {
            double inflow = forcing_data.get_inflow(reach_id, t);
            router.set_lateral_inflow(reach_id, inflow);
        }
        
        // Route
        router.route_timestep();
        
        // Write output
        if (outfile.is_open()) {
            double time_val = (t < forcing_data.times.size()) ? forcing_data.times[t] : t * dt;
            outfile << time_val;
            for (int reach_id : network.topological_order()) {
                outfile << "," << router.get_discharge(reach_id);
            }
            outfile << "\n";
        }
        
        // Progress
        if (t % print_interval == 0) {
            double progress = 100.0 * t / num_timesteps;
            std::cout << "  Progress: " << std::fixed << std::setprecision(0) 
                      << progress << "%\r" << std::flush;
        }
    }
    
    std::cout << "  Progress: 100%    \n";
    
    if (router_config.enable_gradients) {
        router.stop_recording();
        
        // Find outlets and compute gradients
        std::vector<int> outlets;
        for (int reach_id : network.topological_order()) {
            const Reach& r = network.get_reach(reach_id);
            if (r.downstream_junction_id < 0) {
                outlets.push_back(reach_id);
            }
        }
        
        if (!outlets.empty()) {
            std::vector<double> dL_dQ(outlets.size(), 1.0);
            router.compute_gradients(outlets, dL_dQ);
            
            auto grads = router.get_gradients();
            
            // Print summary to console
            std::cout << "\n=== Parameter Gradients (∂Q_outlet/∂θ) ===\n";
            std::cout << "Outlet reach(es): ";
            for (int o : outlets) std::cout << o << " ";
            std::cout << "\n\n";
            
            for (int reach_id : network.topological_order()) {
                std::string prefix = "reach_" + std::to_string(reach_id) + "_";
                double grad_n = grads[prefix + "manning_n"];
                double grad_w = grads[prefix + "width_coef"];
                
                // Only print non-zero gradients
                if (std::abs(grad_n) > 1e-15 || std::abs(grad_w) > 1e-15) {
                    std::cout << "Reach " << std::setw(3) << reach_id << ": ";
                    std::cout << "∂Q/∂n = " << std::setw(12) << std::scientific << grad_n;
                    std::cout << ", ∂Q/∂w = " << std::setw(12) << grad_w;
                    std::cout << "\n";
                }
            }
            
            // Write Jacobian to file if requested
            if (!jacobian_file.empty()) {
                std::ofstream jac_file(jacobian_file);
                if (jac_file.is_open()) {
                    // Header
                    jac_file << "reach_id,manning_n,width_coef,width_exp,depth_coef,depth_exp\n";
                    
                    // Data rows
                    for (int reach_id : network.topological_order()) {
                        std::string prefix = "reach_" + std::to_string(reach_id) + "_";
                        jac_file << reach_id;
                        jac_file << "," << std::scientific << std::setprecision(10) << grads[prefix + "manning_n"];
                        jac_file << "," << grads[prefix + "width_coef"];
                        jac_file << "," << grads[prefix + "width_exp"];
                        jac_file << "," << grads[prefix + "depth_coef"];
                        jac_file << "," << grads[prefix + "depth_exp"];
                        jac_file << "\n";
                    }
                    
                    jac_file.close();
                    std::cout << "\nJacobian written to: " << jacobian_file << "\n";
                } else {
                    std::cerr << "Warning: Could not open Jacobian file: " << jacobian_file << "\n";
                }
            }
        } else {
            std::cout << "Warning: No outlet reaches found - cannot compute gradients\n";
        }
    }
    
    if (outfile.is_open()) {
        outfile.close();
        std::cout << "\nResults written to: " << output_file << "\n";
    }
}

Network create_demo_network() {
    Network net;
    
    // Create a simple Y-shaped network:
    //
    //   reach_0 (headwater 1)
    //                         \
    //                          junction_2 --- reach_3 (main stem) --- outlet
    //                         /
    //   reach_1 (headwater 2)
    //                         |
    //   reach_2 (tributary)  -+
    
    // Headwater reaches
    for (int i = 0; i < 2; ++i) {
        Reach r;
        r.id = i;
        r.name = "headwater_" + std::to_string(i);
        r.length = 3000.0;  // 3 km
        r.slope = 0.002;    // Steeper headwaters
        r.manning_n = Real(0.04);
        r.upstream_junction_id = -1;  // Headwater
        r.downstream_junction_id = 2; // Both flow to junction 2
        net.add_reach(r);
        
        Junction j;
        j.id = i;
        j.is_headwater = true;
        j.downstream_reach_ids = {i};
        net.add_junction(j);
    }
    
    // Tributary
    Reach r2;
    r2.id = 2;
    r2.name = "tributary";
    r2.length = 4000.0;
    r2.slope = 0.0015;
    r2.manning_n = Real(0.035);
    r2.upstream_junction_id = -1;
    r2.downstream_junction_id = 2;
    net.add_reach(r2);
    
    // Confluence junction
    Junction j2;
    j2.id = 2;
    j2.upstream_reach_ids = {0, 1, 2};
    j2.downstream_reach_ids = {3};
    net.add_junction(j2);
    
    // Main stem
    Reach r3;
    r3.id = 3;
    r3.name = "mainstem";
    r3.length = 10000.0;  // 10 km
    r3.slope = 0.0005;    // Gentler slope
    r3.manning_n = Real(0.03);
    r3.upstream_junction_id = 2;
    r3.downstream_junction_id = 3;
    net.add_reach(r3);
    
    // Outlet junction
    Junction j3;
    j3.id = 3;
    j3.upstream_reach_ids = {3};
    j3.is_outlet = true;
    net.add_junction(j3);
    
    net.build_topology();
    return net;
}

std::vector<double> generate_storm_hyetograph(int num_steps, double peak_intensity, 
                                               int peak_step, int duration) {
    std::vector<double> inflow(num_steps, 0.0);
    
    // Simple triangular storm
    int start = peak_step - duration / 2;
    int end = peak_step + duration / 2;
    
    for (int t = start; t <= end && t < num_steps; ++t) {
        if (t < 0) continue;
        double progress;
        if (t <= peak_step) {
            progress = double(t - start) / (peak_step - start);
        } else {
            progress = 1.0 - double(t - peak_step) / (end - peak_step);
        }
        inflow[t] = peak_intensity * progress;
    }
    
    return inflow;
}

void run_demo_simulation() {
    std::cout << "Creating demo network...\n";
    Network net = create_demo_network();
    
    std::cout << "Network created:\n";
    std::cout << "  Reaches: " << net.num_reaches() << "\n";
    std::cout << "  Junctions: " << net.num_junctions() << "\n";
    std::cout << "  Topological order: ";
    for (int id : net.topological_order()) {
        std::cout << id << " ";
    }
    std::cout << "\n\n";
    
    // Configuration
    RouterConfig config;
    config.dt = 900.0;  // 15-minute timestep
    config.enable_gradients = AD_ENABLED;
    
    MuskingumCungeRouter router(net, config);
    
    // Generate storm event (48 hours, peak at 12 hours)
    int num_steps = 48 * 4;  // 48 hours at 15-min intervals
    double peak_intensity = 5.0;  // m³/s per reach
    int peak_step = 12 * 4;
    int storm_duration = 8 * 4;  // 8 hours
    
    auto storm = generate_storm_hyetograph(num_steps, peak_intensity, 
                                            peak_step, storm_duration);
    
    std::cout << "Running " << num_steps << " timesteps (" 
              << num_steps * config.dt / 3600.0 << " hours)...\n\n";
    
    // Start recording for gradients
    if (AD_ENABLED) {
        router.start_recording();
    }
    
    // Run simulation
    std::cout << std::fixed << std::setprecision(3);
    std::cout << "Time (hr)  |  Inflow  |  Outlet Q  |  Velocity  |  K (hr)  |  X\n";
    std::cout << "---------------------------------------------------------------\n";
    
    double max_Q = 0.0;
    
    for (int t = 0; t < num_steps; ++t) {
        // Apply lateral inflow to headwater reaches
        for (int i = 0; i < 3; ++i) {
            router.set_lateral_inflow(i, storm[t]);
        }
        
        // Route
        router.route_timestep();
        
        // Get outlet discharge
        double Q_out = router.get_discharge(3);
        max_Q = std::max(max_Q, Q_out);
        
        // Print every 4 hours
        if (t % 16 == 0 || t == num_steps - 1) {
            const Reach& outlet_reach = net.get_reach(3);
            double hours = router.current_time() / 3600.0;
            std::cout << std::setw(8) << hours << "  |  "
                      << std::setw(7) << storm[t] << "  |  "
                      << std::setw(9) << Q_out << "  |  "
                      << std::setw(9) << to_double(outlet_reach.velocity) << "  |  "
                      << std::setw(7) << to_double(outlet_reach.K) / 3600.0 << "  |  "
                      << std::setw(5) << to_double(outlet_reach.X) << "\n";
        }
    }
    
    std::cout << "\nPeak outlet discharge: " << max_Q << " m³/s\n";
    
    // Compute gradients
    if (AD_ENABLED) {
        router.stop_recording();
        
        std::cout << "\n=== Gradient Computation ===\n";
        
        // Seed with unit gradient at outlet
        std::vector<int> gauge = {3};
        std::vector<double> dL_dQ = {1.0};
        router.compute_gradients(gauge, dL_dQ);
        
        auto grads = router.get_gradients();
        
        std::cout << "Sensitivity of outlet Q to parameters:\n";
        for (int reach_id : net.topological_order()) {
            std::string prefix = "reach_" + std::to_string(reach_id) + "_";
            double dn = grads[prefix + "manning_n"];
            double dw = grads[prefix + "width_coef"];
            
            std::cout << "  Reach " << reach_id << ": "
                      << "∂Q/∂n = " << std::setw(10) << dn << ", "
                      << "∂Q/∂w = " << std::setw(10) << dw << "\n";
        }
    } else {
        std::cout << "\n(Gradient computation disabled - compile with DMC_USE_CODIPACK)\n";
    }
}

void run_bmi_demo() {
    std::cout << "\n=== BMI Interface Demo ===\n\n";
    
    BmiMuskingumCunge model;
    model.Initialize("config.yaml");  // Uses default test network
    
    std::cout << "Component: " << model.GetComponentName() << "\n";
    std::cout << "Time step: " << model.GetTimeStep() << " s\n";
    std::cout << "Input vars: ";
    for (const auto& v : model.GetInputVarNames()) std::cout << v << " ";
    std::cout << "\nOutput vars: ";
    for (const auto& v : model.GetOutputVarNames()) std::cout << v << " ";
    std::cout << "\n\n";
    
    // Run 10 steps
    std::vector<double> lateral(3, 2.0);  // 2 m³/s per reach
    model.SetValue("lateral_inflow", lateral.data());
    
    for (int i = 0; i < 10; ++i) {
        model.Update();
    }
    
    std::vector<double> Q(3);
    model.GetValue("discharge", Q.data());
    
    std::cout << "After 10 timesteps, discharge: ";
    for (double q : Q) std::cout << q << " ";
    std::cout << "m³/s\n";
    
    model.Finalize();
}

int main(int argc, char* argv[]) {
    print_banner();
    
    // Command line options
    std::string network_file;
    std::string forcing_file;
    std::string config_file;
    std::string preset;
    std::string output_file;
    std::string jacobian_file;
    double dt_override = -1;
    bool enable_gradients = false;
    bool run_demo = false;
    
    static struct option long_options[] = {
        {"network",   required_argument, 0, 'n'},
        {"forcing",   required_argument, 0, 'f'},
        {"config",    required_argument, 0, 'c'},
        {"preset",    required_argument, 0, 'p'},
        {"output",    required_argument, 0, 'o'},
        {"jacobian",  required_argument, 0, 'j'},
        {"timestep",  required_argument, 0, 't'},
        {"gradients", no_argument,       0, 'g'},
        {"demo",      no_argument,       0, 'd'},
        {"help",      no_argument,       0, 'h'},
        {0, 0, 0, 0}
    };
    
    int opt;
    while ((opt = getopt_long(argc, argv, "n:f:c:p:o:j:t:gdh", long_options, nullptr)) != -1) {
        switch (opt) {
            case 'n': network_file = optarg; break;
            case 'f': forcing_file = optarg; break;
            case 'c': config_file = optarg; break;
            case 'p': preset = optarg; break;
            case 'o': output_file = optarg; break;
            case 'j': jacobian_file = optarg; enable_gradients = true; break;
            case 't': dt_override = std::stod(optarg); break;
            case 'g': enable_gradients = true; break;
            case 'd': run_demo = true; break;
            case 'h':
                print_usage(argv[0]);
                return 0;
            default:
                print_usage(argv[0]);
                return 1;
        }
    }
    
    std::cout << "AD Status: " << (AD_ENABLED ? "ENABLED (CoDiPack)" : "DISABLED") << "\n";
#ifdef DMC_USE_NETCDF
    std::cout << "NetCDF:    ENABLED\n";
#else
    std::cout << "NetCDF:    DISABLED\n";
#endif
    std::cout << "\n";
    
    try {
<<<<<<< HEAD
=======
        // If config file provided, load it and use values as defaults
        // Command line arguments override config file values
        if (!config_file.empty()) {
            RunoffForcingConfig cfg = RunoffForcingConfig::load_from_control_file(config_file);
            
            // Use config file values if not overridden by command line
            if (network_file.empty()) network_file = cfg.network_file;
            if (forcing_file.empty()) forcing_file = cfg.forcing_file;
            if (output_file.empty()) output_file = cfg.output_file;
            if (jacobian_file.empty()) jacobian_file = cfg.jacobian_file;
            if (!enable_gradients) enable_gradients = cfg.enable_gradients;
            if (dt_override < 0 && cfg.dt > 0) dt_override = cfg.dt;
            
            // Enable gradients if jacobian file specified in config
            if (!jacobian_file.empty()) enable_gradients = true;
        }
        
>>>>>>> b5bae52 (Initial commit of dmc_design codebase)
        if (!network_file.empty() && !forcing_file.empty()) {
            // Run from files
            run_from_files(network_file, forcing_file, config_file, preset,
                          output_file, jacobian_file, dt_override, enable_gradients);
        } else if (run_demo || argc == 1) {
            // Run built-in demo
            run_demo_simulation();
            run_bmi_demo();
<<<<<<< HEAD
        } else {
            std::cerr << "Error: Must specify both --network and --forcing, or use --demo\n\n";
=======
        } else if (!config_file.empty()) {
            std::cerr << "Error: Config file must specify <fname_ntopo> and <fname_qsim>\n\n";
            print_usage(argv[0]);
            return 1;
        } else {
            std::cerr << "Error: Must specify both --network and --forcing, or use --config, or use --demo\n\n";
>>>>>>> b5bae52 (Initial commit of dmc_design codebase)
            print_usage(argv[0]);
            return 1;
        }
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
    
    return 0;
}
