#ifndef DMC_ROUTE_TYPES_HPP
#define DMC_ROUTE_TYPES_HPP

// Configuration-based AD type selection
#ifdef DMC_USE_CODIPACK
    #include <codi.hpp>
    
    namespace dmc {
        // Primary type for gradient computation
        using Real = codi::RealReverse;
        using Tape = typename Real::Tape;
        
        // Alternative types for different use cases
        using RealForward = codi::RealForward;       // Forward mode (few params)
        using RealReverseIdx = codi::RealReverseIndex; // Memory-efficient reverse
        
        constexpr bool AD_ENABLED = true;
        
        inline Tape& get_tape() {
            return Real::getTape();
        }
        
        inline void activate_tape() {
            get_tape().setActive();
        }
        
        inline void deactivate_tape() {
            get_tape().setPassive();
        }
        
        inline void reset_tape() {
            get_tape().reset();
        }
        
        template<typename T>
        inline void register_input(T& val) {
            get_tape().registerInput(val);
        }
        
        template<typename T>
        inline void register_output(T& val) {
            get_tape().registerOutput(val);
        }
        
        template<typename T>
        inline void set_gradient(T& val, double grad) {
            val.setGradient(grad);
        }
        
        template<typename T>
        inline double get_gradient(const T& val) {
            return val.getGradient();
        }
        
        inline void evaluate_tape() {
            get_tape().evaluate();
        }
        
        template<typename T>
        inline double to_double(const T& val) {
            return codi::RealTraits::getPassiveValue(val);
        }
    }
    
#else
    // No AD - use plain doubles
    namespace dmc {
        using Real = double;
        
        constexpr bool AD_ENABLED = false;
        
        inline void activate_tape() {}
        inline void deactivate_tape() {}
        inline void reset_tape() {}
        
        template<typename T>
        inline void register_input(T&) {}
        
        template<typename T>
        inline void register_output(T&) {}
        
        template<typename T>
        inline void set_gradient(T&, double) {}
        
        template<typename T>
        inline double get_gradient(const T&) { return 0.0; }
        
        inline void evaluate_tape() {}
        
        template<typename T>
        inline double to_double(const T& val) { return val; }
    }
#endif

#include <cmath>
#include <algorithm>

namespace dmc {
    // Safe math operations that work with both AD and non-AD types
    template<typename T>
    inline T safe_pow(const T& base, double exp) {
        using std::pow;
        return pow(base, exp);
    }
    
    template<typename T>
    inline T safe_sqrt(const T& x) {
        using std::sqrt;
        return sqrt(x);
    }
    
    template<typename T>
    inline T safe_max(const T& a, const T& b) {
        // Smooth max for AD (to avoid discontinuous gradients)
        // return 0.5 * (a + b + sqrt((a - b) * (a - b) + 1e-10));
        // Or simple max if we don't care about smoothness:
        return (a > b) ? a : b;
    }
    
    template<typename T>
    inline T safe_min(const T& a, const T& b) {
        return (a < b) ? a : b;
    }
    
    template<typename T>
    inline T clamp(const T& val, const T& lo, const T& hi) {
        return safe_max(lo, safe_min(val, hi));
    }
}

#endif // DMC_ROUTE_TYPES_HPP
