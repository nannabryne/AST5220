#ifndef _PERTURBATIONS_HEADER
#define _PERTURBATIONS_HEADER
#ifdef _USEOPENMP
#include <omp.h>
#endif
#include <vector>
#include <fstream>
#include <algorithm>
#include "utils.h"
#include "backgroundcosmology.h"
#include "recombinationhistory.h"

using Vector   = std::vector<double>;
using Vector2D = std::vector<Vector>;

class Perturbations{
  private:

    BackgroundCosmology *cosmo = nullptr;
    RecombinationHistory *rec  = nullptr;
   
    // The scales we integrate over
    const int n_k        = 100;
    const double k_min   = Constants.k_min;
    const double k_max   = Constants.k_max;
    
    // Start and end of the time-integration
    const int n_x        = 10000+1;
    const double x_start = -18.;
    const double x_end   = 0.;

    // Below is a full list of splines you probably need, 
    // but you only need to make the splines you will need

    // Splines of scalar perturbations quantities
    Spline2D delta_c_spline{"delta_c"};
    Spline2D delta_b_spline{"delta_b"};
    Spline2D u_c_spline{"u_c"};
    Spline2D u_b_spline{"u_b"};
    Spline2D Phi_spline{"Phi"};
    Spline2D Pi_spline{"Pi"};
    Spline2D Psi_spline{"Psi"};
   
    // Splines of source functions (ST for temperature; SE for polarization)
    Spline2D ST_spline{"ST"};
    Spline2D SE_spline{"SE"};
    
    // Splines of mulipole quantities
    // NB: If you use there you have to allocate the container first
    // e.g. Theta_spline = std::vector<Spline2D>(n_ell_Theta); before using it
    std::vector<Spline2D> Theta_spline;
    std::vector<Spline2D> Theta_p_spline;
    std::vector<Spline2D> Nu_spline;
    
    //==========================================================
    // [1] Tight coupling ODE system
    //==========================================================

    // Set the initial conditions at the start (which is in tight coupling)
    Vector set_ic(
        const double x, 
        const double k) const;
    
    // Right hand side of the ODE in the tight coupling regime
    int rhs_tight_coupling_ode(double x, double k, const double *y, double *dydx);
    
    // Compute the time when tight coupling ends
    double get_tight_coupling_time(const double k) const;
    
    //==========================================================
    // [2] The full ODE system 
    //==========================================================
    
    // Set initial condition after tight coupling
    Vector set_ic_after_tight_coupling(
        const Vector &y_tight_coupling, 
        const double x, 
        const double k) const;

    // Right hand side of the ODE in the full regime
    int rhs_full_ode(double x, double k, const double *y, double *dydx);
    
    //==========================================================
    // [3] Integrate the full system
    //==========================================================
    
    // Integrate perturbations and spline the result
    void integrate_perturbations();
    
    //==========================================================
    // [4] Compute source functions from the result
    //==========================================================
    
    // Compute source functions and spline the result
    void compute_source_functions();


    double expr_Psi(double x, double k, double Phi, double Theta2) const;
    double expr_Theta2(double x, double k, double Theta1) const;

  public:

    // Constructors
    Perturbations() = default;
    Perturbations(
        BackgroundCosmology *cosmo, 
        RecombinationHistory *rec); 

    // Do all the solving
    void solve();
    
    // Print some useful info about the class
    void info() const;

    // Output info to file
    void output(const double k, const std::string filename) const;

    // Get the quantities we have integrated

    /**
     * @brief Compute the CDM perturbation δ_c(x,k).
     * @param x the time point x = ln(a)
     * @param k the wave number k 
     * @return δ_c(x,k)
    */
    double get_delta_c(const double x, const double k) const;

    /**
     * @brief Compute the baryonic matter perturbation δ_b(x,k).
     * @param x the time point x = ln(a)
     * @param k the wave number k 
     * @return δ_b(x,k)
    */
    double get_delta_b(const double x, const double k) const;

    /**
     * @brief Compute the CDM bulk velocity u_c(x,k).
     * @param x the time point x = ln(a)
     * @param k the wave number k 
     * @return u_c(x,k)
    */
    double get_u_c(const double x, const double k) const;

    /**
     * @brief Compute the baryonic bulk velocity u_b(x,k).
     * @param x the time point x = ln(a)
     * @param k the wave number k 
     * @return u_b(x,k)
    */
    double get_u_b(const double x, const double k) const;

    /**
     * @brief Compute the gravitational potential Φ(x,k).
     * @param x the time point x = ln(a)
     * @param k the wave number k 
     * @return Φ(x,k)
    */
    double get_Phi(const double x, const double k) const;

    /**
     * @brief Compute the XXXXXXX Ψ(x,k).
     * @param x the time point x = ln(a)
     * @param k the wave number k 
     * @return Ψ(x,k)
    */
    double get_Psi(const double x, const double k) const;


    double get_Pi(const double x, const double k) const;

    /**
     * @brief Compute the photon temperature perturbation Θ_ℓ(x,k).
     * @param x the time point x = ln(a)
     * @param k the wave number k 
     * @return Θ_ℓ(x,k)
    */
    double get_Theta(const double x, const double k, const int ell) const;
    double get_Theta_p(const double x, const double k, const int ell) const;
    double get_Nu(const double x, const double k, const int ell) const;
    double get_Source_T(const double x, const double k) const;
    double get_Source_E(const double x, const double k) const;
};

#endif