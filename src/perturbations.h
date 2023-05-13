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

    BackgroundCosmology *cosmo = nullptr;   // background
    RecombinationHistory *rec  = nullptr;   // recombination history
   
    // The scales we integrate over
    const int n_k        = 100;             // #k's to consider
    const double k_min   = Constants.k_min; // minimum wavenumber k
    const double k_max   = Constants.k_max; // maximum wavenumber k
    
    // Start and end of the time-integration
    const int n_x        = int(1e4)+1;    // #points in x-array
    const double x_start = -20.;          // start of x-array
    const double x_end   = 0.;            // end of x-array



    // Splines of scalar perturbations quantities
    Spline2D delta_c_spline{"delta_c"};   // δ_c spline
    Spline2D delta_b_spline{"delta_b"};   // δ_b spline
    Spline2D u_c_spline{"u_c"};           // u_c spline
    Spline2D u_b_spline{"u_b"};           // u_b spline
    Spline2D Phi_spline{"Phi"};           // Φ spline
    Spline2D Psi_spline{"Psi"};           // Ψ spline


   
    // Splines of source functions (ST for temperature; SE for polarization)
    Spline2D ST_spline{"ST"};   // St spline

    Spline2D ST_SW_spline{"ST_Sachs_Wolfe"};
    Spline2D ST_ISW_spline{"ST_integrated_Sachs_Wolfe"};
    Spline2D ST_Doppler_spline{"ST_Doppler"};
    Spline2D ST_pol_spline{"ST_polarisation"};  //??????


    
    //  splines of multipole quantities:
    std::vector<Spline2D> Theta_spline;   //  Θ_ℓ spline, ℓ = 0, 1, 2


    // misc:
    Spline dtau{"abs_dtaudx"};
    Spline dtau_over_Hp{"abs_dtaudx_over_Hp"};



    /**
     * @brief Set initial condition at the start (which is in tight coupling).
     * @param x the time at the start
     * @param k the wavenumber k
     * @return y_init for tight coupling
    */
    Vector set_ic(
        const double x, 
        const double k) const;
    


    /**
     * @brief The right-hand side of the ODEs in the tight coupling regime.
     * @param x the time point x = ln(a)
     * @param k the wavenumber k
     * @param y current (tight coupling-relevant) perturbation quantities 
     * @param dydx the r.h.s. to compute
     * @return GSL_SUCCESS
    */
    int rhs_tight_coupling_ode(double x, double k, const double *y, double *dydx);
    

    /**
     * @brief Compute time at which tight coupling ends.
     * @param k the wavenumber k
     * @return the time of end of tight coupling
    */
    double get_tight_coupling_time(const double k) const;
    


    /**
     * @brief Set initial condition for the full system after tight coupling.
     * @param y_tight_coupling the perturbations at the end of tight coupling
     * @param x the time at the end of tight coupling
     * @param k the wavenumber k
     * @return y at end of tight coupling
    */
    Vector set_ic_after_tight_coupling(
        const Vector &y_tight_coupling, 
        const double x, 
        const double k) const;


    /**
     * @brief The right-hand side of the ODEs in the full system, i.e. after tight coupling.
     * @param x the time point x = ln(a)
     * @param k the wavenumber k
     * @param y current perturbation quantities
     * @param dydx the r.h.s. to compute
     * @return GSL_SUCCESS
    */
    int rhs_full_ode(double x, double k, const double *y, double *dydx);
    

    /**
     * @brief Integrate perturbation equations and spline the result.
     * @param nsteps_x number of steps to integrate over
     * @param npts_k number of k-values
    */
    void integrate_perturbations(int nsteps_x=1e4, int npts_k=100);
    

    /**
     * @brief Compute source function St(x,k) and spline result.
     * @param npts_x number of x-values
     * @param npts_k number of k-values
    */
    void compute_source_functions(int npts_x=1e4+1, int npts_k=100);

    /**
     * @brief Expression for gravitational potential Ψ(x,k).
     * @param x the time point x = ln(a)
     * @param k the wavenumber k 
     * @param Phi current potential Φ(x,k)
     * @param Theta2 current quadrupole Θ_2(x,k)
     * @return Ψ(x,k)
    */
    double expr_Psi(double x, double k, double Phi, double Theta2) const;

    /**
     * @brief Expression for quadrupole Θ_2(x,k) in tight coupling regime.
     * @param x the time point x = ln(a)
     * @param k the wavenumber k 
     * @param Theta1 current dipole Θ_1(x,k)
     * @return Θ_2(x,k)
    */
    double expr_Theta2(double x, double k, double Theta1) const;

    /**
     * @brief Expression for quadrupole Θ_ℓ(x,k), ℓ > 2,  in tight coupling regime.
     * @param x the time point x = ln(a)
     * @param k the wavenumber k 
     * @param ell moment ℓ
     * @param Theta_prev current Θ_{ℓ-1}(x,k)
     * @return Θ_ℓ(x,k)
    */
    double expr_Thetaell(double x, double k, int ell, double Thetaell_prev) const;



  public:

    // Constructors:

    Perturbations() = default;
    /**
     * @brief Create an instance of the Perturbations-class with a defined background and recombination history.
     * @param cosmo instance of the Backgroundcosmology-class, solved 
     * @param rec instance of the RecombinationHistory-class, solved
     * @return Perturbations-object
    */
    Perturbations(
        BackgroundCosmology *cosmo, 
        RecombinationHistory *rec); 

    // Methods:

    /**
     * @brief Print some useful info about the class
    */
    void info() const;


    /**
     * @brief Output some results to a file.
     * @param k the wavenumber k
     * @param filename name of file to be stored in output-directory
    */
    void output(const double k, const std::string filename) const;

    // >> solve-methods:

    /**
     * @brief Do all the solving.
     * @param nsteps_x_perturbations number of steps in x-integration of perturbations
     * @param npts_k_perturbations number of points in k-array for which to integrate perturbations
     * @param npts_x_source number of points in x-array for which to compute source function
     * @param npts_k_source number of points in k-array for which to compute source function
    */
    void solve(int nsteps_x_perturbations=1e4, int npts_x_source=1e4+1, int npts_k_perturbations=100, int npts_k_source=100);


    // >> get-methods:


    /**
     * @brief Compute the CDM perturbation δ_c(x,k).
     * @param x the time point x = ln(a)
     * @param k the wavenumber k 
     * @return δ_c(x,k)
    */
    double get_delta_c(const double x, const double k) const;

    /**
     * @brief Compute the baryonic matter perturbation δ_b(x,k).
     * @param x the time point x = ln(a)
     * @param k the wavenumber k 
     * @return δ_b(x,k)
    */
    double get_delta_b(const double x, const double k) const;

    /**
     * @brief Compute the CDM bulk velocity u_c(x,k).
     * @param x the time point x = ln(a)
     * @param k the wavenumber k 
     * @return u_c(x,k)
    */
    double get_u_c(const double x, const double k) const;

    /**
     * @brief Compute the baryonic bulk velocity u_b(x,k).
     * @param x the time point x = ln(a)
     * @param k the wavenumber k 
     * @return u_b(x,k)
    */
    double get_u_b(const double x, const double k) const;

    /**
     * @brief Compute the gravitational potential Φ(x,k).
     * @param x the time point x = ln(a)
     * @param k the wavenumber k 
     * @return Φ(x,k)
    */
    double get_Phi(const double x, const double k) const;

    /**
     * @brief Compute the gravitational potential Ψ(x,k).
     * @param x the time point x = ln(a)
     * @param k the wavenumber k 
     * @return Ψ(x,k)
    */
    double get_Psi(const double x, const double k) const;

    /**
     * @brief Compute the anisotropic stress Π(x,k).
     * @param x the time point x = ln(a)
     * @param k the wavenumber k 
     * @return Π(x,k)
    */
    double get_Pi(const double x, const double k) const;

    /**
     * @brief Compute the photon temperature perturbation Θ_ℓ(x,k).
     * @param x the time point x = ln(a)
     * @param k the wavenumber k 
     * @param ell the moment ℓ 
     * @return Θ_ℓ(x,k)
    */
    double get_Theta(const double x, const double k, const int ell) const;


    /**
     * @brief Compute the photon temperature source function St(x,k).
     * @param x the time point x = ln(a)
     * @param k the wavenumber k 
     * @return St(x,k)
    */
    double get_Source_T(const double x, const double k) const;

    double get_Source_T(const double x, const double k, const int term) const;

};

#endif