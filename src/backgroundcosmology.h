#ifndef _BACKGROUNDCOSMOLOGY_HEADER
#define _BACKGROUNDCOSMOLOGY_HEADER
#include <iostream>
#include <fstream>
#include "utils.h"

using Vector = std::vector<double>;

class BackgroundCosmology{
  private:
   
    // Cosmological parameters:
    double h0;              // little h0 = H0/(100km/s/Mpc)
    double Omegab0;         // baryon density today
    double Omegac0;         // CDM density today
    double OmegaLambda0;    // dark energy density today
    double Neff;            // effective number of relativistic species (3.046 or 0 if ignoring neutrinos)
    double TCMB0;           // temperature of the CMB today in Kelvin
   
    // Derived parameters:
    double Omegagamma0;     // photon density today (follows from TCMB0)
    double Omeganu0;        // neutrino density today (follows from TCMB0 and Neff)
    double OmegaK0;         // curvature density = 1 - OmegaM0 - OmegaR0 - OmegaLambda0
    double H0;              // the Hubble parameter today H0 = 100h km/s/Mpc
    double OmegaM0;         // total matter density today = Omegab0 + Omegac0
    double OmegaR0;         // total radiation density today = Omegagamma0 + Omeganu0

    // Start and end of x-integration:
    double x_start = Constants.x_start;
    double x_end   = Constants.x_end;
  
    // Splines to be made:
    Spline eta_of_x_spline{"eta"};  // conformal time spline
    Spline t_of_x_spline{"t"};      // cosmic time spline

    // Helper functions:
    /* Ξ_0(x) */
    double Xi0(double x) const; 
    /* Ξ_1(x) */
    double Xi1(double x) const;
    /* Ξ_2(x) */
    double Xi2(double x) const;


    /**
     * @brief Locate milestones and print result.
     * @param nsteps number of steps in x-array
    */
    void milestones(int nsteps);
 
  public:

    // Constructors:

    BackgroundCosmology() = delete;
    /**
     * @brief Create an instance of the Backgroundcosmology-class with a set of parameters describing the Universe today. Compute the remaining parameters based on the input parameters.
     * @param h0 Hubble parameter today (H0) divided by 100 km/s / Mpc
     * @param Omegab0 baryon density today
     * @param Omegac0 cold dark matter density today
     * @param OmegaK0 curvature density today ( = 1 - Σ(Ω_i0) )
     * @param Neff effective number of relativistic species (3.046 or 0 if ignoring neutrinos)
     * @param TCMB0 CMB temperature today [K]
     * @return Backgroundcosmology-object
    */
    BackgroundCosmology(
        double h0, 
        double Omegab0, 
        double Omegac0, 
        double OmegaK0,
        double Neff, 
        double TCMB0
        );
      
    // Methods:

    /**
     * @brief Print some useful information about the class object.
    */
    void info() const;

    /**
     * @brief Output some results to a file.
     * @param filename name of file to be stored in output-directory
    */
    void output(const std::string filename) const;


    // >> solve-methods:

    /**
     * @brief Solve the background by finding eta(x).
     * @param nsteps number of steps in x-array
    */
    void solve_conformal_time(int nsteps=1e4);

    /**
     * @brief Solve the background by finding t(x).
     * @param nsteps number of steps in x-array
    */
    void solve_cosmic_time(int nsteps=1e4);

    /**
     * @brief Solve the background.
     * @param print_milestones whether to print the table showing the time positions of the milestones
     * @param nsteps number of steps in x-array
    */
    void solve(bool print_milestones=true, int nsteps=1e4);


    // >> get-methods:

    /**
     * @brief Compute the redshift z(x).
     * @param x the time point x = ln(a)
     * @return z(x)
    */
    double z_of_x(double x) const;

    /**
     * @brief Compute the conformal time η(x).
     * @param x the time point x = ln(a)
     * @return η(x)
    */
    double eta_of_x(double x) const;

    /**
     * @brief Compute the cosmic time t(x).
     * @param x the time point x = ln(a)
     * @return t(x)
    */
    double t_of_x(double x) const;

    /**
     * @brief Compute the photon radius r(χ).
     * @param Chi the comoving distance χ = η(x=0) - n(x)
     * @return r(χ)
    */
    double r_of_Chi(double Chi) const;

    /**
     * @brief Compute the comoving distance χ(x).
     * @param x the time point x = ln(a)
     * @return χ(x)
    */
    double Chi_of_x(double x) const;

    /**
     * @brief Compute the angular distance d_A(x).
     * @param x the time point x = ln(a)
     * @return d_A(x)
    */
    double dA_of_x(double x) const;

    /**
     * @brief Compute the luminosity distance d_L(x).
     * @param x the time point x = ln(a)
     * @return d_L(x)
    */
    double dL_of_x(double x) const;

    /**
     * @brief Compute the Hubble factor H(x).
     * @param x the time point x = ln(a)
     * @return H(x)
    */
    double H_of_x(double x) const;

    /**
     * @brief Compute the conformal Hubble factor Hp(x).
     * @param x the time point x = ln(a)
     * @return Hp(x)
    */
    double Hp_of_x(double x) const;

    /**
     * @brief Compute the total time derivative of the conformal Hubble factor Hp'(x).
     * @param x the time point x = ln(a)
     * @return d/dx[Hp(x)]
    */
    double dHpdx_of_x(double x) const;

    /**
     * @brief Compute the double total time derivative of the conformal Hubble factor Hp''(x).
     * @param x the time point x = ln(a)
     * @return d^2/dx^2[Hp(x)]
    */
    double ddHpdxx_of_x(double x) const;

    /**
     * @brief Compute the baryon density Ω_b(x).
     * @param x the time point x = ln(a)
     * @return Ω_b(x)
    */
    double get_Omegab(double x = 0.0) const; 

    /**
     * @brief Compute the total matter density Ω_m(x).
     * @param x the time point x = ln(a)
     * @return Ω_m(x)
    */
    double get_OmegaM(double x = 0.0) const; 

    /**
     * @brief Compute the photon density Ω_γ(x).
     * @param x the time point x = ln(a)
     * @return Ω_γ(x)
    */
    double get_Omegagamma(double x = 0.0) const;

    /**
     * @brief Compute the total radiation density Ω_r(x).
     * @param x the time point x = ln(a)
     * @return Ω_r(x)
    */
    double get_OmegaR(double x = 0.0) const; 

    /**
     * @brief Compute the neutrino density Ω_ν(x).
     * @param x the time point x = ln(a)
     * @return Ω_ν(x)
    */
    double get_Omeganu(double x = 0.0) const;

    /**
     * @brief Compute the cold dark matter density Ω_c(x).
     * @param x the time point x = ln(a)
     * @return Ω_c(x)
    */
    double get_Omegac(double x = 0.0) const; 

    /**
     * @brief Compute the dark energy density Ω_Λ(x).
     * @param x the time point x = ln(a)
     * @return Ω_Λ(x)
    */
    double get_OmegaLambda(double x = 0.0) const; 

    /**
     * @brief Compute the curvature density Ω_k(x).
     * @param x the time point x = ln(a)
     * @return Ω_k(x)
    */
    double get_OmegaK(double x = 0.0) const; 

    
    /**
     * @brief Get the Hubble parameter today H(x=0) = H_0.
     * @return H_0
    */
    double get_H0() const;

    /**
     * @brief Get the little Hubble parameter today h = H(x=0) / (100km/s / Mpc).
     * @return h
    */
    double get_h() const;

    /**
     * @brief Get the effective neutrino number N_eff.
     * @return N_eff
    */
    double get_Neff() const;

    /**
     * @brief Compute the CMB temperature T_CMB(x).
     * @param x the time point x = ln(a)
     * @return T_CMB(x)
    */
    double get_TCMB(double x = 0.0) const;

    // >> distance measures (aliases):

    /**
     * @brief Compute the luminosity distance d_L(x).
     * @param x the time point x = ln(a)
     * @return d_L(x)
    */
    double get_luminosity_distance_of_x(double x) const;
    
    /**
     * @brief Compute the comoving distance χ(x).
     * @param x the time point x = ln(a)
     * @return χ(x)
    */
    double get_comoving_distance_of_x(double x) const;


    


};

#endif