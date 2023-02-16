#ifndef _BACKGROUNDCOSMOLOGY_HEADER
#define _BACKGROUNDCOSMOLOGY_HEADER
#include <iostream>
#include <fstream>
#include "utils.h"

using Vector = std::vector<double>;

class BackgroundCosmology{
  private:
   
    // Cosmological parameters
    double h;                       // Little h = H0/(100km/s/Mpc)
    double Omegab0;                  // Baryon density today
    double OmegaCDM0;                // CDM density today
    double OmegaLambda0;             // Dark energy density today
    double Neff;                    // Effective number of relativistic species (3.046 or 0 if ignoring neutrinos)
    double TCMB0;                    // Temperature of the CMB today in Kelvin
   
    // Derived parameters
    double Omegar0;                  // Photon density today (follows from TCMB0)
    double Omeganu0;                 // Neutrino density today (follows from TCMB0 and Neff)
    double Omegak0;                  // Curvature density = 1 - OmegaM - Omegar0 - Omeganu0 - OmegaLambda0
    double H0;                      // The Hubble parameter today H0 = 100h km/s/Mpc

    // Start and end of x-integration (can be changed)
    double x_start = Constants.x_start;
    double x_end   = Constants.x_end;
  

    // Splines to be made
    Spline eta_of_x_spline{"eta"};
    Spline t_of_x_spline{"eta"};
 
  public:

    // Constructors 
    BackgroundCosmology() = delete;
    BackgroundCosmology(
        double h, 
        double Omegab0, 
        double OmegaCDM0, 
        double Omegak0,
        double Neff, 
        double TCMB0
        );

    // Print some useful info about the class
    void info() const;

    // Do all the solving
    void solve();

    // Output some results to file
    void output(const std::string filename) const;

    // Get functions that we must implement
    double eta_of_x(double x) const;
    double H_of_x(double x) const;
    double Hp_of_x(double x) const;
    double dHpdx_of_x(double x) const;
    double ddHpddx_of_x(double x) const;
    double get_Omegab(double x = 0.0) const; 
    double get_Omegam(double x = 0.0) const; 
    double get_Omegar(double x = 0.0) const;
    double get_Omegar_tot(double x = 0.0) const; 
    double get_Omeganu(double x = 0.0) const;
    double get_OmegaCDM(double x = 0.0) const; 
    double get_OmegaLambda(double x = 0.0) const; 
    double get_Omegak(double x = 0.0) const; 
    double get_Omegamnu(double x = 0.0) const; 
    double get_H0() const;
    double get_h() const;
    double get_Neff() const;
    double get_TCMB(double x = 0.0) const;

    double t_of_x(double x) const;
    double r_of_Chi(double Chi) const;
    double Chi_of_x(double x) const;
    double dA_of_x(double x) const;
    double dL_of_x(double x) const;
    void solve_time();

    // Distance measures
    double get_luminosity_distance_of_x(double x) const;
    double get_comoving_distance_of_x(double x) const;


};

#endif