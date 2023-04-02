#ifndef _RECOMBINATION_HISTORY_HEADER
#define _RECOMBINATION_HISTORY_HEADER
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include "utils.h"
#include "backgroundcosmology.h"

using Vector = std::vector<double>;


/*
gt := gtilde
*/

class RecombinationHistory{
  private:

    // The cosmology we use:
    BackgroundCosmology *cosmo = nullptr;
    
    // Helium fraction 
    double Y_P;  // fractional abundance of helium (primordial value)
 
    // The start and end points for recombination arrays
    const double x_start  = -12.;
    const double x_end    = 0.;
  
    double Xe_Saha_limit = 0.99;  // X_e for when to switch between Saha and Peebles



    //===============================================================
    // [1] Computation of Xe (Saha and Peebles equation)
    //===============================================================
    
    /**
     * @brief Compute X_e(x) from the Saha equation.
     * @param x the time point x = ln(a)
     * @return { X_e(x), n_e(x) }
    */
    std::pair<double,double> electron_fraction_from_Saha_equation(double x) const;
    
    /**
     * @brief Compute dX_e/dx from the Peebles equation.
     * @param x the time point x = ln(a)
     * @param Xe current X_e
     * @param dXedx the rhs (dX_e/dx) to compute
     * @return GSL_SUCCESS
    */
    int rhs_Peebles_ode(double x, const double *Xe, double *dXedx);
    
    /**
     * @brief Solve recombination by finding X_e(x) and n_e(x).
     * @param nsteps number of steps in x-array
    */
    void solve_number_density_electrons(int nsteps=80000);
    

    //===============================================================
    // [2] Compute tau and visibility functions
    //===============================================================

    /**
     * @brief Solve ODE to find τ(x) and gt(x).
     * @param nsteps number of steps in integration
    */
    void solve_for_optical_depth_tau(int nsteps=80000);


    //===============================================================
    // [3] Compute sound horizon
    //===============================================================


    /**
     * @brief Solve ODE to find r_s(x).
     * @param nsteps number of steps in integration
    */
    void solve_for_sound_horizon(int nsteps=80000);


    //===============================================================
    // ...
    //===============================================================


    // splines contained in this class:
    Spline log_Xe_of_x_spline{"Xe"};  // electron fraction spline
    Spline tau_of_x_spline{"tau"};    // optical depth spline
    Spline gt_of_x_spline{"gt"};      // visibilty function spline
    Spline rs_of_x_spline{"rs"};      // sound horizon spline

    // "Helper" variables:
    double _8pi = 8*M_PI;                         // 8π
    double _8piG = _8pi*Constants.G;              // 8πG
    double _hbhb = Constants.hbar*Constants.hbar; // hbar^2

    double _H0H0;                                 // H0^2

    /**
     * @brief Locate milestones and print result.
     * @param nsteps number of steps in x-array
    */
    void milestones(int nsteps);
   
    

  public:

    // Construtors:

    RecombinationHistory() = delete;
    /**
     * @brief Create an instance of the RecombinationHistory-class with a defined background and Helium fraction.
     * @param cosmo instance of the Backgroundcosmology-class with a set of parameters describing the Universe today
     * @param Y_P primordial value of the fractional Helium abundance
     * @return RecombinationHistory-object
    */
    RecombinationHistory(
        BackgroundCosmology *cosmo, 
        double Y_P=0);

    /**
     * @brief Change limit for Saha regime.
     * @param Xe_Saha_lower_limit the lowest X_e for which to use Saha equation
    */
    void set_Saha_limit(double Xe_Saha_lower_limit);


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
     * @brief Solve recombination. 
     * @param print_milestones whether to print the table showing the time positions of the milestones
     * @param nsteps_Xe number of steps in x-array for computation of X_e(x) and n_e(x)
     * @param nsteps_tau number of steps in x-array for computation of τ(x) and gt(x)
     * @param nsteps_rs number of steps in x-array for computation of r_s(x)
    */
    void solve(bool print_milestones=false, int nsteps_Xe=8e5, int nsteps_tau=8e5, int nsteps_rs=8e5);
    
    
    // >> get-methods:

    /**
     * @brief Compute the optical depth τ(x).
     * @param x the time point x = ln(a)
     * @return τ(x)
    */
    double tau_of_x(double x) const;

    /**
     * @brief Compute the time derivative of the optical depth dτ(x)/dx.
     * @param x the time point x = ln(a)
     * @return d/dx[τ(x)]
    */
    double dtaudx_of_x(double x) const;

    /**
     * @brief Compute the double time derivative of the optical depth d^2τ(x)/dx^2.
     * @param x the time point x = ln(a)
     * @return d^2/dx^2[τ(x)]
    */
    double ddtaudxx_of_x(double x) const;

    /**
     * @brief Compute the visibility function gt(x).
     * @param x the time point x = ln(a)
     * @return gt(x)
    */
    double gt_of_x(double x) const;

    /**
     * @brief Compute the time derivative of the visibility function dgt(x)/dx.
     * @param x the time point x = ln(a)
     * @return d/dx[gt(x)]
    */
    double dgtdx_of_x(double x) const;

    /**
     * @brief Compute the double time derivative of the visibility function d^2gt(x)/dx^2.
     * @param x the time point x = ln(a)
     * @return d^2/dx^2[gt(x)]
    */
    double ddgtdxx_of_x(double x) const;

    /**
     * @brief Compute the electron fraction X_e(x).
     * @param x the time point x = ln(a)
     * @return X_e(x)
    */
    double Xe_of_x(double x) const;

    /**
     * @brief Compute the electron number density n_e(x).
     * @param x the time point x = ln(a)
     * @return n_e(x)
    */
    double ne_of_x(double x) const;
    
    /**
     * @brief Compute the baryon number density n_b(x).
     * @param x the time point x = ln(a)
     * @return n_b(x)
    */
    double nb_of_x(double x) const;

    /**
     * @brief Compute the sound of speed in the photon-baryon plasma c_s(x). 
     * @param x the time point x = ln(a)
     * @return c_s(x)
    */
    double cs_of_x(double x) const;

    /**
     * @brief Compute the sound horizon s(x). 
     * @param x the time point x = ln(a)
     * @return s(x)
    */
    double rs_of_x(double x) const;


    /**
     * @brief Get fractional abundance of Helium.
     * @return Y_P
    */
    double get_Y_P() const;
    
};

#endif