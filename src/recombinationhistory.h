#ifndef _RECOMBINATION_HISTORY_HEADER
#define _RECOMBINATION_HISTORY_HEADER
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include "utils.h"
#include "backgroundcosmology.h"

using Vector = std::vector<double>;

class RecombinationHistory{
  private:

    // The cosmology we use
    BackgroundCosmology *cosmo = nullptr;
    
    // Helium fraction
    double Yp;
 
    // The start and end points for recombination arrays (can be modified)
    const double x_start  = -14;
    const double x_end    = 2;
    
    // Numbers of points of Xe,ne array (modify as you see fit)
    const int npts_rec_arrays = 4000;
  
    // Xe for when to switch between Saha and Peebles
    const double Xe_saha_limit = 0.99;

    //===============================================================
    // [1] Computation of Xe (Saha and Peebles equation)
    //===============================================================
 
    // Compute Xe from the Saha equation
    std::pair<double,double> electron_fraction_from_saha_equation(double x) const;
    
    // Right hand side of the dXedx Peebles equation
    int rhs_peebles_ode(double x, const double *y, double *dydx);
    
    // Solve for Xe 
    void solve_number_density_electrons();
    
    //===============================================================
    // [2] Compute tau and visibility functions
    //===============================================================

    // The two things we need to solve: Xe/ne and tau
    void solve_for_optical_depth_tau();

    // Splines contained in this class
    // Spline log_Xe_of_x_spline{"Xe"};
    Spline Xe_of_x_spline{"Xe"};
    Spline tau_of_x_spline{"tau"}; 
    Spline gt_of_x_spline{"g"};  

    // "Helper" variables
    double _8pi = 8*M_PI;             // 8π
    double _8piG = _8pi*Constants.G;  // 8πG
    double _hbhb = Constants.hbar*Constants.hbar;   // hbar^2
    
    double _H0H0;   // H0^2
   
    

  public:

    // Construtors
    RecombinationHistory() = delete;
    RecombinationHistory(
        BackgroundCosmology *cosmo, 
        double Yp=0);

    // Do all the solving
    void solve();
    
    // Print some useful info about the class
    void info() const;

    // Output some data to file
    void output(const std::string filename) const;

    // Get functions that we must implement
    double tau_of_x(double x) const;
    double dtaudx_of_x(double x) const;
    double ddtaudxx_of_x(double x) const;
    double gt_of_x(double x) const;
    double dgtdx_of_x(double x) const;
    double ddgtdxx_of_x(double x) const;
    double Xe_of_x(double x) const;
    double ne_of_x(double x) const;
    double get_Yp() const;

    double nb_of_x(double x) const;
};

#endif