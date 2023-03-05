#include "utils.h"
#include "backgroundcosmology.h"
#include "supernovafitting.h"



/* Temporary main-function */
int main(int narg, char **argv){

    //=========================================================================
    // Parameters
    //=========================================================================

    // Background parameters (today)
    double h           = 0.67;
    double Omegab      = 0.05;
    double OmegaCDM    = 0.267;
    double OmegaK      = 0.0;
    double Neff        = 3.046;
    double TCMB        = 2.7255;

    // Recombination parameters
    double Yp          = 0.245;

    // Power-spectrum parameters
    double A_s         = 2.1e-9;
    double n_s         = 0.965;
    double kpivot_mpc  = 0.05;


    //  ----------------------
    //  Milestone I
    //  ----------------------

    // Set up and solve the background

    BackgroundCosmology cosmo(h, Omegab, OmegaCDM, OmegaK, Neff, TCMB);
    cosmo.info();
    cosmo.solve(1e5);
    
    
    // Output background evolution quantities
    cosmo.output("background_cosmology.txt");

    // MCMC analysis

    Utils::StartTiming("MCMC");
    mcmc_fit_to_supernova_data("supernovadata.txt", "mcmc_fitting.txt");
    std::cout << "\n\n";
    Utils::EndTiming("MCMC");
    std::cout << "\n\n";

    // Minimum chi^2 found 29.2803 0.701725 0.25613 0.0765251
    // Elapsed time for [MCMC]: 215.309 sec


    return 0;
}









