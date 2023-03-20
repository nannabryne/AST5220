#include "utils.h"
#include "backgroundcosmology.h"
#include "supernovafitting.h"
#include "recombinationhistory.h"



void m1_MCMC(){
    Utils::StartTiming("MCMC");
    BackgroundCosmology new_cosmo = mcmc_fit_to_supernova_data("supernovadata.txt", "mcmc_fitting.txt");
    std::cout << "\n\n";
    Utils::EndTiming("MCMC");
    std::cout << "\n\n";

    new_cosmo.info();
    new_cosmo.solve(1e5, false);

    new_cosmo.output("revised_background_cosmology.txt");
}


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
    double Yp          = 0.0;

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
    cosmo.solve(1e5, false);
    
    
    // Output background evolution quantities
    cosmo.output("background_cosmology.txt");

    // m1_MCMC();

    // if(narg>1){
    //     // MCMC analysis
    //     if(argv[1]=="MCMC")
    //         m1_MCMC();
    // }
    
    // //  ----------------------
    // //  Milestone II
    // //  ----------------------


    // // Solve the recombination history
    // RecombinationHistory rec(&cosmo, Yp);
    // rec.solve();
    // rec.info();

    // // Output recombination quantities
    // rec.output("recombination.txt");


    return 0;
}









