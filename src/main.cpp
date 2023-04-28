#include "utils.h"
#include "backgroundcosmology.h"
#include "supernovafitting.h"
#include "recombinationhistory.h"
#include "perturbations.h"



void m1_MCMC(){
    Utils::StartTiming("MCMC");
    BackgroundCosmology new_cosmo = mcmc_fit_to_supernova_data("supernovadata.txt", "mcmc_fitting.txt");
    std::cout << "\n\n";
    Utils::EndTiming("MCMC");
    std::cout << "\n\n";

    new_cosmo.info();
    new_cosmo.solve(false, 1e5);

    new_cosmo.output("revised_background_cosmology.txt");
}


void m2_Saha(BackgroundCosmology cosmo){
    // Compare with solution using only Saha
    RecombinationHistory rec_noPee(&cosmo, 0);
    rec_noPee.set_Saha_limit(0.0);
    rec_noPee.solve(false);

    // Output recombination quantities
    rec_noPee.output("recombination_Saha_only.txt");
}


/* Temporary main-function */
int main(int narg, char **argv){

    //=========================================================================
    // Parameters
    //=========================================================================

    // Background parameters (today)
    double h           = 0.67;
    double Omegab      = 0.05;
    double Omegac    = 0.267;
    double OmegaK      = 0.0;
    double Neff        = 3.046;
    double TCMB        = 2.7255;

    // Recombination parameters
    double Yp          = 0.0;

    // Power-spectrum parameters
    double A_s         = 2.1e-9;
    double n_s         = 0.965;
    double kpivot_mpc  = 0.05;

    // Print tables?
    bool print = false;


    //  ----------------------
    //  Milestone I
    //  ----------------------

    // Set up and solve the background

    // BackgroundCosmology cosmo(h, Omegab, Omegac, OmegaK, Neff, TCMB);
    // cosmo.info();
    // cosmo.solve(print, 1e5);
    
    
    // // Output background evolution quantities
    // cosmo.output("background_cosmology.txt");

    // // m1_MCMC();
    

    // //  ----------------------
    // //  Milestone II
    // //  ----------------------


    // // Solve the recombination history
    // RecombinationHistory rec(&cosmo, Yp);
    // rec.info();
    // rec.solve(print);

    // // Output recombination quantities
    // rec.output("recombination.txt");

    // m2_Saha(cosmo);


    //  ----------------------
    //  Milestone III
    //  ----------------------

    BackgroundCosmology cosmo2(h, Omegab, Omegac, OmegaK, 0, TCMB);     // Same with Neff=0
    cosmo2.solve(false, 1e5);
    RecombinationHistory rec2(&cosmo2, Yp);
    rec2.solve(false);

    // Solve the perturbations
    Perturbations pert(&cosmo2, &rec2);
    pert.solve();
    pert.info();
    
    // Output perturbation quantities
    double kvalue;
    kvalue = 0.001 / Constants.Mpc;
    pert.output(kvalue, "perturbations_k0.001.txt");
    kvalue = 0.01 / Constants.Mpc;
    pert.output(kvalue, "perturbations_k0.01.txt");
    kvalue = 0.1 / Constants.Mpc;
    pert.output(kvalue, "perturbations_k0.1.txt");


    /*
    Sett Neff = 0, så får du de samme resultatene som Hans!
    */
        


    return 0;
}









