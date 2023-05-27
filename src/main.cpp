#include "utils.h"
#include "backgroundcosmology.h"
#include "supernovafitting.h"
#include "recombinationhistory.h"
#include "perturbations.h"
#include "powerspectrum.h"



//  Functions to run extra time-demanding stuff:

void m1_MCMC();
void m2_Saha(BackgroundCosmology cosmo);
void m4_decomp(PowerSpectrum power);





int main(int narg, char **argv){

    //=========================================================================
    // Parameters
    //=========================================================================

    // Background parameters (today)
    double h           = 0.67;
    double Omegab      = 0.05;
    double Omegac      = 0.267;
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

    
    // Run which milestones?
    bool milestone1 = false;
    bool milestone2 = false;
    bool milestone3 = false;
    bool milestone4 = true;



    Utils::StartTiming("CMB");


    //  ----------------------
    //  Milestone I
    //  ----------------------

    if(milestone1){

        // Set up and solve the background
        BackgroundCosmology cosmo(h, Omegab, Omegac, OmegaK, Neff, TCMB);
        cosmo.info();
        cosmo.solve(print, 1e5);
        
        // Output background evolution quantities
        cosmo.output("background_cosmology.txt");

        // m1_MCMC();
        
    }


    //  ----------------------
    //  Milestone II
    //  ----------------------

    if(milestone2){

        BackgroundCosmology cosmo(h, Omegab, Omegac, OmegaK, Neff, TCMB);
        cosmo.solve(false, 1e5);

        // Solve the recombination history
        RecombinationHistory rec(&cosmo, Yp);
        rec.info();
        rec.solve(print);

        // Output recombination quantities
        rec.output("recombination.txt");

        // m2_Saha(cosmo);
    }
    


    //  ----------------------
    //  Milestone III
    //  ----------------------

    if(milestone3){

        BackgroundCosmology cosmo2(h, Omegab, Omegac, OmegaK, 0, TCMB);     // Same with Neff=0
        cosmo2.solve(false, 1e5);
        RecombinationHistory rec2(&cosmo2, Yp);
        rec2.solve(false, 1e5, 1e5, 1e5);

        // Solve the perturbations
        Perturbations pert(&cosmo2, &rec2);
        pert.info();
        pert.solve();
        
        // Output perturbation quantities
        double kvalue;
        kvalue = 0.001 / Constants.Mpc;
        pert.output(kvalue, "perturbations_k0.001.txt");
        kvalue = 0.01 / Constants.Mpc;
        pert.output(kvalue, "perturbations_k0.01.txt");
        kvalue = 0.1 / Constants.Mpc;
        pert.output(kvalue, "perturbations_k0.1.txt");


    }

    //  ----------------------
    //  Milestone IV
    //  ----------------------

    if(milestone4){

        BackgroundCosmology cosmo(h, Omegab, Omegac, OmegaK, 0, TCMB);     // Same with Neff=0
        cosmo.solve(false, 1e5);
        RecombinationHistory rec(&cosmo, Yp);
        rec.solve(false, 1e5, 1e5, 1e5);
        Perturbations pert(&cosmo, &rec);
        pert.solve(5e3, 5e3, 100, 100);


        PowerSpectrum power(&cosmo, &rec, &pert, A_s, n_s, kpivot_mpc);
        power.solve();
        power.output("dells.txt");

        std::vector<int> ellvalues{6, 100, 200, 500, 1000};
        power.output("power.txt", ellvalues);



    }

    // TESTING


    // BackgroundCosmology cosmo(h, Omegab, Omegac, OmegaK, 0, TCMB);     // Same with Neff=0
    // // cosmo.info();
    // cosmo.solve(true, 1e4);
    // RecombinationHistory rec(&cosmo, Yp);
    // rec.solve(true, 1e5, 1e5, 1e5);



    // std::cout << "k/Mpc = 0.1:   " << cosmo.find_horizon_entry(0.1, -13, -8) << std::endl;
    // std::cout << "k/Mpc = 0.01:  " << cosmo.find_horizon_entry(0.01, -11, -6) << std::endl;
    // std::cout << "k/Mpc = 0.001: " << cosmo.find_horizon_entry(0.001, -8, -3) << std::endl;



    // // Solve the perturbations
    // Perturbations pert(&cosmo, &rec);
    // // pert.info();
    // pert.solve();



    // PowerSpectrum power(&cosmo, &rec, &pert, A_s, n_s, kpivot_mpc);
    // // power.info();
    // // power.solve();


    std::cout << "-------------------------------------------" << std::endl;
    Utils::EndTiming("CMB");

    return 0;
}
















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


void m4_decomp(PowerSpectrum power){
    // Decompose Cell/Dell into different contributions
    power.solve_decomposed();
    power.output_decomposed("dells_decomp.txt");

}
