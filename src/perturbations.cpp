#include "perturbations.h"

//  ----------------------------------
//  CONSTRUCT
//  ----------------------------------

Perturbations::Perturbations(
    BackgroundCosmology *cosmo, 
    RecombinationHistory *rec) : 
  cosmo(cosmo), 
  rec(rec)
{}



//  ----------------------------------
//  SOLVE
//  ----------------------------------


void Perturbations::solve(int nsteps_x_perturbations, int npts_x_source, int npts_k_perturbations, int npts_k_source){

  //  integrate all the perturbation equation and spline the result:
  Utils::StartTiming("integrateperturbation");
  integrate_perturbations(nsteps_x_perturbations, npts_k_perturbations);
  Utils::EndTiming("integrateperturbation");

  //  compute source functions and spline the result:
  Utils::StartTiming("source");
  compute_source_functions(npts_x_source, npts_k_source);
  Utils::EndTiming("source");
}


void Perturbations::integrate_perturbations(int nsteps_x, int npts_k){


  const gsl_odeiv2_step_type * stepper = gsl_odeiv2_step_rk4;
  const double hstart = 1e-6;
  const double abserr = 1e-8;
  const double relerr = 1e-8;

  //  set up k-array for the k's we are to integrate over:
  const int n_k = npts_k;
  Vector k_array = exp(Utils::linspace(log(k_min), log(k_max), n_k));
  
  const int n_x = nsteps_x + 1;
  Vector x_array = Utils::linspace(x_start, x_end, n_x);

  // // create helper splines:
  // Vector dtau(n_x);
  // Vector dtau_H(n_x);
  // double x;
  // for(int i; i<n_x; i++){
  //   x = x_array[i]
  //   dtau[i] = - rec->tau_of_x(x);
  //   dtau_H[i] = dtau[i]/cosmo->Hp_of_x(x);


  // }

  const int nn = n_x*n_k;

  Vector delta_c_array(nn);
  Vector delta_b_array(nn);
  Vector u_c_array(nn);
  Vector u_b_array(nn);
  Vector Phi_array(nn);
  Vector Theta0_array(nn);
  Vector Theta1_array(nn);
  Vector Theta2_array(nn);
  Vector Theta3_array(nn);
  Vector Psi_array(nn);


  //  loop over all wavenumbers:   
  #pragma omp parallel for schedule(dynamic, 1)
  for(int ik=0; ik<n_k; ik++){

    // //  progress bar... (remove to improve performance)
    // if( (10*ik) / n_k != (10*ik+10) / n_k ) {
    //   std::cout << (100*ik+100)/n_k << "% " << std::flush;
    //   if(ik == n_k-1) std::cout << std::endl;
    // }

    double k = k_array[ik];   // current value of k

    //  Find value to integrate to

    double x_end_tight = get_tight_coupling_time(k);
    int n_x_tight = int(n_x*(x_end_tight-x_start)/(x_end-x_start));
    int idx_end_tight = n_x_tight - 1;
    int n_x_full = n_x-n_x_tight + 1;

    Vector x_array_tc = Vector(x_array.begin(), x_array.end()-(n_x_full-1));
    Vector x_array_full = Vector(x_array.begin()+idx_end_tight, x_array.end());
    x_end_tight = x_array_full[0];

  

    //  (1) TIGHT COUPLING INTEGRATION
    

    //  set up initial conditions in the tight coupling regime:
    auto y_tight_coupling_ini = set_ic(x_start, k);

    //  the tight coupling ODE system:
    ODEFunction dydx_tight_coupling = [&](double x, const double *y, double *dydx){
      return rhs_tight_coupling_ode(x, k, y, dydx);
    };

   
    //  integrate from x_start -> x_end_tight:
    ODESolver ode_tight;
    ode_tight.set_accuracy(hstart, abserr, relerr);
    ode_tight.solve(dydx_tight_coupling, x_array_tc, y_tight_coupling_ini, stepper);
    Vector y_end_tight = ode_tight.get_final_data();
    Vector2D y_sol_tc = ode_tight.get_data();



    //  (2) FULL SYSTEM INTEGRATION

    //  set up initial conditions (y_end_tight is the solution at the end of tight coupling):
    auto y_full_ini = set_ic_after_tight_coupling(y_end_tight, x_end_tight, k);

    //  the full ODE system:
    ODEFunction dydx_full = [&](double x, const double *y, double *dydx){
      return rhs_full_ode(x, k, y, dydx);
    };

    //  integrate from x_end_tight -> x_end:
    ODESolver ode_full;
    ode_full.set_accuracy(hstart, abserr, relerr);
    ode_full.solve(dydx_full, x_array_full, y_full_ini, stepper);
    Vector2D y_sol_full = ode_full.get_data();
    

    //  Store the data 


    int d = n_x*ik;


    // Vector2D y_sol_tc2 = ode_tight.get_data_transpose();
    // Vector2D y_sol_full2 = ode_full.get_data_transpose();

    // int d_tc = d;
    // int d_full = d_tc + idx_end_tight;

    // std::copy(y_sol_tc2[Constants.ind_deltac_tc].begin(), y_sol_tc2[Constants.ind_deltac_tc].end()-1, delta_c_array.begin()+d_tc);
    // std::copy(y_sol_full2[Constants.ind_deltac].begin(),  y_sol_full2[Constants.ind_deltac].end(),    delta_c_array.begin()+d_full);

    // std::copy(y_sol_tc2[Constants.ind_deltab_tc].begin(), y_sol_tc2[Constants.ind_deltab_tc].end()-1, delta_b_array.begin()+d_tc);
    // std::copy(y_sol_full2[Constants.ind_deltab].begin(),  y_sol_full2[Constants.ind_deltab].end(),    delta_b_array.begin()+d_full);

    // std::copy(y_sol_tc2[Constants.ind_uc_tc].begin(), y_sol_tc2[Constants.ind_uc_tc].end()-1, u_c_array.begin()+d_tc);
    // std::copy(y_sol_full2[Constants.ind_uc].begin(),  y_sol_full2[Constants.ind_uc].end(),    u_c_array.begin()+d_full);

    // std::copy(y_sol_tc2[Constants.ind_ub_tc].begin(), y_sol_tc2[Constants.ind_ub_tc].end()-1, u_b_array.begin()+d_tc);
    // std::copy(y_sol_full2[Constants.ind_ub].begin(),  y_sol_full2[Constants.ind_ub].end(),    u_b_array.begin()+d_full);



    // std::copy(y_sol_tc2[Constants.ind_Phi_tc].begin(), y_sol_tc2[Constants.ind_Phi_tc].end()-1, Phi_array.begin()+d_tc);
    // std::copy(y_sol_full2[Constants.ind_Phi].begin(),  y_sol_full2[Constants.ind_Phi].end(),    Phi_array.begin()+d_full);

    // std::copy(y_sol_tc2[Constants.ind_start_Theta_tc].begin(), y_sol_tc2[Constants.ind_start_Theta_tc].end()-1, Theta0_array.begin()+d_tc);
    // std::copy(y_sol_full2[Constants.ind_start_Theta].begin(),  y_sol_full2[Constants.ind_start_Theta].end(),    Theta0_array.begin()+d_full);

    // std::copy(y_sol_tc2[Constants.ind_start_Theta_tc+1].begin(), y_sol_tc2[Constants.ind_start_Theta_tc+1].end()-1, Theta1_array.begin()+d_tc);
    // std::copy(y_sol_full2[Constants.ind_start_Theta+1].begin(),  y_sol_full2[Constants.ind_start_Theta+1].end(),    Theta1_array.begin()+d_full);

    double x;
    int i, j;
    Vector y_sol_curr;


    //  during tight coupling:
    j = 0;
    for(int ix=0; ix<idx_end_tight; ix++){
      i = ix + d;
      x = x_array[ix];
      y_sol_curr = y_sol_tc[j];

      delta_c_array[i]    = y_sol_curr[Constants.ind_deltac_tc];
      delta_b_array[i]    = y_sol_curr[Constants.ind_deltab_tc];
      u_c_array[i]        = y_sol_curr[Constants.ind_uc_tc];
      u_b_array[i]        = y_sol_curr[Constants.ind_ub_tc];
      Theta0_array[i]     = y_sol_curr[Constants.ind_start_Theta_tc];
      Theta1_array[i]     = y_sol_curr[Constants.ind_start_Theta_tc+1];
      Theta2_array[i]     = expr_Theta2(x, k, Theta1_array[i]);
      Theta3_array[i]     = expr_Thetaell(x, k, 3, Theta2_array[i]);
      Phi_array[i]        = y_sol_curr[Constants.ind_Phi_tc];
      Psi_array[i]        = expr_Psi(x, k, Phi_array[i], Theta2_array[i]);

      j += 1;
    }
    //  after tight coupling:
    j = 0;
    for(int ix=idx_end_tight; ix<n_x; ix++){
      i = ix + d;
      // j = ix - n_x_tight;
      x = x_array[ix];
      y_sol_curr = y_sol_full[j];

      delta_c_array[i]    = y_sol_curr[Constants.ind_deltac];
      delta_b_array[i]    = y_sol_curr[Constants.ind_deltab];
      u_c_array[i]        = y_sol_curr[Constants.ind_uc];
      u_b_array[i]        = y_sol_curr[Constants.ind_ub];
      Theta0_array[i]     = y_sol_curr[Constants.ind_start_Theta];
      Theta1_array[i]     = y_sol_curr[Constants.ind_start_Theta+1];
      Theta2_array[i]     = y_sol_curr[Constants.ind_start_Theta+2];
      Theta3_array[i]     = y_sol_curr[Constants.ind_start_Theta+3];
      Phi_array[i]        = y_sol_curr[Constants.ind_Phi];
      Psi_array[i]        = expr_Psi(x, k, Phi_array[i], Theta2_array[i]);

      j += 1;
    }
    


  }
  

  //  make all the splines we need

  delta_c_spline.create(x_array, k_array, delta_c_array, "delta_c");
  delta_b_spline.create(x_array, k_array, delta_b_array, "delta_b");
  u_c_spline.create(x_array, k_array, u_c_array, "u_c");
  u_b_spline.create(x_array, k_array, u_b_array, "u_b");

  Phi_spline.create(x_array, k_array, Phi_array, "Phi");
  Psi_spline.create(x_array, k_array, Psi_array, "Psi");

  Theta_spline = std::vector<Spline2D>(4);
  Theta_spline[0].create(x_array, k_array, Theta0_array, "Theta0");
  Theta_spline[1].create(x_array, k_array, Theta1_array, "Theta1");
  Theta_spline[2].create(x_array, k_array, Theta2_array, "Theta2");
  Theta_spline[3].create(x_array, k_array, Theta3_array, "Theta3");
  

}




//  ----------------------------------
//  INITIAL CONDITIONS ++
//  ----------------------------------


Vector Perturbations::set_ic(const double x, const double k) const{

  Vector y_tc(Constants.n_ell_tot_tc);  // the vector we are going to fill
  
  // For integration of perturbations in tight coupling regime (Only 2 photon multipoles + neutrinos needed)


  const int n_ell_Theta_tc      = Constants.n_ell_Theta_tc;
  const int n_ell_neutrinos_tc  = Constants.n_ell_neutrinos_tc;
  const int n_ell_tot_tc        = Constants.n_ell_tot_tc;
  const bool polarization       = Constants.polarization;
  const bool neutrinos          = Constants.neutrinos;

  //  reference to the tight coupling quantities:
  double &delta_c      =  y_tc[Constants.ind_deltac_tc];        // δ_c(x,k) 
  double &delta_b      =  y_tc[Constants.ind_deltab_tc];        // δ_b(x,k)
  double &u_c          =  y_tc[Constants.ind_uc_tc];            // u_c(x,k)
  double &u_b          =  y_tc[Constants.ind_ub_tc];            // u_b(x,k)
  double &Phi          =  y_tc[Constants.ind_Phi_tc];           // Φ(x,k) 
  double *Theta        = &y_tc[Constants.ind_start_Theta_tc];   // Θ_ℓ(x,k)

  
  //  fetch cosmological parameters and variables:
  double Hp     = cosmo->Hp_of_x(x);        // Hp(x)
  double dHpdx  = cosmo->dHpdx_of_x(x);     // d/dx[Hp(x)]

  //  fetch recombination variables:
  double R          = rec->R_of_x(x);         // R(x)
  double dtaudx     = rec->dtaudx_of_x(x);    // d/dx[τ(x)]
  double ddtaudxx   = rec->ddtaudxx_of_x(x);  // d^2/dx^2[τ(x)]

  //  fetch perturbation variables:
  double Psi = -2./3.;    // Ψ(x=x_init,k)

  //  define useful factors:
  double ckH      = Constants.c*k/Hp;   // ck/Hp(x)


  //  Set initial conditions in the tight coupling regime

  //  compute scalar quantities (Φ, δ, u):

  Phi = - Psi; 

  delta_c = -3./2. * Psi;
  delta_b = delta_c;

  u_c = -0.5*ckH * Psi;
  u_b = u_c;

  // compute photon temperature perturbations (Θ_ℓ):

  Theta[0] = -0.5 * Psi;
  Theta[1] = ckH/6. * Psi;

  return y_tc;
}



Vector Perturbations::set_ic_after_tight_coupling(
    const Vector &y_tc, 
    const double x, 
    const double k) const{

  Vector y(Constants.n_ell_tot_full);   // vector to fill
  

  // Number of multipoles we have in the full regime


  const int n_ell_Theta         = Constants.n_ell_Theta;
  const int n_ell_Thetap        = Constants.n_ell_Thetap;
  const int n_ell_neutrinos     = Constants.n_ell_neutrinos;
  const bool polarization       = Constants.polarization;
  const bool neutrinos          = Constants.neutrinos;

  // Number of multipoles we have in the tight coupling regime


  const int n_ell_Theta_tc      = Constants.n_ell_Theta_tc;
  const int n_ell_neutrinos_tc  = Constants.n_ell_neutrinos_tc;

  //  reference to the tight coupling quantities:
  const double &delta_c_tc      =  y_tc[Constants.ind_deltac_tc];       // δ_c(x,k) from tight coupling
  const double &delta_b_tc      =  y_tc[Constants.ind_deltab_tc];       // δ_b(x,k) from tight coupling
  const double &u_c_tc          =  y_tc[Constants.ind_uc_tc];           // u_c(x,k) from tight coupling
  const double &u_b_tc          =  y_tc[Constants.ind_ub_tc];           // u_b(x,k) from tight coupling
  const double &Phi_tc          =  y_tc[Constants.ind_Phi_tc];          // Φ(x,k) from tight coupling
  const double *Theta_tc        = &y_tc[Constants.ind_start_Theta_tc];  // Θ_ℓ(x,k) from tight coupling

  //  reference to the quantities we are going to set:
  double &delta_c         =  y[Constants.ind_deltac];        // δ_c(x,k) after tight coupling
  double &delta_b         =  y[Constants.ind_deltab];        // δ_b(x,k) after tight coupling
  double &u_c             =  y[Constants.ind_uc];            // u_c(x,k) after tight coupling
  double &u_b             =  y[Constants.ind_ub];            // u_b(x,k) after tight coupling
  double &Phi             =  y[Constants.ind_Phi];           // Φ(x,k) after tight coupling
  double *Theta           = &y[Constants.ind_start_Theta];   // Θ_ℓ(x,k) after tight coupling


  //  fetch constants and variables:
  const double c              = Constants.c;              // c
  const double H0             = cosmo->get_H0();          // H_0


  //  compute scalar quantities (Φ, δ, u):
  Phi = Phi_tc;

  delta_c = delta_c_tc;
  delta_b = delta_b_tc;

  u_c = u_c_tc;
  u_b = u_b_tc;


  // compute photon temperature perturbations (Θ_ℓ):

  Theta[0] = Theta_tc[0];
  Theta[1] = Theta_tc[1];

  double U = c*k/(cosmo->Hp_of_x(x)*rec->dtaudx_of_x(x)) ;
  // Theta[2] = -8./15.*U * Theta[1]; //idk
  Theta[2] = expr_Theta2(x, k, Theta[1]);

  for(int ell=3; ell<n_ell_Theta; ell++){
    Theta[ell] = - (double)ell/(2.*ell + 1.) * U * Theta[ell-1];
  }

  return y;
}


double Perturbations::get_tight_coupling_time(const double k) const{

  // FIX THIS!


  double x_tc_end = 1000.0;
  int  idx_tc_end = n_x*10;

  Vector x_array = Utils::linspace(x_start, x_end, n_x);


  for (int i=0; i<n_x; i++){
    double x          = x_array[i];
    double ckH        = Constants.c * k / cosmo->Hp_of_x(x);
    double dtaudx     = -rec->dtaudx_of_x(x);

    bool cond1 = dtaudx < 10.0;
    bool cond2 = dtaudx < 10.0 * ckH;
    bool cond3 = x > -8.3;  

    if (cond1 || cond2 || cond3) {
      idx_tc_end = i - 1;
      x_tc_end = x_array[idx_tc_end];
      break; 
    }
  }


  // std::cout << "k = " << k*Constants.Mpc << " -> x_tc_end = " << x_tc_end << std::endl;
  return x_tc_end;
}



void Perturbations::compute_source_functions(int npts_x, int npts_k){

  const int n_k = npts_k;
  const int n_x = npts_x;
  
  Vector k_array = exp(Utils::linspace(log(k_min), log(k_max), n_k));
  Vector x_array = Utils::linspace(x_start, x_end, n_x);

  // Make storage for the source functions (in 1D array to be able to pass it to the spline)

  const int nn = n_k*n_x;

  Vector ST_array(nn);   // temperature source

  Vector SW_array(nn);
  Vector ISW_array(nn);
  Vector Doppler_array(nn);
  Vector pol_array(nn);

  //  fetch constants and variables:
  const double c              = Constants.c;              // c
  const double H0             = cosmo->get_H0();          // H_0
  const double H0H0           = H0*H0;                    // H_0^2
  const double Omegagamma0    = cosmo->get_Omegagamma();  // Ω_γ0
  const double Omegac0        = cosmo->get_Omegac();      // Ω_c0
  const double Omegab0        = cosmo->get_Omegab();      // Ω_b0

  //  compute source functions:
  for(auto ix = 0; ix < x_array.size(); ix++){
    const double x = x_array[ix];
    double a = exp(x);

    //  fetch all the things we need...:
    const double Hp       = cosmo->Hp_of_x(x);
    const double dHpdx    = cosmo->dHpdx_of_x(x);
    const double ddHpdxx  = cosmo->ddHpdxx_of_x(x);
    const double tau      = rec->tau_of_x(x);
    const double dtaudx   = rec->dtaudx_of_x(x);
    const double gt       = rec->gt_of_x(x);
    const double dgtdx    = rec->dgtdx_of_x(x);
    const double ddgtdxx  = rec->ddgtdxx_of_x(x);

    

    for(auto ik = 0; ik < k_array.size(); ik++){
      const double k = k_array[ik];

      // NB: This is the format the data needs to be stored 
      // in a 1D array for the 2D spline routine source(ix,ik) -> S_array[ix + nx * ik]
      const int index = ix + n_x * ik;

      //=============================================================================
      // TODO: Compute the source functions
      //=============================================================================
      
      //  fetch the rest of the things we need...:
      double ckH            = c*k/Hp;
      double ck_inv         = 1./(c*k);
      const double Psi      = get_Psi(x, k);        // Ψ(x,k)
      const double Phi      = get_Phi(x, k);        // Φ(x,k)
      const double Theta0   = get_Theta(x, k, 0);   // Θ_0(x,k)
      const double Theta1   = get_Theta(x, k, 1);   // Θ_1(x,k)
      const double Theta2   = get_Theta(x, k, 2);   // Θ_2(x,k)
      const double Theta3   = get_Theta(x, k, 3);   // Θ_3(x,k)
      const double u_b      = get_u_b(x, k);        // u_b(x,k)

      const double dPhidx       = Phi_spline.deriv_x(x, k);         // d/dx[Φ(x,k)]
      const double dPsidx       = Psi_spline.deriv_x(x, k);         // d/dx[Ψ(x,k)]
      const double du_bdx       = u_b_spline.deriv_x(x, k);         // d/dx[u_b(x,k)]
      const double dTheta2dx    = Theta_spline[2].deriv_x(x, k);    // d/dx[Θ_2(x,k)]
      const double ddTheta2dxx  = Theta_spline[2].deriv_xx(x, k);   // fixme!

      double ISW_term, SW_term, Doppler_term, pol_term;

      SW_term = gt * (Theta0 + Psi + 0.25*Theta2);

      ISW_term = exp(-tau)* ( dPsidx - dPhidx );

      Doppler_term = - ck_inv * ( dHpdx*gt*u_b + Hp*( dgtdx*u_b + gt*du_bdx ) );

      pol_term = .75*ck_inv*ck_inv * 
              ( 3*Hp*dHpdx * ( dgtdx*Theta2 + gt*dTheta2dx )
              + Hp*Hp * ( ddgtdxx*Theta2 + 2*dgtdx*dTheta2dx + gt*ddTheta2dxx ) 
              + ( dHpdx*dHpdx + Hp*ddHpdxx ) * gt*Theta2 );

      // Temperature source
      ST_array[index] = SW_term + ISW_term + Doppler_term + pol_term;

      SW_array[index]       = SW_term;
      ISW_array[index]      = ISW_term;
      Doppler_array[index]  = Doppler_term;
      pol_array[index]      = pol_term;

    }
  }

  // Spline the source functions
  ST_spline.create(x_array, k_array, ST_array, "Source_Temp_x_k");

  //  spline the different parts:
  ST_SW_spline.create(x_array, k_array, SW_array, "Source_Temp_x_k_SW");
  ST_ISW_spline.create(x_array, k_array, ISW_array, "Source_Temp_x_k_ISW");
  ST_Doppler_spline.create(x_array, k_array, Doppler_array, "Source_Temp_x_k_Doppler");
  ST_pol_spline.create(x_array, k_array, pol_array, "Source_Temp_x_k_polarisation");


}





//  ----------------------------------
//  R.H.S. OF ODES
//  ----------------------------------


int Perturbations::rhs_tight_coupling_ode(double x, double k, const double *y, double *dydx){

  
  // For integration of perturbations in tight coupling regime (Only 2 photon multipoles + neutrinos needed)

  const int n_ell_Theta_tc      = Constants.n_ell_Theta_tc;
  const int n_ell_neutrinos_tc  = Constants.n_ell_neutrinos_tc;
  const bool neutrinos          = Constants.neutrinos;

  //  reference to the different quantities in the y-array:
  const double &delta_c         =  y[Constants.ind_deltac_tc];      // δ_c(x,k)
  const double &delta_b         =  y[Constants.ind_deltab_tc];      // δ_b(x,k)
  const double &u_c             =  y[Constants.ind_uc_tc];          // u_c(x,k)
  const double &u_b             =  y[Constants.ind_ub_tc];          // u_b(x,k)
  const double &Phi             =  y[Constants.ind_Phi_tc];         // Φ(x,k)
  const double *Theta           = &y[Constants.ind_start_Theta_tc]; // Θ_ℓ(x,k)

  //  reference to the quantities we are going to set in the dydx array:
  double &ddelta_cdx      =  dydx[Constants.ind_deltac_tc];       // d/dx[δ_c(x,k)]
  double &ddelta_bdx      =  dydx[Constants.ind_deltab_tc];       // d/dx[δ_b(x,k)]
  double &du_cdx          =  dydx[Constants.ind_uc_tc];           // d/dx[u_c(x,k)]
  double &du_bdx          =  dydx[Constants.ind_ub_tc];           // d/dx[u_b(x,k)]
  double &dPhidx          =  dydx[Constants.ind_Phi_tc];          // d/dx[Φ(x,k)]
  double *dThetadx        = &dydx[Constants.ind_start_Theta_tc];  // d/dx[Θ_ℓ(x,k)]  

  //  fetch cosmological parameters and variables:
  double Hp             = cosmo->Hp_of_x(x);          // Hp(x)
  double dHpdx          = cosmo->dHpdx_of_x(x);       // d/dx[Hp(x)]
  double H0             = cosmo->get_H0();            // H_0
  double Omegagamma0    = cosmo->get_Omegagamma();    // Ω_γ0
  double Omegac0        = cosmo->get_Omegac();        // Ω_c0
  double Omegab0        = cosmo->get_Omegab();        // Ω_b0 
  double eta            = cosmo->eta_of_x(x);         // η(x)

  //  fetch recombination variables::
  double R          = rec->R_of_x(x);           // R(x)
  double dtaudx     = rec->dtaudx_of_x(x);      // d/dx[τ(x)]
  double ddtaudxx   = rec->ddtaudxx_of_x(x);    // d^2/dx^2[τ(x)]

  

  // useful factors:
  double c      = Constants.c;  // c
  double ckH    = c*k/Hp;       // ck/Hp
  double Rplus  = (1.+R);       // ( 1 + R )
  double Rminus = (1.-R);       // ( 1 - R )
  double H0H0   = H0*H0;        // H_0^2
  double a      = exp(x);       // e^x

  double expr;

  //  fetch perturbation variables:
  const double Theta2 = expr_Theta2(x, k, Theta[1]);    // Θ_2(x,k)
  const double Psi    = expr_Psi(x, k, Phi, Theta2);    // Ψ(x,k)

  
  // reocurring terms:
  double U1 = ( 1 - dHpdx/Hp );               // ( 1 - 1/Hp*d/dx[Hp] )
  double U2 = ( -Theta[0] + 2.*Theta2 );      // ( - Θ_0 + 2Θ_2 )


  //  compute scalar quantities (Φ, δ, u):

  expr = (Omegac0*delta_c + Omegab0*delta_b)/a + 4.*Omegagamma0*Theta[0]/(a*a);
  dPhidx = Psi - 1./3*ckH*ckH*Phi + 0.5*(H0H0/(Hp*Hp)) * expr;
  
  expr = 3.*dPhidx;
  ddelta_cdx = ckH*u_c - expr;
  ddelta_bdx = ckH*u_b - expr;
  
  expr = ckH*Psi;
  du_cdx = - u_c - expr;

  //  compute photon temperature perturbations (Θ_0):

  dThetadx[0] = -ckH * Theta[1] - dPhidx;

  //  compute scalar quantities (u_b):

  double denom = ( dtaudx*Rplus - R*U1 );
  double term1 = - ( ddtaudxx*Rplus - dtaudx*Rminus ) * (3*Theta[1] + u_b);
  double term2 = - ckH*R * ( Psi - U1*U2 + dThetadx[0]);
  double q = (term1+term2)/denom;

  du_bdx = ( q - R*u_b + ckH*U2)/Rplus - expr;

  //  compute photon temperature perturbations (Θ_1):
  
  dThetadx[1] = 1/3.*( q - du_bdx );

  return GSL_SUCCESS;
}


int Perturbations::rhs_full_ode(double x, double k, const double *y, double *dydx){

  //  index and number of the different quantities:
  const int n_ell_Theta         = Constants.n_ell_Theta;
  const int n_ell_Thetap        = Constants.n_ell_Thetap;
  const int n_ell_neutrinos     = Constants.n_ell_neutrinos;
  const bool polarization       = Constants.polarization;
  const bool neutrinos          = Constants.neutrinos;

  //  reference to the different quantities in the y array:
  const double &delta_c         =  y[Constants.ind_deltac];       // δ_c(x,k)
  const double &delta_b         =  y[Constants.ind_deltab];       // δ_b(x,k)
  const double &u_c             =  y[Constants.ind_uc];           // u_c(x,k)
  const double &u_b             =  y[Constants.ind_ub];           // u_c(x,k)
  const double &Phi             =  y[Constants.ind_Phi];          // Φ(x,k)
  const double *Theta           = &y[Constants.ind_start_Theta];  // Θ_ℓ(x,k) 

  //  reference to the quantities we are going to set in the dydx array:
  double &ddelta_cdx      =  dydx[Constants.ind_deltac];        // d/dx[δ_c(x,k)]
  double &ddelta_bdx      =  dydx[Constants.ind_deltab];        // d/dx[δ_b(x,k)]
  double &du_cdx          =  dydx[Constants.ind_uc];            // d/dx[u_c(x,k)]
  double &du_bdx          =  dydx[Constants.ind_ub];            // d/dx[u_b(x,k)]
  double &dPhidx          =  dydx[Constants.ind_Phi];           // d/dx[Φ(x,k)]
  double *dThetadx        = &dydx[Constants.ind_start_Theta];   // d/dx[Θ_ℓ(x,k)]

  //  fetch cosmological parameters and variables:
  const double Hp           = cosmo->Hp_of_x(x);          // Hp(x)
  const double dHpdx        = cosmo->dHpdx_of_x(x);       // d/dx[Hp(x)]
  const double eta          = cosmo->eta_of_x(x);         // η(x)
  const double H0           = cosmo->get_H0();            // H_0
  const double Omegagamma0  = cosmo->get_Omegagamma();    // Ω_γ0
  const double Omegac0      = cosmo->get_Omegac();        // Ω_c0
  const double Omegab0      = cosmo->get_Omegab();        // Ω_b0

  //  fetch recombination variables:
  const double R          = rec->R_of_x(x);           // R(x)
  const double dtaudx     = rec->dtaudx_of_x(x);      // d/dx[τ(x)]

  //  perturbation variables:
  const double Psi = expr_Psi(x, k, Phi, Theta[2]);  // Ψ(x,k)

  
  // useful factors:
  const double c      = Constants.c;    // c
  const double ckH    = c*k/Hp;         // ck/Hp
  const double H0H0   = H0*H0;          // (H_0)^2
  const double a      = exp(x);         // e^x

  double expr;

  //  compute derivatives of scalar quantities (Φ, δ, u):

  expr = ( Omegac0*delta_c + Omegab0*delta_b) / a + 4.*Omegagamma0*Theta[0]/(a*a);
  dPhidx = Psi - ckH*ckH*Phi/3. + 0.5*(H0H0/(Hp*Hp)) * expr;
  
  expr = 3.*dPhidx;
  ddelta_cdx = ckH*u_c - expr;
  ddelta_bdx = ckH*u_b - expr;

  expr = ckH*Psi;
  du_cdx = - u_c - expr;
  du_bdx = - u_b - expr + dtaudx * ( 3.*Theta[1] + u_b ) / R;

  //  compute derivatives of photon multipoles (Θ_ℓ):

  int ell           = 0;                // ℓ
  const int ell_max = n_ell_Theta - 1;  // ℓ_max

  dThetadx[0] = -ckH * Theta[1] - dPhidx;
  dThetadx[1] = ckH/3. * Theta[0] - 2.*ckH/3.*Theta[2] + ckH/3.*Psi + dtaudx*( Theta[1] + u_b/3. );
  
  for(int ell=2; ell<ell_max; ell++){
    expr = ckH/(2.*ell + 1.) * ( ell*Theta[ell-1] - (ell+1.)*Theta[ell+1]);
    dThetadx[ell] = expr + dtaudx * Theta[ell];
  }
  // dThetadx[1] += ckH*Psi + 1./3*dtaudx*u_b;
  ell = 2;
  dThetadx[ell] -= dtaudx*Theta[ell]/10.;
  ell = ell_max;
  dThetadx[ell] = ckH*Theta[ell-1] - c*(ell+1.) / (Hp*eta) * Theta[ell] + dtaudx*Theta[ell];

  return GSL_SUCCESS;
}




//  ----------------------------------
//  REOCURRING EXPRESSIONS
//  ----------------------------------

double Perturbations::expr_Psi(double x, double k, double Phi, double Theta2) const{
  double expr = cosmo->get_H0()/ ( Constants.c*k);
  return - Phi - 12.*expr*expr*exp(-2*x) * cosmo->get_Omegagamma()*Theta2;
}

double Perturbations::expr_Theta2(double x, double k, double Theta1) const{
  double fac =  (20. * Constants.c*k ) / ( 45. * cosmo->Hp_of_x(x) * rec->dtaudx_of_x(x) );
  return -fac*Theta1;
}

double Perturbations::expr_Thetaell(double x, double k, int ell, double Thetaell_prev) const{
  double fac =  ell/(2.*ell + 1.)*( Constants.c*k ) / (cosmo->Hp_of_x(x) * rec->dtaudx_of_x(x) );
  return -fac*Thetaell_prev;
}



//  ----------------------------------
//  GET PERTURBATION QUANTITIES
//  ----------------------------------


double Perturbations::get_delta_c(const double x, const double k) const{
  return delta_c_spline(x,k);
}
double Perturbations::get_delta_b(const double x, const double k) const{
  return delta_b_spline(x,k);
}
double Perturbations::get_u_c(const double x, const double k) const{
  return u_c_spline(x,k);
}
double Perturbations::get_u_b(const double x, const double k) const{
  return u_b_spline(x,k);
}
double Perturbations::get_Phi(const double x, const double k) const{
  return Phi_spline(x,k);
}
double Perturbations::get_Psi(const double x, const double k) const{
  return Psi_spline(x,k);
}
double Perturbations::get_Pi(const double x, const double k) const{
  return Theta_spline[2](x,k);
}
double Perturbations::get_Source_T(const double x, const double k) const{
  return ST_spline(x,k);
}
double Perturbations::get_Source_T(const double x, const double k, const int term) const{
  if(term==1)
    return ST_SW_spline(x,k);
  else if(term==2)
    return ST_ISW_spline(x,k);
  else if(term==3)
    return ST_Doppler_spline(x,k);
  else if(term==4)
    return ST_pol_spline(x,k);
  else
    return ST_spline(x,k);
}


double Perturbations::get_Theta(const double x, const double k, const int ell) const{
  return Theta_spline[ell](x,k);
}




//  ----------------------------------
//  HANDLE I/O
//  ----------------------------------



void Perturbations::info() const{
  std::cout << "\n";
  std::cout << "Info about perturbations class:\n";
  std::cout << "x_start:       " << x_start                << "\n";
  std::cout << "x_end:         " << x_end                  << "\n";
  std::cout << "n_x:     " << n_x              << "\n";
  std::cout << "k_min (1/Mpc): " << k_min * Constants.Mpc  << "\n";
  std::cout << "k_max (1/Mpc): " << k_max * Constants.Mpc  << "\n";
  std::cout << "n_k:     " << n_k              << "\n";
  if(Constants.polarization)
    std::cout << "We include polarization\n";
  else
    std::cout << "We do not include polarization\n";
  if(Constants.neutrinos)
    std::cout << "We include neutrinos\n";
  else
    std::cout << "We do not include neutrinos\n";

  std::cout << "Information about the perturbation system:\n";
  std::cout << "ind_deltac:         " << Constants.ind_deltac           << "\n";
  std::cout << "ind_deltab:         " << Constants.ind_deltab           << "\n";
  std::cout << "ind_u_c:            " << Constants.ind_uc               << "\n";
  std::cout << "ind_u_b:            " << Constants.ind_ub               << "\n";
  std::cout << "ind_Phi:            " << Constants.ind_Phi              << "\n";
  std::cout << "ind_start_Theta:    " << Constants.ind_start_Theta      << "\n";
  std::cout << "n_ell_Theta:        " << Constants.n_ell_Theta          << "\n";
  if(Constants.polarization){
    std::cout << "ind_start_Thetap:   " << Constants.ind_start_Thetap   << "\n";
    std::cout << "n_ell_Thetap:       " << Constants.n_ell_Thetap       << "\n";
  }
  if(Constants.neutrinos){
    std::cout << "ind_start_nu:       " << Constants.ind_start_nu       << "\n";
    std::cout << "n_ell_neutrinos     " << Constants.n_ell_neutrinos    << "\n";
  }
  std::cout << "n_ell_tot_full:     " << Constants.n_ell_tot_full       << "\n";

  std::cout << "Information about the perturbation system in tight coupling:\n";
  std::cout << "ind_deltac:         " << Constants.ind_deltac_tc        << "\n";
  std::cout << "ind_deltab:         " << Constants.ind_deltab_tc        << "\n";
  std::cout << "ind_u_c:            " << Constants.ind_uc_tc            << "\n";
  std::cout << "ind_u_b:            " << Constants.ind_ub_tc            << "\n";
  std::cout << "ind_Phi:            " << Constants.ind_Phi_tc           << "\n";
  std::cout << "ind_start_Theta:    " << Constants.ind_start_Theta_tc   << "\n";
  std::cout << "n_ell_Theta:        " << Constants.n_ell_Theta_tc       << "\n";
  if(Constants.neutrinos){
    std::cout << "ind_start_nu:       " << Constants.ind_start_nu_tc    << "\n";
    std::cout << "n_ell_neutrinos     " << Constants.n_ell_neutrinos_tc << "\n";
  }
  std::cout << "n_ell_tot_tc:       " << Constants.n_ell_tot_tc         << "\n";
  std::cout << std::endl;
}



void Perturbations::output(const double k, const std::string filename) const{
  std::ofstream fp(OUTPUT_PATH + filename.c_str());
  const int npts = 5000+1;
  auto x_array = Utils::linspace(x_start>=-18 ? x_start:-18, x_end, npts);
  auto print_data = [&] (const double x) {
    double arg = k * (cosmo->eta_of_x(0.0) - cosmo->eta_of_x(x));
    fp << x                   << " ";     // 0
    fp << get_delta_c(x,k)    << " ";     // 1
    fp << get_delta_b(x,k)    << " ";     // 2
    fp << get_u_c(x,k)        << " ";     // 3
    fp << get_u_b(x,k)        << " ";     // 4
    fp << get_Theta(x,k,0)    << " ";     // 5 
    fp << get_Theta(x,k,1)    << " ";     // 6
    fp << get_Theta(x,k,2)    << " ";     // 7
    fp << get_Phi(x,k)        << " ";     // 8
    fp << get_Psi(x,k)        << " ";     // 9
    fp << get_Pi(x,k)         << " ";     // 10
    fp << get_Source_T(x,k)  << " ";
    fp << get_Source_T(x,k) * Utils::j_ell(5,   arg)           << " ";
    fp << get_Source_T(x,k) * Utils::j_ell(50,  arg)           << " ";
    fp << get_Source_T(x,k) * Utils::j_ell(500, arg)           << " ";
    fp << "\n";
  };
  std::for_each(x_array.begin(), x_array.end(), print_data);
}




/*

\Psi : Ψ
\Phi : Φ
\Theta : Θ
\Omega : Ω



\gamma : γ
\delta : δ
\eta : η
\tau : τ



\ell : ℓ

*/