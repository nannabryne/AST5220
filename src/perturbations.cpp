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

void Perturbations::solve(){

  //  integrate all the perturbation equation and spline the result:
  Utils::StartTiming("integrateperturbation");
  integrate_perturbations();
  Utils::EndTiming("integrateperturbation");

  //  compute source functions and spline the result:
  Utils::StartTiming("source");
  compute_source_functions();
  Utils::EndTiming("source");
}

//====================================================
// The main work: integrate all the perturbations
// and spline the results
//====================================================

void Perturbations::integrate_perturbations(){
  

  //===================================================================
  // TODO: Set up the k-array for the k's we are going to integrate over
  // Start at k_min end at k_max with n_k points with either a
  // quadratic or a logarithmic spacing
  //===================================================================

  const gsl_odeiv2_step_type * stepper = gsl_odeiv2_step_rkf45;

  //  set up k-array for the k's we are to integrate over:
  Vector k_array = exp(Utils::linspace(log(k_min), log(k_max), n_k));
  

  Vector x_array = Utils::linspace(x_start, x_end, n_x);

  Vector delta_c_array(n_x*n_k);
  Vector delta_b_array(n_x*n_k);
  Vector u_c_array(n_x*n_k);
  Vector u_b_array(n_x*n_k);
  Vector Phi_array(n_x*n_k);
  Vector Theta0_array(n_x*n_k);
  Vector Theta1_array(n_x*n_k);
  Vector Theta2_array(n_x*n_k);
  Vector Psi_array(n_x*n_k);

  //  loop over all wavenumbers:  
  for(int ik=0; ik<n_k; ik++){

    //  progress bar...
    if( (10*ik) / n_k != (10*ik+10) / n_k ) {
      std::cout << (100*ik+100)/n_k << "% " << std::flush;
      if(ik == n_k-1) std::cout << std::endl;
    }

    double k = k_array[ik];   // current value of k

    //  find value to integrate to:
    double x_end_tight = get_tight_coupling_time(k);
    double idk = n_x*(x_end_tight-x_start)/(x_end-x_start);
    int n_x_tight = int(idk)+1;
    int n_x_full = n_x-n_x_tight;
    // Vector x_array_tc = Utils::linspace(x_start, x_end_tight, n_x_tight);
    // Vector x_array_full = Utils::linspace(x_end_tight, x_end, n_x_full);

    Vector x_array_tc = Vector(x_array.begin(), x_array.begin()+n_x_tight);
    Vector x_array_full = Vector(x_array.begin()+n_x_tight-1, x_array.end());



    //===================================================================
    // TODO: Tight coupling integration
    // Remember to implement the routines:
    // set_ic : The IC at the start
    // rhs_tight_coupling_ode : The dydx for our coupled ODE system
    //===================================================================
    

    //  set up initial conditions in the tight coupling regime:
    auto y_tight_coupling_ini = set_ic(x_start, k);

    // The tight coupling ODE system

    ODEFunction dydx_tight_coupling = [&](double x, const double *y, double *dydx){
      return rhs_tight_coupling_ode(x, k, y, dydx);
    };

   
    // Integrate from x_start -> x_end_tight
    // ...
    // ...
    // ...
    // ...
    // ...
    // std::cout << y_tight_coupling_ini[Constants.ind_start_Theta_tc+1] << std::endl;

    ODESolver ode_tight;
    ode_tight.set_accuracy(1e-5, 1e-8, 1e-8);
    ode_tight.solve(dydx_tight_coupling, x_array_tc, y_tight_coupling_ini, stepper);
    Vector y_tight_coupling = ode_tight.get_final_data();
    Vector2D y_sol_tc = ode_tight.get_data();

    // std::cout << y_tight_coupling[0] << ", " << y_tight_coupling[2] << std::endl;


    //====i===============================================================
    // TODO: Full equation integration
    // Remember to implement the routines:
    // set_ic_after_tight_coupling : The IC after tight coupling ends
    // rhs_full_ode : The dydx for our coupled ODE system
    //===================================================================

    // Set up initial conditions (y_tight_coupling is the solution at the end of tight coupling)
    auto y_full_ini = set_ic_after_tight_coupling(y_tight_coupling, x_end_tight, k);

    // The full ODE system
    ODEFunction dydx_full = [&](double x, const double *y, double *dydx){
      return rhs_full_ode(x, k, y, dydx);
    };

    // Integrate from x_end_tight -> x_end
    ODESolver ode_full;
    ode_full.set_accuracy(1e-5, 1e-8, 1e-8);
    ode_full.solve(dydx_full, x_array_full, y_full_ini, stepper);
    Vector2D y_sol_full = ode_full.get_data();
    

    //  Store the data 


    int d = n_x*ik;

    double x;
    int i, j;
    Vector y_sol_curr;


    //  during tight coupling:
    for(int ix=0; ix<n_x_tight; ix++){
      i = ix + d;
      x = x_array_tc[ix];
      y_sol_curr = y_sol_tc[ix];

      delta_c_array[i]    = y_sol_curr[Constants.ind_deltac_tc];
      delta_b_array[i]    = y_sol_curr[Constants.ind_deltab_tc];
      u_c_array[i]        = y_sol_curr[Constants.ind_uc_tc];
      u_b_array[i]        = y_sol_curr[Constants.ind_ub_tc];
      Theta0_array[i]     = y_sol_curr[Constants.ind_start_Theta_tc];
      Theta1_array[i]     = y_sol_curr[Constants.ind_start_Theta_tc+1];
      Theta2_array[i]     = expr_Theta2(x, k, Theta1_array[i]);
      Phi_array[i]        = y_sol_curr[Constants.ind_Phi_tc];
      Psi_array[i]        = expr_Psi(x, k, Phi_array[i], Theta2_array[i]);
    }
    //  after tight coupling:
    for(int ix=n_x_tight; ix<n_x; ix++){
      i = ix + d;
      j = ix - n_x_tight;
      x = x_array_full[j];
      y_sol_curr = y_sol_full[j];

      delta_c_array[i]    = y_sol_curr[Constants.ind_deltac];
      delta_b_array[i]    = y_sol_curr[Constants.ind_deltab];
      u_c_array[i]        = y_sol_curr[Constants.ind_uc];
      u_b_array[i]        = y_sol_curr[Constants.ind_ub];
      Theta0_array[i]     = y_sol_curr[Constants.ind_start_Theta];
      Theta1_array[i]     = y_sol_curr[Constants.ind_start_Theta+1];
      Theta2_array[i]     = y_sol_curr[Constants.ind_start_Theta+2];
      Phi_array[i]        = y_sol_curr[Constants.ind_Phi];
      Psi_array[i]        = expr_Psi(x, k, Phi_array[i], Theta2_array[i]);
    }
    


  }
  

  //=============================================================================
  // TODO: Make all splines needed: Theta0,Theta1,Theta2,Phi,Psi,...
  //=============================================================================
  // ...
  // ...
  // ...

  delta_c_spline.create(x_array, k_array, delta_c_array, "delta_c");
  delta_b_spline.create(x_array, k_array, delta_b_array, "delta_b");
  u_c_spline.create(x_array, k_array, u_c_array, "u_c");
  u_b_spline.create(x_array, k_array, u_b_array, "u_b");

  Phi_spline.create(x_array, k_array, Phi_array, "Phi");
  Psi_spline.create(x_array, k_array, Psi_array, "Psi");

  Theta_spline = std::vector<Spline2D>(3);
  Theta_spline[0].create(x_array, k_array, Theta0_array, "Theta0");
  Theta_spline[1].create(x_array, k_array, Theta1_array, "Theta1");
  Theta_spline[2].create(x_array, k_array, Theta2_array, "Theta2");
  

}



//====================================================
// Set IC at the start of the run (this is in the
// tight coupling regime)
//====================================================
Vector Perturbations::set_ic(const double x, const double k) const{

  Vector y_tc(Constants.n_ell_tot_tc);  // the vector we are going to fill

  //=============================================================================
  // Compute where in the y_tc array each component belongs
  // This is just an example of how to do it to make it easier
  // Feel free to organize the component any way you like
  //=============================================================================
  
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
  // double *Nu           = &y_tc[Constants.ind_start_nu_tc];

  
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


//====================================================
// Set IC for the full ODE system after tight coupling 
// regime ends
//====================================================

Vector Perturbations::set_ic_after_tight_coupling(
    const Vector &y_tc, 
    const double x, 
    const double k) const{

  Vector y(Constants.n_ell_tot_full);   // vector to fill
  
  //=============================================================================
  // Compute where in the y array each component belongs and where corresponding
  // components are located in the y_tc array
  // This is just an example of how to do it to make it easier
  // Feel free to organize the component any way you like
  //=============================================================================

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
  // const double *Nu_tc           = &y_tc[Constants.ind_start_nu_tc];

  //  reference to the quantities we are going to set:
  double &delta_c         =  y[Constants.ind_deltac];        // δ_c(x,k) after tight coupling
  double &delta_b         =  y[Constants.ind_deltab];        // δ_b(x,k) after tight coupling
  double &u_c             =  y[Constants.ind_uc];            // u_c(x,k) after tight coupling
  double &u_b             =  y[Constants.ind_ub];            // u_b(x,k) after tight coupling
  double &Phi             =  y[Constants.ind_Phi];           // Φ(x,k) after tight coupling
  double *Theta           = &y[Constants.ind_start_Theta];   // Θ_ℓ(x,k) after tight coupling
  // double *Theta_p         = &y[Constants.ind_start_Thetap];
  // double *Nu              = &y[Constants.ind_start_nu];


  //  fetch constants and variables:
  const double c              = Constants.c;              // c
  const double H0             = cosmo->get_H0();          // H_0
  // const double H0H0           = H0*H0;                    // H_0^2
  // const double Omegagamma0    = cosmo->get_Omegagamma();  // Ω_γ0


  //  compute scalar quantities (Φ, δ, u):
  Phi = Phi_tc;

  delta_c = delta_c_tc;
  delta_b = delta_b_tc;

  u_c = u_c_tc;
  u_b = u_b_tc;


  // compute photon temperature perturbations (Θ_ℓ):

  //FIXME
  int ell       = 0;                // ℓ
  int ell_max   = n_ell_Theta - 1;  // ℓ_max

  for(ell=0; ell<2; ell++){
    Theta[ell] = Theta_tc[ell];
  }

  double U = c*k/(cosmo->Hp_of_x(x)*rec->dtaudx_of_x(x)) ;
  // Theta[2] = -8./15.*U * Theta[1]; //idk
  Theta[2] = expr_Theta2(x, k, Theta[1]);

  for(ell=3; ell<=ell_max; ell++){
    Theta[ell] = - (double)ell/(2.*ell + 1.) * U * Theta[ell-1];
  }

  return y;
}

//====================================================
// The time when tight coupling end
//====================================================

double Perturbations::get_tight_coupling_time(const double k) const{
  

  //=============================================================================
  // TODO: compute and return x for when tight coupling ends
  // Remember all the three conditions in Callin
  //=============================================================================
  // ...
  // ...

  //  Find x for when the tight coupling is over

  double x_rec = -8.3;
  int npts = int((1e6+1)*(x_rec-x_start)/(x_end-x_start));
  Vector x_array = Utils::linspace(x_start, x_rec, npts);
  
  bool one, two, three;
  double c = Constants.c;
  double x, val, lim;

  double x_tight_coupling_end = x_rec;
  
  for(int i=0; i<npts; i++){
    x = x_array[i];
    
    lim = c*k/cosmo->Hp_of_x(x);
    if(lim>1)
      lim = 1;
    lim *= 10.;

    val = abs(rec->dtaudx_of_x(x));


    if(val<=10 || val<=10*c*k/cosmo->Hp_of_x(x)){
      x_tight_coupling_end = x;
      break;
    }

    if(val<=lim){
      x_tight_coupling_end = x;
      break;
    }

  }

  // std::cout << "k = " << k*Constants.Mpc << " -> x_tc_end = " << x_tight_coupling_end << std::endl;
  return x_tight_coupling_end;
}

//====================================================
// After integrsating the perturbation compute the
// source function(s)
//====================================================
void Perturbations::compute_source_functions(){
  

  //=============================================================================
  // TODO: Make the x and k arrays to evaluate over and use to make the splines
  //=============================================================================
  // ...
  // ...


  Vector k_array = exp(Utils::linspace(log(k_min), log(k_max), n_k));
  Vector x_array = Utils::linspace(x_start, x_end, n_x);

  // Make storage for the source functions (in 1D array to be able to pass it to the spline)

  Vector ST_array(k_array.size() * x_array.size());   // temperature source
  Vector SE_array(k_array.size() * x_array.size());

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
      const double Phi      = get_Phi(x, k);        // 
      const double Theta0   = get_Theta(x, k, 0);   //
      const double Theta2   = get_Theta(x, k, 2);   //
      const double u_b      = get_u_b(x, k);        //

      const double dPhidx       = Phi_spline.deriv_x(x, k);
      const double dPsidx       = Psi_spline.deriv_x(x, k);
      const double du_bdx       = u_b_spline.deriv_x(x, k);
      const double dTheta2dx    = Theta_spline[2].deriv_x(x, k);
      const double ddTheta2dxx  = Theta_spline[2].deriv_xx(x, k);

      double term1, term2, term3, term4;

      term1 = gt * (Theta0 + Psi + 0.25*Theta2);

      term2 = exp(-tau)* ( dPsidx - dPhidx );

      term3 = - ck_inv * ( dHpdx*gt*u_b + Hp*( dgtdx*u_b + gt*du_bdx ) );

      term4 = .75*ck_inv*ck_inv * 
              ( 3*Hp*dHpdx * ( dgtdx*Theta2 + gt*dTheta2dx )
              + Hp*Hp * ( ddgtdxx*Theta2 + 2*dgtdx*dTheta2dx + gt*ddTheta2dxx ) 
              + ( dHpdx*dHpdx + Hp*ddHpdxx ) * gt*Theta2 );


      // Temperature source
      ST_array[index] = term1 + term2 + term3 + term4;

      // Polarization source
      if(Constants.polarization){
        SE_array[index] = 0.0;
      }
    }
  }

  // Spline the source functions
  ST_spline.create (x_array, k_array, ST_array, "Source_Temp_x_k");
  if(Constants.polarization){
    SE_spline.create (x_array, k_array, SE_array, "Source_Pol_x_k");
  }
}





//====================================================
// The right hand side of the perturbations ODE
// in the tight coupling regime
//====================================================

// Derivatives in the tight coupling regime
int Perturbations::rhs_tight_coupling_ode(double x, double k, const double *y, double *dydx){

  //=============================================================================
  // Compute where in the y / dydx array each component belongs
  // This is just an example of how to do it to make it easier
  // Feel free to organize the component any way you like
  //=============================================================================
  
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
  const double Theta2 = expr_Theta2(x, k, Theta[1]);  // Θ_2(x,k)
  const double Psi = expr_Psi(x, k, Phi, Theta2);     // Ψ(x,k)

  
  // reocurring terms:
  double U1 = ( 1 - dHpdx/Hp );             // ( 1 - 1/Hp*d/dx[Hp] )
  double U2 = ( -Theta[0] + 2*Theta2 );     // ( - Θ_0 + 2Θ_2 )


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
  
  dThetadx[1] = ( q - du_bdx )/3.;

  return GSL_SUCCESS;
}



//====================================================
// The right hand side of the full ODE
//====================================================

int Perturbations::rhs_full_ode(double x, double k, const double *y, double *dydx){
  
  //=============================================================================
  // Compute where in the y / dydx array each component belongs
  // This is just an example of how to do it to make it easier
  // Feel free to organize the component any way you like
  //=============================================================================

  // Index and number of the different quantities
  const int n_ell_Theta         = Constants.n_ell_Theta;
  const int n_ell_Thetap        = Constants.n_ell_Thetap;
  const int n_ell_neutrinos     = Constants.n_ell_neutrinos;
  const bool polarization       = Constants.polarization;
  const bool neutrinos          = Constants.neutrinos;

  // The different quantities in the y array
  const double &delta_c         =  y[Constants.ind_deltac];       // δ_c(x,k)
  const double &delta_b         =  y[Constants.ind_deltab];       // δ_b(x,k)
  const double &u_c             =  y[Constants.ind_uc];           // u_c(x,k)
  const double &u_b             =  y[Constants.ind_ub];           // u_c(x,k)
  const double &Phi             =  y[Constants.ind_Phi];          // Φ(x,k)
  const double *Theta           = &y[Constants.ind_start_Theta];  // Θ_ℓ(x,k) 

  // References to the quantities we are going to set in the dydx array
  double &ddelta_cdx      =  dydx[Constants.ind_deltac];        // d/dx[δ_c(x,k)]
  double &ddelta_bdx      =  dydx[Constants.ind_deltab];        // d/dx[δ_b(x,k)]
  double &du_cdx          =  dydx[Constants.ind_uc];            // d/dx[u_c(x,k)]
  double &du_bdx          =  dydx[Constants.ind_ub];            // d/dx[u_b(x,k)]
  double &dPhidx          =  dydx[Constants.ind_Phi];           // d/dx[Φ(x,k)]
  double *dThetadx        = &dydx[Constants.ind_start_Theta];   // d/dx[Θ_ℓ(x,k)]

  // Cosmological parameters and variables 
  const double Hp           = cosmo->Hp_of_x(x);          // Hp(x)
  const double dHpdx        = cosmo->dHpdx_of_x(x);       // d/dx[Hp(x)]
  const double eta          = cosmo->eta_of_x(x);         // η(x)
  const double H0           = cosmo->get_H0();            // H_0
  const double Omegagamma0  = cosmo->get_Omegagamma();    // Ω_γ0
  const double Omegac0      = cosmo->get_Omegac();        // Ω_c0
  const double Omegab0      = cosmo->get_Omegab();        // Ω_b0

  // Recombination variables
  const double R          = rec->R_of_x(x);           // R(x)
  const double dtaudx     = rec->dtaudx_of_x(x);      // d/dx[τ(x)]
  // const double ddtaudxx   = rec->ddtaudxx_of_x(x);    // d^2/dx^2[τ(x)]

  
  //  useful factors
  const double c      = Constants.c;    // c
  const double ckH    = c*k/Hp;         // ck/Hp
  // const double Rplus  = (1.+R);         // ( 1 + R )
  // const double Rminus = (1.-R);         // ( 1 - R ) 
  const double H0H0   = H0*H0;          // (H_0)^2
  const double a      = exp(x);         // e^x

  double expr;

  // Perturbation variables
  const double Psi = expr_Psi(x, k, Phi, Theta[2]);  // Ψ(x,k)

  //  compute derivatives of scalar quantities (Φ, δ, u):

  expr = Omegac0*delta_c/a + Omegab0*delta_b/a + 4.*Omegagamma0*Theta[0]/(a*a);
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
  
  double dell;  // ℓ (double)
  double denom;
  for(ell=2; ell<ell_max; ell++){
    // dell = 1.*ell;
    denom = (2.*ell + 1.);
    expr = ell*ckH*Theta[ell-1]/denom - (ell+1.)*ckH*Theta[ell+1]/denom;

    // expr = ckH  *  dell*Theta[ell-1]/ (2.*dell + 1.) - ckH*(dell+1.)*Theta[ell+1] / (2.*dell + 1.);
    
    dThetadx[ell] = expr + dtaudx * Theta[ell];
  }
  // dThetadx[1] += ckH*Psi + 1./3*dtaudx*u_b;
  dThetadx[2] -= dtaudx*Theta[2]/10.;
  
  ell = ell_max;
  dThetadx[ell] = ckH*dThetadx[ell-1] + -c*(ell+1.)/(Hp*eta)*Theta[ell] + dtaudx*Theta[ell];


  return GSL_SUCCESS;
}






double Perturbations::expr_Psi(double x, double k, double Phi, double Theta2) const{
  double expr = cosmo->get_H0()/ ( Constants.c*k * exp(x));
  return - Phi - 12.*expr*expr * cosmo->get_Omegagamma()*Theta2;
}


double Perturbations::expr_Theta2(double x, double k, double Theta1) const{
  double fac =  (20. * Constants.c*k ) / ( 45. * cosmo->Hp_of_x(x) * rec->dtaudx_of_x(x) );
  return -fac*Theta1;
}





//====================================================
// Get methods
//====================================================



//  ----------------------------------
//  GET METHODS
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
double Perturbations::get_Source_E(const double x, const double k) const{
  return SE_spline(x,k);
}
double Perturbations::get_Theta(const double x, const double k, const int ell) const{
  return Theta_spline[ell](x,k);
}
double Perturbations::get_Theta_p(const double x, const double k, const int ell) const{
  return Theta_p_spline[ell](x,k);
}
double Perturbations::get_Nu(const double x, const double k, const int ell) const{
  return Nu_spline[ell](x,k);
}

//====================================================
// Print some useful info about the class
//====================================================

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

//====================================================
// Output some results to file for a given value of k
//====================================================

void Perturbations::output(const double k, const std::string filename) const{
  std::ofstream fp(OUTPUT_PATH + filename.c_str());
  const int npts = 5000;
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