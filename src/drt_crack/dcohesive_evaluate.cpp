#include "dcohesive.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_inpar/inpar_parameterlist_utils.H"

#include "../drt_lib/drt_globalproblem.H" //can be removed???

int DRT::ELEMENTS::Dcohesive::Evaluate(Teuchos::ParameterList&    params,
                                    DRT::Discretization&      discretization,
                                    std::vector<int>&         lm,
                                    Epetra_SerialDenseMatrix& elemat1,
                                    Epetra_SerialDenseMatrix& elemat2,
                                    Epetra_SerialDenseVector& elevec1,
                                    Epetra_SerialDenseVector& elevec2,
                                    Epetra_SerialDenseVector& elevec3)
{

  // the spring is already failed which means there should be no connectivity
  // between the two crack surfaces
  // no need to do any calculation
  if( failNorm_ or failTang_ )
    return 0;

  /*******************************************************************/ //blockkk
  if( elevec1 == Teuchos::null )
    return 0;
  /*******************************************************************/

  // initialize this elements if not done before
  if( coheStrength_.size() == 0 )
    InitializeElement();

  /*DRT::ELEMENTS::Dcohesive::ActionType act = Dcohesive::calc_none;
  // get the action required
  std::string action = params.get<std::string>("action","calc_none");
  if (action == "calc_none") dserror("No action supplied");
  else if (action=="calc_struct_linstiff")      act = Dcohesive::calc_struct_linstiff;
  else if (action=="calc_struct_nlnstiff")      act = Dcohesive::calc_struct_nlnstiff;
  else if (action=="calc_struct_internalforce") return 0;//act = Dcohesive::calc_struct_internalforce;
  else if (action=="calc_struct_linstiffmass")  return 0;//act = Dcohesive::calc_struct_linstiffmass;
  else if (action=="calc_struct_nlnstiffmass")  return 0;//act = Dcohesive::calc_struct_nlnstiffmass;
  else if (action=="calc_struct_nlnstifflmass") return 0;//act = Dcohesive::calc_struct_nlnstifflmass;
  else if (action=="calc_struct_stress")        act = Dcohesive::calc_struct_stress;
  else if (action=="calc_struct_eleload")       return 0;//act = Dcohesive::calc_struct_eleload;
  else if (action=="calc_struct_fsiload")       return 0;//act = Dcohesive::calc_struct_fsiload;
  else if (action=="calc_struct_update_istep")  return 0;//act = Dcohesive::calc_struct_update_istep;
  else if (action=="calc_struct_reset_istep")   return 0;//act = Dcohesive::calc_struct_reset_istep;
  else if (action=="calc_struct_ptcstiff")      return 0;//act = Dcohesive::calc_struct_ptcstiff;
  else dserror("Unknown type of action for Dcohesive");*/


  // calling appropriate crack modeling procedure
  switch ( model_ )
  {
  case INPAR::CRACK::crack_dczm:
  {
    // ------------------------- Discrete Cohesive Zone Modeling ----------------------------
    EvaluateDCZM( params, discretization, lm, elemat1, elemat2, elevec1, elevec2, elevec3 );
    break;
  }

  case INPAR::CRACK::crack_ddzm:
  {
    // ------------------------- Discrete Damage Zone Modeling ----------------------------
    EvaluateDDZM( params, discretization, lm, elemat1, elemat2, elevec1, elevec2, elevec3 );
    break;
  }
  default:
    dserror("check the crack model in your input file\n");
    break;
  }

  return 0;
}

int DRT::ELEMENTS::Dcohesive::EvaluateDDZM( Teuchos::ParameterList&    params,
                                            DRT::Discretization&      discretization,
                                            std::vector<int>&         lm,
                                            Epetra_SerialDenseMatrix& elemat1,
                                            Epetra_SerialDenseMatrix& elemat2,
                                            Epetra_SerialDenseVector& elevec1,
                                            Epetra_SerialDenseVector& elevec2,
                                            Epetra_SerialDenseVector& elevec3)
{
  std::cout<<"!--------------in Discrete Damage Zone Modeling---------------------!\n";

  //--------------------------------------------------------------------------
  // STEP 1 : compure required constants (B, initial stiffness, critical force)
  //--------------------------------------------------------------------------
  double b[2], initStiff[2], critForce[2], critDisp[2];
  computeConstants( b, initStiff, critForce, critDisp );

  //--------------------------------------------------------------------------
  // STEP 2 : calculate deflections of the spring elements in local coodinates
  //--------------------------------------------------------------------------
  std::vector<double> deflec = calculateDeflections( discretization, lm );

  //--------------------------------------------------------------------------
  // STEP 3 : compute damage in the spring elements and check if it is failed
  //--------------------------------------------------------------------------
  std::vector<double> forceSpr(2);
  computeDamage( b, initStiff, critForce, critDisp, deflec, forceSpr );
  if( failNorm_ and failTang_ )
    return 0;

  //--------------------------------------------------------------------------
  // STEP 4 : evaluate stiffness matrix
  //--------------------------------------------------------------------------
  elevec1(0) = forceSpr[0];
  elevec1(3) = -forceSpr[0];
  double stiff = (1.0-damage_(0,0))*initStiff[0];
  elemat1(0,0) = elemat1(3,3) = stiff;
  elemat1(0,3) = elemat1(3,0) = -stiff;
  //elemat1.Print(std::cout);
  //elevec1.Print(std::cout);
  //dserror("done");//blockkk

  //--------------------------------------------------------------------------
  // STEP 5 : copy new damage into old time step
  //--------------------------------------------------------------------------
  damagePrev_ = damage_;


  return 0;
}

int DRT::ELEMENTS::Dcohesive::EvaluateDCZM( Teuchos::ParameterList&    params,
                                            DRT::Discretization&      discretization,
                                            std::vector<int>&         lm,
                                            Epetra_SerialDenseMatrix& elemat1,
                                            Epetra_SerialDenseMatrix& elemat2,
                                            Epetra_SerialDenseVector& elevec1,
                                            Epetra_SerialDenseVector& elevec2,
                                            Epetra_SerialDenseVector& elevec3)
{
  //--------------------------------------------------------------------------
  // STEP 1 : calculate deflections of the spring elements in local coodinates
  //--------------------------------------------------------------------------
  std::vector<double> deflec = calculateDeflections( discretization, lm );

  //--------------------------------------------------------------------------
  // STEP 2 : evaluate spring forces and stiffness
  //--------------------------------------------------------------------------
  std::vector<double> forceSpr(2),stifSpr(2);
  ForceStiffnessSpring( deflec, forceSpr, stifSpr );

  if( failNorm_ and failTang_ )
  {
    return 0;
  }

  //--------------------------------------------------------------------------
  // STEP 3 : evaluate stiffness matrix
  //--------------------------------------------------------------------------
  elevec1(0) = forceSpr[0];
  elevec1(3) = -forceSpr[0];
  elemat1(0,0) = stifSpr[0];
  elemat1(0,3) = -stifSpr[0];
  elemat1(3,3) = stifSpr[0];
  elemat1(3,0) = -stifSpr[0];
  elemat1.Print(std::cout);
  elevec1.Print(std::cout);
  //dserror("done");//blockkk

  //double a;
  //std::cin>>a;

  return 0;
}


int DRT::ELEMENTS::Dcohesive::EvaluateNeumann(Teuchos::ParameterList& params,
                                           DRT::Discretization&      discretization,
                                           DRT::Condition&           condition,
                                           std::vector<int>&         lm,
                                           Epetra_SerialDenseVector& elevec1,
                                           Epetra_SerialDenseMatrix* elemat1)
{
  return 0;
}

/*------------------------------------------------------------------------------------------*
 * calculate the deflection of spring in local coordinate of this element         sudhakar 04/13
 *------------------------------------------------------------------------------------------*/
std::vector<double> DRT::ELEMENTS::Dcohesive::calculateDeflections( DRT::Discretization& discret,
                                                                    std::vector<int>& lm )
{
  LINALG::Matrix<3,1> xyzDef;

  /*******************************************************/
  // get nodal coordinates
  /*DRT::Node** n = this->Nodes();
  const double* co1 = n[0]->X();
  const double* co2 = n[1]->X();

  // calculate deflections in x-y-z coordinates

  for( unsigned i=0;i<3;i++ )
    xyzDef[i] = co2[i]-co1[i];*/
  /*******************************************************/

  const int numdf = 3;
  RCP<const Epetra_Vector> disp = discret.GetState("displacement");
  if (disp==Teuchos::null) dserror("Cannot get state vectors 'displacement'");
  std::vector<double> mydisp(lm.size());
  DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);
  for( unsigned i=0;i<3;i++ )
  {
    xyzDef(i,0) = Nodes()[1]->X()[i] + mydisp[numdf+i]-
                  Nodes()[0]->X()[i] + mydisp[i];
  }

  // project them in local coordinates of this element
  //TODO: coordinate transformation is not yet implemented
  std::vector<double>deflec(3);

  //normal deflection = xyzDef.refNormal_
  deflec[0] = 0.0;
  for( unsigned i=0;i<3;i++ )
    deflec[0] += refNormal_(i,0)*xyzDef(i,0);

  //tangential deflection = xyzDef - (xyzDef.refNormal_)refNormal_
  deflec[1] = 0.0;
  for( unsigned i=0;i<3;i++ )
    deflec[1] += pow((xyzDef(i,0)-deflec[0]),2);
  deflec[1] = sqrt(deflec[1]);

  // dummy variable
  deflec[2] = 0.0;
  //deflec[0] = fabs(xyzDef[0]);
  //deflec[1] = fabs(xyzDef[1]);
  //deflec[2] = fabs(xyzDef[2]);

  std::cout<<"deflection = "<<deflec[0]<<"\n";//blockkk

  //std::cout<<deflec[0]<<"\t"<<deflec[1]<<"\t"<<deflec[2]<<"\n";//blockkk
  //dserror("okay");//blockkk

  //std::cout<<"DEFLECTION OF SPRING = "<<deflec[0]<<"\n";//blockkk

  return deflec;
}

// [0] -- corresponding to normal components
// [1] -- corresponding to shear components
void DRT::ELEMENTS::Dcohesive::computeConstants( double * b,
                                                 double * initStiff,
                                                 double * critForce,
                                                 double * critDisp )
{
  // parameters related to mode-I opening and mode II opening)

  switch ( tracSepLaw_ )
  {
  case INPAR::CRACK::linear:
  case INPAR::CRACK::trapezoidal:
  case INPAR::CRACK::sinusoidal:
  case INPAR::CRACK::exponential:
  case INPAR::CRACK::ppr:
    for( unsigned i=0;i<2;i++ )
    {
      critDisp[i] = fracEnergy_[i]/coheStrength_[i]/exp(1);
      b[i] = 1.0/critDisp[i];
      critForce[i] = coheStrength_[i]*area_;
      initStiff[i] = critForce[i]/critDisp[i];

    }
    break;
  }
}

void DRT::ELEMENTS::Dcohesive::computeDamage( double * b,
                                              double * initStiff,
                                              double * critForce,
                                              double * critDisp,
                                              const std::vector<double>& deflec,
                                              std::vector<double>& forceSpr )
{
  double g[2];

  // update damage values and ensure irreversibility
  for( unsigned i=0;i<2;i++ )
  {
    if( deflec[i] > critDisp[i] )
    {
      damage_(i,0) = 1.0-1.0/( exp( b[i]*(deflec[i]-critDisp[i] ) ) );
    }

    // to maintain irreversibility of springs
    // copy old damage values if they are more than present values
    if( damagePrev_(i,0) > damage_(i,0) )
      damage_(i,0) = damagePrev_(i,0);
  }

  // area under the traction-separation curve gives G
  // compure fracture energies (G_I and G_II) with present deflection values
  for( unsigned i=0;i<2;i++ )
  {
    // if displacement is less than critical values, area under the linear part of curve
    if( deflec[i] < critDisp[i] )
    {
      forceSpr[i] = critForce[i]*deflec[i]/critDisp[i];
      g[i] = 0.5*forceSpr[i]*deflec[i];
    }
    // if not total area under linear curve + a part of exponential contribution
    else
    {
      forceSpr[i] = coheStrength_[i]*area_*deflec[i]*exp(1.0-deflec[i]/critDisp[i])/critDisp[i];
      // contribution from linear part of traction-separation curve
      g[i] = 0.5*critForce[i]*critDisp[i];
      // additional contribution from exponential part
     // g1 += ((b[i]*critDisp[i]+1)-(b[i]*deflec[i]+1)*exp(b[i]*(critDisp[i]-deflec[i])))*critForce[i]/(b[i]*b[i]*critDisp[i]);
      g[i] += coheStrength_[i]*area_*(2*critDisp[i] - critDisp[i]*exp(-deflec[i]/critDisp[i])
                                                    - deflec[i]*exp(-deflec[i]/critDisp[i]) );
    }
  }

  // check whether the spring has already failed by now
  double fac = g[0]/fracEnergy_[0] + g[1]/fracEnergy_[1];
  if( fac > 0.0 )
    std::cout<<"the criterion for fracture = "<<fabs(fac)<<"\n";//blockkkk
  if( fabs(0.9-fac) < 1e-15 or (fac-0.9) > 1e-15 )
  {
    failNorm_ = true;
    failTang_ = true;
  }
}

/*--------------------------------------------------------------------------------------*
 * Read all necessary parameters from input file
 * Calculate parameters (ex. max displacement) that do not change with time
 *--------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Dcohesive::InitializeElement()
{
  const Teuchos::ParameterList& params = DRT::Problem::Instance()->CrackParams();
  coheStrength_.push_back(params.get<double>("NORMAL_COHESIVE_STRENGTH"));
  coheStrength_.push_back(params.get<double>("SHEAR_COHESIVE_STRENGTH"));

  fracEnergy_.push_back(params.get<double>("G_I"));
  fracEnergy_.push_back(params.get<double>("G_II"));

  iniSlop_.push_back(params.get<double>("SLOPE_NORMAL"));
  iniSlop_.push_back(params.get<double>("SLOPE_SHEAR"));

  tracSepLaw_ = DRT::INPUT::IntegralValue<INPAR::CRACK::tractionSeparation>(params,"TRACTION_SEPARATION_LAW");
  model_ = DRT::INPUT::IntegralValue<INPAR::CRACK::tractionSeparation>(params,"CRACK_MODEL");

  switch( tracSepLaw_ )
  {

  // linear traction-separation law
  case INPAR::CRACK::linear:
  {
    // area under traction separation law is equivalent to G
    // G = 0.5*maxDisp*cohestrength
    maxDisp_.resize(2);
    maxDisp_[0] = 2.0*fracEnergy_[0]/coheStrength_[0];
    maxDisp_[1] = 2.0*fracEnergy_[1]/coheStrength_[1];

    critDisp_.resize(2);
    critDisp_[0] = coheStrength_[0]/iniSlop_[0];
    critDisp_[1] = coheStrength_[1]/iniSlop_[1];
    break;
  }

  // park-paulino-roesler traction-separation law
  case INPAR::CRACK::ppr:
  {
    alfa_ppr_ = params.get<double>("ALFA_PPR");
    beta_ppr_ = params.get<double>("BETA_PPR");

    //TODO: Find proper definition of initial slope
    iniSlop_[0] = iniSlop_[1] = 0.02;

    // evalute non-dimensional exponenets
    mn_.resize(2);
    mn_[0] = alfa_ppr_*(alfa_ppr_ - 1.0)*iniSlop_[0]*iniSlop_[0]/(1.0-alfa_ppr_*iniSlop_[0]*iniSlop_[0]);
    mn_[1] = beta_ppr_*(beta_ppr_ - 1.0)*iniSlop_[1]*iniSlop_[1]/(1.0-beta_ppr_*iniSlop_[1]*iniSlop_[1]);

    // evaluate maximum displacements
    maxDisp_.resize(2);
    maxDisp_[0] = fracEnergy_[0]/coheStrength_[0]*alfa_ppr_*iniSlop_[0]*pow((1.0-iniSlop_[0]),(alfa_ppr_-1.0))*
                  (alfa_ppr_/mn_[0]+1.0)*pow((alfa_ppr_/mn_[0]*iniSlop_[0]+1.0),(mn_[0]-1.0));
    maxDisp_[1] = fracEnergy_[1]/coheStrength_[1]*beta_ppr_*iniSlop_[1]*pow((1.0-iniSlop_[1]),(beta_ppr_-1.0))*
                      (beta_ppr_/mn_[1]+1.0)*pow((beta_ppr_/mn_[1]*iniSlop_[1]+1.0),(mn_[1]-1.0));

    // evaluate energy constants
    gamma_.resize(2);
    if( fabs(fracEnergy_[0] - fracEnergy_[1]) < 1e-14 )
    {
      gamma_[0] = -1.0*fracEnergy_[0]*pow((alfa_ppr_/mn_[0]),mn_[0]);
      gamma_[1] = pow((beta_ppr_/mn_[1]),mn_[1]);
    }
    else
    {
      double powfac = Macaulay( fracEnergy_[0], fracEnergy_[1] )/(fracEnergy_[0] - fracEnergy_[1]);
      gamma_[0] = pow(-fracEnergy_[0],powfac)*pow((alfa_ppr_/mn_[0]),mn_[0]);
      powfac = Macaulay( fracEnergy_[1], fracEnergy_[0] )/(fracEnergy_[1] - fracEnergy_[0]);
      gamma_[1] = pow(-fracEnergy_[1],powfac)*pow((beta_ppr_/mn_[1]),mn_[1]);
    }
    break;
  }
  default:
  {
    dserror( "Initialization is not implemented for this type of Traction-Separation-Law" );
    break;
  }
  }


  /***********************************************************************************/ //blockkk
  //std::cout<<"cohesive strengths = "<<coheStrength_[0]<<"\t"<<coheStrength_[1]<<"\n";
  //std::cout<<"fracture energies = "<<fracEnergy_[0]<<"\t"<<fracEnergy_[1]<<"\n";
  //std::cout<<"alfa = "<<alfa_ppr_<<"\tbeta = "<<beta_ppr_<<"\n";
  //std::cout<<"lambda = "<<iniSlop_[0]<<"\t"<<iniSlop_[1]<<"\n";
  //std::cout<<"values of m = "<<mn_[0]<<"\n"<<mn_[1]<<"\n";
  //std::cout<<"max disp = "<<maxDisp_[0]<<"\t"<<maxDisp_[1]<<"\n";
  //std::cout<<"crit disp = "<<critDisp_[0]<<"\t"<<critDisp_[1]<<"\n";
  //std::cout<<"gamma = "<<gamma_[0]<<"\t"<<gamma_[1]<<"\n";
  //dserror("check\n");
  /***********************************************************************************/ //blockkk

}

/*--------------------------------------------------------------------------------------*
 * Return <a,b> = (a-b) if (a-b) > 0; otherwise returns a zero
 *--------------------------------------------------------------------------------------*/
double DRT::ELEMENTS::Dcohesive::Macaulay( double a, double b )
{
  double c = a-b;
  if( fabs(c) < 1e-14 or c < 0.0 )
    return 0.0;
  else
    return c;
}


/*--------------------------------------------------------------------------------------*
 * calculates forces and stiffnesses of cohesive elements                   sudhakar 04/13
 * depending on traction-separation law
 *--------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Dcohesive::ForceStiffnessSpring( const std::vector<double>& deflec,
                                                     std::vector<double>& forcSpr,
                                                     std::vector<double>& stifSpr)
{
  // if deflections are more than maximum deflections, spring is failed
  if( deflec[0] > maxDisp_[0] and deflec[1] > maxDisp_[1] )
  {
    failNorm_ = true;
    failTang_ = true;
    return;
  }

  // traction-separation law dictates the form of spring forces and stiffnesses
  switch( tracSepLaw_ )
  {
  case INPAR::CRACK::linear:
  {
    double trac_nor=0.0,trac_shear=0.0;

    if( fabs( deflec[0] - critDisp_[0]) < 1e-8  )       // spring deflection equals critical value
      trac_nor = coheStrength_[0];
    else if( deflec[0] < critDisp_[0] )                 // spring loading regime
      trac_nor = coheStrength_[0]*deflec[0]/critDisp_[0];
    else if( deflec[0] < maxDisp_[0] )                  // spring softening regime
      trac_nor = coheStrength_[0] - coheStrength_[0]* ( deflec[0] - critDisp_[0] ) / ( maxDisp_[0] - critDisp_[0] );
    else                                                // spring already failed
    {
      trac_nor = 0.0;
      failNorm_ = true;
    }

    if( fabs( deflec[1] - critDisp_[1]) < 1e-8  )       // spring deflection equals critical value
      trac_shear = coheStrength_[1];
    else if( deflec[1] < critDisp_[1] )                 // spring loading regime
      trac_shear = coheStrength_[1]*deflec[1]/critDisp_[1];
    else if( deflec[1] < maxDisp_[1] )                  // spring softening regime
      trac_shear = coheStrength_[1] - coheStrength_[1]* ( deflec[1] - critDisp_[1] ) / ( maxDisp_[1] - critDisp_[1] );
    else                                                // spring already failed
    {
      trac_shear = 0.0;
      failTang_ = true;
    }

    forcSpr[0] = trac_nor * area_;
    forcSpr[1] = trac_shear * area_;

    stifSpr[0] = forcSpr[0] / deflec[0];
    stifSpr[1] = forcSpr[1] / deflec[1];

    break;
  }
  case INPAR::CRACK::trapezoidal:
  case INPAR::CRACK::sinusoidal:
  case INPAR::CRACK::exponential:
  {
    dserror("the chosen traction-separation law is not yet implemented\n");
    break;
  }
  //-----------------Park-Paulino-Roesler model----------------------------------//
  /* Park K, Paulino G and Roesler J. A unified potential-based cohesive model of
   * mixed-mode fracture, J. Mech. Phys. Solids, 57 (2009), 891--908      */
  //------------------------------------------------------------------------------//
  case INPAR::CRACK::ppr:
  {
    //double temp11 = (1.0-deflec[0]/maxDisp_[0]);
    //double temp12 = (1.0-fabs(deflec[1])/maxDisp_[1]);

    //double temp21 = mn_[0]/alfa_ppr_+deflec[0]/maxDisp_[0];
    //double temp22 = mn_[1]/beta_ppr_+fabs(deflec[1])/maxDisp_[1];

    //double term1 = gamma_[0]*pow(temp11,alfa_ppr_)*pow(temp21,mn_[0])+Macaulay( fracEnergy_[0], fracEnergy_[1] );
    //double term2 = gamma_[1]*pow(temp12,beta_ppr_)*pow(temp22,mn_[1])+Macaulay( fracEnergy_[1], fracEnergy_[0] );

    double tracN = gamma_[0]/maxDisp_[0];/**(mn_[0]*pow(temp11,alfa_ppr_)*pow(temp21,(mn_[0]-1))
                                         -alfa_ppr_*pow(temp11,(alfa_ppr_-1))*pow(temp21,mn_[0]))*term2;*/

    double tracS = gamma_[1]/maxDisp_[1];/**(mn_[1]*pow(temp12,beta_ppr_)*pow(temp22,(mn_[1]-1))
                                         -beta_ppr_*pow(temp12,(beta_ppr_-1))*pow(temp22,mn_[1]))*term1*deflec[1]/fabs(deflec[1]);*/

    // either spring reaches normal max. displacements
    // or it reached conjugate tangential disp.
    if( tracN < 0.0 or deflec[0] > maxDisp_[0] )
    {
      if ( fabs(tracN) > 1e-14 )
        failNorm_ = true;
      forcSpr[0] = 0.0;
      stifSpr[0] = 0.0;
      std::cout<<"spring failed = "<<tracN<<"\n"<<maxDisp_[0]<<"\n";//blockkk
    }
    else
    {
      forcSpr[0] = tracN * area_;
      /*stifSpr[0] = gamma_[0]/maxDisp_[0]/maxDisp_[0]*(mn_[0]*(mn_[0]-1.0)*pow(temp11,alfa_ppr_)*pow(temp21,(mn_[0]-2.0))
                                                 +alfa_ppr_*(alfa_ppr_-1.0)*pow(temp11,(alfa_ppr_-2.0))*pow(temp21,mn_[0])
                                                 -2.0*alfa_ppr_*mn_[0]*pow(temp11,(alfa_ppr_-1.0))*pow(temp21,(mn_[0]-1.0)));*/
      stifSpr[0] = deflec[0]/forcSpr[0];
    }

    // either spring reaches tangential max. displacements
    // or it reached conjugate normal disp.
    if( tracS < 0.0 or deflec[1] > maxDisp_[1] )
    {
      if ( fabs(tracS) > 1e-14 )
        failTang_ = true;
      forcSpr[1] = 0.0;
      stifSpr[1] = 0.0;
    }
    else
    {
      forcSpr[1] = tracS * area_;
      /*stifSpr[1] = gamma_[1]/maxDisp_[1]/maxDisp_[1]*(mn_[1]*(mn_[1]-1.0)*pow(temp21,beta_ppr_)*pow(temp22,(mn_[1]-2.0))
                                                 +beta_ppr_*(beta_ppr_-1.0)*pow(temp12,(beta_ppr_-2.0))*pow(temp22,mn_[1])
                                                 -2.0*beta_ppr_*mn_[1]*pow(temp12,(beta_ppr_-1.0))*pow(temp22,(mn_[1]-1.0)));*/
      stifSpr[1] = deflec[1]/forcSpr[1];
    }

    /*if( deflec[0] > 0.0 )
    {
      std::cout<<"deflections = "<<deflec[0]<<"\t"<<deflec[1]<<"\t"<<deflec[2]<<"\n";
      std::cout<<"terms = "<<term1<<"\t"<<term2<<"\n";
      std::cout<<"tractions = "<<tracN<<"\t"<<tracS<<"\n";
      std::cout<<"forces = "<<forcSpr[0]<<"\t"<<forcSpr[1]<<"\n";
      std::cout<<"stiffness = "<<stifSpr[0]<<"\t"<<stifSpr[1]<<"\n";
      std::cout<<"critical displacements = "<<iniSlop_[1]*maxDisp_[0]<<"\t"<<iniSlop_[1]*maxDisp_[1]<<"\n";
      dserror("okay");
    }*/
    break;
  }
  }
}
