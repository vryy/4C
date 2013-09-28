/*----------------------------------------------------------------------*/
/*!
\file damage.cpp
\brief Contains the functions to establish local material law /
       stress-strain law for isotropic material following nonlinear isotropic
       von Mises plasticity and a linear elastic material law
       (St.Venant Kirchhoff).

       isotropic hardening
       - describing the nonlinear (piecewise) hardening curve via \f$ sigma_y 
          \ and \f$strainbar_p_ref \ from input file

       geometric linear, for small strains

       example input line:
       MAT 1 MAT_Struct_Damage YOUNG 206.9 NUE 0.29 DENS 0.0 SAMPLENUM 2
         SIGMA_Y 0.45 0.65 EPSBAR_P 0.0 1.0 TOL 1.0e-6

<pre>
Maintainer: Caroline Danowski
            danowski@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15253
</pre>
*/

/*----------------------------------------------------------------------*
 | headers                                                   dano 09/13 |
 *----------------------------------------------------------------------*/
#include "damage.H"
#include "matpar_bundle.H"

#include "../drt_lib/drt_globalproblem.H"


/*----------------------------------------------------------------------*
 | constructor (public)                                      dano 04/11 |
 *----------------------------------------------------------------------*/
MAT::PAR::Damage::Damage(
  Teuchos::RCP<MAT::PAR::Material> matdata
  )
: Parameter(matdata),
  youngs_(matdata->GetDouble("YOUNG")),
  poissonratio_(matdata->GetDouble("NUE")),
  density_(matdata->GetDouble("DENS")),
  sigma_y_(*(matdata->Get<std::vector<double> >("SIGMA_Y"))),
  strainbar_p_ref_(*(matdata->Get<std::vector<double> >("EPSBAR_P"))),
  abstol_(matdata->GetDouble("TOL"))
{
}

/*----------------------------------------------------------------------*
 | is called in Material::Factory from ReadMaterials()       dano 02/12 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<MAT::Material> MAT::PAR::Damage::CreateMaterial()
{
  return Teuchos::rcp(new MAT::Damage(this));
}


MAT::DamageType MAT::DamageType::instance_;


/*----------------------------------------------------------------------*
 | is called in Material::Factory from ReadMaterials()       dano 02/12 |
 *----------------------------------------------------------------------*/
DRT::ParObject* MAT::DamageType::Create( const std::vector<char> & data )
{
  MAT::Damage* plastic = new MAT::Damage();
  plastic->Unpack(data);
  return plastic;
}


/*----------------------------------------------------------------------*
 | constructor (public)                                      dano 04/11 |
 *----------------------------------------------------------------------*/
MAT::Damage::Damage()
: params_(NULL),
  plastic_step_(false)
{
}


/*----------------------------------------------------------------------*
 | copy-constructor (public)                                 dano 04/11 |
 *----------------------------------------------------------------------*/
MAT::Damage::Damage(MAT::PAR::Damage* params)
: params_(params)
{
}


/*----------------------------------------------------------------------*
 | pack (public)                                             dano 04/11 |
 *----------------------------------------------------------------------*/
void MAT::Damage::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm( data );
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);

  // matid
  int matid = -1;
  // in case we are in post-process mode
  if (params_ != NULL) matid = params_->Id();
  AddtoPack(data,matid);

  // pack history data
  int histsize;
  // if material is not initialized, i.e. start simulation, nothing to pack
  if (!Initialized())
  {
    histsize=0;
  }
  else
  {
    // if material is initialized (restart): size equates number of gausspoints
    histsize = strainpllast_->size();
  }
  AddtoPack(data,histsize); // Length of history vector(s)
  for (int var=0; var<histsize; ++var)
  {
    // insert history vectors to AddtoPack
    AddtoPack(data,strainpllast_->at(var));

    AddtoPack(data,strainbarpllast_->at(var));
  }

  AddtoPack(data,plastic_step_);

  return;
} // Pack()


/*----------------------------------------------------------------------*
 | unpack (public)                                           dano 04/11 |
 *----------------------------------------------------------------------*/
void MAT::Damage::Unpack(const std::vector<char>& data)
{
  isinit_=true;
  std::vector<char>::size_type position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position,data,type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");

  // matid and recover params_
  int matid;
  ExtractfromPack(position,data,matid);
  params_ = NULL;
  if (DRT::Problem::Instance()->Materials() != Teuchos::null)
    if (DRT::Problem::Instance()->Materials()->Num() != 0)
    {
      const int probinst = DRT::Problem::Instance()->Materials()->GetReadFromProblem();
      MAT::PAR::Parameter* mat = DRT::Problem::Instance(probinst)->Materials()->ParameterById(matid);
      if (mat->Type() == MaterialType())
        params_ = static_cast<MAT::PAR::Damage*>(mat);
      else
        dserror("Type of parameter material %d does not fit to calling type %d", mat->Type(), MaterialType());
    }

  // history data
  int histsize;
  ExtractfromPack(position,data,histsize);

  // if system is not yet initialized, the history vectors have to be intialized
  if (histsize == 0)
    isinit_ = false;

  // unpack plastic strain vectors
  strainpllast_ = Teuchos::rcp( new std::vector<LINALG::Matrix<NUM_STRESS_3D,1> > );
  strainplcurr_ = Teuchos::rcp( new std::vector<LINALG::Matrix<NUM_STRESS_3D,1> > );

  strainbarpllast_ = Teuchos::rcp( new std::vector<double> );
  strainbarplcurr_ = Teuchos::rcp( new std::vector<double> );

  // initialise
  LINALG::Matrix<NUM_STRESS_3D,1> tmp_vect(true);
  double tmp_scalar = 0.0;
    
  for (int var=0; var<histsize; ++var)
  {
    // vectors of last converged state are unpacked
    ExtractfromPack(position,data,tmp_vect);
    strainpllast_->push_back(tmp_vect);

    // scalar-valued vector of last converged state are unpacked
    ExtractfromPack(position,data,tmp_scalar);
    strainbarpllast_->push_back(tmp_scalar);

    // current vectors have to be initialized
    strainplcurr_->push_back(tmp_vect);

    strainbarplcurr_->push_back(tmp_scalar);
  }

  plastic_step_ = false;
  int plastic_step;
  ExtractfromPack(position,data,plastic_step);

  // if it was already plastic before,
  if (plastic_step != 0) plastic_step_ = true;

  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d",data.size(),position);

  return;

} // Unpack()


/*---------------------------------------------------------------------*
 | initialise / allocate internal stress variables (public)      04/11 |
 *---------------------------------------------------------------------*/
void MAT::Damage::Setup(int numgp, DRT::INPUT::LineDefinition* linedef)
{
  // initialise history variables
  strainpllast_ = Teuchos::rcp(new std::vector<LINALG::Matrix<NUM_STRESS_3D,1> >);
  strainplcurr_ = Teuchos::rcp(new std::vector<LINALG::Matrix<NUM_STRESS_3D,1> >);

  strainbarpllast_ = Teuchos::rcp(new std::vector<double>);
  strainbarplcurr_ = Teuchos::rcp(new std::vector<double>);

  LINALG::Matrix<NUM_STRESS_3D,1> emptymat(true);
  strainpllast_->resize(numgp);
  strainplcurr_->resize(numgp);

  strainbarpllast_->resize(numgp);
  strainbarplcurr_->resize(numgp);

  for (int i=0; i<numgp; i++)
  {
    strainpllast_->at(i) = emptymat;
    strainplcurr_->at(i) = emptymat;

    strainbarpllast_->at(i) = 0.0;
    strainbarplcurr_->at(i) = 0.0;
  }

  isinit_=true;
  return;

}  // Setup()


/*---------------------------------------------------------------------*
 | update internal stress variables (public)                dano 04/11 |
 *---------------------------------------------------------------------*/
void MAT::Damage::Update()
{
  // make current values at time step tlast+1 to values of last step tlast
  strainpllast_ = strainplcurr_;
  
  strainbarpllast_ = strainbarplcurr_;

  // empty vectors of current data
  strainplcurr_ = Teuchos::rcp(new std::vector<LINALG::Matrix<NUM_STRESS_3D,1> >);
  
  strainbarplcurr_ = Teuchos::rcp(new std::vector<double>);

  // get the size of the vector
  // (use the last vector, because it includes latest results, current is empty)
  const int histsize = strainpllast_->size();
  strainplcurr_->resize(histsize);
  
  strainbarplcurr_->resize(histsize);

  const LINALG::Matrix<NUM_STRESS_3D,1> emptyvec(true);
  for (int i=0; i<histsize; i++)
  {
    strainplcurr_->at(i) = emptyvec;
    
    strainbarplcurr_->at(i) = 0.0;
  }

  return;
}  // Update()


/*----------------------------------------------------------------------*
 | evaluate material (public)                                dano 08/11 |
 *----------------------------------------------------------------------*/
void MAT::Damage::Evaluate(
  const LINALG::Matrix<3,3>* defgrd,
  const LINALG::Matrix<NUM_STRESS_3D,1>* linstrain,  // linear strain vector
  Teuchos::ParameterList& params,  // parameter list for communication & HISTORY
  LINALG::Matrix<NUM_STRESS_3D,1>* stress, // 2nd PK-stress
  LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D>* cmat // material stiffness matrix
  )
{
  const int gp = params.get<int>("gp",-1);
  if (gp == -1) dserror("no Gauss point number provided in material");
  const int eleID = params.get<int>("eleID",-1);
  if (eleID == -1) dserror("no element provided in material");
  LINALG::Matrix<MAT::NUM_STRESS_3D,1> plstrain(true);

  // get material parameters
  // Young's modulus
  double young = params_->youngs_;
  // Poisson's ratio
  double nu = params_->poissonratio_;
  
  // initialise material parameters
  // lame constant
  // shear modulus parameter mu == G
  double G = 0.0;
  G = young / ( 2.0 * (1.0 + nu) );
  // bulk modulus kappa = E /( 3 ( 1 - 2 nu) )= lambda + 2/3 * mu
  double kappa = 0.0;
  kappa = young /( 3 * (1 - 2 * nu) );

  // build Cartesian identity 2-tensor I_{AB}
  LINALG::Matrix<NUM_STRESS_3D,1> id2(true);
  for (int i=0; i<3; i++) id2(i) = 1.0;

  // linstrain (in): independent variable passed from the element
  //  strain^p: evolution is determined by the flow rule, history varible
  //  strain^e: definition of additive decomposition:
  //  strain^e = strain - strain^p
  // REMARK: stress-like 6-Voigt vector
  LINALG::Matrix<NUM_STRESS_3D,1> strain(*linstrain);

  //---------------------------------------------------------------------------
  // elastic predictor (trial values)
  //---------------------------------------------------------------------------

  // -------------------------------------------- old plastic strains
  
  // plastic strain vector
  // strain^{p,trial}_{n+1} = strain^p_n
  LINALG::Matrix<NUM_STRESS_3D,1> strain_p(true);
  for (int i=0; i<6; i++)
    strain_p(i,0) = strainpllast_->at(gp)(i,0);

  // get old accumulated/equivalent plastic strain only in case of plastic step
  double strainbar_p = 0.0;
  // accumulated or equivalent plastic strain (scalar-valued)
  // astrain^{p,trial}_{n+1} = astrain^p_n
  strainbar_p = strainbarpllast_->at(gp);
  if (strainbarpllast_->at(gp) < 0.0)
    dserror("accumulated plastic strain has to be equal to or greater than zero!");

  // ----------------------------------------------- physical strains
  // convert engineering shear components into physical components
  // input strain is given in Voigt-notation

  // convert engineering shear component (in) into physical component
  for (int i=3; i<6; ++i) strain(i) /= 2.0;
  for (int i=3; i<6; ++i) strain_p(i) /= 2.0;

  // ----------------------------------------------- elastic trial strain
  // assume load step is elastic
  // strain^{e}_{n+1}
  LINALG::Matrix<NUM_STRESS_3D,1> strain_e(true);

  // strain^{e,trial}_{n+1} = strain_n+1 - strain^p_n
  LINALG::Matrix<NUM_STRESS_3D,1> trialstrain_e(false);
  trialstrain_e.Update( 1.0, strain, 0.0 );
  trialstrain_e.Update( (-1.0), strain_p, 1.0 );

  // volumetric strain
  // trace of strain vector
  double tracestrain = trialstrain_e(0) + trialstrain_e(1) + trialstrain_e(2);
  // volstrain = 1/3 . tr( strain ) . Id
  LINALG::Matrix<NUM_STRESS_3D,1> volumetricstrain(false);
  volumetricstrain.Update( (tracestrain/3.0), id2, 0.0 );

  // deviatoric strain
  // dev^e = strain^e - volstrain^e
  LINALG::Matrix<NUM_STRESS_3D,1> devstrain(false);
  devstrain.Update( 1.0, trialstrain_e, 0.0 );
  devstrain.Update( -1.0, volumetricstrain, 1.0 );

  // ------------------------------------------------------- trial stress
  // pressure = kappa . tr( strain ): saved as scalar
  double p = kappa * tracestrain;

  // deviatoric stress = 2 . G . devstrain
  LINALG::Matrix<NUM_STRESS_3D,1> devstress(false);
  devstress.Update( 2*G, devstrain, 0.0 );
  // be careful for shear stresses (sigma_12)
  // in Voigt-notation the shear strains have to be scaled with 1/2
  // normally done in the material tangent (cf. id4sharp)

  // -------------------------------- elastic trial von Mises effective stress
  // q^(trial)_{n+1} := q(s^(trial)_{n+1}) = \sqrt{ 3 . J2 }
  // J2 = 1/2 (s11^2 + s22^2 + s33^2 + 2 . s12^2 + 2 . s23^2 + 2 . s13^2)
  double J2 = 0.0;
  J2 = 1/2.0 * ( devstress(0)*devstress(0) + devstress(1)*devstress(1)
                 + devstress(2)*devstress(2) ) + devstress(3)*devstress(3) 
       + devstress(4)*devstress(4) + devstress(5)*devstress(5);
                   
  // trial von Mises effective stress
  // q^{trial}_{n+1} := q(s^{trial}_{n+1}) = \sqrt{ 3 . J2 }
  double q = 0.0;
  q = sqrt( 3.0 * J2 );

  //---------------------------------------------------------------------------
  // check plastic admissibility, Phi<=0 is admissble
  //---------------------------------------------------------------------------

  // ------------------------------------------------ trial yield function

  //calculate the current isotropic hardening modulus
  double Hiso = GetIsoHardAtStrainbarnp(strainbar_p);

  // calculate the uniaxial yield stress out of samples
  // TODO extend: 3.argument represents internal hardening variable
  double sigma_y = GetSigmaYAtStrainbarnp(strainbar_p,0);

  // calculate the yield function
  // Phi = \sqrt{ 3.0 . J2 } - sigma_y = q - sigma_y
  // with trial values: Phi_trial = q_trial - sigma_y and Dgamma == 0
  double Phi_trial = 0.0;
  Phi_trial = q - sigma_y;

  // --------------------------------------------------------- initialise

  // if trial state is violated, i.e. it's a plastic load step, there are 2
  // possible states: plastic loading: heaviside = 1, elastic unloading = 0)
  double heaviside = 0.0;
  // incremental plastic multiplier Delta gamma
  double Dgamma = 0.0;

  // unit flow vector Nbar (Prandtl-Reuss)
  // (using the updated relative stress s_n+1, no longer s_n+1^trial)
  // Nbar = ( s^{trial}_{n+1} / || s^{trial}_{n+1} || )
  LINALG::Matrix<NUM_STRESS_3D,1> Nbar(true);

  // flow vector N (Prandtl-Reuss)
  // (using the updated deviatoric stress s_n+1, no longer s_n+1^trial)
  // N = sqrt{3/2} . ( s_{n+1} / || s_{n+1} || )
  LINALG::Matrix<NUM_STRESS_3D,1> N(true);

  //-------------------------------------------------------------------
  // IF consistency condition is violated, i.e. plastic load step
  // ( Phi_trial > 0.0, Dgamma >= 0.0 )
  //-------------------------------------------------------------------
  if (Phi_trial > 1.0e-08)  // if (Phi_trial > 0.0)
  {
    // only first plastic call is output at screen for every processor
    // visualisation of whole plastic behaviour via PLASTIC_STRAIN in postprocessing
    if (plastic_step_ == false)
    {
      if ( (plastic_step_ == false) and (gp == 0) )
        std::cout << "plasticity starts in element = " << eleID << std::endl;

      plastic_step_ = true;
    }

    // --------------------------------------------------------- return-mapping

    // local Newton-Raphson

    // initialise
    const int itermax = 50;  // max. number of iterations
    int itnum = 0;  // iteration counter

    // Res:= residual of Newton iteration == yield function
    double Res = 0.0;
    // calculate residual derivative/tangent
    // ResTan = Phi' = d(Phi)/d(Dgamma)
    double ResTan = 0.0;
    // safety check: set to zero
    Dgamma = 0.0;

    // start iteration with index m for local Newton
    while (true)
    {
      itnum++;
      // check for convergence

      // if not converged m > m_max
      if (itnum > itermax)
      {
        dserror(
          "local Newton iteration did not converge after iteration %3d/%3d with Res=%3d",
          itnum,
          itermax,
          Res
          );
      }
      // else: continue loop m <= m_max

      // Res:= Phi = q^(trial)_{n+1} - sigma_y
      // Res = q - 3 * G * Dgamma - sigma_y;
      // with sigma_y = sigma_y(strainbar_p + Dgamma)
      Res = q - 3 * G * Dgamma - sigma_y;

      // check for convergence
      double norm = abs(Res);
      // check: absolute value of Res has to be smaller than given tolerance
      if (norm < (params_->abstol_))
      {
#ifdef DEBUGMATERIAL
        if (gp == 0)
          printf("Newton method converged after %i iterations; abs(Res)=  %-14.8E\n", itnum, abs(Res));
#endif  // #ifdef DEBUGMATERIAL
        break;
      }

      // plasticity with piecewise linear isotropic hardening
      // ResTan = -3G -Hiso = const.
      ResTan = - 3 * G - Hiso;

      // incremental plastic multiplier Dgamma
      // Dgamma^{m} = Dgamma^{m-1} - Phi / Phi'
      Dgamma += ( -Res )/ResTan;

      // -------------------------- local Newton update of plastic values

      // compute new residual of accumulatd plastic strains
      // astrain^p_{n+1} = astrain^p_n + Dgamma
      // astrain^p_{n+1} = SUM{Dgamma_n} from all time steps n
      // Kuhn-Tucker: Dgamma >= 0.0 --> astrain^p_{n+1} >= 0.0
      strainbar_p = strainbarpllast_->at(gp) + Dgamma;
      if (strainbar_p < 0.0)
      {
        std::cout << "strainbar_p = " << strainbar_p << std::endl;
        dserror("accumulated plastic strain has to be equal or greater than zero");
      }

      // sigma_y = sigma_y(astrain^p_{n+1})
      // TODO damage: 3.argument is R
      sigma_y = GetSigmaYAtStrainbarnp(strainbar_p,0);
      
#ifdef DEBUGMATERIAL
      if (gp == 0)
      {
        std::cout << "am 1.GP: local Newton: Res " << Res << std::endl;
        std::cout << "local Newton: ResTan " << ResTan << std::endl;
        std::cout << "local Newton: Dgamma " << Dgamma << std::endl;
        std::cout << "local Newton: sigma_y " << sigma_y << std::endl;
      }
#endif  // #ifdef DEBUGMATERIAL

    }  // end of local Newton iteration

    // --------------------------------------------------- plastic update

    // ---------------------------------------------- update flow vectors

    
    // deviatoric stress norm || s^{trial}_{n+1} ||
    double devstressnorm = 0.0;
    devstressnorm = sqrt( devstress(0)*devstress(0) + devstress(1)*devstress(1)
                          + devstress(2)*devstress(2) 
                          + 2 * ( devstress(3)*devstress(3) + devstress(4)*devstress(4) 
                                  + devstress(5)*devstress(5) ) );

    // unit flow vector Nbar = s^{trial}_{n+1} / || s^{trial}_{n+1} ||
    Nbar.Update(1.0, devstress, 0.0);
    Nbar.Scale(1 / devstressnorm);

    // flow vector N = sqrt(3/2) . Nbar
    N.Update( (sqrt(3.0 / 2.0)), Nbar, 0.0);

    // deviatoric stress
    // s = s_{n+1}^{trial} - 2 . G . Delta gamma . N
    const double facdevstress = (-2.0) * G * Dgamma;
    devstress.Update( facdevstress, N, 1.0);

    // total stress
    // sigma_{n+1} = s_{n+1} + p_{n+1} . id2
    // pressure/volumetric stress no influence due to plasticity
    Stress( p, devstress, *stress );

    // total strains
    // strain^e_{n+1} = strain^(e,trial)_{n+1} - Dgamma . N
    // compute converged engineering strain components (Voigt-notation)
    strain_e.Update( 1.0, trialstrain_e, 0.0);
    strain_e.Update( (-Dgamma), N, 1.0 );

    // strain^p_{n+1} = strain^p_n + Dgamma . N
    strain_p.Update( Dgamma, N, 1.0 );

    // compute converged engineering strain components (Voigt-notation)
    for (int i=3; i<6; ++i) strain_e(i) *= 2.0;
    for (int i=3; i<6; ++i) strain_p(i) *= 2.0;

    // pass the current plastic strains to the element (for visualisation)
    plstrain.Update(1.0, strain_p, 0.0);

    // --------------------------------------------------- update history
    // plastic strain
    strainplcurr_->at(gp) = strain_p;

    // accumulated plastic strain
    strainbarplcurr_->at(gp) = strainbar_p;

#ifdef DEBUGMATERIAL
    std::cout << "end strain_p\n " << strain_p << std::endl;
    std::cout << "end strainplcurr_->at(gp)\n " << strainplcurr_->at(gp) << std::endl;
#endif //ifdef DEBUGMATERIAL

  }  // plastic corrector

  //------------------------------------------------------------------------
  // ELSE: elastic step (Phi_trial <= 0.0, Dgamma = 0.0)
  //------------------------------------------------------------------------
  else  // (Phi_trial <= 0.0)
  {
    // trial state vectors = result vectors of time step n+1
    // sigma^e_n+1 = sigma^(e,trial)_n+1 = s^{trial}_{n+1} + p . id2
    Stress( p, devstress, *stress );

    // total strains
    // strain^e_{n+1} = strain^(e,trial)_{n+1}
    // compute converged engineering strain components (Voigt-notation)
    strain_e.Update( 1.0, trialstrain_e, 0.0);
    for (int i=3; i<6; ++i) strain_e(i) *= 2.0;

    // no plastic yielding
    Dgamma = 0.0;

    // pass the current plastic strains to the element (for visualisation)
    // compute converged engineering strain components (Voigt-notation)
    for (int i=3; i<6; ++i) strain_p(i) *= 2.0;
    plstrain.Update(1.0, strain_p, 0.0);

    // constant values for
    //  - plastic strains
    //  - accumulated plastic strains
    //  - stress

  }  // elastic step

  //---------------------------------------------------------------------------
  // --------------------------------- consistent elastoplastic tangent modulus
  //---------------------------------------------------------------------------

  // using an associative flow rule: C_ep is symmetric
  // ( generally C_ep is nonsymmetric )
  // if Phi^trial = 0: two tangent stress-strain relations exist
  // plastic loading --> C == C_ep
  if (Dgamma > 0.0)
    heaviside = 1.0;
  // elastic unloading --> C == C_e
  else
    heaviside = 0.0;

  // using an associative flow rule: C_ep is symmetric
  // ( generally C_ep is nonsymmetric )
  SetupCmatElastoPlastic(
    *cmat,
    Dgamma,
    G,
    q,
    Nbar,
    Hiso,
    heaviside
    );

#ifdef DEBUGMATERIAL
  std::cout << "Nach Setup Cep\n" << std::endl;
  std::cout << " Dgamma " << Dgamma << std::endl;
  std::cout << " G " << G << std::endl;
  std::cout << " q " << q << std::endl;
  std::cout << " flow vector " << Nbar << std::endl;
  std::cout << " heaviside " << heaviside << std::endl;
  std::cout << "cmat " << cmat << std::endl;
#endif // #ifdef DEBUGMATERIAL

  // ------------------------------- return plastic strains for post-processing
  params.set<LINALG::Matrix<MAT::NUM_STRESS_3D,1> >("plglstrain",plstrain);

  return;

}  // Evaluate()


/*----------------------------------------------------------------------*
 | computes linear stress tensor                             dano 05/11 |
 *----------------------------------------------------------------------*/
void MAT::Damage::Stress(
  const double p,  // volumetric stress
  const LINALG::Matrix<NUM_STRESS_3D,1>& devstress,  // deviatoric stress tensor
  LINALG::Matrix<NUM_STRESS_3D,1>& stress // 2nd PK-stress
  )
{
  // total stress = deviatoric + hydrostatic pressure . I
  // sigma = s + p . I
  stress.Update(1.0, devstress, 0.0);
  for (int i=0; i<3; ++i) stress(i) += p;

}  // Stress()


/*----------------------------------------------------------------------*
 | computes isotropic elasticity tensor in matrix notion     dano 04/11 |
 | for 3d                                                               |
 *----------------------------------------------------------------------*/
void MAT::Damage::SetupCmat(LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D>& cmat)
{
  // get material parameters
  // Young's modulus (modulus of elasticity)
  double young = params_->youngs_;
  // Poisson's ratio
  double nu = params_->poissonratio_;

  // isotropic elasticity tensor C in Voigt matrix notation, cf. FEscript p.29
  //                       [ 1-nu     nu     nu |          0    0    0 ]
  //                       [        1-nu     nu |          0    0    0 ]
  //           E           [               1-nu |          0    0    0 ]
  //   C = --------------- [ ~~~~   ~~~~   ~~~~   ~~~~~~~~~~  ~~~  ~~~ ]
  //       (1+nu)*(1-2*nu) [                    | (1-2*nu)/2    0    0 ]
  //                       [                    |      (1-2*nu)/2    0 ]
  //                       [ symmetric          |           (1-2*nu)/2 ]
  //
  const double mfac = young / ( (1.0 + nu) * (1.0 - 2.0 * nu) );  // factor

  // clear the material tangent
  cmat.Clear();
  // write non-zero components
  cmat(0,0) = mfac * (1.0 - nu);
  cmat(0,1) = mfac * nu;
  cmat(0,2) = mfac * nu;
  cmat(1,0) = mfac * nu;
  cmat(1,1) = mfac * (1.0 - nu);
  cmat(1,2) = mfac * nu;
  cmat(2,0) = mfac * nu;
  cmat(2,1) = mfac * nu;
  cmat(2,2) = mfac * (1.0 - nu);
  // ~~~
  cmat(3,3) = mfac * 0.5 * (1.0 - 2.0 * nu);
  cmat(4,4) = mfac * 0.5 * (1.0 - 2.0 * nu);
  cmat(5,5) = mfac * 0.5 * (1.0 - 2.0 * nu);

}  // SetupCmat()


/*----------------------------------------------------------------------*
 | computes isotropic elasticity tensor in matrix notion     dano 05/11 |
 | for 3d                                                               |
 *----------------------------------------------------------------------*/
void MAT::Damage::SetupCmatElastoPlastic(
  LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D>& cmat,  // elasto-plastic tangent modulus (out)
  double Dgamma,  // plastic multiplier
  double G,  // shear modulus
  double q,  // elastic trial von Mises effective stress
  LINALG::Matrix<NUM_STRESS_3D,1> Nbar,  // unit flow vector
  double Hiso,  // isotropic hardening modulus
  double heaviside  // Heaviside function
  )
{
  // incremental constitutive function for the stress tensor
  // sigma_n+1 = [ cmat - (Dgamma 6 G^2/q) I_d ] : strain_n+1^{e,trial}
  // consistent tangent operator
  // D^{ep} := dsigma_n+1 / dstrain_n+1^{e,trial}

  // depending on the flow vector Cmat_ep can be a fully-occupied matrix

  // C_ep = C_e - ( H^ . Dgamma . 6 . G^2 ) / q . I_d +
  //        +  H^ . 6 . G^2 ( Dgamma/q - 1/(3 G + Hiso) ) Nbar \otimes Nbar

  // ---------------------------------------------------------- Heaviside
  // if plastic loading:   heaviside = 1.0 --> use C_ep
  // if elastic unloading: heaviside = 0.0 --> use C_e

  // I_d = I_s - 1/3 I . I
  // I_d in Voigt-notation applied to symmetric problem, like stress calculation
  //         [ 2/3   -1/3  -1/3 | 0    0    0  ]
  //         [-1/3    2/3  -1/3 | 0    0    0  ]
  //         [-1/3   -1/3   2/3 | 0    0    0  ]
  //   I_d = [ ~~~~  ~~~~  ~~~~  ~~~  ~~~  ~~~ ]
  //         [                  | 1/2   0   0  ]
  //         [    symmetric     |      1/2  0  ]
  //         [                  |          1/2 ]
  //

  // build Cartesian identity 2-tensor I_{AB}
  LINALG::Matrix<NUM_STRESS_3D,1> id2(true);
  for (int i=0; i<3; i++) id2(i) = 1.0;

  // set Cartesian identity 4-tensor in 6-Voigt matrix notation
  // this is fully 'contra-variant' identity tensor, ie I^{ABCD}
  // REMARK: rows are stress-like 6-Voigt
  //         columns are stress-like 6-Voigt
  LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D> id4sharp(true);
  for (int i=0; i<3; i++) id4sharp(i,i) = 1.0;
  for (int i=3; i<6; i++) id4sharp(i,i) = 0.5;

  // ------------------------------------------------------- elastic term
  // C_ep = C_e
  // add standard isotropic elasticity tensor C_e first
  SetupCmat(cmat);
  
  // ------------------------------------------------------ plastic terms
  
  // if plastic loading:   heaviside = 1.0 --> use C_ep
  // if elastic unloading: heaviside = 0.0 --> use C_e
  
  // ------------------------------------------------- first plastic term
  // - ( H^ . Dgamma . 6 . G^2 ) / q^{trial} . I_d

  double epfac = 0.0;
  double epfac3 = 0.0;
  // elastic trial von Mises effective stress
  if (q != 0.0)
  {
    epfac = (-1.0) * heaviside * Dgamma * 6 * G * G / q;
  }
  // constitutive tensor
  // I_d = id4sharp - 1/3 Id \otimes Id
  // contribution: Id4^#
  cmat.Update(epfac, id4sharp, 1.0);
  // contribution: Id \otimes Id
  double epfac1 = 0.0;
  epfac1 = epfac / (-3.0);
  cmat.MultiplyNT(epfac1, id2, id2, 1.0);

  // ------------------------------------------------ second plastic term

  // loop strains (columns)
  for (int k=0; k<6; ++k)
  {
    // ---------------------------------------------------------- tangent
    // loop stresses (rows)
    for (int i=0; i<6; ++i)
    {
      if (q != 0.0)
      {
        epfac3 =  heaviside * 6 * G * G * ( Dgamma / q - 1.0 / (3 * G + Hiso) );
        cmat(i,k) += epfac3 * Nbar(i) * Nbar(k);
      }  // (q != 0.0)
    }  // end rows, loop i
  }  // end columns, loop k

  // complete material tangent C_ep available

#ifdef DEBUGMATERIAL
  if (Dgamma != 0)
  {
    std::cout << "Ende SetupCmatElastPlast" << std::endl;
    std::cout << "Cep\n" << " Dgamma " << Dgamma << std::endl;
    std::cout << " G " << G << std::endl;
    std::cout << " q " << q << std::endl;
    std::cout << " Nbar " << Nbar << std::endl;
    std::cout << " heaviside " << heaviside << std::endl;
    std::cout << " epfac " << epfac << std::endl;
    std::cout << " epfac1 " << epfac1 << std::endl;
    std::cout << " cmat " << cmat << std::endl;
  }
#endif // #ifdef DEBUGMATERIAL

}  // SetupCmatElastoPlastic()


/*----------------------------------------------------------------------*
 | return derivative of piecewise linear function for the    dano 09/13 |
 | yield stress, i.e. isotropic hardening modulus at current            |
 | accumulated plastic strain                                           |
 *----------------------------------------------------------------------*/
double MAT::Damage::GetIsoHardAtStrainbarnp(
  const double strainbar_p  // current accumulated strain
  )
{
  // Hiso = d sigma_y / d astrain^p_{n+1}
  double Hiso = 0.0;

  // extract vectors of samples
  const std::vector<double> strainbar_p_ref = params_->strainbar_p_ref_;
  const std::vector<double> sigma_y_ref = params_->sigma_y_;
  // how many samples are available
  double samplenumber = sigma_y_ref.size();

  // loop over all samples
  for (int i=1; i<samplenumber; ++i)
  {
    // astrain^{p}_{n+1} > astrain^{p}_ref^[i]
    if (strainbar_p >= strainbar_p_ref[i])
    {
      // astrain^{p}_{n+1} > astrain^{p}_ref^[i] --> sigma_y = sigma_ref[i]
      // --> Hiso = d sigma_y / d astrain^{p}_{n+1} = 0
      Hiso = 0.0;
      break;
    }

    // (strainbar_p < strainbar_p_ref[i])
    else
    {
      // load is still elastic: astrain^{p}_{n+1} < astrain^{p}_ref^{i=0}
      if (i == 0)
      {
        // yield boundary is the initial yield stress (sigma_y^{i=0})
        Hiso = 0.0;
        break;
      }
      // astrain^{p,i-1}_ref < astrain^{p}_{n+1} < astrain^{p,i}_ref
      else
      {
        //         sigma_y_n - sigma_y^{i-1}
        // Hiso =  ---------------------------------------
        //        astrain^{p,i}_ref - astrain^{p,i-1}_ref
        Hiso = (sigma_y_ref[i] - sigma_y_ref[i-1]) / (strainbar_p_ref[i] - strainbar_p_ref[i-1]);
        break;
      }
    }  // load is plastic, hardening can occur
  }  // loop over samples

  // return current isotropic hardening modulus
  return Hiso;

}  // GetIsoHardeningModulus()


/*----------------------------------------------------------------------*
 | compute current yield stress sigma_y(astrain^p)           dano 09/13 |
 | calculate yield stress from (sigma_y-astrain^p)-samples              |
 *----------------------------------------------------------------------*/
double MAT::Damage::GetSigmaYAtStrainbarnp(
  const double strainbar_p,  // current accumulated strain, in case of dependent hardening
  const double R  // isotropic hardening internal variable
  )
{
  // extract vectors of samples
  const std::vector<double> sigma_y_ref = params_->sigma_y_;
  // how many samples are available
  double samplenumber = sigma_y_ref.size();
  // return the yield stress
  double sigma_y_interpol = 0.0;

  // uniaxial yield stress given by piecewise linear curve
  if (samplenumber > 0)
  {
    // get vector astrain^p_ref
    const std::vector<double> strainbar_p_ref = params_->strainbar_p_ref_;
    if (sigma_y_ref.size() != strainbar_p_ref.size())
      dserror("Samples have to fit to each other!");

    // loop over samples
    for (int i=0; i<samplenumber; ++i)
    {
      // astrain^{p}_{n+1} > astrain^{p}_ref^[i]
      if (strainbar_p >= strainbar_p_ref[i])
      {
        sigma_y_interpol = sigma_y_ref[i];
      }
      // current strains are <= strainbar_p_ref_max
      else  // astrain^{p}_{n+1} < astrain^{p}_ref^[i]
      {
        // astrain^{p}_{n+1} < astrain^{p}_ref^{i=0}, i.e. load is still elastic
        if (i == 0)
        {
          // yield boundary is the initial yield stress (sigma_y^{i=0})
          sigma_y_interpol = sigma_y_ref[0];
          break;
        }
        // astrain^{p}_ref^{i=0} < astrain^{p}_{n+1} < astrain^{p}_ref^i
        else
        {
          // astrain^{p,i-1}_ref < astrain^{p}_{n+1} < astrain^{p,i}_ref
          if (strainbar_p < strainbar_p_ref[i])
          {
            // sigma_y_{n+1} = sigma_y^i +
            //                                        sigma_y^i - sigma_y^{i-1}
            // + (astrain^p_{n+1} - astrain^{p,i-1}) ---------------------------------------
            //                                      astrain^{p,i}_ref - astrain^{p,i-1}_ref
            sigma_y_interpol = sigma_y_ref[i-1] + (strainbar_p - strainbar_p_ref[i-1])
                      * (sigma_y_ref[i] - sigma_y_ref[i-1])
                      / (strainbar_p_ref[i] - strainbar_p_ref[i-1]);
          }  // current strains between strain^{i-1} and strain^i
        }  // plastic regime
      }
    }  // loop over all samples
  }  // samplenumber > 1

  // return current yield stress
  return sigma_y_interpol;

}  // GetSigmaYAtStrainbarnp()


/*---------------------------------------------------------------------*
 | return names of visualization data (public)              dano 09/13 |
 *---------------------------------------------------------------------*/
void MAT::Damage::VisNames(std::map<std::string,int>& names)
{
  std::string accumulatedstrain = "accumulatedstrain";
  names[accumulatedstrain] = 1; // scalar
}  // VisNames()

/*---------------------------------------------------------------------*
 | return visualization data (public)                       dano 09/13 |
 *---------------------------------------------------------------------*/
bool MAT::Damage::VisData(const std::string& name, std::vector<double>& data, int numgp)
{
  if (name == "accumulatedstrain")
  {
    if ((int)data.size()!=1) dserror("size mismatch");
    double temp = 0.0;
    for (int iter=0; iter<numgp; iter++)
      temp += AccumulatedStrain(iter);
    data[0] = temp/numgp;
  }
  return true;
}  // VisData()


/*----------------------------------------------------------------------*/

