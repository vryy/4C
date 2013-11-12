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

       geometric linear, for small strains including isotropic ductile damage
       - simplified Lemaitre model only considers isotropic hardening
       
       ductile isotropic damage
       - elasticity-damage coupling, cf. de Souza Neto Compuatational plasticity Chapt. 12
       - Lemaitre's elastoplastic damage theory (simplified version without 
         kinematic hardening)
       - isotropic hardening internal variable \f$ R\
       - damage variable \f$ D\
       
       ductile isotropic damage according to the book of de Souza Neto et al 
       "Computational methods for plasticity", chapter 12

       example input line:
       MAT 1 MAT_Struct_Damage YOUNG 206.9 NUE 0.29 DENS 0.0 SAMPLENUM 2
         SIGMA_Y 0.45 0.65 EPSBAR_P 0.0 1.0 DAMDEN 0.0035 DAMEXP 1.0
         DAMTHRESHOLD 1.0e-06  TOL 1.0e-6

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
  damden_(matdata->GetDouble("DAMDEN")),
  damexp_(matdata->GetDouble("DAMEXP")),
  epsbarD_(matdata->GetDouble("DAMTHRESHOLD")),
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
    AddtoPack(data,isohardvarlast_->at(var));
    AddtoPack(data,damagelast_->at(var));
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
  
  // unpack isotropic hardening variable
  // isohardvarlast_ == isohardvarlast_
  // isohardvarcurr_ == isohardvarcurr_
  isohardvarlast_ = Teuchos::rcp( new std::vector<double> );
  isohardvarcurr_ = Teuchos::rcp( new std::vector<double> );

  // unpack isotropic damage internal state variable
  damagelast_ = Teuchos::rcp( new std::vector<double> );
  damagecurr_ = Teuchos::rcp( new std::vector<double> );

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
    ExtractfromPack(position,data,tmp_scalar);
    isohardvarlast_->push_back(tmp_scalar);
    ExtractfromPack(position,data,tmp_scalar);
    damagelast_->push_back(tmp_scalar);

    // current vectors have to be initialised
    strainplcurr_->push_back(tmp_vect);

    strainbarplcurr_->push_back(tmp_scalar);
    isohardvarcurr_->push_back(tmp_scalar);
    damagecurr_->push_back(tmp_scalar);
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
  strainpllast_ = Teuchos::rcp( new std::vector<LINALG::Matrix<NUM_STRESS_3D,1> > );
  strainplcurr_ = Teuchos::rcp( new std::vector<LINALG::Matrix<NUM_STRESS_3D,1> > );

  strainbarpllast_ = Teuchos::rcp( new std::vector<double> );
  strainbarplcurr_ = Teuchos::rcp( new std::vector<double> );

  isohardvarlast_ = Teuchos::rcp( new std::vector<double> );
  isohardvarcurr_ = Teuchos::rcp( new std::vector<double> );

  damagelast_ = Teuchos::rcp( new std::vector<double> );
  damagecurr_ = Teuchos::rcp( new std::vector<double> );

  // set all history variables to zero
  LINALG::Matrix<NUM_STRESS_3D,1> emptymat(true);
  strainpllast_->resize(numgp);
  strainplcurr_->resize(numgp);

  strainbarpllast_->resize(numgp);
  strainbarplcurr_->resize(numgp);

  isohardvarlast_->resize(numgp);
  isohardvarcurr_->resize(numgp);
  
  damagelast_->resize(numgp);
  damagecurr_->resize(numgp);

  for (int i=0; i<numgp; i++)
  {
    strainpllast_->at(i) = emptymat;
    strainplcurr_->at(i) = emptymat;

    strainbarpllast_->at(i) = 0.0;
    strainbarplcurr_->at(i) = 0.0;

    isohardvarlast_->at(i) = 0.0;
    isohardvarcurr_->at(i) = 0.0;

    damagelast_->at(i) = 0.0;
    damagecurr_->at(i) = 0.0;
  }

  isinit_ = true;
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
  isohardvarlast_ = isohardvarcurr_;
  damagelast_ = damagecurr_;

  // empty vectors of current data
  strainplcurr_ = Teuchos::rcp( new std::vector<LINALG::Matrix<NUM_STRESS_3D,1> > );
  
  strainbarplcurr_ = Teuchos::rcp( new std::vector<double> );
  isohardvarcurr_ = Teuchos::rcp( new std::vector<double> );
  damagecurr_ = Teuchos::rcp( new std::vector<double> );

  // get the size of the vector
  // (use the last vector, because it includes latest results, current is empty)
  const int histsize = strainpllast_->size();
  strainplcurr_->resize(histsize);

  strainbarplcurr_->resize(histsize);
  isohardvarcurr_->resize(histsize);
  damagecurr_->resize(histsize);

  const LINALG::Matrix<NUM_STRESS_3D,1> emptyvec(true);
  for (int i=0; i<histsize; i++)
  {
    strainplcurr_->at(i) = emptyvec;

    strainbarplcurr_->at(i) = 0.0;
    isohardvarcurr_->at(i) = 0.0;
    damagecurr_->at(i) = 0.0;
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
  // damage evolution law denominator r
  double damden = params_->damden_;
  // damage evolution law exponent s
  double damexp = params_->damexp_;
  // damage threshold
  double strainbar_p_D = params_->epsbarD_;

  // calculate some more paramters
  // lame constant / shear modulus parameter mu == G
  double G = 0.0;
  G = young / ( 2.0 * (1.0 + nu) );
  // bulk modulus kappa = E /( 3 ( 1 - 2 nu) )= lambda + 2/3 * mu
  double kappa = 0.0;
  kappa = young / ( 3 * (1 - 2 * nu) );

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
  // elastic predictor (trial values: consider old damage and/or plasticity)
  //---------------------------------------------------------------------------

  // -------------------------------- old damage internal state variables

  // to consider damage in the present material model we introduce the
  // '~'-operator to indicate UNDAMAGED values (e.g. s^{~,trial}_{n+1}) in
  // contrast to DAMAGED values (e.g. s^{trial}_{n+1})

  // set trial damage variable to old one
  // D^{trial}_{n+1} = D_n
  double damage = 0.0;
  damage = damagelast_->at(gp);

    // get old damaged isotropic hardening variable only in case of plastic step
  double Rplast = 0.0;
  // damaged accumulated or equivalent plastic strain (scalar-valued)
  // R^{trial}_{n+1} = R_n
  Rplast = isohardvarlast_->at(gp);
  if (isohardvarlast_->at(gp) < 0.0)
  {
    std::cout << "Rplast am ele = " << eleID << ": " << Rplast << std::endl;
    dserror("damaged isotropic hardening variable has to be equal to or greater than zero!");
  }

  // get old integrity: omega_n = 1 - D_n
  double omegaold = 0.0;
  omegaold = 1.0 - damage;

  // ------------------------------------------------ old plastic strains
  
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

  // --------------------------------------------------- physical strains
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
  // devstrain^e = strain^e - volstrain^e
  LINALG::Matrix<NUM_STRESS_3D,1> devstrain(false);
  devstrain.Update( 1.0, trialstrain_e, 0.0 );
  devstrain.Update( -1.0, volumetricstrain, 1.0 );

  // --------------------------------------------- trial undamaged stress

  // undamaged scalar-valued pressure
  // p^{~} = kappa . tr( strain )
  double p_tilde = kappa * tracestrain;

  // deviatoric stress^{~} = 2 . G . devstrain
  LINALG::Matrix<NUM_STRESS_3D,1> devstress_tilde(false);
  devstress_tilde.Update( 2*G, devstrain, 0.0 );
  // be careful for shear stresses (sigma_12)
  // in Voigt-notation the shear strains have to be scaled with 1/2
  // normally done in the material tangent (cf. id4sharp)

  // --------------- trial (undamaged) elastic von Mises effective stress
  
  // q^{~,trial}_{n+1} := q(s^{trial}_{n+1}) / (1-D_n) = \sqrt{ 3 . J2 } / (1-D_n)
  //                    = sqrt{3/2} . || s^{trial}_{n+1} || / (1-D_n)
  // J2 = 1/2 (s11^2 + s22^2 + s33^2 + 2 . s12^2 + 2 . s23^2 + 2 . s13^2)
  double J2 = 0.0;
  J2 = 1/2.0 * ( devstress_tilde(0)*devstress_tilde(0)
                 + devstress_tilde(1)*devstress_tilde(1)
                 + devstress_tilde(2)*devstress_tilde(2) )
       + devstress_tilde(3)*devstress_tilde(3) + devstress_tilde(4)*devstress_tilde(4)
       + devstress_tilde(5)*devstress_tilde(5);
  double q_tilde = 0.0;
  q_tilde = sqrt( 3.0 * J2 );

  // initialise final (damaged) deviatoric stresses
  LINALG::Matrix<NUM_STRESS_3D,1> devstress(true);

  //---------------------------------------------------------------------------
  // check plastic admissibility, Phi<=0 is admissble
  //---------------------------------------------------------------------------

  // ----------------------------------------------- trial yield function
  
  // initialise yield stress and its derivative
  double Hiso = 0.0;
  double sigma_y = 0.0;

  // --------------------------------------------------- damage threshold

  // calculate damage threshold (according to de Souza Neto, p.483ff)

  // bool which decide if threshold is passed or not
  // current strainbar_p < strainbar_p_D
  if (strainbar_p < strainbar_p_D)
  {
    // --> no damage evolution: damage_{n+1} = damage_n == 0
    // strainbar_p_{n+1}= strainbar_p_n + Dgamma

    // no damage evolution -> omega = 1 -D = 1 - 0 = 1
    damage = 0.0;
    omegaold = 1.0;

    // calculate the isotropic hardening modulus with old plastic strains
    // Hiso = dsigma_y / d astrain^p
    Hiso = GetIsoHardAtStrainbarnp(strainbar_p);

    // calculate the uniaxial yield stress out of samples
    sigma_y = GetSigmaYAtStrainbarnp(strainbar_p);
  }
  else  // current strainbar_p > strainbar_p_D
  {
    // calculate the uniaxial yield stress out of samples
    sigma_y = GetSigmaYAtStrainbarnp(Rplast);
  }

  // calculate the yield function
  // Phi = \sqrt{ 3.0 . J2 } - sigma_y = q - sigma_y
  // with trial values: Phi_trial = q{~,trial} - sigma_y and Dgamma == 0
  double Phi_trial = 0.0;
  Phi_trial = q_tilde - sigma_y;

  // --------------------------------------------------------- initialise

  // if trial state is violated, i.e. it's a plastic load step, there are 2
  // possible states: plastic loading: heaviside = 1, elastic unloading = 0)
  double heaviside = 0.0;
  // incremental plastic multiplier Delta gamma
  double Dgamma = 0.0;
  // damage energy release rate Y
  double energyrelrate = 0.0;
  // calculate derivative of engergy release rate Ytan w.r.t. Dgamma
  double Ytan = 0.0;
  // integrity
  // omega_{n+1} = 1 - D_{n+1}
  double omega = 1.0;  // if not damaged, omega == 1.0
  // flag indicating if damage evolution takes place or not
  bool damevolution = false;

  // unit flow vector Nbar (Prandtl-Reuss)
  // (using s_n+1^trial for undamaged, and s_n+1 for damaged load step)
  // Nbar = ( s^{trial}_{n+1} / || s^{trial}_{n+1} || )
  LINALG::Matrix<NUM_STRESS_3D,1> Nbar(true);

  // flow vector N (Prandtl-Reuss)
  // (using the updated deviatoric stress s_n+1, no longer s_n+1^trial)
  // N = sqrt{3/2}/(1 - D_{n+1}) . ( s_{n+1} / || s_{n+1} || )
  LINALG::Matrix<NUM_STRESS_3D,1> N(true);

  //---------------------------------------------------------------------------
  // IF consistency condition is violated, i.e. plastic load step use return-mapping
  // (Phi_trial > 0.0, Dgamma >= 0.0)
  //---------------------------------------------------------------------------
  if (Phi_trial > 1.0e-08)  // if (Phi_trial > 0.0)
  {
    // ------------------------------------------------ damage threshold

    // calculate damage threshold (according to de Souza Neto, p.483ff)

    // bool which decide if threshold is passed or not
    // current strainbar_p < strainbar_p_D
    if (strainbar_p < strainbar_p_D)
    {
      // --> no damage evolution: damage_{n+1} = damage_n == 0
      // strainbar_p_{n+1}= strainbar_p_n + Dgamma

      // no damage evolution -> omega = 1 -D = 1 - 0 = 1
      omega = 1.0;

      // ------------------------------------------------------- return-mapping

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
        Res = q_tilde - 3 * G * Dgamma - sigma_y;

        // check for convergence
        double norm = abs(Res);
        // check: absolute value of Res has to be smaller than given tolerance
        if (norm < (params_->abstol_))
        {
#ifdef DEBUGMATERIAL
          if (gp == 0)
            printf(
              "Newton method converged after %i iterations; abs(Res)=  %-14.8E\n",
              itnum,
              abs(Res)
              );
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
        strainbar_p = strainbarpllast_->at(gp) + Dgamma/omega;
        if (strainbar_p < 0.0)
        {
          std::cout << "strainbar_p = " << strainbar_p << std::endl;
          dserror("accumulated plastic strain has to be equal or greater than zero");
        }
        
        // isotropic damage variable remains the same using D == 0 (omega=1)
        Rplast = isohardvarlast_->at(gp) + Dgamma;

        // Hiso = dsigma_y / d astrain^p_{n+1}
        Hiso = GetIsoHardAtStrainbarnp(strainbar_p);

        // sigma_y = sigma_y(astrain^p_{n+1})
        sigma_y = GetSigmaYAtStrainbarnp(strainbar_p);

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
      
      // check if newest astrain^{p,m}_{n+1} is still smaller than threshold
      if (strainbar_p < strainbar_p_D)
      {
#ifdef DEBUGMATERIAL
        if (gp == 0)
          std::cout << "No damage! Current load step is admissible!" << std::endl;
#endif  // #ifdef DEBUGMATERIAL

        // no damage occurs, i.e. current solution is correct, no return-map
        // considering damage is necessary
        damevolution = false;

        // --------------------------------------------------- plastic update

        // ---------------------------------------------- update flow vectors

        // deviatoric stress norm || s^{trial}_{n+1} ||
        double devstress_tildenorm = 0.0;
        devstress_tildenorm = sqrt( devstress_tilde(0)*devstress_tilde(0)
                                    + devstress_tilde(1)*devstress_tilde(1)
                                    + devstress_tilde(2)*devstress_tilde(2)
                              + 2 * ( devstress_tilde(3)*devstress_tilde(3)
                                      + devstress_tilde(4)*devstress_tilde(4)
                                      + devstress_tilde(5)*devstress_tilde(5) ) );

        // unit flow vector Nbar = s^{trial}_{n+1} / || s^{trial}_{n+1} ||
        Nbar.Update(1.0, devstress_tilde, 0.0);
        Nbar.Scale(1 / devstress_tildenorm);

        // flow vector N = sqrt(3/2) . Nbar
        N.Update( (sqrt(3.0 / 2.0)), Nbar, 0.0);

        // deviatoric stress
        // s = s_{n+1}^{trial} - 2 . G . Delta gamma . N
        const double facdevstress = (-2.0) * G * Dgamma;
        devstress.Update( 1.0, devstress_tilde, 0.0);
        devstress.Update( facdevstress, N, omega);

        // total stress
        // sigma_{n+1} = s_{n+1} + p_{n+1} . id2
        // pressure/volumetric stress no influence due to plasticity
        Stress( p_tilde, devstress, *stress );

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

        // update damaged isotropic hardening variable R_{n+1}
        isohardvarcurr_->at(gp) = Rplast;

        // update damage variable damage_{n+1}
        damagecurr_->at(gp) = damage;

#ifdef DEBUGMATERIAL
        std::cout << "end strain_p\n " << strain_p << std::endl;
        std::cout << "end strainplcurr_->at(gp)\n " << strainplcurr_->at(gp) << std::endl;
#endif //ifdef DEBUGMATERIAL

      }  // (strainbar_p < strainbar_p_D)
      // updated strainbar_p not valid, recalculate step considering damage
      else  // (strainbar_p > strainbar_p_D)
      {
#ifdef DEBUGMATERIAL
        if (gp == 0)
          std::cout << "New solution strainbar_p^m exceeds damage threshold!"
            "\n Recalculate load step considering damage!" << std::endl;
#endif  // #ifdef DEBUGMATERIAL

        // damage occurs, i.e. current solution is not correct, recalculate
        // load step considering damage
        damevolution = true;
      }  // (strainbar_p > strainbar_p_D)
    }  // no damage: (strainbar_p < strainbar_p_D)

    // ------------------------------------ threshold of trial step is exceeded
    // (strainbar_p > strainbar_p_D), i.e. damage evolves
    else
      damevolution = true;

    //-------------------------------------------------------------------
    // -------------------------------------- return-mapping considering damage
    // damage has to be considered (strainbar_p > strainbar_p_D)
    //-------------------------------------------------------------------
    if (damevolution == true)
    {
      // only first plastic call is output at screen for every processor
      // visualisation of whole plastic behaviour via PLASTIC_STRAIN in postprocessing
      if (plastic_step_ == false)
      {
        if ( (plastic_step_ == false) and (gp == 0) )
          std::cout << "damage starts to evolve in element = " << eleID << std::endl;

        plastic_step_ = true;
      }

#ifdef DEBUGMATERIAL
        std::cout << "Damage has to be considered for current load step and ele = "
          << eleID << ", and gp = " << gp << " ! Threshold exceeded!" << std::endl;
#endif  // #ifdef DEBUGMATERIAL

      // ------------------------------------------------- return-mapping

      // local Newton-Raphson

      // ------------------------------- initial guess for Dgamma (12.49)

      // instead of initial guess Dgamma^{m=0} = 0.0 use perfectly plastic
      // solution with frozen yield surface at beginning of each load increment
      //
      // Dgamma^{m=0} = q^{~,trial} - sigma(Rplast_n) * omega_n / 3G
      Dgamma = omegaold * Phi_trial / (3 * G);

      // -------------------- initialise internal variables with Dgamma^0

      // Rplast = R^{p,m=0}_{n+1} = R^p_n + Dgamma
      Rplast = isohardvarlast_->at(gp) + Dgamma;

      const int itermax = 50;  // max. number of iterations
      int itnum = 0;  // iteration counter

      // Res:= residual of Newton iteration == yield function
      double Res = 0.0;
      // calculate residual derivative/tangent
      // ResTan = Phi' = d(Phi)/d(Dgamma)
      double ResTan = 0.0;

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

        // calculate the uniaxial yield stress out of samples using newest
        // solution of Dgamma^m
        // if m=0: sigma_y = sigma_y(R^{p,m}_{n+1}(Dgamma^{m=0}))
        // if m>0: R^{p,m} was updated at end of last local Newton loop
        // sigma_y = sigma_y(R^{p,m}_{n+1})
        sigma_y = GetSigmaYAtStrainbarnp(Rplast);
        // slope of hardening function
        // Hiso = dsigma_y / d Rplast^m_{n+1}
        Hiso = GetIsoHardAtStrainbarnp(Rplast);
        // plasticity with nonlinear (piecewise linear) isotropic hardening

        // get derivative of energy release rate w.r.t. Dgamma
        // d(-energyrelrate) / dDgamma = - Hiso(Rplast^m) . sigma_y(Rplast^m) / (3 . G)
        Ytan = - Hiso * sigma_y / (3 * G);

        // integrity
        // omega_{n+1} = 3G / (q_tilde - sigma_y) * Dgamma = 1 - D_{n+1}
        omega = 3 * G / (q_tilde - sigma_y) * Dgamma;

        // damage energy release rate only implicitely depending on Dgamma (12.47)
        energyrelrate = - (sigma_y * sigma_y) / (6 * G) - p_tilde * p_tilde / (2 * kappa);

        // compute residual function (12.48)
        // Res := F(Dgamma) = omega(Dgamma) - omega_n
        //                    + Dgamma / omega(Dgamma) . (-Y(Dgamma)/r)^s
        //                      . (q_tilde - sigma_y) / (3 G)
        // here: it is important NOT to use Dgamma^{m=0}=0 --> omega=0 --> '1/0'
        double Res = omega - omegaold + pow( (-energyrelrate/damden), damexp )
                                          / ( (3 * G)/(q_tilde - sigma_y) );

        // check for convergence
        double norm = abs(Res);
        // check: absolute value of Res has to be smaller than given tolerance
        if (norm < (params_->abstol_))
        {
#ifdef DEBUGMATERIAL
          if (gp == 0)
            printf(
              "Newton method converged after %i iterations; abs(Res)=  %-14.8E\n",
              itnum,
              abs(Res)
              );
#endif  // #ifdef DEBUGMATERIAL
          break;
        }

        // if load state is not converged, calculate derivatives w.r.t. Dgamma

        // derviative of residual w.r.t. Dgamma
        // ResTan = (3 . G)/(q_tilde - sigma_y)
        //          + (3 . G)/(q_tilde - sigma_y) . Dgamma . Hiso / (q_tilde - sigma_y)
        //          - Hiso / (3 . G) . (-Y/r)^s
        //          - s . Ytan / ( ( (3 . G)/(q_tilde - sigma_y) ) . r) . (-Y/r)^(s-1)
        ResTan = (3 * G) / (q_tilde - sigma_y)
                 + (3 * G)/(q_tilde - sigma_y) * Dgamma * Hiso / (q_tilde - sigma_y)
                 - Hiso / (3 * G) * pow( (-energyrelrate / damden), damexp)
                 - damexp * Ytan / ( ( (3 * G)/(q_tilde - sigma_y) ) * damden )
                   * pow( (-energyrelrate / damden), (damexp-1) );

        // incremental plastic multiplier Dgamma
        // Dgamma^{m} = Dgamma^{m-1} - Phi / Phi'
        Dgamma += ( -Res )/ResTan;

        // -------------------------- local Newton update of plastic values

        // compute new residual of damaged accumulatd plastic strains
        // R^p_{n+1} = R^p_n + Dgamma
        // R^p_{n+1} = SUM{Dgamma_n} from all time steps n
        // Kuhn-Tucker: Dgamma >= 0.0 --> R^p_{n+1} >= 0.0
        Rplast = isohardvarlast_->at(gp) + Dgamma;

        // compute new residual of accumulatd plastic strains
        // astrain^p_{n+1} = astrain^p_n + Dgamma/omega
        // astrain^p_{n+1} = SUM{Dgamma_n} from all time steps n
        // Kuhn-Tucker: Dgamma >= 0.0 --> astrain^p_{n+1} >= 0.0
        strainbar_p = strainbarpllast_->at(gp) + Dgamma/omega;
        if (strainbar_p < 0.0)
        {
          std::cout << "strainbarpllast_->at(gp) = " << strainbarpllast_->at(gp) << std::endl;
          std::cout << "omega = " << omega << std::endl;
          std::cout << "Dgamma = " << Dgamma << std::endl;
          std::cout << "strainbar_p = " << strainbar_p << std::endl;
          dserror("accumulated plastic strain has to be equal or greater than zero");
        }

#ifdef DEBUGMATERIAL
        if (gp == 0)
        {
          std::cout << "Ende local Newton damage = " << damage << std::endl;
          std::cout << "Ende local Newton strainbar_p = " << strainbar_p << std::endl;
          std::cout << "Ende local Newton Rplast = " << Rplast << std::endl;

          std::cout << "am 1.GP: local Newton: Res " << Res << std::endl;
          std::cout << "local Newton: ResTan " << ResTan << std::endl;
          std::cout << "local Newton: Dgamma " << Dgamma << std::endl;
          std::cout << "local Newton: sigma_y " << sigma_y << std::endl;
        }
#endif  // #ifdef DEBUGMATERIAL

      }  // end of local Newton iteration

      // ------------------------------------------------- update plastic state

      // -------------------------- update hardening and damage variables

      // check if damage variable is acceptable
      // admissible values: (0 <= D < 1) or (1 >= omega > 0)
      // sanity check: omega < 1.0e-20
      if (omega < 1.0e-20)
        dserror("INadmissible value of integrity: omega = %-14.8E in ele = %4d!"
          " \n Omega has to be greater than zero!", eleID, omega);

      // update damage variable damage_{n+1}
      damage = 1.0 - omega;

      // --> damaged isotropic hardening variable has newest value (see L678)

      // final von Mises equivalent stress q_{n+1} = sqrt(3 J)
      // using consistency condition for admissible state
      // Phi_{n+1} == 0 = q / (1-D_{n+1}) - sigma_y(Rplast_{n+1})
      // --> q = sigma_y(Rplast_{n+1}) * omega_{n+1}
      // or alternatively: q = omega_n * q_tilde - 3.0 * G * Dgamma;
      double q = omega * sigma_y;

      // get damaged pressure
      double p = omega * p_tilde;

      // deviatoric stress
      // s_{n+1} = (1 - D_{n+1}) s_{n+1}^{trial}
      //           - 2 . G . Delta gamma . sqrt(3/2) . s^{~,trial}/||s^{~,trial}||
      // or alternative
      // s_{n+1} = s_{n+1}^{~,trial} * q / q_tilde
      // with deviatoric stress^{~} = 2 . G . devstrain
      devstress.Update( (q/q_tilde), devstress_tilde, 0.0);
      // alternatively with identical results:
      // s_{n+1} = [omega - 3G Dgamma / q_tilde] . s_{n+1}^{~,trial}

      // total stress
      // sigma_{n+1} = s_{n+1} + p_{n+1} . id2
      // pressure/volumetric stress no influence due to plasticity
      Stress( p, devstress, *stress );

      // -------------------------------------------- update flow vectors

      // unit flow vector Nbar = s_{n+1} / || s_{n+1} ||
      double devstressnorm = 0.0;
      devstressnorm = sqrt( devstress(0)*devstress(0)
                                  + devstress(1)*devstress(1)
                                  + devstress(2)*devstress(2)
                            + 2 * ( devstress(3)*devstress(3)
                                    + devstress(4)*devstress(4)
                                    + devstress(5)*devstress(5) ) );
      Nbar.Update(1.0, devstress, 0.0);
      Nbar.Scale(1 / devstressnorm);

      // flow vector N = sqrt(3/2) . Nbar . 1/omega (Box 12.3 (iv))
      N.Update( (sqrt(3.0 / 2.0) / omega), Nbar, 0.0 );

      // total strains
      // strain^e_{n+1} = strain^(e,trial)_{n+1} - Dgamma . N
      // or alternatively
      //   strain^e_{n+1} = volstrain^{e,trial} + 1/2G . s_{n+1}
      //     = volstrain^{e,trial} + (1 - 3G . Dgamma / (omega . q_tilde) ) . devstrain
      strain_e.Update( 1.0, trialstrain_e, 0.0);
      strain_e.Update( (-Dgamma), N, 1.0 );

      // strain^p_{n+1} = strain^p_n + Dgamma . N
      // or alternatively
      //   strain^p_{n+1} = strain_{n+1} - strain^e_{n+1}
      strain_p.Update( Dgamma, N, 1.0 );

      // compute converged engineering strain components (Voigt-notation)
      for (int i=3; i<6; ++i) strain_e(i) *= 2.0;
      for (int i=3; i<6; ++i) strain_p(i) *= 2.0;

      // pass the current plastic strains to the element (for visualisation)
      plstrain.Update(1.0, strain_p, 0.0);

      // ------------------------------------------------- update history
      // plastic strain
      strainplcurr_->at(gp) = strain_p;

      // accumulated plastic strain
      strainbarplcurr_->at(gp) = strainbar_p;

      // update damaged isotropic hardening variable
      isohardvarcurr_->at(gp) = Rplast;

      // update damage variable damage_{n+1}
      damagecurr_->at(gp) = damage;

#ifdef DEBUGMATERIAL
      std::cout << "end strain_p\n " << strain_p << std::endl;
      std::cout << "end strainplcurr_->at(gp)\n " << strainplcurr_->at(gp) << std::endl;
#endif //ifdef DEBUGMATERIAL

    }  // damage evolution has to be considered, damage threshold exceeded

  }  // plastic corrector

  //------------------------------------------------------------------------
  // ELSE: elastic step (Phi_trial <= 0.0, Dgamma = 0.0)
  //------------------------------------------------------------------------
  else  // (Phi_trial <= 0.0)
  {
    // -------------------------- update stress using damaged elastic law
    // omega_{n+1} = omega_n
    omega = omegaold;

    // get damaged pressure
    // p = omega_{n+1} . p_tilde
    double p = p_tilde * omega;

    // get damaged deviatoric stresses
    // s_{n+1} = omega_{n+1} . s^{trial}_{n+1}
    devstress.Update(omega, devstress_tilde, 0.0);

    // result vectors of time step n+1 = omega . trial state vectors
    // sigma^e_n+1 = omega . sigma^(e,trial)_n+1
    //             = omega . (s^{trial}_{n+1} + p . id2)
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
    //  - accumulated (un)damaged plastic strains
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
  // (damaged) elastic unloading --> C == 1/(1 - D_{n+1})C_e
  else
    heaviside = 0.0;

  // using an associative flow rule: C_ep is symmetric
  // ( generally C_ep is nonsymmetric )
  SetupCmatElastoPlastic(
    *cmat,
    eleID,
    Dgamma,
    G,
    kappa,
    p_tilde,
    q_tilde,
    energyrelrate,
    Ytan,
    sigma_y,
    Hiso,
    Nbar,
    gp,
    damevolution,
    heaviside
    );

#ifdef DEBUGMATERIAL
  std::cout << "Nach Setup Cep\n" << std::endl;
  std::cout << " Dgamma " << Dgamma << std::endl;
  std::cout << " G " << G << std::endl;
  std::cout << " q " << q << std::endl;
  std::cout << " flow vector " << Nbar << std::endl;
  std::cout << " heaviside " << heaviside << std::endl;
  std::cout << "--> cmat " << cmat << std::endl;
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

  // isotropic elasticity tensor C in Voigt matrix notation
  // e.g. see lecture notes of NiliFEM Script (3.134)
  // C = 2G I + lambda (id2 \otimes id2) with I == I_s
  //   = 2G I_d + kappa (id2 \otimes id2)
  //
  //                       [ 1-nu     nu     nu |          0    0    0 ]
  //                       [        1-nu     nu |          0    0    0 ]
  //           E           [               1-nu |          0    0    0 ]
  //   C = --------------- [ ~~~~   ~~~~   ~~~~   ~~~~~~~~~~  ~~~  ~~~ ]
  //       (1+nu)*(1-2*nu) [                    | (1-2*nu)/2    0    0 ]
  //                       [                    |      (1-2*nu)/2    0 ]
  //                       [ symmetric          |           (1-2*nu)/2 ]
  //
  const double mfac = young / ( (1.0 + nu) * (1.0 - 2.0 * nu) );

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
  int eleID,  // current element ID
  double Dgamma,  // plastic multiplier
  double G,  // shear modulus
  double kappa,  // bulk modulus
  double p_tilde,  // undamaged pressure
  double q_tilde,  // undamaged trial von Mises equivalent stress
  double energyrelrate,  // damage energy release rate
  double Ytan,  // derivative of engergy release rate Ytan w.r.t. Dgamma
  double sigma_y,  // current yield stress
  double Hiso,  // isotropic hardening modulus
  LINALG::Matrix<NUM_STRESS_3D,1> Nbar,  // unit flow vector
  int gp,  // current Gauss point
  bool damevolution,  // flag indicating if damage evolves or not
  double heaviside  // Heaviside function
  )
{
  // damage threshold is still not passed, i.e. use undamaged material tangent
  if (damevolution == false)
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
      if (q_tilde != 0.0)
      {
        epfac = (-1.0) * heaviside * Dgamma * 6 * G * G / q_tilde;
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
          if (q_tilde != 0.0)
          {
            epfac3 =  heaviside * 6 * G * G * ( Dgamma / q_tilde - 1.0 / (3 * G + Hiso) );
            // here: Nbar = s^{trial}_{n+1} / || s^{trial}_{n+1} ||
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
  }  // (damevolution == false)

  // material tangent differs for case damage or not
  // if no damage: use standard tangent of purely plastic behaviour
  else  // (damevolution == true)
  {
#ifdef DEBUGMATERIAL
    if (gp == 0)
      std::cout << "damage evolution takes place in eleID = " << eleID << endl;
#endif // #ifdef DEBUGMATERIAL

    // incremental constitutive function for the stress tensor
    // consistent tangent operator
    // D^{ep} := dsigma_n+1 / dstrain_n+1

    // depending on the flow vector Cmat_ep can be a fully-occupied matrix

    // ---------------------------------------------------------- Heaviside

    // if plastic loading:   heaviside = 1.0 --> use D^{ep}
    // if elastic unloading: heaviside = 0.0 --> use D^e

    // -------------------------------------------------- identity matrices

    // deviatoric projection tensor: I_d = I_s - 1/3 I \otimes I
    // I_d in Voigt-notation applied to symmetric problem, like stress calculation
    //         [ 2/3   -1/3  -1/3 | 0    0    0  ]
    //         [-1/3    2/3  -1/3 | 0    0    0  ]
    //         [-1/3   -1/3   2/3 | 0    0    0  ]
    //   I_d = [ ~~~~  ~~~~  ~~~~  ~~~  ~~~  ~~~ ]
    //         [                  | 1/2   0   0  ]
    //         [    symmetric     |      1/2  0  ]
    //         [                  |          1/2 ]

    // build Cartesian identity 2-tensor I_{AB}
    LINALG::Matrix<NUM_STRESS_3D,1> id2(true);
    for (int i=0; i<3; i++) id2(i) = 1.0;

    // set Cartesian identity 4-tensor in 6-Voigt matrix notation
    // this is fully 'contra-variant' identity tensor, i.e. I^{ABCD}
    // REMARK: rows are stress-like 6-Voigt
    //         columns are stress-like 6-Voigt
    LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D> id4sharp(true);
    for (int i=0; i<3; i++) id4sharp(i,i) = 1.0;
    for (int i=3; i<6; i++) id4sharp(i,i) = 0.5;

    // ------------------------------------- extract current history values
    // integrity omega_{n+1}
    double omega = 1 - damagecurr_->at(gp);

    // ----------------------------------------------------- damaged elastic term
    // D^{ep} = (1 - D_{n+1}) . C^e = omega_{n+1} . C^e
    //        = omega_{n+1} . 2G . I_d + omega_{n+1} . kappa . id2 \otimes id2
    // add standard isotropic elasticity tensor C^e first
    if (heaviside == 0)
    {
      SetupCmat(cmat);
      cmat.Scale(omega);
    }
    else // (heaviside == 1)
    {
      // D^{ep} = a . I_d + b . Nbar \otimes Nbar + c . Nbar \otimes id2
      //          + d . id2 \otimes Nbar + e . id2 \otimes id2

      // ----------------------------------- extract some material parameters

      // get material parameters
      // damage evolution law denominator r
      double damden = params_->damden_;
      // damage evolution law exponent s
      double damexp = params_->damexp_;

      // ---------------------------------------------------------- plastic terms

      // check if omega is admissible
      if (omega < 1.0e-20)
        dserror("Omega has to be greater than zero! omega = %-14.8E\n", omega);

      // be aware: in the 1st equilibrium (Newton) iteration (i=0) D^{ep} is
      // indeterminate due to Dgamma == 0.0
      // a small perturbation Dgamma=1e-08 is used (instead of using a limiting
      // procedure) which was proposed by de Souza Neto
      if (Dgamma == 0)
        Dgamma = 1.0e-08;

      // ------------------------ factors required for the elasto-plastic tangent

      // some constants
      double aux = -energyrelrate / damden;
      double auxb = (q_tilde - sigma_y) / (3 * G);
      // PhiT = qtilde - sigma_y(Rplast_{n+1})
      // be careful: NOT Phi_trial which was calculated using sigma_y(Rplast_n)
      double PhiT = q_tilde - sigma_y;

      // ---------------------------------------- derivatives w.r.t. Dgamma
      // domega/dDgamma
      double Domega = (3.0 * G + omega * Hiso) / PhiT;
      // domega/dq_tilde . dq_tilde/dDgamma
      double DomegaDq_tilde = - omega / PhiT;

      // derviative of residual function dF/dDgamma
      double ResTan = Domega - Hiso / (3 * G) * pow(aux, damexp)
                      - auxb * damexp * Ytan / damden * pow(aux, (damexp -1) );
      // derivative of residual function dF/dp_tilde . dp_tilde/dDgamma
      // DResDp_tilde = s . (q_tilde - sigma_y)/(3G) . pow((-Y/r), s-1) . p_tilde /(r . kappa)
      double DResDp_tilde = damexp * auxb * pow(aux, (damexp -1) ) * p_tilde
                            / (damden * kappa);
      // derivative of residual function F w.r.t. q_tilde
      double DResDq_tilde = DomegaDq_tilde + pow(aux, damexp) / (3 * G);

      // --------- calculate the constants according to de Souza Neto (12.52)

      // a1 = - DResDq_tilde / ResTan
      //    = [ omega / (q_tilde - sigma_y) - pow( (-energyrelrate/damden), damexp)
      //                                      / (3.0 * G) ] / ResTan;
      // or alternatively
      //   a1 = [1 - Dgamma / (omega . omega) . pow( (-Y/r), s)] . omega
      //          / [ (q_tilde - sigma_y) . ResTan]
      double a1 = - DResDq_tilde / ResTan;

      // a2 = - DResDp_tilde / ResTan
      //    = - s . p_tilde . (q_tilde - sigma_y) / (3G . r . kappa . ResTan) . pow((-Y/r), s-1);
      double a2 = 0.0;
      a2 = - DResDp_tilde / ResTan;

      // a3 = a2 * Domega
      double a3 = 0.0;
      a3 = a2 * Domega;

      // a4 = a1 . Domega + DomegaDq_tilde
      //    = a1 . Domega - omega / (q_tilde - sigma_y)
      double a4 = 0.0;
      a4 = a1 * Domega + DomegaDq_tilde;

      // ------- calculate actual coefficients of elasto-plastic material tangent

      // a = 2G . sigma_y . omega / q_tilde
      double a = 0.0;
      a = (2.0 * G * sigma_y * omega) / q_tilde;
      // b = 2G . (a1 . Hiso . omega + a4 . sigma_y - sigma_y . omega / q_tilde
      double b = 0.0;
      b = 2.0 * G * (a1 * Hiso * omega + a4 * sigma_y - sigma_y * omega / q_tilde);
      // c = kappa . (a2 . Hiso . omega + a3 . sigma_y) / sqrt(3/2)
      double c = 0.0;
      c = kappa * (a2 * Hiso * omega + a3 * sigma_y) / sqrt(3.0/2.0);
      // d = p_tilde . 2G . sqrt(3/2) . a4
      double d = 0.0;
      d = ( p_tilde * 2.0 * G * sqrt(3.0/2.0) * a4 );
      // e = kappa . (omega + p_tilde . a3)
      double e = 0.0;
      e = kappa * (omega + p_tilde * a3);

      // ------------------------------- assemble elasto-plastic material tangent

      // empty consistent tangent operator
      cmat.Clear();
      // constitutive tensor
      // I_d = id4sharp - 1/3 Id \otimes Id
      // contribution: Id4^#
      cmat.Update(a, id4sharp, 1.0);
      // contribution: Id \otimes Id
      double epfac = 0.0;
      epfac = a / (-3.0);
      if (q_tilde != 0.0)
      {
        cmat.MultiplyNT(epfac, id2, id2, 1.0);
      }
      cmat.MultiplyNT(b, Nbar, Nbar, 1.0);
      cmat.MultiplyNT(c, Nbar, id2, 1.0);
      cmat.MultiplyNT(d, id2, Nbar, 1.0);
      cmat.MultiplyNT(e, id2, id2, 1.0);

    }  // plastic load step: (heaviside == 1)
    // complete material tangent C_ep available

#ifdef DEBUGMATERIAL
    if (Dgamma != 0)
    {
      std::cout << "Ende SetupCmatElastPlast" << std::endl;
      std::cout << "Cep\n" << " Dgamma " << Dgamma << std::endl;
      std::cout << " G " << G << std::endl;
      std::cout << " q_tilde " << q_tilde << std::endl;
      std::cout << " Nbar " << Nbar << std::endl;
      std::cout << " heaviside " << heaviside << std::endl;
      std::cout << " epfac " << epfac << std::endl;
      std::cout << " epfac1 " << epfac1 << std::endl;
      std::cout << " cmat " << cmat << std::endl;
    }
#endif // #ifdef DEBUGMATERIAL

  }  // damage evolves: (damevolution == true)

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
  for (int i=0; i<samplenumber; ++i)
  {
    // astrain^{p}_{n+1} > astrain^{p}_ref^[i]
    if (strainbar_p >= strainbar_p_ref[i])
    {
      // astrain^{p}_{n+1} > astrain^{p}_ref^[i] --> sigma_y = sigma_ref[i]
      // --> Hiso = d sigma_y / d astrain^{p}_{n+1} = 0
      Hiso = 0.0;
      continue;
    }

    // (strainbar_p < strainbar_p_ref[i])
    else
    {
      // load is still elastic: astrain^{p}_{n+1} < astrain^{p}_ref^{i=0}
      if (i == 0)
      {
        // yield boundary is the initial yield stress (sigma_y^{i=0})
        Hiso = 0.0;
        continue;
      }
      // astrain^{p,i-1}_ref < astrain^{p}_{n+1} < astrain^{p,i}_ref
      else
      {
        //         sigma_y_n - sigma_y^{i-1}
        // Hiso =  ---------------------------------------
        //        astrain^{p,i}_ref - astrain^{p,i-1}_ref
        Hiso = (sigma_y_ref[i] - sigma_y_ref[i-1]) / (strainbar_p_ref[i] - strainbar_p_ref[i-1]);
        continue;
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
  const double strainbar_p  // current accumulated strain, in case of dependent hardening
                             // if damage!=0: isotropic hardening internal variable
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
          continue;
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
  
  std::string isohardeningvar = "isohardeningvar";
  names[isohardeningvar] = 1; // scalar

  std::string damage = "damage";
  names[damage] = 1; // scalar
  
}  // VisNames()


/*---------------------------------------------------------------------*
 | return visualization data (public)                       dano 09/13 |
 *---------------------------------------------------------------------*/
bool MAT::Damage::VisData(
  const std::string& name,
  std::vector<double>& data,
  int numgp
  )
{
  if (name == "accumulatedstrain")
  {
    if ((int)data.size() != 1) dserror("size mismatch");
    double temp = 0.0;
    for (int iter=0; iter<numgp; iter++)
      temp += AccumulatedStrain(iter);
    data[0] = temp/numgp;
  }

  if (name == "isohardeningvar")
  {
    if ((int)data.size() != 1) dserror("size mismatch");
    double temp = 0.0;
    for (int iter=0; iter<numgp; iter++)
      temp += IsotropicHardeningVariable(iter);
    data[0] = temp/numgp;
  }

  if (name == "damage")
  {
    if ((int)data.size() != 1) dserror("size mismatch");
    double temp = 0.0;
    for (int iter=0; iter<numgp; iter++)
      temp += IsotropicDamage(iter);
    data[0] = temp/numgp;
  }

  return true;
}  // VisData()


/*----------------------------------------------------------------------*/

