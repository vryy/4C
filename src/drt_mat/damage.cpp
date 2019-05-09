/*----------------------------------------------------------------------*/
/*!
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
         DAMTHRESHOLD 1.0e-06 KINHARD 17 KINHARD_REC 21 TOL 1.0e-6

\level 2

<pre>
\maintainer Matthias Mayr
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
MAT::PAR::Damage::Damage(Teuchos::RCP<MAT::PAR::Material> matdata)
    : Parameter(matdata),
      youngs_(matdata->GetDouble("YOUNG")),
      poissonratio_(matdata->GetDouble("NUE")),
      density_(matdata->GetDouble("DENS")),
      sigma_y_(*(matdata->Get<std::vector<double>>("SIGMA_Y"))),
      strainbar_p_ref_(*(matdata->Get<std::vector<double>>("EPSBAR_P"))),
      damden_(matdata->GetDouble("DAMDEN")),
      damexp_(matdata->GetDouble("DAMEXP")),
      epsbarD_(matdata->GetDouble("DAMTHRESHOLD")),
      kinhard_(matdata->GetDouble("KINHARD")),
      kinhard_rec_(matdata->GetDouble("KINHARD_REC")),
      sathardening_(matdata->GetDouble("SATHARDENING")),
      hardexpo_(matdata->GetDouble("HARDEXPO")),
      abstol_(matdata->GetDouble("TOL"))
{
  if (hardexpo_ < 0.0) dserror("Nonlinear hardening exponent must be non-negative!");
  if (damden_ == 0.0)
    dserror("Denominator has to be unequal to zero, otherwise floating point exception!");
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
DRT::ParObject* MAT::DamageType::Create(const std::vector<char>& data)
{
  MAT::Damage* plastic = new MAT::Damage();
  plastic->Unpack(data);
  return plastic;
}


/*----------------------------------------------------------------------*
 | constructor (public)                                      dano 04/11 |
 *----------------------------------------------------------------------*/
MAT::Damage::Damage() : params_(NULL), plastic_step_(false) {}


/*----------------------------------------------------------------------*
 | copy-constructor (public)                                 dano 04/11 |
 *----------------------------------------------------------------------*/
MAT::Damage::Damage(MAT::PAR::Damage* params) : params_(params), plastic_step_(false) {}


/*----------------------------------------------------------------------*
 | pack (public)                                             dano 04/11 |
 *----------------------------------------------------------------------*/
void MAT::Damage::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data, type);

  // matid
  int matid = -1;
  // in case we are in post-process mode
  if (params_ != NULL) matid = params_->Id();
  AddtoPack(data, matid);

  // pack history data
  int histsize;
  // if material is not initialised, i.e. start simulation, nothing to pack
  if (!Initialized())
  {
    histsize = 0;
  }
  else
  {
    // if material is initialised (restart): size equates number of gausspoints
    histsize = strainpllast_->size();
  }
  AddtoPack(data, histsize);  // Length of history vector(s)
  for (int var = 0; var < histsize; ++var)
  {
    // insert history vectors to AddtoPack
    AddtoPack(data, strainpllast_->at(var));
    AddtoPack(data, backstresslast_->at(var));

    AddtoPack(data, strainbarpllast_->at(var));
    AddtoPack(data, isohardvarlast_->at(var));
    AddtoPack(data, damagelast_->at(var));
  }

  AddtoPack(data, plastic_step_);

  return;
}  // Pack()


/*----------------------------------------------------------------------*
 | unpack (public)                                           dano 04/11 |
 *----------------------------------------------------------------------*/
void MAT::Damage::Unpack(const std::vector<char>& data)
{
  isinit_ = true;
  std::vector<char>::size_type position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position, data, type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");

  // matid and recover params_
  int matid;
  ExtractfromPack(position, data, matid);
  params_ = NULL;
  if (DRT::Problem::Instance()->Materials() != Teuchos::null)
    if (DRT::Problem::Instance()->Materials()->Num() != 0)
    {
      const int probinst = DRT::Problem::Instance()->Materials()->GetReadFromProblem();
      MAT::PAR::Parameter* mat =
          DRT::Problem::Instance(probinst)->Materials()->ParameterById(matid);
      if (mat->Type() == MaterialType())
        params_ = static_cast<MAT::PAR::Damage*>(mat);
      else
        dserror("Type of parameter material %d does not fit to calling type %d", mat->Type(),
            MaterialType());
    }

  // history data
  int histsize;
  ExtractfromPack(position, data, histsize);

  // if system is not yet initialised, the history vectors have to be intialised
  if (histsize == 0) isinit_ = false;

  // unpack plastic strain vectors
  strainpllast_ = Teuchos::rcp(new std::vector<LINALG::Matrix<NUM_STRESS_3D, 1>>);
  strainplcurr_ = Teuchos::rcp(new std::vector<LINALG::Matrix<NUM_STRESS_3D, 1>>);

  // unpack back stress vectors (for kinematic hardening)
  backstresslast_ = Teuchos::rcp(new std::vector<LINALG::Matrix<NUM_STRESS_3D, 1>>);
  backstresscurr_ = Teuchos::rcp(new std::vector<LINALG::Matrix<NUM_STRESS_3D, 1>>);

  // unpack accumulated plastic strain
  strainbarpllast_ = Teuchos::rcp(new std::vector<double>);
  strainbarplcurr_ = Teuchos::rcp(new std::vector<double>);

  // unpack isotropic hardening variable
  isohardvarlast_ = Teuchos::rcp(new std::vector<double>);
  isohardvarcurr_ = Teuchos::rcp(new std::vector<double>);

  // unpack isotropic damage internal state variable
  damagelast_ = Teuchos::rcp(new std::vector<double>);
  damagecurr_ = Teuchos::rcp(new std::vector<double>);

  // initialise
  LINALG::Matrix<NUM_STRESS_3D, 1> tmp_vect(true);
  double tmp_scalar = 0.0;

  for (int var = 0; var < histsize; ++var)
  {
    // vectors of last converged state are unpacked
    ExtractfromPack(position, data, tmp_vect);
    strainpllast_->push_back(tmp_vect);

    ExtractfromPack(position, data, tmp_vect);
    backstresslast_->push_back(tmp_vect);

    // scalar-valued vector of last converged state are unpacked
    ExtractfromPack(position, data, tmp_scalar);
    strainbarpllast_->push_back(tmp_scalar);

    ExtractfromPack(position, data, tmp_scalar);
    isohardvarlast_->push_back(tmp_scalar);

    ExtractfromPack(position, data, tmp_scalar);
    damagelast_->push_back(tmp_scalar);

    // current vectors have to be initialised
    strainplcurr_->push_back(tmp_vect);
    backstresscurr_->push_back(tmp_vect);

    strainbarplcurr_->push_back(tmp_scalar);
    isohardvarcurr_->push_back(tmp_scalar);
    damagecurr_->push_back(tmp_scalar);
  }

  plastic_step_ = false;
  int plastic_step;
  ExtractfromPack(position, data, plastic_step);

  // if it was already plastic before, set true
  if (plastic_step != 0) plastic_step_ = true;

  if (position != data.size()) dserror("Mismatch in size of data %d <-> %d", data.size(), position);

  return;

}  // Unpack()


/*---------------------------------------------------------------------*
 | initialise / allocate internal stress variables (public)      04/11 |
 *---------------------------------------------------------------------*/
void MAT::Damage::Setup(int numgp, DRT::INPUT::LineDefinition* linedef)
{
  // initialise history variables
  strainpllast_ = Teuchos::rcp(new std::vector<LINALG::Matrix<NUM_STRESS_3D, 1>>);
  strainplcurr_ = Teuchos::rcp(new std::vector<LINALG::Matrix<NUM_STRESS_3D, 1>>);

  backstresslast_ = Teuchos::rcp(new std::vector<LINALG::Matrix<NUM_STRESS_3D, 1>>);
  backstresscurr_ = Teuchos::rcp(new std::vector<LINALG::Matrix<NUM_STRESS_3D, 1>>);

  strainbarpllast_ = Teuchos::rcp(new std::vector<double>);
  strainbarplcurr_ = Teuchos::rcp(new std::vector<double>);

  isohardvarlast_ = Teuchos::rcp(new std::vector<double>);
  isohardvarcurr_ = Teuchos::rcp(new std::vector<double>);

  damagelast_ = Teuchos::rcp(new std::vector<double>);
  damagecurr_ = Teuchos::rcp(new std::vector<double>);

  // set all history variables to zero
  LINALG::Matrix<NUM_STRESS_3D, 1> emptyvect(true);
  strainpllast_->resize(numgp);
  strainplcurr_->resize(numgp);

  backstresslast_->resize(numgp);
  backstresscurr_->resize(numgp);

  strainbarpllast_->resize(numgp);
  strainbarplcurr_->resize(numgp);

  isohardvarlast_->resize(numgp);
  isohardvarcurr_->resize(numgp);

  damagelast_->resize(numgp);
  damagecurr_->resize(numgp);

  for (int i = 0; i < numgp; i++)
  {
    strainpllast_->at(i) = emptyvect;
    strainplcurr_->at(i) = emptyvect;

    backstresslast_->at(i) = emptyvect;
    backstresscurr_->at(i) = emptyvect;

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
  backstresslast_ = backstresscurr_;

  strainbarpllast_ = strainbarplcurr_;
  isohardvarlast_ = isohardvarcurr_;
  damagelast_ = damagecurr_;

  // empty vectors of current data
  strainplcurr_ = Teuchos::rcp(new std::vector<LINALG::Matrix<NUM_STRESS_3D, 1>>);
  backstresscurr_ = Teuchos::rcp(new std::vector<LINALG::Matrix<NUM_STRESS_3D, 1>>);

  strainbarplcurr_ = Teuchos::rcp(new std::vector<double>);
  isohardvarcurr_ = Teuchos::rcp(new std::vector<double>);
  damagecurr_ = Teuchos::rcp(new std::vector<double>);

  // get the size of the vector
  // (use the last vector, because it includes latest results, current is empty)
  const int histsize = strainpllast_->size();
  strainplcurr_->resize(histsize);
  backstresscurr_->resize(histsize);

  strainbarplcurr_->resize(histsize);
  isohardvarcurr_->resize(histsize);
  damagecurr_->resize(histsize);

  const LINALG::Matrix<NUM_STRESS_3D, 1> emptyvec(true);
  for (int i = 0; i < histsize; i++)
  {
    strainplcurr_->at(i) = emptyvec;
    backstresscurr_->at(i) = emptyvec;

    strainbarplcurr_->at(i) = 0.0;
    isohardvarcurr_->at(i) = 0.0;
    damagecurr_->at(i) = 0.0;
  }

  return;
}  // Update()


/*----------------------------------------------------------------------*
 | evaluate material (public)                                dano 08/11 |
 *----------------------------------------------------------------------*/
void MAT::Damage::Evaluate(const LINALG::Matrix<3, 3>* defgrd,
    const LINALG::Matrix<NUM_STRESS_3D, 1>* linstrain,  // linear strain vector
    Teuchos::ParameterList& params,            // parameter list for communication & HISTORY
    LINALG::Matrix<NUM_STRESS_3D, 1>* stress,  // 2nd PK-stress
    LINALG::Matrix<NUM_STRESS_3D, NUM_STRESS_3D>* cmat,  // material stiffness matrix
    const int eleGID)
{
  // in case kinematic hardening is ignored, use implementation according to de
  // Souza Neto, Computational Methods for Plasticity
  if ((params_->kinhard_ == 0.0) and (params_->kinhard_rec_ == 0.0) and (params_->hardexpo_ == 0.0))
    EvaluateSimplifiedLemaitre(defgrd, linstrain, params, stress, cmat, eleGID);
  // in case full Lemaitre material model is considered, i.e. including
  // kinematic hardening, use implementation according to Doghri
  else
    EvaluateFullLemaitre(defgrd, linstrain, params, stress, cmat, eleGID);
}  // Evaluate


/*----------------------------------------------------------------------*
 | evaluate material (public)                                dano 08/11 |
 *----------------------------------------------------------------------*/
void MAT::Damage::EvaluateSimplifiedLemaitre(const LINALG::Matrix<3, 3>* defgrd,
    const LINALG::Matrix<NUM_STRESS_3D, 1>* linstrain,  // linear strain vector
    Teuchos::ParameterList& params,            // parameter list for communication & HISTORY
    LINALG::Matrix<NUM_STRESS_3D, 1>* stress,  // 2nd PK-stress
    LINALG::Matrix<NUM_STRESS_3D, NUM_STRESS_3D>* cmat,  // material stiffness matrix
    const int eleGID)
{
  const int gp = params.get<int>("gp", -1);
  if (gp == -1) dserror("no Gauss point number provided in material");
  if (eleGID == -1) dserror("no element provided in material");
  LINALG::Matrix<MAT::NUM_STRESS_3D, 1> plstrain(true);

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
  G = young / (2.0 * (1.0 + nu));
  // bulk modulus bulk = E /( 3 ( 1 - 2 nu) )= lambda + 2/3 * mu
  double bulk = 0.0;
  bulk = young / (3.0 * (1.0 - 2.0 * nu));

  // build Cartesian identity 2-tensor I_{AB}
  LINALG::Matrix<NUM_STRESS_3D, 1> id2(true);
  for (int i = 0; i < 3; i++) id2(i) = 1.0;

  // linstrain (in): independent variable passed from the element
  //  strain^p: evolution is determined by the flow rule, history varible
  //  strain^e: definition of additive decomposition:
  //  strain^e = strain - strain^p
  // REMARK: stress-like 6-Voigt vector
  LINALG::Matrix<NUM_STRESS_3D, 1> strain(*linstrain);

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
    std::cout << "Rplast am ele = " << eleGID << ": " << Rplast << std::endl;
    dserror("damaged isotropic hardening variable has to be equal to or greater than zero!");
  }

  // get old integrity: omega_n = 1 - D_n
  double omegaold = 0.0;
  omegaold = 1.0 - damage;

  // ------------------------------------------------ old plastic strains

  // plastic strain vector
  // strain^{p,trial}_{n+1} = strain^p_n
  LINALG::Matrix<NUM_STRESS_3D, 1> strain_p(true);
  for (int i = 0; i < 6; i++) strain_p(i, 0) = strainpllast_->at(gp)(i, 0);

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
  for (int i = 3; i < 6; ++i) strain(i) /= 2.0;
  for (int i = 3; i < 6; ++i) strain_p(i) /= 2.0;

  // ----------------------------------------------- elastic trial strain
  // assume load step is elastic
  // strain^{e}_{n+1}
  LINALG::Matrix<NUM_STRESS_3D, 1> strain_e(true);

  // strain^{e,trial}_{n+1} = strain_n+1 - strain^p_n
  LINALG::Matrix<NUM_STRESS_3D, 1> trialstrain_e(false);
  trialstrain_e.Update(1.0, strain, 0.0);
  trialstrain_e.Update((-1.0), strain_p, 1.0);

  // volumetric strain
  // trace of strain vector
  double tracestrain = trialstrain_e(0) + trialstrain_e(1) + trialstrain_e(2);
  // volstrain = 1/3 . tr( strain ) . Id
  LINALG::Matrix<NUM_STRESS_3D, 1> volumetricstrain(false);
  volumetricstrain.Update((tracestrain / 3.0), id2, 0.0);

  // deviatoric strain
  // devstrain^e = strain^e - volstrain^e
  LINALG::Matrix<NUM_STRESS_3D, 1> devstrain(false);
  devstrain.Update(1.0, trialstrain_e, 0.0);
  devstrain.Update(-1.0, volumetricstrain, 1.0);

  // --------------------------------------------- trial undamaged stress

  // undamaged scalar-valued pressure
  // p^{~} = bulk . tr( strain )
  double p_tilde = bulk * tracestrain;

  // deviatoric stress^{~} = 2 . G . devstrain
  LINALG::Matrix<NUM_STRESS_3D, 1> devstress_tilde(false);
  devstress_tilde.Update(2.0 * G, devstrain);
  // be careful for shear stresses (sigma_12)
  // in Voigt-notation the shear strains have to be scaled with 1/2
  // normally done in the material tangent (cf. id4sharp)

  // --------------- trial (undamaged) elastic von Mises effective stress

  // q^{~,trial}_{n+1} := q(s^{trial}_{n+1}) / (1-D_n) = \sqrt{ 3 . J2 } / (1-D_n)
  //                    = sqrt{3/2} . || s^{trial}_{n+1} || / (1-D_n)
  // J2 = 1/2 (s11^2 + s22^2 + s33^2 + 2 . s12^2 + 2 . s23^2 + 2 . s13^2)
  double J2 = 0.0;
  J2 = 1.0 / 2.0 *
           (devstress_tilde(0) * devstress_tilde(0) + devstress_tilde(1) * devstress_tilde(1) +
               devstress_tilde(2) * devstress_tilde(2)) +
       devstress_tilde(3) * devstress_tilde(3) + devstress_tilde(4) * devstress_tilde(4) +
       devstress_tilde(5) * devstress_tilde(5);
  double q_tilde = 0.0;
  q_tilde = sqrt(3.0 * J2);

  // initialise final (damaged) deviatoric stresses
  LINALG::Matrix<NUM_STRESS_3D, 1> devstress(true);

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
  // with trial values: Phi_trial = q_trial^{~} - sigma_y and Dgamma == 0
  double Phi_trial = 0.0;
  Phi_trial = q_tilde - sigma_y;

  // --------------------------------------------------------- initialise

  // if trial state is violated, i.e. it's a plastic load step, there are two
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
  LINALG::Matrix<NUM_STRESS_3D, 1> Nbar(true);

  // flow vector N (Prandtl-Reuss)
  // (using the updated deviatoric stress s_n+1, no longer s_n+1^trial)
  // N = sqrt{3/2} / (1 - D_{n+1}) . ( s_{n+1} / || s_{n+1} || )
  LINALG::Matrix<NUM_STRESS_3D, 1> N(true);

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
      int itnum = 0;           // iteration counter

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

        // if not converged (m > m_max)
        if (itnum > itermax)
        {
          dserror("local Newton iteration did not converge after iteration %3d/%3d with Res=%3d",
              itnum, itermax, Res);
        }
        // else: continue loop m <= m_max

        // Res:= Phi = q^(trial)_{n+1} - sigma_y
        // Res = q - 3 * G * Dgamma - sigma_y;
        // with sigma_y = sigma_y(strainbar_p + Dgamma)
        Res = q_tilde - 3.0 * G * Dgamma - sigma_y;

        // check for convergence
        double norm = abs(Res);
        // check: absolute value of Res has to be smaller than given tolerance
        if (norm < (params_->abstol_))
        {
#ifdef DEBUGMATERIAL
          if (gp == 0)
            printf("Newton method converged after %i iterations; abs(Res)=  %-14.8E\n", itnum,
                abs(Res));
#endif  // #ifdef DEBUGMATERIAL
          break;
        }

        // plasticity with piecewise linear isotropic hardening
        // ResTan = -3G -Hiso = const.
        ResTan = -3.0 * G - Hiso;

        // incremental plastic multiplier Dgamma
        // Dgamma^{m} = Dgamma^{m-1} - Phi / Phi'
        Dgamma += (-Res) / ResTan;

        // -------------------------- local Newton update of plastic values

        // compute new residual of accumulatd plastic strains
        // astrain^p_{n+1} = astrain^p_n + Dgamma
        // astrain^p_{n+1} = SUM{Dgamma_n} from all time steps n
        // Kuhn-Tucker: Dgamma >= 0.0 --> astrain^p_{n+1} >= 0.0
        strainbar_p = strainbarpllast_->at(gp) + Dgamma / omega;
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
        if (gp == 0) std::cout << "No damage! Current load step is admissible!" << std::endl;
#endif  // #ifdef DEBUGMATERIAL

        // no damage occurs, i.e. current solution is correct, no return-map
        // considering damage is necessary
        damevolution = false;

        // --------------------------------------------------- plastic update

        // ---------------------------------------------- update flow vectors

        // deviatoric stress norm || s^{trial}_{n+1} ||
        double devstress_tildenorm = 0.0;
        devstress_tildenorm =
            sqrt(devstress_tilde(0) * devstress_tilde(0) + devstress_tilde(1) * devstress_tilde(1) +
                 devstress_tilde(2) * devstress_tilde(2) +
                 2.0 * (devstress_tilde(3) * devstress_tilde(3) +
                           devstress_tilde(4) * devstress_tilde(4) +
                           devstress_tilde(5) * devstress_tilde(5)));

        // unit flow vector Nbar = s^{trial}_{n+1} / || s^{trial}_{n+1} ||
        Nbar.Update((1.0 / devstress_tildenorm), devstress_tilde);

        // flow vector N = sqrt(3/2) . Nbar
        N.Update((sqrt(3.0 / 2.0)), Nbar);

        // deviatoric stress
        // s = s_{n+1}^{trial} - 2 . G . Delta gamma . N
        const double facdevstress = (-2.0) * G * Dgamma;
        devstress.Update(omega, devstress_tilde, facdevstress, N);

        // total stress
        // sigma_{n+1} = s_{n+1} + p_{n+1} . id2
        // pressure/volumetric stress no influence due to plasticity
        Stress(p_tilde, devstress, *stress);

        // total strains
        // strain^e_{n+1} = strain^(e,trial)_{n+1} - Dgamma . N
        // compute converged engineering strain components (Voigt-notation)
        strain_e.Update(1.0, trialstrain_e, (-Dgamma), N);

        // strain^p_{n+1} = strain^p_n + Dgamma . N
        strain_p.Update(Dgamma, N, 1.0);

        // compute converged engineering strain components (Voigt-notation)
        for (int i = 3; i < 6; ++i) strain_e(i) *= 2.0;
        for (int i = 3; i < 6; ++i) strain_p(i) *= 2.0;

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
#endif  // ifdef DEBUGMATERIAL

      }  // (strainbar_p < strainbar_p_D)
      // updated strainbar_p not valid, recalculate step considering damage
      else  // (strainbar_p > strainbar_p_D)
      {
#ifdef DEBUGMATERIAL
        if (gp == 0)
          std::cout << "New solution strainbar_p^m exceeds damage threshold!"
                       "\n Recalculate load step considering damage!"
                    << std::endl;
#endif  // #ifdef DEBUGMATERIAL

        // damage occurs, i.e. current solution is not correct, recalculate
        // load step considering damage
        damevolution = true;
      }  // (strainbar_p > strainbar_p_D)
    }    // no damage: (strainbar_p < strainbar_p_D)

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
#ifdef DEBUGMATERIAL
      // only first plastic call is output at screen for every processor
      // visualisation of whole plastic behaviour via PLASTIC_STRAIN in postprocessing
      if ((plastic_step_ == false) and (gp == 0))
      {
        std::cout << "damage starts to evolve in element = " << eleID << std::endl;

        plastic_step_ = true;
      }

      std::cout << "Damage has to be considered for current load step and ele = " << eleID
                << ", and gp = " << gp << " ! Threshold exceeded!" << std::endl;
#endif  // #ifdef DEBUGMATERIAL

      // ------------------------------------------------- return-mapping

      // local Newton-Raphson

      // ------------------------------- initial guess for Dgamma (12.49)

      // instead of initial guess Dgamma^{m=0} = 0.0 use perfectly plastic
      // solution with frozen yield surface at beginning of each load increment
      //
      // Dgamma^{m=0} = q^{~,trial} - sigma(Rplast_n) * omega_n / 3G
      Dgamma = omegaold * Phi_trial / (3.0 * G);

      // -------------------- initialise internal variables with Dgamma^0

      // Rplast = R^{p,m=0}_{n+1} = R^p_n + Dgamma
      Rplast = isohardvarlast_->at(gp) + Dgamma;

      const int itermax = 50;  // max. number of iterations
      int itnum = 0;           // iteration counter

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
              "local Newton iteration did not converge after iteration %3d/%3d with Res=%3d in "
              "ele=%3d",
              itnum, itermax, Res, eleGID);
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
        Ytan = -Hiso * sigma_y / (3.0 * G);

        // integrity
        // omega_{n+1} = 3G / (q_tilde - sigma_y) * Dgamma = 1 - D_{n+1}
        omega = 3.0 * G / (q_tilde - sigma_y) * Dgamma;

        // damage energy release rate only implicitely depending on Dgamma (12.47)
        energyrelrate = -(sigma_y * sigma_y) / (6.0 * G) - p_tilde * p_tilde / (2.0 * bulk);

        // compute residual function (12.48)
        // Res := F(Dgamma) = omega(Dgamma) - omega_n
        //                    + Dgamma / omega(Dgamma) . (-Y(Dgamma)/r)^s
        //                      . (q_tilde - sigma_y) / (3 G)
        // here: it is important NOT to use Dgamma^{m=0}=0 --> omega=0 --> '1/0'
        double Res =
            omega - omegaold +
            std::pow((-energyrelrate / damden), damexp) / ((3.0 * G) / (q_tilde - sigma_y));

        // check for convergence
        double norm = abs(Res);
        // check: absolute value of Res has to be smaller than given tolerance
        if (norm < (params_->abstol_))
        {
#ifdef DEBUGMATERIAL
          if (gp == 0)
            printf("Newton method converged after %i iterations; abs(Res)=  %-14.8E\n", itnum,
                abs(Res));
#endif  // #ifdef DEBUGMATERIAL
          break;
        }

        // if load state is not converged, calculate derivatives w.r.t. Dgamma

        // derviative of residual w.r.t. Dgamma
        // ResTan = (3 . G)/(q_tilde - sigma_y)
        //          + (3 . G)/(q_tilde - sigma_y) . Dgamma . Hiso / (q_tilde - sigma_y)
        //          - Hiso / (3 . G) . (-Y/r)^s
        //          - s . Ytan / ( ( (3 . G)/(q_tilde - sigma_y) ) . r) . (-Y/r)^(s-1)
        ResTan = (3.0 * G) / (q_tilde - sigma_y) +
                 (3.0 * G) / (q_tilde - sigma_y) * Dgamma * Hiso / (q_tilde - sigma_y) -
                 Hiso / (3.0 * G) * std::pow((-energyrelrate / damden), damexp) -
                 damexp * Ytan / (((3.0 * G) / (q_tilde - sigma_y)) * damden) *
                     std::pow((-energyrelrate / damden), (damexp - 1.0));

        // incremental plastic multiplier Dgamma
        // Dgamma^{m} = Dgamma^{m-1} - Phi / Phi'
        Dgamma += (-Res) / ResTan;

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
        strainbar_p = strainbarpllast_->at(gp) + Dgamma / omega;
        if (strainbar_p < 0.0)
        {
          std::cout << "in element:" << eleGID
                    << ": strainbarpllast_->at(gp) = " << strainbarpllast_->at(gp)
                    << ", omega = " << omega << ", Dgamma = " << Dgamma
                    << ", and strainbar_p = " << strainbar_p << std::endl;
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
        dserror(
            "INadmissible value of integrity: omega = %-14.8E in ele = %4d!"
            " \n Omega has to be greater than zero!",
            omega, eleGID);

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
      devstress.Update((q / q_tilde), devstress_tilde);
      // alternatively with identical results:
      // s_{n+1} = [omega - 3G Dgamma / q_tilde] . s_{n+1}^{~,trial}

      // total stress
      // sigma_{n+1} = s_{n+1} + p_{n+1} . id2
      // pressure/volumetric stress no influence due to plasticity
      Stress(p, devstress, *stress);

      // -------------------------------------------- update flow vectors

      // unit flow vector Nbar = s_{n+1} / || s_{n+1} ||
      double devstressnorm = 0.0;
      devstressnorm = sqrt(devstress(0) * devstress(0) + devstress(1) * devstress(1) +
                           devstress(2) * devstress(2) +
                           2.0 * (devstress(3) * devstress(3) + devstress(4) * devstress(4) +
                                     devstress(5) * devstress(5)));
      Nbar.Update((1.0 / devstressnorm), devstress);

      // flow vector N = sqrt(3/2) . Nbar . 1/omega (Box 12.3 (iv))
      N.Update((sqrt(3.0 / 2.0) / omega), Nbar);

      // total strains
      // strain^e_{n+1} = strain^(e,trial)_{n+1} - Dgamma . N
      // or alternatively
      //   strain^e_{n+1} = volstrain^{e,trial} + 1/2G . s_{n+1}
      //     = volstrain^{e,trial} + (1 - 3G . Dgamma / (omega . q_tilde) ) . devstrain
      strain_e.Update(1.0, trialstrain_e, (-Dgamma), N);

      // strain^p_{n+1} = strain^p_n + Dgamma . N
      // or alternatively
      //   strain^p_{n+1} = strain_{n+1} - strain^e_{n+1}
      strain_p.Update(Dgamma, N, 1.0);

      // compute converged engineering strain components (Voigt-notation)
      for (int i = 3; i < 6; ++i) strain_e(i) *= 2.0;
      for (int i = 3; i < 6; ++i) strain_p(i) *= 2.0;

      // pass the current plastic strains to the element (for visualisation)
      plstrain.Update(strain_p);

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
#endif  // ifdef DEBUGMATERIAL

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
    devstress.Update(omega, devstress_tilde);

    // result vectors of time step n+1 = omega . trial state vectors
    // sigma^e_n+1 = omega . sigma^(e,trial)_n+1
    //             = omega . (s^{trial}_{n+1} + p . id2)
    Stress(p, devstress, *stress);

    // total strains
    // strain^e_{n+1} = strain^(e,trial)_{n+1}
    // compute converged engineering strain components (Voigt-notation)
    strain_e.Update(trialstrain_e);
    for (int i = 3; i < 6; ++i) strain_e(i) *= 2.0;

    // no plastic yielding
    Dgamma = 0.0;

    // pass the current plastic strains to the element (for visualisation)
    // compute converged engineering strain components (Voigt-notation)
    for (int i = 3; i < 6; ++i) strain_p(i) *= 2.0;
    plstrain.Update(strain_p);

    // --------------------------------------------------------- update history
    // constant values for
    //  - plastic strains
    //  - accumulated (un)damaged plastic strains
    //  - stress

    // as current history vectors are set to zero in Update(), the old values
    // need to be set instead, otherwise no constant plastic values are possible
    strainplcurr_->at(gp) = strainpllast_->at(gp);
    strainbarplcurr_->at(gp) = strainbarpllast_->at(gp);
    isohardvarcurr_->at(gp) = isohardvarlast_->at(gp);
    damagecurr_->at(gp) = damagelast_->at(gp);

  }  // elastic step

  //---------------------------------------------------------------------------
  // --------------------------------- consistent elastoplastic tangent modulus
  //---------------------------------------------------------------------------

  // using an associative flow rule: C_ep is symmetric
  // ( generally C_ep is nonsymmetric )
  // if Phi^trial = 0: two tangent stress-strain relations exist
  // plastic loading --> C == C_ep
  if (Dgamma > 0.0) heaviside = 1.0;
  // (damaged) elastic unloading --> C == 1/(1 - D_{n+1})C_e
  else
    heaviside = 0.0;

  // using an associative flow rule: C_ep is symmetric
  // ( generally C_ep is nonsymmetric )
  SetupCmatElastoPlastic(*cmat, eleGID, Dgamma, G, bulk, p_tilde, q_tilde, energyrelrate, Ytan,
      sigma_y, Hiso, Nbar, gp, damevolution, heaviside);

#ifdef DEBUGMATERIAL
  std::cout << "Nach Setup Cep\n" << std::endl;
  std::cout << " Dgamma " << Dgamma << std::endl;
  std::cout << " G " << G << std::endl;
  std::cout << " q " << q << std::endl;
  std::cout << " flow vector " << Nbar << std::endl;
  std::cout << " heaviside " << heaviside << std::endl;
  std::cout << "--> cmat " << cmat << std::endl;
#endif  // #ifdef DEBUGMATERIAL

  // ------------------------------- return plastic strains for post-processing
  params.set<LINALG::Matrix<MAT::NUM_STRESS_3D, 1>>("plglstrain", plstrain);

  return;

}  // EvaluateSimplifiedLemaitre()


/*----------------------------------------------------------------------*
 | evaluate full Lemaitre material model (public)            dano 11/13 |
 *----------------------------------------------------------------------*/
void MAT::Damage::EvaluateFullLemaitre(const LINALG::Matrix<3, 3>* defgrd,
    const LINALG::Matrix<NUM_STRESS_3D, 1>* linstrain,  // linear strain vector
    Teuchos::ParameterList& params,            // parameter list for communication & HISTORY
    LINALG::Matrix<NUM_STRESS_3D, 1>* stress,  // 2nd PK-stress
    LINALG::Matrix<NUM_STRESS_3D, NUM_STRESS_3D>* cmat,  // material stiffness matrix
    int EleGID)
{
  // --------- full Lemaitre material model requires solution of five equations

  // 1) fulfil consistency condition
  // Phi = \sqrt(3 . J2(s_tilde_{n+1} - beta_{n+1})) - sigma_y(R_{n+1}) <= 0
  //
  // 2) determine the stress
  // sigma_{n+1} = (1-D) . C^e : (strain_n+1^{e,trial} - Dgamma . N_{n+1} )
  //
  // 3) determine back stress
  // beta_{n+1} = beta_n + (Hkin . N_tilde_{n+1} - Hkin_rec . beta_n) . Dgamma/(1 + Hkin_rec Dgamma)
  //
  // 4) determine scalar-valued damage
  // D_{n+1} = D_n + 1/(1 - D_{n+1}) . y(s_tilde) . Dgamma / (1 - D)

  const int gp = params.get<int>("gp", -1);
  if (gp == -1) dserror("no Gauss point number provided in material");

  if (EleGID == -1) dserror("no element provided in material");
  LINALG::Matrix<MAT::NUM_STRESS_3D, 1> plstrain(true);

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
  // kinematic hardening modulus, also defined as Hkin = k_1 = a
  double Hkin = params_->kinhard_;
  // scalar-value variable describing recovery, also defined as Hkin_rec = k_2 = b
  double Hkin_rec = params_->kinhard_rec_;
  // hardening exponent
  double hardexpo = params_->hardexpo_;
  // saturation hardening
  double sigma_yinfty = params_->sathardening_;

  // calculate some more paramters
  // lame constant / shear modulus parameter mu == G
  double G = 0.0;
  G = young / (2.0 * (1.0 + nu));
  // bulk modulus bulk = E /( 3 ( 1 - 2 nu) ) = lambda + 2/3 * mu
  double bulk = 0.0;
  bulk = young / (3.0 * (1.0 - 2.0 * nu));

  // build Cartesian identity 2-tensor I_{AB}
  LINALG::Matrix<NUM_STRESS_3D, 1> id2(true);
  for (int i = 0; i < 3; i++) id2(i) = 1.0;

  // linstrain (in): independent variable passed from the element
  //  strain^p: evolution is determined by the flow rule, history varible
  //  strain^e: definition of additive decomposition:
  //  strain^e = strain - strain^p
  // REMARK: stress-like 6-Voigt vector
  LINALG::Matrix<NUM_STRESS_3D, 1> strain(*linstrain);

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
    std::cout << "Rplast am ele = " << EleGID << ": " << Rplast << std::endl;
    dserror("damaged isotropic hardening variable has to be equal to or greater than zero!");
  }

  // get old integrity: omega_n = 1 - D_n
  double omegaold = 0.0;
  omegaold = 1.0 - damage;

  // integrity
  // omega_{n+1} = 1 - D_{n+1}
  double omega = 1.0;  // if not damaged, omega == 1.0

  // ------------------------------------------------ old plastic strains

  // plastic strain vector
  // strain^{p,trial}_{n+1} = strain^p_n
  LINALG::Matrix<NUM_STRESS_3D, 1> strain_p(false);
  strain_p.Update(strainpllast_->at(gp));

  // get old accumulated/equivalent plastic strain only in case of plastic step
  double strainbar_p = 0.0;
  // accumulated or equivalent plastic strain (scalar-valued)
  // astrain^{p,trial}_{n+1} = astrain^p_n
  strainbar_p = strainbarpllast_->at(gp);
  if (strainbarpllast_->at(gp) < 0.0)
    dserror("accumulated plastic strain has to be equal to or greater than zero!");

  // ------------------------------------------------ old back stress
  // beta^{trial}_{n+1} = beta_n
  // beta is a deviatoric tensor
  LINALG::Matrix<NUM_STRESS_3D, 1> beta(false);
  beta.Update(backstresslast_->at(gp));

  // --------------------------------------------------- physical strains
  // convert engineering shear components into physical components
  // input strain is given in Voigt-notation

  // convert engineering shear component (in) into physical component
  for (int i = 3; i < 6; ++i) strain(i) /= 2.0;
  for (int i = 3; i < 6; ++i) strain_p(i) /= 2.0;

  // ----------------------------------------------- elastic trial strain
  // assume load step is elastic
  // strain^{e}_{n+1}
  LINALG::Matrix<NUM_STRESS_3D, 1> strain_e(true);

  // strain^{e,trial}_{n+1} = strain_n+1 - strain^p_n
  LINALG::Matrix<NUM_STRESS_3D, 1> trialstrain_e(false);
  trialstrain_e.Update(1.0, strain, (-1.0), strain_p);

  // volumetric strain
  // trace of strain vector
  double tracestrain = trialstrain_e(0) + trialstrain_e(1) + trialstrain_e(2);
  // volstrain = 1/3 . tr( strain ) . Id
  LINALG::Matrix<NUM_STRESS_3D, 1> volumetricstrain(false);
  volumetricstrain.Update((tracestrain / 3.0), id2, 0.0);

  // deviatoric strain
  // devstrain^e = strain^e - volstrain^e
  LINALG::Matrix<NUM_STRESS_3D, 1> devstrain(false);
  devstrain.Update(1.0, trialstrain_e, (-1.0), volumetricstrain);

  // --------------------------------------------- trial undamaged stress

  // undamaged scalar-valued pressure
  // p^{~} = bulk . tr( strain )
  double p_tilde = bulk * tracestrain;

  // undamaged deviatoric stress^{~} = 2 . G . devstrain
  LINALG::Matrix<NUM_STRESS_3D, 1> devstress_tilde(false);
  devstress_tilde.Update((2.0 * G), devstrain);
  // be careful for shear stresses (sigma_12)
  // in Voigt-notation the shear strains have to be scaled with 1/2
  // normally done in the material tangent (cf. id4sharp)

  // -----------------------------------------  relative effective stress
  // eta^{~,trial}_{n+1} = s^{~,trial}_{n+1} - beta^{trial}_{n+1}
  LINALG::Matrix<NUM_STRESS_3D, 1> eta_tilde(true);
  RelStress(devstress_tilde, beta, eta_tilde);

  // --------------- trial (undamaged) elastic von Mises effective stress

  // q^{~,trial}_{n+1} := q(s^{trial}_{n+1}) / (1 - D_n) = \sqrt{ 3 . J2 } / (1-D_n)
  //                    = sqrt{3/2} . || s^{trial}_{n+1} || / (1 - D_n)
  // J2 = 1/2 (s11^2 + s22^2 + s33^2 + 2 . s12^2 + 2 . s23^2 + 2 . s13^2)
  double J2 = 0.0;
  J2 = 1.0 / 2.0 *
           (devstress_tilde(0) * devstress_tilde(0) + devstress_tilde(1) * devstress_tilde(1) +
               devstress_tilde(2) * devstress_tilde(2)) +
       +devstress_tilde(3) * devstress_tilde(3) + devstress_tilde(4) * devstress_tilde(4) +
       devstress_tilde(5) * devstress_tilde(5);
  double q_tilde = 0.0;
  q_tilde = sqrt(3.0 * J2);

  // ------ trial (undamaged) relative elastic von Mises effective stress

  // qbar^{~,trial}_{n+1} := sqrt{3/2} . || s^{~,trial}_{n+1} - beta^{trial}_{n+1} ||
  //                       = sqrt{3/2} . || eta^{~,trial}_{n+1} ||
  //                       = qbar(eta^{trial}_{n+1}) / (1 - D_n)
  //                       = \sqrt{ 3 . J2bar(eta) } / (1 - D_n)

  // J2bar = 1/2 (eta11^2 + eta22^2 + eta33^2 + 2 . eta12^2 + 2 . eta23^2 + 2 . eta13^2)
  double J2bar = 0.0;
  J2bar = 1.0 / 2.0 *
              (eta_tilde(0) * eta_tilde(0) + eta_tilde(1) * eta_tilde(1) +
                  eta_tilde(2) * eta_tilde(2)) +
          eta_tilde(3) * eta_tilde(3) + eta_tilde(4) * eta_tilde(4) + eta_tilde(5) * eta_tilde(5);
  double qbar_tilde = 0.0;
  qbar_tilde = sqrt(3.0 * J2bar);

#ifdef DEBUGMATERIAL
  if (gp == 0)
  {
    std::cout << ": devstress_tilde\n " << devstress_tilde << std::endl;
    std::cout << ": devstrain\n " << devstrain << std::endl;
    std::cout << ": beta\n " << beta << std::endl;
    std::cout << ": eta_tilde\n " << eta_tilde << std::endl;
    std::cout << "plastic load: strainbarplcurr_->at(gp)\n " << strainbarplcurr_->at(gp)
              << std::endl;
    std::cout << "plastic load: strainbarpllast_->at(gp)\n " << strainbarpllast_->at(gp)
              << std::endl;
    std::cout << "plastic load: strain_p\n " << strain_p << std::endl;
    std::cout << "plastic load: strainplcurr_->at(gp)\n " << strainplcurr_->at(gp) << std::endl;
    std::cout << "plastic load: strainpllast\n " << strainpllast_->at(gp) << std::endl;
    std::cout << "elastic load: backstresscurr_->at(gp)\n " << backstresscurr_->at(gp) << std::endl;
    std::cout << "elastic load: backstresslast_->at(gp)\n " << backstresslast_->at(gp) << std::endl;
    std::cout << ": q_tilde\n " << q_tilde << std::endl;
    std::cout << ": qbar_tilde\n " << qbar_tilde << std::endl;
  }
#endif  // DEBUGMATERIAL

  // initialise final (damaged) deviatoric stresses
  LINALG::Matrix<NUM_STRESS_3D, 1> devstress(true);

  // ----------- initialise variables required due to kinematic hardening

  // initialise variables calculated within local Newton and required also for
  // the damaged elastoplastic material tangent C^{ep}
  // initialise yield stress and its derivative
  double sigma_y = 0.0;
  // isotropic thermodynamical force
  double kappa = 0.0;
  // dkappa_dR = sigma_infty . exp (-delta . R) . delta
  double dkappa_dR = 0.0;
  // factor
  double g = 0.0;
  // denominator
  double h_alg = 0.0;
  // s_N = 2 G . (Dgamma / omega)^2 . dy/ds_tilde : N_tilde
  double s_N = 0.0;
  // Nbetaold = N_tilde : beta_n
  double Nbetaold = 0.0;
  // initialise energy-release rate y := (-Y/r)^s
  double y = 0.0;
  // dy/dsigma_tilde
  LINALG::Matrix<NUM_STRESS_3D, 1> dy_dsigma_tilde(true);
  // b_NbetaoldN = beta_n - 2/3 . (N_tilde : beta_n) . N_tilde
  LINALG::Matrix<NUM_STRESS_3D, 1> b_NbetaoldN(true);
  // correction for difference (back stress - deviatoric stress)
  // c_s - c_beta (49)
  LINALG::Matrix<NUM_STRESS_3D, 1> c_s_b(true);

  //---------------------------------------------------------------------------
  // check plastic admissibility, Phi<=0 is admissble
  //---------------------------------------------------------------------------

  // ----------------------------------------------- trial yield function

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
    omegaold = omega = 1.0;

    // calculate the uniaxial yield stress out of samples
    double sigma_y0 = GetSigmaYAtStrainbarnp(strainbar_p);

    // kappa = sigma_yinfty . (1 - exp (-delta . astrain))
    kappa = sigma_yinfty * (1.0 - exp(-hardexpo * strainbar_p));
    // yield stress sigma_y = sigma_y0 + kappa(strainbar_p)
    sigma_y = sigma_y0 + kappa;

    // compute the derivative of the isotropic hardening kappa(R^p) w.r.t. R^p

    // dkappa_dR = sigma_infty . (- exp (-delta . R) ) . (-delta)
    //           = sigma_infty . exp (-delta . R) . delta
    dkappa_dR = sigma_yinfty * exp(-hardexpo * strainbar_p) * (-hardexpo);

  }     // sigma_y(astrain_n)
  else  // current strainbar_p > strainbar_p_D
  {
    // calculate the uniaxial yield stress out of samples
    double sigma_y0 = GetSigmaYAtStrainbarnp(Rplast);

    // kappa = sigma_yinfty . (1 - exp (-delta . astrain))
    kappa = sigma_yinfty * (1.0 - exp(-hardexpo * Rplast));

    // yield stress sigma_y = sigma_y0 + kappa(strainbar_p)
    sigma_y = sigma_y0 + kappa;

    // compute the derivative of the isotropic hardening kappa(R^p) w.r.t. R^p
    // dkappa_dR = sigma_infty . (- exp (-delta . R) ) . (-delta)
    //           = sigma_infty . exp (-delta . R) . delta
    dkappa_dR = sigma_yinfty * exp(-hardexpo * Rplast) * (-hardexpo);
  }  // sigma_y(R_n)

  // calculate the yield function
  // Phi = \sqrt{ 3.0 . J2 } - sigma_y = qbar_tilde - sigma_y
  // with trial values: Phi_trial = qbar{~,trial} - sigma_y and Dgamma == 0
  double Phi_trial = 0.0;
  Phi_trial = qbar_tilde - sigma_y;

#ifdef DEBUGMATERIAL
  if (gp == 0)
  {
    std::cout << ": Phi_trial  " << Phi_trial << std::endl;
    std::cout << ": qbar_tilde  " << qbar_tilde << std::endl;
    std::cout << ": sigma_y  " << sigma_y << std::endl;
    std::cout << ": kappa " << kappa << std::endl;
    std::cout << ": dkappa_dR " << dkappa_dR << std::endl;
  }
#endif  // DEBUGMATERIAL

  // --------------------------------------------------------- initialise

  // if trial state is violated, i.e. it's a plastic load step, there are 2
  // possible states: plastic loading: heaviside = 1, elastic unloading = 0)
  double heaviside = 0.0;
  // incremental plastic multiplier Delta gamma
  double Dgamma = 0.0;

  // unit flow vector Nbar (Prandtl-Reuss)
  // (using s_n+1^trial for undamaged, and s_n+1 for damaged load step)
  // Nbar = ( eta^{~,trial}_{n+1} / || eta^{~,trial}_{n+1} || )
  LINALG::Matrix<NUM_STRESS_3D, 1> Nbar(true);

  // flow vector N_tilde according to Doghri and his explicit updating
  // N_tilde = 3/2 . (s_tilde - beta) / qbar_tilde
  LINALG::Matrix<NUM_STRESS_3D, 1> N_tilde(true);

  //---------------------------------------------------------------------------
  // IF consistency condition is violated, i.e. plastic load step use return-mapping
  // (Phi_trial > 0.0, Dgamma >= 0.0)
  //---------------------------------------------------------------------------
  if (Phi_trial > 1.0e-08)  // if (Phi_trial > 0.0)
  {
    if (plastic_step_ == false)
    {
      if ((plastic_step_ == false) and (gp == 0))
        std::cout << "plasticity starts in element = " << EleGID << std::endl;

      plastic_step_ = true;
    }

    // ------------------------------------------------------- damage threshold

    // calculate damage threshold (according to de Souza Neto, p.483ff)

    // bool which decide if threshold is passed or not
    // current strainbar_p < strainbar_p_D
    if (strainbar_p < strainbar_p_D)
    {
      // damage threshold not exceeded
      omega = omegaold = 1.0;
      damage = 0.0;
    }  // no damage: (strainbar_p < strainbar_p_D)
    // else (strainbar_p > strainbar_p_D), i.e. damage evolves
    // use history values

    // only first plastic call is output at screen for every processor
    // visualisation of whole plastic behaviour via PLASTIC_STRAIN in postprocessing
    if (plastic_step_ == false)
    {
      if ((plastic_step_ == false) and (gp == 0))
        std::cout << "damage starts to evolve in element = " << EleGID << std::endl;

      plastic_step_ = true;
    }

#ifdef DEBUGMATERIAL
    std::cout << "Damage has to be considered for current load step and ele = " << eleID
              << ", and gp = " << gp << " ! Threshold exceeded!" << std::endl;
#endif  // #ifdef DEBUGMATERIAL

    // eta_tilde = s_tilde - beta

    // deviatoric stress norm || eta^{~}_{n+1} ||
    double eta_tildenorm = 0.0;
    eta_tildenorm = sqrt(eta_tilde(0) * eta_tilde(0) + eta_tilde(1) * eta_tilde(1) +
                         eta_tilde(2) * eta_tilde(2) +
                         2.0 * (eta_tilde(3) * eta_tilde(3) + eta_tilde(4) * eta_tilde(4) +
                                   eta_tilde(5) * eta_tilde(5)));

    // undamaged unit flow vector
    // Nbar = eta^{trial}_{n+1} / || eta^{trial}_{n+1} || = eta_{n+1} / || eta_{n+1} ||
    // using the fact that eta_{n+1} is proportional to eta^{trial}_{n+1}
    Nbar.Update((1.0 / eta_tildenorm), eta_tilde);

    // undamaged flow vector N_tilde according to Doghri
    // N_tilde^{trial} = 3/2 . (s_tilde^{trial} - beta_n) / sqrt{ 3/2 .
    //                   (s_tilde^{trial} - beta_n) : (s_tilde^{trial} - beta_n) }
    //                 = 3/2 . (s_tilde^{trial} - beta_n)
    //                   / ( sqrt{ 3/2 } . || (s_tilde^{trial} - beta_n) || )
    //                 = sqrt{ 3/2 } . (s_tilde^{trial} - beta_n)
    //                   /  || (s_tilde^{trial} - beta_n) || = sqrt{ 3/2 } . Nbar
    N_tilde.Update((sqrt(3.0 / 2.0)), Nbar);

    // ------------------------------------------------ calculate corrections

    // 1.step) compute the corrections over the elastic trial predictor
    // 2.step) if residuals > tolerance criteria compute
    //         --> calculate the corrections in a loop over plastic iterations
    //             m until criteria is fulfilled, i.e. r_i < tolerance criteria

    // -------------------- calculate corrections for elastic trial predictor

    // calculate N : beta_n
    double N_tildebeta = 0.0;
    for (int i = 0; i < 6; ++i)
      N_tildebeta += N_tilde(i, 0) * backstresslast_->at(gp)(i, 0);  // (6x1)(6x1)

    // calculate trial value of hardening modulus (27)
    // h_trial = 3 . G + (1 - D_n) . (dkappa/dR + 3/2 . a - b . N : beta_n)
    // dkappa/dR = Hiso with R = R_n, Dgamma = 0
    double h_trial = 0.0;
    h_trial = 3.0 * G + omegaold * (dkappa_dR + 3.0 / 2.0 * Hkin - Hkin_rec * N_tildebeta);

    // correction for the internal hardening variable R (57)
    // R' = c_R/Dt:  c_R = (1-D_{n}) . Phi^{trial}_{n+1} / h^{trial}_{n+1}
    double c_R = 0.0;
    c_R = omegaold * Phi_trial / h_trial;

    // correction for accumulated plastic strain (58a)
    // astrain' = c_astrain/Dt: c_astrain = c_R / (1 - D_n) = c_R / omega_n
    double c_strainbar = 0.0;
    c_strainbar = c_R / omegaold;

    // corrections for effective, undamaged deviatoric stress (58b)
    // s = c_s/Dt: c_s = -2G . N_tilde^{trial} . c_astrain
    LINALG::Matrix<NUM_STRESS_3D, 1> c_s(false);
    double fac_cs = -2.0 * G * c_strainbar;
    c_s.Update(fac_cs, N_tilde);

    // corrections for back stress (58c)
    // beta' = c_beta/Dt: c_beta = (Hkin . N_tilde^{trial} - Hkin_rec . beta_n) . c_R
    LINALG::Matrix<NUM_STRESS_3D, 1> c_beta(false);
    c_beta.Update(Hkin, N_tilde);
    c_beta.Update((-Hkin_rec), beta, 1.0);
    c_beta.Scale(c_R);

    // corrections for damage variable
    // damage energy release rate
    // Y = -1/2(1-D)^2 sigma : C^e : sigma
    //   = - q^2 / [6 G . (1 - D)^2] - p^2 / [2 bulk . (1 - D)^2]
    // with q^{~}_{n+1} = q_{n+1} / (1 - D_{n+1})
    // and p^{~}_{n+1} = p^{~, trial} = p_{n+1}  / (1 - D_{n+1})
    double Y = 0.0;
    Y = -q_tilde * q_tilde / (6.0 * G) - p_tilde * p_tilde / (2.0 * bulk);
    // y := (Y^trial_{n+1} / r)^s
    y = std::pow((-Y / damden), damexp);

    // c_D = (Y^trial_{n+1} / r)^s . c_astrain
    double c_D = 0.0;
    c_D = y * c_strainbar;

#ifdef DEBUGMATERIAL
    if (gp == 0)
    {
      std::cout << ": Phi_trial  " << Phi_trial << std::endl;
      std::cout << ": qbar_tilde  " << qbar_tilde << std::endl;
      std::cout << ": Y  " << Y << std::endl;
      std::cout << ": y " << y << std::endl;
      std::cout << ": c_D " << c_D << std::endl;
    }
#endif  // DEBUGMATERIAL

    // ------------------------------------------------- return-mapping

    // local Newton-Raphson

    const int itermax = 50;  // max. number of iterations
    int itnum = 0;           // iteration counter

    // start iteration with index m for local Newton
    while (true)
    {
      itnum++;
      // check for convergence

      // if not converged m > m_max
      if (itnum > itermax)
      {
        dserror("local Newton iteration did not converge after iteration %3d/%3d", itnum, itermax);
      }
      // else: continue loop m <= m_max

      // ------------------------- update variables within local Newton
      // c_(.)/dt = ()'
      // --> c_(.) == (.)^{m+1} - (.)^m
      //     --> (.)^{m+1} = (.)^m + c_(.)

      // isotropic hardening variable
      // R^{p,m+1} = R^{p,m}_{n+1} = R^p_n + c_R
      Rplast = isohardvarlast_->at(gp) + c_R;

      // (undamaged) accumulated plastic strain
      // astrain^{m+1}_{n+1} = astrain^m_n + c_astrain
      strainbar_p = strainbarpllast_->at(gp) + c_strainbar;

      // (undamaged, effective) deviatoric stress
      // s^{~,m+1}_{n+1} = s^{~,trial/m}_{n+1} + c_s
      devstress_tilde.Update(1.0, c_s, 1.0);

      // back stress
      // beta^{m+1}_{n+1} = beta^{trial/m}_{n+1} + c_beta
      beta.Update(1.0, c_beta, 1.0);

      // check damage threshold
      if (strainbar_p < strainbar_p_D)
      {
        damage = 0.0;
        omega = 1.0;
      }
      else
      {
        // damage
        // D^{m+1}_{n+1} = D^m_n + c_D
        damage = damagelast_->at(gp) + c_D;
        // update omega
        omega = 1.0 - damage;
      }

      // update plastic multiplier
      // Dgamma = c_R (gamma' = R' = c_R/Dt)
      Dgamma = c_R;

      // update q^{~}_{n+1}
      // q^{~}_{n+1} := q(s^{~}_{n+1}) = \sqrt{ 3 . J2(s^{~} }
      //              = sqrt{3/2} . || s^{~}_{n+1} ||
      // J2 = 1/2 (s11^2 + s22^2 + s33^2 + 2 . s12^2 + 2 . s23^2 + 2 . s13^2)
      J2 = 1.0 / 2.0 *
               (devstress_tilde(0) * devstress_tilde(0) + devstress_tilde(1) * devstress_tilde(1) +
                   devstress_tilde(2) * devstress_tilde(2)) +
           +devstress_tilde(3) * devstress_tilde(3) + devstress_tilde(4) * devstress_tilde(4) +
           devstress_tilde(5) * devstress_tilde(5);
      q_tilde = sqrt(3.0 * J2);

      // update energy release rate
      // here: Y := Y(s_tilde), i.e. depends ONLY on stresses, beta NOT included!
      //
      // Y = - [ \sqrt(3/2) . J2(sigma_tilde) ]^2 / (6 G) - p_tilde^2 / (2 bulk)
      // here we assume that Y depends only on deviatoric part
      // --> q_tilde(sigma_tilde) == q_tilde(s_tilde)
      //
      // Y = - [ q_tilde ]^2 / (6 G) - p_tilde^2 / (2 bulk)
      Y = -q_tilde * q_tilde / (6.0 * G) - p_tilde * p_tilde / (2.0 * bulk);
      // y == (-Y(s^{~,m}_{n+1}) / r)^s
      y = std::pow((-Y / damden), damexp);

      // ----------------- calculate derivatives of energy-release rate

      // d(-Y)/ds_tilde
      //   = d{[ q_tilde(s_tilde) ]^2 / (6 G) + p_tilde^2 / (2 bulk)}/ds_tilde
      //   = d{[ q_tilde(s_tilde) ]^2 / (6 G) + 0 }/ds_tilde
      //   = d{[ q_tilde(s_tilde) ]^2 }/ds_tilde / (6 G)
      //   = 2 . q_tilde(s_tilde) . sqrt(3/2) . s_tilde/|| s_tilde || / (6 G)
      //   = s_tilde / (2 G)

      // dy_ds_tilde = s . (-Y / r)^{s-1} . (1 / r) . s_tilde / (2 G)
      //             = s / (2 G . r) . (-Y / r)^{s-1} .s_tilde
      LINALG::Matrix<NUM_STRESS_3D, 1> dy_dstilde(false);
      double fac_dy_dstilde = damexp / (damden * 2.0 * G) * std::pow((-Y / damden), (damexp - 1.0));
      dy_dstilde.Update(fac_dy_dstilde, devstress_tilde);

      // d(-Y)/dsigma_tilde
      //   = d{[ q_tilde(s_tilde) ]^2 }/ds_tilde / (6 G) +
      //     + d{ p_tilde^2 / (2 bulk)}/d (tr sigma_tilde)
      //   = s_tilde / (2 G) + p_tilde / bulk . I
      LINALG::Matrix<NUM_STRESS_3D, 1> dYneg_dsigma_tilde(false);
      dYneg_dsigma_tilde.Update((1.0 / (2.0 * G)), devstress_tilde);
      for (int i = 0; i < 3; ++i) dYneg_dsigma_tilde(i) += bulk * p_tilde;

      // dy_dsigma_tilde = s . (-Y / r)^{s-1} . (1 / r)
      //                   . [ s_tilde / (2 G) + p_tilde / bulk . I ]
      //
      // dy_dsigma_tilde = s . (-Y / r)^{s-1} . (1 / r) . d(-Y)/dsigma_tilde
      double fac_dy_dsigma_tilde = damexp / damden * std::pow((-Y / damden), (damexp - 1.0));
      dy_dsigma_tilde.Update(fac_dy_dsigma_tilde, dYneg_dsigma_tilde);

      // update effective stress
      // eta^{~}_{n+1} = s^{~}_{n+1} - beta_{n+1}
      LINALG::Matrix<NUM_STRESS_3D, 1> eta_tilde(true);
      RelStress(devstress_tilde, beta, eta_tilde);
      // update the invariant of the stress deviator
      J2bar = 1.0 / 2.0 *
                  (eta_tilde(0) * eta_tilde(0) + eta_tilde(1) * eta_tilde(1) +
                      eta_tilde(1) * eta_tilde(1) + eta_tilde(2) * eta_tilde(2)) +
              eta_tilde(3) * eta_tilde(3) + eta_tilde(4) * eta_tilde(4) +
              eta_tilde(5) * eta_tilde(5);
      // qbar_tilde = sqrt( 3/2 eta : eta ) = sqrt( 3/2 ) . || eta ||
      qbar_tilde = sqrt(3.0 * J2bar);

      // calculate the uniaxial yield stress
      // in contrast to implementation of de Souza Neto, no pairs are used
      // sigma_y = sigma_y0 + kappa(R^p)
      // kappa(R^p) = sigma_infty . (1 - exp(-delta . R^p) )
      // initial yield stress sigma_y0 for R^p == 0;
      double sigma_y0 = GetSigmaYAtStrainbarnp(0);
      kappa = sigma_yinfty * (1.0 - exp(-hardexpo * Rplast));
      sigma_y = sigma_y0 + kappa;

      // update flow vectors: deviatoric stress norm || eta^{~}_{n+1} ||
      double eta_tildenorm = 0.0;
      eta_tildenorm = sqrt(eta_tilde(0) * eta_tilde(0) + eta_tilde(1) * eta_tilde(1) +
                           eta_tilde(2) * eta_tilde(2) +
                           2.0 * (eta_tilde(3) * eta_tilde(3) + eta_tilde(4) * eta_tilde(4) +
                                     eta_tilde(5) * eta_tilde(5)));

      // update undamaged unit flow vector
      // Nbar = eta^{~}_{n+1} / || eta^{~}_{n+1} || = eta_{n+1} / || eta_{n+1} ||
      // using the fact that due eta_{n+1} is proportional to eta^{trial}_{n+1}
      Nbar.Update((1.0 / eta_tildenorm), eta_tilde);

      // update undamaged flow vector N_tilde
      // N_tilde = 3/2 . (s^{~} - beta_n) / sqrt{ 3/2 . || s^{~} - beta_n || }
      N_tilde.Update((sqrt(3.0 / 2.0)), Nbar);

      // -------------------------------------------------- compute residuals

      // ---------------- compute residual of deviatoric stresses (42a)

      // k_s_tilde = s^{~}_{n+1} - s^{~,trial} + 2G . N_tilde . Dgamma/omega
      // with s_tilde^{trial} = 2G . devstrain
      // k_s_tilde --> = 0
      LINALG::Matrix<NUM_STRESS_3D, 1> k_s_tilde(false);
      k_s_tilde.Update((-2.0 * G), devstrain);
      k_s_tilde.Update(1.0, devstress_tilde, 1.0);
      double fac_ks_tilde = 0.0;
      fac_ks_tilde = 2.0 * G * Dgamma / omega;
      k_s_tilde.Update(fac_ks_tilde, N_tilde, 1.0);

      // -------------- compute residual of consistency condition (42b)

      // Phi = qbar_tilde - sigma_y(R^p)
      double Phi = qbar_tilde - sigma_y;

      // ------------------------ compute residual of back stress (42c)

      // k_b = beta_{n+1} - beta_n - (Hkin . N_tilde - Hkin_rec . beta_n)
      //                    . Dgamma / (1 + Hkin_rec . Dgamma)
      // k_b = beta_{n+1} - beta_n
      //       - (Hkin . Dgamma / (1 + Hkin_rec . Dgamma) . N_tilde
      //       + (Hkin_rec . Dgamma / (1 + Hkin_rec . Dgamma) . beta_n
      double fac_k_b1 = 0.0;
      fac_k_b1 = (-1.0) * Hkin * Dgamma / (1.0 + Hkin_rec * Dgamma);
      double fac_k_b2 = 0.0;
      fac_k_b2 = Hkin_rec * Dgamma / (1.0 + Hkin_rec * Dgamma);

      LINALG::Matrix<NUM_STRESS_3D, 1> k_b(false);
      k_b.Update((-1.0), backstresslast_->at(gp));
      k_b.Update(1.0, beta, 1.0);
      k_b.Update(fac_k_b1, N_tilde, 1.0);
      k_b.Update(fac_k_b2, backstresslast_->at(gp), 1.0);

      // ----------------------------- compute residual of damage (42d)

      // k_D = D_{n+1} - D_n - y . Dgamma / (1 - D_{n+1})
      double k_D = 0.0;
      k_D = damage - damagelast_->at(gp) - y * Dgamma / omega;

      // ------------------------------------- calculate residual norms

      // calculate norm of deviatoric residual k_s_tilde
      double norm_k_s_tilde = 0.0;
      norm_k_s_tilde = k_s_tilde.Norm2();

      // calculate norm of consistency condition Phi
      double norm_Phi = 0.0;
      norm_Phi = abs(Phi_trial);

      // calculate norm of back stress evolution equation
      double norm_k_beta = 0.0;
      norm_k_beta = k_b.Norm2();

      // calculate norm of damage evolution equation
      double norm_k_D = 0.0;
      norm_k_D = abs(k_D);

      // ---------------------------------------------- check for convergence

      // load step is convergence ONLY if all equations are convergenced!
      if ((norm_k_s_tilde < params_->abstol_) and (norm_Phi < params_->abstol_) and
          (norm_k_beta < params_->abstol_) and (norm_k_D < params_->abstol_))
      {
        if (omega < 0.0)
        {
          dserror("Damage variable has converged to unacceptable value");
        }
        break;
      }
      else
#ifdef DEBUGMATERIAL
          if (gp == 0)
        printf(
            "Newton method converged after %i iterations; "
            "norm_k_s = %-14.8E, "
            "norm_k_Phi = %-14.8E, "
            "norm_k_b = %-14.8E, "
            "norm_k_D = %-14.8E \n",
            itnum, norm_k_s_tilde, norm_Phi, norm_k_beta, norm_k_D);
#endif  // #ifdef DEBUGMATERIAL

      // --------------- else: load step NOT converged, calculate corrections

      // the sequence of corrections is decisive due to the dependency of single
      // terms on the former ones

      // ------------------------------------- 1.) calculating c_R (54)
      // c_R = num / h_alg

      // calculate numerator num

      // num = [ (1 - D_{n+1}) - y . Delta_astrain
      //         + 2 G . (Delta_astrain)^2 . dy/ds_tilde : N_tilde ]
      //       . [Phi - N_tilde : (k_s_tilde - k_beta) ] +
      //       + 3 G . [ Delta_astrain . k_D + (Delta_astrain)^2 .
      //                 . dy/ds_tilde : k_s_tilde ] -
      //       - 3 G . (Delta_astrain)^2 . 3 G / qbar_tilde . Delta_astrain
      //         / (1 + 3/2 . g) . dy/ds_tilde :
      // [ k_s_tilde - k_b - 2/3 . ( N_tilde : (k_s_tilde - k_b) ) . N_tilde ]
      //
      // num = [ omega - y . (Dgamma / omega) + 2 G . (Dgamma / omega)^2
      //         . dy/ds_tilde : N_tilde ] . [ Phi - N_tilde : (k_s_tilde - k_b) ]
      //       + 3 G . [ (Dgamma / omega) . k_D + (Dgamma / omega)^2 .
      //                 dy/ds_tilde : k_s_tilde ]
      //       - 3 G . (Dgamma / omega)^2 . 3 G / qbar_tilde . (Dgamma / omega)
      //        / (1 + 3/2 . g) . dy/ds_tilde :
      // [ k_s_tilde - k_b - 2/3 . ( N_tilde : (k_s_tilde - k_b) ) . N_tilde ]

      // calculate the factor g included in third term of num (47)
      // g = [ 2 G . Delta_astrain + Hkin . Dgamma / (1 + Hkin_rec . Dgamma) ] / qbar_tilde
      //   = [ 2 G .( Dgamma / omega) + Hkin . Dgamma / (1 + Hkin_rec . Dgamma) ] / qbar_tilde
      if (qbar_tilde != 0.0)
      {
        g = (2.0 * G * Dgamma / omega + Hkin * Dgamma / (1.0 + Hkin_rec * Dgamma)) / qbar_tilde;
      }
      else
        dserror("do not divide by zero!");

      if (g <= 0.0) dserror("factor g has to be greater zero! g =  %-14.8E", g);

      // h_alg = 3 G + [ (1 - D_{n+1}) - y . Delta_astrain +
      //         + 2 G . (Delta_astrain)^2 . dy/ds_tilde : N_tilde ]
      //           . [ dkappa/dR + (3/2 . Hkin - Hkin_rec . (N_tilde : beta_n) )
      //                             / (1 + Hkin_rec . Dgamma)^2 ]
      //        - 3 G . (Delta_astrain)^2 . 3 G / qbar_tilde . Delta_astrain
      //          / (1 + 3/2 . g) . Hkin_rec / (1 + Hkin_rec . Dgamma)^2
      //          . dy/ds_tilde : [ beta_n - 2/3 . (N_tilde : beta_n) . N^tilde ]
      //
      // h_alg = 3 G + [ omega - y . (Dgamma / omega) +
      //         + 2 G . (Dgamma / omega)^2 . dy/ds_tilde : N_tilde ]
      //           . [ dkappa/dR + (3/2 . Hkin - Hkin_rec . (N_tilde : beta_n) )
      //                             / (1 + Hkin_rec . Dgamma)^2 ]
      //        - 3 G . (Dgamma / omega)^2 . 3 G / qbar_tilde . (Dgamma / omega)
      //          / (1 + 3/2 . g) . Hkin_rec / (1 + Hkin_rec . Dgamma)^2
      //          . dy/ds_tilde : [ beta_n - 2/3 . (N_tilde : beta_n) . N^tilde ]

      // dy/ds_tilde : N_tilde
      double dydstilde_Ntilde = 0.0;
      for (int i = 0; i < 6; ++i)
        dydstilde_Ntilde += dy_dstilde(i, 0) * N_tilde(i, 0);  // (6x1)(6x1)
      // s_N = 2 G . (Dgamma / omega)^2 . dy/ds_tilde : N_tilde
      s_N = 2.0 * G * dydstilde_Ntilde;

      // Nbetaold = N_tilde : beta_n
      for (int i = 0; i < 6; ++i)
        Nbetaold += N_tilde(i, 0) * backstresslast_->at(gp)(i, 0);  // (6x1)(6x1)

      // b_NbetaoldN = beta_n - 2/3 . (N_tilde : beta_n) . N_tilde
      b_NbetaoldN.Update(1.0, backstresslast_->at(gp), (-2.0 / 3.0 * Nbetaold), N_tilde);

      // dydsb_NbetaoldN = dy_dstilde : [ beta_n - 2/3 . (N_tilde : beta_n) . N_tilde ]
      //                  = dy/ds_tilde . b_NbetaoldN
      double dydsb_NbetaoldN = 0.0;
      for (int i = 0; i < 6; ++i) dydsb_NbetaoldN += dy_dstilde(i, 0) * b_NbetaoldN(i, 0);

      // N_ksb = N_tilde : (k_s - k_b)
      double N_ksb = 0.0;
      for (int i = 0; i < 6; ++i) N_ksb += N_tilde(i) * (k_s_tilde(i) - k_b(i));

      // s_k = dy/ds_tilde : k_s
      double s_k = 0.0;
      for (int i = 0; i < 6; ++i) s_k += dy_dstilde(i, 0) * k_s_tilde(i, 0);

      // k_NkN = [ k_s - k_b - 2/3 . ( N_tilde : (k_s - k_b) ) . N_tilde ]
      LINALG::Matrix<NUM_STRESS_3D, 1> k_NkN(false);
      k_NkN.Update(1.0, k_s_tilde, (-1.0), k_b);
      double fac_k_NkN = 0.0;
      fac_k_NkN = -2.0 / 3.0 * N_ksb;
      k_NkN.Update(fac_k_NkN, N_tilde, 1.0);

      // dy_dstildekN = dy/ds_tilde . [ k_s - k_b - 2/3 . ( N_tilde : (k_s - k_b) ) . N_tilde ]
      //              = dy_dstilde . k_NkN
      double dy_dstildek_NkN = 0.0;
      for (int i = 0; i < 6; ++i) dy_dstildek_NkN += dy_dstilde(i, 0) * k_NkN(i, 0);

      // ------------------------------ calculate numerator num for c_R

      // num = [ omega - y . (Dgamma / omega) + 2 G . (Dgamma / omega)^2
      //         . dy/ds_tilde : N_tilde ] .
      //         . [ Phi - N_tilde : (k_s_tilde - k_b) ]
      //       + 3 G . [ (Dgamma / omega) . k_D + (Dgamma / omega)^2 .
      //                 dy/ds_tilde : k_s_tilde ]
      //       - 3 G . (Dgamma / omega)^2 . 3 G / qbar_tilde . (Dgamma / omega)
      //         / (1 + 3/2 . g) . dy/ds_tilde :
      //         : [k_s_tilde - k_b - 2/3 . ( N_tilde : (k_s_tilde - k_b) ) . N_{n+1} ]
      // num = [ omega - y . (Dgamma / omega) + 2 G . (Dgamma / omega)^2
      //         . dy/ds_tilde : N_tilde ] .
      //         . [ Phi - N_tilde : (k_s_tilde - k_b) ]
      //       + 3 G . [ (Dgamma / omega) . k_D + (Dgamma / omega)^2 . s_k ]
      //       - 3 G . (Dgamma / omega)^2 . 3 G / qbar_tilde . (Dgamma / omega)
      //         / (1 + 3/2 . g) . dy_dstildek_NkN
      double num = 0.0;
      num = (omega - y * (Dgamma / omega) + s_N) * (Phi - N_ksb) +
            3.0 * G * ((Dgamma / omega) * k_D + (Dgamma / omega) * (Dgamma / omega) * s_k) -
            3.0 * G * (Dgamma / omega) * (Dgamma / omega) * 3.0 * G / qbar_tilde *
                (Dgamma / omega) / (1.0 + 3.0 / 2.0 * g) * dy_dstildek_NkN;

      // -------------------------- calculate denominator h_alg for c_R

      // h_alg = 3 G + [ omega - y . Dgamma / omega + s_N ] . [ dkappa/dR +
      //                 . (3/2 . Hkin - Hkin_rec . Nbetaold )
      //                 / (1 + Hkin_rec . Dgamma)^2 ]
      //        - 3 G (Dgamma / omega)^2 . 3 G / qbar_tilde . Dgamma / omega
      //          / (1 + 3/2 . g) . Hkin_rec / (1 + Hkin_rec . Dgamma)^2
      //          . dydsb_NbetaoldN
      h_alg = 3.0 * G +
              (omega - y * Dgamma / omega + s_N) *
                  (dkappa_dR + (3.0 / 2.0 * Hkin - Hkin_rec * Nbetaold) /
                                   ((1.0 + Hkin_rec * Dgamma) * (1.0 + Hkin_rec * Dgamma))) -
              3.0 * G * ((Dgamma / omega) * (Dgamma / omega)) * 3.0 * G / qbar_tilde *
                  (Dgamma / omega) / (1.0 + 3.0 / 2.0 * g) * Hkin_rec /
                  ((1.0 + Hkin_rec * Dgamma) * (1.0 + Hkin_rec * Dgamma)) * dydsb_NbetaoldN;

      // ------------------------------------- 1.) calculating c_R (54)

      // c_R = num / h_alg
      c_R = num / h_alg;

      // ------------------------------- 2.) calculate c_strainbar (53)

      // c_strainbar = 1/(3 G) . { Phi - N_tilde : (k_s -k_b)
      //           - [ dkappa_dR + (3/2 . Hkin - Hkin_rec . N_tilde : beta_n)
      //               / (1 + Hkin_rec . Dgamma)^2 ] . c_R }
      c_strainbar = (Phi - N_ksb -
                        (dkappa_dR + (3.0 / 2.0 * Hkin - Hkin_rec * Nbetaold) /
                                         ((1 + Hkin_rec * Dgamma) * (1 + Hkin_rec * Dgamma))) *
                            c_R) /
                    (3.0 * G);

      // --------------------------------- 3.) calculate c_s_tilde (51)

      // dN_tilde/ds_tilde : (c_s - c_b) := dNds_csb (50)

      LINALG::Matrix<NUM_STRESS_3D, 1> bracket(false);
      // bracket = [ k_s - k_b - 2/3 . N_tilde : (k_s - k_b) . N_tilde
      //            - Hkin_rec / (1 + Hkin_rec . Dgamma)^2 .
      //              ( beta_n - 2/3 . (N_tilde : beta_n) . N_tilde ) . c_R ]
      // with k_NkN = [ k_s - k_b - 2/3 . ( N_tilde : (k_s - k_b) ) . N_tilde ]
      //
      // bracket = [ k_NkN + fac_bracket . ( beta_n - 2/3 . Nbetaold . N_tilde ) ]
      double fac_bracket = -Hkin_rec * c_R / ((1 + Hkin_rec * Dgamma) * (1 + Hkin_rec * Dgamma));
      double fac_bracket_1 = fac_bracket * (-2 / 3.0) * Nbetaold;
      bracket.Update(1.0, k_NkN, fac_bracket, backstresslast_->at(gp));
      bracket.Update(fac_bracket_1, N_tilde, 1.0);

      // dNds_csb = - 3/2 . 1/qbar_tilde . 1/(1 + 3/2 . g) . bracket (50)
      LINALG::Matrix<NUM_STRESS_3D, 1> dNds_csb(false);
      double fac_dNds_csb = -3.0 / 2.0 * (1.0 / qbar_tilde) * (1.0 / (1.0 + 3.0 / 2.0 * g));
      dNds_csb.Update(fac_dNds_csb, bracket);

      // c_s_tilde = c_s =  - k_s - 2 G . N_tilde . c_strainbar
      //       + (3 G / qbar_tilde) . ( (Dgamma / omega) / (1 + 3/2 . g)) . dNds_csb
      double fac_c_s1 = (-2.0) * G * c_strainbar;
      c_s.Update((-1.0), k_s_tilde, fac_c_s1, N_tilde);
      double fac_c_s2 = 3.0 * G / qbar_tilde * (Dgamma / omega) / (1.0 + 3.0 / 2.0 * g);
      c_s.Update(fac_c_s2, dNds_csb, 1.0);

      // ---------------------- 4.) calculate (c_s_tilde - c_beta) (49)

      // c_s - c_b = - 1/(1 + 3/2 . g) . [ k_s - k_b + g . N_tilde : (k_s - k_b)
      //               . N_tilde ] -
      //             - [ 2 G . c_strainbar + Hkin_rec /
      //               (1 + Hkin_rec . Dgamma)^2 . c_R ] . N_tilde +
      //             + 1/(1 + 3/2 . g) . Hkin_rec / (1 + Hkin_rec . Dgamma)^2
      //               . [ beta_n + g . (N_tilde : beta_n) . N_tilde ] . c_R
      //
      // c_s - c_b = - fac_csb . [ k_s - k_b + g . N_ksb . N_tilde ]
      //             - fac_csb_1 . N_tilde +
      //             + fac_csb_2 . [ beta_n + g . Nbetaold . N_tilde ]
      double fac_csb = 1.0 / (1.0 + 3.0 / 2.0 * g);
      double fac_csb_1 = fac_csb * 2.0 * G * c_strainbar +
                         Hkin_rec * c_R / ((1.0 + Hkin_rec * Dgamma) * (1.0 + Hkin_rec * Dgamma));
      double fac_csb_2 =
          fac_csb * Hkin_rec * c_R / ((1.0 + Hkin_rec * Dgamma) * (1.0 + Hkin_rec * Dgamma));
      c_s_b.Update((-fac_csb), k_s_tilde, fac_csb, k_b);
      c_s_b.Update((-fac_csb * g * N_ksb), N_tilde, 1.0);
      c_s_b.Update((-fac_csb_1), N_tilde, 1.0);
      c_s_b.Update(fac_csb_2, backstresslast_->at(gp), 1.0);
      c_s_b.Update((fac_csb_2 * g * Nbetaold), N_tilde, 1.0);

      // ------ 5.) calculate c_beta = c_s_tilde - (c_s_tilde - c_beta)

      // c_b = c_s - (c_s - c_b) = c_s - c_s_b
      c_beta.Update(1.0, c_s, (-1.0), c_s_b);

      // --------------------------------------- 6.) calculate c_D (52)

      // c_D = - k_D + y . c_strainbar + (Dgamma / omega) .
      //         dy/ds_tilde : [ - k_s - 2 G . N_tilde . c_strainbar
      //          + 3 G / qbar_tilde . (Dgamma / omega) / (1 + 3/2 . g) . dNds_csb ]
      //
      // c_D = - k_D + y . c_strainbar + (Dgamma / omega) . dy/ds_tilde : matrix_c_D

      // matrix_c_D = [ - k_s - 2 G . N_tilde . c_strainbar
      //                + 3 G / qbar_tilde . (Dgamma / omega) / (1 + 3/2 . g) . dNds_csb ]
      LINALG::Matrix<NUM_STRESS_3D, 1> matrix_c_D(false);
      matrix_c_D.Update((-1.0), k_s_tilde, (-2.0 * G * c_strainbar), N_tilde);
      double fac_c_D = 3.0 * G / qbar_tilde * (Dgamma / omega) / (1.0 + 3.0 / 2.0 * g);
      matrix_c_D.Update(fac_c_D, dNds_csb, 1.0);

      // dyds_matrixcD = dy_ds_tilde : matrix_c_D
      double dyds_matrixcD = 0.0;
      for (int i = 0; i < 6; ++i) dyds_matrixcD += dy_dstilde(i, 0) * matrix_c_D(i, 0);

      // c_D = - k_D + y . c_strainbar + (Dgamma / omega) dy_ds_matrixcD
      c_D = -k_D + y * c_strainbar + (Dgamma / omega) * dyds_matrixcD;

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
      dserror(
          "INadmissible value of integrity: omega = %-14.8E in ele = %4d!"
          " \n Omega has to be greater than zero!",
          omega, EleGID);

    // update damage variable damage_{n+1}
    damage = 1.0 - omega;

    // --> damaged isotropic hardening variable has newest value (see L678)

    // get damaged pressure
    double p = omega * p_tilde;

    // deviatoric stress
    // s_{n+1} = omega_{n+1} . s^tilde_{n+1}
    devstress.Update(omega, devstress_tilde);

    // total stress
    // sigma_{n+1} = s_{n+1} + p_{n+1} . id2
    // pressure/volumetric stress no influence due to plasticity
    Stress(p, devstress, *stress);

    // ------------------------------------------------- update strains

    // strain^e_{n+1} = strain^(e,trial)_{n+1} - Dgamma . N_tilde
    strain_e.Update(1.0, trialstrain_e, (-Dgamma), N_tilde);

    // strain^p_{n+1} = strain^p_n + Dgamma . N_tilde
    // or alternatively
    //   strain^p_{n+1} = strain_{n+1} - strain^e_{n+1}
    strain_p.Update(Dgamma, N_tilde, 1.0);

    // compute converged engineering strain components (Voigt-notation)
    for (int i = 3; i < 6; ++i) strain_e(i) *= 2.0;
    for (int i = 3; i < 6; ++i) strain_p(i) *= 2.0;

    // pass the current plastic strains to the element (for visualisation)
    plstrain.Update(strain_p);

    // ------------------------------------------------- update history
    // plastic strain
    strainplcurr_->at(gp) = strain_p;

    // backstress
    backstresscurr_->at(gp) = beta;

    // accumulated plastic strain
    strainbarplcurr_->at(gp) = strainbar_p;

    // update damaged isotropic hardening variable
    isohardvarcurr_->at(gp) = Rplast;

    // update damage variable damage_{n+1}
    damagecurr_->at(gp) = damage;

#ifdef DEBUGMATERIAL
    std::cout << "end strain_p\n " << strain_p << std::endl;
    std::cout << "end strainplcurr_->at(gp)\n " << strainplcurr_->at(gp) << std::endl;
#endif  // ifdef DEBUGMATERIAL

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
    devstress.Update(omega, devstress_tilde);

    // result vectors of time step n+1 = omega . trial state vectors
    // sigma^e_n+1 = omega . sigma^(e,trial)_n+1
    //             = omega . (s^{trial}_{n+1} + p . id2)
    Stress(p, devstress, *stress);

    // total strains
    // strain^e_{n+1} = strain^(e,trial)_{n+1}
    // compute converged engineering strain components (Voigt-notation)
    strain_e.Update(trialstrain_e);
    for (int i = 3; i < 6; ++i) strain_e(i) *= 2.0;

    // no plastic yielding
    Dgamma = 0.0;

    // pass the current plastic strains to the element (for visualisation)
    // compute converged engineering strain components (Voigt-notation)
    for (int i = 3; i < 6; ++i) strain_p(i) *= 2.0;
    plstrain.Update(strain_p);

    // --------------------------------------------------------- update history
    // constant values for
    //  - plastic strains
    //  - accumulated (un)damaged plastic strains
    //  - back stress
    //    (--> relative stress)
    //  - stress
    // as current history vectors are set to zero in Update(), the old values
    // need to be set instead, otherwise no constant plastic values are possible
    strainplcurr_->at(gp) = strainpllast_->at(gp);
    strainbarplcurr_->at(gp) = strainbarpllast_->at(gp);
    backstresscurr_->at(gp) = backstresslast_->at(gp);
    isohardvarcurr_->at(gp) = isohardvarlast_->at(gp);
    damagecurr_->at(gp) = damagelast_->at(gp);

  }  // elastic step

  //---------------------------------------------------------------------------
  // --------------------------------- consistent elastoplastic tangent modulus
  //---------------------------------------------------------------------------

  // using an associative flow rule: C_ep is symmetric
  // ( generally C_ep is nonsymmetric )
  // if Phi^trial = 0: two tangent stress-strain relations exist
  // plastic loading --> C == C_ep
  if (Dgamma > 0.0) heaviside = 1.0;
  // (damaged) elastic unloading --> C == 1/(1 - D_{n+1})C_e
  else
    heaviside = 0.0;

  // using an associative flow rule: C_ep is symmetric
  // ( generally C_ep is nonsymmetric )
  SetupCmatElastoPlasticFullLemaitre(*cmat, N_tilde, *stress, heaviside, Dgamma, s_N, g, h_alg, G,
      dkappa_dR, bulk, Hkin, Hkin_rec, Nbetaold, gp, qbar_tilde, y, dy_dsigma_tilde, b_NbetaoldN);

#ifdef DEBUGMATERIAL
  std::cout << "Nach Setup Cep\n" << std::endl;
  std::cout << " Dgamma " << Dgamma << std::endl;
  std::cout << " G " << G << std::endl;
  std::cout << " q " << q << std::endl;
  std::cout << " flow vector " << Nbar << std::endl;
  std::cout << " heaviside " << heaviside << std::endl;
  std::cout << "--> cmat " << cmat << std::endl;
#endif  // #ifdef DEBUGMATERIAL

  // ------------------------------- return plastic strains for post-processing
  params.set<LINALG::Matrix<MAT::NUM_STRESS_3D, 1>>("plglstrain", plstrain);

  return;

}  // EvaluateFullLemaitre()


/*----------------------------------------------------------------------*
 | computes stress tensor                                    dano 11/13 |
 *----------------------------------------------------------------------*/
void MAT::Damage::Stress(const double p,                // volumetric stress
    const LINALG::Matrix<NUM_STRESS_3D, 1>& devstress,  // deviatoric stress tensor
    LINALG::Matrix<NUM_STRESS_3D, 1>& stress            // 2nd PK-stress
)
{
  // total stress = deviatoric + hydrostatic pressure . I
  // sigma = s + p . I
  stress.Update(devstress);
  for (int i = 0; i < 3; ++i) stress(i) += p;

}  // Stress()


/*----------------------------------------------------------------------*
 | compute relative deviatoric stress tensor                 dano 08/11 |
 *----------------------------------------------------------------------*/
void MAT::Damage::RelStress(
    const LINALG::Matrix<NUM_STRESS_3D, 1>& devstress,  // deviatoric stress tensor
    const LINALG::Matrix<NUM_STRESS_3D, 1>& beta,       // back stress tensor
    LINALG::Matrix<NUM_STRESS_3D, 1>& eta               // relative stress
)
{
  // relative stress = deviatoric - back stress
  // eta = s - beta
  eta.Update(1.0, devstress, (-1.0), beta);

}  // RelStress()


/*----------------------------------------------------------------------*
 | computes isotropic elasticity tensor in matrix notion     dano 04/11 |
 | for 3d                                                               |
 *----------------------------------------------------------------------*/
void MAT::Damage::SetupCmat(LINALG::Matrix<NUM_STRESS_3D, NUM_STRESS_3D>& cmat)
{
  // get material parameters
  // Young's modulus (modulus of elasticity)
  double young = params_->youngs_;
  // Poisson's ratio
  double nu = params_->poissonratio_;

  // isotropic elasticity tensor C in Voigt matrix notation
  // e.g. see lecture notes of NiliFEM Script (3.134)
  // C = 2G I + lambda (id2 \otimes id2) with I == I_s
  //   = 2G I_d + bulk (id2 \otimes id2)
  //
  //                       [ 1-nu     nu     nu |          0    0    0 ]
  //                       [        1-nu     nu |          0    0    0 ]
  //           E           [               1-nu |          0    0    0 ]
  //   C = --------------- [ ~~~~   ~~~~   ~~~~   ~~~~~~~~~~  ~~~  ~~~ ]
  //       (1+nu)*(1-2*nu) [                    | (1-2*nu)/2    0    0 ]
  //                       [                    |      (1-2*nu)/2    0 ]
  //                       [ symmetric          |           (1-2*nu)/2 ]
  //
  const double mfac = young / ((1.0 + nu) * (1.0 - 2.0 * nu));

  // clear the material tangent
  cmat.Clear();
  // write non-zero components
  cmat(0, 0) = mfac * (1.0 - nu);
  cmat(0, 1) = mfac * nu;
  cmat(0, 2) = mfac * nu;
  cmat(1, 0) = mfac * nu;
  cmat(1, 1) = mfac * (1.0 - nu);
  cmat(1, 2) = mfac * nu;
  cmat(2, 0) = mfac * nu;
  cmat(2, 1) = mfac * nu;
  cmat(2, 2) = mfac * (1.0 - nu);
  // ~~~
  cmat(3, 3) = mfac * 0.5 * (1.0 - 2.0 * nu);
  cmat(4, 4) = mfac * 0.5 * (1.0 - 2.0 * nu);
  cmat(5, 5) = mfac * 0.5 * (1.0 - 2.0 * nu);

}  // SetupCmat()


/*----------------------------------------------------------------------*
 | computes isotropic damaged elastoplastic tensor in        dano 05/11 |
 | matrix notion for 3d                                                 |
 *----------------------------------------------------------------------*/
void MAT::Damage::SetupCmatElastoPlastic(
    LINALG::Matrix<NUM_STRESS_3D, NUM_STRESS_3D>& cmat,  // elasto-plastic tangent modulus (out)
    int eleID,                                           // current element ID
    double Dgamma,                                       // plastic multiplier
    double G,                                            // shear modulus
    double bulk,                                         // bulk modulus
    double p_tilde,                                      // undamaged pressure
    double q_tilde,                         // undamaged trial von Mises equivalent stress
    double energyrelrate,                   // damage energy release rate
    double Ytan,                            // derivative of engergy release rate Ytan w.r.t. Dgamma
    double sigma_y,                         // current yield stress
    double Hiso,                            // isotropic hardening modulus
    LINALG::Matrix<NUM_STRESS_3D, 1> Nbar,  // unit flow vector
    int gp,                                 // current Gauss point
    bool damevolution,                      // flag indicating if damage evolves or not
    double heaviside                        // Heaviside function
)
{
  // damage threshold is still not passed, i.e. use undamaged material tangent
  if (damevolution == false)
  {
    // incremental constitutive function for the stress tensor
    // sigma_n+1 = [ cmat - (Dgamma 6 G^2/q) I_d ] : strain_n+1^{e,trial}
    // consistent tangent operator
    // C^{ep} := dsigma_n+1 / dstrain_n+1^{e,trial}

    // depending on the flow vector Cmat_ep can be a fully-occupied matrix

    // C^{ep} = C^e - ( H^ . Dgamma . 6 . G^2 ) / q . I_d +
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
    LINALG::Matrix<NUM_STRESS_3D, 1> id2(true);
    for (int i = 0; i < 3; i++) id2(i) = 1.0;

    // set Cartesian identity 4-tensor in 6-Voigt matrix notation
    // this is fully 'contra-variant' identity tensor, ie I^{ABCD}
    // REMARK: rows are stress-like 6-Voigt
    //         columns are stress-like 6-Voigt
    LINALG::Matrix<NUM_STRESS_3D, NUM_STRESS_3D> id4sharp(true);
    for (int i = 0; i < 3; i++) id4sharp(i, i) = 1.0;
    for (int i = 3; i < 6; i++) id4sharp(i, i) = 0.5;

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
      epfac = (-1.0) * heaviside * Dgamma * 6.0 * G * G / q_tilde;
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

    if (q_tilde != 0.0)
    {
      // loop strains (columns)
      for (int k = 0; k < 6; ++k)
      {
        // ---------------------------------------------------------- tangent
        // loop stresses (rows)
        for (int i = 0; i < 6; ++i)
        {
          epfac3 = heaviside * 6.0 * G * G * (Dgamma / q_tilde - 1.0 / (3.0 * G + Hiso));
          // here: Nbar = s^{trial}_{n+1} / || s^{trial}_{n+1} ||
          cmat(i, k) += epfac3 * Nbar(i) * Nbar(k);
        }  // end rows, loop i
      }    // end columns, loop k
    }      // (q != 0.0)

    // complete material tangent C_ep available

#ifdef DEBUGMATERIAL
    if (Dgamma != 0)
    {
      std::cout << "Ende SetupCmatElastPlast" << std::endl;
      std::cout << "Cep\n"
                << " Dgamma " << Dgamma << std::endl;
      std::cout << " G " << G << std::endl;
      std::cout << " q " << q << std::endl;
      std::cout << " Nbar " << Nbar << std::endl;
      std::cout << " heaviside " << heaviside << std::endl;
      std::cout << " epfac " << epfac << std::endl;
      std::cout << " epfac1 " << epfac1 << std::endl;
      std::cout << " cmat " << cmat << std::endl;
    }
#endif  // #ifdef DEBUGMATERIAL
  }     // (damevolution == false)

  // material tangent differs for case damage or not
  // if no damage: use standard tangent of purely plastic behaviour
  else  // (damevolution == true)
  {
#ifdef DEBUGMATERIAL
    if (gp == 0) std::cout << "damage evolution takes place in eleID = " << eleID << endl;
#endif  // #ifdef DEBUGMATERIAL

    // incremental constitutive function for the stress tensor
    // consistent tangent operator
    // C^{ep} := dsigma_n+1 / dstrain_n+1

    // depending on the flow vector Cmat_ep can be a fully-occupied matrix

    // ---------------------------------------------------------- Heaviside

    // if plastic loading:   heaviside = 1.0 --> use C^{ep}
    // if elastic unloading: heaviside = 0.0 --> use C^e

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
    LINALG::Matrix<NUM_STRESS_3D, 1> id2(true);
    for (int i = 0; i < 3; i++) id2(i) = 1.0;

    // set Cartesian identity 4-tensor in 6-Voigt matrix notation
    // this is fully 'contra-variant' identity tensor, i.e. I^{ABCD}
    // REMARK: rows are stress-like 6-Voigt
    //         columns are stress-like 6-Voigt
    LINALG::Matrix<NUM_STRESS_3D, NUM_STRESS_3D> id4sharp(true);
    for (int i = 0; i < 3; i++) id4sharp(i, i) = 1.0;
    for (int i = 3; i < 6; i++) id4sharp(i, i) = 0.5;

    // ------------------------------------- extract current history values
    // integrity omega_{n+1}
    double omega = 1.0 - damagecurr_->at(gp);

    // ----------------------------------------------------- damaged elastic term
    // C^{ep} = (1 - D_{n+1}) . C^e = omega_{n+1} . C^e
    //        = omega_{n+1} . 2G . I_d + omega_{n+1} . bulk . id2 \otimes id2
    // add standard isotropic elasticity tensor C^e first
    if (heaviside == 0)
    {
      SetupCmat(cmat);
      cmat.Scale(omega);
    }
    else  // (heaviside == 1)
    {
      // C^{ep} = a . I_d + b . Nbar \otimes Nbar + c . Nbar \otimes id2
      //          + d . id2 \otimes Nbar + e . id2 \otimes id2

      // ----------------------------------- extract some material parameters

      // get material parameters
      // damage evolution law denominator r
      double damden = params_->damden_;
      // damage evolution law exponent s
      double damexp = params_->damexp_;

      // ---------------------------------------------------------- plastic terms

      // check if omega is admissible
      if (omega < 1.0e-20) dserror("Omega has to be greater than zero! omega = %-14.8E\n", omega);

      // be aware: in the 1st equilibrium (Newton) iteration (i=0) C^{ep} is
      // indeterminate due to Dgamma == 0.0
      // a small perturbation Dgamma=1e-08 is used (instead of using a limiting
      // procedure) which was proposed by de Souza Neto
      if (Dgamma == 0) Dgamma = 1.0e-08;

      // ------------------------ factors required for the elasto-plastic tangent

      // some constants
      double aux = -energyrelrate / damden;
      double auxb = (q_tilde - sigma_y) / (3.0 * G);
      // PhiT = qtilde - sigma_y(Rplast_{n+1})
      // be careful: NOT Phi_trial which was calculated using sigma_y(Rplast_n)
      double PhiT = q_tilde - sigma_y;

      // ---------------------------------------- derivatives w.r.t. Dgamma
      // domega/dDgamma
      double Domega = (3.0 * G + omega * Hiso) / PhiT;
      // domega/dq_tilde . dq_tilde/dDgamma
      double DomegaDq_tilde = -omega / PhiT;

      // derviative of residual function dF/dDgamma
      double ResTan = Domega - Hiso / (3.0 * G) * std::pow(aux, damexp) -
                      auxb * damexp * Ytan / damden * std::pow(aux, (damexp - 1));
      // derivative of residual function dF/dp_tilde . dp_tilde/dDgamma
      // DResDp_tilde = s . (q_tilde - sigma_y)/(3G) . std::pow((-Y/r), s-1) . p_tilde /(r . bulk)
      double DResDp_tilde =
          damexp * auxb * std::pow(aux, (damexp - 1.0)) * p_tilde / (damden * bulk);
      // derivative of residual function F w.r.t. q_tilde
      double DResDq_tilde = DomegaDq_tilde + std::pow(aux, damexp) / (3.0 * G);

      // --------- calculate the constants according to de Souza Neto (12.52)

      // a1 = - DResDq_tilde / ResTan
      //    = [ omega / (q_tilde - sigma_y) - (-energyrelrate/damden)^damexp
      //                                      / (3.0 * G) ] / ResTan;
      // or alternatively
      //   a1 = [1 - Dgamma / (omega . omega) . (-Y/r)^s] . omega
      //          / [ (q_tilde - sigma_y) . ResTan]
      double a1 = -DResDq_tilde / ResTan;

      // a2 = - DResDp_tilde / ResTan
      //    = - s . p_tilde . (q_tilde - sigma_y) / (3G . r . bulk . ResTan) .(-Y/r)^{s-1};
      double a2 = 0.0;
      a2 = -DResDp_tilde / ResTan;

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
      // b = 2G . (a1 . Hiso . omega + a4 . sigma_y - sigma_y . omega / q_tilde
      double b = 0.0;

      if (q_tilde != 0.0)
      {
        a = (2.0 * G * sigma_y * omega) / q_tilde;
        b = 2.0 * G * (a1 * Hiso * omega + a4 * sigma_y - sigma_y * omega / q_tilde);
      }
      // c = bulk . (a2 . Hiso . omega + a3 . sigma_y) / sqrt(3/2)
      double c = 0.0;
      c = bulk * (a2 * Hiso * omega + a3 * sigma_y) / sqrt(3.0 / 2.0);
      // d = p_tilde . 2G . sqrt(3/2) . a4
      double d = 0.0;
      d = (p_tilde * 2.0 * G * sqrt(3.0 / 2.0) * a4);
      // e = bulk . (omega + p_tilde . a3)
      double e = 0.0;
      e = bulk * (omega + p_tilde * a3);

      // ------------------------------- assemble elasto-plastic material tangent

      // empty consistent tangent operator
      cmat.Clear();
      // constitutive tensor
      // I_d = id4sharp - 1/3 Id \otimes Id
      // contribution: Id4^#
      cmat.Update(a, id4sharp);
      // contribution: Id \otimes Id
      cmat.MultiplyNT((a / (-3.0)), id2, id2, 1.0);
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
      std::cout << "Cep\n"
                << " Dgamma " << Dgamma << std::endl;
      std::cout << " G " << G << std::endl;
      std::cout << " q_tilde " << q_tilde << std::endl;
      std::cout << " Nbar " << Nbar << std::endl;
      std::cout << " heaviside " << heaviside << std::endl;
      std::cout << " epfac " << epfac << std::endl;
      std::cout << " epfac1 " << epfac1 << std::endl;
      std::cout << " cmat " << cmat << std::endl;
    }
#endif  // #ifdef DEBUGMATERIAL

  }  // damage evolves: (damevolution == true)

}  // SetupCmatElastoPlastic()


/*----------------------------------------------------------------------*
 | computes isotropic elastoplastic tensor for full Lemaitre dano 11/13 |
 | model in matrix notion for 3d including ductile damage               |
 *----------------------------------------------------------------------*/
void MAT::Damage::SetupCmatElastoPlasticFullLemaitre(
    LINALG::Matrix<NUM_STRESS_3D, NUM_STRESS_3D>& cmat,  // elasto-plastic tangent modulus (out)
    LINALG::Matrix<NUM_STRESS_3D, 1> N_tilde,            // flow vector
    LINALG::Matrix<NUM_STRESS_3D, 1>& stress,            // stress
    double heaviside,                                    // Heaviside-function
    double Dgamma,                                       // plastic multiplier
    double s_N,  // s_N = 2 G . (Dgamma / omega)^2 . dy/ds_tilde : N_tilde
    double g,    // g = (2 G . (Dgamma / omega) + Hkin . Dgamma / (1 + Hkin_rec . Dgamma) / q_tilde
    double h_alg,       // see definition of expression above
    double G,           // shear modulus
    double dkappa_dR,   // derivative of the hardening curve
    double bulk,        // bulk modulus
    double Hkin,        // kinematic hardening variable 1 (kinematic hardening modulus)
    double Hkin_rec,    // kinematic hardening variable 2 (saturation effect)
    double Nbetaold,    // N_beta = Hkin_rec . N_tilde : beta_n, CARO:Nbetaold =N_tilde : beta_n
    int gp,             // current Gauss point
    double qbar_tilde,  // effective trial stress ^tilde
    double y,           // (-Y / r)^s
    LINALG::Matrix<NUM_STRESS_3D, 1> dy_dsigma_tilde,
    LINALG::Matrix<NUM_STRESS_3D, 1> b_NbetaoldN  // beta_n - 2/3 . (N_tilde : beta_n) . N_tilde
)
{
  // incremental constitutive function for the stress tensor
  // consistent tangent operator
  // C^{ep} := dsigma_n+1 / dstrain_n+1

  // depending on the flow vector Cmat_ep can be a fully-occupied matrix

  // ---------------------------------------------------------- Heaviside

  // if plastic loading:   heaviside = 1.0 --> use C^{ep}
  // if elastic unloading: heaviside = 0.0 --> use C^e

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
  LINALG::Matrix<NUM_STRESS_3D, 1> id2(true);
  for (int i = 0; i < 3; i++) id2(i) = 1.0;

  // set Cartesian identity 4-tensor in 6-Voigt matrix notation
  // this is fully 'contra-variant' identity tensor, i.e. I^{ABCD}
  // REMARK: rows are stress-like 6-Voigt
  //         columns are stress-like 6-Voigt
  LINALG::Matrix<NUM_STRESS_3D, NUM_STRESS_3D> id4sharp(true);
  for (int i = 0; i < 3; i++) id4sharp(i, i) = 1.0;
  for (int i = 3; i < 6; i++) id4sharp(i, i) = 0.5;

  // ------------------------------------- extract current history values
  // integrity omega_{n+1}
  double omega = 1.0 - damagecurr_->at(gp);

  // ------------------------------------------ elastic undamaged tangent
  // C^e = 2G . I_d + bulk . id2 \otimes id2
  SetupCmat(cmat);

  if (heaviside == 0)
  {
    // ------------------------------------------------ damaged elastic tangent

    // C^e_D = (1 - D_{n+1}) . C^e = omega_{n+1} . C^e
    //        = omega_{n+1} . 2G . I_d + omega_{n+1} . bulk . id2 \otimes id2
    // add standard isotropic elasticity tensor C^e first
    cmat.Scale(omega);
  }
  else  // (heaviside == 1)
  {
    // ------------------------------------------ damaged elastoplastic tangent

    // consistent/algorithmic damamged elastoplastic tangent operator (70)
    // C^{ep} = omega . C_tilde^{ep} - stress_tilde \otimes .
    //          . { Dgamma / omega . C_tilde^{ep} : dy/dsigma_tilde
    //              + 2/3 . y . N_tilde
    //              - y / 3 G . [ dkappa_dR + (3/2 . Hkin - Hkin_rec .
    //                . N_tilde : beta_n ) / (1 + Hkin_rec . Dgamma)^2 ]
    //                 . n_alg / h_alg }

    // consistent tangent operator is achieved due to variation of effective,
    // undamaged stress sigma_tilde (68)
    // C^{~,ep} = C^e - (2 G)^2 . (Dgamma / omega) / (1 + 3/2 . g) . dPhi_dsigma_tilde_square
    //            - 2/3 . N_tilde \otimes (2 G . N_tilde)
    //            + { 2/3 . [ dkappa_dR + 3/2 . Hkin - Hkin_rec . (N_tilde : beta_n)
    //                        / (1 + Hkin_rec . Dgamma)^2 ] . N_tilde
    //                 - 3 G / qbar_tilde . (Dgamma / omega) / (1 + 3/2 . g) .
    //                     Hkin_rec / (1 + Hkin_rec . Dgamma)^2 .
    //                     [ beta_n - 2/3 . (N_tilde : beta_n) . N_tilde ]
    //                 } \otimes (n_alg / h_alg)

    // stress_tilde = stress / omega = sigma / omega
    LINALG::Matrix<NUM_STRESS_3D, 1> stress_tilde(false);
    stress_tilde.Update((1.0 / omega), stress);

    // --------------------------- linearisation of delta_R = delta_gamma

    // Lin[ R ] = (1 - D) . Lin[ strainbar ] - Dgamma / omega . Lin[ D ]  (60e)
    //          = (n^alg : Lin[ strain ]) / h^alg                         (64)
    //
    // n_alg = [ omega - y . Dgamma / omega + 2 G . (Dgamma / omega)^2 .
    //           dy_ds_tilde : N_tilde ] . 2 G . N_tilde
    //         - 3 G . (Dgamma / omega)^2 . (C^e - (2 G)^2 . Dgamma / omega
    //           / (1 + 3/2 . g) . dPhi_dsigma_tilde_square ] : dy_dsigma_tilde
    //
    // n_alg = [ omega - y . Dgamma / omega + s_N ] . 2 G . N_tilde
    //         - 3 G . (Dgamma / omega)^2 . (C^e - (2 G)^2 . Dgamma / omega
    //           / (1 + 3/2 . g) . dPhi_dsigma_tilde_square ] : dy_dsigma_tilde

    // dPhi_dsigma_tilde_square = (3/2 . I_d - N_tilde \otimes N_tilde) / qbar_tilde
    LINALG::Matrix<NUM_STRESS_3D, NUM_STRESS_3D> dPhi_dsigma_tilde_square(true);
    double fac_dPhidsigma = 0.0;
    if (qbar_tilde != 0)
    {
      fac_dPhidsigma = 3.0 / 2.0 / qbar_tilde;
      // I_d = id4sharp - 1/3 Id \otimes Id
      // contribution: Id4^#
      dPhi_dsigma_tilde_square.Update(fac_dPhidsigma, id4sharp);
      // contribution: Id \otimes Id
      double fac_dPhidsigma_1 = fac_dPhidsigma / (-3.0);
      dPhi_dsigma_tilde_square.MultiplyNT(fac_dPhidsigma_1, id2, id2, 1.0);
      // loop strains (columns)
      for (int k = 0; k < 6; ++k)
      {
        // ---------------------------------------------------------- tangent
        // loop stresses (rows)
        for (int i = 0; i < 6; ++i)
        {
          // (- 1 / qbar_tilde) . N_tilde \otimes N_tilde
          dPhi_dsigma_tilde_square(i, k) += -1.0 / qbar_tilde * N_tilde(i) * N_tilde(k);
        }  // end rows, loop i
      }    // end columns, loop k
    }      // (qbar_tilde != 0.0)

    // n_alg = [ omega - y . Dgamma / omega + s_N ] . 2 G . N_tilde
    //         - 3 G . (Dgamma / omega)^2 . [ C^e - (2 G)^2 . Dgamma / omega
    //           / (1 + 3/2 . g) . dPhi_dsigma_tilde_square ] : dy_dsigma_tilde

    // n_alg = [ omega - y . Dgamma / omega + s_N ] . 2 G . N_tilde
    //         - 3 G . (Dgamma / omega)^2 . C^e : dy_dsigma_tilde
    //         - 3 G . (Dgamma / omega)^2 . - (2 G)^2 . Dgamma / omega
    //           / (1 + 3/2 . g) . dPhi_dsigma_tilde_square : dy_dsigma_tilde

    // n_alg = [ omega - y . Dgamma / omega + s_N ] . 2 G . N_tilde
    //
    // with s_N = = 2 G . (Dgamma / omega)^2 . dy/ds_tilde : N_tilde
    double fac_n_alg_1 = (omega - y * (Dgamma / omega) + s_N) * 2.0 * G;
    LINALG::Matrix<NUM_STRESS_3D, 1> n_alg(false);
    n_alg.Update(fac_n_alg_1, N_tilde);

    // n_alg += - 3 G . (Dgamma / omega)^2 . C^e : dy_dsigma_tilde
    // calculate C^e : dy_dsigma_tilde
    //          (6x6)  (6x1)
    LINALG::Matrix<NUM_STRESS_3D, 1> ce_dydsigmatilde(false);
    ce_dydsigmatilde.Multiply(cmat, dy_dsigma_tilde);
    double fac_n_alg_2 = -3.0 * G * (Dgamma / omega) * (Dgamma / omega);
    n_alg.Update(fac_n_alg_2, ce_dydsigmatilde, 1.0);

    // n_alg += 3 G . (Dgamma / omega)^3 . (2 G)^2
    //           / (1 + 3/2 . g) . dPhi_dsigma_tilde_square : dy_dsigma_tilde
    // dPhi_dsigma_tilde_square : dy_dsigma_tilde
    LINALG::Matrix<NUM_STRESS_3D, 1> dPhidsigma_dydsigma(false);
    dPhidsigma_dydsigma.Multiply(dPhi_dsigma_tilde_square, dy_dsigma_tilde);
    double fac_n_alg_3 =
        3.0 * G * std::pow((Dgamma / omega), 3) * std::pow((2 * G), 2) / (1 + 3.0 / 2.0 * g);
    n_alg.Update(fac_n_alg_3, dPhidsigma_dydsigma, 1.0);

    // dkappa_HHNb = dkappa_dR + [ 3/2 . Hkin - Hkin_rec . (N_tilde : beta_n) ]
    //                             / (1 + Hkin_rec . Dgamma)^2
    double dkappa_HHNb = dkappa_dR + (3.0 / 2.0 * Hkin - Hkin_rec * Nbetaold) /
                                         ((1.0 + Hkin_rec * Dgamma) * (1.0 + Hkin_rec * Dgamma));

    // b_NbetaoldN = beta_n - 2/3 . (N_tilde : beta_n) . N_tilde

    // C^{~,ep} = C^e
    LINALG::Matrix<NUM_STRESS_3D, NUM_STRESS_3D> Cep_tilde(false);
    Cep_tilde.Update(cmat);

    // C^{~,ep} += - (2 G)^2 . (Dgamma / omega) / (1 + 3/2 . g) . dPhi_dsigma_tilde_square
    //          += fac_Cep_tilde1 . dPhi_dsigma_tilde_square
    // - (2 G)^2 . (Dgamma / omega) / (1 + 3/2 . g)
    double fac_Cep_tilde1 = -(2.0 * G) * (2.0 * G) * (Dgamma / omega) / (1.0 + 3.0 / 2.0 * g);
    Cep_tilde.Update(fac_Cep_tilde1, dPhi_dsigma_tilde_square, 1.0);

    // C^{~,ep} += - 2/3 . 2 G . N_tilde \otimes  N_tilde)
    //          += fac_Cep_tilde2 . N_tilde \otimes  N_tilde)
    double fac_Cep_tilde2 = (-2.0) / 3.0 * 2.0 * G;
    Cep_tilde.MultiplyNT(fac_Cep_tilde2, N_tilde, N_tilde, 1.0);

    // prematrix_n_alg = 2/3 . dkappa_HHNb . N_tilde
    //                   - 3 G / qbar_tilde . (Dgamma / omega) / (1 + 3/2 . g) .
    //                  . Hkin_rec / (1 + Hkin_rec . Dgamma)^2 . b_NbetaoldN
    LINALG::Matrix<NUM_STRESS_3D, 1> prematrix_n_alg(false);

    // fac_prematrix_n_alg = - 3 G / qbar_tilde . (Dgamma / omega) / (1 + 3/2 . g)
    //                       . Hkin_rec / (1 + Hkin_rec . Dgamma)^2
    double fac_prematrix_n_alg = -3.0 * G / qbar_tilde * (Dgamma / omega) / (1.0 + 3.0 / 2.0 * g) *
                                 Hkin_rec / ((1.0 + Hkin_rec * Dgamma) * (1.0 + Hkin_rec * Dgamma));
    // prematrix_n_alg = 2/3 . dkappa_HHNb . N_tilde
    //                   + fac_prematrix_n_alg . b_NbetaoldN
    prematrix_n_alg.Update((2.0 / 3.0 * dkappa_HHNb), N_tilde);
    prematrix_n_alg.Update(fac_prematrix_n_alg, b_NbetaoldN, 1.0);

    // C^{~,ep} += + prematrix_n_alg \otimes (n_alg / h_alg)
    double fac_Cep_tilde3 = 1.0 / h_alg;
    Cep_tilde.MultiplyNT(fac_Cep_tilde3, prematrix_n_alg, n_alg, 1.0);

    // ------------- assemble consistent damaged elastoplastic material tangent

    // C^{ep} = omega . C^{~,ep}
    //          - stress_tilde \otimes cep_tilde_dydsigma_tilde
    //          - stress_tilde \otimes . 2/3 . y . N_tilde
    //          - stress_tilde \otimes . fac_cep_nalg . n_alg

    // C^{ep} = omega . C^{~,ep}
    cmat.Update(omega, Cep_tilde);

    // cep_tilde_dydsigma_tilde = (Dgamma / omega) . C^{~,ep} : dy/dsigma_tilde
    LINALG::Matrix<NUM_STRESS_3D, 1> cep_tilde_dydsigma_tilde(false);
    cep_tilde_dydsigma_tilde.Multiply((Dgamma / omega), Cep_tilde, dy_dsigma_tilde, 0.0);
    // C^{ep} += - stress_tilde \otimes cep_tilde_dydsigma_tilde
    cmat.MultiplyNT((-1.0), stress_tilde, cep_tilde_dydsigma_tilde, 1.0);

    // C^{ep} += - stress_tilde \otimes . (2/3 . y) . N_tilde }
    cmat.MultiplyNT((-2.0 / 3.0 * y), stress_tilde, N_tilde, 1.0);

    // C^{ep} += - stress_tilde \otimes (fac_cep_nalg . n_alg)
    double fac_cep_nalg = -y / (3 * G) * dkappa_HHNb / h_alg;
    cmat.MultiplyNT((-fac_cep_nalg), stress_tilde, n_alg, 1.0);

  }  // (heaviside == 1)

}  // SetupCmatElastoPlasticFullLemaitre()


/*----------------------------------------------------------------------*
 | return derivative of piecewise linear function for the    dano 09/13 |
 | yield stress, i.e. isotropic hardening modulus at current            |
 | accumulated plastic strain                                           |
 *----------------------------------------------------------------------*/
double MAT::Damage::GetIsoHardAtStrainbarnp(const double strainbar_p  // current accumulated strain
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
  for (int i = 0; i < samplenumber; ++i)
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
        Hiso =
            (sigma_y_ref[i] - sigma_y_ref[i - 1]) / (strainbar_p_ref[i] - strainbar_p_ref[i - 1]);
        continue;
      }
    }  // load is plastic, hardening can occur
  }    // loop over samples

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
    if (sigma_y_ref.size() != strainbar_p_ref.size()) dserror("Samples have to fit to each other!");

    // loop over samples
    for (int i = 0; i < samplenumber; ++i)
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
            sigma_y_interpol =
                sigma_y_ref[i - 1] + (strainbar_p - strainbar_p_ref[i - 1]) *
                                         (sigma_y_ref[i] - sigma_y_ref[i - 1]) /
                                         (strainbar_p_ref[i] - strainbar_p_ref[i - 1]);
          }  // current strains between strain^{i-1} and strain^i
        }    // plastic regime
      }
    }  // loop over all samples
  }    // samplenumber > 1

  // return current yield stress
  return sigma_y_interpol;

}  // GetSigmaYAtStrainbarnp()


/*---------------------------------------------------------------------*
 | return names of visualization data (public)              dano 09/13 |
 *---------------------------------------------------------------------*/
void MAT::Damage::VisNames(std::map<std::string, int>& names)
{
  std::string accumulatedstrain = "accumulatedstrain";
  names[accumulatedstrain] = 1;  // scalar

  std::string isohardeningvar = "isohardeningvar";
  names[isohardeningvar] = 1;  // scalar

  std::string damage = "damage";
  names[damage] = 1;  // scalar

}  // VisNames()


/*---------------------------------------------------------------------*
 | return visualization data (public)                       dano 09/13 |
 *---------------------------------------------------------------------*/
bool MAT::Damage::VisData(const std::string& name, std::vector<double>& data, int numgp, int eleID)
{
  if (name == "accumulatedstrain")
  {
    if ((int)data.size() != 1) dserror("size mismatch");
    double temp = 0.0;
    for (int iter = 0; iter < numgp; iter++) temp += AccumulatedStrain(iter);
    data[0] = temp / numgp;
  }

  if (name == "isohardeningvar")
  {
    if ((int)data.size() != 1) dserror("size mismatch");
    double temp = 0.0;
    for (int iter = 0; iter < numgp; iter++) temp += IsotropicHardeningVariable(iter);
    data[0] = temp / numgp;
  }

  if (name == "damage")
  {
    if ((int)data.size() != 1) dserror("size mismatch");
    double temp = 0.0;
    for (int iter = 0; iter < numgp; iter++) temp += IsotropicDamage(iter);
    data[0] = temp / numgp;
  }

  return true;
}  // VisData()


/*----------------------------------------------------------------------*/
