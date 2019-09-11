/*----------------------------------------------------------------------*/
/*! \file
\brief Contains the functions to establish local material law /
       stress-strain law for isotropic material for a 3D hex element
       following an associative von Mises plasticity and a linear thermoelastic
       material law (St.Venant Kirchhoff).

       Combination of the two material model
       - PlastiLinElastic
       - ThermoStVenantKirchhoff

        Thermo-elasto-plastic material with mixed hardening
        - so far:
          - linear isotropic hardening with constant hardening modulus Hiso
            - yield stress sigma_y = sigma_y0 + Hiso . strainplbar
          - linear kinematic hardening with constant hardening modulus Hkin
            - 'linear Armstrong-Frederick kinematic hardening'
            - describing the Bauschinger effect via \f Hkin \,=\, const.\f

       small strains including temperature dependency

       strain-energy potential
       \f \rho \psi \,=\,\rho (\psi^e \,+\, \psi^p)
          \,=\, mu strain^e : strain^e + 1/2 . lambda (strain^e : I)^2 +
                + m (T - T_0) strain^e : I - rho C_V (T ln T/T_o - (T - T_0))
                + 1/2 strain^p : strain^p . (2/3 Hkin)
                + 1/2 strainbar^p : strainbar^p . Hiso \f

       example input line:
       MAT 1 MAT_Struct_ThrPlasticLinElast YOUNG 206.9 NUE 0.29 DENS 0.0
         THEXPANS 1.72e-5 INITTEMP 293.15 YIELD 0.45 ISOHARD 0.0 KINHARD 0.0
         SAMPLENUM 2 SIGMA_Y 0.1 0.2 EPSBAR_P 0.0 1.0 TOL 1.0e-6

\level 2

\maintainer Jan Schnabel
*/
/*----------------------------------------------------------------------*
 | headers                                                   dano 08/11 |
 *----------------------------------------------------------------------*/
#include "thermoplasticlinelast.H"
#include "matpar_bundle.H"

#include "../drt_lib/drt_globalproblem.H"

#include "../drt_tsi/tsi_defines.H"


/*----------------------------------------------------------------------*
 | constructor (public)                                      dano 08/11 |
 *----------------------------------------------------------------------*/
MAT::PAR::ThermoPlasticLinElast::ThermoPlasticLinElast(Teuchos::RCP<MAT::PAR::Material> matdata)
    : Parameter(matdata),
      youngs_(matdata->GetDouble("YOUNG")),
      poissonratio_(matdata->GetDouble("NUE")),
      density_(matdata->GetDouble("DENS")),
      thermexpans_(matdata->GetDouble("THEXPANS")),
      thetainit_(matdata->GetDouble("INITTEMP")),
      yield_(matdata->GetDouble("YIELD")),
      isohard_(matdata->GetDouble("ISOHARD")),
      kinhard_(matdata->GetDouble("KINHARD")),
      sigma_y_(*(matdata->Get<std::vector<double>>("SIGMA_Y"))),
      strainbar_p_ref_(*(matdata->Get<std::vector<double>>("EPSBAR_P"))),
      abstol_(matdata->GetDouble("TOL"))
{
}


/*----------------------------------------------------------------------*
 | is called in Material::Factory from ReadMaterials()       dano 08/11 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<MAT::Material> MAT::PAR::ThermoPlasticLinElast::CreateMaterial()
{
  return Teuchos::rcp(new MAT::ThermoPlasticLinElast(this));
}


MAT::ThermoPlasticLinElastType MAT::ThermoPlasticLinElastType::instance_;


/*----------------------------------------------------------------------*
 | is called in Material::Factory from ReadMaterials()       dano 08/11 |
 *----------------------------------------------------------------------*/
DRT::ParObject* MAT::ThermoPlasticLinElastType::Create(const std::vector<char>& data)
{
  MAT::ThermoPlasticLinElast* plastic = new MAT::ThermoPlasticLinElast();
  plastic->Unpack(data);
  return plastic;
}


/*----------------------------------------------------------------------*
 | constructor (public)                                      dano 08/11 |
 *----------------------------------------------------------------------*/
MAT::ThermoPlasticLinElast::ThermoPlasticLinElast() : params_(NULL) {}


/*----------------------------------------------------------------------*
 | constructor (public)                                      dano 08/11 |
 | called in ReadMaterials --> CreateMaterial                           |
 *----------------------------------------------------------------------*/
MAT::ThermoPlasticLinElast::ThermoPlasticLinElast(MAT::PAR::ThermoPlasticLinElast* params)
    : params_(params), plastic_step_(false)
{
}


/*----------------------------------------------------------------------*
 | pack (public)                                             dano 08/11 |
 *----------------------------------------------------------------------*/
void MAT::ThermoPlasticLinElast::Pack(DRT::PackBuffer& data) const
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

    AddtoPack(data, Dmech_->at(var));
    AddtoPack(data, Dmech_d_->at(var));

    AddtoPack(data, Incstrainpl_->at(var));
    AddtoPack(data, strainelrate_->at(var));
  }

  AddtoPack(data, plastic_step_);

  return;

}  // Pack()


/*----------------------------------------------------------------------*
 | unpack (public)                                           dano 08/11 |
 *----------------------------------------------------------------------*/
void MAT::ThermoPlasticLinElast::Unpack(const std::vector<char>& data)
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
        params_ = static_cast<MAT::PAR::ThermoPlasticLinElast*>(mat);
      else
        dserror("Type of parameter material %d does not fit to calling type %d", mat->Type(),
            MaterialType());
    }

  // history data
  int histsize;
  ExtractfromPack(position, data, histsize);

  // if system is not yet initialised, the history vectors have to be intialized
  if (histsize == 0) isinit_ = false;

  // unpack plastic history vectors
  strainpllast_ = Teuchos::rcp(new std::vector<LINALG::Matrix<NUM_STRESS_3D, 1>>);
  strainplcurr_ = Teuchos::rcp(new std::vector<LINALG::Matrix<NUM_STRESS_3D, 1>>);

  // unpack back stress vectors (for kinematic hardening)
  backstresslast_ = Teuchos::rcp(new std::vector<LINALG::Matrix<NUM_STRESS_3D, 1>>);
  backstresscurr_ = Teuchos::rcp(new std::vector<LINALG::Matrix<NUM_STRESS_3D, 1>>);

  strainbarpllast_ = Teuchos::rcp(new std::vector<double>);
  strainbarplcurr_ = Teuchos::rcp(new std::vector<double>);

  // unpack dissipation stuff
  Dmech_ = Teuchos::rcp(new std::vector<double>);
  Dmech_d_ = Teuchos::rcp(new std::vector<LINALG::Matrix<NUM_STRESS_3D, 1>>);

  Incstrainpl_ = Teuchos::rcp(new std::vector<LINALG::Matrix<NUM_STRESS_3D, 1>>);
  strainelrate_ = Teuchos::rcp(new std::vector<LINALG::Matrix<NUM_STRESS_3D, 1>>);

  for (int var = 0; var < histsize; ++var)
  {
    LINALG::Matrix<NUM_STRESS_3D, 1> tmp_vect(true);
    double tmp_scalar = 0.0;
    // vectors of last converged state are unpacked
    ExtractfromPack(position, data, tmp_vect);
    strainpllast_->push_back(tmp_vect);
    ExtractfromPack(position, data, tmp_vect);
    backstresslast_->push_back(tmp_vect);
    // scalar-valued vector of last converged state are unpacked
    ExtractfromPack(position, data, tmp_scalar);
    strainbarpllast_->push_back(tmp_scalar);
    ExtractfromPack(position, data, tmp_scalar);
    Dmech_->push_back(tmp_scalar);
    ExtractfromPack(position, data, tmp_vect);
    Dmech_d_->push_back(tmp_vect);
    ExtractfromPack(position, data, tmp_vect);
    Incstrainpl_->push_back(tmp_vect);
    ExtractfromPack(position, data, tmp_vect);
    strainelrate_->push_back(tmp_vect);

    // current vectors have to be initialised
    strainplcurr_->push_back(tmp_vect);
    backstresscurr_->push_back(tmp_vect);
    strainbarplcurr_->push_back(tmp_scalar);
  }

  plastic_step_ = false;
  int plastic_step;
  ExtractfromPack(position, data, plastic_step);

  // if it was already plastic before
  if (plastic_step != 0) plastic_step_ = true;

  if (position != data.size()) dserror("Mismatch in size of data %d <-> %d", data.size(), position);

  return;

}  // Unpack()


/*---------------------------------------------------------------------*
 | initialise / allocate internal stress variables (public) dano 08/11 |
 *---------------------------------------------------------------------*/
void MAT::ThermoPlasticLinElast::Setup(int numgp, DRT::INPUT::LineDefinition* linedef)
{
  // initialise history variables
  strainpllast_ = Teuchos::rcp(new std::vector<LINALG::Matrix<NUM_STRESS_3D, 1>>);
  strainplcurr_ = Teuchos::rcp(new std::vector<LINALG::Matrix<NUM_STRESS_3D, 1>>);

  backstresslast_ = Teuchos::rcp(new std::vector<LINALG::Matrix<NUM_STRESS_3D, 1>>);
  backstresscurr_ = Teuchos::rcp(new std::vector<LINALG::Matrix<NUM_STRESS_3D, 1>>);

  strainbarpllast_ = Teuchos::rcp(new std::vector<double>);
  strainbarplcurr_ = Teuchos::rcp(new std::vector<double>);

  Dmech_ = Teuchos::rcp(new std::vector<double>);
  Dmech_d_ = Teuchos::rcp(new std::vector<LINALG::Matrix<NUM_STRESS_3D, 1>>);

  Incstrainpl_ = Teuchos::rcp(new std::vector<LINALG::Matrix<NUM_STRESS_3D, 1>>);
  strainelrate_ = Teuchos::rcp(new std::vector<LINALG::Matrix<NUM_STRESS_3D, 1>>);

  LINALG::Matrix<NUM_STRESS_3D, 1> emptymat_vect(true);
  strainpllast_->resize(numgp);
  strainplcurr_->resize(numgp);

  backstresslast_->resize(numgp);
  backstresscurr_->resize(numgp);

  strainbarpllast_->resize(numgp);
  strainbarplcurr_->resize(numgp);

  Dmech_->resize(numgp);
  Dmech_d_->resize(numgp);

  Incstrainpl_->resize(numgp);
  strainelrate_->resize(numgp);

  for (int i = 0; i < numgp; i++)
  {
    strainpllast_->at(i) = emptymat_vect;
    strainplcurr_->at(i) = emptymat_vect;

    backstresslast_->at(i) = emptymat_vect;
    backstresscurr_->at(i) = emptymat_vect;

    strainbarpllast_->at(i) = 0.0;
    strainbarplcurr_->at(i) = 0.0;

    Dmech_->at(i) = 0.0;
    Dmech_d_->at(i) = emptymat_vect;

    Incstrainpl_->at(i) = emptymat_vect;
    strainelrate_->at(i) = emptymat_vect;
  }

  isinit_ = true;

  return;

}  // Setup()


/*---------------------------------------------------------------------*
 | update internal stress variables (public)                dano 08/11 |
 *---------------------------------------------------------------------*/
void MAT::ThermoPlasticLinElast::Update()
{
  // make current values at time step t_n+1 to values of last step t_n
  strainpllast_ = strainplcurr_;
  backstresslast_ = backstresscurr_;

  strainbarpllast_ = strainbarplcurr_;

  // empty vectors of current data
  strainplcurr_ = Teuchos::rcp(new std::vector<LINALG::Matrix<NUM_STRESS_3D, 1>>);
  backstresscurr_ = Teuchos::rcp(new std::vector<LINALG::Matrix<NUM_STRESS_3D, 1>>);

  strainbarplcurr_ = Teuchos::rcp(new std::vector<double>);

  // get the size of the vector
  // (use the last vector, because it includes latest results, current is empty)
  const int histsize = strainpllast_->size();
  strainplcurr_->resize(histsize);
  backstresscurr_->resize(histsize);

  strainbarplcurr_->resize(histsize);

  const LINALG::Matrix<NUM_STRESS_3D, 1> emptyvec(true);
  for (int i = 0; i < histsize; i++)
  {
    strainplcurr_->at(i) = emptyvec;
    backstresscurr_->at(i) = emptyvec;

    strainbarplcurr_->at(i) = 0.0;
  }

  return;
}  // Update()


/*----------------------------------------------------------------------*
 | evaluate material (public)                                dano 08/11 |
 *----------------------------------------------------------------------*/
void MAT::ThermoPlasticLinElast::Evaluate(const LINALG::Matrix<3, 3>* defgrd,
    const LINALG::Matrix<NUM_STRESS_3D, 1>* linstrain,
    Teuchos::ParameterList& params,            // parameter list for communication & HISTORY
    LINALG::Matrix<NUM_STRESS_3D, 1>* stress,  // 2nd PK-stress
    LINALG::Matrix<NUM_STRESS_3D, NUM_STRESS_3D>* cmat,  // material stiffness matrix
    const int eleGID)
{
  const int gp = params.get<int>("gp", -1);
  if (gp == -1) dserror("no Gauss point number provided in material");
  LINALG::Matrix<MAT::NUM_STRESS_3D, 1> plstrain(true);
  if (eleGID == -1) dserror("no element provided in material");

  // get material parameters
  // Young's modulus
  double young = params_->youngs_;
  // Poisson's ratio
  double nu = params_->poissonratio_;
  // initial yield stress
  double sigma_y0 = params_->yield_;
  // linear isotropic hardening modulus
  double Hiso = params_->isohard_;
  // linear kinematic hardening modulus
  double Hkin = params_->kinhard_;

  // initialise scalars
  // lame constant
  // shear modulus parameter mu == G
  double G = 0.0;
  G = young / (2.0 * (1.0 + nu));
  // bulk modulus kappa = E /( 3 ( 1 - 2 nu) )= lambda + 2/3 * mu
  double kappa = 0.0;
  kappa = young / (3.0 * (1.0 - 2.0 * nu));

  // build Cartesian identity 2-tensor I_{AB}
  LINALG::Matrix<6, 1> id2(true);
  for (int i = 0; i < 3; i++) id2(i) = 1.0;

  // linstrain (in): independent variable passed from the element
  //  strain^p: evolution is determined by the flow rule, history varible
  //  strain^e: definition of additive decomposition:
  //  strain^e = strain - strain^p
  // REMARK: stress-like 6-Voigt vector
  LINALG::Matrix<NUM_STRESS_3D, 1> strain(*linstrain);

  //---------------------------------------------------------------------------
  // elastic predictor (trial values)
  //---------------------------------------------------------------------------

  // ------------------------------------------------- old plastic strain
  // strain^{p,trial}_{n+1} = strain^p_n
  LINALG::Matrix<NUM_STRESS_3D, 1> strain_p(false);
  for (int i = 0; i < 6; i++) strain_p(i, 0) = strainpllast_->at(gp)(i, 0);

  // get old equivalent plastic strain only in case of plastic step
  double strainbar_p = 0.0;
  // accumulated or equivalent plastic strain (scalar-valued)
  // astrain^p,trial}_{n+1} = astrain^p_n
  strainbar_p = (strainbarpllast_->at(gp));
  if (strainbarpllast_->at(gp) < 0.0)
    dserror("accumulated plastic strain has to be equal to or greater than zero!");

  // ---------------------------------------------------- old back stress
  // beta^{trial}_{n+1} = beta_n
  LINALG::Matrix<NUM_STRESS_3D, 1> beta(false);
  for (int i = 0; i < 6; i++) beta(i, 0) = backstresslast_->at(gp)(i, 0);

  // --------------------------------------------------------- physical strains
  // convert engineering shear components into physical components
  // input strain is given in Voigt-notation

  // convert engineering shear component (in) into physical component
  for (int i = 3; i < 6; ++i) strain(i) /= 2.0;
  for (int i = 3; i < 6; ++i) strain_p(i) /= 2.0;

  // ----------------------------------------------- elastic trial strain
  // assume load step is elastic
  // strain^e_{n+1}
  LINALG::Matrix<NUM_STRESS_3D, 1> strain_e(true);

  // strain^{e,trial}_{n+1} = strain_{n+1} - strain^p_n
  LINALG::Matrix<NUM_STRESS_3D, 1> trialstrain_e(false);
  trialstrain_e.Update(1.0, strain, (-1.0), strain_p);

  // volumetric strain
  // trace of strain vector
  double tracestrain = trialstrain_e(0) + trialstrain_e(1) + trialstrain_e(2);
  // volstrain = 1/3 . tr( strain ) . Id
  LINALG::Matrix<NUM_STRESS_3D, 1> volumetricstrain(false);
  volumetricstrain.Update((tracestrain / 3.0), id2);

  // deviatoric strain
  // devstrain^e = strain^e - volstrain^e
  LINALG::Matrix<NUM_STRESS_3D, 1> devstrain(false);
  devstrain.Update(1.0, trialstrain_e, (-1.0), volumetricstrain);

  // ------------------------------------------------------- trial stress
  // pressure = kappa . tr( strain ): saved as scalar
  double p = kappa * tracestrain;

  // deviatoric stress = 2 . G . devstrain
  LINALG::Matrix<NUM_STRESS_3D, 1> devstress(false);
  devstress.Update((2.0 * G), devstrain);
  // be careful for shear stresses (sigma_12)
  // in Voigt-notation the shear strains have to be scaled with 1/2
  // normally done in the material tangent (cf. id4sharp)

  // ------------------------------------------ relative effective stress
  // eta^{trial}_{n+1} = s^{trial}_{n+1} - beta^{trial}_{n+1}
  LINALG::Matrix<NUM_STRESS_3D, 1> eta(true);
  RelDevStress(devstress, beta, eta);

  // J2 = 1/2 ( (eta11^{trial})^2 + (eta22^{trial})^2 + (eta33^{trial})^2
  //      + 2 . (eta12^{trial})^2 + 2 . (eta23^{trial})^2 + 2 . (eta13^{trial})^2)
  double J2 = 0.0;
  J2 = 1.0 / 2.0 * (eta(0) * eta(0) + eta(1) * eta(1) + eta(2) * eta(2)) + +eta(3) * eta(3) +
       eta(4) * eta(4) + eta(5) * eta(5);
  double etanorm = 0.0;
  etanorm = sqrt(eta(0) * eta(0) + eta(1) * eta(1) + eta(2) * eta(2) +
                 +2.0 * (eta(3) * eta(3) + eta(4) * eta(4) + eta(5) * eta(5)));

  // trial effective relative stress
  // qbar^{trial}_{n+1} := qbar(eta^{trial}_{n+1}) = \sqrt{ 3 . J2 }
  double qbar = 0.0;
  qbar = sqrt(3.0 * J2);

  // initialise the isotropic work hardening von Mises stress
  // sigma_yiso:= kappa = kappa(strainbar^p)
  double sigma_yiso = 0.0;

  //---------------------------------------------------------------------------
  // check plastic admissibility, Phi<=0 is admissble
  //---------------------------------------------------------------------------

  // ----------------------------------------------- trial yield function

  // calculate the yield stress
  // sigma_y = sigma_y0 + kappa(strainbar^p)
  // kappa == sigma_yiso, kappa is already used as material parameter
  double sigma_y = 0.0;

  bool bool_linisotrophard = false;
  // check if constant Hiso is given in input file
  if (((params_->sigma_y_).size()) == 0)
  {
    // so far: with linear isotropic hardening
    //         = sigma_y0 + Hiso . strainbar^p_{n}
    //         = sigma_y0 + Hiso . strainbar^{p, trial}_{n+1}
    sigma_yiso = Hiso * strainbar_p;
    sigma_y = sigma_y0 + sigma_yiso;

    bool_linisotrophard = true;
  }
  // calculate the isotropic hardening modulus and yield stress out of samples
  else
  {
    // calculate the isotropic hardening modulus with old plastic strains
    // Hiso = dsigma_y / d astrain^p
    Hiso = GetIsoHardAtStrainbarnp(strainbar_p);

    // calculate the uniaxial yield stress out of samples
    sigma_y = GetSigmaYAtStrainbarnp(strainbar_p);
  }

  // calculate the yield function with Dgamma = 0
  // Phi = \sqrt{ 3.0 . J2 } - sigma_y = q - sigma_y
  // with trial values: Phi_trial = q_trial - sigma_y
  double Phi_trial = 0.0;
  Phi_trial = qbar - sigma_y;

  // --------------------------------------------------------- initialise

  // if trial state is violated, i.e. it's a plastic load step, there are 2
  // possible states: plastic loading: heaviside = 1, elastic unloading = 0)
  double heaviside = 0.0;
  // incremental plastic multiplier Delta gamma
  double Dgamma = 0.0;

  // kinematic hardening curve of current time step and old time step
  // betabar = Hkin . strainbar_p
  // linear kinematic hardening: Hkin = const., else: Hkin = Hkin(strainnbar_p)
  double betabarold = 0.0;
  double betabar = 0.0;

  // unit flow vector Nbar (Prandtl-Reuss)
  // (using the updated relative stress eta_{n+1}, no longer eta_{n+1}^trial)
  // Nbar = ( eta^{trial}_{n+1} / || eta^{trial}_{n+1} || )
  LINALG::Matrix<NUM_STRESS_3D, 1> Nbar(true);

  // flow vector N (Prandtl-Reuss)
  // (using the updated relative stress eta_{n+1}, no longer eta^{trial}_{n+1})
  // N = sqrt{3/2} . ( eta^{trial}_{n+1} / || eta^{trial}_{n+1} || )
  //   = sqrt{3/2} . Nbar
  LINALG::Matrix<NUM_STRESS_3D, 1> N(true);

  //---------------------------------------------------------------------------
  // IF consistency condition is violated, i.e. plastic load step
  // ( Phi^{trial} > 0.0, Dgamma >= 0.0 )
  //---------------------------------------------------------------------------
  if (Phi_trial > 1.0e-08)  // if (Phi^{trial} > 0.0)
  {
    // only first plastic call is output at screen for every processor
    // visualisation of whole plastic behaviour via PLASTIC_STRAIN in postprocessing
    if (plastic_step_ == false)
    {
      plastic_eleID_ = eleGID;

      if ((plastic_step_ == false) and (eleGID == plastic_eleID_) and (gp == 0))
        std::cout << "plasticity starts in element = " << plastic_eleID_ << std::endl;

      plastic_step_ = true;
    }

    // calculate kinematic hardening stress of old time step
    // beta_{n} = Hkin . astrain^p_{n} = Hkin . strainbar^p, trial}_{n+1}
    betabarold = Hkin * strainbar_p;

    // --------------------------------------------------------- return-mapping

    // local Newton-Raphson

    // initialise
    const int itermax = 50;  // max. number of iterations
    int itnum = 0;           // iteration counter

    // Res:= residual of Newton iteration == yield function Phi
    double Res = 0.0;
    // calculate derivative of residual or tangent
    // ResTan = Phi' = d(Phi)/d(Dgamma)
    double ResTan = 0.0;
    // safety check: set to zero
    Dgamma = 0.0;
    betabar = 0.0;

    // start iteration with index m for local Newton
    while (true)
    {
      itnum++;
      // check for convergence

      // if not converged and (m > m_max)
      if (itnum > itermax)
      {
        dserror("local Newton iteration did not converge after iteration %3d/%3d with Res=%3f",
            itnum, itermax, Res);
      }
      // else: continue loop (m <= m_max)

      // Res := Phi = qbar^{trial}_{n+1}
      //             - Delta gamma (3 . G + Hkin(strainbar_p + Dgamma) )
      //             - sigma_y(strainbar_p + Dgamma)
      // - Hiso is introduced in residual via sigma_y(strainbar_p + Dgamma)
      // - Hkin: using the relation: Delta gamma . Hkin := betabar - betabarold
      //
      // --> Res = qbar^{trial} - 3 * G * Dgamma
      //           - betabar(strainbar_p + Dgamma) + betabarold
      //           - sigma_y(strainbar_p + Dgamma)
      Res = qbar - 3.0 * G * Dgamma - betabar + betabarold - sigma_y;

      // check for convergence
      double norm = abs(Res);
      // check: absolute value of Res has to be smaller than given tolerance
      if (norm < (params_->abstol_))
      {
#ifdef DEBUGMATERIAL
        if (gp == 0)
          printf(
              "Newton method converged after %i iterations; abs(Res)=  %-14.8E\n", itnum, abs(Res));
#endif  // #ifdef DEBUGMATERIAL
        break;
      }

      // plasticity with linear kinematic hardening
      // ResTan = -3G - Hkin(strainbar_p + Dgamma^{m-1}) - Hiso(strainbar_p + Dgamma^{m-1})
      // with Hiso = const. when considering LINEAR isotropic hardening
      ResTan = -3.0 * G - Hkin - Hiso;

      // incremental plastic multiplier Dgamma
      // Dgamma^{m} = Dgamma^{m-1} - Phi / Phi'
      Dgamma += (-Res) / ResTan;

      // -------------------------- local Newton update of plastic values

      // compute new residual of accumulatd plastic strains
      // astrain^p_{n+1} = astrain^p_n + Dgamma
      // astrain^p_{n+1} = SUM{Dgamma_n} from all time steps n
      // Kuhn-Tucker: Dgamma >= 0.0 --> astrain^p_{n+1} >= 0.0
      strainbar_p = strainbarpllast_->at(gp) + Dgamma;
      if (strainbar_p < 0.0)
        dserror("accumulated plastic strain has to be equal or greater than zero");

      // Prager's linear kinemativ hardening rule
      // kinematic hardening stress betabar (scalar-valued)
      // beta_{n+1} = Hkin * astrain^p_{n+1}
      betabar = Hkin * strainbar_p;

      if (bool_linisotrophard == true)
      {
        // linear isotropic hardening
        // sigma = sigma_y0 + sigma_yiso(strainbar^p_{n+1})
        sigma_yiso = Hiso * strainbar_p;
        sigma_y = sigma_y0 + sigma_yiso;
      }
      else  // constant_Hiso == false
      {
        // Hiso = dsigma_y / d astrain^p_{n+1}
        Hiso = GetIsoHardAtStrainbarnp(strainbar_p);
        // sigma_y = sigma_y(astrain^p_{n+1})
        sigma_y = GetSigmaYAtStrainbarnp(strainbar_p);
      }

#ifdef DEBUGMATERIAL
      if (gp == 0)
      {
        std::cout << "am 1.GP: local Newton: Res " << Res << std::endl;
        std::cout << "local Newton: ResTan " << ResTan << std::endl;
        std::cout << "local Newton: Dgamma " << Dgamma << std::endl;
        std::cout << "local Newton: betabarold " << betabarold << std::endl;
        std::cout << "local Newton: betabar " << betabar << "\n" << std::endl;
      }
#endif  // #ifdef DEBUGMATERIAL

    }  // end of local Newton iteration

    // --------------------------------------------------- plastic update

    // ---------------------------------------------- update flow vectors
    // unit flow vector Nbar = eta_{n+1}^{trial} / || eta_{n+1}^{trial} ||
    Nbar.Update(eta);
    Nbar.Scale(1.0 / etanorm);

    // flow vector N = sqrt(3/2) eta_{n+1}^{trial} / || eta_{n+1}^{trial} ||
    N.Update((sqrt(3.0 / 2.0)), Nbar);

    // update relative stress eta_{n+1}, cf. (7.193)
    // eta = ( 1 - (Delta gamma / qbar_{n+1}^{trial}) . [ 3 . G + Hkin] ) eta_{n+1}^{trial}
    // H_iso is not needed for update of the stress
    const double etafac = 1.0 - ((Dgamma / qbar) * (3.0 * G + Hkin));
    eta.Scale(etafac);

    // update back stress, cf. (7.197)
    // beta_{n+1} = beta_n . sqrt(2/3) . (betabar - betabarold) . eta / etanorm;
    // sqrt(2/3) N =  2/3 . ( sqrt(3/2) eta / etanorm)
    const double facbeta = 2.0 / 3.0 * (betabar - betabarold);
    beta.Update(facbeta, N, 1.0);

    // deviatoric stress
    // s = s_{n+1}^{trial} - 2 . G . Delta gamma . N
    const double facdevstress = (-2.0) * G * Dgamma;
    devstress.Update(facdevstress, N, 1.0);

    // total stress
    // sigma_{n+1} = s_{n+1} + p_{n+1} . id2
    // pressure/volumetric stress no influence due to plasticity
    Stress(p, devstress, *stress);

    // total strains
    // strain^e_{n+1} = strain^(e,trial)_{n+1} - Dgamma . N
    // compute converged engineering strain components (Voigt-notation)
    strain_e.Update(1.0, trialstrain_e, (-Dgamma), N);

    // strain^p_{n+1} = strain^p_n + Dgamma . N
    strain_p.Update(Dgamma, N, 1.0);

    // compute converged engineering strain components (Voigt-notation)
    for (int i = 3; i < 6; ++i) strain_e(i) *= 2.0;
    for (int i = 3; i < 6; ++i) strain_p(i) *= 2.0;

    // --------------------------------------------------- update history
    // plastic strain
    strainplcurr_->at(gp) = strain_p;

    // accumulated plastic strain
    strainbarplcurr_->at(gp) = strainbar_p;

    // back stress
    backstresscurr_->at(gp) = beta;

#ifdef DEBUGMATERIAL
    if (gp == 0)
    {
      std::cout << "LAST values\nplastic load: strainbarpllast_->at(gp) = "
                << strainbarpllast_->at(gp) << std::endl;
      std::cout << "plastic load: strainpllast->at(gp)\n " << strainpllast_->at(gp) << std::endl;
      std::cout << "plastic load: backstresslast_->at(gp)\n " << backstresslast_->at(gp)
                << std::endl;
      std::cout << "CURRENT values\n plastic load: strainbar_p = " << strainbar_p << std::endl;
      std::cout << "CURRENT values\nplastic load: strainbarplcurr_->at(gp) = "
                << strainbarplcurr_->at(gp) << std::endl;
      std::cout << "plastic load: strain_p\n " << strain_p << std::endl;
      std::cout << "plastic load: strainplcurr_->at(gp)\n " << strainplcurr_->at(gp) << std::endl;
      std::cout << "plastic load: backstresscurr_->at(gp)\n " << backstresscurr_->at(gp)
                << std::endl;
    }
#endif  // ifdef DEBUGMATERIAL

  }  // plastic corrector

  //---------------------------------------------------------------------------
  // ELSE: elastic step (Phi_trial <= 0.0, Dgamma = 0.0)
  //---------------------------------------------------------------------------
  else  // (Phi_trial <= 0.0)
  {
    // trial state vectors = result vectors of time step n+1
    // sigma^e_{n+1} = sigma^{e,trial}_{n+1} = s^{trial}_{n+1} + p . id2
    Stress(p, devstress, *stress);

    // total strains
    // strain^e_{n+1} = strain^(e,trial)_{n+1}
    // compute converged engineering strain components (Voigt-notation)
    strain_e.Update(trialstrain_e);
    for (int i = 3; i < 6; ++i) strain_e(i) *= 2.0;

    // no plastic yielding
    Dgamma = 0.0;
    // kinematic hardening curve of current time step and old time step
    // betabar = Hkin * strainbar_p
    // linear kinematic hardening: Hkin = const., else: Hkin = Hkin(strainnbar_p)
    betabarold = 0.0;
    betabar = 0.0;

    // pass the current plastic strains to the element (for visualisation)
    // compute converged engineering strain components (Voigt-notation)
    for (int i = 3; i < 6; ++i) strain_p(i) *= 2.0;

    // --------------------------------------------------------- update history
    // constant values for
    //  - plastic strains
    //  - accumulated plastic strains
    //  - back stress
    //    (--> relative stress)

    // as current history vectors are set to zero in Update(), the old values
    // need to be set instead, otherwise no constant plastic values are possible
    strainplcurr_->at(gp) = strainpllast_->at(gp);
    strainbarplcurr_->at(gp) = strainbarpllast_->at(gp);
    backstresscurr_->at(gp) = backstresslast_->at(gp);

#ifdef DEBUGMATERIAL
    if (gp == 0)
    {
      std::cout << "LAST values\nelastic load: strainbarpllast_->at(gp) = "
                << strainbarpllast_->at(gp) << std::endl;
      std::cout << "elastic load: strainpllast->at(gp)\n " << strainpllast_->at(gp) << std::endl;
      std::cout << "elastic load: backstresslast_->at(gp)\n " << backstresslast_->at(gp)
                << std::endl;
      std::cout << "CURRENT values\n elastic load: strainbar_p = " << strainbar_p << std::endl;
      std::cout << "CURRENT values\nelastic load: strainbarplcurr_->at(gp) = "
                << strainbarplcurr_->at(gp) << std::endl;
      std::cout << "elastic load: strain_p\n " << strain_p << std::endl;
      std::cout << "elastic load: strainplcurr_->at(gp)\n " << strainplcurr_->at(gp) << std::endl;
      std::cout << "elastic load: backstresscurr_->at(gp)\n " << backstresscurr_->at(gp)
                << std::endl;
    }
#endif  // DEBUGMATERIAL

  }  // elastic step

  //---------------------------------------------------------------------------
  // --------------------------------- consistent elastoplastic tangent modulus
  //---------------------------------------------------------------------------

  // using an associative flow rule: C_ep is symmetric
  // ( generally C_ep is nonsymmetric )
  // if Phi^trial = 0: two tangent stress-strain relations exist
  // plastic loading --> C == C_ep
  if (Dgamma > 0.0) heaviside = 1.0;
  // elastic unloading --> C == C_e
  else
    heaviside = 0.0;

  // using an associative flow rule: C_ep is symmetric
  // ( generally C_ep is nonsymmetric )
  SetupCmatElastoPlastic(*cmat, Dgamma, G, qbar, N, Nbar, heaviside, Hiso, Hkin);

#ifdef DEBUGMATERIAL
  std::cout << "Nach Setup Cep\n" << std::endl;
  std::cout << " Dgamma " << Dgamma << std::endl;
  std::cout << " G " << G << std::endl;
  std::cout << " qbar " << qbar << std::endl;
  std::cout << " flow vector " << N << std::endl;
  std::cout << " heaviside " << heaviside << std::endl;
  std::cout << " Kinematic hardening modul " << Hkin << std::endl;

  // build the elasto-plastic tangent modulus
  LINALG::Matrix<6, 6> cmatFD(true);

  std::cout << "cmat " << *cmat << std::endl;
#endif  // #ifdef DEBUGMATERIAL

  //---------------------------------------------------------------------------
  // ------------------------------------------ internal/mechanical dissipation
  //---------------------------------------------------------------------------

  // ----------------------------------- compute plastic strain increment
  // strain^p_{n+1}' = gamma' . N
  // with implicit Euler:
  // (strain^p_{n+1}-strain^p_n)/dt = Dgamma/dt . N
  // Inc_strain^p_{n+1} := (strain^p_{n+1}-strain^p_n)
  //                     = Dgamma . N = Dgamma . sqrt{3/2} eta_{n+1} / || eta_{n+1}||
  for (int i = 0; i < 6; i++) Incstrainpl_->at(gp)(i, 0) = Dgamma * N(i, 0);
  // --> plastic strain rate: strain^p_{n+1}' = Incstrainpl_/dt
  // scale with dt in StrainRateSplit()

  // ------------------------------------------------ dissipation for r_T
  // calculate mechanical dissipation required for thermo balance equation
  Dissipation(gp, sigma_yiso, Dgamma, N, *stress);

  // --------------------------------------- kinematic hardening for k_TT
  // temperature-dependent dissipated mechanical power
  // if (tr(strain^p) == 0) and (sigma_T(i,i)=const.) --> dot product of both is zero
  // safety check:
  double tracestrainp = 0.0;
  tracestrainp = strain_p(0) + strain_p(1) + strain_p(2);
  if (tracestrainp > 1.0E-8) dserror("trace of plastic strains is not equal to zero!");

  // ----------------------------------- linearisation of D_mech for k_Td
  DissipationCouplCond(*cmat, gp, G, Hiso, Hkin, heaviside, etanorm, Dgamma, N, *stress);

  // ----------------------------------------------------- postprocessing

  // plastic strain
  plstrain = strainplcurr_->at(gp);
  // save the plastic strain for postprocessing
  params.set<LINALG::Matrix<MAT::NUM_STRESS_3D, 1>>("plglstrain", plstrain);

  return;

}  // Evaluate()


/*----------------------------------------------------------------------*
 | computes linear stress tensor                             dano 05/11 |
 *----------------------------------------------------------------------*/
void MAT::ThermoPlasticLinElast::Stress(const double p,  // volumetric stress
    const LINALG::Matrix<NUM_STRESS_3D, 1>& devstress,   // deviatoric stress tensor
    LINALG::Matrix<NUM_STRESS_3D, 1>& stress             // 2nd PK-stress
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
void MAT::ThermoPlasticLinElast::RelDevStress(
    const LINALG::Matrix<NUM_STRESS_3D, 1>& devstress,  // deviatoric stress tensor
    const LINALG::Matrix<NUM_STRESS_3D, 1>& beta,       // back stress tensor
    LINALG::Matrix<NUM_STRESS_3D, 1>& eta               // relative stress
)
{
  // relative stress = deviatoric - back stress
  // eta = s - beta
  eta.Update(1.0, devstress, (-1.0), beta);

}  // RelDevStress()


/*----------------------------------------------------------------------*
 | computes isotropic elasticity tensor in matrix notion     dano 08/11 |
 | for 3d                                                               |
 *----------------------------------------------------------------------*/
void MAT::ThermoPlasticLinElast::SetupCmat(LINALG::Matrix<NUM_STRESS_3D, NUM_STRESS_3D>& cmat)
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
  const double mfac = young / ((1.0 + nu) * (1.0 - 2.0 * nu));  // factor

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
 | computes isotropic elasticity tensor in matrix notion     dano 05/11 |
 | for 3d                                                               |
 *----------------------------------------------------------------------*/
void MAT::ThermoPlasticLinElast::SetupCmatElastoPlastic(
    LINALG::Matrix<NUM_STRESS_3D, NUM_STRESS_3D>& cmat,  // elasto-plastic tangent modulus (out)
    double Dgamma,                                       // plastic multiplier
    double G,                                            // shear modulus
    double q,                                            // elastic trial von Mises effective stress
    LINALG::Matrix<NUM_STRESS_3D, 1> flowvector,         // flow vector
    LINALG::Matrix<NUM_STRESS_3D, 1> Nbar,               // unit flow vector
    double heaviside,                                    // Heaviside function
    double Hiso,                                         // isotropic hardening modulus
    double Hkin                                          // kinematic hardening modulus
)
{
  // incremental constitutive function for the stress tensor
  // sigma_{n+1} = [ cmat - (Dgamma 6 G^2/q) I_d ] : strain^{e,trial}_{n+1}
  // consistent tangent operator
  // D^{ep} := dsigma_{n+1} / dstrain^{e,trial}_{n+1}

  // depending on the flow vector Cmat_ep can be a fully-occupied matrix

  // C_ep = C_e - ( H^ . Dgamma . 6 . G^2 ) / qbar^{trial} . I_d +
  //        +  H^ . 6 . G^2 ( Dgamma/qbar^{trial} - 1/(3 G + Hkin + Hiso) ) Nbar \otimes Nbar

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

  // unit flow vector Nbar (cf. de Souza Neto (7.117)/(7.210) )
  // Nbar = eta^{trial}_{n+1} / || eta^{trial}_{n+1} || = sqrt(2/3) . N
  double flowfac = sqrt(2.0 / 3.0);
  flowvector.Scale(flowfac);

  // ------------------------------------------------------- elastic term
  // C_ep = C_e
  // add standard isotropic elasticity tensor C_e first
  SetupCmat(cmat);

  // ------------------------------------------------------ plastic terms

  // ------------------------------------------------- first plastic term
  // - ( H^ . Dgamma . 6 . G^2 ) / qbar^{trial} . I_d
  double epfac = 0.0;
  // elastic trial von Mises effective stress
  if (q != 0.0)
  {
    epfac = (-1.0) * heaviside * Dgamma * 6.0 * G * G / q;
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
  // +  H^ . 6 . G^2 ( Dgamma/qbar^{trial} - 1/(3 G + Hkin + Hiso) ) Nbar \otimes Nbar

  // unit flow vector (using co-linearity between trial and end state of eta)
  // Nbar = eta_{n+1} / || eta_{n+1} ||
  //      = eta^{trial}_{n+1} / || eta^{trial}_{n+1} ||

  // ------------------------------------------------------------ tangent

  if (q != 0.0)
  {
    double epfac2 = 0.0;
    // loop strains (columns)
    for (int k = 0; k < 6; ++k)
    {
      // loop stresses (rows)
      for (int i = 0; i < 6; ++i)
      {
        epfac2 = heaviside * 6.0 * G * G * (Dgamma / q - 1.0 / (3.0 * G + Hkin + Hiso));
        cmat(i, k) += epfac2 * Nbar(i) * Nbar(k);
      }  // end rows, loop i
    }    // end columns, loop k
  }      // (q != 0.0)

#ifdef DEBUGMATERIAL
  std::cout << "Ende SetupCmatElastPlast" << std::endl;
  std::cout << "Cep\n"
            << " Dgamma " << Dgamma << std::endl;
  std::cout << " G " << G << std::endl;
  std::cout << " q " << q << std::endl;
  std::cout << " flowvector " << flowvector << std::endl;
  std::cout << " heaviside " << heaviside << std::endl;
  std::cout << " epfac " << epfac << std::endl;
  std::cout << " epfac1 " << epfac1 << std::endl;
  std::cout << " epfac2 " << epfac2 << std::endl;
  std::cout << " cmat " << cmat << std::endl;
#endif  // #ifdef DEBUGMATERIAL

}  // SetupCmatElastoPlastic()


/*----------------------------------------------------------------------*
 | split given strain rate into elastic and plastic term     dano 08/11 |
 *----------------------------------------------------------------------*/
void MAT::ThermoPlasticLinElast::StrainRateSplit(int gp,  // current Gauss point
    const double stepsize,                                // step size
    LINALG::Matrix<NUM_STRESS_3D, 1>& strainrate          // total strain rate, i.e. B d'
)
{
  // elastic strain rate strain^e'
  // strain^e' = strain' - strain^p'
  for (int i = 0; i < NUM_STRESS_3D; i++)
  {
    // strain^e' = strain' - strain^p'

    // with strain^p' = Inc_strain^p / dt: use implicit Euler scheme
    strainelrate_->at(gp)(i) = strainrate(i, 0) - (1.0 / stepsize) * Incstrainpl_->at(gp)(i);
  }

  return;

}  // StrainRateSplit


/*----------------------------------------------------------------------*
 | compute internal dissipation term                         dano 04/13 |
 *----------------------------------------------------------------------*/
void MAT::ThermoPlasticLinElast::Dissipation(int gp,  // current Gauss point
    double sigma_yiso,                                // isotropic work hardening von Mises stress
    double Dgamma,                                    // plastic multiplier/increment
    LINALG::Matrix<NUM_STRESS_3D, 1> N,               // flow vector
    LINALG::Matrix<NUM_STRESS_3D, 1> stress           // total mechanical stress
)
{
  // D_mech = stress : strain^p' + beta : d(Phi)/d(beta) Dgamma
  //       = (stress - beta) : strain^p' = (stress - beta) : gamma' . N
  //       = eta : Dgamma/Dt . N

  // ----------------------------- dissipation due to kinematic hardening
  // (stress - beta) : strain^p_{n+1}'
  // with total stress: stress_d + stress_T

  // --------------------------------------- kinematic hardening for fint
  // stressdiff = stress_d_{n+1} - beta_{n+1} = s_{n+1} + p_{n+1} . I - beta_{n+1}
  LINALG::Matrix<NUM_STRESS_3D, 1> stressdiff(false);
  stressdiff.Update(1.0, stress, (-1.0), (backstresscurr_->at(gp)));

  // Dmech = (stress_d + sigma_T - beta) : Inc_strain^p_{n+1}
  double stressIncstrainpl = 0.0;
  stressIncstrainpl =
      Incstrainpl_->at(gp)(0) * stressdiff(0) + Incstrainpl_->at(gp)(1) * stressdiff(1) +
      Incstrainpl_->at(gp)(2) * stressdiff(2) + Incstrainpl_->at(gp)(3) * stressdiff(3) +
      Incstrainpl_->at(gp)(4) * stressdiff(4) + Incstrainpl_->at(gp)(5) * stressdiff(5);

  // --------------------------------------- isotropic hardening for fint
  // kappa(strainbar^p) . strainbar^p' = sigma_yiso . Dgamma/dt
  double isotropicdis = sigma_yiso * Dgamma;

  // return mechanical disspiation term due to mixed hardening hardening
  Dmech_->at(gp) = -stressIncstrainpl + isotropicdis;
  // time step not yet considered, i.e., Dmech_ is an energy, not a power
  // accumulated plastic strain

}  // Dissipation()


/*----------------------------------------------------------------------*
 | compute linearisation of internal dissipation for k_Td    dano 04/13 |
 *----------------------------------------------------------------------*/
void MAT::ThermoPlasticLinElast::DissipationCouplCond(
    LINALG::Matrix<NUM_STRESS_3D, NUM_STRESS_3D>& cmat,  // elasto-plastic tangent modulus (out)
    int gp,                                              // current Gauss point
    double G,                                            // shear modulus
    double Hiso,                                         // isotropic hardening modulus
    double Hkin,                                         // kinematic hardening modulus
    double heaviside,                                    // Heaviside function
    double etanorm,                                      // norm of eta^{trial}_{n+1}
    double Dgamma,                                       // plastic multiplier
    LINALG::Matrix<NUM_STRESS_3D, 1>& N,                 // flow vector
    LINALG::Matrix<NUM_STRESS_3D, 1>& stress             // flow vector
)
{
  // ----------------------------------- linearisation of D_mech for k_Td

  // calculate the derivation of the dissipation w.r.t. to the strains
  // (dD_mech/dstrain)
  // = N_T^T . (- dDmech_kin/ dstrain + dDmech_iso/ dstrain )
  //
  // = - N_T^T . (d [ (sigma_{d,T} - beta) . strain^p' ]/ dstrain)
  //   + N_T^T . (d [ kappa(strainbar^p) . strainbar^p' ]/ dstrain)

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
  // build Cartesian identity 2-tensor I_{AB}
  LINALG::Matrix<6, 1> id2(true);
  for (int i = 0; i < 3; i++) id2(i) = 1.0;

  // set Cartesian identity 4-tensor in 6-Voigt matrix notation
  // this is fully 'contra-variant' identity tensor, ie I^{ABCD}
  // REMARK: rows are stress-like 6-Voigt
  //         columns are stress-like 6-Voigt
  LINALG::Matrix<NUM_STRESS_3D, NUM_STRESS_3D> id4sharp(true);
  for (int i = 0; i < 3; i++) id4sharp(i, i) = 1.0;
  for (int i = 3; i < 6; i++) id4sharp(i, i) = 0.5;

  // ---------------------- linearisation of KINEMATIC hardening for k_Td

  // (dD_mech_kin/dstrain) = (d [ (sigma_{d,T} - beta) . strain^p' ]/ dstrain)
  //
  // d[ (sigma_{d,T} - beta) . strain^p' ]/ dstrain
  // = d(sigma_{d,T} - beta)/dstrain . strain^p'
  //   + (sigma_{d,T} - beta) . (dstrain^p')/ dstrain)
  //
  // sigma_T is independent of deformation, i.e. strains: dsigma_T/dstrain = 0
  //
  // = d(sigma_d - beta)/dstrain . strain^p'
  //   + (sigma_{d,T} - beta) . [(dstrain^p')/ dstrain]

  // d(sigma_d - beta)/dstrain = dstress_d/dstrain = C_ep
  // calculate C_ep . Inc_strain^p_{n+1}
  LINALG::Matrix<6, 1> cmatstrainpinc(false);
  cmatstrainpinc.Multiply(cmat, Incstrainpl_->at(gp));
  // --> divide by dt in thermo_ele

  // (sigma_d - beta) . [(dstrain^p')/ dstrain] = eta_{n+1} . [(dstrain^p')/ dstrain]

  // ---------------------- linearisation of plastic strain
  // [(dstrain^p')/ dstrain] = 1/Dt . [(dstrain^p_{n+1})/ dstrain_{n+1}^{e,trial}]
  // = 2G/(3 G + Hkin + Hiso) . N \otimes N
  //   + Dgamma . 2G / || eta^{trial}_{n+1} || [sqrt(3/2) I_d - N \otimes N]

  LINALG::Matrix<NUM_STRESS_3D, NUM_STRESS_3D> Dmech_kin_d(false);
  double fac_kinlin = 0.0;
  if (etanorm != 0.0)
  {
    fac_kinlin = heaviside * Dgamma * 2.0 * G / etanorm;
  }

  // I_d = id4sharp - 1/3 Id \otimes Id
  double fac_kinlin1 = 0.0;
  fac_kinlin1 = sqrt(3.0 / 2.0) * fac_kinlin;
  double fac_kinlin2 = 0.0;
  fac_kinlin2 = fac_kinlin1 / (-3.0);
  // contribution: Id4^#
  Dmech_kin_d.Update(fac_kinlin1, id4sharp);
  // contribution: Id \otimes Id
  Dmech_kin_d.MultiplyNT(fac_kinlin2, id2, id2, 1.0);

  double fac_lin_3 = 0.0;
  double fac_lin_4 = 0.0;
  fac_lin_4 = 3.0 * G * Hkin * Hiso;
  if (fac_lin_4 != 0) fac_lin_3 = heaviside * 2.0 * G / fac_lin_4;
  double fac_kinlin_flowvect = 0.0;
  fac_kinlin_flowvect = fac_lin_3 - fac_kinlin;

  // loop strains (columns)
  for (int k = 0; k < 6; ++k)
  {
    // loop stresses (rows)
    for (int i = 0; i < 6; ++i)
    {
      Dmech_kin_d(i, k) += fac_kinlin_flowvect * N(i) * N(k);
    }  // end rows, loop i
  }    // end columns, loop k

  // ---------------------- linearisation of ISOTROPIC hardening for k_Td

  // ----------------------------------------linearisation of Dmech_iso
  // dD_mech_iso/dstrain = (d [ kappa(strainbar^p) . strainbar^p' ]/ dstrain)
  //                     = Hiso . (d[ strainbar^p . strainbar^p' ]/dstrain)
  //
  // with linear isotropic hardening, i.e. sigma_yiso := kappa = Hiso . strainbar^p
  //
  // d[ sigma_yiso . strainbar^p' ]/dstrain^{trial}_{n+1}
  // = dkappa/dstrain^{trial}_{n+1} . strainbar^p' + kappa . dstrainbar^p'/dstrain^{trial}_{n+1}
  // = Hiso . dstrainbar^p/dstrain^{trial}_{n+1} . strainbar^p'
  //   + sigma_yiso . dstrainbar^p'/dstrain^{trial}_{n+1}
  //
  // dstrainbar^p/dstrain^{trial}_{n+1} = dDgamma/dstrain^{trial}_{n+1}
  //                                    = 2G/(3G + Hkin + Hiso) N_{n+1}
  // dstrainbar^p'/dstrain^{trial}_{n+1} = 1/Dt . 2G/(3G + Hkin + Hiso) N_{n+1}
  // --> as usual: calculation with dt is done in thermo_ele

  // dD_mech_iso/dstrain
  // = 2G/(3G + Hkin + Hiso) . (N_{n+1} . strainbar^p' + sigma_yiso . 1/Dt N_{n+1})
  // = 2G/(3G + Hkin + Hiso) . 1/Dt . (Inc_strainbar^p_{n+1} + sigma_yiso) N_{n+1}
  // = 2G/(3G + Hkin + Hiso) . 1/Dt . (strainbar^p_{n+1} - strainbar^p_n +
  //   + Hiso . strainbar^p_{n+1}) N_{n+1}
  // = 2G/(3G + Hkin + Hiso) . 1/Dt . (strainbar^p_{n+1} (1+Hiso) - strainbar^p_n) N_{n+1}
  double fac_liniso = 0.0;
  fac_liniso = fac_lin_3 * ((1.0 + Hiso) * strainbarplcurr_->at(gp) - -strainbarpllast_->at(gp));

  // ------------------------------------------------------ term for k_Td
  // add the linearisation term to D_mech_d
  LINALG::Matrix<NUM_STRESS_3D, 1> D_mech_d(false);
  D_mech_d.Multiply(Dmech_kin_d, stress);
  D_mech_d.Update((-1.0), cmatstrainpinc, (-1.0));
  D_mech_d.Update((fac_liniso), N, 1.0);
  // update history
  Dmech_d_->at(gp) = D_mech_d;

}  // DissipationCouplCond()


/*----------------------------------------------------------------------*
 | calculate stresses by evaluating the temperature tangent  dano 08/11 |
 *----------------------------------------------------------------------*/
void MAT::ThermoPlasticLinElast::Evaluate(
    const LINALG::Matrix<1, 1>& Ntemp,  // shapefcts . temperatures
    LINALG::Matrix<6, 1>& ctemp, LINALG::Matrix<6, 1>& stresstemp)
{
  SetupCthermo(ctemp);

  // calculate the temperature difference
  LINALG::Matrix<1, 1> init(false);
  init(0, 0) = (params_->thetainit_);
  // Delta T = T - T_0
  LINALG::Matrix<1, 1> deltaT(false);
  deltaT.Update(1.0, Ntemp, (-1.0), init);

  // temperature dependent stress
  // sigma = C_theta * Delta T = (m*I) * Delta T
  stresstemp.MultiplyNN(ctemp, deltaT);

  // if stresstemp(i,i)=const.: (sigma_T : strainp' == 0), because (tr(strainp') == 0)
  // for different thermal stresses, term has to be considered!!!
  //  // calculate temperature-dependent stress term for dissipation
  //  double tempstressIncstrainpl = 0.0;
  //  tempstressIncstrainpl = stresstemp(0) * Incstrainpl_->at(gp)(0)
  //                          + stresstemp(1) * Incstrainpl_->at(gp)(1)
  //                          + stresstemp(2) * Incstrainpl_->at(gp)(2)
  //                          + stresstemp(3) * Incstrainpl_->at(gp)(3)
  //                          + stresstemp(4) * Incstrainpl_->at(gp)(4)
  //                          + stresstemp(5) * Incstrainpl_->at(gp)(5);
  //  // term enters as negative value into the balance equation

}  // Evaluate


/*----------------------------------------------------------------------*
 | computes temperature dependent isotropic                  dano 05/10 |
 | elasticity tensor in matrix notion for 3d, second(!) order tensor    |
 *----------------------------------------------------------------------*/
void MAT::ThermoPlasticLinElast::SetupCthermo(LINALG::Matrix<NUM_STRESS_3D, 1>& ctemp)
{
  double m = STModulus();

  // isotropic elasticity tensor C_temp in Voigt matrix notation C_temp = m I
  //
  // Matrix-notation for 3D case
  //              [ m      0      0 ]
  //   C_temp =   [ 0      m      0 ]
  //              [ 0      0      m ]
  //
  //  in Vector notation
  //   C_temp =   [m, m, m, 0, 0, 0]^T
  //
  // write non-zero components

  // clear the material tangent
  ctemp.Clear();

  // loop over the element nodes
  for (int i = 0; i < 3; ++i) ctemp(i, 0) = m;  // non-zero entries only in main directions
  // remaining terms zero

}  // SetupCthermo()


/*----------------------------------------------------------------------*
 | calculates stress-temperature modulus                     dano 08/11 |
 *----------------------------------------------------------------------*/
double MAT::ThermoPlasticLinElast::STModulus()
{
  // initialise the parameters for the lame constants
  const double ym = params_->youngs_;
  const double pv = params_->poissonratio_;

  // initialise the thermal expansion coefficient
  const double thermexpans = params_->thermexpans_;

  // plane strain, rotational symmetry
  // E / (1+nu)
  const double c1 = ym / (1.0 + pv);
  // (E*nu) / ((1+nu)(1-2nu))
  const double b1 = c1 * pv / (1.0 - 2.0 * pv);

  // build the lame constants
  //            E
  //   mu = --------
  //        2*(1+nu)
  //                  E*nu
  //   lambda = ----------------
  //            (1+nu)*(1-2*nu)
  //
  //  \f \mu =  \frac{E}{2(1+\nu)} \f
  const double mu = 0.5 * c1;
  // lambda
  // \f \frac{E\,\nu}{(1-2\nu)(1+\nu)} \f
  const double lambda = b1;

  // stress-temperature modulus
  // \f m\, = \, -(2\,\cdot \mu \, +\, 3\cdot\lambda)\cdot\varalpha_T \f
  const double stmodulus = (-1.0) * (2.0 * mu + 3.0 * lambda) * thermexpans;

  return stmodulus;

}  // STModulus()


/*----------------------------------------------------------------------*
 | return derivative of piecewise linear function for the    dano 02/14 |
 | yield stress, i.e. isotropic hardening modulus at current            |
 | accumulated plastic strain                                           |
 *----------------------------------------------------------------------*/
double MAT::ThermoPlasticLinElast::GetIsoHardAtStrainbarnp(
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
 | compute current yield stress sigma_y(astrain^p)           dano 02/14 |
 | calculate yield stress from (sigma_y-astrain^p)-samples              |
 *----------------------------------------------------------------------*/
double MAT::ThermoPlasticLinElast::GetSigmaYAtStrainbarnp(
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
 | return names of visualization data (public)              dano 03/13 |
 *---------------------------------------------------------------------*/
void MAT::ThermoPlasticLinElast::VisNames(std::map<std::string, int>& names)
{
  std::string accumulatedstrain = "accumulatedstrain";
  names[accumulatedstrain] = 1;  // scalar
}  // VisNames()


/*---------------------------------------------------------------------*
 | return visualization data (public)                       dano 03/13 |
 *---------------------------------------------------------------------*/
bool MAT::ThermoPlasticLinElast::VisData(
    const std::string& name, std::vector<double>& data, int numgp, int eleID)
{
  if (name == "accumulatedstrain")
  {
    if ((int)data.size() != 1) dserror("size mismatch");
    double temp = 0.0;
    for (int iter = 0; iter < numgp; iter++) temp += AccumulatedStrain(iter);
    data[0] = temp / numgp;
  }
  return true;
}  // VisData()


/*----------------------------------------------------------------------*/
