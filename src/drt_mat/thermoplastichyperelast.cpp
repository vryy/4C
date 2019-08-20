/*----------------------------------------------------------------------*/
/*! \file
\brief Contains the functions to establish local material law /
       stress-strain law for an isotropic material following finite strain
       von-Mises plasticity with nonlinear isotropic hardening and general
       hyperelasticity (for the time being: NeoHooke).

       implementation is based on
       Simo and Miehe: "Associative coupled thermoplasticity at finite strains:
       Formulation, numerical analysis and implementation", in Computer Methods
       in Applied Mechanics and Engineering, 98:41–104, 1992.

       geometrically nonlinear, finite strains, rate-independent, thermo-plasticity

       example input line:
       [mm,ms,kg,K,GPa]
       MAT 1 MAT_Struct_ThrPlasticHyperElast YOUNG 206.9 NUE 0.29 DENS 7.8e-6
         CTE 1e-5 INITTEMP 293 YIELD 0.45 ISOHARD 0.12924 SATHARDENING 0.715
         HARDEXPO 16.93 YIELDSOFT 0.002 HARDSOFT 0.002 DISSFACT 0.9 TOL 1.0e-06

         Seitz 11/16: There are still some linearizations for the plastic
         heating terms missing, use with caution.

\level 3

\maintainer Matthias Mayr
*/
/*----------------------------------------------------------------------*
 | headers                                                              |
 *----------------------------------------------------------------------*/
#include "thermoplastichyperelast.H"
#include "matpar_bundle.H"
#include "material_service.H"

#include "../drt_lib/drt_globalproblem.H"


/*----------------------------------------------------------------------*
 | constructor (public)                                      dano 03/13 |
 *----------------------------------------------------------------------*/
MAT::PAR::ThermoPlasticHyperElast::ThermoPlasticHyperElast(Teuchos::RCP<MAT::PAR::Material> matdata)
    : Parameter(matdata),
      youngs_(matdata->GetDouble("YOUNG")),
      poissonratio_(matdata->GetDouble("NUE")),
      density_(matdata->GetDouble("DENS")),
      cte_(matdata->GetDouble("CTE")),
      inittemp_(matdata->GetDouble("INITTEMP")),
      yield_(matdata->GetDouble("YIELD")),
      isohard_(matdata->GetDouble("ISOHARD")),
      sathardening_(matdata->GetDouble("SATHARDENING")),
      hardexpo_(matdata->GetDouble("HARDEXPO")),
      yieldsoft_(matdata->GetDouble("YIELDSOFT")),
      hardsoft_(matdata->GetDouble("HARDSOFT")),
      abstol_(matdata->GetDouble("TOL"))
{
  if (sathardening_ < yield_)
    dserror("Saturation hardening must not be less than initial yield stress!");
  if (hardexpo_ < 0.0) dserror("Nonlinear hardening exponent must be non-negative!");
}


/*----------------------------------------------------------------------*
 | is called in Material::Factory from ReadMaterials()       dano 03/13 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<MAT::Material> MAT::PAR::ThermoPlasticHyperElast::CreateMaterial()
{
  return Teuchos::rcp(new MAT::ThermoPlasticHyperElast(this));
}

MAT::ThermoPlasticHyperElastType MAT::ThermoPlasticHyperElastType::instance_;


/*----------------------------------------------------------------------*
 | is called in Material::Factory from ReadMaterials()       dano 03/13 |
 *----------------------------------------------------------------------*/
DRT::ParObject* MAT::ThermoPlasticHyperElastType::Create(const std::vector<char>& data)
{
  MAT::ThermoPlasticHyperElast* thrplhyper = new MAT::ThermoPlasticHyperElast();
  thrplhyper->Unpack(data);
  return thrplhyper;
}


/*----------------------------------------------------------------------*
 | constructor (public)                                      dano 03/13 |
 *----------------------------------------------------------------------*/
MAT::ThermoPlasticHyperElast::ThermoPlasticHyperElast() : params_(NULL) {}


/*----------------------------------------------------------------------*
 | copy-constructor (public)                                 dano 03/13 |
 *----------------------------------------------------------------------*/
MAT::ThermoPlasticHyperElast::ThermoPlasticHyperElast(MAT::PAR::ThermoPlasticHyperElast* params)
    : params_(params), plastic_step_(false)
{
}


/*----------------------------------------------------------------------*
 | pack (public)                                             dano 03/13 |
 *----------------------------------------------------------------------*/
void MAT::ThermoPlasticHyperElast::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data, type);

  // matid
  int matid = -1;
  if (params_ != NULL) matid = params_->Id();  // in case we are in post-process mode
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
    histsize = defgrdlast_->size();
  }

  AddtoPack(data, histsize);  // length of history vector(s)
  for (int var = 0; var < histsize; ++var)
  {
    // insert history vectors to AddtoPack
    AddtoPack(data, defgrdlast_->at(var));
    AddtoPack(data, bebarlast_->at(var));
    AddtoPack(data, accplstrainlast_->at(var));

    // variables corresponding to temperture-dependency
    AddtoPack(data, mechdiss_->at(var));
    AddtoPack(data, mechdiss_kTT_->at(var));
    AddtoPack(data, mechdiss_kTd_->at(var));
    AddtoPack(data, Cmat_kdT_->at(var));
    AddtoPack(data, thrplheat_->at(var));
    AddtoPack(data, thrplheat_kTT_->at(var));
    AddtoPack(data, thrplheat_kTd_->at(var));
  }

  AddtoPack(data, plastic_step_);

  return;
}  // Pack()


/*----------------------------------------------------------------------*
 | unpack (public)                                           dano 03/13 |
 *----------------------------------------------------------------------*/
void MAT::ThermoPlasticHyperElast::Unpack(const std::vector<char>& data)
{
  isinit_ = true;
  std::vector<char>::size_type position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position, data, type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");

  // matid
  int matid;
  ExtractfromPack(position, data, matid);
  params_ = NULL;
  if (DRT::Problem::Instance()->Materials() != Teuchos::null)
  {
    if (DRT::Problem::Instance()->Materials()->Num() != 0)
    {
      const int probinst = DRT::Problem::Instance()->Materials()->GetReadFromProblem();
      MAT::PAR::Parameter* mat =
          DRT::Problem::Instance(probinst)->Materials()->ParameterById(matid);
      if (mat->Type() == MaterialType())
        params_ = static_cast<MAT::PAR::ThermoPlasticHyperElast*>(mat);
      else
        dserror("Type of parameter material %d does not fit to calling type %d", mat->Type(),
            MaterialType());
    }
  }

  // history data
  int histsize;
  ExtractfromPack(position, data, histsize);

  // if system is not yet initialised, the history vectors have to be intialized
  if (histsize == 0) isinit_ = false;

  defgrdlast_ = Teuchos::rcp(new std::vector<LINALG::Matrix<3, 3>>);
  defgrdcurr_ = Teuchos::rcp(new std::vector<LINALG::Matrix<3, 3>>);

  bebarlast_ = Teuchos::rcp(new std::vector<LINALG::Matrix<3, 3>>);
  bebarcurr_ = Teuchos::rcp(new std::vector<LINALG::Matrix<3, 3>>);

  accplstrainlast_ = Teuchos::rcp(new std::vector<double>);
  accplstraincurr_ = Teuchos::rcp(new std::vector<double>);

  mechdiss_ = Teuchos::rcp(new std::vector<double>);
  mechdiss_kTT_ = Teuchos::rcp(new std::vector<double>);
  mechdiss_kTd_ = Teuchos::rcp(new std::vector<LINALG::Matrix<NUM_STRESS_3D, 1>>);
  Cmat_kdT_ = Teuchos::rcp(new std::vector<LINALG::Matrix<NUM_STRESS_3D, 1>>);
  thrplheat_ = Teuchos::rcp(new std::vector<double>);
  thrplheat_kTT_ = Teuchos::rcp(new std::vector<double>);
  thrplheat_kTd_ = Teuchos::rcp(new std::vector<LINALG::Matrix<NUM_STRESS_3D, 1>>);

  for (int var = 0; var < histsize; ++var)
  {
    // initialise
    LINALG::Matrix<3, 3> tmp_matrix(true);
    LINALG::Matrix<NUM_STRESS_3D, 1> tmp_vect(true);
    double tmp_scalar = 0.0;

    ExtractfromPack(position, data, tmp_matrix);
    defgrdlast_->push_back(tmp_matrix);

    ExtractfromPack(position, data, tmp_matrix);
    bebarlast_->push_back(tmp_matrix);

    ExtractfromPack(position, data, tmp_scalar);
    accplstrainlast_->push_back(tmp_scalar);

    ExtractfromPack(position, data, tmp_scalar);
    mechdiss_->push_back(tmp_scalar);

    ExtractfromPack(position, data, tmp_scalar);
    mechdiss_kTT_->push_back(tmp_scalar);

    ExtractfromPack(position, data, tmp_vect);
    mechdiss_kTd_->push_back(tmp_vect);

    ExtractfromPack(position, data, tmp_vect);
    Cmat_kdT_->push_back(tmp_vect);

    ExtractfromPack(position, data, tmp_scalar);
    thrplheat_->push_back(tmp_scalar);

    ExtractfromPack(position, data, tmp_scalar);
    thrplheat_kTT_->push_back(tmp_scalar);

    ExtractfromPack(position, data, tmp_vect);
    thrplheat_kTd_->push_back(tmp_vect);

    // current vectors have to be initialised
    defgrdcurr_->push_back(tmp_matrix);
    bebarcurr_->push_back(tmp_matrix);
    accplstraincurr_->push_back(tmp_scalar);
  }

  plastic_step_ = false;
  int plastic_step;
  ExtractfromPack(position, data, plastic_step);

  // if it was already plastic before, set true
  if (plastic_step != 0) plastic_step_ = true;

  if (position != data.size()) dserror("Mismatch in size of data %d <-> %d", data.size(), position);

  return;

}  // Unpack


/*---------------------------------------------------------------------*
 | initialise / allocate internal variables (public)                   |
 *---------------------------------------------------------------------*/
void MAT::ThermoPlasticHyperElast::Setup(int numgp, DRT::INPUT::LineDefinition* linedef)
{
  // initialise hist variables
  defgrdlast_ = Teuchos::rcp(new std::vector<LINALG::Matrix<3, 3>>);
  defgrdcurr_ = Teuchos::rcp(new std::vector<LINALG::Matrix<3, 3>>);

  bebarlast_ = Teuchos::rcp(new std::vector<LINALG::Matrix<3, 3>>);
  bebarcurr_ = Teuchos::rcp(new std::vector<LINALG::Matrix<3, 3>>);

  accplstrainlast_ = Teuchos::rcp(new std::vector<double>);
  accplstraincurr_ = Teuchos::rcp(new std::vector<double>);

  mechdiss_ = Teuchos::rcp(new std::vector<double>);
  mechdiss_kTT_ = Teuchos::rcp(new std::vector<double>);
  mechdiss_kTd_ = Teuchos::rcp(new std::vector<LINALG::Matrix<6, 1>>);
  Cmat_kdT_ = Teuchos::rcp(new std::vector<LINALG::Matrix<6, 1>>);
  thrplheat_ = Teuchos::rcp(new std::vector<double>);
  thrplheat_kTT_ = Teuchos::rcp(new std::vector<double>);
  thrplheat_kTd_ = Teuchos::rcp(new std::vector<LINALG::Matrix<6, 1>>);

  defgrdlast_->resize(numgp);
  defgrdcurr_->resize(numgp);

  bebarlast_->resize(numgp);
  bebarcurr_->resize(numgp);

  accplstrainlast_->resize(numgp);
  accplstraincurr_->resize(numgp);

  mechdiss_->resize(numgp);
  mechdiss_kTT_->resize(numgp);
  mechdiss_kTd_->resize(numgp);
  Cmat_kdT_->resize(numgp);
  thrplheat_->resize(numgp);
  thrplheat_kTT_->resize(numgp);
  thrplheat_kTd_->resize(numgp);

  LINALG::Matrix<3, 3> emptymat(true);
  for (int i = 0; i < 3; i++) emptymat(i, i) = 1.0;
  LINALG::Matrix<6, 1> emptyvect(true);

  for (int i = 0; i < numgp; i++)
  {
    defgrdlast_->at(i) = emptymat;
    defgrdcurr_->at(i) = emptymat;

    bebarlast_->at(i) = emptymat;
    bebarcurr_->at(i) = emptymat;

    accplstrainlast_->at(i) = 0.0;
    accplstraincurr_->at(i) = 0.0;

    mechdiss_->at(i) = 0.0;
    mechdiss_kTT_->at(i) = 0.0;
    mechdiss_kTd_->at(i) = emptyvect;
    Cmat_kdT_->at(i) = emptyvect;
    thrplheat_->at(i) = 0.0;
    thrplheat_kTT_->at(i) = 0.0;
    thrplheat_kTd_->at(i) = emptyvect;
  }

  isinit_ = true;
  return;

}  // Setup()


/*----------------------------------------------------------------------*
 | update internal variables                                 dano 03/13 |
 *----------------------------------------------------------------------*/
void MAT::ThermoPlasticHyperElast::Update()
{
  // make current values at time step t_n+1 to values of last step t_n
  defgrdlast_ = defgrdcurr_;
  bebarlast_ = bebarcurr_;
  accplstrainlast_ = accplstraincurr_;

  // empty vectors of current data
  defgrdcurr_ = Teuchos::rcp(new std::vector<LINALG::Matrix<3, 3>>);
  bebarcurr_ = Teuchos::rcp(new std::vector<LINALG::Matrix<3, 3>>);
  accplstraincurr_ = Teuchos::rcp(new std::vector<double>);

  // get the size of the vector
  // (use the last vector, because it includes latest results, current is empty)
  const int histsize = defgrdlast_->size();
  defgrdcurr_->resize(histsize);
  bebarcurr_->resize(histsize);
  accplstraincurr_->resize(histsize);

  LINALG::Matrix<3, 3> emptymat(true);
  for (int i = 0; i < histsize; i++)
  {
    defgrdcurr_->at(i) = emptymat;
    bebarcurr_->at(i) = emptymat;
    accplstraincurr_->at(i) = 0.0;
  }

  return;
}  // Update()


/*----------------------------------------------------------------------*
 | calculate stress and constitutive tensor                  dano 03/13 |
 *----------------------------------------------------------------------*/
void MAT::ThermoPlasticHyperElast::Evaluate(const LINALG::Matrix<3, 3>* defgrd,
    const LINALG::Matrix<6, 1>* glstrain, Teuchos::ParameterList& params,
    LINALG::Matrix<6, 1>* stress, LINALG::Matrix<6, 6>* cmat, const int eleGID)
{
  // extract the gauss points from the parameter list
  const int gp = params.get<int>("gp", -1);
  if (gp == -1) dserror("no Gauss point number provided in material");
  if (eleGID == -1) dserror("no element provided in material");

  // elastic material data
  // -------------------------------------------------- get material parameters
  // Young's modulus
  const double ym = params_->youngs_;
  // Poisson's ratio
  const double nu = params_->poissonratio_;
  // shear modulus, mu=G
  const double G = ym / (2.0 * (1.0 + nu));
  // bulk modulus
  const double bulk = ym / (3.0 * (1.0 - 2.0 * nu));
  // linear isotropic hardening modulus
  double Hiso = params_->isohard_;
  // initial yield stress
  double sigma_y0 = params_->yield_;
  // yield stress softening
  double omega_0 = params_->yieldsoft_;
  // hardening softening
  double omega_h = params_->hardsoft_;
  // saturation hardening
  double sigma_y0infty = params_->sathardening_;
  // hardening exponent
  double hardexpo = params_->hardexpo_;
  // reference temperature
  double inittemp = params_->inittemp_;

  // 3x3 2nd-order identity matrix
  LINALG::Matrix<3, 3> id2(true);
  for (int i = 0; i < 3; i++) id2(i, i) = 1.0;

  // start with current deformation
  defgrdcurr_->at(gp) = *defgrd;
  // get the inverse F^{-1}
  LINALG::Matrix<3, 3> invdefgrdcurr(*defgrd);
  invdefgrdcurr.Invert();
  // calculate the Jacobi-determinant J = det(F_{n+1})
  double J = defgrd->Determinant();
  // determinant has to be >= 0
  // in case of too large dt, determinant is often negative --> reduce dt
  if (J < 0) dserror("Jacobi determinant is not allowed to be less than zero!");

  // ------------------------------------------------ multiplicative split of F
  // split deformation gradient in elastic and plastic part
  // F_{n+1} = F_{n+1}^e . F_{n+1}^p

  // relative deformation gradient
  // f_{n+1} = F_{n+1} . (F_n)^-1
  LINALG::Matrix<3, 3> defgrddelta(false);
  LINALG::Matrix<3, 3> invdefgrdlast(defgrdlast_->at(gp));
  invdefgrdlast.Invert();
  defgrddelta.Multiply(*defgrd, invdefgrdlast);

  // isochoric part of relative deformation gradient
  // fbar_{n+1} = Fbar_{n+1} . Fbar_n^{-1} = (J_{n+1}/J_n)^{-1/3}) . f_{n+1}
  // with J_{n+1}/J_n = det(fbar_)
  LINALG::Matrix<3, 3> defgrddeltabar(defgrddelta);
  defgrddeltabar.Scale(pow(defgrddelta.Determinant(), -1.0 / 3.0));

  // --------------------------------------------------------------------------
  // elastic predictor (trial values)
  // --------------------------------------------------------------------------

  // ----------------------------------------------------- elastic trial strain

  // assume load step is elastic
  // elastic left Cauchy-Green (LCG) trial state (isochoric) (9.3.13)
  // bbar_{n+1}^{e,trial} = Fbar_{n+1} (Cbar_{n}^{p-1}) . Fbar_{n+1}^T
  //                      = fbar_{n+1} (bbar_{n} . fbar_{n+1}^T
  // with history variable Cbar_{n+1}^{p-1})^{trial} = Cbar_{n}^{p-1}
  LINALG::Matrix<3, 3> bebar_trial(false);
  LINALG::Matrix<3, 3> tmp(false);
  tmp.Multiply(defgrddeltabar, bebarlast_->at(gp));
  bebar_trial.MultiplyNT(tmp, defgrddeltabar);
  // trace of strain vector
  double tracebebar = bebar_trial(0, 0) + bebar_trial(1, 1) + bebar_trial(2, 2);

  // ------------------------------------------------------- trial stress

  // trial Kirchhoff stress deviator (9.3.9)
  // s_{n+1)^{trial} = G . dev_bebar_{n+1}^{e,trial}
  // dev_bebar_trial = bebar_trial - volstrain^e
  //                 = bebar_trial - 1/3 . tr( bebar_trial ) . id2
  LINALG::Matrix<3, 3> devtau_trial(bebar_trial);
  for (int i = 0; i < 3; i++) devtau_trial(i, i) -= 1.0 / 3.0 * tracebebar;
  devtau_trial.Scale(G);

  // trial equivalent von Mises stress
  // q^{trial} = sqrt(s^{trial}_ij . s^{trial}_ij)
  double q_trial = 0.0;
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) q_trial += devtau_trial(i, j) * devtau_trial(i, j);
  q_trial = sqrt(q_trial);

  // -------------------------- extract scalar-valued element temperature
  // initialise temperature
  double scalartemp = 0.0;
  if (params.getEntryPtr("scalartemp") != NULL)
  {
    // TSI, i.e. temperature is available --> use this temperature
    scalartemp = params.get<double>("scalartemp", -1.0);
    if (scalartemp < 0.0) dserror("INadmissible value for the temperature: T=%3d", scalartemp);
  }
  // in case of purely structural analysis, i.e. isothermal: T = T_0, DeltaT = 0
  else
    scalartemp = inittemp;

  // ------------------------ temperature-dependent yield stress function

  // temperature-dependent isotropic hardening modulus
  double Hiso_temp = 0.0;
  Hiso_temp = Hiso * (1.0 - omega_h * (scalartemp - inittemp));

  // temperature-dependent yield stress
  double sigma_y0_temp = 0.0;
  sigma_y0_temp = sigma_y0 * (1.0 - omega_0 * (scalartemp - inittemp));

  // temperature-dependent saturation hardening
  double sigma_y0infty_temp = 0.0;
  sigma_y0infty_temp = sigma_y0infty * (1.0 - omega_h * (scalartemp - inittemp));

  // get old accumulated or equivalent plastic strain accplstrainlast
  double alpha = 0.0;
  alpha = accplstrainlast_->at(gp);
  if (accplstrainlast_->at(gp) < 0.0)
    dserror("accumulated plastic strain has to be equal to or greater than 0!");

  // thermodynamic force / hardening flux describing nonlinear isotropic hardening
  // sigma_iso = rho (d psi^_p) / (d accplstrain)
  // rho psi^_p = 1/2 . h(T) . alpha^2
  //              + [sigma_y0infty(T) - sigma_y(T)] . (alpha - (1- exp(-delta . alpha))
  // --> sigma_iso = h(T) . alpha + [sigma_y0infty(T) - sigma_y0(T)] . (1 - exp(-delta . alpha)
  double sigma_iso =
      Hiso_temp * alpha + (sigma_y0infty_temp - sigma_y0_temp) * (1.0 - exp(-hardexpo * alpha));

  // complete yield stress
  // sigma_y = sigma_y0_temp + sigma_iso
  double sigma_y = sigma_y0_temp + sigma_iso;

  // calculate yield function at trial state
  // Phi = || s_{n+1}^{trial} || - sqrt(2/3) . sigma_y
  double Phi_trial = 0.0;
  Phi_trial = q_trial - sqrt(2.0 / 3.0) * sigma_y;

  // stress variables tau = J_{n+1} . p_{n+1} . I + s_{n+1}
  LINALG::Matrix<3, 3> devtau(false);

  // some computations
  // mubar = 1/3 mu tr(bebar_{n+1}^{e,trial})
  double mubar = G * 1.0 / 3.0 * tracebebar;

  // initialise incremental plastic multiplier Delta gamma
  double Dgamma = 0.0;

  // unit spatial flow vector
  // n = s^{trial}_{n+1} / || s^{trial}_{n+1} || = s^{trial}_{n+1} / q^{trial}
  LINALG::Matrix<3, 3> n(devtau_trial);
  if (q_trial != 0.0) n.Scale(1.0 / q_trial);

  //-------------------------------------------------------------------
  // IF: elastic step (Phi_trial <= 0.0, Dgamma = 0.0)
  //-------------------------------------------------------------------
  if (Phi_trial <= 0.0)
  {
    // trial state vectors = result vectors of time step n+1

    // _n^{trial} --> _{n+1}
    bebarcurr_->at(gp) = bebar_trial;
    accplstraincurr_->at(gp) = accplstrainlast_->at(gp);
    devtau = devtau_trial;
    Dgamma = 0.0;

    // elastic load --> values are zero
    mechdiss_->at(gp) = 0.0;
    mechdiss_kTT_->at(gp) = 0.0;
    mechdiss_kTd_->at(gp).Scale(0.0);
    Cmat_kdT_->at(gp).Scale(0.0);
    thrplheat_->at(gp) = 0.0;
    thrplheat_kTT_->at(gp) = 0.0;
    thrplheat_kTd_->at(gp).Scale(0.0);
  }  // end if (Phi_trial <= 0.0), i.e. elastic step

  //-------------------------------------------------------------------
  // ELSE consistency condition is violated, i.e. plastic load step
  // (Phi_trial > 0.0, Dgamma >= 0.0)
  //-------------------------------------------------------------------
  else
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

    // ------------ local Newton Raphson to determine plastic multiplier Dgamma

    // initialise
    const int itermax = 50;  // max. number of iterations
    int itnum = 0;           // iteration counter

    // Res:= residual of Newton iteration == yield function Phi
    double Res = 0.0;
    // calculate derivative of residual or tangent
    // ResTan = Phi' = d(Phi)/d(Dgamma)
    double ResTan = 0.0;

    // start iteration with index m for local Newton
    while (true)
    {
      itnum++;
      // check for convergence
      if (itnum > itermax)
      {
        dserror("local Newton iteration did not converge after iteration %3d/%3d with Res=%3d",
            itnum, itermax, Res);
      }  // itnum > itermax
      // else: continue loop, i.e. m <= m_max

      // Res := Phi = qbar^{trial}_{n+1} - sqrt{2/3} (sigma_y0 + sigma_iso(alpha^{trial} )
      //      = qbar^{trial}_{n+1} - sqrt{2/3} (sigma_y0 + sigma_iso(alpha_n )
      Res = q_trial - 2.0 * mubar * Dgamma -
            sqrt(2.0 / 3.0) *
                (sigma_y0_temp + Hiso_temp * alpha +
                    (sigma_y0infty_temp - sigma_y0_temp) * (1.0 - exp(-hardexpo * alpha)));

      // check for convergence
      double norm = abs(Res);

      // check: absolute value of Res has to be smaller than given tolerance
      if (norm < (params_->abstol_))
      {
        break;
      }

      // tangent
      // ResTan = - 2 . mubar - sqrt(2/3) . [ H_iso . sqrt(2/3) +
      // (sigma_y0infty - sigma_y0) . (-exp(-delta . alpha)) . (-delta . sqrt(2/3))
      ResTan = -2.0 * mubar - (2.0 / 3.0) * Hiso_temp -
               2.0 / 3.0 * (sigma_y0infty_temp - sigma_y0_temp) * exp(-hardexpo * alpha) * hardexpo;

      // incremental plastic multiplier Dgamma
      // Dgamma^{m} = Dgamma^{m-1} - Phi / Phi'
      Dgamma += (-Res) / ResTan;
      // update accumulated plastic strain
      // alpha_{n+1} = alpha_n + sqrt{2/3} . Dgamma
      alpha = accplstrainlast_->at(gp) + (sqrt(2.0 / 3.0) * Dgamma);

      // q_trial and bebar (mubar) maintain constant during iteration

    }  // end of local Newton Raphson iteration

    // ------------------------------------------------ update Kirchhoff stress

    // update deviatoric Kirchhoff stress
    // s_{n+1} = s_{n+1}^{trial} - 2 . mubar . Delta gamma . n
    devtau.Update(devtau_trial);
    devtau.Update(-2.0 * mubar * Dgamma, n, 1.0);

    // --------------------------------------------------------- update history

    // update accumulated plastic strain
    accplstraincurr_->at(gp) = alpha;

    // update elastic LCG
    // bbar_{n+1}^e = bbar_{n+1}^{e,trial} - 2/3 . Dgamma . tr(bbar_{n+1}^{e,trial}) . n
    // see e.g. Simo, Comp.Inelaticity (9.3.7)
    // for nightly test case use present update procedure
    bebarcurr_->at(gp) = (bebar_trial);
    bebarcurr_->at(gp).Update((-2.0 / 3.0 * Dgamma * tracebebar), n, 1.0);

    // plastic step update
    if (Dgamma != 0.0)
    {
      // ----------------------------------------- preliminary calculations

      // ------------------------------------- calculate plastic directions

      // pull-back of the spatial flow vector n to N
      // N = F^{-1} . n . F^{-T}
      LINALG::Matrix<3, 3> N(false);
      tmp.Scale(0.0);  // reuse tmp, but reset first
      tmp.Multiply(invdefgrdcurr, n);
      N.MultiplyNT(tmp, invdefgrdcurr);

      // pull-back of the deviatoric part of n^2 to dev[N^2]
      // dev (n^2)
      LINALG::Matrix<3, 3> devnsquare(false);
      devnsquare.Multiply(n, n);
      double tracensquare = (devnsquare(0, 0) + devnsquare(1, 1) + devnsquare(2, 2));
      for (int i = 0; i < 3; i++) devnsquare(i, i) -= 1.0 / 3.0 * tracensquare;
      // dev (N^2) = F^{-1} . dev(n^2) . F^{-T}
      LINALG::Matrix<3, 3> devNsquare(false);
      tmp.Scale(0.0);  // reuse tmp, but reset first
      tmp.Multiply(invdefgrdcurr, devnsquare);
      devNsquare.MultiplyNT(tmp, invdefgrdcurr);

      // ------------------------------------------------------- linearisations
      // dsigma_y0_temp/dT_{n+1} = - omega_0 . sigma_y0 . N_T
      double dsigma_y0_temp_dT = sigma_y0 * (-omega_0);

      // dkappa_temp/dT_{n+1} = [- omega_h . Hiso . astrain^p +
      //                         + (- omega_h . sigma_y0infty + omega_0 . sigma_y0)
      //                            . (1 - exp(-delta . astrain^p)) ]  . N_T
      double dkappaT_dT =
          -omega_h * Hiso * alpha +
          (sigma_y0infty * (-omega_h) - sigma_y0 * (-omega_0)) * (1.0 - exp(-hardexpo * alpha));

      // dkappa_temp/dastrain_{n+1} = Hiso_temp +
      //      + (-delta) . (sigma_y0infty_temp - sigma_y0_temp) . [-exp(-delta . astrain)]
      double dkappa_dastrain = Hiso_temp + (sigma_y0infty_temp - sigma_y0_temp) *
                                               (-exp(-hardexpo * alpha)) * (-hardexpo);

      // beta_0 = 1 + 1/(3 mubar) . dkappa/dastrain^p
      double beta0 = 1.0 + 1.0 / (3.0 * mubar) * dkappa_dastrain;

      // dDgamma/dT_{n+1} = -sqrt(2/3) . dsigma_y/dT / [ 2 . mubar . beta0 ]
      double dDgamma_dT =
          -sqrt(2.0 / 3.0) * (dsigma_y0_temp_dT + dkappaT_dT) / (2.0 * mubar * beta0);

      // dkappa/dT = dkappa(T)/dT + dkappa(astrain)/dT . dastrain/dDgamma . dDgamma/dT
      //           = dkappaT_dT + dkappa_dastrain . \sqrt{2/3} . dDgamma_dT
      double dkappa_dT = dkappaT_dT + dkappa_dastrain * sqrt(2.0 / 3.0) * dDgamma_dT;

      // -------------------------------- calculate second derivatives of kappa
      // 2nd derivatives is required for K_TT/K_Td

      // d/dT(dkappa/dT) = d/dT(dkappa_T/dT) + d/dT(dkappa_dastrain . dastrain/dDgamma . dDgamma/dT)

      // purely thermal derivatives d/dT(dkappa_T/dT)
      // ddkappa_T_dTT = d/dT(dkappa_T/dT) = d/dT(0 + dkappaT_dastrain . \sqrt{2/3} . dDgamma_dT)
      //               = dT_dkappaT_dastrain . \sqrt{2/3} . dDgamma_dT)

      // implicit thermal derivatives d/dT(dkappa_dastrain . dastrain/dDgamma . dDgamma/dT)
      // d/dT(dkappa_dastrain . dastrain/dDgamma . dDgamma/dT)
      // = dT_dkappa_dastrain . sqrt(2/3) * dDgamma_dT + dkappa_dastrain . dastrain/dDgamma .
      // dT_dDgamma/dT

      // d/dT(dkappa/dT) = 2 . ddkappa_dastraindT . dastrain/dDgamma . dDgamma/dT
      double ddkappa_dTT = sqrt(2.0 / 3.0) * dDgamma_dT * 2.0 *
                           (Hiso * (-omega_h) + (sigma_y0infty * (-omega_h) + sigma_y0 * omega_0) *
                                                    (-exp(-hardexpo * alpha)) * (-hardexpo));
      // TODO in case of bad convergence: add dkappa_dastrain . dastrain/dDgamma . dT_dDgamma/dT

      double dkappa_dTdastrain = Hiso * (-omega_h) +
                                 (sigma_y0infty * (-omega_h) + sigma_y0 * omega_0) *
                                     (-exp(-hardexpo * alpha)) * (-hardexpo) +
                                 sqrt(2.0 / 3.0) * dDgamma_dT *
                                     (sigma_y0infty_temp - sigma_y0_temp) *
                                     (-exp(-hardexpo * alpha)) * hardexpo * hardexpo;

      // ------------------------------ calculate derivative of Dgamma w.r.t. E

      // spatial description:
      // 2 . dDgamma/dg = 1/beta0 . [ (1 - 2/3 || s || Dgamma / mubar ) . n
      //                          + || s || / mubar . dev[n^2] ]
      // material description:
      // dDgamma/dE = F^{-1} . (2 . dDgamma/dg) . F^{-T}
      //            = 1/beta0 . [ (1 - 2/3 . || s || . Dgamma / mubar ) . N
      //                          + || s || / mubar . dev[N^2] ]
      LINALG::Matrix<3, 3> dDgamma_dg(false);
      dDgamma_dg.Update((1.0 - 2.0 / 3.0 * q_trial * Dgamma / mubar), N);
      dDgamma_dg.Update((q_trial / mubar), devNsquare, 1.0);
      dDgamma_dg.Scale(1.0 / beta0);

      // -------------------------------------- internal/mechanical dissipation

      // --------------------------------------------------------- D_mech
      // D_mech = tau . [-1/2  Lie(b_e)] . b_e^{-1} - kappa . alpha'
      //        = sqrt(2/3) . sigma_y0(T_{n+1}) . Dgamma/Dt
      // be aware: multiplication with Dt is done in thermo element
      double mechdiss = sqrt(2.0 / 3.0) * Dgamma * sigma_y0_temp;
      mechdiss_->at(gp) = mechdiss;

      // -------------------------------------------- dD_mech/dT for k_TT

      // k_TT += ... dD_mech/dT_{n+1} ...
      // with dD_mech/dT_{n+1} = 1/Dt . sqrt(2/3) . [ dDgamma/dT_{n+1} . sigma_y0_temp
      //                         + Dgamma . dsigma_y0_temp/dT_{n+1} ]

      // with sigma_y0_temp = sigma_y0 . (1.0 - omega_0 . (scalartemp - inittemp) )
      // calculate the derivativ of sigma_y0(T_{n+1}) w.r.t. T_{n+1}
      // derivative of mechanical Dissipation w.r.t. temperatures
      mechdiss_kTT_->at(gp) =
          sqrt(2.0 / 3.0) * (dDgamma_dT * sigma_y0_temp + Dgamma * dsigma_y0_temp_dT);

      // -------------------------------------------- dD_mech/dd for k_Td
      // k_Td += dD_mech/dd_{n+1}
      //      += 1/Dt . sqrt(2/3) . [ sigma_y0_temp . dDgamma/dE ]
      LINALG::Matrix<3, 3> mechdiss_kTd_matrix(false);
      mechdiss_kTd_matrix.Update((sqrt(2.0 / 3.0) * sigma_y0_temp), dDgamma_dg);
      // Voigt notation
      LINALG::Matrix<6, 1> mechdiss_kTd_vct(false);
      mechdiss_kTd_vct(0) = mechdiss_kTd_matrix(0, 0);
      mechdiss_kTd_vct(1) = mechdiss_kTd_matrix(1, 1);
      mechdiss_kTd_vct(2) = mechdiss_kTd_matrix(2, 2);
      mechdiss_kTd_vct(3) = 0.5 * (mechdiss_kTd_matrix(0, 1) + mechdiss_kTd_matrix(1, 0));
      mechdiss_kTd_vct(4) = 0.5 * (mechdiss_kTd_matrix(1, 2) + mechdiss_kTd_matrix(2, 1));
      mechdiss_kTd_vct(5) = 0.5 * (mechdiss_kTd_matrix(0, 2) + mechdiss_kTd_matrix(2, 0));
      mechdiss_kTd_->at(gp).Update(mechdiss_kTd_vct);

      // ------------------------------------------- thermoplastic heating term

      // ------------------------------------------------------------ H_p

      // H_p = T . dkappa/dT . astrain'
      // with astrain' = sqrt(2/3) . Dgamma/Dt
      // thrplheat_ = dkappa/dT . sqrt(2/3) . Dgamma
      thrplheat_->at(gp) = dkappa_dT * sqrt(2.0 / 3.0) * Dgamma;
      // H_p = T . thrplheat_ . 1/Dt
      // be aware: multiplication with 1/Dt and (T = N_T . T) is done in thermo
      //           element

      // ----------------------------------------------- dH_p/dT for k_TT

      // k_TT += dH_p/dT_{n+1}
      //       = dkappa/dT . astrain^p'
      //         + 1/Dt . T . [ d/dT(dkappa/dT) . astrain^p'
      //                        + dkappa/dT . dastrain^p'/dT ]
      //
      //       = [ thrplheat_ . 1/Dt + T . thrplheat_kTT_ . 1/Dt ] . N_T
      // be aware: multiplication with 1/Dt and T is done in thermo element
      thrplheat_kTT_->at(gp) =
          ddkappa_dTT * sqrt(2.0 / 3.0) * Dgamma + dkappa_dT * sqrt(2.0 / 3.0) * dDgamma_dT;

      // ----------------------------------------------- dH_p/dd for k_Td
      // k_Td += ... dH_p/dE . dE/dd ...
      //
      // dH_p/dE = 1/Dt . [ ddkappa/(dT dE) . sqrt{2/3} . dDgamma + dkappa/dT. sqrt{2/3} .
      // dDgamma/dE
      //
      // with
      // ddkappa/(dT dE) = ddkappa/(dT dastrain) . dastrain/dDgamma . dDgamma/dE
      //                 = dkappa_dTdastrain . \sqrt(2/3) . dDgamma/dE
      //
      // dH_p/dE = [ dkappa_dTdastrain . (2/3) . dDgamma + dkappa/dT. sqrt{2/3} ]
      //           . 1/Dt . dDgamma/dE
      // be aware: multiplication with 1/Dt and dE/dd is done in thermo element
      double fac_thrpl_kTd = dkappa_dTdastrain * (2.0 / 3.0) * Dgamma + dkappa_dT * sqrt(2.0 / 3.0);
      LINALG::Matrix<3, 3> thrplheat_kTd_matrix(false);
      thrplheat_kTd_matrix.Update(fac_thrpl_kTd, dDgamma_dg);
      // Voigt notation
      LINALG::Matrix<6, 1> thrplheat_kTd_vct(false);
      thrplheat_kTd_vct(0) = thrplheat_kTd_matrix(0, 0);
      thrplheat_kTd_vct(1) = thrplheat_kTd_matrix(1, 1);
      thrplheat_kTd_vct(2) = thrplheat_kTd_matrix(2, 2);
      thrplheat_kTd_vct(3) = 0.5 * (thrplheat_kTd_matrix(0, 1) + thrplheat_kTd_matrix(1, 0));
      thrplheat_kTd_vct(4) = 0.5 * (thrplheat_kTd_matrix(1, 2) + thrplheat_kTd_matrix(2, 1));
      thrplheat_kTd_vct(5) = 0.5 * (thrplheat_kTd_matrix(0, 2) + thrplheat_kTd_matrix(2, 0));
      thrplheat_kTd_->at(gp).Update(thrplheat_kTd_vct);

      //------ linearisation of mechanical material tangent w.r.t. temperatures

      // ---------------------------------------------- dCmat_dT for k_dT

      // dCmat_dT += F^{-1} . ds_{n+1}/dT_{n+1} . F^{-T}
      // with ds_{n+1}/dT_{n+1} = - 2 . mubar . dDgamma/dT . n
      //                        = + sqrt(2/3) . dsigma_y_dT . 1/beta0 . N := beta5 . N
      LINALG::Matrix<3, 3> Cmat_kdT_matrix(false);
      Cmat_kdT_matrix.Update((-2.0 * mubar * dDgamma_dT), N);
      // Voigt notation
      LINALG::Matrix<6, 1> Cmat_kdT_vct(false);
      Cmat_kdT_vct(0) = Cmat_kdT_matrix(0, 0);
      Cmat_kdT_vct(1) = Cmat_kdT_matrix(1, 1);
      Cmat_kdT_vct(2) = Cmat_kdT_matrix(2, 2);
      Cmat_kdT_vct(3) = 0.5 * (Cmat_kdT_matrix(0, 1) + Cmat_kdT_matrix(1, 0));
      Cmat_kdT_vct(4) = 0.5 * (Cmat_kdT_matrix(1, 2) + Cmat_kdT_matrix(2, 1));
      Cmat_kdT_vct(5) = 0.5 * (Cmat_kdT_matrix(0, 2) + Cmat_kdT_matrix(2, 0));
      Cmat_kdT_->at(gp).Update(Cmat_kdT_vct);

    }  // (Dgamma != 0.0)

#ifdef DEBUGMATERIAL
    std::cout << "dsigma_y0_temp_dT = " << dsigma_y0_temp_dT << std::endl;
    std::cout << "dkappaT_dT = " << dkappaT_dT << std::endl;
    std::cout << "beta0 = " << beta0 << std::endl;
    std::cout << "- 2 * mubar * dDgamma_d = " << -2.0 * mubar * dDgamma_dT << std::endl;
    std::cout << "Cmat_kdT_vct = " << Cmat_kdT_vct << std::endl;
    std::cout << "mubar = " << mubar << std::endl;
    std::cout << "dDgamma_dT = " << dDgamma_dT << std::endl;
    std::cout << "N = " << N << std::endl;
    std::cout << "mechdiss_kTT_ = " << mechdiss_kTT_->at(gp) << std::endl;
    std::cout << "mechdiss_ = " << mechdiss_->at(gp) << std::endl;
    std::cout << "mechdiss_kTd_->at(gp) = " << mechdiss_kTd_->at(gp) << std::endl;
    std::cout << "Cmat_kdT_vct = " << *Cmat_kdT_vct << std::endl;
#endif  // DEBUGMATERIAL

  }  // end plastic step

  // ------------------------------------------------------ update final stress
  // add mean stress to gain Kirchhoff stress tau (9.2.6)
  // tau = J . p . I + devtau
  // with p := U'(J) = 1/2 . bulk . (J^2 -1 ) / J
  // --> tau = 1/2 . bulk ( (J^2 -1 ) . I + devtau
  // different to Miehe (2.37): p = bulk/2 (J^2 - 1) / J
  double p = bulk / 2.0 * (J * J - 1.0) / J;
  LINALG::Matrix<3, 3> tau(devtau);
  for (int i = 0; i < 3; i++) tau(i, i) += J * p;

  // transform Kirchhoff stress to 2.PK-stress
  // PK2 = F^{-1} . tau . F^{-T}
  LINALG::Matrix<3, 3> PK2;
  tmp.Scale(0.0);  // reuse tmp, but reset first
  tmp.Multiply(invdefgrdcurr, tau);
  PK2.MultiplyNT(tmp, invdefgrdcurr);

  // output PK2-stress in Voigt-notation
  (*stress)(0) = PK2(0, 0);
  (*stress)(1) = PK2(1, 1);
  (*stress)(2) = PK2(2, 2);
  (*stress)(3) = 0.5 * (PK2(0, 1) + PK2(1, 0));
  (*stress)(4) = 0.5 * (PK2(1, 2) + PK2(2, 1));
  (*stress)(5) = 0.5 * (PK2(0, 2) + PK2(2, 0));

  // ----------------------- consistent elastoplastic tangent modulus (Box 9.2)
  SetupCmatElastoPlastic(*cmat,  // (o) elasto-plastic tangent modulus
      Dgamma,                    // plastic multiplier
      Hiso_temp,                 // H(T)
      sigma_y0infty_temp,        // trial value of saturation yield stress
      sigma_y0_temp,             // trial value of initial yield stress
      mubar, q_trial,
      *defgrd,        // F
      invdefgrdcurr,  // F^{-1}
      n,              // spatial flow vector
      bulk,           // isotropic thermodynamic force
      gp              // current Gauss point
  );

  return;

}  // Evaluate()


/*----------------------------------------------------------------------*
 | Calculation of consistent elastoplastic tangent modulus   dano 09/13 |
 *----------------------------------------------------------------------*/
void MAT::ThermoPlasticHyperElast::SetupCmatElastoPlastic(
    LINALG::Matrix<NUM_STRESS_3D, NUM_STRESS_3D>& cmat,  // elasto-plastic tangent modulus (out)
    double Dgamma, double Hiso_temp, double sigma_y0infty_temp, double sigma_y0_temp, double mubar,
    double q_trial,                      // || s_{n+1}^{trial} ||
    const LINALG::Matrix<3, 3>& defgrd,  // F
    LINALG::Matrix<3, 3> invdefgrdcurr, LINALG::Matrix<3, 3> n,
    double bulk,  // bulk modulus
    int gp        // current Gauss point
)
{
  // ---------------------------------------------- intialise material tangents
  LINALG::Matrix<6, 6> Cmat(true);
  LINALG::Matrix<6, 6> Cbar_trialMaterial(true);

  // Cmat = C_ep = C_e + Cbar_trial + Cbar_p

  // Cbar_ep = (1 - beta1) . Cbar_{n+1}^{trial} - 2 mubar beta3 n \otimes n
  //           - 2 mubar beta4 sym[ n \otimes dev[n^2] ]^s
  // with Cbar_{n+1}^{trial} = 2 mubar I_d
  //              - 2/3 (I \otimes s_{n+1}^{trial} + s_{n+1}^{trial} \otimes I)

  // ---------------------------------------------------------- calculate terms
  // initialise some variables
  LINALG::Matrix<3, 3> tmp(false);
  // hardening exponent
  double hardexpo = params_->hardexpo_;
  // determinant of the deformation gradient
  double J = defgrd.Determinant();

  // calculate the right Cauchy Green (RCG) deformation tensor and its inverse
  LINALG::Matrix<3, 3> RCG(false);
  RCG.MultiplyTN(defgrd, defgrd);
  LINALG::Matrix<3, 3> invRCG;
  invRCG.Invert(RCG);

  // --------------------------------------- calculate plastic directions
  // pull-back of the spatial flow vector n to N
  // N = F^{-1} n F^{-T}
  LINALG::Matrix<3, 3> N(false);
  tmp.Multiply(invdefgrdcurr, n);
  N.MultiplyNT(tmp, invdefgrdcurr);

  // dev (n^2)
  LINALG::Matrix<3, 3> devnsquare(false);
  devnsquare.Multiply(n, n);
  double tracensquare = (devnsquare(0, 0) + devnsquare(1, 1) + devnsquare(2, 2));
  for (int i = 0; i < 3; i++) devnsquare(i, i) -= 1.0 / 3.0 * tracensquare;

  // pull-back of dev(n^2)
  LINALG::Matrix<3, 3> devNsquare(false);
  tmp.Scale(0.0);  // reuse tmp, but reset first
  tmp.Multiply(invdefgrdcurr, devnsquare);
  devNsquare.MultiplyNT(tmp, invdefgrdcurr);

  // ----------------------------------------------------------- calculate Cmat
  // ------------------------------------------------ isochoric part Cbar
  // Cbar
  // spatial: cbar_trial = 2 . mubar . I_d - 2/3 qbar [n \otimes I + I \otimes n]
  // with I_d = I_s - 1/3 . I \otimes I
  // pull-back of I --> invRCG
  // Cbar += Cbar_trial = 2 . mubar . pullback_I_d
  ElastSymTensor_o_Multiply(Cbar_trialMaterial, 2.0 * mubar, invRCG, invRCG, 1.0);
  ElastSymTensorMultiply(Cbar_trialMaterial, -2.0 / 3.0 * mubar, invRCG, invRCG, 1.0);
  // Cbar += - 2/3 qbar [N \otimes C^{-1} + C^{-1} \otimes N]
  ElastSymTensorMultiplyAddSym(Cbar_trialMaterial, -2.0 / 3.0 * q_trial, N, invRCG, 1.0);

  // ------------------------------------------------ volumetric part C_e
  // spatial c_e = (J . U')' . J . I \otimes I - 2 J U' I4
  // with U'(J) = bulk/2 . (J^2 -1)  / J
  // C_e = bulk . J^2 [C^{-1} \otimes C^{-1}] - bulk ( J^2 -1 ) [C^{-1} \otimes C^{-1}]
  // with - bulk ( J^2 -1 ) [C^{-1} \otimes C^{-1}] = - bulk ( J^2 -1 ) [C^{-1} boeppel C^{-1}]
  ElastSymTensorMultiply(Cmat, bulk * J * J, invRCG, invRCG, 1.0);
  ElastSymTensor_o_Multiply(Cmat, -1.0 * bulk * (J * J - 1.0), invRCG, invRCG, 1.0);
  Cmat.Update(1.0, Cbar_trialMaterial, 1.0);

  // plastic step update
  if (Dgamma != 0.0)
  {
    // ------------------------------------ scaling factors for spatial tangent

    // beta_0 = 1 + 1/(3 mubar) . dkappa_{n+1}/dastrain_{n+1}
    // with dkappa_{n+1}/dastrain_{n+1} = Hiso_temp + (sigma_y0infty_temp - sigma_y0_temp)
    //                                    . [-exp(-delta . astrain)] . (-delta)
    double beta0 = 0.0;
    beta0 = 1.0 + (Hiso_temp + (sigma_y0infty_temp - sigma_y0_temp) *
                                   exp(-hardexpo * accplstraincurr_->at(gp)) * hardexpo) /
                      (3.0 * mubar);

    // beta_1 = 2 . mubar . Dgamma / || s_{n+1}^{trial} ||
    double beta1 = 0.0;
    beta1 = 2.0 * mubar * Dgamma / q_trial;

    // beta_2 = (1 - 1/beta_0) . 2/3 . Dgamma / mubar . || s_{n+1}^{trial} ||
    double beta2 = 0.0;
    beta2 = (1.0 - 1.0 / beta0) * 2.0 / 3.0 * Dgamma / mubar * q_trial;

    // beta_3 = 1/beta_0 - beta_1 + beta_2
    double beta3 = 0.0;
    beta3 = 1.0 / beta0 - beta1 + beta2;

    // beta_4 = (1/beta_0 - beta_1) . || s_{n+1}^{trial} || / mubar
    double beta4 = 0.0;
    beta4 = (1.0 / beta0 - beta1) * q_trial / mubar;

    // this is nonlinear mechanics
    Cmat.Update((-1.0 * beta1), Cbar_trialMaterial, 1.0);
    ElastSymTensorMultiply(Cmat, (-2.0 * mubar * beta3), N, N, 1.0);
    ElastSymTensorMultiply(Cmat, (-2.0 * mubar * beta4), N, devNsquare, 1.0);
  }  // Dgamma != 0

  // update material tangent
  // cmat = C_ep = C_e + Cbar_trial + Cbar_p
  cmat = Cmat;

}  // SetupCmatElastoPlastic()


/*----------------------------------------------------------------------*
 | calculate final isochoric elastic LCG bbar^e_{n+1}        dano 09/13 |
 *----------------------------------------------------------------------*/
void MAT::ThermoPlasticHyperElast::CalculateCurrentBebar(
    const LINALG::Matrix<3, 3>& devtau,  // s_{n+1}
    double G,                            // shear modulus
    const LINALG::Matrix<3, 3>& id2,     // second-order identity
    int gp                               // current Gauss-point
)
{
  // calculate equivalent von Mises stress || s_{n+1} ||
  double q = 0.0;
  for (int i = 0; i < 3; i++)
  {
    for (int j = 0; j < 3; j++)
    {
      q += devtau(i, j) * devtau(i, j);
    }
  }

  // ------------------------------------------- calculate isochoric invariants
  // calculate isochoric second invariant of dev[bebar] J2_bar
  double J2_bar = 0.0;
  J2_bar = 1.0 / 2.0 * q * q / (G * G);

  // calculate isochoric third invariant of dev[bebar] (i.e. determinant) J3_bar
  double J3_bar = 0.0;
  J3_bar = 1.0 / (G * G * G) * devtau.Determinant();

  // --------------------------------------------------- calculate coefficients
  // coefficient q
  // q = (1 - J3_bar) / 2
  double q_coeff = 0.0;
  q_coeff = 1.0 / 2.0 * (1.0 - J3_bar);

  // coefficient d
  // d = - J2_bar^3 / 27 + q^2
  double d_coeff = 0.0;
  d_coeff = -std::pow(J2_bar, 3.0) / 27.0 + std::pow(q_coeff, 2.0);

  // -------------------------------------------------- calculate updated bebar

  // calculate scaling factor for updated bebar
  // 1/3 tr(bebar) = 1/3 . Ibar_1
  double const third_Ibar_1 = std::pow((q_coeff + sqrt(d_coeff)), 1.0 / 3.0) +
                              std::pow((q_coeff - sqrt(d_coeff)), 1.0 / 3.0);

  // bebar_{n+1} = s_{n+1}/mu + 1/3 Ibar_1 . id2
  bebarcurr_->at(gp) = (devtau);
  bebarcurr_->at(gp).Scale(1 / G);
  bebarcurr_->at(gp).Update(third_Ibar_1, id2, 1.0);

}  // CalculateCurrentBebar()


/*----------------------------------------------------------------------*
 | calculate temperature-dependent stresses                  dano 09/13 |
 | is called from so3thermo element                                     |
 *----------------------------------------------------------------------*/
void MAT::ThermoPlasticHyperElast::Evaluate(const LINALG::Matrix<1, 1>& Ntemp,
    LINALG::Matrix<6, 1>& ctemp, LINALG::Matrix<6, 6>& cmat_T, LINALG::Matrix<6, 1>& stresstemp,
    Teuchos::ParameterList& params)
{
  // calculate the temperature difference
  LINALG::Matrix<1, 1> init(false);
  init(0, 0) = params_->inittemp_;
  // Delta T = T - T_0
  LINALG::Matrix<1, 1> deltaT(false);
  deltaT.Update(1.0, Ntemp, (-1.0), init);

  // get the temperature-dependent material tangent
  SetupCthermo(ctemp, params);

  // get the temperature-dependent mechanical material tangent
  SetupCmatThermo(Ntemp, cmat_T, params);

  // calculate thermal stresses
  // tau = ctemp_AK . Delta T = m_0/2.0 . (J + 1/J) . I . Delta T
  // pull-back of Kirchhoff-stresses to PK2-stresses
  // PK2 = F^{-1} . tau . F^{-T}
  // --> PK2 = ctemp . Delta T = m_0/2.0 . (J + 1/J). Cinv . Delta T
  stresstemp.MultiplyNN(ctemp, deltaT);

#ifdef DEBUGMATERIAL
  // ------------- FDcheck of temperature-dependent mechanical material tangent

  // in case we want to test the material tangent without Delta T in the FD Check
  //  stresstemp.Update(ctemp);

  // build the elasto-plastic tangent modulus
  LINALG::Matrix<NUM_STRESS_3D, NUM_STRESS_3D> cmat_TFD(true);
  FDCheck(stresstemp, cmat_T, cmat_TFD, Ntemp, params);
  std::cout << "cmat_T " << cmat_T << std::endl;
  std::cout << "cmat_TFD " << cmat_TFD << std::endl;

  std::cout << "Evaluate Material: Ntemp = " << Ntemp << std::endl;
  std::cout << "Evaluate Material: deltaT = " << deltaT << std::endl;
  std::cout << "Evaluate Material: ctemp\n" << ctemp << std::endl;
  std::cout << "Evaluate Material: cmat_T\n" << cmat_T << std::endl;
  std::cout << "Evaluate Material: thermal stress stresstemp\n" << stresstemp << std::endl;
#endif  // DEBUGMATERIAL

}  // THREvaluate()

/*----------------------------------------------------------------------*
 | computes temperature-dependent isotropic                  dano 09/13 |
 | elasticity tensor in matrix notion for 3d, second(!) order tensor    |
 *----------------------------------------------------------------------*/
void MAT::ThermoPlasticHyperElast::SetupCmatThermo(
    const LINALG::Matrix<1, 1>& Ntemp, LINALG::Matrix<6, 6>& cmat_T, Teuchos::ParameterList& params)
{
  // temperature-dependent material tangent
  // cmat_T = cmat_vol,dT = dstresstemp/dE = 2 dstresstemp/dC
  //        = (T - T_0) . m_0/2.0 . (J - 1/J) (C^{-1} \otimes C^{-1}) -
  //          - (T - T_0) . m_0 . (J + 1/J) ( Cinv boeppel Cinv )

  // calculate the temperature difference
  // Delta T = T - T_0
  const double deltaT = Ntemp(0, 0) - (params_->inittemp_);

  // extract F and Cinv from params
  LINALG::Matrix<3, 3> defgrd = params.get<LINALG::Matrix<3, 3>>("defgrd");

  // get stress-temperature modulus
  double m_0 = STModulus();
  // get Jacobi
  double J = defgrd.Determinant();
  // calculate the right Cauchy Green (RCG) deformation tensor and its inverse
  LINALG::Matrix<3, 3> RCG(false);
  RCG.MultiplyTN(defgrd, defgrd);
  LINALG::Matrix<3, 3> invRCG;
  invRCG.Invert(RCG);

  // clear the material tangent
  cmat_T.Clear();

  // cmat_T = 2 . dS_vol,dT/dd
  //        = (T - T_0) . m_0/2 . (J - 1/J) (C^{-1} \otimes C^{-1})
  //          - (T - T_0) . m_0 . (J + 1/J) ( Cinv boeppel Cinv )
  ElastSymTensorMultiply(cmat_T, (deltaT * m_0 / 2.0 * (J - 1 / J)), invRCG, invRCG, 1.0);
  ElastSymTensor_o_Multiply(cmat_T, (-deltaT * m_0 * (J + 1 / J)), invRCG, invRCG, 1.0);

#ifdef DEBUGMATERIAL
  std::cout << "SetupCmatThermo(): Jacobi determinant J = " << J << std::endl;
  std::cout << "SetupCmatThermo(): 1.0 * (deltaT * m_0/2 * (J + 1/J)) = "
            << 1.0 * (deltaT * m_0 / 2.0 * (J + 1 / J)) << std::endl;
  std::cout << "SetupCmatThermo(): deltaT = " << deltaT << std::endl;
  std::cout << "SetupCmatThermo(): Ntemp = " << Ntemp << std::endl;
  std::cout << "SetupCmatThermo(): inittemp = " << inittemp << std::endl;
#endif  // DEBUGMATERIAL

}  // SetupCmatThermo()


/*----------------------------------------------------------------------*
 | computes temperature-dependent isotropic                  dano 09/13 |
 | elasticity tensor in matrix notion for 3d, second(!) order tensor    |
 *----------------------------------------------------------------------*/
void MAT::ThermoPlasticHyperElast::SetupCthermo(
    LINALG::Matrix<6, 1>& ctemp, Teuchos::ParameterList& params)
{
  // temperature-dependent material tangent
  // C_T = m_0/2.0 . (J + 1/J) . Cinv

  // extract F and Cinv from params
  LINALG::Matrix<3, 3> defgrd = params.get<LINALG::Matrix<3, 3>>("defgrd");
  LINALG::Matrix<6, 1> Cinv = params.get<LINALG::Matrix<6, 1>>("Cinv_vct");

  // temperature-dependent stress temperature modulus
  // m = m(J) = m_0 .(J+1)/J = m_0 . (J + 1/J)
  const double m_0 = STModulus();
  const double J = defgrd.Determinant();
  const double m = m_0 * (J + 1.0 / J);

  // clear the material tangent
  ctemp.Clear();

  // C_T = m_0/2.0 . (J + 1/J) . Cinv
  ctemp.Update((m / 2.0), Cinv);

}  // SetupCthermo()


/*----------------------------------------------------------------------*
 | calculates stress-temperature modulus m_0                 dano 09/13 |
 *----------------------------------------------------------------------*/
double MAT::ThermoPlasticHyperElast::STModulus()
{
  // m_0 := -(2 . mu + 3 . lambda) . alpha_T = - 3 . bulk . alpha_T

  // initialise the parameters for the lame constants
  const double ym = params_->youngs_;                 // Young's modulus
  const double nu = params_->poissonratio_;           // Poisson's ratio
  const double bulk = ym / (3.0 * (1.0 - 2.0 * nu));  // bulk modulus
  const double cte = params_->cte_;

  // stress-temperature modulus
  const double stmodulus = (-1.0) * 3.0 * bulk * cte;

  return stmodulus;

}  // STModulus()


/*---------------------------------------------------------------------*
 | return names of visualization data (public)                         |
 *---------------------------------------------------------------------*/
void MAT::ThermoPlasticHyperElast::VisNames(std::map<std::string, int>& names)
{
  std::string accumulatedstrain = "accumulatedstrain";
  names[accumulatedstrain] = 1;  // scalar

  std::string mechdiss = "mechdiss";
  names[mechdiss] = 1;  // scalar

  std::string thrplheating = "thrplheating";
  names[thrplheating] = 1;  // scalar
}  // VisNames()


/*---------------------------------------------------------------------*
 | return visualization data (public)                                  |
 *---------------------------------------------------------------------*/
bool MAT::ThermoPlasticHyperElast::VisData(
    const std::string& name, std::vector<double>& data, int numgp, int eleID)
{
  // accumulated strain
  if (name == "accumulatedstrain")
  {
    if ((int)data.size() != 1) dserror("size mismatch");
    double temp = 0.0;
    for (int iter = 0; iter < numgp; iter++) temp += AccumulatedStrain(iter);
    data[0] = temp / numgp;
  }

  // mechanical dissipation
  if (name == "mechdiss")
  {
    if ((int)data.size() != 1) dserror("size mismatch");
    double temp = 0.0;
    for (int iter = 0; iter < numgp; iter++) temp += MechDiss(iter);
    data[0] = temp / numgp;
  }

  // thermoplastic heating term
  if (name == "thrplheating")
  {
    if ((int)data.size() != 1) dserror("size mismatch");
    double temp = 0.0;
    for (int iter = 0; iter < numgp; iter++) temp += ThermoPlastHeating(iter);
    data[0] = temp / numgp;
  }

  return true;
}  // VisData()


/*---------------------------------------------------------------------*
 | finite difference check for the material tangent.        dano 12/13 |
 | Meant for debugging only! (public)                                  |
 *---------------------------------------------------------------------*/
void MAT::ThermoPlasticHyperElast::FDCheck(
    LINALG::Matrix<NUM_STRESS_3D, 1>& stress,  // updated stress sigma_n+1
    LINALG::Matrix<NUM_STRESS_3D, NUM_STRESS_3D>&
        cmat,  // material tangent calculated with FD of stresses
    LINALG::Matrix<NUM_STRESS_3D, NUM_STRESS_3D>&
        cmatFD,  // material tangent calculated with FD of stresses
    const LINALG::Matrix<1, 1>& Ntemp, Teuchos::ParameterList& params)
{
  // teste 2dS/dC see linearisation according to Holzapfel
  // extract F and Cinv from params
  LINALG::Matrix<3, 3> defgrd = params.get<LINALG::Matrix<3, 3>>("defgrd");
  LINALG::Matrix<6, 1> Cinv_vct = params.get<LINALG::Matrix<6, 1>>("Cinv_vct");

  // calculate the right Cauchy Green (RCG) deformation tensor and its inverse
  LINALG::Matrix<3, 3> RCG_disturb(false);
  RCG_disturb.MultiplyTN(defgrd, defgrd);

  // value of disturbance
  const double delta = 1.0e-8;
  // disturb the respective strain quantities
  for (int i = 0; i < 3; ++i)
  {
    for (int k = 0; k < 3; ++k)
    {
      printf("-------------------------------------\n");
      printf("-------------------------------------\n");
      printf("STRAIN term %d\n", k);

      RCG_disturb(i, k) += delta / 2.0;
      RCG_disturb(k, i) += delta / 2.0;

      // calculate Jacobi-determinant of disturbed RCG
      double detRCG_disturb = RCG_disturb.Determinant();
      double J_disturb = sqrt(detRCG_disturb);
      LINALG::Matrix<3, 3> invRCG_disturb;
      invRCG_disturb.Invert(RCG_disturb);
      // use vector-notation
      LINALG::Matrix<6, 1> disturb_Cinv_vct(false);
      disturb_Cinv_vct(0) = invRCG_disturb(0, 0);
      disturb_Cinv_vct(1) = invRCG_disturb(1, 1);
      disturb_Cinv_vct(2) = invRCG_disturb(2, 2);
      disturb_Cinv_vct(3) = invRCG_disturb(0, 1);
      disturb_Cinv_vct(4) = invRCG_disturb(1, 2);
      disturb_Cinv_vct(5) = invRCG_disturb(2, 0);

      // calculate the temperature difference
      LINALG::Matrix<1, 1> init(false);
      init(0, 0) = params_->inittemp_;
      // Delta T = T - T_0
      LINALG::Matrix<1, 1> deltaT(false);
      deltaT.Update(1.0, Ntemp, (-1.0), init);

      // temperature-dependent stress temperature modulus
      // m = m(J) = m_0 .(J+1)/J = m_0 . (J + 1/J)
      double m_0 = STModulus();
      double m = m_0 * (J_disturb + 1.0 / J_disturb);
      // in case of testing only the dJ/dd, use undisturbed m, but disturb_Cinv_vct
      // double J = defgrd.Determinant();
      // double m = m_0 * (J + 1.0 / J);
      // clear the material tangent
      LINALG::Matrix<NUM_STRESS_3D, 1> disturb_ctemp(true);
      disturb_ctemp.Clear();
      // C_T = m_0/2.0 . (J + 1/J) . Cinv
      disturb_ctemp.Update((m / 2.0), disturb_Cinv_vct);
      // in case of testing only factor, use undisturbed Cinv_vct
      // disturb_ctemp.Update( (m/2.0), Cinv_vct, 0.0);

      // calculate thermal stresses
      // tau = ctemp_AK . Delta T = m_0/2 . (J + 1/J) . I . Delta T
      // pull-back of Kirchhoff-stresses to PK2-stresses
      // PK2 = F^{-1} . tau . F^{-T}
      // --> PK2 = ctemp . Delta T = m_0/2 . (J + 1/J). Cinv . Delta T
      // initialise disturbed total stresses
      LINALG::Matrix<NUM_STRESS_3D, 1> disturb_stresstemp(true);
      disturb_stresstemp.MultiplyNN(disturb_ctemp, deltaT);
      // in case of testing only disturb_ctemp, ignore deltaT
      // disturb_stresstemp.Update(disturb_ctemp);

#ifdef DEBUGMATERIAL
      std::cout << std::scientific;
      std::cout << "Cinv_vct\n " << Cinv_vct << std::endl;
      std::cout << "disturb_Cinv_vct\n " << disturb_Cinv_vct << std::endl;
      std::cout << "Jacobi determinant disturb = " << J_disturb << std::endl;
      std::cout << "Jacobi determinant = " << J << std::endl;
      std::cout << "deltaT = " << deltaT << std::endl;
      std::cout << "Ntemp = " << Ntemp << std::endl;
      std::cout << "inittemp = " << inittemp << std::endl;
      std::cout << "m DT = " << m * deltaT(0, 0) << std::endl;
      std::cout << "disturb_ctemp\n " << disturb_ctemp << std::endl;
      std::cout << "stress\n " << stress << std::endl;
      std::cout << "disturb_stresstemp\n " << disturb_stresstemp << std::endl;
#endif  // DEBUGMATERIAL

      // be careful we save the disturbed RCG in tensor notation, i.e. (3x3)
      // to insert the corresponding terms in cmat (6x6) copy the terms to
      // their correct position using array VOIGT3X3SYM_ as is done, e.g. in
      // neohooke or so3_plast
      double array[3][3] = {{0, 3, 5}, {3, 1, 4}, {5, 4, 2}};

      for (int stress_comp = 0; stress_comp < 6; stress_comp++)
      {
        // build the finite difference tangent
        cmatFD(stress_comp, array[i][k]) = 0.0;
        // scale with factor 2 due to comparison with cmat_T and cmat_T = 2 dSvol/dC
        cmatFD(stress_comp, array[i][k]) +=
            2.0 * ((disturb_stresstemp(stress_comp) / (delta)-stress(stress_comp) / (delta)));

        std::cout << i << k << stress_comp << "fd: "
                  << 2.0 *
                         ((disturb_stresstemp(stress_comp) / (delta)-stress(stress_comp) / (delta)))
                  << "ref: " << cmat(stress_comp, array[i][k]) << std::endl;
      }
      // undisturb the respective strain quantities (disturbstrain=strain)
      RCG_disturb(i, k) -= delta / 2.0;
      RCG_disturb(k, i) -= delta / 2.0;
    }  // loop stresses

  }  // loop strains

}  // FDCheck()


/*----------------------------------------------------------------------*/
