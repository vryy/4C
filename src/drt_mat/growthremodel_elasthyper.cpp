/*----------------------------------------------------------------------*/
/*!
\brief This file is used to manage the homogenized constraint mixture during growth and remodeling

The input line should read
MAT 1 MAT_GrowthRemodel_ElastHyper NUMMATRF 2 MATIDSRF 11 21 NUMMATEL3D 1 MATIDSEL3D 31 NUMMATEL2D 1
MATIDSEL2D 41 MATIDELPENALTY 51 ELMASSFRAC 0.23 DENS 1050 PRESTRETCHELASTINCIR 1.34
PRESTRETCHELASTINAX 1.25 THICKNESS 1.4099822832e-3 MEANPRESSURE 13332.2668 RADIUS 1e-2 DAMAGE 1
GROWTHTYPE 1 LOCTIMEINT 1 MEMBRANE 0

\level 3

\maintainer Fabian Braeu
*/
/*----------------------------------------------------------------------*/
/* headers */
#include "growthremodel_elasthyper.H"
#include "../drt_lib/standardtypes_cpp.H"
#include "../drt_matelast/elast_summand.H"
#include "../drt_matelast/elast_remodelfiber.H"
#include "../drt_lib/drt_linedefinition.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_mat/matpar_bundle.H"
#include "../drt_mat/material_service.H"
#include "Epetra_SerialDenseSolver.h"
#include "../drt_matelast/elast_isoneohooke.H"
#include "../drt_matelast/elast_volsussmanbathe.H"


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::PAR::GrowthRemodel_ElastHyper::GrowthRemodel_ElastHyper(
    Teuchos::RCP<MAT::PAR::Material> matdata)
    : Parameter(matdata),
      nummat_remodelfiber_(matdata->GetInt("NUMMATRF")),
      nummat_elastiniso_(matdata->GetInt("NUMMATEL3D")),
      nummat_elastinmem_(matdata->GetInt("NUMMATEL2D")),
      matids_remodelfiber_(matdata->Get<std::vector<int>>("MATIDSRF")),
      matids_elastiniso_(matdata->Get<std::vector<int>>("MATIDSEL3D")),
      matids_elastinmem_(matdata->Get<std::vector<int>>("MATIDSEL2D")),
      matid_penalty_(matdata->GetInt("MATIDELPENALTY")),
      init_w_el_(matdata->GetDouble("ELMASSFRAC")),
      density_(matdata->GetDouble("DENS")),
      lamb_prestretch_cir_(matdata->GetDouble("PRESTRETCHELASTINCIR")),
      lamb_prestretch_ax_(matdata->GetDouble("PRESTRETCHELASTINAX")),
      t_ref_(matdata->GetDouble("THICKNESS")),
      p_mean_(matdata->GetDouble("MEANPRESSURE")),
      ri_(matdata->GetDouble("RADIUS")),
      damage_(matdata->GetInt("DAMAGE")),
      growthtype_(matdata->GetInt("GROWTHTYPE")),
      loctimeint_(matdata->GetInt("LOCTIMEINT")),
      membrane_(matdata->GetInt("MEMBRANE")),
      cylinder_(matdata->GetInt("CYLINDER"))
{
  // check if sizes fit
  if (nummat_remodelfiber_ != (int)matids_remodelfiber_->size())
    dserror(
        "number of remodelfiber materials %d does not fit to size of remodelfiber material vector "
        "%d",
        nummat_remodelfiber_, matids_remodelfiber_->size());

  if (nummat_elastinmem_ != (int)matids_elastinmem_->size())
    dserror("number of elastin materials %d does not fit to size of elastin material vector %d",
        nummat_elastinmem_, matids_elastinmem_->size());

  if (membrane_ == 1)
  {
    if ((growthtype_ != 1) || (loctimeint_ != 0))
      dserror(
          "using membrane elements you can only simulate anisotropic growth in thickness direction"
          "and solve the local evolution equations with a Forward Euler Method");
  }
  else
  {
    if (nummat_elastiniso_ != (int)matids_elastiniso_->size())
      dserror("number of elastin materials %d does not fit to size of elastin material vector %d",
          nummat_elastiniso_, matids_elastiniso_->size());
    if (nummat_elastiniso_ == 0) dserror("you have to set a 3D elastin material");
    if (matid_penalty_ == -1) dserror("you have to set a volumetric penalty material");
    if ((p_mean_ == -1) || (ri_ == -1) || (t_ref_ == -1))
      dserror(
          "you have to set the mean blood pressure, the inner radius of the vessel and thickness "
          "of the vessel wall");
  }

  if (cylinder_ != -1 && cylinder_ != 1 && cylinder_ != 2 && cylinder_ != 3)
    dserror(
        "The parameter CYLINDER has to be either 1, 2 or 3. If you have defined a fiber direction "
        "in your"
        ".dat file then just skip this parameter");
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<MAT::Material> MAT::PAR::GrowthRemodel_ElastHyper::CreateMaterial()
{
  return Teuchos::rcp(new MAT::GrowthRemodel_ElastHyper(this));
}


MAT::GrowthRemodel_ElastHyperType MAT::GrowthRemodel_ElastHyperType::instance_;


DRT::ParObject* MAT::GrowthRemodel_ElastHyperType::Create(const std::vector<char>& data)
{
  MAT::GrowthRemodel_ElastHyper* gr_elhy = new MAT::GrowthRemodel_ElastHyper();
  gr_elhy->Unpack(data);

  return gr_elhy;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::GrowthRemodel_ElastHyper::GrowthRemodel_ElastHyper()
    : params_(NULL),
      potsumrf_(0),
      potsumeliso_(0),
      potsumelmem_(0),
      potsumelpenalty_(0),
      t_tot_(0),
      nr_rf_tot_(0)
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::GrowthRemodel_ElastHyper::GrowthRemodel_ElastHyper(MAT::PAR::GrowthRemodel_ElastHyper* params)
    : params_(params),
      potsumrf_(0),
      potsumeliso_(0),
      potsumelmem_(0),
      potsumelpenalty_(0),
      t_tot_(0),
      nr_rf_tot_(0)
{
  // make sure the referenced materials in material list have quick access parameters
  std::vector<int>::const_iterator m;

  // RemodelFiber
  for (m = params_->matids_remodelfiber_->begin(); m != params_->matids_remodelfiber_->end(); ++m)
  {
    const int matid = *m;
    Teuchos::RCP<MAT::ELASTIC::RemodelFiber> sum =
        Teuchos::rcp_static_cast<MAT::ELASTIC::RemodelFiber>(MAT::ELASTIC::Summand::Factory(matid));
    if (sum == Teuchos::null) dserror("Failed to allocate");
    potsumrf_.push_back(sum);
  }

  // 2d Elastin matrix
  for (m = params_->matids_elastinmem_->begin(); m != params_->matids_elastinmem_->end(); ++m)
  {
    const int matid = *m;
    Teuchos::RCP<MAT::ELASTIC::Summand> sum = MAT::ELASTIC::Summand::Factory(matid);
    if (sum == Teuchos::null) dserror("Failed to allocate");
    if (sum->MaterialType() != INPAR::MAT::mes_isoneohooke)
      dserror(
          "2D Elastin Material: So far, you have to use a IsoNeoHooke material as the "
          "prestretching algorithm needs it. "
          "The prestretching algorithm can easily be expanded to other materials!");
    potsumelmem_.push_back(sum);
  }

  if (params_->membrane_ != 1)
  {
    // 3d Elastin matrix
    for (m = params_->matids_elastiniso_->begin(); m != params_->matids_elastiniso_->end(); ++m)
    {
      const int matid = *m;
      Teuchos::RCP<MAT::ELASTIC::Summand> sum = MAT::ELASTIC::Summand::Factory(matid);
      if (sum == Teuchos::null) dserror("Failed to allocate");
      if (sum->MaterialType() != INPAR::MAT::mes_isoneohooke)
        dserror(
            "3D Elastin Material: So far, you have to use an IsoNeoHooke material as the "
            "prestretching algorithm needs it"
            "The prestretching algorithm can easily be expanded to other materials!");
      potsumeliso_.push_back(sum);
    }

    // VolPenalty
    Teuchos::RCP<MAT::ELASTIC::Summand> sum =
        MAT::ELASTIC::Summand::Factory(params_->matid_penalty_);
    if (sum == Teuchos::null) dserror("Failed to allocate");
    if (sum->MaterialType() != INPAR::MAT::mes_volsussmanbathe)
      dserror(
          "Volumetric Penalty Material: So far, you have to use a CoupNeoHooke material as the "
          "prestretching algorithm needs it. "
          "This can easily be expanded to other materials!");
    potsumelpenalty_ = sum;
  }

  // initialize total simulation time
  t_tot_ = 0.0;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::GrowthRemodel_ElastHyper::Pack(DRT::PackBuffer& data) const
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


  // mass fraction of elastin
  AddtoPack(data, cur_rho_el_);
  AddtoPack(data, init_rho_el_);

  AddtoPack(data, v_);
  AddtoPack(data, gp_ax_);
  AddtoPack(data, gp_rad_);
  AddtoPack(data, AcirM_);
  AddtoPack(data, AaxM_);
  AddtoPack(data, AradM_);
  AddtoPack(data, Aradv_);
  AddtoPack(data, AplM_);
  AddtoPack(data, AgM_);
  AddtoPack(data, GM_);
  AddtoPack(data, radaxicirc_);
  AddtoPack(data, mue_frac_);
  AddtoPack(data, setup_);
  AddtoPack(data, nr_rf_tot_);

  if (params_ != NULL)  // summands are not accessible in postprocessing mode
  {
    // loop map of associated potential summands
    for (unsigned int p = 0; p < potsumrf_.size(); ++p) potsumrf_[p]->PackSummand(data);

    // loop map of associated potential summands
    for (unsigned int p = 0; p < potsumelmem_.size(); ++p) potsumelmem_[p]->PackSummand(data);

    if (params_->membrane_ != 1)
    {
      // loop map of associated potential summands
      for (unsigned int p = 0; p < potsumeliso_.size(); ++p) potsumeliso_[p]->PackSummand(data);

      potsumelpenalty_->PackSummand(data);
    }
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::GrowthRemodel_ElastHyper::Unpack(const std::vector<char>& data)
{
  // make sure we have a pristine material
  params_ = NULL;
  potsumrf_.clear();
  potsumeliso_.clear();
  potsumelmem_.clear();

  std::vector<char>::size_type position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position, data, type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");

  // matid and recover params_
  int matid;
  ExtractfromPack(position, data, matid);
  if (DRT::Problem::Instance()->Materials() != Teuchos::null)
  {
    if (DRT::Problem::Instance()->Materials()->Num() != 0)
    {
      const unsigned int probinst = DRT::Problem::Instance()->Materials()->GetReadFromProblem();
      MAT::PAR::Parameter* mat =
          DRT::Problem::Instance(probinst)->Materials()->ParameterById(matid);
      if (mat->Type() == MaterialType())
        params_ = static_cast<MAT::PAR::GrowthRemodel_ElastHyper*>(mat);
      else
        dserror("Type of parameter material %d does not fit to calling type %d", mat->Type(),
            MaterialType());
    }
  }


  // mass fractions of elastin and ground matrix
  ExtractfromPack(position, data, cur_rho_el_);
  ExtractfromPack(position, data, init_rho_el_);

  ExtractfromPack(position, data, v_);
  ExtractfromPack(position, data, gp_ax_);
  ExtractfromPack(position, data, gp_rad_);
  ExtractfromPack(position, data, AcirM_);
  ExtractfromPack(position, data, AaxM_);
  ExtractfromPack(position, data, AradM_);
  ExtractfromPack(position, data, Aradv_);
  ExtractfromPack(position, data, AplM_);
  ExtractfromPack(position, data, AgM_);
  ExtractfromPack(position, data, GM_);
  ExtractfromPack(position, data, radaxicirc_);
  ExtractfromPack(position, data, mue_frac_);
  ExtractfromPack(position, data, setup_);
  ExtractfromPack(position, data, nr_rf_tot_);

  if (params_ != NULL)  // summands are not accessible in postprocessing mode
  {
    // make sure the referenced materials in material list have quick access parameters
    std::vector<int>::const_iterator m;

    // RemodelFiber
    for (m = params_->matids_remodelfiber_->begin(); m != params_->matids_remodelfiber_->end(); ++m)
    {
      const int matid = *m;
      Teuchos::RCP<MAT::ELASTIC::RemodelFiber> sum =
          Teuchos::rcp_static_cast<MAT::ELASTIC::RemodelFiber>(
              MAT::ELASTIC::Summand::Factory(matid));
      if (sum == Teuchos::null) dserror("Failed to allocate");
      potsumrf_.push_back(sum);
    }
    // loop map of associated potential summands
    for (unsigned int p = 0; p < potsumrf_.size(); ++p) potsumrf_[p]->UnpackSummand(data, position);

    // 2D Elastin matrix
    for (m = params_->matids_elastinmem_->begin(); m != params_->matids_elastinmem_->end(); ++m)
    {
      const int matid = *m;
      Teuchos::RCP<MAT::ELASTIC::Summand> sum = MAT::ELASTIC::Summand::Factory(matid);
      if (sum == Teuchos::null) dserror("Failed to allocate");
      if (sum->MaterialType() != INPAR::MAT::mes_isoneohooke)
        dserror(
            "2D Elastin Material: So far, you have to use a IsoNeoHooke material as the "
            "prestretching algorithm needs it. "
            "This can easily be expanded to other materials!");
      potsumelmem_.push_back(sum);
    }
    // loop map of associated potential summands
    for (unsigned int p = 0; p < potsumelmem_.size(); ++p)
      potsumelmem_[p]->UnpackSummand(data, position);

    if (params_->membrane_ != 1)
    {
      // 3D Elastin matrix
      for (m = params_->matids_elastiniso_->begin(); m != params_->matids_elastiniso_->end(); ++m)
      {
        const int matid = *m;
        Teuchos::RCP<MAT::ELASTIC::Summand> sum = MAT::ELASTIC::Summand::Factory(matid);
        if (sum == Teuchos::null) dserror("Failed to allocate");
        if (sum->MaterialType() != INPAR::MAT::mes_isoneohooke)
          dserror(
              "3D Elastin Material: So far, you have to use an IsoNeoHooke material as the "
              "prestretching algorithm needs it"
              "This can easily be expanded to other materials!");
        potsumeliso_.push_back(sum);
      }
      // loop map of associated potential summands
      for (unsigned int p = 0; p < potsumeliso_.size(); ++p)
        potsumeliso_[p]->UnpackSummand(data, position);

      // VolPenalty
      Teuchos::RCP<MAT::ELASTIC::Summand> sum =
          MAT::ELASTIC::Summand::Factory(params_->matid_penalty_);
      if (sum == Teuchos::null) dserror("Failed to allocate");
      if (sum->MaterialType() != INPAR::MAT::mes_volsussmanbathe)
        dserror(
            "Volumetric Penalty Material: So far, you have to use a CoupNeoHooke material as the "
            "prestretching algorithm needs it. "
            "This can easily be expanded to other materials!");
      potsumelpenalty_ = sum;
      potsumelpenalty_->UnpackSummand(data, position);
    }

    // in the postprocessing mode, we do not unpack everything we have packed
    // -> position check cannot be done in this case
    if (position != data.size())
      dserror("Mismatch in size of data %d <-> %d", data.size(), position);
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::GrowthRemodel_ElastHyper::Setup(int numgp, DRT::INPUT::LineDefinition* linedef)
{
  // Initialize some variables
  v_.resize(numgp, 1.0);
  gp_ax_.resize(numgp, 0.0);
  gp_rad_.resize(numgp, 0.0);
  cur_rho_el_.resize(numgp);
  init_rho_el_.resize(numgp);
  GM_.resize(numgp, LINALG::Matrix<3, 3>(true));
  setup_.resize(numgp, 1);

  for (int gp = 0; gp < numgp; ++gp)
  {
    init_rho_el_[gp] = params_->density_ * params_->init_w_el_;
    cur_rho_el_[gp] = init_rho_el_[gp];
  }

  // Setup summands
  // remodelfiber
  for (unsigned int p = 0; p < potsumrf_.size(); ++p)
    potsumrf_[p]->Setup(numgp, params_->density_, linedef);

  // 2D elastin matrix
  for (unsigned int p = 0; p < potsumelmem_.size(); ++p) potsumelmem_[p]->Setup(linedef);

  if (params_->membrane_ != 1)
  {
    // 3D elastin matrix
    for (unsigned int p = 0; p < potsumeliso_.size(); ++p) potsumeliso_[p]->Setup(linedef);

    // volpenalty
    potsumelpenalty_->Setup(linedef);
  }

  // Setup circumferential, radial and axial structural tensor
  SetupAxiCirRadStructuralTensor(linedef);

  // Setup growth tensors in the case of anisotropic growth (--> growth in thickness direction)
  if (params_->growthtype_ == 1) SetupAnisoGrowthTensors();

  // variable which is multiplied with the material parameter of the elastin sheets
  mue_frac_.resize(numgp, 1.0);

  // total number of remodel fibers
  for (unsigned p = 0; p < potsumrf_.size(); ++p) nr_rf_tot_ += potsumrf_[p]->GetNumFibers();

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::GrowthRemodel_ElastHyper::SetupAxiCirRadStructuralTensor(
    DRT::INPUT::LineDefinition* linedef)
{
  // CIR-AXI-RAD nomenclature
  if (linedef->HaveNamed("RAD") and linedef->HaveNamed("AXI") and linedef->HaveNamed("CIR"))
  {
    // Read in of data
    // read local (cylindrical) cosy-directions at current element
    LINALG::Matrix<3, 1> dir;

    // Axial direction
    ReadDir(linedef, "AXI", dir);
    for (int i = 0; i < 3; ++i) radaxicirc_(i, 1) = dir(i);
    AaxM_.MultiplyNT(1.0, dir, dir, 0.0);

    // Circumferential direction
    ReadDir(linedef, "CIR", dir);
    for (int i = 0; i < 3; ++i) radaxicirc_(i, 2) = dir(i);
    AcirM_.MultiplyNT(1.0, dir, dir, 0.0);

    // Radial direction
    ReadDir(linedef, "RAD", dir);
    for (int i = 0; i < 3; ++i) radaxicirc_(i, 0) = dir(i);
    AradM_.MultiplyNT(1.0, dir, dir, 0.0);

    // radial structural tensor in "stress-like" Voigt notation
    MatrixtoStressLikeVoigtNotation(AradM_, Aradv_);
  }
  // Cylinder flag defined in .dat file, dummy AXI, CIR and RAD directions are written till
  // location of element center in reference configuration is available
  else if (params_->cylinder_ == 1 || params_->cylinder_ == 2 || params_->cylinder_ == 3)
  {
    radaxicirc_(0, 0) = radaxicirc_(1, 1) = radaxicirc_(2, 2) = 1.0;
    AaxM_(0, 0) = 1.0;
    AcirM_(1, 1) = 1.0;
    AradM_(2, 2) = 1.0;
    Aradv_(2) = 1.0;
  }
  // No AXI CIR RAD-direction defined in .dat file and additionally no cylinder flag was set
  else
    dserror(
        "Homogenized Constrained Mixture Model can so far only be used by defining AXI-, CIR- and "
        "RAD-direction in the .dat file or by defining the Cylinder flag!");

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::GrowthRemodel_ElastHyper::SetupAnisoGrowthTensors()
{
  AgM_.Update(1.0, AradM_, 0.0);

  // structural tensor of the plane in which all fibers are located
  LINALG::Matrix<3, 3> id(true);
  for (int i = 0; i < 3; ++i) id(i, i) = 1.0;
  AplM_.Update(1.0, AradM_, 0.0);
  AplM_.Update(1.0, id, -1.0);

  return;
}


/*----------------------------------------------------------------------*
 * Function which reads in the AXI CIR RAD directions
 *----------------------------------------------------------------------*/
void MAT::GrowthRemodel_ElastHyper::ReadDir(
    DRT::INPUT::LineDefinition* linedef, std::string specifier, LINALG::Matrix<3, 1>& dir)
{
  std::vector<double> fiber;
  linedef->ExtractDoubleVector(specifier, fiber);
  double fnorm = 0.;
  // normalization
  for (int i = 0; i < 3; ++i)
  {
    fnorm += fiber[i] * fiber[i];
  }
  fnorm = sqrt(fnorm);

  // fill final normalized vector
  for (int i = 0; i < 3; ++i) dir(i) = fiber[i] / fnorm;

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::GrowthRemodel_ElastHyper::SetupAxiCirRadCylinder(
    LINALG::Matrix<1, 3> elecenter, double const dt)
{
  // Clear dummy directions
  radaxicirc_.Clear();
  LINALG::Matrix<3, 1> axdir;
  LINALG::Matrix<3, 1> raddir;
  LINALG::Matrix<3, 1> cirdir;

  // Axial direction
  radaxicirc_(params_->cylinder_ - 1, 1) = 1.0;
  for (int i = 0; i < 3; ++i) axdir(i) = radaxicirc_(i, 1);
  AaxM_.MultiplyNT(1.0, axdir, axdir, 0.0);

  // Radial direction
  elecenter(0, params_->cylinder_ - 1) = 0.0;
  raddir.UpdateT(1. / elecenter.Norm2(), elecenter, 0.0);
  for (int i = 0; i < 3; ++i) radaxicirc_(i, 0) = raddir(i);
  AradM_.MultiplyNT(1.0, raddir, raddir, 0.0);

  // Circumferential direction
  cirdir.CrossProduct(axdir, raddir);
  for (int i = 0; i < 3; ++i) radaxicirc_(i, 2) = cirdir(i);
  AcirM_.MultiplyNT(1.0, cirdir, cirdir, 0.0);

  // radial structural tensor in "stress-like" Voigt notation
  MatrixtoStressLikeVoigtNotation(AradM_, Aradv_);

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::GrowthRemodel_ElastHyper::Update(LINALG::Matrix<3, 3> const& defgrd, int const gp,
    Teuchos::ParameterList& params, int const eleGID)
{
  // Update individual volume of elastin
  static double dt_pre = params.get<double>("delta time");
  if (params_->damage_ == 1) EvaluateElastinDamage(dt_pre);

  // loop map of associated potential summands
  // remodelfiber
  for (unsigned int p = 0; p < potsumrf_.size(); ++p) potsumrf_[p]->Update();

  // 2D elastin matrix
  for (unsigned int p = 0; p < potsumelmem_.size(); ++p) potsumelmem_[p]->Update();

  if (params_->membrane_ != 1)
  {
    // 3D elastin matrix
    for (unsigned int p = 0; p < potsumeliso_.size(); ++p) potsumeliso_[p]->Update();

    // volpenalty
    potsumelpenalty_->Update();
  }


  // build inelastic growth deformation gradient
  LINALG::Matrix<3, 3> FgM(true);
  LINALG::Matrix<3, 3> iFgM(true);
  LINALG::Matrix<3, 3> dFgdrhoM(true);
  LINALG::Matrix<3, 3> diFgdrhoM(true);
  EvaluateGrowthDefGrad(FgM, iFgM, dFgdrhoM, diFgdrhoM, gp);

  // time step size
  double dt = params.get<double>("delta time");

  switch (params_->loctimeint_)
  {
    case 0:
      for (unsigned int p = 0; p < potsumrf_.size(); ++p)
        potsumrf_[p]->EvaluateGrowthAndRemodelingExpl(defgrd, dt, iFgM, gp, eleGID);
      break;
    case 1:  // do nothing
      break;
    default:
      dserror("LOCTIMEINT has to be either 1 (Backward Euler) or 0 (Forward Euler)");
      break;
  }

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::GrowthRemodel_ElastHyper::EvaluatePrestretch(
    LINALG::Matrix<3, 3> const* const defgrd, int const gp, int const eleGID)
{
  // setup prestretch of elastin in axial and circumferential direction
  GM_[gp].Update(params_->lamb_prestretch_cir_, AcirM_, 0.0);
  GM_[gp].Update(params_->lamb_prestretch_ax_, AaxM_, 1.0);

  // safety check
  if (potsumeliso_.size() != 1)
    dserror(
        "So far, the prestretching routine does only work with ONE 3D isochoric elastin material");

  Teuchos::RCP<MAT::ELASTIC::IsoNeoHooke> matiso;
  Teuchos::RCP<MAT::ELASTIC::VolSussmanBathe> matvol;
  matiso = Teuchos::rcp_dynamic_cast<MAT::ELASTIC::IsoNeoHooke>(potsumeliso_[0]);
  matvol = Teuchos::rcp_dynamic_cast<MAT::ELASTIC::VolSussmanBathe>(potsumelpenalty_);
  double R = 1.0;
  double dRdlamb_pre = 0.0;
  double lamb_pre = 1. / (params_->lamb_prestretch_cir_ * params_->lamb_prestretch_ax_);
  while (fabs(R) > 1.0e-10)
  {
    R = matiso->Mue() * init_rho_el_[gp] *
            std::pow(params_->lamb_prestretch_cir_ * params_->lamb_prestretch_cir_ *
                         params_->lamb_prestretch_ax_ * params_->lamb_prestretch_ax_ * lamb_pre *
                         lamb_pre,
                -2. / 3.) *
            (lamb_pre * lamb_pre -
                (1. / 3.) * (params_->lamb_prestretch_cir_ * params_->lamb_prestretch_cir_ +
                                params_->lamb_prestretch_ax_ * params_->lamb_prestretch_ax_ +
                                lamb_pre * lamb_pre)) +
        matvol->Kappa() * init_rho_el_[gp] *
            ((params_->lamb_prestretch_cir_ * params_->lamb_prestretch_ax_ * lamb_pre) *
                    (params_->lamb_prestretch_cir_ * params_->lamb_prestretch_ax_ * lamb_pre) -
                (params_->lamb_prestretch_cir_ * params_->lamb_prestretch_ax_ * lamb_pre)) +
        ((1.0 - (gp_rad_[gp] - 10.0e-3) / params_->t_ref_) * params_->p_mean_);

    dRdlamb_pre =
        matiso->Mue() * init_rho_el_[gp] *
            (-(4. / 3.) *
                std::pow(params_->lamb_prestretch_cir_ * params_->lamb_prestretch_ax_ * lamb_pre,
                    -7. / 3.) *
                params_->lamb_prestretch_cir_ * params_->lamb_prestretch_ax_) *
            (lamb_pre * lamb_pre -
                (1. / 3.) * (params_->lamb_prestretch_cir_ * params_->lamb_prestretch_cir_ +
                                params_->lamb_prestretch_ax_ * params_->lamb_prestretch_ax_ +
                                lamb_pre * lamb_pre)) +
        matiso->Mue() * init_rho_el_[gp] *
            std::pow(params_->lamb_prestretch_cir_ * params_->lamb_prestretch_cir_ *
                         params_->lamb_prestretch_ax_ * params_->lamb_prestretch_ax_ * lamb_pre *
                         lamb_pre,
                -2. / 3.) *
            (2.0 * lamb_pre - (1. / 3.) * (2.0 * lamb_pre)) +
        matvol->Kappa() * init_rho_el_[gp] *
            (2.0 * (params_->lamb_prestretch_cir_ * params_->lamb_prestretch_ax_ * lamb_pre) *
                    params_->lamb_prestretch_cir_ * params_->lamb_prestretch_ax_ -
                params_->lamb_prestretch_cir_ * params_->lamb_prestretch_ax_);

    lamb_pre = lamb_pre + (-R / dRdlamb_pre);
  }

  // update radial prestretch of elastin
  GM_[gp].Update(lamb_pre, AradM_, 1.0);


  // calculate circumferential residual stress
  double sig = 0.0;
  LINALG::Matrix<6, 1> Acir_strain(true);
  MatrixtoStrainLikeVoigtNotation(AcirM_, Acir_strain);

  // total circumferential Cauchy stress ("membrane stress")
  double sig_tot = (params_->p_mean_ * params_->ri_) / params_->t_ref_;

  // evaluate anisotropic remodel fibers
  LINALG::Matrix<3, 3> CM(true);
  CM.MultiplyTN(1.0, *defgrd, *defgrd, 0.0);
  LINALG::Matrix<6, 1> stressaniso(true);
  LINALG::Matrix<6, 6> cmataniso(true);
  LINALG::Matrix<3, 3> id(true);
  for (int i = 0; i < 3; ++i) id(i, i) = 1.0;
  for (unsigned p = 0; p < potsumrf_.size(); ++p)
  {
    potsumrf_[p]->EvaluateAnisotropicStressCmat(CM, id, cmataniso, stressaniso, gp, 1.0, eleGID);
    sig += stressaniso.Dot(Acir_strain);
  }

  // build stress response and elasticity tensor of 3D material
  LINALG::Matrix<NUM_STRESS_3D, 1> stressiso(true);
  LINALG::Matrix<NUM_STRESS_3D, NUM_STRESS_3D> cmatiso(true);
  LINALG::Matrix<NUM_STRESS_3D, 9> dSdiFgiso(true);
  EvaluateStressCmatIso(defgrd, GM_[gp], stressiso, cmatiso, dSdiFgiso, gp, eleGID);
  sig += stressiso.Dot(Acir_strain);

  // residual stress
  double sig_res = sig_tot - sig;

  // safety check
  if (potsumelmem_.size() != 1)
    dserror("So far, only ONE CoupNeoHooke material can be chosen for the membrane \"material\"");

  // build stress response and elasticity tensor of membrane material
  LINALG::Matrix<6, 1> stressmem(true);
  LINALG::Matrix<6, 6> cmatmem(true);
  LINALG::Matrix<NUM_STRESS_3D, 9> dSdiFgmem(true);
  EvaluateStressCmatMembrane(CM, id, stressmem, cmatmem, dSdiFgmem, gp, eleGID);

  mue_frac_[gp] = sig_res / stressmem.Dot(Acir_strain);

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::GrowthRemodel_ElastHyper::SetupGR3D(LINALG::Matrix<3, 3> const* const defgrd,
    Teuchos::ParameterList& params, const double dt, const int gp, const int eleGID)
{
  LINALG::Matrix<1, 3> gprefecoord(true);  // gp coordinates in reference configuration
  LINALG::Matrix<1, 3> axdir(true);
  LINALG::Matrix<1, 3> raddir(true);
  gprefecoord = params.get<LINALG::Matrix<1, 3>>("gprefecoord");

  if ((params_->cylinder_ == 1 || params_->cylinder_ == 2 || params_->cylinder_ == 3) &&
      (setup_[0] == 1))
  {
    LINALG::Matrix<1, 3> elerefecoord(true);  // element center coordinates in reference system
    elerefecoord = params.get<LINALG::Matrix<1, 3>>("elecenter");
    SetupAxiCirRadCylinder(elerefecoord, dt);
    if (params_->growthtype_ == 1) SetupAnisoGrowthTensors();

    // Update fiber directions with new local coordinate system (radaxicirc_)
    for (unsigned k = 0; k < potsumrf_.size(); ++k) potsumrf_[k]->UpdateFiberDirs(radaxicirc_, dt);
  }

  // Evaluate radial and axial distance between origin and current Gauß-Point
  for (int i = 0; i < 3; ++i)
  {
    axdir(0, i) = radaxicirc_(i, 1);
    raddir(0, i) = radaxicirc_(i, 0);
  }
  gp_ax_[gp] = axdir.Dot(gprefecoord);
  gp_rad_[gp] = raddir.Dot(gprefecoord);

  // TODO: BE CAREFULL! So far, this prestretching procedure is only valid for certain materials and
  // a cylindrical geometry.
  //       The principle of the prestretching routine can easily be adapted to other materials or
  //       general geometries!!!
  EvaluatePrestretch(defgrd, gp, eleGID);

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::GrowthRemodel_ElastHyper::Evaluate(const LINALG::Matrix<3, 3>* defgrd,
    const LINALG::Matrix<6, 1>* glstrain, Teuchos::ParameterList& params,
    LINALG::Matrix<6, 1>* stress, LINALG::Matrix<6, 6>* cmat, const int eleGID)
{
  // save current simulation time (used for the evaluation of elastin degradation)
  t_tot_ = params.get<double>("total time");

  // current gp
  int gp = params.get<int>("gp");

  // time step size
  double dt = params.get<double>("delta time");

  if (setup_[gp] == 1)
  {
    SetupGR3D(defgrd, params, dt, gp, eleGID);
    setup_[gp] = 0;
  }

  // blank resulting quantities
  // ... even if it is an implicit law that cmat is zero upon input
  stress->Clear();
  cmat->Clear();

  // some static variables
  static LINALG::Matrix<NUM_STRESS_3D, 1> stressaniso(true);
  static LINALG::Matrix<NUM_STRESS_3D, NUM_STRESS_3D> cmataniso(true);
  static LINALG::Matrix<NUM_STRESS_3D, NUM_STRESS_3D> cmatanisoadd(true);
  static LINALG::Matrix<3, 3> iFinM(true);
  static LINALG::Matrix<3, 3> CM(true);
  CM.MultiplyTN(1.0, *defgrd, *defgrd, 0.0);
  static LINALG::Matrix<6, 1> stressmem(true);
  static LINALG::Matrix<6, 6> cmatmem(true);
  static LINALG::Matrix<6, 9> dSdiFgmem(true);
  static LINALG::Matrix<NUM_STRESS_3D, 1> stressiso(true);
  static LINALG::Matrix<NUM_STRESS_3D, NUM_STRESS_3D> cmatiso(true);
  static LINALG::Matrix<6, 9> dSdiFgiso(true);
  static LINALG::Matrix<NUM_STRESS_3D, NUM_STRESS_3D> cmatadd(true);

  // build growth deformation gradient
  static LINALG::Matrix<3, 3> FgM(true);
  static LINALG::Matrix<3, 3> iFgM(true);
  static LINALG::Matrix<3, 3> dFgdrhoM(true);
  static LINALG::Matrix<3, 3> diFgdrhoM(true);
  EvaluateGrowthDefGrad(FgM, iFgM, dFgdrhoM, diFgdrhoM, gp);

  // some initialization
  // only do at first time step
  if (t_tot_ == dt)
  {
    // evaluate anisotropic remodel fibers
    for (unsigned p = 0; p < potsumrf_.size(); ++p)
    {
      potsumrf_[p]->EvaluateAnisotropicStressCmat(CM, iFgM, cmataniso, stressaniso, gp, dt, eleGID);
      stress->Update(1.0, stressaniso, 1.0);
      cmat->Update(1.0, cmataniso, 1.0);
    }

    // build stress response and elasticity tensor of 3D material
    iFinM.MultiplyNN(1.0, iFgM, GM_[gp], 0.0);
    EvaluateStressCmatIso(defgrd, iFinM, stressiso, cmatiso, dSdiFgiso, gp, eleGID);
    stress->Update(1.0, stressiso, 1.0);
    cmat->Update(1.0, cmatiso, 1.0);

    // build stress response and elasticity tensor of membrane material
    EvaluateStressCmatMembrane(CM, iFgM, stressmem, cmatmem, dSdiFgmem, gp, eleGID);
    stress->Update(1.0, stressmem, 1.0);
    cmat->Update(1.0, cmatmem, 1.0);

    return;
  }
  // Growth and Remodeling
  else
  {
    switch (params_->loctimeint_)
    {
      case 0:
      {
        // evaluate anisotropic remodel fibers
        for (unsigned p = 0; p < potsumrf_.size(); ++p)
        {
          potsumrf_[p]->EvaluateAnisotropicStressCmat(
              CM, iFgM, cmataniso, stressaniso, gp, dt, eleGID);
          stress->Update(1.0, stressaniso, 1.0);
          cmat->Update(1.0, cmataniso, 1.0);
        }

        // build stress response and elasticity tensor of 3D material
        iFinM.MultiplyNN(1.0, iFgM, GM_[gp], 0.0);
        EvaluateStressCmatIso(defgrd, iFinM, stressiso, cmatiso, dSdiFgiso, gp, eleGID);
        stress->Update(1.0, stressiso, 1.0);
        cmat->Update(1.0, cmatiso, 1.0);

        // build stress response and elasticity tensor of membrane material
        EvaluateStressCmatMembrane(CM, iFgM, stressmem, cmatmem, dSdiFgmem, gp, eleGID);
        stress->Update(1.0, stressmem, 1.0);
        cmat->Update(1.0, cmatmem, 1.0);

        break;
      }
      case 1:
      {
        static LINALG::SerialDenseMatrix K_T(2 * nr_rf_tot_, 2 * nr_rf_tot_, true);
        SolveForRhoLambr(K_T, FgM, iFgM, dFgdrhoM, diFgdrhoM, defgrd, dt, gp, eleGID);

        // split solution vector in dlambda_r/dC and drho/dC
        static std::vector<LINALG::Matrix<1, 6>> drhodC(nr_rf_tot_, LINALG::Matrix<1, 6>(true));
        static std::vector<LINALG::Matrix<1, 6>> dlambrdC(nr_rf_tot_, LINALG::Matrix<1, 6>(true));
        static LINALG::Matrix<1, 6> sum_drhodC(true);
        SolveFordrhodCdlambrdC(drhodC, dlambrdC, sum_drhodC, K_T, iFgM, defgrd, dt, gp, eleGID);

        // Evaluate anisotropic remodel fibers
        int nr_grf_proc = 0;
        for (unsigned p = 0; p < potsumrf_.size(); ++p)
        {
          potsumrf_[p]->EvaluateAnisotropicStressCmat(
              CM, iFgM, cmataniso, stressaniso, gp, dt, eleGID);
          stress->Update(1.0, stressaniso, 1.0);
          cmat->Update(1.0, cmataniso, 1.0);
          potsumrf_[p]->EvaluateAdditionalGrowthRemodelCmat(
              defgrd, nr_grf_proc, iFgM, diFgdrhoM, drhodC, dlambrdC, cmatanisoadd, gp, eleGID);
          cmat->Update(1.0, cmatanisoadd, 1.0);
          nr_grf_proc += potsumrf_[p]->GetNumFibers();
        }

        // Evaluate 3D elastin material
        iFinM.MultiplyNN(1.0, iFgM, GM_[gp], 0.0);
        EvaluateStressCmatIso(defgrd, iFinM, stressiso, cmatiso, dSdiFgiso, gp, eleGID);
        stress->Update(1.0, stressiso, 1.0);
        cmat->Update(1.0, cmatiso, 1.0);

        // Build stress response and elasticity tensor of membrane material
        EvaluateStressCmatMembrane(CM, iFgM, stressmem, cmatmem, dSdiFgmem, gp, eleGID);
        stress->Update(1.0, stressmem, 1.0);
        cmat->Update(1.0, cmatmem, 1.0);

        // Evaluate additional terms for the elasticity tensor
        static LINALG::Matrix<6, 9> dSdiFg_sum(true);
        dSdiFg_sum.Update(1.0, dSdiFgiso, 0.0);
        dSdiFg_sum.Update(1.0, dSdiFgmem, 1.0);
        EvaluateAdditionalCmat(cmatadd, diFgdrhoM, sum_drhodC, dSdiFg_sum, gp);
        cmat->Update(1.0, cmatadd, 1.0);

        break;
      }
      default:
        dserror("LOCTIMEINT has to be 1 (Backward Euler Method) or 0 (Forward Euler Method)");
        break;
    }
  }

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::GrowthRemodel_ElastHyper::SolveForRhoLambr(LINALG::SerialDenseMatrix& K_T,
    LINALG::Matrix<3, 3>& FgM, LINALG::Matrix<3, 3>& iFgM, LINALG::Matrix<3, 3>& dFgdrhoM,
    LINALG::Matrix<3, 3>& diFgdrhoM, LINALG::Matrix<3, 3> const* const defgrd, double const& dt,
    int const gp, int const eleGID)
{
  // allocate some variables (allocate variables only once, therefore static)
  static std::vector<std::vector<double>> dWdrho(nr_rf_tot_, std::vector<double>(nr_rf_tot_, 0.0));
  static std::vector<std::vector<double>> dWdlambr(
      nr_rf_tot_, std::vector<double>(nr_rf_tot_, 0.0));
  static std::vector<double> W(nr_rf_tot_, 0.0);
  static std::vector<std::vector<double>> dEdrho(nr_rf_tot_, std::vector<double>(nr_rf_tot_, 0.0));
  static std::vector<std::vector<double>> dEdlambr(
      nr_rf_tot_, std::vector<double>(nr_rf_tot_, 0.0));
  static std::vector<double> E(nr_rf_tot_, 0.0);
  static Epetra_SerialDenseSolver solver;
  // residual vector of assembled system of equation
  static Epetra_SerialDenseMatrix R(2 * nr_rf_tot_, 1);
  for (unsigned i = 0; i < 2 * nr_rf_tot_; ++i) R(i, 0) = 1.0;
  // solution vector of assembled system of equation
  static Epetra_SerialDenseMatrix dsol(2 * nr_rf_tot_, 1);

  int nr_grf_proc = 0;
  int iter = 0;
  int l = 0;
  while (R.NormOne() > 1.0e-10)
  {
    if (iter != 0)
    {
      // Solve linearized system of equations
      solver.SetMatrix(K_T);
      solver.SetVectors(dsol, R);
      solver.SolveToRefinedSolution(true);
      solver.ApplyRefinement();
      solver.Solve();

      l = 0;
      for (unsigned p = 0; p < potsumrf_.size(); ++p)
      {
        for (unsigned k = 0; k < potsumrf_[p]->GetNumFibers(); ++k)
        {
          // update inelastic fiber stretch and current mass density
          potsumrf_[p]->UpdateGrowthRemodelParameter(dsol(l, 0), dsol(nr_rf_tot_ + l, 0), k, gp);
          ++l;
        }
      }
    }
    // Update growth deformation gradient
    EvaluateGrowthDefGrad(FgM, iFgM, dFgdrhoM, diFgdrhoM, gp);

    // global number of fibers which were already processed
    nr_grf_proc = 0;
    for (unsigned p = 0; p < potsumrf_.size(); ++p)
    {
      potsumrf_[p]->EvaluateDerivativesInternalNewton(defgrd, nr_grf_proc, nr_rf_tot_, gp, dt,
          eleGID, iFgM, dFgdrhoM, diFgdrhoM, dWdrho, dWdlambr, W, dEdrho, dEdlambr, E);
      nr_grf_proc += potsumrf_[p]->GetNumFibers();
    }

    // Assembly
    for (unsigned i = 0; i < nr_rf_tot_; ++i)
    {
      R(i, 0) = -W[i];
      R(nr_rf_tot_ + i, 0) = -E[i];
      for (unsigned j = 0; j < nr_rf_tot_; ++j)
      {
        K_T(i, j) = dWdrho[i][j];
        K_T(i, nr_rf_tot_ + j) = dWdlambr[i][j];
        K_T(nr_rf_tot_ + i, j) = dEdrho[i][j];
        K_T(nr_rf_tot_ + i, nr_rf_tot_ + j) = dEdlambr[i][j];
      }
    }

    if (iter > 95)
    {
      std::cout << "iteration:  " << iter << std::endl
                << std::endl
                << "Vector of residuals: " << std::endl
                << R << std::endl
                << "Matrix of derivatives:" << std::endl
                << K_T;
      std::cout << "=================================================================" << std::endl
                << std::endl;
      if (iter > 100) dserror("Internal Newton (at Gauß-Point) does not converge!");
    }
    ++iter;
  }

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::GrowthRemodel_ElastHyper::SolveFordrhodCdlambrdC(
    std::vector<LINALG::Matrix<1, 6>>& drhodC, std::vector<LINALG::Matrix<1, 6>>& dlambrdC,
    LINALG::Matrix<1, 6>& sum_drhodC, LINALG::SerialDenseMatrix& K_T,
    LINALG::Matrix<3, 3> const& iFgM, LINALG::Matrix<3, 3> const* const defgrd, double const& dt,
    int const gp, int const eleGID) const
{
  static std::vector<LINALG::Matrix<1, 6>> dWdC(nr_rf_tot_, LINALG::Matrix<1, 6>(true));
  static std::vector<LINALG::Matrix<1, 6>> dEdC(nr_rf_tot_, LINALG::Matrix<1, 6>(true));
  int nr_grf_proc = 0;
  for (unsigned p = 0; p < potsumrf_.size(); ++p)
  {
    potsumrf_[p]->EvaluateDerivativesCauchyGreen(
        defgrd, nr_grf_proc, gp, dt, iFgM, dWdC, dEdC, eleGID);
    nr_grf_proc += potsumrf_[p]->GetNumFibers();
  }

  // Assembly
  // Rcmat
  static LINALG::SerialDenseMatrix Rcmat(2 * nr_rf_tot_, 6);
  for (unsigned i = 0; i < nr_rf_tot_; ++i)
    for (unsigned j = 0; j < 6; ++j)
    {
      Rcmat(i, j) = -dWdC[i](0, j);
      Rcmat(nr_rf_tot_ + i, j) = -dEdC[i](0, j);
    }

  // Solve
  static Epetra_SerialDenseSolver solver;
  static LINALG::SerialDenseMatrix dsolcmat(2 * nr_rf_tot_, 6);
  solver.SetMatrix(K_T);
  solver.SetVectors(dsolcmat, Rcmat);
  solver.SolveToRefinedSolution(true);
  solver.ApplyRefinement();
  solver.Solve();

  sum_drhodC.Clear();
  for (unsigned i = 0; i < nr_rf_tot_; ++i)
  {
    for (unsigned j = 0; j < 6; ++j)
    {
      drhodC[i](0, j) = dsolcmat(i, j);
      dlambrdC[i](0, j) = dsolcmat(nr_rf_tot_ + i, j);
    }
    sum_drhodC.Update(1.0, drhodC[i], 1.0);
  }

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::GrowthRemodel_ElastHyper::EvaluateStressCmatIso(LINALG::Matrix<3, 3> const* const defgrd,
    LINALG::Matrix<3, 3> const& iFinM, LINALG::Matrix<6, 1>& stressiso,
    LINALG::Matrix<6, 6>& cmatiso, LINALG::Matrix<6, 9>& dSdiFg, int const gp,
    int const eleGID) const
{
  // clear some variables
  stressiso.Clear();
  cmatiso.Clear();
  dSdiFg.Clear();

  // Evaluate elastin matrix
  // some variables
  static LINALG::Matrix<6, 1> iCinv(true);
  static LINALG::Matrix<6, 1> iCinCiCinv(true);
  static LINALG::Matrix<6, 1> iCv(true);
  static LINALG::Matrix<3, 1> prinv(true);
  static LINALG::Matrix<3, 3> iCinCM(true);
  static LINALG::Matrix<3, 3> iFinCeM(true);
  static LINALG::Matrix<9, 1> CiFin9x1(true);
  LINALG::Matrix<9, 1> CiFinCe9x1(true);
  LINALG::Matrix<9, 1> CiFiniCe9x1(true);

  EvaluateKinQuantElast(defgrd, iFinM, iCinv, iCinCiCinv, iCv, iCinCM, iFinCeM, CiFin9x1,
      CiFinCe9x1, CiFiniCe9x1, prinv, gp);

  LINALG::Matrix<3, 1> dPIe(true);
  LINALG::Matrix<6, 1> ddPIIe(true);
  EvaluateInvariantDerivatives(prinv, dPIe, ddPIIe, gp, eleGID);

  // 2nd Piola Kirchhoff stress factors (according to Holzapfel-Nonlinear Solid Mechanics p. 216)
  static LINALG::Matrix<3, 1> gamma(true);
  // constitutive tensor factors (according to Holzapfel-Nonlinear Solid Mechanics p. 261)
  static LINALG::Matrix<8, 1> delta(true);

  // compose coefficients
  CalculateGammaDelta(gamma, delta, prinv, dPIe, ddPIIe);

  EvaluateIsotropicPrincElast(stressiso, cmatiso, iCinv, iCinCiCinv, iCv, gamma, delta);

  EvaluatedSdiFg(dSdiFg, gamma, delta, iFinM, iCinCM, iCinv, CiFin9x1, CiFinCe9x1, iCinCiCinv,
      CiFiniCe9x1, iCv, iFinCeM, gp);

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::GrowthRemodel_ElastHyper::EvaluateKinQuantElast(LINALG::Matrix<3, 3> const* const defgrd,
    LINALG::Matrix<3, 3> const& iFinM, LINALG::Matrix<6, 1>& iCinv,
    LINALG::Matrix<6, 1>& iCinCiCinv, LINALG::Matrix<6, 1>& iCv, LINALG::Matrix<3, 3>& iCinCM,
    LINALG::Matrix<3, 3>& iFinCeM, LINALG::Matrix<9, 1>& CiFin9x1, LINALG::Matrix<9, 1>& CiFinCe9x1,
    LINALG::Matrix<9, 1>& CiFiniCe9x1, LINALG::Matrix<3, 1>& prinv, const int gp) const
{
  // inverse inelastic right Cauchy-Green
  static LINALG::Matrix<3, 3> iCinM(true);
  iCinM.MultiplyNT(1.0, iFinM, iFinM, 0.0);
  MatrixtoStressLikeVoigtNotation(iCinM, iCinv);

  // inverse right Cauchy-Green
  static LINALG::Matrix<3, 3> iCM(true);
  static LINALG::Matrix<3, 3> CM(true);
  CM.MultiplyTN(1.0, *defgrd, *defgrd, 0.0);
  iCM.Invert(CM);
  MatrixtoStressLikeVoigtNotation(iCM, iCv);

  // C_{in}^{-1} * C * C_{in}^{-1}
  static LINALG::Matrix<3, 3> tmp(true);
  static LINALG::Matrix<3, 3> iCinCiCinM;
  tmp.MultiplyNN(1.0, iCinM, CM, 0.0);
  iCinCiCinM.MultiplyNN(1.0, tmp, iCinM, 0.0);
  MatrixtoStressLikeVoigtNotation(iCinCiCinM, iCinCiCinv);

  // elastic right Cauchy-Green in strain-like Voigt notation.
  tmp.MultiplyNN(1.0, *defgrd, iFinM, 0.0);
  static LINALG::Matrix<3, 3> CeM(true);
  CeM.MultiplyTN(1.0, tmp, tmp, 0.0);
  static LINALG::Matrix<6, 1> Ce_strain(true);
  MatrixtoStrainLikeVoigtNotation(CeM, Ce_strain);

  // principal invariants of elastic right Cauchy-Green strain
  InvariantsPrincipal<MAT::VoigtNotation::strain>(prinv, Ce_strain);

  // C_{in}^{-1} * C
  iCinCM.MultiplyNN(1.0, iCinM, CM, 0.0);

  // F_{in}^{-1} * C_e
  iFinCeM.MultiplyNN(1.0, iFinM, CeM, 0.0);

  // C * F_{in}^{-1}
  static LINALG::Matrix<3, 3> CiFinM(true);
  CiFinM.MultiplyNN(1.0, CM, iFinM, 0.0);
  Matrix3x3to9x1(CiFinM, CiFin9x1);

  // C * F_{in}^{-1} * C_e
  static LINALG::Matrix<3, 3> CiFinCeM(true);
  tmp.MultiplyNN(1.0, CM, iFinM, 0.0);
  CiFinCeM.MultiplyNN(1.0, tmp, CeM, 0.0);
  Matrix3x3to9x1(CiFinCeM, CiFinCe9x1);

  // C * F_{in}^{-1} * C_e^{-1}
  static LINALG::Matrix<3, 3> CiFiniCeM(true);
  static LINALG::Matrix<3, 3> iCeM(true);
  iCeM.Invert(CeM);
  tmp.MultiplyNN(1.0, CM, iFinM, 0.0);
  CiFiniCeM.MultiplyNN(1.0, tmp, iCeM, 0.0);
  Matrix3x3to9x1(CiFiniCeM, CiFiniCe9x1);

  return;
}


/*----------------------------------------------------------------------/
 * Reads derivatives with respect to invariants and modified invariants
 * from all materials of the elasthyper-toolbox                         */
/*----------------------------------------------------------------------*/
void MAT::GrowthRemodel_ElastHyper::EvaluateInvariantDerivatives(LINALG::Matrix<3, 1> const& prinv,
    LINALG::Matrix<3, 1>& dPIw, LINALG::Matrix<6, 1>& ddPIIw, int const gp, int const eleGID) const
{
  // derivatives of principal materials weighted with their mass fraction in the constraint mixture
  static LINALG::Matrix<3, 1> dPgrowthI(true);
  static LINALG::Matrix<6, 1> ddPgrowthII(true);

  // loop map of associated potential summands
  // derivatives of strain energy function w.r.t. principal invariants
  for (unsigned p = 0; p < potsumeliso_.size(); ++p)
  {
    dPgrowthI.Clear();
    ddPgrowthII.Clear();
    potsumeliso_[p]->AddDerivativesPrincipal(dPgrowthI, ddPgrowthII, prinv, eleGID);
    dPIw.Update(cur_rho_el_[gp], dPgrowthI, 1.0);
    ddPIIw.Update(cur_rho_el_[gp], ddPgrowthII, 1.0);
  }

  // derivatives of decoupled (volumetric or isochoric) materials weighted with their mass fraction
  // in the constraint mixture
  static LINALG::Matrix<3, 1> modinv(true);
  InvariantsModified(modinv, prinv);
  LINALG::Matrix<3, 1> dPmodI(true);
  LINALG::Matrix<6, 1> ddPmodII(true);
  for (unsigned p = 0; p < potsumeliso_.size(); ++p)
  {
    dPgrowthI.Clear();
    ddPgrowthII.Clear();
    potsumeliso_[p]->AddDerivativesModified(dPgrowthI, ddPgrowthII, modinv, eleGID);
    dPmodI.Update(cur_rho_el_[gp], dPgrowthI, 1.0);
    ddPmodII.Update(cur_rho_el_[gp], ddPgrowthII, 1.0);
  }

  // volpenalty
  dPgrowthI.Clear();
  ddPgrowthII.Clear();
  potsumelpenalty_->AddDerivativesModified(dPgrowthI, ddPgrowthII, modinv, 0);
  dPmodI.Update(cur_rho_el_[gp], dPgrowthI, 1.0);
  ddPmodII.Update(cur_rho_el_[gp], ddPgrowthII, 1.0);

  // convert decoupled derivatives to principal derivatives
  ConvertModToPrinc(prinv, dPmodI, ddPmodII, dPIw, ddPIIw);

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::GrowthRemodel_ElastHyper::InvariantsModified(
    LINALG::Matrix<3, 1>& modinv, LINALG::Matrix<3, 1> const& prinv) const
{
  // 1st invariant, trace
  modinv(0) = prinv(0) * std::pow(prinv(2), -1. / 3.);
  // 2nd invariant
  modinv(1) = prinv(1) * std::pow(prinv(2), -2. / 3.);
  // J
  modinv(2) = std::pow(prinv(2), 1. / 2.);

  return;
}


/*----------------------------------------------------------------------/
 * Converts derivatives with respect to modified invariants in derivatives
 * with respect to principal invariants                                 */
/*----------------------------------------------------------------------*/
void MAT::GrowthRemodel_ElastHyper::ConvertModToPrinc(LINALG::Matrix<3, 1> const& prinv,
    LINALG::Matrix<3, 1> const& dPmodI, LINALG::Matrix<6, 1> const& ddPmodII,
    LINALG::Matrix<3, 1>& dPI, LINALG::Matrix<6, 1>& ddPII) const
{
  // Conversions to dPI
  dPI(0) += std::pow(prinv(2), -1. / 3.) * dPmodI(0);
  dPI(1) += std::pow(prinv(2), -2. / 3.) * dPmodI(1);
  dPI(2) += 0.5 * std::pow(prinv(2), -0.5) * dPmodI(2) -
            1. / 3. * prinv(0) * std::pow(prinv(2), -4. / 3.) * dPmodI(0) -
            2. / 3. * prinv(1) * std::pow(prinv(2), -5. / 3.) * dPmodI(1);

  // Conversions to ddPII
  ddPII(0) += std::pow(prinv(2), -2. / 3.) * ddPmodII(0);
  ddPII(1) += std::pow(prinv(2), -4. / 3.) * ddPmodII(1);
  ddPII(2) += (1. / 9.) * std::pow(prinv(2), -8. / 3.) * prinv(0) * prinv(0) * ddPmodII(0) +
              (4. / 9.) * prinv(0) * prinv(1) * std::pow(prinv(2), -3.) * ddPmodII(5) -
              (1. / 3.) * std::pow(prinv(2), -11. / 6.) * prinv(0) * ddPmodII(4) +
              (4. / 9.) * std::pow(prinv(2), -7. / 3.) * prinv(0) * dPmodI(0) +
              (4. / 9.) * std::pow(prinv(2), -10. / 3.) * prinv(1) * prinv(1) * ddPmodII(1) -
              (2. / 3.) * std::pow(prinv(2), -13. / 6.) * prinv(1) * ddPmodII(3) +
              (10. / 9.) * std::pow(prinv(2), -8. / 3.) * prinv(1) * dPmodI(1) +
              0.25 * std::pow(prinv(2), -1.) * ddPmodII(2) -
              0.25 * std::pow(prinv(2), -1.5) * dPmodI(2);
  ddPII(3) += -(1. / 3.) * std::pow(prinv(2), -2.) * prinv(0) * ddPmodII(5) -
              (2. / 3.) * std::pow(prinv(2), -7. / 3.) * prinv(1) * ddPmodII(1) +
              0.5 * std::pow(prinv(2), -7. / 6.) * ddPmodII(3) -
              (2. / 3.) * std::pow(prinv(2), -5. / 3.) * dPmodI(1);
  ddPII(4) += -(1. / 3.) * std::pow(prinv(2), -5. / 3.) * prinv(0) * ddPmodII(0) -
              (2. / 3.) * std::pow(prinv(2), -2.) * prinv(1) * ddPmodII(5) +
              0.5 * std::pow(prinv(2), -5. / 6.) * ddPmodII(4) -
              (1. / 3.) * std::pow(prinv(2), -4. / 3.) * dPmodI(0);
  ddPII(5) += std::pow(prinv(2), -1.) * ddPmodII(5);

  return;
}


/*----------------------------------------------------------------------*
 * Calculate the coefficients gamma and delta from the partial          *
 * derivatives w.r.t. invariants                                        *
 *----------------------------------------------------------------------*/
void MAT::GrowthRemodel_ElastHyper::CalculateGammaDelta(LINALG::Matrix<3, 1>& gamma,
    LINALG::Matrix<8, 1>& delta, LINALG::Matrix<3, 1> const& prinv, LINALG::Matrix<3, 1> const& dPI,
    LINALG::Matrix<6, 1> const& ddPII) const
{
  // according to Holzapfel-Nonlinear Solid Mechanics p. 216
  gamma(0) = 2. * (dPI(0) + prinv(0) * dPI(1));
  gamma(1) = -2. * dPI(1);
  gamma(2) = 2. * prinv(2) * dPI(2);

  // according to Holzapfel-Nonlinear Solid Mechanics p. 261
  delta(0) = 4. * (ddPII(0) + 2. * prinv(0) * ddPII(5) + dPI(1) + prinv(0) * prinv(0) * ddPII(1));
  delta(1) = -4. * (ddPII(5) + prinv(0) * ddPII(1));
  delta(2) = 4. * (prinv(2) * ddPII(4) + prinv(0) * prinv(2) * ddPII(3));
  delta(3) = 4. * ddPII(1);
  delta(4) = -4. * prinv(2) * ddPII(3);
  delta(5) = 4. * (prinv(2) * dPI(2) + prinv(2) * prinv(2) * ddPII(2));
  delta(6) = -4. * prinv(2) * dPI(2);
  delta(7) = -4. * dPI(1);

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::GrowthRemodel_ElastHyper::EvaluateIsotropicPrincElast(
    LINALG::Matrix<6, 1>& stressisoprinc, LINALG::Matrix<6, 6>& cmatisoprinc,
    LINALG::Matrix<6, 1> const& iCinv, LINALG::Matrix<6, 1> const& iCinCiCinv,
    LINALG::Matrix<6, 1> const& iCv, LINALG::Matrix<3, 1> const& gamma,
    LINALG::Matrix<8, 1> const& delta) const
{
  // 2nd Piola Kirchhoff stresses
  stressisoprinc.Update(gamma(0), iCinv, 1.0);
  stressisoprinc.Update(gamma(1), iCinCiCinv, 1.0);
  stressisoprinc.Update(gamma(2), iCv, 1.0);

  // constitutive tensor
  cmatisoprinc.MultiplyNT(delta(0), iCinv, iCinv, 1.);
  cmatisoprinc.MultiplyNT(delta(1), iCinCiCinv, iCinv, 1.);
  cmatisoprinc.MultiplyNT(delta(1), iCinv, iCinCiCinv, 1.);
  cmatisoprinc.MultiplyNT(delta(2), iCinv, iCv, 1.);
  cmatisoprinc.MultiplyNT(delta(2), iCv, iCinv, 1.);
  cmatisoprinc.MultiplyNT(delta(3), iCinCiCinv, iCinCiCinv, 1.);
  cmatisoprinc.MultiplyNT(delta(4), iCinCiCinv, iCv, 1.);
  cmatisoprinc.MultiplyNT(delta(4), iCv, iCinCiCinv, 1.);
  cmatisoprinc.MultiplyNT(delta(5), iCv, iCv, 1.);
  AddtoCmatHolzapfelProduct(cmatisoprinc, iCv, delta(6));
  AddtoCmatHolzapfelProduct(cmatisoprinc, iCinv, delta(7));

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::GrowthRemodel_ElastHyper::EvaluatedSdiFg(LINALG::Matrix<6, 9>& dSdiFg,
    LINALG::Matrix<3, 1> const& gamma, LINALG::Matrix<8, 1> const& delta,
    LINALG::Matrix<3, 3> const& iFinM, LINALG::Matrix<3, 3> const& iCinCM,
    LINALG::Matrix<6, 1> const& iCinv, LINALG::Matrix<9, 1> const& CiFin9x1,
    LINALG::Matrix<9, 1> const& CiFinCe9x1, LINALG::Matrix<6, 1> const& iCinCiCinv,
    LINALG::Matrix<9, 1> const& CiFiniCe9x1, LINALG::Matrix<6, 1> const& iCv,
    LINALG::Matrix<3, 3> const& iFinCeM, int const gp) const
{
  // clear some variables
  dSdiFg.Clear();

  static LINALG::Matrix<3, 3> id(true);
  for (int i = 0; i < 3; ++i) id(i, i) = 1.0;
  static LINALG::Matrix<6, 9> dSdiFin(true);
  dSdiFin.Clear();

  // derivative of second Piola Kirchhoff stress w.r.t. inverse growth deformation gradient
  MAT::AddRightNonSymmetricHolzapfelProduct(dSdiFin, id, iFinM, gamma(0));
  MAT::AddRightNonSymmetricHolzapfelProduct(dSdiFin, iCinCM, iFinM, gamma(1));
  dSdiFin.MultiplyNT(delta(0), iCinv, CiFin9x1, 1.);
  dSdiFin.MultiplyNT(delta(1), iCinv, CiFinCe9x1, 1.);
  dSdiFin.MultiplyNT(delta(1), iCinCiCinv, CiFin9x1, 1.);
  dSdiFin.MultiplyNT(delta(2), iCinv, CiFiniCe9x1, 1.);
  dSdiFin.MultiplyNT(delta(2), iCv, CiFin9x1, 1.);
  dSdiFin.MultiplyNT(delta(3), iCinCiCinv, CiFinCe9x1, 1.);
  dSdiFin.MultiplyNT(delta(4), iCinCiCinv, CiFiniCe9x1, 1.);
  dSdiFin.MultiplyNT(delta(4), iCv, CiFinCe9x1, 1.);
  dSdiFin.MultiplyNT(delta(5), iCv, CiFiniCe9x1, 1.);
  MAT::AddRightNonSymmetricHolzapfelProduct(dSdiFin, id, iFinCeM, 0.5 * delta(7));

  // diFin/diFg
  static LINALG::Matrix<9, 9> diFindiFg(true);
  diFindiFg.Clear();
  MAT::AddNonSymmetricProduct(1.0, id, GM_[gp], diFindiFg);

  dSdiFg.MultiplyNN(1.0, dSdiFin, diFindiFg, 0.0);

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::GrowthRemodel_ElastHyper::EvaluateAdditionalCmat(LINALG::Matrix<6, 6>& cmatadd,
    LINALG::Matrix<3, 3> const& diFgdrhoM, LINALG::Matrix<1, 6> const& sum_drhodC,
    LINALG::Matrix<6, 9> const& dSdiFg, int const gp) const
{
  // clear some variables
  cmatadd.Clear();

  // diFg/dC
  LINALG::Matrix<9, 6> diFgdC(true);
  LINALG::Matrix<9, 1> diFgdrho9x1(true);
  Matrix3x3to9x1(diFgdrhoM, diFgdrho9x1);
  diFgdC.MultiplyNN(1.0, diFgdrho9x1, sum_drhodC, 0.0);

  // update elasticity tensor
  cmatadd.MultiplyNN(2.0, dSdiFg, diFgdC, 0.0);

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::GrowthRemodel_ElastHyper::EvaluateStressCmatMembrane(LINALG::Matrix<3, 3> const& CM,
    LINALG::Matrix<3, 3> const& iFgM, LINALG::Matrix<6, 1>& stress, LINALG::Matrix<6, 6>& cmat,
    LINALG::Matrix<6, 9>& dSdiFg, const int gp, const int eleGID) const
{
  // clear some variables
  stress.Clear();
  cmat.Clear();
  dSdiFg.Clear();

  // 2nd Piola Kirchhoff stress
  LINALG::FADMatrix<3, 3> iFgM_fad(true);
  iFgM_fad = iFgM;
  iFgM_fad.diff(0, 9);
  static LINALG::FADMatrix<3, 3> CM_fad(true);
  CM_fad = CM;
  static LINALG::FADMatrix<3, 3> FgM_fad(true);
  FgM_fad.Invert(iFgM_fad);
  static LINALG::FADMatrix<3, 3> GM_fad(true);
  GM_fad = GM_[gp];
  static LINALG::FADMatrix<3, 3> iFinM_fad(true);
  iFinM_fad.MultiplyNN(1.0, iFgM_fad, GM_fad, 0.0);
  static LINALG::FADMatrix<3, 3> FinM_fad(true);
  FinM_fad.Invert(iFinM_fad);
  LINALG::FADMatrix<3, 3> CinM_fad(true);
  CinM_fad.MultiplyTN(1.0, FinM_fad, FinM_fad, 0.0);

  static LINALG::FADMatrix<3, 3> CeM_fad(true);
  static LINALG::FADMatrix<3, 3> tmp_fad(true);
  tmp_fad.MultiplyNN(1.0, CM_fad, iFinM_fad, 0.0);
  CeM_fad.MultiplyTN(1.0, iFinM_fad, tmp_fad, 0.0);

  static LINALG::FADMatrix<3, 3> id_fad(true);
  for (int i = 0; i < 3; ++i) id_fad(i, i) = 1.0;
  static LINALG::FADMatrix<3, 3> AradM_fad(true);
  AradM_fad = AradM_;
  static LINALG::FADMatrix<3, 3> AradgrM_fad(true);
  static LINALG::FADMatrix<3, 3> AorthgrM_fad(true);
  tmp_fad.MultiplyNN(1.0, FinM_fad, AradM_fad, 0.0);
  AradgrM_fad.MultiplyNT(1.0 / CinM_fad.Dot(AradM_fad), tmp_fad, FinM_fad, 0.0);
  AorthgrM_fad.Update(1.0, id_fad, 0.0);
  AorthgrM_fad.Update(-1.0, AradgrM_fad, 1.0);

  // X = Aorthgr*Ce*Aorthgr + Aradgr
  static LINALG::FADMatrix<3, 3> XM_fad(true);
  tmp_fad.MultiplyNN(1.0, AorthgrM_fad, CeM_fad, 0.0);
  XM_fad.MultiplyNN(1.0, tmp_fad, AorthgrM_fad, 0.0);
  XM_fad.Update(1.0, AradgrM_fad, 1.0);

  static LINALG::FADMatrix<3, 3> iFinAorthgriFinTM_fad(true);
  tmp_fad.MultiplyNN(1.0, iFinM_fad, AorthgrM_fad, 0.0);
  iFinAorthgriFinTM_fad.MultiplyNT(1.0, tmp_fad, iFinM_fad, 0.0);

  double mue_el_mem = 0.0;
  Teuchos::RCP<MAT::ELASTIC::IsoNeoHooke> matmem;
  for (unsigned k = 0; k < potsumelmem_.size(); ++k)
  {
    matmem = Teuchos::rcp_dynamic_cast<MAT::ELASTIC::IsoNeoHooke>(potsumelmem_[k]);
    mue_el_mem += matmem->Mue();
  }

  FAD X_det = XM_fad(0, 0) * (XM_fad(1, 1) * XM_fad(2, 2) - XM_fad(1, 2) * XM_fad(2, 1)) -
              XM_fad(0, 1) * (XM_fad(1, 0) * XM_fad(2, 2) - XM_fad(1, 2) * XM_fad(2, 0)) +
              XM_fad(0, 2) * (XM_fad(1, 0) * XM_fad(2, 1) - XM_fad(1, 1) * XM_fad(2, 0));

  static LINALG::FADMatrix<3, 3> iXM_fad(true);
  iXM_fad.Invert(XM_fad);
  static LINALG::FADMatrix<3, 3> ZM_fad(
      true);  // Z = F_{in}^{-T}*A_{gr}^{T}*X^{-1}*A_{gr}^{T}*F_{in}^{-1}
  static LINALG::FADMatrix<3, 3> ZTM_fad(true);
  static LINALG::FADMatrix<3, 3> AorthgriFinM_fad(true);
  AorthgriFinM_fad.MultiplyNN(1.0, AorthgrM_fad, iFinM_fad, 0.0);
  tmp_fad.MultiplyTT(1.0, AorthgriFinM_fad, iXM_fad, 0.0);
  ZTM_fad.MultiplyNN(1.0, tmp_fad, AorthgriFinM_fad, 0.0);
  ZM_fad.UpdateT(1.0, ZTM_fad, 0.0);

  static LINALG::FADMatrix<3, 3> stress_fad(true);
  stress_fad.Update(mue_el_mem * cur_rho_el_[gp] * mue_frac_[gp], iFinAorthgriFinTM_fad, 0.0);
  stress_fad.Update(-0.5 * mue_el_mem * cur_rho_el_[gp] * mue_frac_[gp] / X_det, ZM_fad, 1.0);
  stress_fad.Update(-0.5 * mue_el_mem * cur_rho_el_[gp] * mue_frac_[gp] / X_det, ZTM_fad, 1.0);

  static LINALG::Matrix<3, 3> stressM(true);
  stressM = stress_fad.ConverttoDouble();
  MatrixtoStressLikeVoigtNotation(stressM, stress);

  // Derivative of 2nd Piola Kirchhoff stress w.r.t. the inverse growth deformation gradient
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 9; ++j) dSdiFg(i, j) = stress_fad(i, i).dx(j);
  for (int j = 0; j < 9; ++j)
  {
    dSdiFg(3, j) = 0.5 * (stress_fad(0, 1).dx(j) + stress_fad(1, 0).dx(j));
    dSdiFg(4, j) = 0.5 * (stress_fad(1, 2).dx(j) + stress_fad(2, 1).dx(j));
    dSdiFg(5, j) = 0.5 * (stress_fad(0, 2).dx(j) + stress_fad(2, 0).dx(j));
  }

  // Elasticity tensor
  static LINALG::Matrix<3, 3> ZM(true);
  static LINALG::Matrix<3, 3> ZTM(true);
  static LINALG::Matrix<3, 3> XM(true);
  ZM = ZM_fad.ConverttoDouble();
  ZTM = ZTM_fad.ConverttoDouble();
  XM = XM_fad.ConverttoDouble();
  // Y = Z^T + Z
  static LINALG::Matrix<3, 3> YM(true);
  YM.Update(1.0, ZM, 0.0);
  YM.Update(1.0, ZTM, 1.0);
  static LINALG::Matrix<6, 1> Yv(true);
  MatrixtoStressLikeVoigtNotation(YM, Yv);
  static LINALG::Matrix<6, 6> dYdC(true);
  dYdC.Clear();
  MAT::AddtoCmatHolzapfelProduct(dYdC, Yv, -0.5);

  cmat.MultiplyNT(0.5 * mue_el_mem * cur_rho_el_[gp] * mue_frac_[gp] / X_det.val(), Yv, Yv, 0.0);
  cmat.Update(-mue_el_mem * cur_rho_el_[gp] * mue_frac_[gp] / X_det.val(), dYdC, 1.0);

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::GrowthRemodel_ElastHyper::EvaluateElastinDamage(double const& dt_pre)
{
  double T_el = 101 * 365.25;  // in days
  double t_dam = 40.0;

  for (unsigned gp = 0; gp < cur_rho_el_.size(); ++gp)
    cur_rho_el_[gp] = init_rho_el_[gp] * (exp(-(t_tot_ - dt_pre) / T_el) -
                                             0.5 * (1.0 - exp(-(t_tot_ - dt_pre) / t_dam)) *
                                                 exp(-0.5 * (100.0 * (fabs(gp_ax_[gp]) - 0.09)) *
                                                     (100.0 * (fabs(gp_ax_[gp]) - 0.09))));

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::GrowthRemodel_ElastHyper::EvaluateGrowthDefGrad(LINALG::Matrix<3, 3>& FgM,
    LINALG::Matrix<3, 3>& iFgM, LINALG::Matrix<3, 3>& dFgdrhoM, LINALG::Matrix<3, 3>& diFgdrhoM,
    const int gp)
{
  double rho_col_sum = 0.0;
  // evaluate volume change
  for (unsigned p = 0; p < potsumrf_.size(); ++p)
    for (unsigned k = 0; k < potsumrf_[p]->GetNumFibers(); ++k)
      rho_col_sum += potsumrf_[p]->GetCurMassDensity(k, gp);
  v_[gp] = (rho_col_sum + cur_rho_el_[gp]) / params_->density_;

  switch (params_->growthtype_)
  {
    case 1:
      // build growth and inverse growth deformation gradient (and its derivative) for anisotropic
      // growth
      FgM.Update(1.0, AplM_, 0.0);
      FgM.Update(v_[gp], AgM_, 1.0);
      iFgM.Update(1.0, AplM_, 0.0);
      iFgM.Update(1.0 / v_[gp], AgM_, 1.0);
      diFgdrhoM.Update(-(1. / (v_[gp] * v_[gp])) * (1. / params_->density_), AgM_, 0.0);
      dFgdrhoM.Update(1.0 / params_->density_, AgM_, 0.0);
      break;
    case 0:
      // build growth and inverse growth deformation gradient (and its derivative) for isotropic
      // growth
      FgM.Update(std::pow(v_[gp], 1. / 3.), AcirM_, 0.0);
      FgM.Update(std::pow(v_[gp], 1. / 3.), AradM_, 1.0);
      FgM.Update(std::pow(v_[gp], 1. / 3.), AaxM_, 1.0);
      iFgM.Update(std::pow(v_[gp], -1. / 3.), AcirM_, 0.0);
      iFgM.Update(std::pow(v_[gp], -1. / 3.), AradM_, 1.0);
      iFgM.Update(std::pow(v_[gp], -1. / 3.), AaxM_, 1.0);
      diFgdrhoM.Update(
          -1. / 3. * std::pow(v_[gp], -4. / 3.) * (1. / params_->density_), AcirM_, 0.0);
      diFgdrhoM.Update(
          -1. / 3. * std::pow(v_[gp], -4. / 3.) * (1. / params_->density_), AradM_, 1.0);
      diFgdrhoM.Update(
          -1. / 3. * std::pow(v_[gp], -4. / 3.) * (1. / params_->density_), AaxM_, 1.0);
      dFgdrhoM.Update(1. / 3. * std::pow(v_[gp], -2. / 3.) * (1. / params_->density_), AcirM_, 0.0);
      dFgdrhoM.Update(1. / 3. * std::pow(v_[gp], -2. / 3.) * (1. / params_->density_), AradM_, 1.0);
      dFgdrhoM.Update(1. / 3. * std::pow(v_[gp], -2. / 3.) * (1. / params_->density_), AaxM_, 1.0);
      break;
    default:
      dserror("growthtype has to be either 1: anisotropic growth or 0: isotropic growth");
      break;
  }

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::GrowthRemodel_ElastHyper::SetupGR2D(
    Teuchos::ParameterList& params, const double dt, const int gp)
{
  LINALG::Matrix<1, 3> gprefecoord(true);  // gp coordinates in reference configuration
  LINALG::Matrix<1, 3> axdir(true);
  gprefecoord = params.get<LINALG::Matrix<1, 3>>("gprefecoord");

  if ((params_->cylinder_ == 1 || params_->cylinder_ == 2 || params_->cylinder_ == 3) &&
      (setup_[0] == 1))
  {
    LINALG::Matrix<1, 3> elerefecoord(true);  // element center coordinates in reference system
    elerefecoord = params.get<LINALG::Matrix<1, 3>>("elecenter");
    SetupAxiCirRadCylinder(elerefecoord, dt);
    if (params_->growthtype_ == 1) SetupAnisoGrowthTensors();

    // Update fiber directions with new local coordinate system (radaxicirc_)
    for (unsigned k = 0; k < potsumrf_.size(); ++k) potsumrf_[k]->UpdateFiberDirs(radaxicirc_, dt);
  }

  for (int i = 0; i < 3; ++i) axdir(0, i) = radaxicirc_(i, 1);
  gp_ax_[gp] = axdir.Dot(gprefecoord);

  // update elastin prestretch in radial direction
  GM_[gp].Update(params_->lamb_prestretch_cir_, AcirM_, 0.0);
  GM_[gp].Update(params_->lamb_prestretch_ax_, AaxM_, 1.0);
  GM_[gp].Update(1. / (params_->lamb_prestretch_ax_ * params_->lamb_prestretch_cir_), AradM_, 1.0);

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::GrowthRemodel_ElastHyper::EvaluateMembrane(LINALG::Matrix<3, 3> const& defgrd_glob,
    double& rcg33, Teuchos::ParameterList& params, LINALG::Matrix<3, 3>& pk2M_glob,
    LINALG::Matrix<6, 6>& cmat_glob, const int eleGID)
{
  // blank resulting quantities
  pk2M_glob.Clear();
  cmat_glob.Clear();

  // save current simulation time (used for the evaluation of elastin degradation)
  t_tot_ = params.get<double>("total time");

  // current gauss point
  int gp = params.get<int>("gp");

  // time step size
  double dt = params.get<double>("delta time");

  if (setup_[gp] == 1)
  {
    SetupGR2D(params, dt, gp);
    setup_[gp] = 0;
  }

  // Evaluate growth deformation gradient
  // growth deformation gradient
  LINALG::Matrix<3, 3> FgM(true);
  LINALG::Matrix<3, 3> iFgM(true);
  LINALG::Matrix<3, 3> dFgdrhoM(true);
  LINALG::Matrix<3, 3> diFgdrhoM(true);
  EvaluateGrowthDefGrad(FgM, iFgM, dFgdrhoM, diFgdrhoM, gp);

  // Elastic deformation gradient
  LINALG::Matrix<3, 3> CM(true);
  CM.MultiplyTN(1.0, defgrd_glob, defgrd_glob, 0.0);
  LINALG::Matrix<3, 3> CeM(true);
  LINALG::Matrix<3, 3> iFinM(true);
  LINALG::Matrix<3, 3> tmp(true);
  iFinM.MultiplyNN(1.0, iFgM, GM_[gp], 0.0);
  tmp.MultiplyTN(1.0, iFinM, CM, 0.0);
  CeM.MultiplyNN(1.0, tmp, iFinM, 0.0);

  // Calculate adapted right Cauchy Green stretch in thickness direction
  LINALG::Matrix<3, 3> AplCeAplAradM(true);
  LINALG::Matrix<3, 3> CinM(true);
  tmp.Multiply(1.0, AplM_, CeM, 0.0);
  AplCeAplAradM.MultiplyNN(1.0, tmp, AplM_, 0.0);
  AplCeAplAradM.Update(1.0, AradM_, 1.0);
  tmp.MultiplyNT(1.0, iFinM, iFinM, 0.0);
  CinM.Invert(tmp);
  rcg33 = CinM.Dot(AradM_) / AplCeAplAradM.Determinant();

  // Update right Cauchy Green tensor
  double rcg33_old = CM.Dot(AradM_);
  CM.Update(-rcg33_old, AradM_, 1.0);
  CM.Update(rcg33, AradM_, 1.0);


  // Evaluate anisotropic remodel fibers
  static LINALG::Matrix<6, 1> pk2v_glob(true);
  pk2v_glob.Clear();
  static LINALG::Matrix<NUM_STRESS_3D, 1> stressaniso(true);
  static LINALG::Matrix<NUM_STRESS_3D, NUM_STRESS_3D> cmataniso(true);
  for (unsigned p = 0; p < potsumrf_.size(); ++p)
  {
    potsumrf_[p]->EvaluateAnisotropicStressCmat(CM, iFgM, cmataniso, stressaniso, gp, dt, eleGID);
    pk2v_glob.Update(1.0, stressaniso, 1.0);
    cmat_glob.Update(1.0, cmataniso, 1.0);
  }

  // Build stress response and elasticity tensor of membrane material
  static LINALG::Matrix<NUM_STRESS_3D, 1> stressmem(true);
  static LINALG::Matrix<NUM_STRESS_3D, NUM_STRESS_3D> cmatmem(true);
  static LINALG::Matrix<6, 9> dummy(true);
  EvaluateStressCmatMembrane(CM, iFgM, stressmem, cmatmem, dummy, gp, eleGID);
  pk2v_glob.Update(1.0, stressmem, 1.0);
  cmat_glob.Update(1.0, cmatmem, 1.0);
  StressLikeVoigtNotationtoMatrix(pk2v_glob, pk2M_glob);

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::GrowthRemodel_ElastHyper::VisNames(std::map<std::string, int>& names)
{
  std::string result_mass_fraction_el;
  result_mass_fraction_el = "mass_fraction_el";

  names[result_mass_fraction_el] = 1;


  std::string result_v_growth;
  result_v_growth = "v_growth";

  names[result_v_growth] = 1;


  // loop map of associated potential summands
  // remodelfiber
  for (unsigned int p = 0; p < potsumrf_.size(); ++p) potsumrf_[p]->VisNames(names, p);

  // 2D elastin matrix
  for (unsigned int p = 0; p < potsumelmem_.size(); ++p) potsumelmem_[p]->VisNames(names);

  if (params_->membrane_ != 1)
  {
    // 3D elastin matrix
    for (unsigned int p = 0; p < potsumeliso_.size(); ++p) potsumeliso_[p]->VisNames(names);

    // volpenalty
    potsumelpenalty_->VisNames(names);
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool MAT::GrowthRemodel_ElastHyper::VisData(
    const std::string& name, std::vector<double>& data, int numgp, int eleID)
{
  int return_val = 0;

  if (name == "mass_fraction_el")
  {
    if (data.size() != 1) dserror("size mismatch");
    for (unsigned i = 0; i < cur_rho_el_.size(); ++i)
    {
      data[0] += cur_rho_el_[i];
    }
    data[0] = data[0] / cur_rho_el_.size();

    return true;
  }

  if (name == "v_growth")
  {
    if (data.size() != 1) dserror("size mismatch");
    for (unsigned i = 0; i < v_.size(); ++i)
    {
      data[0] += v_[i];
    }
    data[0] = data[0] / v_.size();

    return true;
  }

  // loop map of associated potential summands
  // remodelfiber
  if (name.at(name.length() - 3) == '0')
    return_val += potsumrf_[0]->VisData(name, data, numgp, eleID);
  if (name.at(name.length() - 3) == '1')
    return_val += potsumrf_[1]->VisData(name, data, numgp, eleID);

  // 2D elastin matrix
  for (unsigned int p = 0; p < potsumelmem_.size(); ++p)
    return_val += potsumelmem_[p]->VisData(name, data, numgp, eleID);

  if (params_->membrane_ != 1)
  {
    // 3D elastin matrix
    for (unsigned int p = 0; p < potsumeliso_.size(); ++p)
      return_val += potsumeliso_[p]->VisData(name, data, numgp, eleID);

    // volpenalty
    return_val += potsumelpenalty_->VisData(name, data, numgp, eleID);
  }

  return (bool)return_val;
}
