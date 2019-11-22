/*----------------------------------------------------------------------*/
/*! \file
\brief
This file contains the viscohyperelastic material toolbox.
The available viscous models can be applied to any hyperelastic law
of the Elasthyper toolbox.
The viscous part is summed up with the hyperelastic laws
to build a viscohyperelastic material model.

The input line should read
MAT 0   MAT_ViscoElastHyper   NUMMAT 2 MATIDS 1 2 DENS 0

\level 1

\maintainer Amadeus Gebauer
 */

/*----------------------------------------------------------------------*/
#include "viscoelasthyper.H"
#include "../drt_lib/standardtypes_cpp.H"
#include "../drt_matelast/elast_summand.H"
#include "../linalg/linalg_utils.H"
#include "../drt_lib/drt_linedefinition.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_mat/matpar_bundle.H"
#include "../drt_mat/material_service.H"
#include "../drt_matelast/visco_generalizedgenmax.H"
#include "../drt_lib/voigt_notation.H"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::PAR::ViscoElastHyper::ViscoElastHyper(Teuchos::RCP<MAT::PAR::Material> matdata)
    : MAT::PAR::ElastHyper(matdata)
{
  // polyconvexity check is just implemented for isotropic hyperlastic materials
  if (polyconvex_)
    dserror(
        "This polyconvexity-check is just implemented for isotropic "
        "hyperelastic-materials (do not use for viscoelastic materials).");
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<MAT::Material> MAT::PAR::ViscoElastHyper::CreateMaterial()
{
  return Teuchos::rcp(new MAT::ViscoElastHyper(this));
}


MAT::ViscoElastHyperType MAT::ViscoElastHyperType::instance_;


DRT::ParObject* MAT::ViscoElastHyperType::Create(const std::vector<char>& data)
{
  MAT::ViscoElastHyper* elhy = new MAT::ViscoElastHyper();
  elhy->Unpack(data);

  return elhy;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::ViscoElastHyper::ViscoElastHyper() : MAT::ElastHyper()
{
  isovisco_ = false;
  viscogenmax_ = false;
  viscogeneralizedgenmax_ = false;
  viscofract_ = false;

  // history data 09/13
  isinitvis_ = false;
  histscgcurr_ = Teuchos::null;
  histscglast_ = Teuchos::null;
  histmodrcgcurr_ = Teuchos::null;
  histmodrcglast_ = Teuchos::null;
  histstresscurr_ = Teuchos::null;
  histstresslast_ = Teuchos::null;

  // genmax
  histartstresscurr_ = Teuchos::null;
  histartstresslast_ = Teuchos::null;

  // generalized genmax
  histbranchstresscurr_ = Teuchos::null;
  histbranchstresslast_ = Teuchos::null;
  histbranchelaststresscurr_ = Teuchos::null;
  histbranchelaststresslast_ = Teuchos::null;

  // fract
  histfractartstresscurr_ = Teuchos::null;
  histfractartstresslastall_ = Teuchos::null;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::ViscoElastHyper::ViscoElastHyper(MAT::PAR::ViscoElastHyper* params)
    : MAT::ElastHyper(params),
      isovisco_(false),
      viscogenmax_(false),
      viscogeneralizedgenmax_(false),
      viscofract_(false),
      isinitvis_(false)
{
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ViscoElastHyper::Pack(DRT::PackBuffer& data) const
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
  summandProperties_.Pack(data);
  AddtoPack(data, isovisco_);
  AddtoPack(data, viscogenmax_);
  AddtoPack(data, viscogeneralizedgenmax_);
  AddtoPack(data, viscofract_);

  if (params_ != NULL)  // summands are not accessible in postprocessing mode
  {
    // loop map of associated potential summands
    for (unsigned int p = 0; p < potsum_.size(); ++p)
    {
      potsum_[p]->PackSummand(data);
    }
  }

  //  pack history data 09/13
  int histsize;
  if (!Initialized())
  {
    histsize = 0;
  }
  else
  {
    histsize = histscglast_->size();
  }
  AddtoPack(data, histsize);  // Length of history vector(s)
  for (int var = 0; var < histsize; ++var)
  {
    AddtoPack(data, histscglast_->at(var));
    AddtoPack(data, histmodrcglast_->at(var));
    AddtoPack(data, histstresslast_->at(var));
    AddtoPack(data, histartstresslast_->at(var));
  }

  if (viscogeneralizedgenmax_)
  {
    for (int var = 0; var < histsize; ++var)
    {
      AddtoPack(data, histbranchstresslast_->at(var));
      AddtoPack(data, histbranchelaststresslast_->at(var));
    }
  }

  // pack history of FSLS-model
  if (viscofract_)
  {
    // check if history exists
    AddtoPack(data, (int)(histfractartstresslastall_ != Teuchos::null));
    if (!(int)(histfractartstresslastall_ != Teuchos::null))
      dserror("Something got wrong with your history data.");

    // pack stepsize
    AddtoPack(data, (int)histfractartstresslastall_->at(0).size());
    // pack history values
    for (int gp = 0; gp < (int)histfractartstresslastall_->size(); ++gp)
      for (int step = 0; step < (int)histfractartstresslastall_->at(gp).size(); ++step)
        AddtoPack(data, histfractartstresslastall_->at(gp).at(step));
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ViscoElastHyper::Unpack(const std::vector<char>& data)
{
  // make sure we have a pristine material
  params_ = NULL;
  potsum_.clear();

  summandProperties_.Clear();
  isovisco_ = false;
  viscogenmax_ = false;
  viscogeneralizedgenmax_ = false;
  viscofract_ = false;

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
        params_ = static_cast<MAT::PAR::ViscoElastHyper*>(mat);
      else
        dserror("Type of parameter material %d does not fit to calling type %d", mat->Type(),
            MaterialType());
    }
  }

  summandProperties_.Unpack(position, data);
  isovisco_ = (bool)ExtractInt(position, data);
  viscogenmax_ = (bool)ExtractInt(position, data);
  viscogeneralizedgenmax_ = (bool)ExtractInt(position, data);
  viscofract_ = (bool)ExtractInt(position, data);


  if (params_ != NULL)  // summands are not accessible in postprocessing mode
  {
    // make sure the referenced materials in material list have quick access parameters
    std::vector<int>::const_iterator m;
    for (m = params_->matids_->begin(); m != params_->matids_->end(); ++m)
    {
      const int matid = *m;
      Teuchos::RCP<MAT::ELASTIC::Summand> sum = MAT::ELASTIC::Summand::Factory(matid);
      if (sum == Teuchos::null) dserror("Failed to allocate");
      potsum_.push_back(sum);
    }

    // loop map of associated potential summands
    for (unsigned int p = 0; p < potsum_.size(); ++p)
    {
      potsum_[p]->UnpackSummand(data, position);
    }

    // history data 09/13
    isinitvis_ = true;
    int histsize;
    ExtractfromPack(position, data, histsize);

    if (histsize == 0) isinitvis_ = false;

    // initialize current variables
    histscgcurr_ = Teuchos::rcp(new std::vector<LINALG::Matrix<NUM_STRESS_3D, 1>>(histsize));
    histmodrcgcurr_ = Teuchos::rcp(new std::vector<LINALG::Matrix<NUM_STRESS_3D, 1>>(histsize));
    histstresscurr_ = Teuchos::rcp(new std::vector<LINALG::Matrix<NUM_STRESS_3D, 1>>(histsize));
    histartstresscurr_ = Teuchos::rcp(new std::vector<LINALG::Matrix<NUM_STRESS_3D, 1>>(histsize));

    // initialize last variables
    histscglast_ = Teuchos::rcp(new std::vector<LINALG::Matrix<NUM_STRESS_3D, 1>>(histsize));
    histmodrcglast_ = Teuchos::rcp(new std::vector<LINALG::Matrix<NUM_STRESS_3D, 1>>(histsize));
    histstresslast_ = Teuchos::rcp(new std::vector<LINALG::Matrix<NUM_STRESS_3D, 1>>(histsize));
    histartstresslast_ = Teuchos::rcp(new std::vector<LINALG::Matrix<NUM_STRESS_3D, 1>>(histsize));


    for (int gp = 0; gp < histsize; ++gp)
    {
      ExtractfromPack(position, data, histscglast_->at(gp));
      ExtractfromPack(position, data, histmodrcglast_->at(gp));
      ExtractfromPack(position, data, histstresslast_->at(gp));
      ExtractfromPack(position, data, histartstresslast_->at(gp));
    }

    if (viscogeneralizedgenmax_)
    {
      histbranchstresscurr_ =
          Teuchos::rcp(new std::vector<std::vector<LINALG::Matrix<NUM_STRESS_3D, 1>>>(histsize));
      histbranchelaststresscurr_ =
          Teuchos::rcp(new std::vector<std::vector<LINALG::Matrix<NUM_STRESS_3D, 1>>>(histsize));
      histbranchstresslast_ =
          Teuchos::rcp(new std::vector<std::vector<LINALG::Matrix<NUM_STRESS_3D, 1>>>(histsize));
      histbranchelaststresslast_ =
          Teuchos::rcp(new std::vector<std::vector<LINALG::Matrix<NUM_STRESS_3D, 1>>>(histsize));

      for (int gp = 0; gp < histsize; ++gp)
      {
        ExtractfromPack(position, data, histbranchstresslast_->at(gp));
        ExtractfromPack(position, data, histbranchelaststresslast_->at(gp));
      }
    }

    // for FSLS-model
    if (viscofract_)
    {
      // check if history data is saved
      bool have_historyalldata = (bool)ExtractInt(position, data);
      if (!have_historyalldata) dserror("Something got wrong with your history data.");

      histfractartstresscurr_ =
          Teuchos::rcp(new std::vector<LINALG::Matrix<NUM_STRESS_3D, 1>>(histsize));

      int histfractartstressall_stepsize = ExtractInt(position, data);
      histfractartstresslastall_ = Teuchos::rcp(new std::vector<std::vector<LINALG::Matrix<6, 1>>>(
          histsize, std::vector<LINALG::Matrix<6, 1>>(histfractartstressall_stepsize)));
      for (int gp = 0; gp < histsize; ++gp)
        for (int step = 0; step < histfractartstressall_stepsize; ++step)
          ExtractfromPack(position, data, histfractartstresslastall_->at(gp).at(step));
    }

    // in the postprocessing mode, we do not unpack everything we have packed
    // -> position check cannot be done in this case
    if (position != data.size())
      dserror("Mismatch in size of data %d <-> %d", data.size(), position);
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ViscoElastHyper::Setup(int numgp, DRT::INPUT::LineDefinition* linedef)
{
  // Setup summands
  for (unsigned int p = 0; p < potsum_.size(); ++p) potsum_[p]->Setup(linedef);

  // find out which formulations are used
  isovisco_ = false;
  viscogenmax_ = false;
  viscogeneralizedgenmax_ = false;
  viscofract_ = false;

  summandProperties_.Clear();
  ElastHyperProperties(potsum_, summandProperties_);


  if (summandProperties_.viscoGeneral)
  {
    for (unsigned int p = 0; p < potsum_.size(); ++p)
    {
      potsum_[p]->SpecifyViscoFormulation(
          isovisco_, viscogenmax_, viscogeneralizedgenmax_, viscofract_);
    }
  }

  // Initialise/allocate history variables 09/13
  const LINALG::Matrix<6, 1> emptyvec(true);
  LINALG::Matrix<6, 1> idvec(true);
  for (int i = 0; i < 3; ++i) idvec(i) = 1.;

  histscgcurr_ = Teuchos::rcp(new std::vector<LINALG::Matrix<6, 1>>(numgp, idvec));
  histscglast_ = Teuchos::rcp(new std::vector<LINALG::Matrix<6, 1>>(numgp, idvec));
  histmodrcgcurr_ = Teuchos::rcp(new std::vector<LINALG::Matrix<6, 1>>(numgp, idvec));
  histmodrcglast_ = Teuchos::rcp(new std::vector<LINALG::Matrix<6, 1>>(numgp, idvec));
  histstresscurr_ =
      Teuchos::rcp(new std::vector<LINALG::Matrix<NUM_STRESS_3D, 1>>(numgp, emptyvec));
  histstresslast_ =
      Teuchos::rcp(new std::vector<LINALG::Matrix<NUM_STRESS_3D, 1>>(numgp, emptyvec));
  histartstresscurr_ =
      Teuchos::rcp(new std::vector<LINALG::Matrix<NUM_STRESS_3D, 1>>(numgp, emptyvec));
  histartstresslast_ =
      Teuchos::rcp(new std::vector<LINALG::Matrix<NUM_STRESS_3D, 1>>(numgp, emptyvec));

  if (viscogeneralizedgenmax_)
  {
    const std::vector<LINALG::Matrix<6, 1>> emptybigvec(true);
    histbranchstresscurr_ = Teuchos::rcp(
        new std::vector<std::vector<LINALG::Matrix<NUM_STRESS_3D, 1>>>(numgp, emptybigvec));
    histbranchstresslast_ = Teuchos::rcp(
        new std::vector<std::vector<LINALG::Matrix<NUM_STRESS_3D, 1>>>(numgp, emptybigvec));
    histbranchelaststresscurr_ = Teuchos::rcp(
        new std::vector<std::vector<LINALG::Matrix<NUM_STRESS_3D, 1>>>(numgp, emptybigvec));
    histbranchelaststresslast_ = Teuchos::rcp(
        new std::vector<std::vector<LINALG::Matrix<NUM_STRESS_3D, 1>>>(numgp, emptybigvec));
  }

  // in case of FSLS-model
  if (viscofract_)
  {
    histfractartstresscurr_ =
        Teuchos::rcp(new std::vector<LINALG::Matrix<NUM_STRESS_3D, 1>>(numgp, emptyvec));
    // set true that history size is known
    histfractartstresslastall_ =
        Teuchos::rcp(new std::vector<std::vector<LINALG::Matrix<NUM_STRESS_3D, 1>>>(
            numgp, std::vector<LINALG::Matrix<NUM_STRESS_3D, 1>>(true)));
  }

  isinitvis_ = true;

  return;
}

/*------------------------------------------------------------------------------------------*
|  Setup internal stress variables - special for the inverse analysis (public)         09/13|
 *-------------------------------------------------------------------------------------------*/
void MAT::ViscoElastHyper::ResetAll(const int numgp)
{
  // Initialise/allocate history variables 09/13
  const LINALG::Matrix<6, 1> emptyvec(true);
  LINALG::Matrix<6, 1> idvec(true);
  for (int i = 0; i < 3; ++i) idvec(i) = 1.;

  histscgcurr_ = Teuchos::rcp(new std::vector<LINALG::Matrix<6, 1>>(numgp, idvec));
  histscglast_ = Teuchos::rcp(new std::vector<LINALG::Matrix<6, 1>>(numgp, idvec));
  histmodrcgcurr_ = Teuchos::rcp(new std::vector<LINALG::Matrix<6, 1>>(numgp, idvec));
  histmodrcglast_ = Teuchos::rcp(new std::vector<LINALG::Matrix<6, 1>>(numgp, idvec));
  histstresscurr_ =
      Teuchos::rcp(new std::vector<LINALG::Matrix<NUM_STRESS_3D, 1>>(numgp, emptyvec));
  histstresslast_ =
      Teuchos::rcp(new std::vector<LINALG::Matrix<NUM_STRESS_3D, 1>>(numgp, emptyvec));
  histartstresscurr_ =
      Teuchos::rcp(new std::vector<LINALG::Matrix<NUM_STRESS_3D, 1>>(numgp, emptyvec));
  histartstresslast_ =
      Teuchos::rcp(new std::vector<LINALG::Matrix<NUM_STRESS_3D, 1>>(numgp, emptyvec));


  if (viscogeneralizedgenmax_)
  {
    const std::vector<LINALG::Matrix<6, 1>> emptybigvec(true);
    histbranchstresscurr_ = Teuchos::rcp(
        new std::vector<std::vector<LINALG::Matrix<NUM_STRESS_3D, 1>>>(numgp, emptybigvec));
    histbranchstresslast_ = Teuchos::rcp(
        new std::vector<std::vector<LINALG::Matrix<NUM_STRESS_3D, 1>>>(numgp, emptybigvec));
    histbranchelaststresscurr_ = Teuchos::rcp(
        new std::vector<std::vector<LINALG::Matrix<NUM_STRESS_3D, 1>>>(numgp, emptybigvec));
    histbranchelaststresslast_ = Teuchos::rcp(
        new std::vector<std::vector<LINALG::Matrix<NUM_STRESS_3D, 1>>>(numgp, emptybigvec));
  }

  // in case of FSLS-model
  if (viscofract_)
  {
    histfractartstresscurr_ =
        Teuchos::rcp(new std::vector<LINALG::Matrix<NUM_STRESS_3D, 1>>(numgp, emptyvec));
    histfractartstresslastall_ =
        Teuchos::rcp(new std::vector<std::vector<LINALG::Matrix<NUM_STRESS_3D, 1>>>(
            numgp, std::vector<LINALG::Matrix<NUM_STRESS_3D, 1>>(true)));
  }

  isinitvis_ = true;

  return;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void MAT::ViscoElastHyper::Update()
{
  MAT::ElastHyper::Update();

  // Update history values 09/13
  histscglast_ = histscgcurr_;
  histmodrcglast_ = histmodrcgcurr_;
  histstresslast_ = histstresscurr_;
  histartstresslast_ = histartstresscurr_;

  // for FSLS-model
  if (viscofract_)
  {
    // To calculate the fractional derivative the history of all previous timesteps is saved in each
    // gauß-point

    // numsteps
    const Teuchos::ParameterList& sdyn = DRT::Problem::Instance()->StructuralDynamicParams();
    const int numsteps = sdyn.get<int>("NUMSTEP");
    // maximal size of history (in time steps)
    const unsigned int max_hist = numsteps + 1;

    // for each gp
    for (int gp = 0; gp < (int)histfractartstresslastall_->size(); ++gp)
    {
      // add current stress (all gp) at the end of history data
      histfractartstresslastall_->at(gp).push_back(histfractartstresscurr_->at(gp));

      // this is in the moment not used, because all history data is saved
      // if maximal size of history-vector is reached
      if (histfractartstresslastall_->at(gp).size() > max_hist)
      {
        // Hint: If you want to use this functionaltiy, reimplement it with a pointer (abirzle
        // 12/17)

        // save current history data in tmp_vec, but skip the first data
        // (= oldest data, which should deleted)
        std::vector<LINALG::Matrix<6, 1>> tmp_vec(
            ++histfractartstresslastall_->at(gp).begin(), histfractartstresslastall_->at(gp).end());
        // save data back to history-vector
        histfractartstresslastall_->at(gp) = tmp_vec;
      }
    }
  }

  // initialize current data
  const LINALG::Matrix<6, 1> emptyvec(true);

  LINALG::Matrix<6, 1> idvec(true);
  for (int i = 0; i < 3; ++i) idvec(i) = 1.;
  const int numgp = histscglast_->size();

  histscgcurr_ = Teuchos::rcp(new std::vector<LINALG::Matrix<6, 1>>(numgp, idvec));
  histmodrcgcurr_ = Teuchos::rcp(new std::vector<LINALG::Matrix<6, 1>>(numgp, idvec));
  histstresscurr_ =
      Teuchos::rcp(new std::vector<LINALG::Matrix<NUM_STRESS_3D, 1>>(numgp, emptyvec));
  histartstresscurr_ =
      Teuchos::rcp(new std::vector<LINALG::Matrix<NUM_STRESS_3D, 1>>(numgp, emptyvec));
  histfractartstresscurr_ =
      Teuchos::rcp(new std::vector<LINALG::Matrix<NUM_STRESS_3D, 1>>(numgp, emptyvec));

  if (viscogeneralizedgenmax_)
  {
    histbranchstresslast_ = histbranchstresscurr_;
    histbranchelaststresslast_ = histbranchelaststresscurr_;

    const std::vector<LINALG::Matrix<6, 1>> emptybigvec(true);
    histbranchstresscurr_ = Teuchos::rcp(
        new std::vector<std::vector<LINALG::Matrix<NUM_STRESS_3D, 1>>>(numgp, emptybigvec));
    histbranchelaststresscurr_ = Teuchos::rcp(
        new std::vector<std::vector<LINALG::Matrix<NUM_STRESS_3D, 1>>>(numgp, emptybigvec));
  }

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ViscoElastHyper::Evaluate(const LINALG::Matrix<3, 3>* defgrd,
    const LINALG::Matrix<6, 1>* glstrain, Teuchos::ParameterList& params,
    LINALG::Matrix<6, 1>* stress, LINALG::Matrix<6, 6>* cmat, const int eleGID)
{
  LINALG::Matrix<6, 1> C_strain(true);
  LINALG::Matrix<6, 1> C_stress(true);
  LINALG::Matrix<6, 1> iC_strain(true);
  LINALG::Matrix<6, 1> iC_stress(true);
  LINALG::Matrix<6, 1> modC_strain(true);
  LINALG::Matrix<6, 1> id2(true);
  LINALG::Matrix<6, 1> modrcg(true);
  LINALG::Matrix<6, 6> id4(true);
  LINALG::Matrix<6, 6> id4sharp(true);
  LINALG::Matrix<3, 1> prinv(true);
  LINALG::Matrix<3, 1> modinv(true);
  LINALG::Matrix<7, 1> rateinv(true);
  LINALG::Matrix<7, 1> modrateinv(true);

  LINALG::Matrix<3, 1> dPI(true);
  LINALG::Matrix<6, 1> ddPII(true);

  LINALG::Matrix<6, 1> scgrate(true);
  LINALG::Matrix<6, 1> modrcgrate(true);

  // for extension: LINALG::Matrix<6,1> modicgrate(true);
  LINALG::Matrix<8, 1> mu(true);
  LINALG::Matrix<8, 1> modmu(true);
  LINALG::Matrix<33, 1> xi(true);
  LINALG::Matrix<33, 1> modxi(true);

  EvaluateRightCauchyGreenStrainLikeVoigt(*glstrain, C_strain);
  VStrainUtils::InverseTensor(C_strain, iC_strain);
  VStrainUtils::ToStressLike(iC_strain, iC_stress);
  VStrainUtils::ToStressLike(C_strain, C_stress);
  VStrainUtils::InvariantsPrincipal(prinv, C_strain);


  UTILS::VOIGT::IdentityMatrix(id2);

  using VoigtNotation = UTILS::VOIGT::NotationType;
  UTILS::VOIGT::FourthOrderIdentityMatrix<VoigtNotation::stress, VoigtNotation::stress>(id4sharp);
  UTILS::VOIGT::FourthOrderIdentityMatrix<VoigtNotation::stress, VoigtNotation::strain>(id4);

  ElastHyperEvaluateInvariantDerivatives(prinv, dPI, ddPII, potsum_, summandProperties_, eleGID);

  if (isovisco_)
  {
    if (summandProperties_.isomod)
    {
      // calculate modified invariants
      InvariantsModified(modinv, prinv);
    }
    // calculate viscous quantities
    EvaluateKinQuantVis(C_strain, C_stress, iC_stress, prinv, rateinv, modC_strain, params, scgrate,
        modrcgrate, modrateinv);
    EvaluateMuXi(prinv, modinv, mu, modmu, xi, modxi, rateinv, modrateinv, params, eleGID);
  }

  // blank resulting quantities
  // ... even if it is an implicit law that cmat is zero upon input
  stress->Clear();
  cmat->Clear();

  // add isotropic part
  ElastHyperAddIsotropicStressCmat(*stress, *cmat, C_strain, iC_strain, prinv, dPI, ddPII);


  // add viscous part
  if (isovisco_)
  {
    if (summandProperties_.isomod)
    {
      // add viscous part decoupled
      LINALG::Matrix<NUM_STRESS_3D, 1> stressisomodisovisco(true);
      LINALG::Matrix<NUM_STRESS_3D, NUM_STRESS_3D> cmatisomodisovisco(true);
      LINALG::Matrix<NUM_STRESS_3D, 1> stressisomodvolvisco(true);
      LINALG::Matrix<NUM_STRESS_3D, NUM_STRESS_3D> cmatisomodvolvisco(true);
      EvaluateIsoViscoModified(stressisomodisovisco, stressisomodvolvisco, cmatisomodisovisco,
          cmatisomodvolvisco, prinv, modinv, modmu, modxi, C_strain, id2, iC_stress, id4,
          modrcgrate);
      stress->Update(1.0, stressisomodisovisco, 1.0);
      stress->Update(1.0, stressisomodvolvisco, 1.0);
      cmat->Update(1.0, cmatisomodisovisco, 1.0);
      cmat->Update(1.0, cmatisomodvolvisco, 1.0);
    }

    if (summandProperties_.isoprinc)
    {
      // add viscous part coupled
      LINALG::Matrix<NUM_STRESS_3D, 1> stressisovisco(true);
      LINALG::Matrix<NUM_STRESS_3D, NUM_STRESS_3D> cmatisovisco(true);
      EvaluateIsoViscoPrincipal(stressisovisco, cmatisovisco, mu, xi, id4sharp, scgrate);
      stress->Update(1.0, stressisovisco, 1.0);
      cmat->Update(1.0, cmatisovisco, 1.0);
    }
  }

  // add contribution of viscogenmax-material
  if (viscogenmax_)
  {
    LINALG::Matrix<NUM_STRESS_3D, 1> Q(true);  // artificial viscous stress
    LINALG::Matrix<NUM_STRESS_3D, NUM_STRESS_3D> cmatq(true);
    EvaluateViscoGenMax(stress, cmat, Q, cmatq, params);
    stress->Update(1.0, Q, 1.0);
    cmat->Update(1.0, cmatq, 1.0);
  }

  // add contribution of generalized Maxwell model
  if (viscogeneralizedgenmax_)
  {
    LINALG::Matrix<NUM_STRESS_3D, 1> Q(true);
    LINALG::Matrix<NUM_STRESS_3D, NUM_STRESS_3D> cmatq(true);
    EvaluateViscoGeneralizedGenMax(Q, cmatq, params, glstrain, eleGID);
    stress->Update(1.0, Q, 1.0);
    cmat->Update(1.0, cmatq, 1.0);
  }

  // add contribution of viscofract-material
  if (viscofract_)
  {
    LINALG::Matrix<NUM_STRESS_3D, 1> Q(true);  // artificial viscous stress
    LINALG::Matrix<NUM_STRESS_3D, NUM_STRESS_3D> cmatq(true);
    EvaluateViscoFract(*stress, *cmat, Q, cmatq, params);
    stress->Update(1.0, Q, 1.);
    cmat->Update(1.0, cmatq, 1.);
  }


  /*----------------------------------------------------------------------*/
  // coefficients in principal stretches
  if (summandProperties_.coeffStretchesPrinc || summandProperties_.coeffStretchesMod)
  {
    ElastHyperAddResponseStretches(*cmat, *stress, C_strain, potsum_, summandProperties_, eleGID);
  }

  /*----------------------------------------------------------------------*/
  // Do all the anisotropic stuff!
  if (summandProperties_.anisoprinc)
  {
    ElastHyperAddAnisotropicPrinc(*stress, *cmat, C_strain, params, eleGID, potsum_);
  }

  if (summandProperties_.anisomod)
  {
    ElastHyperAddAnisotropicMod(*stress, *cmat, C_strain, iC_strain, prinv, eleGID, potsum_);
  }
  return;
}

/*----------------------------------------------------------------------*/
/* Evaluate Quantities for viscous Part                           09/13 */
/*----------------------------------------------------------------------*/
void MAT::ViscoElastHyper::EvaluateKinQuantVis(LINALG::Matrix<6, 1>& rcg, LINALG::Matrix<6, 1>& scg,
    LINALG::Matrix<6, 1>& icg, LINALG::Matrix<3, 1>& prinv, LINALG::Matrix<7, 1>& rateinv,
    LINALG::Matrix<6, 1>& modrcg, Teuchos::ParameterList& params, LINALG::Matrix<6, 1>& scgrate,
    LINALG::Matrix<6, 1>& modrcgrate, LINALG::Matrix<7, 1>& modrateinv)
{
  // time derivative
  // -------------------------------------------------------------------
  // get gauss point number of this element
  const int gp = params.get<int>("gp", -1);

  // get time algorithmic parameters
  double dt = params.get<double>("delta time");

  // modrcg : \overline{C} = J^{-\frac{2}{3}} C
  const double modscale = std::pow(prinv(2), -1. / 3.);
  modrcg.Update(modscale, rcg);

  // read history
  LINALG::Matrix<6, 1> scglast(histscglast_->at(gp));
  LINALG::Matrix<6, 1> modrcglast(histmodrcglast_->at(gp));

  // Update history of Cauchy-Green Tensor
  histscgcurr_->at(gp) = scg;        // principal material: store C^{n}
  histmodrcgcurr_->at(gp) = modrcg;  // decoupled material: store \overline{C}^{n}

  // rate of Cauchy-Green Tensor
  // REMARK: strain-like 6-Voigt vector
  scgrate.Update(1.0, scg, 1.0);  // principal material: \dot{C} = \frac{C^n - C^{n-1}}{\Delta t}
  scgrate.Update(-1.0, scglast, 1.0);
  scgrate.Scale(1 / dt);

  modrcgrate.Update(1.0, modrcg, 1.0);  // decoupled material: \overline{\dot{C}} =
                                        // \frac{\overline{C}^n - \overline{C}^{n-1}}{\Delta t}
  modrcgrate.Update(-1.0, modrcglast, 1.0);
  modrcgrate.Scale(1 / dt);

  // invariants
  // -------------------------------------------------------------------
  // Second Invariant of modrcgrate \bar{J}_2 = \frac{1}{2} \tr (\dot{\overline{C^2}}
  modrateinv(1) =
      0.5 * (modrcgrate(0) * modrcgrate(0) + modrcgrate(1) * modrcgrate(1) +
                modrcgrate(2) * modrcgrate(2) + .5 * modrcgrate(3) * modrcgrate(3) +
                .5 * modrcgrate(4) * modrcgrate(4) + .5 * modrcgrate(5) * modrcgrate(5));


  // For further extension of material law (not necassary at the moment)
  /*
  // necassary transfer variable: LINALG::Matrix<6,1>& modicgrate
  // \overline{J}_3 = determinant of modified rate of right Cauchy-Green-Tensor
  modrateinv(2) = modrcgrate(0)*modrcgrate(1)*modrcgrate(2)
      + 0.25 * modrcgrate(3)*modrcgrate(4)*modrcgrate(5)
      - 0.25 * modrcgrate(1)*modrcgrate(5)*modrcgrate(5)
      - 0.25 * modrcgrate(2)*modrcgrate(3)*modrcgrate(3)
      - 0.25 * modrcgrate(0)*modrcgrate(4)*modrcgrate(4);

  // invert modified rate of right Cauchy-Green tensor
  // REMARK: stress-like 6-Voigt vector
  {
    modicgrate(0) = ( modrcgrate(1)*modrcgrate(2) - 0.25*modrcgrate(4)*modrcgrate(4) ) /
  modrateinv(2); modicgrate(1) = ( modrcgrate(0)*modrcgrate(2) - 0.25*modrcgrate(5)*modrcgrate(5) )
  / modrateinv(2); modicgrate(2) = ( modrcgrate(0)*modrcgrate(1) - 0.25*modrcgrate(3)*modrcgrate(3)
  ) / modrateinv(2); modicgrate(3) = ( 0.25*modrcgrate(5)*modrcgrate(4) -
  0.5*modrcgrate(3)*modrcgrate(2) ) / modrateinv(2); modicgrate(4) = (
  0.25*modrcgrate(3)*modrcgrate(5) - 0.5*modrcgrate(0)*modrcgrate(4) ) / modrateinv(2);
    modicgrate(5) = ( 0.25*modrcgrate(3)*modrcgrate(4) - 0.5*modrcgrate(5)*modrcgrate(1) ) /
  modrateinv(2);
  }
   */
}

/*----------------------------------------------------------------------*/
/* Evaluate Factors for viscous Quantities                        09/13 */
/*----------------------------------------------------------------------*/
void MAT::ViscoElastHyper::EvaluateMuXi(LINALG::Matrix<3, 1>& prinv, LINALG::Matrix<3, 1>& modinv,
    LINALG::Matrix<8, 1>& mu, LINALG::Matrix<8, 1>& modmu, LINALG::Matrix<33, 1>& xi,
    LINALG::Matrix<33, 1>& modxi, LINALG::Matrix<7, 1>& rateinv, LINALG::Matrix<7, 1>& modrateinv,
    Teuchos::ParameterList& params, const int eleGID)
{
  // principal materials
  if (summandProperties_.isoprinc)
  {
    // loop map of associated potential summands
    for (unsigned int p = 0; p < potsum_.size(); ++p)
    {
      potsum_[p]->AddCoefficientsViscoPrincipal(prinv, mu, xi, rateinv, params, eleGID);
    }
  }

  // decoupled (volumetric or isochoric) materials
  if (summandProperties_.isomod)
  {
    // loop map of associated potential summands
    for (unsigned int p = 0; p < potsum_.size(); ++p)
    {
      potsum_[p]->AddCoefficientsViscoModified(modinv, modmu, modxi, modrateinv, params, eleGID);
    }
  }
}

/*----------------------------------------------------------------------*/
/* stress and constitutive tensor for principal viscous materials       */
/*                                                        pfaller May15 */
/*----------------------------------------------------------------------*/
void MAT::ViscoElastHyper::EvaluateIsoViscoPrincipal(LINALG::Matrix<6, 1>& stress,
    LINALG::Matrix<6, 6>& cmat, LINALG::Matrix<8, 1>& mu, LINALG::Matrix<33, 1>& xi,
    LINALG::Matrix<6, 6>& id4sharp, LINALG::Matrix<6, 1>& scgrate)
{
  // contribution: \dot{C}
  stress.Update(mu(2), scgrate, 1.0);

  // contribution: id4sharp_{ijkl} = 1/2 (\delta_{ik}\delta_{jl} + \delta_{il}\delta_{jk})
  cmat.Update(xi(2), id4sharp, 1.0);

  return;
}

/*----------------------------------------------------------------------*/
/* stress and constitutive tensor for decoupled viscous materials 09/13 */
/*----------------------------------------------------------------------*/
void MAT::ViscoElastHyper::EvaluateIsoViscoModified(LINALG::Matrix<6, 1>& stressisomodisovisco,
    LINALG::Matrix<6, 1>& stressisomodvolvisco, LINALG::Matrix<6, 6>& cmatisomodisovisco,
    LINALG::Matrix<6, 6>& cmatisomodvolvisco, LINALG::Matrix<3, 1>& prinv,
    LINALG::Matrix<3, 1>& modinv, LINALG::Matrix<8, 1>& modmu, LINALG::Matrix<33, 1>& modxi,
    LINALG::Matrix<6, 1>& rcg, LINALG::Matrix<6, 1>& id2, LINALG::Matrix<6, 1>& icg,
    LINALG::Matrix<6, 6>& id4, LINALG::Matrix<6, 1>& modrcgrate)
{
  // define necessary variables
  const double modscale = std::pow(prinv(2), -1. / 3.);

  // 2nd Piola Kirchhoff stresses

  // isochoric contribution
  LINALG::Matrix<6, 1> modstress(true);
  modstress.Update(modmu(1), id2);
  modstress.Update(modmu(2), modrcgrate, 1.0);
  // build 4-tensor for projection as 6x6 tensor
  LINALG::Matrix<6, 6> Projection;
  Projection.MultiplyNT(1. / 3., icg, rcg);
  Projection.Update(1.0, id4, -1.0);
  // isochoric stress
  stressisomodisovisco.MultiplyNN(modscale, Projection, modstress, 1.0);

  // volumetric contribution:
  // with visco_isoratedep: no volumetric part added --> always 0


  // Constitutive Tensor

  // isochoric contribution
  // modified constitutive tensor
  LINALG::Matrix<6, 6> modcmat(true);
  LINALG::Matrix<6, 6> modcmat2(true);
  // contribution:  Id \otimes \overline{\dot{C}} + \overline{\dot{C}} \otimes Id
  modcmat.MultiplyNT(modxi(1), id2, modrcgrate);
  modcmat.MultiplyNT(modxi(1), modrcgrate, id2, 1.0);
  // contribution: Id4
  modcmat.Update(modxi(2), id4, 1.0);
  // scaling
  modcmat.Scale(std::pow(modinv(2), -4. / 3.));
  // contribution: P:\overline{C}:P
  modcmat2.MultiplyNN(Projection, modcmat);
  cmatisomodisovisco.MultiplyNT(1.0, modcmat2, Projection, 1.0);
  // contribution: 2/3*Tr(J^(-2/3)modstress) (Cinv \odot Cinv - 1/3 Cinv \otimes Cinv)
  modcmat.Clear();
  modcmat.MultiplyNT(-1.0 / 3.0, icg, icg);
  AddtoCmatHolzapfelProduct(modcmat, icg, 1.0);
  LINALG::Matrix<1, 1> tracemat;
  tracemat.MultiplyTN(2. / 3. * std::pow(modinv(2), -2. / 3.), modstress, rcg);
  cmatisomodisovisco.Update(tracemat(0, 0), modcmat, 1.0);
  // contribution: -2/3 (Cinv \otimes S_iso^v + S_iso^v \otimes Cinv)
  cmatisomodisovisco.MultiplyNT(-2. / 3., icg, stressisomodisovisco, 1.0);
  cmatisomodisovisco.MultiplyNT(-2. / 3., stressisomodisovisco, icg, 1.0);

  // volumetric contribution:
  // with visco_isoratedep: no volumetric part added --> always 0

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ViscoElastHyper::EvaluateViscoGenMax(LINALG::Matrix<6, 1>* stress,
    LINALG::Matrix<6, 6>* cmat, LINALG::Matrix<6, 1>& Q, LINALG::Matrix<6, 6>& cmatq,
    Teuchos::ParameterList& params)
{
  // initialize material parameters
  double tau = -1.0;
  double beta = -1.0;
  std::string solve = "";
  // alpha not used in viscogenmax (viscofract with alpha=1 results in viscogenmax-material)
  double alpha(true);

  // read material parameters of viscogenmax material
  for (unsigned int p = 0; p < potsum_.size(); ++p)
    potsum_[p]->ReadMaterialParametersVisco(tau, beta, alpha, solve);

  if (solve == "OST")
  {
    // initialize scalars
    double lambdascalar1(true);
    double lambdascalar2(true);
    double deltascalar(true);
    double theta = 0.5;

    // get theta of global time integration scheme to use it here
    // if global time integration scheme is not ONESTEPTHETA, theta is by default = 0.5 (abirzle
    // 09/14)
    std::string dyntype =
        DRT::Problem::Instance()->StructuralDynamicParams().get<std::string>("DYNAMICTYP");
    if (dyntype == "OneStepTheta")
      theta = DRT::Problem::Instance()
                  ->StructuralDynamicParams()
                  .sublist("ONESTEPTHETA")
                  .get<double>("THETA");

    // get time algorithmic parameters
    // NOTE: dt can be zero (in restart of STI) for Generalized Maxwell model
    // there is no special treatment required. Adaptation for Kelvin-Voigt were necessary.
    double dt = params.get<double>("delta time");  // TIMESTEP in the .dat file

    // evaluate scalars to compute
    // Q^(n+1) = tau/(tau+theta*dt) [(tau-dt+theta*dt)/tau Q + beta(S^(n+1) - S^n)]
    lambdascalar1 = tau / (tau + theta * dt);
    lambdascalar2 = (tau - dt + theta * dt) / tau;

    // factor to calculate visco stiffness matrix from elastic stiffness matrix
    // old Version: scalarvisco =
    // 1+beta_isoprinc*exp(-dt/(2*tau_isoprinc));//+alpha1*tau/(tau+theta*dt); Alines' version:
    // scalarvisco = beta_isoprinc*exp(-dt/(2*tau_isoprinc));//+alpha1*tau/(tau+theta*dt); Scalar
    // consistent to derivation of Q with one-step-theta-schema (abirzle 09/14):
    deltascalar = beta * lambdascalar1;

    // read history
    const int gp = params.get<int>("gp", -1);
    if (gp == -1) dserror("no Gauss point number provided in material");
    LINALG::Matrix<NUM_STRESS_3D, 1> S_n(histstresslast_->at(gp));
    LINALG::Matrix<NUM_STRESS_3D, 1> Q_n(histartstresslast_->at(gp));

    // calculate artificial viscous stresses Q
    Q.Update(lambdascalar2, Q_n, 1.0);
    Q.Update(beta, *stress, 1.0);
    Q.Update(-beta, S_n, 1.0);
    Q.Scale(lambdascalar1);  // Q^(n+1) = lambdascalar1* [lambdascalar2* Q^n + beta*(S^(n+1) - S^n)]


    // update history
    histstresscurr_->at(gp) = *stress;
    histartstresscurr_->at(gp) = Q;

    // viscous constitutive tensor
    // contribution : Cmat_vis = Cmat_inf*deltascalar
    cmatq.Update(deltascalar, *cmat, 1.0);
  }
  else if (solve == "CONVOL")
  {
    // initialize scalars
    double xiscalar1(true);
    double xiscalar2(true);

    // get time algorithmic parameters
    // NOTE: dt can be zero (in restart of STI) for Generalized Maxwell model
    // there is no special treatment required. Adaptation for Kelvin-Voigt were necessary.
    double dt = params.get<double>("delta time");  // TIMESTEP in the .dat file
    // evaluate scalars to compute
    // Q_(n+1)=\exp(2*(-dt/(2*tau)))*Q_n+exp(-dt/(2*tau))*beta*(S_(n+1)-S_n)

    xiscalar1 = exp(-dt / tau);
    xiscalar2 = exp(-dt / (2 * tau)) * beta;

    // read history
    const int gp = params.get<int>("gp", -1);
    if (gp == -1) dserror("no Gauss point number provided in material");
    LINALG::Matrix<NUM_STRESS_3D, 1> S_n(histstresslast_->at(gp));
    LINALG::Matrix<NUM_STRESS_3D, 1> Q_n(histartstresslast_->at(gp));

    // calculate artificial viscous stresses Q
    Q.Update(xiscalar1, Q_n, 1.0);
    Q.Update(xiscalar2, *stress, 1.0);
    Q.Update(-xiscalar2, S_n, 1.0);  // Q_(n+1) = xiscalar1* Q_n + xiscalar2*(S_(n+1) - S_n)

    // update history
    histstresscurr_->at(gp) = *stress;
    histartstresscurr_->at(gp) = Q;

    // viscous constitutive tensor
    // contribution : Cmat_vis = Cmat_inf*beta*exp(-dt/(2*tau))
    cmatq.Update(xiscalar2, *cmat, 1.0);
  }
  else
    dserror("Invalid input. Try valid input OST or CONVOL");
  return;
}  // end EvaluateViscoGenMax

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ViscoElastHyper::EvaluateViscoGeneralizedGenMax(LINALG::Matrix<6, 1>& Q,
    LINALG::Matrix<6, 6>& cmatq, Teuchos::ParameterList& params,
    const LINALG::Matrix<6, 1>* glstrain, const int eleGID)
{
  int numbranch = -1;
  double tau = -1.0;
  std::string solve = "";
  const std::vector<int>* matids = NULL;
  std::vector<Teuchos::RCP<MAT::ELASTIC::Summand>> branchpotsum(
      0);  // vector of summands in one branch
  std::vector<std::vector<Teuchos::RCP<MAT::ELASTIC::Summand>>> branchespotsum(
      0);  // vector for each branch of vectors of summands in each branch

  // get parameters of ViscoGeneralizedGenMax
  for (unsigned int p = 0; p < potsum_.size(); ++p)
  {
    Teuchos::RCP<MAT::ELASTIC::GeneralizedGenMax> GeneralizedGenMax =
        Teuchos::rcp_dynamic_cast<MAT::ELASTIC::GeneralizedGenMax>(potsum_[p]);

    if (GeneralizedGenMax != Teuchos::null)
    {
      GeneralizedGenMax->ReadMaterialParameters(numbranch, matids, solve);
      branchespotsum = GeneralizedGenMax->GetBranchespotsum();
    }
  }

  LINALG::Matrix<6, 6> cmatqbranch(true);
  std::vector<LINALG::Matrix<6, 1>> S(numbranch);
  std::vector<LINALG::Matrix<6, 1>> Qbranch(numbranch);
  std::vector<LINALG::Matrix<6, 1>> S_n(numbranch);
  std::vector<LINALG::Matrix<6, 1>> Q_n(numbranch);

  // read history
  const int gp = params.get<int>("gp", -1);
  if (gp == -1) dserror("no Gauss point number provided in material");
  S_n = histbranchelaststresslast_->at(gp);
  Q_n = histbranchstresslast_->at(gp);


  // check size, correct if needed (only needed in the first time step)
  if (S_n.size() != (unsigned int)(abs(numbranch)))
  {
    std::vector<LINALG::Matrix<6, 1>> S_n(numbranch);
    for (int branch = 0; branch < numbranch; branch++) S_n.at(branch).Clear();
    histbranchelaststresslast_->at(gp) = S_n;
  }

  if (Q_n.size() != (unsigned int)(abs(numbranch)))
  {
    std::vector<LINALG::Matrix<6, 1>> Q_n(numbranch);
    for (int branch = 0; branch < numbranch; branch++) Q_n.at(branch).Clear();
    histbranchstresslast_->at(gp) = Q_n;
  }

  S_n = histbranchelaststresslast_->at(gp);
  Q_n = histbranchstresslast_->at(gp);

  // save switches
  SummandProperties branchProperties;

  /////////////////////////////////////////////////
  // Loop over all viscoelastic Maxwell branches //
  /////////////////////////////////////////////////
  for (int i = 0; i < numbranch; ++i)
  {
    // get parameter of each visco branch
    branchpotsum = branchespotsum[i];

    branchProperties.Clear();
    ElastHyperProperties(branchpotsum, branchProperties);

    if (isovisco_)
      dserror("case isovisco for branch in generalized Maxwell model not yet considered!");
    if (branchProperties.anisoprinc)
      dserror("case anisoprinc for branch in generalized Maxwell model not yet considered!");
    if (branchProperties.anisomod)
      dserror("case anisomod for branch in generalized Maxwell model not yet considered!");

    LINALG::Matrix<6, 1> C_strain(true);
    LINALG::Matrix<6, 1> iC_strain(true);
    LINALG::Matrix<6, 1> modrcg(true);
    LINALG::Matrix<3, 1> prinv(true);
    LINALG::Matrix<3, 1> dPI(true);
    LINALG::Matrix<6, 1> ddPII(true);

    EvaluateRightCauchyGreenStrainLikeVoigt(*glstrain, C_strain);
    VStrainUtils::InverseTensor(C_strain, iC_strain);

    VStrainUtils::InvariantsPrincipal(prinv, C_strain);
    ElastHyperEvaluateInvariantDerivatives(
        prinv, dPI, ddPII, branchpotsum, branchProperties, eleGID);

    // blank resulting quantities
    // ... even if it is an implicit law that cmat is zero upon input
    S.at(i).Clear();
    cmatqbranch.Clear();

    // build stress response and elasticity tensor
    LINALG::Matrix<NUM_STRESS_3D, 1> stressiso(true);
    LINALG::Matrix<NUM_STRESS_3D, NUM_STRESS_3D> cmatiso(true);
    ElastHyperAddIsotropicStressCmat(stressiso, cmatiso, C_strain, iC_strain, prinv, dPI, ddPII);
    S.at(i).Update(1.0, stressiso, 1.0);
    cmatqbranch.Update(1.0, cmatiso, 1.0);

    // get Parameter of ViscoGeneralizedGenMax
    for (unsigned int q = 0; q < branchpotsum.size(); ++q)
    {
      Teuchos::RCP<MAT::ELASTIC::ViscoPart> ViscoPart =
          Teuchos::rcp_dynamic_cast<MAT::ELASTIC::ViscoPart>(branchpotsum[q]);
      if (ViscoPart != Teuchos::null) ViscoPart->ReadMaterialParameters(tau);
    }

    // make sure Qbranch in this branch is empty
    Qbranch.at(i).Clear();
    double deltascalar = 1.0;
    if (solve == "OST")
    {
      // initialize scalars
      double lambdascalar1(true);
      double lambdascalar2(true);
      double theta = 0.5;

      // get theta of global time integration scheme to use it here
      // if global time integration scheme is not ONESTEPTHETA, theta is by default = 0.5 (abirzle
      // 09/14)
      std::string dyntype =
          DRT::Problem::Instance()->StructuralDynamicParams().get<std::string>("DYNAMICTYP");
      if (dyntype == "OneStepTheta")
        theta = DRT::Problem::Instance()
                    ->StructuralDynamicParams()
                    .sublist("ONESTEPTHETA")
                    .get<double>("THETA");

      // get time algorithmic parameters
      // NOTE: dt can be zero (in restart of STI) for Generalized Maxwell model
      // there is no special treatment required. Adaptation for Kelvin-Voigt were necessary.
      double dt = params.get<double>("delta time");  // TIMESTEP in the .dat file

      // evaluate scalars to compute
      // Q^(n+1) = tau/(tau+theta*dt) [(tau-dt+theta*dt)/tau Q + beta(S^(n+1) - S^n)]
      lambdascalar1 = tau / (tau + theta * dt);
      lambdascalar2 = (tau - dt + theta * dt) / tau;

      // see EvaluateViscoGenMax
      deltascalar = lambdascalar1;

      // calculate artificial viscous stresses Q
      // Q_(n+1) = lambdascalar1*[lamdascalar2* Q_n + (Sa_(n+1) - Sa_n)]
      Qbranch.at(i).Update(lambdascalar2, Q_n.at(i), 1.0);
      Qbranch.at(i).Update(1.0, S.at(i), 1.0);
      Qbranch.at(i).Update(-1.0, S_n.at(i), 1.0);
      Qbranch.at(i).Scale(lambdascalar1);
    }
    else if (solve == "CONVOL")
    {
      // initialize scalars
      double xiscalar1(true);
      double xiscalar2(true);
      double dt = params.get<double>("delta time");

      xiscalar1 = exp(-dt / tau);
      xiscalar2 = exp(-dt / (2 * tau));

      deltascalar = xiscalar2;

      // calculate artificial stresses Q
      // Q_(n+1) = xiscalar1* Q_n + xiscalar2*(Sa_(n+1) - Sa_n)
      Qbranch.at(i).Update(xiscalar1, Q_n.at(i), 1.0);
      Qbranch.at(i).Update(xiscalar2, S.at(i), 1.0);
      Qbranch.at(i).Update(-xiscalar2, S_n.at(i), 1.0);
    }
    else
      dserror("Invalid input. Try valid input THETA or CONVOL");

    // sum up branches
    Q.Update(1.0, Qbranch.at(i), 1.0);
    cmatq.Update(deltascalar, cmatqbranch, 1.0);

  }  // end for loop over branches

  // update history
  histbranchelaststresscurr_->at(gp) = S;
  histbranchstresscurr_->at(gp) = Qbranch;
}  // EvaluateViscoGeneralizedGenMax


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ViscoElastHyper::EvaluateViscoFract(LINALG::Matrix<6, 1> stress,
    LINALG::Matrix<6, 6> cmat, LINALG::Matrix<6, 1>& Q, LINALG::Matrix<6, 6>& cmatq,
    Teuchos::ParameterList& params)
{
  // initialize parameters
  double tau(true);
  double alpha(true);
  double beta(true);
  // string not used in viscofract
  std::string solve = "";

  // read material parameters of viscofract-material
  for (unsigned int p = 0; p < potsum_.size(); ++p)
  {
    potsum_[p]->ReadMaterialParametersVisco(tau, beta, alpha, solve);
  }

  // safety checks for alpha
  if (alpha == 1) dserror("Alpha cannot be 1 in fractional viscoelasticity. Use Genmax instead.");

  if (alpha < 0) dserror("Alpha has to be between 0 and 1 in fractional viscoelasticity.");

  // get time algorithmic parameters
  // NOTE: dt can be zero (in restart of STI) for this model
  // there is no special treatment required.
  double dt = params.get<double>("delta time");


  // read history of last time step at gp
  // -> Q_n and history size
  const int gp = params.get<int>("gp", -1);
  if (gp == -1) dserror("No Gauss point number provided in material");

  int hs = histfractartstresslastall_->at(0).size();  // history size
  LINALG::Matrix<NUM_STRESS_3D, 1> Q_n(histfractartstresslastall_->at(gp).at(hs - 1));


  // calculate artificial history stress Qq with weights b_j
  // Qq = sum[j=1 up to j=n][b_j*Q_(n+1-j)] (short: b*Qj)

  // b_j = gamma(j-alpha)/[gamma(-alpha)*gamma(j+1)]
  // with recurstion formula for gamma functions b_j shortens to
  // b_j = (j-1-alpha)/j * b_(j-1)
  // with b_0 = 1 and b_1 = -alpha ...
  double bj = 1.;   // b_0=1
  double fac = 1.;  // pre-factor (j-1-alpha)/j  for calculation of b
  LINALG::Matrix<NUM_STRESS_3D, 1> Qq(true);

  // j=1...n, hs=n
  for (int j = 1; j <= hs; j++)
  {
    fac = (j - 1. - alpha) / j;
    bj = bj * fac;

    LINALG::Matrix<NUM_STRESS_3D, 1> Qj(
        histfractartstresslastall_->at(gp).at(hs - j));  //(hs-1) correlates with last entry = n
    Qq.Update(bj, Qj, 1.0);
  }


  // calculate artificial stress Q

  // Version 1: As in Adolfson and Enelund (2003): Fractional Derivative Visocelasticity at Large
  // Deformations
  //  // initialize and evaluate scalars to compute
  //  // Q^(n+1) = [((dt/tau)^alpha)/(1+theta*(dt/tau)^alpha)]*[theta*S^(n+1)+(1-theta)(S^n-Q^n)]-
  //  //           [1/(1+theta*(dt/tau)^alpha)]*Qq^n


  // Version 2: Anna's Version of calculation
  // Difference:  1.) No one-step theta schema necessary
  //              2.) Introduce beta
  // Q^(n+1) = (dt^alpha / (dt^alpha + tau^alpha))*S^(n+1) - (tau^alpha / (dt^alpha +
  // tau^alpha))*Qq^n
  double dtalpha = std::pow(dt, alpha);
  double taualpha = std::pow(tau, alpha);
  double lambdascalar1 = dtalpha / (dtalpha + taualpha);
  double lambdascalar2 = -1. * taualpha / (dtalpha + taualpha);

  Q.Update(lambdascalar1 * beta, stress, 0.);
  Q.Update(lambdascalar2, Qq, 1.);


  // update history for next step
  histfractartstresscurr_->at(gp) = Q;  // Q_n+1


  // calculate final stress here and in Evaluate
  // S = elastic stress of Psi
  // S_2 = S ; S_1 = beta*S ; Q = Q(S1) = Q(beta*S)
  // S_final = S + beta*S - Q(beta*S)
  Q.Update(beta, stress, -1.);

  // viscos constitutive tensor
  cmatq.Update(lambdascalar1 * beta, cmat, 0.);  // contribution of Q
  cmatq.Update(beta, cmat, -1.);


  return;
}
