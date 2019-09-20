/*----------------------------------------------------------------------*/
/*! \file
\brief Active stress material

\maintainer Amadeus Gebauer

\level 2
*/
/*---------------------------------------------------------------------*/

/*----------------------------------------------------------------------*/
/* headers */
#include "elast_anisoactivestress_evolution.H"

#include "../drt_mat/matpar_material.H"
#include "../drt_mat/material.H"
#include "../drt_lib/standardtypes_cpp.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_mat/material_service.H"

/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
MAT::ELASTIC::PAR::AnisoActiveStress_Evolution::AnisoActiveStress_Evolution(
    const Teuchos::RCP<MAT::PAR::Material>& matdata)
    : ParameterAniso(matdata),
      sigma_(matdata->GetDouble("SIGMA")),
      tauc0_(matdata->GetDouble("TAUC0")),
      maxactiv_(matdata->GetDouble("MAX_ACTIVATION")),
      minactiv_(matdata->GetDouble("MIN_ACTIVATION")),
      activationthreshold_(matdata->GetDouble("ACTIVATION_THRES")),
      sourceactiv_(matdata->GetInt("SOURCE_ACTIVATION")),
      strain_dep_(matdata->GetInt("STRAIN_DEPENDENCY")),
      lambda_lower_(matdata->GetDouble("LAMBDA_LOWER")),
      lambda_upper_(matdata->GetDouble("LAMBDA_UPPER")),
      gamma_(matdata->GetDouble("GAMMA")),
      theta_(matdata->GetDouble("THETA")),
      init_(matdata->GetInt("INIT")),
      adapt_angle_(matdata->GetInt("ADAPT_ANGLE"))
{
}

/*----------------------------------------------------------------------*
 *         Constructor Material Parameter Class                         *
 *----------------------------------------------------------------------*/
MAT::ELASTIC::AnisoActiveStress_Evolution::AnisoActiveStress_Evolution(
    MAT::ELASTIC::PAR::AnisoActiveStress_Evolution* params)
    : Anisotropy(1), params_(params), tauc_np_(0.0), tauc_n_(0.0)
{
  Anisotropy::Initialize(params->init_, params->StructuralTensorStrategy());
}

/*----------------------------------------------------------------------*
 *            Constructor Material Class                               *
 *----------------------------------------------------------------------*/
void MAT::ELASTIC::AnisoActiveStress_Evolution::PackSummand(DRT::PackBuffer& data) const
{
  AddtoPack(data, tauc_n_);

  PackAnisotropy(data);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::AnisoActiveStress_Evolution::UnpackSummand(
    const std::vector<char>& data, std::vector<char>::size_type& position)
{
  ExtractfromPack(position, data, tauc_n_);

  UnpackAnisotropy(data, position);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::AnisoActiveStress_Evolution::Setup(
    int numgp, DRT::INPUT::LineDefinition* linedef)
{
  Anisotropy::SetNumberOfGaussPoints(numgp);
  Anisotropy::ReadAnisotropyFromElement(linedef);

  // Setup of active stress model
  tauc_n_ = params_->tauc0_;
  tauc_np_ = params_->tauc0_;

  // reasonability check...
  if (params_->strain_dep_ != 0)
  {
    if (params_->lambda_lower_ >= params_->lambda_upper_)
    {
      dserror(
          "LAMBDA_LOWER should be lower than LAMBDA_UPPER. Seems reasonable, doesn't it? Dude...");
    }
  }
}

void MAT::ELASTIC::AnisoActiveStress_Evolution::PostSetup(Teuchos::ParameterList& params)
{
  Summand::PostSetup(params);

  Anisotropy::ReadAnisotropyFromParameterList(params);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::AnisoActiveStress_Evolution::AddStressAnisoPrincipal(
    const LINALG::Matrix<6, 1>& rcg, LINALG::Matrix<6, 6>& cmat, LINALG::Matrix<6, 1>& stress,
    Teuchos::ParameterList& params, const int eleGID)
{
  // Virtual GP (is zero for element fibers, otherwise it is the current GP)
  int virtualgp = GetGPId(params);

  LINALG::Matrix<6, 1> A = GetStructuralTensor_stress(virtualgp, 0);

  double dt = params.get("delta time", -1.0);
  if (dt < 0.0)
  {
    dserror("Parameter 'delta time' could not be read!");
  }

  double activationFunction = 0.0;
  Teuchos::RCP<MAT::Material> scatramat;
  if (params_->sourceactiv_ == 0)
  {
    activationFunction = params.get<double>("scalar", 0.0);
  }
  else if (params_->sourceactiv_ == -1)
  {
    if (params.isParameter("gp_conc"))
    {
      // get pointer to vector containing the scalar values at the Gauss points
      Teuchos::RCP<std::vector<std::vector<double>>> conc =
          params.get<Teuchos::RCP<std::vector<std::vector<double>>>>("gp_conc");
      // safety check
      if (conc == Teuchos::null)
      {
        dserror("No concentration from scatra provided for the activation");
      }
      // get Gauss point number
      const int gp = params.get<int>("gp", -1);
      // safety check
      if (gp == -1)
      {
        dserror("No Gauss point number provided in material");
      }
      // check if activation is above the given threshold
      if (conc->at(gp)[0] >= params_->activationthreshold_)
      {
        activationFunction = 1.0;
      }
    }
  }
  else
  {
    double totaltime = params.get<double>("total time", -1.0);
    if (totaltime < 0.0)
    {
      dserror("Parameter 'total time' could not be read!");
    }
    Teuchos::RCP<std::vector<double>> pos_ =
        params.get<Teuchos::RCP<std::vector<double>>>("position");
    const double* coordgpref_ = &(*pos_)[0];
    activationFunction = DRT::Problem::Instance()
                             ->Funct(params_->sourceactiv_ - 1)
                             .Evaluate(0, coordgpref_, totaltime);
  }

  double lambda = 0.0;
  double scale = 0.0;
  double n0 = 1.0;
  if (params_->strain_dep_ != 0)
  {
    // squared stretch along fiber direction
    double I4;
    I4 = A(0) * rcg(0) + A(1) * rcg(1) + A(2) * rcg(2) + A(3) * rcg(3) + A(4) * rcg(4) +
         A(5) * rcg(5);
    // stretch along fiber direction
    lambda = sqrt(I4);
    scale = -4. / (pow((params_->lambda_lower_ - params_->lambda_upper_), 2.0));
    // Frank-Starling factor n0 for strain-dependent contractility
    if (lambda >= params_->lambda_lower_ and lambda <= params_->lambda_upper_)
    {
      n0 = scale * (lambda - params_->lambda_lower_) * (lambda - params_->lambda_upper_);
    }
    else
    {
      n0 = 0.0;
    }
  }

  activationFunction =
      activationFunction * (params_->maxactiv_ - params_->minactiv_) + params_->minactiv_;
  double abs_u_ = abs(activationFunction);
  double absplus_u_ = abs_u_ * static_cast<double>(activationFunction > 0.0);
  tauc_np_ = (tauc_n_ / dt + n0 * params_->sigma_ * absplus_u_) / (1 / dt + abs_u_);
  stress.Update(tauc_np_, A, 1.0);

  // only contribution to cmat if we have strain dependency!
  if (params_->strain_dep_ != 0)
  {
    double dtauc_np_dC;
    // Cmat_active = 2 * dS_active / dC = 2 * dtau(t) / dC \otimes (a_0 \otimes a_0)
    if (lambda >= params_->lambda_lower_ and lambda <= params_->lambda_upper_)
    {
      dtauc_np_dC = 2. * scale *
                    (1. - (params_->lambda_upper_ + params_->lambda_lower_) / (2. * lambda)) *
                    params_->sigma_ * absplus_u_ / (1 / dt + abs_u_);
    }
    else
    {
      dtauc_np_dC = 0.0;
    }
    cmat.MultiplyNT(dtauc_np_dC, A, A, 1.0);
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::AnisoActiveStress_Evolution::GetFiberVecs(
    std::vector<LINALG::Matrix<3, 1>>& fibervecs  ///< vector of all fiber vectors
)
{
  if (params_->init_ == INIT_MODE_NODAL_FIBERS)
  {
    // This method expects constant fibers within this element but the init mode is such that
    // fibers are defined on the Gauss points
    // We therefore cannot return sth here.

    // ToDo: This may needs improvements later on if needed!
    return;
  }

  fibervecs.push_back(GetFiber(GPDEFAULT, 0));
}

/*----------------------------------------------------------------------*
 |  Update internal stress variables              (public)         05/08|
 *----------------------------------------------------------------------*/
void MAT::ELASTIC::AnisoActiveStress_Evolution::Update() { tauc_n_ = tauc_np_; }

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::AnisoActiveStress_Evolution::SetFiberVecs(
    const double newgamma, const LINALG::Matrix<3, 3>& locsys, const LINALG::Matrix<3, 3>& defgrd)
{
  if ((params_->gamma_ < -90) || (params_->gamma_ > 90))
  {
    dserror("Fiber angle not in [-90,90]");
  }
  // convert
  double gamma = (params_->gamma_ * PI) / 180.;

  if (params_->adapt_angle_ != 0 && newgamma != -1.0)
  {
    if (gamma * newgamma < 0.0)
    {
      gamma = -1.0 * newgamma;
    }
    else
    {
      gamma = newgamma;
    }
  }

  LINALG::Matrix<3, 1> ca(true);
  for (int i = 0; i < 3; ++i)
  {
    // a = cos gamma e3 + sin gamma e2
    ca(i) = cos(gamma) * locsys(i, 2) + sin(gamma) * locsys(i, 1);
  }
  // pull back in reference configuration
  LINALG::Matrix<3, 1> a_0(true);
  LINALG::Matrix<3, 3> idefgrd(true);
  idefgrd.Invert(defgrd);

  a_0.Multiply(idefgrd, ca);
  a_0.Scale(1.0 / a_0.Norm2());

  SetFibers(GPDEFAULT, std::vector<LINALG::Matrix<3, 1>>(1, a_0));
}

void MAT::ELASTIC::AnisoActiveStress_Evolution::SetupFiberByCosy(LINALG::Matrix<3, 3>& locsys)
{
  LINALG::Matrix<3, 3> Id(false);
  MAT::IdentityMatrix(Id);

  // final setup of fiber data
  SetFiberVecs(0.0, locsys, Id);
}

void MAT::ELASTIC::AnisoActiveStress_Evolution::DoFiberInitialization()
{
  // fibers aligned in YZ-plane with gamma around Z in global cartesian cosy
  LINALG::Matrix<3, 3> Id(true);
  LINALG::Matrix<3, 3> locsys(true);
  MAT::IdentityMatrix(Id);

  // To realize a full rotated fiber orientation and to keep the general structure of
  // SetFiberVec() the input of locsys has to be adapted if one sets
  //               1           0                 sin(theta_)
  //  locsys := [  0       sin(theta_)               0                  ]
  //               0  sin(gamma)*cos(theta_)    cos(gamma_)*cos(theta_)
  // The call of SetFiberVec() will leed to the following fiber direction
  // a = cos(gamma_)*locsys(:,2) + sin(gamma_)*locsys(:,1)
  //         cos(gamma_)*sin(theta_)               0 cos(gamma_)*sin(theta_)
  //   = [             0              ] + [  sin(gamma_)*sin(theta_)  ] = [
  //   sin(gamma_)*sin(theta_) ] =: sperical coordinates
  //         cos(gamma)^2*cos(theta_)        sin(gamma_)^2*cos(theta_)               cos(theta_)
  //
  {
    // Local initialization of spherical angles
    double theta = (params_->theta_);
    double gamma = (params_->gamma_);
    if (gamma < 0.0 || gamma > 180.0 || abs(theta) > 180.0)
    {
      dserror(
          "Wrong choice of sherical coodinates. Correct domain is gamma in [0,180], theta in "
          "[-180, 180]");
    }
    // conversion to radian measure
    theta = (theta * PI) / 180.0;
    gamma = (gamma * PI) / 180.0;
    locsys(1, 1) = sin(theta);
    locsys(2, 1) = sin(gamma) * cos(theta);
    locsys(0, 2) = sin(theta);
    locsys(2, 2) = cos(gamma) * cos(theta);
  }
  SetFiberVecs(-1.0, locsys, Id);
}
