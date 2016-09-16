/*----------------------------------------------------------------------*/
/*!
\file elast_coupanisoneohooke_ActiveStress.cpp

\brief Active stress material

\maintainer Marc Hirschvogel

\level 2
*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*/
/* headers */
#include "elast_coupanisoneohooke_ActiveStress.H"
#include "../drt_mat/matpar_material.H"
#include "../drt_mat/material.H"
#include "../drt_lib/standardtypes_cpp.H"
#include "../drt_lib/drt_linedefinition.H"
#include "../drt_lib/drt_globalproblem.H"

/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
MAT::ELASTIC::PAR::CoupAnisoNeoHooke_ActiveStress::CoupAnisoNeoHooke_ActiveStress(
  Teuchos::RCP<MAT::PAR::Material> matdata
  )
: Parameter(matdata),
  sigma_(matdata->GetDouble("SIGMA")),
  tauc0_(matdata->GetDouble("TAUC0")),
  maxactiv_(matdata->GetDouble("MAX_ACTIVATION")),
  minactiv_(matdata->GetDouble("MIN_ACTIVATION")),
  activationthreshold_(matdata->GetDouble("ACTIVATION_THRES")),
  sourceactiv_(matdata->GetInt("SOURCE_ACTIVATION")),
  gamma_(matdata->GetDouble("GAMMA")),
  theta_(matdata->GetDouble("THETA")),
  init_(matdata->GetInt("INIT")),
  adapt_angle_(matdata->GetInt("ADAPT_ANGLE"))
{
}

/*----------------------------------------------------------------------*
 *         Constructor Material Parameter Class                         *
 *----------------------------------------------------------------------*/
MAT::ELASTIC::CoupAnisoNeoHooke_ActiveStress::CoupAnisoNeoHooke_ActiveStress(MAT::ELASTIC::PAR::CoupAnisoNeoHooke_ActiveStress* params)
  : params_(params)
{
}

/*----------------------------------------------------------------------*
 *            Constructor Material Class                               *
 *----------------------------------------------------------------------*/
void MAT::ELASTIC::CoupAnisoNeoHooke_ActiveStress::PackSummand(DRT::PackBuffer& data) const
{
  AddtoPack(data,a_);
  AddtoPack(data,A_);
  AddtoPack(data,tauc_last_);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::CoupAnisoNeoHooke_ActiveStress::UnpackSummand(
  const std::vector<char>& data,
  std::vector<char>::size_type& position
  )
{
  ExtractfromPack(position,data,a_);
  ExtractfromPack(position,data,A_);
  ExtractfromPack(position,data,tauc_last_);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::CoupAnisoNeoHooke_ActiveStress::Setup(DRT::INPUT::LineDefinition* linedef)
{
  // Setup of active stress model
  tauc_last_ = params_->tauc0_;
  tauc_ = params_->tauc0_;

  // path if fibers aren't given in .dat file
  if (params_->init_ == 0)
  {
    // fibers aligned in YZ-plane with gamma around Z in global cartesian cosy
    LINALG::Matrix<3,3> Id(true);
    LINALG::Matrix<3,3> locsys(true);
    for (int i=0; i<3; i++)
      Id(i,i) = 1.0;

    // To realize a full rotated fiber orientation and to keep the general structure of SetFiberVec() the
    // input of locsys has to be adapted
    // if one sets
    //               1           0                 sin(theta_)
    //  locsys := [  0       sin(theta_)               0                  ]
    //               0  sin(gamma)*cos(theta_)    cos(gamma_)*cos(theta_)
    // The call of SetFiberVec() will leed to the following fiber direction
    // a = cos(gamma_)*locsys(:,2) + sin(gamma_)*locsys(:,1)
    //         cos(gamma_)*sin(theta_)               0                          cos(gamma_)*sin(theta_)
    //   = [             0              ] + [  sin(gamma_)*sin(theta_)  ] = [   sin(gamma_)*sin(theta_) ] =: sperical coordinates
    //         cos(gamma)^2*cos(theta_)        sin(gamma_)^2*cos(theta_)               cos(theta_)
    //
    {
      // Local initialization of spherical angles
      double theta = (params_->theta_);
      double gamma = (params_->gamma_);
      if ( gamma < 0.0 || gamma > 180.0 || abs(theta) > 180.0 )
      {
        dserror("Wrong choice of sherical coodinates. Correct domain is gamma in [0,180], theta in [-180, 180]");
      }
      // conversion to radian measure
      theta = (theta*PI)/180.0;
      gamma = (gamma*PI)/180.0;
      locsys(1,1) = sin(theta);
      locsys(2,1) = sin(gamma)*cos(theta);
      locsys(0,2) = sin(theta);
      locsys(2,2) = cos(gamma)*cos(theta);
    }
    SetFiberVecs(-1.0,locsys,Id);
  }

  // path if fibers are given in .dat file
  else if (params_->init_ == 1)
  {

    // CIR-AXI-RAD nomenclature
    if (linedef->HaveNamed("RAD") and
        linedef->HaveNamed("AXI") and
        linedef->HaveNamed("CIR"))
    {
      // Read in of data
      LINALG::Matrix<3,3> locsys(true);
      ReadRadAxiCir(linedef, locsys);
      LINALG::Matrix<3,3> Id(true);
      for (int i=0; i<3; i++)
        Id(i,i) = 1.0;
      // final setup of fiber data
      SetFiberVecs(0.0,locsys,Id);
    }

    // FIBER1 nomenclature
    else if ( linedef->HaveNamed("FIBER1") )
    {
      // Read in of fiber data and setting fiber data
      ReadFiber(linedef, "FIBER1", a_);
      SetupStructuralTensor(a_,A_);
    }

    // error path
    else
    {
      dserror("Reading of element local cosy for anisotropic materials failed");
    }

  }
  else
    dserror("INIT mode not implemented");
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::CoupAnisoNeoHooke_ActiveStress::AddStressAnisoPrincipal(
    const LINALG::Matrix<6,1> rcg,
    LINALG::Matrix<6,6>& cmat,
    LINALG::Matrix<6,1>& stress,
    Teuchos::ParameterList& params,
    const int eleGID
)
{
  double dt = params.get("delta time",-1.0);
  if(dt<0.0) dserror("Parameter 'delta time' could not be read!");

  double activationFunction = 0.0;
  Teuchos::RCP<MAT::Material> scatramat;
  if(params_->sourceactiv_==0)
    activationFunction=params.get<double>("scalar",0.0);
  else if(params_->sourceactiv_==-1)
  {
    if (params.isParameter("gp_conc"))
    {
      // get pointer to vector containing the scalar values at the Gauss points
      Teuchos::RCP<std::vector<std::vector<double> > > conc =
          params.get< Teuchos::RCP<std::vector<std::vector<double> > > >("gp_conc");
      // safety check
      if (conc == Teuchos::null)
        dserror("No concentration from scatra provided for the activation");
      // get Gauss point number
      const int gp = params.get<int>("gp",-1);
      // safety check
      if (gp == -1)
        dserror("No Gauss point number provided in material");
      // check if activation is above the given threshold
      if(conc->at(gp)[0] >= params_->activationthreshold_)
        activationFunction = 1.0;
    }
  }
  else
  {
    double totaltime = params.get<double>("total time",-1.0);
    if(totaltime<0.0) dserror("Parameter 'total time' could not be read!");
    Teuchos::RCP<std::vector<double> >  pos_ = params.get<Teuchos::RCP<std::vector<double> > >("position");
    const double* coordgpref_ = &(*pos_)[0];
    activationFunction = DRT::Problem::Instance()->Funct(params_->sourceactiv_-1).Evaluate(0,coordgpref_,totaltime,NULL);
  }
  activationFunction = activationFunction*(params_->maxactiv_-params_->minactiv_)+params_->minactiv_;
  double abs_u_ = abs(activationFunction);
  double absplus_u_ = abs_u_*(activationFunction>0.0);
  tauc_ =  (tauc_last_/dt + params_->sigma_*absplus_u_)/(1/dt + abs_u_);
  stress.Update(tauc_, A_, 1.0);

   // no contribution to cmat
  // double delta = 0.0;
  // cmat.MultiplyNT(delta, A_, A_, 1.0);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::CoupAnisoNeoHooke_ActiveStress::GetFiberVecs(
    std::vector<LINALG::Matrix<3,1> >& fibervecs ///< vector of all fiber vectors
)
{
  fibervecs.push_back(a_);
}

/*----------------------------------------------------------------------*
 |  Update internal stress variables              (public)         05/08|
 *----------------------------------------------------------------------*/
void MAT::ELASTIC::CoupAnisoNeoHooke_ActiveStress::Update()
{
  tauc_last_ = tauc_;
  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::CoupAnisoNeoHooke_ActiveStress::SetFiberVecs(
    const double newgamma,
    const LINALG::Matrix<3,3> locsys,
    const LINALG::Matrix<3,3> defgrd
)
{
  if ((params_->gamma_<-90) || (params_->gamma_ >90)) dserror("Fiber angle not in [-90,90]");
  //convert
  double gamma = (params_->gamma_*PI)/180.;

  if (params_->adapt_angle_ && newgamma != -1.0)
  {
    if (gamma*newgamma < 0.0)
      gamma = -1.0 * newgamma;
    else
      gamma = newgamma;
  }

  LINALG::Matrix<3,1> ca(true);
  for (int i = 0; i < 3; ++i)
  {
    // a = cos gamma e3 + sin gamma e2
    ca(i) = cos(gamma)*locsys(i,2) + sin(gamma)*locsys(i,1);
  }
  // pull back in reference configuration
  LINALG::Matrix<3,1> a_0(true);
  LINALG::Matrix<3,3> idefgrd(true);
  idefgrd.Invert(defgrd);

  a_0.Multiply(idefgrd,ca);
  a_.Update(1./a_0.Norm2(),a_0);

  SetupStructuralTensor(a_,A_);
}
