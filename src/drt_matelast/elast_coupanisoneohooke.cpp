/*----------------------------------------------------------------------*/
/*!
\file elast_coupanisoneohooke.cpp
\brief


the input line should read
  MAT 1 CoupAnisoNeoHooke C 100 GAMMA 35.0 INIT 0 ADAPT_ANGLE 0

<pre>
Maintainer: Susanna Tinkl
            tinkl@lnm.mw.tum.de
            089/289 15265
</pre>
*/

/*----------------------------------------------------------------------*/
/* headers */
#include "elast_coupanisoneohooke.H"
#include "../drt_mat/matpar_material.H"
#include "../drt_lib/standardtypes_cpp.H"
#include "../drt_lib/drt_linedefinition.H"

/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
MAT::ELASTIC::PAR::CoupAnisoNeoHooke::CoupAnisoNeoHooke(
  Teuchos::RCP<MAT::PAR::Material> matdata
  )
: Parameter(matdata),
  c_(matdata->GetDouble("C")),
  gamma_(matdata->GetDouble("GAMMA")),
  init_(matdata->GetInt("INIT")),
  adapt_angle_(matdata->GetInt("ADAPT_ANGLE"))
{
}


Teuchos::RCP<MAT::Material> MAT::ELASTIC::PAR::CoupAnisoNeoHooke::CreateMaterial()
{
  return Teuchos::null;
  //return Teuchos::rcp( new MAT::ELASTIC::CoupAnisoNeoHooke( this ) );
}


/*----------------------------------------------------------------------*
 |  Constructor                                   (public)  bborn 04/09 |
 *----------------------------------------------------------------------*/
MAT::ELASTIC::CoupAnisoNeoHooke::CoupAnisoNeoHooke()
  : Summand(),
    params_(NULL)
{
}


/*----------------------------------------------------------------------*
 |  Constructor                             (public)   bborn 04/09 |
 *----------------------------------------------------------------------*/
MAT::ELASTIC::CoupAnisoNeoHooke::CoupAnisoNeoHooke(MAT::ELASTIC::PAR::CoupAnisoNeoHooke* params)
  : params_(params)
{
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::CoupAnisoNeoHooke::PackSummand(DRT::PackBuffer& data) const
{
  AddtoPack(data,a_);
  AddtoPack(data,A_);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::CoupAnisoNeoHooke::UnpackSummand(const std::vector<char>& data,
																										std::vector<char>::size_type& position)
{
  ExtractfromPack(position,data,a_);
  ExtractfromPack(position,data,A_);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::CoupAnisoNeoHooke::Setup(DRT::INPUT::LineDefinition* linedef)
{
  // path if fibers aren't given in path file
  if (params_->init_ == 0)
  {
    // fibers aligned in YZ-plane with gamma around Z in global cartesian cosy
    LINALG::Matrix<3,3> Id(true);
    for (int i=0; i<3; i++)
      Id(i,i) = 1.0;
    SetFiberVecs(-1.0,Id,Id);
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
      ReadFiber1(linedef);
    }

    // error path
    else
    {
      dserror("Reading of element local cosy for anisotropic materials failed");
    }

    // Setup of structural tensors
    for (int i = 0; i < 3; ++i)
      A_(i) = a_(i)*a_(i);
    A_(3) = a_(0)*a_(1); A_(4) = a_(1)*a_(2); A_(5) = a_(0)*a_(2);

  }
  else
    dserror("INIT mode not implemented");
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::CoupAnisoNeoHooke::AddStressAnisoPrincipal(
    const LINALG::Matrix<6,1> rcg,
    LINALG::Matrix<6,6>& cmat,
    LINALG::Matrix<6,1>& stress,
    Teuchos::ParameterList& params
)
{
  double c=params_->c_;

  double gamma = 2.*c;
  stress.Update(gamma, A_, 1.0);

  // no contribution to cmat
  // double delta = 0.0;
  // cmat.MultiplyNT(delta, A_, A_, 1.0);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::CoupAnisoNeoHooke::GetFiberVecs(
    std::vector<LINALG::Matrix<3,1> >& fibervecs ///< vector of all fiber vectors
)
{
  fibervecs.push_back(a_);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
// Function which reads in the given fiber value due to the FIBER1 nomenclature
void MAT::ELASTIC::CoupAnisoNeoHooke::ReadFiber1(
    DRT::INPUT::LineDefinition* linedef
)
{
  std::vector<double> fiber1;
  linedef->ExtractDoubleVector("FIBER1",fiber1);
  double f1norm=0.;
  //normalization
  for (int i = 0; i < 3; ++i)
  {
    f1norm += fiber1[i]*fiber1[i];
  }
  f1norm = sqrt(f1norm);

  // fill final fiber vector
  for (int i = 0; i < 3; ++i)
    a_(i) = fiber1[i]/f1norm;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
// Function which reads in the given fiber value due to the CIR-AXI-RAD nomenclature
void MAT::ELASTIC::CoupAnisoNeoHooke::ReadRadAxiCir(
    DRT::INPUT::LineDefinition* linedef,
    LINALG::Matrix<3,3>& locsys
)
{
  // fibers aligned in local element cosy with gamma around circumferential direction
  // -> check whether element supports local element cosy
  if (linedef->HaveNamed("RAD") and
      linedef->HaveNamed("AXI") and
      linedef->HaveNamed("CIR"))
  {
    // read local (cylindrical) cosy-directions at current element
    // basis is local cosy with third vec e3 = circumferential dir and e2 = axial dir
    std::vector<double> rad;
    std::vector<double> axi;
    std::vector<double> cir;

    linedef->ExtractDoubleVector("RAD",rad);
    linedef->ExtractDoubleVector("AXI",axi);
    linedef->ExtractDoubleVector("CIR",cir);
    double radnorm=0.; double axinorm=0.; double cirnorm=0.;

    for (int i = 0; i < 3; ++i)
    {
      radnorm += rad[i]*rad[i]; axinorm += axi[i]*axi[i]; cirnorm += cir[i]*cir[i];
    }
    radnorm = sqrt(radnorm); axinorm = sqrt(axinorm); cirnorm = sqrt(cirnorm);

    for (int i=0; i<3; ++i)
    {
      locsys(i,0) = rad[i]/radnorm;
      locsys(i,1) = axi[i]/axinorm;
      locsys(i,2) = cir[i]/cirnorm;
    }
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::CoupAnisoNeoHooke::SetFiberVecs(
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
}
