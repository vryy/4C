/*----------------------------------------------------------------------*/
/*!
\file elast_coupanisoneohooke.cpp
\brief


the input line should read
  MAT 1 CoupAnisoNeoHooke C 100 GAMMA 35.0

<pre>
Maintainer: Susanna Tinkl
            tinkl@lnm.mw.tum.de
            089/289 15265
</pre>
*/

/*----------------------------------------------------------------------*/
/* headers */
#include "elast_coupanisoneohooke.H"
#include "../drt_lib/drt_linedefinition.H"

/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
MAT::ELASTIC::PAR::CoupAnisoNeoHooke::CoupAnisoNeoHooke(
  Teuchos::RCP<MAT::PAR::Material> matdata
  )
: Parameter(matdata),
  c_(matdata->GetDouble("C")),
  gamma_(matdata->GetDouble("GAMMA"))
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
                                                  vector<char>::size_type& position)
{
  ExtractfromPack(position,data,a_);
  ExtractfromPack(position,data,A_);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::CoupAnisoNeoHooke::Setup(DRT::INPUT::LineDefinition* linedef)
{
  // fibers aligned in local element cosy with gamma around circumferential direction
  // -> check whether element supports local element cosy
  vector<double> rad;
  vector<double> axi;
  vector<double> cir;
  if (linedef->HaveNamed("RAD") and
      linedef->HaveNamed("AXI") and
      linedef->HaveNamed("CIR"))
  {
    // read local (cylindrical) cosy-directions at current element
    // basis is local cosy with third vec e3 = circumferential dir and e2 = axial dir
    LINALG::Matrix<3,3> locsys(true);
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

    SetFiberVecs(locsys);
  }
  else
  {
    dserror("Reading of element local cosy for anisotropic materials failed");
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::CoupAnisoNeoHooke::AddStressAnisoPrincipal(
    const LINALG::Matrix<6,1> rcg,
    LINALG::Matrix<6,6>& cmat,
    LINALG::Matrix<6,1>& stress
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
void MAT::ELASTIC::CoupAnisoNeoHooke::SetFiberVecs(
    LINALG::Matrix<3,3> locsys
)
{
  if ((params_->gamma_<-90) || (params_->gamma_ >90)) dserror("Fiber angle not in [-90,90]");
  //convert
  const double gamma = (params_->gamma_*PI)/180.;

  for (int i = 0; i < 3; ++i)
  {
    // a = cos gamma e3 + sin gamma e2
    a_(i) = cos(gamma)*locsys(i,2) + sin(gamma)*locsys(i,1);
  }

  for (int i = 0; i < 3; ++i)
    A_(i) = a_(i)*a_(i);

  A_(3) = a_(0)*a_(1); A_(4) = a_(1)*a_(2); A_(5) = a_(0)*a_(2);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::CoupAnisoNeoHooke::GetFiberVecs(
    std::vector<LINALG::Matrix<3,1> >& fibervecs ///< vector of all fiber vectors
)
{
  fibervecs.push_back(a_);
}
