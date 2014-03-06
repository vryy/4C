/*----------------------------------------------------------------------*/
/*!
 \file growth_scd.cpp

 \brief

 This file contains routines for an integration point based growth law.
It is derived from the growth law implemented in growth_ip.cpp with
additional scalar dependency (scd), e.g. nutrients.

See also Bachelor thesis Jasper Rieser (2013)

 <pre>
   Maintainer: Anh-Tu Vuong
               vuong@lnm.mw.tum.de
               http://www.lnm.mw.tum.de
               089 - 289-15264
 </pre>
 *----------------------------------------------------------------------*/


#include "growth_scd.H"
#include "growth_law.H"
#include "../drt_mat/matpar_bundle.H"
#include "../drt_lib/drt_globalproblem.H"   // for function Factory in Unpack
#include "../drt_lib/drt_utils_factory.H"   // for function Factory in Unpack


/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
MAT::PAR::GrowthScd::GrowthScd(
  Teuchos::RCP<MAT::PAR::Material> matdata
  )
: Growth(matdata),
  satcoeff_(matdata->GetDouble("SATCOEFF"))
{
}


Teuchos::RCP<MAT::Material> MAT::PAR::GrowthScd::CreateMaterial()
{
  return Teuchos::rcp(new MAT::GrowthScd(this));
}


MAT::GrowthScdType MAT::GrowthScdType::instance_;

DRT::ParObject* MAT::GrowthScdType::Create( const std::vector<char> & data )
{
  MAT::GrowthScd* grow = new MAT::GrowthScd();
  grow->Unpack(data);
  return grow;
}

/*----------------------------------------------------------------------*
 |  Constructor                                   (public)         02/10|
 *----------------------------------------------------------------------*/
MAT::GrowthScd::GrowthScd()
  : Growth(),
    params_(NULL)
{
}


/*----------------------------------------------------------------------*
 |  Copy-Constructor                             (public)          02/10|
 *----------------------------------------------------------------------*/
MAT::GrowthScd::GrowthScd(MAT::PAR::GrowthScd* params)
  : Growth(params),
    params_(params)
{
}


/*----------------------------------------------------------------------*
 |  Pack                                          (public)         02/10|
 *----------------------------------------------------------------------*/
void MAT::GrowthScd::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm( data );
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);
  // matid
  int matid = -1;
  if (params_ != NULL) matid = params_->Id();  // in case we are in post-process mode
  AddtoPack(data,matid);

  int numgp=0;
  if (isinit_)
  {
    numgp = dtheta_->size();   // size is number of gausspoints
  }
  AddtoPack(data,numgp);
  // Pack internal variables
  for (int gp = 0; gp < numgp; ++gp)
  {
    AddtoPack(data,detFe_->at(gp));
    AddtoPack(data,dtheta_->at(gp));
  }

  // Pack base class material
  Growth::Pack(data);

  return;
}


/*----------------------------------------------------------------------*
 |  Unpack                                        (public)         02/10|
 *----------------------------------------------------------------------*/
void MAT::GrowthScd::Unpack(const std::vector<char>& data)
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
        params_ = static_cast<MAT::PAR::GrowthScd*>(mat);
      else
        dserror("Type of parameter material %d does not fit to calling type %d", mat->Type(), MaterialType());
    }

  int numgp;
  ExtractfromPack(position,data,numgp);
  if (numgp == 0){ // no history data to unpack
    isinit_=false;
    if (position != data.size())
      dserror("Mismatch in size of data %d <-> %d",data.size(),position);
    return;
  }

  // unpack growth internal variables
  detFe_ = Teuchos::rcp(new std::vector<double> (numgp));
  dtheta_ = Teuchos::rcp(new std::vector<double> (numgp));
  for (int gp = 0; gp < numgp; ++gp) {
    double a;
    ExtractfromPack(position,data,a);
    detFe_->at(gp) = a;
    ExtractfromPack(position,data,a);
    dtheta_->at(gp) = a;
  }

  // extract base class material
  std::vector<char> basedata(0);
  Growth::ExtractfromPack(position,data,basedata);
  Growth::Unpack(basedata);


  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d",data.size(),position);

  return;
}

/*----------------------------------------------------------------------*
 |  Setup                                         (public)         02/10|
 *----------------------------------------------------------------------*/
void MAT::GrowthScd::Setup(int numgp, DRT::INPUT::LineDefinition* linedef)
{
  detFe_ = Teuchos::rcp(new std::vector<double> (numgp,1.0));
  dtheta_ = Teuchos::rcp(new std::vector<double> (numgp,0.0));

  //setup base class
  Growth::Setup(numgp, linedef);
  return;
}

/*----------------------------------------------------------------------*
 |  Update internal growth variables              (public)         02/10|
 *----------------------------------------------------------------------*/
void MAT::GrowthScd::Update()
{
//  const int histsize = theta_->size();
//  for (int i=0; i<histsize; i++)
//  {
//    thetaold_->at(i) = theta_->at(i);
//  }

  //update base class material
  Growth::Update();

  return;
}

/*----------------------------------------------------------------------*
 |  Evaluate Material                             (public)         02/10|
 *----------------------------------------------------------------------*
 The deformation gradient is decomposed into an elastic and growth part:
     F = Felastic * F_g
 Only the elastic part contributes to the stresses, thus we have to
 compute the elastic Cauchy Green Tensor Cdach and elastic 2PK stress Sdach.
 */
void MAT::GrowthScd::Evaluate
(
    const LINALG::Matrix<3,3>* defgrd,
    const LINALG::Matrix<6,1>* glstrain,
    Teuchos::ParameterList& params,
    LINALG::Matrix<6,1>* stress,
    LINALG::Matrix<6,6>* cmat,
    const int eleGID
)
{
  // get gauss point number
  const int gp = params.get<int>("gp",-1);
  if (gp == -1) dserror("no Gauss point number provided in material");

  concentration_ = params.get<double>("scalar", 0.0);

  //if(concentration_==-1.0) dserror("no scalar concentration provided for concentration dependent growth law");

  Growth::Evaluate(defgrd,glstrain,params,stress,cmat,eleGID);

  // build identity tensor I
  LINALG::Matrix<NUM_STRESS_3D,1> Id(true);
  for (int i = 0; i < 3; i++) Id(i) = 1.0;
  // right Cauchy-Green Tensor  C = 2 * E + I
  LINALG::Matrix<NUM_STRESS_3D,1> C(*glstrain);
  C.Scale(2.0);
  C += Id;

  const double theta = theta_->at(gp);
  // elastic right Cauchy-Green Tensor Cdach = F_g^-T C F_g^-1
  LINALG::Matrix<NUM_STRESS_3D,1> Cdach(C);
  Cdach.Scale(1.0/theta/theta);

  // determinate of F_e (necessary for scatra/nutrientconsumption)
  detFe_->at(gp) =   pow( + Cdach(0)*(Cdach(1)*Cdach(2)-Cdach(4)*Cdach(4))
                          - Cdach(3)*(Cdach(3)*Cdach(2)-Cdach(5)*Cdach(4))
                          + Cdach(5)*(Cdach(3)*Cdach(4)-Cdach(5)*Cdach(1)),
                          1/2);

  //store dtheta
  const double dt = params.get<double>("delta time",-1.0);
  for (unsigned i=0; i<theta_->size(); i++)
    (*dtheta_)[i] = ((*theta_)[i]-(*thetaold_)[i])/dt;
}

/*----------------------------------------------------------------------*
 |  Evaluate growth function                           (protected)        02/10|
 *----------------------------------------------------------------------*/
void MAT::GrowthScd::EvaluateGrowthFunction
(
    double & growthfunc, // (o)
    double traceM,       // (i)
    double theta         // (i)
)
{
  //call stress based growth law
  Growth::EvaluateGrowthFunction(growthfunc,traceM,theta);

  // parameters
  const double satcoeff = params_->satcoeff_;

  // scale with concentration dependent factor
 growthfunc = concentration_/(satcoeff+concentration_) * growthfunc;

  return;
}

/*----------------------------------------------------------------------*
 |  Evaluate derivative of growth function       (protected)        02/10|
 *----------------------------------------------------------------------*/
void MAT::GrowthScd::EvaluateGrowthFunctionDerivTheta
(
    double & dgrowthfunctheta,
    double traceM,
    double theta,
    const LINALG::Matrix<NUM_STRESS_3D, 1>& Cdach,
    const LINALG::Matrix<NUM_STRESS_3D, NUM_STRESS_3D>& cmatelastic
)
{
  //call stress based growth law
  Growth::EvaluateGrowthFunctionDerivTheta(dgrowthfunctheta,traceM,theta,Cdach,cmatelastic);

  // parameters
  const double satcoeff = params_->satcoeff_;

  // scale with concentration dependent factor
  dgrowthfunctheta = concentration_/(satcoeff+concentration_) * dgrowthfunctheta;

  return;
}

/*----------------------------------------------------------------------*
 |  Evaluate derivative of growth function       (protected)        02/10|
 *----------------------------------------------------------------------*/
void MAT::GrowthScd::EvaluateGrowthFunctionDerivC
(
    LINALG::Matrix<NUM_STRESS_3D, 1>& dgrowthfuncdC,
    double traceM,
    double theta,
    const LINALG::Matrix<NUM_STRESS_3D, 1>& C,
    const LINALG::Matrix<NUM_STRESS_3D, 1>& S,
    const LINALG::Matrix<NUM_STRESS_3D, NUM_STRESS_3D>& cmat
)
{
  //call stress based growth law
  Growth::EvaluateGrowthFunctionDerivC(dgrowthfuncdC,traceM,theta,C,S,cmat);

  // parameters
  const double satcoeff = params_->satcoeff_;

  // scale with concentration dependent factor
  dgrowthfuncdC.Scale(concentration_/(satcoeff+concentration_));

  return;
}

