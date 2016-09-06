/*!----------------------------------------------------------------------
\file myocard.cpp

\brief myocard material

<pre>
\level 3

\maintainer Julia Hoermann
            hoermann@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15264
</pre>
*/


/*----------------------------------------------------------------------*
 |  headers                                                  cbert 09/12 |
 *----------------------------------------------------------------------*/

#include <vector>
#include "myocard.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_mat/matpar_bundle.H"
#include "../drt_lib/drt_linedefinition.H"
#include <fstream> // For plotting ion concentrations
// Possible Cell models
#include "myocard_minimal.H"
#include "myocard_fitzhugh_nagumo.H"
#include "myocard_inada.H"
#include "myocard_tentusscher.H"
#include "myocard_san_garny.H"


/*----------------------------------------------------------------------*
 |                                                          cbert 09/12 |
 *----------------------------------------------------------------------*/
MAT::PAR::Myocard::Myocard( Teuchos::RCP<MAT::PAR::Material> matdata )
: Parameter(matdata),
  diff1(matdata->GetDouble("DIFF1")),
  diff2(matdata->GetDouble("DIFF2")),
  diff3(0.0),
  dt_deriv(matdata->GetDouble("PERTUBATION_DERIV")),
  model(matdata->Get<std::string>("MODEL")),
  tissue(matdata->Get<std::string>("TISSUE")),
  time_scale(matdata->GetDouble("TIME_SCALE")),
  num_gp(matdata->GetInt("GP"))
  {
  }

Teuchos::RCP<MAT::Material> MAT::PAR::Myocard::CreateMaterial()
  {
    return Teuchos::rcp(new MAT::Myocard(this));
  }


MAT::MyocardType MAT::MyocardType::instance_;


DRT::ParObject* MAT::MyocardType::Create( const std::vector<char> & data )
  {
    MAT::Myocard* myocard = new MAT::Myocard();
    myocard->Unpack(data);
    return myocard;
  }


/*----------------------------------------------------------------------*
 |  Constructor                                    (public)  cbert 08/13 |
 *----------------------------------------------------------------------*/
MAT::Myocard::Myocard()
  : params_(NULL),
  difftensor_(true),
  nb_state_variables_(0),
  myocard_mat_(Teuchos::null)
{

}


/*----------------------------------------------------------------------*
 |  Constructor                                   (public)   cbert 08/13 |
 *----------------------------------------------------------------------*/
MAT::Myocard::Myocard(MAT::PAR::Myocard* params)
  : params_(params),
  difftensor_(true)
{
  Initialize();
  nb_state_variables_ = myocard_mat_->GetNumberOfInternalStateVariables();

}


/*----------------------------------------------------------------------*
 |  Pack                                           (public)  cbert 09/12 |
 *----------------------------------------------------------------------*/
void MAT::Myocard::Pack(DRT::PackBuffer& data) const
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
  AddtoPack(data,nb_state_variables_);
  AddtoPack(data, difftensor_);

  // pack history data
  if (myocard_mat_ != Teuchos::null)
    {
      for(int k=-1; k<nb_state_variables_;++k) // Starting from -1 for mechanical activation
        {
          AddtoPack(data, myocard_mat_->GetInternalState(k));
        }
    }

  return;
}

/*----------------------------------------------------------------------*
 |  Unpack                                         (public)  cbert 09/12 |
 *----------------------------------------------------------------------*/
void MAT::Myocard::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  // extract type
  int type = 0;
  ExtractfromPack(position,data,type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");

    // matid
  int matid;
  ExtractfromPack(position,data,matid);
  ExtractfromPack(position, data, nb_state_variables_);
  ExtractfromPack(position, data, difftensor_);
  params_ = NULL;
  if (DRT::Problem::Instance()->Materials() != Teuchos::null)  // it does not enter here in postprocessing
  {
    if (DRT::Problem::Instance()->Materials()->Num() != 0)
    {
      const int probinst = DRT::Problem::Instance()->Materials()->GetReadFromProblem();
      MAT::PAR::Parameter* mat = DRT::Problem::Instance(probinst)->Materials()->ParameterById(matid);
      if (mat->Type() == MaterialType())
      {
        params_ = static_cast<MAT::PAR::Myocard*>(mat);
      }
      else
        dserror("Type of parameter material %d does not fit to calling type %d", mat->Type(), MaterialType());

      // Initialize material
      Initialize();

      // unpack history data
     double val;
     for(int k=-1; k<nb_state_variables_;++k) // Starting from -1 for mechanical activation
     {
       ExtractfromPack(position, data, val);
       myocard_mat_->SetInternalState(k,val);
     }

     if (position != data.size())
     dserror("Mismatch in size of data %d <-> %d",data.size(),position);

    }
  }


}


/*----------------------------------------------------------------------*
 |  Setup conductivity tensor                                cbert 02/13 |
 *----------------------------------------------------------------------*/
void MAT::Myocard::Setup(const LINALG::Matrix<3,1>& fiber1)
{
  SetupDiffusionTensor(fiber1);
}

void MAT::Myocard::Setup(const LINALG::Matrix<2,1>& fiber1)
{
  SetupDiffusionTensor(fiber1);
}

void MAT::Myocard::Setup(DRT::INPUT::LineDefinition* linedef)
{
  std::vector<double> fiber1(3);
  linedef->ExtractDoubleVector("FIBER1",fiber1);
  SetupDiffusionTensor(fiber1);
}


void MAT::Myocard::SetupDiffusionTensor(const std::vector<double> &fiber1)
  {

  // Normalize fiber1
  double fiber1normS = fiber1[0]*fiber1[0]+fiber1[1]*fiber1[1]+fiber1[2]*fiber1[2];

  // get conductivity values of main fiber direction and perpendicular to fiber direction (rot symmetry)
   const double diff1 = params_->diff1;
   const double diff2  = params_->diff2;

   // ******** SETUP ORTHOTROPIC DIFFUSION TENSOR: diff2*Id + (diff1-diff2)*fiber1*fiber1'
   // first row
   difftensor_(0,0)=diff2 + (diff1-diff2)*fiber1[0]*fiber1[0]/fiber1normS;
   difftensor_(0,1)=(diff1-diff2)*fiber1[0]*fiber1[1]/fiber1normS;
   difftensor_(0,2)=(diff1-diff2)*fiber1[0]*fiber1[2]/fiber1normS;
   // second row
   difftensor_(1,0)=(diff1-diff2)*fiber1[1]*fiber1[0]/fiber1normS;
   difftensor_(1,1)=diff2 + (diff1-diff2)*fiber1[1]*fiber1[1]/fiber1normS;
   difftensor_(1,2)=(diff1-diff2)*fiber1[1]*fiber1[2]/fiber1normS;
   // third row
   difftensor_(2,0)=(diff1-diff2)*fiber1[2]*fiber1[0]/fiber1normS;
   difftensor_(2,1)=(diff1-diff2)*fiber1[2]*fiber1[1]/fiber1normS;
   difftensor_(2,2)=diff2 + (diff1-diff2)*fiber1[2]*fiber1[2]/fiber1normS;
   // done
   return;

  }

void MAT::Myocard::SetupDiffusionTensor(const LINALG::Matrix<3,1>& fiber1)
{

// Normalize fiber1
double fiber1normS = fiber1(0)*fiber1(0)+fiber1(1)*fiber1(1)+fiber1(2)*fiber1(2);

// get conductivity values of main fiber direction and perpendicular to fiber direction (rot symmetry)
 const double diff1 = params_->diff1;
 const double diff2  = params_->diff2;

 // ******** SETUP ORTHOTROPIC DIFFUSION TENSOR: diff2*Id + (diff1-diff2)*fiber1*fiber1'
 // first row
 difftensor_(0,0)=diff2 + (diff1-diff2)*fiber1(0)*fiber1(0)/fiber1normS;
 difftensor_(0,1)=(diff1-diff2)*fiber1(0)*fiber1(1)/fiber1normS;
 difftensor_(0,2)=(diff1-diff2)*fiber1(0)*fiber1(2)/fiber1normS;
 // second row
 difftensor_(1,0)=(diff1-diff2)*fiber1(1)*fiber1(0)/fiber1normS;
 difftensor_(1,1)=diff2 + (diff1-diff2)*fiber1(1)*fiber1(1)/fiber1normS;
 difftensor_(1,2)=(diff1-diff2)*fiber1(1)*fiber1(2)/fiber1normS;
 // third row
 difftensor_(2,0)=(diff1-diff2)*fiber1(2)*fiber1(0)/fiber1normS;
 difftensor_(2,1)=(diff1-diff2)*fiber1(2)*fiber1(1)/fiber1normS;
 difftensor_(2,2)=diff2 + (diff1-diff2)*fiber1(2)*fiber1(2)/fiber1normS;
 // done
 return;

}

void MAT::Myocard::SetupDiffusionTensor(const LINALG::Matrix<2,1>& fiber1)
{

// Normalize fiber1
double fiber1normS = fiber1(0)*fiber1(0)+fiber1(1)*fiber1(1);

// get conductivity values of main fiber direction and perpendicular to fiber direction (rot symmetry)
 const double diff1 = params_->diff1;
 const double diff2  = params_->diff2;

 // ******** SETUP ORTHOTROPIC DIFFUSION TENSOR: diff2*Id + (diff1-diff2)*fiber1*fiber1'
 // first row
 difftensor_(0,0)=diff2 + (diff1-diff2)*fiber1(0)*fiber1(0)/fiber1normS;
 difftensor_(0,1)=(diff1-diff2)*fiber1(0)*fiber1(1)/fiber1normS;
 // second row
 difftensor_(1,0)=(diff1-diff2)*fiber1(1)*fiber1(0)/fiber1normS;
 difftensor_(1,1)=diff2 + (diff1-diff2)*fiber1(1)*fiber1(1)/fiber1normS;
 // done
 return;

}

void MAT::Myocard::Diffusivity(LINALG::Matrix<1,1>& diffus3) const
{
  diffus3(0,0) = difftensor_(0,0); return;
}

void MAT::Myocard::Diffusivity(LINALG::Matrix<2,2>& diffus3) const
{
  for (int i=0; i<2; i++){
    for (int j=0; j<2; j++){
      diffus3(i,j) = difftensor_(i,j);
    }
  }

  return;
}

void MAT::Myocard::Diffusivity(LINALG::Matrix<3,3>& diffus3) const
{
  for (int i=0; i<3; i++){
    for (int j=0; j<3; j++){
      diffus3(i,j) = difftensor_(i,j);
    }
  }
  return;
}

/*----------------------------------------------------------------------*
 |                                                           cbert 09/13 |
 *----------------------------------------------------------------------*/
double MAT::Myocard::ReaCoeff(const double phi, const double dt) const
{
  double reacoeff = params_->time_scale;
  reacoeff *= myocard_mat_->ReaCoeff(phi,dt*params_->time_scale);

  return reacoeff;
}

/*----------------------------------------------------------------------*
 |                                                       hoermann 09/15 |
 *----------------------------------------------------------------------*/
double MAT::Myocard::ReaCoeff(const double phi, const double dt, int gp) const
{
  double reacoeff = params_->time_scale;
  if(gp == -1)
    reacoeff *= myocard_mat_->ReaCoeff(phi,dt*params_->time_scale);
  else
    reacoeff *= myocard_mat_->ReaCoeff(phi,dt*params_->time_scale, gp);

  return reacoeff;
}


/*----------------------------------------------------------------------*
 |                                                           cbert 09/13 |
 *----------------------------------------------------------------------*/
double MAT::Myocard::ReaCoeffDeriv(const double phi, const double dt) const
{
  double ReaCoeffDeriv=0.0;
  if(params_->dt_deriv != 0.0)
  {
    double ReaCoeff_t2 = ReaCoeff((phi+params_->dt_deriv), dt);
    double ReaCoeff_t1 = ReaCoeff(phi, dt);
    ReaCoeffDeriv = (ReaCoeff_t2 - ReaCoeff_t1)/(params_->dt_deriv);
  }
  return ReaCoeffDeriv;
}


/*----------------------------------------------------------------------*
 |                                                       hoermann 09/15 |
 *----------------------------------------------------------------------*/
double MAT::Myocard::ReaCoeffDeriv(const double phi, const double dt, int gp) const
{
  double ReaCoeffDeriv=0.0;
  if(params_->dt_deriv != 0.0)
  {
    double ReaCoeff_t2 = ReaCoeff((phi+params_->dt_deriv), dt, gp);
    double ReaCoeff_t1 = ReaCoeff(phi, dt, gp);
    ReaCoeffDeriv = (ReaCoeff_t2 - ReaCoeff_t1)/(params_->dt_deriv);
  }
  return ReaCoeffDeriv;
}


/*----------------------------------------------------------------------*
 |  returns number of internal state variables              cbert 08/13 |
 *----------------------------------------------------------------------*/
int MAT::Myocard::GetNumberOfInternalStateVariables() const
{
  int val=0;
  val = myocard_mat_->GetNumberOfInternalStateVariables();
  return val;
}

/*----------------------------------------------------------------------*
 |  returns current internal state of the material          cbert 08/13 |
 *----------------------------------------------------------------------*/
double MAT::Myocard::GetInternalState(const int k) const
{
  double val=0.0;
  val = myocard_mat_->GetInternalState(k);
  return val;
}

/*----------------------------------------------------------------------*
 |  returns current internal state of the material       hoermann 09/15 |
 |  for multiple points per element                                     |
 *----------------------------------------------------------------------*/
double MAT::Myocard::GetInternalState(const int k, int gp) const
{
  double val=0.0;
  val = myocard_mat_->GetInternalState(k, gp);
  return val;
}

/*----------------------------------------------------------------------*
 |  returns current internal state of the material          cbert 08/13 |
 *----------------------------------------------------------------------*/
void MAT::Myocard::SetInternalState(const int k, const double val)
{
  myocard_mat_->SetInternalState(k,val);
}

/*----------------------------------------------------------------------*
 |  returns current internal state of the material       hoermann 09/15 |
 |  for multiple points per element                                     |
 *----------------------------------------------------------------------*/
void MAT::Myocard::SetInternalState(const int k, const double val, int gp)
{
  myocard_mat_->SetInternalState(k,val,gp);
}

/*----------------------------------------------------------------------*
 |  returns number of internal state variables of the material  cbert 08/13 |
 *----------------------------------------------------------------------*/
int MAT::Myocard::GetNumberOfIonicCurrents() const
{
  int val=0;
  val = myocard_mat_->GetNumberOfIonicCurrents();
  return val;
}

/*----------------------------------------------------------------------*
 |  returns current internal currents                    hoermann 09/15 |
 |  for multiple points per element                                     |
 *----------------------------------------------------------------------*/
double MAT::Myocard::GetIonicCurrents(const int k, int gp) const
{
  double val=0.0;
  val = myocard_mat_->GetIonicCurrents(k, gp);
  return val;
}

/*----------------------------------------------------------------------*
 |  returns current internal currents                       cbert 08/13 |
 *----------------------------------------------------------------------*/
double MAT::Myocard::GetIonicCurrents(const int k) const
{
  double val=0.0;
  val = myocard_mat_->GetIonicCurrents(k);
  return val;
}


/*----------------------------------------------------------------------*
 |  initialize internal variables (called by constructors)   cbert 09/12 |
 *----------------------------------------------------------------------*/
void MAT::Myocard::Initialize()
{
  if (*(params_->model) == "MV") myocard_mat_ = Teuchos::rcp(new Myocard_Minimal(params_->dt_deriv,*(params_->tissue),params_->num_gp));
  else if (*(params_->model) == "FHN") myocard_mat_ = Teuchos::rcp(new Myocard_Fitzhugh_Nagumo(params_->dt_deriv,*(params_->tissue),params_->num_gp));
  else if (*(params_->model) == "INADA") myocard_mat_ = Teuchos::rcp(new Myocard_Inada(params_->dt_deriv,*(params_->tissue)));
  else if (*(params_->model) == "TNNP") myocard_mat_ = Teuchos::rcp(new Myocard_TenTusscher(params_->dt_deriv,*(params_->tissue)));
  else if (*(params_->model) == "SAN") myocard_mat_ = Teuchos::rcp(new Myocard_SAN_Garny(params_->dt_deriv,*(params_->tissue)));
  else dserror("Myocard Material type is not supported! (for the moment only MV,FHN,INADA and TNNP)");

  return;
}

/*----------------------------------------------------------------------*
 |  update of material at the end of a time step             ljag 07/12 |
 *----------------------------------------------------------------------*/
void MAT::Myocard::Update(const double phi, const double dt)
{
  myocard_mat_->Update(phi, dt*(params_->time_scale));

  return;
}

