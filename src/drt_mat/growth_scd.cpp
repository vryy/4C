/*----------------------------------------------------------------------*/
/*!
 \file growth_scd.cpp

 \brief

This file contains routines for an integration point based and scalar dependend volumetric growth law.
It is derived from the growth laws implemented in growth_ip.cpp and additional adds the scalar dependency, e.g. nutrients.

 <pre>
   Maintainer: Moritz Thon
               thon@lnm.mw.tum.de
               http://www.mhpc.mw.tum.de
               089 - 289-10364
 </pre>
 *----------------------------------------------------------------------*/


#include "growth_scd.H"
#include "growth_law.H"
#include "../drt_mat/matpar_bundle.H"
#include "../drt_lib/drt_globalproblem.H"   // for function Factory in Unpack
#include "../drt_lib/drt_utils_factory.H"   // for function Factory in Unpack


/*----------------------------------------------------------------------*
 |                                                           vuong 06/11  |
 *----------------------------------------------------------------------*/
MAT::PAR::GrowthScd::GrowthScd(
  Teuchos::RCP<MAT::PAR::Material> matdata
  )
: Growth(matdata),
  rearate_(matdata->GetDouble("REARATE")),
  satcoeff_(matdata->GetDouble("SATCOEFF"))
{
  if (growthlaw_->MaterialType() == INPAR::MAT::m_growth_linear or growthlaw_->MaterialType() == INPAR::MAT::m_growth_exponential)
  {
    if (rearate_ <= 0)
        dserror("You need to choose a positive reaction rate!");
    if (satcoeff_ < 0)
        dserror("You need to choose a non-negative saturation coefficient!");
  }
}


Teuchos::RCP<MAT::Material> MAT::PAR::GrowthScd::CreateMaterial()
{
  Teuchos::RCP<MAT::Material> mat;

  switch (growthlaw_->MaterialType())
  {
  case INPAR::MAT::m_growth_linear:
  case INPAR::MAT::m_growth_exponential:
  case INPAR::MAT::m_growth_biofilm:
    mat = Teuchos::rcp(new MAT::GrowthScd(this));
    break;
  case INPAR::MAT::m_growth_ac:
    mat = Teuchos::rcp(new MAT::GrowthScdAC(this));
    break;
  case INPAR::MAT::m_growth_ac_radial:
  case INPAR::MAT::m_growth_ac_radial_refconc:
    mat = Teuchos::rcp(new MAT::GrowthScdACRadial(this));
    break;
  default:
    dserror("The growth law you have chosen is not valid one!");
    mat = Teuchos::null;
    break;
  }

  return mat;
}


MAT::GrowthScdType MAT::GrowthScdType::instance_;

DRT::ParObject* MAT::GrowthScdType::Create( const std::vector<char> & data )
{
  MAT::GrowthScd* grow = new MAT::GrowthScd();
  grow->Unpack(data);
  return grow;
}

/*----------------------------------------------------------------------*
 |  Constructor                                   (public)  vuong 02/10|
 *----------------------------------------------------------------------*/
MAT::GrowthScd::GrowthScd()
  : GrowthMandel(),
    detFe_(Teuchos::null), //initialized in GrowthScd::Unpack
    dtheta_(Teuchos::null), //initialized in GrowthScd::Unpack
    paramsScd_(NULL)
{
}


/*----------------------------------------------------------------------*
 |  Copy-Constructor                             (public)    vuong  02/10|
 *----------------------------------------------------------------------*/
MAT::GrowthScd::GrowthScd(MAT::PAR::GrowthScd* params)
  : GrowthMandel(params),
    detFe_(Teuchos::null), //initialized in GrowthScd::Setup
    dtheta_(Teuchos::null), //initialized in GrowthScd::Setup
    paramsScd_(params)
{
}

/*----------------------------------------------------------------------*
 |  ResetAll                                      (public)         11/12|
 *----------------------------------------------------------------------*/
void MAT::GrowthScd::ResetAll(const int numgp)
{
  for (int j=0; j<numgp; ++j)
  {
    detFe_->at(j) = 1.0;
    dtheta_->at(j) = 0.0;
  }

  MAT::GrowthMandel::ResetAll(numgp);
}

/*----------------------------------------------------------------------*
 |  Pack                                          (public)    vuong  02/10|
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

  MAT::PAR::GrowthScd* params=Parameter();
  if (params != NULL) matid = params->Id();  // in case we are in post-process mode
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
  GrowthMandel::Pack(data);

  return;
}


/*----------------------------------------------------------------------*
 |  Unpack                                        (public)   vuong 02/10|
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

  paramsScd_ = NULL;
  if (DRT::Problem::Instance()->Materials() != Teuchos::null)
    if (DRT::Problem::Instance()->Materials()->Num() != 0)
    {
      const int probinst = DRT::Problem::Instance()->Materials()->GetReadFromProblem();
      MAT::PAR::Parameter* mat = DRT::Problem::Instance(probinst)->Materials()->ParameterById(matid);
      if (mat->Type() == MaterialType())
        paramsScd_ = static_cast<MAT::PAR::GrowthScd*>(mat);
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
  GrowthMandel::ExtractfromPack(position,data,basedata);
  GrowthMandel::Unpack(basedata);

  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d",data.size(),position);

  return;
}

/*----------------------------------------------------------------------*
 |  Setup                                         (public) vuong  02/10|
 *----------------------------------------------------------------------*/
void MAT::GrowthScd::Setup(int numgp, DRT::INPUT::LineDefinition* linedef)
{
  detFe_ = Teuchos::rcp(new std::vector<double> (numgp,1.0));
  dtheta_ = Teuchos::rcp(new std::vector<double> (numgp,0.0));

  //setup base class
  GrowthMandel::Setup(numgp, linedef);
  return;
}


/*----------------------------------------------------------------------*
 |  Evaluate Material                             (public)  vuong   02/10|
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

  GrowthMandel::Evaluate(defgrd,glstrain,params,stress,cmat,eleGID);

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
    (*dtheta_)[i] = ((*theta_)[i]-(*ThetaOld())[i])/dt;
}

// evaluate the volumetric growth factor theta
void MAT::GrowthScd::EvaluateGrowth(double* theta,
                      LINALG::Matrix<6,1>* dthetadC,
                      const LINALG::Matrix<3,3>* defgrd,
                      const LINALG::Matrix<6,1>* glstrain,
                      Teuchos::ParameterList& params,
                      const int eleGID )
{

  // parameters
  MAT::PAR::GrowthScd* matparams=Parameter();
  const double rearate = matparams->rearate_;
  const double satcoeff = matparams->satcoeff_;

  // get gauss point number
  const int gp = params.get<int>("gp",-1);
  if (gp == -1)
    dserror("no Gauss point number provided in material");

  //get pointer to vector containing the scalar values at the gauß points
  Teuchos::RCP<std::vector<std::vector<double> > > concentrations=
      params.get< Teuchos::RCP<std::vector<std::vector<double> > > >("gp_conc");

  //NOTE: in this growth law we always assume the first scalar to induce growth!
  double concentration_first=(concentrations->at(gp)).at(0);

  // concentration dependent factor
  double fac = rearate * concentration_first/(satcoeff+concentration_first);

  Parameter()->growthlaw_->SetFactor(fac);
  GrowthMandel::EvaluateGrowth(theta,dthetadC,defgrd,glstrain,params,eleGID);

  return;
}

MAT::GrowthScdACType MAT::GrowthScdACType::instance_;

DRT::ParObject* MAT::GrowthScdACType::Create( const std::vector<char> & data )
{
  MAT::GrowthScdAC* grow = new MAT::GrowthScdAC();
  grow->Unpack(data);
  return grow;
}

/*----------------------------------------------------------------------*
 |  Constructor                                               Thon 11/14|
 *----------------------------------------------------------------------*/
MAT::GrowthScdAC::GrowthScdAC()
  : GrowthMandel(),
//    concentrations_(Teuchos::null),
    paramsScdAC_(NULL)
{
}


/*----------------------------------------------------------------------*
 |  Copy-Constructor                                          Thon 11/14|
 *----------------------------------------------------------------------*/
MAT::GrowthScdAC::GrowthScdAC(MAT::PAR::GrowthScd* params)
  : GrowthMandel(params),
    paramsScdAC_(params)
{
}

/*----------------------------------------------------------------------*
 |  Pack                                                      Thon 11/14|
 *----------------------------------------------------------------------*/
void MAT::GrowthScdAC::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm( data );
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);
  // matid
  int matid = -1;

  MAT::PAR::GrowthScd* params=Parameter();
  if (params != NULL) matid = params->Id();  // in case we are in post-process mode
  AddtoPack(data,matid);

  // Pack base class material
  GrowthMandel::Pack(data);

  return;
}

/*----------------------------------------------------------------------*
 |  Unpack                                                    Thon 11/14|
 *----------------------------------------------------------------------*/
void MAT::GrowthScdAC::Unpack(const std::vector<char>& data)
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

  paramsScdAC_=NULL;
  if (DRT::Problem::Instance()->Materials() != Teuchos::null)
    if (DRT::Problem::Instance()->Materials()->Num() != 0)
    {
      const int probinst = DRT::Problem::Instance()->Materials()->GetReadFromProblem();
      MAT::PAR::Parameter* mat = DRT::Problem::Instance(probinst)->Materials()->ParameterById(matid);
      if (mat->Type() == MaterialType())
        paramsScdAC_ = dynamic_cast<MAT::PAR::GrowthScd*>(mat);
      else
        dserror("Type of parameter material %d does not fit to calling type %d", mat->Type(), MaterialType());
    }

  // extract base class material
  std::vector<char> basedata(0);
  GrowthMandel::ExtractfromPack(position,data,basedata);
  GrowthMandel::Unpack(basedata);


  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d",data.size(),position);

  return;
}

/*----------------------------------------------------------------------*
 | Evaluate the volumetric growth factor theta               Thon 11/14|
 *----------------------------------------------------------------------*/
void MAT::GrowthScdAC::EvaluateGrowth(double* theta,
                      LINALG::Matrix<6,1>* dthetadC,
                      const LINALG::Matrix<3,3>* defgrd,
                      const LINALG::Matrix<6,1>* glstrain,
                      Teuchos::ParameterList& params,
                      const int eleGID )
{
  // get gauss point number
  const int gp = params.get<int>("gp",-1);
  if (gp == -1)
    dserror("no Gauss point number provided in material");

  //get pointer to vector containing the scalar values at the gauß points
  Teuchos::RCP<std::vector<std::vector<double> > > concentrations=
      params.get< Teuchos::RCP<std::vector<std::vector<double> > > >("gp_conc");

  Parameter()->growthlaw_->SetFactor( concentrations->at(gp) );

  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////

  GrowthMandel::EvaluateGrowth(theta,dthetadC,defgrd,glstrain,params,eleGID);

  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
}

/*----------------------------------------------------------------------------*/
void MAT::GrowthScdAC::VisNames(std::map<std::string,int>& names)
{
  std::string fiber = "theta";
  names[fiber] = 1;

  Matelastic()->VisNames(names);
}

/*----------------------------------------------------------------------------*/
bool MAT::GrowthScdAC::VisData(const std::string& name, std::vector<double>& data, int numgp , int eleID)
{
  if (name == "theta")
  {
    if ((int)data.size()!=1)
      dserror("size mismatch");
    double temp = 0.0;
    for (int iter=0; iter<numgp; iter++)
      temp += theta_->at(iter);
    data[0] = temp/numgp;
  }
  else
  {
    return Matelastic()->VisData(name, data, numgp, eleID);
  }
  return true;
}

MAT::GrowthScdACRadialType MAT::GrowthScdACRadialType::instance_;

DRT::ParObject* MAT::GrowthScdACRadialType::Create( const std::vector<char> & data )
{
  MAT::GrowthScdACRadial* grow = new MAT::GrowthScdACRadial();
  grow->Unpack(data);
  return grow;
}

/*----------------------------------------------------------------------*
 |  Constructor                                               Thon 11/14|
 *----------------------------------------------------------------------*/
MAT::GrowthScdACRadial::GrowthScdACRadial()
  : GrowthScdAC(),
    N_(true),
    n_(Teuchos::null),
    n_for_update_(Teuchos::null),
    F_g_hist_(Teuchos::null)
{

}


/*----------------------------------------------------------------------*
 |  Copy-Constructor                                          Thon 11/14|
 *----------------------------------------------------------------------*/
MAT::GrowthScdACRadial::GrowthScdACRadial(MAT::PAR::GrowthScd* params)
  : GrowthScdAC(params),
    N_(true),
    n_(Teuchos::null),
    n_for_update_(Teuchos::null),
    F_g_hist_(Teuchos::null)
{

}

/*----------------------------------------------------------------------*
 |  Setup                                                     Thon 11/14|
 *----------------------------------------------------------------------*/
void MAT::GrowthScdACRadial::Setup(int numgp, DRT::INPUT::LineDefinition* linedef)
{
  // CIR-AXI-RAD nomenclature
  if (not (linedef->HaveNamed("RAD")))
    dserror("If you want growth into the radial direction you need to specify RAD in your input file!");

  ReadFiber(linedef,"RAD",N_);
  n_=std::vector<LINALG::Matrix<3,1> >(numgp,N_);
  n_for_update_=std::vector<LINALG::Matrix<3,1> >(numgp,N_);

  LINALG::Matrix<3,3> Id(true);
  for (int i = 0; i < 3; i++)
    Id(i,i) = 1.0;
  F_g_hist_ = std::vector<LINALG::Matrix<3,3> >(numgp,Id);

  GrowthScdAC::Setup(numgp, linedef);
  return;
}

/*----------------------------------------------------------------------*
 | Function which reads in the given fiber value due to the             |
 | FIBER1 nomenclature                                       Thon 01/15 |
 *----------------------------------------------------------------------*/
void MAT::GrowthScdACRadial::ReadFiber(
    DRT::INPUT::LineDefinition* linedef,
    std::string specifier,
    LINALG::Matrix<3,1> &fiber_vector
)
{
  std::vector<double> fiber1;
  linedef->ExtractDoubleVector(specifier,fiber1);
  double f1norm=0.;
  //normalization
  for (int i = 0; i < 3; ++i)
  {
    f1norm += fiber1[i]*fiber1[i];
  }
  f1norm = sqrt(f1norm);

  // fill final fiber vector
  for (int i = 0; i < 3; ++i)
    fiber_vector(i) = fiber1[i]/f1norm;
}

/*----------------------------------------------------------------------*
 |  Pack                                                      Thon 11/14|
 *----------------------------------------------------------------------*/
void MAT::GrowthScdACRadial::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm( data );
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);

  // Pack internal variables
  for (int i = 0; i < 3; ++i)
  {
    AddtoPack(data,N_(i,0));
  }

  int numgp=0;
  if (isinit_)
  {
    numgp = theta_->size();;   // size is number of gausspoints
  }
  AddtoPack(data,numgp);

  for (int gp=0; gp<numgp; gp++)
  {
    LINALG::Matrix<3,3> F_g_hist = F_g_hist_.at(gp);
    LINALG::Matrix<3,1> n = n_.at(gp);
    LINALG::Matrix<3,1> n_for_update = n_for_update_.at(gp);

    for (int i = 0; i < 3; ++i)
    {
      for (int j = 0; j < 3; ++j)
      {
        AddtoPack(data,F_g_hist(i,j));
      }
      AddtoPack(data,n(i,0));
      AddtoPack(data,n_for_update(i,0));
    }
  }

  // Pack base class material
  GrowthScdAC::Pack(data);

  return;
}

/*----------------------------------------------------------------------*
 |  Unpack                                                    Thon 11/14|
 *----------------------------------------------------------------------*/
void MAT::GrowthScdACRadial::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position,data,type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");

  // unpack growth internal variables
  // Pack internal variables
  for (int i = 0; i < 3; ++i)
  {
    double Ni;
    ExtractfromPack(position,data,Ni);
    N_(i,0) =Ni;
  }

  int numgp;
  ExtractfromPack(position,data,numgp);
  if (not (numgp == 0))
  {
    F_g_hist_ = std::vector<LINALG::Matrix<3,3> >(numgp,LINALG::Matrix<3,3>(true));
    n_ = std::vector<LINALG::Matrix<3,1> >(numgp,LINALG::Matrix<3,1>(true));
    n_for_update_ = std::vector<LINALG::Matrix<3,1> >(numgp,LINALG::Matrix<3,1>(true));

    for (int gp=0; gp<numgp; gp++)
    {
      LINALG::Matrix<3,3> F_g_hist(true);
      LINALG::Matrix<3,1> n(true);
      LINALG::Matrix<3,1> n_for_update(true);

      for (int i = 0; i < 3; ++i)
      {
        for (int j = 0; j < 3; ++j)
        {
          double F_g_hist_ij;
          ExtractfromPack(position,data,F_g_hist_ij);
          F_g_hist(i,j) =F_g_hist_ij;
        }
        double n_i;
        ExtractfromPack(position,data,n_i);
        n(i,0) =n_i;

        double n_for_update_i;
        ExtractfromPack(position,data,n_for_update_i);
        n_for_update(i,0) =n_for_update_i;
      }
      F_g_hist_.at(gp)=F_g_hist;
      n_.at(gp)=n;
      n_for_update_.at(gp)=n_for_update;
    }
  }

  // extract base class material
  std::vector<char> basedata(0);
  GrowthScdAC::ExtractfromPack(position,data,basedata);
  GrowthScdAC::Unpack(basedata);


  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d",data.size(),position);

  return;
}

/*----------------------------------------------------------------------*
|  Time update of material                                   Thon 01/15|
*----------------------------------------------------------------------*/
void MAT::GrowthScdACRadial::Update()
{
  const int numgp = theta_->size();

  for (int gp=0; gp<numgp; gp++)
  {
    LINALG::Matrix<3,3> F_g_hist_new(true);

    CalcFg(theta_->at(gp),gp,F_g_hist_new);

    F_g_hist_.at(gp) = F_g_hist_new;

    n_.at(gp)=n_for_update_.at(gp);
  }

  GrowthScdAC::Update();
}

/*----------------------------------------------------------------------*
 |  Evaluate Material                                         Thon 01/15|
 *----------------------------------------------------------------------*/
void MAT::GrowthScdACRadial::Evaluate
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
  const int gp = params.get<int>("gp", -1);
  if (gp == -1)
    dserror("no Gauss point number provided in material");

  double time = params.get<double>("total time", -1.0);
  if( abs(time+1.0) < 1e-14 ) dserror("no time step or no total time given for growth material!");
  std::string action = params.get<std::string>("action", "none");
  bool output = false;
  if (action == "calc_struct_stress")
    output = true;

  const double eps = 1.0e-14;
  MAT::PAR::Growth* growth_params = Parameter();
  const double endtime = growth_params->endtime_;
  const double starttime = growth_params->starttime_;
  // when stress output is calculated the final parameters already exist
  // we should not do another local Newton iteration, which uses eventually a wrong thetaold
  if (output)
    time = endtime + 1.0;

  if (time > starttime + eps && time <= endtime + eps) //iff growth is active
  {
    //push-forward of N
    LINALG::Matrix<3,3> defgrdinv(true);
    defgrdinv.Invert(*defgrd);
    LINALG::Matrix<3,1> n_for_update(true);
    n_for_update.MultiplyTN(defgrd->Determinant(),defgrdinv,N_);
    //scale n to length of one
    n_for_update.Scale(1.0/n_for_update.Norm2());
    //save for time update
    n_for_update_.at(gp)=n_for_update;

//    n_.at(gp)=n_for_update;

    //--------------------------------------------------------------------------------------
    // calculate \theta and  \frac{\partial \theta}{\partial C}
    //--------------------------------------------------------------------------------------
    double theta;
    LINALG::Matrix<6,1> dthetadCvec(true);
    EvaluateGrowth(&theta,&dthetadCvec,defgrd,glstrain,params,eleGID);

    LINALG::Matrix<6,1> S(true);
    LINALG::Matrix<6,6> cmatdach(true);

    GetSAndCmatdach(theta,defgrd,&S,&cmatdach,params,eleGID);

    *stress=S;

    //--------------------------------------------------------------------------------------
    // calculate elastic material stiffness matrix = \frac{\partial S}{\partial E}
    //--------------------------------------------------------------------------------------

    // calculate growth part F_g of the deformation gradient F
    LINALG::Matrix<3,3> F_g(true);
    CalcFg(theta,gp,F_g);

    // calculate F_g^(-1)
    LINALG::Matrix<3,3> F_ginv(true);
    F_ginv.Invert(F_g);

    // constitutive matrix including growth cmat = F_g^-1 F_g^-1 cmatdach F_g^-T F_g^-T
    LINALG::Matrix<NUM_STRESS_3D, NUM_STRESS_3D> cmatelast(true);

    cmatelast = PullBack4Tensor(F_ginv,cmatdach);

    if (params.isParameter("dconc_zero_ratio_dC"))
    {
      //if the growth law shall be proportional to the scalar in the
      //spatial configuration on has to a lot more :-(

      //--------------------------------------------------------------------------------------
      // calculate \frac{\partial S}{\partial \theta}
      //--------------------------------------------------------------------------------------
      // NOTE: we do this by a FD approximation, which is really cheap here due to the fact
      // that theta is a scalar value (hence only one more material evaluation is needed!)

      const double espilon = 1.0e-8;

      LINALG::Matrix<6,1> SEps(true);
      LINALG::Matrix<6,6> cmatdachEps(true);

      GetSAndCmatdach(theta+espilon,defgrd,&SEps,&cmatdachEps,params,eleGID);

      //--------------------------------------------------------------------------------------
      // calculate \frac{d S}{d E} = \frac{\partial S}{\partial E} +
      //    + 2* \left( \frac{\partial S}{\partial \theta} \times \frac{\partial \theta}{\partial C} \right)
      //--------------------------------------------------------------------------------------
      for (int i = 0; i < 6; i++)
      {
        for (int j = 0; j < 6; j++)
        {

          (*cmat)(i, j) =  cmatelast(i, j) + 2.0 *(SEps(i) - S(i))/espilon * dthetadCvec(j);
        }
      }

      //clean up for safety reasons
      params.remove("dconc_zero_ratio_dC");
      params.remove("conc_zero_ratio");
    }
    else //dthetadCvec is zero
    {
      //Note: we want to end up here only if our theta does not depend on structural variables, hence
      // iff dthetadCvec is zero. We better check:
      if (dthetadCvec.NormInf() > 1e-14)
        dserror("dthetadCvec should be zero in this case!");

      //we are done!
      *cmat=cmatelast;

      //clean up for safety reasons
      params.remove("conc_zero_ratio");
    }

    // store theta
    theta_->at(gp) = theta;
  }
  else if (time > endtime + eps)
  {
    double theta = theta_->at(gp);

    LINALG::Matrix<6,1> S(true);
    LINALG::Matrix<6,6> cmatdach(true);

    GetSAndCmatdach(theta,defgrd,&S,&cmatdach,params,eleGID);

    *stress=S;

    //--------------------------------------------------------------------------------------
    // calculate material stiffness matrix = dS/dE
    //--------------------------------------------------------------------------------------

    // calculate growth part F_g of the deformation gradient F
    LINALG::Matrix<3,3> F_g(true);
    CalcFg(theta,gp,F_g);

    // calculate F_g^(-1)
    LINALG::Matrix<3,3> F_ginv(true);
    F_ginv.Invert(F_g);

    // constitutive matrix including growth cmat = F_g^-1 F_g^-1 cmatdach F_g^-T F_g^-T
    LINALG::Matrix<NUM_STRESS_3D, NUM_STRESS_3D> cmatelast(true);

    cmatelast = PullBack4Tensor(F_ginv,cmatdach);

    *cmat = cmatelast;
  }
  else //no growth has happened jet
  {
    Matelastic()->Evaluate(defgrd, glstrain, params, stress, cmat, eleGID); //evaluate the standard material
  }

}

///*----------------------------------------------------------------------*
// | calculate growth part of deformation gradient              Thon 01/16|
// *----------------------------------------------------------------------*/
void MAT::GrowthScdACRadial::CalcFg( const double& theta, const int& gp, LINALG::Matrix<3,3>& F_g )
{
  LINALG::Matrix<3,3> ndn(true);
  ndn.MultiplyNT(n_.at(gp),n_.at(gp));

  LINALG::Matrix<3,3> F_g_incr(true);
  for (int i = 0; i < 3; i++)
    F_g_incr(i,i) = 1.0;
  F_g_incr.Update((theta-ThetaOld()->at(gp))/ThetaOld()->at(gp),ndn,1.0);

  F_g.MultiplyNN(F_g_incr,F_g_hist_.at(gp));
}

///*----------------------------------------------------------------------*
// | calculate stresses and elastic material tangent                      |
// | (both in Voigt notation)                                   Thon 01/16|
// *----------------------------------------------------------------------*/
void MAT::GrowthScdACRadial::GetSAndCmatdach
(
    const double theta,
    const LINALG::Matrix<3,3>* defgrd,
    LINALG::Matrix<6,1>* stress,
    LINALG::Matrix<6,6>* cmatdach,
    Teuchos::ParameterList& params,
    const int eleGID
)
{
   // get gauss point number
    const int gp = params.get<int>("gp", -1);
    if (gp == -1)
      dserror("no Gauss point number provided in material");

  // calculate growth part F_g of the deformation gradient F
  LINALG::Matrix<3,3> F_g(true);
  CalcFg(theta,gp,F_g);

  // calculate F_g^(-1)
  LINALG::Matrix<3,3> F_ginv(true);
  F_ginv.Invert(F_g);

  //elastic deformation gradient F_e = F * F_g^(-1)
  LINALG::Matrix<3, 3> defgrddach(true);//*defgrd);
  defgrddach.MultiplyNN(*defgrd,F_ginv); // Scale(1.0 / theta);

  // elastic right Cauchy-Green Tensor Cdach = F_e^T * F_e (= F_g^-T C F_g^-1)
  LINALG::Matrix<3, 3> Cdach(true);
  Cdach.MultiplyTN(defgrddach,defgrddach);

  //transform Cdach into a vector
  LINALG::Matrix<6, 1> Cdachvec(true);
  MatrixToVector(Cdach,Cdachvec,MAT::voigt_strain);

  //--------------------------------------------------------------------------------------
  // call material law with elastic part of defgr and elastic part of glstrain
  //--------------------------------------------------------------------------------------
  // build identity tensor I
  LINALG::Matrix<6, 1> Id(true);
  for (int i = 0; i < 3; i++)
    Id(i) = 1.0;

  LINALG::Matrix<6, 1> glstraindachvec(Cdachvec);
  glstraindachvec -= Id;
  glstraindachvec.Scale(0.5);

  LINALG::Matrix<6, 1> Sdachvec(true);
  // elastic 2 PK stress and constitutive matrix
  Matelastic()->Evaluate(&defgrddach,
                        &glstraindachvec,
                        params,
                        &Sdachvec,
                        cmatdach,
                        eleGID);

  // calculate stress
  // 2PK stress S = F_g^-1 Sdach F_g^-T
  LINALG::Matrix<3, 3> Sdach(true);
  VectorToMatrix(Sdach,Sdachvec,MAT::voigt_stress);

  LINALG::Matrix<3, 3> tmp(true);
  tmp.MultiplyNT(Sdach,F_ginv);
  LINALG::Matrix<3, 3> S(true);
  S.MultiplyNN(F_ginv,tmp);

  MatrixToVector(S,*stress,MAT::voigt_stress);
}


/*----------------------------------------------------------------------*
 |     ///transform vector in voigt notation into symmetric             |
 | 2Tensor (in matrix notation)                              Thon 01/15 |
 *----------------------------------------------------------------------*/
void MAT::GrowthScdACRadial::VectorToMatrix(LINALG::Matrix<3,3>& Matrix,
                    const LINALG::Matrix<6,1>& Vector,
                    const MAT::VoigtType Type )
{
  double alpha;

  switch (Type)
  {
  case MAT::voigt_stress:
    alpha=1.0;
    break;
  case MAT::voigt_strain:
    alpha=0.5;
    break;
  default:
    dserror("No supported VoigtType!");
    alpha=1.0;
    break;
  }

  Matrix(0,0)=Vector(0);
  Matrix(0,1)=alpha*Vector(3);
  Matrix(0,2)=alpha*Vector(5);
  Matrix(1,0)=Matrix(0,1);//alpha*Vector(3);
  Matrix(1,1)=Vector(1);
  Matrix(1,2)=alpha*Vector(4);
  Matrix(2,0)=Matrix(0,2);//alpha*Vector(5);
  Matrix(2,1)=Matrix(1,2);//alpha*Vector(4);
  Matrix(2,2)=Vector(2);

}

/*----------------------------------------------------------------------*
 | ///transform symmetric 2Tensor(in matrix notation) into voigt        |
 | notation ( e.g. vector notation)                          Thon 01/15 |
 *----------------------------------------------------------------------*/
void MAT::GrowthScdACRadial::MatrixToVector(const LINALG::Matrix<3,3>& Matrix,
                    LINALG::Matrix<6,1>& Vector,
                    const MAT::VoigtType Type)
{
  double alpha;
  switch (Type)
  {
  case MAT::voigt_stress:
    alpha=1.0;
    break;
  case MAT::voigt_strain:
    alpha=2.0;
    break;
  default:
    dserror("No supported VoigtType!");
    alpha=1.0;
    break;
  }
  Vector(0)=Matrix(0,0);
  Vector(1)=Matrix(1,1);
  Vector(2)=Matrix(2,2);
  Vector(3)=alpha*Matrix(0,1);
  Vector(4)=alpha*Matrix(1,2);
  Vector(5)=alpha*Matrix(0,2);
}

/*----------------------------------------------------------------------*
 | pull back of a symmertic elastic 4th order tensor (in matrix/voigt   |
 | notation) via the 2th order deformation gradient (also in matrix     |
 | notation)                                                  thon 01/15|
 *----------------------------------------------------------------------*/
LINALG::Matrix<6,6> MAT::GrowthScdACRadial::PullBack4Tensor(const LINALG::Matrix<3,3>& defgr,
                    const LINALG::Matrix<6,6>& Cmat)
{
  double CMAT[3][3][3][3] = {{{{0.0}}}};
  Setup4Tensor(CMAT,Cmat);
//  PrintFourTensor(CMAT);

//This would be the long way....
//
//  double tmp[3][3][3][3] = {{{{0.0}}}};
//  for(int i=0;i<3;++i)
//    for(int j=0;j<3;++j)
//      for(int k=0;k<3;++k)
//        for(int l=0;l<3;++l)
//          for(int A=0;A<3;++A)
//            for(int B=0;B<3;++B)
//              for(int C=0;C<3;++C)
//                for(int D=0;D<3;++D)
//                  tmp[i][j][k][l] += defgr(i,A)*defgr(j,B)*defgr(k,C)*defgr(l,D)*CMAT[A][B][C][D];
//        }
//  PrintFourTensor(tmp);
//  LINALG::Matrix<6,6> CResult1(true);
//  Setup6x6VoigtMatrix(CResult1,tmp1);

  //But we can use the fact that CResult(i,j,k,l)=CResult(k,l,i,j) iff we have a hyperelatic material
  LINALG::Matrix<6,6> CResult(true);

  CResult(0,0)=PullBack4Tensorijkl(defgr,CMAT,0,0,0,0);
  CResult(0,1)=PullBack4Tensorijkl(defgr,CMAT,0,0,1,1);
  CResult(0,2)=PullBack4Tensorijkl(defgr,CMAT,0,0,2,2);
  CResult(0,3)=PullBack4Tensorijkl(defgr,CMAT,0,0,0,1);
  CResult(0,4)=PullBack4Tensorijkl(defgr,CMAT,0,0,1,2);
  CResult(0,5)=PullBack4Tensorijkl(defgr,CMAT,0,0,0,2);
  CResult(1,0)=CResult(0,1);
  CResult(1,1)=PullBack4Tensorijkl(defgr,CMAT,1,1,1,1);
  CResult(1,2)=PullBack4Tensorijkl(defgr,CMAT,1,1,2,2);
  CResult(1,3)=PullBack4Tensorijkl(defgr,CMAT,1,1,0,1);
  CResult(1,4)=PullBack4Tensorijkl(defgr,CMAT,1,1,1,2);
  CResult(1,5)=PullBack4Tensorijkl(defgr,CMAT,1,1,0,2);
  CResult(2,0)=CResult(0,2);
  CResult(2,1)=CResult(1,2);
  CResult(2,2)=PullBack4Tensorijkl(defgr,CMAT,2,2,2,2);
  CResult(2,3)=PullBack4Tensorijkl(defgr,CMAT,2,2,0,1);
  CResult(2,4)=PullBack4Tensorijkl(defgr,CMAT,2,2,1,2);
  CResult(2,5)=PullBack4Tensorijkl(defgr,CMAT,2,2,0,2);
  CResult(3,0)=CResult(0,3);
  CResult(3,1)=CResult(1,3);
  CResult(3,2)=CResult(2,3);
  CResult(3,3)=PullBack4Tensorijkl(defgr,CMAT,0,1,0,1);
  CResult(3,4)=PullBack4Tensorijkl(defgr,CMAT,0,1,1,2);
  CResult(3,5)=PullBack4Tensorijkl(defgr,CMAT,0,1,0,2);
  CResult(4,0)=CResult(0,4);
  CResult(4,1)=CResult(1,4);
  CResult(4,2)=CResult(2,4);
  CResult(4,3)=CResult(3,4);
  CResult(4,4)=PullBack4Tensorijkl(defgr,CMAT,1,2,1,2);
  CResult(4,5)=PullBack4Tensorijkl(defgr,CMAT,1,2,0,2);
  CResult(5,0)=CResult(0,5);
  CResult(5,1)=CResult(1,5);
  CResult(5,2)=CResult(2,5);
  CResult(5,3)=CResult(3,5);
  CResult(5,4)=CResult(4,5);
  CResult(5,5)=PullBack4Tensorijkl(defgr,CMAT,0,2,0,2);

  return CResult;
}

/*-------------------------------------------------------------------------------------*
 |  pull back the ijkl-th entry of a symmetric elastic 4th order                       |
 | tensor (in matrix/voigt notation) /// via the 2th order deformation                 |
 | gradient (also in matrix notation)                                      Thon  01/15 |
 *-------------------------------------------------------------------------------------*/
double MAT::GrowthScdACRadial::PullBack4Tensorijkl(const LINALG::Matrix<3,3>& defgr,
                      const double (&FourTensor)[3][3][3][3],
                      const double& i,
                      const double& j,
                      const double& k,
                      const double& l)
{
  double CResult_ijkl=0;

  for(int A=0;A<3;++A)
    for(int B=0;B<3;++B)
      for(int C=0;C<3;++C)
        for(int D=0;D<3;++D)
          CResult_ijkl += defgr(i,A)*defgr(j,B)*defgr(k,C)*defgr(l,D)*FourTensor[A][B][C][D];

  return CResult_ijkl;
}

/*-------------------------------------------------------------------------------------*
 |  Setup 4-Tensor from 6x6 Voigt notation                                 thon  01/15 |
 *-------------------------------------------------------------------------------------*/
void MAT::GrowthScdACRadial::Setup4Tensor(
    double (&FourTensor)[3][3][3][3],
    const LINALG::Matrix<6,6>& VoigtMatrix
)
{
  //Clear4Tensor(FourTensor);
  // Setup 4-Tensor from 6x6 Voigt matrix (which has to be the representative of a 4 tensor with at least minor symmetries)
  FourTensor[0][0][0][0] = VoigtMatrix(0,0);//C1111
  FourTensor[0][0][1][1] = VoigtMatrix(0,1);//C1122
  FourTensor[0][0][2][2] = VoigtMatrix(0,2);//C1133
  FourTensor[0][0][0][1] = VoigtMatrix(0,3);//C1112
  FourTensor[0][0][1][0] = VoigtMatrix(0,3);//C1121
  FourTensor[0][0][1][2] = VoigtMatrix(0,4);//C1123
  FourTensor[0][0][2][1] = VoigtMatrix(0,4);//C1132
  FourTensor[0][0][0][2] = VoigtMatrix(0,5);//C1113
  FourTensor[0][0][2][0] = VoigtMatrix(0,5);//C1131

  FourTensor[1][1][0][0] = VoigtMatrix(1,0);//C2211
  FourTensor[1][1][1][1] = VoigtMatrix(1,1);//C2222
  FourTensor[1][1][2][2] = VoigtMatrix(1,2);//C2233
  FourTensor[1][1][0][1] = VoigtMatrix(1,3);//C2212
  FourTensor[1][1][1][0] = VoigtMatrix(1,3);//C2221
  FourTensor[1][1][1][2] = VoigtMatrix(1,4);//C2223
  FourTensor[1][1][2][1] = VoigtMatrix(1,4);//C2232
  FourTensor[1][1][0][2] = VoigtMatrix(1,5);//C2213
  FourTensor[1][1][2][0] = VoigtMatrix(1,5);//C2231

  FourTensor[2][2][0][0] = VoigtMatrix(2,0);//C3311
  FourTensor[2][2][1][1] = VoigtMatrix(2,1);//C3322
  FourTensor[2][2][2][2] = VoigtMatrix(2,2);//C3333
  FourTensor[2][2][0][1] = VoigtMatrix(2,3);//C3312
  FourTensor[2][2][1][0] = VoigtMatrix(2,3);//C3321
  FourTensor[2][2][1][2] = VoigtMatrix(2,4);//C3323
  FourTensor[2][2][2][1] = VoigtMatrix(2,4);//C3332
  FourTensor[2][2][0][2] = VoigtMatrix(2,5);//C3313
  FourTensor[2][2][2][0] = VoigtMatrix(2,5);//C3331

  FourTensor[0][1][0][0] = VoigtMatrix(3,0);
  FourTensor[1][0][0][0] = VoigtMatrix(3,0);//C1211 = C2111
  FourTensor[0][1][1][1] = VoigtMatrix(3,1);
  FourTensor[1][0][1][1] = VoigtMatrix(3,1);//C1222 = C2122
  FourTensor[0][1][2][2] = VoigtMatrix(3,2);
  FourTensor[1][0][2][2] = VoigtMatrix(3,2);//C1233 = C2133
  FourTensor[0][1][0][1] = VoigtMatrix(3,3);
  FourTensor[1][0][0][1] = VoigtMatrix(3,3);//C1212 = C2112
  FourTensor[0][1][1][0] = VoigtMatrix(3,3);
  FourTensor[1][0][1][0] = VoigtMatrix(3,3);//C1221 = C2121
  FourTensor[0][1][1][2] = VoigtMatrix(3,4);
  FourTensor[1][0][1][2] = VoigtMatrix(3,4);//C1223 = C2123
  FourTensor[0][1][2][1] = VoigtMatrix(3,4);
  FourTensor[1][0][2][1] = VoigtMatrix(3,4);//C1232 = C2132
  FourTensor[0][1][0][2] = VoigtMatrix(3,5);
  FourTensor[1][0][0][2] = VoigtMatrix(3,5);//C1213 = C2113
  FourTensor[0][1][2][0] = VoigtMatrix(3,5);
  FourTensor[1][0][2][0] = VoigtMatrix(3,5);//C1231 = C2131

  FourTensor[1][2][0][0] = VoigtMatrix(4,0);
  FourTensor[2][1][0][0] = VoigtMatrix(4,0);//C2311 = C3211
  FourTensor[1][2][1][1] = VoigtMatrix(4,1);
  FourTensor[2][1][1][1] = VoigtMatrix(4,1);//C2322 = C3222
  FourTensor[1][2][2][2] = VoigtMatrix(4,2);
  FourTensor[2][1][2][2] = VoigtMatrix(4,2);//C2333 = C3233
  FourTensor[1][2][0][1] = VoigtMatrix(4,3);
  FourTensor[2][1][0][1] = VoigtMatrix(4,3);//C2312 = C3212
  FourTensor[1][2][1][0] = VoigtMatrix(4,3);
  FourTensor[2][1][1][0] = VoigtMatrix(4,3);//C2321 = C3221
  FourTensor[1][2][1][2] = VoigtMatrix(4,4);
  FourTensor[2][1][1][2] = VoigtMatrix(4,4);//C2323 = C3223
  FourTensor[1][2][2][1] = VoigtMatrix(4,4);
  FourTensor[2][1][2][1] = VoigtMatrix(4,4);//C2332 = C3232
  FourTensor[1][2][0][2] = VoigtMatrix(4,5);
  FourTensor[2][1][0][2] = VoigtMatrix(4,5);//C2313 = C3213
  FourTensor[1][2][2][0] = VoigtMatrix(4,5);
  FourTensor[2][1][2][0] = VoigtMatrix(4,5);//C2331 = C3231

  FourTensor[0][2][0][0] = VoigtMatrix(5,0);
  FourTensor[2][0][0][0] = VoigtMatrix(5,0);//C1311 = C3111
  FourTensor[0][2][1][1] = VoigtMatrix(5,1);
  FourTensor[2][0][1][1] = VoigtMatrix(5,1);//C1322 = C3122
  FourTensor[0][2][2][2] = VoigtMatrix(5,2);
  FourTensor[2][0][2][2] = VoigtMatrix(5,2);//C1333 = C3133
  FourTensor[0][2][0][1] = VoigtMatrix(5,3);
  FourTensor[2][0][0][1] = VoigtMatrix(5,3);//C1312 = C3112
  FourTensor[0][2][1][0] = VoigtMatrix(5,3);
  FourTensor[2][0][1][0] = VoigtMatrix(5,3);//C1321 = C3121
  FourTensor[0][2][1][2] = VoigtMatrix(5,4);
  FourTensor[2][0][1][2] = VoigtMatrix(5,4);//C1323 = C3123
  FourTensor[0][2][2][1] = VoigtMatrix(5,4);
  FourTensor[2][0][2][1] = VoigtMatrix(5,4);//C1332 = C3132
  FourTensor[0][2][0][2] = VoigtMatrix(5,5);
  FourTensor[2][0][0][2] = VoigtMatrix(5,5);//C1313 = C3113
  FourTensor[0][2][2][0] = VoigtMatrix(5,5);
  FourTensor[2][0][2][0] = VoigtMatrix(5,5);//C1331 = C3131

}  // Setup4Tensor()


///*------------------------------------------------------------------------------------------*
// |  Setup 6x6 matrix in Voigt notation from 4-Tensor                            thon  01/15 |
// *------------------------------------------------------------------------------------------*/
//void MAT::GrowthScdACRadial::Setup6x6VoigtMatrix(
//    LINALG::Matrix<6,6>& VoigtMatrix,
//    const double (&FourTensor)[3][3][3][3]
//)
//{
/////*  [      C1111                 C1122                C1133                0.5*(C1112+C1121)               0.5*(C1123+C1132)                0.5*(C1113+C1131)      ]
////    [      C2211                 C2222                C2233                0.5*(C2212+C2221)               0.5*(C2223+C2232)                0.5*(C2213+C2231)      ]
////    [      C3311                 C3322                C3333                0.5*(C3312+C3321)               0.5*(C3323+C3332)                0.5*(C3313+C3331)      ]
////    [0.5*(C1211+C2111)    0.5*(C1222+C2122)    0.5*(C1233+C2133)    0.5*(C1212+C2112+C1221+C2121)    0.5*(C1223+C2123+C1232+C2132)    0.5*(C1213+C2113+C1231+C2131)]
////    [0.5*(C2311+C3211)    0.5*(C2322+C3222)    0.5*(C2333+C3233)    0.5*(C2312+C3212+C2321+C3221)    0.5*(C2323+C3223+C2332+C3232)    0.5*(C2313+C3213+C2331+C3231)]
////    [0.5*(C1322+C3122)    0.5*(C1322+C3122)    0.5*(C1333+C3133)    0.5*(C1312+C3112+C1321+C3121)    0.5*(C1323+C3123+C1332+C3132)    0.5*(C1313+C3113+C1331+C3131)] */
//
//  // Setup 4-Tensor from 6x6 Voigt matrix
//  VoigtMatrix(0,0) = FourTensor[0][0][0][0]; //C1111
//  VoigtMatrix(0,1) = FourTensor[0][0][1][1]; //C1122
//  VoigtMatrix(0,2) = FourTensor[0][0][2][2]; //C1133
//  VoigtMatrix(0,3) = 0.5 * (FourTensor[0][0][0][1] + FourTensor[0][0][1][0]); //0.5*(C1112+C1121)
//  VoigtMatrix(0,4) = 0.5 * (FourTensor[0][0][1][2] + FourTensor[0][0][2][1]); //0.5*(C1123+C1132)
//  VoigtMatrix(0,5) = 0.5 * (FourTensor[0][0][0][2] + FourTensor[0][0][2][0]); //0.5*(C1113+C1131)
//
//  VoigtMatrix(1,0) = FourTensor[1][1][0][0]; //C2211
//  VoigtMatrix(1,1) = FourTensor[1][1][1][1]; //C2222
//  VoigtMatrix(1,2) = FourTensor[1][1][2][2]; //C2233
//  VoigtMatrix(1,3) = 0.5 * (FourTensor[1][1][0][1] + FourTensor[1][1][1][0]); //0.5*(C2212+C2221)
//  VoigtMatrix(1,4) = 0.5 * (FourTensor[1][1][1][2] + FourTensor[1][1][2][1]); //0.5*(C2223+C2232)
//  VoigtMatrix(1,5) = 0.5 * (FourTensor[1][1][0][2] + FourTensor[1][1][2][0]); //0.5*(C2213+C2231)
//
//  VoigtMatrix(2,0) = FourTensor[2][2][0][0]; //C3311
//  VoigtMatrix(2,1) = FourTensor[2][2][1][1]; //C3322
//  VoigtMatrix(2,2) = FourTensor[2][2][2][2]; //C3333
//  VoigtMatrix(2,3) = 0.5 * (FourTensor[2][2][0][1] + FourTensor[2][2][1][0]); //0.5*(C3312+C3321)
//  VoigtMatrix(2,4) = 0.5 * (FourTensor[2][2][1][2] + FourTensor[2][2][2][1]); //0.5*(C3323+C3332)
//  VoigtMatrix(2,5) = 0.5 * (FourTensor[2][2][0][2] + FourTensor[2][2][2][0]); //0.5*(C3313+C3331)
//
//  VoigtMatrix(3,0) = 0.5 * (FourTensor[0][1][0][0] + FourTensor[1][0][0][0]); //0.5*(C1211+C2111)
//  VoigtMatrix(3,1) = 0.5 * (FourTensor[0][1][1][1] + FourTensor[1][0][1][1]); //0.5*(C1222+C2122)
//  VoigtMatrix(3,2) = 0.5 * (FourTensor[0][1][2][2] + FourTensor[1][0][2][2]); //0.5*(C1233+C2133)
//  VoigtMatrix(3,3) = 0.25 * (FourTensor[0][1][0][1] + FourTensor[1][0][0][1] + FourTensor[0][1][1][0] + FourTensor[1][0][1][0]); //0.5*(C1212+C2112+C1221+C2121)
//  VoigtMatrix(3,4) = 0.25 * (FourTensor[0][1][1][2] + FourTensor[1][0][1][2] + FourTensor[0][1][2][1] + FourTensor[1][0][2][1]); //0.5*(C1223+C2123+C1232+C2132)
//  VoigtMatrix(3,5) = 0.25 * (FourTensor[0][1][0][2] + FourTensor[1][0][0][2] + FourTensor[0][1][2][0] + FourTensor[1][0][2][0]); //0.5*(C1213+C2113+C1231+C2131)
//
//  VoigtMatrix(4,0) = 0.5 * (FourTensor[1][2][0][0] + FourTensor[2][1][0][0]); //0.5*(C2311+C3211)
//  VoigtMatrix(4,1) = 0.5 * (FourTensor[1][2][1][1] + FourTensor[2][1][1][1]); //0.5*(C2322+C3222)
//  VoigtMatrix(4,2) = 0.5 * (FourTensor[1][2][2][2] + FourTensor[2][1][2][2]); //0.5*(C2333+C3233)
//  VoigtMatrix(4,3) = 0.25 * (FourTensor[1][2][0][1] + FourTensor[2][1][0][1] + FourTensor[1][2][1][0] + FourTensor[2][1][1][0]); //0.5*(C2312+C3212+C2321+C3221)
//  VoigtMatrix(4,4) = 0.25 * (FourTensor[1][2][1][2] + FourTensor[2][1][1][2] + FourTensor[1][2][2][1] + FourTensor[2][1][2][1]); //0.5*(C2323+C3223+C2332+C3232)
//  VoigtMatrix(4,5) = 0.25 * (FourTensor[1][2][0][2] + FourTensor[2][1][0][2] + FourTensor[1][2][2][0] + FourTensor[2][1][2][0]); //0.5*(C2313+C3213+C2331+C3231)
//
//  VoigtMatrix(5,0) = 0.5 * (FourTensor[0][2][0][0] + FourTensor[2][0][0][0]); //0.5*(C1311+C3111)
//  VoigtMatrix(5,1) = 0.5 * (FourTensor[0][2][1][1] + FourTensor[2][0][1][1]); //0.5*(C1322+C3122)
//  VoigtMatrix(5,2) = 0.5 * (FourTensor[0][2][2][2] + FourTensor[2][0][2][2]); //0.5*(C1333+C3133)
//  VoigtMatrix(5,3) = 0.25 * (FourTensor[0][2][0][1] + FourTensor[2][0][0][1] + FourTensor[0][2][1][0] + FourTensor[2][0][1][0]); //0.5*(C1312+C3112+C1321+C3121)
//  VoigtMatrix(5,4) = 0.25 * (FourTensor[0][2][1][2] + FourTensor[2][0][1][2] + FourTensor[0][2][2][1] + FourTensor[2][0][2][1]); //0.5*(C1323+C3123+C1332+C3132)
//  VoigtMatrix(5,5) = 0.25 * (FourTensor[0][2][0][2] + FourTensor[2][0][0][2] + FourTensor[0][2][2][0] + FourTensor[2][0][2][0]); //0.5*(C1313+C3113+C1331+C3131)
//
//}  // Setup6x6VoigtMatrix()

/*------------------------------------------------------------------------------------------*
 |  Print Four Tensor                                                           thon  01/15 |
 *------------------------------------------------------------------------------------------*/
void MAT::GrowthScdACRadial::PrintFourTensor(
    double (&FourTensor)[3][3][3][3]
                                            )
{
  std::cout<<"-----------------Print Four Tensor--------------"<<std::endl;

  for(int i=0;i<3;++i)
    for(int j=0;j<3;++j)
      for(int k=0;k<3;++k)
        for(int l=0;l<3;++l)
          std::cout<<"ELEMENT "<<i<<j<<k<<l<<" : "<<FourTensor[i][j][k][l]<<std::endl;

  std::cout<<"------------------------------------------------"<<std::endl;
  return;
}

/*----------------------------------------------------------------------------*/
void MAT::GrowthScdACRadial::VisNames(std::map<std::string,int>& names)
{
  std::string fiber = "growth_direction";
  names[fiber] = 3; // 3-dim vector

  GrowthScdAC::VisNames(names);
}

/*----------------------------------------------------------------------------*/
bool MAT::GrowthScdACRadial::VisData(const std::string& name, std::vector<double>& data, int numgp , int eleID)
{
  if (name == "growth_direction")
  {
    if ((int)data.size()!=3)
      dserror("size mismatch");

    LINALG::Matrix<3,1> temp(true);
    for (int gp=0; gp<numgp; gp++)
    {
      if ( (abs((n_. at(gp)).Norm2())-1.0)>1e-14 )
      {
        dserror("length of n_ is not one. This should not be possible!");
      }

      temp.Update(1.0/(double)numgp,n_.at(gp),1.0);
    }
    // we fix the length, since the mean of all n's must not have length one, even if
    // all n's individually have length one! So:
    temp.Scale(1.0/temp.Norm2());

    data[0] = temp(0,0);
    data[1] = temp(1,0);
    data[2] = temp(2,0);
  }

  GrowthScdAC::VisData(name,data,numgp,eleID);

  return true;
}
