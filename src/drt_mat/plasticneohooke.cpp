/*----------------------------------------------------------------------*/
/*!
\file plasticneohooke.cpp 
\brief Contains the functions to establish local material law /
       stress-strain law for isotropic material for a 3D hex element
       following rate-independent von Mises plasticity with combined
       (nonlinear) isotropic and kinematic hardening and a nearly
       incompressible Neo Hookean material law.
       
       Refer also to the master thesis of A.A. van der Stelt (2009).

<pre>
Maintainer: Alexander Popp
            popp@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15264
</pre>
*/
/*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include <vector>
#include <Epetra_SerialDenseMatrix.h>
#include <Epetra_SerialDenseVector.h>
#include "plasticneohooke.H"
#include "../drt_lib/linalg_utils.H"

//#include <iostream>

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::PAR::PlasticNeoHooke::PlasticNeoHooke(
  Teuchos::RCP<MAT::PAR::Material> matdata
  )
: Parameter(matdata),
  youngs_(matdata->GetDouble("YOUNG")),
  poissonratio_(matdata->GetDouble("NUE")),
  density_(matdata->GetDouble("DENS")),
  isohard_(matdata->GetDouble("ISOHARD")),
  yield_(matdata->GetDouble("YIELD")),
  infyield_(matdata->GetDouble("INFYIELD")),
  exp_(matdata->GetDouble("EXP")),
  kinhard_(matdata->GetDouble("KINHARD"))
{
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::PlasticNeoHooke::PlasticNeoHooke()
  : params_(NULL)
{
  isinit_=false;
  histplasticrcgcurr_=rcp(new vector<LINALG::Matrix<NUM_STRESS_3D,1> >);
  histeplasticscurr_=rcp(new vector<LINALG::Matrix<1,1> >); 
  histplasticrcglast_=rcp(new vector<LINALG::Matrix<NUM_STRESS_3D,1> >);
  histeplasticslast_=rcp(new vector<LINALG::Matrix<1,1> >);
  histplasticbstresscurr_=rcp(new vector<LINALG::Matrix<3,1> >); /*added*/
  histplasticbstresslast_=rcp(new vector<LINALG::Matrix<3,1> >); /*added*/
  lamdacurr_=rcp(new vector<LINALG::Matrix<3,1> >); /*added*/
  lamdalast_=rcp(new vector<LINALG::Matrix<3,1> >); /*added*/
  nsetcurr_=rcp(new vector<LINALG::Matrix<3,3> >); /*added*/
  nsetlast_=rcp(new vector<LINALG::Matrix<3,3> >); /*added*/
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::PlasticNeoHooke::PlasticNeoHooke(MAT::PAR::PlasticNeoHooke* params)
  : params_(params)
{
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::PlasticNeoHooke::Pack(vector<char>& data) const
{  
  data.resize(0);

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);
  
  // matid
  int matid = -1;
  if (params_ != NULL) matid = params_->Id();  // in case we are in post-process mode
  AddtoPack(data,matid);

  //Pack history data
  int histsize;
  if (!Initialized())
  {
    histsize=0;
  }
  else
  {
    histsize = histplasticrcglast_->size();
  }

  AddtoPack(data,2*histsize); // Length of history vector(s)
  for (int var = 0; var < histsize; ++var)
  {
    AddtoPack(data,histplasticrcglast_->at(var));
    AddtoPack(data,histeplasticslast_->at(var));
    AddtoPack(data,histplasticbstresslast_->at(var)); /*added*/
    AddtoPack(data,lamdalast_->at(var)); /*added*/
    AddtoPack(data,nsetlast_->at(var)); /*added*/
  } 
  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::PlasticNeoHooke::Unpack(const vector<char>& data)
{
  isinit_=true;
  int position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position,data,type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");

  // matid and recover params_
  int matid;
  ExtractfromPack(position,data,matid);
  // in post-process mode we do not have any instance of DRT::Problem
  if (DRT::Problem::NumInstances() > 0)
  {
    const int probinst = DRT::Problem::Instance()->Materials()->GetReadFromProblem();
  MAT::PAR::Parameter* mat = DRT::Problem::Instance(probinst)->Materials()->ParameterById(matid);
  if (mat->Type() == MaterialType())
    params_ = static_cast<MAT::PAR::PlasticNeoHooke*>(mat);
  else
    dserror("Type of parameter material %d does not fit to calling type %d", mat->Type(), MaterialType());
  }
  else
  {
    params_ = NULL;
  }

  // history data
  int twicehistsize;
  ExtractfromPack(position,data,twicehistsize);

  if (twicehistsize == 0) isinit_=false;

  histplasticrcgcurr_=rcp(new vector<LINALG::Matrix<NUM_STRESS_3D,1> >);
  histeplasticscurr_=rcp(new vector<LINALG::Matrix<1,1> >); 
  histplasticrcglast_=rcp(new vector<LINALG::Matrix<NUM_STRESS_3D,1> >);
  histeplasticslast_=rcp(new vector<LINALG::Matrix<1,1> >);
  histplasticbstresscurr_=rcp(new vector<LINALG::Matrix<3,1> >); /*added*/
  histplasticbstresslast_=rcp(new vector<LINALG::Matrix<3,1> >); /*added*/
  lamdacurr_=rcp(new vector<LINALG::Matrix<3,1> >); /*added*/
  lamdalast_=rcp(new vector<LINALG::Matrix<3,1> >); /*added*/
  nsetcurr_=rcp(new vector<LINALG::Matrix<3,3> >); /*added*/
  nsetlast_=rcp(new vector<LINALG::Matrix<3,3> >); /*added*/

  for (int var=0; var<twicehistsize; var+=2) 
  {
    LINALG::Matrix<NUM_STRESS_3D,1> tmp1(true);
    LINALG::Matrix<1,1> tmp2(true);
    LINALG::Matrix<3,1> tmp3(true); /*added*/
    LINALG::Matrix<3,3> tmp4(true); /*added*/
    histplasticrcgcurr_->push_back(tmp1);
    histeplasticscurr_->push_back(tmp2);
    histplasticbstresscurr_->push_back(tmp3); /*added*/
    lamdacurr_->push_back(tmp3); /*added*/
    nsetcurr_->push_back(tmp4); /*added*/
    ExtractfromPack(position,data,tmp1);
    histplasticrcglast_->push_back(tmp1);
    ExtractfromPack(position,data,tmp2);
    histeplasticslast_->push_back(tmp2);
    ExtractfromPack(position,data,tmp3); /*added*/
    histplasticbstresslast_->push_back(tmp3); /*added*/
    ExtractfromPack(position,data,tmp3); /*added*/
    lamdalast_->push_back(tmp3); /*added*/
    ExtractfromPack(position,data,tmp4); /*added*/
    nsetlast_->push_back(tmp4); /*added*/
  }

  if (position != (int)data.size())
    dserror("Mismatch in size of data %d <-> %d",(int)data.size(),position);

  return;
}

/*---------------------------------------------------------------------*
|  Initialise/allocate internal stress variables (public)         11/09|
*----------------------------------------------------------------------*/
void MAT::PlasticNeoHooke::Setup(const int numgp)
{
  histplasticrcgcurr_=rcp(new vector<LINALG::Matrix<NUM_STRESS_3D,1> >);
  histeplasticscurr_=rcp(new vector<LINALG::Matrix<1,1> >);
  histplasticbstresscurr_=rcp(new vector<LINALG::Matrix<3,1> >); /*added*/
  lamdacurr_=rcp(new vector<LINALG::Matrix<3,1> >); /*added*/
  nsetcurr_=rcp(new vector<LINALG::Matrix<3,3> >); /*added*/
  histplasticrcglast_=rcp(new vector<LINALG::Matrix<NUM_STRESS_3D,1> >);
  histeplasticslast_=rcp(new vector<LINALG::Matrix<1,1> >);
  histplasticbstresslast_=rcp(new vector<LINALG::Matrix<3,1> >); /*added*/
  lamdalast_=rcp(new vector<LINALG::Matrix<3,1> >); /*added*/
  nsetlast_=rcp(new vector<LINALG::Matrix<3,3> >); /*added*/
  LINALG::Matrix<NUM_STRESS_3D,1> emptyvec1(true);
  for (int j=0;j<3;j++) {emptyvec1(j,0)=1.0;} //create Identity matrix
  const LINALG::Matrix<1,1> emptyvec2(true);
  const LINALG::Matrix<3,1> emptyvec3(true); /*added*/
  LINALG::Matrix<3,1> emptyvec4(true); /*added*/
  for (int j=0;j<3;j++) {emptyvec4(j,0)=1.0;} /*added*/
  const LINALG::Matrix<3,3> emptyvec5(true); /*added*/

  histplasticrcgcurr_->resize(numgp);
  histeplasticscurr_->resize(numgp);
  histplasticbstresscurr_->resize(numgp); /*added*/
  lamdacurr_->resize(numgp); /*added*/
  nsetcurr_->resize(numgp); /*added*/
  histplasticrcglast_->resize(numgp);
  histeplasticslast_->resize(numgp);
  histplasticbstresslast_->resize(numgp); /*added*/
  lamdalast_->resize(numgp); /*added*/
  nsetlast_->resize(numgp); /*added*/


  for (int j=0; j<numgp;++j)  {
    histplasticrcgcurr_->at(j) = emptyvec1;
    histeplasticscurr_->at(j) = emptyvec2;
    histplasticrcglast_->at(j) = emptyvec1;
    histeplasticslast_->at(j) = emptyvec2;
    histplasticbstresscurr_->at(j) = emptyvec3; /*added*/ 
    histplasticbstresslast_->at(j) = emptyvec3; /*added*/
    lamdacurr_->at(j) = emptyvec4; /*added*/ 
    lamdalast_->at(j) = emptyvec4; /*added*/
    nsetcurr_->at(j) = emptyvec5; /*added*/ 
    nsetlast_->at(j) = emptyvec5; /*added*/
 
  }
  isinit_=true;
  return;
}

/*---------------------------------------------------------------------*
|  Update internal stress variables (public)                      11/09|
*----------------------------------------------------------------------*/
void MAT::PlasticNeoHooke::Update()
{
  histplasticrcglast_ = histplasticrcgcurr_;
  histeplasticslast_ = histeplasticscurr_;
  histplasticbstresslast_ = histplasticbstresscurr_; /*added*/
  lamdalast_ = lamdacurr_; /*added*/
  nsetlast_= nsetcurr_; /*added*/
  const LINALG::Matrix<NUM_STRESS_3D,1> emptyvec1(true);
  const LINALG::Matrix<1,1> emptyvec2(true);
  const LINALG::Matrix<3,1> emptyvec3(true); /*added*/
  const LINALG::Matrix<3,3> emptyvec4(true); /*added*/
  histplasticrcgcurr_=rcp(new vector<LINALG::Matrix<NUM_STRESS_3D,1> >);
  histeplasticscurr_=rcp(new vector<LINALG::Matrix<1,1> >); 
  histplasticbstresscurr_=rcp(new vector<LINALG::Matrix<3,1> >);/*added*/
  lamdacurr_=rcp(new vector<LINALG::Matrix<3,1> >);/*added*/
  nsetcurr_=rcp(new vector<LINALG::Matrix<3,3> >); /*added*/
  const int numgp=histplasticrcglast_->size();
  histplasticrcgcurr_->resize(numgp);
  histeplasticscurr_->resize(numgp);
  histplasticbstresscurr_->resize(numgp); /*added*/
  lamdacurr_->resize(numgp); /*added*/
  nsetcurr_->resize(numgp); /*added*/
  for (int j=0; j<numgp; ++j)  {
    histplasticrcgcurr_->at(j) = emptyvec1;
    histeplasticscurr_->at(j) = emptyvec2;
    histplasticbstresscurr_->at(j) = emptyvec3; /*added*/
    lamdacurr_->at(j) = emptyvec3; /*added*/
    nsetcurr_->at(j) = emptyvec4; /*added*/ 
  }

  return;
}

/*---------------------------------------------------------------------*
|  Reset internal stress variables (public)                       11/09|
*----------------------------------------------------------------------*/
void MAT::PlasticNeoHooke::Reset()
{
  // do nothing,
  // because #histplasticrcgcurr_ and #histeplasticscurr_ are recomputed
  // anyway at every iteration based upon #histplasticrcglast_ and
  // #histeplasticslast_ untouched within time step
  
  return;
}

/*----------------------------------------------------------------------*
 |  Evaluate material (public)                                     11/09|
 *----------------------------------------------------------------------*/
void MAT::PlasticNeoHooke::Evaluate
  (
    const LINALG::Matrix<3,3>* defgrd, ///<deformation gradient
    const int gp, ///< current Gauss point
    Teuchos::ParameterList& params,  ///parameter list for communication & HISTORY
    LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D>* cmat, ///<material stiffness matrix
    LINALG::Matrix<NUM_STRESS_3D,1>* stress ///< 2nd PK-stress
  )
{
  // get material parameters
  double ym =       params_->youngs_;       /// Young's modulus //
  double nu =       params_->poissonratio_; /// Poisson's ratio //
  double H =        params_->isohard_;      /// isotropic hardening parameter //
  double KH_y0 =    params_->yield_;        /// initial Kirchhoff yield stress //
  double sigmainf = params_->infyield_;     /// saturation model (nonlin. iso. hard.) 0.715 //
  double deltad =   params_->exp_;          /// saturation model (nonlin. iso. hard.) 16.93 //
  double Hkin =     params_->kinhard_;      /// kinematic hardening parameter //

  //cout << "\nisohardening " << H << " kinhardening " << Hkin << "\n";

  // initialize scalars
  double lambda = 0.0;
  double mue = 0.0;
  double kappa = 0.0;
  double f = 0.0;
  double dev_tau_abs = 0.0;
  double dj = 0.0;
  double cab = 0.0;
  double cbb = 0.0;
  double c1 = 0.0;

  // initialize vectors
  LINALG::Matrix<3,1> Le_trial;           // store trial eigenvalues
  LINALG::Matrix<3,1> Le;                 // store new eigenvalues
  LINALG::Matrix<3,1> dev_KH_trial;       // prinicple trial Kirchhoff deviatoric stresses
  LINALG::Matrix<3,1> dev_KH;             // principle Kirchhoff deviatoric stresses
  LINALG::Matrix<3,1> dev_PK2_principle;  // principle PK2 deviatoric stresses
  LINALG::Matrix<3,1> PK2_principle;      // prinple PK2 stresses
  LINALG::Matrix<3,3> PK2;                // PK2 stresses
  LINALG::Matrix<3,1> Vec;                // direction vector return mapping

  // initialize matrices
  LINALG::Matrix<3,3> be_trial;           // left-Cauchy Green tensor trial
  LINALG::Matrix<3,3> n_trial;
  LINALG::Matrix<3,3> be;
  LINALG::Matrix<3,3> N_trial; 

  // add value
  lambda = nu*ym/((1.0+nu)*(1.0-2.0*nu)); // Lame's first parameter lambda
  mue = ym/(2.0*(1.0+nu));                // Shear modulus parameter mu
  kappa = lambda+2.0/3.0*mue;             // Bulk modulus

  // build identity tensor I & C
  LINALG::Matrix<NUM_STRESS_3D,1> Id;
  for (int i = 0; i < 3; i++) Id(i) = 1.0;
  for (int i =3; i<6;i++) Id(i)=0.0;

  LINALG::Matrix<3,3> Chelp;
  Chelp.MultiplyTN((*defgrd),(*defgrd));
  
  LINALG::Matrix<NUM_STRESS_3D,1> C;
  
  C(0) = Chelp(0,0);
  C(1) = Chelp(1,1);
  C(2) = Chelp(2,2);
  C(3) = Chelp(0,1)+Chelp(1,0);
  C(4) = Chelp(1,2)+Chelp(2,1);
  C(5) = Chelp(0,2)+Chelp(2,0);
  
  // invariants
  // const double I1 = C(0)+C(1)+C(2); //1st invariant, trace
  double I3 = C(0)*C(1)*C(2)
     + 0.25 * C(3)*C(4)*C(5)
     - 0.25 * C(1)*C(5)*C(5)
     - 0.25 * C(2)*C(3)*C(3)
     - 0.25 * C(0)*C(4)*C(4); // 3rd invariant

  // invert C
  LINALG::Matrix<NUM_STRESS_3D,1> Cinv;
  Cinv(0) = C(1)*C(2) - 0.25*C(4)*C(4);
  Cinv(1) = C(0)*C(2) - 0.25*C(5)*C(5);
  Cinv(2) = C(0)*C(1) - 0.25*C(3)*C(3);
  Cinv(3) = 0.25*C(5)*C(4) - 0.5*C(3)*C(2);
  Cinv(4) = 0.25*C(3)*C(5) - 0.5*C(0)*C(4);
  Cinv(5) = 0.25*C(3)*C(4) - 0.5*C(5)*C(1);
  Cinv.Scale(1.0/I3);
  
  // Inverse plastic Right Cauchy-Green Tensor
  LINALG::Matrix<3,3> C_Plastic_Inv;
  C_Plastic_Inv(0,0) = (histplasticrcglast_->at(gp))(0,0);
  C_Plastic_Inv(1,1) = (histplasticrcglast_->at(gp))(1,0);
  C_Plastic_Inv(2,2) = (histplasticrcglast_->at(gp))(2,0);
  C_Plastic_Inv(0,1) = (histplasticrcglast_->at(gp))(3,0);
  C_Plastic_Inv(1,0) = (histplasticrcglast_->at(gp))(3,0);
  C_Plastic_Inv(1,2) = (histplasticrcglast_->at(gp))(4,0);
  C_Plastic_Inv(2,1) = (histplasticrcglast_->at(gp))(4,0);
  C_Plastic_Inv(0,2) = (histplasticrcglast_->at(gp))(5,0);
  C_Plastic_Inv(2,0) = (histplasticrcglast_->at(gp))(5,0);

  // Oldbackstress
  LINALG::Matrix<3,1> BSold; /*added*/
  BSold(0,0) = (histplasticbstresslast_->at(gp))(0,0); /*added*/
  BSold(1,0) = (histplasticbstresslast_->at(gp))(1,0); /*added*/
  BSold(2,0) = (histplasticbstresslast_->at(gp))(2,0); /*added*/

  LINALG::Matrix<3,1> Blamda; /*added*/
  Blamda(0,0) = (lamdalast_->at(gp))(0,0); /*added*/
  Blamda(1,0) = (lamdalast_->at(gp))(1,0); /*added*/
  Blamda(2,0) = (lamdalast_->at(gp))(2,0); /*added*/

  // Equivalent plastic strain
  double epsilon_p = (histeplasticslast_->at(gp))(0,0);
  
  // Energy elastic part: NeoHooke
  // NeoHooke with penalty W = W^dev(C)+U(J)
  // W = 1/2 mue (^I1-3) + 1/2 kappa (ln J)^2
  // ^I1 = tr(C^) 
  
  // Jacobian
  double J = (*defgrd).Determinant(); /// Jacobian
  
  // Pressure
  double p = kappa*(log(J)/J); /// Pressure
  
  // be_trial
  LINALG::Matrix<3,3> Matrix1;
  Matrix1.MultiplyNT(C_Plastic_Inv,*defgrd);
  be_trial.MultiplyNN(*defgrd,Matrix1); /// trial left Cauchy-Green tensor
  
  // Decompose be_trial
  // evaluate eigenproblem
  Epetra_SerialDenseVector Le_trial_E(3);
  Epetra_SerialDenseMatrix be_trial_E(3,3);
  for (int i = 0; i < 3; ++i)
  {
    for (int j = 0; j < 3; ++j)
    {
      be_trial_E(i,j) = be_trial(i,j);
    }
  }
  
  // watch out! be_trial_ matrix will temporarily hold eigenvectors!
  LINALG::SymmetricEigenProblem(be_trial_E,Le_trial_E); /// Eigenvectors and eigenvalues  
  for (int i = 0; i < 3; ++i)
  {
    Le_trial(i,0) = Le_trial_E(i,0);
    for (int j = 0; j < 3; ++j)
    {
      n_trial(i,j) = be_trial_E(i,j);
    }
  }

  /*
  (nsetcurr_->at(gp))(0,0) = n_trial(0,0);
  (nsetcurr_->at(gp))(0,1) = n_trial(0,1);
  (nsetcurr_->at(gp))(0,2) = n_trial(0,2);
  (nsetcurr_->at(gp))(1,0) = n_trial(1,0);
  (nsetcurr_->at(gp))(1,1) = n_trial(1,1);
  (nsetcurr_->at(gp))(1,2) = n_trial(1,2);
  (nsetcurr_->at(gp))(2,0) = n_trial(2,0);
  (nsetcurr_->at(gp))(2,1) = n_trial(2,1); 
  (nsetcurr_->at(gp))(2,2) = n_trial(2,2);
  */

  // Compute principal dev_stresses
  for (int i = 0; i<3; i++) dev_KH_trial(i,0) = mue*( log(Le_trial(i,0))-((2.0/3.0)*log(J)) );

  LINALG::Matrix<3,3> AB2;
  AB2(0,0) = (nsetlast_->at(gp))(0,0); AB2(0,1) = (nsetlast_->at(gp))(0,1); AB2(0,2) = (nsetlast_->at(gp))(0,2);
  AB2(1,0) = (nsetlast_->at(gp))(1,0); AB2(1,1) = (nsetlast_->at(gp))(1,1); AB2(1,2) = (nsetlast_->at(gp))(1,2);
  AB2(2,0) = (nsetlast_->at(gp))(2,0); AB2(2,1) = (nsetlast_->at(gp))(2,1); AB2(2,2) = (nsetlast_->at(gp))(2,2);
  
  //add//
  LINALG::Matrix<3,3> AB;
  LINALG::Matrix<3,1> ABhelp;
  LINALG::Matrix<3,1> VABtemp1;
  LINALG::Matrix<3,1> VABtemp2;

  for (int i = 0; i < 3; ++i) 
  {   
    for (int j = 0; j < 3; ++j) 
    { 
      VABtemp1(0,0) = n_trial(0,i); // check if eigenvectors are arranged like this in be_trial
      VABtemp1(1,0) = n_trial(1,i);
      VABtemp1(2,0) = n_trial(2,i);
  
      VABtemp2(0,0) = n_trial(0,j); // check if eigenvectors are arranged like this in be_trial
      VABtemp2(1,0) = n_trial(1,j);
      VABtemp2(2,0) = n_trial(2,j);
      
      ABhelp.MultiplyNN(AB2,VABtemp1) ;
      AB(i,j) = ABhelp(0,0)*VABtemp2(0,0)+ABhelp(1,0)*VABtemp2(1,0)+ABhelp(2,0)*VABtemp2(2,0);
    }
  }

  BSold(0,0) = AB(0,0); BSold(1,0) = AB(1,1);  BSold(2,0) = AB(2,2); 
  
  // Check for yielding
  dev_tau_abs = 0.0;
  for (int i = 0; i < 3; ++i) dev_tau_abs += (dev_KH_trial(i,0)-BSold(i,0) )*(dev_KH_trial(i,0)-BSold(i,0));
  dev_tau_abs = sqrt(dev_tau_abs);
  
  f = dev_tau_abs - sqrt(2.0/3.0) * (H*epsilon_p + KH_y0 + (sigmainf - KH_y0) * (1.0 - exp(-deltad*epsilon_p)));
  
  // Radial return algorithm
  //**********************************************************************
  if (f > 0.0) // Plastic
  //**********************************************************************
  {
    LINALG::Matrix<3,1> tempdev;
    tempdev(0,0) = dev_KH_trial(0,0)-BSold(0,0);
    tempdev(1,0) = dev_KH_trial(1,0)-BSold(1,0);
    tempdev(2,0) = dev_KH_trial(2,0)-BSold(2,0);
    Vec.Update(1.0/(dev_tau_abs),tempdev);
    
    dj = 0.0;
    double Rdj = 1.0;
    double Rtan;
    double ddj;
    double teller = 0.0;
    
    for (int poep = 0; poep<1000; poep++)
    {
      teller = teller + 1.0;
      
      Rdj = f - 2.0*mue*dj - 2.0/3.0*H*dj - 2.0/3.0*Hkin*dj - sqrt(2.0/3.0)*(sigmainf-KH_y0)
          * (-exp(-deltad*(epsilon_p+sqrt(2.0/3.0)*dj))+exp(-deltad*epsilon_p));
      
      Rtan = -2.0*mue-2.0/3.0*H- 2.0/3.0*Hkin-2.0/3.0*deltad*(sigmainf-KH_y0)
          * (exp(-deltad*(epsilon_p+sqrt(2.0/3.0)*dj)));
      
      ddj = -Rdj/Rtan;
      
      dj = dj+ddj;
      
      Rdj = f - 2.0*mue*dj - 2.0/3.0*H*dj - 2.0/3.0*Hkin*dj - sqrt(2.0/3.0)*(sigmainf-KH_y0)
          * (-exp(-deltad*(epsilon_p+sqrt(2.0/3.0)*dj))+exp(-deltad*epsilon_p));
    }
    
    for (int i = 0; i<3; i++) Le(i,0) = exp(log(Le_trial(i,0)) - 2.0*dj*Vec(i,0));
    for (int i = 0; i<3; i++) dev_KH(i,0) = dev_KH_trial(i,0) -2.0*mue*dj*Vec(i,0);
  }
  //**********************************************************************
  else // Elastic
  //**********************************************************************
  {
    for (int i = 0; i<3; i++)
    {
      Le(i,0) = Le_trial(i,0);
      dev_KH(i,0) = dev_KH_trial(i,0);
      dj=0.0;
    }
  }
  
  // Update inverse of elastic left Cauchy-Green tensor
  for (int i = 0; i<3; i++)
  {
    for (int j = 0; j<3; j++)
    {
      for (int k = 0; k<3; k++)
      {
        be(i,j) = 0.0 ; /// Determine PK2
      }
    }
  }
  
  for (int i = 0; i<3; i++)
  {
    for (int j = 0; j<3; j++)
    {
      for (int k = 0; k<3; k++)
      {
        be(i,j) += Le(k,0)*n_trial(i,k)*n_trial(j,k) ; /// Determine PK2
      }
    }
  }
  
  //  Update stress
  
  /// New deviatoric PK2 stresses in principle directions
  dev_PK2_principle.EDivide(dev_KH,Le);
  
  /// New PK2 stresses in principle directions
  for (int i = 0; i<3; i++) PK2_principle(i,0) = dev_PK2_principle(i,0) + p*J/Le(i,0); 
  
  // Calculation of F^-1 (Finv)
  LINALG::Matrix<3,3> Finv(*defgrd);
  Finv.Invert();
  
  LINALG::Matrix<3,1> Vtemp1;
  LINALG::Matrix<3,1> Vtemp2;

  for (int i = 0; i<3; i++) 
  {
    // check if eigenvectors are arranged like this in be_trial
    Vtemp1(0,0) = n_trial(0,i); 
    Vtemp1(1,0) = n_trial(1,i);
    Vtemp1(2,0) = n_trial(2,i);
      
    LINALG::Matrix<3,1> Ntemp;
    Ntemp.Multiply(sqrt(Le(i,0)),Finv,Vtemp1) ;
    
    N_trial(0,i) = Ntemp(0,0);
    N_trial(1,i) = Ntemp(1,0);
    N_trial(2,i) = Ntemp(2,0);
  }
  
  /// Determine PK2
  for (int i = 0; i<3; i++)
  {
    for (int j = 0; j<3; j++)
    {
      for (int k = 0; k<3; k++)
      {
        PK2(i,j) = 0.0;
      }
    }
  }
  
  for (int i = 0; i<3; i++)
  {
    for (int j = 0; j<3; j++)
    {
      for (int k = 0; k<3; k++)
      {
        PK2(i,j) += PK2_principle(k,0)*N_trial(i,k)*N_trial(j,k);
      }
    }
  }

  // Transfer PK2 tensor to stress vector
  (*stress)(0) = PK2(0,0);
  (*stress)(1) = PK2(1,1);
  (*stress)(2) = PK2(2,2);
  (*stress)(3) = 0.5*(PK2(0,1)+PK2(1,0));
  (*stress)(4) = 0.5*(PK2(1,2)+PK2(2,1));
  (*stress)(5) = 0.5*(PK2(0,2)+PK2(2,0));
  
  LINALG::Matrix<3,3> temp1;
  LINALG::Matrix<3,3> temp2;
  
  // Update tangent modulus
  LINALG::Matrix<6,6> ET(false);
  for (int a = 0; a<6; a++)
  {
    for (int b = 0; b<6; b++)
    {
      ET(a,b) =0.0;
    }
  }

  for (int a = 0; a < 3; a++)
  {
    for (int b = 0; b < 3; b++)
    {
      if (f>0.0)
      {
        cab = 1.0/(Le(a,0)*Le(b,0))*((1.0-2.0*mue*dj/(dev_tau_abs))*(2.0*mue*(a==b)-2.0/3.0*mue)
             -2.0*mue*Vec(a,0)*Vec(b,0)*(2.0*mue/(2.0*mue+2.0/3.0*H+2.0/3.0*Hkin+2.0/3.0*deltad
             *(sigmainf-KH_y0)*(exp(-deltad*(epsilon_p+sqrt(2.0/3.0)*dj))) )-2.0*mue*dj/dev_tau_abs));
        cbb = 1.0/(Le(a,0)*Le(b,0))*((1.0-2.0*mue*dj/(dev_tau_abs))*(2.0*mue*(b==b)-2.0/3.0*mue)
             -2.0*mue*Vec(b,0)*Vec(b,0)*(2.0*mue/(2.0*mue+2.0/3.0*H+2.0/3.0*Hkin+2.0/3.0*deltad
             *(sigmainf-KH_y0)*(exp(-deltad*(epsilon_p+sqrt(2.0/3.0)*dj))) )-2.0*mue*dj/dev_tau_abs));
      }
        
      else
      {
        cab = 1.0/(Le(a,0)*Le(b,0))*(2.0*mue*(a==b)-2.0/3.0*mue);
        cbb = 1.0/(Le(a,0)*Le(b,0))*(2.0*mue*(b==b)-2.0/3.0*mue);
      }
      
      // check if eigenvectors are arranged like this in be_trial
      Vtemp1(0,0) = N_trial(0,a); 
      Vtemp1(1,0) = N_trial(1,a);
      Vtemp1(2,0) = N_trial(2,a);
      
      // check if eigenvectors are arranged like this in be_trial
      Vtemp2(0,0) = N_trial(0,b);
      Vtemp2(1,0) = N_trial(1,b);
      Vtemp2(2,0) = N_trial(2,b);
                  
      temp1.MultiplyNT(Vtemp1,Vtemp1);
      temp2.MultiplyNT(Vtemp2,Vtemp2);
      
      ElastSymTensorMultiply(ET, cab, temp1, temp2, 1.0);
  
      if (a==b)
      {
        ElastSymTensorMultiply(ET, -2.0*dev_PK2_principle(a,0)/Le(a,0), temp1, temp1, 1.0);
      }
      
      if (a!=b)
      {
        if (Le_trial(a,0) != Le_trial(b,0))
        {       
          c1 = (dev_PK2_principle(a,0) * Le_trial(b,0)/Le(b,0) - dev_PK2_principle(b,0)
              * Le_trial(a,0)/Le(a,0)) / (Le_trial(a,0) -Le_trial(b,0));
        }
        
        if (Le_trial(a,0) == Le_trial(b,0))
        {
          c1 = 0.5*cbb*(Le_trial(a,0)/Le_trial(b,0))-0.5*cab-dev_PK2_principle(a,0)/Le(b,0);
        }
        
        temp1.MultiplyNT(Vtemp1,Vtemp2);
        temp2.MultiplyNT(Vtemp2,Vtemp1);

        ElastSymTensorMultiply(ET, c1, temp1, temp1, 1.0);
        ElastSymTensorMultiply(ET, c1, temp1, temp2, 1.0);
      }
    }
  }
  
  double c2 = kappa*(1.0-log(J));
  double c3 = p*J;
  
  // add scalar Cinv x Cinv
  for (int i=0; i<6; ++i)
  {
    for (int j=0; j<6; ++j)
    {
      (ET)(i,j) += c2 * Cinv(i) * Cinv(j) +  c3 * Cinv(i) * Cinv(j);
    }
  }
  
  AddtoCmatHolzapfelProduct(ET,Cinv,-c3*2.0);
  
  for (int i=0; i<6; ++i)
  {
    for (int j=0; j<6; ++j)
    {
      (*cmat)(i,j) =ET(i,j);
    }
  }

  // Update state variables
  LINALG::Matrix<3,3> temp;
  temp.MultiplyNT(be,Finv);
  LINALG::Matrix<3,3> temp4;
  temp4.MultiplyNN(Finv, temp); /// Update plastic right Cauchy-Green tensor
  
  if  (f>0.0) 
  {
    //cout << " plastic ";
    (histplasticrcgcurr_->at(gp))(0,0) = temp4(0,0);
    (histplasticrcgcurr_->at(gp))(1,0) = temp4(1,1);
    (histplasticrcgcurr_->at(gp))(2,0) = temp4(2,2);
    (histplasticrcgcurr_->at(gp))(3,0) = 0.5*(temp4(0,1)+temp4(1,0));
    (histplasticrcgcurr_->at(gp))(4,0) = 0.5*(temp4(2,1)+temp4(1,2));
    (histplasticrcgcurr_->at(gp))(5,0) = 0.5*(temp4(2,0)+temp4(0,2));
  
    (histplasticbstresscurr_->at(gp))(0,0) = BSold(0,0)+2.0/3.0*dj*Hkin*Vec(0,0);
    (histplasticbstresscurr_->at(gp))(1,0) = BSold(1,0)+2.0/3.0*dj*Hkin*Vec(1,0);
    (histplasticbstresscurr_->at(gp))(2,0) = BSold(2,0)+2.0/3.0*dj*Hkin*Vec(2,0);
  
    (lamdacurr_->at(gp))(0,0) = Le_trial(0,0); /*added*/
    (lamdacurr_->at(gp))(1,0) = Le_trial(1,0); /*added*/
    (lamdacurr_->at(gp))(2,0) = Le_trial(2,0); /*added*/
  }

  else
  {
    //cout << " elastic ";
    (histplasticrcgcurr_->at(gp))(0,0) = (histplasticrcglast_->at(gp))(0,0);
    (histplasticrcgcurr_->at(gp))(1,0) = (histplasticrcglast_->at(gp))(1,0);
    (histplasticrcgcurr_->at(gp))(2,0) = (histplasticrcglast_->at(gp))(2,0);
    (histplasticrcgcurr_->at(gp))(3,0) = (histplasticrcglast_->at(gp))(3,0);
    (histplasticrcgcurr_->at(gp))(4,0) = (histplasticrcglast_->at(gp))(4,0);
    (histplasticrcgcurr_->at(gp))(5,0) = (histplasticrcglast_->at(gp))(5,0);
  
    (histplasticbstresscurr_->at(gp))(0,0) = BSold(0,0);
    (histplasticbstresscurr_->at(gp))(1,0) = BSold(1,0);
    (histplasticbstresscurr_->at(gp))(2,0) = BSold(2,0);
  
    (lamdacurr_->at(gp))(0,0) = Le_trial(0,0); /*added*/
    (lamdacurr_->at(gp))(1,0) = Le_trial(1,0); /*added*/
    (lamdacurr_->at(gp))(2,0) = Le_trial(2,0); /*added*/
  }

  // back stress
  LINALG::Matrix<3,3> Backstress;

  for (int i = 0; i<3; i++)
  {
    for (int j = 0; j<3; j++)
    {
      for (int k = 0; k<3; k++)
      {
        Backstress(i,j) = 0.0 ; /// Determine PK2
      }
    }
  }

  for (int i = 0; i<3; i++)
  {
    for (int j = 0; j<3; j++)
    {
      for (int k = 0; k<3; k++)
      {
        Backstress(i,j)+= (histplasticbstresscurr_->at(gp))(k,0)*n_trial(i,k)*n_trial(j,k) ; /// Determine PK2
      }
    }
  }

  (nsetcurr_->at(gp))(0,0) = Backstress(0,0);
  (nsetcurr_->at(gp))(0,1) = Backstress(0,1); 
  (nsetcurr_->at(gp))(0,2) = Backstress(0,2);
  (nsetcurr_->at(gp))(1,0) = Backstress(1,0);
  (nsetcurr_->at(gp))(1,1) = Backstress(1,1); 
  (nsetcurr_->at(gp))(1,2) = Backstress(1,2);
  (nsetcurr_->at(gp))(2,0) = Backstress(2,0);
  (nsetcurr_->at(gp))(2,1) = Backstress(2,1); 
  (nsetcurr_->at(gp))(2,2) = Backstress(2,2);

  // equivalent plastic strain
  (histeplasticscurr_->at(gp))(0,0) = epsilon_p + sqrt(2.0/3.0)*dj;

  return;
}

#endif
