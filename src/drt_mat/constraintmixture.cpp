/*!----------------------------------------------------------------------
\file constraintmixture.cpp
\brief
This file contains routines for constraint mixture growth and remodeling.
example input line
MAT 1 MAT_ConstraintMixture DENS 0.001 MUE 1.0E3 PHIE 0.08 PREELA 1.0
K1 1.0 K2 1.0 PRECOLL 1.0 KAPPA 1.0E4 LIFETIME 5.0
HOMSTR 0.0 GROWTHFAC 0.0 STARTTIME 100.0 EXPLICIT 1 TOL 1.0E-8

Here an approach for growth and remodeling of an artery is modeled.
For a detailed description see:
- Humphrey, J. & Rajagopal, K.: A constrained mixture model for arterial
  adaptations to a sustained step change in blood flow,
  Biomechanics and Modeling in Mechanobiology, 2003, 2, 109-126

<pre>
Maintainer: Susanna Tinkl
            tinkl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15265
</pre>
*----------------------------------------------------------------------*/

#ifdef CCADISCRET

#include <vector>
#include "constraintmixture.H"
#include "../drt_lib/drt_globalproblem.H"
#include "matpar_bundle.H"
#include "../drt_lib/drt_linedefinition.H"
#include "../drt_io/io_gmsh.H" // for debug plotting with gmsh
#include "../drt_io/io_control.H" // for debug plotting with gmsh
#include "contchainnetw.H" // for debug plotting with gmsh

/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
MAT::PAR::ConstraintMixture::ConstraintMixture(
  Teuchos::RCP<MAT::PAR::Material> matdata
  )
: Parameter(matdata),
  density_(matdata->GetDouble("DENS")),
  mue_(matdata->GetDouble("MUE")),
  phielastin_(matdata->GetDouble("PHIE")),
  prestretchelastin_(matdata->GetDouble("PREELA")),
  k1_(matdata->GetDouble("K1")),
  k2_(matdata->GetDouble("K2")),
  prestretchcollagen_(matdata->GetDouble("PRECOLL")),
  kappa_(matdata->GetDouble("KAPPA")),
  lifetime_(matdata->GetDouble("LIFETIME")),
  homstress_(matdata->GetDouble("HOMSTR")),
  growthfactor_(matdata->GetDouble("GROWTHFAC")),
  starttime_(matdata->GetDouble("STARTTIME")),
  explicit_(matdata->GetInt("EXPLICIT")),
  abstol_(matdata->GetDouble("TOL"))
{
}


Teuchos::RCP<MAT::Material> MAT::PAR::ConstraintMixture::CreateMaterial()
{
  return Teuchos::rcp(new MAT::ConstraintMixture(this));
}


MAT::ConstraintMixtureType MAT::ConstraintMixtureType::instance_;

DRT::ParObject* MAT::ConstraintMixtureType::Create( const std::vector<char> & data )
{
  MAT::ConstraintMixture* comix = new MAT::ConstraintMixture();
  comix->Unpack(data);
  return comix;
}

/*----------------------------------------------------------------------*
 |  Constructor                                   (public)         12/10|
 *----------------------------------------------------------------------*/
MAT::ConstraintMixture::ConstraintMixture()
  : params_(NULL)
{
}


/*----------------------------------------------------------------------*
 |  Copy-Constructor                             (public)          12/10|
 *----------------------------------------------------------------------*/
MAT::ConstraintMixture::ConstraintMixture(MAT::PAR::ConstraintMixture* params)
  : params_(params)
{
}

/*----------------------------------------------------------------------*
 |  Pack                                          (public)         12/10|
 *----------------------------------------------------------------------*/
void MAT::ConstraintMixture::Pack(vector<char>& data) const
{
  data.resize(0);

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);

  // matid
  int matid = -1;
  if (params_ != NULL) matid = params_->Id();  // in case we are in post-process mode
  AddtoPack(data,matid);

  int numgp;
  if (!isinit_)
  {
    numgp = 0; // not initialized -> nothing to pack
  }
  else
  {
    numgp = a1_->size();   // size is number of gausspoints
  }
  AddtoPack(data,numgp);
  // Pack internal variables
  for (int gp = 0; gp < numgp; ++gp)
  {
    AddtoPack(data,a1_->at(gp));
    AddtoPack(data,a2_->at(gp));
    AddtoPack(data,a3_->at(gp));
    AddtoPack(data,a4_->at(gp));
    AddtoPack(data,collagenstretch1_->at(gp));
    AddtoPack(data,collagenstretch2_->at(gp));
    AddtoPack(data,collagenstretch3_->at(gp));
    AddtoPack(data,collagenstretch4_->at(gp));
    AddtoPack(data,massprod1_->at(gp));
    AddtoPack(data,massprod2_->at(gp));
    AddtoPack(data,massprod3_->at(gp));
    AddtoPack(data,massprod4_->at(gp));
    AddtoPack(data,vismassstress_->at(gp));
  }
  if (numgp > 0) AddtoPack(data,massprodbasal_);

  return;
}


/*----------------------------------------------------------------------*
 |  Unpack                                        (public)         12/10|
 *----------------------------------------------------------------------*/
void MAT::ConstraintMixture::Unpack(const vector<char>& data)
{
  isinit_=true;
  vector<char>::size_type position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position,data,type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");

  // matid and recover params_
  int matid;
  ExtractfromPack(position,data,matid);
  if (DRT::Problem::Instance()->Materials() != Teuchos::null)
  {
    const int probinst = DRT::Problem::Instance()->Materials()->GetReadFromProblem();
    MAT::PAR::Parameter* mat = DRT::Problem::Instance(probinst)->Materials()->ParameterById(matid);
    if (mat->Type() == MaterialType())
      params_ = static_cast<MAT::PAR::ConstraintMixture*>(mat);
    else
      dserror("Type of parameter material %d does not fit to calling type %d", mat->Type(), MaterialType());
  }
  else
  {
    params_ = NULL;
  }

  int numgp;
  ExtractfromPack(position,data,numgp);
  if (numgp == 0){ // no history data to unpack
    isinit_=false;
    if (position != data.size())
      dserror("Mismatch in size of data %d <-> %d",data.size(),position);
    return;
  }
  // unpack fiber internal variables
  a1_ = Teuchos::rcp(new vector<LINALG::Matrix<3,1> >(numgp));
  a2_ = Teuchos::rcp(new vector<LINALG::Matrix<3,1> >(numgp));
  a3_ = Teuchos::rcp(new vector<LINALG::Matrix<3,1> >(numgp));
  a4_ = Teuchos::rcp(new vector<LINALG::Matrix<3,1> >(numgp));
  collagenstretch1_ = Teuchos::rcp(new vector<vector<double> > (numgp));
  collagenstretch2_ = Teuchos::rcp(new vector<vector<double> > (numgp));
  collagenstretch3_ = Teuchos::rcp(new vector<vector<double> > (numgp));
  collagenstretch4_ = Teuchos::rcp(new vector<vector<double> > (numgp));
  massprod1_ = Teuchos::rcp(new vector<vector<double> > (numgp));
  massprod2_ = Teuchos::rcp(new vector<vector<double> > (numgp));
  massprod3_ = Teuchos::rcp(new vector<vector<double> > (numgp));
  massprod4_ = Teuchos::rcp(new vector<vector<double> > (numgp));
  vismassstress_ = Teuchos::rcp(new vector<double> (numgp));

  for (int gp = 0; gp < numgp; ++gp) {
    LINALG::Matrix<3,1> alin;
    ExtractfromPack(position,data,alin);
    a1_->at(gp) = alin;
    ExtractfromPack(position,data,alin);
    a2_->at(gp) = alin;
    ExtractfromPack(position,data,alin);
    a3_->at(gp) = alin;
    ExtractfromPack(position,data,alin);
    a4_->at(gp) = alin;
    vector<double> a;
    ExtractfromPack(position,data,a);
    collagenstretch1_->at(gp) = a;
    ExtractfromPack(position,data,a);
    collagenstretch2_->at(gp) = a;
    ExtractfromPack(position,data,a);
    collagenstretch3_->at(gp) = a;
    ExtractfromPack(position,data,a);
    collagenstretch4_->at(gp) = a;
    ExtractfromPack(position,data,a);
    massprod1_->at(gp) = a;
    ExtractfromPack(position,data,a);
    massprod2_->at(gp) = a;
    ExtractfromPack(position,data,a);
    massprod3_->at(gp) = a;
    ExtractfromPack(position,data,a);
    massprod4_->at(gp) = a;
    double b;
    ExtractfromPack(position,data,b);
    vismassstress_->at(gp) = b;
  }
  numpast_ = massprod1_->at(0).size();
  double basal;
  ExtractfromPack(position,data,basal);
  massprodbasal_=basal;

  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d",data.size(),position);

  return;
}

/*----------------------------------------------------------------------*
 |  Setup                                         (public)         12/10|
 *----------------------------------------------------------------------*/
void MAT::ConstraintMixture::Setup (const int numgp, DRT::INPUT::LineDefinition* linedef)
{
  // history variables
  collagenstretch1_ = Teuchos::rcp(new vector<vector<double> > (numgp));
  collagenstretch2_ = Teuchos::rcp(new vector<vector<double> > (numgp));
  collagenstretch3_ = Teuchos::rcp(new vector<vector<double> > (numgp));
  collagenstretch4_ = Teuchos::rcp(new vector<vector<double> > (numgp));
  massprod1_ = Teuchos::rcp(new vector<vector<double> > (numgp));
  massprod2_ = Teuchos::rcp(new vector<vector<double> > (numgp));
  massprod3_ = Teuchos::rcp(new vector<vector<double> > (numgp));
  massprod4_ = Teuchos::rcp(new vector<vector<double> > (numgp));
  vismassstress_ = Teuchos::rcp(new vector<double> (numgp));

  // I don't like this way of getting dt!!!
  const Teuchos::ParameterList& timeintegr = DRT::Problem::Instance()->StructuralDynamicParams();
  double dt = timeintegr.get<double>("TIMESTEP");
  numpast_ = static_cast<int>(round(params_->lifetime_ / dt)) + params_->explicit_;
  // basal mass production rate determined by DENS, PHIE and LIFETIME
  massprodbasal_ = (1 - params_->phielastin_) * params_->density_ / 4.0 / params_->lifetime_;
  for (int gp = 0; gp < numgp; gp++)
  {
    collagenstretch1_->at(gp).resize(numpast_);
    collagenstretch2_->at(gp).resize(numpast_);
    collagenstretch3_->at(gp).resize(numpast_);
    collagenstretch4_->at(gp).resize(numpast_);
    massprod1_->at(gp).resize(numpast_);
    massprod2_->at(gp).resize(numpast_);
    massprod3_->at(gp).resize(numpast_);
    massprod4_->at(gp).resize(numpast_);
    for (int idpast = 0; idpast < numpast_; idpast++)
    {
      collagenstretch1_->at(gp)[idpast] = 1.0;
      collagenstretch2_->at(gp)[idpast] = 1.0;
      collagenstretch3_->at(gp)[idpast] = 1.0;
      collagenstretch4_->at(gp)[idpast] = 1.0;
      massprod1_->at(gp)[idpast] = massprodbasal_;
      massprod2_->at(gp)[idpast] = massprodbasal_;
      massprod3_->at(gp)[idpast] = massprodbasal_;
      massprod4_->at(gp)[idpast] = massprodbasal_;
    }
    vismassstress_->at(gp) = 0.0;
  }

  // fiber vectors
  a1_ = Teuchos::rcp(new vector<LINALG::Matrix<3,1> > (numgp));
  a2_ = Teuchos::rcp(new vector<LINALG::Matrix<3,1> > (numgp));
  a3_ = Teuchos::rcp(new vector<LINALG::Matrix<3,1> > (numgp));
  a4_ = Teuchos::rcp(new vector<LINALG::Matrix<3,1> > (numgp));

  // read local (cylindrical) cosy-directions at current element
  vector<double> rad;
  vector<double> axi;
  vector<double> cir;
  linedef->ExtractDoubleVector("RAD",rad);
  linedef->ExtractDoubleVector("AXI",axi);
  linedef->ExtractDoubleVector("CIR",cir);

  LINALG::Matrix<3,3> locsys;
  // basis is local cosy with third vec e3 = circumferential dir and e2 = axial dir
  double radnorm=0.; double axinorm=0.; double cirnorm=0.;
  for (int i = 0; i < 3; i++) {
  radnorm += rad[i]*rad[i]; axinorm += axi[i]*axi[i]; cirnorm += cir[i]*cir[i];
  }
  radnorm = sqrt(radnorm); axinorm = sqrt(axinorm); cirnorm = sqrt(cirnorm);
  for (int i=0; i<3; i++){
    locsys(i,0) = rad[i]/radnorm;
    locsys(i,1) = axi[i]/axinorm;
    locsys(i,2) = cir[i]/cirnorm;
  }

  const double gamma = (45*PI)/180.; //angle for diagonal fibers

  for (int gp = 0; gp < numgp; gp++)
  {
    for (int i = 0; i < 3; i++) {
      // a1 = e3, circumferential direction
      a1_->at(gp)(i) = locsys(i,2);
      // a2 = e2, axial direction
      a2_->at(gp)(i) = locsys(i,1);
      // a3 = cos gamma e3 + sin gamma e2
      a3_->at(gp)(i) = cos(gamma)*locsys(i,2) + sin(gamma)*locsys(i,1);
      // a4 = cos gamma e3 - sin gamma e2
      a4_->at(gp)(i) = cos(gamma)*locsys(i,2) - sin(gamma)*locsys(i,1);
    }
  }

  isinit_ = true;
}

/*----------------------------------------------------------------------*
 |  Update internal variables                     (public)         12/10|
 *----------------------------------------------------------------------*/
void MAT::ConstraintMixture::Update()
{
  int numpast = numpast_;
  int numgp = massprod1_->size();
  for (int gp = 0; gp < numgp; gp++)
  {
    for (int idpast =  numpast-1; idpast > 0; idpast--)
    {
       collagenstretch1_->at(gp)[idpast] = collagenstretch1_->at(gp)[idpast-1];
       collagenstretch2_->at(gp)[idpast] = collagenstretch2_->at(gp)[idpast-1];
       collagenstretch3_->at(gp)[idpast] = collagenstretch3_->at(gp)[idpast-1];
       collagenstretch4_->at(gp)[idpast] = collagenstretch4_->at(gp)[idpast-1];
       massprod1_->at(gp)[idpast] = massprod1_->at(gp)[idpast-1];
       massprod2_->at(gp)[idpast] = massprod2_->at(gp)[idpast-1];
       massprod3_->at(gp)[idpast] = massprod3_->at(gp)[idpast-1];
       massprod4_->at(gp)[idpast] = massprod4_->at(gp)[idpast-1];
       //cout << idpast << endl;
    }
    // nothing to do for collagenstretch as the actual fiber stretch is already stored there
    massprod1_->at(gp)[0] = massprodbasal_;
    massprod2_->at(gp)[0] = massprodbasal_;
    massprod3_->at(gp)[0] = massprodbasal_;
    massprod4_->at(gp)[0] = massprodbasal_;
  }
}

/*----------------------------------------------------------------------*
 |  Evaluate                                      (public)         12/10|
 *----------------------------------------------------------------------*/
void MAT::ConstraintMixture::Evaluate
(
  const LINALG::Matrix<NUM_STRESS_3D,1>* glstrain,
  const int gp,
  LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D>* cmat,
  LINALG::Matrix<NUM_STRESS_3D,1>* stress,
  double dt,
  double time,
  bool output
)
{
  // differentiate between explicit and implicit description
  bool explintegration = params_->explicit_;
//  int numsteps;
  int firstiter;
  if (explintegration == 1)
  {
//    numsteps = numpast -1;
    firstiter = 1;
  } else {
//    numsteps = numpast;
    firstiter = 0;
  }

  //--------------------------------------------------------------------------------------
  // build identity tensor I
  LINALG::Matrix<NUM_STRESS_3D,1> Id(true);
  for (int i = 0; i < 3; i++) Id(i) = 1.0;
  // right Cauchy-Green Tensor  C = 2 * E + I
  LINALG::Matrix<NUM_STRESS_3D,1> C(*glstrain);
  C.Scale(2.0);
  C += Id;

  //--------------------------------------------------------------------------------------
  // compute actual collagen stretches and store them
  double eps = 1.0e-12;
  if (time > params_->starttime_ + eps)
  {
    double actcollagenstretch =  a1_->at(gp)(0)*a1_->at(gp)(0)*C(0) + a1_->at(gp)(1)*a1_->at(gp)(1)*C(1)
          + a1_->at(gp)(2)*a1_->at(gp)(2)*C(2) + a1_->at(gp)(0)*a1_->at(gp)(1)*C(3)
          + a1_->at(gp)(1)*a1_->at(gp)(2)*C(4) + a1_->at(gp)(0)*a1_->at(gp)(2)*C(5); // = trace(A1:C)
    actcollagenstretch = sqrt(actcollagenstretch);
    collagenstretch1_->at(gp)[0] = actcollagenstretch;
    actcollagenstretch =  a2_->at(gp)(0)*a2_->at(gp)(0)*C(0) + a2_->at(gp)(1)*a2_->at(gp)(1)*C(1)
          + a2_->at(gp)(2)*a2_->at(gp)(2)*C(2) + a2_->at(gp)(0)*a2_->at(gp)(1)*C(3)
          + a2_->at(gp)(1)*a2_->at(gp)(2)*C(4) + a2_->at(gp)(0)*a2_->at(gp)(2)*C(5); // = trace(A2:C)
    actcollagenstretch = sqrt(actcollagenstretch);
    collagenstretch2_->at(gp)[0] = actcollagenstretch;
    actcollagenstretch =  a3_->at(gp)(0)*a3_->at(gp)(0)*C(0) + a3_->at(gp)(1)*a3_->at(gp)(1)*C(1)
          + a3_->at(gp)(2)*a3_->at(gp)(2)*C(2) + a3_->at(gp)(0)*a3_->at(gp)(1)*C(3)
          + a3_->at(gp)(1)*a3_->at(gp)(2)*C(4) + a3_->at(gp)(0)*a3_->at(gp)(2)*C(5); // = trace(A3:C)
    actcollagenstretch = sqrt(actcollagenstretch);
    collagenstretch3_->at(gp)[0] = actcollagenstretch;
    actcollagenstretch =  a4_->at(gp)(0)*a4_->at(gp)(0)*C(0) + a4_->at(gp)(1)*a4_->at(gp)(1)*C(1)
          + a4_->at(gp)(2)*a4_->at(gp)(2)*C(2) + a4_->at(gp)(0)*a4_->at(gp)(1)*C(3)
          + a4_->at(gp)(1)*a4_->at(gp)(2)*C(4) + a4_->at(gp)(0)*a4_->at(gp)(2)*C(5); // = trace(A4:C)
    actcollagenstretch = sqrt(actcollagenstretch);
    collagenstretch4_->at(gp)[0] = actcollagenstretch;
  }

  //--------------------------------------------------------------------------------------
  // some variables
  //int numpast = numpast_;
  //double prestretchelastin = params_->prestretchelastin_;
  double prestretchcollagen = params_->prestretchcollagen_;
  double density = params_->density_;
  //double currmassdens = 0.0;

  EvaluateStress(glstrain, gp, cmat, stress, firstiter, dt, time);

  if (time > params_->starttime_ + eps && params_->growthfactor_ != 0.0)
  {
    //--------------------------------------------------------------------------------------
    // compute new deposition rates
    LINALG::Matrix<3,3> Cmatrix;
    LINALG::Matrix<3,3> Smatrix;
    LINALG::Matrix<4,1> massstress;
    LINALG::Matrix<4,1> massprodcomp;
    MassProduction(gp,C,&Cmatrix,*stress,&Smatrix,&massstress,&massprodcomp);

    // start values for local Newton iteration are computed
    // distinguish between explicit and implicit integration
    if (explintegration == 1.0)
    {
      massprod1_->at(gp)[0] = massprodcomp(0);
      massprod2_->at(gp)[0] = massprodcomp(1);
      massprod3_->at(gp)[0] = massprodcomp(2);
      massprod4_->at(gp)[0] = massprodcomp(3);
      vismassstress_->at(gp) = massstress(1);
    } else {
      LINALG::Matrix<10,1> Residual(true);
      // first six components are zero as the start value for stress is the computed value
      Residual(6) = massprod1_->at(gp)[0] - massprodcomp(0);
      Residual(7) = massprod2_->at(gp)[0] - massprodcomp(1);
      Residual(8) = massprod3_->at(gp)[0] - massprodcomp(2);
      Residual(9) = massprod4_->at(gp)[0] - massprodcomp(3);

      LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D> cmatelastic(*cmat);
      // isotropic invariant
      const double I3 = C(0)*C(1)*C(2)
            + 0.25 * C(3)*C(4)*C(5)
            - 0.25 * C(1)*C(5)*C(5)
            - 0.25 * C(2)*C(3)*C(3)
            - 0.25 * C(0)*C(4)*C(4);    // 3rd invariant, determinant
      const double J = sqrt(I3);     // determinant of F
      //-------------------------------------
      // invert C
      LINALG::Matrix<NUM_STRESS_3D,1> Cinv(6);
      Cinv(0) = C(1)*C(2) - 0.25*C(4)*C(4);
      Cinv(1) = C(0)*C(2) - 0.25*C(5)*C(5);
      Cinv(2) = C(0)*C(1) - 0.25*C(3)*C(3);
      Cinv(3) = 0.25*C(5)*C(4) - 0.5*C(3)*C(2);
      Cinv(4) = 0.25*C(3)*C(5) - 0.5*C(0)*C(4);
      Cinv(5) = 0.25*C(3)*C(4) - 0.5*C(5)*C(1);
      Cinv.Scale(1.0/I3);
      double qdegrad = 0.0;
      Degradation(0.0, &qdegrad);
      // include the case homstress = 0.0!!
      double homstressfixed = 1.0;
      if (params_->homstress_ != 0.0) homstressfixed = params_->homstress_;

      //--------------------------------------------------------------------------------------
      // local Newton iteration
      int localistep = 0;
      int maxstep = 50;
      while (Residual.Norm2() > params_->abstol_ && localistep < maxstep)
      {
        localistep += 1;
        //--------------------------------------------------------------------------------------
        // derivative of residual
        LINALG::Matrix<10,10> DResidual(true);
        for (int id = 0; id < 10; id++) DResidual(id,id) = 1.0;

        // derivative of stress part with respect to mass production
        // do it for all 4 fiber families
        double stretch = prestretchcollagen/collagenstretch1_->at(gp)[0];
        LINALG::Matrix<NUM_STRESS_3D,1> Saniso1(true);
        EvaluateFiber(glstrain, NULL, & Saniso1, a1_->at(gp), stretch);
        double facS = pow(stretch,2) * qdegrad / density * dt;
        Saniso1.Scale(facS);
        for (int id = 0; id < 6; id++) DResidual(id,6) = - Saniso1(id) + dt / density * qdegrad * params_->kappa_ * J * Cinv(id);

        stretch = prestretchcollagen/collagenstretch2_->at(gp)[0];
        LINALG::Matrix<NUM_STRESS_3D,1> Saniso2(true);
        EvaluateFiber(glstrain, NULL, & Saniso2, a2_->at(gp), stretch);
        facS = pow(stretch,2) * qdegrad / density * dt;
        Saniso2.Scale(facS);
        for (int id = 0; id < 6; id++) DResidual(id,7) = - Saniso2(id) + dt / density * qdegrad * params_->kappa_ * J * Cinv(id);

        stretch = prestretchcollagen/collagenstretch3_->at(gp)[0];
        LINALG::Matrix<NUM_STRESS_3D,1> Saniso3(true);
        EvaluateFiber(glstrain, NULL, & Saniso3, a3_->at(gp), stretch);
        facS = pow(stretch,2) * qdegrad / density * dt;
        Saniso3.Scale(facS);
        for (int id = 0; id < 6; id++) DResidual(id,8) = - Saniso3(id) + dt / density * qdegrad * params_->kappa_ * J * Cinv(id);

        stretch = prestretchcollagen/collagenstretch4_->at(gp)[0];
        LINALG::Matrix<NUM_STRESS_3D,1> Saniso4(true);
        EvaluateFiber(glstrain, NULL, & Saniso4, a4_->at(gp), stretch);
        facS = pow(stretch,2) * qdegrad / density * dt;
        Saniso4.Scale(facS);
        for (int id = 0; id < 6; id++) DResidual(id,9) = - Saniso4(id) + dt / density * qdegrad * params_->kappa_ * J * Cinv(id);

        // derivative of mass production with respect to stress
        LINALG::Matrix<3,3> CSC(true);
        LINALG::Matrix<3,3> temp1(true);
        temp1.Multiply(Smatrix,Cmatrix);
        CSC.Multiply(Cmatrix,temp1);

        LINALG::Matrix<3,1> Ca(true);
        Ca.Multiply(Cmatrix,a1_->at(gp));
        LINALG::Matrix<3,1> CSCa(true);
        CSCa.Multiply(CSC,a1_->at(gp));
        double fac = - massprodbasal_ * params_->growthfactor_ / homstressfixed /pow(J,2)
                   / massstress(0) / pow(collagenstretch1_->at(gp)[0],2);
        DResidual(6,0) = fac * Ca(0) * CSCa(0);
        DResidual(6,1) = fac * Ca(1) * CSCa(1);
        DResidual(6,2) = fac * Ca(2) * CSCa(2);
        DResidual(6,3) = fac * 0.5 * (Ca(0) * CSCa(1) + Ca(1) * CSCa(0));
        DResidual(6,4) = fac * 0.5 * (Ca(1) * CSCa(2) + Ca(2) * CSCa(1));
        DResidual(6,5) = fac * 0.5 * (Ca(0) * CSCa(2) + Ca(2) * CSCa(0));

        Ca.Scale(0.0);
        Ca.Multiply(Cmatrix,a2_->at(gp));
        CSCa.Scale(0.0);
        CSCa.Multiply(CSC,a2_->at(gp));
        fac = - massprodbasal_ * params_->growthfactor_ / homstressfixed /pow(J,2)
            / massstress(1) / pow(collagenstretch2_->at(gp)[0],2);
        DResidual(7,0) = fac * Ca(0) * CSCa(0);
        DResidual(7,1) = fac * Ca(1) * CSCa(1);
        DResidual(7,2) = fac * Ca(2) * CSCa(2);
        DResidual(7,3) = fac * 0.5 * (Ca(0) * CSCa(1) + Ca(1) * CSCa(0));
        DResidual(7,4) = fac * 0.5 * (Ca(1) * CSCa(2) + Ca(2) * CSCa(1));
        DResidual(7,5) = fac * 0.5 * (Ca(0) * CSCa(2) + Ca(2) * CSCa(0));

        Ca.Scale(0.0);
        Ca.Multiply(Cmatrix,a3_->at(gp));
        CSCa.Scale(0.0);
        CSCa.Multiply(CSC,a3_->at(gp));
        fac = - massprodbasal_ * params_->growthfactor_ / homstressfixed /pow(J,2)
            / massstress(2) / pow(collagenstretch3_->at(gp)[0],2);
        DResidual(8,0) = fac * Ca(0) * CSCa(0);
        DResidual(8,1) = fac * Ca(1) * CSCa(1);
        DResidual(8,2) = fac * Ca(2) * CSCa(2);
        DResidual(8,3) = fac * 0.5 * (Ca(0) * CSCa(1) + Ca(1) * CSCa(0));
        DResidual(8,4) = fac * 0.5 * (Ca(1) * CSCa(2) + Ca(2) * CSCa(1));
        DResidual(8,5) = fac * 0.5 * (Ca(0) * CSCa(2) + Ca(2) * CSCa(0));

        Ca.Scale(0.0);
        Ca.Multiply(Cmatrix,a4_->at(gp));
        CSCa.Scale(0.0);
        CSCa.Multiply(CSC,a4_->at(gp));
        fac = - massprodbasal_ * params_->growthfactor_ / homstressfixed /pow(J,2)
            / massstress(3) / pow(collagenstretch4_->at(gp)[0],2);
        DResidual(9,0) = fac * Ca(0) * CSCa(0);
        DResidual(9,1) = fac * Ca(1) * CSCa(1);
        DResidual(9,2) = fac * Ca(2) * CSCa(2);
        DResidual(9,3) = fac * 0.5 * (Ca(0) * CSCa(1) + Ca(1) * CSCa(0));
        DResidual(9,4) = fac * 0.5 * (Ca(1) * CSCa(2) + Ca(2) * CSCa(1));
        DResidual(9,5) = fac * 0.5 * (Ca(0) * CSCa(2) + Ca(2) * CSCa(0));

        //----------------------------------------------------
        // solve linear system of equations: gradF * incr = -F
        //----------------------------------------------------
        // F = F*-1.0
        Residual.Scale(-1.0);
        LINALG::Matrix<10,1> increment(true);
        // solve A.X=B
        LINALG::FixedSizeSerialDenseSolver<10,10,1> solver;
        solver.SetMatrix(DResidual);              // set A=DResidual
        solver.SetVectors(increment, Residual);           // set X=increment, B=Residual
        solver.FactorWithEquilibration(true); // "some easy type of preconditioning" (Michael)
        int err2 = solver.Factor();           // ?
        int err = solver.Solve();             // X = A^-1 B
        if ((err!=0) || (err2!=0))
          dserror("solving linear system in Newton-Raphson method for implicit integration failed");

        // damping strategy
        double omega = 2.0;
        LINALG::Matrix<6,1> stepstress(true);
        LINALG::Matrix<4,1> oldmass(true);
        oldmass(0) = massprod1_->at(gp)[0];
        oldmass(1) = massprod2_->at(gp)[0];
        oldmass(2) = massprod3_->at(gp)[0];
        oldmass(3) = massprod4_->at(gp)[0];
        LINALG::Matrix<10,1> Residualtemp(Residual);
        double omegamin = 1.0/64.0;
        while (Residualtemp.Norm2() > (1.0-0.5*omega)*Residual.Norm2() && omega > omegamin)
        {
          // update of theta
          omega = omega/2.0;
          for (int id = 0; id < 6; id++) stepstress(id) = (*stress)(id) + omega*increment(id);
          massprod1_->at(gp)[0] = oldmass(0) + omega*increment(6);
          massprod2_->at(gp)[0] = oldmass(1) + omega*increment(7);
          massprod3_->at(gp)[0] = oldmass(2) + omega*increment(8);
          massprod4_->at(gp)[0] = oldmass(3) + omega*increment(9);

          LINALG::Matrix<NUM_STRESS_3D,1> Scomp(true);
          cmatelastic.Scale(0.0);
          EvaluateStress(glstrain, gp, &cmatelastic, &Scomp, firstiter, dt, time);
          MassProduction(gp,C,&Cmatrix,stepstress,&Smatrix,&massstress,&massprodcomp);

          for (int id = 0; id < 6; id++) Residualtemp(id) = stepstress(id) - Scomp(id);
          Residualtemp(6) = massprod1_->at(gp)[0] - massprodcomp(0);
          Residualtemp(7) = massprod2_->at(gp)[0] - massprodcomp(1);
          Residualtemp(8) = massprod3_->at(gp)[0] - massprodcomp(2);
          Residualtemp(9) = massprod4_->at(gp)[0] - massprodcomp(3);
        }
        if (omega <= omegamin && Residualtemp.Norm2() > (1.0-0.5*omega)*Residual.Norm2()) dserror("no damping coefficient found");

        for (int id = 0; id < 6; id++) (*stress)(id) = stepstress(id);
        Residual = Residualtemp;

//        for (int id = 0; id < 6; id++) (*stress)(id) += increment(id);
        Smatrix(0,0) = (*stress)(0);
        Smatrix(0,1) = (*stress)(3);
        Smatrix(0,2) = (*stress)(5);
        Smatrix(1,0) = Smatrix(0,1);
        Smatrix(1,1) = (*stress)(1);
        Smatrix(1,2) = (*stress)(4);
        Smatrix(2,0) = Smatrix(0,2);
        Smatrix(2,1) = Smatrix(1,2);
        Smatrix(2,2) = (*stress)(2);
//        massprod1_->at(gp)[0] += increment(6);
//        massprod2_->at(gp)[0] += increment(7);
//        massprod3_->at(gp)[0] += increment(8);
//        massprod4_->at(gp)[0] += increment(9);

        if ((massprod1_->at(gp)[0] < 0.0) || (massprod2_->at(gp)[0] < 0.0) ||
            (massprod3_->at(gp)[0] < 0.0) || (massprod4_->at(gp)[0] < 0.0))
        {
          cout << "1: " << massprod1_->at(gp)[0] << endl;
          cout << "2: " << massprod2_->at(gp)[0] << endl;
          cout << "3: " << massprod3_->at(gp)[0] << endl;
          cout << "4: " << massprod4_->at(gp)[0] << endl;
          dserror("negative mass production computed for at least one collagen fiber family!");
        }

//        LINALG::Matrix<NUM_STRESS_3D,1> Scomp(true);
//        cmatelastic.Scale(0.0);
//        EvaluateStress(glstrain, gp, &cmatelastic, &Scomp, firstiter, dt, time);
//        MassProduction(gp,C,&Cmatrix,*stress,&Smatrix,&massstress,&massprodcomp);
//
//        for (int id = 0; id < 6; id++) Residual(id) = (*stress)(id) - Scomp(id);
//        Residual(6) = massprod1_->at(gp)[0] - massprodcomp(0);
//        Residual(7) = massprod2_->at(gp)[0] - massprodcomp(1);
//        Residual(8) = massprod3_->at(gp)[0] - massprodcomp(2);
//        Residual(9) = massprod4_->at(gp)[0] - massprodcomp(3);

      } //while loop
      if (localistep == maxstep && Residual.Norm2() > params_->abstol_) dserror("local Newton iteration did not converge %e", Residual.Norm2());

      vismassstress_->at(gp) = massstress(1);

      //--------------------------------------------------------------------------------------
      // compute cmat
      // formulas are given in the corresponding pdf file, they are to long for here
      LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D> RHS(cmatelastic);
      LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D> LM(true);
      for (int id = 0; id < 6; id++) LM(id,id) = 1.0;
      LINALG::Matrix<3,3> SC(true);
      LINALG::Matrix<3,3> CSC(true);
      LINALG::Matrix<3,3> SCSC(true);
      SC.Multiply(Smatrix,Cmatrix);
      CSC.Multiply(Cmatrix,SC);
      SCSC.Multiply(Smatrix,CSC);

      // Fiber1
      double stretch = prestretchcollagen/collagenstretch1_->at(gp)[0];
      LINALG::Matrix<NUM_STRESS_3D,1> left(true);
      EvaluateFiber(glstrain, NULL, & left, a1_->at(gp), stretch);
      double facS = pow(stretch,2) * qdegrad / density * dt;
      left.Scale(facS);
      left.Update(-dt/density*qdegrad*params_->kappa_*J,Cinv,1.0);
      LINALG::Matrix<NUM_STRESS_3D,1> right(true);
      right(0) = a1_->at(gp)(0)*a1_->at(gp)(0); right(1) = a1_->at(gp)(1)*a1_->at(gp)(1);
      right(2) = a1_->at(gp)(2)*a1_->at(gp)(2); right(3) = a1_->at(gp)(0)*a1_->at(gp)(1);
      right(4) = a1_->at(gp)(1)*a1_->at(gp)(2); right(5) = a1_->at(gp)(0)*a1_->at(gp)(2);
      right.Update(1.0,Cinv,1.0/pow(collagenstretch1_->at(gp)[0],2));
      right.Scale(- 1.0 * massstress(0));
      LINALG::Matrix<3,1> SCa(true);
      SCa.Multiply(SC,a1_->at(gp));
      LINALG::Matrix<3,1> SCSCa(true);
      SCSCa.Multiply(SCSC,a1_->at(gp));
      double fac = 1.0 / massstress(0) / pow(collagenstretch1_->at(gp)[0],2) / pow(J,2);
      right(0) += fac * (2.0 * a1_->at(gp)(0) * SCSCa(0) + SCa(0) * SCa(0));
      right(1) += fac * (2.0 * a1_->at(gp)(1) * SCSCa(1) + SCa(1) * SCa(1));
      right(2) += fac * (2.0 * a1_->at(gp)(2) * SCSCa(2) + SCa(2) * SCa(2));
      right(3) += fac * (a1_->at(gp)(0) * SCSCa(1) + a1_->at(gp)(1) * SCSCa(0) + SCa(0) * SCa(1));
      right(4) += fac * (a1_->at(gp)(1) * SCSCa(2) + a1_->at(gp)(2) * SCSCa(1) + SCa(1) * SCa(2));
      right(5) += fac * (a1_->at(gp)(0) * SCSCa(2) + a1_->at(gp)(2) * SCSCa(0) + SCa(0) * SCa(2));
      right.Scale(0.5*massprodbasal_*params_->growthfactor_/homstressfixed);
      RHS.MultiplyNT(2.0,left,right,1.0);
      LINALG::Matrix<3,1> Ca(true);
      Ca.Multiply(Cmatrix,a1_->at(gp));
      LINALG::Matrix<3,1> CSCa(true);
      CSCa.Multiply(CSC,a1_->at(gp));
      fac = massprodbasal_ * params_->growthfactor_ / homstressfixed /pow(J,2)
          / massstress(0) / pow(collagenstretch1_->at(gp)[0],2);
      right(0) = fac * Ca(0) * CSCa(0);
      right(1) = fac * Ca(1) * CSCa(1);
      right(2) = fac * Ca(2) * CSCa(2);
      right(3) = fac * 0.5 * (Ca(0) * CSCa(1) + Ca(1) * CSCa(0));
      right(4) = fac * 0.5 * (Ca(1) * CSCa(2) + Ca(2) * CSCa(1));
      right(5) = fac * 0.5 * (Ca(0) * CSCa(2) + Ca(2) * CSCa(0));
      LM.MultiplyNT(-1.0,left,right,1.0);

      // Fiber2
      stretch = prestretchcollagen/collagenstretch2_->at(gp)[0];
      left.Scale(0.0);
      EvaluateFiber(glstrain, NULL, & left, a2_->at(gp), stretch);
      facS = pow(stretch,2) * qdegrad / density * dt;
      left.Scale(facS);
      left.Update(-dt/density*qdegrad*params_->kappa_*J,Cinv,1.0);
      right(0) = a2_->at(gp)(0)*a2_->at(gp)(0); right(1) = a2_->at(gp)(1)*a2_->at(gp)(1);
      right(2) = a2_->at(gp)(2)*a2_->at(gp)(2); right(3) = a2_->at(gp)(0)*a2_->at(gp)(1);
      right(4) = a2_->at(gp)(1)*a2_->at(gp)(2); right(5) = a2_->at(gp)(0)*a2_->at(gp)(2);
      right.Update(1.0,Cinv,1.0/pow(collagenstretch2_->at(gp)[0],2));
      right.Scale(- 1.0 * massstress(1));
      SCa.Scale(0.0);
      SCa.Multiply(SC,a2_->at(gp));
      SCSCa.Scale(0.0);
      SCSCa.Multiply(SCSC,a2_->at(gp));
      fac = 1.0 / massstress(1) / pow(collagenstretch2_->at(gp)[0],2) / pow(J,2);
      right(0) += fac * (2.0 * a2_->at(gp)(0) * SCSCa(0) + SCa(0) * SCa(0));
      right(1) += fac * (2.0 * a2_->at(gp)(1) * SCSCa(1) + SCa(1) * SCa(1));
      right(2) += fac * (2.0 * a2_->at(gp)(2) * SCSCa(2) + SCa(2) * SCa(2));
      right(3) += fac * (a2_->at(gp)(0) * SCSCa(1) + a2_->at(gp)(1) * SCSCa(0) + SCa(0) * SCa(1));
      right(4) += fac * (a2_->at(gp)(1) * SCSCa(2) + a2_->at(gp)(2) * SCSCa(1) + SCa(1) * SCa(2));
      right(5) += fac * (a2_->at(gp)(0) * SCSCa(2) + a2_->at(gp)(2) * SCSCa(0) + SCa(0) * SCa(2));
      right.Scale(0.5*massprodbasal_*params_->growthfactor_/homstressfixed);
      RHS.MultiplyNT(2.0,left,right,1.0);
      Ca.Scale(0.0);
      Ca.Multiply(Cmatrix,a2_->at(gp));
      CSCa.Scale(0.0);
      CSCa.Multiply(CSC,a2_->at(gp));
      fac = massprodbasal_ * params_->growthfactor_ / homstressfixed /pow(J,2)
          / massstress(1) / pow(collagenstretch2_->at(gp)[0],2);
      right(0) = fac * Ca(0) * CSCa(0);
      right(1) = fac * Ca(1) * CSCa(1);
      right(2) = fac * Ca(2) * CSCa(2);
      right(3) = fac * 0.5 * (Ca(0) * CSCa(1) + Ca(1) * CSCa(0));
      right(4) = fac * 0.5 * (Ca(1) * CSCa(2) + Ca(2) * CSCa(1));
      right(5) = fac * 0.5 * (Ca(0) * CSCa(2) + Ca(2) * CSCa(0));
      LM.MultiplyNT(-1.0,left,right,1.0);

      // Fiber3
      stretch = prestretchcollagen/collagenstretch3_->at(gp)[0];
      left.Scale(0.0);
      EvaluateFiber(glstrain, NULL, & left, a3_->at(gp), stretch);
      facS = pow(stretch,2) * qdegrad / density * dt;
      left.Scale(facS);
      left.Update(-dt/density*qdegrad*params_->kappa_*J,Cinv,1.0);
      right(0) = a3_->at(gp)(0)*a3_->at(gp)(0); right(1) = a3_->at(gp)(1)*a3_->at(gp)(1);
      right(2) = a3_->at(gp)(2)*a3_->at(gp)(2); right(3) = a3_->at(gp)(0)*a3_->at(gp)(1);
      right(4) = a3_->at(gp)(1)*a3_->at(gp)(2); right(5) = a3_->at(gp)(0)*a3_->at(gp)(2);
      right.Update(1.0,Cinv,1.0/pow(collagenstretch3_->at(gp)[0],2));
      right.Scale(- 1.0 * massstress(2));
      SCa.Scale(0.0);
      SCa.Multiply(SC,a3_->at(gp));
      SCSCa.Scale(0.0);
      SCSCa.Multiply(SCSC,a3_->at(gp));
      fac = 1.0 / massstress(2) / pow(collagenstretch3_->at(gp)[0],2) / pow(J,2);
      right(0) += fac * (2.0 * a3_->at(gp)(0) * SCSCa(0) + SCa(0) * SCa(0));
      right(1) += fac * (2.0 * a3_->at(gp)(1) * SCSCa(1) + SCa(1) * SCa(1));
      right(2) += fac * (2.0 * a3_->at(gp)(2) * SCSCa(2) + SCa(2) * SCa(2));
      right(3) += fac * (a3_->at(gp)(0) * SCSCa(1) + a3_->at(gp)(1) * SCSCa(0) + SCa(0) * SCa(1));
      right(4) += fac * (a3_->at(gp)(1) * SCSCa(2) + a3_->at(gp)(2) * SCSCa(1) + SCa(1) * SCa(2));
      right(5) += fac * (a3_->at(gp)(0) * SCSCa(2) + a3_->at(gp)(2) * SCSCa(0) + SCa(0) * SCa(2));
      right.Scale(0.5*massprodbasal_*params_->growthfactor_/homstressfixed);
      RHS.MultiplyNT(2.0,left,right,1.0);
      Ca.Scale(0.0);
      Ca.Multiply(Cmatrix,a3_->at(gp));
      CSCa.Scale(0.0);
      CSCa.Multiply(CSC,a3_->at(gp));
      fac = massprodbasal_ * params_->growthfactor_ / homstressfixed /pow(J,2)
          / massstress(2) / pow(collagenstretch3_->at(gp)[0],2);
      right(0) = fac * Ca(0) * CSCa(0);
      right(1) = fac * Ca(1) * CSCa(1);
      right(2) = fac * Ca(2) * CSCa(2);
      right(3) = fac * 0.5 * (Ca(0) * CSCa(1) + Ca(1) * CSCa(0));
      right(4) = fac * 0.5 * (Ca(1) * CSCa(2) + Ca(2) * CSCa(1));
      right(5) = fac * 0.5 * (Ca(0) * CSCa(2) + Ca(2) * CSCa(0));
      LM.MultiplyNT(-1.0,left,right,1.0);

      // Fiber4
      stretch = prestretchcollagen/collagenstretch4_->at(gp)[0];
      left.Scale(0.0);
      EvaluateFiber(glstrain, NULL, & left, a4_->at(gp), stretch);
      facS = pow(stretch,2) * qdegrad / density * dt;
      left.Scale(facS);
      left.Update(-dt/density*qdegrad*params_->kappa_*J,Cinv,1.0);
      right(0) = a4_->at(gp)(0)*a4_->at(gp)(0); right(1) = a4_->at(gp)(1)*a4_->at(gp)(1);
      right(2) = a4_->at(gp)(2)*a4_->at(gp)(2); right(3) = a4_->at(gp)(0)*a4_->at(gp)(1);
      right(4) = a4_->at(gp)(1)*a4_->at(gp)(2); right(5) = a4_->at(gp)(0)*a4_->at(gp)(2);
      right.Update(1.0,Cinv,1.0/pow(collagenstretch4_->at(gp)[0],2));
      right.Scale(- 1.0 * massstress(3));
      SCa.Scale(0.0);
      SCa.Multiply(SC,a4_->at(gp));
      SCSCa.Scale(0.0);
      SCSCa.Multiply(SCSC,a4_->at(gp));
      fac = 1.0 / massstress(3) / pow(collagenstretch4_->at(gp)[0],2) / pow(J,2);
      right(0) += fac * (2.0 * a4_->at(gp)(0) * SCSCa(0) + SCa(0) * SCa(0));
      right(1) += fac * (2.0 * a4_->at(gp)(1) * SCSCa(1) + SCa(1) * SCa(1));
      right(2) += fac * (2.0 * a4_->at(gp)(2) * SCSCa(2) + SCa(2) * SCa(2));
      right(3) += fac * (a4_->at(gp)(0) * SCSCa(1) + a4_->at(gp)(1) * SCSCa(0) + SCa(0) * SCa(1));
      right(4) += fac * (a4_->at(gp)(1) * SCSCa(2) + a4_->at(gp)(2) * SCSCa(1) + SCa(1) * SCa(2));
      right(5) += fac * (a4_->at(gp)(0) * SCSCa(2) + a4_->at(gp)(2) * SCSCa(0) + SCa(0) * SCa(2));
      right.Scale(0.5*massprodbasal_*params_->growthfactor_/homstressfixed);
      RHS.MultiplyNT(2.0,left,right,1.0);
      Ca.Scale(0.0);
      Ca.Multiply(Cmatrix,a4_->at(gp));
      CSCa.Scale(0.0);
      CSCa.Multiply(CSC,a4_->at(gp));
      fac = massprodbasal_ * params_->growthfactor_ / homstressfixed /pow(J,2)
          / massstress(3) / pow(collagenstretch4_->at(gp)[0],2);
      right(0) = fac * Ca(0) * CSCa(0);
      right(1) = fac * Ca(1) * CSCa(1);
      right(2) = fac * Ca(2) * CSCa(2);
      right(3) = fac * 0.5 * (Ca(0) * CSCa(1) + Ca(1) * CSCa(0));
      right(4) = fac * 0.5 * (Ca(1) * CSCa(2) + Ca(2) * CSCa(1));
      right(5) = fac * 0.5 * (Ca(0) * CSCa(2) + Ca(2) * CSCa(0));
      LM.MultiplyNT(-1.0,left,right,1.0);

      //----------------------------------------------------
      // solve linear system of equations: A.X=B
      //----------------------------------------------------
      LINALG::FixedSizeSerialDenseSolver<6,6,6> solver;
      solver.SetMatrix(LM);              // set A=LM
      solver.SetVectors(*cmat, RHS);           // set X=increment, B=RHS
      solver.FactorWithEquilibration(true); // "some easy type of preconditioning" (Michael)
      int err2 = solver.Factor();           // ?
      int err = solver.Solve();             // X = A^-1 B
      if ((err!=0) || (err2!=0))
        dserror("solving linear system for cmat failed");

    } // implicit

  } else {
    // Visualization of massstresss also for time < starttime
    LINALG::Matrix<3,3> temp1(true);
    LINALG::Matrix<3,3> temp2(true);
    LINALG::Matrix<4,1> temp3(true);
    LINALG::Matrix<4,1> massstress(true);
    MassProduction(gp,C,&temp1,(*stress),&temp2,&massstress,&temp3);
    vismassstress_->at(gp) = massstress(1);
  }// time > starttime
} // Evaluate

/*----------------------------------------------------------------------*
 |  EvaluateStress                                (private)        01/11|
 *----------------------------------------------------------------------*/
void MAT::ConstraintMixture::EvaluateStress
(
  const LINALG::Matrix<NUM_STRESS_3D,1>* glstrain,
  const int gp,
  LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D>* cmat,
  LINALG::Matrix<NUM_STRESS_3D,1>* stress,
  const int firstiter,
  double dt,
  double time
)
{
  //--------------------------------------------------------------------------------------
  // some variables
  int numpast = numpast_;
  double prestretchelastin = params_->prestretchelastin_;
  double prestretchcollagen = params_->prestretchcollagen_;
  double density = params_->density_;
  double currmassdens = 0.0;

  //--------------------------------------------------------------------------------------
  // calculate stress and elasticity matrix

  // 1st step: elastin
  //==========================
  double refmassdenselastin = params_->phielastin_ * density;
  if (refmassdenselastin < 0.0 || refmassdenselastin > density) dserror("mass fraction of elastin not in [0;1]");
  // account for isotropic prestretch of elastin
  LINALG::Matrix<NUM_STRESS_3D,1> glstrainiso(*glstrain);
  glstrainiso.Scale(pow(prestretchelastin,2));
  LINALG::Matrix<NUM_STRESS_3D,1> Siso(true);
  LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D> cmatiso(true);
  EvaluateElastin(& glstrainiso, & cmatiso, & Siso);
  Siso.Scale(refmassdenselastin/density*pow(prestretchelastin,2));
  (*stress) = Siso;
  cmatiso.Scale(refmassdenselastin/density*pow(prestretchelastin,4));
  (*cmat) = cmatiso;
  currmassdens += refmassdenselastin;

  // 2nd step: collagen1
  //==========================
  for (int idpast = firstiter; idpast < numpast; idpast++)
  {
    double stretch = prestretchcollagen/collagenstretch1_->at(gp)[idpast];
    LINALG::Matrix<NUM_STRESS_3D,1> Saniso1(true);
    LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D> cmataniso1(true);
    EvaluateFiber(glstrain, & cmataniso1, & Saniso1, a1_->at(gp), stretch);
    double qdegrad = 0.0;
    Degradation(dt*idpast, &qdegrad);
    double facS = pow(stretch,2) * qdegrad * massprod1_->at(gp)[idpast] / density * dt;
    Saniso1.Scale(facS);
    (*stress) += Saniso1;
    double faccmat = pow(stretch,4) * qdegrad * massprod1_->at(gp)[idpast] / density * dt;
    cmataniso1.Scale(faccmat);
    (*cmat) += cmataniso1;
    currmassdens += qdegrad * massprod1_->at(gp)[idpast] * dt;
  }

  // 3rd step: collagen2
  //==========================
  for (int idpast = firstiter; idpast < numpast; idpast++)
  {
    double stretch = prestretchcollagen/collagenstretch2_->at(gp)[idpast];
    LINALG::Matrix<NUM_STRESS_3D,1> Saniso2(true);
    LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D> cmataniso2(true);
    EvaluateFiber(glstrain, & cmataniso2, & Saniso2, a2_->at(gp), stretch);
    double qdegrad = 0.0;
    Degradation(dt*idpast, &qdegrad);
    double facS = pow(stretch,2) * qdegrad * massprod2_->at(gp)[idpast] / density * dt;
    Saniso2.Scale(facS);
    (*stress) += Saniso2;
    double faccmat = pow(stretch,4) * qdegrad * massprod2_->at(gp)[idpast] / density * dt;
    cmataniso2.Scale(faccmat);
    (*cmat) += cmataniso2;
    currmassdens += qdegrad * massprod2_->at(gp)[idpast] * dt;
  }

  // 4th step: collagen3
  //==========================
  for (int idpast = firstiter; idpast < numpast; idpast++)
  {
    double stretch = prestretchcollagen/collagenstretch3_->at(gp)[idpast];
    LINALG::Matrix<NUM_STRESS_3D,1> Saniso3(true);
    LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D> cmataniso3(true);
    EvaluateFiber(glstrain, & cmataniso3, & Saniso3, a3_->at(gp), stretch);
    double qdegrad = 0.0;
    Degradation(dt*idpast, &qdegrad);
    double facS = pow(stretch,2) * qdegrad * massprod3_->at(gp)[idpast] / density * dt;
    Saniso3.Scale(facS);
    (*stress) += Saniso3;
    double faccmat = pow(stretch,4) * qdegrad * massprod3_->at(gp)[idpast] / density * dt;
    cmataniso3.Scale(faccmat);
    (*cmat) += cmataniso3;
    currmassdens += qdegrad * massprod3_->at(gp)[idpast] * dt;
  }

  // 5th step: collagen4
  //==========================
  for (int idpast = firstiter; idpast < numpast; idpast++)
  {
    double stretch = prestretchcollagen/collagenstretch4_->at(gp)[idpast];
    LINALG::Matrix<NUM_STRESS_3D,1> Saniso4(true);
    LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D> cmataniso4(true);
    EvaluateFiber(glstrain, & cmataniso4, & Saniso4, a4_->at(gp), stretch);
    double qdegrad = 0.0;
    Degradation(dt*idpast, &qdegrad);
    double facS = pow(stretch,2) * qdegrad * massprod4_->at(gp)[idpast] / density * dt;
    Saniso4.Scale(facS);
    (*stress) += Saniso4;
    double faccmat = pow(stretch,4) * qdegrad * massprod4_->at(gp)[idpast] / density * dt;
    cmataniso4.Scale(faccmat);
    (*cmat) += cmataniso4;
    currmassdens += qdegrad * massprod4_->at(gp)[idpast] * dt;
  }

  // 6th step: volumetric part
  //==========================
  LINALG::Matrix<NUM_STRESS_3D,1> Svol(true);
  LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D> cmatvol(true);
  EvaluateVolumetric(glstrain, & cmatvol, & Svol, currmassdens, density);
  (*stress) += Svol;
  (*cmat) += cmatvol;
}
/*----------------------------------------------------------------------*
 |  EvaluateFiber                                 (private)        12/10|
 *----------------------------------------------------------------------*
 strain energy function

 W    = k1/(2.0*k2)*(exp(k2*(I_4 - 1.0)^2-1.0))

 I_4 .. invariant accounting for the fiber direction
 */
void MAT::ConstraintMixture::EvaluateFiber
(
  const LINALG::Matrix<NUM_STRESS_3D,1>* glstrain,
  LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D>* cmat,
  LINALG::Matrix<NUM_STRESS_3D,1>* stress,
  LINALG::Matrix<3,1> a,
  double stretch
)
{
  const double k1 = params_->k1_;
  const double k2 = params_->k2_;

  //--------------------------------------------------------------------------------------
  // build identity tensor I
  LINALG::Matrix<NUM_STRESS_3D,1> Id(true);
  for (int i = 0; i < 3; i++) Id(i) = 1.0;
  // right Cauchy-Green Tensor  C = 2 * E + I
  LINALG::Matrix<NUM_STRESS_3D,1> C(*glstrain);
  C.Scale(2.0);
  C += Id;

  //--------------------------------------------------------------------------------------
  // structural tensors in voigt notation
  // A = a x a
  LINALG::Matrix<NUM_STRESS_3D,1>  A;
  for (int i = 0; i < 3; i++) {
    A(i) = a(i)*a(i);
  }
  A(3) = a(0)*a(1); A(4) = a(1)*a(2); A(5) = a(0)*a(2);

  double I4 =  A(0)*C(0) + A(1)*C(1) + A(2)*C(2)
             + 1.*(A(3)*C(3) + A(4)*C(4) + A(5)*C(5)); //I4 = trace(A C)
  I4 = I4 * pow(stretch,2);  // account for prestretch and stretch at deposition time

  //--------------------------------------------------------------------------------------
  // fibers can only stretch/compress down to a minimal value
  // look this up, I am not sure. There is nothing stated explicitly.
  double fib1_tension = 1.0;
  if (I4 < 1.0)
  {
    I4 = 1.0;
    fib1_tension = 0.0;
  }

  //--- determine 2nd Piola Kirchhoff stresses S -----------------------------------------
  LINALG::Matrix<NUM_STRESS_3D,1> Saniso(A); // first compute S = 2 dW/dI4 A
  const double exp1 = exp(k2*(I4-1.)*(I4-1.));
  const double fib1 = 2.*(k1*(I4-1.)*exp1);  // 2 dW/dI4
  Saniso.Scale(fib1);  //S

  *stress = Saniso;

  //--- do elasticity matrix -------------------------------------------------------------
  if (cmat != NULL)
  {
    const double delta7 = fib1_tension* 4.0*(k1*exp1 + 2.0*k1*k2*(I4-1.0)*(I4-1.0)*exp1); // 4 d^2Wf/dI4dI4
    LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D> cmataniso; // isochoric elastic C from Fib
    for (int i = 0; i < 6; ++i) {
      for (int j = 0; j < 6; ++j) {
          cmataniso(i,j) = delta7 * A(i) * A(j);  // delta7 A x A
      }
    }
    *cmat = cmataniso;
  }
}

/*----------------------------------------------------------------------*
 |  EvaluateElastin                               (private)         12/10|
 *----------------------------------------------------------------------*
 strain energy function

 W    = 1/2 mue (I1-3)

*/
void MAT::ConstraintMixture::EvaluateElastin
(
  const LINALG::Matrix<NUM_STRESS_3D,1>* glstrain,
  LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D>* cmat,
  LINALG::Matrix<NUM_STRESS_3D,1>* stress
)
{
  double mue = params_->mue_;
  LINALG::Matrix<NUM_STRESS_3D,1> Siso(true);
  for (int i = 0; i < 3; i++) Siso(i) = mue;

  *stress = Siso;
}

/*----------------------------------------------------------------------*
 |  EvaluateVolumetric                            (private)        12/10|
 *----------------------------------------------------------------------*
 strain energy function

 W    = 1/2 kappa (J-M(t)/M(0))^2

*/
void MAT::ConstraintMixture::EvaluateVolumetric
(
  const LINALG::Matrix<NUM_STRESS_3D,1>* glstrain,
  LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D>* cmat,
  LINALG::Matrix<NUM_STRESS_3D,1>* stress,
  double currMassDens,
  double refMassDens
)
{
  double kappa = params_->kappa_;

  //--------------------------------------------------------------------------------------
  // build identity tensor I
  LINALG::Matrix<NUM_STRESS_3D,1> Id(true);
  for (int i = 0; i < 3; i++) Id(i) = 1.0;
  // right Cauchy-Green Tensor  C = 2 * E + I
  LINALG::Matrix<NUM_STRESS_3D,1> C(*glstrain);
  C.Scale(2.0);
  C += Id;

  // isotropic invariant
  const double I3 = C(0)*C(1)*C(2)
        + 0.25 * C(3)*C(4)*C(5)
        - 0.25 * C(1)*C(5)*C(5)
        - 0.25 * C(2)*C(3)*C(3)
        - 0.25 * C(0)*C(4)*C(4);    // 3rd invariant, determinant
  const double J = sqrt(I3);     // determinant of F
  const double p = kappa*(J- currMassDens/refMassDens); // dW_vol/dJ

  //--------------------------------------------------------------------------------------
  // invert C
  LINALG::Matrix<NUM_STRESS_3D,1> Cinv(6);
  Cinv(0) = C(1)*C(2) - 0.25*C(4)*C(4);
  Cinv(1) = C(0)*C(2) - 0.25*C(5)*C(5);
  Cinv(2) = C(0)*C(1) - 0.25*C(3)*C(3);
  Cinv(3) = 0.25*C(5)*C(4) - 0.5*C(3)*C(2);
  Cinv(4) = 0.25*C(3)*C(5) - 0.5*C(0)*C(4);
  Cinv(5) = 0.25*C(3)*C(4) - 0.5*C(5)*C(1);
  Cinv.Scale(1.0/I3);

  //--- determine 2nd Piola Kirchhoff stresses S -----------------------------------------
  LINALG::Matrix<NUM_STRESS_3D,1> Svol(true);
  // volumetric part J*kappa*(J-M(t)/M(0))*Cinv
  for (int i = 0; i < 6; i++) Svol(i) = J*p * Cinv(i);

  *stress = Svol;

  //--- do elasticity matrix -------------------------------------------------------------
  // cmatvol = J(p + J dp/dJ) Cinv x Cinv  -  2 J p Cinv o Cinv
  AddtoCmatHolzapfelProduct((*cmat),Cinv,(-2*J*p));
  for (int i = 0; i < 6; i++) {
    for (int j = 0; j < 6; j++) {
      (*cmat)(i,j) += J*(p+J*kappa) * Cinv(i) * Cinv(j);
    }
  }
}

/*----------------------------------------------------------------------*
 |  MassProduction                                (private)        01/11|
 *----------------------------------------------------------------------*
compute new deposition rates
I derived a strange form for the driving factor just with S and C:

m^k = mbasal * (1 + K * ( sqrt(a_i^T*C*S*C*S*C*a_i/detC) /lambda_k(t)/homeo - 1))

therefore we need C and S as matrices
*/
void MAT::ConstraintMixture::MassProduction
(
  const int gp,
  LINALG::Matrix<NUM_STRESS_3D,1> C,
  LINALG::Matrix<3,3>* Cmat,
  LINALG::Matrix<NUM_STRESS_3D,1> S,
  LINALG::Matrix<3,3>* Smat,
  LINALG::Matrix<4,1>* massstress,
  LINALG::Matrix<4,1>* massprodcomp
)
{
  LINALG::Matrix<3,3> Smatrix(true);
  Smatrix(0,0) = S(0);
  Smatrix(0,1) = S(3);
  Smatrix(0,2) = S(5);
  Smatrix(1,0) = Smatrix(0,1);
  Smatrix(1,1) = S(1);
  Smatrix(1,2) = S(4);
  Smatrix(2,0) = Smatrix(0,2);
  Smatrix(2,1) = Smatrix(1,2);
  Smatrix(2,2) = S(2);
  (*Smat) = Smatrix;
  LINALG::Matrix<3,3> Cmatrix(true);
  Cmatrix(0,0) = C(0);
  Cmatrix(0,1) = 0.5 * C(3);
  Cmatrix(0,2) = 0.5 * C(5);
  Cmatrix(1,0) = Cmatrix(0,1);
  Cmatrix(1,1) = C(1);
  Cmatrix(1,2) = 0.5 * C(4);
  Cmatrix(2,0) = Cmatrix(0,2);
  Cmatrix(2,1) = Cmatrix(1,2);
  Cmatrix(2,2) = C(2);
  (*Cmat) = Cmatrix;
  double detC = C(0)*C(1)*C(2) + 0.25 * C(3)*C(4)*C(5) - 0.25 * C(1)*C(5)*C(5)
        - 0.25 * C(2)*C(3)*C(3) - 0.25 * C(0)*C(4)*C(4);
  LINALG::Matrix<3,3> temp1(true);
  LINALG::Matrix<3,3> temp2(true);
  temp1.Multiply(Smatrix,Cmatrix);
  temp2.Multiply(Cmatrix,temp1);
  temp1.Scale(0.0);
  temp1.Multiply(Smatrix,temp2);
  temp2.Scale(0.0);
  temp2.Multiply(1.0/detC,Cmatrix,temp1);

  double massstress1 = 0.0;
  for (int i = 0; i < 3; i++){
    for (int j = 0; j < 3; j++){
      massstress1 += a1_->at(gp)(i) * temp2(i,j) * a1_->at(gp)(j);
    }
  }
  massstress1 = sqrt(massstress1) / collagenstretch1_->at(gp)[0];
  if (params_->homstress_ != 0.0)
    (*massprodcomp)(0) = massprodbasal_ * (1.0 + params_->growthfactor_ * (massstress1 / params_->homstress_ - 1.0));
  else
    (*massprodcomp)(0) = massprodbasal_ * (1.0 + params_->growthfactor_ * massstress1);
  //cout << "1: " << massstress1 << endl;

  double massstress2 = 0.0;
  for (int i=0; i < 3; i++){
    for (int j=0; j < 3; j++){
      massstress2 += a2_->at(gp)(i) * temp2(i,j) * a2_->at(gp)(j);
    }
  }
  massstress2 = sqrt(massstress2) / collagenstretch2_->at(gp)[0];
  if (params_->homstress_ != 0.0)
    (*massprodcomp)(1) = massprodbasal_ * (1.0 + params_->growthfactor_ * (massstress2 / params_->homstress_ - 1.0));
  else
    (*massprodcomp)(1) = massprodbasal_ * (1.0 + params_->growthfactor_ * massstress2);
  //cout << "2: " << massstress2 << endl;

  double massstress3 = 0.0;
  for (int i=0; i < 3; i++){
    for (int j=0; j < 3; j++){
      massstress3 += a3_->at(gp)(i) * temp2(i,j) * a3_->at(gp)(j);
    }
  }
  massstress3 = sqrt(massstress3) / collagenstretch3_->at(gp)[0];
  if (params_->homstress_ != 0.0)
    (*massprodcomp)(2) = massprodbasal_ * (1.0 + params_->growthfactor_ * (massstress3 / params_->homstress_ - 1.0));
  else
    (*massprodcomp)(2) = massprodbasal_ * (1.0 + params_->growthfactor_ * massstress3);
  //cout << "3: " << massstress3 << endl;

  double massstress4 = 0.0;
  for (int i=0; i < 3; i++){
    for (int j=0; j < 3; j++){
      massstress4 += a4_->at(gp)(i) * temp2(i,j) * a4_->at(gp)(j);
    }
  }
  massstress4 = sqrt(massstress4) / collagenstretch4_->at(gp)[0];
  if (params_->homstress_ != 0.0)
    (*massprodcomp)(3) = massprodbasal_ * (1.0 + params_->growthfactor_ * (massstress4 / params_->homstress_ - 1.0));
  else
    (*massprodcomp)(3) = massprodbasal_ * (1.0 + params_->growthfactor_ * massstress4);
  //cout << "4: " << massstress4 << endl;

  (*massstress)(0) = massstress1;
  (*massstress)(1) = massstress2;
  (*massstress)(2) = massstress3;
  (*massstress)(3) = massstress4;

  if (((*massprodcomp)(0) < 0) || ((*massprodcomp)(1) < 0) ||
      ((*massprodcomp)(2) < 0) || ((*massprodcomp)(3) < 0))
  {
    cout << "1: " << (*massprodcomp)(0) << endl;
    cout << "2: " << (*massprodcomp)(1) << endl;
    cout << "3: " << (*massprodcomp)(2) << endl;
    cout << "4: " << (*massprodcomp)(3) << endl;
    dserror("negative mass production computed for at least one collagen fiber family!");
  }
}

/*----------------------------------------------------------------------*
 |  EvaluateFiberVecs                             (public)         02/11|
 *----------------------------------------------------------------------*/
void MAT::ConstraintMixture::EvaluateFiberVecs
(const int gp, const LINALG::Matrix<3,3>& locsys, const LINALG::Matrix<3,3>& defgrd)
{
  // locsys holds the principal directions
  // The deformation gradient (defgrd) is needed in remodeling as then locsys is given in the
  // spatial configuration and thus the fiber vectors have to be pulled back in the reference
  // configuration as the material is evaluated there.
  // If this function is called during Setup defgrd should be replaced by the Identity.

  const double gamma = (45*PI)/180.; //angle for diagonal fibers
  LINALG::Matrix<3,1> ca1(true);
  LINALG::Matrix<3,1> ca2(true);
  LINALG::Matrix<3,1> ca3(true);
  LINALG::Matrix<3,1> ca4(true);

  for (int i = 0; i < 3; i++) {
    // store old fiber directions
//    olda1_->at(gp)[i] = a1_->at(gp)[i];
//    olda2_->at(gp)[i] = a2_->at(gp)[i];
//    olda3_->at(gp)[i] = a3_->at(gp)[i];
//    olda4_->at(gp)[i] = a4_->at(gp)[i];
    // a1 = e3, circumferential direction, used for collagen and smooth muscle
    ca1(i) = locsys(i,2);
    // a2 = e2
    ca2(i) = locsys(i,1);
    // a3 = cos gamma e3 + sin gamma e2
    ca3(i) = cos(gamma)*locsys(i,2) + sin(gamma)*locsys(i,1);
    // a4 = cos gamma e3 - sin gamma e2
    ca4(i) = cos(gamma)*locsys(i,2) - sin(gamma)*locsys(i,1);
  }

  // pull back in reference configuration
  LINALG::Matrix<3,1> a1_0(true);
  LINALG::Matrix<3,1> a2_0(true);
  LINALG::Matrix<3,1> a3_0(true);
  LINALG::Matrix<3,1> a4_0(true);
  LINALG::Matrix<3,3> idefgrd(false);
  idefgrd.Invert(defgrd);
  for (int i = 0; i < 3; i++) {
    a1_0(i) = idefgrd(i,0)*ca1(0) + idefgrd(i,1)*ca1(1) + idefgrd(i,2)*ca1(2);
    a2_0(i) = idefgrd(i,0)*ca2(0) + idefgrd(i,1)*ca2(1) + idefgrd(i,2)*ca2(2);
    a3_0(i) = idefgrd(i,0)*ca3(0) + idefgrd(i,1)*ca3(1) + idefgrd(i,2)*ca3(2);
    a4_0(i) = idefgrd(i,0)*ca4(0) + idefgrd(i,1)*ca4(1) + idefgrd(i,2)*ca4(2);
  }
  double a1_0norm = sqrt(a1_0(0)*a1_0(0) + a1_0(1)*a1_0(1) + a1_0(2)*a1_0(2));
  double a2_0norm = sqrt(a2_0(0)*a2_0(0) + a2_0(1)*a2_0(1) + a2_0(2)*a2_0(2));
  double a3_0norm = sqrt(a3_0(0)*a3_0(0) + a3_0(1)*a3_0(1) + a3_0(2)*a3_0(2));
  double a4_0norm = sqrt(a4_0(0)*a4_0(0) + a4_0(1)*a4_0(1) + a4_0(2)*a4_0(2));
  for (int i = 0; i < 3; i++) {
    a1_->at(gp)(i) = a1_0(i)/a1_0norm;
    a2_->at(gp)(i) = a2_0(i)/a2_0norm;
    a3_->at(gp)(i) = a3_0(i)/a3_0norm;
    a4_->at(gp)(i) = a4_0(i)/a4_0norm;
  }

  return;
}

/*----------------------------------------------------------------------*
 |  Debug output to gmsh-file                                      01/11|
 *----------------------------------------------------------------------*
 this needs to be copied to STR::TimInt::OutputStep() to enable debug output
 {
   discret_->SetState("displacement",Dis());
   MAT::ConstraintMixtureOutputToGmsh(discret_, GetStep(), 1);
 }
 just works with strtimint!
 don't forget to include constraintmixture.H */
void MAT::ConstraintMixtureOutputToGmsh
(
  const Teuchos::RCP<DRT::Discretization> dis,
  const int timestep,
  const int iter
)
{
  const std::string filebase = DRT::Problem::Instance()->OutputControlFile()->FileName();
  // file for stress
  std::stringstream filename;
  filename << filebase << "_visual" << std::setw(3) << setfill('0') << timestep << std::setw(2) << setfill('0') << iter << ".pos";
  std::ofstream f_system(filename.str().c_str());
  stringstream gmshfilecontent;
  gmshfilecontent << "View \" Time: " << timestep << " Iter: " << iter << " \" {" << endl;


  for (int iele=0; iele<dis->NumMyColElements(); ++iele)
  {
    const DRT::Element* actele = dis->lColElement(iele);

    // build current configuration
    vector<int> lm;
    vector<int> lmowner;
    vector<int> lmstride;
    actele->LocationVector(*dis,lm,lmowner,lmstride);
    Teuchos::RCP<const Epetra_Vector> disp = dis->GetState("displacement");
    vector<double> mydisp(lm.size(),0);

    Teuchos::RefCountPtr<MAT::Material> mat = actele->Material();
    MAT::ConstraintMixture* grow = static_cast <MAT::ConstraintMixture*>(mat.get());
    // Teuchos::RCP<vector<double> > mandel = grow->GetVis();

    // material plot at gauss points
    int ngp = 8; //very bad, does not work for other elements than HEX8!!!
    for (int gp = 0; gp < ngp; ++gp){
      vector<double> point = MAT::MatPointCoords(actele,mydisp,gp); //defined in contchainnetw

      // write mandel stress
      double mandelgp = grow->GetVis(gp);
      gmshfilecontent << "SP(" << scientific << point[0] << ",";
      gmshfilecontent << scientific << point[1] << ",";
      gmshfilecontent << scientific << point[2] << ")";
      gmshfilecontent << "{" << scientific
      << mandelgp
      << "};" << endl;

    }
  }

  gmshfilecontent << "};" << endl;
  f_system << gmshfilecontent.str();
  f_system.close();

  return;
}

#endif // CCADISCRET
