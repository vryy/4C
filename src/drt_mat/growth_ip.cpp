/*!----------------------------------------------------------------------
\file growth_ip.cpp
\brief
This file contains routines for an integration point based growth law
example input line
MAT 1 MAT_GROWTH DENS 1.0 IDMATELASTIC 2 STARTTIME 0.2 ENDTIME 100.0 KPLUS 0.5 MPLUS 4.0 KMINUS 0.25 MMINUS 5.0

Here a kinematic integration point based approach of growth is modeled.
For a detailed description see:
- Lubarda, V. & Hoger, A., On the mechanics of solids with a growing mass,
  International Journal of Solids and Structures, 2002, 39, 4627-4664
- Himpel, G.; Kuhl, E.; Menzel, A. & Steinmann, P., Computational modelling
  of isotropic multiplicative growth, Computer Modeling in Engineering
  and Sciences, 2005, 8, 119-134

<pre>
Maintainer: Susanna Tinkl
            tinkl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15265
</pre>
*----------------------------------------------------------------------*/

#ifdef CCADISCRET

#include <vector>
#include "growth_ip.H"
#include "logneohooke.H"
#include "holzapfelcardiovascular.H"
#include "humphreycardiovascular.H"
#include "aaaneohooke.H"
#include "neohooke.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_mat/matpar_bundle.H"
#include "../drt_lib/drt_utils.H"  // for function Factory in Unpack
#include "../drt_io/io_gmsh.H" // for debug plotting with gmsh
#include "../drt_io/io_control.H" // for debug plotting with gmsh
#include "contchainnetw.H" // for debug plotting with gmsh


/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
MAT::PAR::Growth::Growth(
  Teuchos::RCP<MAT::PAR::Material> matdata
  )
: Parameter(matdata),
  density_(matdata->GetDouble("DENS")),
  idmatelastic_(matdata->GetInt("IDMATELASTIC")),
  starttime_(matdata->GetDouble("STARTTIME")),
  endtime_(matdata->GetDouble("ENDTIME")),
  abstol_(matdata->GetDouble("TOL")),
  kthetaplus_(matdata->GetDouble("KPLUS")),
  mthetaplus_(matdata->GetDouble("MPLUS")),
  kthetaminus_(matdata->GetDouble("KMINUS")),
  mthetaminus_(matdata->GetDouble("MMINUS"))
{
}


Teuchos::RCP<MAT::Material> MAT::PAR::Growth::CreateMaterial()
{
  return Teuchos::rcp(new MAT::Growth(this));
}


MAT::GrowthType MAT::GrowthType::instance_;

DRT::ParObject* MAT::GrowthType::Create( const std::vector<char> & data )
{
  MAT::Growth* grow = new MAT::Growth();
  grow->Unpack(data);
  return grow;
}

/*----------------------------------------------------------------------*
 |  Constructor                                   (public)         02/10|
 *----------------------------------------------------------------------*/
MAT::Growth::Growth()
  : params_(NULL)
{
}


/*----------------------------------------------------------------------*
 |  Copy-Constructor                             (public)          02/10|
 *----------------------------------------------------------------------*/
MAT::Growth::Growth(MAT::PAR::Growth* params)
  : params_(params)
{
}


/*----------------------------------------------------------------------*
 |  Pack                                          (public)         02/10|
 *----------------------------------------------------------------------*/
void MAT::Growth::Pack(vector<char>& data) const
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
    numgp = theta_->size();   // size is number of gausspoints
  }
  AddtoPack(data,numgp);
  // Pack internal variables
  for (int gp = 0; gp < numgp; ++gp)
  {
    AddtoPack(data,theta_->at(gp));
    AddtoPack(data,thetaold_->at(gp));
    AddtoPack(data,mandel_->at(gp));
  }

  // Pack data of elastic material
  vector<char> dataelastic;
  if (matelastic_!=Teuchos::null) {
  matelastic_->Pack(dataelastic);
  }
  AddtoPack(data,dataelastic);

  return;
}


/*----------------------------------------------------------------------*
 |  Unpack                                        (public)         02/10|
 *----------------------------------------------------------------------*/
void MAT::Growth::Unpack(const vector<char>& data)
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
  // in post-process mode we do not have any instance of DRT::Problem
  if (DRT::Problem::NumInstances() > 0)
  {
    const int probinst = DRT::Problem::Instance()->Materials()->GetReadFromProblem();
    MAT::PAR::Parameter* mat = DRT::Problem::Instance(probinst)->Materials()->ParameterById(matid);
    if (mat->Type() == MaterialType())
      params_ = static_cast<MAT::PAR::Growth*>(mat);
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

  // unpack growth internal variables
  theta_ = rcp(new vector<double> (numgp));
  thetaold_ = rcp(new vector<double> (numgp));
  mandel_ = rcp(new vector<double> (numgp));
  for (int gp = 0; gp < numgp; ++gp) {
    double a;
    ExtractfromPack(position,data,a);
    theta_->at(gp) = a;
    ExtractfromPack(position,data,a);
    thetaold_->at(gp) = a;
    ExtractfromPack(position,data,a);
    mandel_->at(gp) = a;
  }

  // Unpack data of elastic material (these lines are copied from drt_element.cpp)
  vector<char> dataelastic;
  ExtractfromPack(position,data,dataelastic);
  if (dataelastic.size()>0) {
    DRT::ParObject* o = DRT::UTILS::Factory(dataelastic);  // Unpack is done here
    MAT::Material* matel = dynamic_cast<MAT::Material*>(o);
    if (matel==NULL)
      dserror("failed to unpack elastic material");
    matelastic_ = rcp(matel);
  } else matelastic_ = Teuchos::null;

  // alternative way to unpack, but not in postprocessing
  // if (params_!=NULL) {
  //   matel_ = MAT::Material::Factory(params_->matelastic_);
  //   matel_->Unpack(dataelastic);
  // }

  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d",data.size(),position);

  return;
}

/*----------------------------------------------------------------------*
 |  Setup                                         (public)         02/10|
 *----------------------------------------------------------------------*/
void MAT::Growth::Setup(const int numgp, DRT::INPUT::LineDefinition* linedef)
{
  theta_ = rcp(new vector<double> (numgp));
  thetaold_ = rcp(new vector<double> (numgp));
  mandel_ = rcp(new vector<double> (numgp));
  for (int j=0; j<numgp; ++j)
  {
    theta_->at(j) = 1.0;
    thetaold_->at(j) = 1.0;
    mandel_->at(j) = 0.0;
  }

  // Setup of elastic material
  matelastic_ = MAT::Material::Factory(params_->idmatelastic_);
  if (matelastic_->MaterialType() == INPAR::MAT::m_holzapfelcardiovascular) {
    MAT::HolzapfelCardio* holz = static_cast <MAT::HolzapfelCardio*>(matelastic_.get());
    holz->Setup(numgp, linedef);
  } else if (matelastic_->MaterialType() == INPAR::MAT::m_humphreycardiovascular) {
    MAT::HumphreyCardio* hum = static_cast <MAT::HumphreyCardio*>(matelastic_.get());
    hum->Setup(numgp, linedef);
  }

  isinit_ = true;
  return;
}

/*----------------------------------------------------------------------*
 |  Update internal growth variables              (public)         02/10|
 *----------------------------------------------------------------------*/
void MAT::Growth::Update()
{
  const int histsize = theta_->size();
  for (int i=0; i<histsize; i++)
  {
    thetaold_->at(i) = theta_->at(i);
  }

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
void MAT::Growth::Evaluate
(
  const LINALG::Matrix<NUM_STRESS_3D,1>* glstrain,
  const int gp,
  LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D> * cmat,
  LINALG::Matrix<NUM_STRESS_3D,1> * stress,
  double dt,
  double time,
  bool output
)
{
  double eps = 1.0e-12;
  double endtime = params_->endtime_;

  // when stress output is calculated the final parameters already exist
  // we should not do another local Newton iteration, which uses eventually a wrong thetaold
  if (output) time = endtime + dt;

  if (time > params_->starttime_ + eps && time <= endtime + eps) {
    LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D> cmatelastic(true);
    LINALG::Matrix<NUM_STRESS_3D,1> Sdach(true);
    const double thetaold = thetaold_->at(gp);
    //double theta = theta_->at(gp);
    double theta = thetaold;
    double mandelcrit = 0.0; //1.0E-6;
    double signmandel = 1.0; // adjusts sign of mandelcrit to sign of mandel

    // check wether starttime is divisible by dt, if not adapt dt in first growth step
    if (time < params_->starttime_ + dt - eps) dt = time - params_->starttime_;

    //--------------------------------------------------------------------------------------
    // build identity tensor I
    LINALG::Matrix<NUM_STRESS_3D,1> Id(true);
    for (int i = 0; i < 3; i++) Id(i) = 1.0;

    // right Cauchy-Green Tensor  C = 2 * E + I
    LINALG::Matrix<NUM_STRESS_3D,1> C(*glstrain);
    C.Scale(2.0);
    C += Id;

    // elastic right Cauchy-Green Tensor Cdach = F_g^-T C F_g^-1
    LINALG::Matrix<NUM_STRESS_3D,1> Cdach(C);
    Cdach.Scale(1.0/theta/theta);
    // elastic Green Lagrange strain
    LINALG::Matrix<NUM_STRESS_3D,1> glstraindach(Cdach);
    glstraindach -= Id;
    glstraindach.Scale(0.5);
    // elastic 2 PK stress and constitutive matrix
    EvaluateElastic(&glstraindach,gp,&cmatelastic,&Sdach,output);

    // trace of elastic Mandel stress Mdach = Cdach Sdach
    double mandel = Cdach(0)*Sdach(0) + Cdach(1)*Sdach(1) + Cdach(2)*Sdach(2) +
                    Cdach(3)*Sdach(3) + Cdach(4)*Sdach(4) + Cdach(5)*Sdach(5);
    if (signmandel*mandel < 0) signmandel = -1.0*signmandel;

    // Evaluate growth law
    double ktheta = 0.0;
    double dktheta = 0.0;
    EvaluateGrowthLaw(&ktheta, &dktheta, mandel, theta, mandelcrit);
    double residual = thetaold - theta + ktheta*(mandel-signmandel*mandelcrit)*dt;

    int localistep = 0;
    double thetaquer = 0.0;
    int maxstep = 30;
    double abstol = params_->abstol_;

    // local Newton iteration to obtain exact theta
    while (abs(residual) > abstol && localistep < maxstep)
    {
      localistep += 1;
      double temp = 0.0;
      for (int i=0; i<6; i++)
      {
        temp += Cdach(i)*(cmatelastic(i,0)*Cdach(0) + cmatelastic(i,1)*Cdach(1) + cmatelastic(i,2)*Cdach(2)
              + cmatelastic(i,3)*Cdach(3) + cmatelastic(i,4)*Cdach(4) + cmatelastic(i,5)*Cdach(5));
      }
      thetaquer = 1.0 - (dktheta*(mandel-signmandel*mandelcrit) - ktheta/theta*(2.0*mandel+temp))*dt;

      // damping strategy
      double omega = 2.0;
      double thetatemp = theta;
      double residualtemp = residual;
      double omegamin = 1.0/64.0;
      while (abs(residualtemp) > (1.0-0.5*omega)*abs(residual) && omega > omegamin)
      {
        // update of theta
        omega = omega/2.0;
        thetatemp = theta + omega*residual/thetaquer;

        // update elastic variables
        Cdach = C;
        Cdach.Scale(1.0/thetatemp/thetatemp);
        glstraindach = Cdach;
        glstraindach -= Id;
        glstraindach.Scale(0.5);
        cmatelastic.Scale(0.0);
        Sdach.Scale(0.0);
        EvaluateElastic(&glstraindach,gp,&cmatelastic,&Sdach,output);

        // trace of mandel stress
        mandel = Cdach(0)*Sdach(0) + Cdach(1)*Sdach(1) + Cdach(2)*Sdach(2) +
                 Cdach(3)*Sdach(3) + Cdach(4)*Sdach(4) + Cdach(5)*Sdach(5);
        if (signmandel*mandel < 0) signmandel = -1.0*signmandel;

        ktheta = 0.0;
        dktheta = 0.0;
        EvaluateGrowthLaw(&ktheta, &dktheta, mandel, thetatemp, mandelcrit);
        residualtemp = thetaold - thetatemp + ktheta*(mandel-signmandel*mandelcrit)*dt;
      } // end of damping loop
      residual = residualtemp;
      theta = thetatemp;
      if (omega <= omegamin && abs(residualtemp) > (1.0-0.5*omega)*abs(residual))
      {
        cout << gp << ": Theta " << thetatemp << " residual " << residualtemp << " stress " << mandel << endl;
        //dserror("no damping coefficient found");
      }

    } // end of local Newton iteration
    if (localistep == maxstep && abs(residual) > abstol) dserror("local Newton iteration did not converge %e %f %f %e", residual, thetaold, theta, mandel);

    double temp = 0.0;
    for (int i=0; i<6; i++)
    {
      temp += Cdach(i)*(cmatelastic(i,0)*Cdach(0) + cmatelastic(i,1)*Cdach(1) + cmatelastic(i,2)*Cdach(2)
            + cmatelastic(i,3)*Cdach(3) + cmatelastic(i,4)*Cdach(4) + cmatelastic(i,5)*Cdach(5));
    }
    thetaquer = 1.0 - (dktheta*(mandel-signmandel*mandelcrit) - ktheta/theta*(2.0*mandel+temp))*dt;

    // 2PK stress S = F_g^-1 Sdach F_g^-T
    LINALG::Matrix<NUM_STRESS_3D,1> S(Sdach);
    S.Scale(1.0/theta/theta);
    *stress = S;

    // constitutive matrix including growth
    cmatelastic.Scale(1.0/theta/theta/theta/theta);
    for (int i=0; i<6; i++) {
      for (int j=0; j<6; j++) {
        double cmatelasCi = cmatelastic(i,0)*C(0) + cmatelastic(i,1)*C(1) + cmatelastic(i,2)*C(2)
                            + cmatelastic(i,3)*C(3) + cmatelastic(i,4)*C(4) + cmatelastic(i,5)*C(5);
        double Ccmatelasj = cmatelastic(0,j)*C(0) + cmatelastic(1,j)*C(1) + cmatelastic(2,j)*C(2)
                            + cmatelastic(3,j)*C(3) + cmatelastic(4,j)*C(4) + cmatelastic(5,j)*C(5);
        (*cmat)(i,j) = cmatelastic(i,j) - 2.0/theta/thetaquer*ktheta*dt* (2.0*S(i)+cmatelasCi) *(S(j)+0.5*Ccmatelasj);
      }
    }

    // store theta
    theta_->at(gp) = theta;
    mandel_->at(gp) = mandel;

  } else if (time > endtime + eps) {  // turn off growth or calculate stresses for output
    LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D> cmatelastic(true);
    LINALG::Matrix<NUM_STRESS_3D,1> Sdach(true);
    double theta = theta_->at(gp);

    //--------------------------------------------------------------------------------------
    // build identity tensor I
    LINALG::Matrix<NUM_STRESS_3D,1> Id(true);
    for (int i = 0; i < 3; i++) Id(i) = 1.0;

    // right Cauchy-Green Tensor  C = 2 * E + I
    LINALG::Matrix<NUM_STRESS_3D,1> C(*glstrain);
    C.Scale(2.0);
    C += Id;

    // elastic right Cauchy-Green Tensor Cdach = F_g^-T C F_g^-1
    LINALG::Matrix<NUM_STRESS_3D,1> Cdach(C);
    Cdach.Scale(1.0/theta/theta);
    // elastic Green Lagrange strain
    LINALG::Matrix<NUM_STRESS_3D,1> glstraindach(Cdach);
    glstraindach -= Id;
    glstraindach.Scale(0.5);
    // elastic 2 PK stress and constitutive matrix
    EvaluateElastic(&glstraindach,gp,&cmatelastic,&Sdach,output);

    // 2PK stress S = F_g^-1 Sdach F_g^-T
    LINALG::Matrix<NUM_STRESS_3D,1> S(Sdach);
    S.Scale(1.0/theta/theta);
    *stress = S;

    // constitutive matrix including growth
    cmatelastic.Scale(1.0/theta/theta/theta/theta);
    *cmat = cmatelastic;

    // trace of elastic Mandel stress Mdach = Cdach Sdach
    double mandel = Cdach(0)*Sdach(0) + Cdach(1)*Sdach(1) + Cdach(2)*Sdach(2) +
                    Cdach(3)*Sdach(3) + Cdach(4)*Sdach(4) + Cdach(5)*Sdach(5);
    mandel_->at(gp) = mandel;

  } else {
    EvaluateElastic(glstrain,gp,cmat,stress,output);
    // build identity tensor I
    LINALG::Matrix<NUM_STRESS_3D,1> Id(true);
    for (int i = 0; i < 3; i++) Id(i) = 1.0;
    // right Cauchy-Green Tensor  C = 2 * E + I
    LINALG::Matrix<NUM_STRESS_3D,1> C(*glstrain);
    C.Scale(2.0);
    C += Id;
    LINALG::Matrix<NUM_STRESS_3D,1> S(true);
    S = *stress;
    mandel_->at(gp) = C(0)*S(0) + C(1)*S(1) + C(2)*S(2) + C(3)*S(3) + C(4)*S(4) + C(5)*S(5);
  }
}

/*----------------------------------------------------------------------*
 |  Evaluate elastic Material                     (private)        02/10|
 *----------------------------------------------------------------------*/
void MAT::Growth::EvaluateElastic
(
  const LINALG::Matrix<NUM_STRESS_3D,1>* glstrain,
  const int gp,
  LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D> * cmat,
  LINALG::Matrix<NUM_STRESS_3D,1> * stress,
  bool output
)
{
  if (matelastic_->MaterialType() == INPAR::MAT::m_logneohooke) {
    MAT::LogNeoHooke* log = static_cast <MAT::LogNeoHooke*>(matelastic_.get());
    log->Evaluate(*glstrain, *cmat, *stress);
  } else if (matelastic_->MaterialType() == INPAR::MAT::m_holzapfelcardiovascular) {
    MAT::HolzapfelCardio* holz = static_cast <MAT::HolzapfelCardio*>(matelastic_.get());
    holz->Evaluate(glstrain, gp, cmat, stress, output);
  } else if (matelastic_->MaterialType() == INPAR::MAT::m_humphreycardiovascular) {
    MAT::HumphreyCardio* hum = static_cast <MAT::HumphreyCardio*>(matelastic_.get());
    hum->Evaluate(glstrain, gp, cmat, stress, output);
  } else if (matelastic_->MaterialType() == INPAR::MAT::m_aaaneohooke){
    MAT::AAAneohooke* aaaneo = static_cast <MAT::AAAneohooke*>(matelastic_.get());
    aaaneo->Evaluate(*glstrain, *cmat, *stress);
  } else if (matelastic_->MaterialType() == INPAR::MAT::m_neohooke){
    MAT::NeoHooke* neo = static_cast <MAT::NeoHooke*>(matelastic_.get());
    neo->Evaluate(*glstrain, *cmat, *stress);
  } else dserror("material not implemented for growth");

}

/*----------------------------------------------------------------------*
 |  Evaluate growth law                           (private)        02/10|
 *----------------------------------------------------------------------*/
void MAT::Growth::EvaluateGrowthLaw
(
  double * ktheta,
  double * dktheta,
  double traceM,
  double theta,
  double mandelcrit
)
{
  // parameters
  double kthetaplus = params_->kthetaplus_; //1.0; 0.5;
  double mthetaplus = params_->mthetaplus_; //2.0; 4.0;
  double thetaplus = 1.5;
  double kthetaminus = params_->kthetaminus_; //2.0; 0.25;
  double mthetaminus = params_->mthetaminus_; //3.0; 5.0;
  double thetaminus = 0.5;
  // ktheta and dktheta should be zero!
  if (traceM > mandelcrit) {
    *ktheta=kthetaplus*pow((thetaplus-theta)/(thetaplus-1.0),mthetaplus);
    *dktheta=kthetaplus*pow((thetaplus-theta)/(thetaplus-1.0),mthetaplus-1.0)/(1.0-thetaplus);
  } else if (traceM < -1.0*mandelcrit) {
    *ktheta=kthetaminus*pow((theta-thetaminus)/(1.0-thetaminus),mthetaminus);
    *dktheta=kthetaminus*pow((theta-thetaminus)/(1.0-thetaminus),mthetaminus-1.0)/(1.0-thetaminus);
  }

}

/*----------------------------------------------------------------------*
 |  Debug output to gmsh-file                                      10/10|
 *----------------------------------------------------------------------*
 this needs to be copied to STR::TimInt::OutputStep() to enable debug output
 {
   discret_->SetState("displacement",Dis());
   MAT::GrowthOutputToGmsh(discret_, GetStep(), 1);
 }
 don't forget to include growth_ip.H */
void MAT::GrowthOutputToGmsh
(
  const Teuchos::RCP<DRT::Discretization> dis,
  const int timestep,
  const int iter
)
{
  const std::string filebase = DRT::Problem::Instance()->OutputControlFile()->FileName();
  // file for mandel stress
  std::stringstream filename_mandel;
  filename_mandel << filebase << "_mandel" << std::setw(3) << setfill('0') << timestep << std::setw(2) << setfill('0') << iter << ".pos";
  std::ofstream f_system_mandel(filename_mandel.str().c_str());
  stringstream gmshfilecontent_mandel;
  gmshfilecontent_mandel << "View \" Time: " << timestep << " Iter: " << iter << " \" {" << endl;

  // file for theta
  std::stringstream filename_theta;
  filename_theta << filebase << "_theta" << std::setw(3) << setfill('0') << timestep << std::setw(2) << setfill('0') << iter << ".pos";
  std::ofstream f_system_theta(filename_theta.str().c_str());
  stringstream gmshfilecontent_theta;
  gmshfilecontent_theta << "View \" Time: " << timestep << " Iter: " << iter << " \" {" << endl;

  for (int iele=0; iele<dis->NumMyColElements(); ++iele)
  {
    const DRT::Element* actele = dis->lColElement(iele);

    // build current configuration
    vector<int> lm;
    vector<int> lmowner;
    actele->LocationVector(*dis,lm,lmowner);
    RCP<const Epetra_Vector> disp = dis->GetState("displacement");
    vector<double> mydisp(lm.size(),0);

    RefCountPtr<MAT::Material> mat = actele->Material();
    MAT::Growth* grow = static_cast <MAT::Growth*>(mat.get());
    RCP<vector<double> > mandel = grow->Getmandel();
    RCP<vector<double> > theta = grow->Gettheta();

    // material plot at gauss points
    int ngp = grow->Gettheta()->size();
    for (int gp = 0; gp < ngp; ++gp){
      vector<double> point = MAT::MatPointCoords(actele,mydisp,gp); //defined in contchainnetw

      // write mandel stress
      double mandelgp = mandel->at(gp);
      gmshfilecontent_mandel << "SP(" << scientific << point[0] << ",";
      gmshfilecontent_mandel << scientific << point[1] << ",";
      gmshfilecontent_mandel << scientific << point[2] << ")";
      gmshfilecontent_mandel << "{" << scientific
      << mandelgp
      << "};" << endl;

      // write theta
      double thetagp = theta->at(gp);
	  gmshfilecontent_theta << "SP(" << scientific << point[0] << ",";
	  gmshfilecontent_theta << scientific << point[1] << ",";
	  gmshfilecontent_theta << scientific << point[2] << ")";
	  gmshfilecontent_theta << "{" << scientific
	  << thetagp
      << "};" << endl;
    }
  }
  gmshfilecontent_mandel << "};" << endl;
  f_system_mandel << gmshfilecontent_mandel.str();
  f_system_mandel.close();

  gmshfilecontent_theta << "};" << endl;
  f_system_theta << gmshfilecontent_theta.str();
  f_system_theta.close();

  return;
}

#endif // CCADISCRET
