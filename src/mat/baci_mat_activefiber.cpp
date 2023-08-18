/*----------------------------------------------------------------------*/
/*! \file
\brief

 This file contains routines for a local material law that contains active
 fiber formation and orientation for the modeling of living cells.

 example input line
 MAT 1 MAT_ACTIVEFIBER DENS 1.0 DECAY 720.0 IDMATPASSIVE 2 KFOR 10.0 KBACK 1.0 KVAR 10.0
SIGMAX 3.9E+03 EPSNULL 2.8E-04

 For a detailed description of the model see:

 - Deshpande, V., McMeeking, R.M., Evans, A.G., 2007. A model for the
   contractility of the cytoskeleton including the effects of stress-fibre
   formation and dissociation, Proceedings of the Royal Society A:
   Mathematical, Physical and Engineering Sciences 463, 787-815.

   ________________________________________________________________________________________
  | !!! ATTENTION !!! Many major mistakes in literature which are corrected in here !!!!   |
  |                                                                                        |
  | SMOOTHED implementation using sigmoid function for improved convergence and stability. |
  |________________________________________________________________________________________|



\level 3



*----------------------------------------------------------------------*/

#include "baci_mat_activefiber.H"

#include "baci_discretization_fem_general_utils_integration.H"
#include "baci_io_control.H"
#include "baci_lib_globalproblem.H"
#include "baci_linalg_four_tensor.H"
#include "baci_linalg_utils_densematrix_eigen.H"
#include "baci_linalg_utils_densematrix_svd.H"
#include "baci_mat_par_bundle.H"
#include "baci_mat_service.H"

/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
MAT::PAR::ActiveFiber::ActiveFiber(Teuchos::RCP<MAT::PAR::Material> matdata)
    : Parameter(matdata),
      density_(matdata->GetDouble("DENS")),
      decayconst_(matdata->GetDouble("DECAY")),
      idmatpassive_(matdata->GetInt("IDMATPASSIVE")),
      kforwards_(matdata->GetDouble("KFOR")),
      kbackwards_(matdata->GetDouble("KBACK")),
      kvariance_(matdata->GetDouble("KVAR")),
      sigmamax_(matdata->GetDouble("SIGMAX")),
      epsilonnull_(matdata->GetDouble("EPSNULL"))
{
  analyticalmaterialtangent_ =
      DRT::INPUT::IntegralValue<int>(
          DRT::Problem::Instance()->StructuralDynamicParams(), "MATERIALTANGENT") == 0;
}


Teuchos::RCP<MAT::Material> MAT::PAR::ActiveFiber::CreateMaterial()
{
  return Teuchos::rcp(new MAT::ActiveFiber(this));
}


/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
MAT::ActiveFiberType MAT::ActiveFiberType::instance_;


/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
DRT::ParObject* MAT::ActiveFiberType::Create(const std::vector<char>& data)
{
  auto* actfiber = new MAT::ActiveFiber();
  actfiber->Unpack(data);
  return actfiber;
}


/*----------------------------------------------------------------------*
 |  Constructor                                             rauch  07/14|
 *----------------------------------------------------------------------*/
MAT::ActiveFiber::ActiveFiber() : params_(nullptr) {}


/*----------------------------------------------------------------------*
 |  Copy-Constructor                                        rauch  07/14|
 *----------------------------------------------------------------------*/
MAT::ActiveFiber::ActiveFiber(MAT::PAR::ActiveFiber* params) : params_(params) {}


/*----------------------------------------------------------------------*
 |  Pack                                                    rauch  07/14|
 *----------------------------------------------------------------------*/
void MAT::ActiveFiber::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data, type);

  // analytical or fd material tangent
  // AddtoPack(data,analyticalmaterialtangent_);

  // matid
  int matid = -1;
  if (params_ != nullptr) matid = params_->Id();  // in case we are in post-process mode
  AddtoPack(data, matid);

  // pack history data
  int histsize;
  // if material is not initialized, i.e. start simulation, nothing to pack
  if (!Initialized())
  {
    histsize = 0;
  }
  else
  {
    // if material is initialized (restart): size equates number of gausspoints
    histsize = histdefgrdlast_->size();
  }

  AddtoPack(data, histsize);  // length of history vector(s)
  for (int var = 0; var < histsize; ++var)
  {
    // insert history vectors to AddtoPack
    AddtoPack(data, histdefgrdlast_->at(var));
    AddtoPack(data, etalast_->at(var));
    AddtoPack(data, sigmaomegaphilast_->at(var));
    AddtoPack(data, etahat_->at(var));
    AddtoPack(data, etahor_->at(var));
    AddtoPack(data, etaver_->at(var));
    AddtoPack(data, etadiag_->at(var));
    //    AddtoPack(data,dxx_->at(var));
    //    AddtoPack(data,dyy_->at(var));
    //    AddtoPack(data,dzz_->at(var));
    //    AddtoPack(data,dxy_->at(var));
    //    AddtoPack(data,dyz_->at(var));
    //    AddtoPack(data,dxz_->at(var));
  }

  // Pack data of passive elastic material
  if (matpassive_ != Teuchos::null)
  {
    matpassive_->Pack(data);
  }
}  // Pack()


/*----------------------------------------------------------------------*
 |  Unpack                                                  rauch  07/14|
 *----------------------------------------------------------------------*/
void MAT::ActiveFiber::Unpack(const std::vector<char>& data)
{
  isinit_ = true;
  std::vector<char>::size_type position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position, data, type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");

  // analyticalmaterialtangent_ = ExtractInt(position,data);

  // matid and recover params_
  int matid;
  ExtractfromPack(position, data, matid);
  params_ = nullptr;
  if (DRT::Problem::Instance()->Materials() != Teuchos::null)
  {
    if (DRT::Problem::Instance()->Materials()->Num() != 0)
    {
      const int probinst = DRT::Problem::Instance()->Materials()->GetReadFromProblem();
      MAT::PAR::Parameter* mat =
          DRT::Problem::Instance(probinst)->Materials()->ParameterById(matid);
      if (mat->Type() == MaterialType())
        params_ = static_cast<MAT::PAR::ActiveFiber*>(mat);
      else
        dserror("Type of parameter material %d does not fit to calling type %d", mat->Type(),
            MaterialType());
    }
  }

  // history data
  int histsize;
  ExtractfromPack(position, data, histsize);

  // if system is not yet initialized, the history vectors have to be initialized
  if (histsize == 0)
  {
    isinit_ = false;
    if (position != data.size())
      dserror("Mismatch in size of data %d <-> %d", data.size(), position);
    return;
  }

  // unpack deformation gradient matrices
  histdefgrdlast_ = Teuchos::rcp(new std::vector<CORE::LINALG::Matrix<3, 3>>);
  histdefgrdcurr_ = Teuchos::rcp(new std::vector<CORE::LINALG::Matrix<3, 3>>);

  // unpack formation level eta for n=240 different fiber orientations
  etalast_ = Teuchos::rcp(new std::vector<CORE::LINALG::Matrix<numbgp, twice>>);
  etacurr_ = Teuchos::rcp(new std::vector<CORE::LINALG::Matrix<numbgp, twice>>);

  // unpack fiber stress sigma for n=240 different fiber orientations
  sigmaomegaphilast_ = Teuchos::rcp(new std::vector<CORE::LINALG::Matrix<numbgp, twice>>);
  sigmaomegaphicurr_ = Teuchos::rcp(new std::vector<CORE::LINALG::Matrix<numbgp, twice>>);

  // unpack average intensity level at every point in the cytoplasm
  etahat_ = Teuchos::rcp(new std::vector<double>);
  etahor_ = Teuchos::rcp(new std::vector<double>);
  etaver_ = Teuchos::rcp(new std::vector<double>);
  etadiag_ = Teuchos::rcp(new std::vector<double>);

  //  dxx_  = Teuchos::rcp( new std::vector<double> );
  //  dyy_  = Teuchos::rcp( new std::vector<double> );
  //  dzz_  = Teuchos::rcp( new std::vector<double> );
  //  dxy_  = Teuchos::rcp( new std::vector<double> );
  //  dyz_  = Teuchos::rcp( new std::vector<double> );
  //  dxz_  = Teuchos::rcp( new std::vector<double> );

  for (int var = 0; var < histsize; ++var)
  {
    // initialize
    CORE::LINALG::Matrix<3, 3> tmp_matrix3x3(true);
    CORE::LINALG::Matrix<numbgp, twice> tmp_matrix(true);
    double tmp_scalar = 0.0;

    // matrices of last converged state are unpacked
    ExtractfromPack(position, data, tmp_matrix3x3);
    histdefgrdlast_->push_back(tmp_matrix3x3);

    // scalar-valued vector of last converged state are unpacked
    ExtractfromPack(position, data, tmp_matrix);
    etalast_->push_back(tmp_matrix);
    ExtractfromPack(position, data, tmp_matrix);
    sigmaomegaphilast_->push_back(tmp_matrix);

    ExtractfromPack(position, data, tmp_scalar);
    etahat_->push_back(tmp_scalar);
    ExtractfromPack(position, data, tmp_scalar);
    etahor_->push_back(tmp_scalar);
    ExtractfromPack(position, data, tmp_scalar);
    etaver_->push_back(tmp_scalar);
    ExtractfromPack(position, data, tmp_scalar);
    etadiag_->push_back(tmp_scalar);
    //    ExtractfromPack(position, data, tmp_scalar);
    //    dxx_ ->push_back(tmp_scalar);
    //    ExtractfromPack(position, data, tmp_scalar);
    //    dyy_ ->push_back(tmp_scalar);
    //    ExtractfromPack(position, data, tmp_scalar);
    //    dzz_ ->push_back(tmp_scalar);
    //    ExtractfromPack(position, data, tmp_scalar);
    //    dxy_ ->push_back(tmp_scalar);
    //    ExtractfromPack(position, data, tmp_scalar);
    //    dyz_ ->push_back(tmp_scalar);
    //    ExtractfromPack(position, data, tmp_scalar);
    //    dxz_ ->push_back(tmp_scalar);

    // current vectors have to be initialized
    histdefgrdcurr_->push_back(tmp_matrix3x3);
    etacurr_->push_back(tmp_matrix);
    sigmaomegaphicurr_->push_back(tmp_matrix);
  }

  // Unpack data of passive elastic material (these lines are copied from element.cpp)
  std::vector<char> dataelastic;
  ExtractfromPack(position, data, dataelastic);
  if (dataelastic.size() > 0)
  {
    DRT::ParObject* o = DRT::UTILS::Factory(dataelastic);  // Unpack is done here
    auto* matel = dynamic_cast<MAT::So3Material*>(o);
    if (matel == nullptr) dserror("failed to unpack passive material");
    matpassive_ = Teuchos::rcp(matel);
  }
  else
    matpassive_ = Teuchos::null;


  if (position != data.size()) dserror("Mismatch in size of data %d <-> %d", data.size(), position);

}  // Unpack()

/*----------------------------------------------------------------------*
 | initialize / allocate internal variables (public)        rauch  07/14|
 *----------------------------------------------------------------------*/
void MAT::ActiveFiber::Setup(int numgp, DRT::INPUT::LineDefinition* linedef)
{
  // initialize history variables
  histdefgrdlast_ = Teuchos::rcp(new std::vector<CORE::LINALG::Matrix<3, 3>>);
  histdefgrdcurr_ = Teuchos::rcp(new std::vector<CORE::LINALG::Matrix<3, 3>>);

  etalast_ = Teuchos::rcp(new std::vector<CORE::LINALG::Matrix<numbgp, twice>>);
  etacurr_ = Teuchos::rcp(new std::vector<CORE::LINALG::Matrix<numbgp, twice>>);

  sigmaomegaphilast_ = Teuchos::rcp(new std::vector<CORE::LINALG::Matrix<numbgp, twice>>);
  sigmaomegaphicurr_ = Teuchos::rcp(new std::vector<CORE::LINALG::Matrix<numbgp, twice>>);

  etahat_ = Teuchos::rcp(new std::vector<double>);
  etahor_ = Teuchos::rcp(new std::vector<double>);
  etaver_ = Teuchos::rcp(new std::vector<double>);
  etadiag_ = Teuchos::rcp(new std::vector<double>);

  //  dxx_  = Teuchos::rcp( new std::vector<double> );
  //  dyy_  = Teuchos::rcp( new std::vector<double> );
  //  dzz_  = Teuchos::rcp( new std::vector<double> );
  //  dxy_  = Teuchos::rcp( new std::vector<double> );
  //  dyz_  = Teuchos::rcp( new std::vector<double> );
  //  dxz_  = Teuchos::rcp( new std::vector<double> );

  // set all history variables to zero
  histdefgrdlast_->resize(numgp);
  histdefgrdcurr_->resize(numgp);

  etalast_->resize(numgp);
  etacurr_->resize(numgp);

  sigmaomegaphilast_->resize(numgp);
  sigmaomegaphicurr_->resize(numgp);

  etahat_->resize(numgp);
  etahor_->resize(numgp);
  etaver_->resize(numgp);
  etadiag_->resize(numgp);

  //  dxx_->resize(numgp);
  //  dyy_->resize(numgp);
  //  dzz_->resize(numgp);
  //  dxy_->resize(numgp);
  //  dyz_->resize(numgp);
  //  dxz_->resize(numgp);

  CORE::LINALG::Matrix<numbgp, twice> emptymat(true);

  CORE::LINALG::Matrix<3, 3> emptymat3x3(true);
  for (int i = 0; i < 3; i++) emptymat3x3(i, i) = 1.0;

  for (int i = 0; i < numgp; i++)
  {
    histdefgrdlast_->at(i) = emptymat3x3;
    histdefgrdcurr_->at(i) = emptymat3x3;

    etalast_->at(i) = emptymat;
    etacurr_->at(i) = emptymat;

    sigmaomegaphilast_->at(i) = emptymat;
    sigmaomegaphicurr_->at(i) = emptymat;

    etahat_->at(i) = 0.0;
    etahor_->at(i) = 0.0;
    etaver_->at(i) = 0.0;
    etadiag_->at(i) = 0.0;

    //    dxx_->at(i) = 0.0;
    //    dyy_->at(i) = 0.0;
    //    dzz_->at(i) = 0.0;
    //    dxy_->at(i) = 0.0;
    //    dyz_->at(i) = 0.0;
    //    dxz_->at(i) = 0.0;
  }

  // Setup of passive material
  matpassive_ =
      Teuchos::rcp_dynamic_cast<MAT::So3Material>(MAT::Material::Factory(params_->idmatpassive_));
  matpassive_->Setup(numgp, linedef);

  isinit_ = true;

}  // Setup()

/*----------------------------------------------------------------------*
 |  ResetAll                                                rauch  07/14|
 *----------------------------------------------------------------------*/
void MAT::ActiveFiber::ResetAll(const int numgp)
{
  // initialize history variables
  histdefgrdlast_ = Teuchos::rcp(new std::vector<CORE::LINALG::Matrix<3, 3>>);
  histdefgrdcurr_ = Teuchos::rcp(new std::vector<CORE::LINALG::Matrix<3, 3>>);

  etalast_ = Teuchos::rcp(new std::vector<CORE::LINALG::Matrix<numbgp, twice>>);
  etacurr_ = Teuchos::rcp(new std::vector<CORE::LINALG::Matrix<numbgp, twice>>);

  sigmaomegaphilast_ = Teuchos::rcp(new std::vector<CORE::LINALG::Matrix<numbgp, twice>>);
  sigmaomegaphicurr_ = Teuchos::rcp(new std::vector<CORE::LINALG::Matrix<numbgp, twice>>);

  etahat_ = Teuchos::rcp(new std::vector<double>);
  etahor_ = Teuchos::rcp(new std::vector<double>);
  etaver_ = Teuchos::rcp(new std::vector<double>);
  etadiag_ = Teuchos::rcp(new std::vector<double>);

  //  dxx_  = Teuchos::rcp( new std::vector<double> );
  //  dyy_  = Teuchos::rcp( new std::vector<double> );
  //  dzz_  = Teuchos::rcp( new std::vector<double> );
  //  dxy_  = Teuchos::rcp( new std::vector<double> );
  //  dyz_  = Teuchos::rcp( new std::vector<double> );
  //  dxz_  = Teuchos::rcp( new std::vector<double> );

  // set all history variables to zero
  histdefgrdlast_->resize(numgp);
  histdefgrdcurr_->resize(numgp);

  etalast_->resize(numgp);
  etacurr_->resize(numgp);

  sigmaomegaphilast_->resize(numgp);
  sigmaomegaphicurr_->resize(numgp);

  etahat_->resize(numgp);
  etahor_->resize(numgp);
  etaver_->resize(numgp);
  etadiag_->resize(numgp);

  CORE::LINALG::Matrix<numbgp, twice> emptymat(true);
  CORE::LINALG::Matrix<3, 3> emptymat3x3(true);
  for (int i = 0; i < 3; i++) emptymat3x3(i, i) = 1.0;

  for (int i = 0; i < numgp; i++)
  {
    histdefgrdlast_->at(i) = emptymat3x3;
    histdefgrdcurr_->at(i) = emptymat3x3;

    etalast_->at(i) = emptymat;
    etacurr_->at(i) = emptymat;

    sigmaomegaphilast_->at(i) = emptymat;
    sigmaomegaphicurr_->at(i) = emptymat;

    etahat_->at(i) = 0.0;
    etahor_->at(i) = 0.0;
    etaver_->at(i) = 0.0;
    etadiag_->at(i) = 0.0;

    //    dxx_->at(i) = 0.0;
    //    dyy_->at(i) = 0.0;
    //    dzz_->at(i) = 0.0;
    //    dxy_->at(i) = 0.0;
    //    dyz_->at(i) = 0.0;
    //    dxz_->at(i) = 0.0;
  }

  matpassive_->ResetAll(numgp);

  isinit_ = true;

}  // ResetAll()

/*----------------------------------------------------------------------*
 |  Update internal variables                               rauch  07/14|
 *----------------------------------------------------------------------*/
void MAT::ActiveFiber::Update()
{
  // make current values at time step t_last+1 to values of last step t_last
  histdefgrdlast_ = histdefgrdcurr_;

  etalast_ = etacurr_;
  sigmaomegaphilast_ = sigmaomegaphicurr_;

  // empty vectors of current data
  histdefgrdcurr_ = Teuchos::rcp(new std::vector<CORE::LINALG::Matrix<3, 3>>);

  etacurr_ = Teuchos::rcp(new std::vector<CORE::LINALG::Matrix<numbgp, twice>>);
  sigmaomegaphicurr_ = Teuchos::rcp(new std::vector<CORE::LINALG::Matrix<numbgp, twice>>);

  // get the size of the vector
  // (use the last vector, because it includes latest results, current is empty)
  const int histsize = histdefgrdlast_->size();

  histdefgrdcurr_->resize(histsize);
  etacurr_->resize(histsize);
  sigmaomegaphicurr_->resize(histsize);

  CORE::LINALG::Matrix<numbgp, twice> emptymat(true);
  CORE::LINALG::Matrix<3, 3> emptymat3x3(true);
  for (int i = 0; i < 3; i++) emptymat3x3(i, i) = 1.0;
  for (int i = 0; i < histsize; i++)
  {
    histdefgrdcurr_->at(i) = emptymat3x3;
    etacurr_->at(i) = emptymat;
    sigmaomegaphicurr_->at(i) = emptymat;
  }

  matpassive_->Update();

}  // Update()

/*----------------------------------------------------------------------*
 |  Reset internal variables                                rauch  07/14|
 *----------------------------------------------------------------------*/
void MAT::ActiveFiber::ResetStep() { matpassive_->ResetStep(); }  // ResetStep()


/*----------------------------------------------------------------------*
 |  Evaluate Material                                       rauch  07/14|
 *----------------------------------------------------------------------*
 The stress response is decomposed into a passive and an active part:
     \sigma = \sigma_{passive} + \sigma_{active}
 */
void MAT::ActiveFiber::Evaluate(const CORE::LINALG::Matrix<3, 3>* defgrd,
    const CORE::LINALG::Matrix<6, 1>* glstrain, Teuchos::ParameterList& params,
    CORE::LINALG::Matrix<6, 1>* stress, CORE::LINALG::Matrix<6, 6>* cmat, const int gp,
    const int eleGID)
{
  //
  //          C1111 C1122 C1133 C1123 C1113 C1112
  //          C2211 C2222 C2233 C2223 C2213 C2212
  //          C3311 C3322 C3333 C3323 C3313 C3312     d S      d S_ij
  // Cmat =   C2311 C2322 C2333 C2323 C2313 C2312 = 2 ---- = 2 -------
  //          C1311 C1322 C1333 C1323 C1313 C1312     d C      d C_kl
  //          C1211 C1222 C1233 C1223 C1213 C1212
  //
  //

  //  std::cout<<" \n EVALUATE ACTIVE STRESS LAW FOR ELE WITH GID: "<<eleGID<<"\n"<<std::endl;

  // Get time algorithmic parameters
  double dt = params.get<double>("delta time", -1.0);

#ifdef DEBUG
  if (dt == -1.0) dserror("no time step size provided in material");
#endif

  //******************
  // PASSIVE PART
  //******************
  // Initialize passive stress and elasticity tensor
  CORE::LINALG::Matrix<NUM_STRESS_3D, NUM_STRESS_3D> cmatpassive(true);
  CORE::LINALG::Matrix<NUM_STRESS_3D, 1> Spassive(true);

  if (cmat != nullptr)
  {
    // Evaluate passive PK2 stress Spassive and passive elasticity tensor cmatpassive
    matpassive_->Evaluate(defgrd, glstrain, params, &Spassive, &cmatpassive, gp, eleGID);
  }

  //******************
  // ACTIVE PART
  //******************
  // Parameters for active constitutive model
  // Decay constant of activation signal
  double decayconst = params_->decayconst_;
  // Non-dimensional parameter governing the rate of formation of stress fibers
  double kforwards = params_->kforwards_;
  // Non-dimensional parameter governing the rate of dissociation of stress fibers
  double kbackwards = params_->kbackwards_;
  // Non-dimensional fiber rate sensitivity
  double kvariance = params_->kvariance_;
  // Maximum tension exerted by stress fibers
  double sigmamax = params_->sigmamax_;
  // Reference strain rate of cross-bridge dynamics law
  double epsilonnull = params_->epsilonnull_;
  // Use of analytical or finite difference material tangent?
  bool analyticalmaterialtangent = params_->analyticalmaterialtangent_;

  // Setup inverse of deformation gradient
  CORE::LINALG::Matrix<3, 3> invdefgrd(*defgrd);
  invdefgrd.Invert();

  // Setup deformation gradient rate, rotation tensor, strain rate and rotation rate
  // \dot{F} = \frac {F^n - F^{n-1}} {\Delta t}
  CORE::LINALG::Matrix<3, 3> defgrdrate(true);
  // R = F * U^{-1}
  CORE::LINALG::Matrix<3, 3> R(true);
  // \dot{\epsilon} = d = 0.5 * (\dot{F}F^{-1} + (\dot{F}F^{-1})^{T}
  CORE::LINALG::Matrix<6, 1> strainrate(true);
  // \dot{R} = \frac {R^n - R^{n-1}} {\Delta t}
  CORE::LINALG::Matrix<3, 3> rotationrate(true);

  //  // e = F^{-T} * E * F^{-1}
  //  CORE::LINALG::Matrix<3,3> eastrain(true);
  //  GLtoEA(*glstrain,invdefgrd,eastrain);

  SetupRates(*defgrd, invdefgrd, params, defgrdrate, R, strainrate, rotationrate, gp, dt);

  // Evaluate active stress
  CORE::LINALG::Matrix<NUM_STRESS_3D, 1> sigma(true);  // 6x1
  CORE::LINALG::Matrix<3, 3> cauchystress(true);       // 3x3

  // Setup history variables
  CORE::LINALG::Matrix<numbgp, twice> emptymat(true);

  const CORE::LINALG::Matrix<numbgp, twice> etalast = etalast_->at(gp);
  CORE::LINALG::Matrix<numbgp, twice> etanew;

  double etahat = 0.0;
  double etahor = 0.0;
  double etaver = 0.0;
  double etadiag = 0.0;

  const CORE::LINALG::Matrix<numbgp, twice> sigmaomegaphilast = sigmaomegaphilast_->at(gp);
  CORE::LINALG::Matrix<numbgp, twice> sigmaomegaphinew(true);

  // Calculation of activation signal at current and last time step
  // C = \exp(frac{-t}{\theta})
  double Csignal = 0.0;
  double Csignalold = 0.0;
  CalcActivationSignal(&Csignal, params, &Csignalold);

  double theta = 1.0;

  // Setting up gauss quadrature (do not forget to adapt defines "numbgp" and "twice" in
  // activefiber.H)
  const CORE::DRT::UTILS::IntegrationPoints1D gausspoints(
      CORE::DRT::UTILS::GaussRule1D::line_10point);


  double tol = 1e-12;
  int maxiter = 6;

  // parameters for sigmoid function
  //
  // 1/(1+exp(-(ax+b)))
  //
  double onebyeps0 = 1.0 / epsilonnull;

  double s = 0.96;
  double onebys = 1.0 / s;
  double aa = 2.0 * log(s / (1.0 - s)) * kvariance * onebyeps0;
  double b = -1.0 * log(1.0 / s - 1.0);
  double exponent;
  double efunct;
  double denom;
  double onebydenom;

  // frequently used products and ratios
  double etafac = decayconst + theta * dt * (kforwards * Csignal + kbackwards);
  double onebyetafac = 1.0 / etafac;
  double auxfac = theta * dt * kforwards * Csignal;
  double auxfac2 = theta * dt * kbackwards;
  double onebysigmax = 1.0 / sigmamax;

  double deta_dsigma = onebyetafac * auxfac2 * onebysigmax * onebys;

  //////////////////////////////////////////////////
  // gauss quadrature loop over spherical surface //
  //////////////////////////////////////////////////
  for (int j = 0; j < 2 * gausspoints.nquad; j++)
  {
    for (int i = 0; i < gausspoints.nquad; i++)
    {
      // angles according to Atkinson_1982
      double omega = acos(gausspoints.qxg[i][0]);
      double phi = (((double)(j + 1)) * M_PI) / ((double)gausspoints.nquad);
      // sigma at \dot eps = 0 = s*eta*sigmax
      double eta_m =
          (1.0 - theta) * dt *
          ((1.0 - etalast(i, j)) * kforwards * Csignalold -
              kbackwards * (etalast(i, j) - (sigmaomegaphilast(i, j) * onebys * onebysigmax)));

      CORE::LINALG::Matrix<3, 1> m(true);
      m(0) = sin(omega) * cos(phi);
      m(1) = sin(omega) * sin(phi);
      m(2) = cos(omega);
      CORE::LINALG::Matrix<numbgp, twice> epsomegaphi(true);

      // Transform strain rate at each point to fiber strain rate in (omega,phi) direction
      // \dot{\epsilon} = \dot{\epsilon}_{ij} m_{i} m_{j}
      epsomegaphi(i, j) = strainrate(0) * m(0) * m(0) + strainrate(1) * m(1) * m(1) +
                          strainrate(2) * m(2) * m(2) + 2. * strainrate(3) * m(0) * m(1) +
                          2. * strainrate(4) * m(1) * m(2) + 2. * strainrate(5) * m(0) * m(2);

      if (abs(epsomegaphi(i, j)) < 1e-12) epsomegaphi(i, j) = 0.0;

      ///////////////////////////////////////////////
      // local newton to find sigma_(n+1)
      ///////////////////////////////////////////////
      double residual = 1.0;

      // predict (sigma at \dot eps = 0 = s*eta*sigmax)
      etanew(i, j) =
          onebyetafac * (etalast(i, j) * decayconst + auxfac +
                            auxfac2 * sigmaomegaphinew(i, j) * onebys * onebysigmax + eta_m);

      exponent = (1.0 / etanew(i, j)) * aa * epsomegaphi(i, j) + b;

      if (exponent > 31.0)
      {
        efunct = 0.0;
        denom = 1.0;
        onebydenom = 1.0;

        residual = sigmaomegaphinew(i, j) - (etanew(i, j) * sigmamax);
      }
      else if (exponent < -30.0)
      {
        efunct = 1e13;
        denom = 1e13;
        onebydenom = 0.0;

        sigmaomegaphinew(i, j) = 0.0;
        etanew(i, j) =
            onebyetafac * (etalast(i, j) * decayconst + auxfac +
                              auxfac2 * sigmaomegaphinew(i, j) * onebys * onebysigmax + eta_m);
        residual = 0.0;
      }
      else
      {
        efunct = exp(-exponent);
        denom = 1.0 + efunct;
        onebydenom = 1. / denom;

        residual = sigmaomegaphinew(i, j) - (etanew(i, j) * sigmamax) * onebydenom;
      }

      /////////////////
      // newton loop //
      /////////////////
      int newtonstep = 0;
      double inc = 0.0;
      while (abs(residual) > tol)
      {
        // do newton step sig_i+1 = sig_i - f/f'
        if (exponent < -30.0)
        {
          inc = residual;
          sigmaomegaphinew(i, j) = 0.0;
          etanew(i, j) =
              onebyetafac * (etalast(i, j) * decayconst + auxfac +
                                auxfac2 * sigmaomegaphinew(i, j) * onebys * onebysigmax + eta_m);
          residual = 0.0;
          break;
        }
        else if (exponent > 31.0)
          inc = residual / (1.0 - (sigmamax * deta_dsigma));
        else
        {
          inc = residual / (1.0 - ((denom * sigmamax) - (1.0 / etanew(i, j) * sigmamax * efunct *
                                                            aa * epsomegaphi(i, j))) *
                                      deta_dsigma * onebydenom * onebydenom);
        }

        // increment stress
        sigmaomegaphinew(i, j) = sigmaomegaphinew(i, j) - inc;

        // evaluate corresponding fiber activation
        etanew(i, j) =
            onebyetafac * (etalast(i, j) * decayconst + auxfac +
                              auxfac2 * sigmaomegaphinew(i, j) * onebys * onebysigmax + eta_m);

        // update
        exponent = 1.0 / etanew(i, j) * aa * epsomegaphi(i, j) + b;

        if (exponent > 31.0)
        {
          efunct = 0.0;
          denom = 1.0;
          onebydenom = 1.0;

          // evaluate residual
          residual = sigmaomegaphinew(i, j) - (etanew(i, j) * sigmamax);
        }
        else if (exponent < -30.0)
        {
          efunct = 1e13;
          denom = 1e13;
          onebydenom = 0.0;

          sigmaomegaphinew(i, j) = 0.0;
          etanew(i, j) =
              onebyetafac * (etalast(i, j) * decayconst + auxfac +
                                auxfac2 * sigmaomegaphinew(i, j) * onebys * onebysigmax + eta_m);

          // set residual
          residual = 0.0;
        }
        else
        {
          efunct = exp(-exponent);
          denom = 1.0 + efunct;
          onebydenom = 1. / denom;

          // evaluate residual
          residual = sigmaomegaphinew(i, j) - (etanew(i, j) * sigmamax) * onebydenom;
        }

        newtonstep++;

        if (newtonstep == maxiter)
        {
          std::cout << "velocity gradient: \n" << std::setprecision(16) << strainrate << std::endl;
          std::cout << "direction vector m: \n" << std::setprecision(16) << m << std::endl;
          printf(
              "local newton iteration did not converge after %d steps. gid= %d , gp= %d , i= %d , "
              "j= %d ,omega= %f ,phi= %f , Residual= %e ,efunct= %e ,eta= %e ,sig= %e ,eps= %e "
              ",detadsig= %e ,eta_m= %e ,etan= %e ,etafac= %e ,sigmax= %e ,aa= %f ,b = %f ,eps0= "
              "%e ,kv= %f,Csig= %f,theta= %f , dt= %f ,decay= %f ,s= %f ,exponent= %f",
              newtonstep, eleGID, gp, i, j, omega, phi, residual, efunct, etanew(i, j),
              sigmaomegaphinew(i, j), epsomegaphi(i, j), deta_dsigma, eta_m, etalast(i, j), etafac,
              sigmamax, aa, b, epsilonnull, kvariance, Csignal, theta, dt, decayconst, s, exponent);
          residual = 0.0;
        }

      }  // newton loop

      // i = phi index  j = theta index
      // etahat is eta of new timestep integrated over sphere
      etahat += etanew(i, j) * (gausspoints.qwgt[i] / (4.0 * (double)gausspoints.nquad));

      // average stress tensor in voigt notation
      sigma(0) += sigmaomegaphinew(i, j) * gausspoints.qwgt[i] * m(0) * m(0);  // sigma_11
      sigma(1) += sigmaomegaphinew(i, j) * gausspoints.qwgt[i] * m(1) * m(1);  // sigma_22
      sigma(2) += sigmaomegaphinew(i, j) * gausspoints.qwgt[i] * m(2) * m(2);  // sigma_33
      sigma(3) += sigmaomegaphinew(i, j) * gausspoints.qwgt[i] * m(0) * m(1);  // sigma_13
      sigma(4) += sigmaomegaphinew(i, j) * gausspoints.qwgt[i] * m(1) * m(2);  // sigma_23
      sigma(5) += sigmaomegaphinew(i, j) * gausspoints.qwgt[i] * m(0) * m(2);  // sigma_13

    }  // loop over i
  }    // loop over j

  // eta in horizontal direction (omega =81 (90) degree, phi = 360 degree)
  etahor = (etanew(5, 19) + etanew(4, 19)) * 0.5;
  // eta in vertical direction (omega =81 (90) degree, phi = 90 degree)
  etaver = (etanew(5, 4) + etanew(4, 4)) * 0.5;
  // eta in diagonal direction (omega =81 (90) degree, phi = 36 (45) degree)
  etadiag = (etanew(5, 1) + etanew(4, 2)) * 0.5;

  // proper scaling of sigma with respect to gauss quadrature rule
  sigma.Scale(3. / (8. * gausspoints.nquad));

  // Update history
  etacurr_->at(gp) = etanew;
  sigmaomegaphicurr_->at(gp) = sigmaomegaphinew;
  // Store etahat
  etahat_->at(gp) = etahat;
  etahor_->at(gp) = etahor;
  etaver_->at(gp) = etaver;
  etadiag_->at(gp) = etadiag;

  // Transform Cauchy stress to PK2 stress
  // S = J * F^{-1} sigma F^{-T}
  CORE::LINALG::Matrix<NUM_STRESS_3D, 1> Sactive(true);  // 6x1
  CauchytoPK2(Sactive, cauchystress, *defgrd, invdefgrd, sigma);

  // Stress including active and passive part
  if (params.get<int>("iostress") == 0)
  {
    stress->Update(1.0, Spassive, 0.0);
    stress->Update(1.0, Sactive, 1.0);
  }
  else
  {
    stress->Update(1.0, Sactive, 1.0);  // write only active stresses as output
  }

#ifndef MATERIALFDCHECK
  if (cmat != nullptr and analyticalmaterialtangent)
  {
    // Setup active elasticity tensor cmatactive
    CORE::LINALG::Matrix<NUM_STRESS_3D, NUM_STRESS_3D> cmatactive(true);
    SetupCmatActive(cmatactive, rotationrate, strainrate, *defgrd, defgrdrate, R, invdefgrd, etanew,
        sigmaomegaphinew, cauchystress, params, theta, Csignal);

    // constitutive matrix including active and passive part
    cmat->Update(1.0, cmatpassive, 0.0);
    cmat->Update(1.0, cmatactive, 1.0);
  }
#else
  if (cmat != nullptr)
  {
    // Setup active elasticity tensor cmatactive
    CORE::LINALG::Matrix<NUM_STRESS_3D, NUM_STRESS_3D> cmatactive(true);
    SetupCmatActive(cmatactive, rotationrate, strainrate, *defgrd, defgrdrate, R, invdefgrd, etanew,
        sigmaomegaphinew, cauchystress, params, theta, Csignal);

    // constitutive matrix including active and passive part
    cmat->Update(1.0, cmatpassive, 0.0);
    cmat->Update(1.0, cmatactive, 1.0);
  }
#endif
}

/*----------------------------------------------------------------------------------*
 | Calculation of deformation gradient rate, rotation tensor, strain rate and       |
 | rotation rate (finite difference scheme)                            rauch  07/14 |
 *----------------------------------------------------------------------------------*/
// see also: viscoelasthyper.cpp line 781 ff
void MAT::ActiveFiber::SetupRates(const CORE::LINALG::Matrix<3, 3>& defgrd,
    const CORE::LINALG::Matrix<3, 3>& invdefgrd, Teuchos::ParameterList& params,
    CORE::LINALG::Matrix<3, 3>& defgrdrate, CORE::LINALG::Matrix<3, 3>& R,
    CORE::LINALG::Matrix<6, 1>& strainrate, CORE::LINALG::Matrix<3, 3>& rotationrate, const int& gp,
    const double& dt)
{
  // Read history
  CORE::LINALG::Matrix<3, 3> defgrdlast = histdefgrdlast_->at(gp);

  //  CORE::LINALG::Matrix<3,3> scaleddefgrd(defgrd);
  //  CORE::LINALG::Matrix<3,3> scaleddefgrdinv(scaleddefgrd.Invert());
  //  scaleddefgrd.Scale(params_->epsilonnull_);

  // Rate of deformation gradient: \dot{F} = \frac {F^n - F^{n-1}} {\Delta t}
  defgrdrate.Update(1.0, defgrd, 0.0);
  defgrdrate.Update(-1.0, defgrdlast, 1.0);
  defgrdrate.Scale(1.0 / dt);

  // Calculate velocity gradient l = \dot{F}.F^{-1}
  CORE::LINALG::Matrix<3, 3> velgradient(true);
  velgradient.MultiplyNN(defgrdrate, invdefgrd);

  // Rate of strain/symmetric part of velocity gradient
  // d = 0.5 * (l + l^{T}) = 0.5 * (\dot{F}F^{-1} + (\dot{F}F^{-1})^{T})
  // Remark: strain-like 6-Voigt vector
  strainrate(0) = velgradient(0, 0) + velgradient(0, 0);
  strainrate(1) = velgradient(1, 1) + velgradient(1, 1);
  strainrate(2) = velgradient(2, 2) + velgradient(2, 2);
  strainrate(3) = velgradient(0, 1) + velgradient(1, 0);
  strainrate(4) = velgradient(1, 2) + velgradient(2, 1);
  strainrate(5) = velgradient(0, 2) + velgradient(2, 0);
  strainrate.Scale(0.5);

  //  dxx_->at(gp)=strainrate(0);
  //  dyy_->at(gp)=strainrate(1);
  //  dzz_->at(gp)=strainrate(2);
  //  dxy_->at(gp)=strainrate(3);
  //  dyz_->at(gp)=strainrate(4);
  //  dxz_->at(gp)=strainrate(5);

  // Rate of rotation tensor (!= skew part w of velocity gradient l, see Holzapfel S.99)
  // Determine rotation tensor R from F (F=R*U) -> polar decomposition of displacement based F
  CORE::LINALG::Matrix<3, 3> Q(true);
  CORE::LINALG::Matrix<3, 3> S(true);
  CORE::LINALG::Matrix<3, 3> VT(true);

  // Calculate rotcurr from defgrd
  CORE::LINALG::SVD<3, 3>(defgrd, Q, S,
      VT);  // Singular Value Decomposition,analogously to micromaterial_evaluate.cpp lines 81ff
  R.MultiplyNN(Q, VT);

  // Update history of deformation gradient
  histdefgrdcurr_->at(gp) = defgrd;

}  // Evaluate()


/*---------------------------------------------------------------------------------------------*
 | Calculation of activation signal of current and last time step                  rauch  07/14|
 *---------------------------------------------------------------------------------------------*/
void MAT::ActiveFiber::CalcActivationSignal(
    double* Csignal, Teuchos::ParameterList& params, double* Csignalold)
{
  // Get constitutive law specific parameters
  double decayconst = params_->decayconst_;  // 720.0;
  // Get time algorithmic parameters
  double dt = params.get<double>("delta time", -1.0);
  double time = params.get<double>("total time", -1.0);

  if (time == -1.0) dserror("time step or total time not available");

  // Calculate activation signal (0 <= signal <= 1)
  // current
  *Csignal = exp(-1. * time / decayconst);

  // from last time step
  if (exp(-1. * (time - dt) / decayconst) < 1.0)
    *Csignalold = exp(-1. * (time - dt) / decayconst);
  else
    *Csignalold = 1.0;

}  // CalcActivationSignal()

/*----------------------------------------------------------------------*
 | pull back of spatial stresses                           rauch  07/14 |
 *----------------------------------------------------------------------*/
void MAT::ActiveFiber::CauchytoPK2(CORE::LINALG::Matrix<6, 1>& Sactive,
    CORE::LINALG::Matrix<3, 3>& cauchystress, const CORE::LINALG::Matrix<3, 3>& defgrd,
    const CORE::LINALG::Matrix<3, 3>& invdefgrd, CORE::LINALG::Matrix<6, 1> sigma)
{
  // calculate the Jacobi-determinant
  const double detF = defgrd.Determinant();  // const???

  // Convert stress like 6x1-Voigt vector to 3x3 matrix
  cauchystress(0, 0) = sigma(0);
  cauchystress(0, 1) = sigma(3);
  cauchystress(0, 2) = sigma(5);
  cauchystress(1, 0) = cauchystress(0, 1);
  cauchystress(1, 1) = sigma(1);
  cauchystress(1, 2) = sigma(4);
  cauchystress(2, 0) = cauchystress(0, 2);
  cauchystress(2, 1) = cauchystress(1, 2);
  cauchystress(2, 2) = sigma(2);

  // S = J * F^{-1} * sigma * F^{-T}
  CORE::LINALG::Matrix<3, 3> temp(true);
  CORE::LINALG::Matrix<3, 3> S(true);
  temp.MultiplyNN(invdefgrd, cauchystress);
  S.MultiplyNT(temp, invdefgrd);
  S.Scale(detF);

#ifdef DEBUG
  if (abs(S(1, 2) - S(2, 1)) > 1e-13 or abs(S(0, 2) - S(2, 0)) > 1e-13 or
      abs(S(0, 1) - S(1, 0)) > 1e-13)
  {
    std::cout << S << std::endl;
    dserror("PK2 not symmetric!!");
  }
#endif

  // Sactive is stress like 6x1-Voigt vector
  Sactive(0) = S(0, 0);
  Sactive(1) = S(1, 1);
  Sactive(2) = S(2, 2);
  Sactive(3) = S(0, 1);
  Sactive(4) = S(1, 2);
  Sactive(5) = S(0, 2);

}  // PK2toCauchy()

///*----------------------------------------------------------------------*
// | pull back of spatial stresses                           rauch  07/14 |
// *----------------------------------------------------------------------*/
// void MAT::ActiveFiber::GLtoEA(
//  CORE::LINALG::Matrix<6,1> glstrain,
//  CORE::LINALG::Matrix<3,3> invdefgrd,
//  CORE::LINALG::Matrix<3,3>& eastrain)    // muss noch initialisiert werden
//{
//  // Convert strain like 6x1-Voigt vector to 3x3 matrix
//  CORE::LINALG::Matrix<3,3> greenlagrange(true);
//  greenlagrange(0,0) = glstrain(0);
//  greenlagrange(0,1) = 0.5*glstrain(3);            // 0.5???
//  greenlagrange(0,2) = 0.5*glstrain(5);            // 0.5???
//  greenlagrange(1,0) = greenlagrange(0,1);
//  greenlagrange(1,1) = glstrain(1);
//  greenlagrange(1,2) = 0.5*glstrain(4);            //0.5???
//  greenlagrange(2,0) = greenlagrange(0,2);
//  greenlagrange(2,1) = greenlagrange(1,2);
//  greenlagrange(2,2) = glstrain(2);
//
//
//  // e = F^{-T} * E * F^{-1}
//  CORE::LINALG::Matrix<3,3> temp(true);
//  temp.MultiplyTN(invdefgrd,greenlagrange);
//  eastrain.MultiplyNT(temp,invdefgrd);
//
//}  // GLtoEA

/*---------------------------------------------------------------*
 |  matrix root                                     rauch  07/14 |
 *---------------------------------------------------------------*/
void MAT::ActiveFiber::MatrixRoot3x3(CORE::LINALG::Matrix<3, 3>& MatrixInOut)
{
  double Norm = MatrixInOut.Norm2();
  // direct calculation for zero-matrix
  if (Norm == 0.)
  {
    MatrixInOut.Clear();
    return;
  }
  else
  {
    CORE::LINALG::Matrix<3, 3> EV(MatrixInOut);
    CORE::LINALG::Matrix<3, 3> EW;

    // MatixInOut = EV * EW * EV^{-1}
    CORE::LINALG::SYEV(EV, EW, EV);

    // sqrt(MatrixInOut) = EV * sqrt(EW) * EVT
    // loop over all eigenvalues
    for (int a = 0; a < 3; a++) EW(a, a) = sqrt(EW(a, a));

    MatrixInOut.Clear();

    // temp = sqrt(EW) * EVT
    CORE::LINALG::Matrix<3, 3> temp;
    temp.MultiplyNT(EW, EV);

    // sqrt(MatrixInOut) = EV * sqrt(EW) * EV^{-1} = EV * temp
    MatrixInOut.MultiplyNN(EV, temp);
  }

}  // MatrixRoot3x3()

/*-------------------------------------------------------------------------------------*
 |  matrix root derivative of a symmetric 3x3 matrix                       rauch  07/14|
 *-------------------------------------------------------------------------------------*/
void MAT::ActiveFiber::MatrixRootDerivativeSym3x3(
    const CORE::LINALG::Matrix<3, 3>& MatrixIn, CORE::LINALG::Matrix<6, 6>& MatrixRootDeriv)
{
  double Norm = MatrixIn.Norm2();

  CORE::LINALG::Matrix<6, 6> id4sharp(true);  // souza S.31 eq. (2.110)???
  for (int i = 0; i < 3; i++) id4sharp(i, i) = 1.0;
  for (int i = 3; i < 6; i++) id4sharp(i, i) = 0.5;

  // direct calculation for zero-matrix
  if (Norm == 0.)
  {
    MatrixRootDeriv = id4sharp;
    dserror("d sqrt(C)/ d C not defined for C==0");
    return;
  }

  else
  {
    double EWtolerance = 1.e-12;  // see souza S.737 Remark A.2

    CORE::LINALG::Matrix<3, 3> EV(MatrixIn);
    CORE::LINALG::Matrix<3, 3> EW;
    CORE::LINALG::SYEV(EV, EW, EV);

    MatrixRootDeriv.Clear();
    // souza eq. (A.52)
    // note: EW stored in ascending order

    //  d X^2 / d X  =  1/2 * (  delta_jk X_lj + delta_il X_kj
    //                         + delta_jl X_ik + delta_kj X_il )    souza eq. (A.46)
    //
    // y_i = sqrt(x_i)
    // dy_i / dx_j = delta_ij 1/(2*sqrt(x_i))

    CORE::LINALG::Matrix<3, 3> id2(true);
    for (int i = 0; i < 3; i++) id2(i, i) = 1.0;
    //  // --------------------------------- switch by number of equal eigenvalues

    if (abs(EW(0, 0) - EW(1, 1)) < EWtolerance &&
        abs(EW(1, 1) - EW(2, 2)) < EWtolerance)  // ------------------ x_a == x_b == x_c
    {
      // calculate derivative
      MatrixRootDeriv = id4sharp;
      MatrixRootDeriv.Scale(1.0 / (2.0 * sqrt(EW(0, 0))));
    }

    else if ((abs(EW(0, 0) - EW(1, 1)) < EWtolerance && abs(EW(1, 1) - EW(2, 2)) > EWtolerance) ||
             (abs(EW(0, 0) - EW(1, 1)) > EWtolerance &&
                 abs(EW(1, 1) - EW(2, 2)) <
                     EWtolerance))  // ---- x_a != x_b == x_c or x_a == x_b != x_c
    {
      // scalar factors
      double s1 = 0.0;
      double s2 = 0.0;
      double s3 = 0.0;
      double s4 = 0.0;
      double s5 = 0.0;
      double s6 = 0.0;

      int a = 0;
      int c = 0;

      // switch which two EW are equal
      if (abs(EW(0, 0) - EW(1, 1)) < EWtolerance &&
          abs(EW(1, 1) - EW(2, 2)) > EWtolerance)  // ----------------------- x_a == x_b != x_c
      {
        a = 2;
        c = 0;
      }
      else if (abs(EW(0, 0) - EW(1, 1)) > EWtolerance &&
               abs(EW(1, 1) - EW(2, 2)) < EWtolerance)  // ------------------ x_a != x_b == x_c
      {
        a = 0;
        c = 2;
      }
      else
        dserror("you should not end up here");

      // in souza eq. (A.53):
      s1 = (sqrt(EW(a, a)) - sqrt(EW(c, c))) / (pow(EW(a, a) - EW(c, c), 2.0)) -
           (1.0 / (2.0 * sqrt(EW(c, c)))) / (EW(a, a) - EW(c, c));
      s2 = 2.0 * EW(c, c) * (sqrt(EW(a, a)) - sqrt(EW(c, c))) / (pow(EW(a, a) - EW(c, c), 2.0)) -
           (EW(a, a) + EW(c, c)) / (EW(a, a) - EW(c, c)) * (1.0 / (2.0 * sqrt(EW(c, c))));
      s3 = 2.0 * (sqrt(EW(a, a)) - sqrt(EW(c, c))) / (pow(EW(a, a) - EW(c, c), 3.0)) -
           ((1.0 / (2.0 * sqrt(EW(a, a)))) + (1.0 / (2.0 * sqrt(EW(c, c))))) /
               (pow(EW(a, a) - EW(c, c), 2.0));
      s4 = EW(c, c) * s3;
      s5 = s4;
      s6 = EW(c, c) * EW(c, c) * s3;

      // calculate derivative
      // + s_1 (d X^2 / d X)
      MAT::AddToCmatDerivTensorSquare(MatrixRootDeriv, s1, MatrixIn, 1.);
      // - s_2 I_s
      MatrixRootDeriv.Update(-s2, id4sharp, 1.);
      // - s_3 (X \dyad X)
      MAT::ElastSymTensorMultiply(MatrixRootDeriv, -1. * s3, MatrixIn, MatrixIn, 1.);
      // + s_4 (X \dyad I)
      MAT::ElastSymTensorMultiply(MatrixRootDeriv, s4, MatrixIn, id2, 1.);
      // + s_5 (I \dyad X)
      MAT::ElastSymTensorMultiply(MatrixRootDeriv, s5, id2, MatrixIn, 1.);
      // - s_6 (I \dyad I)
      MAT::ElastSymTensorMultiply(MatrixRootDeriv, -1. * s6, id2, id2, 1.);
    }

    else if (abs(EW(0, 0) - EW(1, 1)) > EWtolerance &&
             abs(EW(1, 1) - EW(2, 2)) > EWtolerance)  // ----------------- x_a != x_b != x_c
    {
      for (int a = 0; a < 3; a++)  // loop over all eigenvalues
      {
        // a=0 || a=1 || a=2
        int b = (a + 1) % 3;  // b=1 || b=2 || b=0     even (cyclic) permutations of (a,b,c)
        int c = (a + 2) % 3;  // c=2 || c=0 || c=1

        CORE::LINALG::Matrix<3, 1> ea;
        CORE::LINALG::Matrix<3, 1> eb;
        CORE::LINALG::Matrix<3, 1> ec;
        for (int i = 0; i < 3; i++)
        {
          ea(i) = EV(i, a);
          eb(i) = EV(i, b);
          ec(i) = EV(i, c);
        }
        CORE::LINALG::Matrix<3, 3> Ea;
        Ea.MultiplyNT(ea, ea);  // souza S.26 eq. (2.63)
        CORE::LINALG::Matrix<3, 3> Eb;
        Eb.MultiplyNT(eb, eb);
        CORE::LINALG::Matrix<3, 3> Ec;
        Ec.MultiplyNT(ec, ec);

        double fac = sqrt(EW(a, a)) / ((EW(a, a) - EW(b, b)) * (EW(a, a) - EW(c, c)));

        // calculate derivative
        // + d X^2 / d X
        MAT::AddToCmatDerivTensorSquare(MatrixRootDeriv, fac, MatrixIn, 1.);
        // - (x_b + x_c) I_s
        MatrixRootDeriv.Update(-1. * (EW(b, b) + EW(c, c)) * fac, id4sharp, 1.);
        // - [(x_a - x_b) + (x_a - x_c)] (E_a \dyad E_a)
        MAT::ElastSymTensorMultiply(MatrixRootDeriv,
            -1. * fac * ((EW(a, a) - EW(b, b)) + (EW(a, a) - EW(c, c))), Ea, Ea, 1.);
        // - (x_b - x_c) (E_b \dyad E_b)
        MAT::ElastSymTensorMultiply(MatrixRootDeriv, -1. * fac * (EW(b, b) - EW(c, c)), Eb, Eb, 1.);
        // + (x_b - x_c) (E_c \dyad E_c)
        MAT::ElastSymTensorMultiply(MatrixRootDeriv, fac * (EW(b, b) - EW(c, c)), Ec, Ec, 1.);
        // dy / dx_a (E_a \dyad E_a)
        MAT::ElastSymTensorMultiply(MatrixRootDeriv, 1.0 / (2.0 * sqrt(EW(a, a))), Ea, Ea, 1.);
      }  // end loop over all eigenvalues
    }

    else
      dserror("you should not end up here.");
  }

}  // MatrixRootDerivativeSym3x3()

/*-------------------------------------------------------------------------------------------*
 |  Computes active elasticity tensor in 6x6-Voigt notation                      rauch  07/14|
 *-------------------------------------------------------------------------------------------*/
void MAT::ActiveFiber::SetupCmatActive(CORE::LINALG::Matrix<6, 6>& cmatactive,
    const CORE::LINALG::Matrix<3, 3>& rotationrate, CORE::LINALG::Matrix<6, 1> strainrate,
    const CORE::LINALG::Matrix<3, 3>& defgrd, const CORE::LINALG::Matrix<3, 3>& defgrdrate,
    const CORE::LINALG::Matrix<3, 3>& R, const CORE::LINALG::Matrix<3, 3>& invdefgrd,
    CORE::LINALG::Matrix<numbgp, twice> etanew,
    const CORE::LINALG::Matrix<numbgp, twice>& sigmaomegaphicurr,
    const CORE::LINALG::Matrix<3, 3>& cauchystress, Teuchos::ParameterList& params, double theta,
    double Csignal)
{
  // Parameters of constitutive law
  double kforwards = params_->kforwards_;
  double kbackwards = params_->kbackwards_;
  double kvariance = params_->kvariance_;
  double sigmamax = params_->sigmamax_;
  double epsilonnull = params_->epsilonnull_;
  double decayconst = params_->decayconst_;

  // Time algorithmic parameter
  double dt = params.get<double>("delta time", -1.0);

  // Setup integration rule
  const CORE::DRT::UTILS::IntegrationPoints1D gausspoints(
      CORE::DRT::UTILS::GaussRule1D::line_10point);

  // Jacobi Determinant
  const double detF = defgrd.Determinant();

  // Right-Cauchy-Green tensor(3x3): C = F^{T} * F
  CORE::LINALG::Matrix<3, 3> C(true);
  C.MultiplyTN(defgrd, defgrd);
  // Inverse of C: C^{-1}
  CORE::LINALG::Matrix<3, 3> Cinv(C);
  Cinv.Invert();
  // Root of C: \sqrt{C}
  CORE::LINALG::Matrix<3, 3> RootC(C);
  MatrixRoot3x3(RootC);

  // Inverse of sqrt(C): \sqrt(C)^{-1}
  CORE::LINALG::Matrix<3, 3> RootCInv(RootC);
  RootCInv.Invert();


  // Derivative of \sqrt{C} with respect to C: DerviC =  d sqrt(C) / d C
  // CORE::LINALG::Matrix<6,6> DerivC(true);// 6x6 Voigt matrix
  // MatrixRootDerivativeSym3x3(C,DerivC);
  CORE::LINALG::FourTensor<3> TensorDerivC(true);  // initialize 81 4-Tensor values
  // Setup4Tensor(TensorDerivC,DerivC);// 3x3x3x3 Tensor

  // FDCHECK for d sqrt{C}/ d C
  double delta = 1e-8;
  CORE::LINALG::Matrix<3, 3> Ccopy(C);
  CORE::LINALG::Matrix<3, 3> RootCcopy(true);
  CORE::LINALG::Matrix<3, 3> InvRootCcopy(true);

  CORE::LINALG::FourTensor<3> RootCInvDerivCRootCInv(true);
  for (int k = 0; k < 3; ++k)
  {
    for (int l = 0; l < 3; ++l)
    {
      // toggle C_ij
      Ccopy(k, l) += delta / 2;
      Ccopy(l, k) += delta / 2;
      // calc root of toggled C
      RootCcopy.Update(Ccopy);
      MatrixRoot3x3(RootCcopy);
      InvRootCcopy.Update(RootCcopy);
      InvRootCcopy.Invert();
      // finite difference
      RootCcopy.Update(-1.0, RootC, 1.0);
      InvRootCcopy.Update(-1.0, RootCInv, 1.0);
      RootCcopy.Scale(1e8);
      InvRootCcopy.Scale(1e8);
      // fill 4-Tensor
      for (int i = 0; i < 3; ++i)
      {
        for (int j = 0; j < 3; ++j)
        {
          TensorDerivC(i, j, k, l) = RootCcopy(i, j);
          RootCInvDerivCRootCInv(i, j, k, l) = InvRootCcopy(i, j);
        }
      }
      // reset Ccopy
      Ccopy(k, l) -= delta / 2;
      Ccopy(l, k) -= delta / 2;
    }
  }
  //  // compare
  //  for (int i=0; i<3; i++)
  //    for (int j=0; j<3; j++)
  //      for (int k=0; k<3; k++)
  //        for (int l=0; l<3; l++)
  //          if(1)//abs(TensorDerivC[i][j][k][l])-abs(TensorDerivC_fd[i][j][k][l])>1e-6)
  //          {
  //            std::cout<<"dsqrt(C)/dC = "<<TensorDerivC[i][j][k][l]<<" Approx =
  //            "<<TensorDerivC_fd[i][j][k][l]<<" at ijkl="<<i<<j<<k<<l<<std::endl;
  //          }

  // Setup transposed matrices
  CORE::LINALG::Matrix<3, 3> Rtrans(true);
  Rtrans.UpdateT(R);
  CORE::LINALG::Matrix<3, 3> rotationratetrans(true);
  rotationratetrans.UpdateT(rotationrate);
  CORE::LINALG::Matrix<3, 3> defgrdratetrans(true);
  defgrdratetrans.UpdateT(defgrdrate);
  CORE::LINALG::Matrix<3, 3> invdefgrdtrans(true);
  invdefgrdtrans.UpdateT(invdefgrd);

  // 3x3x3x3 Tensor auxiliary variables
  CORE::LINALG::FourTensor<3> temptens1(true);
  CORE::LINALG::FourTensor<3> temptens3(true);
  CORE::LINALG::FourTensor<3> temptens4(true);
  CORE::LINALG::FourTensor<3> temptens5(true);
  CORE::LINALG::FourTensor<3> temptens7(true);
  CORE::LINALG::FourTensor<3> temptens8(true);
  CORE::LINALG::FourTensor<3> temptens10(true);
  CORE::LINALG::FourTensor<3> temptensgauss(true);
  // 3x3 matrix auxiliary variables
  CORE::LINALG::Matrix<3, 3> tempmat1(true);
  CORE::LINALG::Matrix<3, 3> tempmat2(true);
  // 6x6 matrix auxiliary variables
  CORE::LINALG::Matrix<6, 6> tempvoigt1(true);
  CORE::LINALG::Matrix<6, 6> tempvoigt2(true);
  CORE::LINALG::Matrix<6, 6> tempvoigtgauss(true);


  ///////////////////////////////////
  // Calculate constitutive tensor
  //////////////////////////////////
  // [F^-1 * sigma * d F^-T d C] = F^-1 * sigma * R * d sqrt(C)^-1/d C
  // (^T12 in assembly by switching i and j)
  tempmat1.MultiplyNN(invdefgrd, cauchystress);
  tempmat2.MultiplyNN(tempmat1, R);
  MAT::MultiplyMatrixFourTensor<3>(temptens1, tempmat2, RootCInvDerivCRootCInv, false);

  //  // F^-1 * sigma * R * d sqrt(C)^-1/d C
  //  tempmat1.MultiplyNN(invdefgrd,cauchystress);
  //  tempmat2.MultiplyNN(tempmat1,R);
  //  MultMatrixFourTensor(temptens2,tempmat2,RootCInvDerivCRootCInv,true);

  // (F^-1*sigma*F^-T) x dJ/dC : + F^{-1} \sigma F^{-T} dyad 0.5*C^{-1} (dyadic product to obtain
  // 4-Tensor)
  tempmat1.MultiplyNN(invdefgrd, cauchystress);
  tempmat2.MultiplyNT(tempmat1, invdefgrd);
  MAT::ElastSymTensorMultiply(tempvoigt1, 0.5, tempmat2, Cinv, 0.0);
  MAT::SetupFourTensor(temptens3, tempvoigt1);

  /////////////////////////////////////////////////////////////
  // velocity gradient with respect to right cauchy-green
  // d d / d C
  /////////////////////////////////////////////////////////////

  // F^-T * [R * d sqrt(C)/d C]^T12
  tempmat1.Update(invdefgrdtrans);
  tempmat1.Scale(1. / (theta * dt));
  MAT::MultiplyMatrixFourTensor<3>(temptens8, R, TensorDerivC, false);
  TransposeFourTensor12(temptens10, temptens8);
  MAT::MultiplyMatrixFourTensor<3>(temptens7, tempmat1, temptens10, false);

  // F^dot * [R * d sqrt(C)^-1 dC]^T12
  tempmat1.Update(defgrdrate);
  MAT::MultiplyMatrixFourTensor<3>(temptens8, R, RootCInvDerivCRootCInv, true);
  TransposeFourTensor12(temptens10, temptens8);
  MAT::MultiplyMatrixFourTensor<3>(temptens5, tempmat1, temptens10, false);

  for (int i = 0; i < 3; i++)
  {
    for (int j = 0; j < 3; j++)
    {
      for (int k = 0; k < 3; k++)
      {
        for (int l = 0; l < 3; l++)  // T12                                         // T12
          temptens8(i, j, k, l) = temptens7(j, i, k, l) + temptens5(i, j, k, l) +
                                  temptens5(j, i, k, l) + temptens7(i, j, k, l);
      }
    }
  }


  double onebyeps0 = 1. / epsilonnull;

  // parameters for sigmoid function
  //
  // 1/(1+exp(-(axct+b)))
  //
  double s = 0.96;
  double onebys = 1.0 / s;
  double aa = 2.0 * log(s / (1.0 - s)) * kvariance * onebyeps0;
  double b = -1.0 * log(1 / s - 1.0);

  // frequently used products and ratios
  double etafac = decayconst + theta * dt * (kforwards * Csignal + kbackwards);
  double onebyetafac = 1. / etafac;
  double auxfac2 = theta * dt * kbackwards;

  double onebysigmax = 1. / sigmamax;

  double deta_dsigma = onebyetafac * auxfac2 * onebysigmax * onebys;

  double exponent = 0.0;
  double efunct = 0.0;
  double denom = 0.0;
  ///////////////////////////////////////////////////////////
  // gauss integration over surface of a sphere
  ///////////////////////////////////////////////////////////
  for (int j = 0; j < 2 * gausspoints.nquad; j++)
  {
    for (int i = 0; i < gausspoints.nquad; i++)
    {
      double omega = acos(gausspoints.qxg[i][0]);
      double phi = ((double)(j + 1) * M_PI) / ((double)gausspoints.nquad);

      CORE::LINALG::Matrix<3, 1> m(true);
      m(0) = sin(omega) * cos(phi);
      m(1) = sin(omega) * sin(phi);
      m(2) = cos(omega);
      CORE::LINALG::Matrix<3, 3> M(true);
      M.MultiplyNT(m, m);
      CORE::LINALG::Matrix<numbgp, twice> epsomegaphi(true);

      // Transform strain rate at each point to fiber strain rate in (omega,phi) direction
      // \dot{\epsilon} = \dot{\epsilon}_{ij} m_{i} m_{j}
      epsomegaphi(i, j) = strainrate(0) * M(0, 0) + strainrate(1) * M(1, 1) +
                          strainrate(2) * M(2, 2) + 2. * strainrate(3) * M(0, 1) +
                          2. * strainrate(4) * M(1, 2) + 2. * strainrate(5) * M(0, 2);

      if (epsomegaphi(i, j) < 1e-12) epsomegaphi(i, j) = 0.0;

      exponent = 1.0 / etanew(i, j) * aa * epsomegaphi(i, j) + b;

      if (exponent <= 31.0 and exponent >= -30.0)
      {
        efunct = exp(-exponent);
        denom = 1.0 + efunct;
      }


      double factor = 0.0;

      if (exponent < -30.0 or exponent > 31.0)
        factor = 0.0;
      else
      {
        double temp = -denom * sigmamax * deta_dsigma +
                      sigmamax * efunct * aa * epsomegaphi(i, j) * deta_dsigma / etanew(i, j);
        factor = (3.0 / (16.0 * ((double)gausspoints.nquad))) *
                 ((sigmamax * efunct * aa) / (denom * denom + temp));
      }

      // Gauss quadrature
      for (int mm = 0; mm < 3; mm++)
      {
        for (int n = 0; n < 3; n++)
        {
          for (int o = 0; o < 3; o++)
          {
            for (int p = 0; p < 3; p++)
            {
              for (int q = 0; q < 3; q++)
              {
                for (int r = 0; r < 3; r++)
                  temptensgauss(mm, n, o, p) += temptens8(q, r, o, p) * m(q) * m(r) * factor *
                                                gausspoints.qwgt[i] * m(mm) * m(n);
              }
            }
          }
        }
      }

    }  // loop over i
  }    // loop over j

  // F^-1 * [F^-1 * sigma]^T12
  MAT::MultiplyMatrixFourTensor<3>(temptens8, invdefgrd, temptensgauss, true);
  TransposeFourTensor12(temptensgauss, temptens8);
  MAT::MultiplyMatrixFourTensor<3>(temptens4, invdefgrd, temptensgauss, false);


  // Put together active constitutive tensor
  for (int i = 0; i < 3; i++)
  {
    for (int j = 0; j < 3; j++)
    {
      for (int k = 0; k < 3; k++)
      {
        for (int l = 0; l < 3; l++)  // T12                                           // T12
          temptens5(i, j, k, l) = temptens1(i, j, k, l) + temptens1(j, i, k, l) +
                                  temptens3(i, j, k, l) + temptens4(j, i, k, l);
      }
    }
  }

  MAT::Setup6x6VoigtMatrix(cmatactive, temptens5);
  cmatactive.Scale(2.0 * detF);

}  // SetupCmatActive()

/*----------------------------------------------------------------------*
 |  Names of gp data to be visualized                       rauch  07/14|
 *----------------------------------------------------------------------*/
void MAT::ActiveFiber::VisNames(std::map<std::string, int>& names)
{
  std::string fiber = "Etahat";
  names[fiber] = 1;  // scalar
  fiber = "Etahor";
  names[fiber] = 1;  // scalar
  fiber = "Etaver";
  names[fiber] = 1;  // scalar
  fiber = "Etadiag";
  names[fiber] = 1;  // scalar
                     //  fiber = "dxx";
                     //  names[fiber] = 1; // scalar
                     //  fiber = "dyy";
                     //  names[fiber] = 1; // scalar
                     //  fiber = "dzz";
                     //  names[fiber] = 1; // scalar
                     //  fiber = "dxy";
                     //  names[fiber] = 1; // scalar
                     //  fiber = "dyz";
                     //  names[fiber] = 1; // scalar
                     //  fiber = "dxz";
                     //  names[fiber] = 1; // scalar
  matpassive_->VisNames(names);

}  // VisNames()

/*----------------------------------------------------------------------*
 |  gp data to be visualized                                rauch  07/14|
 *----------------------------------------------------------------------*/
bool MAT::ActiveFiber::VisData(
    const std::string& name, std::vector<double>& data, int numgp, int eleID)
{
  if (name == "Etahat")
  {
    if ((int)data.size() != 1) dserror("size mismatch");
    double temp = 0.0;
    for (int iter = 0; iter < numgp; iter++) temp += etahat_()->at(iter);
    data[0] = temp / numgp;
  }
  else if (name == "Etahor")
  {
    if ((int)data.size() != 1) dserror("size mismatch");
    double temp = 0.0;
    for (int iter = 0; iter < numgp; iter++) temp += etahor_()->at(iter);
    data[0] = temp / numgp;
  }
  else if (name == "Etaver")
  {
    if ((int)data.size() != 1) dserror("size mismatch");
    double temp = 0.0;
    for (int iter = 0; iter < numgp; iter++) temp += etaver_()->at(iter);
    data[0] = temp / numgp;
  }
  else if (name == "Etadiag")
  {
    if ((int)data.size() != 1) dserror("size mismatch");
    double temp = 0.0;
    for (int iter = 0; iter < numgp; iter++) temp += etadiag_()->at(iter);
    data[0] = temp / numgp;
  }
  //  else if (name == "dxx")
  //  {
  //    if ((int)data.size()!=1)
  //      dserror("size mismatch");
  //    double temp = 0.0;
  //    for (int iter=0; iter<numgp; iter++)
  //      temp += dxx_()->at(iter);
  //    data[0] = temp/numgp;
  //  }
  //  else if (name == "dyy")
  //  {
  //    if ((int)data.size()!=1)
  //      dserror("size mismatch");
  //    double temp = 0.0;
  //    for (int iter=0; iter<numgp; iter++)
  //      temp += dyy_()->at(iter);
  //    data[0] = temp/numgp;
  //  }
  //  else if (name == "dzz")
  //  {
  //    if ((int)data.size()!=1)
  //      dserror("size mismatch");
  //    double temp = 0.0;
  //    for (int iter=0; iter<numgp; iter++)
  //      temp += dzz_()->at(iter);
  //    data[0] = temp/numgp;
  //  }
  //  else if (name == "dxy")
  //  {
  //    if ((int)data.size()!=1)
  //      dserror("size mismatch");
  //    double temp = 0.0;
  //    for (int iter=0; iter<numgp; iter++)
  //      temp += dxy_()->at(iter);
  //    data[0] = temp/numgp;
  //  }
  //  else if (name == "dyz")
  //  {
  //    if ((int)data.size()!=1)
  //      dserror("size mismatch");
  //    double temp = 0.0;
  //    for (int iter=0; iter<numgp; iter++)
  //      temp += dyz_()->at(iter);
  //    data[0] = temp/numgp;
  //  }
  //  else if (name == "dxz")
  //  {
  //    if ((int)data.size()!=1)
  //      dserror("size mismatch");
  //    double temp = 0.0;
  //    for (int iter=0; iter<numgp; iter++)
  //      temp += dxz_()->at(iter);
  //    data[0] = temp/numgp;
  //  }
  else
  {
    return matpassive_->VisData(name, data, numgp, eleID);
  }
  return true;
}
