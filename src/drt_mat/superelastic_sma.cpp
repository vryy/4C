/*----------------------------------------------------------------------*/
/*!
\file superelastic_sma.cpp
\brief Material law for superelastic isotropic material
       following finite strain with linear and exponential
       flow rules.
       The implementation follows
       Auricchio F. and Taylor R.: "Shape-memory alloys: modelling and numerical
       simulations of finite-strain superelastic behavior."
       Computer Methods in Applied Mechanics and Engineering, 1997
           and
       Auricchio F.: A robust integration-algorithm for a finite-strain
       shape-memory-alloy superelastic model."
       International Journal of Plasticity, 2001


       geometrically nonlinear, finite strains, rate-independent, isothermal

       example input line for the exponential model:
       MAT 1 MAT_Struct_SuperElastSMA YOUNG 60000 DENS 1.0 NUE 0.3 EPSILON_L 0.075
       T_AS_s 0 T_AS_f 0 T_SA_s 0 T_SA_f 0 C_AS 0 C_SA 0 SIGMA_AS_s 520 SIGMA_AS_f 750
       SIGMA_SA_s 550 SIGMA_SA_f 200 ALPHA 0.15 MODEL 1 BETA_AS 250 BETA_SA 20

       example input line for the linear model:
       MAT 1 MAT_Struct_SuperElastSMA YOUNG 60000 DENS 1.0 NUE 0.3 EPSILON_L 0.075
       T_AS_s 0 T_AS_f 0 T_SA_s 0 T_SA_f 0 C_AS 0 C_SA 0 SIGMA_AS_s 520 SIGMA_AS_f 600
       SIGMA_SA_s 300 SIGMA_SA_f 200 ALPHA 0.15 MODEL 2 BETA_AS 0 BETA_SA 0

\level 3

\maintainer André Hemmler
            hemmler@mhpc.mw.tum.de
            http://www.mhpc.mw.tum.de
            089 - 289-1038


*/
/*----------------------------------------------------------------------*
 | headers                                                hemmler 09/16 |
 *----------------------------------------------------------------------*/
#include "superelastic_sma.H"
#include "matpar_bundle.H"
#include "material_service.H"

#include "../linalg/linalg_utils.H"
#include "../drt_lib/drt_globalproblem.H"


/*----------------------------------------------------------------------*
 | constructor (public)                                   hemmler 09/16 |
 *----------------------------------------------------------------------*/
MAT::PAR::SuperElasticSMA::SuperElasticSMA(
  Teuchos::RCP<MAT::PAR::Material> matdata
  )
: Parameter(matdata),
  density_(matdata->GetDouble("DENS")),
  youngs_(matdata->GetDouble("YOUNG")),
  poissonratio_(matdata->GetDouble("NUE")),
  epsilon_L_(matdata->GetDouble("EPSILON_L")),
  T_AS_s_(matdata->GetDouble("T_AS_s")),
  T_AS_f_(matdata->GetDouble("T_AS_f")),
  T_SA_s_(matdata->GetDouble("T_SA_s")),
  T_SA_f_(matdata->GetDouble("T_SA_f")),
  C_AS_(matdata->GetDouble("C_AS")),
  C_SA_(matdata->GetDouble("C_SA")),
  sigma_AS_s_(matdata->GetDouble("SIGMA_AS_s")),
  sigma_AS_f_(matdata->GetDouble("SIGMA_AS_f")),
  sigma_SA_s_(matdata->GetDouble("SIGMA_SA_s")),
  sigma_SA_f_(matdata->GetDouble("SIGMA_SA_f")),
  alpha_(matdata->GetDouble("ALPHA")),
  model_(matdata->GetInt("MODEL")),
  beta_AS_(matdata->GetDouble("BETA_AS")),
  beta_SA_(matdata->GetDouble("BETA_SA"))
{
}


/*----------------------------------------------------------------------*
 | struct with all material parameters                    hemmler 09/16 |
 *----------------------------------------------------------------------*/
struct MAT::SuperElasticSMA::material {
  // elastic material data
  double youngs;
  double poisson;
  double shear;
  double G1;
  double bulk;

  // superelastic material data
  double T_AS_s;
  double T_AS_f;
  double T_SA_s;
  double T_SA_f;
  double C_AS;
  double C_SA;
  double sigma_AS_s;
  double sigma_AS_f;
  double sigma_SA_s;
  double sigma_SA_f;
  double alpha;
  double epsilon_L;
  double model;
  double beta_AS;
  double beta_SA;
  double R_AS_s;
  double R_AS_f;
  double R_SA_s;
  double R_SA_f;
  double temperature;
};

struct MAT::SuperElasticSMA::loadingData {
  double drucker_prager;
  double drucker_prager_last;
  double drucker_prager_AS;
  double drucker_prager_AS_last;
  double drucker_prager_SA;
  double drucker_prager_SA_last;
  double F_AS_f;
  double F_SA_f;
  int H_AS;
  int H_SA;
};


/*----------------------------------------------------------------------*
 | is called in Material::Factory from ReadMaterials()    hemmler 09/16 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<MAT::Material> MAT::PAR::SuperElasticSMA::CreateMaterial()
{
  return Teuchos::rcp(new MAT::SuperElasticSMA(this));
}


MAT::SuperElasticSMAType MAT::SuperElasticSMAType::instance_;


/*----------------------------------------------------------------------*
 | is called in Material::Factory from ReadMaterials()    hemmler 09/16 |
 *----------------------------------------------------------------------*/
DRT::ParObject* MAT::SuperElasticSMAType::Create(const std::vector<char>& data)
{
  MAT::SuperElasticSMA* superelast = new MAT::SuperElasticSMA();
  superelast->Unpack(data);
  return superelast;
}


/*----------------------------------------------------------------------*
 | constructor (public)                                   hemmler 09/16 |
 *----------------------------------------------------------------------*/
MAT::SuperElasticSMA::SuperElasticSMA()
  : params_(NULL)
{
}


/*----------------------------------------------------------------------*
 | copy-constructor (public)                              hemmler 09/16 |
 *----------------------------------------------------------------------*/
MAT::SuperElasticSMA::SuperElasticSMA(MAT::PAR::SuperElasticSMA* params)
  : params_(params)
{
}


/*----------------------------------------------------------------------*
 | pack (public)                                          hemmler 09/16 |
 *----------------------------------------------------------------------*/
void MAT::SuperElasticSMA::Pack(DRT::PackBuffer& data) const
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

  // pack history data
  int histsize;
  // if material is not initialized, i.e. start simulation, nothing to pack
  if (!Initialized())
  {
    histsize=0;
  }
  else
  {
    // if material is initialized (restart): size equates number of gausspoints
    histsize = xi_s_curr_->size();
  }
  AddtoPack(data,histsize); // Length of history vector(s)
  for (int var=0; var<histsize; ++var)
  {
    // insert history vectors to AddtoPack
    AddtoPack(data, druckerpragerloadingcurr_->at(var));
    AddtoPack(data, druckerpragerloadinglast_->at(var));
    AddtoPack(data, xi_s_curr_->at(var));
    AddtoPack(data, xi_s_last_->at(var));
  }

  return;
}  // Pack()


/*----------------------------------------------------------------------*
 | unpack (public)                                        hemmler 09/16 |
 *----------------------------------------------------------------------*/
void MAT::SuperElasticSMA::Unpack(const std::vector<char>& data)
{
  isinit_=true;
  std::vector<char>::size_type position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position,data,type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");

  // matid
  int matid;
  ExtractfromPack(position,data,matid);
  params_ = NULL;
  if (DRT::Problem::Instance()->Materials() != Teuchos::null)
    if (DRT::Problem::Instance()->Materials()->Num() != 0)
    {
      const int probinst = DRT::Problem::Instance()->Materials()->GetReadFromProblem();
      MAT::PAR::Parameter* mat = DRT::Problem::Instance(probinst)->Materials()->ParameterById(matid);
      if (mat->Type() == MaterialType())
        params_ = static_cast<MAT::PAR::SuperElasticSMA*>(mat);
      else
        dserror("Type of parameter material %d does not fit to calling type %d", mat->Type(), MaterialType());
    }

  // history data
  int histsize;
  ExtractfromPack(position,data,histsize);

  // if system is not yet initialized, the history vectors have to be intialized
  if (histsize == 0)
    isinit_ = false;
  druckerpragerloadinglast_ = Teuchos::rcp(new std::vector<double>);
  druckerpragerloadingcurr_ = Teuchos::rcp(new std::vector<double>);
  xi_s_curr_ = Teuchos::rcp(new std::vector<double>);
  xi_s_last_ = Teuchos::rcp(new std::vector<double>);


  for (int var=0; var<histsize; ++var)
  {

    double tmpDouble;

    // vectors of last converged state are unpacked

    ExtractfromPack(position,data,tmpDouble);
    druckerpragerloadingcurr_->push_back(tmpDouble);

    ExtractfromPack(position,data,tmpDouble);
    druckerpragerloadinglast_->push_back(tmpDouble);

    ExtractfromPack(position,data,tmpDouble);
    xi_s_curr_->push_back(tmpDouble);

    ExtractfromPack(position,data,tmpDouble);
    xi_s_last_->push_back(tmpDouble);
  }

  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d",data.size(),position);

  return;

}  // Unpack()


/*---------------------------------------------------------------------*
 | initialize / allocate internal variables (public)     hemmler 09/16 |
 *---------------------------------------------------------------------*/
void MAT::SuperElasticSMA::Setup(int numgp, DRT::INPUT::LineDefinition* linedef)
{
  druckerpragerloadingcurr_ = Teuchos::rcp(new std::vector<double>);
  druckerpragerloadinglast_ = Teuchos::rcp(new std::vector<double>);
  xi_s_curr_ = Teuchos::rcp(new std::vector<double>);
  xi_s_last_ = Teuchos::rcp(new std::vector<double>);

  druckerpragerloadingcurr_->resize(numgp);
  druckerpragerloadinglast_->resize(numgp);
  xi_s_curr_->resize(numgp);
  xi_s_last_->resize(numgp);

  for (int i=0; i<numgp; i++)
  {
    druckerpragerloadinglast_->at(i) = 0.0;
    druckerpragerloadingcurr_->at(i) = 0.0;

    xi_s_curr_->at(i) = 0.0; // Start with zero single variant martensitic fraction
    xi_s_last_->at(i) = 0.0; // Start with zero single variant martensitic fraction
  }

  isinit_ = true;
  return;

}  // Setup()


/*----------------------------------------------------------------------*
 | update internal variables                              hemmler 09/16 |
 *----------------------------------------------------------------------*/
void MAT::SuperElasticSMA::Update()
{
  druckerpragerloadinglast_ = druckerpragerloadingcurr_;
  xi_s_last_ = xi_s_curr_;

  druckerpragerloadingcurr_ = Teuchos::rcp(new std::vector<double >);
  xi_s_curr_ = Teuchos::rcp(new std::vector<double >);

  // Empty vectors of current data
  const int histsize = xi_s_last_->size();
  druckerpragerloadingcurr_->resize(histsize);
  xi_s_curr_->resize(histsize);

  for (int i=0;i<histsize;i++) {
    druckerpragerloadingcurr_->at(i) = 0.0;
    xi_s_curr_->at(i) = 0.0;
  }
  return;
}  // Update()


/*----------------------------------------------------------------------*
 | calculate stress and constitutive tensor               hemmler 09/16 |
 *----------------------------------------------------------------------*/
void MAT::SuperElasticSMA::Evaluate(
  const LINALG::Matrix<3,3>* defgrd,
  const LINALG::Matrix<6,1>* glstrain,
  Teuchos::ParameterList& params,
  LINALG::Matrix<6,1>* stress,
  LINALG::Matrix<6,6>* cmat,
  const int eleGID
  )
{
  // extract the gauss points from the parameter list
  const int gp = params.get<int>("gp",-1);
  if (gp == -1) dserror("no Gauss point number provided in material");


  /*
   **********************************************************
   * Step 0.1 * Load material data from the input file        *
   **********************************************************
   */
  material matdata;
  // elastic material data
  matdata.youngs = params_->youngs_;  // Young's modulus
  matdata.poisson = params_->poissonratio_;  // Poisson's ratio
  matdata.shear = matdata.youngs / (2.0 * (1.0 + matdata.poisson));  // shear modulus
  matdata.bulk = matdata.youngs / (3.0 * (1.0 - 2.0 * matdata.poisson));  // bulk modulus


  // superelastic material data
  matdata.T_AS_s = params_->T_AS_s_;
  matdata.T_AS_f = params_->T_AS_f_;
  matdata.T_SA_s = params_->T_SA_s_;
  matdata.T_SA_f = params_->T_SA_f_;
  matdata.C_AS = params_->C_AS_;
  matdata.C_SA = params_->C_SA_;
  matdata.sigma_AS_s = params_->sigma_AS_s_;
  matdata.sigma_AS_f = params_->sigma_AS_f_;
  matdata.sigma_SA_s = params_->sigma_SA_s_;
  matdata.sigma_SA_f = params_->sigma_SA_f_;
  matdata.alpha = params_->alpha_;
  matdata.epsilon_L = params_->epsilon_L_ / (std::sqrt(2.0/3.0)+matdata.alpha);
  matdata.model = params_->model_;
  matdata.temperature = 0.0;


  // material data that is only needed for the exponential model
  matdata.beta_AS = params_->beta_AS_ * (std::sqrt(2.0/3.0)+matdata.alpha); // Eqn. 44
  matdata.beta_SA = params_->beta_SA_ * (std::sqrt(2.0/3.0)+matdata.alpha); // Eqn. 45
  matdata.G1 = (2.0 * matdata.shear + 9.0 * matdata.bulk * std::pow(matdata.alpha,2.0))*matdata.epsilon_L; // Eqn. 80, 81

  // Compute material parameters
  matdata.R_AS_s = matdata.sigma_AS_s*(std::sqrt(2.0/3.0)+matdata.alpha)-matdata.C_AS * matdata.T_AS_s; // Eqn. 6
  matdata.R_AS_f = matdata.sigma_AS_f*(std::sqrt(2.0/3.0)+matdata.alpha)-matdata.C_AS * matdata.T_AS_f; // Eqn. 7

  matdata.R_SA_s = matdata.sigma_SA_s*(std::sqrt(2.0/3.0)+matdata.alpha)-matdata.C_SA * matdata.T_SA_s; // Eqn. 6
  matdata.R_SA_f = matdata.sigma_SA_f*(std::sqrt(2.0/3.0)+matdata.alpha)-matdata.C_SA * matdata.T_SA_f; // Eqn. 7

  // Check whether a model is given (exponential or linear)
  if(matdata.model != 1 and matdata.model != 2) {
    // No proper model is given --> throw error
    dserror("No sma-model given. Use 1 for the exponential model and 2 for the linear model.");
  }

  LINALG::Matrix<3,3> deformation_gradient_invert(*defgrd);
  deformation_gradient_invert.Invert();

  /*
   **********************************************************
   * Step 0.2 * Read state of the last time-step            *
   **********************************************************
   */
  double xi_S = xi_s_last_->at(gp);
  double druckerpragerloadinglast = druckerpragerloadinglast_->at(gp);


  /*
   **********************************************************
   * Step 1.1 * Compute the Cauchy-Green tensor and its     *
   *          * eigenvalues and -vectors.                   *
   **********************************************************
   */
  // b = F * F^T
  LINALG::Matrix<3,3> cauchy_green_tensor(true);
  cauchy_green_tensor.MultiplyNT(*defgrd, *defgrd);

  // To compute the spectral decomposition, the Cauchy-Green tensor
  // must be converted to Epetra format
  Epetra_SerialDenseMatrix cauchy_green_eigenvectors(3,3);
  Epetra_SerialDenseVector cauchy_green_eigenvalues(3);

  for (int i=0;i<3;i++)
    for (int j=0;j<3;j++)
      cauchy_green_eigenvectors(i,j) = cauchy_green_tensor(i,j);

  // Solve the eigen-problem
  LINALG::SymmetricEigenProblem(cauchy_green_eigenvectors, cauchy_green_eigenvalues);

  /*
   **********************************************************
   * Step 1.2 * Compute logarithmic strains and scaled      *
   *          * principal directions                        *
   **********************************************************
   */

  double logarithmic_strain_volumetric;
  LINALG::Matrix<3,3> logarithmic_strain_deviatoric_tensor(true);
  LINALG::Matrix<3,3> material_scaled_load_deviatoric_tensor(true);
  double logarithmic_strain_deviatoric_norm = 0.0;
  std::vector< LINALG::Matrix<3,1> > spatial_principal_directions(true);
  spatial_principal_directions.resize(3);
  std::vector< LINALG::Matrix<3,3> > spatial_principal_matrices(true);
  std::vector< LINALG::Matrix<3,3> > material_principal_matrices(true);
  spatial_principal_matrices.resize(3);
  material_principal_matrices.resize(3);

  for (int i=0;i<3;i++) {

    for (int j=0;j<3;j++) {
      spatial_principal_directions.at(i)(j) = cauchy_green_eigenvectors(j,i);
    }

    LINALG::Matrix<3,1> material_principal_direction(true);
    material_principal_direction.MultiplyNN(deformation_gradient_invert, spatial_principal_directions.at(i));
    material_principal_direction.Scale(cauchy_green_eigenvalues(i));

    spatial_principal_matrices.at(i).MultiplyNT(spatial_principal_directions.at(i), spatial_principal_directions.at(i));
    material_principal_matrices.at(i).MultiplyNT(material_principal_direction, material_principal_direction);
  }

  // Compute Jacobian
  double cauchy_green_jacobian = std::sqrt(cauchy_green_eigenvalues(0)*cauchy_green_eigenvalues(1)*cauchy_green_eigenvalues(2));


  // Compute logrithmic strains
  logarithmic_strain_volumetric = std::log(cauchy_green_jacobian);

  for (int i=0;i<3;i++) {
    // Compute the deviatoric part of the eigenvalues of the CauchyGreenTensor
    double cauchy_green_principal_deviatoric = std::pow( cauchy_green_jacobian , -1.0 / 3.0 ) * std::sqrt( cauchy_green_eigenvalues(i) );

    double logarithmic_strain_principal_deviatoric = std::log( cauchy_green_principal_deviatoric );

    logarithmic_strain_deviatoric_norm += std::pow( logarithmic_strain_principal_deviatoric , 2.0 );

    // Compute the deviatoric logarithmic strain matrix
    logarithmic_strain_deviatoric_tensor.Update(logarithmic_strain_principal_deviatoric, spatial_principal_matrices.at(i), 1.0);
    material_scaled_load_deviatoric_tensor.Update(logarithmic_strain_principal_deviatoric, spatial_principal_matrices.at(i), 1.0);
  }
  logarithmic_strain_deviatoric_norm = std::sqrt(logarithmic_strain_deviatoric_norm);

  LINALG::Matrix<3,3> scaled_load_deviatoric_tensor(logarithmic_strain_deviatoric_tensor);

  if(logarithmic_strain_deviatoric_norm != 0.0) {
    scaled_load_deviatoric_tensor.Scale( 1.0 / logarithmic_strain_deviatoric_norm );
    material_scaled_load_deviatoric_tensor.Scale( 1.0 / logarithmic_strain_deviatoric_norm );
  }

  /*
   **********************************************************
   * Step 2.1 * Compute trial state of the Drucker-Prager-  *
   *          * type loading function                       *
   **********************************************************
   */
  double kirchhoff_trial_stress_volumetric = matdata.bulk * ( logarithmic_strain_volumetric - 3.0 * matdata.alpha * matdata.epsilon_L * xi_S );
  double kirchhoff_trial_stress_deviatoric_norm = 2.0 * matdata.shear * std::abs( logarithmic_strain_deviatoric_norm - matdata.epsilon_L * xi_S );

  double drucker_prager_loading_trial = kirchhoff_trial_stress_deviatoric_norm + 3.0 * matdata.alpha * kirchhoff_trial_stress_volumetric;

  double drucker_prager_loading_AS_trial = drucker_prager_loading_trial - matdata.C_AS * matdata.temperature;
  double drucker_prager_loading_SA_trial = drucker_prager_loading_trial - matdata.C_SA * matdata.temperature;

  double drucker_prager_loading_AS_last = druckerpragerloadinglast - matdata.C_AS * matdata.temperature;
  double drucker_prager_loading_SA_last = druckerpragerloadinglast - matdata.C_SA * matdata.temperature;

  /*
   **********************************************************
   * Step 2.2 - 4 * Check phase transformations             *
   **********************************************************
   */
  double lambda_AS = 0.0;
  double lambda_SA = 0.0;
  int H_AS = 0;
  int H_SA = 0;
  bool solution_found = false;

  // Check solution xi_S = 0
  double kirchhoff_stress_volumetric_zero = matdata.bulk * logarithmic_strain_volumetric;
  double kirchhoff_stress_deviatoric_norm_zero = 2.0 * matdata.shear * logarithmic_strain_deviatoric_norm;
  double drucker_prager_loading_zero_SA = kirchhoff_stress_deviatoric_norm_zero + 3.0 * matdata.alpha * kirchhoff_stress_volumetric_zero - matdata.C_SA * matdata.temperature;
  if (drucker_prager_loading_zero_SA < matdata.R_SA_f) {
    // xi_S = 0 is the appropriate solution
    xi_S = 0;
    H_SA = 1;
    lambda_SA = xi_S;
    solution_found = true;
  }

  // Check solution xi_S = 1
  double kirchhoff_stress_volumetric_full = matdata.bulk * ( logarithmic_strain_volumetric - 3 * matdata.alpha * matdata.epsilon_L );
  double kirchhoff_stress_deviatoric_norm_full = 2.0 * matdata.shear * ( logarithmic_strain_deviatoric_norm - matdata.epsilon_L );
  double drucker_prager_loading_full_AS = kirchhoff_stress_deviatoric_norm_full + 3.0 * matdata.alpha * kirchhoff_stress_volumetric_full - matdata.C_AS * matdata.temperature;
  if (drucker_prager_loading_full_AS > matdata.R_AS_f) {
    // xi_S = 0 is the appropriate solution
    xi_S = 1;
    H_AS = 1;
    lambda_AS = 1.0 - xi_S;
    solution_found = true;
  }

  if (!solution_found) {

    // Check A -> S phase transformation
    if(drucker_prager_loading_AS_trial > matdata.R_AS_s && drucker_prager_loading_AS_trial > drucker_prager_loading_AS_last && xi_S < 1.0) {
      H_AS = 1;
    }

    // Check S -> A phase transformation
    if(drucker_prager_loading_SA_trial < matdata.R_SA_s && drucker_prager_loading_SA_trial < drucker_prager_loading_SA_last && xi_S > 0.0) {
      H_SA = 1;
    }

  }

  /*
   **********************************************************
   * Step 5 * Compute martensitic evolution, if at least    *
   *        * one phase transformation is active            *
   **********************************************************
   */

  if(!solution_found && (H_AS == 1 || H_SA == 1)) {
    // Return back and update martensitic fraction

    LINALG::Matrix<2,1> lambda_S(true);


    // Parameter for the damped Newton
    int innerdamp_iter_max = 5;
    int damp_iter = 0;
    int damp_iter_max = 100;
    double damping = 1.0;

    bool converged = false;
    int iter = 0;
    int maxiter = 1000;
    double tol = 1e-9;
    double fp;
    double fk;


    while(damp_iter < damp_iter_max) {
      iter = 0;

      LINALG::Matrix<2,1> res_vec(true);
      LINALG::Matrix<2,1> R_k(true);
      LINALG::Matrix<2,1> R_p(true);
      LINALG::Matrix<2,2> d_R_d_lambda(true);
      LINALG::Matrix<2,2> d_R_d_lambda_inv(true);
      lambda_S(0) = lambda_AS;
      lambda_S(1) = lambda_SA;
      LINALG::Matrix<2,1> lambda_S_p(lambda_S);



      while (iter < maxiter) {
        int innerdamp_iter = 0;
        double gamma = 1.0;
        double xi_s_tmp = xi_S + lambda_S(0) + lambda_S(1);

        bool accept = false;

        // Compute local Newton step
        loadingData loading = ComputeLocalNewtonLoading(xi_s_tmp, logarithmic_strain_volumetric, logarithmic_strain_deviatoric_norm, matdata);
        loading.drucker_prager_AS_last = drucker_prager_loading_AS_last;
        loading.drucker_prager_SA_last = drucker_prager_loading_SA_last;
        loading.drucker_prager_last = druckerpragerloadinglast;
        loading.H_AS = H_AS;
        loading.H_SA = H_SA;
        R_k = ComputeLocalNewtonResidual(lambda_S, xi_s_tmp, loading, matdata);
        d_R_d_lambda = ComputeLocalNewtonJacobian(lambda_S, xi_s_tmp, loading, matdata);
        d_R_d_lambda_inv.Invert(d_R_d_lambda);
        res_vec.MultiplyNN(d_R_d_lambda_inv, R_k); // Eqn. 57

        // Compute function, that has to be minimized
        fk = std::pow(R_k(0), 2.0) + std::pow(R_k(1), 2.0);




        do {
          lambda_S_p.Update(1.0, lambda_S, -gamma * damping, res_vec); // Eqn. 57
          double xi_s_p = xi_S + lambda_S_p(0) + lambda_S_p(1);
          loadingData loading_p = ComputeLocalNewtonLoading(xi_s_p, logarithmic_strain_volumetric, logarithmic_strain_deviatoric_norm, matdata);
          loading_p.drucker_prager_AS_last = drucker_prager_loading_AS_last;
          loading_p.drucker_prager_SA_last = drucker_prager_loading_SA_last;
          loading_p.drucker_prager_last = druckerpragerloadinglast;
          loading_p.H_AS = H_AS;
          loading_p.H_SA = H_SA;
          R_p = ComputeLocalNewtonResidual(lambda_S_p, xi_s_p, loading_p, matdata);

          fp = std::pow(R_p(0), 2.0) + std::pow(R_p(1), 2.0);

          // Check condition for acceptance of the Newton-step
          if (fp <= (1 - gamma / 2.0)*fk) {
            accept = true;
          } else {
            gamma *= 0.5;
          }
          innerdamp_iter++;
        } while (!accept && innerdamp_iter < innerdamp_iter_max);

        // Updagte lambda_S with the one calculated with Armijo
        lambda_S.Update(1.0, lambda_S_p); // Eqn. 57

        if(fp < tol) {
          converged = true;
          break;
        }

        iter++;
      }
      damp_iter++;
      if(!converged) {
        damping *= 0.5;
      }

    }


    if(!converged) {
      // local newton method unconverged
      bool error_tol = false;
      if(params.isParameter("tolerate_errors")) {
        error_tol = params.get<bool>("tolerate_errors");
      }


      if(error_tol) {
        params.set<bool>("eval_error", true);
      } else {
        dserror("Local Newton iteration unconverged in %i iterations. Residual: %e Tol: %e", maxiter, fp, tol);
      }
    }

    lambda_AS = lambda_S(0);
    lambda_SA = lambda_S(1);

    xi_S = xi_S + lambda_AS + lambda_SA;
  }
  /*
   **********************************************************
   * Step 6 * Compute stresses with the new martensitic     *
   *        * fraction                                      *
   **********************************************************
   */


  double kirchhoff_stress_volumetric = matdata.bulk * ( logarithmic_strain_volumetric - 3.0 * matdata.alpha * matdata.epsilon_L * xi_S );
  double kirchhoff_stress_deviatoric_norm = 2.0 * matdata.shear * ( logarithmic_strain_deviatoric_norm - matdata.epsilon_L * xi_S );
  LINALG::Matrix<3,3> kirchhoff_stress_deviatoric(scaled_load_deviatoric_tensor);
  kirchhoff_stress_deviatoric.Scale( kirchhoff_stress_deviatoric_norm );

  LINALG::Matrix<3,3> kirchhoff_stress(kirchhoff_stress_deviatoric);
  for (int i=0;i<3;i++)
    kirchhoff_stress(i,i) += kirchhoff_stress_volumetric;

  // Convert Kirchhoff stress tensor in the second Piola-Kirchhoff tensor
  LINALG::Matrix<3,3> PK2(true);

  PK2.MultiplyNN(deformation_gradient_invert, kirchhoff_stress);
  PK2.MultiplyNT(PK2, deformation_gradient_invert);

  // Convert PK2 into Voigt-Notation
  (*stress)(0) = PK2(0,0);
  (*stress)(1) = PK2(1,1);
  (*stress)(2) = PK2(2,2);
  (*stress)(3) = 0.5 * (PK2(0,1) + PK2(1,0));
  (*stress)(4) = 0.5 * (PK2(1,2) + PK2(2,1));
  (*stress)(5) = 0.5 * (PK2(2,0) + PK2(0,2));

  /*
   **********************************************************
   * Step 7 * Compute algorithmic tangent for the global    *
   *        * Newton method                                 *
   **********************************************************
   */
  double drucker_prager_loading = kirchhoff_stress_deviatoric_norm + 3.0 * matdata.alpha * kirchhoff_stress_volumetric;

  double drucker_prager_loading_AS = drucker_prager_loading - matdata.C_AS * matdata.temperature;
  double drucker_prager_loading_SA = drucker_prager_loading - matdata.C_SA * matdata.temperature;

  double F_AS_f = drucker_prager_loading_AS - matdata.R_AS_f;
  double F_SA_f = drucker_prager_loading_SA - matdata.R_SA_f;

  double a;
  double b;
  double c;
  double d;

  double A_AS;
  double A_SA;
  if(matdata.model == 1) {
    // exponential model
    a = H_AS * matdata.beta_AS * (drucker_prager_loading_AS - drucker_prager_loading_AS_last) + matdata.G1 * ( H_AS * matdata.beta_AS * ( 1.0 - xi_S ) - 2.0 * lambda_AS * F_AS_f ) + std::pow(F_AS_f,2.0); // Eqn. 88
    b = H_AS * matdata.beta_AS * (drucker_prager_loading_AS - drucker_prager_loading_AS_last) + matdata.G1 * ( H_AS * matdata.beta_AS * ( 1.0 - xi_S ) - 2.0 * lambda_AS * F_AS_f ); // Eqn. 89
    c = -H_SA * matdata.beta_SA * (drucker_prager_loading_SA - drucker_prager_loading_SA_last) +  matdata.G1 * ( H_SA * matdata.beta_SA * xi_S - 2.0 *  lambda_SA * F_SA_f); // Eqn. 90
    d = -H_SA * matdata.beta_SA * (drucker_prager_loading_SA - drucker_prager_loading_SA_last) +  matdata.G1 * ( H_SA * matdata.beta_SA * xi_S - 2.0 *  lambda_SA * F_SA_f) + std::pow(F_SA_f,2.0); // Eqn. 91

    A_AS = -2.0 * lambda_AS * F_AS_f +  H_AS * matdata.beta_AS * (1.0-xi_S); // Eqn. 92
    A_SA = -2.0 * lambda_SA * F_SA_f +  H_SA * matdata.beta_SA * xi_S; // Eqn. 93
  } else {
    // linear model
    a = -matdata.G1 * lambda_AS - H_AS * ( ( drucker_prager_loading_AS - drucker_prager_loading_AS_last ) + (1.0 - xi_S) * matdata.G1 ) + F_AS_f; // Eqn. 82
    b = -matdata.G1 * lambda_AS - H_AS * ( ( drucker_prager_loading_AS - drucker_prager_loading_AS_last ) + (1.0 - xi_S) * matdata.G1 ); // Eqn. 83
    c = -matdata.G1 * lambda_SA - H_SA * ( ( drucker_prager_loading_SA - drucker_prager_loading_SA_last ) - xi_S * matdata.G1); // Eqn. 84
    d = -matdata.G1 * lambda_SA - H_SA * ( ( drucker_prager_loading_SA - drucker_prager_loading_SA_last ) - xi_S * matdata.G1) + F_SA_f; // Eqn. 85

    A_AS = -lambda_AS - H_AS * ( 1.0 - xi_S ); // Eqn. 86
    A_SA = -lambda_SA + H_SA * xi_S; // Eqn. 87
  }


  double det = a*d-c*b;
  double B = d/det; // Eqn. 58
  double C = -b/det; // Eqn. 58
  double D = -c/det; // Eqn. 58
  double E = a/det; // Eqn. 58

  double T_AS_1 = 2.0 * matdata.shear * (B*A_AS+C*A_SA); // Eqn. 68
  double T_AS_2 = 3.0 * matdata.bulk * matdata.alpha * ( B * A_AS + C * A_SA ); // Eqn. 69
  double T_SA_1 = 2.0 * matdata.shear * ( D * A_AS + E * A_SA ); // Eqn. 70
  double T_SA_2 = 3.0 * matdata.bulk * matdata.alpha * ( D * A_AS + E * A_SA ); // Eqn. 71

  double kappa_star = matdata.bulk * ( 1.0 - 3.0 * matdata.epsilon_L * matdata.alpha * (T_AS_2 + T_SA_2)); // Eqn. 76
  double G_star;
  double M1_star;

  if(xi_S > 0) {

    G_star = matdata.shear * ( 1.0 - matdata.epsilon_L * xi_S / logarithmic_strain_deviatoric_norm ); // Eqn. 77
    M1_star = 2.0 * matdata.shear * matdata.epsilon_L * ( xi_S / logarithmic_strain_deviatoric_norm - (T_AS_1 + T_SA_1) ); // Eqn. 78

  } else {

    G_star = matdata.shear; // Eqn. 77
    M1_star = -2.0 * matdata.shear * matdata.epsilon_L * (T_AS_1 + T_SA_1); // Eqn. 78

  }
  double M2_star = 2.0 * matdata.shear * matdata.epsilon_L * (T_AS_2 + T_SA_2); // Eqn. 79

  // D = [ K * (I (x) I) + 2G*I_dev + M1* (n (x) n) - M2*(n(x)I+1(x)n)]
  // Compute the algorithmic tangent incl. conversion into Voigt-notation



  if (xi_S>0.0) {
    LINALG::Matrix<3,3> eye(true);
    for (int i=0;i<3;i++)
      eye(i,i) = 1.0;
    LINALG::Matrix<6,6> cmat_eul(true);
    LINALG::Matrix<6,6> cmat_eul_tmp(true);
    LINALG::Matrix<6,6> cmat_eul_1(true);
    LINALG::Matrix<6,6> cmat_eul_2(true);
    LINALG::Matrix<6,6> cmat_eul_3(true);
    LINALG::Matrix<6,6> cmat_eul_4(true);
    LINALG::Matrix<6,6> cmat_eul_4_tmp(true);

    // Build up (1 (x) 1) and scale with K*
    for (int i=0;i<3;i++)
      for (int j=0;j<3;j++)
        cmat_eul_1(i,j) = 1.0;
    cmat_eul_1.Scale(kappa_star);

    // Build up I_dev and scale with 2G*
    for (int i=0;i<3;i++)
      for (int j=0;j<3;j++)
        if(i == j)
          cmat_eul_2(i,j) = 2.0 / 3.0;
        else
          cmat_eul_2(i,j) = -1.0 / 3.0;
    for (int i=3;i<6;i++)
      cmat_eul_2(i,i) = 0.5;
    cmat_eul_2.Scale( 2.0 * G_star );

    // Build up (n (x) n) and scale with M1*
    ElastSymTensorMultiply(cmat_eul_3, M1_star, scaled_load_deviatoric_tensor, scaled_load_deviatoric_tensor, 0.0);

    // Build up (n (x) 1 + 1 (x) n) and scale with M2*
    ElastSymTensorMultiply(cmat_eul_4, -M2_star, scaled_load_deviatoric_tensor, eye, 0.0);
    ElastSymTensorMultiply(cmat_eul_4_tmp, -M2_star, eye, scaled_load_deviatoric_tensor, 0.0);
    cmat_eul_4.Update(1.0, cmat_eul_4_tmp, 1.0);


    cmat_eul.Update(1.0, cmat_eul_1, 1.0);
    cmat_eul.Update(1.0, cmat_eul_2, 1.0);
    cmat_eul.Update(1.0, cmat_eul_3, 1.0);
    cmat_eul.Update(1.0, cmat_eul_4, 1.0);


    Pullback4thTensorVoigt(cauchy_green_jacobian, deformation_gradient_invert, cmat_eul, cmat);
  } else {
    // matrices for temporary stuff
     LINALG::Matrix<3,3> tmp1;
     LINALG::Matrix<3,3> tmp2;

     // 3x3 2nd-order identity matrix
     LINALG::Matrix<3,3> id2(true);
     LINALG::Matrix<3,3> Idev;
     for (int i=0; i<3; i++)
     {
       id2(i,i)=1.0;
       for (int j=0; j<3; j++)
       {
         if (i==j) Idev(i,j) = 2.0/3.0;
         else Idev(i,j) = -1.0/3.0;
       }
     }

     // linear elasticity tensor in principal directions
     LINALG::Matrix<3,3> D_ep_principal(Idev);
     D_ep_principal.Scale(2.0 * matdata.shear);
     for (int i=0; i<3; i++)
       for (int j=0; j<3; j++)
         D_ep_principal(i,j) += matdata.bulk;

     LINALG::Matrix<3,1> dev_KH(true);
     Epetra_SerialDenseVector lambda_trial_square(3);
     std::vector<LINALG::Matrix<3,1> > material_principal_directions;
     material_principal_directions.resize(3);

     double pressure=0.0;
     pressure=kirchhoff_stress_volumetric;

     for (int i=0; i<3; i++)
     {
       material_principal_directions.at(i).Multiply(deformation_gradient_invert,spatial_principal_directions.at(i));
       double logStrainPrincipal = std::log(std::pow(cauchy_green_jacobian, -1.0 / 3.0) * std::sqrt(cauchy_green_eigenvalues(i)));
       double logLoadingDirection = 0.0;
       if(logarithmic_strain_deviatoric_norm > 0) {
         logLoadingDirection = logStrainPrincipal / logarithmic_strain_deviatoric_norm;
       }
       dev_KH(i)=2.0 * matdata.shear * ( logStrainPrincipal - matdata.epsilon_L * xi_S * logLoadingDirection );
       lambda_trial_square(i)=cauchy_green_eigenvalues(i);
     }

     // express coefficents of tangent in Kirchhoff stresses
     cmat->Clear();
     for (int a=0; a<3; a++)
     {
       // - sum_1^3 (2 * tau N_aaaa)
       tmp1.MultiplyNT(material_principal_directions.at(a),material_principal_directions.at(a));  // N_{aa}
       ElastSymTensorMultiply(*cmat, -2.0 * (dev_KH(a) + pressure), tmp1, tmp1, 1.0);

       for (int b=0; b<3; b++)
       {
         // c_ab N_aabb
         // result of return mapping of deviatoric component c_ab
         tmp1.MultiplyNT(material_principal_directions.at(a),material_principal_directions.at(a));  // N_{aa}
         tmp2.MultiplyNT(material_principal_directions.at(b),material_principal_directions.at(b));  // N_{bb}
         ElastSymTensorMultiply(*cmat, D_ep_principal(a,b), tmp1, tmp2, 1.0);

         if (a != b)
         {
           double fac = 0.0;
           if (lambda_trial_square(a) != lambda_trial_square(b))
           {
             // (tau_aa * lambda_b^2 - tau_bb * lambda_a^2) / (lambda_a^2 - lambda_b^2)
             fac = ( (dev_KH(a) + pressure) * lambda_trial_square(b)
                    -(dev_KH(b) + pressure) * lambda_trial_square(a) )
                           / (lambda_trial_square(a) - lambda_trial_square(b) );
           } // end if lambda_a != lambda_b
           else // lambda_a = lambda_b
           {
             // 1/2 [(d^2 Psi)/(d ln lambda_b * d ln lambda_b) - (d^2 Psi)/(d ln lambda_a * d ln lambda_b)]
             // - tau_bb, cf. (6.91)
             fac = 0.5 * (D_ep_principal(b,b) - D_ep_principal(a,b)) - (dev_KH(b) + pressure);
           } // end lambda_a = lambda_b
           tmp1.MultiplyNT(material_principal_directions.at(a),material_principal_directions.at(b));  // N_{ab}
           tmp2.MultiplyNT(material_principal_directions.at(b),material_principal_directions.at(a));  // N_{ba}
           ElastSymTensorMultiply(*cmat,fac,tmp1,tmp1,1.0);  // N_{abab}
           ElastSymTensorMultiply(*cmat,fac,tmp1,tmp2,1.0);  // N_{abba}

         } // end if (a!=b)
       } // end loop b
     } // end loop a
   }



  /*
   **********************************************************
   * Step 8 * Save xi_S and the Drucker-Prager loading      *
   *        * function for the next time step.              *
   **********************************************************
   */
  xi_s_curr_->at(gp) = xi_S;
  druckerpragerloadingcurr_->at(gp) = drucker_prager_loading;

  return;

}  // Evaluate()

/*---------------------------------------------------------------------*
 | return names of visualization data (public)           hemmler 09/16 |
 *---------------------------------------------------------------------*/
void MAT::SuperElasticSMA::VisNames(std::map<std::string,int>& names)
{
  names["martensiticfraction"] = 1; // scalar
  names["druckerprager"] = 1; // scalar
}  // VisNames()

LINALG::Matrix<2,1> MAT::SuperElasticSMA::ComputeLocalNewtonResidual(LINALG::Matrix<2,1> lambda_s, double xi_s, loadingData loading, material mat_data) {
  LINALG::Matrix<2,1> R;
  if(mat_data.model == 1) {
    // Exponential model
    R(0) = std::pow( loading.F_AS_f , 2.0) * lambda_s(0) - loading.H_AS * mat_data.beta_AS * ( 1.0 - xi_s ) * ( loading.drucker_prager_AS - loading.drucker_prager_AS_last ); // Eqn. 50
    R(1) = std::pow( loading.F_SA_f , 2.0) * lambda_s(1) - loading.H_SA * mat_data.beta_SA * xi_s * ( loading.drucker_prager_SA - loading.drucker_prager_SA_last ); // Eqn. 51
  } else {
    // Linear model
    R(0) = loading.F_AS_f * lambda_s(0) + loading.H_AS * ( 1.0 - xi_s ) * ( loading.drucker_prager_AS - loading.drucker_prager_AS_last ); // Eqn. 52
    R(1) = loading.F_SA_f * lambda_s(1) - loading.H_SA * xi_s * ( loading.drucker_prager_SA - loading.drucker_prager_SA_last ); // Eqn. 53
  }
  return R;
}

LINALG::Matrix<2,2> MAT::SuperElasticSMA::ComputeLocalNewtonJacobian(LINALG::Matrix<2,1> lambda_s, double xi_s, loadingData loading, material mat_data) {
  LINALG::Matrix<2,2> d_R_d_lambda;
  if(mat_data.model == 1) {
    // exponential model
    d_R_d_lambda(0,0) = loading.H_AS * mat_data.beta_AS * (loading.drucker_prager_AS - loading.drucker_prager_AS_last) + mat_data.G1 * ( loading.H_AS * mat_data.beta_AS * ( 1.0 - xi_s ) - 2.0 * lambda_s(0) * loading.F_AS_f ) + std::pow(loading.F_AS_f,2.0); // Eqn. 88
    d_R_d_lambda(0,1)  = loading.H_AS * mat_data.beta_AS * (loading.drucker_prager_AS - loading.drucker_prager_AS_last) + mat_data.G1 * ( loading.H_AS * mat_data.beta_AS * ( 1.0 - xi_s ) - 2.0 * lambda_s(0) * loading.F_AS_f ); // Eqn. 89
    d_R_d_lambda(1,0)  = -loading.H_SA * mat_data.beta_SA * (loading.drucker_prager_SA - loading.drucker_prager_SA_last) +  mat_data.G1 * ( loading.H_SA * mat_data.beta_SA * xi_s - 2.0 *  lambda_s(1) * loading.F_SA_f); // Eqn. 90
    d_R_d_lambda(1,1)  = -loading.H_SA * mat_data.beta_SA * (loading.drucker_prager_SA - loading.drucker_prager_SA_last) +  mat_data.G1 * ( loading.H_SA * mat_data.beta_SA * xi_s - 2.0 *  lambda_s(1) * loading.F_SA_f) + std::pow(loading.F_SA_f,2.0); // Eqn. 91
  } else {
    // linear model
    d_R_d_lambda(0,0) = -mat_data.G1 * lambda_s(0) - loading.H_AS * ( ( loading.drucker_prager_AS - loading.drucker_prager_AS_last ) + (1.0 - xi_s) * mat_data.G1 ) + loading.F_AS_f; // Eqn. 82
    d_R_d_lambda(0,1) = -mat_data.G1 * lambda_s(0) - loading.H_AS * ( ( loading.drucker_prager_AS - loading.drucker_prager_AS_last ) + (1.0 - xi_s) * mat_data.G1 ); // Eqn. 83
    d_R_d_lambda(1,0) = -mat_data.G1 * lambda_s(1) - loading.H_SA * ( ( loading.drucker_prager_SA - loading.drucker_prager_SA_last ) - xi_s * mat_data.G1); // Eqn. 84
    d_R_d_lambda(1,1) = -mat_data.G1 * lambda_s(1) - loading.H_SA * ( ( loading.drucker_prager_SA - loading.drucker_prager_SA_last ) - xi_s * mat_data.G1) + loading.F_SA_f; // Eqn. 85
  }

  return d_R_d_lambda;
}

MAT::SuperElasticSMA::loadingData MAT::SuperElasticSMA::ComputeLocalNewtonLoading(double xi_S, double log_strain_vol, double log_strain_dev_norm, material mat_data) {
  loadingData loading;
  double kirchhoff_stress_volumetric_tmp = mat_data.bulk * ( log_strain_vol - 3.0 * mat_data.alpha * mat_data.epsilon_L * xi_S );
  double kirchhoff_stress_deviatoric_norm_tmp = 2.0 * mat_data.shear * ( log_strain_dev_norm - mat_data.epsilon_L * xi_S );

  loading.drucker_prager = kirchhoff_stress_deviatoric_norm_tmp + 3.0 * mat_data.alpha * kirchhoff_stress_volumetric_tmp;

  loading.drucker_prager_AS = loading.drucker_prager - mat_data.C_AS * mat_data.temperature;
  loading.drucker_prager_SA = loading.drucker_prager - mat_data.C_SA * mat_data.temperature;

  loading.F_AS_f = loading.drucker_prager_AS - mat_data.R_AS_f;
  loading.F_SA_f = loading.drucker_prager_SA - mat_data.R_SA_f;

  return loading;
}


/*---------------------------------------------------------------------*
 | return visualization data (public)                    hemmler 09/16 |
 *---------------------------------------------------------------------*/
bool MAT::SuperElasticSMA::VisData(
  const std::string& name,
  std::vector<double>& data,
  int numgp,
  int eleID
  )
{
  if (name == "martensiticfraction")
  {
    if ((int)data.size()!=1) dserror("size mismatch");
    double temp = 0.0;
    for (int iter=0;iter<numgp;iter++)
      temp += xi_s_last_->at(iter);
    data[0] = temp/numgp;
  } else if (name == "druckerprager")
  {
    if ((int)data.size()!=1) dserror("size mismatch");
    double F = 0.0;
    for (int iter=0;iter<numgp;iter++)
      if(F < druckerpragerloadinglast_->at(iter))
        F = druckerpragerloadinglast_->at(iter);
    data[0] = F;
  }
  return false;

}  // VisData()

/*---------------------------------------------------------------------*
 | return Kronecker delta (public)                    hemmler 09/16 |
 *---------------------------------------------------------------------*/
int MAT::SuperElasticSMA::Kron(int i, int j) {
  if(i==j) return 1;
  return 0;
} // kron()

/*---------------------------------------------------------------------*
 | return fourth order deviatoric identity tensor (public) hemmler 09/16 |
 *---------------------------------------------------------------------*/
double MAT::SuperElasticSMA::Idev(int i, int j, int k, int l) {
  return Kron(i,j) * ( Kron(i,k) * Kron(k,l)-1.0/3.0 * Kron(k,l) )  + 0.5 * ( 1 - Kron(i,j) ) * Kron(i,k) * Kron(j,l);
  //return kron(i,j) * ( kron(i,k) * kron(k,l)-1.0/3.0 * kron(k,l) )  + ( 1 - kron(i,j) ) * kron(i,k) * kron(j,l);
} // Idev()

/*---------------------------------------------------------------------*
 | transform voigt to tensor notation (public)           hemmler 09/16 |
 *---------------------------------------------------------------------*/
void MAT::SuperElasticSMA::VoigtIndex(int p, int *i, int *j) {
  switch (p) {
  case 0:
    *i = 0;
    *j = 0;
    break;
  case 1:
    *i = 1;
    *j = 1;
    break;
  case 2:
    *i = 2;
    *j = 2;
    break;
  case 3:
    *i = 0;
    *j = 1;
    break;
  case 4:
    *i = 1;
    *j = 2;
    break;
  case 5:
    *i = 0;
    *j = 2;
    break;
  }
} // voigtindex()

/*---------------------------------------------------------------------*
 | transform tensor to voigt notation (public)           hemmler 09/16 |
 *---------------------------------------------------------------------*/
void MAT::SuperElasticSMA::VoigtIndex(int i, int j, int *p) {
  if(i == 0 && j == 0)
    *p = 0;
  else if(i == 1 && j == 1)
    *p = 1;
  else if(i == 2 && j == 2)
    *p = 2;
  else if((i == 0 && j == 1) || (i == 1 && j == 0))
    *p = 3;
  else if((i == 1 && j == 2) || (i == 2 && j == 1))
    *p = 4;
  else if((i == 0 && j == 2) || (i == 2 && j == 0))
    *p = 5;
} // voigtindex()

/*---------------------------------------------------------------------*
 | Pullback of material tangent                          hemmler 09/16 |
 *---------------------------------------------------------------------*/
void MAT::SuperElasticSMA::Pullback4thTensorVoigt(const double jacobian, const LINALG::Matrix<3,3>& defgrdinv, const LINALG::Matrix<6,6>& cmatEul, LINALG::Matrix<6,6>* cmatLag) {
  int i;
  int j;
  int k;
  int l;
  for (int p=0;p<6;p++) {
    for (int q=0;q<6;q++) {
      int M;
      int N;
      VoigtIndex(p, &i, &j);
      VoigtIndex(q, &k, &l);

      for (int A=0;A<3;A++) {
        for (int B=0;B<3;B++) {
          for (int C=0;C<3;C++) {
            for (int D=0;D<3;D++) {

              VoigtIndex(A,B,&M);
              VoigtIndex(C,D,&N);

              (*cmatLag)(p,q) +=  1.0 / jacobian * defgrdinv(i,A) * defgrdinv(j,B) * defgrdinv(k,C) * defgrdinv(l,D) * cmatEul(M,N);


            }
          }
        }
      }
    }
  }

} // pullback4thTensorVoigt()

/*----------------------------------------------------------------------*/

