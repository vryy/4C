/*----------------------------------------------------------------------*/
/*! \file
\brief Contains the functions to establish local material law /
       stress-strain law for isotropic material following finite strain
       von-Mises plasticity with linear isotropic hardening
       and logarithmic hyperelastic material (i.e. linear relation
       between Kirchhoff-stress and logarithmic strain; aka Hencky
       material model).
       The principal stress based implementation follows
       Bonet and Wood: "Nonlinear continuum mechanics for finite element analysis."
       Cambridge University Press, Cambridge, 2008

       geometrically nonlinear, finite strains, rate-independent

       example input line:
       MAT 1 MAT_Struct_PlasticNlnLogNeoHooke YOUNG 206.9 NUE 0.29 DENS 0.0
         YIELD 0.45 ISOHARD 0.12924 SATHARDENING 0.715 HARDEXPO 16.93 VISC 1.0 RATE_DEPENDENCY 0.1

\level 2

\maintainer Jan Schnabel
*/
/*----------------------------------------------------------------------*
 | headers                                                              |
 *----------------------------------------------------------------------*/
#include "plasticnlnlogneohooke.H"
#include "matpar_bundle.H"
#include "material_service.H"

#include "../linalg/linalg_utils.H"
#include "../drt_lib/drt_globalproblem.H"


/*----------------------------------------------------------------------*
 | constructor (public)                                                 |
 *----------------------------------------------------------------------*/
MAT::PAR::PlasticNlnLogNeoHooke::PlasticNlnLogNeoHooke(Teuchos::RCP<MAT::PAR::Material> matdata)
    : Parameter(matdata),
      youngs_(matdata->GetDouble("YOUNG")),
      poissonratio_(matdata->GetDouble("NUE")),
      density_(matdata->GetDouble("DENS")),
      yield_(matdata->GetDouble("YIELD")),
      isohard_(matdata->GetDouble("ISOHARD")),
      infyield_(matdata->GetDouble("SATHARDENING")),
      hardexp_(matdata->GetDouble("HARDEXPO")),
      visc_(matdata->GetDouble("VISC")),
      rate_dependency_(matdata->GetDouble("RATE_DEPENDENCY"))
{
}


/*----------------------------------------------------------------------*
 | is called in Material::Factory from ReadMaterials()                  |
 *----------------------------------------------------------------------*/
Teuchos::RCP<MAT::Material> MAT::PAR::PlasticNlnLogNeoHooke::CreateMaterial()
{
  return Teuchos::rcp(new MAT::PlasticNlnLogNeoHooke(this));
}


MAT::PlasticNlnLogNeoHookeType MAT::PlasticNlnLogNeoHookeType::instance_;


/*----------------------------------------------------------------------*
 | is called in Material::Factory from ReadMaterials()                  |
 *----------------------------------------------------------------------*/
DRT::ParObject* MAT::PlasticNlnLogNeoHookeType::Create(const std::vector<char>& data)
{
  MAT::PlasticNlnLogNeoHooke* plasticneo = new MAT::PlasticNlnLogNeoHooke();
  plasticneo->Unpack(data);
  return plasticneo;
}


/*----------------------------------------------------------------------*
 | constructor (public)                                                 |
 *----------------------------------------------------------------------*/
MAT::PlasticNlnLogNeoHooke::PlasticNlnLogNeoHooke() : params_(NULL) {}


/*----------------------------------------------------------------------*
 | copy-constructor (public)                                            |
 *----------------------------------------------------------------------*/
MAT::PlasticNlnLogNeoHooke::PlasticNlnLogNeoHooke(MAT::PAR::PlasticNlnLogNeoHooke* params)
    : params_(params)
{
}


/*----------------------------------------------------------------------*
 | pack (public)                                                        |
 *----------------------------------------------------------------------*/
void MAT::PlasticNlnLogNeoHooke::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data, type);

  // matid
  int matid = -1;
  if (params_ != NULL) matid = params_->Id();  // in case we are in post-process mode
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
    histsize = accplstrainlast_->size();
  }
  AddtoPack(data, histsize);  // Length of history vector(s)
  for (int var = 0; var < histsize; ++var)
  {
    // insert history vectors to AddtoPack
    AddtoPack(data, accplstrainlast_->at(var));
    AddtoPack(data, invplrcglast_->at(var));
  }

  return;
}  // Pack()


/*----------------------------------------------------------------------*
 | unpack (public)                                                      |
 *----------------------------------------------------------------------*/
void MAT::PlasticNlnLogNeoHooke::Unpack(const std::vector<char>& data)
{
  isinit_ = true;
  std::vector<char>::size_type position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position, data, type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");

  // matid
  int matid;
  ExtractfromPack(position, data, matid);
  params_ = NULL;
  if (DRT::Problem::Instance()->Materials() != Teuchos::null)
    if (DRT::Problem::Instance()->Materials()->Num() != 0)
    {
      const int probinst = DRT::Problem::Instance()->Materials()->GetReadFromProblem();
      MAT::PAR::Parameter* mat =
          DRT::Problem::Instance(probinst)->Materials()->ParameterById(matid);
      if (mat->Type() == MaterialType())
        params_ = static_cast<MAT::PAR::PlasticNlnLogNeoHooke*>(mat);
      else
        dserror("Type of parameter material %d does not fit to calling type %d", mat->Type(),
            MaterialType());
    }

  // history data
  int histsize;
  ExtractfromPack(position, data, histsize);

  // if system is not yet initialized, the history vectors have to be intialized
  if (histsize == 0) isinit_ = false;

  invplrcglast_ = Teuchos::rcp(new std::vector<LINALG::Matrix<3, 3>>);
  invplrcgcurr_ = Teuchos::rcp(new std::vector<LINALG::Matrix<3, 3>>);

  accplstrainlast_ = Teuchos::rcp(new std::vector<double>);
  accplstraincurr_ = Teuchos::rcp(new std::vector<double>);

  for (int var = 0; var < histsize; ++var)
  {
    double tmp1 = 0.0;
    // scalar-valued vector of last converged state are unpacked
    ExtractfromPack(position, data, tmp1);
    accplstrainlast_->push_back(tmp1);

    LINALG::Matrix<3, 3> tmp(true);
    // vectors of last converged state are unpacked
    ExtractfromPack(position, data, tmp);
    invplrcglast_->push_back(tmp);

    // current vectors have to be initialized
    accplstraincurr_->push_back(tmp1);
    invplrcgcurr_->push_back(tmp);
  }

  if (position != data.size()) dserror("Mismatch in size of data %d <-> %d", data.size(), position);

  return;

}  // Unpack()


/*---------------------------------------------------------------------*
 | initialise / allocate internal variables (public)                   |
 *---------------------------------------------------------------------*/
void MAT::PlasticNlnLogNeoHooke::Setup(int numgp, DRT::INPUT::LineDefinition* linedef)
{
  invplrcglast_ = Teuchos::rcp(new std::vector<LINALG::Matrix<3, 3>>);
  invplrcgcurr_ = Teuchos::rcp(new std::vector<LINALG::Matrix<3, 3>>);

  accplstrainlast_ = Teuchos::rcp(new std::vector<double>);
  accplstraincurr_ = Teuchos::rcp(new std::vector<double>);

  invplrcglast_->resize(numgp);
  invplrcgcurr_->resize(numgp);

  accplstrainlast_->resize(numgp);
  accplstraincurr_->resize(numgp);

  LINALG::Matrix<3, 3> emptymat(true);
  for (int i = 0; i < 3; i++) emptymat(i, i) = 1.0;

  for (int i = 0; i < numgp; i++)
  {
    invplrcglast_->at(i) = emptymat;
    invplrcgcurr_->at(i) = emptymat;

    accplstrainlast_->at(i) = 0.0;
    accplstraincurr_->at(i) = 0.0;
  }

  isinit_ = true;
  return;

}  // Setup()


/*----------------------------------------------------------------------*
 | update internal variables                                            |
 *----------------------------------------------------------------------*/
void MAT::PlasticNlnLogNeoHooke::Update()
{
  // make current values at time step t_n+1 to values of last step t_n
  invplrcglast_ = invplrcgcurr_;
  accplstrainlast_ = accplstraincurr_;

  // empty vectors of current data
  invplrcgcurr_ = Teuchos::rcp(new std::vector<LINALG::Matrix<3, 3>>);
  accplstraincurr_ = Teuchos::rcp(new std::vector<double>);

  // get the size of the vector
  // (use the last vector, because it includes latest results, current is empty)
  const int histsize = invplrcglast_->size();
  invplrcgcurr_->resize(histsize);
  accplstraincurr_->resize(histsize);

  LINALG::Matrix<3, 3> emptymat(true);
  for (int i = 0; i < 3; i++) emptymat(i, i) = 1.0;

  for (int i = 0; i < histsize; i++)
  {
    invplrcgcurr_->at(i) = emptymat;
    accplstraincurr_->at(i) = 0.0;
  }
  return;
}  // Update()


/*----------------------------------------------------------------------*
 | calculate stress and constitutive tensor                             |
 *----------------------------------------------------------------------*/
void MAT::PlasticNlnLogNeoHooke::Evaluate(const LINALG::Matrix<3, 3>* defgrd,
    const LINALG::Matrix<6, 1>* glstrain, Teuchos::ParameterList& params,
    LINALG::Matrix<6, 1>* stress, LINALG::Matrix<6, 6>* cmat, const int eleGID)
{
  // extract the gauss points from the parameter list
  const int gp = params.get<int>("gp", -1);
  if (gp == -1) dserror("no Gauss point number provided in material");

  // elastic material data
  // get material parameters
  const double ym = params_->youngs_;                  // Young's modulus
  const double nu = params_->poissonratio_;            // Poisson's ratio
  const double G = ym / (2.0 * (1.0 + nu));            // shear modulus, mu=G
  const double kappa = ym / (3.0 * (1.0 - 2.0 * nu));  // bulk modulus

  // plastic material data
  const double yield = params_->yield_;          // initial yield stress
  const double isohard = params_->isohard_;      // linear isotropic hardening
  const double infyield = params_->infyield_;    // saturation yield stress
  const double hardexp = params_->hardexp_;      // nonlinear hardening exponent
  const double visc = params_->visc_;            // viscosity
  const double eps = params_->rate_dependency_;  // rate dependency

  const double detF = defgrd->Determinant();

  const double dt = params.get<double>("delta time");

  LINALG::Matrix<3, 3> invdefgrd(*defgrd);
  invdefgrd.Invert();

  // matrices for temporary stuff
  LINALG::Matrix<3, 3> tmp1;
  LINALG::Matrix<3, 3> tmp2;

  // 3x3 2nd-order identity matrix
  LINALG::Matrix<3, 3> id2(true);
  // 3x3 2nd-order deviatoric identity matrix in principal directions
  LINALG::Matrix<3, 3> Idev;
  for (int i = 0; i < 3; i++)
  {
    id2(i, i) = 1.0;
    for (int j = 0; j < 3; j++)
    {
      if (i == j)
        Idev(i, j) = 2.0 / 3.0;
      else
        Idev(i, j) = -1.0 / 3.0;
    }
  }

  // linear elasticity tensor in principal directions
  LINALG::Matrix<3, 3> D_ep_principal(Idev);
  D_ep_principal.Scale(2.0 * G);
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) D_ep_principal(i, j) += kappa;

  // plastic increment
  double gamma = 0.0;

  // ------------------------------------------------------- trial strain
  // elastic left Cauchy-Green deformation tensor (LCG) at trial state
  // Be_trial = b_{e,n+1}^{trial} = F_{n+1} C_{p,n}^{-1} F_{n+1}
  LINALG::Matrix<3, 3> Be_trial;
  tmp1.Multiply(*defgrd, invplrcglast_->at(gp));
  Be_trial.MultiplyNT(tmp1, *defgrd);

  // elastic LCG at final state
  LINALG::Matrix<3, 3> Be;

  // ***************************************************
  // Here we start the principal stress based algorithm
  // ***************************************************

  // convert to epetra format and solve eigenvalue problem
  // this matrix contains spatial eigen directions (the second
  // index corresponds to the eigenvalue)
  Epetra_SerialDenseMatrix n(3, 3);
  Epetra_SerialDenseVector lambda_trial_square(3);

  // convert Input Matrix in Epetra format
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) n(i, j) = Be_trial(i, j);

  // calculate eigenvectors and eigenvalues
  LINALG::SymmetricEigenProblem(n, lambda_trial_square);
  // eigenvectors are stored in n, i.e. original matrix inserted in method is
  // destroyed.
  // eigenvalues correspond to the square of the stretches

  // spatial principal direction n_alpha (in LINALG format)
  std::vector<LINALG::Matrix<3, 1>> spatial_principal_directions;
  spatial_principal_directions.resize(3);
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) (spatial_principal_directions.at(i))(j) = n(j, i);

  // material principal directions N_alpha
  // note, that for convenience in the later programming this is NOT
  // based on the correct pull-back N_alpha = lambda_alpha * F^-1 * n_alpha
  // Instead we just do N_alpha = F^{-1} * n_alpha
  std::vector<LINALG::Matrix<3, 1>> material_principal_directions;
  material_principal_directions.resize(3);
  for (int i = 0; i < 3; i++)
    material_principal_directions.at(i).Multiply(invdefgrd, spatial_principal_directions.at(i));

  // deviatoric Kirchhoff stress at trial state
  // tau^{trial} = 2 * G * log(sqrt(lambda^2)) - 2/3 * G * log(detF)
  LINALG::Matrix<3, 1> dev_KH_trial;
  for (int i = 0; i < 3; i++)
    dev_KH_trial(i) = G * std::log(lambda_trial_square(i)) - 2.0 / 3.0 * G * std::log(detF);

  // deviatoric Kirchhoff stress at final state
  // For now we store the trial states. In case of plastic
  // material reaction, we will update it later
  LINALG::Matrix<3, 1> dev_KH;

  // pressure (equal at trial state and final state) for the use with tau
  // tau = tau'_aa + p = tau'_aa + kappa * log(detF)
  // later used for S = F^{-1} . tau . F^{-T}, no scaling with detF required
  double pressure = kappa * std::log(detF);

  // ----------------------------------------------------- yield function
  // calculate von Mises equivalent stress at trial state
  double abs_dev_KH_trial = 0.0;
  for (int i = 0; i < 3; i++) abs_dev_KH_trial += dev_KH_trial(i) * dev_KH_trial(i);
  abs_dev_KH_trial = std::sqrt(abs_dev_KH_trial);

  double y_d = (yield + isohard * accplstrainlast_->at(gp) +
                   (infyield - yield) * (1. - exp(-hardexp * accplstrainlast_->at(gp)))) *
               pow(visc * gamma / dt + 1., eps);

  double f_trial = std::sqrt(3.0 / 2.0) * abs_dev_KH_trial - y_d;

  // switch elastic (f <= 0) and plastic (f > 0) state
  if (f_trial <= 0.0)  // ----------------------------------- elastic step
  {
    // trial state variables are final variables
    dev_KH.Update(dev_KH_trial);
    Be.Update(Be_trial);

    // plastic increment is zero
    gamma = 0.0;
  }

  else  // -------------------------------------------------- plastic step
  {
    // local Newton iteration for nonlinear isotropic hardening
    int maxiter = 10;
    double tol = 1.e-12;
    double res = 0.;
    double tan = 0.;
    int iter = 0;
    double dy_d_dgamma = 0.;

    y_d = (yield + isohard * (accplstrainlast_->at(gp) + gamma) +
              (infyield - yield) * (1. - exp(-hardexp * (accplstrainlast_->at(gp) + gamma)))) *
          pow(visc * gamma / dt + 1., eps);
    dy_d_dgamma =
        (isohard +
            (infyield - yield) * hardexp * exp(-hardexp * (accplstrainlast_->at(gp) + gamma))) *
            pow(visc * gamma / dt + 1., eps) +
        (yield + isohard * (accplstrainlast_->at(gp) + gamma) +
            (infyield - yield) * (1. - exp(-hardexp * (accplstrainlast_->at(gp) + gamma)))) *
            pow(visc * gamma / dt + 1., eps - 1.) * eps * visc / dt;

    for (iter = 0; iter < maxiter; ++iter)
    {
      res = sqrt(3.0 / 2.0) * abs_dev_KH_trial - 3. * G * gamma - y_d;

      if (abs(res) < tol) break;

      tan = -3. * G - dy_d_dgamma;
      gamma -= res / tan;

      y_d = (yield + isohard * (accplstrainlast_->at(gp) + gamma) +
                (infyield - yield) * (1. - exp(-hardexp * (accplstrainlast_->at(gp) + gamma)))) *
            pow(visc * gamma / dt + 1., eps);
      dy_d_dgamma =
          (isohard +
              (infyield - yield) * hardexp * exp(-hardexp * (accplstrainlast_->at(gp) + gamma))) *
              pow(visc * gamma / dt + 1., eps) +
          (yield + isohard * (accplstrainlast_->at(gp) + gamma) +
              (infyield - yield) * (1. - exp(-hardexp * (accplstrainlast_->at(gp) + gamma)))) *
              pow(visc * gamma / dt + 1., eps - 1.) * eps * visc / dt;
    }

    if (iter == maxiter)
    {
      // check, if errors are tolerated or should throw a dserror
      bool error_tol = false;
      if (params.isParameter("tolerate_errors")) error_tol = params.get<bool>("tolerate_errors");
      if (error_tol)
      {
        params.set<bool>("eval_error", true);
        return;
      }
      else
        dserror("Local Newton iteration unconverged in %i iterations. Residual: %d", iter, res);
    }
    // flow vector (7.54a)
    LINALG::Matrix<3, 1> flow_vector(dev_KH_trial);
    flow_vector.Scale(1.0 / (sqrt(2.0 / 3.0) * abs_dev_KH_trial));

    // stress return mapping (7.60)
    double fac_dev_KH = 1.0 - 2.0 * G * gamma / (std::sqrt(2.0 / 3.0) * abs_dev_KH_trial);
    dev_KH.Update(fac_dev_KH, dev_KH_trial, 0.0);

    // strain return mapping
    // b_{e,n+1} = sum_i^3 ( lambda_{n+1}^2 . n \otimes n )
    // lambda_{n+1} = lambda_{n+1}^{trial} / exp(gamma . flow_vector)
    Be.Clear();
    for (int i = 0; i < 3; i++)
    {
      tmp1.MultiplyNT(spatial_principal_directions.at(i), spatial_principal_directions.at(i));
      Be.Update((lambda_trial_square(i) / exp(2.0 * flow_vector(i) * gamma)), tmp1, 1.0);
    }

    // update tangent modulus
    double fac_D_ep_1 = -4.0 * G * G * gamma / (sqrt(2.0 / 3.0) * abs_dev_KH_trial);
    D_ep_principal.Update(fac_D_ep_1, Idev, 1.0);
    tmp1.MultiplyNT(flow_vector, flow_vector);
    double fac_D_ep_2 =
        4.0 * G * G * (sqrt(2.0 / 3.0) * gamma / abs_dev_KH_trial - 1.0 / (3.0 * G + dy_d_dgamma));
    D_ep_principal.Update(fac_D_ep_2, tmp1, 1.0);
  }

  // -------------------------------------------------- output PK2 stress
  // tau_ij = tau_ii * n_i \otimes n_i
  // or tau_ij = tau'_ii + kappa ln(J)
  // S_ij   = F^-1 tau_ij F^-T
  //        = tau_ij (F^-1 n_i) \otimes (F^-1 n_i)
  LINALG::Matrix<3, 3> PK2(true);
  for (int i = 0; i < 3; i++)
  {
    tmp1.MultiplyNT(material_principal_directions.at(i), material_principal_directions.at(i));
    PK2.Update((dev_KH(i) + pressure), tmp1, 1.0);
  }

  // output stress in Voigt notation
  (*stress)(0) = PK2(0, 0);
  (*stress)(1) = PK2(1, 1);
  (*stress)(2) = PK2(2, 2);
  (*stress)(3) = 0.5 * (PK2(0, 1) + PK2(1, 0));
  (*stress)(4) = 0.5 * (PK2(1, 2) + PK2(2, 1));
  (*stress)(5) = 0.5 * (PK2(2, 0) + PK2(0, 2));

  // ---------------------------------------------------- tangent modulus
  // express coefficents of tangent in Kirchhoff stresses
  cmat->Clear();
  for (int a = 0; a < 3; a++)
  {
    // - sum_1^3 (2 * tau N_aaaa)
    tmp1.MultiplyNT(
        material_principal_directions.at(a), material_principal_directions.at(a));  // N_{aa}
    ElastSymTensorMultiply(*cmat, -2.0 * (dev_KH(a) + pressure), tmp1, tmp1, 1.0);

    for (int b = 0; b < 3; b++)
    {
      // c_ab N_aabb
      // result of return mapping of deviatoric component c_ab
      tmp1.MultiplyNT(
          material_principal_directions.at(a), material_principal_directions.at(a));  // N_{aa}
      tmp2.MultiplyNT(
          material_principal_directions.at(b), material_principal_directions.at(b));  // N_{bb}
      ElastSymTensorMultiply(*cmat, D_ep_principal(a, b), tmp1, tmp2, 1.0);

      if (a != b)
      {
        double fac = 0.0;
        if (lambda_trial_square(a) != lambda_trial_square(b))
        {
          // (tau_aa * lambda_b^2 - tau_bb * lambda_a^2) / (lambda_a^2 - lambda_b^2)
          fac = ((dev_KH(a) + pressure) * lambda_trial_square(b) -
                    (dev_KH(b) + pressure) * lambda_trial_square(a)) /
                (lambda_trial_square(a) - lambda_trial_square(b));
        }     // end if lambda_a != lambda_b
        else  // lambda_a = lambda_b
        {
          // 1/2 [(d^2 Psi)/(d ln lambda_b * d ln lambda_b) - (d^2 Psi)/(d ln lambda_a * d ln
          // lambda_b)]
          // - tau_bb, cf. (6.91)
          fac = 0.5 * (D_ep_principal(b, b) - D_ep_principal(a, b)) - (dev_KH(b) + pressure);
        }  // end lambda_a = lambda_b
        tmp1.MultiplyNT(
            material_principal_directions.at(a), material_principal_directions.at(b));  // N_{ab}
        tmp2.MultiplyNT(
            material_principal_directions.at(b), material_principal_directions.at(a));  // N_{ba}
        ElastSymTensorMultiply(*cmat, fac, tmp1, tmp1, 1.0);                            // N_{abab}
        ElastSymTensorMultiply(*cmat, fac, tmp1, tmp2, 1.0);                            // N_{abba}

      }  // end if (a!=b)
    }    // end loop b
  }      // end loop a

  // --------------------------------------------- update plastic history
  // plastic inverse of right Cauchy-Green deformation tensor (RCG)
  tmp1.Multiply(invdefgrd, Be);
  invplrcgcurr_->at(gp).MultiplyNT(tmp1, invdefgrd);

  // accumulated plastic strain
  accplstraincurr_->at(gp) = accplstrainlast_->at(gp) + gamma;

  // Green-Lagrange plastic strains can be easily calculated, in contrast
  // Euler-Almansi requires special treatment, which is not yet considered in the
  // element formulation

  return;

}  // Evaluate()


/*---------------------------------------------------------------------*
 | return names of visualization data (public)                         |
 *---------------------------------------------------------------------*/
void MAT::PlasticNlnLogNeoHooke::VisNames(std::map<std::string, int>& names)
{
  std::string accumulatedstrain = "accumulatedstrain";
  names[accumulatedstrain] = 1;  // scalar

}  // VisNames()


/*---------------------------------------------------------------------*
 | return visualization data (public)                                  |
 *---------------------------------------------------------------------*/
bool MAT::PlasticNlnLogNeoHooke::VisData(
    const std::string& name, std::vector<double>& data, int numgp, int eleID)
{
  if (name == "accumulatedstrain")
  {
    if ((int)data.size() != 1) dserror("size mismatch");
    double temp = 0.0;
    for (int iter = 0; iter < numgp; iter++) temp += AccumulatedStrain(iter);
    data[0] = temp / numgp;
  }
  return false;

}  // VisData()


/*----------------------------------------------------------------------*/
