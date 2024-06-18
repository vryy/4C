/*----------------------------------------------------------------------*/
/*! \file
 \brief Active strain membrane material for gastric electromechanics

 The input line should read
 MAT 0 MAT_Membrane_ActiveStrain MATIDPASSIVE 1 SCALIDVOLTAGE 3 BETA1 1.2 BETA2 2.3 VOLTHRESH 0.5
 ALPHA1 1.0 ALPHA2 1.25

 \level 3


 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | headers                                         brandstaeter 05/2018 |
 *----------------------------------------------------------------------*/
#include "4C_mat_membrane_active_strain.hpp"

#include "4C_comm_utils_factory.hpp"
#include "4C_global_data.hpp"
#include "4C_linalg_utils_densematrix_inverse.hpp"
#include "4C_mat_membrane_elasthyper.hpp"
#include "4C_mat_par_bundle.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 | constructor                                     brandstaeter 05/2018 |
 *----------------------------------------------------------------------*/
Mat::PAR::MembraneActiveStrain::MembraneActiveStrain(const Core::Mat::PAR::Parameter::Data& matdata)
    : Parameter(matdata),
      matid_passive_(matdata.parameters.get<int>("MATIDPASSIVE")),
      scalid_voltage_(matdata.parameters.get<int>("SCALIDVOLTAGE")),
      density_(matdata.parameters.get<double>("DENS")),
      beta1_(matdata.parameters.get<double>("BETA1")),
      beta2_(matdata.parameters.get<double>("BETA2")),
      voltage_threshold_(matdata.parameters.get<double>("VOLTHRESH")),
      alpha1_(matdata.parameters.get<double>("ALPHA1")),
      alpha2_(matdata.parameters.get<double>("ALPHA2"))
{
  return;
}  // Mat::PAR::MembraneActiveStrain::MembraneActiveStrain

/*----------------------------------------------------------------------*
 |                                                 brandstaeter 05/2018 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<Core::Mat::Material> Mat::PAR::MembraneActiveStrain::create_material()
{
  return Teuchos::rcp(new Mat::MembraneActiveStrain(this));
}  // Mat::PAR::MembraneActiveStrain::create_material

Mat::MembraneActiveStrainType Mat::MembraneActiveStrainType::instance_;

Core::Communication::ParObject* Mat::MembraneActiveStrainType::Create(const std::vector<char>& data)
{
  Mat::MembraneActiveStrain* membrane_activestrain = new Mat::MembraneActiveStrain();
  membrane_activestrain->Unpack(data);

  return membrane_activestrain;
}  // Mat::Membrane_ActiveStrainType::Create

/*----------------------------------------------------------------------*
 |                                                 brandstaeter 05/2018 |
 *----------------------------------------------------------------------*/
Mat::MembraneActiveStrain::MembraneActiveStrain()
    : params_(nullptr),
      matpassive_(nullptr),
      voltage_(Teuchos::null),
      activation_(Teuchos::null),
      isinit_(false),
      fibervecs_(false)
{
  return;
}  // Mat::MembraneActiveStrain::MembraneActiveStrain()

/*----------------------------------------------------------------------*
 |                                                 brandstaeter 05/2018 |
 *----------------------------------------------------------------------*/
Mat::MembraneActiveStrain::MembraneActiveStrain(Mat::PAR::MembraneActiveStrain* params)
    : params_(params),
      matpassive_(nullptr),
      voltage_(Teuchos::null),
      activation_(Teuchos::null),
      isinit_(false),
      fibervecs_(false)
{
  return;
}  // Mat::MembraneActiveStrain::MembraneActiveStrain()

/*----------------------------------------------------------------------*
 |                                                 brandstaeter 05/2018 |
 *----------------------------------------------------------------------*/
void Mat::MembraneActiveStrain::Pack(Core::Communication::PackBuffer& data) const
{
  Core::Communication::PackBuffer::SizeMarker sm(data);

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  add_to_pack(data, type);

  // matid
  int matid = -1;
  if (params_ != nullptr) matid = params_->Id();  // in case we are in post-process mode
  add_to_pack(data, matid);

  // fiber vectors: Fiber1, Fiber2, Normal
  add_to_pack(data, fibervecs_);

  // data of passive elastic material
  if (matpassive_ != Teuchos::null)
  {
    matpassive_->Pack(data);
  }

  // pack internal variables
  int numgp;
  // if material is not initialized, i.e. start simulation, nothing to pack
  if (!isinit_)
  {
    numgp = 0;
  }
  else
  {
    // if material is initialized, i.e. restart of simulation, size equates number of gausspoints
    numgp = voltage_->size();
  }
  // Length of internal vector(s)
  add_to_pack(data, numgp);
  for (int gp = 0; gp < numgp; ++gp)
  {
    // insert internal vectors to add_to_pack
    add_to_pack(data, voltage_->at(gp));
    add_to_pack(data, activation_->at(gp));
  }

  return;
}  // Mat::MembraneActiveStrain::Pack()

/*----------------------------------------------------------------------*
 |                                                 brandstaeter 05/2018 |
 *----------------------------------------------------------------------*/
void Mat::MembraneActiveStrain::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  Core::Communication::ExtractAndAssertId(position, data, UniqueParObjectId());

  // matid and recover params_
  int matid = -1;
  extract_from_pack(position, data, matid);
  if (Global::Problem::Instance()->Materials() != Teuchos::null)
    if (Global::Problem::Instance()->Materials()->Num() != 0)
    {
      const int probinst = Global::Problem::Instance()->Materials()->GetReadFromProblem();
      Core::Mat::PAR::Parameter* mat =
          Global::Problem::Instance(probinst)->Materials()->ParameterById(matid);
      if (mat->Type() == MaterialType())
        params_ = static_cast<Mat::PAR::MembraneActiveStrain*>(mat);
      else
        FOUR_C_THROW("Type of parameter material %d does not fit to calling type %d", mat->Type(),
            MaterialType());
    }

  // fiber vectors: Fiber1, Fiber2, Normal
  extract_from_pack(position, data, fibervecs_);

  // unpack data of passive material
  std::vector<char> matpassive_data;
  extract_from_pack(position, data, matpassive_data);
  if (matpassive_data.size() > 0)
  {
    Core::Communication::ParObject* o =
        Core::Communication::Factory(matpassive_data);  // Unpack is done here
    Mat::So3Material* matpassive = dynamic_cast<Mat::So3Material*>(o);
    if (matpassive == nullptr) FOUR_C_THROW("failed to unpack passive material");

    matpassive_ = Teuchos::rcp(matpassive);
  }
  else
  {
    matpassive_ = Teuchos::null;
  }

  int numgp;
  extract_from_pack(position, data, numgp);
  isinit_ = true;
  if (numgp == 0)  // no internal data to unpack
  {
    isinit_ = false;
    if (position != data.size())
      FOUR_C_THROW("Mismatch in size of data %d <-> %d", data.size(), position);
    return;
  }

  // unpack internal variables
  voltage_ = Teuchos::rcp(new std::vector<double>(numgp));
  activation_ = Teuchos::rcp(new std::vector<double>(numgp));
  double voltage_gp;
  double activation_gp;
  for (int gp = 0; gp < numgp; ++gp)
  {
    extract_from_pack(position, data, voltage_gp);
    voltage_->at(gp) = voltage_gp;
    extract_from_pack(position, data, activation_gp);
    activation_->at(gp) = activation_gp;
  }
  return;
}  // Mat::MembraneActiveStrain::Unpack()

/*----------------------------------------------------------------------*
 |                                                 brandstaeter 05/2018 |
 *----------------------------------------------------------------------*/
void Mat::MembraneActiveStrain::setup(int numgp, Input::LineDefinition* linedef)
{
  // setup fibervectors
  setup_fiber_vectors(numgp, linedef);

  // setup of passive material
  matpassive_ =
      Teuchos::rcp_dynamic_cast<Mat::So3Material>(Mat::Factory(params_->matid_passive_), true);
  matpassive_->setup(numgp, linedef);

  // setup internal variables
  voltage_ = Teuchos::rcp(new std::vector<double>);
  voltage_->resize(numgp);

  activation_ = Teuchos::rcp(new std::vector<double>);
  activation_->resize(numgp);

  for (int gp = 0; gp < numgp; ++gp)
  {
    voltage_->at(gp) = 0.0;
    activation_->at(gp) = 0.0;
  }

  isinit_ = true;
  return;
}  // Mat::MembraneActiveStrain::setup()

/*----------------------------------------------------------------------*
 | active strain and hyperelastic stress response plus elasticity tensor|
 |                                                 brandstaeter 05/2018 |
 *----------------------------------------------------------------------*/
void Mat::MembraneActiveStrain::EvaluateMembrane(const Core::LinAlg::Matrix<3, 3>& defgrd,
    const Core::LinAlg::Matrix<3, 3>& cauchygreen, Teuchos::ParameterList& params,
    const Core::LinAlg::Matrix<3, 3>& Q_trafo, Core::LinAlg::Matrix<3, 1>& stress,
    Core::LinAlg::Matrix<3, 3>& cmat, const int gp, const int eleGID)
{
  // blank resulting quantities
  stress.Clear();
  cmat.Clear();

  // get pointer to vector containing the scalar states at the gauss points
  Teuchos::RCP<std::vector<std::vector<double>>> gpscalar =
      params.get<Teuchos::RCP<std::vector<std::vector<double>>>>("gp_scalar",
          Teuchos::rcp(new std::vector<std::vector<double>>(4, std::vector<double>(4, 0.0))));

  const unsigned int scalarid_voltage = params_->scalid_voltage_;

  if (scalarid_voltage >= gpscalar->at(0).size())
    FOUR_C_THROW("Mismatch in requested scalar and number of supplied scalars.");

  // voltage at current gp
  double gpvoltage = gpscalar->at(gp).at(scalarid_voltage);

  // save voltage for visualization
  voltage_->at(gp) = gpvoltage;

  // structural tensor in local coordinates
  std::vector<Core::LinAlg::Matrix<3, 3>> structural_tensors_loc;

  // loop over all fiber vectors
  Core::LinAlg::Matrix<3, 1> fibervector(true);
  Core::LinAlg::Matrix<3, 3> structuraltensor(true);
  for (unsigned int p = 0; p < 3; ++p)
  {
    fibervector.MultiplyTN(1.0, Q_trafo, fibervecs_[p], 0.0);
    structuraltensor.MultiplyNT(1.0, fibervector, fibervector, 0.0);
    structural_tensors_loc.push_back(structuraltensor);
  }

  //******************
  // ACTIVE deformation gradient in local coordinates
  //******************
  Core::LinAlg::Matrix<3, 3> defgrd_active_inv_loc(true);

  // set defgrd_active to identity tensor
  for (int i = 0; i < 3; i++) defgrd_active_inv_loc(i, i) = 1.0;

  // create full active def-grd
  double voltage_theshold = params_->voltage_threshold_;
  double beta1 = params_->beta1_;
  double beta2 = params_->beta2_;

  double gamma = 0;
  if (gpvoltage > voltage_theshold)
  {
    gamma = (1 - std::exp(-beta1 * (gpvoltage - voltage_theshold))) *
            (1 - std::exp(-beta2 * (gpvoltage - voltage_theshold)));
  }

  activation_->at(gp) = gamma;

  double gamma1 = gamma * params_->alpha1_;
  double gamma2 = gamma * params_->alpha2_;
  double gammaNormal =
      (1 - (1 - gamma1) * (1 - gamma2)) /
      ((1 - gamma1) * (1 - gamma2));  // compute gamma_n such that active material is incompressible

  defgrd_active_inv_loc.Update(gamma1 / (1.0 - gamma1), structural_tensors_loc.at(0), 1.0);
  defgrd_active_inv_loc.Update(gamma2 / (1.0 - gamma2), structural_tensors_loc.at(1), 1.0);
  defgrd_active_inv_loc.Update(
      -gammaNormal / (1.0 + gammaNormal), structural_tensors_loc.at(2), 1.0);

  //******************
  // PASSIVE cauchy green in local coordinates
  //******************
  Core::LinAlg::Matrix<3, 3> cauchygreen_passive_local(true);
  Core::LinAlg::Matrix<3, 3> defgrd_passive_local(true);
  defgrd_passive_local.MultiplyNN(1.0, defgrd, defgrd_active_inv_loc, 0.0);
  cauchygreen_passive_local.MultiplyTN(1.0, defgrd_passive_local, defgrd_passive_local, 0.0);

  // compute passive green lagrange strain
  Core::LinAlg::Matrix<3, 3> cmatpassive_loc(true);
  Core::LinAlg::Matrix<3, 1> S_passive_loc_voigt(true);
  Teuchos::rcp_dynamic_cast<Mat::MembraneElastHyper>(matpassive_, true)
      ->EvaluateMembrane(defgrd_passive_local, cauchygreen_passive_local, params, Q_trafo,
          S_passive_loc_voigt, cmatpassive_loc, gp, eleGID);

  //******************
  // FULL PART
  //******************
  Core::LinAlg::Matrix<2, 2> S_tot(true);
  Core::LinAlg::Matrix<2, 2> S_passive_loc(true);
  S_passive_loc(0, 0) = S_passive_loc_voigt(0);
  S_passive_loc(1, 1) = S_passive_loc_voigt(1);
  S_passive_loc(1, 0) = S_passive_loc_voigt(2);
  S_passive_loc(0, 1) = S_passive_loc_voigt(2);

  Core::LinAlg::Matrix<2, 2> defgrd_active_inv_loc_red(true);
  defgrd_active_inv_loc_red(0, 0) = defgrd_active_inv_loc(0, 0);
  defgrd_active_inv_loc_red(1, 0) = defgrd_active_inv_loc(1, 0);
  defgrd_active_inv_loc_red(0, 1) = defgrd_active_inv_loc(0, 1);
  defgrd_active_inv_loc_red(1, 1) = defgrd_active_inv_loc(1, 1);

  Core::LinAlg::Matrix<2, 2> temp2(true);
  temp2.MultiplyNT(1.0, S_passive_loc, defgrd_active_inv_loc_red, 0.0);
  S_tot.MultiplyNN(1.0, defgrd_active_inv_loc_red, temp2, 0.0);

  stress(0) = S_tot(0, 0);
  stress(1) = S_tot(1, 1);
  stress(2) = 0.5 * (S_tot(1, 0) + S_tot(0, 1));

  // pullback of the linearization
  pullback4th_tensor_voigt(defgrd_active_inv_loc_red, cmatpassive_loc, cmat);

  return;
}  // Mat::MembraneActiveStrain::Evaluate

/*----------------------------------------------------------------------*
 | Update internal variables                       brandstaeter 05/2018 |
 *----------------------------------------------------------------------*/
void Mat::MembraneActiveStrain::Update()
{
  matpassive_->Update();
}  // Mat::MembraneActiveStrain::Update

/*----------------------------------------------------------------------*
 | Reset internal variables                        brandstaeter 05/2018 |
 *----------------------------------------------------------------------*/
void Mat::MembraneActiveStrain::reset_step()
{
  matpassive_->reset_step();
}  // Mat::MembraneActiveStrain::reset_step

/*----------------------------------------------------------------------*
 |                                                 brandstaeter 05/2018 |
 *----------------------------------------------------------------------*/
void Mat::MembraneActiveStrain::VisNames(std::map<std::string, int>& names)
{
  matpassive_->VisNames(names);
  names["voltage"] = 1;     // scalar
  names["activation"] = 1;  // scalar
}  // Mat::MembraneActiveStrain::VisNames

/*----------------------------------------------------------------------*
 |                                                 brandstaeter 05/2018 |
 *----------------------------------------------------------------------*/
bool Mat::MembraneActiveStrain::VisData(
    const std::string& name, std::vector<double>& data, int numgp, int eleID)
{
  if (name == "voltage")
  {
    if ((int)data.size() != 1) FOUR_C_THROW("size mismatch");

    for (int gp = 0; gp < numgp; gp++) data[0] += voltage_->at(gp);

    data[0] = data[0] / numgp;
    return true;
  }
  else if (name == "activation")
  {
    if ((int)data.size() != 1) FOUR_C_THROW("size mismatch");

    for (int gp = 0; gp < numgp; gp++) data[0] += activation_->at(gp);

    data[0] = data[0] / numgp;
    return true;
  }

  return matpassive_->VisData(name, data, numgp, eleID);
}  // Mat::MembraneActiveStrain::VisData

/*----------------------------------------------------------------------*
 | setup fiber vectors                                                  |
 *----------------------------------------------------------------------*/
void Mat::MembraneActiveStrain::setup_fiber_vectors(int numgp, Input::LineDefinition* linedef)
{
  Core::LinAlg::Matrix<3, 1> dir;

  // CIR-AXI-RAD nomenclature
  if (linedef->has_named("RAD") and linedef->has_named("AXI") and linedef->has_named("CIR"))
  {
    // Axial direction
    read_dir(linedef, "AXI", dir);
    fibervecs_.push_back(dir);

    // Circumferential direction
    read_dir(linedef, "CIR", dir);
    fibervecs_.push_back(dir);

    // Radial direction
    read_dir(linedef, "RAD", dir);
    fibervecs_.push_back(dir);
  }
  // FIBER nomenclature
  else if (linedef->has_named("FIBER1") and linedef->has_named("FIBER2"))
  {
    for (int i = 1; i < 3; ++i)
    {
      std::ostringstream ss;
      ss << i;
      std::string fibername = "FIBER" + ss.str();  // FIBER Name
                                                   // FiberN direction
      read_dir(linedef, fibername, dir);
      fibervecs_.push_back(dir);
    }

    setup_normal_direction();
  }
  else
  {
    FOUR_C_THROW("Either use Fiber or CIR-AXI-RAD nomenclature to set fiber directions");
  }

  // Check orthonormal basis
  if (fibervecs_.size() != 3)
    FOUR_C_THROW(
        "Wrong number of fiber vectors. This material need three, it is %i", fibervecs_.size());

  double eps = 1e-12;
  if (std::abs(fibervecs_[0].Dot(fibervecs_[1])) > eps or
      std::abs(fibervecs_[1].Dot(fibervecs_[2])) > eps or
      std::abs(fibervecs_[0].Dot(fibervecs_[2])) > eps)
  {
    std::cout << std::endl;
    std::cout << "\tWARNING: fiber vectors do NOT build orthonormal basis!" << std::endl;
    std::cout << std::endl;
    FOUR_C_THROW(
        "Fiber vectors are not orthonormal: while this is not necessary in general, for now we "
        "limit ourselves to the orthonomal case!\n"
        "In particular the calculation of the inverse active deformation gradient depends on this "
        "assumption!");
  }

}  // Mat::MembraneActiveStrain::setup_fiber_vectors

/*----------------------------------------------------------------------*
 * Function which reads in the fiber direction
 *----------------------------------------------------------------------*/
void Mat::MembraneActiveStrain::read_dir(
    Input::LineDefinition* linedef, std::string specifier, Core::LinAlg::Matrix<3, 1>& dir)
{
  std::vector<double> fiber;
  linedef->extract_double_vector(specifier, fiber);
  double fnorm = 0.;
  // normalization
  for (int i = 0; i < 3; ++i)
  {
    fnorm += fiber[i] * fiber[i];
  }
  fnorm = sqrt(fnorm);

  // fill final normalized vector
  for (int i = 0; i < 3; ++i) dir(i) = fiber[i] / fnorm;

  return;
}  // Mat::MembraneActiveStrain::read_dir

void Mat::MembraneActiveStrain::setup_normal_direction()
{
  if (fibervecs_.size() != 2)
  {
    FOUR_C_THROW("Wrong number of fiber vectors to calculate a normal direction.");
  }

  Core::LinAlg::Matrix<3, 1> dir1 = fibervecs_[0];
  Core::LinAlg::Matrix<3, 1> dir2 = fibervecs_[1];
  Core::LinAlg::Matrix<3, 1> normaldir;

  normaldir(0) = dir1(1) * dir2(2) - dir1(2) * dir2(1);
  normaldir(1) = dir1(2) * dir2(0) - dir1(0) * dir2(2);
  normaldir(2) = dir1(0) * dir2(1) - dir1(1) * dir2(0);

  // normalization
  double norm = normaldir.Norm2();
  normaldir.Scale(1 / norm);

  fibervecs_.push_back(normaldir);
}  // Mat::MembraneActiveStrain::setup_normal_direction

/*---------------------------------------------------------------------*
 | Pullback of the tangent from intermediate to reference configuration|
 *---------------------------------------------------------------------*/
void Mat::MembraneActiveStrain::pullback4th_tensor_voigt(
    const Core::LinAlg::Matrix<2, 2>& defgrd_active_inv_red,
    const Core::LinAlg::Matrix<3, 3>& cmat_passive_intermediate,
    Core::LinAlg::Matrix<3, 3>& cmat_reference)
{
  int i;
  int j;
  int k;
  int l;
  for (int p = 0; p < 3; ++p)
  {
    for (int q = 0; q < 3; ++q)
    {
      int M;
      int N;
      tensor2x2_indices(p, &i, &j);
      tensor2x2_indices(q, &k, &l);

      for (int A = 0; A < 2; ++A)
      {
        for (int B = 0; B < 2; ++B)
        {
          for (int C = 0; C < 2; ++C)
          {
            for (int D = 0; D < 2; ++D)
            {
              voigt3_index(A, B, &M);
              voigt3_index(C, D, &N);

              cmat_reference(p, q) += defgrd_active_inv_red(i, A) * defgrd_active_inv_red(j, B) *
                                      defgrd_active_inv_red(k, C) * defgrd_active_inv_red(l, D) *
                                      cmat_passive_intermediate(M, N);
            }
          }
        }
      }
    }
  }

}  // Mat::MembraneActiveStrain::pullback4th_tensor_voigt

/*---------------------------------------------------------------------*
 | transform voigt to tensor notation                                  |
 *---------------------------------------------------------------------*/
void Mat::MembraneActiveStrain::tensor2x2_indices(int p, int* i, int* j)
{
  switch (p)
  {
    case 0:
      *i = 0;
      *j = 0;
      break;
    case 1:
      *i = 1;
      *j = 1;
      break;
    case 2:
      *i = 0;
      *j = 1;
      break;
  }
}  // Mat::MembraneActiveStrain::voigt3_index

/*---------------------------------------------------------------------*
 | transform tensor to voigt notation (public)                         |
 *---------------------------------------------------------------------*/
void Mat::MembraneActiveStrain::voigt3_index(int i, int j, int* p)
{
  if (i == 0 && j == 0)
    *p = 0;
  else if (i == 1 && j == 1)
    *p = 1;
  else if ((i == 0 && j == 1) || (i == 1 && j == 0))
    *p = 2;
}  // Mat::MembraneActiveStrain::voigt3_index

FOUR_C_NAMESPACE_CLOSE
