/*!----------------------------------------------------------------------
 \file membrane_active_strain.cpp

 \brief Active strain membrane material for gastric electromechanics

 The input line should read
 MAT 0 MAT_Membrane_ActiveStrain MATIDPASSIVE 1 SCALIDVOLTAGE 3 BETA1 1.2 BETA2 2.3 VOLTHRESH 0.5
 ALPHA1 1.0 ALPHA2 1.25

 \level 3

 \maintainer Sebastian Brandstaeter
             brandstaeter@lnm.mw.tum.de
             089 - 289-15276

 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | headers                                         brandstaeter 05/2018 |
 *----------------------------------------------------------------------*/
#include "membrane_active_strain.H"
#include "../drt_mat/matpar_bundle.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../linalg/linalg_utils.H"
#include "membrane_elasthyper.H"

/*----------------------------------------------------------------------*
 | constructor                                     brandstaeter 05/2018 |
 *----------------------------------------------------------------------*/
MAT::PAR::Membrane_ActiveStrain::Membrane_ActiveStrain(Teuchos::RCP<MAT::PAR::Material> matdata)
    : Parameter(matdata),
      matid_passive_(matdata->GetInt("MATIDPASSIVE")),
      scalid_voltage_(matdata->GetInt("SCALIDVOLTAGE")),
      density_(matdata->GetDouble("DENS")),
      beta1_(matdata->GetDouble("BETA1")),
      beta2_(matdata->GetDouble("BETA2")),
      voltage_threshold_(matdata->GetDouble("VOLTHRESH")),
      alpha1_(matdata->GetDouble("ALPHA1")),
      alpha2_(matdata->GetDouble("ALPHA2"))
{
  return;
}  // MAT::PAR::Membrane_ActiveStrain::Membrane_ActiveStrain

/*----------------------------------------------------------------------*
 |                                                 brandstaeter 05/2018 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<MAT::Material> MAT::PAR::Membrane_ActiveStrain::CreateMaterial()
{
  return Teuchos::rcp(new MAT::Membrane_ActiveStrain(this));
}  // MAT::PAR::Membrane_ActiveStrain::CreateMaterial

MAT::Membrane_ActiveStrainType MAT::Membrane_ActiveStrainType::instance_;

DRT::ParObject* MAT::Membrane_ActiveStrainType::Create(const std::vector<char>& data)
{
  MAT::Membrane_ActiveStrain* membrane_activestrain = new MAT::Membrane_ActiveStrain();
  membrane_activestrain->Unpack(data);

  return membrane_activestrain;
}  // MAT::Membrane_ActiveStrainType::Create

/*----------------------------------------------------------------------*
 |                                                 brandstaeter 05/2018 |
 *----------------------------------------------------------------------*/
MAT::Membrane_ActiveStrain::Membrane_ActiveStrain()
    : params_(NULL),
      matpassive_(0),
      voltage_(Teuchos::null),
      activation_(Teuchos::null),
      isinit_(false),
      fibervecs_(false)
{
  return;
}  // MAT::Membrane_ActiveStrain::Membrane_ActiveStrain()

/*----------------------------------------------------------------------*
 |                                                 brandstaeter 05/2018 |
 *----------------------------------------------------------------------*/
MAT::Membrane_ActiveStrain::Membrane_ActiveStrain(MAT::PAR::Membrane_ActiveStrain* params)
    : params_(params),
      matpassive_(0),
      voltage_(Teuchos::null),
      activation_(Teuchos::null),
      isinit_(false),
      fibervecs_(false)
{
  return;
}  // MAT::Membrane_ActiveStrain::Membrane_ActiveStrain()

/*----------------------------------------------------------------------*
 |                                                 brandstaeter 05/2018 |
 *----------------------------------------------------------------------*/
void MAT::Membrane_ActiveStrain::Pack(DRT::PackBuffer& data) const
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

  // fiber vectors: Fiber1, Fiber2, Normal
  AddtoPack(data, fibervecs_);

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
  AddtoPack(data, numgp);
  for (int gp = 0; gp < numgp; ++gp)
  {
    // insert internal vectors to AddtoPack
    AddtoPack(data, voltage_->at(gp));
    AddtoPack(data, activation_->at(gp));
  }

  return;
}  // MAT::Membrane_ActiveStrain::Pack()

/*----------------------------------------------------------------------*
 |                                                 brandstaeter 05/2018 |
 *----------------------------------------------------------------------*/
void MAT::Membrane_ActiveStrain::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position, data, type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");

  // matid and recover params_
  int matid = -1;
  ExtractfromPack(position, data, matid);
  if (DRT::Problem::Instance()->Materials() != Teuchos::null)
    if (DRT::Problem::Instance()->Materials()->Num() != 0)
    {
      const int probinst = DRT::Problem::Instance()->Materials()->GetReadFromProblem();
      MAT::PAR::Parameter* mat =
          DRT::Problem::Instance(probinst)->Materials()->ParameterById(matid);
      if (mat->Type() == MaterialType())
        params_ = static_cast<MAT::PAR::Membrane_ActiveStrain*>(mat);
      else
        dserror("Type of parameter material %d does not fit to calling type %d", mat->Type(),
            MaterialType());
    }

  // fiber vectors: Fiber1, Fiber2, Normal
  ExtractfromPack(position, data, fibervecs_);

  // unpack data of passive material
  std::vector<char> matpassive_data;
  ExtractfromPack(position, data, matpassive_data);
  if (matpassive_data.size() > 0)
  {
    DRT::ParObject* o = DRT::UTILS::Factory(matpassive_data);  // Unpack is done here
    MAT::So3Material* matpassive = dynamic_cast<MAT::So3Material*>(o);
    if (matpassive == NULL) dserror("failed to unpack passive material");

    matpassive_ = Teuchos::rcp(matpassive);
  }
  else
  {
    matpassive_ = Teuchos::null;
  }

  int numgp;
  ExtractfromPack(position, data, numgp);
  isinit_ = true;
  if (numgp == 0)  // no internal data to unpack
  {
    isinit_ = false;
    if (position != data.size())
      dserror("Mismatch in size of data %d <-> %d", data.size(), position);
    return;
  }

  // unpack internal variables
  voltage_ = Teuchos::rcp(new std::vector<double>(numgp));
  activation_ = Teuchos::rcp(new std::vector<double>(numgp));
  double voltage_gp;
  double activation_gp;
  for (int gp = 0; gp < numgp; ++gp)
  {
    ExtractfromPack(position, data, voltage_gp);
    voltage_->at(gp) = voltage_gp;
    ExtractfromPack(position, data, activation_gp);
    activation_->at(gp) = activation_gp;
  }
  return;
}  // MAT::Membrane_ActiveStrain::Unpack()

/*----------------------------------------------------------------------*
 |                                                 brandstaeter 05/2018 |
 *----------------------------------------------------------------------*/
void MAT::Membrane_ActiveStrain::Setup(int numgp, DRT::INPUT::LineDefinition* linedef)
{
  // setup fibervectors
  SetupFiberVectors(numgp, linedef);

  // setup of passive material
  matpassive_ = Teuchos::rcp_dynamic_cast<MAT::So3Material>(
      MAT::Material::Factory(params_->matid_passive_), true);
  matpassive_->Setup(numgp, linedef);

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
}  // MAT::Membrane_ActiveStrain::Setup()

/*----------------------------------------------------------------------*
 | active strain and hyperelastic stress response plus elasticity tensor|
 |                                                 brandstaeter 05/2018 |
 *----------------------------------------------------------------------*/
void MAT::Membrane_ActiveStrain::Evaluate(
    LINALG::Matrix<3, 3>& cauchygreen,  ///< right Cauchy-Green tensor
    Teuchos::ParameterList& params,     ///< Container for additional information
    LINALG::Matrix<3, 3>&
        Q_trafo,  ///< Trafo from local membrane orthonormal coordinates to global coordinates
    LINALG::Matrix<3, 1>* stress,  ///< 2nd Piola-Kirchhoff stresses in stress-like voigt notation
    LINALG::Matrix<3, 3>* cmat,    ///< Constitutive matrix
    const int eleGID               ///< Element GID
)
{
  // blank resulting quantities
  stress->Clear();
  cmat->Clear();

  // current gp
  const int gp = params.get<int>("gp");

  // get pointer to vector containing the scalar states at the gauss points
  Teuchos::RCP<std::vector<std::vector<double>>> gpscalar =
      params.get<Teuchos::RCP<std::vector<std::vector<double>>>>("gp_scalar",
          Teuchos::rcp(new std::vector<std::vector<double>>(4, std::vector<double>(4, 0.0))));

  const unsigned int scalarid_voltage = params_->scalid_voltage_;

  if (scalarid_voltage >= gpscalar->at(0).size())
    dserror("Mismatch in requested scalar and number of supplied scalars.");

  // voltage at current gp
  double gpvoltage = gpscalar->at(gp).at(scalarid_voltage);

  // save voltage for visualization
  voltage_->at(gp) = gpvoltage;

  // structural tensor in local coordinates
  std::vector<LINALG::Matrix<3, 3>> structural_tensors_loc;

  // loop over all fiber vectors
  LINALG::Matrix<3, 1> fibervector(true);
  LINALG::Matrix<3, 3> structuraltensor(true);
  for (unsigned int p = 0; p < 3; ++p)
  {
    fibervector.MultiplyTN(1.0, Q_trafo, fibervecs_[p], 0.0);
    structuraltensor.MultiplyNT(1.0, fibervector, fibervector, 0.0);
    structural_tensors_loc.push_back(structuraltensor);
  }

  //******************
  // ACTIVE deformation gradient in local coordinates
  //******************
  LINALG::Matrix<3, 3> defgrd_active_inv_loc(true);

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
  LINALG::Matrix<3, 3> cauchygreen_passive_local(true);
  LINALG::Matrix<3, 3> temp(true);
  temp.MultiplyTN(1.0, defgrd_active_inv_loc, cauchygreen, 0.0);
  cauchygreen_passive_local.MultiplyNN(1.0, temp, defgrd_active_inv_loc, 0.0);

  // compute passive green lagrange strain
  LINALG::Matrix<3, 3> cmatpassive_loc(true);
  LINALG::Matrix<3, 1> S_passive_loc_voigt(true);
  Teuchos::rcp_dynamic_cast<MAT::Membrane_ElastHyper>(matpassive_, true)
      ->Evaluate(cauchygreen_passive_local, params, Q_trafo, &S_passive_loc_voigt, &cmatpassive_loc,
          eleGID);

  //******************
  // FULL PART
  //******************
  LINALG::Matrix<2, 2> S_tot(true);
  LINALG::Matrix<2, 2> S_passive_loc(true);
  S_passive_loc(0, 0) = S_passive_loc_voigt(0);
  S_passive_loc(1, 1) = S_passive_loc_voigt(1);
  S_passive_loc(1, 0) = S_passive_loc_voigt(2);
  S_passive_loc(0, 1) = S_passive_loc_voigt(2);

  LINALG::Matrix<2, 2> defgrd_active_inv_loc_red(true);
  defgrd_active_inv_loc_red(0, 0) = defgrd_active_inv_loc(0, 0);
  defgrd_active_inv_loc_red(1, 0) = defgrd_active_inv_loc(1, 0);
  defgrd_active_inv_loc_red(0, 1) = defgrd_active_inv_loc(0, 1);
  defgrd_active_inv_loc_red(1, 1) = defgrd_active_inv_loc(1, 1);

  LINALG::Matrix<2, 2> temp2(true);
  temp2.MultiplyNT(1.0, S_passive_loc, defgrd_active_inv_loc_red, 0.0);
  S_tot.MultiplyNN(1.0, defgrd_active_inv_loc_red, temp2, 0.0);

  (*stress)(0) = S_tot(0, 0);
  (*stress)(1) = S_tot(1, 1);
  (*stress)(2) = 0.5 * (S_tot(1, 0) + S_tot(0, 1));

  // pullback of the linearization
  Pullback4thTensorVoigt(defgrd_active_inv_loc_red, cmatpassive_loc, cmat);

  return;
}  // MAT::Membrane_ActiveStrain::Evaluate

/*----------------------------------------------------------------------*
 | Update internal variables                       brandstaeter 05/2018 |
 *----------------------------------------------------------------------*/
void MAT::Membrane_ActiveStrain::Update()
{
  matpassive_->Update();
}  // MAT::Membrane_ActiveStrain::Update

/*----------------------------------------------------------------------*
 | Reset internal variables                        brandstaeter 05/2018 |
 *----------------------------------------------------------------------*/
void MAT::Membrane_ActiveStrain::ResetStep()
{
  matpassive_->ResetStep();
}  // MAT::Membrane_ActiveStrain::ResetStep

/*----------------------------------------------------------------------*
 |  ResetAll                                       brandstaeter 05/2018 |
 *----------------------------------------------------------------------*/
void MAT::Membrane_ActiveStrain::ResetAll(const int numgp)
{
  matpassive_->ResetAll(numgp);

  for (int gp = 0; gp < numgp; ++gp)
  {
    voltage_->at(gp) = 0.0;
    activation_->at(gp) = 0.0;
  }
  return;

}  // MAT::Membrane_ActiveStrain::ResetAll()

/*----------------------------------------------------------------------*
 |                                                 brandstaeter 05/2018 |
 *----------------------------------------------------------------------*/
void MAT::Membrane_ActiveStrain::VisNames(std::map<std::string, int>& names)
{
  matpassive_->VisNames(names);
  names["voltage"] = 1;     // scalar
  names["activation"] = 1;  // scalar
}  // MAT::Membrane_ActiveStrain::VisNames

/*----------------------------------------------------------------------*
 |                                                 brandstaeter 05/2018 |
 *----------------------------------------------------------------------*/
bool MAT::Membrane_ActiveStrain::VisData(
    const std::string& name, std::vector<double>& data, int numgp, int eleID)
{
  if (name == "voltage")
  {
    if ((int)data.size() != 1) dserror("size mismatch");

    for (int gp = 0; gp < numgp; gp++) data[0] += voltage_->at(gp);

    data[0] = data[0] / numgp;
    return true;
  }
  else if (name == "activation")
  {
    if ((int)data.size() != 1) dserror("size mismatch");

    for (int gp = 0; gp < numgp; gp++) data[0] += activation_->at(gp);

    data[0] = data[0] / numgp;
    return true;
  }

  return matpassive_->VisData(name, data, numgp, eleID);
}  // MAT::Membrane_ActiveStrain::VisData

/*----------------------------------------------------------------------*
 | setup fiber vectors                                                  |
 *----------------------------------------------------------------------*/
void MAT::Membrane_ActiveStrain::SetupFiberVectors(int numgp, DRT::INPUT::LineDefinition* linedef)
{
  LINALG::Matrix<3, 1> dir;

  // CIR-AXI-RAD nomenclature
  if (linedef->HaveNamed("RAD") and linedef->HaveNamed("AXI") and linedef->HaveNamed("CIR"))
  {
    // Axial direction
    ReadDir(linedef, "AXI", dir);
    fibervecs_.push_back(dir);

    // Circumferential direction
    ReadDir(linedef, "CIR", dir);
    fibervecs_.push_back(dir);

    // Radial direction
    ReadDir(linedef, "RAD", dir);
    fibervecs_.push_back(dir);
  }
  // FIBER nomenclature
  else if (linedef->HaveNamed("FIBER1") and linedef->HaveNamed("FIBER2"))
  {
    for (int i = 1; i < 3; ++i)
    {
      std::ostringstream ss;
      ss << i;
      std::string fibername = "FIBER" + ss.str();  // FIBER Name
                                                   // FiberN direction
      ReadDir(linedef, fibername, dir);
      fibervecs_.push_back(dir);
    }

    SetupNormalDirection();
  }
  else
  {
    dserror("Either use Fiber or CIR-AXI-RAD nomenclature to set fiber directions");
  }

  // Check orthonormal basis
  if (fibervecs_.size() != 3)
    dserror("Wrong number of fiber vectors. This material need three, it is %i", fibervecs_.size());

  double eps = 1e-12;
  if (std::abs(fibervecs_[0].Dot(fibervecs_[1])) > eps or
      std::abs(fibervecs_[1].Dot(fibervecs_[2])) > eps or
      std::abs(fibervecs_[0].Dot(fibervecs_[2])) > eps)
  {
    std::cout << std::endl;
    std::cout << "\tWARNING: fiber vectors do NOT build orthonormal basis!" << std::endl;
    std::cout << std::endl;
    dserror(
        "Fiber vectors are not orthonormal: while this is not necessary in general, for now we "
        "limit ourselves to the orthonomal case!\n"
        "In particular the calculation of the inverse active deformation gradient depends on this "
        "assumption!");
  }

}  // MAT::Membrane_ActiveStrain::SetupFiberVectors

/*----------------------------------------------------------------------*
 * Function which reads in the fiber direction
 *----------------------------------------------------------------------*/
void MAT::Membrane_ActiveStrain::ReadDir(
    DRT::INPUT::LineDefinition* linedef, std::string specifier, LINALG::Matrix<3, 1>& dir)
{
  std::vector<double> fiber;
  linedef->ExtractDoubleVector(specifier, fiber);
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
}  // MAT::Membrane_ActiveStrain::ReadDir

void MAT::Membrane_ActiveStrain::SetupNormalDirection()
{
  if (fibervecs_.size() != 2)
  {
    dserror("Wrong number of fiber vectors to calculate a normal direction.");
  }

  LINALG::Matrix<3, 1> dir1 = fibervecs_[0];
  LINALG::Matrix<3, 1> dir2 = fibervecs_[1];
  LINALG::Matrix<3, 1> normaldir;

  normaldir(0) = dir1(1) * dir2(2) - dir1(2) * dir2(1);
  normaldir(1) = dir1(2) * dir2(0) - dir1(0) * dir2(2);
  normaldir(2) = dir1(0) * dir2(1) - dir1(1) * dir2(0);

  // normalization
  double norm = normaldir.Norm2();
  normaldir.Scale(1 / norm);

  fibervecs_.push_back(normaldir);
}  // MAT::Membrane_ActiveStrain::SetupNormalDirection

/*---------------------------------------------------------------------*
 | Pullback of the tangent from intermediate to reference configuration|
 *---------------------------------------------------------------------*/
void MAT::Membrane_ActiveStrain::Pullback4thTensorVoigt(
    const LINALG::Matrix<2, 2>& defgrd_active_inv_red,
    const LINALG::Matrix<3, 3>& cmat_passive_intermediate, LINALG::Matrix<3, 3>* cmat_reference)
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
      TensorIndices(p, &i, &j);
      TensorIndices(q, &k, &l);

      for (int A = 0; A < 2; ++A)
      {
        for (int B = 0; B < 2; ++B)
        {
          for (int C = 0; C < 2; ++C)
          {
            for (int D = 0; D < 2; ++D)
            {
              VoigtIndex(A, B, &M);
              VoigtIndex(C, D, &N);

              (*cmat_reference)(p, q) += defgrd_active_inv_red(i, A) * defgrd_active_inv_red(j, B) *
                                         defgrd_active_inv_red(k, C) * defgrd_active_inv_red(l, D) *
                                         cmat_passive_intermediate(M, N);
            }
          }
        }
      }
    }
  }

}  // MAT::Membrane_ActiveStrain::Pullback4thTensorVoigt

/*---------------------------------------------------------------------*
 | transform voigt to tensor notation                                  |
 *---------------------------------------------------------------------*/
void MAT::Membrane_ActiveStrain::TensorIndices(int p, int* i, int* j)
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
}  // MAT::Membrane_ActiveStrain::VoigtIndex

/*---------------------------------------------------------------------*
 | transform tensor to voigt notation (public)                         |
 *---------------------------------------------------------------------*/
void MAT::Membrane_ActiveStrain::VoigtIndex(int i, int j, int* p)
{
  if (i == 0 && j == 0)
    *p = 0;
  else if (i == 1 && j == 1)
    *p = 1;
  else if ((i == 0 && j == 1) || (i == 1 && j == 0))
    *p = 2;
}  // MAT::Membrane_ActiveStrain::VoigtIndex
