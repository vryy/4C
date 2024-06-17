/*----------------------------------------------------------------------*/
/*! \file
 \brief wrapper for structure material of porous media


\level 2
 *-----------------------------------------------------------------------*/

#include "4C_mat_structporo.hpp"

#include "4C_comm_utils_factory.hpp"  // for function Factory in Unpack
#include "4C_global_data.hpp"
#include "4C_mat_par_bundle.hpp"
#include "4C_mat_poro_law.hpp"
#include "4C_mat_so3_material.hpp"

#include <vector>

FOUR_C_NAMESPACE_OPEN

Mat::PAR::StructPoro::StructPoro(const Core::Mat::PAR::Parameter::Data& matdata)
    : Parameter(matdata),
      matid_(matdata.parameters.Get<int>("MATID")),
      poro_law_ID_(matdata.parameters.Get<int>("POROLAWID")),
      init_porosity_(matdata.parameters.Get<double>("INITPOROSITY"))
{
  // retrieve problem instance to read from
  const int probinst = Global::Problem::Instance()->Materials()->GetReadFromProblem();

  // for the sake of safety
  if (Global::Problem::Instance(probinst)->Materials() == Teuchos::null)
    FOUR_C_THROW("List of materials cannot be accessed in the global problem instance.");
  // yet another safety check
  if (Global::Problem::Instance(probinst)->Materials()->Num() == 0)
    FOUR_C_THROW("List of materials in the global problem instance is empty.");

  // retrieve validated input line of material ID in question
  auto* curmat = Global::Problem::Instance(probinst)->Materials()->ParameterById(poro_law_ID_);

  switch (curmat->Type())
  {
    case Core::Materials::m_poro_law_linear:
    {
      poro_law_ = static_cast<Mat::PAR::PoroLaw*>(curmat);
      break;
    }
    case Core::Materials::m_poro_law_constant:
    {
      poro_law_ = static_cast<Mat::PAR::PoroLaw*>(curmat);
      break;
    }
    case Core::Materials::m_poro_law_logNeoHooke_Penalty:
    {
      poro_law_ = static_cast<Mat::PAR::PoroLaw*>(curmat);
      break;
    }
    case Core::Materials::m_poro_law_incompr_skeleton:
    {
      poro_law_ = static_cast<Mat::PAR::PoroLawIncompSkeleton*>(curmat);
      break;
    }
    case Core::Materials::m_poro_law_linear_biot:
    {
      poro_law_ = static_cast<Mat::PAR::PoroLawLinBiot*>(curmat);
      break;
    }
    case Core::Materials::m_poro_law_density_dependent:
    {
      poro_law_ = static_cast<Mat::PAR::PoroLawDensityDependent*>(curmat);
      break;
    }
    default:
      FOUR_C_THROW("invalid material for porosity law %d", curmat->Type());
      break;
  }
}

Teuchos::RCP<Core::Mat::Material> Mat::PAR::StructPoro::create_material()
{
  return Teuchos::rcp(new Mat::StructPoro(this));
}

Mat::StructPoroType Mat::StructPoroType::instance_;

Core::Communication::ParObject* Mat::StructPoroType::Create(const std::vector<char>& data)
{
  auto* struct_poro = new Mat::StructPoro();
  struct_poro->Unpack(data);
  return struct_poro;
}

Mat::StructPoro::StructPoro()
    : params_(nullptr),
      mat_(Teuchos::null),
      porosity_(Teuchos::null),
      surf_porosity_(Teuchos::null),
      is_initialized_(false)
{
}

Mat::StructPoro::StructPoro(Mat::PAR::StructPoro* params)
    : params_(params),
      porosity_(Teuchos::null),
      surf_porosity_(Teuchos::null),
      is_initialized_(false)
{
  mat_ = Teuchos::rcp_dynamic_cast<Mat::So3Material>(Mat::Factory(params_->matid_));
  if (mat_ == Teuchos::null)
    FOUR_C_THROW("Mat::StructPoro: underlying material should be of type Mat::So3Material");
}

void Mat::StructPoro::poro_setup(int numgp, Input::LineDefinition* linedef)
{
  porosity_ = Teuchos::rcp(new std::vector<double>(numgp, params_->init_porosity_));
  surf_porosity_ = Teuchos::rcp(new std::map<int, std::vector<double>>);

  is_initialized_ = true;
}

inline Core::Materials::MaterialType Mat::StructPoro::PoroLawType() const
{
  return params_->poro_law_->Type();
}

double Mat::StructPoro::InvBulkModulus() const { return params_->poro_law_->InvBulkModulus(); }

double Mat::StructPoro::Density() const
{
  if (params_->init_porosity_ == 1.0)
    return mat_->Density();
  else
    return ((1.0 - params_->init_porosity_) * mat_->Density());
}

double Mat::StructPoro::DensitySolidPhase() const { return mat_->Density(); }

void Mat::StructPoro::Pack(Core::Communication::PackBuffer& data) const
{
  if (not is_initialized_) FOUR_C_THROW("poro material not initialized. Not a poro element?");

  Core::Communication::PackBuffer::SizeMarker sm(data);

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  add_to_pack(data, type);

  // matid
  int matid = -1;
  if (params_ != nullptr) matid = params_->Id();  // in case we are in post-process mode
  add_to_pack(data, matid);

  // porosity_
  int size = static_cast<int>(porosity_->size());
  add_to_pack(data, size);
  for (int i = 0; i < size; ++i)
  {
    add_to_pack(data, (*porosity_)[i]);
  }

  // surf_porosity_
  size = static_cast<int>(surf_porosity_->size());
  add_to_pack(data, size);
  // iterator
  std::map<int, std::vector<double>>::const_iterator iter;
  for (iter = surf_porosity_->begin(); iter != surf_porosity_->end(); ++iter)
  {
    add_to_pack(data, iter->first);
    add_to_pack(data, iter->second);
  }

  // Pack data of underlying material
  if (mat_ != Teuchos::null) mat_->Pack(data);
}

void Mat::StructPoro::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  Core::Communication::ExtractAndAssertId(position, data, UniqueParObjectId());

  // matid
  int matid;
  extract_from_pack(position, data, matid);
  params_ = nullptr;
  if (Global::Problem::Instance()->Materials() != Teuchos::null)
  {
    if (Global::Problem::Instance()->Materials()->Num() != 0)
    {
      const int probinst = Global::Problem::Instance()->Materials()->GetReadFromProblem();
      Core::Mat::PAR::Parameter* mat =
          Global::Problem::Instance(probinst)->Materials()->ParameterById(matid);
      if (mat->Type() == MaterialType())
        params_ = static_cast<Mat::PAR::StructPoro*>(mat);
      else
        FOUR_C_THROW("Type of parameter material %d does not fit to calling type %d", mat->Type(),
            MaterialType());
    }
  }

  // porosity_
  int size = 0;
  extract_from_pack(position, data, size);
  porosity_ = Teuchos::rcp(new std::vector<double>);
  double tmp = 0.0;
  for (int i = 0; i < size; ++i)
  {
    extract_from_pack(position, data, tmp);
    porosity_->push_back(tmp);
  }

  // surface porosity (i think it is not necessary to pack/unpack this...)
  extract_from_pack(position, data, size);
  surf_porosity_ = Teuchos::rcp(new std::map<int, std::vector<double>>);
  for (int i = 0; i < size; i++)
  {
    int dof;
    std::vector<double> value;
    extract_from_pack(position, data, dof);
    extract_from_pack(position, data, value);

    // add to map
    surf_porosity_->insert(std::pair<int, std::vector<double>>(dof, value));
  }

  // Unpack data of sub material (these lines are copied from element.cpp)
  std::vector<char> datamat;
  extract_from_pack(position, data, datamat);
  if (datamat.size() > 0)
  {
    Core::Communication::ParObject* o =
        Core::Communication::Factory(datamat);  // Unpack is done here
    auto* mat = dynamic_cast<Mat::So3Material*>(o);
    if (mat == nullptr) FOUR_C_THROW("failed to unpack elastic material");
    mat_ = Teuchos::rcp(mat);
  }
  else
    mat_ = Teuchos::null;

  is_initialized_ = true;
}

void Mat::StructPoro::compute_porosity(const double& refporosity, const double& press,
    const double& J, const int& gp, double& porosity, double* dphi_dp, double* dphi_dJ,
    double* dphi_dJdp, double* dphi_dJJ, double* dphi_dpp, double* dphi_dphiref, bool save)
{
  if (refporosity == 1.0)
  {
    // this is pure fluid. The porosity does not change

    porosity = refporosity;
    if (dphi_dp) *dphi_dp = 0.0;
    if (dphi_dJ) *dphi_dJ = 0.0;
    if (dphi_dJdp) *dphi_dJdp = 0.0;
    if (dphi_dJJ) *dphi_dJJ = 0.0;
    if (dphi_dpp) *dphi_dpp = 0.0;
    if (dphi_dphiref) *dphi_dphiref = 0.0;
    return;
  }

  params_->poro_law_->compute_porosity(refporosity, press, J, gp, porosity, dphi_dp, dphi_dJ,
      dphi_dJdp, dphi_dJJ, dphi_dpp, dphi_dphiref);

  // save porosity
  if (save) porosity_->at(gp) = porosity;
}

void Mat::StructPoro::compute_porosity(Teuchos::ParameterList& params, double press, double J,
    int gp, double& porosity, double* dphi_dp, double* dphi_dJ, double* dphi_dJdp, double* dphi_dJJ,
    double* dphi_dpp, bool save)
{
  compute_porosity(
      params_
          ->init_porosity_,  // reference porosity equals initial porosity for non reactive material
      press, J, gp, porosity, dphi_dp, dphi_dJ, dphi_dJdp, dphi_dJJ, dphi_dpp,
      nullptr,  // reference porosity is constant (non reactive) -> derivative not needed
      save);
}

void Mat::StructPoro::compute_porosity(
    Teuchos::ParameterList& params, double press, double J, int gp, double& porosity, bool save)
{
  compute_porosity(
      params, press, J, gp, porosity, nullptr, nullptr, nullptr, nullptr, nullptr, save);
}

void Mat::StructPoro::ComputeSurfPorosity(Teuchos::ParameterList& params, double press, double J,
    const int surfnum, int gp, double& porosity, double* dphi_dp, double* dphi_dJ,
    double* dphi_dJdp, double* dphi_dJJ, double* dphi_dpp, bool save)
{
  compute_porosity(
      params, press, J, gp, porosity, dphi_dp, dphi_dJ, dphi_dJdp, dphi_dJJ, dphi_dpp, save);

  if (save)
  {
    if (gp == 0)  // it's a new iteration, so old values are not needed any more
      ((*surf_porosity_)[surfnum]).clear();

    ((*surf_porosity_)[surfnum]).push_back(porosity);
  }
}

void Mat::StructPoro::ComputeSurfPorosity(Teuchos::ParameterList& params, double press, double J,
    const int surfnum, int gp, double& porosity, bool save)
{
  ComputeSurfPorosity(
      params, press, J, surfnum, gp, porosity, nullptr, nullptr, nullptr, nullptr, nullptr, save);
}


double Mat::StructPoro::PorosityAv() const
{
  double porosityav = 0.0;

  std::vector<double>::const_iterator m;
  for (m = porosity_->begin(); m != porosity_->end(); ++m)
  {
    porosityav += *m;
  }
  porosityav = porosityav / (porosity_->size());

  return porosityav;
}

void Mat::StructPoro::CouplStress(const Core::LinAlg::Matrix<3, 3>& defgrd, const double& press,
    Core::LinAlg::Matrix<6, 1>& couplstress) const
{
  const double J = defgrd.Determinant();

  // Right Cauchy-Green tensor = F^T * F
  Core::LinAlg::Matrix<3, 3> cauchygreen;
  cauchygreen.MultiplyTN(defgrd, defgrd);

  // inverse Right Cauchy-Green tensor
  Core::LinAlg::Matrix<3, 3> C_inv;
  C_inv.Invert(cauchygreen);

  // inverse Right Cauchy-Green tensor as vector
  Core::LinAlg::Matrix<6, 1> C_inv_vec;
  for (int i = 0, k = 0; i < 3; i++)
    for (int j = 0; j < 3 - i; j++, k++) C_inv_vec(k) = C_inv(i + j, j);

  for (int i = 0; i < 6; i++) couplstress(i) = -1.0 * J * press * C_inv_vec(i);
}

void Mat::StructPoro::CouplStress(const Core::LinAlg::Matrix<2, 2>& defgrd, const double& press,
    Core::LinAlg::Matrix<4, 1>& couplstress) const
{
  const double J = defgrd.Determinant();

  // Right Cauchy-Green tensor = F^T * F
  Core::LinAlg::Matrix<2, 2> cauchygreen;
  cauchygreen.MultiplyTN(defgrd, defgrd);

  // inverse Right Cauchy-Green tensor
  Core::LinAlg::Matrix<2, 2> C_inv;
  C_inv.Invert(cauchygreen);

  // inverse Right Cauchy-Green tensor as vector
  Core::LinAlg::Matrix<3, 1> C_inv_vec;
  for (int i = 0, k = 0; i < 2; i++)
    for (int j = 0; j < 2 - i; j++, k++) C_inv_vec(k) = C_inv(i + j, j);

  couplstress(0) = -1.0 * J * press * C_inv_vec(0);
  couplstress(1) = -1.0 * J * press * C_inv_vec(1);
  couplstress(2) =
      0.0;  // this is needed to be compatible with the implementation of the wall element
  couplstress(3) = -1.0 * J * press * C_inv_vec(2);
}

void Mat::StructPoro::constitutive_derivatives(Teuchos::ParameterList& params, double press,
    double J, double porosity, double* dW_dp, double* dW_dphi, double* dW_dJ, double* dW_dphiref,
    double* W)
{
  if (porosity == 0.0) FOUR_C_THROW("porosity equals zero!! Wrong initial porosity?");

  constitutive_derivatives(
      params, press, J, porosity, params_->init_porosity_, dW_dp, dW_dphi, dW_dJ, dW_dphiref, W);
}

void Mat::StructPoro::constitutive_derivatives(Teuchos::ParameterList& params, double press,
    double J, double porosity, double refporosity, double* dW_dp, double* dW_dphi, double* dW_dJ,
    double* dW_dphiref, double* W)
{
  params_->poro_law_->constitutive_derivatives(
      params, press, J, porosity, refporosity, dW_dp, dW_dphi, dW_dJ, dW_dphiref, W);
}

void Mat::StructPoro::VisNames(std::map<std::string, int>& names)
{
  mat_->VisNames(names);
  std::string porosity = "porosity";
  names[porosity] = 1;  // scalar
}

bool Mat::StructPoro::VisData(
    const std::string& name, std::vector<double>& data, int numgp, int eleID)
{
  if (mat_->VisData(name, data, numgp)) return true;
  if (name == "porosity")
  {
    if ((int)data.size() != 1) FOUR_C_THROW("size mismatch");
    data[0] = PorosityAv();
    return true;
  }
  return false;
}

FOUR_C_NAMESPACE_CLOSE
