/*----------------------------------------------------------------------*/
/*! \file
 \brief wrapper for structure material of porous media


\level 2
 *-----------------------------------------------------------------------*/

#include "baci_mat_structporo.H"

#include "baci_comm_utils_factory.H"  // for function Factory in Unpack
#include "baci_lib_globalproblem.H"
#include "baci_mat_par_bundle.H"
#include "baci_mat_poro_law.H"
#include "baci_mat_so3_material.H"

#include <vector>

BACI_NAMESPACE_OPEN

MAT::PAR::StructPoro::StructPoro(Teuchos::RCP<MAT::PAR::Material> matdata)
    : Parameter(matdata),
      matid_(matdata->GetInt("MATID")),
      poro_law_ID_(matdata->GetInt("POROLAWID")),
      init_porosity_(matdata->GetDouble("INITPOROSITY"))
{
  // retrieve problem instance to read from
  const int probinst = DRT::Problem::Instance()->Materials()->GetReadFromProblem();

  // for the sake of safety
  if (DRT::Problem::Instance(probinst)->Materials() == Teuchos::null)
    dserror("List of materials cannot be accessed in the global problem instance.");
  // yet another safety check
  if (DRT::Problem::Instance(probinst)->Materials()->Num() == 0)
    dserror("List of materials in the global problem instance is empty.");

  // retrieve validated input line of material ID in question
  Teuchos::RCP<MAT::PAR::Material> curmat =
      DRT::Problem::Instance(probinst)->Materials()->ById(poro_law_ID_);

  switch (curmat->Type())
  {
    case INPAR::MAT::m_poro_law_linear:
    {
      if (curmat->Parameter() == nullptr) curmat->SetParameter(new MAT::PAR::PoroLawLinear(curmat));
      poro_law_ = static_cast<MAT::PAR::PoroLaw*>(curmat->Parameter());
      break;
    }
    case INPAR::MAT::m_poro_law_constant:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::PoroLawConstant(curmat));
      poro_law_ = static_cast<MAT::PAR::PoroLaw*>(curmat->Parameter());
      break;
    }
    case INPAR::MAT::m_poro_law_logNeoHooke_Penalty:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::PoroLawNeoHooke(curmat));
      poro_law_ = static_cast<MAT::PAR::PoroLaw*>(curmat->Parameter());
      break;
    }
    case INPAR::MAT::m_poro_law_incompr_skeleton:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::PoroLawIncompSkeleton(curmat));
      poro_law_ = static_cast<MAT::PAR::PoroLawIncompSkeleton*>(curmat->Parameter());
      break;
    }
    case INPAR::MAT::m_poro_law_linear_biot:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::PoroLawLinBiot(curmat));
      poro_law_ = static_cast<MAT::PAR::PoroLawLinBiot*>(curmat->Parameter());
      break;
    }
    case INPAR::MAT::m_poro_law_density_dependent:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::PoroLawDensityDependent(curmat));
      poro_law_ = static_cast<MAT::PAR::PoroLawDensityDependent*>(curmat->Parameter());
      break;
    }
    default:
      dserror("invalid material for porosity law %d", curmat->Type());
      break;
  }
}

Teuchos::RCP<MAT::Material> MAT::PAR::StructPoro::CreateMaterial()
{
  return Teuchos::rcp(new MAT::StructPoro(this));
}

MAT::StructPoroType MAT::StructPoroType::instance_;

CORE::COMM::ParObject* MAT::StructPoroType::Create(const std::vector<char>& data)
{
  auto* struct_poro = new MAT::StructPoro();
  struct_poro->Unpack(data);
  return struct_poro;
}

MAT::StructPoro::StructPoro()
    : params_(nullptr),
      mat_(Teuchos::null),
      porosity_(Teuchos::null),
      surf_porosity_(Teuchos::null),
      is_initialized_(false)
{
}

MAT::StructPoro::StructPoro(MAT::PAR::StructPoro* params)
    : params_(params),
      porosity_(Teuchos::null),
      surf_porosity_(Teuchos::null),
      is_initialized_(false)
{
  mat_ = Teuchos::rcp_dynamic_cast<MAT::So3Material>(MAT::Material::Factory(params_->matid_));
  if (mat_ == Teuchos::null)
    dserror("MAT::StructPoro: underlying material should be of type MAT::So3Material");
}

void MAT::StructPoro::PoroSetup(int numgp, DRT::INPUT::LineDefinition* linedef)
{
  porosity_ = Teuchos::rcp(new std::vector<double>(numgp, params_->init_porosity_));
  surf_porosity_ = Teuchos::rcp(new std::map<int, std::vector<double>>);

  is_initialized_ = true;
}

inline INPAR::MAT::MaterialType MAT::StructPoro::PoroLawType() const
{
  return params_->poro_law_->Type();
}

double MAT::StructPoro::InvBulkModulus() const { return params_->poro_law_->InvBulkModulus(); }

double MAT::StructPoro::Density() const
{
  if (params_->init_porosity_ == 1.0)
    return mat_->Density();
  else
    return ((1.0 - params_->init_porosity_) * mat_->Density());
}

double MAT::StructPoro::DensitySolidPhase() const { return mat_->Density(); }

void MAT::StructPoro::Pack(CORE::COMM::PackBuffer& data) const
{
  if (not is_initialized_) dserror("poro material not initialized. Not a poro element?");

  CORE::COMM::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data, type);

  // matid
  int matid = -1;
  if (params_ != nullptr) matid = params_->Id();  // in case we are in post-process mode
  AddtoPack(data, matid);

  // porosity_
  int size = static_cast<int>(porosity_->size());
  AddtoPack(data, size);
  for (int i = 0; i < size; ++i)
  {
    AddtoPack(data, (*porosity_)[i]);
  }

  // surf_porosity_
  size = static_cast<int>(surf_porosity_->size());
  AddtoPack(data, size);
  // iterator
  std::map<int, std::vector<double>>::const_iterator iter;
  for (iter = surf_porosity_->begin(); iter != surf_porosity_->end(); ++iter)
  {
    AddtoPack(data, iter->first);
    AddtoPack(data, iter->second);
  }

  // Pack data of underlying material
  if (mat_ != Teuchos::null) mat_->Pack(data);
}

void MAT::StructPoro::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  CORE::COMM::ExtractAndAssertId(position, data, UniqueParObjectId());

  // matid
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
        params_ = static_cast<MAT::PAR::StructPoro*>(mat);
      else
        dserror("Type of parameter material %d does not fit to calling type %d", mat->Type(),
            MaterialType());
    }
  }

  // porosity_
  int size = 0;
  ExtractfromPack(position, data, size);
  porosity_ = Teuchos::rcp(new std::vector<double>);
  double tmp = 0.0;
  for (int i = 0; i < size; ++i)
  {
    ExtractfromPack(position, data, tmp);
    porosity_->push_back(tmp);
  }

  // surface porosity (i think it is not necessary to pack/unpack this...)
  ExtractfromPack(position, data, size);
  surf_porosity_ = Teuchos::rcp(new std::map<int, std::vector<double>>);
  for (int i = 0; i < size; i++)
  {
    int dof;
    std::vector<double> value;
    ExtractfromPack(position, data, dof);
    ExtractfromPack(position, data, value);

    // add to map
    surf_porosity_->insert(std::pair<int, std::vector<double>>(dof, value));
  }

  // Unpack data of sub material (these lines are copied from element.cpp)
  std::vector<char> datamat;
  ExtractfromPack(position, data, datamat);
  if (datamat.size() > 0)
  {
    CORE::COMM::ParObject* o = CORE::COMM::Factory(datamat);  // Unpack is done here
    auto* mat = dynamic_cast<MAT::So3Material*>(o);
    if (mat == nullptr) dserror("failed to unpack elastic material");
    mat_ = Teuchos::rcp(mat);
  }
  else
    mat_ = Teuchos::null;

  is_initialized_ = true;
}

void MAT::StructPoro::ComputePorosity(const double& refporosity, const double& press,
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

  params_->poro_law_->ComputePorosity(refporosity, press, J, gp, porosity, dphi_dp, dphi_dJ,
      dphi_dJdp, dphi_dJJ, dphi_dpp, dphi_dphiref);

  // save porosity
  if (save) porosity_->at(gp) = porosity;
}

void MAT::StructPoro::ComputePorosity(Teuchos::ParameterList& params, double press, double J,
    int gp, double& porosity, double* dphi_dp, double* dphi_dJ, double* dphi_dJdp, double* dphi_dJJ,
    double* dphi_dpp, bool save)
{
  ComputePorosity(
      params_
          ->init_porosity_,  // reference porosity equals initial porosity for non reactive material
      press, J, gp, porosity, dphi_dp, dphi_dJ, dphi_dJdp, dphi_dJJ, dphi_dpp,
      nullptr,  // reference porosity is constant (non reactive) -> derivative not needed
      save);
}

void MAT::StructPoro::ComputePorosity(
    Teuchos::ParameterList& params, double press, double J, int gp, double& porosity, bool save)
{
  ComputePorosity(
      params, press, J, gp, porosity, nullptr, nullptr, nullptr, nullptr, nullptr, save);
}

void MAT::StructPoro::ComputeSurfPorosity(Teuchos::ParameterList& params, double press, double J,
    const int surfnum, int gp, double& porosity, double* dphi_dp, double* dphi_dJ,
    double* dphi_dJdp, double* dphi_dJJ, double* dphi_dpp, bool save)
{
  ComputePorosity(
      params, press, J, gp, porosity, dphi_dp, dphi_dJ, dphi_dJdp, dphi_dJJ, dphi_dpp, save);

  if (save)
  {
    if (gp == 0)  // it's a new iteration, so old values are not needed any more
      ((*surf_porosity_)[surfnum]).clear();

    ((*surf_porosity_)[surfnum]).push_back(porosity);
  }
}

void MAT::StructPoro::ComputeSurfPorosity(Teuchos::ParameterList& params, double press, double J,
    const int surfnum, int gp, double& porosity, bool save)
{
  ComputeSurfPorosity(
      params, press, J, surfnum, gp, porosity, nullptr, nullptr, nullptr, nullptr, nullptr, save);
}


double MAT::StructPoro::PorosityAv() const
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

void MAT::StructPoro::CouplStress(const CORE::LINALG::Matrix<3, 3>& defgrd, const double& press,
    CORE::LINALG::Matrix<6, 1>& couplstress) const
{
  const double J = defgrd.Determinant();

  // Right Cauchy-Green tensor = F^T * F
  CORE::LINALG::Matrix<3, 3> cauchygreen;
  cauchygreen.MultiplyTN(defgrd, defgrd);

  // inverse Right Cauchy-Green tensor
  CORE::LINALG::Matrix<3, 3> C_inv;
  C_inv.Invert(cauchygreen);

  // inverse Right Cauchy-Green tensor as vector
  CORE::LINALG::Matrix<6, 1> C_inv_vec;
  for (int i = 0, k = 0; i < 3; i++)
    for (int j = 0; j < 3 - i; j++, k++) C_inv_vec(k) = C_inv(i + j, j);

  for (int i = 0; i < 6; i++) couplstress(i) = -1.0 * J * press * C_inv_vec(i);
}

void MAT::StructPoro::CouplStress(const CORE::LINALG::Matrix<2, 2>& defgrd, const double& press,
    CORE::LINALG::Matrix<4, 1>& couplstress) const
{
  const double J = defgrd.Determinant();

  // Right Cauchy-Green tensor = F^T * F
  CORE::LINALG::Matrix<2, 2> cauchygreen;
  cauchygreen.MultiplyTN(defgrd, defgrd);

  // inverse Right Cauchy-Green tensor
  CORE::LINALG::Matrix<2, 2> C_inv;
  C_inv.Invert(cauchygreen);

  // inverse Right Cauchy-Green tensor as vector
  CORE::LINALG::Matrix<3, 1> C_inv_vec;
  for (int i = 0, k = 0; i < 2; i++)
    for (int j = 0; j < 2 - i; j++, k++) C_inv_vec(k) = C_inv(i + j, j);

  couplstress(0) = -1.0 * J * press * C_inv_vec(0);
  couplstress(1) = -1.0 * J * press * C_inv_vec(1);
  couplstress(2) =
      0.0;  // this is needed to be compatible with the implementation of the wall element
  couplstress(3) = -1.0 * J * press * C_inv_vec(2);
}

void MAT::StructPoro::ConstitutiveDerivatives(Teuchos::ParameterList& params, double press,
    double J, double porosity, double* dW_dp, double* dW_dphi, double* dW_dJ, double* dW_dphiref,
    double* W)
{
  if (porosity == 0.0) dserror("porosity equals zero!! Wrong initial porosity?");

  ConstitutiveDerivatives(
      params, press, J, porosity, params_->init_porosity_, dW_dp, dW_dphi, dW_dJ, dW_dphiref, W);
}

void MAT::StructPoro::ConstitutiveDerivatives(Teuchos::ParameterList& params, double press,
    double J, double porosity, double refporosity, double* dW_dp, double* dW_dphi, double* dW_dJ,
    double* dW_dphiref, double* W)
{
  params_->poro_law_->ConstitutiveDerivatives(
      params, press, J, porosity, refporosity, dW_dp, dW_dphi, dW_dJ, dW_dphiref, W);
}

void MAT::StructPoro::VisNames(std::map<std::string, int>& names)
{
  mat_->VisNames(names);
  std::string porosity = "porosity";
  names[porosity] = 1;  // scalar
}

bool MAT::StructPoro::VisData(
    const std::string& name, std::vector<double>& data, int numgp, int eleID)
{
  if (mat_->VisData(name, data, numgp)) return true;
  if (name == "porosity")
  {
    if ((int)data.size() != 1) dserror("size mismatch");
    data[0] = PorosityAv();
    return true;
  }
  return false;
}

BACI_NAMESPACE_CLOSE
