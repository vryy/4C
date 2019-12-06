/*----------------------------------------------------------------------*/
/*! \file

\brief Implementation of common functionality for anisotropic materials

\level 3

\maintainer Amadeus Gebauer

*/
/*----------------------------------------------------------------------*/
#include "anisotropy.H"
#include "material_service.H"
#include "../drt_lib/drt_parobject.H"
#include "../drt_matelast/elast_aniso_structuraltensor_strategy.H"
#include "../drt_fiber/nodal_fiber_holder.H"
#include "../drt_lib/voigt_notation.H"

MAT::Anisotropy::Anisotropy(int number_fibers)
    : number_fibers_(number_fibers),
      fibers_initialized_(false),
      definitionMode_(DefinitionMode::undefined)
{
  // empty
}

MAT::Anisotropy::Anisotropy(int number_fibers, int init,
    const Teuchos::RCP<ELASTIC::StructuralTensorStrategyBase>& structuralTensorStrategy)
    : number_fibers_(number_fibers),
      init_mode_(init),
      fibers_initialized_(false),
      definitionMode_(DefinitionMode::undefined),
      structuralTensorStrategy_(structuralTensorStrategy)
{
  // empty
}

void MAT::Anisotropy::PackAnisotropy(DRT::PackBuffer& data) const
{
  DRT::ParObject::AddtoPack(data, numgp_);
  PackFiberVector<LINALG::Matrix<3, 1>>(data, fibers_);
  PackFiberVector<LINALG::Matrix<6, 1>>(data, structuralTensors_stress_);
  PackFiberVector<LINALG::Matrix<3, 3>>(data, structuralTensors_);
  DRT::ParObject::AddtoPack(data, static_cast<const int>(fibers_initialized_));
  DRT::ParObject::AddtoPack(data, definitionMode_);
}

void MAT::Anisotropy::UnpackAnisotropy(
    const std::vector<char>& data, std::vector<char>::size_type& position)
{
  DRT::ParObject::ExtractfromPack(position, data, numgp_);
  UnpackFiberVector<LINALG::Matrix<3, 1>>(position, data, fibers_);
  UnpackFiberVector<LINALG::Matrix<6, 1>>(position, data, structuralTensors_stress_);
  UnpackFiberVector<LINALG::Matrix<3, 3>>(position, data, structuralTensors_);
  fibers_initialized_ = (bool)DRT::ParObject::ExtractInt(position, data);
  DRT::ParObject::ExtractfromPack(position, data, definitionMode_);
}

void MAT::Anisotropy::Initialize(
    int init, const Teuchos::RCP<ELASTIC::StructuralTensorStrategyBase>& structuralTensorStrategy)
{
  init_mode_ = init;
  structuralTensorStrategy_ = structuralTensorStrategy;
}

void MAT::Anisotropy::SetNumberOfGaussPoints(int numgp)
{
  numgp_ = numgp;

  // As we now know the number of Gauss points, we can resize the fiber_ vector
  fibers_.resize(numgp);
  for (int gp = 0; gp < numgp; ++gp)
  {
    fibers_[gp].resize(number_fibers_);
  }
}

void MAT::Anisotropy::ReadAnisotropyFromElement(DRT::INPUT::LineDefinition* linedef)
{
  if (init_mode_ == INIT_MODE_EXTERNAL)
  {
    // Call external fiber setup
    DoFiberInitialization();
  }
  else if (init_mode_ == INIT_MODE_ELEMENT_FIBERS)
  {
    // CIR-AXI-RAD nomenclature
    if (linedef->HaveNamed("RAD") and linedef->HaveNamed("AXI") and linedef->HaveNamed("CIR"))
    {
      // Read in of data
      // read local (cylindrical) cosy-directions at current element
      LINALG::Matrix<3, 3> locsys;

      // Single coordinate system unite vectors
      LINALG::Matrix<3, 1> dir_rad;  // radial direction
      LINALG::Matrix<3, 1> dir_axi;  // axial direction
      LINALG::Matrix<3, 1> dir_cir;  // circular direction

      // Read fibers
      ReadAnisotropyFiber(linedef, "RAD", dir_rad);
      ReadAnisotropyFiber(linedef, "AXI", dir_axi);
      ReadAnisotropyFiber(linedef, "CIR", dir_cir);

      // Build local coordinate system
      for (unsigned i = 0; i < 3; ++i)
      {
        locsys(i, 0) = dir_rad(i);
        locsys(i, 1) = dir_axi(i);
        locsys(i, 2) = dir_cir(i);
      }

      // Call to setup fibers by coordinate system
      SetupFiberByCosy(locsys);

      definitionMode_ = DefinitionMode::radaxicir;
    }

    // FIBERi nomenclature
    for (unsigned i = 0; i < number_fibers_; ++i)
    {
      if (linedef->HaveNamed("FIBER" + std::to_string(i + 1)))
      {
        // Read in of fiber data and setting fiber data
        ReadAnisotropyFiber(linedef, "FIBER" + std::to_string(i + 1), fibers_[GPDEFAULT][i]);

        definitionMode_ = DefinitionMode::fiberi;
      }
    }

    // check if all read fibers are of unit length
    for (unsigned i = 0; i < number_fibers_; ++i)
    {
      if (std::abs(fibers_[GPDEFAULT][i].Norm2() - 1.0) > TOLERANCE)
      {
        dserror("Fiber %i is not an unit vector, which is fatal.", i);
      }
    }
    fibers_initialized_ = true;

    // Call the notifier method
    OnFibersInitialized();
  }
  // fibers defined on nodes
  else if (init_mode_ == INIT_MODE_NODAL_FIBERS)
  {
    // nothing to do here. gp fibers are passed directly from the element via
    // MAT::Anisotropy::SetFiber(int, const std::vector<LINALG::Matrix<3, 1>>&)
  }
  else
  {
    dserror("Initialization mode for fibers (INIT=%i) is not implemented!", init_mode_);
  }
}

void MAT::Anisotropy::ReadAnisotropyFiber(
    DRT::INPUT::LineDefinition* linedef, std::string specifier, LINALG::Matrix<3, 1>& fiber_vector)
{
  std::vector<double> fiber;
  linedef->ExtractDoubleVector(std::move(specifier), fiber);
  double f1norm = 0.;
  // normalization
  for (std::vector<double>::size_type i = 0; i < 3; ++i)
  {
    f1norm += fiber[i] * fiber[i];
  }
  f1norm = sqrt(f1norm);

  // fill final fiber vector
  for (std::vector<double>::size_type i = 0; i < 3; ++i)
  {
    fiber_vector(i) = fiber[i] / f1norm;
  }
}

void MAT::Anisotropy::ReadAnisotropyFromParameterList(Teuchos::ParameterList& params)
{
  if (init_mode_ == INIT_MODE_NODAL_FIBERS)
  {
    if (params.isParameter("fiberholder"))
    {
      const auto& fiberHolder = params.get<DRT::FIBER::NodalFiberHolder>("fiberholder");

      std::vector<LINALG::Matrix<3, 1>> gpfiber1 =
          fiberHolder.GetFiber(DRT::FIBER::FiberType::Fiber1);
      std::vector<LINALG::Matrix<3, 1>> gpfiber2 =
          fiberHolder.GetFiber(DRT::FIBER::FiberType::Fiber2);

      std::vector<LINALG::Matrix<3, 1>> fibers(0);
      for (int gp = 0; gp < numgp_; ++gp)
      {
        fibers.emplace_back(gpfiber1[gp]);
        fibers.emplace_back(gpfiber2[gp]);
        SetFibers(gp, fibers);
        fibers.clear();
      }
    }
    else
    {
      dserror(
          "The fibers are not set in the ParameterList. Those fibers have to be computed "
          "in the MaterialPostSetup() routine of the Element.");
    }
  }
}

void MAT::Anisotropy::SetFibers(int gp, const std::vector<LINALG::Matrix<3, 1>>& fibers)
{
  // This method should only be called if fiber initialization mode is external or with nodal fibers
  if (gp != GPDEFAULT and init_mode_ != INIT_MODE_NODAL_FIBERS and init_mode_ != INIT_MODE_EXTERNAL)
  {
    dserror(
        "The initialization mode (INIT=%d) does not allow to set nodal fibers. Nodal fibers "
        "are only allowed with nodal fibers and external fiber initialization",
        init_mode_);
  }
  definitionMode_ = DefinitionMode::external;

  // first check, whether the size is correct
  if (fibers.size() < number_fibers_)
  {
    dserror(
        "The number of given fibers (%i) is less than the needed fibers for this material (%i).",
        fibers.size(), number_fibers_);
  }

  if (gp >= numgp_)
  {
    dserror("The current Gauss point %i is out of range of the expected number of Gauss points %i.",
        gp, numgp_);
  }

  // Store fibers
  for (std::vector<LINALG::Matrix<3, 1>>::size_type i = 0; i < number_fibers_; ++i)
  {
    fibers_[gp][i].Update(fibers[i]);
  }
  fibers_initialized_ = true;

  // Call the notifier method
  OnFibersInitialized();
}

void MAT::Anisotropy::ComputeStructuralTensors_stress()
{
  // Need to compute the stuctural tensors
  if (Teuchos::is_null(structuralTensorStrategy_))
  {
    dserror("Structural tensor strategy is not initialized. Do it in Initialize(...).");
  }

  structuralTensors_stress_.resize(fibers_.size());
  for (std::vector<std::vector<LINALG::Matrix<3, 1>>>::size_type gp = 0; gp < fibers_.size(); ++gp)
  {
    structuralTensors_stress_[gp].resize(fibers_[gp].size());
    for (std::vector<LINALG::Matrix<3, 1>>::size_type i = 0; i < fibers_[gp].size(); ++i)
    {
      LINALG::Matrix<6, 1> A_stress(false);
      structuralTensorStrategy_->SetupStructuralTensor(fibers_[gp][i], A_stress);

      structuralTensors_stress_[gp][i].Update(A_stress);
    }
  }
}

void MAT::Anisotropy::ComputeStructuralTensors()
{
  // Need to compute the stuctural tensors
  if (Teuchos::is_null(structuralTensorStrategy_))
  {
    dserror("Structural tensor strategy is not initialized. Do it in Initialize(...).");
  }

  structuralTensors_.resize(fibers_.size());
  for (std::vector<std::vector<LINALG::Matrix<3, 1>>>::size_type gp = 0; gp < fibers_.size(); ++gp)
  {
    structuralTensors_[gp].resize(fibers_[gp].size());
    for (std::vector<LINALG::Matrix<3, 1>>::size_type i = 0; i < fibers_[gp].size(); ++i)
    {
      LINALG::Matrix<6, 1> A_stress(false);
      LINALG::Matrix<3, 3> A(false);
      structuralTensorStrategy_->SetupStructuralTensor(fibers_[gp][i], A_stress);

      // Convert from stress like Voigt notation to matrix notation
      UTILS::VOIGT::Stresses::VectorToMatrix(A_stress, A);

      structuralTensors_[gp][i].Update(A);
    }
  }
}

LINALG::Matrix<6, 1>& MAT::Anisotropy::GetStructuralTensor_stress(int gp, int i)
{
  if (!fibers_initialized_)
  {
    dserror("The fibers are not yet initialized");
  }
  if (structuralTensors_stress_.empty())
  {
    ComputeStructuralTensors_stress();
  }

  return structuralTensors_stress_[gp][i];
}

LINALG::Matrix<3, 3>& MAT::Anisotropy::GetStructuralTensor(int gp, int i)
{
  if (!fibers_initialized_)
  {
    dserror("The fibers are not yet initialized");
  }
  if (structuralTensors_.empty())
  {
    ComputeStructuralTensors();
  }

  return structuralTensors_[gp][i];
}

LINALG::Matrix<3, 1>& MAT::Anisotropy::GetFiber(int gp, int i)
{
  if (!fibers_initialized_)
  {
    dserror("The fibers are not yet initialized");
  }
  return fibers_[gp][i];
}

int MAT::Anisotropy::GetFibersPerElement()
{
  if (init_mode_ == INIT_MODE_NODAL_FIBERS)
  {
    // Nodal fibers
    // Number of fibers is number of Integration points
    return numgp_;
  }

  // Element fibers
  return 1;
}

int MAT::Anisotropy::GetGPId(Teuchos::ParameterList& params)
{
  int gp = GPDEFAULT;
  if (init_mode_ == INIT_MODE_NODAL_FIBERS)
  {
    // get Gauss point from Parameter list
    gp = params.get<int>("gp", -1);

    if (gp < 0)
    {
      dserror("Somehow the material does not write the Gauss point into the parameter list.");
    }
  }

  return gp;
}

template <typename T>
void MAT::Anisotropy::PackFiberVector(
    DRT::PackBuffer& buffer, const std::vector<std::vector<T>>& vct) const
{
  DRT::ParObject::AddtoPack(buffer, static_cast<int>(vct.size()));

  for (const auto& list : vct)
  {
    DRT::ParObject::AddtoPack(buffer, list);
  }
}

template <typename T>
void MAT::Anisotropy::UnpackFiberVector(std::vector<char>::size_type& position,
    const std::vector<char>& data, std::vector<std::vector<T>>& vct) const
{
  vct.clear();
  int numfibs = DRT::ParObject::ExtractInt(position, data);
  for (int i = 0; i < numfibs; ++i)
  {
    std::vector<T> mat(0);
    DRT::ParObject::ExtractfromPack(position, data, mat);
    vct.emplace_back(mat);
  }
}

/*----------------------------------------------------------------------------*/
// explicit instantiation of template functions
template void MAT::Anisotropy::PackFiberVector(
    DRT::PackBuffer& buffer, const std::vector<std::vector<LINALG::Matrix<3, 3>>>& vct) const;
template void MAT::Anisotropy::UnpackFiberVector(std::vector<char>::size_type& position,
    const std::vector<char>& data, std::vector<std::vector<LINALG::Matrix<3, 3>>>& vct) const;
template void MAT::Anisotropy::PackFiberVector(
    DRT::PackBuffer& buffer, const std::vector<std::vector<LINALG::Matrix<3, 1>>>& vct) const;
template void MAT::Anisotropy::UnpackFiberVector(std::vector<char>::size_type& position,
    const std::vector<char>& data, std::vector<std::vector<LINALG::Matrix<3, 1>>>& vct) const;
template void MAT::Anisotropy::PackFiberVector(
    DRT::PackBuffer& buffer, const std::vector<std::vector<LINALG::Matrix<6, 1>>>& vct) const;
template void MAT::Anisotropy::UnpackFiberVector(std::vector<char>::size_type& position,
    const std::vector<char>& data, std::vector<std::vector<LINALG::Matrix<6, 1>>>& vct) const;