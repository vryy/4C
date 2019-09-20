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

MAT::Anisotropy::Anisotropy(int number_fibers) : number_fibers_(number_fibers)
{
  // empty
}

MAT::Anisotropy::Anisotropy(int number_fibers, int init,
    const Teuchos::RCP<ELASTIC::StructuralTensorStrategyBase>& structuralTensorStrategy)
    : number_fibers_(number_fibers),
      init_mode_(init),
      structuralTensorStrategy_(structuralTensorStrategy)
{
  // empty
}

void MAT::Anisotropy::PackAnisotropy(DRT::PackBuffer& data) const
{
  DRT::ParObject::AddtoPack(data, numgp_);
  DRT::ParObject::AddtoPack(data, fibers_);
  DRT::ParObject::AddtoPack(data, structuralTensors_stress_);
  DRT::ParObject::AddtoPack(data, structuralTensors_);
  DRT::ParObject::AddtoPack(data, static_cast<const int>(fibers_initialized_));
  DRT::ParObject::AddtoPack(data, definitionMode_);
}

void MAT::Anisotropy::UnpackAnisotropy(
    const std::vector<char>& data, std::vector<char>::size_type& position)
{
  DRT::ParObject::ExtractfromPack(position, data, numgp_);
  DRT::ParObject::ExtractfromPack(position, data, fibers_);
  DRT::ParObject::ExtractfromPack(position, data, structuralTensors_stress_);
  DRT::ParObject::ExtractfromPack(position, data, structuralTensors_);
  fibers_initialized_ = (bool)DRT::ParObject::ExtractInt(position, data);
  DRT::ParObject::ExtractfromPack(position, data, definitionMode_);
}

void MAT::Anisotropy::Initialize(
    int init, const Teuchos::RCP<ELASTIC::StructuralTensorStrategyBase>& structuralTensorStrategy)
{
  init_mode_ = init;
  structuralTensorStrategy_ = structuralTensorStrategy;
}

void MAT::Anisotropy::ReadAnisotropyFromElement(
    int const numgp, DRT::INPUT::LineDefinition* linedef)
{
  numgp_ = numgp;
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
      for (int i = 0; i < 3; ++i)
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
        ReadAnisotropyFiber(linedef, "FIBER" + std::to_string(i), fibers_[GPDEFAULT][i]);

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
  for (int i = 0; i < 3; ++i)
  {
    f1norm += fiber[i] * fiber[i];
  }
  f1norm = sqrt(f1norm);

  // fill final fiber vector
  for (int i = 0; i < 3; ++i)
  {
    fiber_vector(i) = fiber[i] / f1norm;
  }
}

void MAT::Anisotropy::SetFibers(int gp, const std::vector<LINALG::Matrix<3, 1>>& fibers)
{
  // only save fibers if initialization mode is for nodes
  if (init_mode_ != INIT_MODE_NODAL_FIBERS and init_mode_ != INIT_MODE_EXTERNAL)
  {
    // Do nothing in this case
    return;
  }
  definitionMode_ = DefinitionMode::external;

  // first check, whether the size is correct (only in DEBUG mode)
#ifdef DEBUG
  if (fibers.size() != number_fibers_)
  {
    dserror("The number of fibers does not match. Given %i fibers but expected %i.", fibers.size(),
        number_fibers_);
  }

  if (gp >= numgp_)
  {
    dserror("The current Gauß point %i is out of range of the expected number of Gauß points %i.",
        gp, numgp_);
  }
#endif

  // Store fibers
  for (unsigned i = 0; i < number_fibers_; ++i)
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
  for (unsigned long gp = 0; gp < fibers_.size(); ++gp)
  {
    structuralTensors_stress_[0].resize(fibers_[gp].size());
    for (unsigned long i = 0; i < fibers_[gp].size(); ++i)
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
  for (unsigned long gp = 0; gp < fibers_.size(); ++gp)
  {
    structuralTensors_[gp].resize(fibers_[gp].size());
    for (unsigned long i = 0; i < fibers_[gp].size(); ++i)
    {
      LINALG::Matrix<6, 1> A_stress(false);
      LINALG::Matrix<3, 3> A(false);
      structuralTensorStrategy_->SetupStructuralTensor(fibers_[gp][i], A_stress);

      // Convert from stress like Voigt notation to matrix notation
      VStressUtils::ToTensor(A_stress, A);

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
    gp = params.get<int>("gp", -1);

    if (gp < 0)
    {
      dserror("Somehow the material does not write the Gauß point into the parameter list");
    }
  }

  return gp;
}
