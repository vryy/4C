/*----------------------------------------------------------------------*/
/*! \file

\brief

\level 3

\maintainer Amadeus Gebauer

*/
/*----------------------------------------------------------------------*/
#include "anisotropy_extension_default.H"
#include "material_service.H"
#include "../headers/standardtypes.h"
#include "anisotropy_extension.H"
#include "../drt_lib/drt_parobject.H"
#include <algorithm>

template <unsigned int numfib>
MAT::DefaultAnisotropyExtension<numfib>::DefaultAnisotropyExtension(const int init_mode,
    const double gamma, const bool adapt_angle,
    const Teuchos::RCP<ELASTIC::StructuralTensorStrategyBase>& stucturalTensorStrategy,
    std::array<int, numfib> fiber_ids)
    : FiberAnisotropyExtension<numfib>(stucturalTensorStrategy),
      init_mode_(init_mode),
      gamma_(gamma),
      adapt_angle_(adapt_angle),
      fiber_ids_(fiber_ids)
{
  if (init_mode_ == INIT_MODE_NODAL_FIBERS || init_mode_ == INIT_MODE_NODAL_EXTERNAL)
  {
    this->SetFiberLocation(FiberLocation::GPFibers);
  }
  else
  {
    this->SetFiberLocation(FiberLocation::ElementFibers);
  }
}

template <unsigned int numfib>
void MAT::DefaultAnisotropyExtension<numfib>::PackAnisotropy(DRT::PackBuffer& data) const
{
  // Call base packing
  MAT::FiberAnisotropyExtension<numfib>::PackAnisotropy(data);

  DRT::ParObject::AddtoPack(data, static_cast<int>(initialized_));
}

template <unsigned int numfib>
void MAT::DefaultAnisotropyExtension<numfib>::UnpackAnisotropy(
    const std::vector<char>& data, std::vector<char>::size_type& position)
{
  // Call base unpacking
  MAT::FiberAnisotropyExtension<numfib>::UnpackAnisotropy(data, position);

  initialized_ = static_cast<bool>(DRT::ParObject::ExtractInt(position, data));
}

template <unsigned int numfib>
void MAT::DefaultAnisotropyExtension<numfib>::SetFiberVecs(
    const double newgamma, const LINALG::Matrix<3, 3>& locsys, const LINALG::Matrix<3, 3>& defgrd)
{
  LINALG::Matrix<3, 1> ca1(true);
  LINALG::Matrix<3, 1> ca2(true);

  // Fiber direction derived from local cosy
  if (init_mode_ == INIT_MODE_ELEMENT_EXTERNAL || init_mode_ == INIT_MODE_ELEMENT_FIBERS)
  {
    // alignment angles gamma_i are read from first entry of then unnecessary vectors a1 and a2
    if ((gamma_ < -90) || (gamma_ > 90)) dserror("Fiber angle not in [-90,90]");
    // convert
    double gamma = (gamma_ * PI) / 180.;

    if (adapt_angle_ && newgamma != -1.0)
    {
      if (gamma * newgamma < 0.0)
      {
        gamma = -1.0 * newgamma;
      }
      else
      {
        gamma = newgamma;
      }
    }

    for (int i = 0; i < 3; ++i)
    {
      // a1 = cos gamma e3 + sin gamma e2
      ca1(i) = std::cos(gamma) * locsys(i, 2) + std::sin(gamma) * locsys(i, 1);
      // a2 = cos gamma e3 - sin gamma e2
      ca2(i) = std::cos(gamma) * locsys(i, 2) - std::sin(gamma) * locsys(i, 1);
    }
  }
  else
  {
    dserror(
        "Setting the fiber vectors is only possible for external element fibers mode or using a "
        "coordinate system.");
  }

  // pull back in reference configuration
  LINALG::Matrix<3, 1> a1_0(true);
  LINALG::Matrix<3, 1> a2_0(true);
  LINALG::Matrix<3, 3> idefgrd(true);
  idefgrd.Invert(defgrd);


  std::array<LINALG::Matrix<3, 1>, numfib> fibers;

  if (numfib >= 1)
  {
    fibers[0].Multiply(idefgrd, ca1);
    fibers[0].Scale(1.0 / fibers[0].Norm2());
  }
  if (numfib >= 2)
  {
    fibers[1].Multiply(idefgrd, ca2);
    fibers[1].Scale(1.0 / fibers[1].Norm2());
  }
  if (numfib >= 3)
  {
    dserror(
        "This kind of initialization method is not implemented for materials that need more than 2 "
        "fibers.");
  }

  this->SetFibers(BaseAnisotropyExtension::GPDEFAULT, fibers);
}

template <unsigned int numfib>
void MAT::DefaultAnisotropyExtension<numfib>::SetFiberVecs(const LINALG::Matrix<3, 1>& fibervec)
{
  std::array<LINALG::Matrix<3, 1>, numfib> fibers;
  fibers[0].Update(fibervec);

  if (numfib >= 2)
  {
    dserror("This method can only be called for materials with one fiber!");
  }

  this->SetFibers(BaseAnisotropyExtension::GPDEFAULT, fibers);
}

template <unsigned int numfib>
bool MAT::DefaultAnisotropyExtension<numfib>::DoElementFiberInitialization()
{
  switch (init_mode_)
  {
    case INIT_MODE_ELEMENT_EXTERNAL:
      DoExternalFiberInitialization();
      return true;
    case INIT_MODE_ELEMENT_FIBERS:

      // check, whether a coordinate system is given
      if (this->GetAnisotropy()->HasElementCylinderCoordinateSystem())
      {
        // initialize fiber vector with local coordinate system
        LINALG::Matrix<3, 3> locsys(true);
        LINALG::Matrix<3, 3> Id(true);
        MAT::IdentityMatrix(Id);
        this->GetAnisotropy()->GetElementCylinderCoordinateSystem().EvaluateLocalCoordinateSystem(
            locsys);

        this->SetFiberVecs(-1.0, locsys, Id);
      }
      else if (this->GetAnisotropy()->GetNumberOfElementFibers() > 0)
      {
        // initialize fibers from global given fibers
        std::array<LINALG::Matrix<3, 1>, numfib> fibers;
        for (unsigned int i = 0; i < numfib; ++i)
        {
          fibers.at(i) = this->GetAnisotropy()->GetElementFibers()[fiber_ids_.at(i)];
        }
        this->SetFibers(BaseAnisotropyExtension::GPDEFAULT, fibers);
      }
      else
      {
        dserror("Could not find element coordinate system or element fibers!");
      }

      return true;
    default:
      return false;
  }
}

template <unsigned int numfib>
bool MAT::DefaultAnisotropyExtension<numfib>::DoGPFiberInitialization()
{
  switch (init_mode_)
  {
    case INIT_MODE_NODAL_EXTERNAL:
      DoExternalFiberInitialization();
      return true;
    case INIT_MODE_NODAL_FIBERS:

      // check, whether a coordinate system is given
      if (this->GetAnisotropy()->HasGPCylinderCoordinateSystem())
      {
        dserror(
            "Gauss-point fibers defined via Gauss-point cylinder coordinate systems is not yet "
            "defined");
      }
      else if (this->GetAnisotropy()->GetNumberOfGPFibers() > 0)
      {
        // initialize fibers from global given fibers
        int gp = 0;
        for (const auto& fiberList : this->GetAnisotropy()->GetGPFibers())
        {
          std::array<LINALG::Matrix<3, 1>, numfib> fibers;

          int i = 0;
          for (int id : fiber_ids_)
          {
            fibers.at(i) = fiberList[id];
            ++i;
          }
          this->SetFibers(gp, fibers);
          ++gp;
        }
      }
      else
      {
        dserror("Could not find Gauss-point coordinate systems or Gauss-point fibers!");
      }

      return true;
    default:
      return false;
  }
}

template <unsigned int numfib>
void MAT::DefaultAnisotropyExtension<numfib>::DoExternalFiberInitialization()
{
  LINALG::Matrix<3, 3> Id(false);
  MAT::IdentityMatrix(Id);
  SetFiberVecs(-1.0, Id, Id);
}


// explicit instatiations of template classes
template class MAT::DefaultAnisotropyExtension<1u>;
template class MAT::DefaultAnisotropyExtension<2u>;