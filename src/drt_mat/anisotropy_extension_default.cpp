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

MAT::DefaultAnisotropyExtension::DefaultAnisotropyExtension(const int init_mode, const double gamma,
    const bool adapt_angle,
    const Teuchos::RCP<ELASTIC::StructuralTensorStrategyBase>& stucturalTensorStrategy)
    : FiberAnisotropyExtension(stucturalTensorStrategy),
      init_mode_(init_mode),
      gamma_(gamma),
      adapt_angle_(adapt_angle)
{
  if (init_mode_ == INIT_MODE_NODAL_FIBERS || init_mode_ == INIT_MODE_NODAL_EXTERNAL)
  {
    SetFiberLocation(FiberLocation::GPFibers);
  }
  else
  {
    SetFiberLocation(FiberLocation::ElementFibers);
  }
}

void MAT::DefaultAnisotropyExtension::SetFiberVecs(
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

  a1_0.Multiply(idefgrd, ca1);
  a1_0.Scale(1.0 / a1_0.Norm2());

  a2_0.Multiply(idefgrd, ca2);
  a2_0.Scale(1.0 / a2_0.Norm2());

  std::vector<LINALG::Matrix<3, 1>> fibers(0);
  fibers.emplace_back(a1_0);
  fibers.emplace_back(a2_0);

  SetFibers(BaseAnisotropyExtension::GPDEFAULT, fibers);
}

void MAT::DefaultAnisotropyExtension::SetFiberVecs(const LINALG::Matrix<3, 1>& fibervec)
{
  std::vector<LINALG::Matrix<3, 1>> fibers(0);
  fibers.emplace_back(fibervec);

  SetFibers(BaseAnisotropyExtension::GPDEFAULT, fibers);
}

bool MAT::DefaultAnisotropyExtension::DoElementFiberInitialization()
{
  switch (init_mode_)
  {
    case INIT_MODE_ELEMENT_EXTERNAL:
      DoExternalFiberInitialization();
      return true;
    case INIT_MODE_ELEMENT_FIBERS:

      // check, whether a coordinate system is given
      if (GetAnisotropy()->HasElementCylinderCoordinateSystem())
      {
        // initialize fiber vector with local coordinate system
        LINALG::Matrix<3, 3> locsys(true);
        LINALG::Matrix<3, 3> Id(true);
        MAT::IdentityMatrix(Id);
        GetAnisotropy()->GetElementCylinderCoordinateSystem().EvaluateLocalCoordinateSystem(locsys);

        SetFiberVecs(-1.0, locsys, Id);
      }
      else if (GetAnisotropy()->GetNumberOfElementFibers() > 0)
      {
        // initialize fibers from global given fibers
        SetFibers(GPDEFAULT, GetAnisotropy()->GetElementFibers());
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

bool MAT::DefaultAnisotropyExtension::DoGPFiberInitialization()
{
  switch (init_mode_)
  {
    case INIT_MODE_NODAL_EXTERNAL:
      DoExternalFiberInitialization();
      return true;
    case INIT_MODE_NODAL_FIBERS:

      // check, whether a coordinate system is given
      if (GetAnisotropy()->HasGPCylinderCoordinateSystem())
      {
        dserror(
            "Gauss-point fibers defined via Gauss-point cylinder coordinate systems is not yet "
            "defined");
      }
      else if (GetAnisotropy()->GetNumberOfGPFibers() > 0)
      {
        // initialize fibers from global given fibers
        SetFibers(GetAnisotropy()->GetGPFibers());
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

void MAT::DefaultAnisotropyExtension::DoExternalFiberInitialization()
{
  LINALG::Matrix<3, 3> Id(false);
  MAT::IdentityMatrix(Id);
  SetFiberVecs(-1.0, Id, Id);
}
