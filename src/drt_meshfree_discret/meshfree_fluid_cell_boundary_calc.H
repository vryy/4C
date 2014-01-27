/*----------------------------------------------------------------------*/
/*!
\file meshfree_fluid_cell_boundary_calc.H

\brief evaluation of meshfree fluid terms at integration points

<pre>
Maintainer: Keijo Nissen
            nissen@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15253
</pre>

*---------------------------------------------------------------------------*/

#ifndef MESHFREE_FLUID_CELL_BOUNDARY_CALC_H
#define MESHFREE_FLUID_CELL_BOUNDARY_CALC_H

#include "meshfree_fluid_cell_boundary_interface.H"
#include "drt_meshfree_cell_utils.H"
#include "../linalg/linalg_fixedsizematrix.H"
#include "../linalg/linalg_serialdensematrix.H"
#include "../linalg/linalg_serialdensevector.H"
#include "../drt_fem_general/drt_utils_local_connectivity_matrices.H" // definition of GetGaussPointsAtX()
#include "../drt_lib/drt_element.H"                                   // Element::DiscretizationType

namespace DRT
{
  class Condition;

namespace MESHFREE
{
  class MeshfreeDiscretization;
}

namespace ELEMENTS
{
  class MeshfreeFluidBoundary;
  class FluidEleParameter;
  class FluidEleParameterTimInt;

  /// Internal FluidBoundary element implementation
  /*!
    This internal class keeps all the working arrays needed to
    calculate the FluidBoundary element. Additionally the method Sysmat()
    provides a clean and fast element implementation.

    <h3>Purpose</h3>

    The idea is to separate the element maintenance (class FluidBoundary)
    from the mathematical contents (this class). Of course there are
    different implementations of the FluidBoundary element, this is just one
    such implementation.

    The FluidBoundary element will allocate exactly one object of this class
    for all FluidBoundary elements with the same number of nodes in the mesh.
    This allows us to use exactly matching working arrays (and keep them
    around.)

    The code is meant to be as clean as possible. This is the only way
    to keep it fast. The number of working arrays has to be reduced to
    a minimum so that the element fits into the cache. (There might be
    room for improvements.)

    <h3>History</h3>

    The implementation here is the standard convection-diffusion element
    capable of dealing with systems of transported scalars.

    Right now we do not read any stabilization parameters from the
    input file but have a fixed version.

    \author nis
    \date Nov13
  */
template<DRT::Element::DiscretizationType distype>
class MeshfreeFluidBoundaryCalc: public MeshfreeFluidBoundaryInterface
{
  /*========================================================================*/
  //! @name friends
  /*========================================================================*/

  friend class FluidEleParameter;

  /*========================================================================*/
  //! @name Static member variables
  /*========================================================================*/

  //! number of spatial dimensions of boundary element
  static const int bdrynsd_ = DRT::UTILS::DisTypeToDim<distype>::dim;

  //! number of points of boundary element
  static const int bdrynek_ = DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;

  //! number of spatial dimensions of parent element
  static const int nsd_ = bdrynsd_+1;

  //! number of Gauss points
  static const int ngp_ = DRT::MESHFREE::DisTypeToNumGaussPoints<distype>::ngp;

  //! number of degrees of freedom per node
  static const int numdofpernode_ = nsd_ + 1;

public:
  /*========================================================================*/
  //! @name Constructors and destructors and related methods
  /*========================================================================*/
  MeshfreeFluidBoundaryCalc();

  /*========================================================================*/
  //! @name Evaluate methods
  /*========================================================================*/

  //! Evaluate something at the Neumann boundary
  void EvaluateAction(
    DRT::ELEMENTS::MeshfreeFluidBoundary* ele1,
    Teuchos::ParameterList&               params,
    DRT::Discretization&                  discretization,
    std::vector<int>&                     lm,
    Epetra_SerialDenseMatrix&             elemat1,
    Epetra_SerialDenseMatrix&             elemat2,
    Epetra_SerialDenseVector&             elevec1,
    Epetra_SerialDenseVector&             elevec2,
    Epetra_SerialDenseVector&             elevec3);

  //! Evaluate a Neumann boundary condition
  virtual int EvaluateNeumann(
    DRT::ELEMENTS::MeshfreeFluidBoundary* ele,
    Teuchos::ParameterList&               params,
    DRT::Discretization&                  discretization,
    DRT::Condition&                       condition,
    std::vector<int>&                     lm,
    Epetra_SerialDenseVector&             elevec1_epetra,
    Epetra_SerialDenseMatrix*             elemat1);

  /*========================================================================*/
  //! @name Query methods
  /*========================================================================*/

  //! get density
  void GetDensity(
    Teuchos::RCP<const MAT::Material> material);  ///< reference pointer to material

protected:

  /*========================================================================*/
  //! @name Dofs and nodes
  /*========================================================================*/

  //! number of boundary element nodes
  int bdrynen_;

  //! reference to discretisation
  DRT::MESHFREE::MeshfreeDiscretization* discret_;

  /*========================================================================*/
  //! @name other - te be kept clean!!!
  /*========================================================================*/

  //! pointer to parameter list for time integration
  DRT::ELEMENTS::FluidEleParameterTimInt* fldparatimint_;
  //! pointer to parameter list
  DRT::ELEMENTS::FluidEleParameter* fldpara_;

  //! node coordinates for boundary cell
  LINALG::SerialDenseMatrix nxyz_;
  //! point coordinates for boundary cell
  LINALG::SerialDenseMatrix kxyz_;
  //! Gauss point coordinates for boundary cell
  LINALG::SerialDenseMatrix gxyz_;
  //! Gauss point weights for boundary cell
  LINALG::SerialDenseVector gw_;
  //! array for shape functions for boundary element
  LINALG::SerialDenseVector funct_;
  //! array for shape function derivatives for boundary element
  LINALG::SerialDenseMatrix deriv_;
  //! infinitesimal area element drs
  double drs_;
  //! integration factor
  double fac_;
  //! physical viscosity
  double visc_;
  //! density at t_(n+alpha_F) or t_(n+1)
  double densaf_;
  //! density factor for Neumann boundary conditions (set according to problem type)
  double densfac_;

  }; // class MeshfreeFluidBoundaryCalc

} // namespace ELEMENTS
} // namespace DRT

#endif

