/*!-------------------------------------------------------------------------*\
 * \file meshfree_fluid_cell_boundary_interface.H
 *
 * \brief Interface of meshfree fluid cells
 *
 * <pre>
 * Maintainer: Keijo Nissen (nis)
 *             nissen@lnm.mw.tum.de
 *             http://www.lnm.mw.tum.de
 *             089 - 289-15253
 * </pre>
 *
\*--------------------------------------------------------------------------*/

#ifndef MESHFREE_FLUID_CELL_BOUNDARY_INTERFACE_H
#define MESHFREE_FLUID_CELL_BOUNDARY_INTERFACE_H

#include <vector>
#include "Teuchos_ParameterList.hpp"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_SerialDenseMatrix.h"

#include "../drt_lib/drt_singletondestruction.H"

namespace MAT
{
  class Material;
}

namespace DRT
{
  class Discretization;
  class Condition;

namespace ELEMENTS
{
  class MeshfreeFluidBoundary;

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Interface base class for MeshfreeFluidBoundaryCalc
   *
   * This pure virtual class exists to provide a common interface from which
   * all template versions of MeshfreeFluidBoundaryCell can be derived - so
   * far only the standard one. This is done to stay consistent in structure
   * with the finite element fluid code.
   *
   * \author Keijo Nissen (nissen@lnm.mw.tum.de)
   *
   * \date November, 2013
   */
  /*------------------------------------------------------------------------*/
  class MeshfreeFluidBoundaryInterface : DRT::SingletonDestruction
  {

  public:

    /// Empty constructor
    MeshfreeFluidBoundaryInterface() {};

    /// Empty destructor
    virtual ~MeshfreeFluidBoundaryInterface() {};

    virtual void EvaluateAction(
      DRT::ELEMENTS::MeshfreeFluidBoundary* ele1,
      Teuchos::ParameterList&               params,
      DRT::Discretization&                  discretization,
      std::vector<int>&                     lm,
      Epetra_SerialDenseMatrix&             elemat1,
      Epetra_SerialDenseMatrix&             elemat2,
      Epetra_SerialDenseVector&             elevec1,
      Epetra_SerialDenseVector&             elevec2,
      Epetra_SerialDenseVector&             elevec3) = 0;

    virtual int EvaluateNeumann(
      DRT::ELEMENTS::MeshfreeFluidBoundary* ele,
      Teuchos::ParameterList&               params,
      DRT::Discretization&                  discretization,
      DRT::Condition&                       condition,
      std::vector<int>&                     lm,
      Epetra_SerialDenseVector&             elevec1_epetra,
      Epetra_SerialDenseMatrix*             elemat1
    ) = 0;

  };


} //namespace Elements
} //namespace DRT

#endif // MESHFREE_FLUID_CELL_INTERFACE_H

