/*!-------------------------------------------------------------------------*\
 * \file meshfree_fluid_cell_interface.H
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

#ifndef MESHFREE_FLUID_CELL_INTERFACE_H
#define MESHFREE_FLUID_CELL_INTERFACE_H

#include <vector>
#include <map>
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_SerialDenseMatrix.h"

#include "../drt_lib/drt_singletondestruction.H"
#include "../drt_fem_general/drt_utils_gausspoints.H"

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
  class MeshfreeFluid;
  class MeshfreeFluidBoundary;

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Interface base class for MeshfreeFluidCellCalc
   *
   * This pure virtual class exists to provide a common interface form which
   * all template versions of MeshfreeFluidCell can be derived - so far only
   * the standard one. This is done to stay consistent in structure with the
   * finite element fluid code.
   *
   * \author Keijo Nissen (nissen@lnm.mw.tum.de)
   *
   * \date November, 2013
   */
  /*------------------------------------------------------------------------*/
  class MeshfreeFluidCellInterface : public DRT::SingletonDestruction
  {
  public:

    /*======================================================================*/
    //! @name Constructors and destructors and related methods
    /*======================================================================*/

    /*----------------------------------------------------------------------*/
    /*!
     * \brief Empty constructor
     *///                                                 (public) nis Nov13
    /*----------------------------------------------------------------------*/
    MeshfreeFluidCellInterface() {}

    /*======================================================================*/
    //! @name Query methods
    /*======================================================================*/

    /*----------------------------------------------------------------------*/
    /*!
     * \brief Evaluate the cell - virtual function
     *
     * This class does not provide a definition for this function, it must
     * be defined in MeshfreeFluidCellCalc.
     *///                                                 (public) nis Nov13
    /*----------------------------------------------------------------------*/
    virtual int Evaluate(
      DRT::ELEMENTS::MeshfreeFluid* ele,
      DRT::Discretization &         discretization,
      const std::vector<int>&       lm,
      Teuchos::ParameterList&       params,
      Teuchos::RCP<MAT::Material>&  mat,
      Epetra_SerialDenseMatrix&     elemat1_epetra,
      Epetra_SerialDenseMatrix&     elemat2_epetra,
      Epetra_SerialDenseVector&     elevec1_epetra,
      Epetra_SerialDenseVector&     elevec2_epetra,
      Epetra_SerialDenseVector&     elevec3_epetra,
      bool                          offdiag = false
      ) = 0;

    /// Evaluate supporting methods of the element
    virtual int EvaluateService(
      DRT::ELEMENTS::MeshfreeFluid* ele,
      Teuchos::ParameterList&       params,
      Teuchos::RCP<MAT::Material>&  mat,
      DRT::Discretization&          discretization,
      std::vector<int>&             lm,
      Epetra_SerialDenseMatrix&     elemat1,
      Epetra_SerialDenseMatrix&     elemat2,
      Epetra_SerialDenseVector&     elevec1,
      Epetra_SerialDenseVector&     elevec2,
      Epetra_SerialDenseVector&     elevec3
      ) = 0;

    virtual int ComputeError(
      DRT::ELEMENTS::MeshfreeFluid*       ele,
      Teuchos::ParameterList&             params,
      Teuchos::RCP<MAT::Material>&        mat,
      DRT::Discretization&                discretization,
      std::vector<int>&                   lm,
      Epetra_SerialDenseVector&           elevec1,
      const DRT::UTILS::GaussIntegration& intpoints2
      ) = 0;

  };

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Interface base class for MeshfreeFluidCellCalc
   *
   * This pure virtual class exists to provide a common interface form which
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

