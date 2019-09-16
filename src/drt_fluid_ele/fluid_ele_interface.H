/*----------------------------------------------------------------------*/
/*! \file

\brief Interface class for the evaluation routines of the fluid element

\maintainer Martin Kronbichler

\level 1

*/
/*----------------------------------------------------------------------*/

#ifndef FLUID_ELE_INTERFACE_H
#define FLUID_ELE_INTERFACE_H

#include <vector>
#include <map>
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Epetra_Vector.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_SerialDenseMatrix.h"

#include "../drt_lib/drt_singletondestruction.H"
#include "../drt_fem_general/drt_utils_gausspoints.H"

#include "../drt_cut/cut_utils.H"

namespace MAT
{
  class Material;
}

namespace GEO
{
  namespace CUT
  {
    class BoundaryCell;
    class VolumeCell;
  }  // namespace CUT
}  // namespace GEO

namespace XFEM
{
  class ConditionManager;
  class MeshCouplingFluidFluid;
}  // namespace XFEM

namespace DRT
{
  class Discretization;

  namespace ELEMENTS
  {
    class Fluid;

    /// Interface base class for FluidEleCalc
    /*!
      This class exists to provide a common interface for all template
      versions of FluidEleCalc. The only function this class actually defines
      is Ele, which returns a pointer to the appropriate version of FluidEleCalc.
     */
    class FluidEleInterface : public DRT::SingletonDestruction
    {
     public:
      /// Empty constructor
      FluidEleInterface() {}

      /// Evaluate the element
      /*!
        This class does not provide a definition for this function; it
        must be defined in FluidEleCalc.
       */
      virtual int Evaluate(DRT::ELEMENTS::Fluid* ele, DRT::Discretization& discretization,
          const std::vector<int>& lm, Teuchos::ParameterList& params,
          Teuchos::RCP<MAT::Material>& mat, Epetra_SerialDenseMatrix& elemat1_epetra,
          Epetra_SerialDenseMatrix& elemat2_epetra, Epetra_SerialDenseVector& elevec1_epetra,
          Epetra_SerialDenseVector& elevec2_epetra, Epetra_SerialDenseVector& elevec3_epetra,
          bool offdiag = false) = 0;

      /// evaluate element at specified Gauss points
      virtual int Evaluate(DRT::ELEMENTS::Fluid* ele, DRT::Discretization& discretization,
          const std::vector<int>& lm, Teuchos::ParameterList& params,
          Teuchos::RCP<MAT::Material>& mat, Epetra_SerialDenseMatrix& elemat1_epetra,
          Epetra_SerialDenseMatrix& elemat2_epetra, Epetra_SerialDenseVector& elevec1_epetra,
          Epetra_SerialDenseVector& elevec2_epetra, Epetra_SerialDenseVector& elevec3_epetra,
          const DRT::UTILS::GaussIntegration& intpoints, bool offdiag = false) = 0;

      /// Evaluate the XFEM cut element
      virtual int EvaluateXFEM(DRT::ELEMENTS::Fluid* ele, DRT::Discretization& discretization,
          const std::vector<int>& lm, Teuchos::ParameterList& params,
          Teuchos::RCP<MAT::Material>& mat, Epetra_SerialDenseMatrix& elemat1_epetra,
          Epetra_SerialDenseMatrix& elemat2_epetra, Epetra_SerialDenseVector& elevec1_epetra,
          Epetra_SerialDenseVector& elevec2_epetra, Epetra_SerialDenseVector& elevec3_epetra,
          const std::vector<DRT::UTILS::GaussIntegration>& intpoints,
          const GEO::CUT::plain_volumecell_set& cells, bool offdiag = false) = 0;

      virtual int IntegrateShapeFunction(DRT::ELEMENTS::Fluid* ele,
          DRT::Discretization& discretization, const std::vector<int>& lm,
          Epetra_SerialDenseVector& elevec1, const DRT::UTILS::GaussIntegration& intpoints) = 0;

      virtual int IntegrateShapeFunctionXFEM(DRT::ELEMENTS::Fluid* ele,
          DRT::Discretization& discretization, const std::vector<int>& lm,
          Epetra_SerialDenseVector& elevec1,
          const std::vector<DRT::UTILS::GaussIntegration>& intpoints,
          const GEO::CUT::plain_volumecell_set& cells) = 0;

      /// Evaluate supporting methods of the element
      virtual int EvaluateService(DRT::ELEMENTS::Fluid* ele, Teuchos::ParameterList& params,
          Teuchos::RCP<MAT::Material>& mat, DRT::Discretization& discretization,
          std::vector<int>& lm, Epetra_SerialDenseMatrix& elemat1,
          Epetra_SerialDenseMatrix& elemat2, Epetra_SerialDenseVector& elevec1,
          Epetra_SerialDenseVector& elevec2, Epetra_SerialDenseVector& elevec3) = 0;

      virtual int ComputeError(DRT::ELEMENTS::Fluid* ele, Teuchos::ParameterList& params,
          Teuchos::RCP<MAT::Material>& mat, DRT::Discretization& discretization,
          std::vector<int>& lm, Epetra_SerialDenseVector& elevec1,
          const DRT::UTILS::GaussIntegration& intpoints2) = 0;

      virtual int ComputeErrorInterface(DRT::ELEMENTS::Fluid* ele,   ///< fluid element
          DRT::Discretization& dis,                                  ///< background discretization
          const std::vector<int>& lm,                                ///< element local map
          const Teuchos::RCP<XFEM::ConditionManager>& cond_manager,  ///< XFEM condition manager
          Teuchos::RCP<MAT::Material>& mat,                          ///< material
          Epetra_SerialDenseVector& ele_interf_norms,  /// squared element interface norms
          const std::map<int, std::vector<GEO::CUT::BoundaryCell*>>& bcells,  ///< boundary cells
          const std::map<int, std::vector<DRT::UTILS::GaussIntegration>>&
              bintpoints,                               ///< boundary integration points
          const GEO::CUT::plain_volumecell_set& vcSet,  ///< set of plain volume cells
          Teuchos::ParameterList& params                ///< parameter list
          ) = 0;

      virtual void ElementXfemInterfaceHybridLM(DRT::ELEMENTS::Fluid* ele,  ///< fluid element
          DRT::Discretization& dis,                                  ///< background discretization
          const std::vector<int>& lm,                                ///< element local map
          const Teuchos::RCP<XFEM::ConditionManager>& cond_manager,  ///< XFEM condition manager
          const std::vector<DRT::UTILS::GaussIntegration>& intpoints,  ///< element gauss points
          const std::map<int, std::vector<GEO::CUT::BoundaryCell*>>& bcells,  ///< boundary cells
          const std::map<int, std::vector<DRT::UTILS::GaussIntegration>>&
              bintpoints,  ///< boundary integration points
          const std::map<int, std::vector<int>>&
              patchcouplm,  ///< lm vectors for coupling elements, key= global coupling side-Id
          std::map<int, std::vector<Epetra_SerialDenseMatrix>>&
              side_coupling,                         ///< side coupling matrices
          Teuchos::ParameterList& params,            ///< parameter list
          Teuchos::RCP<MAT::Material>& mat,          ///< material
          Epetra_SerialDenseMatrix& elemat1_epetra,  ///< local system matrix of intersected element
          Epetra_SerialDenseVector&
              elevec1_epetra,               ///< local element vector of intersected element
          Epetra_SerialDenseMatrix& Cuiui,  ///< coupling matrix of a side with itself
          const GEO::CUT::plain_volumecell_set& vcSet  ///< set of plain volume cells
          ) = 0;

      /// add interface condition at cut to element matrix and rhs (two-sided Nitsche coupling)
      virtual void ElementXfemInterfaceNIT(DRT::ELEMENTS::Fluid* ele,  ///< fluid element
          DRT::Discretization& dis,                                  ///< background discretization
          const std::vector<int>& lm,                                ///< element local map
          const Teuchos::RCP<XFEM::ConditionManager>& cond_manager,  ///< XFEM condition manager
          const std::map<int, std::vector<GEO::CUT::BoundaryCell*>>& bcells,  ///< boundary cells
          const std::map<int, std::vector<DRT::UTILS::GaussIntegration>>&
              bintpoints,  ///< boundary integration points
          const std::map<int, std::vector<int>>& patchcouplm,
          Teuchos::ParameterList& params,               ///< parameter list
          Teuchos::RCP<MAT::Material>& mat_master,      ///< material for the coupled side
          Teuchos::RCP<MAT::Material>& mat_slave,       ///< material for the coupled side
          Epetra_SerialDenseMatrix& elemat1_epetra,     ///< element matrix
          Epetra_SerialDenseVector& elevec1_epetra,     ///< element vector
          const GEO::CUT::plain_volumecell_set& vcSet,  ///< volumecell sets in this element
          std::map<int, std::vector<Epetra_SerialDenseMatrix>>&
              side_coupling,               ///< side coupling matrices
          Epetra_SerialDenseMatrix& Cuiui  ///< ui-ui coupling matrix
          ) = 0;

      virtual void CalculateContinuityXFEM(DRT::ELEMENTS::Fluid* ele, DRT::Discretization& dis,
          const std::vector<int>& lm, Epetra_SerialDenseVector& elevec1_epetra,
          const DRT::UTILS::GaussIntegration& intpoints) = 0;

      virtual void CalculateContinuityXFEM(DRT::ELEMENTS::Fluid* ele, DRT::Discretization& dis,
          const std::vector<int>& lm, Epetra_SerialDenseVector& elevec1_epetra) = 0;
    };

  }  // namespace ELEMENTS
}  // namespace DRT

#endif
