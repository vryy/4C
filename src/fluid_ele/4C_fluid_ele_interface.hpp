/*----------------------------------------------------------------------*/
/*! \file

\brief Interface class for the evaluation routines of the fluid element


\level 1

*/
/*----------------------------------------------------------------------*/

#ifndef FOUR_C_FLUID_ELE_INTERFACE_HPP
#define FOUR_C_FLUID_ELE_INTERFACE_HPP

#include "4C_config.hpp"

#include "4C_cut_utils.hpp"
#include "4C_discretization_fem_general_utils_gausspoints.hpp"
#include "4C_linalg_serialdensematrix.hpp"
#include "4C_linalg_serialdensevector.hpp"

#include <Epetra_Vector.h>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>

#include <map>
#include <vector>

FOUR_C_NAMESPACE_OPEN

namespace Mat
{
  class Material;
}

namespace Core::Geo
{
  namespace Cut
  {
    class BoundaryCell;
    class VolumeCell;
  }  // namespace Cut
}  // namespace Core::Geo

namespace XFEM
{
  class ConditionManager;
  class MeshCouplingFluidFluid;
}  // namespace XFEM

namespace Discret
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
    class FluidEleInterface
    {
     public:
      /**
       * Virtual destructor.
       */
      virtual ~FluidEleInterface() = default;

      /// Empty constructor
      FluidEleInterface() = default;

      /// Evaluate the element
      /*!
        This class does not provide a definition for this function; it
        must be defined in FluidEleCalc.
       */
      virtual int Evaluate(Discret::ELEMENTS::Fluid* ele, Discret::Discretization& discretization,
          const std::vector<int>& lm, Teuchos::ParameterList& params,
          Teuchos::RCP<Core::Mat::Material>& mat, Core::LinAlg::SerialDenseMatrix& elemat1_epetra,
          Core::LinAlg::SerialDenseMatrix& elemat2_epetra,
          Core::LinAlg::SerialDenseVector& elevec1_epetra,
          Core::LinAlg::SerialDenseVector& elevec2_epetra,
          Core::LinAlg::SerialDenseVector& elevec3_epetra, bool offdiag = false) = 0;

      /// evaluate element at specified Gauss points
      virtual int Evaluate(Discret::ELEMENTS::Fluid* ele, Discret::Discretization& discretization,
          const std::vector<int>& lm, Teuchos::ParameterList& params,
          Teuchos::RCP<Core::Mat::Material>& mat, Core::LinAlg::SerialDenseMatrix& elemat1_epetra,
          Core::LinAlg::SerialDenseMatrix& elemat2_epetra,
          Core::LinAlg::SerialDenseVector& elevec1_epetra,
          Core::LinAlg::SerialDenseVector& elevec2_epetra,
          Core::LinAlg::SerialDenseVector& elevec3_epetra,
          const Core::FE::GaussIntegration& intpoints, bool offdiag = false) = 0;

      /// Evaluate the XFEM cut element
      virtual int EvaluateXFEM(Discret::ELEMENTS::Fluid* ele,
          Discret::Discretization& discretization, const std::vector<int>& lm,
          Teuchos::ParameterList& params, Teuchos::RCP<Core::Mat::Material>& mat,
          Core::LinAlg::SerialDenseMatrix& elemat1_epetra,
          Core::LinAlg::SerialDenseMatrix& elemat2_epetra,
          Core::LinAlg::SerialDenseVector& elevec1_epetra,
          Core::LinAlg::SerialDenseVector& elevec2_epetra,
          Core::LinAlg::SerialDenseVector& elevec3_epetra,
          const std::vector<Core::FE::GaussIntegration>& intpoints,
          const Core::Geo::Cut::plain_volumecell_set& cells, bool offdiag = false) = 0;

      virtual int integrate_shape_function(Discret::ELEMENTS::Fluid* ele,
          Discret::Discretization& discretization, const std::vector<int>& lm,
          Core::LinAlg::SerialDenseVector& elevec1,
          const Core::FE::GaussIntegration& intpoints) = 0;

      virtual int integrate_shape_function_xfem(Discret::ELEMENTS::Fluid* ele,
          Discret::Discretization& discretization, const std::vector<int>& lm,
          Core::LinAlg::SerialDenseVector& elevec1,
          const std::vector<Core::FE::GaussIntegration>& intpoints,
          const Core::Geo::Cut::plain_volumecell_set& cells) = 0;

      /// Evaluate supporting methods of the element
      virtual int EvaluateService(Discret::ELEMENTS::Fluid* ele, Teuchos::ParameterList& params,
          Teuchos::RCP<Core::Mat::Material>& mat, Discret::Discretization& discretization,
          std::vector<int>& lm, Core::LinAlg::SerialDenseMatrix& elemat1,
          Core::LinAlg::SerialDenseMatrix& elemat2, Core::LinAlg::SerialDenseVector& elevec1,
          Core::LinAlg::SerialDenseVector& elevec2, Core::LinAlg::SerialDenseVector& elevec3) = 0;

      virtual int compute_error(Discret::ELEMENTS::Fluid* ele, Teuchos::ParameterList& params,
          Teuchos::RCP<Core::Mat::Material>& mat, Discret::Discretization& discretization,
          std::vector<int>& lm, Core::LinAlg::SerialDenseVector& elevec1,
          const Core::FE::GaussIntegration& intpoints2) = 0;

      virtual int compute_error_interface(Discret::ELEMENTS::Fluid* ele,  ///< fluid element
          Discret::Discretization& dis,                              ///< background discretization
          const std::vector<int>& lm,                                ///< element local map
          const Teuchos::RCP<XFEM::ConditionManager>& cond_manager,  ///< XFEM condition manager
          Teuchos::RCP<Core::Mat::Material>& mat,                    ///< material
          Core::LinAlg::SerialDenseVector& ele_interf_norms,  /// squared element interface norms
          const std::map<int, std::vector<Core::Geo::Cut::BoundaryCell*>>&
              bcells,  ///< boundary cells
          const std::map<int, std::vector<Core::FE::GaussIntegration>>&
              bintpoints,                                     ///< boundary integration points
          const Core::Geo::Cut::plain_volumecell_set& vcSet,  ///< set of plain volume cells
          Teuchos::ParameterList& params                      ///< parameter list
          ) = 0;

      virtual void element_xfem_interface_hybrid_lm(
          Discret::ELEMENTS::Fluid* ele,                             ///< fluid element
          Discret::Discretization& dis,                              ///< background discretization
          const std::vector<int>& lm,                                ///< element local map
          const Teuchos::RCP<XFEM::ConditionManager>& cond_manager,  ///< XFEM condition manager
          const std::vector<Core::FE::GaussIntegration>& intpoints,  ///< element gauss points
          const std::map<int, std::vector<Core::Geo::Cut::BoundaryCell*>>&
              bcells,  ///< boundary cells
          const std::map<int, std::vector<Core::FE::GaussIntegration>>&
              bintpoints,  ///< boundary integration points
          const std::map<int, std::vector<int>>&
              patchcouplm,  ///< lm vectors for coupling elements, key= global coupling side-Id
          std::map<int, std::vector<Core::LinAlg::SerialDenseMatrix>>&
              side_coupling,                       ///< side coupling matrices
          Teuchos::ParameterList& params,          ///< parameter list
          Teuchos::RCP<Core::Mat::Material>& mat,  ///< material
          Core::LinAlg::SerialDenseMatrix&
              elemat1_epetra,  ///< local system matrix of intersected element
          Core::LinAlg::SerialDenseVector&
              elevec1_epetra,                      ///< local element vector of intersected element
          Core::LinAlg::SerialDenseMatrix& Cuiui,  ///< coupling matrix of a side with itself
          const Core::Geo::Cut::plain_volumecell_set& vcSet  ///< set of plain volume cells
          ) = 0;

      /// add interface condition at cut to element matrix and rhs (two-sided Nitsche coupling)
      virtual void element_xfem_interface_nit(Discret::ELEMENTS::Fluid* ele,  ///< fluid element
          Discret::Discretization& dis,                              ///< background discretization
          const std::vector<int>& lm,                                ///< element local map
          const Teuchos::RCP<XFEM::ConditionManager>& cond_manager,  ///< XFEM condition manager
          const std::map<int, std::vector<Core::Geo::Cut::BoundaryCell*>>&
              bcells,  ///< boundary cells
          const std::map<int, std::vector<Core::FE::GaussIntegration>>&
              bintpoints,  ///< boundary integration points
          const std::map<int, std::vector<int>>& patchcouplm,
          Teuchos::ParameterList& params,                     ///< parameter list
          Teuchos::RCP<Core::Mat::Material>& mat_master,      ///< material for the coupled side
          Teuchos::RCP<Core::Mat::Material>& mat_slave,       ///< material for the coupled side
          Core::LinAlg::SerialDenseMatrix& elemat1_epetra,    ///< element matrix
          Core::LinAlg::SerialDenseVector& elevec1_epetra,    ///< element vector
          const Core::Geo::Cut::plain_volumecell_set& vcSet,  ///< volumecell sets in this element
          std::map<int, std::vector<Core::LinAlg::SerialDenseMatrix>>&
              side_coupling,                       ///< side coupling matrices
          Core::LinAlg::SerialDenseMatrix& Cuiui,  ///< ui-ui coupling matrix
          bool evaluated_cut  ///< the CUT was updated before this evaluation is called
          ) = 0;

      virtual void calculate_continuity_xfem(Discret::ELEMENTS::Fluid* ele,
          Discret::Discretization& dis, const std::vector<int>& lm,
          Core::LinAlg::SerialDenseVector& elevec1_epetra,
          const Core::FE::GaussIntegration& intpoints) = 0;

      virtual void calculate_continuity_xfem(Discret::ELEMENTS::Fluid* ele,
          Discret::Discretization& dis, const std::vector<int>& lm,
          Core::LinAlg::SerialDenseVector& elevec1_epetra) = 0;
    };

  }  // namespace ELEMENTS
}  // namespace Discret

FOUR_C_NAMESPACE_CLOSE

#endif
