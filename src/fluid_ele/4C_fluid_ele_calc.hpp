/*----------------------------------------------------------------------*/
/*! \file

\brief main file containing routines for calculation of fluid element

\level 1


*/
/*----------------------------------------------------------------------*/

#ifndef FOUR_C_FLUID_ELE_CALC_HPP
#define FOUR_C_FLUID_ELE_CALC_HPP

#include "4C_config.hpp"

#include "4C_fluid_ele.hpp"
#include "4C_fluid_ele_interface.hpp"
#include "4C_inpar_fluid.hpp"

FOUR_C_NAMESPACE_OPEN

// forward declarations
namespace FLD
{
  template <Core::FE::CellType distype, int numdofpernode,
      Discret::ELEMENTS::Fluid::EnrichmentType enrtype = Discret::ELEMENTS::Fluid::none>
  class RotationallySymmetricPeriodicBC;

  class TDSEleData;
}  // namespace FLD

namespace Mat
{
  class Material;
}

namespace Discret
{
  namespace ELEMENTS
  {
    class FluidEleParameter;
    class FluidEleParameterTimInt;

    /// Fluid element implementation
    /*!
      This internal class keeps all the working arrays needed to
      calculate the fluid element. Additionally, the method Sysmat()
      provides a clean and fast element implementation.

      <h3>Purpose</h3>

      The idea is to separate the element maintenance (class Fluid) from the
      mathematical contents (this class). There are different
      implementations of the fluid element, this is just one such
      implementation.

      The fluid element will allocate exactly one object of this class for all
      fluid elements with the same number of nodes in the mesh. This
      allows us to use exactly matching working arrays (and keep them
      around.)

      The code is meant to be as clean as possible. This is the only way
      to keep it fast. The number of working arrays has to be reduced to
      a minimum so that the element fits into the cache. (There might be
      room for improvements.)

      <h3>Usability</h3>

      The calculations are done by the evaluate() method. There are two
      version. The virtual method that is inherited from FluidEleInterface
      (and called from Fluid) and the non-virtual one that does the actual
      work. The non-virtual evaluate() method must be callable without an actual
      Fluid object.

      \author u.kue
      \date 07/07
    */

    template <Core::FE::CellType distype,
        Discret::ELEMENTS::Fluid::EnrichmentType enrtype = Discret::ELEMENTS::Fluid::none>
    class FluidEleCalc : public FluidEleInterface
    {
     public:
      //! nen_: number of element nodes (T. Hughes: The Finite Element Method)
      static constexpr int nen_ =
          Discret::ELEMENTS::MultipleNumNode<enrtype>::multipleNode * Core::FE::num_nodes<distype>;

      //! number of space dimensions
      static constexpr int nsd_ = Core::FE::dim<distype>;

      static constexpr int numdofpernode_ = nsd_ + 1;

      virtual int integrate_shape_function(Discret::ELEMENTS::Fluid* ele,
          Core::FE::Discretization& discretization, const std::vector<int>& lm,
          Core::LinAlg::SerialDenseVector& elevec1);


      int integrate_shape_function(Discret::ELEMENTS::Fluid* ele,
          Core::FE::Discretization& discretization, const std::vector<int>& lm,
          Core::LinAlg::SerialDenseVector& elevec1,
          const Core::FE::GaussIntegration& intpoints) override;


      int integrate_shape_function_xfem(Discret::ELEMENTS::Fluid* ele,
          Core::FE::Discretization& discretization, const std::vector<int>& lm,
          Core::LinAlg::SerialDenseVector& elevec1,
          const std::vector<Core::FE::GaussIntegration>& intpoints,
          const Core::Geo::Cut::plain_volumecell_set& cells) override
      {
        FOUR_C_THROW("Implemented in derived xfem class!");
        return 1;
      };


      /// Evaluate supporting methods of the element
      /*!
        Interface function for supporting methods of the element
       */
      int EvaluateService(Discret::ELEMENTS::Fluid* ele, Teuchos::ParameterList& params,
          Teuchos::RCP<Core::Mat::Material>& mat, Core::FE::Discretization& discretization,
          std::vector<int>& lm, Core::LinAlg::SerialDenseMatrix& elemat1,
          Core::LinAlg::SerialDenseMatrix& elemat2, Core::LinAlg::SerialDenseVector& elevec1,
          Core::LinAlg::SerialDenseVector& elevec2,
          Core::LinAlg::SerialDenseVector& elevec3) override;

      /*! \brief Calculate a integrated divergence operator in vector form
       *
       *   The vector valued operator \f$B\f$ is constructed such that
       *   \f$\int_\Omega div (u) \,\mathrm{d}\Omega = B^T u = 0\f$
       *
       *   \author mayr.mt
       *   \date   04/2012
       */
      virtual int CalcDivOp(Discret::ELEMENTS::Fluid* ele,  //< current fluid element
          Core::FE::Discretization& discretization,         //< fluid discretization
          std::vector<int>& lm,                             //< some DOF management
          Core::LinAlg::SerialDenseVector& elevec1  //< reference to element vector to be filled
      );

      /*! \brief Calculate element mass matrix
       *
       *  \author mayr.mt \date 05/2014
       */
      virtual int CalcMassMatrix(Discret::ELEMENTS::Fluid* ele,
          //    Teuchos::ParameterList&              params,
          Core::FE::Discretization& discretization, const std::vector<int>& lm,
          Teuchos::RCP<Core::Mat::Material>& mat, Core::LinAlg::SerialDenseMatrix& elemat1_epetra
          //    const Core::FE::GaussIntegration & intpoints
      );

      /*! \brief Interpolate velocity gradient and pressure to given point
       *
       *  \author rauch \date 05/2014
       */
      int interpolate_velocity_gradient_and_pressure(Discret::ELEMENTS::Fluid* ele,
          Core::FE::Discretization& discretization, const std::vector<int>& lm,
          Core::LinAlg::SerialDenseVector& elevec1_epetra,
          Core::LinAlg::SerialDenseVector& elevec2_epetra);

      /*! \brief Interpolate velocity to given point
       *
       *  \author rauch \date 05/2014
       */
      int interpolate_velocity_to_node(Teuchos::ParameterList& params,
          Discret::ELEMENTS::Fluid* ele, Core::FE::Discretization& discretization,
          const std::vector<int>& lm, Core::LinAlg::SerialDenseVector& elevec1_epetra,
          Core::LinAlg::SerialDenseVector& elevec2_epetra);

      /*! \brief Interpolate velocity to given point
       *
       *  \author rauch \date 07/2015
       */
      int correct_immersed_bound_velocities(Teuchos::ParameterList& params,
          Discret::ELEMENTS::Fluid* ele, Core::FE::Discretization& discretization,
          const std::vector<int>& lm, Teuchos::RCP<Core::Mat::Material>& mat,
          Core::LinAlg::SerialDenseVector& elevec1_epetra,
          Core::LinAlg::SerialDenseVector& elevec2_epetra);

      /*---------------------------------------------------------------------*
       | Action type: interpolate_velocity_to_given_point                    |
       | calculate velocity at given point                       ghamm 12/15 |
       *---------------------------------------------------------------------*/
      int interpolate_velocity_to_point(Discret::ELEMENTS::Fluid* ele,
          Teuchos::ParameterList& params, Core::FE::Discretization& discretization,
          std::vector<int>& lm, Core::LinAlg::SerialDenseVector& elevec1,
          Core::LinAlg::SerialDenseVector& elevec2);

      /*! \brief Interpolate pressure to given point
       *
       *  \author ghamm \date 06/2015
       */
      int interpolate_pressure_to_point(Discret::ELEMENTS::Fluid* ele,
          Teuchos::ParameterList& params, Core::FE::Discretization& discretization,
          std::vector<int>& lm, Core::LinAlg::SerialDenseVector& elevec1_epetra);

      /*! \brief Reset debug output of immersed element
       *
       *  \author rauch \date 05/2014
       */
      int ResetImmersedEle(Discret::ELEMENTS::Fluid* ele, Teuchos::ParameterList& params);

      /*! \brief Calculate coordinates and velocities and element center
       *
       *  \author bk \date 01/2015
       */
      virtual int calc_vel_gradient_ele_center(Discret::ELEMENTS::Fluid* ele,
          Core::FE::Discretization& discretization, const std::vector<int>& lm,
          Core::LinAlg::SerialDenseVector& elevec1, Core::LinAlg::SerialDenseVector& elevec2);


      /*! \brief Calculate properties for adaptive time step based on CFL number
       *
       *  \author bk \date 08/2014
       */
      virtual int CalcTimeStep(Discret::ELEMENTS::Fluid* ele,
          Core::FE::Discretization& discretization, const std::vector<int>& lm,
          Core::LinAlg::SerialDenseVector& elevec1);

      /*! \brief Calculate channel statistics
       *
       *  \author bk \date 05/2014
       */
      virtual int calc_channel_statistics(Discret::ELEMENTS::Fluid* ele,
          Teuchos::ParameterList& params, Core::FE::Discretization& discretization,
          const std::vector<int>& lm, Teuchos::RCP<Core::Mat::Material>& mat);

      /*! \brief Calculate mass flow for periodic hill
       *
       *  \author bk \date 12/2014
       */
      virtual int calc_mass_flow_periodic_hill(Discret::ELEMENTS::Fluid* ele,
          Teuchos::ParameterList& params, Core::FE::Discretization& discretization,
          const std::vector<int>& lm, Core::LinAlg::SerialDenseVector& elevec1,
          Teuchos::RCP<Core::Mat::Material>& mat);

      /*! \brief Project velocity gradient to nodal level
       *
       *   \author ghamm
       *   \date   06/2014
       */
      virtual int vel_gradient_projection(Discret::ELEMENTS::Fluid* ele,
          Teuchos::ParameterList& params, Core::FE::Discretization& discretization,
          std::vector<int>& lm, Core::LinAlg::SerialDenseMatrix& elemat1,
          Core::LinAlg::SerialDenseMatrix& elemat2);


      /*! \brief Project pressure gradient to nodal level
       *
       *   \author mwinter
       *   \date   09/2015
       */
      virtual int pres_gradient_projection(Discret::ELEMENTS::Fluid* ele,
          Teuchos::ParameterList& params, Core::FE::Discretization& discretization,
          std::vector<int>& lm, Core::LinAlg::SerialDenseMatrix& elemat1,
          Core::LinAlg::SerialDenseMatrix& elemat2);

      /*! \brief Calculate a divergence of velocity at the element center
       *
       *   \author ehrl
       *   \date   12/2012
       */
      virtual int ComputeDivU(Discret::ELEMENTS::Fluid* ele,  //< current fluid element
          Core::FE::Discretization& discretization,           //< fluid discretization
          std::vector<int>& lm,                               //< location vector for DOF management
          Core::LinAlg::SerialDenseVector& elevec1  //< reference to element vector to be filled
      );

      /// Evaluate element ERROR
      /*!
          general function to compute the error (analytical solution) for particular problem type
       */
      virtual int compute_error(Discret::ELEMENTS::Fluid* ele, Teuchos::ParameterList& params,
          Teuchos::RCP<Core::Mat::Material>& mat, Core::FE::Discretization& discretization,
          std::vector<int>& lm, Core::LinAlg::SerialDenseVector& elevec);

      int compute_error(Discret::ELEMENTS::Fluid* ele, Teuchos::ParameterList& params,
          Teuchos::RCP<Core::Mat::Material>& mat, Core::FE::Discretization& discretization,
          std::vector<int>& lm, Core::LinAlg::SerialDenseVector& elevec1,
          const Core::FE::GaussIntegration& intpoints2) override;

      /*!
       \brief Evaluates the analytic solution in the given point
       */
      static void evaluate_analytic_solution_point(const Core::LinAlg::Matrix<nsd_, 1>& xyzint,
          const double t, const Inpar::FLUID::CalcError calcerr, const int calcerrfunctno,
          const Teuchos::RCP<Core::Mat::Material>& mat, Core::LinAlg::Matrix<nsd_, 1>& u, double& p,
          Core::LinAlg::Matrix<nsd_, nsd_>& dervel, bool isFullImplPressure = false,
          double deltat = 0.0);

      /// Evaluate the element
      /*!
        Generic virtual interface function. Called via base pointer.
       */
      int evaluate(Discret::ELEMENTS::Fluid* ele, Core::FE::Discretization& discretization,
          const std::vector<int>& lm, Teuchos::ParameterList& params,
          Teuchos::RCP<Core::Mat::Material>& mat, Core::LinAlg::SerialDenseMatrix& elemat1_epetra,
          Core::LinAlg::SerialDenseMatrix& elemat2_epetra,
          Core::LinAlg::SerialDenseVector& elevec1_epetra,
          Core::LinAlg::SerialDenseVector& elevec2_epetra,
          Core::LinAlg::SerialDenseVector& elevec3_epetra, bool offdiag = false) override;

      /// Evaluate the element at specified gauss points
      int evaluate(Discret::ELEMENTS::Fluid* ele, Core::FE::Discretization& discretization,
          const std::vector<int>& lm, Teuchos::ParameterList& params,
          Teuchos::RCP<Core::Mat::Material>& mat, Core::LinAlg::SerialDenseMatrix& elemat1_epetra,
          Core::LinAlg::SerialDenseMatrix& elemat2_epetra,
          Core::LinAlg::SerialDenseVector& elevec1_epetra,
          Core::LinAlg::SerialDenseVector& elevec2_epetra,
          Core::LinAlg::SerialDenseVector& elevec3_epetra,
          const Core::FE::GaussIntegration& intpoints, bool offdiag = false) override;

      int compute_error_interface(Discret::ELEMENTS::Fluid* ele,     ///< fluid element
          Core::FE::Discretization& dis,                             ///< background discretization
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
          ) override
      {
        FOUR_C_THROW("Implemented in derived xfem class!");
        return 1;
      };

      /// Evaluate the XFEM cut element
      int EvaluateXFEM(Discret::ELEMENTS::Fluid* ele, Core::FE::Discretization& discretization,
          const std::vector<int>& lm, Teuchos::ParameterList& params,
          Teuchos::RCP<Core::Mat::Material>& mat, Core::LinAlg::SerialDenseMatrix& elemat1_epetra,
          Core::LinAlg::SerialDenseMatrix& elemat2_epetra,
          Core::LinAlg::SerialDenseVector& elevec1_epetra,
          Core::LinAlg::SerialDenseVector& elevec2_epetra,
          Core::LinAlg::SerialDenseVector& elevec3_epetra,
          const std::vector<Core::FE::GaussIntegration>& intpoints,
          const Core::Geo::Cut::plain_volumecell_set& cells, bool offdiag = false) override
      {
        FOUR_C_THROW("Implemented in derived xfem class!");
        return 1;
      }

      /*!
        \brief calculate dissipation of various terms (evaluation of turbulence models)
      */
      virtual int calc_dissipation(Fluid* ele, Teuchos::ParameterList& params,
          Core::FE::Discretization& discretization, std::vector<int>& lm,
          Teuchos::RCP<Core::Mat::Material> mat);

      /*!
        \brief finite difference check for debugging
      */
      virtual void FDcheck(const Core::LinAlg::Matrix<nsd_, nen_>& evelaf,
          const Core::LinAlg::Matrix<nsd_, nen_>& eveln,
          const Core::LinAlg::Matrix<nsd_, nen_>& fsevelaf,
          const Core::LinAlg::Matrix<nen_, 1>& epreaf,
          const Core::LinAlg::Matrix<nsd_, nen_>& eaccam,
          const Core::LinAlg::Matrix<nen_, 1>& escaaf, const Core::LinAlg::Matrix<nen_, 1>& escaam,
          const Core::LinAlg::Matrix<nen_, 1>& escadtam,
          const Core::LinAlg::Matrix<nsd_, nen_>& emhist,
          const Core::LinAlg::Matrix<nsd_, nen_>& edispnp,
          const Core::LinAlg::Matrix<nsd_, nen_>& egridv,
          const Core::LinAlg::Matrix<(nsd_ + 1) * nen_, (nsd_ + 1) * nen_>& estif,
          const Core::LinAlg::Matrix<(nsd_ + 1) * nen_, (nsd_ + 1) * nen_>& emesh,
          const Core::LinAlg::Matrix<(nsd_ + 1) * nen_, 1>& eforce, const double thermpressaf,
          const double thermpressam, const double thermpressdtaf, const double thermpressdtam,
          const Teuchos::RCP<const Core::Mat::Material> material, const double timefac,
          const double& Cs, const double& Cs_delta_sq, const double& l_tau);


      void element_xfem_interface_hybrid_lm(Discret::ELEMENTS::Fluid* ele,  ///< fluid element
          Core::FE::Discretization& dis,                             ///< background discretization
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
          ) override
      {
        FOUR_C_THROW("Implemented in derived xfem class!");
        return;
      }


      void element_xfem_interface_nit(Discret::ELEMENTS::Fluid* ele,  ///< fluid element
          Core::FE::Discretization& dis,                              ///< background discretization
          const std::vector<int>& lm,                                 ///< element local map
          const Teuchos::RCP<XFEM::ConditionManager>& cond_manager,   ///< XFEM condition manager
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
          ) override
      {
        FOUR_C_THROW("Implemented in derived xfem class!");
        return;
      }

      void calculate_continuity_xfem(Discret::ELEMENTS::Fluid* ele, Core::FE::Discretization& dis,
          const std::vector<int>& lm, Core::LinAlg::SerialDenseVector& elevec1_epetra,
          const Core::FE::GaussIntegration& intpoints) override
      {
        FOUR_C_THROW("Implemented in derived xfem class!");
        return;
      }

      void calculate_continuity_xfem(Discret::ELEMENTS::Fluid* ele, Core::FE::Discretization& dis,
          const std::vector<int>& lm, Core::LinAlg::SerialDenseVector& elevec1_epetra) override
      {
        FOUR_C_THROW("Implemented in derived xfem class!");
        return;
      }

      /// calculate body force from nodal conditions. Static function interface to allow
      /// for its use in FluidEleCalcHDG.
      static void BodyForce(Discret::ELEMENTS::Fluid* ele,  //< pointer to element
          const double time,                                //< current time
          const Inpar::FLUID::PhysicalType physicaltype,    //< physical type
          Core::LinAlg::Matrix<nsd_, nen_>& ebofoaf,        //< body force at nodes
          Core::LinAlg::Matrix<nsd_, nen_>&
              eprescpgaf,  //< prescribed pressure gradient (required for turbulent channel flow!)
          Core::LinAlg::Matrix<nen_, 1>& escabofoaf  //< scatra body force at nodes
      );

      /// calculate correction term at nodes
      static void CorrectionTerm(Discret::ELEMENTS::Fluid* ele,  //< pointer to element
          Core::LinAlg::Matrix<1, nen_>& ecorrectionterm         //<correction term at nodes
      );


     protected:
      /// private Constructor since we are a Singleton.
      FluidEleCalc();

      /*!
        \brief evaluate function for fluid element

        Specific evaluate function without any knowledge about DRT objects. This
        way the element evaluation is independent of the specific mesh storage.
       */
      virtual int evaluate(Teuchos::ParameterList& params,
          const Core::LinAlg::Matrix<nsd_, nen_>& ebofoaf,
          const Core::LinAlg::Matrix<nsd_, nen_>& eprescpgaf,
          const Core::LinAlg::Matrix<nsd_, nen_>& ebofon,
          const Core::LinAlg::Matrix<nsd_, nen_>& eprescpgn,
          Core::LinAlg::Matrix<(nsd_ + 1) * nen_, (nsd_ + 1) * nen_>& elemat1,
          Core::LinAlg::Matrix<(nsd_ + 1) * nen_, (nsd_ + 1) * nen_>& elemat2,
          Core::LinAlg::Matrix<(nsd_ + 1) * nen_, 1>& elevec1,
          const Core::LinAlg::Matrix<nsd_, nen_>& evelaf,
          const Core::LinAlg::Matrix<nen_, 1>& epreaf,
          const Core::LinAlg::Matrix<nsd_, nen_>& evelam,
          const Core::LinAlg::Matrix<nen_, 1>& epream, const Core::LinAlg::Matrix<nen_, 1>& eprenp,
          const Core::LinAlg::Matrix<nsd_, nen_>& evelnp,
          const Core::LinAlg::Matrix<nen_, 1>& escaaf,
          const Core::LinAlg::Matrix<nsd_, nen_>& emhist,
          const Core::LinAlg::Matrix<nsd_, nen_>& eaccam,
          const Core::LinAlg::Matrix<nen_, 1>& escadtam,
          const Core::LinAlg::Matrix<nsd_, nen_>& eveldtam,
          const Core::LinAlg::Matrix<nen_, 1>& epredtam,
          const Core::LinAlg::Matrix<nen_, 1>& escabofoaf,
          const Core::LinAlg::Matrix<nen_, 1>& escabofon,
          const Core::LinAlg::Matrix<nsd_, nen_>& eveln, const Core::LinAlg::Matrix<nen_, 1>& epren,
          const Core::LinAlg::Matrix<nen_, 1>& escaam,
          const Core::LinAlg::Matrix<nsd_, nen_>& edispnp,
          const Core::LinAlg::Matrix<nsd_, nen_>& egridv,
          const Core::LinAlg::Matrix<nsd_, nen_>& egridvn,
          const Core::LinAlg::Matrix<nsd_, nen_>& fsevelaf,
          const Core::LinAlg::Matrix<nen_, 1>& fsescaaf,
          const Core::LinAlg::Matrix<nsd_, nen_>& evel_hat,
          const Core::LinAlg::Matrix<nsd_ * nsd_, nen_>& ereynoldsstress_hat,
          const Core::LinAlg::Matrix<nen_, 1>& eporo,
          const Core::LinAlg::Matrix<nsd_, 2 * nen_>& egradphi,
          const Core::LinAlg::Matrix<nen_, 2 * 1>& ecurvature,
          Teuchos::RCP<Core::Mat::Material> mat, bool isale, bool isowned, double CsDeltaSq,
          double CiDeltaSq, double* saccn, double* sveln, double* svelnp,
          const Core::FE::GaussIntegration& intpoints, bool offdiag);

      /*!
        \brief calculate element matrix and rhs

        @param ebofoaf           (i) body force at n+alpha_F/n+1
        @param eprescpgaf        (i) prescribed pressure gradient at n+alpha_F/n+1 (required for
        turbulent channel flow)
        @param ebofon            (i) body force at n
        @param eprescpgaf        (i) prescribed pressure gradient at n (required for turbulent
        channel flow)
        @param evelaf            (i) nodal velocities at n+alpha_F/n+1
        @param evelam           (i) nodal velocities at n+alpha_M/n
        @param eveln            (i) nodal velocities at n
        @param evelnp           (i) nodal velocities at n+1 (np_genalpha)
        @param fsevelaf         (i) fine-scale nodal velocities at n+alpha_F/n+1
        @param fsescaaf         (i) fine-scale nodal scalar at n+alpha_F/n+1
        @param epreaf           (i) nodal pressure at n+alpha_F/n+1
        @param epream           (i) nodal pressure at n+alpha_M/n
        @param eprenp           (i) nodal pressure at n+1
        @param eaccam           (i) nodal accelerations at n+alpha_M
        @param escaaf           (i) nodal scalar at n+alpha_F/n+1
        @param escaam           (i) nodal scalar at n+alpha_M/n
        @param escadtam         (i) nodal scalar derivatives at n+alpha_M/n+1
        @param eveldtam         (i) nodal velocity derivatives at n+alpha_M/n+1
        @param epredtam         (i) nodal pressure derivatives at n+alpha_M/n+1
        @param emhist           (i) time rhs for momentum equation
        @param edispnp          (i) nodal displacements (on moving mesh)
        @param egridv           (i) grid velocity (on moving mesh)
        @param estif            (o) element matrix to calculate
        @param emesh            (o) linearization wrt mesh motion
        @param eforce           (o) element rhs to calculate
        @param egradphi         (i) gradient of nodal scalar at nodes
        @param ecurvature       (i) curvature of scalar at nodes
        @param thermpressaf     (i) thermodynamic pressure at n+alpha_F/n+1
        @param thermpressam     (i) thermodynamic pressure at n+alpha_M/n
        @param thermpressdtaf   (i) thermodynamic pressure derivative at n+alpha_F/n+1
        @param thermpressdtam   (i) thermodynamic pressure derivative at n+alpha_M/n+1
        @param material         (i) fluid material
        @param Cs_delta_sq      (i) parameter for dynamic Smagorinsky model (Cs*h*h)
        @param isale            (i) ALE flag
        @param intpoints        (i) Gaussian integration points

        */
      virtual void sysmat(const Core::LinAlg::Matrix<nsd_, nen_>& ebofoaf,
          const Core::LinAlg::Matrix<nsd_, nen_>& eprescpgaf,
          const Core::LinAlg::Matrix<nsd_, nen_>& ebofon,
          const Core::LinAlg::Matrix<nsd_, nen_>& eprescpgn,
          const Core::LinAlg::Matrix<nsd_, nen_>& evelaf,
          const Core::LinAlg::Matrix<nsd_, nen_>& evelam,
          const Core::LinAlg::Matrix<nsd_, nen_>& eveln,
          const Core::LinAlg::Matrix<nsd_, nen_>& evelnp,
          const Core::LinAlg::Matrix<nsd_, nen_>& fsevelaf,
          const Core::LinAlg::Matrix<nen_, 1>& fsescaaf,
          const Core::LinAlg::Matrix<nsd_, nen_>& evel_hat,
          const Core::LinAlg::Matrix<nsd_ * nsd_, nen_>& ereynoldsstress_hat,
          const Core::LinAlg::Matrix<nen_, 1>& epreaf, const Core::LinAlg::Matrix<nen_, 1>& epream,
          const Core::LinAlg::Matrix<nen_, 1>& epren, const Core::LinAlg::Matrix<nen_, 1>& eprenp,
          const Core::LinAlg::Matrix<nsd_, nen_>& eaccam,
          const Core::LinAlg::Matrix<nen_, 1>& escaaf, const Core::LinAlg::Matrix<nen_, 1>& escaam,
          const Core::LinAlg::Matrix<nen_, 1>& escadtam,
          const Core::LinAlg::Matrix<nsd_, nen_>& eveldtam,
          const Core::LinAlg::Matrix<nen_, 1>& epredtam,
          const Core::LinAlg::Matrix<nen_, 1>& escabofoaf,
          const Core::LinAlg::Matrix<nen_, 1>& escabofon,
          const Core::LinAlg::Matrix<nsd_, nen_>& emhist,
          const Core::LinAlg::Matrix<nsd_, nen_>& edispnp,
          const Core::LinAlg::Matrix<nsd_, nen_>& egridv,
          Core::LinAlg::Matrix<(nsd_ + 1) * nen_, (nsd_ + 1) * nen_>& estif,
          Core::LinAlg::Matrix<(nsd_ + 1) * nen_, (nsd_ + 1) * nen_>& emesh,
          Core::LinAlg::Matrix<(nsd_ + 1) * nen_, 1>& eforce,
          const Core::LinAlg::Matrix<nen_, 1>& eporo,
          const Core::LinAlg::Matrix<nsd_, 2 * nen_>& egradphi,
          const Core::LinAlg::Matrix<nen_, 2 * 1>& ecurvature, const double thermpressaf,
          const double thermpressam, const double thermpressdtaf, const double thermpressdtam,
          Teuchos::RCP<const Core::Mat::Material> material, double& Cs_delta_sq,
          double& Ci_delta_sq, double& Cv, bool isale, double* saccn, double* sveln, double* svelnp,
          const Core::FE::GaussIntegration& intpoints);


      virtual void sysmat_ost_new(const Core::LinAlg::Matrix<nsd_, nen_>& ebofoaf,
          const Core::LinAlg::Matrix<nsd_, nen_>& eprescpgaf,
          const Core::LinAlg::Matrix<nsd_, nen_>& ebofon,
          const Core::LinAlg::Matrix<nsd_, nen_>& eprescpgn,
          const Core::LinAlg::Matrix<nsd_, nen_>& evelaf,
          const Core::LinAlg::Matrix<nsd_, nen_>& evelam,
          const Core::LinAlg::Matrix<nsd_, nen_>& eveln,
          const Core::LinAlg::Matrix<nsd_, nen_>& evelnp,
          const Core::LinAlg::Matrix<nsd_, nen_>& fsevelaf,
          const Core::LinAlg::Matrix<nen_, 1>& fsescaaf,
          const Core::LinAlg::Matrix<nsd_, nen_>& evel_hat,
          const Core::LinAlg::Matrix<nsd_ * nsd_, nen_>& ereynoldsstress_hat,
          const Core::LinAlg::Matrix<nen_, 1>& epreaf, const Core::LinAlg::Matrix<nen_, 1>& epream,
          const Core::LinAlg::Matrix<nen_, 1>& epren, const Core::LinAlg::Matrix<nen_, 1>& eprenp,
          const Core::LinAlg::Matrix<nsd_, nen_>& eaccam,
          const Core::LinAlg::Matrix<nen_, 1>& escaaf, const Core::LinAlg::Matrix<nen_, 1>& escaam,
          const Core::LinAlg::Matrix<nen_, 1>& escadtam,
          const Core::LinAlg::Matrix<nen_, 1>& escabofoaf,
          const Core::LinAlg::Matrix<nen_, 1>& escabofon,
          const Core::LinAlg::Matrix<nsd_, nen_>& emhist,
          const Core::LinAlg::Matrix<nsd_, nen_>& edispnp,
          const Core::LinAlg::Matrix<nsd_, nen_>& egridv,
          const Core::LinAlg::Matrix<nsd_, nen_>& egridvn,
          Core::LinAlg::Matrix<(nsd_ + 1) * nen_, (nsd_ + 1) * nen_>& estif,
          Core::LinAlg::Matrix<(nsd_ + 1) * nen_, (nsd_ + 1) * nen_>& emesh,
          Core::LinAlg::Matrix<(nsd_ + 1) * nen_, 1>& eforce,
          const Core::LinAlg::Matrix<nen_, 1>& eporo,
          const Core::LinAlg::Matrix<nsd_, 2 * nen_>& egradphi,
          const Core::LinAlg::Matrix<nen_, 2 * 1>& ecurvature, const double thermpressaf,
          const double thermpressam, const double thermpressdtaf, const double thermpressdtam,
          Teuchos::RCP<const Core::Mat::Material> material, double& Cs_delta_sq,
          double& Ci_delta_sq, double& Cv, bool isale, double* saccn, double* sveln, double* svelnp,
          const Core::FE::GaussIntegration& intpoints);



      //! number of components necessary to store second derivatives
      /*!
       1 component  for nsd=1:  (N,xx)

       3 components for nsd=2:  (N,xx ; N,yy ; N,xy)

       6 components for nsd=3:  (N,xx ; N,yy ; N,zz ; N,xy ; N,xz ; N,yz)
      */
      static constexpr int numderiv2_ = Core::FE::DisTypeToNumDeriv2<distype>::numderiv2;

      //! calculate body force from nodal conditions
      void body_force(Discret::ELEMENTS::Fluid* ele,  //< pointer to element
          Core::LinAlg::Matrix<nsd_, nen_>& ebofoaf,  //< body force at nodes
          Core::LinAlg::Matrix<nsd_, nen_>&
              eprescpgaf,  //< prescribed pressure gradient (required for turbulent channel flow!)
          Core::LinAlg::Matrix<nen_, 1>& escabofoaf  //< scatra body force at nodes
      );

      //! calculate body force contribution from surface tension
      void add_surface_tension_force(
          const Core::LinAlg::Matrix<nen_, 1>& escaaf,  ///< scalar at time n+alpha_f / n+1
          const Core::LinAlg::Matrix<nen_, 1>& escaam,  ///< scalar at time n+alpha_m / n
          const Core::LinAlg::Matrix<nsd_, 2 * nen_>&
              egradphi,  //<gradient of scalar function at nodes
          const Core::LinAlg::Matrix<nen_, 2 * 1>&
              ecurvature  //<curvature of scalar function at nodes
      );


      //! evaluate shape functions and their derivatives at element center
      virtual void eval_shape_func_and_derivs_at_ele_center();

      //! brief evaluate shape functions and their derivatives at integration point
      virtual void eval_shape_func_and_derivs_at_int_point(
          const double* gpcoords,  ///< actual integration point (coords)
          double gpweight          ///< actual integration point (weight)
      );

      //! get ALE grid displacements and grid velocity for element
      void get_grid_disp_vel_ale(Core::FE::Discretization& discretization,
          const std::vector<int>& lm, Core::LinAlg::Matrix<nsd_, nen_>& edispnp,
          Core::LinAlg::Matrix<nsd_, nen_>& egridv);

      //! get ALE grid displacements and grid velocity for element
      void get_grid_disp_vel_aleost_new(Core::FE::Discretization& discretization,
          const std::vector<int>& lm, Core::LinAlg::Matrix<nsd_, nen_>& edispnp,
          Core::LinAlg::Matrix<nsd_, nen_>& egridvnp, Core::LinAlg::Matrix<nsd_, nen_>& egridvn);

      //! get ALE grid displacements for element
      virtual void get_grid_disp_ale(Core::FE::Discretization& discretization,
          const std::vector<int>& lm, Core::LinAlg::Matrix<nsd_, nen_>& edispnp);

      //! set the (relative) convective velocity at integration point for various physical types
      void set_convective_velint(const bool isale);

      //! set the (relative) convective velocity at integration point for various physical types
      void set_convective_velint_n(const bool isale);

      //! set element advective field for Oseen problems
      void set_advective_vel_oseen(Discret::ELEMENTS::Fluid* ele);

      //! get material parameters
      void get_material_params(
          Teuchos::RCP<const Core::Mat::Material> material,  ///< reference pointer to material
          const Core::LinAlg::Matrix<nsd_, nen_>& evelaf,    ///< velocity at time n+alpha_f / n+1
          const Core::LinAlg::Matrix<nen_, 1>& epreaf,       ///< pressure at time n+alpha_f / n+1
          const Core::LinAlg::Matrix<nen_, 1>& epream,       ///< pressure at time n+alpha_m / n
          const Core::LinAlg::Matrix<nen_, 1>& escaaf,       ///< scalar at time n+alpha_f / n+1
          const Core::LinAlg::Matrix<nen_, 1>& escaam,       ///< scalar at time n+alpha_m / n
          const Core::LinAlg::Matrix<nen_, 1>&
              escabofoaf,               ///< body force for scalar transport at time n+alpha_f / n+1
          const double thermpressaf,    ///< thermodynamic pressure at time n+alpha_f / n+1
          const double thermpressam,    ///< thermodynamic pressure at time n+alpha_m / n
          const double thermpressdtaf,  ///< time derivative of thermodynamic pressure at time
                                        ///< n+alpha_f / n+1
          const double thermpressdtam,  ///< time derivative of thermodynamic pressure at time
                                        ///< n+alpha_m / n+1
          const double vol              ///< element volume
      );

      //! get material parameters
      void get_material_params(Teuchos::RCP<const Core::Mat::Material> material,
          const Core::LinAlg::Matrix<nsd_, nen_>& evelaf,
          const Core::LinAlg::Matrix<nen_, 1>& epreaf, const Core::LinAlg::Matrix<nen_, 1>& epream,
          const Core::LinAlg::Matrix<nen_, 1>& escaaf, const Core::LinAlg::Matrix<nen_, 1>& escaam,
          const Core::LinAlg::Matrix<nen_, 1>& escabofoaf, const double thermpressaf,
          const double thermpressam, const double thermpressdtaf, const double thermpressdtam,
          const double vol, double& densam, double& densaf, double& densn, double& visc,
          double& viscn, double& gamma);

      //! return constant mk for stabilization parameters
      virtual double get_mk();

      //! calculate stabilization parameter
      void calc_stab_parameter(const double vol);  ///< volume

      //! calculate characteristic element length
      void calc_char_ele_length(const double vol,  ///< volume
          const double vel_norm,                   ///< norm of velocity vector
          double& h_u,                             ///< length for tau_Mu
          double& h_p);                            ///< length for tau_Mp/tau_C


      //! calculate div(epsilon(u))
      void calc_div_eps(
          const Core::LinAlg::Matrix<nsd_, nen_>& evelaf,  ///< velocity at time n+alpha_f / n+1
          const Core::LinAlg::Matrix<nsd_, nen_>& eveln);  ///< velocity at time n+alpha_f / n+1

      //! compute residual of momentum equation and subgrid-scale velocity
      virtual void compute_subgrid_scale_velocity(
          const Core::LinAlg::Matrix<nsd_, nen_>& eaccam,  ///< acceleration at time n+alpha_M
          double& fac1,                                    ///< factor for old s.-s. velocities
          double& fac2,                                    ///< factor for old s.-s. accelerations
          double& fac3,     ///< factor for residual in current s.-s. velocities
          double& facMtau,  ///< facMtau = modified tau_M (see code)
          int iquad,        ///< integration point
          double* saccn,    ///< s.-s. acceleration at time n+alpha_a / n
          double* sveln,    ///< s.-s. velocity at time n+alpha_a / n
          double* svelnp    ///< s.-s. velocity at time n+alpha_f / n+1
      );


      //! compute residual of momentum equation and subgrid-scale velocity
      void compute_subgrid_scale_velocity_ost_new(
          const Core::LinAlg::Matrix<nsd_, nen_>& eaccam,  ///< acceleration at time n+alpha_M
          double& fac1,                                    ///< factor for old s.-s. velocities
          double& fac2,                                    ///< factor for old s.-s. accelerations
          double& fac3,     ///< factor for residual in current s.-s. velocities
          double& facMtau,  ///< facMtau = modified tau_M (see code)
          int iquad,        ///< integration point
          double* saccn,    ///< s.-s. acceleration at time n+alpha_a / n
          double* sveln,    ///< s.-s. velocity at time n+alpha_a / n
          double* svelnp    ///< s.-s. velocity at time n+alpha_f / n+1
      );

      //! Provide linearization of Garlerkin momentum residual with respect to the velocities
      virtual void lin_gal_mom_res_u(
          Core::LinAlg::Matrix<nsd_ * nsd_, nen_>&
              lin_resM_Du,          ///< linearisation of the Garlerkin momentum residual
          const double& timefacfac  ///< = timefac x fac
      );

      //! Provide linearization of Garlerkin momentum residual with respect to the velocities
      void lin_gal_mom_res_uost_new(
          Core::LinAlg::Matrix<nsd_ * nsd_, nen_>&
              lin_resM_Du,          ///< linearisation of the Garlerkin momentum residual
          const double& timefacfac  ///< = timefac x fac
      );

      //! Provide linearization of Garlerkin momentum residual with respect to the velocities in the
      //! case if subscales
      void lin_gal_mom_res_u_subscales(Core::LinAlg::Matrix<nen_ * nsd_, nen_>&
                                           estif_p_v,  ///< block (weighting function v x pressure)
          Core::LinAlg::Matrix<nsd_ * nsd_, nen_>&
              lin_resM_Du,  ///< linearisation of the Garlerkin momentum residual
          Core::LinAlg::Matrix<nsd_, 1>& resM_Du,  ///< residual of the fluid momentum equation
          const double& timefacfac,                ///< (time factor) x (integration factor)
          const double& facMtau                    ///< facMtau = modified tau_M (see code)
      );

      //! Compute element matrix and rhs entries: inertia, convective andyn
      //! reactive terms of the Galerkin part
      void inertia_convection_reaction_gal_part(Core::LinAlg::Matrix<nen_ * nsd_, nen_ * nsd_>&
                                                    estif_u,  ///< block (weighting function v x u)
          Core::LinAlg::Matrix<nsd_, nen_>& velforce,         ///< rhs forces velocity
          Core::LinAlg::Matrix<nsd_ * nsd_, nen_>&
              lin_resM_Du,  ///< linearisation of the Garlerkin momentum residual
          Core::LinAlg::Matrix<nsd_, 1>&
              resM_Du,           ///< linearisation of the Garlerkin momentum residual
          const double& rhsfac,  ///< right-hand-side factor
          const double& rhsfacn  ///< right-hand-side factor time step n
      );

      //! Compute element matrix entries: for the viscous terms of the Galerkin part
      void viscous_gal_part(Core::LinAlg::Matrix<nen_ * nsd_, nen_ * nsd_>&
                                estif_u,                 ///< block (weighting function v x u)
          Core::LinAlg::Matrix<nsd_, nen_>& velforce,    ///< rhs forces velocity
          Core::LinAlg::Matrix<nsd_, nsd_>& viscstress,  ///< viscous stresses
          const double& timefacfac,                      ///< = timefac x fac
          const double& rhsfac,                          ///< right-hand-side factor
          const double& rhsfacn                          ///< right-hand-side factor time step n
      );

      //! Compute element matrix entries: div-grad stabilization and the rhs of the viscous term
      void cont_stab(Core::LinAlg::Matrix<nen_ * nsd_, nen_ * nsd_>&
                         estif_u,                      ///< block (weighting function v x u)
          Core::LinAlg::Matrix<nsd_, nen_>& velforce,  ///< rhs forces velocity
          const double& timefac,                       ///< time factor
          const double& timefacfac,                    ///< = timefac x fac
          const double& timefacfacpre,                 ///< = timefac x fac
          const double& rhsfac,                        ///< right-hand-side factor
          const double& rhsfacn                        ///< right-hand-side factor time step n
      );

      //! Compute element matrix entries: pressure terms of the Garlerkin part and rhs
      void pressure_gal_part(Core::LinAlg::Matrix<nen_ * nsd_, nen_>&
                                 estif_p_v,            ///< block (weighting function v x pressure)
          Core::LinAlg::Matrix<nsd_, nen_>& velforce,  ///< rhs forces velocity
          const double& timefacfac,                    ///< = timefac x fac
          const double& timefacfacpre,
          const double& rhsfac,   ///< right-hand-side factor
          const double& rhsfacn,  ///< right-hand-side factor time step n
          const double& press,    ///< pressure at integration point
          const double& pressn    ///< pressure at integration point
      );

      //! Compute element matrix entries: continuity terms of the Garlerkin part and rhs
      void continuity_gal_part(
          Core::LinAlg::Matrix<nen_, nen_ * nsd_>& estif_q_u,  ///< block (weighting function q x u)
          Core::LinAlg::Matrix<nen_, 1>& preforce,             ///< rhs forces pressure
          const double& timefacfac,                            ///< = timefac x fac
          const double& timefacfacpre,
          const double& rhsfac,  ///< right-hand-side factor
          const double& rhsfacn  ///< right-hand-side factor at time step n
      );

      //! Compute element matrix entries: pressure projection terms
      void pressure_projection(Core::LinAlg::Matrix<nen_, nen_>& ppmat);

      //! Finalize pressure projection terms
      void pressure_projection_finalize(Core::LinAlg::Matrix<nen_, nen_>& ppmat,
          Core::LinAlg::Matrix<nen_, 1>& preforce, const Core::LinAlg::Matrix<nen_, 1>& eprenp);

      //! Compute element matrix entries: body force terms on rhs
      void body_force_rhs_term(Core::LinAlg::Matrix<nsd_, nen_>& velforce,  ///< rhs forces velocity
          const double& rhsfac,  ///< right-hand-side factor for residuals
          const double rhsfacn   //= 0.0 ///< right-hand-side factor for residuals
      );

      //! Compute element matrix entries: conservative formulation
      virtual void conservative_formulation(Core::LinAlg::Matrix<nen_ * nsd_, nen_ * nsd_>&
                                                estif_u,  ///< block (weighting function v x u)
          Core::LinAlg::Matrix<nsd_, nen_>& velforce,     ///< rhs forces velocity
          const double& timefacfac,                       ///< = timefac x fac
          const double& rhsfac                            ///< right-hand-side factor
      );

      //! Provide linearization of stabilization residual with respect to the velocities
      void stab_lin_gal_mom_res_u(Core::LinAlg::Matrix<nsd_ * nsd_, nen_>&
                                      lin_resM_Du,  ///< linearisation of the stabilization residual
          const double& timefacfac                  ///< = timefac x fac
      );

      //! Compute element matrix entries: PSPG
      void pspg(
          Core::LinAlg::Matrix<nen_, nen_ * nsd_>& estif_q_u,  ///< block (weighting function q x u)
          Core::LinAlg::Matrix<nen_, nen_>& ppmat,             ///< block (weighting function q x p)
          Core::LinAlg::Matrix<nen_, 1>& preforce,             ///< rhs forces pressure
          Core::LinAlg::Matrix<nsd_ * nsd_, nen_>&
              lin_resM_Du,           ///< linearisation of the stabilization residual
          const double& fac3,        ///< factor for residual in current subgrid velocities
          const double& timefacfac,  ///< = timefac x fac
          const double& timefacfacpre,
          const double& rhsfac,  ///< right-hand-side factor for residuals
          const int iquad        ///< index of current integration point
      );

      //! Compute element matrix entries: PSPG
      void pspgost_new(
          Core::LinAlg::Matrix<nen_, nen_ * nsd_>& estif_q_u,  ///< block (weighting function q x u)
          Core::LinAlg::Matrix<nen_, nen_>& ppmat,             ///< block (weighting function q x p)
          Core::LinAlg::Matrix<nen_, 1>& preforce,             ///< rhs forces pressure
          Core::LinAlg::Matrix<nsd_ * nsd_, nen_>&
              lin_resM_Du,           ///< linearisation of the stabilization residual
          const double& fac3,        ///< factor for residual in current subgrid velocities
          const double& timefacfac,  ///< = timefac x fac
          const double& timefacfacpre,
          const double& rhsfac,  ///< right-hand-side factor for residuals
          const int iquad        ///< index of current integration point
      );

     protected:
      //! Compute element matrix entries: SUPG
      void supg(Core::LinAlg::Matrix<nen_ * nsd_, nen_ * nsd_>&
                    estif_u,                                   ///< block (weighting function v x u)
          Core::LinAlg::Matrix<nen_ * nsd_, nen_>& estif_p_v,  ///< block (weighting function v x p)
          Core::LinAlg::Matrix<nsd_, nen_>& velforce,          ///< rhs forces velocity
          Core::LinAlg::Matrix<nen_, 1>& preforce,             ///< rhs forces pressure
          Core::LinAlg::Matrix<nsd_ * nsd_, nen_>&
              lin_resM_Du,           ///< linearisation of the stabilization residual
          const double& fac3,        ///< factor for residual in current subgrid velocities
          const double& timefacfac,  ///< = timefac x fac
          const double& timefacfacpre,
          const double& rhsfac  ///< right-hand-side factor for residuals
      );

      //! Compute element matrix entries: SUPG
      void supgost_new(Core::LinAlg::Matrix<nen_ * nsd_, nen_ * nsd_>&
                           estif_u,                            ///< block (weighting function v x u)
          Core::LinAlg::Matrix<nen_ * nsd_, nen_>& estif_p_v,  ///< block (weighting function v x p)
          Core::LinAlg::Matrix<nsd_, nen_>& velforce,          ///< rhs forces velocity
          Core::LinAlg::Matrix<nen_, 1>& preforce,             ///< rhs forces pressure
          Core::LinAlg::Matrix<nsd_ * nsd_, nen_>&
              lin_resM_Du,           ///< linearisation of the stabilization residual
          const double& fac3,        ///< factor for residual in current subgrid velocities
          const double& timefacfac,  ///< = timefac x fac
          const double& timefacfacpre,
          const double& rhsfac  ///< right-hand-side factor for residuals
      );

      //! Compute element matrix entries: reactive stabilization
      void reac_stab(Core::LinAlg::Matrix<nen_ * nsd_, nen_ * nsd_>&
                         estif_u,                              ///< block (weighting function v x u)
          Core::LinAlg::Matrix<nen_ * nsd_, nen_>& estif_p_v,  ///< block (weighting function v x p)
          Core::LinAlg::Matrix<nsd_, nen_>& velforce,          ///< rhs forces velocity
          Core::LinAlg::Matrix<nsd_ * nsd_, nen_>&
              lin_resM_Du,              ///< linearisation of the stabilization residual
          const double& timefacfac,     ///< = timefac x fac
          const double& timefacfacpre,  ///< = timefacpre x fac
          const double& rhsfac,         ///< right-hand-side factor for residuals
          const double& fac3            ///< factor for residual in current subgrid velocities
      );

      //! Compute element matrix entries: viscous stabilization
      void visc_stab(Core::LinAlg::Matrix<nen_ * nsd_, nen_ * nsd_>&
                         estif_u,                              ///< block (weighting function v x u)
          Core::LinAlg::Matrix<nen_ * nsd_, nen_>& estif_p_v,  ///< block (weighting function v x p)
          Core::LinAlg::Matrix<nsd_, nen_>& velforce,          ///< rhs forces velocity
          Core::LinAlg::Matrix<nsd_ * nsd_, nen_>&
              lin_resM_Du,              ///< linearisation of the stabilization residual
          const double& timefacfac,     ///< = timefac x fac
          const double& timefacfacpre,  ///< = timefac x fac
          const double& rhsfac,         ///< right-hand-side factor for residuals
          const double& fac3            ///< factor for residual in current subgrid velocities
      );

      //! Compute element matrix entries: convective divergence stabilization for XFEM
      void conv_div_stab(Core::LinAlg::Matrix<nen_ * nsd_, nen_ * nsd_>&
                             estif_u,                  ///< block (weighting function v x u)
          Core::LinAlg::Matrix<nsd_, nen_>& velforce,  ///< rhs forces velocity
          const double& timefacfac,                    ///< = timefac x fac
          const double& rhsfac                         ///< right-hand-side factor for residuals
      );

      //! Compute element matrix entries: cross stress stabilization
      void cross_stress_stab(Core::LinAlg::Matrix<nen_ * nsd_, nen_ * nsd_>&
                                 estif_u,                      ///< block (weighting function v x u)
          Core::LinAlg::Matrix<nen_ * nsd_, nen_>& estif_p_v,  ///< block (weighting function v x p)
          Core::LinAlg::Matrix<nsd_, nen_>& velforce,          ///< rhs forces velocity
          Core::LinAlg::Matrix<nsd_ * nsd_, nen_>&
              lin_resM_Du,           ///< linearisation of the stabilization residual
          const double& timefacfac,  ///< = timefac x fac
          const double& timefacfacpre,
          const double& rhsfac,  ///< right-hand-side factor for residuals
          const double& fac3     ///< factor for residual in current subgrid velocities
      );

      //! Compute element matrix entries: Reynolds stress stabilization
      void reynolds_stress_stab(Core::LinAlg::Matrix<nen_ * nsd_, nen_ * nsd_>&
                                    estif_u,                   ///< block (weighting function v x u)
          Core::LinAlg::Matrix<nen_ * nsd_, nen_>& estif_p_v,  ///< block (weighting function v x p)
          Core::LinAlg::Matrix<nsd_ * nsd_, nen_>&
              lin_resM_Du,           ///< linearisation of the stabilization residual
          const double& timefacfac,  ///< timefac x fac
          const double& timefacfacpre,
          const double& fac3  ///< factor for residual in current subgrid velocities
      );

      //! turbulence related methods
      //! definition in fluid_impl_turbulence_service.cpp

      //! get parameters for multifractal subgrid scales
      void prepare_multifractal_subgr_scales(
          Core::LinAlg::Matrix<nsd_, 1>&
              B_mfs,      ///< coefficient multifractal subgrid scales velocity
          double& D_mfs,  ///< coefficient multifractal subgrid scales scalar (loma only)
          const Core::LinAlg::Matrix<nsd_, nen_>& evelaf,  ///< velocity at time n+alpha_f / n+1
          const Core::LinAlg::Matrix<nsd_, nen_>&
              fsevelaf,     ///< fine scale velocity at time n+alpha_f / n+1
          const double vol  ///< volume
      );

      //! get turbulence parameter
      void get_turbulence_params(
          Teuchos::ParameterList& turbmodelparams,  ///< pointer general turbulence parameter list
          double& Cs_delta_sq,                      ///< parameter CS in dynamic Smagorinsky
          double& Ci_delta_sq,  ///< parameter CI in dynamic Smagorinsky for loma
          int& nlayer,  ///< number of layers for computation of parameter CS in dynamic Smagorinsky
          double CsDeltaSq,  ///< parameter CS in dynamic Smagorinsky computed in DynSmagFilter()
          double CiDeltaSq   ///< parameter CI in dynamic Smagorinsky for loma computed in
                             ///< DynSmagFilter()
      );

      //! calculate (all-scale) subgrid viscosity
      void calc_subgr_visc(
          const Core::LinAlg::Matrix<nsd_, nen_>& evelaf,  ///< velocity at time n+alpha_f / n+1
          const double vol,                                ///< volume
          double&
              Cs_delta_sq,     ///< parameter CS in dynamic Smagorinsky // !or Cv in dynamic Vreman!
          double& Ci_delta_sq  ///< parameter CS in dynamic Smagorinsky for loma
      );

      //! calculate fine-scale subgrid viscosity
      void calc_fine_scale_subgr_visc(
          const Core::LinAlg::Matrix<nsd_, nen_>& evelaf,  ///< velocity at time n+alpha_f / n+1
          const Core::LinAlg::Matrix<nsd_, nen_>&
              fsevelaf,     ///< fine scale velocity at time n+alpha_f / n+1
          const double vol  ///< volume
      );

      //! get coefficient for multifractal subgrid scales (velocity)
      void calc_multi_frac_subgrid_vel_coef(const double Csgs,  ///< model coefficient
          const double alpha,                                   ///< filter width ratio
          const std::vector<double> Nvel,                       ///< number of casacade steps
          Core::LinAlg::Matrix<nsd_, 1>& B_mfs                  ///< final coefficient
      );

      //! get coefficient for multifractal subgrid scales (scalar) (loma only)
      void calc_multi_frac_subgrid_sca_coef(const double Csgs,  ///< model coefficient
          const double alpha,                                   ///< model coefficient
          const double Pr,                                      ///< Prandtl number
          const double Pr_limit,  ///< Prandtl number to distinguish between low and high Prandtl
                                  ///< number regime
          const std::vector<double> Nvel,  ///< number of casacade steps (velocity)
          double Nphi,                     ///< number of casacade steps (scalar)
          double& D_mfs                    ///< final coefficient
      );

      //! Compute element matrix entries: fine scale subgrid viscousity rhs term
      void fine_scale_sub_grid_viscosity_term(
          Core::LinAlg::Matrix<nsd_, nen_>& velforce,  ///< rhs forces velocity
          const double& fssgviscfac  ///< = (fine scale subgrid viscousity) x timefacfac
      );

      void multfrac_sub_grid_scales_cross(Core::LinAlg::Matrix<nen_ * nsd_, nen_ * nsd_>& estif_u,
          Core::LinAlg::Matrix<nsd_, nen_>& velforce, const double& timefacfac,
          const double& rhsfac);

      void multfrac_sub_grid_scales_reynolds(
          Core::LinAlg::Matrix<nen_ * nsd_, nen_ * nsd_>& estif_u,
          Core::LinAlg::Matrix<nsd_, nen_>& velforce, const double& timefacfac,
          const double& rhsfac);

      void multfrac_sub_grid_scales_consistent_residual();

      //! loma related methods
      //! definition in fluid_impl_loma_service.cpp

      //! update material parameters including subgrid-scale part of scalar
      void update_material_params(
          Teuchos::RCP<const Core::Mat::Material> material,  ///< reference pointer to material
          const Core::LinAlg::Matrix<nsd_, nen_>& evelaf,    ///< velocity at time n+alpha_f / n+1
          const Core::LinAlg::Matrix<nen_, 1>& epreaf,       ///< pressure at time n+alpha_f / n+1
          const Core::LinAlg::Matrix<nen_, 1>& epream,       ///< pressure at time n+alpha_m / n
          const Core::LinAlg::Matrix<nen_, 1>& escaaf,       ///< scalar at time n+alpha_f / n+1
          const Core::LinAlg::Matrix<nen_, 1>& escaam,       ///< scalar at time n+alpha_m / n
          const double thermpressaf,  ///< thermodynamic pressure at time n+alpha_f / n+1
          const double thermpressam,  ///< thermodynamic pressure at time n+alpha_m / n
          const double sgsca          ///< subgrid scalar at integration point
      );

      //! compute additional Galerkin terms on right-hand side of continuity equation
      //! (only required for variable-density flow at low Mach number)
      void compute_gal_rhs_cont_eq(
          const Core::LinAlg::Matrix<nsd_, nen_>& eveln,  ///< velocity at time n
          const Core::LinAlg::Matrix<nen_, 1>& escaaf,    ///< scalar at time n+alpha_F/n+1
          const Core::LinAlg::Matrix<nen_, 1>& escaam,    ///< scalar at time n+alpha_M/n
          const Core::LinAlg::Matrix<nen_, 1>& escadtam,  ///< acceleration at time n+alpha_M/n
          bool isale                                      ///< flag for ALE case
      );

      //! compute additional Galerkin terms on right-hand side of continuity equation
      //! (only required for weakly compressibility)
      void compute_gal_rhs_cont_eq_weak_comp(
          const Core::LinAlg::Matrix<nen_, 1>& epreaf,  ///< pressure at time n+alpha_F/n+1
          const Core::LinAlg::Matrix<nen_, 1>&
              epredtam,  ///< derivative of pressure at time n+alpha_M/n
          bool isale     ///< flag for ALE case
      );

      //! compute additional Galerkin terms on right-hand side of continuity equation
      //! (only required for artificial compressibility)
      void compute_gal_rhs_cont_eq_art_comp(
          const Core::LinAlg::Matrix<nen_, 1>& epreaf,   ///< pressure at time n+alpha_F/n+1
          const Core::LinAlg::Matrix<nen_, 1>& epren,    ///< pressure at time n
          const Core::LinAlg::Matrix<nen_, 1>& escadtam  ///< acceleration at time n+alpha_M/n
      );

      //! compute residual of scalar equation and subgrid-scale part of scalar
      //! (only required for variable-density flow at low Mach number)
      void compute_subgrid_scale_scalar(
          const Core::LinAlg::Matrix<nen_, 1>& escaaf,  ///< scalar at time n+alpha_F/n+1
          const Core::LinAlg::Matrix<nen_, 1>& escaam   ///< scalar at time n+alpha_M/n
      );

      //! recompute Galerkin terms based on updated material parameters
      //! including s.-s. part of scalar and compute cross-stress term on
      //! right-hand side of continuity equation
      //! (only required for variable-density flow at low Mach number)
      void recompute_gal_and_compute_cross_rhs_cont_eq();

      //! Compute element matrix entries: LOMA
      void loma_gal_part(
          Core::LinAlg::Matrix<nen_, nen_ * nsd_>& estif_q_u,  ///< block (weighting function q x u)
          Core::LinAlg::Matrix<nen_, 1>& preforce,             ///< rhs forces pressure
          const double& timefacfac,                            ///< = timefac x fac
          const double& rhsfac  ///< right-hand-side factor for residuals
      );

      //! Compute element matrix entries: artificial compressibility
      void art_comp_pressure_inertia_gal_partand_cont_stab(
          Core::LinAlg::Matrix<nen_ * nsd_, nen_>& estif_p_v,  ///< block (weighting function v x p)
          Core::LinAlg::Matrix<nen_, nen_>& ppmat              ///< block (weighting function q x p)
      );

      //! Compute element matrix entries: weak compressibility
      void weak_comp_pressure_inertia_gal_part(
          Core::LinAlg::Matrix<nen_ * nsd_, nen_>& estif_p_v,  ///< block (weighting function v x p)
          Core::LinAlg::Matrix<nen_, nen_>& ppmat              ///< block (weighting function q x p)
      );

      //! ale related methods
      //! definition in fluid_impl_ale_service.cpp

      //! linearisation in the case of mesh motion 2-D
      virtual void lin_mesh_motion_2_d(
          Core::LinAlg::Matrix<(nsd_ + 1) * nen_, (nsd_ + 1) * nen_>& emesh,  ///< mesh motion
          const Core::LinAlg::Matrix<nsd_, nen_>& evelaf,  ///< velocity at time n+alpha_f / n+1
          const double& press,                             ///< pressure at integration point
          const double& timefac,                           ///< time factor
          const double& timefacfac                         ///< = timefac x fac
      );

      //! linearisation in the case of mesh motion 3-D
      virtual void lin_mesh_motion_3_d(
          Core::LinAlg::Matrix<(nsd_ + 1) * nen_, (nsd_ + 1) * nen_>& emesh,  ///< mesh motion
          const Core::LinAlg::Matrix<nsd_, nen_>& evelaf,  ///< velocity at time n+alpha_f / n+1
          const double& press,                             ///< pressure at integration point
          const double& timefac,                           ///< time factor
          const double& timefacfac                         ///< = timefac x fac
      );

      void get_porosity_at_gp(const Core::LinAlg::Matrix<nen_, 1>& eporo);

      /*!
        \brief calculate rate of strain of (fine-scale) velocity

      @param evel       (i) nodal velocity values
      @param derxy      (i) shape function derivatives
      @param velderxy   (o) velocity derivatives


        \return computed rate of strain
      */
      double get_strain_rate(const Core::LinAlg::Matrix<nsd_, nen_>& evel)
      {
        double rateofstrain = 0.0;

        // velderxy is computed here since the evaluation of the strain rate can be performed
        // at the element center before the gauss loop

        // get velocity derivatives at integration point
        //
        //              +-----  dN (x)
        //   dvel (x)    \        k
        //   -------- =   +     ------ * vel
        //      dx       /        dx        k
        //        j     +-----      j
        //              node k
        //
        // j : direction of derivative x/y/z
        //
        Core::LinAlg::Matrix<nsd_, nsd_> velderxy;
        velderxy.multiply_nt(evel, derxy_);

        // compute (resolved) rate of strain
        //
        //          +-                                 -+ 1
        //          |          /   \           /   \    | -
        //          | 2 * eps | vel |   * eps | vel |   | 2
        //          |          \   / ij        \   / ij |
        //          +-                                 -+
        //
        Core::LinAlg::Matrix<nsd_, nsd_> two_epsilon;
        for (int rr = 0; rr < nsd_; ++rr)
        {
          for (int mm = 0; mm < nsd_; ++mm)
          {
            two_epsilon(rr, mm) = velderxy(rr, mm) + velderxy(mm, rr);
          }
        }

        for (int rr = 0; rr < nsd_; ++rr)
        {
          for (int mm = 0; mm < nsd_; ++mm)
          {
            rateofstrain += two_epsilon(rr, mm) * two_epsilon(mm, rr);
          }
        }

        // sqrt(two_epsilon(rr,mm)*two_epsilon(mm,rr)/4.0*2.0)

        return (sqrt(rateofstrain / 2.0));
      }

      //! output values of Cs, visceff and Cs_delta_sq for statistics
      void store_model_parameters_for_output(const double Cs_delta_sq, const double Ci_delta_sq,
          const int nlayer, const bool isowned, Teuchos::ParameterList& turbmodelparams);

      /*!
       * \brief fill elment matrix and vectors with the global values
       */
      void extract_values_from_global_vector(
          const Core::FE::Discretization& discretization,  ///< discretization
          const std::vector<int>& lm,                      ///<
          FLD::RotationallySymmetricPeriodicBC<distype, nsd_ + 1, enrtype>& rotsymmpbc,  ///<
          Core::LinAlg::Matrix<nsd_, nen_>* matrixtofill,  ///< vector field
          Core::LinAlg::Matrix<nen_, 1>* vectortofill,     ///< scalar field
          const std::string state);                        ///< state of the global vector

      //! identify elements of inflow section
      void inflow_element(Core::Elements::Element* ele);

      // FLD::RotationallySymmetricPeriodicBC<distype> & rotsymmpbc, ///<
      //  {
      //    // get state of the global vector
      //    Teuchos::RCP<const Epetra_Vector> matrix_state = discretization.GetState(state);
      //    if(matrix_state == Teuchos::null)
      //      FOUR_C_THROW("Cannot get state vector %s", state.c_str());
      //
      //    // extract local values of the global vectors
      //    std::vector<double> mymatrix(lm.size());
      //    Core::FE::ExtractMyValues(*matrix_state,mymatrix,lm);
      //
      //    // rotate the vector field in the case of rotationally symmetric boundary conditions
      //    if(matrixtofill != nullptr)
      //      rotsymmpbc.rotate_my_values_if_necessary(mymatrix);
      //
      //    for (int inode=0; inode<nen_; ++inode)  // number of nodes
      //    {
      //      // fill a vector field via a pointer
      //      if (matrixtofill != nullptr)
      //      {
      //        for(int idim=0; idim<nsd_; ++idim) // number of dimensions
      //        {
      //          (*matrixtofill)(idim,inode) = mymatrix[idim+(inode*numdofpernode_)];
      //        }  // end for(idim)
      //      }
      //      // fill a scalar field via a pointer
      //      if (vectortofill != nullptr)
      //        (*vectortofill)(inode,0) = mymatrix[nsd_+(inode*numdofpernode_)];
      //    }
      //  }



      //==================================================================================
      // OLD FLUID ELE CALC ROUTINES BEFORE OST-HIST MIGRATION.

      //! calculate div(epsilon(u))
      void calc_div_eps(
          const Core::LinAlg::Matrix<nsd_, nen_>& evelaf);  ///< velocity at time n+alpha_f / n+1

      //! Compute element matrix and rhs entries: inertia, convective andyn
      //! reactive terms of the Galerkin part
      virtual void inertia_convection_reaction_gal_part(
          Core::LinAlg::Matrix<nen_ * nsd_, nen_ * nsd_>&
              estif_u,                                 ///< block (weighting function v x u)
          Core::LinAlg::Matrix<nsd_, nen_>& velforce,  ///< rhs forces velocity
          Core::LinAlg::Matrix<nsd_ * nsd_, nen_>&
              lin_resM_Du,  ///< linearisation of the Garlerkin momentum residual
          Core::LinAlg::Matrix<nsd_, 1>&
              resM_Du,          ///< linearisation of the Garlerkin momentum residual
          const double& rhsfac  ///< right-hand-side factor
      );

      //! Compute element matrix entries: for the viscous terms of the Galerkin part
      void viscous_gal_part(Core::LinAlg::Matrix<nen_ * nsd_, nen_ * nsd_>&
                                estif_u,                 ///< block (weighting function v x u)
          Core::LinAlg::Matrix<nsd_, nen_>& velforce,    ///< rhs forces velocity
          Core::LinAlg::Matrix<nsd_, nsd_>& viscstress,  ///< viscous stresses
          const double& timefacfac,                      ///< = timefac x fac
          const double& rhsfac                           ///< right-hand-side factor
      );

      //! Compute element matrix entries: div-grad stabilization and the rhs of the viscous term
      void cont_stab(Core::LinAlg::Matrix<nen_ * nsd_, nen_ * nsd_>&
                         estif_u,                      ///< block (weighting function v x u)
          Core::LinAlg::Matrix<nsd_, nen_>& velforce,  ///< rhs forces velocity
          const double& timefac,                       ///< time factor
          const double& timefacfac,                    ///< = timefac x fac
          const double& timefacfacpre,                 ///< = timefac x fac
          const double& rhsfac                         ///< right-hand-side factor
      );

      //! Compute element matrix entries: pressure terms of the Garlerkin part and rhs
      void pressure_gal_part(Core::LinAlg::Matrix<nen_ * nsd_, nen_>&
                                 estif_p_v,            ///< block (weighting function v x pressure)
          Core::LinAlg::Matrix<nsd_, nen_>& velforce,  ///< rhs forces velocity
          const double& timefacfac,                    ///< = timefac x fac
          const double& timefacfacpre,
          const double& rhsfac,  ///< right-hand-side factor
          const double& press    ///< pressure at integration point
      );

      //! Compute element matrix entries: continuity terms of the Garlerkin part and rhs
      virtual void continuity_gal_part(
          Core::LinAlg::Matrix<nen_, nen_ * nsd_>& estif_q_u,  ///< block (weighting function q x u)
          Core::LinAlg::Matrix<nen_, 1>& preforce,             ///< rhs forces pressure
          const double& timefacfac,                            ///< = timefac x fac
          const double& timefacfacpre,
          const double& rhsfac  ///< right-hand-side factor
      );

      //! Compute element matrix entries: body force terms on rhs
      void body_force_rhs_term(Core::LinAlg::Matrix<nsd_, nen_>& velforce,  ///< rhs forces velocity
          const double& rhsfac);

      //==================================================================================



      //! for the handling of rotationally symmetric periodic boundary conditions
      Teuchos::RCP<FLD::RotationallySymmetricPeriodicBC<distype, nsd_ + 1, enrtype>> rotsymmpbc_;
      //! element id
      int eid_;
      //! Flag to (de)activate higher order elements
      //! elements with only mixed second order derivatives are not counted as higher order elements
      //! (see definition of higher order elements in fluid3_ele_impl_utils.cpp)
      bool is_higher_order_ele_;
      //! pointer to parameter lists
      Discret::ELEMENTS::FluidEleParameter* fldpara_;
      //! pointer to parameter list for time integration
      Discret::ELEMENTS::FluidEleParameterTimInt* fldparatimint_;
      //! element type: nurbs
      bool isNurbs_;
      //! weights for nurbs elements
      Core::LinAlg::Matrix<nen_, 1> weights_;
      //! knot vector for nurbs elements
      std::vector<Core::LinAlg::SerialDenseVector> myknots_;
      //! Gaussian integration points
      Core::FE::GaussIntegration intpoints_;
      //! identify elements of inflow section
      //! required for turbulence modeling
      bool is_inflow_ele_;

      //========================================================

      //! element stiffness block velocity-velocity
      Core::LinAlg::Matrix<nen_ * nsd_, nen_ * nsd_> estif_u_;
      //! element stiffness block pressure-velocity
      Core::LinAlg::Matrix<nen_ * nsd_, nen_> estif_p_v_;
      //! element stiffness block velocity-pressure
      Core::LinAlg::Matrix<nen_, nen_ * nsd_> estif_q_u_;
      //! element stiffness block pressure-pressure
      Core::LinAlg::Matrix<nen_, nen_> ppmat_;

      // definition of vectors
      //! element rhs blocks pressure
      Core::LinAlg::Matrix<nen_, 1> preforce_;
      //! element rhs blocks velocity
      Core::LinAlg::Matrix<nsd_, nen_> velforce_;

      //! definition of velocity-based momentum residual vectors
      Core::LinAlg::Matrix<nsd_ * nsd_, nen_> lin_resM_Du_;
      Core::LinAlg::Matrix<nsd_, 1> resM_Du_;


      //========================================================

      //! nodal based quantities for an element
      Core::LinAlg::Matrix<nsd_, nen_> ebofoaf_;
      Core::LinAlg::Matrix<nsd_, nen_> eprescpgaf_;
      Core::LinAlg::Matrix<nen_, 1> escabofoaf_;

      Core::LinAlg::Matrix<nsd_, nen_> ebofon_;
      Core::LinAlg::Matrix<nsd_, nen_> eprescpgn_;
      Core::LinAlg::Matrix<nen_, 1> escabofon_;

      Core::LinAlg::Matrix<nsd_, nen_> evelaf_;
      Core::LinAlg::Matrix<nen_, 1> epreaf_;
      Core::LinAlg::Matrix<nsd_, nen_> evelam_;
      Core::LinAlg::Matrix<nen_, 1> epream_;
      Core::LinAlg::Matrix<nsd_, nen_> evelnp_;
      Core::LinAlg::Matrix<nen_, 1> eprenp_;
      Core::LinAlg::Matrix<nsd_, nen_> eveln_;
      Core::LinAlg::Matrix<nen_, 1> epren_;
      Core::LinAlg::Matrix<nsd_, nen_> eaccam_;
      Core::LinAlg::Matrix<nen_, 1> escadtam_;
      Core::LinAlg::Matrix<nsd_, nen_> eveldtam_;
      Core::LinAlg::Matrix<nen_, 1> epredtam_;
      Core::LinAlg::Matrix<nen_, 1> escaaf_;
      Core::LinAlg::Matrix<nen_, 1> escaam_;
      Core::LinAlg::Matrix<nsd_, nen_> emhist_;
      Core::LinAlg::Matrix<nen_, 1> eporo_;

      Core::LinAlg::Matrix<nsd_, nen_> gradphiele_;
      Core::LinAlg::Matrix<nen_, 1> curvatureele_;
      Core::LinAlg::Matrix<nsd_, nen_> gradphielen_;
      Core::LinAlg::Matrix<nen_, 1> curvatureelen_;

      Core::LinAlg::Matrix<nsd_, 2 * nen_> gradphieletot_;
      Core::LinAlg::Matrix<nen_, 2> curvatureeletot_;

      Core::LinAlg::Matrix<nsd_, nen_> edispnp_;
      Core::LinAlg::Matrix<nsd_, nen_> egridv_;
      Core::LinAlg::Matrix<nsd_, nen_> egridvn_;

      Core::LinAlg::Matrix<nsd_, nen_> fsevelaf_;
      Core::LinAlg::Matrix<nen_, 1> fsescaaf_;

      Core::LinAlg::Matrix<nsd_, nen_> evel_hat_;
      Core::LinAlg::Matrix<nsd_ * nsd_, nen_> ereynoldsstress_hat_;


      //! node coordinates
      Core::LinAlg::Matrix<nsd_, nen_> xyze_;
      //! array for shape functions
      Core::LinAlg::Matrix<nen_, 1> funct_;
      //! array for shape function derivatives w.r.t r,s,t
      Core::LinAlg::Matrix<nsd_, nen_> deriv_;
      //! array for second derivatives of shape function w.r.t r,s,t
      Core::LinAlg::Matrix<numderiv2_, nen_> deriv2_;
      //! transposed jacobian "dx/ds"
      Core::LinAlg::Matrix<nsd_, nsd_> xjm_;
      //! inverse of transposed jacobian "ds/dx"
      Core::LinAlg::Matrix<nsd_, nsd_> xji_;
      //! global velocity derivatives in gausspoint w.r.t x,y,z
      Core::LinAlg::Matrix<nsd_, nsd_> vderxy_;
      //! global derivatives of shape functions w.r.t x,y,z
      Core::LinAlg::Matrix<nsd_, nen_> derxy_;
      //! global second derivatives of shape functions w.r.t x,y,z
      Core::LinAlg::Matrix<numderiv2_, nen_> derxy2_;
      //! bodyforce in gausspoint
      Core::LinAlg::Matrix<nsd_, 1> bodyforce_;
      // New One Step Theta variables
      //======================================================
      //! denisty multiplied with instationary term for OST implementations
      double dens_theta_;
      //! bodyforce in gausspoint (n)
      Core::LinAlg::Matrix<nsd_, 1> bodyforcen_;
      //! (u^n_old*nabla)u^n_old
      Core::LinAlg::Matrix<nsd_, 1> conv_oldn_;
      //! div epsilon(u^n_old)
      Core::LinAlg::Matrix<nsd_, 1> visc_oldn_;
      //! pressure gradient in gausspoint (n)
      Core::LinAlg::Matrix<nsd_, 1> gradpn_;
      //! velocity vector in gausspoint (n)
      Core::LinAlg::Matrix<nsd_, 1> velintn_;
      //! physical viscosity (n)
      double viscn_;
      //! old residual of continuity equation (n)
      double conres_oldn_;
      //! prescribed pressure gradient (required for turbulent channel flow!)
      Core::LinAlg::Matrix<nsd_, 1> generalbodyforcen_;
      //======================================================
      //! prescribed pressure gradient (required for turbulent channel flow!)
      Core::LinAlg::Matrix<nsd_, 1> generalbodyforce_;
      //! vector containing all values from previous timelevel n for momentum equation
      Core::LinAlg::Matrix<nsd_, 1> histmom_;
      //! velocity vector in gausspoint
      Core::LinAlg::Matrix<nsd_, 1> velint_;
      //! subgrid-scale velocity vector in gausspoint
      Core::LinAlg::Matrix<nsd_, 1> sgvelint_;
      //! grid velocity u_G at integration point
      Core::LinAlg::Matrix<nsd_, 1> gridvelint_;
      //! grid velocity u_G at integration point for new ost
      Core::LinAlg::Matrix<nsd_, 1> gridvelintn_;
      //! ale convective velocity c=u-u_G at integration point
      Core::LinAlg::Matrix<nsd_, 1> convvelint_;
      //! Oseen advective velocity at element nodes
      Core::LinAlg::Matrix<nsd_, nen_> eadvvel_;
      //! acceleration vector in gausspoint
      Core::LinAlg::Matrix<nsd_, 1> accint_;
      //! pressure gradient in gausspoint
      Core::LinAlg::Matrix<nsd_, 1> gradp_;
      //! the stabilisation parameters -> it is a (3,1) vector for 2D and 3D
      Core::LinAlg::Matrix<3, 1> tau_;
      //! viscous term including 2nd derivatives
      //! (This array once had three dimensions, now the first two are combined to one.)
      Core::LinAlg::Matrix<nsd_ * nsd_, nen_> viscs2_;
      //! linearisation of convection, convective part
      Core::LinAlg::Matrix<nen_, 1> conv_c_;
      //! linearisation of subgrid-scale convection, convective part
      Core::LinAlg::Matrix<nen_, 1> sgconv_c_;
      //! velocity divergenceat at t_(n+alpha_F) or t_(n+1)
      double vdiv_;
      //! total right hand side terms at int.-point for momentum equation
      Core::LinAlg::Matrix<nsd_, 1> rhsmom_;
      //! (u_old*nabla)u_old
      Core::LinAlg::Matrix<nsd_, 1> conv_old_;
      //! div epsilon(u_old)
      Core::LinAlg::Matrix<nsd_, 1> visc_old_;
      //! old residual of momentum equation
      Core::LinAlg::Matrix<nsd_, 1> momres_old_;
      //! old residual of continuity equation
      double conres_old_;
      //! 2nd derivatives of coord.-functions w.r.t r,s,t
      Core::LinAlg::Matrix<numderiv2_, nsd_> xder2_;
      //! global velocity second derivatives in gausspoint w.r.t local coordinates
      Core::LinAlg::Matrix<nsd_, nsd_> vderiv_;
      //! coordinates of current integration point in reference coordinates
      Core::LinAlg::Matrix<nsd_, 1> xsi_;
      //! Jacobian determinant
      double det_;
      //! integration factor
      double fac_;
      //! physical viscosity
      double visc_;
      //! effective viscosity = physical viscosity + (all-scale) subgrid viscosity
      double visceff_;
      //! reaction coefficient
      double reacoeff_;

      //! Two-phase specific variables:
      //! surface tension
      double gamma_;

      //! LOMA-specific variables:
      //! physical diffusivity of scalar equation
      double diffus_;
      //! right-hand-side term at int.-point for continuity equation
      double rhscon_;
      //! density at t_(n+alpha_F) or t_(n+1)
      double densaf_;
      //! density at t_(n+alpha_M)
      double densam_;
      //! density at t_(n)
      double densn_;
      //! delta density for Boussinesq Approximation
      double deltadens_;
      //! factor for scalar time derivative
      double scadtfac_;
      //! factor for convective scalar term at t_(n+alpha_F) or t_(n+1)
      double scaconvfacaf_;
      //! factor for convective scalar term at t_(n)
      double scaconvfacn_;
      //! addition to continuity equation due to thermodynamic pressure
      double thermpressadd_;
      //! convective velocity vector in gausspoint at t_(n)
      Core::LinAlg::Matrix<nsd_, 1> convvelintn_;
      //! global velocity derivatives in gausspoint w.r.t x,y,z at t_(n)
      Core::LinAlg::Matrix<nsd_, nsd_> vderxyn_;
      //! velocity divergence at at t_(n)
      double vdivn_;
      //! scalar gradient at t_(n+alpha_F) or t_(n+1)
      Core::LinAlg::Matrix<nsd_, 1> grad_scaaf_;
      //! scalar gradient at t_(n)
      Core::LinAlg::Matrix<nsd_, 1> grad_scan_;
      //! scalar at t_(n+alpha_F) or t_(n+1)
      double scaaf_;
      //! scalar at t_(n)
      double scan_;
      //! time derivative of scalar term (only required for generalized-alpha scheme)
      double tder_sca_;
      //! convective scalar term at t_(n+alpha_F) or t_(n+1)
      double conv_scaaf_;
      //! convective scalar term at t_(n)
      double conv_scan_;
      //! right-hand side of scalar equation
      double scarhs_;
      //! subgrid-scale part of scalar at integration point
      double sgscaint_;

      //! weakly_compressible-specific variables:
      //! pressure at t_(n+alpha_F) or t_(n+1)
      double preaf_;
      //! pressure at t_(n+alpha_M) or t_(n)
      double pream_;
      //! factor for convective pressure term at t_(n+alpha_F) or t_(n+1)
      double preconvfacaf_;
      //! time derivative of pressure
      double tder_pre_;
      //! factor for pressure time derivative
      double predtfac_;
      //! pressure gradient at t_(n+alpha_F) or t_(n+1)
      Core::LinAlg::Matrix<nsd_, 1> grad_preaf_;
      //! convective pressure term at t_(n+alpha_F) or t_(n+1)
      double conv_preaf_;
      //! // element correction term
      Core::LinAlg::Matrix<1, nen_> ecorrectionterm_;

      //! turbulence-specific variables:
      //! fine-scale velocity vector in gausspoint
      Core::LinAlg::Matrix<nsd_, 1> fsvelint_;
      //! fine-scale velocity vector in gausspoint for multifractal subgrid-scale modeling
      Core::LinAlg::Matrix<nsd_, 1> mffsvelint_;
      //! fine-scale global velocity derivatives in gausspoint w.r.t x,y,z
      Core::LinAlg::Matrix<nsd_, nsd_> fsvderxy_;
      //! fine-scale global velocity derivatives in gausspoint w.r.t x,y,z for multifractal
      //! subgrid-scale modeling
      Core::LinAlg::Matrix<nsd_, nsd_> mffsvderxy_;
      //! fine scale velocity divergence for multifractal subgrid-scale modeling
      double mffsvdiv_;
      //! (all-scale) subgrid viscosity
      double sgvisc_;
      //! fine-scale subgrid viscosity
      double fssgvisc_;
      //! model parameter for isotropic part of subgrid-stress tensor (dyn Smag for loma)
      double q_sq_;
      //! multifractal subgrid-scale part of scalar at integration point
      double mfssgscaint_;
      //! gradient of multifractal subgrid-scale scalar (for loma)
      Core::LinAlg::Matrix<nsd_, 1> grad_fsscaaf_;

      //! norm of velocity at integration point at time t^{n+1}
      double vel_normnp_;
      //! time-dependent subgrid-scales (pointer to element-specific data)
      Teuchos::RCP<FLD::TDSEleData> tds_;

      // polynomial pressure projection matrices
      double D_;
      Core::LinAlg::Matrix<nen_, 1> E_;

      Core::LinAlg::Matrix<nsd_ * nsd_, nen_> evelafgrad_;
      Core::LinAlg::Matrix<nsd_ * nsd_, nen_> evelngrad_;


      // protected:
      // static std::map<int,std::map<int,Discret::ELEMENTS::FluidEleCalc<distype>* >* > instances_;
    };
    // template<Core::FE::CellType distype, Discret::ELEMENTS::Fluid::EnrichmentType
    // enrtype = Discret::ELEMENTS::Fluid::none>
    // std::map<int,std::map<int,Discret::ELEMENTS::FluidEleCalc<distype>* >* >
    // FluidEleCalc<distype>::instances_;
  }  // namespace ELEMENTS
}  // namespace Discret

FOUR_C_NAMESPACE_CLOSE

#endif
