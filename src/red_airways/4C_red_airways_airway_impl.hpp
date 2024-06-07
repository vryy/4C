/*---------------------------------------------------------------------*/
/*! \file

\brief Internal implementation of airway_impl element


\level 3

*/
/*---------------------------------------------------------------------*/



#ifndef FOUR_C_RED_AIRWAYS_AIRWAY_IMPL_HPP
#define FOUR_C_RED_AIRWAYS_AIRWAY_IMPL_HPP

#include "4C_config.hpp"

#include "4C_fem_general_utils_local_connectivity_matrices.hpp"
#include "4C_linalg_fixedsizematrix.hpp"
#include "4C_red_airways_elementbase.hpp"

FOUR_C_NAMESPACE_OPEN


namespace Discret
{
  namespace ELEMENTS
  {
    /// Interface base class for airway_impl
    /*!
      This class exists to provide a common interface for all template
      versions of airway_impl. The only function
      this class actually defines is Impl, which returns a pointer to
      the appropriate version of airway_impl.
     */
    class RedAirwayImplInterface
    {
     public:
      /// Empty constructor
      RedAirwayImplInterface() {}
      /// Empty destructor
      virtual ~RedAirwayImplInterface() = default;  /// Evaluate the element
      virtual int Evaluate(RedAirway* ele, Teuchos::ParameterList& params,
          Discret::Discretization& discretization, std::vector<int>& lm,
          Core::LinAlg::SerialDenseMatrix& elemat1_epetra,
          Core::LinAlg::SerialDenseMatrix& elemat2_epetra,
          Core::LinAlg::SerialDenseVector& elevec1_epetra,
          Core::LinAlg::SerialDenseVector& elevec2_epetra,
          Core::LinAlg::SerialDenseVector& elevec3_epetra,
          Teuchos::RCP<Core::Mat::Material> mat) = 0;

      virtual void Initial(RedAirway* ele, Teuchos::ParameterList& params,
          Discret::Discretization& discretization, std::vector<int>& lm,
          Core::LinAlg::SerialDenseVector& radii_in, Core::LinAlg::SerialDenseVector& radii_out,
          Teuchos::RCP<const Core::Mat::Material> material) = 0;

      virtual void EvaluateTerminalBC(RedAirway* ele, Teuchos::ParameterList& params,
          Discret::Discretization& discretization, std::vector<int>& lm,
          Core::LinAlg::SerialDenseVector& elevec1_epetra,
          Teuchos::RCP<Core::Mat::Material> mat) = 0;

      virtual void CalcFlowRates(RedAirway* ele, Teuchos::ParameterList& params,
          Discret::Discretization& discretization, std::vector<int>& lm,
          Teuchos::RCP<Core::Mat::Material> mat) = 0;

      virtual void CalcElemVolume(RedAirway* ele, Teuchos::ParameterList& params,
          Discret::Discretization& discretization, std::vector<int>& lm,
          Teuchos::RCP<Core::Mat::Material> mat) = 0;

      virtual void GetCoupledValues(RedAirway* ele, Teuchos::ParameterList& params,
          Discret::Discretization& discretization, std::vector<int>& lm,
          Teuchos::RCP<Core::Mat::Material> material) = 0;

      virtual void get_junction_volume_mix(RedAirway* ele, Teuchos::ParameterList& params,
          Discret::Discretization& discretization, Core::LinAlg::SerialDenseVector& volumeMix_np,
          std::vector<int>& lm, Teuchos::RCP<Core::Mat::Material> material) = 0;

      virtual void CalcCFL(RedAirway* ele, Teuchos::ParameterList& params,
          Discret::Discretization& discretization, std::vector<int>& lm,
          Teuchos::RCP<Core::Mat::Material> material) = 0;

      virtual void update_scatra(RedAirway* ele, Teuchos::ParameterList& params,
          Discret::Discretization& discretization, std::vector<int>& lm,
          Teuchos::RCP<Core::Mat::Material> material) = 0;

      virtual void UpdateElem12Scatra(RedAirway* ele, Teuchos::ParameterList& params,
          Discret::Discretization& discretization, std::vector<int>& lm,
          Teuchos::RCP<Core::Mat::Material> material) = 0;

      virtual void eval_nodal_essential_values(RedAirway* ele, Teuchos::ParameterList& params,
          Discret::Discretization& discretization, Core::LinAlg::SerialDenseVector& nodal_surface,
          Core::LinAlg::SerialDenseVector& nodal_volume,
          Core::LinAlg::SerialDenseVector& nodal_flow, std::vector<int>& lm,
          Teuchos::RCP<Core::Mat::Material> material) = 0;


      /// Internal implementation class for airway element
      static RedAirwayImplInterface* Impl(Discret::ELEMENTS::RedAirway* airway);
    };


    /// Internal airway implementation
    /*!
      This internal class keeps all the working arrays needed to
      calculate the airway element. Additionally the method Sysmat()
      provides a clean and fast element implementation.

      <h3>Purpose</h3>

      \author ismail
      \date 01/09
    */

    template <Core::FE::CellType distype>
    class AirwayImpl : public RedAirwayImplInterface
    {
     public:
      /// Default constructor. Will need to query data from ParameterList.
      AirwayImpl() = default;

      //! number of nodes
      static constexpr int iel = Core::FE::num_nodes<distype>;


      /// Evaluate
      /*!
        The evaluate function for the general airway case.
       */
      int Evaluate(RedAirway* ele, Teuchos::ParameterList& params,
          Discret::Discretization& discretization, std::vector<int>& lm,
          Core::LinAlg::SerialDenseMatrix& elemat1_epetra,
          Core::LinAlg::SerialDenseMatrix& elemat2_epetra,
          Core::LinAlg::SerialDenseVector& elevec1_epetra,
          Core::LinAlg::SerialDenseVector& elevec2_epetra,
          Core::LinAlg::SerialDenseVector& elevec3_epetra,
          Teuchos::RCP<Core::Mat::Material> mat) override;

      /*!
        \brief Calculate virtual trajectory xnp and state of airway (open/closed)

        \param ele              (i) the element those matrix is calculated
        \param epn              (i) nodal pressure at n
        \param dt               (i) timestep
        */
      void EvaluateCollapse(RedAirway* ele, Core::LinAlg::SerialDenseVector& epn,
          Teuchos::ParameterList& params, double dt);
      /*!
        \brief Calculate Pextn and Pextnp

        \param ele              (i) the element those matrix is calculated
        \param epn              (i) nodal pressure at n
        \param dt               (i) timestep
        */
      void ComputePext(RedAirway* ele, Teuchos::RCP<const Epetra_Vector> pn,
          Teuchos::RCP<const Epetra_Vector> pnp, Teuchos::ParameterList& params);


      void EvaluateTerminalBC(RedAirway* ele, Teuchos::ParameterList& params,
          Discret::Discretization& discretization, std::vector<int>& lm,
          Core::LinAlg::SerialDenseVector& rhs, Teuchos::RCP<Core::Mat::Material> mat) override;

      /*!
        \brief get the initial values of the degrees of freedome at the node

        \param ele              (i) the element those matrix is calculated
        \param eqnp             (i) nodal volumetric flow rate at n+1
        \param evelnp           (i) nodal velocity at n+1
        \param eareanp          (i) nodal cross-sectional area at n+1
        \param eprenp           (i) nodal pressure at n+1
        \param estif            (o) element matrix to calculate
        \param eforce           (o) element rhs to calculate
        \param material         (i) airway material/dimesion
        \param time             (i) current simulation time
        \param dt               (i) timestep
        */
      void Initial(RedAirway* ele, Teuchos::ParameterList& params,
          Discret::Discretization& discretization, std::vector<int>& lm,
          Core::LinAlg::SerialDenseVector& radii_in, Core::LinAlg::SerialDenseVector& radii_out,
          Teuchos::RCP<const Core::Mat::Material> material) override;

      /*!
       \Essential functions to compute the results of essentail matrices
      */
      void CalcFlowRates(RedAirway* ele, Teuchos::ParameterList& params,
          Discret::Discretization& discretization, std::vector<int>& lm,
          Teuchos::RCP<Core::Mat::Material> mat) override;

      /*!
       \Essential functions to compute the volume of elements
      */
      void CalcElemVolume(RedAirway* ele, Teuchos::ParameterList& params,
          Discret::Discretization& discretization, std::vector<int>& lm,
          Teuchos::RCP<Core::Mat::Material> mat) override;

      /*!
       \Essential functions to compute the volume mixing and  flowing into a junction
      */
      void get_junction_volume_mix(RedAirway* ele, Teuchos::ParameterList& params,
          Discret::Discretization& discretization,
          Core::LinAlg::SerialDenseVector& junctionVolumeMix_np, std::vector<int>& lm,
          Teuchos::RCP<Core::Mat::Material> material) override;

      /*!
       \Essential functions to evaluate the coupled results
      */
      void GetCoupledValues(RedAirway* ele, Teuchos::ParameterList& params,
          Discret::Discretization& discretization, std::vector<int>& lm,
          Teuchos::RCP<Core::Mat::Material> material) override;

      void CalcCFL(RedAirway* ele, Teuchos::ParameterList& params,
          Discret::Discretization& discretization, std::vector<int>& lm,
          Teuchos::RCP<Core::Mat::Material> material) override;

      void update_scatra(RedAirway* ele, Teuchos::ParameterList& params,
          Discret::Discretization& discretization, std::vector<int>& lm,
          Teuchos::RCP<Core::Mat::Material> material) override;

      void UpdateElem12Scatra(RedAirway* ele, Teuchos::ParameterList& params,
          Discret::Discretization& discretization, std::vector<int>& lm,
          Teuchos::RCP<Core::Mat::Material> material) override;

      void eval_nodal_essential_values(RedAirway* ele, Teuchos::ParameterList& params,
          Discret::Discretization& discretization, Core::LinAlg::SerialDenseVector& nodal_surface,
          Core::LinAlg::SerialDenseVector& nodal_volume,
          Core::LinAlg::SerialDenseVector& nodal_avg_scatra, std::vector<int>& lm,
          Teuchos::RCP<Core::Mat::Material> material) override;


     private:
    };

  }  // namespace ELEMENTS
}  // namespace Discret

FOUR_C_NAMESPACE_CLOSE

#endif
