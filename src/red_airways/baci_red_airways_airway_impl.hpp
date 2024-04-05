/*---------------------------------------------------------------------*/
/*! \file

\brief Internal implementation of airway_impl element


\level 3

*/
/*---------------------------------------------------------------------*/



#ifndef FOUR_C_RED_AIRWAYS_AIRWAY_IMPL_HPP
#define FOUR_C_RED_AIRWAYS_AIRWAY_IMPL_HPP

#include "baci_config.hpp"

#include "baci_discretization_fem_general_utils_local_connectivity_matrices.hpp"
#include "baci_linalg_fixedsizematrix.hpp"
#include "baci_red_airways_elementbase.hpp"

BACI_NAMESPACE_OPEN


namespace DRT
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
          DRT::Discretization& discretization, std::vector<int>& lm,
          CORE::LINALG::SerialDenseMatrix& elemat1_epetra,
          CORE::LINALG::SerialDenseMatrix& elemat2_epetra,
          CORE::LINALG::SerialDenseVector& elevec1_epetra,
          CORE::LINALG::SerialDenseVector& elevec2_epetra,
          CORE::LINALG::SerialDenseVector& elevec3_epetra, Teuchos::RCP<MAT::Material> mat) = 0;

      virtual void Initial(RedAirway* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, std::vector<int>& lm,
          CORE::LINALG::SerialDenseVector& radii_in, CORE::LINALG::SerialDenseVector& radii_out,
          Teuchos::RCP<const MAT::Material> material) = 0;

      virtual void EvaluateTerminalBC(RedAirway* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, std::vector<int>& lm,
          CORE::LINALG::SerialDenseVector& elevec1_epetra, Teuchos::RCP<MAT::Material> mat) = 0;

      virtual void CalcFlowRates(RedAirway* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, std::vector<int>& lm,
          Teuchos::RCP<MAT::Material> mat) = 0;

      virtual void CalcElemVolume(RedAirway* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, std::vector<int>& lm,
          Teuchos::RCP<MAT::Material> mat) = 0;

      virtual void GetCoupledValues(RedAirway* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, std::vector<int>& lm,
          Teuchos::RCP<MAT::Material> material) = 0;

      virtual void GetJunctionVolumeMix(RedAirway* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, CORE::LINALG::SerialDenseVector& volumeMix_np,
          std::vector<int>& lm, Teuchos::RCP<MAT::Material> material) = 0;

      virtual void SolveScatra(RedAirway* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, CORE::LINALG::SerialDenseVector& scatra_np,
          CORE::LINALG::SerialDenseVector& volumeMix_np, std::vector<int>& lm,
          Teuchos::RCP<MAT::Material> material) = 0;

      virtual void SolveScatraBifurcations(RedAirway* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, CORE::LINALG::SerialDenseVector& scatra_np,
          CORE::LINALG::SerialDenseVector& volumeMix_np, std::vector<int>& lm,
          Teuchos::RCP<MAT::Material> material) = 0;

      virtual void CalcCFL(RedAirway* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, std::vector<int>& lm,
          Teuchos::RCP<MAT::Material> material) = 0;

      virtual void UpdateScatra(RedAirway* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, std::vector<int>& lm,
          Teuchos::RCP<MAT::Material> material) = 0;

      virtual void UpdateElem12Scatra(RedAirway* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, std::vector<int>& lm,
          Teuchos::RCP<MAT::Material> material) = 0;

      virtual void EvalPO2FromScatra(RedAirway* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, std::vector<int>& lm,
          Teuchos::RCP<MAT::Material> material) = 0;

      virtual void EvalNodalEssentialValues(RedAirway* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, CORE::LINALG::SerialDenseVector& nodal_surface,
          CORE::LINALG::SerialDenseVector& nodal_volume,
          CORE::LINALG::SerialDenseVector& nodal_flow, std::vector<int>& lm,
          Teuchos::RCP<MAT::Material> material) = 0;


      /// Internal implementation class for airway element
      static RedAirwayImplInterface* Impl(DRT::ELEMENTS::RedAirway* airway);
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

    template <CORE::FE::CellType distype>
    class AirwayImpl : public RedAirwayImplInterface
    {
     public:
      /// Default constructor. Will need to query data from ParameterList.
      AirwayImpl() = default;

      //! number of nodes
      static constexpr int iel = CORE::FE::num_nodes<distype>;


      /// Evaluate
      /*!
        The evaluate function for the general airway case.
       */
      int Evaluate(RedAirway* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, std::vector<int>& lm,
          CORE::LINALG::SerialDenseMatrix& elemat1_epetra,
          CORE::LINALG::SerialDenseMatrix& elemat2_epetra,
          CORE::LINALG::SerialDenseVector& elevec1_epetra,
          CORE::LINALG::SerialDenseVector& elevec2_epetra,
          CORE::LINALG::SerialDenseVector& elevec3_epetra,
          Teuchos::RCP<MAT::Material> mat) override;

      /*!
        \brief Calculate virtual trajectory xnp and state of airway (open/closed)

        \param ele              (i) the element those matrix is calculated
        \param epn              (i) nodal pressure at n
        \param dt               (i) timestep
        */
      void EvaluateCollapse(RedAirway* ele, CORE::LINALG::SerialDenseVector& epn,
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
          DRT::Discretization& discretization, std::vector<int>& lm,
          CORE::LINALG::SerialDenseVector& rhs, Teuchos::RCP<MAT::Material> mat) override;

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
          DRT::Discretization& discretization, std::vector<int>& lm,
          CORE::LINALG::SerialDenseVector& radii_in, CORE::LINALG::SerialDenseVector& radii_out,
          Teuchos::RCP<const MAT::Material> material) override;

      /*!
       \Essential functions to compute the results of essentail matrices
      */
      void CalcFlowRates(RedAirway* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, std::vector<int>& lm,
          Teuchos::RCP<MAT::Material> mat) override;

      /*!
       \Essential functions to compute the volume of elements
      */
      void CalcElemVolume(RedAirway* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, std::vector<int>& lm,
          Teuchos::RCP<MAT::Material> mat) override;

      /*!
       \Essential functions to compute the volume mixing and  flowing into a junction
      */
      void GetJunctionVolumeMix(RedAirway* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization,
          CORE::LINALG::SerialDenseVector& junctionVolumeMix_np, std::vector<int>& lm,
          Teuchos::RCP<MAT::Material> material) override;

      /*!
       \Essential functions to compute the volume mixing and  flowing into a junction
      */
      void SolveScatra(RedAirway* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, CORE::LINALG::SerialDenseVector& scatranp,
          CORE::LINALG::SerialDenseVector& volumeMix_np, std::vector<int>& lm,
          Teuchos::RCP<MAT::Material> material) override;

      /*!
       \Essential functions to compute the volume mixing and  flowing into a junction
      */
      void SolveScatraBifurcations(RedAirway* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, CORE::LINALG::SerialDenseVector& scatra_np,
          CORE::LINALG::SerialDenseVector& volumeMix_np, std::vector<int>& lm,
          Teuchos::RCP<MAT::Material> material) override;
      /*!
       \Essential functions to evaluate the coupled results
      */
      void GetCoupledValues(RedAirway* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, std::vector<int>& lm,
          Teuchos::RCP<MAT::Material> material) override;

      void CalcCFL(RedAirway* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, std::vector<int>& lm,
          Teuchos::RCP<MAT::Material> material) override;

      void UpdateScatra(RedAirway* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, std::vector<int>& lm,
          Teuchos::RCP<MAT::Material> material) override;

      void UpdateElem12Scatra(RedAirway* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, std::vector<int>& lm,
          Teuchos::RCP<MAT::Material> material) override;

      void EvalPO2FromScatra(RedAirway* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, std::vector<int>& lm,
          Teuchos::RCP<MAT::Material> material) override;

      void EvalNodalEssentialValues(RedAirway* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, CORE::LINALG::SerialDenseVector& nodal_surface,
          CORE::LINALG::SerialDenseVector& nodal_volume,
          CORE::LINALG::SerialDenseVector& nodal_avg_scatra, std::vector<int>& lm,
          Teuchos::RCP<MAT::Material> material) override;


     private:
    };

  }  // namespace ELEMENTS
}  // namespace DRT

BACI_NAMESPACE_CLOSE

#endif
