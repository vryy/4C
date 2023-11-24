/*---------------------------------------------------------------------*/
/*! \file

\brief Internal implementation of acinus_impl element


\level 3

*/
/*---------------------------------------------------------------------*/



#ifndef BACI_RED_AIRWAYS_ACINUS_IMPL_H
#define BACI_RED_AIRWAYS_ACINUS_IMPL_H

#include "baci_discretization_fem_general_utils_local_connectivity_matrices.H"
#include "baci_linalg_fixedsizematrix.H"
#include "baci_red_airways_elementbase.H"


namespace DRT
{
  namespace ELEMENTS
  {
    /// Interface base class for acinus_impl
    /*!
      This class exists to provide a common interface for all template
      versions of acinus_impl. The only function this class actually
      defines is Impl, which returns a pointer to the appropriate version
      of acinus_impl.
    */
    class RedAcinusImplInterface
    {
     public:
      /// Empty constructor
      RedAcinusImplInterface() {}
      /// Empty destructor
      virtual ~RedAcinusImplInterface() = default;  /// Evaluate the element
      virtual int Evaluate(RedAcinus* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, std::vector<int>& lm,
          CORE::LINALG::SerialDenseMatrix& elemat1_epetra,
          CORE::LINALG::SerialDenseMatrix& elemat2_epetra,
          CORE::LINALG::SerialDenseVector& elevec1_epetra,
          CORE::LINALG::SerialDenseVector& elevec2_epetra,
          CORE::LINALG::SerialDenseVector& elevec3_epetra, Teuchos::RCP<MAT::Material> mat) = 0;

      virtual void Initial(RedAcinus* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, std::vector<int>& lm,
          Teuchos::RCP<const MAT::Material> material) = 0;

      virtual void EvaluateTerminalBC(RedAcinus* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, std::vector<int>& lm,
          CORE::LINALG::SerialDenseVector& elevec1_epetra, Teuchos::RCP<MAT::Material> mat) = 0;

      virtual void CalcFlowRates(RedAcinus* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, std::vector<int>& lm,
          Teuchos::RCP<MAT::Material> mat) = 0;

      virtual void CalcElemVolume(RedAcinus* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, std::vector<int>& lm,
          Teuchos::RCP<MAT::Material> mat) = 0;

      virtual void GetCoupledValues(RedAcinus* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, std::vector<int>& lm,
          Teuchos::RCP<MAT::Material> material) = 0;

      virtual void GetJunctionVolumeMix(RedAcinus* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, CORE::LINALG::SerialDenseVector& volumeMix_np,
          std::vector<int>& lm, Teuchos::RCP<MAT::Material> material) = 0;

      virtual void SolveScatra(RedAcinus* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, CORE::LINALG::SerialDenseVector& scatra_np,
          CORE::LINALG::SerialDenseVector& volumeMix_np, std::vector<int>& lm,
          Teuchos::RCP<MAT::Material> material) = 0;

      virtual void SolveScatraBifurcations(RedAcinus* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, CORE::LINALG::SerialDenseVector& scatranp,
          CORE::LINALG::SerialDenseVector& volumeMix_np, std::vector<int>& lm,
          Teuchos::RCP<MAT::Material> material) = 0;

      virtual void UpdateScatra(RedAcinus* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, std::vector<int>& lm,
          Teuchos::RCP<MAT::Material> material) = 0;

      virtual void UpdateElem12Scatra(RedAcinus* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, std::vector<int>& lm,
          Teuchos::RCP<MAT::Material> material) = 0;

      virtual void EvalPO2FromScatra(RedAcinus* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, std::vector<int>& lm,
          Teuchos::RCP<MAT::Material> material) = 0;

      virtual void EvalNodalEssentialValues(RedAcinus* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, CORE::LINALG::SerialDenseVector& nodal_surface,
          CORE::LINALG::SerialDenseVector& nodal_volume,
          CORE::LINALG::SerialDenseVector& nodal_flow, std::vector<int>& lm,
          Teuchos::RCP<MAT::Material> material) = 0;

      /// Internal implementation class for acinus element
      static RedAcinusImplInterface* Impl(DRT::ELEMENTS::RedAcinus* acinus);
    };


    /// Internal acinus implementation
    /*!
      This internal class keeps all the working arrays needed to
      calculate the acinus element. Additionally the method Sysmat()
      provides a clean and fast element implementation.

      <h3>Purpose</h3>

      \author ismail
      \date 01/13
    */

    template <CORE::FE::CellType distype>
    class AcinusImpl : public RedAcinusImplInterface
    {
     public:
      /// Constructor
      explicit AcinusImpl();

      //! number of nodes
      static const int iel = CORE::FE::num_nodes<distype>;


      /// Evaluate
      /*!
        The evaluate function for the general acinus case.
      */
      int Evaluate(RedAcinus* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, std::vector<int>& lm,
          CORE::LINALG::SerialDenseMatrix& elemat1_epetra,
          CORE::LINALG::SerialDenseMatrix& elemat2_epetra,
          CORE::LINALG::SerialDenseVector& elevec1_epetra,
          CORE::LINALG::SerialDenseVector& elevec2_epetra,
          CORE::LINALG::SerialDenseVector& elevec3_epetra,
          Teuchos::RCP<MAT::Material> mat) override;

      void EvaluateTerminalBC(RedAcinus* ele, Teuchos::ParameterList& params,
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
        \param material         (i) acinus material/dimesion
        \param time             (i) current simulation time
        \param dt               (i) timestep
      */
      void Initial(RedAcinus* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, std::vector<int>& lm,
          Teuchos::RCP<const MAT::Material> material) override;

      /*!
        \Essential functions to compute the results of essentail matrices
      */
      void CalcFlowRates(RedAcinus* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, std::vector<int>& lm,
          Teuchos::RCP<MAT::Material> mat) override;

      /*!
        \Essential functions to compute the volume of an element
      */
      void CalcElemVolume(RedAcinus* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, std::vector<int>& lm,
          Teuchos::RCP<MAT::Material> mat) override;

      /*!
        \Essential functions to evaluate the coupled results
      */
      void GetCoupledValues(RedAcinus* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, std::vector<int>& lm,
          Teuchos::RCP<MAT::Material> material) override;

      /*!
        \Essential functions to evaluate mixed volume flowing into a junction
      */
      void GetJunctionVolumeMix(RedAcinus* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, CORE::LINALG::SerialDenseVector& volumeMix_np,
          std::vector<int>& lm, Teuchos::RCP<MAT::Material> material) override;

      /*!
        \Essential functions to solve the forward scatra flow
      */
      void SolveScatra(RedAcinus* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, CORE::LINALG::SerialDenseVector& scatranp,
          CORE::LINALG::SerialDenseVector& volumeMix_np, std::vector<int>& lm,
          Teuchos::RCP<MAT::Material> material) override;

      void SolveScatraBifurcations(RedAcinus* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, CORE::LINALG::SerialDenseVector& scatranp,
          CORE::LINALG::SerialDenseVector& volumeMix_np, std::vector<int>& lm,
          Teuchos::RCP<MAT::Material> material) override;

      void UpdateScatra(RedAcinus* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, std::vector<int>& lm,
          Teuchos::RCP<MAT::Material> material) override;

      void UpdateElem12Scatra(RedAcinus* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, std::vector<int>& lm,
          Teuchos::RCP<MAT::Material> material) override;

      void EvalPO2FromScatra(RedAcinus* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, std::vector<int>& lm,
          Teuchos::RCP<MAT::Material> material) override;

      void EvalNodalEssentialValues(RedAcinus* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, CORE::LINALG::SerialDenseVector& nodal_surface,
          CORE::LINALG::SerialDenseVector& nodal_volume,
          CORE::LINALG::SerialDenseVector& nodal_avg_scatra, std::vector<int>& lm,
          Teuchos::RCP<MAT::Material> material) override;


     private:
    };
  }  // namespace ELEMENTS
}  // namespace DRT


#endif
