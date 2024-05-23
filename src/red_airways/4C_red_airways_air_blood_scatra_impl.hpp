/*---------------------------------------------------------------------*/
/*! \file

\brief Incomplete! - Purpose: Internal implementation of RedAirBloodScatra element


\level 3

*/
/*---------------------------------------------------------------------*/



#ifndef FOUR_C_RED_AIRWAYS_AIR_BLOOD_SCATRA_IMPL_HPP
#define FOUR_C_RED_AIRWAYS_AIR_BLOOD_SCATRA_IMPL_HPP

#include "4C_config.hpp"

#include "4C_discretization_fem_general_utils_local_connectivity_matrices.hpp"
#include "4C_linalg_fixedsizematrix.hpp"
#include "4C_red_airways_elementbase.hpp"

FOUR_C_NAMESPACE_OPEN


namespace DRT
{
  namespace ELEMENTS
  {
    /// Interface base class for red_air_blood_scatra_impl
    /*!
      This class exists to provide a common interface for all template
      versions of red_air_blood_scatra_impl. The only function
      this class actually defines is Impl, which returns a pointer to
      the appropriate version of ired_air_blood_scatra_impl.
     */
    class RedAirBloodScatraImplInterface
    {
     public:
      /// Empty constructor
      RedAirBloodScatraImplInterface() {}
      /// Empty destructor
      virtual ~RedAirBloodScatraImplInterface() = default;  /// Evaluate the element
      virtual int Evaluate(RedAirBloodScatra* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, std::vector<int>& lm,
          CORE::LINALG::SerialDenseMatrix& elemat1_epetra,
          CORE::LINALG::SerialDenseMatrix& elemat2_epetra,
          CORE::LINALG::SerialDenseVector& elevec1_epetra,
          CORE::LINALG::SerialDenseVector& elevec2_epetra,
          CORE::LINALG::SerialDenseVector& elevec3_epetra,
          Teuchos::RCP<CORE::MAT::Material> mat) = 0;

      virtual void Initial(RedAirBloodScatra* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, std::vector<int>& lm,
          Teuchos::RCP<const CORE::MAT::Material> material) = 0;

      virtual void EvaluateTerminalBC(RedAirBloodScatra* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, std::vector<int>& lm,
          CORE::LINALG::SerialDenseVector& elevec1_epetra,
          Teuchos::RCP<CORE::MAT::Material> mat) = 0;

      virtual void CalcFlowRates(RedAirBloodScatra* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, CORE::LINALG::SerialDenseVector& a_volumen,
          CORE::LINALG::SerialDenseVector& a_volumenp, std::vector<int>& lm,
          Teuchos::RCP<CORE::MAT::Material> mat) = 0;

      virtual void GetCoupledValues(RedAirBloodScatra* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, std::vector<int>& lm,
          Teuchos::RCP<CORE::MAT::Material> material) = 0;

      virtual void solve_blood_air_transport(RedAirBloodScatra* ele,
          CORE::LINALG::SerialDenseVector& dscatra, CORE::LINALG::SerialDenseVector& dvO2,
          CORE::LINALG::SerialDenseVector& scatraAcinus, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, std::vector<int>& lm,
          Teuchos::RCP<CORE::MAT::Material> material) = 0;

      /// Internal implementation class for acinus element
      static RedAirBloodScatraImplInterface* Impl(DRT::ELEMENTS::RedAirBloodScatra* red_acinus);
    };


    /// Internal acinus implementation
    /*!
      This internal class keeps all the working arrays needed to
      calculate the acinus element. Additionally the method Sysmat()
      provides a clean and fast element implementation.

      <h3>Purpose</h3>

      \author ismail
      \date 01/09
    */

    template <CORE::FE::CellType distype>
    class RedAirBloodScatraImpl : public RedAirBloodScatraImplInterface
    {
     public:
      /// Constructor
      explicit RedAirBloodScatraImpl();

      //! number of nodes
      static constexpr int iel = CORE::FE::num_nodes<distype>;


      /// Evaluate
      /*!
        The evaluate function for the general acinus case.
       */
      int Evaluate(RedAirBloodScatra* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, std::vector<int>& lm,
          CORE::LINALG::SerialDenseMatrix& elemat1_epetra,
          CORE::LINALG::SerialDenseMatrix& elemat2_epetra,
          CORE::LINALG::SerialDenseVector& elevec1_epetra,
          CORE::LINALG::SerialDenseVector& elevec2_epetra,
          CORE::LINALG::SerialDenseVector& elevec3_epetra,
          Teuchos::RCP<CORE::MAT::Material> mat) override;

      /*!
        \brief calculate element matrix and rhs

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
      void Sysmat(RedAirBloodScatra* ele, CORE::LINALG::SerialDenseVector& epnp,
          CORE::LINALG::SerialDenseVector& epn, CORE::LINALG::SerialDenseVector& epnm,
          CORE::LINALG::SerialDenseMatrix& estif, CORE::LINALG::SerialDenseVector& eforce,
          Teuchos::RCP<const CORE::MAT::Material> material, Teuchos::ParameterList& params,
          double time, double dt){};


      void EvaluateTerminalBC(RedAirBloodScatra* ele, Teuchos::ParameterList& params,
          DRT::Discretization& disctretization, std::vector<int>& lm,
          CORE::LINALG::SerialDenseVector& elevec1_epetra,
          Teuchos::RCP<CORE::MAT::Material> mat) override{};

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
      void Initial(RedAirBloodScatra* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, std::vector<int>& lm,
          Teuchos::RCP<const CORE::MAT::Material> material) override;

      /*!
       \Essential functions to compute the results of essentail matrices
      */
      void CalcFlowRates(RedAirBloodScatra* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, CORE::LINALG::SerialDenseVector& a_volumen_strain_np,
          CORE::LINALG::SerialDenseVector& a_volumenp, std::vector<int>& lm,
          Teuchos::RCP<CORE::MAT::Material> mat) override{};

      /*!
       \Essential functions to evaluate the coupled results
      */
      void GetCoupledValues(RedAirBloodScatra* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, std::vector<int>& lm,
          Teuchos::RCP<CORE::MAT::Material> material) override;

      void solve_blood_air_transport(RedAirBloodScatra* ele,
          CORE::LINALG::SerialDenseVector& dscatra, CORE::LINALG::SerialDenseVector& dvO2,
          CORE::LINALG::SerialDenseVector& scatra_acinus, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, std::vector<int>& lm,
          Teuchos::RCP<CORE::MAT::Material> material) override;

     private:
    };
  }  // namespace ELEMENTS
}  // namespace DRT

FOUR_C_NAMESPACE_CLOSE

#endif
