/*---------------------------------------------------------------------*/
/*! \file

\brief Incomplete! - Purpose: Internal implementation of RedAirBloodScatraLine3 element


\level 3

*/
/*---------------------------------------------------------------------*/



#ifndef FOUR_C_RED_AIRWAYS_AIR_BLOOD_SCATRALINE3_IMPL_HPP
#define FOUR_C_RED_AIRWAYS_AIR_BLOOD_SCATRALINE3_IMPL_HPP

#include "4C_config.hpp"

#include "4C_fem_general_utils_local_connectivity_matrices.hpp"
#include "4C_linalg_fixedsizematrix.hpp"
#include "4C_red_airways_elementbase.hpp"

FOUR_C_NAMESPACE_OPEN


namespace Discret
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
    class RedAirBloodScatraLine3ImplInterface
    {
     public:
      /// Empty constructor
      RedAirBloodScatraLine3ImplInterface() {}
      /// Empty destructor
      virtual ~RedAirBloodScatraLine3ImplInterface() = default;  /// Evaluate the element
      virtual int Evaluate(RedAirBloodScatraLine3* ele, Teuchos::ParameterList& params,
          Core::FE::Discretization& discretization, std::vector<int>& lm,
          Core::LinAlg::SerialDenseMatrix& elemat1_epetra,
          Core::LinAlg::SerialDenseMatrix& elemat2_epetra,
          Core::LinAlg::SerialDenseVector& elevec1_epetra,
          Core::LinAlg::SerialDenseVector& elevec2_epetra,
          Core::LinAlg::SerialDenseVector& elevec3_epetra,
          Teuchos::RCP<Core::Mat::Material> mat) = 0;

      virtual void Initial(RedAirBloodScatraLine3* ele, Teuchos::ParameterList& params,
          Core::FE::Discretization& discretization, std::vector<int>& lm,
          Teuchos::RCP<const Core::Mat::Material> material) = 0;

      virtual void EvaluateTerminalBC(RedAirBloodScatraLine3* ele, Teuchos::ParameterList& params,
          Core::FE::Discretization& discretization, std::vector<int>& lm,
          Core::LinAlg::SerialDenseVector& elevec1_epetra,
          Teuchos::RCP<Core::Mat::Material> mat) = 0;

      virtual void CalcFlowRates(RedAirBloodScatraLine3* ele, Teuchos::ParameterList& params,
          Core::FE::Discretization& discretization, Core::LinAlg::SerialDenseVector& a_volumen,
          Core::LinAlg::SerialDenseVector& a_volumenp, std::vector<int>& lm,
          Teuchos::RCP<Core::Mat::Material> mat) = 0;

      virtual void GetCoupledValues(RedAirBloodScatraLine3* ele, Teuchos::ParameterList& params,
          Core::FE::Discretization& discretization, std::vector<int>& lm,
          Teuchos::RCP<Core::Mat::Material> material) = 0;

      virtual void solve_blood_air_transport(RedAirBloodScatraLine3* ele,
          Core::LinAlg::SerialDenseVector& dscatra, Core::LinAlg::SerialDenseVector& dvo2,
          Teuchos::ParameterList& params, Core::FE::Discretization& discretization,
          std::vector<int>& lm, Teuchos::RCP<Core::Mat::Material> material) = 0;

      /// Internal implementation class for acinus element
      static RedAirBloodScatraLine3ImplInterface* Impl(
          Discret::ELEMENTS::RedAirBloodScatraLine3* red_acinus);
    };


    /// Internal acinus implementation
    /*!
      This internal class keeps all the working arrays needed to
      calculate the acinus element. Additionally the method Sysmat()
      provides a clean and fast element implementation.

      <h3>Purpose</h3>

      \author ismail
      \date 6/13
    */

    template <Core::FE::CellType distype>
    class RedAirBloodScatraLine3Impl : public RedAirBloodScatraLine3ImplInterface
    {
     public:
      /// Constructor
      explicit RedAirBloodScatraLine3Impl();

      //! number of nodes
      static constexpr int iel = Core::FE::num_nodes<distype>;


      /// Evaluate
      /*!
        The evaluate function for the general acinus case.
       */
      int Evaluate(RedAirBloodScatraLine3* ele, Teuchos::ParameterList& params,
          Core::FE::Discretization& discretization, std::vector<int>& lm,
          Core::LinAlg::SerialDenseMatrix& elemat1_epetra,
          Core::LinAlg::SerialDenseMatrix& elemat2_epetra,
          Core::LinAlg::SerialDenseVector& elevec1_epetra,
          Core::LinAlg::SerialDenseVector& elevec2_epetra,
          Core::LinAlg::SerialDenseVector& elevec3_epetra,
          Teuchos::RCP<Core::Mat::Material> mat) override;

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
      void Sysmat(RedAirBloodScatraLine3* ele, Core::LinAlg::SerialDenseVector& epnp,
          Core::LinAlg::SerialDenseVector& epn, Core::LinAlg::SerialDenseVector& epnm,
          Core::LinAlg::SerialDenseMatrix& estif, Core::LinAlg::SerialDenseVector& eforce,
          Teuchos::RCP<const Core::Mat::Material> material, Teuchos::ParameterList& params,
          double time, double dt){};


      void EvaluateTerminalBC(RedAirBloodScatraLine3* ele, Teuchos::ParameterList& params,
          Core::FE::Discretization& disctretization, std::vector<int>& lm,
          Core::LinAlg::SerialDenseVector& elevec1_epetra,
          Teuchos::RCP<Core::Mat::Material> mat) override{};

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
      void Initial(RedAirBloodScatraLine3* ele, Teuchos::ParameterList& params,
          Core::FE::Discretization& discretization, std::vector<int>& lm,
          Teuchos::RCP<const Core::Mat::Material> material) override;

      /*!
       \Essential functions to compute the results of essentail matrices
      */
      void CalcFlowRates(RedAirBloodScatraLine3* ele, Teuchos::ParameterList& params,
          Core::FE::Discretization& discretization,
          Core::LinAlg::SerialDenseVector& a_volumen_strain_np,
          Core::LinAlg::SerialDenseVector& a_volumenp, std::vector<int>& lm,
          Teuchos::RCP<Core::Mat::Material> mat) override{};

      /*!
       \Essential functions to evaluate the coupled results
      */
      void GetCoupledValues(RedAirBloodScatraLine3* ele, Teuchos::ParameterList& params,
          Core::FE::Discretization& discretization, std::vector<int>& lm,
          Teuchos::RCP<Core::Mat::Material> material) override;

      void solve_blood_air_transport(RedAirBloodScatraLine3* ele,
          Core::LinAlg::SerialDenseVector& dscatra, Core::LinAlg::SerialDenseVector& dvo2,
          Teuchos::ParameterList& params, Core::FE::Discretization& discretization,
          std::vector<int>& lm, Teuchos::RCP<Core::Mat::Material> material) override{};

     private:
    };
  }  // namespace ELEMENTS
}  // namespace Discret

FOUR_C_NAMESPACE_CLOSE

#endif
