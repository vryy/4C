/*---------------------------------------------------------------------*/
/*! \file

\brief Internal implementation of inter_acinar_dep_impl element


\level 3

*/
/*---------------------------------------------------------------------*/



#ifndef FOUR_C_RED_AIRWAYS_INTERACINARDEP_IMPL_HPP
#define FOUR_C_RED_AIRWAYS_INTERACINARDEP_IMPL_HPP

#include "4C_config.hpp"

#include "4C_discretization_fem_general_utils_local_connectivity_matrices.hpp"
#include "4C_linalg_fixedsizematrix.hpp"
#include "4C_red_airways_elementbase.hpp"

FOUR_C_NAMESPACE_OPEN


namespace DRT
{
  namespace ELEMENTS
  {
    /// Interface base class for inter_acinar_dep_impl
    /*!
      This class exists to provide a common interface for all template
      versions of inter_acinar_dep_impl. The only function
      this class actually defines is Impl, which returns a pointer to
      the appropriate version of inter_acinar_dep_impl.
     */
    class RedInterAcinarDepImplInterface
    {
     public:
      /// Empty constructor
      RedInterAcinarDepImplInterface() {}
      /// Empty destructor
      virtual ~RedInterAcinarDepImplInterface() = default;  /// Evaluate the element
      virtual int Evaluate(RedInterAcinarDep* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, std::vector<int>& lm,
          CORE::LINALG::SerialDenseMatrix& elemat1_epetra,
          CORE::LINALG::SerialDenseMatrix& elemat2_epetra,
          CORE::LINALG::SerialDenseVector& elevec1_epetra,
          CORE::LINALG::SerialDenseVector& elevec2_epetra,
          CORE::LINALG::SerialDenseVector& elevec3_epetra,
          Teuchos::RCP<CORE::MAT::Material> mat) = 0;

      virtual void Initial(RedInterAcinarDep* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, std::vector<int>& lm,
          CORE::LINALG::SerialDenseVector& elevec3_epetra,
          Teuchos::RCP<const CORE::MAT::Material> material) = 0;

      virtual void EvaluateTerminalBC(RedInterAcinarDep* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, std::vector<int>& lm,
          CORE::LINALG::SerialDenseVector& elevec1_epetra,
          Teuchos::RCP<CORE::MAT::Material> mat) = 0;

      virtual void CalcFlowRates(RedInterAcinarDep* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, CORE::LINALG::SerialDenseVector& a_volumen,
          CORE::LINALG::SerialDenseVector& a_volumenp, std::vector<int>& lm,
          Teuchos::RCP<CORE::MAT::Material> mat) = 0;

      virtual void GetCoupledValues(RedInterAcinarDep* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, std::vector<int>& lm,
          Teuchos::RCP<CORE::MAT::Material> material) = 0;

      /// Internal implementation class for inter-acinar linker element
      static RedInterAcinarDepImplInterface* Impl(DRT::ELEMENTS::RedInterAcinarDep* acinus);
    };


    /// Internal inter-acinar linker implementation
    /*!
      This internal class keeps all the working arrays needed to
      calculate the inter-acinar linker element. Additionally the
      method Sysmat() provides a clean and fast element implementation.

      <h3>Purpose</h3>

      \author ismail
      \date 01/09
    */

    template <CORE::FE::CellType distype>
    class InterAcinarDepImpl : public RedInterAcinarDepImplInterface
    {
     public:
      /// Constructor
      explicit InterAcinarDepImpl();

      //! number of nodes
      static constexpr int iel = CORE::FE::num_nodes<distype>;


      /// Evaluate
      /*!
        The evaluate function for the general inter-acinar linker case.
       */
      int Evaluate(RedInterAcinarDep* ele, Teuchos::ParameterList& params,
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
      void Sysmat(std::vector<double>& ial, CORE::LINALG::SerialDenseMatrix& sysmat,
          CORE::LINALG::SerialDenseVector& rhs);


      void EvaluateTerminalBC(RedInterAcinarDep* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, std::vector<int>& lm,
          CORE::LINALG::SerialDenseVector& rhs,
          Teuchos::RCP<CORE::MAT::Material> material) override;

      /*!
        \brief get the initial values of the degrees of freedom at the node

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
      void Initial(RedInterAcinarDep* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, std::vector<int>& lm,
          CORE::LINALG::SerialDenseVector& n_intr_acn_l,
          Teuchos::RCP<const CORE::MAT::Material> material) override;

      /*!
       \Essential functions to compute the results of essential matrices
      */
      void CalcFlowRates(RedInterAcinarDep* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, CORE::LINALG::SerialDenseVector& a_volumen_strain_np,
          CORE::LINALG::SerialDenseVector& a_volumenp, std::vector<int>& lm,
          Teuchos::RCP<CORE::MAT::Material> mat) override{};

      /*!
       \Essential functions to evaluate the coupled results
      */
      void GetCoupledValues(RedInterAcinarDep* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, std::vector<int>& lm,
          Teuchos::RCP<CORE::MAT::Material> material) override{};

     private:
    };

  }  // namespace ELEMENTS
}  // namespace DRT

FOUR_C_NAMESPACE_CLOSE

#endif
