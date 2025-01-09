// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_RED_AIRWAYS_ACINUS_IMPL_HPP
#define FOUR_C_RED_AIRWAYS_ACINUS_IMPL_HPP

#include "4C_config.hpp"

#include "4C_fem_general_utils_local_connectivity_matrices.hpp"
#include "4C_linalg_fixedsizematrix.hpp"
#include "4C_red_airways_elementbase.hpp"

FOUR_C_NAMESPACE_OPEN


namespace Discret
{
  namespace Elements
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
      virtual int evaluate(RedAcinus* ele, Teuchos::ParameterList& params,
          Core::FE::Discretization& discretization, std::vector<int>& lm,
          Core::LinAlg::SerialDenseMatrix& elemat1_epetra,
          Core::LinAlg::SerialDenseMatrix& elemat2_epetra,
          Core::LinAlg::SerialDenseVector& elevec1_epetra,
          Core::LinAlg::SerialDenseVector& elevec2_epetra,
          Core::LinAlg::SerialDenseVector& elevec3_epetra,
          std::shared_ptr<Core::Mat::Material> mat) = 0;

      virtual void initial(RedAcinus* ele, Teuchos::ParameterList& params,
          Core::FE::Discretization& discretization, std::vector<int>& lm,
          std::shared_ptr<const Core::Mat::Material> material) = 0;

      virtual void evaluate_terminal_bc(RedAcinus* ele, Teuchos::ParameterList& params,
          Core::FE::Discretization& discretization, std::vector<int>& lm,
          Core::LinAlg::SerialDenseVector& elevec1_epetra,
          std::shared_ptr<Core::Mat::Material> mat) = 0;

      virtual void calc_flow_rates(RedAcinus* ele, Teuchos::ParameterList& params,
          Core::FE::Discretization& discretization, std::vector<int>& lm,
          std::shared_ptr<Core::Mat::Material> mat) = 0;

      virtual void calc_elem_volume(RedAcinus* ele, Teuchos::ParameterList& params,
          Core::FE::Discretization& discretization, std::vector<int>& lm,
          std::shared_ptr<Core::Mat::Material> mat) = 0;

      virtual void get_coupled_values(RedAcinus* ele, Teuchos::ParameterList& params,
          Core::FE::Discretization& discretization, std::vector<int>& lm,
          std::shared_ptr<Core::Mat::Material> material) = 0;

      /// Internal implementation class for acinus element
      static RedAcinusImplInterface* impl(Discret::Elements::RedAcinus* acinus);
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

    template <Core::FE::CellType distype>
    class AcinusImpl : public RedAcinusImplInterface
    {
     public:
      /// Constructor
      explicit AcinusImpl();

      //! number of nodes
      static constexpr int iel = Core::FE::num_nodes<distype>;


      /// Evaluate
      /*!
        The evaluate function for the general acinus case.
      */
      int evaluate(RedAcinus* ele, Teuchos::ParameterList& params,
          Core::FE::Discretization& discretization, std::vector<int>& lm,
          Core::LinAlg::SerialDenseMatrix& elemat1_epetra,
          Core::LinAlg::SerialDenseMatrix& elemat2_epetra,
          Core::LinAlg::SerialDenseVector& elevec1_epetra,
          Core::LinAlg::SerialDenseVector& elevec2_epetra,
          Core::LinAlg::SerialDenseVector& elevec3_epetra,
          std::shared_ptr<Core::Mat::Material> mat) override;

      void evaluate_terminal_bc(RedAcinus* ele, Teuchos::ParameterList& params,
          Core::FE::Discretization& discretization, std::vector<int>& lm,
          Core::LinAlg::SerialDenseVector& rhs, std::shared_ptr<Core::Mat::Material> mat) override;

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
      void initial(RedAcinus* ele, Teuchos::ParameterList& params,
          Core::FE::Discretization& discretization, std::vector<int>& lm,
          std::shared_ptr<const Core::Mat::Material> material) override;

      /*!
        \Essential functions to compute the results of essential matrices
      */
      void calc_flow_rates(RedAcinus* ele, Teuchos::ParameterList& params,
          Core::FE::Discretization& discretization, std::vector<int>& lm,
          std::shared_ptr<Core::Mat::Material> mat) override;

      /*!
        \Essential functions to compute the volume of an element
      */
      void calc_elem_volume(RedAcinus* ele, Teuchos::ParameterList& params,
          Core::FE::Discretization& discretization, std::vector<int>& lm,
          std::shared_ptr<Core::Mat::Material> mat) override;

      /*!
        \Essential functions to evaluate the coupled results
      */
      void get_coupled_values(RedAcinus* ele, Teuchos::ParameterList& params,
          Core::FE::Discretization& discretization, std::vector<int>& lm,
          std::shared_ptr<Core::Mat::Material> material) override;

     private:
    };
  }  // namespace Elements
}  // namespace Discret


FOUR_C_NAMESPACE_CLOSE

#endif
