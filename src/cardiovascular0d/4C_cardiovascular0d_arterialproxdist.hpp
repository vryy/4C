// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_CARDIOVASCULAR0D_ARTERIALPROXDIST_HPP
#define FOUR_C_CARDIOVASCULAR0D_ARTERIALPROXDIST_HPP

#include "4C_config.hpp"

#include "4C_cardiovascular0d.hpp"
#include "4C_fem_general_utils_integration.hpp"
#include "4C_inpar_cardiovascular0d.hpp"
#include "4C_linalg_vector.hpp"
#include "4C_utils_parameter_list.fwd.hpp"

#include <Epetra_FECrsMatrix.h>
#include <Epetra_Operator.h>
#include <Epetra_RowMatrix.h>

#include <memory>

FOUR_C_NAMESPACE_OPEN

// forward declarations
namespace Core::FE
{
  class Discretization;
}  // namespace Core::FE

namespace Core::LinAlg
{
  class SparseMatrix;
  class SparseOperator;
}  // namespace Core::LinAlg

namespace Utils
{
  class Cardiovascular0DArterialProxDist : public Cardiovascular0D

  {
   public:
    /*!
    \brief Constructor of a Cardiovascular0D based on conditions with a given name. It also
    takes care of the Cardiovascular0D IDs.
    */

    Cardiovascular0DArterialProxDist(std::shared_ptr<Core::FE::Discretization>
                                         discr,  ///< discretization where Cardiovascular0D lives on
        const std::string& conditionname,  ///< Name of condition to create Cardiovascular0D from
        std::vector<int>& curID            ///< current ID
    );



    /// initialization routine called by the manager ctor to get correct reference base values and
    /// activating the right conditions at the beginning
    void initialize(
        Teuchos::ParameterList&
            params,  ///< parameter list to communicate between elements and discretization
        std::shared_ptr<Core::LinAlg::Vector<double>>
            sysvec1,  ///< distributed vector that may be filled by
                      ///< assembly of element contributions
        std::shared_ptr<Core::LinAlg::Vector<double>>
            sysvec2  ///< distributed vector that may be filled by assembly of element contributions
        ) override;

    //! Evaluate routine to call from outside. In here the right action is determined and the
    //! #EvaluateCardiovascular0D routine is called
    void evaluate(
        Teuchos::ParameterList&
            params,  ///< parameter list to communicate between elements and discretization
        std::shared_ptr<Core::LinAlg::SparseMatrix> sysmat1,  ///< Cardiovascular0D stiffness matrix
        std::shared_ptr<Core::LinAlg::SparseOperator>
            sysmat2,  ///< Cardiovascular0D offdiagonal matrix dV/dd
        std::shared_ptr<Core::LinAlg::SparseOperator>
            sysmat3,  ///< Cardiovascular0D offdiagonal matrix dfext/dp
        std::shared_ptr<Core::LinAlg::Vector<double>>
            sysvec1,  ///< distributed vectors that may be filled by
                      ///< assembly of element contributions
        std::shared_ptr<Core::LinAlg::Vector<double>> sysvec2,
        std::shared_ptr<Core::LinAlg::Vector<double>> sysvec3,
        const std::shared_ptr<Core::LinAlg::Vector<double>> sysvec4,
        std::shared_ptr<Core::LinAlg::Vector<double>> sysvec5) override;

   private:
    // don't want = operator, cctor and destructor

    Cardiovascular0DArterialProxDist operator=(const Cardiovascular0DArterialProxDist& old);
    Cardiovascular0DArterialProxDist(const Cardiovascular0DArterialProxDist& old);



  };  // class
}  // namespace Utils

FOUR_C_NAMESPACE_CLOSE

#endif
