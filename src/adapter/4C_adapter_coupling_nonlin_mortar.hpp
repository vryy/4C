// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_ADAPTER_COUPLING_NONLIN_MORTAR_HPP
#define FOUR_C_ADAPTER_COUPLING_NONLIN_MORTAR_HPP

/*---------------------------------------------------------------------*
 | headers                                                 farah 10/14 |
 *---------------------------------------------------------------------*/
#include "4C_config.hpp"

#include "4C_coupling_adapter_mortar.hpp"
#include "4C_fem_condition.hpp"
#include "4C_linalg_vector.hpp"
#include "4C_utils_exceptions.hpp"

#include <Epetra_Map.h>
#include <Teuchos_ParameterListAcceptorDefaultBase.hpp>

#include <memory>

FOUR_C_NAMESPACE_OPEN

/*---------------------------------------------------------------------*
 | forward declarations                                    farah 10/14 |
 *---------------------------------------------------------------------*/
namespace Core::FE
{
  class Discretization;
}  // namespace Core::FE

namespace Core::Elements
{
  class Element;
}

namespace Core::Nodes
{
  class Node;
}

namespace CONTACT
{
  class Interface;
}

namespace Core::LinAlg
{
  class SparseMatrix;
}

namespace Adapter
{
  class CouplingNonLinMortar : public Coupling::Adapter::CouplingMortar
  {
   public:
    /**
     * Construct nonlinear coupling with basic parameters. The remaining information is passed in
     * setup().
     */
    CouplingNonLinMortar(int spatial_dimension, Teuchos::ParameterList mortar_coupling_params,
        Teuchos::ParameterList contact_dynamic_params,
        Core::FE::ShapeFunctionType shape_function_type);

    /*!
    \brief initialize routine

    */
    virtual void setup(std::shared_ptr<Core::FE::Discretization> masterdis,
        std::shared_ptr<Core::FE::Discretization> slavedis, std::vector<int> coupleddof,
        const std::string& couplingcond);

    virtual void setup_spring_dashpot(std::shared_ptr<Core::FE::Discretization> masterdis,
        std::shared_ptr<Core::FE::Discretization> slavedis,
        std::shared_ptr<Core::Conditions::Condition> spring, const int coupling_id, MPI_Comm comm);

    virtual void integrate_lin_d(const std::string& statename,
        const std::shared_ptr<Core::LinAlg::Vector<double>> vec,
        const std::shared_ptr<Core::LinAlg::Vector<double>> veclm);

    virtual void integrate_lin_dm(const std::string& statename,
        const std::shared_ptr<Core::LinAlg::Vector<double>> vec,
        const std::shared_ptr<Core::LinAlg::Vector<double>> veclm);

    virtual void integrate_all(const std::string& statename,
        const std::shared_ptr<Core::LinAlg::Vector<double>> vec,
        const std::shared_ptr<Core::LinAlg::Vector<double>> veclm);

    virtual void evaluate_sliding(const std::string& statename,
        const std::shared_ptr<Core::LinAlg::Vector<double>> vec,
        const std::shared_ptr<Core::LinAlg::Vector<double>> veclm);

    virtual void print_interface(std::ostream& os);

    virtual std::shared_ptr<Core::LinAlg::SparseMatrix> d_lin_matrix()
    {
      if (DLin_ == nullptr) FOUR_C_THROW("ERROR: DLin Matrix is null pointer!");
      return DLin_;
    };

    virtual std::shared_ptr<Core::LinAlg::SparseMatrix> m_lin_matrix()
    {
      if (MLin_ == nullptr) FOUR_C_THROW("ERROR: MLin Matrix is null pointer!");
      return MLin_;
    };

    virtual std::shared_ptr<Core::LinAlg::SparseMatrix> h_matrix()
    {
      if (H_ == nullptr) FOUR_C_THROW("ERROR: H Matrix is null pointer!");
      return H_;
    };

    virtual std::shared_ptr<Core::LinAlg::SparseMatrix> t_matrix()
    {
      if (T_ == nullptr) FOUR_C_THROW("ERROR: T Matrix is null pointer!");
      return T_;
    };

    virtual std::shared_ptr<Core::LinAlg::SparseMatrix> n_matrix()
    {
      if (N_ == nullptr) FOUR_C_THROW("ERROR: N Matrix is null pointer!");
      return N_;
    };

    // create projection operator Dinv*M
    void create_p() override;

    virtual std::shared_ptr<Core::LinAlg::Vector<double>> gap()
    {
      if (gap_ == nullptr) FOUR_C_THROW("ERROR: gap vector is null pointer!");
      return gap_;
    };

    /// the mortar interface itself
    std::shared_ptr<CONTACT::Interface> interface() const { return interface_; }

   protected:
    /*!
    \brief Read Mortar Condition

    */
    virtual void read_mortar_condition(std::shared_ptr<Core::FE::Discretization> masterdis,
        std::shared_ptr<Core::FE::Discretization> slavedis, std::vector<int> coupleddof,
        const std::string& couplingcond, Teuchos::ParameterList& input,
        std::map<int, Core::Nodes::Node*>& mastergnodes,
        std::map<int, Core::Nodes::Node*>& slavegnodes,
        std::map<int, std::shared_ptr<Core::Elements::Element>>& masterelements,
        std::map<int, std::shared_ptr<Core::Elements::Element>>& slaveelements);

    /*!
    \brief Add Mortar Nodes

    */
    virtual void add_mortar_nodes(std::shared_ptr<Core::FE::Discretization> masterdis,
        std::shared_ptr<Core::FE::Discretization> slavedis, std::vector<int> coupleddof,
        Teuchos::ParameterList& input, std::map<int, Core::Nodes::Node*>& mastergnodes,
        std::map<int, Core::Nodes::Node*>& slavegnodes,
        std::map<int, std::shared_ptr<Core::Elements::Element>>& masterelements,
        std::map<int, std::shared_ptr<Core::Elements::Element>>& slaveelements,
        std::shared_ptr<CONTACT::Interface>& interface, int numcoupleddof);

    /*!
    \brief Add Mortar Elements

    */
    virtual void add_mortar_elements(std::shared_ptr<Core::FE::Discretization> masterdis,
        std::shared_ptr<Core::FE::Discretization> slavedis, Teuchos::ParameterList& input,
        std::map<int, std::shared_ptr<Core::Elements::Element>>& masterelements,
        std::map<int, std::shared_ptr<Core::Elements::Element>>& slaveelements,
        std::shared_ptr<CONTACT::Interface>& interface, int numcoupleddof);

    /*!
    \brief complete interface, store as internal variable
           store maps as internal variable and do parallel redist.

    */
    virtual void complete_interface(std::shared_ptr<Core::FE::Discretization> masterdis,
        std::shared_ptr<CONTACT::Interface>& interface);

    /*!
    \brief initialize matrices (interla variables)

    */
    virtual void init_matrices();

    /*!
    \brief create strategy object if required

    */
    virtual void create_strategy(std::shared_ptr<Core::FE::Discretization> masterdis,
        std::shared_ptr<Core::FE::Discretization> slavedis, Teuchos::ParameterList& input,
        int numcoupleddof);

    /*!
    \brief transform back to initial parallel distribution

    */
    virtual void matrix_row_col_transform();

    /// check setup call
    const bool& is_setup() const { return issetup_; };

    /// check init and setup call
    void check_setup() const override
    {
      if (!is_setup()) FOUR_C_THROW("ERROR: Call setup() first!");
    }

   protected:
    bool issetup_;   ///< check for setup
    MPI_Comm comm_;  ///< communicator
    int myrank_;     ///< my proc id

    std::shared_ptr<Epetra_Map>
        slavenoderowmap_;  ///< map of slave row nodes (after parallel redist.)
    std::shared_ptr<Epetra_Map>
        pslavenoderowmap_;  ///< map of slave row nodes (before parallel redist.)
    std::shared_ptr<Epetra_Map>
        smdofrowmap_;  ///< map of sm merged row dofs (after parallel redist.)
    std::shared_ptr<Epetra_Map>
        psmdofrowmap_;  ///< map of sm merged row dofs (before parallel redist.)

    std::shared_ptr<Core::LinAlg::SparseMatrix> DLin_;  ///< linearization of D matrix
    std::shared_ptr<Core::LinAlg::SparseMatrix> MLin_;  ///< linearization of M matrix

    std::shared_ptr<Core::LinAlg::SparseMatrix>
        H_;  ///< Matrix containing the tangent derivatives with respect to slave dofs
    std::shared_ptr<Core::LinAlg::SparseMatrix>
        T_;  ///< Matrix containing the tangent vectors of the slave nodes
    std::shared_ptr<Core::LinAlg::SparseMatrix>
        N_;  ///< Matrix containing the (weighted) gap derivatives
             ///< with respect to master and slave dofs
    std::shared_ptr<Core::LinAlg::Vector<double>> gap_;  ///< gap vector

    std::shared_ptr<CONTACT::Interface> interface_;  ///< interface
  };
}  // namespace Adapter

FOUR_C_NAMESPACE_CLOSE

#endif
