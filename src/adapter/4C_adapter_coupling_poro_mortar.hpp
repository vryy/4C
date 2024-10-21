// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_ADAPTER_COUPLING_PORO_MORTAR_HPP
#define FOUR_C_ADAPTER_COUPLING_PORO_MORTAR_HPP

/*---------------------------------------------------------------------*
 | headers                                                  ager 10/15 |
 *---------------------------------------------------------------------*/
#include "4C_config.hpp"

#include "4C_adapter_coupling_nonlin_mortar.hpp"
#include "4C_contact_lagrange_strategy_poro.hpp"
#include "4C_coupling_adapter.hpp"
#include "4C_linalg_vector.hpp"

#include <Epetra_Comm.h>
#include <Epetra_Map.h>
#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN
/*---------------------------------------------------------------------*
 | forward declarations                                     ager 10/15 |
 *---------------------------------------------------------------------*/
namespace Core::FE
{
  class Discretization;
}  // namespace Core::FE

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
  class CouplingPoroMortar : public CouplingNonLinMortar
  {
   public:
    /*!
    \brief Empty constructor

    */
    CouplingPoroMortar(int spatial_dimension, Teuchos::ParameterList mortar_coupling_params,
        Teuchos::ParameterList contact_dynamic_params,
        Core::FE::ShapeFunctionType shape_function_type);


    virtual void evaluate_poro_mt(Teuchos::RCP<Core::LinAlg::Vector<double>> fvel,
        Teuchos::RCP<Core::LinAlg::Vector<double>> svel,
        Teuchos::RCP<Core::LinAlg::Vector<double>> fpres,
        Teuchos::RCP<Core::LinAlg::Vector<double>> sdisp,
        const Teuchos::RCP<Core::FE::Discretization> sdis,
        Teuchos::RCP<Core::LinAlg::SparseMatrix>& f, Teuchos::RCP<Core::LinAlg::SparseMatrix>& k_fs,
        Teuchos::RCP<Core::LinAlg::Vector<double>>& frhs, Coupling::Adapter::Coupling& coupfs,
        Teuchos::RCP<const Epetra_Map> fdofrowmap);

    void update_poro_mt();

    void recover_fluid_lm_poro_mt(Teuchos::RCP<Core::LinAlg::Vector<double>> disi,
        Teuchos::RCP<Core::LinAlg::Vector<double>> veli);  // h.Willmann

    // return the used poro lagrange strategy
    Teuchos::RCP<CONTACT::LagrangeStrategyPoro> get_poro_strategy()
    {
      if (porolagstrategy_ == Teuchos::null) FOUR_C_THROW("GetPoroStrategy(): No strategy set!");
      return porolagstrategy_;
    };

   protected:
    /*!
    \brief Read Mortar Condition

    */
    void read_mortar_condition(Teuchos::RCP<Core::FE::Discretization> masterdis,
        Teuchos::RCP<Core::FE::Discretization> slavedis, std::vector<int> coupleddof,
        const std::string& couplingcond, Teuchos::ParameterList& input,
        std::map<int, Core::Nodes::Node*>& mastergnodes,
        std::map<int, Core::Nodes::Node*>& slavegnodes,
        std::map<int, Teuchos::RCP<Core::Elements::Element>>& masterelements,
        std::map<int, Teuchos::RCP<Core::Elements::Element>>& slaveelements) override;

    /*!
    \brief Add Mortar Elments

    */
    void add_mortar_elements(Teuchos::RCP<Core::FE::Discretization> masterdis,
        Teuchos::RCP<Core::FE::Discretization> slavedis, Teuchos::ParameterList& input,
        std::map<int, Teuchos::RCP<Core::Elements::Element>>& masterelements,
        std::map<int, Teuchos::RCP<Core::Elements::Element>>& slaveelements,
        Teuchos::RCP<CONTACT::Interface>& interface, int numcoupleddof) override;

    /*!
    \brief complete interface, store as internal variable
           store maps as internal variable and do parallel redist.

    */
    void complete_interface(Teuchos::RCP<Core::FE::Discretization> masterdis,
        Teuchos::RCP<CONTACT::Interface>& interface) override;

    /*!
    \brief create strategy object if required

    */
    void create_strategy(Teuchos::RCP<Core::FE::Discretization> masterdis,
        Teuchos::RCP<Core::FE::Discretization> slavedis, Teuchos::ParameterList& input,
        int numcoupleddof) override;

   private:
    // poro lagrange strategy
    Teuchos::RCP<CONTACT::LagrangeStrategyPoro> porolagstrategy_;

    // firstinit
    bool firstinit_;

    int slavetype_;   // 1 poro, 0 struct, -1 default
    int mastertype_;  // 1 poro, 0 struct, -1 default
  };
}  // namespace Adapter

FOUR_C_NAMESPACE_CLOSE

#endif
