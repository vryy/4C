/*-----------------------------------------------------------------------*/
/*! \file
// Masterthesis of h.Willmann under supervision of Anh-Tu Vuong and Matthias Mayr
// Originates from contact_poro_lagrange_strategy

\brief Meshtying of porous media using Lagrange multipliers


\level 3
*/
/*-----------------------------------------------------------------------*/
#ifndef FOUR_C_CONTACT_MESHTYING_PORO_LAGRANGE_STRATEGY_HPP
#define FOUR_C_CONTACT_MESHTYING_PORO_LAGRANGE_STRATEGY_HPP

#include "4C_config.hpp"

#include "4C_contact_meshtying_lagrange_strategy.hpp"
#include "4C_utils_exceptions.hpp"

FOUR_C_NAMESPACE_OPEN

namespace CONTACT
{
  class PoroMtLagrangeStrategy : public MtLagrangeStrategy
  {
   public:
    /*!
    \brief Standard Constructor

    */
    PoroMtLagrangeStrategy(const Epetra_Map* dof_row_map, const Epetra_Map* NodeRowMap,
        Teuchos::ParameterList params, std::vector<Teuchos::RCP<Mortar::Interface>> interface,
        int dim, Teuchos::RCP<Epetra_Comm> comm, double alphaf, int maxdof);


    /*!
    /brief initial Poro Meshtying calculations
    */
    void InitializePoroMt(Teuchos::RCP<Core::LinAlg::SparseMatrix>& kteffoffdiag);

    /*!
    /brief modify off diagonal system matrix for structural displacement meshtying
     */
    void evaluate_meshtying_poro_off_diag(Teuchos::RCP<Core::LinAlg::SparseMatrix>& kteffoffdiag);

    /*!
    \brief Recovery method
    This method recovers the Langrangemultiplier correctly for the fluid-structure coupling
    matrix block. Complete Recovery only correct if the structural part is recovered elsewhere
    additionally. This is not needed for solving the problem in pure elast-poroelast meshtying
    cases. It is needed for postprocessing though.*/

    void recover_coupling_matrix_partof_lmp(Teuchos::RCP<Epetra_Vector> veli);

    Teuchos::RCP<Core::LinAlg::SparseMatrix> cs_;  // slave matrix block row (needed for LM)

    Teuchos::RCP<Epetra_Map> fvelrow_;  // fluid row map (needed for splitting)

  };  // class POROLagrangeStrategy
}  // namespace CONTACT
FOUR_C_NAMESPACE_CLOSE

#endif
