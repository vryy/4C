/*-----------------------------------------------------------------------*/
/*! \file
// Masterthesis of h.Willmann under supervision of Anh-Tu Vuong and Matthias Mayr
// Originates from contact_poro_lagrange_strategy

\brief Meshtying of porous media using Lagrange multipliers

\maintainer Matthias Mayr

\level 3
*/
/*-----------------------------------------------------------------------*/
#ifndef MESHTYING_PORO_LAGRANGE_STRATEGY_H
#define MESHTYING_PORO_LAGRANGE_STRATEGY_H

#include "meshtying_lagrange_strategy.H"
#include "../drt_lib/drt_dserror.H"

namespace CONTACT
{
  class PoroMtLagrangeStrategy : public MtLagrangeStrategy
  {
   public:
    /*!
    \brief Standard Constructor

    */
    PoroMtLagrangeStrategy(const Epetra_Map* DofRowMap, const Epetra_Map* NodeRowMap,
        Teuchos::ParameterList params, std::vector<Teuchos::RCP<MORTAR::MortarInterface>> interface,
        int dim, Teuchos::RCP<Epetra_Comm> comm, double alphaf, int maxdof);

    /*!
    \brief Destructor

    */
    virtual ~PoroMtLagrangeStrategy(){};

    /*!
    /brief initial Poro Meshtying calculations
    */
    void InitializePoroMt(Teuchos::RCP<LINALG::SparseMatrix>& kteffoffdiag);

    /*!
    /brief modify off diagonal system matrix for structural displacement meshtying
     */
    void EvaluateMeshtyingPoroOffDiag(Teuchos::RCP<LINALG::SparseMatrix>& kteffoffdiag);

    /*!
    \brief Recovery method
    This method recovers the Langrangemultiplier correctly for the fluid-structure coupling
    matrix block. Complete Recovery only correct if the structural part is recovered elsewhere
    additionally. This is not needed for solving the problem in pure elast-poroelast meshtying
    cases. It is needed for postprocessing though.*/

    void RecoverCouplingMatrixPartofLMP(Teuchos::RCP<Epetra_Vector> veli);

    Teuchos::RCP<LINALG::SparseMatrix> cs_;  // slave matrix block row (needed for LM)

    Teuchos::RCP<Epetra_Map> fvelrow_;  // fluid row map (needed for splitting)

  };  // class POROLagrangeStrategy
}  // namespace CONTACT
#endif  // #ifndef CONTACT_PORO_LAGRANGE_STRATEGY_H
