/*----------------------------------------------------------------------*/
/*! \file

\brief Fluid field adapter for fsi airway simulations with attached
parenchyma balloon

\level 2


*/
/*----------------------------------------------------------------------*/


#ifndef BACI_ADAPTER_FLD_LUNG_H
#define BACI_ADAPTER_FLD_LUNG_H

#include "baci_adapter_fld_fluid_fsi.H"

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>

#include <set>

BACI_NAMESPACE_OPEN


// forward declarations
namespace DRT
{
  class Condition;
}

namespace CORE::LINALG
{
  class MapExtractor;
}


namespace ADAPTER
{
  class FluidLung : public FluidFSI
  {
   public:
    /// Constructor
    FluidLung(Teuchos::RCP<Fluid> fluid, Teuchos::RCP<DRT::Discretization> dis,
        Teuchos::RCP<CORE::LINALG::Solver> solver, Teuchos::RCP<Teuchos::ParameterList> params,
        Teuchos::RCP<IO::DiscretizationWriter> output, bool isale, bool dirichletcond);

    /// initialize algorithm
    void Init() override;

    /// List of fluid-structure volume constraints
    void ListLungVolCons(std::set<int>& LungVolConIDs, int& MinLungVolConID);

    /// Initialize fluid part of lung volume constraint
    void InitializeVolCon(Teuchos::RCP<Epetra_Vector> initflowrate, const int offsetID);

    /// Evaluate fluid/ale part of lung volume constraint
    void EvaluateVolCon(Teuchos::RCP<CORE::LINALG::BlockSparseMatrixBase> FluidShapeDerivMatrix,
        Teuchos::RCP<CORE::LINALG::SparseMatrix> FluidConstrMatrix,
        Teuchos::RCP<CORE::LINALG::SparseMatrix> ConstrFluidMatrix,
        Teuchos::RCP<CORE::LINALG::BlockSparseMatrixBase> AleConstrMatrix,
        Teuchos::RCP<CORE::LINALG::BlockSparseMatrixBase> ConstrAleMatrix,
        Teuchos::RCP<Epetra_Vector> FluidRHS, Teuchos::RCP<Epetra_Vector> CurrFlowRates,
        Teuchos::RCP<Epetra_Vector> lagrMultVecRed, const int offsetID, const double dttheta);

    /// Write additional forces due to volume constraint
    void OutputForces(Teuchos::RCP<Epetra_Vector> Forces);

    /// Get map extractor for fsi <-> full map
    Teuchos::RCP<CORE::LINALG::MapExtractor> FSIInterface() { return fsiinterface_; }

    /// Get map extractor for asi, other <-> full inner map
    Teuchos::RCP<CORE::LINALG::MapExtractor> InnerSplit() { return innersplit_; }

   private:
    /// conditions, that define the lung volume constraints
    std::vector<DRT::Condition*> constrcond_;

    /// map extractor for fsi <-> full map
    /// this is needed since otherwise "OtherMap" contains only dofs
    /// which are not part of a condition. however, asi dofs are of
    /// course also "inner" dofs for the fluid field.
    Teuchos::RCP<CORE::LINALG::MapExtractor> fsiinterface_;

    /// map extractor for asi, other <-> full inner map
    Teuchos::RCP<CORE::LINALG::MapExtractor> innersplit_;

    /// map extractor for outflow fsi <-> full map
    Teuchos::RCP<CORE::LINALG::MapExtractor> outflowfsiinterface_;
  };
}  // namespace ADAPTER

BACI_NAMESPACE_CLOSE

#endif
