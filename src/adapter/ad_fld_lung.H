/*----------------------------------------------------------------------*/
/*! \file

\brief Fluid field adapter for fsi airway simulations with attached
parenchyma balloon

\level 2


*/
/*----------------------------------------------------------------------*/


#ifndef ADAPTER_FLUID_LUNG_H
#define ADAPTER_FLUID_LUNG_H

#include <set>

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>

#include "ad_fld_fluid_fsi.H"


// forward declarations
namespace DRT
{
  class Condition;
}

namespace LINALG
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
        Teuchos::RCP<LINALG::Solver> solver, Teuchos::RCP<Teuchos::ParameterList> params,
        Teuchos::RCP<IO::DiscretizationWriter> output, bool isale, bool dirichletcond);

    /// initialize algorithm
    virtual void Init();

    /// List of fluid-structure volume constraints
    void ListLungVolCons(std::set<int>& LungVolConIDs, int& MinLungVolConID);

    /// Initialize fluid part of lung volume constraint
    void InitializeVolCon(Teuchos::RCP<Epetra_Vector> initflowrate, const int offsetID);

    /// Evaluate fluid/ale part of lung volume constraint
    void EvaluateVolCon(Teuchos::RCP<LINALG::BlockSparseMatrixBase> FluidShapeDerivMatrix,
        Teuchos::RCP<LINALG::SparseMatrix> FluidConstrMatrix,
        Teuchos::RCP<LINALG::SparseMatrix> ConstrFluidMatrix,
        Teuchos::RCP<LINALG::BlockSparseMatrixBase> AleConstrMatrix,
        Teuchos::RCP<LINALG::BlockSparseMatrixBase> ConstrAleMatrix,
        Teuchos::RCP<Epetra_Vector> FluidRHS, Teuchos::RCP<Epetra_Vector> CurrFlowRates,
        Teuchos::RCP<Epetra_Vector> lagrMultVecRed, const int offsetID, const double dttheta);

    /// Write additional forces due to volume constraint
    void OutputForces(Teuchos::RCP<Epetra_Vector> Forces);

    /// Get map extractor for fsi <-> full map
    Teuchos::RCP<LINALG::MapExtractor> FSIInterface() { return fsiinterface_; }

    /// Get map extractor for asi, other <-> full inner map
    Teuchos::RCP<LINALG::MapExtractor> InnerSplit() { return innersplit_; }

   private:
    /// conditions, that define the lung volume constraints
    std::vector<DRT::Condition*> constrcond_;

    /// map extractor for fsi <-> full map
    /// this is needed since otherwise "OtherMap" contains only dofs
    /// which are not part of a condition. however, asi dofs are of
    /// course also "inner" dofs for the fluid field.
    Teuchos::RCP<LINALG::MapExtractor> fsiinterface_;

    /// map extractor for asi, other <-> full inner map
    Teuchos::RCP<LINALG::MapExtractor> innersplit_;

    /// map extractor for outflow fsi <-> full map
    Teuchos::RCP<LINALG::MapExtractor> outflowfsiinterface_;
  };
}  // namespace ADAPTER

#endif
