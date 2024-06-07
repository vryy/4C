/*----------------------------------------------------------------------*/
/*! \file

\brief Fluid field adapter for fsi airway simulations with attached
parenchyma balloon

\level 2


*/
/*----------------------------------------------------------------------*/


#ifndef FOUR_C_ADAPTER_FLD_LUNG_HPP
#define FOUR_C_ADAPTER_FLD_LUNG_HPP

#include "4C_config.hpp"

#include "4C_adapter_fld_fluid_fsi.hpp"
#include "4C_fem_condition.hpp"

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>

#include <set>

FOUR_C_NAMESPACE_OPEN


// forward declarations

namespace Core::LinAlg
{
  class MapExtractor;
}


namespace Adapter
{
  class FluidLung : public FluidFSI
  {
   public:
    /// Constructor
    FluidLung(Teuchos::RCP<Fluid> fluid, Teuchos::RCP<Discret::Discretization> dis,
        Teuchos::RCP<Core::LinAlg::Solver> solver, Teuchos::RCP<Teuchos::ParameterList> params,
        Teuchos::RCP<Core::IO::DiscretizationWriter> output, bool isale, bool dirichletcond);

    /// initialize algorithm
    void Init() override;

    /// List of fluid-structure volume constraints
    void ListLungVolCons(std::set<int>& LungVolConIDs, int& MinLungVolConID);

    /// Initialize fluid part of lung volume constraint
    void InitializeVolCon(Teuchos::RCP<Epetra_Vector> initflowrate, const int offsetID);

    /// Evaluate fluid/ale part of lung volume constraint
    void EvaluateVolCon(Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase> FluidShapeDerivMatrix,
        Teuchos::RCP<Core::LinAlg::SparseMatrix> FluidConstrMatrix,
        Teuchos::RCP<Core::LinAlg::SparseMatrix> ConstrFluidMatrix,
        Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase> AleConstrMatrix,
        Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase> ConstrAleMatrix,
        Teuchos::RCP<Epetra_Vector> FluidRHS, Teuchos::RCP<Epetra_Vector> CurrFlowRates,
        Teuchos::RCP<Epetra_Vector> lagrMultVecRed, const int offsetID, const double dttheta);

    /// Write additional forces due to volume constraint
    void OutputForces(Teuchos::RCP<Epetra_Vector> Forces);

    /// Get map extractor for fsi <-> full map
    Teuchos::RCP<Core::LinAlg::MapExtractor> FSIInterface() { return fsiinterface_; }

    /// Get map extractor for asi, other <-> full inner map
    Teuchos::RCP<Core::LinAlg::MapExtractor> InnerSplit() { return innersplit_; }

   private:
    /// conditions, that define the lung volume constraints
    std::vector<Core::Conditions::Condition*> constrcond_;

    /// map extractor for fsi <-> full map
    /// this is needed since otherwise "OtherMap" contains only dofs
    /// which are not part of a condition. however, asi dofs are of
    /// course also "inner" dofs for the fluid field.
    Teuchos::RCP<Core::LinAlg::MapExtractor> fsiinterface_;

    /// map extractor for asi, other <-> full inner map
    Teuchos::RCP<Core::LinAlg::MapExtractor> innersplit_;

    /// map extractor for outflow fsi <-> full map
    Teuchos::RCP<Core::LinAlg::MapExtractor> outflowfsiinterface_;
  };
}  // namespace Adapter

FOUR_C_NAMESPACE_CLOSE

#endif
