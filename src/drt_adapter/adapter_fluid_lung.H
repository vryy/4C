/*----------------------------------------------------------------------*/
/*!
\file adapter_fluid_lung.H

\brief Fluid field adapter for fsi airway simulations with attached
parenchyma balloon

<pre>
Maintainer: Lena Wiechert
            wiechert@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15303
</pre>
*/
/*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#ifndef ADAPTER_FLUID_LUNG_H
#define ADAPTER_FLUID_LUNG_H

#include <Teuchos_RefCountPtr.hpp>
#include <Teuchos_ParameterList.hpp>

#include "adapter_fluid_wrapper.H"


namespace ADAPTER {

  class FluidLung : public FluidWrapper
  {
  public:

    /// Constructor
    FluidLung(Teuchos::RCP<Fluid> fluid);

    /// List of fluid-structure volume constraints
    void ListLungVolCons(std::set<int>& LungVolConIDs,
                         int& MinLungVolConID);

    /// Initialize fluid part of lung volume constraint
    void InitializeVolCon(Teuchos::RCP<Epetra_Vector> initflowrate,
                          const int offsetID);

    /// Evaluate fluid/ale part of lung volume constraint
    void EvaluateVolCon(Teuchos::RCP<LINALG::BlockSparseMatrixBase> FluidShapeDerivMatrix,
                        Teuchos::RCP<LINALG::SparseMatrix> FluidConstrMatrix,
                        Teuchos::RCP<LINALG::SparseMatrix> ConstrFluidMatrix,
                        Teuchos::RCP<LINALG::BlockSparseMatrixBase> AleConstrMatrix,
                        Teuchos::RCP<LINALG::BlockSparseMatrixBase> ConstrAleMatrix,
                        Teuchos::RCP<Epetra_Vector> FluidRHS,
                        Teuchos::RCP<Epetra_Vector> CurrFlowRates,
                        Teuchos::RCP<Epetra_Vector> lagrMultVecRed,
                        const int offsetID,
                        const double dttheta);

  private:

    /// conditions, that define the lung volume constraints
    std::vector<DRT::Condition*> constrcond_;

  };

}

#endif
#endif
