/*!----------------------------------------------------------------------
\file scatra_timint_meshtying_strategy_std_elch.cpp

\brief Standard solution strategy for electrochemistry problems (without meshtying)

<pre>
Maintainer: Rui Fang
            fang@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/
            089 - 289-15251
</pre>

*----------------------------------------------------------------------*/

#include "../drt_fluid/fluid_utils.H"

#include "../drt_scatra/scatra_utils_splitstrategy.H"

#include "scatra_timint_meshtying_strategy_std_elch.H"

/*----------------------------------------------------------------------*
 | constructor                                               fang 12/14 |
 *----------------------------------------------------------------------*/
SCATRA::MeshtyingStrategyStdElch::MeshtyingStrategyStdElch(
    SCATRA::ScaTraTimIntElch* elchtimint
    ) :
MeshtyingStrategyStd(elchtimint)
{
  return;
} // SCATRA::MeshtyingStrategyStdElch::MeshtyingStrategyStdElch


/*----------------------------------------------------------------------*
 | initialize system matrix for electrochemistry problems    fang 12/14 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<LINALG::SparseOperator> SCATRA::MeshtyingStrategyStdElch::InitSystemMatrix() const
{
  Teuchos::RCP<LINALG::SparseOperator> systemmatrix;

  if(DRT::INPUT::IntegralValue<int>(*(scatratimint_->ScatraParameterList()),"BLOCKPRECOND"))
  {
    // safety checks
    if(ElchTimInt()->EquPot() == INPAR::ELCH::equpot_undefined)
      dserror("Type of closing equation for electric potential not correctly set!");
    if(ElchTimInt()->EquPot() != INPAR::ELCH::equpot_enc)
      dserror("Special ELCH assemble strategy for block-matrix will not assemble A_11 block!");
    if(scatratimint_->NumScal() < 1)
      dserror("Number of transported scalars not correctly set!");

    // initial guess for non-zeros per row: 27 neighboring nodes for hex8
    // this is enough! A higher guess would require too much memory!
    // A more precise guess for every submatrix would read:
    // A_00: 27*1,  A_01: 27*1,  A_10: 27*numscal due to electroneutrality, A_11: EMPTY matrix !!!!!
    // usage of a split strategy that makes use of the ELCH-specific sparsity pattern
    LINALG::MapExtractor splitter;
    FLD::UTILS::SetupFluidSplit(*(scatratimint_->Discretization()),scatratimint_->NumScal(),splitter);
    systemmatrix = Teuchos::rcp(new LINALG::BlockSparseMatrix<SCATRA::SplitStrategy>(splitter,splitter,27,false,true));
    Teuchos::rcp_dynamic_cast<LINALG::BlockSparseMatrix<SCATRA::SplitStrategy> >(systemmatrix)->SetNumScal(scatratimint_->NumScal());
  }

  else
    systemmatrix = SCATRA::MeshtyingStrategyStd::InitSystemMatrix();

  return systemmatrix;
} // SCATRA::MeshtyingStrategyStdElch::InitSystemMatrix
