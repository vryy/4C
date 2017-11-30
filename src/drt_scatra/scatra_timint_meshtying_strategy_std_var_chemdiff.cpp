/*!----------------------------------------------------------------------
\file scatra_timint_meshtying_strategy_std_var_chemdiff.cpp

\brief Standard solution strategy for variational chemical diffusion
     problems (without meshtying)

\level 3

<pre>
\maintainer Jorge De Anda Salazar
            deanda@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/
            089-289-15251
</pre>

*----------------------------------------------------------------------*/
#include "scatra_timint_meshtying_strategy_std_var_chemdiff.H"
#include "../drt_fluid/fluid_utils.H"
#include "../drt_scatra/scatra_utils_splitstrategy.H"
#include "../linalg/linalg_solver.H"

/*----------------------------------------------------------------------*
 | constructor                                              deanda 10/17 |
 *----------------------------------------------------------------------*/
SCATRA::MeshtyingStrategyStdVar::MeshtyingStrategyStdVar(
    SCATRA::TimIntVariational* vartimint
    ) :
    MeshtyingStrategyStd(vartimint)
{
  return;
} // SCATRA::MeshtyingStrategyStdVar::MeshtyingStrategyStdVar

/*----------------------------------------------------------------------*
 | initialize system matrix for variational problems       deanda 10/17 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<LINALG::SparseOperator> SCATRA::MeshtyingStrategyStdVar::InitSystemMatrix() const
{
  Teuchos::RCP<LINALG::SparseOperator> systemmatrix;

  if(DRT::INPUT::IntegralValue<int>((VarTimInt()->VarParameterList())->sublist("VARIATIONAL"),"BLOCKPRECOND"))
  {
    systemmatrix = Teuchos::rcp(new LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>(*scatratimint_->Splitter(),*scatratimint_->Splitter(),27,false,true));
  }

  else
    systemmatrix = SCATRA::MeshtyingStrategyStd::InitSystemMatrix();

  return systemmatrix;
} // SCATRA::MeshtyingStrategyStdVar::InitSystemMatrix

/*----------------------------------------------------------------------*
 | setup meshtying objects                                  deanda 10/17 |
 *----------------------------------------------------------------------*/
void SCATRA::MeshtyingStrategyStdVar::SetupMeshtying()
{

  // feed AMGnxn block preconditioner with null space information for each block of global block system matrix
  BuildBlockNullSpaces();

  return;
}

/*-------------------------------------------------------------------------------*
 | build null spaces associated with blocks of global system matrix  deanda 10/17 |
 *-------------------------------------------------------------------------------*/
void SCATRA::MeshtyingStrategyStdVar::BuildBlockNullSpaces() const
{
  // loop over blocks of global system matrix
  for(int iblock=0; iblock< scatratimint_->Splitter()->NumMaps(); ++iblock)
  {
    // store number of current block as string, starting from 1
    std::ostringstream iblockstr;

    iblockstr << iblock+1;

    // equip smoother for current matrix block with empty parameter sublists to trigger null space computation
    Teuchos::ParameterList& blocksmootherparams = scatratimint_->Solver()->Params().sublist("Inverse"+iblockstr.str());

    blocksmootherparams.sublist("Aztec Parameters");
    blocksmootherparams.sublist("MueLu Parameters");


    Teuchos::ParameterList& mllist = blocksmootherparams.sublist("MueLu Parameters",true);
    mllist.set("PDE equations",1);
    mllist.set("null space: dimension",1);
    mllist.set("null space: type","pre-computed");
    mllist.set("null space: add default vectors",false);
    const Teuchos::RCP<std::vector<double> > ns = Teuchos::rcp(new std::vector<double>(scatratimint_->Discretization()->DofRowMap()->NumMyElements(),1.));
    mllist.set<Teuchos::RCP<std::vector<double> > >("nullspace",ns);
    mllist.set("null space: vectors",&((*ns)[0]));
    mllist.set<bool>("ML validate parameter list",false);


    // equip smoother for current matrix block with null space associated with all degrees of freedom on discretization
    scatratimint_->Discretization()->ComputeNullSpaceIfNecessary(blocksmootherparams);

    // reduce full null space to match degrees of freedom associated with current matrix block
    LINALG::Solver::FixMLNullspace("Block "+iblockstr.str(),*scatratimint_->Discretization()->DofRowMap(),*scatratimint_->Splitter()->Map(iblock),blocksmootherparams);
  }

  return;
} // SCATRA::MeshtyingStrategyStdVar::BuildBlockNullSpaces
