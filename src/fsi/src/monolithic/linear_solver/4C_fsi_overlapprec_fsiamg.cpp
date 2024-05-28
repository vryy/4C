/*----------------------------------------------------------------------*/
/*! \file

\brief Special version of block matrix that includes the FSI block preconditioner


\level 1
*/
/*----------------------------------------------------------------------*/

#include "4C_fsi_overlapprec_fsiamg.hpp"

#include "4C_adapter_fld_fluid.hpp"
#include "4C_adapter_str_fsiwrapper.hpp"
#include "4C_comm_utils.hpp"
#include "4C_fsi_overlapprec_hybrid.hpp"
#include "4C_global_data.hpp"
#include "4C_linear_solver_method_linalg.hpp"
#include "4C_linear_solver_preconditioner_linalg.hpp"

#include <EpetraExt_SolverMap_CrsMatrix.h>
#include <ml_MultiLevelPreconditioner.h>
#include <MLAPI_CompObject.h>
#include <MLAPI_Expressions.h>
#include <MLAPI_LoadBalanceInverseOperator.h>
#include <MLAPI_LoadBalanceOperator.h>
#include <MLAPI_MultiVector.h>
#include <MLAPI_Operator_Utils.h>
#include <MLAPI_Workspace.h>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
FSI::OverlappingBlockMatrixFSIAMG::OverlappingBlockMatrixFSIAMG(
    const CORE::LINALG::MultiMapExtractor& maps, ADAPTER::FSIStructureWrapper& structure,
    ADAPTER::Fluid& fluid, ADAPTER::AleFsiWrapper& ale, bool structuresplit, int symmetric,
    std::vector<std::string>& blocksmoother, std::vector<double>& schuromega,
    std::vector<double>& omega, std::vector<int>& iterations, std::vector<double>& somega,
    std::vector<int>& siterations, std::vector<double>& fomega, std::vector<int>& fiterations,
    std::vector<double>& aomega, std::vector<int>& aiterations, int analyze,
    INPAR::FSI::LinearBlockSolver strategy, INPAR::FSI::Verbosity verbosity,
    OverlappingBlockMatrixHybridSchwarz* hybridPrec)
    : OverlappingBlockMatrix(Teuchos::null, maps, structure, fluid, ale, structuresplit, symmetric,
          omega[0], iterations[0], somega[0],
          siterations[0] - 1,  // base class counts iterations starting from 0
          fomega[0], fiterations[0] - 1, aomega[0], aiterations[0] - 1),
      sisml_(false),
      fisml_(false),
      aisml_(false),
      srun_(0),
      frun_(0),
      arun_(0),
      minnlevel_(0),
      analyze_(analyze),
      strategy_(strategy),
      blocksmoother_(blocksmoother),
      schuromega_(schuromega),
      pcomega_(omega),
      pciter_(iterations),
      somega_(somega),
      siterations_(siterations),
      fomega_(fomega),
      fiterations_(fiterations),
      aomega_(aomega),
      aiterations_(aiterations),
      verbosity_(verbosity),
      hybridPrec_(hybridPrec)
{
  if (strategy_ != INPAR::FSI::PreconditionedKrylov && strategy_ != INPAR::FSI::LinalgSolver)
    FOUR_C_THROW("Type of LINEARBLOCKSOLVER parameter not recognized by this class");
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void FSI::OverlappingBlockMatrixFSIAMG::SetupPreconditioner()
{
  TEUCHOS_FUNC_TIME_MONITOR("FSI::OverlappingBlockMatrixFSIAMG::SetupPreconditioner");

#ifdef BLOCKMATRIXMERGE
  FOUR_C_THROW("class OverlappingBlockMatrixFSIAMG does not support #define BLOCKMATRIXMERGE");
#endif

  if (GLOBAL::Problem::Instance()->GetCommunicators()->NumGroups() != 1)
    FOUR_C_THROW(
        "No nested parallelism for AMG FSI. See comments in "
        "FSI::OverlappingBlockMatrixFSIAMG::SetupPreconditioner()!");
  // Attention: No nested parallelism for AMG FSI due to MLAPI incompatibility
  // MLAPI::Space::Reshape constructs an ML_RowMatrix object using a hard coded MPI_COMM_WORLD in
  // MLAPI_Operator.h Fixing this needs major changes in Trilinos/MLAPI which is not desirable.

  // MLAPI::Init() without arguments uses internally MPI_COMM_WOLRD
  MLAPI::Init();
  const int myrank = (hybridPrec_ == nullptr) ? Matrix(0, 0).Comm().MyPID()
                                              : hybridPrec_->Matrix(0, 0).Comm().MyPID();

  const CORE::LINALG::SparseMatrix& structInnerOp =
      (hybridPrec_ == nullptr) ? Matrix(0, 0) : hybridPrec_->Matrix(0, 0);
  const CORE::LINALG::SparseMatrix& fluidInnerOp =
      (hybridPrec_ == nullptr) ? Matrix(1, 1) : hybridPrec_->Matrix(1, 1);
  const CORE::LINALG::SparseMatrix& aleInnerOp =
      (hybridPrec_ == nullptr) ? Matrix(2, 2) : hybridPrec_->Matrix(2, 2);

  Teuchos::RCP<CORE::LINALG::MapExtractor> fsidofmapex = Teuchos::null;
  Teuchos::RCP<Epetra_Map> irownodes = Teuchos::null;

  // build AMG hierarchies
  structuresolver_->Setup(structInnerOp.EpetraMatrix());
  fluidsolver_->Setup(fluidInnerOp.EpetraMatrix(), fsidofmapex, fluid_.discretization(), irownodes,
      structuresplit_);
  if (constalesolver_ == Teuchos::null) alesolver_->Setup(aleInnerOp.EpetraMatrix());

  // get the ml_MultiLevelPreconditioner class from within struct/fluid/ale solver
  Teuchos::RCP<Epetra_Operator> sprec = structuresolver_->EpetraOperator();
  Teuchos::RCP<Epetra_Operator> fprec = fluidsolver_->EpetraOperator();
  Teuchos::RCP<Epetra_Operator> aprec = alesolver_->EpetraOperator();

  // get ML preconditioner class
  ML_Epetra::MultiLevelPreconditioner* smlclass =
      dynamic_cast<ML_Epetra::MultiLevelPreconditioner*>(sprec.get());
  ML_Epetra::MultiLevelPreconditioner* fmlclass =
      dynamic_cast<ML_Epetra::MultiLevelPreconditioner*>(fprec.get());
  ML_Epetra::MultiLevelPreconditioner* amlclass =
      dynamic_cast<ML_Epetra::MultiLevelPreconditioner*>(aprec.get());

  if (smlclass) sisml_ = true;
  if (fmlclass) fisml_ = true;
  if (amlclass) aisml_ = true;

  // get copy of ML parameter list
  if (sis_amg()) sparams_ = structuresolver_->Params().sublist("ML Parameters");
  if (fis_amg()) fparams_ = fluidsolver_->Params().sublist("ML Parameters");
  if (ais_amg()) aparams_ = alesolver_->Params().sublist("ML Parameters");
  // cout << "====================Structure Params\n" << sparams_;
  // cout << "========================Fluid Params\n" << fparams_;
  // cout << "========================Ale Params\n" << aparams_;
  // find for which field using nonsymmetric AMG
  bool SisPetrovGalerkin = sparams_.get("energy minimization: enable", false);
  bool FisPetrovGalerkin = fparams_.get("energy minimization: enable", false);
  bool AisPetrovGalerkin = aparams_.get("energy minimization: enable", false);

  // get ML handle
  ML* sml = nullptr;
  ML* fml = nullptr;
  ML* aml = nullptr;
  if (sis_amg()) sml = const_cast<ML*>(smlclass->GetML());
  if (fis_amg()) fml = const_cast<ML*>(fmlclass->GetML());
  if (ais_amg()) aml = const_cast<ML*>(amlclass->GetML());
  if ((!sml && sis_amg()) || (!fml && fis_amg()) || (!aml && ais_amg()))
    FOUR_C_THROW("Not using ML for Fluid, Structure or Ale");

  // number of grids in structure, fluid and ale
  if (sis_amg())
    snlevel_ = sml->ML_num_actual_levels;
  else
    snlevel_ = 1;
  if (fis_amg())
    fnlevel_ = fml->ML_num_actual_levels;
  else
    fnlevel_ = 1;
  if (ais_amg())
    anlevel_ = aml->ML_num_actual_levels;
  else
    anlevel_ = 1;

  minnlevel_ = 1;
  maxnlevel_ = 1;

  if (!myrank && verbosity_ == INPAR::FSI::verbosity_full)
  {
    printf("       -----------------------------------------------------------------------\n");
    switch (strategy_)
    {
      case INPAR::FSI::PreconditionedKrylov:
      {
        printf(
            "       Setting up BGS(AMG): snlevel %d fnlevel %d anlevel %d minnlevel %d maxnlevel "
            "%d\n",
            snlevel_, fnlevel_, anlevel_, minnlevel_, maxnlevel_);
        fflush(stdout);
        break;
      }
      default:
        break;
    }
  }

  // check whether we have enough iteration and damping factors
  if ((int)pciter_.size() < 3 || (int)pcomega_.size() < 3)
    FOUR_C_THROW("You need at least 3 values of PCITER and PCOMEGA in input file");
  if ((int)siterations_.size() < snlevel_ || (int)somega_.size() < snlevel_)
    FOUR_C_THROW(
        "You need at least %d values of STRUCTPCITER and STRUCTPCOMEGA in input file", snlevel_);
  if ((int)fiterations_.size() < fnlevel_ || (int)fomega_.size() < fnlevel_)
    FOUR_C_THROW(
        "You need at least %d values of FLUIDPCITER and FLUIDPCOMEGA in input file", fnlevel_);
  if ((int)aiterations_.size() < anlevel_ || (int)aomega_.size() < anlevel_)
    FOUR_C_THROW("You need at least %d values of ALEPCITER and ALEPCOMEGA in input file", anlevel_);
  if ((int)blocksmoother_.size() < 1)
    FOUR_C_THROW("You need at least 1 value of BLOCKSMOOTHER in input file");
  if (blocksmoother_[0] == "Schur") FOUR_C_THROW("Schur(AMG) not implemented");

  Ass_.resize(std::max(maxnlevel_, snlevel_));
  Pss_.resize(std::max(maxnlevel_, snlevel_) - 1);
  Rss_.resize(std::max(maxnlevel_, snlevel_) - 1);
  Sss_.resize(std::max(maxnlevel_, snlevel_));

  Aff_.resize(std::max(maxnlevel_, fnlevel_));
  Pff_.resize(std::max(maxnlevel_, fnlevel_) - 1);
  Rff_.resize(std::max(maxnlevel_, fnlevel_) - 1);
  Sff_.resize(std::max(maxnlevel_, fnlevel_));
  Schurff_.resize(std::max(maxnlevel_, fnlevel_));

  Aaa_.resize(std::max(maxnlevel_, anlevel_));
  Paa_.resize(std::max(maxnlevel_, anlevel_) - 1);
  Raa_.resize(std::max(maxnlevel_, anlevel_) - 1);
  Saa_.resize(std::max(maxnlevel_, anlevel_));

  ASF_.resize(maxnlevel_);
  AFS_.resize(maxnlevel_);
  AFA_.resize(maxnlevel_);
  AAF_.resize(maxnlevel_);

  //---------------------------------------------------------- timing
  Teuchos::Time etime("", true);
  //------------------------------------------------------- Structure
  {
    // fine space matching Epetra objects
    MLAPI::Space finespace(structInnerOp.RowMap());
    if (!myrank && verbosity_ == INPAR::FSI::verbosity_full)
      printf("       Structure: NumGlobalElements fine level %d\n",
          structInnerOp.RowMap().NumGlobalElements());

    // extract transfer operator P,R from ML
    MLAPI::Space fspace;
    MLAPI::Space cspace;
    MLAPI::Operator P;
    MLAPI::Operator R;
    for (int i = 1; i < snlevel_; ++i)
    {
      ML_Operator* Pml = &(sml->Pmat[i]);
      if (i == 1)
        fspace = finespace;
      else
        fspace.Reshape(-1, Pml->outvec_leng);
      cspace.Reshape(-1, Pml->invec_leng);
      P.Reshape(cspace, fspace, Pml, false);
      if (SisPetrovGalerkin)
      {
        ML_Operator* Rml = &(sml->Rmat[i - 1]);
        R.Reshape(fspace, cspace, Rml, false);
      }
      else
        R = MLAPI::GetTranspose(P);
      Pss_[i - 1] = P;
      Rss_[i - 1] = R;
    }
    // extract matrix A from ML
    MLAPI::Space space;
    for (int i = 0; i < snlevel_; ++i)
    {
      ML_Operator* Aml = nullptr;
      if (sis_amg()) Aml = &(sml->Amat[i]);
      if (i == 0)
        space = finespace;
      else
        space.Reshape(-1, Aml->invec_leng);
      if (sis_amg())
      {
        MLAPI::Operator A(space, space, Aml, false);
        Ass_[i] = A;
      }
      else
      {
        MLAPI::Operator A(space, space, structInnerOp.EpetraMatrix().get(), false);
        Ass_[i] = A;
      }
    }
    for (int i = snlevel_; i < maxnlevel_; ++i)
    {
      Ass_[i] = Ass_[i - 1];
      fspace = Ass_[i - 1].GetRangeSpace();
      Pss_[i - 1] = GetIdentity(fspace, fspace);
      Rss_[i - 1] = Pss_[i - 1];
    }
  }

  //------------------------------------------------------- Fluid
  {
    // fine space matching Epetra objects
    MLAPI::Space finespace(fluidInnerOp.RowMap());
    if (!myrank && verbosity_ == INPAR::FSI::verbosity_full)
      printf("       Fluid    : NumGlobalElements fine level %d\n",
          fluidInnerOp.RowMap().NumGlobalElements());

    // extract transfer operator P,R from ML
    MLAPI::Space fspace;
    MLAPI::Space cspace;
    MLAPI::Operator P;
    MLAPI::Operator R;
    for (int i = 1; i < fnlevel_; ++i)
    {
      ML_Operator* Pml = &(fml->Pmat[i]);
      if (i == 1)
        fspace = finespace;
      else
        fspace.Reshape(-1, Pml->outvec_leng);
      cspace.Reshape(-1, Pml->invec_leng);
      P.Reshape(cspace, fspace, Pml, false);
      if (FisPetrovGalerkin)
      {
        ML_Operator* Rml = &(fml->Rmat[i - 1]);
        R.Reshape(fspace, cspace, Rml, false);
      }
      else
      {
        R = MLAPI::GetTranspose(P);
      }
      Pff_[i - 1] = P;
      Rff_[i - 1] = R;
    }
    // extract matrix A from ML
    MLAPI::Space space;
    for (int i = 0; i < fnlevel_; ++i)
    {
      ML_Operator* Aml = nullptr;
      if (fis_amg()) Aml = &(fml->Amat[i]);
      if (i == 0)
        space = finespace;
      else
        space.Reshape(-1, Aml->invec_leng);
      if (fis_amg())
      {
        MLAPI::Operator A(space, space, Aml, false);
        Aff_[i] = A;
      }
      else
      {
        MLAPI::Operator A(space, space, fluidInnerOp.EpetraMatrix().get(), false);
        Aff_[i] = A;
      }
    }
    for (int i = fnlevel_; i < maxnlevel_; ++i)
    {
      Aff_[i] = Aff_[i - 1];
      fspace = Aff_[i - 1].GetRangeSpace();
      Pff_[i - 1] = GetIdentity(fspace, fspace);
      Rff_[i - 1] = Pff_[i - 1];
    }
  }

  //------------------------------------------------------- Ale
  {
    // fine space matching Epetra objects
    MLAPI::Space finespace(aleInnerOp.RowMap());
    if (!myrank && verbosity_ == INPAR::FSI::verbosity_full)
      printf("       Ale      : NumGlobalElements fine level %d\n",
          aleInnerOp.RowMap().NumGlobalElements());

    // extract transfer operator P,R from ML
    MLAPI::Space fspace;
    MLAPI::Space cspace;
    MLAPI::Operator P;
    MLAPI::Operator R;
    for (int i = 1; i < anlevel_; ++i)
    {
      ML_Operator* Pml = &(aml->Pmat[i]);
      if (i == 1)
        fspace = finespace;
      else
        fspace.Reshape(-1, Pml->outvec_leng);
      cspace.Reshape(-1, Pml->invec_leng);
      P.Reshape(cspace, fspace, Pml, false);
      if (AisPetrovGalerkin)
      {
        ML_Operator* Rml = &(aml->Rmat[i - 1]);
        R.Reshape(fspace, cspace, Rml, false);
      }
      else
        R = MLAPI::GetTranspose(P);

      Paa_[i - 1] = P;
      Raa_[i - 1] = R;
    }
    // extract matrix A from ML
    MLAPI::Space space;
    for (int i = 0; i < anlevel_; ++i)
    {
      ML_Operator* Aml = nullptr;
      if (ais_amg()) Aml = &(aml->Amat[i]);
      if (i == 0)
        space = finespace;
      else
        space.Reshape(-1, Aml->invec_leng);
      if (ais_amg())
      {
        MLAPI::Operator A(space, space, Aml, false);
        Aaa_[i] = A;
      }
      else
      {
        MLAPI::Operator A(space, space, aleInnerOp.EpetraMatrix().get(), false);
        Aaa_[i] = A;
      }
    }
    for (int i = anlevel_; i < maxnlevel_; ++i)
    {
      Aaa_[i] = Aaa_[i - 1];
      fspace = Aaa_[i - 1].GetRangeSpace();
      Paa_[i - 1] = GetIdentity(fspace, fspace);
      Raa_[i - 1] = Paa_[i - 1];
    }
  }

  //-----------------------------------------------------------------
  // wrap the off-diagonal matrix blocks into MLAPI operators
  {
    const CORE::LINALG::SparseMatrix& Matrix01 =
        (hybridPrec_ == nullptr) ? Matrix(0, 1) : hybridPrec_->Matrix(0, 1);
    MLAPI::Space dspace(Matrix01.EpetraMatrix()->DomainMap());
    MLAPI::Space rspace(Matrix01.EpetraMatrix()->RangeMap());
    Asf_.Reshape(dspace, rspace, Matrix01.EpetraMatrix().get(), false);
    ASF_[0] = Asf_;
  }
  {
    const CORE::LINALG::SparseMatrix& Matrix10 =
        (hybridPrec_ == nullptr) ? Matrix(1, 0) : hybridPrec_->Matrix(1, 0);
    MLAPI::Space dspace(Matrix10.EpetraMatrix()->DomainMap());
    MLAPI::Space rspace(Matrix10.EpetraMatrix()->RangeMap());
    Afs_.Reshape(dspace, rspace, Matrix10.EpetraMatrix().get(), false);
    AFS_[0] = Afs_;
  }
  {
    const CORE::LINALG::SparseMatrix& Matrix12 =
        (hybridPrec_ == nullptr) ? Matrix(1, 2) : hybridPrec_->Matrix(1, 2);
    MLAPI::Space dspace(Matrix12.EpetraMatrix()->DomainMap());
    MLAPI::Space rspace(Matrix12.EpetraMatrix()->RangeMap());
    Afa_.Reshape(dspace, rspace, Matrix12.EpetraMatrix().get(), false);
    AFA_[0] = Afa_;
  }

  if (structuresplit_)
  {
    const CORE::LINALG::SparseMatrix& Matrix21 =
        (hybridPrec_ == nullptr) ? Matrix(2, 1) : hybridPrec_->Matrix(2, 1);
    MLAPI::Space dspace(Matrix21.EpetraMatrix()->DomainMap());
    MLAPI::Space rspace(Matrix21.EpetraMatrix()->RangeMap());
    Aaf_.Reshape(dspace, rspace, Matrix21.EpetraMatrix().get(), false);
    AAF_[0] = Aaf_;
  }
  else
  {
    const CORE::LINALG::SparseMatrix& Matrix20 =
        (hybridPrec_ == nullptr) ? Matrix(2, 0) : hybridPrec_->Matrix(2, 0);
    MLAPI::Space dspace(Matrix20.EpetraMatrix()->DomainMap());
    MLAPI::Space rspace(Matrix20.EpetraMatrix()->RangeMap());
    Aaf_.Reshape(dspace, rspace, Matrix20.EpetraMatrix().get(), false);
    AAF_[0] = Aaf_;
  }

  //================set up MLAPI smoothers for structure, fluid, ale on each level
  Teuchos::RCP<MLAPI::InverseOperator> S;
  Teuchos::RCP<MLAPI::LoadBalanceInverseOperator> lbS;
  if (sis_amg())
  {
    // structure
    for (int i = 0; i < snlevel_ - 1; ++i)
    {
      Teuchos::ParameterList p;
      Teuchos::ParameterList pushlist(sparams_.sublist("smoother: ifpack list"));
      char levelstr[19];
      sprintf(levelstr, "(level %d)", i);
      Teuchos::ParameterList& subp = sparams_.sublist("smoother: list " + std::string(levelstr));
      std::string type = "";
      select_mlapi_smoother(type, i, subp, p, pushlist);
      if (type == "ILU")
      {
        lbS = Teuchos::rcp(new MLAPI::LoadBalanceInverseOperator());
        wrap_ilu_smoother(sml, Ass_[i], *lbS, i);
        Sss_[i] = lbS;
      }
      else
      {
        S = Teuchos::rcp(new MLAPI::InverseOperator());
        S->Reshape(Ass_[i], type, p, &pushlist);
        Sss_[i] = S;
      }
    }

    // structure coarse grid:
    S = Teuchos::rcp(new MLAPI::InverseOperator());
    S->Reshape(Ass_[snlevel_ - 1], "Amesos-KLU");
    Sss_[snlevel_ - 1] = S;

    // dummy coarser then coarse grids
    for (int i = snlevel_; i < maxnlevel_; ++i)
    {
      Sss_[i] = Sss_[i - 1];
    }
  }
  else
  {
    // setup direct solver/ILU prec and do a dummy solve to create factorization/preconditioner
    const CORE::LINALG::SparseMatrix& Matrix00 =
        (hybridPrec_ == nullptr) ? Matrix(0, 0) : hybridPrec_->Matrix(0, 0);
    structuresolver_->Setup(Matrix00.EpetraMatrix());
    Teuchos::RCP<Epetra_Vector> b = Teuchos::rcp(new Epetra_Vector(Matrix00.RangeMap(), true));
    Teuchos::RCP<Epetra_Vector> x = Teuchos::rcp(new Epetra_Vector(Matrix00.DomainMap(), true));
    structuresolver_->Solve(Matrix00.EpetraMatrix(), x, b, true, true);
    srun_ = 1;  // a first solve has been performed
  }

  if (fis_amg())
  {
    // fluid
    for (int i = 0; i < fnlevel_ - 1; ++i)
    {
      Teuchos::ParameterList p;
      Teuchos::ParameterList pushlist(fparams_.sublist("smoother: ifpack list"));
      char levelstr[19];
      sprintf(levelstr, "(level %d)", i);
      Teuchos::ParameterList& subp = fparams_.sublist("smoother: list " + std::string(levelstr));
      std::string type = "";
      select_mlapi_smoother(type, i, subp, p, pushlist);
      if (blocksmoother_[i] == "Schur")  // Schur Complement block smoother
      {
        schur_complement_operator(Schurff_[i], Ass_[i], Aff_[i], Aaa_[i], ASF_[i], AFS_[i], AFA_[i],
            AAF_[i], schuromega_[i], structuresplit_);
        S = Teuchos::rcp(new MLAPI::InverseOperator());
        S->Reshape(Schurff_[i], type, p, &pushlist);
        Sff_[i] = S;
      }
      else  // BGS Smoother
      {
        if (type == "ILU")  // wrap existing ILU decomp from ML hierarchy
        {
          lbS = Teuchos::rcp(new MLAPI::LoadBalanceInverseOperator());
          wrap_ilu_smoother(fml, Aff_[i], *lbS, i);
          Sff_[i] = lbS;
        }
        else  // build new smoother
        {
          S = Teuchos::rcp(new MLAPI::InverseOperator());
          S->Reshape(Aff_[i], type, p, &pushlist);
          Sff_[i] = S;
        }
      }
    }
    // fluid coarse grid:
    S = Teuchos::rcp(new MLAPI::InverseOperator());
    if (blocksmoother_[fnlevel_ - 1] == "Schur")
    {
      const int i = fnlevel_ - 1;
      schur_complement_operator(Schurff_[i], Ass_[i], Aff_[i], Aaa_[i], ASF_[i], AFS_[i], AFA_[i],
          AAF_[i], schuromega_[i], structuresplit_);
      S->Reshape(Schurff_[fnlevel_ - 1], "Amesos-KLU");
    }
    else
    {
      S->Reshape(Aff_[fnlevel_ - 1], "Amesos-KLU");
    }
    Sff_[fnlevel_ - 1] = S;
    // dummy coarser then coarse grids
    for (int i = fnlevel_; i < maxnlevel_; ++i)
    {
      Sff_[i] = Sff_[i - 1];
    }
  }
  else
  {
    // setup direct solver/ILU prec and do a dummy solve to create factorization/preconditioner
    Teuchos::RCP<CORE::LINALG::MapExtractor> fsidofmapex = Teuchos::null;
    Teuchos::RCP<Epetra_Map> irownodes = Teuchos::null;
    const CORE::LINALG::SparseMatrix& Matrix11 =
        (hybridPrec_ == nullptr) ? Matrix(1, 1) : hybridPrec_->Matrix(1, 1);
    fluidsolver_->Setup(
        Matrix11.EpetraMatrix(), fsidofmapex, fluid_.discretization(), irownodes, structuresplit_);
    Teuchos::RCP<Epetra_Vector> b = Teuchos::rcp(new Epetra_Vector(Matrix11.RangeMap(), true));
    Teuchos::RCP<Epetra_Vector> x = Teuchos::rcp(new Epetra_Vector(Matrix11.DomainMap(), true));
    fluidsolver_->Solve(Matrix11.EpetraMatrix(), x, b, true, true);
    frun_ = 1;  // a first solve has been performed
  }

  if (ais_amg())
  {
    // ale
    for (int i = 0; i < anlevel_ - 1; ++i)
    {
      Teuchos::ParameterList p;
      Teuchos::ParameterList pushlist(aparams_.sublist("smoother: ifpack list"));
      char levelstr[19];
      sprintf(levelstr, "(level %d)", i);
      Teuchos::ParameterList& subp = aparams_.sublist("smoother: list " + std::string(levelstr));
      std::string type = "";
      select_mlapi_smoother(type, i, subp, p, pushlist);
      if (type == "ILU")
      {
        lbS = Teuchos::rcp(new MLAPI::LoadBalanceInverseOperator());
        wrap_ilu_smoother(aml, Aaa_[i], *lbS, i);
        Saa_[i] = lbS;
      }
      else
      {
        S = Teuchos::rcp(new MLAPI::InverseOperator());
        S->Reshape(Aaa_[i], type, p, &pushlist);
        Saa_[i] = S;
      }
    }
    // ale coarse grid:
    S = Teuchos::rcp(new MLAPI::InverseOperator());
    S->Reshape(Aaa_[anlevel_ - 1], "Amesos-KLU");
    Saa_[anlevel_ - 1] = S;
    // dummy coarser then coarse grids
    for (int i = anlevel_; i < maxnlevel_; ++i)
    {
      Saa_[i] = Saa_[i - 1];
    }
  }
  else
  {
    // setup direct solver/ILU prec and do a dummy solve to create factorization/preconditioner
    const CORE::LINALG::SparseMatrix& Matrix22 =
        (hybridPrec_ == nullptr) ? Matrix(2, 2) : hybridPrec_->Matrix(2, 2);
    alesolver_->Setup(Matrix22.EpetraMatrix());
    Teuchos::RCP<Epetra_Vector> b = Teuchos::rcp(new Epetra_Vector(Matrix22.RangeMap(), true));
    Teuchos::RCP<Epetra_Vector> x = Teuchos::rcp(new Epetra_Vector(Matrix22.DomainMap(), true));
    alesolver_->Solve(Matrix22.EpetraMatrix(), x, b, true, true);
    arun_ = 1;  // a first solve has been performed
  }

  //---------------------------------------------------------------------
  // in case we do FSIAMG now switch the minnlevel to maxnlevel since
  // we have set up dummy coarser then coarse grids for the smaller hierarchies
  minnlevel_ = maxnlevel_;

  //-------------------------------------------------------------- timing
  if (!myrank && verbosity_ == INPAR::FSI::verbosity_full)
  {
    printf("       -----------------------------------------------------------------------\n");
    printf("       Additional AMG(BGS/Schur)/ BGS(AMG) setup time %10.5e [s]\n",
        etime.totalElapsedTime(true));
  }

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void FSI::OverlappingBlockMatrixFSIAMG::schur_complement_operator(MLAPI::Operator& Schur,
    MLAPI::Operator& Ass, MLAPI::Operator& Aff, MLAPI::Operator& Aaa, MLAPI::Operator& Asf,
    MLAPI::Operator& Afs, MLAPI::Operator& Afa, MLAPI::Operator& Aaf, const double omega,
    const bool structuresplit)
{
  MLAPI::Operator FSSinvSF;  // structure contribution
  Epetra_CrsMatrix* diagSinv;
  const Epetra_Map& rmaps = Ass.GetRCPRowMatrix()->OperatorRangeMap();
  const Epetra_Map& dmaps = Ass.GetRCPRowMatrix()->OperatorDomainMap();
  {
    Epetra_Vector invdiagS(Ass.GetRCPRowMatrix()->OperatorRangeMap());
    Ass.GetRCPRowMatrix()->ExtractDiagonalCopy(invdiagS);
    int err = invdiagS.Reciprocal(invdiagS);
    if (err) FOUR_C_THROW("Inverse of diagonal of S returned %d", err);
    diagSinv = new Epetra_CrsMatrix(::Copy, rmaps, 1, true);
    for (int j = 0; j < rmaps.NumMyElements(); ++j)
    {
      int gid = rmaps.GID(j);
      err = diagSinv->InsertGlobalValues(gid, 1, &invdiagS[j], &gid);
      if (err) FOUR_C_THROW("Epetra_CrsMatrix::InsertGlobalValues returned %d", err);
    }
    diagSinv->FillComplete();
  }

  MLAPI::Operator FAAinvAF;  // ale contribution
  Epetra_CrsMatrix* diagAinv;
  const Epetra_Map& rmapa = Aaa.GetRCPRowMatrix()->OperatorRangeMap();
  const Epetra_Map& dmapa = Aaa.GetRCPRowMatrix()->OperatorDomainMap();
  {
    Epetra_Vector invdiagA(Aaa.GetRCPRowMatrix()->OperatorRangeMap());
    Aaa.GetRCPRowMatrix()->ExtractDiagonalCopy(invdiagA);
    int err = invdiagA.Reciprocal(invdiagA);
    if (err) FOUR_C_THROW("Inverse of diagonal of A returned %d", err);
    diagAinv = new Epetra_CrsMatrix(::Copy, rmapa, 1, true);
    for (int j = 0; j < rmapa.NumMyElements(); ++j)
    {
      int gid = rmapa.GID(j);
      err = diagAinv->InsertGlobalValues(gid, 1, &invdiagA[j], &gid);
      if (err) FOUR_C_THROW("Epetra_CrsMatrix::InsertGlobalValues returned %d", err);
    }
    diagAinv->FillComplete();
  }

  if (structuresplit)
  {
    {
      diagSinv->Scale(omega);
      MLAPI::Space rspace(rmaps);
      MLAPI::Space dspace(dmaps);
      MLAPI::Operator S(dspace, rspace, diagSinv);
      MLAPI::Operator tmp;
      tmp = S * Asf;
      FSSinvSF = Afs * tmp;
    }
    {
      diagAinv->Scale(omega);
      MLAPI::Space rspace(rmapa);
      MLAPI::Space dspace(dmapa);
      MLAPI::Operator A(dspace, rspace, diagAinv);
      MLAPI::Operator tmp;
      tmp = A * Aaf;
      FAAinvAF = Afa * tmp;
    }
    Schur = Aff - FSSinvSF;
    Schur = Schur - FAAinvAF;
  }
  else  // fluidsplit
  {
    {
      MLAPI::Space rspace_s(rmaps);
      MLAPI::Space dspace_s(dmaps);
      MLAPI::Operator S(dspace_s, rspace_s, diagSinv, false);  // need this again

      diagAinv->Scale(omega);
      MLAPI::Space rspace_a(rmapa);
      MLAPI::Space dspace_a(dmapa);
      MLAPI::Operator A(dspace_a, rspace_a, diagAinv);

      MLAPI::Operator tmp1;
      tmp1 = S * Asf;

      MLAPI::Operator tmp2;
      tmp2 = A * Aaf;

      MLAPI::Operator tmp3;
      tmp3 = tmp2 * tmp1;

      FAAinvAF = Afa * tmp3;
    }
    {
      diagSinv->Scale(omega);
      MLAPI::Space rspace_s(rmaps);
      MLAPI::Space dspace_s(dmaps);
      MLAPI::Operator S(dspace_s, rspace_s, diagSinv);
      MLAPI::Operator tmp;
      tmp = S * Asf;
      FSSinvSF = Afs * tmp;
    }
    Schur = Aff - FSSinvSF;
    Schur = Schur + FAAinvAF;
  }

  return;
}



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void FSI::OverlappingBlockMatrixFSIAMG::ra_poffdiagonals()
{
  Teuchos::RCP<Epetra_CrsMatrix> Matrix01 = (hybridPrec_ == nullptr)
                                                ? Matrix(0, 1).EpetraMatrix()
                                                : hybridPrec_->Matrix(0, 1).EpetraMatrix();
  Teuchos::RCP<Epetra_CrsMatrix> Matrix10 = (hybridPrec_ == nullptr)
                                                ? Matrix(1, 0).EpetraMatrix()
                                                : hybridPrec_->Matrix(1, 0).EpetraMatrix();
  Teuchos::RCP<Epetra_CrsMatrix> Matrix12 = (hybridPrec_ == nullptr)
                                                ? Matrix(1, 2).EpetraMatrix()
                                                : hybridPrec_->Matrix(1, 2).EpetraMatrix();
  Teuchos::RCP<Epetra_CrsMatrix> Matrix21 = (hybridPrec_ == nullptr)
                                                ? Matrix(2, 1).EpetraMatrix()
                                                : hybridPrec_->Matrix(2, 1).EpetraMatrix();
  Teuchos::RCP<Epetra_CrsMatrix> Matrix20 = (hybridPrec_ == nullptr)
                                                ? Matrix(2, 0).EpetraMatrix()
                                                : hybridPrec_->Matrix(2, 0).EpetraMatrix();

  for (int i = 0; i < maxnlevel_ - 1; ++i)
  {
    //------ Asf (trouble maker)
    if (!i)
      ra_pfine(ASF_[i + 1], Rss_[i], Matrix01, Pff_[i]);
    else
      ra_pcoarse(ASF_[i + 1], Rss_[i], ASF_[i], Pff_[i]);
    //------ Afs (trouble maker)
    if (!i)
      ra_pfine(AFS_[i + 1], Rff_[i], Matrix10, Pss_[i]);
    else
      ra_pcoarse(AFS_[i + 1], Rff_[i], AFS_[i], Pss_[i]);
    //------ Afa
    if (!i)
      ra_pfine(AFA_[i + 1], Rff_[i], Matrix12, Paa_[i]);
    else
      ra_pcoarse(AFA_[i + 1], Rff_[i], AFA_[i], Paa_[i]);
    //------ Aaf
    if (structuresplit_)
    {
      if (!i)
        ra_pfine(AAF_[i + 1], Raa_[i], Matrix21, Pff_[i]);
      else
        ra_pcoarse(AAF_[i + 1], Raa_[i], AAF_[i], Pff_[i]);
    }
    else
    {
      if (!i)
        ra_pfine(AAF_[i + 1], Raa_[i], Matrix20, Pss_[i]);
      else
        ra_pcoarse(AAF_[i + 1], Raa_[i], AAF_[i], Pss_[i]);
    }
  }
  return;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void FSI::OverlappingBlockMatrixFSIAMG::ra_pfine(MLAPI::Operator& RAP, const MLAPI::Operator& R,
    Teuchos::RCP<Epetra_CrsMatrix> A, const MLAPI::Operator& P)
{
  // this epetraext thingy patches the inherent ml<->epetra conflict
  EpetraExt::CrsMatrix_SolverMap transform;
  Epetra_CrsMatrix* Btrans = &(transform(*A));

  // down to the salt mines of ML....
  ML_Operator* mlB = ML_Operator_Create(MLAPI::GetML_Comm());
  ML_Operator_WrapEpetraMatrix(Btrans, mlB);
  ML_Operator* mlBP = ML_Operator_Create(MLAPI::GetML_Comm());
  ML_2matmult(mlB, P.GetML_Operator(), mlBP, ML_CSR_MATRIX);

  ML_Operator* mlRBP = ML_Operator_Create(MLAPI::GetML_Comm());
  ML_2matmult(R.GetML_Operator(), mlBP, mlRBP, ML_CSR_MATRIX);

  ML_Operator_Destroy(&mlB);
  ML_Operator_Destroy(&mlBP);

  // take ownership of coarse operator
  RAP.Reshape(P.GetDomainSpace(), R.GetRangeSpace(), mlRBP, true);

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void FSI::OverlappingBlockMatrixFSIAMG::ra_pcoarse(MLAPI::Operator& RAP, const MLAPI::Operator& R,
    const MLAPI::Operator& A, const MLAPI::Operator& P)
{
  MLAPI::Operator AP;
  // Intentionally do not use MLAPI's build in RAP product
  AP = A * P;
  RAP = R * AP;
  return;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void FSI::OverlappingBlockMatrixFSIAMG::wrap_ilu_smoother(
    ML* ml, MLAPI::Operator& A, MLAPI::LoadBalanceInverseOperator& S, const int level)
{
  // we use pre == post smoother here, so get postsmoother from ml
  void* data = ml->pre_smoother[level].smoother->data;
  Ifpack_Handle_Struct* handle = static_cast<Ifpack_Handle_Struct*>(data);
  bool takepart = false;
  if (handle->A_Base) takepart = true;
  Ifpack_Preconditioner* prec = nullptr;
  Epetra_RowMatrix* mat = nullptr;
  if (takepart)
  {
    prec = static_cast<Ifpack_Preconditioner*>(handle->A_Base);
    mat = const_cast<Epetra_RowMatrix*>(&(prec->Matrix()));
  }


  int NumMyElements = 0;
  int* MyGlobalElements = nullptr;
  if (takepart)
  {
    NumMyElements = mat->OperatorRangeMap().NumMyElements();
    MyGlobalElements = mat->OperatorRangeMap().MyGlobalElements();
  }
  MLAPI::Space range(-1, NumMyElements, MyGlobalElements);
  if (takepart)
  {
    NumMyElements = mat->OperatorDomainMap().NumMyElements();
    MyGlobalElements = mat->OperatorDomainMap().MyGlobalElements();
  }
  MLAPI::Space domain(-1, NumMyElements, MyGlobalElements);

  MLAPI::LoadBalanceOperator tmp;
  tmp.Reshape(domain, range, mat, false);
  S.Reshape(prec, tmp, false);

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void FSI::OverlappingBlockMatrixFSIAMG::select_mlapi_smoother(std::string& type, const int level,
    Teuchos::ParameterList& subp, Teuchos::ParameterList& p, Teuchos::ParameterList& pushlist)
{
  type = subp.get("smoother: type", "none");
  if (type == "none") FOUR_C_THROW("Cannot find msoother type");
  if (type == "symmetric Gauss-Seidel" || type == "Gauss-Seidel")
  {
    const int sweeps = subp.get("smoother: sweeps", 1);
    const double damping = subp.get("smoother: damping factor", 1.0);
    p.set("smoother: sweeps", sweeps);
    p.set("smoother: damping factor", damping);
  }
  else if (type == "IFPACK")
  {
    type = subp.get("smoother: ifpack type", "ILU");
    const double lof = subp.get<double>("smoother: ifpack level-of-fill", 0);
    const double damping = subp.get("smoother: damping factor", 1.0);
    p.set("smoother: ilu fill", (int)lof);
    p.set("smoother: damping factor", damping);
    p.set("schwarz: reordering type", "rcm");
    pushlist.set("ILU: sweeps", (int)lof);
    pushlist.set("fact: absolute threshold", 0.0);
    pushlist.set("fact: ict level-of-fill", lof);
    pushlist.set("fact: ilut level-of-fill", lof);
    pushlist.set("schwarz: reordering type", "rcm");
  }
  else if (type == "MLS")
  {
    const int poly = subp.get("smoother: MLS polynomial order", 3);
    p.set("smoother: MLS polynomial order", poly);
  }
  else if (type == "Amesos-KLU")
    ;  // nothing to do
  else if (type == "Amesos-Superludist")
    FOUR_C_THROW("No SuperLUDist support in MLAPI");
  else
    FOUR_C_THROW("Smoother not recognized");
  p.set("relaxation: zero starting solution", false);
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void FSI::OverlappingBlockMatrixFSIAMG::explicit_block_vcycle(const int level, const int nlevel,
    MLAPI::MultiVector& mlsy, MLAPI::MultiVector& mlfy, MLAPI::MultiVector& mlay,
    const MLAPI::MultiVector& mlsx, const MLAPI::MultiVector& mlfx,
    const MLAPI::MultiVector& mlax) const
{
  mlsy = 0.0;
  mlay = 0.0;
  mlfy = 0.0;

  // coarsest common level
  if (level == nlevel - 1)
  {
    // on the coarsest common level, use a BlockRichardson that
    // uses "the leftover peak" of the individual multigrid hierachies instead of
    // the simple smoothing schemes within the level. In case a field does not
    // have a remaining "leftover peak", direct solve will be called
    // automatically.
    explicit_block_gauss_seidel_smoother(level, mlsy, mlfy, mlay, mlsx, mlfx, mlax, true);
    return;
  }

  //-------------------------- presmoothing block Gauss Seidel
  explicit_block_gauss_seidel_smoother(level, mlsy, mlfy, mlay, mlsx, mlfx, mlax, false);

  //----------------------------------- coarse level residuals
  // structure
  // sxc = Rss_[level] * ( mlsx - Ass_[level] * mlsy - Asf[level] * mlfy)
  MLAPI::MultiVector sxc;
  {
    MLAPI::MultiVector tmp;
    tmp = mlsx - Ass_[level] * mlsy;
    tmp = tmp - ASF_[level] * mlfy;
    sxc.Reshape(Rss_[level].GetOperatorRangeSpace());
    sxc = Rss_[level] * tmp;
  }

  // ale
  // axc = Raa_[level] * ( mlax - Aaa_[level] * mlay - Aaf[level] * mlfy)
  MLAPI::MultiVector axc;
  {
    MLAPI::MultiVector tmp;
    tmp = mlax - Aaa_[level] * mlay;
    if (structuresplit_)
      tmp = tmp - AAF_[level] * mlfy;
    else
      tmp = tmp - AAF_[level] * mlsy;
    axc.Reshape(Raa_[level].GetOperatorRangeSpace());
    axc = Raa_[level] * tmp;
  }

  // fluid
  // fxc = Rff_[level] * ( mlfx - Aff_[level] * mlfy - Afs[level] * mlsy - Afa[level] * mlay)
  MLAPI::MultiVector fxc;
  {
    MLAPI::MultiVector tmp;
    tmp = mlfx - Aff_[level] * mlfy;
    tmp = tmp - AFS_[level] * mlsy;
    tmp = tmp - AFA_[level] * mlay;
    fxc.Reshape(Rff_[level].GetOperatorRangeSpace());
    fxc = Rff_[level] * tmp;
  }

  //----------------------------------- coarse level corrections
  MLAPI::MultiVector syc(sxc.GetVectorSpace(), 1, false);
  MLAPI::MultiVector ayc(axc.GetVectorSpace(), 1, false);
  MLAPI::MultiVector fyc(fxc.GetVectorSpace(), 1, false);

  //--------------------------------------- solve coarse problem
  explicit_block_vcycle(level + 1, nlevel, syc, fyc, ayc, sxc, fxc, axc);

  //------------------------------- prolongate coarse correction
  {
    MLAPI::MultiVector tmp;
    tmp = Pss_[level] * syc;
    mlsy.Update(1.0, tmp, 1.0);
    tmp = Paa_[level] * ayc;
    mlay.Update(1.0, tmp, 1.0);
    tmp = Pff_[level] * fyc;
    mlfy.Update(1.0, tmp, 1.0);
  }

  //---------------------------- postsmoothing block Gauss Seidel
  // (do NOT zero initial guess)
  explicit_block_gauss_seidel_smoother(level, mlsy, mlfy, mlay, mlsx, mlfx, mlax, false);

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void FSI::OverlappingBlockMatrixFSIAMG::explicit_block_gauss_seidel_smoother(const int level,
    MLAPI::MultiVector& mlsy, MLAPI::MultiVector& mlfy, MLAPI::MultiVector& mlay,
    const MLAPI::MultiVector& mlsx, const MLAPI::MultiVector& mlfx, const MLAPI::MultiVector& mlax,
    const bool amgsolve) const
{
  MLAPI::MultiVector sx;
  MLAPI::MultiVector fx;
  MLAPI::MultiVector ax;
  MLAPI::MultiVector sz(mlsy.GetVectorSpace(), 1, false);
  MLAPI::MultiVector fz(mlfy.GetVectorSpace(), 1, false);
  MLAPI::MultiVector az(mlay.GetVectorSpace(), 1, false);

  for (int run = 0; run < pciter_[level]; ++run)
  {
    //-------------- structure block
    {
      // compute ( r - A y ) for structure row
      sx = mlsx - Ass_[level] * mlsy;
      sx = sx - ASF_[level] * mlfy;
      // zero initial guess
      sz = 0.0;
      local_block_richardson(siterations_[level], somega_[level], level, amgsolve, snlevel_, sz, sx,
          Ass_, Sss_, Pss_, Rss_);
      // if (!amgsolve) Sss_[level].Apply(sx,sz);
      // else           Vcycle(level,snlevel_,sz,sx,Ass_,Sss_,Pss_,Rss_);
      mlsy.Update(pcomega_[level], sz, 1.0);
    }

    //-------------------- ale block
    {
      // compute ( r - A y ) for ale row
      ax = mlax - Aaa_[level] * mlay;
      if (structuresplit_)
        ax = ax - AAF_[level] * mlfy;
      else
        ax = ax - AAF_[level] * mlsy;
      // zero initial guess
      az = 0.0;
      local_block_richardson(aiterations_[level], aomega_[level], level, amgsolve, anlevel_, az, ax,
          Aaa_, Saa_, Paa_, Raa_);
      // if (!amgsolve) Saa_[level].Apply(ax,az);
      // else           Vcycle(level,anlevel_,az,ax,Aaa_,Saa_,Paa_,Raa_);
      mlay.Update(pcomega_[level], az, 1.0);
    }

    //------------------ fluid block
    {
      // compute ( r - A y ) for fluid row
      fx = mlfx - Aff_[level] * mlfy;
      fx = fx - AFS_[level] * mlsy;
      fx = fx - AFA_[level] * mlay;
      // zero initial guess
      fz = 0.0;
      local_block_richardson(fiterations_[level], fomega_[level], level, amgsolve, fnlevel_, fz, fx,
          Aff_, Sff_, Pff_, Rff_);
      // if (!amgsolve) Sff_[level].Apply(fx,fz);
      // else           Vcycle(level,fnlevel_,fz,fx,Aff_,Sff_,Pff_,Rff_);
      mlfy.Update(pcomega_[level], fz, 1.0);
    }

  }  // for (int run=0; run<pciter_[level]; ++run)

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void FSI::OverlappingBlockMatrixFSIAMG::local_block_richardson(const int iterations,
    const double omega, const int level, const bool amgsolve, const int nlevel,
    MLAPI::MultiVector& z, const MLAPI::MultiVector& b, const std::vector<MLAPI::Operator>& A,
    const std::vector<Teuchos::RCP<MLAPI::InverseOperator>>& S,
    const std::vector<MLAPI::Operator>& P, const std::vector<MLAPI::Operator>& R) const
{
  // do one iteration in any case (start counting iterations from zero)
  if (!amgsolve)
    S[level]->Apply(b, z);
  else
    vcycle(level, nlevel, z, b, A, S, P, R);

  if (iterations > 0)
  {
    MLAPI::MultiVector tmpz(z.GetVectorSpace(), 1, false);
    MLAPI::MultiVector tmpb(b.GetVectorSpace(), 1, false);
    z.Scale(omega);
    for (int i = 0; i < iterations; ++i)
    {
      tmpb = b - A[level] * z;
      tmpz = 0.0;
      if (!amgsolve)
        S[level]->Apply(tmpb, tmpz);
      else
        vcycle(level, nlevel, tmpz, tmpb, A, S, P, R);
      z = z + omega * tmpz;
    }
  }

  return;
}

/*----------------------------------------------------------------------*
strongly coupled AMG-Block-Gauss-Seidel
 *----------------------------------------------------------------------*/
void FSI::OverlappingBlockMatrixFSIAMG::sgs(
    const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  if (symmetric_) FOUR_C_THROW("FSIAMG symmetric Block Gauss-Seidel not impl.");

  // rewrap the matrix every time as it is killed irrespective
  // of whether the precond is reused or not.
  {
    Teuchos::RCP<Epetra_CrsMatrix> Matrix00 = (hybridPrec_ == nullptr)
                                                  ? Matrix(0, 0).EpetraMatrix()
                                                  : hybridPrec_->Matrix(0, 0).EpetraMatrix();
    MLAPI::Space dspace(Matrix00->DomainMap());
    MLAPI::Space rspace(Matrix00->RangeMap());
    Ass_[0].Reshape(dspace, rspace, Matrix00.get(), false);
  }
  {
    Teuchos::RCP<Epetra_CrsMatrix> Matrix01 = (hybridPrec_ == nullptr)
                                                  ? Matrix(0, 1).EpetraMatrix()
                                                  : hybridPrec_->Matrix(0, 1).EpetraMatrix();
    MLAPI::Space dspace(Matrix01->DomainMap());
    MLAPI::Space rspace(Matrix01->RangeMap());
    Asf_.Reshape(dspace, rspace, Matrix01.get(), false);
    ASF_[0] = Asf_;
  }
  {
    Teuchos::RCP<Epetra_CrsMatrix> Matrix10 = (hybridPrec_ == nullptr)
                                                  ? Matrix(1, 0).EpetraMatrix()
                                                  : hybridPrec_->Matrix(1, 0).EpetraMatrix();
    MLAPI::Space dspace(Matrix10->DomainMap());
    MLAPI::Space rspace(Matrix10->RangeMap());
    Afs_.Reshape(dspace, rspace, Matrix10.get(), false);
    AFS_[0] = Afs_;
  }
  {
    Teuchos::RCP<Epetra_CrsMatrix> Matrix11 = (hybridPrec_ == nullptr)
                                                  ? Matrix(1, 1).EpetraMatrix()
                                                  : hybridPrec_->Matrix(1, 1).EpetraMatrix();
    MLAPI::Space dspace(Matrix11->DomainMap());
    MLAPI::Space rspace(Matrix11->RangeMap());
    Aff_[0].Reshape(dspace, rspace, Matrix11.get(), false);
  }
  {
    Teuchos::RCP<Epetra_CrsMatrix> Matrix12 = (hybridPrec_ == nullptr)
                                                  ? Matrix(1, 2).EpetraMatrix()
                                                  : hybridPrec_->Matrix(1, 2).EpetraMatrix();
    MLAPI::Space dspace(Matrix12->DomainMap());
    MLAPI::Space rspace(Matrix12->RangeMap());
    Afa_.Reshape(dspace, rspace, Matrix12.get(), false);
    AFA_[0] = Afa_;
  }
  if (structuresplit_)
  {
    Teuchos::RCP<Epetra_CrsMatrix> Matrix21 = (hybridPrec_ == nullptr)
                                                  ? Matrix(2, 1).EpetraMatrix()
                                                  : hybridPrec_->Matrix(2, 1).EpetraMatrix();
    MLAPI::Space dspace(Matrix21->DomainMap());
    MLAPI::Space rspace(Matrix21->RangeMap());
    Aaf_.Reshape(dspace, rspace, Matrix21.get(), false);
    AAF_[0] = Aaf_;
  }
  else
  {
    Teuchos::RCP<Epetra_CrsMatrix> Matrix20 = (hybridPrec_ == nullptr)
                                                  ? Matrix(2, 0).EpetraMatrix()
                                                  : hybridPrec_->Matrix(2, 0).EpetraMatrix();
    MLAPI::Space dspace(Matrix20->DomainMap());
    MLAPI::Space rspace(Matrix20->RangeMap());
    Aaf_.Reshape(dspace, rspace, Matrix20.get(), false);
    AAF_[0] = Aaf_;
  }
  {
    Teuchos::RCP<Epetra_CrsMatrix> Matrix22 = (hybridPrec_ == nullptr)
                                                  ? Matrix(2, 2).EpetraMatrix()
                                                  : hybridPrec_->Matrix(2, 2).EpetraMatrix();
    MLAPI::Space dspace(Matrix22->DomainMap());
    MLAPI::Space rspace(Matrix22->RangeMap());
    Aaa_[0].Reshape(dspace, rspace, Matrix22.get(), false);
  }

  const Epetra_Vector& x = Teuchos::dyn_cast<const Epetra_Vector>(X);

  // various range and domain spaces
  const CORE::LINALG::SparseMatrix& Matrix00 =
      (hybridPrec_ == nullptr) ? Matrix(0, 0) : hybridPrec_->Matrix(0, 0);
  const CORE::LINALG::SparseMatrix& Matrix11 =
      (hybridPrec_ == nullptr) ? Matrix(1, 1) : hybridPrec_->Matrix(1, 1);
  const CORE::LINALG::SparseMatrix& Matrix22 =
      (hybridPrec_ == nullptr) ? Matrix(2, 2) : hybridPrec_->Matrix(2, 2);
  MLAPI::Space rsspace(Matrix00.RangeMap());
  MLAPI::Space rfspace(Matrix11.RangeMap());
  MLAPI::Space raspace(Matrix22.RangeMap());

  MLAPI::Space dsspace(Matrix00.DomainMap());
  MLAPI::Space dfspace(Matrix11.DomainMap());
  MLAPI::Space daspace(Matrix22.DomainMap());

  // initial guess has to be zero!
  Epetra_Vector& y = Teuchos::dyn_cast<Epetra_Vector>(Y);

  Teuchos::RCP<Epetra_Vector> sy = (hybridPrec_ == nullptr)
                                       ? RangeExtractor().ExtractVector(y, 0)
                                       : hybridPrec_->RangeExtractor().ExtractVector(y, 0);
  Teuchos::RCP<Epetra_Vector> fy = (hybridPrec_ == nullptr)
                                       ? RangeExtractor().ExtractVector(y, 1)
                                       : hybridPrec_->RangeExtractor().ExtractVector(y, 1);
  Teuchos::RCP<Epetra_Vector> ay = (hybridPrec_ == nullptr)
                                       ? RangeExtractor().ExtractVector(y, 2)
                                       : hybridPrec_->RangeExtractor().ExtractVector(y, 2);
  MLAPI::MultiVector mlsy(rsspace, sy->Pointers());
  MLAPI::MultiVector mlfy(rfspace, fy->Pointers());
  MLAPI::MultiVector mlay(raspace, ay->Pointers());
  mlsy = 0.0;
  mlfy = 0.0;
  mlay = 0.0;

  // rhs
  Teuchos::RCP<Epetra_Vector> sx = (hybridPrec_ == nullptr)
                                       ? DomainExtractor().ExtractVector(x, 0)
                                       : hybridPrec_->DomainExtractor().ExtractVector(x, 0);
  Teuchos::RCP<Epetra_Vector> fx = (hybridPrec_ == nullptr)
                                       ? DomainExtractor().ExtractVector(x, 1)
                                       : hybridPrec_->DomainExtractor().ExtractVector(x, 1);
  Teuchos::RCP<Epetra_Vector> ax = (hybridPrec_ == nullptr)
                                       ? DomainExtractor().ExtractVector(x, 2)
                                       : hybridPrec_->DomainExtractor().ExtractVector(x, 2);
  MLAPI::MultiVector mlsx(dsspace, sx->Pointers());
  MLAPI::MultiVector mlfx(dfspace, fx->Pointers());
  MLAPI::MultiVector mlax(daspace, ax->Pointers());


  // run FSIAMG
  switch (strategy_)
  {
    case INPAR::FSI::PreconditionedKrylov:
    {
      int myrank = X.Comm().MyPID();
      std::vector<int> Vsweeps(3, 1);
      std::vector<double> Vdamps(3, 1.0);
      Vsweeps[0] = pciter_[0];
      Vdamps[0] = pcomega_[0];
      Vsweeps[1] = pciter_[1];
      Vdamps[1] = pcomega_[1];
      Vsweeps[2] = pciter_[2];
      Vdamps[2] = pcomega_[2];

      AnalyzeBest shierachy(snlevel_);
      for (int i = 0; i < snlevel_; ++i)
      {
        shierachy.S()[i] = Sss_[i];
        shierachy.Type()[i] = "structure smoother";
        shierachy.Sweeps()[i] = siterations_[i];
        shierachy.Damp()[i] = somega_[i];
      }
      AnalyzeBest fhierachy(fnlevel_);
      for (int i = 0; i < fnlevel_; ++i)
      {
        fhierachy.S()[i] = Sff_[i];
        fhierachy.Type()[i] = "fluid smoother";
        fhierachy.Sweeps()[i] = fiterations_[i];
        fhierachy.Damp()[i] = fomega_[i];
      }
      AnalyzeBest ahierachy(anlevel_);
      for (int i = 0; i < anlevel_; ++i)
      {
        ahierachy.S()[i] = Saa_[i];
        ahierachy.Type()[i] = "ale smoother";
        ahierachy.Sweeps()[i] = aiterations_[i];
        ahierachy.Damp()[i] = aomega_[i];
      }

      if (sis_amg() && fis_amg() && ais_amg())
        richardson_bgs_v(myrank, 1, 1.0, Vsweeps, Vdamps, shierachy, fhierachy, ahierachy, mlsy,
            mlfy, mlay, mlsx, mlfx, mlax, Ass_, const_cast<std::vector<MLAPI::Operator>&>(Pss_),
            const_cast<std::vector<MLAPI::Operator>&>(Rss_), Aff_,
            const_cast<std::vector<MLAPI::Operator>&>(Pff_),
            const_cast<std::vector<MLAPI::Operator>&>(Rff_), Aaa_,
            const_cast<std::vector<MLAPI::Operator>&>(Paa_),
            const_cast<std::vector<MLAPI::Operator>&>(Raa_), ASF_, AFS_, AFA_, AAF_, true, false,
            true);
      else
        richardson_bgs_mixed(myrank, 1, 1.0, Vsweeps, Vdamps, sis_amg(), fis_amg(), ais_amg(),
            shierachy, fhierachy, ahierachy, mlsy, mlfy, mlay, mlsx, mlfx, mlax, Ass_,
            const_cast<std::vector<MLAPI::Operator>&>(Pss_),
            const_cast<std::vector<MLAPI::Operator>&>(Rss_), Aff_,
            const_cast<std::vector<MLAPI::Operator>&>(Pff_),
            const_cast<std::vector<MLAPI::Operator>&>(Rff_), Aaa_,
            const_cast<std::vector<MLAPI::Operator>&>(Paa_),
            const_cast<std::vector<MLAPI::Operator>&>(Raa_), ASF_, AFS_, AFA_, AAF_, true, false,
            true);

      break;
    }
    case INPAR::FSI::LinalgSolver:
    {
      // Do nothing. Will be done by CORE::LINALG::Solver internally.
      break;
    }
    default:
    {
      FOUR_C_THROW("Unknown type of preconditioner choice");
      break;
    }
  }

  // Note that mlsy, mlfy, mlay are views of sy, fy, ay, respectively.
  if (hybridPrec_ == nullptr)
  {
    RangeExtractor().InsertVector(*sy, 0, y);
    RangeExtractor().InsertVector(*fy, 1, y);
    RangeExtractor().InsertVector(*ay, 2, y);
  }
  else
  {
    hybridPrec_->RangeExtractor().InsertVector(*sy, 0, y);
    hybridPrec_->RangeExtractor().InsertVector(*fy, 1, y);
    hybridPrec_->RangeExtractor().InsertVector(*ay, 2, y);
  }

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void FSI::OverlappingBlockMatrixFSIAMG::vcycle(const int level, const int nlevel,
    MLAPI::MultiVector& z, const MLAPI::MultiVector& b, const std::vector<MLAPI::Operator>& A,
    const std::vector<Teuchos::RCP<MLAPI::InverseOperator>>& S,
    const std::vector<MLAPI::Operator>& P, const std::vector<MLAPI::Operator>& R,
    const bool trigger) const
{
  // in presmoothing, the initial guess has to be zero, we do this manually here.
  // in postsmoothing, the initial guess has to be nonzero. This is tricky, as
  // SGS smoothers assume nonzero initial guess, but ILU smoothers ALWAYS assume
  // zero guess. We circumvent this by reformulating the postsmoothing step (see below)
  // such that the initial guess can be zero by hand.

  // coarse solve
  if (level == nlevel - 1)
  {
    z = *(S[level]) * b;
    return;
  }

  // presmoothing (initial guess = 0)
  z = 0.0;
  S[level]->Apply(b, z);

  // coarse level residual and correction
  MLAPI::MultiVector bc;
  MLAPI::MultiVector zc(P[level].GetDomainSpace(), 1, true);

  // compute residual and restrict to coarser level
  bc = R[level] * (b - A[level] * z);

  // solve coarse problem
  vcycle(level + 1, nlevel, zc, bc, A, S, P, R);

  // prolongate correction
  z = z + P[level] * zc;


  // postsmoothing (initial guess != 0 !!)
  MLAPI::MultiVector r(b.GetVectorSpace(), true);
  MLAPI::MultiVector dz(b.GetVectorSpace(), true);
  r = A[level] * z;
  r.Update(1.0, b, -1.0);
  S[level]->Apply(r, dz);
  z.Update(1.0, dz, 1.0);


  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
const char* FSI::OverlappingBlockMatrixFSIAMG::Label() const
{
  if (strategy_ == INPAR::FSI::PreconditionedKrylov)
  {
    if (structuresplit_)
      return "structuresplit BGS(AMG) / PreconditionedKrylov";
    else
      return "fluidsplit BGS(AMG) / PreconditionedKrylov";
  }
  return "Unknown strategy in FSI::OverlappingBlockMatrixFSIAMG";
}

FOUR_C_NAMESPACE_CLOSE
