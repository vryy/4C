/*--------------------------------------------------------------------------*/
/*!
\file ssi_monolithic.cpp

\brief monolithic scalar-structure interaction

\level 2

<pre>
\maintainer Rui Fang & Christoph Schmidt
            {fang,schmidt}@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/
            089-289-15251
</pre>
*/
/*--------------------------------------------------------------------------*/
#include "ssi_monolithic.H"

#include <Epetra_Time.h>

#include "../drt_adapter/ad_str_ssiwrapper.H"
#include "../drt_adapter/adapter_scatra_base_algorithm.H"

#include "../drt_inpar/inpar_scatra.H"
#include "../drt_inpar/inpar_ssi.H"

#include "../drt_io/io_control.H"

#include "../drt_lib/drt_assemblestrategy.H"
#include "../drt_lib/drt_globalproblem.H"

#include "../drt_scatra/scatra_timint_implicit.H"

#include "../drt_scatra_ele/scatra_ele_action.H"

#include "../linalg/linalg_mapextractor.H"
#include "../linalg/linalg_solver.H"
#include "../linalg/linalg_utils.H"

/*--------------------------------------------------------------------------*
 | constructor                                                   fang 08/17 |
 *--------------------------------------------------------------------------*/
SSI::SSI_Mono::SSI_Mono(
    const Epetra_Comm&              comm,              //!< communicator
    const Teuchos::ParameterList&   globaltimeparams   //!< parameter list for time integration
    )
    // call base class constructor
  : SSI_Base(comm, globaltimeparams),

    // initialize member variables
    dtele_(0.),
    dtsolve_(0.),
    iter_(0),
    itermax_(globaltimeparams.get<int>("ITEMAX")),
    itertol_(globaltimeparams.sublist("MONOLITHIC").get<double>("CONVTOL")),
    maps_(Teuchos::null),
    matrixtype_(DRT::INPUT::IntegralValue<INPAR::SSI::MatrixType>(globaltimeparams.sublist("MONOLITHIC"),"MATRIXTYPE")),
    residual_(Teuchos::null),
    restol_(globaltimeparams.sublist("MONOLITHIC").get<double>("ABSTOLRES")),
    scatradofnorm_(0.),
    scatraincnorm_(0.),
    scatraresnorm_(0.),
    scatrastructureblock_(Teuchos::null),
    solver_(Teuchos::rcp(new LINALG::Solver(
        DRT::Problem::Instance()->SolverParams(globaltimeparams.sublist("MONOLITHIC").get<int>("LINEAR_SOLVER")),
        comm,
        DRT::Problem::Instance()->ErrorFile()->Handle()
        ))),
    structuredofnorm_(0.),
    structureresnorm_(0.),
    structureincnorm_(0.),
    structurescatrablock_(Teuchos::null),
    systemmatrix_(Teuchos::null),
    timer_(Teuchos::rcp(new Epetra_Time(comm)))
{
  return;
}


/*--------------------------------------------------------------------------*
 | assemble global system of equations                           fang 08/17 |
 *--------------------------------------------------------------------------*/
void SSI::SSI_Mono::AssembleMatAndRHS()
{
  // pass scalar transport degrees of freedom to structural discretization
  SetScatraSolution(scatra_->ScaTraField()->Phinp());

  // build system matrix and residual for structure field
  if(iter_ == 1)
    structure_->Evaluate();
  else
    structure_->Evaluate(maps_->ExtractVector(increment_,1));

  // apply Dirichlet conditions on structural system matrix
  structure_->SystemMatrix()->ApplyDirichlet(*structure_->GetDBCMapExtractor()->CondMap(),true);

  // pass structural degrees of freedom to scalar transport discretization
  SetStructSolution(structure_->Dispnp(),structure_->Velnp());

  // build system matrix and residual for scalar transport field
  scatra_->ScaTraField()->PrepareLinearSolve();

  // assemble off-diagonal scatra-structure block of global system matrix
  AssembleODBlockScatraStructure();

  //! assemble off-diagonal structure-scatra block of global system matrix
  AssembleODBlockStructureScatra();

  // build global system matrix
  switch(matrixtype_)
  {
    case INPAR::SSI::matrix_sparse:
    {
      // check global system matrix
      Teuchos::RCP<LINALG::SparseMatrix> systemmatrix = Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(systemmatrix_);
      if(systemmatrix == Teuchos::null)
        dserror("System matrix is not a sparse matrix!");

      // construct global system matrix by adding matrix blocks
      systemmatrix->Add(*scatra_->ScaTraField()->SystemMatrix(),false,1.,0.);
      systemmatrix->Add(*scatrastructureblock_,false,1.,1.);
      systemmatrix->Add(*structurescatrablock_,false,1.,1.);
      systemmatrix->Add(*structure_->SystemMatrix(),false,1.,1.);

      break;
    }

    default:
    {
      dserror("Type of global system matrix for scalar-structure interaction not recognized!");
      break;
    }
  }

  // finalize global system matrix
  systemmatrix_->Complete();

  // create monolithic right-hand side vector
  residual_->PutScalar(0.);
  maps_->InsertVector(scatra_->ScaTraField()->Residual(),0,residual_);
  maps_->AddVector(structure_->RHS(),1,residual_,-1.);

  return;
}


/*-----------------------------------------------------------------------------------*
 | assemble off-diagonal scatra-structure block of global system matrix   fang 08/17 |
 *-----------------------------------------------------------------------------------*/
void SSI::SSI_Mono::AssembleODBlockScatraStructure() const
{
  // initialize scatra-structure matrix block
  scatrastructureblock_->Zero();

  // create parameter list for element evaluation
  Teuchos::ParameterList eleparams;

  // action for elements
  eleparams.set<int>("action",SCATRA::calc_scatra_mono_odblock_mesh);

  // number of dofset associated with displacement-related dofs on scalar transport discretization
  eleparams.set<int>("ndsdisp",1);

  // number of dofset associated with velocity-related dofs on scalar transport discretization
  eleparams.set<int>("ndsvel",1);

  // remove state vectors from scalar transport discretization
  scatra_->ScaTraField()->Discretization()->ClearState();

  // add state vectors to scalar transport discretization
  scatra_->ScaTraField()->AddTimeIntegrationSpecificVectors();

  // create strategy for assembly of scatra-structure matrix block
  DRT::AssembleStrategy strategyscatrastructure(
      0,                       // row assembly based on number of dofset associated with scalar transport dofs on scalar transport discretization
      1,                       // column assembly based on number of dofset associated with structural dofs on scalar transport discretization
      scatrastructureblock_,   // scatra-structure matrix block
      Teuchos::null,           // no additional matrices or vectors
      Teuchos::null,
      Teuchos::null,
      Teuchos::null
      );

  // assemble scatra-structure matrix block
  scatra_->ScaTraField()->Discretization()->Evaluate(eleparams,strategyscatrastructure);

  // finalize scatra-structure matrix block
  scatrastructureblock_->Complete(*maps_->Map(1),*maps_->Map(0));

  // apply Dirichlet boundary conditions to scatra-structure matrix block
  scatrastructureblock_->ApplyDirichlet(*scatra_->ScaTraField()->DirichMaps()->CondMap(),false);

  // remove state vectors from scalar transport discretization
  scatra_->ScaTraField()->Discretization()->ClearState();

  return;
}


/*-----------------------------------------------------------------------------------*
 | assemble off-diagonal structure-scatra block of global system matrix   fang 08/17 |
 *-----------------------------------------------------------------------------------*/
void SSI::SSI_Mono::AssembleODBlockStructureScatra() const
{
  // nothing here at the moment
  structurescatrablock_->Complete(*maps_->Map(0),*maps_->Map(1));

  return;
}


/*--------------------------------------------------------------------------*
 | return global map of degrees of freedom                       fang 08/17 |
 *--------------------------------------------------------------------------*/
const Teuchos::RCP<const Epetra_Map>& SSI::SSI_Mono::DofRowMap() const
{
  return maps_->FullMap();
}


/*---------------------------------------------------------------------------------------------------*
 | pass structural degrees of freedom to scalar transport discretization and vice versa   fang 08/17 |
 *---------------------------------------------------------------------------------------------------*/
void SSI::SSI_Mono::ExchangeStateVectors()
{
  // pass structural degrees of freedom to scalar transport discretization and vice versa
  SetStructSolution(structure_->Dispnp(),structure_->Velnp());
  SetScatraSolution(scatra_->ScaTraField()->Phinp());

  return;
}


/*--------------------------------------------------------------------------*
 | check termination criterion for Newton-Raphson iteration      fang 08/17 |
 *--------------------------------------------------------------------------*/
bool SSI::SSI_Mono::ExitNewtonRaphson()
{
  // initialize exit flag
  bool exit(false);

  // compute vector norms for convergence check
  scatra_->ScaTraField()->Phinp()->Norm2(&scatradofnorm_);
  maps_->ExtractVector(residual_,0)->Norm2(&scatraresnorm_);
  maps_->ExtractVector(increment_,0)->Norm2(&scatraincnorm_);
  structure_->Dispnp()->Norm2(&structuredofnorm_);
  maps_->ExtractVector(residual_,1)->Norm2(&structureresnorm_);
  maps_->ExtractVector(increment_,1)->Norm2(&structureincnorm_);

  // safety checks
  if(std::isnan(scatradofnorm_) or
     std::isnan(scatraresnorm_) or
     std::isnan(scatraincnorm_) or
     std::isnan(structuredofnorm_) or
     std::isnan(structureresnorm_) or
     std::isnan(structureincnorm_))
    dserror("Vector norm is not a number!");
  if(std::isinf(scatradofnorm_) or
     std::isinf(scatraresnorm_) or
     std::isinf(scatraincnorm_) or
     std::isinf(structuredofnorm_) or
     std::isinf(structureresnorm_) or
     std::isinf(structureincnorm_))
    dserror("Vector norm is infinity!");

  // prevent division by zero
  if(scatradofnorm_ < 1.e-10)
    scatradofnorm_ = 1.e-10;
  if(structuredofnorm_ < 1.e-10)
    scatradofnorm_ = 1.e-10;

  // first Newton-Raphson iteration
  if(iter_ == 1)
  {
    // print first line of convergence table to screen
    // solution increment not yet available during first Newton-Raphson iteration
    if(Comm().MyPID() == 0)
      std::cout << "|  " << std::setw(3) << iter_ << "/" << std::setw(3) << itermax_ << "   | "
                << std::setw(10) << std::setprecision(3) << std::scientific << itertol_ << "[L_2 ]  | "
                << std::setw(10) << std::setprecision(3) << std::scientific << scatraresnorm_
                << "   |      --      | "
                << std::setw(10) << std::setprecision(3) << std::scientific << structureresnorm_
                << "   |      --      | "
                << "(       --      , te = "
                << std::setw(10) << std::setprecision(3) << dtele_ << ")" << std::endl;
  }

  // subsequent Newton-Raphson iterations
  else
  {
    // print current line of convergence table to screen
    if(Comm().MyPID() == 0)
      std::cout << "|  " << std::setw(3) << iter_ << "/" << std::setw(3) << itermax_ << "   | "
                << std::setw(10) << std::setprecision(3) << std::scientific << itertol_ << "[L_2 ]  | "
                << std::setw(10) << std::setprecision(3) << std::scientific << scatraresnorm_ << "   | "
                << std::setw(10) << std::setprecision(3) << std::scientific << scatraincnorm_/scatradofnorm_ << "   | "
                << std::setw(10) << std::setprecision(3) << std::scientific << structureresnorm_ << "   | "
                << std::setw(10) << std::setprecision(3) << std::scientific << structureincnorm_/structuredofnorm_ << "   | (ts = "
                << std::setw(10) << std::setprecision(3) << dtsolve_ << ", te = "
                << std::setw(10) << std::setprecision(3) << dtele_ << ")" << std::endl;

    // convergence check
    if(scatraresnorm_ <= itertol_ and
       structureresnorm_ <= itertol_ and
       scatraincnorm_/scatradofnorm_ <= itertol_ and
       structureincnorm_/structuredofnorm_ <= itertol_)
      // exit Newton-Raphson iteration upon convergence
      exit = true;
  }

  // exit Newton-Raphson iteration when residuals are small enough to prevent unnecessary additional solver calls
  if(scatraresnorm_ < restol_ and structureresnorm_ < restol_)
    exit = true;

  // print warning to screen if maximum number of Newton-Raphson iterations is reached without convergence
  if(iter_ == itermax_)
  {
    if(Comm().MyPID() == 0)
    {
      std::cout << "+------------+-------------------+--------------+--------------+--------------+--------------+" << std::endl;
      std::cout << "|      Newton-Raphson method has not converged after a maximum number of " << std::setw(2) << itermax_ << " iterations!      |" << std::endl;
    }

    // proceed to next time step
    exit = true;
  }

  // print finish line of convergence table to screen
  if(exit and Comm().MyPID() == 0)
    std::cout << "+------------+-------------------+--------------+--------------+--------------+--------------+" << std::endl;

  return exit;
}


/*--------------------------------------------------------------------------*
 | finite difference check for global system matrix              fang 08/17 |
 *--------------------------------------------------------------------------*/
void SSI::SSI_Mono::FDCheck()
{
  // initial screen output
  if(Comm().MyPID() == 0)
    std::cout << std::endl << "FINITE DIFFERENCE CHECK FOR SSI SYSTEM MATRIX" << std::endl;

  // create global state vector
  Teuchos::RCP<Epetra_Vector> statenp(LINALG::CreateVector(*DofRowMap(),true));
  maps_->InsertVector(scatra_->ScaTraField()->Phinp(),0,statenp);
  maps_->InsertVector(structure_->Dispnp(),1,statenp);

  // make a copy of global state vector to undo perturbations later
  Teuchos::RCP<Epetra_Vector> statenp_original = Teuchos::rcp(new Epetra_Vector(*statenp));

  // make a copy of system matrix as Epetra_CrsMatrix
  Teuchos::RCP<Epetra_CrsMatrix> sysmat_original = Teuchos::null;
  if(Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(systemmatrix_) != Teuchos::null)
    sysmat_original = (new LINALG::SparseMatrix(*(Teuchos::rcp_static_cast<LINALG::SparseMatrix>(systemmatrix_))))->EpetraMatrix();
  else if(Teuchos::rcp_dynamic_cast<LINALG::BlockSparseMatrixBase>(systemmatrix_) != Teuchos::null)
    sysmat_original = (new LINALG::SparseMatrix(*(Teuchos::rcp_static_cast<LINALG::BlockSparseMatrixBase>(systemmatrix_)->Merge())))->EpetraMatrix();
  else
    dserror("Type of system matrix unknown!");
  sysmat_original->FillComplete();

  // make a copy of system right-hand side vector
  Teuchos::RCP<Epetra_Vector> rhs_original = Teuchos::rcp(new Epetra_Vector(*residual_));

  // copy and zero out system increment vector if necessary
  Teuchos::RCP<Epetra_Vector> increment_original(Teuchos::null);
  if(iter_ != 1)
  {
    increment_original = Teuchos::rcp(new Epetra_Vector(*increment_));
    increment_->PutScalar(0.);
  }

  // initialize counter for system matrix entries with failing finite difference check
  int counter(0);

  // initialize tracking variable for maximum absolute and relative errors
  double maxabserr(0.);
  double maxrelerr(0.);

  for (int colgid=0; colgid<=sysmat_original->ColMap().MaxAllGID(); ++colgid)
  {
    // continue loop if current column index is not a valid global column index
    int collid(sysmat_original->ColMap().LID(colgid));
    int maxcollid(-1);
    Comm().MaxAll(&collid,&maxcollid,1);
    if(maxcollid < 0)
      continue;

    // fill global state vector with original state variables
    statenp->Update(1.,*statenp_original,0.);

    // impose perturbation
    if(statenp->Map().MyGID(colgid))
      if(statenp->SumIntoGlobalValue(colgid,0,scatra_->ScaTraField()->FDCheckEps()))
        dserror("Perturbation could not be imposed on state vector for finite difference check!");
    scatra_->ScaTraField()->Phinp()->Update(1.,*maps_->ExtractVector(statenp,0),0.);
    structure_->SetState(maps_->ExtractVector(statenp,1));

    // calculate element right-hand side vector for perturbed state
    AssembleMatAndRHS();

    // Now we compare the difference between the current entries in the system matrix
    // and their finite difference approximations according to
    // entries ?= (residual_perturbed - residual_original) / epsilon

    // Note that the residual_ vector actually denotes the right-hand side of the linear
    // system of equations, i.e., the negative system residual.
    // To account for errors due to numerical cancellation, we additionally consider
    // entries + residual_original / epsilon ?= residual_perturbed / epsilon

    // Note that we still need to evaluate the first comparison as well. For small entries in the system
    // matrix, the second comparison might yield good agreement in spite of the entries being wrong!
    for(int rowlid=0; rowlid<DofRowMap()->NumMyElements(); ++rowlid)
    {
      // get global index of current matrix row
      const int rowgid = sysmat_original->RowMap().GID(rowlid);
      if(rowgid < 0)
        dserror("Invalid global ID of matrix row!");

      // skip matrix rows associated with Dirichlet boundary conditions
      if(scatra_->ScaTraField()->DirichMaps()->CondMap()->MyGID(rowgid) or structure_->GetDBCMapExtractor()->CondMap()->MyGID(rowgid))
        continue;

      // get relevant entry in current row of original system matrix
      double entry(0.);
      int length = sysmat_original->NumMyEntries(rowlid);
      int numentries;
      std::vector<double> values(length);
      std::vector<int> indices(length);
      sysmat_original->ExtractMyRowCopy(rowlid,length,numentries,&values[0],&indices[0]);
      for(int ientry=0; ientry<length; ++ientry)
      {
        if(sysmat_original->ColMap().GID(indices[ientry]) == colgid)
        {
          entry = values[ientry];
          break;
        }
      }

      // finite difference suggestion (first divide by epsilon and then add for better conditioning)
      const double fdval = -(*residual_)[rowlid] / scatra_->ScaTraField()->FDCheckEps() + (*rhs_original)[rowlid] / scatra_->ScaTraField()->FDCheckEps();

      // confirm accuracy of first comparison
      if(abs(fdval) > 1.e-17 and abs(fdval) < 1.e-15)
        dserror("Finite difference check involves values too close to numerical zero!");

      // absolute and relative errors in first comparison
      const double abserr1 = entry - fdval;
      if(abs(abserr1) > maxabserr)
        maxabserr = abs(abserr1);
      double relerr1(0.);
      if(abs(entry) > 1.e-17)
        relerr1 = abserr1 / abs(entry);
      else if(abs(fdval) > 1.e-17)
        relerr1 = abserr1 / abs(fdval);
      if(abs(relerr1) > maxrelerr)
        maxrelerr = abs(relerr1);

      // evaluate first comparison
      if(abs(relerr1) > scatra_->ScaTraField()->FDCheckTol())
      {
        std::cout << "sysmat[" << rowgid << "," << colgid << "]:  " << entry << "   ";
        std::cout << "finite difference suggestion:  " << fdval << "   ";
        std::cout << "absolute error:  " << abserr1 << "   ";
        std::cout << "relative error:  " << relerr1 << std::endl;

        counter++;
      }

      // first comparison OK
      else
      {
        // left-hand side in second comparison
        const double left  = entry - (*rhs_original)[rowlid] / scatra_->ScaTraField()->FDCheckEps();

        // right-hand side in second comparison
        const double right = -(*residual_)[rowlid] / scatra_->ScaTraField()->FDCheckEps();

        // confirm accuracy of second comparison
        if(abs(right) > 1.e-17 and abs(right) < 1.e-15)
          dserror("Finite difference check involves values too close to numerical zero!");

        // absolute and relative errors in second comparison
        const double abserr2 = left - right;
        if(abs(abserr2) > maxabserr)
          maxabserr = abs(abserr2);
        double relerr2(0.);
        if(abs(left) > 1.e-17)
          relerr2 = abserr2 / abs(left);
        else if(abs(right) > 1.e-17)
          relerr2 = abserr2 / abs(right);
        if(abs(relerr2) > maxrelerr)
          maxrelerr = abs(relerr2);

        // evaluate second comparison
        if(abs(relerr2) > scatra_->ScaTraField()->FDCheckTol())
        {
          std::cout << "sysmat[" << rowgid << "," << colgid << "]-rhs[" << rowgid << "]/eps:  " << left << "   ";
          std::cout << "-rhs_perturbed[" << rowgid << "]/eps:  " << right << "   ";
          std::cout << "absolute error:  " << abserr2 << "   ";
          std::cout << "relative error:  " << relerr2 << std::endl;

          counter++;
        }
      }
    }
  }

  // communicate tracking variables
  int counterglobal(0);
  Comm().SumAll(&counter,&counterglobal,1);
  double maxabserrglobal(0.);
  Comm().MaxAll(&maxabserr,&maxabserrglobal,1);
  double maxrelerrglobal(0.);
  Comm().MaxAll(&maxrelerr,&maxrelerrglobal,1);

  // final screen output
  if(Comm().MyPID() == 0)
  {
    if(counterglobal)
    {
      printf("--> FAILED AS LISTED ABOVE WITH %d CRITICAL MATRIX ENTRIES IN TOTAL\n\n",counterglobal);
      dserror("Finite difference check failed for SSI system matrix!");
    }
    else
      printf("--> PASSED WITH MAXIMUM ABSOLUTE ERROR %+12.5e AND MAXIMUM RELATIVE ERROR %+12.5e\n\n",maxabserrglobal,maxrelerrglobal);
  }

  // undo perturbations of state variables
  scatra_->ScaTraField()->Phinp()->Update(1.,*maps_->ExtractVector(statenp_original,0),0.);
  structure_->SetState(maps_->ExtractVector(statenp_original,1));

  // recompute system matrix and right-hand side vector based on original state variables
  AssembleMatAndRHS();

  // restore system increment vector if necessary
  if(increment_original != Teuchos::null)
    increment_ = increment_original;

  return;
}


/*--------------------------------------------------------------------------*
 | initialize monolithic algorithm                               fang 08/17 |
 *--------------------------------------------------------------------------*/
int SSI::SSI_Mono::Init(
    const Epetra_Comm&              comm,               //!< communicator
    const Teuchos::ParameterList&   globaltimeparams,   //!< parameter list for time integration
    const Teuchos::ParameterList&   scatraparams,       //!< parameter list for scalar transport
    const Teuchos::ParameterList&   structparams,       //!< parameter list for structure
    const std::string               struct_disname,     //!< name of structural discretization
    const std::string               scatra_disname,     //!< name of scalar transport discretization
    bool                            isAle               //!< flag for ALE
    )
{
  // check input parameters for scalar transport field
  if(DRT::INPUT::IntegralValue<INPAR::SCATRA::VelocityField>(scatraparams,"VELOCITYFIELD") != INPAR::SCATRA::velocity_Navier_Stokes)
    dserror("Invalid type of velocity field for scalar-structure interaction!");

  // call base class routine
  return SSI_Base::Init(comm,globaltimeparams,scatraparams,structparams,struct_disname,scatra_disname,isAle);
}


/*--------------------------------------------------------------------------*
 | output solution to screen and files                           fang 08/17 |
 *--------------------------------------------------------------------------*/
void SSI::SSI_Mono::Output()
{
  // output scalar transport field
  scatra_->ScaTraField()->Output();

  // output structure field
  structure_->Output();

  return;
}


/*--------------------------------------------------------------------------*
 | prepare time step                                             fang 08/17 |
 *--------------------------------------------------------------------------*/
void SSI::SSI_Mono::PrepareTimeStep()
{
  // update time and time step
  IncrementTimeAndStep();

  // pass structural degrees of freedom to scalar transport discretization and vice versa
  ExchangeStateVectors();

  // prepare time step for scalar transport field
  scatra_->ScaTraField()->PrepareTimeStep();

  // prepare time step for structural field
  structure_->PrepareTimeStep();

  // print time step information to screen
  scatra_->ScaTraField()->PrintTimeStepInfo();

  return;
}


/*--------------------------------------------------------------------------*
 | setup monolithic algorithm                                    fang 08/17 |
 *--------------------------------------------------------------------------*/
void SSI::SSI_Mono::Setup()
{
  // call base class routine
  SSI_Base::Setup();

  // check maps from scalar transport and structure discretizations
  if(scatra_->ScaTraField()->DofRowMap()->NumGlobalElements() == 0)
    dserror("Scalar transport discretization does not have any degrees of freedom!");
  if(structure_->DofRowMap()->NumGlobalElements() == 0)
    dserror("Structure discretization does not have any degrees of freedom!");

  // additional safety checks
  if(!scatra_->ScaTraField()->IsIncremental())
    dserror("Must have incremental solution approach for scalar-structure interaction!");
  if(scatra_->ScaTraField()->SystemMatrix() == Teuchos::null)
    dserror("System matrix associated with scalar transport field must be a sparse matrix!");

  return;
}


/*--------------------------------------------------------------------------*
 | setup global system of equations                              fang 08/17 |
 *--------------------------------------------------------------------------*/
void SSI::SSI_Mono::SetupSystem()
{
  // initialize global map extractor
  maps_ = Teuchos::rcp(new LINALG::MapExtractor(
      *LINALG::MergeMap(
          *scatra_->ScaTraField()->DofRowMap(),
          *structure_->DofRowMap(),
          false
          ),
      structure_->DofRowMap(),
      scatra_->ScaTraField()->DofRowMap()
      ));

  // check global map extractor
  maps_->CheckForValidMapExtractor();

  // initialize global increment vector for Newton-Raphson iteration
  increment_ = LINALG::CreateVector(
      *DofRowMap(),
      true
      );

  // initialize global residual vector
  residual_ = LINALG::CreateVector(
      *DofRowMap(),
      true
      );

  // perform initializations associated with global system matrix
  switch(matrixtype_)
  {
    case INPAR::SSI::matrix_sparse:
    {
      // safety check
      if(scatra_->ScaTraField()->SystemMatrix() == Teuchos::null)
        dserror("Incompatible matrix type associated with scalar transport field!");

      // initialize global system matrix
      systemmatrix_ = Teuchos::rcp(new LINALG::SparseMatrix(*DofRowMap(),27,false,true));

      break;
    }

    default:
    {
      dserror("Type of global system matrix for scalar-structure interaction not recognized!");
      break;
    }
  }

  // initialize scatra-structure block
  scatrastructureblock_ = Teuchos::rcp(new LINALG::SparseMatrix(
      *scatra_->ScaTraField()->DofRowMap(),
      27,
      false,
      true
      ));

  // initialize structure-scatra block
  structurescatrablock_ = Teuchos::rcp(new LINALG::SparseMatrix(
      *structure_->DofRowMap(),
      27,
      false,
      true
      ));

  return;
}


/*--------------------------------------------------------------------------*
 | evaluate time step using Newton-Raphson iteration             fang 08/17 |
 *--------------------------------------------------------------------------*/
void SSI::SSI_Mono::Solve()
{
  // initialize counter for Newton-Raphson iteration
  iter_ = 0;

  // print header of convergence table to screen
  if(Comm().MyPID() == 0)
  {
    std::cout << "+------------+-------------------+--------------+--------------+--------------+--------------+" << std::endl;
    std::cout << "|- step/max -|- tolerance[norm] -|- scatra-res -|- scatra-inc -|- struct-res -|- struct-inc -|" << std::endl;
  }

  // start Newton-Raphson iteration
  while(true)
  {
    // update iteration counter
    ++iter_;

    // reset timer
    timer_->ResetStartTime();

    // store time before evaluating elements and assembling global system of equations
    double time = timer_->WallTime();

    // assemble global system of equations
    AssembleMatAndRHS();

    // determine time needed for evaluating elements and assembling global system of equations,
    // and take maximum over all processors via communication
    double mydtele = timer_->WallTime()-time;
    Comm().MaxAll(&mydtele,&dtele_,1);

    // safety check
    if(!systemmatrix_->Filled())
      dserror("Complete() has not been called on global system matrix yet!");

    // perform finite difference check on time integrator level
    if(scatra_->ScaTraField()->FDCheckType() == INPAR::SCATRA::fdcheck_global)
      FDCheck();

    // check termination criterion for Newton-Raphson iteration
    if(ExitNewtonRaphson())
      break;

    // initialize global increment vector
    increment_->PutScalar(0.);

    // store time before solving global system of equations
    time = timer_->WallTime();

    // solve global system of equations
    // Dirichlet boundary conditions have already been applied to global system of equations
    solver_->Solve(
        systemmatrix_->EpetraOperator(),
        increment_,
        residual_,
        true,
        iter_==1
        );

    // determine time needed for solving global system of equations,
    // and take maximum over all processors via communication
    double mydtsolve = timer_->WallTime()-time;
    Comm().MaxAll(&mydtsolve,&dtsolve_,1);

    // update scalar transport field
    scatra_->ScaTraField()->UpdateIter(maps_->ExtractVector(increment_,0));
    scatra_->ScaTraField()->ComputeIntermediateValues();

    // structure field is updated during the next Newton-Raphson iteration step
  } // Newton-Raphson iteration

  return;
}


/*--------------------------------------------------------------------------*
 | time loop                                                     fang 08/17 |
 *--------------------------------------------------------------------------*/
void SSI::SSI_Mono::Timeloop()
{
  // output initial solution to screen and files
  if(Step() == 0)
  {
    SetStructSolution(structure_->Dispnp(),structure_->Velnp());
    Output();
  }

  // time loop
  while(NotFinished())
  {
    // prepare time step
    PrepareTimeStep();

    // evaluate time step
    Solve();

    // update scalar transport and structure fields
    Update();

    // output solution to screen and files
    Output();
  } // while(NotFinished())

  return;
}


/*--------------------------------------------------------------------------------------*
 | update scalar transport and structure fields after time step evaluation   fang 08/17 |
 *--------------------------------------------------------------------------------------*/
void SSI::SSI_Mono::Update()
{
  // update scalar transport field
  scatra_->ScaTraField()->Update();

  // update structure field
  structure_->Update();

  return;
}
