/*----------------------------------------------------------------------*/
/*!
 \file poro_partitioned.cpp

 \brief  Partitioned poroelasticity algorithm

\level 2

\maintainer Ager Christoph
            ager@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289 15249
 *-----------------------------------------------------------------------*/

#include "poro_partitioned.H"

#include "../drt_adapter/ad_str_fpsiwrapper.H"
#include "../drt_adapter/ad_fld_poro.H"
#include "../drt_structure/stru_aux.H"

#include "../drt_lib/drt_globalproblem.H"

#include "../linalg/linalg_utils.H"

/*----------------------------------------------------------------------*
 | constructor (public)                                    vuong 01/12  |
 *----------------------------------------------------------------------*/
POROELAST::Partitioned::Partitioned(
    const Epetra_Comm& comm, const Teuchos::ParameterList& timeparams)
    : PoroBase(comm, timeparams),
      fluidincnp_(Teuchos::rcp(new Epetra_Vector(*(FluidField()->Velnp())))),
      structincnp_(Teuchos::rcp(new Epetra_Vector(*(StructureField()->Dispnp())))),
      del_(Teuchos::null),
      delhist_(Teuchos::null),
      omegan_(0.0),
      omeganp_(0.0)
{
  const Teuchos::ParameterList& porodyn = DRT::Problem::Instance()->PoroelastDynamicParams();
  // Get the parameters for the ConvergenceCheck
  itmax_ = porodyn.get<int>("ITEMAX");     // default: =10
  ittol_ = porodyn.get<double>("INCTOL");  // default: =1e-6

  fluidveln_ = LINALG::CreateVector(*(FluidField()->DofRowMap()), true);
  fluidveln_->PutScalar(0.0);
}

/*----------------------------------------------------------------------*
                                                           vuong 01/12  |
*----------------------------------------------------------------------*/
void POROELAST::Partitioned::DoTimeStep()
{
  PrepareTimeStep();

  Solve();

  UpdateAndOutput();
}

/*----------------------------------------------------------------------*
                                                           vuong 01/12  |
*----------------------------------------------------------------------*/
void POROELAST::Partitioned::SetupSystem() {}  // SetupSystem()

/*----------------------------------------------------------------------*
                                                           vuong 01/12  |
*----------------------------------------------------------------------*/
void POROELAST::Partitioned::UpdateAndOutput()
{
  PrepareOutput();

  Update();

  Output();
}  // UpdateAndOutput()

/*----------------------------------------------------------------------*
                                                           vuong 01/12  |
*----------------------------------------------------------------------*/
void POROELAST::Partitioned::Solve()
{
  int itnum = 0;
  bool stopnonliniter = false;

  if (Comm().MyPID() == 0)
  {
    std::cout << "\n****************************************\n          OUTER ITERATION "
                 "LOOP\n****************************************\n";
  }

  // initially solve coupled scalar transport equation system
  // DoScatraStep();

  if (Step() == 1)
  {
    fluidveln_->Update(1.0, *(FluidField()->Veln()), 0.0);
  }

  while (stopnonliniter == false)
  {
    itnum++;

    // store increment from first solution for convergence check (like in
    // elch_algorithm: use current values)
    fluidincnp_->Update(1.0, *FluidField()->Velnp(), 0.0);
    structincnp_->Update(1.0, *StructureField()->Dispnp(), 0.0);

    // get current fluid velocities due to solve fluid step, like predictor in FSI
    // 1. iteration: get velocities of old time step (T_n)
    if (itnum == 1)
    {
      fluidveln_->Update(1.0, *(FluidField()->Veln()), 0.0);
    }
    else  // itnum > 1
    {
      // save velocity solution of old iteration step T_{n+1}^i
      fluidveln_->Update(1.0, *(FluidField()->Velnp()), 0.0);
    }

    // set fluid- and structure-based scalar transport values required in FSI
    SetFluidSolution();

    if (itnum != 1) StructureField()->PreparePartitionStep();
    // solve structural system
    DoStructStep();

    // set mesh displacement and velocity fields
    SetStructSolution();

    // solve scalar transport equation
    DoFluidStep();
    // ScatraEvaluateSolveIterUpdate();

    // check convergence for all fields and stop iteration loop if
    // convergence is achieved overall
    stopnonliniter = ConvergenceCheck(itnum);

    // AitkenRelax();
  }

  // initial guess for next time step n+1
  // use maximum between omeganp_ and 1.0 as start value for omega_n+1^{i=0}
  // in case of doubt use 1.0, meaning that direction of new solution vector
  // dispnp better old one
  omegan_ = std::max(omeganp_, 1.0);

  return;
}  // OuterLoop()

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void POROELAST::Partitioned::DoStructStep()
{
  if (Comm().MyPID() == 0)
  {
    std::cout << "\n***********************\n STRUCTURE SOLVER \n***********************\n";
  }

  // Newton-Raphson iteration
  StructureField()->Solve();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void POROELAST::Partitioned::DoFluidStep()
{
  if (Comm().MyPID() == 0)
  {
    std::cout << "\n***********************\n FLUID SOLVER \n***********************\n";
  }

  // FluidField()->PrepareSolve();
  // Newton-Raphson iteration
  FluidField()->Solve();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void POROELAST::Partitioned::PrepareTimeStep()
{
  IncrementTimeAndStep();
  PrintHeader();

  StructureField()->PrepareTimeStep();
  SetStructSolution();
  FluidField()->PrepareTimeStep();
  SetFluidSolution();
}

/*----------------------------------------------------------------------*
 | convergence check for both fields (fluid & structure)
 *----------------------------------------------------------------------*/
bool POROELAST::Partitioned::ConvergenceCheck(int itnum)
{
  // convergence check based on the temperature increment
  bool stopnonliniter = false;

  //    | temperature increment |_2
  //  -------------------------------- < Tolerance
  //     | temperature_n+1 |_2

  // variables to save different L2 - Norms
  // define L2-norm of incremental temperature and temperature
  // here: only the temperature field is checked for convergence!!!
  double fluidincnorm_L2(0.0);
  double fluidnorm_L2(0.0);
  double dispincnorm_L2(0.0);
  double structnorm_L2(0.0);

  // build the current temperature increment Inc T^{i+1}
  // \f Delta T^{k+1} = Inc T^{k+1} = T^{k+1} - T^{k}  \f
  fluidincnp_->Update(1.0, *(FluidField()->Velnp()), -1.0);
  structincnp_->Update(1.0, *(StructureField()->Dispnp()), -1.0);

  // build the L2-norm of the temperature increment and the temperature
  fluidincnp_->Norm2(&fluidincnorm_L2);
  FluidField()->Velnp()->Norm2(&fluidnorm_L2);
  structincnp_->Norm2(&dispincnorm_L2);
  StructureField()->Dispnp()->Norm2(&structnorm_L2);

  // care for the case that there is (almost) zero temperature
  // (usually not required for temperature)
  if (fluidnorm_L2 < 1e-6) fluidnorm_L2 = 1.0;
  if (structnorm_L2 < 1e-6) structnorm_L2 = 1.0;

  // print the incremental based convergence check to the screen
  if (Comm().MyPID() == 0)  // and PrintScreenEvry() and (Step()%PrintScreenEvry()==0))
  {
    std::cout << "\n";
    std::cout
        << "***********************************************************************************\n";
    std::cout << "    OUTER ITERATION STEP    \n";
    std::cout
        << "***********************************************************************************\n";
    printf("+--------------+------------------------+--------------------+--------------------+\n");
    printf(
        "|-  step/max  -|-  tol      [norm]     -|--  fluid-inc      --|--  disp-inc      --|\n");
    printf("|   %3d/%3d    |  %10.3E[L_2 ]      | %10.3E         | %10.3E         |", itnum, itmax_,
        ittol_, fluidincnorm_L2 / fluidnorm_L2, dispincnorm_L2 / structnorm_L2);
    printf("\n");
    printf("+--------------+------------------------+--------------------+--------------------+\n");
  }

  // converged
  if ((fluidincnorm_L2 / fluidnorm_L2 <= ittol_) && (dispincnorm_L2 / structnorm_L2 <= ittol_))
  {
    stopnonliniter = true;
    if (Comm().MyPID() == 0)  // and PrintScreenEvry() and (Step()%PrintScreenEvry()==0))
    {
      printf("\n");
      printf(
          "|  Outer Iteration loop converged after iteration %3d/%3d !                       |\n",
          itnum, itmax_);
      printf(
          "+--------------+------------------------+--------------------+--------------------+\n");
    }
  }

  // warn if itemax is reached without convergence, but proceed to next
  // timestep
  if ((itnum == itmax_) and
      ((fluidincnorm_L2 / fluidnorm_L2 > ittol_) || (dispincnorm_L2 / structnorm_L2 > ittol_)))
  {
    stopnonliniter = true;
    if ((Comm().MyPID() == 0))  // and PrintScreenEvry() and (Step()%PrintScreenEvry()==0))
    {
      printf(
          "|     >>>>>> not converged in itemax steps!                                       |\n");
      printf(
          "+--------------+------------------------+--------------------+--------------------+\n");
      printf("\n");
      printf("\n");
    }
  }

  return stopnonliniter;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void POROELAST::Partitioned::AitkenRelax()
{
  dserror("Aitken Relaxation not yet implemented");
  // ------------------------------------------------------------
  // if r_{i+1} not converged in ConvergenceCheck()
  // --> apply Aitken relaxation to displacements
  // as implemented in FSI by Irons and Tuck (1969)

  // Aitken factor
  // relaxation can start when two vectors are available, starting at i=2
  // (also stated in Diss Uli, p. 82ff)
  // we need vectors from iteration step 0,1,2 to calculate relaxed step
  // with omega^i_0 = 0.0
  // Irons & Tuck:
  // increment r := old - new

  // BACI:
  // increment inc := new - old
  // omega^{i+1}_{n+1} = omega^i_{n+1} +
  //
  //                         ( r^i_{n+1} - r^{i+1}_{n+1} )^T . ( r^{i+1}_{n+1} )
  // + (omega^i_{n+1} - 1) . -------------------------------------------------
  //                                 | r^i_{n+1} - r^{i+1}_{n+1} |^2
  //
  // | r^i_{n+1} - r^{i+1}_{n+1} |^2
  // = | (-1)^2*(r^{i+1}_{n+1} - r^i_{n+1} |^2
  // = | (r^{i+1}_{n+1} - r^i_{n+1} |^2

  // initialise increment vector with solution of last iteration (i)
  // update del_ with current residual vector
  // difference of last two solutions
  if (del_ == Teuchos::null)
  {
    del_ = LINALG::CreateVector(*(FluidField()->DofRowMap(0)), true);
    delhist_ = LINALG::CreateVector(*(FluidField()->DofRowMap(0)), true);
    del_->PutScalar(1.0e20);
    delhist_->PutScalar(0.0);
  }

  // calculate difference of current (i+1) and old (i) residual vector
  // del = r^{i+1}_{n+1}
  del_->Update(1.0, *fluidincnp_, 0.0);
  // delhist = ( r^{i+1}_{n+1} - r^i_{n+1} )
  delhist_->Update(1.0, *del_, (-1.0));
  double normdel = 0.0;
  double dot = 0.0;
  delhist_->Norm2(&normdel);
  // calculate dot product
  // dot = delhist_ . del_ = ( r^{i+1}_{n+1} - r^i_{n+1} )^T . r^{i+1}_{n+1}
  del_->Dot(*delhist_, &dot);

  // Aikten factor
  // omega^{i+1}_{n+1} == omeganp_
  // omega^{i}_{n+1} == omegan_
  // ome^{i+1} = ome^i + (ome^i -1) . (r^i - r^{i+1})^T . r^{i+1} / |r^{i+1} - r^{i}|^2
  omeganp_ = omegan_ + (omegan_ - 1.0) * (-dot) / (normdel * normdel);

  // relaxation parameter
  // omega_relax^{i+1} = 1- omega^{i+1}
  double relax = 1.0 - omeganp_;

  if (Comm().MyPID() == 0)
    std::cout << "Aitken relaxation with omega_relax = " << relax << std::endl;

  // relax displacement solution for next iteration step
  // overwrite temp_ with relaxed solution vector
  // T^{i+1} = relax . T^{i+1} + (1- relax^{i+1}) T^i
  //         = T^i + relax^{i+1} * ( T^{i+1} - T^i )
  fluidveln_->Update(relax, *del_, 1.0);

  // update Aitken parameter omega^{i+1}_{n+1}
  omegan_ = omeganp_;

  // update history vector with residual displacement of old iteration step
  delhist_->Update(1.0, *del_, 0.0);

  // end Aitken relaxation
  // ------------------------------------------------------------
}

/*----------------------------------------------------------------------*
                                                           vuong 01/12  |
*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> POROELAST::Partitioned::DofRowMapStructure()
{
  return StructureField()->DofRowMap();
}

/*----------------------------------------------------------------------*
                                                           vuong 01/12  |
*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> POROELAST::Partitioned::DofRowMapFluid()
{
  return FluidField()->DofRowMap();
}
