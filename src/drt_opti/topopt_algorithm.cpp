/*---------------------------------------------------------------------*/
/*! \file

\brief base algorithm for topology optimization of fluid domains

\maintainer Martin Kronbichler

\level 3

*/
/*---------------------------------------------------------------------*/

#include "topopt_algorithm.H"
#include "topopt_fluidAdjoint_timeint.H"
#include "topopt_optimizer.H"
#include "../drt_fluid/fluid_timint_topopt.H"

#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_globalproblem.H"

#include <Teuchos_TimeMonitor.hpp>


/*------------------------------------------------------------------------------------------------*
 | constructor                                                                    winklmaier 12/11 |
 *------------------------------------------------------------------------------------------------*/
TOPOPT::Algorithm::Algorithm(const Epetra_Comm& comm, const Teuchos::ParameterList& topopt)
    : FluidTopOptCouplingAlgorithm(comm, topopt),
      topopt_(topopt),
      optimizer_(Optimizer()),
      doGradient_(true),
      gradienttype_(
          DRT::INPUT::IntegralValue<INPAR::TOPOPT::GradientType>(topopt_, "GRADIENT_TYPE")),
      restarttype_(INPAR::TOPOPT::no_restart)
{
  return;
}


void TOPOPT::Algorithm::TimeLoop() { dserror("No time loop in main optimization routine!"); }


/*------------------------------------------------------------------------------------------------*
 | public: unused time loop of the main algorithm                                winklmaier 12/11 |
 *------------------------------------------------------------------------------------------------*/
void TOPOPT::Algorithm::OptimizationLoop()
{
  PrepareOptimization();

  // prepare data for new optimization iteration
  PrepareFluidField();

  // solve the primary field
  DoFluidField();

  // stop if only primal equations are tested
  if (DRT::INPUT::IntegralValue<INPAR::TOPOPT::AdjointCase>(
          AlgoParameters().sublist("TOPOLOGY ADJOINT FLUID"), "TESTCASE") ==
      INPAR::TOPOPT::adjointtest_primal)
    return;

  // required when restart is called, otherwise only objective value required
  FinishOptimizationStep();

  // optimization process has not yet finished
  while (OptimizationFinished() == false)
  {
    if (doGradient_)
    {
      switch (gradienttype_)
      {
        case INPAR::TOPOPT::gradientByAdjoints:
        {
          // Transfer data from primary field to adjoint field
          PrepareAdjointField();

          // solve the adjoint equations
          DoAdjointField();

          break;
        }
        case INPAR::TOPOPT::gradientByFD1:
        {
          FDGradient(1);
          break;
        }
        case INPAR::TOPOPT::gradientByFD2:
        {
          FDGradient(2);
          break;
        }
        default:
        {
          dserror("unknown approach for computation of objective gradient");
          break;
        }
      }
    }

    if (DRT::INPUT::IntegralValue<INPAR::TOPOPT::AdjointCase>(
            AlgoParameters().sublist("TOPOLOGY ADJOINT FLUID"), "TESTCASE") !=
        INPAR::TOPOPT::adjointtest_no)
      break;  // stop here if adjoints are tested

    // compute the gradient of the objective function
    PrepareOptimizationStep();

    // update objective due to optimization approach
    DoOptimizationStep();

    // write output of optimization step
    Output();

    // Update data for new optimization step
    Update();

    // prepare data for new optimization iteration
    PrepareFluidField();

    // solve the primary field
    DoFluidField();

    // handle inner optimization stuff for which fluid field is required
    FinishOptimizationStep();
  }
}



/*------------------------------------------------------------------------------------------------*
 | protected: prepare the optimization routine                                   winklmaier 12/11 |
 *------------------------------------------------------------------------------------------------*/
void TOPOPT::Algorithm::PrepareOptimization()
{
  // set parameter required by optimizer, belonging to fluid field(s)
  optimizer_->ImportFlowParams(AdjointFluidField()->AdjointParams());

  return;
}



/*------------------------------------------------------------------------------------------------*
 | protected: check if optimization process has finished                         winklmaier 12/11 |
 *------------------------------------------------------------------------------------------------*/
bool TOPOPT::Algorithm::OptimizationFinished()
{
  INPAR::TOPOPT::AdjointCase testcase = DRT::INPUT::IntegralValue<INPAR::TOPOPT::AdjointCase>(
      AlgoParameters().sublist("TOPOLOGY ADJOINT FLUID"), "TESTCASE");

  if (testcase != INPAR::TOPOPT::adjointtest_no and optimizer_->Iter() >= 1)
    return true;  // special test cases for adjoint equations -> no optimization

  if (FluidField()->Discretization()->Comm().MyPID() == 0)
  {
    printf("\n");
    printf(
        "          +----------------------------------------------------------------------+        "
        "  \n");
    printf(
        "          |                          Optimization Step                           |        "
        "  \n");
    printf(
        "          +----------------------------------------------------------------------+        "
        "  \n");
  }

  return optimizer_->Converged(doGradient_);
}


/*------------------------------------------------------------------------------------------------*
 | protected: prepare one (stationary or non-stationary) fluid solution          winklmaier 12/11 |
 *------------------------------------------------------------------------------------------------*/
void TOPOPT::Algorithm::PrepareFluidField()
{
  Teuchos::rcp_dynamic_cast<FLD::TimIntTopOpt>(FluidField())
      ->SetTopOptData(Optimizer()->Density(), optimizer_);

  // reset initial field
  const Teuchos::ParameterList& fdyn = DRT::Problem::Instance()->FluidDynamicParams();
  INPAR::FLUID::InitialField initfield =
      DRT::INPUT::IntegralValue<INPAR::FLUID::InitialField>(fdyn, "INITIALFIELD");
  if (initfield != INPAR::FLUID::initfield_zero_field)
  {
    int startfuncno = fdyn.get<int>("STARTFUNCNO");
    if (initfield != INPAR::FLUID::initfield_field_by_function and
        initfield != INPAR::FLUID::initfield_disturbed_field_from_function)
      startfuncno = -1;
    FluidField()->SetInitialFlowField(initfield, startfuncno);
  }

  return;
}



/*------------------------------------------------------------------------------------------------*
 | protected: perform one (stationary or non-stationary) fluid solution          winklmaier 12/11 |
 *------------------------------------------------------------------------------------------------*/
void TOPOPT::Algorithm::DoFluidField()
{
  if (FluidField()->Discretization()->Comm().MyPID() == 0)
  {
    printf("\n");
    printf(
        "          +----------------------------------------------------------------------------+  "
        "        \n");
    printf(
        "          |                          Solving Fluid Equations                           |  "
        "        \n");
    printf(
        "          +----------------------------------------------------------------------------+  "
        "        \n");
  }

  FluidField()->Integrate();
  return;
}



/*------------------------------------------------------------------------------------------------*
 | protected: prepare one (stationary or non-stationary) fluid solution          winklmaier 12/11 |
 *------------------------------------------------------------------------------------------------*/
void TOPOPT::Algorithm::PrepareAdjointField()
{
  AdjointFluidField()->SetTopOptData(
      Optimizer()->ExportFluidData(), Optimizer()->Density(), optimizer_);

  // reset initial field
  const Teuchos::ParameterList& adjointfdyn =
      DRT::Problem::Instance()->OptimizationControlParams().sublist("TOPOLOGY ADJOINT FLUID");
  INPAR::TOPOPT::InitialAdjointField initfield =
      DRT::INPUT::IntegralValue<INPAR::TOPOPT::InitialAdjointField>(adjointfdyn, "INITIALFIELD");
  int startfuncno = adjointfdyn.get<int>("INITFUNCNO");
  if (initfield != INPAR::TOPOPT::initadjointfield_field_by_function) startfuncno = -1;
  AdjointFluidField()->SetInitialAdjointField(initfield, startfuncno);

  return;
}



/*------------------------------------------------------------------------------------------------*
 | protected: perform one (stationary or non-stationary) adjoint solution        winklmaier 12/11 |
 *------------------------------------------------------------------------------------------------*/
void TOPOPT::Algorithm::DoAdjointField()
{
  if (AdjointFluidField()->Discretization()->Comm().MyPID() == 0)
  {
    printf("\n");
    printf(
        "          "
        "+------------------------------------------------------------------------------+          "
        "\n");
    printf(
        "          |                          Solving Adjoint Equations                           "
        "|          \n");
    printf(
        "          "
        "+------------------------------------------------------------------------------+          "
        "\n");
  }

  AdjointFluidField()->Integrate();
}



/*------------------------------------------------------------------------------------------------*
 | compute gradient by finite differences                                        winklmaier 06/13 |
 *------------------------------------------------------------------------------------------------*/
void TOPOPT::Algorithm::FDGradient(const int numFDPoints)
{
  Teuchos::RCP<const Epetra_Vector> density = Optimizer()->Density();

  if (density->GlobalLength() > 5000)
    dserror("really that much fluid solutions for gradient computation???");

  double c =
      1.0e-5;  // good step size for gradient approximation (low discretization and round-off error)

  for (int i = 0; i < Optimizer()->RowMap()->NumGlobalElements(); i++)
  {
    if (Optimizer()->OptiDis()->Comm().MyPID() == 0)
    {
      printf("+----------------------------------------------------------------------+\n");
      printf("|- computing finite difference gradient in direction %6d of %6d -|\n", i + 1,
          density->GlobalLength());
      printf("+----------------------------------------------------------------------+\n");
    }

    for (int j = 0; j < numFDPoints; j++)
    {
      if (j == 1) c = -c;  // first direction is +c, second -c

      // change optimization variable in current direction
      Optimizer()->AdoptDensityForFD(c, i);

      // adapt porosity field and clear no more required fluid data
      Update();

      // data of primal (fluid) problem no more required
      FluidField()->Reset(
          true, AlgoParameters().get<int>("NUM_OUTPUT_STEPS"), Optimizer()->Iter() + 2);

      // prepare data for new optimization iteration
      PrepareFluidField();

      // solve the primary field
      DoFluidField();

      // compute gradient in this direction
      Optimizer()->ComputeGradientDirectionForFD(c, i, j, numFDPoints);

      // remove change in order to recover original optimization variable
      Optimizer()->AdoptDensityForFD(-c, i);
    }
  }
}



/*------------------------------------------------------------------------------------------------*
 | protected: evaluate the gradient of the optimization objectives               winklmaier 12/11 |
 *------------------------------------------------------------------------------------------------*/
void TOPOPT::Algorithm::PrepareOptimizationStep()
{
  if (doGradient_ and gradienttype_ == INPAR::TOPOPT::gradientByAdjoints) optimizer_->Gradients();

  // data of primal (fluid) problem no more required
  FluidField()->Reset(true, AlgoParameters().get<int>("NUM_OUTPUT_STEPS"), Optimizer()->Iter() + 2);

  // hack to get new optimization parameters in fluid without new routines...
  FluidField()->Init();

  return;
}



/*------------------------------------------------------------------------------------------------*
 | protected: perform one optimization step for the topology optimization        winklmaier 12/11 |
 *------------------------------------------------------------------------------------------------*/
void TOPOPT::Algorithm::DoOptimizationStep()
{
  optimizer_->Iterate(doGradient_);
  return;
}



/*------------------------------------------------------------------------------------------------*
 | protected: finish one optimization step                                       winklmaier 12/11 |
 *------------------------------------------------------------------------------------------------*/
void TOPOPT::Algorithm::FinishOptimizationStep()
{
  optimizer_->Values();

  optimizer_->FinishIteration(doGradient_);

  // data of dual (adjoint) problem no more required
  if (doGradient_)
    AdjointFluidField()->Reset(
        AlgoParameters().get<int>("NUM_OUTPUT_STEPS"), Optimizer()->Iter() + 1);

  return;
}



/*------------------------------------------------------------------------------------------------*
 | protected: output                                                             winklmaier 12/11 |
 *------------------------------------------------------------------------------------------------*/
void TOPOPT::Algorithm::Output()
{
  // write only output of optimizer here!
  // all other data has been written before by the respective routines
  Optimizer()->Output();
  return;
}



/*------------------------------------------------------------------------------------------------*
 | protected: update after one optimization step                                 winklmaier 12/11 |
 *------------------------------------------------------------------------------------------------*/
void TOPOPT::Algorithm::Update()
{
  // clear the field data of the primal and the dual equations
  optimizer_->ClearFieldData();

  Optimizer()->UpdateOptimizationParameter();

  return;
}



/* -----------------------------------------------------------------------------------------------*
 * Restart topology optimization at a specifiec point                            winklmaier 12/11 |
 * -----------------------------------------------------------------------------------------------*/
void TOPOPT::Algorithm::Restart(const int step, const INPAR::TOPOPT::Restart type)
{
  restarttype_ = type;
  switch (type)
  {
    case INPAR::TOPOPT::adjoint:
    {
      dserror("currently not implemented restart type");
      AdjointFluidField()->ReadRestart(step);
      break;
    }
    case INPAR::TOPOPT::fluid:
    {
      dserror("currently not implemented restart type");
      FluidField()->ReadRestart(step);
      break;
    }
    case INPAR::TOPOPT::gradient:
    {
      dserror("currently not implemented restart type");
      break;
    }
    case INPAR::TOPOPT::opti_step:
    {
      if (gradienttype_ != INPAR::TOPOPT::gradientByAdjoints)
        dserror("restart only implemented for adjoints approach!");

      Optimizer()->ReadRestart(step);
      break;
    }
    default:
    {
      dserror("unknown restart type");
      break;
    }
  }

  return;
}
