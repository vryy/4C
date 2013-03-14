/*!------------------------------------------------------------------------------------------------*
\file topopt_algorithm.cpp

\brief base algorithm for topology optimization of fluid domains

<pre>
Maintainer: Martin Winklmaier
            winklmaier@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>
 *------------------------------------------------------------------------------------------------*/


#include "topopt_algorithm.H"
#include "topopt_fluidAdjoint_timeint.H"
#include "topopt_optimizer.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_mat/matpar_bundle.H"
#include "../drt_mat/optimization_density.H"
#include "../linalg/linalg_utils.H"

#include <Teuchos_TimeMonitor.hpp>


/*------------------------------------------------------------------------------------------------*
 | constructor                                                                    winklmaier 12/11 |
 *------------------------------------------------------------------------------------------------*/
TOPOPT::Algorithm::Algorithm(
    const Epetra_Comm& comm,
    const Teuchos::ParameterList& topopt
)
: FluidTopOptCouplingAlgorithm(comm, topopt),
  topopt_(topopt),
  optimizer_(Optimizer()),
  doGradient_(true)
{
  // initialize system vector (without values)
  poro_ = Teuchos::rcp(new Epetra_Vector(*optimizer_->OptiDis()->NodeColMap(),false));

  return;
}


void TOPOPT::Algorithm::TimeLoop()
{
  dserror("No time loop in main optimization routine!");
}


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

  // optimization process has not yet finished
  while (OptimizationFinished() == false)
  {
    if (doGradient_)
    {
      // Transfer data from primary field to adjoint field
      PrepareAdjointField();

      // solve the adjoint equations
      DoAdjointField();
    }

    if (DRT::INPUT::IntegralValue<INPAR::TOPOPT::AdjointCase>(AlgoParameters().sublist("TOPOLOGY ADJOINT FLUID"),"TESTCASE")!=INPAR::TOPOPT::adjointtest_no)
      break; // stop here if adjoints are tested

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

  UpdatePorosity();

  return;
}



/*------------------------------------------------------------------------------------------------*
 | protected: check if optimization process has finished                         winklmaier 12/11 |
 *------------------------------------------------------------------------------------------------*/
bool TOPOPT::Algorithm::OptimizationFinished()
{
  INPAR::TOPOPT::AdjointCase testcase =
      DRT::INPUT::IntegralValue<INPAR::TOPOPT::AdjointCase>(
          AlgoParameters().sublist("TOPOLOGY ADJOINT FLUID"),"TESTCASE");

  if (testcase==INPAR::TOPOPT::adjointtest_primal)
    return true; // test of primal equations
  else if (testcase!=INPAR::TOPOPT::adjointtest_no and optimizer_->Iter()>=1)
    return true; // special test cases for adjoint equations -> no optimization

  if (FluidField().Discretization()->Comm().MyPID()==0)
  {
    printf("\n");
    printf("          +----------------------------------------------------------------------+          \n");
    printf("          |                          Optimization Step                           |          \n");
    printf("          +----------------------------------------------------------------------+          \n");
  }

  return optimizer_->Converged(doGradient_);
}


/*------------------------------------------------------------------------------------------------*
 | protected: prepare one (stationary or non-stationary) fluid solution          winklmaier 12/11 |
 *------------------------------------------------------------------------------------------------*/
void TOPOPT::Algorithm::PrepareFluidField()
{
  FluidField().SetTopOptData(poro_,optimizer_);

  // reset initial field
  const Teuchos::ParameterList& fdyn = DRT::Problem::Instance()->FluidDynamicParams();
  INPAR::FLUID::InitialField initfield = DRT::INPUT::IntegralValue<INPAR::FLUID::InitialField>(fdyn,"INITIALFIELD");
  if(initfield != INPAR::FLUID::initfield_zero_field)
  {
    int startfuncno = fdyn.get<int>("STARTFUNCNO");
    if (initfield != INPAR::FLUID::initfield_field_by_function and
        initfield != INPAR::FLUID::initfield_disturbed_field_from_function)
      startfuncno=-1;
    FluidField().SetInitialFlowField(initfield,startfuncno);
  }

  return;
}



/*------------------------------------------------------------------------------------------------*
 | protected: perform one (stationary or non-stationary) fluid solution          winklmaier 12/11 |
 *------------------------------------------------------------------------------------------------*/
void TOPOPT::Algorithm::DoFluidField()
{
  if (FluidField().Discretization()->Comm().MyPID()==0)
  {
    printf("\n");
    printf("          +----------------------------------------------------------------------------+          \n");
    printf("          |                          Solving Fluid Equations                           |          \n");
    printf("          +----------------------------------------------------------------------------+          \n");
  }

  FluidField().Integrate();
  return;
}



/*------------------------------------------------------------------------------------------------*
 | protected: prepare one (stationary or non-stationary) fluid solution          winklmaier 12/11 |
 *------------------------------------------------------------------------------------------------*/
void TOPOPT::Algorithm::PrepareAdjointField()
{
  AdjointFluidField()->SetTopOptData(optimizer_->ExportFluidData(),poro_,optimizer_);

  // reset initial field
  const Teuchos::ParameterList& adjointfdyn = DRT::Problem::Instance()->OptimizationControlParams().sublist("TOPOLOGY ADJOINT FLUID");
  INPAR::TOPOPT::InitialAdjointField initfield = DRT::INPUT::IntegralValue<INPAR::TOPOPT::InitialAdjointField>(adjointfdyn,"INITIALFIELD");
  int startfuncno = adjointfdyn.get<int>("INITFUNCNO");
  if (initfield != INPAR::TOPOPT::initadjointfield_field_by_function)
    startfuncno=-1;
  AdjointFluidField()->SetInitialAdjointField(initfield,startfuncno);

  return;
}



/*------------------------------------------------------------------------------------------------*
 | protected: perform one (stationary or non-stationary) adjoint solution        winklmaier 12/11 |
 *------------------------------------------------------------------------------------------------*/
void TOPOPT::Algorithm::DoAdjointField()
{
  if (AdjointFluidField()->Discretization()->Comm().MyPID()==0)
  {
    printf("\n");
    printf("          +------------------------------------------------------------------------------+          \n");
    printf("          |                          Solving Adjoint Equations                           |          \n");
    printf("          +------------------------------------------------------------------------------+          \n");
  }

  AdjointFluidField()->Integrate();
}



/*------------------------------------------------------------------------------------------------*
 | protected: evaluate the gradient of the optimization objectives               winklmaier 12/11 |
 *------------------------------------------------------------------------------------------------*/
void TOPOPT::Algorithm::PrepareOptimizationStep()
{
  if (optimizer_->Iter()==0)
    optimizer_->ComputeValues();

  if (doGradient_)
    optimizer_->ComputeGradients();

  FluidField().Reset(
      true,
      DRT::INPUT::IntegralValue<bool>(AlgoParameters(),"OUTPUT_EVERY_ITER"),
      Optimizer()->Iter());

  AdjointFluidField()->Reset(
      DRT::INPUT::IntegralValue<bool>(AlgoParameters(),"OUTPUT_EVERY_ITER"),
      Optimizer()->Iter());

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
  optimizer_->ComputeValues();

  optimizer_->FinishIteration(doGradient_);
  return;
}

/*------------------------------------------------------------------------------------------------*
 | protected: update the porosity field                                          winklmaier 12/11 |
 *------------------------------------------------------------------------------------------------*/
void TOPOPT::Algorithm::UpdatePorosity()
{
  /* The density field is a row vector with own discretization (possibly different
   * than the fluid discretization.
   *
   * Currently the discretizations have to fit. The fluid requires the porosity
   * in dofcolmap. Thus the following steps are required in general:
   *
   * 1) opti-density -> opti-porosity (row maps)
   * 2) opti-porosity -> fluid-porosity (row maps)
   * 3) row fluid-porosity -> col fluid-porosity
   *
   * According to the above description, step 2 is not required currently
   * due to matching discretizations.
   *
   * winklmaier 12/11
   * */
  RCP<const Epetra_Vector> density = optimizer_->Density();

  // get the optimization material
  const MAT::PAR::TopOptDens* mat = NULL;
  const int nummat = DRT::Problem::Instance()->Materials()->Num();
  for (int id = 1; id-1 < nummat; ++id)
  {
    Teuchos::RCP<const MAT::PAR::Material> imat = DRT::Problem::Instance()->Materials()->ById(id);

    if (imat == Teuchos::null)
      dserror("Could not find material Id %d", id);
    else
    {
      if (imat->Type() == INPAR::MAT::m_opti_dens)
      {
        const MAT::PAR::Parameter* matparam = imat->Parameter();
        mat = static_cast<const MAT::PAR::TopOptDens* >(matparam);
        break;
      }
    }
  }
  if (mat==NULL)
    dserror("optimization material not found");

  // and then the material parameters
  const double poro_bd_down = mat->poro_bd_down_;
  const double poro_bd_up = mat->poro_bd_up_;
  const double fac = mat->smear_fac_;

  const double poro_diff = poro_bd_down - poro_bd_up;

  // evaluate porosity due to density
  Epetra_Vector poro_row(density->Map(),false);
  for (int lnodeid=0; lnodeid<optimizer_->OptiDis()->NumMyRowNodes(); lnodeid++)  // loop over processor nodes
  {
    const double dens = (*density)[lnodeid];
    (poro_row)[lnodeid] = poro_bd_up + poro_diff*dens*(1+fac)/(dens+fac);
  }

  // density always in rowmap -> porosity for fluid evaluation in col map required
  LINALG::Export(poro_row,*poro_);
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
  INPAR::TOPOPT::AdjointCase testcase =
      DRT::INPUT::IntegralValue<INPAR::TOPOPT::AdjointCase>(
          AlgoParameters().sublist("TOPOLOGY ADJOINT FLUID"),"TESTCASE");

  if (testcase!=INPAR::TOPOPT::adjointtest_no)
    return; // don't delete results for test cases -> resulttest

  // clear the field data of the primal and the dual equations
  optimizer_->ClearFieldData();

  // update porosity field
  UpdatePorosity();


  return;
}



/* -----------------------------------------------------------------------------------------------*
 * Restart topology optimization at a specifiec point                            winklmaier 12/11 |
 * -----------------------------------------------------------------------------------------------*/
void TOPOPT::Algorithm::Restart(int step, const int type)
{dserror("not implemented");
  return;
}
