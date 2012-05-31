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
  /* TODO description (winklmaier)
  */
TOPOPT::Algorithm::Algorithm(
    const Epetra_Comm& comm,
    const Teuchos::ParameterList& topopt
)
: FluidTopOptCouplingAlgorithm(comm, topopt),
  topopt_(topopt),
  myrank_(comm.MyPID()),
  iter_(0),
  max_iter_(topopt_.sublist("TOPOLOGY OPTIMIZER").get<int>("MAX_ITER")),
  res_tol_(topopt_.get<double>("RESTOL")),
  inc_tol_(topopt_.get<double>("INCTOL")),
  conv_check_type_(DRT::INPUT::IntegralValue<INPAR::TOPOPT::ConvCheck>(topopt_,"CONV_CHECK_TYPE")),
  optimizer_(Optimizer())
{
  // initialize system vector (without values)
  poro_ = rcp(new Epetra_Vector(*optimizer_->OptiDis()->NodeColMap(),false));

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

  // optimization process has not yet finished
  while (OptimizationNotFinished())
  {
    // prepare data for new optimization iteration
    PrepareFluidField();

    // solve the primary field
    DoFluidField();

    // Transfer data from primary field to adjoint field
    PrepareAdjointField();

    // solve the adjoint equations
    DoAdjointField();

    // compute the gradient of the objective function
    PrepareOptimizationStep();
//
//    // update objective due to optimization approach
//    DoOptimizationStep();

//    write output of optimization step
//    Output();

    // Update data for new optimization step
    Update();
  }


  return;
}



/*------------------------------------------------------------------------------------------------*
 | protected: prepare the optimization routine                                   winklmaier 12/11 |
 *------------------------------------------------------------------------------------------------*/
void TOPOPT::Algorithm::PrepareOptimization()
{
  optimizer_->ImportFlowParams(AdjointFluidField()->AdjointParams());
  return;
}



/*------------------------------------------------------------------------------------------------*
 | protected: check if optimization process has finished                         winklmaier 12/11 |
 *------------------------------------------------------------------------------------------------*/
bool TOPOPT::Algorithm::OptimizationNotFinished()
{
  /*
   * Usually it notfinished is initially set true and it is checked if
   * convergence is reached
   *
   * Since here in optimization usually only the objective function is checked
   * and the density increment usually not, it is done here the other way:
   * setting notfinished false and checking whether convergence is not reached
   *
   * winklmaier 12/11
   */
  iter_++;
  if (max_iter_<1) dserror("maximal number of optimization steps smaller than 1");

  bool notFinished = false;

  if (iter_==1)
  {
    notFinished = true;
    // TODO initial output
  }
  else
  {
    if (iter_>max_iter_)
    {
      notFinished = false;
      // TODO output
    }
    else // TODO check convergence also if maxiter reached?
    {
      if ((conv_check_type_==INPAR::TOPOPT::inc) or
          (conv_check_type_==INPAR::TOPOPT::inc_and_res))
      {
        Epetra_Vector inc(*optimizer_->RowMap(),false);
        inc.Update(1.0,*optimizer_->DensityIp(),-1.0,*optimizer_->DensityI(),0.0);

        double incvelnorm;
        inc.Norm2(&incvelnorm);

        if (incvelnorm>inc_tol_)
          notFinished = true;
      }

      if ((conv_check_type_==INPAR::TOPOPT::res) or
          (conv_check_type_==INPAR::TOPOPT::inc_and_res))
      {
        if (fabs(objective_ip_-objective_i_)>res_tol_)
          notFinished = true;
      }
      // TODO output
    }
  }


  if (notFinished)
  {
    UpdatePorosity();
  }

  return notFinished;
}


/*------------------------------------------------------------------------------------------------*
 | protected: prepare one (stationary or non-stationary) fluid solution          winklmaier 12/11 |
 *------------------------------------------------------------------------------------------------*/
void TOPOPT::Algorithm::PrepareFluidField()
{
  FluidField().SetTopOptData(poro_,optimizer_);
  return;
}



/*------------------------------------------------------------------------------------------------*
 | protected: perform one (stationary or non-stationary) fluid solution          winklmaier 12/11 |
 *------------------------------------------------------------------------------------------------*/
void TOPOPT::Algorithm::DoFluidField()
{
  FluidField().Integrate();
  return;
}



/*------------------------------------------------------------------------------------------------*
 | protected: prepare one (stationary or non-stationary) fluid solution          winklmaier 12/11 |
 *------------------------------------------------------------------------------------------------*/
void TOPOPT::Algorithm::PrepareAdjointField()
{
  AdjointFluidField()->SetTopOptData(optimizer_->ExportFluidData(),poro_,optimizer_);
  return;
}



/*------------------------------------------------------------------------------------------------*
 | protected: perform one (stationary or non-stationary) adjoint solution        winklmaier 12/11 |
 *------------------------------------------------------------------------------------------------*/
void TOPOPT::Algorithm::DoAdjointField()
{
  AdjointFluidField()->Integrate();
}



/*------------------------------------------------------------------------------------------------*
 | protected: evaluate the gradient of the optimization objectives               winklmaier 12/11 |
 *------------------------------------------------------------------------------------------------*/
void TOPOPT::Algorithm::PrepareOptimizationStep()
{
  optimizer_->ComputeObjective();
  optimizer_->ComputeGradient();
  return;
}



/*------------------------------------------------------------------------------------------------*
 | protected: perform one optimization step for the topology optimization        winklmaier 12/11 |
 *------------------------------------------------------------------------------------------------*/
void TOPOPT::Algorithm::DoOptimizationStep()
{
  dserror("global topology optimization has no time loop");
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
  RCP<const Epetra_Vector> density = optimizer_->DensityIp();

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
void TOPOPT::Algorithm::Output() const
{
  dserror("This function is currently unused!\n"
      "Output of subfields is written by their own output-functions!\n"
      "Maybe, this will be used for optimization variable(s) or for overview");
  return;
}



/*------------------------------------------------------------------------------------------------*
 | protected: update after one optimization step                                 winklmaier 12/11 |
 *------------------------------------------------------------------------------------------------*/
void TOPOPT::Algorithm::Update()
{
  // clear the field data of the primal and the dual equations
  optimizer_->ClearFieldData();
  return;
}



/* -----------------------------------------------------------------------------------------------*
 * Restart topology optimization at a specifiec point                            winklmaier 12/11 |
 * -----------------------------------------------------------------------------------------------*/
void TOPOPT::Algorithm::Restart(int step, const int type)
{dserror("not implemented");
  return;
}
