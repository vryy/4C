/*!----------------------------------------------------------------------
\file xfluid_dyn_nln_drt.cpp
\brief Control routine for fluid time integration. Includes

     o Singele step one-step-theta time integration

     o Two step BDF2 Gear's methode with one-step-theta start step



<pre>
Maintainer: Axel Gerstenberger
            gerstenberger@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include <ctime>
#include <cstdlib>
#include <iostream>
#include <set>
#include <Teuchos_TimeMonitor.hpp>

#ifdef PARALLEL
#include <mpi.h>
#endif

#include "xfluid_dyn_nln_drt.H"
#include "xfluidimplicitintegration.H"
#include "../drt_lib/drt_resulttest.H"
#include "xfluidresulttest.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_xfem/intersection.H"
#include "../drt_xfem/integrationcell.H"
#include "../drt_xfem/dof_management.H"
#include "../drt_xfem/gmsh.H"

extern struct _FIELD    *field;
extern struct _GENPROB  genprob;
extern struct _FILES  	allfiles;
extern struct _IO_FLAGS ioflags;
extern struct _SOLVAR   *solv;
extern         ALLDYNA  *alldyn;
extern struct _CURVE  	*curve;


/*----------------------------------------------------------------------*
 * Time integration loop for fluid.
 *
 *        o One-step-theta
 *        o BDF2
 *
 *----------------------------------------------------------------------*/
void xdyn_fluid_drt()
{
  cout << "Hallo, ich bin ein Fluid-XFEM Problem" << endl;
  flush(cout);
  // -------------------------------------------------------------------
  // access the discretization
  // -------------------------------------------------------------------
  RefCountPtr<DRT::Discretization> fluiddis = null;
  fluiddis = DRT::Problem::Instance()->Dis(genprob.numff,0);
  RefCountPtr<DRT::Discretization> soliddis = null;
  soliddis = DRT::Problem::Instance()->Dis(genprob.numsf,0);
 
  const int fmyrank = fluiddis->Comm().MyPID();
  cout << "FluidProc: " << fmyrank << endl;
  flush(cout);
  const int smyrank = soliddis->Comm().MyPID();
  cout << "SolidProc: " << smyrank << endl;
  flush(cout);
  //cout << *soliddis;
  
  // -------------------------------------------------------------------
  // set degrees of freedom in the discretization
  // -------------------------------------------------------------------
  if (!fluiddis->Filled()) fluiddis->FillComplete();
  if (!soliddis->Filled()) soliddis->FillComplete();


  // -------------------------------------------------------------------
  // context for output and restart
  // -------------------------------------------------------------------
  IO::DiscretizationWriter output(fluiddis);
  output.WriteMesh(0,0.0);

  // Intersection
  XFEM::Intersection is;
  map<int, XFEM::DomainIntCells > elementalDomainIntCells;
  map<int, XFEM::BoundaryIntCells > elementalBoundaryIntCells;
  is.computeIntersection(fluiddis,soliddis,elementalDomainIntCells,elementalBoundaryIntCells);

  // debug: write both meshes to file in Gmsh format
  ofstream f_system;
  f_system.open ("elements_coupled_system.pos");
  f_system << GMSH::disToString("Fluid", 0.0, fluiddis, elementalDomainIntCells);
  f_system << GMSH::disToString("Solid", 1.0, soliddis);
  f_system << GMSH::getConfigString(2);
  f_system.close();

  // apply enrichments
  DofManager dofmanager(fluiddis, elementalDomainIntCells);
  
  // debug: print enrichments to screen
  cout << dofmanager.toString() << endl;

  
  // -------------------------------------------------------------------
  // set some pointers and variables
  // -------------------------------------------------------------------
  SOLVAR        *fluidsolv  = &solv[genprob.numff];
  SOLVAR        *solidsolv  = &solv[genprob.numsf];

  dsassert(genprob.timetyp==time_dynamic, "alldyn not allocated!!!\n For stationary computations, choose TIMETYP Dynamic and switch time integration for fluid to stationary instead");
  FLUID_DYNAMIC *fdyn     = alldyn[genprob.numff].fdyn;
  STRUCT_DYNAMIC *sdyn    = alldyn[genprob.numsf].sdyn;
  
  fdyn->step              =   0;
  fdyn->acttime           = 0.0;
  
  // -------------------------------------------------------------------
  // create a solver
  // -------------------------------------------------------------------
  RefCountPtr<ParameterList> solveparams = rcp(new ParameterList());
  LINALG::Solver solver(solveparams,fluiddis->Comm(),allfiles.out_err);
  solver.TranslateSolverParameters(*solveparams,fluidsolv);
  fluiddis->ComputeNullSpaceIfNecessary(*solveparams);

  if(fdyn->iop == timeint_stationary
     ||
     fdyn->iop == timeint_one_step_theta
     ||
     fdyn->iop == timeint_bdf2
    )
  {
    // -------------------------------------------------------------------
    // create a fluid nonlinear time integrator
    // -------------------------------------------------------------------
    ParameterList fluidtimeparams;
    XFluidImplicitTimeInt::SetDefaults(fluidtimeparams);

    fluidtimeparams.set<int>              ("number of velocity degrees of freedom" ,genprob.ndim);
    fluidtimeparams.set<double>           ("time step size"           ,fdyn->dt);
    fluidtimeparams.set<double>           ("total time"               ,fdyn->maxtime);
    fluidtimeparams.set<double>           ("theta"                    ,fdyn->theta);
    fluidtimeparams.set<FLUID_TIMEINTTYPE>("time int algo"            ,fdyn->iop);
    fluidtimeparams.set<int>              ("max number timesteps"     ,fdyn->nstep);
    fluidtimeparams.set<int>              ("number of start steps"    ,fdyn->nums);
    fluidtimeparams.set<double>           ("start theta"              ,fdyn->thetas);


    // ---------------------------------------------- nonlinear iteration
    // set linearisation scheme
    if(fdyn->ite==2)
    {
      fluidtimeparams.set<bool>("Use reaction terms for linearisation",true);
    }
    else
    {
      fluidtimeparams.set<bool>("Use reaction terms for linearisation",false);
    }
    fluidtimeparams.set<int>             ("max nonlin iter steps"     ,fdyn->itemax);
    // stop nonlinear iteration when both incr-norms are below this bound
    fluidtimeparams.set<double>          ("tolerance for nonlin iter" ,fdyn->ittol);

    // ----------------------------------------------- restart and output
    fluidtimeparams.set                  ("write restart every"       ,fdyn->uprestart);
    fluidtimeparams.set                  ("write solution every"      ,fdyn->upres);
    fluidtimeparams.set                  ("write stresses"            ,ioflags.fluid_stress);

    //--------------------------------------------------
    fluidtimeparams.set                  ("eval err for analyt sol"   ,fdyn->init);


    //--------------------------------------------------
    // create all vectors and variables associated with the time
    // integration (call the constructor)
    // the only parameter from the list required here is the number of
    // velocity degrees of freedom
    XFluidImplicitTimeInt fluidimplicit(fluiddis,
                                       solver,
                                       fluidtimeparams,
                                       output);

    //--------------------------------------------------
    if (genprob.restart)
    {
      // read the restart information, set vectors and variables
      fluidimplicit.ReadRestart(genprob.restart);
    }
    else
    {
      // set initial field for analytical test problems etc
      if(fdyn->init>0)
      {
        fluidimplicit.SetInitialFlowField(fdyn->init,fdyn->startfuncno);
      }
    }

    //--------------------------------------------------
    // do the time integration (start algo and standard algo)
    fluidimplicit.Integrate();



    //--------------------------------------------------
    // do the result test
#ifdef RESULTTEST
    DRT::ResultTestManager testmanager(fluiddis->Comm());
    testmanager.AddFieldTest(rcp(new XFluidResultTest(fluidimplicit)));
    testmanager.TestAll();
#endif
  }
  else
  {
    dserror("Unknown time type for drt xfluid");
  }

  return;

} // end of xdyn_fluid_drt()

#endif  // #ifdef CCADISCRET
