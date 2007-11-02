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
#ifdef D_FLUID

#include <ctime>
#include <cstdlib>
#include <iostream>
#include <set>
#include <Teuchos_TimeMonitor.hpp>

#ifdef PARALLEL
#include <mpi.h>
#endif

#include "xfluid_dyn_nln_drt.H"
#include "fluidimplicitintegration.H"
#include "fluid_genalpha_integration.H"
#include "../drt_lib/drt_resulttest.H"
#include "fluidresulttest.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_xfem/intersection.H"
#include "../drt_xfem/integrationcell.H"
#include "../drt_xfem/dof_management.H"
#include "../drt_xfem/gmsh.H"

/*----------------------------------------------------------------------*
  |                                                       m.gee 06/01    |
  | vector of numfld FIELDs, defined in global_control.c                 |
 *----------------------------------------------------------------------*/
extern struct _FIELD      *field;

/*----------------------------------------------------------------------*
  |                                                       m.gee 06/01    |
  | general problem data                                                 |
  | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;

/*!----------------------------------------------------------------------
\brief file pointers

<pre>                                                         m.gee 8/00
This structure struct _FILES allfiles is defined in input_control_global.c
and the type is in standardtypes.h
It holds all file pointers and some variables needed for the FRSYSTEM
</pre>
*----------------------------------------------------------------------*/
extern struct _FILES  allfiles;

/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | structure of flags to control output                                 |
 | defined in out_global.c                                              |
 *----------------------------------------------------------------------*/
extern struct _IO_FLAGS     ioflags;

/*----------------------------------------------------------------------*
 | global variable *solv, vector of lenght numfld of structures SOLVAR  |
 | defined in solver_control.c                                          |
 |                                                                      |
 |                                                       m.gee 11/00    |
 *----------------------------------------------------------------------*/
extern struct _SOLVAR  *solv;

/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | pointer to allocate dynamic variables if needed                      |
 | dedfined in global_control.c                                         |
 | ALLDYNA               *alldyn;                                       |
 *----------------------------------------------------------------------*/
extern ALLDYNA      *alldyn;

/*----------------------------------------------------------------------*
 |                                                       m.gee 02/02    |
 | number of load curves numcurve                                       |
 | vector of structures of curves                                       |
 | defined in input_curves.c                                            |
 | INT                   numcurve;                                      |
 | struct _CURVE      *curve;                                           |
 *----------------------------------------------------------------------*/
extern INT            numcurve;
extern struct _CURVE *curve;




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
  map<int, XFEM::DomainIntCells > elementDomainIntCellsMap;
  is.computeIntersection(fluiddis,soliddis,elementDomainIntCellsMap);

  // debug: write both meshes to file in Gmsh format
  ofstream f_system;
  f_system.open ("elements_coupled_system.pos");
  f_system << GMSH::disToString("Fluid", 0.0, fluiddis, elementDomainIntCellsMap);
  f_system << GMSH::disToString("Solid", 1.0, soliddis);
  f_system << GMSH::getConfigString(2);
  f_system.close();

  // apply enrichments
  map<int, const set <XFEM::EnrPhysVar> > nodalDofMap = XFEM::createNodalDofMap(fluiddis, elementDomainIntCellsMap);
  const XFEM::Enrichment enr_std(0, XFEM::Enrichment::typeStandard);
  
  // debug: print enrichments to screen
  for (int i=0; i<fluiddis->NumMyRowNodes(); ++i)
  {
    const int gid = fluiddis->lRowNode(i)->Id();
    const set <XFEM::EnrPhysVar> actset = nodalDofMap[gid];
    for ( set<XFEM::EnrPhysVar>::const_iterator var = actset.begin(); var != actset.end(); ++var )
    {
      cout << "Node: " << gid << ", " << var->toString() << endl;
    };
  };

  
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
    FluidImplicitTimeInt::SetDefaults(fluidtimeparams);

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
    FluidImplicitTimeInt fluidimplicit(fluiddis,
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
    testmanager.AddFieldTest(rcp(new FluidResultTest(fluidimplicit)));
    testmanager.TestAll();
#endif
  }
  else if (fdyn->iop == timeint_gen_alpha)
  {
    // -------------------------------------------------------------------
    // create a generalised alpha time integrator for fluid problems
    // -------------------------------------------------------------------
    // ------------------ set the parameter list
    ParameterList fluidtimeparams;

    fluidtimeparams.set<int>              ("number of velocity degrees of freedom" ,genprob.ndim);
    fluidtimeparams.set<double>           ("time step size"           ,fdyn->dt);
    fluidtimeparams.set<double>           ("total time"               ,fdyn->maxtime);
    fluidtimeparams.set<double>           ("alpha_M"                  ,fdyn->alpha_m);
    fluidtimeparams.set<double>           ("alpha_F"                  ,fdyn->alpha_f);
    fluidtimeparams.set<int>              ("max number timesteps"     ,fdyn->nstep);

    // ---------------------------------------------- nonlinear iteration
    //     // set linearisation scheme
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

    //------------evaluate error for test flows with analytical solutions
    fluidtimeparams.set                  ("eval err for analyt sol"   ,fdyn->init);

    //------------compute statistical data for turbulent channel LES
    if(fdyn->turbu==4)
    {
      fluidtimeparams.set("normal to hom. planes in channel",fdyn->planenormal);
      fluidtimeparams.set("evaluate turbulence statistic",true);
      fluidtimeparams.set("statistics outfile",allfiles.outputfile_kenner);
    }
    else
    {
      fluidtimeparams.set("evaluate turbulence statistic",false);
    }

    //--------------------------------------------------
    // create all vectors and variables associated with the time
    // integration (call the constructor)
    // the only parameter from the list required here is the number of
    // velocity degrees of freedom
    FluidGenAlphaIntegration genalphaint(fluiddis,
                                         solver,
                                         fluidtimeparams,
                                         output);


    //------------- initialise the field from input or restart
    if (genprob.restart)
    {
      // read the restart information, set vectors and variables
      genalphaint.ReadRestart(genprob.restart);
    }
    else
    {
      // set initial field for analytical test problems etc
      if(fdyn->init>0)
      {
        genalphaint.SetInitialFlowField(fdyn->init,fdyn->startfuncno);
      }
    }

    //------------------------- do timeintegration till maxtime
    genalphaint.GenAlphaIntegrateTo(fdyn->nstep,fdyn->maxtime);

    //--------------------------------------------------
    // do the result test
#ifdef RESULTTEST
    DRT::ResultTestManager testmanager(fluiddis->Comm());
    testmanager.AddFieldTest(rcp(new FluidResultTest(genalphaint)));
    testmanager.TestAll();
#endif

  }
  else
  {
    dserror("Unknown time type for drt fluid");
  }

  return;
} // end of dyn_fluid_drt()

#endif
#endif  // #ifdef CCADISCRET
