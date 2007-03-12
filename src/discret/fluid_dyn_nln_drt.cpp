/*!----------------------------------------------------------------------
\file
\brief Control routine for fluid time integration. Includes

     o Singele step one-step-theta time integration

     o Two step BDF2 Gear's methode with one-step-theta start step


     
<pre>
Maintainer: Peter Gamnitzer
            gamnitzer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE
#ifdef D_FLUID

#include <ctime>
#include <cstdlib>
#include <iostream>

#ifdef PARALLEL
#include <mpi.h>
#endif

#include "fluid_dyn_nln_drt.H"


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



//#include "fluid_dyn_nln_drt.H"


/*----------------------------------------------------------------------*
 * Time integration loop for fluid.
 * 
 *        o One-step-theta
 *        o BDF2
 *
 *----------------------------------------------------------------------*/

void dyn_fluid_drt()
{


 DSTraceHelper dst("dyn_fluid_drt");


  // -------------------------------------------------------------------
  // access the discretization
  // -------------------------------------------------------------------
  RefCountPtr<DRT::Discretization> actdis = null;
  {
    vector<RefCountPtr<DRT::Discretization> >* fool =
              (vector<RefCountPtr<DRT::Discretization> >*)field[0].ccadis;
    actdis = (*fool)[0];
  }
  // set degrees of freedom in the discretization
  if (!actdis->Filled()) actdis->FillComplete();

  
  // -------------------------------------------------------------------
  // get a communicator and myrank
  // -------------------------------------------------------------------
  const Epetra_Comm& Comm = actdis->Comm();
  const int myrank  = Comm.MyPID();

  
  // -------------------------------------------------------------------
  // set some pointers and variables
  // -------------------------------------------------------------------
  SOLVAR        *actsolv  = &solv[0];
  FLUID_DYNAMIC *fdyn     = alldyn[0].fdyn;

  
  //-----------------------------------------------------create a solver
  RefCountPtr<ParameterList> solveparams = rcp(new ParameterList());
  LINALG::Solver solver(solveparams,actdis->Comm(),allfiles.out_err);
  solver.TranslateSolverParameters(*solveparams,actsolv);
  actdis->ComputeNullSpaceIfNecessary(*solveparams);


  // -------------------------------------------------------------------
  // get a vector layout from the discretization to construct matching
  // vectors and matrices
  //                 local <-> global dof numbering
  // -------------------------------------------------------------------
  const Epetra_Map* dofrowmap = actdis->DofRowMap();


  // -------------------------------------------------------------------
  // create empty system matrix --- stiffness and mass are assembled in
  // one system matrix!
  // -------------------------------------------------------------------
  RefCountPtr<Epetra_CrsMatrix> sys_mat
      = LINALG::CreateMatrix(*dofrowmap,81);


  // -------------------------------------------------------------------
  // create empty right hand side vector for linear system
  // -------------------------------------------------------------------
  RefCountPtr<Epetra_Vector> rhs  = LINALG::CreateVector(*dofrowmap,true);

  
  // -------------------------------------------------------------------
  // create empty vectors 
  // -------------------------------------------------------------------

  // Vectors passed to the element
  // -----------------------------
  
  // accelerations at time n and n-1
  RefCountPtr<Epetra_Vector> accn  = LINALG::CreateVector(*dofrowmap,true);
  RefCountPtr<Epetra_Vector> accnm = LINALG::CreateVector(*dofrowmap,true);

  // velocities and pressures at time n+1, n and n-1
  RefCountPtr<Epetra_Vector> velnp = LINALG::CreateVector(*dofrowmap,true);
  RefCountPtr<Epetra_Vector> veln  = LINALG::CreateVector(*dofrowmap,true);
  RefCountPtr<Epetra_Vector> velnm = LINALG::CreateVector(*dofrowmap,true);

  // histvector --- a linear combination of velnm, veln (BDF)
  //                or veln, accn (One-Step-Theta)             
  RefCountPtr<Epetra_Vector> hist  = LINALG::CreateVector(*dofrowmap,true);
  
  // Vectors associated to boundary conditions
  // -----------------------------------------
  
  // toggle vector indicating which dofs have Dirichlet BCs
  RefCountPtr<Epetra_Vector> dirichtoggle = LINALG::CreateVector(*dofrowmap,true);

  // a vector of zeros to be used to enforce zero dirichlet boundary conditions
  RefCountPtr<Epetra_Vector> zeros = LINALG::CreateVector(*dofrowmap,true);


  // the vector containing body and surface forces
  RefCountPtr<Epetra_Vector> neumann_loads = LINALG::CreateVector(*dofrowmap,true);

  
  // Vectors used for solution process 
  // ---------------------------------

  // The residual vector --- more or less the rhs for the incremental
  // formulation!!!
  RefCountPtr<Epetra_Vector> residual = LINALG::CreateVector(*dofrowmap,true);

  // Nonlinear iteration increment vector
  RefCountPtr<Epetra_Vector> incvel = LINALG::CreateVector(*dofrowmap,true);

  /*------------------------------------------- set initial step and time */
  fdyn->step    =   0;
  fdyn->acttime = 0.0;
  
  /*--------------------------------------- init all applied time curves -*/
  for (int actcurve=0; actcurve<numcurve; actcurve++)
  {
   /* the last three parameters are obsolete!!! */  
   dyn_init_curve(actcurve,fdyn->step,fdyn->dt,fdyn->maxtime);
  }

  /*------------------------------------------------------- printout head */
  if (myrank==0)
  {
      
  } /* end if (myrank==0) */


  // -------------------------------------------------------------------
  // context for output and restart
  // -------------------------------------------------------------------
  DiscretizationWriter output(actdis);
  output.WriteMesh(0,0.0);

  //------------------------------------------------- output initial state
  output.NewStep(fdyn->step, fdyn->acttime);
  output.WriteVector("vel_and_pres", velnp);

  
  {
   // save all fluid-dynamic info which will be overwritten by startingalgo
   FLUID_TIMEINTTYPE iop_s    = fdyn->iop;
   double theta_s  = fdyn->theta;
   if (theta_s < EPS5)
   {
     theta_s = 1.0;
   }
   if (iop_s == timeint_bdf2) /* BDF2 */
   {
     theta_s = 1.0;
   }
      
   bool stop_timeloop=false;
   /*<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>*/
   /*<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>*/
   /*<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>*/
   /*                               TIMELOOP                             */
   /*<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>*/
   /*<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>*/
   /*<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>*/
   while (stop_timeloop==false)
   {
    // increase counters and time
    fdyn->step++;
    fdyn->dta  = fdyn->dt;                   /* constant time step size */
    fdyn->acttime +=fdyn->dta;

    // check (starting) algorithm 
    if (fdyn->step<=fdyn->nums) /* set parameter for starting phase */
    {
     fdyn->iop    = timeint_one_step_theta;
     fdyn->theta  = fdyn->thetas;
    }
    else if (fdyn->step==(fdyn->nums+1)) /* set original parameter */
    {
     fdyn->iop    = iop_s;
     fdyn->theta  = theta_s;
    }
    
    // out to screen
    if (myrank==0)
    {
	fluid_algoout();
    }
    
    // set time specific parameters
    switch (fdyn->iop)
    {
	case timeint_one_step_theta: /* One step Theta time integration */
	    fdyn->thsl = fdyn->dta*fdyn->theta;
	    fdyn->thpl = fdyn->thsl;
	    fdyn->thsr = (1.0 - fdyn->theta)*fdyn->dta;
	    fdyn->thpr = fdyn->thsr;
	    break;
	case timeint_bdf2:	/* 2nd order backward differencing BDF2	*/
	    fdyn->thsl = fdyn->dta*2.0/3.0;
	    fdyn->thpl = fdyn->thsl;
	    fdyn->thsr = 0.0;
	    fdyn->thpr = fdyn->thsr;
	    break;
	default:
	    dserror("Time integration scheme unknown for mass rhs!");	    
    }

    
    // part of the residual vector belonging to the old timestep

    /*

         One-step-Theta:

                  vel(n) + dt*(1-Theta)*acc(n)


         BDF2: for constant time step:

                  4/3 vel(n) - 1/3 vel(n-1)
                  
    */

    
    switch (fdyn->iop)
    {
	case timeint_one_step_theta: /* One step Theta time integration */
	    hist->Update(1.0,*veln,fdyn->dta*(1.0-fdyn->theta),*accn,0.0);
	    break;
	case timeint_bdf2:	/* 2nd order backward differencing BDF2	*/
	    hist->Update(4./3.,*veln,-1./3.,*velnm,0.0);
	    break;
	default:
	    dserror("Time integration scheme unknown for mass rhs!");	    
    }

    /* do explicit predictor step to start iteration from better value */
    /*           velnp is still containing veln, the first trial value */
    if(fdyn->step>1)
    {
	velnp->Update(fdyn->dta*(1.0+fdyn->dta/fdyn->dtp),*accn,
		      DSQR(fdyn->dta/fdyn->dtp),*velnm,
		      -DSQR(fdyn->dta/fdyn->dtp));
    }
    
    //-------- set dirichlet boundary conditions 
    //------------------------------------- evaluate Neumann and Dirichlet BCs
    {
     ParameterList params;
     // action for elements
     params.set("action","calc_fluid_eleload");
     // choose what to assemble
     params.set("assemble matrix 1",false);
     params.set("assemble matrix 2",false);
     params.set("assemble vector 1",true);
     params.set("assemble vector 2",false);
     params.set("assemble vector 3",false);
     // other parameters needed by the elements
     params.set("total time",fdyn->acttime);
     params.set("delta time",fdyn->dt);
     // set vector values needed by elements
     actdis->ClearState();
     actdis->SetState("u and p at time n+1 (trial)",velnp);
     // predicted dirichlet values
     // velnp then also holds prescribed new dirichlet values
     actdis->EvaluateDirichlet(params,*velnp,*dirichtoggle);
     actdis->ClearState();

     // evaluate Neumann conditions
     neumann_loads->PutScalar(0.0);
     actdis->EvaluateNeumann(params,*neumann_loads);
     actdis->ClearState();
    }
    
    {
    /*<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>*/
    /*<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>*/
    /*<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>*/
    /*                NONLINEAR ITERATION (FLUID)                       */
    /*<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>*/
    /*<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>*/
    /*<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>*/
     int  n_itnum=0;
     bool stop_nonliniter=false;
     if(myrank == 0)
     {
      printf("-------------------------------------------------\n");
     }
     while (stop_nonliniter==false)
     {
      n_itnum++;
      // -------------------------------------------------------------------
      // call elements to calculate system matrix
      // -------------------------------------------------------------------
      {
       // zero out the stiffness matrix
       sys_mat = LINALG::CreateMatrix(*dofrowmap,81);
       // zero out residual
       residual->PutScalar(0.0);
	  
       // create the parameters for the discretization
       ParameterList params;

       // action for elements
       params.set("action","calc_fluid_systemmat_and_residual");
       // choose what to assemble
       params.set("assemble matrix 1",true);
       params.set("assemble matrix 2",false);
       params.set("assemble vector 1",true);
       params.set("assemble vector 2",false);
       params.set("assemble vector 3",false);
       // other parameters that might be needed by the elements
       params.set("total time",fdyn->acttime);
       params.set("delta time",fdyn->dt);
       // set vector values needed by elements
       actdis->ClearState();
       actdis->SetState("u and p at time n+1 (trial)",velnp);
       actdis->SetState("old solution data for rhs"  ,hist );

       // call loop over elements
       actdis->Evaluate(params,sys_mat,null,residual,null,null);
       actdis->ClearState();

       // finalize the system matrix
       LINALG::Complete(*sys_mat);
      }

      //--------- Apply dirichlet boundary conditions to system of equations
      //          residual discplacements are supposed to be zero at
      //          boundary conditions
      LINALG::ApplyDirichlettoSystem(sys_mat,incvel,residual,
				     zeros,dirichtoggle);
      
      
      
      //-------solve for residual displacements to correct incremental displacements
      incvel->PutScalar(0.0);
      //solver.Solve(sys_mat,incvel,residual,true,false);

      //------------------------------------------------ update (u,p) trial
      velnp->Update(1.0,*incvel,1.0);
      
      // check convergence
      {
       double incvelnorm_L2;
       incvel->Norm2(&incvelnorm_L2);
       
       double velnorm_L2;
       velnp->Norm2(&velnorm_L2);
       
       
       if (velnorm_L2<EPS5)
       {
	velnorm_L2 = 1.0;
       }
       
       if(myrank == 0)
       {
	printf("|  %3d/%3d   | %10.3E[L_2 ]  | %10.3E   |\n",
	       n_itnum,fdyn->itemax,fdyn->ittol,incvelnorm_L2/velnorm_L2);
       }
       
       if(incvelnorm_L2/velnorm_L2<fdyn->ittol
	  ||
	  n_itnum == fdyn->itemax)
       {
	stop_nonliniter=true;
	if(myrank == 0)
	{
	 printf("-------------------------------------------------\n");
	}
       }
      }
     }

     
    /*<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>*/
    /*<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>*/
    /*<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>*/
    /*               END NONLINEAR ITERATION (FLUID)                    */
    /*<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>*/
    /*<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>*/
    /*<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>*/
    }


    //-------------------------------------------- output of solution
    output.NewStep(fdyn->step, fdyn->acttime);
    output.WriteVector("vel_and_pres", velnp);
    
    // update acceleration

    if (fdyn->step == 1)
    {
	// do just a linear interpolation within the first timestep
	accn->Update(1.0/fdyn->dta,*velnp,
		     1.0/fdyn->dta,*veln ,
		     0.0);
	
	// ???
	accnm->Update(1.0,*accn,0.0);
		
    }
    else
    { 
	// prev. acceleration becomes (n-1)-accel. of next time step
	accnm->Update(1.0,*accn,0.0);

	/*

        One-step-Theta:

        acc(n+1) = (vel(n+1)-vel(n)) / (Theta * dt(n)) - (1/Theta -1) * acc(n)


        BDF2:

   	              2*dt(n)+dt(n-1)		     dt(n)+dt(n-1)
        acc(n+1) = --------------------- vel(n+1) - --------------- vel(n)
	           dt(n)*[dt(n)+dt(n-1)]	     dt(n)*dt(n-1)

          		    dt(n)
	         + ----------------------- vel(n+1)
 	           dt(n-1)*[dt(n)+dt(n-1)]

        */

	switch (fdyn->iop)
	{
	    case timeint_one_step_theta: /* One step Theta time integration */
		accn->Update(1.0/(fdyn->theta*fdyn->dta),*velnp,
			     1.0/(fdyn->theta*fdyn->dta),*veln ,
			     (-1.0/fdyn->theta-1));
		break;
	    case timeint_bdf2:	/* 2nd order backward differencing BDF2	*/
	        {
		 double dta = fdyn->dta;
		 double dtp = fdyn->dtp;
		 if (dta*dtp < EPS15)
		     dserror("Zero time step size!!!!!");
		 double sum = dta + dtp;
		 
		 accn->Update((2.0*dta+dtp)/(dta*sum) + dta/(dtp*sum),*velnp,
			      -sum / (dta*dtp),*veln ,
			      0.0);
		}
		break;
	    default:
		dserror("Time integration scheme unknown for mass rhs!");	    
	}
    }

    // solution of this step becomes most recent solution of the last step
    velnm->Update(1.0,*veln ,0.0);
    veln ->Update(1.0,*velnp,0.0);

    // update time step sizes
    fdyn->dtp = fdyn->dta;
    
    
    // check steady state, maxiter and maxtime
    if(fdyn->step==fdyn->nstep||fdyn->acttime>=fdyn->maxtime)
    {
	stop_timeloop=true;
    }
   }
  
  }
  /*<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>*/
  /*<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>*/
  /*<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>*/
  /*                            END TIMELOOP                            */
  /*<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>*/
  /*<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>*/
  /*<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>*/
  
  
  //---------- this is the end. Beautiful friend. My only friend, The end.
  // thanks to RefCountPtr<> we do not need to delete anything here!
  
  return;

} // end of dyn_fluid_drt()

#endif  // #ifdef D_FLUID
#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
