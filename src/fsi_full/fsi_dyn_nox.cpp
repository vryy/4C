
#ifdef TRILINOS_PACKAGE
#ifdef D_FSI

#ifdef PARALLEL
#include <mpi.h>
#endif

#include "fsi_dyn_nox.H"
#include "fsi_nox_aitken.H"
#include "fsi_nox_fixpoint.H"
#include "../discret/dstrc.H"

#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <algorithm>

#include <EpetraExt_MapColoring.h>
#include <EpetraExt_MapColoringIndex.h>
#include <Epetra_CrsGraph.h>
#include <Epetra_IntVector.h>
#include <Epetra_MapColoring.h>
#include <Epetra_SerialComm.h>
#include <NOX_Abstract_Group.H>
#include <NOX_Epetra_BroydenOperator.H>
#include <NOX_Epetra_FiniteDifferenceColoring.H>
#include <NOX_Epetra_Vector.H>
#include <NOX_Solver_Manager.H>

#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_Time.hpp>

#include "fsi_nox_interface_helper.H"


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
extern Teuchos::RefCountPtr<Teuchos::ParameterList> globalparameterlist;


/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | structure of flags to control output                                 |
 | defined in out_global.c                                              |
 *----------------------------------------------------------------------*/
extern struct _IO_FLAGS     ioflags;

/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | pointer to allocate dynamic variables if needed                      |
 | dedfined in global_control.c                                         |
 | ALLDYNA               *alldyn;                                       |
 *----------------------------------------------------------------------*/
extern ALLDYNA      *alldyn;

/*!----------------------------------------------------------------------
\brief file pointers

<pre>                                                         m.gee 8/00
This structure struct _FILES allfiles is defined in input_control_global.c
and the type is in standardtypes.h
It holds all file pointers and some variables needed for the FRSYSTEM
</pre>
*----------------------------------------------------------------------*/
extern struct _FILES  allfiles;

/*!----------------------------------------------------------------------
\brief ranks and communicators

<pre>                                                         m.gee 8/00
This structure struct _PAR par; is defined in main_ccarat.c
and the type is in partition.h
</pre>

*----------------------------------------------------------------------*/
extern struct _PAR   par;


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
extern "C" void debug_out_data(FIELD *actfield, CHAR* n, NODE_ARRAY array, INT pos)
{
  static INT call = 0;
  INT i;
  FILE* f;
  CHAR name[50];

  if (getenv("DEBUG")==NULL || par.myrank!=0)
    return;

  sprintf(name, "plot/%s%d.plot", n, call);
  printf("write '" YELLOW_LIGHT "%s" END_COLOR "'\n", name);

  f = fopen(name, "w");

  for (i=0; i<actfield->dis[0].numele; ++i)
  {
    INT j;
    ELEMENT* actele = &(actfield->dis[0].element[i]);
    fprintf(f, "# element %d\n", actele->Id);
    for (j=0; j<=actele->numnp; ++j)
    {
      INT k;
      NODE* actnode = actele->node[j%actele->numnp];
      fprintf(f, "%e %e ",
              actnode->x[0], actnode->x[1]);
      for (k=0; k<actnode->numdf; ++k)
      {
        switch (array)
        {
        case node_array_sol:
          fprintf(f, "%20.20e ",
                  actnode->sol.a.da[pos][k]);
          break;
        case node_array_sol_increment:
          fprintf(f, "%20.20e ",
                  actnode->sol_increment.a.da[pos][k]);
          break;
        case node_array_sol_mf:
          fprintf(f, "%20.20e ",
                  actnode->sol_mf.a.da[pos][k]);
          break;
        default:
          dserror("node array %d unsupported", array);
        }
      }
      fprintf(f, "\t\t# node %d\n", actnode->Id);
    }
    fprintf(f, "\n\n");
  }

  fclose(f);
  call++;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FSI_InterfaceProblem::FSI_InterfaceProblem(Epetra_Comm& Comm,
                                           FSI_STRUCT_WORK* struct_work,
                                           FSI_FLUID_WORK* fluid_work,
                                           FSI_ALE_WORK* ale_work)
  : struct_work_(struct_work),
    fluid_work_(fluid_work),
    ale_work_(ale_work),
    Comm_(Comm),
    counter_(7)
{
  fsidyn_ = alldyn[3].fsidyn;

  // We need to use the normal structural dof numbers for the global
  // ids of the interface vector. That's ok, trilinos lets us use any
  // set of distinct values for global ids.

  FIELD* structfield = struct_work->struct_field;
  int disnum = 0;
  int numnp_total = structfield->dis[disnum].numnp;
  int numaf       = genprob.numaf;

  std::vector<int> gid;
  std::vector<int> redundant_gid;

  // Find all interface dofs.

  for (int i=0;i<numnp_total;i++)
  {
    NODE* actsnode  = &(structfield->dis[disnum].node[i]);
    GNODE* actsgnode = actsnode->gnode;

    NODE* actanode  = actsgnode->mfcpnode[numaf];
    if (actanode == NULL) continue;
    GNODE* actagnode = actanode->gnode;

    /* check for coupling nodes */
    if (actagnode->dirich == NULL)
      dserror("no dirich condition for coupled ALE node #%d",actanode->Id);

    if (actagnode->dirich->dirich_type != DIRICH_CONDITION::dirich_FSI) continue;

    for (int j=0;j<actanode->numdf;j++)
    {
      if (actsnode->proc==Comm.MyPID())
      {
        gid.push_back(actsnode->dof[j]);
      }
      redundant_gid.push_back(actsnode->dof[j]);
    }
  }

  // We need to write the results to all interface nodes on all processors.
  redundantmap_ = Teuchos::rcp(new Epetra_Map(fsidyn_->numsid,
                                              redundant_gid.size(), &redundant_gid[0],
                                              0, Comm));
  redundantsol_ = Teuchos::rcp(new Epetra_Vector(*redundantmap_));
  redundantf_   = Teuchos::rcp(new Epetra_Vector(*redundantmap_));

  StandardMap_ = Teuchos::rcp(new Epetra_Map(fsidyn_->numsid,
                                             gid.size(), &gid[0],
                                             0, Comm));
  soln_ = Teuchos::rcp(new Epetra_Vector(*StandardMap_));

  // We need to make the solution vector redundant before we can put
  // the results to the nodes. We have an unique source map, so use
  // the importer here.
  if (!StandardMap_->UniqueGIDs())
    dserror("source map not unique");
  importredundant_ = Teuchos::rcp(new Epetra_Import(*redundantmap_,*StandardMap_));

  // create connection graph of interface elements
  rawGraph_ = Teuchos::rcp(new Epetra_CrsGraph(Copy,*StandardMap_,12));

  int numele_total = structfield->dis[disnum].numele;
  for (int i=0; i<numele_total; ++i)
  {
    ELEMENT* actele = &(structfield->dis[disnum].element[i]);
    for (int j=0; j<actele->numnp; ++j)
    {
      NODE* actsnode = actele->node[j];
      GNODE* actsgnode = actsnode->gnode;

      if (actsnode->proc!=Comm.MyPID())
        continue;

      NODE* actanode  = actsgnode->mfcpnode[numaf];
      if (actanode != NULL)
      {
        GNODE* actagnode = actanode->gnode;
        if ((actagnode->dirich != NULL) &&
            (actagnode->dirich->dirich_type == DIRICH_CONDITION::dirich_FSI))
        {
          // So this is a node at the FSI interface. Couple it with
          // all the nodes of this element that are at the interface
          // as well.

          for (int l=0; l<actele->numnp; ++l)
          {
            NODE* snode = actele->node[l];
            NODE* anode = snode->gnode->mfcpnode[numaf];
            if ((anode != NULL) &&
                (anode->gnode->dirich != NULL) &&
                (anode->gnode->dirich->dirich_type == DIRICH_CONDITION::dirich_FSI))
            {
              for (int dof1=0; dof1<actanode->numdf; ++dof1)
              {
                // The ale node known the number of dofs to be
                // coupled, but the dofs of the structural node are
                // actually used.
                int err = rawGraph_->InsertGlobalIndices(actsnode->dof[dof1],anode->numdf,snode->dof);
                if (err != 0)
                  dserror("Epetra_CrsGraph::InsertGlobalIndices returned %d", err);
              }
            }
          }
        }
      }
    }
  }
  rawGraph_->FillComplete();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool FSI_InterfaceProblem::computeF(const Epetra_Vector& x,
                                    Epetra_Vector& F,
                                    const FillType fillFlag)
{
  char* flags[] = { "Residual", "Jac", "Prec", "FD_Res", "MF_Res", "MF_Jac", "User", NULL };

  if (Comm_.MyPID()==0)
    cout << "\n==================================================================================================\n"
         << "FSI_InterfaceProblem::computeF: fillFlag = " RED << flags[fillFlag] << END_COLOR "\n\n";

  // we count the number of times the residuum is build
  counter_[fillFlag] += 1;

  if (!x.Map().UniqueGIDs())
    dserror("source map not unique");

  FIELD* structfield = struct_work_->struct_field;
  FIELD* fluidfield  = fluid_work_ ->fluid_field;
  FIELD* alefield    = ale_work_   ->ale_field;

  int disnum = 0;

  // subdivision is not supported.
  INT a_disnum_calc = 0;
  INT a_disnum_io   = 0;
  INT f_disnum_calc = 0;
  INT f_disnum_io   = 0;
  INT s_disnum_calc = 0;
  INT s_disnum_io   = 0;

  // set new interface displacement x
  ARRAY_POSITION *ipos = &(structfield->dis[disnum].ipos);

  // no consecutive iterations here. We have to care about the fields
  // backup from the outside.
  INT itnum = 1;

  // restore structfield state
  solserv_sol_copy(structfield, s_disnum_calc, node_array_sol_mf, node_array_sol_mf, 8,ipos->mf_dispnp);
  solserv_sol_copy(structfield, s_disnum_calc, node_array_sol_mf, node_array_sol_mf, 9,ipos->mf_reldisp);
  if (itnum > 0)
    solserv_sol_copy(structfield,s_disnum_calc,node_array_sol,node_array_sol,11, 9);
  if (itnum == 0)
    solserv_sol_copy(structfield,s_disnum_calc,node_array_sol,node_array_sol,12, 1);

#ifdef DEBUG
  if (getenv("DEBUG")!=NULL && Comm_.NumProc()==1)
  {
    static int in_counter;
    std::ostringstream filename;
    filename << "plot/x_" << in_counter << ".plot";
    std::cout << "write '" YELLOW_LIGHT << filename.str() << END_COLOR "'\n";
    std::ofstream out(filename.str().c_str());
    for (int i=0; i<=60; i+=2)
    {
      out << i << " " << x[i] << " " << x[i+1] << "\n";
    }
    in_counter += 1;
  }
#endif

  redundantsol_->PutScalar(0);

  int err = redundantsol_->Import(x,*importredundant_,Insert);
  if (err!=0)
    dserror("Import failed with err=%d",err);

  DistributeDisplacements dd(*redundantsol_,ipos->mf_dispnp);
  loop_interface(structfield,dd,*redundantsol_);

  // Calculate new interface displacements starting from the given
  // ones.

  /*------------------------------- CMD -------------------------------*/
  perf_begin(44);
  fsi_ale_calc(ale_work_,alefield,a_disnum_calc,a_disnum_io,structfield,s_disnum_calc);
  perf_begin(44);

//   debug_out_data(alefield, "ale_disp_incr", node_array_sol_increment, alefield->dis[a_disnum_calc].ipos.dispnp);
  debug_out_data(alefield, "ale_disp", node_array_sol_mf, alefield->dis[a_disnum_calc].ipos.mf_dispnp);

  /*------------------------------- CFD -------------------------------*/
  perf_begin(42);
  fsi_fluid_calc(fluid_work_,fluidfield,f_disnum_calc,f_disnum_io,alefield,a_disnum_calc);
  perf_end(42);

  debug_out_data(fluidfield, "fluid_vel", node_array_sol_increment, fluidfield->dis[0].ipos.velnp);

  /*------------------------------- CSD -------------------------------*/
  perf_begin(43);
  fsi_struct_calc(struct_work_,structfield,s_disnum_calc,s_disnum_io,itnum,fluidfield,f_disnum_calc);
  perf_end(43);

  debug_out_data(structfield, "struct_disp", node_array_sol, 0);

  // Fill the distributed interface displacement vector F. We
  // don't have a direct mapping, so loop all structural nodes and
  // select the dofs that belong to the local part of the interface.

  GatherResiduum gr(F,ipos->mf_dispnp,ipos->mf_reldisp);
  loop_interface(structfield,gr,F);

#ifdef DEBUG
  if (getenv("DEBUG")!=NULL && Comm_.NumProc()==1)
  {
    static int out_counter;
    std::ostringstream filename;
    filename << "plot/interface_" << out_counter << ".plot";
    std::cout << "write '" YELLOW_LIGHT << filename.str() << END_COLOR "'\n";
    std::ofstream out(filename.str().c_str());
    OutputDisplacements od(out,ipos->mf_dispnp,ipos->mf_reldisp);
    loop_interface(structfield,od,F);
    out_counter += 1;
  }
#endif

  return true;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RefCountPtr<Epetra_Vector> FSI_InterfaceProblem::soln()
{
  FIELD *structfield = struct_work_->struct_field;

  int disnum = 0;
  ARRAY_POSITION *ipos = &(structfield->dis[disnum].ipos);

  GatherDisplacements gd(*soln_,ipos->mf_dispnp);
  loop_interface(structfield,gd,*soln_);

  return soln_;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI_InterfaceProblem::timeloop(const Teuchos::RefCountPtr<NOX::Epetra::Interface::Required>& interface)
{
  int resstep=0;                /* counter for output control  */
  int restartstep=0;            /* counter for restart control */
  FILE *out = allfiles.out_out;

  FIELD* structfield = struct_work_->struct_field;
  FIELD* fluidfield  = fluid_work_ ->fluid_field;
  FIELD* alefield    = ale_work_   ->ale_field;

  int disnum = 0;
  ARRAY_POSITION *ipos = &(structfield->dis[disnum].ipos);

  // subdivision is not supported.
  INT a_disnum_calc = 0;
  INT a_disnum_io   = 0;
  INT f_disnum_calc = 0;
  INT f_disnum_io   = 0;
  INT s_disnum_calc = 0;
  INT s_disnum_io   = 0;

  FLUID_DYNAMIC  *fdyn;
  FSI_DYNAMIC    *fsidyn;
  STRUCT_DYNAMIC *sdyn;
  ALE_DYNAMIC    *adyn;

  sdyn= alldyn[0].sdyn;
  fdyn= alldyn[1].fdyn;
  adyn= alldyn[2].adyn;
  fsidyn= alldyn[3].fsidyn;

  /* print the mesh to GiD */
  if (ioflags.output_gid==1 && par.myrank==0)
  {
    out_gid_msh();
    out_gid_sol_fsi(fluidfield,structfield,f_disnum_io,s_disnum_io);
  }

  /*
   * Binary output has to be done by the algorithms because the
   * contexts are there. */
  if (ioflags.output_bin==1)
  {
    fsi_ale_output(ale_work_,alefield,a_disnum_calc,a_disnum_io);
    fsi_fluid_output(fluid_work_,fluidfield,f_disnum_calc,f_disnum_io);
    fsi_struct_output(struct_work_,structfield,s_disnum_calc,s_disnum_io);
  }

  /* write general data to .out */
  if (par.myrank==0)
  {
    fprintf(out,"max. values:\n");
    fprintf(out,"============\n");

    /* table head */
    fprintf(out," time |            |field|fluid| fluid error in ");

    switch(fdyn->itnorm)
    {
    case FLUID_DYNAMIC::fncc_Linf: /* infinity norm */
      fprintf(out,"inf-norm");
      break;
    case FLUID_DYNAMIC::fncc_L1: /* L_1 norm */
      fprintf(out,"L_1-norm");
      break;
    case FLUID_DYNAMIC::fncc_L2: /* L_2 norm */
      fprintf(out,"L_2-norm");
      break;
    default:
      dserror("Norm for nonlin. convergence check unknown!!\n");
    }

    fprintf(out," |struc| convergence| relaxation |   total    |\n");

    fprintf(out," step |  sim. time | ite | ite |     vel.   |     pre.   | ite | over fields|  parameter | calc. time |\n");
    fprintf(out,"-------------------------------------------------------------------------------------------------------\n");

    /* max values */
    fprintf(out,"%5d | %10.3f | %3d | %3d |        %10.3E       | %3d | %10.3E |            |            |\n",
            fdyn->nstep,
            fdyn->maxtime,
            fsidyn->itemax,
            fdyn->itemax,
            fdyn->ittol,
            sdyn->maxiter,
            fsidyn->convtol);
    fprintf(out,"-------------------------------------------------------------------------------------------------------\n");

    fprintf(out,"\n\ntimeloop:  ");
    fprintf(out,"Iterative Staggered Scheme based on NOX\n");
    fprintf(out,"=========\n");



    /* table head */
    fprintf(out," time |            |field|fluid| fluid error in ");

    switch(fdyn->itnorm)
    {
      case FLUID_DYNAMIC::fncc_Linf: /* infinity norm */
        fprintf(out,"inf-norm");
        break;
      case FLUID_DYNAMIC::fncc_L1: /* L_1 norm */
        fprintf(out,"L_1-norm");
        break;
      case FLUID_DYNAMIC::fncc_L2: /* L_2 norm */
        fprintf(out,"L_2-norm");
        break;
      default:
        dserror("Norm for nonlin. convergence check unknown!!\n");
    }

    fprintf(out," |struc| convergence| relaxation |   total    |\n");

    fprintf(out," step |  sim. time | ite | ite |     vel.   |     pre.   | ite | over fields|  parameter | calc. time |\n");
    fprintf(out,"-------------------------------------------------------------------------------------------------------\n");
  }

  fflush(out);

  // Get the top level parameter list
  Teuchos::ParameterList& nlParams = *globalparameterlist;

  // sublists

  //Teuchos::ParameterList& searchParams = nlParams.sublist("Line Search");
  Teuchos::ParameterList& printParams = nlParams.sublist("Printing");

  Teuchos::ParameterList& dirParams = nlParams.sublist("Direction");
  Teuchos::ParameterList& newtonParams = dirParams.sublist(dirParams.get("Method","Newton"));
  Teuchos::ParameterList& lsParams = newtonParams.sublist("Linear Solver");

  printParams.set("MyPID", par.myrank);

  // Create printing utilities
  Teuchos::RefCountPtr<NOX::Utils> utils = Teuchos::rcp(new NOX::Utils(printParams));

  // Set user defined aitken line search object.
  if (nlParams.sublist("Line Search").get("Method","Full Step")=="Aitken")
  {
    // insert user defined aitken relaxation
    Teuchos::ParameterList& linesearch = nlParams.sublist("Line Search");
    Teuchos::RefCountPtr<NOX::LineSearch::Generic> aitken = Teuchos::rcp(new AitkenRelaxation(utils,linesearch));

    // We change the method here.
    linesearch.set("Method","User Defined");
    linesearch.set("User Defined Line Search",aitken);
  }

  // ==================================================================

  // log solver iterations

  Teuchos::RefCountPtr<std::ofstream> log;
  if (par.myrank==0)
  {
    std::string s = allfiles.outputfile_kenner;
    s.append(".iteration");
    log = Teuchos::rcp(new std::ofstream(s.c_str()));
    (*log) << "# num procs      = " << par.nprocs << "\n"
           << "# Method         = " << dirParams.get("Method","Newton") << "\n"
           << "# Jacobian       = " << nlParams.get("Jacobian", "Matrix Free") << "\n"
           << "# Preconditioner = " << nlParams.get("Preconditioner","None") << "\n"
           << "# Line Search    = " << nlParams.sublist("Line Search").get("Method","Full Step") << "\n"
           << "#\n"
           << "# step  time/step  #nliter  #liter  Residual  Jac  Prec  FD_Res  MF_Res  MF_Jac  User\n"
      ;
  }

  Teuchos::Time timer("time step timer");

  // Create the Epetra_RowMatrix using Finite Difference with Coloring
  // Here we take the graph of the dof connections at the FSI
  // interface and create a coloring with the EpetraExt helpers. This
  // is needed to setup the Finite Difference with Coloring class
  // later on.

  EpetraExt::CrsGraph_MapColoring::ColoringAlgorithm algType = EpetraExt::CrsGraph_MapColoring::GREEDY;
  EpetraExt::CrsGraph_MapColoring tmpMapColoring(algType);
  Teuchos::RefCountPtr<Epetra_MapColoring> colorMap = Teuchos::rcp(&tmpMapColoring(*rawGraph_));
  EpetraExt::CrsGraph_MapColoringIndex colorMapIndex(*colorMap);
  Teuchos::RefCountPtr< vector<Epetra_IntVector> > columns = Teuchos::rcp(&colorMapIndex(*rawGraph_));

  // A distance 1 coloring can be used to build just the diagonal of
  // the Jacobian. This is much faster and might be a suitable
  // preconditioner.

  EpetraExt::CrsGraph_MapColoring distance1MapColoring(algType,0,true);
  Teuchos::RefCountPtr<Epetra_MapColoring> distance1ColorMap = Teuchos::rcp(&distance1MapColoring(*rawGraph_));
  EpetraExt::CrsGraph_MapColoringIndex distance1ColorMapIndex(*distance1ColorMap);
  Teuchos::RefCountPtr< vector<Epetra_IntVector> > distance1Columns = Teuchos::rcp(&distance1ColorMapIndex(*rawGraph_));

  while (fsidyn->step < fsidyn->nstep && fsidyn->time <= fsidyn->maxtime)
  {
    double t2=ds_cputime();

    fsidyn->step++;
    fsidyn->time += fsidyn->dt;

    fdyn->step=fsidyn->step;
    sdyn->step=fsidyn->step;
    adyn->step=fsidyn->step;

    fdyn->acttime = fsidyn->time;
    sdyn->time = fsidyn->time;
    adyn->time = fsidyn->time;

    if (par.myrank==0)
      fsi_algoout(0);

    // reset all counters
    std::fill(counter_.begin(),counter_.end(),0);
    lsParams.sublist("Output").set("Total Number of Linear Iterations",0);

    // start time measurement
    Teuchos::RefCountPtr<Teuchos::TimeMonitor> timemonitor = rcp(new Teuchos::TimeMonitor(timer,true));

    // aus dem Strukturl�er:

      /*======================================================================*
       *                      S O L U T I O N    P H A S E                    *
       *======================================================================*
       * nodal solution history structural field:                             *
       * sol[0][j]           ... total displacements at time (t)              *
       * sol[1][j]           ... velocities at time (t)                       *
       * sol[2][j]           ... accels at time (t)                           *
       * sol[3][j]           ... prescribed displacements at time (t-dt)      *
       * sol[4][j]           ... prescribed displacements at time (t)         *
       * sol[5][j]           ... place 4 - place 3                            *
       * sol[6][j]           ... the  velocities of prescribed dofs           *
       * sol[7][j]           ... the  accels of prescribed dofs               *
       * sol[8][j]           ... working space                                *
       * sol[9][j]           ... total displacements at time (t-dt)           *
       * sol[10][j]          ... velocities at time (t-dt)                    *
       * sol_mf[0][j]        ... latest struct-displacements                  *
       * sol_mf[1][j]        ... (relaxed) displ. of the last iteration step  *
       * sol_mf[2][j]        ... converged relaxed displ. at time (t-dt)      *
       * sol_mf[3][j]        ... actual dispi                                 *
       * sol_mf[4][j]        ... FSI coupl.-forces at the end of the timestep *
       * sol_mf[5][j]        ... FSI coupl.-forces at beginning of the timest.*
       * sol_mf[6][j]        ... used in fsi_gradient.c                       *
       *======================================================================*/

    // no consecutive iterations here. We have to care about the fields
    // backup from the outside.
    INT itnum = 1;

    // backup current structural values
    //solserv_sol_copy(structfield,s_disnum_calc,node_array_sol,node_array_sol,0,9);
    //solserv_sol_copy(structfield,s_disnum_calc,node_array_sol,node_array_sol,1,10);
    /* backup for Newton via finite differences */
    {
      ARRAY_POSITION *ipos;
      ipos = &(structfield->dis[s_disnum_calc].ipos);
      solserv_sol_copy(structfield, s_disnum_calc, node_array_sol_mf, node_array_sol_mf, ipos->mf_dispnp, 8);
      solserv_sol_copy(structfield, s_disnum_calc, node_array_sol_mf, node_array_sol_mf, ipos->mf_reldisp, 9);
      if (itnum > 0)
        solserv_sol_copy(structfield,s_disnum_calc,node_array_sol,node_array_sol,9,11);
      if (itnum == 0)
        solserv_sol_copy(structfield,s_disnum_calc,node_array_sol,node_array_sol,1,12);
    }

    // Begin Nonlinear Solver ************************************

    NOX::Epetra::Vector noxSoln(soln(), NOX::Epetra::Vector::CreateView);

#ifdef DEBUG
    if (getenv("DEBUG")!=NULL && Comm_.NumProc()==1)
    {
      static int out_counter;
      std::ostringstream filename;
      filename << "plot/noxSoln_" << out_counter << ".plot";
      std::cout << "write '" YELLOW_LIGHT << filename.str() << END_COLOR "'\n";
      std::ofstream out(filename.str().c_str());
      //noxSoln.print(out);
      const Epetra_Vector& s = noxSoln.getEpetraVector();
      for (int i=0; i<=60; i+=2)
      {
        out << i << " " << s[i] << " " << s[i+1] << "\n";
      }
      out_counter += 1;
    }
#endif


    // Create the Epetra_RowMatrix.
    // 1. User supplied (Epetra_RowMatrix)
    //Teuchos::RefCountPtr<Epetra_RowMatrix> Analytic = Problem.getJacobian();
    // 2. Matrix-Free (Epetra_Operator)
    //Teuchos::RefCountPtr<NOX::Epetra::MatrixFree> MF = Teuchos::rcp(new NOX::Epetra::MatrixFree(printParams, interface, noxSoln));
    // 3. Finite Difference (Epetra_RowMatrix)
    //Teuchos::RefCountPtr<NOX::Epetra::FiniteDifference> FD = Teuchos::rcp(new NOX::Epetra::FiniteDifference(printParams, interface, noxSoln));

    // Create the linear system
    Teuchos::RefCountPtr<NOX::Epetra::Interface::Required> iReq = interface;

    // Ok. The variables.

    Teuchos::RefCountPtr<NOX::Epetra::MatrixFree> MF;
    Teuchos::RefCountPtr<NOX::Epetra::FiniteDifferenceColoring> FDC;
    Teuchos::RefCountPtr<NOX::Epetra::FiniteDifferenceColoring> FDC1;
    Teuchos::RefCountPtr<NOX::Epetra::BroydenOperator> B;

    Teuchos::RefCountPtr<NOX::Epetra::Interface::Jacobian> iJac;
    Teuchos::RefCountPtr<NOX::Epetra::Interface::Preconditioner> iPrec;

    Teuchos::RefCountPtr<Epetra_Operator> J;
    Teuchos::RefCountPtr<Epetra_Operator> M;

    Teuchos::RefCountPtr<NOX::Epetra::LinearSystemAztecOO> linSys;

    // ==================================================================
    // decide on Jacobian and preconditioner
    // We migh want to use no preconditioner at all. Some kind of
    // Jacobian has to be provided, otherwise the linear system uses
    // plain finite differences.

    std::string jacobian = nlParams.get("Jacobian", "MatrixFree");
    std::string preconditioner = nlParams.get("Preconditioner", "None");

    // Matrix Free Newton Krylov. This is supposed to be the most
    // appropiate choice.
    if (jacobian=="Matrix Free")
    {
      // MatrixFree seems to be the most interessting choice. But you
      // must set a rather low tolerance for the linear solver.

      MF = Teuchos::rcp(new NOX::Epetra::MatrixFree(printParams, interface, noxSoln));
      //MF->setLambda(1e-2);
      iJac = MF;
      J = MF;
    }

    // No Jacobian at all. Do a fix point iteration. This is a user
    // extension, so we have to modify the parameter list here.
    else if (jacobian=="None")
    {
      Teuchos::RefCountPtr<NOX::Direction::Generic> fixpoint = Teuchos::rcp(new FixPoint(utils,nlParams));
      dirParams.set("Method","User Defined");
      dirParams.set("User Defined Direction",fixpoint);
      lsParams.set("Preconditioner","None");
      preconditioner="None";
    }

    // Finite Difference with coloring. Build a Jacobian as good as it
    // gets. This is the closest we get to the true Jacobian.
    else if (jacobian=="Finite Difference")
    {
      Teuchos::ParameterList& fdParams = nlParams.sublist("Finite Difference");
      double alpha = fdParams.get("alpha", 1.0e-4);
      double beta  = fdParams.get("beta",  1.0e-6);
      std::string dt = fdParams.get("Difference Type","Forward");
      NOX::Epetra::FiniteDifferenceColoring::DifferenceType dtype = NOX::Epetra::FiniteDifferenceColoring::Forward;
      if (dt=="Forward")
        dtype = NOX::Epetra::FiniteDifferenceColoring::Forward;
      else if (dt=="Backward")
        dtype = NOX::Epetra::FiniteDifferenceColoring::Backward;
      else if (dt=="Centered")
        dtype = NOX::Epetra::FiniteDifferenceColoring::Centered;
      else
        dserror("unsupported difference type '%s'",dt.c_str());

      FDC = Teuchos::rcp(new NOX::Epetra::FiniteDifferenceColoring(printParams, interface, noxSoln, rawGraph_, colorMap, columns, true, false, beta, alpha));
      FDC->setDifferenceMethod(dtype);

      iJac = FDC;
      J = FDC;
    }

    // Finite Difference with distance 1 coloring. Build just a
    // diagonal Jacobian.
    else if (jacobian=="Finite Difference 1")
    {
      Teuchos::ParameterList& fdParams = nlParams.sublist("Finite Difference");
      double alpha = fdParams.get("alpha", 1.0e-4);
      double beta  = fdParams.get("beta",  1.0e-6);
      std::string dt = fdParams.get("Difference Type","Forward");
      NOX::Epetra::FiniteDifferenceColoring::DifferenceType dtype = NOX::Epetra::FiniteDifferenceColoring::Forward;
      if (dt=="Forward")
        dtype = NOX::Epetra::FiniteDifferenceColoring::Forward;
      else if (dt=="Backward")
        dtype = NOX::Epetra::FiniteDifferenceColoring::Backward;
      else if (dt=="Centered")
        dtype = NOX::Epetra::FiniteDifferenceColoring::Centered;
      else
        dserror("unsupported difference type '%s'",dt.c_str());

      FDC1 = Teuchos::rcp(new NOX::Epetra::FiniteDifferenceColoring(printParams, interface, noxSoln, rawGraph_, distance1ColorMap, distance1Columns, true, true, beta, alpha));
      FDC1->setDifferenceMethod(dtype);

      iJac = FDC1;
      J = FDC1;
    }

    // Broyden. A strange idea that starts with a real Jacobian and
    // does modify it along the steps of the nonlinear solve.
    else if (jacobian=="Broyden")
    {
#if 0
      // we need a Jacobian to start with.
      if (is_null(FDC))
      {
        FDC = Teuchos::rcp(new NOX::Epetra::FiniteDifferenceColoring(printParams, interface, noxSoln, rawGraph_, colorMap, columns, true));
      }
      FDC->computeJacobian(*soln());

      Teuchos::RefCountPtr<Epetra_CrsMatrix> mat = Teuchos::rcp(new Epetra_CrsMatrix(FDC->getUnderlyingMatrix()));
      //B = Teuchos::rcp(new NOX::Epetra::BroydenOperator(nlParams, *soln(), mat));
      B = Teuchos::rcp(new NOX::Epetra::BroydenOperator(nlParams, utils, *soln(), mat));

      iJac = B;
      J = B;
#else
      dserror("Broyden Operator still unfinished");
#endif
    }
    else
    {
      dserror("unsupported Jacobian '%s'",jacobian.c_str());
    }

    // ==================================================================

    // No preconditioning at all. This might work. But on large
    // systems it probably won't.
    if (preconditioner=="None")
    {
      if (lsParams.get("Preconditioner", "None")!="None")
      {
        if (par.myrank==0)
          utils->out() << RED "Warning: Preconditioner turned on in linear solver settings.\n"
                       << "Jacobian operator will be used for preconditioning as well." END_COLOR "\n";
      }

      if (Teuchos::is_null(iJac))
      {
        // if no Jacobian has been set this better be the fix point
        // method.
        if (dirParams.get("Method","Newton")!="User Defined")
        {
          if (par.myrank==0)
            utils->out() << RED "Warning: No Jacobian for solver " << dirParams.get("Method","Newton") << END_COLOR "\n";
        }
        linSys = Teuchos::rcp(new NOX::Epetra::LinearSystemAztecOO(printParams, lsParams, iReq, noxSoln));
      }
      else
      {
        linSys = Teuchos::rcp(new NOX::Epetra::LinearSystemAztecOO(printParams, lsParams, iReq, iJac, J, noxSoln));
      }
    }

    // Finite Difference with coloring. The best (and most expensive)
    // Jacobian we can get. It might do good on really hard problems.
    else if (preconditioner=="Finite Difference")
    {
      if (lsParams.get("Preconditioner", "None")=="None")
      {
        if (par.myrank==0)
          utils->out() << RED "Warning: Preconditioner turned off in linear solver settings." END_COLOR "\n";
      }

      // A real (approximated) matrix for preconditioning
#if 0
      if (is_null(FDC))
      {
        // FiniteDifferenceColoring might be a good preconditioner for
        // really hard problems. But you should construct it only once
        // a time step. Probably.
        FDC = Teuchos::rcp(new NOX::Epetra::FiniteDifferenceColoring(printParams, interface, noxSoln, rawGraph_, colorMap, columns, true));
      }
      iPrec = FDC;
      M = FDC;
#else
      Teuchos::RefCountPtr<NOX::Epetra::FiniteDifferenceColoring> precFDC =
        Teuchos::rcp(new NOX::Epetra::FiniteDifferenceColoring(printParams, interface, noxSoln, rawGraph_, colorMap, columns, true));
      iPrec = precFDC;
      M = precFDC;
#endif

      linSys = Teuchos::rcp(new NOX::Epetra::LinearSystemAztecOO(printParams, lsParams, iJac, J, iPrec, M, noxSoln));
    }

    // Finite Difference with distance 1 coloring. Supposed to be a
    // suitable (and cheap) preconditioner
    else if (preconditioner=="Finite Difference 1")
    {
      if (lsParams.get("Preconditioner", "None")=="None")
      {
        if (par.myrank==0)
          utils->out() << RED "Warning: Preconditioner truned off in linear solver settings." END_COLOR "\n";
      }
#if 0
      if (is_null(FDC1))
      {
        FDC1 = Teuchos::rcp(new NOX::Epetra::FiniteDifferenceColoring(printParams, interface, noxSoln, rawGraph_, distance1ColorMap, distance1Columns, true, true));
      }
      iPrec = FDC1;
      M = FDC1;
#else
      Teuchos::RefCountPtr<NOX::Epetra::FiniteDifferenceColoring> precFDC =
        Teuchos::rcp(new NOX::Epetra::FiniteDifferenceColoring(printParams, interface, noxSoln, rawGraph_, distance1ColorMap, distance1Columns, true, true));
      iPrec = precFDC;
      M = precFDC;
#endif

      linSys = Teuchos::rcp(new NOX::Epetra::LinearSystemAztecOO(printParams, lsParams, iJac, J, iPrec, M, noxSoln));
    }
    else
    {
      dserror("unsupported preconditioner '%s'",preconditioner.c_str());
    }

    // ==================================================================
    // Convergence Tests

    // Create the Group
    Teuchos::RefCountPtr<NOX::Epetra::Group> grp = Teuchos::rcp(new NOX::Epetra::Group(printParams, iReq, noxSoln, linSys));

    // Create the convergence tests
    Teuchos::RefCountPtr<NOX::StatusTest::Combo> combo       = Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR));
    Teuchos::RefCountPtr<NOX::StatusTest::Combo> converged   = Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::AND));

    Teuchos::RefCountPtr<NOX::StatusTest::MaxIters> maxiters = Teuchos::rcp(new NOX::StatusTest::MaxIters(nlParams.get("Max Iterations", 20)));
    Teuchos::RefCountPtr<NOX::StatusTest::FiniteValue> fv    = Teuchos::rcp(new NOX::StatusTest::FiniteValue);

    combo->addStatusTest(fv);
    combo->addStatusTest(converged);
    combo->addStatusTest(maxiters);

    Teuchos::RefCountPtr<NOX::StatusTest::NormF> absresid = Teuchos::rcp(new NOX::StatusTest::NormF(nlParams.get("Norm abs F", 1.0e-6)));
    converged->addStatusTest(absresid);

    if (nlParams.isParameter("Norm Update"))
    {
      Teuchos::RefCountPtr<NOX::StatusTest::NormUpdate> update = Teuchos::rcp(new NOX::StatusTest::NormUpdate(nlParams.get("Norm Update", 1.0e-5)));
      converged->addStatusTest(update);
    }

    if (nlParams.isParameter("Norm rel F"))
    {
      Teuchos::RefCountPtr<NOX::StatusTest::NormF> relresid = Teuchos::rcp(new NOX::StatusTest::NormF(*grp.get(), nlParams.get("Norm rel F", 1.0e-2)));
      converged->addStatusTest(relresid);
    }

    //Teuchos::RefCountPtr<NOX::StatusTest::NormWRMS> wrms     = Teuchos::rcp(new NOX::StatusTest::NormWRMS(1.0e-2, 1.0e-8));
    //converged->addStatusTest(wrms);

    //////////////////////////////////////////////////////////////////

    // Create the method
    NOX::Solver::Manager solver(grp, combo, globalparameterlist);
    NOX::StatusTest::StatusType status = solver.solve();

    if (status != NOX::StatusTest::Converged)
      if (par.myrank==0)
        utils->out() << "Nonlinear solver failed to converge!" << endl;

    // Get the Epetra_Vector with the final solution from the solver
    const NOX::Epetra::Group& finalGroup = dynamic_cast<const NOX::Epetra::Group&>(solver.getSolutionGroup());
    const Epetra_Vector& finalSolution = (dynamic_cast<const NOX::Epetra::Vector&>(finalGroup.getX())).getEpetraVector();
    const Epetra_Vector& finalF        = (dynamic_cast<const NOX::Epetra::Vector&>(finalGroup.getF())).getEpetraVector();

    // End Nonlinear Solver **************************************

    // Output the parameter list
    if (utils->isPrintType(NOX::Utils::Parameters))
      if (par.myrank==0)
      {
        utils->out() << endl
                     << "Final Parameters" << endl
                     << "****************" << endl;
        solver.getList().print(utils->out());
        utils->out() << endl;
      }

#ifdef DEBUG
    if (getenv("DEBUG")!=NULL && Comm_.NumProc()==1)
    {
      static int out_counter;
      std::ostringstream filename;
      filename << "plot/finalSolution_" << out_counter << ".plot";
      std::cout << "write '" YELLOW_LIGHT << filename.str() << END_COLOR "'\n";
      std::ofstream out(filename.str().c_str());
      for (int i=0; i<=60; i+=2)
      {
        out << i << " " << finalSolution[i]+finalF[i] << " " << finalSolution[i+1]+finalF[i+1] << "\n";
      }
      //finalGroup.getX().print(out);
      out_counter += 1;
    }
#endif

    // ==================================================================
    // return results

    redundantsol_->PutScalar(0);
    int err = redundantsol_->Import(finalSolution,*importredundant_,Insert);
    if (err!=0)
      dserror("Import failed with err=%d",err);

    // The last F has not yet been added to the solution. This F has
    // satisfied the norms, so it is very small. But never the less...
#if 0
    // We can add it here because the source is known to be unique.
    err = redundantsol_->Import(finalF,*importredundant_,Add);
    if (err!=0)
      dserror("Import failed with err=%d",err);
#endif
    redundantf_->PutScalar(0);
    err = redundantf_->Import(finalF,*importredundant_,Insert);
    if (err!=0)
      dserror("Import failed with err=%d",err);

    //DistributeDisplacements dd(*redundantsol_,ipos->mf_dispnp);
    DistributeSolution dd(*redundantsol_,*redundantf_,ipos->mf_dispnp);
    loop_interface(structfield,dd,*redundantsol_);

    // ==================================================================

    // stop time measurement
    timemonitor = Teuchos::null;

    if (par.myrank==0)
    {
      (*log) << fsidyn->step
             << " " << timer.totalElapsedTime()
             << " " << nlParams.sublist("Output").get("Nonlinear Iterations",0)
             << " " << nlParams.sublist("Output").get("2-Norm of Residual", 0.)
             << " " << lsParams.sublist("Output").get("Total Number of Linear Iterations",0)
        ;
      for (std::vector<int>::size_type i=0; i<counter_.size(); ++i)
      {
        (*log) << " " << counter_[i];
      }
      (*log) << std::endl;
      log->flush();
    }

    // ==================================================================

    /*--------------------- update MESH data -------------------------*/
    perf_begin(44);
    fsi_ale_final(ale_work_,alefield,a_disnum_calc,a_disnum_io);
    perf_end(44);

    /*-------------------- update FLUID data -------------------------*/
    perf_begin(42);
    fsi_fluid_final(fluid_work_,fluidfield,f_disnum_calc,f_disnum_io);
    perf_end(42);

    /*------------------ update STRUCTURE data -----------------------*/
    perf_begin(43);
    fsi_struct_final(struct_work_,structfield,s_disnum_calc,s_disnum_io);
    perf_end(43);


    double tt=ds_cputime()-t2;
    if (par.myrank==0)
    {
      fprintf(out," %10.3f |\n",tt);
      fprintf(out,"-------------------------------------------------------------------------------------------------------\n");
    }
    fflush(out);

    /* write current solution */
    resstep++;
    restartstep++;

    if (resstep == fsidyn->upres)
    {
      resstep=0;

      /* print out solution to GiD */
      if (ioflags.output_gid==1 && par.myrank==0)
      {
        out_checkfilesize(1);
        out_gid_sol_fsi(fluidfield,structfield,f_disnum_io,s_disnum_io);
      }


#ifdef BINIO
      /*
       * Binary output has to be done by the algorithms because the
       * contexts are there. */
      if (ioflags.output_bin==1)
      {
        fsi_ale_output(ale_work_,alefield,a_disnum_calc,a_disnum_io);
        fsi_fluid_output(fluid_work_,fluidfield,f_disnum_calc,f_disnum_io);
        fsi_struct_output(struct_work_,structfield,s_disnum_calc,s_disnum_io);
      }
#endif
    }

    /* write fsi-restart data */
    if (restartstep==fsidyn->uprestart)
    {
      restartstep=0;
#ifdef BINIO
      restart_write_bin_fsidyn(fsidyn);
#else
      restart_write_fsidyn(fsidyn);
#endif
    }

    /* energy check */
    if (fsidyn->ichecke>0) fsi_energycheck();
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
extern "C"
void dyn_fsi_nox(FSI_STRUCT_WORK* struct_work,
                 FSI_FLUID_WORK* fluid_work,
                 FSI_ALE_WORK* ale_work)
{
  // Create a communicator for Epetra objects
#ifdef PARALLEL
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  // Create the interface between the test problem and the nonlinear solver
  Teuchos::RefCountPtr<FSI_InterfaceProblem> interface = Teuchos::rcp(new FSI_InterfaceProblem(Comm,struct_work,fluid_work,ale_work));
  interface->timeloop(interface);
}

#endif
#endif
