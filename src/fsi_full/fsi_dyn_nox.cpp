
#ifdef TRILINOS_PACKAGE
#ifdef D_FSI

#include "fsi_dyn_nox.H"
#include "fsi_aitken_nox.H"
#include "../discret/dstrc.H"

#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <vector>

#include <NOX_Solver_Manager.H>
#include <NOX_Epetra_Vector.H>
#include <NOX_Epetra_FiniteDifferenceColoring.H>

#include <Teuchos_XMLParameterListHelpers.hpp>
#include <Epetra_SerialComm.h>
#include <EpetraExt_MapColoring.h>
#include <EpetraExt_MapColoringIndex.h>
#include <Epetra_CrsGraph.h>
#include <Epetra_MapColoring.h>
#include <Epetra_IntVector.h>

#include "fsi_nox_interface_helper.H"

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

  if (getenv("DEBUG")==NULL)
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
    ale_work_(ale_work)
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

  // Find all interface dofs.

  for (int i=0;i<numnp_total;i++)
  {
    NODE* actsnode  = &(structfield->dis[disnum].node[i]);
    GNODE* actsgnode = actsnode->gnode;

    NODE* actanode  = actsgnode->mfcpnode[numaf];
    if (actanode == NULL) continue;
    GNODE* actagnode = actanode->gnode;

#ifdef PARALLEL
    dserror("test ownership!");
#endif

    /* check for coupling nodes */
    if (actagnode->dirich == NULL)
      dserror("no dirich condition for coupled ALE node #%d",actanode->Id);

    if (actagnode->dirich->dirich_type != DIRICH_CONDITION::dirich_FSI) continue;

    for (int j=0;j<actanode->numdf;j++)
    {
      gid.push_back(actsnode->dof[j]);
    }
  }

  StandardMap_ = Teuchos::rcp(new Epetra_Map(fsidyn_->numsid,
                                             gid.size(), &gid[0],
                                             0, Comm));
  soln_ = Teuchos::rcp(new Epetra_Vector(*StandardMap_));

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

  cout << "\n==================================================================================================\n"
       << "FSI_InterfaceProblem::computeF: fillFlag = " RED << flags[fillFlag] << END_COLOR "\n\n";

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
  if (getenv("DEBUG")!=NULL)
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

  DistributeDisplacements dd(x,ipos->mf_dispnp);
  loop_interface(structfield,dd,x);

  // Calculate new interface displacements starting from the given
  // ones.

  // reset ale field
  // sollte nicht nötig sein...
//   solserv_sol_copy(alefield,a_disnum_calc,
//                    node_array_sol_increment,
//                    node_array_sol_increment,
//                    alefield->dis[a_disnum_calc].ipos.dispn,
//                    alefield->dis[a_disnum_calc].ipos.dispnp);

//   debug_out_data(alefield, "ale_disp_incr", node_array_sol_increment, alefield->dis[a_disnum_calc].ipos.dispnp);

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
  if (getenv("DEBUG")!=NULL)
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
  DSTraceHelper("FSI_InterfaceProblem::soln");
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
  DSTraceHelper("FSI_InterfaceProblem::timeloop");
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

    fsi_algoout(0);

    // aus dem Strukturlöser:

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

#if 0
#ifdef PARALLEL
    int MyPID = Comm.MyPID();
#else
    int MyPID = 0;
#endif
#endif

    NOX::Epetra::Vector noxSoln(soln(), NOX::Epetra::Vector::CreateView);

#ifdef DEBUG
    if (getenv("DEBUG")!=NULL)
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

    // Create the top level parameter list
    Teuchos::RefCountPtr<Teuchos::ParameterList> nlParamsPtr = Teuchos::rcp(new Teuchos::ParameterList);

    // Read parameters. Hack. This must be merged with our input file.
    Teuchos::updateParametersFromXmlFile("input.xml", nlParamsPtr.get());

    Teuchos::ParameterList& nlParams = *nlParamsPtr.get();

    // sublists

    //Teuchos::ParameterList& searchParams = nlParams.sublist("Line Search");
    Teuchos::ParameterList& printParams = nlParams.sublist("Printing");

    Teuchos::ParameterList& dirParams = nlParams.sublist("Direction");
    Teuchos::ParameterList& newtonParams = dirParams.sublist("Newton");
    Teuchos::ParameterList& lsParams = newtonParams.sublist("Linear Solver");

    //cout << lsParams;

    // Create printing utilities
    NOX::Utils utils(printParams);

#if 0
    // use user defined aitken relaxation
    Teuchos::RefCountPtr<NOX::Direction::Generic> aitken = Teuchos::rcp(new Aitken());
    dirParams.set("Method", "User Defined");
    dirParams.set("User Defined Direction", aitken);
#endif

    // Create the Epetra_RowMatrix.  Uncomment one or more of the following:
    // 1. User supplied (Epetra_RowMatrix)
    //Teuchos::RefCountPtr<Epetra_RowMatrix> Analytic = Problem.getJacobian();
    // 2. Matrix-Free (Epetra_Operator)
    //Teuchos::RefCountPtr<NOX::Epetra::MatrixFree> MF = Teuchos::rcp(new NOX::Epetra::MatrixFree(printParams, interface, noxSoln));
    // 3. Finite Difference (Epetra_RowMatrix)
    //Teuchos::RefCountPtr<NOX::Epetra::FiniteDifference> FD = Teuchos::rcp(new NOX::Epetra::FiniteDifference(printParams, interface, noxSoln));

    Teuchos::RefCountPtr<NOX::Epetra::FiniteDifferenceColoring> FDC =
      Teuchos::rcp(new NOX::Epetra::FiniteDifferenceColoring(printParams, interface, noxSoln, rawGraph_, colorMap, columns));
    Teuchos::RefCountPtr<NOX::Epetra::Interface::Jacobian> iJac_FDC = FDC;

    // Create the linear system
    Teuchos::RefCountPtr<NOX::Epetra::Interface::Required> iReq = interface;
    //Teuchos::RefCountPtr<NOX::Epetra::Interface::Jacobian> iJac = MF;
    //Teuchos::RefCountPtr<NOX::Epetra::Interface::Preconditioner> iPrec = interface;

#if 0
    Teuchos::RefCountPtr<NOX::Epetra::LinearSystemAztecOO> linSys = Teuchos::rcp(new NOX::Epetra::LinearSystemAztecOO(printParams, lsParams,
                                                                                                                      iReq, iJac, MF,
                                                                                                                      noxSoln));
#endif
#if 0
    Teuchos::RefCountPtr<NOX::Epetra::LinearSystemAztecOO> linSys = Teuchos::rcp(new NOX::Epetra::LinearSystemAztecOO(printParams, lsParams,
                                                                                                                      iReq, noxSoln));
#endif
    Teuchos::RefCountPtr<NOX::Epetra::LinearSystemAztecOO> linSys = Teuchos::rcp(new NOX::Epetra::LinearSystemAztecOO(printParams, lsParams,
                                                                                                                      iReq, iJac_FDC, FDC, noxSoln));

    // Create the Group
    Teuchos::RefCountPtr<NOX::Epetra::Group> grp = Teuchos::rcp(new NOX::Epetra::Group(printParams, iReq, noxSoln,
                                                                                       linSys));

    // Create the convergence tests
    Teuchos::RefCountPtr<NOX::StatusTest::NormF> absresid    = Teuchos::rcp(new NOX::StatusTest::NormF(1.0e-6));
    //Teuchos::RefCountPtr<NOX::StatusTest::NormF> relresid    = Teuchos::rcp(new NOX::StatusTest::NormF(*grp.get(), 1.0e-2));
    Teuchos::RefCountPtr<NOX::StatusTest::NormUpdate> update = Teuchos::rcp(new NOX::StatusTest::NormUpdate(1.0e-5));
    Teuchos::RefCountPtr<NOX::StatusTest::NormWRMS> wrms     = Teuchos::rcp(new NOX::StatusTest::NormWRMS(1.0e-2, 1.0e-8));
    Teuchos::RefCountPtr<NOX::StatusTest::Combo> converged   = Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::AND));
    converged->addStatusTest(absresid);
    //converged->addStatusTest(relresid);
    //converged->addStatusTest(wrms);
    //converged->addStatusTest(update);
    Teuchos::RefCountPtr<NOX::StatusTest::MaxIters> maxiters = Teuchos::rcp(new NOX::StatusTest::MaxIters(20));
    Teuchos::RefCountPtr<NOX::StatusTest::FiniteValue> fv    = Teuchos::rcp(new NOX::StatusTest::FiniteValue);
    Teuchos::RefCountPtr<NOX::StatusTest::Combo> combo       = Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR));
    combo->addStatusTest(fv);
    combo->addStatusTest(converged);
    combo->addStatusTest(maxiters);

    // Create the method
    NOX::Solver::Manager solver(grp, combo, nlParamsPtr);
    NOX::StatusTest::StatusType status = solver.solve();

    if (status != NOX::StatusTest::Converged)
      if (par.myrank==0)
        utils.out() << "Nonlinear solver failed to converge!" << endl;

    // Get the Epetra_Vector with the final solution from the solver
    const NOX::Epetra::Group& finalGroup = dynamic_cast<const NOX::Epetra::Group&>(solver.getSolutionGroup());
    const Epetra_Vector& finalSolution = (dynamic_cast<const NOX::Epetra::Vector&>(finalGroup.getX())).getEpetraVector();
    const Epetra_Vector& finalF        = (dynamic_cast<const NOX::Epetra::Vector&>(finalGroup.getF())).getEpetraVector();

    // End Nonlinear Solver **************************************

    // Output the parameter list
    if (utils.isPrintType(NOX::Utils::Parameters))
    {
      utils.out() << endl
                  << "Final Parameters" << endl
                  << "****************" << endl;
      solver.getList().print(utils.out());
      utils.out() << endl;
    }

#ifdef DEBUG
    if (getenv("DEBUG")!=NULL)
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

    DistributeDisplacements dd(finalSolution,ipos->mf_dispnp);
    // Use the last step we calculated.
    // No need to do this. The last but one solution is (by
    // definition) the same as the last one.
    //DistributeSolution dd(finalSolution,finalF,ipos->mf_dispnp);
    loop_interface(structfield,dd,finalSolution);

    // Wenn der letzte Lï¿½eraufruf nicht mit der gefundenen Lï¿½ung
    // erfolgte, mssen die Felder erst einmal neu gelï¿½t werden.

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
  DSTraceHelper("dyn_fsi_nox");
  // Create a communicator for Epetra objects
#if 0
#ifdef PARALLEL
  Epetra_MpiComm Comm(mpicomm);
#else
  Epetra_SerialComm Comm;
#endif
#endif

  // right now the interface is redundant anyway
  Epetra_SerialComm Comm;

  // Create the interface between the test problem and the nonlinear solver
  Teuchos::RefCountPtr<FSI_InterfaceProblem> interface = Teuchos::rcp(new FSI_InterfaceProblem(Comm,struct_work,fluid_work,ale_work));
  interface->timeloop(interface);
}

#endif
#endif
