
#ifdef TRILINOS_PACKAGE
#ifdef D_FSI

#include "fsi_dyn_nox.H"
#include "fsi_aitken_nox.H"
#include "../discret/dstrc.H"

#include <fstream>
#include <iostream>
#include <string>
#include <sstream>

#include <NOX_Solver_Manager.H>
#include <NOX_Epetra_Vector.H>

#include <Teuchos_XMLParameterListHelpers.hpp>

/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | structure of flags to control output                                 |
 | defined in out_global.c                                              |
 *----------------------------------------------------------------------*/
extern struct _IO_FLAGS     ioflags;

/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | struct _GENPROB       genprob; defined in global_control.c           |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;

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


void debug_out_data(FIELD *actfield, CHAR* n, NODE_ARRAY array, INT pos)
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

#if 0
  /* dann noch schnell die python Variante */
  sprintf(name, "plot/%s%d.py", n, call);
  f = fopen(name, "w");
  fprintf(f, "#from Matrix import Matrix\n");
  fprintf(f, "array = [\n");
  for (i=0; i<actfield->dis[0].numnp; ++i)
  {
    INT k;
    NODE* actnode = &(actfield->dis[0].node[i]);

    fprintf(f, "[");

    for (k=0; k<actnode->numdf; ++k)
    {
      switch (array)
      {
      case node_array_sol:
        fprintf(f, "%20.20e, ",
                actnode->sol.a.da[pos][k]);
        break;
      case node_array_sol_increment:
        fprintf(f, "%20.20e, ",
                actnode->sol_increment.a.da[pos][k]);
        break;
      case node_array_sol_mf:
        fprintf(f, "%20.20e, ",
                actnode->sol_mf.a.da[pos][k]);
        break;
      default:
        dserror("node array %d unsupported", array);
      }
    }

    fprintf(f, "],\n");
  }
  fprintf(f, "]\n\n");
  fclose(f);
#endif
  call++;
}


/*----------------------------------------------------------------------*/
//! \brief loop the interface and call the operator for all dofs of
//         this processor.
/*----------------------------------------------------------------------*/
template <class Op>
void loop_interface(FIELD* structfield, Op& op, const Epetra_Vector& x)
{
  DSTraceHelper("loop_interface");

  int disnum = 0;
  int numnp_total = structfield->dis[disnum].numnp;
  int numaf       = genprob.numaf;

  // We don't have a direct mapping, so loop all structural nodes and
  // select the dofs that belong to the local part of the interface.

  int count = 0;
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
      if (x.Map().MyGID(count))
      {
        op(actsnode,j,x.Map().LID(count));
      }
      count += 1;
    }
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
class DistributeDisplacements
{
  const Epetra_Vector& x_;
  int pos_;

public:
  DistributeDisplacements(const Epetra_Vector& x,int pos)
    : x_(x), pos_(pos) {}

  void operator()(NODE* actsnode, int dof, int LID)
    { actsnode->sol_mf.a.da[pos_][dof] = x_[LID]; }
};


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
class DistributeSolution
{
  const Epetra_Vector& x_;
  const Epetra_Vector& F_;
  int pos_;

public:
  DistributeSolution(const Epetra_Vector& x,const Epetra_Vector& F,int pos)
    : x_(x), F_(F), pos_(pos) {}

  void operator()(NODE* actsnode, int dof, int LID)
    { actsnode->sol_mf.a.da[pos_][dof] = F_[LID]+x_[LID]; }
};


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
class GatherResiduum
{
  Epetra_Vector& F_;
  int dispnp_;
  int dispn_;

public:
  GatherResiduum(Epetra_Vector& F,int dispnp,int dispn)
    : F_(F), dispnp_(dispnp), dispn_(dispn) {}

  void operator()(NODE* actsnode, int dof, int LID)
    {
      F_[LID] = actsnode->sol_mf.a.da[dispnp_][dof] -
                actsnode->sol_mf.a.da[dispn_ ][dof];
    }
};


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
class GatherDisplacements
{
  Epetra_Vector& soln_;
  int pos_;

public:
  GatherDisplacements(Epetra_Vector& soln,int pos)
    : soln_(soln), pos_(pos) {}

  void operator()(NODE* actsnode, int dof, int LID)
    { soln_[LID] = actsnode->sol_mf.a.da[pos_][dof]; }
};


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
class OutputDisplacements
{
  std::ostream& out_;
  int dispnp_;
  int dispn_;

public:
  OutputDisplacements(std::ostream& out,int dispnp,int dispn)
    : out_(out), dispnp_(dispnp), dispn_(dispn) {}

  // a hack
  void operator()(NODE* actsnode, int dof, int LID)
    { if (dof==0)
        out_ << LID << " "
             << actsnode->sol_mf.a.da[dispnp_][0] << " "
             << actsnode->sol_mf.a.da[dispnp_][1] << " "
             << actsnode->sol_mf.a.da[dispn_][0] << " "
             << actsnode->sol_mf.a.da[dispn_][1] << " "
             << actsnode->sol_mf.a.da[dispnp_][0]-actsnode->sol_mf.a.da[dispn_][0] << " "
             << actsnode->sol_mf.a.da[dispnp_][1]-actsnode->sol_mf.a.da[dispn_][1] << "\n"; }
};





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
  DSTraceHelper("FSI_InterfaceProblem::FSI_InterfaceProblem");
  fsidyn_ = alldyn[3].fsidyn;

  StandardMap_ = Teuchos::rcp(new Epetra_Map(fsidyn_->numsid, 0, Comm));
  soln_ = Teuchos::rcp(new Epetra_Vector(*StandardMap_));
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool FSI_InterfaceProblem::computeF(const Epetra_Vector& x,
                                    Epetra_Vector& F,
                                    const FillType fillFlag)
{
  DSTraceHelper("FSI_InterfaceProblem::computeF");
  FIELD* structfield = struct_work_->struct_field;
  FIELD* fluidfield  = fluid_work_ ->fluid_field;
  FIELD* alefield    = ale_work_   ->ale_field;

  int disnum = 0;

  // set new interface displacement x
  ARRAY_POSITION *ipos = &(structfield->dis[disnum].ipos);

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

  // subdivision is not supported.
  INT a_disnum_calc = 0;
  INT a_disnum_io   = 0;
  INT f_disnum_calc = 0;
  INT f_disnum_io   = 0;
  INT s_disnum_calc = 0;
  INT s_disnum_io   = 0;

  // no consecutive iterations here. We have to care about the fields
  // backup from the outside.
  INT itnum = 1;

  /*------------------------------- CMD -------------------------------*/
  perf_begin(44);
  fsi_ale_calc(ale_work_,alefield,a_disnum_calc,a_disnum_io,structfield,s_disnum_calc);
  perf_begin(44);

  //debug_out_data(alefield, "ale_disp", node_array_sol_mf, alefield->dis[0].ipos.mf_dispnp);

  /*------------------------------- CFD -------------------------------*/
  perf_begin(42);
  fsi_fluid_calc(fluid_work_,fluidfield,f_disnum_calc,f_disnum_io,alefield,a_disnum_calc);
  perf_end(42);

  //debug_out_data(fluidfield, "fluid_vel", node_array_sol_increment, fluidfield->dis[0].ipos.velnp);

  /*------------------------------- CSD -------------------------------*/
  perf_begin(43);
  fsi_struct_calc(struct_work_,structfield,s_disnum_calc,s_disnum_io,itnum,fluidfield,f_disnum_calc);
  perf_end(43);

  //debug_out_data(structfield, "struct_disp", node_array_sol, 0);

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

    // backup current structural values
    solserv_sol_copy(structfield,s_disnum_calc,node_array_sol,node_array_sol,0,9);
    solserv_sol_copy(structfield,s_disnum_calc,node_array_sol,node_array_sol,1,10);

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
    Teuchos::RefCountPtr<NOX::Epetra::MatrixFree> MF = Teuchos::rcp(new NOX::Epetra::MatrixFree(printParams, interface, noxSoln));
    // 3. Finite Difference (Epetra_RowMatrix)
    //Teuchos::RefCountPtr<NOX::Epetra::FiniteDifference> FD = Teuchos::rcp(new NOX::Epetra::FiniteDifference(printParams, interface, noxSoln));

    // Create the linear system
    Teuchos::RefCountPtr<NOX::Epetra::Interface::Required> iReq = interface;
    Teuchos::RefCountPtr<NOX::Epetra::Interface::Jacobian> iJac = MF;
    //Teuchos::RefCountPtr<NOX::Epetra::Interface::Preconditioner> iPrec = interface;
    Teuchos::RefCountPtr<NOX::Epetra::LinearSystemAztecOO> linSys = Teuchos::rcp(new NOX::Epetra::LinearSystemAztecOO(printParams, lsParams,
                                                                                                                      iReq, iJac, MF,
                                                                                                                      noxSoln));

    // Create the Group
    Teuchos::RefCountPtr<NOX::Epetra::Group> grp = Teuchos::rcp(new NOX::Epetra::Group(printParams, iReq, noxSoln,
                                                                                       linSys));

    // Create the convergence tests
    Teuchos::RefCountPtr<NOX::StatusTest::NormF> absresid    = Teuchos::rcp(new NOX::StatusTest::NormF(1.0e-8));
    //Teuchos::RefCountPtr<NOX::StatusTest::NormF> relresid    = Teuchos::rcp(new NOX::StatusTest::NormF(*grp.get(), 1.0e-2));
    Teuchos::RefCountPtr<NOX::StatusTest::NormUpdate> update = Teuchos::rcp(new NOX::StatusTest::NormUpdate(1.0e-5));
    Teuchos::RefCountPtr<NOX::StatusTest::NormWRMS> wrms     = Teuchos::rcp(new NOX::StatusTest::NormWRMS(1.0e-2, 1.0e-8));
    Teuchos::RefCountPtr<NOX::StatusTest::Combo> converged   = Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::AND));
    converged->addStatusTest(absresid);
    //converged->addStatusTest(relresid);
    //converged->addStatusTest(wrms);
    converged->addStatusTest(update);
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

    // Wenn der letzte Löseraufruf nicht mit der gefundenen Lösung
    // erfolgte, müssen die Felder erst einmal neu gelöst werden.

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
#ifdef PARALLEL
  Epetra_MpiComm Comm(mpicomm);
#else
  Epetra_SerialComm Comm;
#endif

  // Create the interface between the test problem and the nonlinear solver
  Teuchos::RefCountPtr<FSI_InterfaceProblem> interface = Teuchos::rcp(new FSI_InterfaceProblem(Comm,struct_work,fluid_work,ale_work));
  interface->timeloop(interface);
}

#endif
#endif
