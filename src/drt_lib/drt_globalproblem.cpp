/*----------------------------------------------------------------------*/
/*!
\file drt_globalproblem.cpp

\brief global list of problems

<pre>
Maintainer: Ulrich Kuettler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15238
</pre>
*/
/*----------------------------------------------------------------------*/

#ifdef CCADISCRET

#include <Teuchos_StandardParameterEntryValidators.hpp>
#include <Teuchos_ParameterListExceptions.hpp>

#include <Epetra_Time.h>

#include "drt_conditiondefinition.H"
#include "drt_function.H"
#include "drt_globalproblem.H"
#include "drt_inputreader.H"
#include "drt_timecurve.H"
#include "drt_utils.H"
#include "drt_validconditions.H"
#include "drt_validparameters.H"

#ifdef PARALLEL
#include <Epetra_MpiComm.h>
#endif

#include <Epetra_SerialComm.h>



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

/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | structure of flags to control output                                 |
 | defined in out_global.c                                              |
 *----------------------------------------------------------------------*/
extern struct _IO_FLAGS     ioflags;

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
 | global variable *solv, vector of lenght numfld of structures SOLVAR  |
 | defined in solver_control.c                                          |
 |                                                                      |
 |                                                       m.gee 11/00    |
 *----------------------------------------------------------------------*/
extern struct _SOLVAR  *solv;

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
extern struct _MATERIAL *mat;


/*----------------------------------------------------------------------*/
// Lena said: do it the easy way.
/*----------------------------------------------------------------------*/
extern "C"
void drt_problem_done()
{
  DRT::Problem::Done();
}


/*----------------------------------------------------------------------*/
// the instances
/*----------------------------------------------------------------------*/
vector<RefCountPtr<DRT::Problem> > DRT::Problem::instances_;


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
RefCountPtr<DRT::Problem> DRT::Problem::Instance(int num)
{
  if (num > static_cast<int>(instances_.size())-1)
  {
    instances_.resize(num+1);
    instances_[num] = rcp(new Problem());
  }
  return instances_[num];
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::Problem::Done()
{
  // This is called at the very end of a baci run.
  //
  // It removes all global problem objects. Therefore all
  // discretizations as well and everything inside those.
  //
  // There is a whole lot going on here...
  instances_.clear();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::Problem::ReadParameter(DRT::INPUT::DatFileReader& reader)
{
  RCP<ParameterList> list = rcp(new ParameterList("DAT FILE"));

  reader.ReadGidSection("--DISCRETISATION", *list);
  reader.ReadGidSection("--PROBLEM SIZE", *list);
  reader.ReadGidSection("--PROBLEM TYP", *list);
  reader.ReadGidSection("--IO", *list);
  reader.ReadGidSection("--DESIGN DESCRIPTION", *list);
  //reader.ReadGidSection("--STATIC", *list);
  //reader.ReadGidSection("--EIGENVALUE ANALYSIS", *list);
  reader.ReadGidSection("--STRUCTURAL DYNAMIC", *list);
  reader.ReadGidSection("--FLUID DYNAMIC", *list);
  reader.ReadGidSection("--FLUID DYNAMIC/STABILIZATION", *list);
  reader.ReadGidSection("--FLUID DYNAMIC/TURBULENCE MODEL", *list);
  reader.ReadGidSection("--ALE DYNAMIC", *list);
  reader.ReadGidSection("--FSI DYNAMIC", *list);

  reader.ReadGidSection("--FLUID SOLVER", *list);
  reader.ReadGidSection("--STRUCT SOLVER", *list);
  reader.ReadGidSection("--ALE SOLVER", *list);
  reader.ReadGidSection("--THERMAL SOLVER", *list);

  // a special section for condition names that contains a list of key-integer
  // pairs but is not validated since the keys are arbitrary.
  reader.ReadGidSection("--CONDITION NAMES", *list);
  list->sublist("CONDITION NAMES").disableRecursiveValidation();

  setParameterList(list);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::Problem::InputControl()
{
  input_ReadGlobalParameterList();

  // Play it save and fill the old C structures here.
  // We have to get rid of them eventually.

  const Teuchos::ParameterList& size = ProblemSizeParams();

  genprob.nele  = size.get<int>("ELEMENTS");
  genprob.nnode = size.get<int>("NODES");
  genprob.ndim  = size.get<int>("DIM");
  genprob.nmat  = size.get<int>("MATERIALS");
  genprob.numdf = size.get<int>("NUMDF");

  if (genprob.nmat<=0)
    dserror("No material defined!");

  const Teuchos::ParameterList& type = ProblemTypeParams();

  genprob.probtyp        = Teuchos::getIntegralValue<PROBLEM_TYP>(type,"PROBLEMTYP");
  genprob.timetyp        = Teuchos::getIntegralValue<TIME_TYP>(type,"TIMETYP");
  genprob.restart        = type.get<int>("RESTART");
  genprob.numfld         = type.get<int>("NUMFIELD");
  genprob.multisc_struct = type.get<int>("MULTISC_STRUCT");

  // set field numbers depending on problem type and numfld
  switch (genprob.probtyp)
  {
  case prb_fsi:
  {
    genprob.numsf=0;
    genprob.numff=1;
    genprob.numaf=2;
    break;
  }
  case prb_fsi_xfem:
  {
    genprob.numsf=0;
    genprob.numff=1;
    if (genprob.numfld==3)
    {
      dserror("ale support tuned of for fsi xfem problems. choose numfld = 2!");
      genprob.numaf=2;
    }
    break;
  }
  case prb_fluid:
  {
    genprob.numff=0;
    if (genprob.numfld==2)
      genprob.numaf=1;
    break;
  }
  case prb_fluid_xfem:
  {
    genprob.numsf=0;
    genprob.numff=1;
    if (genprob.numfld==3)
      genprob.numaf=2;
    break;
  }
  case prb_condif:
    genprob.numff=0;
    break;
  case prb_ale:
    genprob.numaf=0;
    break;
  case prb_structure:
    genprob.numsf=0;
    break;
  case prb_struct_multi:
  {
    genprob.numsf=0;
    break;
  }
  case prb_tsi:
  {
    genprob.numsf = 0;  /* structural field index */
    genprob.numtf = 1;  /* thermal field index */
    break;
  }
  default:
    dserror("problem type %d unknown", genprob.probtyp);
  }

  // input / output file choices

  // anachronisms all around...
  memset(&ioflags, 0, sizeof(ioflags));

  const Teuchos::ParameterList& io = IOParams();

  ioflags.output_out = Teuchos::getIntegralValue<int>(io,"OUTPUT_OUT");
  ioflags.output_gid = Teuchos::getIntegralValue<int>(io,"OUTPUT_GID");
  ioflags.output_bin = Teuchos::getIntegralValue<int>(io,"OUTPUT_BIN");
  ioflags.struct_disp = Teuchos::getIntegralValue<int>(io,"STRUCT_DISP");
  ioflags.struct_stress = Teuchos::getIntegralValue<int>(io,"STRUCT_STRESS");
  ioflags.struct_stress_smo = Teuchos::getIntegralValue<int>(io,"STRUCT_STRESS_SMO");
  ioflags.struct_sm_disp = Teuchos::getIntegralValue<int>(io,"STRUCT_SM_DISP");
  ioflags.struct_sm_stress = Teuchos::getIntegralValue<int>(io,"STRUCT_SM_STRESS");
  ioflags.fluid_sol = Teuchos::getIntegralValue<int>(io,"FLUID_SOL");
  ioflags.fluid_stress = Teuchos::getIntegralValue<int>(io,"FLUID_STRESS");
  ioflags.fluid_vis = Teuchos::getIntegralValue<int>(io,"FLUID_VIS");
  ioflags.ale_disp = Teuchos::getIntegralValue<int>(io,"ALE_DISP");
  ioflags.therm_temper = Teuchos::getIntegralValue<int>(io,"THERM_TEMPERATURE");
  ioflags.therm_heatflux = Teuchos::getIntegralValue<int>(io,"THERM_HEATFLUX");

  ioflags.steps_per_file = io.get<int>("FILESTEPS");

  // input of general solver data

  switch (genprob.probtyp)
  {
  /* for FSI */
  case prb_fsi:
  {
    if (genprob.numfld!=3) dserror("numfld != 3 for FSI");

    solver_.resize(genprob.numfld);
    solv = &solver_[0];

    solv[genprob.numsf].fieldtyp = structure;
    InputSolverControl("STRUCT SOLVER",&(solv[genprob.numsf]));

    solv[genprob.numff].fieldtyp = fluid;
    InputSolverControl("FLUID SOLVER",&(solv[genprob.numff]));

    solv[genprob.numaf].fieldtyp = ale;
    InputSolverControl("ALE SOLVER",&(solv[genprob.numaf]));
    break;
  }
  case prb_fsi_xfem:
  {
    if (genprob.numfld!=2) dserror("numfld != 2 for FSI XFEM problem");

    solver_.resize(genprob.numfld);
    solv = &solver_[0];

    solv[genprob.numsf].fieldtyp = structure;
    InputSolverControl("STRUCT SOLVER",&(solv[genprob.numsf]));

    solv[genprob.numff].fieldtyp = fluid;
    InputSolverControl("FLUID SOLVER",&(solv[genprob.numff]));
    break;
  }
  /* for structure */
  case prb_structure:
  {
    if (genprob.numfld!=1) dserror("numfld != 1 for Structural Problem");

    solver_.resize(genprob.numfld);
    solv = &solver_[0];

    solv[genprob.numsf].fieldtyp = structure;
    InputSolverControl("STRUCT SOLVER",&(solv[genprob.numsf]));
    break;
  }
  /* for fluid */
  case prb_fluid:
  {
    solver_.resize(genprob.numfld);
    solv = &solver_[0];

    solv[genprob.numff].fieldtyp = fluid;
    InputSolverControl("FLUID SOLVER",&(solv[genprob.numff]));

    if (genprob.numfld==2)
    {
      solv[genprob.numaf].fieldtyp = ale;
      InputSolverControl("ALE SOLVER",&(solv[genprob.numaf]));
    }
    break;
  }
  case prb_fluid_xfem:
  {
    dsassert(genprob.numfld == 2, "numfld != 2 for Fluid XFEM problem");

    solver_.resize(genprob.numfld);
    solv = &solver_[0];

    solv[genprob.numsf].fieldtyp = structure;
    InputSolverControl("STRUCT SOLVER",&(solv[genprob.numsf]));

    solv[genprob.numff].fieldtyp = fluid;
    InputSolverControl("FLUID SOLVER",&(solv[genprob.numff]));
    break;
  }
  /* for convection-diffusion */
  case prb_condif:
  {
    solver_.resize(genprob.numfld);
    solv = &solver_[0];

    solv[genprob.numff].fieldtyp = fluid;
    InputSolverControl("FLUID SOLVER",&(solv[genprob.numff]));
    break;
  }
  /* for (plain) ALE */
  case prb_ale:
  {
    if (genprob.numfld!=1) dserror("numfld != 1 for Ale Problem");

    solver_.resize(genprob.numfld);
    solv = &solver_[0];

    solv[genprob.numaf].fieldtyp = ale;
    InputSolverControl("ALE SOLVER",&(solv[genprob.numaf]));
    break;
  }
#ifdef D_TSI
  /* for TSI */
  case prb_tsi:
  {
    if (genprob.numfld != 2)
    {
      dserror("numfld != 2 for TSI: impossible");
    }

    /* allocate memory for solver entries */
    solver_.resize(genprob.numfld);
    solv = &solver_[0];

    /* set solver details of structure solver */
    solv[genprob.numsf].fieldtyp = structure;
    InputSolverControl("STRUCT SOLVER",&(solv[genprob.numsf]));

    /* set solver details of thermal solver */
    solv[genprob.numtf].fieldtyp = thermal;
    InputSolverControl("THERMAL SOLVER",&(solv[genprob.numtf]));
    break;
  }
#endif

  /* for multi-scale (structure) */
  case prb_struct_multi:
  {
    if (genprob.numfld!=1) dserror("numfld != 1 for Structural Multi-Scale Problem");

    solver_.resize(genprob.numfld);
    solv = &solver_[0];

    solv[genprob.numsf].fieldtyp = structure;
    InputSolverControl("STRUCT SOLVER",&(solv[genprob.numsf]));
    break;
  }

  default:
    dserror("problem type %d unknown", genprob.probtyp);
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::Problem::ReadConditions(const DRT::INPUT::DatFileReader& reader)
{
  Epetra_Time time(*reader.Comm());
  if (reader.Comm()->MyPID()==0)
  {
    std::cout << "Read conditions                          in....";
    std::cout.flush();
  }

  /*---------------------------------------------- input of time curves */
  DRT::UTILS::TimeCurveManager::Instance().ReadInput();
  /*---------------------------------------- input of spatial functions */
  DRT::UTILS::FunctionManager::Instance().ReadInput();
  //------------------------------- read number of design objects we have
  // this currently serves to determine how many node sets we might have
  const Teuchos::ParameterList& design = DesignDescriptionParams();
  int ndnode = design.get<int>("NDPOINT");
  int ndline = design.get<int>("NDLINE");
  int ndsurf = design.get<int>("NDSURF");
  int ndvol  = design.get<int>("NDVOL");

  //--------------------------------------------- read generic node sets
  // read design nodes <-> nodes
  vector<vector<int> > dnode_fenode(ndnode);
  reader.ReadDesign("DNODE",dnode_fenode);

  // read design lines <-> nodes
  vector<vector<int> > dline_fenode(ndline);
  reader.ReadDesign("DLINE",dline_fenode);

  // read design surfaces <-> nodes
  vector<vector<int> > dsurf_fenode(ndsurf);
  reader.ReadDesign("DSURF",dsurf_fenode);

  // read design volumes <-> nodes
  vector<vector<int> > dvol_fenode(ndvol);
  reader.ReadDesign("DVOL",dvol_fenode);

  // create list of known conditions
  Teuchos::RCP<std::vector<Teuchos::RCP<DRT::INPUT::ConditionDefinition> > > vc = DRT::INPUT::ValidConditions();
  std::vector<Teuchos::RCP<DRT::INPUT::ConditionDefinition> >& condlist = *vc;

  // test for each condition definition (input file condition section)
  // - read all conditions that match the definition
  // - add the nodal clouds to the conditions
  // - add the conditions to the appropiate discretizations
  //
  // Note that this will reset (un-FillComplete) the discretizations.

  for (unsigned c=0; c<condlist.size(); ++c)
  {
    std::multimap<int,Teuchos::RCP<DRT::Condition> > cond;

    // read conditions from dat file
    condlist[c]->Read(*this,reader,cond);

    // add nodes to conditions
    multimap<int,RefCountPtr<DRT::Condition> >::const_iterator curr;
    for (curr=cond.begin(); curr!=cond.end(); ++curr)
    {
      switch (curr->second->GType())
      {
      case Condition::Point:   curr->second->Add("Node Ids",dnode_fenode[curr->first]); break;
      case Condition::Line:    curr->second->Add("Node Ids",dline_fenode[curr->first]); break;
      case Condition::Surface: curr->second->Add("Node Ids",dsurf_fenode[curr->first]); break;
      case Condition::Volume:  curr->second->Add("Node Ids",dvol_fenode [curr->first]); break;
      default:
        dserror("geometry type unspecified");
      }

      // Iterate through all discretizations and sort the appropiate condition
      // into the correct discretization it applies to
      for (unsigned i=0; i<NumFields(); ++i)
      {
        for (unsigned j=0; j<NumDis(i); ++j)
        {
          Teuchos::RCP<DRT::Discretization> actdis = Dis(i,j);

          const vector<int>* nodes = curr->second->Nodes();
          if (nodes->size()==0)
            dserror("%s condition %d has no nodal cloud",
                    condlist[c]->Description().c_str(),
                    curr->second->Id());

          const int firstnode = (*nodes)[0];
          int foundit = actdis->HaveGlobalNode(firstnode);
          int found=0;
          actdis->Comm().SumAll(&foundit,&found,1);
          if (found)
            actdis->SetCondition(condlist[c]->Name(),curr->second);
        }
      }
    }
  }

  // debug
#if 0
  for (unsigned i=0; i<NumFields(); ++i)
  {
    for (unsigned j=0; j<NumDis(i); ++j)
    {
      for (unsigned c=0; c<condlist.size(); ++c)
      {
        Teuchos::RCP<DRT::Discretization> actdis = Dis(i,j);
        condlist[c]->Print(cout,&*actdis,true);
      }
    }
  }
#endif

  if (reader.Comm()->MyPID()==0)
  {
    std::cout << time.ElapsedTime() << " secs\n";
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::Problem::InputSolverControl(std::string section, _SOLVAR* solv)
{
  const Teuchos::ParameterList& solverparams = getParameterList()->sublist(section);

  solv->matrixtyp = matrix_none;
  solv->solvertyp = Teuchos::getIntegralValue<_SOLVER_TYP>(solverparams,"SOLVER");

  switch (solv->solvertyp)
  {
  case aztec_msr:
  case vm3:
    solv->azvar = new AZVAR;

    solv->azvar->azsolvertyp = Teuchos::getIntegralValue<_AZSOLVERTYP>(solverparams,"AZSOLVE");
    solv->azvar->azprectyp   = Teuchos::getIntegralValue<_AZPRECTYP>(solverparams,"AZPREC");
    solv->azvar->azconv      = Teuchos::getIntegralValue<int>(solverparams,"AZCONV");
    solv->azvar->azscal      = Teuchos::getIntegralValue<int>(solverparams,"AZSCAL");

    solv->azvar->azomega   = solverparams.get<double>("AZOMEGA");
    solv->azvar->azdrop    = solverparams.get<double>("AZDROP" );
    solv->azvar->azfill    = solverparams.get<double>("AZFILL" );
    solv->azvar->aztol     = solverparams.get<double>("AZTOL"  );

    solv->azvar->azoutput  = solverparams.get<int>("AZOUTPUT"  );
    solv->azvar->azreuse   = solverparams.get<int>("AZREUSE"   );
    solv->azvar->azgfill   = solverparams.get<int>("AZGFILL"   );
    solv->azvar->aziter    = solverparams.get<int>("AZITER"    );
    solv->azvar->azsub     = solverparams.get<int>("AZSUB"     );
    solv->azvar->azgraph   = solverparams.get<int>("AZGRAPH"   );
    solv->azvar->azpoly    = solverparams.get<int>("AZPOLY"    );
    solv->azvar->blockdiag = solverparams.get<int>("AZBDIAG"   );
    solv->azvar->azoverlap = solverparams.get<int>("AZOVERLAP" );

    // parameters of ML preconditioner

    solv->azvar->mlprint        = solverparams.get<int>("ML_PRINT");
    solv->azvar->mlcsize        = solverparams.get<int>("ML_MAXCOARSESIZE");
    solv->azvar->mlmaxlevel     = solverparams.get<int>("ML_MAXLEVEL");
    solv->azvar->mlaggsize      = solverparams.get<int>("ML_AGG_SIZE");

    solv->azvar->mldamp_fine    = solverparams.get<double>("ML_DAMPFINE");
    solv->azvar->mldamp_med     = solverparams.get<double>("ML_DAMPMED");
    solv->azvar->mldamp_coarse  = solverparams.get<double>("ML_DAMPCOARSE");
    solv->azvar->mldamp_prolong = solverparams.get<double>("ML_PROLONG_SMO");
    solv->azvar->ml_threshold   = solverparams.get<double>("ML_PROLONG_THRES");

    if (solv->azvar->mlmaxlevel)
    {
      std::string mlsmotimes = solverparams.get<std::string>("ML_SMOTIMES");
      std::istringstream t;
      t.str(mlsmotimes);
      for (int i=0; i<solv->azvar->mlmaxlevel; ++i)
      {
        if (i >= 15)
          dserror("up to 15 levels supported in ML binding");
        t >> solv->azvar->mlsmotimes[i];
      }
    }

    solv->azvar->mlcoarsentype    = Teuchos::getIntegralValue<int>(solverparams,"ML_COARSEN");
    solv->azvar->mlsmotype_fine   = Teuchos::getIntegralValue<int>(solverparams,"ML_SMOOTHERFINE");
    solv->azvar->mlsmotype_med    = Teuchos::getIntegralValue<int>(solverparams,"ML_SMOOTHERMED");
    solv->azvar->mlsmotype_coarse = Teuchos::getIntegralValue<int>(solverparams,"ML_SMOOTHERCOARSE");

    break;
  case parsuperlu:
    solv->psuperluvars = new PSUPERLUVARS;
    break;
  case lapack_sym:
  case lapack_nonsym:
    solv->lapackvars = new LAPACKVARS;
    break;
  default:
    // nothing special to allocate
    break;
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::Problem::ReadFields(DRT::INPUT::DatFileReader& reader)
{
  genprob.create_dis = 0;
  genprob.create_ale = 0;
  genprob.maxnode    = 0;
  genprob.nodeshift  = genprob.nnode;

  // read elements the first time to create graph object
  // row distribution of nodes
  // column distribution of nodes
  // graph of problem
  RefCountPtr<Epetra_Map> rownodes   = null;
  RefCountPtr<Epetra_Map> colnodes   = null;
  RefCountPtr<Epetra_Map> roweles    = null;
  RefCountPtr<Epetra_Map> coleles    = null;
  RefCountPtr<Epetra_CrsGraph> graph = null;

  RefCountPtr<DRT::Discretization> structdis = null;
  RefCountPtr<DRT::Discretization> fluiddis  = null;
  RefCountPtr<DRT::Discretization> aledis    = null;
  RefCountPtr<DRT::Discretization> structdis_macro = null;
  RefCountPtr<DRT::Discretization> structdis_micro = null;


  switch (genprob.probtyp){
  case prb_fsi:
  {
    // allocate and input general old stuff....
    if (genprob.numfld!=3) dserror("numfld != 3 for fsi problem");
    field = (FIELD*)CCACALLOC(genprob.numfld,sizeof(FIELD));
    field[genprob.numsf].fieldtyp = structure;
    field[genprob.numff].fieldtyp = fluid;
    field[genprob.numaf].fieldtyp = ale;

    // obsolete
    field[genprob.numsf].ndis = DiscretisationParams().get<int>("NUMSTRUCDIS");
    field[genprob.numff].ndis = DiscretisationParams().get<int>("NUMFLUIDDIS");
    field[genprob.numaf].ndis = DiscretisationParams().get<int>("NUMALEDIS");

    structdis = rcp(new DRT::Discretization("Structure",reader.Comm()));
    fluiddis = rcp(new DRT::Discretization("Fluid",reader.Comm()));
    aledis = rcp(new DRT::Discretization("Ale",reader.Comm()));

    AddDis(genprob.numsf, structdis);
    AddDis(genprob.numff, fluiddis);
    AddDis(genprob.numaf, aledis);

    DRT::INPUT::NodeReader nodereader(reader, "--NODE COORDS");

    nodereader.AddElementReader(rcp(new DRT::INPUT::ElementReader(structdis, reader, "--STRUCTURE ELEMENTS")));
    nodereader.AddElementReader(rcp(new DRT::INPUT::ElementReader(fluiddis, reader, "--FLUID ELEMENTS")));
    nodereader.AddElementReader(rcp(new DRT::INPUT::ElementReader(aledis, reader, "--ALE ELEMENTS")));

    nodereader.Read();
    break;
  }
  case prb_fsi_xfem:
  {
    // allocate and input general old stuff....
    dsassert(genprob.numfld==2, "numfld != 2 for fluid problem with XFEM interfaces");
    field = (FIELD*)CCACALLOC(genprob.numfld,sizeof(FIELD));
    field[genprob.numsf].fieldtyp = structure;
    field[genprob.numff].fieldtyp = fluid;

    // obsolete
    field[genprob.numsf].ndis = DiscretisationParams().get<int>("NUMSTRUCDIS");
    field[genprob.numff].ndis = DiscretisationParams().get<int>("NUMFLUIDDIS");

    structdis = rcp(new DRT::Discretization("Structure",reader.Comm()));
    fluiddis = rcp(new DRT::Discretization("Fluid",reader.Comm()));

    AddDis(genprob.numsf, structdis);
    AddDis(genprob.numff, fluiddis);

    DRT::INPUT::NodeReader nodereader(reader, "--NODE COORDS");

    nodereader.AddElementReader(rcp(new DRT::INPUT::ElementReader(structdis, reader, "--STRUCTURE ELEMENTS")));
    nodereader.AddElementReader(rcp(new DRT::INPUT::ElementReader(fluiddis, reader, "--FLUID ELEMENTS")));

    nodereader.Read();
    break;

  }
  case prb_ale:
  {
    // allocate and input general old stuff....
    if (genprob.numfld!=1) dserror("numfld != 1 for ale problem");
    field = (FIELD*)CCACALLOC(genprob.numfld,sizeof(FIELD));
    field[genprob.numaf].fieldtyp = ale;

    // obsolete
    field[genprob.numaf].ndis = DiscretisationParams().get<int>("NUMALEDIS");

    aledis = rcp(new DRT::Discretization("Ale",reader.Comm()));
    AddDis(genprob.numaf, aledis);

    DRT::INPUT::NodeReader nodereader(reader, "--NODE COORDS");
    nodereader.AddElementReader(rcp(new DRT::INPUT::ElementReader(aledis, reader, "--ALE ELEMENTS")));
    nodereader.Read();
    break;
  }
  case prb_fluid: case prb_condif:
  {
    // allocate and input general old stuff....
    if (genprob.numfld!=1) dserror("numfld != 1 for fluid problem");
    field = (FIELD*)CCACALLOC(genprob.numfld,sizeof(FIELD));
    field[genprob.numff].fieldtyp = fluid;

    // obsolete
    field[genprob.numff].ndis = DiscretisationParams().get<int>("NUMFLUIDDIS");

    fluiddis = rcp(new DRT::Discretization("Fluid",reader.Comm()));
    AddDis(genprob.numff, fluiddis);

    DRT::INPUT::NodeReader nodereader(reader, "--NODE COORDS");
    nodereader.AddElementReader(rcp(new DRT::INPUT::ElementReader(fluiddis, reader, "--FLUID ELEMENTS")));
    nodereader.Read();
    break;
  }
  case prb_fluid_xfem:
  {
    // allocate and input general old stuff....
    dsassert(genprob.numfld==2, "numfld != 2 for fluid problem with XFEM interfaces");
    field = (FIELD*)CCACALLOC(genprob.numfld,sizeof(FIELD));
    field[genprob.numsf].fieldtyp = structure;
    field[genprob.numff].fieldtyp = fluid;

    // obsolete
    field[genprob.numsf].ndis = DiscretisationParams().get<int>("NUMSTRUCDIS");
    field[genprob.numff].ndis = DiscretisationParams().get<int>("NUMFLUIDDIS");

    structdis = rcp(new DRT::Discretization("Structure",reader.Comm()));
    fluiddis = rcp(new DRT::Discretization("Fluid",reader.Comm()));

    AddDis(genprob.numsf, structdis);
    AddDis(genprob.numff, fluiddis);

    DRT::INPUT::NodeReader nodereader(reader, "--NODE COORDS");

    nodereader.AddElementReader(rcp(new DRT::INPUT::ElementReader(structdis, reader, "--STRUCTURE ELEMENTS")));
    nodereader.AddElementReader(rcp(new DRT::INPUT::ElementReader(fluiddis, reader, "--FLUID ELEMENTS")));

    nodereader.Read();
    break;
  }
  case prb_fluid_pm:
  {
    dserror("prb_fluid_pm not yet impl.");
    break;
  }
  case prb_tsi:
  {
    dserror("prb_tsi not yet impl.");
    break;
  }
  case prb_structure:
  {
    // allocate and input general old stuff....
    if (genprob.numfld!=1) dserror("numfld != 1 for structural problem");
    field = (FIELD*)CCACALLOC(genprob.numfld,sizeof(FIELD));
    field[genprob.numsf].fieldtyp = structure;

    // obsolete
    field[genprob.numsf].ndis = DiscretisationParams().get<int>("NUMSTRUCDIS");

    structdis = rcp(new DRT::Discretization("Structure",reader.Comm()));
    AddDis(genprob.numsf, structdis);

    DRT::INPUT::NodeReader nodereader(reader, "--NODE COORDS");
    nodereader.AddElementReader(rcp(new DRT::INPUT::ElementReader(structdis, reader, "--STRUCTURE ELEMENTS")));
    nodereader.Read();
    break;
  } // end of else if (genprob.probtyp==prb_structure)
  case prb_struct_multi:
  {
    // allocate and input general old stuff....

    if (genprob.numfld!=1) dserror("numfld != 1 for structural multi-scale problem");
    field = (FIELD*)CCACALLOC(genprob.numfld,sizeof(FIELD));
    field[genprob.numsf].fieldtyp = structure;

    // obsolete
    field[genprob.numsf].ndis = DiscretisationParams().get<int>("NUMSTRUCDIS");

    // read macroscale fields from main inputfile

    structdis_macro = rcp(new DRT::Discretization("Macro Structure",reader.Comm()));
    AddDis(genprob.numsf, structdis_macro);

    DRT::INPUT::NodeReader nodereader(reader, "--NODE COORDS");
    nodereader.AddElementReader(rcp(new DRT::INPUT::ElementReader(structdis_macro, reader, "--STRUCTURE ELEMENTS")));
    nodereader.Read();

    // read microscale fields from second inputfile

    RefCountPtr<DRT::Problem> micro_problem = DRT::Problem::Instance(1);
    RefCountPtr<Epetra_SerialComm> serialcomm = rcp(new Epetra_SerialComm());

    string micro_inputfile_name = mat->m.struct_multiscale->micro_inputfile_name;

    if (micro_inputfile_name[0]!='/')
    {
      string filename = reader.MyInputfileName();
      string::size_type pos = filename.rfind('/');
      if (pos!=string::npos)
      {
        string path = filename.substr(0,pos+1);
        micro_inputfile_name.insert(micro_inputfile_name.begin(), path.begin(), path.end());
      }
    }

    if (!structdis_macro->Comm().MyPID())
        cout << "input for microscale is read from        " << micro_inputfile_name << "\n";

    DRT::INPUT::DatFileReader micro_reader(micro_inputfile_name, serialcomm, 1);
    micro_reader.Activate();

    structdis_micro = rcp(new DRT::Discretization("Micro Structure", micro_reader.Comm()));
    micro_problem->AddDis(genprob.numsf, structdis_micro);

    micro_problem->ReadParameter(micro_reader);

    // read materials of microscale
    // CAUTION: materials for microscale can not be read until
    // micro_reader is activated, since else materials will again be
    // read from macroscale inputfile. Besides, materials MUST be read
    // before elements are read since elements establish a connection
    // to the corresponding material! Thus do not change position of
    // function calls!

    micro_problem->ReadMaterial();


    DRT::INPUT::NodeReader micronodereader(micro_reader, "--NODE COORDS");
    micronodereader.AddElementReader(rcp(new DRT::INPUT::ElementReader(structdis_micro, micro_reader, "--STRUCTURE ELEMENTS")));
    micronodereader.Read();

    // read conditions of microscale -> note that no time curves and
    // spatial functions can be read!

    micro_problem->ReadConditions(micro_reader);

    // At this point, everything for the microscale is read,
    // subsequent reading is only for macroscale
    structdis_micro->FillComplete();


    // reactivate reader of macroscale as well as macroscale material

    reader.Activate();
    ActivateMaterial();

    if (0)
    {
    for (int i=0; i<structdis_macro->Comm().NumProc(); ++i)
    {
      cout << "\n";
      if (structdis_macro->Comm().MyPID()==i)
      {
        cout << "Proc " << i << "meine serielle Diskretisierung:\n";
        cout << *structdis_micro;
      }
      structdis_macro->Comm().Barrier();
    }
    cout << "\n";
    }
    //exit(0);
    break;
  } // end of else if (genprob.probtyp==prb_struct_multi)
  default:
    dserror("Type of problem unknown");
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::Problem::setParameterList(Teuchos::RCP< Teuchos::ParameterList > const &paramList)
{
  try {

    // Test parameter list against valid parameters, set default values
    // and set validator objects to extract numerical values for string
    // parameters.
    paramList->validateParametersAndSetDefaults(*this->getValidParameters());
  }
  catch (Teuchos::Exceptions::InvalidParameter& err)
  {
    std::cerr << "\n\n" << err.what();
    dserror("Input parameter validation failed. Fix your input file.");
  }

  // yes, it is my list
  setMyParamList(paramList);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Teuchos::ParameterList> DRT::Problem::getValidParameters() const
{
  // call the external method to get the valid parameters
  // this way the parameter configuration is separate from the source
  return DRT::INPUT::ValidParameters();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
RCP<DRT::Discretization> DRT::Problem::Dis(int fieldnum, int disnum) const
{
  return discretizations_[fieldnum][disnum];
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::Problem::AddDis(int fieldnum, RefCountPtr<Discretization> dis)
{
  if (fieldnum > static_cast<int>(discretizations_.size())-1)
  {
    discretizations_.resize(fieldnum+1);
  }
  discretizations_[fieldnum].push_back(dis);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::Problem::AddMaterial(const _MATERIAL& m)
{
  if (m.m.fluid==NULL)
    dserror("invalid material added");
  material_.push_back(m); // copy!
  ActivateMaterial();     // always reset!
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::Problem::SetDis(int fieldnum, int disnum, RefCountPtr<Discretization> dis)
{
  if (fieldnum > static_cast<int>(discretizations_.size())-1)
  {
    discretizations_.resize(fieldnum+1);
  }
  if (disnum > static_cast<int>(discretizations_[fieldnum].size()-1))
  {
    discretizations_[fieldnum].resize(disnum+1);
  }
  discretizations_[fieldnum][disnum] = dis;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::Problem::ActivateMaterial()
{
  mat = &material_[0];
}


#endif
