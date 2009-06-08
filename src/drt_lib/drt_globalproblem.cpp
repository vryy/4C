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
#include "drt_materialdefinition.H"
#include "drt_function.H"
#include "drt_globalproblem.H"
#include "drt_inputreader.H"
#include "drt_timecurve.H"
#include "drt_utils.H"
#include "../drt_mat/material.H"
#include "../drt_mat/matpar_bundle.H"
#include "../drt_inpar/drt_validconditions.H"
#include "../drt_inpar/drt_validparameters.H"
#include "../drt_inpar/drt_validmaterials.H"
#include "../drt_mat/micromaterial.H"

#include "../drt_io/io_control.H"

#ifdef PARALLEL
#include <Epetra_MpiComm.h>
#endif

#include <Epetra_SerialComm.h>



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
\brief ranks and communicators

<pre>                                                         m.gee 8/00
This structure struct _PAR par; is defined in main_ccarat.c
and the type is in partition.h
</pre>

*----------------------------------------------------------------------*/
extern struct _PAR   par;




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
std::string DRT::Problem::ProblemType() const
{
  static const char* problemnames[] = PROBLEMNAMES;
  return problemnames[genprob.probtyp];
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::string DRT::Problem::SpatialApproximation() const
{
  // decide which kind of spatial representation is required
  const Teuchos::ParameterList& ptype = ProblemTypeParams();

  std::string basis_fct_type = ptype.get<std::string>("SHAPEFCT");

  return basis_fct_type;
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
  reader.ReadGidSection("--STRUCTURAL DYNAMIC/TIMEADAPTIVITY", *list);
  reader.ReadGidSection("--STRUCTURAL DYNAMIC/GENALPHA", *list);
  reader.ReadGidSection("--STRUCTURAL DYNAMIC/ONESTEPTHETA", *list);
  reader.ReadGidSection("--STRUCTURAL DYNAMIC/GEMM", *list);
  reader.ReadGidSection("--INVERSE ANALYSIS", *list);
  reader.ReadGidSection("--STRUCTURAL CONTACT", *list);
  reader.ReadGidSection("--BROWNIAN MOTION", *list);
  reader.ReadGidSection("--STATISTICAL MECHANICS", *list);
  reader.ReadGidSection("--FLUID DYNAMIC", *list);
  reader.ReadGidSection("--FLUID DYNAMIC/STABILIZATION", *list);
  reader.ReadGidSection("--FLUID DYNAMIC/TURBULENCE MODEL", *list);
  reader.ReadGidSection("--COMBUSTION CONTROL", *list);
  reader.ReadGidSection("--COMBUSTION CONTROL/COMBUSTION FLUID", *list);
  reader.ReadGidSection("--COMBUSTION CONTROL/COMBUSTION GFUNCTION", *list);
  reader.ReadGidSection("--SCALAR TRANSPORT DYNAMIC", *list);
  reader.ReadGidSection("--SCALAR TRANSPORT DYNAMIC/NONLINEAR", *list);
  reader.ReadGidSection("--SCALAR TRANSPORT DYNAMIC/STABILIZATION", *list);
  reader.ReadGidSection("--ALE DYNAMIC", *list);
  reader.ReadGidSection("--FSI DYNAMIC", *list);
  reader.ReadGidSection("--ARTERIAL DYNAMIC", *list);
  reader.ReadGidSection("--XFEM GENERAL", *list);
  reader.ReadGidSection("--LOMA CONTROL", *list);
  reader.ReadGidSection("--ELCH CONTROL", *list);

  reader.ReadSection("--STRUCT NOX", *list);
  reader.ReadSection("--STRUCT NOX/Direction", *list);
  reader.ReadSection("--STRUCT NOX/Direction/Newton", *list);
  reader.ReadSection("--STRUCT NOX/Direction/Steepest Descent", *list);
  reader.ReadSection("--STRUCT NOX/Line Search", *list);
  reader.ReadSection("--STRUCT NOX/Line Search/Full Step", *list);
  reader.ReadSection("--STRUCT NOX/Line Search/Backtrack", *list);
  reader.ReadSection("--STRUCT NOX/Line Search/Polynomial", *list);
  reader.ReadSection("--STRUCT NOX/Line Search/More'-Thuente", *list);
  reader.ReadSection("--STRUCT NOX/Trust Region", *list);
  reader.ReadSection("--STRUCT NOX/Printing", *list);

  reader.ReadGidSection("--FLUID SOLVER", *list);
  reader.ReadGidSection("--FLUID PRESSURE SOLVER", *list);
  reader.ReadGidSection("--STRUCT SOLVER", *list);
  reader.ReadGidSection("--ALE SOLVER", *list);
  reader.ReadGidSection("--THERMAL SOLVER", *list);
  reader.ReadGidSection("--SCALAR TRANSPORT SOLVER", *list);
  reader.ReadGidSection("--SCALAR TRANSPORT ELECTRIC POTENTIAL SOLVER", *list);
  reader.ReadGidSection("--ARTERY NETWORK SOLVER", *list);


  // a special section for condition names that contains a list of key-integer
  // pairs but is not validated since the keys are arbitrary.
  reader.ReadGidSection("--CONDITION NAMES", *list);

  setParameterList(list);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::Problem::InputControl()
{
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
  case prb_fluid_xfem:
  {
    genprob.numsf=0;
    genprob.numff=1;
    genprob.numaf=2;
    break;
  }
  case prb_fluid:
  {
    genprob.numff=0;
    if (genprob.numfld==2)
      genprob.numaf=1;
    break;
  }
  case prb_fluid_ale:
  {
    genprob.numff=0;
    genprob.numaf=1;
    break;
  }
  case prb_freesurf:
  {
    genprob.numff=0;
    genprob.numaf=1;
    break;
  }
  case prb_scatra:
    genprob.numff = 0;  /* fluid field index */
    genprob.numscatra = 1; /* scalar transport field index */
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
  case prb_loma:
  {
    genprob.numff = 0;  /* fluid field index */
    genprob.numscatra = 1; /* scalar transport field index */
    break;
  }
  case prb_elch:
  {
    genprob.numff = 0;  /* fluid field index */
    genprob.numscatra = 1; /* scalar transport field index */
    genprob.numaf=2; /* ALE field index */
    break;
  }
  case prb_combust:
  {
    genprob.numff = 0;  /* fluid field index */
    genprob.numscatra = 1; /* scalar transport field (=G-function) index */
    break;
  }
  case prb_art_net:
  {
    genprob.numartf = 0;  /* arterial network field index */
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
//   ioflags.struct_stress = Teuchos::getIntegralValue<int>(io,"STRUCT_STRESS");
  ioflags.struct_sm_disp = Teuchos::getIntegralValue<int>(io,"STRUCT_SM_DISP");
  ioflags.fluid_sol = Teuchos::getIntegralValue<int>(io,"FLUID_SOL");
  ioflags.fluid_stress = Teuchos::getIntegralValue<int>(io,"FLUID_STRESS");
  ioflags.fluid_vis = Teuchos::getIntegralValue<int>(io,"FLUID_VIS");
  ioflags.ale_disp = Teuchos::getIntegralValue<int>(io,"ALE_DISP");
  ioflags.therm_temper = Teuchos::getIntegralValue<int>(io,"THERM_TEMPERATURE");
  ioflags.therm_heatflux = Teuchos::getIntegralValue<int>(io,"THERM_HEATFLUX");

  ioflags.steps_per_file = io.get<int>("FILESTEPS");

}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::Problem::ReadMaterials(const DRT::INPUT::DatFileReader& reader)
{
  if (materials_ != Teuchos::null)
    dserror("Trying to create again the material bundle");
  materials_ = Teuchos::rcp(new MAT::PAR::Bundle());

  // create list of known materials
  Teuchos::RCP<std::vector<Teuchos::RCP<DRT::INPUT::MaterialDefinition> > > vm = DRT::INPUT::ValidMaterials();
  std::vector<Teuchos::RCP<DRT::INPUT::MaterialDefinition> >& matlist = *vm;

  // test for each material definition (input file --MATERIALS section)
  // and store in #matmap_
  for (unsigned m=0; m<matlist.size(); ++m)
  {
    // read material from DAT file of type #matlist[m]
    matlist[m]->Read(*this,reader,materials_);
  }

  // check if every material was identified
  const std::string name = "--MATERIALS";
  std::vector<const char*> section = reader.Section(name);
  int nummat = 0;
  if (section.size() > 0)
  {
    for (std::vector<const char*>::iterator i=section.begin();
         i!=section.end();
         ++i)
    {
      Teuchos::RCP<std::stringstream> condline = Teuchos::rcp(new std::stringstream(*i));

      std::string mat;
      std::string number;
      std::string name;
      (*condline) >> mat >> number >> name;
      if ( (not (*condline)) or (mat != "MAT") )
        dserror("invalid material line in '%s'",name.c_str());

      // extract material ID
      int matid = -1;
      {
        char* ptr;
        matid = strtol(number.c_str(),&ptr,10);
        if (ptr == number.c_str())
          dserror("failed to read material object number '%s'",
                  number.c_str());
      }

      // processed?
      if (materials_->Find(matid) == -1)
        dserror("Material 'MAT %d' with name '%s' could not be identified", matid, name.c_str());

      // count number of materials provided in file
      nummat += 1;
    }
  }

  // make fast access parameters
  materials_->MakeParameters();


  // inform user
  //std::cout << "Number of successfully read materials is " << nummat << std::endl;
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
    multimap<int,RCP<DRT::Condition> >::const_iterator curr;
    for (curr=cond.begin(); curr!=cond.end(); ++curr)
    {
      switch (curr->second->GType())
      {
      case Condition::Point:
        if (curr->first < 0 or static_cast<unsigned>(curr->first) >= dnode_fenode.size())
          dserror("DPoint %d not in range [0:%d[\n"
                  "DPoint condition on non existent DPoint?",
                  curr->first,dnode_fenode.size());
        curr->second->Add("Node Ids",dnode_fenode[curr->first]);
        break;
      case Condition::Line:
        if (curr->first < 0 or static_cast<unsigned>(curr->first) >= dline_fenode.size())
          dserror("DLine %d not in range [0:%d[\n"
                  "DLine condition on non existent DLine?",
                  curr->first,dline_fenode.size());
        curr->second->Add("Node Ids",dline_fenode[curr->first]);
        break;
      case Condition::Surface:
        if (curr->first < 0 or static_cast<unsigned>(curr->first) >= dsurf_fenode.size())
          dserror("DSurface %d not in range [0:%d[\n"
                  "DSurface condition on non existent DSurface?",
                  curr->first,dsurf_fenode.size());
        curr->second->Add("Node Ids",dsurf_fenode[curr->first]);
        break;
      case Condition::Volume:
        if (curr->first < 0 or static_cast<unsigned>(curr->first) >= dvol_fenode.size())
          dserror("DVolume %d not in range [0:%d[\n"
                  "DVolume condition on non existent DVolume?",
                  curr->first,dvol_fenode.size());
        curr->second->Add("Node Ids",dvol_fenode [curr->first]);
        break;
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

          int foundit = 0;
          for (unsigned i=0; i<nodes->size(); ++i)
          {
            const int node = (*nodes)[i];
            foundit = actdis->HaveGlobalNode(node);
            if (foundit)
              break;
          }
          int found=0;
          actdis->Comm().SumAll(&foundit,&found,1);
          if (found)
	    // Insert a copy since we might insert the same condition in many discretizations.
            actdis->SetCondition(condlist[c]->Name(),Teuchos::rcp(new Condition(*curr->second)));
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
void DRT::Problem::ReadKnots(const DRT::INPUT::DatFileReader& reader)
{
  // decide which kind of spatial representation is required
  const Teuchos::ParameterList& ptype = ProblemTypeParams();
  // get information on the spatial approximation --- we only read knots
  // in the nurbs case
  std::string distype = ptype.get<std::string>("SHAPEFCT");

  // get problem dimension
  const Teuchos::ParameterList& psize= ProblemSizeParams();

  int dim = psize.get<int>("DIM");

  // Iterate through all discretizations and sort the appropiate condition
  // into the correct discretization it applies to
  
  for (unsigned i=0; i<NumFields(); ++i)
  {
    for (unsigned j=0; j<NumDis(i); ++j)
    {
      Teuchos::RCP<DRT::Discretization> actdis = Dis(i,j);

      if(distype == "Nurbs")
      {
        // cast discretisation to nurbs variant to be able
        // to add the knotvector
        DRT::NURBS::NurbsDiscretization* nurbsdis
          =
          dynamic_cast<DRT::NURBS::NurbsDiscretization*>(&(*actdis));

        // define an empty knot vector object
        Teuchos::RCP<DRT::NURBS::Knotvector> disknots=Teuchos::null;

        // read the knotvector data from the input
        reader.ReadKnots(dim,actdis->Name(),disknots);

        if(disknots==Teuchos::null)
        {
          dserror("Knotvector read failed in Nurbs discretisation\n");
        }
        
        // make sure atdis is fillcompleted, to be able to call 
        // ElementRowMap() on it
        if(!actdis->Filled())
        {
          actdis->FillComplete();
        }

        // the smallest gid in the discretisation determines the access
        // pattern via the element offset
        int smallest_gid_in_dis=actdis->ElementRowMap()->MinAllGID();

        // consistency checks
        disknots->FinishKnots(smallest_gid_in_dis);

        // add knots to discretisation
        nurbsdis->SetKnotVector(disknots);
      }
    } //loop discretisations of field
  } //loop fields

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::Problem::OpenControlFile(const Epetra_Comm& comm, std::string inputfile, std::string prefix)
{
  if (genprob.restart)
    inputcontrol_ = Teuchos::rcp(new IO::InputControl(prefix));

  outputcontrol_ = Teuchos::rcp(new IO::OutputControl(comm,
                                                      ProblemType(),
                                                      SpatialApproximation(),
                                                      inputfile,
                                                      prefix,
                                                      genprob.ndim,
                                                      genprob.restart,
                                                      IOParams().get<int>("FILESTEPS")));

}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::Problem::OpenErrorFile(const Epetra_Comm& comm, std::string prefix)
{
  errorfilecontrol_ = Teuchos::rcp(new IO::ErrorFileControl(comm, prefix));
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
  RCP<Epetra_Map> rownodes   = null;
  RCP<Epetra_Map> colnodes   = null;
  RCP<Epetra_Map> roweles    = null;
  RCP<Epetra_Map> coleles    = null;
  RCP<Epetra_CrsGraph> graph = null;

  RCP<DRT::Discretization> structdis       = null;
  RCP<DRT::Discretization> fluiddis        = null;
  RCP<DRT::Discretization> aledis          = null;
  RCP<DRT::Discretization> boundarydis     = null;
  RCP<DRT::Discretization> scatradis       = null;
  RCP<DRT::Discretization> structdis_macro = null;
  RCP<DRT::Discretization> structdis_micro = null;
  RCP<DRT::Discretization> arterydis       = null; //_1D_ARTERY_

  // decide which kind of spatial representation is required
  const Teuchos::ParameterList& ptype = ProblemTypeParams();

  switch (genprob.probtyp){
  case prb_fsi:
  {
    // allocate and input general old stuff....
    if (genprob.numfld!=3) dserror("numfld != 3 for fsi problem");

    std::string distype = ptype.get<std::string>("SHAPEFCT");

    if(distype == "Nurbs")
    {
      structdis = rcp(new DRT::NURBS::NurbsDiscretization("structure",reader.Comm()));
      fluiddis  = rcp(new DRT::NURBS::NurbsDiscretization("fluid"    ,reader.Comm()));
      aledis    = rcp(new DRT::NURBS::NurbsDiscretization("ale"      ,reader.Comm()));
    }
    else
    {
      structdis = rcp(new DRT::Discretization("structure",reader.Comm()));
      fluiddis  = rcp(new DRT::Discretization("fluid"    ,reader.Comm()));
      aledis    = rcp(new DRT::Discretization("ale"      ,reader.Comm()));
    }


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
  case prb_fluid_xfem:
  {
    structdis = rcp(new DRT::Discretization("structure",reader.Comm()));
    fluiddis = rcp(new DRT::Discretization("fluid",reader.Comm()));

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

    std::string distype = ptype.get<std::string>("SHAPEFCT");

    if(distype == "Nurbs")
    {
      aledis = rcp(new DRT::NURBS::NurbsDiscretization("ale",reader.Comm()));
    }
    else
    {
      aledis = rcp(new DRT::Discretization("ale",reader.Comm()));
    }

    AddDis(genprob.numaf, aledis);

    DRT::INPUT::NodeReader nodereader(reader, "--NODE COORDS");
    nodereader.AddElementReader(rcp(new DRT::INPUT::ElementReader(aledis, reader, "--ALE ELEMENTS")));
    nodereader.Read();
    break;
  }
  case prb_fluid:
  {
    // allocate and input general old stuff....
    if (genprob.numfld!=1) dserror("numfld != 1 for fluid problem");

    std::string distype = ptype.get<std::string>("SHAPEFCT");

    if(distype == "Nurbs")
    {
      fluiddis = rcp(new DRT::NURBS::NurbsDiscretization("fluid",reader.Comm()));
    }
    else
    {
      fluiddis = rcp(new DRT::Discretization("fluid",reader.Comm()));
    }

    AddDis(genprob.numff, fluiddis);

    DRT::INPUT::NodeReader nodereader(reader, "--NODE COORDS");
    nodereader.AddElementReader(rcp(new DRT::INPUT::ElementReader(fluiddis, reader, "--FLUID ELEMENTS")));
    nodereader.Read();
    break;
  }
  case prb_scatra:
  {
    // allocate and input general old stuff....
    if (genprob.numfld!=2) dserror("numfld != 2 for scalar transport problem");

    // create empty discretizations
    fluiddis = rcp(new DRT::Discretization("fluid",reader.Comm()));
    AddDis(genprob.numff, fluiddis);

    scatradis = rcp(new DRT::Discretization("scatra",reader.Comm()));
    AddDis(genprob.numscatra, scatradis);

    DRT::INPUT::NodeReader nodereader(reader, "--NODE COORDS");
    nodereader.AddElementReader(rcp(new DRT::INPUT::ElementReader(fluiddis, reader, "--FLUID ELEMENTS")));
    nodereader.AddElementReader(rcp(new DRT::INPUT::ElementReader(scatradis, reader, "--TRANSPORT ELEMENTS")));
    nodereader.Read();
    break;
  }
  case prb_fluid_ale:
  case prb_freesurf:
  {
    // allocate and input general old stuff....
    if (genprob.numfld!=2) dserror("numfld != 2 for fluid problem on ale");

    std::string distype = ptype.get<std::string>("SHAPEFCT");

    if(distype == "Nurbs")
    {
      fluiddis  = rcp(new DRT::NURBS::NurbsDiscretization("fluid"    ,reader.Comm()));
      aledis    = rcp(new DRT::NURBS::NurbsDiscretization("ale"      ,reader.Comm()));
    }
    else
    {
      fluiddis  = rcp(new DRT::Discretization("fluid"    ,reader.Comm()));
      aledis    = rcp(new DRT::Discretization("ale"      ,reader.Comm()));
    }

    AddDis(genprob.numff, fluiddis);
    AddDis(genprob.numaf, aledis);

    DRT::INPUT::NodeReader nodereader(reader, "--NODE COORDS");

    nodereader.AddElementReader(rcp(new DRT::INPUT::ElementReader(fluiddis, reader, "--FLUID ELEMENTS")));
    nodereader.AddElementReader(rcp(new DRT::INPUT::ElementReader(aledis, reader, "--ALE ELEMENTS")));

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

    std::string distype = ptype.get<std::string>("SHAPEFCT");

    if(distype == "Nurbs")
    {
      structdis = rcp(new DRT::NURBS::NurbsDiscretization("structure",reader.Comm()));
    }
    else
    {
      structdis = rcp(new DRT::Discretization("structure",reader.Comm()));
    }

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

    // read macroscale fields from main inputfile

    structdis_macro = rcp(new DRT::Discretization("structure",reader.Comm()));
    AddDis(genprob.numsf, structdis_macro);

    DRT::INPUT::NodeReader nodereader(reader, "--NODE COORDS");
    nodereader.AddElementReader(rcp(new DRT::INPUT::ElementReader(structdis_macro, reader, "--STRUCTURE ELEMENTS")));
    nodereader.Read();

    // read microscale fields from second, third, ... inputfile

    for (std::map<int,Teuchos::RCP<MAT::PAR::Material> >::const_iterator i=materials_->Map()->begin();
         i!=materials_->Map()->end();
         ++i)
    {
      int matid = i->first;
      Teuchos::RCP<MAT::PAR::Material> material = i->second;
      if (material->Type() == INPAR::MAT::m_struct_multiscale)
      {
        Teuchos::RCP<MAT::Material> mat = MAT::Material::Factory(matid);
        MAT::MicroMaterial* micromat = static_cast<MAT::MicroMaterial*>(mat.get());
        int microdisnum = micromat->MicroDisNum();

        RCP<DRT::Problem> micro_problem = DRT::Problem::Instance(microdisnum);
#ifdef PARALLEL
        RCP<Epetra_MpiComm> serialcomm = rcp(new Epetra_MpiComm(MPI_COMM_SELF));
#else
        RCP<Epetra_SerialComm> serialcomm = rcp(new Epetra_SerialComm());
#endif

        string micro_inputfile_name = micromat->MicroInputFileName();

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

        structdis_micro = rcp(new DRT::Discretization("structure", micro_reader.Comm()));
        micro_problem->AddDis(genprob.numsf, structdis_micro);

        micro_problem->ReadParameter(micro_reader);

        /* input of not mesh or time based problem data  */
        micro_problem->InputControl();

        // read materials of microscale
        // CAUTION: materials for microscale can not be read until
        // micro_reader is activated, since else materials will again be
        // read from macroscale inputfile. Besides, materials MUST be read
        // before elements are read since elements establish a connection
        // to the corresponding material! Thus do not change position of
        // function calls!
        materials_->SetReadFromProblem(microdisnum);

        micro_problem->ReadMaterials(micro_reader);

        DRT::INPUT::NodeReader micronodereader(micro_reader, "--NODE COORDS");
        micronodereader.AddElementReader(rcp(new DRT::INPUT::ElementReader(structdis_micro, micro_reader, "--STRUCTURE ELEMENTS")));
        micronodereader.Read();

        // read conditions of microscale
        // -> note that no time curves and spatial functions can be read!

        micro_problem->ReadConditions(micro_reader);

        // At this point, everything for the microscale is read,
        // subsequent reading is only for macroscale
        structdis_micro->FillComplete();
      }
    }

    // reactivate reader of macroscale as well as macroscale material

    reader.Activate();
//    ActivateMaterial();
    materials_->ResetReadFromProblem();

    break;
  } // end of else if (genprob.probtyp==prb_struct_multi)

  case prb_loma:
  {
    // allocate and input general old stuff....
    if (genprob.numfld!=2) dserror("numfld != 2 for low-Mach-number flow problem");

    // create empty discretizations
    fluiddis = rcp(new DRT::Discretization("fluid",reader.Comm()));
    AddDis(genprob.numff, fluiddis);

    scatradis = rcp(new DRT::Discretization("scatra",reader.Comm()));
    AddDis(genprob.numscatra, scatradis);

    DRT::INPUT::NodeReader nodereader(reader, "--NODE COORDS");
    nodereader.AddElementReader(rcp(new DRT::INPUT::ElementReader(fluiddis, reader, "--FLUID ELEMENTS")));
    nodereader.AddElementReader(rcp(new DRT::INPUT::ElementReader(scatradis, reader, "--TRANSPORT ELEMENTS")));
    nodereader.Read();

    break;
  } // end of else if (genprob.probtyp==prb_loma)

  case prb_elch:
  {
    // allocate and input general old stuff....
    if (genprob.numfld>3) dserror("numfld != 2 for Electrochemistry problem");

    // create empty discretizations
    fluiddis = rcp(new DRT::Discretization("fluid",reader.Comm()));
    AddDis(genprob.numff, fluiddis);

    scatradis = rcp(new DRT::Discretization("scatra",reader.Comm()));
    AddDis(genprob.numscatra, scatradis);

    aledis = rcp(new DRT::Discretization("ale",reader.Comm()));
    AddDis(genprob.numaf, aledis);

    DRT::INPUT::NodeReader nodereader(reader, "--NODE COORDS");
    nodereader.AddElementReader(rcp(new DRT::INPUT::ElementReader(fluiddis, reader, "--FLUID ELEMENTS")));
    nodereader.AddElementReader(rcp(new DRT::INPUT::ElementReader(scatradis,reader, "--TRANSPORT ELEMENTS")));
    nodereader.AddElementReader(rcp(new DRT::INPUT::ElementReader(aledis,   reader, "--ALE ELEMENTS")));
    nodereader.Read();

    break;
  } // end of else if (genprob.probtyp==prb_elch)

  case prb_combust:
  {
    // allocate and input general old stuff....
    if (genprob.numfld!=2) dserror("numfld != 2 for combustion problem");

    // create empty discretizations
    fluiddis = rcp(new DRT::Discretization("fluid",reader.Comm()));
    AddDis(genprob.numff, fluiddis);

    scatradis = rcp(new DRT::Discretization("scatra",reader.Comm()));
    AddDis(genprob.numscatra, scatradis);

    DRT::INPUT::NodeReader nodereader(reader, "--NODE COORDS");
    nodereader.AddElementReader(rcp(new DRT::INPUT::ElementReader(fluiddis, reader, "--FLUID ELEMENTS")));
    nodereader.Read();

    break;
  } // end of else if (genprob.probtyp==prb_combust)

  case prb_art_net: // _1D_ARTERY_
  {
    // allocate and input general old stuff....
    if (genprob.numfld!=1) dserror("numfld != 1 for arterial network problem");

    // create empty discretizations
    arterydis = rcp(new DRT::Discretization("artery",reader.Comm()));
    AddDis(genprob.numartf, arterydis);


    DRT::INPUT::NodeReader nodereader(reader, "--NODE COORDS");
    nodereader.AddElementReader(rcp(new DRT::INPUT::ElementReader(arterydis, reader, "--ARTERY ELEMENTS")));
    nodereader.Read();

    break;
  } // end of else if (genprob.probtyp==prb_art_net)
  default:
    dserror("Type of problem unknown");
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::Problem::setParameterList(Teuchos::RCP< Teuchos::ParameterList > const &paramList)
{
  try 
  {

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
void DRT::Problem::AddDis(int fieldnum, RCP<Discretization> dis)
{
  if (fieldnum > static_cast<int>(discretizations_.size())-1)
  {
    discretizations_.resize(fieldnum+1);
  }
  discretizations_[fieldnum].push_back(dis);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::Problem::SetDis(int fieldnum, int disnum, RCP<Discretization> dis)
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



#endif
