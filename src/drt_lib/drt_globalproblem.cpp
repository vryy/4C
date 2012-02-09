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
#include "drt_singletondestruction.H"
#include "drt_globalproblem.H"
#include "drt_inputreader.H"
#include "drt_elementreader.H"
#include "drt_nodereader.H"
#include "drt_timecurve.H"
#include "drt_utils.H"
#include "../drt_mat/material.H"
#include "../drt_mat/matpar_bundle.H"
#include "../drt_inpar/drt_validconditions.H"
#include "../drt_inpar/drt_validparameters.H"
#include "../drt_inpar/drt_validmaterials.H"
#include "../drt_mat/micromaterial.H"
#include "../drt_nurbs_discret/drt_nurbs_discret.H"

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

/*!----------------------------------------------------------------------
\brief ranks and communicators

<pre>                                                         m.gee 8/00
This structure struct _PAR par; is defined in main_ccarat.c
and the type is in partition.h
</pre>

*----------------------------------------------------------------------*/
extern struct _PAR   par;


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
  // destroy singleton objects when the problem object is still alive
  for ( vector<RCP<Problem> >::iterator i=instances_.begin(); i!=instances_.end(); ++i )
  {
    Problem * p = &**i;
    for (vector<DRT::SingletonDestruction *>::iterator i=p->sds_.begin(); i!=p->sds_.end(); ++i)
    {
      DRT::SingletonDestruction * sd = *i;
      sd->Done();
    }
    p->sds_.clear();
  }

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
DRT::Problem::Problem()
{
  materials_ = Teuchos::rcp(new MAT::PAR::Bundle());
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
DRT::Problem::~Problem()
{
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
  reader.ReadGidSection("--PATIENT SPECIFIC", *list);
  reader.ReadGidSection("--STRUCTURAL DYNAMIC", *list);
  reader.ReadGidSection("--STRUCTURAL DYNAMIC/TIMEADAPTIVITY", *list);
  reader.ReadGidSection("--STRUCTURAL DYNAMIC/GENALPHA", *list);
  reader.ReadGidSection("--STRUCTURAL DYNAMIC/ONESTEPTHETA", *list);
  reader.ReadGidSection("--STRUCTURAL DYNAMIC/GEMM", *list);
  reader.ReadGidSection("--INVERSE ANALYSIS", *list);
  reader.ReadGidSection("--MULTI LEVEL MONTE CARLO", *list);
  reader.ReadGidSection("--MESHTYING AND CONTACT", *list);
  reader.ReadGidSection("--INTERACTION POTENTIAL", *list);
  reader.ReadGidSection("--FLUCTUATING HYDRODYNAMICS", *list);
  reader.ReadGidSection("--STATISTICAL MECHANICS", *list);
  reader.ReadGidSection("--THERMAL DYNAMIC", *list);
  reader.ReadGidSection("--THERMAL DYNAMIC/GENALPHA", *list);
  reader.ReadGidSection("--THERMAL DYNAMIC/ONESTEPTHETA", *list);
  reader.ReadGidSection("--TSI DYNAMIC", *list);
  reader.ReadGidSection("--POROELASTICITY DYNAMIC", *list);
  reader.ReadGidSection("--FLUID DYNAMIC", *list);
  reader.ReadGidSection("--FLUID DYNAMIC/STABILIZATION", *list);
  reader.ReadGidSection("--FLUID DYNAMIC/TURBULENCE MODEL", *list);
  reader.ReadGidSection("--FLUID DYNAMIC/SUBGRID VISCOSITY", *list);
  reader.ReadGidSection("--FLUID DYNAMIC/MULTIFRACTAL SUBGRID SCALES", *list);
  reader.ReadGidSection("--FLUID DYNAMIC/TURBULENT INFLOW", *list);
  reader.ReadGidSection("--COMBUSTION CONTROL", *list);
  reader.ReadGidSection("--COMBUSTION CONTROL/COMBUSTION FLUID", *list);
  reader.ReadGidSection("--COMBUSTION CONTROL/COMBUSTION GFUNCTION", *list);
  reader.ReadGidSection("--COMBUSTION CONTROL/COMBUSTION PDE REINITIALIZATION", *list);
  reader.ReadGidSection("--SCALAR TRANSPORT DYNAMIC", *list);
  reader.ReadGidSection("--SCALAR TRANSPORT DYNAMIC/NONLINEAR", *list);
  reader.ReadGidSection("--SCALAR TRANSPORT DYNAMIC/STABILIZATION", *list);
  //reader.ReadGidSection("--SCALAR TRANSPORT DYNAMIC/LEVELSET", *list);
  reader.ReadGidSection("--FS3I CONTROL", *list);
  reader.ReadGidSection("--ALE DYNAMIC", *list);
  reader.ReadGidSection("--FSI DYNAMIC", *list);
  reader.ReadGidSection("--FSI DYNAMIC/CONSTRAINT", *list);
  reader.ReadGidSection("--ARTERIAL DYNAMIC", *list);
  reader.ReadGidSection("--REDUCED DIMENSIONAL AIRWAYS DYNAMIC", *list);
  reader.ReadGidSection("--SEARCH TREE", *list);
  reader.ReadGidSection("--XFEM GENERAL", *list);
  reader.ReadGidSection("--LOMA CONTROL", *list);
  reader.ReadGidSection("--ELCH CONTROL", *list);
  reader.ReadGidSection("--BIOFILM CONTROL", *list);
  reader.ReadGidSection("--TOPOLOGY OPTIMIZATION CONTROL", *list);
  reader.ReadGidSection("--TOPOLOGY OPTIMIZATION CONTROL/TOPOLOGY OPTIMIZER", *list);

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
  reader.ReadGidSection("--FLUID PRESSURE SOLVER", *list);      // TODO: remove me. use CONSTRAINT SOLVER block instead
  reader.ReadGidSection("--XFLUID PROJECTION SOLVER", *list);   // what about me? only used by Axel in xfluidimplicitintegration
  reader.ReadGidSection("--STRUCT SOLVER", *list);
  reader.ReadGidSection("--ALE SOLVER", *list);
  reader.ReadGidSection("--THERMAL SOLVER", *list);
  reader.ReadGidSection("--FLUID SCALAR TRANSPORT SOLVER", *list);
  reader.ReadGidSection("--STRUCTURE SCALAR TRANSPORT SOLVER", *list);
  reader.ReadGidSection("--SCALAR TRANSPORT ELECTRIC POTENTIAL SOLVER", *list);
  reader.ReadGidSection("--COUPLED FLUID AND SCALAR TRANSPORT SOLVER", *list);
  reader.ReadGidSection("--ARTERY NETWORK SOLVER", *list);
  reader.ReadGidSection("--REDUCED DIMENSIONAL AIRWAYS SOLVER", *list);
  reader.ReadGidSection("--TSI MONOLITHIC SOLVER", *list);
  reader.ReadGidSection("--MESHTYING SOLVER", *list);             // MESHTYING SOLVER for structure/fluid meshtying
  reader.ReadGidSection("--POROELASTICITY MONOLITHIC SOLVER", *list);
  reader.ReadGidSection("--CONTACT SOLVER", *list);               // CONTACT SOLVER for contact problems (stores all special parameters for contact preconditioner)
  reader.ReadGidSection("--CONTACT CONSTRAINT SOLVER", *list);    // only used for constraint block in a saddle point problem (for contact/meshtying)

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

  genprob.ndim  = size.get<int>("DIM");

  const Teuchos::ParameterList& type = ProblemTypeParams();

  genprob.probtyp        = DRT::INPUT::IntegralValue<PROBLEM_TYP>(type,"PROBLEMTYP");

  // If there is a restart flag on the command line, ignore the input file.
  if ( genprob.restart==0 )
  {
    genprob.restart        = type.get<int>("RESTART");
  }
  else
  {
    // If there is a non-zero restart flag on the command line, the
    // RESTART flag in the input file should be zero or have the same value!
    const int restartflaginfile = type.get<int>("RESTART");
    if ((restartflaginfile > 0) and (restartflaginfile != genprob.restart))
      dserror("Restart flags in input file and command line are non-zero and different!");
  }

  // If we have an adaptive mesh, things are totally different.
  genprob.adaptive       = DRT::INPUT::IntegralValue<int>(type,"ADAPTIVE");

  // set field numbers depending on problem type and numfld
  switch (genprob.probtyp)
  {
  case prb_fsi:
  case prb_fsi_lung:
  {
    genprob.numsf=0;
    genprob.numff=1;
    genprob.numaf=2;
    break;
  }
  case prb_gas_fsi:
  case prb_biofilm_fsi:
  case prb_thermo_fsi:
  {
    genprob.numsf=0;
    genprob.numff=1;
    genprob.numaf=2;
    genprob.numscatra=3;
    break;
  }
  case prb_fsi_xfem:
  case prb_fluid_xfem:
  case prb_fluid_xfem2:
  {
    genprob.numsf=0;
    genprob.numff=1;
    genprob.numaf=2;
    break;
  }
  case prb_fluid_fluid_ale:
  {
    genprob.numff=0;
    genprob.numaf=1;
    break;
  }
  case prb_fluid_fluid_fsi:
  {
    genprob.numff=0;
    genprob.numaf=1;
    genprob.numsf=2;
    break;
  }
  case prb_fluid_fluid:
  {
    genprob.numff=0;
    break;
  }
  case prb_fluid:
  {
    genprob.numff=0;
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
  case prb_tsi:
  case prb_tfsi_aero:
  {
    genprob.numsf = 0;  /* structural field index */
    genprob.numtf = 1;  /* thermal field index */
    break;
  }
  case prb_thermo:
  {
    genprob.numtf = 0;  /* thermal field index */
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
  case prb_red_airways:
  {
    genprob.numawf = 5;  /* reduced airway network field index */
    break;
  }
  case prb_struct_ale:
  {
    genprob.numsf = 0;  /* structural field index */
    genprob.numaf = 1;  /* ale field index */
    break;
  }
  case prb_fluid_topopt:
  {
    genprob.numff = 0; /* fluid field index */
    break;
  }
  case prb_poroelast:
  {
    genprob.numsf=0;  /* structural field index */
    genprob.numff=1;  /* fluid field index */
    break;
  }
  default:
    dserror("problem type %d unknown", genprob.probtyp);
  }

  // set field ARTNET and RED_AIRWAY numbers
  // this is the numbering used for such fields coupled to higher dimensional fields
  switch (genprob.probtyp)
  {
  case prb_fsi:
  case prb_fsi_lung:
  case prb_fluid_ale:
  case prb_fluid:
  case prb_scatra:
  case prb_loma:
  case prb_elch:
  {
#ifdef D_ARTNET
    //    genprob.numartf = 3;
    genprob.numartf = 5;
#endif
#ifdef D_RED_AIRWAYS
    //    genprob.numawf  = 4;
    genprob.numawf  = 6;
#endif
  }
  default:
    break;
  }

}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::Problem::ReadMaterials(DRT::INPUT::DatFileReader& reader)
{
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
void DRT::Problem::ReadClonedMaterials(DRT::INPUT::DatFileReader& reader)
{
  const std::string name = "--CLONING MATERIAL MAP";
  std::vector<const char*> section = reader.Section(name);
  if (section.size() > 0)
  {
    for (std::vector<const char*>::iterator i=section.begin();
         i!=section.end();
         ++i)
    {
      Teuchos::RCP<std::stringstream> condline = Teuchos::rcp(new std::stringstream(*i));

      std::string field_sep1;
      std::string src_field;
      std::string num_sep1;
      std::string src_num;
      std::string field_sep2;
      std::string tar_field;
      std::string num_sep2;
      std::string tar_num;

      (*condline) >> field_sep1 >> src_field >> num_sep1 >> src_num >>
        field_sep2 >> tar_field >> num_sep2 >> tar_num;
      if ( (not (*condline)) or (field_sep1 != "SRC_FIELD" and
                                 field_sep2 != "TAR_FIELD" and
                                 num_sep1 != "SRC_MAT" and
                                 num_sep2 != "TAR_MAT") )
        dserror("invalid material line in '%s'",name.c_str());

      std::pair<string,string> fields(src_field,tar_field);

      // extract material IDs
      int src_matid = -1;
      int tar_matid = -1;
      {
        char* src_ptr;
        src_matid = strtol(src_num.c_str(),&src_ptr,10);
        if (src_ptr == src_num.c_str())
          dserror("failed to read material object number '%s'",
                  src_num.c_str());

        char* tar_ptr;
        tar_matid = strtol(tar_num.c_str(),&tar_ptr,10);
        if (tar_ptr == tar_num.c_str())
          dserror("failed to read material object number '%s'",
                  tar_num.c_str());
      }

      // processed?
      if (materials_->Find(src_matid) == -1)
        dserror("Source material 'MAT %d' could not be identified", src_matid);
      if (materials_->Find(tar_matid) == -1)
        dserror("Target material 'MAT %d' could not be identified", src_matid);

      pair<int,int> matmap(src_matid,tar_matid);
      clonefieldmatmap_[fields].insert(matmap);
    }
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::Problem::ReadConditions(DRT::INPUT::DatFileReader& reader)
{
  Epetra_Time time(*reader.Comm());
  if (reader.Comm()->MyPID()==0)
  {
    std::cout << "Read conditions                          in....";
    std::cout.flush();
  }

  /*---------------------------------------------- input of time curves */
  timecurvemanager_.ReadInput(reader);
  /*---------------------------------------- input of spatial functions */
  functionmanager_.ReadInput(reader);
  /*-------------------------------------- input of result descriptions */
  resulttest_.ReadInput(reader);
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
          {
	    // Insert a copy since we might insert the same condition in many discretizations.
            actdis->SetCondition(condlist[c]->Name(),Teuchos::rcp(new Condition(*curr->second)));
          }
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
void DRT::Problem::ReadKnots(DRT::INPUT::DatFileReader& reader)
{
  // get information on the spatial approximation --- we only read knots
  // in the nurbs case
  std::string distype = SpatialApproximation();

  // get problem dimension
  const Teuchos::ParameterList& psize= ProblemSizeParams();
  int dim = psize.get<int>("DIM");

  // Iterate through all discretizations and sort the appropriate condition
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

        if (nurbsdis==NULL)
          dserror("Discretization %s is not a NurbsDiscretization! Panic.", actdis->Name().c_str());

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
        // do not initialize elements, since this would require knot
        // vector values
        if(!actdis->Filled())
        {
          actdis->FillComplete(false,false,false);
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
    inputcontrol_ = Teuchos::rcp(new IO::InputControl(prefix, comm));

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
void DRT::Problem::WriteInputParameters()
{
  std::string s = OutputControlFile()->FileName();
  s.append(".parameter");
  std::ofstream stream(s.c_str());
  DRT::INPUT::PrintDatHeader(stream,*getParameterList(),"",false,false);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::Problem::ReadFields(DRT::INPUT::DatFileReader& reader, const bool readmesh)
{
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
  RCP<DRT::Discretization> xfluiddis       = null;
  RCP<DRT::Discretization> aledis          = null;
  RCP<DRT::Discretization> boundarydis     = null;
  RCP<DRT::Discretization> thermdis        = null;
  RCP<DRT::Discretization> scatradis       = null;
  RCP<DRT::Discretization> fluidscatradis  = null;
  RCP<DRT::Discretization> structscatradis = null;
  RCP<DRT::Discretization> arterydis       = null; //_1D_ARTERY_
  RCP<DRT::Discretization> airwaydis       = null; //
  RCP<DRT::Discretization> optidis         = null;

  // decide which kind of spatial representation is required
  std::string distype = SpatialApproximation();

  // the basic node reader. now add desired element readers to it!
  DRT::INPUT::NodeReader nodereader(reader, "--NODE COORDS");

  switch (genprob.probtyp)
  {
  case prb_fsi:
  case prb_fsi_lung:
  {
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
      xfluiddis = rcp(new DRT::Discretization("xfluid"   ,reader.Comm()));
      aledis    = rcp(new DRT::Discretization("ale"      ,reader.Comm()));
    }

    AddDis(genprob.numsf, structdis);
    AddDis(genprob.numff, fluiddis);
    if (xfluiddis!=Teuchos::null)
      AddDis(genprob.numff, xfluiddis); // xfem discretization on slot 1
    AddDis(genprob.numaf, aledis);

    std::set<std::string> fluidelementtypes;
    fluidelementtypes.insert("FLUID");
    fluidelementtypes.insert("FLUID2");
    fluidelementtypes.insert("FLUID3");

    nodereader.AddElementReader(rcp(new DRT::INPUT::ElementReader(structdis, reader, "--STRUCTURE ELEMENTS")));
    nodereader.AddElementReader(rcp(new DRT::INPUT::ElementReader(fluiddis, reader, "--FLUID ELEMENTS", fluidelementtypes)));
    if (xfluiddis!=Teuchos::null)
      nodereader.AddElementReader(rcp(new DRT::INPUT::ElementReader(xfluiddis, reader, "--FLUID ELEMENTS", "XFLUID3")));
    nodereader.AddElementReader(rcp(new DRT::INPUT::ElementReader(aledis, reader, "--ALE ELEMENTS")));

    break;
  }
  case prb_fluid_fluid_ale:
  {
    fluiddis  = rcp(new DRT::Discretization("fluid"    ,reader.Comm()));
    xfluiddis = rcp(new DRT::Discretization("xfluid"   ,reader.Comm()));
    aledis    = rcp(new DRT::Discretization("ale"      ,reader.Comm()));

    AddDis(genprob.numff, fluiddis);
    if (xfluiddis!=Teuchos::null)
      AddDis(genprob.numff, xfluiddis); // xfem discretization on slot 1
    AddDis(genprob.numaf, aledis);

    nodereader.AddElementReader(rcp(new DRT::INPUT::ElementReader(fluiddis, reader, "--FLUID ELEMENTS", "FLUID3")));
    nodereader.AddElementReader(rcp(new DRT::INPUT::ElementReader(xfluiddis, reader, "--FLUID ELEMENTS", "FLUID3")));
    nodereader.AddElementReader(rcp(new DRT::INPUT::ElementReader(aledis, reader, "--ALE ELEMENTS")));

    break;
  }
  case prb_fluid_fluid_fsi:
  {
    fluiddis  = rcp(new DRT::Discretization("fluid"    ,reader.Comm()));
    xfluiddis = rcp(new DRT::Discretization("xfluid"   ,reader.Comm()));
    aledis    = rcp(new DRT::Discretization("ale"      ,reader.Comm()));
    structdis = rcp(new DRT::Discretization("structure",reader.Comm()));

    AddDis(genprob.numff, fluiddis);
    if (xfluiddis!=Teuchos::null)
      AddDis(genprob.numff, xfluiddis); // xfem discretization on slot 1
    AddDis(genprob.numaf, aledis);
    AddDis(genprob.numsf, structdis);

    nodereader.AddElementReader(rcp(new DRT::INPUT::ElementReader(fluiddis, reader, "--FLUID ELEMENTS", "FLUID3")));
    nodereader.AddElementReader(rcp(new DRT::INPUT::ElementReader(xfluiddis, reader, "--FLUID ELEMENTS", "FLUID3")));
    nodereader.AddElementReader(rcp(new DRT::INPUT::ElementReader(structdis, reader, "--STRUCTURE ELEMENTS")));
    nodereader.AddElementReader(rcp(new DRT::INPUT::ElementReader(aledis, reader, "--ALE ELEMENTS")));

    break;
  }
  case prb_fluid_fluid:
  {
    fluiddis  = rcp(new DRT::Discretization("fluid"    ,reader.Comm()));
    xfluiddis = rcp(new DRT::Discretization("xfluid"   ,reader.Comm()));

    AddDis(genprob.numff, fluiddis);
    if (xfluiddis!=Teuchos::null)
      AddDis(genprob.numff, xfluiddis); // xfem discretization on slot 1

    nodereader.AddElementReader(rcp(new DRT::INPUT::ElementReader(fluiddis, reader, "--FLUID ELEMENTS", "FLUID3")));
    nodereader.AddElementReader(rcp(new DRT::INPUT::ElementReader(xfluiddis, reader, "--FLUID ELEMENTS", "FLUID3")));

    break;
  }
  case prb_gas_fsi:
  case prb_biofilm_fsi:
  case prb_thermo_fsi:
  {
    if(distype == "Nurbs")
    {
      dserror("Nurbs discretization not possible for lung gas exchange!");
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

    std::set<std::string> fluidelementtypes;
    fluidelementtypes.insert("FLUID");
    fluidelementtypes.insert("FLUID2");
    fluidelementtypes.insert("FLUID3");

    nodereader.AddElementReader(rcp(new DRT::INPUT::ElementReader(structdis, reader, "--STRUCTURE ELEMENTS")));
    nodereader.AddElementReader(rcp(new DRT::INPUT::ElementReader(fluiddis, reader, "--FLUID ELEMENTS", fluidelementtypes)));
    nodereader.AddElementReader(rcp(new DRT::INPUT::ElementReader(aledis, reader, "--ALE ELEMENTS")));

#ifdef EXTENDEDPARALLELOVERLAP
    structdis->CreateExtendedOverlap(false,false,false);
#endif

    // fluid scatra field
    fluidscatradis = rcp(new DRT::Discretization("scatra1",reader.Comm()));
    AddDis(genprob.numscatra, fluidscatradis);

    // structure scatra field
    structscatradis = rcp(new DRT::Discretization("scatra2",reader.Comm()));
    AddDis(genprob.numscatra, structscatradis);

    break;
  }
  case prb_fsi_xfem:
  case prb_fluid_xfem:
  case prb_fluid_xfem2:
  {
    structdis = rcp(new DRT::Discretization("structure",reader.Comm()));
    fluiddis  = rcp(new DRT::Discretization("fluid"    ,reader.Comm()));
    aledis    = rcp(new DRT::Discretization("ale"      ,reader.Comm()));

    AddDis(genprob.numsf, structdis);
    AddDis(genprob.numff, fluiddis);
    AddDis(genprob.numaf, aledis);

    nodereader.AddElementReader(rcp(new DRT::INPUT::ElementReader(structdis, reader, "--STRUCTURE ELEMENTS")));
    nodereader.AddElementReader(rcp(new DRT::INPUT::ElementReader(fluiddis, reader, "--FLUID ELEMENTS")));
    nodereader.AddElementReader(rcp(new DRT::INPUT::ElementReader(aledis, reader, "--ALE ELEMENTS")));

    break;
  }
  case prb_ale:
  {
    if(distype == "Nurbs")
    {
      aledis = rcp(new DRT::NURBS::NurbsDiscretization("ale",reader.Comm()));
    }
    else
    {
      aledis = rcp(new DRT::Discretization("ale",reader.Comm()));
    }

    AddDis(genprob.numaf, aledis);

    nodereader.AddElementReader(rcp(new DRT::INPUT::ElementReader(aledis, reader, "--ALE ELEMENTS")));

    break;
  }
  case prb_fluid:
  {
    if(distype == "Nurbs")
    {
      fluiddis = rcp(new DRT::NURBS::NurbsDiscretization("fluid",reader.Comm()));
    }
    else
    {
      fluiddis  = rcp(new DRT::Discretization("fluid",reader.Comm()));
      xfluiddis = rcp(new DRT::Discretization("xfluid",reader.Comm()));
    }

    AddDis(genprob.numff, fluiddis);
    if ( xfluiddis!=Teuchos::null )
      AddDis(genprob.numff, xfluiddis); // xfem discretization on slot 1

    std::set<std::string> fluidelementtypes;
    fluidelementtypes.insert("FLUID");
    fluidelementtypes.insert("FLUID2");
    fluidelementtypes.insert("FLUID3");

    nodereader.AddElementReader(rcp(new DRT::INPUT::ElementReader(fluiddis, reader, "--FLUID ELEMENTS", fluidelementtypes)));

    if (xfluiddis!=Teuchos::null)
      nodereader.AddElementReader(rcp(new DRT::INPUT::ElementReader(xfluiddis, reader, "--FLUID ELEMENTS", "XFLUID3")));

    break;
  }
  case prb_scatra:
  {
    // create empty discretizations
    if(distype == "Nurbs")
    {
      fluiddis = rcp(new DRT::NURBS::NurbsDiscretization("fluid",reader.Comm()));
      scatradis = rcp(new DRT::NURBS::NurbsDiscretization("scatra",reader.Comm()));
    }
    else
    {
      fluiddis = rcp(new DRT::Discretization("fluid",reader.Comm()));
      scatradis = rcp(new DRT::Discretization("scatra",reader.Comm()));
    }

    AddDis(genprob.numff, fluiddis);
    AddDis(genprob.numscatra, scatradis);

    nodereader.AddElementReader(rcp(new DRT::INPUT::ElementReader(fluiddis, reader, "--FLUID ELEMENTS")));
    nodereader.AddElementReader(rcp(new DRT::INPUT::ElementReader(scatradis, reader, "--TRANSPORT ELEMENTS")));

    break;
  }
  case prb_fluid_ale:
  case prb_freesurf:
  {
    if(distype == "Nurbs")
    {
      fluiddis  = rcp(new DRT::NURBS::NurbsDiscretization("fluid"    ,reader.Comm()));
      aledis    = rcp(new DRT::NURBS::NurbsDiscretization("ale"      ,reader.Comm()));
    }
    else
    {
      fluiddis  = rcp(new DRT::Discretization("fluid"    ,reader.Comm()));
      xfluiddis = rcp(new DRT::Discretization("xfluid"   ,reader.Comm()));
      aledis    = rcp(new DRT::Discretization("ale"      ,reader.Comm()));
    }

    AddDis(genprob.numff, fluiddis);
    if (xfluiddis!=Teuchos::null)
      AddDis(genprob.numff, xfluiddis); // xfem discretization on slot 1
    AddDis(genprob.numaf, aledis);

    std::set<std::string> fluidelementtypes;
    fluidelementtypes.insert("FLUID");
    fluidelementtypes.insert("FLUID2");
    fluidelementtypes.insert("FLUID3");

    nodereader.AddElementReader(rcp(new DRT::INPUT::ElementReader(fluiddis, reader, "--FLUID ELEMENTS",fluidelementtypes)));
    if (xfluiddis!=Teuchos::null)
      nodereader.AddElementReader(rcp(new DRT::INPUT::ElementReader(xfluiddis, reader, "--FLUID ELEMENTS", "XFLUID3")));
    nodereader.AddElementReader(rcp(new DRT::INPUT::ElementReader(aledis, reader, "--ALE ELEMENTS")));

    break;
  }
  case prb_fluid_pm:
  {
    dserror("prb_fluid_pm not supported by BACI. just use projection method as fluid-solver for fluid problems!");
    break;
  }
  case prb_tsi:
  case prb_tfsi_aero:
  {
    structdis = Teuchos::rcp(new DRT::Discretization("structure",reader.Comm()));
    thermdis  = Teuchos::rcp(new DRT::Discretization("thermo"   ,reader.Comm()));

    AddDis(genprob.numsf, structdis);
    AddDis(genprob.numtf, thermdis);

    nodereader.AddElementReader(Teuchos::rcp(new DRT::INPUT::ElementReader(structdis, reader, "--STRUCTURE ELEMENTS")));
//    nodereader.AddElementReader(Teuchos::rcp(new DRT::INPUT::ElementReader(thermdis, reader, "--THERMO ELEMENTS")));

    break;
  }
  case prb_thermo:
  {
    {
      thermdis = Teuchos::rcp(new DRT::Discretization("thermo",reader.Comm()));
      AddDis(genprob.numtf, thermdis);

      nodereader.AddElementReader(Teuchos::rcp(new DRT::INPUT::ElementReader(thermdis, reader, "--THERMO ELEMENTS")));

      break;
    } // end of else if (genprob.probtyp==prb_thermo)
  }

  case prb_structure:
  {
    if(distype == "Nurbs")
    {
      structdis = rcp(new DRT::NURBS::NurbsDiscretization("structure",reader.Comm()));
    }
    else
    {
      structdis = rcp(new DRT::Discretization("structure",reader.Comm()));
    }

    AddDis(genprob.numsf, structdis);

    nodereader.AddElementReader(rcp(new DRT::INPUT::ElementReader(structdis, reader, "--STRUCTURE ELEMENTS")));

    break;
  } // end of else if (genprob.probtyp==prb_structure)

  case prb_loma:
  {
    // create empty discretizations
    fluiddis = rcp(new DRT::Discretization("fluid",reader.Comm()));
    AddDis(genprob.numff, fluiddis);

    scatradis = rcp(new DRT::Discretization("scatra",reader.Comm()));
    AddDis(genprob.numscatra, scatradis);

    nodereader.AddElementReader(rcp(new DRT::INPUT::ElementReader(fluiddis, reader, "--FLUID ELEMENTS")));
    nodereader.AddElementReader(rcp(new DRT::INPUT::ElementReader(scatradis, reader, "--TRANSPORT ELEMENTS")));

    break;
  } // end of else if (genprob.probtyp==prb_loma)

  case prb_elch:
  {
    // create empty discretizations
    if(distype == "Nurbs")
    {
      fluiddis  = rcp(new DRT::NURBS::NurbsDiscretization("fluid",reader.Comm()));
      scatradis = rcp(new DRT::NURBS::NurbsDiscretization("scatra",reader.Comm()));
      aledis    = rcp(new DRT::NURBS::NurbsDiscretization("ale",reader.Comm()));
    }
    else
    {
      fluiddis  = rcp(new DRT::Discretization("fluid",reader.Comm()));
      scatradis = rcp(new DRT::Discretization("scatra",reader.Comm()));
      aledis    = rcp(new DRT::Discretization("ale",reader.Comm()));
    }

    AddDis(genprob.numff, fluiddis);
    AddDis(genprob.numscatra, scatradis);
    AddDis(genprob.numaf, aledis);

    nodereader.AddElementReader(rcp(new DRT::INPUT::ElementReader(fluiddis, reader, "--FLUID ELEMENTS")));
    nodereader.AddElementReader(rcp(new DRT::INPUT::ElementReader(scatradis,reader, "--TRANSPORT ELEMENTS")));
    nodereader.AddElementReader(rcp(new DRT::INPUT::ElementReader(aledis,   reader, "--ALE ELEMENTS")));

    break;
  } // end of else if (genprob.probtyp==prb_elch)

  case prb_combust:
  {
    // create empty discretizations
    fluiddis = rcp(new DRT::Discretization("fluid",reader.Comm()));
    AddDis(genprob.numff, fluiddis);

    scatradis = rcp(new DRT::Discretization("scatra",reader.Comm()));
    AddDis(genprob.numscatra, scatradis);

    nodereader.AddElementReader(rcp(new DRT::INPUT::ElementReader(fluiddis, reader, "--FLUID ELEMENTS")));

    break;
  } // end of else if (genprob.probtyp==prb_combust)

  case prb_art_net: // _1D_ARTERY_
  {
    // create empty discretizations
    arterydis = rcp(new DRT::Discretization("artery",reader.Comm()));
    AddDis(genprob.numartf, arterydis);

    nodereader.AddElementReader(rcp(new DRT::INPUT::ElementReader(arterydis, reader, "--ARTERY ELEMENTS")));

    break;
  } // end of else if (genprob.probtyp==prb_art_net)
  case prb_red_airways: // _reduced D airways
  {
    // create empty discretizations
    airwaydis = rcp(new DRT::Discretization("red_airway",reader.Comm()));
    AddDis(genprob.numawf, airwaydis);

    nodereader.AddElementReader(rcp(new DRT::INPUT::ElementReader(airwaydis, reader, "--REDUCED D AIRWAYS ELEMENTS")));

    break;
  } // end of else if (genprob.probtyp==prb_red_airways)
  case prb_struct_ale: // structure with ale
  {
    structdis = Teuchos::rcp(new DRT::Discretization("structure",reader.Comm()));
    aledis  = Teuchos::rcp(new DRT::Discretization("ale"   ,reader.Comm()));

    AddDis(genprob.numsf, structdis);
    AddDis(genprob.numaf, aledis);

    nodereader.AddElementReader(Teuchos::rcp(new DRT::INPUT::ElementReader(structdis, reader, "--STRUCTURE ELEMENTS")));

    break;
  } // end of else if (genprob.probtyp==prb_struct_ale)
  case prb_fluid_topopt:
  {
    // create empty discretizations
    fluiddis = rcp(new DRT::Discretization("fluid",reader.Comm()));
    AddDis(genprob.numff, fluiddis);

    nodereader.AddElementReader(rcp(new DRT::INPUT::ElementReader(fluiddis, reader, "--FLUID ELEMENTS")));

    break;
  } // end of else if (genprob.probtyp==prb_combust)
  case prb_poroelast:
  {
    // create empty discretizations
	  structdis = Teuchos::rcp(new DRT::Discretization("structure",reader.Comm()));
	  fluiddis  = Teuchos::rcp(new DRT::Discretization("fluid"   ,reader.Comm()));

	  AddDis(genprob.numsf, structdis);
	  AddDis(genprob.numff, fluiddis);

	  nodereader.AddElementReader(Teuchos::rcp(new DRT::INPUT::ElementReader(structdis, reader, "--STRUCTURE ELEMENTS")));

	  break;
  }// end of else if (genprob.probtyp==prb_poroelast)

  default:
    dserror("Unknown problem type: %d",genprob.probtyp);
  }

  // add artery or airways discretizations only for the following problem types
  switch (genprob.probtyp)
  {
  case prb_fsi:
  case prb_fsi_lung:
  case prb_fluid_ale:
  case prb_fluid:
  case prb_scatra:
  case prb_loma:
  case prb_elch:
  {
#ifdef D_ARTNET
    if(distype != "Nurbs")
    {
      // create empty discretizations
      arterydis = rcp(new DRT::Discretization("artery",reader.Comm()));
      AddDis(genprob.numartf, arterydis);
      nodereader.AddElementReader(rcp(new DRT::INPUT::ElementReader(arterydis, reader, "--ARTERY ELEMENTS")));
    }
#endif
#ifdef D_RED_AIRWAYS
    if(distype != "Nurbs")
    {
      airwaydis = rcp(new DRT::Discretization("red_airway",reader.Comm()));
      AddDis(genprob.numawf, airwaydis);
      nodereader.AddElementReader(rcp(new DRT::INPUT::ElementReader(airwaydis, reader, "--REDUCED D AIRWAYS ELEMENTS")));
    }
#endif
  }
  default:
    break;
  }

  if (readmesh) // now read and allocate!
  {
    // we read nodes and elements for the desired fields as specified above
    nodereader.Read();

    // care for special applications
    switch (genprob.probtyp)
    {
    case prb_fsi:
    case prb_fsi_lung:
    {
      // read microscale fields from second, third, ... inputfile if necessary
      // (in case of multi-scale material models in structure field)
      ReadMicroFields(reader);
      break;
    }
    case prb_structure:
    {
      // read microscale fields from second, third, ... inputfile if necessary
      // (in case of multi-scale material models)
      ReadMicroFields(reader);

      // Read in another discretization for MultiLevel Monte Carlo use
      ReadMultiLevelDiscretization(reader);
      break;
    }
    default:
      break;
    }
  } // if(readmesh)

}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::Problem::ReadMicroFields(DRT::INPUT::DatFileReader& reader)
{
  // do we have micro materials at all?

  bool foundmicromat = false;

  for (std::map<int,Teuchos::RCP<MAT::PAR::Material> >::const_iterator i=materials_->Map()->begin();
       i!=materials_->Map()->end();
       ++i)
  {
    Teuchos::RCP<MAT::PAR::Material> material = i->second;
    if (material->Type() == INPAR::MAT::m_struct_multiscale)
    {
      foundmicromat = true;
      break;
    }
  }

  if (foundmicromat)
  {
    // make sure that we read the micro discretizations only on the processors on
    // which elements with the corresponding micro material are evaluated

    RCP<DRT::Problem> macro_problem = DRT::Problem::Instance();
    RCP<DRT::Discretization> macro_dis = macro_problem->Dis(genprob.numsf,0);

    std::set<int> my_multimat_IDs;

    // take care also of ghosted elements! -> ElementColMap!
    for (int i=0; i<macro_dis->ElementColMap()->NumMyElements(); ++i)
    {
      DRT::Element* actele = macro_dis->lColElement(i);
      RCP<MAT::Material> actmat = actele->Material();

      if (actmat->MaterialType() == INPAR::MAT::m_struct_multiscale)
      {
        MAT::PAR::Parameter* actparams = actmat->Parameter();
        my_multimat_IDs.insert(actparams->Id());
      }
    }

    for (std::map<int,Teuchos::RCP<MAT::PAR::Material> >::const_iterator i=materials_->Map()->begin();
         i!=materials_->Map()->end();
         ++i)
    {
      int matid = i->first;

      if (my_multimat_IDs.find(matid)!=my_multimat_IDs.end())
      {
        Teuchos::RCP<MAT::Material> mat = MAT::Material::Factory(matid);
        MAT::MicroMaterial* micromat = static_cast<MAT::MicroMaterial*>(mat.get());
        int microdisnum = micromat->MicroDisNum();

        RCP<DRT::Problem> micro_problem = DRT::Problem::Instance(microdisnum);

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

        DRT::INPUT::DatFileReader micro_reader(micro_inputfile_name, reader.Comm(), 1);

        RCP<DRT::Discretization> structdis_micro = rcp(new DRT::Discretization("structure", micro_reader.Comm()));
        micro_problem->AddDis(genprob.numsf, structdis_micro);

        micro_problem->ReadParameter(micro_reader);

        /* input of not mesh or time based problem data  */
        micro_problem->InputControl();

        // read materials of microscale
        // CAUTION: materials for microscale cannot be read until
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

        // set the problem number from which to call materials again to zero
        // (i.e. macro problem), cf. MAT::Material::Factory!
        materials_->ResetReadFromProblem();
      }
    }
    materials_->ResetReadFromProblem();
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::Problem::ReadMultiLevelDiscretization(DRT::INPUT::DatFileReader& reader)
{
  // check whether multilvel monte carlo is on
  const Teuchos::ParameterList& mlmcp = DRT::Problem::Instance()->MultiLevelMonteCarloParams();
  // should not read in second discretization if not needed
   if(Teuchos::getIntegralValue<int>(mlmcp,"MLMC")!= false && Teuchos::getIntegralValue<int>(mlmcp,"PROLONGATERES")!=false)
  //if(Teuchos::getIntegralValue<int>(mlmcp,"MLMC")!= false)
  {
    string second_input_file = mlmcp.get<std::string>("DISCRETIZATION_FOR_PROLONGATION");

    RCP<DRT::Problem> multilevel_problem = DRT::Problem::Instance(1);

    // Read in other level
    DRT::INPUT::DatFileReader multilevel_reader(second_input_file, reader.Comm(), 1);

    RCP<DRT::Discretization> structdis_multilevel = rcp(new DRT::Discretization("structure", multilevel_reader.Comm()));

    multilevel_problem->AddDis(genprob.numsf, structdis_multilevel);
    multilevel_problem->ReadParameter(multilevel_reader);
    /* input of not mesh or time based problem data  */
    //multilevel_problem->InputControl();
    // Read in Materials
    DRT::Problem::Instance()->materials_->SetReadFromProblem(1);
    multilevel_problem->ReadMaterials(multilevel_reader);
    // Read in Nodes and Elements
    DRT::INPUT::NodeReader multilevelnodereader(multilevel_reader, "--NODE COORDS");
    multilevelnodereader.AddElementReader(rcp(new DRT::INPUT::ElementReader(structdis_multilevel, multilevel_reader, "--STRUCTURE ELEMENTS")));
    multilevelnodereader.Read();

    // read conditions of other levels
    // -> note that no time curves and spatial functions can be read!
    multilevel_problem->ReadConditions(multilevel_reader);

    // At this point, everything for the other levels is read,
    // subsequent reading is only for level 0
    structdis_multilevel->FillComplete();

    // set the problem number from which to call materials again to zero
    materials_->SetReadFromProblem(0);

    materials_->ResetReadFromProblem();
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


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
DRT::UTILS::Function& DRT::Problem::Funct(int num)
{
  return functionmanager_.Funct(num);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
DRT::UTILS::TimeCurve& DRT::Problem::Curve(int num)
{
  return timecurvemanager_.Curve(num);
}


#endif

