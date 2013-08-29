/*!----------------------------------------------------------------------
\file statmech_manager.cpp
\brief management and auxiliary functions for statistical mechanics

<pre>
Maintainer: Kei Müller
            mueller@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15276
</pre>

*----------------------------------------------------------------------*/

#include "statmech_manager.H"

#include <Teuchos_Time.hpp>

#include "../drt_inpar/inpar_statmech.H"
#include "../drt_lib/drt_discret.H"
#include "../linalg/linalg_utils.H"
#include "../drt_lib/drt_element.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_condition_utils.H"
#include "../linalg/linalg_sparsematrix.H"
#include "../linalg/linalg_serialdensematrix.H"
#include "../drt_fem_general/largerotations.H"
#include "../drt_beamcontact/beam3contact_manager.H"
#include "../drt_beamcontact/beam3contact_octtree.H"
#include "../drt_io/io.H"
#include "../drt_io/io_control.H"

#include "../drt_beam3/beam3.H"
#include "../drt_beam3ii/beam3ii.H"
#include "../drt_beam3eb/beam3eb.H"
#include "../drt_beam3cl/beam3cl.H"
#include "../drt_truss3/truss3.H"
#include "../drt_torsion3/torsion3.H"

// defines flags for debugging and optimization purposes
//#define MEASURETIME
//#define DEBUGCOUT

/*----------------------------------------------------------------------*
 |  ctor (public)                                             cyron 09/08|
 *----------------------------------------------------------------------*/
STATMECH::StatMechManager::StatMechManager(Teuchos::RCP<DRT::Discretization> discret):
statmechparams_( DRT::Problem::Instance()->StatisticalMechanicsParams() ),
unconvergedsteps_(0),
starttimeoutput_(-1.0),
timeintervalstep_(0),
endtoendref_(0.0),
istart_(0),
basisnodes_(discret->NumGlobalNodes()),
basiselements_(discret->NumGlobalElements()),
outputfilenumber_(-1),
discret_(discret),
sumsquareincpar_(0.0),
sumsquareincort_(0.0),
sumrotmiddle_(0.0),
sumsquareincmid_(0.0),
sumsquareincrot_(0.0),
useinitdbcset_(false)
{
  Teuchos::ParameterList parameters = DRT::Problem::Instance()->StructuralDynamicParams();

  //initialize random generators
  SeedRandomGenerators(0);

  // retrieve output root path
  BuildStatMechRootPath();

  // retrieve the dimensions of the periodic boundary box and
  // set spatial resolution for search algorithm binding spots x crosslinkers
  InitializeStatMechValues();

  // Generate fully overlapping node column map.
  // Reason: Necessary for nearest neighbor search of filament/filament and binding site/crosslinker pairs
  CreateFullyOverlappingNodeMap();

  // initial clearing of dbc management vectors
  dbcnodesets_.clear();

  // filament number conditions: quick reference which node belongs to which filament
  InitializeFilamentNumbers();

  // force sensor conditions: for viscoelastic and creep simulations, forces need to be meassured...
  InitializeForceSensors();

  // read points in time from input file where certain actions take place
  SetStartStep(parameters);

  // Initialize crosslinker positions (point-like particles while not doubly-bound)
  CrosslinkerMoleculeInit();

  return;
}// StatMechManager::StatMechManager


/*----------------------------------------------------------------------*
 | seed all random generators of this object with fixed seed if given and|
 | with system time otherwise; seedparameter is used only in the first   |
 | case to calculate the actual seed variable based on some given fixed  |
 | seed value; note that seedparameter may be any integer, but has to be |
 | been set in a deterministic way so that it for a certain call of this |
 | method at a certain point in the program always the same number       |
 | whenever the program is used                               cyron 11/10|
 *----------------------------------------------------------------------*/
void STATMECH::StatMechManager::SeedRandomGenerators(const int seedparameter)
{
  //integer for seeding all random generators
  int seedvariable = 0;

  /*if input flag FIXEDSEED == YES: use same random numbers in each program start;
   *to this end compute seedvariable from given parameter FIXEDSEED and some other
   *deterministic parameter seedparameter given to this method at runtime*/
  if(DRT::INPUT::IntegralValue<int>(statmechparams_,"FIXEDSEED"))
  {
    seedvariable = (statmechparams_.get<int>("INITIALSEED", 0) + seedparameter)*(discret_->Comm().MyPID() + 1);

    randomnumbergen_.seed((unsigned int)seedvariable);
  }
   /*else set seed according to system time and different for each processor
   *(pseudo-random seed) if seedparameter == 0 (this allows for conveniently
   *using a random seed only at certain points in the program, e.g. only once
   *in the beginning; one has just to make sure that seedparameter == 0 does
   *not happen at any other point in the program*/
  else if(seedparameter == 0)
  {
    seedvariable = time(0)*(discret_->Comm().MyPID() + 1);

    randomnumbergen_.seed((unsigned int)seedvariable);
  }

#if (BOOST_MAJOR_VERSION == 1) && (BOOST_MINOR_VERSION >= 47)
  uniformgen_ = Teuchos::rcp(new boost::uniform_01<randnumgen&>(randomnumbergen_));
#else
  boost::uniform_01<>           uniformdist;
  uniformgen_ = Teuchos::rcp(new boost::variate_generator<randnumgen&,boost::uniform_01<> >(randomnumbergen_,uniformdist));
#endif
  boost::normal_distribution<>  normaldist(0.0,1.0);
  normalgen_ = Teuchos::rcp(new boost::variate_generator<randnumgen&,boost::normal_distribution<> >(randomnumbergen_,normaldist));

  return;
}

/*----------------------------------------------------------------------*
 | Set Period Length and Search Resolution       mueller (public)  04/12|
 *----------------------------------------------------------------------*/
void STATMECH::StatMechManager::InitializeStatMechValues()
{
  if(!discret_->Comm().MyPID())
    std::cout<<"====== StatMechManager Init ======="<<std::endl;

  // determine network type / simulation type
  switch(DRT::INPUT::IntegralValue<INPAR::STATMECH::NetworkType>(statmechparams_, "NETWORKTYPE"))
  {
    case INPAR::STATMECH::networktype_std:
    {
      networktype_ = statmech_network_std;
      if(!discret_->Comm().MyPID())
        std::cout<<"-- standard random network simulation"<<std::endl;
    }
    break;
    case INPAR::STATMECH::networktype_casimir:
    {
      networktype_ = statmech_network_casimir;
      if(!discret_->Comm().MyPID())
        std::cout<<"-- casimir-force simulation setup..."<<std::endl;
    }
    break;
    case INPAR::STATMECH::networktype_motassay:
    {
      networktype_ = statmech_network_motassay;
      if(!discret_->Comm().MyPID())
        std::cout<<"-- motility assay simulation setup..."<<std::endl;
    }
    break;
    default:
      dserror("Unknown network type %i", DRT::INPUT::IntegralValue<INPAR::STATMECH::FilamentModel>(statmechparams_, "NETWORKTYPE"));
  }

  // determine filament model
  switch(DRT::INPUT::IntegralValue<INPAR::STATMECH::FilamentModel>(statmechparams_, "FILAMENTMODEL"))
  {
    case INPAR::STATMECH::filamentmodel_std:
    {
      filamentmodel_ = statmech_filament_std;
      if(!discret_->Comm().MyPID())
        std::cout<<"-- standard filaments"<<std::endl;
    }
    break;
    case INPAR::STATMECH::filamentmodel_helical:
    {
      filamentmodel_ = statmech_filament_helical;
      if(!discret_->Comm().MyPID())
        std::cout<<"-- filaments with helical binding site orientation..."<<std::endl;
    }
    break;
    default:
      dserror("Unknown filament model %i", DRT::INPUT::IntegralValue<INPAR::STATMECH::FilamentModel>(statmechparams_, "FILAMENTMODEL"));
  }

  //determine linker model
  switch(DRT::INPUT::IntegralValue<INPAR::STATMECH::LinkerModel>(statmechparams_, "LINKERMODEL"))
  {
    case INPAR::STATMECH::linkermodel_none:
      linkermodel_ = statmech_linker_none;
      if(!discret_->Comm().MyPID())
        std::cout<<"-- no linkers"<<std::endl;
      break;
    case INPAR::STATMECH::linkermodel_std:
      linkermodel_ = statmech_linker_std;
      if(!discret_->Comm().MyPID())
        std::cout<<"-- standard linkers"<<std::endl;
      break;
    case INPAR::STATMECH::linkermodel_stdintpol:
      linkermodel_ = statmech_linker_stdintpol;
      if(!discret_->Comm().MyPID())
        std::cout<<"-- standard linkers with interpolated binding site positions"<<std::endl;
      break;
    case INPAR::STATMECH::linkermodel_bellseq:
      linkermodel_ = statmech_linker_bellseq;
      if(!discret_->Comm().MyPID())
        std::cout<<"-- linkers accounting for Bell's equation"<<std::endl;
      break;
    case INPAR::STATMECH::linkermodel_bellseqintpol:
      linkermodel_ = statmech_linker_bellseqintpol;
      if(!discret_->Comm().MyPID())
        std::cout<<"-- linkers accounting for Bell's equation (interpolated binding site positions)"<<std::endl;
      break;
    case INPAR::STATMECH::linkermodel_active:
      linkermodel_ = statmech_linker_active;
      if(!discret_->Comm().MyPID())
        std::cout<<"-- active linkers"<<std::endl;
      break;
    case INPAR::STATMECH::linkermodel_activeintpol:
      linkermodel_ = statmech_linker_activeintpol;
      if(!discret_->Comm().MyPID())
        std::cout<<"-- active linkers (interpolated binding site positions)"<<std::endl;
      break;
    default:
      dserror("Unknown linker model %i", DRT::INPUT::IntegralValue<INPAR::STATMECH::LinkerModel>(statmechparams_, "LINKERMODEL"));
  }

  periodlength_ = Teuchos::rcp(new std::vector<double>);
  periodlength_->clear();
  {
    std::istringstream PL(Teuchos::getNumericStringParameter(statmechparams_,"PERIODLENGTH"));
    std::string word;
    char* input;
    while (PL >> word)
      periodlength_->push_back(std::strtod(word.c_str(), &input));
  }

  if((int)periodlength_->size()<3)
    dserror("You only gave %d values for PERIODLENGTH! Check your input file.", (int)periodlength_->size());
  for(int i=0; i<(int)periodlength_->size(); i++)
    if(periodlength_->at(i)<0.0)
      dserror("PERIODLENGTH( %d ) = %4.2f < 0.0 does not make any sense! Check your input file.", i, periodlength_->at(i));

  // set spatial resolution for search algorithm binding spots x crosslinkers
  if(statmechparams_.get<int>("SEARCHRES",1)<1)
    dserror("Please give a plausible value (>0) for SEARCHRES!");

  // determine search resolution in each spatial direction according to the periodlength_ vector
  searchres_ = Teuchos::rcp(new std::vector<int>(3,statmechparams_.get<int>("SEARCHRES",1)));
  // in case of a non-cubic periodic volume
  if(fabs(pow(periodlength_->at(0)*periodlength_->at(1)*periodlength_->at(2), 1.0/3.0)-periodlength_->at(0))>1e-4)
  {
    double Hmax = std::max(periodlength_->at(0), periodlength_->at(1));
    Hmax = std::max(Hmax, periodlength_->at(2));
    for(int i=0; i<(int)searchres_->size(); i++)
      searchres_->at(i) = (int)(floor((periodlength_->at(i)/Hmax) * (double)(statmechparams_.get<int>("SEARCHRES",1))));
  }

  // read times for actions and the corresponding step sizes from input file
  actiontime_ = Teuchos::rcp(new std::vector<double>);
  actiontime_->clear();
  {
    std::istringstream TIME(Teuchos::getNumericStringParameter(statmechparams_,"ACTIONTIME"));
    std::string word;
    char* input;
    while (TIME >> word)
      actiontime_->push_back(std::strtod(word.c_str(), &input));
  }
  timestepsizes_ = Teuchos::rcp(new std::vector<double>);
  timestepsizes_->clear();
  {
    std::istringstream DT(Teuchos::getNumericStringParameter(statmechparams_,"ACTIONDT"));
    std::string word;
    char* input;
    while (DT >> word)
      timestepsizes_->push_back(std::strtod(word.c_str(), &input));
  }

  // sanity checks for time points and time step sizes
  Teuchos::ParameterList sdynparams = DRT::Problem::Instance()->StructuralDynamicParams();
  double iodt = sdynparams.get<double>("TIMESTEP", 1.0e9);
  double maxtime = sdynparams.get<double>("MAXTIME",-1.0);

  for(int i=0; i<(int)actiontime_->size()-1; i++)
    if(actiontime_->at(i)>actiontime_->at(i+1))
      dserror("ACTIONTIME values must be monotonously increasing!");
  if((int)actiontime_->size()!=(int)timestepsizes_->size())
    dserror("ACTIONTIME and ACTIONDT have to be equal in number in the input file!");
  if((int)actiontime_->size()<2)
  {
    if(timestepsizes_->at(0)>0.0) // unwanted stuff
      dserror("Give at least two values for ACTIONTIME AND ACTIONDT");
    else if(timestepsizes_->at(0)==0.0) //non-sensical stuff
      dserror("Given ACTIONDT is 0.0!");
    else // handled as default case from drt_validparameters.cpp and the rest of the non-sensical stuff
    {
      timestepsizes_->at(0) = iodt;
      timestepsizes_->push_back(iodt);
      actiontime_->at(0) = maxtime;
      actiontime_->push_back(maxtime);
    }
  }

  // set dbctimeindex_ (position in actiontime_ where DBCs start being applied
  dbctimeindex_ = statmechparams_.get<int>("DBCTIMEINDEX", -1);

  // DBC sanity checks
  DBCSanityCheck();

  if(!discret_->Comm().MyPID())
  {
    std::cout<<"t_equilib  = t(0) = "<<std::setprecision(10)<<actiontime_->at(0)<<" @ dt = "<<timestepsizes_->at(0)<<std::endl;
    std::cout<<"t_ktswitch = t(1) = "<<std::setprecision(10)<<actiontime_->at(1)<<" @ dt = "<<timestepsizes_->at(1)<<std::endl;
    if(dbctimeindex_>-1)
      std::cout<<"t_dbc     = t("<<dbctimeindex_<<") ="<<std::setprecision(10)<<actiontime_->at(dbctimeindex_)<<" @ dt = "<<timestepsizes_->at(dbctimeindex_)<<std::endl;
    if(actiontime_->size()>3)
      std::cout<<"other action times..."<<std::endl;
    for(int i=0; i<(int)actiontime_->size(); i++)
      if(i>1 && i!=dbctimeindex_)
        std::cout<<"t("<<i<<")        = "<<std::setprecision(10)<<actiontime_->at(i)<<" @ dt = "<<timestepsizes_->at(i)<<std::endl;
    std::cout<<"==================================="<<std::endl;
  }
  return;
}
/*----------------------------------------------------------------------*
 | Create fully overlapping node map             (private) mueller 08/13|
 *----------------------------------------------------------------------*/
void STATMECH::StatMechManager::CreateFullyOverlappingNodeMap()
{
  const Epetra_Map noderowmap = *(discret_->NodeRowMap());

  // fill my own row node ids into vector sdata
  std::vector<int> sdata(noderowmap.NumMyElements());
  for (int i=0; i<noderowmap.NumMyElements(); ++i)
    sdata[i] = noderowmap.GID(i);

  //if current processor has elements it writes its number into stproc
  std::vector<int> stproc(0);


  if (noderowmap.NumMyElements())
    stproc.push_back(discret_->Comm().MyPID());


  //information how many processors work at all
  std::vector<int> allproc(discret_->Comm().NumProc());

  //in case of n processors allproc becomes a std::vector with entries (0,1,...,n-1)
  for (int i=0; i<discret_->Comm().NumProc(); ++i) allproc[i] = i;

  //declaring new variable into which the information of stproc on all processors is gathered
  std::vector<int> rtproc(0);

  /*gathers information of stproc and writes it into rtproc; in the end rtproc is a vector which
   * contains the numbers of all processors which have elements*/
  LINALG::Gather<int>(stproc,rtproc,discret_->Comm().NumProc(),&allproc[0],discret_->Comm());

  /*in analogy to stproc and rtproc the variable rdata gathers all the element numbers which are
   * stored on different processors in their own variables sdata; thereby each processor gets
   * the information about all the nodes numbers existing in this problem*/
  std::vector<int> rdata;

  // gather all gids of nodes redundantly from sdata into rdata
  LINALG::Gather<int>(sdata,rdata,(int)rtproc.size(),&rtproc[0],discret_->Comm());

  // build completely overlapping map (on participating processors)
  Teuchos::RCP<Epetra_Map> newnodecolmap = Teuchos::rcp(new Epetra_Map(-1,(int)rdata.size(),&rdata[0],0,discret_->Comm()));
  sdata.clear();
  stproc.clear();
  rdata.clear();
  allproc.clear();
  // rtproc still in use

  //pass new fully overlapping column map to discretization
  discret_->ExportColumnNodes(*newnodecolmap);
  /*rebuild discretization based on the new column map so that each processor creates new ghost elements
   * if necessary; after the following line we have a discretization, where each processor has a fully
   * overlapping column map regardlesse of how overlapping was managed when starting BACI; having ensured
   * this allows convenient and correct (albeit not necessarily efficient) use of search algorithms and
   * crosslinkers in parallel computing*/
  discret_->FillComplete(true,false,false);
  return;
}

/*----------------------------------------------------------------------*
 | intialize filament numbers in ctor            (private) mueller 08/13|
 *----------------------------------------------------------------------*/
void STATMECH::StatMechManager::InitializeFilamentNumbers()
{
  filamentnumber_ = Teuchos::rcp( new Epetra_Vector(*(discret_->NodeColMap())) );
  filamentnumber_->PutScalar(-1.0);

  //getting a vector consisting of pointers to all filament number conditions set
  std::vector<DRT::Condition*> filamentnumberconditions(0);
  discret_->GetCondition("FilamentNumber",filamentnumberconditions);

  //next all the pointers to all the different conditions are looped
  for (int i=0; i<(int)filamentnumberconditions.size(); ++i)
  {
    //get filament number described by the current condition
    int filamentnumber = filamentnumberconditions[i]->GetInt("Filament Number");

    //get a pointer to nodal cloud covered by the current condition
    const std::vector<int>* nodeids = filamentnumberconditions[i]->Nodes();

    //loop through all the nodes of the nodal cloud
    for(int j=0; j<(int)nodeids->size() ; j++)
    {
      //get the node's global id
      int nodenumber = (*nodeids)[j];

      //turning global id into local one
      nodenumber = discret_->NodeColMap()->LID(nodenumber);

      //if the node does not belong to current processor Id is set by LID() to -1 and the node is ignored
      if(nodenumber > -1)
        (*filamentnumber_)[nodenumber] = filamentnumber;
    }
  }
  return;
}

/*----------------------------------------------------------------------*
 | initialize force sensors in ctor              (private) mueller 08/13|
 *----------------------------------------------------------------------*/
void STATMECH::StatMechManager::InitializeForceSensors()
{
  forcesensor_ = Teuchos::rcp( new Epetra_Vector(*(discret_->DofColMap())) );
  forcesensor_->PutScalar(-1.0);

  //gettin a vector consisting of pointers to all filament number conditions set
  std::vector<DRT::Condition*> forcesensorconditions(0);
  discret_->GetCondition("ForceSensor",forcesensorconditions);

  //next all the pointers to all the different conditions are looped
  for (int i=0; i<(int)forcesensorconditions.size(); ++i)
  {
    //get number of nodal dof with respect to which force is to be measured; note: numbering starts with zero
    int nodedofnumber = forcesensorconditions[i]->GetInt("DOF Number") ;

    //get a pointer to nodal cloud covered by the current condition
    const std::vector<int>* nodeids = forcesensorconditions[i]->Nodes();

    //loop through all the nodes of the nodal cloud
    for (int j=0; j<(int)nodeids->size(); j++)
    {
      //get the node's global id
      int nodenumber = (*nodeids)[j];

      //testing whether current nodedofnumber makes sense for current node
      if (nodedofnumber < 0 || nodedofnumber >= discret_->NumDof(discret_->gNode(nodenumber)))
        dserror("ForceSensor condition applied with improper local dof number");

        //global id of degree of freedom at which force is to be measured
        int dofnumber = discret_->Dof( discret_->gNode(nodenumber), nodedofnumber-1 );

        /*if the node does not belong to current processor Id is set by LID() to -1 and the node is ignored; otherwise the degrees of
         * freedom affected by this condition are marked in the vector *forcesensor_ by a one entry*/
        if(nodenumber > -1)
          (*forcesensor_)[dofnumber] = 1.0;
    }
  }
  return;
}

/*----------------------------------------------------------------------*
 | Set start step                                (private) mueller 08/13|
 *----------------------------------------------------------------------*/
void STATMECH::StatMechManager::SetStartStep(Teuchos::ParameterList& parameters)
{
  // set istart_, the number of steps after which we start writing output
  double starttimeout = statmechparams_.get<double>("STARTTIMEOUT", 0.0);
  // intermediate time intervals
  for(int i=0; i<(int)actiontime_->size(); i++)
  {
    // first time interval
    if(i==0)
    {
      if(starttimeout>=actiontime_->front())
        istart_ = (int)(round(actiontime_->front()/parameters.get<double>("delta time",0.01)));
      else
      {
        istart_ = (int)(round((starttimeout)/parameters.get<double>("delta time",0.01)));
        break;
      }
    }
    // last time interval
    else if(i==(int)actiontime_->size()-1)
    {
      if(starttimeout>=actiontime_->at(i))
        istart_ += (int)(round((actiontime_->back()-actiontime_->at(i-1))/timestepsizes_->at(i-1)) + (int)((starttimeout-actiontime_->back())/timestepsizes_->back()));
      else
        istart_ += (int)(round((starttimeout-actiontime_->at(i-1)) / timestepsizes_->at(i-1)));
    }
    // intermediate time intervals
    else
    {
      if(starttimeout>=actiontime_->at(i))
        istart_ += (int)(round((actiontime_->at(i)-actiontime_->at(i-1)) / timestepsizes_->at(i-1)));
      else
      {
        istart_ += (int)(round((starttimeout-actiontime_->at(i-1)) / timestepsizes_->at(i-1)));
        break;
      }
    }
  }
  return;
}

/*----------------------------------------------------------------------*
 | write special output for statistical mechanics (public)   cyron 09/08|
 *----------------------------------------------------------------------*/
void STATMECH::StatMechManager::Update(const int&                                      istep,
                                       const double&                                   timen,
                                       const double&                                   dt,
                                       Epetra_Vector&                                  disrow,
                                       Teuchos::RCP<LINALG::SparseOperator>&           stiff,
                                       int& ndim, Teuchos::RCP<CONTACT::Beam3cmanager> beamcmanager,
                                       bool                                            rebuildoctree,
                                       bool                                            printscreen)
{
#ifdef MEASURETIME
  const double t0 = Teuchos::Time::wallTime();
#endif // #ifdef MEASURETIME
  /* first we modify the displacement vector so that current nodal position at the end of current time step complies with
   * periodic boundary conditions, i.e. no node lies outside a cube of edge length periodlength_*/

  //if dynamic crosslinkers are used update comprises adding and deleting crosslinkers
  if (linkermodel_ != statmech_linker_none)
  {
#ifdef MEASURETIME
    const double t1 = Teuchos::Time::wallTime();
#endif

    //column maps cointaining binding spot positions and their spatial orientation
    Teuchos::RCP<Epetra_MultiVector> bspotpositions = Teuchos::rcp(new Epetra_MultiVector(*bspotcolmap_,3,true));
    Teuchos::RCP<Epetra_MultiVector> bspotrotations = Teuchos::rcp(new Epetra_MultiVector(*bspotcolmap_,3,true));

    /* note: access by ExtractMyValues requires column map vector, whereas displacements on level of time integration are
     * handled as row map vector*/
    Epetra_Vector discol(*discret_->DofColMap(), true);
    LINALG::Export(disrow, discol);
    GetBindingSpotPositions(discol, bspotpositions, bspotrotations);

#ifdef MEASURETIME
    const double t2 = Teuchos::Time::wallTime();
#endif

     // crosslink molecule diffusion
    double standarddev = sqrt(statmechparams_.get<double> ("KT", 0.0) / (2*M_PI * statmechparams_.get<double> ("ETA", 0.0) * statmechparams_.get<double> ("R_LINK", 0.0)) * dt);
    CrosslinkerDiffusion(*bspotpositions, 0.0, standarddev, dt);

#ifdef MEASURETIME
    const double t3 = Teuchos::Time::wallTime();
    double t4 = -1.0e9;
    double t5 = -1.0e9;
#endif

    // set crosslinkers, i.e. considering crosslink molecule diffusion after filaments had time to equilibrate
    if(timen>=actiontime_->front() || fabs(timen-actiontime_->front())<(dt/1e3))
    {
      if(beamcmanager!=Teuchos::null && rebuildoctree)
      {
        std::map<int, LINALG::Matrix<3, 1> > currentpositions;
        std::map<int, LINALG::Matrix<3, 1> > currentrotations;
        GetNodePositionsFromDisVec(discol, currentpositions, currentrotations, true);
        beamcmanager->OcTree()->OctTreeSearch(currentpositions);
      }
      SearchAndSetCrosslinkers(istep, timen, dt, bspotpositions, bspotrotations, discol, beamcmanager, printscreen);

#ifdef DEBUGCOUT
      for(int proc=0; proc<discret_->Comm().NumProc(); proc++)
      {
        if(proc==discret_->Comm().MyPID())
        {
          std::cout<<"\n====Proc "<<proc<<" ===="<<std::endl;
          std::cout<<"numbond: "<<std::endl;
          for(int i=0; i<numbond_->MyLength(); i++)
            if((*numbond_)[i]>0.1)
              std::cout<<"i="<<i<<": "<<(*numbond_)[i]<<std::endl;
          std::cout<<"crosslinkerbond:"<<std::endl;
          for(int i=0; i<crosslinkerbond_->MyLength(); i++)
            if((*numbond_)[i]>0.1)
              std::cout<<"i="<<i<<": "<<(*crosslinkerbond_)[0][i]<<", "<<(*crosslinkerbond_)[1][i]<<std::endl;
          std::cout<<"crosslinkerpositions:"<<std::endl;
          for(int i=0; i<crosslinkerpositions_->MyLength(); i++)
            if((*numbond_)[i]>0.1)
              std::cout<<"i="<<i<<": "<<(*crosslinkerpositions_)[0][i]<<", "<<(*crosslinkerpositions_)[1][i]<<", "<<(*crosslinkerpositions_)[2][i]<<std::endl;
        }
        discret_->Comm().Barrier();
      }
#endif

#ifdef MEASURETIME
      t4 = Teuchos::Time::wallTime();
#endif

      if(beamcmanager!=Teuchos::null && rebuildoctree)
        beamcmanager->ResetPairs();
      
      SearchAndDeleteCrosslinkers(istep, timen, dt, bspotpositions, bspotrotations, discol,beamcmanager,printscreen);

#ifdef MEASURETIME
      t5 = Teuchos::Time::wallTime();
#endif
    }

    // reset thermal energy to new value (simple value change for now, maybe Load Curve later on)
    if(fabs(timen-actiontime_->at(1))<(dt/1e3))
      statmechparams_.set("KT",statmechparams_.get<double>("KTACT",statmechparams_.get<double>("KT",0.0)));

    // // Set a certain number of double-bonded crosslinkers free
    if(fabs(timen-dt-actiontime_->back())<(dt/1e3) && statmechparams_.get<int>("REDUCECROSSLINKSBY",0)>0)
      ReduceNumOfCrosslinkersBy(statmechparams_.get<int>("REDUCECROSSLINKSBY",0), beamcmanager);

    // change length of active crosslinkers
    if(linkermodel_ == statmech_linker_active || linkermodel_ == statmech_linker_activeintpol)
      ChangeActiveLinkerLength(timen, dt, crosslinkeractlengthconv_, crosslinkeractlength_);

    /*settling administrative stuff in order to make the discretization ready for the next time step:
     * done in SearchAndeDeleteCrosslinkers():
     * synchronize the Filled() state on all processors after having added or deleted elements by CheckFilledGlobally();
     * then build new element maps and call FillComplete();
     * done here: finally Crs matrices stiff_ has to be deleted completely and made ready
     * for new assembly since their graph was changed*/
    stiff->Reset();

#ifdef MEASURETIME
    if(!discret_->Comm().MyPID())
    {
      std::cout<<"\n=================Time  Measurement================"<<std::endl;
      std::cout<<"StatMechManager::Update"<<std::endl;
      std::cout<<"Binding Spot Positions :\t"<<std::setprecision(4)<<t2-t1<<"\ts"<<std::endl;
      std::cout<<"Crosslinker Diffusion  :\t"<<std::setprecision(4)<<t3-t2<<"\ts"<<std::endl;
      std::cout<<"Set Crosslinkers       :\t"<<std::setprecision(4)<<t4-t3<<"\ts"<<std::endl;
      std::cout<<"Delete Crosslinkers    :\t"<<std::setprecision(4)<<t5-t4<<"\ts"<<std::endl;
      std::cout<<"=================================================="<<std::endl;
      std::cout<<"Update total time      :\t"<<std::setprecision(4)<<Teuchos::Time::wallTime()-t0<<"\ts"<<std::endl;
    }
#endif // #ifdef MEASURETIME
  }//linkermodel_ != statmech_linker_none
  return;
} // StatMechManager::Update()

/*----------------------------------------------------------------------*
 | Update time step size in time integration       (public)mueller 06/12|
 *----------------------------------------------------------------------*/
void STATMECH::StatMechManager::UpdateTimeAndStepSize(double& dt,
                                                      double& timeconverged,
                                                      bool    initialset)
{
  double eps = 2.0e-11;
  if(initialset)
  {
    Teuchos::ParameterList sdynparams = DRT::Problem::Instance()->StructuralDynamicParams();
    double maxtime = sdynparams.get<double>("MAXTIME",-1.0);
    timeintervalstep_ = -1;
    for(int i=0; i<(int)actiontime_->size(); i++)
    {
      if(timeconverged>=actiontime_->at(i) || fabs(timeconverged-actiontime_->at(i))<eps)
        timeintervalstep_++;
      else
      {
        if(actiontime_->at(i)==maxtime)
          timeintervalstep_ = i;
        break;
      }
    }
    dt = timestepsizes_->at(timeintervalstep_);
    // "++" needed so that after initialization, no additional update occurs
    timeintervalstep_++;
  }
  else
  {
    if(timeintervalstep_<(int)timestepsizes_->size())
    {
      double nexttimethreshold = actiontime_->at(timeintervalstep_);
      double dtnew = timestepsizes_->at(timeintervalstep_);
      do
      {
        // update time step
        double dtnew = timestepsizes_->at(timeintervalstep_);
        // update step size
        nexttimethreshold = actiontime_->at(timeintervalstep_);
        if((timeconverged>=nexttimethreshold || fabs(timeconverged-nexttimethreshold)<eps) && dtnew>0.0)
        {
          if(!discret_->Comm().MyPID())
            std::cout<<"---time step size switched from dt = "<<dt<<" to dt = "<<dtnew<<std::endl;
          dt = dtnew;
          timeintervalstep_++;
        }
      }
      while((timeconverged>=nexttimethreshold || fabs(timeconverged-nexttimethreshold)<eps) && dtnew>0.0 && timeintervalstep_<(int)actiontime_->size());
    }
  }
  return;
}

/*----------------------------------------------------------------------*
 | Retrieve nodal positions                     (public)   mueller 09/11|
 *----------------------------------------------------------------------*/
void STATMECH::StatMechManager::GetNodePositionsFromDisVec(Epetra_Vector&                      discol,
                                                           std::map<int,LINALG::Matrix<3,1> >& currentpositions,
                                                           std::map<int,LINALG::Matrix<3,1> >& currentrotations,
                                                           bool                                positionsonly)
{
  /*in preparation for later decision whether a crosslink should be established between two nodes (binding spots) we first store the
   * current positions of all column map nodes (column map binding spots) in the map currentpositions; additionally we store the rotational displacements
   * analogously in a map currentrotations for later use in setting up reference geometry of crosslinkers; the maps
   * currentpositions and currentrotations relate positions and rotations to a local column map node Id, respectively*/
  currentpositions.clear();
  if(!positionsonly)
    currentrotations.clear();

    for (int i=0; i<discret_->NumMyColNodes(); ++i)
    {

      //get pointer at a node
      const DRT::Node* node = discret_->lColNode(i);

      //get GIDs of this node's degrees of freedom
      std::vector<int> dofnode = discret_->Dof(node);

      LINALG::Matrix<3, 1> currpos;
      LINALG::Matrix<3, 1> currrot;

      currpos(0) = node->X()[0] + discol[discret_->DofColMap()->LID(dofnode[0])];
      currpos(1) = node->X()[1] + discol[discret_->DofColMap()->LID(dofnode[1])];
      currpos(2) = node->X()[2] + discol[discret_->DofColMap()->LID(dofnode[2])];
      //if node has also rotational degrees of freedom
      if (discret_->NumDof(node) == 6 && !positionsonly)
      {
        currrot(0) = discol[discret_->DofColMap()->LID(dofnode[3])];
        currrot(1) = discol[discret_->DofColMap()->LID(dofnode[4])];
        currrot(2) = discol[discret_->DofColMap()->LID(dofnode[5])];
      }

      currentpositions[node->LID()] = currpos;
      if(!positionsonly)
        currentrotations[node->LID()] = currrot;
    }
  return;
}

/*----------------------------------------------------------------------*
 | Retrieve binding spot positions              (public)   mueller 10/12|
 *----------------------------------------------------------------------*/
void STATMECH::StatMechManager::GetBindingSpotPositions(Epetra_Vector& discol, Teuchos::RCP<Epetra_MultiVector> bspotpositions, Teuchos::RCP<Epetra_MultiVector> bspotrotations)
{
  /*in preparation for later decision whether a crosslink should be established between two nodes (binding spots) we first store the
   * current positions of all column map nodes (column map binding spots) in the map currentpositions; additionally we store the rotational displacements
   * analogously in a map currentrotations for later use in setting up reference geometry of crosslinkers; the maps
   * currentpositions and currentrotations relate positions and rotations to a local column map node Id, respectively*/

  // in case the four-noded crosslinker beam element is used, currentpositions and currentrotations have to be set up another way
  if(linkermodel_ == statmech_linker_stdintpol || linkermodel_ == statmech_linker_activeintpol || linkermodel_ == statmech_linker_bellseqintpol)
    GetInterpolatedBindingSpotPositions(discol, bspotpositions, bspotrotations);
  else  // conventional crosslinker beam element, i.e. binding spots coincide with nodes
    GetNodalBindingSpotPositionsFromDisVec(discol, bspotpositions, bspotrotations);
  return;
}

/*----------------------------------------------------------------------*
 | Retrieve binding spot positions (nodes)      (public)   mueller 10/12|
 *----------------------------------------------------------------------*/
void STATMECH::StatMechManager::GetNodalBindingSpotPositionsFromDisVec(const Epetra_Vector&              discol,
                                                                       Teuchos::RCP<Epetra_MultiVector>  bspotpositions,
                                                                       Teuchos::RCP<Epetra_MultiVector>  bspotrotations)
{
  /*in preparation for later decision whether a crosslink should be established between two nodes (binding spots) we first store the
   * current positions of all column map nodes (column map binding spots) in the map bspotpositions; additionally we store the rotational displacements
   * analogously in a map bspotrotations for later use in setting up reference geometry of crosslinkers; the maps
   * bspotpositions and bspotrotations relate positions and rotations to a local column map node Id, respectively*/
  bool getpositions = true;
  bool getrotations = true;
  if(bspotpositions==Teuchos::null)
    getpositions = false;
  if(bspotrotations==Teuchos::null)
    getrotations = false;

  //sanity checks
  if(!getpositions && !getrotations)
    dserror("Both vectors are null!");
  if(getpositions)
  {
    if((bspotpositions->MyLength() != discret_->NodeColMap()->NumMyElements()) &&
        (bspotpositions->MyLength() != discret_->NodeRowMap()->NumMyElements()))
      dserror("Positions vector is neither of node column map format nor of node row map format!");
  }
  if(getrotations)
  {
    if((bspotrotations->MyLength() != discret_->NodeColMap()->NumMyElements()) &&
        (bspotrotations->MyLength() != discret_->NodeRowMap()->NumMyElements()))
      dserror("Rotations vector is neither of node column map format nor of node row map format!");
  }
  if(getpositions && getrotations)
  {
    if(bspotpositions->MyLength() != bspotrotations->MyLength())
      dserror("Positions vector and rotations vector have different lengths!");
  }

  // determine map format
  bool colmapformat = true;
  if(getpositions)
  {
    if(bspotpositions->MyLength()==discret_->NodeRowMap()->NumMyElements() && discret_->Comm().NumProc()>1)
      colmapformat = false;
  }
  else
  {
    if(bspotrotations->MyLength()==discret_->NodeRowMap()->NumMyElements() && discret_->Comm().NumProc()>1)
      colmapformat = false;
  }

  //first get triads at all row nodes (works, since here: bspotrowmap_ = noderowmap )
  Teuchos::RCP<Epetra_MultiVector> bspotpositionsrow = Teuchos::null;
  Teuchos::RCP<Epetra_MultiVector> bspotrotationsrow = Teuchos::null;
  if(bspotpositions!=Teuchos::null)
    bspotpositionsrow = Teuchos::rcp(new Epetra_MultiVector(*(discret_->NodeRowMap()),3));
  if(bspotrotations!=Teuchos::null)
    bspotrotationsrow = Teuchos::rcp(new Epetra_MultiVector(*(discret_->NodeRowMap()),3));

  //update nodaltriads_
  for (int i=0; i<discret_->NodeRowMap()->NumMyElements(); i++)
  {
    //get pointer at a node
    const DRT::Node* node = discret_->lRowNode(i);
    //get GIDs of this node's degrees of freedom
    std::vector<int> dofnode = discret_->Dof(node);

    if(getpositions)
      for(int j=0; j<bspotpositionsrow->NumVectors(); j++)
        (*bspotpositionsrow)[j][i] = node->X()[j] + discol[discret_->DofColMap()->LID(dofnode[j])];

    //if node has also rotational degrees of freedom
    if (discret_->NumDof(node) == 6 && getrotations)
      for(int j=0; j<bspotrotationsrow->NumVectors(); j++)
        (*bspotrotationsrow)[j][i] = discol[discret_->DofColMap()->LID(dofnode[j+3])];
  }

  // Export
  // column map format
  if(colmapformat)
  {
    if(getpositions)
      CommunicateMultiVector(bspotpositionsrow, bspotpositions, false, true);
    if(getrotations)
      CommunicateMultiVector(bspotrotationsrow, bspotrotations, false, true);
  }
  //row map : Why does "bspotrotations = Teuchos::rcp(new Epetra_MultiVector(*bspotrotationsrow));" only work until the end of the function?
  else
  {
    if(getpositions)
      //bspotpositions = Teuchos::rcp(new Epetra_MultiVector(*bspotpositionsrow),true);
      for(int i=0; i<bspotpositions->NumVectors(); i++)
        for(int j=0; j<bspotpositions->MyLength(); j++)
          (*bspotpositions)[i][j] = (*bspotpositionsrow)[i][j];
    if(getrotations)
      //bspotrotations = Teuchos::rcp(new Epetra_MultiVector(*bspotrotationsrow),true);
      for(int i=0; i<bspotrotations->NumVectors(); i++)
        for(int j=0; j<bspotrotations->MyLength(); j++)
          (*bspotrotations)[i][j] = (*bspotrotationsrow)[i][j];
  }
  return;
} // GetNodalBindingSpotPositionsFromDisVec()

/*----------------------------------------------------------------------*
 | Updates positions and rotations of binding spots        mueller 10/12|
 *----------------------------------------------------------------------*/
void STATMECH::StatMechManager::GetInterpolatedBindingSpotPositions(const Epetra_Vector&              discol,
                                                                    Teuchos::RCP<Epetra_MultiVector>  bspotpositions,
                                                                    Teuchos::RCP<Epetra_MultiVector>  bspotrotations)
{
  Teuchos::RCP<Epetra_MultiVector> bspotpositionsrow = Teuchos::rcp(new Epetra_MultiVector(*bspotrowmap_,3));
  Teuchos::RCP<Epetra_MultiVector> bspotrotationsrow = Teuchos::rcp(new Epetra_MultiVector(*bspotrowmap_,3));

  DRT::Element* filelement = NULL;
  DRT::Node* node0 = NULL;
  DRT::Node* node1 = NULL;

  // get current positions and rotations of this element's nodes
  std::vector<int> dofnode0;
  std::vector<int> dofnode1;

  LINALG::Matrix<3, 1> currpos;
  LINALG::Matrix<3, 1> currrot;
  LINALG::Matrix<1,2>  Ibp;
  LINALG::Matrix<3,1> ThetaXi;
  // Vector containing Nodal Positions. Only translational DOFs are evaluated
  std::vector<double> position(6);
  //Quaternions at relevant filament nodes
  LINALG::Matrix<4,1>  bQ;
  std::vector<LINALG::Matrix<4,1> >  nQ(2);
  DRT::ELEMENTS::Beam3ii* filele = NULL;

  // loop over row binding spots CHANGED
  int prevelegid = -1;
  for(int i=0; i<bspotrowmap_->NumMyElements(); i++)
  {
    // retrieve the element the binding spot belongs to and retrieve its nodes
    int elegid = (int)(*bspot2element_)[i];
    // only recalculate nodal positions and rotations if the element GID changed compared to the previous binding spot
    if(elegid!=prevelegid)
    {
      filelement = discret_->gElement(elegid);
      filele = dynamic_cast<DRT::ELEMENTS::Beam3ii*> (discret_->gElement(elegid));

      node0 = discret_->gNode(filelement->NodeIds()[0]);
      node1 = discret_->gNode(filelement->NodeIds()[1]);
      dofnode0 = discret_->Dof(node0);
      dofnode1 = discret_->Dof(node1);

      position[0] = node0->X()[0] + discol[discret_->DofColMap()->LID(dofnode0[0])];
      position[1] = node0->X()[1] + discol[discret_->DofColMap()->LID(dofnode0[1])];
      position[2] = node0->X()[2] + discol[discret_->DofColMap()->LID(dofnode0[2])];
      position[3] = node1->X()[0] + discol[discret_->DofColMap()->LID(dofnode1[0])];
      position[4] = node1->X()[1] + discol[discret_->DofColMap()->LID(dofnode1[1])];
      position[5] = node1->X()[2] + discol[discret_->DofColMap()->LID(dofnode1[2])];

      // NodeShift
      UnshiftPositions(position);

      // To interpolate rotations we need to get the two nodal Triads (Quaternions) of the Element
      for(int j=0;j<2;j++)
        for(int k=0;k<4;k++)
          nQ[j](k) = (filele->Quaternion()[j])(k);

      prevelegid = elegid;
    }

    // Interpolation of Positions
    DRT::UTILS::shape_function_1D(Ibp,(double)(*bspotxi_)[bspotrowmap_->GID(i)],filelement->Shape());
    currpos.PutScalar(0);
    for(int j=0;j<3;j++)
      currpos(j)= Ibp(0)*position[j]+Ibp(1)*position[j+3];

    /*if bspot currently has coordinate value greater than statmechparams_.get<double>("PeriodLength",0.0),or smaller than 0
     *it is shifted by -statmechparams_.get<double>("PeriodLength",0.0) sufficiently often to lie again in the domain*/
    for(int j=0;j<(int)(periodlength_->size());j++)
    {  if(currpos(j) > periodlength_->at(j))
          currpos(j) -= periodlength_->at(j)*floor(currpos(j)/periodlength_->at(j));

       if(currpos(j) < 0.0)
         currpos(j) -= periodlength_->at(j)*floor(currpos(j)/periodlength_->at(j));
    }

    //Interpolation of Rotations
    currrot.PutScalar(0);
    InterpolateTriadonBindingSpot(Ibp,nQ,currrot,bQ);

    for(int j=0;j<3;j++)
    {
      (*bspotpositionsrow)[j][i] = currpos(j);
      (*bspotrotationsrow)[j][i] = currrot(j);
    }
  }
  // Export row map entries to column maps
  CommunicateMultiVector(bspotpositionsrow, bspotpositions, false, true);
  CommunicateMultiVector(bspotrotationsrow, bspotrotations, false, true);

  // DEBUGGING CHECK wether bspotpositions make sense
/*
  if(discret_->Comm().MyPID()==0)
  {
  std::cout<< " **************** CHECK WETHER BSPOTPOSITIONS MAKE SENSE DISTANCES  ****************** " << std::endl;
  double ll;
  std::cout << " Anzahl BSPOTS : " << bspotpositions.MyLength() << std::endl;
  for(int ii=1;ii<bspotpositions.MyLength();ii++)
  {
      ll= sqrt(pow(bspotpositions[2][ii]-bspotpositions[2][ii-1],2)+pow(bspotpositions[1][ii]-bspotpositions[1][ii-1],2)+pow(bspotpositions[0][ii]-bspotpositions[0][ii-1],2));
      std::cout<< " BSPOTDISTANCE "<< ll<< std::endl;
  }
  for(int ii=0;ii<bspotpositions.MyLength();ii++)
    std::cout<< "BSPOTXI " << (*bspotxi_)[ii] << " BSPOT2NODES " << (*bspot2nodes_)[0][ii] << " "<<(*bspot2nodes_)[1][ii] <<   std::endl;

  std::cout<< " **************** CHECK WETHER BSPOTPOSITIONS MAKE SENSE POSITIONS  ****************** " << std::endl;
  std::cout<< bspotpositions << std::endl;
  }
  */
}

/*------------------------------------------------------------------------------------*
 | (private) update internodal triads at binding positions                            |
 | This is done by a transfer of the current rotations into Quaternions Mueller 11/12 |
 *------------------------------------------------------------------------------------*/
void STATMECH::StatMechManager::GetBindingSpotTriads(const Teuchos::RCP<Epetra_MultiVector> bspotrotations,
                                                     Teuchos::RCP<Epetra_MultiVector>       bspottriads)
{
  if (DRT::INPUT::IntegralValue<int>(statmechparams_, "CHECKORIENT") ||
      filamentmodel_ == statmech_filament_helical ||
      linkermodel_ == statmech_linker_stdintpol ||
      linkermodel_ == statmech_linker_activeintpol ||
      linkermodel_ == statmech_linker_bellseqintpol)
  {
    if (linkermodel_ == statmech_linker_stdintpol || linkermodel_ == statmech_linker_activeintpol || linkermodel_ == statmech_linker_bellseqintpol)
      GetInterpolatedBindingSpotTriads(bspotrotations, bspottriads);
    else
      GetElementBindingSpotTriads(bspottriads);
  }
  return;
}

/*----------------------------------------------------------------------*
 | (private) update nodal triads                           mueller 1/11 |
 *----------------------------------------------------------------------*/
void STATMECH::StatMechManager::GetElementBindingSpotTriads(Teuchos::RCP<Epetra_MultiVector> nodetriads)
{
  //first get triads at all row nodes
  Teuchos::RCP<Epetra_MultiVector> nodetriadsrow = Teuchos::rcp(new Epetra_MultiVector(*(discret_->NodeRowMap()), 4, true));
  //update nodaltriads_
  for (int i=0; i<discret_->NodeRowMap()->NumMyElements(); i++)
  {
    //lowest GID of any connected element (the related element cannot be a crosslinker, but has to belong to the actual filament discretization)
    int lowestid(((discret_->lRowNode(i)->Elements())[0])->Id());
    int lowestidele(0);
    for (int j=0; j<discret_->lRowNode(i)->NumElement(); j++)
      if (((discret_->lRowNode(i)->Elements())[j])->Id() < lowestid)
      {
        lowestid = ((discret_->lRowNode(i)->Elements())[j])->Id();
        lowestidele = j;
      }

    //check type of element (orientation triads are not for all elements available in the same way
    DRT::ElementType & eot = ((discret_->lRowNode(i)->Elements())[lowestidele])->ElementType();
    //if element is of type beam3ii get nodal triad
    if (eot == DRT::ELEMENTS::Beam3iiType::Instance())
    {
      DRT::ELEMENTS::Beam3ii* filele = NULL;
      filele = dynamic_cast<DRT::ELEMENTS::Beam3ii*> (discret_->lRowNode(i)->Elements()[lowestidele]);

      //check whether crosslinker is connected to first or second node of that element
      int nodenumber = 0;
      if(discret_->lRowNode(i)->Id() == ((filele->Nodes())[1])->Id() )
        nodenumber = 1;

      //save nodal triad of this node in nodaltriadrow
      for(int j=0; j<4; j++)
        (*nodetriadsrow)[j][i] = ((filele->Qnew())[nodenumber])(j);
    }
    else if (eot == DRT::ELEMENTS::Beam3Type::Instance())
    {
      DRT::ELEMENTS::Beam3* filele = NULL;
      filele = dynamic_cast<DRT::ELEMENTS::Beam3*> (discret_->lRowNode(i)->Elements()[lowestidele]);

      //approximate nodal triad by triad at the central element Gauss point (assuming 2-noded beam elements)
      for(int j=0; j<4; j++)
        (*nodetriadsrow)[j][i] = ((filele->Qnew())[0])(j);
    }
    else
      dserror("Filaments have to be discretized with beam3ii elements for orientation check!!!");
  }
  // communicate the appropriate vector
  if(nodetriads->MyLength()==discret_->NodeColMap()->NumMyElements())
    CommunicateMultiVector(nodetriadsrow, nodetriads, false, true, false);
  else
    nodetriads = Teuchos::rcp(new Epetra_MultiVector(*nodetriadsrow));
}//StatMechManager::GetNodalTriads

/*------------------------------------------------------------------------------------*
 | (private) update internodal triads at binding positions                            |
 | This is done by a transfer of the current rotations into Quaternions Mueller 10/12 |
 *------------------------------------------------------------------------------------*/
void STATMECH::StatMechManager::GetInterpolatedBindingSpotTriads(const Teuchos::RCP<Epetra_MultiVector> bspotrotations,
                                                                 Teuchos::RCP<Epetra_MultiVector>       bspottriads)
{
  Teuchos::RCP<Epetra_MultiVector>  bspottriadsrow = Teuchos::rcp(new Epetra_MultiVector(*bspotrowmap_, 4, true));
  Teuchos::RCP<Epetra_MultiVector>  bspotrotationsrow = Teuchos::rcp(new Epetra_MultiVector(*bspotrowmap_,3));
  //Get my bspotrotations in row format
  CommunicateMultiVector(bspotrotationsrow, bspotrotations,true,false,false,true);

  LINALG::Matrix<4,1> tempQ;
  LINALG::Matrix<3,1> tempM;
  for(int i=0; i<bspotrowmap_->NumMyElements(); i++)
  {
    for(int k=0;k<3;k++)
      tempM(k)=(*bspotrotationsrow)[k][i];
    LARGEROTATIONS::angletoquaternion(tempM,tempQ);
    for(int j=0;j<4;j++)
      (*bspottriadsrow)[j][i]=tempQ(j);
  }
  //export nodaltriadsrow to col map
  if(bspottriads->MyLength()==bspotcolmap_->NumMyElements())
    CommunicateMultiVector(bspottriadsrow, bspottriads, false, true);
  else
    bspottriads = Teuchos::rcp(new Epetra_MultiVector(*bspottriadsrow));
}

/*----------------------------------------------------------------------------------------------*
 | Shift vector in accordance to periodic BCs                                    mueller 10/12  |
 *----------------------------------------------------------------------------------------------*/
void STATMECH::StatMechManager::ShiftPositions(std::vector<double>& pos, const int& ndim, const int& nnode)
{
  return;
}


/*----------------------------------------------------------------------------------------------*
 | Reverse the effect of periodic BCs on a vector                                mueller 10/12  |
 *----------------------------------------------------------------------------------------------*/
void STATMECH::StatMechManager::UnshiftPositions(std::vector<double>& pos, const int& ndim, const int& nnode)
{
  if (periodlength_->at(0) > 0.0)
  {
    for(int i=1;i<nnode;i++)
    {
      for(int dof= ndim - 1; dof > -1; dof--)
      {
        if( fabs( pos[3*i+dof] + periodlength_->at(dof) - pos[dof] ) < fabs( pos[3*i+dof] - pos[dof] ) )
          pos[3*i+dof] += periodlength_->at(dof);
        if( fabs( pos[3*i+dof] - periodlength_->at(dof) - pos[dof]) < fabs( pos[3*i+dof] - pos[dof] ) )
          pos[3*i+dof] -= periodlength_->at(dof);
      }
    }
  }
  return;
}

/*------------------------------------------------------------------------------*
 | Interpolation of the binding spot triads between nodes mueller 10/12 |
 *------------------------------------------------------------------------------*/
void STATMECH::StatMechManager::InterpolateTriadonBindingSpot(const LINALG::Matrix<1, 2>&        Ibp,
                                                              std::vector<LINALG::Matrix<4,1> >& bQnew,
                                                              LINALG::Matrix<3, 1>&              ThetaXi,
                                                              LINALG::Matrix<4,1>&               bQxi)
{
  LINALG::Matrix<3,1> rphiIJ;
  std::vector<LINALG::Matrix<3,1> >  rPsili;
  rPsili.resize(2);
  LINALG::Matrix<3,1>  rPsil;
  LINALG::Matrix<4,1>  rQr;
  LINALG::Matrix<4,1> rQIJ;
  LINALG::Matrix<4,1> rQIJhalf;
  LINALG::Matrix<4,1>  rQl;
  LINALG::Matrix<4,1>  rQli;
  //rotation quaternion at Binding Position
  LINALG::Matrix<4,1>  rQxi;

  LARGEROTATIONS::quaternionproduct(bQnew[1],LARGEROTATIONS::inversequaternion(bQnew[0]),rQIJ);
  LARGEROTATIONS::quaterniontoangle(rQIJ,rphiIJ);
  rphiIJ.Scale(0.5);
  LARGEROTATIONS::angletoquaternion(rphiIJ,rQIJhalf);
  LARGEROTATIONS::quaternionproduct(rQIJhalf,bQnew[0],rQr);

  for (int node=0; node<2; ++node)
  {
    LARGEROTATIONS::quaternionproduct(bQnew[node],LARGEROTATIONS::inversequaternion(rQr),rQli);
    LARGEROTATIONS::quaterniontoangle(rQli,rPsili[node]);
  }
  rPsil.PutScalar(0);
    for(int node=0;node<2;node++)
      for(int i=0;i<3;i++)
        rPsil(i)+= Ibp(node)*rPsili[node](i);

  LARGEROTATIONS::angletoquaternion(rPsil,rQl);
  LARGEROTATIONS::quaternionproduct(rQl,rQr,rQxi);
  LARGEROTATIONS::quaterniontoangle(rQxi,ThetaXi);
  bQxi=rQxi;
}


/*----------------------------------------------------------------------*
 | Shifts current position of nodes so that they comply with periodic   |
 | boundary conditions                                       cyron 04/10|
 *----------------------------------------------------------------------*/
void STATMECH::StatMechManager::PeriodicBoundaryShift(Epetra_Vector& disrow,
                                                      int            ndim,
                                                      const double&  timen,
                                                      const double&  dt)
{
  double starttime = actiontime_->at((int)(actiontime_->size()-1));
  double shearamplitude = statmechparams_.get<double> ("SHEARAMPLITUDE", 0.0);
  int curvenumber = statmechparams_.get<int> ("CURVENUMBER", -1)-1;
  int oscilldir = statmechparams_.get<int> ("DBCDISPDIR", -1)-1;


  //only if period length >0 has been defined periodic boundary conditions are swithced on
  if (periodlength_->at(0) > 0.0)
  {
    for (int i=0; i<discret_->NumMyRowNodes(); i++)
    {
      //get a pointer at i-th row node
      const DRT::Node* node = discret_->lRowNode(i);

      //get GIDs of this node's degrees of freedom
      std::vector<int> dofnode = discret_->Dof(node);

      for (int j=ndim-1; j>-1; j--)
      {
        //current coordinate value
        double xcurr = node->X()[j] + disrow[discret_->DofRowMap()->LID(dofnode[j])];

        /*if node currently has coordinate value greater than periodlength,
         *it is shifted by -periodlength sufficiently often to lie again in the domain*/
        if (xcurr > (*periodlength_)[j])
        {
          disrow[discret_->DofRowMap()->LID(dofnode[j])] -= (*periodlength_)[j]*floor(xcurr/(*periodlength_)[j]);

          /*the upper domain surface orthogonal to the z-direction may be subject to shear Dirichlet boundary condition; the lower surface
           *may fixed by DBC. To avoid problems when nodes exit the domain through the upper z-surface and reenter through the lower
           *z-surface, the shear has to be substracted from nodal coordinates in that case */
          if (j == 2 && curvenumber >= 0 && timen > starttime && fabs(timen-starttime)>dt/1e4)
            disrow[discret_->DofRowMap()->LID(dofnode[oscilldir])] -= shearamplitude * DRT::Problem::Instance()->Curve(curvenumber).f(timen);
        }
        /*if node currently has coordinate value smaller than zero, it is shifted by periodlength sufficiently often
         *to lie again in the domain*/
        if (node->X()[j] + disrow[discret_->DofRowMap()->LID(dofnode[j])] < 0.0)
        {
          disrow[discret_->DofRowMap()->LID(dofnode[j])] -= (*periodlength_)[j]*floor(xcurr/(*periodlength_)[j]);

          /*the upper domain surface orthogonal to the z-direction may be subject to shear Dirichlet boundary condition; the lower surface
           *may be fixed by DBC. To avoid problems when nodes exit the domain through the lower z-surface and reenter through the upper
           *z-surface, the shear has to be added to nodal coordinates in that case */
          if (j == 2 && curvenumber >= 0 && timen > starttime && fabs(timen-starttime)>dt/1e4)
            disrow[discret_->DofRowMap()->LID(dofnode[oscilldir])] += shearamplitude * DRT::Problem::Instance()->Curve(curvenumber).f(timen);
        }
      }
    }
  }
  return;
}

/*------------------------------------------------------------------------*
 | This function loops through all the elements of the discretization and |
 | tests whether beam3 are broken by periodic boundary conditions in the  |
 | reference configuration; if yes initial values of curvature and jacobi |
 | determinants are adapted in a proper way                    cyron 02/10|
 *-----------------------------------------------------------------------*/
void STATMECH::StatMechManager::PeriodicBoundaryBeam3Init(DRT::Element* element)
{
  DRT::ELEMENTS::Beam3* beam = dynamic_cast<DRT::ELEMENTS::Beam3*> (element);

  //3D beam elements are embedded into R^3:
  const int ndim = 3;

  /*get reference configuration of beam3 element in proper format for later call of SetUpReferenceGeometry;
   * note that rotrefe for beam3 elements is related to the entry in the global total Lagrange displacement
   * vector related to a certain rotational degree of freedom; as the displacement is initially zero also
   * rotrefe is set to zero here*/
  std::vector<double> xrefe(beam->NumNode() * ndim, 0);
  std::vector<double> rotrefe(beam->NumNode() * ndim, 0);

  for (int i=0; i<beam->NumNode(); i++)
    for (int dof = 0; dof < ndim; dof++)
    {
      xrefe[3* i + dof] = beam->Nodes()[i]->X()[dof];
      rotrefe[3* i + dof] = 0.0;
    }

  /*loop through all nodes except for the first node which remains fixed as reference node; all other nodes are
   * shifted due to periodic boundary conditions if required*/
  for (int i=1; i<beam->NumNode(); i++)
  {
    for (int dof = 0; dof < ndim; dof++)
    {
      /*if the distance in some coordinate direction between some node and the first node becomes smaller by adding or subtracting
       * the period length, the respective node has obviously been shifted due to periodic boundary conditions and should be shifted
       * back for evaluation of element matrices and vectors; this way of detecting shifted nodes works as long as the element length
       * is smaller than half the periodic length*/
      if (fabs((beam->Nodes()[i]->X()[dof]) + periodlength_->at(dof) - (beam->Nodes()[0]->X()[dof])) < fabs((beam->Nodes()[i]->X()[dof]) - (beam->Nodes()[0]->X()[dof])))
        xrefe[3* i + dof] += periodlength_->at(dof);

      if (fabs((beam->Nodes()[i]->X()[dof]) - periodlength_->at(dof) - (beam->Nodes()[0]->X()[dof])) < fabs((beam->Nodes()[i]->X()[dof]) - (beam->Nodes()[0]->X()[dof])))
        xrefe[3* i + dof] -= periodlength_->at(dof);
    }
  }

  /*SetUpReferenceGeometry is a templated function; note that the third argument "true" is necessary as all beam elements
   * have already been initialized once upon reading input file*/
  switch (beam->NumNode())
  {
    case 2:
    {
      beam->SetUpReferenceGeometry<2> (xrefe, rotrefe, true);
      break;
    }
    case 3:
    {
      beam->SetUpReferenceGeometry<3> (xrefe, rotrefe, true);
      break;
    }
    case 4:
    {
      beam->SetUpReferenceGeometry<4> (xrefe, rotrefe, true);
      break;
    }
    case 5:
    {
      beam->SetUpReferenceGeometry<5> (xrefe, rotrefe, true);
      break;
    }
    default:
      dserror("Only Line2, Line3, Line4 and Line5 Elements implemented.");
    break;
  }
  return;
}

/*------------------------------------------------------------------------*
 | Beam3ii initialization when periodic BCs are applied        cyron 02/10|
 *-----------------------------------------------------------------------*/
void STATMECH::StatMechManager::PeriodicBoundaryBeam3iiInit(DRT::Element* element)
{
  // note: in analogy to PeriodicBoundaryBeam3Init()

  DRT::ELEMENTS::Beam3ii* beam = dynamic_cast<DRT::ELEMENTS::Beam3ii*>(element);

  const int ndim = 3;

  std::vector<double> xrefe(beam->NumNode()*ndim,0);
  std::vector<double> rotrefe(beam->NumNode()*ndim,0);

  for(int i=0;i<beam->NumNode();i++)
    for(int dof=0; dof<ndim; dof++)
    {
      xrefe[3*i+dof] = beam->Nodes()[i]->X()[dof];
      rotrefe[3*i+dof] = 0.0;
    }

  for(int i=1;i<beam->NumNode();i++)
  {
    for(int dof=0; dof<ndim; dof++)
    {
      if( fabs( (beam->Nodes()[i]->X()[dof]) + periodlength_->at(dof) - (beam->Nodes()[0]->X()[dof]) ) < fabs( (beam->Nodes()[i]->X()[dof]) - (beam->Nodes()[0]->X()[dof]) ) )
        xrefe[3*i+dof] += periodlength_->at(dof);
      if( fabs( (beam->Nodes()[i]->X()[dof]) - periodlength_->at(dof) - (beam->Nodes()[0]->X()[dof]) ) < fabs( (beam->Nodes()[i]->X()[dof]) - (beam->Nodes()[0]->X()[dof]) ) )
        xrefe[3*i+dof] -= periodlength_->at(dof);
    }
  }

  switch(beam->NumNode())
  {
    case 2:
    {
      beam->SetUpReferenceGeometry<2>(xrefe,rotrefe,true);
      break;
    }
    case 3:
    {
      beam->SetUpReferenceGeometry<3>(xrefe,rotrefe,true);
      break;
    }
    case 4:
    {
      beam->SetUpReferenceGeometry<4>(xrefe,rotrefe,true);
      break;
    }
    case 5:
    {
      beam->SetUpReferenceGeometry<5>(xrefe,rotrefe,true);
      break;
    }
    default:
      dserror("Only Line2, Line3, Line4 and Line5 Elements implemented.");
    break;
  }
  return;
}

/*------------------------------------------------------------------------*
 | BeamCL initialization when periodic BCs are applied       mueller 10/12|
 *-----------------------------------------------------------------------*/
void STATMECH::StatMechManager::PeriodicBoundaryBeamCLInit(DRT::Element* element)
{
  // note: in analogy to PeriodicBoundaryBeam3Init()

  DRT::ELEMENTS::BeamCL* beam = dynamic_cast<DRT::ELEMENTS::BeamCL*>(element);
  if(beam->NumNode()!= 4)
    dserror("PeriodicBoundaryBeamCLInit() only implemented for 4-noded BeamCL");
  const int ndim = 3;

  std::vector<double> xrefe(6,0);
  std::vector<double> rotrefe(6,0);

  xrefe=beam->XRef();
  for(int i=0;i<6;i++)
   rotrefe[i]=0.0;
  for(int dof=0; dof<ndim; dof++)
  {
    if( fabs( beam->XRef()[3+dof] + periodlength_->at(dof) - beam->XRef()[dof] ) < fabs( beam->XRef()[3+dof] - beam->XRef()[dof] ) )
      xrefe[3+dof] += periodlength_->at(dof);
    if( fabs( beam->XRef()[3+dof] - periodlength_->at(dof) - beam->XRef()[dof] ) < fabs( beam->XRef()[3+dof] - beam->XRef()[dof] ) )
      xrefe[3+dof] -= periodlength_->at(dof);
  }

  beam->SetUpReferenceGeometry<2>(xrefe,rotrefe,true);
}

/*------------------------------------------------------------------------*
 | Beam3eb initialization when periodic BCs are applied      mueller 03/13|
 *-----------------------------------------------------------------------*/
void STATMECH::StatMechManager::PeriodicBoundaryBeam3ebInit(DRT::Element* element)
{
  // note: in analogy to PeriodicBoundaryBeam3Init()

//  DRT::ELEMENTS::Beam3eb* beam = dynamic_cast<DRT::ELEMENTS::Beam3eb*>(element);
//  const int ndim = 3;
//  std::vector<double> xrefe(beam->NumNode()*ndim,0);
//
//  for(int i=0;i<beam->NumNode();i++)
//    for(int dof=0; dof<ndim; dof++)
//      xrefe[3*i+dof] = beam->Nodes()[i]->X()[dof];
//
//  for(int i=1;i<beam->NumNode();i++)
//  {
//    for(int dof=0; dof<ndim; dof++)
//    {
//      if( fabs( (beam->Nodes()[i]->X()[dof]) + periodlength_->at(dof) - (beam->Nodes()[0]->X()[dof]) ) < fabs( (beam->Nodes()[i]->X()[dof]) - (beam->Nodes()[0]->X()[dof]) ) )
//        xrefe[3*i+dof] += periodlength_->at(dof);
//      if( fabs( (beam->Nodes()[i]->X()[dof]) - periodlength_->at(dof) - (beam->Nodes()[0]->X()[dof]) ) < fabs( (beam->Nodes()[i]->X()[dof]) - (beam->Nodes()[0]->X()[dof]) ) )
//        xrefe[3*i+dof] -= periodlength_->at(dof);
//    }
//  }
//  beam->SetUpReferenceGeometry(xrefe,true);
}

/*------------------------------------------------------------------------*
 | This function loops through all the elements of the discretization and |
 | tests whether truss3 are broken by periodic boundary conditions in the |
 | reference configuration; if yes initial values of jacobi determinants  |
 | are adapted in a proper way                                 cyron 03/10|
 *-----------------------------------------------------------------------*/
void STATMECH::StatMechManager::PeriodicBoundaryTruss3Init(DRT::Element* element)
{
  DRT::ELEMENTS::Truss3* truss = dynamic_cast<DRT::ELEMENTS::Truss3*> (element);

  //3D truss elements are embeddet into R^3:
  const int ndim = 3;

  /*get reference configuration of truss3 element in proper format for later call of SetUpReferenceGeometry*/
  std::vector<double> xrefe(truss->NumNode() * ndim, 0);

  for (int i=0; i<truss->NumNode(); i++)
    for (int dof = 0; dof < ndim; dof++)
      xrefe[3* i + dof] = truss->Nodes()[i]->X()[dof];

  /*loop through all nodes except for the first node which remains fixed as reference node; all other nodes are
   * shifted due to periodic boundary conditions if required*/
  for (int i=1; i<truss->NumNode(); i++)
  {
    for (int dof = 0; dof < ndim; dof++)
    {
      /*if the distance in some coordinate direction between some node and the first node becomes smaller by adding or subtracting
       * the period length, the respective node has obviously been shifted due to periodic boundary conditions and should be shifted
       * back for evaluation of element matrices and vectors; this way of detecting shifted nodes works as long as the element length
       * is smaller than half the periodic length*/
      if (fabs((truss->Nodes()[i]->X()[dof]) + periodlength_->at(dof) - (truss->Nodes()[0]->X()[dof])) < fabs((truss->Nodes()[i]->X()[dof]) - (truss->Nodes()[0]->X()[dof])))
        xrefe[3* i + dof] += periodlength_->at(dof);

      if (fabs((truss->Nodes()[i]->X()[dof]) - periodlength_->at(dof) - (truss->Nodes()[0]->X()[dof])) < fabs((truss->Nodes()[i]->X()[dof]) - (truss->Nodes()[0]->X()[dof])))
        xrefe[3* i + dof] -= periodlength_->at(dof);
    }
  }

  /*note that the third argument "true" is necessary as all truss elements have already been initialized once upon reading input file*/
  truss->SetUpReferenceGeometry(xrefe, true);
}

/*----------------------------------------------------------------------*
 | Assign crosslink molecules and nodes to volume partitions            |
 |                                                (public) mueller 08/10|
 *----------------------------------------------------------------------*/
void STATMECH::StatMechManager::PartitioningAndSearch(const Epetra_MultiVector&         bspotpositions,
                                                      Epetra_MultiVector&               bspottriadscol,
                                                      Teuchos::RCP<Epetra_MultiVector>& neighbourslid)
{
  std::vector<double> limit(6,0.0);
  if(periodlength_->at(0)==0.0)
  {
    // initialize
    for(int i=0; i<(int)limit.size(); i++)
    {
      if(i%2==0)
        limit[i] = 1e9;
      else
        limit[i] = -1e9;
    }
    // extreme values
    // binding spots
    for(int i=0; i<(int)bspotpositions.NumVectors(); i++)
    {
      for(int j=0; j<(int)bspotpositions.MyLength(); j++)
      {
        // max
        if(limit[2*i+1]<bspotpositions[i][j])
          limit[2*i+1] = bspotpositions[i][j];
        // min
        if(limit[2*i]>bspotpositions[i][j])
          limit[2*i] = bspotpositions[i][j];
      }
    }
    // linker positions
    for(int i=0; i<(int)crosslinkerpositions_->NumVectors(); i++)
    {
      for(int j=0; j<(int)crosslinkerpositions_->MyLength(); j++)
      {
        // max
        if(limit[2*i+1]<(*crosslinkerpositions_)[i][j])
          limit[2*i+1] = (*crosslinkerpositions_)[i][j];
        // min
        if(limit[2*i]>(*crosslinkerpositions_)[i][j])
          limit[2*i] = (*crosslinkerpositions_)[i][j];
      }
    }
  }
  else
  {
    for(int i=0; i<(int)limit.size(); i++)
      if(i%2==0)
        limit[i] = 0.0;
      else
        limit[i] = (*periodlength_)[(i-1)/2];
  }

  //filament binding spots and crosslink molecules are indexed according to their positions within the boundary box volume
  std::vector<std::vector<std::vector<int> > > bspotinpartition;
  for(int i=0; i<(int)searchres_->size(); i++)
    bspotinpartition.push_back(std::vector<std::vector<int> >((*searchres_)[i], std::vector<int>()));

  /*nodes*/
  // loop over node positions to map their column map LIDs to partitions
  for(int i=0; i<(int)bspotinpartition.size(); i++) // note: bspotinpartition.size==3
  {
    for(int j=0;j<bspotpositions.MyLength();j++)
    {
      int partition = (int)std::floor((bspotpositions[i][j]-limit[2*i])/(limit[2*i+1]-limit[2*i])*(double)(*searchres_)[i]);
      if(partition==(int)(*searchres_)[i])
        partition--;
      bspotinpartition[i][partition].push_back(bspotcolmap_->LID(j)); //column lid
    }
  }
  /*crosslink molecules*/
  // Export crosslinkerpositions_ to transfermap_ format (row map format for crosslink molecules)
  Teuchos::RCP<Epetra_MultiVector> crosslinkerpositionstrans = Teuchos::rcp(new Epetra_MultiVector(*transfermap_, 3, true));
  Teuchos::RCP<Epetra_Vector> numbondtrans = Teuchos::rcp(new Epetra_Vector(*transfermap_, true));
  Teuchos::RCP<Epetra_MultiVector> crosslinkpartitiontrans = Teuchos::rcp(new Epetra_MultiVector(*transfermap_, 3, false));

  CommunicateVector(numbondtrans, numbond_, true, false, false, true);
  CommunicateMultiVector(crosslinkerpositionstrans, crosslinkerpositions_, true, false, false, true);

  for(int i=0; i<crosslinkpartitiontrans->NumVectors(); i++)
  {
    for(int j=0; j<crosslinkpartitiontrans->MyLength(); j++)
    {
      // mark entries with double-bonded crosslink molecules
      if((*numbondtrans)[j]>1.9)
      {
        (*crosslinkpartitiontrans)[i][j] = -1.0;
        continue;
      }
      else
      {
        int partition = (int)std::floor(((*crosslinkerpositionstrans)[i][j]-limit[2*i])/(limit[2*i+1]-limit[2*i])*(double)(*searchres_)[i]);
        if(partition==(*searchres_)[i])
          partition--;
        (*crosslinkpartitiontrans)[i][j] = partition;
      }
    }
  }

  // detection of nodes within search proximity of the crosslink molecules
  DetectNeighbourNodes(bspotpositions, bspotinpartition, *numbondtrans, *crosslinkerpositionstrans, *crosslinkpartitiontrans, bspottriadscol, neighbourslid);
  return;
}//void StatMechManager::PartitioningAndSearch

/*----------------------------------------------------------------------*
 | detect neighbour nodes to crosslink molecules (public) mueller (7/10)|
 *----------------------------------------------------------------------*/
void STATMECH::StatMechManager::DetectNeighbourNodes(const Epetra_MultiVector&                           bspotpositions,
                                                     const std::vector<std::vector<std::vector<int> > >& bspotinpartition,
                                                     const Epetra_Vector&                                numbond,
                                                     const Epetra_MultiVector&                           crosslinkerpositions,
                                                     const Epetra_MultiVector&                           crosslinkpartitions,
                                                     const Epetra_MultiVector&                           bspottriadscol,
                                                     Teuchos::RCP<Epetra_MultiVector>&                   neighbourslid)
  {
  /* Description:
   * A vector containing the volume partitions of the crosslink molecules is handed over to this method. We loop over the partition
   * number of the first (x) component and its two neighbouring partition numbers, hence gathering information on all nodes
   * within these three layers.
   *
   * Now, we check, whether or not an LID of the first component matches one of the LIDs in the second component's partition layers
   * (once again, we check the given partition number and its immediate neighbours).
   * If so, we head to the third component and repeat this procedure.
   * If a match is found for all components, we can be sure that the crosslink molecule in question lies within
   * the 27 partitions encompassing the crosslink molecule partition.
   *
   * Eventually, the distance between the crosslink molecule and the filament node is calculated.
   * If the distance lies beneath the boundaries given by the crosslinker lengths rmin and rmax,
   * the node's LID is added to a storage vector for further use.
   *
   * After having found a match in the next component, we exit the loop to avoid unnecessary computational cost
   */
  // distribute information of crosslinkerbond_ to processorspecific maps
  Teuchos::RCP<Epetra_MultiVector> crosslinkerbondtrans = Teuchos::rcp(new Epetra_MultiVector(*transfermap_, 2, true));
  CommunicateMultiVector(crosslinkerbondtrans, crosslinkerbond_,true,false,false,true);

  std::vector<std::vector<int> > neighbournodes(crosslinkpartitions.MyLength(), std::vector<int>());

  int maxneighbourslocal = 0;
  int maxneighboursglobal = 0;

  for(int part=0; part<crosslinkpartitions.MyLength(); part++)
  {
    // i.e. numbond!=2.0
    if(crosslinkpartitions[0][part]>-0.9)
    {
      double rmin = (statmechparams_.get<double>("R_LINK", 0.0)-statmechparams_.get<double>("DeltaR_LINK", 0.0)) / 2.0;
      double rmax = (statmechparams_.get<double>("R_LINK", 0.0)+statmechparams_.get<double>("DeltaR_LINK", 0.0)) / 2.0;
      // determine search radius in accordance to bonding status (only without helical binding structure)
      // reason for this: we set the crosslinker position of molecules with numbond>0 on a node, i.e. on a binding spot of a filement/linker
      // So, the maximal radius is now R_LINK-DeltaR_LINK assuming the linker can bond in any direction.
      // When numbond==0 -> crosslinker position is thought to be in the center of a sphere with the adjusted search radii rmin and rmax.
      // In case of a helical binding spot structure, we exactly know the position of the singly bound linker. Therefore, no multiplication with 2.
      if ((int)numbond[part]>0.9)
      {
        rmin *= 2.0;
        rmax *= 2.0;
      }

      // first component
      for(int ilayer=(int)crosslinkpartitions[0][part]-1; ilayer<(int)crosslinkpartitions[0][part]+2; ilayer++)
      {
        if(ilayer>-1 && ilayer<(*searchres_)[0])
        {
          for(int i=0; i<(int)bspotinpartition[0][ilayer].size(); i++)
          {
            int tmplid = (int)bspotinpartition[0][ilayer][i];
            // second component
            for(int jlayer=(int)crosslinkpartitions[1][part]-1; jlayer<(int)crosslinkpartitions[1][part]+2; jlayer++)
            {
              if(jlayer>-1 && jlayer<(*searchres_)[1])
              {
                for(int j=0; j<(int)bspotinpartition[1][jlayer].size(); j++)
                {
                  if(bspotinpartition[1][jlayer][j]==tmplid)
                  {
                    //third component
                    for(int klayer=(int)crosslinkpartitions[2][part]-1; klayer<(int)crosslinkpartitions[2][part]+2; klayer++)
                    {
                      if(klayer>-1 && klayer<(*searchres_)[2])
                      {
                        for(int k=0; k<(int)bspotinpartition[2][klayer].size(); k++)
                        {
                          if(bspotinpartition[2][klayer][k]==tmplid)
                          {
                            // calculate distance crosslinker-node
                            LINALG::Matrix<3, 1> difference;
                            for (int l=0; l<(int)difference.M(); l++)
                              difference(l) = crosslinkerpositions[l][part]-bspotpositions[l][bspotcolmap_->GID(tmplid)];
                            // 1. criterion: distance between linker and binding spot within given interval
                            if(difference.Norm2()<rmax && difference.Norm2()>rmin)
                            {
                              // further calculations in case of helical binding spot geometry and singly bound crosslinkers
                              if(filamentmodel_ == statmech_filament_helical)
                              {
                                bool bspot1=false;
                                bool bspot2=false;
                                // 2. criterion: linker has to lie in the oriented cone with peak "binding spot location"
                                // first and second triad vectors (tangent and normal)
                                LINALG::Matrix<3,1> firstdir;
                                LINALG::Matrix<3,1> bspotvec;
                                // retrieve tangential and normal vector from binding spot quaternions
                                  LINALG::Matrix<3,3> bspottriad;
                                  // auxiliary variable for storing a triad in quaternion form
                                  LINALG::Matrix<4, 1> qnode;
                                  // triad of node on first filament which is affected by the new crosslinker
                                  for (int l=0; l<4; l++)
                                    qnode(l) = bspottriadscol[l][bspotcolmap_->GID(tmplid)];
                                  LARGEROTATIONS::quaterniontotriad(qnode, bspottriad);
                                  for (int l=0; l<(int)bspottriad.M(); l++)
                                  {
                                    firstdir(l) = bspottriad(l,0);
                                    bspotvec(l) = bspottriad(l,1);
                                  }

                                // rotation matrix around tangential vector by given angle
                                RotationAroundFixedAxis(firstdir,bspotvec,(*bspotorientations_)[bspotcolmap_->GID(tmplid)]);

                                // linker position
                                LINALG::Matrix<3,1> crossbspotdiff;
                                for(int l=0; l<(int)crossbspotdiff.M(); l++)
                                  crossbspotdiff(l) = crosslinkerpositions[l][part]-bspotpositions[l][bspotcolmap_->GID(tmplid)];
                                // line parameter of the intersection point of the line through the binding spot with the orientation of the binding spot.
                                // a)lambda must be larger than zero in order to lie within the cone
                                double lambda = (crossbspotdiff.Dot(bspotvec))/(bspotvec.Dot(bspotvec));
                                if(lambda>0)
                                {
                                  difference.Scale(1.0/difference.Norm2()); // normalized vector
                                  // b) the angle between the binding spot orientation and the vector between bspot and crosslinker must be smaller than PHIBSPOT
                                  //    Only then does the crosslinker lie within the cone
                                  if(acos(fabs(bspotvec.Dot(difference))) < statmechparams_.get<double>("PHIBSPOT", 0.524))
                                    bspot1=true;
                                  if((int)numbond[part]<0.1)
                                    bspot2=true;
                                  else //check wether first Bspot lies fullfills orientation criterium
                                  {       //id of bindingspot that was bound ealier
                                    int bspotID=0;
                                    if((int)(*crosslinkerbondtrans)[0][part] < -0.9)
                                      bspotID=(int)(*crosslinkerbondtrans)[1][part];
                                    else if((int)(*crosslinkerbondtrans)[1][part] < -0.9)
                                      bspotID=(int)(*crosslinkerbondtrans)[0][part];
                                    else //
                                      dserror("Error in crosslinker management and/or search!");
                                      //turn difference Vectors direction
                                    for(int l=0; l<(int)crossbspotdiff.M(); l++)
                                     crossbspotdiff(l) = -crossbspotdiff(l);

                                    for (int l=0; l<4; l++)
                                       qnode(l) = bspottriadscol[l][bspotID];
                                    LARGEROTATIONS::quaterniontotriad(qnode, bspottriad);
                                    for (int l=0; l<(int)bspottriad.M(); l++)
                                    {
                                      firstdir(l) = bspottriad(l,0);
                                      bspotvec(l) = bspottriad(l,1);
                                    }

                                    RotationAroundFixedAxis(firstdir,bspotvec,(*bspotorientations_)[bspotID]);

                                    lambda = (crossbspotdiff.Dot(bspotvec))/(bspotvec.Dot(bspotvec));
                                    if(lambda>0)
                                    {
                                     difference.Scale(1.0/difference.Norm2()); // normalized vector
                                    // b) the angle between the binding spot orientation and the vector between bspot and crosslinker must be smaller than PHIBSPOT
                                    //    Only then does the crosslinker lie within the cone
                                    if(acos(fabs(bspotvec.Dot(difference))) < statmechparams_.get<double>("PHIBSPOT", 0.524))
                                     bspot2=true;
                                    }
                                  }

                                  if(bspot1 && bspot2)
                                    neighbournodes[part].push_back(tmplid);
                                }
                              }
                              else // only difference criterion applied
                                neighbournodes[part].push_back(tmplid);
                            }
                            // exit loop immediately
                            break;
                          }
                        }
                      }
                    }
                    break;
                  }
                }
              }
            }
          }
        }
      }
      // "-1" indicates the possibility of a crosslink molecule becoming passive, i.e. hypothetically bonding to the same filament
      if((int)numbond[part]==1)
        neighbournodes[part].push_back(-1);
      // store local maximal number of LIDs per molecule in order to determine neighbourslid->NumVectors()
      maxneighbourslocal = std::max(maxneighbourslocal, (int)neighbournodes[part].size());
    }
  }

  // get global maximal number of LIDs per molecule
  discret_->Comm().MaxAll(&maxneighbourslocal, &maxneighboursglobal, 1);
  if(maxneighboursglobal==0)
    maxneighboursglobal = 1;
  // copy information to Epetra_MultiVector for communication
  Teuchos::RCP<Epetra_MultiVector> neighbourslidtrans = Teuchos::rcp(new Epetra_MultiVector(*transfermap_, maxneighboursglobal, false));

  /* assign "-2" (also to distinguish between 'empty' and passive crosslink molecule "-1"
   * to be able to determine entries which remain "empty" due to number of LIDs < maxneighboursglobal*/
  neighbourslidtrans->PutScalar(-2.0);
  for(int i=0; i<neighbourslidtrans->MyLength(); i++)
    for(int j=0; j<(int)neighbournodes[i].size(); j++)
      (*neighbourslidtrans)[j][i] = (double)neighbournodes[i][j];

  // make information redundant on all Procs
  neighbourslid = Teuchos::rcp(new Epetra_MultiVector(*crosslinkermap_,maxneighboursglobal,false));
  CommunicateMultiVector(neighbourslidtrans, neighbourslid, false, true);

  return;
}// StatMechManager::DetectNeighbourNodes

/*----------------------------------------------------------------------*
 | rotation around a fixed axis by angle phirot          mueller (11/11)|
 *----------------------------------------------------------------------*/
void STATMECH::StatMechManager::RotationAroundFixedAxis(LINALG::Matrix<3,1>& axis, LINALG::Matrix<3,1>& vector, const double& phirot)
{
  // unit length
  axis.Scale(1.0/axis.Norm2());
  vector.Scale(1.0/vector.Norm2());

  LINALG::Matrix<3, 3> RotMat;

  for (int j=0; j<3; j++)
    RotMat(j, j) = cos(phirot) + axis(j) * axis(j) * (1 - cos(phirot));
  RotMat(0, 1) = axis(0) * axis(1) * (1 - cos(phirot)) - axis(2) * sin(phirot);
  RotMat(0, 2) = axis(0) * axis(2) * (1 - cos(phirot)) + axis(1) * sin(phirot);
  RotMat(1, 0) = axis(1) * axis(0) * (1 - cos(phirot)) + axis(2) * sin(phirot);
  RotMat(1, 2) = axis(1) * axis(2) * (1 - cos(phirot)) - axis(0) * sin(phirot);
  RotMat(2, 0) = axis(2) * axis(0) * (1 - cos(phirot)) - axis(1) * sin(phirot);
  RotMat(2, 1) = axis(2) * axis(1) * (1 - cos(phirot)) + axis(0) * sin(phirot);

  LINALG::Matrix<3,1> tmpvec = vector;
  vector.Multiply(RotMat, tmpvec);

  return;
}


/*----------------------------------------------------------------------*
 | Searches for crosslink molecule-filament node pairs and adds actual  |
 | crosslinker elements once linking conditions are met.                |
 | (private)                                              mueller (7/10)|
 *----------------------------------------------------------------------*/
void STATMECH::StatMechManager::SearchAndSetCrosslinkers(const int&                                       istep,
                                                         const double&                                    timen,
                                                         const double&                                    dt,
                                                         const Teuchos::RCP<Epetra_MultiVector>           bspotpositions,
                                                         const Teuchos::RCP<Epetra_MultiVector>           bspotrotations,
                                                         const Epetra_Vector&                             discol,
                                                         Teuchos::RCP<CONTACT::Beam3cmanager>             beamcmanager,
                                                         bool                                             printscreen)
{
  double t0 = Teuchos::Time::wallTime();
  /*preliminaries*/
  // BINDING SPOT TRIAD UPDATE
  Teuchos::RCP<Epetra_MultiVector> bspottriadscol = Teuchos::rcp(new Epetra_MultiVector(*bspotcolmap_,4,true));
  GetBindingSpotTriads(bspotrotations, bspottriadscol);

  //get current on-rate for crosslinkers
  double kon = 0;
  double starttime = actiontime_->at(1);

  if(timen <= starttime || (timen>starttime && fabs(timen-starttime) < dt/1e4))
    kon = statmechparams_.get<double>("K_ON_start",0.0);
  else
    kon = statmechparams_.get<double>("K_ON_end",0.0);

  //probability with which a crosslinker is established between crosslink molecule and neighbour node
  double plink = 1.0 - exp( -dt*kon );
  //probability with which a crosslinker blocks two binding spots on the same filament (self-binding)
  double pself = 1.0 - exp( -dt*statmechparams_.get<double>("K_ON_SELF", 0.0) );

  //Volume partitioning, assignment of nodes and molecules to partitions, search for neighbours
  // map filament (column map, i.e. entire node set on each proc) node positions to volume partitions every SEARCHINTERVAL timesteps

  Teuchos::RCP<Epetra_MultiVector> neighbourslid;
  if(statmechparams_.get<int>("SEARCHRES",1)>0)
    PartitioningAndSearch(*bspotpositions,*bspottriadscol, neighbourslid);

#ifdef MEASURETIME
  double t1 = Teuchos::Time::wallTime();
  double t2 = -1.0e9;
#endif
  /*the following part of the code is executed on processor 0 only and leads to the decision, at which nodes crosslinks shall be set
   *processor 0 goes through all the crosslink molecules and checks whether a crosslink is to be set; this works precisely as follows:
   *(1) the crosslink molecules are looped through in a random order
   *(2) if a node has not yet reached its maximal number of crosslinks, a crosslink may be set
   *(3) a crosslink is set if and only if the node passes a probability check
   *note: a fully overlapping node map on processor 0 is imperative for the following part of the code to work correctly!*/

  // a vector indicating the crosslink molecule which is going to constitute a crosslinker element
  Teuchos::RCP<Epetra_Vector> addcrosselement = Teuchos::rcp(new Epetra_Vector(*crosslinkermap_, true));

  int numsetelements = 0;
  if(discret_->Comm().MyPID()==0)
  {
    // obtain a random order in which the crosslinkers are addressed
    std::vector<int> order = Permutation(crosslinkermap_->NumMyElements());

    for(int i=0; i<numbond_->MyLength(); i++)
    {
      int irandom = order[i];

      // only consider crosslink molecule with less than two nodes or if it is actively searching for bonding partners
      if((*numbond_)[irandom]<1.9)
      {
        // PROBABILITY CHECK & PREPARATION OF VECTORS
        // obtain a random order of neighboursLID indices
        std::vector<int> neighbourorder = Permutation(neighbourslid->NumVectors());

        // loop over neighbouring binding spots
        for(int j=0; j<neighbourslid->NumVectors(); j++)
        {
          // random index
          int index = neighbourorder[j];
          // skip this, if neighbourslid entry is '-2', meaning empty or binding spot is unavailable due to input file specs
          if((*neighbourslid)[index][irandom] > -1.9)
          {
            // check the space between already set linkers and the about to be set linker
            bool setlink = true;
            int bspotinterval = statmechparams_.get<int>("BSPOTINTERVAL", 1);
            if(bspotinterval>1 && (*neighbourslid)[index][irandom]>-0.9)
            {
              int bspotlid = (*neighbourslid)[index][irandom];
              //getting a vector consisting of pointers to all filament number conditions set
              std::vector<DRT::Condition*> filaments(0);
              discret_->GetCondition("FilamentNumber",filaments);
              int filnum = (*filamentnumber_)[bspotlid];
              int startlid = bspotlid - int(ceil(double(bspotinterval/2)));
              int endlid = bspotlid + int(ceil(double(bspotinterval/2)));
              if(bspotcolmap_->GID(startlid)<(*filaments[filnum]->Nodes())[0])
                startlid = bspotcolmap_->LID((*filaments[filnum]->Nodes())[0]);
              if(bspotcolmap_->GID(endlid)>(*filaments[filnum]->Nodes())[(int)filaments[filnum]->Nodes()->size()-1])
                endlid = bspotcolmap_->LID((*filaments[filnum]->Nodes())[(int)filaments[filnum]->Nodes()->size()-1]);

              for(int k=startlid; k<endlid+1; k++)
              {
                if((*bspotstatus_)[k]>-1)
                {
                  setlink = false;
                  break;
                }
              }
            }

            if(setlink)
            {
              // current neighbour LID
              int bspotLID = (int)(*neighbourslid)[index][irandom];

              //skip the rest of this iteration of the i-loop in case the binding spot is occupied
              if(bspotLID>-1)
                if((*bspotstatus_)[bspotLID]>-0.1)
                  continue;

              // flag indicating loop break after first new bond has been established between i-th crosslink molecule and j-th neighbour node
              bool bondestablished = false;
              // necessary condition to be fulfilled in order to set a crosslinker
              double probability = 0.0;
              // switch between probability for regular inter-filament binding and self-binding crosslinker
              if(bspotLID>=0)
                probability = plink;
              else
                probability = pself;

              if( (*uniformgen_)() < probability )
              {
                int free = 0;
                int occupied = 0;
                bool set = false;

                // check both entries of crosslinkerbond_
                for(int k=0; k<crosslinkerbond_->NumVectors(); k++)
                  if((*crosslinkerbond_)[k][irandom]<-0.9 && !set)
                  {
                    // store free bond position
                    free = k;
                    set = true;
                  }
                  else
                    occupied = k;

                switch((int)(*numbond_)[irandom])
                {
                  // free crosslink molecule creating a bond
                  case 0:
                  {
                    // crosslinkers only allowed on the horizontal filament when looking at a loom setup
                    if(networktype_ == statmech_network_casimir && (*filamentnumber_)[bspotLID]!=0)
                      break;
                    else // standard case
                    {
                      // update of crosslink molecule positions
                      LINALG::SerialDenseMatrix LID(1,1,true);
                      LID(0,0) = bspotLID;
                      (*crosslinkerbond_)[free][irandom] = bspotcolmap_->GID(bspotLID);
                      // update status of binding spot by putting in the crosslinker id
                      ((*bspotstatus_)[bspotLID]) = irandom;
                      // increment the number of bonds of this crosslinker
                      ((*numbond_)[irandom]) = 1.0;
                      CrosslinkerIntermediateUpdate(*bspotpositions, LID, irandom);
                      bondestablished = true;
                    }
                  }
                  break;
                  // one bond already exists -> establish a second bond/passive crosslink molecule
                  // Potentially, add an element (later)
                  case 1:
                  {
                    //Col map LIDs of nodes to be crosslinked
                     Epetra_SerialDenseMatrix LID(2,1);
                     LID(0,0) = bspotcolmap_->LID((int)(*crosslinkerbond_)[occupied][irandom]);
                    // distinguish between a real bspotLID and the entry '-1', which indicates a passive crosslink molecule
                    if(bspotLID>=0)
                      LID(1,0) = bspotLID;
                    else // choose the neighbor node on the same filament as bspotLID as second entry and take basisnodes_ into account
                    {
                      int currfilament = (int)(*filamentnumber_)[(int)LID(0,0)];
                      if((int)LID(0,0)<basisnodes_-1)
                      {
                        if((int)(*filamentnumber_)[(int)LID(0,0)+1]==currfilament)
                          LID(1,0) = LID(0,0) + 1.0;
                        else
                          LID(1,0) = LID(0,0) - 1.0;
                      }
                      if((int)LID(0,0)==basisnodes_-1)
                        if((int)(*filamentnumber_)[(int)LID(0,0)-1]==currfilament)
                          LID(1,0) = LID(0,0) - 1.0;
                    }

                    // do not do anything if a crosslinker is about to occupy two binding spots on the same filament and K_ON_SELF is zero
                    if(linkermodel_ == statmech_linker_stdintpol || linkermodel_ == statmech_linker_activeintpol || linkermodel_ == statmech_linker_bellseqintpol)
                    {
                      if((*filamentnumber_)[discret_->NodeColMap()->LID((int)(*bspot2nodes_)[0][bspotcolmap_->GID((int)LID(0,0))])] == (*filamentnumber_)[discret_->NodeColMap()->LID((int)(*bspot2nodes_)[0][bspotcolmap_->GID((int)LID(1,0))])] && statmechparams_.get<double>("K_ON_SELF",0.0)==0.0)
                        break;
                    }
                    else
                    {
                      if((*filamentnumber_)[(int)LID(0,0)]==(*filamentnumber_)[(int)LID(1,0)] && statmechparams_.get<double>("K_ON_SELF",0.0)==0.0)
                        break;
                    }

                    // do not do anything in case of a Loom set up if conditions hereafter are not met
                    //CHECK make this more elegant later
                    if(networktype_ == statmech_network_casimir)
                    {
                      if(!SetCrosslinkerLoom(LID, bspotpositions, bspottriadscol))
                        break;
                    }

                    //unit direction vector between currently considered two nodes
                    LINALG::Matrix<3,1> direction;
                    for(int j=0;j<3;j++)
                      direction(j) = (*bspotpositions)[j][(int)LID(0,0)]-(*bspotpositions)[j][(int)LID(1,0)];

                    // skip binding process, if cycle time is still detached and
                    // skip binding process, if polarity of active linker is not fulfilled
                    bool polaritycriterion = true;
                    if(linkermodel_ == statmech_linker_active || linkermodel_ == statmech_linker_activeintpol)
                      polaritycriterion = LinkerPolarityCheckAttach(bspottriadscol, LID, direction);

                    if(polaritycriterion)
                    {
                      /* In case of beam contact, evaluate whether a crosslinker would intersect with any exisiting elements
                       * if it were to be set between the two nodes considered*/
                      bool intersection = false;
                      if(beamcmanager!=Teuchos::null)
                      {
                        Epetra_SerialDenseMatrix nodecoords(3,2);
                        for(int k=0; k<nodecoords.M(); k++)
                        {
                          nodecoords(k,0) = (*bspotpositions)[k][bspotcolmap_->GID((int)LID(0,0))];
                          nodecoords(k,1) = (*bspotpositions)[k][bspotcolmap_->GID((int)LID(1,0))];
                        }
                        for(int k=0; k<(int)LID.M(); k++)
                        {
                          intersection = beamcmanager->OcTree()->IntersectBBoxesWith(nodecoords, LID);
                          if(intersection)
                            break;
                        }
                      }

                      /*Orientation Check: only if the two nodes in question pass a check for their
                       * orientation, a marker, indicating an element to be added, will be set
                       * In addition to that, if beam contact is enabled, a linker is only set if
                       * it does not intersect with other existing elements*/
                      if(CheckOrientation(direction,*bspottriadscol,LID) && !intersection)
                      {
                        numsetelements++;

                        ((*numbond_)[irandom]) = 2.0;
                        // actually existing two bonds
                        if(bspotLID>=0)
                        {
                           (*crosslinkerbond_)[free][irandom] = bspotcolmap_->GID(bspotLID);
                           ((*bspotstatus_)[bspotLID]) = irandom;
                          // set flag at irandom-th crosslink molecule that an element is to be added
                          (*addcrosselement)[irandom] = 1.0;
                          // update molecule positions
                          CrosslinkerIntermediateUpdate(*bspotpositions, LID, irandom);

                          // consider crosslinkers covering two binding spots of one filament
                          if(linkermodel_ == statmech_linker_stdintpol || linkermodel_ == statmech_linker_activeintpol || linkermodel_ == statmech_linker_bellseqintpol)
                          {
                            if((*filamentnumber_)[discret_->NodeColMap()->LID((int)(*bspot2nodes_)[0][bspotcolmap_->GID((int)LID(0,0))])] == (*filamentnumber_)[discret_->NodeColMap()->LID((int)(*bspot2nodes_)[0][bspotcolmap_->GID((int)LID(1,0))])] && statmechparams_.get<double>("K_ON_SELF",0.0)>0.0)
                              (*crosslinkonsamefilament_)[irandom] = 1.0;
                          }
                          else
                          {
                            if((*filamentnumber_)[(int)LID(0,0)]==(*filamentnumber_)[(int)LID(1,0)] && statmechparams_.get<double>("K_ON_SELF",0.0)>0.0)
                              (*crosslinkonsamefilament_)[irandom] = 1.0;
                          }
                        }
                        else // passive crosslink molecule
                        {
                          (*searchforneighbours_)[irandom] = 0.0;
                          LINALG::SerialDenseMatrix oneLID(1,1,true);
                          oneLID(0,0) = LID(0,0);
                          CrosslinkerIntermediateUpdate(*bspotpositions, oneLID, irandom);
                        }

                        if(DRT::INPUT::IntegralValue<INPAR::STATMECH::StatOutput>(statmechparams_, "SPECIAL_OUTPUT")==INPAR::STATMECH::statout_structanaly)
                        {
                          // set t_on and reset t_off-t_on
                          (*crosslinkunbindingtimes_)[0][irandom] = timen;
                          (*crosslinkunbindingtimes_)[1][irandom] = -1.0;
                        }

                        bondestablished = true;
                      }
                    }
                  }
                  break;
                }// switch((int)(*numbond_)[irandom])
                // for now, break j-loop after a new bond was established, i.e crosslinker elements cannot be established starting from zero bonds
                if(bondestablished)
                  break;
              }//if( (*uniformgen_)() < probability )
            }//
          }// if(neighbourslid...)
        }// for(int j=0; j<(int)neighboursLID.size(); j++)
      }//if((*numbond_)[irandom]<1.9)
    }// for(int i=0; i<numbond_->MyLength(); i++)
#ifdef MEASURETIME
    t2 = Teuchos::Time::wallTime();
#else
    if(printscreen)
      std::cout << "\nsearch + linker admin time: " << Teuchos::Time::wallTime() - t0<< " seconds";
#endif
  }// if(discret_->Comm().MypPID==0)

  /* note: searchforneighbours_ and crosslinkonsamefilament_ are not being communicated
   * to the other Procs because their information is of concern to Proc 0 only.*/
  //synchronize information about number of bonded filament nodes by exporting it to row map format and then reimporting it to column map format
  Teuchos::RCP<Epetra_Vector> bspotstatusrow = Teuchos::rcp(new Epetra_Vector(*bspotrowmap_,true));
  CommunicateVector(bspotstatusrow, bspotstatus_);

  // transfer vectors
  Teuchos::RCP<Epetra_MultiVector> crosslinkerpositionstrans = Teuchos::rcp(new Epetra_MultiVector(*transfermap_,3,true));
  Teuchos::RCP<Epetra_MultiVector> crosslinkerbondtrans = Teuchos::rcp(new Epetra_MultiVector(*transfermap_,2,true));
  Teuchos::RCP<Epetra_Vector> numbondtrans = Teuchos::rcp(new Epetra_Vector(*transfermap_, true));
  Teuchos::RCP<Epetra_Vector> addcrosselementtrans = Teuchos::rcp(new Epetra_Vector(*transfermap_, true));
  // exports and reimports
  CommunicateMultiVector(crosslinkerpositionstrans, crosslinkerpositions_);
  CommunicateMultiVector(crosslinkerbondtrans, crosslinkerbond_);
  CommunicateVector(numbondtrans, numbond_);
  CommunicateVector(addcrosselementtrans, addcrosselement);

#ifdef MEASURETIME
  double t3 = Teuchos::Time::wallTime();
#endif

  // ADDING ELEMENTS
  // the node row map is needed in order to adequately assing ownership of elements to processors
  const Epetra_Map noderowmap(*(discret_->NodeRowMap()));

  int numsetelementsall = -1;
  discret_->Comm().MaxAll(&numsetelements, &numsetelementsall, 1);

  if(numsetelementsall>0)
  {
    // rotations from displacement vector
    Teuchos::RCP<Epetra_MultiVector> nodalrotations = Teuchos::null;
    // NODAL quaternions from filament elements
    Teuchos::RCP<Epetra_MultiVector> nodalquaternions = Teuchos::null;
    if(linkermodel_ == statmech_linker_stdintpol || linkermodel_ == statmech_linker_activeintpol || linkermodel_ == statmech_linker_bellseqintpol)
    {
      nodalrotations = Teuchos::rcp(new Epetra_MultiVector(*(discret_->NodeColMap()),3));
      GetNodalBindingSpotPositionsFromDisVec(discol, Teuchos::null, nodalrotations);
      nodalquaternions = Teuchos::rcp(new Epetra_MultiVector(*(discret_->NodeColMap()),4));
      GetElementBindingSpotTriads(nodalquaternions);
    }

    // add elements to problem discretization (processor specific)
    for(int i=0; i<addcrosselement->MyLength(); i++)
    {
      if((*addcrosselement)[i]>0.9)
      {
        /*there is the problem of how to assign to a crosslinker element a GID which is certainly not used by any
         * other element; we know that for a network at the beginning (without crosslinkers) each node has a
         * connectivity of 1 or 2 so that the number of all elements used to discretize the filaments is smaller
         * than basisnodes_. Thus we choose crosslinker GIDs >= basisnodes_ for the additionally added crosslinker
         * elements. To make sure that two crosslinkers never have the same GID we add to basisnodes_ the value
         * basisnodes_*GID1, where GID1 is the GID of the first node of the crosslinker element. Then we add GID2
         * (note: GID2 < basisnodes_), which is the GID of the second node of the crosslinker element.
         * Hence basisnodes_ + GID1*basisnodes_ + GID2 always gives a GID which cannot be used by any other element;
         * the first node of the crosslinker element is always assumed to be the one with the greater GID; the owner
         * of the first node is assumed to be the owner of the crosslinker element*/

        // copied binding spot map (in order to handle the addition of elements based upon an invariant map just in case)
        Epetra_Map bspotcolmap(*bspotcolmap_);

        // obtain binding spot GID
        std::vector<int> bspotgid(2,0);
        // determine smaller and larger of the GIDs
        bspotgid.at(1) = std::min((int)(*crosslinkerbond_)[0][i],(int)(*crosslinkerbond_)[1][i]);
        bspotgid.at(0) = std::max((int)(*crosslinkerbond_)[0][i],(int)(*crosslinkerbond_)[1][i]);

        // different sizes due to different linker elements
        Teuchos::RCP<std::vector<int> > globalnodeids;
        if(linkermodel_ == statmech_linker_stdintpol || linkermodel_ == statmech_linker_activeintpol || linkermodel_ == statmech_linker_bellseqintpol)
        {
          globalnodeids = Teuchos::rcp(new std::vector<int>(4,0));
          for(int l=0; l<(int)bspot2nodes_->NumVectors(); l++)
          {
            (*globalnodeids)[l]=(int)(*bspot2nodes_)[l][bspotgid[0]];
            (*globalnodeids)[l+2]=(int)(*bspot2nodes_)[l][bspotgid[1]];
          }
        }
        else
        {
          globalnodeids = Teuchos::rcp(new std::vector<int>(2,0));
          (*globalnodeids)[0] = bspotgid[0];
          (*globalnodeids)[1] = bspotgid[1];
        }

        // calculate element GID
        int newcrosslinkerGID = (bspotgid[0] + 1)*basisnodes_ + bspotgid[1];

        /* Correction of the crosslinker GID if there happens to be another crosslinker with the same ID. This might occur
         * if two different crosslink molecules bind to the same two nodes. Then, the calculated crosslinker GID will be the
         * same in both cases, leading to problems within the discretization.
         * As long as an unused GID cannot be found, the crosslinker GID keeps getting incremented by 1.*/
        discret_->Comm().Barrier();
        while(1)
        {
          int gidexists = 1;
          // query existance of node on this Proc
          int gidonproc = (int)(discret_->HaveGlobalElement(newcrosslinkerGID));
          // sum over all processors
          discret_->Comm().MaxAll(&gidonproc, &gidexists, 1);
          // calculate new GID if necessary by shifting the initial GID
          if(gidexists>0)
            newcrosslinkerGID++;
          else
            break;
        }

        /* Create mapping from crosslink molecule to crosslinker element GID
         * Note: No need for the usual procedure of exporting and reimporting to make things redundant
         * because info IS already redundant by design here.*/
        (*crosslink2element_)[i] = newcrosslinkerGID;

        //save positions of nodes between which a crosslinker has to be established in variables xrefe and rotrefe:
        std::vector<double> rotrefe(6);
        std::vector<double> xrefe(6);
        // resize in case of interpolated crosslinker element
        if (linkermodel_ == statmech_linker_stdintpol || linkermodel_ == statmech_linker_activeintpol || linkermodel_ == statmech_linker_bellseqintpol)
          rotrefe.resize(12);

        for(int k=0; k<3; k++)
        {
          xrefe[k ] = (*bspotpositions)[k][bspotcolmap.LID(bspotgid.at(0))];
          xrefe[k+3] = (*bspotpositions)[k][bspotcolmap.LID(bspotgid.at(1))];

          //set nodal rotations (not true ones, only those given in the displacement vector)
          if (linkermodel_ == statmech_linker_stdintpol || linkermodel_ == statmech_linker_activeintpol || linkermodel_ == statmech_linker_bellseqintpol)
          {
            rotrefe[k]   = (*nodalrotations)[k][bspotcolmap.LID((*globalnodeids)[0])];
            rotrefe[k+3] = (*nodalrotations)[k][bspotcolmap.LID((*globalnodeids)[1])];
            rotrefe[k+6] = (*nodalrotations)[k][bspotcolmap.LID((*globalnodeids)[2])];
            rotrefe[k+9] = (*nodalrotations)[k][bspotcolmap.LID((*globalnodeids)[3])];
          }
          else
          {
            rotrefe[k] = (*bspotrotations)[k][bspotcolmap.LID((*globalnodeids)[0])];
            rotrefe[k+3] = (*bspotrotations)[k][bspotcolmap.LID((*globalnodeids)[1])];
          }
        }

        /*a crosslinker is added on each processor which is row node owner of at least one of its nodes;
         *the processor which is row map owner of the node with the larger GID will be the owner of the new crosslinker element; on the other processors the new
         *crosslinker element will be present as a ghost element only*/
        bool hasrownode = false;
        for(int k=0; k<(int)globalnodeids->size(); k++)
          if(noderowmap.LID(globalnodeids->at(k)) > -1)
          {
            hasrownode = true;
            break;
          }

        if(hasrownode)
          AddNewCrosslinkerElement(newcrosslinkerGID, globalnodeids,bspotgid, xrefe,rotrefe,*discret_,nodalquaternions);

        // add all new elements to contact discretization on all Procs
        if(beamcmanager!=Teuchos::null)
          AddNewCrosslinkerElement(newcrosslinkerGID, globalnodeids,bspotgid, xrefe,rotrefe,beamcmanager->ContactDiscret(),nodalquaternions);

        // call of FillComplete() necessary after each added crosslinker because different linkers may share the same nodes (only BeamCL)
//        if(i<addcrosselement.MyLength()-1 && (linkermodel_ == statmech_linker_stdintpol || linkermodel_ == statmech_linker_activeintpol || linkermodel_ == statmech_linker_bellseqintpol))
//        {
//          discret_->CheckFilledGlobally();
//          discret_->FillComplete(true, false, false);
//
//          if(beamcmanager!=Teuchos::null)
//          {
//            beamcmanager->ContactDiscret().CheckFilledGlobally();
//            beamcmanager->ContactDiscret().FillComplete(true, false, false);
//          }
//        }
      }
    }

    // synchronization for problem discretization
    discret_->CheckFilledGlobally();
    discret_->FillComplete(true, false, false);

    // synchronization for contact discretization
    if(beamcmanager!=Teuchos::null)
    {
      beamcmanager->ContactDiscret().CheckFilledGlobally();
      beamcmanager->ContactDiscret().FillComplete(true, false, false);
    }
  }
#ifdef MEASURETIME
  if(!discret_->Comm().MyPID())
  {
    std::cout<<"\n=================Time  Measurement================"<<std::endl;
    std::cout<<"StatMechManager::SearchAndSetCrosslinkers"<<std::endl;
    std::cout<<"Partitioning + Search  :\t"<<std::setprecision(4)<<t1-t0<<"\ts"<<std::endl;
    std::cout<<"Linker administration  :\t"<<std::setprecision(4)<<t2-t1<<"\ts"<<std::endl;
    std::cout<<"Vector communication   :\t"<<std::setprecision(4)<<t3-t2<<"\ts"<<std::endl;
    std::cout<<"Addition of elements   :\t"<<std::setprecision(4)<<Teuchos::Time::wallTime()-t3<<"\ts"<<std::endl;
    std::cout<<"=================================================="<<std::endl;
    std::cout<<"total time             :\t"<<std::setprecision(4)<<Teuchos::Time::wallTime()-t0<<"\ts"<<std::endl;
  }
#endif

  //std::couts
  if(!discret_->Comm().MyPID() && printscreen)
    std::cout<<"\n\n"<<numsetelements<<" crosslinker element(s) added!"<<std::endl;
}//void StatMechManager::SearchAndSetCrosslinkers

/*----------------------------------------------------------------------*
 | create a new crosslinker element and add it to your discretization of|
 | choice (public)                                       mueller (11/11)|
 *----------------------------------------------------------------------*/
void STATMECH::StatMechManager::AddNewCrosslinkerElement(const int&                       crossgid,
                                                         Teuchos::RCP<std::vector<int> >  globalnodeids,
                                                         const std::vector<int>&          bspotgid,
                                                         const std::vector<double>&       xrefe,
                                                         const std::vector<double>&       rotrefe,
                                                         DRT::Discretization&             mydiscret,
                                                         Teuchos::RCP<Epetra_MultiVector> nodalquaternions,
                                                         bool                             addinitlinks)
{
  // get the nodes from the discretization (redundant on both the problem as well as the contact discretization
  // crosslinker is a 4 noded beam Element
  if(statmechparams_.get<double>("ILINK",0.0) > 0.0)
  {
    switch((int)globalnodeids->size())
    {
      case 2:
      {
        // globalnodeids[0] is the larger of the two node GIDs
        Teuchos::RCP<DRT::ELEMENTS::Beam3> newcrosslinker = Teuchos::rcp(new DRT::ELEMENTS::Beam3(crossgid,(mydiscret.gNode(globalnodeids->at(0)))->Owner() ) );

        DRT::Node* nodes[2] = {mydiscret.gNode( globalnodeids->at(0) ), mydiscret.gNode( globalnodeids->at(1) ) };

        newcrosslinker->SetNodeIds((int)globalnodeids->size(),&(globalnodeids->at(0)));
        newcrosslinker->BuildNodalPointers(&nodes[0]);

        //setting up crosslinker element parameters and reference geometry
        newcrosslinker->SetCrossSec(statmechparams_.get<double>("ALINK",0.0));
        newcrosslinker->SetCrossSecShear(1.1*statmechparams_.get<double>("ALINK",0.0));
        newcrosslinker->SetIyy(statmechparams_.get<double>("ILINK",0.0));
        newcrosslinker->SetIzz(statmechparams_.get<double>("ILINK",0.0));
        newcrosslinker->SetIrr(statmechparams_.get<double>("IPLINK",0.0));
        newcrosslinker->SetMaterial(2);
        //set up reference configuration of crosslinker
        newcrosslinker->SetUpReferenceGeometry<2>(xrefe,rotrefe);
        //add element to discretization
        if(!addinitlinks)
        {
          // beam contact discretization
          if(mydiscret.Name()=="beam contact")
          addedcelements_.push_back(crossgid);
          else
            addedelements_.push_back(crossgid);
        }

        mydiscret.AddElement(newcrosslinker);
      }
      break;
      case 4:
      {
        Teuchos::RCP<DRT::ELEMENTS::BeamCL> newcrosslinker = Teuchos::rcp(new DRT::ELEMENTS::BeamCL(crossgid,(mydiscret.gNode(globalnodeids->at(0)))->Owner() ) );
        DRT::Node* nodes[4] = {mydiscret.gNode( globalnodeids->at(0) ),
                               mydiscret.gNode( globalnodeids->at(1) ),
                               mydiscret.gNode( globalnodeids->at(2) ),
                               mydiscret.gNode( globalnodeids->at(3) )};

        newcrosslinker->SetNodeIds((int)globalnodeids->size(),&(globalnodeids->at(0)));
        newcrosslinker->BuildNodalPointers(&nodes[0]);

        //setting up crosslinker element parameters and reference geometry
        newcrosslinker->SetCrossSec(statmechparams_.get<double>("ALINK",0.0));
        newcrosslinker->SetCrossSecShear(1.1*statmechparams_.get<double>("ALINK",0.0));
        newcrosslinker->SetIyy(statmechparams_.get<double>("ILINK",0.0));
        newcrosslinker->SetIzz(statmechparams_.get<double>("ILINK",0.0));
        newcrosslinker->SetIrr(statmechparams_.get<double>("IPLINK",0.0));
        newcrosslinker->SetMaterial(2);
        newcrosslinker->SetBindingPosition((double)(*bspotxi_)[bspotgid[0]],(double)(*bspotxi_)[bspotgid[1]]);

        //set up reference configuration of crosslinker
        newcrosslinker->SetUpReferenceGeometry<2>(xrefe,rotrefe);

        // set initial triads/quaternions
        if(nodalquaternions==Teuchos::null)
          dserror("No nodal quaternions delivered to this method!");

        std::vector<LINALG::Matrix<4,1> > nodequat((int)globalnodeids->size(),LINALG::Matrix<4, 1>(true));
        for(int i=0; i<(int)nodequat.size(); i++)
          for(int j=0; j<(int)nodequat[i].M(); j++)
            (nodequat[i])(j) =(*nodalquaternions)[j][nodalquaternions->Map().LID((*globalnodeids)[i])];

        newcrosslinker->SetInitialQuaternions(nodequat);

        //add element to discretization
        if(!addinitlinks)
        {
          // beam contact discretization
          if(mydiscret.Name()=="beam contact")
            addedcelements_.push_back(crossgid);
          else
            addedelements_.push_back(crossgid);
        }
        mydiscret.AddElement(newcrosslinker);
      }
      break;
      default:
        dserror("Received %i global node IDs! Implemented: 2 or 4!", (int)globalnodeids->size());
    }
  }
  else
  {
    Teuchos::RCP<DRT::ELEMENTS::Truss3> newcrosslinker = Teuchos::rcp(new DRT::ELEMENTS::Truss3(crossgid, (mydiscret.gNode(globalnodeids->at(0)))->Owner() ) );
    DRT::Node* nodes[2] = {mydiscret.gNode( globalnodeids->at(0) ), mydiscret.gNode( globalnodeids->at(1)) };

    newcrosslinker->SetNodeIds(2,&(*globalnodeids)[0]);
    newcrosslinker->BuildNodalPointers(&nodes[0]);

    //setting up crosslinker element parameters
    newcrosslinker ->SetCrossSec(statmechparams_.get<double>("ALINK",0.0));
    newcrosslinker->SetMaterial(2);

    //correct reference configuration data is computed for the new crosslinker element;
    newcrosslinker->SetUpReferenceGeometry(xrefe);

    //add element to discretization
    if(!addinitlinks)
      addedelements_.push_back(crossgid);
    mydiscret.AddElement(newcrosslinker);
  }

  return;
} //void StatMechManager::AddNewCrosslinkerElement()

/*----------------------------------------------------------------------*
 | setting of crosslinkers for loom network setup                       |
 | (private)                                              mueller (2/12)|
 *----------------------------------------------------------------------*/
bool STATMECH::StatMechManager::SetCrosslinkerLoom(Epetra_SerialDenseMatrix& LID, const Teuchos::RCP<Epetra_MultiVector> bspotpositions, const Teuchos::RCP<Epetra_MultiVector> bspottriadscol)
{
  if(networktype_ == statmech_network_casimir)
  {
    bool setcrosslinker = true;
    // 1) if a crosslink between two vertical filaments is considered, do not set a linker
    if((*filamentnumber_)[(int)LID(0,0)]!= 0 && (*filamentnumber_)[(int)LID(1,0)]!=0)
      setcrosslinker = false;
    else
    {
      // 2) consider direction of the crosslinker
      LINALG::Matrix<3,1> crossdir;
      for(int j=0; j<(int)crossdir.M(); j++)
        crossdir(j) = (*bspotpositions)[j][(int)LID(1,0)]-(*bspotpositions)[j][(int)LID(0,0)];
      crossdir.Scale(1.0/crossdir.Norm2());

      // retrieve tangential vector from binding spot quaternions
      std::vector<double> alpha((int)LID.M(),0.0);
      for(int i=0; i<LID.M(); i++)
      {
        LINALG::Matrix<3,1> firstdir;
        LINALG::Matrix<3,3> bspottriad;
        // auxiliary variable for storing a triad in quaternion form
        LINALG::Matrix<4, 1> qnode;
        // triad of node on first filament which is affected by the new crosslinker
        for (int l=0; l<4; l++)
          qnode(l) = (*bspottriadscol)[l][(int)LID(i,0)];
        LARGEROTATIONS::quaterniontotriad(qnode, bspottriad);

        for (int l=0; l<(int)bspottriad.M(); l++)
         firstdir(l) = bspottriad(l,0);
        firstdir.Scale(1.0/firstdir.Norm2());

        alpha.at(i) = acos(fabs(firstdir.Dot(crossdir)));
      }
      // angle: impose orthogonality
      if(alpha.at(0) > 5.0/12.0*M_PI && alpha.at(1) > 5.0/12.0*M_PI)
      {
        // 3) determine which node LID entry is the non-zero entry
        int currfilnumber = -1;
        for(int k=0; k<(int)LID.M(); k++)
        {
          if((int)(*filamentnumber_)[(int)LID(k,0)]!= 0)
          {
            currfilnumber = (int)(*filamentnumber_)[(int)LID(k,0)];
            break;
          }
        }
        for(int k=0; k<filamentnumber_->MyLength(); k++)
        {
          // skip the rest of the nodes when having reached the end of the current filament
          if((*filamentnumber_)[k]==currfilnumber+1)
            break;
          // check if any of the filament's nodes has a linker attached
          // make sure that the binding spot is occupied (i.e. the entry is the linker ID, not "-1")
          if((*bspotstatus_)[k]>-0.1)
          {
            if((*filamentnumber_)[k]==currfilnumber && (*numbond_)[(int)(*bspotstatus_)[k]]>1.9)
            {
              setcrosslinker = false;
              break;
            }
          }
        }
      }
      else
        setcrosslinker = false;
    }

    return setcrosslinker;
  }
  else
    return true;
}

/*----------------------------------------------------------------------*
 |  filament polarity check for active linkers upon bond establishment  |
 |                                             (private)   mueller 08/13|
 *----------------------------------------------------------------------*/
bool STATMECH::StatMechManager::LinkerPolarityCheckAttach(Teuchos::RCP<Epetra_MultiVector> bspottriadscol,
                                                          Epetra_SerialDenseMatrix&      LID,
                                                          LINALG::Matrix<3,1>&           direction)
{
  // cycle time check:
  // time during which an active linker is in recovery conformation (default values from myosin cycle)
  double recoverytime = (statmechparams_.get<double>("ACTIVELINKERCYCLE",0.04))*(statmechparams_.get<double>("ACTIVERECOVERYFRACTION",0.95));
  int crosslid = (*bspotstatus_)[(int)LID(0,0)];
  if((*crosslinkeractcycletime_)[crosslid] >= recoverytime)
  {
    double dirlength = direction.Norm2();
    direction.Scale(1.0/direction.Norm2());

    // 1. polarity check:
    // crosslinker moves to one end (myosin to +end) of the filament.
    // direction vector from substrate filament node to filament node must be in positive direction of the filament

    // direction vector, calculated before:
    // = current position (free filament node) - current position (substrate filament node)
    LINALG::Matrix<3,1> linkdirection(true);
    linkdirection -= direction;

    // unit tangential direction of the filament:
    // first triad vector = tangent
    LINALG::Matrix<3,1> firstdir(true);
    // retrieve tangential and normal vector from binding spot quaternions
    LINALG::Matrix<3,3> bspottriad(true);
    // auxiliary variable for storing a triad in quaternion form
    LINALG::Matrix<4, 1> qnode(true);
    // triad of node on free filament
    for (int l=0; l<4; l++)
      qnode(l) = (*bspottriadscol)[l][(int)LID(1,0)];
    LARGEROTATIONS::quaterniontotriad(qnode, bspottriad);
    for (int l=0; l<(int)bspottriad.M(); l++)
      firstdir(l) = bspottriad(l,0);

    // calculate orientation of direction vector to tangential direction
    double tangentialscalefactor = linkdirection.Dot(firstdir);

    // polarity criterion not fulfilled, if director is not parallel to the free filament
    if (tangentialscalefactor >= 0)
    {
      // 2. polarity criterion (2D): filaments are not allowed to link from to far away
      // alpha = angle between new crosslinker and vertical line --> has to be smaller than pi/2
      double rlink =  statmechparams_.get<double> ("R_LINK", 0.0);
      double alpha = acos(rlink/dirlength);
      if (alpha > statmechparams_.get<double>("PHIBSPOT", 1.0472))
        return false;
    }
    else
      return false;
  }
  else
    return false;

  // all tests passed -> good to go!
  return true;
}

/*----------------------------------------------------------------------*
 | searches crosslinkers and deletes them if probability check is passed|
 | (public)                                               mueller (7/10)|
 *----------------------------------------------------------------------*/
void STATMECH::StatMechManager::SearchAndDeleteCrosslinkers(const int&                                 istep,
                                                            const double&                              timen,
                                                            const double&                              dt,
                                                            const Teuchos::RCP<Epetra_MultiVector>     bspotpositions,
                                                            const Teuchos::RCP<Epetra_MultiVector>     bspotrotations,
                                                            const Epetra_Vector&                       discol,
                                                            Teuchos::RCP<CONTACT::Beam3cmanager>       beamcmanager,
                                                            bool                                       printscreen)
{
//get current off-rate for crosslinkers
  double koff = 0;
  double starttime = actiontime_->at(1);

  if (timen <= starttime || (timen>starttime && fabs(starttime)<dt/1e4))
    koff = statmechparams_.get<double> ("K_OFF_start", 0.0);
  else
    koff = statmechparams_.get<double> ("K_OFF_end", 0.0);

  //probability with which a crosslink breaks up in the current time step
  double p = 1.0 - exp(-dt * koff);

  Teuchos::RCP<Epetra_MultiVector> punlink = Teuchos::rcp(new Epetra_MultiVector(*crosslinkermap_,2));
  punlink->PutScalar(p);

  // binding spot triads
  Teuchos::RCP<Epetra_MultiVector> bspottriadscol = Teuchos::rcp(new Epetra_MultiVector(*bspotcolmap_,4));
  bspottriadscol->PutScalar(0.0);
  GetBindingSpotTriads(bspotrotations, bspottriadscol);

  // if off-rate is also dependent on the forces acting within the crosslinker
  if (statmechparams_.get<double>("DELTABELLSEQ", 0.0) != 0.0)
    ForceDependentOffRate(dt, koff, discol, punlink);

  int numdetachpolarity = 0;
  if(linkermodel_ == statmech_linker_active || linkermodel_ == statmech_linker_activeintpol)
    numdetachpolarity = LinkerPolarityCheckDetach(punlink, bspotpositions, bspottriadscol);

  // a vector indicating the upcoming deletion of crosslinker elements
  Teuchos::RCP<Epetra_Vector> delcrosselement = Teuchos::rcp(new Epetra_Vector(*crosslinkermap_, true));
  delcrosselement->PutScalar(-1.0);

  // SEARCH
  // search and setup for the deletion of elements is done by Proc 0
  int numdelelements = 0;
  if (discret_->Comm().MyPID()==0)
  {
    // create random order of indices
    std::vector<int> order = Permutation(numbond_->MyLength());
    // loop over crosslinkers
    for (int i=0; i<numbond_->MyLength(); i++)
    {
      int irandom = order[i];
      // take action according to the number of bonds of a crosslink molecule
      switch ((int)(*numbond_)[irandom])
      {
        case 0:
          continue;
        break;
        // crosslink molecule with one bond
        case 1:
        {
          // singly bound crosslinker in motility assay is just bond on the substrate filament --> don't delete this bond
          if(networktype_ == statmech_network_motassay)
            continue;

          for (int j=0; j<crosslinkerbond_->NumVectors(); j++)
            if ((*crosslinkerbond_)[j][irandom]>-0.9)
              if ((*uniformgen_)() < (*punlink)[j][irandom])
              {
                // obtain LID and reset crosslinkerbond_ at this position
                int bspotLID = bspotcolmap_->LID((int) (*crosslinkerbond_)[j][irandom]);
                ((*bspotstatus_)[bspotLID]) = -1.0;
                (*numbond_)[irandom] = 0.0;
                (*crosslinkerbond_)[j][irandom] = -1.0;

                LINALG::SerialDenseMatrix LID(1, 1, true);
                LID(0, 0) = bspotLID;
                CrosslinkerIntermediateUpdate(*bspotpositions, LID, irandom, false);
              }
        }
        break;
        // crosslinker element
        case 2:
        {
          // currently, if an element is deleted, a one-bonded crosslink molecule remains (no simultaneous cut of both bonds)
          std::vector<int> jorder = Permutation(crosslinkerbond_->NumVectors());

          int J = crosslinkerbond_->NumVectors();
          if(networktype_ == statmech_network_motassay)
            J = 1;
          for (int j=0; j<J; j++)
          {
            // random pick of one of the crosslinkerbond entries
            int jrandom = jorder[j];
            // in motility assays we have deterministic randomness ;-)
            if(networktype_ == statmech_network_motassay)
              jrandom = 1;
            if ((*uniformgen_)() < (*punlink)[jrandom][irandom])
            {
              (*numbond_)[irandom] = 1.0;

              // an actual crosslinker element exists
              if((*searchforneighbours_)[irandom]>0.9)
              {
                numdelelements++;
                // get the nodal LIDs
                int bspotLID = bspotcolmap_->LID((int)(*crosslinkerbond_)[jrandom][irandom]);

                // in case of a casimir setup, we want the vertical filaments to be free of linkers. Hence, we choose the node on the vertical filament to let go of the linker
                // note: jrandom is either 0 or 1. Hence, "!".
                if(networktype_ == statmech_network_casimir && (*filamentnumber_)[bspotLID]==0)
                {
                  jrandom = (!jrandom);
                  bspotLID = bspotcolmap_->LID((int)(*crosslinkerbond_)[jrandom][irandom]);
                }

                // enter the crosslinker element ID into the deletion list
                (*delcrosselement)[irandom] = (*crosslink2element_)[irandom];

                // vector updates
                (*crosslink2element_)[irandom] = -1.0;
                (*bspotstatus_)[bspotLID] = -1.0;
                (*crosslinkerbond_)[jrandom][irandom] = -1.0;

                // reset linker to long conformation upon unbinding and reset the cycle time
                if(linkermodel_ == statmech_linker_active || linkermodel_ == statmech_linker_activeintpol)
                {
                  (*crosslinkeractlength_)[irandom] = 0.0;
                  (*crosslinkeractcycletime_)[irandom]=0.0;
                }

                // if the linker to be removed occupies to binding spots on the same filament
                if((*crosslinkonsamefilament_)[irandom] > 0.9)
                  (*crosslinkonsamefilament_)[irandom] = 0.0;
              }
              else  // passive crosslink molecule
                (*searchforneighbours_)[irandom] = 1.0;

              for (int k=0; k<crosslinkerbond_->NumVectors(); k++)
                if ((*crosslinkerbond_)[k][irandom]>-0.9)
                {
                  LINALG::SerialDenseMatrix LID(1, 1, true);
                  LID(0,0) = bspotcolmap_->LID((int)(*crosslinkerbond_)[k][irandom]);
                  CrosslinkerIntermediateUpdate(*bspotpositions, LID, irandom);
                  break;
                }

              if(DRT::INPUT::IntegralValue<INPAR::STATMECH::StatOutput>(statmechparams_, "SPECIAL_OUTPUT")==INPAR::STATMECH::statout_structanaly)
              {
                double t_on = (*crosslinkunbindingtimes_)[0][irandom];
                (*crosslinkunbindingtimes_)[1][irandom] = timen-t_on;
                (*crosslinkunbindingtimes_)[0][irandom] = -1.0;
              }

              // leave j-loop after first bond is dissolved
              break;
            }
          }
        }
        break;
      }// switch ((int)(*numbond_)[irandom])
    }// for (int i=0; i<numbond_->MyLength(); i++)
  }// if(discret_->Comm().MyPID()==0)

  // synchronize information about number of bonded filament nodes by exporting it to row map format and then reimporting it to column map format
  // note: searchforneighbours_ and crosslinkonsamefilament_ are not communicated
  // transfer vectors
  Teuchos::RCP<Epetra_Vector> bspotstatusrow = Teuchos::rcp(new Epetra_Vector(*bspotrowmap_, true));
  Teuchos::RCP<Epetra_MultiVector> crosslinkerpositionstrans = Teuchos::rcp(new Epetra_MultiVector(*transfermap_, 3, true));
  Teuchos::RCP<Epetra_MultiVector> crosslinkerbondtrans = Teuchos::rcp(new Epetra_MultiVector(*transfermap_, 2, true));
  Teuchos::RCP<Epetra_Vector> numbondtrans = Teuchos::rcp(new Epetra_Vector(*transfermap_, true));
  Teuchos::RCP<Epetra_Vector> crosslink2elementtrans = Teuchos::rcp(new Epetra_Vector(*transfermap_, true));
  Teuchos::RCP<Epetra_Vector> delcrosselementtrans = Teuchos::rcp(new Epetra_Vector(*transfermap_, true));

  // export and reimport
  CommunicateVector(bspotstatusrow, bspotstatus_);
  CommunicateMultiVector(crosslinkerpositionstrans, crosslinkerpositions_);
  CommunicateMultiVector(crosslinkerbondtrans, crosslinkerbond_);
  CommunicateVector(numbondtrans, numbond_);
  CommunicateVector(crosslink2elementtrans, crosslink2element_);
  CommunicateVector(delcrosselementtrans, delcrosselement);
  if(linkermodel_==statmech_linker_active || linkermodel_ == statmech_linker_activeintpol)
  {
    Teuchos::RCP<Epetra_Vector> crosslinkeractcycletimetrans = Teuchos::rcp(new Epetra_Vector(*transfermap_,true));
    CommunicateVector(crosslinkeractcycletimetrans,crosslinkeractcycletime_);
  }

  // DELETION OF ELEMENTS
  RemoveCrosslinkerElements(*discret_,*delcrosselement,&deletedelements_);
  // contact elements
  if(beamcmanager!=Teuchos::null)
    RemoveCrosslinkerElements(beamcmanager->ContactDiscret(),*delcrosselement,&deletedcelements_);

  if(!discret_->Comm().MyPID() && printscreen)
  {
    std::cout<<numdelelements<<" crosslinker element(s) deleted!"<<std::endl;
    if((linkermodel_==statmech_linker_active || linkermodel_ == statmech_linker_activeintpol) && numdetachpolarity>0)
      std::cout<<"  - thereof "<<numdetachpolarity<<" due to filament polarity!"<<std::endl;
  }
  return;
} //StatMechManager::SearchAndDeleteCrosslinkers()

/*----------------------------------------------------------------------*
 | (private) remove crosslinker element                    mueller 1/11 |
 *----------------------------------------------------------------------*/
void STATMECH::StatMechManager::RemoveCrosslinkerElements(DRT::Discretization& mydiscret, Epetra_Vector& delcrosselement, std::vector<std::vector<char> >* deletedelements)
{
  unsigned startsize = deletedelements->size();
  std::vector<DRT::PackBuffer> vdata( startsize );

  int delel = 0;
  for (int i=0; i<delcrosselement.MyLength(); i++)
  {
    if (mydiscret.HaveGlobalElement((int)delcrosselement[i]))
    {
      //save the element by packing before elimination to make it restorable in case that needed
      vdata.push_back( DRT::PackBuffer() );
      mydiscret.gElement((int)delcrosselement[i])->Pack( vdata.back() );
    }
  }
  for ( unsigned i=startsize; i<vdata.size(); ++i )
    vdata[i].StartPacking();
  for (int i=0; i<delcrosselement.MyLength(); i++)
    if (mydiscret.HaveGlobalElement((int)delcrosselement[i]))
    {
      //save the element by packing before elimination to make it restorable in case that needed
      deletedelements->push_back( std::vector<char>() );
      mydiscret.gElement((int)delcrosselement[i])->Pack( vdata[deletedelements->size()-1] );
      mydiscret.DeleteElement( (int)delcrosselement[i] );
      delel++;
    }
  for ( unsigned i=startsize; i<vdata.size(); ++i )
    swap( (*deletedelements)[i], vdata[i]() );

  /*synchronize
   *the Filled() state on all processors after having added or deleted elements by ChekcFilledGlobally(); then build
   *new element maps and call FillComplete();*/
  mydiscret.CheckFilledGlobally();
  mydiscret.FillComplete(true, false, false);

  return;
}

/*----------------------------------------------------------------------*
 | (private) force dependent off-rate                      mueller 1/11 |
 *----------------------------------------------------------------------*/
void STATMECH::StatMechManager::ForceDependentOffRate(const double&                    dt,
                                                      const double&                    koff0,
                                                      const Epetra_Vector&             discol,
                                                      Teuchos::RCP<Epetra_MultiVector> punlink)
{
  // thermal energy
  double kt = statmechparams_.get<double>("KT",0.00404531);
  // characteristic vond length (=nodereldis)
  // from B. Gui and W. Guilford: Mechanics of actomyosin bonds in different nucleotide states are tuned to muscle contraction
  // Fig 2: slip pathway, ADP, delta = 0.0004;
  // Note: delta<0 -> catch bond, delta>0 -> bond-weakening
  double delta = statmechparams_.get<double>("DELTABELLSEQ", 0.0);

  // update element2crosslink_
  element2crosslink_ = Teuchos::rcp(new Epetra_Vector(*(discret_->ElementRowMap())));
  element2crosslink_->PutScalar(-1.0);

  for(int i=0; i<crosslinkermap_->NumMyElements(); i++)
  {
    if((*crosslink2element_)[i]>-0.9) // there exists a linker element
    {
      // reverse mapping
      int elelid = discret_->ElementRowMap()->LID((int)(*crosslink2element_)[i]);
      if(elelid>-0.9)
        (*element2crosslink_)[elelid] = i;
    }
  }

  for(int i=0; i<discret_->ElementRowMap()->NumMyElements(); i++)
  {
    if((*element2crosslink_)[i]>-0.9)
    {
      DRT::Element* crosslinker = discret_->lRowElement(i);
      int crosslid = crosslinkermap_->LID((int)(*element2crosslink_)[i]);
      int bspotgid0 = (int)(*crosslinkerbond_)[0][crosslid];
      int bspotgid1 = (int)(*crosslinkerbond_)[1][crosslid];

      // GID (be it a node ID or a binding spot ID)
      int checkgid = -1;

      // element internal force vector
      Epetra_SerialDenseVector force;
      // normal strain
      double eps = 0.0;

      const DRT::ElementType &eot = crosslinker->ElementType();
      // retrieve internal force vector
      if(eot == DRT::ELEMENTS::Beam3Type::Instance())
      {
        force.Resize(crosslinker->NumNode()*6);
        force = (dynamic_cast<DRT::ELEMENTS::Beam3*>(crosslinker))->InternalForces();
        eps = (dynamic_cast<DRT::ELEMENTS::Beam3*>(crosslinker))->EpsilonSgn();
        checkgid = crosslinker->Nodes()[0]->Id();
      }
      else if(eot == DRT::ELEMENTS::BeamCLType::Instance())
      {
        if(crosslinker->NumNode()!=4)
          dserror("Currently only implemented for BEAM3CL with four nodes.");
        force.Resize(crosslinker->NumNode()/2*6);
        force = (dynamic_cast<DRT::ELEMENTS::BeamCL*>(crosslinker))->InternalForces();
        eps = (dynamic_cast<DRT::ELEMENTS::BeamCL*>(crosslinker))->EpsilonSgn();

        // mapping from binding spot to internal forces
        // convention (see SearchAndSetCrosslinkers): two subsequent node IDs belong to the same filament
        int crossnodegid0 = crosslinker->Nodes()[0]->Id();
        int crossnodegid1 = crosslinker->Nodes()[1]->Id();

        // node IDs from management vector
        int bspot0nodegid0 = (*bspot2nodes_)[0][bspotcolmap_->LID(bspotgid0)];

        if(bspot0nodegid0 == crossnodegid0 || bspot0nodegid0 == crossnodegid1)
          checkgid = bspotgid0;
        else
          checkgid = bspotgid1;
      }

      // currently, only forces (not moments) considered
      // nodal forces
      LINALG::Matrix<3,1> f0;
      LINALG::Matrix<3,1> f1;

      for(int j=0; j<(int)f0.M(); j++)
      {
        f0(j) = force[j];
        f1(j) = force[6+j];
      }

      // pick the larger absolute force value among the nodes
      double Fbspot0 = f0.Norm2();
      double Fbspot1 = f1.Norm2();
      double sgn = 0.0;
      if (eps<0) sgn = -1.0;
      else if (eps>0) sgn = 1.0;

      // adjusted off-rate according to Bell's equation (Howard, eq 5.10, p.89)
      double koffbspot0 = koff0 * exp(sgn * Fbspot0 * delta/kt);
      double koffbspot1 = koff0 * exp(sgn * Fbspot1 * delta/kt);

      // relate off-rates tp binding spots
      if(bspotgid0 == checkgid)
      {
        (*punlink)[0][crosslid] = 1 - exp(-dt * koffbspot0);
        (*punlink)[1][crosslid] = 1 - exp(-dt * koffbspot1);

        if (unbindingprobability_ != Teuchos::null)
        {
          (*unbindingprobability_)[0][crosslid] = 1 - exp(-dt * koffbspot0);
          (*unbindingprobability_)[1][crosslid] = 1 - exp(-dt * koffbspot1);
        }
      }
      else
      {
        (*punlink)[1][crosslid] = 1 - exp(-dt * koffbspot0);
        (*punlink)[0][crosslid] = 1 - exp(-dt * koffbspot1);
        if (unbindingprobability_ != Teuchos::null)
        {
          (*unbindingprobability_)[0][crosslid] = 1 - exp(-dt * koffbspot0);
          (*unbindingprobability_)[1][crosslid] = 1 - exp(-dt * koffbspot1);
        }
      }

//      FILE* fp = NULL;
//      std::ostringstream filename;
//      filename << outputrootpath_<<"/StatMechOutput/LinkerInternalForces.dat";
//      fp = fopen(filename.str().c_str(), "a");
//      std::stringstream internalforces;
//      internalforces << sgn <<" " << Fnode0 << " " << Fnode1 <<std::endl;
//      fprintf(fp, internalforces.str().c_str());
//      fclose(fp);
//      std::cout<<"\n\nelement "<<crosslinker->Id()<<":"<<std::endl;
//      std::cout<<" bspot "<<bspotgid0<<std::setprecision(8)<<": F = "<<Fbspot0<<", eps = "<<eps<<", k_off = "<<koffbspot0<<", p = "<<(*punlink)[0][crosslid]<<std::endl;
//      std::cout<<" bspot "<<bspotgid1<<std::setprecision(8)<<": F = "<<Fbspot1<<", eps = "<<eps<<", k_off = "<<koffbspot1<<", p = "<<(*punlink)[1][crosslid]<<std::endl;
    }
  }

  // import -> redundancy on all Procs
  Teuchos::RCP<Epetra_MultiVector> punlinktrans = Teuchos::rcp(new Epetra_MultiVector(*transfermap_, true));
  CommunicateMultiVector(punlinktrans, punlink, true, true, false, false);
  if (unbindingprobability_ != Teuchos::null)
  {
    punlinktrans->PutScalar(0.0);
    CommunicateMultiVector(punlinktrans, unbindingprobability_, true, true,false, false);
  }

  return;
} // StatMechManager::ForceDependentOffRate

/*----------------------------------------------------------------------*
 | (private) cahnge reference length of crosslinker elements            |
 |                                                         mueller 08/13|
 *----------------------------------------------------------------------*/
void STATMECH::StatMechManager::ChangeActiveLinkerLength(const double&               timen,
                                                         const double&               dt,
                                                         Teuchos::RCP<Epetra_Vector> actlinklengthin,
                                                         Teuchos::RCP<Epetra_Vector> actlinklengthout,
                                                         bool                        revertchanges)
{
  int numprobshort = 0;
  int numproblong = 0;
  if(!revertchanges)
  {
    // get current hydrolysis-rate for crosslinkers to change it's length to long or short
    double kactlinklong = 0;
    double kactlinkshort = 0;
    double starttime = actiontime_->at(1);

    if(timen <= starttime || (timen>starttime && fabs(timen-starttime) < dt/1e4))
    {
      kactlinklong = statmechparams_.get<double>("K_ACT_LONG_start",0.0);
      kactlinkshort = statmechparams_.get<double>("K_ACT_SHORT_start",0.0);
    }
    else
    {
      kactlinklong = statmechparams_.get<double>("K_ACT_LONG_end",0.0);
      kactlinkshort = statmechparams_.get<double>("K_ACT_SHORT_end",0.0);
    }

    // probability with which a crosslinker change it's length to long
    double plong = 1.0 - exp( -dt*kactlinklong );
    // probability with which a crosslinker change it's length to short
    double pshort = 1.0 - exp( -dt*kactlinkshort );

    /*the following part of the code leads to the decision, which active crosslinker should change it's length.
     * This works precisely as follows:
     *(1) the crosslink molecules are looped through
     *(2) long crosslinkers are checked, if probability to change to short is fulfilled --> if:
     *(3)   beam3-pointer on this crosslinker to change linker length of this crosslinker from long to short
     *(4) do also with short crosslinkers
     *note: a fully overlapping node map on processor 0 is imperative for the following part of the code to work correctly!*/

    Teuchos::RCP<Epetra_Vector> actlinklengthtrans = Teuchos::rcp(new Epetra_Vector(*transfermap_, true));
    CommunicateVector(actlinklengthtrans,actlinklengthout,true,false);

    // probability check, if length of active linker should be changed
    for(int i=0; i<actlinklengthtrans->MyLength(); i++)
    {
      // only for doubly-bound crosslinkers
      if((*crosslink2element_)[i]>-0.9)
      {
        switch((int)(*actlinklengthtrans)[i])
        {
          case 0:
            // probability check for long to short
            if((*uniformgen_)() < pshort)
            {
              numprobshort++;
              (*actlinklengthtrans)[i] = 1.0;
            }
          break;
          case 1:
            // probability check for short to long
            if((*uniformgen_)() < plong)
            {
              numproblong++;
              (*actlinklengthtrans)[i] = 0.0;
            }
          break;
          default: dserror("Wrong status %d in actlinklengthtrans", (int)(*actlinklengthtrans)[i]);
        }
      }
    }
    CommunicateVector(actlinklengthtrans, actlinklengthout, false, true);
  }

  int numprobshortglob = 0;
  int numproblongglob = 0;

  discret_->Comm().SumAll(&numprobshort, &numprobshortglob, 1);
  discret_->Comm().SumAll(&numproblong, &numproblongglob, 1);


  // store linker length change
  int toshort=0;
  int tolong=0;
  // change reference length in actual elements
  for(int i=0; i<actlinklengthout->MyLength(); i++)
  {
    // if length of the active linker changed in the probability check
    if(fabs((*actlinklengthout)[i]-(*actlinklengthin)[i])>1e-6 && (int)(*crosslink2element_)[i]>-0.9)
    {
      // just on the processor, were element is located
      // scaling-factor = 0.75 (myosin II)
      int rowlid = discret_->ElementRowMap()->LID((int)(*crosslink2element_)[i]);
      double sca = statmechparams_.get<double>("LINKERSCALEFACTOR", 0.8);

      switch((int)(*actlinklengthout)[i])
      {
        // change linker length to short
        case 1:
        {
          if (rowlid != -1)
          {
            switch(linkermodel_)
            {
              case statmech_linker_active:
                (dynamic_cast<DRT::ELEMENTS::Beam3*>(discret_->lRowElement(rowlid)))->SetReferenceLength(sca);
              break;
              case statmech_linker_activeintpol:
                (dynamic_cast<DRT::ELEMENTS::BeamCL*>(discret_->lRowElement(rowlid)))->SetReferenceLength(sca);
              break;
              default: dserror("Unknown active linker beam element!");
            }
            toshort++;
          }
        }
        break;
        // change linker length to long
        case 0:
        {
          if (rowlid != -1)
          {
            switch(linkermodel_)
            {
              case statmech_linker_active:
                (dynamic_cast<DRT::ELEMENTS::Beam3*>(discret_->lRowElement(rowlid)))->SetReferenceLength(1.0/sca);
              break;
              case statmech_linker_activeintpol:
                (dynamic_cast<DRT::ELEMENTS::BeamCL*>(discret_->lRowElement(rowlid)))->SetReferenceLength(1.0/sca);
              break;
              default: dserror("Unknown active linker beam element!");
            }
            tolong++;
          }
        }
        break;
        default: dserror("Wrong status %d in actlinklength_", (int)(*crosslinkeractlength_)[i]);
      }
    }
  }

  int toshortglob = 0;
  int tolongglob = 0;

  discret_->Comm().SumAll(&toshort, &toshortglob, 1);
  discret_->Comm().SumAll(&tolong, &tolongglob,1);

  // number of changed crosslinker lengths
  if(!discret_->Comm().MyPID())
  {
    if(revertchanges)
      std::cout<<"\nRevert reference lengths..."<<std::endl;
    else
    {
      if((toshort+tolong)>0)
      {
        std::cout<<"\n"<<"--Crosslinker length changes--"<<std::endl;
        std::cout<<" - long to short: "<<toshort<<std::endl;
        std::cout<<"   - due to prob: "<<numprobshortglob<<std::endl;
        std::cout<<" - short to long: "<<tolong<<std::endl;
        std::cout<<"   - due to prob: "<<numproblongglob<<std::endl;
        std::cout<<"------------------------------"<<std::endl;
      }
    }

    // update cycle time
    for(int i=0; i<crosslinkeractcycletime_->MyLength(); i++)
      (*crosslinkeractcycletime_)[i] += dt;
  }

  Teuchos::RCP<Epetra_Vector> crosslinkeractcycletimetrans = Teuchos::rcp(new Epetra_Vector(*transfermap_,true));
  CommunicateVector(crosslinkeractcycletimetrans,crosslinkeractcycletime_);

  return;
}

/*----------------------------------------------------------------------*
 |  filament polarity check for active linkers upon bond breaking       |
 |                                               (private) mueller 08/13|
 *----------------------------------------------------------------------*/
int STATMECH::StatMechManager::LinkerPolarityCheckDetach(Teuchos::RCP<Epetra_MultiVector>       punlink,
                                                          const Teuchos::RCP<Epetra_MultiVector> bspotpositions,
                                                          const Teuchos::RCP<Epetra_MultiVector> bspottriadscol)
{

  int numdetach = 0;

  if(!discret_->Comm().MyPID())
  {
    for(int i=0; i<crosslinkermap_->NumMyElements(); i++)
    {
      // there exists a crosslinker
      if((*crosslink2element_)[i]>-0.9)
      {
        // polarity check:
        // crosslinker goes to one end (myosin to +end) of the filament
        // direction vector from substrate filament node to filament node must be in positive direction
        // of the filament

        //Col map LIDs of nodes, which are crosslinked
        int bspotlid0 = bspotcolmap_->LID((int)(*crosslinkerbond_)[0][i]);
        int bspotlid1 = bspotcolmap_->LID((int)(*crosslinkerbond_)[1][i]);

        //direction vector between currently considered two nodes
        // = current position (free filament node) - current position (substrate filament node)
        LINALG::Matrix<3,1> linkdirection(true);
        for(int j=0; j<(int)linkdirection.M(); j++)
          linkdirection(j) = (*bspotpositions)[j][bspotlid1] - (*bspotpositions)[j][bspotlid0];

        // unit tangential direction of the filament:
        // first triad vector = tangent
        LINALG::Matrix<3,1> firstdir(true);
        // retrieve tangential and normal vector from binding spot quaternions
        LINALG::Matrix<3,3> bspottriad(true);
        // auxiliary variable for storing a triad in quaternion form
        LINALG::Matrix<4, 1> qnode(true);
        // triad of node on free filament
        for (int l=0; l<4; l++)
          qnode(l) = (*bspottriadscol)[l][bspotlid1];
        LARGEROTATIONS::quaterniontotriad(qnode, bspottriad);
        for (int l=0; l<(int)bspottriad.M(); l++)
          firstdir(l) = bspottriad(l,0);

        // calculate orientation of direction vector to tangential direction
        double tangentialscalefactor = linkdirection.Dot(firstdir);
        // polarity not fulfilled, if direction vector is not in the direction of the free filament
        //TODO MULTIVECTOR -> i.e. accounting for both binding sites!!!
        // For now, manipulate j==1 since this is the free filament binding site
        if (tangentialscalefactor < 0.0)
        {
          numdetach++;
          (*punlink)[1][i]=1.0;
        }
      }
    }
  }

  Teuchos::RCP<Epetra_MultiVector> punlinktrans = Teuchos::rcp(new Epetra_MultiVector(*transfermap_,2,true));
  CommunicateMultiVector(punlinktrans,punlink);

  return numdetach;
}// StatMechManager::LinkerPolarityCheckDetach()

/*----------------------------------------------------------------------*
 | (private) Reduce currently existing crosslinkers by a certain        |
 | percentage.                                              mueller 1/11|
 *----------------------------------------------------------------------*/
void STATMECH::StatMechManager::ReduceNumOfCrosslinkersBy(const int                            numtoreduce,
                                                          Teuchos::RCP<CONTACT::Beam3cmanager> beamcmanager)
{
  if(!discret_->Comm().MyPID())
  {
    std::cout<<"\n\n==========================================================="<<std::endl;
    std::cout<<"-- "<<crosslinkermap_->NumMyElements()<<" crosslink molecules in volume"<<std::endl;
    std::cout<<"-- removing "<<statmechparams_.get<int>("REDUCECROSSLINKSBY",0)<<" crosslinkers...\n"<<std::endl;
  }

  int ncrosslink = statmechparams_.get<int>("N_crosslink", 0);

  // check for the correctness of the given input value
  if(numtoreduce>ncrosslink)
    dserror("REDUCECROSSLINKSBY is greater than N_crosslink. Please check your input file!");

  // a vector indicating the upcoming deletion of crosslinker elements
  Teuchos::RCP<Epetra_Vector> delcrosselement = Teuchos::rcp(new Epetra_Vector(*crosslinkermap_));
  Teuchos::RCP<Epetra_Vector> delcrossmolecules = Teuchos::rcp(new Epetra_Vector(*crosslinkermap_, true));
  delcrosselement->PutScalar(-1.0);

  // SEARCH
  // search and setup for the deletion of elements is done by Proc 0
  int numdelelements = 0;
  int numdelmolecules = 0;
  if (discret_->Comm().MyPID()==0)
  {
    //create random order in which crosslinkers are addressed
    std::vector<int> randomorder = Permutation(ncrosslink);
    // number of deleted crosslinkers
    int numdelcrosslinks = 0;
    for (int i=0; i<numbond_->MyLength(); i++)
    {
      int irandom = randomorder[i];
      if(numdelcrosslinks<numtoreduce)
      {
        // take action according to the number of bonds of a crosslink molecule
        switch ((int)(*numbond_)[irandom])
        {
          case 0:
          {
            numdelmolecules++;
            (*delcrossmolecules)[irandom] = 1.0;
          }
          break;
          // crosslink molecule with one bond
          case 1:
          {
            numdelmolecules++;
            for (int j=0; j<crosslinkerbond_->NumVectors(); j++)
              if ((*crosslinkerbond_)[j][irandom]>-0.9)
              {
                // obtain LID and reset crosslinkerbond_ at this position
                int bspotLID = bspotcolmap_->LID((int) (*crosslinkerbond_)[j][irandom]);
                (*bspotstatus_)[bspotLID] = -1.0;
                (*delcrossmolecules)[irandom] = 1.0;
              }
          }
          break;
          // crosslinker element
          case 2:
          {
            // an actual crosslinker element exists
            if((*searchforneighbours_)[irandom]>0.9)
            {
              numdelelements++;
              // enter the crosslinker element ID into the deletion list
              (*delcrosselement)[irandom] = (*crosslink2element_)[irandom];
              (*delcrossmolecules)[irandom] = 1.0;

              ((*bspotstatus_)[(int)(*crosslinkerbond_)[0][irandom]]) = -1.0;
              ((*bspotstatus_)[(int)(*crosslinkerbond_)[1][irandom]]) = -1.0;
            }
            else  // passive crosslink molecule
            {
              numdelmolecules++;
              for (int j=0; j<crosslinkerbond_->NumVectors(); j++)
                if ((*crosslinkerbond_)[j][irandom]>-0.9)
                {
                  // obtain LID and reset crosslinkerbond_ at this position
                  int bspotLID = bspotcolmap_->LID((int)(*crosslinkerbond_)[j][irandom]);
                  ((*bspotstatus_)[bspotLID]) = -1.0;
                  (*delcrossmolecules)[irandom] = 1.0;
                }
            }
          }
          break;
        }// switch ((int)(*numbond_)[i])
        numdelcrosslinks++;
      }// if(numdelcrosslinks<numtoreduce)
      else // upon reaching the given number of deleted crosslinks, exit the loop
        break;
    }// for (int i=0; i<numbond_->MyLength(); i++)
  }// if(discret_->Comm().MyPID()==0)

  //synchronize information about number of bonded filament nodes by exporting it to row map format and then reimporting it to column map format
  // transfer vector
  Teuchos::RCP<Epetra_Vector> bspotstatusrow = Teuchos::rcp(new Epetra_Vector(*bspotrowmap_, true));
  Teuchos::RCP<Epetra_Vector> delcrosselementtrans = Teuchos::rcp(new Epetra_Vector(*transfermap_, true));
  Teuchos::RCP<Epetra_Vector> delcrossmoleculestrans = Teuchos::rcp(new Epetra_Vector(*transfermap_, true));
  //export and reimport
  CommunicateVector(bspotstatusrow, bspotstatus_);
  CommunicateVector(delcrosselementtrans, delcrosselement);
  CommunicateVector(delcrossmoleculestrans, delcrossmolecules);

  //PREPARE REBUILDING OF CROSSLINKER MAPS AND CORRESPONDING VECTORS
  // std::vectors for transfer of data to new maps and vectors
  std::vector<std::vector<double> > newcrosslinkerpositions;
  std::vector<std::vector<double> > newcrosslinkerbond;
  std::vector<std::vector<double> > newvisualizepositions;
  std::vector<std::vector<double> > newcrosslinkunbindingtimes;
  std::vector<int> newcrosslinkonsamefilament;
  std::vector<int> newsearchforneighbours;
  std::vector<int> newnumbond;
  std::vector<int> newcrosslink2element;
  // Build temporary vectors
  for(int i=0; i<delcrossmolecules->MyLength(); i++)
  {
    // if crosslinker is kept
    if((*delcrossmolecules)[i]<0.1)
    {
      // add the crosslinkers which are kept, to temporary vectors
      std::vector<double> crosspos;
      std::vector<double> visualpos;
      std::vector<double> crossbond;
      for(int j=0; j<crosslinkerpositions_->NumVectors(); j++)
        crosspos.push_back((*crosslinkerpositions_)[j][i]);
      for(int j=0; j<visualizepositions_->NumVectors(); j++)
        visualpos.push_back((*visualizepositions_)[j][i]);
      for(int j=0; j<crosslinkerbond_->NumVectors(); j++)
        crossbond.push_back((*crosslinkerbond_)[j][i]);

      if(DRT::INPUT::IntegralValue<INPAR::STATMECH::StatOutput>(statmechparams_, "SPECIAL_OUTPUT")==INPAR::STATMECH::statout_structanaly)
      {
        std::vector<double> unbind;
        for(int j=0; j<crosslinkunbindingtimes_->NumVectors(); j++)
          unbind.push_back((*crosslinkunbindingtimes_)[j][i]);
        newcrosslinkunbindingtimes.push_back(unbind);
      }

      newcrosslinkerpositions.push_back(crosspos);
      newvisualizepositions.push_back(visualpos);
      newcrosslinkerbond.push_back(crossbond);
      newcrosslinkonsamefilament.push_back((int)(*crosslinkonsamefilament_)[i]);
      newsearchforneighbours.push_back((int)(*searchforneighbours_)[i]);
      newnumbond.push_back((int)(*numbond_)[i]);
      newcrosslink2element.push_back((int)(*crosslink2element_)[i]);
    }
  }

  unsigned startsize = deletedelements_.size();
  std::vector<DRT::PackBuffer> vdata( startsize );

  // DELETION OF ELEMENTS
  for (int i=0; i<delcrosselement->MyLength(); i++)
    if (discret_->HaveGlobalElement((int)(*delcrosselement)[i]))
    {
      //save the element by packing before elimination to make it restorable in case that needed
      //std::cout<<"Proc "<<discret_->Comm().MyPID()<<": deleting element
      //"<<(int)deletecrosslinkerelements[i]<<std::endl;
      vdata.push_back( DRT::PackBuffer() );
      discret_->gElement((int)(*delcrosselement)[i])->Pack( vdata.back() );
    }
  for ( unsigned i=startsize; i<vdata.size(); ++i )
    vdata[i].StartPacking();
  for (int i=0; i<delcrosselement->MyLength(); i++)
    if (discret_->HaveGlobalElement((int)(*delcrosselement)[i]))
    {
      //save the element by packing before elimination to make it restorable in case that needed
      //std::cout<<"Proc "<<discret_->Comm().MyPID()<<": deleting element "<<(int)deletecrosslinkerelements[i]<<std::endl;
      deletedelements_.push_back( std::vector<char>() );
      discret_->gElement((int)(*delcrosselement)[i])->Pack( vdata[deletedelements_.size()-1] );
      discret_->DeleteElement( (int)(*delcrosselement)[i]);
    }
  for ( unsigned i=startsize; i<vdata.size(); ++i )
    swap( deletedelements_[i], vdata[i]() );
  /*synchronize
  *the Filled() state on all processors after having added or deleted elements by ChekcFilledGlobally(); then build
  *new element maps and call FillComplete();*/
  discret_->CheckFilledGlobally();
  discret_->FillComplete(true, false, false);

  // SET UP OF NEW CROSSLINKER MAPS AND VECTORS
  // create crosslinker maps
  std::vector<int> newgids;
  for (int i=0; i<(int)newcrosslinkerpositions.size(); i++)
    newgids.push_back(i);

  // crosslinker column and row maps
  crosslinkermap_ = Teuchos::rcp(new Epetra_Map(-1, (int)newgids.size(), &newgids[0], 0, discret_->Comm()));
  transfermap_    = Teuchos::rcp(new Epetra_Map((int)newgids.size(), 0, discret_->Comm()));
  //vectors
  crosslinkerpositions_ = Teuchos::rcp(new Epetra_MultiVector(*crosslinkermap_, 3));
  visualizepositions_ = Teuchos::rcp(new Epetra_MultiVector(*crosslinkermap_, 3));
  crosslinkerbond_ = Teuchos::rcp(new Epetra_MultiVector(*crosslinkermap_, 2));
  crosslinkonsamefilament_ = Teuchos::rcp(new Epetra_Vector(*crosslinkermap_));
  searchforneighbours_ = Teuchos::rcp(new Epetra_Vector(*crosslinkermap_));
  numbond_ = Teuchos::rcp(new Epetra_Vector(*crosslinkermap_));
  crosslink2element_ = Teuchos::rcp(new Epetra_Vector(*crosslinkermap_));
  if(DRT::INPUT::IntegralValue<INPAR::STATMECH::StatOutput>(statmechparams_, "SPECIAL_OUTPUT")==INPAR::STATMECH::statout_structanaly)
    crosslinkunbindingtimes_ = Teuchos::rcp(new Epetra_MultiVector(*crosslinkermap_,3));

  //copy information from the temporary vectors to the adjusted crosslinker vectors
  for(int i=0; i<crosslinkerpositions_->MyLength(); i++)
  {
    for(int j=0; j<crosslinkerpositions_->NumVectors(); j++)
      (*crosslinkerpositions_)[j][i] = newcrosslinkerpositions[i][j];
    for(int j=0; j<visualizepositions_->NumVectors(); j++)
      (*visualizepositions_)[j][i] = (double)newvisualizepositions[i][j];
    for(int j=0; j<crosslinkerbond_->NumVectors(); j++)
      (*crosslinkerbond_)[j][i] = (double)newcrosslinkerbond[i][j];

    if(DRT::INPUT::IntegralValue<INPAR::STATMECH::StatOutput>(statmechparams_, "SPECIAL_OUTPUT")==INPAR::STATMECH::statout_structanaly)
      for(int j=0; j<crosslinkerbond_->NumVectors(); j++)
        (*crosslinkunbindingtimes_)[j][i] = (double)newcrosslinkunbindingtimes[i][j];

    (*crosslinkonsamefilament_)[i] = (double)newcrosslinkonsamefilament[i];
    (*searchforneighbours_)[i] = (double)newsearchforneighbours[i];
    (*numbond_)[i] = (double)newnumbond[i];
    (*crosslink2element_)[i] = (double)newcrosslink2element[i];
  }
  if(!discret_->Comm().MyPID())
  {
    std::cout<<"-- "<<numdelelements<<" crosslinker elements removed"<<std::endl;
    std::cout<<"-- "<<numdelmolecules<<" free/one-bonded crosslinker molecules removed"<<std::endl;
    std::cout<<"------------------------------------------------------"<<std::endl;
    std::cout<<"-- "<<numdelelements+numdelmolecules<<" crosslinkers removed"<<std::endl;
    std::cout<<"\n-- "<<crosslinkermap_->NumMyElements()<<" crosslink molecules left in volume"<<std::endl;
    std::cout<<"===========================================================\n"<<std::endl;
  }
  // ReduceNumberOfCrosslinkersBy() is called at the beginning of the time step. Hence, we can just call WriteConv() again to update the converged state
  WriteConv(beamcmanager);

  return;
}//ReduceNumberOfCrosslinkersBy()


/*----------------------------------------------------------------------*
 | (public) generate gaussian randomnumbers with mean "meanvalue" and   |
 | standarddeviation "standarddeviation" for parallel use     cyron10/09|
 *----------------------------------------------------------------------*/
void STATMECH::StatMechManager::GenerateGaussianRandomNumbers(Teuchos::RCP<Epetra_MultiVector> randomnumbers, const double meanvalue, const double standarddeviation)
{
  randomnumbers->PutScalar(0.0);

  //multivector for stochastic forces evaluated by each element based on row map
  Teuchos::RCP<Epetra_MultiVector> randomnumbersrow = Teuchos::rcp(new Epetra_MultiVector(*(discret_->ElementRowMap()), randomnumbers->NumVectors()));

  for (int i=0; i<randomnumbersrow->MyLength(); i++)
    for (int j=0; j<randomnumbersrow->NumVectors(); j++)
      (*randomnumbersrow)[j][i] = standarddeviation*(*normalgen_)() + meanvalue;

  //export stochastic forces from row map to column map (unusual CommunicateMultiVector() call but does the job!)
  CommunicateMultiVector(randomnumbers,randomnumbersrow,true,false,false);

#ifdef DEBUGCOUT
  for(int pid=0; pid<discret_->Comm().NumProc(); pid++)
  {
    if(pid==discret_->Comm().MyPID())
      std::cout<<"\n\nProc "<<pid<<": Row\n\n"<<*randomnumbersrow<<std::endl;
    discret_->Comm().Barrier();
  }

  discret_->Comm().Barrier();

  for(int pid=0; pid<discret_->Comm().NumProc(); pid++)
  {
    if(pid==discret_->Comm().MyPID())
      std::cout<<"\n\nProc "<<pid<<": Column\n\n"<<*randomnumbers<<std::endl;
    discret_->Comm().Barrier();
  }
#endif

  return;
} // StatMechManager::SynchronizeRandomForces()

/*----------------------------------------------------------------------*
 | (public) writing restart information for manager objects   cyron 12/08|
 *----------------------------------------------------------------------*/
void STATMECH::StatMechManager::WriteRestart(Teuchos::RCP<IO::DiscretizationWriter> output, double& dt)
{
  output->WriteInt("istart", istart_);
  output->WriteInt("unconvergedsteps", unconvergedsteps_);
  output->WriteDouble("starttimeoutput", starttimeoutput_);
  output->WriteDouble("endtoendref", endtoendref_);
  output->WriteInt("basisnodes", basisnodes_);
  output->WriteInt("outputfilenumber", outputfilenumber_);
  //note: beginold_,endold_,sumdispmiddle_ not considered; related methods not restartable
  output->WriteInt("basiselements", basiselements_);
  output->WriteDouble("sumsquareincpar", sumsquareincpar_);
  output->WriteDouble("sumsquareincort", sumsquareincort_);
  output->WriteDouble("sumrotmiddle", sumrotmiddle_);
  output->WriteDouble("sumsquareincmid", sumsquareincmid_);
  output->WriteDouble("sumsquareincrot", sumsquareincrot_);
  output->WriteDouble("timestepsize", dt);
  /*note: crosslinkermap_, transfermap_, ddcorrrowmap_, ddcorrcolmap_,
   * filamentnumber_, forcesensor_, not considered, because generated in constructor*/

  //note: Using WriteVector and ReadMultiVector requires unique map of MultiVector thus export/import for restart to/from row map
  Epetra_Export exporter(*bspotcolmap_,*bspotrowmap_);
  Teuchos::RCP<Epetra_Vector> bspotstatusrow = Teuchos::rcp(new Epetra_Vector(*bspotrowmap_,true));
  bspotstatusrow->Export(*bspotstatus_,exporter,Insert);
  output->WriteVector("bspotstatus",bspotstatusrow,IO::DiscretizationWriter::nodevector);

  output->WriteRedundantDoubleVector("startindex",startindex_);

  WriteRestartRedundantMultivector(output,"crosslinkerbond",crosslinkerbond_);
  WriteRestartRedundantMultivector(output,"crosslinkerpositions",crosslinkerpositions_);
  WriteRestartRedundantMultivector(output,"numbond",numbond_);
  WriteRestartRedundantMultivector(output,"crosslink2element",crosslink2element_);
  WriteRestartRedundantMultivector(output,"crosslinkonsamefilament",crosslinkonsamefilament_);
  WriteRestartRedundantMultivector(output,"visualizepositions",visualizepositions_);
  WriteRestartRedundantMultivector(output,"searchforneighbours",searchforneighbours_);

  if(crosslinkunbindingtimes_!=Teuchos::null)
    WriteRestartRedundantMultivector(output,"crosslinkunbindingtimes",crosslinkunbindingtimes_);



  return;
} // StatMechManager::WriteRestart()

/*----------------------------------------------------------------------------*
 | (public) write restart information for fully redundant   Epetra_Multivector|
 | with name "name"                                                cyron 11/10|
 *----------------------------------------------------------------------------*/
void STATMECH::StatMechManager::WriteRestartRedundantMultivector(Teuchos::RCP<IO::DiscretizationWriter> output, const std::string name, RCP<Epetra_MultiVector> multivector)
{
  //create stl vector to store information in multivector
  Teuchos::RCP<std::vector<double> > stlvector = Teuchos::rcp(new std::vector<double>);
  stlvector->resize(multivector->MyLength()*multivector->NumVectors());

  for (int i=0; i<multivector->MyLength(); i++)
    for (int j=0; j<multivector->NumVectors(); j++)
      (*stlvector)[i + j*multivector->MyLength()] = (*multivector)[j][i];

  //write information to output file; note that WriteRedundantDoubleVector is active on proc 0 only
  output->WriteRedundantDoubleVector(name,stlvector);

  return;
} // StatMechManager::WriteRestartRedundantMultivector()



/*----------------------------------------------------------------------*
 |read restart information for statistical mechanics (public)cyron 12/08|
 *----------------------------------------------------------------------*/
void STATMECH::StatMechManager::ReadRestart(IO::DiscretizationReader& reader, double& dt)
{

  // read restart information for statistical mechanics
  istart_ = reader.ReadInt("istart");
  unconvergedsteps_ = reader.ReadInt("unconvergedsteps");
  starttimeoutput_ = reader.ReadDouble("starttimeoutput");
  endtoendref_ = reader.ReadDouble("endtoendref");
  basisnodes_ = reader.ReadInt("basisnodes");
  outputfilenumber_ = reader.ReadInt("outputfilenumber");
  //note: beginold_,endold_,sumdispmiddle_ not considered; related methods not restartable
  basiselements_ = reader.ReadInt("basiselements");
  sumsquareincpar_ = reader.ReadDouble("sumsquareincpar");
  sumsquareincort_ = reader.ReadDouble("sumsquareincort");
  sumrotmiddle_ = reader.ReadDouble("sumrotmiddle");
  sumsquareincmid_ = reader.ReadDouble("sumsquareincmid");
  sumsquareincrot_ = reader.ReadDouble("sumsquareincrot");
  dt = reader.ReadDouble("timestepsize");
  /*note: crosslinkermap_, transfermap_, ddcorrrowmap_, ddcorrcolmap_,
   * filamentnumber_, forcesensor_, not considered, because generated in constructor*/

  //note: Using WriteVector and ReadMultiVector requires uniquen map of MultiVector thus export/import for restart to/from row map
  Epetra_Import importer(*bspotcolmap_,*bspotrowmap_);
  Teuchos::RCP<Epetra_Vector> bspotstatusrow = Teuchos::rcp(new Epetra_Vector(*bspotrowmap_),true);
  reader.ReadVector(bspotstatusrow,"bspotstatus");
  bspotstatus_->Import(*bspotstatusrow,importer,Insert);


  //Read redundant Epetra_Multivectors and STL vectors
  reader.ReadRedundantDoubleVector(startindex_,"startindex");
  ReadRestartRedundantMultivector(reader,"crosslinkerbond",crosslinkerbond_);
  ReadRestartRedundantMultivector(reader,"crosslinkerpositions",crosslinkerpositions_);
  ReadRestartRedundantMultivector(reader,"numbond",numbond_);
  ReadRestartRedundantMultivector(reader,"crosslink2element",crosslink2element_);
  ReadRestartRedundantMultivector(reader,"crosslinkonsamefilament",crosslinkonsamefilament_);
  ReadRestartRedundantMultivector(reader,"visualizepositions",visualizepositions_);
  ReadRestartRedundantMultivector(reader,"searchforneighbours",searchforneighbours_);
  if(DRT::INPUT::IntegralValue<INPAR::STATMECH::StatOutput>(statmechparams_, "SPECIAL_OUTPUT")==INPAR::STATMECH::statout_structanaly)
    ReadRestartRedundantMultivector(reader,"crosslinkunbindingtimes",crosslinkunbindingtimes_);


  return;
}// StatMechManager::ReadRestart()

/*----------------------------------------------------------------------------*
 | (public) read restart information for fully redundant Epetra_Multivector   |
 | with name "name"                                                cyron 11/10|
 *----------------------------------------------------------------------------*/
void STATMECH::StatMechManager::ReadRestartRedundantMultivector(IO::DiscretizationReader& reader, const std::string name, Teuchos::RCP<Epetra_MultiVector> multivector)
{
    //we assume that information was stored like for a redundant stl vector
    Teuchos::RCP<std::vector<double> > stlvector = Teuchos::rcp(new std::vector<double>);
    stlvector->resize(multivector->MyLength()*multivector->NumVectors());

    /*ReadRedundantDoubleVector reads information of stlvector on proc 0 and distributes
     *this information then to all other processors so that after the following line
     *the variable stlvector has the same information on all procs*/
    reader.ReadRedundantDoubleVector(stlvector,name);

    //transfer data from stlvector to Epetra_Multivector
    for (int i=0; i<multivector->MyLength(); i++)
      for (int j=0; j<multivector->NumVectors(); j++)
         (*multivector)[j][i] = (*stlvector)[i + j*multivector->MyLength()];

  return;
} // StatMechManager::WriteRestartRedundantMultivector()

/*-----------------------------------------------------------------------*
 | (public) saves all relevant variables *_ as *conv_ to allow  for      |
 | returning to the beginning of a time step                 cyron 11/10 |
 *-----------------------------------------------------------------------*/
void STATMECH::StatMechManager::WriteConv(Teuchos::RCP<CONTACT::Beam3cmanager> beamcmanager)
{
  //save relevant class variables at the very end of the time step
  crosslinkerbondconv_ = Teuchos::rcp(new Epetra_MultiVector(*crosslinkerbond_));
  crosslinkerpositionsconv_ = Teuchos::rcp(new Epetra_MultiVector(*crosslinkerpositions_));
  bspotstatusconv_ = Teuchos::rcp(new Epetra_Vector(*bspotstatus_));
  numbondconv_ = Teuchos::rcp(new Epetra_Vector(*numbond_));
  crosslinkonsamefilamentconv_ = Teuchos::rcp(new Epetra_Vector(*crosslinkonsamefilament_));
  crosslink2elementconv_ = Teuchos::rcp(new Epetra_Vector(*crosslink2element_));
  searchforneighboursconv_ = Teuchos::rcp(new Epetra_Vector(*searchforneighbours_));
  if(linkermodel_ == statmech_linker_active || linkermodel_ == statmech_linker_activeintpol)
    crosslinkeractlengthconv_ = Teuchos::rcp(new Epetra_Vector(*crosslinkeractlength_));

  //set addedelements_, deletedelements_ empty vectors
  addedelements_.clear();
  deletedelements_.clear();
  if(beamcmanager!=Teuchos::null)
  {
    addedcelements_.clear();
    deletedcelements_.clear();
  }

  return;
} // StatMechManager::WriteConv()

/*-----------------------------------------------------------------------*
 | communicate Vector to all Processors                    mueller 11/11 |
 *-----------------------------------------------------------------------*/
void STATMECH::StatMechManager::CommunicateVector(Teuchos::RCP<Epetra_Vector> InVec,
                                                  Teuchos::RCP<Epetra_Vector> OutVec,
                                                  bool                        doexport,
                                                  bool                        doimport,
                                                  bool                        zerofy,
                                                  bool                        exportinsert)
{
  /* zerofy InVec at the beginning of each search except for Proc 0
   * for subsequent export and reimport. This way, we guarantee redundant information
   * on all processors. */

  const Epetra_BlockMap& InMap = InVec->Map();
  const Epetra_BlockMap& OutMap = OutVec->Map();

  // first, export the values of OutVec on Proc 0 to InVecs of all participating processors
  Epetra_Export exporter(OutMap, InMap);
  Epetra_Import importer(OutMap, InMap);
  if(doexport)
  {
    // zero out all vectors which are not Proc 0. Then, export Proc 0 data to InVec map.
    if(discret_->Comm().MyPID()!=0 && zerofy)
      OutVec->PutScalar(0.0);
    if(exportinsert)
      InVec->Export(*OutVec, exporter, Insert);
    else
      InVec->Export(*OutVec, exporter, Add);
  }
  if(doimport)
    OutVec->Import(*InVec,importer,Insert);
  return;
}

/*-----------------------------------------------------------------------*
 | communicate MultiVector to all Processors               mueller 11/11 |
 *-----------------------------------------------------------------------*/
void STATMECH::StatMechManager::CommunicateMultiVector(Teuchos::RCP<Epetra_MultiVector> InVec,
                                                       Teuchos::RCP<Epetra_MultiVector> OutVec,
                                                       bool                             doexport,
                                                       bool                             doimport,
                                                       bool                             zerofy,
                                                       bool                             exportinsert)
{
  // first, export the values of OutVec on Proc 0 to InVecs of all participating processors
  Epetra_Export exporter(OutVec->Map(), InVec->Map());
  Epetra_Import importer(OutVec->Map(), InVec->Map());

  if(doexport)
  {
    // zero out all vectors which are not Proc 0. Then, export Proc 0 data to InVec map.
    if(discret_->Comm().MyPID()!=0 && zerofy)
      OutVec->PutScalar(0.0);
    if(exportinsert)
      InVec->Export(*OutVec, exporter, Insert);
    else
      InVec->Export(*OutVec, exporter, Add);
  }
  if(doimport)
    OutVec->Import(*InVec,importer,Insert);
  return;
}

/*-----------------------------------------------------------------------*
 | (public) restore state at the beginning of this time step cyron 11/10 |
 *-----------------------------------------------------------------------*/
void STATMECH::StatMechManager::RestoreConv(Teuchos::RCP<LINALG::SparseOperator>& stiff, Teuchos::RCP<CONTACT::Beam3cmanager> beamcmanager)
{
  //restore state at the beginning of time step for relevant class variables
  crosslinkerbond_ = Teuchos::rcp(new Epetra_MultiVector(*crosslinkerbondconv_));
  crosslinkerpositions_ = Teuchos::rcp(new Epetra_MultiVector(*crosslinkerpositionsconv_));
  bspotstatus_ = Teuchos::rcp(new Epetra_Vector(*bspotstatusconv_));
  numbond_ = Teuchos::rcp(new Epetra_Vector(*numbondconv_));
  crosslinkonsamefilament_ = Teuchos::rcp(new Epetra_Vector(*crosslinkonsamefilamentconv_));
  crosslink2element_ = Teuchos::rcp(new Epetra_Vector(*crosslink2elementconv_));
  searchforneighbours_ = Teuchos::rcp(new Epetra_Vector(*searchforneighboursconv_));
  if(linkermodel_ == statmech_linker_active || linkermodel_ == statmech_linker_activeintpol)
  {
    // reverts reference lengths of elements back to the way they were before... (has to be done prior to restoring actlinklength_. Otherwise, fatal loss of information)
    ChangeActiveLinkerLength(1e9, 1e9, crosslinkeractlength_, crosslinkeractlengthconv_,true);
    crosslinkeractlength_ = Teuchos::rcp(new Epetra_Vector(*crosslinkeractlengthconv_));
  }

  /*restore state of the discretization at the beginning of this time step; note that to this and
   *adding and deleting crosslinker element has to be undone exactly vice-versa compared to the way
   *it was done first in order to handle also those crosslinkers correctly added and deleted in one
   *and the same time step*/

  //loop through all elements deleted in this time step and restore them in the discretization
  for(int i=0; i<(int)deletedelements_.size(); i++)
  {
    std::vector<char> tmp;
    std::vector<char>::size_type position = 0;
    DRT::ParObject::ExtractfromPack(position,deletedelements_[i],tmp);
    DRT::ParObject* o = DRT::UTILS::Factory(tmp);
    DRT::Element* ele = dynamic_cast<DRT::Element*>(o);
    if (ele == NULL)
      dserror("Failed to build an element from the element data");
    discret_->AddElement(Teuchos::rcp(ele));
  }
  deletedelements_.clear();

  //loop through addedelements_, delete all these elements and then set addedelements_ an empty vector
  for(int i=0; i<(int)addedelements_.size(); i++)
    discret_->DeleteElement(addedelements_[i]);
  addedelements_.clear();

  /*settling administrative stuff in order to make the discretization ready for the next time step: synchronize
   *the Filled() state on all processors after having added or deleted elements by ChekcFilledGlobally(); then build
   *new element maps and call FillComplete(); finally Crs matrices stiff_ has to be deleted completely and made ready
   *for new assembly since their graph was changed*/
  discret_->CheckFilledGlobally();
  discret_->FillComplete(true, false, false);
  stiff->Reset();

  // same procedure for contact discretization
  if(beamcmanager!=Teuchos::null)
  {
    for(int i=0; i<(int)deletedcelements_.size(); i++)
    {
      std::vector<char> tmp;
      std::vector<char>::size_type position = 0;
      DRT::ParObject::ExtractfromPack(position,deletedcelements_[i],tmp);
      DRT::ParObject* o = DRT::UTILS::Factory(tmp);
      DRT::Element* ele = dynamic_cast<DRT::Element*>(o);
      if (ele == NULL)
        dserror("Failed to build an element from the element data");
      beamcmanager->ContactDiscret().AddElement(Teuchos::rcp(ele));
    }
    deletedcelements_.clear();

    for(int i=0; i<(int)addedcelements_.size(); i++)
      beamcmanager->ContactDiscret().DeleteElement(addedcelements_[i]);
    addedcelements_.clear();

    // contact discretization
    beamcmanager->ContactDiscret().CheckFilledGlobally();
    beamcmanager->ContactDiscret().FillComplete();
  }

  return;
} // StatMechManager::RestoreConv()

/*----------------------------------------------------------------------*
 | check for broken element                        (public) mueller 3/10|
 *----------------------------------------------------------------------*/
bool STATMECH::StatMechManager::CheckForBrokenElement(LINALG::SerialDenseMatrix& coord, LINALG::SerialDenseMatrix& cut)
{
  // empty cut just in case it was handed over non-empty
  cut.Zero();
  bool broken = false;
  int ndim = coord.M();
  // flag "broken" signals a broken element, flag "cut" hints at location of nodes 0 and 1 (of the two-noded beam)
  for (int dof = 0; dof < ndim; dof++)
    //loop through columns of "cut"
    for (int n = 0; n < cut.N(); n++)
    {
      // broken element with node_n close to "0.0"-boundary
      if (fabs(coord(dof, n + 1) - periodlength_->at(dof) - coord(dof, n)) < fabs(coord(dof, n + 1) - coord(dof, n)))
      {
        broken = true;
        // set value for the spatial component in question at n-th cut
        cut(dof, n) = 1.0;
      }
      else if (fabs(coord(dof, n + 1) + periodlength_->at(dof) - coord(dof, n)) < fabs(coord(dof, n + 1) - coord(dof, n)))
      {
        broken = true;
        cut(dof, n) = 2.0;
      }
    }
  return broken;
}// StatMechManager::CheckForBrokenElement

/*----------------------------------------------------------------------*
 | get a matrix with node coordinates and their DOF LIDs                |
 |                                                 (public) mueller 3/10|
 *----------------------------------------------------------------------*/
void STATMECH::StatMechManager::GetElementNodeCoords(DRT::Element* element, Teuchos::RCP<Epetra_Vector> discol, LINALG::SerialDenseMatrix& coord, std::vector<int>* lids)
{
  // clear LID vector just in case it was handed over non-empty
  lids->clear();
  for (int j=0; j<element->NumNode(); j++)
  {
    int nodegid = element->NodeIds()[j];
    DRT::Node* node = discret_->lColNode(discret_->NodeColMap()->LID(nodegid));
    for (int k=0; k<3; k++)
    {
      // obtain k-th spatial component of the reference position of the j-th node
      double referenceposition = node->X()[k];
      // get the GIDs of the node's DOFs
      std::vector<int> dofnode = discret_->Dof(node);
      // store the displacement of the k-th spatial component
      double displacement = (*discol)[discret_->DofColMap()->LID(dofnode[k])];
      // write updated components into coord (only translational
      coord(k, j) = referenceposition + displacement;
      // store current lid(s) (3 translational DOFs per node)
      if (lids != NULL)
        lids->push_back(discret_->DofRowMap()->LID(dofnode[k]));
    }
  }
  return;
} // StatMechManager::GetElementNodeCoords

/*----------------------------------------------------------------------*
 | update force sensor locations                   (public) mueller 3/10|
 *----------------------------------------------------------------------*/
void STATMECH::StatMechManager::UpdateForceSensors(std::vector<int>& sensornodes, int dbcdispdir)
{
  // reinitialize forcesensor_
  forcesensor_->PutScalar(-1.0);

  // loop over DOFs subjected to oscillation (by DBC)
  for (int i=0; i<(int)sensornodes.size(); i++)
  {
    // check if node is row node (this ensures processor-wise unique forcesensor_ vectors)
    if(discret_->NodeRowMap()->LID(sensornodes[i])>-1)
    {
      // get the node
      DRT::Node* actnode = discret_->gNode(sensornodes.at(i));
      // get the GID of the DOF of the oscillatory motion
      int dofgid = discret_->Dof(0, actnode)[dbcdispdir];
      // now, get the LID
      int lid = discret_->DofColMap()->LID(dofgid);
      // activate force sensor at lid-th position
      (*forcesensor_)[lid] = 1.0;
    }
  }
} // StatMechManager::UpdateForceSensors

/*----------------------------------------------------------------------*
 | update number of unconverged steps              (public) mueller 3/10|
 *----------------------------------------------------------------------*/
void STATMECH::StatMechManager::UpdateNumberOfUnconvergedSteps(bool flag)
{
  if(flag)
    unconvergedsteps_++;
  else
    unconvergedsteps_--;
  return;
}

/*----------------------------------------------------------------------*
 | add parameters to given parameter list          (public) mueller 3/10|
 *----------------------------------------------------------------------*/
void STATMECH::StatMechManager::AddStatMechParamsTo(Teuchos::ParameterList& params, Teuchos::RCP<Epetra_MultiVector> randomnumbers)
{
  params.set("ETA",statmechparams_.get<double>("ETA",0.0));
  params.set("THERMALBATH",DRT::INPUT::IntegralValue<INPAR::STATMECH::ThermalBathType>(statmechparams_,"THERMALBATH"));
  params.set<int>("FRICTION_MODEL",DRT::INPUT::IntegralValue<INPAR::STATMECH::FrictionModel>(statmechparams_,"FRICTION_MODEL"));
  params.set<INPAR::STATMECH::DBCType>("DBCTYPE", DRT::INPUT::IntegralValue<INPAR::STATMECH::DBCType>(statmechparams_,"DBCTYPE"));
  if(randomnumbers!=Teuchos::null)
    params.set("RandomNumbers",randomnumbers);
  params.set("SHEARAMPLITUDE",statmechparams_.get<double>("SHEARAMPLITUDE",0.0));
  // note CURVENUMBER and DBCDISPDIR follow the counting convention 1,2,3,... (not C++). They are internally decremented
  params.set("CURVENUMBER",statmechparams_.get<int>("CURVENUMBER",-1));
  params.set("DBCDISPDIR",statmechparams_.get<int>("DBCDISPDIR",-1));
  params.set("STARTTIMEACT",actiontime_->at(dbctimeindex_));
  if(DRT::INPUT::IntegralValue<INPAR::STATMECH::DBCType>(statmechparams_,"DBCTYPE")==INPAR::STATMECH::dbctype_affineshear)
  params.set("STARTTIMEACT",actiontime_->at(dbctimeindex_));
  params.set("DELTA_T_NEW",timestepsizes_->at(dbctimeindex_));
  params.set("PERIODLENGTH",GetPeriodLength());
  if(linkermodel_ == statmech_linker_bellseq || linkermodel_ == statmech_linker_bellseqintpol ||
     linkermodel_ == statmech_linker_active || linkermodel_ == statmech_linker_activeintpol ||
     networktype_ == statmech_network_casimir)
    params.set<std::string>("internalforces","yes");
  return;
}

/*----------------------------------------------------------------------*
 | generates a vector with a random permutation of the numbers between 0|
 | and N - 1                                        (public) cyron 06/10|
 *----------------------------------------------------------------------*/
std::vector<int> STATMECH::StatMechManager::Permutation(const int& N)
{
  //auxiliary variable
  int j = 0;

  //result vector initialized with ordered numbers from 0 to N-1
  std::vector<int> result(N, 0);
  for (int i=0; i<(int)result.size(); i++)
    result[i] = i;

  for (int i=0; i<N; ++i)
  {
    //generate random number between 0 and i
    j = (int)floor((i + 1.0)*(*uniformgen_)());

    /*exchange values at positions i and j (note: value at position i is i due to above initialization
     *and because so far only positions <=i have been changed*/
    result[i] = result[j];
    result[j] = i;
  }

  return result;
} // StatMechManager::Permutation

/*----------------------------------------------------------------------*
 | Computes current internal energy of discret_ (public)     cyron 12/10|
 *----------------------------------------------------------------------*/
void STATMECH::StatMechManager::ComputeInternalEnergy(const Teuchos::RCP<Epetra_Vector> dis,
                                                      std::vector<double>&              energy,
                                                      const double&                     dt,
                                                      const std::ostringstream&         filename,
                                                      bool                              fillzeros,
                                                      bool                              writefile)
{
  ParameterList p;
  p.set("action", "calc_struct_energy");

  //add statistical vector to parameter list for statistical forces and damping matrix computation
  p.set("delta time",dt);
  p.set("ETA",statmechparams_.get<double>("ETA",0.0));
  p.set("THERMALBATH",DRT::INPUT::IntegralValue<INPAR::STATMECH::ThermalBathType>(statmechparams_,"THERMALBATH"));
  p.set<int>("FRICTION_MODEL",DRT::INPUT::IntegralValue<INPAR::STATMECH::FrictionModel>(statmechparams_,"FRICTION_MODEL"));
  p.set("SHEARAMPLITUDE",statmechparams_.get<double>("SHEARAMPLITUDE",0.0));
  p.set("CURVENUMBER",statmechparams_.get<int>("CURVENUMBER",-1));
  p.set("DBCDISPDIR",statmechparams_.get<int>("DBCDISPDIR",-1));
  p.set("PERIODLENGTH", periodlength_);

  discret_->ClearState();
  discret_->SetState("displacement", dis);
  Teuchos::RCP<Epetra_SerialDenseVector> energies = Teuchos::rcp(new Epetra_SerialDenseVector(1));
  energies->Scale(0.0);
  //filaments
  p.set("energyoftype", "beam3ii");
  discret_->EvaluateScalars(p, energies);
  energy.push_back((*energies)(0));
  //crosslinkers
  energies->Scale(0.0);
  p.remove("energyoftype", true);
  if(linkermodel_ == statmech_linker_stdintpol  || linkermodel_ == statmech_linker_activeintpol || linkermodel_ == statmech_linker_bellseqintpol)
    p.set("energyoftype", "beam3cl");
  else
    p.set("energyoftype", "beam3");
  discret_->EvaluateScalars(p, energies);
  energy.push_back((*energies)(0));
  discret_->ClearState();

  if(!discret_->Comm().MyPID() && writefile)
  {
    FILE* fp = NULL;
    fp = fopen(filename.str().c_str(), "a");
    std::stringstream writetofile;
    for(int i=0; i<(int)energy.size(); i++)
      writetofile<<energy[i]<<"    ";
    if(fillzeros)
      for(int i=0; i<17-(int)energy.size(); i++)
        writetofile<<"    "<<0;
    writetofile<<std::endl;
    fprintf(fp, writetofile.str().c_str());
    fclose(fp);
  }
  return;
} // StatMechManager::ComputeInternalEnergy

/*----------------------------------------------------------------------*
 | checks orientation of crosslinker relative to linked filaments       |
 |                                                  (public) cyron 06/10|
 *----------------------------------------------------------------------*/

bool STATMECH::StatMechManager::CheckOrientation(const LINALG::Matrix<3, 1> direction, const Epetra_MultiVector& nodaltriadscol, const Epetra_SerialDenseMatrix& LID, Teuchos::RCP<double> phifil)
{

  //if orientation is not to be checked explicitly, this function always returns true
  if (!DRT::INPUT::IntegralValue<int>(statmechparams_, "CHECKORIENT"))
    return true;

  if(statmechparams_.get<double>("KT") == 0.0)
    dserror("Set KT to a non-zero value.");

  if(LID.M()!=2 && LID.N()!=1)
    dserror("LID has wrong dimensions %d x %d", LID.M(),LID.M());

  //triads on filaments at the two nodes connected by crosslinkers
  LINALG::Matrix<3, 3> T1;
  LINALG::Matrix<3, 3> T2;

  //auxiliary variable for storing a triad in quaternion form
  LINALG::Matrix<4, 1> qnode;

  //angle between filament axes at crosslinked points, respectively
  double Phi;

  //Deltaphi = Phi - Phi0, where Phi0 is the angle between crosslinked filaments with zero potential energy (i.e. the most likely one)
  double DeltaPhi;

  //triad of node on first filament which is affected by the new crosslinker
  for (int j=0; j<4; j++)
    qnode(j) = nodaltriadscol[j][(int) LID(0,0)];
  LARGEROTATIONS::quaterniontotriad(qnode, T1);

  //triad of node on second filament which is affected by the new crosslinker
  for (int j=0; j<4; j++)
    qnode(j) = nodaltriadscol[j][(int) LID(1,0)];
  LARGEROTATIONS::quaterniontotriad(qnode, T2);

  //auxiliary variable
  double scalarproduct = T1(0, 0)*T2(0,0) + T1(1,0)*T2(1,0) + T1(2,0)*T2(2,0);

  Phi = acos(scalarproduct);
  //Phi should be the acute angle between 0° and 90° between the filament axes
  if(Phi > M_PI/2.0)
    Phi = M_PI - Phi;

  if(phifil!=Teuchos::null)
    *phifil = Phi;

  DeltaPhi = Phi - statmechparams_.get<double> ("PHI0",0.0);

  //assuming bending and torsion potentials 0.5*EI*Deltaphi^2 and a Boltzmann distribution for the different states of the crosslinker we get
  double pPhi = exp(-0.5 * statmechparams_.get<double> ("CORIENT", 0.0) * DeltaPhi * DeltaPhi / statmechparams_.get<double> ("KT", 0.0));


  //pPhi = 0.0 if DeltaPhi is outside allowed range
  if(Phi < statmechparams_.get<double>("PHI0",0.0) - 0.5*statmechparams_.get<double>("PHI0DEV",6.28) ||
     Phi > statmechparams_.get<double>("PHI0",0.0) + 0.5*statmechparams_.get<double>("PHI0DEV",6.28))
     pPhi = 0.0;

  //crosslinker has to pass three probability checks with respect to orientation
  return((*uniformgen_)() < pPhi);
} // StatMechManager::CheckOrientation

/*----------------------------------------------------------------------*
| simulation of crosslinker diffusion           (public) mueller 07/10|
*----------------------------------------------------------------------*/
void STATMECH::StatMechManager::CrosslinkerDiffusion(const Epetra_MultiVector& bspotpositions,
                                                     double                    mean,
                                                     double                    standarddev,
                                                     const double&             dt)
{
  /* Here, the diffusion of crosslink molecules is handled.
   * Depending on the number of occupied binding spots of the molecule, its motion
   * is calculated differently.
   */
#ifdef DEBUGCOUT
  Teuchos::RCP<Epetra_MultiVector> crosslinkerdeltatrans = Teuchos::rcp(new Epetra_MultiVector(*transfermap_, 3, true));
#endif

  Teuchos::RCP<Epetra_MultiVector> crosslinkerpositionstrans = Teuchos::rcp(new Epetra_MultiVector(*transfermap_, 3, true));
  CommunicateMultiVector(crosslinkerpositionstrans, crosslinkerpositions_, true, false,true);

  // bonding cases
  for (int i=0; i<crosslinkerpositionstrans->MyLength(); i++)
  {
    int crosslid = crosslinkermap_->LID(transfermap_->GID(i));
    switch ((int)(*numbond_)[crosslid])
    {
      // bonding case 1:  no bonds, diffusion
      case 0:
        for (int j=0; j<crosslinkerpositionstrans->NumVectors(); j++)
        {
          (*crosslinkerpositionstrans)[j][i] += standarddev*(*normalgen_)() + mean;
#ifdef DEBUGCOUT
          (*crosslinkerdeltatrans)[j][i] = standarddev*(*normalgen_)() + mean;
#endif
        }
      break;
      // bonding case 2: crosslink molecule attached to one filament
      case 1:
      {
        int bspotLID = bspotcolmap_->LID(std::max((int)(*crosslinkerbond_)[0][crosslid],(int)(*crosslinkerbond_)[1][crosslid]));
        for (int j=0; j<crosslinkerpositionstrans->NumVectors(); j++)
          (*crosslinkerpositionstrans)[j][i] = bspotpositions[j][bspotLID];
      }
      break;
      // bonding case 3: actual crosslinker has been established or passified molecule
      case 2:
      {
        int largerbspotLID = bspotcolmap_->LID(std::max((int)(*crosslinkerbond_)[0][crosslid],(int)(*crosslinkerbond_)[1][crosslid]));
        int smallerbspotLID = bspotcolmap_->LID(std::min((int)(*crosslinkerbond_)[0][crosslid],(int)(*crosslinkerbond_)[1][crosslid]));
        if(smallerbspotLID>-1)
          for (int j=0; j<crosslinkerpositionstrans->NumVectors(); j++)
            (*crosslinkerpositionstrans)[j][i] = (bspotpositions[j][largerbspotLID]+bspotpositions[j][smallerbspotLID])/2.0;
        else
          for (int j=0; j<crosslinkerpositionstrans->NumVectors(); j++)
            (*crosslinkerpositionstrans)[j][i] = bspotpositions[j][largerbspotLID];
      }
      break;
    }
  }
  // check for compliance with periodic boundary conditions if existent
  if (periodlength_->at(0) > 0.0)
    CrosslinkerPeriodicBoundaryShift(crosslinkerpositionstrans);

  // Update by Broadcast: make this information redundant on all procs
  CommunicateMultiVector(crosslinkerpositionstrans, crosslinkerpositions_, false, true);

#ifdef DEBUGCOUT
  Teuchos::RCP<Epetra_MultiVector> crosslinkerdelta = Teuchos::rcp(new Epetra_MultiVector(*crosslinkermap_, 3, true));
  CommunicateMultiVector(crosslinkerdeltatrans, crosslinkerdelta, false,true);

  for(int i=0; i<crosslinkerpositionstrans->NumVectors(); i++)
    for(int j=0; j<crosslinkerpositionstrans->MyLength(); j++)
      if((*crosslinkerpositionstrans)[i][j]>periodlength_->at(i) ||(*crosslinkerpositionstrans)[i][j]<0.0)
        dserror("Crosslinker %i outside of the periodic box: %d, %d, %d ", j, (*crosslinkerpositionstrans)[0][j], (*crosslinkerpositionstrans)[1][j], (*crosslinkerpositionstrans)[2][j]);

  if(!discret_->Comm().MyPID())
  {
    // check increments of the stochastic process
    std::stringstream crossposfile;
    crossposfile<<"./crosspos.txt";
    FILE* fpcross = NULL;
    fpcross = fopen(crossposfile.str().c_str(), "w");

    std::stringstream crosspos;
    for(int i=0; i<crosslinkerdelta->NumVectors(); i++)
      for(int j=0; j<crosslinkerdelta->MyLength(); j++)
        crosspos<<(*crosslinkerdelta)[i][j]<<" "<<(*crosslinkerdelta)[i][j]<<" "<<(*crosslinkerdelta)[i][j]<<std::endl;
    // move temporary std::stringstream to file and close it
    fprintf(fpcross, crosspos.str().c_str());
    fclose(fpcross);

    // check paths of the stochastic process
    if(crosslinkerpositions_->MyLength()>1)
      dserror("Check only for a single crosslink molecule!");
    std::ostringstream filename;
    filename << outputrootpath_ << "/crosslinkermotion.dat";
    FILE *fp = fopen(filename.str().c_str(), "a");

    std::stringstream linkerposition;
    linkerposition<<(*crosslinkerpositions_)[0][0]<<" "<<(*crosslinkerpositions_)[1][0]<<" "<<(*crosslinkerpositions_)[2][0]<<std::endl;
    fprintf(fp, linkerposition.str().c_str());
    fclose(fp);
  }
#endif
  return;
}// StatMechManager::CrosslinkerDiffusion

/*----------------------------------------------------------------------*
 | update crosslink molecule positions                                  |
 |                                                (public) mueller 08/10|
 *----------------------------------------------------------------------*/
void STATMECH::StatMechManager::CrosslinkerIntermediateUpdate(const Epetra_MultiVector& bspotpositions,
                                                              const LINALG::SerialDenseMatrix& LID,
                                                              const int& crosslinkernumber,
                                                              bool coupledmovement)
{
  // case: one-bonded crosslink molecule (i.e. two cases: +1 bond (starting from 0 bonds) or -1 bond (molecule is free to diffuse again)
  if (LID.M()==1 && LID.N()==1)
  {
    // set molecule position to node position
    if (coupledmovement)
    {
      for (int i=0; i < crosslinkerpositions_->NumVectors(); i++)
        (*crosslinkerpositions_)[i][crosslinkernumber] = bspotpositions[i][(int)LID(0,0)];
    }
    else
    {
      // generate vector in random direction of length R_LINK to "reset" crosslink molecule position:
      // it may now reenter or leave the bonding proximity
      LINALG::Matrix<3, 1> deltapos;
      for (int i=0; i<(int)deltapos.M(); i++)
        deltapos(i) = (*uniformgen_)();
      deltapos.Scale(statmechparams_.get<double> ("R_LINK", 0.0) / deltapos.Norm2());
      for (int i=0; i<crosslinkerpositions_->NumVectors(); i++)
        (*crosslinkerpositions_)[i][crosslinkernumber] += deltapos(i);
    }
  }
  // case: crosslinker element
  if (LID.M()==2 && LID.N()==1)
  {
    int bspotlid = std::max((int)LID(0,0), (int)LID(1,0));
    // in case of a loom network, we want the crosslinker position to lie on the horizontal filament
    if(networktype_ == statmech_network_casimir)
      for(int i=0; i<LID.M(); i++)
        if((*filamentnumber_)[(int)LID(i,0)] == 0)
        {
          bspotlid = (int)LID(i,0);
          break;
        }

    for (int i=0; i<crosslinkerpositions_->NumVectors(); i++)
      (*crosslinkerpositions_)[i][crosslinkernumber] = bspotpositions[i][bspotcolmap_->GID(bspotlid)];
  }
  return;
}// StatMechManager::CrosslinkerIntermediateUpdate

/*----------------------------------------------------------------------*
 | Initialize crosslinker positions               (public) mueller 07/10|
 *----------------------------------------------------------------------*/
void STATMECH::StatMechManager::CrosslinkerMoleculeInit()
{
  int ncrosslink = statmechparams_.get<int> ("N_crosslink", 0);
  int numbins = statmechparams_.get<int>("HISTOGRAMBINS", 1);

  // create crosslinker maps
  std::vector<int> gids;
  for (int i=0; i<ncrosslink; i++)
    gids.push_back(i);
  // crosslinker column and row map
  crosslinkermap_ = Teuchos::rcp(new Epetra_Map(-1, ncrosslink, &gids[0], 0, discret_->Comm()));
  transfermap_    = Teuchos::rcp(new Epetra_Map(ncrosslink, 0, discret_->Comm()));
  startindex_ = Teuchos::rcp(new std::vector<double>);

  // create maps for binding spots in case of helical binding spot geometry of the actin filament
  if(filamentmodel_ == statmech_filament_helical)
  {
    // rise and rotation per binding spot (default values: 2.77nm and -166.15deg (left-handed, one-start helix))
    double riseperbspot = statmechparams_.get<double>("RISEPERBSPOT",0.00277);
    double rotperbspot = statmechparams_.get<double>("ROTPERBSPOT", -2.8999);

    //getting a vector consisting of pointers to all filament number conditions set
    std::vector<DRT::Condition*> filaments(0);
    discret_->GetCondition("FilamentNumber",filaments);

    // apply new beam element with intermediate binding spot positions
    if(linkermodel_ == statmech_linker_stdintpol  || linkermodel_ == statmech_linker_activeintpol || linkermodel_ == statmech_linker_bellseqintpol)
    {
      // temporary vectors
      // binding spot orientations relative to second triad vector
      std::vector<double> bspotorientations;
      // value of binding spot in parameter space
      std::vector<double> bspotxi;
      // maps binding spot to element
      std::vector<int>    bspot2element;
      // vector holding all bspot gids
      std::vector<int>    bspotgids;
      // vector signaling the existence of a binding spot on this proc
      std::vector<int>     bspotonproc;
      // CHECK vector
      std::vector<int>   filnodes;


      //loop over filaments to create redundant column map
      for (int i=0; i<(int)filaments.size(); i++)
      {
        //get a pointer to nodal cloud covered by the current condition
        const std::vector<int>* nodeids = filaments[i]->Nodes();
        // retrieve filament length using first and last node of the current filament
        DRT::Node* node0 = discret_->lColNode(discret_->NodeColMap()->LID((*nodeids)[0]));
        double lfil = 0.0;
        std::vector<double> pos(6);
        // loop through all Nodes of the filaments, consider breakage by Periodic Boundary and calculate length of filament i
        for(int j=1;j<(int)nodeids->size();j++)
        {
          DRT::Node* tmpnode = discret_->lColNode(discret_->NodeColMap()->LID((*nodeids)[j]));
          for(int l=0;l<3;l++)
            pos[l+3]= tmpnode->X()[l];

          if(j==1)
          { for(int l=0;l<3;l++)
             pos[l]= node0->X()[l];
          }

          //shift position in case of breakage due to Periodic Boundary
          for(int k=0;k<3;k++)
          {
            if( fabs( pos[k+3] + periodlength_->at(k) - pos[k] ) < fabs( pos[k+3] - pos[k] ) )
            {
              pos[k+3] += periodlength_->at(k);
            }
            if( fabs( pos[k+3] - periodlength_->at(k) - pos[k] ) < fabs( pos[k+3] - pos[k] ) )
            {
              pos[k+3] -= periodlength_->at(k);
            }
          }

          lfil += sqrt((pos[3]-pos[0])*(pos[3]-pos[0])+
                       (pos[4]-pos[1])*(pos[4]-pos[1])+
                       (pos[5]-pos[2])*(pos[5]-pos[2]));

          for(int l=0;l<3;l++)
            pos[l] = pos[l+3];
        }
        // element length (assuming equal node spacing)
        double elelength = lfil/(double)(filaments[i]->Nodes()->size()-1);

        // add as many binding spots to the i-th filament as possible
        int bspot = 0;
        while((double)bspot*riseperbspot <= lfil)
        {
          // put binding spot id into the vector in order to create a fully overlapping binding spot map afterwards
          bspotgids.push_back((int)bspotgids.size());
          //std::cout<< " BSPOTGIDS "<< (int)bspotgids << std::endl;


          // calculate the orientation of the binding spot as well as the respective curve parameter
          /* determine element current element GID:
           * As set up by the Inputfile-Generator, the filament element GID can be retrieved by knowing the filament number and the lesser node number.
           * Attention: If that property is changed, one should of course refer to node->Elements(). In the meanwhile, this method seems faster*/
          //int elegid = node0->Id() - i + (int)(floor((double)bspot*riseperbspot/lfil));
          int ival;
          // we need to make sure that bspot belongs to existing Filaments. In case a bspot is located on a node we always choose Element with the lower GID
          ival = (int)(floor((double)bspot*riseperbspot/elelength));
          // There might be a more elegant way
          if((int)(ceil((double)bspot*riseperbspot/elelength))==ival && bspot != 0)
            ival-=1;
          int elegid = node0->Id() - i + ival;

          // add element GID if this element is a row map element on this processor
          if(discret_->ElementRowMap()->LID(elegid)>-1)
          {
            // mark binding spot as being on this proc
            bspotonproc.push_back(1);
            bspot2element.push_back(elegid);

            // orientation [0;2pi] and parameter xi_bs of the binding spot
            double phibs = 0.0;
            double xibs = -1;
            if(bspot==0)
            {
              bspotorientations.push_back(phibs);
              bspotxi.push_back(xibs);
            }
            else
            {
              phibs = (((double)(bspot)*rotperbspot)/(2.0*M_PI)-floor(((double)(bspot-1)*rotperbspot)/(2.0*M_PI)))*2.0*M_PI;
              ival = (int)floor((double)bspot*riseperbspot/elelength);
              if((int)ceil((double)bspot*riseperbspot/elelength)==ival)
                  ival-=1;
              xibs = 2*((double)bspot*riseperbspot/elelength - ival - 0.5);
              bspotorientations.push_back(phibs);
              bspotxi.push_back(xibs);
            }

          }
          else // only bspotonproc stores the information about off-proc elements since the structure can be easily rebuilt for the other vectors
            bspotonproc.push_back(0);
          // next binding spot on filament
          bspot++;
        }
      }

      // maps
      // create redundant binding spot column map based upon the bspotids
      bspotcolmap_ = Teuchos::rcp(new Epetra_Map(-1, (int)bspotgids.size(), &bspotgids[0], 0, discret_->Comm()));
      // create processor-specific row maps
      std::vector<int> bspotrowgids;
      for(int i=0; i<(int)bspotonproc.size(); i++)
        if(bspotonproc[i]==1)
          bspotrowgids.push_back(i);  // note: since column map is fully overlapping: i=col. LID = GID
      bspotrowmap_ = Teuchos::rcp(new Epetra_Map((int)bspotgids.size(), (int)bspotrowgids.size(), &bspotrowgids[0], 0, discret_->Comm()));

      // vectors
      // initialize class vectors
      bspotstatus_ = Teuchos::rcp(new Epetra_Vector(*bspotcolmap_));
      bspotstatus_->PutScalar(-1.0);
      bspotorientations_ = rcp(new Epetra_Vector(*bspotcolmap_,true));
      bspotxi_ = rcp(new Epetra_Vector(*bspotcolmap_,true));
      bspot2element_ = rcp(new Epetra_Vector(*bspotrowmap_, true));
      bspot2nodes_=rcp(new Epetra_MultiVector(*bspotcolmap_,2,true));

      Epetra_Vector bspotorientationsrow(*bspotrowmap_);
      Epetra_Vector bspotxirow(*bspotrowmap_);
      Epetra_MultiVector bspot2nodesrow(*bspotrowmap_,2);
      for(int i=0; i<bspotrowmap_->NumMyElements(); i++)
      {
        bspotorientationsrow[i] = bspotorientations[i];
        bspotxirow[i] = bspotxi[i];
        (*bspot2element_)[i] = bspot2element[i];
        DRT::Element *ele = discret_->gElement(bspot2element[i]);
        for(int l=0;l<ele->NumNode();l++)
          bspot2nodesrow[l][i]=ele->NodeIds()[l];
      }
      // make information on orientations and curve parameters redundant on all procs
      Epetra_Import bspotimporter(*bspotcolmap_, *bspotrowmap_);
      bspotorientations_->Import(bspotorientationsrow,bspotimporter,Insert);
      bspotxi_->Import(bspotxirow,bspotimporter,Insert);
      bspot2nodes_->Import(bspot2nodesrow,bspotimporter,Insert);
    }
    else // beam3/beam3ii element: new binding spot maps are equivalent to node maps
    {
      // maps
      bspotcolmap_ = Teuchos::rcp(new Epetra_Map(*(discret_->NodeColMap())));
      bspotrowmap_ = Teuchos::rcp(new Epetra_Map(*(discret_->NodeRowMap())));
      // initialize class vectors (we do not need nspotxi_ here since the nodes are the binding spots)
      bspotstatus_ = Teuchos::rcp(new Epetra_Vector(*bspotcolmap_));
      bspotstatus_->PutScalar(-1.0);
      bspotorientations_ = Teuchos::rcp(new Epetra_Vector(*bspotcolmap_));
      bspot2element_ = Teuchos::rcp(new Epetra_Vector(*bspotrowmap_));

      // orientations and vectors
      int bspotgid = 0;
      for(int i=0; i<(int)filaments.size(); i++)
      {
        //get a pointer to nodal cloud covered by the current condition
        const std::vector<int>* nodeids = filaments[i]->Nodes();
        // retrieve filament length using first and last node of the current filament
        DRT::Node* node0 = discret_->lColNode(discret_->NodeColMap()->LID((*nodeids)[0]));
        //DRT::Node* node1 = discret_->lColNode(discret_->NodeColMap()->LID((*nodeids)[((int)nodeids->size())-1]));
       double lfil= 0.0;
       std::vector<double> pos(6);
        // loop through all Nodes of the filaments, consider breakage by Periodic Boundary and calculate length of filament i
        for(int j=1;j<(int)nodeids->size();j++)
        {
          DRT::Node* tmpnode = discret_->lColNode(discret_->NodeColMap()->LID((*nodeids)[j]));
          for(int l=0;l<3;l++)
            pos[l+3]= tmpnode->X()[l];

          if(j==1)
          {
            for(int l=0;l<3;l++)
             pos[l]= node0->X()[l];
          }
          //shift position in case of breakage due to Periodic Boundary
          for(int k=0;k<3;k++)
          { 
          	if( fabs( pos[k+3] + periodlength_->at(k) - pos[k] ) < fabs( pos[k+3] - pos[k] ) )
            {
              pos[k+3] += periodlength_->at(k);
            }
           	if( fabs( pos[k+3] - periodlength_->at(k) - pos[k] ) < fabs( pos[k+3] - pos[k] ) )
            {
              pos[k+3] -= periodlength_->at(k);
            }
          }

          lfil += sqrt(pow(pos[3]-pos[0],2)+pow(pos[4]-pos[1],2)+pow(pos[5]-pos[2],2));

          for(int l=0;l<3;l++)
           pos[l] = pos[l+3];
        }

        // element length (assuming equal node spacing)
        double elelength = lfil/(double)(filaments[i]->Nodes()->size()-1);

        for(int j=0; j<(int)filaments[i]->Nodes()->size(); j++)
        {
          double phibs = 0.0;
          // first node
          if(j==0)
            (*bspotorientations_)[bspotcolmap_->LID(bspotgid)] = phibs;
          else
          {
            // determine orientation (relative to second triad vector)
            phibs = (((double)(j)*(elelength/riseperbspot)*rotperbspot)/(2.0*M_PI)-floor(((double)(j-1)*(elelength/riseperbspot)*rotperbspot)/(2.0*M_PI)))*2.0*M_PI;
            (*bspotorientations_)[bspotcolmap_->LID(bspotgid)] = phibs;
          }
          // assumption:
          if(bspotrowmap_->LID(bspotgid)>-1)
            (*bspot2element_)[bspotrowmap_->LID(bspotgid)] = discret_->lRowNode(discret_->NodeRowMap()->LID(bspotgid))->Elements()[0]->Id();
          bspotgid++;
        }
      }
    }
  }
  else // without helical orientation
  {
    if(linkermodel_ == statmech_linker_stdintpol || linkermodel_ == statmech_linker_activeintpol || linkermodel_ == statmech_linker_bellseqintpol)
    {
      // rise per binding spot (default values: 2.77nm)
      double riseperbspot = statmechparams_.get<double>("RISEPERBSPOT",0.00277);

      //getting a vector consisting of pointers to all filament number conditions set
      std::vector<DRT::Condition*> filaments(0);
      discret_->GetCondition("FilamentNumber",filaments);

      // value of binding spot in parameter space
      std::vector<double> bspotxi;
      // maps binding spot to element
      std::vector<int>    bspot2element;
      // vector holding all bspot gids
      std::vector<int>    bspotgids;
      // vector signaling the existence of a binding spot on this proc
      std::vector<int>    bspotonproc;
      // CHECK vector
      std::vector<int>   filnodes;

      //loop over filaments to create redundant column map
      for (int i=0; i<(int)filaments.size(); i++)
      {
        // filament length
        double lfil = 0;
        //get a pointer to nodal cloud covered by the current condition
        const std::vector<int>* nodeids = filaments[i]->Nodes();
        // retrieve filament length using first and last node of the current filament
        DRT::Node* node0 = discret_->lColNode(discret_->NodeColMap()->LID((*nodeids)[0]));
        //DRT::Node* node1 = discret_->lColNode(discret_->NodeColMap()->LID((*nodeids)[((int)nodeids->size())-1]));

        std::vector<double> pos(6);
        // loop through all Nodes of the filaments, consider breakage by Periodic Boundary and calculate length of filament i
        for(int j=1;j<(int)nodeids->size();j++)
        {
          DRT::Node* tmpnode = discret_->lColNode(discret_->NodeColMap()->LID((*nodeids)[j]));
          for(int l=0;l<3;l++)
            pos[l+3]= tmpnode->X()[l];

          if(j==1)
          { for(int l=0;l<3;l++)
             pos[l]= node0->X()[l];
          }
          //shift position in case of breakage due to Periodic Boundary
          for(int k=0;k<3;k++)
          { 
          	if( fabs( pos[k+3] + periodlength_->at(k) - pos[k] ) < fabs( pos[k+3] - pos[k] ) )
            {
              pos[k+3] += periodlength_->at(k);
            }
           	if( fabs( pos[k+3] - periodlength_->at(k) - pos[k] ) < fabs( pos[k+3] - pos[k] ) )
            {
              pos[k+3] -= periodlength_->at(k);
            }
          }

          lfil += sqrt((pos[3]-pos[0])*(pos[3]-pos[0])+
                       (pos[4]-pos[1])*(pos[4]-pos[1])+
                      (pos[5]-pos[2])*(pos[5]-pos[2]));

          for(int l=0;l<3;l++)
          	pos[l] = pos[l+3];
        }
        // element length (assuming equal node spacing)
        double elelength = lfil/(double)(filaments[i]->Nodes()->size()-1);

        // add as many binding spots to the i-th filament as possible
        int bspot = 0;
        //more elegant way for tolerance?
        double ltol = 1e-7;
        while((double)bspot*riseperbspot < lfil+ltol)
        {
          // put binding spot id into the vector in order to create a fully overlapping binding spot map afterwards
          bspotgids.push_back((int)bspotgids.size());
          // calculate the orientation of the binding spot as well as the respective curve parameter
          /* determine element current element GID:
           * As set up by the Inputfile-Generator, the filament element GID can be retrieved by knowing the filament number and the lesser node number.
           * Attention: If that property is changed, one should of course refer to node->Elements(). In the meanwhile, this method seems faster*/
          //int elegid = node0->Id() - i + (int)(floor((double)bspot*riseperbspot/lfil));
          int ival;
          // we need to make sure that bspot belongs to existing Filaments. In case a bspot is located on a node we always choose Element with the lower GID
          ival = (int)(floor((double)bspot*riseperbspot/elelength));
          // There might be a more elegant way
          if((int)(ceil((double)bspot*riseperbspot/elelength))==ival && bspot != 0)
            ival-=1;
          int elegid = node0->Id() - i + ival;

          // add element GID if this element is a row map element on this processor
          if(discret_->ElementRowMap()->LID(elegid)>-1)
          {
            // mark binding spot as being on this proc
            bspotonproc.push_back(1);
            bspot2element.push_back(elegid);


            // orientation [0;2pi] and parameter xi_bs of the binding spot
            double xibs = -1; // 0
            if(bspot==0)
            {
              bspotxi.push_back(xibs);
            }
            else
            {
              ival = (int)floor((double)bspot*riseperbspot/elelength);
              if((int)ceil((double)bspot*riseperbspot/elelength)==ival)
                  ival-=1;
              xibs = 2*((double)bspot*riseperbspot/elelength - ival - 0.5);
              bspotxi.push_back(xibs);
            }

          }
          else // only bspotonproc stores the information about off-proc elements since the structure can be easily rebuilt for the other vectors
            bspotonproc.push_back(0);
          // next binding spot on filament
          bspot++;
        }
      }
      // maps
      // create redundant binding spot column map based upon the bspotids
      bspotcolmap_ = rcp(new Epetra_Map(-1, (int)bspotgids.size(), &bspotgids[0], 0, discret_->Comm()));
      // create processor-specific row maps
      std::vector<int> bspotrowgids;
      for(int i=0; i<(int)bspotonproc.size(); i++)
        if(bspotonproc[i]==1)
          bspotrowgids.push_back(i);  // note: since column map is fully overlapping: i=col. LID = GID

      bspotrowmap_ = rcp(new Epetra_Map((int)bspotgids.size(), (int)bspotrowgids.size(), &bspotrowgids[0], 0, discret_->Comm()));

      // vectors
      // initialize class vectors
      bspotstatus_ = rcp(new Epetra_Vector(*bspotcolmap_, true));
      bspotstatus_->PutScalar(-1.0);
      bspotxi_ = rcp(new Epetra_Vector(*bspotcolmap_,true));
      bspot2element_ = rcp(new Epetra_Vector(*bspotrowmap_, true));
      bspot2nodes_=rcp(new Epetra_MultiVector(*bspotcolmap_,2,true));

      Epetra_Vector bspotxirow(*bspotrowmap_);
      Epetra_MultiVector bspot2nodesrow(*bspotrowmap_,2);
      for(int i=0; i<bspotrowmap_->NumMyElements(); i++)
      {
        bspotxirow[i] = bspotxi[i];
        (*bspot2element_)[i] = bspot2element[i];
        DRT::Element *ele = discret_->gElement(bspot2element[i]);
        for(int l=0;l<2;l++)
          bspot2nodesrow[l][i]=ele->NodeIds()[l];
      }
      // make information on orientations and curve parameters redundant on all procs
      Epetra_Import bspotimporter(*bspotcolmap_, *bspotrowmap_);
      bspotxi_->Import(bspotxirow,bspotimporter,Insert);
      bspot2nodes_->Import(bspot2nodesrow,bspotimporter,Insert);
    }
    else
    {
      bspotcolmap_ = Teuchos::rcp(new Epetra_Map(*(discret_->NodeColMap())));
      bspotrowmap_ = Teuchos::rcp(new Epetra_Map(*(discret_->NodeRowMap())));
      bspotstatus_ = Teuchos::rcp(new Epetra_Vector(*bspotcolmap_));
      bspotstatus_->PutScalar(-1.0);
    }
  }

  // create density-density-correlation-function map with
  if(DRT::INPUT::IntegralValue<INPAR::STATMECH::StatOutput>(statmechparams_, "SPECIAL_OUTPUT")==INPAR::STATMECH::statout_structanaly ||
     DRT::INPUT::IntegralValue<int>(statmechparams_, "GMSHOUTPUT"))
  {
    std::vector<int> bins;
    for(int i=0; i<discret_->Comm().NumProc()*numbins; i++)
      bins.push_back(i);
    ddcorrcolmap_ = Teuchos::rcp(new Epetra_Map(-1, discret_->Comm().NumProc()*numbins, &bins[0], 0, discret_->Comm()));
    // create processor-specific density-density-correlation-function map
    ddcorrrowmap_ = Teuchos::rcp(new Epetra_Map(discret_->Comm().NumProc()*numbins, 0, discret_->Comm()));
    // create new trafo matrix (for later use in DDCorr Function where we evaluate in layer directions), initialize with identity matrix
    trafo_ = Teuchos::rcp(new LINALG::SerialDenseMatrix(3,3,true));
    for(int i=0; i<trafo_->M(); i++)
      (*trafo_)(i,i) = 1.0;
    for(int i=0; i<(int)cog_.M(); i++)
      cog_(i) = periodlength_->at(i)/2.0;

    // start indices for parallel handling of orientation correlation etc in NumLinkerSpotsAndOrientation()
    // calculation of start indices for each processor
    // number of overall independent combinations
    int numnodes = discret_->NodeColMap()->NumMyElements();
    int numcombinations = (numnodes*numnodes-numnodes)/2;
    // combinations on each processor
    int combinationsperproc = (int)floor((double)numcombinations/(double)discret_->Comm().NumProc());
    int remainder = numcombinations%combinationsperproc;

    // get starting index tuples for later use
    startindex_->assign(2*discret_->Comm().NumProc(), 0.0);

    for(int mypid=0; mypid<discret_->Comm().NumProc()-1; mypid++)
    {
      std::vector<int> start(2,0);
      bool continueloop = false;
      bool quitloop = false;
      int counter = 0;
      int appendix = 0;
      if(mypid==discret_->Comm().NumProc()-1)
        appendix = remainder;

      // loop over crosslinker pairs
      for(int i=0; i<discret_->NodeColMap()->NumMyElements(); i++)
      {
        for(int j=0; j<discret_->NodeColMap()->NumMyElements(); j++)
        {
          if(i==(*startindex_)[2*mypid] && j==(*startindex_)[2*mypid+1])
            continueloop = true;
          if(j>i && continueloop)
          {
            if(counter<combinationsperproc+appendix)
              counter++;
            else
            {
              // new start index j
              if(j==discret_->NodeColMap()->NumMyElements()-1)
                start[1] = 0;
              else
                start[1] = j;
              quitloop = true;
              break;
            }
          }
        }
        if(quitloop)
        {
          // new start index i
          if(start[1]==0)
            start[0] = i+1;
          else
            start[0] = i;
          // new start tuple
          (*startindex_)[2*(mypid+1)] = (double)(start[0]);
          (*startindex_)[2*(mypid+1)+1] = (double)(start[1]);
          break;
        }
      }
    }

    if(DRT::INPUT::IntegralValue<INPAR::STATMECH::StatOutput>(statmechparams_, "SPECIAL_OUTPUT")==INPAR::STATMECH::statout_structanaly)
    {
      // initialize unbinding times vector to "-1"
      crosslinkunbindingtimes_ = Teuchos::rcp(new Epetra_MultiVector(*crosslinkermap_,2));
      crosslinkunbindingtimes_->PutScalar(-1.0);
    }
  }

  // in case the internal forces of the crosslinker affect the off-rate
  if(linkermodel_ == statmech_linker_bellseq || linkermodel_ == statmech_linker_bellseqintpol ||
     linkermodel_ == statmech_linker_active || linkermodel_ == statmech_linker_activeintpol)
  {
    element2crosslink_ = Teuchos::rcp(new Epetra_Vector(*(discret_->ElementRowMap())));
    element2crosslink_->PutScalar(-1.0);
    unbindingprobability_ = Teuchos::rcp(new Epetra_MultiVector(*crosslinkermap_,2,true));
  }
  else
  {
    element2crosslink_ = Teuchos::null;
    unbindingprobability_ = Teuchos::null;
  }

  std::vector<double> upperbound = *periodlength_;
  // handling both cases: with and without periodic boundary conditions
  if (periodlength_->at(0) == 0.0)
    for(int i=0; i<(int)upperbound.size(); i++)
      upperbound.at(i) = statmechparams_.get<double> ("MaxRandValue", 0.0);

  crosslinkerpositions_ = Teuchos::rcp(new Epetra_MultiVector(*crosslinkermap_, 3, true));

  Teuchos::RCP<Epetra_MultiVector> crosslinkerpositionstrans = Teuchos::rcp(new Epetra_MultiVector(*transfermap_,3,true));
  for (int i=0; i<crosslinkerpositionstrans->MyLength(); i++)
    for (int j=0; j<crosslinkerpositionstrans->NumVectors(); j++)
      (*crosslinkerpositionstrans)[j][i] = upperbound.at(j) * (*uniformgen_)();
  CommunicateMultiVector(crosslinkerpositionstrans, crosslinkerpositions_,false,true);

  // initial bonding status is set (no bonds)
  crosslinkerbond_ = Teuchos::rcp(new Epetra_MultiVector(*crosslinkermap_, 2));
  crosslinkerbond_->PutScalar(-1.0);
  // initial bond counter is set (no bonds)
  numbond_ = Teuchos::rcp(new Epetra_Vector(*crosslinkermap_, true));

  crosslinkonsamefilament_ = Teuchos::rcp(new Epetra_Vector(*crosslinkermap_, true));

  // crosslinker element IDs of the crosslink molecules
  crosslink2element_ = Teuchos::rcp(new Epetra_Vector(*crosslinkermap_));
  crosslink2element_->PutScalar(-1.0);

  // initialize the beautiful visuals vector (aka beevee-vector)
  visualizepositions_ = Teuchos::rcp(new Epetra_MultiVector(*crosslinkermap_, 3, false));
  for (int i=0; i<visualizepositions_->MyLength(); i++)
    for (int j=0; j<visualizepositions_->NumVectors(); j++)
      (*visualizepositions_)[j][i] = (*crosslinkerpositions_)[j][i];

  searchforneighbours_ = Teuchos::rcp(new Epetra_Vector(*crosslinkermap_, false));
  searchforneighbours_->PutScalar(1.0);

  if(linkermodel_ == statmech_linker_active || linkermodel_ == statmech_linker_activeintpol)
  {
    // set management vector to store the actual length of an active crosslinker: 0=long, 1=short
    // initial set long (=0)
    crosslinkeractlength_ = Teuchos::rcp(new Epetra_Vector(*crosslinkermap_));
    //set management vector to store the time since the last time the crosslinker was linked
    crosslinkeractcycletime_ = Teuchos::rcp(new Epetra_Vector(*crosslinkermap_));
    crosslinkeractcycletime_->PutScalar(statmechparams_.get<double>("ACTIVELINKERCYCLE",0.04));
  }

  // calculate  number of total binding spot depending on BSPOTINTERVAL
  numbspots_ = 0;
  int bspotinterval = statmechparams_.get<int>("BSPOTINTERVAL",1);
  if(bspotinterval>1)
  {
    for(int i=0; i<bspotstatus_->MyLength(); i++)
      if(i%bspotinterval==0)
        numbspots_++;
  }
  else if(bspotinterval==1)
    numbspots_ = bspotstatus_->MyLength();
  else
    dserror("Check your input file parameter BSPOTINTERVAL! It's %i !", bspotinterval);

#ifdef DEBUGCOUT
  std::cout<<"Proc "<<discret_->Comm().MyPID()<<": total number of binding spots = "<<numbspots_<<std::endl;
#endif

  return;
}//StatMechManager::CrosslinkerMoleculeInit

/*----------------------------------------------------------------------*
| Set crosslinkers wherever possible before the first time step        |
|                                                (private) mueller 11/11|
*----------------------------------------------------------------------*/
void STATMECH::StatMechManager::SetInitialCrosslinkers(Teuchos::RCP<CONTACT::Beam3cmanager> beamcmanager)
{
  if(linkermodel_ != statmech_linker_none)
  {
    if(beamcmanager!=Teuchos::null && (linkermodel_ == statmech_linker_stdintpol || linkermodel_ == statmech_linker_activeintpol || linkermodel_ == statmech_linker_bellseqintpol))
      dserror("Beam contact not implemented yet for the interpolated crosslinker element");

    const Epetra_Map noderowmap = *discret_->NodeRowMap();
    const Epetra_Map nodecolmap = *discret_->NodeColMap();
    // TODO
//==NEW
    // new node positions and rotations
    Teuchos::RCP<Epetra_MultiVector> bspotpositions = rcp(new Epetra_MultiVector(*bspotcolmap_,3,true));
    Teuchos::RCP<Epetra_MultiVector> bspotrotations = rcp(new Epetra_MultiVector(*bspotcolmap_,3,true));
    Epetra_Vector disrow(*discret_->DofRowMap(), true);
    Epetra_Vector discol(*discret_->DofColMap(), true);
    GetBindingSpotPositions(discol, bspotpositions, bspotrotations);

    Teuchos::RCP<Epetra_MultiVector> bspottriadscol = Teuchos::rcp(new Epetra_MultiVector(*bspotcolmap_,4,true));
    GetBindingSpotTriads(bspotrotations, bspottriadscol);
//==ENDNEW
//==OLD
    //node positions and rotations (binding spot positions and rotations when applying 4-noded beam element)
    std::map<int, LINALG::Matrix<3, 1> > currentpositions;

    {
      std::map<int, LINALG::Matrix<3, 1> > currentrotations;
      // Vectors hold zero->ok, since no displacement at this stage of the simulation
      GetNodePositionsFromDisVec(discol, currentpositions, currentrotations, true);
      // do initial octree build
      if(beamcmanager!=Teuchos::null)
        beamcmanager->OcTree()->OctTreeSearch(currentpositions);
    }
//==ENDOLD

    // generate a random order of binding spots and crosslinkers (only on Proc 0)
    std::vector<int> randbspot;
    std::vector<int> randlink;

    // consider only every BSPOTINTERVAL-th binding spot for binding
    int bspotinterval = statmechparams_.get<int>("BSPOTINTERVAL",1);

    //1. set singly bound linkers
    int numsinglybound = 0;
    if(discret_->Comm().MyPID()==0)
    {
      randbspot = Permutation(bspotcolmap_->NumMyElements());
      randlink = Permutation(statmechparams_.get<int>("N_crosslink", 0));
      int numbspots = statmechparams_.get<int>("INITOCCUPIEDBSPOTS",0);

      if(linkermodel_ != statmech_linker_none && numbspots>bspotcolmap_->NumMyElements())
        dserror("Given number of initially occupied binding spots (%i) exceeds the total binding spot count (%i)! Check your input file!",numbspots, bspotcolmap_->NumMyElements());

      // attach linkers with one binding domain to filaments
      if(networktype_ == statmech_network_motassay)
      {
        // attach linkers to substrate filaments
        int subfil = statmechparams_.get<int>("NUMSUBSTRATEFIL",0);
        int crosslinknum=0;
        for(int i=0; i<numbspots; i++)
        {
          // depending on whether binding spots are real nodes or virtual binding sites along a filament, we have to map the binding spot LID to its corresponding node LID
          int nodelid = i;
          if(linkermodel_ == statmech_linker_stdintpol || linkermodel_ == statmech_linker_activeintpol || linkermodel_ == statmech_linker_bellseqintpol)
            nodelid = (*bspot2nodes_)[0][i];
          if (((int)(*filamentnumber_)[nodelid])<subfil) //vector starts with zero element
          {
            // if this binding spot is still unoccupied
            if((*bspotstatus_)[i]<-0.9)
            {
              // attach it to the first binding spot (i.e. update of relevant class vectors)
              (*bspotstatus_)[i] = crosslinknum;
              (*numbond_)[crosslinknum] = 1.0;
              for(int j=0; j<crosslinkerbond_->NumVectors(); j++)
                if((*crosslinkerbond_)[j][crosslinknum]<-0.1)
                {
                  (*crosslinkerbond_)[j][crosslinknum] = bspotcolmap_->GID(i);
                  numsinglybound++;
                  break;
                }

              //update crosslinker position
              Epetra_SerialDenseMatrix LID(1,1);
              LID(0,0) = i;
              CrosslinkerIntermediateUpdate(*bspotpositions, LID, crosslinknum);
              // update visualization
              for (int j=0; j<visualizepositions_->NumVectors(); j++)
                (*visualizepositions_)[j][crosslinknum] = (*crosslinkerpositions_)[j][crosslinknum];
              crosslinknum++;
            }
          }
          else
            break;
        }
      }
      else
      {
        // first, establish specified number of singly bound crosslinkers
        for(int i=0; i<numbspots; i++)
        {
          int firstbspot = randbspot[i];
          // if this binding spot is still unoccupied
          if((*bspotstatus_)[firstbspot]<-0.9 && bspotcolmap_->GID(firstbspot)%bspotinterval==0)
          {
            // get the ilink-th random crosslinker
            int currlink = randlink[i];
            // attach it to the first binding spot (i.e. update of relevant class vectors)
            (*bspotstatus_)[firstbspot] = currlink;
            (*numbond_)[currlink] = 1.0;
            for(int j=0; j<crosslinkerbond_->NumVectors(); j++)
              if((*crosslinkerbond_)[j][currlink]<0.1)
              {
                (*crosslinkerbond_)[j][currlink] = bspotcolmap_->GID(firstbspot);
                numsinglybound++;
                break;
              }

            //update crosslinker position
            Epetra_SerialDenseMatrix LID(1,1);
            LID(0,0) = firstbspot;
            CrosslinkerIntermediateUpdate(*bspotpositions, LID, currlink);
            // update visualization
            for (int j=0; j<visualizepositions_->NumVectors(); j++)
              (*visualizepositions_)[j][currlink] = (*crosslinkerpositions_)[j][currlink];
          }
        }
      }
    }

    // Communication, first stage (not sure if we need to communicate all of these vectors here (or later)
    // transfer vectors
    Teuchos::RCP<Epetra_Vector> bspotstatusrow = Teuchos::rcp(new Epetra_Vector(*discret_->NodeRowMap(),true));
    Teuchos::RCP<Epetra_Vector> numbondtrans = Teuchos::rcp(new Epetra_Vector(*transfermap_, true));
    Teuchos::RCP<Epetra_MultiVector> crosslinkerpositionstrans = Teuchos::rcp(new Epetra_MultiVector(*transfermap_,3,true));
    Teuchos::RCP<Epetra_MultiVector> visualizepositionstrans = Teuchos::rcp(new Epetra_MultiVector(*transfermap_,3,true));
    Teuchos::RCP<Epetra_MultiVector> crosslinkerbondtrans = Teuchos::rcp(new Epetra_MultiVector(*transfermap_,2,true));

    CommunicateVector(bspotstatusrow, bspotstatus_);
    CommunicateVector(numbondtrans, numbond_);
    CommunicateMultiVector(crosslinkerpositionstrans, crosslinkerpositions_);
    CommunicateMultiVector(visualizepositionstrans, visualizepositions_);
    CommunicateMultiVector(crosslinkerbondtrans, crosslinkerbond_);

    // this section creates crosslinker finite elements
    // Commented because of performance issues  when adding a large number of elements at once.
    // 2. Now, parallely search for neighbour nodes
    Teuchos::RCP<Epetra_MultiVector> neighbourslid = Teuchos::null;
    if(statmechparams_.get<int>("SEARCHRES",1)>0)
      PartitioningAndSearch(*bspotpositions,*bspottriadscol, neighbourslid);

    // 3. create double bonds
    // a vector indicating the crosslink molecule which is going to constitute a crosslinker element
    RCP<Epetra_Vector> addcrosselement = Teuchos::rcp(new Epetra_Vector(*crosslinkermap_, true));
    int numsetelements = 0;

    if(discret_->Comm().MyPID()==0 && networktype_ != statmech_network_motassay)
    {
      int numbspots = statmechparams_.get<int>("INITOCCUPIEDBSPOTS",0);
      for(int i=0; i<numbspots; i++)
      {
        // get current linker
        int currlink = randlink[i];
        // obtain a random order of neighboursLID indices
        std::vector<int> neighbourorder = Permutation(neighbourslid->NumVectors());
        // if there is a potential second binding spot in the vicinity of the crosslinker (sufficient to check the first entry)
        for(int j=0; j<neighbourslid->NumVectors(); j++)
        {
          int currneighbour = neighbourorder[j];
          int secondbspot = (int)(*neighbourslid)[currneighbour][currlink];

          // if second binding binding spot exists and spot is unoccupied
          if((*neighbourslid)[currneighbour][currlink]>-0.1 && (*bspotstatus_)[secondbspot] < -0.9)// && bspotcolmap_->GID(secondbspot)%bspotinterval==0)
          {
            Epetra_SerialDenseMatrix LID(2,1);
            for(int k=0; k<crosslinkerbond_->NumVectors(); k++)
            {
              if((*crosslinkerbond_)[k][currlink]<-0.9)
              {
                LID(k,0) = secondbspot;
                (*crosslinkerbond_)[k][currlink] = bspotcolmap_->GID(secondbspot);
              }
              else
                LID(k,0) = (*crosslinkerbond_)[k][currlink];
            }

            // possibility of occupation of two binding spots on the same filament depends on K_ON_SELF
            if(statmechparams_.get<double>("K_ON_SELF",0.0)<1e-8)
            {
              if(linkermodel_ == statmech_linker_stdintpol  || linkermodel_ == statmech_linker_activeintpol || linkermodel_ == statmech_linker_bellseqintpol)
              {
                int filno1 = (*filamentnumber_)[discret_->NodeColMap()->LID((int)(*bspot2nodes_)[0][bspotcolmap_->GID((int)LID(0,0))])];
                int filno2 = (*filamentnumber_)[discret_->NodeColMap()->LID((int)(*bspot2nodes_)[0][bspotcolmap_->GID((int)LID(1,0))])];
                if( filno1 == filno2)
                  break;
              }
              else
              {
                if((*filamentnumber_)[(int)LID(0,0)]==(*filamentnumber_)[(int)LID(1,0)])
                  break;
              }
            }

            // check for intersection in case of beam contact (CURRENTLY ONLY FOR CONVENTIONAL BEAM3 ELEMENT)
            bool intersection = false;
            if(beamcmanager!=Teuchos::null)
            {
              Epetra_SerialDenseMatrix nodecoords(3,2);
              for(int k=0; k<nodecoords.M(); k++)
              {
                nodecoords(k,0) = (*bspotpositions)[k][(int)LID(0,0)];
                nodecoords(k,1) = (*bspotpositions)[k][(int)LID(1,0)];
              }
              for(int k=0; k<(int)LID.M(); k++)
              {
                intersection = beamcmanager->OcTree()->IntersectBBoxesWith(nodecoords, LID);
                if(intersection)
                  break;
              }
            }

            // orientation check
            //unit direction vector between currently considered two nodes
            LINALG::Matrix<3,1> direction;
            for(int j=0;j<3;j++)
              direction(j) = (*bspotpositions)[j][(int)LID(0,0)]-(*bspotpositions)[j][(int)LID(1,0)];
            direction.Scale(1.0/direction.Norm2());

            if(CheckOrientation(direction, *bspottriadscol,LID) && !intersection)
            {
              numsetelements++;
              (*addcrosselement)[currlink] = 1.0;
              // establish double bond to the first given neighbour
              // attach it to the second binding spot
              (*bspotstatus_)[secondbspot] = currlink;
              (*numbond_)[currlink] = 2.0;

              if(DRT::INPUT::IntegralValue<INPAR::STATMECH::StatOutput>(statmechparams_, "SPECIAL_OUTPUT")==INPAR::STATMECH::statout_structanaly)
                for(int k=0; k<crosslinkunbindingtimes_->NumVectors(); k++)
                (*crosslinkunbindingtimes_)[k][currlink] = 0.0;

              //update crosslinker position
              CrosslinkerIntermediateUpdate(*bspotpositions, LID, currlink);
              // update visualization
              for (int k=0; k<visualizepositions_->NumVectors(); k++)
                (*visualizepositions_)[k][currlink] = (*crosslinkerpositions_)[k][currlink];
              break;
            }
          }
        }
      }
    }

    // Communication, second stage
    RCP<Epetra_Vector> addcrosselementtrans = Teuchos::rcp(new Epetra_Vector(*transfermap_,true));
    bspotstatusrow->PutScalar(0.0);
    crosslinkerpositionstrans->PutScalar(0.0);
    visualizepositionstrans->PutScalar(0.0);
    crosslinkerbondtrans->PutScalar(0.0);
    numbondtrans->PutScalar(0.0);

    // export and reimport
    CommunicateVector(addcrosselementtrans, addcrosselement);
    CommunicateVector(bspotstatusrow, bspotstatus_);
    CommunicateVector(numbondtrans, numbond_);
    CommunicateMultiVector(crosslinkerpositionstrans, crosslinkerpositions_);
    CommunicateMultiVector(visualizepositionstrans, visualizepositions_);
    CommunicateMultiVector(crosslinkerbondtrans, crosslinkerbond_);

    int numsetelementsall = -1;
    discret_->Comm().MaxAll(&numsetelements, &numsetelementsall, 1);

    if(numsetelementsall>0)
    {
      // rotations from displacement vector
      Teuchos::RCP<Epetra_MultiVector> nodalrotations = Teuchos::null;
      // NODAL quaternions from filament elements
      Teuchos::RCP<Epetra_MultiVector> nodalquaternions = Teuchos::null;
      if(linkermodel_ == statmech_linker_stdintpol  || linkermodel_ == statmech_linker_activeintpol || linkermodel_ == statmech_linker_bellseqintpol)
      {
        nodalrotations = Teuchos::rcp(new Epetra_MultiVector(*(discret_->NodeColMap()),3));
        GetNodalBindingSpotPositionsFromDisVec(discol, Teuchos::null, nodalrotations);
        nodalquaternions = Teuchos::rcp(new Epetra_MultiVector(*(discret_->NodeColMap()),4));
        GetElementBindingSpotTriads(nodalquaternions);
      }

      // add elements to problem discretization (processor specific)
      for(int i=0; i<addcrosselement->MyLength(); i++)
      {
        if((*addcrosselement)[i]>0.9)
        {
          /*there is the problem of how to assign to a crosslinker element a GID which is certainly not used by any
           * other element; we know that for a network at the beginning (without crosslinkers) each node has a
           * connectivity of 1 or 2 so that the number of all elements used to discretize the filaments is smaller
           * than basisnodes_. Thus we choose crosslinker GIDs >= basisnodes_ for the additionally added crosslinker
           * elements. To make sure that two crosslinkers never have the same GID we add to basisnodes_ the value
           * basisnodes_*GID1, where GID1 is the GID of the first node of the crosslinker element. Then we add GID2
           * (note: GID2 < basisnodes_), which is the GID of the second node of the crosslinker element.
           * Hence basisnodes_ + GID1*basisnodes_ + GID2 always gives a GID which cannot be used by any other element;
           * the first node of the crosslinker element is always assumed to be the one with the greater GID; the owner
           * of the first node is assumed to be the owner of the crosslinker element*/

          // copied binding spot map (in order to handle the addition of elements based upon an invariant map just in case)
          Epetra_Map bspotcolmap(*bspotcolmap_);

          // obtain binding spot GID
          std::vector<int> bspotgid(2);
          // determine smaller and larger of the GIDs
          bspotgid.at(1) = std::min((int)(*crosslinkerbond_)[0][i],(int)(*crosslinkerbond_)[1][i]);
          bspotgid.at(0) = std::max((int)(*crosslinkerbond_)[0][i],(int)(*crosslinkerbond_)[1][i]);

          // different sizes due to different linker elements
          Teuchos::RCP<std::vector<int> > globalnodeids;
          if(linkermodel_ == statmech_linker_stdintpol || linkermodel_ == statmech_linker_activeintpol || linkermodel_ == statmech_linker_bellseqintpol)
          {
            globalnodeids = Teuchos::rcp(new std::vector<int>(4,0));
            for(int l=0;l<2;l++)
            {
              (*globalnodeids)[l]=(int)(*bspot2nodes_)[l][bspotgid[0]];
              (*globalnodeids)[l+2]=(int)(*bspot2nodes_)[l][bspotgid[1]];
            }
          }
          else
          {
            globalnodeids = Teuchos::rcp(new std::vector<int>(2,0));
            (*globalnodeids)[0] = bspotgid[0];
            (*globalnodeids)[1] = bspotgid[1];
          }

          // calculate element GID
          int newcrosslinkerGID = (bspotgid[0] + 1)*basisnodes_ + bspotgid[1];

          /* Correction of the crosslinker GID if there happens to be another crosslinker with the same ID. This might occur
           * if two different crosslink molecules bind to the same two nodes. Then, the calculated crosslinker GID will be the
           * same in both cases, leading to problems within the discretization.
           * As long as an unused GID cannot be found, the crosslinker GID keeps getting incremented by 1.*/
          discret_->Comm().Barrier();
          while(1)
          {
            int gidexists = 1;
            // query existance of node on this Proc
            int gidonproc = (int)(discret_->HaveGlobalElement(newcrosslinkerGID));
            // sum over all processors
            discret_->Comm().MaxAll(&gidonproc, &gidexists, 1);
            // calculate new GID if necessary by shifting the initial GID
            if(gidexists>0)
              newcrosslinkerGID++;
            else
              break;
          }
          /* Create mapping from crosslink molecule to crosslinker element GID
           * Note: No need for the usual procedure of exporting and reimporting to make things redundant
           * because info IS already redundant by design here.*/
          (*crosslink2element_)[i] = newcrosslinkerGID;

          //save positions of nodes between which a crosslinker has to be established in variables xrefe and rotrefe:
          std::vector<double> rotrefe(6);
          std::vector<double> xrefe(6);
          // resize in case of interpolated crosslinker element
          if (linkermodel_ == statmech_linker_stdintpol  || linkermodel_ == statmech_linker_activeintpol || linkermodel_ == statmech_linker_bellseqintpol)
            rotrefe.resize(12);

          for(int k=0; k<3; k++)
          {
            xrefe[k ] = (*bspotpositions)[k][bspotcolmap.LID(bspotgid.at(0))];
            xrefe[k+3] = (*bspotpositions)[k][bspotcolmap.LID(bspotgid.at(1))];

            //set nodal rotations (not true ones, only those given in the displacement vector)
            if (linkermodel_ == statmech_linker_stdintpol  || linkermodel_ == statmech_linker_activeintpol || linkermodel_ == statmech_linker_bellseqintpol)
            {
              rotrefe[k]   = (*nodalrotations)[k][bspotcolmap.LID((*globalnodeids)[0])];
              rotrefe[k+3] = (*nodalrotations)[k][bspotcolmap.LID((*globalnodeids)[1])];
              rotrefe[k+6] = (*nodalrotations)[k][bspotcolmap.LID((*globalnodeids)[2])];
              rotrefe[k+9] = (*nodalrotations)[k][bspotcolmap.LID((*globalnodeids)[3])];
            }
            else
            {
              rotrefe[k] = (*bspotrotations)[k][bspotcolmap.LID((*globalnodeids)[0])];
              rotrefe[k+3] = (*bspotrotations)[k][bspotcolmap.LID((*globalnodeids)[1])];
            }
          }

          /*a crosslinker is added on each processor which is row node owner of at least one of its nodes;
           *the processor which is row map owner of the node with the larger GID will be the owner of the new crosslinker element; on the other processors the new
           *crosslinker element will be present as a ghost element only*/
          bool hasrownode = false;
          for(int k=0; k<(int)globalnodeids->size(); k++)
            if(noderowmap.LID(globalnodeids->at(k)) > -1)
            {
              hasrownode = true;
              break;
            }

          if(hasrownode)
            AddNewCrosslinkerElement(newcrosslinkerGID, globalnodeids,bspotgid, xrefe,rotrefe,*discret_,nodalquaternions);

          // add all new elements to contact discretization on all Procs
          if(beamcmanager!=Teuchos::null)
            AddNewCrosslinkerElement(newcrosslinkerGID, globalnodeids,bspotgid, xrefe,rotrefe,beamcmanager->ContactDiscret(),nodalquaternions);
        }
      }
    }

    // synchronization for problem discretization
    discret_->CheckFilledGlobally();
    discret_->FillComplete(true, false, false);

    if(beamcmanager!=Teuchos::null)
    {
      beamcmanager->ContactDiscret().CheckFilledGlobally();
      beamcmanager->ContactDiscret().FillComplete(true, false, false);
    }

    // reduce number of contact pairs to zero again to avoid unnecessary computations
    if(beamcmanager!=Teuchos::null && !(DRT::INPUT::IntegralValue<int>(DRT::Problem::Instance()->ContactDynamicParams(), "BEAMS_NEWGAP")))
      beamcmanager->ResetPairs();

    //Gmsh output
    if(DRT::INPUT::IntegralValue<int>(statmechparams_,"GMSHOUTPUT"))
    {
      std::ostringstream filename;
      filename << StatMechRootPath() <<"/GmshOutput/InitLinks.pos";
      Epetra_Vector disrow(*discret_->DofRowMap(), true);
      GmshOutput(disrow,filename,0);
    }
    if(beamcmanager!=Teuchos::null && DRT::INPUT::IntegralValue<INPAR::STATMECH::StatOutput>(statmechparams_, "SPECIAL_OUTPUT")==INPAR::STATMECH::statout_octree)
    {
      // "-2" for initial octree output
      beamcmanager->OcTree()->OctTreeSearch(currentpositions,-2);
      beamcmanager->ResetPairs();
    }
    //std::couts
    if(!discret_->Comm().MyPID())
    {
      std::cout<<"====setting initial crosslinkers===="<<std::endl;
      std::cout<<"singly bound: "<<numsinglybound-numsetelements<<std::endl;
      std::cout<<"doubly bound: "<<numsetelements<<std::endl;
      std::cout<<"===================================="<<std::endl;
    }
  }
  return;
} // SetInitialCrosslinkers()

/*----------------------------------------------------------------------*
 |  Periodic Boundary Shift for crosslinker diffusion simulation        |
 |                                                (public) mueller 07/10|
 *----------------------------------------------------------------------*/
void STATMECH::StatMechManager::CrosslinkerPeriodicBoundaryShift(Teuchos::RCP<Epetra_MultiVector> crosslinkerpositions)
{
  for (int i=0; i<crosslinkerpositions->NumVectors(); i++)
    for (int j=0; j<crosslinkerpositions->MyLength(); j++)
    {
      if ((*crosslinkerpositions)[i][j] > (*periodlength_)[i])
        (*crosslinkerpositions)[i][j] -= (floor((*crosslinkerpositions)[i][j]/(*periodlength_)[i]))*(*periodlength_)[i];
      if ((*crosslinkerpositions)[i][j] < 0.0)
        (*crosslinkerpositions)[i][j] -= (floor((*crosslinkerpositions)[i][j]/(*periodlength_)[i]))*(*periodlength_)[i];
    }
  return;
}// StatMechManager::CrosslinkerPeriodicBoundaryShift

/*----------------------------------------------------------------------*
 |  Evaluate DBCs in case of periodic BCs (public)         mueller 03/12|
 *----------------------------------------------------------------------*/
void STATMECH::StatMechManager::EvaluateDirichletStatMech(Teuchos::ParameterList&            params,
                                                          Teuchos::RCP<Epetra_Vector>        dis,
                                                          Teuchos::RCP<LINALG::MapExtractor> dbcmapextractor)
{
  if(dbcmapextractor!=Teuchos::null)
  {
    INPAR::STATMECH::DBCType dbctype = DRT::INPUT::IntegralValue<INPAR::STATMECH::DBCType>(statmechparams_,"DBCTYPE");

    switch(dbctype)
    {
      // standard DBC application
      case INPAR::STATMECH::dbctype_std:
        discret_->EvaluateDirichlet(params, dis, Teuchos::null, Teuchos::null, Teuchos::null, dbcmapextractor);
      break;
      // default: everything involving periodic boundary conditions
      default:
        EvaluateDirichletPeriodic(params, dis, dbcmapextractor);
      break;
    }
  }
  else
    dserror("Only new DBC application method implemented using the map extractor! Old version using toggle vectors dicontinued!");

  return;
}
/*----------------------------------------------------------------------*
 |  Evaluate DBCs in case of periodic BCs (private         mueller 02/10|
 *----------------------------------------------------------------------*/
void STATMECH::StatMechManager::EvaluateDirichletPeriodic(Teuchos::ParameterList& params,
                                                          Teuchos::RCP<Epetra_Vector> dis,
                                                          Teuchos::RCP<LINALG::MapExtractor> dbcmapextractor)
{
#ifdef MEASURETIME
  const double t_start = Teuchos::Time::wallTime();
#endif // #ifdef MEASURETIME

  // some preliminary checks
  if (!(discret_->Filled())) dserror("FillComplete() was not called");
  if (!(discret_->HaveDofs())) dserror("AssignDegreesOfFreedom() was not called");

  // retrieve DBC type
  INPAR::STATMECH::DBCType dbctype = DRT::INPUT::IntegralValue<INPAR::STATMECH::DBCType>(statmechparams_,"DBCTYPE");
  // vector holding Dirichlet increments
  Teuchos::RCP<Epetra_Vector> deltadbc = Teuchos::rcp(new Epetra_Vector(*(discret_->DofRowMap()),true));
  // vector of DOF-IDs which are Dirichlet BCs
  Teuchos::RCP<std::set<int> > dbcgids = Teuchos::rcp(new std::set<int>());

  // take action according to DBC type in the input file
  switch(dbctype)
  {
    // shear with a fixed Dirichlet node set
    case INPAR::STATMECH::dbctype_shearfixed:
    {
      // if point in time is reached at which the application of Dirichlet values starts
      if(!DBCStart(params))
        return;
      DBCOscillatoryMotion(params,dis,deltadbc);
      useinitdbcset_ = true;
    }
    break;
    // shear with a fixed Dirichlet node set and setting axial stiffness of dirsupted elements close to zero
    case INPAR::STATMECH::dbctype_shearfixeddel:
    {
      // if point in time is reached at which the application of Dirichlet values starts
      if(!DBCStart(params))
        return;
      DBCOscillatoryMotion(params,dis,deltadbc);
      useinitdbcset_ = true;
    }
    break;
    // shear with an updated Dirichlet node set (only DOF in direction of oscillation is subject to BC, others free)
    case INPAR::STATMECH::dbctype_sheartrans:
    {
      if(!DBCStart(params))
        return;
      DBCOscillatoryMotion(params,dis,deltadbc);
    }
    break;
    // pin down and release individual nodes
    case INPAR::STATMECH::dbctype_pinnodes:
    {
      // apply DBCs from input file first
      DBCGetPredefinedConditions(params,dis,dbcgids);
      // then, determine nodes that get inhibited by Crosslinkers
      DBCPinNodes();
    }
    break;
    // apply affine shear displacement to all nodes
    case INPAR::STATMECH::dbctype_affineshear:
    {
      if(!DBCStart(params))
        return;
      DBCAffineShear(params,dis,deltadbc);
      // used here to store the node set that remains under DBCs after the initial affine displacement
      useinitdbcset_ = true;
    }
    break;
    // apply affine shear displacement to all nodes
    case INPAR::STATMECH::dbctype_affinesheardel:
    {
      if(!DBCStart(params))
        return;
      DBCAffineShear(params,dis,deltadbc);
      // used here to store the node set that remains under DBCs after the initial affine displacement
      useinitdbcset_ = true;
    }
    break;
    default:
      return;
  }
//------------------------------------set Dirichlet values
  DBCSetValues(dis, deltadbc, dbcgids);

//----create DBC and free map and build their common extractor
  DBCCreateMap(dbcgids, dbcmapextractor);

#ifdef MEASURETIME
  const double t_end = Teuchos::Time::wallTime();
  if(!discret_->Comm().MyPID())
    std::cout<<"DBC Evaluation time: "<<t_end-t_start<<std::endl;
#endif // #ifdef MEASURETIME

  return;
} //EvaluateDirichletPeriodic()

/*----------------------------------------------------------------------*
 | Determine if application of DBCs starts at a given time  mueller 5/12|
 *----------------------------------------------------------------------*/
bool STATMECH::StatMechManager::DBCStart(Teuchos::ParameterList& params)
{
  double eps = 2.0e-11;
  // get the current time
  double time = params.get<double>("total time", 0.0);
  double starttime = actiontime_->at(dbctimeindex_);
  //double dt = params.get<double>("delta time", 0.01);
  if (time<0.0) dserror("t = %f ! Something is utterly wrong here. The absolute time should be positive!", time);

  if(time < starttime && fabs(starttime-time)>eps)
    return false;
  else
    return true;
}

/*----------------------------------------------------------------------*
 | Some simple sanity checks                   (private)  mueller 01/13 |
 *----------------------------------------------------------------------*/
void STATMECH::StatMechManager::DBCSanityCheck()
{
  if(dbctimeindex_<0) // default
    dbctimeindex_ = (int)actiontime_->size()-1;
  else if(dbctimeindex_==0)
    dserror("Given index DBCTIMEINDEX = %i ! Start counting at 1!", dbctimeindex_);
  else if(dbctimeindex_>(int)actiontime_->size())
    dserror("Given index DBCTIMEINDEX = %i lies outside the ACTIONTIME vector! Check your input file!", dbctimeindex_);
  else
    dbctimeindex_--;

  INPAR::STATMECH::DBCType dbctype = DRT::INPUT::IntegralValue<INPAR::STATMECH::DBCType>(statmechparams_,"DBCTYPE");
  switch(dbctype)
  {
    case INPAR::STATMECH::dbctype_std:
      break;
    case INPAR::STATMECH::dbctype_pinnodes:
      break;
    case INPAR::STATMECH::dbctype_none:
      break;
    // default: everything involving periodic boundary conditions and rheology
    default:
    {
      // "-1" because of counting convention
      int curvenumber = statmechparams_.get<int>("CURVENUMBER",0)-1;
      int displacementdir = statmechparams_.get<int>("DBCDISPDIR",-1)-1;

      if(periodlength_->at(0)<=0.0)
        dserror("PERIODLENGTH = %f! This method works only in conjunction with periodic boundary conditions!", periodlength_->at(0));
      if(periodlength_->at(0)!=periodlength_->at(1) || periodlength_->at(0)!=periodlength_->at(2) || periodlength_->at(1)!=periodlength_->at(2))
        dserror("Check PERIODLENGTH in your input file! This method only works for a cubic boundary volume");
      if(displacementdir>2 || displacementdir<0 || curvenumber<0)
        dserror("In case of imposed DBC values, please define the StatMech parameters DBCDISPDIR={1,2,3} and/or CURVENUMBER correctly");
      if(dbctype==INPAR::STATMECH::dbctype_affineshear)
      {
        if((int)actiontime_->size()<3)
          dserror("For affine deformation, give three time values for ACTIONTIME! 1. equilibration time 2. dbc application time 3. release of DBC nodes other than nodes of elements shifted in z-direction!");
        else if(dbctimeindex_!=1)
          dserror("DBCTIMEINDEX must be set to 2 for affine shear displacement!");
      }
    }
  }
  return;
}

/*----------------------------------------------------------------------*
 | Gather information on where DBCs are to be applied in the case of    |
 | viscoelastic measurements                   (private)  mueller 05/12 |
 *----------------------------------------------------------------------*/
void STATMECH::StatMechManager::DBCOscillatoryMotion(Teuchos::ParameterList& params,
                                                     Teuchos::RCP<Epetra_Vector> dis,
                                                     Teuchos::RCP<Epetra_Vector> deltadbc)
{
  INPAR::STATMECH::DBCType dbctype = DRT::INPUT::IntegralValue<INPAR::STATMECH::DBCType>(statmechparams_,"DBCTYPE");
  // get the absolute time, time step size, time curve number and oscillation direction
  const double time = params.get<double>("total time", 0.0); // target time (i.e. timen_)
  double dt = params.get<double>("delta time", 0.01);
  int curvenumber = statmechparams_.get<int>("CURVENUMBER",0)-1;
  int oscdir = statmechparams_.get<int>("DBCDISPDIR",-1)-1;
  double amp = statmechparams_.get<double>("SHEARAMPLITUDE",0.0);

//------------------------------------------------------gather dbc node set(s)
  // after the the very first entry into the following part, reenter only if the Dirichlet Node Set is dynamically updated
  if(!useinitdbcset_)
  {
    // note: we need column map displacements as we might need access to ghost nodes
    Teuchos::RCP<Epetra_Vector> discol = Teuchos::rcp(new Epetra_Vector(*discret_->DofColMap(), true));
    LINALG::Export(*dis, *discol);

    // clear node sets, as they will be filled for the first time or updated with new node sets
    dbcnodesets_.clear();
    // vectors to manipulate DBC
    std::vector<int> oscillnodes;
    std::vector<int> fixednodes;
    oscillnodes.clear();
    fixednodes.clear();
    // toggle vector
    std::vector<bool> nodedbcstatus(discret_->NumMyColNodes(), false);

    for(int lid=0; lid<discret_->NumMyRowElements(); lid++)
    {
      // An element used to browse through local Row Elements
      DRT::Element* element = discret_->lRowElement(lid);
      // skip element if it is a crosslinker element or in addition, in case of the Bead Spring model, Torsion3 elements
      if(element->Id() <= basisnodes_ && element->Id() < statmechparams_.get<int>("NUM_EVAL_ELEMENTS", basiselements_))
      {
        // number of translational DOFs
        int numdof = 3;
        // positions of nodes of an element with n nodes
        LINALG::SerialDenseMatrix coord(numdof,(int)discret_->lRowElement(lid)->NumNode(), true);
        // indicates location, direction and component of a broken element with n nodes->n-1 possible cuts
        LINALG::SerialDenseMatrix cut(numdof,(int)discret_->lRowElement(lid)->NumNode()-1,true);
//-------------------------------- obtain nodal coordinates of the current element
        std::vector<int> doflids;
        GetElementNodeCoords(element, discol, coord, &doflids);
//-----------------------detect broken/fixed/free elements and fill position vector
        // determine existence and location of broken element
        if(CheckForBrokenElement(coord,cut))
        {
          // reduce the axial stiffness of the element drastically, close to zero in order to take this element out
          // only for Beam3ii case
          if(element->ElementType()==DRT::ELEMENTS::Beam3iiType::Instance() && dbctype == INPAR::STATMECH::dbctype_shearfixeddel)
          {
            for(int n=0; n<cut.N(); n++)
            {
              if(cut(2,n)>0.0)
              {
                //std::cout<<"Element "<<element->Id()<<" now has a reduced cross section"<<std::endl;
                dynamic_cast<DRT::ELEMENTS::Beam3ii*>(element)->SetCrossSec(1.0e-9);
                dynamic_cast<DRT::ELEMENTS::Beam3ii*>(element)->SetCrossSecShear(1.1e-9);
                break;
              }
            }
          }
          // loop over number of cuts (columns)
          for(int n=0; n<cut.N(); n++)
          {
            int nodelidn = discret_->NodeColMap()->LID(element->Nodes()[n]->Id());
            int nodelidnp = discret_->NodeColMap()->LID(element->Nodes()[n+1]->Id());
            int noderowlidn = discret_->NodeRowMap()->LID(element->Nodes()[n]->Id());
            int noderowlidnp = discret_->NodeRowMap()->LID(element->Nodes()[n+1]->Id());

            // only consider this node pair if it is not already subject to Dirichlet BCs
            if(nodedbcstatus.at(nodelidn)==false && nodedbcstatus.at(nodelidnp)==false)
            {
              // case 1: broken element (in z-dir); node_n+1 oscillates, node_n is fixed in dir. of oscillation
              if(cut(2,n)==1.0)
              {
                // update dbc status for these nodes (dofs)
                nodedbcstatus.at(nodelidn) = true;
                nodedbcstatus.at(nodelidnp) = true;
                // add GID of fixed node to fixed-nodes-vector (to be added to condition later)
                if(noderowlidn>-1)
                  fixednodes.push_back(element->Nodes()[n]->Id());
                // add GID of oscillating node to osc.-nodes-vector
                if(noderowlidnp>-1)
                {
                  oscillnodes.push_back(element->Nodes()[n+1]->Id());

                  // incremental Dirichlet displacement for an oscillating node (all DOFs except oscdir = 0.0)
                  // time curve increment
                  double tcincrement = 0.0;
                  if(curvenumber>-1)
                    tcincrement = DRT::Problem::Instance()->Curve(curvenumber).f(time) -
                                  DRT::Problem::Instance()->Curve(curvenumber).f(time-dt);
                  (*deltadbc)[doflids.at(numdof*(n+1)+oscdir)] = amp*tcincrement;
                }
              }// end of case 1
              if(cut(2,n)==2.0)// case 2: broken element (in z-dir); node_n oscillates, node_n+1 is fixed in dir. of oscillation
              {
                nodedbcstatus.at(nodelidn) = true;
                nodedbcstatus.at(nodelidnp) = true;
                if(noderowlidnp>-1)
                  fixednodes.push_back(element->Nodes()[n+1]->Id());
                // oscillating node
                if(noderowlidn>-1)
                {
                  oscillnodes.push_back(element->Nodes()[n]->Id());
                  double tcincrement = 0.0;
                  if(curvenumber>-1)
                    tcincrement = DRT::Problem::Instance()->Curve(curvenumber).f(time) -
                                  DRT::Problem::Instance()->Curve(curvenumber).f(time-dt);
                  (*deltadbc)[doflids.at(numdof*n+oscdir)] = amp*tcincrement;
                }
              } // end of case 2
            }
          }
        }
      }
    }
//-----------------------------------------update node sets
    dbcnodesets_.push_back(oscillnodes);
    dbcnodesets_.push_back(fixednodes);
//---------check/set force sensors anew for each time step
    UpdateForceSensors(dbcnodesets_[0], oscdir);
  }
  else // only ever use the fixed set of nodes after the first time dirichlet values were imposed
  {
    // dofs of fixednodes_ remain untouched since fixednodes_ dofs are 0.0
    for(int i=0; i<(int)dbcnodesets_[0].size(); i++)
    {
      int nodelid = discret_->NodeRowMap()->LID(dbcnodesets_[0][i]);
      DRT::Node* oscnode = discret_->lRowNode(nodelid);
      std::vector<int> dofnode = discret_->Dof(oscnode);
      // oscillating node
      double tcincrement = 0.0;
      if(curvenumber>-1)
        tcincrement = DRT::Problem::Instance()->Curve(curvenumber).f(time) -
                      DRT::Problem::Instance()->Curve(curvenumber).f(time-dt);
      (*deltadbc)[discret_->DofRowMap()->LID(dofnode[oscdir])] = amp*tcincrement;
    }
  }
  return;
}

/*----------------------------------------------------------------------*
 | pins down a particular given node by Dirichlet BCs                   |
 |                                           (private)     mueller 05/12|
 *----------------------------------------------------------------------*/
void STATMECH::StatMechManager::DBCPinNodes()
{
  // gather node sets
  // end nodes of horizontal filaments
  std::vector<int> hfilfixednodes;
  // end nodes of vertical filaments
  std::vector<int> vfilfixednodes;
  // regular nodes
  std::vector<int> fixednodes;

  // clear sets initially
  dbcnodesets_.clear();
  hfilfixednodes.clear();
  vfilfixednodes.clear();
  fixednodes.clear();

  int maxnodegid = discret_->NumGlobalNodes()-1;

  for(int i=0; i<bspotstatus_->MyLength(); i++)
  {
    if((*bspotstatus_)[i]>-0.1)
    {
      int nodegid = bspotcolmap_->GID(i);
      if(nodegid==0)
      {
        hfilfixednodes.push_back(nodegid);
        continue;
      }
      // the very last node
      if(nodegid==maxnodegid)
      {
        if((*filamentnumber_)[bspotcolmap_->LID(maxnodegid)]==0)
          hfilfixednodes.push_back(nodegid);
        else
          vfilfixednodes.push_back(nodegid);
        continue;
      }
      // in between end nodes and regular nodes
      if(nodegid<maxnodegid && nodegid>0)
      {
        // look for end nodes
        if((*filamentnumber_)[i-1]==(*filamentnumber_)[i] && (*filamentnumber_)[i+1]!=(*filamentnumber_)[i])
        {
          if((*filamentnumber_)[i]==0)
            hfilfixednodes.push_back(bspotcolmap_->GID(i));
          else
            vfilfixednodes.push_back(bspotcolmap_->GID(i));
        }
        else
          fixednodes.push_back(bspotcolmap_->GID(i));
      }
    }
  }

  // add to node sets
  dbcnodesets_.push_back(hfilfixednodes);
  dbcnodesets_.push_back(vfilfixednodes);
  dbcnodesets_.push_back(fixednodes);

  return;
}

/*----------------------------------------------------------------------*
 | Apply affine shear displacement              (private)  mueller 05/12|
 *----------------------------------------------------------------------*/
void STATMECH::StatMechManager::DBCAffineShear(Teuchos::ParameterList&     params,
                                               Teuchos::RCP<Epetra_Vector> dis,
                                               Teuchos::RCP<Epetra_Vector> deltadbc)
{
  /*====the following is mostly taken from DBCOscillatoryMotion()====
   * we apply an affine deformation to the whole network, i.e. all nodes
   * are displaced according to their vertical position in order to
   * shear the network by an angle gamma*/

  const double time = params.get<double>("total time", 0.0);
  double dt = params.get<double>("delta time", 0.01);
  int curvenumber = statmechparams_.get<int>("CURVENUMBER",0)-1;
  int displacementdir = statmechparams_.get<int>("DBCDISPDIR",-1)-1;
  double amp = statmechparams_.get<double>("SHEARAMPLITUDE",0.0);
  double tcincrement = DRT::Problem::Instance()->Curve(curvenumber).f(time) -
                       DRT::Problem::Instance()->Curve(curvenumber).f(time-dt);


  double tol = 1e-10;

  // upon hitting the third time threshold, erase the last node set and thus free these nodes of DBCs
  if(fabs(time-dt-actiontime_->at(2))<tol && (int) dbcnodesets_.size()==3)
    dbcnodesets_.erase(dbcnodesets_.end());

  // We need column map displacements as we might need access to ghost nodes
  Teuchos::RCP<Epetra_Vector> discol = Teuchos::rcp(new Epetra_Vector(*discret_->DofColMap(), true));
  LINALG::Export(*dis, *discol);

  // reenter once the third time threshold is hit in order to release the
  if(!useinitdbcset_)
  {
    dbcnodesets_.clear();
    std::vector<int> sensornodes;
    std::vector<int> fixednodes;
    std::vector<int> affineshearnodes;
    sensornodes.clear();
    fixednodes.clear();
    affineshearnodes.clear();
    std::vector<bool> nodedbcstatus(discret_->NumMyRowNodes(), false);

    // displaced and fixed nodes
    for(int lid=0; lid<discret_->NumMyRowElements(); lid++)
    {
      DRT::Element* element = discret_->lRowElement(lid);
      // skip element if it is a crosslinker element or in addition, in case of the Bead Spring model, Torsion3 elements
      if(element->Id() <= basisnodes_ && element->Id() < statmechparams_.get<int>("NUM_EVAL_ELEMENTS", basiselements_))
      {
        // number of translational DOFs
        int numdof = 3;
        LINALG::SerialDenseMatrix coord(numdof,(int)discret_->lRowElement(lid)->NumNode(), true);
        // indicates location, direction and component of a broken element with n nodes->n-1 possible cuts
        LINALG::SerialDenseMatrix cut(numdof,(int)discret_->lRowElement(lid)->NumNode()-1,true);

        std::vector<int> doflids;
        GetElementNodeCoords(element, discol, coord, &doflids);

        // shifted elements
        if(CheckForBrokenElement(coord,cut))
        {
          // reduce the axial stiffness of the element drastically, close to zero in order to take this element out
          // only for Beam3ii case
          for(int n=0; n<cut.N(); n++)
            if(element->ElementType()==DRT::ELEMENTS::Beam3iiType::Instance() && cut(2,n)>0.0)
            {
              dynamic_cast<DRT::ELEMENTS::Beam3ii*>(element)->SetCrossSec(1.0e-12);
              dynamic_cast<DRT::ELEMENTS::Beam3ii*>(element)->SetCrossSecShear(1.1e-12);
              break;
            }

          for(int n=0; n<cut.N(); n++)
          {
            int nodelidn = discret_->NodeRowMap()->LID(element->Nodes()[n]->Id());
            //std::cout<<"Proc "<<discret_->Comm().MyPID()<<": nodelidn = "<<nodelidn<<"/"<<discret_->NodeRowMap()->GID(nodelidn)<<std::endl;
            if(nodelidn>-1)
            {
              //std::cout<<"Proc "<<discret_->Comm().MyPID()<<": nodelidn cut = "<<(int)cut(2,n)<<std::endl;
              if(nodedbcstatus.at(nodelidn)==false)
              {
                switch((int)(round(cut(2,n))))
                {
                  // case 1: broken element (in z-dir); node_n+1 displaced, node_n is fixed in dir. of displacement
                  case 1:
                  {
                    //std::cout<<"Proc "<<discret_->Comm().MyPID()<<": case 1 n"<<std::endl;
                    nodedbcstatus.at(nodelidn) = true;
                    fixednodes.push_back(element->Nodes()[n]->Id());
                  }
                  break;
                  // case 2: broken element (in z-dir); node_n displaced, node_n+1 is fixed in dir. of displacement
                  case 2:
                  {
                    //std::cout<<"Proc "<<discret_->Comm().MyPID()<<": case 2 n"<<std::endl;
                    nodedbcstatus.at(nodelidn) = true;
                    sensornodes.push_back(element->Nodes()[n]->Id());
                    (*deltadbc)[doflids.at(numdof*n+displacementdir)] =  (coord(2,n)/(*periodlength_)[2])*amp*tcincrement;
                  }
                  break;
                }
              }
            }
            int nodelidnp = discret_->NodeRowMap()->LID(element->Nodes()[n+1]->Id());
            //std::cout<<"Proc "<<discret_->Comm().MyPID()<<": nodelidnp = "<<nodelidnp<<"/"<<discret_->NodeRowMap()->GID(nodelidnp)<<std::endl;
            if(nodelidnp>-1)
            {
              //std::cout<<"Proc "<<discret_->Comm().MyPID()<<": nodelidnp cut = "<<(int)cut(2,n)<<std::endl;
              if(nodedbcstatus.at(nodelidnp)==false)
              {
                switch((int)(round(cut(2,n))))
                {
                  // case 1: broken element (in z-dir); node_n+1 displaced, node_n is fixed in dir. of displacement
                  case 1:
                  {
                    //std::cout<<"Proc "<<discret_->Comm().MyPID()<<": case 1 np"<<std::endl;
                    nodedbcstatus.at(nodelidnp) = true;
                    sensornodes.push_back(element->Nodes()[n+1]->Id());
                    (*deltadbc)[doflids.at(numdof*(n+1)+displacementdir)] = (coord(2,n+1)/(*periodlength_)[2])*amp*tcincrement;
                  }
                  break;
                  // case 2: broken element (in z-dir); node_n displaced, node_n+1 is fixed in dir. of displacement
                  case 2:
                  {
                   // std::cout<<"Proc "<<discret_->Comm().MyPID()<<": case 2 np"<<std::endl;
                    nodedbcstatus.at(nodelidnp) = true;
                    fixednodes.push_back(element->Nodes()[n+1]->Id());
                  }
                  break;
                }
              }
            }
          }
        }
      }
    }
    // affine shear nodes (in analogy to prior DBC nodes
    for(int lid=0; lid<discret_->NumMyRowElements(); lid++)
    {
      DRT::Element* element = discret_->lRowElement(lid);
      if(element->Id() <= basisnodes_ && element->Id() < statmechparams_.get<int>("NUM_EVAL_ELEMENTS", basiselements_))
      {
        int numdof = 3;
        LINALG::SerialDenseMatrix coord(numdof,(int)discret_->lRowElement(lid)->NumNode(), true);
        LINALG::SerialDenseMatrix cut(numdof,(int)discret_->lRowElement(lid)->NumNode()-1,true);
        std::vector<int> doflids;
        GetElementNodeCoords(element, dis, coord, &doflids);
        for(int n=0; n<cut.N(); n++)
        {
          for(int m=n; m<n+2; m++)
          {
            int nodelidm = discret_->NodeRowMap()->LID(element->Nodes()[m]->Id());
            if(nodelidm>-1)
            {
              if(nodedbcstatus.at(nodelidm)==false)
              {
                nodedbcstatus.at(nodelidm) = true;
                affineshearnodes.push_back(discret_->NodeRowMap()->GID(nodelidm));
                (*deltadbc)[doflids.at(numdof*m+displacementdir)]     =  (coord(2,m)/(*periodlength_)[2])*amp*tcincrement;
              }
            }
          }
        }
      }
    }
    dbcnodesets_.push_back(sensornodes);
    dbcnodesets_.push_back(fixednodes);
    // without if-statement, this node set might added again after restart beyond actiontime_->at(2)
    if(time<=actiontime_->at(2))
      dbcnodesets_.push_back(affineshearnodes);
    //std::cout<<"A Proc "<<discret_->Comm().MyPID()<<": sizeosc = "<<dbcnodesets_[0].size()<<std::endl;

    UpdateForceSensors(dbcnodesets_[0], displacementdir);
  }
  else // only ever use the fixed set of nodes after the first time dirichlet values were imposed
  {
    // top plate nodes
    for(int i=0; i<(int)dbcnodesets_[0].size(); i++)
    {
      int nodelid = discret_->NodeRowMap()->LID(dbcnodesets_[0][i]);
      DRT::Node* displacednode = discret_->lRowNode(nodelid);
      std::vector<int> dofnode = discret_->Dof(displacednode);

      double znode = displacednode->X()[2] + (*discol)[discret_->DofColMap()->LID(dofnode[2])];
      double tcincrement = DRT::Problem::Instance()->Curve(curvenumber).f(time) -
                           DRT::Problem::Instance()->Curve(curvenumber).f(time-dt);
      (*deltadbc)[discret_->DofRowMap()->LID(dofnode[displacementdir])] = znode/(*periodlength_)[2]*amp*tcincrement;
    }
    // all other network nodes that undergo affine deformation (hard coded for now)
    if(dbcnodesets_.size()==3)
    {
      for(int i=0; i<(int)dbcnodesets_[2].size(); i++)
      {
        int nodelid = discret_->NodeRowMap()->LID(dbcnodesets_[2][i]);
        DRT::Node* affinenode = discret_->lRowNode(nodelid);
        std::vector<int> dofnode = discret_->Dof(affinenode);

        double znode = affinenode->X()[2] + (*dis)[discret_->DofRowMap()->LID(dofnode[2])];
        double tcincrement = DRT::Problem::Instance()->Curve(curvenumber).f(time) -
                             DRT::Problem::Instance()->Curve(curvenumber).f(time-dt);
        (*deltadbc)[discret_->DofRowMap()->LID(dofnode[displacementdir])] = znode/(*periodlength_)[2]*amp*tcincrement;
      }
    }
  }
  return;
}

/*----------------------------------------------------------------------*
 |  evaluate Dirichlet conditions (private), mwgee 01/07, mueller 05/12 |
 *----------------------------------------------------------------------*/
void STATMECH::StatMechManager::DoDirichletConditionPredefined(DRT::Condition&             cond,
                                                               DRT::Discretization&        dis,
                                                               const bool                  usetime,
                                                               const double                time,
                                                               Teuchos::RCP<Epetra_Vector> systemvector,
                                                               Teuchos::RCP<Epetra_Vector> systemvectord,
                                                               Teuchos::RCP<Epetra_Vector> systemvectordd,
                                                               Teuchos::RCP<Epetra_Vector> toggle,
                                                               Teuchos::RCP<std::set<int> > dbcgids)
{
  const std::vector<int>* nodeids = cond.Nodes();
  if (!nodeids) dserror("Dirichlet condition does not have nodal cloud");
  const int nnode = (*nodeids).size();
  const std::vector<int>*    curve  = cond.Get<std::vector<int> >("curve");
  const std::vector<int>*    funct  = cond.Get<std::vector<int> >("funct");
  const std::vector<int>*    onoff  = cond.Get<std::vector<int> >("onoff");
  const std::vector<double>* val    = cond.Get<std::vector<double> >("val");

  // determine highest degree of time derivative
  // and first existent system vector to apply DBC to
  unsigned deg = 0;  // highest degree of requested time derivative
  Teuchos::RCP<Epetra_Vector> systemvectoraux = Teuchos::null;  // auxiliar system vector
  if (systemvector != Teuchos::null)
  {
    deg = 0;
    systemvectoraux = systemvector;
  }
  if (systemvectord != Teuchos::null)
  {
    deg = 1;
    if (systemvectoraux == Teuchos::null)
      systemvectoraux = systemvectord;
  }
  if (systemvectordd != Teuchos::null)
  {
    deg = 2;
    if (systemvectoraux == Teuchos::null)
      systemvectoraux = systemvectordd;
  }
  dsassert(systemvectoraux!=Teuchos::null, "At least one vector must be unequal to null");

  // loop nodes to identify and evaluate load curves and spatial distributions
  // of Dirichlet boundary conditions
  for (int i=0; i<nnode; ++i)
  {
    // do only nodes in my row map
    int nlid = dis.NodeRowMap()->LID((*nodeids)[i]);
    if (nlid < 0) continue;
    DRT::Node* actnode = dis.lRowNode( nlid );

    // call explicitly the main dofset, i.e. the first column
    std::vector<int> dofs = dis.Dof(0,actnode);
    const unsigned total_numdf = dofs.size();

    // Get native number of dofs at this node. There might be multiple dofsets
    // (in xfem cases), thus the size of the dofs vector might be a multiple
    // of this.
    const int numele = actnode->NumElement();
    const DRT::Element * const * myele = actnode->Elements();
    int numdf = 0;
    for (int j=0; j<numele; ++j)
      numdf = std::max(numdf,myele[j]->NumDofPerNode(0,*actnode));

    if ( ( total_numdf % numdf ) != 0 )
      dserror( "illegal dof set number" );

    for (unsigned j=0; j<total_numdf; ++j)
    {
      int onesetj = j % numdf;
      if ((*onoff)[onesetj]==0)
      {
        const int lid = (*systemvectoraux).Map().LID(dofs[j]);
        if (lid<0) dserror("Global id %d not on this proc in system vector",dofs[j]);
        if (toggle!=Teuchos::null)
          (*toggle)[lid] = 0.0;
        // get rid of entry in DBC map - if it exists
        if (dbcgids != Teuchos::null)
          (*dbcgids).erase(dofs[j]);
        continue;
      }
      const int gid = dofs[j];
      std::vector<double> value(deg+1,(*val)[onesetj]);

      // factor given by time curve
      std::vector<double> curvefac(deg+1, 1.0);
      int curvenum = -1;
      if (curve) curvenum = (*curve)[onesetj];
      if (curvenum>=0 && usetime)
        curvefac = DRT::Problem::Instance()->Curve(curvenum).FctDer(time,deg);
      else
        for (unsigned i=1; i<(deg+1); ++i) curvefac[i] = 0.0;

      // factor given by spatial function
      double functfac = 1.0;
      int funct_num = -1;
      if (funct)
      {
        funct_num = (*funct)[onesetj];
        if (funct_num>0)
          functfac =
            DRT::Problem::Instance()->Funct(funct_num-1).Evaluate(onesetj,
                                                                  actnode->X(),
                                                                  time,
                                                                  &dis);
      }

      // apply factors to Dirichlet value
      for (unsigned i=0; i<deg+1; ++i)
      {
        value[i] *= functfac * curvefac[i];
      }

      // assign value
      const int lid = (*systemvectoraux).Map().LID(gid);
      if (lid<0) dserror("Global id %d not on this proc in system vector",gid);
      if (systemvector != Teuchos::null)
        (*systemvector)[lid] = value[0];
      if (systemvectord != Teuchos::null)
        (*systemvectord)[lid] = value[1];
      if (systemvectordd != Teuchos::null)
        (*systemvectordd)[lid] = value[2];
      // set toggle vector
      if (toggle != Teuchos::null)
        (*toggle)[lid] = 1.0;
      // amend vector of DOF-IDs which are Dirichlet BCs
      if (dbcgids != Teuchos::null)
        (*dbcgids).insert(gid);
    }  // loop over nodal DOFs
  }  // loop over nodes
  return;
}

/*----------------------------------------------------------------------*
 |  fill system vector and toggle vector (private)         mueller  3/10|
 *----------------------------------------------------------------------*/
void STATMECH::StatMechManager::DoDirichletConditionPeriodic(std::vector<int>*           nodeids,
                                                             std::vector<int>*           onoff,
                                                             Teuchos::RCP<Epetra_Vector> dis,
                                                             Teuchos::RCP<Epetra_Vector> deltadbc,
                                                             Teuchos::RCP<std::set<int> > dbcgids)
/*
 * This basically does the same thing as DoDirichletCondition() (to be found in drt_discret_evaluate.cpp),
 * but with the slight difference of taking current displacements into account.
 * Time curve values aren't added to the reference position(s) of the discretization as usual,
 * but to the latest known 0-position(s). These positions are calculated using the deltadbc
 * vector holding the latest incremental Dirichlet displacement. It is added to the displacement
 * at the end of the preceding time step.
 */
{
  /*/ "condition output"
  std::cout<<"Node Ids: ";
  for(int i=0; i<(int)nodeids->size(); i++)
    std::cout<<nodeids->at(i)<<" ";
  std::cout<<"onoff: ";
  for(int i=0; i<(int)discret_->Dof(0,discret_->gNode(nodeids->at(0))).size(); i++)
    std::cout<<onoff->at(i)<<" ";
  std::cout<<std::endl;*/

  // some checks for errors
  if (!nodeids) dserror("No Node IDs were handed over!");

  // get the condition properties
  const int nnode = nodeids->size();

  // loop over all nodes in condition
  for (int i=0; i<nnode; ++i)
  {
    // do only nodes in my row map
    if (!discret_->NodeRowMap()->MyGID(nodeids->at(i))) continue;
    DRT::Node* actnode = discret_->gNode(nodeids->at(i));
    if (!actnode) dserror("Cannot find global node %d",nodeids->at(i));

    // call explicitly the main dofset, i.e. the first column
    std::vector<int> dofs = discret_->Dof(0,actnode);
    const unsigned numdf = dofs.size();
    // loop over DOFs
    for (unsigned j=0; j<numdf; ++j)
    {
      // get the GID and the LID of the currently handled DOF
      const int gid = dofs[j];
      const int lid = (*dis).Map().LID(gid);

//---------------------------------------------Dirichlet Value Assignment
      if (onoff->at(j)==1)
      {
        // assign value
        if (lid<0) dserror("Global id %d not on this proc in system vector", dofs[j]);
        if (dis != Teuchos::null)
          (*dis)[lid] += (*deltadbc)[lid];
        // amend vector of DOF-IDs which are Dirichlet BCs
        if (dbcgids != Teuchos::null)
          (*dbcgids).insert(gid);
      }
      else // if DOF in question is not subject to DBCs (anymore)
      {
        if (lid<0)
          dserror("Global id %d not on this proc in system vector",dofs[j]);
        // get rid of entry in DBC map - if it exists
        if (dbcgids != Teuchos::null)
          (*dbcgids).erase(dofs[j]);
      }
    }  // loop over nodal DOFs
  }  // loop over nodes
  return;
}

/*----------------------------------------------------------------------*
 | Get DBCs defined in Input file and ad the DBCs DOFs to dbcgids       |            |
 |                                           (private)     mueller 05/12|
 *----------------------------------------------------------------------*/
void STATMECH::StatMechManager::DBCGetPredefinedConditions(Teuchos::ParameterList&      params,
                                                           Teuchos::RCP<Epetra_Vector>  dis,
                                                           Teuchos::RCP<std::set<int> > dbcgids)
{
  //get absolute time
  double time = params.get<double>("total time", 0.0);
  bool usetime = true;
  if (time<0.0) usetime = false;

  // get the Dirichlet conditions given in the input file
  std::vector<DRT::Condition*> dirichlet;
  discret_->GetCondition("Dirichlet",dirichlet);

  for (int i=0; i<(int)dirichlet.size(); i++)
  {
    if (dirichlet[i]->Type() != DRT::Condition::LineDirichlet) continue;
    DoDirichletConditionPredefined(*dirichlet[i],*discret_,usetime,time,dis,Teuchos::null,Teuchos::null,Teuchos::null,dbcgids);
  }
  for (int i=0; i<(int)dirichlet.size(); i++)
  {
    if (dirichlet[i]->Type() != DRT::Condition::PointDirichlet) continue;
    DoDirichletConditionPredefined(*dirichlet[i],*discret_,usetime,time,dis,Teuchos::null,Teuchos::null,Teuchos::null,dbcgids);
  }

  return;
}

/*----------------------------------------------------------------------*
 |  Set Dirichlet Values                 (private)         mueller  5/12|
 *----------------------------------------------------------------------*/
void STATMECH::StatMechManager::DBCSetValues(Teuchos::RCP<Epetra_Vector> dis,
                                             Teuchos::RCP<Epetra_Vector> deltadbc,
                                             Teuchos::RCP<std::set<int> > dbcgids)
{
  // number of DOFs per node and DOF toggle vectors
  int numdof = (int)discret_->Dof(discret_->gNode(discret_->NodeRowMap()->GID(0))).size();
  // DOF toggle vector
  std::vector<int> onoff(numdof,0);

  // manipulation of onoff vectors according to the different DBC types
  INPAR::STATMECH::DBCType dbctype = DRT::INPUT::IntegralValue<INPAR::STATMECH::DBCType>(statmechparams_,"DBCTYPE");
  switch(dbctype)
  {
    case INPAR::STATMECH::dbctype_shearfixed:
      // inhibit translational degrees of freedom
      for(int i=0; i<3; i++)
        onoff.at(i) = 1;
    break;
    case INPAR::STATMECH::dbctype_shearfixeddel:
      // inhibit translational degrees of freedom
      for(int i=0; i<3; i++)
        onoff.at(i) = 1;
    break;
    case INPAR::STATMECH::dbctype_sheartrans:
    {
      int dbcdispdir = statmechparams_.get<int>("DBCDISPDIR",-1)-1;
      onoff.at(dbcdispdir) = 1;
    }
    break;
    case INPAR::STATMECH::dbctype_affineshear:
      // inhibit translational degrees of freedom
      for(int i=0; i<3; i++)
        onoff.at(i) = 1;
    break;
    case INPAR::STATMECH::dbctype_affinesheardel:
      // inhibit translational degrees of freedom
      for(int i=0; i<3; i++)
        onoff.at(i) = 1;
    break;
    case INPAR::STATMECH::dbctype_pinnodes:
      // DOF toggle vector for end nodes of the horizontal filament
      for(int i=0; i<3; i++)
        onoff.at(i) = 1;
    break;
    default:
      dserror("DBC values cannot be imposed for this DBC type: %d", dbctype);
    break;
  }

  // Apply DBC values
  for(int i=0; i<(int)dbcnodesets_.size(); i++)
    if(!dbcnodesets_[i].empty())
      DoDirichletConditionPeriodic(&dbcnodesets_[i], &onoff, dis, deltadbc, dbcgids);

  return;
}

/*----------------------------------------------------------------------*
 | Create Dirichlet DOF maps                    (private)  mueller 05/12|
 *----------------------------------------------------------------------*/
void STATMECH::StatMechManager::DBCCreateMap(Teuchos::RCP<std::set<int> > dbcgids,
                                             Teuchos::RCP<LINALG::MapExtractor> dbcmapextractor)
{
  if (dbcmapextractor != Teuchos::null)
  {
    // build map of Dirichlet DOFs
    int nummyelements = 0;
    int* myglobalelements = NULL;
    std::vector<int> dbcgidsv;
    if (dbcgids->size() > 0)
    {
      dbcgidsv.reserve(dbcgids->size());
      dbcgidsv.assign(dbcgids->begin(),dbcgids->end());
      nummyelements = dbcgidsv.size();
      myglobalelements = &(dbcgidsv[0]);
    }
    Teuchos::RCP<Epetra_Map> dbcmap = Teuchos::rcp(new Epetra_Map(-1, nummyelements, myglobalelements, discret_->DofRowMap()->IndexBase(), discret_->DofRowMap()->Comm()));
    // build the map extractor of Dirichlet-conditioned and free DOFs
    *dbcmapextractor = LINALG::MapExtractor(*(discret_->DofRowMap()), dbcmap);
  }
  return;
}

