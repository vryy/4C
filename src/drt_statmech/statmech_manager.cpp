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
#include "../drt_truss3/truss3.H"
#include "../drt_torsion3/torsion3.H"

#include <iostream>
#include <iomanip>
#include <cstdio>
#include <math.h>

//MEASURETIME activates measurement of computation time for certain parts of the code
//#define MEASURETIME

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
useinitdbcset_(false)
{
  Teuchos::ParameterList parameters = DRT::Problem::Instance()->StructuralDynamicParams();

  //initialize random generators
  SeedRandomGenerators(0);

  // retrieve the dimensions of the periodic boundary box and
  // set spatial resolution for search algorithm binding spots x crosslinkers
  InitializeStatMechValues();

  /*setting and deleting dynamic crosslinkers needs a special code structure for parallel search trees
   * here we provide a fully overlapping column map which is required by the search tree to look for each
   * node for all the neighbouring nodes within a certain distance; note: we pass this fully overlapping
   * map to the discretization and adapt it accordingly; this is not efficient in parallel use, however, it
   * makes sure that it leads to at least correct results; as a consequence this implementation may be con-
   * sidered an implementation capable for parallel use, but not optimized for it; note: the discretization's
   * column map is turned into a fully overlapping one even in case that dynamic crosslinkers are deactivated;
   * the reason is that also the method for Gmsh output currently relies on a fully overlapping column map
   * in case of parallel computing (otherwise the output is not written correctly*/

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

  /* and filamentnumber_ is generated based on a column map vector as each node has to
   * know about each other node its filament number in order to decide weather a crosslink may be established
   * or not; vectors is initalized with -1, which state is changed if filament numbering is used, only*/
  filamentnumber_ = Teuchos::rcp( new Epetra_Vector(*(discret_->NodeColMap())) );
  filamentnumber_->PutScalar(-1);


  /*force sensors can be applied at any degree of freedom of the discretization the list of force sensors should
   * be based on a column map vector so that each processor has not only the information about each node's
   * displacement, but also about whether this has a force sensor; as a consequence each processor can write the
   * complete information gathered by all force sensors into a file of its own without any additional communication
   * with any other processor; initialization with -1 indicates that so far no forcesensors have been set*/
  forcesensor_ = Teuchos::rcp( new Epetra_Vector(*(discret_->DofColMap())) );
  forcesensor_->PutScalar(-1.0);

  // initial clearing of dbc management vectors
  dbcnodesets_.clear();

  /*since crosslinkers should be established only between different filaments the number of the filament
   * each node is belonging to is stored in the condition FilamentNumber; if no such conditions have been defined
   * the default -1 value is upkept in filamentnumber_ and crosslinkers between nodes belonging to the same filament
   * are allowed*/

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

  /*Young's modulus and loss modulus are to be measured by means of the reaction forces at certain sensor points; example: if for a
   * an actin network between two rheometer plates the stiffness is to be determined this can be done by measuring the forces exerted
   * to the upper plate which is moving forwards and backwards for example; so for measurements of viscoelastic properties of materials
   * whithin a system certain sensor points have to be specified and in order to handle this within BACI these points are marked by
   * means of the condition sensorcondition*/

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

  /*after having generated a search tree and a discretization with fully overlapping column map we initialize the search tree
   * for accelerated search for each nodes neighbouring nodes; note: the tree is based on a bounding box
   * with respect to the reference positions (which are the positions at the beginning of the simulation;
   * in case of large overall deformations of the fiber network such an initialization would have to be carried
   * out in each time step with respect to the current node positions*/

  /*currenpositions is a map which relates each LID of any node on this processor to the nodes coordinates*/
  std::map<int,LINALG::Matrix<3,1> > currentpositions;

  currentpositions.clear();

  for (int lid = 0; lid <discret_->NumMyColNodes(); ++lid)
  {
    const DRT::Node* node = discret_->lColNode(lid);
    LINALG::Matrix<3,1> currpos;
    currpos(0) = node->X()[0];
    currpos(1) = node->X()[1];
    currpos(2) = node->X()[2];
    currentpositions[node->LID()] = currpos;
  }

  /* Initialization of N_CROSSLINK crosslinker molecule REPRESENTATIONS. As long as the molecules do not act as a link
   * between two filaments, only their positions are calculated. Here, the molecules' initial positions are determined.
   * Calculations are made on Proc 0, only.*/
  CrosslinkerMoleculeInit();

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

//  if(!discret_->Comm().MyPID())
//  {
//    for(int i=0; i<(int)periodlength_->size(); i++)
//      cout<<periodlength_->at(i)<<" ";
//    cout<<endl;
//  }

  // set spatial resolution for search algorithm binding spots x crosslinkers
  if(statmechparams_.get<int>("SEARCHRES",1)<1)
    dserror("Please give a plausible value for SEARCHRES!");

  // determine search resolution in each spatial direction according to the periodlength_ vector
  searchres_ = Teuchos::rcp(new std::vector<int>(3,statmechparams_.get<int>("SEARCHRES",1)));
  // in case of a non-cubic periodic volume
  if(fabs(pow(periodlength_->at(0)*periodlength_->at(1)*periodlength_->at(2), 1.0/3.0)-periodlength_->at(0))>1e-4)
  {
    double Hmax = max(periodlength_->at(0), periodlength_->at(1));
    Hmax = max(Hmax, periodlength_->at(2));
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
  for(int i=0; i<(int)actiontime_->size()-1; i++)
    if(actiontime_->at(i)>actiontime_->at(i+1))
      dserror("ACTIONTIME values must be monotonously increasing!");
  if((int)actiontime_->size()!=(int)timestepsizes_->size())
    dserror("ACTIONTIME and ACTIONDT have to be equal in number in the input file!");
  if((int)actiontime_->size()<2)
    dserror("ACTIONTIME has to have at least 2 values");

  // increase the vector position if time values are equal
  for(int i=0; i<(int)actiontime_->size()-1; i++)
    if(fabs(actiontime_->at(i)-actiontime_->at(i+1))<timestepsizes_->at(i)/1e3)
      timeintervalstep_++;

//  if(!discret_->Comm().MyPID())
//  {
//    for(int i=0; i<(int)searchres_->size(); i++)
//      cout<<searchres_->at(i)<<" ";
//    cout<<endl;
//  }
  return;
}

/*----------------------------------------------------------------------*
 | write special output for statistical mechanics (public)    cyron 09/08|
 *----------------------------------------------------------------------*/
void STATMECH::StatMechManager::Update(const int& istep,
                                       const double& timen,
                                       const double& dt,
                                       Epetra_Vector& disrow,
                                       Teuchos::RCP<LINALG::SparseOperator>& stiff,
                                       int& ndim, Teuchos::RCP<CONTACT::Beam3cmanager> beamcmanager,
                                       bool rebuildoctree,
                                       bool printscreen)
{
#ifdef MEASURETIME
  const double t_start = Teuchos::Time::wallTime();
#endif // #ifdef MEASURETIME
  /* first we modify the displacement vector so that current nodal position at the end of current time step complies with
   * periodic boundary conditions, i.e. no node lies outside a cube of edge length periodlength_*/

  //if dynamic crosslinkers are used update comprises adding and deleting crosslinkers
  if (DRT::INPUT::IntegralValue<int>(statmechparams_, "DYN_CROSSLINKERS"))
  {
    // crosslink molecule diffusion
    double standarddev = sqrt(statmechparams_.get<double> ("KT", 0.0) / (2*M_PI * statmechparams_.get<double> ("ETA", 0.0) * statmechparams_.get<double> ("R_LINK", 0.0)) * dt);
    CrosslinkerDiffusion(disrow, 0.0, standarddev, dt);

    /*the following tow rcp pointers are auxiliary variables which are needed in order provide in the very end of the
     * crosslinker administration a node row and column map; these maps have to be taken here before the first modification
     * by deleting and adding elements have been carried out with the discretization since after such modifications the maps
     * cannot be called from the discretization before calling FillComplete() again which should be done only in the very end
     * note: this way of getting node maps is all right only if no nodes are added/ deleted, but elements only*/
    const Epetra_Map noderowmap = *discret_->NodeRowMap();
    const Epetra_Map nodecolmap = *discret_->NodeColMap();

    //node positions and rotations (binding spot positions and rotations when applying 4-noded beam element)
    std::map<int, LINALG::Matrix<3, 1> > currentpositions;
    std::map<int, LINALG::Matrix<3, 1> > currentrotations;

    /*note: access by ExtractMyValues requires column map vector, whereas displacements on level of time integration are
     * handled as row map vector*/
    Epetra_Vector discol(*discret_->DofColMap(), true);
    LINALG::Export(disrow, discol);
    GetNodePositions(discol, currentpositions, currentrotations);

    // set crosslinkers, i.e. considering crosslink molecule diffusion after filaments had time to equilibrate
    if(timen>=actiontime_->front() || fabs(timen-actiontime_->front())<(dt/1e3))
    {
      if(beamcmanager!=Teuchos::null && rebuildoctree)
        beamcmanager->OcTree()->OctTreeSearch(currentpositions);

      SearchAndSetCrosslinkers(istep, timen, dt, noderowmap, nodecolmap, currentpositions,currentrotations,beamcmanager,printscreen);
      
      if(beamcmanager!=Teuchos::null && rebuildoctree)
        beamcmanager->ResetPairs();
      
      SearchAndDeleteCrosslinkers(timen, dt, noderowmap, nodecolmap, currentpositions, discol,beamcmanager,printscreen);
    }

    // reset thermal energy to new value (simple value change for now, maybe Load Curve later on)
    if(fabs(timen-actiontime_->at(1))<(dt/1e3))
      statmechparams_.set("KT",statmechparams_.get<double>("KTACT",statmechparams_.get<double>("KT",0.0)));

    // // Set a certain number of double-bonded crosslinkers free
    if(fabs(timen-dt-actiontime_->back())<(dt/1e3) && statmechparams_.get<int>("REDUCECROSSLINKSBY",0)>0)
      ReduceNumOfCrosslinkersBy(statmechparams_.get<int>("REDUCECROSSLINKSBY",0));

    // force dependent unlinking: store displacement vector of current time step for next one
    if(DRT::INPUT::IntegralValue<int>(statmechparams_, "FORCEDEPUNLINKING"))
      for(int i=0; i<disprev_->MyLength(); i++)
        (*disprev_)[i] = discol[i];

    /*settling administrative stuff in order to make the discretization ready for the next time step:
     * done in SearchAndeDeleteCrosslinkers():
     * synchronize the Filled() state on all processors after having added or deleted elements by CheckFilledGlobally();
     * then build new element maps and call FillComplete();
     * done here: finally Crs matrices stiff_ has to be deleted completely and made ready
     * for new assembly since their graph was changed*/
    stiff->Reset();

#ifdef MEASURETIME
    cout << "\n***\nadministration time: " << Teuchos::Time::wallTime() - t_admin<< " seconds\n***\n";
#endif // #ifdef MEASURETIME
  }//if(DRT::INPUT::IntegralValue<int>(statmechparams_,"DYN_CROSSLINKERS"))

#ifdef MEASURETIME
  const double Delta_t = Teuchos::Time::wallTime()-t_start;
  cout << "\n***\ntotal time: " << Delta_t<< " seconds\n***\n";
#endif // #ifdef MEASURETIME

  return;
} // StatMechManager::Update()

/*----------------------------------------------------------------------*
 | Update time step size in time integration       (public)mueller 06/12|
 *----------------------------------------------------------------------*/
void STATMECH::StatMechManager::UpdateTimeAndStepSize(double& dt, double& timeconverged)
{
  // update time step
  double dtnew = timestepsizes_->at(timeintervalstep_);
  // update step size
  double nexttimethreshold = actiontime_->at(timeintervalstep_);
  double eps = 1.0e-10;
  if((timeconverged>=nexttimethreshold || fabs(timeconverged-nexttimethreshold)<eps) && dtnew>0.0)
  {
    if(!discret_->Comm().MyPID())
      cout<<"dtnew = "<<timestepsizes_->at(timeintervalstep_)<<", nexttimethreshold = "<<nexttimethreshold<<endl;
    dt = dtnew;
    timeintervalstep_++;
  }
  return;
}

/*----------------------------------------------------------------------*
 | Retrieve nodal positions                     (public)   mueller 09/11|
 *----------------------------------------------------------------------*/
void STATMECH::StatMechManager::GetNodePositions(Epetra_Vector& discol,
                                       std::map<int,LINALG::Matrix<3,1> >& currentpositions,
                                       std::map<int,LINALG::Matrix<3,1> >& currentrotations,
                                       bool positionsonly)
{
  /*in preparation for later decision whether a crosslink should be established between two nodes (binding spots) we first store the
   * current positions of all column map nodes (column map binding spots) in the map currentpositions; additionally we store the rotational displacements
   * analogously in a map currentrotations for later use in setting up reference geometry of crosslinkers; the maps
   * currentpositions and currentrotations relate positions and rotations to a local column map node Id, respectively*/
  currentpositions.clear();
  if(!positionsonly)
    currentrotations.clear();

  // in case the four-noded crosslinker beam element is used, currentpositions and currentrotations have to be set up another way
  if(DRT::INPUT::IntegralValue<int>(statmechparams_, "INTERNODALBSPOTS"))
    UpdateBindingSpots(discol, currentpositions);
  else	// conventional crosslinker beam element, i.e. binding spots coincide with nodes
  {
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
  }
  return;
}
/*----------------------------------------------------------------------*
 | Updates positions and rotations of binding spots        mueller 11/11|
 *----------------------------------------------------------------------*/
void STATMECH::StatMechManager::UpdateBindingSpots(const Epetra_Vector& discol,std::map<int,LINALG::Matrix<3,1> >& currentpositions)
{
  Epetra_Vector bspotpositionsrow(*bspotrowmap_);
  Epetra_Vector bspotrotationsrow(*bspotrowmap_);
  Epetra_Import bspotimporter(*bspotcolmap_, *bspotrowmap_);

  DRT::Element* filelement = NULL;
	DRT::Node* node0 = NULL;
	DRT::Node* node1 = NULL;

	//double currelelength = 0.0;

	// get current positions and rotations of this element's nodes
	std::vector<int> dofnode0;
	std::vector<int> dofnode1;

	std::vector<double> currpos0(3,0.0);
	std::vector<double> currpos1(3,0.0);
	std::vector<double> currrot0(3,0.0);
	std::vector<double> currrot1(3,0.0);

  // loop over row binding spots
  int prevelegid = -1;
  for(int i=0; i<bspotrowmap_->NumMyElements(); i++)
  {
    // retrieve the element the binding spot belongs to and retrieve its nodes
    int elegid = (int)(*bspot2element_)[i];
    // only recalculate nodal positions and rotations if the element GID changed compared to the previous binding spot
    if(elegid!=prevelegid)
    {
      filelement = discret_->gElement(elegid);
      node0 = discret_->gNode(filelement->NodeIds()[0]);
      node1 = discret_->gNode(filelement->NodeIds()[1]);

      dofnode0 = discret_->Dof(node0);
      dofnode1 = discret_->Dof(node1);

      // current positions of node0 and node1
      currpos0[0] = node0->X()[0] + discol[discret_->DofColMap()->LID(dofnode0[0])];
      currpos0[1] = node0->X()[1] + discol[discret_->DofColMap()->LID(dofnode0[1])];
      currpos0[2] = node0->X()[2] + discol[discret_->DofColMap()->LID(dofnode0[2])];
      currpos1[0] = node0->X()[0] + discol[discret_->DofColMap()->LID(dofnode1[0])];
      currpos1[1] = node0->X()[1] + discol[discret_->DofColMap()->LID(dofnode1[1])];
      currpos1[2] = node0->X()[2] + discol[discret_->DofColMap()->LID(dofnode1[2])];
      if (discret_->NumDof(node0) == 6)
      {
        currrot0[0] = discol[discret_->DofColMap()->LID(dofnode0[3])];
        currrot0[1] = discol[discret_->DofColMap()->LID(dofnode0[4])];
        currrot0[2] = discol[discret_->DofColMap()->LID(dofnode0[5])];
        currrot1[0] = discol[discret_->DofColMap()->LID(dofnode1[3])];
        currrot1[1] = discol[discret_->DofColMap()->LID(dofnode1[4])];
        currrot1[2] = discol[discret_->DofColMap()->LID(dofnode1[5])];
      }

      // update element length (i.e. distance between the nodes)

      prevelegid = elegid;
    }

    // interpolation of positions and rotations...soon to follow
  }
}

/*----------------------------------------------------------------------*
 | Shifts current position of nodes so that they comply with periodic   |
 | boundary conditions                                       cyron 04/10|
 *----------------------------------------------------------------------*/
void STATMECH::StatMechManager::PeriodicBoundaryShift(Epetra_Vector& disrow, int ndim, const double &timen, const double &dt)
{
  double starttime = actiontime_->at((int)(actiontime_->size()-1));
  double shearamplitude = statmechparams_.get<double> ("SHEARAMPLITUDE", 0.0);
  int curvenumber = statmechparams_.get<int> ("CURVENUMBER", -1);
  int oscilldir = statmechparams_.get<int> ("OSCILLDIR", -1);


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
        if (xcurr > periodlength_->at(j))
        {
          disrow[discret_->DofRowMap()->LID(dofnode[j])] -= periodlength_->at(j)*floor(xcurr/periodlength_->at(j));

          /*the upper domain surface orthogonal to the z-direction may be subject to shear Dirichlet boundary condition; the lower surface
           *may fixed by DBC. To avoid problems when nodes exit the domain through the upper z-surface and reenter through the lower
           *z-surface, the shear has to be substracted from nodal coordinates in that case */
          if (j == 2 && curvenumber >= 1 && timen > starttime && fabs(timen-starttime)>dt/1e4)
            disrow[discret_->DofRowMap()->LID(dofnode[oscilldir])] -= shearamplitude * DRT::Problem::Instance()->Curve(curvenumber - 1).f(timen);
        }
        /*if node currently has coordinate value smaller than zero, it is shifted by periodlength sufficiently often
         *to lie again in the domain*/
        if (xcurr < 0.0)
        {
          disrow[discret_->DofRowMap()->LID(dofnode[j])] -= periodlength_->at(j)*floor(xcurr/periodlength_->at(j));

          /*the upper domain surface orthogonal to the z-direction may be subject to shear Dirichlet boundary condition; the lower surface
           *may be fixed by DBC. To avoid problems when nodes exit the domain through the lower z-surface and reenter through the upper
           *z-surface, the shear has to be added to nodal coordinates in that case */
          if (j == 2 && curvenumber >= 1 && timen > starttime && fabs(timen-starttime)>dt/1e4)
            disrow[discret_->DofRowMap()->LID(dofnode[oscilldir])] += shearamplitude * DRT::Problem::Instance()->Curve(curvenumber - 1).f(timen);
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

  //3D beam elements are embeddet into R^3:
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
}

/*------------------------------------------------------------------------*
 | This function loops through all the elements of the discretization and |
 | tests whether beam3 are broken by periodic boundary conditions in the  |
 | reference configuration; if yes initial values of curvature and jacobi |
 | determinants are adapted in a proper way                    cyron 02/10|
 *-----------------------------------------------------------------------*/
void STATMECH::StatMechManager::PeriodicBoundaryBeam3iiInit(DRT::Element* element)
{
  DRT::ELEMENTS::Beam3ii* beam = dynamic_cast<DRT::ELEMENTS::Beam3ii*>(element);

  //3D beam elements are embeddet into R^3:
  const int ndim = 3;

  /*get reference configuration of beam3ii element in proper format for later call of SetUpReferenceGeometry;
   * note that rotrefe for beam3ii elements is related to the entry in the global total Lagrange displacement
   * vector related to a certain rotational degree of freedom; as the displacement is initially zero also
   * rotrefe is set to zero here*/
  std::vector<double> xrefe(beam->NumNode()*ndim,0);
  std::vector<double> rotrefe(beam->NumNode()*ndim,0);

  for(int i=0;i<beam->NumNode();i++)
  for(int dof=0; dof<ndim; dof++)
  {
    xrefe[3*i+dof] = beam->Nodes()[i]->X()[dof];
    rotrefe[3*i+dof] = 0.0;
  }

  /*loop through all nodes except for the first node which remains fixed as reference node; all other nodes are
   * shifted due to periodic boundary conditions if required*/
  for(int i=1;i<beam->NumNode();i++)
  {
    for(int dof=0; dof<ndim; dof++)
    {
      /*if the distance in some coordinate direction between some node and the first node becomes smaller by adding or subtracting
       * the period length, the respective node has obviously been shifted due to periodic boundary conditions and should be shifted
       * back for evaluation of element matrices and vectors; this way of detecting shifted nodes works as long as the element length
       * is smaller than half the periodic length*/
      if( fabs( (beam->Nodes()[i]->X()[dof]) + periodlength_->at(dof) - (beam->Nodes()[0]->X()[dof]) ) < fabs( (beam->Nodes()[i]->X()[dof]) - (beam->Nodes()[0]->X()[dof]) ) )
        xrefe[3*i+dof] += periodlength_->at(dof);

      if( fabs( (beam->Nodes()[i]->X()[dof]) - periodlength_->at(dof) - (beam->Nodes()[0]->X()[dof]) ) < fabs( (beam->Nodes()[i]->X()[dof]) - (beam->Nodes()[0]->X()[dof]) ) )
        xrefe[3*i+dof] -= periodlength_->at(dof);
    }
  }

  /*SetUpReferenceGeometry is a templated function; note that the third argument "true" is necessary as all beam elements
   * have already been initialized once upon reading input file*/
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
 |																								(public) mueller 08/10|
 *----------------------------------------------------------------------*/
void STATMECH::StatMechManager::PartitioningAndSearch(const std::map<int,LINALG::Matrix<3,1> >& currentpositions, Epetra_MultiVector& bspottriadscol, Teuchos::RCP<Epetra_MultiVector>& neighbourslid)
{
  //filament binding spots and crosslink molecules are indexed according to their positions within the boundary box volume
  std::vector<std::vector<std::vector<int> > > bspotinpartition;
  for(int i=0; i<(int)searchres_->size(); i++)
    bspotinpartition.push_back(std::vector<std::vector<int> >(searchres_->at(i), std::vector<int>()));

  /*nodes*/
  // loop over node positions to map their column map LIDs to partitions
  for (std::map<int, LINALG::Matrix<3, 1> >::const_iterator posi = currentpositions.begin(); posi != currentpositions.end(); posi++)
  {
    for(int j=0; j<(int)bspotinpartition.size(); j++) // bspotinpartition.size==3
    {
      int partition = (int)std::floor((posi->second)(j)/periodlength_->at(j)*(double)(searchres_->at(j)));
      if(partition==(int)searchres_->at(j))
        partition--;
      bspotinpartition[j][partition].push_back((int)(posi->first)); //column lid
    }
  }

  /*crosslink molecules*/
  // Export crosslinkerpositions_ to transfermap_ format (kind of a row map format for crosslink molecules)
  Epetra_MultiVector crosslinkerpositionstrans(*transfermap_, 3, true);
  Epetra_Vector numbondtrans(*transfermap_, true);
  Epetra_MultiVector crosslinkpartitiontrans(*transfermap_, 3, false);

  CommunicateVector(numbondtrans, *numbond_, true, false);
  CommunicateMultiVector(crosslinkerpositionstrans, *crosslinkerpositions_, true, false);

  for(int i=0; i<crosslinkpartitiontrans.MyLength(); i++)
  {
    // mark entries with double-bonded crosslink molecules
    if(numbondtrans[i]>1.9)
    {
      for(int j=0; j<crosslinkpartitiontrans.NumVectors(); j++)
        crosslinkpartitiontrans[j][i] = -1.0;
      continue;
    }
    else
    {
      for(int j=0; j<crosslinkpartitiontrans.NumVectors(); j++)
      {
        int partition = (int)std::floor(crosslinkerpositionstrans[j][i]/periodlength_->at(j)*(double)(searchres_->at(j)));
        if(partition==searchres_->at(j))
          partition--;
        crosslinkpartitiontrans[j][i] = partition;
      }
    }
  }
  // detection of nodes within search proximity of the crosslink molecules
  DetectNeighbourNodes(currentpositions, &bspotinpartition, numbondtrans, crosslinkerpositionstrans, crosslinkpartitiontrans, bspottriadscol, neighbourslid);
  return;
}//void StatMechManager::PartitioningAndSearch

/*----------------------------------------------------------------------*
 | detect neighbour nodes to crosslink molecules (public) mueller (7/10)|
 *----------------------------------------------------------------------*/
void STATMECH::StatMechManager::DetectNeighbourNodes(const std::map<int,LINALG::Matrix<3,1> >& currentpositions,
                                           std::vector<std::vector<std::vector<int> > >* bspotinpartition,
                                           Epetra_Vector& numbond,
                                           Epetra_MultiVector& crosslinkerpositions,
                                           Epetra_MultiVector& crosslinkpartitions,
                                           Epetra_MultiVector& bspottriadscol,
                                           Teuchos::RCP<Epetra_MultiVector>& neighbourslid)
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
  Epetra_MultiVector crosslinkerbond(*transfermap_, 2, true);
  CommunicateMultiVector(crosslinkerbond, *crosslinkerbond_);

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
      if ((int)numbond[part]>0.9 && !DRT::INPUT::IntegralValue<int>(statmechparams_,"HELICALBINDINGSTRUCT"))
      {
        rmin *= 2.0;
        rmax *= 2.0;
      }

      // first component
      for(int ilayer=(int)crosslinkpartitions[0][part]-1; ilayer<(int)crosslinkpartitions[0][part]+2; ilayer++)
        if(ilayer>-1 && ilayer<searchres_->at(0))
          for(int i=0; i<(int)(*bspotinpartition)[0][ilayer].size(); i++)
          {
            int tmplid = (int)(*bspotinpartition)[0][ilayer][i];
            // second component
            for(int jlayer=(int)crosslinkpartitions[1][part]-1; jlayer<(int)crosslinkpartitions[1][part]+2; jlayer++)
              if(jlayer>-1 && jlayer<searchres_->at(1))
                for(int j=0; j<(int)(*bspotinpartition)[1][jlayer].size(); j++)
                  if((*bspotinpartition)[1][jlayer][j]==tmplid)
                  {
                    //third component
                    for(int klayer=(int)crosslinkpartitions[2][part]-1; klayer<(int)crosslinkpartitions[2][part]+2; klayer++)
                      if(klayer>-1 && klayer<searchres_->at(2))
                        for(int k=0; k<(int)(*bspotinpartition)[2][klayer].size(); k++)
                          if((*bspotinpartition)[2][klayer][k]==tmplid)
                          {
                            // get the current node position for the node with LID==tmplid
                            const map<int, LINALG::Matrix<3, 1> >::const_iterator nodepos = currentpositions.find(tmplid);
                            // calculate distance crosslinker-node
                            LINALG::Matrix<3, 1> difference;
                            for (int l=0; l<(int)difference.M(); l++)
                              difference(l) = crosslinkerpositions[l][part]-(nodepos->second)(l);
                            // 1. criterion: distance between linker and binding spot within given interval
                            if(difference.Norm2()<rmax && difference.Norm2()>rmin)
                            {
                              // further calculations in case of helical binding spot geometry and singly bound crosslinkers
                              // note: a singly-bound crosslinker is already oriented if helicality is assumed. Thus below.
                              // Free linkers are much more flexible in their choice of a binding spot since they can translate and
                              // rotate freely. Thus, even if helicality is assumed, the simple proximity criterion holds.
                              if(DRT::INPUT::IntegralValue<int>(statmechparams_, "HELICALBINDINGSTRUCT") && numbond[part]>0.9)
                              {
                                // 2. criterion: linker has to lie in the oriented cone with peak "binding spot location"
                                // first and second triad vectors (tangent and normal)
                                LINALG::Matrix<3,1> firstdir;
                                LINALG::Matrix<3,1> bspotvec;
                                // retrieve tangential and normal vector from binding spot quaternions
                                {
                                  LINALG::Matrix<3,3> bspottriad;
                                  // auxiliary variable for storing a triad in quaternion form
                                  LINALG::Matrix<4, 1> qnode;
                                  // triad of node on first filament which is affected by the new crosslinker
                                  for (int l=0; l<4; l++)
                                    qnode(l) = bspottriadscol[l][tmplid];
                                  LARGEROTATIONS::quaterniontotriad(qnode, bspottriad);
                                  for (int l=0; l<(int)bspottriad.M(); l++)
                                  {
                                    firstdir(l) = bspottriad(l,0);
                                    bspotvec(l) = bspottriad(l,1);
                                  }
                                }
                                // rotation matrix around tangential vector by given angle
                                RotationAroundFixedAxis(firstdir,bspotvec,(*bspotorientations_)[tmplid]);
                                // linker position
                                LINALG::Matrix<3,1> crossbspotdiff;
                                for(int l=0; l<(int)crossbspotdiff.M(); l++)
                                  crossbspotdiff(l) = crosslinkerpositions[l][part]-(nodepos->second)(l);
                                // line parameter of the intersection point of the line through the binding spot with the orientation of the binding spot.
                                // a)lambda must be larger than zero in order to lie within the cone
                                double lambda = -(crossbspotdiff.Dot(bspotvec))/(bspotvec.Dot(bspotvec));
                                if(lambda>0)
                                {
                                  difference.Scale(1.0/difference.Norm2());
                                  // b) the angle between the binding spot orientation and the vector between bspot and crosslinker must be smaller than PHIBSPOT
                                  //    Only then does the crosslinker lie within the cone
                                  if(acos(fabs(bspotvec.Dot(difference)))<statmechparams_.get<double>("PHIBSPOT", 0.524))
                                    neighbournodes[part].push_back(tmplid);
                                }
                              }
                              else // only difference criterion applied
                                neighbournodes[part].push_back(tmplid);
                            }
                            // exit loop immediately
                            break;
                          }
                    break;
                  }
          }
      // "-1" indicates the possibility of a crosslink molecule becoming passive, i.e. hypothetically bonding to the same filament
      if((int)numbond[part]==1)
        neighbournodes[part].push_back(-1);
      // store local maximal number of LIDs per molecule in order to determine neighbourslid->NumVectors()
      maxneighbourslocal = max(maxneighbourslocal, (int)neighbournodes[part].size());
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
  CommunicateMultiVector(*neighbourslidtrans, *neighbourslid, false, true);

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
 | crosslinker elements once linking conditions are met.       					|
 | (private)                                              mueller (7/10)|
 *----------------------------------------------------------------------*/
void STATMECH::StatMechManager::SearchAndSetCrosslinkers(const int& istep, const double& timen, const double& dt, const Epetra_Map& noderowmap, const Epetra_Map& nodecolmap,
                                                         const std::map<int, LINALG::Matrix<3, 1> >& currentpositions,
                                                         const std::map<int,LINALG::Matrix<3, 1> >& currentrotations,
                                                         Teuchos::RCP<CONTACT::Beam3cmanager> beamcmanager,
                                                         bool printscreen)
{
  double t_search = Teuchos::Time::wallTime();
  /*preliminaries*/
  // BINDING SPOT TRIAD UPDATE
  Epetra_MultiVector bspottriadscol(nodecolmap,4,true);
  if (DRT::INPUT::IntegralValue<int>(statmechparams_, "CHECKORIENT") || DRT::INPUT::IntegralValue<int>(statmechparams_, "HELICALBINDINGSTRUCT"))
    GetBindingSpotTriads(&bspottriadscol);

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
    PartitioningAndSearch(currentpositions,bspottriadscol, neighbourslid);

  //cout<<"neighbourslid: "<<neighbourslid->MyLength()<<"x"<<neighbourslid->NumVectors()<<endl;
  //cout<<"\n\nneighbourslid\n"<<*neighbourslid<<endl;
  /*the following part of the code is executed on processor 0 only and leads to the decision, at which nodes crosslinks shall be set
   *processor 0 goes through all the crosslink molecules and checks whether a crosslink is to be set; this works precisely as follows:
   *(1) the crosslink molecules are looped through in a random order
   *(2) if a node has not yet reached its maximal number of crosslinks, a crosslink may be set
   *(3) a crosslink is set if and only if the node passes a probability check
   *note: a fully overlapping node map on processor 0 is imperative for the following part of the code to work correctly!*/

  // a vector indicating the crosslink molecule which is going to constitute a crosslinker element
  Epetra_Vector addcrosselement(*crosslinkermap_, true);

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

        // loop over neighbour nodes
        for(int j=0; j<neighbourslid->NumVectors(); j++)
        {
          // random index
          int index = neighbourorder[j];

          // skip this, if neighbourslid entry is '-2', meaning empty
          if((*neighbourslid)[index][irandom] > -1.9)
          {
            // current neighbour LID
            int nodeLID = (int)(*neighbourslid)[index][irandom];

            //skip the rest of this iteration of the i-loop in case the binding spot is occupied
            if(nodeLID>-1)
              if((*bspotstatus_)[nodeLID]>-0.1)
                continue;

            // flag indicating loop break after first new bond has been established between i-th crosslink molecule and j-th neighbour node
            bool bondestablished = false;
            // necessary condition to be fulfilled in order to set a crosslinker
            double probability = 0.0;
            // switch between probability for regular inter-filament binding and self-binding crosslinker
            if(nodeLID>=0)
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
                  if(DRT::INPUT::IntegralValue<int>(statmechparams_, "LOOMSETUP") && (*filamentnumber_)[nodeLID]!=0)
                    break;
                  else // standard case
                  {
                    // update of crosslink molecule positions
                    LINALG::SerialDenseMatrix LID(1,1,true);
                    LID(0,0) = nodeLID;
                    (*crosslinkerbond_)[free][irandom] = nodecolmap.GID(nodeLID);
                    // update status of binding spot by putting in the crosslinker id
                    ((*bspotstatus_)[nodeLID]) = irandom;
                    // increment the number of bonds of this crosslinker
                    ((*numbond_)[irandom]) = 1.0;
                    CrosslinkerIntermediateUpdate(currentpositions, LID, irandom, bspottriadscol);
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
                  LID(0,0) = nodecolmap.LID((int)(*crosslinkerbond_)[occupied][irandom]);
                  // distinguish between a real nodeLID and the entry '-1', which indicates a passive crosslink molecule
                  if(nodeLID>=0)
                    LID(1,0) = nodeLID;
                  else // choose the neighbor node on the same filament as nodeLID as second entry and take basisnodes_ into account
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
                  if((*filamentnumber_)[(int)LID(0,0)]==(*filamentnumber_)[(int)LID(1,0)] && statmechparams_.get<double>("K_ON_SELF",0.0)==0.0)
                    break;
                  // do not do anything in case of a Loom set up if conditions hereafter are not met
                  if(!SetCrosslinkerLoom(LID, currentpositions, bspottriadscol))
                    break;

                  //unit direction vector between currently considered two nodes
                  LINALG::Matrix<3,1> direction((currentpositions.find((int)LID(0,0)))->second);
                  direction -= (currentpositions.find((int)LID(1,0)))->second;
                  direction.Scale(1.0/direction.Norm2());

                  /* In case of beam contact, evaluate whether a crosslinker would intersect with any exisiting elements
                   * if it were to be set between the two nodes considered*/
                  bool intersection = false;
                  if(DRT::INPUT::IntegralValue<int>(statmechparams_,"BEAMCONTACT"))
                  {
                    Epetra_SerialDenseMatrix nodecoords(3,2);
                    for(int k=0; k<nodecoords.M(); k++)
                    {
                      nodecoords(k,0) = ((currentpositions.find((int)LID(0,0)))->second)(k);
                      nodecoords(k,1) = ((currentpositions.find((int)LID(1,0)))->second)(k);
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
                  if(CheckOrientation(direction,bspottriadscol,LID) && !intersection)
                  {
                    numsetelements++;

                    ((*numbond_)[irandom]) = 2.0;
                    // actually existing two bonds
                    if(nodeLID>=0)
                    {
                      (*crosslinkerbond_)[free][irandom] = nodecolmap.GID(nodeLID);
                      ((*bspotstatus_)[nodeLID]) = irandom;
                      // set flag at irandom-th crosslink molecule that an element is to be added
                      addcrosselement[irandom] = 1.0;
                      // update molecule positions
                      CrosslinkerIntermediateUpdate(currentpositions, LID, irandom, bspottriadscol);

                      // consider crosslinkers covering two binding spots of one filament
                      if((*filamentnumber_)[(int)LID(0,0)]==(*filamentnumber_)[nodeLID])
                        (*crosslinkonsamefilament_)[irandom] = 1.0;
                    }
                    else // passive crosslink molecule
                    {
                      (*searchforneighbours_)[irandom] = 0.0;
                      LINALG::SerialDenseMatrix oneLID(1,1,true);
                      oneLID(0,0) = LID(0,0);
                      CrosslinkerIntermediateUpdate(currentpositions, oneLID, irandom, bspottriadscol);
                    }
                    bondestablished = true;
                  }
                }
                break;
              }// switch((int)(*numbond_)[irandom])

              // for now, break j-loop after a new bond was established, i.e crosslinker elements cannot be established starting from zero bonds
              if(bondestablished)
                break;
            }//if( (*uniformgen_)() < probability )
          }// if(neighbourslid...)
        }// for(int j=0; j<(int)neighboursLID.size(); j++)
      }//if((*numbond_)[irandom]<1.9)
    }// for(int i=0; i<numbond_->MyLength(); i++)
    if(printscreen)
      cout << "\nsearch time: " << Teuchos::Time::wallTime() - t_search<< " seconds";
  }// if(discret_->Comm().MypPID==0)

  /* note: searchforneighbours_ and crosslinkonsamefilament_ are not being communicated
   * to the other Procs because their information is of concern to Proc 0 only.*/
  //synchronize information about number of bonded filament nodes by exporting it to row map format and then reimporting it to column map format
  Epetra_Vector bspotstatusrow(*bspotrowmap_,true);
  CommunicateVector(bspotstatusrow, *bspotstatus_);

  // transfer vectors
  Epetra_MultiVector crosslinkerpositionstrans(*transfermap_,3,true);
  Epetra_MultiVector crosslinkerbondtrans(*transfermap_,2,true);
  Epetra_Vector numbondtrans(*transfermap_, true);
  Epetra_Vector addcrosselementtrans(*transfermap_, true);
  // exports and reimports
  CommunicateMultiVector(crosslinkerpositionstrans, *crosslinkerpositions_);
  CommunicateMultiVector(crosslinkerbondtrans, *crosslinkerbond_);
  CommunicateVector(numbondtrans, *numbond_);
  CommunicateVector(addcrosselementtrans, addcrosselement);

  // ADDING ELEMENTS
  // add elements to problem discretization (processor specific)
  for(int i=0; i<addcrosselement.MyLength(); i++)
  {
    if(addcrosselement[i]>0.9)
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

      // obtain node GID
      int nodeGID[2] = {	(int)(*crosslinkerbond_)[0][i],(int)(*crosslinkerbond_)[1][i]};
      // determine smaller and larger of the GIDs
      int GID2 = min(nodeGID[0],nodeGID[1]);
      int GID1 = max(nodeGID[0],nodeGID[1]);
      int globalnodeids[2] = {GID1, GID2};
      // calculate element GID
      int newcrosslinkerGID = (GID1 + 1)*basisnodes_ + GID2;

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

      //getting current position and rotational status of nodes with GID nodeGID[] (based on problem discretization NodeColmap()
      // node 1
      map< int,LINALG::Matrix<3,1> >::const_iterator pos0 = currentpositions.find( nodecolmap.LID(nodeGID[0]) );
      map< int,LINALG::Matrix<3,1> >::const_iterator rot0 = currentrotations.find( nodecolmap.LID(nodeGID[0]) );
      // node 2
      map< int,LINALG::Matrix<3,1> >::const_iterator pos1 = currentpositions.find( nodecolmap.LID(nodeGID[1]) );
      map< int,LINALG::Matrix<3,1> >::const_iterator rot1 = currentrotations.find( nodecolmap.LID(nodeGID[1]) );

      //save positions of nodes between which a crosslinker has to be established in variables xrefe and rotrefe:
      std::vector<double> rotrefe(6);
      std::vector<double> xrefe(6);

      for(int k=0; k<3; k++)
      {
        //the first three positions in xrefe and rotrefe are for the data of node GID1
        if(nodeGID[0] > nodeGID[1])
        {
          //set nodal positions
          xrefe[k ] = (pos0->second)(k);
          xrefe[k+3] = (pos1->second)(k);

          //set nodal rotations (not true ones, only those given in the displacement vector)
          rotrefe[k ] = (rot0->second)(k);
          rotrefe[k+3] = (rot1->second)(k);
        }
        else
        {
          //set nodal positions
          xrefe[k ] = (pos1->second)(k);
          xrefe[k+3] = (pos0->second)(k);

          //set nodal rotations (not true ones, only those given in the displacement vector)
          rotrefe[k ] = (rot1->second)(k);
          rotrefe[k+3] = (rot0->second)(k);
        }
      }

      /*a crosslinker is added on each processor which is row node owner of at least one of its nodes;
       *the processor which is row map owner of the node with the larger GID will be the owner of the new crosslinker element; on the other processors the new
       *crosslinker element will be present as a ghost element only*/
      if(noderowmap.LID(nodeGID[0]) > -1 || noderowmap.LID(nodeGID[1]) > -1)
        AddNewCrosslinkerElement(newcrosslinkerGID,&globalnodeids[0],xrefe,rotrefe,*discret_);
      // add all new elements to contact discretization on all Procs
      if(DRT::INPUT::IntegralValue<int>(statmechparams_,"BEAMCONTACT"))
        AddNewCrosslinkerElement(newcrosslinkerGID,&globalnodeids[0],xrefe,xrefe,beamcmanager->ContactDiscret());
    }
  }
  // synchronization for problem discretization
  discret_->CheckFilledGlobally();
  discret_->FillComplete(true, false, false);

  // synchronization for contact discretization
  if(DRT::INPUT::IntegralValue<int>(statmechparams_,"BEAMCONTACT"))
  {
    beamcmanager->ContactDiscret().CheckFilledGlobally();
    beamcmanager->ContactDiscret().FillComplete(true, false, false);
  }
  //couts
  if(!discret_->Comm().MyPID() && printscreen)
    cout<<"\n\n"<<numsetelements<<" crosslinker element(s) added!"<<endl;
}//void StatMechManager::SearchAndSetCrosslinkers

/*----------------------------------------------------------------------*
 | create a new crosslinker element and add it to your discretization of|
 | choice (public) 				  														 mueller (11/11)|
 *----------------------------------------------------------------------*/
void STATMECH::StatMechManager::AddNewCrosslinkerElement(const int& crossgid,
                                                         int* globalnodeids,
                                                         const std::vector<double>& xrefe,
                                                         const std::vector<double>& rotrefe,
                                                         DRT::Discretization& mydiscret,
                                                         bool addinitlinks)
{
  // get the nodes from the discretization (redundant on both the problem as well as the contact discretization
  DRT::Node* nodes[2] = {	mydiscret.gNode( globalnodeids[0] ) , mydiscret.gNode( globalnodeids[1] )};
  // crosslinker is a beam element
  if(statmechparams_.get<double>("ILINK",0.0) > 0.0)
  {
    // globalnodeids[0] is the larger of the two node GIDs
    Teuchos::RCP<DRT::ELEMENTS::Beam3> newcrosslinker = Teuchos::rcp(new DRT::ELEMENTS::Beam3(crossgid,(mydiscret.gNode(globalnodeids[0]))->Owner() ) );

    newcrosslinker->SetNodeIds(2,globalnodeids);
    newcrosslinker->BuildNodalPointers(&nodes[0]);

    //setting up crosslinker element parameters and reference geometry
    newcrosslinker ->SetCrossSec(statmechparams_.get<double>("ALINK",0.0));
    newcrosslinker ->SetCrossSecShear(1.1*statmechparams_.get<double>("ALINK",0.0));
    newcrosslinker ->SetIyy(statmechparams_.get<double>("ILINK",0.0));
    newcrosslinker ->SetIzz(statmechparams_.get<double>("ILINK",0.0));
    newcrosslinker ->SetIrr(statmechparams_.get<double>("IPLINK",0.0));
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
  else	// crosslinker is a truss element
  {
    Teuchos::RCP<DRT::ELEMENTS::Truss3> newcrosslinker = Teuchos::rcp(new DRT::ELEMENTS::Truss3(crossgid, (mydiscret.gNode(globalnodeids[0]))->Owner() ) );

    newcrosslinker->SetNodeIds(2,globalnodeids);
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
bool STATMECH::StatMechManager::SetCrosslinkerLoom(Epetra_SerialDenseMatrix& LID, const std::map<int, LINALG::Matrix<3, 1> >& currentpositions, Epetra_MultiVector& bspottriadscol)
{
  if(DRT::INPUT::IntegralValue<int>(statmechparams_, "LOOMSETUP"))
  {
    bool setcrosslinker = true;
    // 1) if a crosslink between two vertical filaments is considered, do not set a linker
    if((*filamentnumber_)[(int)LID(0,0)]!= 0 && (*filamentnumber_)[(int)LID(1,0)]!=0)
      setcrosslinker = false;
    else
    {
      // 2) consider direction of the crosslinker
      LINALG::Matrix<3,1> crossdir;
      std::map<int, LINALG::Matrix<3, 1> >::const_iterator pos0 = currentpositions.find((int)LID(0,0));
      std::map<int, LINALG::Matrix<3, 1> >::const_iterator pos1 = currentpositions.find((int)LID(1,0));
      for(int j=0; j<(int)crossdir.M(); j++)
        crossdir(j) = (pos1->second)(j) - (pos0->second)(j);
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
          qnode(l) = bspottriadscol[l][(int)LID(i,0)];
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
 | searches crosslinkers and deletes them if probability check is passed|
 | (public) 																							mueller (7/10)|
 *----------------------------------------------------------------------*/
void STATMECH::StatMechManager::SearchAndDeleteCrosslinkers(const double& timen,
                                                            const double& dt,
                                                            const Epetra_Map& noderowmap,
                                                            const Epetra_Map& nodecolmap,
                                                            const std::map<int, LINALG::Matrix<3, 1> >& currentpositions,
                                                            Epetra_Vector& discol,
                                                            Teuchos::RCP<CONTACT::Beam3cmanager> beamcmanager,
                                                            bool printscreen)
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

  Epetra_Vector punlink(*crosslinkermap_);
  punlink.PutScalar(p);

  // binding spot triads
  Epetra_MultiVector bspottriadscol(nodecolmap,4,true);
  if (DRT::INPUT::IntegralValue<int>(statmechparams_, "HELICALBINDINGSTRUCT"))
    GetBindingSpotTriads(&bspottriadscol);

  // if off-rate is also dependent on the forces acting within the crosslinker
  if(DRT::INPUT::IntegralValue<int>(statmechparams_, "FORCEDEPUNLINKING"))
    ForceDependentOffRate(dt, koff,&punlink,discol);

  // a vector indicating the upcoming deletion of crosslinker elements
  Epetra_Vector delcrosselement(*crosslinkermap_, true);
  delcrosselement.PutScalar(-1.0);

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
          for (int j=0; j<crosslinkerbond_->NumVectors(); j++)
            if ((*crosslinkerbond_)[j][irandom]>-0.9)
              if ((*uniformgen_)() < punlink[irandom])
              {
                // obtain LID and reset crosslinkerbond_ at this position
                int nodeLID = nodecolmap.LID((int) (*crosslinkerbond_)[j][irandom]);
                ((*bspotstatus_)[nodeLID]) = -1.0;
                (*numbond_)[irandom] = 0.0;
                (*crosslinkerbond_)[j][irandom] = -1.0;

                LINALG::SerialDenseMatrix LID(1, 1, true);
                LID(0, 0) = nodeLID;
                CrosslinkerIntermediateUpdate(currentpositions, LID, irandom, bspottriadscol, false);
              }
        }
        break;
        // crosslinker element
        case 2:
        {
          // currently, if an element is deleted, a one-bonded crosslink molecule remains (no simultaneous cut of both bonds)
          std::vector<int> jorder = Permutation(crosslinkerbond_->NumVectors());
          for (int j=0; j<crosslinkerbond_->NumVectors(); j++)
          {
            if ((*uniformgen_)() < punlink[irandom])
            {
              (*numbond_)[irandom] = 1.0;

              // an actual crosslinker element exists
              if((*searchforneighbours_)[irandom]>0.9)
              {
                numdelelements++;
                // random pick of one of the crosslinkerbond entries
                int jrandom = jorder[j];
                // get the nodal LIDs
                int nodeLID = nodecolmap.LID((int)(*crosslinkerbond_)[jrandom][irandom]);

                // in case of a loom setup, we want the vertical filaments to be free of linkers. Hence, we choose the node on the vertical filament to let go of the linker
                // note: jrandom is either 0 or 1. Hence, "!".
                if(DRT::INPUT::IntegralValue<int>(statmechparams_, "LOOMSETUP") && (*filamentnumber_)[nodeLID]==0)
                {
                  jrandom = (!jrandom);
                  nodeLID = nodecolmap.LID((int)(*crosslinkerbond_)[jrandom][irandom]);
                }

                // enter the crosslinker element ID into the deletion list
                delcrosselement[irandom] = (*crosslink2element_)[irandom];

                // vector updates
                (*crosslink2element_)[irandom] = -1.0;
                (*bspotstatus_)[nodeLID] = -1.0;
                (*crosslinkerbond_)[jrandom][irandom] = -1.0;
                // if the linker to be removed occupies to binding spots on the same filament
                if((*crosslinkonsamefilament_)[irandom] > 0.9)
                  (*crosslinkonsamefilament_)[irandom] = 0.0;
              }
              else	// passive crosslink molecule
                (*searchforneighbours_)[irandom] = 1.0;

              for (int k=0; k<crosslinkerbond_->NumVectors(); k++)
                if ((*crosslinkerbond_)[k][irandom]>-0.9)
                {
                  LINALG::SerialDenseMatrix LID(1, 1, true);
                  LID(0,0) = nodecolmap.LID((int)(*crosslinkerbond_)[k][irandom]);
                  CrosslinkerIntermediateUpdate(currentpositions, LID, irandom, bspottriadscol);
                  break;
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
  Epetra_Vector bspotstatusrow(*bspotrowmap_, true);
  Epetra_MultiVector crosslinkerpositionstrans(*transfermap_, true);
  Epetra_MultiVector crosslinkerbondtrans(*transfermap_, 2, true);
  Epetra_Vector numbondtrans(*transfermap_, true);
  Epetra_Vector crosslink2elementtrans(*transfermap_, true);
  Epetra_Vector delcrosselementtrans(*transfermap_, true);

  // export and reimport
  CommunicateVector(bspotstatusrow, *bspotstatus_);
  CommunicateMultiVector(crosslinkerpositionstrans, *crosslinkerpositions_);
  CommunicateMultiVector(crosslinkerbondtrans, *crosslinkerbond_);
  CommunicateVector(numbondtrans, *numbond_);
  CommunicateVector(crosslink2elementtrans, *crosslink2element_);
  CommunicateVector(delcrosselementtrans, delcrosselement);


  // DELETION OF ELEMENTS
  RemoveCrosslinkerElements(*discret_,delcrosselement,&deletedelements_);
  // contact elements
  if(DRT::INPUT::IntegralValue<int>(statmechparams_,"BEAMCONTACT"))
    RemoveCrosslinkerElements(beamcmanager->ContactDiscret(),delcrosselement,&deletedcelements_);

  if(!discret_->Comm().MyPID() && printscreen)
    cout<<numdelelements<<" crosslinker element(s) deleted!"<<endl;
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
void STATMECH::StatMechManager::ForceDependentOffRate(const double& dt, const double& koff0, Epetra_Vector* punlink, Epetra_Vector& discol)
{
  // update element2crosslink_
  element2crosslink_->PutScalar(-1.0);
  for(int i=0; i<crosslinkermap_->NumMyElements(); i++)
    if((*crosslink2element_)[i]>-0.9) // there exists a linker element
    {
      // reverse mapping
      int elelid = discret_->ElementRowMap()->LID((int)(*crosslink2element_)[i]);
      if(elelid>-0.9)
      {
        // 4-noded crosslinker beam element
        if(DRT::INPUT::IntegralValue<int>(statmechparams_, "INTERNODALBSPOTS"))
        {
          //multiple crosslinkers may be mapped to the same element. Hence, the reverse mapping has to be a MultiVector
        }
        else	//conventional beam3 crosslinker element
          (*element2crosslink_)[elelid] = i;
      }
    }

  for(int i=0; i<discret_->ElementRowMap()->NumMyElements(); i++)
    if((*element2crosslink_)[i]>-0.9)
    {
      Teuchos::RCP<Epetra_SerialDenseVector> force = Teuchos::null;
      DRT::Element* crosslinker = discret_->lRowElement(i);
      const DRT::ElementType & eot = crosslinker->ElementType();
      // retrieve internal force vector
      if(eot == DRT::ELEMENTS::Beam3Type::Instance())
        force = dynamic_cast<DRT::ELEMENTS::Beam3*>(crosslinker)->InternalForces();
      else // nothing yet, new crosslinker element soon to come
      {
        // force value for interpolated position xi given by bspotxi_
      }
      // currently, only forces (not moments) considered
      // nodal forces
      LINALG::Matrix<3,1> f0;
      LINALG::Matrix<3,1> f1;
      for(int j=0; j<(int)f0.M(); j++)
      {
        f0(j) = (*force)[j];
        f1(j) = (*force)[6+j];
      }
      // pick the larger absolute force value among the nodes
      int nodeid = -1;
      if(f0.Norm2()>=f1.Norm2())
        nodeid = crosslinker->NodeIds()[0];
      else
        nodeid = crosslinker->NodeIds()[1];
      DRT::Node* node = discret_->lColNode(discret_->NodeColMap()->LID(nodeid));
      std::vector<int> dofnode = discret_->Dof(node);
      // calculate the relative displacement between t_i and t_(i-1)
      double nodereldis = 0.0;
      for(int j=0; j<3; j++)
        nodereldis += (discol[dofnode[j]]-(*disprev_)[dofnode[j]])*(discol[dofnode[j]]-(*disprev_)[dofnode[j]]);
      nodereldis = sqrt(nodereldis);
      // adjusted off-rate according to Bell's equation (Howard, eq 5.10, p.89)
      double koff = koff0 * exp((max(f0.Norm2(),f1.Norm2())*nodereldis)/statmechparams_.get<double>("KT",0.00404531));
      (*punlink)[(int)(*element2crosslink_)[i]] = 1 - exp(-dt*koff);
    }
    else
      (*punlink)[(int)(*element2crosslink_)[i]] = 0.0;

  // Export and and reimport -> redundancy on all Procs
  Epetra_Vector punlinktrans(*transfermap_,true);
  CommunicateVector(punlinktrans, *punlink);
	return;
}// StatMechManager::ForceDependentOffRate

/*----------------------------------------------------------------------*
 | (private) update nodal triads                           mueller 1/11 |
 *----------------------------------------------------------------------*/
void STATMECH::StatMechManager::GetBindingSpotTriads(Epetra_MultiVector* bspottriadscol)
{
  //first get triads at all row nodes
  Epetra_MultiVector bspottriadsrow(*bspotrowmap_, 4, true);
  //update nodaltriads_
  for (int i=0; i<bspotrowmap_->NumMyElements(); i++)
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
        bspottriadsrow[j][i] = ((filele->Qnew())[nodenumber])(j);
    }
    else if (eot == DRT::ELEMENTS::Beam3Type::Instance())
    {
      DRT::ELEMENTS::Beam3* filele = NULL;
      filele = dynamic_cast<DRT::ELEMENTS::Beam3*> (discret_->lRowNode(i)->Elements()[lowestidele]);

      //approximate nodal triad by triad at the central element Gauss point (assuming 2-noded beam elements)
      for(int j=0; j<4; j++)
        bspottriadsrow[j][i] = ((filele->Qnew())[0])(j);
    }
    else
      dserror("Filaments have to be discretized with beam3ii elements for orientation check!!!");


  }
  //export nodaltriadsrow to col map variable
  CommunicateMultiVector(bspottriadsrow, *bspottriadscol, false, true);
	return;
}//StatMechManager::GetNodalTriads

/*----------------------------------------------------------------------*
 | (private) Reduce currently existing crosslinkers by a certain        |
 | percentage.                                              mueller 1/11|
 *----------------------------------------------------------------------*/
void STATMECH::StatMechManager::ReduceNumOfCrosslinkersBy(const int numtoreduce)
{
  if(!discret_->Comm().MyPID())
  {
    cout<<"\n\n==========================================================="<<endl;
    cout<<"-- "<<crosslinkermap_->NumMyElements()<<" crosslink molecules in volume"<<endl;
    cout<<"-- removing "<<statmechparams_.get<int>("REDUCECROSSLINKSBY",0)<<" crosslinkers...\n"<<endl;
  }
  
  int ncrosslink = statmechparams_.get<int>("N_crosslink", 0);

  // check for the correctness of the given input value
  if(numtoreduce>ncrosslink)
    dserror("REDUCECROSSLINKSBY is greater than N_crosslink. Please check your input file!");

  // a vector indicating the upcoming deletion of crosslinker elements
  Epetra_Vector delcrosselement(*crosslinkermap_);
  Epetra_Vector delcrossmolecules(*crosslinkermap_, true);
  delcrosselement.PutScalar(-1.0);

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
            delcrossmolecules[irandom] = 1.0;
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
                int nodeLID = discret_->NodeColMap()->LID((int) (*crosslinkerbond_)[j][irandom]);
                (*bspotstatus_)[nodeLID] = -1.0;
                delcrossmolecules[irandom] = 1.0;
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
              delcrosselement[irandom] = (*crosslink2element_)[irandom];
              delcrossmolecules[irandom] = 1.0;

              ((*bspotstatus_)[(int)(*crosslinkerbond_)[0][irandom]]) = -1.0;
              ((*bspotstatus_)[(int)(*crosslinkerbond_)[1][irandom]]) = -1.0;
            }
            else	// passive crosslink molecule
            {
              numdelmolecules++;
              for (int j=0; j<crosslinkerbond_->NumVectors(); j++)
                if ((*crosslinkerbond_)[j][irandom]>-0.9)
                {
                  // obtain LID and reset crosslinkerbond_ at this position
                  int nodeLID = discret_->NodeColMap()->LID((int)(*crosslinkerbond_)[j][irandom]);
                  ((*bspotstatus_)[nodeLID]) = -1.0;
                  delcrossmolecules[irandom] = 1.0;
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
  Epetra_Vector bspotstatusrow(*bspotrowmap_, true);
  Epetra_Vector delcrosselementtrans(*transfermap_, true);
  Epetra_Vector delcrossmoleculestrans(*transfermap_, true);
  //export and reimport
  CommunicateVector(bspotstatusrow, *bspotstatus_);
  CommunicateVector(delcrosselementtrans, delcrosselement);
  CommunicateVector(delcrossmoleculestrans, delcrossmolecules);

  //PREPARE REBUILDING OF CROSSLINKER MAPS AND CORRESPONDING VECTORS
  // std::vectors for transfer of data to new maps and vectors
  std::vector<std::vector<double> > newcrosslinkerpositions;
  std::vector<std::vector<double> > newcrosslinkerbond;
  std::vector<std::vector<double> > newvisualizepositions;
  std::vector<int> newcrosslinkonsamefilament;
  std::vector<int> newsearchforneighbours;
  std::vector<int> newnumbond;
  std::vector<int> newcrosslink2element;
  // Build temporary vectors
  for(int i=0; i<delcrossmolecules.MyLength(); i++)
  {
    // if crosslinker is kept
    if(delcrossmolecules[i]<0.1)
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
  for (int i=0; i<delcrosselement.MyLength(); i++)
    if (discret_->HaveGlobalElement((int)delcrosselement[i]))
    {
      //save the element by packing before elimination to make it restorable in case that needed
      //cout<<"Proc "<<discret_->Comm().MyPID()<<": deleting element
      //"<<(int)deletecrosslinkerelements[i]<<endl;
      vdata.push_back( DRT::PackBuffer() );
      discret_->gElement((int)delcrosselement[i])->Pack( vdata.back() );
    }
  for ( unsigned i=startsize; i<vdata.size(); ++i )
    vdata[i].StartPacking();
  for (int i=0; i<delcrosselement.MyLength(); i++)
    if (discret_->HaveGlobalElement((int)delcrosselement[i]))
    {
      //save the element by packing before elimination to make it restorable in case that needed
      //cout<<"Proc "<<discret_->Comm().MyPID()<<": deleting element "<<(int)deletecrosslinkerelements[i]<<endl;
      deletedelements_.push_back( std::vector<char>() );
      discret_->gElement((int)delcrosselement[i])->Pack( vdata[deletedelements_.size()-1] );
      discret_->DeleteElement( (int)delcrosselement[i]);
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

  //copy information from the temporary vectors to the adjusted crosslinker vectors
  for(int i=0; i<crosslinkerpositions_->MyLength(); i++)
  {
    for(int j=0; j<crosslinkerpositions_->NumVectors(); j++)
      (*crosslinkerpositions_)[j][i] = newcrosslinkerpositions[i][j];
    for(int j=0; j<visualizepositions_->NumVectors(); j++)
      (*visualizepositions_)[j][i] = (double)newvisualizepositions[i][j];
    for(int j=0; j<crosslinkerbond_->NumVectors(); j++)
      (*crosslinkerbond_)[j][i] = (double)newcrosslinkerbond[i][j];

    (*crosslinkonsamefilament_)[i] = (double)newcrosslinkonsamefilament[i];
    (*searchforneighbours_)[i] = (double)newsearchforneighbours[i];
    (*numbond_)[i] = (double)newnumbond[i];
    (*crosslink2element_)[i] = (double)newcrosslink2element[i];
  }
  if(!discret_->Comm().MyPID())
  {
    cout<<"-- "<<numdelelements<<" crosslinker elements removed"<<endl;
    cout<<"-- "<<numdelmolecules<<" free/one-bonded crosslinker molecules removed"<<endl;
    cout<<"------------------------------------------------------"<<endl;
    cout<<"-- "<<numdelelements+numdelmolecules<<" crosslinkers removed"<<endl;
    cout<<"\n-- "<<crosslinkermap_->NumMyElements()<<" crosslink molecules left in volume"<<endl;
    cout<<"===========================================================\n"<<endl;
  }
  // ReduceNumberOfCrosslinkersBy() is called at the beginning of the time step. Hence, we can just call WriteConv() again to update the converged state
  WriteConv();
  
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
  Epetra_MultiVector randomnumbersrow(*(discret_->ElementRowMap()), randomnumbers->NumVectors());

  for (int i=0; i<randomnumbersrow.MyLength(); i++)
    for (int j=0; j<randomnumbersrow.NumVectors(); j++)
      randomnumbersrow[j][i] = standarddeviation*(*normalgen_)() + meanvalue;

  //export stochastic forces from row map to column map
  CommunicateMultiVector(*randomnumbers,randomnumbersrow,true,false,false);

  return;
} // StatMechManager::SynchronizeRandomForces()

/*----------------------------------------------------------------------*
 | (public) writing restart information for manager objects   cyron 12/08|
 *----------------------------------------------------------------------*/
void STATMECH::StatMechManager::WriteRestart(Teuchos::RCP<IO::DiscretizationWriter> output, double& dt)
{
  output->WriteInt("istart", istart_);
  output->WriteInt("timeintervalstep", timeintervalstep_);
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


  return;
} // StatMechManager::WriteRestart()

/*----------------------------------------------------------------------------*
 | (public) write restart information for fully redundant   Epetra_Multivector|
 | with name "name"                                                cyron 11/10|
 *----------------------------------------------------------------------------*/
void STATMECH::StatMechManager::WriteRestartRedundantMultivector(Teuchos::RCP<IO::DiscretizationWriter> output, const string name, RCP<Epetra_MultiVector> multivector)
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
  timeintervalstep_ = reader.ReadInt("timeintervalstep");
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


  return;
}// StatMechManager::ReadRestart()

/*----------------------------------------------------------------------------*
 | (public) read restart information for fully redundant Epetra_Multivector   |
 | with name "name"                                                cyron 11/10|
 *----------------------------------------------------------------------------*/
void STATMECH::StatMechManager::ReadRestartRedundantMultivector(IO::DiscretizationReader& reader, const string name, Teuchos::RCP<Epetra_MultiVector> multivector)
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
void STATMECH::StatMechManager::WriteConv()
{
  //save relevant class variables at the very end of the time step
  crosslinkerbondconv_ = Teuchos::rcp(new Epetra_MultiVector(*crosslinkerbond_));
  crosslinkerpositionsconv_ = Teuchos::rcp(new Epetra_MultiVector(*crosslinkerpositions_));
  bspotstatusconv_ = Teuchos::rcp(new Epetra_Vector(*bspotstatus_));
  numbondconv_ = Teuchos::rcp(new Epetra_Vector(*numbond_));
  crosslinkonsamefilamentconv_ = Teuchos::rcp(new Epetra_Vector(*crosslinkonsamefilament_));
  crosslink2elementconv_ = Teuchos::rcp(new Epetra_Vector(*crosslink2element_));
  searchforneighboursconv_ = Teuchos::rcp(new Epetra_Vector(*searchforneighbours_));

  //set addedelements_, deletedelements_ empty vectors
  addedelements_.clear();
  deletedelements_.clear();
  if(DRT::INPUT::IntegralValue<int>(statmechparams_,"BEAMCONTACT"))
  {
    addedcelements_.clear();
    deletedcelements_.clear();
  }

  return;
} // StatMechManager::WriteConv()

/*-----------------------------------------------------------------------*
 | communicate Vector to all Processors                    mueller 11/11 |
 *-----------------------------------------------------------------------*/
void STATMECH::StatMechManager::CommunicateVector(Epetra_Vector& InVec,
                                                  Epetra_Vector& OutVec,
                                                  bool doexport,
                                                  bool doimport,
                                                  bool zerofy)
{
  /* zerofy InVec at the beginning of each search except for Proc 0
   * for subsequent export and reimport. This way, we guarantee redundant information
   * on all processors. */

  // first, export the values of OutVec on Proc 0 to InVecs of all participating processors
  Epetra_Export exporter(OutVec.Map(), InVec.Map());
  Epetra_Import importer(OutVec.Map(), InVec.Map());
  if(doexport)
  {
    // zero out all vectors which are not Proc 0. Then, export Proc 0 data to InVec map.
    if(discret_->Comm().MyPID()!=0 && zerofy)
      OutVec.PutScalar(0.0);
    InVec.Export(OutVec, exporter, Add);
  }
  if(doimport)
    OutVec.Import(InVec,importer,Insert);
  return;
}

/*-----------------------------------------------------------------------*
 | communicate MultiVector to all Processors               mueller 11/11 |
 *-----------------------------------------------------------------------*/
void STATMECH::StatMechManager::CommunicateMultiVector(Epetra_MultiVector& InVec,
                                                       Epetra_MultiVector& OutVec,
                                                       bool doexport,
                                                       bool doimport,
                                                       bool zerofy)
{
  // first, export the values of OutVec on Proc 0 to InVecs of all participating processors
  Epetra_Export exporter(OutVec.Map(), InVec.Map());
  Epetra_Import importer(OutVec.Map(), InVec.Map());
  
  if(doexport)
  {
    // zero out all vectors which are not Proc 0. Then, export Proc 0 data to InVec map.
    if(discret_->Comm().MyPID()!=0 && zerofy)
      OutVec.PutScalar(0.0);
    InVec.Export(OutVec, exporter, Add);
  }
  if(doimport)
    OutVec.Import(InVec,importer,Insert);
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
  if(DRT::INPUT::IntegralValue<int>(statmechparams_,"BEAMCONTACT"))
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
 |																							   (public) mueller 3/10|
 *----------------------------------------------------------------------*/
void STATMECH::StatMechManager::GetElementNodeCoords(DRT::Element* element, Teuchos::RCP<Epetra_Vector> dis, LINALG::SerialDenseMatrix& coord, std::vector<int>* lids)
{
  // clear LID vector just in case it was handed over non-empty
  lids->clear();
  for (int j=0; j<element->NumNode(); j++)
  {
    for (int k=0; k<3; k++)
    {
      // obtain k-th spatial component of the reference position of the j-th node
      double referenceposition = ((element->Nodes())[j])->X()[k];
      // get the GIDs of the node's DOFs
      std::vector<int> dofnode = discret_->Dof((element->Nodes())[j]);
      // store the displacement of the k-th spatial component
      double displacement = (*dis)[discret_->DofRowMap()->LID(dofnode[k])];
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
void STATMECH::StatMechManager::UpdateForceSensors(std::vector<int>& sensornodes, int oscdir)
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
      int dofgid = discret_->Dof(0, actnode)[oscdir];
      // now, get the LID
      int collid = discret_->DofColMap()->LID(dofgid);
      // activate force sensor at lid-th position
      (*forcesensor_)[collid] = 1.0;
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
  params.set("CURVENUMBER",statmechparams_.get<int>("CURVENUMBER",-1));
  params.set("OSCILLDIR",statmechparams_.get<int>("OSCILLDIR",-1));
  params.set("STARTTIMEACT",actiontime_->back());
  params.set("DELTA_T_NEW",timestepsizes_->back());
  params.set("PERIODLENGTH",GetPeriodLength());
  if(DRT::INPUT::IntegralValue<int>(statmechparams_,"FORCEDEPUNLINKING") || DRT::INPUT::IntegralValue<int>(statmechparams_,"LOOMSETUP"))
    params.set<string>("internalforces","yes");
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
void STATMECH::StatMechManager::ComputeInternalEnergy(const Teuchos::RCP<Epetra_Vector> dis, double& energy,const double& dt, const std::ostringstream& filename, bool fillzeros, bool writefile)
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
  p.set("OSCILLDIR",statmechparams_.get<int>("OSCILLDIR",-1));
  p.set("PERIODLENGTH", periodlength_);

  discret_->ClearState();
  discret_->SetState("displacement", dis);
  Teuchos::RCP<Epetra_SerialDenseVector> energies = Teuchos::rcp(new Epetra_SerialDenseVector(1));
  energies->Scale(0.0);
  discret_->EvaluateScalars(p, energies);
  discret_->ClearState();
  energy = (*energies)(0);

  if(!discret_->Comm().MyPID() && writefile)
  {
    FILE* fp = NULL;
    fp = fopen(filename.str().c_str(), "a");
    std::stringstream writetofile;
    writetofile<<energy;
    if(fillzeros)
      for(int i=0; i<16; i++)
        writetofile<<"    "<<0;
    writetofile<<endl;
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
| simulation of crosslinker diffusion		        (public) mueller 07/10|
*----------------------------------------------------------------------*/
void STATMECH::StatMechManager::CrosslinkerDiffusion(const Epetra_Vector& dis, double mean, double standarddev, const double &dt)
{
/* Here, the diffusion of crosslink molecules is handled.
 * Depending on the number of occupied binding spots of the molecule, its motion
 * is calculated differently.
 */

  // export row displacement to column map format
  Epetra_Vector discol(*(discret_->DofColMap()), true);
  LINALG::Export(dis, discol);
  /*In this section, crosslinker positions are updated*/
  // diffusion processed by Proc 0
  if (discret_->Comm().MyPID()==0)
  {
    // bonding cases
    for (int i=0; i<crosslinkerpositions_->MyLength(); i++)
    {
      // number of crosslink molecule bonds to filament nodes, larger GID in case of two bonds
      bool set = false;
      int nodegid0 = -1;
      int nodegid1 = -1;

      // getting node gid for crosslink molecules with one bond (both active and passive, i.e. numbond={1,2})
      for (int j=0; j<crosslinkerbond_->NumVectors(); j++)
        if ((*crosslinkerbond_)[j][i]>-0.9 && !set)
        {
          nodegid0 = (int)(*crosslinkerbond_)[j][i];
          set = true;
        }
      // case: crosslinker beam element (numbond=2 AND both gid entries>-1)
      if ((*crosslinkerbond_)[0][i]>-0.9 && (*crosslinkerbond_)[1][i]>-0.9)
        nodegid1 = (int)(*crosslinkerbond_)[1][i];

      // determine crosslinker position according to bonding status (only crosslinker representations, actual crosslinker elements handled thereafter)
      for (int j=0; j<crosslinkerpositions_->NumVectors(); j++)
      {
        switch ((int)(*numbond_)[i])
        {
          // bonding case 1:  no bonds, diffusion
          case 0:
            (*crosslinkerpositions_)[j][i] += standarddev*(*normalgen_)() + mean;
          break;
          // bonding case 2: crosslink molecule attached to one filament
          case 1:
          {
            DRT::Node *node = discret_->lColNode(discret_->NodeColMap()->LID(nodegid0));
            int dofgid = discret_->Dof(node).at(j);
            // set current crosslink position to coordinates of node which it is attached to
            (*crosslinkerpositions_)[j][i] = node->X()[j] + discol[discret_->DofColMap()->LID(dofgid)];
          }
          break;
          // bonding case 3: an actual crosslinker has been established
          case 2:
          {
            // crosslink molecule position is set to position of the node with the larger GID
            if(nodegid1 == -1)
            {
              DRT::Node *node = discret_->lColNode(discret_->NodeColMap()->LID(nodegid0));
              int dofgid = discret_->Dof(node).at(j);
              (*crosslinkerpositions_)[j][i] = node->X()[j]	+ discol[discret_->DofColMap()->LID(dofgid)];
            }
          }
          break;
        }
      }
      // crosslinker position in case of an actual crosslinker element (mid position)
      if(nodegid1 > -1)
      {
        DRT::Node *node0 = discret_->lColNode(discret_->NodeColMap()->LID(nodegid0));
        DRT::Node *node1 = discret_->lColNode(discret_->NodeColMap()->LID(nodegid1));
        for(int j=0; j<crosslinkerpositions_->NumVectors(); j++)
        {
          int dofgid0 = discret_->Dof(node0).at(j);
          int dofgid1 = discret_->Dof(node1).at(j);
          (*crosslinkerpositions_)[j][i] = (node0->X()[j]+discol[discret_->DofColMap()->LID(dofgid0)] + node1->X()[j]+discol[discret_->DofColMap()->LID(dofgid1)])/2.0;
        }
      }
    }// end of i-loop
    // check for compliance with periodic boundary conditions if existent
    if (periodlength_->at(0) > 0.0)
      CrosslinkerPeriodicBoundaryShift(*crosslinkerpositions_);
  } // if(discret_->Comm().MyPID()==0)
  // Update by Broadcast: copy this information to all processors
  Epetra_MultiVector crosslinkerpositionstrans(*transfermap_, 3, true);
  CommunicateMultiVector(crosslinkerpositionstrans, *crosslinkerpositions_);
  
  return;
}// StatMechManager::CrosslinkerDiffusion

/*----------------------------------------------------------------------*
 | update crosslink molecule positions	                 								|
 |																				        (public) mueller 08/10|
 *----------------------------------------------------------------------*/
void STATMECH::StatMechManager::CrosslinkerIntermediateUpdate(const std::map<int, LINALG::Matrix<3, 1> >& currentpositions,
                                                    const LINALG::SerialDenseMatrix& LID, const int& crosslinkernumber,
                                                    Epetra_MultiVector& bspottriadscol,
                                                    bool coupledmovement)
{
  // case: one-bonded crosslink molecule (i.e. two cases: +1 bond (starting from 0 bonds) or -1 bond (molecule is free to diffuse again)
  if (LID.M()==1 && LID.N()==1)
  {
    // set molecule position to node position
    map<int, LINALG::Matrix<3, 1> >::const_iterator pos0 = currentpositions.find((int)LID(0,0));
    if (coupledmovement)
    {
      // if the binding spots are oriented according to the double helical structure of f-actin
      if(DRT::INPUT::IntegralValue<int>(statmechparams_,"HELICALBINDINGSTRUCT"))
      {
        double ronebond = statmechparams_.get<double> ("R_LINK", 0.0) / 2.0;
        // get binding spot quaternion
        LINALG::Matrix<4,1> qbspot;
        LINALG::Matrix<3,3> R;
        for(int i=0; i<(int)qbspot.M(); i++)
          qbspot(i) = bspottriadscol[i][(int)LID(0,0)];
        LARGEROTATIONS::quaterniontotriad(qbspot,R);

        LINALG::Matrix<3,1> tangent;
        LINALG::Matrix<3,1> normal;
        // retrieve tangential and normal vector from binding spot quaternions
        for (int i=0; i<(int)R.M(); i++)
        {
          tangent(i) = R(i,0);
          normal(i) = R(i,1);
        }
        // rotation matrix around tangential vector by given angle
        RotationAroundFixedAxis(tangent,normal,(*bspotorientations_)[(int)LID(0,0)]);

        // calculation of the visualized point lying in the direction of the rotated normal
        for (int i=0; i<crosslinkerpositions_->NumVectors(); i++)
          (*crosslinkerpositions_)[i][crosslinkernumber] = (pos0->second)(i) + ronebond*normal(i);
      }
      else // simplistic binding spot model (no orientation)
      {
        for (int i=0; i < crosslinkerpositions_->NumVectors(); i++)
          (*crosslinkerpositions_)[i][crosslinkernumber] = (pos0->second)(i);
      }
    }
    else
    {
      if(DRT::INPUT::IntegralValue<int>(statmechparams_,"HELICALBINDINGSTRUCT"))
      {
        // nothing to be done
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
  }
  // case: crosslinker element
  if (LID.M()==2 && LID.N()==1)
  {
    if(DRT::INPUT::IntegralValue<int>(statmechparams_,"HELICALBINDINGSTRUCT"))
    {
      map<int,LINALG::Matrix<3,1> >::const_iterator pos = currentpositions.find((int)LID(0,0));
      for (int i=0; i<crosslinkerpositions_->NumVectors(); i++)
        (*crosslinkerpositions_)[i][crosslinkernumber] = (pos->second)(i);
      pos = currentpositions.find((int)LID(1,0));
      for (int i=0; i<crosslinkerpositions_->NumVectors(); i++)
      {
        (*crosslinkerpositions_)[i][crosslinkernumber] += (pos->second)(i);
        (*crosslinkerpositions_)[i][crosslinkernumber] /= 2.0;
      }
    }
    else
    {
      int chosenlid = max((int)LID(0,0), (int)LID(1,0));
      // in case of a loom network, we want the crosslinker position to lie on the horizontal filament
      if(DRT::INPUT::IntegralValue<int>(statmechparams_, "LOOMSETUP"))
        for(int i=0; i<LID.M(); i++)
          if((*filamentnumber_)[(int)LID(i,0)] == 0)
          {
            chosenlid = (int)LID(i,0);
            break;
          }
      map<int,LINALG::Matrix<3,1> >::const_iterator updatedpos = currentpositions.find(chosenlid);
      for (int i=0; i<crosslinkerpositions_->NumVectors(); i++)
        (*crosslinkerpositions_)[i][crosslinkernumber] = (updatedpos->second)(i);
    }
  }
  return;
}// StatMechManager::CrosslinkerIntermediateUpdate

/*----------------------------------------------------------------------*
 | Initialize crosslinker positions  			        (public) mueller 07/10|
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
  if(DRT::INPUT::IntegralValue<int>(statmechparams_, "HELICALBINDINGSTRUCT"))
  {
    // rise and rotation per binding spot (default values: 2.77nm and -166.15deg (left-handed, one-start helix))
    double riseperbspot = statmechparams_.get<double>("RISEPERBSPOT",0.00277);
    double rotperbspot = statmechparams_.get<double>("ROTPERBSPOT", -2.8999);

    //getting a vector consisting of pointers to all filament number conditions set
    std::vector<DRT::Condition*> filaments(0);
    discret_->GetCondition("FilamentNumber",filaments);

    // apply new beam element with intermediate binding spot positions
    if(DRT::INPUT::IntegralValue<int>(statmechparams_, "INTERNODALBSPOTS"))
    {
      // temporary vectors
      // binding spot orientations relative to second triad vector
      std::vector<double> bspotorientations;
      // value of binding spot in parameter space
      std::vector<double> bspotxi;
      // maps binding spot to element
      std::vector<int> 		bspot2element;
      // vector holding all bspot gids
      std::vector<int> 		bspotgids;
      // vector signaling the existence of a binding spot on this proc
      std::vector<int> 		bspotonproc;

      //loop over filaments to create redundant column map
      for (int i=0; i<(int)filaments.size(); i++)
      {
        //get a pointer to nodal cloud covered by the current condition
        const std::vector<int>* nodeids = filaments[i]->Nodes();
        // retrieve filament length using first and last node of the current filament
        DRT::Node* node0 = discret_->lColNode(discret_->NodeColMap()->LID((*nodeids)[0]));
        DRT::Node* node1 = discret_->lColNode(discret_->NodeColMap()->LID((*nodeids)[((int)nodeids->size())-1]));

        // length of the current filament
        double lfil = 0.0;
        for(int j=0; j<3; j++)
          lfil += (node1->X()[j]-node0->X()[j])*(node1->X()[j]-node0->X()[j]);
        lfil = sqrt(lfil);

        // element length (assuming equal node spacing)
        double elelength = lfil/(double)(filaments[i]->Nodes()->size()-1);
        // add as many binding spots to the i-th filament as possible
        int bspot = 0;
        while((double)bspot*riseperbspot <= lfil)
        {
          // put binding spot id into the vector in order to create a fully overlapping binding spot map afterwards
          bspotgids.push_back((int)bspotgids.size());

          // calculate the orientation of the binding spot as well as the respective curve parameter
          /* determine element current element GID:
           * As set up by the Inputfile-Generator, the filament element GID can be retrieved by knowing the filament number and the lesser node number.
           * Attention: If that property is changed, one should of course refer to node->Elements(). In the meanwhile, this method seems faster*/
          int elegid = node0->Id() - i + (int)(floor((double)bspot*riseperbspot/lfil));
          // add element GID if this element is a row map element on this processor
          if(discret_->ElementRowMap()->LID(elegid)>-1)
          {
            // mark binding spot as being on this proc
            bspotonproc.push_back(1);
            bspot2element.push_back(elegid);

            // orientation [0;2pi] and parameter xi_bs of the binding spot
            double phibs = 0.0;
            double xibs = 0.0;
            if(bspot==0)
            {
              bspotorientations.push_back(phibs);
              bspotxi.push_back(xibs);
            }
            else
            {
              phibs = (((double)(bspot)*rotperbspot)/(2.0*M_PI)-floor(((double)(bspot-1)*rotperbspot)/(2.0*M_PI)))*2.0*M_PI;
              xibs = (double)bspot*riseperbspot/elelength - floor((double)bspot*riseperbspot/elelength) - 0.5;
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
          bspotrowgids.push_back(i);	// note: since column map is fully overlapping: i=col. LID = GID
      bspotrowmap_ = Teuchos::rcp(new Epetra_Map((int)bspotgids.size(), (int)bspotrowgids.size(), &bspotrowgids[0], 0, discret_->Comm()));

      // vectors
      // initialize class vectors
      bspotstatus_ = Teuchos::rcp(new Epetra_Vector(*bspotcolmap_));
      bspotstatus_->PutScalar(-1.0);
      bspotorientations_ = Teuchos::rcp(new Epetra_Vector(*bspotcolmap_,true));
      bspotxi_ = Teuchos::rcp(new Epetra_Vector(*bspotcolmap_,true));
      bspot2element_ = Teuchos::rcp(new Epetra_Vector(*bspotrowmap_, true));

      Epetra_Vector bspotorientationsrow(*bspotrowmap_);
      Epetra_Vector bspotxirow(*bspotrowmap_);
      for(int i=0; i<bspotrowmap_->NumMyElements(); i++)
      {
        bspotorientationsrow[i] = bspotorientations[i];
        bspotxirow[i] = bspotxi[i];
        (*bspot2element_)[i] = bspot2element[i];
      }
      // make information on orientations and curve parameters redundant on all procs
      Epetra_Import bspotimporter(*bspotcolmap_, *bspotrowmap_);
      bspotorientations_->Import(bspotorientationsrow,bspotimporter,Insert);
      bspotxi_->Import(bspotxirow,bspotimporter,Insert);
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
      Epetra_MultiVector bspottriadscol(*bspotcolmap_,4,true);
      GetBindingSpotTriads(&bspottriadscol);

      DRT::Node* node0 = NULL;
      DRT::Node* node1 = NULL;
      LINALG::Matrix<3,1> vectonfil;
      // note default value of pi/2 leads to zero offset from the second direction of the material frame
      double bspotoffset = statmechparams_.get<double>("BSPOTOFFSET",-M_PI/2.0);
      // difference angle to fit position of first nodes of two filaments
      double deltatheta = 0.0;
      int bspotgid = 0;

      for(int i=0; i<(int)filaments.size(); i++)
      {
        // vector between first nodes of consecutive filament PAIRS (0-1; 2-3; ...)
        double angularoffset = 0.0;
        if(i%2==0 && (int)filaments.size()>1)
        {
          // vector from first node of filament i to the first one of i+1
          node0 = discret_->lColNode(discret_->NodeColMap()->LID((*filaments[i]->Nodes())[0]));
          if(i<(int)filaments.size()-1)
            node1 = discret_->lColNode(discret_->NodeColMap()->LID((*filaments[i+1]->Nodes())[0]));
          else
            node1 = discret_->lColNode(discret_->NodeColMap()->LID((*filaments[i-1]->Nodes())[0]));

          for(int j=0; j<(int)vectonfil.M(); j++)
            vectonfil(j) = node1->X()[j]-node0->X()[j];
          vectonfil.Scale(1.0/vectonfil.Norm2());

          // get the second direction of the material triad (since filament is initially straight, all triads of a filament are the same)
          LINALG::Matrix<3,3> bspottriad;
          // auxiliary variable for storing a triad in quaternion form
          LINALG::Matrix<4, 1> qnode;
          LINALG::Matrix<3,1> firstdir;
          LINALG::Matrix<3,1> secdir;
          int nodelid = discret_->NodeColMap()->LID((*filaments[i]->Nodes())[0]);
          for (int l=0; l<4; l++)
            qnode(l) = bspottriadscol[l][nodelid];
          LARGEROTATIONS::quaterniontotriad(qnode, bspottriad);
          for (int l=0; l<(int)bspottriad.M(); l++)
          {
            firstdir(l) = bspottriad(l,0);
            secdir(l) = bspottriad(l,1);
          }
          // relative rotation so that a pair of filaments has facing first binding spots
          deltatheta = acos(secdir.Dot(vectonfil));
          // test rotation to see in which direction to rotate secdir
          RotationAroundFixedAxis(firstdir, secdir, deltatheta);
          // "0.1" just some value >0 and <pi as we want to segregate parallel from antiparallel
          if(acos(secdir.Dot(vectonfil))>0.1)
            deltatheta = -deltatheta;

          angularoffset = deltatheta + fabs(bspotoffset);
        }
        else
          angularoffset  = -(M_PI-deltatheta-bspotoffset);

        for(int j=0; j<(int)filaments[i]->Nodes()->size(); j++)
        {
          // determine orientation (relative to second triad vector). No use of riseperbspot, since somewhat handled by element length
          (*bspotorientations_)[bspotcolmap_->LID(bspotgid)] = j*rotperbspot + angularoffset;
          // assumption:
          if(bspotrowmap_->LID(bspotgid)>-1)
            (*bspot2element_)[bspotrowmap_->LID(bspotgid)] = discret_->lRowNode(discret_->NodeRowMap()->LID(bspotgid))->Elements()[0]->Id();
          // increment binding spot GID
          bspotgid++;
        }
      }
    }
  }
  else // conventional binding spot geometry (bspots are identical to nodes)
  {
    bspotcolmap_ = Teuchos::rcp(new Epetra_Map(*(discret_->NodeColMap())));
    bspotrowmap_ = Teuchos::rcp(new Epetra_Map(*(discret_->NodeRowMap())));
    bspotstatus_ = Teuchos::rcp(new Epetra_Vector(*bspotcolmap_));
    bspotstatus_->PutScalar(-1.0);
  }

  // create density-density-correlation-function map with
  if(DRT::INPUT::IntegralValue<INPAR::STATMECH::StatOutput>(statmechparams_, "SPECIAL_OUTPUT")==INPAR::STATMECH::statout_densitydensitycorr ||
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
  }
  /*cout<<"start indices: ";
  for(int i=0; i<(int)startindex_->size(); i++)
    cout<<(*startindex_)[i]<<" ";
  cout<<endl;*/

  // in case the internal forces of the crosslinker affect the off-rate
  if(DRT::INPUT::IntegralValue<int>(statmechparams_, "FORCEDEPUNLINKING"))
  {
    element2crosslink_ = Teuchos::rcp(new Epetra_Vector(*(discret_->ElementRowMap())));
    element2crosslink_->PutScalar(-1.0);
    disprev_ = Teuchos::rcp(new Epetra_Vector(*(discret_->DofColMap()),true));
  }

  std::vector<double> upperbound = *periodlength_;
  // handling both cases: with and without periodic boundary conditions
  if (periodlength_->at(0) == 0.0)
    for(int i=0; i<(int)upperbound.size(); i++)
      upperbound.at(i) = statmechparams_.get<double> ("MaxRandValue", 0.0);

  crosslinkerpositions_ = Teuchos::rcp(new Epetra_MultiVector(*crosslinkermap_, 3, true));

  for (int i=0; i<crosslinkerpositions_->MyLength(); i++)
    for (int j=0; j<crosslinkerpositions_->NumVectors(); j++)
      (*crosslinkerpositions_)[j][i] = upperbound.at(j) * (*uniformgen_)();

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

  return;
}//StatMechManager::CrosslinkerMoleculeInit

/*----------------------------------------------------------------------*
| Set crosslinkers whereever possible before the first time step        |
|                                                (private) mueller 11/11|
*----------------------------------------------------------------------*/
void STATMECH::StatMechManager::SetInitialCrosslinkers(Teuchos::RCP<CONTACT::Beam3cmanager> beamcmanager)
{
  if(DRT::INPUT::IntegralValue<int>(statmechparams_, "DYN_CROSSLINKERS"))
  {
    Epetra_MultiVector bspottriadscol(*bspotcolmap_,4,true);
    if(DRT::INPUT::IntegralValue<int>(statmechparams_,"HELICALBINDINGSTRUCT"))
      GetBindingSpotTriads(&bspottriadscol);

    const Epetra_Map noderowmap = *discret_->NodeRowMap();
    const Epetra_Map nodecolmap = *discret_->NodeColMap();
    //node positions and rotations (binding spot positions and rotations when applying 4-noded beam element)
    std::map<int, LINALG::Matrix<3, 1> > currentpositions;
    std::map<int, LINALG::Matrix<3, 1> > currentrotations;

    // Vectors hold zero->ok, since no displacement at this stage of the simulation
    Epetra_Vector discol(*discret_->DofColMap(), true);
    GetNodePositions(discol, currentpositions, currentrotations, true);
    // do initial octree build
    if(beamcmanager!=Teuchos::null)
      beamcmanager->OcTree()->OctTreeSearch(currentpositions);

    // generate a random order of binding spots and crosslinkers (only on Proc 0)
    std::vector<int> randbspot;
    std::vector<int> randlink;

    //1. set singly bound linkers
    int numsinglybound = 0;
    if(discret_->Comm().MyPID()==0)
    {
      randbspot = Permutation(bspotcolmap_->NumMyElements());
      randlink = Permutation(statmechparams_.get<int>("N_crosslink", 0));
      int numbspots = statmechparams_.get<int>("INITOCCUPIEDBSPOTS",0);

      if(numbspots>bspotcolmap_->NumMyElements())
        dserror("Given number of initially occupied binding spots (%i) exceeds the total binding spot count (%i)! Check your input file!",numbspots, bspotcolmap_->NumMyElements());

      // first, establish specified number of singly bound crosslinkers
      for(int i=0; i<numbspots; i++)
      {
        int firstbspot = randbspot[i];
        // if this binding spot is still unoccupied
        if((*bspotstatus_)[firstbspot]<-0.9)
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
          CrosslinkerIntermediateUpdate(currentpositions, LID, currlink, bspottriadscol);
          // update visualization
          for (int j=0; j<visualizepositions_->NumVectors(); j++)
            (*visualizepositions_)[j][currlink] = (*crosslinkerpositions_)[j][currlink];
        }
      }
    }

    // Communication, first stage (not sure if we need to communicate all of these vectors here (or later)
    // transfer vectors
    Epetra_Vector bspotstatusrow(*discret_->NodeRowMap(),true);
    Epetra_Vector numbondtrans(*transfermap_, true);
    Epetra_MultiVector crosslinkerpositionstrans(*transfermap_,3,true);
    Epetra_MultiVector visualizepositionstrans(*transfermap_,3,true);
    Epetra_MultiVector crosslinkerbondtrans(*transfermap_,2,true);

    CommunicateVector(bspotstatusrow, *bspotstatus_);
    CommunicateVector(numbondtrans, *numbond_);
    CommunicateMultiVector(crosslinkerpositionstrans, *crosslinkerpositions_);
    CommunicateMultiVector(visualizepositionstrans, *visualizepositions_);
    CommunicateMultiVector(crosslinkerbondtrans, *crosslinkerbond_);

    // this section creates crosslinker finite elements
    // Commented because of performance issues  when adding a large number of elements at once.
    // 2. Now, parallely search for neighbour nodes
    Teuchos::RCP<Epetra_MultiVector> neighbourslid;
    if(statmechparams_.get<int>("SEARCHRES",1)>0)
      PartitioningAndSearch(currentpositions,bspottriadscol, neighbourslid);

    // 3. create double bonds
    // a vector indicating the crosslink molecule which is going to constitute a crosslinker element
    Epetra_Vector addcrosselement(*crosslinkermap_, true);
    int numsetelements = 0;

    if(discret_->Comm().MyPID()==0)
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

          // if second binding exists and spot is unoccupied; mandatory: occupy two binding spots on DIFFERENT filaments
          if((*neighbourslid)[currneighbour][currlink]>-0.1 && (*filamentnumber_)[secondbspot]!=(*filamentnumber_)[randbspot[i]])
          {
            if((*bspotstatus_)[secondbspot] < -0.9)
            {
              numsetelements++;
              addcrosselement[currlink] = 1.0;
              // establish double bond to the first given neighbour
              // attach it to the second binding spot
              (*bspotstatus_)[secondbspot] = currlink;
              (*numbond_)[currlink] = 2.0;
              Epetra_SerialDenseMatrix LID(2,1);
              for(int k=0; k<crosslinkerbond_->NumVectors(); k++)
                if((*crosslinkerbond_)[k][currlink]<0.1)
                {
                  LID(k,0) = secondbspot;
                  (*crosslinkerbond_)[k][currlink] = bspotcolmap_->GID(secondbspot);
                }
                else
                  LID(k,0) = (*crosslinkerbond_)[k][currlink];

              bool intersection = false;
              if(beamcmanager!=Teuchos::null)
              {
                Epetra_SerialDenseMatrix nodecoords(3,2);
                for(int k=0; k<nodecoords.M(); k++)
                {
                  nodecoords(k,0) = ((currentpositions.find((int)LID(0,0)))->second)(k);
                  nodecoords(k,1) = ((currentpositions.find((int)LID(1,0)))->second)(k);
                }
                for(int k=0; k<(int)LID.M(); k++)
                {
                  intersection = beamcmanager->OcTree()->IntersectBBoxesWith(nodecoords, LID);
                  if(intersection)
                    break;
                }
              }

              if(!intersection)
              {
                //update crosslinker position
                CrosslinkerIntermediateUpdate(currentpositions, LID, currlink, bspottriadscol);
                // update visualization
                for (int k=0; k<visualizepositions_->NumVectors(); k++)
                  (*visualizepositions_)[k][currlink] = (*crosslinkerpositions_)[k][currlink];
                break;
              }
            }
          }
        }
      }
    }
    else
    {
      numbond_->PutScalar(0.0);
      crosslinkerbond_->PutScalar(0.0);
      crosslinkerpositions_->PutScalar(0.0);
      visualizepositions_->PutScalar(0.0);
      addcrosselement.PutScalar(0.0);
    }

    // Communication, second stage
    Epetra_Vector addcrosselementtrans(*transfermap_,true);
    bspotstatusrow.PutScalar(0.0);
    crosslinkerpositionstrans.PutScalar(0.0);
    visualizepositionstrans.PutScalar(0.0);
    crosslinkerbondtrans.PutScalar(0.0);
    numbondtrans.PutScalar(0.0);

    // export and reimport
    CommunicateVector(addcrosselementtrans, addcrosselement);
    CommunicateVector(bspotstatusrow, *bspotstatus_);
    CommunicateVector(numbondtrans, *numbond_);
    CommunicateMultiVector(crosslinkerpositionstrans, *crosslinkerpositions_);
    CommunicateMultiVector(visualizepositionstrans, *visualizepositions_);
    CommunicateMultiVector(crosslinkerbondtrans, *crosslinkerbond_);

    // ADDING ELEMENTS
    // add elements to problem discretization (processor specific)
    for(int i=0; i<addcrosselement.MyLength(); i++)
    {
      if(addcrosselement[i]>0.9)
      {
        // obtain node GID
        int nodeGID[2] = {	(int)(*crosslinkerbond_)[0][i],(int)(*crosslinkerbond_)[1][i]};
        // determine smaller and larger of the GIDs
        int GID2 = min(nodeGID[0],nodeGID[1]);
        int GID1 = max(nodeGID[0],nodeGID[1]);
        int globalnodeids[2] = {GID1, GID2};
        // calculate element GID
        int newcrosslinkerGID = (GID1 + 1)*basisnodes_ + GID2;

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

        (*crosslink2element_)[i] = newcrosslinkerGID;

        //getting current position and rotational status of nodes with GID nodeGID[] (based on problem discretization NodeColmap()
        // node 1
        map< int,LINALG::Matrix<3,1> >::const_iterator pos0 = currentpositions.find( nodecolmap.LID(nodeGID[0]) );
        map< int,LINALG::Matrix<3,1> >::const_iterator rot0 = currentrotations.find( nodecolmap.LID(nodeGID[0]) );
        // node 2
        map< int,LINALG::Matrix<3,1> >::const_iterator pos1 = currentpositions.find( nodecolmap.LID(nodeGID[1]) );
        map< int,LINALG::Matrix<3,1> >::const_iterator rot1 = currentrotations.find( nodecolmap.LID(nodeGID[1]) );

        //save positions of nodes between which a crosslinker has to be established in variables xrefe and rotrefe:
        std::vector<double> rotrefe(6);
        std::vector<double> xrefe(6);

        for(int k=0; k<3; k++)
        {
          //the first three positions in xrefe and rotrefe are for the data of node GID1
          if(nodeGID[0] > nodeGID[1])
          {
            //set nodal positions
            xrefe[k ] = (pos0->second)(k);
            xrefe[k+3] = (pos1->second)(k);

            //set nodal rotations (not true ones, only those given in the displacement vector)
            rotrefe[k ] = (rot0->second)(k);
            rotrefe[k+3] = (rot1->second)(k);
          }
          else
          {
            //set nodal positions
            xrefe[k ] = (pos1->second)(k);
            xrefe[k+3] = (pos0->second)(k);

            //set nodal rotations (not true ones, only those given in the displacement vector)
            rotrefe[k ] = (rot1->second)(k);
            rotrefe[k+3] = (rot0->second)(k);
          }
        }
        // add linker element if one of the element's nodes is on the node row map of the processor
        // we do not need to worry about contact yet: if enabled, the contact discretization is initialized
        // after the statmechmanager...
        if(noderowmap.LID(nodeGID[0]) > -1 || noderowmap.LID(nodeGID[1]) > -1)
          AddNewCrosslinkerElement(newcrosslinkerGID,&globalnodeids[0],xrefe,rotrefe,*discret_, true);
        if(DRT::INPUT::IntegralValue<int>(statmechparams_,"BEAMCONTACT"))
          AddNewCrosslinkerElement(newcrosslinkerGID,&globalnodeids[0],xrefe,xrefe,beamcmanager->ContactDiscret(), true);
      }
    }

    // synchronization for problem discretization
    discret_->CheckFilledGlobally();
    discret_->FillComplete(true, false, false);

    if(DRT::INPUT::IntegralValue<int>(statmechparams_,"BEAMCONTACT"))
    {
      beamcmanager->ContactDiscret().CheckFilledGlobally();
      beamcmanager->ContactDiscret().FillComplete(true, false, false);
    }

    // reduce number of contact pairs to zero again to avoid unnecessary computations
    if(beamcmanager!=Teuchos::null)
      beamcmanager->ResetPairs();

    //Gmsh output
    if(DRT::INPUT::IntegralValue<int>(statmechparams_,"GMSHOUTPUT"))
    {
      std::ostringstream filename;
      filename << "./GmshOutput/InitLinks.pos";
      Epetra_Vector disrow(*discret_->DofRowMap(), true);
      GmshOutput(disrow,filename,0);
    }
    if(beamcmanager!=Teuchos::null && DRT::INPUT::IntegralValue<INPAR::STATMECH::StatOutput>(statmechparams_, "SPECIAL_OUTPUT")==INPAR::STATMECH::statout_octree)
    {
      // "-2" for initial octree output
      beamcmanager->OcTree()->OctTreeSearch(currentpositions,-2);
      beamcmanager->ResetPairs();
    }
    //couts
    if(!discret_->Comm().MyPID())
    {
      cout<<"====setting initial crosslinkers===="<<endl;
      cout<<"singly bound: "<<numsinglybound-numsetelements<<endl;
      cout<<"doubly bound: "<<numsetelements<<endl;
      cout<<"===================================="<<endl;
    }
  }
  return;
}

/*----------------------------------------------------------------------*
 |  Periodic Boundary Shift for crosslinker diffusion simulation        |
 |                                                (public) mueller 07/10|
 *----------------------------------------------------------------------*/
void STATMECH::StatMechManager::CrosslinkerPeriodicBoundaryShift(Epetra_MultiVector& crosslinkerpositions)
{
  for (int i=0; i<crosslinkerpositions.MyLength(); i++)
    for (int j=0; j<crosslinkerpositions.NumVectors(); j++)
    {
      if (crosslinkerpositions[j][i] > periodlength_->at(j))
        crosslinkerpositions[j][i] -= periodlength_->at(j);
      if (crosslinkerpositions[j][i] < 0.0)
        crosslinkerpositions[j][i] += periodlength_->at(j);
    }
  return;
}// StatMechManager::CrosslinkerPeriodicBoundaryShift

/*----------------------------------------------------------------------*
 |  Evaluate DBCs in case of periodic BCs (public)         mueller 03/12|
 *----------------------------------------------------------------------*/
void STATMECH::StatMechManager::EvaluateDirichletStatMech(ParameterList&                     params,
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
 |  Evaluate DBCs in case of periodic BCs (public)         mueller 02/10|
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
  	cout<<"DBC Evaluation time: "<<t_end-t_start<<endl;
#endif // #ifdef MEASURETIME

  return;
} //EvaluateDirichletPeriodic()

/*----------------------------------------------------------------------*
 | Determine if application of DBCs starts at a given time  mueller 5/12|
 *----------------------------------------------------------------------*/
bool STATMECH::StatMechManager::DBCStart(Teuchos::ParameterList& params)
{
  // get the current time
  double time = params.get<double>("total time", 0.0);
  double starttime = actiontime_->at((int)(actiontime_->size()-1));
  double dt = params.get<double>("delta time", 0.01);
  if (time<0.0) dserror("t = %f ! Something is utterly wrong here. The total time should be positive!", time);

  if(time < starttime && fabs(starttime-time)>dt/1e4)
    return false;
  else
    return true;
}

/*----------------------------------------------------------------------*
 | Gather information on where DBCs are to be applied in the case of    |
 | viscoelastic measurements                   (private)  mueller 05/12 |
 *----------------------------------------------------------------------*/
void STATMECH::StatMechManager::DBCOscillatoryMotion(Teuchos::ParameterList& params,
                                                     Teuchos::RCP<Epetra_Vector> dis,
                                                     Teuchos::RCP<Epetra_Vector> deltadbc)
{
  // get the absolute time, time step size, time curve number and oscillation direction
  const double time = params.get<double>("total time", 0.0); // target time (i.e. timen_)
  double dt = params.get<double>("delta time", 0.01);
  int curvenumber = statmechparams_.get<int>("CURVENUMBER",0)-1;
  int oscdir = statmechparams_.get<int>("OSCILLDIR",-1);
  double amp = statmechparams_.get<double>("SHEARAMPLITUDE",0.0);

//----------------------------------- sanity checks
  if(periodlength_->at(0)!=periodlength_->at(1) || periodlength_->at(0)!=periodlength_->at(2) || periodlength_->at(1)!=periodlength_->at(2))
    dserror("Check PERIODLENGTH in your input file! This method only works for a cubic boundary volume");
  if(periodlength_->at(0)<=0.0)
    dserror("PERIODLENGTH = %f! This method work only in conjunction with periodic boundary conditions!", periodlength_->at(0));
  if(oscdir>2 || oscdir<0 || curvenumber<0)
    dserror("In case of imposed oscillatory motion, please define the StatMech parameters OSCILLDIR={0,1,2} and/or CURVENUMBER correctly");

//------------------------------------------------------gather dbc node set(s)
  // loop over row elements
  // after the the very first entry into the following part, reenter only if the Dirichlet Node Set is dynamically updated
  if(!useinitdbcset_)
  {
    // clear node sets, as they will be filled for the first time or updated with new node sets
    dbcnodesets_.clear();
    // vectors to manipulate DBC
    std::vector<int> oscillnodes;
    std::vector<int> fixednodes;
    oscillnodes.clear();
    fixednodes.clear();
    // toggle vector
    std::vector<bool> nodedbcstatus(discret_->NumMyRowNodes(), false);

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
        // get nodal coordinates and LIDs of the nodal DOFs
        // store LIDs of element node dofs
        std::vector<int> doflids;
        GetElementNodeCoords(element, dis, coord, &doflids);
//-----------------------detect broken/fixed/free elements and fill position vector
        // determine existence and location of broken element
        if(CheckForBrokenElement(coord,cut))
        {
          // loop over number of cuts (columns)
          for(int n=0; n<cut.N(); n++)
          {
            int nodelidn = discret_->NodeRowMap()->LID(element->Nodes()[n]->Id());
            int nodelidnp = discret_->NodeRowMap()->LID(element->Nodes()[n+1]->Id());

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
                fixednodes.push_back(element->Nodes()[n]->Id());
                // add GID of oscillating node to osc.-nodes-vector
                oscillnodes.push_back(element->Nodes()[n+1]->Id());

                // incremental Dirichlet displacement for an oscillating node (all DOFs except oscdir = 0.0)
                // time curve increment
                double tcincrement = 0.0;
                if(curvenumber>-1)
                  tcincrement = DRT::Problem::Instance()->Curve(curvenumber).f(time) -
                                DRT::Problem::Instance()->Curve(curvenumber).f(time-dt);
                (*deltadbc)[doflids.at(numdof*(n+1)+oscdir)] = amp*tcincrement;
//                cout<<"t="<<time<<", dt="<<dt<<", (*deltadbc)["<<doflids.at(numdof*(n+1)+oscdir)<<"] = "<<amp*tcincrement<<endl;
//                cout<<"      amp="<<amp<<", f(t) = "<<DRT::Problem::Instance()->Curve(curvenumber).f(time)<<", f(t-dt) = "<<DRT::Problem::Instance()->Curve(curvenumber).f(time-dt)<<endl;
              }// end of case 1
              // case 2: broken element (in z-dir); node_n oscillates, node_n+1 is fixed in dir. of oscillation
              if(cut(2,n)==2.0)
              {
                nodedbcstatus.at(nodelidn) = true;
                nodedbcstatus.at(nodelidnp) = true;
                oscillnodes.push_back(element->Nodes()[n]->Id());
                fixednodes.push_back(element->Nodes()[n+1]->Id());
                // oscillating node
                double tcincrement = 0.0;
                if(curvenumber>-1)
                  tcincrement = DRT::Problem::Instance()->Curve(curvenumber).f(time) -
                                DRT::Problem::Instance()->Curve(curvenumber).f(time-dt);
                (*deltadbc)[doflids.at(numdof*n+oscdir)] = amp*tcincrement;
              } // end of case 2

//              if(cut(2,n)>0.0)
//              {
//                cout<<"Element "<<element->Id()<<endl;
//                cout<<"Nodes   "<<element->Nodes()[0]->Id()<<", "<<element->Nodes()[1]->Id()<<endl;
//                cout<<"doflids ";
//                for(int i=0; i<(int)doflids.size(); i++)
//                  cout<<doflids[i]<<" ";
//                cout<<endl;
//              }
            }
          }
        }
      }
    }
//-----------------------------------------update node sets
    dbcnodesets_.push_back(oscillnodes);
    dbcnodesets_.push_back(fixednodes);

//---------check/set force sensors anew for each time step
    // add DOF LID where a force sensor is to be set
    UpdateForceSensors(dbcnodesets_[0], oscdir);
    //if(!discret_->Comm().MyPID())
    //  cout<<"\n=========================================="<<endl;
    //for(int pid=0; pid<discret_->Comm().NumProc();pid++)
    //{
    //  if(pid==discret_->Comm().MyPID())
    //    cout<<"UpdateForceSensors_"<<pid<<": "<<oscillnodes.size()<< " nodes @ t="<<time<<endl;
    //  discret_->Comm().Barrier();
    //}
    //cout<<"==========================================\n"<<endl;
  }
  else // only ever use the fixed set of nodes after the first time dirichlet values were imposed
  {
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
//      cout<<"t="<<time<<", dt="<<dt<<", (*deltadbc)["<<dofnode[oscdir]<<"] = "<<amp*tcincrement<<endl;
//      cout<<"      amp="<<amp<<", f("<<time<<") = "<<DRT::Problem::Instance()->Curve(curvenumber).f(time)<<", f("<<time-dt<<") = "<<DRT::Problem::Instance()->Curve(curvenumber).f(time-dt)<<endl;
    }
    // dofs of fixednodes_ remain untouched since fixednodes_ dofs are 0.0
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
    case INPAR::STATMECH::dbctype_sheartrans:
    {
      int oscdir = statmechparams_.get<int>("OSCILLDIR",-1);
      onoff.at(oscdir) = 1;
    }
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
  const vector<int>* nodeids = cond.Nodes();
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
  cout<<"Node Ids: ";
  for(int i=0; i<(int)nodeids->size(); i++)
    cout<<nodeids->at(i)<<" ";
  cout<<"onoff: ";
  for(int i=0; i<(int)discret_->Dof(0,discret_->gNode(nodeids->at(0))).size(); i++)
    cout<<onoff->at(i)<<" ";
  cout<<endl;*/

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

