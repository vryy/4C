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
#include "statmech_search.H"

#include <Teuchos_Time.hpp>

#include "../drt_inpar/inpar_statmech.H"
#include "../drt_inpar/inpar_structure.H"
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
#include "../drt_spring3/spring3.H"
#include "../drt_beam3r/beam3r.H"
#include "../drt_beam3eb/beam3eb.H"
#include "../drt_beam3cl/beam3cl.H"
#include "../drt_truss3/truss3.H"
#include "../drt_truss3cl/truss3cl.H"
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
numoutputstep_(1),
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
  SeedRandomGenerators(0,-1);

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

  // initial clearing of dbc management vectors
  nbcnodesets_.clear();

  // filament number conditions: quick reference which node belongs to which filament
  InitializeFilamentNumbers();

  // Create randomly assigned endnode sets with barbed and pointed ends
//  CreateBarbedPointedNodeSets();

  // force sensor conditions: for viscoelastic and creep simulations, forces need to be meassured...
  InitializeForceSensors();

  // read points in time from input file where certain actions take place
  SetStartStep(parameters);

  // Initialize crosslinker positions (point-like particles while not doubly-bound)
  CrosslinkerMoleculeInit();

  // Initialize crosslinker positions (point-like particles while not doubly-bound)
//  MonomerInitiatePos();


  return;
}// StatMechManager::StatMechManager


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
  //
  if (linkermodel_ != statmech_linker_none && statmechparams_.get<int>("N_crosslink",0)>0)
  {
#ifdef MEASURETIME
    const double t1 = Teuchos::Time::wallTime();
#endif
    //column maps cointaining binding spot positions and their spatial orientation
    Teuchos::RCP<Epetra_MultiVector> bspotpositions = Teuchos::rcp(new Epetra_MultiVector(*bspotcolmap_,3,true));
    // allocate memory to binding spot rotations only if linkers are not Truss3, Truss3CL elements
    // In case of Truss C/L, bspotrotations indicate the tangent of nodal displacement vector
    Teuchos::RCP<Epetra_MultiVector> bspotrotations = Teuchos::rcp(new Epetra_MultiVector(*bspotcolmap_,3,true));;

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
//    MonomerDiffusion(0.0, standarddev, dt);

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

      SearchAndDeleteCrosslinkers(istep, timen, dt, bspotpositions, bspotrotations, discol,beamcmanager,printscreen);

#ifdef MEASURETIME
      t5 = Teuchos::Time::wallTime();
#endif

      ElementToCrosslinkMapping(element2crosslink_);
    }

    // reset thermal energy to new value (simple value change for now, maybe Load Curve later on)
    if(fabs(timen-actiontime_->at(1))<(dt/1e3))
      statmechparams_.set("KT",statmechparams_.get<double>("KTACT",statmechparams_.get<double>("KT",0.0)));

    // // Set a certain number of double-bonded crosslinkers free
    if(fabs(timen-dt-actiontime_->back())<(dt/1e3) && statmechparams_.get<int>("REDUCECROSSLINKSBY",0)>0)
      ReduceNumOfCrosslinkersBy(statmechparams_.get<int>("REDUCECROSSLINKSBY",0), beamcmanager);

    // change length of active crosslinkers
    if(linkermodel_ == statmech_linker_active || linkermodel_ == statmech_linker_activeintpol || linkermodel_ == statmech_linker_myosinthick)
      ChangeActiveLinkerLength(timen, dt, crosslinkeractlengthconv_, crosslinkeractlength_, printscreen,false,bspotpositions,bspotrotations);

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
                                                      double  timeconverged,
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

/*------------------------------------------------------------------------------------------------*
 | Retrieve binding spot positions                                        (public)   mueller 10/12|
 *------------------------------------------------------------------------------------------------*/
void STATMECH::StatMechManager::GetBindingSpotPositions(Epetra_Vector& discol,
                                                        Teuchos::RCP<Epetra_MultiVector> bspotpositions,
                                                        Teuchos::RCP<Epetra_MultiVector> bspotrotations)
{
  /*in preparation for later decision whether a crosslink should be established between two nodes (binding spots) we first store the
   * current positions of all column map nodes (column map binding spots) in the map currentpositions; additionally we store the rotational displacements
   * analogously in a map currentrotations for later use in setting up reference geometry of crosslinkers; the maps
   * currentpositions and currentrotations relate positions and rotations to a local column map node Id, respectively*/

  // in case the four-noded crosslinker beam or Truss element is used, currentpositions and currentrotations have to be set up another way
  if(linkermodel_ == statmech_linker_stdintpol ||
     linkermodel_ == statmech_linker_activeintpol ||
     linkermodel_ == statmech_linker_bellseqintpol ||
     linkermodel_ == statmech_linker_myosinthick)
    GetInterpolatedBindingSpotPositions(discol, bspotpositions, bspotrotations);
  else  // conventional crosslinker beam3r, Truss 3 element, i.e. binding spots coincide with nodes
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
    DRT::Node* node = discret_->lRowNode(i);
    //get GIDs of this node's degrees of freedom
    std::vector<int> dofnode = discret_->Dof(node);
    //ask the truss element about the first element the first node is connected to
    DRT::Element* Element = node->Elements()[0];

    if(getpositions)
      for(int j=0; j<bspotpositionsrow->NumVectors(); j++)
        (*bspotpositionsrow)[j][i] = node->X()[j] + discol[discret_->DofColMap()->LID(dofnode[j])];

    //if node has also rotational degrees of freedom
    if (discret_->NumDof(node) == 6 && getrotations)
    {
      //Check via dynamic cast, if it's a beam3eb element
      DRT::ELEMENTS::Beam3eb* BeamElement = dynamic_cast<DRT::ELEMENTS::Beam3eb*>(Element);
      //if not, tell the user to use beam3 instead
      if(BeamElement==NULL)
      {
        for(int j=0; j<bspotrotationsrow->NumVectors(); j++)
          (*bspotrotationsrow)[j][i] = discol[discret_->DofColMap()->LID(dofnode[j+3])];
      }
      else
      {
        LINALG::Matrix<3,1> trefNodeAux(true);
        trefNodeAux=BeamElement->Tref()[0];
        for(int j=0; j<bspotrotationsrow->NumVectors(); j++)
          (*bspotrotationsrow)[j][i]=trefNodeAux(j)+discol[discret_->DofColMap()->LID(dofnode[j+3])];
      }
    }
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

/*--------------------------------------------------------------------------------------------------*
 | Updates positions and rotations of binding spots for 4 noded  elements              mueller 10/12|
 *--------------------------------------------------------------------------------------------------*/
void STATMECH::StatMechManager::GetInterpolatedBindingSpotPositions(const Epetra_Vector&              discol,
                                                                    Teuchos::RCP<Epetra_MultiVector>  bspotpositions,
                                                                    Teuchos::RCP<Epetra_MultiVector>  bspotrotations)

{
  bool getpositions = true;
  bool getrotations = true;
  if(bspotpositions==Teuchos::null)
    getpositions = false;
  if(bspotrotations==Teuchos::null)
    getrotations = false;

  //sanity checks
  if(!getpositions && !getrotations)
    dserror("Both vectors are null!");

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
    LINALG::Matrix<1,4>  Ibp_hermite;         //Hermite shape functions for beam3EB element
    std::vector<double> Tcurr_(6); // Vector containing tangent of Nodal Positions.
    // length of first element and second element required for hermite shape function interpolation
    double LengthofElementatRef;
    LINALG::Matrix<3,1> ThetaXi;
    // Vector containing Nodal Positions. Only translational DOFs are evaluated
    std::vector<double> position(6);
    std::vector<double> InitialPosition(6);
    //Quaternions at relevant filament nodes
    LINALG::Matrix<4,1>  bQ;
    std::vector<LINALG::Matrix<4,1> >  nQ(2);
    DRT::ELEMENTS::Beam3r* filele = NULL;

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
        if (filelement->ElementType()==DRT::ELEMENTS::Beam3rType::Instance())
          filele = dynamic_cast<DRT::ELEMENTS::Beam3r*> (discret_->gElement(elegid));

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



        /* Compute the tangential degrees of freedom at current time step.
         * Since filaments are assumed to be straight at reference configuration,
         * initial values are zero.
         */
        if(filelement->ElementType()==DRT::ELEMENTS::Beam3ebType::Instance())
        {
          InitialPosition[0] = node0->X()[0];
          InitialPosition[1] = node0->X()[1];
          InitialPosition[2] = node0->X()[2];
          InitialPosition[3] = node1->X()[0];
          InitialPosition[4] = node1->X()[1];
          InitialPosition[5] = node1->X()[2];

          // NodeShift Initial positions
          UnshiftPositions(InitialPosition);

          //Tangents at reference configuration
          std::vector<LINALG::Matrix<3,1> > Tref_(2);
          double norm2=0;
          for(int node = 0; node<2 ; node++)
          {
            for(int dof = 0; dof< 3 ; dof++ )
            {
              Tref_[node](dof) =  InitialPosition[dof+3]- InitialPosition[dof];
            }
            norm2 = Tref_[node].Norm2();
            LengthofElementatRef=norm2;
            Tref_[node].Scale(1/norm2);
          }

          Tcurr_[0] = Tref_[0](0) + discol[discret_->DofColMap()->LID(dofnode0[3])];
          Tcurr_[1] = Tref_[0](1) + discol[discret_->DofColMap()->LID(dofnode0[4])];
          Tcurr_[2] = Tref_[0](2) + discol[discret_->DofColMap()->LID(dofnode0[5])];
          Tcurr_[3] = Tref_[1](0) + discol[discret_->DofColMap()->LID(dofnode1[3])];
          Tcurr_[4] = Tref_[1](1) + discol[discret_->DofColMap()->LID(dofnode1[4])];
          Tcurr_[5] = Tref_[1](2) + discol[discret_->DofColMap()->LID(dofnode1[5])];
        }

        // NodeShift
        UnshiftPositions(position);

        // To interpolate rotations we need to get the two nodal Triads (Quaternions) of the Element, if no truss element
        if(statmechparams_.get<double>("ILINK",0.0)>0.0)
        {
          for(int j=0;j<2;j++)
            for(int k=0;k<4;k++)
              nQ[j](k) = (filele->Qnew()[j])(k);
        }
        prevelegid = elegid;
      }

      // Interpolation of Positions
      if (filelement->ElementType()==DRT::ELEMENTS::Beam3ebType::Instance())
      {
        //Here we need the position of the internodal Binding Spot at time t=0, when the filament beam elements where initialized
        DRT::UTILS::shape_function_hermite_1D(Ibp_hermite,(double)(*bspotxi_)[bspotrowmap_->GID(i)],LengthofElementatRef,filelement->Shape());
        currpos.PutScalar(0);
        for(int j=0;j<3;j++)
          currpos(j)= Ibp_hermite(0)*position[j]+Ibp_hermite(1)*Tcurr_[j]+Ibp_hermite(2)*position[j+3]+Ibp_hermite(3)*Tcurr_[j+3];
      }
      else
      {
        DRT::UTILS::shape_function_1D(Ibp,(double)(*bspotxi_)[bspotrowmap_->GID(i)],filelement->Shape());
        currpos.PutScalar(0);
        for(int j=0;j<3;j++)
          currpos(j)= Ibp(0)*position[j]+Ibp(1)*position[j+3];
      }

      /*if bspot currently has coordinate value greater than statmechparams_.get<double>("PeriodLength",0.0),or smaller than 0
       *it is shifted by -statmechparams_.get<double>("PeriodLength",0.0) sufficiently often to lie again in the domain*/
      for(int j=0;j<(int)(periodlength_->size());j++)
      {
        if(currpos(j) > periodlength_->at(j))
          currpos(j) -= periodlength_->at(j)*floor(currpos(j)/periodlength_->at(j));

        if(currpos(j) < 0.0)
          currpos(j) -= periodlength_->at(j)*floor(currpos(j)/periodlength_->at(j));
      }

      //Interpolation of Rotations
      if(statmechparams_.get<double>("ILINK",0.0)>0.0)
      {
        currrot.PutScalar(0.0);
        InterpolateTriadonBindingSpot(Ibp,nQ,currrot,bQ);
      }


      for(int j=0;j<3;j++)
      {
        (*bspotpositionsrow)[j][i] = currpos(j);
        if (getrotations)
        (*bspotrotationsrow)[j][i] = currrot(j);
      }
    }
    // Export row map entries to column maps
    if(getpositions)
    CommunicateMultiVector(bspotpositionsrow, bspotpositions, false, true);
    if(getrotations)
      CommunicateMultiVector(bspotrotationsrow, bspotrotations, false, true);






  // DEBUGGING CHECK wether bspotpositions make sense

//  if(discret_->Comm().MyPID()==0)
//  {
//  std::cout<< " **************** CHECK WETHER BSPOTPOSITIONS MAKE SENSE DISTANCES  ****************** " << std::endl;
//  double ll;
//  std::cout << " Anzahl BSPOTS : " << bspotpositions->MyLength() << std::endl;
//  for(int ii=1;ii<bspotpositions->MyLength();ii++)
//  {
//      ll= std::sqrt(std::pow((*bspotpositions)[2][ii]-(*bspotpositions)[2][ii-1],2)+pow((*bspotpositions)[1][ii]-(*bspotpositions)[1][ii-1],2)+pow((*bspotpositions)[0][ii]-(*bspotpositions)[0][ii-1],2));
//      std::cout<< " BSPOTDISTANCE "<< ll<< std::endl;
//  }
//  for(int ii=0;ii<bspotpositions->MyLength();ii++)
//    std::cout<< "BSPOTXI " << (*bspotxi_)[ii] << " BSPOT2NODES " << (*bspot2nodes_)[0][ii] << " "<<(*bspot2nodes_)[1][ii] <<   std::endl;

//  std::cout<< " **************** CHECK WETHER BSPOTPOSITIONS MAKE SENSE POSITIONS  ****************** " << std::endl;
//  std::cout<< *bspotpositions << std::endl;
//  }

}//GetInterpolatedBindingSpotPositions


/*------------------------------------------------------------------------------------*
 | (private) update internodal triads at binding positions                            |
 | This is done by a transfer of the current rotations into Quaternions Mueller 11/12 |
 *------------------------------------------------------------------------------------*/
void STATMECH::StatMechManager::GetBindingSpotTriads(const Teuchos::RCP<Epetra_MultiVector> bspotrotations,
                                                     Teuchos::RCP<Epetra_MultiVector>       bspottriads)
{
  if(bspottriads!=Teuchos::null)
  {
    if (DRT::INPUT::IntegralValue<int>(statmechparams_, "CHECKORIENT") ||
        filamentmodel_ == statmech_filament_helical ||
        linkermodel_ == statmech_linker_activeintpol ||
        linkermodel_ == statmech_linker_bellseqintpol ||
        linkermodel_ == statmech_linker_myosinthick)
    {
      if (linkermodel_ == statmech_linker_stdintpol || linkermodel_ == statmech_linker_activeintpol || linkermodel_ == statmech_linker_bellseqintpol || linkermodel_ ==  statmech_linker_myosinthick)
        GetInterpolatedBindingSpotTriads(bspotrotations, bspottriads);
      else
        GetElementBindingSpotTriads(bspottriads);
    }
  }
  return;
}


/*------------------------------------------------------------------------------------------------*
 | (private) update internodal triads at binding positions for GMSH visualizations                |
 | This is done by a transfer of the current rotations into Quaternions           Mukherjee 12/14 |
 *------------------------------------------------------------------------------------------------*/
void STATMECH::StatMechManager::GetBindingSpotTriadsGMSH(const Teuchos::RCP<Epetra_MultiVector> bspotrotations,
                                                     Teuchos::RCP<Epetra_MultiVector>       bspottriads)
{
  if(bspottriads!=Teuchos::null && bspotrotations!=Teuchos::null)
  {
    if (DRT::INPUT::IntegralValue<int>(statmechparams_, "CHECKORIENT") ||
        filamentmodel_ == statmech_filament_helical ||
        linkermodel_ == statmech_linker_activeintpol ||
        linkermodel_ == statmech_linker_bellseqintpol ||
        linkermodel_ == statmech_linker_myosinthick)
    {
      if (linkermodel_ == statmech_linker_stdintpol || linkermodel_ == statmech_linker_activeintpol || linkermodel_ == statmech_linker_bellseqintpol || linkermodel_ ==  statmech_linker_myosinthick)
        GetInterpolatedBindingSpotTriads(bspotrotations, bspottriads);
      else
        GetElementBindingSpotTriads(bspottriads);
    }
  }
  return;
}

/*----------------------------------------------------------------------*
 | (private) update nodal triads                           mueller 1/11 |
 *----------------------------------------------------------------------*/
void STATMECH::StatMechManager::GetElementBindingSpotTriads(Teuchos::RCP<Epetra_MultiVector> nodetriads)
{
  //first get triads at all row nodes
  // In case of Beam3eb elements, nodal triad returns only the nodal tangent at Reference config.
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
    //if element is of type beam3r get nodal triad
    if (eot == DRT::ELEMENTS::Beam3rType::Instance())
    {
      DRT::ELEMENTS::Beam3r* filele = NULL;
      filele = dynamic_cast<DRT::ELEMENTS::Beam3r*> (discret_->lRowNode(i)->Elements()[lowestidele]);

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
    else if (eot == DRT::ELEMENTS::Beam3ebType::Instance())
    {
      DRT::ELEMENTS::Beam3eb* filele = NULL;
      filele = dynamic_cast<DRT::ELEMENTS::Beam3eb*> (discret_->lRowNode(i)->Elements()[lowestidele]);

      //check whether crosslinker is connected to first or second node of that element
      int nodenumber = 0;
      if(discret_->lRowNode(i)->Id() == ((filele->Nodes())[1])->Id() )
        nodenumber = 1;

      //save nodal triad of this node in nodaltriadrow
      for(int j=0; j<3; j++)
        (*nodetriadsrow)[j][i] = ((filele->Tref())[nodenumber])(j);
      // Unlike quaternion, Triad has only 3 components. Therefore the 4th component is set to zero.
      (*nodetriadsrow)[3][i]=0.0;
    }
    else
      dserror("Filaments have to be discretized with beam3r or Beam3eb elements for orientation check!!!");
  }
  // communicate the appropriate vector
  if(nodetriads->MyLength()==discret_->NodeColMap()->NumMyElements())
    CommunicateMultiVector(nodetriadsrow, nodetriads, false, true, false);
  else
    nodetriads = Teuchos::rcp(new Epetra_MultiVector(*nodetriadsrow));
}//StatMechManager::GetNodalTriads

/*------------------------------------------------------------------------------------*
 | (private) update internodal triads at binding positions                            |
 | This is done by a transfer of the current rotations into Quaternions mueller 10/12 |
 *------------------------------------------------------------------------------------*/
void STATMECH::StatMechManager::GetInterpolatedBindingSpotTriads(const Teuchos::RCP<Epetra_MultiVector> bspotrotations,
                                                                 Teuchos::RCP<Epetra_MultiVector>       bspottriads)
{
  if(bspottriads!=Teuchos::null)
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
  return;
}


/*------------------------------------------------------------------------------*
 | Interpolation of the binding spot triads between nodes         mueller 10/12 |
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
  rPsil.PutScalar(0.0);
  for(int node=0;node<2;node++)
    for(int i=0;i<3;i++)
      rPsil(i)+= Ibp(node)*rPsili[node](i);

  LARGEROTATIONS::angletoquaternion(rPsil,rQl);
  LARGEROTATIONS::quaternionproduct(rQl,rQr,rQxi);
  LARGEROTATIONS::quaterniontoangle(rQxi,ThetaXi);
  bQxi=rQxi;

  return;
}


/*----------------------------------------------------------------------*
 | Assign crosslink molecules and nodes to volume partitions            |
 |                                                (public) mueller 08/10|
 *----------------------------------------------------------------------*/
void STATMECH::StatMechManager::PartitioningAndSearch(const Teuchos::RCP<Epetra_MultiVector> bspotpositions,
                                                      const Teuchos::RCP<Epetra_MultiVector> bspottriadscol,
                                                      Teuchos::RCP<Epetra_MultiVector>&      neighbourslid)
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
    for(int i=0; i<(int)bspotpositions->NumVectors(); i++)
    {
      for(int j=0; j<(int)bspotpositions->MyLength(); j++)
      {
        // max
        if(limit[2*i+1]<(*bspotpositions)[i][j])
          limit[2*i+1] = (*bspotpositions)[i][j];
        // min
        if(limit[2*i]>(*bspotpositions)[i][j])
          limit[2*i] = (*bspotpositions)[i][j];
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
    for(int j=0;j<bspotpositions->MyLength();j++)
    {
      int partition = (int)std::floor(((*bspotpositions)[i][j]-limit[2*i])/(limit[2*i+1]-limit[2*i])*(double)(*searchres_)[i]);
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
  DetectBindingSpots(bspotpositions, bspotinpartition, numbondtrans, crosslinkerpositionstrans, crosslinkpartitiontrans, bspottriadscol, neighbourslid);

  return;
}//void StatMechManager::PartitioningAndSearch

/*----------------------------------------------------------------------*
 | detect neighbour nodes to crosslink molecules (public) mueller (7/10)|
 *----------------------------------------------------------------------*/
void STATMECH::StatMechManager::DetectBindingSpots(const Teuchos::RCP<Epetra_MultiVector>              bspotpositions,
                                                   const std::vector<std::vector<std::vector<int> > >& bspotinpartition,
                                                   const Teuchos::RCP<Epetra_Vector>                   numbondtrans,
                                                   const Teuchos::RCP<Epetra_MultiVector>              crosslinkerpositionstrans,
                                                   const Teuchos::RCP<Epetra_MultiVector>              crosslinkpartitionstrans,
                                                   const Teuchos::RCP<Epetra_MultiVector>              bspottriadscol,
                                                   Teuchos::RCP<Epetra_MultiVector>&                   neighbourslid)
  {
  // distribute information of crosslinkerbond_ to processorspecific maps
  Teuchos::RCP<Epetra_MultiVector> crosslinkerbondtrans = Teuchos::rcp(new Epetra_MultiVector(*transfermap_, 2, true));
  CommunicateMultiVector(crosslinkerbondtrans, crosslinkerbond_,true,false,false,true);

  std::vector<std::vector<int> > neighbournodes(crosslinkpartitionstrans->MyLength(), std::vector<int>());

  int maxneighbourslocal = 0;
  int maxneighboursglobal = 0;

  for(int part=0; part<crosslinkpartitionstrans->MyLength(); part++)
  {
    // i.e. numbond!=2.0
    if((*crosslinkpartitionstrans)[0][part]>-0.9)
    {
      double rmin = (statmechparams_.get<double>("R_LINK", 0.0)-statmechparams_.get<double>("DeltaR_LINK", 0.0)) / 2.0;
      double rmax = (statmechparams_.get<double>("R_LINK", 0.0)+statmechparams_.get<double>("DeltaR_LINK", 0.0)) / 2.0;
      if ((int)(*numbondtrans)[part]>0.9)
      {
        rmin *= 2.0;
        rmax *= 2.0;
      }

      // first component
      for(int ilayer=(int)(*crosslinkpartitionstrans)[0][part]-1; ilayer<(int)(*crosslinkpartitionstrans)[0][part]+2; ilayer++)
      {
        if(ilayer>-1 && ilayer<(*searchres_)[0])
        {
          for(int i=0; i<(int)bspotinpartition[0][ilayer].size(); i++)
          {
            int tmplid = (int)bspotinpartition[0][ilayer][i];
            // second component
            for(int jlayer=(int)(*crosslinkpartitionstrans)[1][part]-1; jlayer<(int)(*crosslinkpartitionstrans)[1][part]+2; jlayer++)
            {
              if(jlayer>-1 && jlayer<(*searchres_)[1])
              {
                for(int j=0; j<(int)bspotinpartition[1][jlayer].size(); j++)
                {
                  if(bspotinpartition[1][jlayer][j]==tmplid)
                  {
                    //third component
                    for(int klayer=(int)(*crosslinkpartitionstrans)[2][part]-1; klayer<(int)(*crosslinkpartitionstrans)[2][part]+2; klayer++)
                    {
                      if(klayer>-1 && klayer<(*searchres_)[2])
                      {
                        for(int k=0; k<(int)bspotinpartition[2][klayer].size(); k++)
                        {
                          if(bspotinpartition[2][klayer][k]==tmplid)
                          {
                            // calculate distance crosslinker-bspot
                            LINALG::Matrix<3, 1> distance;
                            for (int l=0; l<(int)distance.M(); l++)
                              distance(l) = (*crosslinkerpositionstrans)[l][part]-(*bspotpositions)[l][bspotcolmap_->GID(tmplid)];
                            // 1. criterion: distance between linker and binding spot within given interval
                            if(distance.Norm2()<rmax && distance.Norm2()>rmin)
                            {
                              // further calculations in case of helical binding spot geometry and singly bound crosslinkers
                              if(filamentmodel_ == statmech_filament_helical && statmechparams_.get<double>("ILINK",0.0)>0.0)
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
                                  qnode(l) = (*bspottriadscol)[l][bspotcolmap_->GID(tmplid)];
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
                                  crossbspotdiff(l) = (*crosslinkerpositionstrans)[l][part]-(*bspotpositions)[l][bspotcolmap_->GID(tmplid)];
                                // line parameter of the intersection point of the line through the binding spot with the orientation of the binding spot.
                                // a)lambda must be larger than zero in order to lie within the cone
                                double lambda = (crossbspotdiff.Dot(bspotvec))/(bspotvec.Dot(bspotvec));
                                if(lambda>0)
                                {
                                  distance.Scale(1.0/distance.Norm2()); // normalized vector
                                  // b) the angle between the binding spot orientation and the vector between bspot and crosslinker must be smaller than PHIBSPOT
                                  //    Only then does the crosslinker lie within the cone
                                  if(acos(fabs(bspotvec.Dot(distance))) < statmechparams_.get<double>("PHIBSPOT", 0.524))
                                    bspot1=true;
                                  if((int)(*numbondtrans)[part]<0.1)
                                    bspot2=true;
                                  else //check wether first Bspot fullfills orientation criterion
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
                                       qnode(l) = (*bspottriadscol)[l][bspotID];
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
                                      distance.Scale(1.0/distance.Norm2()); // normalized vector
                                    // b) the angle between the binding spot orientation and the vector between bspot and crosslinker must be smaller than PHIBSPOT
                                    //    Only then does the crosslinker lie within the cone
                                    if(acos(fabs(bspotvec.Dot(distance))) < statmechparams_.get<double>("PHIBSPOT", 0.524))
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

  /* assign "-2" -> 'empty'
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
 | Binding spot / linker search by Octree                mueller (09/13)|
 *----------------------------------------------------------------------*/
void STATMECH::StatMechManager::DetectBindingSpotsOctree(const Teuchos::RCP<Epetra_MultiVector> bspotpositions,
                                                         Teuchos::RCP<Epetra_MultiVector>       bspottriadscol,
                                                         Teuchos::RCP<Epetra_MultiVector>&      neighbourslid)
{
  std::vector<std::vector<int> > neighbournodes(crosslinkerpositions_->MyLength(), std::vector<int>());
  int maxneighbourslocal = 0;
  int maxneighboursglobal = 0;

  //Create and build BspotOctree
  //length of bounding boxes arund binding spots (should be larger or equal the crosslinker length + deviation)
  double bindingradius =(statmechparams_.get<double>("R_LINK", 0.0)+statmechparams_.get<double>("DeltaR_LINK", 0.0));
  //Constructor
  Teuchos::RCP<STATMECH::SEARCH::Octree> bspottree = Teuchos::rcp(new STATMECH::SEARCH::Octree(periodlength_,discret_ ,bspotrowmap_,bspotcolmap_, statmechparams_.get<int> ("MAXBSPOTOCTREEDEPTH", 5), statmechparams_.get<int> ("MINBSPOTSINOCT", 5), bindingradius));
  //Build Octree based on the location of the binding spot positions
  bspottree->BuildOctree(bspotpositions);
  //Locate crosslinker positions (i.e. build a vector that maps each crosslinker to one octant)
  bspottree->LocatePositions(crosslinkerpositions_,crosslinkermap_);

 //loop over all octants (this is done in parallel)
  for(int octLID=0;octLID<bspottree->BBoxesInOctRow()->MyLength(); octLID++)
  {
    //retrieve GID of octant from octreemap
    int octGID = (int)bspottree->BBoxesInOctRow()->Map().GID(octLID);
    //loop over all crosslinker in this octant
    for(int i=0; i< bspottree->CrosslinkerInOctants()->NumVectors();i++)
    {
      //break i loop if entry is -9. This means no further crosslinker molecules in this octant
      if((*bspottree->CrosslinkerInOctants())[i][octGID]<-1)
        break;
      // only consider crosslinkers that are unbound (i.e. numbond<2.0)
      if((*numbond_)[bspotcolmap_->LID((int)(*bspottree->CrosslinkerInOctants())[i][octGID])]>1.9)
        continue;
      //temporarily store Crosslinker ID
      int clID=(int)(*bspottree->CrosslinkerInOctants())[i][octGID];
      /*adjust length of linker in case of singly bound crosslinker (not in case of helical bstruct).
       * This is because freely diffusing crosslinkers (numbond=-1) are modelled as point in space located at the center of a crosslinker molecule.
       * If a crosslinker has bound to a binding spot its position is set to the position of this binding spot, i.e. on one of the ends of the crosslinker molecule.
       * Hence the search-radius has to be set to the full crosslinker length*/
      double rmin = (statmechparams_.get<double>("R_LINK", 0.0)-statmechparams_.get<double>("DeltaR_LINK", 0.0)) / 2.0;
      double rmax = (statmechparams_.get<double>("R_LINK", 0.0)+statmechparams_.get<double>("DeltaR_LINK", 0.0)) / 2.0;
      if ((*numbond_)[clID]>0.9)
      {
        rmin *= 2.0;
        rmax *= 2.0;
      }

      //loop over all binding spots within this octant (i.e. possible neighbours)
      for(int k=0; k<bspottree->BBoxesInOctRow()->NumVectors();k++)
      {
        //skip if entry is -9 meaning there are nor further bspots
        if((*bspottree->BBoxesInOctRow())[k][octLID]<-1)
          continue;
        //store binding spot id temporarily
        int tmplid = (int)(*bspottree->BBoxesInOctRow())[k][octLID];
        //Now perform Orientation checks to identify possible binding partners

        // calculate distance crosslinker-node
        LINALG::Matrix<3, 1> difference;
        for (int l=0; l<(int)difference.M(); l++)
          difference(l) = (*crosslinkerpositions_)[l][clID]-(*bspotpositions)[l][bspotcolmap_->GID(tmplid)];
        //difference(l) = (*crosslinkerpositions_)[l][clID]-(*bspotpositions)[l][bspotcolmap_->LID(bspotrowmap_->GID(tmplid))];
        // 1. criterion: distance between linker and binding spot within given interval
        if(difference.Norm2()<rmax && difference.Norm2()>rmin)
        {
          // further calculations in case of helical binding spot geometry and singly bound crosslinkers
          if(filamentmodel_ == statmech_filament_helical)
          {
            bool bspot1=false;
            bool bspot2=false;
            // 2. criterion: linker has to lie in the oriented cone with peak "binding spot location"
            // first and second tria(d vectors (tangent and normal)
            LINALG::Matrix<3,1> firstdir;
            LINALG::Matrix<3,1> bspotvec;
            // retrieve tangential and normal vector from binding spot quaternions
            LINALG::Matrix<3,3> bspottriad;
            // auxiliary variable for storing a triad in quaternion form
            LINALG::Matrix<4, 1> qnode;
            // triad of node on first filament which is affected by the new crosslinker
            for (int l=0; l<4; l++)
              qnode(l) = (*bspottriadscol)[l][bspotcolmap_->GID(tmplid)];
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
              crossbspotdiff(l) = (*crosslinkerpositions_)[l][clID]-(*bspotpositions)[l][bspotcolmap_->GID(tmplid)];
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
              if((int)(*numbond_)[clID]<0.1)
                bspot2=true;
              else //check wether first Bspot lies fullfills orientation criterium
              {       //id of bindingspot that was bound ealier
                int bspotID=0;
                if((int)(*crosslinkerbond_)[0][clID] < -0.9)
                  bspotID=(int)(*crosslinkerbond_)[1][clID];
                else if((int)(*crosslinkerbond_)[1][clID] < -0.9)
                  bspotID=(int)(*crosslinkerbond_)[0][clID];
                else //
                  dserror("Error in crosslinker management and/or search!");
                  //turn difference Vectors direction
                for(int l=0; l<(int)crossbspotdiff.M(); l++)
                 crossbspotdiff(l) = -crossbspotdiff(l);

                for (int l=0; l<4; l++)
                   qnode(l) = (*bspottriadscol)[l][bspotID];
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
                neighbournodes[clID].push_back(tmplid);
            }
          }
          else // only difference criterion applied
            neighbournodes[clID].push_back(tmplid);
        }//end orientation & distance check
      }//end loop over all bspots in octant
      // store local maximal number of LIDs per molecule in order to determine neighbourslid->NumVectors()
      maxneighbourslocal = std::max(maxneighbourslocal, (int)neighbournodes[clID].size());
    }//end loop over all crosslinker in octant
  }//end loop over all octants

  // get global maximal number of LIDs per molecule
  discret_->Comm().MaxAll(&maxneighbourslocal, &maxneighboursglobal, 1);
  if(maxneighboursglobal==0)
    maxneighboursglobal = 1;

  // copy information to Epetra_MultiVector for communication
  Teuchos::RCP<Epetra_MultiVector> neighbourslidtrans = Teuchos::rcp(new Epetra_MultiVector(*transfermap_,maxneighboursglobal, false));
  /* assign "-2" -> 'empty'
   * to be able to determine entries which remain "empty" due to number of LIDs < maxneighboursglobal*/
  neighbourslidtrans->PutScalar(0);
  neighbourslid = Teuchos::rcp(new Epetra_MultiVector(*crosslinkermap_,maxneighboursglobal,true));
  /* here we add 2.0 to every ID that enters neighbourslid avoiding zero-entries resulting from neighboursID=0.
   * Later, after communication we substract 2.0*/
  for(int i=0; i<neighbourslid->MyLength(); i++)
    for(int j=0; j<(int)neighbournodes[i].size(); j++)
      (*neighbourslid)[j][i] = (double)neighbournodes[i][j]+2.0;

  // make information redundant on all Procs
  CommunicateMultiVector(neighbourslidtrans, neighbourslid, true,false,false);
  CommunicateMultiVector(neighbourslidtrans, neighbourslid, false,true);

  //substract 2.0 from all entries in neighbourslid to retrieve right value of ID's.
  //Now, entry=-2 means no neighbour, entry>=0 stands for the neighbour ID' s
  for(int i=0; i<neighbourslid->MyLength(); i++)
    for(int j=0; j<(int)neighbourslid->NumVectors(); j++)
      (*neighbourslid)[j][i]=(*neighbourslid)[j][i]-2.0;
  return;
}


/*----------------------------------------------------------------------*
 | Binding spot / linker search by Binning               mueller (11/13)|
 *----------------------------------------------------------------------*/
void STATMECH::StatMechManager::DetectBindingSpotsBinning(const Teuchos::RCP<Epetra_MultiVector> bspotpositions,
                                                          Teuchos::RCP<Epetra_MultiVector>       bspottriadscol,
                                                          Teuchos::RCP<Epetra_MultiVector>&      neighbourslid)
{
  Teuchos::RCP<STATMECH::SEARCH::BinSearch> binsearch = Teuchos::rcp(new STATMECH::SEARCH::BinSearch(discret_,*periodlength_,*searchres_));

  binsearch->AssignPositionsToBins(bspotpositions);

  // proc-wise neighbor search
  LINALG::Matrix<6,1> rootlimits = binsearch->GetRootLimits();

  Teuchos::RCP<Epetra_MultiVector> crosslinkerpositionstrans = Teuchos::rcp(new Epetra_MultiVector(*transfermap_, 3, true));
  CommunicateMultiVector(crosslinkerpositionstrans, crosslinkerpositions_,true,false,false,true);

  std::vector<std::vector<int> > neighbourbspotstransstd(crosslinkerpositionstrans->MyLength(), std::vector<int>());

  int maxneighbourslocal = -1;

  for(int i=0; i<crosslinkerpositionstrans->MyLength(); i++)
  {
    if((*numbond_)[crosslinkermap_->LID(transfermap_->GID(i))]<1.9)
    {
      std::vector<int> indices(3,-1);
      for(int j=0; j<crosslinkerpositionstrans->NumVectors(); j++)
      {
        indices[j] = (int)std::floor(((*crosslinkerpositionstrans)[j][i]-rootlimits(2*j))/(rootlimits(2*j+1)-rootlimits(2*j))*(double)(*searchres_)[j]);
        if(indices[j]>=(*searchres_)[j] || indices[j]<0)
          dserror("Bin index indices[%i] = %i is outside of the rootlimits! Search resolution: %i", j, indices[j], (*searchres_)[j]);
      }

      int crosscollid = crosslinkermap_->LID(transfermap_->GID(i));
      double rmin = (statmechparams_.get<double>("R_LINK", 0.0)-statmechparams_.get<double>("DeltaR_LINK", 0.0)) / 2.0;
      double rmax = (statmechparams_.get<double>("R_LINK", 0.0)+statmechparams_.get<double>("DeltaR_LINK", 0.0)) / 2.0;
      if ((*numbond_)[crosscollid]>0.9)
      {
        rmin *= 2.0;
        rmax *= 2.0;
      }


      const std::vector<int> bins = binsearch->GetSurroundingBins(indices);

      for(int bin= 0; bin<(int)bins.size(); bin++)
      {
        std::vector<int> bspotsinbin = binsearch->GetBin(bins[bin]).GetBinMembers();

        // loop over binding spots in adjacent bins
        for(int j=0; j<(int)bspotsinbin.size(); j++)
        {
          int bspotlid = bspotcolmap_->LID(bspotsinbin[j]);
          LINALG::Matrix<3, 1> distance;
          for (int k=0; k<(int)distance.M(); k++)
            distance(k) = (*crosslinkerpositionstrans)[k][i]-(*bspotpositions)[k][bspotlid];

          BindingSpotDistanceCriterion(i,crosscollid,bspotlid,rmin,rmax,distance,neighbourbspotstransstd,crosslinkerpositionstrans,bspotpositions,bspottriadscol);
        }
      }
      maxneighbourslocal = std::max(maxneighbourslocal, (int)neighbourbspotstransstd[i].size());
    }
  }
  // communication, gather information into neighbourslid
  int maxneighboursglobal = 0;
  discret_->Comm().MaxAll(&maxneighbourslocal, &maxneighboursglobal, 1);
  if(maxneighboursglobal==0)
    maxneighboursglobal = 1;

  // transfer to communicable format
  Teuchos::RCP<Epetra_MultiVector> neighbourslidtrans = Teuchos::rcp(new Epetra_MultiVector(*transfermap_,maxneighboursglobal));
  neighbourslidtrans->PutScalar(-2.0);
  for(int i=0; i<(int)neighbourslidtrans->MyLength(); i++)
    for(int j=0; j<(int)neighbourbspotstransstd[i].size(); j++)
      (*neighbourslidtrans)[j][i] = neighbourbspotstransstd[i][j];

  neighbourslid = Teuchos::rcp(new Epetra_MultiVector(*crosslinkermap_,maxneighboursglobal));

  CommunicateMultiVector(neighbourslidtrans, neighbourslid,false,true);


  return;
}

/*----------------------------------------------------------------------*
 |Neighbor search binding site distance criterion        mueller (11/13)|
 *----------------------------------------------------------------------*/
void STATMECH::StatMechManager::BindingSpotDistanceCriterion(const int&                             crossrowlid,
                                                             const int&                             crosscollid,
                                                             const int&                             bspotlid,
                                                             const double&                          mindist,
                                                             const double&                          maxdist,
                                                             LINALG::Matrix<3,1>                    distance,
                                                             std::vector<std::vector<int> >&        neighbourbspotstransstd,
                                                             const Teuchos::RCP<Epetra_MultiVector> crosslinkerpositionstrans,
                                                             const Teuchos::RCP<Epetra_MultiVector> bspotpositions,
                                                             const Teuchos::RCP<Epetra_MultiVector> bspottriadscol)
{
  if(distance.Norm2()<maxdist && distance.Norm2()>mindist)
  {
    // further calculations in case of helical binding spot geometry and singly bound crosslinkers
    if(filamentmodel_ == statmech_filament_helical)
    {
      bool bspot1=false;
      bool bspot2=false;
      // 2. criterion: linker has to lie in the oriented cone with peak "binding spot location"
      // first and second tria(d vectors (tangent and normal)
      LINALG::Matrix<3,1> firstdir;
      LINALG::Matrix<3,1> bspotvec;
      // retrieve tangential and normal vector from binding spot quaternions
      LINALG::Matrix<3,3> bspottriad;
      // auxiliary variable for storing a triad in quaternion form
      LINALG::Matrix<4, 1> qnode;
      // triad of node on first filament which is affected by the new crosslinker
      for (int l=0; l<4; l++)
        qnode(l) = (*bspottriadscol)[l][bspotcolmap_->GID(bspotlid)];
      LARGEROTATIONS::quaterniontotriad(qnode, bspottriad);
      for (int l=0; l<(int)bspottriad.M(); l++)
      {
        firstdir(l) = bspottriad(l,0);
        bspotvec(l) = bspottriad(l,1);
      }

      // rotation matrix around tangential vector by given angle
      RotationAroundFixedAxis(firstdir,bspotvec,(*bspotorientations_)[bspotcolmap_->GID(bspotlid)]);

      // linker position
      LINALG::Matrix<3,1> crossbspotdiff;
      for(int l=0; l<(int)crossbspotdiff.M(); l++)
        crossbspotdiff(l) = (*crosslinkerpositionstrans)[l][crossrowlid]-(*bspotpositions)[l][bspotcolmap_->GID(bspotlid)];
      // line parameter of the intersection point of the line through the binding spot with the orientation of the binding spot.
      // a)lambda must be larger than zero in order to lie within the cone
      double lambda = (crossbspotdiff.Dot(bspotvec))/(bspotvec.Dot(bspotvec));
      if(lambda>0)
      {
        distance.Scale(1.0/distance.Norm2()); // normalized vector
        // b) the angle between the binding spot orientation and the vector between bspot and crosslinker must be smaller than PHIBSPOT
        //    Only then does the crosslinker lie within the cone
        if(acos(fabs(bspotvec.Dot(distance))) < statmechparams_.get<double>("PHIBSPOT", 0.524))
          bspot1=true;
        if((int)(*numbond_)[crosscollid]<0.1)
          bspot2=true;
        else //check wether first Bspot lies fullfills orientation criterium
        {       //id of bindingspot that was bound ealier
          int bspotID=0;
          if((int)(*crosslinkerbond_)[0][crosscollid] < -0.9)
            bspotID=(int)(*crosslinkerbond_)[1][crosscollid];
          else if((int)(*crosslinkerbond_)[1][crosscollid] < -0.9)
            bspotID=(int)(*crosslinkerbond_)[0][crosscollid];
          else
          { //
            std::cout<<"0: "<<(*crosslinkerbond_)[0][crosscollid]<<", 1:"<<(*crosslinkerbond_)[1][crosscollid]<<std::endl;
            dserror("Error in crosslinker management and/or search!");
          }
          //turn difference Vectors direction
          for(int l=0; l<(int)crossbspotdiff.M(); l++)
           crossbspotdiff(l) = -crossbspotdiff(l);

          for (int l=0; l<4; l++)
             qnode(l) = (*bspottriadscol)[l][bspotID];
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
            distance.Scale(1.0/distance.Norm2()); // normalized vector
            // b) the angle between the binding spot orientation and the vector between bspot and crosslinker must be smaller than PHIBSPOT
            //    Only then does the crosslinker lie within the cone
            if(acos(fabs(bspotvec.Dot(distance))) < statmechparams_.get<double>("PHIBSPOT", 0.524))
              bspot2=true;
          }
        }
        if(bspot1 && bspot2)
          neighbourbspotstransstd[crossrowlid].push_back(bspotlid);
      }
    }
    else // only difference criterion applied
      neighbourbspotstransstd[crossrowlid].push_back(bspotlid);
  }
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
  // beam elements only. Binding spot triads represent tangential vector for Euler type elements
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

  INPAR::STATMECH::BSSearchType searchtype = DRT::INPUT::IntegralValue<INPAR::STATMECH::BSSearchType>(statmechparams_,"BINDINGSITESEARCH");

  switch(searchtype)
  {
    case INPAR::STATMECH::bsstype_volpart:
    {
      if(statmechparams_.get<int>("SEARCHRES",1)>0)
        PartitioningAndSearch(bspotpositions,bspottriadscol, neighbourslid);
    }
    break;
    case INPAR::STATMECH::bsstype_octree:
      DetectBindingSpotsOctree(bspotpositions, bspottriadscol,neighbourslid);
    break;
    case INPAR::STATMECH::bsstype_binning:
      DetectBindingSpotsBinning(bspotpositions,bspottriadscol,neighbourslid);
    break;
    default: dserror("Unknown binding site search method %d !", searchtype);
  }

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
  int numstdlinkers = 0;
  int numactlinkers = 0;

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
              if((*bspotstatus_)[bspotLID]>-0.1)
                continue;

              // flag indicating loop break after first new bond has been established between i-th crosslink molecule and j-th neighbour node
              bool bondestablished = false;
              // necessary condition to be fulfilled in order to set a crosslinker
              double probability = 0.0;
              // switch between probability for regular inter-filament binding and self-binding crosslinker
              probability = plink;
              // self binding is activated
              if(statmechparams_.get<double>("K_ON_SELF", 0.0)>0.0)
              {
                if((int)(*numbond_)[irandom]==1)
                {
                  int nodelid0 = -1;
                  int nodelid1 = -1;
                  switch(linkermodel_)
                  {
                    case statmech_linker_stdintpol:
                    case statmech_linker_bellseqintpol:
                    case statmech_linker_activeintpol:
                    case statmech_linker_myosinthick:
                    {
                      if(linkermodel_==statmech_linker_myosinthick)
                        std::cout<<"  ===WARNING! Implementation incomplete!!!==="<<std::endl;
                      nodelid1 = discret_->NodeColMap()->LID((int)(*bspot2nodes_)[0][irandom]);
                      for(int k=0; k<crosslinkerbond_->NumVectors(); k++)
                        if((*crosslinkerbond_)[k][irandom]>-0.1)
                        {
                          nodelid0 = discret_->NodeColMap()->LID((int)(*bspot2nodes_)[0][(int)(*crosslinkerbond_)[k][irandom]]);
                          break;
                        }
                    }
                    break;
                    case statmech_linker_std:
                    case statmech_linker_bellseq:
                    case statmech_linker_active:
                    {
                      nodelid1 = bspotLID;
                      for(int k=0; k<crosslinkerbond_->NumVectors(); k++)
                        if((*crosslinkerbond_)[k][irandom]>-0.1)
                        {
                          nodelid0 = discret_->NodeColMap()->LID((int)(*crosslinkerbond_)[k][irandom]);
                          break;
                        }
                    }
                    break;
                    default:
                      dserror("Unknown linker type!");
                  }
                  if((*filamentnumber_)[nodelid0] == (*filamentnumber_)[nodelid1])
                    probability = pself;
                }
              }

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
                  // one bond already exists -> establish a second bond
                  // Potentially, add an element (later)
                  case 1:
                  {
                    //Col map LIDs of nodes to be crosslinked
                    Epetra_SerialDenseMatrix LID(2,1);
                    LID(0,0) = bspotcolmap_->LID((int)(*crosslinkerbond_)[occupied][irandom]);
                    LID(1,0) = bspotLID;

                    // Check if crosslinkers are about to occupy two binding sites on the same filament in case this is NOT wanted
                    if(statmechparams_.get<double>("K_ON_SELF",0.0)==0.0)
                    {
                      bool linkonsamefilament = false;
                      switch(linkermodel_)
                      {
                        case statmech_linker_stdintpol:
                        case statmech_linker_activeintpol:
                        case statmech_linker_bellseqintpol:
                        case statmech_linker_myosinthick:
                          if((*filamentnumber_)[discret_->NodeColMap()->LID((int)(*bspot2nodes_)[0][(int)LID(0,0)])] == (*filamentnumber_)[discret_->NodeColMap()->LID((int)(*bspot2nodes_)[0][(int)LID(1,0)])])
                            linkonsamefilament = true;
                          break;
                        default:
                          if((*filamentnumber_)[(int)LID(0,0)]==(*filamentnumber_)[(int)LID(1,0)])
                            linkonsamefilament = true;
                      }
                      // leave neighourslid-loop
                      if(linkonsamefilament)
                        break;
                    }

                    // casimir-specific handling of crosslinker setting (potentially remove)
                    if(networktype_ == statmech_network_casimir)
                    {
                      if(!SetCrosslinkerLoom(LID, bspotpositions, bspottriadscol))
                        break;
                    }
                    // DBC-specific measures for setting crosslinkers
                    // interpolated crosslinker elements transport mechanical information across the periodic boundary, which
                    // is not intended or sometimes unwanted. Hence, no crosslinks are established across the boundary
                    // Finding such a potential connection results in break -> leaving the current loop step and continuing with next linker
                    if(timen>=actiontime_->at(bctimeindex_))
                    {
                      if(CheckCrossPeriodicBCLink(LID,discol))
                        break;
                    }

                    //direction vector between currently considered two nodes (direction free - occupied)
                    LINALG::Matrix<3,1> direction;
                    for(int j=0;j<bspotpositions->NumVectors();j++)
                      direction(j) = (*bspotpositions)[j][(int)LID(0,0)]-(*bspotpositions)[j][(int)LID(1,0)];

                    // skip binding process, if cycle time is still detached and
                    // skip binding process, if polarity of active linker is not fulfilled
                    bool polaritycriterion = true;
                    if(DRT::INPUT::IntegralValue<int>(statmechparams_,"FILAMENTPOLARITY"))
                    {
                      switch(linkermodel_)
                      {
                        case statmech_linker_active:
                        case statmech_linker_activeintpol:
                        case statmech_linker_myosinthick:
                          polaritycriterion = LinkerPolarityCheckAttach(bspottriadscol, LID, direction);
                          break;
                        default: dserror("Filament polarity not implemented for linkermodel_ %i !", linkermodel_);
                      }
                    }
                    else
                      direction.Scale(1.0/direction.Norm2());

                    if(polaritycriterion)
                    {
                      /* In case of beam contact, evaluate whether a crosslinker would intersect with any exisiting elements
                       * if it were to be set between the two nodes considered*/
                      bool intersection = false;
                      if(beamcmanager!=Teuchos::null)
                      {
                        if(linkermodel_ == statmech_linker_stdintpol ||
                           linkermodel_ == statmech_linker_activeintpol ||
                           linkermodel_ == statmech_linker_bellseqintpol ||
                           linkermodel_ == statmech_linker_myosinthick)
                          dserror("Beam contact not implemented for interpolated elements!");
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
                      if(CheckOrientation(direction,discol,bspottriadscol,LID) && !intersection)
                      {
                        numsetelements++;
                        if(linkermodel_ == statmech_linker_active || linkermodel_ == statmech_linker_activeintpol ||
                           linkermodel_ == statmech_linker_myosinthick || statmechparams_.get<double>("ACTIVELINKERFRACTION",0.0)>0.0)
                        {
                          if((*crosslinkertype_)[irandom]==1.0)
                            numactlinkers++;
                          else
                            numstdlinkers++;
                        }

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
                        }
                        else
                          dserror("LID = %i < 0! Check binding site maps!",bspotLID);

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
  /* note: searchforneighbours_ is not communicated
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
  // the node row map is needed in order to adequately assessing ownership of elements to processors
  const Epetra_Map noderowmap(*(discret_->NodeRowMap()));

  int numsetelementsall = -1;
  discret_->Comm().MaxAll(&numsetelements, &numsetelementsall, 1);

  if(numsetelementsall>0)
  {
    // positions from displacement vector
    Teuchos::RCP<Epetra_MultiVector> nodalpositions = Teuchos::null;
    // rotations from displacement vector
    Teuchos::RCP<Epetra_MultiVector> nodalrotations = Teuchos::null;
    // NODAL quaternions from filament elements
    Teuchos::RCP<Epetra_MultiVector> nodalquaternions = Teuchos::null;

    if((linkermodel_ == statmech_linker_stdintpol ||
       linkermodel_ == statmech_linker_activeintpol ||
       linkermodel_ == statmech_linker_bellseqintpol ||
       linkermodel_ == statmech_linker_myosinthick) && statmechparams_.get<double>("ILINK",0.0)>0.0)
    {

      if(linkermodel_==statmech_linker_myosinthick)
        nodalpositions = Teuchos::rcp(new Epetra_MultiVector(*(discret_->NodeColMap()),3));
      nodalrotations = Teuchos::rcp(new Epetra_MultiVector(*(discret_->NodeColMap()),3));
      GetNodalBindingSpotPositionsFromDisVec(discol, nodalpositions, nodalrotations);
      nodalquaternions = Teuchos::rcp(new Epetra_MultiVector(*(discret_->NodeColMap()),4));
      GetElementBindingSpotTriads(nodalquaternions);
    }
    else if(linkermodel_ == statmech_linker_std  && statmechparams_.get<double>("ILINK",0.0)==0.0 && statmechparams_.get<double>("ALINK",0.0)==0.0)
    {
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

        switch(linkermodel_)
        {
          case statmech_linker_stdintpol:
          case statmech_linker_activeintpol:
          case statmech_linker_bellseqintpol:
          case statmech_linker_myosinthick:
          {
            globalnodeids = Teuchos::rcp(new std::vector<int>(4,0));
            for(int l=0; l<(int)bspot2nodes_->NumVectors(); l++)
            {
              (*globalnodeids)[l]=(int)(*bspot2nodes_)[l][bspotcolmap_->LID(bspotgid[0])];
              (*globalnodeids)[l+2]=(int)(*bspot2nodes_)[l][bspotcolmap_->LID(bspotgid[1])];
            }
            if(linkermodel_==statmech_linker_myosinthick)
              CreateTransverseNodePairs(globalnodeids, nodalpositions);
          }
          break;
          case statmech_linker_std:
          case statmech_linker_bellseq:
          case statmech_linker_active:
          {
            globalnodeids = Teuchos::rcp(new std::vector<int>(2,0));
            (*globalnodeids)[0] = bspotgid[0];
            (*globalnodeids)[1] = bspotgid[1];
          }
          break;
          default: dserror("Unknown linker type %i!",linkermodel_);
        }

        // calculate element GID
        int newcrosslinkerGID = GenerateNewLinkerGID(bspotgid);

        /* Create mapping from crosslink molecule to crosslinker element GID
         * Note: No need for the usual procedure of exporting and reimporting to make things redundant
         * because info IS already redundant by design here.*/
        (*crosslink2element_)[i] = newcrosslinkerGID;

        //save positions of nodes between which a crosslinker has to be established in variables xrefe and rotrefe:
        std::vector<double> rotrefe(6);
        std::vector<double> xrefe(6);
        /* resize in case of interpolated crosslinker element. No need for interpolated Beam3eb element, because rotational
           d.o.f. does not exist */
        if (linkermodel_ == statmech_linker_stdintpol ||
            linkermodel_ == statmech_linker_activeintpol ||
            linkermodel_ == statmech_linker_bellseqintpol ||
            linkermodel_ == statmech_linker_myosinthick)
          rotrefe.resize(12);

        // all beam elements
        if(statmechparams_.get<double>("ILINK",0.0)>0.0)
        {
          for(int k=0; k<3; k++)
          {
            xrefe[k]   = (*bspotpositions)[k][bspotcolmap.LID(bspotgid.at(0))];
            xrefe[k+3] = (*bspotpositions)[k][bspotcolmap.LID(bspotgid.at(1))];

            //set nodal rotations (not true ones, only those given in the displacement vector)
            switch(linkermodel_)
            {
              case statmech_linker_stdintpol:
              case statmech_linker_activeintpol:
              case statmech_linker_bellseqintpol:
              case statmech_linker_myosinthick:
              {
                rotrefe[k]   = (*nodalrotations)[k][bspotcolmap.LID((*globalnodeids)[0])];
                rotrefe[k+3] = (*nodalrotations)[k][bspotcolmap.LID((*globalnodeids)[1])];
                rotrefe[k+6] = (*nodalrotations)[k][bspotcolmap.LID((*globalnodeids)[2])];
                rotrefe[k+9] = (*nodalrotations)[k][bspotcolmap.LID((*globalnodeids)[3])];
              }
              break;
              case statmech_linker_std:
              case statmech_linker_bellseq:
              case statmech_linker_active:
              {
                rotrefe[k] = (*bspotrotations)[k][bspotcolmap.LID((*globalnodeids)[0])];
                rotrefe[k+3] = (*bspotrotations)[k][bspotcolmap.LID((*globalnodeids)[1])];
              }
              break;
              default: dserror("Unknown linkertype_");
            }
          }
        }
        else // truss elements
        {
          for(int k=0; k<3; k++)
          {
            xrefe[k]   = (*bspotpositions)[k][bspotcolmap.LID(bspotgid.at(0))];
            xrefe[k+3] = (*bspotpositions)[k][bspotcolmap.LID(bspotgid.at(1))];
            // Here rotation does not actually mean rotational degrees of freedom, but tangential dof
            rotrefe[k] = (*bspotrotations)[k][bspotcolmap.LID((*globalnodeids)[0])];
            rotrefe[k+3] = (*bspotrotations)[k][bspotcolmap.LID((*globalnodeids)[1])];
          }
        }


        /*a crosslinker is added on each processor which is row node owner of at least one of its nodes;
         *the processor which is row map owner of the node with the larger GID will be the owner of the new crosslinker element; on the other processors the new
         *crosslinker element will be present as a ghost element only*/
        bool hasrownode = CheckIfRowNodes(noderowmap,globalnodeids);

        if(hasrownode)
          AddNewCrosslinkerElement(newcrosslinkerGID, globalnodeids,bspotgid, xrefe,rotrefe,*discret_,nodalquaternions);

        // add all new elements to contact discretization on all Procs
        if(beamcmanager!=Teuchos::null)
          AddNewCrosslinkerElement(newcrosslinkerGID, globalnodeids,bspotgid, xrefe,rotrefe,beamcmanager->BTSolDiscret(),nodalquaternions);

        // add two more beam3 elements for myosin thick filament
        if(linkermodel_==statmech_linker_myosinthick)
          AddSupportElements(noderowmap,i,beamcmanager,nodalpositions,globalnodeids,bspotgid,xrefe,rotrefe);
      }
    }

    // synchronization for problem discretization
    discret_->CheckFilledGlobally();
    discret_->FillComplete(true, false, false);

    // synchronization for contact discretization
    if(beamcmanager!=Teuchos::null)
    {
      beamcmanager->BTSolDiscret().CheckFilledGlobally();
      beamcmanager->BTSolDiscret().FillComplete(true, false, false);
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
  {
    std::cout<<"\n\n"<<numsetelements<<" crosslinker element(s) added!"<<std::endl;
    if(numstdlinkers+numactlinkers>0)
    {
      std::cout<<" - "<<numstdlinkers<<" standard linkers"<<std::endl;
      std::cout<<" - "<<numactlinkers<<" active   linkers"<<std::endl;
    }
  }
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
{ // If "ILINK>0", linker element is a beam, if "ILINK==0" linker element is truss
  if(statmechparams_.get<double>("ILINK",0.0) > 0.0)
  {
  // get the nodes from the discretization (redundant on both the problem as well as the contact discretization
  // crosslinker is a 4 noded beam Element
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
        // globalnodeids[0] is the larger of the two node GIDs
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
        // set interpolated position to mid positions -> xi = 0.0
        if(linkermodel_==statmech_linker_myosinthick)
          newcrosslinker->SetBindingPosition(0.0,0.0);
        else
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
        break;
    }
  }
  // For Truss linkers
  else if(statmechparams_.get<double>("ILINK",0.0) == 0.0 && statmechparams_.get<double>("ALINK",0.0) != 0.0)
  {
    switch((int)globalnodeids->size())
    {
     case 2:
     {
       Teuchos::RCP<DRT::ELEMENTS::Truss3> newcrosslinker = Teuchos::rcp(new DRT::ELEMENTS::Truss3(crossgid, (mydiscret.gNode(globalnodeids->at(0)))->Owner() ) );
       DRT::Node* nodes[2] = {mydiscret.gNode( globalnodeids->at(0) ), mydiscret.gNode( globalnodeids->at(1)) };

       newcrosslinker->SetNodeIds(2,&(*globalnodeids)[0]);
       newcrosslinker->BuildNodalPointers(&nodes[0]);

       //setting up crosslinker element parameters
       newcrosslinker ->SetCrossSec(statmechparams_.get<double>("ALINK",0.0));
       newcrosslinker->SetMaterial(2);

       //correct reference configuration data is computed for the new crosslinker element;
       newcrosslinker->SetUpReferenceGeometry(xrefe,rotrefe);

       //add element to discretization
       if(!addinitlinks)
         addedelements_.push_back(crossgid);
       mydiscret.AddElement(newcrosslinker);
     }
     break;
     case 4:
     {
       Teuchos::RCP<DRT::ELEMENTS::Truss3CL> newcrosslinker = Teuchos::rcp(new DRT::ELEMENTS::Truss3CL(crossgid, (mydiscret.gNode(globalnodeids->at(0)))->Owner() ) );
       DRT::Node* nodes[4] = {mydiscret.gNode( globalnodeids->at(0) ),
                              mydiscret.gNode( globalnodeids->at(1) ),
                              mydiscret.gNode( globalnodeids->at(2) ),
                              mydiscret.gNode( globalnodeids->at(3) )};

       newcrosslinker->SetNodeIds((int)globalnodeids->size(),&(*globalnodeids)[0]);
       newcrosslinker->BuildNodalPointers(&nodes[0]);

       //setting up crosslinker element parameters
       newcrosslinker ->SetCrossSec(statmechparams_.get<double>("ALINK",0.0));
       newcrosslinker->SetMaterial(2);
       newcrosslinker->SetBindingPosition((double)(*bspotxi_)[bspotgid[0]],(double)(*bspotxi_)[bspotgid[1]]);

       //correct reference configuration data is computed for the new crosslinker element;
       newcrosslinker->SetUpReferenceGeometry(xrefe);

       //add element to discretization
       if(!addinitlinks)
         addedelements_.push_back(crossgid);
       mydiscret.AddElement(newcrosslinker);
     }
     break;
     default:
       dserror("Received %i global node IDs! Implemented: 2 or 4!", (int)globalnodeids->size());
       break;
    }
  }
  // For Spring linkers
  else if(statmechparams_.get<double>("ILINK",0.0) == 0.0 && statmechparams_.get<double>("ALINK",0.0) == 0.0)
  {
    switch((int)globalnodeids->size())
    {
     case 2:
     {
       Teuchos::RCP<DRT::ELEMENTS::Spring3> newcrosslinker = Teuchos::rcp(new DRT::ELEMENTS::Spring3(crossgid, (mydiscret.gNode(globalnodeids->at(0)))->Owner() ) );
       DRT::Node* nodes[2] = {mydiscret.gNode( globalnodeids->at(0) ), mydiscret.gNode( globalnodeids->at(1)) };

       newcrosslinker->SetNodeIds(2,&(*globalnodeids)[0]);
       newcrosslinker->BuildNodalPointers(&nodes[0]);

       //setting up crosslinker element parameters
       newcrosslinker->SetMaterial(2);

       //correct reference configuration data is computed for the new crosslinker element;

       DRT::Element * element = discret_->gElement(discret_->lRowElement(0)->Id());
       if (element->ElementType().Name()=="Beam3rType")
       {
         // set initial triads/quaternions
         if(nodalquaternions==Teuchos::null)
           dserror("No nodal quaternions delivered to this method!");

         std::vector<LINALG::Matrix<4,1> > nodequat((int)globalnodeids->size(),LINALG::Matrix<4, 1>(true));
         for(int i=0; i<(int)nodequat.size(); i++)
           for(int j=0; j<(int)nodequat[i].M(); j++)
             (nodequat[i])(j) =(*nodalquaternions)[j][nodalquaternions->Map().LID((*globalnodeids)[i])];

         newcrosslinker->SetInitialTangents(nodequat);
         newcrosslinker->SetUpReferenceGeometry(xrefe,rotrefe,false,true);
       }
       else
         newcrosslinker->SetUpReferenceGeometry(xrefe,rotrefe);


       //add element to discretization
       if(!addinitlinks)
         addedelements_.push_back(crossgid);
       mydiscret.AddElement(newcrosslinker);
     }
     break;
     default:
       dserror("Received %i global node IDs! Implemented: 2 or 4!", (int)globalnodeids->size());
       break;
    }
  }

  return;
} //void StatMechManager::AddNewCrosslinkerElement()


/*----------------------------------------------------------------------*
 | Add supporting elements (myosin thick filament)       mueller (01/14)|
 *----------------------------------------------------------------------*/
void STATMECH::StatMechManager::AddSupportElements(const Epetra_Map&                    noderowmap,
                                                   const int&                           crosslid,
                                                   Teuchos::RCP<CONTACT::Beam3cmanager> beamcmanager,
                                                   Teuchos::RCP<Epetra_MultiVector>     nodalpositions,
                                                   Teuchos::RCP<std::vector<int> >      globalnodeids,
                                                   const std::vector<int>&              bspotgid,
                                                   const std::vector<double>&           xrefe,
                                                   const std::vector<double>&           rotrefe)
{
  if(linkermodel_==statmech_linker_myosinthick)
  {
    // element loop
    for(int i=0; i<2; i++)
    {
      Teuchos::RCP<std::vector<int> > addednodeids = Teuchos::rcp(new std::vector<int>(2,-1));
      (*addednodeids)[0] = globalnodeids->at(2*i);
      (*addednodeids)[1] = globalnodeids->at(2*i+1);

      int additionalgid = GenerateNewLinkerGID(*addednodeids);

      // log additional element IDs for handling later on
      (*additionalcross2ele_)[i][crosslid] = additionalgid;

      // new reference positions
      std::vector<double> addedxrefe(6);
      std::vector<double> addedrotrefe(6);

      //positions
      for(int j=0; j<3; j++)
      {
        addedxrefe[j] = (*nodalpositions)[j][nodalpositions->Map().LID((*addednodeids)[0])];
        addedxrefe[j+3] = (*nodalpositions)[j][nodalpositions->Map().LID((*addednodeids)[1])];
      }
      // rotations
      for(int j=0; j<(int)addedrotrefe.size(); j++)
        addedrotrefe[j] = rotrefe[6*i+j];

      bool hasrownode = CheckIfRowNodes(noderowmap,addednodeids);

      if(hasrownode)
        AddNewCrosslinkerElement(additionalgid, addednodeids,bspotgid, addedxrefe,addedrotrefe,*discret_);

      // add all new elements to contact discretization on all Procs
      if(beamcmanager!=Teuchos::null)
        AddNewCrosslinkerElement(additionalgid, addednodeids,bspotgid, addedxrefe,addedrotrefe,beamcmanager->BTSolDiscret());
    }
  }
  return;
}

/*----------------------------------------------------------------------*
 | Generate a new crosslinker element GID                mueller (01/14)|
 *----------------------------------------------------------------------*/
int STATMECH::StatMechManager::GenerateNewLinkerGID(std::vector<int>& bspotgid)
{
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
  return newcrosslinkerGID;
}

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
 | Extract 3x3 triad from global quaternion vec  (private)mueller (3/14)|
 *----------------------------------------------------------------------*/
LINALG::Matrix<3,3> STATMECH::StatMechManager::ExtractTriadFromGlobalQuaternionVector(const int&                       lid,
                                                                                      Teuchos::RCP<Epetra_MultiVector> quaternionvec)
{
  // ATTENTION, no scaling of vector to unit length! If desired, that has to be done separately
  // (nodal) quaternion
  LINALG::Matrix<4, 1> q(true);
  LINALG::Matrix<3,3> triad(true);
  for (int l=0; l<4; l++)
    q(l) = (*quaternionvec)[l][lid];
  LARGEROTATIONS::quaterniontotriad(q, triad);

  return triad;
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
  Teuchos::RCP<Epetra_MultiVector> bspottriadscol = Teuchos::null;
  if(statmechparams_.get<double>("ILINK",0.0)>0.0)
  {
    bspottriadscol = Teuchos::rcp(new Epetra_MultiVector(*bspotcolmap_,4,true));
    GetBindingSpotTriads(bspotrotations, bspottriadscol);
  }

  // if off-rate is also dependent on the forces acting within the crosslinker
  if (statmechparams_.get<double>("DELTABELLSEQ", 0.0) != 0.0)// ||
//      linkermodel_ == statmech_linker_active ||
//      linkermodel_ == statmech_linker_activeintpol ||
//      linkermodel_ == statmech_linker_myosinthick)
    ForceDependentOffRate(dt, koff, discol, punlink);

  int numdetachpolarity = 0;
  if(DRT::INPUT::IntegralValue<int>(statmechparams_,"FILAMENTPOLARITY"))
//    LinkerPolarityCheckDetach(punlink, bspotpositions, bspottriadscol);
    numdetachpolarity = LinkerPolarityCheckDetach(punlink, bspotpositions, bspottriadscol);

  // a vector indicating the upcoming deletion of crosslinker elements
  Teuchos::RCP<Epetra_Vector> delcrosselement = Teuchos::rcp(new Epetra_Vector(*crosslinkermap_));
  delcrosselement->PutScalar(-1.0);
  Teuchos::RCP<Epetra_MultiVector> additionaldelele = Teuchos::null;
  if(linkermodel_==statmech_linker_myosinthick)
  {
    additionaldelele = Teuchos::rcp(new Epetra_MultiVector(*crosslinkermap_,2));
    additionaldelele->PutScalar(-1.0);
  }


  // SEARCH
  // search and setup for the deletion of elements is done by Proc 0
  int numdelelements = 0;
  int numstdlinkers = 0;
  int numactlinkers = 0;

  if (!discret_->Comm().MyPID())
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
          // singly bound crosslinker in motility assay binds to the substrate filament --> don't delete this bond
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
          for (int j=0; j<crosslinkerbond_->NumVectors(); j++)
          {
            // random pick of one of the crosslinkerbond entries
            int jrandom = jorder[j];
            // in motility assays we want to cut only the second bond (first bond = substrate bond)
            if(networktype_ == statmech_network_motassay)
              jrandom = 1;
            if ((*uniformgen_)() < (*punlink)[jrandom][irandom])
            {
              // after start of DBC application, (in B3CL case, do NOT delete any initial cross boundary elements anymore. Leads to the reverse effect of network softening (don't want that)
              if(timen>=actiontime_->at(bctimeindex_))
              {
                LINALG::SerialDenseMatrix LID(2,1);
                LID(0,0) = bspotcolmap_->LID((int)(*crosslinkerbond_)[0][irandom]);
                LID(1,0) = bspotcolmap_->LID((int)(*crosslinkerbond_)[1][irandom]);
                if(CheckCrossPeriodicBCLink(LID,discol))
                  break;
              }
              (*numbond_)[irandom] = 1.0;

              // an actual crosslinker element exists
              numdelelements++;

              if(statmechparams_.get<double>("ACTIVELINKERFRACTION",0.0)>0.0)
              {
                if((*crosslinkertype_)[irandom]==1.0)
                  numactlinkers++;
                else
                  numstdlinkers++;
              }

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
              // delete additional elements for myosin thick filaments
              if(linkermodel_ == statmech_linker_myosinthick)
              {
                (*additionaldelele)[0][irandom] = (*additionalcross2ele_)[0][irandom];
                (*additionaldelele)[1][irandom] = (*additionalcross2ele_)[1][irandom];
                (*additionalcross2ele_)[0][irandom] = -1.0;
                (*additionalcross2ele_)[1][irandom] = -1.0;
              }

              // vector updates
              (*crosslink2element_)[irandom] = -1.0;
              (*bspotstatus_)[bspotLID] = -1.0;
              (*crosslinkerbond_)[jrandom][irandom] = -1.0;

              // reset linker to long conformation upon unbinding and reset the cycle time
              if(linkermodel_ == statmech_linker_active ||
                 linkermodel_ == statmech_linker_activeintpol ||
                 linkermodel_ == statmech_linker_myosinthick)
              {
                if(crosslinkertype_==Teuchos::null) // all linkers are active
                {
                  (*crosslinkeractlength_)[irandom] = 0.0;
                  (*crosslinkeractcycletime_)[irandom]=0.0;
                }
                else // a fraction of the linkers are active
                {
                  if((*crosslinkertype_)[irandom]==1.0)
                  {
                    (*crosslinkeractlength_)[irandom] = 0.0;
                    (*crosslinkeractcycletime_)[irandom]=0.0;
                  }
                }
              }
              for (int k=0; k<crosslinkerbond_->NumVectors(); k++)
              {
                if ((*crosslinkerbond_)[k][irandom]>-0.9)
                {
                  LINALG::SerialDenseMatrix LID(1, 1, true);
                  LID(0,0) = bspotcolmap_->LID((int)(*crosslinkerbond_)[k][irandom]);
                  CrosslinkerIntermediateUpdate(*bspotpositions, LID, irandom);
                  break;
                }
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
  // note: searchforneighbours_ is not communicated
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
  if(linkermodel_==statmech_linker_active ||
     linkermodel_ == statmech_linker_activeintpol ||
     linkermodel_ == statmech_linker_myosinthick)
  {
    Teuchos::RCP<Epetra_Vector> crosslinkeractcycletimetrans = Teuchos::rcp(new Epetra_Vector(*transfermap_,true));
    CommunicateVector(crosslinkeractcycletimetrans,crosslinkeractcycletime_);
    if(linkermodel_==statmech_linker_myosinthick)
    {
      Teuchos::RCP<Epetra_MultiVector> additionalcross2eletrans = Teuchos::rcp(new Epetra_MultiVector(*transfermap_, 2, true));
      Teuchos::RCP<Epetra_MultiVector> additionaldeleletrans = Teuchos::rcp(new Epetra_MultiVector(*transfermap_, 2, true));
      CommunicateMultiVector(additionaldeleletrans, additionaldelele);
    }
  }

  // DELETION OF ELEMENTS
  RemoveCrosslinkerElements(*discret_,*delcrosselement,&deletedelements_);
  if(linkermodel_ == statmech_linker_myosinthick)
    RemoveMultipleCrosslinkerElements(*discret_,*additionaldelele, &deletedelements_);
  // contact elements
  if(beamcmanager!=Teuchos::null)
  {
    RemoveCrosslinkerElements(beamcmanager->BTSolDiscret(),*delcrosselement,&deletedcelements_);
    if(linkermodel_ == statmech_linker_myosinthick)
      RemoveMultipleCrosslinkerElements(beamcmanager->BTSolDiscret(),*additionaldelele, &deletedelements_);
  }

  // adjust crosslinker-element mapping vector


  if(!discret_->Comm().MyPID() && printscreen)
  {
    std::cout<<numdelelements<<" crosslinker element(s) deleted!"<<std::endl;
    if(numstdlinkers+numactlinkers>0)
    {
      std::cout<<" - "<<numstdlinkers<<" standard linkers"<<std::endl;
      std::cout<<" - "<<numactlinkers<<" active   linkers"<<std::endl;
    }
    if((linkermodel_==statmech_linker_active || linkermodel_ == statmech_linker_activeintpol || linkermodel_ == statmech_linker_myosinthick) && numdetachpolarity>0)
      std::cout<<"  - thereof "<<numdetachpolarity<<" due to filament polarity!"<<std::endl;
  }
  return;
} //StatMechManager::SearchAndDeleteCrosslinkers()

/*----------------------------------------------------------------------*
 | (private) remove crosslinker element                    mueller 1/11 |
 *----------------------------------------------------------------------*/
void STATMECH::StatMechManager::RemoveCrosslinkerElements(DRT::Discretization&             mydiscret,
                                                          Epetra_Vector&                   delcrosselement,
                                                          std::vector<std::vector<char> >* deletedelements)
{
  unsigned startsize = deletedelements->size();
  std::vector<DRT::PackBuffer> vdata( startsize );

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
 | (private) remove crosslinker element (MultiVector)      mueller 1/14 |
 *----------------------------------------------------------------------*/
void STATMECH::StatMechManager::RemoveMultipleCrosslinkerElements(DRT::Discretization&             mydiscret,
                                                                  const Epetra_MultiVector&        delcrosselement,
                                                                  std::vector<std::vector<char> >* deletedelements)
{
  unsigned startsize = deletedelements->size();
  std::vector<DRT::PackBuffer> vdata( startsize );

  for(int h=0; h<delcrosselement.NumVectors(); h++)
  {
    for (int i=0; i<delcrosselement.MyLength(); i++)
    {
      if (mydiscret.HaveGlobalElement((int)delcrosselement[h][i]))
      {
        //save the element by packing before elimination to make it restorable in case that needed
        vdata.push_back( DRT::PackBuffer() );
        mydiscret.gElement((int)delcrosselement[h][i])->Pack( vdata.back() );
      }
    }
    for ( unsigned i=startsize; i<vdata.size(); ++i )
      vdata[i].StartPacking();
    for (int i=0; i<delcrosselement.MyLength(); i++)
    {
      if (mydiscret.HaveGlobalElement((int)delcrosselement[h][i]))
      {
        //save the element by packing before elimination to make it restorable in case that needed
        deletedelements->push_back( std::vector<char>() );
        mydiscret.gElement((int)delcrosselement[h][i])->Pack( vdata[deletedelements->size()-1] );
        mydiscret.DeleteElement( (int)delcrosselement[h][i] );
      }
    }
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
 | (private) Create node pairs for interpolated elements, that do not   |
 | belong to one and the same filament (element)          mueller 01/14 |
 *----------------------------------------------------------------------*/
void STATMECH::StatMechManager::CreateTransverseNodePairs(Teuchos::RCP<std::vector<int> >  globalnodeids,
                                                          Teuchos::RCP<Epetra_MultiVector> nodalpositions)
{
  /* When interpolating not between the nodes of ONE filament element, but between two nodes, that belong to separate elements,
   * this methods comes in handy. The original node order [1 2 3 4] will be shuffled to either [1 3 2 4] or [1 4 2 3] in accordance to
   * the mid position Xbar of the biperpendicular connection of internodal directions v and w
   */

  // distance from calucated mid positions to cog
  std::vector<double> distances;
  distances.clear();

  // COG of the four nodes
  LINALG::Matrix<3,1> cog(true);

  // N1 = first node of element A (first entry in globalnodeids)
  LINALG::Matrix<3,1> Ng(true);
  // N2 = second node of element A (second entry in globalnodeids)
  LINALG::Matrix<3,1> Nh(true);
  // direction nodes 1 3, 2 4
  LINALG::Matrix<3,1> v(true);
  LINALG::Matrix<3,1> w(true);

  LINALG::Matrix<3,1> disttocog(true);

  // calculate COG
  for(int i=0; i<(int)globalnodeids->size(); i++)
    for(int j=0; j<(int)cog.M(); j++)
      cog(j) += (*nodalpositions)[j][nodalpositions->Map().LID(globalnodeids->at(i))];
  cog.Scale(1.0/(double)globalnodeids->size());

  // nodal positions
  for(int i=0; i<(int)Ng.M(); i++)
  {
    Ng(i) = (*nodalpositions)[i][nodalpositions->Map().LID(globalnodeids->at(0))];
    Nh(i) = (*nodalpositions)[i][nodalpositions->Map().LID(globalnodeids->at(1))];
  }

  // permutation vector [2 3 3 2]
  std::vector<int> permut(4,3);
  permut.at(0) = 2;
  permut.at(3) = 2;

  int h_par = -1;

  // look at both possible permutations 1 3, 2 4 and 1 4, 2 3 in order to determine which
  for(int h=0; h<2; h++)
  {
    // get direction vectors
    for(int i=0; i<(int)v.M(); i++)
    {
      v(i) = (*nodalpositions)[i][nodalpositions->Map().LID(globalnodeids->at(permut[2*h]))] - Ng(i);
      w(i) = (*nodalpositions)[i][nodalpositions->Map().LID(globalnodeids->at(permut[2*h+1]))] - Nh(i);
    }

    // n = cross(v,w)
    LINALG::Matrix<3,1> n(true);
    n(0) = v(1)*w(2) - v(2)*w(1);
    n(1) = v(2)*w(0) - v(0)*w(2);
    n(2) = v(0)*w(1) - v(1)*w(0);

    // non-parallel case
    if(n.Norm2()>1e-10)
    {
      n.Scale(1.0/n.Norm2());

      // ng = cross(v,n), nh = cross(w,n)
      LINALG::Matrix<3,1> ng(true);
      LINALG::Matrix<3,1> nh(true);
      ng(0) = v(1)*n(2) - v(2)*n(1);
      ng(1) = v(2)*n(0) - v(0)*n(2);
      ng(2) = v(0)*n(1) - v(1)*n(0);
      nh(0) = w(1)*n(2) - w(2)*n(1);
      nh(1) = w(2)*n(0) - w(0)*n(2);
      nh(2) = w(0)*n(1) - w(1)*n(0);

      // intersection point (h,Eg)
      LINALG::Matrix<3,1> DeltaN(Ng);
      DeltaN -= Nh;
      LINALG::Matrix<3,1> Xh(w);
      Xh.Scale(DeltaN.Dot(ng)/w.Dot(ng));
      Xh += Nh;
      // intersection point (g,Eh)
      DeltaN.Scale(-1.0);
      LINALG::Matrix<3,1> Xg(v);
      Xg.Scale(DeltaN.Dot(nh)/v.Dot(nh));
      Xg += Ng;
      // mid position between intersection points
      LINALG::Matrix<3,1> Xbar(Xg);
      Xbar += Xh;
      Xbar.Scale(0.5);

      LINALG::Matrix<3,1> disttocog(cog);
      disttocog -= Xbar;

      distances.push_back(disttocog.Norm2());
    }
    else // parallelity = 2D = correct set
    {
      h_par = h;
      distances.push_back(0.0);
    }
  }

  std::vector<int> nodeids(*globalnodeids);

  // shuffling node order, non-square case
  if(h_par<0)
  {
    // included the equality case (occurs in special case of exactly same node pair distances 1 3, 2 4, rotated by pi/2)
    // In this case, both alternatives are equivalent also mechanically
    if(distances[0]>=distances[1])
    {
      // node 3 becomes 2 and vice versa
      globalnodeids->at(1) = nodeids.at(2);
      globalnodeids->at(2) = nodeids.at(1);
    }
    else if(distances[1]>distances[0])
    {
      globalnodeids->at(1) = nodeids.at(3);
      globalnodeids->at(2) = nodeids.at(1);
      globalnodeids->at(3) = nodeids.at(2);
    }
  }
  else // square case
  {
    if(h_par==0)
    {
      // node 3 becomes 2 and vice versa
      globalnodeids->at(1) = nodeids.at(2);
      globalnodeids->at(2) = nodeids.at(1);
    }
    else
    {
      globalnodeids->at(1) = nodeids.at(3);
      globalnodeids->at(2) = nodeids.at(1);
      globalnodeids->at(3) = nodeids.at(2);
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 | (private) map element<->crosslinker                    mueller 09/12 |
 *----------------------------------------------------------------------*/
void STATMECH::StatMechManager::ElementToCrosslinkMapping(Teuchos::RCP<Epetra_Vector>& element2crosslink)
{
  // update element2crosslink_ for correct mapping. Reason: Before coming here, elements might have been deleted in
  // SearchAndDeleteCrosslinkers(). The ElementRowMap might have changed -> Rebuild
  if(element2crosslink!=Teuchos::null)
  {
    element2crosslink = Teuchos::rcp(new Epetra_Vector(*(discret_->ElementColMap())));
    element2crosslink->PutScalar(-1.0);

    for(int i=0; i<crosslinkermap_->NumMyElements(); i++)
    {
      if((*crosslink2element_)[i]>-0.9) // there exists a linker element
      {
        // reverse mapping
        int elelid = discret_->ElementColMap()->LID((int)(*crosslink2element_)[i]);
        if(elelid>-0.9)
          (*element2crosslink)[elelid] = i;
      }
    }
  }
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
  if(statmechparams_.get<double>("DELTABELLSEQ", 0.0)!=0.0)
  {
    // thermal energy
    double kt = statmechparams_.get<double>("KT",0.00404531);
    // characteristic vond length (=nodereldis)
    // from B. Gui and W. Guilford: Mechanics of actomyosin bonds in different nucleotide states are tuned to muscle contraction
    // Fig 2: slip pathway, ADP, delta = 0.0004;
    // Note: delta<0 -> catch bond, delta>0 -> bond-weakening
    double delta = statmechparams_.get<double>("DELTABELLSEQ", 0.0);

    // update element2crosslink_
    ElementToCrosslinkMapping(element2crosslink_);

    for(int i=0; i<discret_->ElementColMap()->NumMyElements(); i++)
    {
      if((*element2crosslink_)[i]>-0.9)
      {
        DRT::Element* crosslinker = discret_->lColElement(i);
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
          force = (dynamic_cast<DRT::ELEMENTS::Beam3*>(crosslinker))->InternalForceVector();
          eps = (dynamic_cast<DRT::ELEMENTS::Beam3*>(crosslinker))->EpsilonSgn();
          checkgid = crosslinker->Nodes()[0]->Id();
        }
        else if(eot == DRT::ELEMENTS::BeamCLType::Instance())
        {
          if(crosslinker->NumNode()!=4)
            dserror("Currently only implemented for BEAM3CL with four nodes.");
          force.Resize(crosslinker->NumNode()/2*6);
          force = (dynamic_cast<DRT::ELEMENTS::BeamCL*>(crosslinker))->InternalForceVector();
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
        double koffbspot0 = 0.0;
        double koffbspot1 = 0.0;
        if(kt>1e-16)
        {
          koffbspot0 = koff0 * exp(sgn * Fbspot0 * delta/kt);
          koffbspot1 = koff0 * exp(sgn * Fbspot1 * delta/kt);
        }

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
            (*unbindingprobability_)[1][crosslid] = 1 - exp(-dt * koffbspot0);
            (*unbindingprobability_)[0][crosslid] = 1 - exp(-dt * koffbspot1);
          }
        }
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
  }
  return;
} // StatMechManager::ForceDependentOffRate

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
  if(linkermodel_==statmech_linker_myosinthick)
    dserror("Method not implemented yet for this linker type!");

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
            numdelelements++;
            // enter the crosslinker element ID into the deletion list
            (*delcrosselement)[irandom] = (*crosslink2element_)[irandom];
            (*delcrossmolecules)[irandom] = 1.0;

            ((*bspotstatus_)[(int)(*crosslinkerbond_)[0][irandom]]) = -1.0;
            ((*bspotstatus_)[(int)(*crosslinkerbond_)[1][irandom]]) = -1.0;
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

  //MAXRANDFORCE is a multiple of the standard deviation
  double maxrandforcefac = statmechparams_.get<double>("MAXRANDFORCE",-1.0);

  if(maxrandforcefac==-1.0)
  {
    for (int i=0; i<randomnumbersrow->MyLength(); i++)
      for (int j=0; j<randomnumbersrow->NumVectors(); j++)
        (*randomnumbersrow)[j][i] = standarddeviation*(*normalgen_)() + meanvalue;
  }
  else
  {
    for (int i=0; i<randomnumbersrow->MyLength(); i++)
      for (int j=0; j<randomnumbersrow->NumVectors(); j++)
      {
        (*randomnumbersrow)[j][i] = standarddeviation*(*normalgen_)() + meanvalue;
        if((*randomnumbersrow)[j][i]>maxrandforcefac*standarddeviation + meanvalue)
        {
          (*randomnumbersrow)[j][i]=maxrandforcefac*standarddeviation + meanvalue;
        }
        else if((*randomnumbersrow)[j][i]<-maxrandforcefac*standarddeviation + meanvalue)
        {
          (*randomnumbersrow)[j][i]=-maxrandforcefac*standarddeviation + meanvalue;
        }
      }
  }

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
void STATMECH::StatMechManager::UpdateNumberOfUnconvergedSteps()
{
  unconvergedsteps_++;

  return;
}


/*----------------------------------------------------------------------*
 | checks orientation of crosslinker relative to linked filaments       |
 |                                                  (public) cyron 06/10|
 *----------------------------------------------------------------------*/
bool STATMECH::StatMechManager::CheckOrientation(const LINALG::Matrix<3, 1> direction, const Epetra_Vector& discol, Teuchos::RCP<Epetra_MultiVector> bspottriadscol, const Epetra_SerialDenseMatrix& LID, Teuchos::RCP<double> phifil)
{
  //if orientation is not to be checked explicitly, this function always returns true
  if (statmechparams_.get<double>("ILINK",0.0)>0.0 && (!DRT::INPUT::IntegralValue<int>(statmechparams_, "CHECKORIENT") || bspottriadscol==Teuchos::null))
    return true;
  // Either Spring3 or Truss3 crosslinkers
  else if (statmechparams_.get<double>("ILINK",0.0)==0.0 && (!DRT::INPUT::IntegralValue<int>(statmechparams_, "CHECKORIENT") || bspottriadscol==Teuchos::null))
    return true;

  // check for linkers with both their binding domains on one filament
  if(statmechparams_.get<double>("K_ON_SELF",0.0)>0.0)
  {
    switch(linkermodel_)
    {
      case statmech_linker_stdintpol:
      case statmech_linker_activeintpol:
      case statmech_linker_bellseqintpol:
      case statmech_linker_myosinthick:
        if((*filamentnumber_)[discret_->NodeColMap()->LID((int)(*bspot2nodes_)[0][(int)LID(0,0)])] == (*filamentnumber_)[discret_->NodeColMap()->LID((int)(*bspot2nodes_)[0][(int)LID(1,0)])])
          return true;
        break;
      default:
        if((*filamentnumber_)[(int)LID(0,0)]==(*filamentnumber_)[(int)LID(1,0)])
          return true;
    }
  }

  if(statmechparams_.get<double>("KT") == 0.0)
    dserror("Set KT to a non-zero value.");

  if(LID.M()!=2 && LID.N()!=1)
    dserror("LID has wrong dimensions %d x %d", LID.M(),LID.N());

  //auxiliary variable
  double scalarproduct =0.0;

  //auxiliary variable for storing a triad in quaternion form
  LINALG::Matrix<4, 1> qnode;

  //angle between filament axes at crosslinked points, respectively
  double Phi;

  //Deltaphi = Phi - Phi0, where Phi0 is the angle between crosslinked filaments with zero potential energy (i.e. the most likely one)
  double DeltaPhi;
  if(statmechparams_.get<double>("ILINK",0.0)>0.0 )
  {
    //triads on filaments at the two nodes connected by crosslinkers
    LINALG::Matrix<3, 3> T1;
    LINALG::Matrix<3, 3> T2;

  //triad of binding site on first filament which is affected by the new crosslinker
  for (int j=0; j<4; j++)
    qnode(j) = (*bspottriadscol)[j][(int) LID(0,0)];
  LARGEROTATIONS::quaterniontotriad(qnode, T1);

  //triad of binding site on second filament which is affected by the new crosslinker
  for (int j=0; j<4; j++)
    qnode(j) = (*bspottriadscol)[j][(int) LID(1,0)];
  LARGEROTATIONS::quaterniontotriad(qnode, T2);

  //auxiliary variable
  scalarproduct = T1(0, 0)*T2(0,0) + T1(1,0)*T2(1,0) + T1(2,0)*T2(2,0);
  }

  // Either Spring3 or Truss3 crosslinkers
  else if(statmechparams_.get<double>("ILINK",0.0)==0.0 )
  {
    //triads on filaments at the two nodes connected by crosslinkers
    LINALG::Matrix<3, 1> T(true);
    LINALG::Matrix<6, 1> TrefFil(true);
    LINALG::Matrix<3, 1> T1(true);
    LINALG::Matrix<3, 1> T2(true);

    for (int i=0; i<2; i++)
    {
    //triad of binding site on first filament which is affected by the new crosslinker
      {
        DRT::Node* node = discret_->lColNode((int)LID(i,0));
        std::vector<int> DofNode = discret_->Dof(node);
        //if node has also rotational degrees of freedom
        if (discret_->NumDof(node) == 6)
        {
          //if not, tell the user to use beam3 instead
          for(int j=0; j<3; j++)
          {
            TrefFil(j+3*i) = (*bspottriadscol)[j][(int) LID(i,0)];
            T(j)=TrefFil(j+3*i)+discol[discret_->DofColMap()->LID(DofNode[j+3])];
          }
          if (i==0)
            T1=T;
          else
            T2=T;
        }
      }
    }

    //auxiliary variable
    scalarproduct = T1(0)*T2(0) + T1(1)*T2(1) + T1(2)*T2(2);
  }

  if(statmechparams_.get<double>("ILINK",0.0)>0.0 )
  {
  if(scalarproduct>1.0 && fabs(scalarproduct-1.0)<1e-12)
    scalarproduct = 1.0;
  }
  // Either Spring3 or Truss3 crosslinkers
  if(statmechparams_.get<double>("ILINK",0.0)==0.0 )
  {
    if(scalarproduct>1.0 || scalarproduct<-1.0 )
      return true; // exit
  }

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
      {
        (*crosslinkerpositionstrans)[0][i] += standarddev*(*normalgen_)() + mean;
        (*crosslinkerpositionstrans)[1][i] += standarddev*(*normalgen_)() + mean;
        // in case of two-dimensional linker motion, z-direction is not updated
        if(!DRT::INPUT::IntegralValue<int>(statmechparams_, "PLANELINKERMOTION"))
          (*crosslinkerpositionstrans)[2][i] += standarddev*(*normalgen_)() + mean;
#ifdef DEBUGCOUT
      for (int j=0; j<crosslinkerpositionstrans->NumVectors(); j++)
        (*crosslinkerdeltatrans)[j][i] = standarddev*(*normalgen_)() + mean;
#endif
      }
      break;
      // bonding case 2: crosslink molecule attached to one filament
      case 1:
      {
        int bspotLID = bspotcolmap_->LID(std::max((int)(*crosslinkerbond_)[0][crosslid],(int)(*crosslinkerbond_)[1][crosslid]));
        (*crosslinkerpositionstrans)[0][i] = bspotpositions[0][bspotLID];
        (*crosslinkerpositionstrans)[1][i] = bspotpositions[1][bspotLID];
        if(!DRT::INPUT::IntegralValue<int>(statmechparams_, "PLANELINKERMOTION"))
          (*crosslinkerpositionstrans)[2][i] = bspotpositions[2][bspotLID];
      }
      break;
      // bonding case 3: actual crosslinker has been established
      case 2:
      {
        int largerbspotLID = bspotcolmap_->LID(std::max((int)(*crosslinkerbond_)[0][crosslid],(int)(*crosslinkerbond_)[1][crosslid]));
        int smallerbspotLID = bspotcolmap_->LID(std::min((int)(*crosslinkerbond_)[0][crosslid],(int)(*crosslinkerbond_)[1][crosslid]));
        for (int j=0; j<crosslinkerpositionstrans->NumVectors(); j++)
          (*crosslinkerpositionstrans)[j][i] = (bspotpositions[j][largerbspotLID]+bspotpositions[j][smallerbspotLID])/2.0;
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
      if(DRT::INPUT::IntegralValue<int>(statmechparams_, "PLANELINKERMOTION"))
        deltapos(2) = 0.0;
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
