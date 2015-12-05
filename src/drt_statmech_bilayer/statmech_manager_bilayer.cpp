/*!----------------------------------------------------------------------
\file statmech_manager_bilayer.cpp
\brief management and auxiliary functions for statistical mechanics
       of lipid bilayer

<pre>
Maintainer: Dhrubajyoti Mukherjee
            mukherjee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15276
</pre>

*----------------------------------------------------------------------*/
#include "statmech_manager_bilayer.H"

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
#include "../drt_io/io.H"
#include "../drt_io/io_control.H"
#include "../drt_discsh3/discsh3.H"


/*---------------------------------------------------------------------------*
 |  ctor (public)                                             mukherjee 09/15|
 *---------------------------------------------------------------------------*/
STATMECH::StatMechManagerBilayer::StatMechManagerBilayer(Teuchos::RCP<DRT::Discretization> discret):
statmechBilayerparams_( DRT::Problem::Instance()->StatisticalMechanicsParams() ),
unconvergedsteps_(0),
starttimeoutput_(-1.0),
timeintervalstep_(0),
istart_(0),
numoutputstep_(1),
outputfilenumber_(-1),
discret_(discret)
//facediscret_(Teuchos::null)
{
  Teuchos::ParameterList parameters = DRT::Problem::Instance()->StructuralDynamicParams();

  //initialize random generators
  SeedRandomGenerators(0,-1);

  // retrieve output root path
  BuildStatMechRootPath();

  // retrieve the dimensions of the periodic boundary box and
  // set spatial resolution for search algorithm binding spots x crosslinkers
  InitializeStatMechValues();

  // read points in time from input file where certain actions take place
//  SetStartStep(parameters);


  return;
}// StatMechManagerBilayer::StatMechManagerBilayer


/*---------------------------------------------------------------------------*
 | write special output for statistical mechanics (public)   mukhherjee 09/15|
 *---------------------------------------------------------------------------*/
void STATMECH::StatMechManagerBilayer::Update(const int&                               istep,
                                       const double&                                   timen,
                                       const double&                                   dt,
                                       Epetra_Vector&                                  disrow,
                                       Teuchos::RCP<LINALG::SparseOperator>&           stiff,
                                       bool                                            rebuildoctree,
                                       bool                                            printscreen)
{
#ifdef MEASURETIME
  const double t0 = Teuchos::Time::wallTime();
#endif // #ifdef MEASURETIME
  /* first we modify the displacement vector so that current nodal position at the end of current time step complies with
   * periodic boundary conditions, i.e. no node lies outside a cube of edge length periodlength_*/

    /* note: access by ExtractMyValues requires column map vector, whereas displacements on level of time integration are
     * handled as row map vector*/
    Epetra_Vector discol(*discret_->DofColMap(), true);
    LINALG::Export(disrow, discol);

    if(DRT::INPUT::IntegralValue<int>(statmechBilayerparams_,"BONDFLIP"))
    {
      // Clear NodeIDs of newly formed faces at every time step
      new_faces_NodeID_.clear();
      EvaluateBondFlipsParallel(istep, timen, dt, discol, printscreen);
    }

    // reset thermal energy to new value (simple value change for now, maybe Load Curve later on)
    if(fabs(timen-actiontime_->at(1))<(dt/1e3))
      statmechBilayerparams_.set("KT",statmechBilayerparams_.get<double>("KTACT",statmechBilayerparams_.get<double>("KT",0.0)));

    /*settling administrative stuff in order to make the discretization ready for the next time step:
     * synchronize the Filled() state on all processors after having added or deleted elements by CheckFilledGlobally();
     * then build new element maps and call FillComplete();
     * done here: finally Crs matrices stiff_ has to be deleted completely and made ready
     * for new assembly since their graph was changed*/
    stiff->Reset();

  return;
} // StatMechManagerBilayer::Update()

/*----------------------------------------------------------------------*
 | Update time step size in time integration       (public)mueller 06/12|
 *----------------------------------------------------------------------*/
void STATMECH::StatMechManagerBilayer::UpdateTimeAndStepSize(double& dt,
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
 | update number of unconverged steps              (public) mueller 3/10|
 *----------------------------------------------------------------------*/
void STATMECH::StatMechManagerBilayer::UpdateNumberOfUnconvergedSteps()
{
  unconvergedsteps_++;

  return;
}



/*----------------------------------------------------------------------*
 | Set start step                                (private) mueller 08/13|
 *----------------------------------------------------------------------*/
void STATMECH::StatMechManagerBilayer::SetStartStep(Teuchos::ParameterList& parameters)
{
  // set istart_, the number of steps after which we start writing output
  double starttimeout = statmechBilayerparams_.get<double>("STARTTIMEOUT", 0.0);
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
 | add parameters to given parameter list       (public) mukherjee 09/15|
 *----------------------------------------------------------------------*/
void STATMECH::StatMechManagerBilayer::AddStatMechParamsTo(Teuchos::ParameterList& params, Teuchos::RCP<Epetra_MultiVector> randomnumbers)
{
  params.set("ETA",statmechBilayerparams_.get<double>("ETA",0.0));
  params.set("THERMALBATH",DRT::INPUT::IntegralValue<INPAR::STATMECH::ThermalBathType>(statmechBilayerparams_,"THERMALBATH"));
  params.set<int>("FRICTION_MODEL",DRT::INPUT::IntegralValue<INPAR::STATMECH::FrictionModel>(statmechBilayerparams_,"FRICTION_MODEL"));
  params.set<INPAR::STATMECH::StatOutput>("SPECIAL_OUTPUT", DRT::INPUT::IntegralValue<INPAR::STATMECH::StatOutput>(statmechBilayerparams_,"SPECIAL_OUTPUT"));
  if(randomnumbers!=Teuchos::null)
    params.set("RandomNumbers",randomnumbers);
  // note CURVENUMBER and DBCDISPDIR follow the counting convention 1,2,3,... (not C++). They are internally decremented
  params.set("CURVENUMBER",statmechBilayerparams_.get<int>("CURVENUMBER",-1));
  params.set("DBCDISPDIR",statmechBilayerparams_.get<int>("DBCDISPDIR",-1));
  params.set("STARTTIMEACT",actiontime_->at(bctimeindex_));
  params.set("STARTTIMEACT",actiontime_->at(bctimeindex_));
  params.set("DELTA_T_NEW",timestepsizes_->at(bctimeindex_));

  params.set<std::string>("internalforces","yes");

  // adds information on the used time integration
  params.set<std::string>("DYNAMICTYP", DRT::Problem::Instance()->StructuralDynamicParams().get<std::string>("DYNAMICTYP"));

  return;
}

/*----------------------------------------------------------------------*
 | Set Period Length and Search Resolution       mueller (public)  04/12|
 *----------------------------------------------------------------------*/
void STATMECH::StatMechManagerBilayer::InitializeStatMechValues()
{
  if(!discret_->Comm().MyPID())
    std::cout<<"=================== StatMechManager Init ======================="<<std::endl;

  switch(DRT::INPUT::IntegralValue<INPAR::STATMECH::LinkerModel>(statmechBilayerparams_, "FRICTION_MODEL"))
  {
    case INPAR::STATMECH::frictionmodel_isotropiclumped:
      if(!discret_->Comm().MyPID())
        std::cout<<"  -- friction model: isotropic lumped"<<std::endl;
      break;
    case INPAR::STATMECH::frictionmodel_isotropicconsistent:
      if(!discret_->Comm().MyPID())
        std::cout<<"  -- friction model: isotropic consistent"<<std::endl;
      break;
    case INPAR::STATMECH::frictionmodel_anisotropicconsistent:
      if(!discret_->Comm().MyPID())
        std::cout<<"  -- friction model: anisotropic consistent"<<std::endl;
      break;
    default: dserror("No friction model (i.e. FRICTION_MODEL == none) declared in your input file!");
  }

  periodlength_ = Teuchos::rcp(new std::vector<double>);
  periodlength_->clear();
  {
    std::istringstream PL(Teuchos::getNumericStringParameter(statmechBilayerparams_,"PERIODLENGTH"));
    std::string word;
    char* input;
    while (PL >> word)
      periodlength_->push_back(std::strtod(word.c_str(), &input));
  }

  if((int)periodlength_->size()<3 && (int)periodlength_->size()!=1)
    dserror("You only gave %d values for PERIODLENGTH! Check your input file.", (int)periodlength_->size());
  else if((int)periodlength_->size()==1)
    for(int i=0; i<2; i++)
      periodlength_->push_back(periodlength_->at(0));
  for(int i=0; i<(int)periodlength_->size(); i++)
    if(periodlength_->at(i)<0.0)
      dserror("PERIODLENGTH( %d ) = %4.2f < 0.0 does not make any sense! Check your input file.", i, periodlength_->at(i));

  // read times for actions and the corresponding step sizes from input file
  actiontime_ = Teuchos::rcp(new std::vector<double>);
  actiontime_->clear();
  {
    std::istringstream TIME(Teuchos::getNumericStringParameter(statmechBilayerparams_,"ACTIONTIME"));
    std::string word;
    char* input;
    while (TIME >> word)
      actiontime_->push_back(std::strtod(word.c_str(), &input));
  }
  timestepsizes_ = Teuchos::rcp(new std::vector<double>);
  timestepsizes_->clear();
  {
    std::istringstream DT(Teuchos::getNumericStringParameter(statmechBilayerparams_,"ACTIONDT"));
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

  // set bctimeindex_ (position in actiontime_ where DBCs start being applied
  bctimeindex_ = statmechBilayerparams_.get<int>("BCTIMEINDEX", -1);

  if(!discret_->Comm().MyPID())
  {
    std::cout<<"t_equilib  = t(0) = "<<std::setprecision(10)<<actiontime_->at(0)<<" @ dt = "<<timestepsizes_->at(0)<<std::endl;
    std::cout<<"t_ktswitch = t(1) = "<<std::setprecision(10)<<actiontime_->at(1)<<" @ dt = "<<timestepsizes_->at(1)<<std::endl;
    if(bctimeindex_>-1)
      std::cout<<"t_dbc     = t("<<bctimeindex_<<") ="<<std::setprecision(10)<<actiontime_->at(bctimeindex_)<<" @ dt = "<<timestepsizes_->at(bctimeindex_)<<std::endl;
    if(actiontime_->size()>3)
      std::cout<<"other action times..."<<std::endl;
    for(int i=0; i<(int)actiontime_->size(); i++)
      if(i>1 && i!=bctimeindex_)
        std::cout<<"t("<<i<<")        = "<<std::setprecision(10)<<actiontime_->at(i)<<" @ dt = "<<timestepsizes_->at(i)<<std::endl;
    std::cout<<"================================================================"<<std::endl;

    new_faces_NodeID_.clear();

  }
  return;
}

/*----------------------------------------------------------------------*
 | Evaluate bond-flips                          (public) mukherjee 09/15|
 *----------------------------------------------------------------------*/
void STATMECH::StatMechManagerBilayer::EvaluateBondFlips(const int&         istep,
                                                        const double&      timen,
                                                        const double&      dt,
                                                        Epetra_Vector&     discol,
                                                        bool               printscreen)
{
  facediscret_ = Teuchos::rcp_dynamic_cast<DRT::DiscretizationFaces>(discret_, true);


  // Get number of edges for a particular proc
  int NumFlips=0;              // Number of successful flips
  int NumTrial=0;               // Number of  trials

  std::srand ( unsigned ( std::time(0) ) );
  std::vector<int> RandomFaces;

   // set some values:
   for (int i=1; i<facediscret_->NumMyColFaces(); ++i) RandomFaces.push_back(i); // 1 2 3 4 5 6 7 8 9

   // using built-in random generator:
   std::random_shuffle ( RandomFaces.begin(), RandomFaces.end() );

   int RandEdgeGID = 0;
   for(int i=0; i<facediscret_->NumMyColFaces(); i++)
  {
    // If flip are allowed do a flip
    if(IfBondFlip(RandEdgeGID, discol) && CheckBoundingBoxQuality(RandEdgeGID, discol))
    {
      DoBondFlip(RandEdgeGID);
      NumFlips++;
    }
    else
      RandEdgeGID++;

    NumTrial++;
  }
   std::cout<<"The success rate of number of flips is "<< (double)NumFlips/NumTrial<<std::endl;

  return;
}

/*----------------------------------------------------------------------*
 | Evaluate bond-flip for parallel processors   (public) mukherjee 09/15|
 *----------------------------------------------------------------------*/
void STATMECH::StatMechManagerBilayer::EvaluateBondFlipsParallel(const int&         istep,
                                                        const double&      timen,
                                                        const double&      dt,
                                                        Epetra_Vector&     discol,
                                                        bool               printscreen)
{

  facediscret_ = Teuchos::rcp_dynamic_cast<DRT::DiscretizationFaces>(discret_, true);

  // For GMSH visulization
  std::vector<DRT::Element *> old_faces;
  std::vector<LINALG::Matrix<1,2>  > old_faces_NodeID;
  old_faces.clear();
  {
    for(int i=0; i<facediscret_->NumMyRowFaces(); i++)
    {
      old_faces.push_back(facediscret_->lRowFace(i));
      LINALG::Matrix<1,2> NodeIDs(true);
      NodeIDs(0)=old_faces[i]->NodeIds()[0];
      NodeIDs(1)=old_faces[i]->NodeIds()[1];
      old_faces_NodeID.push_back(NodeIDs);
    }
  }

  int NumFlips=0;              // Number of successful flips
  int NumTrial=0;               // Number of  trials
  int NumTrialAll = 0;
  int NumFlipsAll = 0;

  while(NumFlipsAll<discret_->NumGlobalNodes())
  {
    // Pick random numbers from Number of Edges
    // Set "FIXEDSEED=yes" to generate same set of random numbers for every simulation run

    // loop over the participating processors each of which appends its part of the output to one output file
    for (int proc = 0; proc < discret_->Comm().NumProc(); proc++)
    {
      if (discret_->Comm().MyPID() == proc)
      {
        const int RandEdgeLID =(int)(floor((*uniformgen_)()*(double)(facediscret_->NumMyRowFaces())));

        DRT::Element * ele = facediscret_->lColFace(RandEdgeLID);
        DRT::ELEMENTS::DiscSh3Line* edge = dynamic_cast<DRT::ELEMENTS::DiscSh3Line*>(ele);


        if(IfBondFlipParallel(RandEdgeLID, discol) && CheckBoundingBoxQualityParallel(RandEdgeLID, discol)
            && edge->ParentMasterElement()->Owner()==proc && edge->ParentSlaveElement()->Owner()==proc )
        {
          DoBondFlipParallel(*discret_,*facediscret_,RandEdgeLID);

          NumFlips++;
        }
        NumTrial++;

      }
      discret_->CheckFilledGlobally();
      discret_->FillComplete(true,true,true);
      facediscret_->CheckFilledGlobally();
      facediscret_->FillCompleteFaces(true,true,true,true);
    }

    discret_->Comm().SumAll(&NumFlips, &NumFlipsAll, 1);
    discret_->Comm().SumAll(&NumTrial, &NumTrialAll, 1);

  }
  discret_->CheckFilledGlobally();
  discret_->FillComplete();
  facediscret_->CheckFilledGlobally();
  facediscret_->FillCompleteFaces(true,true,true,true);
  int count=0;
  for(int i=0; i<facediscret_->NumMyRowFaces(); i++)
  {
    bool inoldset=false;
    for(int j=0; j<facediscret_->NumMyRowFaces(); j++)
    {
      if((facediscret_->lRowFace(i)->NodeIds()[0] == old_faces_NodeID[j](0) || facediscret_->lRowFace(i)->NodeIds()[0]==old_faces_NodeID[j](1))
          && (facediscret_->lRowFace(i)->NodeIds()[1] == old_faces_NodeID[j](0) || facediscret_->lRowFace(i)->NodeIds()[1]==old_faces_NodeID[j](1)))
      {
        inoldset=true;
        break;
      }
    }
    if(!inoldset)
    {
      LINALG::Matrix<1,2> NodeIDs(true);
      std::vector<int> NodeIDsNew;
      NodeIDsNew.push_back(facediscret_->lRowFace(i)->NodeIds()[0]);
      NodeIDsNew.push_back(facediscret_->lRowFace(i)->NodeIds()[1]);
      NodeIDs(0)=facediscret_->lRowFace(i)->NodeIds()[0];
      NodeIDs(1)=facediscret_->lRowFace(i)->NodeIds()[1];
      new_faces_NodeID_[count]=NodeIDsNew;
      count++;
    }
  }
  discret_->Comm().Barrier();

  if(!discret_->Comm().MyPID() && printscreen)
  {
    std::cout<<"\n\n"<<"The success rate of number of flips is "<< (double)NumFlipsAll/NumTrialAll<<std::endl;
  }

  return;
}

/*----------------------------------------------------------------------*
 | Chekc if bond-flips allowed                  (public) mukherjee 09/15|
 *----------------------------------------------------------------------*/
bool STATMECH::StatMechManagerBilayer::IfBondFlip(const int&    RandEdgeGID,
                                                  Epetra_Vector&     discol)
{

  bool doflip=false;
  DRT::FaceElement * edge = facediscret_->gFace( RandEdgeGID );

  // Find current length of edge
    LINALG::Matrix <1,6> current_pos(true);
    LINALG::Matrix <1,3> dummy(true);
    for (int dim=0;dim<3;++dim)
    {
      for (int node=0;node<2;++node)
      {
        double referenceposition = ((edge->Nodes())[node])->X()[dim];
        std::vector<int> dofnode = discret_->Dof((edge->Nodes())[node]);
        double displacement = (double)(discol)[discret_->DofColMap()->LID(dofnode[dim])];
        current_pos(3*node+dim) =  referenceposition + displacement;
      }
      dummy(dim)=current_pos(dim)-current_pos(3+dim);
    }

    double length= dummy.Norm2();

    // TODO: import values from input file
    double sigma=0.1;  // diameter of monomer
    // Currently importing dummyvalues
    double l_max=sqrt(100); // Maximum length of tether
    if (length > sigma && length < l_max)
      doflip=true;

    return doflip;
}

/*-------------------------------------------------------------------------------*
 | Chekc if bond-flips allowed for parallel processor    (public) mukherjee 09/15|
 *-------------------------------------------------------------------------------*/
bool STATMECH::StatMechManagerBilayer::IfBondFlipParallel(const int&    RandEdgeLID,
                                                         Epetra_Vector&     discol)
{

  bool doflip=false;
  DRT::Element * ele = facediscret_->lColFace(RandEdgeLID);
  DRT::ELEMENTS::DiscSh3Line* edge = dynamic_cast<DRT::ELEMENTS::DiscSh3Line*>(ele);

  // Find current length of edge
    LINALG::Matrix <1,6> current_pos(true);
    LINALG::Matrix <1,3> dummy(true);
    for (int dim=0;dim<3;++dim)
    {
      for (int node=0;node<2;++node)
      {
        double referenceposition = ((edge->Nodes())[node])->X()[dim];
        std::vector<int> dofnode = discret_->Dof((edge->Nodes())[node]);
        double displacement = (double)(discol)[discret_->DofColMap()->LID(dofnode[dim])];
        current_pos(3*node+dim) =  referenceposition + displacement;
      }
      dummy(dim)=current_pos(dim)-current_pos(3+dim);
    }

    double length= dummy.Norm2();

    // TODO: import values from input file
    double sigma=0.1;  // diameter of monomer
    // Currently importing dummyvalues
    double l_max=sqrt(100); // Maximum length of tether
    if (length > sigma && length < l_max)
      doflip=true;

    return doflip;
}

/*-----------------------------------------------------------------*
 | Check quality of the bounding box consisting of                 |
 | master and slave element                (public) mukherjee 09/15|
 *-----------------------------------------------------------------*/
bool STATMECH::StatMechManagerBilayer::CheckBoundingBoxQuality(const int&    RandEdgeGID,
                                                               Epetra_Vector&     discol)
{
  bool BoundingBoxQuality=true;
  DRT::FaceElement * edge = facediscret_->gFace( RandEdgeGID );

  DRT::Element* pele = edge->ParentMasterElement();  // Master element
  DRT::Element* nele = edge->ParentSlaveElement();   // Slave element

  if (pele == NULL) dserror("pele is NULL");
  if (nele == NULL) dserror("nele is NULL");

  LINALG::Matrix <1,3> current_pos_master(true);
  LINALG::Matrix <1,3> current_pos_slave(true);
  LINALG::Matrix <1,3> edge_vector12(true);
  LINALG::Matrix <1,3> edge_vector21(true);

  // Find the node which doesn't belong to the master element
  for (int i=0; i<nele->NumNode(); i++)
  {
    if(pele->NodeIds()[i]!=edge->NodeIds()[0] && pele->NodeIds()[i]!=edge->NodeIds()[1])
      current_pos_master= GetSpatialPosition((pele->Nodes())[i],discol);
    if(nele->NodeIds()[i]!=edge->NodeIds()[0] && nele->NodeIds()[i]!=edge->NodeIds()[1])
      current_pos_slave= GetSpatialPosition((nele->Nodes())[i],discol);
  }

  // Find current length of edge
  LINALG::Matrix <1,3> edge_node1toNodeMaster(true);
  LINALG::Matrix <1,3> edge_node2toNodeMaster(true);
  LINALG::Matrix <1,3> edge_node1toNodeSlave(true);
  LINALG::Matrix <1,3> edge_node2toNodeSlave(true);

  LINALG::Matrix <1,3> current_pos_edgeNode1 = GetSpatialPosition((edge->Nodes())[0],discol);
  LINALG::Matrix <1,3> current_pos_edgeNode2 = GetSpatialPosition((edge->Nodes())[1],discol);

  for (int dim=0;dim<3;++dim)
  {
    edge_node1toNodeMaster(dim)=current_pos_master(dim)-current_pos_edgeNode1(dim);
    edge_node2toNodeMaster(dim)=current_pos_master(dim)-current_pos_edgeNode2(dim);
    edge_node1toNodeSlave(dim)=current_pos_slave(dim)-current_pos_edgeNode1(dim);
    edge_node2toNodeSlave(dim)=current_pos_slave(dim)-current_pos_edgeNode2(dim);
    edge_vector12(dim)=current_pos_edgeNode2(dim)-current_pos_edgeNode1(dim);
    edge_vector21(dim)=current_pos_edgeNode1(dim)-current_pos_edgeNode2(dim);
  }

  //Nomralize the vectors
  edge_node1toNodeMaster.Scale(1/edge_node1toNodeMaster.Norm2());
  edge_node2toNodeMaster.Scale(1/edge_node2toNodeMaster.Norm2());
  edge_node1toNodeSlave.Scale(1/edge_node1toNodeSlave.Norm2());
  edge_node2toNodeSlave.Scale(1/edge_node2toNodeSlave.Norm2());
  edge_vector12.Scale(1/edge_vector12.Norm2());
  edge_vector21.Scale(1/edge_vector21.Norm2());


  double Theta_edgeNode1AndNodeMaster = acos(edge_vector12.Dot(edge_node1toNodeMaster));
  double Theta_edgeNode2AndNodeMaster = acos(edge_vector21.Dot(edge_node2toNodeMaster));
  double Theta_edgeNode1AndNodeSlave = acos(edge_vector12.Dot(edge_node1toNodeSlave));
  double Theta_edgeNode2AndNodeSlave = acos(edge_vector21.Dot(edge_node2toNodeSlave));

  // Do not flip if these criteria are met

  // Interior angles greater than 120 degrees
  if (Theta_edgeNode1AndNodeMaster>=2.094 ||Theta_edgeNode2AndNodeMaster>=2.094
    || Theta_edgeNode1AndNodeSlave>=2.094 || Theta_edgeNode2AndNodeSlave>=2.094)
    BoundingBoxQuality=false;

  // Interior angles less than 20 degrees
  if (Theta_edgeNode1AndNodeMaster<=0.3490 ||Theta_edgeNode2AndNodeMaster<=0.3490
    || Theta_edgeNode1AndNodeSlave<=0.3490 || Theta_edgeNode2AndNodeSlave<=0.3490)
    BoundingBoxQuality=false;

  // Summation of two greater than 140 degrees
  if (Theta_edgeNode1AndNodeMaster+Theta_edgeNode1AndNodeSlave>=2.4434)
    BoundingBoxQuality=false;

  if (Theta_edgeNode2AndNodeMaster+Theta_edgeNode2AndNodeSlave>=2.4434)
    BoundingBoxQuality=false;


  return BoundingBoxQuality;

}//CheckBoundingBoxQuality


/*--------------------------------------------------------------------------*
 | Check quality of the bounding box consisting of                          |
 | master and slave element for parallel processors (public) mukherjee 09/15|
 *--------------------------------------------------------------------------*/
bool STATMECH::StatMechManagerBilayer::CheckBoundingBoxQualityParallel(const int&    RandEdgeLID,
                                                               Epetra_Vector&     discol)
{
  bool BoundingBoxQuality=true;
  DRT::Element * ele = facediscret_->lColFace(RandEdgeLID);
  DRT::ELEMENTS::DiscSh3Line* edge = dynamic_cast<DRT::ELEMENTS::DiscSh3Line*>(ele);


  DRT::Element* pele = edge->ParentMasterElement();  // Master element
  DRT::Element* nele = edge->ParentSlaveElement();   // Slave element

  if (pele == NULL) dserror("pele is NULL");
  if (nele == NULL) dserror("nele is NULL");

  LINALG::Matrix <1,3> current_pos_master(true);
  LINALG::Matrix <1,3> current_pos_slave(true);
  LINALG::Matrix <1,3> edge_vector12(true);
  LINALG::Matrix <1,3> edge_vector21(true);

  // Find the node which doesn't belong to the master element
  for (int i=0; i<nele->NumNode(); i++)
  {
    if(pele->NodeIds()[i]!=edge->NodeIds()[0] && pele->NodeIds()[i]!=edge->NodeIds()[1])
      current_pos_master= GetSpatialPosition((pele->Nodes())[i],discol);
    if(nele->NodeIds()[i]!=edge->NodeIds()[0] && nele->NodeIds()[i]!=edge->NodeIds()[1])
      current_pos_slave= GetSpatialPosition((nele->Nodes())[i],discol);
  }

  // Find current length of edge
  LINALG::Matrix <1,3> edge_node1toNodeMaster(true);
  LINALG::Matrix <1,3> edge_node2toNodeMaster(true);
  LINALG::Matrix <1,3> edge_node1toNodeSlave(true);
  LINALG::Matrix <1,3> edge_node2toNodeSlave(true);

  LINALG::Matrix <1,3> current_pos_edgeNode1 = GetSpatialPosition((edge->Nodes())[0],discol);
  LINALG::Matrix <1,3> current_pos_edgeNode2 = GetSpatialPosition((edge->Nodes())[1],discol);

  for (int dim=0;dim<3;++dim)
  {
    edge_node1toNodeMaster(dim)=current_pos_master(dim)-current_pos_edgeNode1(dim);
    edge_node2toNodeMaster(dim)=current_pos_master(dim)-current_pos_edgeNode2(dim);
    edge_node1toNodeSlave(dim)=current_pos_slave(dim)-current_pos_edgeNode1(dim);
    edge_node2toNodeSlave(dim)=current_pos_slave(dim)-current_pos_edgeNode2(dim);
    edge_vector12(dim)=current_pos_edgeNode2(dim)-current_pos_edgeNode1(dim);
    edge_vector21(dim)=current_pos_edgeNode1(dim)-current_pos_edgeNode2(dim);
  }

  //Nomralize the vectors
  edge_node1toNodeMaster.Scale(1/edge_node1toNodeMaster.Norm2());
  edge_node2toNodeMaster.Scale(1/edge_node2toNodeMaster.Norm2());
  edge_node1toNodeSlave.Scale(1/edge_node1toNodeSlave.Norm2());
  edge_node2toNodeSlave.Scale(1/edge_node2toNodeSlave.Norm2());
  edge_vector12.Scale(1/edge_vector12.Norm2());
  edge_vector21.Scale(1/edge_vector21.Norm2());


  double Theta_edgeNode1AndNodeMaster = acos(edge_vector12.Dot(edge_node1toNodeMaster));
  double Theta_edgeNode2AndNodeMaster = acos(edge_vector21.Dot(edge_node2toNodeMaster));
  double Theta_edgeNode1AndNodeSlave = acos(edge_vector12.Dot(edge_node1toNodeSlave));
  double Theta_edgeNode2AndNodeSlave = acos(edge_vector21.Dot(edge_node2toNodeSlave));

  // Do not flip if these criteria are met

  // Interior angles greater than 120 degrees
  if (Theta_edgeNode1AndNodeMaster>=2.094 ||Theta_edgeNode2AndNodeMaster>=2.094
    || Theta_edgeNode1AndNodeSlave>=2.094 || Theta_edgeNode2AndNodeSlave>=2.094)
    BoundingBoxQuality=false;

  // Interior angles less than 20 degrees
  if (Theta_edgeNode1AndNodeMaster<=0.3490 ||Theta_edgeNode2AndNodeMaster<=0.3490
    || Theta_edgeNode1AndNodeSlave<=0.3490 || Theta_edgeNode2AndNodeSlave<=0.3490)
    BoundingBoxQuality=false;

  // Summation of two greater than 140 degrees
  if (Theta_edgeNode1AndNodeMaster+Theta_edgeNode1AndNodeSlave>=2.4434)
    BoundingBoxQuality=false;

  if (Theta_edgeNode2AndNodeMaster+Theta_edgeNode2AndNodeSlave>=2.4434)
    BoundingBoxQuality=false;


  return BoundingBoxQuality;

}//CheckBoundingBoxQuality



/*------------------------------------------------------------*
 | Execute bond-flips                 (public) mukherjee 09/15|
 *------------------------------------------------------------*/
void STATMECH::StatMechManagerBilayer::DoBondFlip(const int&    RandEdgeGID)
{
  DRT::FaceElement * edge = facediscret_->gFace( RandEdgeGID );

  DRT::Element* pele = edge->ParentMasterElement();  // Master element
  DRT::Element* nele = edge->ParentSlaveElement();   // Slave element

    if (pele == NULL) dserror("pele is NULL");
    if (nele == NULL) dserror("nele is NULL");


    // map of element ids to be modified with the new node
    // "key" contains the element id, and "value" is dummy here
    // map is used to make sure that one element is stored only once
    std::map<int, int> delEle_master;
    std::map<int, int> delEle_slave;
    //! stores split nodes and their duplicated nodes
    std::map<int,int> oldnew_master;
    std::map<int,int> oldnew_slave;

    int OnlyInSlaveId =0;
    int OnlyInMasterId =0;
    // Find the node which doesn't belong to the master element
    for (int i=0; i<nele->NumNode(); i++)
    {
      if(nele->NodeIds()[i]!=edge->NodeIds()[0] && nele->NodeIds()[i]!=edge->NodeIds()[1])
        OnlyInSlaveId =nele->NodeIds()[i];
      if(pele->NodeIds()[i]!=edge->NodeIds()[0] && pele->NodeIds()[i]!=edge->NodeIds()[1])
        OnlyInMasterId =pele->NodeIds()[i];
    }

    int count=0;
    int MasterToSlaveId=0;
   for (int i=0; i<pele->NumNode(); i++)
   {
     if ((pele->NodeIds()[i]==edge->NodeIds()[0] || pele->NodeIds()[i]==edge->NodeIds()[1]) && count ==0)
     {
       oldnew_master [pele->NodeIds()[i]]=OnlyInSlaveId;
       count++;
     }
     else if ((pele->NodeIds()[i]==edge->NodeIds()[0] || pele->NodeIds()[i]==edge->NodeIds()[1]) && count ==1)
     {
       MasterToSlaveId=pele->NodeIds()[i];
     }
   }

   for (int i=0; i<nele->NumNode(); i++)
   {
     if (nele->NodeIds()[i]==MasterToSlaveId)
     {
       oldnew_slave [nele->NodeIds()[i]]=OnlyInMasterId;
     }
   }

   delEle_master [pele->Id() ] = 0;

   delEle_slave[ nele->Id() ] = 0;


    ModifyElementConnectivity(*discret_,delEle_master, oldnew_master);

    ModifyElementConnectivity(*discret_,delEle_slave, oldnew_slave);

    discret_->FillComplete();
    facediscret_->FillCompleteFaces(true,true,true,true);
    return;
}


/*--------------------------------------------------------------------------*
 | Execute bond-flips for parallel processor        (public) mukherjee 09/15|
 *--------------------------------------------------------------------------*/
void STATMECH::StatMechManagerBilayer::DoBondFlipParallel(DRT::Discretization& mydiscret,
                                                          DRT::DiscretizationFaces& myfacediscret,
                                                          const int&    RandEdgeLID)
{
  DRT::Element * ele = myfacediscret.lColFace(RandEdgeLID);

  DRT::ELEMENTS::DiscSh3Line* edge = dynamic_cast<DRT::ELEMENTS::DiscSh3Line*>(ele);

  DRT::Element* pele = edge->ParentMasterElement();  // Master element
  DRT::Element* nele = edge->ParentSlaveElement();   // Slave element

  if (pele == NULL) dserror("pele is NULL");
  if (nele == NULL) dserror("nele is NULL");


  // map of element ids to be modified with the new node
  // "key" contains the element id, and "value" is dummy here
  // map is used to make sure that one element is stored only once
    std::map<int, int> delEle;
    //! stores split nodes and their duplicated nodes
    std::vector< std::map<int,int> > oldnew;
    std::map<int,int> oldnew_master;
    std::map<int,int> oldnew_slave;
    int OnlyInSlaveId =0;
    int OnlyInMasterId =0;
    // Find the node which doesn't belong to the master element
    for (int i=0; i<nele->NumNode(); i++)
    {
      if(nele->NodeIds()[i]!=edge->NodeIds()[0] && nele->NodeIds()[i]!=edge->NodeIds()[1])
        OnlyInSlaveId =nele->NodeIds()[i];
      if(pele->NodeIds()[i]!=edge->NodeIds()[0] && pele->NodeIds()[i]!=edge->NodeIds()[1])
        OnlyInMasterId =pele->NodeIds()[i];
    }


    int count=0;
    int MasterToSlaveId=0;
   for (int i=0; i<pele->NumNode(); i++)
   {
     if ((pele->NodeIds()[i]==edge->NodeIds()[0] || pele->NodeIds()[i]==edge->NodeIds()[1]) && count ==0)
     {
       oldnew_master [pele->NodeIds()[i]]=OnlyInSlaveId;
       oldnew.push_back(oldnew_master);
       count++;
     }
     else if ((pele->NodeIds()[i]==edge->NodeIds()[0] || pele->NodeIds()[i]==edge->NodeIds()[1]) && count ==1)
     {
       MasterToSlaveId=pele->NodeIds()[i];
     }
   }

   for (int i=0; i<nele->NumNode(); i++)
   {
     if (nele->NodeIds()[i]==MasterToSlaveId)
     {
       oldnew_slave [nele->NodeIds()[i]]=OnlyInMasterId;
       oldnew.push_back(oldnew_slave);
     }
   }

   delEle[pele->Id() ]=0;
   delEle[nele->Id() ]=0;

    ModifyElementConnectivityNew(mydiscret,delEle, oldnew);

    return;
}

/*------------------------------------------------------------------------------------------------------------------------*
 * Modify element connectivity of given discretization
 * discret ---> given discretization
 * delEle  ---> elements for which connectivity is modified (key -- element id, value -- dummy)           mukherjee 09/15
 * oldnew  ---> key -- contains old nodes , val -- contains new nodes
 * This discretization, for delEle elements gets new nodes in place of old nodes
 *-------------------------------------------------------------------------------------------------------------------------*/
void STATMECH::StatMechManagerBilayer::ModifyElementConnectivity(DRT::Discretization& mydiscret,
                                                                 std::map<int,int> & delEle,
                                                                 std::map<int,int> & oldnew )
{
  for( std::map<int,int>::iterator num = delEle.begin(); num != delEle.end(); num++ )
  {
    int eleid = num->first;
    if( mydiscret.HaveGlobalElement( eleid ) )
    {
      bool del = false;
      DRT::Element * ele = mydiscret.gElement( eleid );
      const int * oldnodes = ele->NodeIds();
      std::vector<int> newnodes( ele->NumNode() );

      for( int i = 0; i < ele->NumNode(); i++ )
      {
        std::map<int, int>::iterator delnod = oldnew.find( oldnodes[i] );
        if( delnod != oldnew.end() )
        {
          del = true;
          newnodes[i] = delnod->second;
        }
        else
          newnodes[i] = oldnodes[i];
      }

      if( not del )
        dserror("This element should have atleast one replaceable node\n");

      if( newnodes.size() != static_cast<std::size_t>(ele->NumNode()) )
        dserror("Check the number of new nodes\n");

#if 1
    //-----------------
    // modifying the nodes of an element is easy
    // just modify the node ids of the element
    // when fillcomplete is called the corresponding nodes will be
    //  set through DRT::Element::BuildNodalPointers()
    //-----------------
    ele->SetNodeIds(newnodes.size(), &newnodes[0]);
#endif

    }
  }
}

/*------------------------------------------------------------------------------------------------------------------------*
 * Modify element connectivity of given discretization
 * discret ---> given discretization
 * delEle  ---> elements for which connectivity is modified (key -- element id, value -- dummy)           mukherjee 09/15
 * oldnew  ---> key -- contains old nodes , val -- contains new nodes
 * This discretization, for delEle elements gets new nodes in place of old nodes
 *-------------------------------------------------------------------------------------------------------------------------*/
void STATMECH::StatMechManagerBilayer::ModifyElementConnectivityNew(DRT::Discretization& mydiscret,
                                                                 std::map<int,int> & delEle,
                                                                 std::vector< std::map<int,int> >& oldnew )
{
  int count=0;
  for( std::map<int,int>::iterator num = delEle.begin(); num != delEle.end(); num++ )
  {
    int eleid = num->first;
    if( mydiscret.HaveGlobalElement( eleid ) )
    {
      bool del = false;
      DRT::Element * ele = mydiscret.gElement( eleid );
      const int * oldnodes = ele->NodeIds();
      std::vector<int> newnodes( ele->NumNode() );

      for( int i = 0; i < ele->NumNode(); i++ )
      {
        std::map<int, int>::iterator delnod = oldnew[count].find( oldnodes[i] );
        if( delnod != oldnew[count].end() )
        {
          del = true;
          newnodes[i] = delnod->second;
        }
        else
          newnodes[i] = oldnodes[i];
      }

      if( not del )
        dserror("This element should have atleast one replaceable node\n");

      if( newnodes.size() != static_cast<std::size_t>(ele->NumNode()) )
        dserror("Check the number of new nodes\n");


#if 1
    //-----------------
    // modifying the nodes of an element is easy
    // just modify the node ids of the element
    // when fillcomplete is called the corresponding nodes will be
    //  set through DRT::Element::BuildNodalPointers()
    //-----------------
    ele->SetNodeIds(newnodes.size(), &newnodes[0]);
#endif

    }
    count++;
  }
}

/*----------------------------------------------------------------------*
 | Create Dirichlet DOF maps                  (private)  mukherjee 09/15|
 *----------------------------------------------------------------------*/
void STATMECH::StatMechManagerBilayer::EvaluateNeumannStatMech(Teuchos::ParameterList&              params,
                                                        Teuchos::RCP<Epetra_Vector>          disn,
                                                        Teuchos::RCP<Epetra_Vector>          systemvector,
                                                        Teuchos::RCP<LINALG::SparseOperator> systemmatrix)
{
  // some checks
  if (!discret_->Filled()) dserror("FillComplete() was not called");
  if (!discret_->HaveDofs()) dserror("AssignDegreesOfFreedom() was not called");

  // get load vector
  const Teuchos::ParameterList& sdyn = DRT::Problem::Instance()->StructuralDynamicParams();
  bool loadlin = DRT::INPUT::IntegralValue<INPAR::STR::DynamicType>(sdyn, "LOADLIN");
  INPAR::STATMECH::NBCType nbctype = DRT::INPUT::IntegralValue<INPAR::STATMECH::NBCType>(statmechBilayerparams_,"NBCTYPE");
  switch(nbctype)
  {
    case INPAR::STATMECH::nbctype_std:
    {
      if (!loadlin)
        discret_->EvaluateNeumann(params, systemvector);
      else
      {
        discret_->SetState(0,"displacement new", disn);
        discret_->EvaluateNeumann(params, systemvector, systemmatrix);
      }
    }
    break;
    default: break;
  }
  return;
}

/*----------------------------------------------------------------------*
 |  Evaluate DBCs either in the standard or the "statmech" way          |
 |                                        (public)       mukherjee 09/15|
 *----------------------------------------------------------------------*/
void STATMECH::StatMechManagerBilayer::EvaluateDirichletStatMech(Teuchos::ParameterList&            params,
                                                          Teuchos::RCP<Epetra_Vector>        dis,
                                                          Teuchos::RCP<Epetra_Vector>        vel,
                                                          Teuchos::RCP<LINALG::MapExtractor> dbcmapextractor)
{
  if(dbcmapextractor!=Teuchos::null)
  {
    INPAR::STATMECH::DBCType dbctype = DRT::INPUT::IntegralValue<INPAR::STATMECH::DBCType>(statmechBilayerparams_,"DBCTYPE");

    switch(dbctype)
    {
      // standard DBC application
      case INPAR::STATMECH::dbctype_std:
        discret_->EvaluateDirichlet(params, dis, vel, Teuchos::null, Teuchos::null, dbcmapextractor);
      break;
      // default: everything involving periodic boundary conditions
      default:
      break;
    }
  }
  else
    dserror("Only new DBC application method implemented using the map extractor! Old version using toggle vectors discontinued!");

  return;
}


