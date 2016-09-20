/*----------------------------------------------------------------------*/
/*!
\file beam3contact_manager.cpp

\brief Main class to control beam contact

\level 2

\maintainer Christoph Meier
*/
/*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*/
/* headers */
#include "beam3contact_manager.H"
#include "beam3contact_defines.H"
#include "beam3contact_octtree.H"
#include "beam3contact_utils.H"
#include "../drt_inpar/inpar_beamcontact.H"
#include "../drt_inpar/inpar_beampotential.H"
#include "../drt_inpar/inpar_contact.H"
#include "../drt_inpar/inpar_structure.H"
#include "../linalg/linalg_utils.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_discret.H"
#include "../linalg/linalg_sparsematrix.H"
#include <Teuchos_Time.hpp>

#include "../drt_beam3/beam3.H"
#include "../drt_beam3/beam3r.H"
#include "../drt_beam3/beam3eb.H"
#include "../drt_beam3/beam3k.H"
#include "../drt_rigidsphere/rigidsphere.H"

#include "../drt_inpar/inpar_structure.H"
#include "../drt_contact/contact_element.H"
#include "../drt_contact/contact_node.H"
#include "../drt_io/io.H"

/*----------------------------------------------------------------------*
 |  constructor (public)                                      popp 04/10|
 *----------------------------------------------------------------------*/
CONTACT::Beam3cmanager::Beam3cmanager(DRT::Discretization& discret, double alphaf):
numnodes_(0),
numnodalvalues_(0),
pdiscret_(discret),
pdiscomm_(discret.Comm()),
searchradius_(0.0),
sphericalsearchradius_(0.0),
searchradiuspot_(0.0),
alphaf_(alphaf),
constrnorm_(0.0),
btsolconstrnorm_(0.0),
maxtotalsimgap_(-1000.0),
maxtotalsimgap_cp_(-1000.0),
maxtotalsimgap_gp_(-1000.0),
maxtotalsimgap_ep_(-1000.0),
maxtotalsimrelgap_(-1000.0),
mintotalsimgap_(1000.0),
mintotalsimgap_cp_(1000.0),
mintotalsimgap_gp_(1000.0),
mintotalsimgap_ep_(1000.0),
mintotalsimrelgap_(1000.0),
mintotalsimunconvgap_(1000.0),
totpenaltyenergy_(0.0),
totpenaltywork_(0.0),
maxdeltadisp_(0.0),
totalmaxdeltadisp_(0.0),
firststep_(true),
elementtypeset_(false),
outputcounter_(0),
timen_(0.0),
contactevaluationtime_(0.0),
global_kappa_max_(0.0),
step_(0)
{
  // initialize vectors of contact forces
  fc_ = LINALG::CreateVector(*discret.DofRowMap(), false);
  fcold_ = LINALG::CreateVector(*discret.DofRowMap(), false);
  fc_->PutScalar(0.0);
  fcold_->PutScalar(0.0);

  contactpairmap_.clear();
  oldcontactpairmap_.clear();
  btsolpairmap_.clear();
  oldbtsolpairmap_.clear();
  btsphpairmap_.clear();

  // read parameter lists from DRT::Problem
  sbeamcontact_   = DRT::Problem::Instance()->BeamContactParams();
  sbeampotential_ = DRT::Problem::Instance()->BeamPotentialParams();
  scontact_       = DRT::Problem::Instance()->ContactDynamicParams();
  sstructdynamic_ = DRT::Problem::Instance()->StructuralDynamicParams();

  // indicate if beam-to-solid meshtying, beam-to-solid contact or beam-to-sphere contact is applied
  btsolmt_ = DRT::INPUT::IntegralValue<int>(BeamContactParameters(),"BEAMS_BTSOLMT");
  btsol_   = DRT::INPUT::IntegralValue<int>(BeamContactParameters(),"BEAMS_BTSOL");
  btsph_   = DRT::INPUT::IntegralValue<int>(BeamContactParameters(),"BEAMS_BTSPH");


  InitBeamContactDiscret();

  // check input parameters
  if (sbeamcontact_.get<double>("BEAMS_BTBPENALTYPARAM") < 0.0 ||
      sbeamcontact_.get<double>("BEAMS_BTSPENALTYPARAM") < 0.0 ||
      sbeamcontact_.get<double>("BEAMS_BTSPH_PENALTYPARAM") < 0.0 ||
      sbeamcontact_.get<double>("BEAMS_BTSMTPENALTYPARAM") < 0.0)
  {
    dserror("ERROR: The penalty parameter has to be positive.");
  }

  // initialize beam-to-beam contact element pairs
  pairs_.resize(0);
  oldpairs_.resize(0);

  // initialize beam-to-solid contact element pairs
  btsolpairs_.resize(0);
  oldbtsolpairs_.resize(0);

  // initialize beam-to-solid meshtying element groups
  btsolmtgroups_.resize(0);

  // initialize beam-to-sphere contact element pairs
  btsphpairs_.resize(0);

  // initialize potential-based interaction pairs
  btbpotpairs_.resize(0);
  btsphpotpairs_.resize(0);

  // initialize input parameters
  currentpp_ = sbeamcontact_.get<double>("BEAMS_BTBPENALTYPARAM");
  btspp_     = sbeamcontact_.get<double>("BEAMS_BTSPENALTYPARAM");
  btsmtpp_   = sbeamcontact_.get<double>("BEAMS_BTSMTPENALTYPARAM");

  if (btsph_)
  {
    btsphpp_ = sbeamcontact_.get<double>("BEAMS_BTSPH_PENALTYPARAM");
    if (btsphpp_ == 0.0)
    {
      btsphpp_ = currentpp_;
    }
    else if (btsphpp_ < 0.0)
      dserror("ERROR: The beam-to-sphere penalty parameter has to be positive. Check input file!");
  }

  // initialize Uzawa iteration index
  uzawaiter_ = 0;

  if(!pdiscret_.Comm().MyPID())
  {
    std::cout << "========================= Beam Contact =========================" << std::endl;
    std::cout<<"Elements in discret.   = "<<pdiscret_.NumGlobalElements()<<std::endl;
  }

  //Set maximal and minimal beam/sphere radius occurring in discretization
  SetMinMaxEleRadius();

  //Get search box increment from input file
  searchboxinc_=BEAMCONTACT::DetermineSearchboxInc(sbeamcontact_);

  if(searchboxinc_<0.0)
    dserror("Choose a positive value for the searchbox extrusion factor BEAMS_EXTVAL!");

  // initialize octtree for contact search
  if (DRT::INPUT::IntegralValue<INPAR::BEAMCONTACT::OctreeType>(sbeamcontact_,"BEAMS_OCTREE") != INPAR::BEAMCONTACT::boct_none)
  {
    if (!pdiscret_.Comm().MyPID())
      std::cout << "BTB-CO penalty         = " << currentpp_ << std::endl;

    if (!pdiscret_.Comm().MyPID())
      std::cout << "BTS-CO penalty         = " << btspp_ << std::endl;

    if (!pdiscret_.Comm().MyPID())
      std::cout << "BTS-MT penalty         = " << btsmtpp_ << std::endl;

    tree_ = Teuchos::rcp(new Beam3ContactOctTree(sbeamcontact_,pdiscret_,*btsoldiscret_));
  }
  else
  {
    if(btsol_ || btsolmt_)
      dserror("Beam to solid contact/meshtying are only implemented for the octree contact search!");

    // compute the search radius for searching possible contact pairs
    ComputeSearchRadius();
    tree_ = Teuchos::null;
    if(!pdiscret_.Comm().MyPID())
      std::cout<<"\nBrute Force Search"<<std::endl;
  }

  if(!pdiscret_.Comm().MyPID())
  {
    if (DRT::INPUT::IntegralValue<INPAR::BEAMCONTACT::Strategy>(sbeamcontact_,"BEAMS_STRATEGY") == INPAR::BEAMCONTACT::bstr_penalty )
      std::cout << "Strategy                 Penalty" << std::endl;
    else if (DRT::INPUT::IntegralValue<INPAR::BEAMCONTACT::Strategy>(sbeamcontact_,"BEAMS_STRATEGY") == INPAR::BEAMCONTACT::bstr_uzawa)
    {
      std::cout << "Strategy                 Augmented Lagrange" << std::endl;
      if (DRT::INPUT::IntegralValue<INPAR::BEAMCONTACT::PenaltyLaw>(sbeamcontact_,"BEAMS_PENALTYLAW")!=INPAR::BEAMCONTACT::pl_lp)
        dserror("Augmented Lagrange strategy only implemented for Linear penalty law (LinPen) so far!");
    }
    else if (DRT::INPUT::IntegralValue<INPAR::BEAMCONTACT::Strategy>(sbeamcontact_,"BEAMS_STRATEGY") == INPAR::BEAMCONTACT::bstr_gmshonly )
         std::cout << "Strategy                 Gmsh Only" << std::endl;
    else
      dserror("Unknown strategy for beam contact!");

    switch (DRT::INPUT::IntegralValue<INPAR::BEAMCONTACT::PenaltyLaw>(sbeamcontact_,"BEAMS_PENALTYLAW"))
    {
      case INPAR::BEAMCONTACT::pl_lp:
      {
        std::cout << "Regularization Type      Linear penalty law!" << std::endl;
        break;
      }
      case INPAR::BEAMCONTACT::pl_qp:
      {
        std::cout << "Regularization Type      Quadratic penalty law!" << std::endl;
        break;
      }
      case INPAR::BEAMCONTACT::pl_lnqp:
      {
        std::cout << "Regularization Type      Linear penalty law with quadratic regularization for negative gaps!" << std::endl;
        break;
      }
      case INPAR::BEAMCONTACT::pl_lpqp:
      {
        std::cout << "Regularization Type      Linear penalty law with quadratic regularization for positive gaps!" << std::endl;
        break;
      }
      case INPAR::BEAMCONTACT::pl_lpcp:
      {
        std::cout << "Regularization Type      Linear penalty law with cubic regularization for positive gaps!" << std::endl;
        break;
      }
      case INPAR::BEAMCONTACT::pl_lpdqp:
      {
        std::cout << "Regularization Type      Linear penalty law with double quadratic regularization for positive gaps!" << std::endl;
        break;
      }
      case INPAR::BEAMCONTACT::pl_lpep:
      {
        std::cout << "Regularization Type      Linear penalty law with exponential regularization for positive gaps!" << std::endl;
        break;
      }
    }

    if (DRT::INPUT::IntegralValue<INPAR::BEAMCONTACT::PenaltyLaw>(sbeamcontact_,"BEAMS_PENALTYLAW")!=INPAR::BEAMCONTACT::pl_lp)
    {
      std::cout << "Regularization Params    BEAMS_PENREGPARAM_G0 = " << sbeamcontact_.get<double>("BEAMS_PENREGPARAM_G0",-1.0)\
      << ",  BEAMS_PENREGPARAM_F0 = " << sbeamcontact_.get<double>("BEAMS_PENREGPARAM_F0",-1.0)\
      << ",  BEAMS_PENREGPARAM_C0 = " << sbeamcontact_.get<double>("BEAMS_PENREGPARAM_C0",-1.0)\
      << ",  BEAMS_GAPSHIFTPARAM = " << sbeamcontact_.get<double>("BEAMS_GAPSHIFTPARAM",0.0) << std::endl;
    }

    if (DRT::INPUT::IntegralValue<INPAR::BEAMCONTACT::Damping>(sbeamcontact_,"BEAMS_DAMPING") == INPAR::BEAMCONTACT::bd_no )
      std::cout << "Damping                  No Contact Damping Force Applied!" << std::endl;
    else
    {
      std::cout << "Damping                  BEAMS_DAMPINGPARAM = " << sbeamcontact_.get<double>("BEAMS_DAMPINGPARAM",-1.0)\
      << ",    BEAMS_DAMPREGPARAM1 = " << sbeamcontact_.get<double>("BEAMS_DAMPREGPARAM1",-1.0)\
      << ",   BEAMS_DAMPREGPARAM2 = " << sbeamcontact_.get<double>("BEAMS_DAMPREGPARAM2",-1.0) << std::endl;
    }

    if (sbeamcontact_.get<double>("BEAMS_BASICSTIFFGAP",-1000.0)!=-1000.0)
    {
      std::cout << "Linearization            For gaps < -" << sbeamcontact_.get<double>("BEAMS_BASICSTIFFGAP",-1000.0)\
      << " only the basic part of the contact linearization is applied!"<< std::endl;
    }

    std::cout <<"================================================================\n" << std::endl;
  }

  dis_ = LINALG::CreateVector(*ProblemDiscret().DofRowMap(), true);
  dis_old_ = LINALG::CreateVector(*ProblemDiscret().DofRowMap(), true);



  // read the DLINE conditions specifying charge density of beams
  linechargeconds_.clear();
  ProblemDiscret().GetCondition("BeamPotentialLineCharge", linechargeconds_);

  // initialization stuff for potential-based interaction
  if(linechargeconds_.size()!=0)
  {
    // initialize parameters of applied potential law
    ki_ = Teuchos::rcp(new std::vector<double>);
    mi_ = Teuchos::rcp(new std::vector<double>);
    ki_->clear();
    mi_->clear();
    // read potential law parameters from input and check
    {
      std::istringstream PL(Teuchos::getNumericStringParameter(sbeampotential_,"POT_LAW_EXPONENT"));
      std::string word;
      char* input;
      while (PL >> word)
        mi_->push_back(std::strtod(word.c_str(), &input));
    }
    {
      std::istringstream PL(Teuchos::getNumericStringParameter(sbeampotential_,"POT_LAW_PREFACTOR"));
      std::string word;
      char* input;
      while (PL >> word)
        ki_->push_back(std::strtod(word.c_str(), &input));
    }
    if (!ki_->empty())
    {
      if (ki_->size()!=mi_->size())
        dserror("number of potential law prefactors does not match number of potential law exponents. Check your input file!");

      for (unsigned int i=0; i<mi_->size(); ++i)
        if (mi_->at(i) <= 0)
          dserror("only positive values are allowed for potential law exponent. Check your input file");
    }

    if(!pdiscret_.Comm().MyPID())
    {
      std::cout << "=============== Beam Potential-Based Interaction ===============" << std::endl;

      switch (DRT::INPUT::IntegralValue<INPAR::BEAMPOTENTIAL::BeamPotentialType>(sbeampotential_,"BEAMPOTENTIAL_TYPE"))
      {
        case INPAR::BEAMPOTENTIAL::beampot_surf:
        {
          std::cout << "Potential Type:      Surface" << std::endl;
          break;
        }
        case INPAR::BEAMPOTENTIAL::beampot_vol:
        {
          std::cout << "Potential Type:      Volume" << std::endl;
          break;
        }
      }

      std::cout << "Potential Law:       Phi(r) = ";
      for (unsigned int i=0; i<ki_->size(); ++i)
      {
        if (i>0) std::cout << " + ";
        std::cout << "(" << ki_->at(i) << ") * r^(-" << mi_->at(i) << ")";
      }
      std::cout << std::endl;
    }

    // initialize octtree for search of potential-based interaction pairs
    if (DRT::INPUT::IntegralValue<INPAR::BEAMCONTACT::OctreeType>(sbeampotential_,"BEAMPOT_OCTREE") != INPAR::BEAMCONTACT::boct_none)
    {
      pottree_ = Teuchos::rcp(new Beam3ContactOctTree(sbeampotential_,pdiscret_,*btsoldiscret_));
    }
    else
    {
      // read cutoff radius for search of potential-based interaction pairs
      searchradiuspot_ = sbeampotential_.get<double>("CUTOFFRADIUS",-1.0);
      if (searchradiuspot_<=0)
        dserror("no/invalid value for cutoff radius of potential-based interaction pairs specified. Check your input file!");

      // Compute the search radius for searching possible contact pairs
      ComputeSearchRadius();
      pottree_ = Teuchos::null;
      if(!pdiscret_.Comm().MyPID())
      {
        std::cout << "\nSearch Strategy:     Brute Force Search" << std::endl;
        std::cout << "Search Radius:       " << searchradiuspot_ << std::endl;
      }
    }

    if(!pdiscret_.Comm().MyPID())
    {
      std::cout << "================================================================\n" << std::endl;
    }

    //Parameters to indicate, if beam-to-solid or beam-to-sphere potential-based interaction is applied
    potbtsol_ = DRT::INPUT::IntegralValue<int>(sbeampotential_,"BEAMPOT_BTSOL");
    potbtsph_ = DRT::INPUT::IntegralValue<int>(sbeampotential_,"BEAMPOT_BTSPH");

    // build a map telling us which nodes lie on which DLINE
    dlinenodemap_.clear();
    for (unsigned int i=0; i<linechargeconds_.size(); ++i)
    {
      if (linechargeconds_[i]->Type() != DRT::Condition::BeamPotential_LineChargeDensity)
        dserror("The specified DLINE conditions are not of correct type BeamPotential_LineChargeDensity");

      const std::vector<int>*    node_ids = linechargeconds_[i]->Nodes();
      for (unsigned int j=0; j<node_ids->size(); ++j)
      {
        dlinenodemap_[(*node_ids)[j]] = i;
      }
    }
  } //beam potential line charge condition applied

  return;
}

/*----------------------------------------------------------------------*
 |  print beam3 contact manager (public)                      popp 04/10|
 *----------------------------------------------------------------------*/
void CONTACT::Beam3cmanager::Print(std::ostream& os) const
{
  if (Comm().MyPID()==0)
    os << "Beam3 Contact Discretization:" << std::endl;

  ProblemDiscret().Print(os);

  return;
}

/*----------------------------------------------------------------------*
 |  evaluate contact (public)                                 popp 04/10|
 *----------------------------------------------------------------------*/
void CONTACT::Beam3cmanager::Evaluate(LINALG::SparseMatrix& stiffmatrix,
                                      Epetra_Vector& fres,
                                      const Epetra_Vector& disrow,
                                      Teuchos::ParameterList timeintparams,
                                      bool newsti,
                                      double time)
{
  // get out of here if only interested in gmsh output
  if (DRT::INPUT::IntegralValue<INPAR::BEAMCONTACT::Strategy>(sbeamcontact_,"BEAMS_STRATEGY") == INPAR::BEAMCONTACT::bstr_gmshonly )
    return;

  //set time
  timen_=time;

  //Set class variable
  dis_->Update(1.0,disrow,0.0);

  // map linking node numbers and current node positions
  std::map<int,LINALG::Matrix<3,1> > currentpositions;
  currentpositions.clear();
  //extract fully overlapping displacement vector on contact discretization from
  //displacement vector in row map format on problem discretization
  Epetra_Vector disccol(*BTSolDiscret().DofColMap(),true);
  ShiftDisMap(disrow, disccol);
  // update currentpositions
  SetCurrentPositions(currentpositions,disccol);

  double t_start = 0.0;
  double t_end = 0.0;

  //**********************************************************************
  // SEARCH
  //**********************************************************************
  std::vector<std::vector<DRT::Element* > > elementpairs, elementpairspot;
  elementpairs.clear();
  elementpairspot.clear();
  //**********************************************************************
  // Contact: Octree search
  //**********************************************************************
  if (tree_ != Teuchos::null)
  {
     t_start = Teuchos::Time::wallTime();
     elementpairs = tree_->OctTreeSearch(currentpositions);

     t_end = Teuchos::Time::wallTime() - t_start;
     Teuchos::ParameterList ioparams = DRT::Problem::Instance()->IOParams();
     if(!pdiscret_.Comm().MyPID() && ioparams.get<int>("STDOUTEVRY",0))
       std::cout << "      OctTree Search (Contact): " << t_end << " seconds, found pairs: "<<elementpairs.size()<< std::endl;
  }
  //**********************************************************************
  // Contact: brute-force search
  //**********************************************************************
  else
  {
    t_start = Teuchos::Time::wallTime();

    elementpairs = BruteForceSearch(currentpositions,searchradius_,sphericalsearchradius_);
    t_end = Teuchos::Time::wallTime() - t_start;
    Teuchos::ParameterList ioparams = DRT::Problem::Instance()->IOParams();
    if(!pdiscret_.Comm().MyPID() && ioparams.get<int>("STDOUTEVRY",0))
      std::cout << "      Brute Force Search (Contact): " << t_end << " seconds" << std::endl;
  }

  t_start = Teuchos::Time::wallTime();

  // process the found element pairs and fill the BTB, BTSOL, BTSPH interaction pair vectors
  FillContactPairsVectors(elementpairs);

  if(linechargeconds_.size()!=0)
  {
    //**********************************************************************
    // Potential-based interaction: Octree search
    //**********************************************************************
    if(pottree_ != Teuchos::null)
    {
       double t_start = Teuchos::Time::wallTime();
       elementpairspot = pottree_->OctTreeSearch(currentpositions);

       double t_end = Teuchos::Time::wallTime() - t_start;
       Teuchos::ParameterList ioparams = DRT::Problem::Instance()->IOParams();
       if(!pdiscret_.Comm().MyPID() && ioparams.get<int>("STDOUTEVRY",0))
         std::cout << "           OctTree Search (Potential): " << t_end << " seconds, found pairs: "<<elementpairspot.size()<< std::endl;
    }
    //**********************************************************************
    // Potential-based interaction: brute-force search
    //**********************************************************************
    else
    {
      double t_start = Teuchos::Time::wallTime();

      elementpairspot = BruteForceSearch(currentpositions,searchradiuspot_,searchradiuspot_); // TODO do we need a sphericalsearchradius here as well?
      double t_end = Teuchos::Time::wallTime() - t_start;
      Teuchos::ParameterList ioparams = DRT::Problem::Instance()->IOParams();
      if(!pdiscret_.Comm().MyPID() && ioparams.get<int>("STDOUTEVRY",0))
        std::cout << "           Brute Force Search (Potential): " << t_end << " seconds" << std::endl;
    }

    FillPotentialPairsVectors(elementpairspot);
  }

  // update element state of all pairs with current positions (already calculated in SetCurrentPositions) and current tangents (will be calculated in SetState)
  SetState(currentpositions,disccol);


  // At this point we have all possible contact pairs_ with updated positions
  //**********************************************************************
  // evaluation of contact pairs
  //**********************************************************************
  // every proc that owns one of the nodes of a 'beam3contact' object has
  // to evaluate this object. Fc and Stiffc will be evaluated. Assembly
  // of the additional stiffness will be done by the objects themselves,
  // the assembly of Fc has to be done by the 'beam3cmanager', because the
  // additional force has to be known for the current and the last time
  // step due to generalized alpha time integration. The current contact
  // forces will be stored in 'fc_', the previous ones in 'fcold_'. An
  // update method at the end of each time step manages the data transfer
  // from 'fc_' to 'fcold_'. This update method is called by the time
  // integration class.
  //**********************************************************************
  // initialize global contact force vectors
  fc_->PutScalar(0.0);

  // initialize contact stiffness and uncomplete global stiffness
  stiffc_ = Teuchos::rcp(new LINALG::SparseMatrix(stiffmatrix.RangeMap(),100));
  stiffmatrix.UnComplete();


  t_end = Teuchos::Time::wallTime() - t_start;
  if(!pdiscret_.Comm().MyPID())
    std::cout << "      Pair management: " << t_end << " seconds. "<< std::endl;
  t_start = Teuchos::Time::wallTime();

  // evaluate all element pairs (BTB, BTSOL, BTSPH; Contact and Potential)
  EvaluateAllPairs(timeintparams);

  t_end = Teuchos::Time::wallTime() - t_start;
  if(!pdiscret_.Comm().MyPID())
    std::cout << "      Evaluate Contact Pairs: " << t_end << " seconds. "<< std::endl;
  double sumproc_evaluationtime=0.0;
  Comm().SumAll(&t_end,&sumproc_evaluationtime,1);
  contactevaluationtime_+=sumproc_evaluationtime;
  t_start = Teuchos::Time::wallTime();

  if (DRT::INPUT::IntegralValue<INPAR::STR::MassLin>(sstructdynamic_,"MASSLIN") != INPAR::STR::ml_rotations)
  {
    // assemble contact forces into global fres vector
    fres.Update(1.0-alphaf_,*fc_,1.0);
    fres.Update(alphaf_,*fcold_,1.0);
    // determine contact stiffness matrix scaling factor (new STI)
    // (this is due to the fact that in the new STI, we hand in the
    // already appropriately scaled effective stiffness matrix. Thus,
    // the additional contact stiffness terms must be equally scaled
    // here, as well. In the old STI, the complete scaling operation
    // is done after contact evaluation within the time integrator,
    // therefore no special scaling needs to be applied here.)
    double scalemat = 1.0;
    if (newsti) scalemat = 1.0 - alphaf_;
    // assemble contact stiffness into global stiffness matrix
    stiffc_->Complete();
    stiffmatrix.Add(*stiffc_,false,scalemat,1.0);
    stiffmatrix.Complete();
  }
  else
  {
    // assemble contact forces into global fres vector
    fres.Update(1.0,*fc_,1.0);

    // assemble contact stiffness into global stiffness matrix
    stiffc_->Complete();
    stiffmatrix.Add(*stiffc_,false,1.0,1.0);
    stiffmatrix.Complete();
  }

  //With this line output can be printed every Newton step
  #ifdef OUTPUTEVERYNEWTONSTEP
    ConsoleOutput();
  #endif

  t_end = Teuchos::Time::wallTime() - t_start;
  if(!pdiscret_.Comm().MyPID())
    std::cout << "      Post-manage Pairs: " << t_end << " seconds. "<< std::endl;

  return;
}

/*----------------------------------------------------------------------*
 |  Shift map of displacement vector                         meier 05/14|
 *----------------------------------------------------------------------*/
void CONTACT::Beam3cmanager::ShiftDisMap(const Epetra_Vector& disrow, Epetra_Vector& disccol)
{

  // export displacements into fully overlapping column map format
  Epetra_Vector discrow(*BTSolDiscret().DofRowMap(),true);
  int numbtsdofs = (*BTSolDiscret().DofRowMap()).NumMyElements();

  for (int i=0;i<numbtsdofs;i++)
  {
    int btsolcontact_gid = (*BTSolDiscret().DofRowMap()).GID(i);
    int problem_gid = dofoffsetmap_[btsolcontact_gid];
    double disp = disrow[(*ProblemDiscret().DofRowMap()).LID(problem_gid)];
    discrow.ReplaceGlobalValue(btsolcontact_gid,0,disp);
  }
  LINALG::Export(discrow,disccol);

  return;
}

/*----------------------------------------------------------------------*
 | setup of contact discretization btsoldiscret_             grill 05/16|
 *----------------------------------------------------------------------*/
void CONTACT::Beam3cmanager::InitBeamContactDiscret()
{
  // create new (basically copied) discretization for contact
  // (to ease our search algorithms we afford the luxury of
  // ghosting all nodes and all elements on all procs, i.e.
  // we export the discretization to full overlap. However,
  // we do not want to do this with the actual discretization
  // and thus create a stripped copy here that only contains
  // nodes and elements).
  // Then, within all beam contact specific routines we will
  // NEVER use the underlying problem discretization but always
  // the copied beam contact discretization.

  Teuchos::RCP<Epetra_Comm> comm = Teuchos::rcp(pdiscret_.Comm().Clone());
  btsoldiscret_ = Teuchos::rcp(new DRT::Discretization((std::string)"beam to solid contact",comm));
  dofoffsetmap_.clear();
  std::map<int,std::vector<int> > nodedofs;
  nodedofs.clear();

  // loop over all column nodes of underlying problem discret and add
  for (int i=0;i<(ProblemDiscret().NodeColMap())->NumMyElements();++i)
  {
    DRT::Node* node = ProblemDiscret().lColNode(i);
    if (!node) dserror("Cannot find node with lid %",i);
    Teuchos::RCP<DRT::Node> newnode = Teuchos::rcp(node->Clone());
    if (BEAMCONTACT::BeamNode(*newnode))
    {
      BTSolDiscret().AddNode(newnode);
      nodedofs[node->Id()]=ProblemDiscret().Dof(0,node);
    }
    else if (BEAMCONTACT::RigidsphereNode(*newnode) and btsph_)
    {
      BTSolDiscret().AddNode(newnode);
      nodedofs[node->Id()]=ProblemDiscret().Dof(0,node);
    }
    else
    {
      if (btsol_==false && btsolmt_==false)
        dserror("Only beam elements are allowed as long as the flags btsol_ and btsolmt_ are set to false!");
    }
  }

  int maxproblemid=ProblemDiscret().ElementRowMap()->MaxAllGID();
  // loop over all column elements of underlying problem discret and add
  for (int i=0;i<(ProblemDiscret().ElementColMap())->NumMyElements();++i)
  {
    DRT::Element* ele = ProblemDiscret().lColElement(i);
    if (!ele) dserror("Cannot find element with lid %",i);
    Teuchos::RCP<DRT::Element> newele = Teuchos::rcp(ele->Clone());
    if (BEAMCONTACT::BeamElement(*newele) or BEAMCONTACT::RigidsphereElement(*newele))
    {
      BTSolDiscret().AddElement(newele);
    }
  }

  // begin: determine surface elements and their nodes

  // vector that contains solid-to-solid and beam-to-solid contact pairs
  std::vector<DRT::Condition*> beamandsolidcontactconditions(0);
  ProblemDiscret().GetCondition("Contact", beamandsolidcontactconditions);

  // vector that solely contains beam-to-solid contact pairs
  std::vector<DRT::Condition*> btscontactconditions(0);

  // vector that solely contains beam-to-solid meshtying pairs
  std::vector<DRT::Condition*> btsmeshtyingconditions(0);

  // sort out solid-to-solid contact pairs, since these are treated in the drt_contact framework
  for (int i = 0; i < (int) beamandsolidcontactconditions.size(); ++i)
  {
    if(*(beamandsolidcontactconditions[i]->Get<std::string>("Application"))=="Beamtosolidcontact")
      btscontactconditions.push_back(beamandsolidcontactconditions[i]);
    if(*(beamandsolidcontactconditions[i]->Get<std::string>("Application"))=="Beamtosolidmeshtying")
      btsmeshtyingconditions.push_back(beamandsolidcontactconditions[i]);
  }

  //******************************* BEAM-TO-SOLID CONTACT *******************************
  solcontacteles_.resize(0);
  solcontactnodes_.resize(0);
  int ggsize = 0;

   //-------------------------------------------------- process surface nodes
   for (int j=0;j<(int)btscontactconditions.size();++j)
   {
     // get all nodes and add them
     const std::vector<int>* nodeids = btscontactconditions[j]->Nodes();
     if (!nodeids) dserror("Condition does not have Node Ids");
     for (int k=0; k<(int)(*nodeids).size(); ++k)
     {
       int gid = (*nodeids)[k];
       // do only nodes that I have in my discretization
       if (!ProblemDiscret().NodeColMap()->MyGID(gid)) continue;
       DRT::Node* node = ProblemDiscret().gNode(gid);

       if (!node) dserror("Cannot find node with gid %",gid);

       Teuchos::RCP<CONTACT::CoNode> cnode = Teuchos::rcp(new CONTACT::CoNode(node->Id(),node->X(),
                                                        node->Owner(),
                                                        ProblemDiscret().NumDof(0,node),
                                                        ProblemDiscret().Dof(0,node),
                                                        false,    //all solid elements are master elements
                                                        false));  //no "initially active" decision necessary for beam to solid contact

       // note that we do not have to worry about double entries
       // as the AddNode function can deal with this case!
       // the only problem would have occured for the initial active nodes,
       // as their status could have been overwritten, but is prevented
       // by the "foundinitialactive" block above!
       solcontactnodes_.push_back(cnode);
       BTSolDiscret().AddNode(cnode);
       nodedofs[node->Id()]=ProblemDiscret().Dof(0,node);
     }
   }

  //----------------------------------------------- process surface elements
  for (int j=0;j<(int)btscontactconditions.size();++j)
  {
    // get elements from condition j of current group
    std::map<int,Teuchos::RCP<DRT::Element> >& currele = btscontactconditions[j]->Geometry();

    // elements in a boundary condition have a unique id
    // but ids are not unique among 2 distinct conditions
    // due to the way elements in conditions are build.
    // We therefore have to give the second, third,... set of elements
    // different ids. ids do not have to be continous, we just add a large
    // enough number ggsize to all elements of cond2, cond3,... so they are
    // different from those in cond1!!!
    // note that elements in ele1/ele2 already are in column (overlapping) map
    int lsize = (int)currele.size();
    int gsize = 0;
    Comm().SumAll(&lsize,&gsize,1);

    std::map<int,Teuchos::RCP<DRT::Element> >::iterator fool;
    for (fool=currele.begin(); fool != currele.end(); ++fool)
    {
      //The IDs of the surface elements of each conditions begin with zero. Therefore we have to add ggsize in order to
      //get unique element IDs in the end. Furthermore, only the solid elements are added to the contact discretization btsoldiscret_
      //via the btscontactconditions, whereas all beam elements with their original ID are simply cloned from the problem discretization
      //into the contact discretization. In order to avoid solid element IDs being identical to these beam element IDs within the contact
      //discretization we have to add the additional offset maxproblemid, which is identical to the maximal element ID in the problem
      //discretization.
      Teuchos::RCP<DRT::Element> ele = fool->second;
      Teuchos::RCP<CONTACT::CoElement> cele = Teuchos::rcp(new CONTACT::CoElement(ele->Id()+ggsize+maxproblemid+1,
                                                              ele->Owner(),
                                                              ele->Shape(),
                                                              ele->NumNode(),
                                                              ele->NodeIds(),
                                                              false,    //all solid elements are master elements
                                                              false));  //no nurbs allowed up to now

      solcontacteles_.push_back(cele);
      BTSolDiscret().AddElement(cele);
    } // for (fool=ele1.start(); fool != ele1.end(); ++fool)
    ggsize += gsize; // update global element counter
  }
  //end: determine surface elements and their nodes

  //****************************** BEAM-TO-SOLID MESHTYING ******************************
  solmeshtyingeles_.resize(0);
  solmeshtyingnodes_.resize(0);

   //-------------------------------------------------- process surface nodes
   for (int j=0;j<(int)btsmeshtyingconditions.size();++j)
   {
     // get all nodes and add them
     const std::vector<int>* nodeids = btsmeshtyingconditions[j]->Nodes();
     if (!nodeids) dserror("Condition does not have Node Ids");
     for (int k=0; k<(int)(*nodeids).size(); ++k)
     {
       int gid = (*nodeids)[k];
       // do only nodes that I have in my discretization
       if (!ProblemDiscret().NodeColMap()->MyGID(gid)) continue;
       DRT::Node* node = ProblemDiscret().gNode(gid);

       if (!node) dserror("Cannot find node with gid %",gid);

       Teuchos::RCP<MORTAR::MortarNode> mtnode = Teuchos::rcp(new MORTAR::MortarNode(node->Id(),node->X(),
                                                        node->Owner(),
                                                        ProblemDiscret().NumDof(0,node),
                                                        ProblemDiscret().Dof(0,node),
                                                        false));    //all solid elements are master elements

       // note that we do not have to worry about double entries
       // as the AddNode function can deal with this case!
       // the only problem would have occured for the initial active nodes,
       // as their status could have been overwritten, but is prevented
       // by the "foundinitialactive" block above!
       solmeshtyingnodes_.push_back(mtnode);
       BTSolDiscret().AddNode(mtnode);
       nodedofs[node->Id()]=ProblemDiscret().Dof(0,node);
     }
   }

  //----------------------------------------------- process surface elements
  for (int j=0;j<(int)btsmeshtyingconditions.size();++j)
  {
    // get elements from condition j of current group
    std::map<int,Teuchos::RCP<DRT::Element> >& currele = btsmeshtyingconditions[j]->Geometry();

    // elements in a boundary condition have a unique id
    // but ids are not unique among 2 distinct conditions
    // due to the way elements in conditions are build.
    // We therefore have to give the second, third,... set of elements
    // different ids. ids do not have to be continous, we just add a large
    // enough number ggsize to all elements of cond2, cond3,... so they are
    // different from those in cond1!!!
    // note that elements in ele1/ele2 already are in column (overlapping) map
    int lsize = (int)currele.size();
    int gsize = 0;
    Comm().SumAll(&lsize,&gsize,1);

    std::map<int,Teuchos::RCP<DRT::Element> >::iterator fool;
    for (fool=currele.begin(); fool != currele.end(); ++fool)
    {
      //The IDs of the surface elements of each conditions begin with zero. Therefore we have to add ggsize in order to
      //get unique element IDs in the end. Furthermore, only the solid elements are added to the contact discretization btsoldiscret_
      //via the btscontactconditions, whereas all beam elements with their original ID are simply cloned from the problem discretization
      //into the contact discretization. In order to avoid solid element IDs being identical to these beam element IDs within the contact
      //discretization we have to add the additional offset maxproblemid, which is identical to the maximal element ID in the problem
      //discretization.
      Teuchos::RCP<DRT::Element> ele = fool->second;
      Teuchos::RCP<MORTAR::MortarElement> mtele = Teuchos::rcp(new MORTAR::MortarElement(ele->Id()+ggsize+maxproblemid+1,
                                                              ele->Owner(),
                                                              ele->Shape(),
                                                              ele->NumNode(),
                                                              ele->NodeIds(),
                                                              false,    //all solid elements are master elements
                                                              false));  //no nurbs allowed up to now

      solmeshtyingeles_.push_back(mtele);
      BTSolDiscret().AddElement(mtele);
    } // for (fool=ele1.start(); fool != ele1.end(); ++fool)
    ggsize += gsize; // update global element counter
  }
  //end: determine surface elements and their nodes

  // build maps but do not assign dofs yet, we'll do this below
  // after shuffling around of nodes and elements (saves time)
  BTSolDiscret().FillComplete(false,false,false);

  // store the node and element row and column maps into this manager
  noderowmap_ = Teuchos::rcp(new Epetra_Map(*(BTSolDiscret().NodeRowMap())));
  elerowmap_  = Teuchos::rcp(new Epetra_Map(*(BTSolDiscret().ElementRowMap())));
  nodecolmap_ = Teuchos::rcp(new Epetra_Map(*(BTSolDiscret().NodeColMap())));
  elecolmap_  = Teuchos::rcp(new Epetra_Map(*(BTSolDiscret().ElementColMap())));

  // build fully overlapping node and element maps
  // fill my own row node ids into vector (e)sdata
  std::vector<int> sdata(noderowmap_->NumMyElements());
  std::vector<int> esdata(elerowmap_->NumMyElements());
  for (int i=0; i<noderowmap_->NumMyElements(); ++i)
    sdata[i] = noderowmap_->GID(i);
  for (int i=0; i<elerowmap_->NumMyElements(); ++i)
    esdata[i] = elerowmap_->GID(i);

  // if current proc is participating it writes row IDs into (e)stproc
  std::vector<int> stproc(0);
  std::vector<int> estproc(0);
  if (noderowmap_->NumMyElements())
    stproc.push_back(BTSolDiscret().Comm().MyPID());
  if (elerowmap_->NumMyElements())
    estproc.push_back(BTSolDiscret().Comm().MyPID());

  // information how many processors participate in total
  std::vector<int> allproc(BTSolDiscret().Comm().NumProc());
  for (int i=0;i<BTSolDiscret().Comm().NumProc();++i) allproc[i] = i;

  // declaring new variables into which the info of (e)stproc on all processors is gathered
  std::vector<int> rtproc(0);
  std::vector<int> ertproc(0);

  // gathers information of (e)stproc and writes it into (e)rtproc; in the end (e)rtproc
  // is a vector which contains the numbers of all processors which own nodes/elements.
  LINALG::Gather<int>(stproc,rtproc,BTSolDiscret().Comm().NumProc(),&allproc[0],BTSolDiscret().Comm());
  LINALG::Gather<int>(estproc,ertproc,BTSolDiscret().Comm().NumProc(),&allproc[0],BTSolDiscret().Comm());

  // in analogy to (e)stproc and (e)rtproc the variables (e)rdata gather all the row ID
  // numbers which are  stored on different processors in their own variables (e)sdata; thus,
  // each processor gets the information about all the row ID numbers existing in the problem
  std::vector<int> rdata;
  std::vector<int> erdata;

  // gather all gids of nodes redundantly from (e)sdata into (e)rdata
  LINALG::Gather<int>(sdata,rdata,(int)rtproc.size(),&rtproc[0],BTSolDiscret().Comm());
  LINALG::Gather<int>(esdata,erdata,(int)ertproc.size(),&ertproc[0],BTSolDiscret().Comm());

  // build completely overlapping node map (on participating processors)
  Teuchos::RCP<Epetra_Map> newnodecolmap = Teuchos::rcp(new Epetra_Map(-1,(int)rdata.size(),&rdata[0],0,BTSolDiscret().Comm()));
  sdata.clear();
  stproc.clear();
  rdata.clear();
  allproc.clear();

  // build completely overlapping element map (on participating processors)
  Teuchos::RCP<Epetra_Map> newelecolmap = Teuchos::rcp(new Epetra_Map(-1,(int)erdata.size(),&erdata[0],0,BTSolDiscret().Comm()));
  esdata.clear();
  estproc.clear();
  erdata.clear();

  // store the fully overlapping node and element maps
  nodefullmap_ = Teuchos::rcp(new Epetra_Map(*newnodecolmap));
  elefullmap_  = Teuchos::rcp(new Epetra_Map(*newelecolmap));

  // pass new fully overlapping node and element maps to beam contact discretization
  BTSolDiscret().ExportColumnNodes(*newnodecolmap);
  BTSolDiscret().ExportColumnElements(*newelecolmap);

  // complete beam contact discretization based on the new column maps
  // (this also assign new degrees of freedom what we actually do not
  // want, thus we have to introduce a dof mapping next)
  BTSolDiscret().FillComplete(true,false,false);

  // communicate the map nodedofs to all proccs
  DRT::Exporter ex(*(ProblemDiscret().NodeColMap()),*(BTSolDiscret().NodeColMap()),Comm());
  ex.Export(nodedofs);

  //Determine offset between the IDs of problem discretization and BTSol discretization
  for (int i=0;i<(BTSolDiscret().NodeColMap())->NumMyElements();++i)
  {
    DRT::Node* node = BTSolDiscret().lColNode(i);
    int nodeid = node->Id();
    std::vector<int> btsolnodedofids = BTSolDiscret().Dof(0,node);
    std::vector<int> originalnodedofids = nodedofs[nodeid];

    if (btsolnodedofids.size() != originalnodedofids.size())
      dserror("Number of nodal DoFs does not match! "
          "node (GID %d) originally had %d DoFs, now in BTSOLdiscret %d DoFs!",
          nodeid, originalnodedofids.size(),btsolnodedofids.size());

    for (int j=0;j<(int)btsolnodedofids.size();j++)
    {
      dofoffsetmap_[btsolnodedofids[j]]=originalnodedofids[j];
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 |  Set current displacement state                            popp 04/10|
 *----------------------------------------------------------------------*/
void CONTACT::Beam3cmanager::SetCurrentPositions(std::map<int,LINALG::Matrix<3,1> >& currentpositions,
                                                 const Epetra_Vector& disccol)
{
  //**********************************************************************
  // get positions of all nodes (fully overlapping map) and store
  // the current results into currentpositions
  //**********************************************************************

  // loop over all beam contact nodes
  for (int i=0;i<FullNodes()->NumMyElements();++i)
  {
    // get node pointer
    DRT::Node* node = BTSolDiscret().lColNode(i);

     // TODO maybe this can be done in a more elegant way in the future
    /* check whether node is a beam node which is NOT used for centerline interpolation
     * if so, we simply skip it because it does not have position (and tangent) DoFs */
    if (BEAMCONTACT::BeamNode(*node) and !BEAMCONTACT::BeamCenterlineNode(*node))
      continue;

    // get GIDs of this node's degrees of freedom
    std::vector<int> dofnode = BTSolDiscret().Dof(node);

    // nodal positions
    LINALG::Matrix<3,1> currpos;
    currpos(0) = node->X()[0] + disccol[BTSolDiscret().DofColMap()->LID(dofnode[0])];
    currpos(1) = node->X()[1] + disccol[BTSolDiscret().DofColMap()->LID(dofnode[1])];
    currpos(2) = node->X()[2] + disccol[BTSolDiscret().DofColMap()->LID(dofnode[2])];

    // store into currentpositions
    currentpositions[node->Id()] = currpos;
  }

  return;
}

/*----------------------------------------------------------------------*
 |  Set current displacement state                            popp 04/10|
 *----------------------------------------------------------------------*/
void CONTACT::Beam3cmanager::SetState(std::map<int,LINALG::Matrix<3,1> >& currentpositions,
                                      const Epetra_Vector& disccol)
{

  // map to store the nodal tangent vectors (necessary for Kirchhoff type beams) and adress it with the node ID
  std::map<int,LINALG::Matrix<3,1> > currenttangents;
  currenttangents.clear();

  // Update of nodal tangents for Kirchhoff elements; nodal positions have already been set in SetCurrentPositions
  // loop over all beam contact nodes
  for (int i=0;i<FullNodes()->NumMyElements();++i)
  {
    // get node pointer
    DRT::Node* node = BTSolDiscret().lColNode(i);

    // TODO maybe this can be done in a more elegant way in the future
   /* check whether node is a beam node which is NOT used for centerline interpolation
    * if so, we simply skip it because it does not have position (and tangent) DoFs */
   if (BEAMCONTACT::BeamNode(*node) and !BEAMCONTACT::BeamCenterlineNode(*node))
     continue;


    //get nodal tangents for Kirchhoff elements
    if (numnodalvalues_==2 and BEAMCONTACT::BeamNode(*node))
    {
      // get GIDs of this node's degrees of freedom
      std::vector<int> dofnode = BTSolDiscret().Dof(node);

      LINALG::Matrix<3,1> currtan(true);
      for (int i=0; i<numnodes_;i++)    // TODO for now, use number of centerline nodes numnodes_ (=2) (no matter how many nodes the function call node->Elements()[0]->NumNode() would tell you)
      {
        if (node->Elements()[0]->Nodes()[i]->Id()==node->Id() and  node->Elements()[0]->ElementType() == DRT::ELEMENTS::Beam3ebType::Instance() )
        {
          const DRT::ELEMENTS::Beam3eb* ele = dynamic_cast<const DRT::ELEMENTS::Beam3eb*>(node->Elements()[0]);
          currtan(0)=((ele->Tref())[i])(0) + disccol[BTSolDiscret().DofColMap()->LID(dofnode[3])];
          currtan(1)=((ele->Tref())[i])(1) + disccol[BTSolDiscret().DofColMap()->LID(dofnode[4])];
          currtan(2)=((ele->Tref())[i])(2) + disccol[BTSolDiscret().DofColMap()->LID(dofnode[5])];
        }
        else if (node->Elements()[0]->Nodes()[i]->Id()==node->Id() and node->Elements()[0]->ElementType() == DRT::ELEMENTS::Beam3kType::Instance() )
        {
          const DRT::ELEMENTS::Beam3k* ele = dynamic_cast<const DRT::ELEMENTS::Beam3k*>(node->Elements()[0]);
          currtan(0)=((ele->Tref())[i])(0) + disccol[BTSolDiscret().DofColMap()->LID(dofnode[3])];
          currtan(1)=((ele->Tref())[i])(1) + disccol[BTSolDiscret().DofColMap()->LID(dofnode[4])];
          currtan(2)=((ele->Tref())[i])(2) + disccol[BTSolDiscret().DofColMap()->LID(dofnode[5])];
        }
        else if (node->Elements()[0]->Nodes()[i]->Id()==node->Id() and node->Elements()[0]->ElementType() == DRT::ELEMENTS::Beam3rType::Instance() )
        {
          const DRT::ELEMENTS::Beam3r* ele = dynamic_cast<const DRT::ELEMENTS::Beam3r*>(node->Elements()[0]);
          currtan(0)=((ele->Tref())[i])(0) + disccol[BTSolDiscret().DofColMap()->LID(dofnode[6])];
          currtan(1)=((ele->Tref())[i])(1) + disccol[BTSolDiscret().DofColMap()->LID(dofnode[7])];
          currtan(2)=((ele->Tref())[i])(2) + disccol[BTSolDiscret().DofColMap()->LID(dofnode[8])];
        }
      }
      // store into currenttangents
      currenttangents[node->Id()] = currtan;
    }
    //set currenttangents to zero for Lagrange centerline interpolation
    else
    {
      for (int i=0;i<3;i++)
        currenttangents[node->Id()](i) = 0.0;
    }
  }

  //**********************************************************************
  // update nodal coordinates also in existing contact pairs objects
  //**********************************************************************
  // loop over all pairs

  for(int i=0;i<(int)pairs_.size();++i)
  {
    // temporary matrices to store nodal coordinates of each element
    Epetra_SerialDenseMatrix ele1pos(3*numnodalvalues_,numnodes_);
    Epetra_SerialDenseMatrix ele2pos(3*numnodalvalues_,numnodes_);
    // Positions: Loop over all nodes of element 1
    /* be careful here: beam eles (such as beam3k, beam3r) may have intermediate
     * nodes which are not used for centerline interpolation and thus do not have
     * position or tangent DoFs. we therefore use numnodes_ here rather than the query
     * (btsphpotpairs_[i]->Element1())->NumNode() */
    // TODO do the same for beam-to-solid contact pairs
    for(int m=0;m<numnodes_;m++)
    {
      int tempGID = ((pairs_[i]->Element1())->NodeIds())[m];
      LINALG::Matrix<3,1> temppos = currentpositions[tempGID];

      // store updated nodal coordinates
      for(int n=0;n<3;n++)
        ele1pos(n,m) = temppos(n);
      // store updated nodal tangents
    }
    if (numnodalvalues_==2)
    {
      // Tangents: Loop over all nodes of element 1
      for(int m=0;m<numnodes_;m++)
      {
        int tempGID = ((pairs_[i]->Element1())->NodeIds())[m];
        LINALG::Matrix<3,1> temptan = currenttangents[tempGID];

        // store updated nodal tangents
        for(int n=0;n<3;n++)
          ele1pos(n+3,m) = temptan(n);
      }
    }
    // Positions: Loop over all nodes of element 2
    for(int m=0;m<numnodes_;m++)
    {
      int tempGID = ((pairs_[i]->Element2())->NodeIds())[m];
        LINALG::Matrix<3,1> temppos = currentpositions[tempGID];
      // store updated nodal coordinates
      for(int n=0;n<3;n++)
        ele2pos(n,m) = temppos(n);
    }
    if (numnodalvalues_==2)
    {
      // Tangents: Loop over all nodes of element 2
      for(int m=0;m<numnodes_;m++)
      {
        int tempGID = ((pairs_[i]->Element2())->NodeIds())[m];
          LINALG::Matrix<3,1> temptan = currenttangents[tempGID];

        // store updated nodal tangents
          for(int n=0;n<3;n++)
            ele2pos(n+3,m) = temptan(n);
      }
    }
    // finally update nodal positions in contact pair objects
    pairs_[i]->UpdateElePos(ele1pos,ele2pos);
  }
  // Update also the interpolated tangents if the tangentsmoothing is activated for Reissner beams
  int smoothing = DRT::INPUT::IntegralValue<INPAR::BEAMCONTACT::Smoothing>(sbeamcontact_,"BEAMS_SMOOTHING");
  if (smoothing != INPAR::BEAMCONTACT::bsm_none)
  {
    for(int i=0;i<(int)pairs_.size();++i)
    {
      pairs_[i]->UpdateEleSmoothTangents(currentpositions);
    }
  }

  // Do the same for the beam-to-solid contact pairs
  for(int i=0;i<(int)btsolpairs_.size();++i)
  {
    int numnodessol = ((btsolpairs_[i])->Element2())->NumNode();
    // temporary matrices to store nodal coordinates of each element
    Epetra_SerialDenseMatrix ele1pos(3*numnodalvalues_,numnodes_);
    Epetra_SerialDenseMatrix ele2pos(3,numnodessol);
    // Positions: Loop over all nodes of element 1 (beam element)
    for(int m=0;m<numnodes_;m++)
    {
      int tempGID = ((btsolpairs_[i]->Element1())->NodeIds())[m];
      LINALG::Matrix<3,1> temppos = currentpositions[tempGID];

      // store updated nodal coordinates
      for(int n=0;n<3;n++)
        ele1pos(n,m) = temppos(n);
      // store updated nodal tangents
    }
    if (numnodalvalues_==2)
    {
      // Tangents: Loop over all nodes of element 1
      for(int m=0;m<numnodes_;m++)
      {
        int tempGID = ((btsolpairs_[i]->Element1())->NodeIds())[m];
        LINALG::Matrix<3,1> temptan = currenttangents[tempGID];

        // store updated nodal tangents
        for(int n=0;n<3;n++)
          ele1pos(n+3,m) = temptan(n);
      }
    }
    // Positions: Loop over all nodes of element 2 (solid element)
    for(int m=0;m<(btsolpairs_[i]->Element2())->NumNode();m++)
    {
      int tempGID = ((btsolpairs_[i]->Element2())->NodeIds())[m];
        LINALG::Matrix<3,1> temppos = currentpositions[tempGID];
      // store updated nodal coordinates
      for(int n=0;n<3;n++)
        ele2pos(n,m) = temppos(n);
    }

    // finally update nodal positions in contact pair objects
    btsolpairs_[i]->UpdateElePos(ele1pos,ele2pos);
  }

  // Do the same for the beam-to-solid meshtying groups
  for(int i=0;i<(int)btsolmtgroups_.size();++i)
  {
    int numnodessol = ((btsolmtgroups_[i])->Element2()[0])->NumNode();
    int numele2 = (int)(btsolmtgroups_[i]->Element2().size());
    // temporary matrices to store nodal coordinates of each element
    Epetra_SerialDenseMatrix ele1pos(3*numnodalvalues_,numnodes_);
    std::vector<Epetra_SerialDenseMatrix> ele2pos(numele2);
    for (int e=0;e<numele2;e++)
      ele2pos[e].Shape(3,numnodessol);
    // Positions: Loop over all nodes of element 1 (beam element)
    for(int m=0;m<numnodes_;m++)
    {
      int tempGID = ((btsolmtgroups_[i]->Element1())->NodeIds())[m];
      LINALG::Matrix<3,1> temppos = currentpositions[tempGID];

      // store updated nodal coordinates
      for(int n=0;n<3;n++)
        ele1pos(n,m) = temppos(n);
      // store updated nodal tangents
    }
    if (numnodalvalues_==2)
    {
      // Tangents: Loop over all nodes of element 1
      for(int m=0;m<numnodes_;m++)
      {
        int tempGID = ((btsolmtgroups_[i]->Element1())->NodeIds())[m];
        LINALG::Matrix<3,1> temptan = currenttangents[tempGID];

        // store updated nodal tangents
        for(int n=0;n<3;n++)
          ele1pos(n+3,m) = temptan(n);
      }
    }
    // Positions: Loop over all nodes of all elements 2 (solid surface elements)
    for (int e=0;e<numele2;e++)
    {
      for(int m=0;m<(btsolmtgroups_[i]->Element2()[e])->NumNode();m++)
      {
        int tempGID = ((btsolmtgroups_[i]->Element2()[e])->NodeIds())[m];
        LINALG::Matrix<3,1> temppos = currentpositions[tempGID];
        // store updated nodal coordinates
        for(int n=0;n<3;n++)
          (ele2pos[e])(n,m) = temppos(n);
      }
    }

    // finally update nodal positions in meshtying group objects
    btsolmtgroups_[i]->UpdateElePos(ele1pos,ele2pos);
  }

  // Do the same for the beam-to-sphere contact pairs
  for(int i=0;i<(int)btsphpairs_.size();++i)
  {
    // temporary matrices to store nodal coordinates of each element
    Epetra_SerialDenseMatrix ele1pos(3*numnodalvalues_,numnodes_);
    Epetra_SerialDenseMatrix ele2pos(3,1);
    // Positions: Loop over all nodes of element 1 (beam element)
    for(int m=0;m<numnodes_;m++)
    {
      int tempGID = ((btsphpairs_[i]->Element1())->NodeIds())[m];
      LINALG::Matrix<3,1> temppos = currentpositions[tempGID];

      // store updated nodal coordinates
      for(int n=0;n<3;n++)
        ele1pos(n,m) = temppos(n);
      // store updated nodal tangents
    }
    if (numnodalvalues_==2)
    {
      // Tangents: Loop over all nodes of element 1
      for(int m=0;m<numnodes_;m++)
      {
        int tempGID = ((btsphpairs_[i]->Element1())->NodeIds())[m];
        LINALG::Matrix<3,1> temptan = currenttangents[tempGID];

        // store updated nodal tangents
        for(int n=0;n<3;n++)
          ele1pos(n+3,m) = temptan(n);
      }
    }
    // Positions: (rigid sphere element)
    int tempGID = ((btsphpairs_[i]->Element2())->NodeIds())[0];
      LINALG::Matrix<3,1> temppos = currentpositions[tempGID];
    // store updated nodal coordinates
    for(int n=0;n<3;n++)
      ele2pos(n,0) = temppos(n);

    // finally update nodal positions in contact pair objects
    btsphpairs_[i]->UpdateElePos(ele1pos,ele2pos);
  }

  // loop over all btbpotpairs
  for(int i=0;i<(int)btbpotpairs_.size();++i)
  {
    // temporary matrices to store nodal coordinates of each element
    Epetra_SerialDenseMatrix ele1pos(3*numnodalvalues_,numnodes_);
    Epetra_SerialDenseMatrix ele2pos(3*numnodalvalues_,numnodes_);
    // Positions: Loop over all nodes of element 1
    for(int m=0;m<numnodes_;m++)
    {
      int tempGID = ((btbpotpairs_[i]->Element1())->NodeIds())[m];
      LINALG::Matrix<3,1> temppos = currentpositions[tempGID];

      // store updated nodal coordinates
      for(int n=0;n<3;n++)
        ele1pos(n,m) = temppos(n);
      // store updated nodal tangents
    }
    if (numnodalvalues_==2)
    {
      // Tangents: Loop over all nodes of element 1
      for(int m=0;m<numnodes_;m++)
      {
        int tempGID = ((btbpotpairs_[i]->Element1())->NodeIds())[m];
        LINALG::Matrix<3,1> temptan = currenttangents[tempGID];

        // store updated nodal tangents
        for(int n=0;n<3;n++)
          ele1pos(n+3,m) = temptan(n);
      }
    }
    // Positions: Loop over all nodes of element 2
    for(int m=0;m<numnodes_;m++)
    {
      int tempGID = ((btbpotpairs_[i]->Element2())->NodeIds())[m];
        LINALG::Matrix<3,1> temppos = currentpositions[tempGID];
      // store updated nodal coordinates
      for(int n=0;n<3;n++)
        ele2pos(n,m) = temppos(n);
    }
    if (numnodalvalues_==2)
    {
      // Tangents: Loop over all nodes of element 2
      for(int m=0;m<numnodes_;m++)
      {
        int tempGID = ((btbpotpairs_[i]->Element2())->NodeIds())[m];
          LINALG::Matrix<3,1> temptan = currenttangents[tempGID];

        // store updated nodal tangents
          for(int n=0;n<3;n++)
            ele2pos(n+3,m) = temptan(n);
      }
    }
    // finally update nodal positions in contact pair objects
    btbpotpairs_[i]->UpdateElePos(ele1pos,ele2pos);
  }


  // Do the same for the beam-to-sphere potential pairs
  for(int i=0;i<(int)btsphpotpairs_.size();++i)
  {
    // temporary matrices to store nodal coordinates of each element
    Epetra_SerialDenseMatrix ele1pos(3*numnodalvalues_,numnodes_);
    Epetra_SerialDenseMatrix ele2pos(3,1);
    // Positions: Loop over all nodes of element 1 (beam element)
    for(int m=0;m<numnodes_;m++)
    {
      int tempGID = ((btsphpotpairs_[i]->Element1())->NodeIds())[m];
      LINALG::Matrix<3,1> temppos = currentpositions[tempGID];

      // store updated nodal coordinates
      for(int n=0;n<3;n++)
        ele1pos(n,m) = temppos(n);
      // store updated nodal tangents
    }
    if (numnodalvalues_==2)
    {
      // Tangents: Loop over all nodes of element 1
      for(int m=0;m<numnodes_;m++)
      {
        int tempGID = ((btsphpotpairs_[i]->Element1())->NodeIds())[m];
        LINALG::Matrix<3,1> temptan = currenttangents[tempGID];

        // store updated nodal tangents
        for(int n=0;n<3;n++)
          ele1pos(n+3,m) = temptan(n);
      }
    }
    // Positions: (rigid sphere element)
    int tempGID = ((btsphpotpairs_[i]->Element2())->NodeIds())[0];
      LINALG::Matrix<3,1> temppos = currentpositions[tempGID];
    // store updated nodal coordinates
    for(int n=0;n<3;n++)
      ele2pos(n,0) = temppos(n);

    // finally update nodal positions in contact pair objects
    btsphpotpairs_[i]->UpdateElePos(ele1pos,ele2pos);
  }

  return;
}

/*---------------------------------------------------------------------*
 |  Evaluate all pairs stored in the pair vectors            grill 10/14|
 *----------------------------------------------------------------------*/
void CONTACT::Beam3cmanager::EvaluateAllPairs(Teuchos::ParameterList timeintparams)
{
  //Begin: Determine maximal curvature occuring in complete beam discretization
  double kappa_max=0.0;
  global_kappa_max_=0.0;

  for (int i=0;i<ProblemDiscret().NumMyColElements();i++)
  {
    DRT::Element* element = ProblemDiscret().lColElement(i);
    const DRT::ElementType & eot = element->ElementType();
    if (eot == DRT::ELEMENTS::Beam3ebType::Instance())
    {
      const DRT::ELEMENTS::Beam3eb* beam3ebelement = dynamic_cast<const DRT::ELEMENTS::Beam3eb*>(element);

      if(fabs(beam3ebelement->GetKappaMax())>kappa_max)
        kappa_max=fabs(beam3ebelement->GetKappaMax());
    }
    else if (eot == DRT::ELEMENTS::Beam3rType::Instance())
    {
      const DRT::ELEMENTS::Beam3r* beam3relement = dynamic_cast<const DRT::ELEMENTS::Beam3r*>(element);

      if(fabs(beam3relement->GetKappaMax())>kappa_max)
        kappa_max=fabs(beam3relement->GetKappaMax());
    }
    else
    {
      //std::cout << "Warning: Calculation of kappa_max only implemented for beam3eb elements so far!" << std::endl;
      kappa_max=0.0;
    }
  }

  Comm().MaxAll(&kappa_max,&global_kappa_max_,1);
//  std::cout << "global_kappa_max_: " << global_kappa_max_ << std::endl;
  timeintparams.set("kappa_max", global_kappa_max_);
  //End: Determine maximal curvature occuring in complete beam discretization

  // Loop over all BTB contact pairs
  for (int i=0;i<(int)pairs_.size();++i)
  {
    // only evaluate pair for those procs owning or ghosting at
    // least one node of one of the two elements of the pair
    int firsteleid = (pairs_[i]->Element1())->Id();
    int secondeleid = (pairs_[i]->Element2())->Id();
    bool firstisincolmap = ColElements()->MyGID(firsteleid);
    bool secondisincolmap = ColElements()->MyGID(secondeleid);

    // evaluate additional contact forces and stiffness
    if (firstisincolmap || secondisincolmap)
    {
      pairs_[i]->Evaluate(*stiffc_,*fc_,currentpp_,contactpairmap_,timeintparams);

      //if active, get minimal gap of contact element pair
      if(pairs_[i]->GetContactFlag() == true)
      {
        std::vector<double> pairgaps = pairs_[i]->GetGap();
        for(int i=0;i<(int)pairgaps.size();i++)
        {
          double gap = pairgaps[i];
          if(gap<mintotalsimunconvgap_)
            mintotalsimunconvgap_ = gap;
        }
      }
    }
  }

  // Loop over all BTSOL contact pairs
  for (int i=0;i<(int)btsolpairs_.size();++i)
  {
    // only evaluate pair for those procs owning or ghosting at
    // least one node of one of the two elements of the pair
    int firsteleid = (btsolpairs_[i]->Element1())->Id();
    int secondeleid = (btsolpairs_[i]->Element2())->Id();
    bool firstisincolmap = ColElements()->MyGID(firsteleid);
    bool secondisincolmap = ColElements()->MyGID(secondeleid);
    // evaluate additional contact forces and stiffness
    if (firstisincolmap || secondisincolmap)
    {
      btsolpairs_[i]->Evaluate(*stiffc_,*fc_,btspp_);
    }
  }

  // Loop over all BTSOL meshtying groups
  for (int i=0;i<(int)btsolmtgroups_.size();++i)
  {
    // only evaluate pair for those procs owning or ghosting at
    // least one node of one of the elements of the pair
    int firsteleid = (btsolmtgroups_[i]->Element1())->Id();
    bool firstisincolmap = ColElements()->MyGID(firsteleid);
    bool secondisincolmap = false;
    std::vector<DRT::Element*> ele2 = btsolmtgroups_[i]->Element2();
    for (int m=0;m<(int)ele2.size();++m)
    {
      int secondeleid = ele2[m]->Id();
      secondisincolmap = secondisincolmap || ColElements()->MyGID(secondeleid);
    }
    // evaluate additional meshtying forces and stiffness
    if (firstisincolmap || secondisincolmap)
    {
      btsolmtgroups_[i]->Evaluate(*stiffc_,*fc_,btsmtpp_);
    }
  }

  // Loop over all BTSPH contact pairs
  for (int i=0;i<(int)btsphpairs_.size();++i)
  {
    // only evaluate pair for those procs owning or ghosting at
    // least one node of one of the two elements of the pair
    int firsteleid = (btsphpairs_[i]->Element1())->Id();
    int secondeleid = (btsphpairs_[i]->Element2())->Id();
    bool firstisincolmap = ColElements()->MyGID(firsteleid);
    bool secondisincolmap = ColElements()->MyGID(secondeleid);
    // evaluate additional contact forces and stiffness
    if (firstisincolmap || secondisincolmap)
    {
      btsphpairs_[i]->Evaluate(*stiffc_,*fc_,btsphpp_);
    }
  }

  // Loop over all BTB potential pairs
  for (int i=0;i<(int)btbpotpairs_.size();++i)
  {
    // only evaluate pair for those procs owning or ghosting at
    // least one node of one of the two elements of the pair
    int firsteleid = (btbpotpairs_[i]->Element1())->Id();
    int secondeleid = (btbpotpairs_[i]->Element2())->Id();
    bool firstisincolmap = ColElements()->MyGID(firsteleid);
    bool secondisincolmap = ColElements()->MyGID(secondeleid);
    // evaluate additional contact forces and stiffness
    if (firstisincolmap || secondisincolmap)
    {
      for (unsigned int j=0; j<ki_->size(); ++j)
      {
        if (ki_->at(j)!=0.0)
          btbpotpairs_[i]->Evaluate(*stiffc_,*fc_,ki_->at(j),mi_->at(j));
      }
    }
  }

  // Loop over all BTSPH potential pairs
  for (int i=0;i<(int)btsphpotpairs_.size();++i)
  {
    // only evaluate pair for those procs owning or ghosting at
    // least one node of one of the two elements of the pair
    int firsteleid = (btsphpotpairs_[i]->Element1())->Id();
    int secondeleid = (btsphpotpairs_[i]->Element2())->Id();
    bool firstisincolmap = ColElements()->MyGID(firsteleid);
    bool secondisincolmap = ColElements()->MyGID(secondeleid);
    // evaluate additional contact forces and stiffness
    if (firstisincolmap || secondisincolmap)
    {
      for (unsigned int j=0; j<ki_->size(); ++j)
      {
        if (ki_->at(j)!=0.0)
          btsphpotpairs_[i]->Evaluate(*stiffc_,*fc_,ki_->at(j),mi_->at(j));
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 |  process the found element pairs and fill the corresponding
 |  BTB, BTSOL and BTSPH contact pair vectors                grill 09/14|
 *----------------------------------------------------------------------*/
void CONTACT::Beam3cmanager::FillContactPairsVectors(const std::vector<std::vector<DRT::Element* > > elementpairs)
{
  std::vector<std::vector<DRT::Element* > > formattedelementpairs;
  formattedelementpairs.clear();

  //Besides beam-to-beam contact we also can handle beam-to-solid contact and beam to sphere contact(will be implemented in the future).
  //In all cases element 1 has to be the beam element.

  // ALWAYS START FROM SCRATCH with beam-to-solid meshtying groups
  std::map<int, std::vector<int> > localbtsmtgroups;
  btsolmtgroups_.clear();
  btsolmtgroups_.resize(0);

  //All other element pairs (solid-solid, sphere-solid etc.) will be sorted out later.
  for (int i=0;i<(int)elementpairs.size();i++)
  {
    // if ele1 is a beam element we take the pair directly
    if(BEAMCONTACT::BeamElement(*(elementpairs[i])[0]))
    {
      formattedelementpairs.push_back(elementpairs[i]);
    }
    // if ele1 is no beam element, but ele2 is one, we have to change the order
    else if(BEAMCONTACT::BeamElement(*(elementpairs[i])[1]))
    {
      std::vector<DRT::Element* > elementpairaux;
      elementpairaux.clear();
      elementpairaux.push_back((elementpairs[i])[1]);
      elementpairaux.push_back((elementpairs[i])[0]);
      formattedelementpairs.push_back(elementpairaux);
    }
  }
  //Determine type of applied beam elements and set the corresponding values for the member variables
  //numnodes_ and numnodaldofs_. This has only to be done once in the beginning, since beam contact
  //simulations are only possible when using beam elements of one type!
  if (!elementtypeset_ and formattedelementpairs.size()>0)
  {
    SetElementTypeAndDistype((formattedelementpairs[0])[0]);
    elementtypeset_=true;
  }

  //So far, all beam elements occurring in the pairs_ vector have to be of the same type.
  //This will be checked in the following lines.
  if ((int)formattedelementpairs.size()>0)
  {
    //search for the first beam element in vector pairs
    const DRT::ElementType & pair1_ele1_type = ((formattedelementpairs[0])[0])->ElementType();
    for (int k=0;k<(int)formattedelementpairs.size();++k)
    {
      const DRT::ElementType & ele1_type = ((formattedelementpairs[k])[0])->ElementType();
      const DRT::ElementType & ele2_type = ((formattedelementpairs[k])[1])->ElementType();

      //ele1 and ele2 (in case this is a beam element) have to be of the same type as ele1 of the first pair
      if ( ele1_type!=pair1_ele1_type or (BEAMCONTACT::BeamElement(*(formattedelementpairs[k])[1]) and ele2_type!=pair1_ele1_type))
      {
        dserror("All contacting beam elements have to be of the same type (beam3k, beam3eb, beam3 or beam3r). Check your input file!");
      }
    }
  }
  // Only the element pairs of formattedelementpairs (found in the contact search) which have not been found in the last time step (i.e. which
  // are not in oldpairs_) will be generated as new Beam3contact instances. Pairs which already exist in oldpairs_ will simply be copied to pairs_.
  // This procedure looks a bit circumstantial at first glance: However, it is not possible to solely use formattedelementpairs
  // and to simply delete oldpairs_ at the end of a time step if the new gap function definition is used, since the latter needs
  // history variables of the last time step which are stored in the oldpairs_ vector.
  // Only beam-to-beam contact pairs (not beam-to-solid or beam-to-sphere pairs) need this history information.
  for (int k=0;k<(int)formattedelementpairs.size();k++)
  {
    DRT::Element* ele1 = (formattedelementpairs[k])[0];
    DRT::Element* ele2 = (formattedelementpairs[k])[1];
    int currid1 = ele1->Id();
    int currid2 = ele2->Id();

    //beam-to-beam pair
    if(BEAMCONTACT::BeamElement(*(formattedelementpairs[k])[1]))
    {
      bool foundlasttimestep = false;
      bool isalreadyinpairs = false;

      if(contactpairmap_.find(std::make_pair(currid1,currid2))!=contactpairmap_.end())
        isalreadyinpairs = true;

      if(oldcontactpairmap_.find(std::make_pair(currid1,currid2))!=oldcontactpairmap_.end())
        foundlasttimestep = true;

      if(!isalreadyinpairs and foundlasttimestep)
      {
        pairs_.push_back(oldcontactpairmap_[std::make_pair(currid1,currid2)]);
        if (currid1<currid2)
          contactpairmap_[std::make_pair(currid1,currid2)] = pairs_[pairs_.size()-1];
        else
          dserror("Element 1 has to have the smaller element-ID. Adapt your contact search!");

        isalreadyinpairs = true;
      }

      if(!isalreadyinpairs)
      {
        //Add new contact pair object: The auxiliary_instance of the abstract class Beam3contactinterface is only needed here in order to call
        //the function Impl() which creates an instance of the templated class Beam3contactnew<numnodes, numnodalvalues> !
        pairs_.push_back(CONTACT::Beam3contactinterface::Impl(numnodes_,numnodalvalues_,ProblemDiscret(),BTSolDiscret(),dofoffsetmap_,ele1,ele2, sbeamcontact_));
        if (currid1<=currid2)
          contactpairmap_[std::make_pair(currid1,currid2)] = pairs_[pairs_.size()-1];
        else
          dserror("Element 1 has to have the smaller element-ID. Adapt your contact search!");
      }
    }
    // beam-to-sphere pair
    else if(BEAMCONTACT::RigidsphereElement(*(formattedelementpairs[k])[1]))
    {
      if(btsphpairmap_.find(std::make_pair(currid1,currid2))==btsphpairmap_.end())
      {
        btsphpairs_.push_back(CONTACT::Beam3tospherecontactinterface::Impl(numnodes_,numnodalvalues_,ProblemDiscret(),BTSolDiscret(),dofoffsetmap_,ele1,ele2));
        btsphpairmap_[std::make_pair(currid1,currid2)] = btsphpairs_[btsphpairs_.size()-1];
      }
    }
    // beam-to-solid contact pair
    else if(BEAMCONTACT::SolidContactElement(*(formattedelementpairs[k])[1]))
    {
      bool foundlasttimestep = false;
      bool isalreadyinpairs = false;

      if(btsolpairmap_.find(std::make_pair(currid1,currid2))!=btsolpairmap_.end())
        isalreadyinpairs = true;

      if(oldbtsolpairmap_.find(std::make_pair(currid1,currid2))!=oldbtsolpairmap_.end())
        foundlasttimestep = true;

      if(!isalreadyinpairs and foundlasttimestep)
      {
        btsolpairs_.push_back(oldbtsolpairmap_[std::make_pair(currid1,currid2)]);
        if (currid1<currid2)
          btsolpairmap_[std::make_pair(currid1,currid2)] = btsolpairs_[btsolpairs_.size()-1];
        else
          dserror("Element 1 has to have the smaller element-ID. Adapt your contact search!");

        isalreadyinpairs = true;
      }

      if(!isalreadyinpairs)
      {
        //Add new contact pair object: The auxiliary_instance of the abstract class Beam3contactinterface is only needed here in order to call
        //the function Impl() which creates an instance of the templated class Beam3contactnew<numnodes, numnodalvalues> !
        btsolpairs_.push_back(CONTACT::Beam3tosolidcontactinterface::Impl((formattedelementpairs[k])[1]->NumNode(),numnodes_,numnodalvalues_,ProblemDiscret(),BTSolDiscret(),dofoffsetmap_,ele1,ele2, sbeamcontact_));
        if (currid1<=currid2)
          btsolpairmap_[std::make_pair(currid1,currid2)] = btsolpairs_[btsolpairs_.size()-1];
        else
          dserror("Element 1 has to have the smaller element-ID. Adapt your contact search!");
      }
    }
    // beam-to-solid meyhtying pair
    else if(BEAMCONTACT::SolidMeshtyingElement(*(formattedelementpairs[k])[1]))
    {
     bool isalreadyinpairs = false;

     if(localbtsmtgroups.find(currid1)!=localbtsmtgroups.end())
     {
       std::vector<int> master = localbtsmtgroups[currid1];
       for (int m=0;m<(int)master.size();++m)
       {
         if (master[m]==currid2)
           isalreadyinpairs=true;
       }
     }

     if (!isalreadyinpairs)
     {
       if (currid1<=currid2)
         localbtsmtgroups[currid1].push_back(currid2);
       else
         dserror("Element 1 has to have the smaller element-ID. Adapt your contact search!");
     }
    }
    else
    {
      dserror("ERROR: Unknown element type in beam contact pairs (none of BTB, BTSco, BTSmt, BTSPH)");
    }
  }

  // finish construction of beam-to-solid meshtying groups (this has to be done here after having collected
  // the groups of master (solid surface) elements belonging to one slave (beam) element
  std::map<int,std::vector<int> >::iterator iter;
  for (iter=localbtsmtgroups.begin();iter!=localbtsmtgroups.end();++iter)
  {
    // ID of element 1 (slave = beam)
    int ele1id = iter->first;

    // IDs of list of elements 2 (master = solid surface)
    std::vector<int> ele2id = iter->second;

    // get element 1
    DRT::Element* ele1 = BTSolDiscret().gElement(ele1id);

    // get list of elements 2
    std::vector<DRT::Element*> ele2;
    for (int m=0;m<(int)ele2id.size();++m)
      ele2.push_back(BTSolDiscret().gElement(ele2id[m]));

    // currently all solid elements must be of same type...
    int numnodessol = 0;
    for (int m=0;m<(int)ele2.size();++m)
    {
      if (m==0)
        numnodessol = ele2[m]->NumNode();
      else
      {
        int temp = ele2[m]->NumNode();
        if (temp != numnodessol)
          dserror("ERROR: All BTS meshtying solid elements must still be of same type");
      }
    }

    // add new meyhtying group object: The auxiliary_instance of the abstract class Beam3tosolidmeyhtyinginterface
    // is only needed here in order to call the function Impl() which creates an instance of the templated class
    // Beam3tosolidmeshtying<numnodessol,numnodes, numnodalvalues>s
    btsolmtgroups_.push_back(CONTACT::Beam3tosolidmeshtyinginterface::Impl(numnodessol,numnodes_,numnodalvalues_,
        ProblemDiscret(),BTSolDiscret(),dofoffsetmap_,ele1,ele2, sbeamcontact_));
  }

  if(pdiscret_.Comm().MyPID()==0)
  {
    std::cout << "\t Total number of BTB contact pairs:     " << pairs_.size() << std::endl;
    if (btsph_)   std::cout << "\t Total number of BTSPH contact pairs:    " << btsphpairs_.size() << std::endl;
    if (btsol_)   std::cout << "\t Total number of BTSOL contact pairs:    " << btsolpairs_.size() << std::endl;
    if (btsolmt_) std::cout << "\t Total number of BTSOL meshtying groups: " << btsolmtgroups_.size() << std::endl;
  }
}

/*----------------------------------------------------------------------*
 |  process the found element pairs and fill the corresponding
 |  BTB, BTSOL and BTSPH potential pair vectors              grill 09/14|
 *----------------------------------------------------------------------*/
void CONTACT::Beam3cmanager::FillPotentialPairsVectors(const std::vector<std::vector<DRT::Element* > > elementpairs)
{
  std::vector<std::vector<DRT::Element* > > formattedelementpairs;
  formattedelementpairs.clear();

  //Besides beam-to-beam potentials, we can also handle beam to sphere potentials
  //In all cases element 1 has to be the beam element.

  //All other element pairs (sphere-sphere, solid-solid, sphere-solid etc.) will be sorted out later.
  for (int i=0;i<(int)elementpairs.size();i++)
  {
    // if ele1 is a beam element we take the pair directly
    if(BEAMCONTACT::BeamElement(*(elementpairs[i])[0]))
    {
      formattedelementpairs.push_back(elementpairs[i]);
    }
    // if ele1 is no beam element, but ele2 is one, we have to change the order
    else if(BEAMCONTACT::BeamElement(*(elementpairs[i])[1]))
    {
      std::vector<DRT::Element* > elementpairaux;
      elementpairaux.clear();
      elementpairaux.push_back((elementpairs[i])[1]);
      elementpairaux.push_back((elementpairs[i])[0]);
      formattedelementpairs.push_back(elementpairaux);
    }
  }
  //Determine type of applied beam elements and set the corresponding values for the member variables
  //numnodes_ and numnodaldofs_. This has only to be done once in the beginning, since beam contact
  //simulations are only possible when using beam elements of one type!
  if (numnodalvalues_==0 and formattedelementpairs.size()>0)
    SetElementTypeAndDistype((formattedelementpairs[0])[0]);

  // Create Beam3tobeampotentialinterface instances in the btbpotpairs_ vector
  for (int k=0;k<(int)formattedelementpairs.size();k++)
  {
    DRT::Element* ele1 = (formattedelementpairs[k])[0];
    DRT::Element* ele2 = (formattedelementpairs[k])[1];
    int currid1 = ele1->Id();
    int currid2 = ele2->Id();

    // check the line charge conditions which apply for the nodes of ele1 and ele2
    // find and pass line charge conditions associated with the elements of this pair
    std::vector<DRT::Condition*> currconds;
    currconds.clear();

    bool nocharge = false;
    // for now, we exclude mutual interaction of elements on same beam (same DLINE)
    // TODO read flag from Inputfile whether to do this or not
    bool samedesignline = false;
    // check arbitrary (here: first) node of element for the DLINE, it is part of
    if (dlinenodemap_.find(ele1->NodeIds()[0]) != dlinenodemap_.end() and
        dlinenodemap_.find(ele2->NodeIds()[0]) != dlinenodemap_.end())
    {
      if (dlinenodemap_[ele1->NodeIds()[0]] == dlinenodemap_[ele2->NodeIds()[0]])
      {
        samedesignline = true;
      }
      else
      {
        // find and store line charge condition of first node of element 1/2
        currconds.push_back(linechargeconds_[dlinenodemap_[ele1->NodeIds()[0]]]);
        currconds.push_back(linechargeconds_[dlinenodemap_[ele2->NodeIds()[0]]]);
      }
    }
    else  // none of the elements is "loaded" by a charge -> do not create a btbpotpair
    {
      nocharge = true;
    }

    bool isalreadyinpotpairs = false;
    if(BEAMCONTACT::BeamElement(*(formattedelementpairs[k])[1]))
    {
      // TODO use a potpairmap_ for this query. see contactpairmap
      for (unsigned int i=0; i<btbpotpairs_.size(); ++i)
      {
        int id1 = btbpotpairs_[i]->Element1()->Id();
        int id2 = btbpotpairs_[i]->Element2()->Id();

        if ((id1==currid1 and id2==currid2) or (id1==currid2 and id2==currid1))
          isalreadyinpotpairs = true;
      }
    }
    else if (BEAMCONTACT::RigidsphereElement(*(formattedelementpairs[k])[1]))
    {
      // TODO use a potpairmap_ for this query. see contactpairmap
      for (unsigned int i=0; i<btsphpotpairs_.size(); ++i)
      {
        int id1 = btsphpotpairs_[i]->Element1()->Id();
        int id2 = btsphpotpairs_[i]->Element2()->Id();

        if ((id1==currid1 and id2==currid2) or (id1==currid2 and id2==currid1))
          isalreadyinpotpairs = true;
      }
    }

    if((!isalreadyinpotpairs) and (!samedesignline) and (!nocharge))
    {
      //beam-to-beam pair
      if(BEAMCONTACT::BeamElement(*(formattedelementpairs[k])[1]))
      {
        //Add new potential pair object: The auxiliary_instance of the abstract class Beam3tobeampotentialinterface is only needed here in order to call
        //the function Impl() which creates an instance of the templated class Beam3tobeampotential<numnodes, numnodalvalues> !
        btbpotpairs_.push_back(CONTACT::Beam3tobeampotentialinterface::Impl(numnodes_,numnodalvalues_,ProblemDiscret(),BTSolDiscret(),dofoffsetmap_,ele1,ele2, sbeampotential_,currconds));
      }
      // beam-to-sphere pair
      else if (BEAMCONTACT::RigidsphereElement(*(formattedelementpairs[k])[1]) and potbtsph_)
      {
        //Add new potential pair object: The auxiliary_instance of the abstract class Beam3tobeampotentialinterface is only needed here in order to call
        //the function Impl() which creates an instance of the templated class Beam3tobeampotential<numnodes, numnodalvalues> !
        btsphpotpairs_.push_back(CONTACT::Beam3tospherepotentialinterface::Impl(numnodes_,numnodalvalues_,ProblemDiscret(),BTSolDiscret(),dofoffsetmap_,ele1,ele2, sbeampotential_,currconds));
      }
      // beam-to-solid pair, beam-to-??? pair
      else
      {
         dserror("Only beam-to-beam potential interaction is implemented yet. No other types of elements allowed!");
      }
    }
  }

  if(pdiscret_.Comm().MyPID()==0)
  {
    std::cout << "            Total number of BTB pot pairs: " << btbpotpairs_.size() << std::endl;
    if (potbtsph_)
      std::cout << "            Total number of BTSPH pot pairs: " << btsphpotpairs_.size() << std::endl;
  }
}

/*----------------------------------------------------------------------*
 |  search possible contact element pairs                     popp 04/10|
 *----------------------------------------------------------------------*/
std::vector<std::vector<DRT::Element*> > CONTACT::Beam3cmanager::BruteForceSearch(std::map<int,LINALG::Matrix<3,1> >& currentpositions, const double searchradius, const double sphericalsearchradius)
{
  //**********************************************************************
  // Steps of search for element pairs that MIGHT get into contact:
  //
  // 1) Find non-neighboring node pairs
  // 2) Compute distance between node pairs and compare with search radius
  // 3) Find non-neighboring element pairs based on node pairs
  // 4) Check if new pair already exists. If not we've found a new entry!
  //
  // NOTE: This is a brute force search! The time to search across n nodes
  // goes with n^2, which is not efficient at all....
  //**********************************************************************
  std::vector<std::vector<DRT::Element*> > newpairs;
  newpairs.clear();

  //**********************************************************************
  // LOOP 1: column nodes (overlap = 1)
  // each processor looks for close nodes (directly connect with the corresponding node) for each of these nodes
  //**********************************************************************
  for (int i=0;i<ColNodes()->NumMyElements();++i)
  {
    // get global id, node itself and current position
    int firstgid = ColNodes()->GID(i);
    DRT::Node* firstnode = BTSolDiscret().gNode(firstgid);

    // TODO see also LOOP 2 below
   /* check whether node position has been stored in currentpositions previously;
    * if not, it most likely is a beam node which is NOT used for centerline interpolation
    * if so, we simply skip it because it does not have position (and tangent) DoFs */
    if (currentpositions.find(firstgid) == currentpositions.end())
    {
      if (BEAMCONTACT::BeamNode(*firstnode) and !BEAMCONTACT::BeamCenterlineNode(*firstnode))
      {
        continue;
      }
      else
      {
        dserror("this should not happen!");
      }
    }

    LINALG::Matrix<3,1> firstpos = currentpositions[firstgid];

    // create storage for neighbouring nodes to be excluded.
    std::vector<int> neighbournodeids(0);
    // create storage for near nodes to be identified
    std::vector<int> NearNodesGIDs;

    // get the elements 'firstnode' is linked to
    DRT::Element** neighboureles = firstnode->Elements();
    // loop over all adjacent elements and their nodes
    for (int j=0;j<firstnode->NumElement();++j)
    {
      DRT::Element* thisele = neighboureles[j];
      for (int k=0;k<thisele->NumNode();++k)
      {
        int nodeid = thisele->NodeIds()[k];
        if (nodeid == firstgid) continue;

        // add to neighbouring node vector
        neighbournodeids.push_back(nodeid);
      }
    }
    //**********************************************************************
    // LOOP 2: all nodes (fully overlapping column map)
    // each processor looks for close nodes within these nodes
    //**********************************************************************
    for (int j=0;j<FullNodes()->NumMyElements();++j)
    {
      // get global node id and current position
      int secondgid = FullNodes()->GID(j);

      // TODO see comment above
      if (currentpositions.find(secondgid) == currentpositions.end())
        continue;

      LINALG::Matrix<3,1> secondpos = currentpositions[secondgid];

      // nothing to do for identical pair
      if (firstgid == secondgid) continue;

      // check if second node is neighbouring node
      bool neighbouring = false;
      for (int k=0;k<(int)neighbournodeids.size();++k)
      {
        if (secondgid==neighbournodeids[k])
        {
          neighbouring = true;
          break;
        }
      }
      // compute distance by comparing firstpos <-> secondpos
      if (neighbouring==false)
      {
        LINALG::Matrix<3,1> distance;
        for (int k=0;k<3;k++) distance(k) = secondpos(k)-firstpos(k);

        // nodes are near if distance < search radius
        if (distance.Norm2() < searchradius)
          NearNodesGIDs.push_back(secondgid);
      }
    }
    // AT THIS POINT WE HAVE FOUND AND STORED ALL NODES CLOSE TO FIRSTNODE BESIDES THE DIRECTLY CONNECTED NEIGHBOR NODES!

    //*********************************************************************
    // Up to now we have only considered nodes, but in the end we will need
    // to find element pairs, that might possibly get into contact. To find
    // these element pairs, we combine the elements around 'firstnode' with
    // all elements around each 'NearNodesGIDs'-node. Repetitions of pairs
    // and neighboring pairs will be rejected. For the remaining GIDs,
    // pointers on these elements are be created. With these pointers, the
    // beam3contact objects are set up and stored into the vector pairs_.
    //*********************************************************************
    // vectors of element ids
    std::vector<int> FirstElesGIDs(0);
    std::vector<int> SecondElesGIDs(0);
    // loop over all elements adjacent to firstnode
    for (int j=0;j<firstnode->NumElement();++j)
    {
      //element pointer
      DRT::Element* ele1 = neighboureles[j];

      // insert into element vector
      FirstElesGIDs.push_back(ele1->Id());
    }

    // loop over ALL nodes close to first node
    for (int j=0;j<(int)NearNodesGIDs.size();++j)
    {
      // node pointer
      DRT::Node* tempnode = BTSolDiscret().gNode(NearNodesGIDs[j]);
      // getting the elements tempnode is linked to
      DRT::Element** TempEles = tempnode->Elements();

      // loop over all elements adjacent to tempnode
      for (int k=0;k<tempnode->NumElement();++k)
      {
        //element pointer
        DRT::Element* ele2 = TempEles[k];
        // insert into element vector
        SecondElesGIDs.push_back(ele2->Id());
      }
    }
    // AT THIS POINT WE HAVE FOUND AND STORED ALL ELEMENTS CLOSE TO FIRSTNODE!

    //*********************************************************************
    // The combination and creation of beam3contact objects can now begin.
    // First of all we reject all second element GIDs, that occur twice and
    // generate a new vector 'SecondElesGIDsRej' where each GID occurs only
    // once. This vector will then be used for generating the pair objects.
    //*********************************************************************

    // initialize reduced vector of close elements
    std::vector<int> SecondElesGIDsRej;
    SecondElesGIDsRej.clear();

    // loop over all close elements
    for (int j=0;j<(int)SecondElesGIDs.size();++j)
    {
      // check if this element gid occurs twice
      int tempGID = SecondElesGIDs[j];
      bool twice = false;

      // loop over all close elements again
      for (int k=j+1;k<(int)SecondElesGIDs.size();++k)
      {
        if (tempGID==SecondElesGIDs[k]) twice = true;
      }

      // only insert in reduced vector if not yet there
      if (twice == false) SecondElesGIDsRej.push_back(tempGID);
    }

    // now finally create element pairs via two nested loops
    for (int j=0;j<(int)FirstElesGIDs.size();++j)
    {
      // beam element pointer
      DRT::Element* ele1 = BTSolDiscret().gElement(FirstElesGIDs[j]);

      // node ids adjacent to this element
      const int* NodesEle1 = ele1->NodeIds();

      // loop over all close elements
      for (int k=0;k<(int)SecondElesGIDsRej.size();++k)
      {
        // get and cast a pointer on an element
        DRT::Element* ele2 = BTSolDiscret().gElement(SecondElesGIDsRej[k]);

        // close element id
        const int* NodesEle2 = ele2->NodeIds();

        // check if elements are neighbouring (share one common node)
        bool elements_neighbouring = false;
        for (int m=0;m<ele1->NumNode();++m)
        {
          for (int n=0;n<ele2->NumNode();++n)
          {
            // neighbouring if they share one common node
            if(NodesEle1[m] == NodesEle2[n])
              elements_neighbouring = true;
          }
        }

        // Check if the pair 'jk' already exists in newpairs
        bool foundbefore = false;
        for (int m=0;m<(int)newpairs.size();++m)
        {
          int id1 = (newpairs[m])[0]->Id();
          int id2 = (newpairs[m])[1]->Id();
          int currid1 = FirstElesGIDs[j];
          int currid2 = SecondElesGIDsRej[k];

          // already exists if can be found in newpairs
          if ((id1 == currid1 && id2 == currid2) || (id1 == currid2 && id2 == currid1))
            foundbefore = true;
        }

        // if NOT neighbouring and NOT found before
        // create new beam3contact object and store it into pairs_

        // Here we additonally apply the method CloseMidpointDistance which sorts out all pairs with a midpoint distance
        // larger than sphericalsearchradius. Thus with this additional method the search is based on spherical bounding boxes
        // and node on node distances any longer. The radius of these spheres is sphericalsearchradius/2.0, the center of such
        // a sphere is (r1+r2)/2, with r1 and r2 representing the nodal positions.
        if (!elements_neighbouring && !foundbefore && CloseMidpointDistance(ele1, ele2, currentpositions, sphericalsearchradius))
        {
          std::vector<DRT::Element*> contactelementpair;
          contactelementpair.clear();
          if (ele1->Id()<ele2->Id())
          {
            contactelementpair.push_back(ele1);
            contactelementpair.push_back(ele2);
          }
          else
          {
            contactelementpair.push_back(ele2);
            contactelementpair.push_back(ele1);
          }
          newpairs.push_back(contactelementpair);
        }
      }
    }
  } //LOOP 1
  return newpairs;
}

/*-------------------------------------------------------------------- -*
 |  Compute search radius from discretization data            popp 04/10|
 *----------------------------------------------------------------------*/
void CONTACT::Beam3cmanager::ComputeSearchRadius()
{
  // some local variables
  double charactlength = 0.0;
  double globalcharactlength = 0.0;
  double maxelelength = 0.0;

  // look for maximum element length in the whole discretization
  GetMaxEleLength(maxelelength);

  // select characeteristic length
  if (maxeleradius_ > maxelelength) charactlength = maxeleradius_;
  else                             charactlength = maxelelength;

  // communicate among all procs to find the global maximum
  Comm().MaxAll(&charactlength,&globalcharactlength,1);

  // Compute the search radius. This one is only applied to determine
  // close pairs considering the node-to-node distances.
  double nodalsearchfac = 3.0;
  searchradius_ = nodalsearchfac * (2.0*searchboxinc_+globalcharactlength);

  // In a second step spherical search boxes are applied which consider
  // the midpoint-to-midpiont distance. In the first (nodal-based) search step
  // it has to be ensured that all pairs relevant for this second search step will
  // be found. The most critical case (i.e. the case, where the midpoints are as close as
  // possible but the node distances are as large as possible) is the case where to (straight)
  // beams are perpendicular to each other and the beam midpoint coincide with the closest points
  // between these two beams. One can show, that in this case a value of nodalsearchfac=2.0 is sufficient
  // to find all relevant pairs in the first step. This factor should also be sufficient, if the two beam
  // elements are deformed (maximal assumed deformation of a beam element is a half circle!). To be on
  // the safe side (the number of elements found in the first search step is not very relevant for the
  // overall efficiency), we choose a factor of nodalsearchfac = 3.0.
  sphericalsearchradius_ = 2.0*searchboxinc_ + globalcharactlength;

  // some information for the user
  if (Comm().MyPID()==0)
  {
    std::cout << "Penalty parameter      = " << currentpp_ << std::endl;
    std::cout << "BTS-Penalty parameter  = " << btspp_ << std::endl;
    std::cout << "Maximum element radius = " << maxeleradius_ << std::endl;
    std::cout << "Maximum element length = " << maxelelength << std::endl;
    std::cout << "Search radius          = " << searchradius_  << std::endl << std::endl;
  }

  return;
}

/*-------------------------------------------------------------------- -*
 |  Find maximum element radius in discretization             popp 04/10|
 *----------------------------------------------------------------------*/
void CONTACT::Beam3cmanager::SetMinMaxEleRadius()
{

  mineleradius_=0.0;
  maxeleradius_=0.0;

  bool minbeamradiusinitialized=false;

  // loop over all elements in row map
  for (int i=0;i<RowElements()->NumMyElements();++i)
  {
    // get pointer onto element
    int gid = RowElements()->GID(i);
    DRT::Element* thisele = BTSolDiscret().gElement(gid);

    double eleradius = 0.0;

    if (BEAMCONTACT::BeamElement(*thisele) or BEAMCONTACT::RigidsphereElement(*thisele))
    {  // compute eleradius from moment of inertia
      // (RESTRICTION: CIRCULAR CROSS SECTION !!!)
      eleradius = BEAMCONTACT::CalcEleRadius(thisele);

      // if current radius is larger than maximum radius -> update
      if (eleradius > maxeleradius_) maxeleradius_ = eleradius;

      //Initialize minbeamradius- with the first radius we get; otherwise its value will remain 0.0!
      if (!minbeamradiusinitialized)
      {
        mineleradius_=eleradius;
        minbeamradiusinitialized=true;
      }

      // if current radius is smaller than minimal radius -> update
      if (eleradius < mineleradius_) mineleradius_ = eleradius;
    }
  }

  return;
}

/*-------------------------------------------------------------------- -*
 |  Find maximum element length in discretization             popp 04/10|
 *----------------------------------------------------------------------*/
void CONTACT::Beam3cmanager::GetMaxEleLength(double& maxelelength)
{
  // loop over all elements in row map
  for (int i=0;i<RowElements()->NumMyElements();i++)
  {
    // get pointer onto element
    int gid = RowElements()->GID(i);
    DRT::Element* thisele = BTSolDiscret().gElement(gid);

    double elelength = 0.0;

    if (BEAMCONTACT::BeamElement(*thisele))
    {
      // get global IDs of edge nodes and pointers
      int node0_gid = thisele->NodeIds()[0];
      int node1_gid = thisele->NodeIds()[1];
      DRT::Node* node0 = BTSolDiscret().gNode(node0_gid);
      DRT::Node* node1 = BTSolDiscret().gNode(node1_gid);

      // get coordinates of edge nodes
      std::vector<double> x_n0(3);
      std::vector<double> x_n1(3);
      for (int j=0;j<3;++j)
      {
        x_n0[j] = node0->X()[j];
        x_n1[j] = node1->X()[j];
      }

      // compute distance vector and length
      // (APPROXIMATION FOR HIGHER-ORDER ELEMENTS !!!)
      std::vector<double> dist(3);
      for (int j=0;j<3;++j) dist[j] = x_n0[j] - x_n1[j];
      elelength = sqrt(dist[0]*dist[0]+dist[1]*dist[1]+dist[2]*dist[2]);
    }
    else if (BEAMCONTACT::RigidsphereElement(*thisele))
      continue;         // elelength does not apply for rigid spheres, radius is already considered in MaxEleRadius(), so simply do nothing here
    else
      dserror("The function GetMaxEleLength is only defined for beam elements and rigid sphere elements!");

    // if current length is larger than maximum length -> update
    if (elelength > maxelelength) maxelelength = elelength;
  }

  return;
}

/*-------------------------------------------------------------------- -*
 |  Update contact forces at the end of time step             popp 04/10|
 *----------------------------------------------------------------------*/
void CONTACT::Beam3cmanager::Update(const Epetra_Vector& disrow, const int& timestep,
                                    const int& newtonstep)
{
  // store values of fc_ into fcold_ (generalized alpha)
  fcold_->Update(1.0,*fc_,0.0);

  double disold_L2=0.0;
  dis_old_->Norm2(&disold_L2);

  // calculate (dis_old_-dis_):
  dis_old_->Update(-1.0,*dis_,1.0);
  // calculate inf-norm of (dis_old_-dis_)
  dis_old_->NormInf(&maxdeltadisp_);
  // invert the last step and get again dis_old_
  dis_old_->Update(1.0,*dis_,1.0);
  //update dis_old_ -> dis_
  dis_old_->Update(1.0,*dis_,0.0);

  if(maxdeltadisp_>totalmaxdeltadisp_ and disold_L2>1.0e-12)//don't consider the first time step where disold_L2=0!
    totalmaxdeltadisp_=maxdeltadisp_;

  //If the original gap function definition is applied, the displacement per time is not allowed
  //to be larger than the smalles beam cross section radius occurring in the discretization!
  bool newgapfunction=DRT::INPUT::IntegralValue<int>(BeamContactParameters(),"BEAMS_NEWGAP");
  if(!newgapfunction)
  {
    double maxdeltadisscalefac=sbeamcontact_.get<double>("BEAMS_MAXDELTADISSCALEFAC",1.0);
    //TODO: shall we allow for larger displacements in the first time step in general?
    if(maxdeltadisp_>maxdeltadisscalefac*mineleradius_ and timestep!=1)
    {
      //std::cout << "Minimal element radius: " << mineleradius_ << std::endl;
      //std::cout << "Maximal displacement per time step: " << maxdeltadisp_ << std::endl;
      //dserror("Displacement increment per time step larger than smallest beam element radius, "
      //        "but newgapfunction_ flag is not set. Choose smaller time step!");
    }
  }
  //std::cout << "mineleradius_: " << mineleradius_ << std::endl;
  //std::cout << "totalmaxdeltadisp_: " << totalmaxdeltadisp_ << std::endl;

  // create gmsh output for nice visualization
  // (endoftimestep flag set to TRUE)
#ifdef GMSHTIMESTEPS
  GmshOutput(disrow, timestep, newtonstep, true);
#endif

  // First, we check some restrictions concerning the new gap function definition
  for (int i=0;i<(int)pairs_.size();++i)
  {
    // only relevant if current pair is active
    if (pairs_[i]->GetContactFlag() == true)
    {
      if (pairs_[i]->GetNewGapStatus() == true)
      {
        //Necessary when using the new gap function definition (ngf_=true) for very slender beams in order to avoid crossing:
        //For very low penalty parameters and very slender beams it may happen that in the converged configuration the remaining penetration is
        //larger than the sum of the beam radii (R1+R2), i.e. the beam centerlines remain crossed even in the converged configuration. In this
        //case the sign of the normal vector normal_ has to be inverted at the end of the time step, since this quantity is stored in normal_old_ afterwards.
        //Otherwise the contact force would be applied in the wrong direction and the beams could cross!
        pairs_[i]->InvertNormal();
        Teuchos::ParameterList ioparams = DRT::Problem::Instance()->IOParams();
        if(!pdiscret_.Comm().MyPID() && ioparams.get<int>("STDOUTEVRY",0))
          std::cout << "      Warning: Penetration to large, choose higher penalty parameter!" << std::endl;
      }

      if (pairs_[i]->GetShiftStatus() == true)
      {
      //In case the contact points of two beams are identical (i.e. r1=r2), the nodal coordinates of one beam are shifted by a small
      //predefined value in order to enable evaluation of the contact pair. This makes the Newton scheme more robust. However,
      //in the converged configuration we want to have the real nodal positions for all contact pairs!!!
      dserror("Contact pair with identical contact points (i.e. r1=r2) not possible in the converged configuration!");
      }
    }
  }

  // set normal_old_=normal_ for all contact pairs at the end of the time step
  UpdateAllPairs();

  // print some data to screen
  ConsoleOutput();
  //store pairs_ in oldpairs_ to be available in next time step
  //this is needed for the new gapfunction defintion and also for the output at the end of an time step
  oldpairs_.clear();
  oldpairs_.resize(0);
  oldpairs_ = pairs_;

  oldcontactpairmap_.clear();
  oldcontactpairmap_ = contactpairmap_;

  oldbtsolpairs_.clear();
  oldbtsolpairs_.resize(0);
  oldbtsolpairs_ = btsolpairs_;

  oldbtsolpairmap_.clear();
  oldbtsolpairmap_ = btsolpairmap_;

  // clear potential contact pairs
  contactpairmap_.clear();
  btsphpairmap_.clear();
  btsolpairmap_.clear();

  pairs_.clear();
  pairs_.resize(0);

  btsolpairs_.clear();
  btsolpairs_.resize(0);

  btsolmtgroups_.clear();
  btsolmtgroups_.resize(0);

  btsphpairs_.clear();
  btsphpairs_.resize(0);

  btbpotpairs_.clear();
  btbpotpairs_.resize(0);

  btsphpotpairs_.clear();
  btsphpotpairs_.resize(0);

  return;
}

/*----------------------------------------------------------------------*
 |  Write Gmsh data for current state                          popp 04/10|
 *----------------------------------------------------------------------*/
void CONTACT::Beam3cmanager::GmshOutput(const Epetra_Vector& disrow, const int& timestep,
                                        const int& newtonstep, bool endoftimestep)
{

  //**********************************************************************
  // create filename for ASCII-file for output
  //**********************************************************************

  #if defined(CONTACTPAIRSPECIFICOUTPUT) or defined(DISTINGUISHCONTACTCOLOR)
    //Since all the gmsh-output is written by proc0 (this is necessary in order to keep the correct
    //order of the nodes and intermediate points when visualizing Bezier curves with Blender) the pairs_
    //vector would have to be communicated before writing the output
    if(btsoldiscret_->Comm().NumProc()>1)
      dserror("Contact pair specific gmsh output is not implemented in parallel so far.");
  #endif

//  if (btsol_)
//    dserror("GmshOutput not implemented for beam-to-solid contact so far!");

  #ifdef OUTPUTEVERY
    if(fabs(timen_ -(outputcounter_+1)*OUTPUTEVERY)<1.0e-8)
    {
      outputcounter_++;
    }
    else
    {
      return;
    }
  #endif

  // STEP 1: OUTPUT OF TIME STEP INDEX
  std::ostringstream filename;
  filename << "../o/gmsh_output/";

  #ifdef OUTPUTEVERY
    if (timestep<1000000)
      filename << "beams_t" << std::setw(6) << std::setfill('0') << outputcounter_;
    else /*(timestep>=1000000)*/
      dserror("ERROR: Gmsh output implemented for max 999.999 time steps");
  #else
    if (timestep<1000000)
      filename << "beams_t" << std::setw(6) << std::setfill('0') << timestep;
    else /*(timestep>=1000000)*/
      dserror("ERROR: Gmsh output implemented for max 999.999 time steps");
  #endif

  // STEPS 2/3: OUTPUT OF UZAWA AND NEWTON STEP INDEX
  // (for the end of time step output, omit this)
  int uzawastep = 99;
  if (!endoftimestep)
  {
    // check if Uzawa index is needed or not
    INPAR::BEAMCONTACT::Strategy strategy =
    DRT::INPUT::IntegralValue<INPAR::BEAMCONTACT::Strategy>(BeamContactParameters(),"BEAMS_STRATEGY");

    if (strategy==INPAR::BEAMCONTACT::bstr_uzawa)
    {
      uzawastep = uzawaiter_;
      if (uzawastep<10)
        filename << "_u0";
      else if (uzawastep<100)
        filename << "_u";
      else /*(uzawastep>=100)*/
        dserror("ERROR: Gmsh output implemented for max 99 Uzawa steps");
      filename << uzawastep;
    }

    // Newton index is always needed, of course
    if (newtonstep<10)
      filename << "_n00";
    else if (newtonstep<100)
      filename << "_n0";
    else if (newtonstep<1000)
      filename << "_n";
    else /*(newtonstep>=1000*/
      dserror("ERROR: Gmsh output implemented for max 999 Newton steps");
    filename << newtonstep;
  }

  // finish filename
  filename<<".pos";

  //**********************************************************************
  // gsmh output beam elements as prisms
  //**********************************************************************

  // approximation of the circular cross section with n prisms
  int n=N_CIRCUMFERENTIAL;        // CHOOSE N HERE
  if (n<3) n=3;    // minimum is 3, else no volume is defined.
  int n_axial=N_Axial; //number of devisions of element in axial direction

  //extract fully overlapping displacement vector on contact discretization from
  //displacement vector in row map format on problem discretization
  Epetra_Vector disccol(*BTSolDiscret().DofColMap(),true);
  ShiftDisMap(disrow, disccol);

    // do output to file in c-style
    FILE* fp = NULL;

  //The whole gmsh output is done by proc 0!
  if (btsoldiscret_->Comm().MyPID()==0)
  {
    // open file to write output data into
    fp = fopen(filename.str().c_str(), "w");

    if (fp==NULL) dserror("can't write to specified gmsh output folder!");

    std::stringstream gmshfileheader;
      //write output to temporary std::stringstream;
//    gmshfileheader <<"General.BackgroundGradient = 0;\n";
//    gmshfileheader <<"View.LineType = 1;\n";
//    gmshfileheader <<"View.LineWidth = 1.4;\n";
//    gmshfileheader <<"View.PointType = 1;\n";
//    gmshfileheader <<"View.PointSize = 3;\n";
//    gmshfileheader <<"General.ColorScheme = 1;\n";
//    gmshfileheader <<"General.Color.Background = {255,255,255};\n";
//    gmshfileheader <<"General.Color.Foreground = {255,255,255};\n";
//    //gmshfileheader <<"General.Color.BackgroundGradient = {128,147,255};\n";
//    gmshfileheader <<"General.Color.Foreground = {85,85,85};\n";
//    gmshfileheader <<"General.Color.Text = {0,0,0};\n";
//    gmshfileheader <<"General.Color.Axes = {0,0,0};\n";
//    gmshfileheader <<"General.Color.SmallAxes = {0,0,0};\n";
//    gmshfileheader <<"General.Color.AmbientLight = {25,25,25};\n";
//    gmshfileheader <<"General.Color.DiffuseLight = {255,255,255};\n";
//    gmshfileheader <<"General.Color.SpecularLight = {255,255,255};\n";
//    gmshfileheader <<"View.ColormapAlpha = 1;\n";
//    gmshfileheader <<"View.ColormapAlphaPower = 0;\n";
//    gmshfileheader <<"View.ColormapBeta = 0;\n";
//    gmshfileheader <<"View.ColormapBias = 0;\n";
//    gmshfileheader <<"View.ColormapCurvature = 0;\n";
//    gmshfileheader <<"View.ColormapInvert = 0;\n";
//    gmshfileheader <<"View.ColormapNumber = 2;\n";
//    gmshfileheader <<"View.ColormapRotation = 0;\n";
//    gmshfileheader <<"View.ColormapSwap = 0;\n";
//    gmshfileheader <<"View.ColorTable = {Black,Yellow,Blue,Orange,Red,Cyan,Purple,Brown,Green};\n";

    //general stuff
    gmshfileheader <<"View.Axes = 0;\n";
    gmshfileheader <<"View.LineType = 1;\n";
    gmshfileheader <<"View.LineWidth = 1.5;\n";
    gmshfileheader <<"General.RotationCenterGravity=0;\n";

//    gmshfileheader <<"General.BackgroundGradient = 0;\n";
//    gmshfileheader <<"View.ShowScale = 0;\n";
//    gmshfileheader <<"General.SmallAxes = 0;\n";

    //output: patch test (y-plane view)
//    gmshfileheader <<"General.TrackballQuaternion0 = 4.329637285359677e-17;\n";
//    gmshfileheader <<"General.TrackballQuaternion1 = -0.7071067811865475;\n";
//    gmshfileheader <<"General.TrackballQuaternion2 = -0.7071067811865476;\n";
//    gmshfileheader <<"General.TrackballQuaternion3 = 4.329637285359678e-17;\n";


    //output: beam rotating on arc
    //gmshfileheader <<"View.MaxX = 1.00;\n";
    //gmshfileheader <<"View.MaxY = 0.9;\n";
    //gmshfileheader <<"View.MaxZ = 1.6;\n";
    //gmshfileheader <<"View.MinX = -1.0;\n";
    //gmshfileheader <<"View.MinY = -0.9;\n";
    //gmshfileheader <<"View.MinZ =  0.0;\n";
    //gmshfileheader <<"General.TranslationY = -1.0;\n";


    //output: two straight beams rotating
//    gmshfileheader <<"General.ScaleX = 1.518637403204836;\n";
//    gmshfileheader <<"General.ScaleY = 1.518637403204836;\n";
//    gmshfileheader <<"General.ScaleZ = 1.518637403204836;\n";
//    gmshfileheader <<"General.TrackballQuaternion0 = -0.6055997025224867;\n";
//    gmshfileheader <<"General.TrackballQuaternion1 = -0.7452668955647411;\n";
//    gmshfileheader <<"General.TrackballQuaternion2 = 0.04875456702342353;\n";
//    gmshfileheader <<"General.TrackballQuaternion3 = -0.274680262986493;\n";
//    gmshfileheader <<"General.TranslationX = 0.4366289952435975;\n";
//    gmshfileheader <<"General.TranslationY = -0.9107124357404688;\n";

    //output energy conservation
//    gmshfileheader <<"General.ScaleX = 2.106015884264052;\n";
//    gmshfileheader <<"General.ScaleY = 2.106015884264052;\n";
//    gmshfileheader <<"General.ScaleZ = 2.106015884264052;\n";
//    gmshfileheader <<"General.TrackballQuaternion0 = 0.5237725213884583;\n";
//    gmshfileheader <<"General.TrackballQuaternion1 = -0.3413733097709274;\n";
//    gmshfileheader <<"General.TrackballQuaternion2 = -0.350946206992174;\n";
//    gmshfileheader <<"General.TrackballQuaternion3 = 0.6971107293767856;\n";
//    gmshfileheader <<"General.TranslationX = -0.008622353150367048;\n";
//    gmshfileheader <<"General.TranslationY = -0.01368488040624668;\n";

    //output two beams twisting
//    gmshfileheader <<"General.ScaleX = 0.04642018729706283;\n";
//    gmshfileheader <<"General.ScaleY = 0.04642018729706283;\n";
//    gmshfileheader <<"General.ScaleZ = 0.04642018729706283;\n";
//    gmshfileheader <<"General.TrackballQuaternion0 = 0.108433898949273;\n";
//    gmshfileheader <<"General.TrackballQuaternion1 = -0.003546806330406991;\n";
//    gmshfileheader <<"General.TrackballQuaternion2 = -0.8535041687531515;\n";
//    gmshfileheader <<"General.TrackballQuaternion3 = -0.5096666985830173;\n";
//    gmshfileheader <<"General.TranslationX = 0.5325452732071692;\n";
//    gmshfileheader <<"General.TranslationY = 0.2738989770477544;\n";

//    // General view settings for color table: Set predefined color range (for contact status display with discrete colors)
//    gmshfileheader << "View.RangeType=2;" << std::endl;         // Custom value scale range type
//    gmshfileheader << "View.CustomMin=0;" << std::endl;         // Minimum value to be displayed
//    gmshfileheader << "View.CustomMax=1;" << std::endl;         // Maximum value to be displayed
//    gmshfileheader << "View.IntervalsType=3;" << std::endl;     // Use discrete interval type
//    gmshfileheader << "View.NbIso=3;" << std::endl;             // Use three color ranges
//    gmshfileheader << "View.ColorTable={{0,140,255},{142,255,120},{206,0,0}};" << std::endl; // Define three colors

    //gmshfileheader <<"View.PointSize = 5;\n";

    //no shadow
    //gmshfileheader <<"View.Light = 0;\n";

    //write content into file and close it
    fprintf(fp, gmshfileheader.str().c_str());
    fclose(fp);

    //fp = fopen(filename.str().c_str(), "w");
    fp = fopen(filename.str().c_str(), "a");

    std::stringstream gmshfilecontent;
    gmshfilecontent << "View \" Step T" << timestep;

    // information about Uzawa and Newton step
    if (!endoftimestep)
      gmshfilecontent << " U" << uzawastep << " N" << newtonstep;

    // finish step information
    gmshfilecontent << " \" {" << std::endl;

    // write content into file and close
    fprintf(fp,gmshfilecontent.str().c_str());
    fclose(fp);
  }

#ifdef OUTPUTALLPROCS
  int numoutputloops=btsoldiscret_->Comm().NumProc()
#else
  int numoutputloops=1;
#endif

  // loop over the participating processors each of which appends its part of the output to one output file
  for(int i=0;i<numoutputloops;i++)
  {
    if (btsoldiscret_->Comm().MyPID() == i)
    {
      fp = fopen(filename.str().c_str(), "a");
      std::stringstream gmshfilecontent;

      //loop over fully overlapping column element map of proc 0
      for (int i=0;i<FullElements()->NumMyElements();++i)
      {
        // get pointer onto current beam element
        DRT::Element* element = BTSolDiscret().lColElement(i);

        // write output in specific function
        const DRT::ElementType & eot = element->ElementType();

        // no output for solid elements here
        if (eot != DRT::ELEMENTS::Beam3ebType::Instance() and
            eot != DRT::ELEMENTS::Beam3Type::Instance() and
            eot != DRT::ELEMENTS::Beam3rType::Instance() and
            eot != DRT::ELEMENTS::Beam3kType::Instance()and
            eot != DRT::ELEMENTS::RigidsphereType::Instance())
          continue;

        //**********
        // BEAM3R
        //**********
        // standard procedure for Reissner beams or rigid spheres
        if ( eot == DRT::ELEMENTS::Beam3Type::Instance() or
             eot == DRT::ELEMENTS::Beam3rType::Instance() or
             eot == DRT::ELEMENTS::RigidsphereType::Instance())
        {
          //*******************************************************************
          // special output for BEAM3r with Hermite center line interpolation
          //*******************************************************************
          bool done = false;
          if (eot == DRT::ELEMENTS::Beam3rType::Instance())
          {
            // this cast is necessary in order to use the method ->Tref() and others
            const DRT::ELEMENTS::Beam3r* ele = dynamic_cast<const DRT::ELEMENTS::Beam3r*>(element);

            if (ele->HermiteCenterlineInterpolation())
            {
              // mark as already done
              done = true;

              // this cast is necessary in order to use the method ->Tref()
              const DRT::ELEMENTS::Beam3r* ele = dynamic_cast<const DRT::ELEMENTS::Beam3r*>(element);

              int nnodescl = 2;
              LINALG::SerialDenseMatrix nodalcoords(3,nnodescl);
              LINALG::SerialDenseMatrix nodaltangents(3,nnodescl);
              LINALG::SerialDenseMatrix coord(3,n_axial);

              // compute current nodal positions (center line only)
              for (int i=0;i<3;++i)
              {
                for (int j=0;j<nnodescl;++j)
                {
                  double referenceposition = ((element->Nodes())[j])->X()[i];
                  std::vector<int> dofnode = BTSolDiscret().Dof((element->Nodes())[j]);
                  double displacement = disccol[BTSolDiscret().DofColMap()->LID(dofnode[i])];
                  nodalcoords(i,j) =  referenceposition + displacement;
                  nodaltangents(i,j) =  ((ele->Tref())[j])(i) + disccol[BTSolDiscret().DofColMap()->LID(dofnode[6+i])];
                }
              }

              if (nnodescl==2)
              {
                std::vector<double> disp_totlag;
                disp_totlag.reserve(12);
                for (int n=0;n<2;++n)
                {
                  for (int d=0;d<3;++d)
                    disp_totlag.push_back(nodalcoords(d,n));
                  for (int d=0;d<3;++d)
                    disp_totlag.push_back(nodaltangents(d,n));
                }
                //Calculate axial positions within the element by using the Hermite interpolation
                for (int i=0;i<n_axial;i++)
                {
                  double xi=-1.0 + i*2.0/(n_axial -1); // parameter coordinate of position vector on beam centerline
                  LINALG::Matrix<3,1> r;
                  ele->GetPosAtXi(r,xi,disp_totlag); //position vector on beam centerline

                  for (int j=0;j<3;j++)
                    coord(j,i)=r(j);
                }
              }
              else
              {
                dserror("Only 2-noded center line interpolations (Hermite) possible so far!");
              }
              if(N_CIRCUMFERENTIAL!=0)
                GMSH_N_noded(n,n_axial,coord,element,gmshfilecontent);
              else
                GMSH_N_nodedLine(n,n_axial,coord,element,gmshfilecontent);
            }
          }
          //*******************************************************************
          //*******************************************************************

          // standard case (no Hermite center line interpolation)
          if (!done)
          {
            // prepare storage for nodal coordinates
            int nnodes = element->NumNode();
            LINALG::SerialDenseMatrix coord(3,nnodes);

            // compute current nodal positions
            for (int id=0;id<3;++id)
            {
              for (int jd=0;jd<element->NumNode();++jd)
              {
                double referenceposition = ((element->Nodes())[jd])->X()[id];
                std::vector<int> dofnode = BTSolDiscret().Dof((element->Nodes())[jd]);
                double displacement = disccol[BTSolDiscret().DofColMap()->LID(dofnode[id])];
                coord(id,jd) =  referenceposition + displacement;
              }
            }

            switch (element->NumNode())
            {
              // rigid sphere element (1 node)
              case 1:
              {
                GMSH_sphere(coord,element,gmshfilecontent);
                break;
              }
              // 2-noded beam element (linear interpolation)
              case 2:
              {
                GMSH_2_noded(n,coord,element,gmshfilecontent);
                break;
              }
              // 3-noded beam element (quadratic nterpolation)
              case 3:
              {
                GMSH_3_noded(n,coord,element,gmshfilecontent);
                break;
              }
              // 4-noded beam element (quadratic interpolation)
              case 4:
              {
                GMSH_4_noded(n,coord,element,gmshfilecontent);
                break;
              }
              // 4- or 5-noded beam element (higher-order interpolation)
              default:
              {
                dserror("Gmsh output for %i noded element not yet implemented!", element->NumNode());
                break;
              }
            }
          }
        }

        //**********
        // BEAM3EB
        //**********
        // initially straight torsion-free Kirchhoff beams need a special treatment
        else if (eot == DRT::ELEMENTS::Beam3ebType::Instance())
        {
          // this cast is necessary in order to use the method ->Tref()
          const DRT::ELEMENTS::Beam3eb* ele = dynamic_cast<const DRT::ELEMENTS::Beam3eb*>(element);
          // prepare storage for nodal coordinates
          int nnodes = element->NumNode();
          LINALG::SerialDenseMatrix nodalcoords(3,nnodes);
          LINALG::SerialDenseMatrix nodaltangents(3,nnodes);
          LINALG::SerialDenseMatrix coord(3,n_axial);

          // compute current nodal positions
          for (int i=0;i<3;++i)
          {
            for (int j=0;j<element->NumNode();++j)
            {
              double referenceposition = ((element->Nodes())[j])->X()[i];
              std::vector<int> dofnode = BTSolDiscret().Dof((element->Nodes())[j]);
              double displacement = disccol[BTSolDiscret().DofColMap()->LID(dofnode[i])];
              nodalcoords(i,j) =  referenceposition + displacement;
              nodaltangents(i,j) =  ((ele->Tref())[j])(i) + disccol[BTSolDiscret().DofColMap()->LID(dofnode[3+i])];
            }
          }

          if (nnodes ==2)
          {
            LINALG::Matrix<12,1> disp_totlag(true);
            for (int i=0;i<3;i++)
            {
              disp_totlag(i)=nodalcoords(i,0);
              disp_totlag(i+6)=nodalcoords(i,1);
              disp_totlag(i+3)=nodaltangents(i,0);
              disp_totlag(i+9)=nodaltangents(i,1);
            }
            //Calculate axial positions within the element by using the Hermite interpolation of Kirchhoff beams
            for (int i=0;i<n_axial;i++)
            {
              double xi=-1.0 + i*2.0/(n_axial -1); // parameter coordinate of position vector on beam centerline
              LINALG::Matrix<3,1> r = ele->GetPos(xi, disp_totlag); //position vector on beam centerline

              for (int j=0;j<3;j++)
                coord(j,i)=r(j);
            }
          }
          else
          {
            dserror("Only 2-noded Kirchhoff elements possible so far!");
          }
          if(N_CIRCUMFERENTIAL!=0)
            GMSH_N_noded(n,n_axial,coord,element,gmshfilecontent);
          else
            GMSH_N_nodedLine(n,n_axial,coord,element,gmshfilecontent);
        }

        //**********
        // BEAM3K
        //**********
        // full strong/weak Kirchhoff beams need a special treatment
        else if (eot == DRT::ELEMENTS::Beam3kType::Instance())
        {
          // this cast is necessary in order to use the method ->Tref()
          const DRT::ELEMENTS::Beam3k* ele = dynamic_cast<const DRT::ELEMENTS::Beam3k*>(element);
          // prepare storage for nodal coordinates
          int nnodes = element->NumNode();
          LINALG::SerialDenseMatrix nodalcoords(3,nnodes);
          LINALG::SerialDenseMatrix nodaltangents(3,nnodes);
          LINALG::SerialDenseMatrix coord(3,n_axial);

          // compute current nodal positions
          for (int i=0;i<3;++i)
          {
            for (int j=0;j<element->NumNode();++j)
            {
              double referenceposition = ((element->Nodes())[j])->X()[i];
              std::vector<int> dofnode = BTSolDiscret().Dof((element->Nodes())[j]);
              double displacement = disccol[BTSolDiscret().DofColMap()->LID(dofnode[i])];
              nodalcoords(i,j) =  referenceposition + displacement;
              if(j<2)
              {
                if (ele->RotVec())
                  nodaltangents(i,j) = ((ele->Theta0())[j])(i) + disccol[BTSolDiscret().DofColMap()->LID(dofnode[3+i])];
                else
                  nodaltangents(i,j) =  ((ele->Tref())[j])(i) + disccol[BTSolDiscret().DofColMap()->LID(dofnode[3+i])];
              }
              else
              {
                nodaltangents(i,j) =0.0;
              }
            }
          }

          if (ele->RotVec())
          {
            // compute tangents from rotation vectors
            LINALG::Matrix<3,1> theta;
            LINALG::Matrix<3,3> R;
            for (int j=0;j<element->NumNode();++j)
            {
              theta(0) = nodaltangents(0,j);
              theta(1) = nodaltangents(1,j);
              theta(2) = nodaltangents(2,j);
              R.Clear();
              LARGEROTATIONS::angletotriad(theta,R);
              std::vector<int> dofnode = BTSolDiscret().Dof((element->Nodes())[j]);
              double lt = disccol[BTSolDiscret().DofColMap()->LID(dofnode[6])];
              nodaltangents(0,j) = (1.0+lt) * R(0,0);
              nodaltangents(1,j) = (1.0+lt) * R(1,0);
              nodaltangents(2,j) = (1.0+lt) * R(2,0);
            }
          }

          // remember that the beam3k is a 3-noded element only(!) due to the fact the
          // rotation interpolation needs a third node. With regard to the centerline that is
          // of interest here, the beam3k element is a 2-noded element with Hermite interpolation
          if (nnodes ==3)
          {
            LINALG::Matrix<12,1> disp_totlag(true);
            for (int i=0;i<3;i++)
            {
              disp_totlag(i)=nodalcoords(i,0);
              disp_totlag(i+6)=nodalcoords(i,1);
              disp_totlag(i+3)=nodaltangents(i,0);
              disp_totlag(i+9)=nodaltangents(i,1);
            }
            //Calculate axial positions within the element by using the Hermite interpolation of Kirchhoff beams
            for (int i=0;i<n_axial;i++)
            {
              double xi=-1.0 + i*2.0/(n_axial -1); // parameter coordinate of position vector on beam centerline
              LINALG::Matrix<3,1> r = ele->GetPos(xi, disp_totlag); //position vector on beam centerline

              for (int j=0;j<3;j++)
                coord(j,i)=r(j);
            }
          }
          else
          {
            dserror("Only 3-noded weak Kirchhoff elements possible so far!");
          }
          if(N_CIRCUMFERENTIAL!=0)
            GMSH_N_noded(n,n_axial,coord,element,gmshfilecontent);
          else
            GMSH_N_nodedLine(n,n_axial,coord,element,gmshfilecontent);
        }

        else
        {
          dserror("Your chosen type of beam element is not allowed for beam contact!");
        }
      }//loop over elements

      #ifdef CONTACTPAIRSPECIFICOUTPUT
      //loop over pairs vector in order to print normal vector
      for (int i=0;i<(int)pairs_.size();++i)
      {

        std::vector<LINALG::Matrix<3,1> > r1_vec = pairs_[i]->GetX1();
        std::vector<LINALG::Matrix<3,1> > r2_vec = pairs_[i]->GetX2();
        std::vector<double> contactforce = pairs_[i]->GetContactForce();

        int numcps = pairs_[i]->GetNumCps();

        for(int j=0;j<(int)r1_vec.size();j++)
        {
          LINALG::Matrix<3,1> normal(true);
          LINALG::Matrix<3,1> r1(true);
          LINALG::Matrix<3,1> r2(true);

          double fac=1.0;
          double color = 1.0;

          if(j<numcps)
          {
            fac=0.01;
            color = 0.9;
          }
          else
          {
            fac=0.1;
            color = 0.5;
          }

          for (int k=0;k<3;k++)
          {
            normal(k)=(r2_vec[j](k)-r1_vec[j](k));
          }
          normal.Scale(1.0/normal.Norm2());

          for (int k=0;k<3;k++)
          {
            r1(k)=r1_vec[j](k);
            r2(k)=r1_vec[j](k)+contactforce[j]*fac*normal(k);
          }

          gmshfilecontent << "SL("<< std::scientific;
          gmshfilecontent << r1(0) << "," << r1(1) << "," << r1(2) << ",";
          gmshfilecontent << r2(0) << "," << r2(1) << "," << r2(2);
          gmshfilecontent << "){" << std::scientific;
          gmshfilecontent << color << "," << color << "};"<< std::endl << std::endl;
        }
      }//loop over pairs
      #endif

      //**********
      // SOLIDS
      //**********
      //****************************begin: solid visualization
      #ifdef BTSOLGMSH
      // Loop over all column elements on this processor
      for (int i = 0; i < ProblemDiscret().NumMyColElements(); ++i)
      {
        // Get pointer onto current beam element
        DRT::Element* element = ProblemDiscret().lColElement(i);

        if (!BEAMCONTACT::BeamElement(*element))
        {
          GMSH_Solid(element, disrow, gmshfilecontent);
        }
      }
      // Loop over all column elements on this processor
      for (int i = 0; i < BTSolDiscret().NumMyColElements(); ++i)
      {
        // Get pointer onto current beam element
        DRT::Element* element = BTSolDiscret().lColElement(i);

        if (!BEAMCONTACT::BeamElement(*element))
        {
          GMSH_SolidSurfaceElementNumbers(element, disrow, gmshfilecontent);
        }
      }
      #endif

      #if defined(GMSHDEBUG) || defined(SAVEFORCE)
      // Get information about current problem settings for rendering 2D text in Gmsh and
      // for creating an unique filename if contact forces should be written in a text file

      // Get number of Gauss points used for one contact interval
      const int numgp = DRT::UTILS::IntegrationPoints1D(DRT::UTILS::GAUSSRULE).nquad;

      // Get moment of inertia
      //double Iyy1 = 0.0;
      //if ((int)btsolpairs_.size() > 0)
      //{
      //  DRT::Element* element1 = const_cast<DRT::Element*>(btsolpairs_[0]->Element1());
      //  Iyy1 = (dynamic_cast<DRT::ELEMENTS::Beam3eb*>(element1))->Iyy();
      //}

      // Get number of beam elements and surface elements on contact surface
      int numele1 = 0;
      int numele2 = 0;
      for (int i = 0; i < ColElements()->NumMyElements(); ++i)
      {
        DRT::Element* element = BTSolDiscret().lColElement(i);
        if(BEAMCONTACT::BeamElement(*element))
          numele1++;
        else
          numele2++;
      }

      // Get element type of surface elements
      std::string shape2 = "unknown";
      if ((int)btsolpairs_.size() > 0)
      {
        DRT::Element* element2 = const_cast<DRT::Element*>(btsolpairs_[0]->Element2());

        switch(element2->Shape())
        {
          case DRT::Element::tri3:
            shape2 = "tri3";
            break;
          case DRT::Element::tri6:
            shape2 = "tri6";
            break;
          case DRT::Element::quad4:
            shape2 = "quad4";
            break;
          case DRT::Element::quad8:
            shape2 = "quad8";
            break;
          case DRT::Element::quad9:
            shape2 = "quad9";
            break;
          default:
            shape2 = "unknown";
        }
      }
      #endif

      #ifdef GMSHDEBUG
      // Draw distance vector rD as line and value of gap as 3D text at each Gauss point

      // Initialize variables
      double r1[3];           // Position on beam element centerline
      double r2[3];           // Position on surface element
      double n2[3];           // Surface unit normal vector
      double gap;             // Gap function
      double rt[3];           // Position for text
      int total_numgp = 0;    // Total number of Gauss points
      int total_numgpc = 0;   // Number of Gauss points in contact

      // Loop over all beam to solid contact pairs
      for (int i_pair = 0; i_pair < (int)btsolpairs_.size(); i_pair++)
      {
        std::vector<Beam3tosolidcontactinterface::gmshDebugPoint> gmshDebugPoints = btsolpairs_[i_pair]->GetGmshDebugPoints();

        // Loop over all Gauss points
        for (int i = 0; i < (int)gmshDebugPoints.size(); i++)
        {
          // Declaring variable for color of lines
          double color = 0.5; // Inactive contact

          Beam3tosolidcontactinterface::gmshDebugPoint gmshDebugPoint = gmshDebugPoints[i];

          // Draw only information at Gauss points
          if (gmshDebugPoint.type == 0) // 0: Gauss point
          {
            total_numgp++;
            gap = gmshDebugPoint.gap;

            // Check if current Gauss point has contact
            if (gap < 0)
            {
              // Increase number of Gauss points in contact and change color of line
              total_numgpc++;
              color = 1.0; // Active contact
            }

            // Get positions from gmshDebugPoint and calculate 3D text position
            for (int j = 0; j < 3; j++)
            {
              r1[j] = gmshDebugPoint.r1(j);
              r2[j] = gmshDebugPoint.x2(j);
              n2[j] = gmshDebugPoint.n2(j);
              rt[j] = r1[j] + 0.3 * n2[j];
            }

            // Draw line between beam element centerline and surface element
            gmshfilecontent << "SL(" << std::scientific << r1[0] << "," << r1[1] << "," << r1[2] << "," << r2[0] << "," << r2[1] << "," << r2[2] << ")";
            gmshfilecontent << "{" << std::scientific << color << "," << color << "};" << std::endl;

            // Draw line from Gauss point on beam element centerline in surface unit normal vector direction
            gmshfilecontent << "SL(" << std::scientific << r1[0] << "," << r1[1] << "," << r1[2] << "," << rt[0] << "," << rt[1] << "," << rt[2] << ")";
            gmshfilecontent << "{" << std::scientific << color << "," << color << "};" << std::endl;

            // Draw gap as 3D text
            gmshfilecontent << "T3(" << std::scientific << rt[0] << "," << rt[1] << "," << rt[2] << "," << 17 << ")";
            if (gap < 0)
              gmshfilecontent << "{\"" << " gap = " << gap << "\"};" << std::endl;
            else
              gmshfilecontent << "{\"" << " gap > 0" << "\"};" << std::endl;
          }
        }
      }

      // Plot 2D text with information about current problem settings
      int xpos = 15;
      int ypos = 5;
      const int yspace = 15;
      int textsize = 17;
      gmshfilecontent << "T2(" << std::scientific << xpos << "," << (ypos += yspace) << "," << textsize << ")";
      gmshfilecontent << "{\"" << "Number of Gauss points in contact intervall: " << numgp << "\"};" << std::endl;

      gmshfilecontent << "T2(" << std::scientific << xpos << "," << (ypos += yspace) << "," << textsize << ")";
      gmshfilecontent << "{\"" << "Number of beam elements:                     " << numele1 << "\"};" << std::endl;

      gmshfilecontent << "T2(" << std::scientific << xpos << "," << (ypos += yspace) << "," << textsize << ")";
      gmshfilecontent << "{\"" << "Number of solid elements:                    " << numele2 << "\"};" << std::endl;

      gmshfilecontent << "T2(" << std::scientific << xpos << "," << (ypos += yspace) << "," << textsize << ")";
      gmshfilecontent << "{\"" << "Solid element shape:                         " << shape2.c_str() << "\"};" << std::endl;

      //gmshfilecontent << "T2(" << std::scientific << 15 << "," << 80 << "," << 17 << ")";
      //gmshfilecontent << "{\"" << "Beam moment of inertia:                      " << Iyy1 << "\"};" << std::endl;

      ypos += 5;

      gmshfilecontent << "T2(" << std::scientific << xpos << "," << (ypos += yspace) << "," << textsize << ")";
      gmshfilecontent << "{\"" << "Total number of Gauss points                 " << total_numgp << "\"};" << std::endl;

      gmshfilecontent << "T2(" << std::scientific << xpos << "," << (ypos += yspace) << "," << textsize << ")";
      gmshfilecontent << "{\"" << "Total number of Gauss points in contact:     " << total_numgpc << "\"};" << std::endl;
      #endif
      //***************************end: solid visualization*************************************************************

      // write content into file and close
      fprintf(fp,gmshfilecontent.str().c_str());
      fclose(fp);
    }
    Comm().Barrier();
  }

  Comm().Barrier();
  //Add a white and a black point -> this is necessary in order to get the full color range
  if (btsoldiscret_->Comm().MyPID()==0)
  {
    fp = fopen(filename.str().c_str(), "a");
    std::stringstream gmshfilecontent;

    // finish data section of this view by closing curley brackets (somehow needed to get color)
    gmshfilecontent <<"SP(0.0,0.0,0.0){0.0,0.0};"<<std::endl;
    gmshfilecontent <<"SP(0.0,0.0,0.0){1.0,1.0};"<<std::endl;
    gmshfilecontent << "};" << std::endl;

    // write content into file and close
    fprintf(fp,gmshfilecontent.str().c_str());
    fclose(fp);
  }
  Comm().Barrier();

  //*************************begin: gmsh output of contact forces for solids************************************
  #ifdef GMSHFORCE
  // Draw contact forces acting on solid element using an extra Gmsh output file
  FILE* fp_fc2 = NULL;

  // Create filename
  std::ostringstream filename_fc2;
  filename_fc2 << "o/gmsh_output/";
  if (timestep < 1000000)
    filename_fc2 << "fc2_t" << std::setw(6) << std::setfill('0') << timestep;
  else
    dserror("ERROR: Gmsh output implemented for max 999.999 time steps");

  if (!endoftimestep)
  {
     // Newton index is always needed, of course
     if (newtonstep < 10)
       filename_fc2 << "_n0";
     else if (newtonstep < 100)
       filename_fc2 << "_n";
     else
       dserror("ERROR: Gmsh output implemented for max 99 Newton steps");
     filename_fc2 << newtonstep;
  }

  // Finish filename
  filename_fc2 << ".pos";

  fp_fc2 = fopen(filename_fc2.str().c_str(), "w");

  std::stringstream fc2filecontent;

  // General settings for 3D vector visualization
  fc2filecontent << "General.ArrowHeadRadius=0.08;" << std::endl; // Arrow head radius
  fc2filecontent << "General.ArrowStemLength=0.6;" << std::endl;  // Arrow stem length
  fc2filecontent << "View.CenterGlyphs=2;" << std::endl;          // Locate with position of arrow head
  fc2filecontent << "View.ArrowSizeMax=100;" << std::endl;        // Arrow maximal size

  // Draw solid surface elements with different colors for contact status
  fc2filecontent << "View \"StepT" << timestep << "_Solid\"{" << std::endl;

  // Loop over all column elements on this processor
  for (int i = 0; i < ProblemDiscret().NumMyColElements(); ++i)
  {
    // Get pointer onto current beam element
    DRT::Element* element = ProblemDiscret().lColElement(i);

    if (!BEAMCONTACT::BeamElement(*element))
    {
      GMSH_Solid(element, disrow, fc2filecontent);
    }
  }
  fc2filecontent << "};" << std::endl;

  // General view settings for color table: Set predefined color range (for contact status)
  fc2filecontent << "View.RangeType=2;" << std::endl;                // Custom value scale range type
  fc2filecontent << "View.CustomMin=0;" << std::endl;                // Minimum value to be displayed
  fc2filecontent << "View.CustomMax=1;" << std::endl;                // Maximum value to be displayed
  fc2filecontent << "View.IntervalsType=3;" << std::endl;            // Use discrete interval type
  fc2filecontent << "View.NbIso=1;" << std::endl;                    // Use one color range
  fc2filecontent << "View.ColorTable={{180,180,180}};" << std::endl; // Define color

  // Draw contact forces acting on solid element
  fc2filecontent << "View \"StepT" << timestep << "_fc2\"{" << std::endl;

  // Loop over all beam to solid contact pairs
  for (int i_pair = 0; i_pair < (int)btsolpairs_.size(); i_pair++)
  {
    std::vector<Beam3tosolidcontactinterface::gmshDebugPoint> gmshDebugPoints = btsolpairs_[i_pair]->GetGmshDebugPoints();

    // Loop over all Gauss points
    for (int i = 0; i < (int)gmshDebugPoints.size(); i++)
    {
      Beam3tosolidcontactinterface::gmshDebugPoint gmshDebugPoint = gmshDebugPoints[i];

      // Draw only information at Gauss points
      if (gmshDebugPoint.type == 0) // 0: Gauss point
      {
        double gap = gmshDebugPoint.gap;
        if (gap < 0)
        {
          double r2[3];
          double n2[3];
          double fc2[3];
          for (int j = 0; j < 3; j++)
          {
            r2[j] = gmshDebugPoint.x2(j);
            n2[j] = gmshDebugPoint.n2(j);
            fc2[j] = -gmshDebugPoint.fp * n2[j];
          }

          // Draw 3D vector for contact force
          fc2filecontent << "VP(" << std::scientific
              << r2[0] << "," << r2[1] << "," << r2[2] << "){"
              << fc2[0] << "," << fc2[1] << "," << fc2[2] << "};" << std::endl;
        }
      }
    }
  }
  fc2filecontent << "};" << std::endl;

  // Write content into file and close
  fprintf(fp_fc2, fc2filecontent.str().c_str());
  fclose(fp_fc2);
  #endif
  //*************************end: gmsh output of contact forces for solids************************************

  return;
}

/*----------------------------------------------------------------------*
 |  Compute rotation matrix R                                cyron 01/09|
 *----------------------------------------------------------------------*/
void CONTACT::Beam3cmanager::TransformAngleToTriad(Epetra_SerialDenseVector& theta,
                                                   Epetra_SerialDenseMatrix& R)
{
  // compute spin matrix according to Crisfield Vol. 2, equation (16.8)
  Epetra_SerialDenseMatrix spin(3,3);
  ComputeSpin(spin,theta);

  // nompute norm of theta
  double theta_abs = theta.Norm2();

  // build an identity matrix
  Epetra_SerialDenseMatrix identity(3,3);
  for(int i=0;i<3;i++) identity(i,i) = 1.0;

  // square of spin matrix
  Epetra_SerialDenseMatrix spin2(3,3);
  for (int i=0;i<3;i++)
    for(int j=0;j<3;j++)
      for(int k=0;k<3;k++)
        spin2(i,k) += spin(i,j) * spin(j,k);


  // compute rotation matrix according to Crisfield Vol. 2, equation (16.22)
  for(int i=0;i<3;i++)
    for(int j=0;j<3;j++)
      R(i,j) = identity(i,j) + spin(i,j)*(sin(theta_abs))/theta_abs + (1-(cos(theta_abs)))/(pow(theta_abs,2)) * spin2(i,j);

  return;
}

/*----------------------------------------------------------------------*
 |  Compute spin (private)                                    cyron 01/09|
 *----------------------------------------------------------------------*/
void CONTACT::Beam3cmanager::ComputeSpin(Epetra_SerialDenseMatrix& spin,
                                         Epetra_SerialDenseVector& rotationangle)
{
  // initialization
  const double spinscale = 1.0;
  for (int i=0;i<rotationangle.Length();++i)
    rotationangle[i] *= spinscale;

  // initialize spin with zeros
  for(int i=0;i<3;i++)
    for(int j=0;j<3;j++)
      spin(i,j) = 0;

  // fill spin with values (see Crisfield Vol. 2, equation (16.8))
  spin(0,0) = 0;
  spin(0,1) = -rotationangle[2];
  spin(0,2) = rotationangle[1];
  spin(1,0) = rotationangle[2];
  spin(1,1) = 0;
  spin(1,2) = -rotationangle[0];
  spin(2,0) = -rotationangle[1];
  spin(2,1) = rotationangle[0];
  spin(2,2) = 0;

  return;
}

/*----------------------------------------------------------------------*
 |  intialize second, third, ... Uzawa step                   popp 12/11|
 *----------------------------------------------------------------------*/
void CONTACT::Beam3cmanager::InitializeUzawa(LINALG::SparseMatrix& stiffmatrix,
                                             Epetra_Vector& fres,
                                             const Epetra_Vector& disrow,
                                             Teuchos::ParameterList timeintparams,
                                             bool newsti)
{
  // since we will modify the graph of stiffmatrix by adding additional
  // contact tiffness entries, we have to uncomplete it
  stiffmatrix.UnComplete();

  // determine contact stiffness matrix scaling factor (new STI)
  // (this is due to the fact that in the new STI, we hand in the
  // already appropriately scaled effective stiffness matrix. Thus,
  // the additional contact stiffness terms must be equally scaled
  // here, as well. In the old STI, the complete scaling operation
  // is done after contact evaluation within the time integrator,
  // therefore no special scaling needs to be applied here.)
  double scalemat = 1.0;
  if (newsti) scalemat = 1.0 - alphaf_;

  if (DRT::INPUT::IntegralValue<INPAR::STR::MassLin>(sstructdynamic_,"MASSLIN") != INPAR::STR::ml_rotations)
  {
    // remove contact stiffness terms from stiffmatrix
    stiffmatrix.Add(*stiffc_, false, -scalemat, 1.0);
    // remove old contact force terms from fres
    fres.Update(-(1.0-alphaf_),*fc_,1.0);
    fres.Update(-alphaf_,*fcold_,1.0);
  }
  else
  {

    // remove contact stiffness terms from stiffmatrix
    stiffmatrix.Add(*stiffc_, false, -1.0, 1.0);
    // remove old contact force terms from fres
    fres.Update(-1.0,*fc_,1.0);
  }

  // now redo Evaluate()
  Teuchos::RCP<Epetra_Vector> nullvec = Teuchos::null;
  Evaluate(stiffmatrix,fres,disrow,timeintparams,newsti);

  return;
}

/*----------------------------------------------------------------------*
 |  Reset all Uzawa-based Lagrange multipliers                popp 04/10|
 *----------------------------------------------------------------------*/
void CONTACT::Beam3cmanager::ResetAlllmuzawa()
{
  // loop over all potential contact pairs
  for(int i=0;i<(int)pairs_.size();i++)
    pairs_[i]->Resetlmuzawa();

  for(int i=0;i<(int)btsolpairs_.size();i++)
    btsolpairs_[i]->Resetlmuzawa();

  for(int i=0;i<(int)btsphpairs_.size();i++)
    btsphpairs_[i]->Resetlmuzawa();
  return;
}

/*----------------------------------------------------------------------*
 |  Update contact constraint norm during uzawa iteration    meier 10/14|
 *----------------------------------------------------------------------*/
void CONTACT::Beam3cmanager::UpdateConstrNormUzawa()
{
  //Next we want to find out the maximal and minimal gap
  //We have to distinguish between maximal and minimal gap, since our penalty
  //force law can already become active for positive gaps.
  double maxgap = 0.0;
  double maxallgap = 0.0;
  double mingap = 0.0;
  double minallgap = 0.0;
  double maxrelgap = 0.0;
  double maxallrelgap = 0.0;
  double minrelgap = 0.0;
  double minallrelgap = 0.0;

  // loop over all pairs to find all gaps
  for (int i=0;i<(int)pairs_.size();++i)
  {
    // only relevant if current pair is active
    if (pairs_[i]->GetContactFlag() == true)
    {
      //get smaller radius of the two elements:
      double smallerradius=0.0;
      double radius1=BEAMCONTACT::CalcEleRadius(pairs_[i]->Element1());
      double radius2=BEAMCONTACT::CalcEleRadius(pairs_[i]->Element2());
      if(radius1<radius2)
        smallerradius=radius1;
      else
        smallerradius=radius2;

      std::vector<double> pairgaps = pairs_[i]->GetGap();
      for(int i=0;i<(int)pairgaps.size();i++)
      {
        double gap = pairgaps[i];

        double relgap = gap/smallerradius;

        if (gap>maxgap)
          maxgap = gap;

        if (gap<mingap)
          mingap = gap;

        if (relgap>maxrelgap)
          maxrelgap = relgap;

        if (relgap<minrelgap)
          minrelgap = relgap;
      }
    }
  }

  //So far, we only have the processor local extrema, but we want the extrema of the whole problem
  //As long as the beam contact discretization is full overlapping, all pairs are stored in all procs and
  //don't need this procedure. However, for future applications (i.e. when abstain from a fully overlapping
  //discretization) it might be useful.
  Comm().MaxAll(&maxgap,&maxallgap,1);
  Comm().MinAll(&mingap,&minallgap,1);
  Comm().MaxAll(&maxrelgap,&maxallrelgap,1);
  Comm().MinAll(&minrelgap,&minallrelgap,1);

  //Set class variable
#ifdef RELCONSTRTOL
  constrnorm_=fabs(minallrelgap);
#else
  constrnorm_=fabs(minallgap);
#endif

   // print results to screen
  Teuchos::ParameterList ioparams = DRT::Problem::Instance()->IOParams();
  if (Comm().MyPID()==0 && ioparams.get<int>("STDOUTEVRY",0))
  {
    std::cout << std::endl << "     ************************BTB*************************"<<std::endl;
    std::cout << "      Penalty parameter         = " << currentpp_ << std::endl;
    std::cout << "      Minimal current Gap       = " << minallgap << std::endl;
    std::cout << "      Minimal total unconv. Gap = " << mintotalsimunconvgap_ << std::endl;
    std::cout << "      Minimal current rel. Gap  = " << minallrelgap << std::endl;
    std::cout << "      Current Constraint Norm   = " << constrnorm_ << std::endl;
    std::cout << "      Maximal current Gap       = " << maxallgap << std::endl;
    std::cout << "      Maximal current rel. Gap  = " << maxallrelgap << std::endl;
    if ((int)btsolpairs_.size())
    {
      std::cout << std::endl << "      ************************BTS*************************"<<std::endl;
      std::cout << "      BTS-Penalty parameter     = " << btspp_ << std::endl;
      std::cout << "      Current Constraint Norm   = " << btsolconstrnorm_ << std::endl;
    }
    std::cout<<"      ****************************************************"<<std::endl;
  }

  return;
}

/*----------------------------------------------------------------------*
 |  Update contact constraint norm                            popp 04/10|
 *----------------------------------------------------------------------*/
void CONTACT::Beam3cmanager::UpdateConstrNorm()
{
  //Next we want to find out the maximal and minimal gap
  //We have to distinguish between maximal and minimal gap, since our penalty
  //force law can already become active for positive gaps.
  double maxgap = -1000.0;
  double maxgap_cp = -1000.0;
  double maxgap_gp = -1000.0;
  double maxgap_ep = -1000.0;
  double maxallgap = -1000.0;
  double maxallgap_cp = -1000.0;
  double maxallgap_gp = -1000.0;
  double maxallgap_ep = -1000.0;
  double mingap = 1000.0;
  double mingap_cp = 1000.0;
  double mingap_gp = 1000.0;
  double mingap_ep = 1000.0;
  double minallgap = 1000.0;
  double minallgap_cp = 1000.0;
  double minallgap_gp = 1000.0;
  double minallgap_ep = 1000.0;
  double maxrelgap = -1000.0;
  double maxallrelgap = -1000.0;
  double minrelgap = 1000.0;
  double minallrelgap = 1000.0;

  //Clear class variable
  totpenaltyenergy_=0.0;
  double proclocal_penaltyenergy=0.0;

  //Calculate contact work
  double deltapenaltywork=0.0;
  Epetra_Vector delta_disp(*dis_);
  delta_disp.Update(-1.0,*dis_old_,1.0);

  Epetra_Vector fc_alpha(*fc_);
  fc_alpha.Update(alphaf_,*fcold_,1.0-alphaf_);

  delta_disp.Dot(fc_alpha,&deltapenaltywork);

  deltapenaltywork=-deltapenaltywork;

  totpenaltywork_+=deltapenaltywork;

  // loop over all pairs to find all gaps
  for (int i=0;i<(int)pairs_.size();++i)
  {

    // make sure to evaluate each pair only once
    int firsteleid = (pairs_[i]->Element1())->Id();
    bool firstisinrowmap = RowElements()->MyGID(firsteleid);

    // only relevant if current pair is active
    if (pairs_[i]->GetContactFlag() == true and firstisinrowmap)
    {
      //Update penalty energy
      proclocal_penaltyenergy += pairs_[i]->GetEnergy();

      //get smaller radius of the two elements:
      double smallerradius=0.0;
      double radius1=BEAMCONTACT::CalcEleRadius(pairs_[i]->Element1());
      double radius2=BEAMCONTACT::CalcEleRadius(pairs_[i]->Element2());
      if(radius1<radius2)
        smallerradius=radius1;
      else
        smallerradius=radius2;

      std::vector<double> pairgaps = pairs_[i]->GetGap();
      int numcps=pairs_[i]->GetNumCps();
      int numgps=pairs_[i]->GetNumGps();
      int numeps=pairs_[i]->GetNumEps();

      dsassert(pairgaps.size() == numcps+numgps+numeps, "size mismatch! total"
          " number of gaps unequal sum of individual contact type gaps");

      for(int i=0;i<(int)pairgaps.size();i++)
      {
        double gap = pairgaps[i];

        if (i<numcps)
        {
          if (gap>maxgap_cp)
            maxgap_cp = gap;

          if (gap<mingap_cp)
            mingap_cp = gap;
        }
        else if (i<numcps+numgps)
        {
          if (gap>maxgap_gp)
            maxgap_gp = gap;

          if (gap<mingap_gp)
            mingap_gp = gap;
        }
        else if (i<numcps+numgps+numeps)
        {
          if (gap>maxgap_ep)
            maxgap_ep = gap;

          if (gap<mingap_ep)
            mingap_ep = gap;
        }

        double relgap = gap/smallerradius;

        if (gap>maxgap)
          maxgap = gap;

        if (gap<mingap)
          mingap = gap;

        if (relgap>maxrelgap)
          maxrelgap = relgap;

        if (relgap<minrelgap)
          minrelgap = relgap;
      }
    }
  }

  //So far, we only have the processor local extrema, but we want the extrema of the whole problem
  //As long as the beam contact discretization is full overlapping, all pairs are stored in all procs and
  //don't need this procedure. However, for future applications (i.e. when abstain from a fully overlapping
  //discretization) it might be useful.
  Comm().MaxAll(&maxgap,&maxallgap,1);
  Comm().MaxAll(&maxgap_cp,&maxallgap_cp,1);
  Comm().MaxAll(&maxgap_gp,&maxallgap_gp,1);
  Comm().MaxAll(&maxgap_ep,&maxallgap_ep,1);
  Comm().MinAll(&mingap,&minallgap,1);
  Comm().MinAll(&mingap_cp,&minallgap_cp,1);
  Comm().MinAll(&mingap_gp,&minallgap_gp,1);
  Comm().MinAll(&mingap_ep,&minallgap_ep,1);
  Comm().MaxAll(&maxrelgap,&maxallrelgap,1);
  Comm().MinAll(&minrelgap,&minallrelgap,1);

  Comm().SumAll(&proclocal_penaltyenergy,&totpenaltyenergy_,1);

  //So far, we have determined the extrema of the current time step. Now, we want to determine the
  //extrema of the total simulation:
  if (maxallgap>maxtotalsimgap_)
    maxtotalsimgap_ = maxallgap;
  if (maxallgap_cp>maxtotalsimgap_cp_)
    maxtotalsimgap_cp_ = maxallgap_cp;
  if (maxallgap_gp>maxtotalsimgap_gp_)
    maxtotalsimgap_gp_ = maxallgap_gp;
  if (maxallgap_ep>maxtotalsimgap_ep_)
    maxtotalsimgap_ep_ = maxallgap_ep;

  if (minallgap<mintotalsimgap_)
    mintotalsimgap_ = minallgap;
  if (minallgap_cp<mintotalsimgap_cp_)
    mintotalsimgap_cp_ = minallgap_cp;
  if (minallgap_gp<mintotalsimgap_gp_)
    mintotalsimgap_gp_ = minallgap_gp;
  if (minallgap_ep<mintotalsimgap_ep_)
    mintotalsimgap_ep_ = minallgap_ep;

  if (maxallrelgap>maxtotalsimrelgap_)
    maxtotalsimrelgap_ = maxallrelgap;

  if (minallrelgap<mintotalsimrelgap_)
    mintotalsimrelgap_ = minallrelgap;

  //Set class variable
#ifdef RELCONSTRTOL
  constrnorm_=fabs(minallrelgap);
#else
  constrnorm_=fabs(minallgap);
#endif

  //TODO: Adapt this, as soon as we have a concrete implementation of beam-to-solid contact element pairs
  btsolconstrnorm_ = 0.0;

   // print results to screen
  Teuchos::ParameterList ioparams = DRT::Problem::Instance()->IOParams();
  if (Comm().MyPID()==0 && ioparams.get<int>("STDOUTEVRY",0))
  {
    std::cout << std::endl << "      ***********************************BTB************************************"<<std::endl;
    std::cout << "      Penalty parameter                = " << currentpp_ << std::endl;
    std::cout << "      Minimal current Gap              = " << minallgap << std::endl;
    std::cout << "      Minimal total Point-to-Point Gap = " << mintotalsimgap_cp_ << std::endl;
    std::cout << "      Minimal total Line-to-Line Gap   = " << mintotalsimgap_gp_ << std::endl;
    std::cout << "      Minimal total Endpoint Gap       = " << mintotalsimgap_ep_ << std::endl;
    std::cout << "      Minimal total Gap                = " << mintotalsimgap_ << std::endl;
    std::cout << "      Minimal total unconv. Gap        = " << mintotalsimunconvgap_ << std::endl;
    std::cout << "      Minimal current rel. Gap         = " << minallrelgap << std::endl;
    std::cout << "      Current Constraint Norm          = " << constrnorm_ << std::endl;
    std::cout << "      Minimal total rel. Gap           = " << mintotalsimrelgap_ << std::endl;
    std::cout << "      Maximal current Gap              = " << maxallgap << std::endl;
    std::cout << "      Maximal total Point-to-Point Gap = " << maxtotalsimgap_cp_ << std::endl;
    std::cout << "      Maximal total Line-to-Line Gap   = " << maxtotalsimgap_gp_ << std::endl;
    std::cout << "      Maximal total Endpoint Gap       = " << maxtotalsimgap_ep_ << std::endl;
    std::cout << "      Maximal total Gap                = " << maxtotalsimgap_ << std::endl;
    std::cout << "      Maximal current rel. Gap         = " << maxallrelgap << std::endl;
    std::cout << "      Maximal total rel. Gap           = " << maxtotalsimrelgap_ << std::endl;
    if ((int)btsolpairs_.size())
    {
      std::cout << std::endl << "      ***********************************BTS************************************"<<std::endl;
      std::cout << "      BTS-Penalty parameter     = " << btspp_ << std::endl;
      std::cout << "      Current Constraint Norm   = " << btsolconstrnorm_ << std::endl;
    }
    std::cout<<"      **************************************************************************"<<std::endl;
  }

  return;
}

/*----------------------------------------------------------------------*
 |  Shift normal vector to "normal_old_"                     meier 10/10|
 *----------------------------------------------------------------------*/
void CONTACT::Beam3cmanager::UpdateAllPairs()
{
   // loop over all potential contact pairs
    for (int i=0;i<(int)pairs_.size();++i)
      pairs_[i]->UpdateClassVariablesStep();

    // loop over all potential beam to solid contact pairs
    for (int i = 0; i < (int)btsolpairs_.size(); ++i)
      btsolpairs_[i]->UpdateClassVariablesStep();
}

/*----------------------------------------------------------------------*
 |  Update all Uzawa-based Lagrange multipliers               popp 04/10|
 *----------------------------------------------------------------------*/
void CONTACT::Beam3cmanager::UpdateAlllmuzawa()
{
  // loop over all potential contact pairs
  for (int i=0;i<(int)pairs_.size();++i)
    pairs_[i]->Updatelmuzawa(currentpp_);

  for (int i=0;i<(int)btsolpairs_.size();++i)
    btsolpairs_[i]->Updatelmuzawa(btspp_);

  for (int i=0;i<(int)btsphpairs_.size();++i)
    btsphpairs_[i]->Updatelmuzawa(currentpp_);

  return;
}

/*----------------------------------------------------------------------*
 |  Reset penalty parameter                                   popp 04/10|
 *----------------------------------------------------------------------*/
void CONTACT::Beam3cmanager::ResetCurrentpp()
{
  // get initial value from input file and reset
  currentpp_ = BeamContactParameters().get<double>("BEAMS_BTBPENALTYPARAM");

  // get initial value from input file and reset
  btspp_ = BeamContactParameters().get<double>("BEAMS_BTSPENALTYPARAM");

  return;
}

/*----------------------------------------------------------------------*
 |  Reset Uzawa iteration index                               popp 04/10|
 *----------------------------------------------------------------------*/
void CONTACT::Beam3cmanager::ResetUzawaIter()
{
  // reset index to zero
  uzawaiter_ = 0;

  return;
}

/*----------------------------------------------------------------------*
 |  Update Uzawa iteration index                              popp 04/10|
 *----------------------------------------------------------------------*/
void CONTACT::Beam3cmanager::UpdateUzawaIter()
{
  // increase index by one
  uzawaiter_++;

  return;
}

/*----------------------------------------------------------------------*
 |  Update penalty parameter                                  popp 04/10|
 *----------------------------------------------------------------------*/
bool CONTACT::Beam3cmanager::IncreaseCurrentpp(const double& globnorm)
{
  // check convergence rate of Uzawa iteration
  // if too slow then empirically increase the penalty parameter
  bool update = false;
  if ( (globnorm >= 0.25 * constrnorm_) && (uzawaiter_ >= 2) )
  {
    currentpp_ = currentpp_ * 1.6;
    update = true;
  }
  return update;
}

/*----------------------------------------------------------------------*
 |  Print active set                                          popp 04/10|
 *----------------------------------------------------------------------*/
void CONTACT::Beam3cmanager::ConsoleOutput()
{
  Teuchos::ParameterList ioparams = DRT::Problem::Instance()->IOParams();
  if(ioparams.get<int>("STDOUTEVRY",0))
  {
    // begin output
    if (Comm().MyPID()==0)
    {
      std::cout << "\n      Active contact set------------------------------------------------------------\n";
      printf("      ID1            ID2              T xi       eta      angle   gap         force \n");
    }
    Comm().Barrier();

    double maxangle=0.0;
    double minangle=90.0;
    double mincpgap=1000.0;
    double mingpgap=1000.0;
    double minepgap=1000.0;
    double maxcpgap=-1000.0;
    double maxgpgap=-1000.0;
    double maxepgap=-1000.0;
    int numperpc=0;
    int numparc=0;
    int numepc=0;
    int numperpc_transitions=0;
    int numepc_boundarygausspoint=0;

    #ifdef PRINTGAPFILE
      double error=0.0;
      int gapsize=0;
      double refgap=-0.002;
     #endif

    double perpshiftangle1 = sbeamcontact_.get<double>("BEAMS_PERPSHIFTANGLE1");
    double perpshiftangle2 = sbeamcontact_.get<double>("BEAMS_PERPSHIFTANGLE2");

    #ifdef PRINTGAPSOVERLENGTHFILE
      step_++;
      std::ostringstream filename;
      filename << "gp_gapsandforces_1x49filaments_consistentaxialtension_40elefil_" << step_ <<".txt";

      // do output to file in c-style
      FILE* fp = NULL;
      std::stringstream filecontent;

      //if(firststep_)
      if(true)
      {
        // open file to write output data into
        fp = fopen(filename.str().c_str(), "w");
        filecontent << "ID1 ID2 xi eta gap force angles\n";
        firststep_=false;
      }
      else
      {
        // open file to write output data into
        fp = fopen(filename.str().c_str(), "a");
      }
    #endif

    // loop over all pairs
    for (int i=0;i<(int)pairs_.size();++i)
    {
      // check if this pair is active
      if (pairs_[i]->GetContactFlag())
      {
        // make sure to print each pair only once
        int firsteleid = (pairs_[i]->Element1())->Id();
        bool firstisinrowmap = RowElements()->MyGID(firsteleid);

        // abbreviations
        int id1 = (pairs_[i]->Element1())->Id();
        int id2 = (pairs_[i]->Element2())->Id();
        std::vector<double> gaps = pairs_[i]->GetGap();
        std::vector<int> types = pairs_[i]->GetContactType();
        std::vector<double> forces = pairs_[i]->GetContactForce();
        std::vector<double> angles = pairs_[i]->GetContactAngle();
        std::vector<std::pair<double,double> > closestpoints = pairs_[i]->GetClosestPoint();
        std::pair<int,int> numsegments = pairs_[i]->GetNumSegments();
        std::vector<std::pair<int,int> > segmentids = pairs_[i]->GetSegmentIds();

        // print some output (use printf-method for formatted output)
        if (firstisinrowmap)
        {
          for(int j=0;j<(int)gaps.size();j++)
          {

            #ifdef PRINTGAPFILE
              double gap=gaps[j];
              error+=fabs((gap-refgap)/refgap);
              gapsize++;
            #endif

            if(fabs(angles[j]/M_PI*180.0)>maxangle)
              maxangle=fabs(angles[j]/M_PI*180.0);

            if(fabs(angles[j]/M_PI*180.0)<minangle)
              minangle=fabs(angles[j]/M_PI*180.0);

            if(types[j]==0)
            {
              if(gaps[j]<mincpgap)
                mincpgap=gaps[j];
              if(gaps[j]>maxcpgap)
                maxcpgap=gaps[j];
              numperpc++;

              if(fabs(angles[j]/M_PI*180.0)<perpshiftangle2 and fabs(angles[j]/M_PI*180.0)>perpshiftangle1)
                numperpc_transitions++;
            }

            if(types[j]==1)
            {
              if(gaps[j]<mingpgap)
                mingpgap=gaps[j];
              if(gaps[j]>maxgpgap)
                maxgpgap=gaps[j];
              numparc++;
              if( ((id1+1)%10)==0 and fabs(closestpoints[j].first-0.99)<0.02)
                numepc_boundarygausspoint++;

              #ifdef PRINTGAPSOVERLENGTHFILE
                //if(id1>=10 and id1 <=19 and id2>=20 and id2<=29)
                if(id1>=320 and id1 <=359 and id2>=360 and id2<=399)
                filecontent << id1 << " "<< id2 << " " << closestpoints[j].first  << " "<<  closestpoints[j].second << " " << gaps[j] << " " << forces[j] << " " << angles[j]/M_PI*180.0 << "\n";
              #endif
            }

            if(types[j]==2)
            {
              if(gaps[j]<minepgap)
                minepgap=gaps[j];
              if(gaps[j]>maxepgap)
                maxepgap=gaps[j];
              numepc++;
            }

            printf("      %-6d (%2d/%-2d) %-6d (%2d/%-2d)   %-1d  %-6.2f %-6.2f %-7.2f %-11.8e %-11.2e \n",id1,segmentids[j].first+1,numsegments.first,id2,segmentids[j].second+1,numsegments.second,types[j],closestpoints[j].first,closestpoints[j].second,angles[j]/M_PI*180.0,gaps[j],forces[j]);
            fflush(stdout);
          }
        }
      }
    }

#ifdef PRINTGAPSOVERLENGTHFILE
  //write content into file and close it
  fprintf(fp, filecontent.str().c_str());
  fclose(fp);
#endif

    //Calculate sum over all procs
    double sumpro_maxangle=0.0;
    double sumpro_minangle=0.0;
    double sumpro_mincpgap=0.0;
    double sumpro_mingpgap=0.0;
    double sumpro_minepgap=0.0;
    double sumpro_maxcpgap=0.0;
    double sumpro_maxgpgap=0.0;
    double sumpro_maxepgap=0.0;
    int sumpro_numperpc=0;
    int sumpro_numparc=0;
    int sumpro_numepc=0;
    int sumpro_numperpc_transitions=0;
    int sumpro_numepc_boundarygausspoint=0;

    Comm().MaxAll(&maxangle,&sumpro_maxangle,1);
    Comm().MinAll(&minangle,&sumpro_minangle,1);
    Comm().MinAll(&mincpgap,&sumpro_mincpgap,1);
    Comm().MinAll(&mingpgap,&sumpro_mingpgap,1);
    Comm().MinAll(&minepgap,&sumpro_minepgap,1);
    Comm().MaxAll(&maxcpgap,&sumpro_maxcpgap,1);
    Comm().MaxAll(&maxgpgap,&sumpro_maxgpgap,1);
    Comm().MaxAll(&maxepgap,&sumpro_maxepgap,1);
    Comm().SumAll(&numperpc,&sumpro_numperpc,1);
    Comm().SumAll(&numparc,&sumpro_numparc,1);
    Comm().SumAll(&numepc,&sumpro_numepc,1);
    Comm().SumAll(&numperpc_transitions,&sumpro_numperpc_transitions,1);
    Comm().SumAll(&numepc_boundarygausspoint,&sumpro_numepc_boundarygausspoint,1);

    #ifdef PRINTNUMCONTACTSFILE
    if(Comm().MyPID()==0)
    {
      //TODO
      std::ostringstream filename;
      filename << "activecontacts_minmaxgapangle_statmech_37filaments_noendpoints.txt";

      // do output to file in c-style
      FILE* fp = NULL;
      std::stringstream filecontent;

      if(firststep_)
      {
        // open file to write output data into
        fp = fopen(filename.str().c_str(), "w");
        filecontent << "perpcontacts parcontacts epcontacts pertransitions minangle maxangle minperpgap minpargap minepgap global_kappa_max\n";

        firststep_=false;
      }
      else
      {
        // open file to write output data into
        fp = fopen(filename.str().c_str(), "a");
      }

      filecontent << sumpro_numperpc << " " << sumpro_numparc << " " << sumpro_numepc << " " << sumpro_numperpc_transitions << " ";
      filecontent << sumpro_minangle << " " << sumpro_maxangle << " " << sumpro_mincpgap << " " << sumpro_mingpgap << " " << sumpro_minepgap<< " " << global_kappa_max_ << "\n";

      //write content into file and close it
      fprintf(fp, filecontent.str().c_str());
      fclose(fp);
    }
    #endif

    #ifdef PRINTGAPFILE
      error=error/gapsize;
      std::cout << "error: " << error << std::endl;

      std::ostringstream filename;
      filename << "gaps.txt";

      // do output to file in c-style
      FILE* fp = NULL;
      std::stringstream filecontent;

      if(firststep_)
      {
        // open file to write output data into
        fp = fopen(filename.str().c_str(), "w");
        filecontent << "gaperrors\n";
        firststep_=false;
      }
      else
      {
        // open file to write output data into
        fp = fopen(filename.str().c_str(), "a");
      }

      filecontent << error << "\n";

      //write content into file and close it
      fprintf(fp, filecontent.str().c_str());
      fclose(fp);
    #endif

    if (Comm().MyPID()==0)
    {
      std::cout << "Number of Point-to-Point Contact Pairs: " << sumpro_numperpc << std::endl;
      std::cout << "Number of Point Contacts in Transition Range: " << sumpro_numperpc_transitions << std::endl;
      std::cout << "Number of Line-to-Line Contact Pairs: " << sumpro_numparc << std::endl;
      std::cout << "Number of Boundary Gauss points at left end of bundle: " << sumpro_numepc_boundarygausspoint << std::endl;
      std::cout << "Number of Endpoint Contact Pairs: " << sumpro_numepc << std::endl;

      std::cout << "Minimal contact angle: " << sumpro_minangle << std::endl;
      std::cout << "Maximal contact angle: " << sumpro_maxangle << std::endl;

      std::cout << "Minimal Point-to-Point gap: " << sumpro_mincpgap << std::endl;
      std::cout << "Minimal Line-to-Line gap: " << sumpro_mingpgap << std::endl;
      std::cout << "Minimal Endpoint gap: " << sumpro_minepgap << std::endl;

      std::cout << "Maximal Point-to-Point gap: " << sumpro_maxcpgap << std::endl;
      std::cout << "Maximal Line-to-Line gap: " << sumpro_maxgpgap << std::endl;
      std::cout << "Maximal Endpoint gap: " << sumpro_maxepgap << std::endl;

      std::cout << " global_kappa_max_: " << global_kappa_max_ << std::endl;
      std::cout << " contactevaluationtime_: " << contactevaluationtime_ << std::endl;
    }

        // loop over all btsph pairs
    for (int i=0;i<(int)btsphpairs_.size();++i)
    {
      // check if this pair is active
      if (btsphpairs_[i]->GetContactFlag())
      {
        // get coordinates of contact point of each element
        Epetra_SerialDenseVector x1 = btsphpairs_[i]->GetX1();
        Epetra_SerialDenseVector x2 = btsphpairs_[i]->GetX2();

        // make sure to print each pair only once
        // (TODO: this is not yet enough...)
        int firsteleid = (btsphpairs_[i]->Element1())->Id();
        bool firstisinrowmap = RowElements()->MyGID(firsteleid);

        // abbreviations
        int id1 = (btsphpairs_[i]->Element1())->Id();
        int id2 = (btsphpairs_[i]->Element2())->Id();
        double gap = btsphpairs_[i]->GetGap();
        double lm = btsphpairs_[i]->Getlmuzawa() - currentpp_ * btsphpairs_[i]->GetGap();

        // print some output (use printf-method for formatted output)
        if (firstisinrowmap)
        {
          printf("ACTIVE BTSPH PAIR: %d & %d \t gap: %e \t lm: %e \n",id1,id2,gap,lm);
          fflush(stdout);
        }
      }
    }

    // end output
    Comm().Barrier();
    if (Comm().MyPID()==0) std::cout << std::endl;
  }
  return;
}

/*----------------------------------------------------------------------*
 |  Write reaction forces and moments into a csv-file         popp 04/10|
 *----------------------------------------------------------------------*/
void CONTACT::Beam3cmanager::Reactions(const Epetra_Vector& fint,
                                       const Epetra_Vector& dirichtoggle,
                                       const int& timestep)
{

  dserror("Reaction Forces are not implemented up to now!");

  // we need to address the nodes / dofs via the beam contact
  // discretization, because only this is exported to full overlap
  Epetra_Vector fintbc(fint);
  fintbc.ReplaceMap(*BTSolDiscret().DofRowMap());
  Epetra_Vector dirichtogglebc(dirichtoggle);
  dirichtogglebc.ReplaceMap(*BTSolDiscret().DofRowMap());

  // compute bearing reactions from fint via dirichtoggle
  // Note: dirichtoggle is 1 for DOFs with DBC and 0 elsewise
  Epetra_Vector fbearing(*BTSolDiscret().DofRowMap());
  fbearing.Multiply(1.0,dirichtogglebc,fintbc,0.0);

  // std::stringstream for filename
  std::ostringstream filename;
  filename << "o/gmsh_output/reaction_forces_moments.csv";

  // do output to file in c-style
  FILE* fp = NULL;

  // open file to write output data into
  if (timestep==1) fp = fopen(filename.str().c_str(), "w");
  else             fp = fopen(filename.str().c_str(), "a");

  // std::stringstream for file content
  std::ostringstream CSVcontent;
  CSVcontent << std::endl << timestep <<",";

  // only implemented for one single node
  int i=0;  // CHOOSE YOUR NODE ID HERE!!!
  const DRT::Node* thisnode = BTSolDiscret().gNode(i);
  const std::vector<int> DofGIDs = BTSolDiscret().Dof(thisnode);
  CSVcontent << i;

  // write reaction forces and moments
  for (int j=0;j<6;++j) CSVcontent << "," << fbearing[i*6+j];

  // write content into file and close
  fprintf(fp,CSVcontent.str().c_str());
  fclose(fp);

  return;
}

/*----------------------------------------------------------------------*
 |  Compute gmsh output for 2-noded elements (private)        popp 04/10|
 *----------------------------------------------------------------------*/
void CONTACT::Beam3cmanager::GMSH_2_noded(const int& n,
                                          const Epetra_SerialDenseMatrix& coord,
                                          const DRT::Element* thisele,
                                          std::stringstream& gmshfilecontent)
{
  // some local variables
  Epetra_SerialDenseMatrix prism(3,6);
  Epetra_SerialDenseVector axis(3);
  Epetra_SerialDenseVector radiusvec1(3);
  Epetra_SerialDenseVector radiusvec2(3);
  Epetra_SerialDenseVector auxvec(3);
  Epetra_SerialDenseVector theta(3);
  Epetra_SerialDenseMatrix R(3,3);

  double eleradius = BEAMCONTACT::CalcEleRadius(thisele);

  // declaring variable for color of elements
  double color = 1.0;
#ifdef DISTINGUISHCONTACTCOLOR
  for (int i=0;i<(int)pairs_.size();++i)
  {
    // abbreviations
    int id1 = (pairs_[i]->Element1())->Id();
    int id2 = (pairs_[i]->Element2())->Id();
    bool active = pairs_[i]->GetContactFlag();

    // if element is memeber of an active contact pair, choose different color
    if ( (thisele->Id()==id1 || thisele->Id()==id2) && active) color = 0.5;
  }
  for (int i=0;i<(int)btsphpairs_.size();++i)
  {
    // abbreviations
    int id1 = (btsphpairs_[i]->Element1())->Id();
    int id2 = (btsphpairs_[i]->Element2())->Id();
    bool active = btsphpairs_[i]->GetContactFlag();

    // if element is memeber of an active contact pair, choose different color
    if ( (thisele->Id()==id1 || thisele->Id()==id2) && active) color = 0.5;
  }
#endif

  // compute three dimensional angle theta
  for (int j=0;j<theta.Length();++j)
    axis[j] = coord(j,1) - coord(j,0);
  double norm_axis = axis.Norm2();
  for (int j=0;j<axis.Length();++j)
    theta[j] = axis[j] / norm_axis * 2 * M_PI / n;

  // Compute rotation matirx R
  TransformAngleToTriad(theta,R);

  // Now the first prism will be computed via two radiusvectors, that point from each of
  // the nodes to two points on the beam surface. Further prisms will be computed via a
  // for-loop, where the second node of the previous prism is used as the first node of the
  // next prism, whereas the central points (=nodes) stay  identic for each prism. The
  // second node will be computed by a rotation matrix and a radiusvector.

  // compute radius vector for first surface node of first prims
  for (int j=0;j<3;++j) auxvec[j] = coord(j,0) + norm_axis;

  // radiusvector for point on surface
  radiusvec1[0] = auxvec[1]*axis[2] - auxvec[2]*axis[1];
  radiusvec1[1] = auxvec[2]*axis[0] - auxvec[0]*axis[2];
  radiusvec1[2] = auxvec[0]*axis[1] - auxvec[1]*axis[0];

  // initialize all prism points to nodes
  for (int j=0;j<3;++j)
  {
    prism(j,0) = coord(j,0);
     prism(j,1) = coord(j,0);
    prism(j,2) = coord(j,0);
    prism(j,3) = coord(j,1);
    prism(j,4) = coord(j,1);
    prism(j,5) = coord(j,1);
  }

  // get first point on surface for node1 and node2
  for (int j=0;j<3;++j)
  {
    prism(j,1) += radiusvec1[j] / radiusvec1.Norm2() * eleradius;
    prism(j,4) += radiusvec1[j] / radiusvec1.Norm2() * eleradius;
  }

  // compute radiusvec2 by rotating radiusvec1 with rotation matrix R
  radiusvec2.Multiply('N','N',1,R,radiusvec1,0);

  // get second point on surface for node1 and node2
  for(int j=0;j<3;j++)
  {
    prism(j,2) += radiusvec2[j] / radiusvec2.Norm2() * eleradius;
    prism(j,5) += radiusvec2[j] / radiusvec2.Norm2() * eleradius;
  }

  // now first prism is built -> put coordinates into filecontent-stream
  // Syntax for gmsh input file
  // SI(x,y--,z,  x+.5,y,z,    x,y+.5,z,   x,y,z+.5, x+.5,y,z+.5, x,y+.5,z+.5){1,2,3,3,2,1};
  // SI( coordinates of the six corners ){colors}
  gmshfilecontent << "SI("<< std::scientific;
  gmshfilecontent << prism(0,0) << "," << prism(1,0) << "," << prism(2,0) << ",";
  gmshfilecontent << prism(0,1) << "," << prism(1,1) << "," << prism(2,1) << ",";
  gmshfilecontent << prism(0,2) << "," << prism(1,2) << "," << prism(2,2) << ",";
  gmshfilecontent << prism(0,3) << "," << prism(1,3) << "," << prism(2,3) << ",";
  gmshfilecontent << prism(0,4) << "," << prism(1,4) << "," << prism(2,4) << ",";
  gmshfilecontent << prism(0,5) << "," << prism(1,5) << "," << prism(2,5);
  gmshfilecontent << "){" << std::scientific;
  gmshfilecontent << color << "," << color << "," << color << "," << color << "," << color << "," << color << "};" << std::endl << std::endl;

  // now the other prisms will be computed
  for (int sector=0;sector<n-1;++sector)
  {
    // initialize for next prism
    // some nodes of last prims can be taken also for the new next prism
    for (int j=0;j<3;++j)
    {
      prism(j,1)=prism(j,2);
      prism(j,4)=prism(j,5);
      prism(j,2)=prism(j,0);
      prism(j,5)=prism(j,3);
    }

    // old radiusvec2 is now radiusvec1; radiusvec2 is set to zero
    for (int j=0;j<3;++j)
    {
      radiusvec1[j] = radiusvec2[j];
      radiusvec2[j] = 0.0;
    }

    // compute radiusvec2 by rotating radiusvec1 with rotation matrix R
    radiusvec2.Multiply('N','N',1,R,radiusvec1,0);

    // get second point on surface for node1 and node2
    for (int j=0;j<3;++j)
    {
       prism(j,2) += radiusvec2[j] / radiusvec2.Norm2() * eleradius;
       prism(j,5) += radiusvec2[j] / radiusvec2.Norm2() * eleradius;
    }

    // put coordinates into filecontent-stream
    // Syntax for gmsh input file
    // SI(x,y--,z,  x+.5,y,z,    x,y+.5,z,   x,y,z+.5, x+.5,y,z+.5, x,y+.5,z+.5){1,2,3,3,2,1};
    // SI( coordinates of the six corners ){colors}
    gmshfilecontent << "SI("<< std::scientific;
    gmshfilecontent << prism(0,0) << "," << prism(1,0) << "," << prism(2,0) << ",";
    gmshfilecontent << prism(0,1) << "," << prism(1,1) << "," << prism(2,1) << ",";
    gmshfilecontent << prism(0,2) << "," << prism(1,2) << "," << prism(2,2) << ",";
    gmshfilecontent << prism(0,3) << "," << prism(1,3) << "," << prism(2,3) << ",";
    gmshfilecontent << prism(0,4) << "," << prism(1,4) << "," << prism(2,4) << ",";
    gmshfilecontent << prism(0,5) << "," << prism(1,5) << "," << prism(2,5);
    gmshfilecontent << "){" << std::scientific;
    gmshfilecontent << color << "," << color << "," << color << "," << color << "," << color << "," << color << "};" << std::endl << std::endl;
  }

   // BEISPIEL: output nodal normals
  /*for (int k=0;k<thisele->NumNode();++k)
  {
   // get nodal coordinates
   DRT::Node* currnode = thisele->Nodes()[k];
   double* coord = currnode->X();
   gmshfilecontent << "VP(" << std::scientific << coord[0] << "," << coord[1] << "," << coord[2] << ")";
   gmshfilecontent << "{" << std::scientific << nn[0] << "," << nn[1] << "," << nn[2] << "};" << std::endl;
  }*/

   return;
}

/*----------------------------------------------------------------------*
 |  Compute gmsh output for 3-noded elements (private)        popp 04/10|
 *----------------------------------------------------------------------*/
void CONTACT::Beam3cmanager::GMSH_3_noded(const int& n,
                                          const Epetra_SerialDenseMatrix& allcoord,
                                          const DRT::Element* thisele,
                                          std::stringstream& gmshfilecontent)
{
  // some local variables
  Epetra_SerialDenseMatrix prism(3,6);
  Epetra_SerialDenseVector axis(3);
  Epetra_SerialDenseVector radiusvec1(3);
  Epetra_SerialDenseVector radiusvec2(3);
  Epetra_SerialDenseVector auxvec(3);
  Epetra_SerialDenseVector theta(3);
  Epetra_SerialDenseMatrix R(3,3);
  Epetra_SerialDenseMatrix coord(3,2);

  double eleradius = BEAMCONTACT::CalcEleRadius(thisele);

  // declaring variable for color of elements
  double color = 1.0;
#ifdef DISTINGUISHCONTACTCOLOR
  for (int i=0;i<(int)pairs_.size();++i)
  {
    // abbreviations
    int id1 = (pairs_[i]->Element1())->Id();
    int id2 = (pairs_[i]->Element2())->Id();
    bool active = pairs_[i]->GetContactFlag();

    // if element is memeber of an active contact pair, choose different color
    if ( (thisele->Id()==id1 || thisele->Id()==id2) && active) color = 0.0;
  }
  for (int i=0;i<(int)btsphpairs_.size();++i)
  {
    // abbreviations
    int id1 = (btsphpairs_[i]->Element1())->Id();
    int id2 = (btsphpairs_[i]->Element2())->Id();
    bool active = btsphpairs_[i]->GetContactFlag();

    // if element is memeber of an active contact pair, choose different color
    if ( (thisele->Id()==id1 || thisele->Id()==id2) && active) color = 0.875;
  }
#endif

  // computation of coordinates starts here
  // first, the prisms between node 1 and 3 will be computed.
  // afterwards the prisms between nodes 3 and 2 will follow
  for (int i=0;i<2;++i)
  {
    // prisms between node 1 and node 3
    if (i==0)
    {
      for (int j=0;j<3;++j)
      {
        coord(j,0) = allcoord(j,0);
        coord(j,1) = allcoord(j,2);
      }
    }

    // prisms between node 3 and node 2
    else if (i==1)
    {
      for (int j=0;j<3;++j)
      {
        coord(j,0) = allcoord(j,2);
        coord(j,1) = allcoord(j,1);
      }
    }

    // compute three dimensional angle theta
    for (int j=0;j<theta.Length();++j)
      axis[j] = coord(j,1) - coord(j,0);
    double norm_axis = axis.Norm2();
    for (int j=0;j<axis.Length();++j)
      theta[j] = axis[j] / norm_axis * 2 * M_PI / n;

    // Compute rotation matrix R
    TransformAngleToTriad(theta,R);

    // Now the first prism will be computed via two radiusvectors, that point from each of
    // the nodes to two points on the beam surface. Further prisms will be computed via a
    // for-loop, where the second node of the previous prism is used as the first node of the
    // next prism, whereas the central points (=nodes) stay  identic for each prism. The
    // second node will be computed by a rotation matrix and a radiusvector.

    // compute radius vector for first surface node of first prims
    for (int j=0;j<3;++j) auxvec[j] = coord(j,0) + norm_axis;

    // radiusvector for point on surface
    radiusvec1[0] = auxvec[1]*axis[2] - auxvec[2]*axis[1];
    radiusvec1[1] = auxvec[2]*axis[0] - auxvec[0]*axis[2];
    radiusvec1[2] = auxvec[0]*axis[1] - auxvec[1]*axis[0];

    // initialize all prism points to nodes
    for (int j=0;j<3;++j)
    {
      prism(j,0) = coord(j,0);
      prism(j,1) = coord(j,0);
      prism(j,2) = coord(j,0);
      prism(j,3) = coord(j,1);
      prism(j,4) = coord(j,1);
      prism(j,5) = coord(j,1);
    }

    // get first point on surface for node1 and node2
    for (int j=0;j<3;++j)
    {
      prism(j,1) += radiusvec1[j] / radiusvec1.Norm2() * eleradius;
      prism(j,4) += radiusvec1[j] / radiusvec1.Norm2() * eleradius;
    }

    // compute radiusvec2 by rotating radiusvec1 with rotation matrix R
    radiusvec2.Multiply('N','N',1,R,radiusvec1,0);

    // get second point on surface for node1 and node2
    for(int j=0;j<3;j++)
    {
      prism(j,2) += radiusvec2[j] / radiusvec2.Norm2() * eleradius;
      prism(j,5) += radiusvec2[j] / radiusvec2.Norm2() * eleradius;
    }

    // now first prism is built -> put coordinates into filecontent-stream
    // Syntax for gmsh input file
    // SI(x,y--,z,  x+.5,y,z,    x,y+.5,z,   x,y,z+.5, x+.5,y,z+.5, x,y+.5,z+.5){1,2,3,3,2,1};
    // SI( coordinates of the six corners ){colors}
    gmshfilecontent << "SI("<< std::scientific;
    gmshfilecontent << prism(0,0) << "," << prism(1,0) << "," << prism(2,0) << ",";
    gmshfilecontent << prism(0,1) << "," << prism(1,1) << "," << prism(2,1) << ",";
    gmshfilecontent << prism(0,2) << "," << prism(1,2) << "," << prism(2,2) << ",";
    gmshfilecontent << prism(0,3) << "," << prism(1,3) << "," << prism(2,3) << ",";
    gmshfilecontent << prism(0,4) << "," << prism(1,4) << "," << prism(2,4) << ",";
    gmshfilecontent << prism(0,5) << "," << prism(1,5) << "," << prism(2,5);
    gmshfilecontent << "){" << std::scientific;
    gmshfilecontent << color << "," << color << "," << color << "," << color << "," << color << "," << color << "};" << std::endl << std::endl;

    // now the other prisms will be computed
    for (int sector=0;sector<n-1;++sector)
    {
      // initialize for next prism
      // some nodes of last prims can be taken also for the new next prism
      for (int j=0;j<3;++j)
      {
        prism(j,1)=prism(j,2);
        prism(j,4)=prism(j,5);
        prism(j,2)=prism(j,0);
        prism(j,5)=prism(j,3);
      }

      // old radiusvec2 is now radiusvec1; radiusvec2 is set to zero
      for (int j=0;j<3;++j)
      {
        radiusvec1[j] = radiusvec2[j];
        radiusvec2[j] = 0.0;
      }

      // compute radiusvec2 by rotating radiusvec1 with rotation matrix R
      radiusvec2.Multiply('N','N',1,R,radiusvec1,0);

      // get second point on surface for node1 and node2
      for (int j=0;j<3;++j)
      {
        prism(j,2) += radiusvec2[j] / radiusvec2.Norm2() * eleradius;
        prism(j,5) += radiusvec2[j] / radiusvec2.Norm2() * eleradius;
      }

      // put coordinates into filecontent-stream
      // Syntax for gmsh input file
      // SI(x,y--,z,  x+.5,y,z,    x,y+.5,z,   x,y,z+.5, x+.5,y,z+.5, x,y+.5,z+.5){1,2,3,3,2,1};
      // SI( coordinates of the six corners ){colors}
      gmshfilecontent << "SI("<< std::scientific;
      gmshfilecontent << prism(0,0) << "," << prism(1,0) << "," << prism(2,0) << ",";
      gmshfilecontent << prism(0,1) << "," << prism(1,1) << "," << prism(2,1) << ",";
      gmshfilecontent << prism(0,2) << "," << prism(1,2) << "," << prism(2,2) << ",";
      gmshfilecontent << prism(0,3) << "," << prism(1,3) << "," << prism(2,3) << ",";
      gmshfilecontent << prism(0,4) << "," << prism(1,4) << "," << prism(2,4) << ",";
      gmshfilecontent << prism(0,5) << "," << prism(1,5) << "," << prism(2,5);
      gmshfilecontent << "){" << std::scientific;
      gmshfilecontent << color << "," << color << "," << color << "," << color << "," << color << "," << color << "};" << std::endl << std::endl;
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 |  Compute gmsh output for 3-noded elements (private)       meier 04/14|
 *----------------------------------------------------------------------*/
void CONTACT::Beam3cmanager::GMSH_4_noded(const int& n,
                                          const Epetra_SerialDenseMatrix& allcoord,
                                          const DRT::Element* thisele,
                                          std::stringstream& gmshfilecontent)
{
  // some local variables
  Epetra_SerialDenseMatrix prism(3,6);
  Epetra_SerialDenseVector axis(3);
  Epetra_SerialDenseVector radiusvec1(3);
  Epetra_SerialDenseVector radiusvec2(3);
  Epetra_SerialDenseVector auxvec(3);
  Epetra_SerialDenseVector theta(3);
  Epetra_SerialDenseMatrix R(3,3);
  Epetra_SerialDenseMatrix coord(3,2);

  double eleradius = 0;

  // get radius of element
  const DRT::ElementType & eot = thisele->ElementType();

    if ( eot == DRT::ELEMENTS::Beam3Type::Instance() )
      {
        const DRT::ELEMENTS::Beam3* thisbeam = static_cast<const DRT::ELEMENTS::Beam3*>(thisele);
        eleradius = MANIPULATERADIUSVIS*sqrt(sqrt(4 * (thisbeam->Izz()) / M_PI));
      }
    if ( eot == DRT::ELEMENTS::Beam3rType::Instance() )
      {
        const DRT::ELEMENTS::Beam3r* thisbeam = static_cast<const DRT::ELEMENTS::Beam3r*>(thisele);
        eleradius = MANIPULATERADIUSVIS*sqrt(sqrt(4 * (thisbeam->Izz()) / M_PI));
      }

  // declaring variable for color of elements
  double color = 1.0;
#ifdef DISTINGUISHCONTACTCOLOR
  for (int i=0;i<(int)pairs_.size();++i)
  {
    // abbreviations
    int id1 = (pairs_[i]->Element1())->Id();
    int id2 = (pairs_[i]->Element2())->Id();
    bool active = pairs_[i]->GetContactFlag();

    // if element is memeber of an active contact pair, choose different color
    if ( (thisele->Id()==id1 || thisele->Id()==id2) && active) color = 0.0;
  }
  for (int i=0;i<(int)btsphpairs_.size();++i)
  {
    // abbreviations
    int id1 = (btsphpairs_[i]->Element1())->Id();
    int id2 = (btsphpairs_[i]->Element2())->Id();
    bool active = btsphpairs_[i]->GetContactFlag();

    // if element is memeber of an active contact pair, choose different color
    if ( (thisele->Id()==id1 || thisele->Id()==id2) && active) color = 0.875;
  }
#endif

  // computation of coordinates starts here
  // first, the prisms between node 1 and 3 will be computed.
  // afterwards the prisms between nodes 3 and 2 will follow
  for (int i=0;i<3;++i)
  {
    // prisms between node 1 and node 3
    if (i==0)
    {
      for (int j=0;j<3;++j)
      {
        coord(j,0) = allcoord(j,0);
        coord(j,1) = allcoord(j,2);
      }
    }

    // prisms between node 3 and node 4
    else if (i==1)
    {
      for (int j=0;j<3;++j)
      {
        coord(j,0) = allcoord(j,2);
        coord(j,1) = allcoord(j,3);
      }
    }

    // prisms between node 4 and node 2
    else if (i==2)
    {
      for (int j=0;j<3;++j)
      {
        coord(j,0) = allcoord(j,3);
        coord(j,1) = allcoord(j,1);
      }
    }

    // compute three dimensional angle theta
    for (int j=0;j<theta.Length();++j)
      axis[j] = coord(j,1) - coord(j,0);
    double norm_axis = axis.Norm2();
    for (int j=0;j<axis.Length();++j)
      theta[j] = axis[j] / norm_axis * 2 * M_PI / n;

    // Compute rotation matrix R
    TransformAngleToTriad(theta,R);

    // Now the first prism will be computed via two radiusvectors, that point from each of
    // the nodes to two points on the beam surface. Further prisms will be computed via a
    // for-loop, where the second node of the previous prism is used as the first node of the
    // next prism, whereas the central points (=nodes) stay  identic for each prism. The
    // second node will be computed by a rotation matrix and a radiusvector.

    // compute radius vector for first surface node of first prims
    for (int j=0;j<3;++j) auxvec[j] = coord(j,0) + norm_axis;

    // radiusvector for point on surface
    radiusvec1[0] = auxvec[1]*axis[2] - auxvec[2]*axis[1];
    radiusvec1[1] = auxvec[2]*axis[0] - auxvec[0]*axis[2];
    radiusvec1[2] = auxvec[0]*axis[1] - auxvec[1]*axis[0];

    // initialize all prism points to nodes
    for (int j=0;j<3;++j)
    {
      prism(j,0) = coord(j,0);
      prism(j,1) = coord(j,0);
      prism(j,2) = coord(j,0);
      prism(j,3) = coord(j,1);
      prism(j,4) = coord(j,1);
      prism(j,5) = coord(j,1);
    }

    // get first point on surface for node1 and node2
    for (int j=0;j<3;++j)
    {
      prism(j,1) += radiusvec1[j] / radiusvec1.Norm2() * eleradius;
      prism(j,4) += radiusvec1[j] / radiusvec1.Norm2() * eleradius;
    }

    // compute radiusvec2 by rotating radiusvec1 with rotation matrix R
    radiusvec2.Multiply('N','N',1,R,radiusvec1,0);

    // get second point on surface for node1 and node2
    for(int j=0;j<3;j++)
    {
      prism(j,2) += radiusvec2[j] / radiusvec2.Norm2() * eleradius;
      prism(j,5) += radiusvec2[j] / radiusvec2.Norm2() * eleradius;
    }

    // now first prism is built -> put coordinates into filecontent-stream
    // Syntax for gmsh input file
    // SI(x,y--,z,  x+.5,y,z,    x,y+.5,z,   x,y,z+.5, x+.5,y,z+.5, x,y+.5,z+.5){1,2,3,3,2,1};
    // SI( coordinates of the six corners ){colors}
    gmshfilecontent << "SI("<< std::scientific;
    gmshfilecontent << prism(0,0) << "," << prism(1,0) << "," << prism(2,0) << ",";
    gmshfilecontent << prism(0,1) << "," << prism(1,1) << "," << prism(2,1) << ",";
    gmshfilecontent << prism(0,2) << "," << prism(1,2) << "," << prism(2,2) << ",";
    gmshfilecontent << prism(0,3) << "," << prism(1,3) << "," << prism(2,3) << ",";
    gmshfilecontent << prism(0,4) << "," << prism(1,4) << "," << prism(2,4) << ",";
    gmshfilecontent << prism(0,5) << "," << prism(1,5) << "," << prism(2,5);
    gmshfilecontent << "){" << std::scientific;
    gmshfilecontent << color << "," << color << "," << color << "," << color << "," << color << "," << color << "};" << std::endl << std::endl;

    // now the other prisms will be computed
    for (int sector=0;sector<n-1;++sector)
    {
      // initialize for next prism
      // some nodes of last prims can be taken also for the new next prism
      for (int j=0;j<3;++j)
      {
        prism(j,1)=prism(j,2);
        prism(j,4)=prism(j,5);
        prism(j,2)=prism(j,0);
        prism(j,5)=prism(j,3);
      }

      // old radiusvec2 is now radiusvec1; radiusvec2 is set to zero
      for (int j=0;j<3;++j)
      {
        radiusvec1[j] = radiusvec2[j];
        radiusvec2[j] = 0.0;
      }

      // compute radiusvec2 by rotating radiusvec1 with rotation matrix R
      radiusvec2.Multiply('N','N',1,R,radiusvec1,0);

      // get second point on surface for node1 and node2
      for (int j=0;j<3;++j)
      {
        prism(j,2) += radiusvec2[j] / radiusvec2.Norm2() * eleradius;
        prism(j,5) += radiusvec2[j] / radiusvec2.Norm2() * eleradius;
      }

      // put coordinates into filecontent-stream
      // Syntax for gmsh input file
      // SI(x,y--,z,  x+.5,y,z,    x,y+.5,z,   x,y,z+.5, x+.5,y,z+.5, x,y+.5,z+.5){1,2,3,3,2,1};
      // SI( coordinates of the six corners ){colors}
      gmshfilecontent << "SI("<< std::scientific;
      gmshfilecontent << prism(0,0) << "," << prism(1,0) << "," << prism(2,0) << ",";
      gmshfilecontent << prism(0,1) << "," << prism(1,1) << "," << prism(2,1) << ",";
      gmshfilecontent << prism(0,2) << "," << prism(1,2) << "," << prism(2,2) << ",";
      gmshfilecontent << prism(0,3) << "," << prism(1,3) << "," << prism(2,3) << ",";
      gmshfilecontent << prism(0,4) << "," << prism(1,4) << "," << prism(2,4) << ",";
      gmshfilecontent << prism(0,5) << "," << prism(1,5) << "," << prism(2,5);
      gmshfilecontent << "){" << std::scientific;
      gmshfilecontent << color << "," << color << "," << color << "," << color << "," << color << "," << color << "};" << std::endl << std::endl;
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 |  Compute gmsh output for N-noded elements (private)       meier 02/14|
 *----------------------------------------------------------------------*/
void CONTACT::Beam3cmanager::GMSH_N_noded(const int& n,
                                          int& n_axial,
                                          const Epetra_SerialDenseMatrix& allcoord,
                                          const DRT::Element* thisele,
                                          std::stringstream& gmshfilecontent)
{
  // some local variables
  Epetra_SerialDenseMatrix prism(3,6);
  Epetra_SerialDenseVector axis(3);
  Epetra_SerialDenseVector radiusvec1(3);
  Epetra_SerialDenseVector radiusvec2(3);
  Epetra_SerialDenseVector auxvec(3);
  Epetra_SerialDenseVector theta(3);
  Epetra_SerialDenseMatrix R(3,3);
  Epetra_SerialDenseMatrix coord(3,2);

  double eleradius = 0;

  // only implemented for Kirchhoff beams so far!
  const DRT::ElementType & eot = thisele->ElementType();
  if ( eot == DRT::ELEMENTS::Beam3ebType::Instance() )
  {
    const DRT::ELEMENTS::Beam3eb* thisbeam = static_cast<const DRT::ELEMENTS::Beam3eb*>(thisele);
    eleradius = MANIPULATERADIUSVIS*sqrt(sqrt(4 * (thisbeam->Izz()) / M_PI));
  }
  else if ( eot == DRT::ELEMENTS::Beam3kType::Instance() )
  {
    const DRT::ELEMENTS::Beam3k* thisbeam = static_cast<const DRT::ELEMENTS::Beam3k*>(thisele);
    eleradius = MANIPULATERADIUSVIS*sqrt(sqrt(4 * (thisbeam->Iyy()) / M_PI));
  }
  else if ( eot == DRT::ELEMENTS::Beam3rType::Instance() )
  {
    const DRT::ELEMENTS::Beam3r* thisbeam = static_cast<const DRT::ELEMENTS::Beam3r*>(thisele);
    eleradius = MANIPULATERADIUSVIS*sqrt(sqrt(4 * (thisbeam->Iyy()) / M_PI));
  }
  else
    dserror("ERROR: GSMH_N_noded output not yet implemented for this beam element type");


  // declaring variable for color of elements
  double color = 0.0;

  #ifdef DISTINGUISHCONTACTCOLOR
  for (int i=0;i<(int)pairs_.size();++i)
  {
    // abbreviations
    int id1 = (pairs_[i]->Element1())->Id();
    int id2 = (pairs_[i]->Element2())->Id();
    bool active = pairs_[i]->GetContactFlag();

    // if element is memeber of an active contact pair, choose different color
    if ( (thisele->Id()==id1 || thisele->Id()==id2) && active) color = 1.0;
  }
  for (int i=0;i<(int)btsphpairs_.size();++i)
  {
    // abbreviations
    int id1 = (btsphpairs_[i]->Element1())->Id();
    int id2 = (btsphpairs_[i]->Element2())->Id();
    bool active = btsphpairs_[i]->GetContactFlag();

    // if element is memeber of an active contact pair, choose different color
    if ( (thisele->Id()==id1 || thisele->Id()==id2) && active) color = 0.875;
  }
  for (int i_pair = 0; i_pair < (int)btsolpairs_.size(); ++i_pair)
  {
    // Id of beam element
    int id1 = (btsolpairs_[i_pair]->Element1())->Id();

    // If beam element is memeber of an active beam to solid contact pair, choose different color
    if (thisele->Id() == id1)
    {
      if (btsolpairs_[i_pair]->GetContactFlag())
      {
        // Beam to solid pair is in active contact (color red)
        color = 1.0;
        break;
      }
      else
      {
        // Beam to solid pair might come into active contact (color green)
        color = 0.5;
      }
    }
  }
  #endif

  // computation of coordinates starts here
  for (int i=0;i<n_axial-1;++i)
  {
    // prisms between node i and node i+1
    for (int j=0;j<3;++j)
    {
      coord(j,0) = allcoord(j,i);
      coord(j,1) = allcoord(j,i+1);
    }

    // separate elements by differently colored sub-segments
    double loccolor = color;
    if (i==0) loccolor = 0.25;

    // output of element IDs
    if (i==n_axial/2)
    {
      gmshfilecontent << "T3(" << std::scientific << coord(0,0) << "," << coord(1,0) << "," << coord(2,0) << "," << 17 << ")";
      gmshfilecontent << "{\"" << thisele->Id() << "\"};" << std::endl;
    }

    // compute three dimensional angle theta
    for (int j=0;j<3;++j)
      axis[j] = coord(j,1) - coord(j,0);

    double norm_axis = axis.Norm2();
    for (int j=0;j<3;++j)
      theta[j] = axis[j] / norm_axis * 2 * M_PI / n;

    // Compute rotation matrix R
    TransformAngleToTriad(theta,R);

    // Now the first prism will be computed via two radiusvectors, that point from each of
    // the nodes to two points on the beam surface. Further prisms will be computed via a
    // for-loop, where the second node of the previous prism is used as the first node of the
    // next prism, whereas the central points (=nodes) stay  identic for each prism. The
    // second node will be computed by a rotation matrix and a radiusvector.

    // compute radius vector for first surface node of first prims
    for (int j=0;j<3;++j) auxvec[j] = coord(j,0) + norm_axis;

    // radiusvector for point on surface
    radiusvec1[0] = auxvec[1]*axis[2] - auxvec[2]*axis[1];
    radiusvec1[1] = auxvec[2]*axis[0] - auxvec[0]*axis[2];
    radiusvec1[2] = auxvec[0]*axis[1] - auxvec[1]*axis[0];

    // initialize all prism points to nodes
    for (int j=0;j<3;++j)
    {
      prism(j,0) = coord(j,0);
      prism(j,1) = coord(j,0);
      prism(j,2) = coord(j,0);
      prism(j,3) = coord(j,1);
      prism(j,4) = coord(j,1);
      prism(j,5) = coord(j,1);
    }

    // get first point on surface for node1 and node2
    for (int j=0;j<3;++j)
    {
      prism(j,1) += radiusvec1[j] / radiusvec1.Norm2() * eleradius;
      prism(j,4) += radiusvec1[j] / radiusvec1.Norm2() * eleradius;
    }

    // compute radiusvec2 by rotating radiusvec1 with rotation matrix R
    radiusvec2.Multiply('N','N',1,R,radiusvec1,0);

    // get second point on surface for node1 and node2
    for(int j=0;j<3;j++)
    {
      prism(j,2) += radiusvec2[j] / radiusvec2.Norm2() * eleradius;
      prism(j,5) += radiusvec2[j] / radiusvec2.Norm2() * eleradius;
    }

    // now first prism is built -> put coordinates into filecontent-stream
    // Syntax for gmsh input file
    // SI(x,y--,z,  x+.5,y,z,    x,y+.5,z,   x,y,z+.5, x+.5,y,z+.5, x,y+.5,z+.5){1,2,3,3,2,1};
    // SI( coordinates of the six corners ){colors}
    gmshfilecontent << "SI("<< std::scientific;
    gmshfilecontent << prism(0,0) << "," << prism(1,0) << "," << prism(2,0) << ",";
    gmshfilecontent << prism(0,1) << "," << prism(1,1) << "," << prism(2,1) << ",";
    gmshfilecontent << prism(0,2) << "," << prism(1,2) << "," << prism(2,2) << ",";
    gmshfilecontent << prism(0,3) << "," << prism(1,3) << "," << prism(2,3) << ",";
    gmshfilecontent << prism(0,4) << "," << prism(1,4) << "," << prism(2,4) << ",";
    gmshfilecontent << prism(0,5) << "," << prism(1,5) << "," << prism(2,5);
    gmshfilecontent << "){" << std::scientific;
    gmshfilecontent << loccolor << "," << loccolor << "," << loccolor << "," << loccolor << "," << loccolor << "," << loccolor << "};" << std::endl << std::endl;

    // now the other prisms will be computed
    for (int sector=0;sector<n-1;++sector)
    {
      // initialize for next prism
      // some nodes of last prims can be taken also for the new next prism
      for (int j=0;j<3;++j)
      {
        prism(j,1)=prism(j,2);
        prism(j,4)=prism(j,5);
        prism(j,2)=prism(j,0);
        prism(j,5)=prism(j,3);
      }

      // old radiusvec2 is now radiusvec1; radiusvec2 is set to zero
      for (int j=0;j<3;++j)
      {
        radiusvec1[j] = radiusvec2[j];
        radiusvec2[j] = 0.0;
      }

      // compute radiusvec2 by rotating radiusvec1 with rotation matrix R
      radiusvec2.Multiply('N','N',1,R,radiusvec1,0);

      // get second point on surface for node1 and node2
      for (int j=0;j<3;++j)
      {
        prism(j,2) += radiusvec2[j] / radiusvec2.Norm2() * eleradius;
        prism(j,5) += radiusvec2[j] / radiusvec2.Norm2() * eleradius;
      }

      // put coordinates into filecontent-stream
      // Syntax for gmsh input file
      // SI(x,y--,z,  x+.5,y,z,    x,y+.5,z,   x,y,z+.5, x+.5,y,z+.5, x,y+.5,z+.5){1,2,3,3,2,1};
      // SI( coordinates of the six corners ){colors}
      gmshfilecontent << "SI("<< std::scientific;
      gmshfilecontent << prism(0,0) << "," << prism(1,0) << "," << prism(2,0) << ",";
      gmshfilecontent << prism(0,1) << "," << prism(1,1) << "," << prism(2,1) << ",";
      gmshfilecontent << prism(0,2) << "," << prism(1,2) << "," << prism(2,2) << ",";
      gmshfilecontent << prism(0,3) << "," << prism(1,3) << "," << prism(2,3) << ",";
      gmshfilecontent << prism(0,4) << "," << prism(1,4) << "," << prism(2,4) << ",";
      gmshfilecontent << prism(0,5) << "," << prism(1,5) << "," << prism(2,5);
      gmshfilecontent << "){" << std::scientific;
      gmshfilecontent << loccolor << "," << loccolor << "," << loccolor << "," << loccolor << "," << loccolor << "," << loccolor << "};" << std::endl << std::endl;
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 |  Compute gmsh output for N-noded elements (private)       meier 02/14|
 *----------------------------------------------------------------------*/
void CONTACT::Beam3cmanager::GMSH_N_nodedLine(const int& n,
                                              const int& n_axial,
                                              const Epetra_SerialDenseMatrix& allcoord,
                                              const DRT::Element* thisele,
                                              std::stringstream& gmshfilecontent)
{

  // declaring variable for color of elements
  double color = 1.0;

  #ifdef DISTINGUISHCONTACTCOLOR
    for (int i=0;i<(int)pairs_.size();++i)
    {
      // abbreviations
      int id1 = (pairs_[i]->Element1())->Id();
      int id2 = (pairs_[i]->Element2())->Id();
      bool active = pairs_[i]->GetContactFlag();

      // if element is memeber of an active contact pair, choose different color
      if ( (thisele->Id()==id1 || thisele->Id()==id2) && active) color = 0.0;
    }
    for (int i=0;i<(int)btsphpairs_.size();++i)
    {
      // abbreviations
      int id1 = (btsphpairs_[i]->Element1())->Id();
      int id2 = (btsphpairs_[i]->Element2())->Id();
      bool active = btsphpairs_[i]->GetContactFlag();

      // if element is memeber of an active contact pair, choose different color
      if ( (thisele->Id()==id1 || thisele->Id()==id2) && active) color = 0.875;
    }
  #endif

  for (int i=0;i<n_axial-1;i++)
  {
    // separate elements by black dots
    if (i==0)
      gmshfilecontent <<"SP(" << allcoord(0,i) <<  "," << allcoord(1,i) << "," << allcoord(2,i) << "){0.0,0.0};"<<std::endl;

    // output of element IDs
    if (i==n_axial/2)
    {
      gmshfilecontent << "T3(" << std::scientific << allcoord(0,i) << "," << allcoord(1,i) << "," << allcoord(2,i) << "," << 17 << ")";
      gmshfilecontent << "{\"" << thisele->Id() << "\"};" << std::endl;
    }

    gmshfilecontent << "SL("<< std::scientific;
    gmshfilecontent << allcoord(0,i) << "," << allcoord(1,i) << "," << allcoord(2,i) << ",";
    gmshfilecontent << allcoord(0,i+1) << "," << allcoord(1,i+1) << "," << allcoord(2,i+1);
    gmshfilecontent << "){" << std::scientific;
    gmshfilecontent << color << "," << color << "};" << std::endl;
  }

  return;
}

/*----------------------------------------------------------------------*
 |  Compute gmsh output for sphere elements (private)        grill 03/14|
 *----------------------------------------------------------------------*/
void CONTACT::Beam3cmanager::GMSH_sphere(const Epetra_SerialDenseMatrix& coord,
                                          const DRT::Element* thisele,
                                          std::stringstream& gmshfilecontent)
{
  double eleradius=0.0;
  double color=1.0;

  // get radius of element
  const DRT::ElementType & eot = thisele->ElementType();

  if ( eot == DRT::ELEMENTS::RigidsphereType::Instance() )
  {
    const DRT::ELEMENTS::Rigidsphere* thisparticle = static_cast<const DRT::ELEMENTS::Rigidsphere*>(thisele);
    eleradius = thisparticle->Radius();
  }
  else dserror("GMSH_sphere can only handle elements of Type Rigidsphere!");

  // loop over BTSPH pairs (if any), here we only need to check Element2 (Rigidsphere ele)
  for (int i=0;i<(int)btsphpairs_.size();++i)
  {
    // abbreviations
    int id2 = (btsphpairs_[i]->Element2())->Id();
    bool active = btsphpairs_[i]->GetContactFlag();

    // if element is member of an active contact pair, choose different color
    if ( thisele->Id()==id2 && active) color = 0.875;
  }

  // ********************** Visualization as a point ***********************************************
  // syntax for scalar point:  SP( coordinates x,y,z ){value at point (determines the color)}
//  gmshfilecontent <<"SP(" << std::scientific << coord(0,0) << "," << coord(1,0) << "," << coord(2,0) << "){" << std::scientific  << color << "};"<< std::endl << std::endl;
  // ***********************************************************************************************

  // ********************** Visualization as an icosphere ****************************************

  // for details see http://en.wikipedia.org/wiki/Icosahedron
  // and http://blog.andreaskahler.com/2009/06/creating-icosphere-mesh-in-code.html

  // sphere is visualized as icosphere
  // the basic icosahedron consists of 20 equilateral triangles (12 vertices)
  // further refinement by subdividing the triangles

  // list storing the (x,y,z) coordinates of all vertices
  std::vector< std::vector<double> > vertexlist(12, std::vector<double>(3,0));

  // list storing the indices of the three vertices that define a triangular face
  std::vector< std::vector<int> > facelist(20, std::vector<int>(3,0));

  double normfac = sqrt( 1.0 + 0.25* pow(1+sqrt(5),2) );
  double c = 0.5*(1.0+sqrt(5))/normfac*eleradius;
  double d = 1/normfac*eleradius;

  // compute the final coordinates of the initial 12 vertices
  vertexlist[0][0]+=-d; vertexlist[0][1]+=c; vertexlist[0][2]+=0;
  vertexlist[1][0]+=d; vertexlist[1][1]+=c; vertexlist[1][2]+=0;
  vertexlist[2][0]+=-d; vertexlist[2][1]+=-c; vertexlist[2][2]+=0;
  vertexlist[3][0]+=d; vertexlist[3][1]+=-c; vertexlist[3][2]+=0;

  vertexlist[4][0]+=0; vertexlist[4][1]+=-d; vertexlist[4][2]+=c;
  vertexlist[5][0]+=0; vertexlist[5][1]+=d; vertexlist[5][2]+=c;
  vertexlist[6][0]+=0; vertexlist[6][1]+=-d; vertexlist[6][2]+=-c;
  vertexlist[7][0]+=0; vertexlist[7][1]+=d; vertexlist[7][2]+=-c;

  vertexlist[8][0]+=c; vertexlist[8][1]+=0; vertexlist[8][2]+=-d;
  vertexlist[9][0]+=c; vertexlist[9][1]+=0; vertexlist[9][2]+=d;
  vertexlist[10][0]+=-c; vertexlist[10][1]+=0; vertexlist[10][2]+=-d;
  vertexlist[11][0]+=-c; vertexlist[11][1]+=0; vertexlist[11][2]+=d;


  // fill initial facelist
  facelist[0][0]=0;   facelist[0][1]=11;    facelist[0][2]=5;
  facelist[1][0]=0;   facelist[1][1]=5;     facelist[1][2]=1;
  facelist[2][0]=0;   facelist[2][1]=1;     facelist[2][2]=7;
  facelist[3][0]=0;   facelist[3][1]=7;     facelist[3][2]=10;
  facelist[4][0]=0;   facelist[4][1]=10;    facelist[4][2]=11;

  facelist[5][0]=1;   facelist[5][1]=5;     facelist[5][2]=9;
  facelist[6][0]=5;   facelist[6][1]=11;    facelist[6][2]=4;
  facelist[7][0]=11;  facelist[7][1]=10;    facelist[7][2]=2;
  facelist[8][0]=10;  facelist[8][1]=7;     facelist[8][2]=6;
  facelist[9][0]=7;   facelist[9][1]=1;     facelist[9][2]=8;

  facelist[10][0]=3;  facelist[10][1]=9;    facelist[10][2]=4;
  facelist[11][0]=3;  facelist[11][1]=4;    facelist[11][2]=2;
  facelist[12][0]=3;  facelist[12][1]=2;    facelist[12][2]=6;
  facelist[13][0]=3;  facelist[13][1]=6;    facelist[13][2]=8;
  facelist[14][0]=3;  facelist[14][1]=8;    facelist[14][2]=9;

  facelist[15][0]=4;  facelist[15][1]=9;    facelist[15][2]=5;
  facelist[16][0]=2;  facelist[16][1]=4;    facelist[16][2]=11;
  facelist[17][0]=6;  facelist[17][1]=2;    facelist[17][2]=10;
  facelist[18][0]=8;  facelist[18][1]=6;    facelist[18][2]=7;
  facelist[19][0]=9;  facelist[19][1]=8;    facelist[19][2]=1;

  // level of refinement, num_faces = 20 * 4^(ref_level)
  const int ref_level= 3;
  // refine the icosphere by calling GmshRefineIcosphere
  for (int p=0; p<ref_level; ++p)
  {
    GmshRefineIcosphere(vertexlist,facelist,eleradius);
  }

  const double centercoord[] = {coord(0,0), coord(1,0), coord(2,0)};
  for (unsigned int i=0; i<facelist.size(); ++i)
    PrintGmshTriangleToStream(gmshfilecontent,vertexlist,facelist[i][0],facelist[i][1],facelist[i][2],color,centercoord);

  // ********************* end: visualization as an icosphere ***************************************

  return;
}

/*----------------------------------------------------------------------*
 |  print Gmsh Triangle to stringstream (private)            grill 03/14|
 *----------------------------------------------------------------------*/
void CONTACT::Beam3cmanager::PrintGmshTriangleToStream(std::stringstream& gmshfilecontent,
                                                      const std::vector< std::vector< double > > &vertexlist,
                                                      int i, int j, int k, double color,
                                                      const double centercoord[])
{
  // "ST" is scalar triangle, followed by 3x coordinates (x,y,z) of vertices and color
  gmshfilecontent << "ST("<< std::scientific;
  gmshfilecontent << centercoord[0] + vertexlist[i][0] << "," << centercoord[1] + vertexlist[i][1] << "," << centercoord[2] + vertexlist[i][2] << ",";
  gmshfilecontent << centercoord[0] + vertexlist[j][0] << "," << centercoord[1] + vertexlist[j][1] << "," << centercoord[2] + vertexlist[j][2] << ",";
  gmshfilecontent << centercoord[0] + vertexlist[k][0] << "," << centercoord[1] + vertexlist[k][1] << "," << centercoord[2] + vertexlist[k][2];
  gmshfilecontent << "){" << std::scientific;
  gmshfilecontent << color << "," << color << "," << color  << "};" << std::endl << std::endl;

  return;
}

/*----------------------------------------------------------------------*
 |  Refine Icosphere (private)                                grill 03/14|
 *----------------------------------------------------------------------*/
void CONTACT::Beam3cmanager::GmshRefineIcosphere(std::vector< std::vector<double> > &vertexlist,
                                                  std::vector< std::vector<int> > &facelist,
                                                  double radius)
{
  int num_faces_old = facelist.size();
  std::vector<double> newvertex(3,0.0);
  double scalefac = 0.0;
  std::vector<int> newface(3,0);

  // subdivide each face into four new triangular faces:
  /*               /_\
   *              /_V_\
   */
  for (int i=0; i<num_faces_old; ++i)
  {
    int oldvertices[] = {facelist[i][0],facelist[i][1],facelist[i][2]};

    // compute, normalize and store new vertices in vertexlist
    // subdivide all three edges (all connections of three old vertices)
    for (int j=0; j<3; ++j)
    {
      for (int k=j+1; k<3; ++k)
      {
        newvertex[0]=0.5* (vertexlist[oldvertices[j]][0] + vertexlist[oldvertices[k]][0]);
        newvertex[1]=0.5* (vertexlist[oldvertices[j]][1] + vertexlist[oldvertices[k]][1]);
        newvertex[2]=0.5* (vertexlist[oldvertices[j]][2] + vertexlist[oldvertices[k]][2]);

        // scale new vertex to lie on sphere with given radius
        scalefac = radius / sqrt( pow(newvertex[0],2) + pow(newvertex[1],2) + pow(newvertex[2],2) );
        for (int q=0; q<3; ++q) newvertex[q] *= scalefac;

        vertexlist.push_back(newvertex);
      }
    }

    int len_vertexlist = (int) vertexlist.size();
    // add four new triangles to facelist
    newface[0]=oldvertices[0];  newface[1]=len_vertexlist-3; newface[2]=len_vertexlist-2;
    facelist.push_back(newface);
    newface[0]=oldvertices[1];  newface[1]=len_vertexlist-3; newface[2]=len_vertexlist-1;
    facelist.push_back(newface);
    newface[0]=oldvertices[2];  newface[1]=len_vertexlist-2; newface[2]=len_vertexlist-1;
    facelist.push_back(newface);
    newface[0]=len_vertexlist-3;  newface[1]=len_vertexlist-2; newface[2]=len_vertexlist-1;
    facelist.push_back(newface);
  }

  // erase the old faces
  facelist.erase(facelist.begin(), facelist.begin()+num_faces_old);

  return;
}

/*----------------------------------------------------------------------*
 | Compute Gmsh output for solid elements                               |
 *----------------------------------------------------------------------*/
void CONTACT::Beam3cmanager::GMSH_Solid(
    const DRT::Element* element,
    const Epetra_Vector& disrow,
    std::stringstream& gmshfilecontent)
{
  // Prepare storage for nodal coordinates
  int nnodes = element->NumNode();
  LINALG::SerialDenseMatrix coord(3, nnodes);

  // Compute current nodal positions
  for (int i_dim = 0; i_dim < 3; i_dim++)
  {
    for (int i_node = 0; i_node < element->NumNode(); i_node++)
    {
      double referenceposition = ((element->Nodes())[i_node])->X()[i_dim];
      std::vector<int> dofnode = ProblemDiscret().Dof((element->Nodes())[i_node]);
      double displacement = disrow[ProblemDiscret().DofColMap()->LID(dofnode[i_dim])];
      coord(i_dim,i_node) = referenceposition + displacement;
    }
  }

  // Declaring variable for standard color of solid elements
  double color = 0.2;

  // BACI surface pattern of hexahedral elements hex8, hex20, hex27
  const int surfNodes[6][9] = {{0,3,2,1,11,10,9,8,20}, {0,1,5,4,8,13,16,12,21},
      {1,2,6,5,9,14,17,13,22}, {2,3,7,6,10,15,18,14,23}, {0,4,7,3,12,19,15,11,24}, {4,5,6,7,16,17,18,19,25}};
  double surfColor[6] = {color, color, color, color, color, color};

  switch (element->Shape())
  {
    // 8-noded hexahedron
    case DRT::Element::hex8:
    {
      const int n_surfNodes = 4;
      GMSH_GetSurfColor(element, n_surfNodes, surfNodes, surfColor);

      // 4-noded quadrangle
      for (int i_surf = 0; i_surf < 6; i_surf++)
      {
        double coordsSQ[3][4];
        double colorsSQ[4];
        for (int i_nodeSQ = 0; i_nodeSQ < n_surfNodes; i_nodeSQ++)
        {
          for (int i_dim = 0; i_dim < 3; i_dim++)
          {
            coordsSQ[i_dim][i_nodeSQ] = coord(i_dim,surfNodes[i_surf][i_nodeSQ]);
            colorsSQ[i_nodeSQ] = surfColor[i_surf];
          }
        }
        GMSH_SQ(coordsSQ, colorsSQ, gmshfilecontent);
      }
      break;
    }
    // 20-noded hexahedron
    case DRT::Element::hex20:
    {
      const int n_surfNodes = 8;
      GMSH_GetSurfColor(element, n_surfNodes, surfNodes, surfColor);

      // 8-noded quadrangle (one 4-noded quadrangle and four 3-noded triangles)
      for (int i_surf = 0; i_surf < 6; i_surf++)
      {
        // 4-noded quadrangle
        int nodesSQ[4] = {4,5,6,7};
        double coordsSQ[3][4];
        double colorsSQ[4];
        for (int i_nodeSQ = 0; i_nodeSQ < 4; i_nodeSQ++)
        {
          for (int i_dim = 0; i_dim < 3; i_dim++)
          {
            coordsSQ[i_dim][i_nodeSQ] = coord(i_dim,surfNodes[i_surf][nodesSQ[i_nodeSQ]]);
            colorsSQ[i_nodeSQ] = surfColor[i_surf];
          }
        }
        GMSH_SQ(coordsSQ, colorsSQ, gmshfilecontent);

        // 3-noded triangles
        int nodesST[4][3] = {{0,4,7}, {1,5,4}, {2,6,5}, {3,7,6}};
        double coordsST[3][3];
        double colorsST[3];
        for (int i_ST = 0; i_ST < 4; i_ST++)
        {
          for (int i_nodeST = 0; i_nodeST < 3; i_nodeST++)
          {
            for (int i_dim = 0; i_dim < 3; i_dim++)
            {
              coordsST[i_dim][i_nodeST] = coord(i_dim,surfNodes[i_surf][nodesST[i_ST][i_nodeST]]);
              colorsST[i_nodeST] = surfColor[i_surf];
            }
          }
          GMSH_ST(coordsST, colorsST, gmshfilecontent);
        }
      }
      break;
    }
    default:
    {
      dserror("Gsmh-output: element shape not implemented");
      break;
    }
  }
}

/*----------------------------------------------------------------------*
 | Compute Gmsh output for solid surface element numbers                |
 *----------------------------------------------------------------------*/
void CONTACT::Beam3cmanager::GMSH_SolidSurfaceElementNumbers(
    const DRT::Element* element,
    const Epetra_Vector& disrow,
    std::stringstream& gmshfilecontent)
{
  // Prepare storage for nodal coordinates
  int nnodes = element->NumNode();
  LINALG::SerialDenseMatrix coord(3, nnodes);

  // Compute current nodal positions
  for (int i_dim = 0; i_dim < 3; i_dim++)
  {
    for (int i_node = 0; i_node < element->NumNode(); i_node++)
    {
      double referenceposition = ((element->Nodes())[i_node])->X()[i_dim];
      std::vector<int> dofnode = BTSolDiscret().Dof((element->Nodes())[i_node]);
      std::vector<int> problemdofnode;
      for (int m=0;m<(int)dofnode.size();++m)
        problemdofnode.push_back((dofoffsetmap_.find(dofnode[m]))->second);
      double displacement = disrow[ProblemDiscret().DofColMap()->LID(problemdofnode[i_dim])];
      coord(i_dim,i_node) = referenceposition + displacement;
    }
  }

  // center point of element
  double center[3] = {0.0, 0.0, 0.0};
  for (int k=0;k<nnodes;++k)
  {
    center[0] += coord(0,k);
    center[1] += coord(1,k);
    center[2] += coord(2,k);
  }
  center[0] /= nnodes;
  center[1] /= nnodes;
  center[2] /= nnodes;

  //Output of element IDs
  gmshfilecontent << "T3(" << std::scientific << center[0] << "," << center[1] << "," << center[2] << "," << 17 << ")";
  gmshfilecontent << "{\"" << element->Id() << "\"};" << std::endl;

  return;
}

/*----------------------------------------------------------------------*
 | Get color of solid element surfaces for GMSH-Output                  |
 *----------------------------------------------------------------------*/
void CONTACT::Beam3cmanager::GMSH_GetSurfColor(
    const DRT::Element* element,
    const int& n_surfNodes,
    const int surfNodes[6][9],
    double surfColor[6])
{
  // Loop over all surfaces to get their colors
  for (int i_surf = 0; i_surf < 6; i_surf++)
  {
    // Loop over all beam to solid contact pairs
    for (int i_pair = 0; i_pair < (int)btsolpairs_.size(); i_pair++)
    {
      int found = 0;

      // Get surface element
      DRT::Element* bts_element2 = const_cast<DRT::Element*>(btsolpairs_[i_pair]->Element2());

      // Loop over all nodes of bts surface element
      for (int j = 0; j < bts_element2->NumNode(); j++)
      {
        // Loop over all nodes of solid element surfaces
        for (int k = 0; k < n_surfNodes; k++)
        {
          if ((element->Nodes())[surfNodes[i_surf][k]]->Id() == (bts_element2->Nodes())[j]->Id())
          {
            found++;
          }
        }
      }
      // If the solid element contains all surface element nodes, the surface element is part the solid element
      if (found == bts_element2->NumNode())
      {
        // Get contact flag of beam to solid contact pair
        if (btsolpairs_[i_pair]->GetContactFlag())
        {
          // Beam to solid pair is in active contact
          surfColor[i_surf] = 1.0;

          // Go to the next surface
          break;
        }
        else
        {
          // Beam to solid pair might come into active contact
          surfColor[i_surf] = 0.5;
        }
      }
    }
  }
}

/*----------------------------------------------------------------------*
 | GMSH-Surface-Output for 4-noded quadrangle (SQ)                      |
 *----------------------------------------------------------------------*/
void CONTACT::Beam3cmanager::GMSH_SQ(
    const double coords[3][4],
    const double color[4],
    std::stringstream& gmshfilecontent)
{
  gmshfilecontent << "SQ(" << std::scientific;
  for (int i_node = 0; i_node < 4; i_node++)
  {
    for (int i_dim = 0; i_dim < 3; i_dim++)
    {
      gmshfilecontent << coords[i_dim][i_node];
      if (i_dim < 2)
      {
        // Coordinate x, y, z seperator
        gmshfilecontent << ",";
      }
    }
    if (i_node < 3)
    {
      // Node separator
      gmshfilecontent << ",";
    }
  }
  gmshfilecontent << "){" << std::scientific;
  for (int i_node = 0; i_node < 4; i_node++)
  {
    gmshfilecontent << color[i_node];
    if (i_node < 3)
    {
      // Node separator
      gmshfilecontent << ",";
    }
  }
  gmshfilecontent << "};" << std::endl;
}

/*----------------------------------------------------------------------*
 | GMSH-Surface-Output for 3-noded triangle (ST)                        |
 *----------------------------------------------------------------------*/
void CONTACT::Beam3cmanager::GMSH_ST(
    const double coords[3][3],
    const double color[3],
    std::stringstream& gmshfilecontent)
{
  gmshfilecontent << "ST(" << std::scientific;
  for (int i_node = 0; i_node < 3; i_node++)
  {
    for (int i_dim = 0; i_dim < 3; i_dim++)
    {
      gmshfilecontent << coords[i_dim][i_node];
      if (i_dim < 2)
      {
        // Coordinate x, y, z seperator
        gmshfilecontent << ",";
      }
    }
    if (i_node < 2)
    {
      // Node separator
      gmshfilecontent << ",";
    }
  }
  gmshfilecontent << "){" << std::scientific;
  for (int i_node = 0; i_node < 3; i_node++)
  {
    gmshfilecontent << color[i_node];
    if (i_node < 2)
    {
      // Node separator
      gmshfilecontent << ",";
    }
  }
  gmshfilecontent << "};" << std::endl;
}

/*----------------------------------------------------------------------*
 |  Determine number of nodes and nodal DoFs of element      meier 02/14|
 *----------------------------------------------------------------------*/
void CONTACT::Beam3cmanager::SetElementTypeAndDistype(DRT::Element* ele1)
{
  const DRT::ELEMENTS::Beam3Base* ele = dynamic_cast<const DRT::ELEMENTS::Beam3Base*>(ele1);

  numnodes_ = ele->NumCenterlineNodes();
  numnodalvalues_ = ele->HermiteCenterlineInterpolation() ? 2 : 1;

  return;
}

/*----------------------------------------------------------------------*
 | Is element midpoint distance smaller than search radius?  meier 02/14|
 *----------------------------------------------------------------------*/
bool CONTACT::Beam3cmanager::CloseMidpointDistance(const DRT::Element* ele1, const DRT::Element* ele2, std::map<int,LINALG::Matrix<3,1> >& currentpositions, const double sphericalsearchradius)
{

  LINALG::Matrix<3,1> midpos1(true);
  LINALG::Matrix<3,1> midpos2(true);
  LINALG::Matrix<3,1> diffvector(true);

  // get midpoint position of element 1
  if (ele1->NumNode() == 2)   //2-noded beam element
  {
    const DRT::Node* node1ele1 = ele1->Nodes()[0];
    const DRT::Node* node2ele1 = ele1->Nodes()[1];

    for (int i=0;i<3;++i)
      midpos1(i) = 0.5*((currentpositions[node1ele1->Id()])(i)+(currentpositions[node2ele1->Id()])(i));
  }
  else if (ele1->NumNode() == 1)  //rigidsphere element
  {
    const DRT::Node* node1ele1 = ele1->Nodes()[0];

    for (int i=0;i<3;++i)
      midpos1(i) = (currentpositions[node1ele1->Id()])(i);
  }

  // get midpoint position of element 2
  if (ele2->NumNode() == 2)   //2-noded beam element
  {
    const DRT::Node* node1ele2 = ele2->Nodes()[0];
    const DRT::Node* node2ele2 = ele2->Nodes()[1];

    for (int i=0;i<3;++i)
      midpos2(i) = 0.5*((currentpositions[node1ele2->Id()])(i)+(currentpositions[node2ele2->Id()])(i));
  }
  else if (ele2->NumNode() == 1)  //rigidsphere element
  {
    const DRT::Node* node1ele2 = ele2->Nodes()[0];

    for (int i=0;i<3;++i)
      midpos2(i) = (currentpositions[node1ele2->Id()])(i);
  }

  // compute distance
  for(int i=0;i<3;i++)
    diffvector(i) = midpos1(i) - midpos2(i);

  if(diffvector.Norm2()<=sphericalsearchradius)
    return true;
  else
    return false;
}

/*----------------------------------------------------------------------*
 | read contact force for restart  meier 02/15|
 *----------------------------------------------------------------------*/
void CONTACT::Beam3cmanager::ReadRestart(IO::DiscretizationReader& reader)
{

  reader.ReadVector(fcold_, "fcold");
  reader.ReadVector(dis_old_, "dis_old");
  totpenaltywork_ = reader.ReadDouble("totpenaltywork");

  maxtotalsimgap_ = reader.ReadDouble("maxtotalsimgap");
  maxtotalsimgap_cp_ = reader.ReadDouble("maxtotalsimgap_cp");
  maxtotalsimgap_gp_ = reader.ReadDouble("maxtotalsimgap_gp");
  maxtotalsimgap_ep_ = reader.ReadDouble("maxtotalsimgap_ep");
  maxtotalsimrelgap_ = reader.ReadDouble("maxtotalsimrelgap");
  mintotalsimgap_ = reader.ReadDouble("mintotalsimgap");
  mintotalsimgap_cp_ = reader.ReadDouble("mintotalsimgap_cp");
  mintotalsimgap_gp_ = reader.ReadDouble("mintotalsimgap_gp");
  mintotalsimgap_ep_ = reader.ReadDouble("mintotalsimgap_ep");
  mintotalsimrelgap_ = reader.ReadDouble("mintotalsimrelgap");
  mintotalsimunconvgap_ = reader.ReadDouble("mintotalsimunconvgap");

  outputcounter_ = reader.ReadInt("outputcounter");

  return;
}

/*----------------------------------------------------------------------*
 | write contact force for restart  meier 02/15|
 *----------------------------------------------------------------------*/
void CONTACT::Beam3cmanager::WriteRestart(Teuchos::RCP<IO::DiscretizationWriter> output)
{
  WriteRestart(*output);

  return;
}

/*----------------------------------------------------------------------*
 | write contact force for restart  meier 02/15|
 *----------------------------------------------------------------------*/
void CONTACT::Beam3cmanager::WriteRestart(IO::DiscretizationWriter& output)
{

  output.WriteVector("fcold", fcold_);
  output.WriteVector("dis_old", dis_old_);
  output.WriteDouble("totpenaltywork", totpenaltywork_);

  output.WriteDouble("maxtotalsimgap",maxtotalsimgap_);
  output.WriteDouble("maxtotalsimgap_cp",maxtotalsimgap_cp_);
  output.WriteDouble("maxtotalsimgap_gp",maxtotalsimgap_gp_);
  output.WriteDouble("maxtotalsimgap_ep",maxtotalsimgap_ep_);
  output.WriteDouble("maxtotalsimrelgap",maxtotalsimrelgap_);
  output.WriteDouble("mintotalsimgap",mintotalsimgap_);
  output.WriteDouble("mintotalsimgap_cp",mintotalsimgap_cp_);
  output.WriteDouble("mintotalsimgap_gp",mintotalsimgap_gp_);
  output.WriteDouble("mintotalsimgap_ep",mintotalsimgap_ep_);
  output.WriteDouble("mintotalsimrelgap",mintotalsimrelgap_);
  output.WriteDouble("mintotalsimunconvgap",mintotalsimunconvgap_);

  output.WriteInt("outputcounter", outputcounter_);  // ToDo needed?

  return;
}
