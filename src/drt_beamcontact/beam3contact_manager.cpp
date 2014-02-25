/*!-----------------------------------------------------------------------------------------------------------
\file beam3contact_manager.cpp
\brief Main class to control beam contact

<pre>
Maintainer: Christoph Meier
            meier@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15262
</pre>

*-----------------------------------------------------------------------------------------------------------*/

#include "beam3contact_manager.H"
#include "beam3contact.H"
#include "beam3contact_defines.H"
#include "beam3contact_octtree.H"
#include "../drt_inpar/inpar_beamcontact.H"
#include "../drt_inpar/inpar_contact.H"
#include "../linalg/linalg_utils.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_discret.H"
#include "../linalg/linalg_sparsematrix.H"
#include <Teuchos_Time.hpp>

#include "../drt_beam3/beam3.H"
#include "../drt_beam3ii/beam3ii.H"
#include "../drt_beam3eb/beam3eb.H"

/*----------------------------------------------------------------------*
 |  constructor (public)                                      popp 04/10|
 *----------------------------------------------------------------------*/
CONTACT::Beam3cmanager::Beam3cmanager(DRT::Discretization& discret, double alphaf):
numnodes_(0),
numnodalvalues_(0),
pdiscret_(discret),
pdiscomm_(discret.Comm()),
alphaf_(alphaf),
constrnorm_(0.0)
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
  cdiscret_ = Teuchos::rcp(new DRT::Discretization((std::string)"beam contact",comm));

  // loop over all column nodes of underlying problem discret and add
  for (int i=0;i<(ProblemDiscret().NodeColMap())->NumMyElements();++i)
  {
    DRT::Node* node = ProblemDiscret().lColNode(i);
    if (!node) dserror("Cannot find node with lid %",i);
    Teuchos::RCP<DRT::Node> newnode = Teuchos::rcp(node->Clone());
    ContactDiscret().AddNode(newnode);
  }

  // loop over all column elements of underlying problem discret and add
  for (int i=0;i<(ProblemDiscret().ElementColMap())->NumMyElements();++i)
  {
    DRT::Element* ele = ProblemDiscret().lColElement(i);
    if (!ele) dserror("Cannot find element with lid %",i);
    Teuchos::RCP<DRT::Element> newele = Teuchos::rcp(ele->Clone());
    ContactDiscret().AddElement(newele);
  }

  // build maps but do not assign dofs yet, we'll do this below
  // after shuffling around of nodes and elements (saves time)
  ContactDiscret().FillComplete(false,false,false);

  // store the node and element row and column maps into this manager
  noderowmap_ = Teuchos::rcp(new Epetra_Map(*(ContactDiscret().NodeRowMap())));
  elerowmap_  = Teuchos::rcp(new Epetra_Map(*(ContactDiscret().ElementRowMap())));
  nodecolmap_ = Teuchos::rcp(new Epetra_Map(*(ContactDiscret().NodeColMap())));
  elecolmap_  = Teuchos::rcp(new Epetra_Map(*(ContactDiscret().ElementColMap())));

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
    stproc.push_back(ContactDiscret().Comm().MyPID());
  if (elerowmap_->NumMyElements())
    estproc.push_back(ContactDiscret().Comm().MyPID());

  // information how many processors participate in total
  std::vector<int> allproc(ContactDiscret().Comm().NumProc());
  for (int i=0;i<ContactDiscret().Comm().NumProc();++i) allproc[i] = i;

  // declaring new variables into which the info of (e)stproc on all processors is gathered
  std::vector<int> rtproc(0);
  std::vector<int> ertproc(0);

  // gathers information of (e)stproc and writes it into (e)rtproc; in the end (e)rtproc
  // is a vector which contains the numbers of all processors which own nodes/elements.
  LINALG::Gather<int>(stproc,rtproc,ContactDiscret().Comm().NumProc(),&allproc[0],ContactDiscret().Comm());
  LINALG::Gather<int>(estproc,ertproc,ContactDiscret().Comm().NumProc(),&allproc[0],ContactDiscret().Comm());

  // in analogy to (e)stproc and (e)rtproc the variables (e)rdata gather all the row ID
  // numbers which are  stored on different processors in their own variables (e)sdata; thus,
  // each processor gets the information about all the row ID numbers existing in the problem
  std::vector<int> rdata;
  std::vector<int> erdata;

  // gather all gids of nodes redundantly from (e)sdata into (e)rdata
  LINALG::Gather<int>(sdata,rdata,(int)rtproc.size(),&rtproc[0],ContactDiscret().Comm());
  LINALG::Gather<int>(esdata,erdata,(int)ertproc.size(),&ertproc[0],ContactDiscret().Comm());

  // build completely overlapping node map (on participating processors)
  Teuchos::RCP<Epetra_Map> newnodecolmap = Teuchos::rcp(new Epetra_Map(-1,(int)rdata.size(),&rdata[0],0,ContactDiscret().Comm()));
  sdata.clear();
  stproc.clear();
  rdata.clear();
  allproc.clear();

  // build completely overlapping element map (on participating processors)
  Teuchos::RCP<Epetra_Map> newelecolmap = Teuchos::rcp(new Epetra_Map(-1,(int)erdata.size(),&erdata[0],0,ContactDiscret().Comm()));
  esdata.clear();
  estproc.clear();
  erdata.clear();

  // store the fully overlapping node and element maps
  nodefullmap_ = Teuchos::rcp(new Epetra_Map(*newnodecolmap));
  elefullmap_  = Teuchos::rcp(new Epetra_Map(*newelecolmap));

  // pass new fully overlapping node and element maps to beam contact discretization
  ContactDiscret().ExportColumnNodes(*newnodecolmap);
  ContactDiscret().ExportColumnElements(*newelecolmap);


  // complete beam contact discretization based on the new column maps
  // (this also assign new degrees of freedom what we actually do not
  // want, thus we have to introduce a dof mapping next)
  ContactDiscret().FillComplete(true,false,false);

  // compute dof offset between problem and contact discretization
  const Epetra_Map* pddofs = ProblemDiscret().DofRowMap();
  const Epetra_Map* cddofs = ContactDiscret().DofRowMap();
  if (pddofs->NumGlobalElements() != cddofs->NumGlobalElements())
    dserror("ERROR: Inconsistency in dof mapping between the two discretizations");
  int pdmin = pddofs->MinAllGID();
  int cdmin = cddofs->MinAllGID();
  dofoffset_ = cdmin - pdmin;

  // read parameter lists from DRT::Problem
  sbeamcontact_ = DRT::Problem::Instance()->BeamContactParams();
  scontact_     = DRT::Problem::Instance()->ContactDynamicParams();
  
  // check input parameters
  if (DRT::INPUT::IntegralValue<INPAR::CONTACT::ApplicationType>(scontact_,"APPLICATION") != INPAR::CONTACT::app_beamcontact)
   dserror("ERROR: The given input parameters are not for beam contact");
  if (scontact_.get<double>("PENALTYPARAM") <= 0.0)
   dserror("ERROR: The penalty parameter has to be positive.");  
  
  // initialize contact element pairs
  pairs_.resize(0);
  oldpairs_.resize(0);
  
  // initialize input parameters
  currentpp_ = scontact_.get<double>("PENALTYPARAM");
  
  // initialize Uzawa iteration index
  uzawaiter_ = 0;

  if(!pdiscret_.Comm().MyPID())
  {
    std::cout << "========================= Beam Contact =========================" << std::endl;
    std::cout<<"Elements in discret.   = "<<pdiscret_.NumGlobalElements()<<std::endl;
  }

  // initialize octtree for contact search
  if (DRT::INPUT::IntegralValue<INPAR::BEAMCONTACT::OctreeType>(sbeamcontact_,"BEAMS_OCTREE") != INPAR::BEAMCONTACT::boct_none)
  {
    if (!pdiscret_.Comm().MyPID())
      std::cout << "Penalty parameter      = " << currentpp_ << std::endl;
    tree_ = Teuchos::rcp(new Beam3ContactOctTree(sbeamcontact_,pdiscret_,*cdiscret_,dofoffset_));
  }
  else
  {
    // Compute the search radius for searching possible contact pairs
    ComputeSearchRadius();
    tree_ = Teuchos::null;
    if(!pdiscret_.Comm().MyPID())
      std::cout<<"\nBrute Force Search"<<std::endl;
  }



  if(!pdiscret_.Comm().MyPID())
  {
    if (DRT::INPUT::IntegralValue<INPAR::CONTACT::SolvingStrategy>(scontact_,"STRATEGY") == INPAR::CONTACT::solution_penalty )
      std::cout << "Strategy               = Penalty" << std::endl;
    else if (DRT::INPUT::IntegralValue<INPAR::CONTACT::SolvingStrategy>(scontact_,"STRATEGY") == INPAR::CONTACT::solution_auglag)
      std::cout << "Strategy               = Augmented Lagrange" << std::endl;
    std::cout <<"================================================================\n" << std::endl;
  }

  // initialize newgapfunction_ flag
  newgapfunction_ = DRT::INPUT::IntegralValue<int>(BeamContactParameters(),"BEAMS_NEWGAP");

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
                                      bool newsti)
{
  // map linking node numbers and current node positions
  std::map<int,LINALG::Matrix<3,1> > currentpositions;
  currentpositions.clear();
  
  // update currentpositions
  SetCurrentPositions(currentpositions,disrow);

  std::vector<std::vector<DRT::Element* > > contactelementpairs;
  contactelementpairs.clear();

  //**********************************************************************
  // octtree search (loop over all elements and find closest pairs)
  //**********************************************************************
  if(tree_ != Teuchos::null)
  {
     double t_start = Teuchos::Time::wallTime();

     contactelementpairs = tree_->OctTreeSearch(currentpositions);

     double t_end = Teuchos::Time::wallTime() - t_start;
     Teuchos::ParameterList ioparams = DRT::Problem::Instance()->IOParams();
     if(!pdiscret_.Comm().MyPID() && ioparams.get<int>("STDOUTEVRY",0))
       std::cout << "           Octree Search: " << t_end << " seconds, "<< std::endl;
  }
  //**********************************************************************
  // brute-force search (loop over all elements and find closest pairs)
  //**********************************************************************
  else
  {
    double t_start = Teuchos::Time::wallTime();

    contactelementpairs = SearchPossibleContactPairsNew(currentpositions);
    double t_end = Teuchos::Time::wallTime() - t_start;
    Teuchos::ParameterList ioparams = DRT::Problem::Instance()->IOParams();
    if(!pdiscret_.Comm().MyPID() && ioparams.get<int>("STDOUTEVRY",0))
      std::cout << "\nBrute Force Search: " << t_end << " seconds, "<< std::endl;
  }

  //Determine type of applied beam elements and set the corresponding values for the member variables
  //numnodes_ and numnodaldofs_. This has only to be done once in the beginning, since beam contact
  //simulations are only possible when using beam elements of one type!
  if (oldpairs_.size()==0 and contactelementpairs.size()>0)
    SetElementTypeAndDistype((contactelementpairs[0])[0]);

  //So far, all beam elements occuring in the pairs_ vector have to be of the same type.
  //This will be checked in the following lines.
  if ((int)contactelementpairs.size()>0)
  {
    const DRT::ElementType & pair1_ele1_type = ((contactelementpairs[0])[0])->ElementType();
    for (int k=0;k<(int)contactelementpairs.size();++k)
    {
      const DRT::ElementType & ele1_type = ((contactelementpairs[k])[0])->ElementType();
      const DRT::ElementType & ele2_type = ((contactelementpairs[k])[1])->ElementType();

      if ( ele1_type!=pair1_ele1_type or ele2_type!=pair1_ele1_type)
      {
        dserror("All contacting beam elements have to be of the same type (beam3eb, beam3 or beam3ii). Change your input file!");
      }
    }
  }

  // Only the element pairs of contactelementpairs (found in the contact search) which have not been found in the last time step (i.e. which
  // are not in oldpairs_) will be generated as new Beam3contact instances. Pairs which already exist in oldpairs_ will simply be copied to pairs_.
  // This procedure looks a bit circumstantial at first glance: However, it is not possible to solely use contactelementpairs
  // and to simply delete oldpairs_ at the end of a time step if the new gap function definition is used, since the latter needs
  // history variables of the last time step which are stored in the oldpairs_ vector.
  for (int k=0;k<(int)contactelementpairs.size();k++)
  {
    DRT::Element* ele1 = (contactelementpairs[k])[0];
    DRT::Element* ele2 = (contactelementpairs[k])[1];
    int currid1 = ele1->Id();
    int currid2 = ele2->Id();
    bool foundlasttimestep = false;
    for (int m=0;m<(int)oldpairs_.size();m++)
    {
      int id1 = (oldpairs_[m]->Element1())->Id();
      int id2 = (oldpairs_[m]->Element2())->Id();
      // pair already exists
      if ((id1 == currid1 && id2 == currid2) || (id1 == currid2 && id2 == currid1))
      {
        //Add existing pair of last time step (storded in oldpairs_) to new pairs_ vector
        //So far, this is only done if the new gap function is applied.
        if (!PairAllreadyExisting(currid1,currid2) and newgapfunction_)
        {
          pairs_.push_back(oldpairs_[m]);
          contactpairmap_[std::make_pair(id1,id2)] = pairs_[pairs_.size()-1];
          foundlasttimestep = true;
        }
        break;
      }
    }
    // create new contact pair and add to pairs_ if not found before
    if (!foundlasttimestep && !PairAllreadyExisting(currid1,currid2))
    {
      #ifndef NEWBEAMCONTACT
        //add new contact pair object
        pairs_.push_back(Teuchos::rcp (new Beam3contact(CONTACT::Beam3contact(ProblemDiscret(),ContactDiscret(),DofOffset(),ele1,ele2))));
      #else
        //Add new contact pair objcect: The auxiliary_instance of the abstract class Beam3contactinterface is only needed here in order to call
        //the function Impl() which creates an instance of the templated class Beam3contactnew<numnodes, numnodalvalues> !
        pairs_.push_back(CONTACT::Beam3contactinterface::Impl(numnodes_,numnodalvalues_,ProblemDiscret(),ContactDiscret(),DofOffset(),ele1,ele2, sbeamcontact_));
      #endif
      contactpairmap_[std::make_pair(currid1,currid2)] = pairs_[pairs_.size()-1];
    }
  }

  // update element state of all pairs with current positions (already calculated in SetCurrentPositions) and current tangents (will be calculated in SetState)
  SetState(currentpositions,disrow);
  // At this point we have all possible contact pairs_ with updated positions

  //**********************************************************************
  // evalutation of contact pairs
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
  fc_ = Teuchos::rcp(new Epetra_Vector(fres.Map()));
  if (fcold_==Teuchos::null) fcold_ = Teuchos::rcp(new Epetra_Vector(fres.Map()));
  
  // initialize contact stiffness and uncomplete global stiffness
  stiffc_ = Teuchos::rcp(new LINALG::SparseMatrix(stiffmatrix.RangeMap(),100));
  stiffmatrix.UnComplete();

#ifndef NEWBEAMCONTACT
  const Teuchos::ParameterList& sbeamcontact = DRT::Problem::Instance()->BeamContactParams();
  int beams_smoothing = DRT::INPUT::IntegralValue<INPAR::BEAMCONTACT::Smoothing>(sbeamcontact,"BEAMS_SMOOTHING");
  // Loop over all pairs
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
      pairs_[i]->PreEvaluation(beams_smoothing, currentpositions);
  }

  // Loop over all pairs
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
      bool newgapfunction = DRT::INPUT::IntegralValue<int>(BeamContactParameters(),"BEAMS_NEWGAP");
      pairs_[i]->Evaluate(*stiffc_,*fc_,currentpp_,newgapfunction,contactpairmap_,beams_smoothing);
    }
  }
#else
  // Loop over all pairs
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
      pairs_[i]->Evaluate(*stiffc_,*fc_,currentpp_);
    }
  }
#endif

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

  // debug output
  //std::cout << "fc:    " << std::endl << *fc_ << std::endl;
  //std::cout << "fcold: " << std::endl << *fcold_ << std::endl;
  //std::cout << "fres:  " << std::endl << fres << std::endl;

  return;
}

/*----------------------------------------------------------------------*
 |  Set current displacement state                            popp 04/10|
 *----------------------------------------------------------------------*/
void CONTACT::Beam3cmanager::SetCurrentPositions(std::map<int,LINALG::Matrix<3,1> >& currentpositions,
                                                 const Epetra_Vector& disrow)
{
  //**********************************************************************
  // get positions of all nodes (fully overlapping map) and store
  // the current results into currentpositions
  //**********************************************************************

  // we need to address the nodes / dofs via the beam contact
  // discretization, because only this is exported to full overlap
  Epetra_Vector disbcrow(disrow);
  disbcrow.ReplaceMap(*ContactDiscret().DofRowMap());

  // export displacements into fully overlapping column map format
  Epetra_Vector discol(*ContactDiscret().DofColMap(),true);
  LINALG::Export(disbcrow,discol);

  // loop over all beam contact nodes
  for (int i=0;i<FullNodes()->NumMyElements();++i)
  {
    // get node pointer
    const DRT::Node* node = ContactDiscret().lColNode(i);

    // get GIDs of this node's degrees of freedom
    std::vector<int> dofnode = ContactDiscret().Dof(node);

    // nodal positions
    LINALG::Matrix<3,1> currpos;
    currpos(0) = node->X()[0] + discol[ContactDiscret().DofColMap()->LID(dofnode[0])];
    currpos(1) = node->X()[1] + discol[ContactDiscret().DofColMap()->LID(dofnode[1])];
    currpos(2) = node->X()[2] + discol[ContactDiscret().DofColMap()->LID(dofnode[2])];

    // store into currentpositions
    currentpositions[node->Id()] = currpos;
  }

  return;
}

/*----------------------------------------------------------------------*
 |  Set current displacement state                            popp 04/10|
 *----------------------------------------------------------------------*/
void CONTACT::Beam3cmanager::SetState(std::map<int,LINALG::Matrix<3,1> >& currentpositions,
                                      const Epetra_Vector& disrow)
{
  //**********************************************************************
  // get positions of all nodes (fully overlapping map) and store
  // the current results into currentpositions
  //**********************************************************************

  // we need to address the nodes / dofs via the beam contact
  // discretization, because only this is exported to full overlap
  Epetra_Vector disbcrow(disrow);
  disbcrow.ReplaceMap(*ContactDiscret().DofRowMap());

  // export displacements into fully overlapping column map format
  Epetra_Vector discol(*ContactDiscret().DofColMap(),true);
  LINALG::Export(disbcrow,discol);

  // map to store the nodal tangent vectors (necessary for Kirchhoff type beams) and adress it with the node ID
  std::map<int,LINALG::Matrix<3,1> > currenttangents;
  currenttangents.clear();

  // Update of nodal tangents for Kirchhoff elements; nodal positions have already been set in SetCurrentPositions
  // loop over all beam contact nodes
  for (int i=0;i<FullNodes()->NumMyElements();++i)
  {
    // get node pointer
    const DRT::Node* node = ContactDiscret().lColNode(i);

    //get nodal tangents for Kirchhoff elements
    if (numnodalvalues_==2)
    {
      // get GIDs of this node's degrees of freedom
      std::vector<int> dofnode = ContactDiscret().Dof(node);

      LINALG::Matrix<3,1> currtan(true);
      for (int i=0; i<node->Elements()[0]->NumNode();i++)
      {
        if (node->Elements()[0]->Nodes()[i]->Id()==node->Id())
        {
          const DRT::ELEMENTS::Beam3eb* ele = dynamic_cast<const DRT::ELEMENTS::Beam3eb*>(node->Elements()[0]);
          currtan(0)=((ele->Tref())[i])(0) + discol[ContactDiscret().DofColMap()->LID(dofnode[3])];
          currtan(1)=((ele->Tref())[i])(1) + discol[ContactDiscret().DofColMap()->LID(dofnode[4])];
          currtan(2)=((ele->Tref())[i])(2) + discol[ContactDiscret().DofColMap()->LID(dofnode[5])];
        }
      }
      // store into currenttangents
      currenttangents[node->Id()] = currtan;
    }
    //set currenttangents to zero for Reissner elements
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
    for(int m=0;m<(pairs_[i]->Element1())->NumNode();m++)
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
      for(int m=0;m<(pairs_[i]->Element1())->NumNode();m++)
      {
        int tempGID = ((pairs_[i]->Element1())->NodeIds())[m];
        LINALG::Matrix<3,1> temptan = currenttangents[tempGID];

        // store updated nodal tangents
        for(int n=0;n<3;n++)
          ele1pos(n+3,m) = temptan(n);
      }
    }
    // Positions: Loop over all nodes of element 2
    for(int m=0;m<(pairs_[i]->Element2())->NumNode();m++)
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
      for(int m=0;m<(pairs_[i]->Element2())->NumNode();m++)
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

#ifdef NEWBEAMCONTACT
  // Update also the interpolated tangents if the tangentsmoothing is activated for Reissner beams
  int smoothing = DRT::INPUT::IntegralValue<INPAR::BEAMCONTACT::Smoothing>(sbeamcontact_,"BEAMS_SMOOTHING");
  if (smoothing != INPAR::BEAMCONTACT::bsm_none)
  {
    for(int i=0;i<(int)pairs_.size();++i)
    {
      pairs_[i]->UpdateEleSmoothTangents(currentpositions);
    }
  }
#endif

  return;
}

/*----------------------------------------------------------------------*
 |  search possible contact element pairs                     popp 04/10|
 *----------------------------------------------------------------------*/
std::vector<std::vector<DRT::Element*> > CONTACT::Beam3cmanager::SearchPossibleContactPairsNew(std::map<int,LINALG::Matrix<3,1> >& currentpositions)
{
  //**********************************************************************
  // Steps of search for element pairs that MIGHT get into contact:
  //
  // 1) Find non-neighbourng node pairs
  // 2) Compute distance between node pairs and compare with search radius
  // 3) Find non-neighbouring element pairs based on node pairs
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
    DRT::Node* firstnode = ContactDiscret().gNode(firstgid);
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
       if (distance.Norm2() < searchradius_)
         NearNodesGIDs.push_back(secondgid);
     }
   }

   // AT THIS POINT WE HAVE FOUND AND STROED ALL NODES CLOSE TO FIRSTNODE BESIDES THE DIRECTLY CONNECTED NEIGHBOR NODES!

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
     DRT::Node* tempnode = ContactDiscret().gNode(NearNodesGIDs[j]);

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
     DRT::Element* ele1 = ContactDiscret().gElement(FirstElesGIDs[j]);

     // node ids adjacent to this element
     const int* NodesEle1 = ele1->NodeIds();

     // loop over all close elements
     for (int k=0;k<(int)SecondElesGIDsRej.size();++k)
      {
       // get and cast a pointer on an element
        DRT::Element* ele2 = ContactDiscret().gElement(SecondElesGIDsRej[k]);

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

        // Check if the pair 'jk' already exists in pairs
        bool foundbefore = false;
        for (int m=0;m<(int)newpairs.size();++m)
        {
          int id1 = (newpairs[m])[0]->Id();
          int id2 = (newpairs[m])[1]->Id();
          int currid1 = FirstElesGIDs[j];
          int currid2 = SecondElesGIDsRej[k];

         // already exists if can be found in pairs
          if ((id1 == currid1 && id2 == currid2) || (id1 == currid2 && id2 == currid1))
            foundbefore = true;
        }

        // if NOT neighbouring and NOT found before
        // create new beam3contact object and store it into pairs_
        if (!elements_neighbouring && !foundbefore)
        {
          std::vector<DRT::Element*> contactelementpair;
          contactelementpair.clear();
          contactelementpair.push_back(ele1);
          contactelementpair.push_back(ele2);
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
  // TODO more sophisticated criteria might be added later
  
  // some local variables
  double charactlength = 0.0;
  double globalcharactlength = 0.0;
  double maxeleradius = 0.0;
  double maxelelength = 0.0;
    
  // look for maximum element radius in the whole discretization
  GetMaxEleRadius(maxeleradius);
  
  // look for maximum element length in the whole discretization
  GetMaxEleLength(maxelelength);
  
  // select characeteristic length
  if (maxeleradius > maxelelength) charactlength = maxeleradius;   
  else                             charactlength = maxelelength;  
  
  // communicate among all procs to find the global maximum
  Comm().MaxAll(&charactlength,&globalcharactlength,1);
  
  // compute the search radius
  // the factor is only empiric yet it must be greater than 2
  // in order to account for the circular beam cross sections
  searchradius_ = 2.5 * globalcharactlength;

  // some information for the user
  if (Comm().MyPID()==0)
  {
    std::cout << "Penalty parameter      = " << currentpp_ << std::endl;
    std::cout << "Maximum element radius = " << maxeleradius << std::endl;
    std::cout << "Maximum element length = " << maxelelength << std::endl;
    std::cout << "Search radius          = " << searchradius_  << std::endl << std::endl;
  }
  
  return;
}

/*-------------------------------------------------------------------- -*
 |  Find maximum element radius in discretization             popp 04/10|
 *----------------------------------------------------------------------*/
void CONTACT::Beam3cmanager::GetMaxEleRadius(double& maxeleradius)
{
  // loop over all elements in row map 
  for (int i=0;i<RowElements()->NumMyElements();++i)
  {
    // get pointer onto element
    int gid = RowElements()->GID(i);
    DRT::Element* thisele = ContactDiscret().gElement(gid);
    

    // compute eleradius from moment of inertia
    // (RESTRICTION: CIRCULAR CROSS SECTION !!!)

    double eleradius = 0;

    const DRT::ElementType & eot = thisele->ElementType();

    if ( eot == DRT::ELEMENTS::Beam3Type::Instance() )
    {
      DRT::ELEMENTS::Beam3* thisbeam = static_cast<DRT::ELEMENTS::Beam3*>(thisele);
      eleradius = sqrt(sqrt(4 * (thisbeam->Izz()) / M_PI));
    }
    if ( eot == DRT::ELEMENTS::Beam3iiType::Instance() )
    {
      DRT::ELEMENTS::Beam3ii* thisbeam = static_cast<DRT::ELEMENTS::Beam3ii*>(thisele);
      eleradius = sqrt(sqrt(4 * (thisbeam->Izz()) / M_PI));
    }
    if ( eot == DRT::ELEMENTS::Beam3ebType::Instance() )
    {
      DRT::ELEMENTS::Beam3eb* thisbeam = static_cast<DRT::ELEMENTS::Beam3eb*>(thisele);
      eleradius = sqrt(sqrt(4 * (thisbeam->Iyy()) / M_PI));
    }
    
    // if current radius is larger than maximum radius -> update
    if (eleradius > maxeleradius) maxeleradius = eleradius;
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
    DRT::Element* thisele = ContactDiscret().gElement(gid);
    
    // get global IDs of edge nodes and pointers
    int node0_gid = thisele->NodeIds()[0];
    int node1_gid = thisele->NodeIds()[1];
    DRT::Node* node0 = ContactDiscret().gNode(node0_gid);
    DRT::Node* node1 = ContactDiscret().gNode(node1_gid);
    
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
    double elelength = sqrt(dist[0]*dist[0]+dist[1]*dist[1]+dist[2]*dist[2]);
    
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

  // create gmsh output for nice visualization
  // (endoftimestep flag set to TRUE)
#ifdef GMSHTIMESTEPS
  GmshOutput(disrow, timestep, newtonstep, true);
#endif
  
  // print some data to screen
  ConsoleOutput();

  // set normal_old_=normal_ for all contact pairs at the end of the time step
  if(newgapfunction_) ShiftAllNormal();

  //store pairs_ in oldpairs_ to be available in next time step
  //this is needed for the new gapfunction defintion and also for the output at the end of an time step
  oldpairs_.clear();
  oldpairs_.resize(0);
  oldpairs_ = pairs_;

  // clear potential contact pairs
  pairs_.clear();
  pairs_.resize(0);

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
  
  // STEP 1: OUTPUT OF TIME STEP INDEX
  std::ostringstream filename;
  filename << "o/gmsh_output/";
  if (timestep<1000000)
    filename << "beams_t" << std::setw(6) << std::setfill('0') << timestep;
  else /*(timestep>=1000000)*/
    dserror("ERROR: Gmsh output implemented for max 999.999 time steps");
  
  // STEPS 2/3: OUTPUT OF UZAWA AND NEWTON STEP INDEX
  // (for the end of time step output, omit this)
  int uzawastep = 99;
  if (!endoftimestep)
  {
    // check if Uzawa index is needed or not
    INPAR::CONTACT::SolvingStrategy soltype =
    DRT::INPUT::IntegralValue<INPAR::CONTACT::SolvingStrategy>(GeneralContactParameters(),"STRATEGY");

    if (soltype==INPAR::CONTACT::solution_auglag)
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
      filename << "_n0";
    else if (newtonstep<100)
      filename << "_n";
    else /*(newtonstep>=100*/
      dserror("ERROR: Gmsh output implemented for max 99 Newton steps");
    filename << newtonstep;
  }
  
  // finish filename
  filename<<".pos";

  //**********************************************************************
  // gsmh output beam elements as prisms
  // NOT YET AVAILABLE FOR THE PARALLEL CASE!
  //**********************************************************************
  if (Comm().NumProc() == 1)
  {
    // approximation of the circular cross section with n prisms
    int n=N_CIRCUMFERENTIAL;        // CHOOSE N HERE
    if (n<3) n=3;    // minimum is 3, else no volume is defined.
    int n_axial=N_Axial; //number of devisions of element in axial direction

    // we need to address the nodes / dofs via the beam contact
    // discretization, because only this is exported to full overlap
    Epetra_Vector disbcrow(disrow);
    disbcrow.ReplaceMap(*ContactDiscret().DofRowMap());

    // export displacements into fully overlapping column map format
    Epetra_Vector discol(*ContactDiscret().DofColMap(),true);
    LINALG::Export(disbcrow,discol);

    // do output to file in c-style
    FILE* fp = NULL;
    
    // open file to write output data into
    fp = fopen(filename.str().c_str(), "w");
  
    std::stringstream gmshfileheader;
    // write output to temporary std::stringstream; 
    gmshfileheader <<"General.BackgroundGradient = 0;\n";
    gmshfileheader <<"View.LineType = 1;\n";
    gmshfileheader <<"View.LineWidth = 1.4;\n";
    gmshfileheader <<"View.PointType = 1;\n";
    gmshfileheader <<"View.PointSize = 3;\n";
    gmshfileheader <<"General.ColorScheme = 1;\n";
    gmshfileheader <<"General.Color.Background = {255,255,255};\n";
    gmshfileheader <<"General.Color.Foreground = {255,255,255};\n";
    //gmshfileheader <<"General.Color.BackgroundGradient = {128,147,255};\n";
    gmshfileheader <<"General.Color.Foreground = {85,85,85};\n";
    gmshfileheader <<"General.Color.Text = {0,0,0};\n";
    gmshfileheader <<"General.Color.Axes = {0,0,0};\n";
    gmshfileheader <<"General.Color.SmallAxes = {0,0,0};\n";
    gmshfileheader <<"General.Color.AmbientLight = {25,25,25};\n";
    gmshfileheader <<"General.Color.DiffuseLight = {255,255,255};\n";
    gmshfileheader <<"General.Color.SpecularLight = {255,255,255};\n";
    gmshfileheader <<"View.ColormapAlpha = 1;\n";
    gmshfileheader <<"View.ColormapAlphaPower = 0;\n";
    gmshfileheader <<"View.ColormapBeta = 0;\n";
    gmshfileheader <<"View.ColormapBias = 0;\n";
    gmshfileheader <<"View.ColormapCurvature = 0;\n";
    gmshfileheader <<"View.ColormapInvert = 0;\n";
    gmshfileheader <<"View.ColormapNumber = 2;\n";
    gmshfileheader <<"View.ColormapRotation = 0;\n";
    gmshfileheader <<"View.ColormapSwap = 0;\n";
    gmshfileheader <<"View.ColorTable = {Black,Yellow,Blue,Orange,Red,Cyan,Purple,Brown,Green};\n";

    //write content into file and close it
    fprintf(fp, gmshfileheader.str().c_str());
    fclose(fp);

    fp = fopen(filename.str().c_str(), "w");

    std::stringstream gmshfilecontent;
    gmshfilecontent << "View \" Step T" << timestep;
            
    // information about Uzawa and Newton step
    if (!endoftimestep)
      gmshfilecontent << " U" << uzawastep << " N" << newtonstep;

    // finish step information
    gmshfilecontent << " \" {" << std::endl;
    
    // loop over all column elements on this processor
    for (int i=0;i<ColElements()->NumMyElements();++i)
    {
      // get pointer onto current beam element
      DRT::Element* element = ContactDiscret().lColElement(i);

      // write output in specific function
      const DRT::ElementType & eot = element->ElementType();
      // standard procedure for Reissner beams
      if ( eot == DRT::ELEMENTS::Beam3Type::Instance() or eot == DRT::ELEMENTS::Beam3iiType::Instance() )
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
            std::vector<int> dofnode = ContactDiscret().Dof((element->Nodes())[jd]);
            double displacement = discol[ContactDiscret().DofColMap()->LID(dofnode[id])];
            coord(id,jd) =  referenceposition + displacement;
          }
        }

        switch (element->NumNode())
        {
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
          // 4- or 5-noded beam element (higher-order interpolation)
          default:
          {
            dserror("Gmsh output for %i noded element not yet implemented!", element->NumNode());
            break;
          }
        }
      }
      // Kirchhoff beams need a special treatment
      else if ( eot == DRT::ELEMENTS::Beam3ebType::Instance())
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
            std::vector<int> dofnode = ContactDiscret().Dof((element->Nodes())[j]);
            double displacement = discol[ContactDiscret().DofColMap()->LID(dofnode[i])];
            nodalcoords(i,j) =  referenceposition + displacement;
            nodaltangents(i,j) =  ((ele->Tref())[j])(i) + discol[ContactDiscret().DofColMap()->LID(dofnode[3+i])];
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
        GMSH_N_noded(n,n_axial,coord,element,gmshfilecontent);
      }
      else
      {
        dserror("Your chosen type of beam element is not allowed for beam contact!");
      }
    }
    
    // finish data section of this view by closing curley brackets (somehow needed to get color)
    gmshfilecontent <<"SP(0.0,0.0,0.0){0.0,0.0};"<<std::endl;
    gmshfilecontent << "};" << std::endl;
    
    // write content into file and close
    fprintf(fp,gmshfilecontent.str().c_str());
    fclose(fp); 
  }
 

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

  // remove contact stiffness terms from stiffmatrix
  stiffmatrix.Add(*stiffc_, false, -scalemat, 1.0);

  // remove old contact force terms from fres
  fres.Update(-(1.0-alphaf_),*fc_,1.0);
  fres.Update(-alphaf_,*fcold_,1.0);

  // now redo Evaluate()
  Teuchos::RCP<Epetra_Vector> nullvec = Teuchos::null;
  Evaluate(stiffmatrix,fres,disrow,newsti);

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
  
  return;
}

/*----------------------------------------------------------------------*
 |  Update contact constraint norm                            popp 04/10|
 *----------------------------------------------------------------------*/
void CONTACT::Beam3cmanager::UpdateConstrNorm(double* cnorm)
{
  // some local variables
  int j=0;
  int dim = (int)pairs_.size();
  // vector holding the values of the penetration relative to the radius of the beam with the smaller radius
  Epetra_SerialDenseVector gapvector(dim);
  
  // initalize processor-local and global norms
  double locnorm = 0.0;
  double globnorm = 0.0;
  
  // loop over all pairs to find all gaps
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
        std::cout << "Penetration to large, choose higher penalty parameter!" << std::endl;
      }

      #ifdef NEWBEAMCONTACT
        if (pairs_[i]->GetShiftStatus() == true)
        {
        //In case the contact points of two beams are identical (i.e. r1=r2), the nodal coordinates of one beam are shifted by a small
        //predefined value in order to enable evaluation of the contact pair. This makes the Newton scheme more robust. However,
        //in the converged configuration we want to have the real nodal positions for all contact pairs!!!
        std::cout << "Contact pair with identical contact points (i.e. r1=r2) not possible in the converged configuration!" << std::endl;
        }
      #endif

#ifdef RELCONSTRTOL
      // Retrieve beam radii
      std::vector<double> radii(2,0.0);
      for(int k=0; k<(int)radii.size();k++)
      {
        // get Element type
        DRT::Element* currele = cdiscret_->gElement(pairs_[i]->Element1()->Id());
        if(k>0)
          currele = cdiscret_->gElement(pairs_[i]->Element2()->Id());

        const DRT::ElementType & eot = currele->ElementType();
        if(eot == DRT::ELEMENTS::Beam3Type::Instance())
          radii[k] = sqrt(sqrt(4 * ((dynamic_cast<DRT::ELEMENTS::Beam3*>(currele))->Izz()) / M_PI));
        else if(eot == DRT::ELEMENTS::Beam3iiType::Instance())
          radii[k] = sqrt(sqrt(4 * ((dynamic_cast<DRT::ELEMENTS::Beam3ii*>(currele))->Izz()) / M_PI));
        else if(eot == DRT::ELEMENTS::Beam3ebType::Instance())
                  radii[k] = sqrt(sqrt(4 * ((dynamic_cast<DRT::ELEMENTS::Beam3eb*>(currele))->Iyy()) / M_PI));
        if(radii[k]==0.0)
          dserror("beam radius is 0! Check your element type.");
      }

      double smallerradius = std::min(radii[0], radii[1]);

      gapvector[j] = pairs_[i]->GetGap()/smallerradius;
#else
      gapvector[j] = pairs_[i]->GetGap();
#endif
      j++;
    }
  }
  
  // compute local constraint norm as a L-inf-norm
  locnorm = gapvector.NormInf();
  
  // communicate among all procs to find the global constraint norm
  Comm().MaxAll(&locnorm,&globnorm,1);
  
  // update penalty parameter if necessary
  // (only possible for AUGMENTED LAGRANGE strategy)
  bool updatepp = false;
  INPAR::CONTACT::SolvingStrategy soltype = DRT::INPUT::IntegralValue<INPAR::CONTACT::SolvingStrategy>(GeneralContactParameters(),"STRATEGY");
  if (soltype==INPAR::CONTACT::solution_auglag)
      updatepp = IncreaseCurrentpp(globnorm);
  
   // print results to screen
  Teuchos::ParameterList ioparams = DRT::Problem::Instance()->IOParams();
  if (Comm().MyPID()==0 && cnorm==NULL && ioparams.get<int>("STDOUTEVRY",0))
  {
    std::cout << std::endl << "*********************************************"<<std::endl;
    std::cout << "Global Constraint Norm = " << globnorm << std::endl;
    if (updatepp) std::cout << "Updated penalty parameter = " << currentpp_ << std::endl;
    std::cout<<"*********************************************"<<std::endl;
  }
    
  // update class variable
  if(cnorm==NULL)
    constrnorm_ = globnorm;
  else
    *cnorm = globnorm;
  
  return; 
}

/*----------------------------------------------------------------------*
 |  Shift normal vector to "normal_old_"                     meier 10/10|
 *----------------------------------------------------------------------*/
void CONTACT::Beam3cmanager::ShiftAllNormal()
{
   // loop over all potential contact pairs
    for (int i=0;i<(int)pairs_.size();++i)
      pairs_[i]->ShiftNormal();


}

/*----------------------------------------------------------------------*
 |  Update all Uzawa-based Lagrange multipliers               popp 04/10|
 *----------------------------------------------------------------------*/
void CONTACT::Beam3cmanager::UpdateAlllmuzawa()
{
  // loop over all potential contact pairs
  for (int i=0;i<(int)pairs_.size();++i)
    pairs_[i]->Updatelmuzawa(currentpp_);
    
  return;
}

/*----------------------------------------------------------------------*
 |  Reset penalty parameter                                   popp 04/10|
 *----------------------------------------------------------------------*/
void CONTACT::Beam3cmanager::ResetCurrentpp()
{
  // get initial value from input file and reset
  currentpp_ = GeneralContactParameters().get<double>("PENALTYPARAM");
  
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
 |  Reset pairs vector                                     mueller 11/11|
 *----------------------------------------------------------------------*/
void CONTACT::Beam3cmanager::ResetPairs()
{
  if (!newgapfunction_)
  {
    // Reset pairs vector to size zero
    pairs_.clear();
  }
  else
  {
    dserror("Are you sure you want to delete the pairs_ vector? The newgapfunction_ "
            "flag is on: In this case history variables of the last time step are needed"
            "and the pairs_ vector should not be deleted!");
  }

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
      std::cout << "\nActive contact set--------------------------------------------------------------\n";
    Comm().Barrier();
      
    // loop over all pairs
    for (int i=0;i<(int)pairs_.size();++i)
    {
      // check if this pair is active
      if (pairs_[i]->GetContactFlag())
      {
        // get coordinates of contact point of each element
        Epetra_SerialDenseVector x1 = pairs_[i]->GetX1();
        Epetra_SerialDenseVector x2 = pairs_[i]->GetX2();

        // make sure to print each pair only once
        // (TODO: this is not yet enough...)
        int firsteleid = (pairs_[i]->Element1())->Id();
        bool firstisinrowmap = RowElements()->MyGID(firsteleid);

        // abbreviations
        int id1 = (pairs_[i]->Element1())->Id();
        int id2 = (pairs_[i]->Element2())->Id();
        double gap = pairs_[i]->GetGap();
        double lm = pairs_[i]->Getlmuzawa() - currentpp_ * pairs_[i]->GetGap();

        // print some output (use printf-method for formatted output)
        if (firstisinrowmap)
        {
          printf("ACTIVE PAIR: %d & %d \t gap: %e \t lm: %e \n",id1,id2,gap,lm);
          //printf("x1 = %e %e %e \t x2 = %e %e %e \n",x1[0],x1[1],x1[2],x2[0],x2[1],x2[2]);
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
  // we need to address the nodes / dofs via the beam contact
  // discretization, because only this is exported to full overlap
  Epetra_Vector fintbc(fint);
  fintbc.ReplaceMap(*ContactDiscret().DofRowMap());
  Epetra_Vector dirichtogglebc(dirichtoggle);
  dirichtogglebc.ReplaceMap(*ContactDiscret().DofRowMap());

  // compute bearing reactions from fint via dirichtoggle
  // Note: dirichtoggle is 1 for DOFs with DBC and 0 elsewise
  Epetra_Vector fbearing(*ContactDiscret().DofRowMap());
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
  const DRT::Node* thisnode = ContactDiscret().gNode(i);
  const std::vector<int> DofGIDs = ContactDiscret().Dof(thisnode);
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
        
  double eleradius = 0;

  // get radius of element
  const DRT::ElementType & eot = thisele->ElementType();

  if ( eot == DRT::ELEMENTS::Beam3Type::Instance() )
  {
    const DRT::ELEMENTS::Beam3* thisbeam = static_cast<const DRT::ELEMENTS::Beam3*>(thisele);
    eleradius = MANIPULATERADIUSVIS*sqrt(sqrt(4 * (thisbeam->Izz()) / M_PI));
  }
  if ( eot == DRT::ELEMENTS::Beam3iiType::Instance() )
  {
    const DRT::ELEMENTS::Beam3ii* thisbeam = static_cast<const DRT::ELEMENTS::Beam3ii*>(thisele);
    eleradius = MANIPULATERADIUSVIS*sqrt(sqrt(4 * (thisbeam->Izz()) / M_PI));
  }
  if ( eot == DRT::ELEMENTS::Beam3ebType::Instance() )
  {
    const DRT::ELEMENTS::Beam3eb* thisbeam = static_cast<const DRT::ELEMENTS::Beam3eb*>(thisele);
    eleradius = MANIPULATERADIUSVIS*sqrt(sqrt(4 * (thisbeam->Izz()) / M_PI));
  }
  
  // declaring variable for color of elements
  double color = 1.0;   
  for (int i=0;i<(int)pairs_.size();++i)
  {
    // abbreviations
    int id1 = (pairs_[i]->Element1())->Id();
    int id2 = (pairs_[i]->Element2())->Id();
    bool active = pairs_[i]->GetContactFlag();
    
    // if element is memeber of an active contact pair, choose different color
    if ( (thisele->Id()==id1 || thisele->Id()==id2) && active) color = 0.875;
  }
    
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
    
  double eleradius = 0;

  // get radius of element
  const DRT::ElementType & eot = thisele->ElementType();

    if ( eot == DRT::ELEMENTS::Beam3Type::Instance() )
      {
        const DRT::ELEMENTS::Beam3* thisbeam = static_cast<const DRT::ELEMENTS::Beam3*>(thisele);
        eleradius = MANIPULATERADIUSVIS*sqrt(sqrt(4 * (thisbeam->Izz()) / M_PI));
      }
    if ( eot == DRT::ELEMENTS::Beam3iiType::Instance() )
      {
        const DRT::ELEMENTS::Beam3ii* thisbeam = static_cast<const DRT::ELEMENTS::Beam3ii*>(thisele);
        eleradius = MANIPULATERADIUSVIS*sqrt(sqrt(4 * (thisbeam->Izz()) / M_PI));
      }

  // declaring variable for color of elements
  double color = 1.0;   
  for (int i=0;i<(int)pairs_.size();++i)
  {
    // abbreviations
    int id1 = (pairs_[i]->Element1())->Id();
    int id2 = (pairs_[i]->Element2())->Id();
    bool active = pairs_[i]->GetContactFlag();
    
    // if element is memeber of an active contact pair, choose different color
    if ( (thisele->Id()==id1 || thisele->Id()==id2) && active) color = 0.0; 
  }
  
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
 |  Compute gmsh output for N-noded elements (private)       meier 02/14|
 *----------------------------------------------------------------------*/
void CONTACT::Beam3cmanager::GMSH_N_noded(const int& n,
                                          const int& n_axial,
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

  // declaring variable for color of elements
  double color = 1.0;
  for (int i=0;i<(int)pairs_.size();++i)
  {
    // abbreviations
    int id1 = (pairs_[i]->Element1())->Id();
    int id2 = (pairs_[i]->Element2())->Id();
    bool active = pairs_[i]->GetContactFlag();

    // if element is memeber of an active contact pair, choose different color
    if ( (thisele->Id()==id1 || thisele->Id()==id2) && active) color = 0.0;
  }

  // computation of coordinates starts here
  for (int i=0;i<n_axial-1;++i)
  {
    // prisms between node i and node i+1
    for (int j=0;j<3;++j)
    {
      coord(j,0) = allcoord(j,i);
      coord(j,1) = allcoord(j,i+1);
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
 |  Determine number of nodes and nodal DoFs of element      meier 02/14|
 *----------------------------------------------------------------------*/
void CONTACT::Beam3cmanager::SetElementTypeAndDistype(DRT::Element* ele1)
{
  numnodes_ = ele1->NumNode();

  const DRT::ElementType & ele1_type = ele1->ElementType();
  if (ele1_type ==DRT::ELEMENTS::Beam3Type::Instance() or ele1_type ==DRT::ELEMENTS::Beam3iiType::Instance())
  {
    numnodalvalues_ = 1;
  }
  else if (ele1_type ==DRT::ELEMENTS::Beam3ebType::Instance())
  {
    numnodalvalues_ = 2;
  }
  else
  {
    dserror("Element type not valid: only beam3, beam3ii and beam3eb is possible for beam contact!");
  }

  return;
}

/*----------------------------------------------------------------------*
 |  Check, if pair with given Id is already in pairs_       meier 02/14|
 *----------------------------------------------------------------------*/
bool CONTACT::Beam3cmanager::PairAllreadyExisting(int currid1, int currid2)
{
  bool pair_allready_existing = false;
  for (unsigned int i=0;i<pairs_.size();i++)
  {
    int id1 = (pairs_[i]->Element1())->Id();
    int id2 = (pairs_[i]->Element2())->Id();
    if ((id1==currid1 and id2==currid2) or (id1==currid2 and id2==currid1))
    {
      pair_allready_existing = true;
      break;
    }
  }

  return pair_allready_existing;
}
