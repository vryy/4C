/*!-------------------------------------------------------------------------*\
 * \file drt_meshfree_discret.cpp
 *
 * \brief discretisation with additional reference point vector
 *        (meshfree analysis)
 *
 * <pre>
 * Maintainer: Keijo Nissen (nis)
 *             nissen@lnm.mw.tum.de
 *             http://www.lnm.mw.tum.de
 *             089 - 289-15253
 * </pre>
 *
\*--------------------------------------------------------------------------*/

#include "../drt_lib/drt_globalproblem.H"
#include "drt_meshfree_discret.H"
#include "drt_meshfree_node.H"
#include "drt_meshfree_cell.H"
#include "../drt_fem_general/drt_utils_maxent_basisfunctions.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_utils_factory.H"
#include "../linalg/linalg_solver.H"
#include "../drt_io/io_control.H"

/*--------------------------------------------------------------------------*
 |  ctor                                                 (public) nis Oct11 |
 *--------------------------------------------------------------------------*/
DRT::MESHFREE::MeshfreeDiscretization::MeshfreeDiscretization(
  const std::string         name,
  Teuchos::RCP<Epetra_Comm> comm,
  const Teuchos::ParameterList & params)
  : DRT::Discretization::Discretization(name,comm),
    assigned_(false)
{
  // neighbourhood search
  nodeassigntype_ = DRT::INPUT::IntegralValue<INPAR::MESHFREE::NodeAssignType>(params,"NODEKNOTASSIGNMENT");

  // type of meshfree approximation
  INPAR::MESHFREE::meshfreetype meshfreetype = DRT::INPUT::IntegralValue<INPAR::MESHFREE::meshfreetype>(params,"TYPE");
  switch (meshfreetype)
  {
  case INPAR::MESHFREE::maxent:
  {
    solutionapprox_  = Teuchos::rcp(new MaxEntApprox(params));
    weightingapprox_ = Teuchos::null;
    break;
  }
  default: dserror("No valid meshfree discretization.");
  }

  // create solver for DirichletBC with non-interpolatory basis functions if necessary
  dbcsolver_ = Teuchos::null;
  const int dbcsolvernumber = params.get<int>("DBC_SOLVER");
  if (dbcsolvernumber != (-1))
  {
    dbcsolver_ = Teuchos::rcp(new LINALG::Solver(DRT::Problem::Instance()->SolverParams(dbcsolvernumber),
                                                 this->Comm(),
                                                 DRT::Problem::Instance()->ErrorFile()->Handle()));
    //this->ComputeNullSpaceIfNecessary(dbcsolver_->Params());
  }

  return;
}

/*--------------------------------------------------------------------------*
 |  dtor                                                 (public) nis Oct11 |
 *--------------------------------------------------------------------------*/
DRT::MESHFREE::MeshfreeDiscretization::~MeshfreeDiscretization()
{
  return;
}

/*--------------------------------------------------------------------------*
 |  get knot row map                                     (public) nis Oct11 |
 *--------------------------------------------------------------------------*/
const Epetra_Map* DRT::MESHFREE::MeshfreeDiscretization::KnotRowMap() const
{
#ifdef DEBUG
  if (Filled()) return knotrowmap_.get();
  else dserror("FillComplete() must be called before call to KnotRowMap()");
  return NULL;
#else
  return knotrowmap_.get();
#endif
}

/*--------------------------------------------------------------------------*
 |  get knot column map                                  (public) nis Oct11 |
 *--------------------------------------------------------------------------*/
const Epetra_Map* DRT::MESHFREE::MeshfreeDiscretization::KnotColMap() const
{
#ifdef DEBUG
  if (Filled()) return knotcolmap_.get();
  else dserror("FillComplete() must be called before call to KnotColMap()");
  return NULL;
#else
  return knotcolmap_.get();
#endif
}

/*--------------------------------------------------------------------------*
 | get global no of knots                                (public) nis Jan12 |
 *--------------------------------------------------------------------------*/
int DRT::MESHFREE::MeshfreeDiscretization::NumGlobalKnots() const
{
#ifdef DEBUG
  if (Filled()) return KnotRowMap()->NumGlobalElements();
  else dserror("FillComplete() must be called before call to NumGlobalKnots()");
  return -1;
#else
  return KnotRowMap()->NumGlobalElements();
#endif
}

/*--------------------------------------------------------------------------*
 | get no of my row knots                                (public) nis Jan12 |
 *--------------------------------------------------------------------------*/
int DRT::MESHFREE::MeshfreeDiscretization::NumMyRowKnots() const
{
#ifdef DEBUG
  if (Filled()) return KnotRowMap()->NumMyElements();
  else dserror("FillComplete() must be called before call to NumMyRowKnots()");
  return -1;
#else
  return KnotRowMap()->NumMyElements();
#endif
}

/*--------------------------------------------------------------------------*
 | get no of my column knots                             (public) nis Jan12 |
 *--------------------------------------------------------------------------*/
int DRT::MESHFREE::MeshfreeDiscretization::NumMyColKnots() const
{
  if (Filled()) return KnotColMap()->NumMyElements();
  else return (int)knot_.size();
}

/*--------------------------------------------------------------------------*
 | query existance of node                               (public) nis Jan12 |
 *--------------------------------------------------------------------------*/
bool DRT::MESHFREE::MeshfreeDiscretization::HaveGlobalKnot(int gid) const
{
  std::map<int,Teuchos::RCP<DRT::MESHFREE::MeshfreeNode> >:: const_iterator curr = knot_.find(gid);
  if (curr == knot_.end()) return false;
  else                     return true;
}

/*--------------------------------------------------------------------------*
 | get node with global id                               (public) nis Jan12 |
 *--------------------------------------------------------------------------*/
DRT::MESHFREE::MeshfreeNode* DRT::MESHFREE::MeshfreeDiscretization::gKnot(int gid) const
{
#ifdef DEBUG
  std::map<int,Teuchos::RCP<DRT::MESHFREE::MeshfreeNode> >:: const_iterator curr = knot_.find(gid);
  if (curr == knot_.end()) dserror("Knot with global id gid=%d not stored on this proc",gid);
  else                     return curr->second.get();
  return NULL;
#else
  return knot_.find(gid)->second.get();
#endif
}

/*--------------------------------------------------------------------------*
 | add a knot                                            (public) nis Jan12 |
 *--------------------------------------------------------------------------*/
void DRT::MESHFREE::MeshfreeDiscretization::AddKnot(Teuchos::RCP<DRT::MESHFREE::MeshfreeNode> knot)
{
  knot_[knot->Id()] = knot;
  Reset();
  return;
}

/*--------------------------------------------------------------------------*
 | delete an knot                                        (public) nis Jan12 |
 *--------------------------------------------------------------------------*/
bool DRT::MESHFREE::MeshfreeDiscretization::DeleteKnot(Teuchos::RCP<DRT::MESHFREE::MeshfreeNode> knot)
{
  std::map<int,Teuchos::RCP<DRT::MESHFREE::MeshfreeNode> >::iterator fool = knot_.find(knot->Id());
  if (fool==knot_.end()) return false;
  knot_.erase(fool);
  Reset();
  return true;
}

/*--------------------------------------------------------------------------*
 | delete an knot                                        (public) nis Jan12 |
 *--------------------------------------------------------------------------*/
bool DRT::MESHFREE::MeshfreeDiscretization::DeleteKnot(const int gid)
{
  std::map<int,Teuchos::RCP<DRT::MESHFREE::MeshfreeNode> >::iterator fool = knot_.find(gid);
  if (fool==knot_.end()) return false;
  knot_.erase(fool);
  Reset();
  return true;
}

/*--------------------------------------------------------------------------*
 |  Pack local knots (row map) into buffer               (public) nis Jan12 |
 *--------------------------------------------------------------------------*/
Teuchos::RCP<std::vector<char> > DRT::MESHFREE::MeshfreeDiscretization::PackMyKnots() const
{
  if (!Filled()) dserror("FillComplete was not called on this discretization");

  DRT::PackBuffer buffer;

  for (std::vector<DRT::MESHFREE::MeshfreeNode*>::const_iterator i=knotrowptr_.begin();
       i!=knotrowptr_.end();
       ++i)
  {
    DRT::MESHFREE::MeshfreeNode * k = *i;
    k->Pack(buffer);
  }

  buffer.StartPacking();

  for (std::vector<DRT::MESHFREE::MeshfreeNode*>::const_iterator i=knotrowptr_.begin();
       i!=knotrowptr_.end();
       ++i)
  {
    DRT::MESHFREE::MeshfreeNode * k = *i;
    k->Pack(buffer);
  }

  Teuchos::RCP<std::vector<char> > block = Teuchos::rcp(new std::vector<char>);
  std::swap( *block, buffer() );
  return block;
}

/*--------------------------------------------------------------------------*
 |  Unpack knot buffer and create local knots            (public) nis Jan12 |
 *--------------------------------------------------------------------------*/
void DRT::MESHFREE::MeshfreeDiscretization::UnPackMyKnots(Teuchos::RCP<std::vector<char> > k)
{
  std::vector<char>::size_type index = 0;
  while (index < k->size())
  {
    std::vector<char> data;
    ParObject::ExtractfromPack(index,*k,data);
    DRT::ParObject* o = DRT::UTILS::Factory(data);
    DRT::MESHFREE::MeshfreeNode* knot = dynamic_cast<DRT::MESHFREE::MeshfreeNode*>(o);
    if (knot == NULL)
    {
      dserror("Failed to build a knot from the knot data");
    }
    knot->SetOwner(comm_->MyPID());
    AddKnot(Teuchos::rcp(knot));
  }
  // in case AddKnot forgets...
  Reset();
}


/*--------------------------------------------------------------------------*
 |  initialise assignment node<>{knot,cell}              (public) nis Jan12 |
 *--------------------------------------------------------------------------*/
bool DRT::MESHFREE::MeshfreeDiscretization::InitAssignSingleNode(double const * const & xn,
                                                                 double const & range,
                                                                 int const & numele,
                                                                 DRT::Element** eles,
                                                                 std::set<int>& elegid,
                                                                 std::set<int>& knotgid,
                                                                 int& init_knot)
{
  // initialisiation of auxiliary pointers (and variable)
  DRT::MESHFREE::MeshfreeNode** knots;    // pointer to vector of pointers of knots
  DRT::MESHFREE::MeshfreeNode*  knotcurr; // pointer to current knot
  DRT::MESHFREE::Cell* elecurr;           // pointer to current element; here: meshfree-cell
  bool init(false);  // boolean, whether an initial knot was found for recursive search

  // loop over all elements parsed to this recursion
  int e;
  for (e=0; e<numele; e++){

    // if not yet having dealt with this element, add to set of elegids and...
    if ((elegid.insert(eles[e]->Id())).second){

      // cast element to meshfree cell and get pointer to its knots
      elecurr = dynamic_cast<DRT::MESHFREE::Cell*>(eles[e]);
      knots = elecurr->Knots();

      // loop over all element knots
      for(int k=0; k< elecurr->NumKnot(); k++){
        knotcurr = knots[k];

        // if not yet having dealt with this knot, add to set of knotgids and...
        if ((knotgid.insert(knotcurr->Id())).second){

          // if knot is in range... TODO ? different range for weighting function
          if (knotcurr->DistToPoint(xn)<range){

            // current element to be erased, otherwise it would not be
            // considered in actual search of 'AssignSingleNode([...])'
            elegid.erase(elecurr->Id());

            // set init true, mark current knot as initial one and break
            init = true;
            init_knot = knotcurr->Id();
            break;

          } // end if dist<range
        } // end if new knot
      } // end for all element knots

      // if we already found an inital knot: break
      if (init) break;

    } // end if new element
  } // end for all elements

  return init;
}

/*--------------------------------------------------------------------------*
 |  recursive assignment node<>{knot,cell}               (public) nis Jan12 |
 *--------------------------------------------------------------------------*/
void DRT::MESHFREE::MeshfreeDiscretization::AssignSingleNode(const double* const & xn,
                                                             const int & nodeid,
                                                             DRT::Node* & node,
                                                             double const & range,
                                                             const int & numele,
                                                             DRT::Element** eles,
                                                             std::set<int>& elegid,
                                                             std::set<int>& knotgid,
                                                             const int & myrank,
                                                             std::map<int, std::map<int,int> > & procmap)
{
  // initialisiation of auxiliary pointers
  // (kept to minimum because of recursive calling)
  DRT::MESHFREE::MeshfreeNode** knots;    // pointer to vector of pointers of knots
  DRT::MESHFREE::MeshfreeNode*  knotcurr; // pointer to current knot
  DRT::MESHFREE::Cell* elecurr;           // pointer to current element; here: meshfree-cell

  // loop over all elements parsed to this recursion
  for (int e=0; e<numele; e++){

    // if not yet having dealt with this element, add to set of elegids and...
    if ((elegid.insert(eles[e]->Id())).second){

      // cast element to meshfree cell and get pointer to its knots
      elecurr = dynamic_cast<DRT::MESHFREE::Cell*>(eles[e]);
      knots = elecurr->Knots();

      // loop over all element knots
      for(int k=0; k<elecurr->NumKnot(); k++){
        knotcurr = knots[k];

        // if not yet having dealt with this knot, add to set of knotgids and...
        if ((knotgid.insert(knotcurr->Id())).second){

          // if knot is in range... TODO ? different range for weighting function
          if (knotcurr->DistToPoint(xn)<range){
            // ... and knot is owned by other proc, update procmap.
            //     (every proc only deals with its own knots)
            if (knotcurr->Owner()!=myrank)
              procmap[knotcurr->Owner()].insert(std::make_pair(nodeid,knotcurr->Id()));
            // ... and proc is owner of knot, start recursion over elements of this knot.
            else
              AssignSingleNode(xn,nodeid,node,range,knotcurr->NumElement(),knotcurr->Elements(),elegid,knotgid,myrank,procmap);

            // add node to knot
            knotcurr->AddNodePtr(nodeid, node);

          } // end if dist<range
        } // end if new knot
      } // end for all element knots

      // assign this element<>node combination (after recursive call)
      elecurr->AddNode(nodeid, node);
      node->AddElementPtr(elecurr);

    } // end if new element
  } // end for all elements
  return;
}

/*--------------------------------------------------------------------------*
 |  assignment node<>{knot,cell} by neighbourhood search (public) nis Jan12 |
 *--------------------------------------------------------------------------*/
void DRT::MESHFREE::MeshfreeDiscretization::AssignNodesToKnotsAndCells()
{
  if (knot_.size()==0)
  {
    assigned_ = true;
    return; // dserror("No knot to assign nodes to in meshfree discret.");
  }

  if (assigned_)
    return; // dserror("Already assigned nodes in meshfree discret. You shouldn't be here.");

  // clear all nodal information in all cells - just to make sure
  Unassign();

  switch (nodeassigntype_)
  {
  case INPAR::MESHFREE::procwise:
  {
    std::set<int> elegid;  // set eof element-gids already handled by this proc
    std::set<int> knotgid; // set eof knot-gids already handled by this proc

    bool init;               // boolean, whether an initial knot was found for recursive search
    int init_knotid;         // id of knot to start recursive search with
    MeshfreeNode* init_knot = NULL; // knot to start recursive search with
    int nodeid;              // id of current node
    DRT::Node* node;         // current node

    const int myrank(this->Comm().MyPID()); // id of this processor
    const int numproc(Comm().NumProc());    // total number of processors
    const double range(solutionapprox_->GetRange()); // range of basis functions

    // procmap: <rank <nodeid,knotid> >
    std::map<int, std::map<int,int> > procmap; // maps to ranks a map that maps to nodes the initial knot

    // on proc assignment
    for (int n=0; n<NumMyColNodes(); n++){
      elegid.clear();
      knotgid.clear();
      init_knotid = -1;

      // get new node and nodeid
      node = noderowptr_[n];
      nodeid = node->Id();

      // get initial knot for recursive call on proc:
      // first by searching on elements of node itself and than on proc
      // discretisation (only for own nodes but not for ghosted ones)
      init = InitAssignSingleNode(node->X(),range,node->NumElement(),node->Elements(),elegid,knotgid,init_knotid);
      if (!init){
        if (node->Owner()==myrank) {
          std::cout << "numproccolele = " << NumMyColElements() << std::endl;
          init = InitAssignSingleNode(node->X(),range,NumMyColElements(),&elecolptr_[0],elegid,knotgid,init_knotid);
        }
        else {
          node_.erase(node->Id());
          continue;
        }
      }
      if (!init)
        dserror("There is no knot within the range of node %d on its proc.",nodeid);

      // get initial knot for recursive search
      init_knot = knot_.find(init_knotid)->second.get();

      if (init_knot->Owner()!=myrank){
        procmap[init_knot->Owner()].insert(std::make_pair(nodeid,init_knot->Id()));
      }

      // assignment node<>{knot,cell} on proc by recursive function call
      node->ClearMyElementTopology();

      AssignSingleNode(node->X(),nodeid,node,range,init_knot->NumElement(),init_knot->Elements(),elegid,knotgid,myrank,procmap);
    }

    // off proc communication, if more than one proc
    if (numproc>1){
      dserror("Parallel node assignment not implemented, yet.");
      std::map<int, std::map<int,int> > testmap;
      for (int i=0; i<numproc; i++){
        testmap.clear();
        // make sure every proc is here
        if (i==myrank){
          // send map of nodes to respective procs
        }
        else{
          // wait for map of nodes from proc i and
          // assign recursively and keep track of nodes
        }
        // when test map exists, you have a problem my friend!
        testmap[i].clear();
        for (int j=0; j<numproc; j++)
        if (!testmap[j].empty()){
          // communicate testmap back to proc i
          // let proc i check if any new node have to communicated
          // if so rerun loop
        }
      } // end proc-loop
    } // end if parallel
    break;
  } // end case nodeassigntype_==procwise
  case INPAR::MESHFREE::blockwise:
  {
    dserror("Blockwise assignment of nodes to elements not yet implemented.");
    break;
  } // end case nodeassigntype_==blockwise
  default:
    dserror("Invalid node assignment type.");
  } // end switch (nodeassigntype_)

  assigned_ = true;

  return;
}

/*--------------------------------------------------------------------------*
 |  << operator for meshfree discretization                        nis Jan12|
 *--------------------------------------------------------------------------*/
std::ostream& operator << (std::ostream& os, const DRT::MESHFREE::MeshfreeDiscretization& mdis)
{
  mdis.Print(os);
  return os;
}

/*--------------------------------------------------------------------------*
 |  Print meshfree discretization (public)                         nis Jan12|
 *--------------------------------------------------------------------------*/
void DRT::MESHFREE::MeshfreeDiscretization::Print(std::ostream& os) const
{
  int numglobalelements = 0;
  int numglobalnodes    = 0;
  int numglobalknots    = 0;
  if (Filled())
  {
    numglobalelements = NumGlobalElements();
    numglobalnodes    = NumGlobalNodes();
    numglobalknots    = NumGlobalKnots();
  }
  else
  {
    int nummyele   = 0;
    std::map<int,Teuchos::RCP<DRT::Element> >::const_iterator ecurr;
    for (ecurr=element_.begin(); ecurr != element_.end(); ++ecurr)
      if (ecurr->second->Owner() == Comm().MyPID()) nummyele++;

    int nummynodes = 0;
    std::map<int,Teuchos::RCP<DRT::Node> >::const_iterator ncurr;
    for (ncurr=node_.begin(); ncurr != node_.end(); ++ncurr)
      if (ncurr->second->Owner() == Comm().MyPID()) nummynodes++;

    int nummyknots = 0;
    std::map<int,Teuchos::RCP<DRT::MESHFREE::MeshfreeNode> >::const_iterator kcurr;
    for (kcurr=knot_.begin(); kcurr != knot_.end(); ++kcurr)
      if (kcurr->second->Owner() == Comm().MyPID()) nummyknots++;

    Comm().SumAll(&nummyele,&numglobalelements,1);
    Comm().SumAll(&nummynodes,&numglobalnodes,1);
    Comm().SumAll(&nummyknots,&numglobalknots,1);
  }

  // print head
  if (Comm().MyPID()==0)
  {
    os << "--------------------------------------------------\n";
    os << "Discretization: " << Name() << std::endl;
    os << "--------------------------------------------------\n";
    os << numglobalelements << " Elements " << numglobalnodes << " Nodes (global) " << numglobalknots << " Knots (global)\n";
    os << "--------------------------------------------------\n";
    if (Filled())
    os << "Filled() = true\n";
    else
    os << "Filled() = false\n";
    os << "--------------------------------------------------\n";
  }
  // print elements
  for (int proc=0; proc < Comm().NumProc(); ++proc)
  {
    if (proc == Comm().MyPID())
    {
      if ((int)element_.size())
        os << "-------------------------- Proc " << proc << " :\n";
      std::map<int,Teuchos::RCP<DRT::Element> >:: const_iterator curr;
      for (curr = element_.begin(); curr != element_.end(); ++curr)
      {
        os << *(curr->second);
        if (Filled())
        {
          std::vector<int> dof = Dof(0,&*(curr->second));
          if (dof.size())
          {
            os << " Dofs ";
            for (unsigned i=0; i<dof.size(); ++i) os << std::setw(6) << dof[i] << " ";
          }
        }
        os << std::endl;
      }
      os << std::endl;
    }
    Comm().Barrier();
  }
  // print nodes
  for (int proc=0; proc < Comm().NumProc(); ++proc)
  {
    if (proc == Comm().MyPID())
    {
      if ((int)node_.size())
        os << "-------------------------- Proc " << proc << " :\n";
      std::map<int,Teuchos::RCP<DRT::Node> >:: const_iterator curr;
      for (curr = node_.begin(); curr != node_.end(); ++curr)
      {
        os << *(curr->second);
        if (Filled())
        {
          std::vector<int> dof = Dof(0,&*(curr->second));
          if (dof.size())
          {
            os << " Dofs ";
            for (unsigned i=0; i<dof.size(); ++i) os << std::setw(6) << dof[i] << " ";
          }
        }
        os << std::endl;
      }
      os << std::endl;
    }
    Comm().Barrier();
  }
  // print knots
  for (int proc=0; proc < Comm().NumProc(); ++proc)
  {
    if (proc == Comm().MyPID())
    {
      if ((int)knot_.size())
        os << "-------------------------- Proc " << proc << " :\n";
      std::map<int,Teuchos::RCP<DRT::MESHFREE::MeshfreeNode> >:: const_iterator curr;
      for (curr = knot_.begin(); curr != knot_.end(); ++curr)
      {
        os << *(curr->second);
        if (Filled())
        {
          std::vector<int> dof = Dof(0,&*(curr->second));
          if (dof.size())
          {
            os << " Dofs ";
            for (unsigned i=0; i<dof.size(); ++i) os << std::setw(6) << dof[i] << " ";
          }
        }
        os << std::endl;
      }
      os << std::endl;
    }
    Comm().Barrier();
  }
  // print conditions
  for (int proc=0; proc<Comm().NumProc(); ++proc)
  {
    if (proc==Comm().MyPID())
    {
      int numcond = condition_.size();
      if (numcond)
        os << "-------------------------- Proc " << proc << " :\n";
      if (numcond)
      {
        os << numcond << " Conditions:\n";
        std::map<std::string,Teuchos::RCP<Condition> >::const_iterator curr;
        for (curr=condition_.begin(); curr != condition_.end(); ++curr)
        {
          os << curr->first << " ";
          os << *(curr->second) << std::endl;
        }
      }
      os << std::endl;
    }
    Comm().Barrier();
  }
  return;
}
