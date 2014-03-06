/*!-------------------------------------------------------------------------*\
 * \file drt_meshfree_discret.cpp
 *
 * \brief discretisation for meshfree analysis
 *
 * <pre>
 * Maintainer: Keijo Nissen (nis)
 *             nissen@lnm.mw.tum.de
 *             http://www.lnm.mw.tum.de
 *             089 - 289-15253
 * </pre>
 *
\*--------------------------------------------------------------------------*/

#include "drt_meshfree_discret.H"
#include "drt_meshfree_utils.H"
#include "drt_meshfree_node.H"
#include "drt_meshfree_cell.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_fem_general/drt_utils_maxent_basisfunctions.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_utils_factory.H"
#include "../linalg/linalg_solver.H"
#include "../linalg/linalg_serialdensevector.H"
#include "../drt_io/io_control.H"

/*--------------------------------------------------------------------------*
 |  ctor                                                 (public) nis Oct11 |
 *--------------------------------------------------------------------------*/
DRT::MESHFREE::MeshfreeDiscretization::MeshfreeDiscretization(
  const std::string         name,
  Teuchos::RCP<Epetra_Comm> comm,
  const Teuchos::ParameterList & params)
  : DRT::Discretization::Discretization(name,comm),
    assigned_(false),
    domaintopology_(4)
{
  // neighbourhood search
  nodeassigntype_ = DRT::INPUT::IntegralValue<INPAR::MESHFREE::NodeAssignType>(params,"NODEPOINTASSIGNMENT");

  // type of meshfree approximation
  INPAR::MESHFREE::meshfreetype meshfreetype = DRT::INPUT::IntegralValue<INPAR::MESHFREE::meshfreetype>(params,"TYPE");
  switch (meshfreetype)
  {
  case INPAR::MESHFREE::maxent:
  {
    solutionapprox_  = Teuchos::rcp(new MaxEntApprox(params,0));
    weightingapprox_  = Teuchos::rcp(new MaxEntApprox(params,1));
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
 |  get point row map                                     (public) nis Oct11 |
 *--------------------------------------------------------------------------*/
const Epetra_Map* DRT::MESHFREE::MeshfreeDiscretization::PointRowMap() const
{
#ifdef DEBUG
  if (Filled()) return pointrowmap_.get();
  else dserror("FillComplete() must be called before call to PointRowMap()");
  return NULL;
#else
  return pointrowmap_.get();
#endif
}

/*--------------------------------------------------------------------------*
 |  get point column map                                  (public) nis Oct11 |
 *--------------------------------------------------------------------------*/
const Epetra_Map* DRT::MESHFREE::MeshfreeDiscretization::PointColMap() const
{
#ifdef DEBUG
  if (Filled()) return pointcolmap_.get();
  else dserror("FillComplete() must be called before call to PointColMap()");
  return NULL;
#else
  return pointcolmap_.get();
#endif
}

/*--------------------------------------------------------------------------*
 | get global no of points                                (public) nis Jan12 |
 *--------------------------------------------------------------------------*/
int DRT::MESHFREE::MeshfreeDiscretization::NumGlobalPoints() const
{
#ifdef DEBUG
  if (Filled()) return PointRowMap()->NumGlobalElements();
  else dserror("FillComplete() must be called before call to NumGlobalPoints()");
  return -1;
#else
  return PointRowMap()->NumGlobalElements();
#endif
}

/*--------------------------------------------------------------------------*
 | get no of my row points                                (public) nis Jan12 |
 *--------------------------------------------------------------------------*/
int DRT::MESHFREE::MeshfreeDiscretization::NumMyRowPoints() const
{
#ifdef DEBUG
  if (Filled()) return PointRowMap()->NumMyElements();
  else dserror("FillComplete() must be called before call to NumMyRowPoints()");
  return -1;
#else
  return PointRowMap()->NumMyElements();
#endif
}

/*--------------------------------------------------------------------------*
 | get no of my column points                             (public) nis Jan12 |
 *--------------------------------------------------------------------------*/
int DRT::MESHFREE::MeshfreeDiscretization::NumMyColPoints() const
{
  if (Filled()) return PointColMap()->NumMyElements();
  else return (int)point_.size();
}

/*--------------------------------------------------------------------------*
 | query existance of node                               (public) nis Jan12 |
 *--------------------------------------------------------------------------*/
bool DRT::MESHFREE::MeshfreeDiscretization::HaveGlobalPoint(int gid) const
{
  std::map<int,Teuchos::RCP<DRT::MESHFREE::MeshfreeNode> >:: const_iterator curr = point_.find(gid);
  if (curr == point_.end()) return false;
  else                     return true;
}

/*--------------------------------------------------------------------------*
 | get node with global id                               (public) nis Jan12 |
 *--------------------------------------------------------------------------*/
DRT::MESHFREE::MeshfreeNode* DRT::MESHFREE::MeshfreeDiscretization::gPoint(int gid) const
{
#ifdef DEBUG
  std::map<int,Teuchos::RCP<DRT::MESHFREE::MeshfreeNode> >:: const_iterator curr = point_.find(gid);
  if (curr == point_.end()) dserror("Point with global id gid=%d not stored on this proc",gid);
  else                     return curr->second.get();
  return NULL;
#else
  return point_.find(gid)->second.get();
#endif
}

/*--------------------------------------------------------------------------*
 | add a point                                            (public) nis Jan12 |
 *--------------------------------------------------------------------------*/
void DRT::MESHFREE::MeshfreeDiscretization::AddPoint(Teuchos::RCP<DRT::MESHFREE::MeshfreeNode> point)
{
  point_[point->Id()] = point;
  Reset();
  return;
}

/*--------------------------------------------------------------------------*
 | delete an point                                        (public) nis Jan12 |
 *--------------------------------------------------------------------------*/
bool DRT::MESHFREE::MeshfreeDiscretization::DeletePoint(Teuchos::RCP<DRT::MESHFREE::MeshfreeNode> point)
{
  std::map<int,Teuchos::RCP<DRT::MESHFREE::MeshfreeNode> >::iterator fool = point_.find(point->Id());
  if (fool==point_.end()) return false;
  point_.erase(fool);
  Reset();
  return true;
}

/*--------------------------------------------------------------------------*
 | delete an point                                        (public) nis Jan12 |
 *--------------------------------------------------------------------------*/
bool DRT::MESHFREE::MeshfreeDiscretization::DeletePoint(const int gid)
{
  std::map<int,Teuchos::RCP<DRT::MESHFREE::MeshfreeNode> >::iterator fool = point_.find(gid);
  if (fool==point_.end()) return false;
  point_.erase(fool);
  Reset();
  return true;
}

/*--------------------------------------------------------------------------*
 |  Pack local points (row map) into buffer               (public) nis Jan12 |
 *--------------------------------------------------------------------------*/
Teuchos::RCP<std::vector<char> > DRT::MESHFREE::MeshfreeDiscretization::PackMyPoints() const
{
  if (!Filled()) dserror("FillComplete was not called on this discretization");

  DRT::PackBuffer buffer;

  for (std::vector<DRT::MESHFREE::MeshfreeNode*>::const_iterator i=pointrowptr_.begin();
       i!=pointrowptr_.end();
       ++i)
  {
    DRT::MESHFREE::MeshfreeNode * k = *i;
    k->Pack(buffer);
  }

  buffer.StartPacking();

  for (std::vector<DRT::MESHFREE::MeshfreeNode*>::const_iterator i=pointrowptr_.begin();
       i!=pointrowptr_.end();
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
 |  Unpack point buffer and create local points            (public) nis Jan12 |
 *--------------------------------------------------------------------------*/
void DRT::MESHFREE::MeshfreeDiscretization::UnPackMyPoints(Teuchos::RCP<std::vector<char> > k)
{
  std::vector<char>::size_type index = 0;
  while (index < k->size())
  {
    std::vector<char> data;
    ParObject::ExtractfromPack(index,*k,data);
    DRT::ParObject* o = DRT::UTILS::Factory(data);
    DRT::MESHFREE::MeshfreeNode* point = dynamic_cast<DRT::MESHFREE::MeshfreeNode*>(o);
    if (point == NULL)
    {
      dserror("Failed to build a point from the point data");
    }
    point->SetOwner(comm_->MyPID());
    AddPoint(Teuchos::rcp(point));
  }
  // in case AddPoint forgets...
  Reset();
}

/*----------------------------------------------------------------------------*
 |  Adds node set topology to faces (public)                        nis Jan14 |
 *----------------------------------------------------------------------------*/
void DRT::MESHFREE::MeshfreeDiscretization::AddNodeSetTopology(
  std::vector<std::vector<std::vector<int> >* > nodeset
  )
{
  int dim = DRT::Problem::Instance()->NDim();
  for (int i=(dim+1); i<=3; ++i)
    if (nodeset[i]->size()>0)
      dserror("Can't assign %i-D face in %i-D problem.",i,dim);
  if (nodeset[dim]->size()!=1)
    dserror("There has to be one and only one %i-D face in a %i-D problem.",dim,dim);

  for (int i=0; i<4; ++i)
  {
    const unsigned nface = nodeset[i]->size();
    domaintopology_[i].resize(nface);
    for (unsigned j=0; j<nface; ++j)
      domaintopology_[i][j] = Face((*(nodeset[i]))[j],i,this);
  }

  return;
}

/*--------------------------------------------------------------------------*
 |  initialise assignment node<>{point,cell}              (public) nis Jan12 |
 *--------------------------------------------------------------------------*/
bool DRT::MESHFREE::MeshfreeDiscretization::InitAssignSingleNode(const double * const xn,
                                                                 const double range,
                                                                 const int numcell,
                                                                 DRT::Element** cells,
                                                                 std::set<int>& cellgid,
                                                                 std::set<int>& pointgid,
                                                                 int& init_point)
{
  // initialisiation of auxiliary variables
  bool init(false);  // boolean, whether an initial point was found for recursive search

  // loop over all cell parsed to this recursion
  for (int c=0; c<numcell; c++)
  {
    // get pointer to current cell
    DRT::MESHFREE::Cell* cellcurr = dynamic_cast<DRT::MESHFREE::Cell*>(cells[c]);

    // if not yet having dealt with this cell, add to set of cellgids and...
    if ((cellgid.insert(cellcurr->Id())).second)
    {
      // get pointer to pointer to cell points
      DRT::Node** points = cellcurr->Points();

      // loop over all cell points
      for(int p=0; p< cellcurr->NumPoint(); p++)
      {
        // get pointer to current cell points
        DRT::MESHFREE::MeshfreeNode* pointcurr = dynamic_cast<DRT::MESHFREE::MeshfreeNode*>(points[p]);

        // if not yet having dealt with this point, add to set of pointgids and...
        if ((pointgid.insert(pointcurr->Id())).second)
        {
          // if point is in range... TODO ? different range for weighting function
          if (pointcurr->DistToPoint(xn)<range)
          {
            // current cell to be erased, otherwise it would not be
            // considered in actual search of 'AssignSingleNode([...])'
            cellgid.erase(cellcurr->Id());

            // set init true, mark current point as initial one and break
            init = true;
            init_point = pointcurr->Id();
            break;

          } // end if dist<range
        } // end if new point
      } // end for all cell points

      // if we already found an inital point: break
      if (init) break;

    } // end if new cell

  } // end for all cells

  return init;
}

/*--------------------------------------------------------------------------*
 |  recursive assignment node<>{point,cell}               (public) nis Jan12 |
 *--------------------------------------------------------------------------*/
void DRT::MESHFREE::MeshfreeDiscretization::AssignSingleNode(const double* const xn,
                                                             const int nodeid,
                                                             DRT::Node* node,
                                                             const double range,
                                                             const int numcell,
                                                             DRT::Element** cells,
                                                             std::set<int>& cellgid,
                                                             std::set<int>& pointgid,
                                                             const int myrank,
                                                             std::map<int, std::map<int,int> > & procmap)
{
  // loop over all cellss parsed to this recursion
  for (int c=0; c<numcell; c++)
  {
    // get pointer to current cell
    DRT::MESHFREE::Cell* cellcurr = dynamic_cast<DRT::MESHFREE::Cell*>(cells[c]);

    // if not yet having dealt with this cell, add to set of cellgids and...
    if ((cellgid.insert(cellcurr->Id())).second)
    {
      // get pointer to pointer to cell points
      DRT::Node** points = cellcurr->Points();

      // loop over all cell points
      for(int p=0; p<cellcurr->NumPoint(); p++)
      {
        // get pointer to current cell points
        DRT::MESHFREE::MeshfreeNode* pointcurr = dynamic_cast<DRT::MESHFREE::MeshfreeNode*>(points[p]);

        // if not yet having dealt with this point, add to set of pointgids and...
        if ((pointgid.insert(pointcurr->Id())).second)
        {
         // if point is in range... TODO ? different range for weighting function
          if (pointcurr->DistToPoint(xn)<range)
          {
            // ... and point is owned by other proc, update procmap.
            //     (every proc only deals with its own points)
            if (pointcurr->Owner()!=myrank)
              procmap[pointcurr->Owner()].insert(std::make_pair(nodeid,pointcurr->Id()));
            // ... and proc is owner of point, start recursion over cells of this point.
            else
              AssignSingleNode(xn,nodeid,node,range,pointcurr->NumElement(),pointcurr->Elements(),cellgid,pointgid,myrank,procmap);

            // finally add node to point
            pointcurr->AddNodePtr(nodeid, node);

          } // end if dist<range
        } // end if new point
      } // end for all cell points

      // assign this cell<>node combination (after recursive call)
      cellcurr->AddNode(nodeid, node);
      node->AddElementPtr(cellcurr);

    } // end if new cell
  } // end for all cells

  return;
}

/*--------------------------------------------------------------------------*
 |  assignment node<>{point,cell} by neighbourhood search (public) nis Jan12 |
 *--------------------------------------------------------------------------*/
void DRT::MESHFREE::MeshfreeDiscretization::AssignNodesToPointsAndCells()
{
  if (point_.size()==0)
  {
    assigned_ = true;
    return; // dserror("No point to assign nodes to in meshfree discret.");
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
    std::set<int> pointgid; // set eof point-gids already handled by this proc

    bool init;               // boolean, whether an initial point was found for recursive search
    int init_pointid;         // id of point to start recursive search with
    MeshfreeNode* init_point = NULL; // point to start recursive search with
    int nodeid;              // id of current node
    DRT::Node* node;         // current node

    const int myrank(this->Comm().MyPID()); // id of this processor
    const int numproc(Comm().NumProc());    // total number of processors
    const double range(solutionapprox_->GetRange()); // range of basis functions

    // procmap: <rank <nodeid,pointid> >
    std::map<int, std::map<int,int> > procmap; // maps to ranks a map that maps to nodes the initial point

    // on proc assignment
    for (int n=0; n<NumMyColNodes(); n++){
      elegid.clear();
      pointgid.clear();
      init_pointid = -1;

      // get new node and nodeid
      node = noderowptr_[n];
      nodeid = node->Id();

      // get initial point for recursive call on proc:
      // first by searching on elements of node itself and than on proc
      // discretisation (only for own nodes but not for ghosted ones)
      init = InitAssignSingleNode(node->X(),range,node->NumElement(),node->Elements(),elegid,pointgid,init_pointid);
      if (!init){
        if (node->Owner()==myrank)
          init = InitAssignSingleNode(node->X(),range,NumMyColElements(),elecolptr_.data(),elegid,pointgid,init_pointid);
        else
        {
          node_.erase(node->Id());
          continue;
        }
      }
      if (!init)
        dserror("There is no point within the range of node %d on its proc.",nodeid);

      // get initial point for recursive search
      init_point = point_.find(init_pointid)->second.get();

      if (init_point->Owner()!=myrank){
        procmap[init_point->Owner()].insert(std::make_pair(nodeid,init_point->Id()));
      }

      // assignment node<>{point,cell} on proc by recursive function call
      node->ClearMyElementTopology();
      AssignSingleNode(node->X(),nodeid,node,range,init_point->NumElement(),init_point->Elements(),elegid,pointgid,myrank,procmap);
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
 |  Computes values at nodes from nodal values for non-interpolatory        |
 |  basis function                                                nis Jan14 |
 *--------------------------------------------------------------------------*/
void DRT::MESHFREE::MeshfreeDiscretization::ComputeValuesAtNodes(
  const Teuchos::RCP<const Epetra_Vector>& nodalvalues,
  const Teuchos::RCP<Epetra_Vector>      & valuesatnodes) const
{
  // check if (meshfree!) filled - otherwise no correct node-to-node assignments
  if (!Filled())
    dserror("Fill complete meshfree discretisation before computing values at nodes!");

  // check if we have only single dofset in meshfree discretisation
  if (dofsets_.size()!=1)
    dserror("Up to now values at nodes can only be computes for single dofset meshfree discretisations!");

  // get number of dofs per node from first row node
  const int nummyrownode = NumMyRowNodes();
  if (nummyrownode==0) dserror("can't be!");//return;

  // get number of dofs per node - prerequisite: equal for all nodes
  const int numdofpernode = NumDof(noderowptr_[0]);

  // check if we only have nodal dofs
  if (DofRowMap()->NumGlobalElements()!=numdofpernode*NumGlobalNodes())
    dserror("Either non-nodal dofs exist or number of dofs is not the same for all nodes!");

  // check if nodal dofs are guaranteed te be in a stride
//  if (DRT::Problem::Instance()->BandWidthOpt())
//    dserror("We rely on nodal dofs being in a stride. If this is still guaranted by your BandWidthOpt, remove dserror!");

  // loop over all row nodes
  for (int i=0; i<nummyrownode; ++i)
  {
    //------------------------------------------------------------------------
    // get neighbourhood information
    //------------------------------------------------------------------------
    const MeshfreeNode* const cnode = dynamic_cast<MeshfreeNode*>(noderowptr_[i]);
    const int nneighbour = cnode->NumNode();
    const DRT::Node* const * neighbours = cnode->Nodes();
    const int facedim = cnode->GetFaceDim();

    if (facedim==0)
    {
      //------------------------------------------------------------------------
      // node is a vertex and thus interpolatory: copy nodal values
      //------------------------------------------------------------------------
      std::vector<int>    dofids(numdofpernode);
      std::vector<double> valuesatnode(numdofpernode);
      // loop over all dofs of this node
      for (int idof=0; idof<numdofpernode; ++idof)
      {
        dofids[idof] = Dof(cnode,idof);
        valuesatnode[idof] = (*nodalvalues)[DofRowMap()->LID(Dof(cnode,idof))];
      }

      // write values at nodes into vector
      valuesatnodes->ReplaceGlobalValues(numdofpernode,valuesatnode.data(),dofids.data());
    }
    else
    {
      //------------------------------------------------------------------------
      // construct matrix with distance to neighbours for basis function evaluation
      //------------------------------------------------------------------------
      LINALG::SerialDenseMatrix distnn(3,nneighbour,false);
      const double * cx = cnode->X();
      // loop over all neighbours
      for (int n=0; n<nneighbour; ++n)
      {
        // get position of current neighbour
        const double * xn = neighbours[n]->X();
        // compute distance and write in respective matrix column
        double * ccol = distnn[n];
        for (int k=0; k<3; ++k)
          ccol[k] = cx[k] - xn[k];
      }

      //------------------------------------------------------------------------
      // get solution basis functions of neighbourhood
      //------------------------------------------------------------------------
      Teuchos::RCP<LINALG::SerialDenseVector> funct = Teuchos::rcp(new LINALG::SerialDenseVector(nneighbour,false));
      int err = solutionapprox_->GetMeshfreeBasisFunction(facedim, Teuchos::rcpFromRef(distnn), funct);
      if (err>0)
      {
        std::cout << "When computing the solution basis functions at node " << i << ":" << std::endl;
        DRT::MESHFREE::OutputMeshfreeError(err);
      }

      //------------------------------------------------------------------------
      // compute values at current node
      //------------------------------------------------------------------------
      std::vector<int>    dofids(numdofpernode);
      std::vector<double> valuesatnode(numdofpernode);
      // loop over all dofs of this node
      for (int idof=0; idof<numdofpernode; ++idof)
      {
        dofids[idof] = Dof(cnode,idof);
        double valueatnode = 0.0;
        // loop over all neighbours
        for (int n=0; n<nneighbour; ++n)
          // u(x_i) = sum N(x_i) u_i
          valueatnode += (*funct)[n]*(*nodalvalues)[DofRowMap()->LID(Dof(neighbours[n],idof))];
        valuesatnode[idof] = valueatnode;
      }

      // write values at nodes into vector
      valuesatnodes->ReplaceGlobalValues(numdofpernode,valuesatnode.data(),dofids.data());
    }

  }

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
  int numglobalpoints    = 0;
  if (Filled())
  {
    numglobalelements = NumGlobalElements();
    numglobalnodes    = NumGlobalNodes();
    numglobalpoints   = NumGlobalPoints();
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

    int nummypoints = 0;
    std::map<int,Teuchos::RCP<DRT::MESHFREE::MeshfreeNode> >::const_iterator pcurr;
    for (pcurr=point_.begin(); pcurr != point_.end(); ++pcurr)
      if (pcurr->second->Owner() == Comm().MyPID()) nummypoints++;

    Comm().SumAll(&nummyele,&numglobalelements,1);
    Comm().SumAll(&nummynodes,&numglobalnodes,1);
    Comm().SumAll(&nummypoints,&numglobalpoints,1);
  }

  // print head
  if (Comm().MyPID()==0)
  {
    os << "--------------------------------------------------\n";
    os << "Discretization: " << Name() << std::endl;
    os << "--------------------------------------------------\n";
    os << numglobalelements << " Elements " << numglobalnodes << " Nodes (global) " << numglobalpoints << " Points (global)\n";
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
  // print points
  for (int proc=0; proc < Comm().NumProc(); ++proc)
  {
    if (proc == Comm().MyPID())
    {
      if ((int)point_.size())
        os << "-------------------------- Proc " << proc << " :\n";
      std::map<int,Teuchos::RCP<DRT::MESHFREE::MeshfreeNode> >:: const_iterator curr;
      for (curr = point_.begin(); curr != point_.end(); ++curr)
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
