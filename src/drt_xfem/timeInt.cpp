/*!------------------------------------------------------------------------------------------------*
\file timeInt.cpp

\brief provides the basic time integration classes "TimeInt", "TimeIntStd", "TimeIntEnr"

\level 2

<pre>
\maintainer Benedikt Schott
            schott@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>

\warning this combustion module related file will be deleted within the next time!!!
 *------------------------------------------------------------------------------------------------*/


#include "timeInt.H"
#include "../drt_lib/drt_exporter.H"
#include "../drt_combust/combust_flamefront.H"
#include "../drt_geometry/position_array.H"


/*------------------------------------------------------------------------------------------------*
 * basic XFEM time-integration constructor                                       winklmaier 11/11 *
 *------------------------------------------------------------------------------------------------*/
XFEM::TIMEINT::TIMEINT(
    const Teuchos::RCP<DRT::Discretization> discret,
    const Teuchos::RCP<DofManager> olddofman,
    const Teuchos::RCP<DofManager> newdofman,
    std::vector<Teuchos::RCP<Epetra_Vector> > oldVectors,
    const Teuchos::RCP<COMBUST::FlameFront> flamefront,
    const Epetra_Map& olddofcolmap,
    const Epetra_Map& newdofrowmap,
    const std::map<DofKey, DofGID>& oldNodalDofColDistrib,
    const std::map<DofKey, DofGID>& newNodalDofRowDistrib,
    const Teuchos::RCP<std::map<int,std::vector<int> > > pbcmap
) :
discret_(discret),
olddofman_(olddofman),
newdofman_(newdofman),
olddofcolmap_(olddofcolmap),
newdofrowmap_(newdofrowmap),
oldNodalDofColDistrib_(oldNodalDofColDistrib),
newNodalDofRowDistrib_(newNodalDofRowDistrib),
phin_(flamefront->Phin()),
phinp_(flamefront->Phinp()),
oldinterfacehandle_(flamefront->InterfaceHandleOld()),
newinterfacehandle_(flamefront->InterfaceHandle()),
gradphi_(flamefront->GradPhi()),
oldVectors_(oldVectors),
pbcmap_(pbcmap),
myrank_(discret_->Comm().MyPID()),
numproc_(discret_->Comm().NumProc()),
newton_max_iter_(10),
newton_tol_(1.0e-06)
{
  return;
} // end constructor



/*------------------------------------------------------------------------------------------------*
 * set the computation type with help of the iteration counter                   winklmaier 10/11 *
 *------------------------------------------------------------------------------------------------*/
void XFEM::TIMEINT::type(
    int iter,
    int iterMax
)
{
  if (iter==1)                              FGIType_=FRS1FGI1_;
  else if (iterMax==1 or iter%iterMax==1)   FGIType_=FRS1FGINot1_;
  else                                      FGIType_=FRSNot1_;
} // end function type



/*------------------------------------------------------------------------------------------------*
 * assign the Epetra vectors which shall be computed to the                                       *
 * algorithms data structure                                                     winklmaier 10/11 *
 *------------------------------------------------------------------------------------------------*/
void XFEM::TIMEINT::handleVectors(
    std::vector<Teuchos::RCP<Epetra_Vector> > newRowVectorsn,
    std::vector<Teuchos::RCP<Epetra_Vector> > newRowVectorsnp
)
{
  if (newRowVectorsn.size()!=newRowVectorsnp.size())
    dserror("Number of state-vectors of different times are different!");

  newVectors_ = newRowVectorsn;
  newVectors_.insert(newVectors_.end(),newRowVectorsnp.begin(),newRowVectorsnp.end());

  if (oldVectors_.size() != newVectors_.size())
    dserror("Number of state-vectors at new and old discretization are different!");
} // end function handleVectors



/*------------------------------------------------------------------------------------------------*
 * identify intersection status of an element                                    winklmaier 01/12 *
 *------------------------------------------------------------------------------------------------*/
XFEM::TIMEINT::intersectionType XFEM::TIMEINT::intersectionStatus(
    const DRT::Element* ele,
    bool oldTimeStep
) const
{
  // return status are: 0 = uncut, 1 = bisected/trisected, 2 = numerically intersected

  // initialization
  Teuchos::RCP<COMBUST::InterfaceHandleCombust> ih;
  Teuchos::RCP<Epetra_Vector> phi;

  if (oldTimeStep)
  {
    ih = oldinterfacehandle_;
    phi = phin_;
  }
  else
  {
    ih = newinterfacehandle_;
    phi = phinp_;
  }

  // check if "normally" intersected
  if ((ih->ElementBisected(ele->Id())) or
      (ih->ElementTrisected(ele->Id())))
    return XFEM::TIMEINT::cut_; // bisected or trisected

  // check if numerically intersected:
  // this means, nodes are on different sides of the interface,
  // but the cut status is not bi- or trisected so that some phi-values are nearly zero
  const int* elenodeids = ele->NodeIds();  // nodeids of element
  double phivalue = 0.0;
  bool side = true;
  for (int nodeid=0;nodeid<ele->NumNode();nodeid++) // loop over element nodes
  {
    phivalue = (*phi)[discret_->gNode(elenodeids[nodeid])->LID()];
    if (nodeid == 0)
      side = plusDomain(phivalue);
    else
    {
      if (side!=plusDomain(phivalue)) // different sides
        return XFEM::TIMEINT::numerical_cut_; // numerically intersected
    }
  }

  return XFEM::TIMEINT::uncut_;
} // end function intersectionStatus



/*------------------------------------------------------------------------------------------------*
 * identify interface side of a point in combustion                              winklmaier 10/10 *
 *------------------------------------------------------------------------------------------------*/
int XFEM::TIMEINT::interfaceSide(
    double phi
) const
{
  if (plusDomain(phi) == true) return 1;
  else return -1;
} // end function interfaceSide



/*------------------------------------------------------------------------------------------------*
 * identify interface side of a point in combustion                              winklmaier 10/10 *
 *------------------------------------------------------------------------------------------------*/
int XFEM::TIMEINT::interfaceSide(
    DRT::Element* ele,
    LINALG::Matrix<3,1> x,
    bool newTimeStep
) const
{
  const int nsd = 3; // dimension

  // required interfacehandle
  Teuchos::RCP<COMBUST::InterfaceHandleCombust> ih = newTimeStep ? newinterfacehandle_ : oldinterfacehandle_;
  const GEO::DomainIntCells&  domainIntCells(ih->ElementDomainIntCells(ele->Id()));

  static LINALG::Matrix<nsd,1> eta(true); // local cell coordinates
  bool pointInCell = false; // boolean whether point is in current cell

  // loop over domain integration cells
  for (GEO::DomainIntCells::const_iterator cell = domainIntCells.begin(); cell != domainIntCells.end(); ++cell)
  {
    callXToXiCoords(*cell,x,eta,pointInCell);

    if (pointInCell)
    {
      if (cell->getDomainPlus()==true)        return 1;
      else                                                         return -1;
    }
  }

  // handle special case of tiny cell which was removed by cut:
  // all remaining cells are at same side, values in this ele are all at
  // this side since it is uncut and no more enriched
  bool side = 0;
  for (GEO::DomainIntCells::const_iterator cell = domainIntCells.begin(); cell != domainIntCells.end(); ++cell)
  {
    if (cell == domainIntCells.begin())
      side = cell->getDomainPlus();
    else
      if (side !=cell->getDomainPlus())
        break; // special case of more than once cut element with deleted cell

    if (side==true)              return 1;
    else                            return -1;
  }

  dserror("point in an element should be in one of the elements domain integration cells");
  return 0;
} // end function interfaceSide



/*------------------------------------------------------------------------------------------------*
 * call the computation of local coordinates for an element                      winklmaier 11/11 *
 *------------------------------------------------------------------------------------------------*/
void XFEM::TIMEINT::callXToXiCoords(
    const GEO::DomainIntCell& cell,
    LINALG::Matrix<3,1>& x,
    LINALG::Matrix<3,1>& xi,
    bool& pointInDomain
) const
{
  LINALG::SerialDenseMatrix xyz(cell.CellNodalPosXYZ());
  callXToXiCoords(xyz,cell.Shape(),x,xi,pointInDomain);
} // end function callXToXiCoords



/*------------------------------------------------------------------------------------------------*
 * call the computation of local coordinates for an integration cell             winklmaier 10/10 *
 *------------------------------------------------------------------------------------------------*/
void XFEM::TIMEINT::callXToXiCoords(
    const DRT::Element* ele,
    LINALG::Matrix<3,1>& x,
    LINALG::Matrix<3,1>& xi,
    bool& pointInDomain
) const
{
  LINALG::SerialDenseMatrix xyz(3,ele->NumNode(),true);
  GEO::fillInitialPositionArray(ele, xyz);
  callXToXiCoords(xyz,ele->Shape(),x,xi,pointInDomain);
} // end function callXToXiCoords



/*------------------------------------------------------------------------------------------------*
 * call the computation of local coordinates for a polytop                                        *
 * with corners given by the coordinates                                         winklmaier 10/10 *
 *------------------------------------------------------------------------------------------------*/
void XFEM::TIMEINT::callXToXiCoords(
    LINALG::SerialDenseMatrix& nodecoords,
    DRT::Element::DiscretizationType DISTYPE,
    LINALG::Matrix<3,1>& x,
    LINALG::Matrix<3,1>& xi,
    bool& pointInDomain
) const
{
  switch (DISTYPE)
  {
  case DRT::Element::hex8:
    XToXiCoords<DRT::Element::hex8>(nodecoords,x,xi,pointInDomain);
    break;
  case DRT::Element::hex20:
    XToXiCoords<DRT::Element::hex20>(nodecoords,x,xi,pointInDomain);
    break;
  case DRT::Element::tet4:
    XToXiCoords<DRT::Element::tet4>(nodecoords,x,xi,pointInDomain);
    break;
  default:
  {
    dserror("add your 3D distype and the according transformation!");
    break;
  }
  } // end switch
} // end function callXToXiCoords



/*------------------------------------------------------------------------------------------------*
 * add adjacent elements for a periodic boundary node                            winklmaier 05/11 *
 *------------------------------------------------------------------------------------------------*/
void XFEM::TIMEINT::addPBCelements(
    const DRT::Node* node,
    std::vector<const DRT::Element*>&  eles
) const
{
  const DRT::Element* const* elements = node->Elements(); // element around current node

  for (int iele=0;iele<node->NumElement();iele++)
    eles.push_back(elements[iele]);

  // get pbcnode
  bool pbcnodefound = false; // boolean indicating whether this node is a pbc node
  DRT::Node* pbcnode = NULL;
  findPBCNode(node,pbcnode,pbcnodefound);

  // add elements located around the coupled pbc node
  if (pbcnodefound)
  {
    // get adjacent elements of this node
    const DRT::Element*const* pbcelements = pbcnode->Elements();
    // add elements to list
    for (int iele=0;iele<pbcnode->NumElement();iele++)// = ptToNode->Elements();
    {
      eles.push_back(pbcelements[iele]);
    }
  } // end if pbcnode true
}



void XFEM::TIMEINT::findPBCNode(
    const DRT::Node* node,
    DRT::Node*& pbcnode,
    bool& pbcnodefound
) const
{
  const int nodegid = node->Id(); // global node id
  pbcnodefound = false; // boolean indicating whether this node is a pbc node
  int coupnodegid = -1;
  // loop all nodes with periodic boundary conditions (master nodes)
  for (std::map<int, std::vector<int>  >::const_iterator pbciter= (*pbcmap_).begin(); pbciter != (*pbcmap_).end(); ++pbciter)
  {
    if (pbciter->first == nodegid) // node is a pbc master node
    {
      pbcnodefound = true;
      // coupled node is the (first) slave node
      coupnodegid = pbciter->second[0];
      if (pbciter->second.size()!=1) dserror("this might need to be modified for more than 1 slave per master");
    }
    else
    {
      // loop all slave nodes
      for (size_t islave=0;islave<pbciter->second.size();islave++)
      {
        if (pbciter->second[islave] == nodegid) // node is a pbc slave node
        {
          pbcnodefound = true;
          coupnodegid = pbciter->first; // coupled node is the master node
        }
      } // end loop over slave nodes
    } // end if
  } // end loop over pbc map

  // get pbcnode
  if (pbcnodefound)
    pbcnode = discret_->gNode(coupnodegid);
  } // end if pbcnode true}



/*------------------------------------------------------------------------------------------------*
 * reset a special state with another state in the data class                    winklmaier 10/11 *
 *------------------------------------------------------------------------------------------------*/
void XFEM::TIMEINT::resetState(
    TimeIntData::state oldState,
    TimeIntData::state newState
) const
{
  for (std::vector<TimeIntData>::iterator data=timeIntData_->begin();
      data!=timeIntData_->end(); data++)
  {
    if (data->state_ == oldState)
      data->state_ = newState;
  }
} // end function resetState



/*------------------------------------------------------------------------------------------------*
 * clear the data of all nodes having a special state                            winklmaier 10/11 *
 *------------------------------------------------------------------------------------------------*/
void XFEM::TIMEINT::clearState(
    TimeIntData::state state
) const
{
  std::vector<TimeIntData>::iterator data;
  while(true) // while loop over data to be cleared
  {
    for (data=timeIntData_->begin();
        data!=timeIntData_->end(); data++) // for loop over all data
    {
      if (data->state_==state)
      {
        timeIntData_->erase(data);
        break;
      }
    } // end for loop over all data
    if (data==timeIntData_->end())
      break;
  } // end while loop over data to be cleared
} // end function clear state



#ifdef PARALLEL
/*------------------------------------------------------------------------------------------------*
 * basic function sending data to dest and receiving data from source            winklmaier 10/10 *
 *------------------------------------------------------------------------------------------------*/
void XFEM::TIMEINT::sendData(
    DRT::PackBuffer& dataSend,
    int& dest,
    int& source,
    std::vector<char>& dataRecv
) const
{
  std::vector<int> lengthSend(1,0);
  lengthSend[0] = dataSend().size();
  int size_one = 1;

#ifdef DEBUG
  std::cout << "--- sending "<< lengthSend[0] << " bytes: from proc " << myrank_ << " to proc " << dest << std::endl;
#endif

  // exporter for sending
  DRT::Exporter exporter(discret_->Comm());

  // send length of the data to be received ...
  MPI_Request req_length_data;
  int length_tag = 0;
  exporter.ISend(myrank_, dest, &(lengthSend[0]) , size_one, length_tag, req_length_data);
  // ... and receive length
  std::vector<int> lengthRecv(1,0);
  exporter.Receive(source, length_tag, lengthRecv, size_one);
  exporter.Wait(req_length_data);

  // send actual data ...
  int data_tag = 4;
  MPI_Request req_data;
  exporter.ISend(myrank_, dest, &(dataSend()[0]), lengthSend[0], data_tag, req_data);

  // ... and receive data
  dataRecv.clear(); dataRecv.resize(lengthRecv[0]);
  exporter.ReceiveAny(source, data_tag, dataRecv, lengthRecv[0]);
  exporter.Wait(req_data);

#ifdef DEBUG
  std::cout << "--- receiving "<< lengthRecv[0] << " bytes: to proc " << myrank_ << " from proc " << source << std::endl;
#endif
} // end sendData



/*------------------------------------------------------------------------------------------------*
 * packing a node for parallel communication only with the basic nodal data                       *
 * without an underlying discretization fitting to the node's new prozessor      winklmaier 10/10 *
 *------------------------------------------------------------------------------------------------*/
void XFEM::TIMEINT::packNode(
    DRT::PackBuffer& dataSend,
    DRT::Node& node
) const
{
  const int nsd = 3;
  DRT::ParObject::AddtoPack(dataSend,node.Id());
  DRT::ParObject::AddtoPack(dataSend,LINALG::Matrix<nsd,1>(node.X()));
  DRT::ParObject::AddtoPack(dataSend,node.Owner());
} // end packNode



/*------------------------------------------------------------------------------------------------*
 * unpacking a node after parallel communication only with the basic nodal data                   *
 * without an underlying discretization fitting to the node's new prozessor      winklmaier 10/10 *
 *------------------------------------------------------------------------------------------------*/
void XFEM::TIMEINT::unpackNode(
    std::vector<char>::size_type& posinData,
    std::vector<char>& dataRecv,
    DRT::Node& node
) const
{
  const int nsd = 3; // dimension
  int id; // global id
  LINALG::Matrix<nsd,1> coords; // coordinates
  int owner; // processor

  DRT::ParObject::ExtractfromPack(posinData,dataRecv,id);
  DRT::ParObject::ExtractfromPack(posinData,dataRecv,coords);
  DRT::ParObject::ExtractfromPack(posinData,dataRecv,owner);

  if (owner==myrank_) // real node with all data
  {
    node = *discret_->gNode(id);
  }
  else // just id, coords and owner
  {
    double coordinates[nsd];
    for (int dim=0;dim<nsd;dim++)
      coordinates[dim] = coords(dim);

    DRT::Node newNode(id,coordinates,owner);
    node = newNode;
  } // end if correct processor
} // end function unpackNode
#endif // PARALLEL



/*------------------------------------------------------------------------------------------------*
 * basic XFEM time-integration constructor for standard degrees of freedom       winklmaier 11/11 *
 *------------------------------------------------------------------------------------------------*/
XFEM::STD::STD(
    XFEM::TIMEINT& timeInt,
    INPAR::COMBUST::XFEMTimeIntegration& timeIntType,
    const Teuchos::RCP<Epetra_Vector> veln,
    const double& dt,
    const Teuchos::RCP<COMBUST::FlameFront> flamefront,
    bool initialize
) :
XFEM::TIMEINT::TIMEINT(timeInt),
timeIntType_(timeIntType),
veln_(veln),
dt_(dt),
flamefront_(flamefront)
{
  if (initialize)
  {
    const int nsd = 3; // dimension
    timeIntData_ = Teuchos::rcp(new std::vector<TimeIntData>); // vector containing all data used for computation

    /*------------------------*
     * Initialization         *
     *------------------------*/
    for (int leleid=0; leleid<discret_->NumMyColElements(); leleid++)  // loop over processor nodes
    {
      DRT::Element* iele = discret_->lColElement(leleid);
      if (intersectionStatus(iele)!=TIMEINT::uncut_) // element cut
      {
        const int* nodeGids = iele->NodeIds(); // node gids
        for (int inode=0;inode<iele->NumNode();inode++) // loop over element nodes
        {
          if(olddofcolmap_.MyGID(nodeGids[inode]))
          {}
          oldEnrNodes_.insert(nodeGids[inode]);
        } // end loop over element nodes
      } // end if element cut
    } // end loop over processor nodes

    LINALG::Matrix<nsd,1> dummyStartpoint; // dummy startpoint for comparison
    for (int i=0;i<nsd;i++) dummyStartpoint(i) = 777.777;

    // fill curr_ structure with the data for the nodes which changed interface side
    for (int lnodeid=0; lnodeid<discret_->NumMyColNodes(); lnodeid++)  // loop over processor nodes
    {
      DRT::Node* currnode = discret_->lColNode(lnodeid); // current analysed node

      // node on current processor which changed interface side
      if ((currnode->Owner() == myrank_) &&
          (interfaceSideCompare((*phinp_)[lnodeid],(*phin_)[lnodeid]) == false))
      {
        timeIntData_->push_back(TimeIntData(
            *currnode,
            LINALG::Matrix<nsd,1>(true),
            std::vector<LINALG::Matrix<nsd,nsd> >(oldVectors_.size(),LINALG::Matrix<nsd,nsd>(true)),
            std::vector<LINALG::Matrix<1,nsd> >(oldVectors_.size(),LINALG::Matrix<1,nsd>(true)),
            dummyStartpoint,
            (*phinp_)[lnodeid],
            1,
            0,
            std::vector<int>(1,-1),
            std::vector<int>(1,-1),
            INFINITY,
            TimeIntData::predictor_)); // data created for the node
      }
    } // end loop over processor nodes

    startpoints();

    // test loop if all initial startpoints have been computed
    for (std::vector<TimeIntData>::iterator data=timeIntData_->begin();
        data!=timeIntData_->end(); data++)
    {
      if (data->startpoint_==dummyStartpoint) // startpoint unchanged
        dserror("WARNING! No enriched node on one interface side found!\nThis "
            "indicates that the whole area is at one side of the interface!");
    } // end loop over nodes
  }
  return;
} // end constructor



/*------------------------------------------------------------------------------------------------*
 * initialize data when called in a new FGI                                      winklmaier 10/11 *
 *------------------------------------------------------------------------------------------------*/
void XFEM::STD::importNewFGIData(
    const Teuchos::RCP<DRT::Discretization> discret,
    const Teuchos::RCP<XFEM::DofManager> newdofman,
    const Teuchos::RCP<COMBUST::FlameFront> flamefront,
    const Epetra_Map& newdofrowmap,
    const std::map<DofKey, DofGID>& newNodalDofRowDistrib)
{
  discret_ = discret;
  newdofman_ = newdofman;
  phinpi_ = phinp_;
  phinp_ = flamefront->Phinp();
  flamefront_ = flamefront;
  newinterfacehandle_ = flamefront->InterfaceHandle();
  newdofrowmap_ = newdofrowmap;
  // must loop because XFEM::DofKey does not implement a copy assignment operator
  newNodalDofRowDistrib_.clear();
  for (std::map<XFEM::DofKey, XFEM::DofGID>::const_iterator i=newNodalDofRowDistrib.begin();
       i != newNodalDofRowDistrib.end(); ++i)
    newNodalDofRowDistrib_[i->first] = i->second;
  return;
}



/*------------------------------------------------------------------------------------------------*
 * identify an element containing a point and additional data                    winklmaier 06/10 *
 *------------------------------------------------------------------------------------------------*/
void XFEM::STD::elementSearch(
    DRT::Element*& ele,
    LINALG::Matrix<3,1>& x,
    LINALG::Matrix<3,1>& xi,
    LINALG::Matrix<3,1>& vel,
    double& phi,
    bool& found
) const
{
  DRT::Element* currele = NULL; // current element

  int startid; // local row element id
  if (ele==NULL) startid = 0; // start with first local row element
  else           startid = -1; // pseudo-id so that id+1 will be 0

  //loop over elements
  for (int ieleid = startid;ieleid<discret_->NumMyColElements();ieleid++)
  {
    // if ele != NULL and so initialized,
    // first it should be checked if it is fitting
    if (ieleid == -1)
      currele = ele;
    else
      currele = discret_->lColElement(ieleid);

    //    std::cout << "currently analysed element" << *currele;
    //    std::cout << "startpoint approximation: " << x;
    callXToXiCoords(currele,x,xi,found);
    //    std::cout << "in local xi-coordinates: " << xi;

    if (found)
    {
      ele = currele;
      getGPValues(ele,xi,vel,phi);
      break;
    }
  } // end loop over processor elements
} // end function findElementAndLocalCoords



/*------------------------------------------------------------------------------------------------*
 * interpolate velocity and phi-value for a point in an element                  winklmaier 11/11 *
 *------------------------------------------------------------------------------------------------*/
void XFEM::STD::getGPValues(
    DRT::Element* ele,
    LINALG::Matrix<3,1>& xi,
    LINALG::Matrix<3,1>& vel,
    double& phi
) const
{
  switch (ele->Shape())
  {
  case DRT::Element::hex8:
    getGPValues<DRT::Element::hex8>(ele,xi,vel,phi);
    break;
  case DRT::Element::hex20:
    getGPValues<DRT::Element::hex20>(ele,xi,vel,phi);
    break;
  case DRT::Element::tet4:
    getGPValues<DRT::Element::tet4>(ele,xi,vel,phi);
    break;
  default:
    dserror("add your 3D distype here!");
    break;
  } // end switch
} // end function getGPValues



/*------------------------------------------------------------------------------------------------*
 * Compute starting values for the interface-changing nodes                      winklmaier 06/10 *
 *------------------------------------------------------------------------------------------------*/
void XFEM::STD::startpoints()
{
  //Initialization
  const int nsd = 3; // 3 dimensions for a 3d fluid element
  const double TOL = 1.0e-3; // tolerance

  // loop over processors
  for (int procid=0; procid<numproc_; procid++)
  {
    // loop over nodes which changed interface side
    for (std::vector<TimeIntData>::iterator data=timeIntData_->begin();
        data!=timeIntData_->end(); data++)
    {
      if (data->state_==TimeIntData::basicStd_) // correct state
      {
        LINALG::Matrix<nsd,1> newNodeCoords(data->node_.X()); // coords of endpoint of Lagrangian characteristics

        // loop over intersected elements on processor
        for (std::set<int>::const_iterator enrnode = oldEnrNodes_.begin();
            enrnode != oldEnrNodes_.end();
            enrnode++)
        {
          DRT::Node* lnodeold = discret_->gNode(*enrnode);//elenodeids[elenode]);  // node near interface
          LINALG::Matrix<nsd,1> oldNodeCoords(lnodeold->X());  // coords of potential startpoint

          // just look for points on the same interface side and on the same processor in row map
          if ((interfaceSideCompare(data->phiValue_,(*phin_)[lnodeold->LID()])) and
              (lnodeold->Owner()==myrank_))
          {
            LINALG::Matrix<nsd,1> diff;  // vector from old point at time n to new point at time n+1
            diff.Update(1.0,newNodeCoords,-1.0,oldNodeCoords);

            if (diff.Norm2() < 2*data->dMin_) // possible new nearest node
            {
              // get nodal velocity of a node
              LINALG::Matrix<nsd,1> nodevel(true);  // velocity of "old" node at time n
              const std::set<XFEM::FieldEnr>& fieldenrset(olddofman_->getNodeDofSet(*enrnode));
              for (std::set<XFEM::FieldEnr>::const_iterator fieldenr = fieldenrset.begin();
                  fieldenr != fieldenrset.end();++fieldenr)
              {
                const DofKey olddofkey(*enrnode,*fieldenr);
                const int olddofpos = oldNodalDofColDistrib_.find(olddofkey)->second;

                if (fieldenr->getEnrichment().Type() == XFEM::Enrichment::typeStandard)
                {
                  if (fieldenr->getField() == XFEM::PHYSICS::Velx)
                    nodevel(0) = (*veln_)[olddofcolmap_.LID(olddofpos)];
                  if (fieldenr->getField() == XFEM::PHYSICS::Vely)
                    nodevel(1) = (*veln_)[olddofcolmap_.LID(olddofpos)];
                  if (fieldenr->getField() == XFEM::PHYSICS::Velz)
                    nodevel(2) = (*veln_)[olddofcolmap_.LID(olddofpos)];
                }
              } // end loop over fieldenr

              LINALG::Matrix<1,1> arc(true); // cosinus of angle between dist and vel(x_n+1)*diff.Norm2()
              if (nodevel.Norm2()/diff.Norm2()>1e-2)
                arc.MultiplyTN(1.0/nodevel.Norm2(),diff,nodevel);

              double dist = diff.Norm2() + (diff.Norm2()-arc.Norm2()); // distance, containing arc influence

              if (dist-data->dMin_+TOL*dist < 0) // new nearest node (cosinus shall be near 1!)
              {
                data->startpoint_.Update(1.0,newNodeCoords,-dt_,nodevel);
                data->dMin_ = dist;
                data->startGid_.clear();
                data->startGid_.push_back(*enrnode);
                data->startOwner_.clear();
                data->startOwner_.push_back(myrank_);
              } // end if new best startvalue
              else if ((dist-data->dMin_ > -TOL*dist) and
                  (dist-data->dMin_ < TOL*dist)) // handles special case that two nodes are very similar near
              {
                data->startpoint_.Update(0.5,newNodeCoords,-0.5*dt_,nodevel,0.5); // midpoint is 0.5*(x-dt*vel+old_startpoint)
                data->dMin_ = (dist+data->dMin_)/2.0;
                data->startGid_.push_back(*enrnode);
                data->startOwner_.push_back(myrank_);
              } // end if new equal best startvalue
            } // end if possibly new best startvalue
          } // end if just points on the correct side
        } // end loop over intersected elements
      } // end if correct state
    } // end loop over nodes which changed interface side
#ifdef PARALLEL
    exportStartData();
#endif
  } // end loop over processors
} // end startValuesFinder



/*------------------------------------------------------------------------------------------------*
 * setting the computed data for the standard degrees of freedom into the according               *
 * Epetra Vectors for all handled nodes                                          winklmaier 06/10 *
 *------------------------------------------------------------------------------------------------*/
void XFEM::STD::setFinalData(
)
{
  const int nsd = 3; // 3 dimensions for a 3d fluid element
  std::map<int,int> usedStartpoints; // map containing the used start points
  int numStartpoints; // number of used start points
  double newValue = 0.0; // new value

  // loop over data
  for (std::vector<TimeIntData>::iterator data=timeIntData_->begin();
      data!=timeIntData_->end(); data++)
  {
    if (data->state_==TimeIntData::doneStd_)
    {
      std::vector<LINALG::Matrix<nsd,1> >& velValues(data->velValues_); // velocities of the node
      std::vector<double>& presValues(data->presValues_); // pressures of the node

      const int gnodeid = data->node_.Id(); // global node id

      std::map<int,int>::iterator currstartpoint = usedStartpoints.find(gnodeid); // current start point
      if (currstartpoint==usedStartpoints.end()) // standard case and "standard alternative" case
      {
        usedStartpoints.insert(std::pair<int,int>(gnodeid,1));
        numStartpoints = 1;
      }
      else
      {
        currstartpoint->second +=1;
        numStartpoints = currstartpoint->second;
      }

      // set nodal velocities and pressures with help of the field set of node
      const std::set<XFEM::FieldEnr>& fieldenrset(newdofman_->getNodeDofSet(gnodeid));
      for (std::set<XFEM::FieldEnr>::const_iterator fieldenr = fieldenrset.begin();
          fieldenr != fieldenrset.end();++fieldenr) // loop over field enr set
      {
        const DofKey newdofkey(gnodeid, *fieldenr);
        const int newdofpos = newNodalDofRowDistrib_.find(newdofkey)->second;

        if (fieldenr->getEnrichment().Type() == XFEM::Enrichment::typeStandard)
        {
          /*
        if (fieldenr->getField() == XFEM::PHYSICS::Velx)
          std::cout << (*newVectors_[0])[newdofrowmap_.LID(newdofpos)] << " becomes " << velValues[0](0) << std::endl;
        else if (fieldenr->getField() == XFEM::PHYSICS::Vely)
          std::cout << (*newVectors_[0])[newdofrowmap_.LID(newdofpos)] << " becomes " << velValues[0](1) << std::endl;
        else if (fieldenr->getField() == XFEM::PHYSICS::Velz)
          std::cout << (*newVectors_[0])[newdofrowmap_.LID(newdofpos)] << " becomes " << velValues[0](2) << std::endl;
        else if (fieldenr->getField() == XFEM::PHYSICS::Pres)
          std::cout << (*newVectors_[0])[newdofrowmap_.LID(newdofpos)] << " becomes " << presValues[0] << std::endl;
           */


          // loop over vectors to be set (either all vectors or only the vectors at t^n+1)
          for (size_t index=0;index<vectorSize(data->type_);index++)
          {
            double value = (*newVectors_[index])[newdofrowmap_.LID(newdofpos)];
            if (fieldenr->getField() == XFEM::PHYSICS::Velx)
              newValue = ((numStartpoints-1.0)/numStartpoints)*value
              + velValues[index](0)/numStartpoints;
            else if (fieldenr->getField() == XFEM::PHYSICS::Vely)
              newValue = ((numStartpoints-1.0)/numStartpoints)*value
              + velValues[index](1)/numStartpoints;
            else if (fieldenr->getField() == XFEM::PHYSICS::Velz)
              newValue = ((numStartpoints-1.0)/numStartpoints)*value
              + velValues[index](2)/numStartpoints;
            else if (fieldenr->getField() == XFEM::PHYSICS::Pres)
              newValue = ((numStartpoints-1.0)/numStartpoints)*value
              + presValues[index]/numStartpoints;

            (*newVectors_[index])[newdofrowmap_.LID(newdofpos)] = newValue; // set the value
          }
        }
      } // end loop over fieldenr
      data->type_ = TimeIntData::standard_; // predictor is done, so next time standard
    }
  } // end loop over nodes
} // end setFinalData



#ifdef PARALLEL
/*------------------------------------------------------------------------------------------------*
 * export start data to neighbour proc                                           winklmaier 06/10 *
 *------------------------------------------------------------------------------------------------*/
void XFEM::STD::exportStartData()
{
  const int nsd = 3; // 3 dimensions for a 3d fluid element

  // destination proc (the "next" one)
  int dest = myrank_+1;
  if(myrank_ == (numproc_-1))
    dest = 0;

  // source proc (the "last" one)
  int source = myrank_-1;
  if(myrank_ == 0)
    source = numproc_-1;

  DRT::PackBuffer dataSend; // data to be sent

  // packing the data
  for (std::vector<TimeIntData>::iterator data=timeIntData_->begin();
      data!=timeIntData_->end(); data++)
  {
    packNode(dataSend,data->node_);
    DRT::ParObject::AddtoPack(dataSend,data->vel_);
    DRT::ParObject::AddtoPack(dataSend,data->velDeriv_);
    DRT::ParObject::AddtoPack(dataSend,data->presDeriv_);
    DRT::ParObject::AddtoPack(dataSend,data->startpoint_);
    DRT::ParObject::AddtoPack(dataSend,data->phiValue_);
    DRT::ParObject::AddtoPack(dataSend,data->searchedProcs_);
    DRT::ParObject::AddtoPack(dataSend,data->counter_);
    DRT::ParObject::AddtoPack(dataSend,data->startGid_);
    DRT::ParObject::AddtoPack(dataSend,data->startOwner_);
    DRT::ParObject::AddtoPack(dataSend,data->dMin_);
    DRT::ParObject::AddtoPack(dataSend,(int)data->type_);
  }

  dataSend.StartPacking();

  for (std::vector<TimeIntData>::iterator data=timeIntData_->begin();
      data!=timeIntData_->end(); data++)
  {
    packNode(dataSend,data->node_);
    DRT::ParObject::AddtoPack(dataSend,data->vel_);
    DRT::ParObject::AddtoPack(dataSend,data->velDeriv_);
    DRT::ParObject::AddtoPack(dataSend,data->presDeriv_);
    DRT::ParObject::AddtoPack(dataSend,data->startpoint_);
    DRT::ParObject::AddtoPack(dataSend,data->phiValue_);
    DRT::ParObject::AddtoPack(dataSend,data->searchedProcs_);
    DRT::ParObject::AddtoPack(dataSend,data->counter_);
    DRT::ParObject::AddtoPack(dataSend,data->startGid_);
    DRT::ParObject::AddtoPack(dataSend,data->startOwner_);
    DRT::ParObject::AddtoPack(dataSend,data->dMin_);
    DRT::ParObject::AddtoPack(dataSend,(int)data->type_);
  }

  std::vector<char> dataRecv;
  sendData(dataSend,dest,source,dataRecv);

  // pointer to current position of group of cells in global std::string (counts bytes)
  std::vector<char>::size_type posinData = 0;

  // clear vector that should be filled
  timeIntData_->clear();

  // unpack received data
  while (posinData < dataRecv.size())
  {
    double coords[nsd] = {0.0};
    DRT::Node node(0,(double*)coords,0); // initialize node
    LINALG::Matrix<nsd,1> vel; // velocity at point x
    std::vector<LINALG::Matrix<nsd,nsd> > velDeriv; // derivation of velocity at point x
    std::vector<LINALG::Matrix<1,nsd> > presDeriv; // derivation of pressure at point x
    LINALG::Matrix<nsd,1> startpoint; // startpoint
    double phiValue; // phi-value
    int searchedProcs; // number of searched processors
    int counter; // iteration counter
    std::vector<int> startGid; // global id of first startpoint
    std::vector<int> startOwner; // owner of first startpoint
    double dMin; // minimal distance
    int newtype; // type of the data

    unpackNode(posinData,dataRecv,node);
    DRT::ParObject::ExtractfromPack(posinData,dataRecv,vel);
    DRT::ParObject::ExtractfromPack(posinData,dataRecv,velDeriv);
    DRT::ParObject::ExtractfromPack(posinData,dataRecv,presDeriv);
    DRT::ParObject::ExtractfromPack(posinData,dataRecv,startpoint);
    DRT::ParObject::ExtractfromPack(posinData,dataRecv,phiValue);
    DRT::ParObject::ExtractfromPack(posinData,dataRecv,searchedProcs);
    DRT::ParObject::ExtractfromPack(posinData,dataRecv,counter);
    DRT::ParObject::ExtractfromPack(posinData,dataRecv,startGid);
    DRT::ParObject::ExtractfromPack(posinData,dataRecv,startOwner);
    DRT::ParObject::ExtractfromPack(posinData,dataRecv,dMin);
    DRT::ParObject::ExtractfromPack(posinData,dataRecv,newtype);

    timeIntData_->push_back(TimeIntData(
        node,
        vel,
        velDeriv,
        presDeriv,
        startpoint,
        phiValue,
        searchedProcs,
        counter,
        startGid,
        startOwner,
        dMin,
        (TimeIntData::type)newtype));
  }

  discret_->Comm().Barrier(); // processors wait for each other
} // end exportStartData



/*------------------------------------------------------------------------------------------------*
 * export final data to node's proc                                              winklmaier 06/10 *
 *------------------------------------------------------------------------------------------------*/
void XFEM::STD::exportFinalData()
{
  const int nsd = 3; // 3 dimensions for a 3d fluid element

  // array of vectors which stores data for
  // every processor in one vector
  std::vector<std::vector<TimeIntData> > dataVec(numproc_);

  // fill vectors with the data
  for (std::vector<TimeIntData>::iterator data=timeIntData_->begin();
      data!=timeIntData_->end(); data++)
  {
    if (data->state_==TimeIntData::doneStd_)
      dataVec[data->node_.Owner()].push_back(*data);
  }

  timeIntData_->clear();
  *timeIntData_ = dataVec[myrank_]; // set final data of own processor
  dataVec[myrank_].clear(); // clear data about current proc

  // send data to the processor where the point lies (1. nearest higher neighbour 2. 2nd nearest higher neighbour...)
  for (int dest=(myrank_+1)%numproc_;dest!=myrank_;dest=(dest+1)%numproc_) // dest is the target processor
  {
    // Initialization
    int source = myrank_-(dest-myrank_); // source proc (sends (dest-myrank_) far and gets from (dest-myrank_) earlier)
    if (source<0)                            source+=numproc_;
    else if (source>=numproc_)       source-=numproc_;

    DRT::PackBuffer dataSend;

    // pack data to be sent
    for (std::vector<TimeIntData>::iterator data=dataVec[dest].begin();
        data!=dataVec[dest].end(); data++)
    {
      DRT::ParObject::AddtoPack(dataSend,data->node_.Id());
      DRT::ParObject::AddtoPack(dataSend,data->startpoint_);
      DRT::ParObject::AddtoPack(dataSend,data->phiValue_);
      DRT::ParObject::AddtoPack(dataSend,data->startGid_);
      DRT::ParObject::AddtoPack(dataSend,data->startOwner_);
      DRT::ParObject::AddtoPack(dataSend,data->velValues_);
      DRT::ParObject::AddtoPack(dataSend,data->presValues_);
      DRT::ParObject::AddtoPack(dataSend,data->type_);
    }

    dataSend.StartPacking();

    for (std::vector<TimeIntData>::iterator data=dataVec[dest].begin();
        data!=dataVec[dest].end(); data++)
    {
      DRT::ParObject::AddtoPack(dataSend,data->node_.Id());
      DRT::ParObject::AddtoPack(dataSend,data->startpoint_);
      DRT::ParObject::AddtoPack(dataSend,data->phiValue_);
      DRT::ParObject::AddtoPack(dataSend,data->startGid_);
      DRT::ParObject::AddtoPack(dataSend,data->startOwner_);
      DRT::ParObject::AddtoPack(dataSend,data->velValues_);
      DRT::ParObject::AddtoPack(dataSend,data->presValues_);
      DRT::ParObject::AddtoPack(dataSend,data->type_);
    }

    // clear the no more needed data
    dataVec[dest].clear();

    std::vector<char> dataRecv;
    sendData(dataSend,dest,source,dataRecv);

    // pointer to current position of group of cells in global std::string (counts bytes)
    std::vector<char>::size_type posinData = 0;

    // unpack received data
    while (posinData < dataRecv.size())
    {
      int gid; // global id of node
      LINALG::Matrix<nsd,1> startpoint; // startpoint
      double phiValue; // phi-value
      std::vector<int> startGid; // global id of first startpoint
      std::vector<int> startOwner; // owner of first startpoint
      std::vector<LINALG::Matrix<nsd,1> > velValues; // velocity values
      std::vector<double> presValues; // pressure values
      int newtype; // type of the data

      DRT::ParObject::ExtractfromPack(posinData,dataRecv,gid);
      DRT::ParObject::ExtractfromPack(posinData,dataRecv,startpoint);
      DRT::ParObject::ExtractfromPack(posinData,dataRecv,phiValue);
      DRT::ParObject::ExtractfromPack(posinData,dataRecv,startGid);
      DRT::ParObject::ExtractfromPack(posinData,dataRecv,startOwner);
      DRT::ParObject::ExtractfromPack(posinData,dataRecv,velValues);
      DRT::ParObject::ExtractfromPack(posinData,dataRecv,presValues);
      DRT::ParObject::ExtractfromPack(posinData,dataRecv,newtype);

      timeIntData_->push_back(TimeIntData(
          *discret_->gNode(gid),
          startpoint,
          phiValue,
          startGid,
          startOwner,
          velValues,
          presValues,
          (TimeIntData::type)newtype));
    } // end loop over number of nodes to get

    // processors wait for each other
    discret_->Comm().Barrier();
  } // end loop over processors
} // end exportfinalData
#endif //parallel



/*------------------------------------------------------------------------------------------------*
 * basic XFEM time-integration constructor for enrichment degrees of freedom     winklmaier 11/11 *
 *------------------------------------------------------------------------------------------------*/
XFEM::ENR::ENR(
    XFEM::TIMEINT& timeInt,
    INPAR::COMBUST::XFEMTimeIntegrationEnr& timeIntEnr,
    INPAR::COMBUST::XFEMTimeIntegrationEnrComp& timeIntEnrType
) :
XFEM::TIMEINT::TIMEINT(timeInt),
timeIntEnr_(timeIntEnr),
timeIntEnrType_(timeIntEnrType),
critTol_(1.0e-02)


{
  // Initialization
  if (timeIntEnrType_==INPAR::COMBUST::xfemtimeintenr_standard)
    getCritCutElements();

  timeIntData_ = Teuchos::rcp(new std::vector<TimeIntData>);
  return;
} // end constructor



/*------------------------------------------------------------------------------------------------*
 * out of order!                                                                 winklmaier 10/11 *
 *------------------------------------------------------------------------------------------------*/
void XFEM::ENR::compute(
    std::vector<Teuchos::RCP<Epetra_Vector> > newRowVectorsn,
    std::vector<Teuchos::RCP<Epetra_Vector> > newRowVectorsnp
)
{
  dserror("Unused function! Use a function of the derived classes");
} // end function compute



/*------------------------------------------------------------------------------------------------*
 * initialize data when called in a new FGI                                      winklmaier 10/11 *
 *------------------------------------------------------------------------------------------------*/
void XFEM::ENR::importNewFGIData(
    const Teuchos::RCP<DRT::Discretization> discret,
    const Teuchos::RCP<XFEM::DofManager> newdofman,
    const Teuchos::RCP<COMBUST::FlameFront> flamefront,
    const Epetra_Map& newdofrowmap,
    const std::map<DofKey, DofGID>& newNodalDofRowDistrib,
    const std::map<DofKey, DofGID>& oldNodalDofColDistrib
)
{
  discret_=discret;
  newdofman_=newdofman;
  phinp_=flamefront->Phinp();
  newinterfacehandle_=flamefront->InterfaceHandle();
  newdofrowmap_=newdofrowmap;
  // must loop because XFEM::DofKey does not implement a copy assignment operator
  newNodalDofRowDistrib_.clear();
  for (std::map<XFEM::DofKey, XFEM::DofGID>::const_iterator i=newNodalDofRowDistrib.begin();
       i != newNodalDofRowDistrib.end(); ++i)
    newNodalDofRowDistrib_[i->first] = i->second;
  return;
}



/*------------------------------------------------------------------------------------------------*
 * give back a boolean indicating if a new enrichment value shall be computed.                    *
 * use the FGI FRS type for this analysis as well as data about the node         winklmaier 11/11 *
 *------------------------------------------------------------------------------------------------*/
bool XFEM::ENR::newEnrValueNeeded(
    const DRT::Node* node
) const
{
  switch (FGIType_)
  {
  case FRS1FGI1_:
  {
    // case 1: the node was not enriched in old timestep and therefore needs a new enrichment
    const int gid = node->Id();
    //         std::cout << "here with node " << *node << std::endl;

    const std::set<XFEM::FieldEnr>& fieldenrset(newdofman_->getNodeDofSet(gid)); // field set of node
    for (std::set<XFEM::FieldEnr>::const_iterator fieldenr = fieldenrset.begin();
        fieldenr != fieldenrset.end();++fieldenr)
    {
      const DofKey newdofkey(gid, *fieldenr);

      if (fieldenr->getEnrichment().Type() != XFEM::Enrichment::typeStandard)
      {
        if (timeIntEnrType_==INPAR::COMBUST::xfemtimeintenr_full) // every enriched node gets new enrichment values with this flag
          return true;

        // additional checks for standard or minimal enrichment computation:

        // check if enrichment dofs existed before
        std::map<DofKey, DofGID>::const_iterator olddof = oldNodalDofColDistrib_.find(newdofkey);
        if (olddof == oldNodalDofColDistrib_.end()) // olddof not found -> no enr value before -> enr value has to be set
          return true;

        if (timeIntEnrType_==INPAR::COMBUST::xfemtimeintenr_standard)
          return critCut(node);
        //      std::cout << "bool is " << critCut << std::endl;
      } // end if dof has non standard enrichment
    } // end loop over fieldset
    break;
  } // end case FRS1FGI1
  case FRS1FGINot1_:
  {
    // case 1: the node was not enriched in old timestep and therefore needs a new enrichment
    const int gid = node->Id();
    //         std::cout << "here with node " << *node << std::endl;

    const std::set<XFEM::FieldEnr>& fieldenrset(newdofman_->getNodeDofSet(gid)); // field set of node
    for (std::set<XFEM::FieldEnr>::const_iterator fieldenr = fieldenrset.begin();
        fieldenr != fieldenrset.end();++fieldenr)
    {
      const DofKey newdofkey(gid, *fieldenr);

      if (fieldenr->getEnrichment().Type() != XFEM::Enrichment::typeStandard)
      {
        // check if enrichment dofs existed before
        std::map<DofKey, DofGID>::const_iterator dof_npi = nodalDofColDistrib_npi_.find(newdofkey);
        if (dof_npi == nodalDofColDistrib_npi_.end()) // olddof not found -> no enr value before -> enr value has to be set
        {
          std::map<DofKey, DofGID>::const_iterator dof_n = oldNodalDofColDistrib_.find(newdofkey);
          if (dof_n != oldNodalDofColDistrib_.end()) // dof existed at t^n so use this value if it is not critical
          {
            if ((timeIntEnrType_==INPAR::COMBUST::xfemtimeintenr_standard) and (critCut(node))) // critical cut has to be recomputed
              return true;
            else // use existing value of old solution which was deleted by former FGI
            {
              const int newdofpos = newNodalDofRowDistrib_.find(newdofkey)->second;
              for (size_t index=0;index<newVectors_.size();index++)
                (*newVectors_[index])[newdofrowmap_.LID(newdofpos)] =
                    (*oldVectors_[index])[olddofcolmap_.LID(dof_n->second)];
              return false;
            }
          }
          else // node enriched, was not enriched at t^n+1,i and t^n
            return true;
        }
      } // end if dof has non standard enrichment
    } // end loop over fieldset
    break;
  }
  case FRSNot1_: return false;
  default:
    dserror("undefined type");
    break;
  } // end switch
  return false;
} // end function newEnrValueNeeded



/*------------------------------------------------------------------------------------------------*
 * give back a boolean indicating if all cut elements around the node, which has to be
 * enriched, are critically cut                                                  winklmaier 11/11 *
 *------------------------------------------------------------------------------------------------*/
bool XFEM::ENR::critCut(const DRT::Node* node) const
{
  // case 2: the node was enriched but the according support is very small.
  // Then the enrichment value is potentially much too high and therefore needs a new value.
  bool critCut = false; // true if intersected ele around and all intersected eles have small support for enr shape fcn
  // vector of elements located around this node
  std::vector<const DRT::Element*> eles;
  addPBCelements(node,eles);
  const int numeles=eles.size();

  // determine interface side of node
  bool domainPlus = plusDomain((*phin_)[node->LID()]); // indicator if node is in omega^+

  for (int iele=0;iele<numeles;iele++) // loop over elements around node
  {
    if (intersectionStatus(eles[iele])==XFEM::TIMEINT::cut_) // element intersected
    {
      if (domainPlus) // when node is in plus domain, support of enrichment is the minus part of the element
      {
        std::set<int>::const_iterator tmp = critElesMinus_.find(eles[iele]->Id());
        if (tmp == critElesMinus_.end()) // volumes just saved in critical cases
        {
          critCut = false; // one element has a big enough support for the node -> no problem with values
          break;
        }
        else
          critCut = true;
      } // end if node in plus domain
      else
      {
        std::set<int>::const_iterator tmp = critElesPlus_.find(eles[iele]->Id());
        if (tmp == critElesPlus_.end()) // volumes just saved in critical cases
        {
          critCut = false; // one element has a big enough support for the node -> no problem with values
          break;
        }
        else
          critCut = true;
      } // end if node in minus domain
    } // end if element bisected or touched
  } // end loop over elements around node
  //  if (critCut) std::cout << *node << " is newly set because of critical cut" << std::endl;
  return critCut;
} // end function critCut



/*------------------------------------------------------------------------------------------------*
 * set the critical cut elements with critical part omega^- or omega^+           winklmaier 08/10 *
 *------------------------------------------------------------------------------------------------*/
void XFEM::ENR::getCritCutElements(
)
{
  double plusVol; // volume of integration cells lying in omega^+
  double minusVol; // volume of integration cells lying in omega^-
  double currVol; // current volume of integration cell
  double eleVol; // element volume

  DRT::Element* currEle = NULL; // current element

  for (int iele=0;iele<discret_->NumMyColElements();iele++) // loop over elements
  {
    currEle = discret_->lColElement(iele);

    if (intersectionStatus(currEle)==XFEM::TIMEINT::cut_) // element bisected
    {
      plusVol = 0.0;
      minusVol = 0.0;

      // get domain integration cells for this element
      const GEO::DomainIntCells&  domainIntCells(oldinterfacehandle_->ElementDomainIntCells(currEle->Id())); // domain integration cells of bisected element

      // loop over domain integration cells
      for (GEO::DomainIntCells::const_iterator cell = domainIntCells.begin(); cell != domainIntCells.end(); ++cell)
      {
        currVol = cell->VolumeInPhysicalDomain();

        if (currVol < 0.0)
        {
          std::cout << "negative volume detected and reverted" << std::endl;
          currVol = -currVol;
        }

        if (cell->getDomainPlus())       plusVol += currVol;
        else                                           minusVol += currVol;
      } // end loop over domain integration cells

      if (plusVol == 0.0 || minusVol == 0.0) // no domain integration cell in one subdomain
      {
        std::cout << "element " << *currEle << " shall be bisected" << std::endl;
        std::cout << "element volume in plus domain is " << plusVol << std::endl;
        std::cout << "element volume in minus domain is " << minusVol << std::endl;
        dserror("WARNING!!! Bisected domain shall have integration cells on both interface sides!");
      } // end if no domain integration cell in one subdomain

      eleVol = plusVol + minusVol;

      if (plusVol/eleVol < critTol_)
        critElesPlus_.insert(currEle->Id());

      if (minusVol/eleVol < critTol_)
        critElesMinus_.insert(currEle->Id());
    } // end if element bisected
  } // end loop over col elements
} // end function getCritCutElements



/*------------------------------------------------------------------------------------------------*
 * find a facing flame front patch by projecton of                               winklmaier 08/10 *
 * node into boundary cells of current element                                                    *
 *----------------------------------------------------------------------------------------------- */
bool XFEM::TIMEINT::SignedDistance(
    const DRT::Node* node,
    const int elegid,
    double& dist,
    LINALG::Matrix<3,1>& normal,
    LINALG::Matrix<3,1>& proj,
    bool oldTimeStep
) const
{
  // Initialization
  const int nsd = 3;
  LINALG::Matrix<nsd,1> nodecoord(node->X()); // coordinates of this node

  Teuchos::RCP<COMBUST::InterfaceHandleCombust> ih;
  Teuchos::RCP<Epetra_Vector> phi;

  if (oldTimeStep)
  {
    ih = oldinterfacehandle_;
    phi = phin_;
  }
  else
  {
    ih = newinterfacehandle_;
    phi = phinp_;
  }

  if (intersectionStatus(discret_->gElement(elegid),oldTimeStep)==XFEM::TIMEINT::uncut_)
    dserror("call this function only for intersected elements");

  //-----------------------------------------------------------
  // compute smallest distance to the flame front for this node
  //-----------------------------------------------------------
  // smallest distance to the vertex of a flame front patch
  double vertexdist = 7777.7; // default value
  // smallest distance to the edge of a flame front patch
  double edgedist = 6666.6; // default value
  // smallest distance to flame front
  double mindist = 5555.5; // default value

  LINALG::Matrix<nsd,1> projedge(true);
  LINALG::Matrix<nsd,1> projvertex(true);

  // number of flamefront patches for this element
  const std::vector<GEO::BoundaryIntCell> patches = ih->ElementBoundaryIntCells(elegid);
  const int numpatch = patches.size();

  // loop flame front patches of this element
  for(int ipatch=0; ipatch<numpatch; ++ipatch)
  {
    // get a single patch from group of flamefront patches
    const GEO::BoundaryIntCell patch = patches[ipatch];
    // only triangles and quadrangles are allowed as flame front patches (boundary cells)
    if (!(patch.Shape() == DRT::Element::tri3 or
        patch.Shape() == DRT::Element::quad4))
    {
      dserror("invalid type of boundary integration cell for reinitialization");
    }

    // get coordinates of vertices defining flame front patch
    const LINALG::SerialDenseMatrix& patchcoord = patch.CellNodalPosXYZ();

    // compute normal vector to flame front patch
    normal.Clear();
    ComputeNormalVectorToFlameFront(patch,patchcoord,normal);

    //-----------------------------------------
    // find flame front patches facing the node
    //-----------------------------------------
    // boolean indicating if facing patch was found
    bool facenode = false;
    // distance to the facing patch
    double patchdist = 7777.7; // default value
    // check if this patch faces the node
    FindFacingPatchProjCellSpace(nodecoord,patch,patchcoord,normal,facenode,patchdist);

    // a facing patch was found
    if (facenode == true)
    {
      // overwrite smallest distance if computed patch distance is smaller
      if (fabs(patchdist) < fabs(mindist))
      {
        proj = nodecoord;
        proj.Update(patchdist,normal,1.0);
        mindist = patchdist;
      }
    }

    //-------------------------------------------------------------
    // compute smallest distance to edges of this flame front patch
    //-------------------------------------------------------------
    ComputeDistanceToEdge(nodecoord,patch,patchcoord,projedge,edgedist);

    //----------------------------------------------------------------
    // compute smallest distance to vertices of this flame front patch
    //----------------------------------------------------------------
    ComputeDistanceToPatch(nodecoord,patch,patchcoord,projvertex,vertexdist);
  }

  if (fabs(edgedist) < fabs(mindist)) // case 2a
  {
    // if G-value at the node is negative, the minimal distance has to be negative
    if ((*phi)[node->LID()] < 0.0 )
      mindist = -edgedist;
    else
      mindist = edgedist;

    proj = projedge;
  }

  if (fabs(vertexdist) < fabs(mindist))
  {
    // if the sign has been changed by mistake in ComputeDistanceToPatch(), this has to be corrected here
    if ((*phi)[node->LID()] < 0.0 )
      mindist = -vertexdist;
    else
      mindist = vertexdist;

    proj = projvertex;
  }

  if (mindist == 5555.5) // in touched-minus elements this case can happen
  {
    if (intersectionStatus(discret_->gElement(elegid),oldTimeStep) != TIMEINT::numerical_cut_)
      dserror ("computation of minimal distance failed");

    dist = (*phi)[node->LID()];
    return false;
  }
  else
  {
    dist = mindist;
    return true;
  }
}



/*------------------------------------------------------------------------------------------------*
 | private: find a facing flame front patch by projecton of node into boundary cell space         |
 |                                                                                    henke 12/09 |
 *----------------------------------------------------------------------------------------------- */
void XFEM::TIMEINT::FindFacingPatchProjCellSpace(
    const LINALG::Matrix<3,1>&       node,
    const GEO::BoundaryIntCell&      patch,
    const LINALG::SerialDenseMatrix& patchcoord,
    const LINALG::Matrix<3,1>&       normal,
    bool&                            facenode,
    double&                          patchdist
) const
{
  // indicator
  facenode = false;

  static LINALG::Matrix<2,1> eta(true);
  double alpha = 0.0;

  //-------------------------------------------------------
  // perform Newton-Raphson method to project node on patch
  //-------------------------------------------------------
  bool converged = false;
  switch(patch.Shape())
  {
  case DRT::Element::tri3:
  {
    converged = ProjectNodeOnPatch<DRT::Element::tri3>(node, patch, patchcoord, normal, eta, alpha);
    break;
  }
  case DRT::Element::quad4:
  {
    converged = ProjectNodeOnPatch<DRT::Element::quad4>(node, patch, patchcoord, normal, eta, alpha);
    break;
  }
  default:
    dserror("unknown type of boundary integration cell");
    break;
  }

  // Newton iteration converged
  //  std::cout << "Newton iteration converged in " << iter << " steps!" << std::endl;

  //----------------------------------------------------
  // check if projection lies within boundary cell space
  //----------------------------------------------------
  // remark: - tolerance has to be of same order as the tolerance that coordinates of projected nodes
  //           differ from an exact position on edges of patches (e.g. 1.0E-7 ~ 1.0E-8 -> 1.0E-6)
  //         - if this is not the case, the level set function can become tilted, since valid
  //           patches are ignored
  double TOL= 1e-6;

  switch(patch.Shape())
  {
  case DRT::Element::tri3:
  {
    // criteria for tri3 patch
    if ((eta(0) > -TOL) and (eta(0) < 1.0+TOL) and
        (eta(1) > -TOL) and (eta(1) < 1.0+TOL) and
        (1.0-eta(0)-eta(1) > -TOL) and (1.0-eta(0)-eta(1) < 1.0+TOL) and
        converged)
    {
      facenode = true;
      patchdist = alpha;
      //      std::cout << "facing patch found (tri3 patch)! coordinates eta(0): " << eta(0) << " eta(1) " << eta(1) << std::endl;
    }
    break;
  }
  case DRT::Element::quad4:
  {
    // criteria for quad4 patch
    if ((eta(0) > -1.0-TOL) and (eta(0) < 1.0+TOL) and
        (eta(1) > -1.0-TOL) and (eta(1) < 1.0+TOL) and
        converged)
    {
      facenode = true;
      patchdist = alpha;
      //      std::cout << "facing patch found (quad4 patch)!" << std::endl;
    }
    break;
  }
  default:
    dserror("unknown type of boundary integration cell");
    break;
  }
  //  if (!converged)
  //  {
  //    std::cout << "node x component " << node(0,0) << std::endl;
  //    std::cout << "node y component " << node(1,0) << std::endl;
  //    std::cout << "node z component " << node(2,0) << std::endl;
  //    std::cout << "eta1 " << eta(0) << std::endl;
  //    std::cout << "eta2 " << eta(1) << std::endl;
  //    std::cout << "alpha " << alpha << std::endl;
  //    std::cout << "patch vertices x component " << patchcoord(0,0) << " " << patchcoord(0,1) << " " << patchcoord(0,2) << std::endl;
  //    std::cout << "patch vertices y component " << patchcoord(1,0) << " " << patchcoord(1,1) << " " << patchcoord(1,2) << std::endl;
  //    std::cout << "patch vertices z component " << patchcoord(2,0) << " " << patchcoord(2,1) << " " << patchcoord(2,2) << std::endl;
  //  }

  return;
}



/*------------------------------------------------------------------------------------------------*
 | private: compute distance to edge of patch                                         henke 08/09 |
 *----------------------------------------------------------------------------------------------- */
void XFEM::TIMEINT::ComputeDistanceToEdge(
    const LINALG::Matrix<3,1>&       node,
    const GEO::BoundaryIntCell&      patch,
    const LINALG::SerialDenseMatrix& patchcoord,
    LINALG::Matrix<3,1>&             proj,
    double&                          edgedist
) const
{
  // set temporary edgedist to large value
  double edgedisttmp = edgedist;

  // number of vertices of flame front patch (3 for tri3; 4 for quad4)
  const size_t numvertices = patchcoord.N();

  // current vertex of the patch (first vertex)
  static LINALG::Matrix<3,1> vertex1(true);
  // current next vertex of the patch (second vertex)
  static LINALG::Matrix<3,1> vertex2(true);
  // distance vector from first vertex to node
  static LINALG::Matrix<3,1> vertex1tonode(true);
  // distance vector from first vertex to second vertex
  static LINALG::Matrix<3,1> vertex1tovertex2(true);

  LINALG::Matrix<3,1> tmpproj(true);

  // compute distance to all vertices of patch
  for(size_t ivert = 0; ivert<numvertices; ++ivert)
  {
    // vertex1 of flame front patch
    vertex1(0) = patchcoord(0,ivert);
    vertex1(1) = patchcoord(1,ivert);
    vertex1(2) = patchcoord(2,ivert);

    if (ivert < (numvertices-1))
    {
      vertex2(0) = patchcoord(0,ivert+1);
      vertex2(1) = patchcoord(1,ivert+1);
      vertex2(2) = patchcoord(2,ivert+1);
    }
    else if (ivert == (numvertices-1))
    {
      vertex2(0) = patchcoord(0,0);
      vertex2(1) = patchcoord(1,0);
      vertex2(2) = patchcoord(2,0);
    }

    // compute distance vector from node to current first
    vertex1tonode.Update(1.0, node, -1.0, vertex1);
    // compute distance vector from current second first vertex to current frist vertex (edge)
    vertex1tovertex2.Update(1.0, vertex2, -1.0, vertex1);
    double normvertex1tovertex2 = vertex1tovertex2.Norm2();
    // normalize vector
    vertex1tovertex2.Scale(1.0/normvertex1tovertex2);

    // scalar product of vertex1tonode and the normed vertex1tovertex2
    double lotfusspointdist = vertex1tovertex2.Dot(vertex1tonode);

    if( (lotfusspointdist >= 0.0) and (lotfusspointdist <= normvertex1tovertex2) ) // lotfusspoint on edge
    {
      tmpproj.Update(1.0,vertex1,lotfusspointdist,vertex1tovertex2);
      LINALG::Matrix<3,1> nodetolotfusspoint(true);
      nodetolotfusspoint.Update(1.0,tmpproj,-1.0,node);

      // determine length of vector from node to lot fuss point
      edgedisttmp = nodetolotfusspoint.Norm2();
      if (edgedisttmp < edgedist)
      {
        edgedist = edgedisttmp;
        proj = tmpproj;
      }
    }
  }

  return;
}


/*------------------------------------------------------------------------------------------------*
 | private: compute distance to vertex of patch                                       henke 08/09 |
 *----------------------------------------------------------------------------------------------- */
void XFEM::TIMEINT::ComputeDistanceToPatch(
    const LINALG::Matrix<3,1>&       node,
    const GEO::BoundaryIntCell&      patch,
    const LINALG::SerialDenseMatrix& patchcoord,
    LINALG::Matrix<3,1>&             proj,
    double&                          vertexdist
) const
{
  // set temporary vertexdist to large value
  double vertexdisttmp = vertexdist;

  // number of vertices of flame front patch (3 for tri3; 4 for quad4)
  const size_t numvertices = patchcoord.N();

  // current vertex of the patch
  static LINALG::Matrix<3,1> vertex(true);
  // distance vector from patch to node
  static LINALG::Matrix<3,1> dist(true);

  // compute distance to all vertices of patch
  for(size_t ivert = 0; ivert<numvertices; ++ivert)
  {
    // vertex of flame front patch
    vertex(0) = patchcoord(0,ivert);
    vertex(1) = patchcoord(1,ivert);
    vertex(2) = patchcoord(2,ivert);

    // compute distance vector from flame front to node
    dist.Update(1.0, node, -1.0, vertex);

    // compute L2-norm of distance vector
    vertexdisttmp = dist.Norm2();
    if (vertexdisttmp < vertexdist)
    {
      proj = vertex;
      vertexdist = vertexdisttmp;
    }
  }

  return;
}



/*------------------------------------------------------------------------------------------------*
 | private: compute normal vector to flame front patch                                henke 08/09 |
 *----------------------------------------------------------------------------------------------- */
void XFEM::TIMEINT::ComputeNormalVectorToFlameFront(
    const GEO::BoundaryIntCell&      patch,
    const LINALG::SerialDenseMatrix& patchcoord,
    LINALG::Matrix<3,1>&             normal
) const
{
  // first point of flame front patch
  LINALG::Matrix<3,1> point1;
  point1(0) = patchcoord(0,0);
  point1(1) = patchcoord(1,0);
  point1(2) = patchcoord(2,0);

  // second point of flame front patch
  LINALG::Matrix<3,1> point2;
  point2(0) = patchcoord(0,1);
  point2(1) = patchcoord(1,1);
  point2(2) = patchcoord(2,1);

  // first edge of flame front patch
  LINALG::Matrix<3,1> edge1;
  edge1.Update(1.0, point2, -1.0, point1);

  // third point of flame front patch
  point2(0) = patchcoord(0,2);
  point2(1) = patchcoord(1,2);
  point2(2) = patchcoord(2,2);

  // second edge of flame front patch (if patch is triangle; if not: edge 2 is secant of polygon)
  LINALG::Matrix<3,1> edge2;
  edge2.Update(1.0, point2, -1.0, point1);

  // compute normal vector of patch (cross product: edge1 x edge2)
  // remark: normal vector points into unburnt domain (G<0)
  normal(0) = (edge1(1)*edge2(2) - edge1(2)*edge2(1));
  normal(1) = (edge1(2)*edge2(0) - edge1(0)*edge2(2));
  normal(2) = (edge1(0)*edge2(1) - edge1(1)*edge2(0));
#ifdef COMBUST_2D
  normal(2) = 0.0;
#endif

  //  const Epetra_Comm& comm = scatra_.Discretization()->Comm();
  //  std::cout << "proc " << comm.MyPID() << " normal " <<  normal << std::endl;
  //  std::cout << "proc " << comm.MyPID() << " patch " <<  patchcoord << std::endl;

  // compute unit (normed) normal vector
  double norm = sqrt(normal(0)*normal(0) + normal(1)*normal(1) + normal(2)*normal(2));
  if (norm == 0.0) dserror("norm of normal vector is zero!");
  normal.Scale(1.0/norm);

#ifdef DEBUG
#ifdef COMBUST_2D
  if (!((normal(2) > 0.0-1.0E-8) and (normal(2) < 0.0+1.0E-8)))
  {
    std::cout << "z-component of normal: " << normal(2) << std::endl;
    dserror ("pseudo-3D problem not symmetric anymore!");
  }
#endif
#endif

  return;
}

