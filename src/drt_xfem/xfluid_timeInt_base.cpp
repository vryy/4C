/*!------------------------------------------------------------------------------------------------*
\file timeInt.cpp

\brief provides the basic time integration classes "TimeInt", "TimeIntStd", "TimeIntEnr"

<pre>
Maintainer: Martin Winklmaier
            winklmaier@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>
 *------------------------------------------------------------------------------------------------*/


#include "../drt_lib/drt_exporter.H"
#include "../drt_combust/combust_flamefront.H"
#include "../drt_geometry/position_array.H"
#include "../drt_xfem/xfem_fluidwizard.H"

#include "../drt_cut/cut_volumecell.H"
#include "../drt_cut/cut_elementhandle.H"
#include "../drt_cut/cut_sidehandle.H"
//#include "../drt_cut/cut_intersection.H"

#include "../drt_bele3/bele3.H"
#include "../drt_bele3/bele3_4.H"

#include "../drt_io/io.H"
#include "../drt_io/io_gmsh.H"
#include "../drt_io/io_control.H"

#include "xfluid_timeInt_base.H"


#define DEBUG_TIMINT_STD

/*------------------------------------------------------------------------------------------------*
 * basic XFEM time-integration constructor                                       winklmaier 11/11 *
 *------------------------------------------------------------------------------------------------*/
XFEM::XFLUID_TIMEINT_BASE::XFLUID_TIMEINT_BASE(
    const RCP<DRT::Discretization> discret,
    const RCP<DRT::Discretization> boundarydis,
    Teuchos::RCP<XFEM::FluidWizard>         wizard_old,       /// fluid wizard at t^n
    Teuchos::RCP<XFEM::FluidWizard>         wizard_new,       /// fluid wizard at t^(n+1)
    Teuchos::RCP<XFEM::FluidDofSet>         dofset_old,       /// fluid wizard at t^n
    Teuchos::RCP<XFEM::FluidDofSet>         dofset_new,       /// fluid wizard at t^(n+1)
    std::vector<RCP<Epetra_Vector> > oldVectors,
    const Epetra_Map& olddofcolmap,
    const Epetra_Map& newdofrowmap,
    const RCP<std::map<int,std::vector<int> > > pbcmap
) :
discret_(discret),
boundarydis_(boundarydis),
wizard_old_(wizard_old),       /// fluid wizard at t^n
wizard_new_(wizard_new),       /// fluid wizard at t^(n+1)
dofset_old_(dofset_old),       /// fluid dofset at t^n
dofset_new_(dofset_new),       /// fluid dofset at t^(n+1)
olddofcolmap_(olddofcolmap),
newdofrowmap_(newdofrowmap),
oldVectors_(oldVectors),
pbcmap_(pbcmap),
myrank_(discret_->Comm().MyPID()),
numproc_(discret_->Comm().NumProc()),
newton_max_iter_(10),
newton_tol_(1.0e-10)
{

  return;
} // end constructor



/*------------------------------------------------------------------------------------------------*
 * out of order!                                                                 winklmaier 10/11 *
 *------------------------------------------------------------------------------------------------*/
void XFEM::XFLUID_TIMEINT_BASE::compute(
    std::vector<RCP<Epetra_Vector> > newRowVectorsn,
    std::vector<RCP<Epetra_Vector> > newRowVectorsnp
)
{
  dserror("Unused function! Use a function of the derived classes");
} // end function compute



/*------------------------------------------------------------------------------------------------*
 * set the computation type with help of the iteration counter                       schott 07/12 *
 *------------------------------------------------------------------------------------------------*/
void XFEM::XFLUID_TIMEINT_BASE::type(
    int iter,
    int iterMax
)
{

  if (iter==1)
    FGIType_=FRS1FGI1_;
  else if (iterMax==1 or iter%iterMax==1)
    FGIType_=FRS1FGINot1_;
  else
    FGIType_=FRSNot1_;

  return;
} // end function type



/*------------------------------------------------------------------------------------------------*
 * assign the Epetra vectors which shall be computed to the                                       *
 * algorithms data structure                                                     winklmaier 10/11 *
 *------------------------------------------------------------------------------------------------*/
void XFEM::XFLUID_TIMEINT_BASE::handleVectors(
    std::vector<RCP<Epetra_Vector> >& newRowVectorsn
)
{
//  if (newRowVectorsn.size()!=newRowVectorsnp.size())
//    dserror("Number of state-vectors of different times are different!");

//  newVectors_.insert(newVectors_.end(),newRowVectorsnp.begin(),newRowVectorsnp.end());

  newVectors_ = newRowVectorsn;

  if (oldVectors_.size() != newVectors_.size())
    dserror("Number of state-vectors at new and old discretization are different!");
} // end function handleVectors



/*------------------------------------------------------------------------------------------------*
 * identify intersection status of an element                                    winklmaier 01/12 *
 *------------------------------------------------------------------------------------------------*/
XFEM::XFLUID_TIMEINT_BASE::intersectionType XFEM::XFLUID_TIMEINT_BASE::intersectionStatus(
    const DRT::Element* ele,
    bool oldTimeStep
) const
{
  // return status are: 0 = uncut, 1 = bisected/trisected, 2 = numerically intersected

  // initialization
  RCP<COMBUST::InterfaceHandleCombust> ih;
  RCP<Epetra_Vector> phi;

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
    return XFEM::XFLUID_TIMEINT_BASE::cut_; // bisected or trisected

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
        return XFEM::XFLUID_TIMEINT_BASE::numerical_cut_; // numerically intersected
    }
  }

  return XFEM::XFLUID_TIMEINT_BASE::uncut_;
} // end function intersectionStatus



/*------------------------------------------------------------------------------------------------*
 * identify interface side of a point in combustion                              winklmaier 10/10 *
 *------------------------------------------------------------------------------------------------*/
int XFEM::XFLUID_TIMEINT_BASE::interfaceSide(
    double phi
) const
{
  if (plusDomain(phi) == true) return 1;
  else return -1;
} // end function interfaceSide



/*------------------------------------------------------------------------------------------------*
 * identify interface side of a point in combustion                              winklmaier 10/10 *
 *------------------------------------------------------------------------------------------------*/
int XFEM::XFLUID_TIMEINT_BASE::interfaceSide(
    DRT::Element* ele,
    LINALG::Matrix<3,1> x,
    bool newTimeStep
) const
{
  const int nsd = 3; // dimension

  // required interfacehandle
  RCP<COMBUST::InterfaceHandleCombust> ih = newTimeStep ? newinterfacehandle_ : oldinterfacehandle_;
  const GEO::DomainIntCells&  domainIntCells(ih->ElementDomainIntCells(ele->Id()));

  static LINALG::Matrix<nsd,1> eta(true); // local cell coordinates
  bool pointInCell = false; // boolean whether point is in current cell

  // loop over domain integration cells
  for (GEO::DomainIntCells::const_iterator cell = domainIntCells.begin(); cell != domainIntCells.end(); ++cell)
  {
    callXToXiCoords(*cell,x,eta,pointInCell);

    if (pointInCell)
    {
      if (cell->getDomainPlus()==true) 	return 1;
      else 								return -1;
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

    if (side==true)		return 1;
    else				return -1;
  }

  dserror("point in an element should be in one of the elements domain integration cells");
  return 0;
} // end function interfaceSide


/*------------------------------------------------------------------------------------------------*
 * check if the current point x2 changed the side compared to x1                     schott 07/12 *
 *------------------------------------------------------------------------------------------------*/
bool XFEM::XFLUID_TIMEINT_BASE::changedSideSameTime(
    bool                    newTimeStep,       /// new/old timestep for both points x1 and x2
    DRT::Element*           ele1,              /// first element where x1 lies in
    LINALG::Matrix<3,1>&    x1,                /// global coordinates of point x1
    DRT::Element*           ele2,              /// second element where x2 lies in
    LINALG::Matrix<3,1>&    x2                 /// global coordinates of point x2
) const
{
  //-----------------------------------------------------------------------
  // special case of equal coordinates x1 and x2 -> no line
  LINALG::Matrix<3,1> diff(true);
  diff.Update(1.0, x1, -1.0, x2);

  if(diff.Norm2() < 1.0e-13) return false;

  //-----------------------------------------------------------------------
  // standard case of a real line between x1 and x2

  RCP<XFEM::FluidWizard> wizard = newTimeStep ? wizard_new_ : wizard_old_;

  // REMARK:
  // changing the side of a point at two times (newton steps) with coordinates x1 and x2 is done
  // via changing the cut of this trace between x1 and x2 with surrounding cutting sides
  // if there is a cut between the ray tracer and a cutting side, then the point changed the side
  // while moving from Position x1 to Position x2

  // first of all find involved elements which are possibly passed through by the ray between x1 and x2
  // second we have to find all possible cutting sides

  int mi = 0; // we assume only cutting sides with mesh 0 (meshintersection)

  bool check_allsides=false; // flag to check all sides in case of one common node is a non-row node

  std::set<int> eids;
  std::set<int> common_nodes;

  if(ele1->Id() == ele2->Id()) // line within one element
  {
    eids.insert(ele1->Id()); // just one element to check
  }
  else if( Neighbors(ele1,ele2, common_nodes) )
  {
    // get all elements adjacent to common nodes in ele1 and ele2

    // REMARK: if at least one common node is not a row node, then there is possibly one adjacent element missing,
    // but then there are also involved sides missing -> check all sides in cutdiscret, also if it is very slow!

    for(std::set<int>::iterator it=common_nodes.begin(); it!=common_nodes.end(); it++)
    {
      DRT::Node* n = discret_->gNode(*it);

      if(n->Owner() != myrank_)
      {
        check_allsides=true; // flag to check all sides in case of one common node is a non-row node
        break;
      }

      const int numele = n->NumElement();

      DRT::Element* * elements = n->Elements();

      for(int e_it=0; e_it<numele; e_it++)
      {
        DRT::Element* ele = elements[e_it];

        eids.insert(ele->Id()); cout << "\t\t add element " << ele->Id() << endl;
      }
    }
  }
  else
  {
    //dserror("please check this case: ele1 = %d and ele2 = %d are not Neighbors -> check all sides in cutdiscret", ele1->Id(), ele2->Id());

#ifdef DEBUG_TIMINT_STD
    cout << "\t\t\t all sides have to be check for changing side, ele1 = " << ele1->Id()
         << " and ele2 = " << ele2->Id()
         << " are not Neighbors" << endl;
#endif

    // check all sides ghosted in boundarydis
    check_allsides=true;
  }


  //---------------------------------------------------------------------
  std::set<int> cut_sides;

  // collect all cutting sides
  if(check_allsides)
  {
    cout << " \t\t check all sides" << endl;
    // add all sides of cut_discret to check
    for(int i=0; i< boundarydis_->NumMyColElements(); i++)
    {
      cut_sides.insert(boundarydis_->ElementColMap()->GID(i));
    }
  }
  else
  {
    //loop all found element
    for(std::set<int>::iterator e_it=eids.begin(); e_it!=eids.end(); e_it++)
    {
      DRT::Element* ele = discret_->gElement(*e_it);

      GEO::CUT::ElementHandle* eh = wizard->GetElement(ele);

      // no cutsides within this element, then no side changing possible
      if(eh == NULL)
      {
        continue; // next element
      }

      //--------------------------------------------------------
      // get involved side ids
      GEO::CUT::plain_element_set elements;
      eh->CollectElements( elements );

      // get all side-ids
      for(GEO::CUT::plain_element_set::iterator eles = elements.begin(); eles!=elements.end(); eles++)
      {
        GEO::CUT::Element* sub_ele = *eles;

        GEO::CUT::plain_facet_set facets = sub_ele->Facets();

        for(GEO::CUT::plain_facet_set::const_iterator facet_it = facets.begin(); facet_it!=facets.end(); facet_it++)
        {
          GEO::CUT::Facet* facet = *facet_it;

          GEO::CUT::Side* parent_side = facet->ParentSide();

          bool is_elements_side = sub_ele->OwnedSide(parent_side);

          if(!is_elements_side) // is a cutting side
          {
            cut_sides.insert(parent_side->Id()); cout << "\t\t add side " << parent_side->Id() << endl;
          } // !element's side
        } //facets
      }//sub elements
    }
  }

  //---------------------------------------------------------------------
  // check cut between the line and all cutting sides
  for(std::set<int>::iterator side_it= cut_sides.begin(); side_it!=cut_sides.end(); side_it++)
  {
    // get the side via sidehandle
    GEO::CUT::SideHandle* sh = wizard->GetCutSide(*side_it, mi);

    cout << "changedSideSameTime with help of side" << *side_it << endl;

    if(callSideEdgeIntersection(sh, *side_it, x1, x2))
    {
      cout << " changed side! " << endl;
      return true;
    }
  }

  return false;
}


/*------------------------------------------------------------------------------------------------*
 * check if both element are neighbors sharing at least one common node              schott 07/12 *
 *------------------------------------------------------------------------------------------------*/
bool XFEM::XFLUID_TIMEINT_BASE::Neighbors(
    DRT::Element* ele1,
    DRT::Element* ele2,
    std::set<int> & common_nodes
    ) const
{

  bool is_neighbor = false;

  const int numnode1 =  ele1->NumNode();
  const int numnode2 =  ele2->NumNode();

  const int* nodeids1 = ele1->NodeIds();
  const int* nodeids2 = ele2->NodeIds();

  std::vector<int> nodeids1_vec;
  std::vector<int> nodeids2_vec;

  for(int i=0; i< numnode1; i++)
    nodeids1_vec.push_back(nodeids1[i]);

  for(int i=0; i< numnode2; i++)
    nodeids2_vec.push_back(nodeids2[i]);

  sort(nodeids1_vec.begin(), nodeids1_vec.end());
  sort(nodeids2_vec.begin(), nodeids2_vec.end());

  //------------------------------------------------------
  // find common nodes
  int index1=0;
  int index2=0;
  while(true)
  {
    // if at least one vector of nodes is checked
    if(index1 >= numnode1 or index2 >= numnode2) break;

    if(nodeids1_vec[index1] < nodeids2_vec[index2]) index1++;
    else if(nodeids1_vec[index1] > nodeids2_vec[index2]) index2++;
    else // nodeids1_vec[index1] == nodeids2_vec[index2
    {
      common_nodes.insert(nodeids1_vec[index1]);
    cout << "common node between element " << ele1->Id() << " and element " << ele2->Id() << " is " << nodeids1_vec[index1] << endl;
      is_neighbor = true; // true if at least one node is common

      index1++;
      index2++;
    }
  }

  return is_neighbor;
}


/*------------------------------------------------------------------------------------------------*
 * check if edge between x1 and x2 cuts the side                                     schott 07/12 *
 *------------------------------------------------------------------------------------------------*/
bool XFEM::XFLUID_TIMEINT_BASE::callSideEdgeIntersection(
    GEO::CUT::SideHandle*         sh,   /// side handle
    int                           sid,  /// side id
    LINALG::Matrix<3,1>&          x1,   /// coordinates of edge's start point
    LINALG::Matrix<3,1>&          x2    /// coordinates of edge's end point
    ) const
{

  switch (sh->Shape())
  {
  case DRT::Element::tri3:
  {
    return callSideEdgeIntersectionT<DRT::Element::tri3>(sh,sid,x1,x2);
  }
  case DRT::Element::quad4:
  {
    return callSideEdgeIntersectionT<DRT::Element::quad4>(sh,sid,x1,x2);
  }
  case DRT::Element::quad8:
  {
    return callSideEdgeIntersectionT<DRT::Element::quad8>(sh,sid,x1,x2);
  }
  case DRT::Element::quad9:
  {
    return callSideEdgeIntersectionT<DRT::Element::quad9>(sh,sid,x1,x2);
  }
  default:
  {
    dserror("unknown side shape");
    break;
  }
  }

  return false;
}

/*------------------------------------------------------------------------------------------------*
 * check if edge between x1 and x2 cuts the side (templated)                         schott 07/12 *
 *------------------------------------------------------------------------------------------------*/
template < DRT::Element::DiscretizationType sidetype>
bool XFEM::XFLUID_TIMEINT_BASE::callSideEdgeIntersectionT(
    GEO::CUT::SideHandle*         sh,   /// side handle
    int                           sid,  /// side id
    LINALG::Matrix<3,1>&          x1,   /// coordinates of edge's start point
    LINALG::Matrix<3,1>&          x2    /// coordinates of edge's end point
) const
{

//  cout << "\t\t check changed side via intersection between side " << sid
//       << " and line " << x1 << " and " << x2 << endl;

  const int nsd = 3;
  const int numNodesSurface = DRT::UTILS::DisTypeToNumNodePerEle<sidetype>::numNodePerElement;

  LINALG::Matrix<nsd,2> xyze_lineElement(true);

  for(int i=0; i<nsd; i++)
  {
    xyze_lineElement(i,0) = x1(i);
    xyze_lineElement(i,1) = x2(i);
  }

  Epetra_SerialDenseMatrix xyze_side;
  sh->Coordinates( xyze_side );

  LINALG::Matrix<nsd,numNodesSurface> xyze_surfaceElement(xyze_side);

  LINALG::Matrix<3,1> xsi(true);

  //cout << "\t\tcompute intersection between side and line KERNEL::ComputeIntersection" << endl;

  GEO::CUT::KERNEL::DebugComputeIntersection<DRT::Element::line2, sidetype> ci( xsi );
  //GEO::CUT::KERNEL::ComputeIntersection<DRT::Element::line2, sidetype> ci( xsi );
  if ( ci( xyze_surfaceElement, xyze_lineElement ) )
  {
      return true;
  }

  return false;
}




/*------------------------------------------------------------------------------------------------*
 * call the computation of local coordinates for an element                      winklmaier 11/11 *
 *------------------------------------------------------------------------------------------------*/
void XFEM::XFLUID_TIMEINT_BASE::callXToXiCoords(
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
 * call the computation of local coordinates for an integration cell                 schott 07/12 *
 *------------------------------------------------------------------------------------------------*/
void XFEM::XFLUID_TIMEINT_BASE::callXToXiCoords(
    const DRT::Element*    ele,           /// pointer to element
    LINALG::Matrix<3,1>&   x,             /// global coordinates of point
    LINALG::Matrix<3,1>&   xi,            /// determined local coordinates w.r.t ele
    bool&                  pointInDomain  /// lies point in element ?
) const
{
  LINALG::SerialDenseMatrix xyz(3,ele->NumNode(),true);
  GEO::fillInitialPositionArray(ele, xyz);
  callXToXiCoords(xyz,ele->Shape(),x,xi,pointInDomain);
} // end function callXToXiCoords



/*------------------------------------------------------------------------------------------------*
 * call the computation of local coordinates for a polytop                                        *
 * with corners given by the coordinates                                             schott 07/12 *
 *------------------------------------------------------------------------------------------------*/
void XFEM::XFLUID_TIMEINT_BASE::callXToXiCoords(
    LINALG::SerialDenseMatrix&        nodecoords,    /// node coordinates of element
    DRT::Element::DiscretizationType  DISTYPE,       /// discretization type
    LINALG::Matrix<3,1>&              x,             /// global coordinates of point
    LINALG::Matrix<3,1>&              xi,            /// determined local coordinates w.r.t ele
    bool&                             pointInDomain  /// lies point in element ?
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
  default: dserror("add your 3D distype and the according transformation!");
  } // end switch

  return;
} // end function callXToXiCoords



/*------------------------------------------------------------------------------------------------*
 * add adjacent elements for a periodic boundary node                            winklmaier 05/11 *
 *------------------------------------------------------------------------------------------------*/
void XFEM::XFLUID_TIMEINT_BASE::addPBCelements(
    const DRT::Node* node,
    std::vector<const DRT::Element*>&  eles
) const
{
  dserror("what to do in addPBCelements?");

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



void XFEM::XFLUID_TIMEINT_BASE::findPBCNode(
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
 * reset a special state with another state in the data class                        schott 07/12 *
 *------------------------------------------------------------------------------------------------*/
void XFEM::XFLUID_TIMEINT_BASE::resetState(
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

  return;
} // end function resetState



/*------------------------------------------------------------------------------------------------*
 * clear the data of all nodes having a special state                                schott 07/12 *
 *------------------------------------------------------------------------------------------------*/
void XFEM::XFLUID_TIMEINT_BASE::clearState(
    TimeIntData::state state         /// state of time int to clear
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

  return;
} // end function clear state



#ifdef PARALLEL
/*------------------------------------------------------------------------------------------------*
 * basic function sending data to dest and receiving data from source            winklmaier 10/10 *
 *------------------------------------------------------------------------------------------------*/
void XFEM::XFLUID_TIMEINT_BASE::sendData(
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
  cout << "--- sending "<< lengthSend[0] << " bytes: from proc " << myrank_ << " to proc " << dest << endl;
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
  cout << "--- receiving "<< lengthRecv[0] << " bytes: to proc " << myrank_ << " from proc " << source << endl;
#endif
} // end sendData



/*------------------------------------------------------------------------------------------------*
 * packing a node for parallel communication only with the basic nodal data                       *
 * without an underlying discretization fitting to the node's new prozessor      winklmaier 10/10 *
 *------------------------------------------------------------------------------------------------*/
void XFEM::XFLUID_TIMEINT_BASE::packNode(
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
void XFEM::XFLUID_TIMEINT_BASE::unpackNode(
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
 * basic XFEM time-integration constructor for standard degrees of freedom           schott 07/12 *
 *------------------------------------------------------------------------------------------------*/
XFEM::XFLUID_STD::XFLUID_STD(
    XFEM::XFLUID_TIMEINT_BASE&                                       timeInt,            /// time integration base class object
    const std::map<int, std::vector<INPAR::XFEM::XFluidTimeInt> >&   reconstr_method,    /// reconstruction map for nodes and its dofsets
    INPAR::XFEM::XFluidTimeInt&                                      timeIntType,        /// type of time integration
    const RCP<Epetra_Vector>                                         veln,               /// velocity at time t^n
    const double&                                                    dt,                 /// time step size
    bool                                                             initialize          /// is initialization?
) :
XFEM::XFLUID_TIMEINT_BASE::XFLUID_TIMEINT_BASE(timeInt),
timeIntType_(timeIntType),
veln_(veln),
dt_(dt)
{
  if (initialize)
  {
    const int nsd = 3; // dimension
    timeIntData_ = Teuchos::rcp(new std::vector<TimeIntData>); // vector containing all data used for computation

    /*------------------------*
     * Initialization         *
     *------------------------*/
//    for (int leleid=0; leleid<discret_->NumMyColElements(); leleid++)  // loop over processor nodes
//    {
//      DRT::Element* iele = discret_->lColElement(leleid);
//      if (intersectionStatus(iele)!=XFLUID_TIMEINT_BASE::uncut_) // element cut
//      {
//        const int* nodeGids = iele->NodeIds(); // node gids
//        for (int inode=0;inode<iele->NumNode();inode++) // loop over element nodes
//        {
//          if(olddofcolmap_.MyGID(nodeGids[inode]));
//          oldEnrNodes_.insert(nodeGids[inode]);
//        } // end loop over element nodes
//      } // end if element cut
//    } // end loop over processor nodes

    LINALG::Matrix<nsd,1> dummyStartpoint; // dummy startpoint for comparison
    for (int i=0;i<nsd;i++) dummyStartpoint(i) = 777.777;

    //--------------------------------------------------------------------------------------
    // fill timeIntData_ structure with the data for the nodes which are marked for SEMILAGRANGEAN reconstruction

#ifdef DEBUG_TIMINT_STD
    cout << "\n\tFill timeIntData_ for SEMI-LAGRANGIAN algorithm" << endl;
#endif
    // loop over processor nodes
    for (int lnodeid=0; lnodeid<discret_->NumMyRowNodes(); lnodeid++)
    {
      DRT::Node* currnode = discret_->lRowNode(lnodeid); // current analysed node

      std::map<int, std::vector<INPAR::XFEM::XFluidTimeInt> >::const_iterator it = reconstr_method.find(currnode->Id());
      if(it != reconstr_method.end())
      {
        // reconstruction methods just for marked dofsets
        for(size_t i=0; i<(it->second).size(); i++)
        {
          if((it->second)[i] == INPAR::XFEM::Xf_TimeInt_SemiLagrange )
          {
#ifdef DEBUG_TIMINT_STD
            cout << "\n\tcall Semilagrange Algorithm for node " << currnode->Id() << " and dofset " << i;
#endif

            // constructor for standard computation
            timeIntData_->push_back(TimeIntData(
                *currnode,                                                                             // node
                i,                                                                                     // nds (nodal dofset) at new timestep
                LINALG::Matrix<nsd,1>(true),                                                           // velocity
                vector<LINALG::Matrix<nsd,nsd> >(oldVectors_.size(),LINALG::Matrix<nsd,nsd>(true)),    // vel deriv
                vector<LINALG::Matrix<1,nsd> >(oldVectors_.size(),LINALG::Matrix<1,nsd>(true)),        // pres deriv
                dummyStartpoint,                                                                       // ...
                1,                                                                                     // searchedProcs
                0,                                                                                     // counter
                vector<int>(1,-1),                                                                     // startGid
                vector<int>(1,-1),                                                                     // startOwner
                INFINITY,                                                                              // minimal distance
                TimeIntData::predictor_)); // data created for the node


          } // semi-lagrangian algo
        } // nodaldofsets
      } // some dofsets has to be reconstructed
    } // end loop over processor nodes

    //--------------------------------------------------------------------------------------
    // compute initial points on the right side of the interface to start finding the lagrangian origin
    startpoints();

    //--------------------------------------------------------------------------------------

    // test loop if all initial startpoints have been computed
    for (vector<TimeIntData>::iterator data=timeIntData_->begin();
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
 * out of order!                                                                     schott 07/12 *
 *------------------------------------------------------------------------------------------------*/
void XFEM::XFLUID_STD::compute(
    vector<RCP<Epetra_Vector> >& newRowVectorsn
)
{
  dserror("Unused function! Use a function of the derived classes");
}



///*------------------------------------------------------------------------------------------------*
// * initialize data when called in a new FGI                                      winklmaier 10/11 *
// *------------------------------------------------------------------------------------------------*/
//void XFEM::XFLUID_STD::importNewFGIData(
//    const RCP<DRT::Discretization> discret,
//    const RCP<XFEM::DofManager> newdofman,
//    const RCP<COMBUST::FlameFront> flamefront,
//    const Epetra_Map& newdofrowmap,
//    const map<DofKey, DofGID>& newNodalDofRowDistrib)
//{
//  discret_ = discret;
//  newdofman_ = newdofman;
//  phinpi_ = phinp_;
//  phinp_ = flamefront->Phinp();
//  flamefront_ = flamefront;
//  newinterfacehandle_ = flamefront->InterfaceHandle();
//  newdofrowmap_ = newdofrowmap;
//  newNodalDofRowDistrib_ = newNodalDofRowDistrib;
//  return;
//}



/*------------------------------------------------------------------------------------------------*
 * identify an element containing a point and additional data                        schott 07/12 *
 *------------------------------------------------------------------------------------------------*/
void XFEM::XFLUID_STD::elementSearch(
    DRT::Element*&        ele,   /// pointer to element if point lies in a found element
    LINALG::Matrix<3,1>&  x,     /// global coordiantes of point
    LINALG::Matrix<3,1>&  xi,    /// determined local coordinates w.r.t ele
    bool&                 found  /// is element found?
) const
{

  // REMARK: if ele!= NULL, then check that element first, before loop all row elements

  int startid;                   // local row element id
  if (ele==NULL) startid = 0;    // start with first local row element
  else           startid = -1;   // pseudo-id so that id+1 will be 0 (for additional element check, if ele!=NULL)

  DRT::Element* currele = NULL;  // current element

  //loop over elements
  for (int ieleid = startid;ieleid<discret_->NumMyRowElements();ieleid++)
  {
    // if ele != NULL additional check
    // first it should be checked if it is fitting
    if (ieleid == -1)
    {
      currele = ele;
      ele = NULL; // reset ele, ele will be set if an element is found finally
    }
    else
      currele = discret_->lRowElement(ieleid);

    //    cout << "currently analysed element" << *currele;
    //    cout << "startpoint approximation: " << x;
    callXToXiCoords(currele,x,xi,found);
    //    cout << "in local xi-coordinates: " << xi;

    if (found)
    {
      ele = currele;
      //cout << "coordinates found in element " << ele->Id() << endl;
      break;
    }
  } // end loop over processor elements

  return;
} // end elementSearch



/*------------------------------------------------------------------------------------------------*
 * interpolate velocity and derivatives for a point in an element                    schott 06/12 *
 *------------------------------------------------------------------------------------------------*/
void XFEM::XFLUID_STD::getGPValues(
    DRT::Element*                 ele,            /// pointer to element
    LINALG::Matrix<3,1>&          xi,             /// local coordinates of point w.r.t element
    std::vector<int>&             nds,            /// nodal dofset of point for elemental nodes
    bool                          step_np,        /// computation w.r.t. old or new interface position
    LINALG::Matrix<3,1>&          vel,            /// determine velocity at point
    LINALG::Matrix<3,3>&          vel_deriv,      /// determine velocity derivatives at point
    bool                          compute_deriv   /// shall derivatives be computed?
) const
{

  switch (ele->Shape())
  {
  case DRT::Element::hex8:
    getGPValuesT<DRT::Element::hex8>(ele,xi,nds,step_np,vel,vel_deriv,compute_deriv);
    break;
  case DRT::Element::hex20:
    getGPValuesT<DRT::Element::hex20>(ele,xi,nds,step_np,vel,vel_deriv,compute_deriv);
    break;
  case DRT::Element::tet4:
    getGPValuesT<DRT::Element::tet4>(ele,xi,nds,step_np,vel,vel_deriv,compute_deriv);
    break;
  default: dserror("add your 3D distype here!");
  } // end switch

  return;
} // end function getGPValues






/*------------------------------------------------------------------------------------------------*
 * Compute starting values for nodes marked for semi-lagrangian algo                 schott 07/12 *
 *------------------------------------------------------------------------------------------------*/
void XFEM::XFLUID_STD::startpoints()
{

#ifdef DEBUG_TIMINT_STD
  cout << "\n\t compute initial start points for finding lagrangean origin" << endl;
#endif

  // REMARK: we do not need parallel communication for start point values because
  // structural surface is ghosted on all procs

  // loop over nodes which changed interface side
  for (vector<TimeIntData>::iterator data=timeIntData_->begin();
      data!=timeIntData_->end(); data++)
  {
    if (data->state_==TimeIntData::basicStd_) // correct state
    {
      // find a start approximation for semilagrange backtracking
      // this start approximation can be found by tracking back the projection of current node
      // along structural movement
      ProjectAndTrackback(*data);

    } // end if correct state
  } // end loop over nodes which changed interface side


  return;
} // end startpoints()




/*------------------------------------------------------------------------------------------------*
 * Project the current point onto the structural surface, track it back along structural         *
 * movement and project it back into the fluid domain                                schott 06/12 *
 *------------------------------------------------------------------------------------------------*/
void XFEM::XFLUID_STD::ProjectAndTrackback( TimeIntData& data)
{

#ifdef DEBUG_TIMINT_STD
  cout << "\n\t ProjectAndTrackback for node " << data.node_.Id() << endl;
#endif

  const int nsd = 3;
  LINALG::Matrix<nsd,1> newNodeCoords(data.node_.X()); // coords of endpoint of Lagrangian characteristics

  // determine the smallest distance of node at t^(n+1) to the current side elements

  GEO::CUT::Node* n_new = wizard_new_->GetNode(data.node_.Id());

  if(n_new == NULL) dserror("node %d does not carry a node handle (semilagrange reasonable?)", data.node_.Id());

  //--------------------------------------------------------
  // get involved side ids for projection and distance computation
  //--------------------------------------------------------

  // set of side-ids involved in cutting the current connection of volumecells at t^(n+1)
  std::map<int, std::vector<GEO::CUT::BoundaryCell*> >  bcells_new;

  const std::vector<std::set<GEO::CUT::plain_volumecell_set, GEO::CUT::Cmp > >& dof_cellsets_new = n_new->DofCellSets();

  // the std-set
  const std::set<GEO::CUT::plain_volumecell_set, GEO::CUT::Cmp >& cell_set_new = dof_cellsets_new[data.nds_np_];

  // get all side-ids w.r.t to all volumecells contained in current new set around the current node
  for(std::set<GEO::CUT::plain_volumecell_set>::const_iterator adj_eles = cell_set_new.begin(); adj_eles!=cell_set_new.end(); adj_eles++)
  {
    const GEO::CUT::plain_volumecell_set ele_vc = *adj_eles;

    for(GEO::CUT::plain_volumecell_set::const_iterator vcs=ele_vc.begin(); vcs!=ele_vc.end(); vcs++)
    {
      GEO::CUT::VolumeCell* vc = *vcs;

      // get sides involved in creation boundary cells (map<sideId,bcells>)
      vc->GetBoundaryCells(bcells_new);
    }
  }

  //--------------------------------------------------------
  // distance computation w.r.t surface involved in cutting the adjacent elements(n_new to sides/edges/nodes)
  //--------------------------------------------------------

  //------------------------------------
  // find all involved nodes and sides (edges) for computing distance to node

  std::set<int> points;
  std::set<int> sides;

  // loop bcs and extract sides and nodes
  for(std::map<int, std::vector<GEO::CUT::BoundaryCell*> >::iterator bcells_it= bcells_new.begin();
      bcells_it!=bcells_new.end();
      bcells_it++)
  {
    int sid=bcells_it->first;
    sides.insert(sid);

    // get the side
    DRT::Element* side = boundarydis_->gElement(sid);
    int numnode = side->NumNode();
    const int * side_nodes = side->NodeIds();

    // get the nodes
    for(int i=0; i< numnode; i++)
    {
      points.insert(side_nodes[i]);
    }
  }

  if(points.size() == 0 and sides.size() == 0) dserror("there are no cutting sides around node %d, why Semilagrange here?", data.node_.Id());

  //-------------------------------------
  // initialize data holding information about minimal distance

  // smallest distance
  LINALG::Matrix<3,1> proj_x_np(true);                ///< projected point at t^(n+1)
  LINALG::Matrix<3,1> proj_x_n(true);                 ///< projected point at t^n (tracked back along structural movement)
  LINALG::Matrix<3,1> start_point(true);              ///< final start point for SemiLagrange algo
  int proj_sid = -1;                                  ///< id of side that contains the projected point
  std::map<vector<int>, vector<int> >    proj_lineid;         ///< map< sorted nids, global side IDs >
  // smallest distance w.r.t side
  LINALG::Matrix<2,1> proj_xi_side(true);             ///< local coordinates of projected point if projection w.r.t side
  // smallest distance w.r.t line
  std::map<vector<int>, std::vector<double> > proj_xi_line;        ///< map<sorted nids,local line coordinates w.r.t lines of different sides>
  std::map<vector<int>, vector<int> >    proj_nid_line;       ///< map<sorted nids, sideids>
  // smallest distance w.r.t point
  int proj_nid_np = -1;                               ///< nid of projected point if projection w.r.t point (structural node)

  double min_dist = INFINITY;                         ///< minimal distance

  data.proj_ = TimeIntData::failed_;                  ///< projection method yielding minimal distance

  //----------------------------------------------------------------------
  // check distance to sides and its lines
  //----------------------------------------------------------------------
  for(std::set<int>::iterator it= sides.begin(); it!=sides.end(); it++)
  {

    //-------------------------------------------------

    DRT::Element* side = boundarydis_->gElement(*it);

#ifdef DEBUG_TIMINT_STD
      cout << "\t\t\tcall projection of point for side element " << side->Id() << endl;
#endif

    CallProjectOnSide(side,
                      newNodeCoords,
                      proj_sid,
                      min_dist,
                      proj_x_np,
                      proj_xi_side,
                      data);


    //------------------------------------
    // edges are checked twice (w.r.t both sides)
    //------------------------------------

    // loop all edges of current side
    std::vector<Teuchos::RCP<DRT::Element> > lines = side->Lines();

    int line_count = 0;

    for(std::vector<Teuchos::RCP<DRT::Element> >::iterator line_it = lines.begin(); line_it!= lines.end(); line_it++)
    {

#ifdef DEBUG_TIMINT_STD
      cout << "\t\t\tcall projection of point to line " << line_count << " for side element " << side->Id() << endl;
#endif

      CallProjectOnLine(side,
                        &(**line_it),
                        line_count,           ///< local line id w.r.t side element
                        newNodeCoords,        ///< node coordinates of point that has to be projected
                        min_dist,             ///< minimal distance, potentially updated
                        proj_x_np,            ///< projection of point on this side
                        proj_xi_line,         ///< map<side ID, local coordinates of projection of point w.r.t to this line>
                        proj_lineid,          ///< map<side ID, local line id>
                        proj_nid_line,        ///< map<side ID, vec<line Ids>>
                        proj_sid,             ///< id of side that contains the projected point
                        data                  ///< reference to data
                         );



      line_count++;
    } // loop lines of side
  } // loop sides


  //------------------------------------
  // check minimal distance w.r.t points (structural nodes)
  //------------------------------------
  // loop all points
  for(std::set<int>::iterator it = points.begin(); it!=points.end(); it++)
  {
    // compute distance to following point
    DRT::Node* node = boundarydis_->gNode(*it);

    CallProjectOnPoint(
                       node,                 ///< pointer to node
                       newNodeCoords,        ///< node coordinates of point that has to be projected
                       min_dist,             ///< minimal distance, potentially updated
                       proj_nid_np,          ///< nid id of projected point on surface
                       proj_x_np,            ///< projection of point on this side
                       proj_sid,             ///< id of side that contains the projected point
                       proj_xi_line,         ///< map<side ID, local coordinates of projection of point w.r.t to this line>
                       proj_lineid,          ///< map<side ID, local line id>
                       proj_nid_line,        ///< map<side ID, vec<line Ids>>
                       data                  ///< reference to data
                       );



  } // loop points

  if(data.proj_ == TimeIntData::failed_) dserror("projection of node %d not successful!", n_new->Id());


  // track back the projected points on structural surface from t^(n+1)->t^n
  if(data.proj_ == TimeIntData::onSide_)
  {
    DRT::Element* side = boundarydis_->gElement(proj_sid);

    // side geometry at initial state t^0
    const int numnodes = side->NumNode();
    DRT::Node ** nodes = side->Nodes();
    Epetra_SerialDenseMatrix side_xyze( 3, numnodes );
    for ( int i=0; i<numnodes; ++i )
    {
      const double * x = nodes[i]->X();
      std::copy( x, x+3, &side_xyze( 0, i ) );
    }

    DRT::Element::LocationArray cutla( 1 );
    side->LocationVector(*boundarydis_,cutla,false);

    ComputeStartPoint_Side(side, side_xyze, cutla[0].lm_, proj_xi_side, min_dist, proj_x_n, start_point);


  }
  else if(data.proj_ == TimeIntData::onLine_)
  {
    // check if both lines are identical
    for(std::map<std::vector<int>, std::vector<int> >::iterator line=proj_nid_line.begin(); line!=proj_nid_line.end(); line++)
    {
      std::vector<int>& sides = line->second;

      std::vector<int>&    local_lineIds = proj_lineid.find(line->first)->second;
      std::vector<double>& local_xi_line = proj_xi_line.find(line->first)->second;

      if(sides.size()!= 2) dserror("there must be two sides adjacent to line, but there are %d sides", sides.size());

      //---------------------------------------------------------
      // side 1
      DRT::Element* side_1 = boundarydis_->gElement(sides[0]);

      // side geometry at initial state t^0
      const int numnodes_1 = side_1->NumNode();
      DRT::Node ** nodes_1 = side_1->Nodes();
      Epetra_SerialDenseMatrix side_xyze_1( 3, numnodes_1 );
      for ( int i=0; i<numnodes_1; ++i )
      {
        const double * x = nodes_1[i]->X();
        std::copy( x, x+3, &side_xyze_1( 0, i ) );
      }

      DRT::Element::LocationArray cutla_1( 1 );
      side_1->LocationVector(*boundarydis_,cutla_1,false);

      //---------------------------------------------------------
      // side 2
      DRT::Element* side_2 = boundarydis_->gElement(sides[1]);

      // side geometry at initial state t^0
      const int numnodes_2 = side_2->NumNode();
      DRT::Node ** nodes_2 = side_2->Nodes();
      Epetra_SerialDenseMatrix side_xyze_2( 3, numnodes_2 );
      for ( int i=0; i<numnodes_2; ++i )
      {
        const double * x = nodes_2[i]->X();
        std::copy( x, x+3, &side_xyze_2( 0, i ) );
      }

      DRT::Element::LocationArray cutla_2( 1 );
      side_2->LocationVector(*boundarydis_,cutla_2,false);

      //---------------------------------------------------------

      // line geometry at initial state t^0

      std::vector<Teuchos::RCP<DRT::Element> > lines = side_1->Lines();

      RCP<DRT::Element> line_ele = lines[local_lineIds[0]];
//
//      // line geometry
//      const int numnodes = line->NumNode();
//      DRT::Node ** nodes = line->Nodes();
//      Epetra_SerialDenseMatrix line_xyze( 3, numnodes );
//      for ( int i=0; i<numnodes; ++i )
//      {
//        const double * x = nodes[i]->X();
//        std::copy( x, x+3, &line_xyze( 0, i ) );
//      }
//
//
//      DRT::Element::LocationArray cutla( 1 );
//      line->LocationVector(*boundarydis_,cutla,false);

      call_get_projxn_Line(side_1, &*line_ele, proj_x_n, local_xi_line[0]);


      ComputeStartPoint_Line(
          side_1,          ///< pointer to side element
          side_xyze_1,     ///< side's node coordinates
          side_2,          ///< pointer to side element
          side_xyze_2,     ///< side's node coordinates
          cutla_1[0].lm_,  ///< local map
          cutla_2[0].lm_,  ///< local map
          min_dist,        ///< distance from point to its projection
          proj_x_n,        ///< projected point at t^n
          start_point      ///< final start point
      );

    }
  }
  else if(data.proj_ == TimeIntData::onPoint_)
  {
    dserror("implement startvalue w.r.t point");
  }
  else dserror("unknown projection method for tracking back the structural surface point");



  // track back the projected points from t^(n+1) -> t^n




  if(data.proj_ != TimeIntData::failed_)
  {
    // fill data
    //  data.startpoint_.Update(1.0,newNodeCoords,-dt_,nodevel);
    data.startpoint_.Update(1.0,start_point,0.0);
    data.initialpoint_.Update(1.0,start_point,0.0);
    data.dMin_ = min_dist;
    data.startGid_.clear();
    data.startGid_.push_back(-1);
    data.startOwner_.clear();
    data.startOwner_.push_back(-1);
  }













  // GMSH debug output

//
//  std::string filename="start_value_Semilagrange";
  int myrank = 0; //discret.Comm().MyPID();

  if(myrank==0) std::cout << "\n\t ... writing Gmsh output...\n" << std::flush;

  int step = 0;

  const int step_diff = 100;
  bool screen_out = true; //xfluid_.gmsh_debug_out_screen_;

   // output for Element and Node IDs
   std::ostringstream filename_base_startvalues;
   filename_base_startvalues << "TIMINT"<< "_startvalues_" <<  n_new->Id()<< "_";
   const std::string filename_startvalues = IO::GMSH::GetNewFileNameAndDeleteOldFiles(filename_base_startvalues.str(), step, step_diff, screen_out,myrank);
   if (true) cout << endl;
   std::ofstream gmshfilecontent_startvalues(filename_startvalues.c_str());
   gmshfilecontent_startvalues.setf(std::ios::scientific,std::ios::floatfield);
   gmshfilecontent_startvalues.precision(16);

   gmshfilecontent_startvalues         << "View \"" << "SL " << "node "   << n_new->Id() << "\" {\n";

   {
   // 1. line ( x_np -> x_np projected on side )
   gmshfilecontent_startvalues<<"SL(" <<newNodeCoords(0)<<","<<newNodeCoords(1)<<","<<newNodeCoords(2)<<","
                                      <<proj_x_np(0)<<","<<proj_x_np(1)<<","<<proj_x_np(2)<<")"<<std::endl;
   gmshfilecontent_startvalues<<"{" << 0 << "," << 1 << "};"<<std::endl;
   // 2. line ( x_np projected on side -> x_n projected on side )
   gmshfilecontent_startvalues<<"SL(" <<proj_x_np(0)<<","<<proj_x_np(1)<<","<<proj_x_np(2)<<","
                                      <<proj_x_n(0)<<","<<proj_x_n(1)<<","<<proj_x_n(2)<<")"<<std::endl;
   gmshfilecontent_startvalues<<"{" << 1 << "," << 2 << "};"<<std::endl;
   // 3. line ( x_n projected on side -> start_point)
   gmshfilecontent_startvalues<<"SL(" <<proj_x_n(0)<<","<<proj_x_n(1)<<","<<proj_x_n(2)<<","
                                      <<start_point(0)<<","<<start_point(1)<<","<<start_point(2)<<")"<<std::endl;
   gmshfilecontent_startvalues<<"{" << 2 << "," << 3 << "};"<<std::endl;
   }

   gmshfilecontent_startvalues   << "};\n";

//  int pointno=1,point_begin,point_end,lineno=1;
//  for(plain_facet_set::const_iterator i=facete.begin();i!=facete.end();i++)
//  {
//     point_begin = pointno;
//     Facet *fe = *i;
////the side id of the first facet is used as the file name for every volumecell
//     if(i==facete.begin())
//     {
//       static int sideno = 0;
//       sideno++;
//       std::stringstream out;
//       out <<"side"<<sideno<<".pos";
//       filename = out.str();
//       file.open(filename.c_str());
//     }
//     const std::vector<std::vector<double> > corners = fe->CornerPointsLocal(ParentElement());
//     for(std::vector<std::vector<double> >::const_iterator k=corners.begin();k!=corners.end();k++)
//     {
//       const std::vector<double> coords = *k;
//       file<<"Point("<<pointno<<")={"<<coords[0]<<","<<coords[1]<<","<<coords[2]<<","<<"1"<<"};"<<std::endl;
//       pointno++;
//     }
//     point_end = pointno;
//     for(int i=point_begin;i!=point_end;i++)
//     {
//       if(i!=point_end-1)
//         file<<"Line("<<lineno<<")={"<<i<<","<<i+1<<"};"<<std::endl;
//       else
//         file<<"Line("<<lineno<<")={"<<i<<","<<point_begin<<"};"<<std::endl;
//       lineno++;
//     }
//  }
//
//  for(unsigned i=0;i<gauspts.size();i++)
//  {
//     file<<"Point("<<pointno<<")={"<<gauspts[i][0]<<","<<gauspts[i][1]<<","<<gauspts[i][2]<<","<<"1"<<"};"<<std::endl;
//     pointno++;
//  }
//  file.close();




}


void XFEM::XFLUID_STD::ComputeStartPoint_Side(
    DRT::Element*                side,           ///< pointer to side element
    Epetra_SerialDenseMatrix &   side_xyze,      ///< side's node coordinates
    const std::vector<int>&      lm,             ///< local map
    LINALG::Matrix<2,1> &        xi_side,        ///< local coordinates of projected point w.r.t side
    double &                     dist,           ///< distance from point to its projection
    LINALG::Matrix<3,1>&         proj_x_n,       ///< projected point at t^n
    LINALG::Matrix<3,1>&         start_point     ///< final start point
)
{

  LINALG::Matrix<3,1> normal(true);

  callgetNormalSide_tn(side, normal, side_xyze, lm, proj_x_n, xi_side);

  // map point into fluid domain along normal vector
  start_point.Update(1.0, proj_x_n, dist, normal);

  return;
}

void XFEM::XFLUID_STD::ComputeStartPoint_Line(
    DRT::Element*                side1,           ///< pointer to side element
    Epetra_SerialDenseMatrix &   side1_xyze,      ///< side's node coordinates
    DRT::Element*                side2,           ///< pointer to side element
    Epetra_SerialDenseMatrix &   side2_xyze,      ///< side's node coordinates
    const std::vector<int>&           lm1,             ///< local map
    const std::vector<int>&           lm2,             ///< local map
    double &                     dist,           ///< distance from point to its projection
    LINALG::Matrix<3,1>&         proj_x_n,       ///< projected point at t^n
    LINALG::Matrix<3,1>&         start_point     ///< final start point
)
{

  LINALG::Matrix<3,1> normal_avg(true); // averaged normal vector
  LINALG::Matrix<3,1> normal1(true);
  LINALG::Matrix<3,1> normal2(true);

  LINALG::Matrix<2,1> xi_side1(true);
  LINALG::Matrix<2,1> xi_side2(true);


  LINALG::Matrix<3,1> xi_1_avg(true);
  LINALG::Matrix<3,1> xi_2_avg(true);


  for(int i=0; i< side1->NumNode(); i++)
    xi_1_avg = DRT::UTILS::getNodeCoordinates( i, side1->Shape());

  xi_1_avg.Scale(1.0/ side1->NumNode());

  for(int i=0; i< side2->NumNode(); i++)
    xi_1_avg = DRT::UTILS::getNodeCoordinates( i, side2->Shape());

  xi_2_avg.Scale(1.0/ side2->NumNode());

  xi_side1(0,0) = xi_1_avg(0,0);
  xi_side1(1,0) = xi_1_avg(1,0);

  xi_side2(0,0) = xi_2_avg(0,0);
  xi_side2(1,0) = xi_2_avg(1,0);

  LINALG::Matrix<3,1> proj_x_n_dummy1(true);

  callgetNormalSide_tn(side1, normal1, side1_xyze, lm1, proj_x_n_dummy1, xi_side1);
  callgetNormalSide_tn(side2, normal2, side2_xyze, lm2, proj_x_n_dummy1, xi_side2);

  normal_avg.Update(1.0,normal1, 1.0);
  normal_avg.Update(1.0,normal2, 1.0);

  normal_avg.Scale(1.0/normal_avg.Norm2());

  // map point into fluid domain along normal vector
  start_point.Update(1.0, proj_x_n, dist, normal_avg);

  return;
}

void XFEM::XFLUID_STD::callgetNormalSide_tn(
    DRT::Element*                side,           ///< pointer to side element
    LINALG::Matrix<3,1>&         normal,         ///< normal vector w.r.t side
    Epetra_SerialDenseMatrix &   side_xyze,      ///< side's node coordinates
    const std::vector<int>&      lm,             ///< local map
    LINALG::Matrix<3,1> &        proj_x_n,       ///< projected point on side
    LINALG::Matrix<2,1> &        xi_side         ///< local coordinates of projected point w.r.t side
)
{

  if (side->ElementType() == DRT::ELEMENTS::Bele3Type::Instance()) // three dofs per node, for standard Dirichlet coupling
  {
    const int numdofpernode = 3;

    switch ( side->Shape() )
    {
//      case DRT::Element::tri3:
//      {
//        getNormalSide_tn<DRT::Element::tri3, numdofpernode>(normal, side_xyze, lm, proj_x_n, xi_side);
//        break;
//      }
//      case DRT::Element::tri6:
//      {
//        getNormalSide_tn<DRT::Element::tri6, numdofpernode>(normal, side_xyze, lm, proj_x_n, xi_side);
//        break;
//      }
    case DRT::Element::quad4:
    {
      getNormalSide_tn<DRT::Element::quad4, numdofpernode>(normal, side_xyze, lm, proj_x_n, xi_side);
      break;
    }
    case DRT::Element::quad8:
    {
      getNormalSide_tn<DRT::Element::quad8, numdofpernode>(normal, side_xyze, lm, proj_x_n, xi_side);
      break;
    }
//      case DRT::Element::quad9:
//      {
//        getNormalSide_tn<DRT::Element::quad9, numdofpernode>(normal, side_xyze, lm, proj_x_n, xi_side);
//        break;
//      }
    default:
      dserror( "unsupported side shape %d", side->Shape() );
    }
  }
  else if (side->ElementType() == DRT::ELEMENTS::Bele3_4Type::Instance()) // four dofs per node, for standard Dirichlet coupling
  {
    const int numdofpernode = 4;

    switch ( side->Shape() )
    {
//      case DRT::Element::tri3:
//      {
//        getNormalSide_tn<DRT::Element::tri3, numdofpernode>(normal, side_xyze, lm, proj_x_n, xi_side);
//        break;
//      }
//      case DRT::Element::tri6:
//      {
//        getNormalSide_tn<DRT::Element::tri6, numdofpernode>(normal, side_xyze, lm, proj_x_n, xi_side);
//        break;
//      }
    case DRT::Element::quad4:
    {
      getNormalSide_tn<DRT::Element::quad4, numdofpernode>(normal, side_xyze, lm, proj_x_n, xi_side);
      break;
    }
    case DRT::Element::quad8:
    {
      getNormalSide_tn<DRT::Element::quad8, numdofpernode>(normal, side_xyze, lm, proj_x_n, xi_side);
      break;
    }
//      case DRT::Element::quad9:
//      {
//        getNormalSide_tn<DRT::Element::quad9, numdofpernode>(normal, side_xyze, lm, proj_x_n, xi_side);
//        break;
//      }
    default:
      dserror( "unsupported side shape %d", side->Shape() );
    }
  }
}

template<DRT::Element::DiscretizationType side_distype, const int numdof>
void XFEM::XFLUID_STD::getNormalSide_tn(
    LINALG::Matrix<3,1>&         normal,         ///< normal vector w.r.t side
    Epetra_SerialDenseMatrix &   side_xyze,      ///< side's node coordinates
    const std::vector<int>&      lm,             ///< local map
    LINALG::Matrix<3,1> &        proj_x_n,       ///< projected point on side
    LINALG::Matrix<2,1> &        xi_side         ///< local coordinates of projected point w.r.t side
)
{

  const int side_nen_ = DRT::UTILS::DisTypeToNumNodePerEle<side_distype>::numNodePerElement;

  // add displacements
  addeidisp<side_distype, numdof>(side_xyze, *boundarydis_, "idispn", lm);

  // fill Linalg matrix with current node coordinates including displacements
  LINALG::Matrix<3,side_nen_> xyze_(side_xyze);


  // Initialization
  LINALG::Matrix<side_nen_,1> funct(true);      // shape functions
  LINALG::Matrix<2,side_nen_> deriv(true);      // derivatives dr, ds


  LINALG::Matrix<3,1> x(true);

  LINALG::Matrix<3,2> derxy (true);
  LINALG::Matrix<3,1> dx_dr (true);
  LINALG::Matrix<3,1> dx_ds (true);

  // get current values
  DRT::UTILS::shape_function_2D( funct, xi_side( 0 ), xi_side( 1 ), side_distype );
  DRT::UTILS::shape_function_2D_deriv1( deriv, xi_side( 0 ), xi_side( 1 ), side_distype );

  proj_x_n.Multiply(xyze_, funct);

  derxy.MultiplyNT(xyze_, deriv);


  // set dx_dr and dx_ds
  for (int i=0; i< 3; i++)
  {
    dx_dr(i) = derxy(i,0);
    dx_ds(i) = derxy(i,1);

  }


  normal(0) = dx_dr(1)*dx_ds(2)-dx_ds(1)*dx_dr(2);
  normal(1) = dx_dr(2)*dx_ds(0)-dx_ds(2)*dx_dr(0);
  normal(2) = dx_dr(0)*dx_ds(1)-dx_ds(0)*dx_dr(1);

  normal.Scale(1.0/normal.Norm2());

  return;
}



/*------------------------------------------------------------------------------------------------*
 * call and prepare the projection of point to side                                  schott 07/12 *
 *------------------------------------------------------------------------------------------------*/
void XFEM::XFLUID_STD::call_get_projxn_Line(
    DRT::Element*                    side,                 ///< pointer to structural side element
    DRT::Element*                    line,                 ///< pointer to structural line of side element
    LINALG::Matrix<3,1>&             proj_x_n,
    double &                         xi_line
    )
{

  //-------------------------------------------------

  // line geometry at initial state t^0

  // line geometry
  const int numnodes = line->NumNode();
  DRT::Node ** nodes = line->Nodes();
  Epetra_SerialDenseMatrix line_xyze( 3, numnodes );
  for ( int i=0; i<numnodes; ++i )
  {
    const double * x = nodes[i]->X();
    std::copy( x, x+3, &line_xyze( 0, i ) );
  }


  DRT::Element::LocationArray cutla( 1 );
  line->LocationVector(*boundarydis_,cutla,false);


  if(side->ElementType() == DRT::ELEMENTS::Bele3Type::Instance())
  {
    const int numdofpernode = 3;

    switch ( line->Shape() )
    {
    case DRT::Element::line2:
    {
      get_projxn_Line<DRT::Element::line2, numdofpernode>(line_xyze, cutla[0].lm_, proj_x_n, xi_line);
      break;
    }
    default:
      dserror( "unsupported line shape %d", line->Shape() );
    }
  }
  else if (side->ElementType() == DRT::ELEMENTS::Bele3_4Type::Instance()) // four dofs per node, for standard Dirichlet coupling
  {
    const int numdofpernode = 4;

    switch ( line->Shape() )
    {
    case DRT::Element::line2:
    {
      get_projxn_Line<DRT::Element::line2, numdofpernode>(line_xyze, cutla[0].lm_, proj_x_n, xi_line);
      break;
    }
    default:
      dserror( "unsupported line shape %d", line->Shape() );
    }
  }
}


template<DRT::Element::DiscretizationType line_distype, const int numdof>
void XFEM::XFLUID_STD::get_projxn_Line(
    Epetra_SerialDenseMatrix &   line_xyze,      ///< line's node coordinates
    const std::vector<int>&      lm,             ///< local map
    LINALG::Matrix<3,1> &        proj_x_n,       ///< projected point on side
    double &                     xi_line         ///< local coordinates of projected point w.r.t line
)
{

  const int line_nen_ = DRT::UTILS::DisTypeToNumNodePerEle<line_distype>::numNodePerElement;

   // add displacements
   addeidisp<line_distype, numdof>(line_xyze, *boundarydis_, "idispn", lm);

   // fill LINALG matrix
   LINALG::Matrix<3,line_nen_> xyze_(line_xyze);


  // Initialization
  LINALG::Matrix<line_nen_,1> funct(true);      // shape functions
  LINALG::Matrix<1,line_nen_> deriv(true);      // derivatives dr

  LINALG::Matrix<3,1> x(true);


  // get current values
  DRT::UTILS::shape_function_1D( funct, xi_line, line_distype );

  // projected point tracked back at t^n
  proj_x_n.Multiply(xyze_, funct);


  return;
}


/*--------------------------------------------------------------------------------
 * add side's or line's interface displacements and set current node coordinates
 *--------------------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype, const int numdof>
void XFEM::XFLUID_STD::addeidisp(
    Epetra_SerialDenseMatrix&    xyze,         ///< node coordinates of side or line
    const DRT::Discretization &  cutdis,       ///< cut discretization
    const std::string            state,        ///< state
    const std::vector<int>&      lm            ///< local map
    )
{

  const int nen = DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;

  LINALG::Matrix<3,nen> eidisp(true);


  // get state of the global vector
  Teuchos::RCP<const Epetra_Vector> matrix_state = cutdis.GetState(state);
  if(matrix_state == Teuchos::null)
    dserror("Cannot get state vector %s", state.c_str());

  // extract local values of the global vectors
  std::vector<double> mymatrix(lm.size());
  DRT::UTILS::ExtractMyValues(*matrix_state,mymatrix,lm);

  for (int inode=0; inode<nen; ++inode)  // number of nodes
  {
    for(int idim=0; idim<3; ++idim) // number of dimensions
    {
      (eidisp)(idim,inode) = mymatrix[idim+(inode*numdof)]; // attention! disp state vector has 3+1 dofs for displacement (the same as for (u,p))
    }
  }

  // add the displacement of the interface
  for (int inode = 0; inode < nen; ++inode)
  {
    xyze(0,inode) += eidisp(0, inode);
    xyze(1,inode) += eidisp(1, inode);
    xyze(2,inode) += eidisp(2, inode);
  }

  return;
} // addeidisp


/*------------------------------------------------------------------------------------------------*
 * call and prepare the projection of point to side                                  schott 07/12 *
 *------------------------------------------------------------------------------------------------*/
void XFEM::XFLUID_STD::CallProjectOnSide(
    DRT::Element*          side,                 ///< pointer to structural side element
    LINALG::Matrix<3,1>&   newNodeCoords,        ///< node coordinates of point that has to be projected
    int &                  proj_sid,             ///< id of side when projection lies on side
    double &               min_dist,             ///< minimal distance, potentially updated
    LINALG::Matrix<3,1>&   proj_x_np,            ///< projection of point on this side
    LINALG::Matrix<2,1>&   proj_xi_side,         ///< local coordinates of projection of point w.r.t to this side
    TimeIntData&           data                  ///< reference to data
    )
{
  bool on_side = false;                 ///< lies projection on side?
  double curr_dist = INFINITY;          ///< resulting distance
  LINALG::Matrix<3,1> x_side(true);     ///< resulting projected point
  LINALG::Matrix<2,1> xi_side(true);    ///< local coordinates resulting projected point


  // side geometry at initial state t^0
  const int numnodes = side->NumNode();
  DRT::Node ** nodes = side->Nodes();
  Epetra_SerialDenseMatrix side_xyze( 3, numnodes );
  for ( int i=0; i<numnodes; ++i )
  {
    const double * x = nodes[i]->X();
    std::copy( x, x+3, &side_xyze( 0, i ) );
  }


  DRT::Element::LocationArray cutla( 1 );
  side->LocationVector(*boundarydis_,cutla,false);

  if (side->ElementType() == DRT::ELEMENTS::Bele3Type::Instance()) // three dofs per node, for standard Dirichlet coupling
  {
    const int numdofpernode = 3;

    switch ( side->Shape() )
    {
//      case DRT::Element::tri3:
//      {
//        on_side = ProjectOnSide<DRT::Element::tri3,numdofpernode>(side_xyze, cutla[0].lm_,newNodeCoords,x_side,xi_side);
//        break;
//      }
//      case DRT::Element::tri6:
//      {
//        on_side = ProjectOnSide<DRT::Element::tri6,numdofpernode>(side_xyze, cutla[0].lm_,newNodeCoords,x_side,xi_side);
//        break;
//      }
    case DRT::Element::quad4:
    {
      on_side = ProjectOnSide<DRT::Element::quad4,numdofpernode>(side_xyze, cutla[0].lm_,newNodeCoords,x_side,xi_side, curr_dist);
      break;
    }
    case DRT::Element::quad8:
    {
//        on_side = ProjectOnSide<DRT::Element::quad8,numdofpernode>(side_xyze, cutla[0].lm_,newNodeCoords,x_side,xi_side, curr_dist);
      break;
    }
//      case DRT::Element::quad9:
//      {
//        on_side = ProjectOnSide<DRT::Element::quad9,numdofpernode>(side_xyze, cutla[0].lm_,newNodeCoords,x_side,xi_side, curr_dist);
//        break;
//      }
    default:
      dserror( "unsupported side shape %d", side->Shape() );
    }
  }
  else if (side->ElementType() == DRT::ELEMENTS::Bele3_4Type::Instance()) // four dofs per node, for standard Dirichlet coupling
  {
    const int numdofpernode = 4;

    switch ( side->Shape() )
    {
//      case DRT::Element::tri3:
//      {
//        on_side = ProjectOnSide<DRT::Element::tri3,numdofpernode>(side_xyze, cutla[0].lm_,newNodeCoords,x_side,xi_side, curr_dist);
//        break;
//      }
//      case DRT::Element::tri6:
//      {
//        on_side = ProjectOnSide<DRT::Element::tri6,numdofpernode>(side_xyze, cutla[0].lm_,newNodeCoords,x_side,xi_side, curr_dist);
//        break;
//      }
    case DRT::Element::quad4:
    {
      on_side = ProjectOnSide<DRT::Element::quad4,numdofpernode>(side_xyze, cutla[0].lm_,newNodeCoords,x_side,xi_side, curr_dist);
      break;
    }
    case DRT::Element::quad8:
    {
      //      on_side = ProjectOnSide<DRT::Element::quad8,numdofpernode>(side_xyze, cutla[0].lm_,newNodeCoords,x_side,xi_side, curr_dist);
      break;
    }
//      case DRT::Element::quad9:
//      {
//        //      on_side = ProjectOnSide<DRT::Element::quad9,numdofpernode>(side_xyze, cutla[0].lm_,newNodeCoords,x_side,xi_side, curr_dist);
//        break;
//      }
    default:
      dserror( "unsupported side shape %d", side->Shape() );
    }
  }

  // update minimal distance if possible
  if(on_side)
  {
#ifdef DEBUG_TIMINT_STD
    cout << "\t\t\tprojection of current node lies on side " << side->Id() << " with distance " << curr_dist << endl;
#endif

    if(curr_dist < min_dist && curr_dist > 0)
    {

#ifdef DEBUG_TIMINT_STD
      cout << "\t\t\t\tupdated smallest distance! " << curr_dist << endl;
#endif


      //--------------------
      // set current minimal distance w.r.t side
      data.proj_ = TimeIntData::onSide_;
      // set side that contains the projected point
      proj_sid = side->Id();
      // update minimal distance
      min_dist = curr_dist;
      // update projection point
      proj_x_np.Update(1.0, x_side ,0.0);
      // update local coordinates w.r.t side
      proj_xi_side.Update(1.0, xi_side, 0.0);
      //--------------------

    }
    else if(curr_dist < 0)
    {
#ifdef DEBUG_TIMINT_STD
      cout << "\t\t\t\tnegative distance -> do not update" << endl;
#endif
    }
  }

  return;
}



/*------------------------------------------------------------------------------------------------*
 * call and prepare the projection of point to side                                  schott 07/12 *
 *------------------------------------------------------------------------------------------------*/
void XFEM::XFLUID_STD::CallProjectOnLine(
    DRT::Element*                    side,                 ///< pointer to structural side element
    DRT::Element*                    line,                 ///< pointer to structural line of side element
    int                              line_count,           ///< local line id w.r.t side element
    LINALG::Matrix<3,1>&             newNodeCoords,        ///< node coordinates of point that has to be projected
    double &                         min_dist,             ///< minimal distance, potentially updated
    LINALG::Matrix<3,1>&             proj_x_np,            ///< projection of point on this side
    std::map<vector<int>, std::vector<double> >&  proj_xi_line,         ///< map<sorted nids, local line coordinates of projection of point w.r.t sides >
    std::map<vector<int>, vector<int> >&     proj_lineid,          ///< map<sorted nids, local line id w.r.t sides>
    std::map<vector<int>, vector<int> >&     proj_nid_line,        ///< map<sorted nids, side Ids>
    int &                            proj_sid,             ///< id of side that contains the projected point
    TimeIntData&                     data                  ///< reference to data
    )
{


  bool on_line = false;                ///< is projection on line
  double curr_dist = INFINITY;         ///< resulting distance
  LINALG::Matrix<3,1> x_line(true);    ///< resulting projected point
  double xi_line = INFINITY;           ///< local coordinates resulting projected point

  //-------------------------------------------------

  // line geometry at initial state t^0

  // line geometry
  const int numnodes = line->NumNode();
  DRT::Node ** nodes = line->Nodes();
  Epetra_SerialDenseMatrix line_xyze( 3, numnodes );
  for ( int i=0; i<numnodes; ++i )
  {
    const double * x = nodes[i]->X();
    std::copy( x, x+3, &line_xyze( 0, i ) );
  }


  DRT::Element::LocationArray cutla( 1 );
  line->LocationVector(*boundarydis_,cutla,false);


  if(side->ElementType() == DRT::ELEMENTS::Bele3Type::Instance())
  {
    const int numdofpernode = 3;

    switch ( line->Shape() )
    {
    case DRT::Element::line2:
    {
      on_line = ProjectOnLine<DRT::Element::line2,numdofpernode>(line_xyze, cutla[0].lm_,newNodeCoords,x_line,xi_line, curr_dist);
      break;
    }
    default:
      dserror( "unsupported line shape %d", line->Shape() );
    }
  }
  else if (side->ElementType() == DRT::ELEMENTS::Bele3_4Type::Instance()) // four dofs per node, for standard Dirichlet coupling
  {
    const int numdofpernode = 4;

    switch ( line->Shape() )
    {
    case DRT::Element::line2:
    {
      on_line = ProjectOnLine<DRT::Element::line2,numdofpernode>(line_xyze, cutla[0].lm_,newNodeCoords,x_line,xi_line, curr_dist);
      break;
    }
    default:
      dserror( "unsupported line shape %d", line->Shape() );
    }
  }

  // update minimal distance if possible
  if(on_line)
  {
#ifdef DEBUG_TIMINT_STD
    cout << "\t\t\tprojection lies on line " << " for side "<< side->Id() << endl;
#endif

    const double TOL_dist = 1e-12;

    if(curr_dist < (min_dist+TOL_dist) && curr_dist > 0)
    {
      cout.precision(15);
      cout << "\t\t\t\tupdated smallest distance!" << curr_dist << endl;

      if(curr_dist < (min_dist-TOL_dist)) // smaller distance found
      {
        // another line found, reset already found lines
        proj_lineid.clear();
        proj_xi_line.clear();
        proj_nid_line.clear();
      }
      else
      {
        // add lines, that have the same distance up to TOL_dist
      }

      vector<int> line_nids;
      for ( int i=0; i<numnodes; ++i )
      {
        line_nids.push_back(nodes[i]->Id());
      }

      sort(line_nids.begin(), line_nids.end());

      // reset side id for projection on side
      proj_sid= -1;

      //--------------------
      // set current minimal distance w.r.t line
      data.proj_ = TimeIntData::onLine_;
      // update minimal distance
      min_dist = curr_dist;
      // update projection point
      proj_x_np.Update(1.0, x_line ,0.0);



      std::map<vector<int>, vector<int> >::iterator lines = proj_nid_line.find(line_nids);

      // line already inserted via other side
      if(lines != proj_nid_line.end())
      {

        // set line id w.r.t side->Id() that contains the projected point
        (lines->second).push_back(side->Id());
        // update local line id w.r.t side
        (proj_lineid.find(line_nids))->second.push_back(line_count);
        // update local coordinates w.r.t line
        (proj_xi_line.find(line_nids))->second.push_back(xi_line);
      }
      else
      {

        vector<int> sideids;
        sideids.push_back(side->Id());
        vector<int> locallineids;
        locallineids.push_back(line_count);
        std::vector<double> locallineXiCoords;
        locallineXiCoords.push_back(xi_line);

        // set line id w.r.t side->Id() that contains the projected point
        proj_nid_line.insert(std::pair<vector<int>,vector<int> >(line_nids, sideids));
        // update local line id w.r.t side
        proj_lineid.insert(std::pair<vector<int>,vector<int> >(line_nids,locallineids));
        // update local coordinates w.r.t line
        proj_xi_line.insert(std::pair<vector<int>,vector<double> >(line_nids,locallineXiCoords));
      }

      //--------------------

    }
    else if(curr_dist < 0)
    {
      dserror(" no negative distance to line possible");
    }
  } // on_line

  return;
}

/*------------------------------------------------------------------------------------------------*
 * call and prepare the projection of point to point (distance computation)          schott 07/12 *
 *------------------------------------------------------------------------------------------------*/
void XFEM::XFLUID_STD::CallProjectOnPoint(
    DRT::Node*                     node,                 ///< pointer to node
    LINALG::Matrix<3,1>&           newNodeCoords,        ///< node coordinates of point that has to be projected
    double &                       min_dist,             ///< minimal distance, potentially updated
    int &                          proj_nid_np,          ///< nid id of projected point on surface
    LINALG::Matrix<3,1>&           proj_x_np,            ///< projection of point on this side
    int &                          proj_sid,             ///< id of side that contains the projected point
    std::map<vector<int>, std::vector<double> > proj_xi_line,         ///< map<side ID, local coordinates of projection of point w.r.t to this line>
    std::map<vector<int>, vector<int> >    proj_lineid,          ///< map<side ID, local line id>
    std::map<vector<int>, vector<int> >    proj_nid_line,        ///< map<side ID, vec<line Ids>>
    TimeIntData&                   data                  ///< reference to data
)
{

  double curr_dist = INFINITY;          ///< resulting distance
  LINALG::Matrix<3,1> x_point(true);    ///< resulting projected point


  // its point geometry
  Epetra_SerialDenseMatrix point_xyze( 3, 1 );

  const double * x = node->X();
  std::copy( x, x+3, &point_xyze( 0, 0 ) );

  vector<int> lm = boundarydis_->Dof(0,node);


  // compute distance between two points
  ProjectOnPoint(point_xyze, lm,newNodeCoords, x_point, curr_dist);

  // update minimal distance if possible
  if(curr_dist < min_dist)
  {
#ifdef DEBUG_TIMINT_STD
    cout << "\t\t\t\tupdated smallest distance!" << curr_dist << endl;
#endif
    //--------------------
    // set current minimal distance w.r.t point
    data.proj_ = TimeIntData::onPoint_;
    //reset side id
    proj_sid = -1;
    // reset line id w.r.t side->Id() that contains the projected point
    proj_lineid.clear();
    proj_xi_line.clear();
    proj_nid_line.clear();
    // update minimal distance
    min_dist = curr_dist;
    // update point with minimal distance
    proj_x_np.Update(1.0, x_point, 0.0);
    // update node id of point with minimal distance
    proj_nid_np = node->Id();
    //--------------------

  }

  return;
}


/*--------------------------------------------------------------------------------*
 * project point from in normal direction onto corresponding side    schott 07/12 *
 *--------------------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType side_distype, const int numdof>
bool XFEM::XFLUID_STD::ProjectOnSide(
    Epetra_SerialDenseMatrix &   side_xyze,      ///< side's node coordinates
    const std::vector<int>&      lm,             ///< local map
    LINALG::Matrix<3,1> &        x_gp_lin,       ///< global coordinates of point that has to be projected
    LINALG::Matrix<3,1> &        x_side,         ///< projected point on side
    LINALG::Matrix<2,1> &        xi_side,        ///< local coordinates of projected point w.r.t side
    double &                     dist            ///< distance from point to its projection
)
{
  //  TEUCHOS_FUNC_TIME_MONITOR( "FLD::XFluid::XFluidState::ProjectOnSide" );

  const int side_nen_ = DRT::UTILS::DisTypeToNumNodePerEle<side_distype>::numNodePerElement;

  // add displacements
  addeidisp<side_distype, numdof>(side_xyze, *boundarydis_, "idispnp", lm);

  // fill Linalg matrix with current node coordinates including displacements
  LINALG::Matrix<3,side_nen_> xyze_(side_xyze);


  // Initialization
  LINALG::Matrix<side_nen_,1> funct(true);      // shape functions
  LINALG::Matrix<2,side_nen_> deriv(true);      // derivatives dr, ds
  LINALG::Matrix<3,side_nen_> deriv2(true);     // 2nd derivatives drdr, dsds, drds


  LINALG::Matrix<3,1> x(true);

  LINALG::Matrix<3,2> derxy (true);
  LINALG::Matrix<3,1> dx_dr (true);
  LINALG::Matrix<3,1> dx_ds (true);

  LINALG::Matrix<3,3> derxy2 (true);
  LINALG::Matrix<3,1> dx_drdr (true);
  LINALG::Matrix<3,1> dx_dsds (true);
  LINALG::Matrix<3,1> dx_drds (true);

  LINALG::Matrix<3,1> dx_drdr_times_dx_ds(true);
  LINALG::Matrix<3,1> dx_dr_times_dx_drds(true);
  LINALG::Matrix<3,1> dx_drds_times_dx_ds(true);
  LINALG::Matrix<3,1> dx_dr_times_dx_dsds(true);
  LINALG::Matrix<3,1> dx_dr_times_dx_ds(true);

  LINALG::Matrix<3,1> residuum(true);             // residuum of the newton iteration
  LINALG::Matrix<3,3> sysmat(true);               // matrix for the newton system
  LINALG::Matrix<3,1> incr(true);                 // increment of the newton system

  LINALG::Matrix<3,1> sol(true); // sol carries xi_1, xi_2, d (distance)

  if(side_distype == DRT::Element::tri3 or
     side_distype == DRT::Element::tri6)
  {
    sol(0) = 0.333333333333333;
    sol(1) = 0.333333333333333;
  }
  else if( side_distype == DRT::Element::quad4 or
           side_distype == DRT::Element::quad8 or
           side_distype == DRT::Element::quad9)
  {

    sol(0) = 0.0;
    sol(1) = 0.0;
  }
  else
  {
    dserror("define start side xi-coordinates for unsupported cell type");
  }

  const double absTolIncr = 1.0e-9;   // rel tolerance for the local coordinates increment
  const double relTolRes  = 1.0e-9;   // rel tolerance for the whole residual
  const double absTOLdist = 1.0e-9;   // abs tolerance for distance

  int iter=0;
  const int maxiter = 7;

  bool converged = false;

  while(iter < maxiter && !converged)
  {

    iter++;


    // get current values
    DRT::UTILS::shape_function_2D( funct, sol( 0 ), sol( 1 ), side_distype );
    DRT::UTILS::shape_function_2D_deriv1( deriv, sol( 0 ), sol( 1 ), side_distype );
    DRT::UTILS::shape_function_2D_deriv2( deriv2, sol( 0 ), sol( 1 ), side_distype );

    x.Multiply(xyze_, funct);

    derxy.MultiplyNT(xyze_, deriv);

    derxy2.MultiplyNT(xyze_, deriv2);

    // set dx_dr and dx_ds
    for (int i=0; i< 3; i++)
    {
      dx_dr(i) = derxy(i,0);
      dx_ds(i) = derxy(i,1);

      dx_drdr(i) = derxy2(i,0);
      dx_dsds(i) = derxy2(i,1);
      dx_drds(i) = derxy2(i,2);
    }

    // get vector products

    dx_drdr_times_dx_ds(0) = dx_drdr(1)*dx_ds(2)-dx_ds(1)*dx_drdr(2);
    dx_drdr_times_dx_ds(1) = dx_drdr(2)*dx_ds(0)-dx_ds(2)*dx_drdr(0);
    dx_drdr_times_dx_ds(2) = dx_drdr(0)*dx_ds(1)-dx_ds(0)*dx_drdr(1);

    dx_dr_times_dx_drds(0) = dx_dr(1)*dx_drds(2)-dx_drds(1)*dx_dr(2);
    dx_dr_times_dx_drds(1) = dx_dr(2)*dx_drds(0)-dx_drds(2)*dx_dr(0);
    dx_dr_times_dx_drds(2) = dx_dr(0)*dx_drds(1)-dx_drds(0)*dx_dr(1);

    dx_drds_times_dx_ds(0) = dx_drds(1)*dx_ds(2)-dx_ds(1)*dx_drds(2);
    dx_drds_times_dx_ds(1) = dx_drds(2)*dx_ds(0)-dx_ds(2)*dx_drds(0);
    dx_drds_times_dx_ds(2) = dx_drds(0)*dx_ds(1)-dx_ds(0)*dx_drds(1);

    dx_dr_times_dx_dsds(0) = dx_dr(1)*dx_dsds(2)-dx_dsds(1)*dx_dr(2);
    dx_dr_times_dx_dsds(1) = dx_dr(2)*dx_dsds(0)-dx_dsds(2)*dx_dr(0);
    dx_dr_times_dx_dsds(2) = dx_dr(0)*dx_dsds(1)-dx_dsds(0)*dx_dr(1);

    dx_dr_times_dx_ds(0) = dx_dr(1)*dx_ds(2)-dx_ds(1)*dx_dr(2);
    dx_dr_times_dx_ds(1) = dx_dr(2)*dx_ds(0)-dx_ds(2)*dx_dr(0);
    dx_dr_times_dx_ds(2) = dx_dr(0)*dx_ds(1)-dx_ds(0)*dx_dr(1);

    // define sysmat
    for(int i=0; i< 3; i++)
    {
      // d/dr
      sysmat(i,0) = dx_dr(i) - sol(2) * (dx_drdr_times_dx_ds(i) + dx_dr_times_dx_drds(i));

      // d/ds
      sysmat(i,1) = dx_ds(i) - sol(2) * (dx_drds_times_dx_ds(i) + dx_dr_times_dx_dsds(i));

      // d/d(dist)
      sysmat(i,2) = - dx_dr_times_dx_ds(i);


      // residual
      residuum(i) = x(i) - sol(2) * dx_dr_times_dx_ds(i) - x_gp_lin(i);

    }



    sysmat.Invert();

    //solve Newton iteration
    incr.Clear();
    incr.Multiply(-1.0,sysmat,residuum); // incr = -Systemmatrix^-1 * residuum

    // update solution
    sol.Update(1.0, incr, 1.0);

    if ( (incr.Norm2()/sol.Norm2() < absTolIncr) && (residuum.Norm2()/sol.Norm2() < relTolRes) )
    {
      converged = true;
    }

    // check ° relative criterion for local coordinates (between [-1,1]^2)
    //       ° absolute criterion for distance (-> 0)
    //       ° relative criterion for whole residuum
    if(    //sqrt(incr(0)*incr(0)+incr(1)*incr(1))/sqrt(sol(0)*sol(0)+sol(1)*sol(1)) <  relTolIncr
        sqrt(incr(0)*incr(0)+incr(1)*incr(1)) <  absTolIncr
        && incr(2) < absTOLdist
        && residuum.Norm2()/sol.Norm2() < relTolRes)
    {
      converged = true;
    }

  }

  bool on_side = false;

  if(!converged)
  {

#ifdef DEBUG_TIMINT_STD
    cout.precision(15);

    cout << "increment criterion loc coord "
        //<< sqrt(incr(0)*incr(0)+incr(1)*incr(1))/sqrt(sol(0)*sol(0)+sol(1)*sol(1))
        << sqrt(incr(0)*incr(0)+incr(1)*incr(1))
        << " \tabsTOL: " << absTolIncr
        << endl;
    cout << "absolute criterion for distance "
        << incr(2)
        << " \tabsTOL: " << absTOLdist
        << endl;
    cout << "relative criterion whole residuum "
        << residuum.Norm2()/sol.Norm2()
        << " \trelTOL: " << relTolRes
        << endl;


    cout << "sysmat.Invert" << sysmat << endl;
    cout << "sol-norm " << sol.Norm2() << endl;
    cout << "sol " << sol << endl;
    cout << "x_gp_lin" << x_gp_lin << endl;
    cout << "side " << xyze_ << endl;
#endif

    //dserror( "newton scheme in ProjectOnSide not converged! " );

    xi_side(0) = INFINITY;
    xi_side(1) = INFINITY;

    dist = INFINITY;

    on_side = false;
  }
  else
  {
    LINALG::Matrix<3,1> xsi(true);
    xsi(0) = sol(0);
    xsi(1) = sol(1);
    xsi(2) = 0;

    bool within_limits = WithinLimits<side_distype>(xsi);

    if(within_limits)
    {

      on_side = true;


      // get length of current normal vector
      double normal_length = dx_dr_times_dx_ds.Norm2();


      dist = -sol(2)*normal_length; // negative sol(2)!!! and scaling with normal length

      // evaluate shape function at solution
      DRT::UTILS::shape_function_2D( funct, sol( 0 ), sol( 1 ), side_distype );

      // get projected gauss point
      x_side.Multiply(xyze_, funct);

      // set local coordinates w.r.t side
      xi_side(0) = sol(0);
      xi_side(1) = sol(1);

    }
    else
    {
      on_side = false;

      xi_side(0) = INFINITY;
      xi_side(1) = INFINITY;

      dist = INFINITY;

      x_side(0) = INFINITY;
      x_side(1) = INFINITY;
      x_side(2) = INFINITY;
    }

  }

  return on_side;
} // ProjectOnSide

/*--------------------------------------------------------------------------------*
 * project point in normal direction onto corresponding line         schott 07/12 *
 *--------------------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType line_distype, const int numdof>
bool XFEM::XFLUID_STD::ProjectOnLine(
    Epetra_SerialDenseMatrix &   line_xyze,      ///< line's node coordinates
    const std::vector<int>&      lm,             ///< local map
    LINALG::Matrix<3,1> &        x_point_np,     ///< global coordinates of point that has to be projected
    LINALG::Matrix<3,1> &        x_line,         ///< projected point on line
    double &                     xi_line,        ///< local coordinates of projected point w.r.t line
    double &                     dist            ///< distance from point to its projection
)
{

  bool on_line = false;

  const int line_nen_ = DRT::UTILS::DisTypeToNumNodePerEle<line_distype>::numNodePerElement;

  // add displacements
  addeidisp<line_distype, numdof>(line_xyze, *boundarydis_, "idispnp", lm);

  // fill LINALG matrix
  LINALG::Matrix<3,line_nen_> xyze_(line_xyze);

  // linear line
  if(line_nen_ == 2)
  {
    // get line direction
    LINALG::Matrix<3,1> line_dir(true); // x2-x1
    LINALG::Matrix<3,1> dir_tmp(true);  // x - 0.5(x1 + x2)
    for(int isd=0; isd<3; isd++)
    {
      line_dir(isd) = xyze_(isd,1) - xyze_(isd,0);
      dir_tmp(isd) = x_point_np(isd,0) - 0.5*(xyze_(isd,1) + xyze_(isd,0));
    }

    double line_length = line_dir.Norm2(); // || (x2-x1) ||

    if(line_length < 1e-12) dserror("line has lenght smaller than 1e-12");

    //2.0/|| (x2-x1) ||^2 * < x-0.5(x1+x2),x2-x1 >
    xi_line = 2.0/(line_length*line_length) * line_dir.Dot(dir_tmp);

    // check if projection within line segement between x1 and x2
    LINALG::Matrix<3,1> xsi(true);
    xsi(0) = xi_line;
    xsi(1) = 0;
    xsi(2) = 0;

    bool within_limits = WithinLimits<line_distype>(xsi);

    if(within_limits)
    {
      on_line = true;

      // set projection
      for(int i=0; i<3; i++)
      {
        x_line(i) = 0.5*(1+xi_line)*xyze_(i,1) + 0.5*(1-xi_line)*xyze_(i,0);
      }

      // set distance vector between x and projection
      LINALG::Matrix<3,1> dist_vec(true);
      for(int i=0; i<3; i++)
      {
        dist_vec(i) = x_point_np(i) - x_line(i);
      }

      // compute distance
      dist = dist_vec.Norm2();
    }
    else
    {
      on_line = false;
      xi_line = INFINITY;
      dist = INFINITY;

      x_line(0) = INFINITY;
      x_line(1) = INFINITY;
      x_line(2) = INFINITY;
    }
  }
  else
  {
    dserror("projection on line just for linear lines implemented!");
  }

  return on_line;
} // ProjectOnLine

/*--------------------------------------------------------------------------------*
 * compute distance (project) between two points                     schott 07/12 *
 *--------------------------------------------------------------------------------*/
void XFEM::XFLUID_STD::ProjectOnPoint(
    Epetra_SerialDenseMatrix &   point_xyze,     ///< point's node coordinates
    const std::vector<int>&      lm,             ///< local map
    LINALG::Matrix<3,1> &        x_point_np,     ///< global coordinates of point that has to be projected
    LINALG::Matrix<3,1> &        x_point,        ///< global coordinates of point on surface
    double &                     dist            ///< distance from point to its projection
)
{
  const std::string state = "idispnp";

  Teuchos::RCP<const Epetra_Vector> matrix_state = boundarydis_->GetState( state);
  if(matrix_state == Teuchos::null)
    dserror("Cannot get state vector %s", state.c_str());

  // extract local values of the global vectors
  std::vector<double> mymatrix(lm.size());
  DRT::UTILS::ExtractMyValues(*matrix_state,mymatrix,lm);

  // add the displacement of the interface
  for(int idim=0; idim<3; ++idim) // number of dimensions
  {
    (point_xyze)(idim,0) += mymatrix[idim]; // attention! disp state vector has 3+1 dofs for displacement (the same as for (u,p))
  }

  // fill LINALG matrix
  LINALG::Matrix<3,1> p(point_xyze);

  // compute direction vector between two points
  LINALG::Matrix<3,1> direction(true);
  direction.Update(1.0, x_point_np, -1.0, p);

  // compute distance
  dist = direction.Norm2();

  x_point.Update(1.0, p ,0.0);

  return;
}// ProjectOnPoint


/*--------------------------------------------------------------------------------*
 * check if local coordinates are within limits                      schott 07/12 *
 *--------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType elementtype>
bool XFEM::XFLUID_STD::WithinLimits(LINALG::Matrix<3,1>& xsi_)
{

  switch ( elementtype )
  {
  case DRT::Element::line2:
  case DRT::Element::line3:
    return ( xsi_( 0 ) >= -1-TOLERANCE and xsi_( 0 ) <=  1+TOLERANCE);
  case DRT::Element::tri3:
  case DRT::Element::tri6:
    return ( xsi_( 0 ) >=   -TOLERANCE and xsi_( 1 ) >=   -TOLERANCE and
             xsi_( 0 ) <=  1+TOLERANCE and xsi_( 1 ) <=  1+TOLERANCE and
            (xsi_( 0 ) + xsi_( 1 )) <= 1+TOLERANCE);
  case DRT::Element::quad4:
  case DRT::Element::quad8:
  case DRT::Element::quad9:
    return ( xsi_( 0 ) >= -1-TOLERANCE and xsi_( 1 ) >= -1-TOLERANCE and
             xsi_( 0 ) <=  1+TOLERANCE and xsi_( 1 ) <=  1+TOLERANCE );
  case DRT::Element::hex8:
  case DRT::Element::hex20:
  case DRT::Element::hex27:
    return ( xsi_( 0 ) >= -1-TOLERANCE and xsi_( 1 ) >= -1-TOLERANCE and xsi_( 2 ) >= -1-TOLERANCE and
             xsi_( 0 ) <=  1+TOLERANCE and xsi_( 1 ) <=  1+TOLERANCE and xsi_( 2 ) <=  1+TOLERANCE );
  case DRT::Element::tet4:
  case DRT::Element::tet10:
    return ( xsi_( 0 ) >=   -TOLERANCE and xsi_( 1 ) >=   -TOLERANCE and xsi_( 2 ) >=   -TOLERANCE and
             xsi_( 0 )+xsi_( 1 )+xsi_( 2 ) <= 1+TOLERANCE );
  case DRT::Element::pyramid5:
    return ( xsi_( 0 ) >= -1-TOLERANCE and xsi_( 1 ) >= -1-TOLERANCE and xsi_( 2 ) >=   -TOLERANCE and
             xsi_( 0 ) <=  1+TOLERANCE and xsi_( 1 ) <=  1+TOLERANCE and
             ( fabs( xsi_( 0 ) ) > fabs( xsi_( 1 ) ) ?
               ( xsi_( 2 )+fabs( xsi_( 0 ) ) <= 1+TOLERANCE ) :
               ( xsi_( 2 )+fabs( xsi_( 1 ) ) <= 1+TOLERANCE ) )
      );
  case DRT::Element::wedge6:
  case DRT::Element::wedge15:
    return ( xsi_( 0 ) >=   -TOLERANCE and xsi_( 1 ) >=   -TOLERANCE and xsi_( 2 ) >= -1-TOLERANCE and
             xsi_( 2 ) <=  1+TOLERANCE and
             xsi_( 0 )+xsi_( 1 ) <= 1+TOLERANCE );
  default:
    throw std::runtime_error( "unsupported element type in XFEM::XFLUID_STD::WithinLimits" );
  }
  return false;
}

/*------------------------------------------------------------------------------------------------*
 * setting the computed data for the standard degrees of freedom into the according               *
 * Epetra Vectors for all handled nodes                                              schott 06/10 *
 *------------------------------------------------------------------------------------------------*/
void XFEM::XFLUID_STD::setFinalData(
)
{
  const int nsd = 3; // 3 dimensions for a 3d fluid element

  // loop over data
  for (vector<TimeIntData>::iterator data=timeIntData_->begin();
      data!=timeIntData_->end(); data++)
  {
    if (data->state_!=TimeIntData::doneStd_)
      dserror("when data is set, all computation has to be done");

    vector<LINALG::Matrix<nsd,1> >& velValues(data->velValues_); // velocities of the node
    vector<double>& presValues(data->presValues_); // pressures of the node

    const int gnodeid = data->node_.Id(); // global node id

    //-------------------------------------------------------
    DRT::Node* node = discret_->gNode(gnodeid);

    std::vector<int> dofs;
    dofset_new_->Dof(*node, data->nds_np_, dofs );

    if(dofs.size() != (nsd+1) ) dserror("not the right number of dofs %d for this node ", dofs.size());


    // set velocity dofs
    for(int i=0; i< nsd; i++)
    {
      for (size_t index=0;index<vectorSize(data->type_);index++)
        (*newVectors_[index])[newdofrowmap_.LID(dofs[i])] = velValues[index](i,0); // set the value
    }
    for (size_t index=0;index<vectorSize(data->type_);index++)
      (*newVectors_[index])[newdofrowmap_.LID(dofs[3])] = presValues[index]; // set the value

    data->type_ = TimeIntData::standard_; // predictor is done, so next time standard
  } // end loop over nodes
} // end setFinalData



#ifdef PARALLEL
/*------------------------------------------------------------------------------------------------*
 * export start data to neighbour proc                                           winklmaier 06/10 *
 *------------------------------------------------------------------------------------------------*/
void XFEM::XFLUID_STD::exportStartData()
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
  for (vector<TimeIntData>::iterator data=timeIntData_->begin();
      data!=timeIntData_->end(); data++)
  {
    packNode(dataSend,data->node_);
    DRT::ParObject::AddtoPack(dataSend,data->nds_np_);
    DRT::ParObject::AddtoPack(dataSend,data->vel_);
    DRT::ParObject::AddtoPack(dataSend,data->velDeriv_);
    DRT::ParObject::AddtoPack(dataSend,data->presDeriv_);
    DRT::ParObject::AddtoPack(dataSend,data->startpoint_);
//    DRT::ParObject::AddtoPack(dataSend,data->phiValue_);
    DRT::ParObject::AddtoPack(dataSend,data->searchedProcs_);
    DRT::ParObject::AddtoPack(dataSend,data->counter_);
    DRT::ParObject::AddtoPack(dataSend,data->startGid_);
    DRT::ParObject::AddtoPack(dataSend,data->startOwner_);
    DRT::ParObject::AddtoPack(dataSend,data->dMin_);
    DRT::ParObject::AddtoPack(dataSend,(int)data->type_);
  }

  dataSend.StartPacking();

  for (vector<TimeIntData>::iterator data=timeIntData_->begin();
      data!=timeIntData_->end(); data++)
  {
    packNode(dataSend,data->node_);
    DRT::ParObject::AddtoPack(dataSend,data->nds_np_);
    DRT::ParObject::AddtoPack(dataSend,data->vel_);
    DRT::ParObject::AddtoPack(dataSend,data->velDeriv_);
    DRT::ParObject::AddtoPack(dataSend,data->presDeriv_);
    DRT::ParObject::AddtoPack(dataSend,data->startpoint_);
//    DRT::ParObject::AddtoPack(dataSend,data->phiValue_);
    DRT::ParObject::AddtoPack(dataSend,data->searchedProcs_);
    DRT::ParObject::AddtoPack(dataSend,data->counter_);
    DRT::ParObject::AddtoPack(dataSend,data->startGid_);
    DRT::ParObject::AddtoPack(dataSend,data->startOwner_);
    DRT::ParObject::AddtoPack(dataSend,data->dMin_);
    DRT::ParObject::AddtoPack(dataSend,(int)data->type_);
  }

  vector<char> dataRecv;
  sendData(dataSend,dest,source,dataRecv);

  // pointer to current position of group of cells in global string (counts bytes)
  vector<char>::size_type posinData = 0;

  // clear vector that should be filled
  timeIntData_->clear();

  // unpack received data
  while (posinData < dataRecv.size())
  {
    double coords[nsd] = {0.0};
    DRT::Node node(0,(double*)coords,0); // initialize node
    int nds_np = -1;
    LINALG::Matrix<nsd,1> vel; // velocity at point x
    vector<LINALG::Matrix<nsd,nsd> > velDeriv; // derivation of velocity at point x
    vector<LINALG::Matrix<1,nsd> > presDeriv; // derivation of pressure at point x
    LINALG::Matrix<nsd,1> startpoint; // startpoint
//    double phiValue; // phi-value
    int searchedProcs; // number of searched processors
    int counter; // iteration counter
    vector<int> startGid; // global id of first startpoint
    vector<int> startOwner; // owner of first startpoint
    double dMin; // minimal distance
    int newtype; // type of the data

    unpackNode(posinData,dataRecv,node);
    DRT::ParObject::ExtractfromPack(posinData,dataRecv,nds_np);
    DRT::ParObject::ExtractfromPack(posinData,dataRecv,vel);
    DRT::ParObject::ExtractfromPack(posinData,dataRecv,velDeriv);
    DRT::ParObject::ExtractfromPack(posinData,dataRecv,presDeriv);
    DRT::ParObject::ExtractfromPack(posinData,dataRecv,startpoint);
//    DRT::ParObject::ExtractfromPack(posinData,dataRecv,phiValue);
    DRT::ParObject::ExtractfromPack(posinData,dataRecv,searchedProcs);
    DRT::ParObject::ExtractfromPack(posinData,dataRecv,counter);
    DRT::ParObject::ExtractfromPack(posinData,dataRecv,startGid);
    DRT::ParObject::ExtractfromPack(posinData,dataRecv,startOwner);
    DRT::ParObject::ExtractfromPack(posinData,dataRecv,dMin);
    DRT::ParObject::ExtractfromPack(posinData,dataRecv,newtype);

    timeIntData_->push_back(TimeIntData(
        node,
        nds_np,
        vel,
        velDeriv,
        presDeriv,
        startpoint,
//        phiValue,
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
void XFEM::XFLUID_STD::exportFinalData()
{
  const int nsd = 3; // 3 dimensions for a 3d fluid element

  // array of vectors which stores data for
  // every processor in one vector
  vector<std::vector<TimeIntData> > dataVec(numproc_);

  // fill vectors with the data
  for (vector<TimeIntData>::iterator data=timeIntData_->begin();
      data!=timeIntData_->end(); data++)
  {
    if (data->state_!=TimeIntData::doneStd_)
      dserror("All data should be set here, having status 'done'. Thus something is wrong!");
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
    if (source<0)				source+=numproc_;
    else if (source>=numproc_)	source-=numproc_;

    DRT::PackBuffer dataSend;

    // pack data to be sent
    for (vector<TimeIntData>::iterator data=dataVec[dest].begin();
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

    for (vector<TimeIntData>::iterator data=dataVec[dest].begin();
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

    vector<char> dataRecv;
    sendData(dataSend,dest,source,dataRecv);

    // pointer to current position of group of cells in global string (counts bytes)
    vector<char>::size_type posinData = 0;

    // unpack received data
    while (posinData < dataRecv.size())
    {
      int gid; // global id of node
      LINALG::Matrix<nsd,1> startpoint; // startpoint
      double phiValue; // phi-value
      vector<int> startGid; // global id of first startpoint
      vector<int> startOwner; // owner of first startpoint
      vector<LINALG::Matrix<nsd,1> > velValues; // velocity values
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



