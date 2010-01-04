/*!-----------------------------------------------------------------------------------------------*
 \file combust_reinitializer.cpp

  \brief reinitialization of G-function (level set function)

  detailed description in header file combust_interface.H

<pre>
Maintainer: Florian Henke
            henke@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15265
</pre>
 *------------------------------------------------------------------------------------------------*/
#ifdef CCADISCRET

#include "combust_reinitializer.H"
#include "combust_defines.H"
#include "../drt_fem_general/drt_utils_shapefunctions_service.H"
#include "../drt_lib/linalg_fixedsizematrix.H"

/*------------------------------------------------------------------------------------------------*
 | constructor                                                                        henke 12/09 |
 *------------------------------------------------------------------------------------------------*/
COMBUST::Reinitializer::Reinitializer(
    const Teuchos::ParameterList& combustdyn,
    SCATRA::ScaTraTimIntImpl& scatra,
    const std::map<int,GEO::BoundaryIntCells>& boundaryintcells
    ) :
    combustdyn_(combustdyn),
    scatra_(scatra),
    flamefront_(boundaryintcells),
    reinitaction_(Teuchos::getIntegralValue<INPAR::COMBUST::ReInitialActionGfunc>(combustdyn_.sublist("COMBUSTION GFUNCTION"),"REINITIALIZATION"))
{
  switch (reinitaction_)
  {
    case INPAR::COMBUST::reinitaction_none:
      // do nothing
      break;
    case INPAR::COMBUST::reinitaction_byfunction:
      // read a FUNCTION from the input file and reinitialize the G-function by evaluating it
      scatra_.SetInitialField(INPAR::SCATRA::initfield_field_by_function,combustdyn_.sublist("COMBUSTION GFUNCTION").get<int>("REINITFUNCNO"));
      break;
    case INPAR::COMBUST::reinitaction_signeddistancefunction:
#ifdef PARALLEL
      //dserror("direct computation of signed distance function not available in parallel");
#endif
      SignedDistanceFunction();
      break;
//    case INPAR::COMBUST::sussman:
//      // SCATRA parameters have to be set before in ScaTraFluidCouplingAlgorithm!
//      ScaTraField().Reinitialize();
//      break;
    default:
      dserror ("unknown option to reinitialize the G-function");
  }

  return;
}


/*------------------------------------------------------------------------------------------------*
 | destructor                                                                         henke 12/09 |
 *------------------------------------------------------------------------------------------------*/
COMBUST::Reinitializer::~Reinitializer()
{
  return;
}


/*------------------------------------------------------------------------------------------------*
 | private: build signed distance function                                            henke 06/08 |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::Reinitializer::SignedDistanceFunction()
{
  // get communicator (for output)
  const Epetra_Comm& comm = scatra_.Discretization()->Comm();
  if (comm.MyPID()==0)
    std::cout << "\n--- reinitializing G-function with signed distance function ..." << std::flush;

  // get a pointer to the G-function discretization
  Teuchos::RCP<DRT::Discretization> gfuncdis = scatra_.Discretization();

  const Epetra_Map* dofrowmap = gfuncdis->DofRowMap();

  // vector of coordinates of node, the G-function value of which is to be reinitialized
  static LINALG::Matrix<3,1> nodecoord(true);
  // normal vector of a flame front patch
  static LINALG::Matrix<3,1> normal(true);

  // loop all row nodes on the processor
  for(int lnodeid=0; lnodeid < gfuncdis->NumMyRowNodes(); ++lnodeid)
  {
    // get the processor local node
    const DRT::Node* lnode = gfuncdis->lRowNode(lnodeid);
//cout << "proc " << comm.MyPID() << " reinitialization for node: " << lnode->Id() << endl;

    // the set of degrees of freedom associated with the node
    const vector<int> nodedofset = gfuncdis->Dof(lnode); // this should not be a vector, but a scalar

#ifdef DEBUG
    int numdofs = nodedofset.size(); // this should be 1 (scalar field!)
    if (numdofs != 1)
      dserror("There are more than 1 dof for a node in a scalar field!");
#endif

    const int dofgid = nodedofset[0];
    int doflid = dofrowmap->LID(dofgid);

    if (fabs((*scatra_.Phinp())[doflid]) < 0.15);
    {

    //-----------------------------------------------------------
    // compute smallest distance to the flame front for this node
    //-----------------------------------------------------------
    // smallest distance to the vertex of a flame front patch
    double vertexdist = 7777.7; // default value
    // smallest distance to flame front
    double mindist = 5555.5; // default value

    // get physical coordinates of this node
    nodecoord(0) = lnode->X()[0];
    nodecoord(1) = lnode->X()[1];
    nodecoord(2) = lnode->X()[2];

    // loop groups (vectors) of flamefront patches of all elements
    for(std::map<int,GEO::BoundaryIntCells>::const_iterator elepatches = flamefront_.begin(); elepatches != flamefront_.end(); ++elepatches)
    {
      // number of flamefront patches for this element
      const std::vector<GEO::BoundaryIntCell> patches = elepatches->second;
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
            //cout << "distance to flame front patch: " << mindist << " is overwritten by: " << patchdist << endl;
            mindist = patchdist;
          }
        }

        //----------------------------------------------------------------
        // compute smallest distance to vertices of this flame front patch
        //----------------------------------------------------------------
        ComputeDistanceToPatch(nodecoord,patch,patchcoord,normal,vertexdist);
      }
    }

    //-------------------------------------------
    // determine smallest distance to flame front
    //-------------------------------------------
    // remark: variables have the following meaning here:
    //         - "mindist" is either "patchdist" or still the default value (5555.5)
    //         - "vertexdist" is the distance to the clostest vertex of any flame front patch
    //
    // case 1: a local flame front patch was found for this node -> mindist = patchdist
    // case 2: a flame front patch was found for this node, but it is not local (curved interface);
    //         a local patch was not found, because this node is located in the blind angle of all
    //         local patches -> mindist = vertexdist
    // case 3: something went wrong:  "mindist" still has default value (5555.5)

    //cout << "mindist " << std::setw(18) << std::setprecision(12) << std::scientific << mindist << endl;
    //cout << "vertexdist " << std::setw(18) << std::setprecision(12) << std::scientific << vertexdist << endl;
    if (fabs(vertexdist) < fabs(mindist)) // case 2
    {
      mindist = vertexdist;
      //cout << "node " << lnode->Id() << " does not face a patch (blind angle) or this patch is not local; distance: " << mindist << endl;
    }
    if (mindist == 5555.5) // case 3
      dserror ("G-fuction value of node %d was not reinitialized!", lnode->Id());

    //------------------------------
    // reinitialize G-function field
    //------------------------------
    // assign new values of signed distance function to the G-function field
    scatra_.Phinp()->ReplaceMyValues(1,&mindist,&doflid);
  }
  }
  if (comm.MyPID()==0)
    std::cout << " done" << std::endl;;

  return;
}


/*------------------------------------------------------------------------------------------------*
 | private: find a facing flame front patch by projecton of node into boundary cell space         |
 |                                                                                    henke 12/09 |
 *----------------------------------------------------------------------------------------------- */
void COMBUST::Reinitializer::FindFacingPatchProjCellSpace(
    const LINALG::Matrix<3,1>&       node,
    const GEO::BoundaryIntCell&      patch,
    const LINALG::SerialDenseMatrix& patchcoord,
    const LINALG::Matrix<3,1>&       normal,
    bool&                            facenode,
    double&                          patchdist)
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
  }

  // Newton iteration converged
  //  cout << "Newton iteration converged in " << iter << " steps!" << endl;

  //----------------------------------------------------
  // check if projection lies within boundary cell space
  //----------------------------------------------------
  switch(patch.Shape())
  {
  case DRT::Element::tri3:
  {
    // criteria for tri3 patch
    if ((eta(0) > -1.0E-12) and (eta(0) < 1.0+1.0E-12) and
        (eta(1) > -1.0E-12) and (eta(1) < 1.0+1.0E-12) and
        (1.0-eta(0)-eta(1) > -1.0E-12) and (1.0-eta(0)-eta(1) < 1.0+1.0E-12) and
        converged)
    {
      facenode = true;
      patchdist = alpha;
      //cout << "facing patch found (tri3 patch)! distance: " << alpha << endl;
    }
    break;
  }
  case DRT::Element::quad4:
  {
    // criteria for quad4 patch
    if ((eta(0) > -1.0-1.0E-12) and (eta(0) < 1.0+1.0E-12) and
        (eta(1) > -1.0-1.0E-12) and (eta(1) < 1.0+1.0E-12) and
        converged)
    {
      facenode = true;
      patchdist = alpha;
      cout << "facing patch found (quad4 patch)!" << endl;
    }
    break;
  }
  default:
    dserror("unknown type of boundary integration cell");
  }
//  if (!converged)
//  {
//    cout << "node x component " << node(0,0) << endl;
//    cout << "node y component " << node(1,0) << endl;
//    cout << "node z component " << node(2,0) << endl;
//    cout << "eta1 " << eta(0) << endl;
//    cout << "eta2 " << eta(1) << endl;
//    cout << "alpha " << alpha << endl;
//    cout << "patch vertices x component " << patchcoord(0,0) << " " << patchcoord(0,1) << " " << patchcoord(0,2) << endl;
//    cout << "patch vertices y component " << patchcoord(1,0) << " " << patchcoord(1,1) << " " << patchcoord(1,2) << endl;
//    cout << "patch vertices z component " << patchcoord(2,0) << " " << patchcoord(2,1) << " " << patchcoord(2,2) << endl;
//  }

  return;
}


/*------------------------------------------------------------------------------------------------*
 | private: compute distance to vertex of patch                                       henke 08/09 |
 *----------------------------------------------------------------------------------------------- */
void COMBUST::Reinitializer::ComputeDistanceToPatch(
    const LINALG::Matrix<3,1>&       node,
    const GEO::BoundaryIntCell&      patch,
    const LINALG::SerialDenseMatrix& patchcoord,
    const LINALG::Matrix<3,1>&       normal,
    double&                          vertexdist
)
{
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
    double normdist = sqrt(dist(0)*dist(0) + dist(1)*dist(1) + dist(2)*dist(2));
    if(normdist < fabs(vertexdist))
    {
      // determine sign of distance vector from node to flame front
      double tmp = normal(0)*dist(0) + normal(1)*dist(1) + normal(2)*dist(2);
      // tmp < 0 if node in unburnt domain (G<0) and vice versa
      if (tmp <= 0.0)
        vertexdist = normdist;
      else
        vertexdist = -normdist;
      //cout << "distance to vertex is overwritten by: " << vertexdist << endl;
    }
  }

  return;
}


/*------------------------------------------------------------------------------------------------*
 | private: compute normal vector to flame front patch                                henke 08/09 |
 *----------------------------------------------------------------------------------------------- */
void COMBUST::Reinitializer::ComputeNormalVectorToFlameFront(
      const GEO::BoundaryIntCell&      patch,
      const LINALG::SerialDenseMatrix& patchcoord,
      LINALG::Matrix<3,1>&             normal)
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

//  const Epetra_Comm& comm = scatra_.Discretization()->Comm();
//  cout << "proc " << comm.MyPID() << " normal " <<  normal << endl;
//  cout << "proc " << comm.MyPID() << " patch " <<  patchcoord << endl;

  // compute unit (normed) normal vector
  double norm = sqrt(normal(0)*normal(0) + normal(1)*normal(1) + normal(2)*normal(2));
  if (norm == 0.0) dserror("norm of normal vector is zero!");
  normal.Scale(1.0/norm);

#ifdef DEBUG
//  if (!((normal(2) > 0.0-1.0E-12) and (normal(2) < 0.0+1.0E-12)))
//  {
//    cout << "z-component of normal: " << normal(2) << endl;
//    dserror ("pseudo-3D problem not symmetric anymore!");
//  }
#endif

  return;
}


/*------------------------------------------------------------------------------------------------*
 | private: find facing flame front patch by projection on the patch egdes            henke 08/09 |
 *----------------------------------------------------------------------------------------------- */
void::COMBUST::Reinitializer::FindFacingPatchProjCellEdges(
       const LINALG::Matrix<3,1>&       node,
       const LINALG::SerialDenseMatrix& patch,
       bool&                            facenode,
       double&                          vertexdist)
{
  dserror("Thou shalt not call this function! It is out of order.");
#if 0
  facenode = false;
  // number of column vectors (edges) forming flame front patch
  const int numcol = patch.N();

  // beginning point of flame front patch
  LINALG::Matrix<3,1> vecbeg(true);
  // end point of flame front patch
  LINALG::Matrix<3,1> vecend(true);
  // edge vector of flame front patch
  LINALG::Matrix<3,1> edge(true);
  // distance vector from beginning point to node
  LINALG::Matrix<3,1> dist(true);

  int counter = 0;

  for(int iedge = 0; iedge<numcol; ++iedge)
  {
    // beginning point of flame front patch
    vecbeg(0) = patch(0,iedge);
    vecbeg(1) = patch(1,iedge);
    vecbeg(2) = patch(2,iedge);

    if(iedge == (numcol-1)) // last edge
    {
      // end point of edge is first point of flame front patch
      vecend(0) = patch(0,0);
      vecend(1) = patch(1,0);
      vecend(2) = patch(2,0);
    }
    else // regular edge (not the last one)
    {
      // end point of flame front patch
      vecend(0) = patch(0,iedge+1);
      vecend(1) = patch(1,iedge+1);
      vecend(2) = patch(2,iedge+1);
    }

    // compute edge vector of flame front patch
    edge.Update(1.0, vecend, -1.0, vecbeg);
    // compute distance vector from node to flame front
    dist.Update(1.0, vecend, -1.0, node);

    // compute norms of edge and distance vectors
    double normedge = sqrt(edge(0)*edge(0) + edge(1)*edge(1) + edge(2)*edge(2));
    double normdist = sqrt(dist(0)*dist(0) + dist(1)*dist(1) + dist(2)*dist(2));
    if(normdist < fabs(vertexdist))
    {
      // compute normal vector of flame front patch
      LINALG::Matrix<3,1> normal(true);
      ComputeNormalVectorToFlameFront(patch, normal);

      // determine sign of distance vector from node to flame front
      double tmp = normal(0)*dist(0) + normal(1)*dist(1) + normal(2)*dist(2);
      // tmp < 0 if node in unburnt domain (G<0) and vice versa
      if (tmp <= 0.0)
        vertexdist = -normdist;
      else
        vertexdist = normdist;

      //cout << "distance to vertex is overwritten by: " << vertexdist << endl;
    }

    // divide vectors by their norms to get unit vectors
    if (normedge == 0.0) dserror("length of edge vector is zero");
    for (int icomp=0; icomp<3; ++icomp) edge(icomp) /= normedge;
    if (normdist == 0.0) dserror("length of distance vector is zero");
    for (int icomp=0; icomp<3; ++icomp) dist(icomp) /= normdist;

    // compute projection of distance vector in direction of edge
    double proj = edge(0)*dist(0) + edge(1)*dist(1) + edge(2)*dist(2);

    // if (0 < projection < 1), node faces this edge of the patch
    if((proj >= 0.0-1.0E-13) and (proj <= 1.0+1.0E-13))
      counter += 1;
    // if all projections face their edges, the node faces the patch
    if(counter == numcol)
      facenode = true;
  }
#endif
  return;
}


#endif // #ifdef CCADISCRET
