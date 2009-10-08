/*!-----------------------------------------------------------------------------------------------*
\file combust_algorithm.cpp

\brief base combustion algorithm

    detailed description in header file combust_algorithm.H

<pre>
Maintainer: Florian Henke
            henke@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15265
</pre>
 *------------------------------------------------------------------------------------------------*/
#ifdef CCADISCRET

#include "combust_algorithm.H"
#include "combust_flamefront.H"
#include "combust_utils.H"
#include "../drt_lib/drt_globalproblem.H"
#include "combust_fluidimplicitintegration.H"
#include "../drt_mat/newtonianfluid.H"
#include "../drt_mat/matlist.H"
#include "../drt_geometry/integrationcell.H"
#include "../drt_lib/drt_function.H"

/*------------------------------------------------------------------------------------------------*
 | constructor                                                                        henke 06/08 |
 *------------------------------------------------------------------------------------------------*/
  /* Der constructor sollte den gesamten Algorithmus initialisieren.
   * Hier müssen alle Variablen, die den Einzelfeldern übergeordnet sind, initialisiert werden.
   *
   * das heisst:
   * o G-Funktionsvektor (t und ig+1) auf initial value setzen
   * o Geschwindigkeitsvektor (t und iu+1) auf initial value setzen
   * o alle Zähler auf 0 setzen (step_(0), f_giter_(0), g_iter_(0), f_iter_(0))
   * o alle Normen und Grenzwerte auf 0 setzen
   *
   * Zusammenfassend kann an sagen, dass alles was in der Verbrennungsrechnung vor der Zeitschleife
   * passieren soll, hier passieren muss, weil die combust dyn gleich die Zeitschleife ruft.
   *
   * scalar transport velocity field has been initialized in ScaTraFluidCouplingAlgorithm()
  */
COMBUST::Algorithm::Algorithm(Epetra_Comm& comm, const Teuchos::ParameterList& combustdyn)
: ScaTraFluidCouplingAlgorithm(comm, combustdyn,false),
// initialize member variables
  fgiter_(0),
  fgitermax_(combustdyn.get<int>("ITEMAX")),
  convtol_(combustdyn.get<double>("CONVTOL")),
//  fgvelnormL2_(?),
//  fggfuncnormL2_(?),
/* hier müssen noch mehr Parameter eingelesen werden, nämlich genau die, die für die FGI-Schleife
   nötig sind. Die für die Schleifen über die Einzelfelder existieren in den Zeitintegrationsschemata. */
//  reinitializationaction_(combustdyn.sublist("COMBUSTION GFUNCTION").get<INPAR::COMBUST::ReInitialActionGfunc>("REINITIALIZATION")),
  combusttype_(Teuchos::getIntegralValue<INPAR::COMBUST::CombustionType>(combustdyn.sublist("COMBUSTION FLUID"),"COMBUSTTYPE")),
  reinitializationaction_(Teuchos::getIntegralValue<INPAR::COMBUST::ReInitialActionGfunc>(combustdyn.sublist("COMBUSTION GFUNCTION"),"REINITIALIZATION")),
  reinitinterval_(combustdyn.sublist("COMBUSTION GFUNCTION").get<int>("REINITINTERVAL")),
  combustdyn_(combustdyn),
  interfacehandle_(Teuchos::null),
  flamefront_(Teuchos::null)
  {

  // get pointers to the discretizations from the time integration scheme of each field
  /* remark: fluiddis cannot be of type "const Teuchos::RCP<const DRT::Dis...>", because parent
   * class. InterfaceHandle only accepts "const Teuchos::RCP<DRT::Dis...>"            henke 01/09 */
  const Teuchos::RCP<DRT::Discretization> fluiddis = FluidField().Discretization();
  const Teuchos::RCP<DRT::Discretization> gfuncdis = ScaTraField().Discretization();

  /*----------------------------------------------------------------------------------------------*
   * initialize all data structures needed for the combustion algorithm
   *
   * - capture the flame front and create interface geometry (triangulation)
   * - determine initial enrichment (DofManager wird bereits mit dem Element d.h. Diskretisierung angelegt)
   * - ...
   *----------------------------------------------------------------------------------------------*/
  // construct initial flame front
  flamefront_ = rcp(new COMBUST::FlameFront(fluiddis,gfuncdis));
  flamefront_->ProcessFlameFront(combustdyn_,ScaTraField().Phinp());

  // construct interfacehandle using initial flame front
  interfacehandle_ = rcp(new COMBUST::InterfaceHandleCombust(fluiddis,gfuncdis,flamefront_));
  // get integration cells according to initial flame front
  interfacehandle_->UpdateInterfaceHandle();

  velnpip_ = rcp(new Epetra_Vector(*fluiddis->DofRowMap()),true);
  velnpi_ = rcp(new Epetra_Vector(*fluiddis->DofRowMap()),true);

  phinpip_ = rcp(new Epetra_Vector(*gfuncdis->NodeRowMap()),true);
  phinpi_ = rcp(new Epetra_Vector(*gfuncdis->NodeRowMap()),true);

  if (Comm().MyPID()==0)
  {
    std::cout << "\n-------------------------  Initial state of coupled problem set  -----------------------------" << std::endl;
  }
  // std::cout << "Combustion Algorithm constructor done \n" << endl;
}

/*------------------------------------------------------------------------------------------------*
 | destructor                                                                         henke 06/08 |
 *------------------------------------------------------------------------------------------------*/
COMBUST::Algorithm::~Algorithm()
{
}

/*------------------------------------------------------------------------------------------------*
 | public: time loop of algorithm for dynamic combustion problem                      henke 06/08 |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::Algorithm::TimeLoop()
{

  if (Comm().MyPID()==0)
    std::cout << "Combustion Algorithm Timeloop starting \n" << endl;

//URSULA
  // initial volume of domain minus
  const double volume_start = interfacehandle_->ComputeVolumeMinus();

  // calling Output() here causes an error with the velpressplitterForOutput_ in TransformXFEMToStandardVector()
//  Output();

  // time loop
  while (NotFinished())
  {
    // prepare next time step
    PrepareTimeStep();

    // Fluid-G-function-Interaction loop
    while (NotConvergedFGI())
    {
      // prepare Fluid-G-function iteration
      PrepareFGIteration();

      // solve nonlinear Navier-Stokes system
      DoFluidField();

      // solve linear G-function equation
      DoGfuncField();

      // update field vectors
      UpdateFGIteration();

    } // Fluid-G-function-Interaction loop

    if (ScaTraField().Step() % reinitinterval_ == 0)
    {
      // reinitialize G-function
      ReinitializeGfunc(reinitializationaction_);
    }

    // update all field solvers
    UpdateTimeStep();

    // write output to screen and files
    Output();

  const double volume_current = interfacehandle_->ComputeVolumeMinus();
  // compute mass loss
  double mass_current = fabs(volume_current - volume_start) / volume_start;
  if (Comm().MyPID()==0)
  {
    std::cout << "\n=======================================\n" << endl;
    std::cout << "              Mass calculation         \n" << endl;
    std::cout << "initial mass: " << volume_start << endl;
    std::cout << "final mass:   " << volume_current << endl;
    std::cout << "mass change:  " << mass_current << endl;
    std::cout << "\n=======================================\n" << endl;
  }

  } // time loop

//URSULA
  // final volume of domain minus
  const double volume_end = interfacehandle_->ComputeVolumeMinus();
  // compute mass loss
  double mass_change = fabs(volume_end - volume_start) / volume_start;
  if (Comm().MyPID()==0)
  {
	  std::cout << "\n=======================================\n" << endl;
	  std::cout << "              Mass calculation         \n" << endl;
	  std::cout << "initial mass: " << volume_start << endl;
	  std::cout << "final mass:   " << volume_end << endl;
	  std::cout << "mass change:  " << mass_change << endl;
	  std::cout << "\n=======================================\n" << endl;
  }

  return;
} // TimeLoop()

/*------------------------------------------------------------------------------------------------*
 | public: algorithm for static combustion problem                                    henke 08/08 |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::Algorithm::SolveStationaryProblem()
{
  dserror("Der Algorithmus für statische Verbrennungprobleme kann noch nix!");
  return;
}

/*------------------------------------------------------------------------------------------------*
 | protected: reinitialize G-function                                                 henke 06/08 |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::Algorithm::ReinitializeGfunc(INPAR::COMBUST::ReInitialActionGfunc action)
{
  /* Here, the G-function is reinitialized, because we suspect that the signed distance property has
   * been lost due to numerical inaccuracies. There are various options to reinitialize the
   * G-function. For the time being we use an exact, but expensive, procedure to ensure this
   * property which computes the distance for every point. (SignedDistFunc)           henke 06/08 */
  if (Comm().MyPID()==0)
  {
    cout<<"\n------------------------------------  REINITIALIZE G-FUNCTION  -------------------------------\n\n";
  }

  switch (action)
  {
    case INPAR::COMBUST::reinitaction_none:
      // do nothing
      break;
    case INPAR::COMBUST::reinitaction_byfunction:
      // read a FUNCTION from the input file and reinitialize the G-function bz evaluating it
      ScaTraField().SetInitialField(INPAR::SCATRA::initfield_field_by_function,combustdyn_.sublist("COMBUSTION GFUNCTION").get<int>("REINITFUNCNO"));
      break;
    case INPAR::COMBUST::reinitaction_signeddistancefunction:
#ifdef PARALLEL
  dserror("direct computation of signed distance function not available in parallel");
#endif
      SignedDistFunc();
      break;
//    case INPAR::COMBUST::sussman:
//      // SCATRA parameters have to been set before in ScaTraFluidCouplingAlgorithm!
//      ScaTraField().Reinitialize();
//      break;
    default:
      dserror ("Unknown option to reinitialize the G-function");
  }
  return;
}

/*------------------------------------------------------------------------------------------------*
 | protected: build signed distance function                                          henke 06/08 |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::Algorithm::SignedDistFunc()
{
  /* This member function constructs a G-function field that meets the signed distance property. The
   * algorithm assigns the value of the distance between each node and the surface defined by G=0 as
   * a scalar value to every node in the G-function discretization.*/

  // get a pointer to the G-function discretization
  Teuchos::RCP<DRT::Discretization> gfuncdis = ScaTraField().Discretization();

  const Epetra_Map* dofrowmap = gfuncdis->DofRowMap();

  // loop all row nodes on the processor
  for(int lnodeid=0; lnodeid < gfuncdis->NumMyRowNodes(); ++lnodeid)
  {
    // get the processor local node
    const DRT::Node*  lnode      = gfuncdis->lRowNode(lnodeid);

    // the set of degrees of freedom associated with the node
    const vector<int> nodedofset = gfuncdis->Dof(lnode); // this should not be a vector, but a scalar

#ifdef DEBUG
    int numdofs = nodedofset.size(); // this should be 1 (scalar field!)
    if (numdofs != 1)
      dserror("There are more than 1 dof for a node in a scalar field!");
#endif

    const int dofgid = nodedofset[0];
    int doflid = dofrowmap->LID(dofgid);

    //-------------------------------------------------------
    // compute smallest distance to flame front for this node
    //-------------------------------------------------------
    //cout << "******************************************" << endl;
    //cout << "reinitialization for node: " << lnode->Id() << endl;
    //cout << "******************************************" << endl;
    // smallest distance to the vertex of a flame front patch
    double vertexdist = 5555.5;
    // smallest distance to flame front
    double mindist = 7777.7;

    // get physical coordinates of this node
    LINALG::Matrix<3,1> nodecoord;
    nodecoord(0) = lnode->X()[0];
    nodecoord(1) = lnode->X()[1];
    nodecoord(2) = lnode->X()[2];

    // get boundary integration cells (flame front patches)
    const std::map<int,GEO::BoundaryIntCells>& flamefront = interfacehandle_->GetElementalBoundaryIntCells();

    //-----------------------------------------
    // find flame front patches facing the node
    //-----------------------------------------
    // loop groups (vectors) of flamefront patches of all elements
    for(std::map<int,GEO::BoundaryIntCells>::const_iterator elepatches = flamefront.begin(); elepatches != flamefront.end(); ++elepatches)
    {
      // number of flamefront patches for this element
      const std::vector<GEO::BoundaryIntCell> patches = elepatches->second;
      const int numpatch = patches.size();

      // loop flame front patches of this element
      for(int ipatch=0; ipatch<numpatch; ++ipatch)
      {
        // get coordinates of vertices of flame front patch
        const LINALG::SerialDenseMatrix&  patchcoord = patches[ipatch].CellNodalPosXYZ();

        // check if this patch faces this node
        bool facenode = false;
        FindFacingPatch(nodecoord,patchcoord,facenode,vertexdist);

        // a facing patch was found
        if (facenode == true)
        {
          //-------------------------------------------
          // compute distance to this flame front patch
          //-------------------------------------------
          double patchdist;
          ComputeDistanceToFlameFront(nodecoord,patchcoord,patchdist);

          // overwrite smallest distance if computed patch distance is smaller
          if (fabs(patchdist) < fabs(mindist))
          {
            //cout << "distance to flame front patch: " << mindist << " is overwritten by: " << patchdist << endl;
            mindist = patchdist;
          }
        }
      }
    }

    //-------------------------------------------
    // determine smallest distance to flame front
    //-------------------------------------------
    // remark: variable have the following meaning here:
    //         - "mindist" is either "patchdist" or still the default value (7777.7)
    //         - "vertexdist" is the distance to the clostest vertex of any flame front patch
    //
    // case 1: a local flame front patch was found for this node -> mindist = patchdist
    // case 2: a flame front patch was found for this node, but it is not local (curved interface);
    //         a local patch was not found, because this node is located in the blind angle of all
    //         local patches -> mindist = vertexdist
    // case 3: something went wrong:  "mindist" still has default value (7777.7)

    //cout << "mindist " << std::setw(18) << std::setprecision(12) << std::scientific << mindist << endl;
    //cout << "vertexdist " << std::setw(18) << std::setprecision(12) << std::scientific << vertexdist << endl;
    if (fabs(mindist) > fabs(vertexdist)) // case 2
    {
      mindist = vertexdist;

      //if ((nodecoord(0) == 0.5) or (nodecoord(0) == -0.5) or
      //    (nodecoord(1) == 0.5) or (nodecoord(1) == -0.5))
      //{
      //  cout << "Boundary node: ";
      //}
      cout << "node " << lnode->Id() << " does not face a patch (blind angle) or this patch is not local; distance: " << mindist << endl;
    }
    if (mindist == 7777.7) // case 3
      dserror ("node %d was not reinitialized!", lnode->Id());

    //------------------------------
    // reinitialize G-function field
    //------------------------------
    // assign new values of signed distance function to the G-function field
    ScaTraField().Phinp()->ReplaceMyValues(1,&mindist,&doflid);
  }

  return;
}

/*------------------------------------------------------------------------------------------------*
 | compute normal vector of flame front patch                                         henke 08/09 |
 *----------------------------------------------------------------------------------------------- */
void COMBUST::Algorithm::ComputeNormalVectorToFlameFront(
      const LINALG::SerialDenseMatrix& patch,
      LINALG::Matrix<3,1>&             normal)
{
  // first point of flame front patch
  LINALG::Matrix<3,1> point1;
  point1(0) = patch(0,0);
  point1(1) = patch(1,0);
  point1(2) = patch(2,0);

  // second point of flame front patch
  LINALG::Matrix<3,1> point2;
  point2(0) = patch(0,1);
  point2(1) = patch(1,1);
  point2(2) = patch(2,1);

  // first edge of flame front patch
  LINALG::Matrix<3,1> edge1;
  edge1.Update(1.0, point2, -1.0, point1);

  // third point of flame front patch
  point2(0) = patch(0,2);
  point2(1) = patch(1,2);
  point2(2) = patch(2,2);

  // second edge of flame front patch (if patch is triangle; if not: edge 2 is secant of polygon)
  LINALG::Matrix<3,1> edge2;
  edge2.Update(1.0, point2, -1.0, point1);

  // compute normal vector of patch (cross product: edge1 x edge2)
  // remark: normal vector points into unburnt domain (G<0)
  normal(0) = (edge1(1)*edge2(2) - edge1(2)*edge2(1));
  normal(1) = (edge1(2)*edge2(0) - edge1(0)*edge2(2));
  normal(2) = (edge1(0)*edge2(1) - edge1(1)*edge2(0));

  // compute unit (normed) normal vector
  double norm = sqrt(normal(0)*normal(0) + normal(1)*normal(1) + normal(2)*normal(2));
  if (norm == 0.0) dserror("norm of normal vector is zero!");
  normal.Scale(1.0/norm);

#ifdef DEBUG
  if (!((normal(2) > 0.0-1.0E-12) and (normal(2) < 0.0+1.0E-12)))
  {
    cout << "z-component of normal: " << normal(2) << endl;
    dserror ("pseudo-3D problem not symmetric anymore!");
  }
#endif

  return;
}
/*------------------------------------------------------------------------------------------------*
 | compute distance directly                                                          henke 08/09 |
 *----------------------------------------------------------------------------------------------- */
void COMBUST::Algorithm::ComputeDistanceToFlameFront(
      const LINALG::Matrix<3,1>&       node,
      const LINALG::SerialDenseMatrix& patch,
      double&                          distance)
{
  // first point of flame front patch
  LINALG::Matrix<3,1> point;
  point(0) = patch(0,0);
  point(1) = patch(1,0);
  point(2) = patch(2,0);

  // compute distance vector between node and point on flame front patch (point)
  // remark: distance vector points towards the flame front
  LINALG::Matrix<3,1> dist;
  dist.Update(1.0, point, -1.0, node);

  // compute normal vector of flame front patch
  LINALG::Matrix<3,1> normal(true);
  ComputeNormalVectorToFlameFront(patch, normal);

  // project distance vector in normal direction: compute scalar product (normal * dist)
  // remark: distance < 0 if node in unburnt domain (G<0) and vice versa
  distance = normal(0)*dist(0) + normal(1)*dist(1) + normal(2)*dist(2);

  return;
}

void::COMBUST::Algorithm::FindFacingPatch(
       const LINALG::Matrix<3,1>&       node,
       const LINALG::SerialDenseMatrix& patch,
       bool&                            facenode,
       double&                          vertexdist)
{
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
  return;
}

/*------------------------------------------------------------------------------------------------*
 | protected: overwrite Navier-Stokes velocity                                        henke 08/09 |
 *------------------------------------------------------------------------------------------------*/
const Teuchos::RCP<Epetra_Vector> COMBUST::Algorithm::OverwriteFluidVel()
{
  //----------------------------------------------------------------
  // impose velocity field by function e.g. for level set test cases
  // Navier-Stokes solution velocity field is overwritten
  //----------------------------------------------------------------
  if (Comm().MyPID()==0)
    cout << "\n\n===================== overwriting Navier-Stokes solution ======================\n\n" << endl;

  // get fluid (Navier-Stokes) velocity vector in standard FEM configuration (no XFEM dofs)
  const Teuchos::RCP<Epetra_Vector> convel = FluidField().ExtractInterfaceVeln();
  // velocity function number = 1 (FUNCT1)
  const int velfuncno = 1;

  // loop all nodes on the processor
  for(int lnodeid=0; lnodeid < FluidField().Discretization()->NumMyRowNodes(); ++lnodeid)
  {
    // get the processor local node
    DRT::Node*  lnode = FluidField().Discretization()->lRowNode(lnodeid);
    // get standard dofset from fluid time integration
    vector<int> fluidnodedofs = (*(FluidField().DofSet())).Dof(lnode);
    // determine number of space dimensions (numdof - pressure dof)
    const int numdim = ((int) fluidnodedofs.size()) -1;
    if (numdim != 3) dserror("3 components expected for velocity");

    int err = 0;
    // overwrite velocity dofs only
    for(int index=0; index<numdim; ++index)
    {
      // global and processor's local fluid dof ID
      const int fgid = fluidnodedofs[index];
      int flid = convel->Map().LID(fgid);
      if (flid < 0) dserror("lid not found in map for given gid");

      // get value of corresponding velocity component
      double value = DRT::Problem::Instance()->Funct(velfuncno-1).Evaluate(index,lnode->X(),0.0,NULL);

      // insert velocity value into node-based vector
//        value = 0.0;
      err += convel->ReplaceMyValues(1, &value, &flid);
    }
    if (err != 0) dserror("error overwriting Navier-Stokes solution");
  }

  return convel;
}

/*------------------------------------------------------------------------------------------------*
 | protected: compute flame velocity                                                  henke 07/09 |
 *------------------------------------------------------------------------------------------------*/
const Teuchos::RCP<Epetra_Vector> COMBUST::Algorithm::ComputeFlameVel(const Teuchos::RCP<Epetra_Vector>& convel,
                                                                      const Teuchos::RCP<const DRT::DofSet>& dofset)
{
  // get a pointer to the fluid discretization
  const Teuchos::RCP<const DRT::Discretization> fluiddis = FluidField().Discretization();
  // get G-function value vector on fluid NodeColMap
  const Teuchos::RCP<Epetra_Vector> phinp = flamefront_->Phinp();

#ifdef DEBUG
  // get map of this vector
  const Epetra_BlockMap& phimap = phinp->Map();
  // check, whether this map is still identical with the current node map in the discretization
  if (not phimap.SameAs(*fluiddis->NodeColMap())) dserror("node column map has changed!");
#endif

  // loop over nodes on this processor
  for(int lnodeid=0; lnodeid < fluiddis->NumMyRowNodes(); ++lnodeid)
  {
    // get the processor local node
    DRT::Node*  lnode      = fluiddis->lRowNode(lnodeid);
    // get list of adjacent elements of this node
    DRT::Element** elelist = lnode->Elements();

//    cout << "------------------------------------------------------------" << endl;
//    cout << "run for node: " << lnode->Id() << endl;
//    cout << "------------------------------------------------------------" << endl;

    //--------------------------------------------------------
    // compute "average"/"smoothed" normal vector at this node
    //--------------------------------------------------------
    // vector for average normal vector at this node
    LINALG::Matrix<3,1> avnvec(true);

    // loop over adjacent elements of this node
    for(int iele=0; iele<lnode->NumElement(); ++iele)
    {
      const DRT::Element* ele = elelist[iele];
      const int numnode = ele->NumNode();

//    cout << "------------------------------------------------------------" << endl;
//    cout << "run for element: " << ele->Id() << endl;
//    cout << "------------------------------------------------------------" << endl;

      // extract G-function values for nodes of this element
      Epetra_SerialDenseMatrix myphi(numnode,1);
      DRT::UTILS::ExtractMyNodeBasedValues(ele, myphi, *phinp);

      // get node coordinates of this element
      Epetra_SerialDenseMatrix xyze(3,numnode);
      GEO::fillInitialPositionArray<Epetra_SerialDenseMatrix>(ele, xyze);

      // TODO: function could be templated DISTYPE -> LINALG::Matrix<3,DISTYPE>
      Epetra_SerialDenseMatrix deriv(3,numnode);

#ifdef DEBUG
      bool nodefound = false;
#endif
      // find out which node in the element is my local node lnode
      for (int inode=0; inode<numnode; ++inode)
      {
        if (ele->NodeIds()[inode] == lnode->Id())
        {
          // get local (element) coordinates of this node
          LINALG::Matrix<3,1> coord = DRT::UTILS::getNodeCoordinates(inode,ele->Shape());
          // evaluate derivatives of shape functions at this node
          DRT::UTILS::shape_function_3D_deriv1(deriv,coord(0),coord(1),coord(2),ele->Shape());
#ifdef DEBUG
          nodefound = true;
#endif
        }
      }
#ifdef DEBUG
      if (nodefound==false)
        dserror("node was not found in list of elements");
#endif
      //----------------------------------------------------
      // compute normal vector at this node for this element
      // n = - grad phi / |grad phi|
      //----------------------------------------------------
      // evaluate gradient of G-function field at this node
      // remark: grad phi = sum (grad N_i * phi_i)

      // get transposed of the jacobian matrix d x / d \xi
      // xjm(i,j) = deriv(i,k)*xyze(j,k)
      Epetra_SerialDenseMatrix xjm(3,3);
      xjm.Multiply('N','T',1.0,deriv,xyze,0.0);

      // inverse of jacobian (xjm)
//      Epetra_SerialDenseMatrix xji(3,3);
//      xjm.Invert(xji);

      LINALG::NonSymmetricInverse(xjm,3);

      // compute global derivates
      Epetra_SerialDenseMatrix derxy(3,numnode);
      // derxy(i,j) = xji(i,k) * deriv(k,j)
      derxy.Multiply('N','N',1.0,xjm,deriv,0.0);
//      xji.Multiply(false,deriv,derxy);

      Epetra_SerialDenseMatrix gradphi(3,1);
      derxy.Multiply(false,myphi,gradphi);
      double ngradphi = sqrt(gradphi(0,0)*gradphi(0,0)+gradphi(1,0)*gradphi(1,0)+gradphi(2,0)*gradphi(2,0));

      // normal vector at this node for this element
      LINALG::Matrix<3,1> nvec(true);
      if (ngradphi == 0.0) dserror("length of normal is zero");
      for (int icomp=0; icomp<3; ++icomp)
        nvec(icomp) = -gradphi(icomp,0) / ngradphi;

//      cout << "normal vector for element: " << ele->Id() << " at node: " << lnode->Id() << " is: " << nvec << endl;

      // add normal vector to linear combination (could also be weighted in different ways!)
      for (int icomp=0; icomp<3; ++icomp)
        avnvec(icomp) = avnvec(icomp) + nvec(icomp);
    }
    //---------------------------
    // compute unit normal vector
    //---------------------------
    // compute norm of average normal vector
    double avnorm = sqrt(avnvec(0)*avnvec(0)+avnvec(1)*avnvec(1)+avnvec(2)*avnvec(2));
    // divide vector by its norm to get unit normal vector
    if (avnorm == 0.0) dserror("length of normal is zero");
    for (int icomp=0; icomp<3; ++icomp) avnvec(icomp) /= avnorm;

//    cout << "average normal vector at node: " << lnode->Id() << " is: " << avnvec << endl;

    //------------------------
    // get material parameters
    //------------------------
    // get material from first (arbitrary!) element adjacent to this node
    const Teuchos::RCP<MAT::Material> matlistptr = elelist[0]->Material();
    dsassert(matlistptr->MaterialType() == INPAR::MAT::m_matlist, "material is not of type m_matlist");
    const MAT::MatList* matlist = static_cast<const MAT::MatList*>(matlistptr.get());

    // density burnt domain
    Teuchos::RCP<const MAT::Material> matptrplus = matlist->MaterialById(3);
    dsassert(matptrplus->MaterialType() == INPAR::MAT::m_fluid, "material is not of type m_fluid");
    const MAT::NewtonianFluid* matplus = static_cast<const MAT::NewtonianFluid*>(matptrplus.get());
    const double rhoplus = matplus->Density();

    // density unburnt domain
    Teuchos::RCP<const MAT::Material> matptrminus = matlist->MaterialById(4);
    dsassert(matptrminus->MaterialType() == INPAR::MAT::m_fluid, "material is not of type m_fluid");
    const MAT::NewtonianFluid* matminus = static_cast<const MAT::NewtonianFluid*>(matptrminus.get());
    const double rhominus = matminus->Density();

    // laminar flame speed
    const double sl = combustdyn_.sublist("COMBUSTION FLUID").get<double>("LAMINAR_FLAMESPEED");
    //---------------------------------------------
    // compute relative flame velocity at this node
    //---------------------------------------------
    // get phi value for this node
    const int lid = phinp->Map().LID(lnode->Id());
    const double gfuncval = (*phinp)[lid];

    double speedfac = 0.0;
    if (gfuncval > 0.0) // burnt domain -> burnt material
      // flame speed factor = laminar flame speed * rho_unburnt / rho_burnt
      speedfac = sl * rhominus/rhoplus;
    else // interface or unburnt domain -> unburnt material
      // flame speed factor = laminar flame speed
      speedfac = sl;

    LINALG::Matrix<3,1> flvelrel(true);
    for (int icomp=0; icomp<3; ++icomp)
      flvelrel(icomp) = speedfac * avnvec(icomp);
    //-----------------------------------------------
    // compute (absolute) flame velocity at this node
    //-----------------------------------------------
    LINALG::Matrix<3,1> fluidvel(true);
    // get the set of dof IDs for this node (3 x vel + 1 x pressure) from standard FEM dofset
    const std::vector<int> dofids = (*dofset).Dof(lnode);
    //
    std::vector<int> lids(3);

    // extract velocity values (no pressure!) from global velocity vector
    for (int icomp=0; icomp<3; ++icomp)
    {
      lids[icomp] = convel->Map().LID(dofids[icomp]);
      fluidvel(icomp) = (*convel)[lids[icomp]];
    }

    LINALG::Matrix<3,1> flvelabs(true);
    // add fluid velocity (Navier Stokes solution) and relative flame velocity
    for (int icomp=0; icomp<3; ++icomp)
    {
      flvelabs(icomp) = fluidvel(icomp) + flvelrel(icomp);
      convel->ReplaceMyValues(1,&flvelabs(icomp),&lids[icomp]);
    }
//  cout << "------------------------------------------------------------" << endl;
//  cout << "run for node: " << lnode->Id() << endl;
//  cout << "------------------------------------------------------------" << endl;
//  cout << "convection velocity: " << flvelabs << endl;
  }
  return convel;
}

/*------------------------------------------------------------------------------------------------*
 | protected: FGI iteration converged?                                                henke 06/08 |
 *------------------------------------------------------------------------------------------------*/
bool COMBUST::Algorithm::NotConvergedFGI()
{
  // if (fgiter <= fgitermax and ComputeGfuncNorm() < maxepsg and ComputeFluidNorm() < maxepsf)
  //return (fgiter_ < fgitermax_ and true);

  bool notconverged = true;

  if (fgiter_ <= fgitermax_)
  {
      if (fgiter_ == 0)
      {

      }
      else if (fgiter_ == 1)
      {
    	  if (velnpip_->MyLength() != FluidField().ExtractInterfaceVeln()->MyLength())
    	    dserror("vectors must have the same length");
    	  velnpip_->Update(1.0,*(FluidField().ExtractInterfaceVeln()),0.0);
    	  phinpip_->Update(1.0,*(flamefront_->Phinp()),0.0);

    	  if (fgiter_ == fgitermax_)
    		  notconverged = false;
      }
      else
      {
    	  velnpi_->Update(1.0,*velnpip_,0.0);
    	  phinpi_->Update(1.0,*phinpip_,0.0);
    	  velnpip_->Update(1.0,*(FluidField().ExtractInterfaceVeln()),0.0);
    	  phinpip_->Update(1.0,*(flamefront_->Phinp()),0.0);

    	  fgvelnormL2_ = 1.0;
    	  fggfuncnormL2_ = 1.0;

    	  Teuchos::RCP<Epetra_Vector> incvel = rcp(new Epetra_Vector(velnpip_->Map()),true);
    	  if (incvel->MyLength() != FluidField().ExtractInterfaceVeln()->MyLength())
    	    dserror("vectors must have the same length");
    	  incvel->Update(1.0,*velnpip_,-1.0,*velnpi_,0.0);
    	  incvel->Norm2(&fgvelnormL2_);

    	  Teuchos::RCP<Epetra_Vector> incgfunc = rcp(new Epetra_Vector(*ScaTraField().Discretization()->NodeRowMap()),true);
    	  incgfunc->Update(1.0,*phinpip_,-1.0,*phinpi_,0.0);
    	  incgfunc->Norm2(&fggfuncnormL2_);

    	  if (Comm().MyPID()==0)
    	  {
    	     printf("\n|+---------------------- FGI ----------------------+|");
    	     printf("\n|iter/itermax|----tol-----|-fluid inc--|-g-func inc-|");
    	     printf("\n|   %2d/%2d    | %10.3E | %10.3E | %10.3E |",fgiter_,fgitermax_,convtol_,fgvelnormL2_,fggfuncnormL2_);
    	     printf("\n|+-------------------------------------------------+|\n");
    	  }

    	  if ((fgvelnormL2_ <= convtol_) and (fggfuncnormL2_ <= convtol_))
    	  {
             notconverged = false;
    	  }
          else
          {
             if (fgiter_ == fgitermax_)
             {
            	 notconverged = false;
            	 if (Comm().MyPID()==0)
            	 {
            	    printf("|+---------------- not converged ------------------+|");
            	    printf("\n|+-------------------------------------------------+|");
            	 }
             }
          }
      }
  }
  else
  {

  }

  return notconverged;
}

/*------------------------------------------------------------------------------------------------*
 | protected: prepare time step for combustion algorithm                              henke 06/08 |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::Algorithm::PrepareTimeStep()
{
  IncrementTimeAndStep();
  fgiter_ = 0;
  // fgnormgfunc = large value, determined in Input file
  fgvelnormL2_ = 1.0;
  // fgnormfluid = large value
  fggfuncnormL2_ = 1.0;

  if (Comm().MyPID()==0)
  {
    //cout<<"---------------------------------------  time step  ------------------------------------------\n";
    printf("----------------------Combustion-------  time step %2d ----------------------------------------\n",Step());
    printf("TIME: %11.4E/%11.4E  DT = %11.4E STEP = %4d/%4d \n",Time(),MaxTime(),Dt(),Step(),NStep());
  }

  FluidField().PrepareTimeStep();

  // prepare time step
  /* remark: initial velocity field has been transfered to scalar transport field in constructor of
   * ScaTraFluidCouplingAlgorithm (initialvelset_ == true). Time integration schemes, such as
   * the one-step-theta scheme, are thus initialized correctly.
   */
  ScaTraField().PrepareTimeStep();

  // synchronicity check between combust algorithm and base algorithms
  if (FluidField().Time() != Time())
    dserror("Time in Fluid time integration differs from time in combustion algorithm");
  if (ScaTraField().Time() != Time())
    dserror("Time in ScaTra time integration  differs from time in combustion algorithm");

  return;
}

/*------------------------------------------------------------------------------------------------*
 | protected: prepare time step for combustion algorithm                              henke 06/08 |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::Algorithm::PrepareFGIteration()
{
  fgiter_ += 1;
  if (Comm().MyPID()==0)
  {
    //cout<<"\n---------------------------------------  FGI loop  -------------------------------------------\n";
    printf("\n---------------------------------------  FGI loop: iteration number: %2d ----------------------\n",fgiter_);
  }
}

/*------------------------------------------------------------------------------------------------*
 | protected: perform a fluid time integration step                                   henke 06/08 |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::Algorithm::DoFluidField()
{
  if (Comm().MyPID()==0)
  {
    std:: cout<<"\n---------------------------------------  FLUID SOLVER  ---------------------------------------" << std::endl;
  }

  // export interface information to the fluid time integration
  FluidField().ImportInterface(interfacehandle_);

  // solve nonlinear Navier-Stokes equations
  FluidField().NonlinearSolve();

  return;
}

/*------------------------------------------------------------------------------------------------*
 | protected: perform a G-function time integration step                              henke 06/08 |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::Algorithm::DoGfuncField()
{
  if (Comm().MyPID()==0)
  {
    cout<<"\n---------------------------------------  G-FUNCTION SOLVER  ----------------------------------\n";
  }

  /* Physikalisch hat nur die Geschwindigkeit auf der Flame (G=0) eine Bedeutung für das Feld der
   * G-Function. Tatsächlich kann aber das gesamte Navier-Stokes Fluid Feld als Input genommen werden.
   * Die G-function fragt an jedem Knoten nach der Konvektionsgeschwindigkeit.
   */
  // assign the fluid velocity to the G-function field as convective velocity
  switch(combusttype_)
  {
  case INPAR::COMBUST::combusttype_twophaseflow:
  {
    // for two-phase flow, the fluid velocity field is continuous; it can be directly transferred to
    // the scalar transport field

    ScaTraField().SetVelocityField(
      FluidField().ExtractInterfaceVeln(),
//      OverwriteFluidVel(),
      Teuchos::null,
      FluidField().DofSet(),
      FluidField().Discretization()
    );
    break;
  }
  case INPAR::COMBUST::combusttype_premixedcombustion:
  {
    // for combustion, the velocity field is discontinuous; the relative flame velocity is added
    const Teuchos::RCP<Epetra_Vector> convel = FluidField().ExtractInterfaceVeln();

#if 0
//TEST
//    const Teuchos::RCP<Epetra_Vector> xfemvel = ComputeFlameVel(convel,FluidField().DofSet());
//    if((*convel).MyLength() != (*xfemvel).MyLength())
//      dserror("length is not the same!");
//
//    const int dim = (*convel).MyLength();
//    int counter = 0;
//
//    for(int idof=0; idof < dim ;++idof)
//    {
//      if ((*convel)[idof] == (*xfemvel)[idof])
//        counter++;
//    }
//    cout << "So viele dofs sind gleich: " << counter << endl;
//TEST
#endif

    ScaTraField().SetVelocityField(
      ComputeFlameVel(convel,FluidField().DofSet()),
      Teuchos::null,
      FluidField().DofSet(),
      FluidField().Discretization()
    );
    break;
  }
  default:
    dserror("unknown type of combustion problem");
  }

  // solve nonlinear convection-diffusion equation
  ScaTraField().NonlinearSolve();
  // solve linear convection-diffusion equation
//  ScaTraField().Solve();

  return;
}

/*------------------------------------------------------------------------------------------------*
 | protected: update                                                                  henke 06/08 |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::Algorithm::UpdateFGIteration()
{
  // update flame front according to evolved G-function field
  flamefront_->ProcessFlameFront(combustdyn_,ScaTraField().Phinp());

  // update interfacehandle (get integration cells) according to updated flame front
  interfacehandle_->UpdateInterfaceHandle();

  // update the Fluid and the FGI vector at the end of the FGI loop
  return;
}

/*------------------------------------------------------------------------------------------------*
 | protected: update                                                                  henke 06/08 |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::Algorithm::UpdateTimeStep()
{
  FluidField().Update();
  ScaTraField().Update();
  return;
}

/*------------------------------------------------------------------------------------------------*
 | protected: output                                                                  henke 06/08 |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::Algorithm::Output()
{
  // Note: The order is important here! In here control file entries are
  // written. And these entries define the order in which the filters handle
  // the Discretizations, which in turn defines the dof number ordering of the
  // Discretizations.
  FluidField().Output();
  FluidField().LiftDrag();
  ScaTraField().Output();

  // debug IO
#if 0
  // print out convective velocity vector
  cout<<*velocitynp_<<endl;
#endif

  return;
}

#endif // #ifdef CCADISCRET
