/*!-----------------------------------------------------------------------------------------------*
 \file combust_flamefront.cpp

 \brief

  detailed description in header file combust_interface.H

\level 2

<pre>
\maintainer Benedikt Schott
            schott@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>
 *------------------------------------------------------------------------------------------------*/


#include "combust_flamefront.H"
#include "combust_defines.H"
#include "combust_refinementcell.H"
#include "combust3_utils.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_geometry/tetrahedradecomposition.H"
#include "../drt_io/io_gmsh.H"
#include "../drt_io/io_pstream.H"
#include "../linalg/linalg_utils.H"
#include <Teuchos_TimeMonitor.hpp>
#include "../drt_geometry/element_coordtrafo.H"
#include "../drt_geometry/position_array.H"
#include "../drt_geometry/integrationcell.H"
#include "../drt_geometry/intersection_service_templates.H"
#include "../drt_cut/cut_levelsetintersection.H"
#include "../drt_cut/cut_integrationcell.H"
#include "../drt_cut/cut_volumecell.H"

#ifdef PARALLEL
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

#ifdef DEBUG
#include "combust3.H"
#endif


// this flag turns on the former post-evaluation of the cut algorithm implemented by Florian Henke
// due to the inclusion of ghost-penalty terms in the xfem formulation, this is not necessary anymore
// for the time being, Forian's framework is kept
// however, if problems are not observed with the new framework, it will be removed
//#define TOL_ADAPTION_FOR_CUT
//#define VOLTOL 1e-6
// finally, switched to new cut treatment on 25.2.2014 !!!


/*------------------------------------------------------------------------------------------------*
 | constructor                                                                        henke 10/08 |
 *------------------------------------------------------------------------------------------------*/
COMBUST::FlameFront::FlameFront(
    const Teuchos::RCP<DRT::Discretization> fluiddis,
    const Teuchos::RCP<DRT::Discretization> gfuncdis
) :
fluiddis_(fluiddis),
gfuncdis_(gfuncdis),
phin_(Teuchos::null),
phinp_(Teuchos::null),
gradphi_(Teuchos::null),
gradphi2_(Teuchos::null),
curvature_(Teuchos::null),
refinement_(false),
maxRefinementLevel_(0),
xfeminttype_(INPAR::COMBUST::xfemintegration_cut)
{
  // construct interfacehandles using initial flame front
  interfacehandle_ = Teuchos::rcp(new COMBUST::InterfaceHandleCombust(fluiddis,gfuncdis));
  interfacehandle_old_ = Teuchos::rcp(new COMBUST::InterfaceHandleCombust(fluiddis,gfuncdis));
}


/*------------------------------------------------------------------------------------------------*
 | destructor                                                                         henke 10/08 |
 *------------------------------------------------------------------------------------------------*/
COMBUST::FlameFront::~FlameFront()
{
  return;
}

/*------------------------------------------------------------------------------------------------*
 | public: update the flame front                                                     henke 06/10 |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::FlameFront::UpdateOldInterfaceHandle()
{
  if (interfacehandle_==Teuchos::null)
    dserror ("new interfacehandle has to be set before old interfacehandle is updated!");

  std::map<int, GEO::DomainIntCells> domainIntCells = interfacehandle_->DomainIntCells();
  std::map<int, GEO::BoundaryIntCells> boundaryIntCells = interfacehandle_->BoundaryIntCells();
  std::map<int, COMBUST::InterfaceHandleCombust::CutStatus> cutstatus = interfacehandle_->CutState();

  interfacehandle_old_->UpdateInterfaceHandle(
      domainIntCells,
      boundaryIntCells,
      cutstatus);
}


/*------------------------------------------------------------------------------------------------*
 | public: update the flame front                                                     henke 06/10 |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::FlameFront::UpdateFlameFront(
    const Teuchos::ParameterList& combustdyn,
    const Teuchos::RCP<const Epetra_Vector>& phin,
    const Teuchos::RCP<const Epetra_Vector>& phinp,
    const Teuchos::RCP<const Epetra_Vector> phifld
)
{
  // get
  refinement_ = DRT::INPUT::IntegralValue<int>(combustdyn.sublist("COMBUSTION GFUNCTION"),"REFINEMENT");
  maxRefinementLevel_ = Teuchos::getIntParameter(combustdyn.sublist("COMBUSTION GFUNCTION"),"REFINEMENTLEVEL");
  // get type of cut algorithm
  xfeminttype_ = DRT::INPUT::IntegralValue<INPAR::COMBUST::XFEMIntegration>(combustdyn.sublist("COMBUSTION FLUID"),"XFEMINTEGRATION");

  // rearrange and store phi vectors
  StorePhiVectors(phin, phinp);
  if (phifld!=Teuchos::null)
    SetPhinp(phifld);

  // generate the interface geometry based on the G-function (level set field)
  // remark: must be called after StorePhiVectors, since it relies on phinp_
  ProcessFlameFront(phinp_);

  // compute smoothed gradient of G-function field
  // remark: must be called after ProcessFlameFront, since it relies on the new interface position

  const bool smoothed_boundary_integration = DRT::INPUT::IntegralValue<int>
  (combustdyn.sublist("COMBUSTION FLUID"),"SMOOTHED_BOUNDARY_INTEGRATION");
  const INPAR::COMBUST::SmoothGradPhi smoothgradphi = DRT::INPUT::IntegralValue<INPAR::COMBUST::SmoothGradPhi>
  (combustdyn.sublist("COMBUSTION FLUID"),"SMOOTHGRADPHI");
  if( smoothed_boundary_integration and smoothgradphi!=INPAR::COMBUST::smooth_grad_phi_none )
  {
    CallSmoothGradPhi(combustdyn);

    // there are different requirements, depending on what the curvature is needed for
    switch (DRT::INPUT::IntegralValue<INPAR::COMBUST::CombustionType>(combustdyn.sublist("COMBUSTION FLUID"),"COMBUSTTYPE"))
    {
    case INPAR::COMBUST::combusttype_premixedcombustion:
      ComputeCurvatureForCombustion(combustdyn);
      break;
    case INPAR::COMBUST::combusttype_twophaseflow_surf:
    case INPAR::COMBUST::combusttype_twophaseflowjump:
      ComputeCurvatureForSurfaceTension(combustdyn);
      break;
    default:
      break;
    }
  }

  return;
}

/*------------------------------------------------------------------------------------------------*
 | public: generate the interface geometry based on the G-function (level set field)  henke 10/08 |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::FlameFront::ProcessFlameFront(const Teuchos::RCP<const Epetra_Vector> phi)
{
  /* This function is accessible from outside the FlameFront class. It can be called e.g. by the
   * interface handle constructor. If the FlameFront class will always do the same thing, this
   * function could be spared. Then the content of this function could be added to the constructor
   * of this class.
   *
   * henke 10/08
   */

  if (fluiddis_->Comm().MyPID()==0)
    IO::cout << "---  capturing flame front... " << IO::flush;

  const Teuchos::RCP<Epetra_Vector> phicol = Teuchos::rcp(new Epetra_Vector(*fluiddis_->NodeColMap()));
  if (phicol->MyLength() != phi->MyLength())
    dserror("vector phi needs to be distributed according to fluid node column map");

  // maps are defined
  std::map<int,GEO::DomainIntCells >        elementalDomainIntCells;
  std::map<int,GEO::BoundaryIntCells >      elementalBoundaryIntCells;
  std::map<int,COMBUST::InterfaceHandleCombust::CutStatus > elementcutstatus;

  // loop over fluid (combustion) column elements
  // remark: loop over row elements would be sufficient, but enrichment is done in column loop
  for (int iele=0; iele<fluiddis_->NumMyColElements(); ++iele)
  {
    // get element from discretization
    const DRT::Element *ele = fluiddis_->lColElement(iele);

#ifdef DEBUG
    if(ele->ElementType() != DRT::ELEMENTS::Combust3Type::Instance())
      // this is not compulsory, but combust3 elements are expected here!
      dserror("unexpected element type: this should be of combust3 type!");
#endif

    // create refinement cell from a fluid element -> cell will have same geometry as element!
    const Teuchos::RCP<COMBUST::RefinementCell> rootcell = Teuchos::rcp(new COMBUST::RefinementCell(ele));

    // refinement strategy is turned on
    if (refinement_)
    {
      //IO::cout << "starting refinement ..." << IO::endl;
      //IO::cout << "maximal refinement level " << maxRefinementLevel_ << IO::endl;
      if (maxRefinementLevel_ < 0)
        dserror("maximal refinement level not defined");

      RefineFlameFront(rootcell,phi);
    }
    else // refinement strategy is turned off
      FindFlameFront(rootcell,phi);

    /* jetzt habe ich für jedes Element eine "rootcell", an der entweder (refinement on) ein ganzer
     * Baum von Zellen mit Interfaceinformationen hängt, oder (refinement off) eine einzige Zelle
     * hängt. Aus dieser Information muss jetzt das Oberfläche der Flammenfront berechnet werden.
     */

    // generate flame front (interface) geometry
    CaptureFlameFront(rootcell,
        elementalDomainIntCells,
        elementalBoundaryIntCells,
        elementcutstatus);

    // should not be necessary
#if 0
    // delete all refinement cells of root cell
    if (refinement_ == true)
    {
      rootcell->Clear();
    }
#endif
  }

  interfacehandle_->UpdateInterfaceHandle(
      elementalDomainIntCells,
      elementalBoundaryIntCells,
      elementcutstatus);

  if (fluiddis_->Comm().MyPID()==0)
    IO::cout << "done" << IO::endl;
}

/*------------------------------------------------------------------------------------------------*
 | rearrange and store phi vectors for use in the fluid time integration routine      henke 06/10 |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::FlameFront::StorePhiVectors(
    //const Teuchos::RCP<Epetra_Vector> phinm,
    const Teuchos::RCP<const Epetra_Vector> phin,
    const Teuchos::RCP<const Epetra_Vector> phinp)
{
  /* In the processing of the flame front the fluid discretization is cut by the level set
   * function (G-function). This is the reason why fluid elements need to be able to access the
   * corresponding scalar G-function values for all their nodes. Therefore, the ScaTra dof-based
   * vector phinp has to be rearranged in a parallel environment to represent a Fluid node-based
   * vector. This involves two steps:
   * 1. ScaTra DofRowMap -> Fluid NodeRowMap
   * 2. Fluid NodeRowMap -> Fluid NodeColMap
   *
   * henke 02/09
   */

  // reset vectors
  phin_  = Teuchos::null;
  phinp_ = Teuchos::null;

  //------------------------------------------------------------------------------------------------
  // Rearranging phi vectors from gfuncdis DofRowMap to fluiddis NodeRowMap
  //------------------------------------------------------------------------------------------------
  //const Teuchos::RCP<Epetra_Vector> phinmrow = Teuchos::rcp(new Epetra_Vector(*fluiddis_->NodeRowMap()));
  const Teuchos::RCP<Epetra_Vector> phinrow  = Teuchos::rcp(new Epetra_Vector(*fluiddis_->NodeRowMap()));
  const Teuchos::RCP<Epetra_Vector> phinprow = Teuchos::rcp(new Epetra_Vector(*fluiddis_->NodeRowMap()));

  //----------------------------------------------------------------------------------------------
  // congruent (matching) discretizations (Fluid == G-function)
  //----------------------------------------------------------------------------------------------
  if(true)
  {
    /* The premise of geometrically identical discretizations leads to a significant simplification
     * of the rearrangment process. If nodes of both discretizations are distributed in the same way
     * (which they are) and the dofs belonging to nodes lie on the same processors, then the
     * following holds:
     * G-function NodeMap   identical to   Fluid NodeMap
     * G-function NodeMap     similar to   G-function DofMap   and therefore
     * G-function DofMap      similar to   Fluid NodeMap
     *
     * This means that vector phinp living on the G-function DofRowMap can be directly copied to the
     * Fluid NodeRowMap without any rearrangement for congruent discretizations.
     *
     * Here is what we used to do:
     *
     *  // get the G-function values corresponding to fluid nodes
     *  *phinmrow = *phinm;
     *  *phinrow  = *phin;
     *  *phinprow = *phinp;
     *  // ausgeschrieben, was *phinprow = *phinp in einer Zeile macht:
     *  for(int lnodeid=0;lnodeid<fluiddis_->NumMyRowNodes();lnodeid++)
     *  {
     *    // get the G-function value corresponding to this fluid node
     *    double value = (*phinp)[lnodeid];
     *    phinprow->ReplaceMyValue(lnodeid,value);
     *    }
     *
     * henke 02/09
     *
     * This is not possible anymore if we use periodic boundary conditions. The number of nodes remains
     * the same, but dofs are removed, since master and slave nodes share a dof. Things are not so
     * easy anymore and we have to pick the corresponding G-function dof by looping over all fluid
     * nodes.
     *
     * henke 04/11
     */

    //DRT::UTILS::PrintParallelDistribution(*fluiddis_);
    //DRT::UTILS::PrintParallelDistribution(*gfuncdis_);

    // loop all nodes on the processor
    for(int lfluidnodeid=0;lfluidnodeid<fluiddis_->NumMyRowNodes();lfluidnodeid++)
    {
      // get the processor's scatra node
      // remark: we rely on identical parallel distributions of fluid and G-function discretizations, here
      DRT::Node* gfuncnode = gfuncdis_->lRowNode(lfluidnodeid);
#ifdef DEBUG
      if(gfuncdis_->NumDof(0,gfuncnode)!=1) dserror("G-function node should have 1 dof");
#endif
      // find out the local dof id of the dof at the scatra node
      const int gfuncdofgid = gfuncdis_->Dof(0,gfuncnode,0);

      const int lgfuncdofidn = phin->Map().LID(gfuncdofgid);
      const int lgfuncdofidnp = phinp->Map().LID(gfuncdofgid);

      if (lgfuncdofidn < 0) dserror("local dof id not found in map for given global dof id");
      if (lgfuncdofidnp < 0) dserror("local dof id not found in map for given global dof id");

      // now copy the values
      const double valuen = (*phin)[lgfuncdofidn];
      const int errn = phinrow->ReplaceMyValue(lfluidnodeid,0,valuen);
      if (errn) dserror("error while inserting value into phinrow");

      const double valuenp = (*phinp)[lgfuncdofidnp];
      const int errnp = phinprow->ReplaceMyValue(lfluidnodeid,0,valuenp);
      if (errnp) dserror("error while inserting value into phinprow");
    }

  }
  //----------------------------------------------------------------------------------------------
  // no congruent (matching) discretizations (Fluid != G-function)
  //----------------------------------------------------------------------------------------------
  else
  {
    /* In general this operation is very complex. It involves parallel communication and
     * interpolation of the G-function at fluid nodes, if non-matching discretizations for fluid
     * and G-function are used.
     *
     * henke 02/09
     */
    dserror("non-congruent discretizations cannot be handled yet!");
  }

  //------------------------------------------------------------------------------------------------
  // Export vector phinp from fluiddis NodeRowMap to fluiddis NodeColMap for parallel
  // accessibility.
  // remark: SetState() can not be used here, because it is designed for dof-based vectors only.
  //------------------------------------------------------------------------------------------------
  const Teuchos::RCP<Epetra_Vector> phincol = Teuchos::rcp(new Epetra_Vector(*fluiddis_->NodeColMap()));
  LINALG::Export(*phinrow,*phincol);
  const Teuchos::RCP<Epetra_Vector> phinpcol = Teuchos::rcp(new Epetra_Vector(*fluiddis_->NodeColMap()));
  LINALG::Export(*phinprow,*phinpcol);

  // store vector on fluiddis NodeColMap holding G-function values in member variable
  phin_  = phincol;
  phinp_ = phinpcol;
}


/*------------------------------------------------------------------------------------------------*
 | refine the region around the flame front                                           henke 10/08 |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::FlameFront::RefineFlameFront(const Teuchos::RCP<COMBUST::RefinementCell> cell,
    const Teuchos::RCP<const Epetra_Vector> phi)
{
  // get G-Function values of refinement cell
  FindFlameFront(cell,phi);
  // is maximal refinement level already reached?
  if (cell->RefinementLevel() < maxRefinementLevel_) //->No
  {
    if (cell->Bisected()) // cell is bisected ....
    {
      cell->RefineCell(); // ... and will be refined
      //        IO::cout << "RefineCell done" << IO::endl;
      // loop over all new refinement cells and call of RefineFlameFront
      for (int icell = 0; icell < cell->NumOfChildren(); icell++)
        RefineFlameFront(cell->GetRefinementCell(icell),phi);
    }
  }

  return;
}


/*------------------------------------------------------------------------------------------------*
 | interpolate G-function field for a refinement cell                                 henke 10/08 |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::FlameFront::FindFlameFront(
    const Teuchos::RCP<COMBUST::RefinementCell> cell,
    const Teuchos::RCP<const Epetra_Vector> phi)
{
  // get the element which this cell belongs to
  const DRT::Element* ele = cell->Ele();

  //------------------------------------------------------------------------------------------------
  // congruent (matching) discretizations (Fluid == G-function)
  //------------------------------------------------------------------------------------------------
  if (true)
  {
    //----------------------------------------------------------------------------------------------
    // (refinement level == 0) also means: refinement turned off (fluid element == cell!)
    //----------------------------------------------------------------------------------------------
    if (cell->RefinementLevel() == 0)
    {
      // remark: vector "lm" is neccessary, because ExtractMyValues() only accepts "vector<int>"
      // arguments, but ele->NodeIds delivers an "int*" argument
      std::vector<int> lm(ele->NumNode());
      // get vector of node GIDs for this element
      const int* nodeids = ele->NodeIds();
      for (unsigned inode=0; inode < lm.size(); inode++)
        lm[inode] = nodeids[inode];

      //--------------------------------------------------------------------------------------------
      // extract G-function values for all nodes belonging to this fluid element ( == cell!)
      //--------------------------------------------------------------------------------------------
      // create vector "mygfuncvalues" holding G-function values for this element
      std::vector<double> mygfuncvalues(ele->NumNode(),1000.0);
      // get entries in "gfuncvalues" corresponding to node GIDs "lm" and store them in "mygfuncvalues"
      DRT::UTILS::ExtractMyValues(*phi,mygfuncvalues,lm);
      //TEST
      //if (ele->Id()==0)
      //{
      //  IO::cout<< "Gfunc " << ele->Id() << IO::endl;
      //  for(std::size_t ig=0; ig<mygfuncvalues.size(); ig++)
      //    IO::cout << mygfuncvalues[ig] << IO::endl;
      //}

#ifdef DEBUG
      if (lm.size() != mygfuncvalues.size() )
        dserror("unexpected number of nodes: there should be 8 nodes per element!");
#endif

      //--------------------------------------------------------------------------------------------
      // store G-function values in refinement cell and decide whether this cell is intersected
      //--------------------------------------------------------------------------------------------
      // node numbering in element has to match data structure in cell
      cell->SetGfuncValues(mygfuncvalues);
    }
    //----------------------------------------------------------------------------------------------
    // higher level of refinement
    //----------------------------------------------------------------------------------------------
    else
    {
      // get the vertex coordinates of the refinement cell
      // given in the coordinate system of the element == rootcell
      const std::vector<std::vector<double> > vertexcoord = cell->GetVertexCoord();
      const std::size_t numnodecell = vertexcoord.size();
      // to be filled with the G-Function values of the refinement cell
      std::vector<double> gfunctionvalues(numnodecell);
      // G-Function values of the corresponding element
      if (cell->ReturnRootCell()==NULL)
        dserror("return NULL");
      const std::vector<double> gfuncelement = cell->ReturnRootCell()->GetGfuncValues();
      //TEST
      //IO::cout<< "G-Funct of rootcell "<< IO::endl;
      //for (std::size_t i=0; i<gfuncelement.size(); i++)
      //  IO::cout<< gfuncelement[i] << IO::endl;

      const size_t numnode = cell->Ele()->NumNode();
#ifdef DEBUG
      if (gfuncelement.size() != numnode or numnodecell != numnode)
        dserror("number of G-Function values does not match number of element nodes");
#endif
      // loop over all vertices of the refinement cell
      for (std::size_t ivertex = 0; ivertex < numnodecell; ivertex++)
      {
        // evaluate shape function at the coordinates of the current refinement cell vertex
        Epetra_SerialDenseVector funct(numnode);
        DRT::UTILS::shape_function_3D(funct,vertexcoord[ivertex][0],vertexcoord[ivertex][1],vertexcoord[ivertex][2],DRT::Element::hex8);
        // compute G-Function value
        gfunctionvalues[ivertex] = 0.0;
        for (size_t inode = 0; inode < numnode; inode++)
          gfunctionvalues[ivertex] += gfuncelement[inode] * funct(inode);
      }
      cell->SetGfuncValues(gfunctionvalues);
    }
  }
  //------------------------------------------------------------------------------------------------
  // no congruent (matching) discretizations (Fluid != G-function)
  //------------------------------------------------------------------------------------------------
  else
  {
    dserror("non-congruent discretizations cannot be handled yet!");
    if (true) // (refinement level == 0) also means: refinement turned off (fluid element == cell!)
    {
      dserror("non-congruent discretizations not yet working!");
    }
    else // (refinement level != 0) higher level of refinement
    {
      dserror("refinement for non-congruent discretizations not yet working!");
    }
  }

  //------------------------------------------------------------------------------------------------
  // find intersection points of G-function (zero level set) with refinement cell edges
  //------------------------------------------------------------------------------------------------
  if(cell->Bisected() == true)
  {
    FindIntersectionPoints(cell);
  }
  else
  {
    // done (-> next fluid element)
  }

  return;
}


/*------------------------------------------------------------------------------------------------*
 | find intersection points of G-function (level set zero iso-surface) with refinement cell edges |
 |                                                                                rasthofer 06/09 |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::FlameFront::FindIntersectionPoints(const Teuchos::RCP<COMBUST::RefinementCell> cell)
{
  // get G-function values at vertices from refinement cell
  const std::vector<double>& gfuncvalues = cell->GetGfuncValues();
  // get vertex coordinates (local fluid element coordinates) from refinement cell
  const std::vector<std::vector<double> >& vertexcoord = cell->GetVertexCoord();
  // temporary variable to store intersection points
  std::multimap<int,std::vector<double> > intersectionpoints;
  // for double intersected lines of quadratic element a multimap is necessary
  // std::multimap<int,std::vector<double> > intersectionpoints;

  //-------------------------------------------------
  // get vector of lines with corresponding vertices
  //-------------------------------------------------
  // lines: edge numbers and corresponding vertices (hex8)
  std::vector<std::vector<int> > lines;

  switch(cell->Ele()->Shape())
  {
  case DRT::Element::hex8:
  case DRT::Element::hex20:
  {
    lines = DRT::UTILS::getEleNodeNumberingLines(DRT::Element::hex8);
    // remark: vertices are assumed to be numbered in the same way the nodes are
    //         convention documented in globalreport.pdf
    /* L1:  0 1
     * L2:  1 2
     * L3:  2 3
     * L4:  0 3
     * L5:  0 4
     * L6:  1 5
     * L7:  2 6
     * L8:  3 7
     * L9:  4 5
     * L10: 5 6
     * L11: 6 7
     * L12: 4 7
     */

    //-------------------------------
    // determine intersection points
    //-------------------------------
    // loop edges of refinement cell
    for(std::size_t iline=0; iline<lines.size(); iline++)
    {
      // get G-function value of the two vertices defining an edge
      double gfuncval1 = gfuncvalues[lines[iline][0]];
      double gfuncval2 = gfuncvalues[lines[iline][1]];

      std::vector<double> coordinates(3);

      // check for change of sign along edge
      if (gfuncval1*gfuncval2 < 0.0)
      {
        for (int dim = 0; dim < 3; dim++)
        {
          // vertices have one coordinate in common
          if(vertexcoord[lines[iline][0]][dim] == vertexcoord[lines[iline][1]][dim])
          {
            // intersection point has the same coordinate for that direction
            coordinates[dim] = vertexcoord[lines[iline][0]][dim];
          }
          else // compute intersection point
          {
            // linear interpolation (for hex8)
            // x = x1 + (phi(=0) - phi1)/(phi2 - phi1)*(x2 - x1)
            // store intersection point coordinate (local element coordinates) for every component dim
            coordinates[dim] = vertexcoord[lines[iline][0]][dim] - gfuncval1 / (gfuncval2 - gfuncval1)
                * (vertexcoord[lines[iline][1]][dim] - vertexcoord[lines[iline][0]][dim]);

            // shift intersection point to vertex if it is very close
            //if((fabs(vertexcoord[lines[iline][0]][dim]-coordinates[dim]) < 1.0E-4) or
            //  (fabs(vertexcoord[lines[iline][1]][dim]-coordinates[dim]) < 1.0E-4))
            //{
            //  if(fabs(vertexcoord[lines[iline][0]][dim]-coordinates[dim]) < 1.0E-4)
            //  {
            //    IO::cout << "coordinates shifted to vertex 1: " << coordinates[dim] << IO::endl;
            //    coordinates[dim] = vertexcoord[lines[iline][0]][dim];
            //  }
            //  else if(fabs(vertexcoord[lines[iline][1]][dim]-coordinates[dim]) < 1.0E-4)
            //  {
            //    IO::cout << "coordinates shifted to vertex 2: " << coordinates[dim] << IO::endl;
            //    coordinates[dim] = vertexcoord[lines[iline][1]][dim];
            //  }
            //  else
            //  {
            //    dserror("impossible");
            //  }
            //}
          }
        }

        // store coordinates of intersection point for each line
        intersectionpoints.insert(std::pair<int, std::vector<double> >( iline, coordinates));
        //intersectionpoints[iline] = coordinates;
      }
      else
      {
        //do nothing and go to the next line
      }
    }

    break;
  }
  default:
  {
    dserror("FindIntersectionPoints() does not support this element shape!");
    break;
  }
  }

  //TEST
  //for (std::map<int,std::vector<double> >::const_iterator iter = intersectionpoints.begin(); iter != intersectionpoints.end(); ++iter)
  //{
  //  IO::cout<< iter->first << IO::endl;
  //  std::vector<double> coord = iter->second;
  //  for (std::size_t isd=0; isd<3; isd++)
  //  {
  //    IO::cout<< coord[isd] << IO::endl;
  //  }
  //}

  // store intersection points in refinement cell
  cell->intersectionpoints_ = intersectionpoints;

  // possible Gmsh output to visualize the rootcell and its intersectionpoints
  // cell->RootCellToGmsh();

  return;
}


/*------------------------------------------------------------------------------------------------*
 | capture flame front within one root refinement cell (==element) to create domain and boundary  |
 | integration cells                                                                  henke 08/10 |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::FlameFront::CaptureFlameFront(
    const Teuchos::RCP<const COMBUST::RefinementCell>           rootcell,
    std::map<int,GEO::DomainIntCells >&                         elementDomainIntCells,
    std::map<int,GEO::BoundaryIntCells >&                       elementBoundaryIntCells,
    std::map<int,COMBUST::InterfaceHandleCombust::CutStatus >&  elementcutstatus)
{
  // get all refinement cells attached to this root cell (leaves in tree structure)
  // remark: if refinement is turned off, the root cell (==element) itself is returned
  std::vector<const COMBUST::RefinementCell* > RefinementCells;
  rootcell->SearchRefinementCells(RefinementCells);

  // list of domain integration cells
  GEO::DomainIntCells listDomainIntCellsperEle;
  // list of boundary integration cells
  GEO::BoundaryIntCells listBoundaryIntCellsperEle;

  // set element cut status to default value 'undefined'
  COMBUST::InterfaceHandleCombust::CutStatus cutstat = COMBUST::InterfaceHandleCombust::undefined;

  // get vector of G-function values of root cell
  const std::vector<double>& gfuncvaluesrootcell = rootcell->GetGfuncValues();

  // different (spatial) XFEM integration methodologies
  switch (xfeminttype_)
  {
  case INPAR::COMBUST::xfemintegration_cut:
  {
    // loop over all refinement cells attached to root cell
    // remark: this is a partition of the root cell
    for (std::size_t irefcell = 0; irefcell < RefinementCells.size(); irefcell++)
    {
      //-----------------------
      // use GEO::Cut algorithm
      //-----------------------
      const DRT::Element::DiscretizationType distype_cell = RefinementCells[irefcell]->Shape();
      // check if distype of cell correct (hex8)
      if(distype_cell != DRT::Element::hex8) dserror("hex8 refinement cell expected");

      int numnode_cell = DRT::UTILS::getNumberOfElementNodes(distype_cell);
      std::vector<int> nids;
      nids.reserve(numnode_cell);
      for ( int i=0; i<numnode_cell; ++i )
      {
        nids.push_back(i);
      }
      // get G-function values at vertices from cell
      const std::vector<double>& gfuncvalues = RefinementCells[irefcell]->GetGfuncValues();

      //------------------------
      // call GEO::Cut algorithm
      //------------------------
      GEO::CUT::LevelSetIntersection levelset;

      const std::vector<std::vector<double> >& vertexcoord = RefinementCells[irefcell]->GetVertexCoord();

      LINALG::SerialDenseMatrix cellcoord(3,numnode_cell);
      for (int ivert=0; ivert<numnode_cell; ivert++)
      {
        std::copy(vertexcoord[ivert].begin(), vertexcoord[ivert].end(), &cellcoord(0,ivert));
      }

      // check if this refinement cell is cut, according to its G-function values -> add it to 'levelset'
      levelset.AddElement( 1, nids, cellcoord, distype_cell, &gfuncvalues[0], false );

      try
      {
        levelset.Cut();
      }
      catch ( std::runtime_error & err )
      {
        std::cerr << "failed to cut element\n"
            << "coordinates:\n";
        cellcoord.Print( std::cerr );
        std::cerr << "g-function values:\n"
                  << std::setprecision( 16 );
        std::copy( gfuncvalues.begin(), gfuncvalues.end(), std::ostream_iterator<double>( std::cerr, ", " ) );
        std::cerr << "\n";
        throw;
      }

      //-----------------
      // process cut data
      //-----------------
      GEO::CUT::ElementHandle * ehandle = levelset.GetElement( 1 );

      // cell is in contact with the interface (cut or touched)
      if (ehandle!=NULL)
      {
        // get global coordinates of element this cell belongs to
        const DRT::Element* xfemele = rootcell->Ele();
        const int numnode_ele = xfemele->NumNode();

        LINALG::SerialDenseMatrix xyze(3,numnode_ele);
        for(int inode=0;inode<numnode_ele;inode++)
        {
          const DRT::Node* node = xfemele->Nodes()[inode];
          xyze(0,inode) = node->X()[0];
          xyze(1,inode) = node->X()[1];
          xyze(2,inode) = node->X()[2];
        }

        GEO::CUT::plain_element_set cuteles;
        ehandle->CollectElements( cuteles );

        switch(distype_cell)
        {
        case DRT::Element::hex8:
        {
          if (cuteles.size() != 1)
            dserror("one cut element expected for linear elements");
          break;
        }
        case DRT::Element::hex20:
        case DRT::Element::hex27:
        {
          if (cuteles.size() != 8)
            dserror("eight cut elements expected for quadratic elements");
          break;
        }
        default:
          dserror("distype unknown for level set cut algorithm");
        }

        //--------------------------------------------
        // compute volume of refinement cell (element)
        //--------------------------------------------

        for ( GEO::CUT::plain_element_set::const_iterator icutele=cuteles.begin(); icutele!=cuteles.end(); ++icutele )
        {
          // get pointer to cut element
          GEO::CUT::Element* cutele = *icutele;

          GEO::CUT::plain_volumecell_set volcells;
          volcells = cutele->VolumeCells();
          const int numvolcells = volcells.size();

          //---------------------------------------------------------------
          // SOME IMPORTANT REMARKS ON THE CUT SITUATION IN THE CUT !!!!!
          //---------------------------------------------------------------

          // In the cut algorithm, we are dealing with several tolerances. Two tolerances are particularly relevant for
          // the present level-set intersection:
          //
          // - #define VOLUMETOL   1e-14  // the tolerance when created integration cells are discarded, the boundary cell however will be created
          //   If the intersected element exhibits a volume cell which does not represent a simple geometry, i.e. a hexahedra, a tetrahedra, a wedge,...,
          //   and cannot be tetrahedralize using a simple a priori known scheme (cut library), TetMesh is called to tetrahedralize the volume cell. This
          //   tolerance will prevent TetMesh to tetrahedralize volumes smaller than the given value. Boundary cells obtained from the
          //   adjacent volume cell are kept!
          // - #define TOLERANCE   1e-12
          //   Based on this tolerance for nodal values, the cut algorithm determines whether to start cutting the element.
          //
          // This has the following impact on critical intersection situations:
          //
          // Case 1: a touched element
          //                   |
          //         ----------|--
          //         |         | |
          //         |         | |
          //         |         | |
          //         |         | |
          //         |         | |
          //         ----------|--
          //                   |
          //              straight interface approaching right boundary
          //
          //         element is bisected until interface and, hence, level-set values at the right nodes are within the tolerance
          //         -> two volume cells, one face boundary
          //         then, when the tolerance is reached, the cut returns a touched element
          //         -> one volume cell, one face boundary
          //
          // Case 2: from skew interface to touched element
          //                    |
          //         ----------|--
          //         |        |  |
          //         |       |   |
          //         |      |    |
          //         |     |     |
          //         |    |      |
          //         ----|--------
          //            |
          //         skew interface approaching right boundary at the lower corner and being zero at the upper right corner
          //
          //         element is bisected until interface and, hence, level-set value at the lower right node is within the tolerance
          //         -> two volume cells, one face boundary
          //         then, when the tolerance is reached, the cut returns a touched element
          //         -> one volume cell, one face boundary
          //
          // Note: This smooth transition from bisected to touched is possible, since the related volume cells represent geometries with
          //       simple tetrahedralization. A call of TetMesh is not required.
          //
          // Case 3: small tetrahedra (one edge of a three dimensional domain is cut off) and wedge cuts
          //
          //         -------------
          //         |           |
          //         |           |
          //         |           |
          //         |           |
          //         |          /|
          //         -------------
          //                   interface cuts of small tetrahedra or wedge
          //
          //         For this case, two volume cells are first created. Since the left part of the element does not constitute a simple
          //         geometry anymore, TetMesh is called to mesh the entire geometry (we need matching boundary cells at the interface).
          //         If the right element part is too small, i.e., below the aforementioned volume, TetMesh will only mesh the large volume
          //         cell. The boundary cells at the interface are stored for both volume cells however.
          //
          // Note: While we can directly take the suggestion of the cut algorithm for case 1 and 2, special care for this case is required.

          //-----------------
          // cell is not intersected: CASE 1 and 2
          //-----------------
          if (numvolcells==1)
          {
            //IO::cout << "element " << rootcell->Ele()->Id() << " is uncut " << IO::endl;
            // store domain integration cells
            const size_t numstoredvol = StoreDomainIntegrationCells(cutele,listDomainIntCellsperEle,xyze);
            if (numstoredvol != 1) dserror("No volume stored for touched element!");

            // check for potential boundary cells: corresponds to touched elements
            const GEO::CUT::plain_boundarycell_set bcells = (*volcells.begin())->BoundaryCells();
            if (bcells.size() > 0)
            {
              //IO::cout << "element " << rootcell->Ele()->Id() << " is potentially touched" << IO::endl;
              //element is touched and we have to store the corresponding boundary cell
              const size_t numstoredbound = StoreBoundaryIntegrationCells(cutele,listBoundaryIntCellsperEle,xyze);
              if (numstoredbound == 1)
              {
                //IO::cout << "element " << rootcell->Ele()->Id() << " is touched" << IO::endl;
                cutstat = COMBUST::InterfaceHandleCombust::touched;
              }
              else if (numstoredbound == 0)
                cutstat = COMBUST::InterfaceHandleCombust::uncut;
              else
               dserror("Unexpected number of boundaries!");
            }
            else
              cutstat = COMBUST::InterfaceHandleCombust::uncut;
          }
          //-----------------
          // cell is bisected
          //-----------------
          else if (numvolcells==2)
          {
            //IO::cout << "element " << rootcell->Ele()->Id() << " is really bisected" << IO::endl;
            // store domain integration cells
            const size_t numstoredvol = StoreDomainIntegrationCells(cutele,listDomainIntCellsperEle,xyze);
            // store boundary integration cells
            const size_t numstoredbound =StoreBoundaryIntegrationCells(cutele,listBoundaryIntCellsperEle,xyze);
            if (numstoredbound != 1) dserror("One interface boundary expected for bisected element!");

//            IO::cout << "storedvol " << numstoredvol << IO::endl;
//            IO::cout << "storedbound " << numstoredbound << IO::endl;

            // all volume cells and boundary cells have been stored (this is the usual case)
            if (numstoredvol==2)
            {

              cutstat = COMBUST::InterfaceHandleCombust::bisected;
              //std::cout << "element " << rootcell->Ele()->Id() << " is bisected" << IO::endl;
            }
            // element is intersected but one volume cell is small (boundary is likewise small) and falls below
            // the tolerance for TetMesh such that it does not get any domain integration cells
            // in this case, the element is treated as uncut
            else if (numstoredvol==1)
            {
              // the dofmap creator will consider nodes of touched elements as enriched if phi_node < (0.1*hk_eleDiam)
              // for an element to be treated as touched four nodes have to be below this tolerance
              // this has to be checked here before an element will be denoted as touched

              unsigned nodecount = 0;
              const double hk_eleDiam = COMBUST::getEleDiameter<DRT::Element::hex8>(xyze);
              // enrich only nodes which are touched, that means nodes which have Gfunc ~ 0.0
              for (std::size_t rr = 0; rr<gfuncvalues.size(); ++rr)
              {
                if (std::abs(gfuncvalues[rr]) < 0.05*hk_eleDiam)// this tolerance here must be equal to the one given in xdofmapcreation_combust.cpp
                  nodecount += 1;
              }

              if (nodecount == 4)
                cutstat = COMBUST::InterfaceHandleCombust::touched;
              else if (nodecount < 4)
                  cutstat = COMBUST::InterfaceHandleCombust::uncut;
              else
                dserror("To many nodes with zero level-set for bisected element!");
            }
            // something went wrong
            else
            {
              FlameFrontToGmsh(rootcell->ReturnRootCell(),listBoundaryIntCellsperEle,listDomainIntCellsperEle);
              dserror("flame front for bisected element %d could not be captured",rootcell->Ele()->Id());
            }
          }
          //---------------------------------------
          // cell is trisected
          //---------------------------------------
          else if (numvolcells==3)
          {
            // store domain integration cells
            const size_t numstoredvol = StoreDomainIntegrationCells(cutele,listDomainIntCellsperEle,xyze);
            // store also boundary integration cells
            const size_t numstoredbound = StoreBoundaryIntegrationCells(cutele,listBoundaryIntCellsperEle,xyze);
            if (numstoredbound != 2) dserror("Two interface boundaries expected for trisected element!");

            if (numstoredvol == 3 or numstoredvol == 2)
            {
              // treat element as bisected
              cutstat = COMBUST::InterfaceHandleCombust::bisected;
            }
            // element is trisected, but two volume cells do not contain any domain integration cell
            // thus element is considered uncut
            else if (numstoredvol == 1)
            {
              // update cut status of root cell (element)
              cutstat = COMBUST::InterfaceHandleCombust::uncut;
              // we cannot simply distinguish between not touched, touched and double touched here
              // moreover, everything we are doing here is currently wrong since multiple enrichments are not
              // available in the current implementation
              // therefore, we simply take the easiest way
            }
            else
            {
              FlameFrontToGmsh(rootcell->ReturnRootCell(),listBoundaryIntCellsperEle,listDomainIntCellsperEle);
              dserror("flame front for bisected element %d could not be captured",rootcell->Ele()->Id());
            }
          }
          //---------------------------------------
          // cell is multi-sected
          //---------------------------------------
          else if (numvolcells>3 and numvolcells<6)
          {
            // store domain integration cells
            const size_t numstoredvol = StoreDomainIntegrationCells(cutele,listDomainIntCellsperEle,xyze);
            // store also boundary integration cells
            /* const size_t numstoredbound = */ StoreBoundaryIntegrationCells(cutele,listBoundaryIntCellsperEle,xyze);

            if (numstoredvol > 1)
            {
              // treat element as bisected
              cutstat = COMBUST::InterfaceHandleCombust::bisected;
            }
            // element is multisected, but several volume cells do not contain any domain integration cell
            // thus element is considered uncut
            else if (numstoredvol == 1)
            {
              // update cut status of root cell (element)
              cutstat = COMBUST::InterfaceHandleCombust::uncut;
              // note, everything we are doing here is currently wrong since multiple enrichments are not
              // available in the current implementation
              // therefore, we simply take the easiest way
            }
            else
            {
              FlameFrontToGmsh(rootcell->ReturnRootCell(),listBoundaryIntCellsperEle,listDomainIntCellsperEle);
              dserror("flame front for bisected element %d could not be captured",rootcell->Ele()->Id());
            }
          }
          else dserror("cut status could not be determined");

        }
      }
      //--------------
      // cell is uncut
      //--------------
      else
      {
        // update cut status of root cell (element)
        cutstat = COMBUST::InterfaceHandleCombust::uncut;
        // store refinement cell in list of domain integration cells
        StoreDomainIntegrationCell(RefinementCells[irefcell],listDomainIntCellsperEle);
        // remark: no boundary integration cell required
      }

    } // end loop over all refinement cells
    break;
  }
  //-----------------------------------------------------------------------------------------
  // decompose hexahedral elements into 6 tetrahedra first to create domain integration cells
  //-----------------------------------------------------------------------------------------
  case INPAR::COMBUST::xfemintegration_tetrahedra:
  {
    // loop over all refinement cells attached to root cell
    // remark: this covers the whole domain of the root cell
    for (std::size_t irefcell = 0; irefcell < RefinementCells.size(); irefcell++)
    {
#ifdef DEBUG
      // check if distype of cell correct (hex8)
      const DRT::Element::DiscretizationType distype = RefinementCells[irefcell]->Shape();
      if(distype != DRT::Element::hex8) dserror("hex8 refinement cell expected");
#endif

      if(RefinementCells[irefcell]->Bisected())
      {
        GEO::TetrahedraDecomposition decomposition(RefinementCells[irefcell], listBoundaryIntCellsperEle, listDomainIntCellsperEle);
        // update cut status of root cell (element)
        cutstat = COMBUST::InterfaceHandleCombust::bisected;
      }
      // the cell is touched (interface aligned with a cell surface)
      else if(RefinementCells[irefcell]->Touched())
      {
        // store hex8 domain integration cell
        StoreDomainIntegrationCell(RefinementCells[irefcell],listDomainIntCellsperEle);
        // store quad4 boundary integration cell
        StoreBoundaryIntegrationCell(RefinementCells[irefcell],listBoundaryIntCellsperEle);
        // update cut status of root cell (element)
        cutstat = COMBUST::InterfaceHandleCombust::touched;
      }
      else // (not bisected) and (not touched)
      {
        // store hex8 domain integration cell
        StoreDomainIntegrationCell(RefinementCells[irefcell],listDomainIntCellsperEle);
        // no boundary integration cell required
        // update cut status of root cell (element)
        cutstat = COMBUST::InterfaceHandleCombust::uncut;
      }

    } // end loop over all refinement cells
    break;
  }
  //----------------------------------------------------------------------------
  // use refinement strategy to create small hexahedral domain integration cells
  // which will not be aligned with the interface
  //----------------------------------------------------------------------------
  case INPAR::COMBUST::xfemintegration_hexahedra:
  {
#ifdef DEBUG
    // check if distype of cell correct (hex8)
    const DRT::Element::DiscretizationType distype = rootcell->Shape();
    if(distype != DRT::Element::hex8) dserror("hex8 refinement cell expected");
#endif
    if(rootcell->Bisected())
    {
      // loop over all refinement cells attached to root cell and create hex8 domain integration cells
      // remark: this covers the whole domain of the root cell
      for (std::size_t irefcell = 0; irefcell < RefinementCells.size(); irefcell++)
      {
        // store hex8 domain integration cell
        StoreDomainIntegrationCell(RefinementCells[irefcell],listDomainIntCellsperEle);
      } // end loop over all refinement cells

      // create triangular boundary integration cells for the root cell
      // remark: Creating triangles for every leaf cell would also work, but create many more
      //         boundary integration cells. Building triangles only for the root cell leads to
      //         fewer cells, the size of which is comparable to those obtained for the TetGen
      //         strategy without refinement.
      std::vector<GEO::DomainIntCell> emptyvec;
      // pass empty vector to buildPLC ('listDomainIntCellsperEle' is not needed for this call,
      // since only boundary integration cells are built)
      // - build piecewise linear complex (PLC)
      // - store triangles as boundary integration cells
      //TODO: das ist der richtige Weg, aber dazu darf root cell nicht const uebergeben werden
      //FindIntersectionPoints(rootcell);
      buildPLC(rootcell->ReturnRootCell(),gfuncvaluesrootcell,emptyvec,listBoundaryIntCellsperEle);

      // check, if we really got integration cells on both sides of the interface
      double vol_minus = 0.0;
      double vol_plus = 0.0;
      for (std::size_t intcell = 0; intcell < listDomainIntCellsperEle.size(); intcell++)
      {
        GEO::DomainIntCell iintcell = listDomainIntCellsperEle[intcell];
        if (iintcell.getDomainPlus()==true)
          vol_plus += iintcell.VolumeInXiDomain(*rootcell->Ele());
        else
          vol_minus += iintcell.VolumeInXiDomain(*rootcell->Ele());
      }
      // element is uncut or touched
      if (vol_minus == 0.0 or vol_plus == 0.0)
      {
        //IO::cout << "integration cells only on one side of interface" << IO::endl;
        //IO::cout << "number of integration cells:  " << listDomainIntCellsperEle.size() << IO::endl;
        // delete integration cells
        listDomainIntCellsperEle.clear();
        // and store hex8 domain integration cell instead
        StoreDomainIntegrationCell(rootcell->ReturnRootCell(),listDomainIntCellsperEle);

        // check size of boundary integration cells
        double boundaryarea = 0.0;
        const double areatol = 0.01;
        const double areaele = 4.0;
        for (std::size_t bintcell = 0; bintcell < listBoundaryIntCellsperEle.size(); bintcell++)
        {
          GEO::BoundaryIntCell ibintcell = listBoundaryIntCellsperEle[bintcell];
          dserror("Function AreaInXiDomain() not implemented!");
          //boundaryarea += ibintcell.AreaInXiDomain();
        }
        if (boundaryarea/areaele < areatol)
        {
          // delete boundary integration cells
          listBoundaryIntCellsperEle.clear();
          // adapt cut status
          cutstat = COMBUST::InterfaceHandleCombust::uncut;
        }
        else
        {
          if (listDomainIntCellsperEle[1].getDomainPlus())
          {
            // element is touched
            cutstat = COMBUST::InterfaceHandleCombust::touched;
          }
        }
      }
      // element is bisected
      else
        cutstat = COMBUST::InterfaceHandleCombust::bisected;
    }
    // the rootcell is touched (interface aligned with an element surface)
    // remark: in all other cases (leaf cell is touched) the root cell will be bisected
    //         -> previous case, since domain hexahedra are created for both cases
    else if(rootcell->Touched())
    {
      // store hex8 domain integration cell
      StoreDomainIntegrationCell(rootcell->ReturnRootCell(),listDomainIntCellsperEle);
      // store quad4 boundary integration cell
      StoreBoundaryIntegrationCell(rootcell->ReturnRootCell(),listBoundaryIntCellsperEle);
      // update cutstatus
      cutstat = COMBUST::InterfaceHandleCombust::touched;
    }
    else // (not bisected) and (not touched)
    {
      // store hex8 domain integration cell
      StoreDomainIntegrationCell(rootcell->ReturnRootCell(),listDomainIntCellsperEle);
      // no boundary integration cell required
      // set cutstatus
      cutstat = COMBUST::InterfaceHandleCombust::uncut;
    }
    break;
  }
  default:
  {
    dserror("unknown type of XFEM integration");
    break;
  }
  }

  // - a refined linear element consisting of touched refinement cells is actually bisected, if there
  //   exist refinement cells in both domains (interface exactly aligned with refinement cell boundaries)
  // - a non-refined quadratic element constisting of touched sub-elements is actually bisected, if the
  //   exist refinement cells in both domains (interface exactly aligned with sub-element boundaries)
  if (cutstat == COMBUST::InterfaceHandleCombust::touched)
  {
    if( (refinement_ or rootcell->Ele()->Shape()!=DRT::Element::hex8) )
    {
      // get domain of first domain integration cell
      const bool firstcellplus = listDomainIntCellsperEle[0].getDomainPlus();
      for(GEO::DomainIntCells::const_iterator itercell = listDomainIntCellsperEle.begin(); itercell != listDomainIntCellsperEle.end(); ++itercell )
      {
        // if any domain integration cell belongs to the other domain, the element is bisected
        if (itercell->getDomainPlus() != firstcellplus)
        {
          cutstat = COMBUST::InterfaceHandleCombust::bisected;
          break;
        }
      }
    }
  }

  //------------------------------------------------------
  // store flame front information per element (root cell)
  //------------------------------------------------------
  // store list of integration cells per element in flame front
  const int eleid = rootcell->Ele()->Id();
  elementDomainIntCells[eleid] = listDomainIntCellsperEle;
  // if there exist boundary integration cells
  if(listBoundaryIntCellsperEle.size() > 0)
    elementBoundaryIntCells[eleid] = listBoundaryIntCellsperEle;

  // add cut status for this root cell
  elementcutstatus[eleid] = cutstat;

  //if (cutstat == COMBUST::InterfaceHandleCombusttouched)
  //  FlameFrontToGmsh(rootcell->ReturnRootCell(),listBoundaryIntCellsperEle,listDomainIntCellsperEle);

  return;
}


/*------------------------------------------------------------------------------------------------*
 | identify status of a domain integration cell based on its cut data                 henke 04/11 |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::FlameFront::IdentifyDomainIntegrationCellStatus(
    GEO::DomainIntCell& cell,
    const GEO::CUT::ElementHandle& e)
{
  return;
}


/*------------------------------------------------------------------------------------------------*
 | identify status of a boundary integration cell based on its cut data               henke 04/11 |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::FlameFront::IdentifyBoundaryIntegrationCellStatus(
    GEO::BoundaryIntCell& cell,
    const GEO::CUT::ElementHandle& e)
{
  return;
}


/*------------------------------------------------------------------------------------------------*
 | project midpoint on level set zero iso-surface                                     henke 10/10 |
 *------------------------------------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType DISTYPE,
class V>
bool projectMidpoint(
    V& valuesGcell,                      // values of G-function at vertices of refinement cell
    const std::vector<double>& midpoint, // midpoint vector
    LINALG::Matrix<3,1>&       projpoint // vector of point projected on level set iso-surface
)
{
  // indicator for convergence of Newton-Raphson scheme
  bool converged = false;

  // size of system of equations
  // remark: number space dimensions for 3d combustion problem + auxiliary variable (alpha)
  const size_t nsys = 4;
  // here, a tet4 or hex8 cell is expected
  const size_t numvertices = DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;

  static LINALG::Matrix<numvertices,1> funct(true);  // shape functions
  static LINALG::Matrix<3,numvertices> deriv1(true); // derivatives of shape functions
  static LINALG::Matrix<6,numvertices> deriv2(true); // second derivatives of shape functions
  static LINALG::Matrix<1,1> valueG(true);           // value of G-function
  static LINALG::Matrix<3,1> gradG(true);            // gradient of G-function
  static LINALG::Matrix<6,1> grad2ndG(true);         // second derivatives of G-function

  //----------------------------------
  // start values for iterative scheme
  //----------------------------------
  // start position (projection midpoint = midpoint)
  projpoint(0) = midpoint[0];
  projpoint(1) = midpoint[1];
  projpoint(2) = midpoint[2];

  // auxiliary variable
  // remark: forth unknown to close system of equations; arbitrary value
  double alpha = 0.0;

  // function F (system of equations)
  static LINALG::Matrix<nsys,1> f(true);
  // gradient of function F (dF/deta(0), dF/deta(1), dF/dalpha)
  static LINALG::Matrix<nsys,nsys> gradf(true);
  // increment in Newton iteration (unknown to be solved for)
  static LINALG::Matrix<nsys,1> incr(true);

  // maximum number Newton iterations
  size_t maxiter = 5;
  // convergence tolerance
  double conv = 0.0;

  //------------------------------------------------------
  // Newton-Raphson loop for non-linear projection problem
  //------------------------------------------------------
  for (size_t iter=0;iter<maxiter;++iter)
  {
    // evaluate shape functions in boundary cell space at current position
    funct.Clear();
    DRT::UTILS::shape_function_3D(funct,projpoint(0),projpoint(1),projpoint(2),DISTYPE);
    // evaluate derivatives of shape functions in boundary cell space at current position
    deriv1.Clear();
    DRT::UTILS::shape_function_3D_deriv1(deriv1,projpoint(0),projpoint(1),projpoint(2),DISTYPE);
    deriv2.Clear();
    DRT::UTILS::shape_function_3D_deriv2(deriv2,projpoint(0),projpoint(1),projpoint(2),DISTYPE);

    // evaluate gradient of G-function at current position \xi_1,\xi_2,\xi_3
    // gradG(i,j) = deriv(j,k)*valuesGcell(i,k)
    gradG.Clear();
    gradG.MultiplyNN(deriv1,valuesGcell);

    grad2ndG.Clear();
    grad2ndG.MultiplyNN(deriv2,valuesGcell);

    //---------------------------------------------------
    // build system of equations F and its gradient gradF
    //---------------------------------------------------
    // clear static arrays
    f.Clear();
    gradf.Clear();
    incr.Clear();

    // evaluate function F
    valueG.MultiplyTN(funct,valuesGcell);
    f(0) = valueG(0,0);                                 // G
    f(1) = projpoint(0) - midpoint[0] + gradG(0)*alpha; // xi_1 - midpoint_1 + {d G}/{d xi_1} * alpha
    f(2) = projpoint(1) - midpoint[1] + gradG(1)*alpha; // xi_2 - midpoint_2 + {d G}/{d xi_2} * alpha
    f(3) = projpoint(2) - midpoint[2] + gradG(2)*alpha; // xi_3 - midpoint_3 + {d G}/{d xi_3} * alpha

    // evaluate gradF (Jacobian matrix)
    gradf(0,0) = gradG(0); // {d G}/{dxi_1}
    gradf(0,1) = gradG(1); // {d G}/{dxi_2}
    gradf(0,2) = gradG(2); // {d G}/{dxi_3}
    gradf(0,3) = 0.0;

    gradf(1,0) = 1.0;
    gradf(1,1) = grad2ndG(3)*alpha; // {d^2 G}/{dxi_1*dxi_2} * alpha
    gradf(1,2) = grad2ndG(4)*alpha; // {d^2 G}/{dxi_1*dxi_3} * alpha
    gradf(1,3) = gradG(0);          // {d G}/{dxi_1}

    gradf(2,0) = grad2ndG(3)*alpha; // {d^2 G}/{dxi_2*dxi_1} * alpha
    gradf(2,1) = 1.0;
    gradf(2,2) = grad2ndG(5)*alpha; // {d^2 G}/{dxi_2*dxi_3} * alpha
    gradf(2,3) = gradG(1);          // {d G}/{dxi_2}

    gradf(3,0) = grad2ndG(4)*alpha; // {d^2 G}/{dxi_3*dxi_1} * alpha
    gradf(3,1) = grad2ndG(5)*alpha; // {d^2 G}/{dxi_3*dxi_2} * alpha
    gradf(3,2) = 1.0;
    gradf(3,3) = gradG(2);          // {d G}/{dxi_3}

    // check convergence
    conv = sqrt(f(0)*f(0)+f(1)*f(1)+f(2)*f(2)+f(3)*f(3));
    //std::cout << "iteration " << iter << ": -> |f|=" << conv << IO::endl;

    if (conv <= 1.0E-15) break;
    //----------------------------------------------------
    // solve linear system of equations: gradF * incr = -F
    //----------------------------------------------------
    // F = F*-1.0
    f.Scale(-1.0);
    // solve A.X=B
    LINALG::FixedSizeSerialDenseSolver<nsys,nsys,1> solver;
    solver.SetMatrix(gradf);              // set A=gradF
    solver.SetVectors(incr, f);           // set X=incr, B=F
    solver.FactorWithEquilibration(true); // "some easy type of preconditioning" (Michael)
    int err2 = solver.Factor();           // ?
    int err = solver.Solve();             // incr = gradF^-1.F
    if ((err!=0) || (err2!=0))
      dserror("solving linear system in Newton-Raphson method for projection failed");

    // update projection of point and alpha
    projpoint(0) += incr(0);
    projpoint(1) += incr(1);
    projpoint(2) += incr(2);
    alpha  += incr(3);
  }

  // Newton iteration unconverged
  if (conv > 1.0E-15)
  {
    IO::cout << "Newton-Raphson algorithm for projection of midpoint did not converge" << IO::endl;
  }
  else
  {
    converged = true;
    //LINALG::Matrix<3,1> diff(true);
    //diff(0) = midpoint[0]-projpoint(0);
    //diff(1) = midpoint[1]-projpoint(1);
    //diff(2) = midpoint[2]-projpoint(2);
    //std::cout << diff << IO::endl;
  }
  return converged;
}


/*------------------------------------------------------------------------------------------------*
 | projects midpoint on front and back side of cell for 2D-applications               henke 10/10 |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::FlameFront::projectMidpoint2D(
    std::vector<std::vector<int> >&       trianglelist,
    std::multimap<int,std::vector<int> >& segmentlist,
    std::vector<std::vector<double> >&    pointlist,
    const std::vector<int>&               polypoints,
    const std::vector<double>&            midpoint // midpoint vector
)
{
  if(maxRefinementLevel_>0) dserror("This is not implemented for refinement!");
  if(polypoints.size()-1 != 4) dserror("this is not a 2D problem");

  // determine global third component (pseudo 3D direction)
  // third component of midpoint will be zero in paramenter space (without refinement strategy)
  const double midcomp = 0.0;
  int thirddim = -1;
  for (unsigned idim=0;idim<3;++idim)
  {
    // replace absolute tolerances by relative tolerances
    if ((midpoint[idim] > midcomp-1.0E-8) and (midpoint[idim] < midcomp+1.0E-8))
    {
      thirddim = idim;
      break;
    }
  }
  // create midpoints on front and back side of cell, respectively
  std::vector<double> midpointback(3);
  std::vector<double> midpointfront(3);
  midpointback = midpoint;
  midpointback[thirddim] = -1.0;
  midpointfront = midpoint;
  midpointfront[thirddim] = 1.0;

  //TEST
  //for (std::size_t iter=0; iter<pointlist.size(); ++iter)
  //{
  //  IO::cout<< iter << IO::endl;
  //  std::vector<double> coord = pointlist[iter];
  //  for (std::size_t isd=0; isd<3; isd++)
  //  {
  //    IO::cout<< coord[isd] << IO::endl;
  //  }
  //}

  // find IDs of both intersection points located on the back side
  std::vector<double> point1(3);
  std::vector<double> point2(3);
  size_t ipoint = 777;
  for (ipoint=0;ipoint<4;++ipoint)
  {
    // two consecutive points have to have a "third coordinate" of value -1.0
    point1 = pointlist[polypoints[ipoint]];
    point2 = pointlist[polypoints[ipoint+1]];
    // replace absolute tolerances by relative tolerances
    if ((point1[thirddim] > -1.0-1.0E-8) and (point1[thirddim] < -1.0+1.0E-8) and
        (point2[thirddim] > -1.0-1.0E-8) and (point2[thirddim] < -1.0+1.0E-8))
    {
      // this 'ipoint' is the point
      break;
    }
  }
  // replace absolute tolerances by relative tolerances
  if(!((point1[thirddim] > -1.0-1.0E-8) and (point1[thirddim] < -1.0+1.0E-8) and
       (point2[thirddim] > -1.0-1.0E-8) and (point2[thirddim] < -1.0+1.0E-8)))
    dserror("2D midpoint projection algorithm failed");

  // add frontside and backside mid point to list of points
  size_t midpointback_id = pointlist.size();
  pointlist.push_back(midpointback);
  size_t midpointfront_id = pointlist.size();
  pointlist.push_back(midpointfront);

  // compute IDs of intersection points starting from the first point on the back side
  size_t first_id  = ipoint;       // first point on the back side
  size_t second_id = (ipoint+1)%4; // second point on the front side
  size_t third_id  = (ipoint+2)%4; // first point on the fornt side
  size_t forth_id  = (ipoint+3)%4; // second point on the front side

  // combine 6 points to 4 triangles with correct orientation
  // remark: midpointback is located between the first and second intersection point
  //         midpointfront is located between the third and forth intersection point
  std::vector<int> trianglepoints(3);
  // first triangle
  trianglepoints[0] = polypoints[first_id];
  trianglepoints[1] = midpointfront_id;
  trianglepoints[2] = polypoints[forth_id];
  trianglelist.push_back(trianglepoints);
  // second triangle
  trianglepoints[0] = polypoints[first_id];
  trianglepoints[1] = midpointback_id;
  trianglepoints[2] = midpointfront_id;
  trianglelist.push_back(trianglepoints);
  // third triangle
  trianglepoints[0] = midpointback_id;
  trianglepoints[1] = polypoints[third_id];
  trianglepoints[2] = midpointfront_id;
  trianglelist.push_back(trianglepoints);
  // forth triangle
  trianglepoints[0] = midpointback_id;
  trianglepoints[1] = polypoints[second_id];
  trianglepoints[2] = polypoints[third_id];
  trianglelist.push_back(trianglepoints);

  //TEST
  //IO::cout<<"triangles"<< IO::endl;
  //for (std::size_t itriangle=0; itriangle<trianglelist.size(); itriangle++)
  //{
  //  IO::cout<< "dreieck " << itriangle<< IO::endl;
  //  for (int i=0; i<3; i++)
  //    IO::cout<< trianglelist[itriangle][i]<<IO::endl;
  //}

  // replace segment on back side by two segments connected to midpointback
  size_t segid = 777;
  for (std::map<int,std::vector<int> >::const_iterator iter = segmentlist.begin(); iter != segmentlist.end(); ++iter)
  {
    std::vector<int> points = iter->second;
    if(((points[0]==polypoints[first_id]) or (points[1] == polypoints[first_id])) and
       ((points[0]==polypoints[second_id]) or (points[1] == polypoints[second_id])))
    {
      segid = iter->first;
    }
  }
  segmentlist.erase(segid);
  std::vector<int> segment(2);
  segment[0] = polypoints[first_id];
  segment[1] = midpointback_id;
  segmentlist.insert(std::pair<int,std::vector<int> >(segid,segment));

  segment[0] = midpointback_id;
  segment[1] = polypoints[second_id];
  segmentlist.insert(std::pair<int,std::vector<int> >(segid,segment));

  // replace segment on front side by two segments connected to midpointfront
  for (std::map<int,std::vector<int> >::const_iterator iter = segmentlist.begin(); iter != segmentlist.end(); ++iter)
  {
    std::vector<int> points = iter->second;
    if(((points[0]==polypoints[third_id]) or (points[1] == polypoints[third_id])) and
        ((points[0]==polypoints[forth_id]) or (points[1] == polypoints[forth_id])))
    {
      segid = iter->first;
    }
  }
  segmentlist.erase(segid);
  segment[0] = polypoints[third_id];
  segment[1] = midpointfront_id;
  segmentlist.insert(std::pair<int,std::vector<int> >(segid,segment));
  segment[0] = midpointfront_id;
  segment[1] = polypoints[forth_id];
  segmentlist.insert(std::pair<int,std::vector<int> >(segid,segment));

  // add segment connecting both midpoints (arbitrary ID 6)
  segment[0] = midpointback_id;
  segment[1] = midpointfront_id;
  segmentlist.insert(std::pair<int,std::vector<int> >(6,segment));
}

/*------------------------------------------------------------------------------------------------*
 | triangulate the interface (flame front) inside a refinement cell                   henke 10/08 |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::FlameFront::TriangulateFlameFront(
    const COMBUST::RefinementCell* cell,
    std::vector<std::vector<int> >&       trianglelist,
    std::multimap<int,std::vector<int> >& segmentlist,
    std::vector<std::vector<double> >&    pointlist,
    std::multimap<int,int>&               intersectionpointsids,
    const std::vector<double>&            gfuncvalues
)
{
  // check -> is cell really bisected?
  // This means, not only edges of the cell are contained in segmentlist (key=-1)
  if(segmentlist.count(-1) == segmentlist.size()) dserror("This cell is touched!");

  std::map<int,std::vector<int> > polygonpoints;
  int actpoint;

  // determine start segment for polygon
  // if a surface is intersected twice, the is the start segment
  //  -> so the second segment can not be forgotten
  //  -> and an edge is also not a good beginning
  int numpolygon = 1;
  int beginpolygon=-1;
  for (std::multimap<int,std::vector<int> >::iterator seglistiter = segmentlist.begin(); seglistiter != segmentlist.end(); seglistiter++)
  {
    if(seglistiter->first!=-1)
    {
      if (segmentlist.count(seglistiter->first)==2)
      {
        if (numpolygon==1)
        {
          numpolygon = 2;
          beginpolygon = seglistiter->first;
        }
      }
      else if (segmentlist.count(seglistiter->first)==1)
      {
        if (beginpolygon<0)
        {
          beginpolygon = seglistiter->first;
        }
      }
      else
        dserror("can't build polygon");
    }
  }

  //-----------------------------------------------------------------------
  // build polygon representing the boundary of the interface patch
  //-----------------------------------------------------------------------
  int j=0;
  //loop over all polygons
  for (std::multimap<int,std::vector<int> >::iterator segmentiter=segmentlist.equal_range(beginpolygon).first; segmentiter!=segmentlist.equal_range(beginpolygon).second; segmentiter++)
  {
    // int actkey = segmentiter->first;
    // std::vector<int> actsegmentpoints = segmentiter->second;

    // Its sufficent to determine the orientation of the first segment
    // to get the right direction of rotation (Umlaufsinn) of the polygon.
    IdentifyPolygonOrientation(segmentiter->second, segmentiter->first, intersectionpointsids, gfuncvalues);
    // IdentifyPolygonOrientation(actsegmentpoints, actkey, intersectionpointsids, gfuncvalues);
    std::vector<int> actsegmentpoints = segmentiter->second;
    //TEST
    // IO::cout << "Startsegment " << segmentiter->first << IO::endl;
    // IO::cout << actsegmentpoints[0] << IO::endl;
    // IO::cout << actsegmentpoints[1] << IO::endl;

    //contains the vertices of the polygon
    std::vector<int> polypoints;
    //store points of start segment
    polypoints.push_back(actsegmentpoints[0]);
    actpoint = actsegmentpoints[1];
    polypoints.push_back(actpoint);

    // now, the following segment is searched
    // this means the end point of this segment is equal the end point of the last segment
    // this procedure is repeated until the starting point is reached
    while(actpoint!=polypoints[0])
    {
      for (std::multimap<int,std::vector<int> >::const_iterator iter = segmentlist.begin(); iter != segmentlist.end(); ++iter)
      {
        std::vector<int> it_points = iter->second;
        if (it_points!=actsegmentpoints)
        {
          if (it_points[0]==actpoint)
          {
            actpoint = it_points[1];
            polypoints.push_back(actpoint);
            actsegmentpoints = it_points;
            break; //finished: find next segment
          }
          else if (it_points[1]==actpoint)
          {
            actpoint = it_points[0];
            polypoints.push_back(actpoint);
            actsegmentpoints = it_points;
            break; //finished: find next segment
          }
          else
          {
          }
        }
        else
        {
        }
      }
    }
    //store polypoints in a map containing all polygons
    polygonpoints[j] = polypoints;
    j++;
    //TEST
    //      IO::cout<<"polygonvertices"<< IO::endl;
    //      for (std::size_t ipoly=0; ipoly<polypoints.size(); ipoly++)
    //      {
    //        IO::cout<< polypoints[ipoly] << IO::endl;
    //      }

    // special case: both segments of the double intersected surface belong to the same polygon
    if (numpolygon==2)
    {
      if (polypoints.size()>(segmentlist.size()-3))
      {
        //IO::cout<< "G-function values" << IO::endl;
        //for (std::size_t i=0; i<gfuncvalues.size(); i++)
        //{
        //  IO::cout<< gfuncvalues[i] << IO::endl;
        //}
        //dserror("Unexpected intersection");
        IO::cout << "/!\\ warning: special case of polygon occured" << IO::endl;
        break;
      }
    }
  }

  //--------------------------------------------------------------------------
  // build triangles
  //--------------------------------------------------------------------------
  for (std::size_t ipolygons=0; ipolygons<polygonpoints.size(); ipolygons++)
  {
    std::vector<int> polypoints = polygonpoints[ipolygons];
    if (polypoints.size()<4)
      dserror("TriangulateFlameFront needs at least 3 intersectionpoints");

    //-----------------------------
    // three points form a triangle
    //-----------------------------
    if (polypoints.size()==4) //there are only three different points and three points form a triangle
    {
      std::vector<int> trianglepoints (3);
      for (int i=0; i<3; i++)
      {
        trianglepoints[i] = polypoints[i];
      }
      trianglelist.push_back(trianglepoints);
    }
    //----------------------------------------------------
    // add a midpoint to the pointlist to create triangles
    //----------------------------------------------------
    else
    {
      //-----------------------------
      //calculate midpoint of polygon
      //-----------------------------
      std::size_t numpoints = polypoints.size() - 1;
      std::vector<double> midpoint (3);
      std::vector<double> point1 (3);
      std::vector<double> point2 (3);
      if (numpoints%2==0) // even
      {
        point1 = pointlist[polypoints[0]];
        point2 = pointlist[polypoints[numpoints/2]];
      }
      else // odd
      {
        point1 = pointlist[polypoints[0]];
        point2 = pointlist[polypoints[(numpoints+1)/2]];
      }

      for (int dim=0; dim<3; dim++)
      {
        // compute middle of this coordinate
        midpoint[dim] = (point2[dim] + point1[dim]) * 0.5;
        // remark: modification added by Ursula (31.07.09)
        //         to avoid problems with TetGen, if triangles lie competely on a surface
        if (midpoint[dim]==1 or midpoint[dim]==-1) // midpoint lies on surface of the cell
        {
          for (unsigned idim=0;idim<gfuncvalues.size();idim++)
            IO::cout << "G-function value " << idim << " for this cell: " << gfuncvalues[idim] << IO::endl;
          IO::cout << "/!\\ warning === coodinate of midpoint moved to interior " << midpoint[dim] << IO::endl;
          // TODO: this should not be an absolute factor, but a relative factor (related to the cell size))
          midpoint[dim] = midpoint[dim] * 0.98;  // -> move midpoint into the cell
        }
      }

#if 0
      //--------------------------------------------------------------------------------
      // perform Newton-Raphson method to project midpoint on level set zero iso-surface
      //--------------------------------------------------------------------------------
      bool converged = false;
      // projection of midpoint vector
      LINALG::Matrix<3,1> projpoint(true);

      // get vector of G-function values of root cell
      const std::vector<double>& gfuncvaluesrootcell = cell->ReturnRootCell()->GetGfuncValues();

      switch(cell->Shape())
      {
      case DRT::Element::tet4:
      {
#ifdef DEBUG
        if (gfuncvaluesrootcell.size() != 4)
          dserror("discretization types do not match!");
#endif
        // store G-function values in a fixed size vector
        const size_t numvertices = DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::tet4>::numNodePerElement;
        LINALG::Matrix<numvertices,1> valuesGcell(true);
        for (unsigned irow=0;irow<gfuncvaluesrootcell.size();++irow)
          valuesGcell(irow) = gfuncvaluesrootcell[irow];

        // project midpoint on level set zero iso-surface
        converged = projectMidpoint<DRT::Element::tet4>(valuesGcell, midpoint, projpoint);
        break;
      }
      case DRT::Element::hex8:
      {
#ifdef DEBUG
        if (gfuncvaluesrootcell.size() != 8)
          dserror("discretization types do not match!");
#endif
        // store G-function values in a fixed size vector
        const size_t numvertices = DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::hex8>::numNodePerElement;
        LINALG::Matrix<numvertices,1> valuesGcell(true);
        for (unsigned irow=0;irow<gfuncvaluesrootcell.size();++irow)
          valuesGcell(irow) = gfuncvaluesrootcell[irow];

        // project midpoint on level set zero iso-surface
        converged = projectMidpoint<DRT::Element::hex8>(valuesGcell, midpoint, projpoint);
        break;
      }
      case DRT::Element::hex20:
      {
#ifdef DEBUG
        if (gfuncvaluesrootcell.size() != 20)
          dserror("discretization types do not match!");
#endif
        // store G-function values in a fixed size vector
        const size_t numvertices = DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::hex20>::numNodePerElement;
        LINALG::Matrix<numvertices,1> valuesGcell(true);
        for (unsigned irow=0;irow<gfuncvaluesrootcell.size();++irow)
          valuesGcell(irow) = gfuncvaluesrootcell[irow];

        // project midpoint on level set zero iso-surface
        converged = projectMidpoint<DRT::Element::hex20>(valuesGcell, midpoint, projpoint);
        break;
      }
      default:
        dserror("unknown type of boundary integration cell");
      }
      if (converged)
      {
        // overwrite midpoint with projection of midpoint
        for (int dim=0; dim<3; dim++)
          midpoint[dim] = projpoint(dim);
      }
#endif

//#ifndef COMBUST_2D
      // add midpoint to list of points defining piecewise linear complex (interface surface)
      std::size_t midpoint_id = pointlist.size();//ids start at 0
      //IO::cout<< "pointlistsize " << pointlist.size() << "midpointid " << midpoint_id<<IO::endl;
      pointlist.push_back(midpoint);

      //build triangles
      for (std::size_t j=0; j<polypoints.size()-1; j++)
      {
        std::vector<int> trianglepoints (3);
        trianglepoints[0] = polypoints[j];
        trianglepoints[1] = polypoints[j+1];
        trianglepoints[2] = midpoint_id;
        trianglelist.push_back(trianglepoints);
      }
//#else
//      projectMidpoint2D(trianglelist,segmentlist,pointlist,polypoints,midpoint);
//#endif
    }
  }

  //TEST
  //IO::cout<<"triangles"<< IO::endl;
  //for (std::size_t itriangle=0; itriangle<trianglelist.size(); itriangle++)
  //{
  //  IO::cout<< "dreieck " << itriangle<< IO::endl;
  //  for (int i=0; i<3; i++)
  //    IO::cout<< trianglelist[itriangle][i]<<IO::endl;
  //}

  return;
}


/*------------------------------------------------------------------------------------------------*
 | identify the orientation of the interface polygon, normal vector + -> -                        |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::FlameFront::IdentifyPolygonOrientation(
    std::vector<int>&           segment,
    const int                   surf_id,
    std::multimap<int,int>&     intersectionpointsids,
    const std::vector<double>&  gfuncvalues
)
{
  // array containing the line_id's for each surface
  int surfacelines[6][4] = {{3, 2, 1, 0},
      {0, 5, 8, 4},
      {1, 6, 9, 5},
      {2, 7,10, 6},
      {4,11, 7, 3},
      {8, 9,10,11}};
  // array containing test vertex for each surface
  // that means : test of sign of g-func
  int testnode[6][4] = {{0, 3, 2, 1},
      {0, 1, 5, 4},
      {1, 2, 6, 5},
      {2, 3, 7, 6},
      {0, 4, 7, 3},
      {4, 5, 6, 7}};
  int line_id = -1;
  int testpoint = -1;

  const int numnode = DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::hex8>::numNodePerElement;
  //segment contains intersection point
  if(segment[0]>numnode-1)
  {
    int interpoint = segment[0];

    for(std::map<int,int>::const_iterator iterinterpoint=intersectionpointsids.begin(); iterinterpoint!=intersectionpointsids.end(); iterinterpoint++)
    {
      if(iterinterpoint->second==interpoint)
      {
        line_id = iterinterpoint->first;
      }
    }
    //get point to test, depends on the line_id
    for(int k=0; k<4; k++)
    {
      if(line_id==surfacelines[surf_id][k])
      {
        testpoint = k;
        break;
      }
    }
  }
  else //segment == diagonal
  {
    for(int k=0; k<4; k++)
    {
      if(testnode[surf_id][k]==segment[0])
      {
        if(k==0)
          testpoint = 3;
        else
          testpoint = k-1;
        break;
      }
    }
  }
  if(gfuncvalues[testnode[surf_id][testpoint]]>0) //change! (normal vector in domain plus <0 is needed)
  {
    int temp = segment[0];
    segment[0] = segment[1];
    segment[1] = temp;
  }
  return;
}


/*------------------------------------------------------------------------------------------------*
 | build polygon segments of intersection surface of a refinement cell                            |
 | hex8 only                                                                                      |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::FlameFront::buildFlameFrontSegments(
    std::multimap<int,int>&                  intersectionpointsids,
    std::multimap<int,std::vector<int> >&    segmentlist,
    const std::vector<double>&               gfuncvalues,
    const std::vector<std::vector<double> >& pointlist
)
{
  // array containing the line_id's for each surface
  int surface[6][4] = {{3, 2, 1, 0},
      {0, 5, 8, 4},
      {1, 6, 9, 5},
      {2, 7,10, 6},
      {4,11, 7, 3},
      {8, 9,10,11}};

  //int numsurf = DRT::UTILS::getNumberOfElementSurfaces(DRT::Element::hex8);
  //hex8 only
  for (int i=0; i<6; i++)//loop over all surfaces
  {
    // gets end points of the segment
    std::vector<int> segmentpoints;

    for (int j=0; j<4; j++)//loop over all lines of the surface
    {
      // find line
      std::map<int,int>::const_iterator it_intersectionpoint = intersectionpointsids.find(surface[i][j]);
      if (it_intersectionpoint!=intersectionpointsids.end()) //check: is line intersected?
      {
        //get id of the intersectionpoint and store it in segmentpoints vector
        segmentpoints.push_back(it_intersectionpoint->second);
      }
    }

    //IO::cout << "segmentpoints size" << segmentpoints.size() << IO::endl;

    switch (segmentpoints.size())
    {
    case 0:
    {
      // check, if gfunc==0 at the vertices
      // gets nodes with gfunc==0
      std::vector<int> zeropoints;
      std::vector<std::vector<int> > surfacepointlist = DRT::UTILS::getEleNodeNumberingSurfaces(DRT::Element::hex8);
      for(int k=0; k<4; k++)//loop over nodes
      {
        if (gfuncvalues[surfacepointlist[i][k]]==0)
          zeropoints.push_back(surfacepointlist[i][k]);
      }
      //IO::cout << "zeropoints size" << zeropoints.size() << IO::endl;
      switch (zeropoints.size())
      {
      case 0: //surface not intersected
      case 1: //surface not intersected, but one vertex touchs the interface
      {
        break;
      }
      case 2:
        /*
         * 2 opposite vertices with Phi=0
         * - surface is intersected, iff Phi at the remaining vertices has different sign (segment == diagonal)
         * - otherwise vertiches touch interface
         * 2 neighboring vertices with Phi=0
         * - edge touchs interface, but cell is not intersected
         * - edge limits interface patch -> then, it is needed in TriangulateFlameFront()
         */
      {
        // intersection: opposite vertices with Phi=0
        if(((zeropoints[0]==surfacepointlist[i][0])and(zeropoints[1]==surfacepointlist[i][2])))
        {
          // different sign of Phi at the remaining vertices
          if (gfuncvalues[surfacepointlist[i][1]]*gfuncvalues[surfacepointlist[i][3]]<0)
          {
            //store in segmentlist
            segmentlist.insert(std::pair<int,std::vector<int> >(i,zeropoints));
          }
        }
        else if (((zeropoints[0]==surfacepointlist[i][1])and(zeropoints[1]==surfacepointlist[i][3])))
        {
          if (gfuncvalues[surfacepointlist[i][0]]*gfuncvalues[surfacepointlist[i][2]]<0)
          {
            segmentlist.insert(std::pair<int,std::vector<int> >(i,zeropoints));
          }
        }
        else
        {
          // no intersection
          // possibly segment limits the interface
          // store in segment list (key=-1), if it is not already contained in the segmentlist
          if (zeropoints[0]>zeropoints[1])
          {
            int temp =zeropoints[1];
            zeropoints[1] = zeropoints[0];
            zeropoints[0] = temp;
          }
          //already in segmentlist?
          bool not_in_segmentlist = true;
          for (std::multimap<int,std::vector<int> >::iterator iter=segmentlist.equal_range(-1).first; iter!=segmentlist.equal_range(-1).second; iter++)
          {
            if (iter->second == zeropoints)
              not_in_segmentlist = false;
          }
          if(not_in_segmentlist)
            //store in segmentlist
            segmentlist.insert(std::pair<int,std::vector<int> >(-1,zeropoints));
        }
        break;
      }
      case 3:
      {
        // surface is not intersected, but two following edges are aligned with the interface
        // store in segment list (key=-1), if it is not already contained in the segmentlist

        // zeropoints are assigned to segments
        // checke which vertex is missing
        if(zeropoints[0]==surfacepointlist[i][0])
        {
          if(zeropoints[1]==surfacepointlist[i][1])
          {
            if(zeropoints[2]==surfacepointlist[i][2])
            {
              //vertex 4 is missing -> no problem
            }
            else
            {
              //vertex 3 is missing -> rearrange zeropoints
              int temp = zeropoints[0];
              zeropoints[0] = zeropoints[2];
              zeropoints[2] = zeropoints[1];
              zeropoints[1] = temp;
            }
          }
          else
          {
            //vertex 2 is missing -> rearrange zeropoints
            int temp = zeropoints[0];
            zeropoints[0] = zeropoints[1];
            zeropoints[1] = zeropoints[2];
            zeropoints[2] = temp;
          }
        }
        else
        {
          //vertex 1 is missing -> no problem
        }

        for (std::size_t k=0; k<zeropoints.size()-1; k++)
        {
          std::vector<int> segment (2);
          segment[0] = zeropoints[k];
          segment[1] = zeropoints[k+1];
          if (segment[0]>segment[1])
          {
            int temp =segment[1];
            segment[1] = segment[0];
            segment[0] = temp;
          }
          bool not_in_segmentlist = true;
          for (std::multimap<int,std::vector<int> >::iterator iter=segmentlist.equal_range(-1).first; iter!=segmentlist.equal_range(-1).second; iter++)
          {
            if (iter->second == segment)
              not_in_segmentlist = false;
          }
          if(not_in_segmentlist)
            segmentlist.insert(std::pair<int,std::vector<int> >(-1,segment));
        }
        break;
      }
      case 4:
        // surface is aligned with the interface -> cell is not intersected
        // surface has to be added to BoundaryIntCells
      {

        // as all segments, limiting the interface, are already found in case 2 and 3,
        // nothing remains to do here

        break;
      }
      default:
        dserror("impossible number of zero values");
      }
      break;
    }
    case 1:
    {
      //IO::cout << "one intersection point" << IO::endl;
      //find zeropoint to build segment
      std::vector<int> zeropoints;
      std::vector<std::vector<int> > surfacepointlist = DRT::UTILS::getEleNodeNumberingSurfaces(DRT::Element::hex8);
      for(int k=0; k<4; k++)//loop over nodes
      {
        if (gfuncvalues[surfacepointlist[i][k]]==0.0)
          zeropoints.push_back(surfacepointlist[i][k]);
      }
      if (zeropoints.size()!=1)
      {
        IO::cout << "Surface " << i << IO::endl;
        for (size_t l=0; l < gfuncvalues.size(); l++)
        {
          IO::cout << "Node " << l << " " << gfuncvalues[l] << IO::endl;
        }
        dserror("can't build intersection segment, one intersection point, but not one vertex that has G=0");
      }

      std::vector<int> segment (2);
      segment[0] = segmentpoints[0];
      segment[1] = zeropoints[0];
      segmentlist.insert(std::pair<int,std::vector<int> >(i,segment));

      break;
    }
    case 2:
    {
      //store segment in segmentlist
      segmentlist.insert(std::pair<int,std::vector<int> >(i,segmentpoints));
      break;
    }
    case 3:
    {
      dserror("3 is an impossible number of intersectionpoints for hex8 element surface");
      break;
    }
    case 4:
    {
      //IO::cout << "4 intersection points " << IO::endl;
      // 2 segments per surface are possible for hex8 -> 4 intersection points
      // idea to assign the intersection points to the segments:
      //       look for maximal distance
      //       this combination is not possible
      //       the segments are now obtained
      double distance;
      double maxdist = 0.0;
      int maxdistcounter = 0;
      for(int k=0; k<4; k++)//loop over segmentpoints
      {
        std::vector<double> point1 = pointlist[segmentpoints[k]];
        std::vector<double> point2 (3);
        if(k<3)
        {
          point2 = pointlist[segmentpoints[k+1]];
        }
        else
        {
          point2 = pointlist[segmentpoints[0]];
        }
        distance = sqrt((point1[0]-point2[0])*(point1[0]-point2[0])+(point1[1]-point2[1])*(point1[1]-point2[1])+(point1[2]-point2[2])*(point1[2]-point2[2]));
        if (distance>maxdist)
        {
          maxdist = distance;
          maxdistcounter = k;
        }
      }

      //build segment
      //IO::cout << "Maxdist" << maxdistcounter << IO::endl;
      if (maxdistcounter==0 or maxdistcounter==2)
      {
        std::vector<int> segment1 (2);
        segment1[0] = segmentpoints[3];
        segment1[1] = segmentpoints[0];
        std::vector<int> segment2 (2);
        segment2[0] = segmentpoints[1];
        segment2[1] = segmentpoints[2];
        segmentlist.insert(std::pair<int,std::vector<int> >(i,segment1));
        segmentlist.insert(std::pair<int,std::vector<int> >(i,segment2));
      }
      else
      {
        std::vector<int> segment1 (2);
        segment1[0] = segmentpoints[0];
        segment1[1] = segmentpoints[1];
        std::vector<int> segment2 (2);
        segment2[0] = segmentpoints[2];
        segment2[1] = segmentpoints[3];
        segmentlist.insert(std::pair<int,std::vector<int> >(i,segment1));
        segmentlist.insert(std::pair<int,std::vector<int> >(i,segment2));
      }
      break;
    }
    default:
    {
      dserror("unexpected number of intersectionpoints!");
      break;
    }
    }
  }

  return;
}


/*------------------------------------------------------------------------------------------------*
 | store domain integration cell (only for hex8 cells at present)                     henke 08/10 |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::FlameFront::StoreDomainIntegrationCell(
    const COMBUST::RefinementCell* cell,
    GEO::DomainIntCells& domainintcelllist
)
{
  //if (cell->Shape() != DRT::Element::hex8) dserror("check if this makes sense for non-hex8 elements first");
  // get G-function values at vertices from cell
  const std::vector<double>& gfuncvalues = cell->GetGfuncValues();
  // get vertex coordinates (local fluid element coordinates) from refinement cell
  const std::vector<std::vector<double> >& vertexcoord = cell->GetVertexCoord();

  //---------------------------------
  // get global coordinates of nodes
  //---------------------------------
  const int numnode = cell->Ele()->NumNode();
  LINALG::SerialDenseMatrix xyze(3,numnode);
  for(int inode=0;inode<numnode;inode++)
  {
    xyze(0,inode) = cell->Ele()->Nodes()[inode]->X()[0];
    xyze(1,inode) = cell->Ele()->Nodes()[inode]->X()[1];
    xyze(2,inode) = cell->Ele()->Nodes()[inode]->X()[2];
  }

  //-----------------------------------------
  // get global coordinates of cell vertices
  //-----------------------------------------
  //coordinates of vertices of cell, corresponds to 'vertexcoord')
  LINALG::SerialDenseMatrix cellcoord(3,numnode);
  // if cell == element, globalcellcoord == xyze
  LINALG::SerialDenseMatrix globalcellcoord(3,numnode);

  for(int ivert=0; ivert<numnode; ivert++)
  {
    static LINALG::Matrix<3,1> vertcoord;
    for(int dim=0; dim<3; dim++)
    {
      // cellcoord = transpose of vertexcoord
      cellcoord(dim,ivert) = vertexcoord[ivert][dim];
      // coordinates of a single cell vertex
      vertcoord(dim) = cellcoord(dim,ivert);
    }
    // transform vertex from local (element) coordinates to global (physical) coordinates
    GEO::elementToCurrentCoordinatesInPlace(cell->Ele()->Shape(), xyze, vertcoord);

    // store vertex in array of global cell coordinates
    for(int  dim=0; dim<3; dim++)
      globalcellcoord(dim,ivert) = vertcoord(dim);
#ifdef DEBUG
    // if cell = element -> globalcellcoord = xyze
    //if (DRT::INPUT::IntegralValue<int>(combustdyn.sublist("COMBUSTION GFUNCTION"),"REFINEMENT") == false)
    //  if (globalcellcoord != xyze) dserror("coordinates should be identical");
#endif
  }

  //-------------------------------------------
  // determine which domain the cell belongs to
  //-------------------------------------------
  // compute average G-function value for this refinement cell (= integration cell)
  bool inGplus = GetIntCellDomain(cellcoord,gfuncvalues,cell->Shape());

  //TEST
  //IO::cout << "globalcellcoord " << globalcellcoord(0,3) << globalcellcoord(1,3) << IO::endl;
  //if(inGplus) {
  //  IO::cout << "In G plus" << IO::endl;
  //}
  //else {
  //  IO::cout << "In G minus" << IO::endl;
  //}
  //------------------------
  // create integration cell
  //------------------------
  // create an integration cell and add it to the list of integration cells per element
  // remark: for now, this is restricted to hex8 integration cells
  domainintcelllist.push_back(GEO::DomainIntCell(cell->Shape(), cellcoord, globalcellcoord, inGplus));

  return;
}


/*------------------------------------------------------------------------------------------------*
 | build piecewise linear complex (PLC) in Tetgen format                              henke 10/08 |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::FlameFront::buildPLC(
    const COMBUST::RefinementCell* cell,
    const std::vector<double>& gfuncvaluesrootcell,
    GEO::DomainIntCells& domainintcelllist,
    GEO::BoundaryIntCells& boundaryintcelllist)
{
  /* A PLC contains
   * - points -> vertices, intersection points, center of interface patch
   * - segmentlist -> intersection curve between interface and surfaces of cell: buildFlameFrontSegments()
   * - trianglelist -> triangles approximating the interface: TriangulateFlameFront()
   * A PLC is the input for TetGen -> CallTetGen()
   * The triangles from the triangle list are stored as boundary integration cells. TetGaen also
   * generates triangular boundary integration cells. However, these are not used since their
   * configuaration (vertex numbering) is unknown.
   */

  //-----------------------------------------------------
  // prepare lists, get coordinates and G-function values
  //-----------------------------------------------------
  std::vector<std::vector<double> >    pointlist;
  std::multimap<int,std::vector<int> > segmentlist;
  std::vector<std::vector<int> >       trianglelist;

  const DRT::Element::DiscretizationType distype = cell->Ele()->Shape();
  // global coordinates of element this cell belongs to
  const int numnode = cell->Ele()->NumNode();
  LINALG::SerialDenseMatrix xyze(3,numnode);
  for(int inode=0;inode<numnode;inode++)
  {
    xyze(0,inode) = cell->Ele()->Nodes()[inode]->X()[0];
    xyze(1,inode) = cell->Ele()->Nodes()[inode]->X()[1];
    xyze(2,inode) = cell->Ele()->Nodes()[inode]->X()[2];
  }

  // get vertex coordinates (local fluid element coordinates) from refinement cell
  const std::vector<std::vector<double> >& vertexcoord = cell->GetVertexCoord();

  // get G-function values from refinement cell
  const std::vector<double>& gfuncvalues = cell->GetGfuncValues();

  //TEST
  //IO::cout << cell->Ele()->Id() << IO::endl;
  //IO::cout<< "G-Werte" << IO::endl;
  //for (std::size_t i=0; i<gfuncvalues.size(); i++)
  //{
  //  IO::cout<< gfuncvalues[i] << IO::endl;
  //}

  //---------------------------------------------
  // fill list of points and intersection points
  //---------------------------------------------
  int numofpoints = 0;
  int numvertex = DRT::UTILS::getNumberOfElementCornerNodes(distype);
  // get intersection points from refinement cell

  // corner points
  for (int ivertex=0; ivertex<numvertex; ivertex++)
  {
    // remark: all hex elements/cells have 8 corners;
    //         additional nodes (hex20,hex27) are inner nodes and not corners
    //         here, only the corners are written into the point list
    pointlist.push_back(vertexcoord[ivertex]);
    numofpoints++;
  }

  const std::multimap<int,std::vector<double> >& intersectionpoints = cell->intersectionpoints_;
  // map corresponds to intersection points, but contains the points' IDs instead of their coordinates
  // links 'pointlist' required by TetGen to 'intersectionpointlist'
  std::multimap<int,int> intersectionpointsids;

  // intersection points
  // std::map<ID of cut edge in element, coordinates of intersection point>
  for (std::multimap<int,std::vector<double> >::const_iterator iter = intersectionpoints.begin(); iter != intersectionpoints.end(); ++iter)
  {
    pointlist.push_back(iter->second);
    intersectionpointsids.insert(std::pair<int, int>(iter->first,numofpoints));
    //intersectionpointsids[iter->first] = numofpoints;
    numofpoints++;
  }

  //TEST
  //IO::cout << "number of intersection points " << intersectionpoints.size() << IO::endl;
  //for (std::size_t iter=0; iter<pointlist.size(); ++iter)
  //{
  //  IO::cout<< iter << IO::endl;
  //  std::vector<double> coord = pointlist[iter];
  //  for (std::size_t isd=0; isd<3; isd++)
  //  {
  //    IO::cout<< coord[isd] << IO::endl;
  //  }
  //}

  //---------------------------------------------------------
  // build segments enclosing interface patches within a cell
  //---------------------------------------------------------
  // Segments form the boundary of the interfacepatch in the cell (element, without refinement).
  // Some special cases have to be distinguished
  buildFlameFrontSegments(intersectionpointsids, segmentlist, gfuncvalues, pointlist);


  //TEST
  //IO::cout<<"Segments"<< IO::endl;
  //for (std::map<int,std::vector<int> >::const_iterator iter = segmentlist.begin(); iter != segmentlist.end(); ++iter)
  //{
  //  IO::cout<<"Segment "<< iter->first << IO::endl;
  //  std::vector<int> point = iter->second;
  //  for (std::size_t isd=0; isd<2; isd++)
  //  {
  //    IO::cout<< point[isd] << IO::endl;
  //  }
  //}

  //--------------------------------------
  // triangulate flame front within a cell
  //--------------------------------------
  TriangulateFlameFront(cell, trianglelist, segmentlist, pointlist, intersectionpointsids, gfuncvalues);
  if (trianglelist.size() <= 0) dserror("There are no triangles!");

  //--------------------------------------------
  // store triangular boundary integration cells
  //--------------------------------------------
  // remark: Although Tetgen also delivers triangles as interface patches, these are not used,
  //         because their direction of rotation (Umlaufsinn) is not known. Instead, the triangles
  //         from the 'trianglelist' are used directly.
  for (std::size_t itriangle=0; itriangle<trianglelist.size(); itriangle++)
  {
    LINALG::SerialDenseMatrix trianglecoord(3,3); //3 directions, 3 nodes
    LINALG::SerialDenseMatrix phystrianglecoord(3,3);
    for (int inode=0; inode<3; inode++)
    {
      static LINALG::Matrix<3,1> tcoord;
      for (int dim=0; dim<3; dim++)
      {
        trianglecoord(dim,inode) = pointlist[trianglelist[itriangle][inode]][dim];
        tcoord(dim) = trianglecoord(dim,inode);
      }
      GEO::elementToCurrentCoordinatesInPlace(distype, xyze, tcoord);
      for(int  dim=0; dim<3; dim++)
        phystrianglecoord(dim,inode) = tcoord(dim);
    }

    // store boundary integration cells in boundaryintcelllist
    // remark: boundary integration cells in bisected cells (these have always triangular boundary
    //         integration cells) are defined to belong to the "plus" domain (G>0)
    boundaryintcelllist.push_back(GEO::BoundaryIntCell(DRT::Element::tri3, -1, trianglecoord,
        Teuchos::null, phystrianglecoord, true));
  }

#if 0
  if (xfeminttype_ == INPAR::COMBUST::xfemintegration_cut)
  {
    DRT::Element::DiscretizationType distype_cell = cell->Shape();

    // get vertex coordinates (local fluid element coordinates) from refinement cell
    const std::vector<std::vector<double> >& vertexcoord = cell->GetVertexCoord();

    int numnode = DRT::UTILS::getNumberOfElementNodes( distype_cell );

    std::vector<int> nids;
    nids.reserve(numnode);
    for ( int i=0; i<numnode; ++i )
    {
      nids.push_back(i);
    }

    //-----------------------------------------
    // get global coordinates of cell vertices
    //-----------------------------------------

    LINALG::SerialDenseMatrix globalcellcoord(3,numnode);

    for (int ivert=0; ivert<numnode; ivert++)
    {
      LINALG::Matrix<3,1> vertcoord;

      std::copy(vertexcoord[ivert].begin(),
          vertexcoord[ivert].end(),
          vertcoord.A());

      // transform vertex from local (element) coordinates to global (physical) coordinates
      GEO::elementToCurrentCoordinatesInPlace(distype_cell, xyze, vertcoord);

      std::copy(vertcoord.A(),
          vertcoord.A()+3,
          &globalcellcoord(0,ivert));
    }

    // get G-function values at vertices from cell
    const std::vector<double>& gfuncvalues = cell->GetGfuncValues();

    // GEO::CUT intersection call

    GEO::CUT::LevelSetIntersection levelset;

    levelset.AddElement( 1, nids, globalcellcoord, distype_cell, &gfuncvalues[0] );

    try
    {
      levelset.Cut();
    }
    catch ( std::runtime_error & err )
    {
      std::cerr << "failed to cut element\n"
                << "coordinates:\n"
        ;
      globalcellcoord.Print( std::cerr );
      std::cerr << "g-function values:\n"
                << std::setprecision( 16 );
      std::copy( gfuncvalues.begin(), gfuncvalues.end(), std::ostream_iterator<double>( std::cerr, ", " ) );
      std::cerr << "\n";
      throw;
    }

    GEO::CUT::ElementHandle * e = levelset.GetElement( 1 );

    if ( e!=NULL and e->IsCut() )
    {
      GEO::CUT::plain_integrationcell_set cells;
      e->GetIntegrationCells( cells );

      //std::set<GEO::CUT::BoundaryCell*> bcells;
      //e->GetBoundaryCells( bcells );

      LINALG::Matrix<3,1> physCoordCorner;
      LINALG::Matrix<3,1> eleCoordDomainCorner;

      for ( GEO::CUT::plain_integrationcell_set::iterator i=cells.begin(); i!=cells.end(); ++i )
      {
        GEO::CUT::IntegrationCell * ic = *i;
        DRT::Element::DiscretizationType distype = ic->Shape();
        int numnodes = DRT::UTILS::getNumberOfElementNodes( distype );

        LINALG::SerialDenseMatrix physCoord = ic->Coordinates();
        LINALG::SerialDenseMatrix coord( 3, numnodes );

        for ( int j=0; j<numnodes; ++j )
        {
          std::copy( &physCoord( 0, j ), &physCoord( 0, j ) + 3, physCoordCorner.A() );

          e->LocalCoordinates( physCoordCorner, eleCoordDomainCorner );

          std::copy( eleCoordDomainCorner.A(), eleCoordDomainCorner.A()+3, &coord( 0, j ) );
        }

        // if degenerated don't store
        if ( distype==DRT::Element::tet4 )
          if(GEO::checkDegenerateTet(numnodes, coord, physCoord))
            continue;

        bool inGplus = GetIntCellDomainInElement(coord, gfuncvalues, DRT::Element::hex8, distype);

        //domainintcelllist.push_back(GEO::DomainIntCell(distype, tetrahedroncoord, phystetrahedroncoord, inGplus));
        domainintcelllist.push_back( GEO::DomainIntCell( distype, coord, physCoord, inGplus /*ic->Position()==GEO::CUT::Point::outside*/ ) );
      }
    }
    else
    {
      dserror("cut expected");
    }
  }
#endif
}


/*------------------------------------------------------------------------------------------------*
 | store (quad4) boundary integration cell                                            henke 08/10 |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::FlameFront::StoreBoundaryIntegrationCell(
    const COMBUST::RefinementCell* cell,
    GEO::BoundaryIntCells& boundaryintcelllist
)
{
  // For a given refinement cell, a boundary integartion cell is created, if one surface of the cell
  // is fully aligned with the zero iso-surface of the interface (flame front).
  // For touched refiement cells, boundary integration cells are stored twice, that is for the
  // refinement cells on both sides of the interface with the corresponding boolean 'inGplus'
  // indicating which side they belong to.

  // check cell shape
  DRT::Element::DiscretizationType celldistype = cell->Shape();
#ifdef DEBUG
  if (cell->Shape()!=DRT::Element::hex8)
    dserror("not supported for this cell shape");
  // TODO not sure if this is really a prerequisite for this function
  if (cell->Ele()->Shape()!=DRT::Element::hex8)
    dserror("not supported for this element shape");
#endif

  // get G-function values from refinement cell
  const std::vector<double>& gfuncvaluescell = cell->GetGfuncValues();
  // get vertex coordinates (local fluid element coordinates) from refinement cell
  const std::vector<std::vector<double> >& vertexcoord = cell->GetVertexCoord();
  //---------------------------------
  // get global coordinates of nodes
  //---------------------------------
  const int numnode = cell->Ele()->NumNode();
  LINALG::SerialDenseMatrix xyze(3,numnode);
  for(int inode=0;inode<numnode;inode++)
  {
    xyze(0,inode) = cell->Ele()->Nodes()[inode]->X()[0];
    xyze(1,inode) = cell->Ele()->Nodes()[inode]->X()[1];
    xyze(2,inode) = cell->Ele()->Nodes()[inode]->X()[2];
  }

  // get list of surface IDs of cell
  const std::vector<std::vector<int> >& cellsurfacelist = DRT::UTILS::getEleNodeNumberingSurfaces(celldistype);
  int numsurf = DRT::UTILS::getNumberOfElementSurfaces(celldistype);
  for(int isurf=0; isurf<numsurf; isurf++)
  {
    int numzeronodes = 0;
    // loop over vertices
    for(int k=0; k<4; k++)
    {
      // remark: We do not work with tolerances here, since the G-function values of vertices close
      //         to the interface have been reset to zero previously. Therefore they should relly
      //         be 0.0 exactly.
      if (gfuncvaluescell[cellsurfacelist[isurf][k]]==0.0)
        numzeronodes++;
    }
    // Only if the cell surface is fully aligned with the interface (4 vertices with zero G-function
    // value), a boundary integraton cell is stored.
    if (numzeronodes==4)
    {
      // define coordinate matrices (3 dimensions, 4 nodes)
      LINALG::SerialDenseMatrix quadcoord(3,4);
      LINALG::SerialDenseMatrix physquadcoord(3,4);
      for (int inode=0; inode<4; inode++)
      {
        static LINALG::Matrix<3,1> qcoord;
        for (int dim=0; dim<3; dim++)
        {
          // 'quadcoord' is transpose of 'vertexcoord'
          // remark: because the numbering in 'vertexcoord' coincides with the numbering in
          //         'cellsurfacelist', we can be sure to pick the right nodes, here
          quadcoord(dim,inode) = vertexcoord[cellsurfacelist[isurf][inode]][dim];
          qcoord(dim) = quadcoord(dim,inode);
        }
        GEO::elementToCurrentCoordinatesInPlace(DRT::Element::hex8, xyze, qcoord);
        for(int  dim=0; dim<3; dim++)
          physquadcoord(dim,inode) = qcoord(dim);
      }

      LINALG::SerialDenseMatrix cellcoord(3, numnode);
      for(int inode=0; inode<numnode; inode++)
      {
        for(int dim=0; dim<3; dim++)
          cellcoord(dim,inode) = vertexcoord[inode][dim];
      }

      // Determine which domain the domain integration cell belongs to - the corresponding boundary
      // integration cell will always belong to the same domain.
      bool inGplus = GetIntCellDomain(cellcoord, gfuncvaluescell,celldistype);
      //store boundary integration cells in boundaryintcelllist
      boundaryintcelllist.push_back(GEO::BoundaryIntCell(DRT::Element::quad4, -1, quadcoord,
          Teuchos::null, physquadcoord,inGplus));
    }
  }
}


/*------------------------------------------------------------------------------------------------*
 | store domain integration cell                                                  rasthofer 01/14 |
 *------------------------------------------------------------------------------------------------*/
std::size_t COMBUST::FlameFront::StoreDomainIntegrationCells(
    GEO::CUT::Element*         element,
    GEO::DomainIntCells&       domainintcelllist,
    LINALG::SerialDenseMatrix& xyze
)
{
  // number of stored volumes (not domain integration cells!)
  std::size_t numstoredvol = 0;

  // get volume cells
  GEO::CUT::plain_volumecell_set volcells;
  volcells = element->VolumeCells();

  //-------------------------------------
  // add domain integration cells to list
  //-------------------------------------

  // loop all volume cells
  for ( GEO::CUT::plain_volumecell_set::const_iterator ivolcell=volcells.begin(); ivolcell!=volcells.end(); ++ivolcell )
  {
    GEO::CUT::VolumeCell * volcell = *ivolcell;
    //std::cout << "volumen " << std::setw(24)<< std::setprecision(20) << volcell->Volume() << std::endl;
    //volcell->Print(std::cout);

    // get domain integration cells
    const GEO::CUT::plain_integrationcell_set cells = volcell->IntegrationCells();

    // check, if this volume cell has integration cells
    // this necessary here, since there is a tolerance in the cut algorithm which prevents TetMesh to tetrahedralize
    // volume cells smaller than 10^-14 (boundary cells are created!)
    std::size_t domcellcount = 0;
    if (cells.size()>0)
    {

      // loop domain integration cells
      for ( GEO::CUT::plain_integrationcell_set::const_iterator icell=cells.begin(); icell!=cells.end(); ++icell )
      {
        GEO::CUT::IntegrationCell * ic = *icell;
        DRT::Element::DiscretizationType distype = ic->Shape();
        int numnode = DRT::UTILS::getNumberOfElementNodes( distype );

        LINALG::SerialDenseMatrix physCoord(3,numnode);
        LINALG::SerialDenseMatrix coord = ic->Coordinates();

        for (int ivert=0; ivert<numnode; ivert++)
        {
          // write 'coord' to fixed size matrix format
          LINALG::Matrix<3,1> vertcoord;
          std::copy(&coord(0,ivert), &coord(0,ivert)+3, vertcoord.A());

          // transform vertex from local (element) coordinates to global (physical) coordinates
          GEO::elementToCurrentCoordinatesInPlace(element->Shape(), xyze, vertcoord);

          // write as 'physCoord'
          std::copy(vertcoord.A(), vertcoord.A()+3, &physCoord(0,ivert));
        }

        bool inGplus = false;
        if (ic->Position() == GEO::CUT::Point::outside) // volume cell is in G-plus domain (G>0)
          inGplus = true;

        if (distype != DRT::Element::hex8 and
            distype != DRT::Element::tet4)
        {
          IO::cout << "distype " << distype << IO::endl;
          dserror("unexpected type of domain integration cell");
        }

        double physvolume = GEO::DomainIntCell( distype, coord, physCoord, inGplus ).VolumeInPhysicalDomain();
        if (physvolume < 0.0)
          IO::cout << "negative volume detected for domain integration cell!" << IO::endl;

         domainintcelllist.push_back( GEO::DomainIntCell( distype, coord, physCoord, inGplus ) );
         domcellcount += 1;

      }

      numstoredvol += 1;
      // at least one domain integration cell has been stored
      if (domcellcount == 0) dserror("no domain integration cells for this volume cell");
    }

  }

  return numstoredvol;
}


/*------------------------------------------------------------------------------------------------*
 | store a list of boundary integration cells                                     rasthofer 01/14 |
 *------------------------------------------------------------------------------------------------*/
size_t COMBUST::FlameFront::StoreBoundaryIntegrationCells(
    GEO::CUT::Element*         element,
    GEO::BoundaryIntCells&     boundaryintcelllist,
    LINALG::SerialDenseMatrix& xyze
)
{
  // number of stored boundaries (not boundary integration cells!)
  size_t numstoredbound = 0;

  // get volume cells
  GEO::CUT::plain_volumecell_set volcells;
  volcells = element->VolumeCells();
  const int numvolcells = volcells.size();

  // the volume cells (G-plus or G-minus) we want to loop to determine what boundary cells to store
  GEO::CUT::Point::PointPosition loopvol = GEO::CUT::Point::undecided;

  // regular call of this function
  if(numvolcells!=1)
  {
    // compute G-plus and total volume
    double refvol = 0.0;
    double refvolplus = 0.0;
    for ( GEO::CUT::plain_volumecell_set::const_iterator ivolcell=volcells.begin(); ivolcell!=volcells.end(); ++ivolcell )
    {
      GEO::CUT::VolumeCell * volcell = *ivolcell;

      // compute positive (G>0) volume
      if (volcell->Position() == GEO::CUT::Point::outside)
      {
        refvolplus += volcell->Volume();
      }
      // volume of all volume cells (since we work in parameter space, the result should always be 2.0^3=8.0)
      refvol += volcell->Volume();
    }

    // remark: We want to loop those volume cells (plus or minus) occupying the larger portion of
    //         the element volume.
    // G-plus cells are larger
    if (refvolplus > (refvol-refvolplus))
      loopvol = GEO::CUT::Point::outside; // loop the volume cells in the G-plus domain (G>0)
    // G-minus cells are larger
    else
      loopvol = GEO::CUT::Point::inside; // loop the volume cells in the G-minus domain (G<0)
  }
  else // the touched case
  {
    // loop the G-plus cells
    // remark: if this is not a G-plus cell, the loop will do nothing
    //         there will be a G-plus neighbor element which will contribute its boundary cell so
    //         we create no hole in the interface
    loopvol = GEO::CUT::Point::outside;
  }

  //---------------------------------------
  // add boundary integration cells to list
  //---------------------------------------
  // temporary vectors
  LINALG::Matrix<3,1> physCoordVertex;
  LINALG::Matrix<3,1> eleCoordVertex;

  for ( GEO::CUT::plain_volumecell_set::const_iterator ivolcell=volcells.begin(); ivolcell!=volcells.end(); ++ivolcell )
  {
    GEO::CUT::VolumeCell * volcell = *ivolcell;

    // if the cell belongs to the domain we want to loop
    if (volcell->Position() == loopvol)
    {
      // get boundary integration cells
      const GEO::CUT::plain_boundarycell_set & bcells = volcell->BoundaryCells();

      std::size_t bcellcount = 0;

      for ( GEO::CUT::plain_boundarycell_set::const_iterator ibcell=bcells.begin(); ibcell!=bcells.end(); ++ibcell )
      {
        GEO::CUT::BoundaryCell * bcell = *ibcell;

        DRT::Element::DiscretizationType distype = bcell->Shape();
        int numnode = DRT::UTILS::getNumberOfElementNodes( distype );

        if (distype != DRT::Element::tri3 and
            distype != DRT::Element::quad4)
        {
          IO::cout << "distype " << distype << IO::endl;
          dserror("unexpected type of boundary integration cell");
        }

        // get physical coordinates of this cell
        LINALG::SerialDenseMatrix physCoord(3,numnode);
        // local coordinates of this cell
        LINALG::SerialDenseMatrix coord = bcell->Coordinates();

        // for negative (G<0) volume cells, we have to reverse the order of vertices to reverse
        // the normal vector
        if(loopvol == GEO::CUT::Point::inside)
        {
          LINALG::SerialDenseMatrix tmpcoord = coord;
          for (int ivert=0; ivert<numnode; ivert++)
            std::copy(&tmpcoord(0,ivert), &tmpcoord(0,ivert)+3, &coord(0,numnode-ivert-1));
        }

        for (int ivert=0; ivert<numnode; ivert++)
        {
          // write 'coord' to fixed size matrix format
          LINALG::Matrix<3,1> vertcoord;
          std::copy(&coord(0,ivert), &coord(0,ivert)+3, vertcoord.A());

          // transform vertex from local (element) coordinates to global (physical) coordinates
          GEO::elementToCurrentCoordinatesInPlace(element->Shape(), xyze, vertcoord);

          // write as 'physCoord'
          std::copy(vertcoord.A(), vertcoord.A()+3, &physCoord(0,ivert));
        }

        boundaryintcelllist.push_back(GEO::BoundaryIntCell(distype, -1, coord, Teuchos::null, physCoord, true));
        bcellcount += 1;

      }

      // at least one boundary integration cell has been stored
      if (bcellcount==0) dserror("no boundary integration cells for this volume cell");

      if (numvolcells >= 3)
      {
        // remark: Trisected elements consist of three volume cells and two boundaries in between
        //         them. The middle volume cell has two facets which are boundaries. If the middle
        //         volume cell belongs to the 'outside' the two facets stand for two non-connected
        //         boundaries.
        GEO::CUT::plain_facet_set facets = volcell->Facets();
        for ( GEO::CUT::plain_facet_set::const_iterator ifacet=facets.begin(); ifacet!=facets.end(); ++ifacet )
        {
          if ((*ifacet)->OnCutSide())
            numstoredbound += 1;
         }
      }
      // for bisected and touched elements there will only be one facet which is the boundary
      // remark: this is only a shortcut; we could also loop the facets (there is only one!)
      else
        numstoredbound += 1;

    }
  }

  return numstoredbound;
}


/*------------------------------------------------------------------------------------------------*
 | compute GfuncValue for hexahedral integration cell                              rasthofer 04/10|
 *------------------------------------------------------------------------------------------------*/
bool COMBUST::FlameFront::GetIntCellDomainInElementAtCenter(
    const LINALG::SerialDenseMatrix        IntCellCoord,
    const std::vector<double>&             gfuncvalues_ele,
    const DRT::Element::DiscretizationType xfem_distype,
    const DRT::Element::DiscretizationType distype_cell
)
{
  /*------------------------------------------------------------------------------
   * - compute phi at the midpoint of domain integration cell
   * -> >= 0 in plus  'inGplus'==true
   * -> < 0  in minus 'inGplus'==false
   * ----------------------------------------------------------------------------*/

  bool inGplus = false;
  Epetra_SerialDenseVector  cellcenter(3);

//  int numcellnodes = 0;
//  switch(distype_cell)
//  {
//  //    case DRT::Element::tet4:
//  //    {
//  //      numcellnodes = DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::tet4>::numNodePerElement;
//  //      break;
//  //    }
//  case DRT::Element::hex8:
//  {
//    numcellnodes = DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::hex8>::numNodePerElement;
//    break;
//  }
//  case DRT::Element::hex20: {
//    numcellnodes = DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::hex20>::numNodePerElement;
//    break;
//  }
//  default:
//    dserror("Discretization Type (IntCell) not supported yet!");
//  }

  int numelenodes = 0;
  if (xfem_distype==DRT::Element::hex8) {
    numelenodes = DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::hex8>::numNodePerElement;
  }
  else if (xfem_distype==DRT::Element::hex20) {
    numelenodes = DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::hex20>::numNodePerElement;
  }
  else {
    dserror("Discretization Type (Ele) not supported yet!");
  }

  //compute cell center for hex8 element
  //cell center is midpoint of diagonal of integration cell
  for (int i=0; i<3; i++)
  {
    cellcenter(i) = 0.5 * (IntCellCoord(i,6) - IntCellCoord(i,0)) + IntCellCoord(i,0);
  }

  //compute g-funct value at cell center
  double gvalcellcenter = 0.0;
  //  std::vector<double> gvalcellnodes (numcellnodes);
  //  for (int icellnode=0; icellnode<numcellnodes; icellnode++)
  //    gvalcellnodes[icellnode] = 0.0;

  Epetra_SerialDenseVector  funct(numelenodes);
  DRT::UTILS::shape_function_3D(funct,cellcenter(0),cellcenter(1),cellcenter(2),DRT::Element::hex8);
  for (int ielenode=0; ielenode<numelenodes; ielenode++)
  {
    gvalcellcenter += funct(ielenode) * gfuncvalues_ele[ielenode];
  }

  if (gvalcellcenter >= 0.0)
    inGplus = true;

  return inGplus;
}


/*------------------------------------------------------------------------------------------------*
 | compute average GfuncValue for integration cell                                                |
 *------------------------------------------------------------------------------------------------*/
bool COMBUST::FlameFront::GetIntCellDomain(
    const LINALG::SerialDenseMatrix        IntCellCoord,
    const std::vector<double>&             gfuncvalues_cell,
    const DRT::Element::DiscretizationType distype_cell
)
{
  /*------------------------------------------------------------------------------
   * - compute phi at each node of domain integration cell
   * - compute average G-value for each domain integration cell
   * -> >=0 in plus == true
   * _> <0  in minus == false
   * ----------------------------------------------------------------------------*/

  bool inGplus = false;

  int numcellnodes = 0;
  switch(distype_cell)
  {
  case DRT::Element::tet4:
  {
    numcellnodes = DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::tet4>::numNodePerElement;
    break;
  }
  case DRT::Element::hex8:
  {
    numcellnodes = DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::hex8>::numNodePerElement;
    break;
  }
  case DRT::Element::hex20:
  {
    numcellnodes = DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::hex20>::numNodePerElement;
    break;
  }
  default:
  {
    dserror("Discretization Type (IntCell) not supported yet!");
    break;
  }
  }

  //calculate average Gfunc value
  double averageGvalue = 0.0;

  for (int icellnode=0; icellnode<numcellnodes; icellnode++)
    averageGvalue += gfuncvalues_cell[icellnode];

  if (numcellnodes == 0.0)
    dserror("division by zero: number of vertices per cell is 0!");
  averageGvalue /= numcellnodes;

  // determine DomainIntCell position
  // ">=" because of the approximation of the interface,
  // also averageGvalue=0 is possible; as always approximation errors are made
  if(averageGvalue>=0.0)
    inGplus = true;
  //if(averageGvalue==0)
  //  dserror("can't determine DomainIntCell position");

  return inGplus;
}


/*------------------------------------------------------------------------------------------------*
 | compute average GfuncValue for integration cell                                                |
 *------------------------------------------------------------------------------------------------*/
bool COMBUST::FlameFront::GetIntCellDomainInElement(
    const LINALG::SerialDenseMatrix        IntCellCoord,
    const std::vector<double>&             gfuncvalues_ele,
    const DRT::Element::DiscretizationType xfem_distype,
    const DRT::Element::DiscretizationType distype_cell
)
{
  /*------------------------------------------------------------------------------
   * - compute phi at each node of domain integration cell
   * - compute average G-value for each domain integration cell
   * -> >=0 in plus == true
   * _> <0  in minus == false
   * ----------------------------------------------------------------------------*/

  bool inGplus = false;

  int numcellnodes = 0;
  switch(distype_cell)
  {
  case DRT::Element::tet4:
  {
    numcellnodes = DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::tet4>::numNodePerElement;
    break;
  }
  case DRT::Element::hex8:
  {
    numcellnodes = DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::hex8>::numNodePerElement;
    break;
  }
  case DRT::Element::hex20:
  {
    numcellnodes = DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::hex20>::numNodePerElement;
    break;
  }
  default:
  {
    dserror("Discretization Type (IntCell) not supported yet!");
    break;
  }
  }

  int numelenodes = 0;
  if (xfem_distype==DRT::Element::hex8) {
    numelenodes = DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::hex8>::numNodePerElement;
  }
  else if (xfem_distype==DRT::Element::hex20) {
    numelenodes = DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::hex20>::numNodePerElement;
  }
  else {
    dserror("Discretization Type (Ele) not supported yet!");
  }

  std::vector<double> gvalcellnodes (numcellnodes);
  for (int icellnode=0; icellnode<numcellnodes; icellnode++)
    gvalcellnodes[icellnode] = 0.0;

  for (int icellnode=0; icellnode<numcellnodes; icellnode++)
  {
    //calculate shape function at IntCell node
    Epetra_SerialDenseVector  funct(numelenodes);
    DRT::UTILS::shape_function_3D(funct,IntCellCoord(0,icellnode),IntCellCoord(1,icellnode),IntCellCoord(2,icellnode),DRT::Element::hex8);
    /*
     *         ___
     * gval(x)=\   N(x)*gvali
     *         /
     *         ___
     */
    for (int ielenode=0; ielenode<numelenodes; ielenode++)
    {
      gvalcellnodes[icellnode] += funct(ielenode) * gfuncvalues_ele[ielenode];
    }
  }

  //calculate average Gfunc value
  double averageGvalue = 0.0;

  for (int icellnode=0; icellnode<numcellnodes; icellnode++)
    averageGvalue += gvalcellnodes[icellnode];

  if (numcellnodes == 0.0)
    dserror("division by zero: number of vertices per cell is 0!");
  averageGvalue /= numcellnodes;

  //determine DomainIntCell position
  //URSULA 030809 ----- Modfication
  // >= ; because of the approximation of the interface,
  // also averageGvalue=0 is possible; as always approximation errors are made,
  // this is just one more
  if(averageGvalue>=0.0)
    inGplus = true;
  //  if(averageGvalue==0)
  //    dserror("can't determine DomainIntCell position");

  return inGplus;
}


#ifdef PARALLEL
/*------------------------------------------------------------------------------------------------*
 | export flame front to all processors                                               henke 12/09 |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::FlameFront::ExportFlameFront(std::map<int, GEO::BoundaryIntCells>& myflamefront)
{
  //-------------------------------
  // prepare parallel communication
  //-------------------------------
  // get communicator, e.g. from fluid discretization
  const Epetra_Comm& comm = fluiddis_->Comm();
  const int myrank = comm.MyPID();
  const int numproc = comm.NumProc();

  int size_one = 1;

  DRT::Exporter exporter(comm);

  // destination proc (the "next" one)
  int dest = myrank+1;
  if(myrank == (numproc-1))
    dest = 0;

  // source proc (the "last" one)
  int source = myrank-1;
  if(myrank == 0)
    source = numproc-1;

  //#ifdef DEBUG
  //  IO::cout << "number of cell groups (cut column elements) on proc " << myrank << " " << myflamefront.size() << IO::endl;
  //  for(std::map<int, GEO::BoundaryIntCells>::const_iterator cellgroup=myflamefront.begin(); cellgroup != myflamefront.end(); ++cellgroup)
  //    {
  //      // put ID of cut element here
  //      if (cellgroup->first == 1274)
  //      {
  //        IO::cout << "output for element 1274 packing" << IO::endl;
  //        const size_t numcells = (cellgroup->second).size();
  //        IO::cout << "proc " << myrank << " number of integration cells: " << numcells << IO::endl;
  //        for (size_t icell=0; icell<numcells; ++icell)
  //        {
  //          GEO::BoundaryIntCell cell = cellgroup->second[icell];
  //          IO::cout << "proc " << myrank << " cell " << icell << " vertexcoord " << cell.CellNodalPosXYZ();
  //        }
  //      }
  //    }
  //#endif

#ifdef DEBUG
  IO::cout << "proc " << myrank << " flame front pieces for " << myflamefront.size() << " elements available before export" << IO::endl;
#endif

  DRT::PackBuffer data;
  COMBUST::FlameFront::packBoundaryIntCells(myflamefront, data);
  data.StartPacking();
  COMBUST::FlameFront::packBoundaryIntCells(myflamefront, data);

  //-----------------------------------------------------------------
  // pack data (my boundary integration cell groups) for initial send
  //-----------------------------------------------------------------
  std::vector<char> dataSend;
  swap( dataSend, data() );

  //-----------------------------------------------
  // send data around in a circle to all processors
  //-----------------------------------------------
  // loop over processors
  for(int num = 0; num < numproc-1; num++)
  {
    std::vector<int> lengthSend(1,0);
    lengthSend[0] = dataSend.size();

#ifdef DEBUG
    IO::cout << "--- sending "<< lengthSend[0] << " bytes: from proc " << myrank << " to proc " << dest << IO::endl;
#endif

    // send length of the data to be received ...
    MPI_Request req_length_data;
    int length_tag = 0;
    exporter.ISend(myrank, dest, &(lengthSend[0]) , size_one, length_tag, req_length_data);
    // ... and receive length
    std::vector<int> lengthRecv(1,0);
    exporter.Receive(source, length_tag, lengthRecv, size_one);
    exporter.Wait(req_length_data);

    // send actual data ...
    int data_tag = 4;
    MPI_Request req_data;
    exporter.ISend(myrank, dest, &(dataSend[0]), lengthSend[0], data_tag, req_data);
    // ... and receive data
    std::vector<char> dataRecv(lengthRecv[0]);
    exporter.ReceiveAny(source, data_tag, dataRecv, lengthRecv[0]);
    exporter.Wait(req_data);

#ifdef DEBUG
    IO::cout << "--- receiving "<< lengthRecv[0] << " bytes: to proc " << myrank << " from proc " << source << IO::endl;
#endif

    //-----------------------------------------------
    // unpack data (boundary integration cell groups)
    //-----------------------------------------------
    std::map<int, GEO::BoundaryIntCells> flamefront_recv;

    COMBUST::FlameFront::unpackBoundaryIntCells(dataRecv, flamefront_recv);

    //#ifdef DEBUG
    //    IO::cout << "proc " << myrank << " receiving "<< lengthRecv[0] << " bytes from proc " << source << IO::endl;
    //    IO::cout << "proc " << myrank << " receiving "<< flamefront_recv.size() << " flame front pieces from proc " << source << IO::endl;
    //
    //    for(std::map<int, GEO::BoundaryIntCells>::const_iterator cellgroup=flamefront_recv.begin(); cellgroup != flamefront_recv.end(); ++cellgroup)
    //    {
    //      // put ID of cut element here
    //      if (cellgroup->first == 1274)
    //      {
    //        IO::cout << "output for element 1274 unpacking" << IO::endl;
    //        const size_t numcells = (cellgroup->second).size();
    //        IO::cout << "proc " << myrank << " number of integration cells: " << numcells << IO::endl;
    //        for (size_t icell=0; icell<numcells; ++icell)
    //        {
    //          GEO::BoundaryIntCell cell = cellgroup->second[icell];
    //          IO::cout << "proc " << myrank << " cell " << icell << " vertexcoord " << cell.CellNodalPosXYZ();
    //        }
    //      }
    //    }
    //#endif

    // add group of cells to my flame front map
    // remark: all groups of boundary integration cells (flame front pieces within an element) are collected here
    for (std::map<int, GEO::BoundaryIntCells>::const_iterator cellgroup = flamefront_recv.begin(); cellgroup != flamefront_recv.end(); ++cellgroup)
    {
      myflamefront.insert(*cellgroup);
    }

    // make received data the new 'to be sent' data
    dataSend = dataRecv;

    // processors wait for each other
    comm.Barrier();
  }
#ifdef DEBUG
  IO::cout << "proc " << myrank << " flame front pieces for " << myflamefront.size() << " elements available after export" << IO::endl;
#endif
}
#endif


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void COMBUST::FlameFront::packBoundaryIntCells(
    const std::map<int, GEO::BoundaryIntCells>& intcellmap,
    DRT::PackBuffer&                            dataSend)
{
  // pack data on all processors
  // loop entries of map (groups of boundary integration cells)
  for(std::map<int, GEO::BoundaryIntCells>::const_iterator cellgroup=intcellmap.begin(); cellgroup != intcellmap.end(); ++cellgroup)
  {

    // pack data of all boundary integrations cells belonging to an element
    const int elegid = cellgroup->first;
    DRT::ParObject::AddtoPack(dataSend,elegid);

    const int numcells = (cellgroup->second).size();
    DRT::ParObject::AddtoPack(dataSend,numcells);

    for (int icell=0; icell<numcells; ++icell)
    {
      GEO::BoundaryIntCell cell = cellgroup->second[icell];
      // get all member variables from a single boundary integration cell
      const DRT::Element::DiscretizationType distype = cell.Shape();
      DRT::ParObject::AddtoPack(dataSend,distype);

      // coordinates of cell vertices in (fluid) element parameter space
      //      const Epetra_SerialDenseMatrix& vertices_xi = cell.CellNodalPosXiDomain();
      //      const LINALG::SerialDenseMatrix& vertices_xi = cell.CellNodalPosXiDomain();
      const LINALG::SerialDenseMatrix vertices_xi = cell.CellNodalPosXiDomain();
      DRT::ParObject::AddtoPack(dataSend,vertices_xi);

      // coordinates of cell vertices in physical space
      //      const Epetra_SerialDenseMatrix& vertices_xyz = cell.CellNodalPosXYZ();
      //      const LINALG::SerialDenseMatrix& vertices_xyz = cell.CellNodalPosXYZ();
      const LINALG::SerialDenseMatrix vertices_xyz = cell.CellNodalPosXYZ();
      DRT::ParObject::AddtoPack(dataSend,vertices_xyz);
    }
  }
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void COMBUST::FlameFront::unpackBoundaryIntCells(
    const std::vector<char>&                     data,
    std::map<int, GEO::BoundaryIntCells>&   intcellmap)
{
  // pointer to current position in a group of cells in local std::string (counts bytes)
  std::vector<char>::size_type posingroup = 0;

  while (posingroup < data.size())
  {
    // extract fluid element gid
    int elegid = -1;
    DRT::ParObject::ExtractfromPack(posingroup,data,elegid);
    if (elegid < 0) dserror("extraction of element gid failed");

    //extract number of boundary integration cells for this element
    int numvecs = -1;
    DRT::ParObject::ExtractfromPack(posingroup,data,numvecs);

    // vector holding group of boundary integration cells belonging to this element
    GEO::BoundaryIntCells intcellvector;

    for (int icell=0; icell<numvecs; ++icell)
    {
      //--------------------------------------------------------------------
      // extract all member variables for a single boundary integration cell
      //--------------------------------------------------------------------
      // distype of cell
      DRT::Element::DiscretizationType distype;
      int distypeint = -1;
      DRT::ParObject::ExtractfromPack(posingroup,data,distypeint);
      if (distypeint == 4)
        distype = DRT::Element::tri3;
      else if (distypeint == 1)
        distype = DRT::Element::quad4;
      else
        dserror("unexpected distype %d", distypeint);

      LINALG::SerialDenseMatrix vertices_xi;
      DRT::ParObject::ExtractfromPack(posingroup,data,vertices_xi);

      // coordinates of cell vertices in physical space
      LINALG::SerialDenseMatrix vertices_xyz;
      DRT::ParObject::ExtractfromPack(posingroup,data,vertices_xyz);

      //store boundary integration cells in boundaryintcelllist
      intcellvector.push_back(GEO::BoundaryIntCell(distype, -1, vertices_xi,
          Teuchos::null, vertices_xyz));
    }

    // add group of cells for this element to the map
    intcellmap.insert(std::make_pair(elegid,intcellvector));

  }
  // check correct reading
  if (posingroup != data.size())
    dserror("mismatch in size of data %d <-> %d",(int)data.size(),posingroup);
}

/*----------------------------------------------------------------------------------*
 |output to gmsh                                                     rasthofer 05/10|
 *----------------------------------------------------------------------------------*/
void COMBUST::FlameFront::FlameFrontToGmsh(
    const COMBUST::RefinementCell* cell,
    const GEO::BoundaryIntCells&   boundaryintcelllist,
    const GEO::DomainIntCells&     domainintcelllist
)
{
  const bool screen_out = false;

  const std::string filename = IO::GMSH::GetNewFileNameAndDeleteOldFiles("Flamefront", cell->Ele()->Id(), 200, screen_out, 0);
  std::ofstream gmshfilecontent(filename.c_str());
  {
    gmshfilecontent << "View \" " << "Cell \" {\n";
    const int Id = cell->Ele()->Id();
    const int numnode = cell->Ele()->NumNode();
    const DRT::Element::DiscretizationType distype = cell->Ele()->Shape();
    const std::vector<std::vector<double> > vertexcoord = cell->GetVertexCoord();
    LINALG::SerialDenseMatrix Pos(3,numnode);
    for (int i=0; i<numnode; i++)
    {
      for (size_t k=0; k<3; k++)
      {
        Pos(k,i) = vertexcoord[i][k];
      }
    }
    IO::GMSH::cellWithScalarToStream(distype, Id, Pos, gmshfilecontent);
    gmshfilecontent << "};\n";
  }
  {
    gmshfilecontent << "View \" " << "BoundaryIntCells \" {\n";
    for(std::size_t icell = 0; icell < boundaryintcelllist.size(); icell++)
    {
      GEO::BoundaryIntCell cell = boundaryintcelllist[icell];
      const LINALG::SerialDenseMatrix& cellpos = cell.CellNodalPosXiDomain();//cell.CellNodalPosXYZ();
      const double color = 5;
      gmshfilecontent << IO::GMSH::cellWithScalarToString(cell.Shape(), color, cellpos) << std::endl;
    }
    gmshfilecontent << "};\n";
  }
  {
    gmshfilecontent << "View \" " << "DomainIntCells \" {\n";
    for(std::size_t icell = 0; icell < domainintcelllist.size(); icell++)
    {
      GEO::DomainIntCell cell = domainintcelllist[icell];
      const LINALG::SerialDenseMatrix& cellpos = cell.CellNodalPosXiDomain();//cell.CellNodalPosXYZ();

      int domain_id = 0;
      if (cell.getDomainPlus())
        domain_id = 1;
      const double color = domain_id;
      gmshfilecontent << IO::GMSH::cellWithScalarToString(cell.Shape(), color, cellpos) << std::endl;
    }
    gmshfilecontent << "};\n";
  }
  gmshfilecontent.close();

  return;
}

