/*!-----------------------------------------------------------------------------------------------*
 \file combust_flamefront.cpp

 \brief

  detailed description in header file combust_interface.H

<pre>
Maintainer: Florian Henke
            henke@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15265
</pre>
 *------------------------------------------------------------------------------------------------*/
#ifdef CCADISCRET

#include <iterator>

#include "../drt_lib/standardtypes_cpp.H" // has to be declared first

#include "combust_flamefront.H"
#include "combust_defines.H"
#include "../drt_combust/combust3.H"
#include "../drt_combust/combust3_utils.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_exporter.H"
#include "../drt_lib/drt_parobject.H"
#include "../drt_geometry/integrationcell_coordtrafo.H"
#include "../drt_geometry/tetrahedradecomposition.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_fem_general/drt_utils_integration.H"
#include "../drt_io/io_gmsh.H"
#include "../linalg/linalg_utils.H" // LINALG::Export
#include <Teuchos_StandardParameterEntryValidators.hpp>
#include <Teuchos_TimeMonitor.hpp>
#include "../drt_geometry/position_array.H"
#include "../drt_geometry/intersection_service_templates.H"

#ifdef PARALLEL
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

#include "../drt_cut/cut_levelsetintersection.H"
#include "../drt_cut/cut_meshintersection.H"
#include "../drt_cut/cut_elementhandle.H"
#include "../drt_cut/cut_integrationcell.H"


/*------------------------------------------------------------------------------------------------*
 | constructor                                                                        henke 10/08 |
 *------------------------------------------------------------------------------------------------*/
COMBUST::FlameFront::FlameFront(
    const Teuchos::RCP<const DRT::Discretization> fluiddis,
    const Teuchos::RCP<DRT::Discretization> gfuncdis
) :
fluiddis_(fluiddis),
gfuncdis_(gfuncdis),
phinm_(Teuchos::null),
phin_(Teuchos::null),
phinp_(Teuchos::null),
gradphi_(Teuchos::null),
maxRefinementLevel_(0),
xfeminttype_(INPAR::COMBUST::xfemintegration_tetgen)
{
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
void COMBUST::FlameFront::UpdateFlameFront(
    const Teuchos::ParameterList& combustdyn,
    const Teuchos::RCP<const Epetra_Vector>& phin,
    const Teuchos::RCP<const Epetra_Vector>& phinp,
    bool ReinitModifyPhi
)
{
  // rearrange and store phi vectors
  StorePhiVectors(phin, phinp);
  // modify phi vectors nearly zero
  //  ModifyPhiVector(combustdyn, ReinitModifyPhi);

  // generate the interface geometry based on the G-function (level set field)
  // remark: must be called after StorePhiVectors, since it relies on phinp_
  ProcessFlameFront(combustdyn, phinp_);

  // compute smoothed gradient of G-function field
  // remark: must be called after ProcessFlameFront, since it relies on the new interface position

  const INPAR::COMBUST::SmoothGradPhi smoothgradphi = DRT::INPUT::IntegralValue<INPAR::COMBUST::SmoothGradPhi>
  (combustdyn.sublist("COMBUSTION FLUID"),"SMOOTHGRADPHI");
  if(smoothgradphi!=INPAR::COMBUST::smooth_grad_phi_none)
    CallSmoothGradPhi(combustdyn);

  return;
}

/*------------------------------------------------------------------------------------------------*
 | public: generate the interface geometry based on the G-function (level set field)  henke 10/08 |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::FlameFront::ProcessFlameFront(
    const Teuchos::ParameterList& combustdyn,
    const Teuchos::RCP<const Epetra_Vector> phi)
{
  /* This function is accessible from outside the FlameFront class. It can be called e.g. by the
   * interface handle constructor. If the FlameFront class will always do the same thing, this
   * function could be spared. Then the content of this function could be added to the constructor
   * of this class.
   *
   * henke 10/08
   */

  if (fluiddis_->Comm().MyPID()==0)
    std::cout << "\n---  processing flame front ... " << std::flush;;

  const Teuchos::RCP<Epetra_Vector> phicol = rcp(new Epetra_Vector(*fluiddis_->NodeColMap()));
  if (phicol->MyLength() != phi->MyLength())
    dserror("vector phi needs to be distributed according to fluid node column map");

  // maps are cleared and refilled on every call of this function
  myelementintcells_.clear();
  myboundaryintcells_.clear();

  // get type of integrationcells
  xfeminttype_ = DRT::INPUT::IntegralValue<INPAR::COMBUST::XFEMIntegration>(combustdyn.sublist("COMBUSTION FLUID"),"XFEMINTEGRATION");

  // loop over fluid (combustion) column elements
  // remark: loop over row elements would be sufficient, but enrichment is done in column loop
  for (int iele=0; iele<fluiddis_->NumMyColElements(); ++iele)
  {
    // get element from discretization
    const DRT::Element *ele = fluiddis_->lColElement(iele);

#ifdef DEBUG
#ifdef D_FLUID3
    if(ele->ElementType() != DRT::ELEMENTS::Combust3Type::Instance())
      // this is not compulsory, but combust3 elements are expected here!
      dserror("unexpected element type: this should be of combust3 type!");
#endif
#endif

    // create refinement cell from a fluid element -> cell will have same geometry as element!
    const Teuchos::RCP<COMBUST::RefinementCell> rootcell = rcp(new COMBUST::RefinementCell(ele));

    // refinement strategy is turned on
    if (DRT::INPUT::IntegralValue<int>(combustdyn.sublist("COMBUSTION GFUNCTION"),"REFINEMENT") == true)
    {
      //std::cout << "starting refinement ..." << std::endl;
      maxRefinementLevel_ = Teuchos::getIntParameter(combustdyn.sublist("COMBUSTION GFUNCTION"),"REFINEMENTLEVEL");
      //std::cout << "maximal refinement level " << maxRefinementLevel_ << std::endl;
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
    CaptureFlameFront(rootcell);

    // should not be necessary
#if 0
    // delete all refinement cells of root cell
    if (DRT::INPUT::IntegralValue<int>(combustdyn.sublist("COMBUSTION GFUNCTION"),"REFINEMENT") == true)
    {
      rootcell->Clear();
    }
#endif
  }

  //TEST
  //int cells = 0;
  //for (std::map<int, GEO::DomainIntCells>::iterator it=myelementintcells_.begin(); it!=myelementintcells_.end(); it++)
  //{
  //   GEO::DomainIntCells IntCells = it->second;
  //   //std::cout << IntCells.size() << std::endl;
  //   cells = cells + IntCells.size();
  //}
  //std::cout << "Anzahl Integrationszellen: " << cells << std::endl;

  std::cout << "done" << std::endl;
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
  //const Teuchos::RCP<Epetra_Vector> phinmrow = rcp(new Epetra_Vector(*fluiddis_->NodeRowMap()));
  const Teuchos::RCP<Epetra_Vector> phinrow  = rcp(new Epetra_Vector(*fluiddis_->NodeRowMap()));
  const Teuchos::RCP<Epetra_Vector> phinprow = rcp(new Epetra_Vector(*fluiddis_->NodeRowMap()));

  //if (phinmrow->MyLength() != phinm->MyLength())
  //  dserror("vectors phinmrow and phinm must have the same length");
  if (phinrow->MyLength() != phin->MyLength())
    dserror("vectors phinrow and phin must have the same length");
  if (phinprow->MyLength() != phinp->MyLength())
    dserror("vectors phinprow and phinp must have the same length");

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
     * henke 02/09
     */

    // get the G-function values corresponding to fluid nodes
    //*phinmrow = *phinm;
    *phinrow  = *phin;
    *phinprow = *phinp;
    /* ausgeschrieben, was *phinprow = *phinp in einer Zeile macht:
       for(int lnodeid=0;lnodeid<fluiddis_->NumMyRowNodes();lnodeid++)
       {
         // get the G-function value corresponding to this fluid node
         double value = (*phinp)[lnodeid];
         phinprow->ReplaceMyValue(lnodeid,value);
       }
     */
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
  const Teuchos::RCP<Epetra_Vector> phincol = rcp(new Epetra_Vector(*fluiddis_->NodeColMap()));
  LINALG::Export(*phinrow,*phincol);
  const Teuchos::RCP<Epetra_Vector> phinpcol = rcp(new Epetra_Vector(*fluiddis_->NodeColMap()));
  LINALG::Export(*phinprow,*phinpcol);

  // store vector on fluiddis NodeColMap holding G-function values in member variable
  phin_  = phincol;
  phinp_ = phinpcol;
}

/*------------------------------------------------------------------------------------------------*
 | modify the fluid phi-vectors phinp_ and phin_ but not the originally scatra phi vector         |
 | we set phi-values near to zero (abs(*)< TOL) to 0.0 to avoid small integration cells           |
 | and ill conditioned matrices -> we need this implicit in xdofmapcreation_combust               |
 |                                                                        schott/winklmaier 08/10 |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::FlameFront::ModifyPhiVector(const Teuchos::ParameterList& combustdyn, bool ReinitModifyPhi)
{
  std::cout << "\n---  Modify the fluid phi-vector (G-function) at nodes with small values ... " << std::endl << std::flush ;

  if (ReinitModifyPhi==true)
  {
    // Benedikt:
    // this case shall circumvent the tetgen-problem only!
    // when tetgen is not used any more, this case can be removed

    // get relative tolerance for which we modify G-function values in the fluid part
    const double ModifyTOL = 1e-012;

    cout << "\n \t -> ModifyTOL for G-func values is set to " << ModifyTOL << " (relative to element diameter) to handle problems with tetgen!" << endl;

    // count modified phi node values
    int modified_counter = 0;


    //========================================================================================
    // get an absolute tolerance for which we modify phi values with |phi(node)-0.0| < TOL
    //========================================================================================

    // number space dimensions for 3d combustion element
    const size_t nsd = 3;

    // TODO template this part with switch(DISTYPE)
    const DRT::Element::DiscretizationType DISTYPE = DRT::Element::hex8;

    // number of nodes of this element for interpolation!!!
    const size_t numnode = DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;

    // pointer to phivector at time n with rowmap
    RCP<Epetra_Vector> phinTmp = rcp(new Epetra_Vector(*fluiddis_->NodeRowMap()));
    LINALG::Export(*phin_,*phinTmp);
    // pointer to phivector at time n+1 with rowmap
    RCP<Epetra_Vector> phinpTmp = rcp(new Epetra_Vector(*fluiddis_->NodeRowMap()));
    LINALG::Export(*phinp_,*phinpTmp);


    // phinp_ and phin_ are given in same nodeColMap
    // loop over nodes

    for(int inode=0;inode<fluiddis_->NumMyRowNodes();inode++)
    {
      // get node from fluid discretization
      const DRT::Node *actnode = fluiddis_->lRowNode(inode);

      // this processor has to modify the phi-value for the current node
      // get local processor id of node
      int nodeID = actnode->Id();

      // get all adjacent elements of this row node -> we find these in a
      int numberOfElements = actnode->NumElement();
      const DRT::Element* const* elements = actnode->Elements();

      // initialize maximal element diameter to 0.0
      double maxEleDiam = 0.0;

      //==========================================================================
      // loop over adjacent elements => set max eleDiam
      for(int ele_current=0; ele_current<numberOfElements; ele_current++)
      {
        // get current element
        const DRT::Element* ele_adj = elements[ele_current];

        //-------------------------------------------------------------
        // get element diameter for the current adjacent element
        //-------------------------------------------------------------
        if(ele_adj->Shape() != DRT::Element::hex8) dserror("ModifyPhiVector only implemented for hex8 elements!");

        static LINALG::Matrix<nsd,numnode> xyze_adj;
        GEO::fillInitialPositionArray<DISTYPE>(ele_adj, xyze_adj);

        // calculate element diameter
        const double hk_current = COMBUST::getEleDiameter<DISTYPE>(xyze_adj);


        //-------------------------------------------------------------
        // update maxEleDiam
        //-------------------------------------------------------------
        maxEleDiam = max(hk_current,maxEleDiam);

      } // end loop over adjacent
      //==========================================================================

      // set TOL dependent on maxEleDiam
      const double TOL = maxEleDiam * ModifyTOL;


      // get current phi-value
      // lid is the local processor id of the current node
      // phinp_ is a nodeColMap
      const int lid = (*phinpTmp).Map().LID(nodeID);
      if (lid<0)
        dserror("Proc %d: Cannot find gid=%d in Epetra_Vector",(*phinpTmp).Comm().MyPID(),lid);
      double phinp_current = (*phinpTmp)[lid];
      double phin_current = (*phinTmp)[lid];

      // decide if modification to phi=0.0 is necessary
      // reset small phi values to zero
      if (fabs(phinp_current - 0.0) < TOL)
      {
        cout << "\t !!!warning (ModifyPhiVector)!!! we modify a small phi-value to 0.0" << std::endl;
        (*phinpTmp)[lid] = 0.0;
        modified_counter++;
      }
      if (fabs(phin_current - 0.0) < TOL)
      {
        (*phinTmp)[lid] = 0.0;
        // no additional cout-comment here because this value has been
        // modified in the timestep before and cout was created there
      }

    } // end loop over nodes

    if(modified_counter > 0)
      std::cout << "---  \t number of modified phi-values: " << modified_counter << " ...done" << std::flush;;

    LINALG::Export(*phinTmp,*phin_);
    LINALG::Export(*phinpTmp,*phinp_);
  }
  else // non-ReinitModifyPhi case (standard case)
  {
    //========================================================================================
    // get a tolerance for the main modify criterion: ModifyTOL
    // REMARK:
    // crit_i in [0,1] is a critical value for node i_
    //		crit_i near 0.0 means: "node i is !!!very critical!!! for the adjacent element"
    //		crit_i near 1.0 means: "node i is !!!uncritical!!! for the adjacent element"
    // crit_i := crit(i,1) * ... * crit(i,8) in [0,1] for hex8 elements with 8 nodes
    //           crit(i,j):= critical value for pair of nodes i and j
    // if (crit_i := crit(i,1) * ... * crit(i,8) < ModifyTOL) -> modify the critical node i
    //========================================================================================
    // get relative tolerance for which we modify G-function values ( done just in the fluid part )
    const double ModifyTOL = combustdyn.sublist("COMBUSTION FLUID").get<double>("PHI_MODIFY_TOL");

    //========================================================================================
    // get a tolerance for which we look for nodes which potentially get modified
    // REMARK:
    // SearchTOL means e.g. |phi_i-0.0| < SearchTOL * max(eleDiam, over adjacent elements to node i)
    // thats not a modify criterion!
    // usually we choose SearchTOL = 0.1
    //========================================================================================
    const double SearchTOL = 0.1;
    cout << "\n \t -> SearchTOL for G-func values which get potentially modified is set to " << SearchTOL << " !" << endl;

    // tolerance for which we interpret a value as zero-> these phi-values are set to zero!
    const double TOL_zero = 1e-14;

    cout << "\n \t -> TOL_zero for G-func values which are a numerical null is set to " << TOL_zero << " !" << endl;


    // number space dimensions for 3d combustion element
    const size_t nsd = 3;

    // TODO template this part with switch(DISTYPE)
    const DRT::Element::DiscretizationType DISTYPE = DRT::Element::hex8;

    // number of nodes of this element for interpolation!!!
    const size_t numnode = DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;
    // count modified phi node values
    int modified_counter = 0;


    // pointer to phivector at time n with rowmap
    RCP<Epetra_Vector> phinTmp = rcp(new Epetra_Vector(*fluiddis_->NodeRowMap()));
    LINALG::Export(*phin_,*phinTmp);
    // pointer to phivector at time n+1 with rowmap
    RCP<Epetra_Vector> phinpTmp = rcp(new Epetra_Vector(*fluiddis_->NodeRowMap()));
    LINALG::Export(*phinp_,*phinpTmp);

    // phinp_ and phin_ are given in same nodeColMap
    // loop over nodes

    for(int inode=0;inode<fluiddis_->NumMyRowNodes();inode++)
    {
      // get node from fluid discretization
      const DRT::Node *actnode = fluiddis_->lRowNode(inode);


      // this processor has to modify the phi-value for the current node
      // get local processor id of node
      int nodeID = actnode->Id();


      // get all adjacent elements of this row node
      int numberOfElements = actnode->NumElement();
      const DRT::Element* const* elements = actnode->Elements();


      // element critical value current node
      // initialization with 1.0
      // REMARK: 1.0 means uncritical, 0.0 means highly critical
      Epetra_SerialDenseVector ele_critical_inode_np(numberOfElements);
      //ele_critical_inode_np.Clear();
      Epetra_SerialDenseVector ele_critical_inode_n(numberOfElements);
      //ele_critical_inode_n.Clear();

      const int lid = (*phinpTmp).Map().LID(nodeID);
      if (lid<0)
        dserror("Proc %d: Cannot find gid=%d in Epetra_Vector",(*phinpTmp).Comm().MyPID(),lid);
      double phinp_current = (*phinpTmp)[lid];
      double phin_current = (*phinTmp)[lid];

      // initialize maximal element diameter to 0.0
      double maxEleDiam = 0.0;

      //==========================================================================
      // loop over adjacent elements => set max eleDiam
      for(int ele_current=0; ele_current<numberOfElements; ele_current++)
      {
        // get current element
        const DRT::Element* ele_adj = elements[ele_current];

        //-------------------------------------------------------------
        // get element diameter for the current adjacent element
        //-------------------------------------------------------------
        if(ele_adj->Shape() != DRT::Element::hex8) dserror("ModifyPhiVector only implemented for hex8 elements!");

        static LINALG::Matrix<nsd,numnode> xyze_adj;
        GEO::fillInitialPositionArray<DISTYPE>(ele_adj, xyze_adj);

        //cout << xyze_adj << endl;

        // calculate element diameter
        const double hk_current = COMBUST::getEleDiameter<DISTYPE>(xyze_adj);

        //cout << "ele-diameter for element" << ele_adj->Id() << " is " << hk_current << endl;

        //-------------------------------------------------------------
        // update maxEleDiam
        //-------------------------------------------------------------
        maxEleDiam = max(hk_current,maxEleDiam);

      } // end loop over adjacent
      //==========================================================================
      // decide if the current node "maybe" have to be modified
      // call the main modification check if search-criterion is positive
      //   and phi-value is not 0.0 already
      bool modify_phinp_pot = false;
      bool modify_phin_pot  = false;

      if ( (fabs(phinp_current - 0.0) < SearchTOL * maxEleDiam) and (phinp_current != 0.0) )
      {
        modify_phinp_pot = true;
        cout << "\t --- node "<< nodeID << ": nodal phinp value get potentially modified" << std::endl<< std::flush;
      }

      if ( (fabs(phin_current - 0.0) < SearchTOL * maxEleDiam)  and (phin_current  != 0.0) )
        modify_phin_pot  = true;
      //==========================================================================

      // we decide if we have to show for a need of a modification of the phi-value to zero
      // do nothing if phi-value is zero already or phi-value is large enough
      if( (modify_phinp_pot == true) or (modify_phin_pot == true) )
      {
        // boolian for real modification of phi value
        bool inode_modified_n = false;
        bool inode_modified_np = false;

        // modify numerical zeros to 0.0 directly
        if((fabs(phinp_current - 0.0) < TOL_zero) && (phinp_current!= 0.0))
        {
          inode_modified_np = true;
          cout << "\t\t -> a numerical zero is set to 0.0" << std::endl<< std::flush;
        }
        if((fabs(phin_current - 0.0) < TOL_zero) && (phin_current!= 0.0))  inode_modified_n = true;
        //==========================================================================
        if((inode_modified_np == false) or (inode_modified_n == false))
        {
          // calculate critical node values for all adjacent elements
          for(int ele_current=0; ele_current<numberOfElements; ele_current++)
          {
            // get current element
            const DRT::Element* ele_adj = elements[ele_current];

            if(ele_adj->Shape() != DRT::Element::hex8) dserror("ModifyPhiVector only implemented for hex8 elements!");

            // initialize entry to 1.0
            ele_critical_inode_np(ele_current) = 1.0;
            ele_critical_inode_n(ele_current) = 1.0;

            // loop over all nodes of the current element
            const int* ptToNodeIds_adj = ele_adj->NodeIds();

            // get phi-values of current adjacent element ele_adj
            // create vector "ephinp" holding scalar phi values for this element
            Epetra_SerialDenseVector ephinp(numnode); //local vector phi-values of adjacent element
            Epetra_SerialDenseVector ephin(numnode); //local vector phi-values of adjacent element

            // which node in param space of element ele_adj has actnode
            int ID_param_space = -1;

            // get vector of node GIDs of this adjacent element -> needed for ExtractMyValues
            vector<int> nodeID_adj(numnode);
            for (size_t node=0; node < numnode; node++){
              nodeID_adj[node] = ptToNodeIds_adj[node];
              // get local number of node actnode in ele_adj
              if(actnode->Id() == ptToNodeIds_adj[node]) ID_param_space = node;
            }

            // extract the phi-values of adjacent element with local ids from global vector *phinp
            // get pointer to vector holding G-function values at the fluid nodes
            DRT::UTILS::ExtractMyValues(*phinp_, ephinp, nodeID_adj);
            LINALG::Matrix<numnode,1> ephinp_adj(ephinp);

            DRT::UTILS::ExtractMyValues(*phin_, ephin, nodeID_adj);
            LINALG::Matrix<numnode,1> ephin_adj(ephin);

            //==============================================================================================
            // ------calculate critical values crit_adj_inode(j) for all nodes jnode to current node inode------
            //==============================================================================================
            // REMARK:
            // we define for each node j a critical value to node i in [0,1] if there is a sign change in the
            // corresponding phi-values
            // if the zero point at the line between inode and jnode is near node i, the critical value
            // gets a value near 0 => such zero points yield ill-conditioned entries in matrices
            // with enriched shape functions
            // if the zero point at the line is near node j, the critical value
            // gets a value near 1 => such zero points are not ill for the current node i, but maybe for node j,
            // but this checks the loop iteration for jnode
            //
            // we have to distinguish three types of lines
            // 1. lines between neighbouring nodes -> lines are edges of the element (linear progress of crit_adj_val)
            // 2. lines between nodes within one 2D-face -> lines lie in a 2D face (quadratic progress of crit_adj_val)
            // 3. lines between nodes within the whole 3D-element (3D-diagonals) -> (cubic progress)
            // => we simplify the three cases to case 1.
            // illustration: the intersection points of the resulting interface (boundary integration cells) are connected
            //               linearly. Therefore a quadratic or cubic progress is not necessary.
            //               A estimated linear progress overestimates the quadratic of cubic progress,
            //               so we get a !more critical! value for a linear progress, that's okay!

            // loop over nodes: calculate critical values for each node jnode!=inode to current inode
            LINALG::Matrix<numnode,1> crit_adj_inode_n;
            LINALG::Matrix<numnode,1> crit_adj_inode_np;
            for(int jnode=0;jnode<(int)numnode;jnode++)
            {
              if(jnode != ID_param_space)
              {
                // we take the same formula for all node pairs (inode,jnode), inode=fix, jnode =1..numnode
                // REMARK: see illustration above

                // get phi-values for the to nodes: jnode and inode
                // it holds for a linear progress
                // dx_i denotes the distance between node i and the existing intersection point
                // at the line between inode and jnode holds
                // |dx_i| / (|dx_i| + |dx_j|) = |phi(inode)| / (|phi(inode)| + |phi(jnode)|)
                //
                // we set:
                // crit_adj_inode = max (0, -sign(phi_i)*sign(phi_j) * |phi(inode)| / (|phi(inode)| + |phi(jnode)|) )
                if(modify_phinp_pot == true)
                {
                  double phi_np_i = ephinp_adj(ID_param_space);
                  double phi_np_j = ephinp_adj(jnode);



                  int sign_phi_np_i = (phi_np_i > 0) - (phi_np_i < 0);
                  int sign_phi_np_j = (phi_np_j > 0) - (phi_np_j < 0);


                  // standard cases:
                  // phi_j > 0 and phi_i <~ 0
                  // phi_j < 0 and phi_i >~ 0
                  double relative_distance_np =  fabs(phi_np_i) / (fabs(phi_np_i) + (fabs(phi_np_j)) );


                  if(sign_phi_np_i != sign_phi_np_j)
                  {
                    crit_adj_inode_np(jnode) = relative_distance_np;
                  }
                  else if((sign_phi_np_i == sign_phi_np_j) && (sign_phi_np_i != 0))
                  {
                    // no intersection point between
                    crit_adj_inode_np(jnode) = 1.0;
                  }
                  else
                  {
                    dserror("impossible");
                  }
                } // end if inode_modified_np
                else
                {
                  // set critical value to 1.0
                  crit_adj_inode_np(jnode) = 1.0;
                }
                if(modify_phin_pot ==true)
                {
                  double phi_n_i = ephin_adj(ID_param_space);
                  double phi_n_j = ephin_adj(jnode);

                  int sign_phi_n_i = (phi_n_i > 0) - (phi_n_i < 0);
                  int sign_phi_n_j = (phi_n_j > 0) - (phi_n_j < 0);


                  double relative_distance_n  =  fabs(phi_n_i)  / (fabs(phi_n_i)  + (fabs(phi_n_j))  );
                  if(sign_phi_n_i != sign_phi_n_j)
                  {
                    crit_adj_inode_n(jnode)  = relative_distance_n;
                  }
                  else if((sign_phi_n_i == sign_phi_n_j) && (sign_phi_n_i != 0))
                  {
                    // no intersection point between
                    crit_adj_inode_n(jnode) = 1.0;
                  }
                  else
                  {
                    dserror("impossible");
                  }
                } // end if inode_modified_n
                else
                {
                  // set critical value to 1.0
                  crit_adj_inode_n(jnode)  = 1.0;
                }
              }
              else
              {
                // set critical value for the node inode itself
                crit_adj_inode_np(jnode) = 1.0;
                crit_adj_inode_n(jnode)  = 1.0;
              }
              // update the element critical value for node inode
              ele_critical_inode_np(ele_current) *= crit_adj_inode_np(jnode);
              ele_critical_inode_n(ele_current)  *= crit_adj_inode_n(jnode);
            }// end loop over nodes of current adjacent element
          } // end loop over adjacent

          // ============== decide if the corresponding phi value of node inode gets modified ===========

          double minimal_crit_val_np = 1.0;
          double minimal_crit_val_n = 1.0;

          for(int ele = 0; ele < numberOfElements; ele++)
          {
            if (fabs(ele_critical_inode_np(ele) - 0.0 ) < ModifyTOL)
            {
              inode_modified_np = true;
            }
            if (fabs(ele_critical_inode_n(ele)  - 0.0 ) < ModifyTOL)
            {
              inode_modified_n  = true;
            }
            minimal_crit_val_np = min(minimal_crit_val_np, fabs(ele_critical_inode_np(ele) - 0.0 ));
            minimal_crit_val_n = min(minimal_crit_val_n, fabs(ele_critical_inode_n(ele) - 0.0 ));
          }
          if(inode_modified_np == true)
            cout << "\t \t    minimal_crit_val_np was: " << minimal_crit_val_np << std::endl<< std::flush;
          if(inode_modified_n == true)
            cout << "\t \t    minimal_crit_val_n was: " << minimal_crit_val_n << std::endl<< std::flush;
        } // end set inode_modifed_np status
        //==========================================================================



        // decide if modification to phi=0.0 is necessary
        // reset small phi values to zero
        if (inode_modified_np == true)
        {
          cout << "\t \t -> modified: a critical phi-value was modified to 0.0" << std::endl<< std::flush;
          (*phinpTmp)[lid] = 0.0;
          modified_counter++;
        }

        if (inode_modified_n == true)
        {
          (*phinTmp)[lid] = 0.0;
          // no additional cout-comment here because this value has been
          // modified in the timestep before and cout was created there
        }
      } // end if SearchTOL

      // REMARK: all not row nodes are reconstructed by another processor!!!
    } // end loop over nodes


    if(modified_counter > 0)
      std::cout << "---  \t number of modified phi-values: " << modified_counter << " ...done\n" << std::flush;

    LINALG::Export(*phinTmp,*phin_);
    LINALG::Export(*phinpTmp,*phinp_);

  }

  return;
}


/*------------------------------------------------------------------------------------------------*
 | compute smoothed gradient field of G-function (level set field) for use in fluid time          |
 | integration scheme                                                                schott 06/10 |
 *------------------------------------------------------------------------------------------------*/
template <const size_t nsd_real>
void COMBUST::FlameFront::ComputeSmoothGradPhi(const Teuchos::ParameterList& combustdyn)
{
  // get type of reconstruction
  const INPAR::COMBUST::SmoothGradPhi SmoothGradPhi = DRT::INPUT::IntegralValue<INPAR::COMBUST::SmoothGradPhi>
  (combustdyn.sublist("COMBUSTION FLUID"),"SMOOTHGRADPHI");

  // number space dimensions for 3d combustion element
  const size_t nsd = 3;

  // TODO template the 2nd part with switch(DISTYPE)
  const DRT::Element::DiscretizationType DISTYPE = DRT::Element::hex8;

  // gradphi_ gets a NodeColMap
  gradphi_ = rcp(new Epetra_MultiVector(*fluiddis_->NodeColMap(),nsd));

  // before we can export to NodeColMap we need reconstruction with a NodeRowMap
  const Teuchos::RCP<Epetra_MultiVector> gradphirow = rcp(new Epetra_MultiVector(*fluiddis_->NodeRowMap(),nsd));

  // map of pointers to nodes which must be reconstructed by this processor
  // key is the local id at processor
  std::map<int, const DRT::Node*> nodesToReconstruct;
  nodesToReconstruct.clear();

  //============================================================
  // I.PART:  loop over all elements in column!!! map -> intersected elements at processor boundary
  // must be seen by all processors!!!
  // => find all nodes which must be reconstructed
  //============================================================
  for(int iele=0;iele<fluiddis_->NumMyColElements();iele++)
  {
    // get element from fluid discretization
    const DRT::Element *actele = fluiddis_->lColElement(iele);

    //--------------------------------------------------------------------------------------------
    // find out whether this element is intersected or not by looking at number of integration cells
    //--------------------------------------------------------------------------------------------

    // get global Id of element
    int actele_GID = actele->Id();
    // get vector of all integration cells of this element
    GEO::DomainIntCells IntCells = myelementintcells_[actele_GID];
    // get vector of all boundary integration cells of this element -> also touched elements need reconstruction
    GEO::BoundaryIntCells BoundIntCells = myboundaryintcells_[actele_GID];
    //    cout << "actele->ID()" << actele->Id() << endl;
    //    cout << actele->Id() << "\t" << IntCells.size() << endl;
    //    cout << actele->Id() << "\t" <<BoundIntCells.size() << endl;

    //    if (IntCells.size() > 1 || BoundIntCells.size() > 0) // element is intersected or touched
    //    {
    //cout << "element with ID: " << actele->Id() << "is reconstructed "<< endl;
    // -> assemble all adjacent nodes, these node values must be reconstructed

    // number of nodes of this element (number of vertices)
    // e.g. hex20 elements hasGEO::DomainIntCells IntCells = myelementintcells_[actele_GID]; numnode=20 but numberOfNodes=8
    // if actele is an element at processor boundary, actele will be less than 8
    // TODO:: check this
    const int numberOfNodes = actele->NumNode();

    // get vector of pointers of node (for this element)
    const DRT::Node* const* ele_vecOfPtsToNode = actele->Nodes();

    for(int vec_it = 0; vec_it < numberOfNodes; vec_it++)
    {
      // get owner of the node to compare with my_rank
      int node_owner = (ele_vecOfPtsToNode[vec_it])->Owner();

      // check wheather this node is a row node, compare with actual processor id
      if(node_owner == fluiddis_->Comm().MyPID() )
      {
        // get local processor id of node

        // TODO: check whether insert in map avoids twice-inserting
        int lid = ele_vecOfPtsToNode[vec_it]->LID();
        nodesToReconstruct[lid] = ele_vecOfPtsToNode[vec_it];
      } // end if

      // REMARK: all not row nodes are reconstructed by another processor!!!
    }// end vec iteration (vec of nodes of element)
    //    }// end if element intersected
  }// end for loop over row elements

  typedef std::map<int, const DRT::Node*> Node_Map;

  //===============================================================================================
  // II.PART: reconstruct nodes
  // => row nodes!!! adjacent to intersected elements must be reconstructed by this processor
  //===============================================================================================

  // number of nodes of this element for interpolation!!!
  const size_t numnode = DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;


  //===============loop over all nodes inserted in nodesToReconstruct==============
  // map iterator
  typedef std::map<int, const DRT::Node*> Map_NodesToReconstruct;
  typedef Map_NodesToReconstruct::iterator Reconstruct_iterator;

  for(Reconstruct_iterator it_node = nodesToReconstruct.begin(); it_node!= nodesToReconstruct.end(); it_node++)
  {
    // define vector for smoothed values (Phi, grad_Phi) of current node
    static LINALG::Matrix<nsd_real,1> PHI_SMOOTHED; // for real reconstruction 2D or 3D
    PHI_SMOOTHED.Clear();
    static LINALG::Matrix<nsd,1> PHI_SMOOTHED_3D; // set whole 3D vector also for 2D examples
    PHI_SMOOTHED_3D.Clear();


    // get local processor id of current node and pointer to current node
    //int lid_node = it_node->first;
    const DRT::Node* ptToNode = it_node->second;

    int numberOfElements = ptToNode->NumElement();


    // =================SET TYPE OF RECONSTRUCTION FOR CURRENT NODE =======================================
    // REMARK:
    // although least squares reconstruction is chosen in input file, boundary elements have to be reconstructed
    // by the mean value (average) method because there are not enough adjacent elements to reconstruct
    // the nodal phi gradient
    // minimal number of adjacent elements for least_squares reconstruction is 3 for 3D and 2 for 2D, so nsd_real
    // mean value method is the same for 2D and 3D

    // initialize reconstruction type for current node with none
    INPAR::COMBUST::SmoothGradPhi CurrTypeOfNodeReconst = INPAR::COMBUST::smooth_grad_phi_none;

    // set reconstruction node for current node with special cases (boundary elements)
    if (SmoothGradPhi == INPAR::COMBUST::smooth_grad_phi_none)
      dserror("SmoothGradPhi() should not be called for reconstruction type Grad_phi_None");
    else if (SmoothGradPhi == INPAR::COMBUST::smooth_grad_phi_meanvalue)
      CurrTypeOfNodeReconst = INPAR::COMBUST::smooth_grad_phi_meanvalue;
    else if ( (SmoothGradPhi == INPAR::COMBUST::smooth_grad_phi_leastsquares_3D  && numberOfElements < int(nsd_real)) ||
        (SmoothGradPhi == INPAR::COMBUST::smooth_grad_phi_leastsquares_2Dx && numberOfElements <= int(nsd_real)) ||
        (SmoothGradPhi == INPAR::COMBUST::smooth_grad_phi_leastsquares_2Dy && numberOfElements <= int(nsd_real)) ||
        (SmoothGradPhi == INPAR::COMBUST::smooth_grad_phi_leastsquares_2Dz && numberOfElements <= int(nsd_real))    )
    {
      CurrTypeOfNodeReconst = INPAR::COMBUST::smooth_grad_phi_meanvalue;

      cout << "warning:\n"
          << "\t\tnode\t" << (it_node->second)->Id() << "\tis a boundary node with "<< numberOfElements << " elements and is reconstructed with average (mean value) method" << endl;
    }
    else // type of reconstruction as chosen in input file
      CurrTypeOfNodeReconst = SmoothGradPhi;

    // ====================================================================================================
    // REMARK:
    // special case for intersected elements which have boundary nodes
    // for boundary nodes there are not enough elements to reconstruct node gradient with least squares
    // => build average of gradients of adjacent elements instead of least squares reconstruction
    // at least nsd+1 adjacent elements necessary for least squares reconstruction, we take at least numnode!!!
    if(CurrTypeOfNodeReconst == INPAR::COMBUST::smooth_grad_phi_meanvalue)  //(int)numnode
    {
      //==============================================================================
      //      average (mean value) reconstruction for boundary nodes and as alternative
      //==============================================================================

      // get adjacent elements to current node actnode
      const DRT::Element* const* elements = ptToNode->Elements();

      // loop over adjacent elements
      for(int ele_current=0; ele_current<numberOfElements; ele_current++)
      {
        // get current element
        const DRT::Element* ele_adj = elements[ele_current];

        const int* ptToNodeIds_adj = ele_adj->NodeIds();

        // get phi-values of current adjacent element ele_adj
        // create vector "ephinp" holding scalar phi values for this element
        Epetra_SerialDenseVector ephinp(numnode); //local vector phi-values of adjacent element

        // which node in param space of element ele_adj has actnode
        int ID_param_space = -1;

        // get vector of node GIDs of this adjacent element -> needed for ExtractMyValues
        vector<int> nodeID_adj(numnode);
        for (size_t inode=0; inode < numnode; inode++){
          nodeID_adj[inode] = ptToNodeIds_adj[inode];
          // get local number of node actnode in ele_adj
          if(ptToNode->Id() == ptToNodeIds_adj[inode]) ID_param_space = inode;
        }

        // extract the phi-values of adjacent element with local ids from global vector *phinp
        // get pointer to vector holding G-function values at the fluid nodes
        DRT::UTILS::ExtractMyValues(*phinp_, ephinp, nodeID_adj);
        LINALG::Matrix<numnode,1> ephi_adj(ephinp);

        //===============================================================================
        // calculate gradient of phi at this node

        // get Xi-coordinates of node for current element

        // get derivatives of shapefunctions evaluated at node in XYZ-coordinates
        static LINALG::Matrix<nsd,numnode> deriv3Dele_xyz;
        // get derivatives of shapefunctions evaluates at node in Xi-coordinates
        static LINALG::Matrix<nsd,numnode> deriv3Dele;

        // get Xi-coordinates of current node in current adjacent element
        static LINALG::Matrix<nsd,1> node_Xicoordinates;
        node_Xicoordinates.Clear();
        node_Xicoordinates = DRT::UTILS::getNodeCoordinates(ID_param_space, DISTYPE);


        // get derivatives of shapefunctions at node
        DRT::UTILS::shape_function_3D_deriv1(deriv3Dele,node_Xicoordinates(0),node_Xicoordinates(1),node_Xicoordinates(2),DISTYPE);

        // reconstruct XYZ!-gradient

        // get node coordinates of this element
        static LINALG::Matrix<nsd,numnode> xyze_adj;
        GEO::fillInitialPositionArray<DISTYPE>(ele_adj, xyze_adj);

        // get Jacobi-Matrix for transformation
        static LINALG::Matrix<nsd,nsd> xjm_ele_XiToXYZ;
        xjm_ele_XiToXYZ.MultiplyNT(deriv3Dele,xyze_adj);

        // inverse of jacobian
        static LINALG::Matrix<nsd,nsd> xji_ele_XiToXYZ;
        xji_ele_XiToXYZ.Invert(xjm_ele_XiToXYZ);

        // set XYZ-derivates of shapefunctions
        deriv3Dele_xyz.Multiply(xji_ele_XiToXYZ,deriv3Dele);

        //===============================================================================
        // calculate gradient of phi at node for current element
        static LINALG::Matrix<3,1> nodal_grad_tmp;
        nodal_grad_tmp.Clear();

        // get xyz-gradient
        nodal_grad_tmp.Multiply(deriv3Dele_xyz, ephi_adj);

        //===============================================================================
        // add to vector with smoothed vector
        // we dont need PHI_SMOOTHED(0,0)

        // full 3D reconstruction, also for 2D examples
        PHI_SMOOTHED_3D(0,0) += nodal_grad_tmp(0,0);
        PHI_SMOOTHED_3D(1,0) += nodal_grad_tmp(1,0);
        PHI_SMOOTHED_3D(2,0) += nodal_grad_tmp(2,0);

      }// end loop over all adjacent elements
#if(0)
      // special case for straight_bodyforce
      // this is a node of an intersected element
      // the element must be a boundary element!
      // assume a continued interface across the domain boundary
      if(numberOfElements == 2)
      {
        // set two vectors in y-direction
        PHI_SMOOTHED_3D(0,0) = 1.0;
        PHI_SMOOTHED_3D(1,0) = 1.0;
        PHI_SMOOTHED_3D(2,0) = 0.0;

        cout << "\n\t !!! warning !!! (Rayleigh-Taylor modification) we modify the gradient of phi periodically at the domain boundary";
      }
#endif

#if(0)
      // special case for Rayleigh-Taylor
      // this is a node of an intersected element
      // the element must be a boundary element!
      // assume a continued interface across the domain boundary
      if(numberOfElements == 2)
      {
        // set two vectors in y-direction
        PHI_SMOOTHED_3D(0,0) = 0.0;
        PHI_SMOOTHED_3D(1,0) = 2.0;
        PHI_SMOOTHED_3D(2,0) = 0.0;

        cout << "\n\t !!! warning !!! (Rayleigh-Taylor modification) we modify the gradient of phi periodically at the domain boundary";
      }
#endif

#if(0)
      // special case for flame vortex interaction
      // this is a node of an intersected element
      // the element must be a boundary element!
      // assume a continued interface across the domain boundary
      if(numberOfElements == 2 || numberOfElements == 1)
      {
        // set two vectors in y-direction
        PHI_SMOOTHED_3D(0,0) = 0.0;
        PHI_SMOOTHED_3D(1,0) = 1.0 * numberOfElements;
        PHI_SMOOTHED_3D(2,0) = 0.0;

        cout << "\n\t !!! warning !!! (flame_vortex_interaction modification) we modify the gradient of phi periodically at the domain boundary";
      }
#endif


      // weight sum of nodal_grad_tmp 1/number_of_vectors to get an average value
      PHI_SMOOTHED_3D.Scale(1.0/numberOfElements);
    } // end of average (mean value) reconstruction
    else // least squares reconstruction
    {
      //==============================================================================
      //      standard reconstruction with least squares method, see Merchandise'07
      //==============================================================================

      // we need Epetra-Matrix (numberOfElements is not a valid template parameter)

      Epetra_SerialDenseMatrix RHS_LS(numberOfElements,1);
      Epetra_SerialDenseMatrix MAT_LS(numberOfElements, nsd_real);

      // get adjacent elements to current node actnode
      const DRT::Element* const* elements = ptToNode->Elements();

      int xyz_2D_dim = -1;
      // set one dimension (xyz-dim) to zero for 2D LS-reconstruction
      if(SmoothGradPhi == INPAR::COMBUST::smooth_grad_phi_leastsquares_2Dx) xyz_2D_dim = 0;
      if(SmoothGradPhi == INPAR::COMBUST::smooth_grad_phi_leastsquares_2Dy) xyz_2D_dim = 1;
      if(SmoothGradPhi == INPAR::COMBUST::smooth_grad_phi_leastsquares_2Dz) xyz_2D_dim = 2;


      // loop over adjacent elements
      for(int ele_current=0; ele_current<numberOfElements; ele_current++)
      {
        // get current element
        const DRT::Element* ele_adj = elements[ele_current];

        const int* ptToNodeIds_adj = ele_adj->NodeIds();

        // get phi-values of current adjacent element ele_adj
        // create vector "ephinp" holding scalar phi values for this element
        Epetra_SerialDenseVector ephinp(numnode); //local vector phi-values of adjacent element

        // which node in param space of element ele_adj has actnode
        int ID_param_space = -1;

        // get vector of node GIDs of this adjacent element -> needed for ExtractMyValues
        vector<int> nodeID_adj(numnode);
        for (size_t inode=0; inode < numnode; inode++){
          nodeID_adj[inode] = ptToNodeIds_adj[inode];
          // get local number of node actnode in ele_adj
          if(ptToNode->Id() == ptToNodeIds_adj[inode]) ID_param_space = inode;
        }
        if (ID_param_space == -1) dserror ("node in current adjacent not found!!!");

        // extract the phi-values of adjacent element with local ids from global vector *phinp
        // get pointer to vector holding G-function values at the fluid nodes
        DRT::UTILS::ExtractMyValues(*phinp_, ephinp, nodeID_adj);
        LINALG::Matrix<numnode,1> ephi_adj(ephinp);


        // calculate center of gravity
        static LINALG::Matrix<numnode,1> funct;
        funct.Clear();
        static LINALG::Matrix<nsd,1> centerOfGravXi;
        centerOfGravXi.Clear();

        // xi-coord = 0.0.0 for hex8
        if (DISTYPE != DRT::Element::hex8) dserror("center of gravity implemented only for hex8 elements");
        centerOfGravXi(0) = 0.0;
        centerOfGravXi(1) = 0.0;
        centerOfGravXi(2) = 0.0;



        DRT::UTILS::shape_function_3D(funct,centerOfGravXi(0),centerOfGravXi(1),centerOfGravXi(2),DISTYPE);

        // get node coordinates of this element
        static LINALG::Matrix<nsd,numnode> xyze_adj;
        GEO::fillInitialPositionArray<DISTYPE>(ele_adj, xyze_adj);

        // interpolate center of gravity
        static LINALG::Matrix<nsd,1> centerOfGravXYZ;
        centerOfGravXYZ.Clear();

        centerOfGravXYZ.Multiply(xyze_adj,funct);

        // reconstruct xyz-gradient
        // calculate vector from current node to point of gravity (direction for taylor)
        static LINALG::Matrix<nsd,1> direction;
        direction.Clear();

        direction(0) = centerOfGravXYZ(0) - xyze_adj(0, ID_param_space);
        direction(1) = centerOfGravXYZ(1) - xyze_adj(1, ID_param_space);
        direction(2) = centerOfGravXYZ(2) - xyze_adj(2, ID_param_space);



        // cancel out the 2D direction
        if(xyz_2D_dim != -1) direction(xyz_2D_dim) = 0.0;


        // assemble direction into matrix A

        // calculate ephi at point of gravity via interpolation
        static LINALG::Matrix<1,1> phi_adj;
        phi_adj.Clear();
        phi_adj.MultiplyTN(ephi_adj,funct);

        // set RHS and MAT for least squares method

        if(SmoothGradPhi == INPAR::COMBUST::smooth_grad_phi_leastsquares_3D){
          RHS_LS(ele_current,0) = phi_adj(0,0)-ephi_adj(ID_param_space,0);
          MAT_LS(ele_current,0) = direction(0);
          MAT_LS(ele_current,1) = direction(1);
          MAT_LS(ele_current,2) = direction(2);
        }
        else{ // assemble MAT_LS for 2D cases
          RHS_LS(ele_current,0) = phi_adj(0,0)-ephi_adj(ID_param_space,0);

          // 2D LEAST_SQUARES cases
          size_t MAT_isd = 0;
          for(size_t isd = 0; isd<nsd; isd++)
          {
            if(int(isd)!=xyz_2D_dim){// we have to not!!! assemble the xyz_2D_dim
              MAT_LS(ele_current,MAT_isd) = direction(isd);
              MAT_isd++;
            }
          } // end for
        } // end else

      } // =================================end loop over all adjacent elements===========================

      // the system MAT_LS * phi_smoothed = RHS_LS is only solvable in a least squares manner
      // MAT_LS is not square
      // -> solve MAT_LS^T * MAT_LS * phi_smoothed = MAT_LS^T * RHS_LS
      Epetra_SerialDenseMatrix MAT_tmp(nsd_real,nsd_real);
      Epetra_SerialDenseMatrix RHS_tmp(nsd_real,1);

      MAT_tmp.Multiply('T','N', 1.0, MAT_LS,MAT_LS, 0.0);
      RHS_tmp.Multiply('T','N', 1.0, MAT_LS,RHS_LS, 0.0);

      // set LINALG-Matrix with fixed size
      // initialize element matrix for A^T*A and RHS-vector to solve Normalengleichung

      LINALG::Matrix<nsd_real,nsd_real> MAT(MAT_tmp);
      LINALG::Matrix<nsd_real,1> RHS(RHS_tmp);

      //=================solve A^T*A* phi_smoothed = A^T*RHS (LEAST SQUARES)================

      // solve the system for current node
      LINALG::FixedSizeSerialDenseSolver<nsd_real,nsd_real,1> solver; //1 is dimension of RHS
      solver.SetMatrix(MAT);
      solver.SetVectors(PHI_SMOOTHED,RHS);
      solver.Solve();

      //      cout << " node:" << (it_node->second)->Id() << endl;
      //      cout << "PHI_SMOOTHED" << PHI_SMOOTHED << endl;
      // set full 3D vector especially important for 2D leastsquares reconstruction
      if(CurrTypeOfNodeReconst == INPAR::COMBUST::smooth_grad_phi_leastsquares_3D)
      {
        PHI_SMOOTHED_3D(0,0) = PHI_SMOOTHED(0,0);
        PHI_SMOOTHED_3D(1,0) = PHI_SMOOTHED(1,0);
        PHI_SMOOTHED_3D(2,0) = PHI_SMOOTHED(2,0);
      }
      else // 2D-LeastSquares-cases
      {
        size_t PHI_2D_isd = 0;
        for(size_t isd = 0; isd<nsd; isd++)
        {
          if(int(isd)!=xyz_2D_dim){ // we have to not!!! assemble the xyz_2D_dim
            PHI_SMOOTHED_3D(isd,0) = PHI_SMOOTHED(PHI_2D_isd,0);
            PHI_2D_isd++;
          }
        }
      }

    } // end of standard (Least_squares-) reconstruction case


    //=====================================================================================================
    // set the global vector gradphirow holding the new reconstructed values of gradient of phi in row!!! map
    //=====================================================================================================

    // get the global id for current node
    int GID = (it_node->second)->Id();

    // get local processor id according to global node id
    const int lid = (*gradphirow).Map().LID(GID);

    if (lid<0) dserror("Proc %d: Cannot find gid=%d in Epetra_Vector",(*gradphi_).Comm().MyPID(),GID);

    const int numcol = (*gradphirow).NumVectors();
    if( numcol != (int)nsd) dserror("number of columns in Epetra_MultiVector is not identically to nsd");

    // loop over dimensions (= number of columns in multivector)
    for(int col=0; col< numcol; col++)
    {
      // get columns vector of multivector
      double* globalcolumn = (*gradphirow)[col];


      // set smoothed gradient entry of phi into column of global multivector
      globalcolumn[lid] = PHI_SMOOTHED_3D(col,0);
    }

  }// ====================================== end loop over nodes ==========================================

  // export NodeRowMap to NodeColMap gradphi_
  LINALG::Export(*gradphirow,*gradphi_);

  std::cout << "done" << std::endl;



  //
  //  //===================================== new 2-Ring version
  //  for(Reconstruct_iterator it_node = nodesToReconstruct.begin(); it_node!= nodesToReconstruct.end(); it_node++)
  //  {
  //    // define vector for smoothed values (Phi, grad_Phi) of current node
  //    static LINALG::Matrix<9,1> PHI_SMOOTHED; // for real reconstruction 2D or 3D
  //    PHI_SMOOTHED.Clear();
  //    static LINALG::Matrix<9,1> PHI_SMOOTHED_3D; // set whole 3D vector also for 2D examples
  //    PHI_SMOOTHED_3D.Clear();
  //
  //
  //    // get local processor id of current node and pointer to current node
  //    //int lid_node = it_node->first;
  //    const DRT::Node* ptToNode = it_node->second;
  ////    cout << "node:" << (it_node->second)->Id();
  //
  //    int numberOfElements = ptToNode->NumElement();
  //
  //    if(true) // least squares reconstruction
  //    {
  //        // get all adjacent elements in a 2-Ring
  //        // get adjacent elements to current node actnode (1-Ring)
  //        const DRT::Element* const* elementsOneRing = ptToNode->Elements();
  //
  //
  //        // map of pointers to nodes which must be reconstructed by this processor
  //        // key is the local id at processor
  //        std::map<int, const DRT::Element*> elesinTwoRing;
  //        elesinTwoRing.clear();
  //
  //        // loop over adjacent elements (1-ring elements)
  //        for(int ele_current=0; ele_current<numberOfElements; ele_current++)
  //        {
  //          // get current element
  //          const DRT::Element* ele_adj = elementsOneRing[ele_current];
  //
  //          // these are global Ids
  //          const int* ptToNodeIds_adj = ele_adj->NodeIds();
  //
  //
  //          for(int currentNode = 0; currentNode < (int)numnode; currentNode++)
  //          {
  //        	  // get local ID
  //              int curr_GID = ptToNodeIds_adj[currentNode];
  //
  //              // get local processor id according to global node id
  //              const int curr_lid = (*phinp_).Map().LID(curr_GID);
  //
  //              if (curr_lid<0) dserror("Proc %d: Cannot find gid=%d in Epetra_Vector",(*gradphi_).Comm().MyPID(),curr_GID);
  //
  //
  //        	  // get current Node, the local ID is needed
  //        	  const DRT::Node *actnode = fluiddis_->lColNode(curr_lid);
  //
  //        	  // get pointer to all adjacent elements of this node actnode
  //              const DRT::Element* const* elementsTwoRing = actnode->Elements();
  //
  //              int numberOfElementsTwoRing = actnode->NumElement();
  //
  //        	  // get all elements of the 2-ring and insert these in a vector
  //              // loop over these elements, insert to map
  //              for(int ele_current_two_ring=0; ele_current_two_ring < numberOfElementsTwoRing; ele_current_two_ring++)
  //              {
  //                  // get current element
  //                  const DRT::Element* ele_Two_Ring = elementsTwoRing[ele_current_two_ring];
  //
  //                  // insert in map of all two-ring elements respect to one node
  //                  int lid = elementsTwoRing[ele_current_two_ring]->LID();
  //                  elesinTwoRing[lid] = ele_Two_Ring;
  //              } // end loop over elements in 2-ring
  //          } // end loop over nodes of 1-ring
  //        } // end loop over elements of 1-ring
  //
  //        // get node xyze-coordinates
  //        LINALG::Matrix<nsd,1> xyz_node(ptToNode->X());
  ////        cout << "xyz_node" << xyz_node << endl;
  //
  //
  //    	const int numberOfTwoRingElements = elesinTwoRing.size();
  ////    	cout << "numberOfTwoRingElements" << numberOfTwoRingElements << endl;
  //      Epetra_SerialDenseMatrix RHS_LS(numberOfTwoRingElements,1);
  //      Epetra_SerialDenseMatrix MAT_LS(numberOfTwoRingElements, 9);
  //
  //      // map iterator
  //      typedef std::map<int, const DRT::Element*> Map_EleOfTwoRing;
  //      typedef Map_EleOfTwoRing::iterator TwoRing_iterator;
  //      int element_counter = 0;
  //
  //
  ////      cout << "=====================iterate over elements in two-ring==============" << endl;
  //      for(TwoRing_iterator it_ele = elesinTwoRing.begin(); it_ele!= elesinTwoRing.end(); it_ele++)
  //      {
  //
  //          // get adjacent elements to current node actnode
  //          const DRT::Element* ele_Two_Ring= it_ele->second;
  //
  //          // get pointer to NodeIds of current element (global Ids)
  //          const int* ptToNodeIds_adj = ele_Two_Ring->NodeIds();
  //
  //
  //          // get phi-values of current adjacent element ele_adj
  //          // create vector "ephinp" holding scalar phi values for this element
  //          Epetra_SerialDenseVector ephinp(numnode); //local vector phi-values of adjacent element
  //
  ////          cout << "these are all global node Ids of element" << ele_Two_Ring->Id() << endl;
  //          // get vector of node GIDs of this adjacent element -> needed for ExtractMyValues
  //          vector<int> nodeID_adj(numnode);
  //          for (size_t inode=0; inode < numnode; inode++){
  //            nodeID_adj[inode] = ptToNodeIds_adj[inode];
  ////            cout << nodeID_adj[inode] << endl;
  //          }
  //
  //          // extract the phi-values of adjacent element with local ids from global vector *phinp
  //          // get pointer to vector holding G-function values at the fluid nodes
  //          DRT::UTILS::ExtractMyValues(*phinp_, ephinp, nodeID_adj);
  //          LINALG::Matrix<numnode,1> ephi_adj(ephinp);
  //
  //
  //          // calculate center of gravity
  //          static LINALG::Matrix<numnode,1> funct;
  //          funct.Clear();
  //          static LINALG::Matrix<nsd,1> centerOfGravXi;
  //          centerOfGravXi.Clear();
  //
  //          // xi-coord = 0.0.0 for hex8
  //          if (DISTYPE != DRT::Element::hex8) dserror("center of gravity implemented only for hex8 elements");
  //          centerOfGravXi(0) = 0.0;
  //          centerOfGravXi(1) = 0.0;
  //          centerOfGravXi(2) = 0.0;
  //
  //
  //
  //          DRT::UTILS::shape_function_3D(funct,centerOfGravXi(0),centerOfGravXi(1),centerOfGravXi(2),DISTYPE);
  //
  //          // get node coordinates of this element
  //          static LINALG::Matrix<nsd,numnode> xyze_adj;
  //          GEO::fillInitialPositionArray<DISTYPE>(ele_Two_Ring, xyze_adj);
  //
  //          // interpolate center of gravity
  //          static LINALG::Matrix<nsd,1> centerOfGravXYZ;
  //          centerOfGravXYZ.Clear();
  //
  //          centerOfGravXYZ.Multiply(xyze_adj,funct);
  //
  //          // reconstruct xyz-gradient
  //          // calculate vector from current node to point of gravity (direction for taylor)
  //          static LINALG::Matrix<nsd,1> direction;
  //          direction.Clear();
  //
  //          direction(0) = centerOfGravXYZ(0) - xyz_node(0);
  //          direction(1) = centerOfGravXYZ(1) - xyz_node(1);
  //          direction(2) = centerOfGravXYZ(2) - xyz_node(2);
  //
  ////          cout << "direction" << direction << endl;
  //
  //          // calculate ephi at point of gravity via interpolation
  //          static LINALG::Matrix<1,1> phi_adj;
  //          phi_adj.Clear();
  //          phi_adj.MultiplyTN(ephi_adj,funct);
  //
  ////          cout << "phi_adj" << phi_adj << endl;
  //
  ////          // cancel out the 2D direction
  ////          if(xyz_2D_dim != -1) direction(xyz_2D_dim) = 0.0;
  //          // get the global id for current node
  //          int GID = ptToNode->Id();
  //
  //          // get local processor id according to global node id
  //          const int lid = (*phinp_).Map().LID(GID);
  //
  //          if (lid<0) dserror("Proc %d: Cannot find gid=%d in Epetra_Vector",(*gradphi_).Comm().MyPID(),GID);
  //
  //          double phi_node_value = (*phinp_)[lid];
  ////          cout << "phi_node_value " << phi_node_value << endl;
  //
  //          // set RHS and MAT for least squares method
  //
  //            RHS_LS(element_counter,0) = phi_adj(0,0) -phi_node_value;
  //
  //            MAT_LS(element_counter,0) = direction(0);
  //            MAT_LS(element_counter,1) = direction(1);
  //            MAT_LS(element_counter,2) = direction(2);
  //
  //            MAT_LS(element_counter,3) = 0.5*direction(0)*direction(0);
  //            MAT_LS(element_counter,4) = 0.5*direction(1)*direction(1);
  //            MAT_LS(element_counter,5) = 0.5*direction(2)*direction(2);
  //
  //            MAT_LS(element_counter,6) = direction(0)*direction(1);
  //            MAT_LS(element_counter,7) = direction(0)*direction(2);
  //            MAT_LS(element_counter,8) = direction(1)*direction(2);
  ////            MAT_LS(element_counter,0) = 1.0;
  ////            MAT_LS(element_counter,1) = direction(0);
  ////            MAT_LS(element_counter,2) = direction(1);
  ////            MAT_LS(element_counter,3) = direction(2);
  ////
  ////            MAT_LS(element_counter,4) = 0.5*direction(0)*direction(0);
  ////            MAT_LS(element_counter,5) = 0.5*direction(1)*direction(1);
  ////            MAT_LS(element_counter,6) = 0.5*direction(2)*direction(2);
  ////
  ////            MAT_LS(element_counter,7) = direction(0)*direction(1);
  ////            MAT_LS(element_counter,8) = direction(0)*direction(2);
  ////            MAT_LS(element_counter,9) = direction(1)*direction(2);
  //
  //
  //            element_counter++;
  //      } // end loop over two-ring-elements
  //
  ////      cout << "MATLS" << endl;
  ////      for(int row = 0; row < element_counter; row++)
  ////      {
  ////    	  for(int col = 0; col < 10; col++)
  ////    	  {
  ////    		  cout << MAT_LS(row, col) << "\t";
  ////    	  }
  ////    	  cout << endl;
  ////      }
  ////      cout << endl;
  //
  //      // the system MAT_LS * phi_smoothed = RHS_LS is only solvable in a least squares manner
  //      // MAT_LS is not square
  //      // -> solve MAT_LS^T * MAT_LS * phi_smoothed = MAT_LS^T * RHS_LS
  //      Epetra_SerialDenseMatrix MAT_tmp(9,9);
  //      Epetra_SerialDenseMatrix RHS_tmp(9,1);
  //
  //      MAT_tmp.Multiply('T','N', 1.0, MAT_LS,MAT_LS, 0.0);
  //      RHS_tmp.Multiply('T','N', 1.0, MAT_LS,RHS_LS, 0.0);
  //
  //      // set LINALG-Matrix with fixed size
  //      // initialize element matrix for A^T*A and RHS-vector to solve Normalengleichung
  //
  //      LINALG::Matrix<9,9> MAT(MAT_tmp);
  //      LINALG::Matrix<9,1> RHS(RHS_tmp);
  //
  ////      cout << "MAT" << MAT << endl;
  ////      cout << "RHS" << RHS << endl;
  //      //=================solve A^T*A* phi_smoothed = A^T*RHS (LEAST SQUARES)================
  //
  //      // solve the system for current node
  //      LINALG::FixedSizeSerialDenseSolver<9,9,1> solver; //1 is dimension of RHS
  //      solver.SetMatrix(MAT);
  //      solver.SetVectors(PHI_SMOOTHED,RHS);
  //      solver.Solve();
  //
  ////      cout << " node:" << (it_node->second)->Id() << endl;
  ////      cout << "============PHI_SMOOTHED==========" << PHI_SMOOTHED << endl;
  //      // set full 3D vector especially important for 2D leastsquares reconstruction
  ////        PHI_SMOOTHED_3D(0,0) = PHI_SMOOTHED(1,0);
  ////        PHI_SMOOTHED_3D(1,0) = PHI_SMOOTHED(2,0);
  ////        PHI_SMOOTHED_3D(2,0) = PHI_SMOOTHED(3,0);
  ////        PHI_SMOOTHED_3D(3,0) = PHI_SMOOTHED(4,0);
  ////        PHI_SMOOTHED_3D(4,0) = PHI_SMOOTHED(5,0);
  ////        PHI_SMOOTHED_3D(5,0) = PHI_SMOOTHED(6,0);
  ////        PHI_SMOOTHED_3D(6,0) = PHI_SMOOTHED(7,0);
  ////        PHI_SMOOTHED_3D(7,0) = PHI_SMOOTHED(8,0);
  ////        PHI_SMOOTHED_3D(8,0) = PHI_SMOOTHED(9,0);
  //
  //        PHI_SMOOTHED_3D(0,0) = PHI_SMOOTHED(0,0);
  //        PHI_SMOOTHED_3D(1,0) = PHI_SMOOTHED(1,0);
  //        PHI_SMOOTHED_3D(2,0) = PHI_SMOOTHED(2,0);
  //        PHI_SMOOTHED_3D(3,0) = PHI_SMOOTHED(3,0);
  //        PHI_SMOOTHED_3D(4,0) = PHI_SMOOTHED(4,0);
  //        PHI_SMOOTHED_3D(5,0) = PHI_SMOOTHED(5,0);
  //        PHI_SMOOTHED_3D(6,0) = PHI_SMOOTHED(6,0);
  //        PHI_SMOOTHED_3D(7,0) = PHI_SMOOTHED(7,0);
  //        PHI_SMOOTHED_3D(8,0) = PHI_SMOOTHED(8,0);
  //
  //    } // end of standard (Least_squares-) reconstruction case
  //
  //
  //    //=====================================================================================================
  //    // set the global vector gradphirow holding the new reconstructed values of gradient of phi in row!!! map
  //    //=====================================================================================================
  //
  //    // get the global id for current node
  //    int GID = (it_node->second)->Id();
  //
  //    // get local processor id according to global node id
  //    const int lid = (*gradphirow).Map().LID(GID);
  //
  //    if (lid<0) dserror("Proc %d: Cannot find gid=%d in Epetra_Vector",(*gradphi_).Comm().MyPID(),GID);
  //
  //    const int numcol = (*gradphirow).NumVectors();
  ////    if( numcol != (int)nsd) dserror("number of columns in Epetra_MultiVector is not identically to nsd");
  //
  //    // loop over dimensions (= number of columns in multivector)
  //    for(int col=0; col< numcol; col++)
  //    {
  //      // get columns vector of multivector
  //      double* globalcolumn = (*gradphirow)[col];
  //
  //
  //      // set smoothed gradient entry of phi into column of global multivector
  ////      globalcolumn[lid] = PHI_SMOOTHED_3D(col,0);
  //      globalcolumn[lid] = PHI_SMOOTHED_3D(col,0);
  //    }
  //
  //  }// ====================================== end loop over nodes ==========================================
  //
  //  // export NodeRowMap to NodeColMap gradphi_
  //  LINALG::Export(*gradphirow,*gradphi_);
  //
  //  std::cout << "done" << std::endl;
  //
  //


  return;
}


/*------------------------------------------------------------------------------------------------*
 | call function for different types of computation of smoothed gradient field of G-function      |
 | (level set field) for use in fluid time                                                        |
 | integration scheme                                                                 henke 06/10 |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::FlameFront::CallSmoothGradPhi(const Teuchos::ParameterList& combustdyn)
{
  TEUCHOS_FUNC_TIME_MONITOR("SMOOTHING OF GRAD_PHI");
  if (gfuncdis_->Comm().MyPID()==0)
    std::cout << "\n---  smoothing gradient of phi (G-function) around the interface for surface tension applications ... " << std::flush;

  // get type of reconstruction
  const INPAR::COMBUST::SmoothGradPhi SmoothGradPhi = DRT::INPUT::IntegralValue<INPAR::COMBUST::SmoothGradPhi>
  (combustdyn.sublist("COMBUSTION FLUID"),"SMOOTHGRADPHI");


  // standard dimension is 3, real dimension for least squares reconstruction type is necessary

  size_t nsd_real = 3;
  if( SmoothGradPhi == INPAR::COMBUST::smooth_grad_phi_leastsquares_2Dx ||
      SmoothGradPhi == INPAR::COMBUST::smooth_grad_phi_leastsquares_2Dy ||
      SmoothGradPhi == INPAR::COMBUST::smooth_grad_phi_leastsquares_2Dz ) nsd_real = 2;

  // cout reconstruction type
  if(SmoothGradPhi == INPAR::COMBUST::smooth_grad_phi_leastsquares_2Dx)
    std::cout << "\n---  \t reconstruction with:\t LeastSquares_2Dx... " << std::flush;
  if(SmoothGradPhi == INPAR::COMBUST::smooth_grad_phi_leastsquares_2Dy)
    std::cout << "\n---  \t reconstruction with:\t LeastSquares_2Dy... " << std::flush;
  if(SmoothGradPhi == INPAR::COMBUST::smooth_grad_phi_leastsquares_2Dz)
    std::cout << "\n---  \t reconstruction with:\t LeastSquares_2Dz... " << std::flush;
  if(SmoothGradPhi == INPAR::COMBUST::smooth_grad_phi_leastsquares_3D)
    std::cout << "\n---  \t reconstruction with:\t LeastSquares_3D... " << std::flush;
  if(SmoothGradPhi == INPAR::COMBUST::smooth_grad_phi_meanvalue)
    std::cout << "\n---  \t reconstruction with:\t MeanValue... " << std::flush;

  // switch: real dimension of reconstruction, is real dimension 2D or 3D?
  switch (nsd_real)
  {
  case 3:
    ComputeSmoothGradPhi<3>(combustdyn);
    break;
  case 2:
    ComputeSmoothGradPhi<2>(combustdyn);
    break;
  default:
    dserror("wrong nsd_real");
  };

  if (gfuncdis_->Comm().MyPID()==0)
    std::cout << " done" << std::endl;

  return;
}

/*------------------------------------------------------------------------------------------------*
 | refine the region around the flame front                                           henke 10/08 |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::FlameFront::RefineFlameFront(const Teuchos::RCP<COMBUST::RefinementCell> cell,
    const Teuchos::RCP<const Epetra_Vector> phi)
{
  /* Hier wird der rekursive Verfeinerungsprozess gestartet.

  rufe FindFlameFront() für eine Verfeinerungszelle
  falls Rückgabe = true (Zelle wird geschnitten)
    teile die Verfeinerungszelle
    -> RefineCell() (input: Zelle; output: geteilte Zellen)
    für alle neuen Verfeinerungszellen
      falls maximale Anzahl von Verfeinerungen noch nicht erreicht
        rekursiver Aufruf: RefineFlamefront() (RefineCell())
   */

  // get G-Function values of refinement cell and determine, if it is intersected
  FindFlameFront(cell,phi);
  // is maximal refinement level already reached?
  if (cell->RefinementLevel() < maxRefinementLevel_) //->No
  {
    if (cell->Bisected()) // cell is bisected ....
    {
      cell->RefineCell(); // ... and will be refined
      //        std::cout << "RefineCell done" << std::endl;
      // loop over all new refinement cells and call of RefineFlameFront
      for (int icell = 0; icell < cell->NumOfChildren(); icell++)
        RefineFlameFront(cell->GetRefinementCell(icell),phi);
    }
  }

  return;
}


/*------------------------------------------------------------------------------------------------*
 | find the flame front within a refinement cell according to G-function field        henke 10/08 |
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
      vector<int> lm(ele->NumNode());
      // get vector of node GIDs for this element
      const int* nodeids = ele->NodeIds();
      for (unsigned inode=0; inode < lm.size(); inode++)
        lm[inode] = nodeids[inode];
      /*
  weitere Implementierungsmöglichkeiten
  2.Möglichkeit, führt Schleife für NodeIds() von Hand aus und spart so Umschreiben vin nodeids auf lm
      // get pointer to nodes of this element
      const DRT::Node* const* nodes = ele->Nodes();
      // loop all nodes of this element
      for (int inode = 0; inode < ele->NumNode(); inode++)
      {
        // get GID of this node (fluid node!)
        const int nodeid = nodes[inode]->Id();
        // build vector of all node GIDs of this element
        lm[inode] = nodeid;
      }
  3.Möglichkeit, die funktioniert, wenn phi auf der gfuncdis DofColMap vorliegen würde.
        // get vector of all dof GIDs of this element
        const vector<int> lm = gfuncdis_->Dof(ele);
        // create local vector "myphinp"
        vector<double> myphinp(lm.size()); // vector<double> degeneriert hier zum scalaren double!
       */
      //--------------------------------------------------------------------------------------------
      // extract G-function values for all nodes belonging to this fluid element ( == cell!)
      //--------------------------------------------------------------------------------------------
      // create vector "mygfuncvalues" holding G-function values for this element
      vector<double> mygfuncvalues(ele->NumNode(),1000.0);
      // get entries in "gfuncvalues" corresponding to node GIDs "lm" and store them in "mygfuncvalues"
      DRT::UTILS::ExtractMyValues(*phi,mygfuncvalues,lm);
      //TEST
      //if (ele->Id()==0)
      //{
      //  std::cout<< "Gfunc " << ele->Id() << std::endl;
      //  for(std::size_t ig=0; ig<mygfuncvalues.size(); ig++)
      //    std::cout << mygfuncvalues[ig] << std::endl;
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

      //TEST Einheitswürfel mit verschiedenen Schnitten (Hex8)
      //        vector<double> phis (8);
      //        phis[0]=-1;
      //        phis[1]=-1;
      //        phis[2]=1;
      //        phis[3]=1;
      //        phis[4]=1;
      //        phis[5]=-1;
      //        phis[6]=1;
      //        phis[7]=2;
      //
      //        phis[0]=1;
      //        phis[1]=1;
      //        phis[2]=0;
      //        phis[3]=0;
      //        phis[4]=0;
      //        phis[5]=0;
      //        phis[6]=-1;
      //        phis[7]=-1;
      //
      //        phis[0]=-1.01;
      //        phis[1]=0.99;
      //        phis[2]=0.99;
      //        phis[3]=-1.01;
      //        phis[4]=-1.01;
      //        phis[5]=0.99;
      //        phis[6]=0.99;
      //        phis[7]=-1.01;
      //
      //        phis[0]=-1.5;
      //        phis[1]=0.5;
      //        phis[2]=-1.5;
      //        phis[3]=0.5;
      //        phis[4]=-1.5;
      //        phis[5]=0.5;
      //        phis[6]=-1.5;
      //        phis[7]=0.5;
      //
      //        phis[0]=0;
      //        phis[1]=0;
      //        phis[2]=0;
      //        phis[3]=0;
      //        phis[4]=0;
      //        phis[5]=0;
      //        phis[6]=0;
      //        phis[7]=1;
      //
      //       phis[0]=1.5;
      //       phis[1]=-0.5;
      //       phis[2]=1.5;
      //       phis[3]=2;
      //       phis[4]=1.5;
      //       phis[5]=-0.5;
      //       phis[6]=1.5;
      //       phis[7]=2;
      //
      //       phis[0]=1;
      //       phis[1]=-1;
      //       phis[2]=-1;
      //       phis[3]=1;
      //       phis[4]=1;
      //       phis[5]=-1;
      //       phis[6]=-1;
      //       phis[7]=1;
      //
      //       phis[0]=1.0;
      //       phis[1]=1.5;
      //       phis[2]=-3.0;
      //       phis[3]=-3.0;
      //       phis[4]=1.5;
      //       phis[5]=1;
      //       phis[6]=3.0;
      //       phis[7]=-3.0;
      //
      //             phis[0]=-8.642922e-03;
      //             phis[1]=-8.622150e-03;
      //             phis[2]=1.790659e-02;
      //             phis[3]=1.647029e-02;
      //             phis[4]=9.449617e-03;
      //             phis[5]=8.525568e-03;
      //             phis[6]=-7.705337e-04;
      //             phis[7]=2.363820e-04;
      //              cell->SetGfuncValues(phis);

      //TEST Einheitswürfel mit verschiedenen Schnitten (Hex20)
      //      vector<double> phis (20);
      //      phis[0]=-1;
      //      phis[1]=-1;
      //      phis[2]=3;
      //      phis[3]=3;
      //      phis[4]=-1;
      //      phis[5]=-1;
      //      phis[6]=3;
      //      phis[7]=3;
      //      phis[8]=-1;
      //      phis[9]=1;
      //      phis[10]=3;
      //      phis[11]=1;
      //      phis[12]=-1;
      //      phis[13]=-1;
      //      phis[14]=3;
      //      phis[15]=3;
      //      phis[16]=-1;
      //      phis[17]=1;
      //      phis[18]=3;
      //      phis[19]=1;
      //      cell->SetGfuncValues(phis);

      //Ausgabe
      //      for (int l=0;l<8;l++)
      //      {
      //         std::cout << phis[l] << std::endl;
      //      }

      //      std::cout << "G-Funktionswerte in Zelle verpackt für Element: " << cell->Ele()->Id() << " - geschnitten?: " << cell->Bisected() << endl;
      //      for (unsigned i=0; i < mygfuncvalues.size(); i++)
      //        std::cout << "Wert am Knoten " << i << ": " << mygfuncvalues[i] << endl;
    }
    //----------------------------------------------------------------------------------------------
    // higher level of refinement
    //----------------------------------------------------------------------------------------------
    else
    {
      /* hole GID dieses Fluid Elementes
           // const int gid = ele->Id();
         für alle Eckpunkte der Zelle
           werte G-Funktion an lokalen Koordinaten des Elementes GID aus
             // GID der Elemente beider Dis identisch!
             // lokale Fluid und G-Funktion Koordinaten identisch!
           speichere den Wert am Eckpunkt
       */

      // get the vertex coordinates of the refinement cell
      // given in the coordinate system of the element == rootcell
      const std::vector<std::vector<double> > vertexcoord = cell->GetVertexCoord();
      // to be filled with the G-Function values of the refinement cell
      std::vector<double> gfunctionvalues(8);
      // G-Function values of the corresponding element
      if (cell->ReturnRootCell()==NULL)
        dserror("return NULL");
      const std::vector<double> gfuncelement = cell->ReturnRootCell()->GetGfuncValues();
      //TEST
      //std::cout<< "G-Funct of rootcell "<< std::endl;
      //for (std::size_t i=0; i<gfuncelement.size(); i++)
      //  std::cout<< gfuncelement[i] << std::endl;

      const int numnode = cell->Ele()->NumNode();
#ifdef DEBUG
      if (gfuncelement.size() != static_cast<unsigned>( numnode ))
        dserror("number of G-Function values does not match number of element nodes");
#endif
      // loop over all vertices of the refinement cell
      for (std::size_t ivertex = 0; ivertex < vertexcoord.size(); ivertex++)
      {
        // evaluate shape function at the coordinates of the current refinement cell vertex
        Epetra_SerialDenseVector funct(numnode);
        DRT::UTILS::shape_function_3D(funct,vertexcoord[ivertex][0],vertexcoord[ivertex][1],vertexcoord[ivertex][2],DRT::Element::hex8);
        // compute G-Function value
        gfunctionvalues[ivertex] = 0;
        for (int inode = 0; inode < numnode; inode++)
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
      /* für jeden Knoten des Fluid Elementes
           hole die globalen Koordinaten dieses Knotens
           finde GID des zugehörigen LevelSet Knotens
             // gleiche globale Koordinaten!
           hole Wert der G-Funktion für diese GID aus Lösungsvektor
           speichere den Wert am Eckpunkt
             // hier muss Knoten mit Eckpunkt richtig verbunden werden!
       */
    }
    else // (refinement level != 0) higher level of refinement
    {
      dserror("refinement for non-congruent discretizations not yet working!");
      /* für alle Eckpunkte der Zelle
           berechne globale Koordinaten des Eckpunktes
           finde GID des zugehörigen LevelSet Elements
             // Eckpunkt wird oft auf einer Elementkante liegen!
           transformiere die Eckpunktkoord. in lokale LevelSet Koord.
           werte G-Funktion an lokalen Koordinaten des Elementes GID aus
           speichere den Wert am Eckpunkt
       */
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
  std::map<int,std::vector<double> > intersectionpoints;
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
            //    cout << "coordinates shifted to vertex 1: " << coordinates[dim] << endl;
            //    coordinates[dim] = vertexcoord[lines[iline][0]][dim];
            //  }
            //  else if(fabs(vertexcoord[lines[iline][1]][dim]-coordinates[dim]) < 1.0E-4)
            //  {
            //    cout << "coordinates shifted to vertex 2: " << coordinates[dim] << endl;
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
        intersectionpoints[iline] = coordinates;
      }
      else
      {
        //do nothing and go to the next line
      }
    }

    break;
  }
  case DRT::Element::hex20:
  case DRT::Element::hex27:
  {
    dserror("Hallo Kilian, ich habe noch einige Hinweise!");
    // TODO @ Kilian:Kannst du, bevor du damit deine Level-Set-Beispiele
    //        rechnest noch eine paar Testfaelle dafuer rechnen. Du kannst
    //        von mir dazu ein Inputfile mit einem Hex20-Element haben.
    //        Wenn das passt, dann sollte fuer einfache Faelle auch richtig
    //        geschnitten werden und der Rest durchlaufen. Das hat zumindest
    //        fuer meinen Testfall funktioniert. Dennoch solltest du dir folgende
    //        Stellen nochmal anschauen:
    //        - buildPLC(): * alle vom distype abhängigen Stellen
    //                      * sich daraus ergebende Folgen fuer CallTetGen() und
    //                        TransformIntegrationCell()
    //        - buildFlameFrontSegments(): *bei den Faellen mit 3 bzw 4 Intersectionpoints
    //                                      waere ich mal sehr vorsichtig
    //        - darueberhinaus kann ein Blick in StoreDomainIntegrationCell() und
    //          IdentifyPolygonOrientation() nicht schaden
    //        - außerdem wuerde ich projectmidpoint in TriangulateFlameFront() erstmal rausnehmen
    //        Und noch etwas: Ich habe mal zwei Gmsh-Output-Funktionen geschrieben:
    //        RootCelltoGmsh() und FlamefronttoGmsh(). Vielleicht kannst du sie brauchen.
    //        Bei Fragen einfach vorbeikommen.
    //        Ursula

    lines = DRT::UTILS::getEleNodeNumberingLines(DRT::Element::hex20);
    // remark: vertices are assumed to be numbered in the same way the nodes are
    //         convention documented in globalreport.pdf
    /* L1:  0 1  8
     * L2:  1 2  9
     * L3:  2 3 10
     * L4:  0 3 11
     * L5:  0 4 12
     * L6:  1 5 13
     * L7:  2 6 14
     * L8:  3 7 15
     * L9:  4 5 16
     * L10: 5 6 17
     * L11: 6 7 18
     * L12: 4 7 19
     */

    //-------------------------------
    // determine intersection points
    //-------------------------------
    // loop edges of refinement cell
    for(std::size_t iline=0; iline<lines.size(); iline++)
    {
      // get G-function value of the three nodes
      double gfuncval1 = gfuncvalues[lines[iline][0]]; // left node
      //std::cout << "1: " << gfuncval1 << std::endl;
      double gfuncval2 = gfuncvalues[lines[iline][1]]; // right node
      //std::cout << "2: " << gfuncval2 << std::endl;
      double gfuncval3 = gfuncvalues[lines[iline][2]]; // node in the middle
      //std::cout << "3: " << gfuncval3 << std::endl;

      std::vector<double> coordinates1 (3);
      std::vector<double> coordinates2 (3);

      // ---------------------------------------
      // reformulation of the curve of the g-function along the edge leads
      // a*xi^2 + b*xi + c = 0
      //----------------------------------------
      double a = 0.5*gfuncval1 + 0.5*gfuncval2 - gfuncval3;
      double b = 0.5*gfuncval2 - 0.5*gfuncval1;
      double c = gfuncval3;

      // check weather the curve is really quadratic
      if (a != 0)
      {
        double determinant = b*b -4*a*c;
        //std::cout << "Determinante betraegt: " <<determinant<<" \n" ;

        // exclude complex solutions
        if (determinant >= 0.0)
        {
          // compute intersectionpoints
          double xi_1 = (-b + sqrt(determinant))/(2*a);
          double xi_2 = (-b - sqrt(determinant))/(2*a);
          //Test: Ausgabe der Schnittpunkte
          std::cout << "Schnittpunkt1: " <<xi_1<<" \n" ;
          std::cout << "Schnittpunkt2: " <<xi_2<<" \n" ;

          // check weather intersection points are inside the element, i.e. they have to be in the interval ]-1.0;1.0[
          if ((fabs(xi_1) < 1.0) or (fabs(xi_2) < 1.0))
          {
            // get coordinates
            for (int dim = 0; dim < 3; dim++)
            {
              // vertices have one coordinate in common
              if(vertexcoord[lines[iline][0]][dim] == vertexcoord[lines[iline][1]][dim])
              {
                coordinates1[dim] = vertexcoord[lines[iline][0]][dim];
                coordinates2[dim] = vertexcoord[lines[iline][0]][dim];
              }
              else // store coordinate of intersectionpoint
              {
                coordinates1[dim] = xi_1;
                coordinates2[dim] = xi_2;
              }
            }

            // at the moment the following intersection algorithm can only handle simple intersected elements
            // i.e. it based on the assumption that each edge is intersected only once
            if (fabs(xi_1)<1.0 and fabs(xi_2)<1.0 and (xi_2 != xi_1))
            {
              std::cout << "<!> WARNING: FindIntersectionPoints() has detected double intersected line of hex20-element" << std::endl;
              std::cout << "G-Function-Values of element:" << std::endl;
              for (std::size_t k=0;k<gfuncvalues.size();k++)
                std::cout << gfuncvalues[k] << std::endl;
              std::cout << "Adapt BuildFlameFrontSegements()!" << std::endl;
              dserror("Check this first!");
            }

            // store intersectionpoint if it is located in the interval ]-1.0;1.0[
            if (fabs(xi_1) < 1.0)
              intersectionpoints.insert(pair<int,std::vector<double> >(iline,coordinates1));
            if ((fabs(xi_2) < 1.0) and (xi_2 != xi_1))
              intersectionpoints.insert(pair<int,std::vector<double> >(iline,coordinates2));
          }
        }
      }
      else //linear curve of g-function along the edge
      {
        if (b != 0)
        {
          // compute intersectionpoint
          double xi_1 = -c/b;
          std::cout << "Schnittpunkt1: " <<xi_1<<" \n" ;
          // store intersectionpoint if it is located in the interval ]-1.0;1.0[
          if (fabs(xi_1) < 1.0)
          {
            // get coordinates
            for (int dim = 0; dim < 3; dim++)
            {
              // vertices have one coordinate in common
              if(vertexcoord[lines[iline][0]][dim] == vertexcoord[lines[iline][1]][dim])
              {
                coordinates1[dim] = vertexcoord[lines[iline][0]][dim];
              }
              else // store coordinates of intersectionpoint
              {
                coordinates1[dim] = xi_1;
              }
            }
            intersectionpoints.insert(pair<int,std::vector<double> >(iline,coordinates1));
          }
        }
      }
    }
    break;
  }
  default:
    dserror("FindIntersectionPoints() does not support this element shape!");
  }

  //TEST Ausgabe
  //for (std::map<int,std::vector<double> >::const_iterator iter = intersectionpoints.begin(); iter != intersectionpoints.end(); ++iter)
  //{
  //  std::cout<< iter->first << std::endl;
  //  std::vector<double> coord = iter->second;
  //  for (std::size_t isd=0; isd<3; isd++)
  //  {
  //    std::cout<< coord[isd] << std::endl;
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
void COMBUST::FlameFront::CaptureFlameFront(const Teuchos::RCP<const COMBUST::RefinementCell> rootcell)
{
  // get all refinement cells attached to this root cell (leaves in tree structure)
  // remark: if refinement is turned off, the root cell (==element) itself is returned
  std::vector<const COMBUST::RefinementCell* > RefinementCells;
  rootcell->SearchRefinementCells(RefinementCells);

  // list of domain integration cells
  GEO::DomainIntCells listDomainIntCellsperEle;
  // list of boundary integration cells
  GEO::BoundaryIntCells listBoundaryIntCellsperEle;

  // get vector of G-function values of root cell
  const std::vector<double>& gfuncvaluesrootcell = rootcell->GetGfuncValues();

  // different (spatial) XFEM integration methodologies
  switch (xfeminttype_)
  {
  //----------------------------------------------------------------------------------
  // use Constrained Delaunay Tetrahedralization (CDT) of TetGen to create tetrahedral
  // domain integration cells
  //----------------------------------------------------------------------------------
  case INPAR::COMBUST::xfemintegration_tetgen:
  {
    // loop over all refinement cells attached to root cell
    // remark: this covers the whole domain of the root cell
    for (std::size_t icell = 0; icell < RefinementCells.size(); icell++)
    {
#ifdef DEBUG
      // check if distype of cell correct (hex8)
      const DRT::Element::DiscretizationType distype = RefinementCells[icell]->Shape();
      if(distype != DRT::Element::hex8) dserror("hex8 refinement cell expected");
#endif
      // for cut cells
      if(RefinementCells[icell]->Bisected())
      {
        // - build piecewise linear complex (PLC)
        // - store triangles as boundary integration cells
        // - create tetrahedral domain integration cell using TetGen
        buildPLC(RefinementCells[icell],gfuncvaluesrootcell,listDomainIntCellsperEle,listBoundaryIntCellsperEle);
      }
      // the cell is touched (interface aligned with a cell surface)
      else if(RefinementCells[icell]->Touched())
      {
        // store hex8 domain integration cell
        StoreDomainIntegrationCell(RefinementCells[icell],listDomainIntCellsperEle);
        // store quad4 boundary integration cell
        StoreBoundaryIntegrationCell(RefinementCells[icell],listBoundaryIntCellsperEle);
      }
      else // (not bisected) and (not touched)
      {
        // store hex8 domain integration cell
        StoreDomainIntegrationCell(RefinementCells[icell],listDomainIntCellsperEle);
        // no boundary integration cell required
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
    for (std::size_t icell = 0; icell < RefinementCells.size(); icell++)
    {
#ifdef DEBUG
      // check if distype of cell correct (hex8)
      const DRT::Element::DiscretizationType distype = RefinementCells[icell]->Shape();
      if(distype != DRT::Element::hex8) dserror("hex8 refinement cell expected");
#endif

      if(RefinementCells[icell]->Bisected())
      {
        GEO::TetrahedraDecomposition decomposition(RefinementCells[icell], listBoundaryIntCellsperEle, listDomainIntCellsperEle);
      }
      // the cell is touched (interface aligned with a cell surface)
      else if(RefinementCells[icell]->Touched())
      {
        // store hex8 domain integration cell
        StoreDomainIntegrationCell(RefinementCells[icell],listDomainIntCellsperEle);
        // store quad4 boundary integration cell
        StoreBoundaryIntegrationCell(RefinementCells[icell],listBoundaryIntCellsperEle);
      }
      else // (not bisected) and (not touched)
      {
        // store hex8 domain integration cell
        StoreDomainIntegrationCell(RefinementCells[icell],listDomainIntCellsperEle);
        // no boundary integration cell required
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
      for (std::size_t icell = 0; icell < RefinementCells.size(); icell++)
      {
        // store hex8 domain integration cell
        StoreDomainIntegrationCell(RefinementCells[icell],listDomainIntCellsperEle);
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
      buildPLC(rootcell->ReturnRootCell(),gfuncvaluesrootcell,emptyvec,listBoundaryIntCellsperEle);
    }
    // the rootcell is touched (interface aligned with an element surface)
    // remark: in all other cases (leaf cell is touched) the root cell will be bisected
    //         -> previous case, since domain hexahedra are created for both cases
    if(rootcell->Touched())
    {
      // store hex8 domain integration cell
      StoreDomainIntegrationCell(rootcell->ReturnRootCell(),listDomainIntCellsperEle);
      // store quad4 boundary integration cell
      StoreBoundaryIntegrationCell(rootcell->ReturnRootCell(),listBoundaryIntCellsperEle);
    }
    else // (not bisected) and (not touched)
    {
      // store hex8 domain integration cell
      StoreDomainIntegrationCell(rootcell->ReturnRootCell(),listDomainIntCellsperEle);
      // no boundary integration cell required
    }
    break;
  }
  default: dserror("unknown type of XFEM integration");
  }

  //-----------------------------------------------------------
  // store list of integration cells per element in flame front
  //-----------------------------------------------------------
  myelementintcells_[rootcell->Ele()->Id()] = listDomainIntCellsperEle;
  // if there exist boundary integration cells
  if(listBoundaryIntCellsperEle.size() > 0)
    myboundaryintcells_[rootcell->Ele()->Id()] = listBoundaryIntCellsperEle;

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
  // TODO replace absolute tolerances by relative tolerances
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
    //cout << "iteration " << iter << ": -> |f|=" << conv << endl;

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
    cout << "Newton-Raphson algorithm for projection of midpoint did not converge" << endl;
  }
  else
  {
    converged = true;
    //LINALG::Matrix<3,1> diff(true);
    //diff(0) = midpoint[0]-projpoint(0);
    //diff(1) = midpoint[1]-projpoint(1);
    //diff(2) = midpoint[2]-projpoint(2);
    //cout << diff << endl;
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
  // TODO replace absolute tolerances by relative tolerances
  if(maxRefinementLevel_>0) dserror("This is not implemented for refinement!");
  if(polypoints.size()-1 != 4) dserror("this is not a 2D problem");

  // determine global third component (pseudo 3D direction)
  // third component of midpoint will be zero in paramenter space (without refinement strategy)
  const double midcomp = 0.0;
  int thirddim = -1;
  for (unsigned idim=0;idim<3;++idim)
  {
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
  //  std::cout<< iter << std::endl;
  //  std::vector<double> coord = pointlist[iter];
  //  for (std::size_t isd=0; isd<3; isd++)
  //  {
  //    std::cout<< coord[isd] << std::endl;
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
    if ((point1[thirddim] > -1.0-1.0E-8) and (point1[thirddim] < -1.0+1.0E-8) and
        (point2[thirddim] > -1.0-1.0E-8) and (point2[thirddim] < -1.0+1.0E-8))
    {
      // this 'ipoint' is the point
      break;
    }
  }
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
  //std::cout<<"triangles"<< std::endl;
  //for (std::size_t itriangle=0; itriangle<trianglelist.size(); itriangle++)
  //{
  //  std::cout<< "dreieck " << itriangle<< std::endl;
  //  for (int i=0; i<3; i++)
  //    std::cout<< trianglelist[itriangle][i]<<std::endl;
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
  segmentlist.insert(pair<int,std::vector<int> >(segid,segment));

  segment[0] = midpointback_id;
  segment[1] = polypoints[second_id];
  segmentlist.insert(pair<int,std::vector<int> >(segid,segment));

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
  segmentlist.insert(pair<int,std::vector<int> >(segid,segment));
  segment[0] = midpointfront_id;
  segment[1] = polypoints[forth_id];
  segmentlist.insert(pair<int,std::vector<int> >(segid,segment));

  // add segment connecting both midpoints (arbitrary ID 6)
  segment[0] = midpointback_id;
  segment[1] = midpointfront_id;
  segmentlist.insert(pair<int,std::vector<int> >(6,segment));
}

/*------------------------------------------------------------------------------------------------*
 | triangulate the interface (flame front) inside a refinement cell                   henke 10/08 |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::FlameFront::TriangulateFlameFront(
    const COMBUST::RefinementCell* cell,
    std::vector<std::vector<int> >&       trianglelist,
    std::multimap<int,std::vector<int> >& segmentlist,
    std::vector<std::vector<double> >&    pointlist,
    std::map<int,int>&                    intersectionpointsids,
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
    // std::cout << "Startsegment " << segmentiter->first << std::endl;
    // std::cout << actsegmentpoints[0] << std::endl;
    // std::cout << actsegmentpoints[1] << std::endl;

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
    //      std::cout<<"polygonvertices"<< std::endl;
    //      for (std::size_t ipoly=0; ipoly<polypoints.size(); ipoly++)
    //      {
    //        std::cout<< polypoints[ipoly] << std::endl;
    //      }

    // special case: both segments of the double intersected surface belong to the same polygon
    if (numpolygon==2)
    {
      if (polypoints.size()>(segmentlist.size()-3))
      {
        //std::cout<< "G-function values" << std::endl;
        //for (std::size_t i=0; i<gfuncvalues.size(); i++)
        //{
        //  std::cout<< gfuncvalues[i] << std::endl;
        //}
        //dserror("Unexpected intersection");
        cout << "/!\\ warning: special case of polygon occured" << endl;
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
            cout << "G-function value " << idim << " for this cell: " << gfuncvalues[idim] << endl;
          cout << "/!\\ warning === coodinate of midpoint moved to interior " << midpoint[dim] << endl;
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

#ifndef COMBUST_2D
      // add midpoint to list of points defining piecewise linear complex (interface surface)
      std::size_t midpoint_id = pointlist.size();//ids start at 0
      //std::cout<< "pointlistsize " << pointlist.size() << "midpointid " << midpoint_id<<std::endl;
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
#else
      projectMidpoint2D(trianglelist,segmentlist,pointlist,polypoints,midpoint);
#endif
    }
  }

  //TEST
  //std::cout<<"triangles"<< std::endl;
  //for (std::size_t itriangle=0; itriangle<trianglelist.size(); itriangle++)
  //{
  //  std::cout<< "dreieck " << itriangle<< std::endl;
  //  for (int i=0; i<3; i++)
  //    std::cout<< trianglelist[itriangle][i]<<std::endl;
  //}

  return;
}


/*------------------------------------------------------------------------------------------------*
 | identify the orientation of the interface polygon, normal vector + -> -                        |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::FlameFront::IdentifyPolygonOrientation(
    std::vector<int>&           segment,
    const int                   surf_id,
    std::map<int,int>&          intersectionpointsids,
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
    std::map<int,int>&                       intersectionpointsids,
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

    //std::cout << "segmentpoints size" << segmentpoints.size() << std::endl;

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
      //std::cout << "zeropoints size" << zeropoints.size() << std::endl;
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
            segmentlist.insert(pair<int,std::vector<int> >(i,zeropoints));
          }
        }
        else if (((zeropoints[0]==surfacepointlist[i][1])and(zeropoints[1]==surfacepointlist[i][3])))
        {
          if (gfuncvalues[surfacepointlist[i][0]]*gfuncvalues[surfacepointlist[i][2]]<0)
          {
            segmentlist.insert(pair<int,std::vector<int> >(i,zeropoints));
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
            segmentlist.insert(pair<int,std::vector<int> >(-1,zeropoints));
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
            segmentlist.insert(pair<int,std::vector<int> >(-1,segment));
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
      //std::cout << "one intersection point" << std::endl;
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
        for (size_t l=0; l < gfuncvalues.size(); l++)
        {
          std::cout << gfuncvalues[l] << std::endl;
        }
        dserror("can't build intersection segment");
      }

      std::vector<int> segment (2);
      segment[0] = segmentpoints[0];
      segment[1] = zeropoints[0];
      segmentlist.insert(pair<int,std::vector<int> >(i,segment));

      break;
    }
    case 2:
    {
      //store segment in segmentlist
      segmentlist.insert(pair<int,std::vector<int> >(i,segmentpoints));
      break;
    }
    case 3:
    {
      dserror("impossible number of intersectionpoints for hex8 element surface");
      break;
    }
    case 4:
    {
      //std::cout << "4 intersection points " << std::endl;
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
      //std::cout << "Maxdist" << maxdistcounter << std::endl;
      if (maxdistcounter==0 or maxdistcounter==2)
      {
        std::vector<int> segment1 (2);
        segment1[0] = segmentpoints[3];
        segment1[1] = segmentpoints[0];
        std::vector<int> segment2 (2);
        segment2[0] = segmentpoints[1];
        segment2[1] = segmentpoints[2];
        segmentlist.insert(pair<int,std::vector<int> >(i,segment1));
        segmentlist.insert(pair<int,std::vector<int> >(i,segment2));
      }
      else
      {
        std::vector<int> segment1 (2);
        segment1[0] = segmentpoints[0];
        segment1[1] = segmentpoints[1];
        std::vector<int> segment2 (2);
        segment2[0] = segmentpoints[2];
        segment2[1] = segmentpoints[3];
        segmentlist.insert(pair<int,std::vector<int> >(i,segment1));
        segmentlist.insert(pair<int,std::vector<int> >(i,segment2));
      }
      break;
    }
    default:
      dserror("unexpected number of intersectionpoints!");
    }
  }

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
  //std::cout << cell->Ele()->Id() << std::endl;
  //std::cout<< "G-Werte" << std::endl;
  //for (std::size_t i=0; i<gfuncvalues.size(); i++)
  //{
  //  std::cout<< gfuncvalues[i] << std::endl;
  //}

  //---------------------------------------------
  // fill list of points and intersection points
  //---------------------------------------------
  int numofpoints = 0;
  int numvertex = DRT::UTILS::getNumberOfElementCornerNodes(distype);
  // get intersection points from refinement cell
  const std::map<int,std::vector<double> >& intersectionpoints = cell->intersectionpoints_;

  for (int ivertex=0; ivertex<numvertex; ivertex++)
  {
    // remark: Hex hat 8 Ecken, weitere Koord in vertexcoord bei hex20-27 sind innere Knoten
    //         und keine Ecken!
    pointlist.push_back(vertexcoord[ivertex]);
    numofpoints++;
  }

  // map corresponds to intersection points, but contains the points' IDs instead of their coordinates
  // links 'pointlist' required by TetGen to 'intersectionpointlist'
  std::map<int,int> intersectionpointsids;

  // intersection points
  // map<ID of cut edge in element, coordinates of intersection point>
  for (std::map<int,std::vector<double> >::const_iterator iter = intersectionpoints.begin(); iter != intersectionpoints.end(); ++iter)
  {
    pointlist.push_back(iter->second);
    intersectionpointsids[iter->first] = numofpoints;
    numofpoints++;
  }

  //TEST
  //std::cout << "number of intersection points " << intersectionpoints.size() << endl;
  //for (std::size_t iter=0; iter<pointlist.size(); ++iter)
  //{
  //  std::cout<< iter << std::endl;
  //  std::vector<double> coord = pointlist[iter];
  //  for (std::size_t isd=0; isd<3; isd++)
  //  {
  //    std::cout<< coord[isd] << std::endl;
  //  }
  //}

  //---------------------------------------------------------
  // build segments enclosing interface patches within a cell
  //---------------------------------------------------------
  // Segments form the boundary of the interfacepatch in the cell (element, without refinement).
  // Some special cases have to be distinguished
  buildFlameFrontSegments(intersectionpointsids, segmentlist, gfuncvalues, pointlist);

  //TEST
  //std::cout<<"Segments"<< std::endl;
  //for (std::map<int,std::vector<int> >::const_iterator iter = segmentlist.begin(); iter != segmentlist.end(); ++iter)
  //{
  //  std::cout<<"Segment "<< iter->first << std::endl;
  //  std::vector<int> point = iter->second;
  //  for (std::size_t isd=0; isd<2; isd++)
  //  {
  //    std::cout<< point[isd] << std::endl;
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

#if 1

  if (xfeminttype_ == INPAR::COMBUST::xfemintegration_tetgen)
  {
    DRT::Element::DiscretizationType cell_distype = cell->Shape();

    // get vertex coordinates (local fluid element coordinates) from refinement cell
    const std::vector<std::vector<double> >& vertexcoord = cell->GetVertexCoord();

    int numnode = DRT::UTILS::getNumberOfElementNodes( cell_distype );

    std::vector<int> nids;
    nids.reserve(numnode);
    for ( int i=0; i<numnode; ++i )
    {
      nids.push_back(i);
    }

#if 0

    // Use the existing cut triangles and preserve the surface that way.

    // use element local coordinates since we want to avoid roundoff errors

    GEO::CUT::MeshIntersection intersection;

    LINALG::SerialDenseMatrix trianglecoord(3,3); // 3 directions, 3 nodes

    for ( unsigned i=0; i!=trianglelist.size(); ++i )
    {
      const std::vector<int> & t = trianglelist[i];
      for (int inode=0; inode<3; inode++)
      {
        std::copy( pointlist[t[inode]].begin(),
            pointlist[t[inode]].end(),
            &trianglecoord(0,inode) );
      }

      intersection.AddCutSide( i+1, t, trianglecoord, DRT::Element::tri3, 0 );
    }

    LINALG::SerialDenseMatrix cellcoord(3,numnode);

    for (int ivert=0; ivert<numnode; ivert++)
    {
      std::copy(vertexcoord[ivert].begin(),
          vertexcoord[ivert].end(),
          &cellcoord(0,ivert));
    }

    intersection.AddElement( 1, nids, cellcoord, cell_distype );

    intersection.Cut( true );

    GEO::CUT::ElementHandle * e = intersection.GetElement( 1 );

    if ( e!=NULL and e->IsCut() )
    {
      std::set<GEO::CUT::IntegrationCell*> cells;
      e->GetIntegrationCells( cells );

      //std::set<GEO::CUT::BoundaryCell*> bcells;
      //e->GetBoundaryCells( bcells );

      LINALG::Matrix<3,1> physCoordCorner;
      LINALG::Matrix<3,1> eleCoordDomainCorner;

      for ( std::set<GEO::CUT::IntegrationCell*>::iterator i=cells.begin(); i!=cells.end(); ++i )
      {
        GEO::CUT::IntegrationCell * ic = *i;
        DRT::Element::DiscretizationType distype = ic->Shape();
        int numnodes = DRT::UTILS::getNumberOfElementNodes( distype );

        LINALG::SerialDenseMatrix physCoord( 3, numnodes );
        LINALG::SerialDenseMatrix coord = ic->Coordinates();

        for (int ivert=0; ivert<numnodes; ivert++)
        {
          LINALG::Matrix<3,1> vertcoord;

          std::copy(&coord(0,ivert),
              &coord(0,ivert)+3,
              vertcoord.A());

          // transform vertex from local (element) coordinates to global (physical) coordinates
          GEO::elementToCurrentCoordinatesInPlace(cell_distype, xyze, vertcoord);

          std::copy(vertcoord.A(),
              vertcoord.A()+3,
              &physCoord(0,ivert));
        }

        // if degenerated don't store
        if ( distype==DRT::Element::tet4 )
          if(GEO::checkDegenerateTet(numnodes, coord, physCoord))
            continue;

        bool inGplus = GetIntCellDomainInElement(coord, gfuncvalues, DRT::Element::hex8, distype);

        domainintcelllist.push_back( GEO::DomainIntCell( distype, coord, physCoord, inGplus ) );
      }
    }
    else
    {
      dserror("cut expected");
    }
#else
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
      GEO::elementToCurrentCoordinatesInPlace(cell_distype, xyze, vertcoord);

      std::copy(vertcoord.A(),
          vertcoord.A()+3,
          &globalcellcoord(0,ivert));
    }

    // get G-function values at vertices from cell
    const std::vector<double>& gfuncvalues = cell->GetGfuncValues();

    // GEO::CUT intersection call

    GEO::CUT::LevelSetIntersection levelset;

    levelset.AddElement( 1, nids, globalcellcoord, &gfuncvalues[0], cell_distype );

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
      std::cerr << "g-function values:\n";
      std::copy( gfuncvalues.begin(), gfuncvalues.end(), std::ostream_iterator<double>( std::cerr, ", " ) );
      std::cerr << "\n";
      throw;
    }

    GEO::CUT::ElementHandle * e = levelset.GetElement( 1 );

    if ( e!=NULL and e->IsCut() )
    {
      std::set<GEO::CUT::IntegrationCell*> cells;
      e->GetIntegrationCells( cells );

      //std::set<GEO::CUT::BoundaryCell*> bcells;
      //e->GetBoundaryCells( bcells );

      LINALG::Matrix<3,1> physCoordCorner;
      LINALG::Matrix<3,1> eleCoordDomainCorner;

      for ( std::set<GEO::CUT::IntegrationCell*>::iterator i=cells.begin(); i!=cells.end(); ++i )
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
#endif
  }

#else

  //------------------------------------------------------------------------------------
  // call external program TetGen based on CDT (Constrained Delaunay Tetrahedralization)
  //------------------------------------------------------------------------------------
  // only called if XEFM integratiuon es performed using TetGen
  if (xfeminttype_ == INPAR::COMBUST::xfemintegration_tetgen)
  {
    // the node IDs of each surface are used as point IDs for TetGen
    const std::vector<std::vector<int> >& xfemsurfacelist = DRT::UTILS::getEleNodeNumberingSurfaces(distype);
#ifdef QHULL
    //------------------------------------------------------------------------------
    // create domain integration cells with TetGen for a cell/element and store them
    //------------------------------------------------------------------------------
    // remark: if cell is cut, tetrahedra have to be created serving as integration cells
    CallTetGen(pointlist, segmentlist, xfemsurfacelist, trianglelist, domainintcelllist, xyze,
        gfuncvaluesrootcell);
#else
    dserror("Set QHULL flag to use Tetgen!");
#endif
  }

#endif
  // Gmsh output for flame front
  // FlamefrontToGmsh(cell, pointlist, segmentlist, trianglelist);
}


/*------------------------------------------------------------------------------------------------*
 | store domain integration cell (only for hex8 cells at present)                     henke 08/10 |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::FlameFront::StoreDomainIntegrationCell(
    const COMBUST::RefinementCell* cell,
    GEO::DomainIntCells& domainintcelllist
)
{
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
    // if cell == element, globalcellcoord == xyze!
    // TODO we could use a check here
#endif
  }

  //-------------------------------------------
  // determine which domain the cell belongs to
  //-------------------------------------------
  // compute average G-function value for this refinement cell (= integration cell)
  bool inGplus = GetIntCellDomain(cellcoord,gfuncvalues,cell->Shape());

  //TEST
  //std::cout << "globalcellcoord " << globalcellcoord(0,3) << globalcellcoord(1,3) << std::endl;
  //if(inGplus) {
  //  std::cout << "In G plus" << std::endl;
  //}
  //else {
  //  std::cout << "In G minus" << std::endl;
  //}
  //------------------------
  // create integration cell
  //------------------------
  // create an integration cell and add it to the list of integration cells per element
  // remark: for now, this is restricted to hex8 integration cells
  domainintcelllist.push_back(GEO::DomainIntCell(DRT::Element::hex8, cellcoord, globalcellcoord, inGplus));

  return;
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
      // integration cell will always belont to the same domain.
      bool inGplus = GetIntCellDomain(cellcoord, gfuncvaluescell,celldistype);
      //store boundary integration cells in boundaryintcelllist
      boundaryintcelllist.push_back(GEO::BoundaryIntCell(DRT::Element::quad4, -1, quadcoord,
          Teuchos::null, physquadcoord,inGplus));
    }
  }
}

#if 0
#ifdef QHULL
/*------------------------------------------------------------------------------------------------*
 | calls the CDT to create tetrahedral integration cells in TetGen format             henke 10/08 |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::FlameFront::CallTetGen(
    const std::vector<std::vector<double> >&  pointlist,
    std::multimap<int,std::vector<int> >&     segmentlist,
    const std::vector<std::vector<int> >&     xfemsurfacelist,
    const std::vector<std::vector<int> >&     trianglelist,
    GEO::DomainIntCells&                      domainintcelllist,
    const LINALG::SerialDenseMatrix           xyze,
    const std::vector<double>&                gfuncvalues
)
{
  // tetgenio in is filled with the PLC
  // tetgenio out is filled by TetGen with tetrahedras

  const int dim = 3;
  tetgenio in;
  tetgenio out;
  char switches[] = "pQYY";    //- p     tetrahedralizes a PLC
  //-Q      no terminal output except errors
  // YY     do not generate additional points on surfaces -> fewer cells
  tetgenio::facet *f;
  tetgenio::polygon *p;

  //allocate point list
  in.numberofpoints = pointlist.size();
  in.pointlist = new double[in.numberofpoints * dim];

  // TODO: should scale factor be an input parameter? Is it neccessary?
  const double scalefactor = 1.0E7;

  // fill point list
  int fill = 0;
  for(int i = 0; i <  in.numberofpoints; i++)
  {
    for(int j = 0; j < dim; j++)
    {
      double coord = pointlist[i][j] * scalefactor;
      in.pointlist[fill] = (double) coord;
      fill++;
    }
  }

  in.numberoffacets = xfemsurfacelist.size() + trianglelist.size();
  in.facetlist = new tetgenio::facet[in.numberoffacets];
  in.facetmarkerlist = new int[in.numberoffacets];

  // loop over all element surfaces
  for(std::size_t i=0; i<xfemsurfacelist.size(); i++)
  {
    f = &in.facetlist[i];
    //map
    // int numsegments = 0;
    // std::map<int,std::vector<int> >::const_iterator it_segmentlist = segmentlist.find(i);
    // if (it_segmentlist!=segmentlist.end()) //check: is surface intersected?
    // {
    //   numsegments = 1;
    // }
    // falls zweimal geschnitten : numpolygons = 3
    //multimap
    int numsegments = (int) segmentlist.count(i);
    f->numberofpolygons = 1 + numsegments;
    f->polygonlist = new tetgenio::polygon[f->numberofpolygons];
    f->numberofholes = 0;
    f->holelist = NULL;
    // face of hex
    p = &f->polygonlist[0];
    p->numberofvertices = 4; //hexahedra has 4 verices
    p->vertexlist = new int[p->numberofvertices];
    for(int ivertex = 0; ivertex < p->numberofvertices; ivertex++)
    {
      p->vertexlist[ivertex] = xfemsurfacelist[i][ivertex];
    }
    //store segments if present
    if (numsegments)
    {
      //map
      // for(int j=1; j<(numsegments+1); j++) //vorsicht bei Obergrenze
      // {
      //   p = &f->polygonlist[j];
      //   p->numberofvertices = 2;
      //   p->vertexlist = new int[p->numberofvertices];
      //   for(int k=0; k<2; k++)
      //   {
      //     p->vertexlist[k] = segmentlist[i][k]; //surface i, point k
      //TEST
      //     std::cout << "Surface " << i << std::endl;
      //     std::cout << p->vertexlist[k] << std::endl;
      //   }
      // }
      //multimap
      int j=1;
      for (std::multimap<int,std::vector<int> >::iterator iter=segmentlist.equal_range(i).first; iter!=segmentlist.equal_range(i).second; iter++)
      {
        p = &f->polygonlist[j];
        p->numberofvertices = 2;
        p->vertexlist = new int[p->numberofvertices];
        std::vector<int> segmentpoints = iter->second;
        for(int k=0; k<2; k++)
        {
          p->vertexlist[k] = segmentpoints[k];
        }
        j++;
      }
    }
    //in.facetmarkerlist[i] = 0; nur für BoundaryIntCells
  }
  // store triangles
  std::size_t k = xfemsurfacelist.size();
  for(std::size_t i=0; i<trianglelist.size(); i++)
  {
    f = &in.facetlist[k];
    f->numberofpolygons = 1;
    f->polygonlist = new tetgenio::polygon[f->numberofpolygons];
    f->numberofholes = 0;
    f->holelist = NULL;
    p = &f->polygonlist[0];
    p->numberofvertices = 3; //triangle has 3 vertices
    p->vertexlist = new int[p->numberofvertices];
    //TEST
    //std::cout << "Triangle " << i << std::endl;
    //std::cout << k << std::endl;
    for(int j=0; j<3; j++)
    {
      p->vertexlist[j] = trianglelist[i][j];
      //TEST
      //std::cout << p->vertexlist[j] << std::endl;
    }
    //in.facetmarkerlist[k] = 1; nur für BoundaryIntCells
    k++;
  }

  //std::cout << "----- Start Tetgen -----" << std::endl;
  //in.save_nodes("tetin");
  //in.save_poly("tetin");
  tetrahedralize(switches, &in,  &out);
  //out.save_nodes("tetout");
  //out.save_elements("tetout");
  //out.save_faces("tetout");
  //std::cout << "----- End Tetgen -----" << std::endl;

  // store interface triangles (+ recovery of higher order meshes)
  fill = 0;
  for(int i = 0; i <  out.numberofpoints; i++)
    for(int j = 0; j < dim; j++)
    {
      out.pointlist[fill] = (double) (out.pointlist[fill] * (1.0/scalefactor));
      fill++;
    }

  TransformIntegrationCells(out, domainintcelllist, xyze, gfuncvalues);

  return;
}
#endif
#endif

#if 0
#ifdef QHULL
/*------------------------------------------------------------------------------------------------*
 | transform tetrahedral integration cells from TetGen format to BACI format          henke 10/08 |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::FlameFront::TransformIntegrationCells(
    tetgenio&                       out,
    GEO::DomainIntCells&            domainintcelllist,
    //GEO::BoundaryIntCells&          boundaryintcelllist,
    const LINALG::SerialDenseMatrix xyze,
    const std::vector<double>&      gfuncvalues
) // output: integration cells in baci format
{
  DRT::Element::DiscretizationType distype = DRT::Element::tet4;
  const int numTetNodes = DRT::UTILS::getNumberOfElementNodes(distype);

  //std::cout << "number of Tertrahedra " << out.numberoftetrahedra << std::endl;
  for(int i=0; i<out.numberoftetrahedra; i++ )
  {
    LINALG::SerialDenseMatrix tetrahedroncoord(3, numTetNodes);
    LINALG::SerialDenseMatrix phystetrahedroncoord(3, numTetNodes);
    for(int j = 0; j < numTetNodes; j++)
    {
      static LINALG::Matrix<3,1> tetcoord;
      for(int k = 0; k < 3; k++)
      {
        tetrahedroncoord(k,j) = out.pointlist[out.tetrahedronlist[i*out.numberofcorners+j]*3+k];
        // k: drei aufeinanderfolgende Einträge in der pointlist, entspricht den 3 Richtungen
        // *3: Richtungen je Knoten
        // tetraherdronlist[i*out.numberofcorners+j]: node ID for ith element, jth node
        tetcoord(k) = tetrahedroncoord(k,j);
      }
      // compute physical coordinates
      GEO::elementToCurrentCoordinatesInPlace(DRT::Element::hex8, xyze, tetcoord);
      for(int k = 0; k < 3; k++)
        phystetrahedroncoord(k,j) = tetcoord(k);
    }

    // if degenerated don't store
    if(!GEO::checkDegenerateTet(numTetNodes, tetrahedroncoord, phystetrahedroncoord))
    {
      bool inGplus = GetIntCellDomainInElement(tetrahedroncoord, gfuncvalues, DRT::Element::hex8, distype);

      domainintcelllist.push_back(GEO::DomainIntCell(distype, tetrahedroncoord, phystetrahedroncoord, inGplus));
    }
  }

  return;
}
#endif
#endif

/*------------------------------------------------------------------------------------------------*
 | compute GfuncValue for hexahedral integration cell                              rasthofer 04/10|
 *------------------------------------------------------------------------------------------------*/
bool COMBUST::FlameFront::GetIntCellDomainInElementAtCenter(
    const LINALG::SerialDenseMatrix        IntCellCoord,
    const std::vector<double>&             gfuncvalues_ele,
    const DRT::Element::DiscretizationType xfem_distype,
    const DRT::Element::DiscretizationType cell_distype
)
{
  /*------------------------------------------------------------------------------
   * - compute phi at the midpoint of domain integration cell
   * -> >= 0 in plus  'inGplus'==true
   * -> < 0  in minus 'inGplus'==false
   * ----------------------------------------------------------------------------*/

  bool inGplus = false;
  Epetra_SerialDenseVector  cellcenter(3);

  int numcellnodes = 0;
  switch(cell_distype)
  {
  //    case DRT::Element::tet4:
  //    {
  //      numcellnodes = DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::tet4>::numNodePerElement;
  //      break;
  //    }
  case DRT::Element::hex8:
  {
    numcellnodes = DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::hex8>::numNodePerElement;
    break;
  }
  default:
    dserror("Discretization Type (IntCell) not supported yet!");
  }

  int numelenodes = 0;
  if (xfem_distype==DRT::Element::hex8) {
    numelenodes = DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::hex8>::numNodePerElement;
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
    const DRT::Element::DiscretizationType cell_distype
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
  switch(cell_distype)
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
  default:
    dserror("Discretization Type (IntCell) not supported yet!");
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
    const DRT::Element::DiscretizationType cell_distype
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
  switch(cell_distype)
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
  default:
    dserror("Discretization Type (IntCell) not supported yet!");
  }

  int numelenodes = 0;
  if (xfem_distype==DRT::Element::hex8) {
    numelenodes = DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::hex8>::numNodePerElement;
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
  //  cout << "number of cell groups (cut column elements) on proc " << myrank << " " << myflamefront.size() << endl;
  //  for(std::map<int, GEO::BoundaryIntCells>::const_iterator cellgroup=myflamefront.begin(); cellgroup != myflamefront.end(); ++cellgroup)
  //    {
  //      // put ID of cut element here
  //      if (cellgroup->first == 1274)
  //      {
  //        cout << "output for element 1274 packing" << endl;
  //        const size_t numcells = (cellgroup->second).size();
  //        cout << "proc " << myrank << " number of integration cells: " << numcells << endl;
  //        for (size_t icell=0; icell<numcells; ++icell)
  //        {
  //          GEO::BoundaryIntCell cell = cellgroup->second[icell];
  //          cout << "proc " << myrank << " cell " << icell << " vertexcoord " << cell.CellNodalPosXYZ();
  //        }
  //      }
  //    }
  //#endif

#ifdef DEBUG
  std::cout << "proc " << myrank << " number of flame front pieces available before export " << myflamefront.size() << std::endl;
#endif

  DRT::PackBuffer data;
  COMBUST::FlameFront::packBoundaryIntCells(myflamefront, data);
  data.StartPacking();
  COMBUST::FlameFront::packBoundaryIntCells(myflamefront, data);

  //-----------------------------------------------------------------
  // pack data (my boundary integration cell groups) for initial send
  //-----------------------------------------------------------------
  vector<char> dataSend;
  swap( dataSend, data() );

  //-----------------------------------------------
  // send data around in a circle to all processors
  //-----------------------------------------------
  // loop over processors
  for(int num = 0; num < numproc-1; num++)
  {
    vector<int> lengthSend(1,0);
    lengthSend[0] = dataSend.size();

#ifdef DEBUG
    cout << "--- sending "<< lengthSend[0] << " bytes: from proc " << myrank << " to proc " << dest << endl;
#endif

    // send length of the data to be received ...
    MPI_Request req_length_data;
    int length_tag = 0;
    exporter.ISend(myrank, dest, &(lengthSend[0]) , size_one, length_tag, req_length_data);
    // ... and receive length
    vector<int> lengthRecv(1,0);
    exporter.Receive(source, length_tag, lengthRecv, size_one);
    exporter.Wait(req_length_data);

    // send actual data ...
    int data_tag = 4;
    MPI_Request req_data;
    exporter.ISend(myrank, dest, &(dataSend[0]), lengthSend[0], data_tag, req_data);
    // ... and receive data
    vector<char> dataRecv(lengthRecv[0]);
    exporter.ReceiveAny(source, data_tag, dataRecv, lengthRecv[0]);
    exporter.Wait(req_data);

#ifdef DEBUG
    cout << "--- receiving "<< lengthRecv[0] << " bytes: to proc " << myrank << " from proc " << source << endl;
#endif

    //-----------------------------------------------
    // unpack data (boundary integration cell groups)
    //-----------------------------------------------
    std::map<int, GEO::BoundaryIntCells> flamefront_recv;

    COMBUST::FlameFront::unpackBoundaryIntCells(dataRecv, flamefront_recv);

    //#ifdef DEBUG
    //    cout << "proc " << myrank << " receiving "<< lengthRecv[0] << " bytes from proc " << source << endl;
    //    cout << "proc " << myrank << " receiving "<< flamefront_recv.size() << " flame front pieces from proc " << source << endl;
    //
    //    for(std::map<int, GEO::BoundaryIntCells>::const_iterator cellgroup=flamefront_recv.begin(); cellgroup != flamefront_recv.end(); ++cellgroup)
    //    {
    //      // put ID of cut element here
    //      if (cellgroup->first == 1274)
    //      {
    //        cout << "output for element 1274 unpacking" << endl;
    //        const size_t numcells = (cellgroup->second).size();
    //        cout << "proc " << myrank << " number of integration cells: " << numcells << endl;
    //        for (size_t icell=0; icell<numcells; ++icell)
    //        {
    //          GEO::BoundaryIntCell cell = cellgroup->second[icell];
    //          cout << "proc " << myrank << " cell " << icell << " vertexcoord " << cell.CellNodalPosXYZ();
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
  std::cout << "proc " << myrank << " number of flame front pieces now available " << myflamefront.size() << std::endl;
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
    const vector<char>&                     dataRecv,
    std::map<int, GEO::BoundaryIntCells>&   intcellmap)
{
  // pointer to current position of group of cells in global string (counts bytes)
  vector<char>::size_type posofgroup = 0;

  while (posofgroup < dataRecv.size())
  {
    // pointer to current position in a group of cells in local string (counts bytes)
    vector<char>::size_type posingroup = 0;
    vector<char> data;
    DRT::ParObject::ExtractfromPack(posofgroup, dataRecv, data);

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
        dserror("unexpected distype");

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
    intcellmap.insert(make_pair(elegid,intcellvector));

    // check correct reading
    if (posingroup != data.size())
      dserror("mismatch in size of data %d <-> %d",(int)data.size(),posingroup);
  }
  // check correct reading
  if (posofgroup != dataRecv.size())
    dserror("mismatch in size of data %d <-> %d",(int)dataRecv.size(),posofgroup);
}

/*----------------------------------------------------------------------------------*
 |output to gmsh                                                     rasthofer 05/10|
 *----------------------------------------------------------------------------------*/
void COMBUST::FlameFront::FlamefrontToGmsh(
    const COMBUST::RefinementCell* cell,
    const std::vector<std::vector<double> >& pointlist,
    const std::multimap<int,std::vector<int> >& segmentlist,
    const std::vector<std::vector<int> >& trianglelist)
{
  const bool screen_out = false;

  const std::string filename = IO::GMSH::GetNewFileNameAndDeleteOldFiles("Flamefront", 999, 5, screen_out, 0);
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
    gmshfilecontent << "View \" " << "Intersectionpoints \" {\n";
    const std::map<int,std::vector<double> > intersectionpoints = cell->intersectionpoints_;
    for(std::map<int,std::vector<double> >::const_iterator iter = intersectionpoints.begin(); iter != intersectionpoints.end(); ++iter)
    {
      const int id = iter->first;
      const std::vector<double> point = iter->second;
      LINALG::Matrix<3,1> pos;
      for (size_t k=0; k<point.size(); k++)
        pos(k) = point[k];

      IO::GMSH::cellWithScalarToStream(DRT::Element::point1, id, pos, gmshfilecontent);
    }
    gmshfilecontent << "};\n";
  }
  {
    gmshfilecontent << "View \" " << "Segments \" {\n";
    const std::map<int,std::vector<double> > intersectionpoints = cell->intersectionpoints_;
    for(std::map<int,std::vector<int> >::const_iterator iter = segmentlist.begin(); iter != segmentlist.end(); ++iter)
    {
      const int id = iter->first;
      const std::vector<int> segmentpoints = iter->second;
      LINALG::Matrix<3,2> pos;
      for (int k=0; k<2; k++)
      {
        for (int i=0; i<3; i++)
          pos(i,k) = pointlist[segmentpoints[k]][i];
      }
      IO::GMSH::cellWithScalarToStream(DRT::Element::line2, id, pos, gmshfilecontent);
    }
    gmshfilecontent << "};\n";
  }
  {
    gmshfilecontent << "View \" " << "Tiangles \" {\n";
    const std::map<int,std::vector<double> > intersectionpoints = cell->intersectionpoints_;
    for (size_t l=0; l<trianglelist.size(); l++)
    {
      const std::vector<int> trianglepoints = trianglelist[l];
      LINALG::Matrix<3,3> pos;
      for (int k=0; k<3; k++)
      {
        for (int i=0; i<3; i++)
          pos(i,k) = pointlist[trianglepoints[k]][i];
      }
      IO::GMSH::cellWithScalarToStream(DRT::Element::tri3, 1, pos, gmshfilecontent);
    }
    gmshfilecontent << "};\n";
  }
  gmshfilecontent.close();

  return;
}

#endif // #ifdef CCADISCRET
