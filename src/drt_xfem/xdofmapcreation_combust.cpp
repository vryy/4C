/*!
\file xdofmapcreation_combust.cpp

\brief defines unknowns based on the intersection pattern from the xfem intersection

this is related to the physics of the fluid problem and therefore should not be part of the standard xfem routines

<pre>
Maintainer: Axel Gerstenberger
            gerstenberger@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>
 */
#ifdef CCADISCRET

#include <algorithm>
#include <set>
#include <iterator>


#include "xdofmapcreation.H"
#include "xdofmapcreation_combust.H"
#include "xdofmapcreation_parallel_utils.H"
#include "enrichment_utils.H"
#include "dofkey.H"
#include "../drt_lib/drt_exporter.H"
#include "../drt_lib/drt_parobject.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_combust/combust_defines.H"
#include "../drt_combust/combust3_utils.H"

/*------------------------------------------------------------------------------------------------*
 | original function: XFEM::ApplyNodalEnrichments                                     henke 03/09 |
 *------------------------------------------------------------------------------------------------*/
bool XFEM::ApplyJumpEnrichment(
    const DRT::Element*                       xfemele,
    const std::set<XFEM::PHYSICS::Field>&     fieldset,
    std::map<int, std::set<XFEM::FieldEnr> >& nodeDofMap)
{
  // type of enrichment determined by name of function; label (first argument) = 0
  const XFEM::Enrichment jumpenr(XFEM::Enrichment::typeJump, 0);

  bool skipped_element = false;

  // jump enrichments for all nodes of intersected element
  // remark: jump enrichments are added to already existing standard enrichments
  //         if jump enrichments already exist, nothing is added to set
  const int numnodes = xfemele->NumNode();
  const int* nodeidptrs = xfemele->NodeIds();
  for (int inode = 0; inode<numnodes; ++inode)
  {
    const int nodeid = nodeidptrs[inode];
    //      const bool anothervoidenrichment_in_set = XFEM::EnrichmentInNodalDofSet(node_gid, enrtype, nodalDofSet);
    //      if (not anothervoidenrichment_in_set)
    //      {
    for (std::set<XFEM::PHYSICS::Field>::const_iterator field = fieldset.begin();field != fieldset.end();++field)
    {
#ifdef COMBUST_NORMAL_ENRICHMENT
      if ((*field == XFEM::PHYSICS::Veln) or
          (*field == XFEM::PHYSICS::Pres))
      {
#endif
        nodeDofMap[nodeid].insert(XFEM::FieldEnr(*field, jumpenr));
        //std::cout << "Jump Enrichment applied for node " << nodeid << std::endl;
#ifdef COMBUST_NORMAL_ENRICHMENT
      }
#endif
    }
  };

  return skipped_element;
}


/*------------------------------------------------------------------------------------------------*
 | original function: XFEM::ApplyNodalEnrichments                                    schott 08/10 |
 | apply additional jump enrichments to touched elements                                          |
 *------------------------------------------------------------------------------------------------*/
bool XFEM::ApplyJumpEnrichmentToTouched(
    const COMBUST::InterfaceHandleCombust&    ih,
    const DRT::Element*                       xfemele,
    const std::set<XFEM::PHYSICS::Field>&     fieldset,
    std::map<int, std::set<XFEM::FieldEnr> >& nodeDofMap)
{
  // type of enrichment determined by name of function; label (first argument) = 0
  const XFEM::Enrichment jumpenr(XFEM::Enrichment::typeJump, 0);
  bool skipped_element = false;

    // jump enrichments for all nodes with Gfunc = 0.0 of a touched element
    // remark: jump enrichments are added to already existing standard enrichments
    //         if jump enrichments already exist, nothing is added to set
    const int numnodes = xfemele->NumNode();
    const int* nodeidptrs = xfemele->NodeIds();

    // get element diameter
    const size_t numnode = DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::hex8>::numNodePerElement;
    static LINALG::Matrix<3,numnode> xyze;
    GEO::fillInitialPositionArray<DRT::Element::hex8>(xfemele, xyze);
    const double hk_eleDiam = COMBUST::getEleDiameter<DRT::Element::hex8>(xyze);
    
    // get phi values for nodes
    // get pointer to vector holding G-function values at the fluid nodes
    const Teuchos::RCP<Epetra_Vector> phinp = ih.FlameFront()->Phinp();
    std::vector<double> phinp_;
    
    // extract local (element level) G-function values from global vector
    DRT::UTILS::ExtractMyNodeBasedValues(xfemele, phinp_, *phinp);
    
    if((int)phinp_.size() != numnodes) dserror("phinp_ - vector has not the same length as numnodes");
    
    unsigned nodecount = 0;
    // enrich only nodes which are touched, that means nodes which have Gfunc ~ 0.0
    for (int inode = 0; inode<numnodes; ++inode)
    {
      const int nodeid = nodeidptrs[inode];
      // to avoid phi-values near zero and bad conditioned element matrices, we enrich only
      // the nodes with phi-value approximative 0.0
      if (fabs(phinp_[inode]) < 0.1*hk_eleDiam)//1.0E-10)
      {
        for (std::set<XFEM::PHYSICS::Field>::const_iterator field = fieldset.begin();field != fieldset.end();++field)
        {
          nodeDofMap[nodeid].insert(XFEM::FieldEnr(*field, jumpenr));
          //std::cout << "additional Jump Enrichment applied for node " << nodeid << std::endl;
        }
        nodecount += 1;
      }
    }
    if (nodecount < 4)
      dserror("number of enriched nodes of touched element less than four %d ", nodecount);
    if (nodecount != 4)
      cout << "/!\\ touched element " << xfemele->Id() << " has " << nodecount << " enriched nodes" << endl;

  return skipped_element;
}


bool XFEM::ApplyKinkEnrichment(
    const DRT::Element*                       xfemele,
    const std::set<XFEM::PHYSICS::Field>&     fieldset,
    std::map<int, std::set<XFEM::FieldEnr> >& nodeDofMap)
{
  // type of enrichment determined by name of function; label (first argument) = 0
  // kann ich für mehr als ein Interfacestück je Element diese mit label = 0 bis n-1 durchnummerieren
  const XFEM::Enrichment kinkenr(XFEM::Enrichment::typeKink, 0);

  bool skipped_element = false;

    // kink enrichments for all nodes of intersected element
    // remark: already existing standard enrichments are overwritten, this might be inefficient
    const int numnodes = xfemele->NumNode();
    const int* nodeidptrs = xfemele->NodeIds();
    for (int inode = 0; inode<numnodes; ++inode)
    {
      const int nodeid = nodeidptrs[inode];
//      const bool anothervoidenrichment_in_set = XFEM::EnrichmentInNodalDofSet(node_gid, enrtype, nodalDofSet);
//      if (not anothervoidenrichment_in_set)
//      {
      // jedes Feld bekommet ein KinkEnrichment auch der Druck (im Moment)
      // überschreibt das die alten Werte
        for (std::set<XFEM::PHYSICS::Field>::const_iterator field = fieldset.begin();field != fieldset.end();++field)
        {
//          nodeDofMap[nodeid].insert(XFEM::FieldEnr(*field, kinkenr));
          //cout << "Kink Enrichment applied for all fields" << endl;
//          std::cout << "Kink Enrichment applied for node " << nodeid << std::endl;
          //Falls nur die Geschw ein KinkEnr bekommen
//          if (*field != XFEM::PHYSICS::Pres)
          nodeDofMap[nodeid].insert(XFEM::FieldEnr(*field, kinkenr));
          //cout << "Kink Enrichment applied for velocity fields only" << endl;
        }
    };

  //  }

  return skipped_element;
}


/*
 * applies additional jump enrichments to physical field Pres
 * doesn't apply additional jump enrichments to physical field Velx,Vely,Velz -> standard FEM kinks across the interface automatically
 * needed for two-phase-flow with surface tension -> discontinuous in pressure, continuous (with kinks) in velocity
 */
bool XFEM::ApplyKinkJumpEnrichmentToTouched(
    const COMBUST::InterfaceHandleCombust&    ih,
    const DRT::Element*                       xfemele,
    const std::set<XFEM::PHYSICS::Field>&     fieldset,
    std::map<int, std::set<XFEM::FieldEnr> >& nodeDofMap)
{
  // type of enrichment determined by name of function; label (first argument) = 0
  // kann ich für mehr als ein Interfacestück je Element diese mit label = 0 bis n-1 durchnummerieren
  const XFEM::Enrichment kinkenr(XFEM::Enrichment::typeKink, 0);
  const XFEM::Enrichment jumpenr(XFEM::Enrichment::typeJump, 0);

  bool skipped_element = false;

  // remark: already existing standard enrichments are overwritten, this might be inefficient
  const int numnodes = xfemele->NumNode();
  const int* nodeidptrs = xfemele->NodeIds();

  // get phi values for nodes
  // get pointer to vector holding G-function values at the fluid nodes
  const Teuchos::RCP<Epetra_Vector> phinp = ih.FlameFront()->Phinp();
  std::vector<double> phinp_;

  // extract local (element level) G-function values from global vector
  DRT::UTILS::ExtractMyNodeBasedValues(xfemele, phinp_, *phinp);

  if((int)phinp_.size() != numnodes) dserror("phinp_ - vector has not the same length as numnodes");


  // enrich only nodes which are touched, that means nodes which have Gfunc ~ 0.0
  for (int inode = 0; inode<numnodes; ++inode)
  {
    const int nodeid = nodeidptrs[inode];

    for (std::set<XFEM::PHYSICS::Field>::const_iterator field = fieldset.begin();field != fieldset.end();++field)
    {
      // nodes with abs(Gfunc)< tol * eleDiam => Gfunc is set to zero => touch point of touched face
      if (fabs(phinp_[inode]) < 1e-014)
      {
        // pressure fields gets jump enrichment
        if (*field == XFEM::PHYSICS::Pres){
          nodeDofMap[nodeid].insert(XFEM::FieldEnr(*field, jumpenr));
        }
        else if(*field == XFEM::PHYSICS::Velx ||
                *field == XFEM::PHYSICS::Vely ||
                *field == XFEM::PHYSICS::Velz)
        {
          nodeDofMap[nodeid].insert(XFEM::FieldEnr(*field, kinkenr));
        }
        else {
          dserror ("ApplyKinkJumpEnrichmentToTouched called for wrong XFEM::PHYSICS:: - Field");
        }
      }
    }
  }

  return skipped_element;
}


bool XFEM::ApplyKinkJumpEnrichment(
    const DRT::Element*                       xfemele,
    const std::set<XFEM::PHYSICS::Field>&     fieldset,
    std::map<int, std::set<XFEM::FieldEnr> >& nodeDofMap)
{
  // type of enrichment determined by name of function; label (first argument) = 0
  // kann ich für mehr als ein Interfacestück je Element diese mit label = 0 bis n-1
  // durchnummerieren
  const XFEM::Enrichment kinkenr(XFEM::Enrichment::typeKink, 0);
  const XFEM::Enrichment jumpenr(XFEM::Enrichment::typeJump, 0);

  bool skipped_element = false;

  // kink enrichments for all nodes of intersected element for velcity fields
  // remark: already existing standard enrichments are overwritten, this might be inefficient
  const int numnodes = xfemele->NumNode();
  const int* nodeidptrs = xfemele->NodeIds();
  for (int inode = 0; inode<numnodes; ++inode)
  {
    const int nodeid = nodeidptrs[inode];
    //      const bool anothervoidenrichment_in_set = XFEM::EnrichmentInNodalDofSet(node_gid, enrtype, nodalDofSet);
    //      if (not anothervoidenrichment_in_set)
    //      {
    // jedes Feld Velx,Vely,Velz bekommet ein KinkEnrichment
    // das Druck-Feld bekommt jumps
    for (std::set<XFEM::PHYSICS::Field>::const_iterator field = fieldset.begin();field != fieldset.end();++field)
    {
      // pressure fields gets jump enrichment
      if (*field == XFEM::PHYSICS::Pres){
        nodeDofMap[nodeid].insert(XFEM::FieldEnr(*field, jumpenr));
      }
      else if(*field == XFEM::PHYSICS::Velx ||
              *field == XFEM::PHYSICS::Vely ||
              *field == XFEM::PHYSICS::Velz)
      {
        nodeDofMap[nodeid].insert(XFEM::FieldEnr(*field, kinkenr));
      }
      else {
        dserror ("ApplyKinkJumpEnrichment called for wrong XFEM::PHYSICS:: - Field");
      }
    }
  };
  return skipped_element;
}


/*------------------------------------------------------------------------------------------------*
 | create a DofMap: used for combustion problems only                                 henke 03/09 |
 *------------------------------------------------------------------------------------------------*/
void XFEM::createDofMapCombust(
  const COMBUST::InterfaceHandleCombust&    ih,
  std::map<int, std::set<XFEM::FieldEnr> >& nodeDofMap,
  std::map<int, std::set<XFEM::FieldEnr> >& elementDofMap,
  const std::set<XFEM::PHYSICS::Field>&     fieldset,
  const XFEM::ElementAnsatz&                elementAnsatz,
  const Teuchos::ParameterList&             params
  )
{
  int skipped_node_enr_count = 0;
  int skipped_elem_enr_count = 0;

  // loop my column elements and add enrichments to nodes of each element
  for (int i=0; i < ih.xfemdis()->NumMyColElements(); ++i)
  {
    const DRT::Element* xfemele = ih.xfemdis()->lColElement(i);

    bool skipped_node_enr = false;
#ifdef COMBUST_STRESS_BASED
    bool skipped_elem_enr = false;
#endif

    // create an empty element ansatz map
    map<XFEM::PHYSICS::Field, DRT::Element::DiscretizationType> element_ansatz;

#ifdef COMBUST_STRESS_BASED
    // add discontinuous stress unknowns for this element, if DLM condensation is turned off
    if (not params.get<bool>("DLM_condensation"))
    {
      // ask for appropriate element ansatz (shape functions) for this type of element
      element_ansatz = elementAnsatz.getElementAnsatz(xfemele->Shape());
    }
#endif

    //-------------------------------------------------------
    // apply standard enrichments to every node in every call
    //-------------------------------------------------------
    ApplyStandardEnrichmentCombust(xfemele, fieldset, nodeDofMap);
    //std::cout << "dof initialization: element "<< xfemele->Id() << " gets standard enrichments" << std::endl;

    //--------------------------------------------------------------------
    // apply non-standard enrichments to every node, in every regular call
    //--------------------------------------------------------------------
    // remark: - for any regular call of DofManager the flame front exists (ih.FlameFront() != Teuchos::null),
    //         - for the initialization call of DofManager in CombustFluidImplicitTimeInt constructor the
    //           flame front does not exist yet (ih.FlameFront() == Teuchos::null)
    if (ih.FlameFront() != Teuchos::null)
    {
      //----------------------------------------
      // find out whether an element is bisected
      //----------------------------------------
      if (ih.ElementBisected(xfemele->Id()))
      {
        //std::cout << "Element "<< xfemele->Id() << " ist geschnitten und Knoten werden angereichert" << std::endl;
        const INPAR::COMBUST::CombustionType combusttype = DRT::INPUT::get<INPAR::COMBUST::CombustionType>(params, "combusttype");
        // build a DofMap holding dofs for all nodes including additional dofs of enriched nodes
        switch(combusttype)
        {
        case INPAR::COMBUST::combusttype_premixedcombustion:
        {
          // apply jump enrichments to all nodes of a bisected element
          skipped_node_enr = ApplyJumpEnrichment(xfemele, fieldset, nodeDofMap);

#ifdef COMBUST_STRESS_BASED
          // apply element stress enrichments to a bisected element
          ApplyElementEnrichmentCombust(xfemele, element_ansatz, elementDofMap, ih);
#endif
        }
        break;
        case INPAR::COMBUST::combusttype_twophaseflow:
        {
          // apply kink enrichments to all nodes of a bisected element
          skipped_node_enr = ApplyKinkEnrichment(xfemele, fieldset, nodeDofMap);
        }
        break;
        case INPAR::COMBUST::combusttype_twophaseflow_surf:
        {
          // apply kink enrichments to all nodes for velocity field and jumps to pressure field of a bisected element
          skipped_node_enr = ApplyKinkJumpEnrichment(xfemele, fieldset, nodeDofMap);
        }
        break;
        case INPAR::COMBUST::combusttype_twophaseflowjump:
        {
          // apply jump enrichments to all nodes of a bisected element
          skipped_node_enr = ApplyJumpEnrichment(xfemele, fieldset, nodeDofMap);
        }
        break;
        default:
          dserror("unknown type of combustion problem");
        }
        if (skipped_node_enr) skipped_node_enr_count++;
        if (skipped_node_enr) skipped_node_enr_count++;
        /*
        std::cout << "Enrichments des Elements" << std::endl;
        const int numnodes = xfemele->NumNode();
        const int* nodeidptrs = xfemele->NodeIds();
        for (int inode = 0; inode<numnodes; ++inode)
        {
          const int nodeid = nodeidptrs[inode];
          std::set<XFEM::FieldEnr> nodeenrset = nodeDofMap[nodeid];
          for (std::set<XFEM::FieldEnr>::const_iterator nodeenr = nodeenrset.begin();nodeenr != nodeenrset.end();++nodeenr)
          std::cout << "Angereichertes Feld " << physVarToString(nodeenr->getField()) << "Anreicherungstyp " << nodeenr->getEnrichment().toString() <<  std::endl;
        }
         */
      }
      //---------------------------------------
      // find out whether an element is touched
      //---------------------------------------
      else if( ih.ElementTouched(xfemele->Id()))
      {
        const INPAR::COMBUST::CombustionType combusttype = DRT::INPUT::get<INPAR::COMBUST::CombustionType>(params, "combusttype");
        // build a DofMap holding dofs for all nodes including additional dofs of enriched nodes
        switch(combusttype)
        {
        case INPAR::COMBUST::combusttype_premixedcombustion:
        {
          std::cout << "\n---  element "<< xfemele->Id() << " is touched at a face and nodes with G=0.0 get additionally enriched";
          // TODO: implementation of  additional degrees of freedom for touched elements for premixed combustion STRESS BASED!!!
#ifdef COMBUST_STRESS_BASED
          // apply element stress enrichments to an intersected element
          //skipped_elem_enr = ApplyElementEnrichmentCombust(xfemele, element_ansatz, elementDofMap, ih);
          dserror(" apply enrichments for premixedcombustion with touched elements");
#endif
#ifdef COMBUST_NITSCHE
          skipped_node_enr += ApplyJumpEnrichmentToTouched(ih, xfemele, fieldset, nodeDofMap);
#endif
        }
        break;
        case INPAR::COMBUST::combusttype_twophaseflow:
        {
          // no additional enrichments for touched elements necessary, standard FEM kinks across element faces automatically
        }
        break;
        case INPAR::COMBUST::combusttype_twophaseflow_surf:
        {
          //std::cout << "\n---  element "<< xfemele->Id() << " is touched at a face and nodes with G=0.0 get additionally enriched";
          // apply kink enrichments to all nodes for velocity field and jumps to pressure field of an intersected element
          skipped_node_enr += ApplyKinkJumpEnrichmentToTouched(ih, xfemele, fieldset, nodeDofMap);
        }
        break;
        case INPAR::COMBUST::combusttype_twophaseflowjump:
        {
          std::cout << "\n---  element "<< xfemele->Id() << " is touched at a face and nodes with G=0.0 get additionally enriched";
          // apply additional jump enrichments for pressure to all nodes with Gfunc = 0.0 of a touched element
          skipped_node_enr += ApplyJumpEnrichmentToTouched(ih, xfemele, fieldset, nodeDofMap);
        }
        break;
        default:
          dserror("unknown type of combustion problem");
        }
      }
      else
      {
        // nothing to do, standard element
      }
    }
  }

#ifdef PARALLEL
  syncNodalDofs(ih, nodeDofMap);
#endif

  cout << " skipped "<< skipped_node_enr_count << " node unknowns" << endl;
  cout << " skipped "<< skipped_elem_enr_count << " element unknowns" << endl;
}

/*------------------------------------------------------------------------------------------------*
 | original function: XFEM::ApplyStandardEnrichment                                   henke 03/09 |
 *------------------------------------------------------------------------------------------------*/
void XFEM::ApplyStandardEnrichmentCombust(
    const DRT::Element*                       xfemele,
    const std::set<XFEM::PHYSICS::Field>&     fieldset,
    std::map<int, std::set<XFEM::FieldEnr> >& nodeDofMap
    )
{
  // type of enrichment corresponds to name of function; label (second argument) = 0
  const XFEM::Enrichment stdenr(XFEM::Enrichment::typeStandard, 0);

  // standard enrichments for all nodes of element
  // remark: standard enrichments are added to already existing jump enrichments;
  //         if standard enrichments already exist, nothing is added to set
  const int numnodes = xfemele->NumNode();
  const int* nodeidptrs = xfemele->NodeIds();
  for (int inode = 0; inode<numnodes; ++inode)
  {
    const int nodeid = nodeidptrs[inode];
    //std::cout << "element " << xfemele->Id() << " node ID: " << nodeid << endl;
    for (std::set<XFEM::PHYSICS::Field>::const_iterator fields = fieldset.begin();fields != fieldset.end();++fields)
    {
#ifdef COMBUST_NORMAL_ENRICHMENT
      if (*fields != XFEM::PHYSICS::Veln)
      {
#endif
      nodeDofMap[nodeid].insert(XFEM::FieldEnr(*fields, stdenr));
      //std::cout << "Standard Enrichment applied for node " << nodeid << std::endl;
#ifdef COMBUST_NORMAL_ENRICHMENT
      }
#endif
    }
  };
}


void XFEM::ApplyElementEnrichmentCombust(
    const DRT::Element*                                                xfemele,
    const map<XFEM::PHYSICS::Field, DRT::Element::DiscretizationType>& element_ansatz,
    std::map<int, std::set<XFEM::FieldEnr> >&                          elementDofMap,
    const COMBUST::InterfaceHandleCombust&                             ih
)
{
  // type of enrichment corresponds to name of function; label (second argument) = 0
  const XFEM::Enrichment elementenr1(XFEM::Enrichment::typeStandard, 0);
  const XFEM::Enrichment elementenr2(XFEM::Enrichment::typeJump, 0);

  map<XFEM::PHYSICS::Field, DRT::Element::DiscretizationType>::const_iterator fielditer;
  for (fielditer = element_ansatz.begin();fielditer != element_ansatz.end();++fielditer)
  {
    elementDofMap[xfemele->Id()].insert(XFEM::FieldEnr(fielditer->first, elementenr1));
    elementDofMap[xfemele->Id()].insert(XFEM::FieldEnr(fielditer->first, elementenr2));
    cout << "added element enrichment on crearedofmapCombust" << endl;
  }
}

void XFEM::ApplyElementEnrichmentCombust(
    const DRT::Element*                                                xfemele,
    const map<XFEM::PHYSICS::Field, DRT::Element::DiscretizationType>& element_ansatz,
    std::set<XFEM::FieldEnr>&                                          elementFieldEnrSet,
    const COMBUST::InterfaceHandleCombust&                             ih
)
{
  // type of enrichment corresponds to name of function; label (second argument) = 0
  const XFEM::Enrichment elementenr1(XFEM::Enrichment::typeStandard, 0);
  const XFEM::Enrichment elementenr2(XFEM::Enrichment::typeVoid, 0);

  map<XFEM::PHYSICS::Field, DRT::Element::DiscretizationType>::const_iterator fielditer;
  for (fielditer = element_ansatz.begin();fielditer != element_ansatz.end();++fielditer)
  {
    elementFieldEnrSet.insert(XFEM::FieldEnr(fielditer->first, elementenr1));
    elementFieldEnrSet.insert(XFEM::FieldEnr(fielditer->first, elementenr2));
  }
}


#endif // #ifdef CCADISCRET
