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

/*------------------------------------------------------------------------------------------------*
 | original function: XFEM::ApplyNodalEnrichments                                     henke 03/09 |
 *------------------------------------------------------------------------------------------------*/
bool XFEM::ApplyJumpEnrichment(
    const DRT::Element*                       xfemele,
    const std::set<XFEM::PHYSICS::Field>&     fieldset,
    const double                              volumeRatioLimit,
    std::map<int, std::set<XFEM::FieldEnr> >& nodeDofMap)
{
  // type of enrichment determined by name of function; label (first argument) = 0
  const XFEM::Enrichment jumpenr(XFEM::Enrichment::typeJump, 0);

//  const double volumeratio = XFEM::DomainCoverageRatio(*xfemele,ih);
//  const bool almost_empty_element = (fabs(1.0-volumeratio) < volumeRatioLimit);

  bool skipped_element = false;

//  if ( not almost_empty_element)
//  {

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
          nodeDofMap[nodeid].insert(XFEM::FieldEnr(*field, jumpenr));
//          std::cout << "Jump Enrichment applied for node " << nodeid << std::endl;
        }
//      }
    };

  //  }

  /* almost_empty_element: diesen Fall schließen wir erstmal aus!   henke 03/09
  else
  { // void enrichments only in the fluid domain
    const int nen = xfemele->NumNode();
    const int* nodeidptrs = xfemele->NodeIds();
    for (int inen = 0; inen<nen; ++inen)
    {
      const int node_gid = nodeidptrs[inen];
      const LINALG::Matrix<3,1> nodalpos(ih.xfemdis()->gNode(node_gid)->X());
      const int label = ih.PositionWithinConditionNP(nodalpos);
      const bool in_fluid = (0 == label);

      if (in_fluid)
      {
        const bool anothervoidenrichment_in_set = EnrichmentInNodalDofSet(node_gid, enrtype, nodalDofSet);
        if (not anothervoidenrichment_in_set)
        {
          for (std::set<XFEM::PHYSICS::Field>::const_iterator field = fieldset.begin();field != fieldset.end();++field)
          {
            nodalDofSet[node_gid].insert(XFEM::FieldEnr(*field, voidenr));
          }
        }
      }
    };
    skipped_element = true;
//    cout << "skipped interior void unknowns for element: "<< xfemele->Id() << ", volumeratio limit: " << std::scientific << volumeRatioLimit << ", volumeratio: abs (" << std::scientific << (1.0 - volumeratio) << " )" << endl;
  }*/

  return skipped_element;
}

bool XFEM::ApplyKinkEnrichment(
    const DRT::Element*                       xfemele,
    const std::set<XFEM::PHYSICS::Field>&     fieldset,
    const double                              volumeRatioLimit,
    std::map<int, std::set<XFEM::FieldEnr> >& nodeDofMap)
{
  // type of enrichment determined by name of function; label (first argument) = 0
  // kann ich für mehr als ein Interfacestück je Element diese mit label = 0 bis n-1
  // durchnummerieren
  const XFEM::Enrichment kinkenr(XFEM::Enrichment::typeKink, 0);

//  const double volumeratio = XFEM::DomainCoverageRatio(*xfemele,ih);
//  const bool almost_empty_element = (fabs(1.0-volumeratio) < volumeRatioLimit);

  bool skipped_element = false;

//  if ( not almost_empty_element)
//  {

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
//      }
    };

  //  }

  /* almost_empty_element: diesen Fall schließen wir erstmal aus!   henke 03/09
  else
  { // void enrichments only in the fluid domain
    const int nen = xfemele->NumNode();
    const int* nodeidptrs = xfemele->NodeIds();
    for (int inen = 0; inen<nen; ++inen)
    {
      const int node_gid = nodeidptrs[inen];
      const LINALG::Matrix<3,1> nodalpos(ih.xfemdis()->gNode(node_gid)->X());
      const int label = ih.PositionWithinConditionNP(nodalpos);
      const bool in_fluid = (0 == label);

      if (in_fluid)
      {
        const bool anothervoidenrichment_in_set = EnrichmentInNodalDofSet(node_gid, enrtype, nodalDofSet);
        if (not anothervoidenrichment_in_set)
        {
          for (std::set<XFEM::PHYSICS::Field>::const_iterator field = fieldset.begin();field != fieldset.end();++field)
          {
            nodalDofSet[node_gid].insert(XFEM::FieldEnr(*field, voidenr));
          }
        }
      }
    };
    skipped_element = true;
//    cout << "skipped interior void unknowns for element: "<< xfemele->Id() << ", volumeratio limit: " << std::scientific << volumeRatioLimit << ", volumeratio: abs (" << std::scientific << (1.0 - volumeratio) << " )" << endl;
  }*/

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
  const double volumeRatioLimit = params.get<double>("volumeRatioLimit");
  const double boundaryRatioLimit = params.get<double>("boundaryRatioLimit");

  int skipped_node_enr_count = 0;
  int skipped_elem_enr_count = 0;

  // loop my column elements and add enrichments to nodes of each element
  for (int i=0; i < ih.xfemdis()->NumMyColElements(); ++i)
  {
    const DRT::Element* xfemele = ih.xfemdis()->lColElement(i);

    bool skipped_node_enr = false;
    bool skipped_elem_enr = false;

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
    //std::cout << "dof initialization: element "<< xfemele->Id() << " gets standard enrichments (= 4 dofs)" << std::endl;

    //--------------------------------------------------------------------
    // apply non-standard enrichments to every node, in every regular call
    //--------------------------------------------------------------------
    // remark: - for any regular call of DofManager the flame front exists (ih.FlameFront() != Teuchos::null),
    //         - for the initialization call of DofManager in CombustFluidImplicitTimeInt constructor the
    //           flame front does not exist yet (ih.FlameFront() == Teuchos::null)
    if (ih.FlameFront() != Teuchos::null)
    {
      //--------------------------------------------------------------------------------------------
      // find out whether an element is intersected or not by looking at number of integration cells
      //--------------------------------------------------------------------------------------------
      // get vector of integration cells for this element
      std::size_t numcells= ih.GetNumDomainIntCells(xfemele);
      if (numcells > 1) // element intersected
      {
        //std::cout << "Element "<< xfemele->Id() << " ist geschnitten und Knoten werden angereichert" << std::endl;
        const INPAR::COMBUST::CombustionType combusttype = params.get<INPAR::COMBUST::CombustionType>("combusttype");
        // build a DofMap holding dofs for all nodes including additional dofs of enriched nodes
        switch(combusttype)
        {
          case INPAR::COMBUST::combusttype_premixedcombustion:
          {
            // apply jump enrichments to all nodes of an intersected element
            skipped_node_enr = ApplyJumpEnrichment(xfemele, fieldset, volumeRatioLimit, nodeDofMap);

#ifdef COMBUST_STRESS_BASED
            // apply element stress enrichments to an intersected element
            skipped_elem_enr = ApplyElementEnrichmentCombust(xfemele, element_ansatz, elementDofMap, ih, boundaryRatioLimit);
#endif
          }
          break;
          case INPAR::COMBUST::combusttype_twophaseflow:
          {
            // apply kink enrichments to all nodes of an intersected element
            skipped_node_enr = ApplyKinkEnrichment(xfemele, fieldset, volumeRatioLimit, nodeDofMap);
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
      else if (numcells == 1) // element not intersected
      {
        // nothing to do
        //std::cout << "Element "<< xfemele->Id() << " ist nicht geschnitten" << std::endl;
      }
      else // numcells == 0 or negative number
      {
        // no integration cell -> impossible, something went wrong!
        dserror ("There are no DomainIntCells for element %d ", xfemele->Id());
      }
    }
  }

#ifdef PARALLEL
  syncNodalDofs(ih, nodeDofMap);
#endif

  cout << " skipped "<< skipped_node_enr_count << " node unknowns (volumeratio limit:   " << std::scientific << volumeRatioLimit   << ")" << endl;
  cout << " skipped "<< skipped_elem_enr_count << " element unknowns (boundaryratio limit: " << std::scientific << boundaryRatioLimit << ")" << endl;
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
      nodeDofMap[nodeid].insert(XFEM::FieldEnr(*fields, stdenr));
      //std::cout << "Standard Enrichment applied for node " << nodeid << std::endl;
    }
  };
}


bool XFEM::ApplyElementEnrichmentCombust(
    const DRT::Element*                                                xfemele,
    const map<XFEM::PHYSICS::Field, DRT::Element::DiscretizationType>& element_ansatz,
    std::map<int, std::set<XFEM::FieldEnr> >&                          elementDofMap,
    const COMBUST::InterfaceHandleCombust&                             ih,
    const double                                                       boundaryRatioLimit
)
{
  // check, how much area for integration we have (from BoundaryIntcells)
  const double boundarysize = XFEM::BoundaryCoverageRatio(*xfemele,ih.GetBoundaryIntCells(xfemele->Id()),ih);
  const bool almost_zero_surface = (fabs(boundarysize) < boundaryRatioLimit);

  bool skipped_element = false;

  // type of enrichment corresponds to name of function; label (second argument) = 0
  const XFEM::Enrichment elementenr(XFEM::Enrichment::typeVoid, 0);

  //if (not almost_zero_surface)
  //{
      map<XFEM::PHYSICS::Field, DRT::Element::DiscretizationType>::const_iterator fielditer;
      for (fielditer = element_ansatz.begin();fielditer != element_ansatz.end();++fielditer)
      {
        elementDofMap[xfemele->Id()].insert(XFEM::FieldEnr(fielditer->first, elementenr));
      }
  //}
  //else
  //{
  //  skipped_element = true;
    //cout << "skipped stress unknowns for element: "<< xfemele->Id() << ", boundary size: " << boundarysize << endl;
  //}
  return skipped_element;
}

bool XFEM::ApplyElementEnrichmentCombust(
    const DRT::Element*                                                xfemele,
    const map<XFEM::PHYSICS::Field, DRT::Element::DiscretizationType>& element_ansatz,
    std::set<XFEM::FieldEnr>&                                          elementFieldEnrSet,
    const COMBUST::InterfaceHandleCombust&                             ih,
    const double                                                       boundaryRatioLimit
)
{
  // check, how much area for integration we have (from BoundaryIntcells)
  const double boundarysize = XFEM::BoundaryCoverageRatio(*xfemele,ih.GetBoundaryIntCells(xfemele->Id()),ih);
  const bool almost_zero_surface = (fabs(boundarysize) < boundaryRatioLimit);

  bool skipped_element = false;

  // type of enrichment corresponds to name of function; label (second argument) = 0
  const XFEM::Enrichment elementenr(XFEM::Enrichment::typeVoid, 0);

  //if (not almost_zero_surface)
  //{
      map<XFEM::PHYSICS::Field, DRT::Element::DiscretizationType>::const_iterator fielditer;
      for (fielditer = element_ansatz.begin();fielditer != element_ansatz.end();++fielditer)
      {
        elementFieldEnrSet.insert(XFEM::FieldEnr(fielditer->first, elementenr));
      }
  //}
  //else
  //{
  //  skipped_element = true;
    //cout << "skipped stress unknowns for element: "<< xfemele->Id() << ", boundary size: " << boundarysize << endl;
  //}
  return skipped_element;
}


#endif // #ifdef CCADISCRET
