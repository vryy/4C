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
 | create a DofMap: used for combustion problems only; derived from createDofMap()    henke 03/09 |
 *------------------------------------------------------------------------------------------------*/
void XFEM::createDofMapCombust(
  const COMBUST::InterfaceHandleCombust&    ih,
  std::map<int, std::set<XFEM::FieldEnr> >& nodeDofMap,
  std::map<int, std::set<XFEM::FieldEnr> >& elementDofMap,
  const std::set<XFEM::PHYSICS::Field>&     fieldset,
  const Teuchos::ParameterList&             params
  )
{
  //if (fluiddis->Comm().MyPID() == 0)
  //  std::cout << "Creating DofMap for combustion problem" << std::endl;

  const double volumeRatioLimit = params.get<double>("volumeRatioLimit");

  bool skipped_node_enr = false;

  // loop my column elements and add enrichments to nodes of each element
  for (int i=0; i < ih.xfemdis()->NumMyColElements(); ++i)
  {
    const DRT::Element* xfemele = ih.xfemdis()->lColElement(i);

    /* the following procedure can be abbreviated as soon as the interfacehandle holds domain
     * integration cells. Then the following line will do the job:
     *
     * if (ih.ElementIntersected(xfemele->Id()))
     *
     * Since ElementIntersected() is a member function of the base class InterfaceHandle, then not
     * the specific InterfaceHandleCombust will have to be given to createDofMap(), but a base class
     * instance will do. The following "if" statement is the only reason a InterfaceHandleCombust is
     * needed in here!
     *
     * henke 03/09 */

    // initialization call of DofManager in CombustFluidImplicitTimeInt constructor will end up here
    if (ih.FlameFront() == Teuchos::null)
    {
      // std::cout << "dof initialization: element "<< xfemele->Id() << " gets standard enrichments (= 4 dofs)" << std::endl;
      ApplyStandardEnrichmentCombust(xfemele, fieldset, nodeDofMap);
    }
    // any regular call of DofManager (existing flame front) will end up here
    else
    {
      // apply standard enrichment to every node
      ApplyStandardEnrichmentCombust(xfemele, fieldset, nodeDofMap);

/*
  Check der FlameFrontPatches für das Element ist nur vorübergehende Lösung
  // get the refinement cell belonging to a fluid element
  std::map<int,const Teuchos::RCP<const COMBUST::RefinementCell> >::const_iterator iter = ih.FlameFront()->FlameFrontPatches().find(xfemele->Id());
  if (iter->second->Intersected())
*/
      //---------------------------------------------------
      // find out whether an element is intersected or not
      // by looking at number of integration cells
      //---------------------------------------------------
      // get vector of integration cells for this element
      std::size_t numcells= ih.GetNumDomainIntCells(xfemele);
      if (numcells > 1) // element intersected
      {
        const INPAR::COMBUST::CombustionType combusttype = params.get<INPAR::COMBUST::CombustionType>("combusttype");
        // build a DofMap holding dofs for all nodes including additional dofs of enriched nodes
        switch(combusttype)
        {
          case INPAR::COMBUST::combusttype_premixedcombustion:
          {
            // apply jump enrichments to all nodes of an intersected element
            skipped_node_enr = ApplyJumpEnrichment(xfemele, fieldset, volumeRatioLimit, nodeDofMap);
            // remark: Brauche ich hier tatsächlich einen Rückgabewert, nur um zu checken, ob die
            //         Anreicherung funktioniert hat?
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
//        std::cout << "Element "<< xfemele->Id() << " ist geschnitten und Knoten werden angereichert" << std::endl;
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
        // in case there are enriched element dofs they have to be applied now   henke 04/09
//        ApplyElementEnrichments();
      }
      else if (numcells == 1) // element not intersected
      {
/*
        std::cout << "Element "<< xfemele->Id() << " ist nicht geschnitten" << std::endl;
        // apply standard enrichments to all nodes of an intersected element
        ApplyStandardEnrichmentCombust(xfemele, fieldset, nodeDofMap);
        //
        //  Funktion sollte checken, ob schon JumpEnrichments (andere Enrichments) vorliegen. Sollen
        //  diese überschrieben werden? Kläre die Theorie dazu!
        //
*/
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

}

/*------------------------------------------------------------------------------------------------*
 | original function: XFEM::ApplyStandardEnrichment                                   henke 03/09 |
 *------------------------------------------------------------------------------------------------*/
void XFEM::ApplyStandardEnrichmentCombust(
    const DRT::Element*                           xfemele,
    const std::set<XFEM::PHYSICS::Field>&         fieldset,
    std::map<int, std::set<XFEM::FieldEnr> >&     nodalDofSet)
{
  // type of enrichment determined by name of function; label (first argument) = 0
  const XFEM::Enrichment stdenr(XFEM::Enrichment::typeStandard, 0);

  // standard enrichments for all nodes of element
  // remark: standard enrichments are added to already existing jump enrichments;
  //         if standard enrichments already exist, nothing is added to set
  const int numnodes = xfemele->NumNode();
  const int* nodeidptrs = xfemele->NodeIds();
  for (int inode = 0; inode<numnodes; ++inode)
  {
    const int nodeid = nodeidptrs[inode];
//    std::cout << "element " << xfemele->Id() << " node ID: " << nodeid << endl;
    for (std::set<XFEM::PHYSICS::Field>::const_iterator field = fieldset.begin();field != fieldset.end();++field)
    {
      nodalDofSet[nodeid].insert(XFEM::FieldEnr(*field, stdenr));
//      std::cout << "Standard Enrichment applied for node " << nodeid << std::endl;
    }
  };
}




#endif  // #ifdef CCADISCRET
