/*!
\file xdofmapcreation.cpp

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
#include "xdofmapcreation_fsi.H"
#include "xdofmapcreation_parallel_utils.H"
#include "enrichment_utils.H"
#include "dofkey.H"
#include "../drt_lib/drt_exporter.H"
#include "../drt_lib/drt_parobject.H"
#include "../drt_lib/drt_utils.H"




/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
static bool ApplyElementEnrichments(
    const DRT::Element*                           xfemele,
    const map<XFEM::PHYSICS::Field, DRT::Element::DiscretizationType>&  element_ansatz,
    const XFEM::InterfaceHandleXFSI&              ih,
    const int&                                    label,
    const XFEM::Enrichment::EnrType               enrtype,
    const double                                  boundaryRatioLimit,
    std::set<XFEM::FieldEnr>&                     enrfieldset)
{
  // check, how much area for integration we have (from BoundaryIntcells)
  const double boundarysize = XFEM::BoundaryCoverageRatio(*xfemele,ih.GetBoundaryIntCells(xfemele->Id()),ih);
  const bool almost_zero_surface = (fabs(boundarysize) < boundaryRatioLimit);

  bool skipped_element = false;

  const XFEM::Enrichment enr(enrtype, label);
  if ( not almost_zero_surface)
  {
      map<XFEM::PHYSICS::Field, DRT::Element::DiscretizationType>::const_iterator fielditer;
      for (fielditer = element_ansatz.begin();fielditer != element_ansatz.end();++fielditer)
      {
        enrfieldset.insert(XFEM::FieldEnr(fielditer->first, enr));
      }
  }
  else
  {
    skipped_element = true;
//    cout << "skipped stress unknowns for element: "<< xfemele->Id() << ", boundary size: " << boundarysize << endl;
  }
  return skipped_element;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
static bool ApplyNodalEnrichmentsVoidFSI(
    const DRT::Element*                           xfemele,
    const XFEM::InterfaceHandleXFSI&              ih,
    const XFEM::Enrichment::EnrType               enrtype,
    const std::set<XFEM::PHYSICS::Field>&         fieldset,
    const double                                  volumeRatioLimit,
    std::map<int, std::set<XFEM::FieldEnr> >&     nodalDofSet)
{


  set<int> beles; // all bele that intersect this xfem element
  if (ih.ElementIntersected(xfemele->Id()))
    beles = ih.GetIntersectingBoundaryElementsGID(xfemele->Id());

  map<int, set<int> > elementsByLabel;
  map<int, set<int> > nodesByLabel;
  XFEM::separateByLabel(ih, beles, elementsByLabel, nodesByLabel);

  if (nodesByLabel.size() > 2)
    dserror("needs more generalizashun!");

  bool skipped_element = false;

  if (nodesByLabel.size() == 2)
  {
    //check, whether these patches are connected
    map<int, set<int> >::const_iterator entry1 = nodesByLabel.begin();
    map<int, set<int> >::const_iterator entry2 = entry1++;

//    cout << "(";
//    for (set<int>::const_iterator i = entry1->second.begin(); i != entry1->second.end(); ++i)
//    {
//      if (i != entry1->second.begin())
//        cout << ",";
//      cout << *i ;
//    }
//    cout << ")" << endl;
//
//    cout << "(";
//    for (set<int>::const_iterator i = entry2->second.begin(); i != entry2->second.end(); ++i)
//    {
//      if (i != entry2->second.begin())
//        cout << ",";
//      cout << *i;
//    }
//    cout << ")" << endl;

    const bool connected = XFEM::ConnectedElements(entry1->second, entry2->second);
    if (connected)
    {
      cout << "Nodal VoidFSI: connected patch" << endl;
      const double volumeratio = XFEM::DomainCoverageRatio(*xfemele,ih);
      const bool almost_empty_element = (fabs(1.0-volumeratio) < volumeRatioLimit);

      const XFEM::Enrichment voidenr(XFEM::Enrichment::typeVoidFSI, -1);

      bool skipped_element = false;

      if ( not almost_empty_element)
      { // void enrichments for everybody !!!
        const int nen = xfemele->NumNode();
        const int* nodeidptrs = xfemele->NodeIds();
        for (int inen = 0; inen<nen; ++inen)
        {
          for (std::set<XFEM::PHYSICS::Field>::const_iterator field = fieldset.begin();field != fieldset.end();++field)
          {
            nodalDofSet[nodeidptrs[inen]].insert(XFEM::FieldEnr(*field, voidenr));
          }
        };
      }
      else
      { // void enrichments only in the fluid domain
        const int nen = xfemele->NumNode();
        const int* nodeidptrs = xfemele->NodeIds();
        for (int inen = 0; inen<nen; ++inen)
        {
          const LINALG::Matrix<3,1> nodalpos(xfemele->Nodes()[inen]->X());
          const bool in_fluid = (0 == ih.PositionWithinConditionNP(nodalpos));

          if (in_fluid)
          {
            for (std::set<XFEM::PHYSICS::Field>::const_iterator field = fieldset.begin();field != fieldset.end();++field)
            {
              nodalDofSet[nodeidptrs[inen]].insert(XFEM::FieldEnr(*field, voidenr));
            }
          }
        };
        skipped_element = true;
        //    cout << "skipped interior void unknowns for element: "<< xfemele->Id() << ", volumeratio limit: " << std::scientific << volumeRatioLimit << ", volumeratio: abs (" << std::scientific << (1.0 - volumeratio) << " )" << endl;
      }
    }
  }
  return skipped_element;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
static bool ApplyElementEnrichmentsVoidFSI(
    const DRT::Element*                           xfemele,
    const map<XFEM::PHYSICS::Field, DRT::Element::DiscretizationType>&  element_ansatz,
    const XFEM::InterfaceHandleXFSI&              ih,
    const double                                  boundaryRatioLimit,
    std::set<XFEM::FieldEnr>&                     enrfieldset)
{

  set<int> beles; // alle bele that intersect this xfem element
  if (ih.ElementIntersected(xfemele->Id()))
    beles = ih.GetIntersectingBoundaryElementsGID(xfemele->Id());

  map<int, set<int> > elementsByLabel;
  map<int, set<int> > nodesByLabel;
  XFEM::separateByLabel(ih, beles, elementsByLabel, nodesByLabel);

  if (nodesByLabel.size() > 2)
    dserror("needs more generalizashun!");

  bool skipped_element = false;

  if (nodesByLabel.size() == 2)
  {
    //check, whether these patches are connected
    map<int, set<int> >::const_iterator entry1 = nodesByLabel.begin();
    map<int, set<int> >::const_iterator entry2 = entry1++;
    const bool connected = XFEM::ConnectedElements(entry1->second, entry2->second);
    if (connected)
    {
      cout << "Elemental VoidFSI: connected patch" << endl;
      skipped_element = ApplyElementEnrichments(xfemele, element_ansatz, ih, -1, XFEM::Enrichment::typeVoidFSI,  boundaryRatioLimit, enrfieldset);
    }
    else
    {
      // single enrichments will be done in the Void enrichment processs
    }
  }
  else
  {
    skipped_element = false;
  }
  return skipped_element;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
static bool ApplyNodalEnrichmentsVoid(
    const DRT::Element*                           xfemele,
    const XFEM::InterfaceHandle&                  ih,
    const int&                                    label,
    const XFEM::Enrichment::EnrType               enrtype,
    const std::set<XFEM::PHYSICS::Field>&         fieldset,
    const double                                  volumeRatioLimit,
    std::map<int, std::set<XFEM::FieldEnr> >&     nodalDofSet)
{
  const double volumeratio = XFEM::DomainCoverageRatio(*xfemele,ih);
  const bool almost_empty_element = (fabs(1.0-volumeratio) < volumeRatioLimit);

  const XFEM::Enrichment voidenr(enrtype, label);

  bool skipped_element = false;

  if ( not almost_empty_element)
  { // void enrichments for everybody !!!
    const int nen = xfemele->NumNode();
    const int* nodeidptrs = xfemele->NodeIds();
    for (int inen = 0; inen<nen; ++inen)
    {
      const int node_gid = nodeidptrs[inen];
      const bool VoidFSIenrichment_in_set = EnrichmentInDofSet(XFEM::Enrichment::typeVoidFSI, nodalDofSet[node_gid]);
      if (not VoidFSIenrichment_in_set)
      {
        for (std::set<XFEM::PHYSICS::Field>::const_iterator field = fieldset.begin();field != fieldset.end();++field)
        {
          nodalDofSet[node_gid].insert(XFEM::FieldEnr(*field, voidenr));
        }
      }
    };
  }
  else
  { // void enrichments only in the fluid domain
    const int nen = xfemele->NumNode();
    const int* nodeidptrs = xfemele->NodeIds();
    for (int inen = 0; inen<nen; ++inen)
    {
      const int node_gid = nodeidptrs[inen];
      const bool VoidFSIenrichment_in_set = EnrichmentInDofSet(XFEM::Enrichment::typeVoidFSI, nodalDofSet[node_gid]);
      if (not VoidFSIenrichment_in_set)
      {
        const LINALG::Matrix<3,1> nodalpos(xfemele->Nodes()[inen]->X());
        const bool in_fluid = (0 == ih.PositionWithinConditionNP(nodalpos));

        if (in_fluid)
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
  }
  return skipped_element;
}



#if 0
/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
static bool ApplyNodalEnrichmentsNodeWise(
    const DRT::Element*                           xfemele,
    const XFEM::InterfaceHandle&                  ih,
    const int&                                    label,
    const XFEM::Enrichment::EnrType               enrtype,
    const std::set<XFEM::PHYSICS::Field>&         fieldset,
    const double                                  volumeRatioLimit,
    std::map<int, std::set<XFEM::FieldEnr> >&     nodalDofSet)
{
  const vector<double> ratios = XFEM::DomainCoverageRatioPerNode(*xfemele,ih);

  const XFEM::Enrichment voidenr(label, XFEM::Enrichment::typeVoid);

  bool skipped_element = false;

  const int nen = xfemele->NumNode();
  const int* nodeidptrs = xfemele->NodeIds();
  for (int inen = 0; inen<nen; ++inen)
  {
    const int node_gid = nodeidptrs[inen];

      const bool usefull_contribution = (fabs(ratios[inen]) > volumeRatioLimit);
      if ( usefull_contribution)
      {
        for (std::set<XFEM::PHYSICS::Field>::const_iterator field = fieldset.begin();field != fieldset.end();++field)
        {
          nodalDofSet[node_gid].insert(XFEM::FieldEnr(*field, voidenr));
        }
      }
      else
      {
        skipped_element = true;
        cout << "skipped interior void unknowns for element: "<< xfemele->Id() << ", for node: "<< node_gid << ", volumeratio limit: " << std::scientific << volumeRatioLimit << ", volumeratio: abs (" << std::scientific << fabs(ratios[inen]) << " )" << endl;
        const LINALG::Matrix<3,1> nodalpos(ih.xfemdis()->gNode(node_gid)->X());
        const bool in_fluid = (0 == ih.PositionWithinConditionNP(nodalpos));

        if (in_fluid)
        {
          for (std::set<XFEM::PHYSICS::Field>::const_iterator field = fieldset.begin();field != fieldset.end();++field)
          {
            nodalDofSet[node_gid].insert(XFEM::FieldEnr(*field, voidenr));
          }
        }
      }
  }
  return skipped_element;
}
#endif
//static void findOutMultitude(
//    const XFEM::InterfaceHandleXFSI&              ih,
//    set<int>                                      beles,
//    vector<set<int> >&                            distinct_interfaces_NGid,
//    vector<set<int> >&                            distinct_interfaces_BeleGid
//    )
//{
////  cout << "numbele: " << beles.size() << endl;
//  for (set<int>::const_iterator ibele = beles.begin(); ibele != beles.end(); ++ibele)
//  {
//    const DRT::Element* bele = ih.cutterdis()->gElement(*ibele);
//
//    const int* bnodeids = bele->NodeIds();
//    set<int> current_nodeid_set;
//    for (int inode = 0; inode < bele->NumNode(); ++inode)
//    {
//      current_nodeid_set.insert(bnodeids[inode]);
//    }
//
//    // if first element, create patch
//    if (distinct_interfaces_NGid.empty())
//    {
//      distinct_interfaces_NGid.push_back(current_nodeid_set);
//      distinct_interfaces_BeleGid.push_back(set<int>());
//      distinct_interfaces_BeleGid[0].insert(bele->Id());
//      continue;
//    }
//
//    for (size_t ipatch = 0; ipatch < distinct_interfaces_NGid.size(); ++ipatch)
//    {
//      const bool connected_elements = XFEM::ConnectedElements(distinct_interfaces_NGid[ipatch],current_nodeid_set);
//
//      if (connected_elements)
//      {
//        // add bele to existing set
//        distinct_interfaces_NGid[ipatch].insert(current_nodeid_set.begin(),current_nodeid_set.end());
//        distinct_interfaces_BeleGid[ipatch].insert(bele->Id());
//      }
//      else
//      {
//        // create new patch with bele as first entry
//        distinct_interfaces_NGid.push_back(current_nodeid_set);
//        distinct_interfaces_BeleGid.push_back(set<int>());
//        distinct_interfaces_BeleGid[ipatch+1].insert(bele->Id());
//      }
//    }
//  }
//
//  // sanity check
//  const size_t initnumpatch = distinct_interfaces_NGid.size();
//  if (initnumpatch != distinct_interfaces_BeleGid.size())
//  {
//    dserror("bug!");
//  }
//
//  //check patches with themselves
//  for (size_t ibase=0; ibase < initnumpatch; ++ibase)
//  {
//    if (ibase == 1)
//      break;  // we don't need to merge just one set
//
//    if (ibase == distinct_interfaces_NGid.size())
//      break;  //
//
//    bool restart = false;
//    const set<int>& baseset = distinct_interfaces_NGid[ibase];
//    for (size_t itest=1; itest < distinct_interfaces_NGid.size(); ++itest)
//    {
//      const set<int>& testset = distinct_interfaces_NGid[itest];
//      const bool connected_elements = XFEM::ConnectedElements(baseset,testset);
//      if (connected_elements)
//      {
//        // add things to base sets
//        distinct_interfaces_NGid[ibase].insert(testset.begin(),testset.end());
//        distinct_interfaces_BeleGid[ibase].insert(distinct_interfaces_BeleGid[itest].begin(),distinct_interfaces_BeleGid[itest].end());
//
//        // remove test sets
//        distinct_interfaces_NGid.erase(distinct_interfaces_NGid.begin() + itest);
//        distinct_interfaces_BeleGid.erase(distinct_interfaces_BeleGid.begin() + itest);
//
//        restart = true;
//      }
//      if (restart)
//        break;
//    }
//  }
//
//
////  cout << "distint_interfaces.size() = " << distinct_interfaces_NGid.size() << endl;
//
//}






/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
static bool ApplyElementVoidEnrichments(
    const DRT::Element*                           xfemele,
    const map<XFEM::PHYSICS::Field, DRT::Element::DiscretizationType>&  element_ansatz,
    const XFEM::InterfaceHandleXFSI&              ih,
    const int&                                    label,
    const XFEM::Enrichment::EnrType               enrtype,
    const double                                  boundaryRatioLimit,
    std::set<XFEM::FieldEnr>&                     enrfieldset)
{
  // check, how much area for integration we have (from BoundaryIntcells)
  const double boundarysize = XFEM::BoundaryCoverageRatio(*xfemele,ih.GetBoundaryIntCells(xfemele->Id()),ih);
  const bool almost_zero_surface = (fabs(boundarysize) < boundaryRatioLimit);

  bool skipped_element = false;

  const XFEM::Enrichment enr(enrtype, label);
  //  const XFEM::Enrichment stdenr(0, XFEM::Enrichment::typeStandard);
  if ( not almost_zero_surface)
  {
    const bool VoidFSIenrichment_in_set = EnrichmentInDofSet(XFEM::Enrichment::typeVoidFSI, enrfieldset);
    if (not VoidFSIenrichment_in_set)
    {
      map<XFEM::PHYSICS::Field, DRT::Element::DiscretizationType>::const_iterator fielditer;
      for (fielditer = element_ansatz.begin();fielditer != element_ansatz.end();++fielditer)
      {
        enrfieldset.insert(XFEM::FieldEnr(fielditer->first, enr));
      }
    }
  }
  else
  {
    skipped_element = true;
//    cout << "skipped stress unknowns for element: "<< xfemele->Id() << ", boundary size: " << boundarysize << endl;
  }
  return skipped_element;
}



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
static void processVoidFSIEnrichmentForElementNodes(
    const DRT::Element*                           xfemele,
    const map<XFEM::PHYSICS::Field, DRT::Element::DiscretizationType>&  element_ansatz,
    const XFEM::InterfaceHandleXFSI&              ih,
    const set<int>&                               labels,
    const std::set<XFEM::PHYSICS::Field>&         fieldset,
    const double                                  volumeRatioLimit,
    std::map<int, std::set<XFEM::FieldEnr> >&     nodalDofSet,
    bool&                                         skipped_node_enr
    )
{

  const int xele_gid = xfemele->Id();

  if (ih.ElementIntersected(xele_gid))
  {
    skipped_node_enr = ApplyNodalEnrichmentsVoidFSI(xfemele, ih, XFEM::Enrichment::typeVoidFSI, fieldset, volumeRatioLimit, nodalDofSet);
      //ApplyNodalEnrichmentsNodeWise(xfemele, ih, label, enrtype, fieldset, 2.0e-2, nodalDofSet);
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
static void processVoidEnrichmentForElementNodes(
    const DRT::Element*                           xfemele,
    const map<XFEM::PHYSICS::Field, DRT::Element::DiscretizationType>&  element_ansatz,
    const XFEM::InterfaceHandleXFSI&              ih,
    const set<int>&                               labels,
    const std::set<XFEM::PHYSICS::Field>&         fieldset,
    const double                                  volumeRatioLimit,
    std::map<int, std::set<XFEM::FieldEnr> >&     nodalDofSet,
    bool&                                         skipped_node_enr
    )
{

  const int xele_gid = xfemele->Id();

  if (ih.ElementIntersected(xele_gid))
  {
    for(std::set<int>::const_iterator labeliter = labels.begin(); labeliter!=labels.end(); ++labeliter)
    {
      const int label = *labeliter;

      if (ih.ElementIntersected(xele_gid))
      {
        if (ih.ElementHasLabel(xele_gid, label))
        {
          skipped_node_enr = ApplyNodalEnrichmentsVoid(xfemele, ih, label, XFEM::Enrichment::typeVoid, fieldset, volumeRatioLimit, nodalDofSet);
          //ApplyNodalEnrichmentsNodeWise(xfemele, ih, label, enrtype, fieldset, 2.0e-2, nodalDofSet);
        }
      }
    }
  }
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void XFEM::processVoidEnrichmentForElement(
    const DRT::Element*                           xfemele,
    const map<XFEM::PHYSICS::Field, DRT::Element::DiscretizationType>&  element_ansatz,
    const XFEM::InterfaceHandleXFSI&              ih,
    const set<int>&                               labels,
    const double                                  boundaryRatioLimit,
    std::set<XFEM::FieldEnr>&                     elementDofs,
    bool&                                         skipped_elem_enr
    )
{

  const int xele_gid = xfemele->Id();

  if (ih.ElementIntersected(xele_gid))
  {
    skipped_elem_enr = ApplyElementEnrichmentsVoidFSI(xfemele, element_ansatz, ih, boundaryRatioLimit, elementDofs);
  }

  for(std::set<int>::const_iterator labeliter = labels.begin(); labeliter!=labels.end(); ++labeliter)
  {
    const int label = *labeliter;

    if (ih.ElementIntersected(xele_gid))
    {
      if (ih.ElementHasLabel(xele_gid, label))
      {
        skipped_elem_enr = ApplyElementVoidEnrichments(xfemele, element_ansatz, ih, label, XFEM::Enrichment::typeVoid, boundaryRatioLimit, elementDofs);
      }
    }
  }
}

///*----------------------------------------------------------------------*
// *----------------------------------------------------------------------*/
//static void applyMultitudeNodalVoidEnrichments(
//    const XFEM::InterfaceHandleXFSI&              ih,
//    const std::set<XFEM::PHYSICS::Field>&         fieldset,
//    std::map<int, std::set<XFEM::FieldEnr> >&     nodalDofSet)
//{
//
//  for (int i=0; i<ih.xfemdis()->NumMyColNodes(); ++i)
//  {
//    DRT::Node* node = ih.xfemdis()->lColNode(i);
//
//    DRT::Element** eles = node->Elements();
////    cout << "vorher " << nodalDofSet[node->Id()].size() << endl;
//    vector<set<int> > distinct_interfaces_NGid;
//    vector<set<int> > distinct_interfaces_BeleGid;
//    for (int iele = 0; iele < node->NumElement(); ++iele)
//    {
//      DRT::Element* fluidele = eles[iele];
//
//      set<int> beles;
//      if (ih.ElementIntersected(fluidele->Id()))
//        beles = ih.GetIntersectingBoundaryElementsGID(fluidele->Id());
//
//      findOutMultitude(ih, beles, distinct_interfaces_NGid, distinct_interfaces_BeleGid);
//
//      // loop here
//      if (distinct_interfaces_NGid.size() > 1)
//      {
////        cout << "found multitude" << endl;
//        const LINALG::Matrix<3,1> nodalpos(node->X());
//        const int label = ih.PositionWithinConditionNP(nodalpos);
//
//        const XFEM::Enrichment enr_multi(label, XFEM::Enrichment::typeVoid);
//
//
//        for (std::set<XFEM::PHYSICS::Field>::const_iterator field = fieldset.begin();field != fieldset.end();++field)
//        {
//          nodalDofSet[node->Id()].insert(XFEM::FieldEnr(*field, enr_multi));
//        }
//
//      }
//
//    }
////    cout << "nachher" << nodalDofSet[node->Id()].size() << endl;
//  };
//}




/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
static void processStandardEnrichmentNodalBasedApproach(
    const XFEM::InterfaceHandle&                  ih,
    const std::set<XFEM::PHYSICS::Field>&         fieldset,
    std::map<int, std::set<XFEM::FieldEnr> >&     nodalDofSet)
{
  const XFEM::Enrichment enr_std(XFEM::Enrichment::typeStandard, 0);
  for (int i=0; i<ih.xfemdis()->NumMyColNodes(); ++i)
  {
    const DRT::Node* node = ih.xfemdis()->lColNode(i);
    const LINALG::Matrix<3,1> nodalpos(node->X());

    const bool voidfsienrichment_in_set = EnrichmentInNodalDofSet(node->Id(), XFEM::Enrichment::typeVoidFSI, nodalDofSet);
    const bool voidenrichment_in_set = EnrichmentInNodalDofSet(node->Id(), XFEM::Enrichment::typeVoid, nodalDofSet);

    if (not voidenrichment_in_set and not voidfsienrichment_in_set)
    {
      const bool in_fluid = (0 == ih.PositionWithinConditionNP(nodalpos));

      if (in_fluid)
      {
        for (std::set<XFEM::PHYSICS::Field>::const_iterator field = fieldset.begin();field != fieldset.end();++field)
        {
          nodalDofSet[node->Id()].insert(XFEM::FieldEnr(*field, enr_std));
        }
      }
    }
  };
}


#if 0
/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
static void applyStandardEnrichment(
    const XFEM::InterfaceHandle&                  ih,
    const std::set<XFEM::PHYSICS::Field>&         fieldset,
    std::map<int, std::set<XFEM::FieldEnr> >&     nodalDofSet,
    std::map<int, std::set<XFEM::FieldEnr> >&     elementalDofs)
{
  const int standard_label = 0;
  const XFEM::Enrichment enr_std(standard_label, XFEM::Enrichment::typeStandard);
  for (int i=0; i<ih.xfemdis()->NumMyColElements(); ++i)
  {
    const DRT::Element* xfemele = ih.xfemdis()->lColElement(i);
    if ( not ih.ElementIntersected(xfemele->Id()))
    {
      const int* nodeidptrs = xfemele->NodeIds();
      const LINALG::Matrix<3,1> nodalpos(xfemele->Nodes()[0]->X());

      const int label = ih.PositionWithinConditionNP(nodalpos);
      const bool in_fluid = (0 == label);

      if (in_fluid)
      {
        for (int inen = 0; inen<xfemele->NumNode(); ++inen)
        {
          const int node_gid = nodeidptrs[inen];
          bool voidenrichment_in_set = false;
          //check for void enrichement in a given set, if such set already exists for this node_gid
          std::map<int, std::set<XFEM::FieldEnr> >::const_iterator setiter = nodalDofSet.find(node_gid);
          if (setiter != nodalDofSet.end())
          {

            std::set<XFEM::FieldEnr> fieldenrset = setiter->second;
            for (std::set<XFEM::FieldEnr>::const_iterator fieldenr = fieldenrset.begin(); fieldenr != fieldenrset.end(); ++fieldenr)
            {
              if (fieldenr->getEnrichment().Type() == XFEM::Enrichment::typeVoid)
              {
                voidenrichment_in_set = true;
                break;
              }
            }
          }
          if (not voidenrichment_in_set)
          {
            for (std::set<XFEM::PHYSICS::Field>::const_iterator field = fieldset.begin();field != fieldset.end();++field)
            {
              nodalDofSet[node_gid].insert(XFEM::FieldEnr(*field, enr_std));
            }
          }
        };

        //                // add continuous stress unknowns
        //                const int element_gid = actele->Id();
        //                elementalDofs[element_gid].insert(XFEM::FieldEnr(PHYSICS::Tauxx, enr_std));
        //                elementalDofs[element_gid].insert(XFEM::FieldEnr(PHYSICS::Tauyy, enr_std));
        //                elementalDofs[element_gid].insert(XFEM::FieldEnr(PHYSICS::Tauzz, enr_std));
        //                elementalDofs[element_gid].insert(XFEM::FieldEnr(PHYSICS::Tauxy, enr_std));
        //                elementalDofs[element_gid].insert(XFEM::FieldEnr(PHYSICS::Tauxz, enr_std));
        //                elementalDofs[element_gid].insert(XFEM::FieldEnr(PHYSICS::Tauyz, enr_std));
      }
    }
  };
}
#endif

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void XFEM::createDofMapFSI(
    const XFEM::InterfaceHandleXFSI&                    ih,
    std::map<int, const std::set<XFEM::FieldEnr> >&     nodalDofSetFinal,
    std::map<int, const std::set<XFEM::FieldEnr> >&     elementalDofsFinal,
    const std::set<XFEM::PHYSICS::Field>&               fieldset,
    const XFEM::ElementAnsatz&                          elementAnsatz,
    const Teuchos::ParameterList&                       params
    )
{
  // temporary assembly
  std::map<int, std::set<XFEM::FieldEnr> >  nodalDofSet;
  std::map<int, std::set<XFEM::FieldEnr> >  elementalDofs;

  // get list of coupling label
  const std::set<int> labels = ih.GetAvailableBoundaryLabels();

  const double volumeRatioLimit = params.get<double>("volumeRatioLimit");
  const double boundaryRatioLimit = params.get<double>("boundaryRatioLimit");

  int skipped_node_enr_count = 0;
  int skipped_elem_enr_count = 0;

  // first, the connecting enrichments have to be fully applied
  for (int i=0; i<ih.xfemdis()->NumMyColElements(); ++i)
  {
    const DRT::Element* xfemele = ih.xfemdis()->lColElement(i);
    // add discontinuous stress unknowns
    // the number of each of these parameters will be determined later
    // by using a discretization type and appropriate shape functions
    map<XFEM::PHYSICS::Field, DRT::Element::DiscretizationType> element_ansatz;
    if (not params.get<bool>("DLM_condensation"))
    {
      element_ansatz = elementAnsatz.getElementAnsatz(xfemele->Shape());
    }

    bool skipped_node_enr = false;

    processVoidFSIEnrichmentForElementNodes(
        xfemele, element_ansatz, ih, labels, fieldset,
        volumeRatioLimit,
        nodalDofSet,
        skipped_node_enr);
    if (skipped_node_enr)
      skipped_node_enr_count++;

  };

#ifdef PARALLEL
  // first sync to get VoidFSI enrichments onto all processors, such that
  // Void enrichments can be applied correctly
  syncNodalDofs(ih, nodalDofSet);
#endif

  // now, all nodal void enrichments as well as the element erichments can be applied
  for (int i=0; i<ih.xfemdis()->NumMyColElements(); ++i)
  {
    const DRT::Element* xfemele = ih.xfemdis()->lColElement(i);
    // add discontinuous stress unknowns
    // the number of each of these parameters will be determined later
    // by using a discretization type and appropriate shape functions
    map<XFEM::PHYSICS::Field, DRT::Element::DiscretizationType> element_ansatz;
    if (not params.get<bool>("DLM_condensation"))
    {
      element_ansatz = elementAnsatz.getElementAnsatz(xfemele->Shape());
    }

    bool skipped_node_enr = false;
    bool skipped_elem_enr = false;

    XFEM::processVoidEnrichmentForElement(
        xfemele, element_ansatz, ih, labels,
        boundaryRatioLimit,
        elementalDofs[xfemele->Id()],
        skipped_elem_enr);
    if (skipped_elem_enr)
      skipped_elem_enr_count++;

    processVoidEnrichmentForElementNodes(
        xfemele, element_ansatz, ih, labels, fieldset,
        volumeRatioLimit,
        nodalDofSet,
        skipped_node_enr);
    if (skipped_node_enr)
      skipped_node_enr_count++;

  };


//  for(std::set<int>::const_iterator labeliter = labels.begin(); labeliter!=labels.end(); ++labeliter)
//  {
//    applyMultitudeNodalVoidEnrichments(ih, fieldset, nodalDofSet);
//  };


#ifdef PARALLEL
  syncNodalDofs(ih, nodalDofSet);
#endif

  cout << " skipped node unknowns for "<< skipped_node_enr_count << " elements (volumeratio limit:   " << std::scientific << volumeRatioLimit   << ")" << endl;
  cout << " skipped elem unknowns for "<< skipped_elem_enr_count << " elements (boundaryratio limit: " << std::scientific << boundaryRatioLimit << ")" << endl;

  processStandardEnrichmentNodalBasedApproach(ih, fieldset, nodalDofSet);

//#ifdef PARALLEL
//  syncNodalDofs(ih, nodalDofSet);
//#endif

  // create const sets from standard sets, so the sets cannot be accidently changed
  // could be removed later, if this is a performance bottleneck
  for ( std::map<int, std::set<XFEM::FieldEnr> >::const_iterator oneset = nodalDofSet.begin(); oneset != nodalDofSet.end(); ++oneset )
  {
    nodalDofSetFinal.insert( make_pair(oneset->first, oneset->second));
  };

  for ( std::map<int, std::set<XFEM::FieldEnr> >::const_iterator onevec = elementalDofs.begin(); onevec != elementalDofs.end(); ++onevec )
  {
    elementalDofsFinal.insert( make_pair(onevec->first, onevec->second));
  };
}






#endif  // #ifdef CCADISCRET
