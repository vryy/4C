/*!---------------------------------------------------------------------*
\file combust_utils.cpp

\brief collection of functions in namespace COMBUST

<pre>
Maintainer: Florian Henke
            henke@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15265
</pre>
*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include <iostream>

#include "combust_utils.H"

//#include "xfem_condition.H"
#include "../drt_xfem/enrichment_utils.H"
#include "../drt_xfem/field_enriched.H"
#include "../drt_f3/xfluid3_interpolation.H"



/*------------------------------------------------------------------------------------------------*
 | print COMBUST module logo on screen                                                henke 06/08 |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::printCombustLogo()
{
    std::cout<<"     ___            ___    "<<std::endl;
    std::cout<<"    /   \\          /   \\ "<<std::endl;
    std::cout<<"    \\_   \\        /  __/ "<<std::endl;
    std::cout<<"     _\\   \\      /  /__  "<<" Das ist               "<<std::endl;
    std::cout<<"     \\___  \\____/   __/  "<<" das Verbrennungsmodul "<<std::endl;
    std::cout<<"         \\_       _/      "<<" in BACI               "<<std::endl;
    std::cout<<"           | @ @  \\_      "<<"                       "<<std::endl;
    std::cout<<"           |               "<<" Der Elch wird bald    "<<std::endl;
    std::cout<<"         _/     /\\        "<<" ein feuerspeiender    "<<std::endl;
    std::cout<<"        /o)  (o/\\ \\_     "<<" Drache sein!          "<<std::endl;
    std::cout<<"        \\_____/ /         "<<std::endl;
    std::cout<<"          \\____/          "<<std::endl;
    std::cout<<"                           "<<std::endl;
}


/*------------------------------------------------------------------------------------------------*
 | functions below copied from drt_xfem/xdofmapcreation.cpp                           henke 10/08 |
 |
 | These functions are actually not of general XFEM kind, but combustion (former xfluid) specific!
 | This is why they are "parked" here, in this combust_utils-file.
 |
 | - functionality for elemental dofs was removed
 | - option to condsense the DLMs was removed
 *------------------------------------------------------------------------------------------------*/

/*------------------------------------------------------------------------------------------------*
 | This function is only called in the dofmanager constructor (unused? henke 08/09)   henke 10/08 |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::createDofMap(
    const XFEM::InterfaceHandle&                    ih,
    std::map<int, const std::set<XFEM::FieldEnr> >&     nodalDofSetFinal,
    std::map<int, const std::set<XFEM::FieldEnr> >&     elementalDofsFinal
)
{
  dserror ("unused function COMBUST::createDofMap()");
  return;
//  // temporary assembly
//  std::map<int, std::set<XFEM::FieldEnr> >  nodalDofSet;
//  std::map<int, std::set<XFEM::FieldEnr> >  elementalDofs;
//
//  // get structure (I dont have that!) elements for each coupling label
//  // xfemlabel_1 - element1
//  //             - element2
//  // xfemlabel_2 - element6
//  //             - element9
//  const std::map<int,std::set<int> >& elementsByLabel = *ih.elementsByLabel();
//
//  // invert collection
//  // element1 - xfemlabel_1
//  // element2 - xfemlabel_1
//  // element6 - xfemlabel_2
//  // element9 - xfemlabel_2
//  std::map<int,int> labelPerElementId;
//  XFEM::InvertElementsByLabel(elementsByLabel, labelPerElementId);
//
//  // loop condition labels = loop xfsi bodies = loop level set functions
//  for(std::map<int,std::set<int> >::const_iterator conditer = elementsByLabel.begin(); conditer!=elementsByLabel.end(); ++conditer)
//  {
//    const int label = conditer->first;
//
//    // for surface with label, loop my col elements and add enrichments to each elements member nodes
//    // remark: this implies that the "enrichment band" has width 1 element layer on both sides of interface
//    for (int i=0; i<ih.xfemdis()->NumMyColElements(); ++i)
//    {
//      const DRT::Element* xfemele = ih.xfemdis()->lColElement(i);
//      const int element_gid = xfemele->Id();
//
//      if (ih.ElementIntersected(element_gid)) // Wurde das Element vom Interface geschnitten?
//      {
//        if (ih.ElementHasLabel(element_gid, label)) // Wurde das Element von diesem Interface (label) geschnitten?
//        {
//          ApplyNodalEnrichments(xfemele, ih, label, nodalDofSet);
//        }
//      }
//
//    };
//  };
//
//  applyStandardEnrichmentNodalBasedApproach(ih, nodalDofSet, elementalDofs);
//
//  // create const sets from standard sets, so the sets cannot be accidently changed
//  // could be removed later, if this is a performance bottleneck
//  for ( std::map<int, std::set<XFEM::FieldEnr> >::const_iterator oneset = nodalDofSet.begin(); oneset != nodalDofSet.end(); ++oneset )
//  {
//    nodalDofSetFinal.insert( make_pair(oneset->first, oneset->second));
//  };
//
//  for ( std::map<int, std::set<XFEM::FieldEnr> >::const_iterator onevec = elementalDofs.begin(); onevec != elementalDofs.end(); ++onevec )
//  {
//    elementalDofsFinal.insert( make_pair(onevec->first, onevec->second));
//  };
}

/*------------------------------------------------------------------------------------------------*
 | Diese Funktion beeinhaltet eine Fallunterscheidung, wenn nur ein kleines Eck eines Fluidelementes
 | abgeschnitten wird (siehe Zeichnung von Axel 28.10.08)                             henke 10/08 |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::ApplyNodalEnrichments(
    const DRT::Element*                           xfemele,
    const XFEM::InterfaceHandle&                  ih,
    const int&                                    label,
    std::map<int, std::set<XFEM::FieldEnr> >&     nodalDofSet
)
{
  dserror ("unused function COMBUST::ApplyNodalEnrichments()");
  return;
//  const double volumeratiolimit = 1.0e-4; // Hängt von Auflösung ab -> anpassen?
//
//  const double volumeratio = XFEM::DomainCoverageRatio(*xfemele,ih);
//  // "almost empty" means that the volume ratio is close to 1, meaning that almost the entire fluid
//  // element is covered by the structure element. "Empty" refers to the amount of fluid in a fluid
//  // element domain.
//  const bool almost_empty_element = (fabs(1.0-volumeratio) < volumeratiolimit);
//
//  const XFEM::Enrichment voidenr(label, XFEM::Enrichment::typeVoid);
//
//  if ( not almost_empty_element) // all nodes of that fluid element are enriched
//  { // void enrichments for everybody !!!
//    const int nen = xfemele->NumNode();
//    const int* nodeidptrs = xfemele->NodeIds();
//    for (int inen = 0; inen<nen; ++inen)
//    {
//      const int node_gid = nodeidptrs[inen];
//      nodalDofSet[node_gid].insert(XFEM::FieldEnr(XFEM::PHYSICS::Velx, voidenr));
//      nodalDofSet[node_gid].insert(XFEM::FieldEnr(XFEM::PHYSICS::Vely, voidenr));
//      nodalDofSet[node_gid].insert(XFEM::FieldEnr(XFEM::PHYSICS::Velz, voidenr));
//      nodalDofSet[node_gid].insert(XFEM::FieldEnr(XFEM::PHYSICS::Pres, voidenr));
//    };
//  }
//  else
//  { // void enrichments only in the fluid domain
//    // other nodes connected by intersected element edges will be enriched during the loop over
//    // adjacent elements          henke 10/08
//    const int nen = xfemele->NumNode();  // couldn't part of this stuff around here be moved outside the if/else almost_empty_element?
//    const int* nodeidptrs = xfemele->NodeIds();
//    for (int inen = 0; inen<nen; ++inen)
//    {
//      const int node_gid = nodeidptrs[inen];
//      const BlitzVec3 nodalpos(toBlitzArray(ih.xfemdis()->gNode(node_gid)->X()));
//      const int label = ih.PositionWithinConditionNP(nodalpos);
//      const bool in_fluid = (label == 0); // Abfrage, ob der Knoten im Fluid, oder in der Struktur liegt
//
//      if (in_fluid) // Knoten im Fluid soll ein Enrichment erhalten
//      {
//        nodalDofSet[node_gid].insert(XFEM::FieldEnr(XFEM::PHYSICS::Velx, voidenr));
//        nodalDofSet[node_gid].insert(XFEM::FieldEnr(XFEM::PHYSICS::Vely, voidenr));
//        nodalDofSet[node_gid].insert(XFEM::FieldEnr(XFEM::PHYSICS::Velz, voidenr));
//        nodalDofSet[node_gid].insert(XFEM::FieldEnr(XFEM::PHYSICS::Pres, voidenr));
//      }
//    };
//    cout << "skipped interior void unknowns for element: "<< xfemele->Id() << ", volumeratio limit: " << std::scientific << volumeratiolimit << ", volumeratio: abs ( 1.0 - " << std::scientific << volumeratio << " )" << endl;
//  }
}

/*------------------------------------------------------------------------------------------------*
 | This function is currently not used in Axels XFluid.                          henke 10/08 |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::applyStandardEnrichment(
    const XFEM::InterfaceHandle&              ih,
    std::map<int, std::set<XFEM::FieldEnr> >&     nodalDofSet,
    std::map<int, std::set<XFEM::FieldEnr> >&     elementalDofs
)
{
  // Der Inhalt ist gelöscht, wei die Funktion sowieso momentan nicht genutzt wird. Falls benötigt,
  // ist sie von Axel zu adaptieren.
  return;
}

/*------------------------------------------------------------------------------------------------*
 | This function is the currectly used function of the former                 henke 10/08 |
 | Es werden Fluid Knoten angereichert, die noch kein void Enrichment bekommen haben. Es wird also
 | sozusagen der Standard FEM Fall erzeugt.
 *------------------------------------------------------------------------------------------------*/
void COMBUST::applyStandardEnrichmentNodalBasedApproach(
    const XFEM::InterfaceHandle&              ih,
    std::map<int, std::set<XFEM::FieldEnr> >&     nodalDofSet,
    std::map<int, std::set<XFEM::FieldEnr> >&     elementalDofs
)
{
  return;
//  const int standard_label = 0;
//  const XFEM::Enrichment enr_std(standard_label, XFEM::Enrichment::typeStandard);
//  for (int i=0; i<ih.xfemdis()->NumMyColNodes(); ++i)
//  {
//    const DRT::Node* node = ih.xfemdis()->lColNode(i);
//    const BlitzVec3 nodalpos(toBlitzArray(node->X()));
//
//    const int node_gid = node->Id();
//    bool voidenrichment_in_set = false;
//    //check for void enrichement in a given set, if such set already exists for this node_gid
//    std::map<int, std::set<FieldEnr> >::const_iterator setiter = nodalDofSet.find(node_gid);
//    if (setiter != nodalDofSet.end())
//    {
//      std::set<FieldEnr> fieldenrset = setiter->second;
//      for (std::set<FieldEnr>::const_iterator fieldenr = fieldenrset.begin(); fieldenr != fieldenrset.end(); ++fieldenr)
//      {
//        if (fieldenr->getEnrichment().Type() == Enrichment::typeVoid)
//        {
//          voidenrichment_in_set = true;
//          break;
//        }
//      }
//    }
//    if (not voidenrichment_in_set)
//    {
//      bool in_fluid = false;
//      const int label = ih.PositionWithinConditionNP(nodalpos);
//      if (label == 0)
//      {
//        in_fluid = true;
//      }
//      else
//      {
//        in_fluid = false;
//      }
//      if (in_fluid)
//      {
//        nodalDofSet[node_gid].insert(XFEM::FieldEnr(PHYSICS::Velx, enr_std));
//        nodalDofSet[node_gid].insert(XFEM::FieldEnr(PHYSICS::Vely, enr_std));
//        nodalDofSet[node_gid].insert(XFEM::FieldEnr(PHYSICS::Velz, enr_std));
//        nodalDofSet[node_gid].insert(XFEM::FieldEnr(PHYSICS::Pres, enr_std));
//      }
//    }
//
//  };
}



#endif  // #ifdef CCADISCRET

