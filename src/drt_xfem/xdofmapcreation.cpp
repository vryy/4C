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

#include <blitz/array.h>
#include "xdofmapcreation.H"
#include "xfem_condition.H"
#include "enrichment_utils.H"
#include "../drt_f3/xfluid3_interpolation.H"

void XFEM::ApplyNodalEnrichments(
    const DRT::Element*                           xfemele,
    const XFEM::InterfaceHandle&                  ih,
    const int&                                    label,
    std::map<int, std::set<XFEM::FieldEnr> >&     nodalDofSet
) 
{
  const double volumeratiolimit = 1.0e-2;
  
  const double volumeratio = XFEM::DomainCoverageRatio(*xfemele,ih);
  const bool almost_empty_element = (fabs(1.0-volumeratio) < volumeratiolimit);
  
  const XFEM::Enrichment voidenr(label, XFEM::Enrichment::typeVoid);
  
  if ( not almost_empty_element)  
  { // void enrichments for everybody !!!
    const int nen = xfemele->NumNode();
    const int* nodeidptrs = xfemele->NodeIds();
    for (int inen = 0; inen<nen; ++inen)
    {
      const int node_gid = nodeidptrs[inen];
      nodalDofSet[node_gid].insert(XFEM::FieldEnr(XFEM::PHYSICS::Velx, voidenr));
      nodalDofSet[node_gid].insert(XFEM::FieldEnr(XFEM::PHYSICS::Vely, voidenr));
      nodalDofSet[node_gid].insert(XFEM::FieldEnr(XFEM::PHYSICS::Velz, voidenr));
      nodalDofSet[node_gid].insert(XFEM::FieldEnr(XFEM::PHYSICS::Pres, voidenr));
    };
  }
  else
  { // void enrichments only in the fluid domain
    const int nen = xfemele->NumNode();
    const int* nodeidptrs = xfemele->NodeIds();
    for (int inen = 0; inen<nen; ++inen)
    {
      const int node_gid = nodeidptrs[inen];
      const LINALG::Matrix<3,1> nodalpos(ih.xfemdis()->gNode(node_gid)->X());
      const int label = ih.PositionWithinConditionNP(nodalpos);
      const bool in_fluid = (label == 0);
  
      if (in_fluid)
      {
        nodalDofSet[node_gid].insert(XFEM::FieldEnr(XFEM::PHYSICS::Velx, voidenr));
        nodalDofSet[node_gid].insert(XFEM::FieldEnr(XFEM::PHYSICS::Vely, voidenr));
        nodalDofSet[node_gid].insert(XFEM::FieldEnr(XFEM::PHYSICS::Velz, voidenr));
        nodalDofSet[node_gid].insert(XFEM::FieldEnr(XFEM::PHYSICS::Pres, voidenr));
      }
    };
    cout << "skipped interior void unknowns for element: "<< xfemele->Id() << ", volumeratio limit: " << std::scientific << volumeratiolimit << ", volumeratio: abs (" << std::scientific << (1.0 - volumeratio) << " )" << endl;
  }
}

void XFEM::ApplyElementEnrichments(
    const DRT::Element*                           xfemele,
    const XFEM::InterfaceHandle&                  ih,
    const int&                                    label,
    const bool                                    DLM_condensation,
    std::set<XFEM::FieldEnr>&                     enrfieldset
    )
{
  // check, how much area for integration we have (from BoundaryIntcells)
  const double boundarysize = XFEM::BoundaryCoverageRatio(*xfemele,ih);
  const bool almost_zero_surface = (fabs(boundarysize) < 1.0e-4);
  const XFEM::Enrichment voidenr(label, XFEM::Enrichment::typeVoid);
//  const XFEM::Enrichment stdenr(0, XFEM::Enrichment::typeStandard);
  if ( not almost_zero_surface) 
  {
    // add discontinuous stress unknowns
    // the number of each of these parameters will be determined later
    // by using a discretization type and appropriate shape functions
    map<XFEM::PHYSICS::Field, DRT::Element::DiscretizationType> element_ansatz; 
    if (not DLM_condensation)
    {
      element_ansatz = XFLUID::getElementAnsatz(xfemele->Shape());
    }

    map<XFEM::PHYSICS::Field, DRT::Element::DiscretizationType>::const_iterator fielditer;
    for (fielditer = element_ansatz.begin();fielditer != element_ansatz.end();++fielditer)
    {
      enrfieldset.insert(XFEM::FieldEnr(fielditer->first, voidenr));
//      enrfieldset.insert(XFEM::FieldEnr(fielditer->first, stdenr));
    }
  }
  else
  {
    cout << "skipped stress unknowns for element: "<< xfemele->Id() << ", boundary size: " << boundarysize << endl;
  }
}


void XFEM::ApplyVoidEnrichmentForElement(
    const DRT::Element*                           xfemele,
    const XFEM::InterfaceHandle&                  ih,
    const std::map<int,int>&                      labelPerElementId,
    const int&                                    label,
    const bool                                    DLM_condensation,
    std::map<int, std::set<XFEM::FieldEnr> >&     nodalDofSet,
    std::map<int, std::set<XFEM::FieldEnr> >&     elementalDofs
    )
{
  const int element_gid = xfemele->Id();
  
  if (ih.ElementIntersected(element_gid))
  {
    if (ih.ElementHasLabel(element_gid, label))
    {
      ApplyNodalEnrichments(xfemele, ih, label, nodalDofSet); 
  
      ApplyElementEnrichments(xfemele, ih, label, DLM_condensation, elementalDofs[element_gid]);
    }
  }
}

void XFEM::createDofMap(
    const XFEM::InterfaceHandle&                    ih,
    std::map<int, const std::set<XFEM::FieldEnr> >&     nodalDofSetFinal,
    std::map<int, const std::set<XFEM::FieldEnr> >&     elementalDofsFinal,
    const bool DLM_condensation
)
{
  // temporary assembly
  std::map<int, std::set<XFEM::FieldEnr> >  nodalDofSet;
  std::map<int, std::set<XFEM::FieldEnr> >  elementalDofs;

  // get elements for each coupling label
  const std::map<int,std::set<int> >& elementsByLabel = *ih.elementsByLabel(); 

  // invert collection
  std::map<int,int> labelPerElementId;
  XFEM::InvertElementsByLabel(elementsByLabel, labelPerElementId);

  // loop condition labels
  for(std::map<int,std::set<int> >::const_iterator conditer = elementsByLabel.begin(); conditer!=elementsByLabel.end(); ++conditer)
  {
    const int label = conditer->first;

    // for surface with label, loop my col elements and add enrichments to each elements member nodes
    for (int i=0; i<ih.xfemdis()->NumMyColElements(); ++i)
    {
      const DRT::Element* xfemele = ih.xfemdis()->lColElement(i);

      ApplyVoidEnrichmentForElement(
          xfemele, ih, labelPerElementId, label, DLM_condensation,
          nodalDofSet, elementalDofs);
    };
  };

  applyStandardEnrichmentNodalBasedApproach(ih, nodalDofSet, elementalDofs);

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


void XFEM::applyStandardEnrichment(
    const XFEM::InterfaceHandle&              ih,
    std::map<int, std::set<XFEM::FieldEnr> >&     nodalDofSet,
    std::map<int, std::set<XFEM::FieldEnr> >&     elementalDofs
)
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

      bool in_fluid = false;
      const int label = ih.PositionWithinConditionNP(nodalpos);
      if (label == 0)
      {
        in_fluid = true;
      }
      else
      {
        in_fluid = false;
      }

      if (in_fluid)
      {
        for (int inen = 0; inen<xfemele->NumNode(); ++inen)
        {
          const int node_gid = nodeidptrs[inen];
          bool voidenrichment_in_set = false;
          //check for void enrichement in a given set, if such set already exists for this node_gid
          std::map<int, std::set<FieldEnr> >::const_iterator setiter = nodalDofSet.find(node_gid);
          if (setiter != nodalDofSet.end())
          {

            std::set<FieldEnr> fieldenrset = setiter->second;
            for (std::set<FieldEnr>::const_iterator fieldenr = fieldenrset.begin(); fieldenr != fieldenrset.end(); ++fieldenr)
            {
              if (fieldenr->getEnrichment().Type() == Enrichment::typeVoid)
              {
                voidenrichment_in_set = true;
                break;
              }
            }
          }
          if (not voidenrichment_in_set)
          {
            nodalDofSet[node_gid].insert(XFEM::FieldEnr(PHYSICS::Velx, enr_std));
            nodalDofSet[node_gid].insert(XFEM::FieldEnr(PHYSICS::Vely, enr_std));
            nodalDofSet[node_gid].insert(XFEM::FieldEnr(PHYSICS::Velz, enr_std));
            nodalDofSet[node_gid].insert(XFEM::FieldEnr(PHYSICS::Pres, enr_std));
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

void XFEM::applyStandardEnrichmentNodalBasedApproach(
    const XFEM::InterfaceHandle&              ih,
    std::map<int, std::set<XFEM::FieldEnr> >&     nodalDofSet,
    std::map<int, std::set<XFEM::FieldEnr> >&     elementalDofs
)
{
  const int standard_label = 0;
  const XFEM::Enrichment enr_std(standard_label, XFEM::Enrichment::typeStandard);
  for (int i=0; i<ih.xfemdis()->NumMyColNodes(); ++i)
  {
    const DRT::Node* node = ih.xfemdis()->lColNode(i);
    const LINALG::Matrix<3,1> nodalpos(node->X());

    const int node_gid = node->Id();
    bool voidenrichment_in_set = false;
    //check for void enrichment in a given set, if such set already exists for this node_gid
    std::map<int, std::set<FieldEnr> >::const_iterator setiter = nodalDofSet.find(node_gid);
    if (setiter != nodalDofSet.end())
    {
      std::set<FieldEnr> fieldenrset = setiter->second;
      for (std::set<FieldEnr>::const_iterator fieldenr = fieldenrset.begin(); fieldenr != fieldenrset.end(); ++fieldenr)
      {
        if (fieldenr->getEnrichment().Type() == Enrichment::typeVoid)
        {
          voidenrichment_in_set = true;
          break;
        }
      }
    }
    if (not voidenrichment_in_set)
    {
      bool in_fluid = false;
      const int label = ih.PositionWithinConditionNP(nodalpos);
      if (label == 0)
      {
        in_fluid = true;
      }
      else
      {
        in_fluid = false;
      }
      if (in_fluid)
      {
        nodalDofSet[node_gid].insert(XFEM::FieldEnr(PHYSICS::Velx, enr_std));
        nodalDofSet[node_gid].insert(XFEM::FieldEnr(PHYSICS::Vely, enr_std));
        nodalDofSet[node_gid].insert(XFEM::FieldEnr(PHYSICS::Velz, enr_std));
        nodalDofSet[node_gid].insert(XFEM::FieldEnr(PHYSICS::Pres, enr_std));
      }
    }
    
  };
}




#endif  // #ifdef CCADISCRET
