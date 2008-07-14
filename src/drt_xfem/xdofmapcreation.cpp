/*!
\file xdofmapcreation.cpp

\brief defines unknowns based on the intersection pappern from the xfem intersection

this is related to the physics of the fluid problem and therefore not part of the standard xfem routines

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


/*----------------------------------------------------------------------*
 |  construct dofmap                                            ag 11/07|
 *----------------------------------------------------------------------*/
void XFEM::createDofMap(
    const XFEM::InterfaceHandle&                   ih,
    std::map<int, const set<XFEM::FieldEnr> >&     nodalDofSetFinal,
    std::map<int, const set<XFEM::FieldEnr> >&     elementalDofsFinal
)
{
  // temporary assembly
  std::map<int, set<XFEM::FieldEnr> >  nodalDofSet;
  std::map<int, set<XFEM::FieldEnr> >  elementalDofs;

  // get elements for each coupling label
  const std::map<int,set<int> >& elementsByLabel = *ih.elementsByLabel(); 

  // invert collection
  std::map<int,int> labelPerElementId;
  XFEM::InvertElementsByLabel(elementsByLabel, labelPerElementId);

  for(std::map<int,set<int> >::const_iterator conditer = elementsByLabel.begin(); conditer!=elementsByLabel.end(); ++conditer)
  {
    const int label = conditer->first;

    // for surface with label, loop my col elements and add void enrichments to each elements member nodes
    const XFEM::Enrichment voidenr(label, XFEM::Enrichment::typeVoid);
    for (int i=0; i<ih.xfemdis()->NumMyColElements(); ++i)
    {
      const DRT::Element* xfemele = ih.xfemdis()->lColElement(i);
      const int element_gid = xfemele->Id();

      if (ih.ElementIntersected(element_gid))
      {
        const XFEM::BoundaryIntCells& bcells = ih.elementalBoundaryIntCells()->find(element_gid)->second;
        bool has_label = false;
        for (BoundaryIntCells::const_iterator bcell = bcells.begin(); bcell != bcells.end(); ++bcell)
        {
          const int surface_ele_gid = bcell->GetSurfaceEleGid();
          const int label_for_current_bele = labelPerElementId.find(surface_ele_gid)->second;
          if (label == label_for_current_bele)
          {
            has_label = true;
            break;
          }
        }

        const double ratio = XFEM::AreaRatio(*xfemele,ih);
        const bool almost_empty_element = (fabs(1.0-ratio) < 1.0e-6);
        if (has_label)
        {

          if ( not almost_empty_element)  
          {
            const int nen = xfemele->NumNode();
            const int* nodeidptrs = xfemele->NodeIds();
            for (int inen = 0; inen<nen; ++inen)
            {
              const int node_gid = nodeidptrs[inen];
              nodalDofSet[node_gid].insert(XFEM::FieldEnr(PHYSICS::Velx, voidenr));
              nodalDofSet[node_gid].insert(XFEM::FieldEnr(PHYSICS::Vely, voidenr));
              nodalDofSet[node_gid].insert(XFEM::FieldEnr(PHYSICS::Velz, voidenr));
              nodalDofSet[node_gid].insert(XFEM::FieldEnr(PHYSICS::Pres, voidenr));
            };
          }

          // add discontinuous stress unknowns
          // the number of each of these parameters will be determined later
          // by using a discretization type and appropriate shape functions
          const map<XFEM::PHYSICS::Field, DRT::Element::DiscretizationType> element_ansatz(XFLUID::getElementAnsatz(xfemele->Shape()));

          map<XFEM::PHYSICS::Field, DRT::Element::DiscretizationType>::const_iterator fielditer;
          for (fielditer = element_ansatz.begin();fielditer != element_ansatz.end();++fielditer)
          {
            elementalDofs[element_gid].insert(XFEM::FieldEnr(fielditer->first, voidenr));
          }

        }
      }
    };
  };

  applyStandardEnrichment(ih, nodalDofSet, elementalDofs);

  // create const sets from standard sets, so the sets cannot be accidentily changed
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
    const XFEM::InterfaceHandle&             ih,
    std::map<int, set<XFEM::FieldEnr> >&     nodalDofSet,
    std::map<int, set<XFEM::FieldEnr> >&     elementalDofs
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
      const BlitzVec3 nodalpos(toBlitzArray(xfemele->Nodes()[0]->X()));

      bool in_fluid = false;
      const int label = PositionWithinCondition(nodalpos, ih);
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

#endif  // #ifdef CCADISCRET
