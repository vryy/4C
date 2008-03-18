/*!----------------------------------------------------------------------
\file fluid_utils.cpp
\brief utility functions for fluid problems

<pre>
Maintainer: Axel Gerstenberger
            gerstenberger@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include <stdio.h>

#include "fluid_utils.H"


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void FLUID_UTILS::SetupFluidSplit(const DRT::Discretization& dis,
                                 int ndim,
                                 LINALG::MapExtractor& extractor)
{
  std::set<int> conddofset;
  std::set<int> otherdofset;

  int numrownodes = dis.NumMyRowNodes();
  for (int i=0; i<numrownodes; ++i)
  {
    DRT::Node* node = dis.lRowNode(i);

    std::vector<int> dof = dis.Dof(node);
    for (unsigned j=0; j<dof.size(); ++j)
    {
      // test for dof position
      if (j<static_cast<unsigned>(ndim))
      {
        otherdofset.insert(dof[j]);
      }
      else
      {
        conddofset.insert(dof[j]);
      }
    }
  }

  std::vector<int> conddofmapvec;
  conddofmapvec.reserve(conddofset.size());
  conddofmapvec.assign(conddofset.begin(), conddofset.end());
  conddofset.clear();
  Teuchos::RCP<Epetra_Map> conddofmap =
    Teuchos::rcp(new Epetra_Map(-1,
                                conddofmapvec.size(),
                                &conddofmapvec[0],
                                0,
                                dis.Comm()));
  conddofmapvec.clear();

  std::vector<int> otherdofmapvec;
  otherdofmapvec.reserve(otherdofset.size());
  otherdofmapvec.assign(otherdofset.begin(), otherdofset.end());
  otherdofset.clear();
  Teuchos::RCP<Epetra_Map> otherdofmap =
    Teuchos::rcp(new Epetra_Map(-1,
                                otherdofmapvec.size(),
                                &otherdofmapvec[0],
                                0,
                                dis.Comm()));
  otherdofmapvec.clear();

  extractor.Setup(*dis.DofRowMap(),conddofmap,otherdofmap);
}

// -------------------------------------------------------------------
// -------------------------------------------------------------------
void FLUID_UTILS::SetupXFluidSplit(
        const DRT::Discretization& dis,
        const RCP<XFEM::DofManager> dofman,
        LINALG::MapExtractor& extractor)
{
    // -------------------------------------------------------------------
    // get a vector layout from the discretization for a vector which only
    // contains the velocity dofs and for one vector which only contains
    // pressure degrees of freedom.
    //
    // The maps are designed assuming that every node has pressure and
    // velocity degrees of freedom --- this won't work for inf-sup stable
    // elements at the moment!
    // -------------------------------------------------------------------

    // Allocate integer vectors which will hold the dof number of the
    // velocity or pressure dofs
    vector<int> velmapdata;
    vector<int> premapdata;

    // collect global dofids for velocity and pressure in vectors
    for (int i=0; i<dis.NumMyRowNodes(); ++i) {
        const DRT::Node* node = dis.lRowNode(i);
        const std::set<XFEM::FieldEnr> enrvarset = dofman->getNodeDofSet(node->Id());
        const vector<int> dof = dis.Dof(node);
        dsassert(dof.size() == enrvarset.size(), "mismatch in length!");
        std::set<XFEM::FieldEnr>::const_iterator enrvar;
        unsigned int countdof = 0;
        for (enrvar = enrvarset.begin(); enrvar != enrvarset.end(); ++enrvar) {
            switch (enrvar->getField()) {
                case XFEM::PHYSICS::Velx:
                case XFEM::PHYSICS::Vely:
                case XFEM::PHYSICS::Velz:
                    velmapdata.push_back(dof[countdof]);
                    break;
                case XFEM::PHYSICS::Pres:
                    premapdata.push_back(dof[countdof]);
                    break;
                default:
                    break;
            }
            countdof++;
        }
    }

    // the rowmaps are generated according to the pattern provided by
    // the data vectors
    RCP<Epetra_Map> velrowmap = rcp(new Epetra_Map(-1,
            velmapdata.size(),&velmapdata[0],0,
            dis.Comm()));
    RCP<Epetra_Map> prerowmap = rcp(new Epetra_Map(-1,
            premapdata.size(),&premapdata[0],0,
            dis.Comm()));

    const Epetra_Map* map = dis.DofRowMap();
    extractor.Setup(*map, prerowmap, velrowmap);
}



#endif /* CCADISCRET       */
