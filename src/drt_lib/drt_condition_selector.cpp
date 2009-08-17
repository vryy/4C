
#ifdef CCADISCRET

#include "drt_condition_selector.H"
#include "linalg_utils.H"


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::UTILS::MultiConditionSelector::SetupExtractor(const DRT::Discretization& dis,
                                                        Teuchos::RCP<Epetra_Map> fullmap,
                                                        LINALG::MultiMapExtractor& extractor)
{
  SetupCondDofSets(dis);

  // find all non-conditioned dofs. It is much more clean to do that
  // afterwards.

  std::set<int> otherdofset(fullmap->MyGlobalElements(),
                            fullmap->MyGlobalElements() + fullmap->NumMyElements());

  for (unsigned j=0; j<conddofset_.size(); ++j)
  {
    std::set<int>& conddofset = conddofset_[j];

    std::set<int>::const_iterator conditer;
    for (conditer = conddofset.begin(); conditer != conddofset.end(); ++conditer)
    {
      otherdofset.erase(*conditer);
    }
  }

  // setup all maps. The "other" map goes first so it becomes the zeroth map
  // of the MultiMapExtractor.

  std::vector<RCP<const Epetra_Map> > maps;
  maps.reserve(conddofset_.size()+1);

  maps.push_back(LINALG::CreateMap(otherdofset, dis.Comm()));
  for (unsigned j=0; j<conddofset_.size(); ++j)
  {
    std::set<int>& conddofset = conddofset_[j];
    maps.push_back(LINALG::CreateMap(conddofset, dis.Comm()));
  }

  // MultiMapExtractor setup

  extractor.Setup(*fullmap,maps);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::UTILS::MultiConditionSelector::SetupCondDofSets(const DRT::Discretization& dis)
{
  conddofset_.resize(selectors_.size());

  int numrownodes = dis.NumMyRowNodes();
  for (int i=0; i<numrownodes; ++i)
  {
    DRT::Node* node = dis.lRowNode(i);

    for (unsigned j=0; j<selectors_.size(); ++j)
    {
      ConditionSelector& conds = *selectors_[j];

      // put all conditioned dofs into conddofset
      if (conds.ContainsNode(node->Id()))
      {
        std::vector<int> dof = dis.Dof(node);
        for (unsigned k=0; k<dof.size(); ++k)
        {
          // test for dof position
          if (conds.ContainsDof(dof[k],k))
          {
            conddofset_[j].insert(dof[k]);
          }
        }

        // Node has been found. Done with this one.
        break;
      }
    }
  }
}

#endif
