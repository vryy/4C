/*----------------------------------------------------------------------*/
/*!
\file adapter_utils.cpp

\brief

<pre>
Maintainer: Ulrich Kuettler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15238
</pre>
*/
/*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include "adapter_utils.H"
#include "adapter_coupling_mortar.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_utils.H"

// we need to know all element types for the ale mesh creation
#include "../drt_f2/fluid2.H"
#include "../drt_f3/fluid3.H"

#include "../drt_ale2/ale2.H"
#include "../drt_ale3/ale3.H"

#include <map>
#include <set>
#include <string>
#include <vector>
#include <algorithm>

#ifdef PARALLEL
#include <mpi.h>
#endif

extern "C" /* stuff which is c and is accessed from c++ */
{
#include "../headers/standardtypes.h"
}

/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::UTILS::SetupNDimExtractor(const DRT::Discretization& dis,
                                        std::string condname,
                                        LINALG::MapExtractor& extractor)
{
  std::set<int> conddofset;
  std::set<int> otherdofset;

  std::vector<DRT::Condition*> conds;
  dis.GetCondition(condname, conds);

  int numrownodes = dis.NumMyRowNodes();
  for (int i=0; i<numrownodes; ++i)
  {
    DRT::Node* node = dis.lRowNode(i);

    // test if node is covered by condition
    bool conditioned = false;
    for (unsigned j=0; j<conds.size(); ++j)
    {
      const vector<int>* n = conds[j]->Nodes();

      // DRT::Condition nodes are ordered by design! So we can perform a
      // binary search here.
      if (std::binary_search(n->begin(), n->end(), node->Id()))
      {
        conditioned = true;
        break;
      }
    }

    std::vector<int> dof = dis.Dof(node);
    for (unsigned j=0; j<dof.size(); ++j)
    {
      // test for condition coverage and dof position
      if (conditioned and j<static_cast<unsigned>(genprob.ndim))
      {
        conddofset.insert(dof[j]);
      }
      else
      {
        otherdofset.insert(dof[j]);
      }
    }
  }

  // if there is no such condition, do not waste any more time
  int conddofsetsize = static_cast<int>(conddofset.size());
  int size;
  dis.Comm().SumAll(&conddofsetsize,&size,1);
  if (size==0)
  {
    Teuchos::RCP<Epetra_Map> emptymap =
      Teuchos::rcp(new Epetra_Map(-1,
                                  0,
                                  NULL,
                                  0,
                                  dis.Comm()));
    extractor.Setup(*dis.DofRowMap(),emptymap,Teuchos::rcp(dis.DofRowMap(),false));
  }
  else
  {
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
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Map> ADAPTER::UTILS::ConditionNodeMap(const DRT::Discretization& dis,
                                                          std::string condname)
{
  std::set<int> condnodeset;

  std::vector<DRT::Condition*> conds;
  dis.GetCondition(condname, conds);

  int numrownodes = dis.NumMyRowNodes();
  for (int i=0; i<numrownodes; ++i)
  {
    DRT::Node* node = dis.lRowNode(i);

    // test if node is covered by condition
    for (unsigned j=0; j<conds.size(); ++j)
    {
      const vector<int>* n = conds[j]->Nodes();

      // DRT::Condition nodes are ordered by design! So we can perform a
      // binary search here.
      if (std::binary_search(n->begin(), n->end(), node->Id()))
      {
        condnodeset.insert(node->Id());
        break;
      }
    }
  }

  std::vector<int> condnodemapvec;
  condnodemapvec.reserve(condnodeset.size());
  condnodemapvec.assign(condnodeset.begin(), condnodeset.end());
  condnodeset.clear();
  Teuchos::RCP<Epetra_Map> condnodemap =
    Teuchos::rcp(new Epetra_Map(-1,
                                condnodemapvec.size(),
                                &condnodemapvec[0],
                                0,
                                dis.Comm()));
  condnodemapvec.clear();
  return condnodemap;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<std::set<int> > ADAPTER::UTILS::ConditionElementMap(const DRT::Discretization& dis,
                                                                 std::string condname)
{
  std::vector<DRT::Condition*> conds;
  dis.GetCondition(condname, conds);

  std::set<int> condelementset;
  int nummyelements = dis.NumMyColElements();
  for (int i=0; i<nummyelements; ++i)
  {
    DRT::Element* actele = dis.lColElement(i);
    int numnodes = actele->NumNode();
    DRT::Node** nodes = actele->Nodes();
    for (int n=0; n<numnodes; ++n)
    {
      DRT::Node* actnode = nodes[n];

      // test if node is covered by condition
      for (unsigned j=0; j<conds.size(); ++j)
      {
        const vector<int>* n = conds[j]->Nodes();

        // DRT::Condition nodes are ordered by design! So we can perform a
        // binary search here.
        if (std::binary_search(n->begin(), n->end(), actnode->Id()))
        {
          condelementset.insert(actele->Id());
          break;
        }
      }
    }
  }

  Teuchos::RCP<std::set<int> > condelementmap = Teuchos::rcp(new set<int>());
  swap(*condelementmap,condelementset);
  return condelementmap;
}


#endif
