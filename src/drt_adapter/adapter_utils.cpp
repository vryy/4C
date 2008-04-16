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


#endif
