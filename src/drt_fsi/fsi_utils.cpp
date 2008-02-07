
#ifdef CCADISCRET

#include <map>
#include <set>
#include <string>
#include <vector>
#include <algorithm>

#include "fsi_utils.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/linalg_utils.H"

#include <Epetra_CrsMatrix.h>
#include <EpetraExt_RowMatrixOut.h>
#include <NOX.H>
#include <NOX_Epetra.H>
#include <Epetra_SerialDenseMatrix.h>

#ifdef PARALLEL
#include <mpi.h>
#endif

extern "C" /* stuff which is c and is accessed from c++ */
{
#include "../headers/standardtypes.h"
}

/*!----------------------------------------------------------------------
\brief file pointers

<pre>                                                         m.gee 8/00
This structure struct _FILES allfiles is defined in input_control_global.c
and the type is in standardtypes.h
It holds all file pointers and some variables needed for the FRSYSTEM
</pre>
*----------------------------------------------------------------------*/
extern struct _FILES  allfiles;

/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::UTILS::SetupInterfaceExtractor(const DRT::Discretization& dis,
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


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::UTILS::DumpJacobian(NOX::Epetra::Interface::Required& interface,
                              double alpha,
                              double beta,
                              Teuchos::RefCountPtr<Epetra_Vector> soln,
                              std::string filename)
{
  // that's really stupid again
  const Epetra_BlockMap& bmap = soln->Map();
  Epetra_Map map(bmap.NumGlobalElements(),
                 bmap.NumMyElements(),
                 bmap.MyGlobalElements(),
                 0,
                 bmap.Comm());

  RefCountPtr<Epetra_CrsMatrix> jacobian = rcp(new Epetra_CrsMatrix(Copy, map, map.NumGlobalElements()));

  int nummyelements = map.NumMyElements();
  int mypos = LINALG::FindMyPos(nummyelements, map.Comm());
  double eta = 0.0;

  Epetra_Vector fo(*soln);
  Epetra_Vector fp(*soln);
  Epetra_Vector Jc(*soln);

  // Compute the RHS at the initial solution
  interface.computeF(*soln, fo, NOX::Epetra::Interface::Required::FD_Res);

  Epetra_Vector x_perturb = *soln;

  for (int i=0; i<map.NumGlobalElements(); ++i)
  {
    if (map.Comm().MyPID()==0)
      cout << "calculate column " << i << "\n";

    int proc = 0;
    int idx = 0;
    if (i>=mypos and i<mypos+nummyelements)
    {
      eta = alpha*(*soln)[i-mypos] + beta;
      x_perturb[i-mypos] += eta;
      idx = map.GID(i-mypos);
      proc = map.Comm().MyPID();
    }

    // Find what proc eta is on
    int broadcastProc = 0;
    map.Comm().SumAll(&proc, &broadcastProc, 1);

    // Send the perturbation variable, eta, to all processors
    map.Comm().Broadcast(&eta, 1, broadcastProc);

    map.Comm().Broadcast(&idx, 1, broadcastProc);

    // Compute the perturbed RHS
    interface.computeF(x_perturb, fp, NOX::Epetra::Interface::Required::FD_Res);

    // Compute the column k of the Jacobian
    Jc.Update(1.0, fp, -1.0, fo, 0.0);
    Jc.Scale( 1.0/eta );

    // Insert nonzero column entries into the jacobian
    for (int j = 0; j < map.NumMyElements(); ++j)
    {
      int gid = map.GID(j);
      if (Jc[j] != 0.0)
      {
        int err = jacobian->SumIntoGlobalValues(gid,1,&Jc[j],&idx);
        if (err>0)
        {
          err = jacobian->InsertGlobalValues(gid,1,&Jc[j],&idx);
        }
        if (err != 0)
          dserror("Assembly failed");
      }
    }

    // Unperturb the solution vector
    x_perturb = *soln;
  }

  jacobian->FillComplete();

  EpetraExt::RowMatrixToMatlabFile(filename.c_str(),*jacobian);
}

#endif
