/*!----------------------------------------------------------------------
\file drt_utils_parallel.cpp
\brief A collection of helper methods for namespace DRT

<pre>
-------------------------------------------------------------------------
                 BACI finite element library subsystem
            Copyright (2008) Technical University of Munich

Under terms of contract T004.008.000 there is a non-exclusive license for use
of this work by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library is proprietary software. It must not be published, distributed,
copied or altered in any form or any media without written permission
of the copyright holder. It may be used under terms and conditions of the
above mentioned license by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library may solemnly used in conjunction with the BACI contact library
for purposes described in the above mentioned contract.

This library contains and makes use of software copyrighted by Sandia Corporation
and distributed under LGPL licence. Licensing does not apply to this or any
other third party software used here.

Questions? Contact Dr. Michael W. Gee (gee@lnm.mw.tum.de)
                   or
                   Prof. Dr. Wolfgang A. Wall (wall@lnm.mw.tum.de)

http://www.lnm.mw.tum.de

-------------------------------------------------------------------------
</pre>
<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/


#include "drt_utils_parallel.H"
#include "drt_node.H"
#include "../drt_nurbs_discret/drt_control_point.H"
#include "drt_dofset.H"
#include "drt_discret.H"
#include "../linalg/linalg_utils.H"
#include "../linalg/linalg_mapextractor.H"
#include "drt_dserror.H"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::UTILS::RedistributeWithNewNodalDistribution(
    DRT::Discretization&     dis,
    const Epetra_Map&        noderowmap,
    const Epetra_Map&        nodecolmap
    )
{
  // redistribute nodes to column (ghost) map
  dis.ExportColumnNodes(nodecolmap);

  Teuchos::RCP< Epetra_Map > elerowmap;
  Teuchos::RCP< Epetra_Map > elecolmap;

  // now we have all elements in a linear map roweles
  // build resonable maps for elements from the
  // already valid and final node maps
  dis.BuildElementRowColumn(noderowmap, nodecolmap, elerowmap, elecolmap);

  // we can now export elements to resonable row element distribution
  dis.ExportRowElements(*elerowmap);

  // export to the column map / create ghosting of elements
  dis.ExportColumnElements(*elecolmap);

}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::UTILS::PrintParallelDistribution(const DRT::Discretization& dis)
{
  const int numproc=dis.Comm().NumProc();

  if(numproc>1)
  {
    const int myrank=dis.Comm().MyPID();

    vector<int> my_n_nodes     (numproc,0);
    vector<int>    n_nodes     (numproc,0);
    vector<int> my_n_ghostnodes(numproc,0);
    vector<int>    n_ghostnodes(numproc,0);
    vector<int> my_n_elements  (numproc,0);
    vector<int>    n_elements  (numproc,0);
    vector<int> my_n_ghostele  (numproc,0);
    vector<int>    n_ghostele  (numproc,0);

    my_n_nodes     [myrank]=dis.NumMyRowNodes();
    my_n_ghostnodes[myrank]=dis.NumMyColNodes()-my_n_nodes[myrank];
    my_n_elements  [myrank]=dis.NumMyRowElements();
    my_n_ghostele  [myrank]=dis.NumMyColElements()-my_n_elements[myrank];

    dis.Comm().SumAll(&my_n_nodes     [0],&n_nodes     [0],numproc);
    dis.Comm().SumAll(&my_n_ghostnodes[0],&n_ghostnodes[0],numproc);
    dis.Comm().SumAll(&my_n_elements  [0],&n_elements  [0],numproc);
    dis.Comm().SumAll(&my_n_ghostele  [0],&n_ghostele  [0],numproc);

    if(myrank==0)
    {
      cout << endl;
      cout <<"   Discretization: " << dis.Name() << endl;
      printf("   +-----+---------------+--------------+-----------------+----------------+\n");
      printf("   | PID |  n_rownodes   | n_ghostnodes |  n_rowelements  |   n_ghostele   |\n");
      printf("   +-----+---------------+--------------+-----------------+----------------+\n");
      for(int npid=0;npid<numproc;++npid)
      {
        printf("   | %3d | %13d | %12d | %15d | %14d |\n",npid,n_nodes[npid],n_ghostnodes[npid],n_elements[npid],n_ghostele[npid]);
        printf("   +-----+---------------+--------------+-----------------+----------------+\n");
      }
      cout << endl;
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> DRT::UTILS::GetColVersionOfRowVector(
    const Teuchos::RCP<const DRT::Discretization> dis,
    const Teuchos::RCP<const Epetra_Vector> state)
{
  // note that this routine has the same functionality as SetState,
  // although here we do not store the new vector anywhere
  // maybe this routine can be used in SetState or become a member function of the discretization class

  if (!dis->HaveDofs()) dserror("FillComplete() was not called");
  const Epetra_Map* colmap = dis->DofColMap();
  const Epetra_BlockMap& vecmap = state->Map();

  // if it's already in column map just set a reference
  // This is a rought test, but it might be ok at this place. It is an
  // error anyway to hand in a vector that is not related to our dof
  // maps.
  if (vecmap.PointSameAs(*colmap))
    return state;
  // if it's not in column map export and allocate
  else
  {
    Teuchos::RCP<Epetra_Vector> tmp = LINALG::CreateVector(*colmap,false);
    LINALG::Export(*state,*tmp);
    return tmp;
  }
}

/*----------------------------------------------------------------------*
 |(private)                                                   tk 06/10  |
 |recompute nodecolmap of standard discretization to include all        |
 |nodes as of subdicretization                                          |
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Map> DRT::UTILS::ComputeNodeColMap(
         const RCP<DRT::Discretization> sourcedis,  ///< standard discretization we want to redistribute
         const RCP<DRT::Discretization> subdis ///< subdiscretization prescribing ghosting
         )
{
  const Epetra_Map* oldcolnodemap = sourcedis->NodeColMap();

  vector<int> mycolnodes(oldcolnodemap->NumMyElements());
  oldcolnodemap->MyGlobalElements (&mycolnodes[0]);
  for (int inode = 0; inode != subdis->NumMyColNodes(); ++inode)
  {
      const DRT::Node* newnode = subdis->lColNode(inode);
      const int gid = newnode->Id();
      if (!(sourcedis->HaveGlobalNode(gid)))
      {
          mycolnodes.push_back(gid);
      }
  }

  // now reconstruct the extended colmap
  RCP<Epetra_Map> newcolnodemap = Teuchos::rcp(new Epetra_Map(-1,
                                     mycolnodes.size(),
                                     &mycolnodes[0],
                                     0,
                                     sourcedis->Comm()));
  return newcolnodemap;
}


