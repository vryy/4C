/*!----------------------------------------------------------------------
\file drt_utils.cpp
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

#ifdef CCADISCRET

#ifdef PARALLEL
#ifdef LINUX_MUENCH
typedef int idxtype;
extern "C"
{
  void METIS_PartGraphRecursive(int *, idxtype *, idxtype *,
				idxtype *, idxtype *, int *, int *, int *,
				int *, int *, idxtype *);
  void METIS_PartGraphKway(int *, idxtype *, idxtype *, idxtype *,
			   idxtype *, int *, int *, int *, int *, int *,
			   idxtype *);
}
#endif
#include <mpi.h>
#endif // PARALLEL

#include <algorithm>
#include <numeric>
#include <vector>
#include <set>
#include <map>
#include <string>

#include "../linalg/linalg_utils.H"
#include "../linalg/linalg_mapextractor.H"

#include "drt_utils.H"
#include "drt_node.H"
#include "../drt_nurbs_discret/drt_control_point.H"
#include "drt_dofset.H"
#include "drt_discret.H"

#include "drt_dserror.H"
#include "standardtypes_cpp.H"

#include "drt_parobjectfactory.H"
#include "../drt_s8/shell8.H"

/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;


/*----------------------------------------------------------------------*
 |  allocate an instance of a specific impl. of ParObject (public) mwgee 12/06|
 *----------------------------------------------------------------------*/
DRT::ParObject* DRT::UTILS::Factory(const vector<char>& data)
{
  return ParObjectFactory::Instance().Create( data );
}

/*----------------------------------------------------------------------*
 |  allocate an element of a specific type (public)          mwgee 03|07|
 *----------------------------------------------------------------------*/
Teuchos::RCP<DRT::Element> DRT::UTILS::Factory(const string eletype,
                                              const string eledistype,
                                              const int id,
                                              const int owner)
{
  return ParObjectFactory::Instance().Create( eletype, eledistype, id, owner );
}


/*----------------------------------------------------------------------*
 |  partition a graph using metis  (public)                  mwgee 11/06|
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_CrsGraph> DRT::UTILS::PartGraphUsingMetis(
                                             const Epetra_CrsGraph& graph,
                                             const Epetra_Vector& weights)
{
  const int myrank   = graph.Comm().MyPID();
  const int numproc  = graph.Comm().NumProc();

  if (numproc==1)
  {
    Teuchos::RCP<Epetra_CrsGraph> outgraph = rcp(new Epetra_CrsGraph(graph));
    return outgraph;
  }

  // proc that will do the serial partitioning
  // the graph is collapsed to this proc
  // Normally this would be proc 0 but 0 always has so much to do.... ;-)
  int workrank = 0;
  if (graph.Comm().NumProc()>1) workrank=1;

  // get rowmap of the graph
  const Epetra_BlockMap& tmp = graph.RowMap();
  Epetra_Map rowmap(tmp.NumGlobalElements(),tmp.NumMyElements(),
                    tmp.MyGlobalElements(),0,graph.Comm());

  // build a target map that stores everything on proc workrank
  // We have arbirtary gids here and we do not tell metis about
  // them. So we have to keep rowrecv until the redistributed map is
  // build.


  // rowrecv is a fully redundant vector (size of number of nodes)
  vector<int> rowrecv(rowmap.NumGlobalElements());

  // after Allreduce rowrecv contains
  //
  // *-+-+-    -+-+-*-+-+-    -+-+-*-           -*-+-+-    -+-+-*
  // * | | .... | | * | | .... | | * ..........  * | | .... | | *
  // *-+-+-    -+-+-*-+-+-    -+-+-*-           -*-+-+-    -+-+-*
  //   gids stored     gids stored                  gids stored
  //  on first proc  on second proc                 on last proc
  //
  // the ordering of the gids on the procs is arbitrary (as copied
  // from the map)
  LINALG::AllreduceEMap(rowrecv, rowmap);

  // construct an epetra map from the list of gids
  Epetra_Map tmap(rowmap.NumGlobalElements(),
                  // if ..........    then ............... else
                  (myrank == workrank) ? (int)rowrecv.size() : 0,
                  &rowrecv[0],
                  0,
                  rowmap.Comm());

  // export the graph to tmap
  Epetra_CrsGraph tgraph(Copy,tmap,108,false);
  Epetra_Export exporter(rowmap,tmap);
  int err = tgraph.Export(graph,exporter,Add);
  if (err<0) dserror("Graph export returned err=%d",err);
  tgraph.FillComplete();
  tgraph.OptimizeStorage();

  // export the weights to tmap
  Epetra_Vector tweights(tmap,false);
  err = tweights.Export(weights,exporter,Insert);
  if (err<0) dserror("Vector export returned err=%d",err);

  // do partitioning using metis on workrank
  vector<int> part(tmap.NumMyElements());
  if (myrank==workrank)
  {
    // metis requests indexes. So we need a reverse lookup from gids
    // to indexes.
    map<int,int> idxmap;
    for (unsigned i=0; i<rowrecv.size(); ++i)
    {
      idxmap[rowrecv[i]] = i;
    }

    // xadj points from index i to the index of the
    // first adjacent node
    vector<int> xadj(tmap.NumMyElements()+1);

    // a list of adjacent nodes, adressed using xadj
    vector<int> adjncy(tgraph.NumGlobalNonzeros()); // the size is an upper bound

    // rowrecv(i)       rowrecv(i+1)                      node gids
    //     ^                 ^
    //     |                 |
    //     | idxmap          | idxmap
    //     |                 |
    //     v                 v
    //     i                i+1                       equivalent indices
    //     |                 |
    //     | xadj            | xadj
    //     |                 |
    //     v                 v
    //    +-+-+-+-+-+-+-+-+-+-+                -+-+-+
    //    | | | | | | | | | | | ............... | | |      adjncy
    //    +-+-+-+-+-+-+-+-+-+-+                -+-+-+
    //
    //    |       i's       |    (i+1)'s
    //    |    neighbours   |   neighbours           (numbered by equivalent indices)
    //

    vector<int> vwgt(tweights.MyLength());
    for (int i=0; i<tweights.MyLength(); ++i) vwgt[i] = (int)tweights[i];

    int count=0;
    xadj[0] = 0;
    for (int row=0; row<tgraph.NumMyRows(); ++row)
    {
      //cout << "xadj[" << row << "] = " << xadj[row] << endl;
      int grid = tgraph.RowMap().GID(row);
      int numindices;
      int* lindices;
      int err = tgraph.ExtractMyRowView(row,numindices,lindices);
      if (err) dserror("Epetra_CrsGraph::ExtractMyRowView returned err=%d",err);
      //cout << "adjncy: ";
      for (int col=0; col<numindices; ++col)
      {
        int gcid = tgraph.ColMap().GID(lindices[col]);
        if (gcid==grid) continue;
        adjncy[count] = idxmap[gcid];
        //cout << adjncy[count] << " ";
        ++count;
      }
      //cout << endl;
      xadj[row+1] = count;
    }
    //cout << "xadj[" << xadj.size()-1 << "] = " << xadj[xadj.size()-1] << endl;
    //cout << "tgraph.NumGlobalNonzeros() " << tgraph.NumGlobalNonzeros() << endl
    //     << "tmap.NumMyElements()       " << tmap.NumMyElements() << endl
    //     << "count                      " << count << endl;

    idxmap.clear();

    if (numproc<8) // better for smaller no. of partitions
    {
#ifdef PARALLEL
      int wgtflag=2;
      int numflag=0;
      int npart=numproc;
      int options[5] = { 0,3,1,1,0 };
      int edgecut=0;
      int nummyele = tmap.NumMyElements();
      METIS_PartGraphRecursive(&nummyele,
                               &xadj[0],
                               &adjncy[0],
                               &vwgt[0],
                               NULL,
                               &wgtflag,
                               &numflag,
                               &npart,
                               options,
                               &edgecut,
                               &part[0]);
#endif
    }
    else
    {
#ifdef PARALLEL
      int wgtflag=2;
      int numflag=0;
      int npart=numproc;
      int options[5] = { 0,3,1,1,0 };
      int edgecut=0;
      int nummyele = tmap.NumMyElements();
      METIS_PartGraphKway(&nummyele,
                          &xadj[0],
                          &adjncy[0],
                          &vwgt[0],
                          NULL,
                          &wgtflag,
                          &numflag,
                          &npart,
                          options,
                          &edgecut,
                          &part[0]);
#endif
    }
  } // if (myrank==workrank)

  // broadcast partitioning result
  int size = tmap.NumMyElements();
  tmap.Comm().Broadcast(&size,1,workrank);
  part.resize(size);
  tmap.Comm().Broadcast(&part[0],size,workrank);

  // loop part and count no. of nodes belonging to me
  // (we reuse part to save on memory)
  int count=0;
  for (int i=0; i<size; ++i)
    if (part[i]==myrank)
    {
      part[count] = rowrecv[i];
      ++count;
    }

  rowrecv.clear();

  // create map with new layout
  Epetra_Map newmap(size,count,&part[0],0,graph.Comm());

  // create the output graph and export to it
  Teuchos::RCP<Epetra_CrsGraph> outgraph =
                           rcp(new Epetra_CrsGraph(Copy,newmap,108,false));
  Epetra_Export exporter2(graph.RowMap(),newmap);
  err = outgraph->Export(graph,exporter2,Add);
  if (err<0) dserror("Graph export returned err=%d",err);
  outgraph->FillComplete();
  outgraph->OptimizeStorage();

  return outgraph;
}



/*----------------------------------------------------------------------*
 |  locally extract a subset of values  (public)            mwgee 12/06|
 *----------------------------------------------------------------------*/
void DRT::UTILS::ExtractMyValues(const Epetra_Vector& global,
                                 vector<double>& local,
                                 const vector<int>& lm)
{
  const size_t ldim = lm.size();
  local.resize(ldim);
  for (size_t i=0; i<ldim; ++i)
  {
    const int lid = global.Map().LID(lm[i]);
    if (lid<0)
      dserror("Proc %d: Cannot find gid=%d in Epetra_Vector",global.Comm().MyPID(),lm[i]);
    local[i] = global[lid];
  }
  return;
}


/*----------------------------------------------------------------------*
 |  locally extract a subset of values  (public)             henke 12/09|
 *----------------------------------------------------------------------*/
void DRT::UTILS::ExtractMyValues(const Epetra_Vector&      global,
                                 Epetra_SerialDenseVector& local,
                                 const vector<int>&        lm)
{
  const size_t ldim = lm.size();
  local.Size(ldim);
  for (size_t i=0; i<ldim; ++i)
  {
    const int lid = global.Map().LID(lm[i]);
    if (lid<0) dserror("Proc %d: Cannot find gid=%d in Epetra_Vector",global.Comm().MyPID(),lm[i]);
    local[i] = global[lid];
  }
  return;
}


/*----------------------------------------------------------------------*
 | extract local values from global node-based (multi) vector           |
 |                                                          henke 06/09 |
 *----------------------------------------------------------------------*/
void DRT::UTILS::ExtractMyNodeBasedValues(
    const DRT::Element* ele,
    std::vector<double>& local,
    const Epetra_MultiVector& global)
{
  const int numnode = ele->NumNode();
  const int numcol = global.NumVectors();
  local.resize(numnode*numcol);

  // loop over element nodes
  for (int i=0; i<numnode; ++i)
  {
    const int nodegid = (ele->Nodes()[i])->Id();
    const int lid = global.Map().LID(nodegid);
    if (lid<0) dserror("Proc %d: Cannot find gid=%d in Epetra_Vector",global.Comm().MyPID(),nodegid);

    // loop over multi vector columns (numcol=1 for Epetra_Vector)
    for (int col=0; col<numcol; col++)
    {
      double* globalcolumn = (global)[col];
      local[col+(numcol*i)] = globalcolumn[lid];
    }
  }
  return;
}


/*----------------------------------------------------------------------*
 | extract local values from global node-based multi vector   gjb 08/08 |
 *----------------------------------------------------------------------*/
void DRT::UTILS::ExtractMyNodeBasedValues(
    const DRT::Element* ele,
    Epetra_SerialDenseVector& local,
    const RCP<Epetra_MultiVector>& global,
    const int nsd
    )
{
  if (global==null) dserror("received a TEUCHOS::null pointer");
  if (nsd > global->NumVectors())
    dserror("Requested %d of %d available columns", nsd,global->NumVectors());
  const int iel = ele->NumNode(); // number of nodes
  if (local.Length()!=(iel*nsd)) dserror("vector size mismatch.");

  for (int i=0; i<nsd; i++)
  {
    // access actual component column of multi-vector
    double* globalcolumn = (*global)[i];
    // loop over the element nodes
    for (int j=0;j<iel;j++)
    {
      const int nodegid = (ele->Nodes()[j])->Id();
      const int lid = global->Map().LID(nodegid);
      local(i+(nsd*j))=globalcolumn[lid];
    }
  }
  return;
}


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
  RCP<Epetra_Map> newcolnodemap = rcp(new Epetra_Map(-1,
                                     mycolnodes.size(),
                                     &mycolnodes[0],
                                     0,
                                     sourcedis->Comm()));
  return newcolnodemap;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::UTILS::ComputeStructure3DNullSpace( DRT::Discretization & dis, std::vector<double> & ns, const double * x0, int numdf, int dimns )
{
  const Epetra_Map* rowmap = dis.DofRowMap();
  const int lrows = rowmap->NumMyElements();
  double* mode[6];
  for (int i=0; i<dimns; ++i) mode[i] = &(ns[i*lrows]);

#ifdef D_SHELL8
  // special for shell8
  Epetra_SerialDenseMatrix dir;
  DRT::Element* dwele = dis.lRowElement(0);
  if (dwele->ElementType()==DRT::ELEMENTS::Shell8Type::Instance())
  {
    dir.Shape(dis.NumMyRowNodes(),3);
    for (int i=0; i<dis.NumMyRowNodes(); ++i)
    {
      DRT::Node* actnode = dis.lRowNode(i);
      DRT::ELEMENTS::Shell8* s8 =
        dynamic_cast<DRT::ELEMENTS::Shell8*>(actnode->Elements()[0]);
      if (!s8) dserror("Cannot cast to Shell8");
      int j;
      for (j=0; j<s8->NumNode(); ++j)
        if (s8->Nodes()[j]->Id() == actnode->Id()) break;
      if (j==s8->NumNode()) dserror("Can't find matching node - weird!");
      double h2 = (*s8->GetThickness())[j]/2.0;
      // get director
      const Epetra_SerialDenseMatrix* a3ref = s8->GetDirectors();
      dir(i,0) = (*a3ref)(0,j)*h2;
      dir(i,1) = (*a3ref)(1,j)*h2;
      dir(i,2) = (*a3ref)(2,j)*h2;
    }
  }
#endif

  /* the rigid body modes for structures are:
        xtrans   ytrans  ztrans   xrot       yrot       zrot
        mode[0]  mode[1] mode[2]  mode[3]    mode[4]    mode[5]
  -----------------------------------------------------------
  x   |    1       0       0       0          z-z0      -y+y0
  y   |    0       1       0      -z+z0       0          x-x0
  z   |    0       0       1       y-y0      -x+x0       0
  dx  |    0       0       0       0          a3        -a2
  dy  |    0       0       0      -a3         0          a1
  dz  |    0       0       0       a2        -a1         0
  */

  // works straight for bricks as well
//   if (ele->Type() == DRT::Element::element_shell8 ||
//       ele->Type() == DRT::Element::element_ale3 ||
//       ele->Type() == DRT::Element::element_so_hex8 ||
//       ele->Type() == DRT::Element::element_so_hex20 ||
//       ele->Type() == DRT::Element::element_so_hex27 ||
//       ele->Type() == DRT::Element::element_sosh8 ||
//       ele->Type() == DRT::Element::element_so_tet4 ||
//       ele->Type() == DRT::Element::element_so_tet10 ||
//       ele->Type() == DRT::Element::element_so_weg6 ||
//       ele->Type() == DRT::Element::element_sodisp ||
//       ele->Type() == DRT::Element::element_so_shw6 ||
//       ele->Type() == DRT::Element::element_truss3 ||
//       ele->Type() == DRT::Element::element_torsion3)
  {
    for (int i=0; i<dis.NumMyRowNodes(); ++i)
    {
      DRT::Node* actnode = dis.lRowNode(i);
      const double* x = actnode->X();
      vector<int> dofs = dis.Dof(actnode);
      for (unsigned j=0; j<dofs.size(); ++j)
      {
        const int dof = dofs[j];
        const int lid = rowmap->LID(dof);
        if (lid<0) dserror("Cannot find dof");
        switch (j) // j is degree of freedom
        {
        case 0:
          mode[0][lid] = 1.0;
          mode[1][lid] = 0.0;
          mode[2][lid] = 0.0;
          mode[3][lid] = 0.0;
          mode[4][lid] = x[2] - x0[2];
          mode[5][lid] = -x[1] + x0[1];
        break;
        case 1:
          mode[0][lid] = 0.0;
          mode[1][lid] = 1.0;
          mode[2][lid] = 0.0;
          mode[3][lid] = -x[2] + x0[2];
          mode[4][lid] = 0.0;
          mode[5][lid] = x[0] - x0[0];
        break;
        case 2:
          mode[0][lid] = 0.0;
          mode[1][lid] = 0.0;
          mode[2][lid] = 1.0;
          mode[3][lid] = x[1] - x0[1];
          mode[4][lid] = -x[0] + x0[0];
          mode[5][lid] = 0.0;
        break;
#ifdef D_SHELL8
        case 3:
          mode[0][lid] = 0.0;
          mode[1][lid] = 0.0;
          mode[2][lid] = 0.0;
          mode[3][lid] = 0.0;
          mode[4][lid] = dir(i,2);
          mode[5][lid] = -dir(i,1);
        break;
        case 4:
          mode[0][lid] = 0.0;
          mode[1][lid] = 0.0;
          mode[2][lid] = 0.0;
          mode[3][lid] = -dir(i,2);
          mode[4][lid] = 0.0;
          mode[5][lid] = dir(i,0);
        break;
        case 5:
          mode[0][lid] = 0.0;
          mode[1][lid] = 0.0;
          mode[2][lid] = 0.0;
          mode[3][lid] = dir(i,1);
          mode[4][lid] = -dir(i,0);
          mode[5][lid] = 0.0;
        break;
#endif
        default:
          dserror("Only dofs 0 - 5 supported");
        break;
        } // switch (j)
      } // for (int j=0; j<actnode->Dof().NumDof(); ++j)
    } // for (int i=0; i<NumMyRowNodes(); ++i)
  } // if (ele->Type() == DRT::Element::element_shell8)
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::UTILS::ComputeStructure2DNullSpace( DRT::Discretization & dis, std::vector<double> & ns, const double * x0, int numdf, int dimns )
{
  const Epetra_Map* rowmap = dis.DofRowMap();
  const int lrows = rowmap->NumMyElements();
  double* mode[6];
  for (int i=0; i<dimns; ++i) mode[i] = &(ns[i*lrows]);

//   else if (ele->Type() == DRT::Element::element_wall1 ||
//            ele->Type() == DRT::Element::element_ale2 ||
//            ele->Type() == DRT::Element::element_torsion2)
  {
    for (int i=0; i<dis.NumMyRowNodes(); ++i)
    {
      DRT::Node* actnode = dis.lRowNode(i);
      const double* x = actnode->X();
      vector<int> dofs = dis.Dof(actnode);
      for (unsigned j=0; j<dofs.size(); ++j)
      {
        const int dof = dofs[j];
        const int lid = rowmap->LID(dof);
        if (lid<0) dserror("Cannot find dof");
        switch (j) // j is degree of freedom
        {
        case 0:
          mode[0][lid] = 1.0;
          mode[1][lid] = 0.0;
          mode[2][lid] = -x[1] + x0[1];
        break;
        case 1:
          mode[0][lid] = 0.0;
          mode[1][lid] = 1.0;
          mode[2][lid] = x[0] - x0[0];
        break;
        default:
          dserror("Only dofs 0 - 1 supported");
        break;
        } // switch (j)
      } // for (int j=0; j<actnode->Dof().NumDof(); ++j)
    } // for (int i=0; i<NumMyRowNodes(); ++i)
  } // else if (ele->Type() == DRT::Element::element_wall1)
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::UTILS::ComputeBeam2DNullSpace( DRT::Discretization & dis, std::vector<double> & ns, const double * x0, int numdf, int dimns )
{
  const Epetra_Map* rowmap = dis.DofRowMap();
  const int lrows = rowmap->NumMyElements();
  double* mode[6];
  for (int i=0; i<dimns; ++i) mode[i] = &(ns[i*lrows]);

  /* the variable mode[i][j] describes how the position of a
   * node changes with respect to the j-th degree of freedom
   * in case that the i-th rigid body mode is applied to the
   * structure; the structure altogether always has 3 rigid body
   * modes in R^2 and 6 in R^3; these modes are translation in
   * each coordinate direction, respectively, and rotation
   * around each axis, respectively. This is summed up in the
   * following table where in the left column x,y,z denote
   * translations in x-, y- and z-direction of a node due to
   * the application of a rigid body mode, whereas dx,dy,dz
   * denote increments of the node's rotational degrees of
   * freedom, which relate to a rotation around the x-,y-
   * and z-axis.
   *
        xtrans   ytrans  ztrans   xrot       yrot       zrot
        mode[0]  mode[1] mode[2]  mode[3]    mode[4]    mode[5]
  -----------------------------------------------------------
  x   |    1       0       0       0          z-z0      -y+y0
  y   |    0       1       0      -z+z0       0          x-x0
  z   |    0       0       1       y-y0      -x+x0       0
  dx  |    0       0       0       1          0          0
  dy  |    0       0       0       0          1          0
  dz  |    0       0       0       0          0          1

  for example the first line means: a translation of a node in
  x-direction may be caused either by a x-translation of the whole
  structure (which is rigid body mode 0) or by a rotation either
  around the y-axis or the z-axis. In case of a rotation dtheta around the
  y-axis for example the resulting x-translation is dtheta times the
  lever arm z - z0. Here z0 represents the z-coordinate of the point
  around which the structure is rotated. This point may be chosen
  arbitrarily and by the algorithms underlying to this method it is
  chosen automatically according to some mathematical considerations.
  Note that this holds true for infinitesimal rotations dtehta, only, of
  course.
  On the other hand e.g. the fourth column means that a rigid body
  rotation dtheta around the x-axis entails translations (-z+z0)*dtheta
  in y-direction and (y-y0)*dtheta in z-direction and a rotation
  increment 1*dtheta of the rotational degree of freedom related to the
  x-axis.
  */

  /* for beam2 and beam2r elements the above table reduces to
   *
        xtrans   ytrans    zrot
        mode[0]  mode[1]   mode[2]
  -----------------------------------------------------------
  x   |    1       0       -y+y0
  y   |    0       1       x-x0
  dz  |    0       0       1
  note: for the here employed Timoshenko and Reissner beam elements a rigid
  body rotation entails also an increment of the rotation degree of freedom
  dz which makes the director of the beam move accordingly; only then
  a rotation does not lead to any shear stress and is truely a rigid
  body rotation
  */

  /*two dimensional beam beam elements, where each node has 2 translational
   * degrees of freedom and one rotational degree of freedom*/
//   else if (ele->Type() == DRT::Element::element_beam2 ||
//            ele->Type() == DRT::Element::element_beam2r)
  {
    //looping through all nodes
    for (int i=0; i<dis.NumMyRowNodes(); ++i)
    {
      //getting pointer at current node
      DRT::Node* actnode = dis.lRowNode(i);

      //getting coordinates of current node
      const double* x = actnode->X();

      //getting number of degrees of freedom of current node
      vector<int> dofs = dis.Dof(actnode);

      //looping through all degrees of freedom of a node
      for (unsigned j=0; j<dofs.size(); ++j)
      {
        const int dof = dofs[j];
        const int lid = rowmap->LID(dof);
        if (lid<0) dserror("Cannot find dof");
        // j is degree of freedom; each case refers to one line in the above table
        switch (j)
        {
        case 0:
          mode[0][lid] = 1.0;
          mode[1][lid] = 0.0;
          mode[2][lid] = -x[1] + x0[1];
        break;
        case 1:
          mode[0][lid] = 0.0;
          mode[1][lid] = 1.0;
          mode[2][lid] = x[0] - x0[0];
        break;
        case 2:
          mode[0][lid] = 0.0;
          mode[1][lid] = 0.0;
          mode[2][lid] = 1.0;
        break;
        default:
          dserror("Only dofs 0 - 2 supported");
        break;
        } // switch (j)
      } // for (int j=0; j<actnode->Dof().NumDof(); ++j)
    } // for (int i=0; i<NumMyRowNodes(); ++i)
  } // else if (ele->Type() == DRT::Element::element_beam2 || ele->Type() == DRT::Element::element_beam2r)
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::UTILS::ComputeBeam3DNullSpace( DRT::Discretization & dis, std::vector<double> & ns, const double * x0, int numdf, int dimns )
{
  const Epetra_Map* rowmap = dis.DofRowMap();
  const int lrows = rowmap->NumMyElements();
  double* mode[6];
  for (int i=0; i<dimns; ++i) mode[i] = &(ns[i*lrows]);

  /* for beam3 elements the relation between rigid body modes and
   * increments on the degrees of freedom is non-trivial since
   * rotational increments in 3D are non-additive in general. In
   * general this relation may require calling all the elements.
   * However, in opposition to the SHELL8 element it is not
   * sufficient to just call a director saved in the element.
   * Rather to calculate proper increments for the rotational
   * degrees of freedom due to a rigid body rotation of the
   * complete structure, the triad at each node is required in
   * order to transform non-additive increments into additive ones.
   * However, the beam3 element currently does not save the nodal
   * triads as a class variable, but only the triads at each Gauss
   * point. In the following a wrong (!!!) dummy version is implemneted
   * but commented out. In this dummy version the rotational degrees of
   * freedom are treated identically to the additive translational
   * degrees of freedom. Activating and using this part of the code
   * quickly reveals the problems of such a naive implemnetation.
   * Usually the equation solver simply does not work with this
   * dummy code, i.e. the iterative solution process does not converge.
   * If Algebraic Multigrid methods should be really used for beam3
   * elements, one first has to develop efficient special methods for
   * these elements. Currently trying to use Algebraic multigrid methods
   * for beam3 elements just amounts to an error as no properly working
   * implementation has been available so far*/

//   else if (ele->Type() == DRT::Element::element_beam3 || ele->Type() == DRT::Element::element_beam3ii)
  {
    //looping through all nodes
    for (int i=0; i<dis.NumMyRowNodes(); ++i)
    {
      //getting pointer at current node
      DRT::Node* actnode = dis.lRowNode(i);

      //getting coordinates of current node
      const double* x = actnode->X();

      //getting number of degrees of freedom of current node
      vector<int> dofs = dis.Dof(actnode);

      //looping through all degrees of freedom of a node
      for (unsigned j=0; j<dofs.size(); ++j)
      {
        const int dof = dofs[j];
        const int lid = rowmap->LID(dof);
        if (lid<0) dserror("Cannot find dof");
        // j is degree of freedom; each case refers to one line in the above table
        switch (j)
        {
        case 0:
          mode[0][lid] = 1.0;
          mode[1][lid] = 0.0;
          mode[2][lid] = 0.0;
          mode[3][lid] = 0.0;
          mode[4][lid] =  x[2] - x0[2];
          mode[5][lid] = -x[1] + x0[1];
        break;
        case 1:
          mode[0][lid] = 0.0;
          mode[1][lid] = 1.0;
          mode[2][lid] = 0.0;
          mode[3][lid] = -x[2] + x0[2];
          mode[4][lid] = 0.0;
          mode[5][lid] =  x[0] - x0[0];
        break;
        case 2:
          mode[0][lid] = 0.0;
          mode[1][lid] = 0.0;
          mode[2][lid] = 1.0;
          mode[3][lid] =  x[1] - x0[1];
          mode[4][lid] = -x[0] + x0[0];
          mode[5][lid] = 0.0;
        break;
        case 3:
          mode[0][lid] = 0.0;
          mode[1][lid] = 0.0;
          mode[2][lid] = 0.0;
          mode[3][lid] = 1.0;
          mode[4][lid] = 0.0;
          mode[5][lid] = 0.0;
        break;
        case 4:
          mode[0][lid] = 0.0;
          mode[1][lid] = 0.0;
          mode[2][lid] = 0.0;
          mode[3][lid] = 0.0;
          mode[4][lid] = 1.0;
          mode[5][lid] = 0.0;
        break;
        case 5:
          mode[0][lid] = 0.0;
          mode[1][lid] = 0.0;
          mode[2][lid] = 0.0;
          mode[3][lid] = 0.0;
          mode[4][lid] = 0.0;
          mode[5][lid] = 1.0;
        break;
        default:
          dserror("Only dofs 0 - 5 supported");
        break;
        } // switch (j)
      } // for (int j=0; j<actnode->Dof().NumDof(); ++j)
    } // for (int i=0; i<NumMyRowNodes(); ++i)

  } // else if (ele->Type() == DRT::Element::element_beam3 || ele->Type() == DRT::Element::element_beam3ii)
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::UTILS::ComputeXFluid3DNullSpace( DRT::Discretization & dis, std::vector<double> & ns, const double * x0, int numdf, int dimns )
{
  const Epetra_Map* rowmap = dis.DofRowMap();
  const int lrows = rowmap->NumMyElements();
  double* mode[6];
  for (int i=0; i<dimns; ++i) mode[i] = &(ns[i*lrows]);

  /* the rigid body modes for fluids are:
        xtrans   ytrans  ztrans   pressure
        mode[0]  mode[1] mode[2]  mode[3]
  ----------------------------------------
  x   |    1       0       0       0
  y   |    0       1       0       0
  z   |    0       0       1       0
  p   |    0       0       0       1
  */

//   else if (ele->Type() == DRT::Element::element_xfluid3 ||
//            ele->Type() == DRT::Element::element_combust3 ||
//            ele->Type() == DRT::Element::element_smoothrod ||
//            ele->Type() == DRT::Element::element_sosh8p8)
  {
    for (int i=0; i<dis.NumMyRowNodes(); ++i)
    {
      DRT::Node* actnode = dis.lRowNode(i);
      vector<int> dofs = dis.Dof(actnode);
      for (unsigned j=0; j<dofs.size(); ++j)
      {
        const int dof = dofs[j];
        const int lid = rowmap->LID(dof);
        if (lid<0) dserror("Cannot find dof");
        switch (j) // j is degree of freedom
        {
        case 0:
          mode[0][lid] = 1.0;
          mode[1][lid] = 0.0;
          mode[2][lid] = 0.0;
          mode[3][lid] = 0.0;
        break;
        case 1:
          mode[0][lid] = 0.0;
          mode[1][lid] = 1.0;
          mode[2][lid] = 0.0;
          mode[3][lid] = 0.0;
        break;
        case 2:
          mode[0][lid] = 0.0;
          mode[1][lid] = 0.0;
          mode[2][lid] = 1.0;
          mode[3][lid] = 0.0;
        break;
        case 3:
          mode[0][lid] = 0.0;
          mode[1][lid] = 0.0;
          mode[2][lid] = 0.0;
          mode[3][lid] = 1.0;
        break;
        case 4:
          mode[0][lid] = 0.0;
          mode[1][lid] = 0.0;
          mode[2][lid] = 0.0;
          mode[3][lid] = 0.0;
        break;
        case 5:
          mode[0][lid] = 0.0;
          mode[1][lid] = 0.0;
          mode[2][lid] = 0.0;
          mode[3][lid] = 0.0;
        break;
        case 6:
          mode[0][lid] = 0.0;
          mode[1][lid] = 0.0;
          mode[2][lid] = 0.0;
          mode[3][lid] = 0.0;
        break;
        case 7:
          mode[0][lid] = 0.0;
          mode[1][lid] = 0.0;
          mode[2][lid] = 0.0;
          mode[3][lid] = 0.0;
        break;
        default:
          dserror("Only dofs 0 - 7 supported");
        break;
        } // switch (j)
      } // for (int j=0; j<actnode->Dof().NumDof(); ++j)
    } // for (int i=0; i<NumMyRowNodes(); ++i)
  } // else if (ele->Type() == DRT::Element::element_fluid3 or ele->Type() == DRT::Element::element_xfluid3)
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::UTILS::ComputeFluid2DNullSpace( DRT::Discretization & dis, std::vector<double> & ns, const double * x0, int numdf, int dimns )
{
  const Epetra_Map* rowmap = dis.DofRowMap();
  const int lrows = rowmap->NumMyElements();
  double* mode[6];
  for (int i=0; i<dimns; ++i) mode[i] = &(ns[i*lrows]);

//   else if (ele->Type() == DRT::Element::element_fluid2)
  {
    for (int i=0; i<dis.NumMyRowNodes(); ++i)
    {
      DRT::Node* actnode = dis.lRowNode(i);
      vector<int> dofs = dis.Dof(actnode);
      for (unsigned j=0; j<dofs.size(); ++j)
      {
        const int dof = dofs[j];
        const int lid = rowmap->LID(dof);
        if (lid<0) dserror("Cannot find dof");
        switch (j) // j is degree of freedom
        {
        case 0:
          mode[0][lid] = 1.0;
          mode[1][lid] = 0.0;
          mode[2][lid] = 0.0;
        break;
        case 1:
          mode[0][lid] = 0.0;
          mode[1][lid] = 1.0;
          mode[2][lid] = 0.0;
        break;
        case 2:
          mode[0][lid] = 0.0;
          mode[1][lid] = 0.0;
          mode[2][lid] = 1.0;
        break;
        default:
          dserror("Only dofs 0 - 2 supported");
        break;
        } // switch (j)
      } // for (int j=0; j<actnode->Dof().NumDof(); ++j)
    } // for (int i=0; i<NumMyRowNodes(); ++i)
  } // else if (ele->Type() == DRT::Element::element_fluid2)
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::UTILS::ComputeFluid3DNullSpace( DRT::Discretization & dis, std::vector<double> & ns, const double * x0, int numdf, int dimns )
{
  const Epetra_Map* rowmap = dis.DofRowMap();
  const int lrows = rowmap->NumMyElements();
  double* mode[6];
  for (int i=0; i<dimns; ++i) mode[i] = &(ns[i*lrows]);

//   else if (ele->Type() == DRT::Element::element_transport or
//       ele->Type() == DRT::Element::element_fluid3)
  {
    for (int i=0; i<dis.NumMyRowNodes(); ++i)
    {
      DRT::Node* actnode = dis.lRowNode(i);
      vector<int> dofs = dis.Dof(actnode);
      const unsigned int ndof = dofs.size();
      for (unsigned j=0; j<ndof; ++j)
      {
        const int dof = dofs[j];
        const int lid = rowmap->LID(dof);
        if (lid<0) dserror("Cannot find dof");

        for (unsigned k=0; k<ndof; ++k)
        {
          if (k == j)
            mode[k][lid] = 1.0;
          else
            mode[k][lid] = 0.0;
        }

      } // for (int j=0; j<actnode->Dof().NumDof(); ++j)
    } // for (int i=0; i<NumMyRowNodes(); ++i)
  } // else if (ele->Type() == DRT::Element::element_transport)
}


#endif  // #ifdef CCADISCRET
