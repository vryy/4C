/*!----------------------------------------------------------------------
\file drt_utils_metis.cpp

\brief A collection of helper methods for namespace DRT

<pre>
\level 2

\maintainer Martin Kronbichler
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>

*----------------------------------------------------------------------*/

#include "drt_utils_metis.H"
#include "../linalg/linalg_utils.H"

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
#endif // PARALLEL


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
    Teuchos::RCP<Epetra_CrsGraph> outgraph = Teuchos::rcp(new Epetra_CrsGraph(graph));
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
  std::vector<int> rowrecv(rowmap.NumGlobalElements());

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
  std::vector<int> part(tmap.NumMyElements());
  if (myrank==workrank)
  {
    // metis requests indexes. So we need a reverse lookup from gids
    // to indexes.
    std::map<int,int> idxmap;
    for (unsigned i=0; i<rowrecv.size(); ++i)
    {
      idxmap[rowrecv[i]] = i;
    }

    // xadj points from index i to the index of the
    // first adjacent node
    std::vector<int> xadj(tmap.NumMyElements()+1);

    // a list of adjacent nodes, adressed using xadj
    std::vector<int> adjncy(tgraph.NumGlobalNonzeros()); // the size is an upper bound

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

    std::vector<int> vwgt(tweights.MyLength());
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
                           Teuchos::rcp(new Epetra_CrsGraph(Copy,newmap,108,false));
  Epetra_Export exporter2(graph.RowMap(),newmap);
  err = outgraph->Export(graph,exporter2,Add);
  if (err<0) dserror("Graph export returned err=%d",err);
  outgraph->FillComplete();
  outgraph->OptimizeStorage();

  return outgraph;
}


#if !defined(PARALLEL) || !defined(HAVE_PARMETIS)

#include "drt_inputreader.H"
#include "drt_utils.H"
#include "../linalg/linalg_utils.H"
#include "standardtypes_cpp.H"

#include <Epetra_Time.h>
#include <iterator>

#if 0 // we are not using metis anymore, use parmetis instead
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::UTILS::PartUsingMetis(Teuchos::RCP<Epetra_Map>& rownodes,
                                Teuchos::RCP<Epetra_Map>& colnodes,
                                list<std::vector<int> >& elementnodes,
                                Teuchos::RCP<Epetra_Comm> comm)
{

#if 0
    if (myrank==0)
    {
      cout << "\n\nelementnodes: size=" << elementnodes.size() << endl;
      for (list<std::vector<int> >::iterator i=elementnodes.begin();
           i!=elementnodes.end();
           ++i)
      {
        copy(i->begin(), i->end(), ostream_iterator<int>(cout, " "));
        cout << endl;
      }
    }
#endif

#ifdef WE_DO_NOT_HAVE_MUCH_MEMORY_BUT_A_LOT_OF_TIME

    // construct graph
    Teuchos::RCP<Epetra_CrsGraph> graph = Teuchos::rcp(new Epetra_CrsGraph(Copy,*rownodes,81,false));
    if (myrank==0)
    {
      for (list<std::vector<int> >::iterator i=elementnodes.begin();
           i!=elementnodes.end();
           ++i)
      {
        // get the node ids of this element
        int  numnode = static_cast<int>(i->size());
        int* nodeids = &(*i)[0];

        // loop nodes and add this topology to the row in the graph of every node
        for (int i=0; i<numnode; ++i)
        {
          int err = graph->InsertGlobalIndices(nodeids[i],numnode,nodeids);
          if (err<0) dserror("graph->InsertGlobalIndices returned %d",err);
        }
      }
    }

#else

    // No need to test for myrank==0 as elementnodes is filled on proc 0 only.

    // build the graph ourselves
    std::vector<std::set<int> > localgraph(rownodes->NumMyElements());
    for (list<std::vector<int> >::iterator i=elementnodes.begin();
         i!=elementnodes.end();
         ++i)
    {
      // get the node ids of this element
      int  numnode = static_cast<int>(i->size());
      int* nodeids = &(*i)[0];

      // loop nodes and add this topology to the row in the graph of every node
      for (int n=0; n<numnode; ++n)
      {
        int nodelid = rownodes->LID(nodeids[n]);
        copy(nodeids,
             nodeids+numnode,
             inserter(localgraph[nodelid],
                      localgraph[nodelid].begin()));
      }
    }

    elementnodes.clear();

    // fill exact entries per row vector
    // this will really speed things up for long lines
    std::vector<int> entriesperrow;
    entriesperrow.reserve(rownodes->NumMyElements());

    transform(localgraph.begin(),
              localgraph.end(),
              back_inserter(entriesperrow),
              std::mem_fun_ref(&set<int>::size));

    // construct graph
    Teuchos::RCP<Epetra_CrsGraph> graph = Teuchos::rcp(new Epetra_CrsGraph(Copy,*rownodes,&entriesperrow[0],false));

    entriesperrow.clear();

    for (unsigned i = 0; i<localgraph.size(); ++i)
    {
      std::set<int>& rowset = localgraph[i];
      std::vector<int> row;
      row.reserve(rowset.size());
      row.assign(rowset.begin(),rowset.end());
      rowset.clear();

      int err = graph->InsertGlobalIndices(rownodes->GID(i),row.size(),&row[0]);
      if (err<0) dserror("graph->InsertGlobalIndices returned %d",err);
    }

    localgraph.clear();

#endif

    elementnodes.clear();

    // finalize construction of this graph
    int err = graph->FillComplete(*rownodes,*rownodes);
    if (err) dserror("graph->FillComplete returned %d",err);

    // partition graph using metis
    Epetra_Vector weights(graph->RowMap(),false);
    weights.PutScalar(1.0);
    graph = DRT::UTILS::PartGraphUsingMetis(*graph,weights);

    // replace rownodes, colnodes with row and column maps from the graph
    // do stupid conversion from Epetra_BlockMap to Epetra_Map
    const Epetra_BlockMap& brow = graph->RowMap();
    const Epetra_BlockMap& bcol = graph->ColMap();
    rownodes = Teuchos::rcp(new Epetra_Map(brow.NumGlobalElements(),
                                   brow.NumMyElements(),
                                   brow.MyGlobalElements(),
                                   0,
                                   *comm));
    colnodes = Teuchos::rcp(new Epetra_Map(bcol.NumGlobalElements(),
                                   bcol.NumMyElements(),
                                   bcol.MyGlobalElements(),
                                   0,
                                   *comm));

    // At this point we have a good guess about the system matrix bandwidth.
    //graph->MaxNumIndices();

    graph = Teuchos::null;
}


#endif
#endif
