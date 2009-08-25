
#ifdef CCADISCRET

#ifndef PARMETIS

#include "drt_inputreader.H"
#include "drt_utils.H"
#include "linalg_utils.H"
#include "standardtypes_cpp.H"

#include <Epetra_Time.h>
#include <iterator>


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::UTILS::PartUsingMetis(RCP<Epetra_Map>& rownodes,
                                RCP<Epetra_Map>& colnodes,
                                list<vector<int> >& elementnodes,
                                RCP<Epetra_Comm> comm)
{

#if 0
    if (myrank==0)
    {
      cout << "\n\nelementnodes: size=" << elementnodes.size() << endl;
      for (list<vector<int> >::iterator i=elementnodes.begin();
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
    Teuchos::RCP<Epetra_CrsGraph> graph = rcp(new Epetra_CrsGraph(Copy,*rownodes,81,false));
    if (myrank==0)
    {
      for (list<vector<int> >::iterator i=elementnodes.begin();
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
    vector<set<int> > localgraph(rownodes->NumMyElements());
    for (list<vector<int> >::iterator i=elementnodes.begin();
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
    vector<int> entriesperrow;
    entriesperrow.reserve(rownodes->NumMyElements());

    transform(localgraph.begin(),
              localgraph.end(),
              back_inserter(entriesperrow),
              mem_fun_ref(&set<int>::size));

    // construct graph
    Teuchos::RCP<Epetra_CrsGraph> graph = rcp(new Epetra_CrsGraph(Copy,*rownodes,&entriesperrow[0],false));

    entriesperrow.clear();

    for (unsigned i = 0; i<localgraph.size(); ++i)
    {
      set<int>& rowset = localgraph[i];
      vector<int> row;
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
    rownodes = rcp(new Epetra_Map(brow.NumGlobalElements(),
                                   brow.NumMyElements(),
                                   brow.MyGlobalElements(),
                                   0,
                                   *comm));
    colnodes = rcp(new Epetra_Map(bcol.NumGlobalElements(),
                                   bcol.NumMyElements(),
                                   bcol.MyGlobalElements(),
                                   0,
                                   *comm));

    // At this point we have a good guess about the system matrix bandwidth.
    //graph->MaxNumIndices();

    graph = null;
}

#endif

#endif
