/*!----------------------------------------------------------------------
\file utils.cpp
\brief A collection of helper methods for namespace DRT

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE

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
#endif

#include "drt_utils.H"
#include "drt_node.H"
#include "drt_designnode.H"
#include "drt_designelement.H"
#include "shell8.H"
#include "drt_dserror.H"


/*----------------------------------------------------------------------*
 |  create an instance of ParObject  (public)                mwgee 11/06|
 *----------------------------------------------------------------------*/
DRT::ParObject* DRT::Utils::Factory(const char* data)
{
  // mv ptr behind the size record
  const int* ptr = (const int*)data;
  ptr++;
  
  // get the type
  const int type = *ptr;
  
  switch(type)
  {
    case ParObject_Container:
    {
      DRT::Container* object = new DRT::Container();
      object->Unpack(data);
      return object;
    }
    break;
    case ParObject_Condition:
    {
      DRT::Condition* object = new DRT::Condition(DRT::Condition::condition_none);
      object->Unpack(data);
      return object;
    }
    break;
    case ParObject_Node:
    {
      double dummycoord[3] = {999.,999.,999.};
      DRT::Node* object = new DRT::Node(-1,dummycoord,-1);
      object->Unpack(data);
      return object;
    }
    break;
    case ParObject_DesignNode:
    {
      double dummycoord[3] = {999.,999.,999.};
      DRT::DesignNode* object = new DRT::DesignNode(-1,dummycoord,-1);
      object->Unpack(data);
      return object;
    }
    break;
    case ParObject_Element:
    {
      dserror("DRT::Element is pure virtual, cannot create instance");
    }
    break;
    case ParObject_DesignElement:
    {
      DRT::DesignElement* object = new DRT::DesignElement(-1,DRT::Element::element_none,-1);
      object->Unpack(data);
      return object;
    }
    break;
    case ParObject_Shell8:
    {
      DRT::Shell8* object = new DRT::Shell8(-1,-1);
      object->Unpack(data);
      return object;
    }
    break;
    default:
      dserror("Unknown type of ParObject instance: %d",type);
    break;
  }
  
  return NULL;
}

/*----------------------------------------------------------------------*
 |  partition a graph using metis  (public)                  mwgee 11/06|
 *----------------------------------------------------------------------*/
RefCountPtr<Epetra_CrsGraph> DRT::Utils::PartGraphUsingMetis(
                                             const Epetra_CrsGraph& graph,
                                             const Epetra_Vector& weights)
{
  const int myrank   = graph.Comm().MyPID();
  const int numproc  = graph.Comm().NumProc();
  
  if (numproc==1)
  {
    RefCountPtr<Epetra_CrsGraph> outgraph = rcp(new Epetra_CrsGraph(graph));
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
  vector<int> rowsend(rowmap.NumGlobalElements());
  vector<int> rowrecv(rowmap.NumGlobalElements());
  for (int i=0; i<(int)rowsend.size(); ++i) rowsend[i] = 0;  
  for (int i=0; i<rowmap.NumMyElements(); ++i) 
  {
    const int gid = rowmap.GID(i);
    rowsend[gid] = gid;
  }
  rowmap.Comm().SumAll(&rowsend[0],&rowrecv[0],rowmap.NumGlobalElements());
  rowsend.clear();
  if (myrank != workrank) rowrecv.clear();
  Epetra_Map tmap(rowmap.NumGlobalElements(),(int)rowrecv.size(),&rowrecv[0],0,rowmap.Comm());
  rowrecv.clear();
  
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
    vector<int> xadj(tmap.NumMyElements()+1);
    vector<int> adjncy(tgraph.NumGlobalNonzeros()); // the size is an upper bound
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
        adjncy[count] = gcid;
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
      part[count] = i;
      ++count;
    }
    
  // create map with new layout
  Epetra_Map newmap(size,count,&part[0],0,graph.Comm());

  // create the output graph and export to it
  RefCountPtr<Epetra_CrsGraph> outgraph = 
                           rcp(new Epetra_CrsGraph(Copy,newmap,108,false));
  Epetra_Export exporter2(graph.RowMap(),newmap);
  err = outgraph->Export(graph,exporter2,Add);
  if (err<0) dserror("Graph export returned err=%d",err);
  outgraph->FillComplete();
  outgraph->OptimizeStorage();
  
  return outgraph;
}


#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
