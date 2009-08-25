
#ifdef CCADISCRET

#ifdef PARMETIS

#include "drt_inputreader.H"
#include "drt_utils.H"
#include "linalg_utils.H"
#include "standardtypes_cpp.H"

#include <Epetra_Time.h>
#include <iterator>

namespace DRT {
namespace UTILS {

void PackLocalConnectivity(map<int,set<int> >& lcon, vector<char>& sblock);
void UnpackLocalConnectivity(map<int,set<int> >& lcon, vector<char>& rblock);

}
}

typedef int idxtype;
extern "C"
{
  void ParMETIS_V3_PartKway(idxtype *vtxdist, idxtype *xadj, idxtype *adjncy, idxtype *vwgt,
                            idxtype *adjwgt, int *wgtflag, int *numflag, int *ncon, int *nparts,
                            float *tpwgts, float *ubvec, int *options, int *edgecut, idxtype *part,
                            MPI_Comm *comm);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::UTILS::PackLocalConnectivity(
  map<int,set<int> > & lcon,
  vector<char>       & sblock
  )
{
  map<int,set<int> >::iterator gidinlcon;
  set<int>::iterator           adjacentgid;

  sblock.clear();
  // size (number of nodes we have a connectivity for)
  int size=lcon.size();
  DRT::ParObject::AddtoPack(sblock,size);

  for(gidinlcon=lcon.begin();gidinlcon!=lcon.end();++gidinlcon)
  {
    // add node gid we store the connectivity for
    DRT::ParObject::AddtoPack<int>(sblock,(int)(gidinlcon->first));

    // add number of nodes adjacent to this one
    DRT::ParObject::AddtoPack<int>(sblock,(int)(gidinlcon->second).size());

    // add list of neighbours to this node
    for(adjacentgid =(gidinlcon->second).begin();
        adjacentgid!=(gidinlcon->second).end();
        ++adjacentgid)
    {
      DRT::ParObject::AddtoPack(sblock,*adjacentgid);
    }
  }

  lcon.clear();

  return;
} // end Pack_lcon


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::UTILS::UnpackLocalConnectivity(
  map<int,set<int> > & lcon  ,
  vector<char>       & rblock
  )
{
  if(rblock.empty())
  {
    dserror("trying to extract information from an empty receive block\n");
  }

  lcon.clear();

  // position to extract
  int position = 0;

  // extract size (number of nodes we have a connectivity for)
  int size=0;
  DRT::ParObject::ExtractfromPack<int>(position,rblock,size);

  for(int i=0;i<size;++i)
  {
    // extract node gid we store the connectivity for
    int gid=-1;
    DRT::ParObject::ExtractfromPack<int>(position,rblock,gid);

    if(gid<0)
    {
      dserror("Unable to unpack a proper gid");
    }

    // extract number of adjacent nodes
    int numnb=0;
    DRT::ParObject::ExtractfromPack<int>(position,rblock,numnb);

    if(numnb<1)
    {
      dserror("Everybody should have at least one unpackable neighbour (%d given)",numnb);
    }

    set<int> neighbourset;

    // extract all adjacent nodes and feed them into the set
    for(int j=0;j<numnb;++j)
    {
      int nbgid=0;
      DRT::ParObject::ExtractfromPack(position,rblock,nbgid);
      neighbourset.insert(nbgid);
    }

    // add this node connectivity to local connectivity map
    lcon.insert(pair<int,set<int> >(gid,neighbourset));
  }

  // trash receive block
  rblock.clear();

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::UTILS::PartUsingParMetis(RCP<DRT::Discretization> dis,
                                   RCP<Epetra_Map> roweles,
                                   RCP<Epetra_Map>& rownodes,
                                   RCP<Epetra_Map>& colnodes,
                                   vector<int>& nids,
                                   int nblock,
                                   RCP<Epetra_Comm> comm,
                                   Epetra_Time& time,
                                   bool outflag)
{
    // construct reverse lookup from gids to vertex ids. Again,
    // this is a global, fully redundant map!
    // We need this lookup for the construction of the parmetis
    // adjacency array
    //
    //
    //                         gidtoidx
    //                      ------------->
    //                gid_i                vertexid
    //                      <-------------
    //                           nids
    //
    map<int,int> gidtoidx;
    int numnodes = static_cast<int>(nids.size());
    for (int i=0;i<numnodes;++i)
    {
      gidtoidx[nids[i]]=i;
    }

    int myrank = comm->MyPID();
    int numproc = comm->NumProc();

    if (myrank==0)
    {

      if (outflag)
      cout << time.ElapsedTime() << " secs\n";
      time.ResetStartTime();

      if (outflag)
      {
        cout << "Build initial node graph, call PARMETIS  in....";
        fflush(stdout);
      }
    }

    // --------------------------------------------------
    // - define distributed nodal row map

    // vertices and nodes --- the parmetis call will be based on a
    // consecutive numbering of vertices. We will use a lookup
    // for global node ids which will map the gid to its vertex id.
    // We need this for the construction of the adjacency array.
    //
    //         rownode gid   | (row)vertex id            vtxdist
    //                       |
    //      -----------------+----------------------------+---+
    //   +-        gid0      |      0.....................| 0 |
    //   |         gid1      |      1                     +---+
    //  p|         gid2      |      2                       ^
    //  r|         gid3      |      3                       |
    //  o|         gid4      |      4                       |
    //  c|         gid5      |      5                       |
    //  0|         gid6      |      6                       |
    //   |         gid7      |      7                       |
    //   |         gid8      |      8                       |
    //   +-        gid9      |      9                       v
    //      -----------------+----------------------------+---+
    //  .+-       gid10      |     10.....................| 10|
    //  p|        gid11      |     11                     +---+
    //  r|        gid12      |     12                       ^
    //  o|        gid13      |     13                       |
    //  c|        gid14      |     14                       |
    //  1|        gid15      |     15                       |
    //   +-       gid16      |     16                       v
    //      -----------------+----------------------------+---+
    //  p+-       gid17      |     17.....................| 17|
    //  r|        gid18      |     18                     +---+
    //  o|        gid19      |     19                       ^
    //  c|        gid20      |     20                       |
    //  2|        gid21      |     21                       |
    //   +-       gid22      |     22                       v
    //      ----------------------------------------------+---+
    //      ..............................................| 23|
    //                                                    +---+
    //

    // number of node id junks
    int nbsize = numnodes/nblock;

    // create a simple (pseudo linear) map for nodes
    int mynsize = nbsize;
    if (myrank==numproc-1)
      mynsize = numnodes-(numproc-1)*nbsize;

    // construct the initial linear node rowmap
    RCP<Epetra_Map> lin_noderowmap = rcp(new Epetra_Map(-1,mynsize,&nids[myrank*nbsize],0,*comm));

    // remember my vertex distribution for the later parmetis call
    vector<int> vtxdist(numproc+1);
    for(int np=0;np<numproc;++np)
    {
      vtxdist[np]=np*nbsize;
    }
    vtxdist[numproc]=numnodes;

    // --------------------------------------------------
    // - determine adjacency array (i.e. the infos for the node graph)
    //   using the nodal row distribution and a round robin communication
    //   of element connectivity information

    // this is a set of gids of all nodes for which we have a
    // connectivity information on this proc
    set<int> procnodes;

    // loop all eles on this proc and determine all gids for
    // which we have some connectivity information
    for(int lid=0;lid<roweles->NumMyElements();++lid)
    {
      int gid=roweles->GID(lid);

      DRT::Element* ele=dis->gElement(gid);

      // get the node ids of this element
      const int  numnode = ele->NumNode();
      const int* nodeids = ele->NodeIds();

      // all node gids of this element are inserted into a set of
      // node ids
      copy(nodeids, nodeids+numnode, inserter(procnodes,procnodes.begin()));
    }

    // ---------------------
    // build a processor local node connectivity


    //     node gid contained in one of the elements
    //                    on this proc
    //                         |
    //                         | lcon
    //                         |
    //                         v
    //    set of all gids of adjacent nodes on this proc
    //
    map<int,set<int> > lcon;

    // construct empty local map
    set<int>::iterator procnode;

    for(procnode=procnodes.begin();procnode!=procnodes.end();++procnode)
    {
      lcon.insert(pair<int, set<int> >(*procnode,set<int> ()));
    }

    // loop all eles on this proc and construct the local
    // connectivity information
    map<int,set<int> >::iterator gidinlcon;

    for(int lid=0;lid<roweles->NumMyElements();++lid)
    {
      int gid=roweles->GID(lid);

      DRT::Element* ele=dis->gElement(gid);

      // get the node ids of this element
      const int  numnode = ele->NumNode();
      const int* nodeids = ele->NodeIds();

      // loop nodeids
      for(int rr=0;rr<numnode;++rr)
      {
        // find gid of rr in map
        gidinlcon=lcon.find(nodeids[rr]);

        if(gidinlcon==lcon.end())
        {
          dserror("GID %d should be already contained in the map",nodeids[rr]);
        }

        // add nodeids to local connectivity
        copy(nodeids,
             nodeids+numnode,
             inserter((gidinlcon->second),(gidinlcon->second).begin()));
      }
    }

    //-------------------------
    // (incomplete) round robin loop to build complete node
    //  connectivity for local nodes across all processors

    //       node gid of one row node on this proc
    //                         |
    //                         | gcon
    //                         |
    //                         v
    //    set of all gids of adjacent nodes on all procs

    // prepare empty map
    map<int,set<int> > gcon;
    for(int j=0;j<lin_noderowmap->NumMyElements();++j)
    {
      int gid=lin_noderowmap->GID(j);

      gcon.insert(pair<int, set<int> >(gid,set<int> ()));
    }

    {
#ifdef PARALLEL
      // create an exporter for point to point comunication
      DRT::Exporter exporter(dis->Comm());

      // necessary variables
      MPI_Request request;

      int         tag    =-1;
      int         frompid=-1;
      int         topid  =-1;
      int         length =-1;

      // define send and receive blocks
      vector<char> sblock;
      vector<char> rblock;

#endif

      for (int np=0;np<numproc;++np)
      {
        // in the first step, we cannot receive anything
        if(np >0)
        {
#ifdef PARALLEL
          //----------------------
          // Unpack local graph from the receive block from the
          // last proc

          // make sure that you do not think you received something if
          // you didn't
          if(rblock.empty()==false)
          {
            dserror("rblock not empty");
          }

          // receive from predecessor
          frompid=(myrank+numproc-1)%numproc;
          exporter.ReceiveAny(frompid,tag,rblock,length);

          if(tag!=(myrank+numproc-1)%numproc)
          {
            dserror("received wrong message (ReceiveAny)");
          }

          exporter.Wait(request);

          // for safety
          exporter.Comm().Barrier();
#endif

          // Unpack received block
          UnpackLocalConnectivity(lcon,rblock);

        }

        // -----------------------
        // add local connectivity passing by to global
        // connectivity on this proc

        // loop this procs global connectivity
        map<int,set<int> >::iterator gidingcon;

        for(gidingcon=gcon.begin();gidingcon!=gcon.end();++gidingcon)
        {

          // search in (other) procs local connectivity
          gidinlcon=lcon.find(gidingcon->first);

          if(gidinlcon!=lcon.end())
          {
            // in this case we do have some additional
            // connectivity info from the proc owning lcon

            copy((gidinlcon->second).begin(),
                 (gidinlcon->second).end(),
                 inserter((gidingcon->second),(gidingcon->second).begin()));
          }
        }

        // in the last step, we keep everything on this proc
        if(np < numproc-1)
        {
          //-------------------------
          // Pack local graph into block to send
          PackLocalConnectivity(lcon,sblock);

#ifdef PARALLEL
          // Send block to next proc.

          tag    =myrank;
          frompid=myrank;
          topid  =(myrank+1)%numproc;

          exporter.ISend(frompid,topid,
                         &(sblock[0]),sblock.size(),
                         tag,request);
#endif
        }
      }
    }

    lcon.clear();

    // --------------------------------------------------
    // - do partitioning using parmetis

    //    gid(i)           gid(i+1)                        index gids
    //      |^                 |^
    //  nids||             nids||
    //      ||gidtoidxd        ||gidtoidx
    //      ||                 ||
    //      v|                 v|
    //   nodegid(i)       nodegid(i+1)                      node gids
    //      ^                 ^
    //      |                 |
    //      | lin_noderowmap  | lin_noderowmap
    //      |                 |
    //      v                 v
    //      i                i+1                        local equivalent indices
    //      |                 |
    //      | xadj            | xadj
    //      |                 |
    //      v                 v
    //     +-+-+-+-+-+-+-+-+-+-+                -+-+-+
    //     | | | | | | | | | | | ............... | | |      adjncy
    //     +-+-+-+-+-+-+-+-+-+-+                -+-+-+
    //
    //     |     gid(i)'s    |   gid(i+1)'s
    //     |    neighbours   |   neighbours           (numbered by global equivalent
    //                                                 indices, mapping by nids vector
    //                                                 to global ids)
    //

    // xadj points from vertex index i to the index of the
    // first adjacent vertex.
    vector<int> xadj(lin_noderowmap->NumMyElements()+1);

    // a list of adjacent nodes, adressed using xadj
    vector<int> adjncy;

    int count=0;
    xadj[0] = 0;

    for (int idx=0;idx<lin_noderowmap->NumMyElements();++idx)
    {
      // get global node id of rownode
      int  growid =lin_noderowmap->GID(idx);

      map<int,set<int> >::iterator gidingcon;
      set<int>::iterator           nbgid;

      gidingcon=gcon.find(growid);

      if(gidingcon!=gcon.end())
      {

        for(nbgid=(gidingcon->second).begin();nbgid!=(gidingcon->second).end();++nbgid)
        {

          if(*nbgid!=growid)
          {
            // reverse lookup to determine neighbours index
            int nbidx=gidtoidx[*nbgid];

            adjncy.push_back(nbidx);
            ++count;
          }
        }
      }
      else
      {
        dserror("GID %d not found in global connectivity during setup of adjncy",growid);
      }
      xadj[idx+1] = count;
    }

    // define a vector of nodeweights
    vector<int> vwgt(lin_noderowmap->NumMyElements());

    // at the moment, we use a constant weight distribution at this point
    for (int i=0; i<lin_noderowmap->NumMyElements(); ++i)
    {
      vwgt[i] = 1;
    }

    /*
      This is used to indicate if the graph is weighted.
      wgtflag can take one of four values:
      0  No weights (vwgt and adjwgt are both NULL).
      1  Weights on the edges only (vwgt is NULL).
      2  Weights on the vertices only (adjwgt is NULL).
      3  Weights on both the vertices and edges.
    */
    int wgtflag=2;
    /*
      This is used to indicate the numbering scheme that
      is used for the vtxdist, xadj, adjncy, and part
      arrays. numflag can take one of two values:
      0 C-style numbering that starts from 0.
      1 Fortran-style numbering that starts from 1.
    */
    int numflag=0;
    /*
      This is used to specify the number of weights that
      each vertex has. It is also the number of balance
      constraints that must be satisfied.
    */
    int ncon=1;
    /*
      This is used to specify the number of sub-domains
      that are desired. Note that the number of sub-domains
      is independent of the number of processors that call
      this routine.
    */
    int npart=numproc;
    /*
      This is an array of integers that is used to pass
      additional parameters for the routine. If options[0]=0,
      then the default values are used. If options[0]=1,
      then the remaining two elements of options are
      interpreted as follows:
      options[1]     This specifies the level of information
                     to be returned during the execution of
                     the algorithm. Timing information can be
                     obtained by setting this to 1. Additional
                     options for this parameter can be obtained
                     by looking at the the file defs.h in the
                     ParMETIS-Lib directory. The numerical values
                     there should be added to obtain the correct
                     value. The default value is 0.
      options[2]     This is the random number seed for the routine.
                     The default value is 15.
    */
    int options[3] = { 0,0,15 };
    /*
      Upon successful completion, the number of edges that are cut
      by the partitioning is written to this parameter.
    */
    int edgecut=0;
    /*
      This is an array of size equal to the number of locally-stored
      vertices. Upon successful completion the partition vector of
      the locally-stored vertices is written to this array.
      Note that this vector will not be redistributed by metis, it
      will just contain the information how to redistribute.
    */
    vector<int> part(lin_noderowmap->NumMyElements());
    /*
      An array of size ncon that is used to specify the imbalance
      tolerance for each vertex weight, with 1 being perfect balance
      and nparts being perfect imbalance. A value of 1.05 for each
      of the ncon weights is recommended.
    */
    float ubvec =1.05;
    /*
      An array of size ncon x nparts that is used to specify
      the fraction of vertex weight that should be distributed
      to each sub-domain for each balance constraint. If all
      of the sub-domains are to be of the same size for every
      vertex weight, then each of the ncon x nparts elements
      should be set to a value of 1/nparts. If ncon is greater
      than one, the target sub-domain weights for each sub-domain
      are stored contiguously (similar to the vwgt array). Note
      that the sum of all of the tpwgts for a give vertex weight
      should be one.
    */
    vector<float> tpwgts(npart,1.0/(double)npart);
    /*
      This is a pointer to the MPI communicator of the processes that
      call PARMETIS.
    */
    MPI_Comm mpicomm=(dynamic_cast<const Epetra_MpiComm*>(&(dis->Comm())))->Comm();

    ParMETIS_V3_PartKway(
      &(vtxdist[0]),
      &(xadj   [0]),
      &(adjncy [0]),
      &(vwgt   [0]),
      NULL         ,
      &wgtflag     ,
      &numflag     ,
      &ncon        ,
      &npart       ,
      &(tpwgts[0]) ,
      &ubvec       ,
      &(options[0]),
      &edgecut     ,
      &(part[0])   ,
      &mpicomm);

    // for each vertex j, the proc number to which this vertex
    // belongs was written to part[j]. Note that PARMETIS does
    // not redistribute the graph according to the new partitioning,
    // it simply computes the partitioning and writes it to
    // the part array.

    // construct epetra graph on linear noderowmap
    Teuchos::RCP<Epetra_CrsGraph> graph = rcp(new Epetra_CrsGraph(Copy,*lin_noderowmap,108,false));

    for (map<int,set<int> >::iterator gid=gcon.begin();
         gid!=gcon.end();
         ++gid
      )
    {
      set<int>& rowset = gid->second;
      vector<int> row;
      row.reserve(rowset.size());
      row.assign(rowset.begin(),rowset.end());
      rowset.clear();

      int err = graph->InsertGlobalIndices(gid->first,row.size(),&row[0]);
      if (err<0) dserror("graph->InsertGlobalIndices returned %d",err);
    }

    // trash the stl based graph
    gcon.clear();

    // fill graph and optimize storage
    graph->FillComplete();
    graph->OptimizeStorage();

    // redundant, global vectors to broadcast partition results
    vector<int> lglobalpart(numnodes,0);
    vector<int> globalpart (numnodes,0);

    // insert proc ids to distribute to into global vector
    for (unsigned i=0;i<part.size();++i)
    {
      lglobalpart[gidtoidx[lin_noderowmap->GID(i)]]=part[i];
    }

    // finally broadcast partition results
    comm->SumAll(&(lglobalpart[0]),&(globalpart[0]),numnodes);

    lglobalpart.clear();

    // --------------------------------------------------
    // - build final nodal row map

    // resize part array to new size
    count=0;
    for (unsigned i=0; i<globalpart.size(); ++i)
    {
      if (globalpart[i]==myrank)
      {
        ++count;
      }
    }
    part.resize(count);

    // fill it with the new distribution from the partitioning results
    count=0;
    for (unsigned i=0; i<globalpart.size(); ++i)
    {
      if (globalpart[i]==myrank)
      {
        part[count] = nids[i];
        ++count;
      }
    }

    nids.clear();

    // create map with new layout
    Epetra_Map newmap(globalpart.size(),count,&part[0],0,*comm);

    // create the output graph and export to it
    RefCountPtr<Epetra_CrsGraph> outgraph =
      rcp(new Epetra_CrsGraph(Copy,newmap,108,false));
    Epetra_Export exporter(graph->RowMap(),newmap);
    int err = outgraph->Export(*graph,exporter,Add);
    if (err<0) dserror("Graph export returned err=%d",err);

    //trash old graph
    graph=null;

    // call fill complete and optimize storage
    outgraph->FillComplete();
    outgraph->OptimizeStorage();

    // replace rownodes, colnodes with row and column maps from the graph
    // do stupid conversion from Epetra_BlockMap to Epetra_Map
    const Epetra_BlockMap& brow = outgraph->RowMap();
    const Epetra_BlockMap& bcol = outgraph->ColMap();
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
}

#endif
#endif
