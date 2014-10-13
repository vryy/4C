/*!----------------------------------------------------------------------
\file drt_utils_parmetis.cpp
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

#if defined(PARALLEL) && defined(PARMETIS)

#include "drt_utils_parmetis.H"
#include "../linalg/linalg_utils.H"

#include <Epetra_Time.h>

//#ifdef HAVE_Trilinos_Q1_2013
#if 1
//Include Isorropia_Exception.hpp only because the helper functions at
//the bottom of this file (which create the epetra objects) can
//potentially throw exceptions.
#include <Isorropia_Exception.hpp>

//The Isorropia symbols being demonstrated are declared
//in these headers:
#include <Isorropia_Epetra.hpp>
#include <Isorropia_EpetraRedistributor.hpp>
#include <Isorropia_EpetraPartitioner.hpp>
#include <Isorropia_EpetraCostDescriber.hpp>
#endif

#include <iterator>

namespace DRT {
namespace UTILS {

void PackLocalConnectivity(std::map<int,std::set<int> >& lcon, DRT::PackBuffer& sblock);
void UnpackLocalConnectivity(std::map<int,std::set<int> >& lcon, std::vector<char>& rblock);

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
  std::map<int,std::set<int> > & lcon,
  DRT::PackBuffer    & sblock
  )
{
  std::map<int,std::set<int> >::iterator gidinlcon;
  std::set<int>::iterator           adjacentgid;

  // size (number of nodes we have a connectivity for)
  int size=lcon.size();
  DRT::ParObject::AddtoPack(sblock,size);

  for(gidinlcon=lcon.begin();gidinlcon!=lcon.end();++gidinlcon)
  {
    // add node gid we store the connectivity for
    DRT::ParObject::AddtoPack(sblock,(int)(gidinlcon->first));

    // add number of nodes adjacent to this one
    DRT::ParObject::AddtoPack(sblock,(int)(gidinlcon->second).size());

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
  std::map<int,std::set<int> > & lcon  ,
  std::vector<char>       & rblock
  )
{
  if(rblock.empty())
  {
    dserror("trying to extract information from an empty receive block\n");
  }

  lcon.clear();

  // position to extract
  std::vector<char>::size_type position = 0;

  // extract size (number of nodes we have a connectivity for)
  int size=0;
  DRT::ParObject::ExtractfromPack(position,rblock,size);

  for(int i=0;i<size;++i)
  {
    // extract node gid we store the connectivity for
    int gid=-1;
    DRT::ParObject::ExtractfromPack(position,rblock,gid);

    if(gid<0)
    {
      dserror("Unable to unpack a proper gid");
    }

    // extract number of adjacent nodes
    int numnb=0;
    DRT::ParObject::ExtractfromPack(position,rblock,numnb);

    if(numnb<1)
    {
      dserror("Everybody should have at least one unpackable neighbour (%d given)",numnb);
    }

    std::set<int> neighbourset;

    // extract all adjacent nodes and feed them into the set
    for(int j=0;j<numnb;++j)
    {
      int nbgid=0;
      DRT::ParObject::ExtractfromPack(position,rblock,nbgid);
      neighbourset.insert(nbgid);
    }

    // add this node connectivity to local connectivity map
    lcon.insert(std::pair<int,std::set<int> >(gid,neighbourset));
  }

  // trash receive block
  rblock.clear();

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::UTILS::PartUsingParMetis(Teuchos::RCP<DRT::Discretization> dis,
                                   Teuchos::RCP<Epetra_Map> roweles,
                                   Teuchos::RCP<Epetra_Map>& rownodes,
                                   Teuchos::RCP<Epetra_Map>& colnodes,
                                   Teuchos::RCP<Epetra_Comm> comm,
                                   bool outflag,
                                   int parts)
{
  const int myrank = comm->MyPID();
  Epetra_Time timer(dis->Comm());
  double t1 = timer.ElapsedTime();
  if (!myrank && outflag)
  {
    printf("parmetis:\n");
    fflush(stdout);
  }

  Teuchos::RCP<const Epetra_CrsGraph> graph = BuildGraph(dis,
    roweles,
    rownodes,
    comm,
    outflag);

//#ifdef HAVE_Trilinos_Q1_2013
#if 1

  // use Isorropia

  Teuchos::ParameterList paramlist;
  //No parameters. By default, Isorropia will use Zoltan hypergraph
  //partitioning, treating the graph columns as hyperedges and the
  //graph rows as vertices.

  // if the user wants to use less procs than available (as indicated by
  // the input flag "parts" above) then pass on this information to the
  // parameter list for Zoltan/Isorropia
  paramlist.set("partitioning method", "graph");
  if (parts != -1)
  {
    std::stringstream ss;
    ss << parts;
    std::string s = ss.str();
    paramlist.set("num parts", s);
  }

  //This time, we'll try graph partitioning..
  //Create a parameter sublist for Zoltan parameters.
  /*paramlist.set("PARTITIONING METHOD", "GRAPH");
  paramlist.set("STRUCTURALLY SYMMETRIC", "NO");

  //We could also call ParMetis, if Zoltan was built with ParMetis.
  Teuchos::ParameterList& sublist = paramlist.sublist("Zoltan");
  sublist.set("GRAPH_PACKAGE", "PARMETIS");
  sublist.set("PARMETIS_METHOD", "PARTKWAY");
  sublist.set("CHECK_GRAPH", "2");
  sublist.set("GRAPH_SYMMETRIZE", "TRANSPOSE");
  sublist.set("PARMETIS_OUTPUT_LEVEL", "7");*/

  //Teuchos::ParameterList& sublist = paramlist.sublist("Zoltan");
  //sublist.set("LB_APPROACH", "PARTITION");

  Epetra_CrsGraph *balanced_graph = NULL;
  try {
    balanced_graph =
      Isorropia::Epetra::createBalancedCopy(*graph, paramlist);

  }
  catch(std::exception& exc) {
    std::cout << "Isorropia::createBalancedCopy threw "
         << "exception '" << exc.what() << "' on proc "
         << myrank << std::endl;
    dserror("Error within Isorropia (graph balancing)");
  }

  // obtain the new column and row map
  Teuchos::RCP<Epetra_CrsGraph> rcp_balanced_graph = Teuchos::rcp(balanced_graph);
  rcp_balanced_graph->FillComplete();
  rcp_balanced_graph->OptimizeStorage();
  rownodes = Teuchos::rcp(new Epetra_Map(-1,
      rcp_balanced_graph->RowMap().NumMyElements(),
      rcp_balanced_graph->RowMap().MyGlobalElements(),0,*comm));
  colnodes = Teuchos::rcp(new Epetra_Map(-1,
      rcp_balanced_graph->ColMap().NumMyElements(),
      rcp_balanced_graph->ColMap().MyGlobalElements(),0,*comm));

  double t2 = timer.ElapsedTime();
  if (!myrank && outflag)
  {
    printf("Graph rebalancing:    %10.5e secs\n",t2-t1);
    fflush(stdout);
  }

#else  // use old PARMETIS partition algorithm

  // prepare parmetis input and call parmetis to partition the graph
  std::vector<int> vtxdist(numproc+1,0);
  std::vector<int> xadj(graph->RowMap().NumMyElements()+1);
  std::vector<int> adjncy(0,0);
  {
    Epetra_IntVector ridtoidx(graph->RowMap(),false);
    Epetra_IntVector cidtoidx(graph->ColMap(),false);

    xadj[0] = 0;

    // create vtxdist
    std::vector<int> svtxdist(numproc+1,0);
    svtxdist[myrank+1] = graph->RowMap().NumMyElements();
    comm->SumAll(&svtxdist[0],&vtxdist[0],numproc+1);
    for (int i=0; i<numproc; ++i) vtxdist[i+1] += vtxdist[i];

    // create a Epetra_IntVector in rowmap with new numbering starting from zero and being linear
    for (int i=0; i<ridtoidx.Map().NumMyElements(); ++i) ridtoidx[i] = i + vtxdist[myrank];
    Epetra_Import importer(graph->ColMap(),graph->RowMap());
    int err = cidtoidx.Import(ridtoidx,importer,Insert);
    if (err) dserror("Import using importer returned %d",err);

    // create xadj and adjncy;
    for (int lrid=0; lrid<graph->RowMap().NumMyElements(); ++lrid)
    {
      int metisrid = ridtoidx[lrid];
      int  numindices;
      int* indices;
      int err = graph->ExtractMyRowView(lrid,numindices,indices);
      if (err) dserror("Epetra_CrsGraph::GetMyRowView returned %d",err);
      // switch column indices to parmetis indices and add to vector of columns
      std::vector<int> metisindices(0);
      for (int i=0; i<numindices; ++i)
      {
        if (metisrid == cidtoidx[indices[i]]) continue; // must not contain main diagonal entry
        metisindices.push_back(cidtoidx[indices[i]]);
      }
      //sort(metisindices.begin(),metisindices.end());
      // beginning of new row
      xadj[lrid+1] = xadj[lrid] + (int)metisindices.size();
      for (int i=0; i<(int)metisindices.size(); ++i)
        adjncy.push_back(metisindices[i]);
    }
  }

  comm->Barrier();

  double t2 = timer.ElapsedTime();
  if (!myrank && outflag)
  {
    printf("parmetis setup    %10.5e secs\n",t2-t1);
    fflush(stdout);
  }

  std::vector<int> part;              // output
  // make sure there is memory allocated, even for empty processors
  part.reserve(graph->RowMap().NumMyElements()+1);
  part.resize(graph->RowMap().NumMyElements());
  {
    int wgtflag = 0;             // graph is not weighted
    int numflag = 0;             // numbering start from zero
    int ncon = 1;                // number of weights on each node
    int npart = numproc;         // number of partitions desired
    int options[4] = { 0,0,15,0 }; // use default metis parameters
    int edgecut = 0;             // output, number of edges cut in partitioning
    float ubvec = (float)1.05;
    std::vector<float> tpwgts(npart,static_cast<float>(1.0/(double)npart));
    MPI_Comm mpicomm=(dynamic_cast<const Epetra_MpiComm*>(&(dis->Comm())))->Comm();

    ParMETIS_V3_PartKway(&(vtxdist[0]),&(xadj[0]),&(adjncy[0]),
                         NULL,NULL,&wgtflag,&numflag,&ncon,&npart,&(tpwgts[0]),&ubvec,
                         &(options[0]),&edgecut,&(part[0]),&mpicomm);

  }

  double t3 = timer.ElapsedTime();
  if (!myrank && outflag)
  {
    printf("parmetis call     %10.5e secs\n",t3-t2);
    fflush(stdout);
  }

  // build new row map according to part
  {
    std::vector<int> mygids;
    for (int proc=0; proc<numproc; ++proc)
    {
      int size = 0;
      std::vector<int> sendpart;
      std::vector<int> sendgid;
      if (proc==myrank)
      {
        sendpart = part;
        size = (int)sendpart.size();
        sendgid.resize(size);
        for (int i=0; i<(int)sendgid.size(); ++i) sendgid[i] = graph->RowMap().GID(i);
      }
      comm->Broadcast(&size,1,proc);
      if (proc!=myrank)
      {
        sendpart.resize(size);
        sendgid.resize(size);
      }
      comm->Broadcast(&sendpart[0],size,proc);
      comm->Broadcast(&sendgid[0],size,proc);
      for (int i=0; i<(int)sendpart.size(); ++i)
      {
        if (sendpart[i] != myrank) continue;
        mygids.push_back(sendgid[i]);
      }
      comm->Barrier();
    }
    rownodes = Teuchos::rcp(new Epetra_Map(-1,(int)mygids.size(),&mygids[0],0,*comm));
  }

  // export the graph to the new row map to obtain the new column map
  {
    int bandwith = graph->GlobalMaxNumNonzeros();
    Epetra_CrsGraph finalgraph(Copy,*rownodes,bandwith,false);
    Epetra_Export exporter(graph->RowMap(),*rownodes);
    int err = finalgraph.Export(*graph,exporter,Insert);
    if (err<0) dserror("Graph export returned err=%d",err);
    graph = Teuchos::null;
    finalgraph.FillComplete();
    finalgraph.OptimizeStorage();
    colnodes = Teuchos::rcp(new Epetra_Map(-1,//finalgraph.ColMap().NumGlobalElements(),
                                  finalgraph.ColMap().NumMyElements(),
                                  finalgraph.ColMap().MyGlobalElements(),0,*comm));
  }


  double t4 = timer.ElapsedTime();
  if (!myrank && outflag)
  {
    printf("parmetis cleanup  %10.5e secs\n",t4-t3);
    fflush(stdout);
  }
#endif

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::UTILS::PartUsingZoltanWithWeights(Teuchos::RCP<DRT::Discretization> dis,
                                            Teuchos::RCP<Epetra_Map>& rownodes,
                                            Teuchos::RCP<Epetra_Map>& colnodes,
                                            bool outflag)
{
  const int myrank = dis->Comm().MyPID();
  Epetra_Time timer(dis->Comm());
  double t1 = timer.ElapsedTime();
  if (!myrank && outflag)
  {
    printf("parmetis:\n");
    fflush(stdout);
  }

  // create nodal graph of existing problem
  Teuchos::RCP<Epetra_CrsGraph> initgraph = dis->BuildNodeGraph();

  const Epetra_Map* oldnoderowmap = dis->NodeRowMap();
  // Now we're going to create a Epetra_Vector with vertex weights and a Epetra_CrsMatrix
  // for the edge weights to be used in the partitioning operation.
  // weights must be at least one for zoltan
  Teuchos::RCP<Epetra_CrsMatrix> crs_ge_weights = Teuchos::rcp(new Epetra_CrsMatrix(Copy,*oldnoderowmap, 15));
  Teuchos::RCP<Epetra_Vector> vweights = LINALG::CreateVector(*oldnoderowmap, true);

  // loop all row elements and get their cost of evaluation
  for (int i=0; i<dis->ElementRowMap()->NumMyElements(); ++i)
  {
    DRT::Element* ele = dis->lRowElement(i);
    DRT::Node** nodes = ele->Nodes();
    const int numnode = ele->NumNode();
    std::vector<int> lm(numnode);
    std::vector<int> lmrowowner(numnode);
    for(int n=0; n<numnode; ++n)
    {
      lm[n] = nodes[n]->Id();
      lmrowowner[n] = nodes[n]->Owner();
    }

    // element vector and matrix for weights of nodes and edges
    Epetra_SerialDenseMatrix edgeweigths_ele;
    Epetra_SerialDenseVector nodeweights_ele;
    // evaluate elements to get their evaluation cost
    ele->NodalConnectivity(edgeweigths_ele, nodeweights_ele);

    LINALG::Assemble(*crs_ge_weights, edgeweigths_ele, lm, lmrowowner, lm);
    LINALG::Assemble(*vweights, nodeweights_ele, lm, lmrowowner);
  }

  Teuchos::RCP<Isorropia::Epetra::CostDescriber> costs = Teuchos::rcp(new Isorropia::Epetra::CostDescriber);
  costs->setGraphEdgeWeights(crs_ge_weights);
  costs->setVertexWeights(vweights);

  Teuchos::ParameterList paramlist;
  Teuchos::ParameterList& sublist = paramlist.sublist("Zoltan");
  sublist.set("LB_APPROACH", "PARTITION");

  // Now create the partitioner object
  Teuchos::RCP<Isorropia::Epetra::Partitioner> partitioner =
    Teuchos::rcp(new Isorropia::Epetra::Partitioner(Teuchos::rcp_const_cast<const Epetra_CrsGraph>(initgraph), costs, paramlist));

  Isorropia::Epetra::Redistributor rd(partitioner);

  Teuchos::RCP<Epetra_CrsGraph> balanced_graph = rd.redistribute(*initgraph);

  // extract repartitioned maps
  const Epetra_BlockMap& rntmp = balanced_graph->RowMap();
  rownodes = Teuchos::rcp(new Epetra_Map(-1,rntmp.NumMyElements(),rntmp.MyGlobalElements(),0,dis->Comm()));
  const Epetra_BlockMap& cntmp = balanced_graph->ColMap();
  colnodes = Teuchos::rcp(new Epetra_Map(-1,cntmp.NumMyElements(),cntmp.MyGlobalElements(),0,dis->Comm()));

  double t2 = timer.ElapsedTime();
  if (!myrank && outflag)
  {
    printf("Graph rebalancing:    %10.5e secs\n",t2-t1);
    fflush(stdout);
  }

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_CrsGraph> DRT::UTILS::BuildGraph(
  Teuchos::RCP<DRT::Discretization> dis,
  Teuchos::RCP<Epetra_Map> roweles,
  Teuchos::RCP<Epetra_Map>& rownodes,
  Teuchos::RCP<Epetra_Comm> comm,
  bool outflag)
{
  const int myrank = comm->MyPID();
  const int numproc = comm->NumProc();

  // create a set of all nodes that I have
  std::set<int> mynodes;
  for (int lid=0;lid<roweles->NumMyElements();++lid)
  {
    DRT::Element* ele=dis->gElement(roweles->GID(lid));
    const int  numnode = ele->NumNode();
    const int* nodeids = ele->NodeIds();
    copy(nodeids,nodeids+numnode,inserter(mynodes,mynodes.begin()));
  }

  // build a unique row map from the overlapping sets
  for (int proc=0; proc<numproc; ++proc)
  {
    int size = 0;
    std::vector<int> recvnodes;
    if (proc==myrank)
    {
      recvnodes.clear();
      std::set<int>::iterator fool;
      for (fool = mynodes.begin(); fool != mynodes.end(); ++fool)
        recvnodes.push_back(*fool);
      size=(int)recvnodes.size();
    }
    comm->Broadcast(&size,1,proc);
    if (proc!=myrank) recvnodes.resize(size);
    comm->Broadcast(&recvnodes[0],size,proc);
    if (proc!=myrank)
    {
      for (int i=0; i<size; ++i)
      {
        std::set<int>::iterator fool = mynodes.find(recvnodes[i]);
        if (fool==mynodes.end()) continue;
        else                     mynodes.erase(fool);
      }
    }
    comm->Barrier();
  }


  // copy the set to a vector
  {
    std::vector<int> nodes;
    std::set<int>::iterator fool;
    for (fool = mynodes.begin(); fool != mynodes.end(); ++fool)
      nodes.push_back(*fool);
    mynodes.clear();
    // create a non-overlapping row map
    rownodes = Teuchos::rcp(new Epetra_Map(-1,(int)nodes.size(),&nodes[0],0,*comm));
  }


  // start building the graph object
  std::map<int,std::set<int> > locals;
  std::map<int,std::set<int> > remotes;
  for (int lid=0;lid<roweles->NumMyElements();++lid)
  {
    DRT::Element* ele=dis->gElement(roweles->GID(lid));
    const int  numnode = ele->NumNode();
    const int* nodeids = ele->NodeIds();
    for (int i=0; i<numnode; ++i)
    {
      const int lid = rownodes->LID(nodeids[i]); // am I owner of this gid?
      std::map<int,std::set<int> >* insertmap = NULL;
      if (lid != -1) insertmap = &locals;
      else           insertmap = &remotes;
      // see whether we already have an entry for nodeids[i]
      std::map<int,std::set<int> >::iterator fool = (*insertmap).find(nodeids[i]);
      if (fool==(*insertmap).end()) // no entry in that row yet
      {
        std::set<int> tmp;
        copy(nodeids,nodeids+numnode,inserter(tmp,tmp.begin()));
        (*insertmap)[nodeids[i]] = tmp;
      }
      else
      {
        std::set<int>& imap = fool->second;
        copy(nodeids,nodeids+numnode,inserter(imap,imap.begin()));
      }
    } // for (int i=0; i<numnode; ++i)
  } // for (int lid=0;lid<roweles->NumMyElements();++lid)
  // run through locals and remotes to find the max bandwith
  int maxband = 0;
  {
    int smaxband = 0;
    std::map<int,std::set<int> >::iterator fool;
    for (fool = locals.begin(); fool != locals.end(); ++fool)
      if (smaxband < (int)fool->second.size()) smaxband = (int)fool->second.size();
    for (fool = remotes.begin(); fool != remotes.end(); ++fool)
      if (smaxband < (int)fool->second.size()) smaxband = (int)fool->second.size();
    comm->MaxAll(&smaxband,&maxband,1);
  }
  if (!myrank && outflag)
  {
    printf("parmetis max nodal bandwith %d\n",maxband);
    fflush(stdout);
  }

  Teuchos::RCP<Epetra_CrsGraph> graph = Teuchos::rcp(new Epetra_CrsGraph(Copy,*rownodes,maxband,false));
  comm->Barrier();

  // fill all local entries into the graph
  {
    std::map<int,std::set<int> >::iterator fool = locals.begin();
    for (; fool != locals.end(); ++fool)
    {
      const int grid = fool->first;
      std::vector<int> cols(0,0);
      std::set<int>::iterator setfool = fool->second.begin();
      for (; setfool != fool->second.end(); ++setfool) cols.push_back(*setfool);
      int err = graph->InsertGlobalIndices(grid,(int)cols.size(),&cols[0]);
      if (err<0) dserror("Epetra_CrsGraph::InsertGlobalIndices returned %d for global row %d",err,grid);
    }
    locals.clear();
  }

  comm->Barrier();


  // now we need to communicate and add the remote entries
  for (int proc=0; proc<numproc; ++proc)
  {
    std::vector<int> recvnodes;
    int size = 0;
    if (proc==myrank)
    {
      recvnodes.clear();
      std::map<int,std::set<int> >::iterator mapfool = remotes.begin();
      for (; mapfool != remotes.end(); ++mapfool)
      {
        recvnodes.push_back((int)mapfool->second.size()+1); // length of this entry
        recvnodes.push_back(mapfool->first); // global row id
        std::set<int>::iterator fool = mapfool->second.begin();
        for (; fool!=mapfool->second.end(); ++fool) // global col ids
          recvnodes.push_back(*fool);
      }
      size = (int)recvnodes.size();
    }
    comm->Broadcast(&size,1,proc);
    if (proc!=myrank) recvnodes.resize(size);
    comm->Broadcast(&recvnodes[0],size,proc);
    if (proc!=myrank && size)
    {
      int* ptr = &recvnodes[0];
      while (ptr < &recvnodes[size-1])
      {
        int num  = *ptr;
        int grid = *(ptr+1);
        // see whether I have grid in my row map
        if (rownodes->LID(grid) != -1) // I have it, put stuff in my graph
        {
          int err = graph->InsertGlobalIndices(grid,num-1,(ptr+2));
          if (err<0) dserror("Epetra_CrsGraph::InsertGlobalIndices returned %d",err);
          ptr += (num+1);
        }
        else // I don't have it so I don't care for entries of this row, goto next row
          ptr += (num+1);
      }
    }
    comm->Barrier();
  } //  for (int proc=0; proc<numproc; ++proc)
  remotes.clear();

  comm->Barrier();

  // finish graph
  graph->FillComplete();
  graph->OptimizeStorage();

  comm->Barrier();

  return graph;
}
#endif
