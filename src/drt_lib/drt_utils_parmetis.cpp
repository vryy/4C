/*!----------------------------------------------------------------------
\file drt_utils_parmetis.cpp
\brief Implementation A collection of helper methods for namespace DRT
<pre>
\brief Implementation
\level 0
\maintainer Martin Kronbichler
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>

*----------------------------------------------------------------------*/

#include "drt_utils_parmetis.H"
#include "../linalg/linalg_utils.H"

#include <Epetra_Time.h>

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

#include <iterator>

namespace DRT {
namespace UTILS {

void PackLocalConnectivity(std::map<int,std::set<int> >& lcon, DRT::PackBuffer& sblock);
void UnpackLocalConnectivity(std::map<int,std::set<int> >& lcon, std::vector<char>& rblock);

}
}

/*typedef int idxtype;
extern "C"
{
  void ParMETIS_V3_PartKway(idxtype *vtxdist, idxtype *xadj, idxtype *adjncy, idxtype *vwgt,
                            idxtype *adjwgt, int *wgtflag, int *numflag, int *ncon, int *nparts,
                            float *tpwgts, float *ubvec, int *options, int *edgecut, idxtype *part,
                            MPI_Comm *comm);
}*/


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

  // use Isorropia

  Teuchos::ParameterList paramlist;
  //No parameters. By default, Isorropia will use Zoltan hypergraph
  //partitioning, treating the graph columns as hyperedges and the
  //graph rows as vertices.

  // if the user wants to use less procs than available (as indicated by
  // the input flag "parts" above) then pass on this information to the
  // parameter list for Zoltan/Isorropia
  if (parts != -1)
  {
    std::stringstream ss;
    ss << parts;
    std::string s = ss.str();
    paramlist.set("num parts", s);
  }

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
  for (int proc=numproc-1; proc>=0; --proc)
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


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::UTILS::WeightedRepartitioning(
  Teuchos::RCP<DRT::Discretization> dis,
  bool assigndegreesoffreedom,
  bool initelements,
  bool doboundaryconditions)
{
  // maps to be filled with final distributed node maps
  Teuchos::RCP<Epetra_Map> rownodes;
  Teuchos::RCP<Epetra_Map> colnodes;
  // do weighted repartitioning
  DRT::UTILS::PartUsingZoltanWithWeights(dis, rownodes, colnodes, true);

  // rebuild of the system with new maps
  Teuchos::RCP<Epetra_Map> roweles;
  Teuchos::RCP<Epetra_Map> coleles;
  dis->BuildElementRowColumn(*rownodes,*colnodes,roweles,coleles);
  dis->ExportRowNodes(*rownodes);
  dis->ExportRowElements(*roweles);
  dis->ExportColumnNodes(*colnodes);
  dis->ExportColumnElements(*coleles);

  // do final fill complete
  dis->FillComplete(assigndegreesoffreedom,initelements,doboundaryconditions);

  return;
}
