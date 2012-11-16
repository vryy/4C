/*
 * aggregation_method_uncoupled.cpp
 *
 *  Created on: May 28, 2010
 *      Author: wiesner
 */

#include "aggregation_method_uncoupled.H"

#include <map>
#include <queue>
#include <Epetra_CrsGraph.h>

LINALG::AggregationMethod_Uncoupled::AggregationMethod_Uncoupled(FILE* outfile) :
  nGlobalDirichletBlocks(0),
  map_(Teuchos::null),
  amal_map_(Teuchos::null),
  nVerbose_(0)
{

}


int LINALG::AggregationMethod_Uncoupled::AmalgamateMatrix(const RCP<Epetra_CrsMatrix>& A, ParameterList& params, RCP<Epetra_CrsGraph>& amalA, RCP<Epetra_IntVector>& bdry_array)
{

  int nUnamalgamatedBlockSize = params.get("Unamalgamated BlockSize",3);          // number of real PDE equations (for setting up amalgamated maps)
  int nPDE = params.get("PDE equations",nUnamalgamatedBlockSize);                 // local block size (should be = nUnamalgamated BlockSize in most cases, used for myblockid2blocksize, necessary if there are holes in the RowMatrix of A)

  //////////////// store original map of matrix A
  map_ = Teuchos::rcp(new Epetra_Map(A->RowMap()));

  //////////////// build row map for amalgamated matrix from row map of A
  //std::map<int,int> myrowid2globalamalblockid; // local row id -> global block id
  //std::vector<int>  globalamalblockids;           // vector with global block ids for amalgamated matrix

  int cnt_amalRows = 0;
  for(int i=0; i<A->RowMap().NumMyElements(); i++)
  {
    // calculate global block id (using the GID!)
    int globalblockid = int(A->RowMap().GID(i)/nUnamalgamatedBlockSize); // change me for variable block sizes
    myrowid2globalamalblockid_[i] = globalblockid;
    if(globalamalblockid2myrowid_.count(globalblockid)>0)
    {
      globalamalblockid2myrowid_.find(globalblockid)->second.push_back(i);
    }
    else
    {
      globalamalblockid2myrowid_[globalblockid] = std::vector<int>(1,i);
    }

    if(cnt_amalRows > 0)
    {
      if(globalamalblockids_[cnt_amalRows-1]!=globalblockid)
      {
        globalamalblockids_.push_back(globalblockid);
        cnt_amalRows++;
      }
    }
    else
    {
      globalamalblockids_.push_back(globalblockid);
      cnt_amalRows++;
    }
  }

  int nvblocks = globalamalblockids_.size();   // number of local rows of amalgamated matrix

  ////////////// store information about local block sizes
  std::vector<int> myblockid2blocksize;
  for(int i=0; i<nvblocks; i++)
  {
    myblockid2blocksize.push_back(nPDE); // change me for variable block size
  }

  ////////////// generate row map for amalgamated matrix with same distribution over all procs as row map of A
  //RCP<Epetra_Map> amalA_RowMap = Teuchos::rcp(new Epetra_Map(-1,globalamalblockids_.size(),&globalamalblockids_[0],0,A->Comm()));
  amal_map_ = Teuchos::rcp(new Epetra_Map(-1,globalamalblockids_.size(),&globalamalblockids_[0],0,A->Comm()));
#ifdef DEBUG
  int nRows = A->NumMyRows();  // number of matrix rows for current proc
  if(nvblocks != amal_map_->NumMyElements()) dserror("nvblocks and NumMyElements of amalA_RowMap does not match");
  if(nVerbose_ > 5)
    cout << "PROC " << A->Comm().MyPID() << " " << globalamalblockids_.size() << " amalgamated rows, " << nRows << " normal rows " << endl;
#endif

  ////////////// extract information from overlapping column map
  std::map<int,int> globalcolid2globalamalblockid; // global col id -> local block id
  for(int i=0; i<A->ColMap().NumMyElements(); i++)
  {
    int gcolid = A->ColMap().GID(i);
    int globalblockid = int(A->ColMap().GID(i)/nUnamalgamatedBlockSize); // change me for variable block sizes
    globalcolid2globalamalblockid[gcolid] = globalblockid;
  }

#ifdef DEBUG
  if(nVerbose_ > 5)
  {
    std::map<int,int> cnt_second_vals;
    std::map<int,int>::const_iterator iter;
    for(iter=globalcolid2globalamalblockid.begin(); iter!=globalcolid2globalamalblockid.end(); ++iter)
    {
      cnt_second_vals[iter->second] = 1;
    }
    int cnt_colblocks = cnt_second_vals.size();    // number of different block nids for cols
    cout << "PROC " << A->Comm().MyPID() << " col blocks " << cnt_colblocks << " colmap elements " << A->ColMap().NumMyElements() << " domain map elments " << A->DomainMap().NumMyElements() << endl;
  }
#endif

  ////////////// find boundary blocks
  if(bdry_array != Teuchos::null) bdry_array = Teuchos::null;   // delete input object
  bdry_array = Teuchos::rcp(new Epetra_IntVector(*amal_map_,true));
  bdry_array->PutValue(AGGR_READY);

  //////////////////// create graph for amalgamated matrix
  if(amalA != Teuchos::null) amalA = Teuchos::null;   // delete input object
  amalA = Teuchos::rcp(new Epetra_CrsGraph(Copy,*amal_map_,10));

  //////////////////// fill graph for amalgamated matrix
  // loop over all local block rows
  int nBoundaryBlocks = 0;        // count local boundary blocks
  for (int i=0; i<nvblocks; i++)
  {
    std::vector<int> colblocks; // global column block ids
    bool bIsBoundaryBlock = true;

    // loop over all rows of current block
    for(int bi=0; bi<myblockid2blocksize[i]; bi++)
    {
      // extract local row information
      int row = i*myblockid2blocksize[i] + bi;   // current local row index (change me for variable block size)
      int nnz = A->NumMyEntries(row);
      int* indices;
      double* vals;
      int err = A->ExtractMyRowView(row,nnz,vals,indices);   // extract local row information
      if(err!=0) dserror("error in ExtractMyRowView");

      if(nnz==0) cout << "PROC " << A->Comm().MyPID() << " row" << " OOPS: nnz == 0???"<< endl << endl;
      int realnnz = 0;

      // loop through all local columns
      for(int k=0; k<nnz; k++)
      {
        if(A->MyLCID(indices[k])==false) dserror("problem with columns");
        int gcid = A->GCID(indices[k]);    // global column id
        if(vals[k]!=0.0) // since A is not pruned there are sometimes zeros we can avoid
        {
          colblocks.push_back(globalcolid2globalamalblockid.find(gcid)->second);
          realnnz++;    // count real nnz in matrix row
        }
      }

      // row block cannot be boundary block any more
      if(realnnz!=1) bIsBoundaryBlock = false;

    }

    // eliminate duplicate values in colblocks
    sort( colblocks.begin( ), colblocks.end( ) );
    std::vector< int >::iterator endLocation;
    endLocation = std::unique( colblocks.begin(), colblocks.end() );
    colblocks.erase(endLocation,colblocks.end());

    // colblocks now contains the global column ids of the amalgamated matrix
    // fill matrix graph
    int grid = amal_map_->GID(i);
    if(amalA->MyGRID(grid)==false) dserror("global row id doesn't belong to current proc");
    amalA->InsertGlobalIndices(grid,colblocks.size(),&colblocks[0]);

    ///////////////// mark boundary block row
    if(bIsBoundaryBlock)
    {
      (*bdry_array)[i] = AGGR_BDRY;
      nBoundaryBlocks++;
    }
  }

  amalA->FillComplete();  // it's a quadratic matrix, so DomainMap and RangeMap should be the same as amal_map_
#ifdef DEBUG
  if(amalA->DomainMap().SameAs(*amal_map_)==false) dserror("row map and domain map doesn't match");
  if(amalA->RangeMap().SameAs(*amal_map_)==false) dserror("row map and range map doesn't match");
#endif

  //////////////// store number of dirichlet rows (dirichlet blocks)
  amalA->Comm().SumAll(&nBoundaryBlocks,&nGlobalDirichletBlocks,1);

  return 0;
}

int LINALG::AggregationMethod_Uncoupled::Coarsen(const RCP<Epetra_CrsGraph>& amalA, const RCP<Epetra_IntVector>& bdry_array, ParameterList& params, RCP<Epetra_IntVector>& amal_Aggregates, int& nLocalAggregates)
{
  ///////////// prepare amal_Aggregates vector (contains local aggregate id)
  if (amal_Aggregates==Teuchos::null) amal_Aggregates = Teuchos::rcp(new Epetra_IntVector(amalA->RowMap(),true));
  amal_Aggregates->PutValue(-1);

  ///////////// create internal vector for aggregate status
  // boundary nodes are already marked as AGGR_BDRY
  RCP<Epetra_IntVector> aggr_stat = Teuchos::rcp(new Epetra_IntVector(*bdry_array));

#ifdef DEBUG
  if(nVerbose_>5)
  {
    int m = 0;  // number of boundary nodes on cur proc
    int mg = 0; // global number of boundary nodes
    for(int i=0; i<aggr_stat->MyLength(); i++)
    {
      if ((*aggr_stat)[i] == AGGR_BDRY)
        m++;
    }
    amalA->Comm().SumAll(&m,&mg,1);
    cout << "Aggregation (UC): Phase 0, " << aggr_stat->GlobalLength() - mg << " nodes + " << mg << " boundary nodes = " << aggr_stat->GlobalLength() << endl;
  }
#endif


  // return value
  int ret = -1;

  ////////////// phase 1
  ret = Phase1(amalA,params,amal_Aggregates,aggr_stat,nLocalAggregates);

  ////////////// phase 2
  if(ret > 0) // node attachement only necessary if there are nodes to be attached
  {
    if(params.get("phase 2: node attachement scheme","MaxLink") == "MaxLink")
      ret = Phase2_maxlink(amalA,params,amal_Aggregates,aggr_stat,nLocalAggregates);
    else if(params.get("phase 2: node attachement scheme","MaxLink") == "MinRank")
      ret = Phase2_minrank(amalA,params,amal_Aggregates,aggr_stat,nLocalAggregates);
  }
#ifdef DEBUG
  else if(nVerbose_ > 5) cout << "Aggregation (UC): Phase 1+, PROC " << amalA->Comm().MyPID() << " skip phase 2" << endl;
#endif

  ////////////// phase 3
  if(ret > 0)
    ret = Phase3(amalA,params,amal_Aggregates,aggr_stat,nLocalAggregates);
#ifdef DEBUG
  else if(nVerbose_ > 5) cout << "Aggregation (UC): Phase 2+, PROC " << amalA->Comm().MyPID() << " skip phase 3" << endl;
#endif

  ////////////// phase 4
  if(ret == 0) // there cannot be any non-aggregated nodes anywhere on this proc!
    Phase4(amalA,params,amal_Aggregates,aggr_stat,nLocalAggregates);  // form aggregates for boundary nodes
  else dserror("error in aggregation");

#ifdef DEBUG
  for(int inode=0; inode<amalA->NumMyRows(); inode++)
  {
    if((*aggr_stat)[inode] == AGGR_READY || // <- this is a non-aggregated node
        (*amal_Aggregates)[inode] == -1)     // <- this is not good
    {
      dserror("error in aggregation");
    }
  }
#endif


  return ret;
}


int LINALG::AggregationMethod_Uncoupled::Phase1(const RCP<Epetra_CrsGraph>& amalA, ParameterList& params, RCP<Epetra_IntVector>& amal_Aggregates, RCP<Epetra_IntVector>& aggr_stat, int& nLocalAggregates)
{
  ///////////// set some variables that might be in the parameter list later
  unsigned int min_nodes_per_aggregate = params.get("phase 1: min nodes per aggregate", 9);
  unsigned int max_nonbdry_neighbournodes = params.get("phase 1: max neighbour nodes", 2);
  int ordering = 0;       // ordering scheme (1=canonical, 2=graph)   // TODO: improve this code
  if(params.get("phase 1: node ordering scheme","graph")=="graph")
    ordering = 2;
  else if(params.get("phase 1: node ordering scheme","graph")=="canonical")
    ordering = 1;

  // some internal variables
  int aggr_count = 0;     // number of aggregates
  std::queue<int> graph_ordering_inodes;  // inodes for graph ordering
  int inode2 = 0;
  int inode = 0;
  while(inode2 < amalA->NumMyRows())
  {
    if(ordering == 1)
      inode = inode2++;
    else if(ordering == 2)
    {
      // if there are no nodes for graph ordering scheme
      if(graph_ordering_inodes.size()==0)
      {
        // add exactly one ready node for graph ordering aggregates
        for(int jnode=0; jnode<amalA->NumMyRows(); jnode++)
        {
          if((*aggr_stat)[jnode] == AGGR_READY)
          {
            graph_ordering_inodes.push(jnode);
            break;
          }
        }
      }
      if(graph_ordering_inodes.size()==0) break;  // there's no ready node any more -> end phase 1
      inode = graph_ordering_inodes.front();  // take next node
      graph_ordering_inodes.pop();          // delete this node in list
    }


    // consider further only if the current node is in READY mode
    if((*aggr_stat)[inode] == AGGR_READY)
    {
      // build testwise new aggregate
      Aggregate ag;
      ag.list.push_back(inode);

      // extract column information for current row on current proc
      //int nnz = amalA->NumMyEntries(inode);   // number of nonzero cols in current row
      //int nnz = amalA->NumMyIndices(inode);
      int nnz;
      int* indices;
      amalA->ExtractMyRowView(inode,nnz,indices);

      unsigned int cnt_neighbours = 0;             // number of potential neighbour nodes for current new aggregate

      // build tentative new aggregate
      // count the neighbors of current node that are to be aggregated
      for(int j=0; j<nnz; j++)
      {
        int mycolindex = indices[j];

        // note: this is uncoupled coarsening
        // only column indices of current proc are allowed
        if(amalA->MyLCID(mycolindex))
        {
          // check status of current neighbor node
          if((*aggr_stat)[mycolindex]==AGGR_READY || (*aggr_stat)[mycolindex]==AGGR_NOTSEL)
            ag.list.push_back(mycolindex);  // add neighbour node to current agg (temporarely)
          else if((*aggr_stat)[mycolindex]!=AGGR_BDRY)
            cnt_neighbours++;               // count bdry neighbour nodes
        } // end if, current column index belongs to this proc
      }   // end for, loop over all columns in current row

      // check if new aggregate is acceptable
      if((cnt_neighbours > max_nonbdry_neighbournodes) || // aggregate is to near to boundary nodes
          (ag.list.size() < min_nodes_per_aggregate))   // there are too less nodes in the aggregate
      {
        ///////////////// failed to build acceptable new aggregate
        // ok, the node with id inode is not good for building an aggregate around
        ag.list.clear();                     // clear tentative aggregate info
        (*aggr_stat)[inode] = AGGR_NOTSEL;   // this node isn't good for being supernode, mark him as not selected
        if(ordering == 2) // graph ordering
        {
          // even though the aggregate around inode is not perfect, we try next the nodes where
          // inode is connected to
          for(int j=0; j<nnz; j++)  // loop over all column indices
          {
            int mycolindex = indices[j];

            if((*aggr_stat)[mycolindex]==AGGR_READY)  // if node connected to inode is not aggregated by now
            {
              graph_ordering_inodes.push(mycolindex);
            }
          }
        }
      }
      else
      {
        //////////////// accept new aggregate
        ag.id = aggr_count++;   // aggregate accepted, increment aggregate counter
        (*amal_Aggregates)[inode] = ag.id;  // store aggregate number in result vector
        for(unsigned int k=0; k<ag.list.size(); k++) // loop over all elements in new aggregate
        {
          (*aggr_stat)[ag.list[k]] = AGGR_SELECTED;  // mark all nodes in new aggregate as aggregated
          (*amal_Aggregates)[ag.list[k]] = ag.id;    // store corresponding aggregate number
          if(ordering == 2)
          {
            int colnnz;
            int* colindices;
            amalA->ExtractMyRowView(ag.list[k],colnnz,colindices);  // extract neighboured nodes for k-th node of current new aggregate

            for(int kk=0; kk<colnnz; kk++) // loop over all neighbours of k-th node in current aggregate
            {
              if(aggr_stat->Map().MyLID(colindices[kk]))  // check if kk-th column index belongs to this proc
              {
                if((*aggr_stat)[colindices[kk]] == AGGR_READY)
                {
                  graph_ordering_inodes.push(colindices[kk]);
                }
              }
            }
          }
        }
      } // end if, aggregate accepted

    } // end if aggr_stat[inode] == AGGR_READY

  }   // end loop over all rows

#ifdef DEBUG
  int m = 0;  // number of non-aggregated nodes on cur proc
  int mg = 0; // global number of non-aggregated nodes
  for(int i=0; i<aggr_stat->MyLength(); i++)
  {
    if ((*aggr_stat)[i] == AGGR_READY)
      m++;
  }
  amalA->Comm().SumAll(&m,&mg,1);
  if(mg > 0)
    cout << "Aggregation (UC): Phase 1 (WARNING) " << mg << " unaggregated nodes left (status READY)" << endl;

  if(nVerbose_>5)
  {
    m = 0;  // number of aggregated nodes on cur proc
    mg = 0; // global number of aggregated nodes
    for(int i=0; i<aggr_stat->MyLength(); i++)
    {
      if ((*aggr_stat)[i] == AGGR_SELECTED)
        m++;
    }
    amalA->Comm().SumAll(&m,&mg,1);
    int gk = 0; // number of global aggs
    amalA->Comm().SumAll(&aggr_count,&gk,1);
    cout << "Aggregation (UC): Phase 1, " << mg << " nodes out of " << aggr_stat->GlobalLength() << " nodes aggregated" << endl;
    cout << "Aggregation (UC): Phase 1, PROC " << amalA->Comm().MyPID() << " " << aggr_count << " out of " << gk << " aggregates found" << endl;
  }
#endif



  //////////////// clean up
  if(graph_ordering_inodes.size() > 0)
  {
#ifdef DEBUG
    cout << "graph_ordering_nodes still contains elements" << endl;
#endif
    for(unsigned int k=0; k<graph_ordering_inodes.size(); k++)
      graph_ordering_inodes.pop();
  }

  ////////////// return number of local aggregates
  nLocalAggregates = aggr_count;

  ////////////// collect some local information
  int cntselected = 0;
  int cntbdry = 0;
  int cntnotselected = 0;
  int cntready = 0;
  for (int inode=0; inode<aggr_stat->MyLength(); inode++)
  {
    if((*aggr_stat)[inode] == AGGR_SELECTED)
      cntselected++;
    else if((*aggr_stat)[inode] == AGGR_BDRY)
      cntbdry++;
    else if((*aggr_stat)[inode] == AGGR_NOTSEL)
      cntnotselected++;
    else if((*aggr_stat)[inode] == AGGR_READY)
      cntready++;
  }
#ifdef DEBUG
  if(nVerbose_>5)
    cout << "Aggregation (UC): Phase 1, ~~~~> PROC " << amalA->Comm().MyPID() << " AGGR_READY: " << cntready << " AGGR_SELECTED: " << cntselected << " AGGR_NOTSEL: " << cntnotselected << " AGGR_BDRY: " << cntbdry << " #AGGR: " << aggr_count << endl;
#endif

  /////////// return number of nodes, that are not aggregated (and also no dirichlet boundary nodes!)
  return cntready + cntnotselected;
}

// Phase 2 MAXLINK
// search for a neighboring aggregate that has the most connections to my node
int LINALG::AggregationMethod_Uncoupled::Phase2_maxlink(const RCP<Epetra_CrsGraph>& amalA, ParameterList& params, RCP<Epetra_IntVector>& amal_Aggregates, RCP<Epetra_IntVector>& aggr_stat, int& nLocalAggregates)
{
  ///////////////// initialize local variables
#ifdef DEBUG
  int aggr_count = nLocalAggregates; // again count aggregates
#endif

  // loop over all local rows
  for(int inode=0; inode<amalA->NumMyRows(); inode++)
  {
    if((*aggr_stat)[inode] == AGGR_NOTSEL || // <- this is a non-aggregated node
        (*aggr_stat)[inode] == AGGR_READY)    // <- this is not good
    {

      int selected_aggregate = -1;

      ///////////// extract column information for current row
      int colnnz;
      int* colindices;
      amalA->ExtractMyRowView(inode,colnnz,colindices);  // extract neighboured nodes for k-th node of current new aggregate

      std::map<int,int> aggid2cntconnections;

      for(int kk=0; kk<colnnz; kk++) // loop over all neighbours
      {
        if(aggr_stat->Map().MyLID(colindices[kk]))  // check if kk-th column index belongs to this proc
        {
          if((*aggr_stat)[colindices[kk]] == AGGR_SELECTED)
          {
            // we have found a neighboured aggregate
            int aggid = (*amal_Aggregates)[colindices[kk]];   // get processor local aggregate id
            if(aggid2cntconnections.count(aggid)>0) aggid2cntconnections[aggid] = aggid2cntconnections[aggid] + 1;
            else aggid2cntconnections[aggid] = 1;
          }
        }
      }

      int maxcnt = 0;
      std::map<int,int>::const_iterator iter;
      for(iter=aggid2cntconnections.begin();iter!=aggid2cntconnections.end();++iter)
      {
        if(maxcnt < iter->second)
        {
          maxcnt = iter->second;
          selected_aggregate = iter->first;
        }
      }

      // add nodes to aggregate
      if(selected_aggregate != -1)
      {
        (*aggr_stat)[inode] = AGGR_SELECTED;
        (*amal_Aggregates)[inode] = selected_aggregate;
      }

    }
  }

#ifdef DEBUG
  int m = 0;  // number of non-aggregated nodes on cur proc
  int mg = 0; // global number of non-aggregated nodes
  for(int i=0; i<aggr_stat->MyLength(); i++)
  {
    if ((*aggr_stat)[i] == AGGR_READY)
      m++;
  }
  amalA->Comm().SumAll(&m,&mg,1);
  if(mg > 0)
    cout << "Aggregation (UC): Phase 2 (WARNING) " << mg << " unaggregated nodes left (status READY)" << endl;

  if(nVerbose_>5)
  {
    m = 0;  // number of aggregated nodes on cur proc
    mg = 0; // global number of aggregated nodes
    for(int i=0; i<aggr_stat->MyLength(); i++)
    {
      if ((*aggr_stat)[i] == AGGR_SELECTED)
        m++;
    }
    amalA->Comm().SumAll(&m,&mg,1);
    int gk = 0; // number of global aggs
    amalA->Comm().SumAll(&aggr_count,&gk,1);
    cout << "Aggregation (UC): Phase 2, " << mg << " nodes out of " << aggr_stat->GlobalLength() << " nodes aggregated" << endl;
    cout << "Aggregation (UC): Phase 2, PROC " << amalA->Comm().MyPID() << " " << aggr_count << " out of " << gk << " aggregates found" << endl;
  }
#endif


  ////////////// collect some local information
  int cntselected = 0;
  int cntbdry = 0;
  int cntnotselected = 0;
  int cntready = 0;
  for (int inode=0; inode<aggr_stat->MyLength(); inode++)
  {
    if((*aggr_stat)[inode] == AGGR_SELECTED)
      cntselected++;
    else if((*aggr_stat)[inode] == AGGR_BDRY)
      cntbdry++;
    else if((*aggr_stat)[inode] == AGGR_NOTSEL)
      cntnotselected++;
    else if((*aggr_stat)[inode] == AGGR_READY)
      cntready++;
  }
#ifdef DEBUG
  if(nVerbose_>5)
    cout << "Aggregation (UC): Phase 2, ~~~~> PROC " << amalA->Comm().MyPID() << " AGGR_READY: " << cntready << " AGGR_SELECTED: " << cntselected << " AGGR_NOTSEL: " << cntnotselected << " AGGR_BDRY: " << cntbdry << " #AGGR: " << aggr_count << endl;
#endif

  /////////// return number of nodes, that are not aggregated (and also no dirichlet boundary nodes!)
  return cntready + cntnotselected;
}

// Phase 2 MINRANK
// search for a neighboring aggregate that has the fewest number of nodes
int LINALG::AggregationMethod_Uncoupled::Phase2_minrank(const RCP<Epetra_CrsGraph>& amalA, ParameterList& params, RCP<Epetra_IntVector>& amal_Aggregates, RCP<Epetra_IntVector>& aggr_stat, int& nLocalAggregates)
{
  ///////////////// initialize local variables
#ifdef DEBUG
  int aggr_count = nLocalAggregates; // again count aggregates
#endif

  ///////////////// prepare phase 2


  // store number of nodes for each aggregate in map localaggid2aggsize
  std::map<int,int> localaggid2aggsize;
  for(int inode=0; inode<amalA->NumMyRows(); inode++)
  {
    int aggid = (*amal_Aggregates)[inode];

    if(localaggid2aggsize.count(aggid) > 0)
      localaggid2aggsize[aggid] = localaggid2aggsize.find(aggid)->second + 1;
    else
    {
      localaggid2aggsize[aggid] = 1;
    }
  }

  // loop over all local rows
  for(int inode=0; inode<amalA->NumMyRows(); inode++)
  {
    if((*aggr_stat)[inode] == AGGR_NOTSEL || // <- this is a non-aggregated node
        (*aggr_stat)[inode] == AGGR_READY)    // <- this is not good
    {
      /////////////////////////////////////////////////////////////////////////////
      // search for a neighboring aggregate that has the fewest number of nodes
      int selected_aggregate = -1;
      int mincount = 100000;

      ///////////// extract column information for current row
      int colnnz;
      int* colindices;
      amalA->ExtractMyRowView(inode,colnnz,colindices);  // extract neighboured nodes for k-th node of current new aggregate

      for(int kk=0; kk<colnnz; kk++) // loop over all neighbours
      {
        if(aggr_stat->Map().MyLID(colindices[kk]))  // check if kk-th column index belongs to this proc
        {
          if((*aggr_stat)[colindices[kk]] == AGGR_SELECTED)
          {
            // we have found a neighboured aggregate
            int m = (*amal_Aggregates)[colindices[kk]];   // get processor local aggregate id
            if(localaggid2aggsize.find(m)->second < mincount)
            {
              mincount = localaggid2aggsize.find(m)->second;
              selected_aggregate = m;
            }
          }
        }
      }

      // add nodes to aggregate
      if(selected_aggregate != -1)
      {
        (*aggr_stat)[inode] = AGGR_SELECTED;
        (*amal_Aggregates)[inode] = selected_aggregate;
      }

    }
  }

#ifdef DEBUG
  int m = 0;  // number of non-aggregated nodes on cur proc
  int mg = 0; // global number of non-aggregated nodes
  for(int i=0; i<aggr_stat->MyLength(); i++)
  {
    if ((*aggr_stat)[i] == AGGR_READY)
      m++;
  }
  amalA->Comm().SumAll(&m,&mg,1);
  if(mg > 0)
    cout << "Aggregation (UC): Phase 2 (WARNING) " << mg << " unaggregated nodes left (status READY)" << endl;

  if(nVerbose_>5)
  {
    m = 0;  // number of aggregated nodes on cur proc
    mg = 0; // global number of aggregated nodes
    for(int i=0; i<aggr_stat->MyLength(); i++)
    {
      if ((*aggr_stat)[i] == AGGR_SELECTED)
        m++;
    }
    amalA->Comm().SumAll(&m,&mg,1);
    int gk = 0; // number of global aggs
    amalA->Comm().SumAll(&aggr_count,&gk,1);
    cout << "Aggregation (UC): Phase 2, " << mg << " nodes out of " << aggr_stat->GlobalLength() << " nodes aggregated" << endl;
    cout << "Aggregation (UC): Phase 2, PROC " << amalA->Comm().MyPID() << " " << aggr_count << " out of " << gk << " aggregates found" << endl;
  }
#endif


  ////////////// collect some local information
  int cntselected = 0;
  int cntbdry = 0;
  int cntnotselected = 0;
  int cntready = 0;
  for (int inode=0; inode<aggr_stat->MyLength(); inode++)
  {
    if((*aggr_stat)[inode] == AGGR_SELECTED)
      cntselected++;
    else if((*aggr_stat)[inode] == AGGR_BDRY)
      cntbdry++;
    else if((*aggr_stat)[inode] == AGGR_NOTSEL)
      cntnotselected++;
    else if((*aggr_stat)[inode] == AGGR_READY)
      cntready++;
  }
#ifdef DEBUG
  if(nVerbose_>5)
    cout << "Aggregation (UC): Phase 2, ~~~~> PROC " << amalA->Comm().MyPID() << " AGGR_READY: " << cntready << " AGGR_SELECTED: " << cntselected << " AGGR_NOTSEL: " << cntnotselected << " AGGR_BDRY: " << cntbdry << " #AGGR: " << aggr_count << endl;
#endif

  /////////// return number of nodes, that are not aggregated (and also no dirichlet boundary nodes!)
  return cntready + cntnotselected;
}


int LINALG::AggregationMethod_Uncoupled::Phase3(const RCP<Epetra_CrsGraph>& amalA, ParameterList& params, RCP<Epetra_IntVector>& amal_Aggregates, RCP<Epetra_IntVector>& aggr_stat, int& nLocalAggregates)
{
  ///////////// form new aggregates for non-aggregated nodes
  for(int inode=0; inode<amalA->NumMyRows(); inode++)
  {
    if((*aggr_stat)[inode] == AGGR_READY || // <- this is a non-aggregated node
        (*aggr_stat)[inode] == AGGR_NOTSEL)  // <- this is not good
    {
      /////////////// build new aggregate
      Aggregate ag;
      ag.list.push_back(inode);
      ag.id = nLocalAggregates++;

      /////////////// search for un-aggregated neighbor nodes
      int colnnz;
      int* colindices;
      amalA->ExtractMyRowView(inode,colnnz,colindices);  // extract neighboured nodes for k-th node of current new aggregate

      for(int kk=0; kk<colnnz; kk++) // loop over all neighbours
      {
        if(aggr_stat->Map().MyLID(colindices[kk]))  // check if kk-th column index belongs to this proc
        {
          int colindex = colindices[kk];
          if((*aggr_stat)[colindex] != AGGR_SELECTED && (*aggr_stat)[colindex] != AGGR_BDRY)
          {
            // add node to new aggregate
            ag.list.push_back(colindex);
          }
        }
      }

      //////////////// finalize new aggregate
      for(unsigned int j=0; j<ag.list.size(); j++)
      {
        (*aggr_stat)[ag.list.at(j)] = AGGR_SELECTED;
        (*amal_Aggregates)[ag.list.at(j)] = ag.id;
      }
    }
  }

  ////////////// collect some local information
  int cntselected = 0;
  int cntbdry = 0;
  int cntnotselected = 0;
  int cntready = 0;
  for (int inode=0; inode<aggr_stat->MyLength(); inode++)
  {
    if((*aggr_stat)[inode] == AGGR_SELECTED)
      cntselected++;
    else if((*aggr_stat)[inode] == AGGR_BDRY)
      cntbdry++;
    else if((*aggr_stat)[inode] == AGGR_NOTSEL)
      cntnotselected++;
    else if((*aggr_stat)[inode] == AGGR_READY)
      cntready++;
  }
#ifdef DEBUG
  if(nVerbose_>5)
    cout << "Aggregation (UC): Phase 3, ~~~~> PROC " << amalA->Comm().MyPID() << " AGGR_READY: " << cntready << " AGGR_SELECTED: " << cntselected << " AGGR_NOTSEL: " << cntnotselected << " AGGR_BDRY: " << cntbdry << " #AGGR: " << nLocalAggregates << endl;
#endif

  /////////// return number of nodes, that are not aggregated (and also no dirichlet boundary nodes!)
  return cntready + cntnotselected;
}

int LINALG::AggregationMethod_Uncoupled::Phase4(const RCP<Epetra_CrsGraph>& amalA, ParameterList& params, RCP<Epetra_IntVector>& amal_Aggregates, RCP<Epetra_IntVector>& aggr_stat, int& nLocalAggregates)
{
  ///////////// form new aggregates for boundary nodes
  for(int inode=0; inode<amalA->NumMyRows(); inode++)
  {
    if((*aggr_stat)[inode] == AGGR_BDRY )
    {
      /////////////// build new aggregate
      Aggregate ag;
      ag.list.push_back(inode);
      ag.id = nLocalAggregates++;

      /////////////// search for un-aggregated non-boundary neighbor nodes
      int colnnz;
      int* colindices;
      amalA->ExtractMyRowView(inode,colnnz,colindices);  // extract neighboured nodes for k-th node of current new aggregate

      for(int kk=0; kk<colnnz; kk++) // loop over all neighbours
      {
        if(aggr_stat->Map().MyLID(colindices[kk]))  // check if kk-th column index belongs to this proc
        {
          int colindex = colindices[kk];
          if((*aggr_stat)[colindex] != AGGR_SELECTED && (*aggr_stat)[colindex] != AGGR_BDRY)
          {
            // add node to new aggregate
            ag.list.push_back(colindex);
          }
        }
      }

      //////////////// finalize new aggregate
      for(unsigned int j=0; j<ag.list.size(); j++)
      {
        (*aggr_stat)[ag.list.at(j)] = AGGR_SELECTED;
        (*amal_Aggregates)[ag.list.at(j)] = ag.id;
      }
    }
  }

#ifdef DEBUG
  if(nVerbose_ > 5)
  {
    int cntselected = 0;
    int cntbdry = 0;
    int cntnotselected = 0;
    int cntready = 0;
    for (int inode=0; inode<aggr_stat->MyLength(); inode++)
    {
      if((*aggr_stat)[inode] == AGGR_SELECTED)
        cntselected++;
      else if((*aggr_stat)[inode] == AGGR_BDRY)
        cntbdry++;
      else if((*aggr_stat)[inode] == AGGR_NOTSEL)
        cntnotselected++;
      else if((*aggr_stat)[inode] == AGGR_READY)
        cntready++;
    }
    cout << "Aggregation (UC): Phase 4, ~~~~> PROC " << amalA->Comm().MyPID() << " AGGR_READY: " << cntready << " AGGR_SELECTED: " << cntselected << " AGGR_NOTSEL: " << cntnotselected << " AGGR_BDRY: " << cntbdry << " #AGGR: " << nLocalAggregates << endl;
  }
#endif

  return 0;
}


RCP<Epetra_IntVector> LINALG::AggregationMethod_Uncoupled::UnamalgamateVector(const RCP<Epetra_IntVector>& amal_Aggregates)
{
  /////////////// generate result vector
  RCP<Epetra_IntVector> ret = Teuchos::rcp(new Epetra_IntVector(*map_,true));


  // loop over all amalgamated block ids for current proc
  // should be the same as iterating over all local GIDS of amal_Aggregates->Map
  std::vector<int>::const_iterator iter;
  for(iter=globalamalblockids_.begin(); iter!=globalamalblockids_.end(); ++iter)
  {
    int globalblockid = *iter;  // current block GID on this proc

    // what aggregate belongs to global block id?
    if(!amal_Aggregates->Map().MyGID(globalblockid)) cout << "WARNING: block GID not owned by this proc??" << endl;
    int aggid = (*amal_Aggregates)[amal_Aggregates->Map().LID(globalblockid)];  // determine corresponding aggregate id

    std::vector<int> rowids = globalamalblockid2myrowid_.find(globalblockid)->second; // contains all local rowids (living on map_)
    std::vector<int>::const_iterator rowid_iter;
    for(rowid_iter=rowids.begin(); rowid_iter!=rowids.end(); ++rowid_iter)
    {
      int localrowid = *rowid_iter; // current loacl row id in map_
      (*ret)[localrowid] = aggid;   // store aggregate id
    }
  }

  return ret;
}

int LINALG::AggregationMethod_Uncoupled::PruneMatrixAnisotropic(const RCP<Epetra_CrsMatrix>& A, RCP<Epetra_CrsMatrix>& prunedA, double epsilon, int curlevel)
{
  if(curlevel >= 0)
  {
    epsilon *= std::pow(0.5,curlevel);
  }

  epsilon *= epsilon;     // adapt epsilon


  //double epsilon = 0.02 * pow(0.5,curlevel); //0.03; // TODO: compiled fix in
  //epsilon = epsilon*epsilon;

#ifdef DEBUG
  if(nVerbose_ > 5)
    cout << "PROC " << A->Comm().MyPID() << " nnz of A before pruning " << A->NumMyNonzeros() << endl;
#endif

  int maxEntriesPerRow = 0;
  for (int i=0; i<A->NumMyRows(); i++)
  {
    if(A->NumMyEntries(i)>maxEntriesPerRow)
      maxEntriesPerRow = A->NumMyEntries(i);
  }


#if 0  // THIS WORKS ONLY IN PARALLEL FOR CONTINOUS GIDS STARTING WITH GID 0
  // >>>>>>>>>>>> this is only necessary for anisotropic aggregation
  vector<double> globaldiagonalentries(A->NumGlobalRows(),0.0);
  vector<double> localdiagonalentries(A->NumGlobalRows(),0.0);

  // extract matrix diagonal entries
  for(int i=0; i<A->NumMyRows(); i++)
  {
    int grid = A->GRID(i);

    int colnnz = A->NumGlobalEntries(grid);
    int numEntries;
    vector<int> colindices(colnnz,0);
    std::vector<double> colvals(colnnz,0.0);
#ifdef DEBUG
    int ret = A->ExtractGlobalRowCopy(grid,colnnz,numEntries,&colvals[0],&colindices[0]);
    if(ret !=0 )
    {
      cout << "ExtractGlobalRowCopy returned " << ret << endl;

      if(ret < 0) dserror("Error in PruneMatrix: ExtractGlobalRowCopy < 0");
    }
#else
    A->ExtractGlobalRowCopy(grid,colnnz,numEntries,&colvals[0],&colindices[0]);
#endif

    for(int j=0; j<colnnz; j++)
    {
      int globalcolid = colindices[j];

      if(grid == globalcolid)
      {
        localdiagonalentries[grid] = colvals[j]; // store diagonal entry
      }
    }
  }

  A->Comm().SumAll(&localdiagonalentries[0],&globaldiagonalentries[0],globaldiagonalentries.size());

#ifdef DEBUG
  for(unsigned int i=0; i<globaldiagonalentries.size(); i++)
  {
    if(globaldiagonalentries[i]==0.0)
      dserror("zero on diagonal?");
  }
#endif
  // <<<<<<<<<<<<<<<<<<<<
#else




  // THIS WORKS IN PARALLEL WITHOUT RESTRICTIONS
 vector<int> numrowsofproc(A->Comm().NumProc(),0);
 vector<int> localnumrows(A->Comm().NumProc(),0);
 localnumrows[A->Comm().MyPID()] = A->NumMyRows();
 A->Comm().SumAll(&localnumrows[0],&numrowsofproc[0],A->Comm().NumProc());
 int rowoffset = 0;
 for(int i=0; i<A->Comm().MyPID(); i++) rowoffset+=numrowsofproc[i];

 vector<double> globaldiagonalentries(A->NumGlobalRows(),0.0);
 vector<double> localdiagonalentries(A->NumGlobalRows(),0.0);
 vector<int> globalrowids_allproc(A->NumGlobalRows(),-1);
 vector<int> localrowids_allproc(A->NumGlobalRows(),-1);

 // extract matrix diagonal entries
 for(int i=0; i<A->NumMyRows(); i++)
 {
   int grid = A->GRID(i);
   localrowids_allproc[rowoffset + i] = grid;

   int colnnz = A->NumGlobalEntries(grid);
   int numEntries;
   vector<int> colindices(colnnz,0);
   vector<double> colvals(colnnz,0.0);
#ifdef DEBUG
   int ret = A->ExtractGlobalRowCopy(grid,colnnz,numEntries,&colvals[0],&colindices[0]);
   if(ret !=0 )
   {
     cout << "ExtractGlobalRowCopy returned " << ret << endl;

     if(ret < 0) dserror("Error in PruneMatrix: ExtractGlobalRowCopy < 0");
   }
#else
   A->ExtractGlobalRowCopy(grid,colnnz,numEntries,&colvals[0],&colindices[0]);
#endif

   for(int j=0; j<colnnz; j++)
   {
     int globalcolid = colindices[j];

     if(grid == globalcolid)
     {
       localdiagonalentries[rowoffset + i] = colvals[j]; // store diagonal entry
     }
   }
 }

 A->Comm().SumAll(&localdiagonalentries[0],&globaldiagonalentries[0],globaldiagonalentries.size());
 A->Comm().SumAll(&localrowids_allproc[0],&globalrowids_allproc[0],globalrowids_allproc.size());

#ifdef DEBUG
 for(unsigned int i=0; i<globaldiagonalentries.size(); i++)
 {
   if(globaldiagonalentries[i]==0.0)
     dserror("zero on diagonal?");
 }
 if(globaldiagonalentries.size() != globalrowids_allproc.size()) dserror("size of globalrowids_allproc and globaldiagonalentries doesn't match");
#endif

 // now all diagonal entries of all processors are stored in a locally shrinked vector globaldiagonalentries
 // create map of A->RowMap to diagonal entries
 map<int, double> grid2diagonalentry;
 for(unsigned int i=0; i<globalrowids_allproc.size(); i++)
 {
   grid2diagonalentry[globalrowids_allproc[i]] = globaldiagonalentries[i];
 }

#endif

  // define vectors for storing matrix graph information
  vector<int> globalrowids(A->NumMyRows(),-1);                // local row id -> global row id
  map<int, vector<int> > globalcolids;  // global col ids
  map<int, std::vector<double> > vals;  // vals
  vector<int> numindicesperrow(A->NumMyRows(),0);                                       // number of colindices per local row



  // extract graph information of A + pruning
  for(int i=0; i<A->NumMyRows(); i++)
  {
    int grid = A->GRID(i);    // global row id
    globalrowids[i] = grid;

    int colnnz = A->NumGlobalEntries(grid);
    int numEntries;
    vector<int> colindices(colnnz,0);
    std::vector<double> colvals(colnnz,0.0);
    int ret = A->ExtractGlobalRowCopy(grid,colnnz,numEntries,&colvals[0],&colindices[0]);

    if(ret !=0 )
    {
#ifdef DEBUG
      cout << "ExtractGlobalRowCopy returned " << ret << endl;
#endif
      if(ret < 0) dserror("Error in PruneMatrix: ExtractGlobalRowCopy < 0");
    }

    double diagonalentry = 0.0;   // local variable for diagonal entry of current row
    int diagonalcolid = -1;

    for(int j=0; j<colnnz; j++)
    {
      int globalcolid = colindices[j];

      // anisotropic aggregation
      if(grid != globalcolid)
      {
        // this is an off-diagonal entry
        double dcompare1 = colvals[j] * colvals[j];
        if(dcompare1 > 0.0)  // simple pruning
        {
          double dcompare2 = abs(grid2diagonalentry.find(grid)->second * grid2diagonalentry.find(globalcolid)->second);

          if(dcompare1 >= epsilon * dcompare2)
          {
            globalcolids[grid].push_back(globalcolid);
            vals[grid].push_back(colvals[j]);
          }
        }
      }
      else
      {
        // TODO: what about the ignored values??? -> new variant
        diagonalentry += colvals[j];
        diagonalcolid = globalcolid;
      }
    } // end for (all columns)

    // add diagonal entry
    globalcolids[grid].push_back(diagonalcolid);
    vals[grid].push_back(diagonalentry);


#ifdef DEBUG
    if(globalcolids[i].size() != vals[i].size())  dserror("number of cols in row doesn't match to number of values in same row");
    std::vector<double> tmpvals = vals[grid];
    for(unsigned int j=0; j<tmpvals.size();j++)
    {
      if(tmpvals[j] == 0.0)
        cout << "ERROR  " << endl;
    }
#endif

    // store number of colindices per row
    numindicesperrow[i] = colindices.size();
  }

  // create matrix with new graph
  if(prunedA!=Teuchos::null) prunedA = Teuchos::null;
  prunedA = Teuchos::rcp(new Epetra_CrsMatrix(Copy,A->RowMap(),&numindicesperrow[0],true));

  // fill matrix
  for (int i=0; i<prunedA->NumMyRows(); i++)
  {
    int grid = A->GRID(i);    // global row id

    int ret = prunedA->InsertGlobalValues(grid,globalcolids[grid].size(),&(vals[grid])[0],&(globalcolids[grid])[0]);
    if(ret != 0) cout << "ret != 0" << endl;
  }

  prunedA->FillComplete(true);

#ifdef DEBUG
  if(nVerbose_>5)
    cout << "PROC " << A->Comm().MyPID() << " nnz of A after pruning " << prunedA->NumMyNonzeros() << endl;
#endif
  return 0;
}


int LINALG::AggregationMethod_Uncoupled::PruneMatrixAnisotropicSimple(const RCP<Epetra_CrsMatrix>& A, RCP<Epetra_CrsMatrix>& prunedA, double epsilon, int curlevel)
{
  if(curlevel >= 0)
  {
    epsilon *= std::pow(0.5,curlevel);
  }

  epsilon *= epsilon;     // adapt epsilon

//  double epsilon = 0.03 * pow(0.7,curlevel); //0.03; // TODO: compiled fix in
//  epsilon = epsilon*epsilon;

#ifdef DEBUG
  if(nVerbose_ > 5)
    cout << "PROC " << A->Comm().MyPID() << " nnz of A before pruning " << A->NumMyNonzeros() << endl;
#endif

  int maxEntriesPerRow = 0;
  for (int i=0; i<A->NumMyRows(); i++)
  {
    if(A->NumMyEntries(i)>maxEntriesPerRow)
      maxEntriesPerRow = A->NumMyEntries(i);
  }

  // THIS WORKS IN PARALLEL WITHOUT RESTRICTIONS
 vector<int> numrowsofproc(A->Comm().NumProc(),0);
 vector<int> localnumrows(A->Comm().NumProc(),0);
 localnumrows[A->Comm().MyPID()] = A->NumMyRows();
 A->Comm().SumAll(&localnumrows[0],&numrowsofproc[0],A->Comm().NumProc());
 int rowoffset = 0;
 for(int i=0; i<A->Comm().MyPID(); i++) rowoffset+=numrowsofproc[i];


 // extract matrix diagonal entries
 RCP<Epetra_Vector> diagVec = Teuchos::rcp(new Epetra_Vector(A->RowMap(),true));
 A->ExtractDiagonalCopy(*diagVec);


  // define vectors for storing matrix graph information
  vector<int> globalrowids(A->NumMyRows(),-1);                // local row id -> global row id
  map<int, vector<int> > globalcolids;  // global col ids
  map<int, std::vector<double> > vals;  // vals
  vector<int> numindicesperrow(A->NumMyRows(),0);                                       // number of colindices per local row

  // extract graph information of A + pruning
  for(int i=0; i<A->NumMyRows(); i++)
  {
    int grid = A->GRID(i);    // global row id
    globalrowids[i] = grid;

    int colnnz = A->NumGlobalEntries(grid);
    int numEntries;
    vector<int> colindices(colnnz,0);
    std::vector<double> colvals(colnnz,0.0);
    int ret = A->ExtractGlobalRowCopy(grid,colnnz,numEntries,&colvals[0],&colindices[0]);

    if(ret !=0 )
    {
#ifdef DEBUG
      cout << "ExtractGlobalRowCopy returned " << ret << endl;
#endif
      if(ret < 0) dserror("Error in PruneMatrix: ExtractGlobalRowCopy < 0");
    }

    double diagonalentry = 0.0;   // local variable for diagonal entry of current row
    int diagonalcolid = -1;

    for(int j=0; j<colnnz; j++)
    {
      int globalcolid = colindices[j];

      // check if globalcolid belongs to current proc (so we can access corresponding diagVec entry)
      if(A->RowMap().MyGID(globalcolid) == false) continue;

      // anisotropic aggregation
      if(grid != globalcolid)
      {
        // this is an off-diagonal entry
        double dcompare1 = colvals[j] * colvals[j];
        if(dcompare1 > 0.0)  // simple pruning
        {
          double dcompare2 = abs((*diagVec)[A->RowMap().LID(grid)] * (*diagVec)[A->RowMap().LID(globalcolid)]);

          if(dcompare1 >= epsilon * dcompare2)
          {
            globalcolids[grid].push_back(globalcolid);
            vals[grid].push_back(colvals[j]);
          }
        }
      }
      else
      {
        // TODO: what about the ignored values?
        diagonalentry += colvals[j];
        diagonalcolid = globalcolid;
      }
    } // end for (all columns)

    // add diagonal entry
    globalcolids[grid].push_back(diagonalcolid);
    vals[grid].push_back(diagonalentry);


#ifdef DEBUG
    if(globalcolids[i].size() != vals[i].size())  dserror("number of cols in row doesn't match to number of values in same row");
    std::vector<double> tmpvals = vals[grid];
    for(unsigned int j=0; j<tmpvals.size();j++)
    {
      if(tmpvals[j] == 0.0)
        cout << "ERROR  " << endl;
    }
#endif

    // store number of colindices per row
    numindicesperrow[i] = colindices.size();
  }

  // create matrix with new graph
  if(prunedA!=Teuchos::null) prunedA = Teuchos::null;
  prunedA = Teuchos::rcp(new Epetra_CrsMatrix(Copy,A->RowMap(),&numindicesperrow[0],true));

  // fill matrix
  for (int i=0; i<prunedA->NumMyRows(); i++)
  {
    int grid = A->GRID(i);    // global row id

    int ret = prunedA->InsertGlobalValues(grid,globalcolids[grid].size(),&(vals[grid])[0],&(globalcolids[grid])[0]);
    if(ret != 0) cout << "ret != 0" << endl;
  }

  prunedA->FillComplete(true);

#ifdef DEBUG
  if(nVerbose_>5)
    cout << "PROC " << A->Comm().MyPID() << " nnz of A after pruning " << prunedA->NumMyNonzeros() << endl;
#endif
  return 0;
}

int LINALG::AggregationMethod_Uncoupled::PruneMatrix(const RCP<Epetra_CrsMatrix>& A, RCP<Epetra_CrsMatrix>& prunedA)
{
#ifdef DEBUG
  if(nVerbose_>5)
    cout << "PROC " << A->Comm().MyPID() << " nnz of A before pruning " << A->NumMyNonzeros() << endl;
#endif

  int maxEntriesPerRow = 0;
  for (int i=0; i<A->NumMyRows(); i++)
  {
    if(A->NumMyEntries(i)>maxEntriesPerRow)
      maxEntriesPerRow = A->NumMyEntries(i);
  }

  // define vectors for storing matrix graph information
  vector<int> globalrowids(A->NumMyRows(),-1);                // local row id -> global row id
  map<int, vector<int> > globalcolids;  // global col ids
  map<int, std::vector<double> > vals;  // vals
  vector<int> numindicesperrow(A->NumMyRows(),0);                                       // number of colindices per local row

  // extract graph information of A + pruning
  for(int i=0; i<A->NumMyRows(); i++)
  {
    int grid = A->GRID(i);    // global row id
    globalrowids[i] = grid;

    int colnnz = A->NumGlobalEntries(grid);
    int numEntries;
    vector<int> colindices(colnnz,0);
    std::vector<double> colvals(colnnz,0.0);
    int ret = A->ExtractGlobalRowCopy(grid,colnnz,numEntries,&colvals[0],&colindices[0]);

    if(ret !=0 )
    {
#ifdef DEBUG
      cout << "ExtractGlobalRowCopy returned " << ret << endl;
#endif
      if(ret < 0) dserror("Error in PruneMatrix: ExtractGlobalRowCopy < 0");
    }

    for(int j=0; j<colnnz; j++)
    {
      int globalcolid = colindices[j];

      // simple pruning
      if(colvals[j]!=0.0)   // simple pruning
      {
        globalcolids[grid].push_back(globalcolid);
        vals[grid].push_back(colvals[j]);
      }
    }


#ifdef DEBUG
    if(globalcolids[i].size() != vals[i].size())  dserror("number of cols in row doesn't match to number of values in same row");
    std::vector<double> tmpvals = vals[grid];
    for(unsigned int j=0; j<tmpvals.size();j++)
    {
      if(tmpvals[j] == 0.0)
        cout << "ERROR  " << endl;
    }
#endif

    // store number of colindices per row
    numindicesperrow[i] = colindices.size();
  }

  // create matrix with new graph
  if(prunedA!=Teuchos::null) prunedA = Teuchos::null;
  prunedA = Teuchos::rcp(new Epetra_CrsMatrix(Copy,A->RowMap(),&numindicesperrow[0],true));

  // fill matrix
  for (int i=0; i<prunedA->NumMyRows(); i++)
  {
    int grid = A->GRID(i);    // global row id

    int ret = prunedA->InsertGlobalValues(grid,globalcolids[grid].size(),&(vals[grid])[0],&(globalcolids[grid])[0]);
    if(ret != 0) cout << "ret != 0" << endl;
  }

  prunedA->FillComplete(true);

#ifdef DEBUG
  if(nVerbose_>5)
    cout << "PROC " << A->Comm().MyPID() << " nnz of A after pruning " << prunedA->NumMyNonzeros() << endl;
#endif

#if 0
  if(A->Filled()==false) dserror("A is not filled!!");
  cout << *A << endl;
  if(prunedA!=Teuchos::null) prunedA = Teuchos::null;

  prunedA = Teuchos::rcp(new Epetra_CrsMatrix(Copy,A->RangeMap(),20,false));

  // loop over all local rows of A
  for(int i=0; i<A->NumMyRows(); i++)
  {
    int grid = A->GRID(i);  // global row id

    int colnnz = A->NumGlobalEntries(grid);
    int numEntries;
    vector<int> colindices(colnnz,-1);
    std::vector<double> vals(colnnz,-1.0);
    int ret = A->ExtractGlobalRowCopy(grid,colnnz,numEntries,&vals[0],&colindices[0]);
    if(ret!=0) dserror("error extract global row view");

    std::vector<int> prunedcolids;
    std::vector<double> prunedvals;

    for(int lj=0; lj<colnnz; ++lj) // loop over all column entries
    {

        //if(A->MyLCID(colindices[lj]) && vals[lj]!=0.0)
        {
          //cout << "colindices " << colindices[lj] << endl;
          prunedcolids.push_back(colindices[lj]);
          prunedvals.push_back(vals[lj]);
        }
    } // end loop over all column entries

    if(prunedcolids.size()==0)
      dserror("zero row in pruned matrix! crash");

    ret = prunedA->InsertGlobalValues(grid,prunedcolids.size(),&prunedvals[0],&prunedcolids[0]);
    if(ret != 0)
      cout << "WARNING: InsertGlobalValues != 0 = " << ret << endl;
  }

  cout << A->DomainMap() << endl;
  if(prunedA->Filled()==true) dserror("prunedA already filled");
  prunedA->FillComplete(A->DomainMap(),A->RangeMap(),false);
  if(prunedA->Filled()==false) dserror("prunedA not filled");

  cout << *prunedA << endl;
#endif

  return 0;
}

int LINALG::AggregationMethod_Uncoupled::GetGlobalAggregates(const RCP<Epetra_CrsMatrix>& A, ParameterList& params, RCP<Epetra_IntVector>& aggrinfo, int& naggregates_local, const RCP<Epetra_MultiVector>& ThisNS)
{
  RCP<Epetra_CrsGraph> amal_graph = Teuchos::null; // graph of amalgamated A matrix (lives on amal_map_)
  RCP<Epetra_IntVector> bdry_array = Teuchos::null; // bdry array vector (lives on amal_map_)
  RCP<Epetra_IntVector> amal_aggs = Teuchos::null; // vector with aggregate ids (local aggregate ids)

  nVerbose_ = params.get("ML output",0);  // store verbosity level

  RCP<Epetra_CrsMatrix> prunedA = Teuchos::null;
  string pruningAlgorithm = params.get("aggregation method", "isotropic aggregation");
  if(pruningAlgorithm == "isotropic aggregation")
  {
    PruneMatrix(A,prunedA);
  }
  else if(pruningAlgorithm == "anisotropic aggregation")
  {
    double epsilon = params.get("anisotropic aggregation: epsilon", 0.03);
    int curlevel = params.get("Current Level",0);
    PruneMatrixAnisotropicSimple(A,prunedA,epsilon,curlevel);
  }
  else if(pruningAlgorithm == "anisotropic aggregation (coupled)")
  {
    double epsilon = params.get("anisotropic aggregation: epsilon", 0.03);
    int curlevel = params.get("Current Level",0);
    PruneMatrixAnisotropic(A,prunedA, epsilon,curlevel);
  }
  else
  {
    prunedA = A;  // no pruning at all
  }

  //AmalgamateMatrix(A,params,amal_graph,bdry_array);
  AmalgamateMatrix(prunedA,params,amal_graph,bdry_array);
  Coarsen(amal_graph,bdry_array,params,amal_aggs,naggregates_local);
  RCP<Epetra_IntVector> aggs = UnamalgamateVector(amal_aggs); // stores unamalgamted vector with local agg ids

  // transform local agg ids to global agg ids
  const Epetra_Comm& comm = A->Comm();
  vector<int> local(comm.NumProc());
  vector<int> global(comm.NumProc());
  for(int i=0;i<comm.NumProc();++i) local[i] = 0; // zero out local vector
  local[comm.MyPID()] = naggregates_local;         // store number of local aggregates in vector
  comm.SumAll(&local[0],&global[0],comm.NumProc()); // now all aggregates are known in "global" vector

  // determine offset for global agg ids on current proc
  int offset = 0;
  for(int i=0; i<comm.MyPID(); ++i) offset += global[i];

  // create aggrinfo vector
  if(aggrinfo != Teuchos::null) aggrinfo = Teuchos::null;
  aggrinfo = Teuchos::rcp(new Epetra_IntVector(*map_,true));
#ifdef DEBUG
  if(!aggs->Map().SameAs(aggrinfo->Map())) dserror("maps don't match!");
#endif

  for(int i=0; i<aggs->MyLength(); ++i)
  {
#ifdef DEBUG
    if((*aggs)[i] == -1) dserror("error in aggregation");
#endif
    (*aggrinfo)[i] = (*aggs)[i] + offset;
  }

  // sum up
  int naggregatesglobal = 0;
  for(int i=0; i<comm.NumProc(); ++i)
    naggregatesglobal += global[i];

  return naggregatesglobal;
}


