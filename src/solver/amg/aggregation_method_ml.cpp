/*
 * aggregation_method_ml.cpp
 *
 *  Created on: May 28, 2010
 *      Author: wiesner
 */

#include "aggregation_method_ml.H"

LINALG::AggregationMethod_ML::AggregationMethod_ML(FILE* outfile) :
  nGlobalDirichletBlocks(0)
{

}

int LINALG::AggregationMethod_ML::GetGlobalAggregates(const RCP<Epetra_CrsMatrix>& A, ParameterList& params, RCP<Epetra_IntVector>& aggrinfo, int& naggregates_local, const RCP<Epetra_MultiVector>& ThisNS)
{
  // create aggrinfo vector
  if(aggrinfo != Teuchos::null) aggrinfo = Teuchos::null;
  aggrinfo = Teuchos::rcp(new Epetra_IntVector(A->RowMap(),true));

  int naggregates = GetAggregates(A,params,*ThisNS,aggrinfo);

  const Epetra_Comm& comm = A->Comm();
  vector<int> local(comm.NumProc());
  vector<int> global(comm.NumProc());
  for (int i=0; i<comm.NumProc(); ++i) local[i] = 0;  // zero out local vector
  local[comm.MyPID()] = naggregates;                  // fill in local aggregates
  comm.SumAll(&local[0],&global[0],comm.NumProc());   // now all aggregates are in global
  int offset = 0;
  for (int i=0; i<comm.MyPID(); ++i) offset += global[i];
  for (int i=0; i<aggrinfo->MyLength(); ++i)
    if ((*aggrinfo)[i] < naggregates) (*aggrinfo)[i] += offset; // shift "local" agg id to "global" agg id
    else                              (*aggrinfo)[i] = -1;      // set agg info of all non local dofs to -1

  int naggregatesglobal = 0;
  for (int i=0; i<comm.NumProc(); ++i)    // sum up all number of aggregates over all processors
  {
    naggregatesglobal += global[i];
  }

  naggregates_local = naggregates;  // return local number of aggregates for current processor as reference
  return naggregatesglobal;
}

int LINALG::AggregationMethod_ML::GetAggregates(const RCP<Epetra_CrsMatrix>& A, ParameterList& List, const Epetra_MultiVector& ThisNS, RCP<Epetra_IntVector>& aggrinfo)
{
  if(!A->RowMap().SameAs(aggrinfo->Map())) dserror ("map of aggrinfo must match row map of operator");

  string CoarsenType    = List.get("aggregation: type","Uncoupled");
  double Threshold    = List.get("aggregation: threshold", 0.0);
  int NumPDEEquations   = List.get("PDE equations",1);
  int nsdim         = List.get("null space: dimension", -1);
  if (nsdim==-1)  cout << "dimension of null space not set" << endl;
  int size = A->RowMap().NumMyElements();

  // create ML objects
  ML_Aggregate* agg_object;
  ML_Aggregate_Create(&agg_object);
  ML_Aggregate_KeepInfo(agg_object,1);
  ML_Aggregate_Set_MaxLevels(agg_object,2);
  ML_Aggregate_Set_StartLevel(agg_object,0);
  ML_Aggregate_Set_Threshold(agg_object,Threshold);

  ML_Set_PrintLevel(List.get("ML output", 0));

  // create ML operator
  ML_Operator* ML_Ptent = 0;
  ML_Ptent = ML_Operator_Create(MLAPI::GetML_Comm());

  //if(!thisns) cout << "error: null space is NULL" << endl;
  if (ThisNS.NumVectors() == 0) dserror("error: zero-dimension null space");

  int ns_size = ThisNS.MyLength();

  double* null_vect = 0;
  ML_memory_alloc((void **)(&null_vect), sizeof(double) * ns_size * ThisNS.NumVectors(), (char*)"ns");

  int incr = 1;
  for (int v = 0 ; v < ThisNS.NumVectors() ; ++v)
    DCOPY_F77(&ns_size, (double*)ThisNS[v], &incr,
        null_vect + v * ThisNS.MyLength(), &incr);

  ML_Aggregate_Set_NullSpace(agg_object,NumPDEEquations,nsdim,null_vect,size);

  // set coarsening type
  if(CoarsenType == "Uncoupled")
    agg_object->coarsen_scheme = ML_AGGR_UNCOUPLED;
  else if (CoarsenType == "Uncoupled-MIS")
    agg_object->coarsen_scheme = ML_AGGR_HYBRIDUM;
  else if(CoarsenType == "MIS")
  { // needed for MIS, otherwise it sets the number of equations to the null space dimension
    //agg_object->max_levels = -7; // i don't understand this
    agg_object->coarsen_scheme = ML_AGGR_MIS;
  }
  else if(CoarsenType == "METIS")
    agg_object->coarsen_scheme = ML_AGGR_METIS;
  else
  {
    dserror(std::string("error: requested aggregation scheme (" + CoarsenType + ") not recognized"));
  }

  // create ML_Operator for A
  ML_Operator* ML_A = ML_Operator_Create(MLAPI::GetML_Comm());
  ML_Operator_WrapEpetraMatrix(A.get(),ML_A);

  // run coarsening process
  int NextSize = ML_Aggregate_Coarsen(agg_object, ML_A, &ML_Ptent, MLAPI::GetML_Comm());

  int* aggrmap = NULL;
  ML_Aggregate_Get_AggrMap(agg_object,0,&aggrmap);
  if (!aggrmap) dserror("aggr_info not available");

  assert (NextSize * nsdim != 0);
  for (int i=0; i<size; ++i) (*aggrinfo)[i] = aggrmap[i];

  ML_Aggregate_Destroy(&agg_object);

  ////////////////////////////////
  // -> i think, Michael forgot this
  // since we're only interested in the aggregates we can free the ML_Operators
  // now valgrind isn't complaining any more
  // but there are still two reachable blocks for Uncoupled coarsening scheme (in ml_qr_fix 15 and 20, called by ML_Aggregate_CoarsenUncoupled in line 629, ml_agg_uncoupled.c)
  // i think it is as ML_qr_fix_setNumDeadNod(numDeadNod); is never called???
  ML_Operator_Destroy(&ML_Ptent);
  ML_Operator_Destroy(&ML_A);
  ML_Ptent = NULL;
  ML_A = NULL;
  ML_qr_fix_Destroy();   // <- ok, this is missing in ML_Aggregate_CoarsenUncoupled in line 629, ml_agg_uncoupled.c

  ML_memory_free((void**)(&null_vect));  // temporary vector with null space data
  null_vect = NULL;
  ////////////////////////////////

  return (NextSize/nsdim);
}
