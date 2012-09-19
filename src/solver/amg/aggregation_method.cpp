/*
 * aggregation_method.cpp
 *
 *  Created on: May 20, 2010
 *      Author: wiesner
 */


#include "aggregation_method.H"
#include "aggregation_method_uncoupled.H"
#include "aggregation_method_ml.H"

#include <map>
#include <queue>

LINALG::AggregationMethod::AggregationMethod(FILE* outfile) :
  nGlobalDirichletBlocks(0),
  nVerbose_(0)
{

}


int LINALG::AggregationMethod::GetGlobalAggregates(const RCP<Epetra_CrsMatrix>& A, ParameterList& params, RCP<Epetra_IntVector>& aggrinfo, int& naggregates_local, const RCP<Epetra_MultiVector>& ThisNS)
{
   return 0;
}


void LINALG::AggregationMethod::PrintIntVectorInMatlabFormat(std::string fname, const Epetra_IntVector& V, const bool newfile)
{
  int MyPID = V.Map().Comm().MyPID();
  int NumProc = V.Map().Comm().NumProc();

  std::ofstream os;

  for (int iproc=0; iproc < NumProc; iproc++)
  {
    if (MyPID==iproc)
    {
      // open file for writing
      if ((iproc == 0) && (newfile))
        os.open(fname.c_str(),std::fstream::trunc);
      else
        os.open(fname.c_str(),std::fstream::ate | std::fstream::app);

      int NumMyElements1 = V.Map().NumMyElements();
      int MaxElementSize1 = V.Map().MaxElementSize();
      int * MyGlobalElements1 = V.Map().MyGlobalElements();
      int * FirstPointInElementList1 = NULL;
      if (MaxElementSize1!=1) FirstPointInElementList1 = V.Map().FirstPointInElementList();

      if (MyPID==0)
      {
        os.width(8);
        os <<  "     MyPID"; os << "    ";
        os.width(12);
        if (MaxElementSize1==1)
          os <<  "GID  ";
        else
          os <<  "     GID/Point";
        os.width(20);
        os <<  "Value  ";
        os << endl;
      }
      for (int i=0; i < NumMyElements1; i++)
      {
        for (int ii=0; ii< V.Map().ElementSize(ii); ii++)
        {
          int iii;
          os.width(10);
          os <<  MyPID; os << "    ";
          os.width(10);
          if (MaxElementSize1==1)
          {
            os << MyGlobalElements1[i] << "    ";
            iii = i;
          }
          else
          {
            os <<  MyGlobalElements1[i]<< "/" << ii << "    ";
            iii = FirstPointInElementList1[i]+ii;
          }
          os.width(20);
          os <<  (V.Values())[iii];
          os << endl;
        }
      }
      os << flush;
    }

    // Do a few global ops to give I/O a chance to complete
    V.Map().Comm().Barrier();
    V.Map().Comm().Barrier();
    V.Map().Comm().Barrier();
  }
  return;
}


RCP<LINALG::AggregationMethod> LINALG::AggregationMethodFactory::Create(const string AggregationMethodType, FILE* outfile)
{
  try
  {
    LINALG::AggregationMethod* aggm = NULL;

    if(AggregationMethodType == "ML Uncoupled")
      aggm = new LINALG::AggregationMethod_ML(outfile);
    else if(AggregationMethodType == "Uncoupled")
      aggm = new LINALG::AggregationMethod_Uncoupled(outfile);
    else
      dserror("error");

    return rcp(aggm);
  }
  catch(string str)
  {
    cout << "Error: AggregationMethodFactory::Create: " << str << endl;
    dserror("upps");
  }
  return null;
}
