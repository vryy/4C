/*----------------------------------------------------------------------*/
/*!
\file pre_exodus_validate.cpp

\brief preprocessor for exodusII format 

<pre>
Maintainer: Georg Bauer
            bauer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/bauer
            089 - 289-15252
</pre>

Validate a given BACI input file (after all preprocessing steps)

*/
/*----------------------------------------------------------------------*/
#ifdef D_EXODUS
#ifdef CCADISCRET

#include "pre_exodus_validate.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_inputreader.H"
#include <Teuchos_RefCountPtr.hpp>

#ifdef PARALLEL
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

#include "../headers/standardtypes.h"


/*!----------------------------------------------------------------------
\brief file pointers

<pre>                                                         m.gee 8/00
This structure struct _FILES allfiles is defined in input_control_global.c
and the type is in standardtypes.h
It holds all file pointers and some variables needed for the FRSYSTEM
</pre>
*----------------------------------------------------------------------*/
extern struct _FILES  allfiles; // this CCARAT struct is needed by ReadConditions()


using namespace std;
using namespace Teuchos;


void EXODUS::ValidateInputFile(string datfile)
{
  // read and check the provided header file
  cout << "checking BACI input file       --> "<<datfile<< endl;

  // do some dirty tricks in order to keep ReadConditions() running
  // (compare with ntainp_ccadiscret() )
  char* datfilename = (char*) datfile.c_str();
  allfiles.inputfile_name = datfilename;
  sprintf(allfiles.outputfile_name, "%s.err",datfile.c_str());
  if ((allfiles.out_err = fopen(allfiles.outputfile_name,"w"))==NULL)
  {
    printf("Opening of output file .err failed\n");
  }

  // communication
#ifdef PARALLEL
  int myrank = 0;
  int nproc  = 1;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  if ((nproc>1) && (myrank==0)) dserror("Using more than one processor is not supported.");
  RefCountPtr<Epetra_Comm> comm = rcp(new Epetra_MpiComm(MPI_COMM_WORLD));
#else
  RefCountPtr<Epetra_Comm> comm = rcp(new Epetra_SerialComm());
#endif

  //create a problem instance and a DatFileReader
  Teuchos::RCP<DRT::Problem> problem = DRT::Problem::Instance();
  DRT::INPUT::DatFileReader reader(datfile.c_str(), comm, 0,false);
  reader.Activate();

  // validate dynamic and solver sections
  cout<<"...Read parameters"<<endl;
  problem->ReadParameter(reader);

  // validate all condition definitions
  cout<<"...";
  problem->ReadConditions(reader);

  // materials cannot be checked via problem->ReadMaterial()
  // since the filters use a dummy definition for this method
  
  // do not read the different fields (discretizations) here,
  // since RAM might be a problem for huge problems!

  // the input file seems to be valid
  cout<<"...OK"<<endl<<endl;

  // clean up
  problem->Done();

  return;
}


#endif
#endif
