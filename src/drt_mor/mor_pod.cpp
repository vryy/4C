/*----------------------------------------------------------------------*/
/*! \file

\brief Model Order Reduction (MOR) using Proper Orthogonal Decomposition (POD)

\level 2


 *----------------------------------------------------------------------*/

#include "mor_pod.H"

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_StandardParameterEntryValidators.hpp>
#include <stdio.h>
#include <iostream>

#include "../linalg/linalg_utils_sparse_algebra_math.H"
#include "../linalg/linalg_solver.H"
#include "../linalg/linalg_mapextractor.H"
#include "../linalg/linalg_multiply.H"

#include "../drt_adapter/ad_str_structure.H"
#include "../drt_io/io.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_condition.H"
#include "../drt_io/io_control.H"

/*----------------------------------------------------------------------*
 | constructor                                            pfaller Oct17 |
 *----------------------------------------------------------------------*/
UTILS::MOR::MOR(Teuchos::RCP<DRT::Discretization> discr)
    : actdisc_(discr),
      myrank_(actdisc_->Comm().MyPID()),
      morparams_(DRT::Problem::Instance()->MORParams()),
      havemor_(false),
      projmatrix_(Teuchos::null),
      structmapr_(Teuchos::null),
      redstructmapr_(Teuchos::null),
      structrimpo_(Teuchos::null),
      structrinvimpo_(Teuchos::null)
{
  // check if model order reduction is given in input file
  if (morparams_.get<std::string>("POD_MATRIX") != std::string("none")) havemor_ = true;

  // no mor? -> exit
  if (not havemor_) return;

  // initialize temporary matrix
  Teuchos::RCP<Epetra_MultiVector> tmpmat = Teuchos::null;

  // read projection matrix from binary file
  ReadMatrix(morparams_.get<std::string>("POD_MATRIX"), tmpmat);

  // build an importer
  Teuchos::RCP<Epetra_Import> dofrowimporter =
      Teuchos::rcp(new Epetra_Import(*(actdisc_->DofRowMap()), (tmpmat->Map())));
  projmatrix_ =
      Teuchos::rcp(new Epetra_MultiVector(*(actdisc_->DofRowMap()), tmpmat->NumVectors(), true));
  int err = projmatrix_->Import(*tmpmat, *dofrowimporter, Insert, 0);
  if (err != 0) dserror("POD projection matrix could not be mapped onto the dof map");

  // check row dimension
  if (projmatrix_->GlobalLength() != actdisc_->DofRowMap()->NumGlobalElements())
    if (myrank_ == 0) dserror("Projection matrix does not match discretization.");

  // check orthogonality
  if (not IsOrthogonal(projmatrix_)) dserror("Projection matrix is not orthogonal.");

  // maps for reduced system
  structmapr_ = Teuchos::rcp(new Epetra_Map(projmatrix_->NumVectors(), 0, actdisc_->Comm()));
  redstructmapr_ = Teuchos::rcp(
      new Epetra_Map(projmatrix_->NumVectors(), projmatrix_->NumVectors(), 0, actdisc_->Comm()));
  // LINALG::AllreduceEMap cant't be used here, because NumGlobalElements will be choosen wrong

  // importers for reduced system
  structrimpo_ = Teuchos::rcp(new Epetra_Import(*structmapr_, *redstructmapr_));
  structrinvimpo_ = Teuchos::rcp(new Epetra_Import(*redstructmapr_, *structmapr_));

  return;
}

/*----------------------------------------------------------------------*
 | M_red = V^T * M * V                                    pfaller Oct17 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<LINALG::SparseMatrix> UTILS::MOR::ReduceDiagnoal(Teuchos::RCP<LINALG::SparseMatrix> M)
{
  // right multiply M * V
  Teuchos::RCP<Epetra_MultiVector> M_tmp =
      Teuchos::rcp(new Epetra_MultiVector(M->RowMap(), projmatrix_->NumVectors(), true));
  int err = M->Multiply(false, *projmatrix_, *M_tmp);
  if (err) dserror("Multiplication M * V failed.");

  // left multiply V^T * (M * V)
  Teuchos::RCP<Epetra_MultiVector> M_red_mvec =
      Teuchos::rcp(new Epetra_MultiVector(*structmapr_, M_tmp->NumVectors(), true));
  MultiplyEpetraMultiVectors(
      projmatrix_, 'T', M_tmp, 'N', redstructmapr_, structrimpo_, M_red_mvec);

  // convert Epetra_MultiVector to LINALG::SparseMatrix
  Teuchos::RCP<LINALG::SparseMatrix> M_red =
      Teuchos::rcp(new LINALG::SparseMatrix(*structmapr_, 0, false, true));
  EpetraMultiVectorToLINALGSparseMatrix(M_red_mvec, structmapr_, Teuchos::null, M_red);

  return M_red;
}

/*----------------------------------------------------------------------*
 | M_red = V^T * M                                        pfaller Oct17 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<LINALG::SparseMatrix> UTILS::MOR::ReduceOffDiagonal(
    Teuchos::RCP<LINALG::SparseMatrix> M)
{
  // right multiply M * V
  Teuchos::RCP<Epetra_MultiVector> M_tmp =
      Teuchos::rcp(new Epetra_MultiVector(M->DomainMap(), projmatrix_->NumVectors(), true));
  int err = M->Multiply(true, *projmatrix_, *M_tmp);
  if (err) dserror("Multiplication V^T * M failed.");

  // convert Epetra_MultiVector to LINALG::SparseMatrix
  Teuchos::RCP<Epetra_Map> rangemap = Teuchos::rcp(new Epetra_Map(M->DomainMap()));
  Teuchos::RCP<LINALG::SparseMatrix> M_red =
      Teuchos::rcp(new LINALG::SparseMatrix(*rangemap, 0, false, true));
  EpetraMultiVectorToLINALGSparseMatrix(M_tmp, rangemap, structmapr_, M_red);

  return M_red;
}

/*----------------------------------------------------------------------*
 | v_red = V^T * v                                        pfaller Oct17 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_MultiVector> UTILS::MOR::ReduceRHS(Teuchos::RCP<Epetra_MultiVector> v)
{
  Teuchos::RCP<Epetra_MultiVector> v_red = Teuchos::rcp(new Epetra_Vector(*structmapr_, true));
  MultiplyEpetraMultiVectors(projmatrix_, 'T', v, 'N', redstructmapr_, structrimpo_, v_red);

  return v_red;
}

/*----------------------------------------------------------------------*
 | v_red = V^T * v                                        pfaller Oct17 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> UTILS::MOR::ReduceResidual(Teuchos::RCP<Epetra_Vector> v)
{
  Teuchos::RCP<Epetra_Vector> v_tmp = Teuchos::rcp(new Epetra_Vector(*redstructmapr_));
  int err = v_tmp->Multiply('T', 'N', 1.0, *projmatrix_, *v, 0.0);
  if (err) dserror("Multiplication V^T * v failed.");

  Teuchos::RCP<Epetra_Vector> v_red = Teuchos::rcp(new Epetra_Vector(*structmapr_));
  v_red->Import(*v_tmp, *structrimpo_, Insert, 0);

  return v_red;
}

/*----------------------------------------------------------------------*
 | v = V * v_red                                          pfaller Oct17 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> UTILS::MOR::ExtendSolution(Teuchos::RCP<Epetra_Vector> v_red)
{
  Teuchos::RCP<Epetra_Vector> v_tmp = Teuchos::rcp(new Epetra_Vector(*redstructmapr_, true));
  v_tmp->Import(*v_red, *structrinvimpo_, Insert, 0);
  Teuchos::RCP<Epetra_Vector> v = Teuchos::rcp(new Epetra_Vector(*actdisc_->DofRowMap()));
  int err = v->Multiply('N', 'N', 1.0, *projmatrix_, *v_tmp, 0.0);
  if (err) dserror("Multiplication V * v_red failed.");

  return v;
}

/*----------------------------------------------------------------------*
 | multiply two Epetra MultiVectors                       pfaller Oct17 |
 *----------------------------------------------------------------------*/
void UTILS::MOR::MultiplyEpetraMultiVectors(Teuchos::RCP<Epetra_MultiVector> multivect1,
    char multivect1Trans, Teuchos::RCP<Epetra_MultiVector> multivect2, char multivect2Trans,
    Teuchos::RCP<Epetra_Map> redmap, Teuchos::RCP<Epetra_Import> impo,
    Teuchos::RCP<Epetra_MultiVector> result)
{
  // initialize temporary Epetra_MultiVector (redmap: all procs hold all elements/rows)
  Teuchos::RCP<Epetra_MultiVector> multivect_temp =
      Teuchos::rcp(new Epetra_MultiVector(*redmap, multivect2->NumVectors(), true));

  // do the multiplication: (all procs hold the full result)
  int err = multivect_temp->Multiply(
      multivect1Trans, multivect2Trans, 1.0, *multivect1, *multivect2, 0.0);
  if (err) dserror("Multiplication failed.");

  // import the result to a Epetra_MultiVector whose elements/rows are distributed over all procs
  result->Import(*multivect_temp, *impo, Insert, 0);

  return;
}

/*----------------------------------------------------------------------*
 | Epetra_MultiVector to LINALG::SparseMatrix             pfaller Oct17 |
 *----------------------------------------------------------------------*/
void UTILS::MOR::EpetraMultiVectorToLINALGSparseMatrix(Teuchos::RCP<Epetra_MultiVector> multivect,
    Teuchos::RCP<Epetra_Map> rangemap, Teuchos::RCP<Epetra_Map> domainmap,
    Teuchos::RCP<LINALG::SparseMatrix> sparsemat)
{
  // pointer to values of the Epetra_MultiVector
  double *Values;
  Values = multivect->Values();

  // loop over columns of the Epetra_MultiVector
  for (int i = 0; i < multivect->NumVectors(); i++)
  {
    // loop over rows of the Epetra_MultiVector
    for (int j = 0; j < multivect->MyLength(); j++)
    {
      // assemble the values into the LINALG::SparseMatrix (value by value)
      sparsemat->Assemble(Values[i * multivect->MyLength() + j], rangemap->GID(j), i);
    }
  }

  // Complete the LINALG::SparseMatrix
  if (domainmap == Teuchos::null)
    sparsemat->Complete();
  else
    sparsemat->Complete(*domainmap, *rangemap);
  return;
}

/*----------------------------------------------------------------------*
 | read binary projection matrix from file                pfaller Oct17 |
 |                                                                      |
 | MATLAB code to write projmatrix(DIM x dim):                          |
 |                                                                      |
 |   fid=fopen(filename, 'w');                                          |
 |   fwrite(fid, size(projmatrix, 1), 'int');                           |
 |   fwrite(fid, size(projmatrix, 2), 'int');                           |
 |   fwrite(fid, projmatrix', 'single');                                |
 |   fclose(fid);                                                       |
 |                                                                      |
 *----------------------------------------------------------------------*/
void UTILS::MOR::ReadMatrix(std::string filename, Teuchos::RCP<Epetra_MultiVector> &projmatrix)
{
  // ***************************
  // PART1: Read in Matrix Sizes
  // ***************************

  // insert path to file if necessary
  if (filename[0] != '/')
  {
    std::string pathfilename = DRT::Problem::Instance()->OutputControlFile()->InputFileName();
    std::string::size_type pos = pathfilename.rfind('/');
    if (pos != std::string::npos)
    {
      std::string path = pathfilename.substr(0, pos + 1);
      filename.insert(filename.begin(), path.begin(), path.end());
    }
  }

  // open binary file
  std::ifstream file1(
      filename.c_str(), std::ifstream::in | std::ifstream::binary | std::ifstream::ate);

  if (!file1.good()) dserror("File containing the matrix could not be opened. Check Input-File.");

  // allocation of a memory-block to read the data
  char *sizeblock = new char[8];

  // jump to beginning of file
  file1.seekg(0, std::ifstream::beg);

  // read the first 8 byte into the memory-block
  file1.read(sizeblock, std::streamoff(8));

  // close the file
  file1.close();

  // union for conversion of 4 bytes to int or float
  union CharIntFloat {
    char ValueAsChars[4];
    int ValueAsInt;
    float ValueAsFloat;
  } NumRows, NumCols;

  // number of rows
  for (int k = 0; k < 4; k++) NumRows.ValueAsChars[k] = sizeblock[k];

  // number of columns
  for (int k = 0; k < 4; k++) NumCols.ValueAsChars[k] = sizeblock[k + 4];

  delete[] sizeblock;

  // allocate multivector according to matrix size:
  Teuchos::RCP<Epetra_Map> mymap =
      Teuchos::rcp(new Epetra_Map(NumRows.ValueAsInt, 0, actdisc_->Comm()));
  projmatrix = Teuchos::rcp(new Epetra_MultiVector(*mymap, NumCols.ValueAsInt));


  // ***************************
  // PART2: Read In Matrix
  // ***************************

  // select the format we want to import matrices in
  const int formatfactor = 4;  // 8 is double 4 is single

  // union for conversion of some bytes to int or float
  union NewCharIntFloat {
    char VAsChar[formatfactor];
    double VAsDbl;
    float VAsFlt;
  } Val;

  // calculate a number of bytes that are needed to be reserved in order to fill the multivector
  int mysize = projmatrix->NumVectors() * projmatrix->MyLength() * formatfactor;

  // open binary file (again)
  std::ifstream file2(
      filename.c_str(), std::ifstream::in | std::ifstream::binary | std::ifstream::ate);

  // allocation of a memory-block to read the data
  char *memblock = new char[mysize];

  // calculation of starting points in matrix for each processor
  const Epetra_Comm &comm(actdisc_->Comm());
  const int numproc(comm.NumProc());
  const int mypid(comm.MyPID());
  std::vector<int> localnumbers(numproc, 0);
  std::vector<int> globalnumbers(numproc, 0);
  localnumbers[mypid] = mymap->NumMyElements();
  comm.SumAll(&localnumbers[0], &globalnumbers[0], numproc);

  int factor(0);
  for (int i = 0; i < mypid; i++) factor += globalnumbers[i];

  actdisc_->Comm().Barrier();

  // 64 bit number necessary, as integer can overflow for large matrices
  long long start =
      (long long)factor * (long long)projmatrix->NumVectors() * (long long)formatfactor;

  // leads to a starting point:
  file2.seekg(8 + start, std::ifstream::beg);

  // read into memory
  file2.read(memblock, std::streamoff(mysize));

  // close the file
  file2.close();

  // loop over columns and fill double array
  for (int i = 0; i < projmatrix->NumVectors(); i++)
  {
    // loop over all rows owned by the calling processor
    for (int j = 0; j < projmatrix->MyLength(); j++)
    {
      // current value
      for (int k = 0; k < formatfactor; k++)
        Val.VAsChar[k] = memblock[j * formatfactor * NumCols.ValueAsInt + i * formatfactor + k];

      // write current value to Multivector
      projmatrix->ReplaceMyValue(j, i, double(Val.VAsFlt));
    }
  }

  // delete memory-block
  delete[] memblock;

  // all procs wait until proc number 0 did finish the stuff before
  actdisc_->Comm().Barrier();

  // Inform user
  if (myrank_ == 0) std::cout << " --> Successful\n" << std::endl;

  return;
}

/*----------------------------------------------------------------------*
 | check orthogonality with M^T * M - I == 0              pfaller Oct17 |
 *----------------------------------------------------------------------*/
bool UTILS::MOR::IsOrthogonal(Teuchos::RCP<Epetra_MultiVector> M)
{
  const int n = M->NumVectors();

  // calculate V^T * V (should be identity)
  Epetra_Map map = Epetra_Map(n, n, 0, actdisc_->Comm());
  Epetra_MultiVector identity = Epetra_MultiVector(map, n, true);
  identity.Multiply('T', 'N', 1.0, *M, *M, 0.0);

  // subtract one from diagonal
  for (int i = 0; i < n; ++i) identity.SumIntoGlobalValue(i, i, -1.0);

  // inf norm of columns
  double *norms = new double[n];
  identity.NormInf(norms);

  for (int i = 0; i < n; ++i)
    if (norms[i] > 1.0e-7) return false;

  return true;
}
