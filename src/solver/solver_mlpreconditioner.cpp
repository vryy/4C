/*----------------------------------------------------------------------*/
/*! \file

\brief Interface class to ML preconditoner

\brief Declaration
\level 0
\maintainer Martin Kronbichler
            http://www.lnm.mw.tum.de
            089 - 289-15235
*/

#include "../drt_lib/drt_dserror.H"

#include "ml_common.h"
#include "ml_include.h"
#include "ml_epetra_utils.h"
#include "ml_epetra.h"
#include "ml_epetra_operator.h"
#include "ml_MultiLevelPreconditioner.h"

#include "../linalg/linalg_mlapi_operator.H"  // Michael's MLAPI based ML preconditioner

#include "../linalg/linalg_utils_sparse_algebra_math.H"

#include "solver_mlpreconditioner.H"

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
LINALG::SOLVER::MLPreconditioner::MLPreconditioner(FILE* outfile, Teuchos::ParameterList& mllist)
    : PreconditionerType(outfile), mllist_(mllist)
{
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
void LINALG::SOLVER::MLPreconditioner::Setup(
    bool create, Epetra_Operator* matrix, Epetra_MultiVector* x, Epetra_MultiVector* b)
{
  SetupLinearProblem(matrix, x, b);

  if (create)
  {
    Epetra_CrsMatrix* A = dynamic_cast<Epetra_CrsMatrix*>(matrix);
    if (A == NULL) dserror("CrsMatrix expected");

    // free old matrix first
    P_ = Teuchos::null;
    Pmatrix_ = Teuchos::null;

    // create a copy of the scaled matrix
    // so we can reuse the preconditioner
    Pmatrix_ = Teuchos::rcp(new Epetra_CrsMatrix(*A));

#if 0
    // DO NOT COMMIT THIS STUFF
    LINALG::PrintMatrixInMatlabFormat("Amatrix.out",*A,true);

    std::ofstream os;

    // open file for writing
    os.open("bvector.out",std::fstream::trunc);
    os << "%%MatrixMarket matrix array real general" << std::endl;
    os << b->Map().NumGlobalElements() << " " << 1 << std::endl;

        int NumMyElements1 = b->Map().NumMyElements();
        int MaxElementSize1 = b->Map().MaxElementSize();
        int* MyGlobalElements1 = b->Map().MyGlobalElements();
        int* FirstPointInElementList1(NULL);
        if (MaxElementSize1!=1) FirstPointInElementList1 = b->Map().FirstPointInElementList();
        double ** A_Pointers = b->Pointers();

        for (int i=0; i<NumMyElements1; i++)
        {
            os << std::setw(30) << std::setprecision(16) <<  A_Pointers[0][i];    // print out values of 1. vector (only Epetra_Vector supported, no Multi_Vector)
            os << endl;
        }
        os << flush;

      // close file
      os.close();

    // write out nullspace
    os.open("nspvector.out",std::fstream::trunc);

    Teuchos::RCP<std::vector<double> > nsp = mllist_.get<Teuchos::RCP<std::vector<double> > >("nullspace");
    int nsdim = mllist_.get<int>("null space: dimension");

    os << "%%MatrixMarket matrix array real general" << std::endl;
    os << nsp->size()/nsdim << " " << nsdim << std::endl;
    for(int row = 0; row < nsp->size(); row++) {
        os << std::setw(20) << std::setprecision(16);
        os << (*nsp)[row] << std::endl;
    }

    os << flush;
    os.close();
    // END DO NOT COMMIT THIS STUFF
#endif

    mllist_.remove("init smoother", false);

    // see whether we use standard ml or our own mlapi operator
    const bool domlapioperator = mllist_.get<bool>("LINALG::AMG_Operator", false);
    if (domlapioperator)
    {
      P_ = Teuchos::rcp(new LINALG::AMG_Operator(Pmatrix_, mllist_, true));
    }
    else
    {
      P_ = Teuchos::rcp(new ML_Epetra::MultiLevelPreconditioner(*Pmatrix_, mllist_, true));

      // for debugging ML
      // dynamic_cast<ML_Epetra::MultiLevelPreconditioner&>(*P_).PrintUnused(0);
    }
  }
}
