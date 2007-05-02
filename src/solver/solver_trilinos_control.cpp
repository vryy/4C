/*!----------------------------------------------------------------------
\file
\brief

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/
#ifdef TRILINOS_PACKAGE

#include <ctime>
#include <cstdlib>
#include <iostream>

#ifdef PARALLEL
#include <mpi.h>
#endif

#ifdef PARALLEL
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_CrsMatrix.h"
#include "Epetra_LinearProblem.h"

#include "Amesos_Klu.h"
#include "Amesos_Umfpack.h"
#include "Amesos_Lapack.h"

// Trilinos is configured SuperLUDIST only in the parallel version
#ifdef PARALLEL
#include "Amesos_Superludist.h"
#endif

#include "AztecOO.h"

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RefCountPtr.hpp"

#include "Ifpack.h"
#include "Ifpack_AdditiveSchwarz.h"

#include "ml_common.h"
#include "ml_include.h"
#include "ml_epetra_utils.h"
#include "ml_epetra.h"
#include "ml_epetra_operator.h"
#include "ml_MultiLevelPreconditioner.h"

#include "../solver/solver_trilinos_control.H"
#include "../solver/solver_trilinos_spooles.H"
#include "../solver/solver_trilinos_ml.H"

using namespace std;
using namespace Teuchos;

/*!----------------------------------------------------------------------
\brief file pointers

<pre>                                                         m.gee 8/00
This structure struct _FILES allfiles is defined in input_control_global.c
and the type is in standardtypes.h
It holds all file pointers and some variables needed for the FRSYSTEM
</pre>
*----------------------------------------------------------------------*/
extern struct _FILES  allfiles;

// local methods
static void solve_amesos_klu(TRILINOSMATRIX* tri,
                             DIST_VECTOR*    sol,
                             DIST_VECTOR*    rhs,
                             int             option,
                             bool            sym);
static void solve_amesos_umfpack(TRILINOSMATRIX* tri,
                                 DIST_VECTOR*    sol,
                                 DIST_VECTOR*    rhs,
                                 int             option);
#ifdef PARALLEL
static void solve_amesos_superlu(TRILINOSMATRIX* tri,
                                 DIST_VECTOR*    sol,
                                 DIST_VECTOR*    rhs,
                                 int             option);
#endif
static void solve_amesos_lapack(TRILINOSMATRIX* tri,
                                DIST_VECTOR*    sol,
                                DIST_VECTOR*    rhs,
                                int             option,
                                bool            sym);
static void solve_aztecoo(TRILINOSMATRIX* tri,
                          DIST_VECTOR*    sol,
                          DIST_VECTOR*    rhs,
                          int             option,
                          SOLVAR*         actsolv,
                          FIELD*          actfield,
                          int             disnum);

/*----------------------------------------------------------------------*
  |                                                          m.gee 9/06 |
  |  routine to control all solver calls                                |
  |  in the case where Trilinos matrices are used internally            |
 *----------------------------------------------------------------------*/
void solver_trilinos_control(struct _FIELD          *actfield,
                             int                     disnum,
                             struct _SOLVAR         *actsolv,
                             struct _INTRA          *actintra,
                             enum   _SPARSE_TYP     *sysarray_typ,
                             union  _SPARSE_ARRAY   *sysarray,
                             struct _DIST_VECTOR    *sol,
                             struct _DIST_VECTOR    *rhs,
                             INT                     option)
{
#ifdef DEBUG
  dstrc_enter("solver_trilinos_control");
#endif

// check correct type of matrix
if (*sysarray_typ != trilinos)
  dserror("Matrix is not in Trilinos format");

// get ptr to matrix and epetra object
TRILINOSMATRIX* tri = sysarray->trilinos;

// in init phase, make sure we are clean
if (option==1)
{
  tri->is_factored=0;
  tri->is_init=1;
}

switch (actsolv->solvertyp)
{
  //---------------------------------------------------------------------------
  case superlu:
#ifdef PARALLEL
    solve_amesos_superlu(tri,sol,rhs,option);
#else
    dserror("Superludist only with -DPARALLEL");
#endif
  break;

  //---------------------------------------------------------------------------
  case amesos_klu_sym:
    solve_amesos_klu(tri,sol,rhs,option,true);
  break;

  //---------------------------------------------------------------------------
  case amesos_klu_nonsym:
    solve_amesos_klu(tri,sol,rhs,option,false);
  break;

  //---------------------------------------------------------------------------
  case umfpack:
    solve_amesos_umfpack(tri,sol,rhs,option);
  break;

  //---------------------------------------------------------------------------
  case lapack_sym:
    solve_amesos_lapack(tri,sol,rhs,option,true);
  break;

  //---------------------------------------------------------------------------
  case lapack_nonsym:
    solve_amesos_lapack(tri,sol,rhs,option,false);
  break;

  //---------------------------------------------------------------------------
  case aztec_msr:
    solve_aztecoo(tri,sol,rhs,option,actsolv,actfield,disnum);
  break;

  //---------------------------------------------------------------------------
#ifdef SPOOLES_PACKAGE
  case SPOOLES_nonsym:
    if (option==1) // init phase
    {
      tri->sysarray_typ = (SPARSE_TYP*)  CCACALLOC(1,sizeof(SPARSE_TYP));
      tri->sysarray     = (SPARSE_ARRAY*)CCACALLOC(1,sizeof(SPARSE_ARRAY));
      tri->sysarray_typ[0] = spoolmatrix;
      tri->sysarray[0].spo = (SPOOLMAT*)CCACALLOC(1,sizeof(SPOOLMAT));
      tri->sysarray[0].spo->is_init=1;
      tri->sysarray[0].spo->ncall=0;
      tri->sysarray[0].spo->is_factored=0;
    }
    else
    {
      // Prepare the matrix
      // Note that we do not copy the matrix as we can feed Spooles
      // directly from the Epetra_CrsMatrix
      if (tri->sysarray_typ[0]!=spoolmatrix ||
          !(tri->sysarray)                  ||
          !(tri->sysarray[0].spo)              )
        dserror("Init phase was not called properly");

      Epetra_CrsMatrix* matrix = (Epetra_CrsMatrix*)tri->matrix;
      SPOOLMAT*         spo    = tri->sysarray[0].spo;
      if (spo->is_factored==0)
      {
        spo->nnz         = matrix->NumGlobalNonzeros();
        spo->numeq_total = tri->numeq_total;
        spo->numeq       = tri->numeq;
        spo->is_init     = 1;
        spo->is_factored = 0;
        spo->ncall       = 0;
        if (spo->update.Typ != ARRAY::cca_XX) amdel(&(spo->update));
        am_alloc_copy(&(tri->update),&(spo->update));
      }
      solver_trilinos_spooles(actsolv,actintra,tri,spo,sol,rhs,option);
    }
  break;
#endif

  //---------------------------------------------------------------------------
  default:
     cout << "solvertyp is " << actsolv->solvertyp << ", see headers/enums.h:SOLVER_TYP\n";
     dserror("Unknown solver, it's either not supported by Trilinos algebra or not compiled in");
  break;
}

#ifdef DEBUG
  dstrc_exit();
#endif
  return;
} /* end of solver_trilinos_control */




/*----------------------------------------------------------------------*
  |                                                          m.gee 9/06 |
  |  routine to control solution with AztecOO                           |
  |  in the case where Trilinos matrices are used internally            |
 *----------------------------------------------------------------------*/
void solve_aztecoo(TRILINOSMATRIX* tri,
                   DIST_VECTOR*    sol,
                   DIST_VECTOR*    rhs,
                   int             option,
                   SOLVAR*         actsolv,
                   FIELD*          actfield,
                   int             disnum)
{
#ifdef DEBUG
  dstrc_enter("solve_aztecoo");
#endif
//-----------------------------------------------------------------------
  // get options
  AZVAR* azvar  = actsolv->azvar;
  // get matrix
  if (!tri->matrix) dserror("Matrix is NULL");
  Epetra_CrsMatrix* matrix = (Epetra_CrsMatrix*)tri->matrix;
  //========================================================== init phase
  if (option==1)
  {
    // create the Epetra_LinearProblem
    if (tri->linearproblem)
    {
      Epetra_LinearProblem* lp = (Epetra_LinearProblem*)tri->linearproblem;
      delete lp;
      tri->linearproblem=NULL;
    }
    Epetra_LinearProblem* lp = new Epetra_LinearProblem();
    tri->linearproblem = (void*)lp;

    // delete the Aztec solver
    if (tri->solver)
    {
      AztecOO* solver = (AztecOO*)tri->solver;
      delete solver;
      tri->solver=NULL;
    }

    // delete the nullspace if present
    if (tri->nullspace)
      delete tri->nullspace;
    tri->nullspace = NULL;

    // create a teuchos parameter list to hold sublists for aztec, ml, etc
    if (tri->params)
    {
      ParameterList* list = (ParameterList*)tri->params;
      delete list;
      tri->params = NULL;
    }
    ParameterList* list = new ParameterList();
    tri->params = (void*)list;

    // set aztec parameters from azvar
    ParameterList& azlist = list->sublist("Aztec Parameters");

    // set solver
    switch(azvar->azsolvertyp)
    {
      case azsolv_CG:       azlist.set("AZ_solver",AZ_cg);       break;
      case azsolv_GMRES:    azlist.set("AZ_solver",AZ_gmres);    break;
      case azsolv_CGS:      azlist.set("AZ_solver",AZ_cgs);      break;
      case azsolv_BiCGSTAB: azlist.set("AZ_solver",AZ_bicgstab); break;
      case azsolv_LU:       azlist.set("AZ_solver",AZ_lu);       break;
      case azsolv_TFQMR:    azlist.set("AZ_solver",AZ_tfqmr);    break;
      default: dserror("Unknown solver for AztecOO");            break;
    }
    azlist.set("AZ_kspace",azvar->azsub);
    // set preconditioner
    switch(azvar->azprectyp)
    {
      case azprec_none:
        azlist.set("AZ_precond",AZ_none);
      break;
      case azprec_ILUT:
        // using ifpack, see below
        azlist.set("AZ_precond",AZ_user_precond);
      break;
      case azprec_ILU:
        // using ifpack, see below
        azlist.set("AZ_precond",AZ_user_precond);
      break;
      case azprec_Jacobi:
        azlist.set("AZ_precond",AZ_Jacobi);
      break;
      case azprec_Neumann:
        azlist.set("AZ_precond",AZ_Neumann);
      break;
      case azprec_Least_Squares:
        azlist.set("AZ_precond",AZ_ls);
      break;
      case azprec_SymmGaussSeidel:
        azlist.set("AZ_precond",AZ_sym_GS);
      break;
      case azprec_LU:
        // using ifpack, see below
        azlist.set("AZ_precond",AZ_user_precond);
      break;
      case azprec_RILU:
        azlist.set("AZ_precond",AZ_dom_decomp);
        azlist.set("AZ_subdomain_solve",AZ_rilu);
        azlist.set("AZ_graph_fill",azvar->azgfill);
      break;
      case azprec_ICC:
        // using ifpack, see below
        azlist.set("AZ_precond",AZ_user_precond);
      break;
      case azprec_ML:
      case azprec_MLfluid:
      case azprec_MLfluid2:
#if (!defined(SOLVE_DIRICH)) || (!defined(SOLVE_DIRICH2))
        if (matrix->Comm().MyPID()==0)
        {
          printf("WARNING:\n");
          printf("You should use SOLVE_DIRICH and SOLVE_DIRICH2 with ML\n");
          printf("(Or have NO nodes with partial Dirichlet BCs)\n");
          fflush(stdout);
        }
#endif
        azlist.set("AZ_precond",AZ_user_precond);
      break;
      default:
        dserror("Unknown preconditioner for AztecOO");
      break;
    }

    // set other parameters
    azlist.set("AZ_max_iter",azvar->aziter);
    azlist.set("AZ_overlap",0);
    azlist.set("AZ_type_overlap",AZ_symmetric);
    azlist.set("AZ_poly_ord",azvar->azpoly);
    if (!azvar->azoutput)
      azlist.set("AZ_output",AZ_none);             // AZ_none AZ_all AZ_warnings AZ_last 10
    else
      azlist.set("AZ_output",azvar->azoutput);
    azlist.set("AZ_diagnostics",AZ_none);          // AZ_none AZ_all
    azlist.set("AZ_conv",azvar->azconv);           // AZ_conv from input
    azlist.set("AZ_tol",azvar->aztol);
    azlist.set("AZ_drop",azvar->azdrop);
    azlist.set("AZ_scaling",AZ_none);              // use epetra scaling instead, see below
    azlist.set("AZ_keep_info",1);

    // set flags
    tri->is_init=1;
    tri->ncall=0;
  }
  //========================================================== solution phase
  else if (option==0)
  {
    // scale linear system (be careful with ML though, try without as well)
    bool scaling_infnorm = false;
    bool scaling_symdiag = false;
    if (azvar->azscal==1)
    {
      scaling_infnorm = false;
      scaling_symdiag = true;
    }
    else if (azvar->azscal==2)
    {
      scaling_infnorm = true;
      scaling_symdiag = false;
    }

    // get parameter list
    if (!tri->params) dserror("AztecOO parameters is NULL");
    ParameterList* list = (ParameterList*)tri->params;

    //----------------------------------------- compute a recreation flag
    bool create = true;
    if (tri->ncall == 0)                       create = true;  // first time
    else if (azvar->azreuse==0)                create = true;  // no reuse from input
    else if (tri->ncall % azvar->azreuse == 0) create = true;  // modulo recreate is true
    else                                       create = false; // reuse

    //----------------------------------------- prepare scaled linear system
    // wrap vectors
    Epetra_Vector x(View,matrix->OperatorDomainMap(),sol->vec.a.dv);
    Epetra_Vector b(View,matrix->OperatorRangeMap(),rhs->vec.a.dv);
    // get linear problem
    if (!tri->linearproblem) dserror("linearproblem is NULL");
    Epetra_LinearProblem* lp = (Epetra_LinearProblem*)tri->linearproblem;
    // set vectors into linear problem
    lp->SetLHS(&x);
    lp->SetRHS(&b);
    lp->SetOperator(matrix);
    // do symmetric infnorm scaling (note that ML might hate this!)
    // most matrices experience slightly better performance with this
    // than with symm. diagonal scaling (ML hates that as well).
    // Of course, the 2 scalings could be introduced as options....
    RefCountPtr<Epetra_Vector> rowsum;
    RefCountPtr<Epetra_Vector> colsum;
    if (scaling_infnorm)
    {
      rowsum = rcp(new Epetra_Vector(matrix->RowMap(),false));
      colsum = rcp(new Epetra_Vector(matrix->RowMap(),false));
      matrix->InvRowSums(*rowsum);
      matrix->InvColSums(*colsum);
      lp->LeftScale(*rowsum);
      lp->RightScale(*colsum);
    }
    RefCountPtr<Epetra_Vector> diag;
    if (scaling_symdiag)
    {
      Epetra_Vector invdiag(matrix->RowMap(),false);
      diag = rcp(new Epetra_Vector(matrix->RowMap(),false));
      matrix->ExtractDiagonalCopy(*diag);
      invdiag.Reciprocal(*diag);
      lp->LeftScale(invdiag);
      lp->RightScale(invdiag);
    }

    //---- get solver and recreate it (currently dies when trying to reuse)
    if (tri->solver)
    {
      AztecOO* tmp = (AztecOO*)tri->solver;
      delete tmp;
    }
    AztecOO* solver = new AztecOO();
    tri->solver = (void*)solver;
    // set parameters to be safe
    ParameterList& azlist = list->sublist("Aztec Parameters");
    // set any parameters
    solver->SetParameters(azlist,false);
    // pass linear problem to solver
    solver->SetProblem(*lp);

    //-------------------------- get IFPACK preconditioner and (re)create it
    if (azvar->azprectyp == azprec_ILU  ||
        azvar->azprectyp == azprec_ILUT ||
        azvar->azprectyp == azprec_ICC  ||
        azvar->azprectyp == azprec_LU   )
    {
      if (create)
      {
        if (tri->prec) // destroy preconditioner
        {
          Ifpack_Preconditioner* prec = (Ifpack_Preconditioner*)tri->prec;
          delete prec;
          tri->prec = NULL;
        }
        // create ifpack parameter list
        ParameterList&  ifpacklist = list->sublist("IFPACK Parameters");
        ifpacklist.set("fact: drop tolerance",azvar->azdrop);
        ifpacklist.set("fact: level-of-fill",azvar->azgfill);
        ifpacklist.set("fact: ilut level-of-fill",azvar->azfill);
        //ifpacklist.set("fact: relax value",0.0);
        //ifpacklist.set("fact: absolute threshold",0.1);
        ifpacklist.set("schwarz: combine mode","Add"); // can be "Add", "Zero", "Insert", "InsertAdd", "Average", "AbsMax"
        ifpacklist.set("schwarz: reordering type","rcm"); // "rcm" or "metis"
        ifpacklist.set("amesos: solver type", "Amesos_Klu"); // can be "Amesos_Klu", "Amesos_Umfpack", "Amesos_Superlu"

        // destroy the copy of the matrix we have stored and create a new copy
        // this is needed to reuse the Ifpack operator
        if (tri->precmatrix)
        {
          Epetra_CrsMatrix* tmp = (Epetra_CrsMatrix*)tri->precmatrix;
          delete tmp;
          tri->precmatrix=NULL;
        }
        // create copy of matrix
        Epetra_CrsMatrix* precmatrix = new Epetra_CrsMatrix(*matrix);
        tri->precmatrix = (void*)precmatrix;

        // create the preconditioner
        string prectype = "";
        if      (azvar->azprectyp == azprec_ILU)  prectype = "ILU";
        else if (azvar->azprectyp == azprec_ILUT) prectype = "ILUT";
        else if (azvar->azprectyp == azprec_ICC)  prectype = "IC";
        else if (azvar->azprectyp == azprec_LU)   prectype = "Amesos";
        else dserror("Unknown type of ifpack asm subdomain preconditioner");

        Ifpack Factory;
        Ifpack_Preconditioner* prec =
                Factory.Create(prectype,precmatrix,azvar->azoverlap);
        prec->SetParameters(ifpacklist);
        prec->Initialize();
        prec->Compute();
        tri->prec = (void*)prec;
        solver->SetPrecOperator(prec);
      }
      else // reuse
      {
        if (!tri->prec) dserror("Cannot reuse, prec is NULL");
        Ifpack_Preconditioner* prec =
                             (Ifpack_Preconditioner*)tri->prec;
        solver->SetPrecOperator(prec);
      }
    }
    //------------------------------ get ML preconditioner and (re)create it
    if (azvar->azprectyp == azprec_ML      ||
        azvar->azprectyp == azprec_MLfluid ||
        azvar->azprectyp == azprec_MLfluid2 )
    {
      if (create)
      {
        if (tri->prec) // destroy ML preconditioner
        {
          ML_Epetra::MultiLevelPreconditioner* prec =
           (ML_Epetra::MultiLevelPreconditioner*)tri->prec;
          delete prec;
          tri->prec = NULL;
        }

        // create ML's parameter list and nullspace
        ParameterList& mllist = list->sublist("ML Parameters");
        create_ml_parameterlist(actsolv,mllist,actfield,disnum,*matrix,&(tri->nullspace));

        // destroy the copy of the matrix we have stored and create a new copy
        // this is needed to reuse the ML operator
        if (tri->precmatrix)
        {
          Epetra_CrsMatrix* tmp = (Epetra_CrsMatrix*)tri->precmatrix;
          delete tmp;
          tri->precmatrix=NULL;
        }

        // create copy of matrix
        Epetra_CrsMatrix* precmatrix = new Epetra_CrsMatrix(*matrix);
        tri->precmatrix = (void*)precmatrix;

        ML_Epetra::MultiLevelPreconditioner* prec =
               new ML_Epetra::MultiLevelPreconditioner(*precmatrix,mllist,true);
        tri->prec = (void*)prec;
        solver->SetPrecOperator(prec);
      }
      else // reuse
      {
        if (!tri->prec) dserror("Cannot reuse, prec is NULL");
        ML_Epetra::MultiLevelPreconditioner* prec =
                             (ML_Epetra::MultiLevelPreconditioner*)tri->prec;
        solver->SetPrecOperator(prec);
      }
    }

    //--------------------------------------------------------- iterate
    solver->Iterate(azvar->aziter,azvar->aztol);

    //------------------------------ check status of solve
    const double* status = solver->GetAztecStatus();
    if (status[AZ_why] != AZ_normal)
    {
      bool resolve = false;
      if (status[AZ_why] == AZ_breakdown)
      {
        if (matrix->Comm().MyPID()==0)
          printf("Numerical breakdown in AztecOO, try again with SuperLU\n");
        resolve = true;
      }
      if (status[AZ_why] == AZ_loss)
      {
        if (matrix->Comm().MyPID()==0)
        {
          printf("RANK 0: AztecOO: Numerical loss of precision occured! continue...\n");
          fprintf(allfiles.out_err,"RANK 0: AztecOO: Numerical loss of precision occured, continue...\n");
        }
        resolve = true;
      }
      if (status[AZ_why] == AZ_ill_cond)
      {
        if (matrix->Comm().MyPID()==0)
        {
          printf("RANK 0: AztecOO: Preconditioning ill-conditioned or singular,\n");
          printf("                 solution is least square ! continue...\n");
          fprintf(allfiles.out_err,"RANK 0: AztecOO: Preconditioning ill-conditioned or singular,\n");
          fprintf(allfiles.out_err,"                 solution is least square ! continue...\n");
        }
        resolve = true;
      }
      if (status[AZ_why] == AZ_maxits)
      {
        if (matrix->Comm().MyPID()==0)
        {
         printf("RANK 0: AztecOO: Maximum number of iterations %d reached \n",azvar->aziter);
         fprintf(allfiles.out_err,"RANK 0: AztecOO: Maximum number of iterations %d reached \n",azvar->aziter);
         fflush(allfiles.out_err);
         fflush(stdout);
        }
        resolve = true;
      }
      if (resolve)
      {
#if 0
        // this is not such a good idea after all...
#ifdef PARALLEL
        if (matrix->Comm().MyPID()==0) cout << "Falling back to SuperLU\n";
        Amesos_Superludist superlusolver(*lp);
        int err = superlusolver.SymbolicFactorization();
        if (err) dserror("SuperLU.SymbolicFactorization() returned %d",err);
        err     = superlusolver.NumericFactorization();
        if (err) dserror("SuperLU.NumericFactorization() returned %d",err);
        err     = superlusolver.Solve();
        if (err) dserror("SuperLU.Solve() returned %d",err);
#else
        cout << "Falling back to Amesos_KLU\n";
        Amesos_Klu klusolver(*lp);
        int err = klusolver.SymbolicFactorization();
        if (err) dserror("Amesos_Klu.SymbolicFactorization() returned %d",err);
        err     = klusolver.NumericFactorization();
        if (err) dserror("Amesos_Klu.NumericFactorization() returned %d",err);
        err     = klusolver.Solve();
        if (err) dserror("Amesos_Klu.Solve() returned %d",err);
#endif
#endif
      }
    }
    // print some statistics
    if (matrix->Comm().MyPID()==0)
    {
    if (actsolv->fieldtyp==structure) fprintf(allfiles.out_err,"Structure: ");
    if (actsolv->fieldtyp==fluid)     fprintf(allfiles.out_err,"Fluid: ");
    if (actsolv->fieldtyp==ale)       fprintf(allfiles.out_err,"Ale: ");
    if (actsolv->fieldtyp==pressure)  fprintf(allfiles.out_err,"Pressure: ");
     fprintf(allfiles.out_err,"AztecOO: unknowns/iterations/time %d  %d  %f\n",
             sol->numeq_total,(int)status[AZ_its],status[AZ_solve_time]);
    }

    //------------------------------------ undo scaling of linear problem
    if (scaling_infnorm)
    {
      Epetra_Vector invrowsum(matrix->RowMap(),false);
      invrowsum.Reciprocal(*rowsum);
      rowsum = null;
      Epetra_Vector invcolsum(matrix->RowMap(),false);
      invcolsum.Reciprocal(*colsum);
      colsum = null;
      lp->LeftScale(invrowsum);
      lp->RightScale(invcolsum);
    }
    if (scaling_symdiag)
    {
      lp->LeftScale(*diag);
      lp->RightScale(*diag);
      diag = null;
    }

    // set flags important for reuse
    tri->is_factored=1;
    tri->ncall++;
  } // option==0
  //========================================================== destroy phase
  else if (option==2) // destroy everything associated with this solver/matrix
  {                   // do NOT destroy the matrix, that's done somewhere else
    if (tri->linearproblem)
    {
      Epetra_LinearProblem* lp = (Epetra_LinearProblem*)tri->linearproblem;
      delete lp;
      tri->linearproblem=NULL;
    }
    if (tri->solver)
    {
      AztecOO* solver = (AztecOO*)tri->solver;
      delete solver;
      tri->solver=NULL;
    }
    if (tri->params)
    {
      ParameterList* list = (ParameterList*)tri->params;
      delete list;
      tri->params = NULL;
    }
    if (tri->precmatrix)
    {
          Epetra_CrsMatrix* tmp = (Epetra_CrsMatrix*)tri->precmatrix;
          delete tmp;
          tri->precmatrix=NULL;
    }
    if (tri->nullspace)
    {
      delete tri->nullspace;
      tri->nullspace = NULL;
    }
    if (tri->prec)
    {
      if (azvar->azprectyp == azprec_ML      ||
          azvar->azprectyp == azprec_MLfluid ||
          azvar->azprectyp == azprec_MLfluid2 )
      {
        ML_Epetra::MultiLevelPreconditioner* prec =
         (ML_Epetra::MultiLevelPreconditioner*)tri->prec;
        delete prec;
        tri->prec = NULL;
      }
      else if (azvar->azprectyp == azprec_ILU  ||
               azvar->azprectyp == azprec_ILUT ||
               azvar->azprectyp == azprec_ICC  ||
               azvar->azprectyp == azprec_LU)
      {
        // this might die because Ifpack_Preconditioner is only a base class
        // to the actual derived type
        // (better dynamic_cast from Epetra_Operator?)
        Ifpack_Preconditioner* prec = (Ifpack_Preconditioner*)tri->prec;
        delete prec;
        tri->prec = NULL;
      }
      else dserror("Cannot destroy preconditioner of unknown type");
    }
  }
  else // option other than 1 or 0 or 2
    dserror("Unknown option flag for AztecOO solver");
//-----------------------------------------------------------------------
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
} /* end of solve_aztecoo */


/*----------------------------------------------------------------------*
  |                                                          m.gee 9/06 |
  |  routine to control solution with Amesos_Klu                                |
  |  in the case where Trilinos matrices are used internally            |
 *----------------------------------------------------------------------*/
void solve_amesos_klu(TRILINOSMATRIX* tri,
                      DIST_VECTOR*    sol,
                      DIST_VECTOR*    rhs,
                      int             option,
                      bool            sym)
{
#ifdef DEBUG
  dstrc_enter("solve_amesos_klu");
#endif
//-----------------------------------------------------------------------
  //========================================================== init phase
  if (option==1)
  {
    // create the Epetra_LinearProblem
    if (tri->linearproblem)
    {
      Epetra_LinearProblem* lp = (Epetra_LinearProblem*)tri->linearproblem;
      delete lp;
    }
    Epetra_LinearProblem* lp = new Epetra_LinearProblem();
    tri->linearproblem = (void*)lp;
    // create the Amesos_KLU solver
    if (tri->solver)
    {
      Amesos_Klu* solver = (Amesos_Klu*)tri->solver;
      delete solver;
    }
    Amesos_Klu* solver = new Amesos_Klu(*lp);
    tri->solver = (void*)solver;
    tri->is_init=1;
  }
  //========================================================== solution phase
  else if (option==0)
  {
    // get matrix
    Epetra_CrsMatrix* matrix = (Epetra_CrsMatrix*)tri->matrix;

    // wrap vectors
    Epetra_Vector x(View,matrix->OperatorDomainMap(),sol->vec.a.dv);
    Epetra_Vector b(View,matrix->OperatorDomainMap(),rhs->vec.a.dv);

    // get solver
    Amesos_Klu* solver = (Amesos_Klu*)tri->solver;

    // get linear problem
    Epetra_LinearProblem* lp = const_cast<Epetra_LinearProblem*>(solver->GetProblem());

    // set vectors into linear problem
    lp->SetLHS(&x);
    lp->SetRHS(&b);

    // matrix has not yet been factored
    if (!tri->is_factored)
    {
      int err=0;
      lp->SetOperator(matrix);
      solver->SetUseTranspose(sym);
      err = solver->SymbolicFactorization();
      if (err) dserror("Amesos_Klu::SymbolicFactorization returned an err");
      err = solver->NumericFactorization();
      if (err) dserror("Amesos_Klu::NumericFactorization returned an err");
      err = solver->Solve();
      if (err) dserror("Amesos_Klu::Solve returned an err");
    }
    // matrix has been factored before
    else
    {
      int err = solver->Solve();
      if (err) dserror("Amesos_Klu::Solve returned an err");
    }
    tri->is_factored=1;
  }
  //========================================================== destroy phase
  else if (option==2) // destroy everything associated with this solver/matrix
  {                   // do NOT destroy the matrix, that's done somewhere else
    if (tri->linearproblem)
    {
      Epetra_LinearProblem* lp = (Epetra_LinearProblem*)tri->linearproblem;
      delete lp;
      tri->linearproblem=NULL;
    }
    if (tri->solver)
    {
      Amesos_Klu* solver = (Amesos_Klu*)tri->solver;
      delete solver;
      tri->solver=NULL;
    }
    if (tri->params)
    {
      ParameterList* list = (ParameterList*)tri->params;
      delete list;
      tri->params = NULL;
    }
    if (tri->precmatrix)
    {
          Epetra_CrsMatrix* tmp = (Epetra_CrsMatrix*)tri->precmatrix;
          delete tmp;
          tri->precmatrix=NULL;
    }
    if (tri->nullspace)
    {
      delete tri->nullspace;
      tri->nullspace = NULL;
    }
    if (tri->prec)
    {
      dserror("Cannot destroy preconditioner of unknown type");
    }
  }
  else
    dserror("Unknown option flag");
//-----------------------------------------------------------------------
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
} /* end of solve_amesos_klu */



#ifdef PARALLEL
/*----------------------------------------------------------------------*
  |                                                          m.gee 9/06 |
  |  routine to control solution with Amesos_SuperLU                    |
  |  in the case where Trilinos matrices are used internally            |
 *----------------------------------------------------------------------*/
void solve_amesos_superlu(TRILINOSMATRIX* tri,
                          DIST_VECTOR*    sol,
                          DIST_VECTOR*    rhs,
                          int             option)
{
#ifdef DEBUG
  dstrc_enter("solve_amesos_superlu");
#endif
//-----------------------------------------------------------------------
  //========================================================== init phase
  if (option==1)
  {
    // create the Epetra_LinearProblem
    if (tri->linearproblem)
    {
      Epetra_LinearProblem* lp = (Epetra_LinearProblem*)tri->linearproblem;
      delete lp;
    }
    Epetra_LinearProblem* lp = new Epetra_LinearProblem();
    tri->linearproblem = (void*)lp;
    // create the Amesos_Umfpack solver
    if (tri->solver)
    {
      Amesos_Superludist* solver = (Amesos_Superludist*)tri->solver;
      delete solver;
    }
    Amesos_Superludist* solver = new Amesos_Superludist(*lp);
    tri->solver = (void*)solver;
    tri->is_init=1;
  }
  //========================================================== solution phase
  else if (option==0)
  {
    // get matrix
    Epetra_CrsMatrix* matrix = (Epetra_CrsMatrix*)tri->matrix;

    // wrap vectors
    Epetra_Vector x(View,matrix->OperatorDomainMap(),sol->vec.a.dv);
    Epetra_Vector b(View,matrix->OperatorDomainMap(),rhs->vec.a.dv);

    // get solver
    Amesos_Superludist* solver = (Amesos_Superludist*)tri->solver;

    // get linear problem
    Epetra_LinearProblem* lp = const_cast<Epetra_LinearProblem*>(solver->GetProblem());

    // set vectors into linear problem
    lp->SetLHS(&x);
    lp->SetRHS(&b);

    // matrix has not yet been factored
    if (!tri->is_factored)
    {
      int err=0;
      lp->SetOperator(matrix);
      // in case of umfpack, the amesos docu says one cannot replace the matrix
      // under the hood, so destroy and recreate
      delete solver;
      solver = new Amesos_Superludist(*lp);
      tri->solver = (void*)solver;
      err = solver->SymbolicFactorization();
      if (err) dserror("Amesos_Superludist::SymbolicFactorization returned an err");
      err = solver->NumericFactorization();
      if (err) dserror("Amesos_Superludist::NumericFactorization returned an err");
      err = solver->Solve();
      if (err) dserror("Amesos_Superludist::Solve returned an err");
    }
    // matrix has been factored before
    else
    {
      int err = solver->Solve();
      if (err) dserror("Amesos_Superludist::Solve returned an err");
    }
    tri->is_factored=1;
  }
  //========================================================== destroy phase
  else if (option==2) // destroy everything associated with this solver/matrix
  {                   // do NOT destroy the matrix, that's done somewhere else
    if (tri->linearproblem)
    {
      Epetra_LinearProblem* lp = (Epetra_LinearProblem*)tri->linearproblem;
      delete lp;
      tri->linearproblem=NULL;
    }
    if (tri->solver)
    {
      Amesos_Superludist* solver = (Amesos_Superludist*)tri->solver;
      delete solver;
      tri->solver=NULL;
    }
    if (tri->params)
    {
      ParameterList* list = (ParameterList*)tri->params;
      delete list;
      tri->params = NULL;
    }
    if (tri->precmatrix)
    {
          Epetra_CrsMatrix* tmp = (Epetra_CrsMatrix*)tri->precmatrix;
          delete tmp;
          tri->precmatrix=NULL;
    }
    if (tri->nullspace)
    {
      delete tri->nullspace;
      tri->nullspace = NULL;
    }
    if (tri->prec)
    {
      dserror("Cannot destroy preconditioner of unknown type");
    }
  }
  else
    dserror("Unknown option flag");
//-----------------------------------------------------------------------
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
} /* end of solve_amesos_superlu */
#endif // PARALLEL


/*----------------------------------------------------------------------*
  |                                                          m.gee 9/06 |
  |  routine to control solution with Amesos_Klu                                |
  |  in the case where Trilinos matrices are used internally            |
 *----------------------------------------------------------------------*/
void solve_amesos_umfpack(TRILINOSMATRIX* tri,
                          DIST_VECTOR*    sol,
                          DIST_VECTOR*    rhs,
                          int             option)
{
#ifdef DEBUG
  dstrc_enter("solve_amesos_umfpack");
#endif
//-----------------------------------------------------------------------
  //========================================================== init phase
  if (option==1)
  {
    // create the Epetra_LinearProblem
    if (tri->linearproblem)
    {
      Epetra_LinearProblem* lp = (Epetra_LinearProblem*)tri->linearproblem;
      delete lp;
    }
    Epetra_LinearProblem* lp = new Epetra_LinearProblem();
    tri->linearproblem = (void*)lp;
    // create the Amesos_Umfpack solver
    if (tri->solver)
    {
      Amesos_Umfpack* solver = (Amesos_Umfpack*)tri->solver;
      delete solver;
    }
    Amesos_Umfpack* solver = new Amesos_Umfpack(*lp);
    tri->solver = (void*)solver;
    tri->is_init=1;
  }
  //========================================================== solution phase
  else if (option==0)
  {
    // get matrix
    Epetra_CrsMatrix* matrix = (Epetra_CrsMatrix*)tri->matrix;

    // wrap vectors
    Epetra_Vector x(View,matrix->OperatorDomainMap(),sol->vec.a.dv);
    Epetra_Vector b(View,matrix->OperatorDomainMap(),rhs->vec.a.dv);

    // get solver
    Amesos_Umfpack* solver = (Amesos_Umfpack*)tri->solver;

    // get linear problem
    Epetra_LinearProblem* lp = const_cast<Epetra_LinearProblem*>(solver->GetProblem());

    // set vectors into linear problem
    lp->SetLHS(&x);
    lp->SetRHS(&b);

    // matrix has not yet been factored
    if (!tri->is_factored)
    {
      int err=0;
      lp->SetOperator(matrix);
      // in case of umfpack, the amesos docu says one cannot replace the matrix
      // under the hood, so destroy and recreate
      delete solver;
      solver = new Amesos_Umfpack(*lp);
      tri->solver = (void*)solver;
      err = solver->SymbolicFactorization();
      if (err) dserror("Amesos_Umfpack::SymbolicFactorization returned an err");
      err = solver->NumericFactorization();
      if (err) dserror("Amesos_Umfpack::NumericFactorization returned an err");
      err = solver->Solve();
      if (err) dserror("Amesos_Umfpack::Solve returned an err");
    }
    // matrix has been factored before
    else
    {
      int err = solver->Solve();
      if (err) dserror("Amesos_Umfpack::Solve returned an err");
    }
    tri->is_factored=1;
  }
  //========================================================== destroy phase
  else if (option==2) // destroy everything associated with this solver/matrix
  {                   // do NOT destroy the matrix, that's done somewhere else
    if (tri->linearproblem)
    {
      Epetra_LinearProblem* lp = (Epetra_LinearProblem*)tri->linearproblem;
      delete lp;
      tri->linearproblem=NULL;
    }
    if (tri->solver)
    {
      Amesos_Umfpack* solver = (Amesos_Umfpack*)tri->solver;
      delete solver;
      tri->solver=NULL;
    }
    if (tri->params)
    {
      ParameterList* list = (ParameterList*)tri->params;
      delete list;
      tri->params = NULL;
    }
    if (tri->precmatrix)
    {
          Epetra_CrsMatrix* tmp = (Epetra_CrsMatrix*)tri->precmatrix;
          delete tmp;
          tri->precmatrix=NULL;
    }
    if (tri->nullspace)
    {
      delete tri->nullspace;
      tri->nullspace = NULL;
    }
    if (tri->prec)
    {
      dserror("Cannot destroy preconditioner of unknown type");
    }
  }
  else
    dserror("Unknown option flag");
//-----------------------------------------------------------------------
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
} /* end of solve_amesos_umfpack */


/*----------------------------------------------------------------------*
  |                                                          m.gee 9/06 |
  |  routine to control solution with Amesos_Lapack                                |
  |  in the case where Trilinos matrices are used internally            |
 *----------------------------------------------------------------------*/
void solve_amesos_lapack(TRILINOSMATRIX* tri,
                         DIST_VECTOR*    sol,
                         DIST_VECTOR*    rhs,
                         int             option,
                         bool            sym)
{
#ifdef DEBUG
  dstrc_enter("solve_amesos_lapack");
#endif
//-----------------------------------------------------------------------
  //========================================================== init phase
  if (option==1)
  {
    // create the Epetra_LinearProblem
    if (tri->linearproblem)
    {
      Epetra_LinearProblem* lp = (Epetra_LinearProblem*)tri->linearproblem;
      delete lp;
    }
    Epetra_LinearProblem* lp = new Epetra_LinearProblem();
    tri->linearproblem = (void*)lp;
    // create the Amesos_Lapack solver
    if (tri->solver)
    {
      Amesos_Lapack* solver = (Amesos_Lapack*)tri->solver;
      delete solver;
    }
    Amesos_Lapack* solver = new Amesos_Lapack(*lp);
    tri->solver = (void*)solver;
    tri->is_init=1;
  }
  //========================================================== solution phase
  else if (option==0)
  {
    // get matrix
    Epetra_CrsMatrix* matrix = (Epetra_CrsMatrix*)tri->matrix;

    // wrap vectors
    Epetra_Vector x(View,matrix->OperatorDomainMap(),sol->vec.a.dv);
    Epetra_Vector b(View,matrix->OperatorDomainMap(),rhs->vec.a.dv);

    // get solver
    Amesos_Lapack* solver = (Amesos_Lapack*)tri->solver;

    // get linear problem
    Epetra_LinearProblem* lp = const_cast<Epetra_LinearProblem*>(solver->GetProblem());

    // set vectors into linear problem
    lp->SetLHS(&x);
    lp->SetRHS(&b);

    // matrix has not yet been factored
    if (!tri->is_factored)
    {
      int err=0;
      lp->SetOperator(matrix);
      // in case of umfpack, the amesos docu says one cannot replace the matrix
      // under the hood, so destroy and recreate
      delete solver;
      solver = new Amesos_Lapack(*lp);
      tri->solver = (void*)solver;
      //solver->SetUseTranspose(sym);
      err = solver->SymbolicFactorization();
      if (err) dserror("Amesos_Lapack::SymbolicFactorization returned an err");
      err = solver->NumericFactorization();
      if (err) dserror("Amesos_Lapack::NumericFactorization returned an err");
      err = solver->Solve();
      if (err) dserror("Amesos_Lapack::Solve returned an err");
    }
    // matrix has been factored before
    else
    {
      int err = solver->Solve();
      if (err) dserror("Amesos_Lapack::Solve returned an err");
    }
    tri->is_factored=1;
  }
  //========================================================== destroy phase
  else if (option==2) // destroy everything associated with this solver/matrix
  {                   // do NOT destroy the matrix, that's done somewhere else
    if (tri->linearproblem)
    {
      Epetra_LinearProblem* lp = (Epetra_LinearProblem*)tri->linearproblem;
      delete lp;
      tri->linearproblem=NULL;
    }
    if (tri->solver)
    {
      Amesos_Lapack* solver = (Amesos_Lapack*)tri->solver;
      delete solver;
      tri->solver=NULL;
    }
    if (tri->params)
    {
      ParameterList* list = (ParameterList*)tri->params;
      delete list;
      tri->params = NULL;
    }
    if (tri->precmatrix)
    {
          Epetra_CrsMatrix* tmp = (Epetra_CrsMatrix*)tri->precmatrix;
          delete tmp;
          tri->precmatrix=NULL;
    }
    if (tri->nullspace)
    {
      delete tri->nullspace;
      tri->nullspace = NULL;
    }
    if (tri->prec)
    {
      dserror("Cannot destroy preconditioner of unknown type");
    }
  }
  else
    dserror("Unknown option flag");
//-----------------------------------------------------------------------
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
} /* end of solve_amesos_lapack */

#endif // TRILINOS_PACKAGE
