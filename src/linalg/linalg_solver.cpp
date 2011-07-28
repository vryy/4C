/*!----------------------------------------------------------------------
\file linalg_solver.cpp

\brief

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
#ifdef CCADISCRET

#undef WRITEOUTSTATISTICS
#ifdef WRITEOUTSTATISTICS
#include "Teuchos_Time.hpp"
#endif

#ifdef PARALLEL
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

#include "../drt_inpar/inpar_solver.H"
#include "linalg_solver.H"
#include "linalg_sparsematrix.H"
#include "../drt_lib/standardtypes_cpp.H"

#include "BelosTypes.hpp" // for Belos verbosity codes

// includes for MLMultiply...
#include <ml_MultiLevelPreconditioner.h>
#include <MLAPI_Operator_Utils.h>
#include <MLAPI_Workspace.h>

#include <az_aztec_defs.h> // for translation of parameters

#include <Teuchos_StandardParameterEntryValidators.hpp>
#include <EpetraExt_Transpose_RowMatrix.h>

#include "solver_directsolver.H"
#include "solver_aztecsolver.H"
#ifdef TRILINOS_DEV
#include "solver_stratimikossolver.H"
#endif

#include <Teuchos_TimeMonitor.hpp>

/*----------------------------------------------------------------------*
 |  ctor (public)                                            mwgee 02/07|
 *----------------------------------------------------------------------*/
LINALG::Solver::Solver(RCP<ParameterList> params,
                       const Epetra_Comm& comm, FILE* outfile) :
comm_(comm),
params_(params),
outfile_(outfile)
{
  Setup();
  return;
}

/*----------------------------------------------------------------------*
 |  ctor (public)                                            mwgee 03/08|
 *----------------------------------------------------------------------*/
LINALG::Solver::Solver(const Epetra_Comm& comm, FILE* outfile) :
comm_(comm),
params_(rcp(new ParameterList())),
outfile_(outfile)
{
  // set the default solver
  Params().set("solver","klu");
  Params().set("symmetric",false);

  // set-up
  Setup();

  return;
}

/*----------------------------------------------------------------------*
 |  ctor (public)                                                  11/08|
 *----------------------------------------------------------------------*/
LINALG::Solver::Solver(const Teuchos::ParameterList& inparams,
                       const Epetra_Comm& comm,
                       FILE* outfile) :
comm_(comm),
params_(rcp(new ParameterList())),
outfile_(outfile)
{
  // set solver parameters
  *params_ = TranslateSolverParameters(inparams);

  // set-up
  Setup();

  return;
}

/*----------------------------------------------------------------------*
 |  set-up of stuff common to all constructors                     11/08|
 *----------------------------------------------------------------------*/
void LINALG::Solver::Setup()
{
  solver_ = Teuchos::null;
}

/*----------------------------------------------------------------------*
 |  dtor (public)                                            mwgee 02/07|
 *----------------------------------------------------------------------*/
LINALG::Solver::~Solver()
{
  // destroy in the right order
  Reset();
}

/*----------------------------------------------------------------------*
 |  reset solver (public)                                    mwgee 02/07|
 *----------------------------------------------------------------------*/
void LINALG::Solver::Reset()
{
  solver_ = Teuchos::null;
}

/*----------------------------------------------------------------------*
 |  << operator                                              mwgee 02/07|
 *----------------------------------------------------------------------*/
ostream& operator << (ostream& os, const LINALG::Solver& solver)
{
  solver.Print(os);
  return os;
}

/*----------------------------------------------------------------------*
 |  print solver (public)                                    mwgee 02/07|
 *----------------------------------------------------------------------*/
void LINALG::Solver::Print(ostream& os) const
{
  if (Comm().MyPID()==0)
  {
    os << "============================LINALG::Solver Parameter List\n";
    os << *params_;
    os << "========================end LINALG::Solver Parameter List\n";
  }
  return;
}

/*----------------------------------------------------------------------*
 |  adapt tolerance (public)                                 mwgee 02/08|
 *----------------------------------------------------------------------*/
void LINALG::Solver::AdaptTolerance(const double desirednlnres,
                                    const double currentnlnres,
                                    const double better)
{
  if (!Params().isSublist("Aztec Parameters")) return;
  const int myrank = Comm().MyPID();
  ParameterList& azlist = Params().sublist("Aztec Parameters");
  int output   = azlist.get<int>("AZ_output",1);
  int convtest = azlist.get<int>("AZ_conv",AZ_noscaled);
  if (convtest != AZ_r0) dserror("Using convergence adaptivity: Use AZ_r0 in input file");
  bool havesavedvalue = azlist.isParameter("AZ_tol save");
  if (!havesavedvalue)
  {
    if (!azlist.isParameter("AZ_tol"))
    {
      cout << azlist;
      dserror("No Aztec tolerance in ParameterList");
    }
    azlist.set<double>("AZ_tol save",azlist.get<double>("AZ_tol",1.e-8));
  }
  double tol = azlist.get<double>("AZ_tol save",1.e-8);
  if (!myrank && output)
    printf("                --- Aztec input   relative tolerance %10.3E\n",tol);
  if (currentnlnres*tol < desirednlnres)
  {
    double tolnew = desirednlnres*better/currentnlnres;
    if (tolnew<tol) tolnew = tol;
    if (!myrank && output && tolnew > tol)
      printf("                *** Aztec adapted relative tolerance %10.3E\n",tolnew);
    azlist.set<double>("AZ_tol",tolnew);
  }

  return;
}

/*----------------------------------------------------------------------*
 |  adapt tolerance (public)                                 mwgee 02/08|
 *----------------------------------------------------------------------*/
void LINALG::Solver::ResetTolerance()
{
  if (!Params().isSublist("Aztec Parameters")) return;
  ParameterList& azlist = Params().sublist("Aztec Parameters");
  bool havesavedvalue = azlist.isParameter("AZ_tol save");
  if (!havesavedvalue) return;
  azlist.set<double>("AZ_tol",azlist.get<double>("AZ_tol save",1.e-8));
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void LINALG::Solver::Setup(
  RefCountPtr<Epetra_Operator>     matrix             ,
  RefCountPtr<Epetra_Vector>       x                  ,
  RefCountPtr<Epetra_Vector>       b                  ,
  bool                             refactor           ,
  bool                             reset              ,
  RefCountPtr<Epetra_MultiVector>  weighted_basis_mean,
  RefCountPtr<Epetra_MultiVector>  kernel_c           ,
  bool                             project            )
{
  TEUCHOS_FUNC_TIME_MONITOR("LINALG::Solver:  1)   Setup");

  // reset data flags on demand
   if (reset)
  {
    Reset();
    refactor = true;
  }

  if ( solver_ == Teuchos::null )
  {
    // decide what solver to use
    std::string solvertype = Params().get("solver","none");

    if ("aztec"  ==solvertype)
    {
      solver_ = Teuchos::rcp( new LINALG::SOLVER::AztecSolver( comm_, Params(), outfile_ ) );
    }
    else if ("belos"  ==solvertype)
    {
      //solver_ = Teuchos::rcp( new LINALG::SOLVER::BelosSolver( comm_, Params(), outfile_));
      BuildBelosSolver(comm_, Params(), outfile_);
    }
    else if ("stratimikos" == solvertype)
    {
#ifdef TRILINOS_DEV
      solver_ = Teuchos::rcp( new LINALG::SOLVER::StratimikosSolver(comm_, Params(), outfile_) );
#else
      dserror("Stratimikos not available as solver. Switch TRILINOS_DEV to ON!");
#endif
    }
    else if ("klu"    ==solvertype or
             "umfpack"==solvertype or
             "superlu"==solvertype or
             "lapack" ==solvertype)
    {
      solver_ = Teuchos::rcp( new LINALG::SOLVER::DirectSolver( solvertype ) );
    }
    else
      dserror("Unknown type of solver");
  }

  solver_->Setup( matrix, x, b, refactor, reset, weighted_basis_mean, kernel_c, project );
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void LINALG::Solver::Solve()
{
  TEUCHOS_FUNC_TIME_MONITOR("LINALG::Solver:  2)   Solve");
  solver_->Solve();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int LINALG::Solver::ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  return solver_->ApplyInverse( X, Y );
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void LINALG::Solver::Solve(
  RefCountPtr<Epetra_Operator>     matrix             ,
  RefCountPtr<Epetra_Vector>       x                  ,
  RefCountPtr<Epetra_Vector>       b                  ,
  bool                             refactor           ,
  bool                             reset              ,
  RefCountPtr<Epetra_MultiVector>  weighted_basis_mean,
  RefCountPtr<Epetra_MultiVector>  kernel_c           ,
  bool                             project            )
{
  Setup( matrix, x, b, refactor, reset, weighted_basis_mean, kernel_c, project );
  Solve();
}


/*----------------------------------------------------------------------*
 |  fix an ML nullspace to match a new map  (public)          mwgee 5/11|
 *----------------------------------------------------------------------*/
void LINALG::Solver::FixMLNullspace(char* field,
                                    const Epetra_Map& oldmap,
                                    const Epetra_Map& newmap,
                                    Teuchos::ParameterList& solveparams)
{
  // there is no ML list, do nothing
  if (!solveparams.isSublist("ML Parameters"))
    return;

  Teuchos::ParameterList& params = solveparams.sublist("ML Parameters");

  const int ndim = params.get("null space: dimension",-1);
  if (ndim==-1) dserror("List does not contain nullspace dimension");

  RCP<vector<double> > ns =
        params.get<RCP<vector<double> > >("nullspace",Teuchos::null);
  if (ns==Teuchos::null) dserror("List does not contain nullspace");
  double* ons = &((*ns)[0]);

  const int olength = (int)ns->size() / ndim;
  if (olength != oldmap.NumMyElements())
    dserror("Nullspace does not match old map length");

  const int nlength = newmap.NumMyElements();

  if (olength==nlength) return; // everything should be ok, do nothing

  if (nlength > olength)
    dserror("New problem size larger than old - full rebuild of nullspace neccessary");

  // Allocate a new nullspace and fill it
  RCP<vector<double> > nsnew = rcp(new vector<double>(nlength*ndim,0.0));
  double* nns = &((*nsnew)[0]);

  for (int i=0; i<nlength; ++i)
  {
    int gid = newmap.GID(i);
    int olid = oldmap.LID(gid);
    if (olid==-1) continue;

    // transfer entries for this dof to new nullspace vector
    for (int j=0; j<ndim; ++j)
      nns[j*ndim+i] = ons[j*ndim+olid];
  }

  // put new nullspace in parameter list
  // this print message can go away at some point
  if (!oldmap.Comm().MyPID())
    printf("Fixing %s ML Nullspace\n",field);
  params.set<RCP<vector<double> > >("nullspace",nsnew);
  params.set<double*>("null space: vectors",nns);

  return;
}

/*----------------------------------------------------------------------------------------------*
 |  Translate BACI dat solver block to Stratimikos paramter list  (public)          tawiesn 7/11|
 *----------------------------------------------------------------------------------------------*/
const Teuchos::ParameterList LINALG::Solver::TranslateToStratimikos(const Teuchos::ParameterList& inparams)
{
  Teuchos::ParameterList outparams;

  // switch type of solver
  switch (DRT::INPUT::IntegralValue<INPAR::SOLVER::SolverType>(inparams,"SOLVER"))
  {
  case INPAR::SOLVER::stratimikos_amesos:   //================================================== Lapack
  {
    outparams.set("Linear Solver Type","Amesos");
    ParameterList& alist = outparams.sublist("Linear Solver Types").sublist("Amesos");
    alist = TranslateAmesosToStratimikos(inparams);
  }
  break;
  case INPAR::SOLVER::stratimikos_belos:    //=================================================== Belos
  {
    outparams.set("Linear Solver Type","Belos");
    ParameterList& alist = outparams.sublist("Linear Solver Types").sublist("Belos");
    alist = TranslateBelosToStratimikos(inparams);
  }
  break;
  case INPAR::SOLVER::stratimikos_aztec:    //================================================= AztecOO
  {
    outparams.set("Linear Solver Type","AztecOO");
    ParameterList& alist = outparams.sublist("Linear Solver Types").sublist("AztecOO");
    alist = TranslateAztecToStratimikos(inparams);
  }
  break;
  default:
    dserror("Unsupported type of solver");
  }

  //------------------------------------- set type of preconditioner
  const int prectyp = DRT::INPUT::IntegralValue<INPAR::SOLVER::AzPrecType>(inparams,"AZPREC");
  switch (prectyp)
  {
  case INPAR::SOLVER::azprec_none:
    outparams.set("Preconditioner Type","None");
  break;
  case INPAR::SOLVER::azprec_ILUT:
  case INPAR::SOLVER::azprec_ILU:
  case INPAR::SOLVER::azprec_Jacobi:
  case INPAR::SOLVER::azprec_SymmGaussSeidel:
  case INPAR::SOLVER::azprec_GaussSeidel:
  {
    outparams.set("Preconditioner Type","Ifpack");
    ParameterList& alist = outparams.sublist("Preconditioner Types").sublist("Ifpack");
    alist = TranslateIfpackToStratimikos(inparams);
  }
  break;
  case INPAR::SOLVER::azprec_TekoSIMPLE:
  case INPAR::SOLVER::azprec_BGSnxn:
    // using Teko
    dserror("Teko not supported by Stratimikos by default");
    break;
  case INPAR::SOLVER::azprec_ML:
  case INPAR::SOLVER::azprec_MLfluid:
  case INPAR::SOLVER::azprec_MLAPI:
  case INPAR::SOLVER::azprec_MLfluid2:
  {
    outparams.set("Preconditioner Type","ML");
    ParameterList& alist = outparams.sublist("Preconditioner Types").sublist("ML");
    alist = TranslateMLToStratimikos(inparams);
  }
  break;
  default:
    dserror("Unknown preconditioner for Stratimikos");
  break;
  }

  return outparams;

}

/*-----------------------------------------------------------------------------------------------*
 |  Interpret AMESOS relevant paramters from BACI dat solver block (public)          tawiesn 7/11|
 *-----------------------------------------------------------------------------------------------*/
const Teuchos::ParameterList LINALG::Solver::TranslateAmesosToStratimikos(const Teuchos::ParameterList& inparams)
{
  Teuchos::ParameterList outparams;

  // switch type of solver
  switch (DRT::INPUT::IntegralValue<INPAR::SOLVER::SolverType>(inparams,"SOLVER"))
  {
  case INPAR::SOLVER::stratimikos_amesos://====================================== Tim Davis' KLU
    outparams.set("Solver Type","Klu");
    outparams.set("Refactorization Policy","RepivotOnRefactorization");
    outparams.sublist("Amesos Settings").set("Reindex",true);// in parallel Amesos KLU needs standard maps (ReIndexing necessary)
    break;
  default:
    break;
  }

  // common Amesos parameters
  outparams.sublist("VerboseObject").set("Output File","none");
  outparams.sublist("VerboseObject").set("Verbosity Level","default");
  const int verbosity= inparams.get<int>("AZOUTPUT");
  if (verbosity == 0){ outparams.sublist("VerboseObject").set("Verbosity Level","none");}
  else if (verbosity > 9){ outparams.sublist("VerboseObject").set("Verbosity Level","extreme");}
  else if (verbosity > 6){ outparams.sublist("VerboseObject").set("Verbosity Level","high");}
  else if (verbosity > 3){ outparams.sublist("VerboseObject").set("Verbosity Level","medium");}
  else if (verbosity > 0){ outparams.sublist("VerboseObject").set("Verbosity Level","low");}

  return outparams;
}

/*----------------------------------------------------------------------------------------------*
 |  Interpret AZTEC relevant paramters from BACI dat solver block (public)          tawiesn 7/11|
 *----------------------------------------------------------------------------------------------*/
const Teuchos::ParameterList LINALG::Solver::TranslateAztecToStratimikos(const Teuchos::ParameterList& inparams)
{
  Teuchos::ParameterList outparams;

  outparams.set("Output Every RHS",false);
  ParameterList& fsolve = outparams.sublist("Forward Solve");
  fsolve.set("Max Iterations", inparams.get<int>("AZITER"));
  fsolve.set("Tolerance", inparams.get<double>("AZTOL"));
  ParameterList& aztec = fsolve.sublist("AztecOO Settings");
  //------------------------------------- set type of AZTEC (internal) preconditioner
  const int azprectyp = DRT::INPUT::IntegralValue<INPAR::SOLVER::AzPrecType>(inparams,"AZPREC");
  switch (azprectyp)
  {
  case INPAR::SOLVER::azprec_none:
    aztec.set("Aztec Preconditioner","none");
  break;
  case INPAR::SOLVER::azprec_ILUT:
    aztec.set("Aztec Preconditioner","ilut");
  break;
  case INPAR::SOLVER::azprec_ILU:
    aztec.set("Aztec Preconditioner","ilu");
  break;
  case INPAR::SOLVER::azprec_Jacobi:
    aztec.set("Aztec Preconditioner","Jacobi");
  break;
  case INPAR::SOLVER::azprec_SymmGaussSeidel:
    aztec.set("Aztec Preconditioner","Symmetric Gauss-Seidel");
  break;
  }
  //------------------------------------- set type of AZTEC solver
  switch (DRT::INPUT::IntegralValue<INPAR::SOLVER::AzSolverType>(inparams,"AZSOLVE"))
  {
  case INPAR::SOLVER::azsolv_CG:       aztec.set("Aztec Solver","CG");       break;
  case INPAR::SOLVER::azsolv_GMRES:    aztec.set("Aztec Solver","GMRES");    break;
  case INPAR::SOLVER::azsolv_CGS:      aztec.set("Aztec Solver","CGS");      break;
  case INPAR::SOLVER::azsolv_BiCGSTAB: aztec.set("Aztec Solver","BiCGStab"); break;
  case INPAR::SOLVER::azsolv_LU:       aztec.set("Aztec Solver","LU");       break;
  case INPAR::SOLVER::azsolv_TFQMR:    aztec.set("Aztec Solver","TFQMR");    break;
  default:
  {
    cout << "flag " << DRT::INPUT::IntegralValue<INPAR::SOLVER::AzSolverType>(inparams,"AZSOLVE") << endl;
    dserror("Unknown solver for Stratimikos_Aztec");            break;
  }
  }
  //------------------------------------ other AZTEC parameters
  switch (DRT::INPUT::IntegralValue<int>(inparams,"AZCONV"))
  {
    case AZ_r0:       aztec.set("Convergence Test","r0"); break;
    case AZ_noscaled: aztec.set("Convergence Test","no scaling"); break;
    case AZ_rhs:      aztec.set("Convergence Test","rhs"); break;
    case AZ_Anorm:    aztec.set("Convergence Test","Anorm"); break;
    case AZ_sol:      aztec.set("Convergence Test","sol"); break;
    default: dserror("Convergence Test not supported by Stratimikos_Aztec"); break;
  }
  aztec.set("Output Frequency",inparams.get<int>("AZOUTPUT"));
  aztec.set("Drop Tolerance",inparams.get<double>("AZDROP"));  // for Aztec ILUT // TODO: remove me: we always use ifpack
  aztec.set("Fill Factor",inparams.get<double>("IFPACKFILL")); // for Aztec ILUT // TODO: remove me: we always use ifpack
  aztec.set("Graph Fill",inparams.get<int>("IFPACKGFILL")); // for Aztec ILU    // TODO: remove me: we always use ifpack
  aztec.set("Overlap",inparams.get<int>("IFPACKOVERLAP"));      // TODO: remove me: we always use ifpack
  aztec.set("Polynomial Order",inparams.get<int>("AZPOLY"));    // for Aztec Polynomial prec
  aztec.set("Size of Krylov Subspace",inparams.get<int>("AZSUB")); // size of Krylov subspace

  return outparams;
}

/*----------------------------------------------------------------------------------------------*
 |  Interpret BELOS relevant paramters from BACI dat solver block (public)          tawiesn 7/11|
 *----------------------------------------------------------------------------------------------*/
const Teuchos::ParameterList LINALG::Solver::TranslateBelosToStratimikos(const Teuchos::ParameterList& inparams)
{
  Teuchos::ParameterList outparams;

  //------------------------------------- set type of BELOS solver
  switch (DRT::INPUT::IntegralValue<INPAR::SOLVER::AzSolverType>(inparams,"AZSOLVE"))
  {
  case INPAR::SOLVER::azsolv_CG:
    {
      outparams.set("Solver Type","Block CG");
      ParameterList& stype = outparams.sublist("Solver Types");
      stype.sublist("Block CG").set("Adaptive Block Size",false);  // no adaptive block size, since # rhs always 1
      stype.sublist("Block CG").set("Block Size",1);               // TODO: number of null space vectors?
      stype.sublist("Block CG").set("Convergence Tolerance", inparams.get<double>("AZTOL"));
      stype.sublist("Block CG").set("Maximum Iterations", inparams.get<int>("AZITER"));
      stype.sublist("Block CG").set("Output Frequency", inparams.get<int>("AZOUTPUT"));
    }
    break;
  case INPAR::SOLVER::azsolv_GMRES:
    {
      outparams.set("Solver Type","Block GMRES");
      ParameterList& stype = outparams.sublist("Solver Types");
      stype.sublist("Block GMRES").set("Adaptive Block Size",false);  // no adaptive block size, since # rhs always 1
      stype.sublist("Block GMRES").set("Block Size",1);               // TODO: number of null space vectors?
      stype.sublist("Block GMRES").set("Convergence Tolerance", inparams.get<double>("AZTOL"));
      stype.sublist("Block GMRES").set("Maximum Iterations", inparams.get<int>("AZITER"));
      stype.sublist("Block GMRES").set("Num Blocks", inparams.get<int>("AZSUB")); // number of Krylov space vectors
      stype.sublist("Block GMRES").set("Output Frequency", inparams.get<int>("AZOUTPUT"));

    }
    break;
    default:
    {
      cout << "flag " << DRT::INPUT::IntegralValue<INPAR::SOLVER::AzSolverType>(inparams,"AZSOLVE") << endl;
      dserror("Unknown solver for Stratimikos_Belos");            break;
    }
  }

  return outparams;
}

/*-----------------------------------------------------------------------------------------------*
 |  Interpret IFPACK relevant paramters from BACI dat solver block (public)          tawiesn 7/11|
 *-----------------------------------------------------------------------------------------------*/
const Teuchos::ParameterList LINALG::Solver::TranslateIfpackToStratimikos(const Teuchos::ParameterList& inparams)
{
  Teuchos::ParameterList outparams;

  outparams.set("Overlap", inparams.get<int>("IFPACKOVERLAP"));

  //------------------------------------- set type of preconditioner
  const int prectyp = DRT::INPUT::IntegralValue<INPAR::SOLVER::AzPrecType>(inparams,"AZPREC");
  switch (prectyp)
  {
  case INPAR::SOLVER::azprec_ILUT:    outparams.set("Prec Type", "ILUT"); break;
  case INPAR::SOLVER::azprec_ILU:     outparams.set("Prec Type", "ILU");  break;
  case INPAR::SOLVER::azprec_Jacobi:
  {
      outparams.set("Prec Type","point relaxation");
  }
  break;
  case INPAR::SOLVER::azprec_SymmGaussSeidel:
  {
    outparams.set("Prec Type","point relaxation");
  }
  break;
  case INPAR::SOLVER::azprec_GaussSeidel:
  {
    outparams.set("Prec Type","point relaxation");
  }
  break;
  }


  ParameterList& settings = outparams.sublist("Ifpack Settings");
  settings = TranslateBACIToIfpack(inparams);

  return outparams;
}

/*-------------------------------------------------------------------------------------------*
 |  Interpret ML relevant paramters from BACI dat solver block (public)          tawiesn 7/11|
 *-------------- ----------------------------------------------------------------------------*/
const Teuchos::ParameterList LINALG::Solver::TranslateMLToStratimikos(const Teuchos::ParameterList& inparams)
{
  Teuchos::ParameterList outparams;

  outparams.set("Base Method Defaults", "SA");

  ParameterList& settings = outparams.sublist("ML Settings");
  settings = TranslateBACIToML(inparams,NULL); // no downwinding Gauss Seidel within ML

  return outparams;
}

/*--------------------------------------------------------------------------------------------------------*
 |  Translate IFPACK part of BACI dat solver block to IFPACK parameter list (public)          tawiesn 7/11|
 *--------------------------------------------------------------------------------------------------------*/
const Teuchos::ParameterList LINALG::Solver::TranslateBACIToIfpack(const Teuchos::ParameterList& inparams)
{
  Teuchos::ParameterList ifpacklist;

  ifpacklist.set("relaxation: damping factor",inparams.get<double>("AZOMEGA"));
  ifpacklist.set("fact: drop tolerance",inparams.get<double>("AZDROP"));
  ifpacklist.set("fact: level-of-fill",inparams.get<int>("IFPACKGFILL"));
  ifpacklist.set("fact: ilut level-of-fill",inparams.get<double>("IFPACKFILL"));
  ifpacklist.set("partitioner: overlap",inparams.get<int>("IFPACKOVERLAP"));
  ifpacklist.set("schwarz: combine mode",inparams.get<string>("IFPACKCOMBINE")); // can be "Zero", "Add", "Insert"
  ifpacklist.set("schwarz: reordering type","rcm"); // "rcm" or "metis" or "amd"
  ifpacklist.set("amesos: solver type", "Amesos_Klu"); // can be "Amesos_Klu", "Amesos_Umfpack", "Amesos_Superlu"

  //------------------------------------- set type of preconditioner
  const int prectyp = DRT::INPUT::IntegralValue<INPAR::SOLVER::AzPrecType>(inparams,"AZPREC");
  switch (prectyp)
  {
  case INPAR::SOLVER::azprec_Jacobi:
  {
    ifpacklist.set("relaxation: type","Jacobi");
    ifpacklist.set("relaxation: sweeps",inparams.get<int>("IFPACKGFILL"));
    ifpacklist.set("relaxation: damping factor",inparams.get<double>("AZOMEGA"));
  }
  break;
  case INPAR::SOLVER::azprec_SymmGaussSeidel:
  {
    ifpacklist.set("relaxation: type","symmetric Gauss-Seidel");
    ifpacklist.set("relaxation: sweeps",inparams.get<int>("IFPACKGFILL")); // misuse IFPACKGFILL parameter
  }
  break;
  case INPAR::SOLVER::azprec_GaussSeidel:
  {
    ifpacklist.set("relaxation: type","Gauss-Seidel");
    ifpacklist.set("relaxation: sweeps",inparams.get<int>("IFPACKGFILL")); // misuse IFPACKGFILL parameter
  }
  break;
  case INPAR::SOLVER::azprec_DownwindGaussSeidel:
  {
    // in case of downwinding prevent ifpack from again reordering
    ifpacklist.set("schwarz: reordering type","none");
    ifpacklist.set("relaxation: type","Gauss-Seidel");
    ifpacklist.set("relaxation: sweeps",inparams.get<int>("IFPACKGFILL")); // misuse IFPACKGFILL parameter
  }
  break;
  }

  return ifpacklist;
}

/*------------------------------------------------------------------------------------------------*
 |  Translate ML part of BACI dat solver block to ML parameter list (public)          tawiesn 7/11|
 *------------------------------------------------------------------------------------------------*/
const Teuchos::ParameterList LINALG::Solver::TranslateBACIToML(const Teuchos::ParameterList& inparams, Teuchos::ParameterList* azlist)
{
  Teuchos::ParameterList mllist;

  ML_Epetra::SetDefaults("SA",mllist);
  const int prectyp = DRT::INPUT::IntegralValue<INPAR::SOLVER::AzPrecType>(inparams,"AZPREC");
  switch (prectyp)
  {
  case INPAR::SOLVER::azprec_ML: // do nothing, this is standard
  break;
  case INPAR::SOLVER::azprec_MLAPI: // set flag to use mlapi operator
    mllist.set<bool>("LINALG::AMG_Operator",true);
  break;
  case INPAR::SOLVER::azprec_MLfluid: // unsymmetric, unsmoothed restriction
    mllist.set("aggregation: use tentative restriction",true);
  break;
  case INPAR::SOLVER::azprec_MLfluid2: // full Pretrov-Galerkin unsymmetric smoothed
    mllist.set("energy minimization: enable",true);
    mllist.set("energy minimization: type",3); // 1,2,3 cheap -> expensive
    mllist.set("aggregation: block scaling",false);
  break;
  case INPAR::SOLVER::azprec_AMG:   // set flag to use AMGPreconditioner
    mllist.set<bool>("LINALG::AMGPreconditioner",true);
  break;
  default: dserror("Unknown type of ml preconditioner");
  }
  mllist.set("ML output"                       ,inparams.get<int>("ML_PRINT"));
  if (inparams.get<int>("ML_PRINT")==10)
    mllist.set("print unused"                  ,1);
  else
    mllist.set("print unused"                  ,-2);
  mllist.set("increasing or decreasing"        ,"increasing");
  mllist.set("coarse: max size"                ,inparams.get<int>("ML_MAXCOARSESIZE"));
  mllist.set("max levels"                      ,inparams.get<int>("ML_MAXLEVEL"));
  mllist.set("smoother: pre or post"           ,"both");
  mllist.set("aggregation: threshold"          ,inparams.get<double>("ML_PROLONG_THRES"));
  mllist.set("aggregation: damping factor"     ,inparams.get<double>("ML_PROLONG_SMO"));
  mllist.set("aggregation: nodes per aggregate",inparams.get<int>("ML_AGG_SIZE"));
  // override the default sweeps=2 with a default sweeps=1
  // individual level sweeps are set below
  mllist.set("smoother: sweeps",1);
  // save memory if this is an issue, make ML use single precision
  //mllist.set("low memory usage",true);
  switch (DRT::INPUT::IntegralValue<int>(inparams,"ML_COARSEN"))
  {
    case 0:  mllist.set("aggregation: type","Uncoupled");  break;
    case 1:  mllist.set("aggregation: type","METIS");      break;
    case 2:  mllist.set("aggregation: type","VBMETIS");    break;
    case 3:  mllist.set("aggregation: type","MIS");        break;
    default: dserror("Unknown type of coarsening for ML"); break;
  }

  // set ml smoothers
  const int mlmaxlevel = inparams.get<int>("ML_MAXLEVEL");
  // create vector of integers containing smoothing steps/polynomial order of level
  std::vector<int> mlsmotimessteps;
  {
    std::istringstream mlsmotimes(Teuchos::getNumericStringParameter(inparams,"ML_SMOTIMES"));
    std::string word;
    while (mlsmotimes >> word)
      mlsmotimessteps.push_back(std::atoi(word.c_str()));
  }

  if ((int)mlsmotimessteps.size() < mlmaxlevel)
    dserror("Not enough smoothing steps ML_SMOTIMES=%d, must be larger/equal than ML_MAXLEVEL=%d\n",
            mlsmotimessteps.size(),mlmaxlevel);

  for (int i=0; i<mlmaxlevel-1; ++i)
  {
    char levelstr[11];
    sprintf(levelstr,"(level %d)",i);
    ParameterList& smolevelsublist = mllist.sublist("smoother: list "+(string)levelstr);
    int type;
    double damp;
    if (i==0)
    {
      type = DRT::INPUT::IntegralValue<int>(inparams,"ML_SMOOTHERFINE");
      damp = inparams.get<double>("ML_DAMPFINE");
    }
    else if (i < mlmaxlevel-1)
    {
      type = DRT::INPUT::IntegralValue<int>(inparams,"ML_SMOOTHERMED");
      damp = inparams.get<double>("ML_DAMPMED");
    }
    else
    {
      type = DRT::INPUT::IntegralValue<int>(inparams,"ML_SMOOTHERCOARSE");
      damp = inparams.get<double>("ML_DAMPCOARSE");
    }
    switch (type)
    {
    case 0: // SGS
      smolevelsublist.set("smoother: type"                        ,"symmetric Gauss-Seidel");
      smolevelsublist.set("smoother: sweeps"                      ,mlsmotimessteps[i]);
      smolevelsublist.set("smoother: damping factor"              ,damp);
    break;
    case 7: // GS
      smolevelsublist.set("smoother: type"                        ,"Gauss-Seidel");
      smolevelsublist.set("smoother: sweeps"                      ,mlsmotimessteps[i]);
      smolevelsublist.set("smoother: damping factor"              ,damp);
    break;
    case 8: // DGS
      smolevelsublist.set("smoother: type"                        ,"Gauss-Seidel");
      smolevelsublist.set("smoother: sweeps"                      ,mlsmotimessteps[i]);
      smolevelsublist.set("smoother: damping factor"              ,damp);
      if (azlist!=NULL)
      {
        azlist->set<bool>("downwinding",true);
        azlist->set<double>("downwinding tau",inparams.get<double>("DWINDTAU"));
      }
      else
      {
        cout << "WARNING: cannot set parameters for Downwinding Gauss Seidel" << endl;
      }
      {
        ParameterList& ifpacklist = mllist.sublist("smoother: ifpack list");
        ifpacklist.set("schwarz: reordering type","true");
      }
    break;
    case 1: // Jacobi
      smolevelsublist.set("smoother: type"                        ,"Jacobi");
      smolevelsublist.set("smoother: sweeps"                      ,mlsmotimessteps[i]);
      smolevelsublist.set("smoother: damping factor"              ,damp);
    break;
    case 2: // Chebychev
      smolevelsublist.set("smoother: type"                        ,"MLS");
      smolevelsublist.set("smoother: MLS polynomial order"        ,mlsmotimessteps[i]);
    break;
    case 3: // MLS
      smolevelsublist.set("smoother: type"                        ,"MLS");
      smolevelsublist.set("smoother: MLS polynomial order"        ,-mlsmotimessteps[i]);
    break;
    case 4: // Ifpack's ILU
    {
      smolevelsublist.set("smoother: type"                        ,"IFPACK");
      smolevelsublist.set("smoother: ifpack type"                 ,"ILU");
      smolevelsublist.set("smoother: ifpack overlap"              ,inparams.get<int>("IFPACKOVERLAP"));
      smolevelsublist.set<double>("smoother: ifpack level-of-fill",(double)mlsmotimessteps[i]);
      ParameterList& ifpacklist = mllist.sublist("smoother: ifpack list");
      ifpacklist.set("schwarz: reordering type","rcm"); // "rcm" or "metis" or "amd" or "true"
      ifpacklist.set("schwarz: combine mode",inparams.get<string>("IFPACKCOMBINE")); // can be "Zero", "Insert", "Add"
      ifpacklist.set("partitioner: overlap",inparams.get<int>("IFPACKOVERLAP"));
    }
    break;
    case 5: // Amesos' KLU
      smolevelsublist.set("smoother: type"                        ,"Amesos-KLU");
    break;
    case 9: // Amesos' Umfpack
      smolevelsublist.set("smoother: type"                        ,"Amesos-UMFPACK");
    break;
#ifdef PARALLEL
    case 6: // Amesos' SuperLU_Dist
      smolevelsublist.set("smoother: type"                        ,"Amesos-Superludist");
    break;
#endif
    default: dserror("Unknown type of smoother for ML: tuple %d",type); break;
    } // switch (type)
  } // for (int i=0; i<azvar->mlmaxlevel-1; ++i)

  // set coarse grid solver
  const int coarse = mlmaxlevel-1;
  switch (DRT::INPUT::IntegralValue<int>(inparams,"ML_SMOOTHERCOARSE"))
  {
    case 0:
      mllist.set("coarse: type"          ,"symmetric Gauss-Seidel");
      mllist.set("coarse: sweeps"        ,mlsmotimessteps[coarse]);
      mllist.set("coarse: damping factor",inparams.get<double>("ML_DAMPCOARSE"));
    break;
    case 7:
      mllist.set("coarse: type"          ,"Gauss-Seidel");
      mllist.set("coarse: sweeps"        ,mlsmotimessteps[coarse]);
      mllist.set("coarse: damping factor",inparams.get<double>("ML_DAMPCOARSE"));
    break;
    case 8:
      mllist.set("coarse: type"          ,"Gauss-Seidel");
      mllist.set("coarse: sweeps"        ,mlsmotimessteps[coarse]);
      mllist.set("coarse: damping factor",inparams.get<double>("ML_DAMPCOARSE"));
      if (azlist!=NULL)
      {
        azlist->set<bool>("downwinding",true);
        azlist->set<double>("downwinding tau",inparams.get<double>("DWINDTAU"));
      }
      else
      {
        cout << "WARNING: cannot set parameters for Downwinding Gauss Seidel" << endl;
      }
      {
        ParameterList& ifpacklist = mllist.sublist("smoother: ifpack list");
        ifpacklist.set("schwarz: reordering type","true");
      }
    break;
    case 1:
      mllist.set("coarse: type"          ,"Jacobi");
      mllist.set("coarse: sweeps"        ,mlsmotimessteps[coarse]);
      mllist.set("coarse: damping factor",inparams.get<double>("ML_DAMPCOARSE"));
    break;
    case 2:
      mllist.set("coarse: type"                ,"MLS");
      mllist.set("coarse: MLS polynomial order",mlsmotimessteps[coarse]);
    break;
    case 3:
      mllist.set("coarse: type"                ,"MLS");
      mllist.set("coarse: MLS polynomial order",-mlsmotimessteps[coarse]);
    break;
    case 4:
    {
      mllist.set("coarse: type"          ,"IFPACK");
      mllist.set("coarse: ifpack type"   ,"ILU");
      mllist.set("coarse: ifpack overlap",0);
      mllist.set<double>("coarse: ifpack level-of-fill",(double)mlsmotimessteps[coarse]);
      ParameterList& ifpacklist = mllist.sublist("coarse: ifpack list");
      ifpacklist.set<int>("fact: level-of-fill",mlsmotimessteps[coarse]);
      ifpacklist.set("schwarz: reordering type","rcm");
      ifpacklist.set("schwarz: combine mode",inparams.get<string>("IFPACKCOMBINE")); // can be "Zero", "Insert", "Add"
      ifpacklist.set("partitioner: overlap",inparams.get<int>("IFPACKOVERLAP"));
    }
    break;
    case 5:
      mllist.set("coarse: type","Amesos-KLU");
    break;
    case 9:
      mllist.set("coarse: type","Amesos-UMFPACK");
    break;
    case 6:
      mllist.set("coarse: type","Amesos-Superludist");
    break;
    default: dserror("Unknown type of coarse solver for ML"); break;
  } // switch (azvar->mlsmotype_coarse)
  // default values for nullspace
  mllist.set("PDE equations",1);
  mllist.set("null space: dimension",1);
  mllist.set("null space: type","pre-computed");
  mllist.set("null space: add default vectors",false);
  mllist.set<double*>("null space: vectors",NULL);
#if defined(PARALLEL) && defined(PARMETIS)
  mllist.set("repartition: enable",1);
  mllist.set("repartition: partitioner","ParMETIS");
  mllist.set("repartition: max min ratio",1.3);
  mllist.set("repartition: min per proc",3000);
#endif

  return mllist;
}

const Teuchos::ParameterList LINALG::Solver::TranslateBACIToTeko(const Teuchos::ParameterList& inparams)
{
  Teuchos::ParameterList tekolist;

  const int prectyp = DRT::INPUT::IntegralValue<INPAR::SOLVER::AzPrecType>(inparams,"AZPREC");
  switch (prectyp)
  {
  case INPAR::SOLVER::azprec_TekoSIMPLE:
  {
    tekolist.set("Prec Type","SIMPLE");
    tekolist.set("alpha",inparams.get<double>("SIMPLE_DAMPING"));
  }
  break;
  case INPAR::SOLVER::azprec_BGSnxn:
  {
    tekolist.set("Prec Type","BGS");
  }
  break;
  default:
    dserror("LINALG::Solver::TranslateBACIToTeko: wrong type of preconditioner");
    break;
  }

  return tekolist;
}

/*----------------------------------------------------------------------*
 |  translate solver parameters (public)               mwgee 02/07,11/08|
 *----------------------------------------------------------------------*/
const Teuchos::ParameterList LINALG::Solver::TranslateSolverParameters(const ParameterList& inparams)
{
  TEUCHOS_FUNC_TIME_MONITOR("LINALG::Solver:  0)   TranslateSolverParameters");
  // HINT:
  // input parameter inparams.get<int>("AZGRAPH") is not retrieved

  // make empty output parameters
  Teuchos::ParameterList outparams;

  // switch type of solver
  switch (DRT::INPUT::IntegralValue<INPAR::SOLVER::SolverType>(inparams,"SOLVER"))
  {
#ifdef PARALLEL
  case INPAR::SOLVER::superlu://============================== superlu solver (parallel only)
    outparams.set("solver","superlu");
    outparams.set("symmetric",false);
  break;
#endif
  case INPAR::SOLVER::amesos_klu_sym://====================================== Tim Davis' KLU
    outparams.set("solver","klu");
    outparams.set("symmetric",true);
  break;
  case INPAR::SOLVER::amesos_klu_nonsym://=================================== Tim Davis' KLU
    outparams.set("solver","klu");
    outparams.set("symmetric",false);
  break;
  case INPAR::SOLVER::umfpack://========================================= Tim Davis' Umfpack
    outparams.set("solver","umfpack");
    outparams.set("symmetric",false);
  break;
  case INPAR::SOLVER::lapack_sym://================================================== Lapack
    outparams.set("solver","lapack");
    outparams.set("symmetric",true);
  break;
  case INPAR::SOLVER::lapack_nonsym://=============================================== Lapack
    outparams.set("solver","lapack");
    outparams.set("symmetric",false);
  break;
  case INPAR::SOLVER::stratimikos_amesos: //===================================== Stratimikos
  case INPAR::SOLVER::stratimikos_aztec:  //===================================== Stratimikos
  case INPAR::SOLVER::stratimikos_belos:  //===================================== Stratimikos
  {
    outparams.set("solver","stratimikos");
    outparams.set("xml file",inparams.get<string>("STRATIMIKOS_XMLFILE"));
    ParameterList& stratimikoslist = outparams.sublist("Stratimikos Parameters");
    stratimikoslist = TranslateToStratimikos(inparams);
  }
  break;
  case INPAR::SOLVER::belos:     //=================================================== Belos
  {
    outparams.set("solver","belos");
    outparams.set("symmetric",false);
    ParameterList& beloslist = outparams.sublist("Belos Parameters");
    //--------------------------------------------- set type of solver
    switch (DRT::INPUT::IntegralValue<INPAR::SOLVER::AzSolverType>(inparams,"AZSOLVE"))
    {
    case INPAR::SOLVER::azsolv_CG:       beloslist.set("belos_solver",AZ_cg);       break;
    case INPAR::SOLVER::azsolv_GMRES:    beloslist.set("belos_solver",AZ_gmres);    break;
    default:
    {
      cout << "flag " << DRT::INPUT::IntegralValue<INPAR::SOLVER::AzSolverType>(inparams,"AZSOLVE") << endl;
      dserror("Unknown solver for Belos");            break;
    }
    }
    //------------------------------------- set type of preconditioner
    const int azprectyp = DRT::INPUT::IntegralValue<INPAR::SOLVER::AzPrecType>(inparams,"AZPREC");
    switch (azprectyp)
    {
    case INPAR::SOLVER::azprec_none:
      beloslist.set("Preconditioner Type","none");
    break;
    case INPAR::SOLVER::azprec_ILUT:
      // using ifpack
      beloslist.set("Preconditioner Type","ILUT");
    break;
    case INPAR::SOLVER::azprec_ILU:
      // using ifpack
      beloslist.set("Preconditioner Type","ILU");
    break;
    case INPAR::SOLVER::azprec_TekoSIMPLE:
    case INPAR::SOLVER::azprec_BGSnxn:
      // using Teko
      beloslist.set("Preconditioner Type","Teko");
      break;
    case INPAR::SOLVER::azprec_Jacobi:
    case INPAR::SOLVER::azprec_SymmGaussSeidel:
    case INPAR::SOLVER::azprec_GaussSeidel:
      // using ifpack
      beloslist.set("Preconditioner Type","point relaxation");
    break;
    case INPAR::SOLVER::azprec_DownwindGaussSeidel:
      // using ifpack
      beloslist.set("Preconditioner Type","point relaxation");
      beloslist.set<bool>("downwinding",true);
      beloslist.set<double>("downwinding tau",inparams.get<double>("DWINDTAU"));
      cout << "warning: downwinding not supported?" << endl;
    break;
    case INPAR::SOLVER::azprec_LU:
      // using ifpack
      beloslist.set("Preconditioner Type","Amesos");
    break;
    case INPAR::SOLVER::azprec_ML:
    case INPAR::SOLVER::azprec_MLfluid:
    case INPAR::SOLVER::azprec_MLAPI:
    case INPAR::SOLVER::azprec_MLfluid2:
      beloslist.set("Preconditioner Type","ML");
      break;
    case INPAR::SOLVER::azprec_AMGBS:
    case INPAR::SOLVER::azprec_AMG:
    case INPAR::SOLVER::azprec_BGS2x2:
      beloslist.set("Preconditioner Type","ML");
      break;
    break;
    default:
      dserror("Unknown preconditioner for Belos");
    break;
    }
    //------------------------------------- set other belos parameters
     beloslist.set("Num Blocks",inparams.get<int>("AZSUB"));
     beloslist.set("Maximum Iterations",inparams.get<int>("AZITER"));
     beloslist.set("Adaptive Block Size",true);
     const int verbosity= inparams.get<int>("AZOUTPUT");
     if (verbosity > 0)
       beloslist.set("Verbosity", Belos::Errors);
     else if (verbosity > 1)
       beloslist.set("Verbosity", Belos::Errors + Belos::Warnings);
     else if (verbosity > 3)
       beloslist.set("Verbosity", Belos::Errors + Belos::Warnings + Belos::StatusTestDetails);
     else if (verbosity > 6)
       beloslist.set("Verbosity", Belos::Errors + Belos::Warnings + Belos::StatusTestDetails + Belos::TimingDetails);
     beloslist.set("Convergence Tolerance",inparams.get<double>("AZTOL"));
     //-------------------------------- set parameters for Ifpack if used
     if (azprectyp == INPAR::SOLVER::azprec_ILU  ||
         azprectyp == INPAR::SOLVER::azprec_ILUT ||
         azprectyp == INPAR::SOLVER::azprec_ICC  ||
         azprectyp == INPAR::SOLVER::azprec_LU   ||
         azprectyp == INPAR::SOLVER::azprec_SymmGaussSeidel ||
         azprectyp == INPAR::SOLVER::azprec_GaussSeidel ||
         azprectyp == INPAR::SOLVER::azprec_DownwindGaussSeidel ||
         azprectyp == INPAR::SOLVER::azprec_Jacobi)
     {
       ParameterList& ifpacklist = outparams.sublist("IFPACK Parameters");
       ifpacklist = TranslateBACIToIfpack(inparams);
     }
     //------------------------------------- set parameters for Teko if used
     if (azprectyp == INPAR::SOLVER::azprec_TekoSIMPLE ||
         azprectyp == INPAR::SOLVER::azprec_BGSnxn)
     {
       ParameterList& tekolist = outparams.sublist("Teko Parameters");
       tekolist = TranslateBACIToTeko(inparams);
     }
     //------------------------------------- set parameters for ML if used
     if (azprectyp == INPAR::SOLVER::azprec_ML       ||
         azprectyp == INPAR::SOLVER::azprec_MLfluid  ||
         azprectyp == INPAR::SOLVER::azprec_MLfluid2 ||
         azprectyp == INPAR::SOLVER::azprec_MLAPI ||
         azprectyp == INPAR::SOLVER::azprec_AMG)          // new AMG preconditioner
     {
       ParameterList& mllist = outparams.sublist("ML Parameters");
       mllist = TranslateBACIToML(inparams,&beloslist);
     } // if ml preconditioner

    break;
  }
  case INPAR::SOLVER::aztec_msr://================================================= AztecOO
  {
    outparams.set("solver","aztec");
    outparams.set("symmetric",false);
    ParameterList& azlist = outparams.sublist("Aztec Parameters");
    //--------------------------------- set scaling of linear problem
    const int azscal = DRT::INPUT::IntegralValue<int>(inparams,"AZSCAL");
    if (azscal==1)
      azlist.set("scaling","symmetric");
    else if (azscal==2)
      azlist.set("scaling","infnorm");
    else
      azlist.set("scaling","none");
    //--------------------------------------------- set type of solver
    switch (DRT::INPUT::IntegralValue<INPAR::SOLVER::AzSolverType>(inparams,"AZSOLVE"))
    {
    case INPAR::SOLVER::azsolv_CG:       azlist.set("AZ_solver",AZ_cg);       break;
    case INPAR::SOLVER::azsolv_GMRES:    azlist.set("AZ_solver",AZ_gmres);    break;
    case INPAR::SOLVER::azsolv_GMRESR:   azlist.set("AZ_solver",AZ_GMRESR);   break;
    case INPAR::SOLVER::azsolv_CGS:      azlist.set("AZ_solver",AZ_cgs);      break;
    case INPAR::SOLVER::azsolv_BiCGSTAB: azlist.set("AZ_solver",AZ_bicgstab); break;
    case INPAR::SOLVER::azsolv_LU:       azlist.set("AZ_solver",AZ_lu);       break;
    case INPAR::SOLVER::azsolv_TFQMR:    azlist.set("AZ_solver",AZ_tfqmr);    break;
    default:
    {
      cout << "flag " << DRT::INPUT::IntegralValue<INPAR::SOLVER::AzSolverType>(inparams,"AZSOLVE") << endl;
      dserror("Unknown solver for AztecOO");            break;
    }
    }
    //------------------------------------- set type of preconditioner
    const int azprectyp = DRT::INPUT::IntegralValue<INPAR::SOLVER::AzPrecType>(inparams,"AZPREC");
    switch (azprectyp)
    {
    case INPAR::SOLVER::azprec_none:
      azlist.set("AZ_precond",AZ_none);
      azlist.set("AZ_subdomain_solve",AZ_none);
      azlist.set("AZ_precond",AZ_none);
      azlist.set("Preconditioner Type","none");
    break;
    case INPAR::SOLVER::azprec_ILUT:
      // using ifpack
      azlist.set("AZ_precond",AZ_user_precond);
      azlist.set("Preconditioner Type","ILUT");
    break;
    case INPAR::SOLVER::azprec_ILU:
      // using ifpack
      azlist.set("AZ_precond",AZ_user_precond);
      azlist.set("Preconditioner Type","ILU");
    break;
    case INPAR::SOLVER::azprec_Neumann:
      azlist.set("AZ_precond",AZ_Neumann);
    break;
    case INPAR::SOLVER::azprec_Least_Squares:
      azlist.set("AZ_precond",AZ_ls);
    break;
    case INPAR::SOLVER::azprec_Jacobi:
      // using ifpack
      azlist.set("AZ_precond",AZ_user_precond);
      azlist.set("Preconditioner Type","point relaxation");
    break;
    case INPAR::SOLVER::azprec_SymmGaussSeidel:
      // using ifpack
      azlist.set("AZ_precond",AZ_user_precond);
      azlist.set("Preconditioner Type","point relaxation");
    break;
    case INPAR::SOLVER::azprec_GaussSeidel:
      // using ifpack
      azlist.set("AZ_precond",AZ_user_precond);
      azlist.set("Preconditioner Type","point relaxation");
    break;
    case INPAR::SOLVER::azprec_DownwindGaussSeidel:
      // using ifpack
      azlist.set("AZ_precond",AZ_user_precond);
      azlist.set("Preconditioner Type","point relaxation");
      azlist.set<bool>("downwinding",true);
      azlist.set<double>("downwinding tau",inparams.get<double>("DWINDTAU"));
    break;
    case INPAR::SOLVER::azprec_LU:
      // using ifpack
      azlist.set("AZ_precond",AZ_user_precond);
      azlist.set("Preconditioner Type","Amesos");
    break;
    case INPAR::SOLVER::azprec_RILU:
      azlist.set("AZ_precond",AZ_dom_decomp);
      azlist.set("AZ_subdomain_solve",AZ_rilu);
      azlist.set("AZ_graph_fill",inparams.get<int>("IFPACKGFILL"));
    break;
    case INPAR::SOLVER::azprec_ICC:
      // using ifpack
      azlist.set("AZ_precond",AZ_user_precond);
      azlist.set("Preconditioner Type","IC");
    break;
    case INPAR::SOLVER::azprec_ML:
    case INPAR::SOLVER::azprec_MLfluid:
    case INPAR::SOLVER::azprec_MLAPI:
    case INPAR::SOLVER::azprec_MLfluid2:
    case INPAR::SOLVER::azprec_AMGBS:
    case INPAR::SOLVER::azprec_AMG:
    case INPAR::SOLVER::azprec_BGS2x2:
    case INPAR::SOLVER::azprec_BGSnxn:
    case INPAR::SOLVER::azprec_TekoSIMPLE:
      azlist.set("AZ_precond",AZ_user_precond);
    break;
    default:
      dserror("Unknown preconditioner for AztecOO");
    break;
    }
    //------------------------------------- set other aztec parameters
    azlist.set("AZ_kspace",inparams.get<int>("AZSUB"));
    azlist.set("AZ_max_iter",inparams.get<int>("AZITER"));
    azlist.set("AZ_overlap",inparams.get<int>("IFPACKOVERLAP"));
    azlist.set("AZ_type_overlap",AZ_symmetric);
    azlist.set("AZ_poly_ord",inparams.get<int>("AZPOLY"));
    const int azoutput = inparams.get<int>("AZOUTPUT");
    if (!azoutput)
      azlist.set("AZ_output",AZ_none);             // AZ_none AZ_all AZ_warnings AZ_last 10
    else
      azlist.set("AZ_output",azoutput);
    azlist.set("AZ_diagnostics",inparams.get<int>("AZBDIAG"));          // AZ_none AZ_all
    azlist.set("AZ_conv",DRT::INPUT::IntegralValue<int>(inparams,"AZCONV"));
    azlist.set("AZ_tol",inparams.get<double>("AZTOL"));
    azlist.set("AZ_drop",inparams.get<double>("AZDROP"));
    azlist.set("AZ_scaling",AZ_none);
    azlist.set("AZ_keep_info",0);
    // set reuse parameters
    azlist.set("ncall",0);                         // counting number of solver calls
    azlist.set("reuse",inparams.get<int>("AZREUSE"));            // reuse info for n solver calls
    //-------------------------------- set parameters for Ifpack if used
    if (azprectyp == INPAR::SOLVER::azprec_ILU  ||
        azprectyp == INPAR::SOLVER::azprec_ILUT ||
        azprectyp == INPAR::SOLVER::azprec_ICC  ||
        azprectyp == INPAR::SOLVER::azprec_LU   ||
        azprectyp == INPAR::SOLVER::azprec_SymmGaussSeidel ||
        azprectyp == INPAR::SOLVER::azprec_GaussSeidel ||
        azprectyp == INPAR::SOLVER::azprec_DownwindGaussSeidel ||
        azprectyp == INPAR::SOLVER::azprec_Jacobi)
    {
      ParameterList& ifpacklist = outparams.sublist("IFPACK Parameters");
      ifpacklist = TranslateBACIToIfpack(inparams);

      /*ifpacklist.set("relaxation: damping factor",inparams.get<double>("AZOMEGA"));
      ifpacklist.set("fact: drop tolerance",inparams.get<double>("AZDROP"));
      ifpacklist.set("fact: level-of-fill",inparams.get<int>("IFPACKGFILL"));
      ifpacklist.set("fact: ilut level-of-fill",inparams.get<double>("IFPACKFILL"));
      ifpacklist.set("partitioner: overlap",inparams.get<int>("IFPACKOVERLAP"));
      ifpacklist.set("schwarz: combine mode",inparams.get<string>("IFPACKCOMBINE")); // can be "Zero", "Add", "Insert"
      ifpacklist.set("schwarz: reordering type","rcm"); // "rcm" or "metis" or "amd"
      ifpacklist.set("amesos: solver type", "Amesos_Klu"); // can be "Amesos_Klu", "Amesos_Umfpack", "Amesos_Superlu"
      if (azprectyp == INPAR::SOLVER::azprec_SymmGaussSeidel)
      {
        ifpacklist.set("relaxation: type","symmetric Gauss-Seidel");
        ifpacklist.set("relaxation: sweeps",inparams.get<int>("IFPACKGFILL"));
        ifpacklist.set("relaxation: damping factor",inparams.get<double>("AZOMEGA"));
      }
      if (azprectyp == INPAR::SOLVER::azprec_GaussSeidel)
      {
        ifpacklist.set("relaxation: type","Gauss-Seidel");
        ifpacklist.set("relaxation: sweeps",inparams.get<int>("IFPACKGFILL"));
        ifpacklist.set("relaxation: damping factor",inparams.get<double>("AZOMEGA"));
      }
      if (azprectyp == INPAR::SOLVER::azprec_DownwindGaussSeidel)
      {
        // in case of downwinding prevent ifpack from again reordering
        ifpacklist.set("schwarz: reordering type","none");
        ifpacklist.set("relaxation: type","Gauss-Seidel");
        ifpacklist.set("relaxation: sweeps",inparams.get<int>("IFPACKGFILL"));
        ifpacklist.set("relaxation: damping factor",inparams.get<double>("AZOMEGA"));
      }
      if (azprectyp == INPAR::SOLVER::azprec_Jacobi)
      {
        ifpacklist.set("relaxation: type","Jacobi");
        ifpacklist.set("relaxation: sweeps",inparams.get<int>("IFPACKGFILL"));
        ifpacklist.set("relaxation: damping factor",inparams.get<double>("AZOMEGA"));
      }*/
    }
    //------------------------------------- set parameters for Teko if used
    if (azprectyp == INPAR::SOLVER::azprec_TekoSIMPLE ||
        azprectyp == INPAR::SOLVER::azprec_BGSnxn)
    {
      ParameterList& tekolist = outparams.sublist("Teko Parameters");
      tekolist = TranslateBACIToTeko(inparams);
      // TODO: Teko preconditioners -> Stratimikos interface? -> ask Eric about his plans
    }
    //------------------------------------- set parameters for ML if used
    if (azprectyp == INPAR::SOLVER::azprec_ML       ||
        azprectyp == INPAR::SOLVER::azprec_MLfluid  ||
        azprectyp == INPAR::SOLVER::azprec_MLfluid2 ||
        azprectyp == INPAR::SOLVER::azprec_MLAPI ||
        azprectyp == INPAR::SOLVER::azprec_AMG)          // new AMG preconditioner
    {
      ParameterList& mllist = outparams.sublist("ML Parameters");
      mllist = TranslateBACIToML(inparams,&azlist);

      /*ML_Epetra::SetDefaults("SA",mllist);
      switch (azprectyp)
      {
      case INPAR::SOLVER::azprec_ML: // do nothing, this is standard
      break;
      case INPAR::SOLVER::azprec_MLAPI: // set flag to use mlapi operator
        mllist.set<bool>("LINALG::AMG_Operator",true);
      break;
      case INPAR::SOLVER::azprec_MLfluid: // unsymmetric, unsmoothed restruction
        mllist.set("aggregation: use tentative restriction",true);
      break;
      case INPAR::SOLVER::azprec_MLfluid2: // full Pretrov-Galerkin unsymmetric smoothed
        mllist.set("energy minimization: enable",true);
        mllist.set("energy minimization: type",3); // 1,2,3 cheap -> expensive
        mllist.set("aggregation: block scaling",false);
      break;
      case INPAR::SOLVER::azprec_AMG:   // set flag to use AMGPreconditioner
        mllist.set<bool>("LINALG::AMGPreconditioner",true);
      break;
      default: dserror("Unknown type of ml preconditioner");
      }
      mllist.set("ML output"                       ,inparams.get<int>("ML_PRINT"));
      if (inparams.get<int>("ML_PRINT")==10)
        mllist.set("print unused"                  ,1);
      else
        mllist.set("print unused"                  ,-2);
      mllist.set("increasing or decreasing"        ,"increasing");
      mllist.set("coarse: max size"                ,inparams.get<int>("ML_MAXCOARSESIZE"));
      mllist.set("max levels"                      ,inparams.get<int>("ML_MAXLEVEL"));
      mllist.set("smoother: pre or post"           ,"both");
      mllist.set("aggregation: threshold"          ,inparams.get<double>("ML_PROLONG_THRES"));
      mllist.set("aggregation: damping factor"     ,inparams.get<double>("ML_PROLONG_SMO"));
      mllist.set("aggregation: nodes per aggregate",inparams.get<int>("ML_AGG_SIZE"));
      // override the default sweeps=2 with a default sweeps=1
      // individual level sweeps are set below
      mllist.set("smoother: sweeps",1);
      // save memory if this is an issue, make ML use single precision
      //mllist.set("low memory usage",true);
      switch (DRT::INPUT::IntegralValue<int>(inparams,"ML_COARSEN"))
      {
        case 0:  mllist.set("aggregation: type","Uncoupled");  break;
        case 1:  mllist.set("aggregation: type","METIS");      break;
        case 2:  mllist.set("aggregation: type","VBMETIS");    break;
        case 3:  mllist.set("aggregation: type","MIS");        break;
        default: dserror("Unknown type of coarsening for ML"); break;
      }

      // set ml smoothers
      const int mlmaxlevel = inparams.get<int>("ML_MAXLEVEL");
      // create vector of integers containing smoothing steps/polynomial order of level
      std::vector<int> mlsmotimessteps;
      {
        std::istringstream mlsmotimes(Teuchos::getNumericStringParameter(inparams,"ML_SMOTIMES"));
        std::string word;
        while (mlsmotimes >> word)
          mlsmotimessteps.push_back(std::atoi(word.c_str()));
      }

      if ((int)mlsmotimessteps.size() < mlmaxlevel)
        dserror("Not enough smoothing steps ML_SMOTIMES=%d, must be larger/equal than ML_MAXLEVEL=%d\n",
                mlsmotimessteps.size(),mlmaxlevel);

      for (int i=0; i<mlmaxlevel-1; ++i)
      {
        char levelstr[11];
        sprintf(levelstr,"(level %d)",i);
        ParameterList& smolevelsublist = mllist.sublist("smoother: list "+(string)levelstr);
        int type;
        double damp;
        if (i==0)
        {
          type = DRT::INPUT::IntegralValue<int>(inparams,"ML_SMOOTHERFINE");
          damp = inparams.get<double>("ML_DAMPFINE");
        }
        else if (i < mlmaxlevel-1)
        {
          type = DRT::INPUT::IntegralValue<int>(inparams,"ML_SMOOTHERMED");
          damp = inparams.get<double>("ML_DAMPMED");
        }
        else
        {
          type = DRT::INPUT::IntegralValue<int>(inparams,"ML_SMOOTHERCOARSE");
          damp = inparams.get<double>("ML_DAMPCOARSE");
        }
        switch (type)
        {
        case 0: // SGS
          smolevelsublist.set("smoother: type"                        ,"symmetric Gauss-Seidel");
          smolevelsublist.set("smoother: sweeps"                      ,mlsmotimessteps[i]);
          smolevelsublist.set("smoother: damping factor"              ,damp);
        break;
        case 7: // GS
          smolevelsublist.set("smoother: type"                        ,"Gauss-Seidel");
          smolevelsublist.set("smoother: sweeps"                      ,mlsmotimessteps[i]);
          smolevelsublist.set("smoother: damping factor"              ,damp);
        break;
        case 8: // DGS
          smolevelsublist.set("smoother: type"                        ,"Gauss-Seidel");
          smolevelsublist.set("smoother: sweeps"                      ,mlsmotimessteps[i]);
          smolevelsublist.set("smoother: damping factor"              ,damp);
          azlist.set<bool>("downwinding",true);
          azlist.set<double>("downwinding tau",inparams.get<double>("DWINDTAU"));
          {
            ParameterList& ifpacklist = mllist.sublist("smoother: ifpack list");
            ifpacklist.set("schwarz: reordering type","true");
          }
        break;
        case 1: // Jacobi
          smolevelsublist.set("smoother: type"                        ,"Jacobi");
          smolevelsublist.set("smoother: sweeps"                      ,mlsmotimessteps[i]);
          smolevelsublist.set("smoother: damping factor"              ,damp);
        break;
        case 2: // Chebychev
          smolevelsublist.set("smoother: type"                        ,"MLS");
          smolevelsublist.set("smoother: MLS polynomial order"        ,mlsmotimessteps[i]);
        break;
        case 3: // MLS
          smolevelsublist.set("smoother: type"                        ,"MLS");
          smolevelsublist.set("smoother: MLS polynomial order"        ,-mlsmotimessteps[i]);
        break;
        case 4: // Ifpack's ILU
        {
          smolevelsublist.set("smoother: type"                        ,"IFPACK");
          smolevelsublist.set("smoother: ifpack type"                 ,"ILU");
          smolevelsublist.set("smoother: ifpack overlap"              ,inparams.get<int>("IFPACKOVERLAP"));
          smolevelsublist.set<double>("smoother: ifpack level-of-fill",(double)mlsmotimessteps[i]);
          ParameterList& ifpacklist = mllist.sublist("smoother: ifpack list");
          ifpacklist.set("schwarz: reordering type","rcm"); // "rcm" or "metis" or "amd" or "true"
          ifpacklist.set("schwarz: combine mode",inparams.get<string>("IFPACKCOMBINE")); // can be "Zero", "Insert", "Add"
          ifpacklist.set("partitioner: overlap",inparams.get<int>("IFPACKOVERLAP"));
        }
        break;
        case 5: // Amesos' KLU
          smolevelsublist.set("smoother: type"                        ,"Amesos-KLU");
        break;
        case 9: // Amesos' Umfpack
          smolevelsublist.set("smoother: type"                        ,"Amesos-UMFPACK");
        break;
#ifdef PARALLEL
        case 6: // Amesos' SuperLU_Dist
          smolevelsublist.set("smoother: type"                        ,"Amesos-Superludist");
        break;
#endif
        default: dserror("Unknown type of smoother for ML: tuple %d",type); break;
        } // switch (type)
      } // for (int i=0; i<azvar->mlmaxlevel-1; ++i)

      // set coarse grid solver
      const int coarse = mlmaxlevel-1;
      switch (DRT::INPUT::IntegralValue<int>(inparams,"ML_SMOOTHERCOARSE"))
      {
        case 0:
          mllist.set("coarse: type"          ,"symmetric Gauss-Seidel");
          mllist.set("coarse: sweeps"        ,mlsmotimessteps[coarse]);
          mllist.set("coarse: damping factor",inparams.get<double>("ML_DAMPCOARSE"));
        break;
        case 7:
          mllist.set("coarse: type"          ,"Gauss-Seidel");
          mllist.set("coarse: sweeps"        ,mlsmotimessteps[coarse]);
          mllist.set("coarse: damping factor",inparams.get<double>("ML_DAMPCOARSE"));
        break;
        case 8:
          mllist.set("coarse: type"          ,"Gauss-Seidel");
          mllist.set("coarse: sweeps"        ,mlsmotimessteps[coarse]);
          mllist.set("coarse: damping factor",inparams.get<double>("ML_DAMPCOARSE"));
          azlist.set<bool>("downwinding",true);
          azlist.set<double>("downwinding tau",inparams.get<double>("DWINDTAU"));
          {
            ParameterList& ifpacklist = mllist.sublist("smoother: ifpack list");
            ifpacklist.set("schwarz: reordering type","true");
          }
        break;
        case 1:
          mllist.set("coarse: type"          ,"Jacobi");
          mllist.set("coarse: sweeps"        ,mlsmotimessteps[coarse]);
          mllist.set("coarse: damping factor",inparams.get<double>("ML_DAMPCOARSE"));
        break;
        case 2:
          mllist.set("coarse: type"                ,"MLS");
          mllist.set("coarse: MLS polynomial order",mlsmotimessteps[coarse]);
        break;
        case 3:
          mllist.set("coarse: type"                ,"MLS");
          mllist.set("coarse: MLS polynomial order",-mlsmotimessteps[coarse]);
        break;
        case 4:
        {
          mllist.set("coarse: type"          ,"IFPACK");
          mllist.set("coarse: ifpack type"   ,"ILU");
          mllist.set("coarse: ifpack overlap",0);
          mllist.set<double>("coarse: ifpack level-of-fill",(double)mlsmotimessteps[coarse]);
          ParameterList& ifpacklist = mllist.sublist("coarse: ifpack list");
          ifpacklist.set<int>("fact: level-of-fill",mlsmotimessteps[coarse]);
          ifpacklist.set("schwarz: reordering type","rcm");
          ifpacklist.set("schwarz: combine mode",inparams.get<string>("IFPACKCOMBINE")); // can be "Zero", "Insert", "Add"
          ifpacklist.set("partitioner: overlap",inparams.get<int>("IFPACKOVERLAP"));
        }
        break;
        case 5:
          mllist.set("coarse: type","Amesos-KLU");
        break;
        case 9:
          mllist.set("coarse: type","Amesos-UMFPACK");
        break;
        case 6:
          mllist.set("coarse: type","Amesos-Superludist");
        break;
        default: dserror("Unknown type of coarse solver for ML"); break;
      } // switch (azvar->mlsmotype_coarse)
      // default values for nullspace
      mllist.set("PDE equations",1);
      mllist.set("null space: dimension",1);
      mllist.set("null space: type","pre-computed");
      mllist.set("null space: add default vectors",false);
      mllist.set<double*>("null space: vectors",NULL);
#if defined(PARALLEL) && defined(PARMETIS)
      mllist.set("repartition: enable",1);
      mllist.set("repartition: partitioner","ParMETIS");
      mllist.set("repartition: max min ratio",1.3);
      mllist.set("repartition: min per proc",3000);
#endif
      //cout << mllist << endl << endl << endl; fflush(stdout);*/
    } // if ml preconditioner
    //------------------------------------- set parameters for AMGBS if used
    if (azprectyp == INPAR::SOLVER::azprec_AMGBS      )
    {
      ParameterList& amglist = outparams.sublist("AMGBS Parameters");
      ML_Epetra::SetDefaults("SA",amglist);
      amglist.set("amgbs: smoother: pre or post"    ,"both");
      amglist.set("amgbs: prolongator smoother (vel)",inparams.get<string>("AMGBS_PSMOOTHER_VEL"));
      amglist.set("amgbs: prolongator smoother (pre)",inparams.get<string>("AMGBS_PSMOOTHER_PRE"));

      amglist.set("output"                          ,inparams.get<int>("ML_PRINT"));
      amglist.set("coarse: max size"                ,inparams.get<int>("ML_MAXCOARSESIZE"));
      amglist.set("max levels"                      ,inparams.get<int>("ML_MAXLEVEL"));
      amglist.set("aggregation: threshold"          ,inparams.get<double>("ML_PROLONG_THRES"));
      amglist.set("aggregation: damping factor"     ,inparams.get<double>("ML_PROLONG_SMO"));
      amglist.set("aggregation: nodes per aggregate",inparams.get<int>("ML_AGG_SIZE"));
      // override the default sweeps=2 with a default sweeps=1
      // individual level sweeps are set below
      amglist.set("smoother: sweeps",1);
      switch (DRT::INPUT::IntegralValue<int>(inparams,"ML_COARSEN"))
      {
        case 0:  amglist.set("aggregation: type","Uncoupled");  break;
        case 1:  amglist.set("aggregation: type","METIS");      break;
        case 2:  amglist.set("aggregation: type","VBMETIS");    break;
        case 3:  amglist.set("aggregation: type","MIS");        break;
        default: dserror("Unknown type of coarsening for ML"); break;
      }

      //////////////////// set braess-sarazin smoothers
      const int mlmaxlevel = inparams.get<int>("ML_MAXLEVEL");
      // create vector of integers containing smoothing steps with braess-sarazin
      std::vector<int> mlsmotimessteps;
      {
        std::istringstream mlsmotimes(Teuchos::getNumericStringParameter(inparams,"ML_SMOTIMES"));
        std::string word;
        while (mlsmotimes >> word)
          mlsmotimessteps.push_back(std::atoi(word.c_str()));
      }

      if ((int)mlsmotimessteps.size() < mlmaxlevel)
        dserror("Not enough smoothing steps ML_SMOTIMES=%d, must be larger/equal than ML_MAXLEVEL=%d\n",
                mlsmotimessteps.size(),mlmaxlevel);

      //////////////////// read in damping parameters for Braess Sarazin
      std::vector<double> bsdamping;   // damping parameters for Braess-Sarazin
      {
        double word;
        std::istringstream bsdampingstream(Teuchos::getNumericStringParameter(inparams,"AMGBS_BS_DAMPING"));
        while (bsdampingstream >> word)
          bsdamping.push_back(word);
      }
      if ((int)bsdamping.size() < mlmaxlevel)
        dserror("Not enough damping factors AMGBS_BS_DAMPING=%d, must be larger/equal than ML_MAXLEVEL=%d\n",
                bsdamping.size(),mlmaxlevel);

      //////////////////// read in smoothing time steps for pressure correction equation with jacobi/sgs iteration
      std::vector<int> bspcsweeps;
      {
        std::istringstream bsdampingstream(Teuchos::getNumericStringParameter(inparams,"AMGBS_BS_PCSWEEPS"));
        std::string word;
        while (bsdampingstream >> word)
          bspcsweeps.push_back(std::atoi(word.c_str()));
      }
      if ((int)bspcsweeps.size() < mlmaxlevel)
        dserror("Not enough integers for number of smoothing iterations of pressure correction equation within Braess-Sarazin. AMGBS_BS_PCSWEEPS=%d, must be larger/equal than ML_MAXLEVEL=%d\n",
            bspcsweeps.size(),mlmaxlevel);

      //////////////////// read in damping parameters for jacobi/sgs iteration in pressure correction solver
      std::vector<double> bspcdamping;
      {
        double word;
        std::istringstream bspcdampingstream(Teuchos::getNumericStringParameter(inparams,"AMGBS_BS_PCDAMPING"));
        while (bspcdampingstream >> word)
          bspcdamping.push_back(word);
      }
      if ((int)bspcdamping.size() < mlmaxlevel)
        dserror("Not enough damping factors AMGBS_BS_PCDAMPING=%d, must be larger/equal than ML_MAXLEVEL=%d\n",
                bspcdamping.size(),mlmaxlevel);

      // Parameters for presssure correction smoother
      for (int i=0; i<mlmaxlevel; ++i)
      {
        char levelstr[11];
        sprintf(levelstr,"(level %d)",i);
        ParameterList& smolevelsublist = amglist.sublist("braess-sarazin: list "+(string)levelstr);

          smolevelsublist.set("braess-sarazin: sweeps"                      ,mlsmotimessteps[i]);
          smolevelsublist.set("braess-sarazin: damping factor"              ,bsdamping[i]);
          smolevelsublist.set("pressure correction approx: type"      ,"ILU");   // TODO choose Umfpack, KLU, ILU, Jacobi, Gauss-Seidel, symmetric Gauss-Seidel
          smolevelsublist.set("Ifpack overlap",inparams.get<int>("IFPACKOVERLAP"));

          switch (DRT::INPUT::IntegralValue<int>(inparams,"AMGBS_BS_PCCOARSE"))
          {
            case 0:
            {
              smolevelsublist.set("coarse: type","Umfpack");
            }
            break;
            case 1:
            {
              smolevelsublist.set("coarse: type","KLU");
            }
            break;
            case 2: // ILU is chosen
            {
              smolevelsublist.set("coarse: type"          ,"ILU");
              ParameterList& ifpacklist = smolevelsublist.sublist("IFPACK Parameters coarse");
              ifpacklist.set("relaxation: damping factor",inparams.get<double>("AZOMEGA"));
              ifpacklist.set("fact: drop tolerance",inparams.get<double>("AZDROP"));
              ifpacklist.set("fact: level-of-fill",inparams.get<int>("IFPACKGFILL"));
              ifpacklist.set("fact: ilut level-of-fill",inparams.get<double>("IFPACKFILL"));
              ifpacklist.set("schwarz: combine mode",inparams.get<string>("IFPACKCOMBINE")); // can be "Add", "Zero", "Insert", "InsertAdd", "Average", "AbsMax"
              ifpacklist.set("schwarz: reordering type","rcm"); // "rcm" or "metis"
              ifpacklist.set("partitioner: overlap",inparams.get<int>("IFPACKOVERLAP"));
            }
            break;
            case 3:
            case 6:
            {
              smolevelsublist.set("coarse: type","Jacobi");
              ParameterList& ifpacklist = smolevelsublist.sublist("IFPACK Parameters coarse");
              ifpacklist.set("fact: drop tolerance",inparams.get<double>("AZDROP"));
              ifpacklist.set("partitioner: overlap",inparams.get<int>("IFPACKOVERLAP"));
              ifpacklist.set("schwarz: combine mode",inparams.get<string>("IFPACKCOMBINE")); // can be "Add", "Zero", "Insert", "InsertAdd", "Average", "AbsMax"
              ifpacklist.set("schwarz: reordering type","rcm"); // "rcm" or "metis"
              ifpacklist.set("relaxation: type","Jacobi");
              ifpacklist.set("relaxation: sweeps",bspcsweeps[i]);
              ifpacklist.set("relaxation: damping factor",bspcdamping[i]);
              break;
            }
            case 4:
            case 7:
            {
              smolevelsublist.set("coarse: type","Gauss-Seidel");
              ParameterList& ifpacklist = smolevelsublist.sublist("IFPACK Parameters coarse");
              ifpacklist.set("fact: drop tolerance",inparams.get<double>("AZDROP"));
              ifpacklist.set("partitioner: overlap",inparams.get<int>("IFPACKOVERLAP"));
              ifpacklist.set("schwarz: combine mode",inparams.get<string>("IFPACKCOMBINE")); // can be "Add", "Zero", "Insert", "InsertAdd", "Average", "AbsMax"
              ifpacklist.set("schwarz: reordering type","rcm"); // "rcm" or "metis"
              ifpacklist.set("relaxation: type","Gauss-Seidel");
              ifpacklist.set("relaxation: sweeps",bspcsweeps[i]);
              ifpacklist.set("relaxation: damping factor",bspcdamping[i]);
              break;
            }
            case 5:
            case 8:
            {
              smolevelsublist.set("coarse: type","symmetric Gauss-Seidel");
              ParameterList& ifpacklist = smolevelsublist.sublist("IFPACK Parameters coarse");
              ifpacklist.set("fact: drop tolerance",inparams.get<double>("AZDROP"));
              ifpacklist.set("partitioner: overlap",inparams.get<int>("IFPACKOVERLAP"));
              ifpacklist.set("schwarz: combine mode",inparams.get<string>("IFPACKCOMBINE")); // can be "Add", "Zero", "Insert", "InsertAdd", "Average", "AbsMax"
              ifpacklist.set("schwarz: reordering type","rcm"); // "rcm" or "metis"
              ifpacklist.set("relaxation: type","symmetric Gauss-Seidel");
              ifpacklist.set("relaxation: sweeps",bspcsweeps[i]);
              ifpacklist.set("relaxation: damping factor",bspcdamping[i]);
              break;
            }
            default: dserror("Unknown type of coarse solver for pressure correction equation"); break;
          } // switch (azvar->mlsmotype_coarse)

          switch (DRT::INPUT::IntegralValue<int>(inparams,"AMGBS_BS_PCMEDIUM"))
          {
            case 0:
            {
              smolevelsublist.set("medium: type","Umfpack");
            }
            break;
            case 1:
            {
              smolevelsublist.set("medium: type","KLU");
            }
            break;
            case 2:
            {
              smolevelsublist.set("medium: type"          ,"ILU");
              ParameterList& ifpacklist = smolevelsublist.sublist("IFPACK Parameters medium");
              ifpacklist.set("relaxation: damping factor",inparams.get<double>("AZOMEGA"));
              ifpacklist.set("fact: drop tolerance",inparams.get<double>("AZDROP"));
              ifpacklist.set("fact: level-of-fill",inparams.get<int>("IFPACKGFILL"));
              ifpacklist.set("fact: ilut level-of-fill",inparams.get<double>("IFPACKFILL"));
              ifpacklist.set("schwarz: combine mode",inparams.get<string>("IFPACKCOMBINE")); // can be "Add", "Zero", "Insert", "InsertAdd", "Average", "AbsMax"
              ifpacklist.set("schwarz: reordering type","rcm"); // "rcm" or "metis"
              ifpacklist.set("partitioner: overlap",inparams.get<int>("IFPACKOVERLAP"));
            }
            break;
            case 3:
            case 6:
            {
              smolevelsublist.set("medium: type","Jacobi");
              ParameterList& ifpacklist = smolevelsublist.sublist("IFPACK Parameters medium");
              ifpacklist.set("fact: drop tolerance",inparams.get<double>("AZDROP"));
              ifpacklist.set("partitioner: overlap",inparams.get<int>("IFPACKOVERLAP"));
              ifpacklist.set("schwarz: combine mode",inparams.get<string>("IFPACKCOMBINE")); // can be "Add", "Zero", "Insert", "InsertAdd", "Average", "AbsMax"
              ifpacklist.set("schwarz: reordering type","rcm"); // "rcm" or "metis"
              ifpacklist.set("relaxation: type","Jacobi");
              ifpacklist.set("relaxation: sweeps",bspcsweeps[i]);
              ifpacklist.set("relaxation: damping factor",bspcdamping[i]);
              break;
            }
            case 4:
            case 7:
            {
              smolevelsublist.set("medium: type","Gauss-Seidel");
              ParameterList& ifpacklist = smolevelsublist.sublist("IFPACK Parameters medium");
              ifpacklist.set("fact: drop tolerance",inparams.get<double>("AZDROP"));
              ifpacklist.set("partitioner: overlap",inparams.get<int>("IFPACKOVERLAP"));
              ifpacklist.set("schwarz: combine mode",inparams.get<string>("IFPACKCOMBINE")); // can be "Add", "Zero", "Insert", "InsertAdd", "Average", "AbsMax"
              ifpacklist.set("schwarz: reordering type","rcm"); // "rcm" or "metis"
              ifpacklist.set("relaxation: type","Gauss-Seidel");
              ifpacklist.set("relaxation: sweeps",bspcsweeps[i]);
              ifpacklist.set("relaxation: damping factor",bspcdamping[i]);
              break;
            }
            case 5:
            case 8:
            {
              smolevelsublist.set("medium: type","symmetric Gauss-Seidel");
              ParameterList& ifpacklist = smolevelsublist.sublist("IFPACK Parameters medium");
              ifpacklist.set("fact: drop tolerance",inparams.get<double>("AZDROP"));
              ifpacklist.set("partitioner: overlap",inparams.get<int>("IFPACKOVERLAP"));
              ifpacklist.set("schwarz: combine mode",inparams.get<string>("IFPACKCOMBINE")); // can be "Add", "Zero", "Insert", "InsertAdd", "Average", "AbsMax"
              ifpacklist.set("schwarz: reordering type","rcm"); // "rcm" or "metis"
              ifpacklist.set("relaxation: type","symmetric Gauss-Seidel");
              ifpacklist.set("relaxation: sweeps",bspcsweeps[i]);
              ifpacklist.set("relaxation: damping factor",bspcdamping[i]);
              break;
            }
            default: dserror("Unknown type of medium level solver for pressure correction equation"); break;
          }

          switch (DRT::INPUT::IntegralValue<int>(inparams,"AMGBS_BS_PCFINE"))
          {
            case 0:
            {
              smolevelsublist.set("fine: type","Umfpack");
            }
            break;
            case 1:
            {
              smolevelsublist.set("fine: type","KLU");
            }
            break;
            case 2:
            {
              smolevelsublist.set("fine: type"          ,"ILU");
              ParameterList& ifpacklist = smolevelsublist.sublist("IFPACK Parameters fine");
              ifpacklist.set("relaxation: damping factor",inparams.get<double>("AZOMEGA"));
              ifpacklist.set("fact: drop tolerance",inparams.get<double>("AZDROP"));
              ifpacklist.set("fact: level-of-fill",inparams.get<int>("IFPACKGFILL"));
              ifpacklist.set("fact: ilut level-of-fill",inparams.get<double>("IFPACKFILL"));
              ifpacklist.set("schwarz: combine mode",inparams.get<string>("IFPACKCOMBINE")); // can be "Add", "Zero", "Insert", "InsertAdd", "Average", "AbsMax"
              ifpacklist.set("schwarz: reordering type","rcm"); // "rcm" or "metis"
              ifpacklist.set("partitioner: overlap",inparams.get<int>("IFPACKOVERLAP"));
            }
            break;
            case 3:
            case 6:
            {
              smolevelsublist.set("fine: type","Jacobi");
              ParameterList& ifpacklist = smolevelsublist.sublist("IFPACK Parameters fine");
              ifpacklist.set("fact: drop tolerance",inparams.get<double>("AZDROP"));
              ifpacklist.set("partitioner: overlap",inparams.get<int>("IFPACKOVERLAP"));
              ifpacklist.set("schwarz: combine mode",inparams.get<string>("IFPACKCOMBINE")); // can be "Add", "Zero", "Insert", "InsertAdd", "Average", "AbsMax"
              ifpacklist.set("schwarz: reordering type","rcm"); // "rcm" or "metis"
              ifpacklist.set("relaxation: type","Jacobi");
              ifpacklist.set("relaxation: sweeps",bspcsweeps[i]);
              ifpacklist.set("relaxation: damping factor",bspcdamping[i]);
              break;
            }
            case 4:
            case 7:
            {
              smolevelsublist.set("fine: type","Gauss-Seidel");
              ParameterList& ifpacklist = smolevelsublist.sublist("IFPACK Parameters fine");
              ifpacklist.set("fact: drop tolerance",inparams.get<double>("AZDROP"));
              ifpacklist.set("partitioner: overlap",inparams.get<int>("IFPACKOVERLAP"));
              ifpacklist.set("schwarz: combine mode",inparams.get<string>("IFPACKCOMBINE")); // can be "Add", "Zero", "Insert", "InsertAdd", "Average", "AbsMax"
              ifpacklist.set("schwarz: reordering type","rcm"); // "rcm" or "metis"
              ifpacklist.set("relaxation: type","Gauss-Seidel");
              ifpacklist.set("relaxation: sweeps",bspcsweeps[i]);
              ifpacklist.set("relaxation: damping factor",bspcdamping[i]);
              break;
            }
            case 5:
            case 8:
            {
              smolevelsublist.set("fine: type","symmetric Gauss-Seidel");
              ParameterList& ifpacklist = smolevelsublist.sublist("IFPACK Parameters fine");
              ifpacklist.set("fact: drop tolerance",inparams.get<double>("AZDROP"));
              ifpacklist.set("partitioner: overlap",inparams.get<int>("IFPACKOVERLAP"));
              ifpacklist.set("schwarz: combine mode",inparams.get<string>("IFPACKCOMBINE")); // can be "Add", "Zero", "Insert", "InsertAdd", "Average", "AbsMax"
              ifpacklist.set("schwarz: reordering type","rcm"); // "rcm" or "metis"
              ifpacklist.set("relaxation: type","symmetric Gauss-Seidel");
              ifpacklist.set("relaxation: sweeps",bspcsweeps[i]);
              ifpacklist.set("relaxation: damping factor",bspcdamping[i]);
              break;
            }
            default: dserror("Unknown type of coarse solver for pressure correction equation"); break;
          } // switch (azvar->mlsmotype_coarse)


      } // for (int i=0; i<azvar->mlmaxlevel-1; ++i)



      amglist.set("PDE equations",1);
      amglist.set("null space: dimension",1);
      amglist.set("null space: type","pre-computed");
      amglist.set("null space: add default vectors",false);
      amglist.set<double*>("null space: vectors",NULL);



      cout << amglist << endl; // TODO delete me
    } // if AMGBS preconditioner
    if (azprectyp == INPAR::SOLVER::azprec_BGS2x2)
    {
      ParameterList& bgslist = outparams.sublist("BGS Parameters");
      bgslist.set("numblocks",2);

      // currently, the number of Gauss-Seidel iterations and the relaxation
      // parameter on the global level are set to 1 and 1.0, respectively
      bgslist.set("global_iter",1);
      bgslist.set("global_omega",1.0);

      // currently, the order of blocks in the given EpetraOperator is not changed
      // in the Gauss-Seidel procedure
      bgslist.set("fliporder",false);

      // currently, the number of Richardson iteratios and the relaxation
      // parameter on the individual block level are set to 1 and 1.0, respectively
      bgslist.set("block1_iter",1);
      bgslist.set("block1_omega",1.0);
      bgslist.set("block2_iter",1);
      bgslist.set("block2_omega",1.0);
    }
  }
  break;
  default:
    dserror("Unsupported type of solver");
  break;
  }

  //================================================================== deliver
  return outparams;
}


/*----------------------------------------------------------------------*
 | Multiply matrices A*B                                     mwgee 02/08|
 *----------------------------------------------------------------------*/
RCP<LINALG::SparseMatrix> LINALG::MLMultiply(const LINALG::SparseMatrix& A,
                                             const LINALG::SparseMatrix& B,
                                             bool complete)
{
  return MLMultiply(*A.EpetraMatrix(),*B.EpetraMatrix(),
                    A.explicitdirichlet_,A.savegraph_,complete);
}

/*----------------------------------------------------------------------*
 | Multiply matrices A*B                                     mwgee 02/08|
 *----------------------------------------------------------------------*/
RCP<LINALG::SparseMatrix> LINALG::MLMultiply(const LINALG::SparseMatrix& A,
                                             const LINALG::SparseMatrix& B,
                                             bool explicitdirichlet,
                                             bool savegraph,
                                             bool complete)
{
  return MLMultiply(*A.EpetraMatrix(),*B.EpetraMatrix(),
                    explicitdirichlet,savegraph,complete);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<LINALG::SparseMatrix> LINALG::MLMultiply(const LINALG::SparseMatrix& A,
                                                      bool transA,
                                                      const LINALG::SparseMatrix& B,
                                                      bool transB,
                                                      bool explicitdirichlet,
                                                      bool savegraph,
                                                      bool completeoutput)
{
  // make sure FillComplete was called on the matrices
  if (!A.Filled()) dserror("A has to be FillComplete");
  if (!B.Filled()) dserror("B has to be FillComplete");

  EpetraExt::RowMatrix_Transpose transposera(true,NULL,false);
  EpetraExt::RowMatrix_Transpose transposerb(true,NULL,false);
  Epetra_CrsMatrix* Atrans = NULL;
  Epetra_CrsMatrix* Btrans = NULL;
  if (transA)
    Atrans = &(dynamic_cast<Epetra_CrsMatrix&>(transposera(*A.EpetraMatrix())));
  else
    Atrans = A.EpetraMatrix().get();
  if (transB)
    Btrans = &(dynamic_cast<Epetra_CrsMatrix&>(transposerb(*B.EpetraMatrix())));
  else
    Btrans = B.EpetraMatrix().get();

  Teuchos::RCP<LINALG::SparseMatrix> C;
  C = LINALG::MLMultiply(*Atrans,*Btrans,explicitdirichlet,savegraph,completeoutput);

  return C;
}

/*----------------------------------------------------------------------*
 | Multiply matrices A*B                                     mwgee 02/08|
 *----------------------------------------------------------------------*/
//static void CopySortDeleteZeros(const Epetra_CrsMatrix& A, Epetra_CrsMatrix& As);
RCP<LINALG::SparseMatrix> LINALG::MLMultiply(const Epetra_CrsMatrix& Aorig,
                                             const Epetra_CrsMatrix& Borig,
                                             bool explicitdirichlet,
                                             bool savegraph,
                                             bool complete)
{
  EpetraExt::CrsMatrix_SolverMap Atransform;
  EpetraExt::CrsMatrix_SolverMap Btransform;
  const Epetra_CrsMatrix& A = Atransform(const_cast<Epetra_CrsMatrix&>(Aorig));
  const Epetra_CrsMatrix& B = Btransform(const_cast<Epetra_CrsMatrix&>(Borig));

  // make sure FillComplete was called on the matrices
  if (!A.Filled()) dserror("A has to be FillComplete");
  if (!B.Filled()) dserror("B has to be FillComplete");

  // For debugging, it might be helpful when all columns are
  // sorted and all zero values are wiped from the input:
  //RCP<Epetra_CrsMatrix> As = CreateMatrix(A.RowMap(),A.MaxNumEntries());
  //RCP<Epetra_CrsMatrix> Bs = CreateMatrix(B.RowMap(),B.MaxNumEntries());
  //CopySortDeleteZeros(A,*As);
  //CopySortDeleteZeros(B,*Bs);
  ML_Operator* ml_As = ML_Operator_Create(MLAPI::GetML_Comm());
  ML_Operator* ml_Bs = ML_Operator_Create(MLAPI::GetML_Comm());
  //ML_Operator_WrapEpetraMatrix(As.get(),ml_As);
  //ML_Operator_WrapEpetraMatrix(Bs.get(),ml_Bs);
  ML_Operator_WrapEpetraMatrix(const_cast<Epetra_CrsMatrix*>(&A),ml_As);
  ML_Operator_WrapEpetraMatrix(const_cast<Epetra_CrsMatrix*>(&B),ml_Bs);
  ML_Operator* ml_AtimesB = ML_Operator_Create(MLAPI::GetML_Comm());
  ML_2matmult(ml_As,ml_Bs,ml_AtimesB,ML_CSR_MATRIX); // do NOT use ML_EpetraCRS_MATRIX !!
  ML_Operator_Destroy(&ml_As);
  ML_Operator_Destroy(&ml_Bs);
  // For ml_AtimesB we have to reconstruct the column map in global indexing,
  // The following is going down to the salt-mines of ML ...
  int N_local = ml_AtimesB->invec_leng;
  ML_CommInfoOP* getrow_comm = ml_AtimesB->getrow->pre_comm;
  if (!getrow_comm) dserror("ML_Operator does not have CommInfo");
  ML_Comm* comm = ml_AtimesB->comm;
  if (N_local != B.DomainMap().NumMyElements())
    dserror("Mismatch in local row dimension between ML and Epetra");
  int N_rcvd  = 0;
  int N_send  = 0;
  int flag    = 0;
  for (int i=0; i<getrow_comm->N_neighbors; i++)
  {
    N_rcvd += (getrow_comm->neighbors)[i].N_rcv;
    N_send += (getrow_comm->neighbors)[i].N_send;
    if (  ((getrow_comm->neighbors)[i].N_rcv != 0) &&
       ((getrow_comm->neighbors)[i].rcv_list != NULL) )  flag = 1;
  }
  // For some unknown reason, ML likes to have stuff one larger than
  // neccessary...
  vector<double> dtemp(N_local+N_rcvd+1);
  vector<int>    cmap(N_local+N_rcvd+1);
  for (int i=0; i<N_local; ++i)
  {
    cmap[i] = B.DomainMap().GID(i);
    dtemp[i] = (double)cmap[i];
  }
  ML_cheap_exchange_bdry(&dtemp[0],getrow_comm,N_local,N_send,comm);
  if (flag)
  {
    int count = N_local;
    const int neighbors = getrow_comm->N_neighbors;
    for (int i=0; i<neighbors; i++)
    {
      const int nrcv = getrow_comm->neighbors[i].N_rcv;
      for (int j=0; j<nrcv; j++)
        cmap[getrow_comm->neighbors[i].rcv_list[j]] = (int)dtemp[count++];
    }
  }
  else
    for (int i=0; i<N_local+N_rcvd; ++i) cmap[i] = (int)dtemp[i];
  dtemp.clear();

  // we can now determine a matching column map for the result
  Epetra_Map gcmap(-1,N_local+N_rcvd,&cmap[0],0,A.Comm());

  int allocated=0;
  int rowlength;
  double* val=NULL;
  int* bindx=NULL;
  const int myrowlength = A.RowMap().NumMyElements();
  const Epetra_Map& rowmap = A.RowMap();

  // determine the maximum bandwith for the result matrix.
  // replaces the old, very(!) memory-consuming guess:
  // int guessnpr = A.MaxNumEntries()*B.MaxNumEntries();
  int educatedguess = 0;
  for (int i=0; i<myrowlength; ++i)
  {
    // get local row
    ML_get_matrix_row(ml_AtimesB,1,&i,&allocated,&bindx,&val,&rowlength,0);
    if (rowlength>educatedguess) educatedguess = rowlength;
  }

  // allocate our result matrix and fill it
  RCP<Epetra_CrsMatrix> result
    = rcp(new Epetra_CrsMatrix(::Copy,A.RangeMap(),gcmap,educatedguess,false));

  vector<int> gcid(educatedguess);
  for (int i=0; i<myrowlength; ++i)
  {
    const int grid = rowmap.GID(i);
    // get local row
    ML_get_matrix_row(ml_AtimesB,1,&i,&allocated,&bindx,&val,&rowlength,0);
    if (!rowlength) continue;
    if ((int)gcid.size() < rowlength) gcid.resize(rowlength);
    for (int j=0; j<rowlength; ++j)
    {
      gcid[j] = gcmap.GID(bindx[j]);
#ifdef DEBUG
      if (gcid[j]<0) dserror("This is really bad... cannot find gcid");
#endif
    }
#ifdef DEBUG
    int err = result->InsertGlobalValues(grid,rowlength,val,&gcid[0]);
    if (err!=0 && err!=1) dserror("Epetra_CrsMatrix::InsertGlobalValues returned err=%d",err);
#else
    result->InsertGlobalValues(grid,rowlength,val,&gcid[0]);
#endif
  }
  if (bindx) ML_free(bindx);
  if (val) ML_free(val);
  ML_Operator_Destroy(&ml_AtimesB);
  if (complete)
  {
    int err = result->FillComplete(B.DomainMap(),A.RangeMap());
    if (err) dserror("Epetra_CrsMatrix::FillComplete returned err=%d",err);

#if 0 // the current status is that we don't need this (mwgee)
    EpetraExt::CrsMatrix_SolverMap ABtransform;
    const Epetra_CrsMatrix& tmp = ABtransform(*result);
    RCP<Epetra_CrsMatrix> finalresult = rcp(new Epetra_CrsMatrix(*result));
    if (!finalresult->Filled())
    {
      finalresult->FillComplete(B.DomainMap(),A.RangeMap());
      finalresult->OptimizeStorage();
    }
    result = null;
    return rcp(new SparseMatrix(finalresult,explicitdirichlet,savegraph));
#endif
  }
  return rcp(new SparseMatrix(result,explicitdirichlet,savegraph));
}
/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
/*
  static void CopySortDeleteZeros(const Epetra_CrsMatrix& A, Epetra_CrsMatrix& As)
  {
  vector<int>    scindices(A.MaxNumEntries());
  vector<double> scvalues(A.MaxNumEntries());
  for (int i=0; i<A.NumMyRows(); ++i)
  {
  int grid = A.RowMap().GID(i);
  int numentries;
  double* values;
  int* indices;
  A.ExtractMyRowView(i,numentries,values,indices);
  int snumentries=0;
  for (int j=0; j<numentries; ++j)
  {
  if (values[j]==0.0) continue;
  scindices[snumentries] = A.ColMap().GID(indices[j]);
  scvalues[snumentries] = values[j];
  snumentries++;
  }
  ML_az_sort(&scindices[0],snumentries,NULL,&scvalues[0]);
  int err = As.InsertGlobalValues(grid,snumentries,&scvalues[0],&scindices[0]);
  if (err) dserror("Epetra_CrsMatrix::InsertGlobalValues returned err=%d",err);
  }
  if (A.Filled()) As.FillComplete(A.DomainMap(),A.RangeMap(),true);
  return;
  }
*/


#endif  // #ifdef CCADISCRET
