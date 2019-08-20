/*----------------------------------------------------------------------*/
/*! \file

\brief Implementation

\level 2

\maintainer Martin Kronbichler
            http://www.lnm.mw.tum.de
            089 - 289-15235

*----------------------------------------------------------------------*/
#include "linalg_mlapi_operator.H"
#include "Epetra_Vector.h"

/*----------------------------------------------------------------------*
 |  ctor (public)                                            mwgee 10/07|
 *----------------------------------------------------------------------*/
LINALG::AMG_Operator::AMG_Operator(
    Teuchos::RCP<Epetra_RowMatrix> A, Teuchos::ParameterList& params, const bool compute)
    : Epetra_Operator(), label_("LINALG::AMG_Operator"), params_(params), nlevel_(0)
{
  Epetra_CrsMatrix* tmp = dynamic_cast<Epetra_CrsMatrix*>(A.get());
  if (!tmp) dserror("Expected Epetra_CrsMatrix");
  Ainput_ = Teuchos::rcp(tmp, false);
  SetupNonSymStab();
  return;
}


/*----------------------------------------------------------------------*
 | v cycle (private)                                         mwgee 10/07|
 *----------------------------------------------------------------------*/
void LINALG::AMG_Operator::Vcycle(
    const MLAPI::MultiVector& b_f, MLAPI::MultiVector& x_f, const int level) const
{
  // coarse grid solve
  if (level == NLevel() - 1)
  {
    x_f = S(level) * b_f;
    return;
  }

  MLAPI::MultiVector r_c(P(level).GetDomainSpace(), 1, false);
  MLAPI::MultiVector x_c(P(level).GetDomainSpace(), 1, true);

  // presmoothing
  x_f = 0.0;
  S(level).Apply(b_f, x_f);

  // compute residual and restrict it to next coarser grid
  r_c = R(level) * (b_f - A(level) * x_f);

  // solve coarser problem
  Vcycle(r_c, x_c, level + 1);

  // prolongate correction
  x_f = x_f + P(level) * x_c;

  // postsmooth
  S(level).Apply(b_f, x_f);

  return;
}

/*----------------------------------------------------------------------*
 |  setup phase (private)                                    mwgee 10/07|
 *----------------------------------------------------------------------*/
void LINALG::AMG_Operator::SetupNonSymStab()
{
  MLAPI::Init();

  //------------------------------------------------- get some parameters
  int maxlevels = Params().get<int>("max levels", 10);
  int maxcoarsesize = Params().get<int>("coarse: max size", 1);
  double* nullspace = Params().get<double*>("null space: vectors", (double*)NULL);
  if (!nullspace) dserror("No nullspace supplied in parameter list");
  int nsdim = Params().get<int>("null space: dimension", 1);
  double damping = Params().get<double>("aggregation: damping factor", 1.33);
  std::string eigenanalysis = Params().get("eigen-analysis: type", "cg");

  //--------------------------------------------- get input matrix wrapped
  MLAPI::Space space(A()->RowMatrixRowMap());
  MLAPI::Operator mlapiA(space, space, A().get(), false);

  //------------------------------------------------ get nullspace wrapped
  MLAPI::MultiVector NS;
  NS.Reshape(mlapiA.GetRangeSpace(), nsdim);
  for (int i = 0; i < nsdim; ++i)
  {
    const int length = NS.GetMyLength();
    for (int j = 0; j < length; ++j) NS(j, i) = nullspace[i * length + j];
  }

  //--------------------------------------- start constructing the hierarchy
  int level = 0;
  A_.resize(maxlevels);
  Astab_.resize(maxlevels);
  Anonstab_.resize(maxlevels);
  R_.resize(maxlevels - 1);
  P_.resize(maxlevels - 1);
  S_.resize(maxlevels);
  A_[0] = mlapiA;
  Astab_[0] = mlapiA;
  Anonstab_[0] = mlapiA;
  MLAPI::InverseOperator S;
  MLAPI::Operator Ptent;
  MLAPI::Operator Rtent;
  MLAPI::Operator P;
  MLAPI::Operator R;
  MLAPI::MultiVector NextNS;
  MLAPI::Operator Ac;
  // temporary parameter list for ml (using most of Params() but changing some)
  Teuchos::ParameterList params(Params());
  for (level = 0; level < maxlevels - 1; ++level)
  {
    //---------------------------------------stabilized current grid matrix
    mlapiA = A_[level];
    if (level) params.set("PDE equations", NS.GetNumVectors());
    params.set("workspace: current level", level);

    //------------------------------------------------------- build smoother
    Teuchos::ParameterList p;
    char levelstr[11];
    sprintf(levelstr, "(level %d)", level);
    std::string type =
        Params().get("smoother: type " + (std::string)levelstr, "symmetric Gauss-Seidel");
    if (type == "Jacobi" || type == "symmetric Gauss-Seidel")
    {
      int sweeps = Params().get("smoother: sweeps " + (std::string)levelstr, 1);
      double damping = Params().get("smoother: damping factor " + (std::string)levelstr, 1.0);
      p.set<int>("smoother: sweeps", sweeps);
      p.set<double>("smoother: damping factor", damping);
    }
    else if (type == "IFPACK")
    {
      type = Params().get("smoother: ifpack type " + (std::string)levelstr, "ILU");
      int lof = Params().sublist("smoother: ifpack list").get("fact: level-of-fill", 0);
      p.set<int>("fact: level-of-fill", lof);
    }
    else if (type == "Amesos-KLU")
      ;  // nothing to do
    else if (type == "Amesos-Superludist")
      dserror("Amesos-Superludist not supported by MLAPI");
    else if (type == "MLS")
    {
      int sweeps = Params().get("smoother: MLS polynomial order " + (std::string)levelstr, 1);
      p.set("smoother: MLS polynomial order", sweeps);
    }
    S.Reshape(mlapiA, type, p);

    //---------------------------------- build Ptent and next level nullspace
    GetPtent(mlapiA, params, NS, Ptent, NextNS);
    Ptent = -1.0 * Ptent;
    NS = -1.0 * NextNS;

    //----------------------------- build symmetric and nonsymmetric part of A
    MLAPI::Operator Asym = GetTranspose(mlapiA);
    Asym = 0.5 * (Asym + mlapiA);
    MLAPI::Operator Askew = mlapiA - Asym;

    //---------------------------------------build row norms of Asym and Askew
    MLAPI::MultiVector skewnorm = Row1Norm(Askew);
    MLAPI::MultiVector symnorm = Row1Norm(Asym);

    //------------------------------- build ratio between skewnorm and symnorm
    MLAPI::MultiVector ratio(mlapiA.GetRangeSpace(), 1, false);
    for (int i = 0; i < ratio.GetMyLength(); ++i) ratio(i) = skewnorm(i) / symnorm(i);
    // cout << "skewnorm\n" << skewnorm << "symnorm\n" << symnorm << "ratio\n" << ratio;
    // cout << "fineratio\n" << ratio << endl;

    //------------------------ build prolongation smoothing operators
    if (damping)
    {
      double lambdamax = 0.0;
      type = Params().get("eigen-analysis: type", "cg");
      if (type == "Anorm")
        lambdamax = MaxEigAnorm(mlapiA, true);
      else if (type == "cg")
        lambdamax = MaxEigCG(mlapiA, true);
      else if (type == "power-method")
        lambdamax = MaxEigPowerMethod(mlapiA, true);
      else
        dserror("Unknown type of eigenanalysis for MLAPI");
      // cout << "Max eigenvalue = " << lambdamax << endl;
      MLAPI::Operator I = GetIdentity(mlapiA.GetDomainSpace(), mlapiA.GetRangeSpace());
      MLAPI::MultiVector Diag = GetDiagonal(mlapiA);
      Diag.Reciprocal();
      Diag.Scale(damping / lambdamax);
      MLAPI::Operator Dinv = GetDiagonal(Diag);
      MLAPI::Operator IminuswDinvA = I - Dinv * mlapiA;
      MLAPI::Operator IminuswDinvAT = I - Dinv * GetTranspose(mlapiA);
      P = IminuswDinvA * Ptent;
      R = IminuswDinvAT * Ptent;
      R = GetTranspose(R);
      Rtent = GetTranspose(Ptent);
    }
    else
    {
      P = Ptent;
      Rtent = GetTranspose(Ptent);
      R = Rtent;
    }

    //------- coarsen matrix and build skew/sym part of coarse matrix
    //----------- use non stab. matrix to build nonstab coarse matrix
    Ac = GetRAP(R, Anonstab_[level], P);
    MLAPI::Operator Acsym = GetTranspose(Ac);
    Acsym = 0.5 * (Acsym + Ac);
    MLAPI::Operator Acskew = Ac - Acsym;

    //--------------------------- build row norms of Acsym and Acskew
    MLAPI::MultiVector cskewnorm = Row1Norm(Acskew);
    MLAPI::MultiVector csymnorm = Row1Norm(Acsym);

    //---------------------- ratio between coarse skew/sym norms (ist)
    MLAPI::MultiVector cratioist(Ac.GetRangeSpace(), 1, false);
    for (int i = 0; i < cratioist.GetMyLength(); ++i) cratioist(i) = cskewnorm(i) / csymnorm(i);
    // cout << "cskewnorm\n" << cskewnorm << "csymnorm\n" << csymnorm << "cratioist\n" << cratioist;
    // cout << "cratioist\n" << cratioist;

    //---------------------------- coarsen ratio from fine level (soll)
    MLAPI::MultiVector cratiosoll;
    cratiosoll = Rtent * ratio;
    // cout << "cratiosoll\n" << cratiosoll;

    //---------------------------------- create vector of boost factors
    MLAPI::MultiVector boost(Ac.GetRangeSpace(), 1, false);
    bool test = true;
    int printit = 0;
#if 1  // boost stuff above 1.0 only
    for (int i = 0; i < boost.GetMyLength(); ++i)
    {
      double nom = abs(cratioist(i));
      double denom = abs(cratiosoll(i));
      if (denom <= 1.0e-10) denom = 1.0e-08;
      boost(i) = nom / denom;
      if (boost(i) <= 1.0)
      {
        boost(i) = 0.0;
        continue;
      }
      if (boost(i) >= 2.0) boost(i) = 2.0;
      boost(i) = sqrt(boost(i) - 1.0);
      if (boost(i) != 0.0 && test)
      {
        test = false;
        printit = 1;
      }
    }
#endif
#if 0  // boost everything
    for (int i=0; i<boost.GetMyLength(); ++i)
    {
      double nom   = abs(cratioist(i));
      double denom = abs(cratiosoll(i));
      if (denom<=1.0e-10) denom = 1.0e-08;
      boost(i) = nom/denom;
      if (boost(i)!=0.0  && test)
      {
        test = false;
        printit  = 0;
      }
    }
#endif

    int printout = 0;
    Comm().SumAll(&printit, &printout, 1);
    test = true;

#if 0  // print nonzero boost vector
    if (printout)
    {
      if (!Comm().MyPID()) cout << "boost\n";
      fflush(stdout);
      cout << boost;
    }
#endif

    //---------------------- build the boosting matrix B = diag(boost)
    MLAPI::Operator B = GetDiagonal(boost);

    //--------------------------------- boost the symmetric part of Ac
    MLAPI::Operator Acstab = B * Acsym * B;
#if 0
    if (printout)
    {
      if (!Comm().MyPID()) cout << "boost matrix\n";
      fflush(stdout);
      cout << Acstab;
    }
#endif
    MLAPI::Operator Acnonstab = Ac;
    Ac = Ac + Acstab;

    //-------------------------------------- store values in hierarchy
    // maybe store the B*Acsym*B as well later
    P_[level] = P;
    R_[level] = R;
    S_[level] = S;
    A_[level + 1] = Ac;
    Astab_[level + 1] = Acstab;
    Anonstab_[level + 1] = Acnonstab;

    // break if coarsest level is below specified size
    if (Ac.GetRangeSpace().GetNumGlobalElements() <= maxcoarsesize)
    {
      ++level;
      break;
    }

  }  // for (level=0; level<maxlevels-1; ++level)


  //------------------------------------------------ setup coarse solver
  {
    Teuchos::ParameterList p;
    S_[level].Reshape(A_[level], "Amesos-KLU", p);
  }

  //---------------------------------------------- store number of levels
  nlevel_ = level + 1;

  return;
}

/*----------------------------------------------------------------------*
 |  apply operator (public)                                  mwgee 10/07|
 *----------------------------------------------------------------------*/
int LINALG::AMG_Operator::ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  // do not do anything for testing
  // Y.Update(1.0,X,0.0);
  // return 0;

  // create a space
  const Epetra_BlockMap& bmap = X.Map();
  MLAPI::Space space;
  space.Reshape(bmap.NumGlobalElements(), bmap.NumMyElements(), bmap.MyGlobalElements());

  // wrap incoming and outgoing vectors as MLAPI::MultiVector
  // note: Aztec might pass X and Y as physically identical objects,
  // so we deep copy here
  MLAPI::MultiVector in(space, X.Pointers(), 1);
  MLAPI::MultiVector out(space, Y.Pointers(), 1);
  MLAPI::MultiVector b(space, 1, false);
  MLAPI::MultiVector x(space, 1, true);
  const int mylength = X.Map().NumMyElements();
  for (int i = 0; i < mylength; ++i) b(i) = in(i);

  // call V cycle multigrid
  Vcycle(b, x, 0);

  // copy the solution back
  for (int i = 0; i < mylength; ++i) out(i) = x(i);


  return 0;
}

/*----------------------------------------------------------------------*
 |  build row wise 1-norm (private)                          mwgee 10/07|
 *----------------------------------------------------------------------*/
MLAPI::MultiVector LINALG::AMG_Operator::Row1Norm(MLAPI::Operator& A)
{
  MLAPI::MultiVector norm(A.GetRangeSpace(), 1, true);
  const Epetra_RowMatrix* B = A.GetRowMatrix();
  if (!B) dserror("Cannot get Epetra_RowMatrix from MLAPI Operator");
  int nummyrows = B->NumMyRows();
  int maxnumentries = B->MaxNumEntries() * 2;
  std::vector<double> vals(maxnumentries);
  std::vector<int> indices(maxnumentries);
  for (int i = 0; i < nummyrows; ++i)
  {
    int numentries;
    int err = B->ExtractMyRowCopy(i, maxnumentries, numentries, &vals[0], &indices[0]);
    if (err) dserror("Epetra_RowMatrix::ExtractMyRowCopy returned err=%d", err);
    double sum = 0.0;
    for (int j = 0; j < numentries; ++j) sum += abs(vals[j]);
    norm(i) = sum;
  }
  return norm;
}
