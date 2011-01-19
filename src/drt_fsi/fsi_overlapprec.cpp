#ifdef CCADISCRET

#include "fsi_overlapprec.H"
#include "fsi_debugwriter.H"
#include <Epetra_Time.h>
#include <EpetraExt_MatrixMatrix.h>

extern struct _GENPROB     genprob;
/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
FSI::BlockPreconditioningMatrix::BlockPreconditioningMatrix(Teuchos::RCP<UTILS::MonolithicDebugWriter> pcdbg,
                                                            const LINALG::MultiMapExtractor& maps,
                                                            ADAPTER::Structure& structure,
                                                            ADAPTER::Fluid& fluid,
                                                            ALE::Ale& ale,
                                                            int symmetric,
                                                            double omega,
                                                            int iterations,
                                                            double somega,
                                                            int siterations,
                                                            double fomega,
                                                            int fiterations,
                                                            double aomega,
                                                            int aiterations,
                                                            FILE* err)
  : LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>(maps,maps,81,false,true),
    symmetric_(symmetric),
    omega_(omega),
    iterations_(iterations),
    somega_(somega),
    siterations_(siterations),
    fomega_(fomega),
    fiterations_(fiterations),
    aomega_(aomega),
    aiterations_(aiterations),
    err_(err),
    pcdbg_(pcdbg)
{
  fluidsolver_ = Teuchos::rcp(new LINALG::Preconditioner(fluid.LinearSolver()));

#ifndef BLOCKMATRIXMERGE
  structuresolver_ = Teuchos::rcp(new LINALG::Preconditioner(structure.LinearSolver()));

  constalesolver_ = ale.ConstPreconditioner();
  if (constalesolver_==Teuchos::null)
    alesolver_ = Teuchos::rcp(new LINALG::Preconditioner(ale.LinearSolver()));
  else
    alesolver_ = constalesolver_;
#endif
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int FSI::BlockPreconditioningMatrix::ApplyInverse(const Epetra_MultiVector &X, Epetra_MultiVector &Y) const
{
  if (UseTranspose())
    dserror("no transpose preconditioning");

  Teuchos::RCP<Epetra_Vector> r;

  if (pcdbg_!=Teuchos::null)
  {
    pcdbg_->NewIteration();

    // X and Y are the same at this point (if we have been called by aztec!)
    Epetra_Vector &y = Teuchos::dyn_cast<Epetra_Vector>(Y);
    pcdbg_->WriteVector("x",Teuchos::rcp(&y,false));

    r = Teuchos::rcp(new Epetra_Vector(y.Map()));
    Apply(X,*r);
  }

#ifdef BLOCKMATRIXMERGE
  MergeSolve(X, Y);
#else
  SGS(X, Y);
#endif

  if (pcdbg_!=Teuchos::null)
  {
    Epetra_Vector &y = Teuchos::dyn_cast<Epetra_Vector>(Y);
    pcdbg_->WriteVector("y",Teuchos::rcp(&y,false));
    r->Update(-1,y,1);
    pcdbg_->WriteVector("r",r);
  }

  return 0;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void FSI::BlockPreconditioningMatrix::MergeSolve(const Epetra_MultiVector &X, Epetra_MultiVector &Y) const
{
#ifdef BLOCKMATRIXMERGE
  const Epetra_Vector &x = Teuchos::dyn_cast<const Epetra_Vector>(X);
  Epetra_Vector &y = Teuchos::dyn_cast<Epetra_Vector>(Y);

#if 0
  const std::string fname = "mergedmatrix.mtl";
  cout << sparse_->RowMap()<<endl;
  cout << sparse_->RangeMap()<<endl;
  cout << sparse_->DomainMap()<<endl;
  cout << x.Map()<<endl;
#endif


  fluidsolver_->Solve(sparse_->EpetraMatrix(),
                      Teuchos::rcp(&y,false),
                      Teuchos::rcp(new Epetra_Vector(x)),
                      true);
#endif
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void FSI::BlockPreconditioningMatrix::SetupPreconditioner()
{
#ifdef BLOCKMATRIXMERGE

  // this is really evil :)
  sparse_ = Merge();
  fluidsolver_->Setup(sparse_->EpetraMatrix());

#if 0
  Matrix(0,0).Dump("dump-struct");
  Matrix(1,1).Dump("dump-fluid");
  Matrix(2,2).Dump("dump-ale");

  static int count;
  count++;
  std::stringstream s;
  s << "dump-" << count;
  cout << "write: " << s.str() << "\n";
  sparse_->Dump(s.str());

//   Epetra_Vector diagonal(sparse_->RowMap());
//   int err = sparse_->ExtractDiagonalCopy(diagonal);
//   diagonal.Print(cout);
#endif

#else
  const LINALG::SparseMatrix& structInnerOp = Matrix(0,0);
  const LINALG::SparseMatrix& fluidInnerOp  = Matrix(1,1);
  const LINALG::SparseMatrix& aleInnerOp    = Matrix(2,2);

  structuresolver_->Setup(structInnerOp.EpetraMatrix());
  fluidsolver_    ->Setup(fluidInnerOp .EpetraMatrix());
  if (constalesolver_==Teuchos::null)
    alesolver_    ->Setup(aleInnerOp   .EpetraMatrix());
#endif
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void FSI::BlockPreconditioningMatrix::LocalBlockRichardson(Teuchos::RCP<LINALG::Preconditioner> solver,
                                                           const LINALG::SparseMatrix& innerOp,
                                                           Teuchos::RCP<Epetra_Vector> x,
                                                           Teuchos::RCP<Epetra_Vector> y,
                                                           Teuchos::RCP<Epetra_Vector> tmpx,
                                                           int iterations,
                                                           double omega,
                                                           FILE* err,
                                                           const Epetra_Comm& comm)
{
  if (iterations > 0)
  {
    y->Scale(omega);
    Teuchos::RCP<Epetra_Vector> tmpy = Teuchos::rcp(new Epetra_Vector(y->Map()));
    if (err!=NULL)
      if (comm.MyPID()==0)
        fprintf(err,"    fluid richardson (%d,%f):",iterations,omega);
    for (int i=0; i<iterations; ++i)
    {
      innerOp.EpetraMatrix()->Multiply(false,*y,*tmpx);
      tmpx->Update(1.0,*x,-1.0);

      if (err!=NULL)
      {
        double n;
        tmpx->Norm2(&n);
        if (comm.MyPID()==0)
          fprintf(err," %e",n);
      }

      solver->Solve(innerOp.EpetraMatrix(),tmpy,tmpx,false);
      y->Update(omega,*tmpy,1.0);
    }
    if (err!=NULL)
      if (comm.MyPID()==0)
        fprintf(err,"\n");
  }
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
FSI::OverlappingBlockMatrix::OverlappingBlockMatrix(Teuchos::RCP<UTILS::MonolithicDebugWriter> pcdbg,
                                                    const LINALG::MultiMapExtractor& maps,
                                                    ADAPTER::Structure& structure,
                                                    ADAPTER::Fluid& fluid,
                                                    ALE::Ale& ale,
                                                    bool structuresplit,
                                                    int symmetric,
                                                    double omega,
                                                    int iterations,
                                                    double somega,
                                                    int siterations,
                                                    double fomega,
                                                    int fiterations,
                                                    double aomega,
                                                    int aiterations,
                                                    FILE* err)
  : BlockPreconditioningMatrix(pcdbg,
                               maps,
                               structure,
                               fluid,
                               ale,
                               symmetric,
                               omega,
                               iterations,
                               somega,
                               siterations,
                               fomega,
                               fiterations,
                               aomega,
                               aiterations,
                               err),
    structuresplit_(structuresplit),
    structure_(structure),
    fluid_(fluid),
    ale_(ale)
{
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void FSI::OverlappingBlockMatrix::SetupPreconditioner()
{
#ifdef BLOCKMATRIXMERGE

  // this is really evil :)
  sparse_ = Merge();
  fluidsolver_->Setup(sparse_->EpetraMatrix());

#if 0
  Matrix(0,0).Dump("dump-struct");
  Matrix(1,1).Dump("dump-fluid");
  Matrix(2,2).Dump("dump-ale");

  static int count;
  count++;
  std::stringstream s;
  s << "dump-" << count;
  cout << "write: " << s.str() << "\n";
  sparse_->Dump(s.str());

//   Epetra_Vector diagonal(sparse_->RowMap());
//   int err = sparse_->ExtractDiagonalCopy(diagonal);
//   diagonal.Print(cout);
#endif

#else
  const LINALG::SparseMatrix& structInnerOp = Matrix(0,0);
  const LINALG::SparseMatrix& fluidInnerOp  = Matrix(1,1);
  const LINALG::SparseMatrix& aleInnerOp    = Matrix(2,2);

  RCP<LINALG::MapExtractor> fsidofmapex = null;
  RCP<Epetra_Map>           irownodes = null;
#if 0
  if (structuresplit_)
  {
    RCP<DRT::Discretization> dis      = fluid_.Discretization();
    const Epetra_Map*        nrowmap  = dis->NodeRowMap();
    DRT::Condition*          fsi      = dis->GetCondition("FSICoupling");
    if (!fsi) dserror("there should be an FSI interface condition in here...");
    const vector<int>*       fsinodes = fsi->Nodes();
    vector<int>              fsidofs;
    vector<int>              fsirownodes;
    for (int i=0; i<(int)fsinodes->size(); ++i)
      if (nrowmap->MyGID((*fsinodes)[i]))
      {
        fsirownodes.push_back((*fsinodes)[i]);
        DRT::Node* node  = dis->gNode((*fsinodes)[i]);
        const int numdof = dis->NumDof(node);
        for (int j=0; j<numdof; ++j)
          fsidofs.push_back(dis->Dof(node,j));
      }
    RCP<Epetra_Map> fsimap =
      rcp(new Epetra_Map(-1,(int)fsidofs.size(),&fsidofs[0],0,dis->DofRowMap()->Comm()));
    irownodes = rcp(new Epetra_Map(-1,(int)fsirownodes.size(),&(fsirownodes[0]),0,fsimap->Comm()));
    fsidofmapex = rcp(new LINALG::MapExtractor(*(dis->DofRowMap()),fsimap));
  }
#endif

  structuresolver_->Setup(structInnerOp.EpetraMatrix());
  fluidsolver_->Setup(fluidInnerOp.EpetraMatrix(),
                      fsidofmapex,
                      fluid_.Discretization(),
                      irownodes,
                      structuresplit_);
  if (constalesolver_==Teuchos::null)
    alesolver_->Setup(aleInnerOp.EpetraMatrix());
#endif
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void FSI::OverlappingBlockMatrix::SGS(const Epetra_MultiVector &X, Epetra_MultiVector &Y) const
{
  const LINALG::SparseMatrix& structInnerOp = Matrix(0,0);
  const LINALG::SparseMatrix& fluidInnerOp  = Matrix(1,1);
  const LINALG::SparseMatrix& fluidMeshOp   = Matrix(1,2);
  const LINALG::SparseMatrix& fluidBoundOp  = Matrix(1,0);
  const LINALG::SparseMatrix& aleInnerOp    = Matrix(2,2);

  // Extract vector blocks
  // RHS

  const Epetra_Vector &x = Teuchos::dyn_cast<const Epetra_Vector>(X);

  // initial guess

  Epetra_Vector &y = Teuchos::dyn_cast<Epetra_Vector>(Y);

  Teuchos::RCP<Epetra_Vector> sy = RangeExtractor().ExtractVector(y,0);
  Teuchos::RCP<Epetra_Vector> fy = RangeExtractor().ExtractVector(y,1);
  Teuchos::RCP<Epetra_Vector> ay = RangeExtractor().ExtractVector(y,2);

  Teuchos::RCP<Epetra_Vector> sz = Teuchos::rcp(new Epetra_Vector(sy->Map()));
  Teuchos::RCP<Epetra_Vector> fz = Teuchos::rcp(new Epetra_Vector(fy->Map()));
  Teuchos::RCP<Epetra_Vector> az = Teuchos::rcp(new Epetra_Vector(ay->Map()));

  Teuchos::RCP<Epetra_Vector> tmpsx = Teuchos::rcp(new Epetra_Vector(DomainMap(0)));
  Teuchos::RCP<Epetra_Vector> tmpfx = Teuchos::rcp(new Epetra_Vector(DomainMap(1)));
  Teuchos::RCP<Epetra_Vector> tmpax = Teuchos::rcp(new Epetra_Vector(DomainMap(2)));

  // outer Richardson loop
  for (int run=0; run<iterations_; ++run)
  {
    Teuchos::RCP<Epetra_Vector> sx = DomainExtractor().ExtractVector(x,0);
    Teuchos::RCP<Epetra_Vector> fx = DomainExtractor().ExtractVector(x,1);
    Teuchos::RCP<Epetra_Vector> ax = DomainExtractor().ExtractVector(x,2);

    // ----------------------------------------------------------------
    // lower GS

    {
      if (run>0)
      {
        const LINALG::SparseMatrix& structBoundOp  = Matrix(0,1);

        structInnerOp.Multiply(false,*sy,*tmpsx);
        sx->Update(-1.0,*tmpsx,1.0);
        structBoundOp.Multiply(false,*fy,*tmpsx);
        sx->Update(-1.0,*tmpsx,1.0);
      }

      // Solve structure equations for sy with the rhs sx
      structuresolver_->Solve(structInnerOp.EpetraMatrix(),sz,sx,true);
      // do Richardson iteration
      LocalBlockRichardson(structuresolver_,structInnerOp,sx,sz,tmpsx,siterations_,somega_,err_,Comm());

      if (run>0)
      {
        sy->Update(omega_,*sz,1.0);
      }
      else
      {
        sy->Update(omega_,*sz,0.0);
      }
    }

    {
      // Solve ale equations for ay with the rhs ax - A(I,Gamma) sy

      if (run>0)
      {
        aleInnerOp.Multiply(false,*ay,*tmpax);
        ax->Update(-1.0,*tmpax,1.0);
      }

      if (structuresplit_)
      {
        if (run>0)
        {
          const LINALG::SparseMatrix& aleBoundOp    = Matrix(2,1);
          aleBoundOp.Multiply(false,*fy,*tmpax);
          ax->Update(-1.0,*tmpax,1.0);
        }
      }
      else
      {
        const LINALG::SparseMatrix& aleBoundOp    = Matrix(2,0);
        aleBoundOp.Multiply(false,*sy,*tmpax);
        ax->Update(-1.0,*tmpax,1.0);
      }

      alesolver_->Solve(aleInnerOp.EpetraMatrix(),az,ax,true);
      // do Richardson iteration
      LocalBlockRichardson(alesolver_,aleInnerOp,ax,az,tmpax,aiterations_,aomega_,err_,Comm());

      if (run>0)
      {
        ay->Update(omega_,*az,1.0);
      }
      else
      {
        ay->Update(omega_,*az,0.0);
      }
    }

    {
      // Solve fluid equations for fy with the rhs fx - F(I,Gamma) sy - F(Mesh) ay

      if (run>0)
      {
        fluidInnerOp.Multiply(false,*fy,*tmpfx);
        fx->Update(-1.0,*tmpfx,1.0);
      }

      fluidBoundOp.Multiply(false,*sy,*tmpfx);
      fx->Update(-1.0,*tmpfx,1.0);
      fluidMeshOp.Multiply(false,*ay,*tmpfx);
      fx->Update(-1.0,*tmpfx,1.0);
      fluidsolver_->Solve(fluidInnerOp.EpetraMatrix(),fz,fx,true);

      LocalBlockRichardson(fluidsolver_,fluidInnerOp,fx,fz,tmpfx,fiterations_,fomega_,err_,Comm());

      if (run>0)
      {
        fy->Update(omega_,*fz,1.0);
      }
      else
      {
        fy->Update(omega_,*fz,0.0);
      }
    }

    // ----------------------------------------------------------------
    // the symmetric part of the pc can be skipped

    if (symmetric_)
    {

      sx = DomainExtractor().ExtractVector(x,0);
      fx = DomainExtractor().ExtractVector(x,1);
      ax = DomainExtractor().ExtractVector(x,2);

      // ----------------------------------------------------------------
      // upper GS

      {
        fluidInnerOp.Multiply(false,*fy,*tmpfx);
        fx->Update(-1.0,*tmpfx,1.0);
        fluidBoundOp.Multiply(false,*sy,*tmpfx);
        fx->Update(-1.0,*tmpfx,1.0);
        fluidMeshOp.Multiply(false,*ay,*tmpfx);
        fx->Update(-1.0,*tmpfx,1.0);

        fluidsolver_->Solve(fluidInnerOp.EpetraMatrix(),fz,fx,true);

        LocalBlockRichardson(fluidsolver_,fluidInnerOp,fx,fz,tmpfx,fiterations_,fomega_,err_,Comm());
        fy->Update(omega_,*fz,1.0);
      }

      {
        aleInnerOp.Multiply(false,*ay,*tmpax);
        ax->Update(-1.0,*tmpax,1.0);
        if (structuresplit_)
        {
          const LINALG::SparseMatrix& aleBoundOp    = Matrix(2,1);
          aleBoundOp.Multiply(false,*fy,*tmpax);
          ax->Update(-1.0,*tmpax,1.0);
        }
        else
        {
          const LINALG::SparseMatrix& aleBoundOp    = Matrix(2,0);
          aleBoundOp.Multiply(false,*sy,*tmpax);
          ax->Update(-1.0,*tmpax,1.0);
        }
        alesolver_->Solve(aleInnerOp.EpetraMatrix(),az,ax,true);

        LocalBlockRichardson(alesolver_,aleInnerOp,ax,az,tmpax,aiterations_,aomega_,err_,Comm());
        ay->Update(omega_,*az,1.0);
      }

      {
        const LINALG::SparseMatrix& structBoundOp  = Matrix(0,1);

        structInnerOp.Multiply(false,*sy,*tmpsx);
        sx->Update(-1.0,*tmpsx,1.0);
        structBoundOp.Multiply(false,*fy,*tmpsx);
        sx->Update(-1.0,*tmpsx,1.0);

        structuresolver_->Solve(structInnerOp.EpetraMatrix(),sz,sx,true);

        LocalBlockRichardson(structuresolver_,structInnerOp,sx,sz,tmpsx,siterations_,somega_,err_,Comm());
        sy->Update(omega_,*sz,1.0);
      }
    }
  }

  RangeExtractor().InsertVector(*sy,0,y);
  RangeExtractor().InsertVector(*fy,1,y);
  RangeExtractor().InsertVector(*ay,2,y);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
const char* FSI::OverlappingBlockMatrix::Label() const
{
  return "FSI::OverlappingBlockMatrix";
}



// /*----------------------------------------------------------------------*
//  *----------------------------------------------------------------------*/
FSI::LungOverlappingBlockMatrix::LungOverlappingBlockMatrix(const LINALG::MultiMapExtractor& maps,
                                                            ADAPTER::Structure& structure,
                                                            ADAPTER::Fluid& fluid,
                                                            ALE::Ale& ale,
                                                            bool structuresplit,
                                                            int symmetric,
                                                            double omega,
                                                            int iterations,
                                                            double somega,
                                                            int siterations,
                                                            double fomega,
                                                            int fiterations,
                                                            double aomega,
                                                            int aiterations,
                                                            FILE* err)
  : OverlappingBlockMatrix(Teuchos::null,
                           maps,
                           structure,
                           fluid,
                           ale,
                           structuresplit,
                           symmetric,
                           omega,
                           iterations,
                           somega,
                           siterations,
                           fomega,
                           fiterations,
                           aomega,
                           aiterations,
                           err)
{
  // determine map of all dofs not related to constraint

  std::vector<Teuchos::RCP<const Epetra_Map> > fsimaps;
  fsimaps.push_back(maps.Map(0));
  fsimaps.push_back(maps.Map(1));
  fsimaps.push_back(maps.Map(2));
  overallfsimap_ = LINALG::MultiMapExtractor::MergeMaps(fsimaps);
  fsiextractor_ = LINALG::MultiMapExtractor(*overallfsimap_, fsimaps);

  StructSchur_ = Teuchos::rcp(new LungSchurComplement());
  FluidSchur_ = Teuchos::rcp(new LungSchurComplement());

  // stuff needed for SIMPLE preconditioner
  const Teuchos::ParameterList& fsidyn   = DRT::Problem::Instance()->FSIDynamicParams();
  alpha_ = fsidyn.sublist("CONSTRAINT").get<double>("ALPHA");
  simpleiter_ = fsidyn.sublist("CONSTRAINT").get<int>("SIMPLEITER");
  prec_ = Teuchos::getIntegralValue<INPAR::FSI::PrecConstr>(fsidyn.sublist("CONSTRAINT"),"PRECONDITIONER");

  RCP<Teuchos::ParameterList> constrsolvparams = rcp(new Teuchos::ParameterList);
  constrsolvparams->set("solver","umfpack");
  constraintsolver_ = rcp(new LINALG::Solver(constrsolvparams,
                                             maps.Map(0)->Comm(),
                                             DRT::Problem::Instance()->ErrorFile()->Handle()));
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void FSI::LungOverlappingBlockMatrix::SetupPreconditioner()
{
#ifdef BLOCKMATRIXMERGE

  // this is really evil :)
  // note that this only works with UMFPACK as fluid solver (saddle
  // point problem!)

  sparse_ = Merge();
  fluidsolver_->Setup(sparse_->EpetraMatrix());

#else

  const LINALG::SparseMatrix& structInnerOp = Matrix(0,0);
  const LINALG::SparseMatrix& fluidInnerOp  = Matrix(1,1);
  const LINALG::SparseMatrix& aleInnerOp    = Matrix(2,2);

  RCP<LINALG::MapExtractor> fsidofmapex = null;
  RCP<Epetra_Map>           irownodes = null;

  structuresolver_->Setup(structInnerOp.EpetraMatrix());
  fluidsolver_->Setup(fluidInnerOp.EpetraMatrix(),
                      fsidofmapex,
                      fluid_.Discretization(),
                      irownodes,
                      structuresplit_);
  if (constalesolver_==Teuchos::null)
    alesolver_->Setup(aleInnerOp.EpetraMatrix());

#endif
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void FSI::LungOverlappingBlockMatrix::SGS(const Epetra_MultiVector &X, Epetra_MultiVector &Y) const
{
  const LINALG::SparseMatrix& StructInnerOp  = Matrix(0,0);
  const LINALG::SparseMatrix& StructBoundOp  = Matrix(0,1);
  const LINALG::SparseMatrix& StructConOp    = Matrix(0,3);
  const LINALG::SparseMatrix& FluidBoundOp   = Matrix(1,0);
  const LINALG::SparseMatrix& FluidInnerOp   = Matrix(1,1);
  const LINALG::SparseMatrix& FluidMeshOp    = Matrix(1,2);
  const LINALG::SparseMatrix& FluidConOp     = Matrix(1,3);
  const LINALG::SparseMatrix& AleInnerOp     = Matrix(2,2);
  const LINALG::SparseMatrix& ConStructOp    = Matrix(3,0);
  const LINALG::SparseMatrix& ConFluidOp     = Matrix(3,1);


  // Extract vector blocks
  // RHS

  const Epetra_Vector &x = Teuchos::dyn_cast<const Epetra_Vector>(X);

  // initial guess

  Epetra_Vector &y = Teuchos::dyn_cast<Epetra_Vector>(Y);

  Teuchos::RCP<Epetra_Vector> sy = RangeExtractor().ExtractVector(y,0);
  Teuchos::RCP<Epetra_Vector> fy = RangeExtractor().ExtractVector(y,1);
  Teuchos::RCP<Epetra_Vector> ay = RangeExtractor().ExtractVector(y,2);
  Teuchos::RCP<Epetra_Vector> cy = RangeExtractor().ExtractVector(y,3);

  Teuchos::RCP<Epetra_Vector> sz = Teuchos::rcp(new Epetra_Vector(sy->Map()));
  Teuchos::RCP<Epetra_Vector> fz = Teuchos::rcp(new Epetra_Vector(fy->Map()));
  Teuchos::RCP<Epetra_Vector> az = Teuchos::rcp(new Epetra_Vector(ay->Map()));

  Teuchos::RCP<Epetra_Vector> tmpsx = Teuchos::rcp(new Epetra_Vector(DomainMap(0)));
  Teuchos::RCP<Epetra_Vector> tmpfx = Teuchos::rcp(new Epetra_Vector(DomainMap(1)));
  Teuchos::RCP<Epetra_Vector> tmpax = Teuchos::rcp(new Epetra_Vector(DomainMap(2)));

  for (int outerrun=0; outerrun<simpleiter_; ++outerrun)
  {

    // -------------------------------------------------------------------
    // intermediate fsi dofs: u_(n+1/2) = F^(-1)(f-B^T*lambda_n)
    // -------------------------------------------------------------------

    // inner Richardson loop (FSI block)

    for (int run=0; run<iterations_; ++run)
    {
      Teuchos::RCP<Epetra_Vector> sx = DomainExtractor().ExtractVector(x,0);
      Teuchos::RCP<Epetra_Vector> fx = DomainExtractor().ExtractVector(x,1);
      Teuchos::RCP<Epetra_Vector> ax = DomainExtractor().ExtractVector(x,2);

      // ----------------------------------------------------------------
      // lower GS


      // Structure
      {
        StructConOp.Multiply(false,*cy,*tmpsx);
        sx->Update(-1.0, *tmpsx, 1.0);

        if (run>0 or outerrun>0)
        {
          StructInnerOp.Multiply(false,*sy,*tmpsx);
          sx->Update(-1.0,*tmpsx,1.0);
          StructBoundOp.Multiply(false,*fy,*tmpsx);
          sx->Update(-1.0,*tmpsx,1.0);
        }

        // Solve structure equations for sy with the rhs sx
        structuresolver_->Solve(StructInnerOp.EpetraMatrix(),sz,sx,true);

        // do Richardson iteration
        LocalBlockRichardson(structuresolver_,StructInnerOp,sx,sz,tmpsx,siterations_,somega_,err_,Comm());

        if (run>0 or outerrun>0)
        {
          sy->Update(omega_,*sz,1.0);
        }
        else
        {
          sy->Update(omega_,*sz,0.0);
        }
      }

      // Ale
      {
        // Solve ale equations for ay with the rhs ax - A(I,Gamma) sy

        if (run>0 or outerrun>0)
        {
          AleInnerOp.Multiply(false,*ay,*tmpax);
          ax->Update(-1.0,*tmpax,1.0);
        }
        if (structuresplit_)
        {
          if (run>0 or outerrun>0)
          {
            const LINALG::SparseMatrix& aleFBoundOp = Matrix(2,1);
            aleFBoundOp.Multiply(false,*fy,*tmpax);
            ax->Update(-1.0,*tmpax,1.0);
          }

          const LINALG::SparseMatrix& aleSBoundOp   = Matrix(2,0);
          aleSBoundOp.Multiply(false,*sy,*tmpax);
          ax->Update(-1.0,*tmpax,1.0);
        }
        else
        {
          const LINALG::SparseMatrix& aleFBoundOp   = Matrix(2,0);
          aleFBoundOp.Multiply(false,*sy,*tmpax);
          ax->Update(-1.0,*tmpax,1.0);
        }

        alesolver_->Solve(AleInnerOp.EpetraMatrix(),az,ax,true);

        // do Richardson iteration
        LocalBlockRichardson(alesolver_,AleInnerOp,ax,az,tmpax,aiterations_,aomega_,err_,Comm());

        if (run>0 or outerrun>0)
        {
          ay->Update(omega_,*az,1.0);
        }
        else
        {
          ay->Update(omega_,*az,0.0);
        }
      }

      // Fluid
      {
        // Solve fluid equations for fy with the rhs fx - F(I,Gamma) sy - F(Mesh) ay - F(Constr) cy

        if (run>0 or outerrun>0)
        {
          FluidInnerOp.Multiply(false,*fy,*tmpfx);
          fx->Update(-1.0,*tmpfx,1.0);
        }

        FluidBoundOp.Multiply(false,*sy,*tmpfx);
        fx->Update(-1.0,*tmpfx,1.0);
        FluidMeshOp.Multiply(false,*ay,*tmpfx);
        fx->Update(-1.0,*tmpfx,1.0);
        FluidConOp.Multiply(false,*cy, *tmpfx);
        fx->Update(-1.0, *tmpfx, 1.0);

        fluidsolver_->Solve(FluidInnerOp.EpetraMatrix(),fz,fx,true);

        LocalBlockRichardson(fluidsolver_,FluidInnerOp,fx,fz,tmpfx,fiterations_,fomega_,err_,Comm());

        if (run>0 or outerrun>0)
        {
          fy->Update(omega_,*fz,1.0);
        }
        else
        {
          fy->Update(omega_,*fz,0.0);
        }
      }

      // ----------------------------------------------------------------
      // the symmetric part of the pc can be skipped

      if (symmetric_)
      {
        sx = DomainExtractor().ExtractVector(x,0);
        fx = DomainExtractor().ExtractVector(x,1);
        ax = DomainExtractor().ExtractVector(x,2);

        // ----------------------------------------------------------------
        // upper GS

        {
          FluidInnerOp.Multiply(false,*fy,*tmpfx);
          fx->Update(-1.0,*tmpfx,1.0);
          FluidBoundOp.Multiply(false,*sy,*tmpfx);
          fx->Update(-1.0,*tmpfx,1.0);
          FluidMeshOp.Multiply(false,*ay,*tmpfx);
          fx->Update(-1.0,*tmpfx,1.0);
          FluidConOp.Multiply(false,*cy,*tmpfx);
          fx->Update(-1.0, *tmpfx, 1.0);

          fluidsolver_->Solve(FluidInnerOp.EpetraMatrix(),fz,fx,true);

          LocalBlockRichardson(fluidsolver_,FluidInnerOp,fx,fz,tmpfx,fiterations_,fomega_,err_,Comm());
          fy->Update(omega_,*fz,1.0);
        }

        {
          AleInnerOp.Multiply(false,*ay,*tmpax);
          ax->Update(-1.0,*tmpax,1.0);

          if (structuresplit_)
          {
            const LINALG::SparseMatrix& aleFBoundOp = Matrix(2,1);
            aleFBoundOp.Multiply(false,*fy,*tmpax);
            ax->Update(-1.0,*tmpax,1.0);

            const LINALG::SparseMatrix& aleSBoundOp = Matrix(2,0);
            aleSBoundOp.Multiply(false,*sy,*tmpax);
            ax->Update(-1.0,*tmpax,1.0);
          }
          else
          {
            const LINALG::SparseMatrix& aleFBoundOp = Matrix(2,0);
            aleFBoundOp.Multiply(false,*sy,*tmpax);
            ax->Update(-1.0,*tmpax,1.0);
          }

          alesolver_->Solve(AleInnerOp.EpetraMatrix(),az,ax,true);

          LocalBlockRichardson(alesolver_,AleInnerOp,ax,az,tmpax,aiterations_,aomega_,err_,Comm());

          ay->Update(omega_,*az,1.0);
        }

        {
          StructInnerOp.Multiply(false,*sy,*tmpsx);
          sx->Update(-1.0,*tmpsx,1.0);
          StructBoundOp.Multiply(false,*fy,*tmpsx);
          sx->Update(-1.0,*tmpsx,1.0);
          StructConOp.Multiply(false,*cy,*tmpsx);
          sx->Update(-1.0, *tmpsx, 1.0);


          structuresolver_->Solve(StructInnerOp.EpetraMatrix(),sz,sx,true);

          LocalBlockRichardson(structuresolver_,StructInnerOp,sx,sz,tmpsx,siterations_,somega_,err_,Comm());
          sy->Update(omega_,*sz,1.0);
        }
      }
    }

    // -----------------------------------------------------------------------
    // intermediate constraint dofs: Dlambda~ = S^(-1) * (cx - B^ * u_(n+1/2))
    // -----------------------------------------------------------------------

    // cx - B^ * u_(n+1/2)

    Teuchos::RCP<Epetra_Vector> cx = DomainExtractor().ExtractVector(x,3);

    Epetra_Vector inter(cx->Map());
    ConStructOp.Multiply(false,*sy, inter);
    cx->Update(-1.0, inter, 1.0);
    ConFluidOp.Multiply(false,*fy, inter);
    cx->Update(-1.0, inter, 1.0);


    LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy> invDiag(fsiextractor_, fsiextractor_, 1, false, true);

    if (prec_ == INPAR::FSI::Simple)
    {
      // D^{-1} = diag(A(0,0))^{-1}

      Epetra_Vector structDiagVec(StructInnerOp.RowMap(),false);
      StructInnerOp.ExtractDiagonalCopy(structDiagVec);
      int err = structDiagVec.Reciprocal(structDiagVec);
      if (err) dserror("Epetra_MultiVector::Reciprocal (structure matrix) returned %d",err);
      LINALG::SparseMatrix invstructDiag(structDiagVec);

      Epetra_Vector fluidDiagVec(FluidInnerOp.RowMap(),false);
      FluidInnerOp.ExtractDiagonalCopy(fluidDiagVec);
      err = fluidDiagVec.Reciprocal(fluidDiagVec);
      if (err) dserror("Epetra_MultiVector::Reciprocal (fluid matrix) returned %d",err);
      LINALG::SparseMatrix invfluidDiag(fluidDiagVec);

      Epetra_Vector aleDiagVec(AleInnerOp.RowMap(),false);
      AleInnerOp.ExtractDiagonalCopy(aleDiagVec);
      err = aleDiagVec.Reciprocal(aleDiagVec);
      if (err) dserror("Epetra_MultiVector::Reciprocal (ale matrix) returned %d",err);
      LINALG::SparseMatrix invaleDiag(aleDiagVec);

      invDiag.Assign(0,0,View,invstructDiag);
      invDiag.Assign(1,1,View,invfluidDiag);
      invDiag.Assign(2,2,View,invaleDiag);
    }
    else if (prec_ == INPAR::FSI::Simplec)
    {
      // D^{-1} = sum(abs(A(0,0)))^{-1}

      Epetra_Vector structDiagVec(StructInnerOp.RowMap(),false);
      int err = StructInnerOp.EpetraMatrix()->InvRowSums(structDiagVec);
      if (err) dserror("Epetra_CrsMatrix::InvRowSums (structure matrix) returned %d",err);
      LINALG::SparseMatrix invstructDiag(structDiagVec);

      Epetra_Vector fluidDiagVec(FluidInnerOp.RowMap(),false);
      err = FluidInnerOp.EpetraMatrix()->InvRowSums(fluidDiagVec);
      if (err) dserror("Epetra_CrsMatrix::InvRowSums (fluid matrix) returned %d",err);
      LINALG::SparseMatrix invfluidDiag(fluidDiagVec);

      Epetra_Vector aleDiagVec(AleInnerOp.RowMap(),false);
      err =  AleInnerOp.EpetraMatrix()->InvRowSums(aleDiagVec);
      if (err) dserror("Epetra_CrsMatrix::InvRowSums (ale matrix) returned %d",err);
      LINALG::SparseMatrix invaleDiag(aleDiagVec);

      invDiag.Assign(0,0,View,invstructDiag);
      invDiag.Assign(1,1,View,invfluidDiag);
      invDiag.Assign(2,2,View,invaleDiag);
    }
    else
      dserror("Unknown type of preconditioner for constraint fsi system");

    invDiag.Complete();


    // S = - B^ * D^{-1} * B^T

    Teuchos::RCP<LINALG::SparseMatrix> interconA = StructSchur_->CalculateSchur(Matrix(3,0),invDiag.Matrix(0,0),Matrix(0,3));
    Teuchos::RCP<LINALG::SparseMatrix> temp = FluidSchur_->CalculateSchur(Matrix(3,1),invDiag.Matrix(1,1),Matrix(1,3));

    interconA->Add(*temp,false,1.0,1.0);

    interconA->Complete(StructConOp.DomainMap(),ConStructOp.RangeMap());
    interconA->Scale(-1.0);

    Teuchos::RCP<Epetra_Vector> interconsol = rcp(new Epetra_Vector(ConStructOp.RangeMap()));
    constraintsolver_->Solve(interconA->EpetraOperator(),interconsol,cx,true,true);

    // -------------------------------------------------------------------
    // update of all dofs
    // -------------------------------------------------------------------

    cy->Update(alpha_, *interconsol, 1.0);

    Teuchos::RCP<Epetra_Vector> temp1;
    Teuchos::RCP<Epetra_Vector> temp2;

    temp1 = rcp(new Epetra_Vector(sy->Map()));
    temp2 = rcp(new Epetra_Vector(sy->Map()));
    StructConOp.Multiply(false,*interconsol, *temp1);
    temp1->Scale(alpha_);
    invDiag.Matrix(0,0).Multiply(false,*temp1, *temp2);
    sy->Update(-1.0, *temp2, 1.0);

    temp1 = rcp(new Epetra_Vector(fy->Map()));
    temp2 = rcp(new Epetra_Vector(fy->Map()));
    FluidConOp.Multiply(false,*interconsol, *temp1);
    temp1->Scale(alpha_);
    invDiag.Matrix(1,1).Multiply(false,*temp1, *temp2);
    fy->Update(-1.0, *temp2, 1.0);
  }

  RangeExtractor().InsertVector(*sy,0,y);
  RangeExtractor().InsertVector(*fy,1,y);
  RangeExtractor().InsertVector(*ay,2,y);
  RangeExtractor().InsertVector(*cy,3,y);
}



Teuchos::RCP<LINALG::SparseMatrix> FSI::LungSchurComplement::CalculateSchur(const LINALG::SparseMatrix& A,
                                                                            const LINALG::SparseMatrix& B,
                                                                            const LINALG::SparseMatrix& C)
{
  // make sure FillComplete was called on the matrices
  if (!A.Filled()) dserror("A has to be FillComplete");
  if (!B.Filled()) dserror("B has to be FillComplete");
  if (!C.Filled()) dserror("C has to be FillComplete");

  if (temp_ == Teuchos::null)
  {
    const int npr = max(A.MaxNumEntries(),B.MaxNumEntries());
    temp_ = Teuchos::rcp(new LINALG::SparseMatrix(A.RangeMap(),npr,false,true));
  }
  int err = EpetraExt::MatrixMatrix::Multiply(*A.EpetraMatrix(),false,*B.EpetraMatrix(),false,*(temp_->EpetraMatrix()),true);
  if (err) dserror("EpetraExt::MatrixMatrix::Multiply returned err = %d",err);

  if (res_ == Teuchos::null)
  {
    const int npr = max(temp_->MaxNumEntries(),C.MaxNumEntries());
    res_ = Teuchos::rcp(new LINALG::SparseMatrix(temp_->RangeMap(),npr,false,true));
  }
  err = EpetraExt::MatrixMatrix::Multiply(*temp_->EpetraMatrix(),false,*C.EpetraMatrix(),false,*(res_->EpetraMatrix()),true);
  if (err) dserror("EpetraExt::MatrixMatrix::Multiply returned err = %d",err);

  return res_;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
// void FSI::LungOverlappingBlockMatrix::MyMatrixMatrixMultiply(const LINALG::SparseMatrix& A,
//                                                              const LINALG::SparseMatrix& B,
//                                                              Teuchos::RCP<LINALG::SparseMatrix> result) const
// {
//   // make sure FillComplete was called on the matrices
//   if (!A.Filled()) dserror("A has to be FillComplete");
//   if (!B.Filled()) dserror("B has to be FillComplete");

//   if (result == Teuchos::null)
//   {
//     const int npr = max(A.MaxNumEntries(),B.MaxNumEntries());
//     result = Teuchos::rcp(new LINALG::SparseMatrix(A.RangeMap(),npr,false,true));
//   }
//   int err = EpetraExt::MatrixMatrix::Multiply(*A.EpetraMatrix(),false,*B.EpetraMatrix(),false,*(result->EpetraMatrix()),true);
//   if (err) dserror("EpetraExt::MatrixMatrix::Multiply returned err = %d",err);
// }


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
const char* FSI::LungOverlappingBlockMatrix::Label() const
{
  return "FSI::LungOverlappingBlockMatrix";
}

// /*----------------------------------------------------------------------*
//  *----------------------------------------------------------------------*/
FSI::ConstrOverlappingBlockMatrix::ConstrOverlappingBlockMatrix(const LINALG::MultiMapExtractor& maps,
                                                            ADAPTER::Structure& structure,
                                                            ADAPTER::Fluid& fluid,
                                                            ALE::Ale& ale,
                                                            bool structuresplit,
                                                            int symmetric,
                                                            double omega,
                                                            int iterations,
                                                            double somega,
                                                            int siterations,
                                                            double fomega,
                                                            int fiterations,
                                                            double aomega,
                                                            int aiterations,
                                                            FILE* err)
  : OverlappingBlockMatrix(Teuchos::null,
                           maps,
                           structure,
                           fluid,
                           ale,
                           structuresplit,
                           symmetric,
                           omega,
                           iterations,
                           somega,
                           siterations,
                           fomega,
                           fiterations,
                           aomega,
                           aiterations,
                           err)
{
  // determine map of all dofs not related to constraint

  std::vector<Teuchos::RCP<const Epetra_Map> > fsimaps;
  fsimaps.push_back(maps.Map(0));
  fsimaps.push_back(maps.Map(1));
  fsimaps.push_back(maps.Map(2));
  overallfsimap_ = LINALG::MultiMapExtractor::MergeMaps(fsimaps);
  fsiextractor_ = LINALG::MultiMapExtractor(*overallfsimap_, fsimaps);

  // stuff needed for SIMPLE preconditioner -> this needs to be read
  // in from the input file one day!
  const Teuchos::ParameterList& fsidyn   = DRT::Problem::Instance()->FSIDynamicParams();
  alpha_ = fsidyn.sublist("CONSTRAINT").get<double>("ALPHA");
  simpleiter_ = fsidyn.sublist("CONSTRAINT").get<int>("SIMPLEITER");
  prec_ = Teuchos::getIntegralValue<INPAR::FSI::PrecConstr>(fsidyn.sublist("CONSTRAINT"),"PRECONDITIONER");
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void FSI::ConstrOverlappingBlockMatrix::SetupPreconditioner()
{
#ifdef BLOCKMATRIXMERGE

  // this is really evil
  // note that this only works with UMFPACK as fluid solver (saddle
  // point problem!)

  sparse_ = Merge();
  fluidsolver_->Setup(sparse_->EpetraMatrix());

#else

  const LINALG::SparseMatrix& structInnerOp = Matrix(0,0);
  const LINALG::SparseMatrix& fluidInnerOp  = Matrix(1,1);
  const LINALG::SparseMatrix& aleInnerOp    = Matrix(2,2);

  RCP<LINALG::MapExtractor> fsidofmapex = null;
  RCP<Epetra_Map>           irownodes = null;

  structuresolver_->Setup(structInnerOp.EpetraMatrix());
  fluidsolver_->Setup(fluidInnerOp.EpetraMatrix(),
                      fsidofmapex,
                      fluid_.Discretization(),
                      irownodes,
                      structuresplit_);
  if (constalesolver_==Teuchos::null)
    alesolver_->Setup(aleInnerOp.EpetraMatrix());

#endif
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void FSI::ConstrOverlappingBlockMatrix::SGS(const Epetra_MultiVector &X, Epetra_MultiVector &Y) const
{
  const LINALG::SparseMatrix& StructInnerOp  = Matrix(0,0);
  const LINALG::SparseMatrix& StructBoundOp  = Matrix(0,1);
  const LINALG::SparseMatrix& StructConOp    = Matrix(0,3);
  const LINALG::SparseMatrix& FluidBoundOp   = Matrix(1,0);
  const LINALG::SparseMatrix& FluidInnerOp   = Matrix(1,1);
  const LINALG::SparseMatrix& FluidMeshOp    = Matrix(1,2);
  const LINALG::SparseMatrix& FluidConOp     = Matrix(1,3);
  const LINALG::SparseMatrix& AleInnerOp     = Matrix(2,2);
  const LINALG::SparseMatrix& ConStructOp    = Matrix(3,0);
  const LINALG::SparseMatrix& ConFluidOp     = Matrix(3,1);


  // Extract vector blocks
  // RHS

  const Epetra_Vector &x = Teuchos::dyn_cast<const Epetra_Vector>(X);

  // initial guess

  Epetra_Vector &y = Teuchos::dyn_cast<Epetra_Vector>(Y);

  Teuchos::RCP<Epetra_Vector> sy = RangeExtractor().ExtractVector(y,0);
  Teuchos::RCP<Epetra_Vector> fy = RangeExtractor().ExtractVector(y,1);
  Teuchos::RCP<Epetra_Vector> ay = RangeExtractor().ExtractVector(y,2);
  Teuchos::RCP<Epetra_Vector> cy = RangeExtractor().ExtractVector(y,3);

  Teuchos::RCP<Epetra_Vector> sz = Teuchos::rcp(new Epetra_Vector(sy->Map()));
  Teuchos::RCP<Epetra_Vector> fz = Teuchos::rcp(new Epetra_Vector(fy->Map()));
  Teuchos::RCP<Epetra_Vector> az = Teuchos::rcp(new Epetra_Vector(ay->Map()));

  Teuchos::RCP<Epetra_Vector> tmpsx = Teuchos::rcp(new Epetra_Vector(DomainMap(0)));
  Teuchos::RCP<Epetra_Vector> tmpfx = Teuchos::rcp(new Epetra_Vector(DomainMap(1)));
  Teuchos::RCP<Epetra_Vector> tmpax = Teuchos::rcp(new Epetra_Vector(DomainMap(2)));

  for (int outerrun=0; outerrun<simpleiter_; ++outerrun)
  {

    // -------------------------------------------------------------------
    // intermediate fsi dofs: u_(n+1/2) = F^(-1)(f-B^T*lambda_n)
    // -------------------------------------------------------------------

    // inner Richardson loop (FSI block)

    for (int run=0; run<iterations_; ++run)
    {
      Teuchos::RCP<Epetra_Vector> sx = DomainExtractor().ExtractVector(x,0);
      Teuchos::RCP<Epetra_Vector> fx = DomainExtractor().ExtractVector(x,1);
      Teuchos::RCP<Epetra_Vector> ax = DomainExtractor().ExtractVector(x,2);

      // ----------------------------------------------------------------
      // lower GS


      // Structure
      {

        if (run>0 or outerrun>0)
        {
          StructConOp.Multiply(false,*cy,*tmpsx);
          sx->Update(-1.0, *tmpsx, 1.0);

          StructInnerOp.Multiply(false,*sy,*tmpsx);
          sx->Update(-1.0,*tmpsx,1.0);
          StructBoundOp.Multiply(false,*fy,*tmpsx);
          sx->Update(-1.0,*tmpsx,1.0);
        }

        // Solve structure equations for sy with the rhs sx
        structuresolver_->Solve(StructInnerOp.EpetraMatrix(),sz,sx,true);

        // do Richardson iteration
        LocalBlockRichardson(structuresolver_,StructInnerOp,sx,sz,tmpsx,siterations_,somega_,err_,Comm());

        if (run>0 or outerrun>0)
        {
          sy->Update(omega_,*sz,1.0);
        }
        else
        {
          sy->Update(omega_,*sz,0.0);
        }
      }

      // Ale
      {
        // Solve ale equations for ay with the rhs ax - A(I,Gamma) sy

        if (run>0 or outerrun>0)
        {
          AleInnerOp.Multiply(false,*ay,*tmpax);
          ax->Update(-1.0,*tmpax,1.0);
        }
        if (structuresplit_)
        {
          if (run>0 or outerrun>0)
          {
            const LINALG::SparseMatrix& aleBoundOp = Matrix(2,1);
            aleBoundOp.Multiply(false,*fy,*tmpax);
            ax->Update(-1.0,*tmpax,1.0);
          }

        }
        else
        {
          const LINALG::SparseMatrix& aleBoundOp   = Matrix(2,0);
          aleBoundOp.Multiply(false,*sy,*tmpax);
          ax->Update(-1.0,*tmpax,1.0);
        }

        alesolver_->Solve(AleInnerOp.EpetraMatrix(),az,ax,true);

        // do Richardson iteration
        LocalBlockRichardson(alesolver_,AleInnerOp,ax,az,tmpax,aiterations_,aomega_,err_,Comm());

        if (run>0 or outerrun>0)
        {
          ay->Update(omega_,*az,1.0);
        }
        else
        {
          ay->Update(omega_,*az,0.0);
        }
      }

      // Fluid
      {
        // Solve fluid equations for fy with the rhs fx - F(I,Gamma) sy - F(Mesh) ay - F(Constr) cy

        if (run>0 or outerrun>0)
        {
          FluidInnerOp.Multiply(false,*fy,*tmpfx);
          fx->Update(-1.0,*tmpfx,1.0);
        }

        FluidBoundOp.Multiply(false,*sy,*tmpfx);
        fx->Update(-1.0,*tmpfx,1.0);
        FluidMeshOp.Multiply(false,*ay,*tmpfx);
        fx->Update(-1.0,*tmpfx,1.0);
        FluidConOp.Multiply(false,*cy, *tmpfx);
        fx->Update(-1.0, *tmpfx, 1.0);

        fluidsolver_->Solve(FluidInnerOp.EpetraMatrix(),fz,fx,true);

        LocalBlockRichardson(fluidsolver_,FluidInnerOp,fx,fz,tmpfx,fiterations_,fomega_,err_,Comm());

        if (run>0 or outerrun>0)
        {
          fy->Update(omega_,*fz,1.0);
        }
        else
        {
          fy->Update(omega_,*fz,0.0);
        }
      }

      // ----------------------------------------------------------------
      // the symmetric part of the pc can be skipped

      if (symmetric_)
      {
        sx = DomainExtractor().ExtractVector(x,0);
        fx = DomainExtractor().ExtractVector(x,1);
        ax = DomainExtractor().ExtractVector(x,2);

        // ----------------------------------------------------------------
        // upper GS

        {
          FluidInnerOp.Multiply(false,*fy,*tmpfx);
          fx->Update(-1.0,*tmpfx,1.0);
          FluidBoundOp.Multiply(false,*sy,*tmpfx);
          fx->Update(-1.0,*tmpfx,1.0);
          FluidMeshOp.Multiply(false,*ay,*tmpfx);
          fx->Update(-1.0,*tmpfx,1.0);
          FluidConOp.Multiply(false,*cy,*tmpfx);
          fx->Update(-1.0, *tmpfx, 1.0);

          fluidsolver_->Solve(FluidInnerOp.EpetraMatrix(),fz,fx,true);

          LocalBlockRichardson(fluidsolver_,FluidInnerOp,fx,fz,tmpfx,fiterations_,fomega_,err_,Comm());
          fy->Update(omega_,*fz,1.0);
        }

        {
          AleInnerOp.Multiply(false,*ay,*tmpax);
          ax->Update(-1.0,*tmpax,1.0);

          if (structuresplit_)
          {
            const LINALG::SparseMatrix& aleBoundOp = Matrix(2,1);
            aleBoundOp.Multiply(false,*fy,*tmpax);
            ax->Update(-1.0,*tmpax,1.0);

          }
          else
          {
            const LINALG::SparseMatrix& aleBoundOp = Matrix(2,0);
            aleBoundOp.Multiply(false,*sy,*tmpax);
            ax->Update(-1.0,*tmpax,1.0);
          }

          alesolver_->Solve(AleInnerOp.EpetraMatrix(),az,ax,true);

          LocalBlockRichardson(alesolver_,AleInnerOp,ax,az,tmpax,aiterations_,aomega_,err_,Comm());

          ay->Update(omega_,*az,1.0);
        }

        {
          StructInnerOp.Multiply(false,*sy,*tmpsx);
          sx->Update(-1.0,*tmpsx,1.0);
          StructBoundOp.Multiply(false,*fy,*tmpsx);
          sx->Update(-1.0,*tmpsx,1.0);
          StructConOp.Multiply(false,*cy,*tmpsx);
          sx->Update(-1.0, *tmpsx, 1.0);


          structuresolver_->Solve(StructInnerOp.EpetraMatrix(),sz,sx,true);

          LocalBlockRichardson(structuresolver_,StructInnerOp,sx,sz,tmpsx,siterations_,somega_,err_,Comm());
          sy->Update(omega_,*sz,1.0);
        }
      }
    }

    // -----------------------------------------------------------------------
    // intermediate constraint dofs: Dlambda~ = S^(-1) * (cx - B^ * u_(n+1/2))
    // -----------------------------------------------------------------------

    // cx - B^ * u_(n+1/2)

    Teuchos::RCP<Epetra_Vector> cx = DomainExtractor().ExtractVector(x,3);

    Epetra_Vector inter(cx->Map());
    ConStructOp.Multiply(false,*sy, inter);
    cx->Update(-1.0, inter, 1.0);
    ConFluidOp.Multiply(false,*fy, inter);
    cx->Update(-1.0, inter, 1.0);


    LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy> invDiag(fsiextractor_, fsiextractor_, 1, false, true);

    if (prec_ == INPAR::FSI::Simple)
    {
      // D^{-1} = diag(A(0,0))^{-1}

      Epetra_Vector structDiagVec(StructInnerOp.RowMap(),false);
      StructInnerOp.ExtractDiagonalCopy(structDiagVec);
      int err = structDiagVec.Reciprocal(structDiagVec);
      if (err) dserror("Epetra_MultiVector::Reciprocal (structure matrix) returned %d",err);
      LINALG::SparseMatrix invstructDiag(structDiagVec);

      Epetra_Vector fluidDiagVec(FluidInnerOp.RowMap(),false);
      FluidInnerOp.ExtractDiagonalCopy(fluidDiagVec);
      err = fluidDiagVec.Reciprocal(fluidDiagVec);
      if (err) dserror("Epetra_MultiVector::Reciprocal (fluid matrix) returned %d",err);
      LINALG::SparseMatrix invfluidDiag(fluidDiagVec);

      Epetra_Vector aleDiagVec(AleInnerOp.RowMap(),false);
      AleInnerOp.ExtractDiagonalCopy(aleDiagVec);
      err = aleDiagVec.Reciprocal(aleDiagVec);
      if (err) dserror("Epetra_MultiVector::Reciprocal (ale matrix) returned %d",err);
      LINALG::SparseMatrix invaleDiag(aleDiagVec);

      invDiag.Assign(0,0,View,invstructDiag);
      invDiag.Assign(1,1,View,invfluidDiag);
      invDiag.Assign(2,2,View,invaleDiag);
    }
    else if (prec_ == INPAR::FSI::Simplec)
    {
      // D^{-1} = sum(abs(A(0,0)))^{-1}

      Epetra_Vector structDiagVec(StructInnerOp.RowMap(),false);
      int err = StructInnerOp.EpetraMatrix()->InvRowSums(structDiagVec);
      if (err) dserror("Epetra_CrsMatrix::InvRowSums (structure matrix) returned %d",err);
      LINALG::SparseMatrix invstructDiag(structDiagVec);

      Epetra_Vector fluidDiagVec(FluidInnerOp.RowMap(),false);
      err = FluidInnerOp.EpetraMatrix()->InvRowSums(fluidDiagVec);
      if (err) dserror("Epetra_CrsMatrix::InvRowSums (fluid matrix) returned %d",err);
      LINALG::SparseMatrix invfluidDiag(fluidDiagVec);

      Epetra_Vector aleDiagVec(AleInnerOp.RowMap(),false);
      err =  AleInnerOp.EpetraMatrix()->InvRowSums(aleDiagVec);
      if (err) dserror("Epetra_CrsMatrix::InvRowSums (ale matrix) returned %d",err);
      LINALG::SparseMatrix invaleDiag(aleDiagVec);

      invDiag.Assign(0,0,View,invstructDiag);
      invDiag.Assign(1,1,View,invfluidDiag);
      invDiag.Assign(2,2,View,invaleDiag);
    }
    else
      dserror("Unknown type of preconditioner for constraint fsi system");

    invDiag.Complete();


    // S = - B^ * D^{-1} * B^T

    Teuchos::RCP<LINALG::SparseMatrix> temps = LINALG::Multiply(ConStructOp, false, invDiag.Matrix(0,0), false);
    Teuchos::RCP<LINALG::SparseMatrix> interconA = LINALG::Multiply(*temps, false, StructConOp, false, false);

    temps = LINALG::Multiply(ConFluidOp, false, invDiag.Matrix(1,1), false);
    Teuchos::RCP<LINALG::SparseMatrix> tempss = LINALG::Multiply(*temps, false, FluidConOp, false);
    interconA->Add(*tempss, false, 1.0, 1.0);
    interconA->Complete(StructConOp.DomainMap(),ConStructOp.RangeMap());
    interconA->Scale(-1.0);

    RCP<Teuchos::ParameterList> constrsolvparams = rcp(new Teuchos::ParameterList);
    constrsolvparams->set("solver","umfpack");
    RCP<Epetra_Vector> interconsol = rcp(new Epetra_Vector(ConStructOp.RangeMap()));
    RCP<LINALG::Solver> ConstraintSolver = rcp(new LINALG::Solver(constrsolvparams,
                                                                  interconA->Comm(),
                                                                  DRT::Problem::Instance()->ErrorFile()->Handle()));
    ConstraintSolver->Solve(interconA->EpetraOperator(),interconsol,cx,true,true);

    // -------------------------------------------------------------------
    // update of all dofs
    // -------------------------------------------------------------------

    cy->Update(alpha_, *interconsol, 1.0);

    Teuchos::RCP<Epetra_Vector> temp1;
    Teuchos::RCP<Epetra_Vector> temp2;

    temp1 = rcp(new Epetra_Vector(sy->Map()));
    temp2 = rcp(new Epetra_Vector(sy->Map()));
    StructConOp.Multiply(false,*interconsol, *temp1);
    temp1->Scale(alpha_);
    invDiag.Matrix(0,0).Multiply(false,*temp1, *temp2);
    sy->Update(-1.0, *temp2, 1.0);

    temp1 = rcp(new Epetra_Vector(fy->Map()));
    temp2 = rcp(new Epetra_Vector(fy->Map()));
    FluidConOp.Multiply(false,*interconsol, *temp1);
    temp1->Scale(alpha_);
    invDiag.Matrix(1,1).Multiply(false,*temp1, *temp2);
    fy->Update(-1.0, *temp2, 1.0);
  }

  RangeExtractor().InsertVector(*sy,0,y);
  RangeExtractor().InsertVector(*fy,1,y);
  RangeExtractor().InsertVector(*ay,2,y);
  RangeExtractor().InsertVector(*cy,3,y);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
const char* FSI::ConstrOverlappingBlockMatrix::Label() const
{
  return "FSI::ConstrOverlappingBlockMatrix";
}

#endif
