#ifdef CCADISCRET

#include "fsi_overlapprec_fsiamg.H"
#include <Epetra_Time.h>
#include <ml_MultiLevelPreconditioner.h>
#include "MLAPI_Operator_Utils.h"
//#include "MLAPI_Error.h"
//#include "MLAPI_CompObject.h"
//#include "MLAPI_TimeObject.h"
//#include "MLAPI_MultiVector.h"
//#include "MLAPI_Expressions.h"
//#include "MLAPI_BaseOperator.h"
//#include "MLAPI_Workspace.h"
//#include "MLAPI_Aggregation.h"
//#include "MLAPI_Eig.h"



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
FSI::OverlappingBlockMatrixFSIAMG::OverlappingBlockMatrixFSIAMG(
                                                    const LINALG::MultiMapExtractor& maps,
                                                    ADAPTER::Structure& structure,
                                                    ADAPTER::Fluid& fluid,
                                                    ADAPTER::Ale& ale,
                                                    bool structuresplit,
                                                    int symmetric,
                                                    double omega,
                                                    int iterations,
                                                    double somega,
                                                    int siterations,
                                                    double fomega,
                                                    int fiterations,
                                                    FILE* err)
  : OverlappingBlockMatrix(maps,
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
                           err)
{
  MLAPI::Init();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void FSI::OverlappingBlockMatrixFSIAMG::SetupPreconditioner()
{
  const int myrank = Matrix(0,0).Comm().MyPID();

  const LINALG::SparseMatrix& structInnerOp = Matrix(0,0);
  const LINALG::SparseMatrix& fluidInnerOp  = Matrix(1,1);
  const LINALG::SparseMatrix& aleInnerOp    = Matrix(2,2);
  
  RCP<LINALG::MapExtractor> fsidofmapex = null;
  RCP<Epetra_Map>           irownodes = null;
  
  // build AMG hierarchies
  structuresolver_->Setup(structInnerOp.EpetraMatrix());
  fluidsolver_->Setup(fluidInnerOp.EpetraMatrix(),
                      fsidofmapex,
                      fluid_.Discretization(),
                      irownodes,
                      structuresplit_);
  if (constalesolver_==Teuchos::null)
    alesolver_->Setup(aleInnerOp.EpetraMatrix());

  // get the ml_MultiLevelPreconditioner class from within struct/fluid solver
  RCP<Epetra_Operator> sprec = structuresolver_->EpetraOperator();
  RCP<Epetra_Operator> fprec = fluidsolver_->EpetraOperator();
  RCP<Epetra_Operator> aprec = alesolver_->EpetraOperator();
  
  // get ML preconditioner class
  ML_Epetra::MultiLevelPreconditioner* smlclass = 
    dynamic_cast<ML_Epetra::MultiLevelPreconditioner*>(sprec.get());
  ML_Epetra::MultiLevelPreconditioner* fmlclass = 
    dynamic_cast<ML_Epetra::MultiLevelPreconditioner*>(fprec.get());
  ML_Epetra::MultiLevelPreconditioner* amlclass = 
    dynamic_cast<ML_Epetra::MultiLevelPreconditioner*>(aprec.get());
  if (!smlclass || !fmlclass || !amlclass) dserror("Not using ML for Fluid, Structure or Ale");

  // get copy of ML parameter list
  sparams_ = structuresolver_->Params().sublist("ML Parameters");
  fparams_ = fluidsolver_->Params().sublist("ML Parameters");
  aparams_ = alesolver_->Params().sublist("ML Parameters");
  //cout << "====================Structure Params\n" << sparams_;
  //cout << "========================Fluid Params\n" << fparams_;
  //cout << "========================Ale Params\n" << aparams_;
  // find for which field using nonsymmetric AMG
  bool SisPetrovGalerkin = sparams_.get("energy minimization: enable",false);
  bool FisPetrovGalerkin = fparams_.get("energy minimization: enable",false);
  bool AisPetrovGalerkin = aparams_.get("energy minimization: enable",false);
  if (AisPetrovGalerkin) dserror("MLFLUID2 hierarchy not supported for Ale, choose ML");

  // get ML handle
  ML* sml = const_cast<ML*>(smlclass->GetML());
  ML* fml = const_cast<ML*>(fmlclass->GetML());
  ML* aml = const_cast<ML*>(amlclass->GetML());
  if (!sml || !fml || !aml) dserror("Not using ML for Fluid, Structure or Ale");

  // number of grids in structure, fluid and ale
  snlevel_ = sml->ML_num_actual_levels;
  fnlevel_ = fml->ML_num_actual_levels;
  anlevel_ = aml->ML_num_actual_levels;
  // min of number of grids over fields
  minnlevel_ = min(snlevel_,fnlevel_);
  minnlevel_ = min(minnlevel_,anlevel_);
  
  if (!myrank)
  {
    printf("Setting up FSIAMG: snlevel %d fnlevel %d anlevel %d minnlevel %d\n",
           snlevel_,fnlevel_,anlevel_,minnlevel_);
  }

  Ass_.resize(snlevel_);
  Pss_.resize(snlevel_-1);
  Rss_.resize(snlevel_-1);
  Sss_.resize(snlevel_);

  Aff_.resize(fnlevel_);
  Pff_.resize(fnlevel_-1);
  Rff_.resize(fnlevel_-1);
  Sff_.resize(fnlevel_);
  
  Aaa_.resize(anlevel_);
  Paa_.resize(anlevel_-1);
  Saa_.resize(anlevel_);

  //------------------------------------------------------- Structure
  {
    // fine space matching Epetra objects
    MLAPI::Space finespace(structInnerOp.RowMap());
    if (!myrank) printf("Structure: Epetra NumMyElements fine space %d\n",structInnerOp.RowMap().NumMyElements());
    
    // extract transfer operator P,R from ML
    MLAPI::Space fspace;
    MLAPI::Space cspace;
    MLAPI::Operator P;
    MLAPI::Operator R;
    for (int i=1; i<snlevel_; ++i)
    {
      ML_Operator* Pml = &(sml->Pmat[i]);
      if (!myrank) printf("Structure P: level %d domain %d range %d\n",i,Pml->invec_leng,Pml->outvec_leng);
      if (i==1) fspace = finespace;
      else      fspace.Reshape(-1,Pml->outvec_leng);
      cspace.Reshape(-1,Pml->invec_leng);
      P.Reshape(cspace,fspace,Pml,false);
      if (SisPetrovGalerkin)
      {
        ML_Operator* Rml = &(sml->Rmat[i-1]);
        if (!myrank) printf("Structure R: level %d domain %d range %d\n",i,Rml->invec_leng,Rml->outvec_leng);
        R.Reshape(fspace,cspace,Rml,false);
      }
      else 
      {
        R = MLAPI::GetTranspose(P);
      }
      Pss_[i-1] = P;
      Rss_[i-1] = R;
    }
    // extract matrix A from ML
    MLAPI::Space space;
    for (int i=0; i<snlevel_; ++i)
    {
      ML_Operator* Aml = &(sml->Amat[i]);
      if (i==0) space = finespace;
      else      space.Reshape(-1,Aml->invec_leng);
      MLAPI::Operator A(space,space,Aml,false);
      Ass_[i] = A;
    }
  }

  //------------------------------------------------------- Fluid
  {
    // fine space matching Epetra objects
    MLAPI::Space finespace(fluidInnerOp.RowMap());
    if (!myrank) printf("Fluid: Epetra NumMyElements fine space %d\n",fluidInnerOp.RowMap().NumMyElements());
    
    // extract transfer operator P,R from ML
    MLAPI::Space fspace;
    MLAPI::Space cspace;
    MLAPI::Operator P;
    MLAPI::Operator R;
    for (int i=1; i<fnlevel_; ++i)
    {
      ML_Operator* Pml = &(fml->Pmat[i]);
      if (!myrank) printf("Fluid P: level %d domain %d range %d\n",i,Pml->invec_leng,Pml->outvec_leng);
      if (i==1) fspace = finespace;
      else      fspace.Reshape(-1,Pml->outvec_leng);
      cspace.Reshape(-1,Pml->invec_leng);
      P.Reshape(cspace,fspace,Pml,false);
      if (FisPetrovGalerkin)
      {
        ML_Operator* Rml = &(fml->Rmat[i-1]);
        if (!myrank) printf("Fluid R: level %d domain %d range %d\n",i,Rml->invec_leng,Rml->outvec_leng);
        R.Reshape(fspace,cspace,Rml,false);
      }
      else 
      {
        R = MLAPI::GetTranspose(P);
      }
      Pff_[i-1] = P;
      Rff_[i-1] = R;
    }
    // extract matrix A from ML
    MLAPI::Space space;
    for (int i=0; i<fnlevel_; ++i)
    {
      ML_Operator* Aml = &(fml->Amat[i]);
      if (i==0) space = finespace;
      else      space.Reshape(-1,Aml->invec_leng);
      MLAPI::Operator A(space,space,Aml,false);
      Aff_[i] = A;
    }
  }

  //------------------------------------------------------- Ale
  {
    // fine space matching Epetra objects
    MLAPI::Space finespace(aleInnerOp.RowMap());
    if (!myrank) printf("Ale: Epetra NumMyElements fine space %d\n",structInnerOp.RowMap().NumMyElements());
    
    // extract transfer operator P,R from ML
    MLAPI::Space fspace;
    MLAPI::Space cspace;
    MLAPI::Operator P;
    for (int i=1; i<anlevel_; ++i)
    {
      ML_Operator* Pml = &(aml->Pmat[i]);
      if (!myrank) printf("Ale P: level %d domain %d range %d\n",i,Pml->invec_leng,Pml->outvec_leng);
      if (i==1) fspace = finespace;
      else      fspace.Reshape(-1,Pml->outvec_leng);
      cspace.Reshape(-1,Pml->invec_leng);
      P.Reshape(cspace,fspace,Pml,false);
      Paa_[i-1] = P;
    }
    // extract matrix A from ML
    MLAPI::Space space;
    for (int i=0; i<anlevel_; ++i)
    {
      ML_Operator* Aml = &(aml->Amat[i]);
      if (i==0) space = finespace;
      else      space.Reshape(-1,Aml->invec_leng);
      MLAPI::Operator A(space,space,Aml,false);
      Aaa_[i] = A;
    }
  }

  //================set up MLAPI smoothers for structure, fluid, ale on each level
  MLAPI::InverseOperator S;
  // structure
  for (int i=0; i<snlevel_-1; ++i)
  {
    Teuchos::ParameterList p;
    char levelstr[11];
    sprintf(levelstr,"(level %d)",i);
    Teuchos::ParameterList& subp = sparams_.sublist("smoother: list "+(string)levelstr);
    //cout << "*********************level " << i << " sublist:\n" << subp;
    string type = "";
    SelectMLAPISmoother(type,subp,p);
    S.Reshape(Ass_[i],type,p);
    Sss_[i] = S;
  }
  // coarse grid:
  S.Reshape(Ass_[snlevel_-1],"Amesos-KLU");
  Sss_[snlevel_-1] = S;
  // fluid
  for (int i=0; i<fnlevel_-1; ++i)
  {
    Teuchos::ParameterList p;
    char levelstr[11];
    sprintf(levelstr,"(level %d)",i);
    Teuchos::ParameterList& subp = fparams_.sublist("smoother: list "+(string)levelstr);
    //cout << "*********************level " << i << " sublist:\n" << subp;
    string type = "";
    SelectMLAPISmoother(type,subp,p);
    S.Reshape(Ass_[i],type,p);
    Sff_[i] = S;
  }
  // coarse grid:
  S.Reshape(Aff_[fnlevel_-1],"Amesos-KLU");
  Sff_[fnlevel_-1] = S;
  // ale
  for (int i=0; i<anlevel_-1; ++i)
  {
    Teuchos::ParameterList p;
    char levelstr[11];
    sprintf(levelstr,"(level %d)",i);
    Teuchos::ParameterList& subp = aparams_.sublist("smoother: list "+(string)levelstr);
    //cout << "*********************level " << i << " sublist:\n" << subp;
    string type = "";
    SelectMLAPISmoother(type,subp,p);
    S.Reshape(Aaa_[i],type,p);
    Saa_[i] = S;
  }
  // coarse grid:
  S.Reshape(Aaa_[anlevel_-1],"Amesos-KLU");
  Saa_[anlevel_-1] = S;
  

  // wrap the off-diagonal matrix blocks into MLAPI operators
  const LINALG::SparseMatrix& structBoundOp = Matrix(0,1); // Asf
  const LINALG::SparseMatrix& fluidBoundOp  = Matrix(1,0); // Afs
  const LINALG::SparseMatrix& fluidMeshOp   = Matrix(1,2); // Afa
  const LINALG::SparseMatrix& aleBoundOp    = Matrix(2,1); // Aaf

  MLAPI::Operator  Asf;
  {
    MLAPI::Space dspace(structBoundOp.EpetraMatrix()->DomainMap());
    MLAPI::Space rspace(structBoundOp.EpetraMatrix()->RangeMap());
    Asf.Reshape(dspace,rspace,structBoundOp.EpetraMatrix().get(),false);
  }
  MLAPI::Operator  Afs;
  {
    MLAPI::Space dspace(fluidBoundOp.EpetraMatrix()->DomainMap());
    MLAPI::Space rspace(fluidBoundOp.EpetraMatrix()->RangeMap());
    Afs.Reshape(dspace,rspace,fluidBoundOp.EpetraMatrix().get(),false);
  }
  MLAPI::Operator  Afa;
  {
    MLAPI::Space dspace(fluidMeshOp.EpetraMatrix()->DomainMap());
    MLAPI::Space rspace(fluidMeshOp.EpetraMatrix()->RangeMap());
    Afa.Reshape(dspace,rspace,fluidMeshOp.EpetraMatrix().get(),false);
  }
  MLAPI::Operator  Aaf;
  {
    MLAPI::Space dspace(aleBoundOp.EpetraMatrix()->DomainMap());
    MLAPI::Space rspace(aleBoundOp.EpetraMatrix()->RangeMap());
    Aaf.Reshape(dspace,rspace,aleBoundOp.EpetraMatrix().get(),false);
  }
  exit(0);
  

  return;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void FSI::OverlappingBlockMatrixFSIAMG::SelectMLAPISmoother(std::string& type,
                                                            Teuchos::ParameterList& subp,
                                                            Teuchos::ParameterList& p)
{
  type = subp.get("smoother: type","none");
  if (type=="none") dserror("Cannot find msoother type");
  if (type=="symmetric Gauss-Seidel" || type=="Gauss-Seidel")
  {
    const int    sweeps  = subp.get("smoother: sweeps",1);
    const double damping = subp.get("smoother: damping factor",1.0);
    p.set("smoother: sweeps",sweeps);
    p.set("smoother: damping factor",damping); 
  }
  else if (type=="IFPACK")
  {
    type              = subp.get("smoother: ifpack type","ILU");
    const int lof     = subp.get("smoother: ifpack level-of-fill",0);
    p.set("fact: level-of-fill",lof);
  }
  else if (type=="MLS")
  {
    const int poly = subp.get("smoother: MLS polynomial order",3);
    p.set("smoother: MLS polynomial order",poly);
  }
  else if (type=="Amesos-KLU"); // nothing to do
  else if (type=="Amesos-Superludist") dserror("No SuperLUDist support in MLAPI");
  else dserror("Smoother not recognized");
  return;
}                                                            


#if 1
/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void FSI::OverlappingBlockMatrixFSIAMG::SGS(
                 const Epetra_MultiVector &X, Epetra_MultiVector &Y) const
{
  if (!structuresplit_) dserror("FSIAMG for structuresplit monoFSI only");
  if (symmetric_)       dserror("FSIAMG symmetric Block Gauss-Seidel not impl.");
  
  const LINALG::SparseMatrix& structInnerOp = Matrix(0,0); // Ass
  const LINALG::SparseMatrix& structBoundOp = Matrix(0,1); // Asf
  const LINALG::SparseMatrix& fluidBoundOp  = Matrix(1,0); // Afs
  const LINALG::SparseMatrix& fluidInnerOp  = Matrix(1,1); // Aff
  const LINALG::SparseMatrix& fluidMeshOp   = Matrix(1,2); // Afa
  const LINALG::SparseMatrix& aleBoundOp    = Matrix(2,1); // Aaf
  const LINALG::SparseMatrix& aleInnerOp    = Matrix(2,2); // Aaa

  const MLAPI::Operator& Ass = Ass_[0];
  const MLAPI::Operator& Aff = Aff_[0];
  const MLAPI::Operator& Aaa = Aaa_[0];

  // Extract vector blocks
  // RHS
  const Epetra_Vector &x = Teuchos::dyn_cast<const Epetra_Vector>(X);

  // initial guess
  Epetra_Vector& y = Teuchos::dyn_cast<Epetra_Vector>(Y);

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
    // rhs
    Teuchos::RCP<Epetra_Vector> sx = DomainExtractor().ExtractVector(x,0);
    Teuchos::RCP<Epetra_Vector> fx = DomainExtractor().ExtractVector(x,1);
    Teuchos::RCP<Epetra_Vector> ax = DomainExtractor().ExtractVector(x,2);

    // ----------------------------------------------------------------
    // lower GS
    {
      if (run>0)
      {
        structInnerOp.Multiply(false,*sy,*tmpsx);
        sx->Update(-1.0,*tmpsx,1.0);
        structBoundOp.Multiply(false,*fy,*tmpsx);
        sx->Update(-1.0,*tmpsx,1.0);
      }

      // Solve structure equations for sy with the rhs sx
      structuresolver_->Solve(structInnerOp.EpetraMatrix(),sz,sx,true);

      // do Richardson iteration
      LocalBlockRichardson(structuresolver_,structInnerOp,sx,sz,tmpsx,siterations_,somega_,err_,Comm());

      if (!run) sy->Update(omega_,*sz,0.0);
      else      sy->Update(omega_,*sz,1.0);
    }

    {
      // Solve ale equations for ay with the rhs ax - A(I,Gamma) sy
      if (run>0)
      {
        aleInnerOp.Multiply(false,*ay,*tmpax);
        ax->Update(-1.0,*tmpax,1.0);

        aleBoundOp.Multiply(false,*fy,*tmpax);
        ax->Update(-1.0,*tmpax,1.0);
      }

      alesolver_->Solve(aleInnerOp.EpetraMatrix(),az,ax,true);

      if (!run) ay->Update(omega_,*az,0.0);
      else      ay->Update(omega_,*az,1.0);
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

      if (!run) fy->Update(omega_,*fz,0.0);
      else      fy->Update(omega_,*fz,1.0);
    }

  }

  RangeExtractor().InsertVector(*sy,0,y);
  RangeExtractor().InsertVector(*fy,1,y);
  RangeExtractor().InsertVector(*ay,2,y);
}
#endif

#if 0
/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void FSI::OverlappingBlockMatrixFSIAMG::SGS(
                 const Epetra_MultiVector &X, Epetra_MultiVector &Y) const
{
  if (!structuresplit_) dserror("FSIAMG for structuresplit monoFSI only");
  if (symmetric_)       dserror("FSIAMG symmetric Block Gauss-Seidel not impl.");
  
  const LINALG::SparseMatrix& structInnerOp = Matrix(0,0); // Ass
  const LINALG::SparseMatrix& structBoundOp = Matrix(0,1); // Asf
  const LINALG::SparseMatrix& fluidBoundOp  = Matrix(1,0); // Afs
  const LINALG::SparseMatrix& fluidInnerOp  = Matrix(1,1); // Aff
  const LINALG::SparseMatrix& fluidMeshOp   = Matrix(1,2); // Afa
  const LINALG::SparseMatrix& aleBoundOp    = Matrix(2,1); // Aaf
  const LINALG::SparseMatrix& aleInnerOp    = Matrix(2,2); // Aaa


  // Extract vector blocks
  // RHS
  const Epetra_Vector &x = Teuchos::dyn_cast<const Epetra_Vector>(X);

  // initial guess
  Epetra_Vector& y = Teuchos::dyn_cast<Epetra_Vector>(Y);

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
    // rhs
    Teuchos::RCP<Epetra_Vector> sx = DomainExtractor().ExtractVector(x,0);
    Teuchos::RCP<Epetra_Vector> fx = DomainExtractor().ExtractVector(x,1);
    Teuchos::RCP<Epetra_Vector> ax = DomainExtractor().ExtractVector(x,2);

    // ----------------------------------------------------------------
    // lower GS
    {
      if (run>0)
      {
        structInnerOp.Multiply(false,*sy,*tmpsx);
        sx->Update(-1.0,*tmpsx,1.0);
        structBoundOp.Multiply(false,*fy,*tmpsx);
        sx->Update(-1.0,*tmpsx,1.0);
      }

      // Solve structure equations for sy with the rhs sx
      structuresolver_->Solve(structInnerOp.EpetraMatrix(),sz,sx,true);

      // do Richardson iteration
      LocalBlockRichardson(structuresolver_,structInnerOp,sx,sz,tmpsx,siterations_,somega_,err_,Comm());

      if (!run) sy->Update(omega_,*sz,0.0);
      else      sy->Update(omega_,*sz,1.0);
    }

    {
      // Solve ale equations for ay with the rhs ax - A(I,Gamma) sy
      if (run>0)
      {
        aleInnerOp.Multiply(false,*ay,*tmpax);
        ax->Update(-1.0,*tmpax,1.0);

        aleBoundOp.Multiply(false,*fy,*tmpax);
        ax->Update(-1.0,*tmpax,1.0);
      }

      alesolver_->Solve(aleInnerOp.EpetraMatrix(),az,ax,true);

      if (!run) ay->Update(omega_,*az,0.0);
      else      ay->Update(omega_,*az,1.0);
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

      if (!run) fy->Update(omega_,*fz,0.0);
      else      fy->Update(omega_,*fz,1.0);
    }

  }

  RangeExtractor().InsertVector(*sy,0,y);
  RangeExtractor().InsertVector(*fy,1,y);
  RangeExtractor().InsertVector(*ay,2,y);
}
#endif

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
const char* FSI::OverlappingBlockMatrixFSIAMG::Label() const
{
  return "FSI::OverlappingBlockMatrix_FSIAMG";
}

#endif
