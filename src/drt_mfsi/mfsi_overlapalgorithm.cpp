/*----------------------------------------------------------------------*/
/*!
\file mfsi_overlapalgorithm.cpp

\brief Solve monolithic overlapping block coupled FSI problems

<pre>
Maintainer: Ulrich Kuettler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15238
</pre>
*/
/*----------------------------------------------------------------------*/

#ifdef CCADISCRET

#include "mfsi_overlapalgorithm.H"

#include "mfsi_algorithm.H"
#include "mfsi_couplingoperator.H"
#include "mfsi_nox_thyra_group.H"
#include "mfsi_overlappreccondfactory.H"
#include "mfsi_statustest.H"
#include "mfsi_overlappreccondfactory.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_colors.H"
#include "../drt_lib/linalg_utils.H"

#include <Teuchos_StandardParameterEntryValidators.hpp>

#include <NOX.H>

#include <Thyra_EpetraThyraWrappers.hpp>
#include <Thyra_EpetraLinearOp.hpp>
#include <Thyra_get_Epetra_Operator.hpp>

#include <Thyra_VectorStdOps.hpp>
#include <Thyra_DefaultIdentityLinearOp.hpp>
#include <Thyra_DefaultAddedLinearOp.hpp>
#include <Thyra_DefaultScaledAdjointLinearOp.hpp>

// fix clashes between ccarat and Thyra::AmesosLinearOpWithSolveFactory
#ifdef UMFPACK
#undef UMFPACK
#endif

#include <Thyra_DefaultRealLinearSolverBuilder.hpp>
#include <Thyra_AmesosLinearOpWithSolveFactory.hpp>
#include <Thyra_AztecOOLinearOpWithSolveFactory.hpp>

#ifdef PARALLEL
#include <mpi.h>
#endif

extern "C"
{
#include "../headers/standardtypes.h"
}


/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;

/*----------------------------------------------------------------------*
 | global variable *solv, vector of lenght numfld of structures SOLVAR  |
 | defined in solver_control.c                                          |
 |                                                                      |
 |                                                       m.gee 11/00    |
 *----------------------------------------------------------------------*/
extern struct _SOLVAR  *solv;

/*!----------------------------------------------------------------------
\brief file pointers

<pre>                                                         m.gee 8/00
This structure struct _FILES allfiles is defined in input_control_global.c
and the type is in standardtypes.h
It holds all file pointers and some variables needed for the FRSYSTEM
</pre>
*----------------------------------------------------------------------*/
extern struct _FILES  allfiles;

/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | structure of flags to control output                                 |
 | defined in out_global.c                                              |
 *----------------------------------------------------------------------*/
extern struct _IO_FLAGS     ioflags;


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
extern Teuchos::RCP<Teuchos::ParameterList> globalparameterlist;


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MFSI::OverlapAlgorithm::OverlapAlgorithm(Epetra_Comm& comm)
  : Algorithm(comm)
{
  // right now we use matching meshes at the interface

  Coupling& coupsf = StructureFluidCoupling();
  Coupling& coupsa = StructureAleCoupling();
  Coupling& coupfa = FluidAleCoupling();

  // structure to fluid

  coupsf.SetupConditionCoupling(StructureField()->Interface(),
                                FluidField()->Interface());

  // structure to ale

  coupsa.SetupConditionCoupling(StructureField()->Interface(),
                                AleField()->Interface());

  // In the following we assume that both couplings find the same dof
  // map at the structural side. This enables us to use just one
  // interface dof map for all fields and have just one transfer
  // operator from the interface map to the full field map.
  if (not coupsf.MasterDofMap()->SameAs(*coupsa.MasterDofMap()))
    dserror("structure interface dof maps do not match");

  if (coupsf.MasterDofMap()->NumGlobalElements()==0)
    dserror("No nodes in matching FSI interface. Empty FSI coupling condition?");

  // init transfer from interface to field
  StructureField()->SetInterfaceMap(coupsf.MasterDofMap());
  FluidField()    ->SetInterfaceMap(coupsf.SlaveDofMap());
  AleField()      ->SetInterfaceMap(coupsa.SlaveDofMap());

  // the fluid-ale coupling always matches
  const Epetra_Map* fluidnodemap = FluidField()->Discretization()->NodeRowMap();
  const Epetra_Map* alenodemap   = AleField()->Discretization()->NodeRowMap();

  coupfa.SetupCoupling(*FluidField()->Discretization(),
                       *AleField()->Discretization(),
                       *fluidnodemap,
                       *alenodemap);

  FluidField()->SetMeshMap(coupfa.MasterDofMap());

  // create combined map

  int numBlocks = 3;
  std::vector<Teuchos::RCP<const Thyra::VectorSpaceBase<double> > > vecSpaces(numBlocks);

  //simap_ = Thyra::create_VectorSpace(StructureField()->Interface().OtherDofMap());
  //sgmap_ = Thyra::create_VectorSpace(StructureField()->Interface().CondDofMap());
  smap_ =  Thyra::create_VectorSpace(StructureField()->DofRowMap());
  fimap_ = Thyra::create_VectorSpace(FluidField()    ->Interface().OtherDofMap());
  aimap_ = Thyra::create_VectorSpace(AleField()      ->Interface().OtherDofMap());

  //vecSpaces[0] = simap_;
  //vecSpaces[1] = sgmap_;
  vecSpaces[0] = smap_;
  vecSpaces[1] = fimap_;
  vecSpaces[2] = aimap_;

  SetDofRowMap(Teuchos::rcp(new Thyra::DefaultProductVectorSpace<double>(numBlocks, &vecSpaces[0])));

  // field solvers used within the block preconditioner
  structsolverfactory_ = Teuchos::rcp(new Thyra::AmesosLinearOpWithSolveFactory(Thyra::Amesos::UMFPACK));
  //interfacesolverfactory_ = Teuchos::rcp(new Thyra::AztecOOLinearOpWithSolveFactory());
  fluidsolverfactory_ = Teuchos::rcp(new Thyra::AmesosLinearOpWithSolveFactory(Thyra::Amesos::UMFPACK));
  alesolverfactory_ = Teuchos::rcp(new Thyra::AmesosLinearOpWithSolveFactory(Thyra::Amesos::UMFPACK));

  // the factory to create the special block preconditioner
  SetPreconditionerFactory(
    Teuchos::rcp(new MFSI::OverlappingPCFactory(structsolverfactory_,
                                                //interfacesolverfactory_,
                                                fluidsolverfactory_,
                                                alesolverfactory_)));

  // Lets use aztec for now. This a about the only choice we have got.
  SetSolverFactory(Teuchos::rcp(new Thyra::AztecOOLinearOpWithSolveFactory()));
  SolverFactory()->setPreconditionerFactory(PreconditionerFactory(), "FSI block preconditioner");
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MFSI::OverlapAlgorithm::InitialGuess(Thyra::DefaultProductVector<double>& ig)
{
  SetupVector(ig,
              StructureField()->InitialGuess(),
              FluidField()->InitialGuess(),
              AleField()->InitialGuess(),
              0.0);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MFSI::OverlapAlgorithm::SetupRHS(Thyra::DefaultProductVector<double> &f) const
{
  SetupVector(f,
              StructureField()->RHS(),
              FluidField()->RHS(),
              AleField()->RHS(),
              FluidField()->ResidualScaling());

  // NOX expects a different sign here.
  Thyra::scale(-1., &f);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MFSI::OverlapAlgorithm::SetupVector(Thyra::DefaultProductVector<double> &f,
                                         Teuchos::RCP<const Epetra_Vector> sv,
                                         Teuchos::RCP<const Epetra_Vector> fv,
                                         Teuchos::RCP<const Epetra_Vector> av,
                                         double fluidscale) const
{

  // extract the inner and boundary dofs of all three fields

  //Teuchos::RCP<Epetra_Vector> sov = StructureField()->Interface().ExtractOtherVector(sv);
  Teuchos::RCP<Epetra_Vector> fov = FluidField()    ->Interface().ExtractOtherVector(fv);
  Teuchos::RCP<Epetra_Vector> aov = AleField()      ->Interface().ExtractOtherVector(av);

  //Teuchos::RCP<Epetra_Vector> scv = StructureField()->Interface().ExtractCondVector(sv);
  Teuchos::RCP<Epetra_Vector> fcv = FluidField()    ->Interface().ExtractCondVector(fv);

  //scv->Update(fluidscale, *FluidToStruct(fcv), 1.0);

  Teuchos::RCP<Epetra_Vector> modsv = StructureField()->Interface().InsertCondVector(FluidToStruct(fcv));
  modsv->Update(1.0, *sv, fluidscale);

  int numBlocks = 3;
  std::vector<Teuchos::RCP<Thyra::VectorBase<double> > > vec(numBlocks);

  //vec[0] = Thyra::create_Vector(sov,simap_);
  //vec[1] = Thyra::create_Vector(scv,sgmap_);
  vec[0] = Thyra::create_Vector(modsv,smap_);
  vec[1] = Thyra::create_Vector(fov,fimap_);
  vec[2] = Thyra::create_Vector(aov,aimap_);

  f.initialize(DofRowMap(),&vec[0]);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MFSI::OverlapAlgorithm::ExtractFieldVectors(Teuchos::RCP<const Thyra::DefaultProductVector<double> > x,
                                                 Teuchos::RCP<const Epetra_Vector>& sx,
                                                 Teuchos::RCP<const Epetra_Vector>& fx,
                                                 Teuchos::RCP<const Epetra_Vector>& ax) const
{
  sx = Thyra::get_Epetra_Vector(*StructureField()->DofRowMap(), x->getVectorBlock(0));

  Teuchos::RCP<const Epetra_Vector> fox = Thyra::get_Epetra_Vector(*FluidField()    ->Interface().OtherDofMap(), x->getVectorBlock(1));
  Teuchos::RCP<const Epetra_Vector> aox = Thyra::get_Epetra_Vector(*AleField()      ->Interface().OtherDofMap(), x->getVectorBlock(2));

  Teuchos::RCP<const Epetra_Vector> scx = StructureField()->Interface().ExtractCondVector(sx);

  Teuchos::RCP<Epetra_Vector> fcx = StructToFluid(scx);
  Teuchos::RCP<Epetra_Vector> acx = StructToAle(scx);

  // We convert Delta d to Delta u here.
  fcx->Scale(1./Dt());

  Teuchos::RCP<Epetra_Vector> f = FluidField()    ->Interface().InsertOtherVector(fox);
  Teuchos::RCP<Epetra_Vector> a = AleField()      ->Interface().InsertOtherVector(aox);

  FluidField()    ->Interface().InsertCondVector(fcx, f);
  AleField()      ->Interface().InsertCondVector(acx, a);

  fx = f;
  ax = a;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MFSI::OverlapAlgorithm::SetupSysMat(Thyra::DefaultBlockedLinearOp<double>& mat) const
{
  // extract Jacobian matrices and put them into composite system
  // matrix W

  const Coupling& coupsf = StructureFluidCoupling();
  const Coupling& coupsa = StructureAleCoupling();

  Teuchos::RCP<Epetra_CrsMatrix> s = Teuchos::rcp_dynamic_cast<Epetra_CrsMatrix>(StructureField()->SysMat());

  // split fluid matrix

  Teuchos::RCP<Epetra_CrsMatrix> fii;
  Teuchos::RCP<Epetra_CrsMatrix> fig;
  Teuchos::RCP<Epetra_CrsMatrix> fgi;
  Teuchos::RCP<Epetra_CrsMatrix> fgg;

  Teuchos::RCP<Epetra_Map> fgmap = FluidField()->Interface().CondDofMap();
  Teuchos::RCP<Epetra_Map> fimap = FluidField()->Interface().OtherDofMap();

  Teuchos::RCP<Epetra_CrsMatrix> f = Teuchos::rcp_dynamic_cast<Epetra_CrsMatrix>(FluidField()->SysMat());

  // map between fluid interface and structure column map

  if (fluidstructcolmap_.size()==0)
  {
    DRT::Exporter ex(f->RowMap(),f->ColMap(),f->Comm());
    coupsf.FillSlaveToMasterMap(fluidstructcolmap_);
    ex.Export(fluidstructcolmap_);
  }

  if (not LINALG::SplitMatrix2x2(f,fimap,fgmap,fii,fig,fgi,fgg))
  {
    dserror("failed to split fluid matrix");
  }

  double scale = FluidField()->ResidualScaling();

  // transform fluid interface matrix to structure interface and add it to
  // structure matrix

  AddFluidInterface(scale,fgg,s);

  // split ale matrix

  Teuchos::RCP<Epetra_CrsMatrix> aii;
  Teuchos::RCP<Epetra_CrsMatrix> aig;
  Teuchos::RCP<Epetra_CrsMatrix> agi;
  Teuchos::RCP<Epetra_CrsMatrix> agg;

  Teuchos::RCP<Epetra_Map> agmap = AleField()->Interface().CondDofMap();
  Teuchos::RCP<Epetra_Map> aimap = AleField()->Interface().OtherDofMap();

  Teuchos::RCP<Epetra_CrsMatrix> a = Teuchos::rcp_dynamic_cast<Epetra_CrsMatrix>(AleField()->SysMat());

  // map between ale interface and structure column map

  if (alestructcolmap_.size()==0)
  {
    DRT::Exporter ex(a->RowMap(),a->ColMap(),a->Comm());
    coupsa.FillSlaveToMasterMap(alestructcolmap_);
    ex.Export(alestructcolmap_);
  }

  if (not LINALG::SplitMatrix2x2(a,aimap,agmap,aii,aig,agi,agg))
  {
    dserror("failed to split ale matrix");
  }

  // build block matrix

  mat.beginBlockFill(DofRowMap(), DofRowMap());

  mat.setBlock(0,0,Thyra::nonconstEpetraLinearOp(s));
  mat.setBlock(0,1,Thyra::nonconstEpetraLinearOp(ConvertFgiRowmap(fgi,scale,s->RowMap())));

//   mat.setBlock(1,2,
//                Thyra::nonconstScale<double>(
//                  scale,
//                  nonconstCouplingOp(
//                    fgi,
//                    Teuchos::null,
//                    Teuchos::rcp(&coupsf,false))
//                  ));
//   mat.setBlock(2,1,nonconstCouplingOp(fig,Teuchos::rcp(&coupsf,false)));
  mat.setBlock(1,0,Thyra::nonconstEpetraLinearOp(ConvertFigColmap(fig,
                                                                  fluidstructcolmap_,
                                                                  s->DomainMap())));

  mat.setBlock(1,1,Thyra::nonconstEpetraLinearOp(fii));

  //mat.setBlock(3,1,nonconstCouplingOp(aig,Teuchos::rcp(&coupsa,false)));
  mat.setBlock(2,0,Thyra::nonconstEpetraLinearOp(ConvertFigColmap(aig,
                                                                  alestructcolmap_,
                                                                  s->DomainMap())));
  mat.setBlock(2,2,Thyra::nonconstEpetraLinearOp(aii));

  mat.endBlockFill();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MFSI::OverlapAlgorithm::AddFluidInterface(double scale,
                                               Teuchos::RCP<Epetra_CrsMatrix> fgg,
                                               Teuchos::RCP<Epetra_CrsMatrix> s) const
{
  const Coupling& coupsf = StructureFluidCoupling();

  Teuchos::RCP<Epetra_CrsMatrix> perm_fgg = coupsf.SlaveToPermSlave(fgg);

  const Epetra_Map& perm_rowmap = perm_fgg->RowMap();
  const Epetra_Map& perm_colmap = perm_fgg->ColMap();

  if (permfluidstructcolmap_.size()==0)
  {
    // fgg->RowMap()
    DRT::Exporter ex(*coupsf.SlaveDofMap(),perm_colmap,perm_fgg->Comm());
    coupsf.FillSlaveToMasterMap(permfluidstructcolmap_);
    ex.Export(permfluidstructcolmap_);
  }

  int rows = perm_fgg->NumMyRows();
  for (int i=0; i<rows; ++i)
  {
    int NumEntries;
    double *Values;
    int *Indices;
    int err = perm_fgg->ExtractMyRowView(i, NumEntries, Values, Indices);
    if (err!=0)
      dserror("ExtractMyRowView error: %d", err);

    int gid = perm_rowmap.GID(i);
    std::map<int,int>::iterator iter = permfluidstructcolmap_.find(gid);
    if (iter==permfluidstructcolmap_.end())
      dserror("gid %d not found", gid);
    int structgid = iter->second;

    for (int j=0; j<NumEntries; ++j)
    {
      int perm_gid = perm_colmap.GID(Indices[j]);
      iter = permfluidstructcolmap_.find(perm_gid);
      if (iter==permfluidstructcolmap_.end())
        dserror("gid %d not found", perm_gid);
      int index = iter->second;
      double value = Values[j]*scale;
      err = s->SumIntoGlobalValues(structgid, 1, &value, &index);
      if (err)
        dserror("SumIntoGlobalValues error: %d", err);
    }
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_CrsMatrix>
MFSI::OverlapAlgorithm::ConvertFigColmap(Teuchos::RCP<Epetra_CrsMatrix> fig,
                                         const std::map<int,int>& convcolmap,
                                         const Epetra_Map& domainmap) const
{
  const Epetra_Map& rowmap = fig->RowMap();
  const Epetra_Map& colmap = fig->ColMap();

  Teuchos::RCP<Epetra_CrsMatrix> convfig = Teuchos::rcp(new Epetra_CrsMatrix(Copy,rowmap,
                                                                             fig->MaxNumEntries()));

  int rows = fig->NumMyRows();
  for (int i=0; i<rows; ++i)
  {
    int NumEntries;
    double *Values;
    int *Indices;
    int err = fig->ExtractMyRowView(i, NumEntries, Values, Indices);
    if (err!=0)
      dserror("ExtractMyRowView error: %d", err);

    for (int j=0; j<NumEntries; ++j)
    {
      int gid = colmap.GID(Indices[j]);
      std::map<int,int>::const_iterator iter = convcolmap.find(gid);
      if (iter==convcolmap.end())
        dserror("gid %d not found in map", gid);
      Indices[j] = iter->second;
    }

    err = convfig->InsertGlobalValues(rowmap.GID(i), NumEntries, Values, Indices);
    if (err)
      dserror("InsertGlobalValues error: %d", err);
  }

  convfig->FillComplete(domainmap,rowmap);

  return convfig;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_CrsMatrix>
MFSI::OverlapAlgorithm::ConvertFgiRowmap(Teuchos::RCP<Epetra_CrsMatrix> fgi,
                                         double scale,
                                         const Epetra_Map& structrowmap) const
{
  const Coupling& coupsf = StructureFluidCoupling();

  // redistribute fluid matrix to match distribution of structure matrix

  fgi = coupsf.SlaveToPermSlave(fgi);

  // create new matrix with structure row map and fill it

  Teuchos::RCP<Epetra_CrsMatrix> sfgi = Teuchos::rcp(new Epetra_CrsMatrix(Copy,structrowmap,fgi->MaxNumEntries()));

  //const Epetra_Map& rowmap = fgi->RowMap();
  const Epetra_Map& colmap = fgi->ColMap();

  int rows = fgi->NumMyRows();
  for (int i=0; i<rows; ++i)
  {
    int NumEntries;
    double *Values;
    int *Indices;
    int err = fgi->ExtractMyRowView(i, NumEntries, Values, Indices);
    if (err!=0)
      dserror("ExtractMyRowView error: %d", err);

    // pull indices back to global
    for (int j=0; j<NumEntries; ++j)
    {
      Indices[j] = colmap.GID(Indices[j]);
    }

    // put row into matrix with structure row map
    err = sfgi->InsertGlobalValues(coupsf.MasterDofMap()->GID(i), NumEntries, Values, Indices);
    if (err)
      dserror("InsertGlobalValues error: %d", err);
  }

  sfgi->FillComplete(fgi->DomainMap(),structrowmap);
  sfgi->Scale(scale);
  return sfgi;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<NOX::StatusTest::Combo>
MFSI::OverlapAlgorithm::CreateStatusTest(Teuchos::ParameterList& nlParams,
                                         Teuchos::RCP<NOX::Thyra::Group> grp)
{
  // Create the convergence tests
  Teuchos::RCP<NOX::StatusTest::Combo> combo       = Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR));
  Teuchos::RCP<NOX::StatusTest::Combo> converged   = Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::AND));

  Teuchos::RCP<NOX::StatusTest::MaxIters> maxiters = Teuchos::rcp(new NOX::StatusTest::MaxIters(nlParams.get("Max Iterations", 100)));
  Teuchos::RCP<NOX::StatusTest::FiniteValue> fv    = Teuchos::rcp(new NOX::StatusTest::FiniteValue);

  combo->addStatusTest(fv);
  combo->addStatusTest(converged);
  combo->addStatusTest(maxiters);

  // Very simple absolute test. We need to do something more
  // sophisticated here.

  Teuchos::RCP<NOX::StatusTest::NormF> absresid =
    Teuchos::rcp(new NOX::StatusTest::NormF(nlParams.get("Norm abs F", 1.0e-6)));
  converged->addStatusTest(absresid);

  return combo;
}


#endif
