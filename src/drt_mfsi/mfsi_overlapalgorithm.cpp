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
#include "mfsi_preconditionerfactory.H"
#include "mfsi_statustest.H"

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
  //
  // We could use the fluid interface map for the coupling
  // equations. We assume that most of the time the structural matrix
  // will be smaller that the fluid one, so there would less work to
  // do to transfer the stuff.
  //
  // On the other hand, we have no direct interface only Fluid-Ale
  // coupling. So for a start it will be simpler to use the structural
  // interface map.

  int numBlocks = 4;
  std::vector<Teuchos::RCP<const Thyra::VectorSpaceBase<double> > > vecSpaces(numBlocks);

  simap_ = Thyra::create_VectorSpace(StructureField()->Interface().OtherDofMap());
  sgmap_ = Thyra::create_VectorSpace(StructureField()->Interface().CondDofMap());
  fimap_ = Thyra::create_VectorSpace(FluidField()    ->Interface().OtherDofMap());
  aimap_ = Thyra::create_VectorSpace(AleField()      ->Interface().OtherDofMap());

  vecSpaces[0] = simap_;
  vecSpaces[1] = sgmap_;
  vecSpaces[2] = fimap_;
  vecSpaces[3] = aimap_;

  SetDofRowMap(Teuchos::rcp(new Thyra::DefaultProductVectorSpace<double>(numBlocks, &vecSpaces[0])));

}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MFSI::OverlapAlgorithm::InitialGuess(Thyra::DefaultProductVector<double>& ig)
{
  SetupVector(ig,
              StructureField()->InitialGuess(),
              FluidField()->InitialGuess(),
              AleField()->InitialGuess(),
              1.0);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MFSI::OverlapAlgorithm::SetupRHS(Thyra::DefaultProductVector<double> &f) const
{
  SetupVector(f,
              StructureField()->RHS(),
              FluidField()->RHS(),
              AleField()->RHS(),
              1.0);

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

  Teuchos::RCP<Epetra_Vector> sov = StructureField()->Interface().ExtractOtherVector(sv);
  Teuchos::RCP<Epetra_Vector> fov = FluidField()    ->Interface().ExtractOtherVector(fv);
  Teuchos::RCP<Epetra_Vector> aov = AleField()      ->Interface().ExtractOtherVector(av);

  Teuchos::RCP<Epetra_Vector> scv = StructureField()->Interface().ExtractCondVector(sv);
  Teuchos::RCP<Epetra_Vector> fcv = FluidField()    ->Interface().ExtractCondVector(fv);

  scv->Update(fluidscale, *FluidToStruct(fcv), 1.0);

  int numBlocks = 4;
  std::vector<Teuchos::RCP<Thyra::VectorBase<double> > > vec(numBlocks);

  vec[0] = Thyra::create_Vector(sov,simap_);
  vec[1] = Thyra::create_Vector(scv,sgmap_);
  vec[2] = Thyra::create_Vector(fov,fimap_);
  vec[3] = Thyra::create_Vector(aov,aimap_);

  f.initialize(DofRowMap(),&vec[0]);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MFSI::OverlapAlgorithm::SetupSysMat(Thyra::DefaultBlockedLinearOp<double>& mat) const
{
  // extract Jacobian matrices and put them into composite system
  // matrix W

  const Coupling& coupsf = StructureFluidCoupling();
  const Coupling& coupsa = StructureAleCoupling();

  // split structural matrix

  Teuchos::RCP<Epetra_CrsMatrix> sii;
  Teuchos::RCP<Epetra_CrsMatrix> sig;
  Teuchos::RCP<Epetra_CrsMatrix> sgi;
  Teuchos::RCP<Epetra_CrsMatrix> sgg;

  Teuchos::RCP<Epetra_Map> sgmap = StructureField()->Interface().CondDofMap();
  Teuchos::RCP<Epetra_Map> simap = StructureField()->Interface().OtherDofMap();

  Teuchos::RCP<Epetra_CrsMatrix> s = Teuchos::rcp_dynamic_cast<Epetra_CrsMatrix>(StructureField()->SysMat());

  if (not LINALG::SplitMatrix2x2(s,simap,sgmap,sii,sig,sgi,sgg))
  {
    dserror("failed to split structural matrix");
  }

  // split fluid matrix

  Teuchos::RCP<Epetra_CrsMatrix> fii;
  Teuchos::RCP<Epetra_CrsMatrix> fig;
  Teuchos::RCP<Epetra_CrsMatrix> fgi;
  Teuchos::RCP<Epetra_CrsMatrix> fgg;

  Teuchos::RCP<Epetra_Map> fgmap = FluidField()->Interface().CondDofMap();
  Teuchos::RCP<Epetra_Map> fimap = FluidField()->Interface().OtherDofMap();

  Teuchos::RCP<Epetra_CrsMatrix> f = Teuchos::rcp_dynamic_cast<Epetra_CrsMatrix>(FluidField()->SysMat());

  if (not LINALG::SplitMatrix2x2(f,fimap,fgmap,fii,fig,fgi,fgg))
  {
    dserror("failed to split fluid matrix");
  }

  // split ale matrix

  Teuchos::RCP<Epetra_CrsMatrix> aii;
  Teuchos::RCP<Epetra_CrsMatrix> aig;
  Teuchos::RCP<Epetra_CrsMatrix> agi;
  Teuchos::RCP<Epetra_CrsMatrix> agg;

  Teuchos::RCP<Epetra_Map> agmap = AleField()->Interface().CondDofMap();
  Teuchos::RCP<Epetra_Map> aimap = AleField()->Interface().OtherDofMap();

  Teuchos::RCP<Epetra_CrsMatrix> a = Teuchos::rcp_dynamic_cast<Epetra_CrsMatrix>(AleField()->SysMat());

  if (not LINALG::SplitMatrix2x2(a,aimap,agmap,aii,aig,agi,agg))
  {
    dserror("failed to split ale matrix");
  }

  // build block matrix

  mat.beginBlockFill(DofRowMap(), DofRowMap());

  mat.setBlock(0,0,Thyra::nonconstEpetraLinearOp(sii));
  mat.setBlock(0,1,Thyra::nonconstEpetraLinearOp(sig));
  mat.setBlock(1,0,Thyra::nonconstEpetraLinearOp(sgi));

  double scale = 1.0;

  mat.setBlock(1,1,
               Thyra::nonconstAdd<double>(
                 Thyra::nonconstEpetraLinearOp(sgg),
                 Thyra::nonconstScale<double>(
                   scale,
                   nonconstCouplingOp(
                     fgg,
                     Teuchos::rcp(&coupsf,false),
                     Teuchos::rcp(&coupsf,false)))
                 ));

  mat.setBlock(1,2,nonconstCouplingOp(fgi,Teuchos::null,Teuchos::rcp(&coupsf,false)));
  mat.setBlock(2,1,nonconstCouplingOp(fig,Teuchos::rcp(&coupsf,false)));
  mat.setBlock(2,2,Thyra::nonconstEpetraLinearOp(fii));

  mat.setBlock(3,1,nonconstCouplingOp(aig,Teuchos::rcp(&coupsa,false)));
  mat.setBlock(3,3,Thyra::nonconstEpetraLinearOp(aii));

  mat.endBlockFill();
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

  // ToDo

  return combo;
}


#endif
