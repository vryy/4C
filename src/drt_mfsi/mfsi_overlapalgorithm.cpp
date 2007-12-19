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

  coupsf.SetupConditionCoupling(*StructureField()->Discretization(),
                                 *FluidField()->Discretization(),
                                 "FSICoupling");

  coupsf.SetupInnerDofMaps(*StructureField()->Discretization()->DofRowMap(),
                           *FluidField()->Discretization()->DofRowMap());

  // structure to ale

  coupsa.SetupConditionCoupling(*StructureField()->Discretization(),
                                 *AleField()->Discretization(),
                                 "FSICoupling");

  coupsa.SetupInnerDofMaps(*StructureField()->Discretization()->DofRowMap(),
                           *AleField()->Discretization()->DofRowMap());

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

  vecSpaces[0] = Thyra::create_VectorSpace(coupsf.MasterInnerDofMap());
  vecSpaces[1] = Thyra::create_VectorSpace(coupsf.MasterDofMap());
  vecSpaces[2] = Thyra::create_VectorSpace(coupsf.SlaveInnerDofMap());
  vecSpaces[3] = Thyra::create_VectorSpace(coupsa.SlaveInnerDofMap());

  SetDofRowMap(Teuchos::rcp(new Thyra::DefaultProductVectorSpace<double>(numBlocks, &vecSpaces[0])));

}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MFSI::OverlapAlgorithm::InitialGuess(Thyra::DefaultProductVector<double>& ig)
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MFSI::OverlapAlgorithm::SetupRHS(Thyra::DefaultProductVector<double> &f) const
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MFSI::OverlapAlgorithm::SetupSysMat(Thyra::DefaultBlockedLinearOp<double>& mat) const
{
  // extract Jacobian matrices and put them into composite system
  // matrix W
  //Teuchos::RCP<Thyra::LinearOpBase<double> > smat = Teuchos::rcp(new Thyra::EpetraLinearOp(StructureField()->SysMat()));
  //Teuchos::RCP<Thyra::LinearOpBase<double> > fmat = Teuchos::rcp(new Thyra::EpetraLinearOp(FluidField()->SysMat()));
  //Teuchos::RCP<Thyra::LinearOpBase<double> > amat = Teuchos::rcp(new Thyra::EpetraLinearOp(AleField()->SysMat()));

  const Coupling& coupsf = StructureFluidCoupling();
  const Coupling& coupsa = StructureAleCoupling();

  // split structural matrix

  Teuchos::RCP<Epetra_CrsMatrix> sii;
  Teuchos::RCP<Epetra_CrsMatrix> sig;
  Teuchos::RCP<Epetra_CrsMatrix> sgi;
  Teuchos::RCP<Epetra_CrsMatrix> sgg;

  Teuchos::RCP<Epetra_Map> sgmap = coupsf.MasterDofMap();
  Teuchos::RCP<Epetra_Map> simap = coupsf.MasterInnerDofMap();

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

  Teuchos::RCP<Epetra_Map> fgmap = coupsf.SlaveDofMap();
  Teuchos::RCP<Epetra_Map> fimap = coupsf.SlaveInnerDofMap();

  Teuchos::RCP<Epetra_CrsMatrix> f = Teuchos::rcp_dynamic_cast<Epetra_CrsMatrix>(FluidField()->SysMat());

  if (not LINALG::SplitMatrix2x2(f,fimap,fgmap,fii,fig,fgi,fgg))
  {
    dserror("failed to split fluid matrix");
  }

//   fgg = coupsf.RCSlaveToMaster(fgg);
//   fgi = coupsf.RSlaveToMaster(fgi);
//   fig = coupsf.CSlaveToMaster(fig);

  // split ale matrix

  Teuchos::RCP<Epetra_CrsMatrix> aii;
  Teuchos::RCP<Epetra_CrsMatrix> aig;
  Teuchos::RCP<Epetra_CrsMatrix> agi;
  Teuchos::RCP<Epetra_CrsMatrix> agg;

  Teuchos::RCP<Epetra_Map> agmap = coupsa.SlaveDofMap();
  Teuchos::RCP<Epetra_Map> aimap = coupsa.SlaveInnerDofMap();

  Teuchos::RCP<Epetra_CrsMatrix> a = Teuchos::rcp_dynamic_cast<Epetra_CrsMatrix>(AleField()->SysMat());

  if (not LINALG::SplitMatrix2x2(a,aimap,agmap,aii,aig,agi,agg))
  {
    dserror("failed to split ale matrix");
  }

//   aig = coupsa.CSlaveToMaster(aig);

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
