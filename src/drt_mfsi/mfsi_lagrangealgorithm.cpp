/*----------------------------------------------------------------------*/
/*!
\file mfsi_algorithm.cpp

\brief Solve monolithic Lagrangian block coupled FSI problems

<pre>
Maintainer: Ulrich Kuettler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15238
</pre>
*/
/*----------------------------------------------------------------------*/

#ifdef CCADISCRET

#include "mfsi_lagrangealgorithm.H"

#include "mfsi_algorithm.H"
#include "mfsi_nox_thyra_group.H"
#include "mfsi_preconditionerfactory.H"
#include "mfsi_statustest.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_colors.H"

#include <Teuchos_StandardParameterEntryValidators.hpp>

#include <NOX.H>

#include <Thyra_EpetraThyraWrappers.hpp>
#include <Thyra_EpetraLinearOp.hpp>
#include <Thyra_get_Epetra_Operator.hpp>

#include <Thyra_VectorStdOps.hpp>
#include <Thyra_DefaultIdentityLinearOp.hpp>

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
MFSI::LagrangeAlgorithm::LagrangeAlgorithm(Epetra_Comm& comm)
  : Algorithm(comm)
{
  // right now we use matching meshes at the interface

  Coupling& coupsf = StructureFluidCoupling();
  Coupling& coupsa = StructureAleCoupling();
  Coupling& coupfa = FluidAleCoupling();

  // structure to fluid

  coupsf.SetupConditionCoupling(*StructureField()->Discretization(),
                                StructureField()->Interface(),
                                *FluidField()->Discretization(),
                                FluidField()->Interface(),
                                "FSICoupling");

  // setup the very simple structure to fluid coupling
  // u(n+1)*dt = d(n+1) - d(n)
  // but we solve both fields for the increments
  // the structure is even solved for the middle increment
  //
  // du(n+1)*dt = 1/(1-alpha_f)*dd(n+m)

  coupsf.SetupCouplingMatrices(*StructureField()->Discretization()->DofRowMap(),
                               *FluidField()->Discretization()->DofRowMap());

//   coupsf.MasterToMasterMat()->Scale(StructureField()->DispIncrFactor());
//   coupsf.MasterToMasterMatTrans()->Scale(StructureField()->DispIncrFactor());
  coupsf.SlaveToMasterMat()->Scale(-Dt());
  coupsf.SlaveToMasterMatTrans()->Scale(-Dt());

  // structure to ale

  coupsa.SetupConditionCoupling(*StructureField()->Discretization(),
                                StructureField()->Interface(),
                                *AleField()->Discretization(),
                                AleField()->Interface(),
                                "FSICoupling");

  // setup structure to ale coupling
  //
  // dd(G,n+1) = 1/(1-alpha_f)*dd(n+m)

  coupsa.SetupCouplingMatrices(*StructureField()->Discretization()->DofRowMap(),
                                *AleField()->Discretization()->DofRowMap());

//   coupsa.MasterToMasterMat()->Scale(StructureField()->DispIncrFactor());
//   coupsa.MasterToMasterMatTrans()->Scale(StructureField()->DispIncrFactor());
  coupsa.SlaveToMasterMat()->Scale(-1.);
  coupsa.SlaveToMasterMatTrans()->Scale(-1.);

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

  ifstruct_ = Teuchos::rcp(new Epetra_Vector(*StructureField()->InterfaceMap()));
  iastruct_ = Teuchos::rcp(new Epetra_Vector(*StructureField()->InterfaceMap()));

  // the fluid-ale coupling always matches
  const Epetra_Map* fluidnodemap = FluidField()->Discretization()->NodeRowMap();
  const Epetra_Map* alenodemap   = AleField()->Discretization()->NodeRowMap();

  coupfa.SetupCoupling(*FluidField()->Discretization(),
                        *AleField()->Discretization(),
                        *fluidnodemap,
                        *alenodemap);

  FluidField()->SetMeshMap(coupfa.MasterDofMap());

  // create Thyra vector spaces from Epetra maps for all fields
  smap_ = Thyra::create_VectorSpace(StructureField()->DofRowMap());
  fmap_ = Thyra::create_VectorSpace(FluidField()->DofRowMap());
  amap_ = Thyra::create_VectorSpace(AleField()->DofRowMap());

  sfmap_ = Thyra::create_VectorSpace(coupsf.MasterDofMap());

  // stack vector spaces to build a composite that contains them all
  int numBlocks = 5;
  std::vector<Teuchos::RCP<const Thyra::VectorSpaceBase<double> > > vecSpaces(numBlocks);
  vecSpaces[0] = smap_;
  vecSpaces[1] = fmap_;
  vecSpaces[2] = amap_;
  vecSpaces[3] = sfmap_;
  vecSpaces[4] = sfmap_;
  SetDofRowMap(Teuchos::rcp(new Thyra::DefaultProductVectorSpace<double>(numBlocks, &vecSpaces[0])));

#if 0
  Teuchos::RCP<Teuchos::FancyOStream> out = Teuchos::VerboseObjectBase::getDefaultOStream();

  // use a solver builder for the standard field solvers so we get all the default behaviour
  Thyra::DefaultRealLinearSolverBuilder linearSolverBuilder;

  if (Comm().MyPID()==0)
  {
    std::cout << *linearSolverBuilder.getValidParameters() << std::endl;

    //linearSolverBuilder.getValidParameters()->print(std::cout,
    //                                                Teuchos::ParameterList::PrintOptions().showDoc(true).indent(2).showTypes(true));
  }

  Teuchos::RCP< Teuchos::ParameterList > paramList = Teuchos::rcp(new ParameterList());
  linearSolverBuilder.setParameterList(paramList);

  if (Comm().MyPID()==0)
  {
    paramList->print(std::cout);
  }

  structsolverfactory_ = linearSolverBuilder.createLinearSolveStrategy("");
  structsolverfactory_->setOStream(out);
  structsolverfactory_->setVerbLevel(Teuchos::VERB_LOW);

  fluidsolverfactory_ = linearSolverBuilder.createLinearSolveStrategy("");
  fluidsolverfactory_->setOStream(out);
  fluidsolverfactory_->setVerbLevel(Teuchos::VERB_LOW);

  alesolverfactory_ = linearSolverBuilder.createLinearSolveStrategy("");
  alesolverfactory_->setOStream(out);
  alesolverfactory_->setVerbLevel(Teuchos::VERB_LOW);
#else
  // field solvers used within the block preconditioner
  structsolverfactory_ = Teuchos::rcp(new Thyra::AmesosLinearOpWithSolveFactory(Thyra::Amesos::KLU));
  fluidsolverfactory_ = Teuchos::rcp(new Thyra::AmesosLinearOpWithSolveFactory(Thyra::Amesos::UMFPACK));
  alesolverfactory_ = Teuchos::rcp(new Thyra::AmesosLinearOpWithSolveFactory(Thyra::Amesos::KLU));
#endif

  sfidentity_ = Teuchos::rcp(new Thyra::DefaultIdentityLinearOp<double>(Thyra::create_VectorSpace(StructureField()->InterfaceMap())));

  //Thyra::ConstLinearOperator<double> sfihandle = sfidentity;

  // the factory to create the special block preconditioner
  SetPreconditionerFactory(
    Teuchos::rcp(new MFSI::PreconditionerFactory(structsolverfactory_,
                                                 fluidsolverfactory_,
                                                 alesolverfactory_,
                                                 sfidentity_)));

  // Lets use aztec for now. This a about the only choice we have got.
  SetSolverFactory(Teuchos::rcp(new Thyra::AztecOOLinearOpWithSolveFactory()));
  SolverFactory()->setPreconditionerFactory(PreconditionerFactory(), "FSI block preconditioner");

  // We cannot use Amesos, since it expects the unterlying matrix to
  // be a Epetra_Operator.
  //solverfactory_ = Teuchos::rcp(new Thyra::AmesosLinearOpWithSolveFactory());
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MFSI::LagrangeAlgorithm::InitialGuess(Thyra::DefaultProductVector<double>& ig)
{
  // the linear field systems must be setup before the initial guess
  // is known
  Teuchos::RCP< const Thyra::VectorBase< double > > sig = Thyra::create_Vector(StructureField()->InitialGuess(), smap_);
  Teuchos::RCP< const Thyra::VectorBase< double > > fig = Thyra::create_Vector(FluidField()->InitialGuess(), fmap_);
  Teuchos::RCP< const Thyra::VectorBase< double > > aig = Thyra::create_Vector(AleField()->InitialGuess(), amap_);

  Teuchos::RCP< const Thyra::VectorBase< double > > sfig =
    Thyra::create_Vector(ifstruct_, Thyra::create_VectorSpace(StructureField()->InterfaceMap()));

  Teuchos::RCP< const Thyra::VectorBase< double > > saig =
    Thyra::create_Vector(iastruct_, Thyra::create_VectorSpace(StructureField()->InterfaceMap()));

  int numBlocks = 5;
  std::vector<Teuchos::RCP<const Thyra::VectorBase<double> > > vec(numBlocks);
  vec[0] = sig;
  vec[1] = fig;
  vec[2] = aig;
  vec[3] = sfig;
  vec[4] = saig;
  ig.initialize(DofRowMap(),&vec[0]);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MFSI::LagrangeAlgorithm::SetupRHS(Thyra::DefaultProductVector<double> &f) const
{
  // Extract RHS and put it into f

  // make a local copy so that we can modify the rhs vectors
  Teuchos::RCP<Epetra_Vector> sv = rcp(new Epetra_Vector(*StructureField()->RHS()));
  Teuchos::RCP<Epetra_Vector> fv = rcp(new Epetra_Vector(*FluidField()->RHS()));
  Teuchos::RCP<Epetra_Vector> av = rcp(new Epetra_Vector(*AleField()->RHS()));

//   debug_.DumpVector("sf",*StructureField()->Discretization(),*sv);
//   debug_.DumpVector("ff",*FluidField()->Discretization(),*fv);
//   debug_.DumpVector("af",*AleField()->Discretization(),*av);

  // wrap epetra vectors in thyra
  Teuchos::RCP< Thyra::VectorBase< double > > srhs = Thyra::create_Vector(sv, smap_);
  Teuchos::RCP< Thyra::VectorBase< double > > frhs = Thyra::create_Vector(fv, fmap_);
  Teuchos::RCP< Thyra::VectorBase< double > > arhs = Thyra::create_Vector(av, amap_);

  // We couple absolute vectors, no increments. So we have a nonzero rhs.

  Teuchos::RCP< Epetra_Vector > ifstruct = StructureField()->FluidCondRHS();
  ifstruct->Update(1.0*Dt(),*FluidToStruct(FluidField()->StructCondRHS()),-1.0);
  Teuchos::RCP< Thyra::VectorBase< double > > sfrhs =
    Thyra::create_Vector(ifstruct, Thyra::create_VectorSpace(StructureField()->InterfaceMap()));

  Teuchos::RCP< Epetra_Vector > iastruct = StructureField()->MeshCondRHS();
  iastruct->Update(1.0,*AleToStruct(AleField()->StructCondRHS()),-1.0);
  Teuchos::RCP< Thyra::VectorBase< double > > sarhs =
    Thyra::create_Vector(iastruct, Thyra::create_VectorSpace(StructureField()->InterfaceMap()));

  // create block vector
  int numBlocks = 5;
  std::vector<Teuchos::RCP<Thyra::VectorBase<double> > > vec(numBlocks);
  vec[0] = srhs;
  vec[1] = frhs;
  vec[2] = arhs;
  vec[3] = sfrhs;
  vec[4] = sarhs;
  f.initialize(DofRowMap(),&vec[0]);

  // NOX expects a different sign here.
  Thyra::scale(-1., &f);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MFSI::LagrangeAlgorithm::SetupSysMat(Thyra::DefaultBlockedLinearOp<double>& mat) const
{
  // extract Jacobian matrices and put them into composite system
  // matrix W
  Teuchos::RCP<Thyra::LinearOpBase<double> > smat = Teuchos::rcp(new Thyra::EpetraLinearOp(StructureField()->SysMat()));
  Teuchos::RCP<Thyra::LinearOpBase<double> > fmat = Teuchos::rcp(new Thyra::EpetraLinearOp(FluidField()->SysMat()));
  Teuchos::RCP<Thyra::LinearOpBase<double> > amat = Teuchos::rcp(new Thyra::EpetraLinearOp(AleField()->SysMat()));

  mat.beginBlockFill(DofRowMap(), DofRowMap());
  mat.setBlock(0,0,smat);
  mat.setBlock(1,1,fmat);
  mat.setBlock(2,2,amat);

  const Coupling& coupsf = StructureFluidCoupling();
  const Coupling& coupsa = StructureAleCoupling();

  // structure to fluid coupling
  mat.setBlock(3,0,Teuchos::rcp(new Thyra::EpetraLinearOp(coupsf.MasterToMasterMat())));
  mat.setBlock(3,1,Teuchos::rcp(new Thyra::EpetraLinearOp(coupsf.SlaveToMasterMat())));
  mat.setBlock(0,3,Teuchos::rcp(new Thyra::EpetraLinearOp(coupsf.MasterToMasterMatTrans())));
  mat.setBlock(1,3,Teuchos::rcp(new Thyra::EpetraLinearOp(coupsf.SlaveToMasterMatTrans())));
  //mat.setBlock(3,3,sfidentity_);

  // structure to ale coupling
  // note there is no ale effect on the structural equations
  mat.setBlock(4,0,Teuchos::rcp(new Thyra::EpetraLinearOp(coupsa.MasterToMasterMat())));
  mat.setBlock(4,2,Teuchos::rcp(new Thyra::EpetraLinearOp(coupsa.SlaveToMasterMat())));
  //mat.setBlock(0,4,Teuchos::rcp(new Thyra::EpetraLinearOp(coupsa.MasterToMasterMatTrans())));
  mat.setBlock(2,4,Teuchos::rcp(new Thyra::EpetraLinearOp(coupsa.SlaveToMasterMatTrans())));
  //mat.setBlock(4,4,sfidentity_);

  mat.endBlockFill();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MFSI::LagrangeAlgorithm::ExtractFieldVectors(Teuchos::RCP<const Thyra::DefaultProductVector<double> > x,
                                                  Teuchos::RCP<const Epetra_Vector>& sx,
                                                  Teuchos::RCP<const Epetra_Vector>& fx,
                                                  Teuchos::RCP<const Epetra_Vector>& ax) const
{
  sx = Thyra::get_Epetra_Vector(*StructureField()->DofRowMap(), x->getVectorBlock(0));
  fx = Thyra::get_Epetra_Vector(*FluidField()    ->DofRowMap(), x->getVectorBlock(1));
  ax = Thyra::get_Epetra_Vector(*AleField()      ->DofRowMap(), x->getVectorBlock(2));
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<NOX::StatusTest::Combo>
MFSI::LagrangeAlgorithm::CreateStatusTest(Teuchos::ParameterList& nlParams,
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

  Teuchos::RCP<PartialNormF> structureDisp =
    Teuchos::rcp(new PartialNormF("displacement",
                                  0,
                                  *StructureField()->DofRowMap(),
                                  *StructureField()->InnerDisplacementRowMap(),
                                  nlParams.get("Norm abs disp", 1.0e-6),
                                  PartialNormF::Scaled));
  converged->addStatusTest(structureDisp);

  Teuchos::RCP<PartialNormF> innerFluidVel =
    Teuchos::rcp(new PartialNormF("velocity",
                                  1,
                                  *FluidField()->DofRowMap(),
                                  *FluidField()->InnerVelocityRowMap(),
                                  nlParams.get("Norm abs vel", 1.0e-6),
                                  PartialNormF::Scaled));
  converged->addStatusTest(innerFluidVel);

  Teuchos::RCP<PartialNormF> fluidPress =
    Teuchos::rcp(new PartialNormF("pressure",
                                  1,
                                  *FluidField()->DofRowMap(),
                                  *FluidField()->PressureRowMap(),
                                  nlParams.get("Norm abs pres", 1.0e-6),
                                  PartialNormF::Scaled));
  converged->addStatusTest(fluidPress);

  const Coupling& coupsf = StructureFluidCoupling();

  Teuchos::RCP<InterfaceNormF> interface =
    Teuchos::rcp(new InterfaceNormF(1.,
                                    *StructureField()->DofRowMap(),
                                    *coupsf.MasterDofMap(),
                                    FluidField()->ResidualScaling(),
                                    *FluidField()->DofRowMap(),
                                    *coupsf.SlaveDofMap(),
                                    coupsf,
                                    nlParams.get("Norm abs interface", 1.0e-6),
                                    PartialNormF::Scaled));
  converged->addStatusTest(interface);

#if 0
  if (nlParams.isParameter("Norm Update"))
  {
    Teuchos::RCP<NOX::StatusTest::NormUpdate> update =
      Teuchos::rcp(new NOX::StatusTest::NormUpdate(nlParams.get("Norm Update", 1.0e-5)));
    converged->addStatusTest(update);
  }

  if (nlParams.isParameter("Norm rel F"))
  {
    Teuchos::RCP<NOX::StatusTest::NormF> relresid =
      Teuchos::rcp(new NOX::StatusTest::NormF(*grp.get(), nlParams.get("Norm rel F", 1.0e-2)));
    converged->addStatusTest(relresid);
  }
#endif

  //Teuchos::RCP<NOX::StatusTest::NormWRMS> wrms     = Teuchos::rcp(new NOX::StatusTest::NormWRMS(1.0e-2, 1.0e-8));
  //converged->addStatusTest(wrms);

  return combo;
}


#endif
