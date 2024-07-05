/*----------------------------------------------------------------------------*/
/*! \file


\brief Solve FSI problem with volume constraints

\level 2
*/

/*----------------------------------------------------------------------------*/

#include "4C_fsi_constrmonolithic.hpp"

#include "4C_adapter_ale_fsi.hpp"
#include "4C_adapter_fld_fluid_fsi.hpp"
#include "4C_adapter_str_fsiwrapper.hpp"
#include "4C_adapter_str_structure.hpp"
#include "4C_ale_utils_mapextractor.hpp"
#include "4C_constraint_manager.hpp"
#include "4C_coupling_adapter.hpp"
#include "4C_fem_discretization.hpp"
#include "4C_fluid_utils_mapextractor.hpp"
#include "4C_fsi_overlapprec_fsiamg.hpp"
#include "4C_fsi_statustest.hpp"
#include "4C_global_data.hpp"
#include "4C_inpar_fsi.hpp"
#include "4C_io_control.hpp"
#include "4C_linalg_sparsematrix.hpp"
#include "4C_structure_aux.hpp"

#include <NOX_Epetra_LinearSystem.H>
#include <NOX_Epetra_LinearSystem_AztecOO.H>

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FSI::ConstrMonolithic::ConstrMonolithic(
    const Epetra_Comm& comm, const Teuchos::ParameterList& timeparams)
    : BlockMonolithic(comm, timeparams), conman_(structure_field()->get_constraint_manager())
{
  icoupfa_ = Teuchos::rcp(new Core::Adapter::Coupling());
  coupsaout_ = Teuchos::rcp(new Core::Adapter::Coupling());
  coupfsout_ = Teuchos::rcp(new Core::Adapter::Coupling());
  coupfaout_ = Teuchos::rcp(new Core::Adapter::Coupling());

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::ConstrMonolithic::GeneralSetup()
{
  const Teuchos::ParameterList& fsidyn = Global::Problem::Instance()->FSIDynamicParams();
  const Teuchos::ParameterList& fsimono = fsidyn.sublist("MONOLITHIC SOLVER");
  linearsolverstrategy_ =
      Core::UTILS::IntegralValue<Inpar::FSI::LinearBlockSolver>(fsimono, "LINEARBLOCKSOLVER");

  set_default_parameters(fsidyn, nox_parameter_list());

  // ToDo: Set more detailed convergence tolerances like in standard FSI
  // additionally set tolerance for volume constraint
  nox_parameter_list().set("Norm abs vol constr", fsimono.get<double>("CONVTOL"));

  //-----------------------------------------------------------------------------
  // ordinary fsi coupling
  //-----------------------------------------------------------------------------

  // right now we use matching meshes at the interface

  Core::Adapter::Coupling& coupsf = structure_fluid_coupling();
  Core::Adapter::Coupling& coupsa = structure_ale_coupling();
  Core::Adapter::Coupling& coupfa = fluid_ale_coupling();

  // structure to fluid
  const int ndim = Global::Problem::Instance()->NDim();
  coupsf.setup_condition_coupling(*structure_field()->discretization(),
      structure_field()->Interface()->fsi_cond_map(), *fluid_field()->discretization(),
      fluid_field()->Interface()->fsi_cond_map(), "FSICoupling", ndim);

  // structure to ale

  coupsa.setup_condition_coupling(*structure_field()->discretization(),
      structure_field()->Interface()->fsi_cond_map(), *ale_field()->discretization(),
      ale_field()->Interface()->fsi_cond_map(), "FSICoupling", ndim);

  // fluid to ale at the interface

  icoupfa_->setup_condition_coupling(*fluid_field()->discretization(),
      fluid_field()->Interface()->fsi_cond_map(), *ale_field()->discretization(),
      ale_field()->Interface()->fsi_cond_map(), "FSICoupling", ndim);

  // In the following we assume that both couplings find the same dof
  // map at the structural side. This enables us to use just one
  // interface dof map for all fields and have just one transfer
  // operator from the interface map to the full field map.
  if (not coupsf.MasterDofMap()->SameAs(*coupsa.MasterDofMap()))
    FOUR_C_THROW("structure interface dof maps do not match");

  if (coupsf.MasterDofMap()->NumGlobalElements() == 0)
    FOUR_C_THROW("No nodes in matching FSI interface. Empty FSI coupling condition?");

  // the fluid-ale coupling always matches
  const Epetra_Map* fluidnodemap = fluid_field()->discretization()->NodeRowMap();
  const Epetra_Map* alenodemap = ale_field()->discretization()->NodeRowMap();

  coupfa.setup_coupling(*fluid_field()->discretization(), *ale_field()->discretization(),
      *fluidnodemap, *alenodemap, ndim);

  fluid_field()->SetMeshMap(coupfa.MasterDofMap());

  aleresidual_ = Teuchos::rcp(new Epetra_Vector(*ale_field()->Interface()->Map(0)));

  // ---------------------------------------------------------------------------
  // Build the global Dirichlet map extractor
  //
  // Dirichlet maps for structure and fluid do not intersect with interface map.
  // ALE Dirichlet map might intersect with interface map, but ALE interface DOFs
  // are not part of the final system of equations. Hence, we just need the
  // intersection of inner ALE DOFs with Dirichlet ALE DOFs.
  std::vector<Teuchos::RCP<const Epetra_Map>> aleintersectionmaps;
  aleintersectionmaps.push_back(ale_field()->get_dbc_map_extractor()->cond_map());
  aleintersectionmaps.push_back(ale_field()->Interface()->other_map());
  Teuchos::RCP<Epetra_Map> aleintersectionmap =
      Core::LinAlg::MultiMapExtractor::intersect_maps(aleintersectionmaps);

  // Merge Dirichlet maps of structure, fluid and ALE to global FSI Dirichlet map
  std::vector<Teuchos::RCP<const Epetra_Map>> dbcmaps;
  dbcmaps.push_back(structure_field()->get_dbc_map_extractor()->cond_map());
  dbcmaps.push_back(fluid_field()->get_dbc_map_extractor()->cond_map());
  dbcmaps.push_back(aleintersectionmap);
  Teuchos::RCP<const Epetra_Map> dbcmap = Core::LinAlg::MultiMapExtractor::merge_maps(dbcmaps);

  // Finally, create the global FSI Dirichlet map extractor
  dbcmaps_ = Teuchos::rcp(new Core::LinAlg::MapExtractor(*dof_row_map(), dbcmap, true));
  if (dbcmaps_ == Teuchos::null)
  {
    FOUR_C_THROW("Creation of FSI Dirichlet map extractor failed.");
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::ConstrMonolithic::evaluate(Teuchos::RCP<const Epetra_Vector> step_increment)
{
  //-----------------------------------------------------------------------------
  // Increment lagrange multiplier
  //-----------------------------------------------------------------------------
  if (step_increment != Teuchos::null)
  {
    Teuchos::RCP<Epetra_Vector> lagrincr = extractor().extract_vector(step_increment, 3);
    conman_->UpdateTotLagrMult(lagrincr);
  }

  //-----------------------------------------------------------------------------
  // evaluation of all fields; constraints are evaluated by strucuture
  //-----------------------------------------------------------------------------
  FSI::Monolithic::evaluate(step_increment);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::ConstrMonolithic::scale_system(Core::LinAlg::BlockSparseMatrixBase& mat, Epetra_Vector& b)
{
  // should we scale the system?
  const Teuchos::ParameterList& fsidyn = Global::Problem::Instance()->FSIDynamicParams();
  const Teuchos::ParameterList& fsimono = fsidyn.sublist("MONOLITHIC SOLVER");
  const bool scaling_infnorm = (bool)Core::UTILS::IntegralValue<int>(fsimono, "INFNORMSCALING");

  if (scaling_infnorm)
  {
    // The matrices are modified here. Do we have to change them back later on?

    Teuchos::RCP<Epetra_CrsMatrix> A = mat.Matrix(0, 0).EpetraMatrix();
    srowsum_ = Teuchos::rcp(new Epetra_Vector(A->RowMap(), false));
    scolsum_ = Teuchos::rcp(new Epetra_Vector(A->RowMap(), false));
    A->InvRowSums(*srowsum_);
    A->InvColSums(*scolsum_);
    if (A->LeftScale(*srowsum_) or A->RightScale(*scolsum_) or
        mat.Matrix(0, 1).EpetraMatrix()->LeftScale(*srowsum_) or
        mat.Matrix(0, 2).EpetraMatrix()->LeftScale(*srowsum_) or
        mat.Matrix(0, 3).EpetraMatrix()->LeftScale(*srowsum_) or
        mat.Matrix(1, 0).EpetraMatrix()->RightScale(*scolsum_) or
        mat.Matrix(2, 0).EpetraMatrix()->RightScale(*scolsum_) or
        mat.Matrix(3, 0).EpetraMatrix()->RightScale(*scolsum_))
      FOUR_C_THROW("structure scaling failed");

    A = mat.Matrix(2, 2).EpetraMatrix();
    arowsum_ = Teuchos::rcp(new Epetra_Vector(A->RowMap(), false));
    acolsum_ = Teuchos::rcp(new Epetra_Vector(A->RowMap(), false));
    A->InvRowSums(*arowsum_);
    A->InvColSums(*acolsum_);
    if (A->LeftScale(*arowsum_) or A->RightScale(*acolsum_) or
        mat.Matrix(2, 0).EpetraMatrix()->LeftScale(*arowsum_) or
        mat.Matrix(2, 1).EpetraMatrix()->LeftScale(*arowsum_) or
        mat.Matrix(0, 2).EpetraMatrix()->RightScale(*acolsum_) or
        mat.Matrix(1, 2).EpetraMatrix()->RightScale(*acolsum_))
      FOUR_C_THROW("ale scaling failed");

    Teuchos::RCP<Epetra_Vector> sx = extractor().extract_vector(b, 0);
    Teuchos::RCP<Epetra_Vector> ax = extractor().extract_vector(b, 2);

    if (sx->Multiply(1.0, *srowsum_, *sx, 0.0)) FOUR_C_THROW("structure scaling failed");
    if (ax->Multiply(1.0, *arowsum_, *ax, 0.0)) FOUR_C_THROW("ale scaling failed");

    extractor().insert_vector(*sx, 0, b);
    extractor().insert_vector(*ax, 2, b);
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::ConstrMonolithic::unscale_solution(
    Core::LinAlg::BlockSparseMatrixBase& mat, Epetra_Vector& x, Epetra_Vector& b)
{
  const Teuchos::ParameterList& fsidyn = Global::Problem::Instance()->FSIDynamicParams();
  const Teuchos::ParameterList& fsimono = fsidyn.sublist("MONOLITHIC SOLVER");
  const bool scaling_infnorm = (bool)Core::UTILS::IntegralValue<int>(fsimono, "INFNORMSCALING");

  if (scaling_infnorm)
  {
    Teuchos::RCP<Epetra_Vector> sy = extractor().extract_vector(x, 0);
    Teuchos::RCP<Epetra_Vector> ay = extractor().extract_vector(x, 2);

    if (sy->Multiply(1.0, *scolsum_, *sy, 0.0)) FOUR_C_THROW("structure scaling failed");
    if (ay->Multiply(1.0, *acolsum_, *ay, 0.0)) FOUR_C_THROW("ale scaling failed");

    extractor().insert_vector(*sy, 0, x);
    extractor().insert_vector(*ay, 2, x);

    Teuchos::RCP<Epetra_Vector> sx = extractor().extract_vector(b, 0);
    Teuchos::RCP<Epetra_Vector> ax = extractor().extract_vector(b, 2);

    if (sx->ReciprocalMultiply(1.0, *srowsum_, *sx, 0.0)) FOUR_C_THROW("structure scaling failed");
    if (ax->ReciprocalMultiply(1.0, *arowsum_, *ax, 0.0)) FOUR_C_THROW("ale scaling failed");

    extractor().insert_vector(*sx, 0, b);
    extractor().insert_vector(*ax, 2, b);

    Teuchos::RCP<Epetra_CrsMatrix> A = mat.Matrix(0, 0).EpetraMatrix();
    srowsum_->Reciprocal(*srowsum_);
    scolsum_->Reciprocal(*scolsum_);
    if (A->LeftScale(*srowsum_) or A->RightScale(*scolsum_) or
        mat.Matrix(0, 1).EpetraMatrix()->LeftScale(*srowsum_) or
        mat.Matrix(0, 2).EpetraMatrix()->LeftScale(*srowsum_) or
        mat.Matrix(0, 3).EpetraMatrix()->LeftScale(*srowsum_) or
        mat.Matrix(1, 0).EpetraMatrix()->RightScale(*scolsum_) or
        mat.Matrix(2, 0).EpetraMatrix()->RightScale(*scolsum_) or
        mat.Matrix(3, 0).EpetraMatrix()->RightScale(*scolsum_))
      FOUR_C_THROW("structure scaling failed");

    A = mat.Matrix(2, 2).EpetraMatrix();
    arowsum_->Reciprocal(*arowsum_);
    acolsum_->Reciprocal(*acolsum_);
    if (A->LeftScale(*arowsum_) or A->RightScale(*acolsum_) or
        mat.Matrix(2, 0).EpetraMatrix()->LeftScale(*arowsum_) or
        mat.Matrix(2, 1).EpetraMatrix()->LeftScale(*arowsum_) or
        mat.Matrix(0, 2).EpetraMatrix()->RightScale(*acolsum_) or
        mat.Matrix(1, 2).EpetraMatrix()->RightScale(*acolsum_))
      FOUR_C_THROW("ale scaling failed");
  }

  // very simple hack just to see the linear solution

  Epetra_Vector r(b.Map());
  mat.Apply(x, r);
  r.Update(1., b, 1.);

  Teuchos::RCP<Epetra_Vector> sr = extractor().extract_vector(r, 0);
  Teuchos::RCP<Epetra_Vector> fr = extractor().extract_vector(r, 1);
  Teuchos::RCP<Epetra_Vector> ar = extractor().extract_vector(r, 2);

  // increment additional ale residual
  aleresidual_->Update(-1., *ar, 0.);

  std::ios_base::fmtflags flags = utils()->out().flags();

  double n, ns, nf, na;
  r.Norm2(&n);
  sr->Norm2(&ns);
  fr->Norm2(&nf);
  ar->Norm2(&na);
  utils()->out() << std::scientific << "\nlinear solver quality:\n"
                 << "L_2-norms:\n"
                 << "   |r|=" << n << "   |rs|=" << ns << "   |rf|=" << nf << "   |ra|=" << na
                 << "\n";
  r.NormInf(&n);
  sr->NormInf(&ns);
  fr->NormInf(&nf);
  ar->NormInf(&na);
  utils()->out() << "L_inf-norms:\n"
                 << "   |r|=" << n << "   |rs|=" << ns << "   |rf|=" << nf << "   |ra|=" << na
                 << "\n";

  utils()->out().flags(flags);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<::NOX::Epetra::LinearSystem> FSI::ConstrMonolithic::create_linear_system(
    Teuchos::ParameterList& nlParams, ::NOX::Epetra::Vector& noxSoln,
    Teuchos::RCP<::NOX::Utils> utils)
{
  Teuchos::RCP<::NOX::Epetra::LinearSystem> linSys;

  Teuchos::ParameterList& printParams = nlParams.sublist("Printing");
  Teuchos::ParameterList& dirParams = nlParams.sublist("Direction");
  Teuchos::ParameterList& newtonParams = dirParams.sublist("Newton");
  Teuchos::ParameterList* lsParams = nullptr;

  // in case of nonlinCG the linear solver list is somewhere else
  if (dirParams.get("Method", "User Defined") == "User Defined")
    lsParams = &(newtonParams.sublist("Linear Solver"));
  else if (dirParams.get("Method", "User Defined") == "NonlinearCG")
    lsParams = &(dirParams.sublist("Nonlinear CG").sublist("Linear Solver"));
  else
    FOUR_C_THROW("Unknown nonlinear method");

  ::NOX::Epetra::Interface::Jacobian* iJac = this;
  ::NOX::Epetra::Interface::Preconditioner* iPrec = this;
  const Teuchos::RCP<Epetra_Operator> J = systemmatrix_;
  const Teuchos::RCP<Epetra_Operator> M = systemmatrix_;

  switch (linearsolverstrategy_)
  {
    case Inpar::FSI::PreconditionedKrylov:
      linSys = Teuchos::rcp(new ::NOX::Epetra::LinearSystemAztecOO(printParams, *lsParams,
          Teuchos::rcp(iJac, false), J, Teuchos::rcp(iPrec, false), M, noxSoln));
      break;
    default:
      FOUR_C_THROW("Unsupported type of monolithic solver/preconditioner!");
      break;
  }

  return linSys;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<::NOX::StatusTest::Combo> FSI::ConstrMonolithic::create_status_test(
    Teuchos::ParameterList& nlParams, Teuchos::RCP<::NOX::Epetra::Group> grp)
{
  // Create the convergence tests
  Teuchos::RCP<::NOX::StatusTest::Combo> combo =
      Teuchos::rcp(new ::NOX::StatusTest::Combo(::NOX::StatusTest::Combo::OR));
  Teuchos::RCP<::NOX::StatusTest::Combo> converged =
      Teuchos::rcp(new ::NOX::StatusTest::Combo(::NOX::StatusTest::Combo::AND));

  Teuchos::RCP<::NOX::StatusTest::MaxIters> maxiters =
      Teuchos::rcp(new ::NOX::StatusTest::MaxIters(nlParams.get("Max Iterations", 100)));
  Teuchos::RCP<::NOX::StatusTest::FiniteValue> fv =
      Teuchos::rcp(new ::NOX::StatusTest::FiniteValue);

  combo->addStatusTest(fv);
  combo->addStatusTest(converged);
  combo->addStatusTest(maxiters);

  // require one solve
  converged->addStatusTest(Teuchos::rcp(new NOX::FSI::MinIters(1)));

  // setup tests for structural displacements

  Teuchos::RCP<::NOX::StatusTest::Combo> structcombo =
      Teuchos::rcp(new ::NOX::StatusTest::Combo(::NOX::StatusTest::Combo::OR));

  Teuchos::RCP<NOX::FSI::PartialNormF> structureDisp = Teuchos::rcp(new NOX::FSI::PartialNormF(
      "displacement", extractor(), 0, nlParams.get<double>("Norm abs disp"),
      ::NOX::Abstract::Vector::TwoNorm, NOX::FSI::PartialNormF::Scaled));
  Teuchos::RCP<NOX::FSI::PartialNormUpdate> structureDispUpdate =
      Teuchos::rcp(new NOX::FSI::PartialNormUpdate("displacement update", extractor(), 0,
          nlParams.get<double>("Norm abs disp"), NOX::FSI::PartialNormUpdate::Scaled));

  add_status_test(structureDisp);
  structcombo->addStatusTest(structureDisp);

  converged->addStatusTest(structcombo);

  // setup tests for fluid velocities

  std::vector<Teuchos::RCP<const Epetra_Map>> fluidvel;
  fluidvel.push_back(fluid_field()->InnerVelocityRowMap());
  fluidvel.push_back(Teuchos::null);
  Core::LinAlg::MultiMapExtractor fluidvelextract(*dof_row_map(), fluidvel);

  Teuchos::RCP<::NOX::StatusTest::Combo> fluidvelcombo =
      Teuchos::rcp(new ::NOX::StatusTest::Combo(::NOX::StatusTest::Combo::OR));

  Teuchos::RCP<NOX::FSI::PartialNormF> innerFluidVel = Teuchos::rcp(new NOX::FSI::PartialNormF(
      "velocity", fluidvelextract, 0, nlParams.get<double>("Norm abs vel"),
      ::NOX::Abstract::Vector::TwoNorm, NOX::FSI::PartialNormF::Scaled));
  Teuchos::RCP<NOX::FSI::PartialNormUpdate> innerFluidVelUpdate =
      Teuchos::rcp(new NOX::FSI::PartialNormUpdate("velocity update", fluidvelextract, 0,
          nlParams.get<double>("Norm abs vel"), NOX::FSI::PartialNormUpdate::Scaled));

  add_status_test(innerFluidVel);
  fluidvelcombo->addStatusTest(innerFluidVel);

  converged->addStatusTest(fluidvelcombo);

  // setup tests for fluid pressure

  std::vector<Teuchos::RCP<const Epetra_Map>> fluidpress;
  fluidpress.push_back(fluid_field()->PressureRowMap());
  fluidpress.push_back(Teuchos::null);
  Core::LinAlg::MultiMapExtractor fluidpressextract(*dof_row_map(), fluidpress);

  Teuchos::RCP<::NOX::StatusTest::Combo> fluidpresscombo =
      Teuchos::rcp(new ::NOX::StatusTest::Combo(::NOX::StatusTest::Combo::OR));

  Teuchos::RCP<NOX::FSI::PartialNormF> fluidPress = Teuchos::rcp(new NOX::FSI::PartialNormF(
      "pressure", fluidpressextract, 0, nlParams.get<double>("Norm abs pres"),
      ::NOX::Abstract::Vector::TwoNorm, NOX::FSI::PartialNormF::Scaled));
  Teuchos::RCP<NOX::FSI::PartialNormUpdate> fluidPressUpdate =
      Teuchos::rcp(new NOX::FSI::PartialNormUpdate("pressure update", fluidpressextract, 0,
          nlParams.get<double>("Norm abs pres"), NOX::FSI::PartialNormUpdate::Scaled));

  add_status_test(fluidPress);
  fluidpresscombo->addStatusTest(fluidPress);

  // setup tests for volume constraint

  std::vector<Teuchos::RCP<const Epetra_Map>> volconstr;
  volconstr.push_back(conman_->GetConstraintMap());
  volconstr.push_back(Teuchos::null);
  Core::LinAlg::MultiMapExtractor volconstrextract(*dof_row_map(), volconstr);

  Teuchos::RCP<::NOX::StatusTest::Combo> volconstrcombo =
      Teuchos::rcp(new ::NOX::StatusTest::Combo(::NOX::StatusTest::Combo::OR));

  Teuchos::RCP<NOX::FSI::PartialNormF> VolConstr = Teuchos::rcp(new NOX::FSI::PartialNormF(
      "constraints", volconstrextract, 0, nlParams.get<double>("Norm abs vol constr"),
      ::NOX::Abstract::Vector::TwoNorm, NOX::FSI::PartialNormF::Scaled));
  Teuchos::RCP<NOX::FSI::PartialNormUpdate> VolConstrUpdate =
      Teuchos::rcp(new NOX::FSI::PartialNormUpdate("constraints update", volconstrextract, 0,
          nlParams.get<double>("Norm abs vol constr"), NOX::FSI::PartialNormUpdate::Scaled));

  add_status_test(VolConstr);
  volconstrcombo->addStatusTest(VolConstr);

  converged->addStatusTest(volconstrcombo);

  return combo;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::ConstrMonolithic::create_system_matrix(bool structuresplit)
{
  const Teuchos::ParameterList& fsidyn = Global::Problem::Instance()->FSIDynamicParams();
  const Teuchos::ParameterList& fsimono = fsidyn.sublist("MONOLITHIC SOLVER");

  // get the PCITER from inputfile
  std::vector<int> pciter;
  std::vector<double> pcomega;
  std::vector<int> spciter;
  std::vector<double> spcomega;
  std::vector<int> fpciter;
  std::vector<double> fpcomega;
  std::vector<int> apciter;
  std::vector<double> apcomega;
  {
    std::string word1;
    std::string word2;
    {
      std::istringstream pciterstream(Teuchos::getNumericStringParameter(fsimono, "PCITER"));
      std::istringstream pcomegastream(Teuchos::getNumericStringParameter(fsimono, "PCOMEGA"));
      while (pciterstream >> word1) pciter.push_back(std::atoi(word1.c_str()));
      while (pcomegastream >> word2) pcomega.push_back(std::atof(word2.c_str()));
    }
    {
      std::istringstream pciterstream(Teuchos::getNumericStringParameter(fsimono, "STRUCTPCITER"));
      std::istringstream pcomegastream(
          Teuchos::getNumericStringParameter(fsimono, "STRUCTPCOMEGA"));
      while (pciterstream >> word1) spciter.push_back(std::atoi(word1.c_str()));
      while (pcomegastream >> word2) spcomega.push_back(std::atof(word2.c_str()));
    }
    {
      std::istringstream pciterstream(Teuchos::getNumericStringParameter(fsimono, "FLUIDPCITER"));
      std::istringstream pcomegastream(Teuchos::getNumericStringParameter(fsimono, "FLUIDPCOMEGA"));
      while (pciterstream >> word1) fpciter.push_back(std::atoi(word1.c_str()));
      while (pcomegastream >> word2) fpcomega.push_back(std::atof(word2.c_str()));
    }
    {
      std::istringstream pciterstream(Teuchos::getNumericStringParameter(fsimono, "ALEPCITER"));
      std::istringstream pcomegastream(Teuchos::getNumericStringParameter(fsimono, "ALEPCOMEGA"));
      while (pciterstream >> word1) apciter.push_back(std::atoi(word1.c_str()));
      while (pcomegastream >> word2) apcomega.push_back(std::atof(word2.c_str()));
    }
  }

  //-----------------------------------------------------------------------------
  // create block system matrix
  //-----------------------------------------------------------------------------

  switch (linearsolverstrategy_)
  {
    case Inpar::FSI::PreconditionedKrylov:
      systemmatrix_ = Teuchos::rcp(new ConstrOverlappingBlockMatrix(extractor(), *structure_field(),
          *fluid_field(), *ale_field(), structuresplit,
          Core::UTILS::IntegralValue<int>(fsimono, "SYMMETRICPRECOND"), pcomega[0], pciter[0],
          spcomega[0], spciter[0], fpcomega[0], fpciter[0], apcomega[0], apciter[0]));

      break;
    default:
      FOUR_C_THROW("Unsupported type of monolithic solver/preconditioner!");
      break;
  }
}

FOUR_C_NAMESPACE_CLOSE
