/*----------------------------------------------------------------------*/
/*! \file

 \brief base for porous media monolithic meshtying method

// Masterthesis of h.Willmann under supervision of Anh-Tu Vuong and Johannes Kremheller
// Originates from poro_monolithic

\level 2

 *----------------------------------------------------------------------*/

#include "4C_poroelast_monolithicmeshtying.hpp"

#include "4C_adapter_coupling_poro_mortar.hpp"
#include "4C_adapter_fld_poro.hpp"
#include "4C_adapter_str_fpsiwrapper.hpp"
#include "4C_coupling_adapter.hpp"
#include "4C_global_data.hpp"
#include "4C_linalg_utils_sparse_algebra_manipulation.hpp"
#include "4C_structure_aux.hpp"

FOUR_C_NAMESPACE_OPEN

PoroElast::MonolithicMeshtying::MonolithicMeshtying(const Epetra_Comm& comm,
    const Teuchos::ParameterList& timeparams,
    Teuchos::RCP<Core::LinAlg::MapExtractor> porosity_splitter)
    : Monolithic(comm, timeparams, porosity_splitter), normrhsfactiven_(0.0), tolfres_ncoup_(0.0)
{
  // Initialize mortar adapter for meshtying interface
  mortar_adapter_ = Teuchos::rcp(new Adapter::CouplingPoroMortar(
      Global::Problem::Instance()->NDim(), Global::Problem::Instance()->mortar_coupling_params(),
      Global::Problem::Instance()->contact_dynamic_params(),
      Global::Problem::Instance()->spatial_approximation_type()));

  const int ndim = Global::Problem::Instance()->NDim();
  std::vector<int> coupleddof(ndim, 1);  // 1,1,1 should be in coupleddof
  // coupleddof[ndim]=0; // not necessary because structural discretization is used
  mortar_adapter_->Setup(structure_field()->discretization(), structure_field()->discretization(),
      coupleddof, "Mortar");

  fvelactiverowdofmap_ = Teuchos::rcp(new Core::LinAlg::MultiMapExtractor);

  // mesh tying not yet works for non-matching structure and fluid discretizations
  if (not matchinggrid_)
  {
    FOUR_C_THROW(
        "The coupling algorithm 'poro_monolithicmeshtying' does not yet work for non-matching "
        "discretizations!");
  }
}

void PoroElast::MonolithicMeshtying::SetupSystem() { Monolithic::SetupSystem(); }

void PoroElast::MonolithicMeshtying::evaluate(
    Teuchos::RCP<const Epetra_Vector> iterinc, bool firstiter)
{
  // evaluate monolithic system for newton iterations
  Monolithic::evaluate(iterinc, firstiter);

  // get state vectors to store in contact data container
  Teuchos::RCP<Epetra_Vector> fvel = fluid_structure_coupling().SlaveToMaster(
      fluid_field()->ExtractVelocityPart(fluid_field()->Velnp()));

  // modified pressure vector modfpres is used to get pressure values to mortar/contact integrator.
  // the pressure values will be written on first displacement DOF

  // extract fluid pressures from full fluid state vector
  Teuchos::RCP<const Epetra_Vector> fpres =
      fluid_field()->ExtractPressurePart(fluid_field()->Velnp());
  // initialize modified pressure vector with fluid velocity dof map
  Teuchos::RCP<Epetra_Vector> modfpres =
      Teuchos::rcp(new Epetra_Vector(*fluid_field()->VelocityRowMap(), true));

  const int ndim = Global::Problem::Instance()->NDim();
  int* mygids = fpres->Map().MyGlobalElements();
  double* val = fpres->Values();
  for (int i = 0; i < fpres->MyLength(); ++i)
  {
    int gid = mygids[i] - ndim;
    // copy pressure value into first velocity DOF of the same node
    modfpres->ReplaceGlobalValues(1, &val[i], &gid);
  }
  // convert velocity map to structure displacement map
  modfpres = fluid_structure_coupling().SlaveToMaster(modfpres);

  // for the set_state() methods in EvaluatePoroMt() non const state vectores are needed
  // ->WriteAccess... methods are used (even though we will not change the states ...)
  Teuchos::RCP<Epetra_Vector> svel = structure_field()->WriteAccessVelnp();
  Teuchos::RCP<Epetra_Vector> sdisp = structure_field()->WriteAccessDispnp();

  // for the EvaluatePoroMt() method RCPs on the matrices are needed...
  Teuchos::RCP<Core::LinAlg::SparseMatrix> f =
      Teuchos::rcpFromRef<Core::LinAlg::SparseMatrix>(systemmatrix_->Matrix(1, 1));
  Teuchos::RCP<Core::LinAlg::SparseMatrix> k_fs =
      Teuchos::rcpFromRef<Core::LinAlg::SparseMatrix>(systemmatrix_->Matrix(1, 0));

  Teuchos::RCP<Epetra_Vector> frhs = Extractor()->ExtractVector(rhs_, 1);

  // modify system matrix and rhs for meshtying
  mortar_adapter_->EvaluatePoroMt(fvel, svel, modfpres, sdisp, structure_field()->discretization(),
      f, k_fs, frhs, fluid_structure_coupling(), fluid_field()->dof_row_map());

  // assign modified parts of system matrix into full system matrix
  systemmatrix_->Assign(1, 1, Core::LinAlg::View, *f);
  systemmatrix_->Assign(1, 0, Core::LinAlg::View, *k_fs);

  // assign modified part of RHS vector into full RHS vector
  Extractor()->InsertVector(*frhs, 1, *rhs_);

  // because the mesh tying interface stays the same, the map extractors for a separate convergence
  // check of the mesh tying fluid coupling condition is only build once
  if ((iter_ == 1) and (Step() == 1)) SetupExtractor();
}

void PoroElast::MonolithicMeshtying::Update()
{
  Monolithic::Update();
  mortar_adapter_->UpdatePoroMt();
}

void PoroElast::MonolithicMeshtying::recover_lagrange_multiplier_after_newton_step(
    Teuchos::RCP<const Epetra_Vector> iterinc)
{
  Monolithic::recover_lagrange_multiplier_after_newton_step(iterinc);

  // displacement and fluid velocity & pressure incremental vector
  Teuchos::RCP<const Epetra_Vector> s_iterinc;
  Teuchos::RCP<const Epetra_Vector> f_iterinc;
  extract_field_vectors(iterinc, s_iterinc, f_iterinc);

  // RecoverStructuralLM
  Teuchos::RCP<Epetra_Vector> tmpsx = Teuchos::rcp<Epetra_Vector>(new Epetra_Vector(*s_iterinc));
  Teuchos::RCP<Epetra_Vector> tmpfx = Teuchos::rcp<Epetra_Vector>(new Epetra_Vector(*f_iterinc));

  mortar_adapter_->recover_fluid_lm_poro_mt(tmpsx, tmpfx);
}

void PoroElast::MonolithicMeshtying::build_convergence_norms()
{
  //------------------------------------------------------------ build residual force norms
  normrhs_ = UTILS::calculate_vector_norm(vectornormfres_, rhs_);
  Teuchos::RCP<const Epetra_Vector> rhs_s;

  Teuchos::RCP<const Epetra_Vector> rhs_f;
  Teuchos::RCP<const Epetra_Vector> rhs_fvel;
  // split velocity into part of normal coupling & tangential condition and velocities
  Teuchos::RCP<const Epetra_Vector> rhs_fvel_activen;
  Teuchos::RCP<const Epetra_Vector> rhs_fvel_other;
  Teuchos::RCP<const Epetra_Vector> rhs_fpres;

  // process structure unknowns of the first field (structure)
  rhs_s = Extractor()->ExtractVector(rhs_, 0);

  // process fluid unknowns of the second field
  rhs_f = Extractor()->ExtractVector(rhs_, 1);
  rhs_fvel = fluid_field()->ExtractVelocityPart(rhs_f);
  // now split it
  rhs_fvel_activen = fluid_vel_active_dof_extractor()->ExtractVector(rhs_fvel, 0);
  rhs_fvel_other = fluid_vel_active_dof_extractor()->ExtractVector(rhs_fvel, 1);
  // pressure is treated separately anyway
  rhs_fpres = fluid_field()->ExtractPressurePart(rhs_f);

  if (porosity_dof_)
  {
    FOUR_C_THROW("porosity dof not implemented for poro_monolithicmeshtying");
    // consult method of mother class for further hints how to do this
  }
  else
  {
    normrhsstruct_ = UTILS::calculate_vector_norm(vectornormfres_, rhs_s);
  }

  normrhsfluid_ = UTILS::calculate_vector_norm(vectornormfres_, rhs_f);
  normrhsfluidvel_ = UTILS::calculate_vector_norm(vectornormfres_, rhs_fvel_other);
  // residual norm of normal coupling condition on poro-fluid
  normrhsfactiven_ = UTILS::calculate_vector_norm(vectornormfres_, rhs_fvel_activen);

  normrhsfluidpres_ = UTILS::calculate_vector_norm(vectornormfres_, rhs_fpres);


  //------------------------------------------------------------- build residual increment norms
  // can stay exactly the same because a monolithic scheme with the same increments as without
  // meshtying is used
  iterinc_->Norm2(&norminc_);

  // displacement and fluid velocity & pressure incremental vector
  Teuchos::RCP<const Epetra_Vector> interincs;
  Teuchos::RCP<const Epetra_Vector> interincf;
  Teuchos::RCP<const Epetra_Vector> interincfvel;
  Teuchos::RCP<const Epetra_Vector> interincfpres;
  // process structure unknowns of the first field
  interincs = Extractor()->ExtractVector(iterinc_, 0);
  // process fluid unknowns of the second field
  interincf = Extractor()->ExtractVector(iterinc_, 1);
  interincfvel = fluid_field()->ExtractVelocityPart(interincf);
  interincfpres = fluid_field()->ExtractPressurePart(interincf);

  normincstruct_ = UTILS::calculate_vector_norm(vectornorminc_, interincs);
  normincfluid_ = UTILS::calculate_vector_norm(vectornorminc_, interincf);
  normincfluidvel_ = UTILS::calculate_vector_norm(vectornorminc_, interincfvel);
  normincfluidpres_ = UTILS::calculate_vector_norm(vectornorminc_, interincfpres);
}

void PoroElast::MonolithicMeshtying::SetupExtractor()
{
  // some maps and vectors
  Teuchos::RCP<Epetra_Map> factivenmap;
  Teuchos::RCP<Epetra_Map> factivenmapcomplement;
  std::vector<Teuchos::RCP<const Epetra_Map>> fluidveldofmapvec;

  // get activemap from poro lagrange strategy of the adapter
  factivenmap = mortar_adapter_->GetPoroStrategy()->FluidActiveNDofMap();

  // build the complement part of the map
  factivenmapcomplement = Core::LinAlg::SplitMap(*fluid_field()->VelocityRowMap(), *factivenmap);

  // write things into the vector for ->Setup
  fluidveldofmapvec.emplace_back(factivenmap);
  fluidveldofmapvec.emplace_back(factivenmapcomplement);

  fvelactiverowdofmap_->Setup(*fluid_field()->VelocityRowMap(), fluidveldofmapvec);
}

bool PoroElast::MonolithicMeshtying::Converged()
{
  // check for single norms
  bool convinc = false;
  bool convfres = false;

  // convinc can stay the same because the increments are the same as without meshtying
  // residual increments
  switch (normtypeinc_)
  {
    case Inpar::PoroElast::convnorm_abs_global:
      convinc = norminc_ < tolinc_;
      break;
    case Inpar::PoroElast::convnorm_abs_singlefields:
      convinc = (normincstruct_ < tolinc_struct_ and normincfluidvel_ < tolinc_velocity_ and
                 normincfluidpres_ < tolinc_pressure_ and normincporo_ < tolinc_porosity_);
      break;
    default:
      FOUR_C_THROW("Cannot check for convergence of residual values!");
      break;
  }

  // residual forces
  switch (normtypefres_)
  {
    case Inpar::PoroElast::convnorm_abs_global:
      convfres = normrhs_ < tolfres_;
      break;
    case Inpar::PoroElast::convnorm_abs_singlefields:
      convfres = (normrhsstruct_ < tolfres_struct_ and normrhsfluidvel_ < tolfres_velocity_ and
                  normrhsfluidpres_ < tolfres_pressure_ and normrhsporo_ < tolfres_porosity_ and
                  normrhsfactiven_ < tolfres_ncoup_);
      break;
    default:
      FOUR_C_THROW("Cannot check for convergence of residual forces!");
      break;
  }

  // combine increments and forces
  bool conv = false;
  if (combincfres_ == Inpar::PoroElast::bop_and)
    conv = convinc and convfres;
  else if (combincfres_ == Inpar::PoroElast::bop_or)
    conv = convinc or convfres;
  else
    FOUR_C_THROW("Something went terribly wrong with binary operator!");

  return conv;
}

bool PoroElast::MonolithicMeshtying::SetupSolver()
{
  Monolithic::SetupSolver();

  // get dynamic section of poroelasticity
  const Teuchos::ParameterList& poroelastdyn =
      Global::Problem::Instance()->poroelast_dynamic_params();

  tolfres_ncoup_ = poroelastdyn.get<double>("TOLRES_NCOUP");

  return true;
}

void PoroElast::MonolithicMeshtying::print_newton_iter_header_stream(std::ostringstream& oss)
{
  oss << "------------------------------------------------------------" << std::endl;
  oss << "                   Newton-Raphson Scheme                    " << std::endl;
  oss << "                NormRES " << VectorNormString(vectornormfres_);
  oss << "     NormINC " << VectorNormString(vectornorminc_) << "                    " << std::endl;
  oss << "------------------------------------------------------------" << std::endl;

  // enter converged state etc
  oss << "numiter";

  // different style due relative or absolute error checking

  // --------------------------------------------------------global system test
  // residual forces
  switch (normtypefres_)
  {
    case Inpar::PoroElast::convnorm_abs_global:
      oss << std::setw(15) << "abs-res"
          << "(" << std::setw(5) << std::setprecision(2) << tolfres_ << ")";
      break;
    case Inpar::PoroElast::convnorm_abs_singlefields:
      oss << std::setw(15) << "abs-s-res"
          << "(" << std::setw(5) << std::setprecision(2) << tolfres_struct_ << ")";
      if (porosity_dof_)
        oss << std::setw(15) << "abs-poro-res"
            << "(" << std::setw(5) << std::setprecision(2) << tolfres_porosity_ << ")";
      oss << std::setw(15) << "abs-fvel-res"
          << "(" << std::setw(5) << std::setprecision(2) << tolfres_velocity_ << ")";
      oss << std::setw(15) << "abs-fpres-res"
          << "(" << std::setw(5) << std::setprecision(2) << tolfres_pressure_ << ")";
      oss << std::setw(15) << "abs-fncoup-res"
          << "(" << std::setw(5) << std::setprecision(2) << tolfres_ncoup_ << ")";
      break;
    default:
      FOUR_C_THROW("Unknown or undefined convergence form for residual.");
      break;
  }

  switch (normtypeinc_)
  {
    case Inpar::PoroElast::convnorm_abs_global:
      oss << std::setw(15) << "abs-inc"
          << "(" << std::setw(5) << std::setprecision(2) << tolinc_ << ")";
      break;
    case Inpar::PoroElast::convnorm_abs_singlefields:
      oss << std::setw(15) << "abs-s-inc"
          << "(" << std::setw(5) << std::setprecision(2) << tolinc_struct_ << ")";
      if (porosity_dof_)
        oss << std::setw(15) << "abs-poro-inc"
            << "(" << std::setw(5) << std::setprecision(2) << tolinc_porosity_ << ")";
      oss << std::setw(15) << "abs-fvel-inc"
          << "(" << std::setw(5) << std::setprecision(2) << tolinc_velocity_ << ")";
      oss << std::setw(15) << "abs-fpres-inc"
          << "(" << std::setw(5) << std::setprecision(2) << tolinc_pressure_ << ")";
      break;
    default:
      FOUR_C_THROW("Unknown or undefined convergence form for increment.");
      break;
  }
}

void PoroElast::MonolithicMeshtying::print_newton_iter_text_stream(std::ostringstream& oss)
{
  // enter converged state etc
  oss << std::setw(7) << iter_;

  // different style due relative or absolute error checking

  // --------------------------------------------------------global system test
  // residual forces
  switch (normtypefres_)
  {
    case Inpar::PoroElast::convnorm_abs_global:
      oss << std::setw(22) << std::setprecision(5) << std::scientific << normrhs_;
      break;
    case Inpar::PoroElast::convnorm_abs_singlefields:
      break;
    default:
      FOUR_C_THROW("Unknown or undefined convergence form for global residual.");
      break;
  }
  // increments
  switch (normtypeinc_)
  {
    case Inpar::PoroElast::convnorm_abs_global:
      oss << std::setw(22) << std::setprecision(5) << std::scientific << norminc_;
      break;
    case Inpar::PoroElast::convnorm_abs_singlefields:
      break;
    default:
      FOUR_C_THROW("Unknown or undefined convergence form for global increment.");
      break;
  }

  // --------------------------------------------------------single field test
  switch (normtypefres_)
  {
    case Inpar::PoroElast::convnorm_abs_singlefields:
      oss << std::setw(22) << std::setprecision(5) << std::scientific << normrhsstruct_;
      if (porosity_dof_)
        oss << std::setw(22) << std::setprecision(5) << std::scientific << normrhsporo_;
      oss << std::setw(22) << std::setprecision(5) << std::scientific << normrhsfluidvel_;
      oss << std::setw(22) << std::setprecision(5) << std::scientific << normrhsfluidpres_;
      oss << std::setw(22) << std::setprecision(5) << std::scientific << normrhsfactiven_;
      break;
    case Inpar::PoroElast::convnorm_abs_global:
      break;
    default:
      FOUR_C_THROW("Unknown or undefined convergence form for single field residual.");
      break;
  }

  switch (normtypeinc_)
  {
    case Inpar::PoroElast::convnorm_abs_singlefields:
      oss << std::setw(22) << std::setprecision(5) << std::scientific << normincstruct_;
      if (porosity_dof_)
        oss << std::setw(22) << std::setprecision(5) << std::scientific << normincporo_;
      oss << std::setw(22) << std::setprecision(5) << std::scientific << normincfluidvel_;
      oss << std::setw(22) << std::setprecision(5) << std::scientific << normincfluidpres_;
      break;
    case Inpar::PoroElast::convnorm_abs_global:
      break;
    default:
      FOUR_C_THROW("Unknown or undefined convergence form for single field increment.");
      break;
  }
}

FOUR_C_NAMESPACE_CLOSE
