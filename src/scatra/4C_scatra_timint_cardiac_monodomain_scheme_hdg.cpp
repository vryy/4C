/*----------------------------------------------------------------------*/
/*! \file
\brief time-integration scheme for HDG with extensions for
       cardiac monodomain problems

\level 3


*/
/*----------------------------------------------------------------------*/

#include "4C_scatra_timint_cardiac_monodomain_scheme_hdg.hpp"

#include "4C_global_data.hpp"
#include "4C_io.hpp"
#include "4C_scatra_ele_action.hpp"
#include "4C_scatra_ele_calc_hdg.hpp"
#include "4C_utils_parameter_list.hpp"

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
ScaTra::TimIntCardiacMonodomainHDG::TimIntCardiacMonodomainHDG(
    Teuchos::RCP<Core::FE::Discretization> actdis, Teuchos::RCP<Core::LinAlg::Solver> solver,
    Teuchos::RCP<Teuchos::ParameterList> params,
    Teuchos::RCP<Teuchos::ParameterList> sctratimintparams,
    Teuchos::RCP<Teuchos::ParameterList> extraparams,
    Teuchos::RCP<Core::IO::DiscretizationWriter> output)
    : ScaTraTimIntImpl(actdis, solver, sctratimintparams, extraparams, output),
      TimIntCardiacMonodomain(actdis, solver, params, sctratimintparams, extraparams, output),
      TimIntHDG(actdis, solver, sctratimintparams, extraparams, output)
{
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ScaTra::TimIntCardiacMonodomainHDG::setup()
{
  // call init()-functions of base classes
  // note: this order is important
  TimIntHDG::setup();
  TimIntCardiacMonodomain::setup();

  // Activation time at time n+1
  activation_time_interpol_.reset(new Epetra_Vector(*discret_->NodeRowMap()));
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ScaTra::TimIntCardiacMonodomainHDG::Update()
{
  // time update of myocard material
  element_material_time_update();

  // Standard Update
  TimIntHDG::Update();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ScaTra::TimIntCardiacMonodomainHDG::element_material_time_update()
{
  discret_->ClearState(true);

  Teuchos::ParameterList eleparams;
  Core::UTILS::AddEnumClassToParameterList<ScaTra::Action>(
      "action", ScaTra::Action::time_update_material, eleparams);

  discret_->set_state("phiaf", phinp_);
  discret_->set_state(nds_intvar_, "intphin", intphin_);
  discret_->set_state(0, "phin", phin_);


  Core::LinAlg::SerialDenseMatrix dummyMat;
  Core::LinAlg::SerialDenseVector dummyVec;
  Core::Elements::Element::LocationArray la(discret_->NumDofSets());


  for (int iele = 0; iele < discret_->NumMyColElements(); ++iele)
  {
    Core::Elements::Element *ele = discret_->lColElement(iele);
    ele->LocationVector(*discret_, la, false);

    ele->evaluate(eleparams, *discret_, la, dummyMat, dummyMat, dummyVec, dummyVec, dummyVec);
  }

  discret_->ClearState(true);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ScaTra::TimIntCardiacMonodomainHDG::output_state()
{
  // Call function from base class
  ScaTra::TimIntHDG::output_state();

  if (nb_max_mat_int_state_vars_)
  {
    material_internal_state_np_->PutScalar(0.0);
    Teuchos::ParameterList params;
    Core::UTILS::AddEnumClassToParameterList<ScaTra::Action>(
        "action", ScaTra::Action::get_material_internal_state, params);
    params.set<Teuchos::RCP<Epetra_MultiVector>>(
        "material_internal_state", material_internal_state_np_);
    discret_->evaluate(params);
    material_internal_state_np_ =
        params.get<Teuchos::RCP<Epetra_MultiVector>>("material_internal_state");
    if (material_internal_state_np_ == Teuchos::null)
      FOUR_C_THROW("Cannot get state vector material internal state");

    output_->write_vector("ionic_currents_hdg", material_internal_state_np_);

    for (int k = 0; k < material_internal_state_np_->NumVectors(); ++k)
    {
      std::ostringstream temp;
      temp << k + 1;
      material_internal_state_np_component_ =
          Teuchos::rcp((*material_internal_state_np_)(k), false);
      output_->write_vector("mat_int_state_hdg" + temp.str(), material_internal_state_np_component_,
          Core::IO::elementvector);
    }
  }


  // copy values from node to dof vector
  Teuchos::RCP<Epetra_Vector> dofphi = Core::LinAlg::CreateVector(*discret_->NodeRowMap());

  for (int i = 0; i < dofphi->MyLength(); ++i)
  {
    int dofgid = discret_->NodeRowMap()->GID(i);
    dofphi->ReplaceMyValue(discret_->NodeRowMap()->LID(dofgid), 0, (*interpolatedPhinp_)[i]);
  }
  output_->write_vector("phinp", dofphi);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ScaTra::TimIntCardiacMonodomainHDG::write_problem_specific_output(
    Teuchos::RCP<Epetra_Vector> interpolatedPhi)
{
  // Compute and write activation time
  if (activation_time_interpol_ != Teuchos::null)
  {
    for (int k = 0; k < interpolatedPhi->MyLength(); k++)
    {
      if ((*interpolatedPhi)[k] >= activation_threshold_ &&
          (*activation_time_interpol_)[k] <= dta_ * 0.9)
        (*activation_time_interpol_)[k] = time_;
    }
    output_->write_vector(
        "activation_time_np_hdg", activation_time_interpol_, Core::IO::nodevector);
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ScaTra::TimIntCardiacMonodomainHDG::pack_material()
{
  Core::Communication::PackBuffer buffer;

  // loop over elements
  for (int iele = 0; iele < discret_->NumMyColElements(); ++iele)
  {
    auto *hdgele = dynamic_cast<Discret::ELEMENTS::ScaTraHDG *>(discret_->lColElement(iele));
    hdgele->pack_material(buffer);
  }

  Teuchos::RCP<std::vector<char>> block = Teuchos::rcp(new std::vector<char>);
  std::swap(*block, buffer());
  data_ = block;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ScaTra::TimIntCardiacMonodomainHDG::unpack_material()
{
  std::vector<char>::size_type index = 0;
  // loop over elements
  for (int iele = 0; iele < discret_->NumMyColElements(); ++iele)
  {
    auto *hdgele = dynamic_cast<Discret::ELEMENTS::ScaTraHDG *>(discret_->lColElement(iele));
    std::vector<char> data;
    hdgele->extract_from_pack(index, *data_, data);
    hdgele->unpack_material(data);
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ScaTra::TimIntCardiacMonodomainHDG::project_material()
{
  discret_->ClearState(true);
  // set action
  Teuchos::ParameterList eleparams;
  Core::UTILS::AddEnumClassToParameterList<ScaTra::Action>(
      "action", ScaTra::Action::project_material_field, eleparams);

  Core::LinAlg::SerialDenseMatrix dummyMat;
  Core::LinAlg::SerialDenseVector dummyVec;
  Core::Elements::Element::LocationArray dummy(1);

  for (int iele = 0; iele < discret_->NumMyColElements(); ++iele)
  {
    Core::Elements::Element *ele = discret_->lColElement(iele);

    // call routine on elements to project material field
    ele->evaluate(eleparams, *discret_, dummy, dummyMat, dummyMat, dummyVec, dummyVec, dummyVec);
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ScaTra::TimIntCardiacMonodomainHDG::read_restart(
    const int step, Teuchos::RCP<Core::IO::InputControl> input)
{
  // Call function from base class
  ScaTra::TimIntHDG::read_restart(step, input);

  activation_time_interpol_.reset(new Epetra_Vector(*discret_->NodeRowMap()));
}

FOUR_C_NAMESPACE_CLOSE
