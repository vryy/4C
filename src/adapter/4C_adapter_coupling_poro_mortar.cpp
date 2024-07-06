/*---------------------------------------------------------------------*/
/*! \file

\brief coupling adapter for porous meshtying

Masterthesis of h.Willmann under supervision of Anh-Tu Vuong
and Johannes Kremheller, Originates from Adapter::CouplingNonLinMortar

\level 3


*/
/*---------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 |  includes                                                  ager 10/15|
 *----------------------------------------------------------------------*/
// lib
#include "4C_adapter_coupling_poro_mortar.hpp"

#include "4C_contact_element.hpp"
#include "4C_contact_interface.hpp"
#include "4C_contact_node.hpp"
#include "4C_coupling_adapter.hpp"
#include "4C_fem_discretization.hpp"
#include "4C_fem_nurbs_discretization.hpp"
#include "4C_fem_nurbs_discretization_control_point.hpp"
#include "4C_fem_nurbs_discretization_knotvector.hpp"
#include "4C_global_data.hpp"
#include "4C_inpar_contact.hpp"
#include "4C_inpar_structure.hpp"
#include "4C_linalg_sparsematrix.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 |  ctor                                                      ager 10/15|
 *----------------------------------------------------------------------*/
Adapter::CouplingPoroMortar::CouplingPoroMortar(int spatial_dimension,
    Teuchos::ParameterList mortar_coupling_params, Teuchos::ParameterList contact_dynamic_params,
    Core::FE::ShapeFunctionType shape_function_type)
    : CouplingNonLinMortar(
          spatial_dimension, mortar_coupling_params, contact_dynamic_params, shape_function_type),
      firstinit_(false),
      slavetype_(-1),
      mastertype_(-1)
{
  // empty...
}

/*----------------------------------------------------------------------*
 |  Read Mortar Condition                                     ager 10/15|
 *----------------------------------------------------------------------*/
void Adapter::CouplingPoroMortar::read_mortar_condition(
    Teuchos::RCP<Core::FE::Discretization> masterdis,
    Teuchos::RCP<Core::FE::Discretization> slavedis, std::vector<int> coupleddof,
    const std::string& couplingcond, Teuchos::ParameterList& input,
    std::map<int, Core::Nodes::Node*>& mastergnodes, std::map<int, Core::Nodes::Node*>& slavegnodes,
    std::map<int, Teuchos::RCP<Core::Elements::Element>>& masterelements,
    std::map<int, Teuchos::RCP<Core::Elements::Element>>& slaveelements)
{
  // Call Base Class
  CouplingNonLinMortar::read_mortar_condition(masterdis, slavedis, coupleddof, couplingcond, input,
      mastergnodes, slavegnodes, masterelements, slaveelements);

  // Set Problem Type to Poro
  switch (Global::Problem::instance()->get_problem_type())
  {
    case Core::ProblemType::poroelast:
      input.set<int>("PROBTYPE", Inpar::CONTACT::poroelast);
      break;
    case Core::ProblemType::poroscatra:
      input.set<int>("PROBTYPE", Inpar::CONTACT::poroscatra);
      break;
    default:
      FOUR_C_THROW("Invalid poro problem is specified");
      break;
  }

  // porotimefac = 1/(theta*dt) --- required for derivation of structural displacements!
  const Teuchos::ParameterList& stru = Global::Problem::instance()->structural_dynamic_params();
  double porotimefac =
      1 / (stru.sublist("ONESTEPTHETA").get<double>("THETA") * stru.get<double>("TIMESTEP"));
  input.set<double>("porotimefac", porotimefac);
  const Teuchos::ParameterList& porodyn = Global::Problem::instance()->poroelast_dynamic_params();
  input.set<bool>("CONTACTNOPEN",
      Core::UTILS::IntegralValue<int>(porodyn, "CONTACTNOPEN"));  // used in the integrator
  if (!Core::UTILS::IntegralValue<int>(porodyn, "CONTACTNOPEN"))
    FOUR_C_THROW("Set CONTACTNOPEN for Poroelastic meshtying!");
}

/*----------------------------------------------------------------------*
 |  Add Mortar Elements                                        ager 10/15|
 *----------------------------------------------------------------------*/
void Adapter::CouplingPoroMortar::add_mortar_elements(
    Teuchos::RCP<Core::FE::Discretization> masterdis,
    Teuchos::RCP<Core::FE::Discretization> slavedis, Teuchos::ParameterList& input,
    std::map<int, Teuchos::RCP<Core::Elements::Element>>& masterelements,
    std::map<int, Teuchos::RCP<Core::Elements::Element>>& slaveelements,
    Teuchos::RCP<CONTACT::Interface>& interface, int numcoupleddof)
{
  bool isnurbs = input.get<bool>("NURBS");

  // get problem dimension (2D or 3D) and create (Mortar::Interface)
  const int dim = Global::Problem::instance()->n_dim();

  // We need to determine an element offset to start the numbering of the slave
  // mortar elements AFTER the master mortar elements in order to ensure unique
  // eleIDs in the interface discretization. The element offset equals the
  // overall number of master mortar elements (which is not equal to the number
  // of elements in the field that is chosen as master side).
  //
  // If masterdis==slavedis, the element numbering is right without offset
  int eleoffset = 0;
  if (masterdis.get() != slavedis.get())
  {
    int nummastermtreles = masterelements.size();
    comm_->SumAll(&nummastermtreles, &eleoffset, 1);
  }

  // feeding master elements to the interface
  std::map<int, Teuchos::RCP<Core::Elements::Element>>::const_iterator elemiter;
  for (elemiter = masterelements.begin(); elemiter != masterelements.end(); ++elemiter)
  {
    Teuchos::RCP<Core::Elements::Element> ele = elemiter->second;
    Teuchos::RCP<CONTACT::Element> cele = Teuchos::rcp(new CONTACT::Element(
        ele->id(), ele->owner(), ele->shape(), ele->num_node(), ele->node_ids(), false, isnurbs));

    Teuchos::RCP<Core::Elements::FaceElement> faceele =
        Teuchos::rcp_dynamic_cast<Core::Elements::FaceElement>(ele, true);
    if (faceele == Teuchos::null) FOUR_C_THROW("Cast to FaceElement failed!");
    cele->phys_type() = Mortar::Element::other;

    std::vector<Teuchos::RCP<Core::Conditions::Condition>> porocondvec;
    masterdis->get_condition("PoroCoupling", porocondvec);
    for (unsigned int i = 0; i < porocondvec.size(); ++i)
    {
      std::map<int, Teuchos::RCP<Core::Elements::Element>>::const_iterator eleitergeometry;
      for (eleitergeometry = porocondvec[i]->geometry().begin();
           eleitergeometry != porocondvec[i]->geometry().end(); ++eleitergeometry)
      {
        if (faceele->parent_element()->id() == eleitergeometry->second->id())
        {
          if (mastertype_ == 0)
            FOUR_C_THROW(
                "struct and poro master elements on the same processor - no mixed interface "
                "supported");
          cele->phys_type() = Mortar::Element::poro;
          mastertype_ = 1;
          break;
        }
      }
    }
    if (cele->phys_type() == Mortar::Element::other)
    {
      if (mastertype_ == 1)
        FOUR_C_THROW(
            "struct and poro master elements on the same processor - no mixed interface supported");
      cele->phys_type() = Mortar::Element::structure;
      mastertype_ = 0;
    }

    cele->set_parent_master_element(faceele->parent_element(), faceele->face_parent_number());

    if (isnurbs)
    {
      Teuchos::RCP<Core::FE::Nurbs::NurbsDiscretization> nurbsdis =
          Teuchos::rcp_dynamic_cast<Core::FE::Nurbs::NurbsDiscretization>(masterdis);

      Teuchos::RCP<Core::FE::Nurbs::Knotvector> knots = (*nurbsdis).get_knot_vector();
      std::vector<Core::LinAlg::SerialDenseVector> parentknots(dim);
      std::vector<Core::LinAlg::SerialDenseVector> mortarknots(dim - 1);

      Teuchos::RCP<Core::Elements::FaceElement> faceele =
          Teuchos::rcp_dynamic_cast<Core::Elements::FaceElement>(ele, true);
      double normalfac = 0.0;
      bool zero_size = knots->get_boundary_ele_and_parent_knots(parentknots, mortarknots, normalfac,
          faceele->parent_master_element()->id(), faceele->face_master_number());

      // store nurbs specific data to node
      cele->zero_sized() = zero_size;
      cele->knots() = mortarknots;
      cele->normal_fac() = normalfac;
    }

    interface->add_element(cele);
  }

  // feeding slave elements to the interface
  for (elemiter = slaveelements.begin(); elemiter != slaveelements.end(); ++elemiter)
  {
    Teuchos::RCP<Core::Elements::Element> ele = elemiter->second;

    Teuchos::RCP<CONTACT::Element> cele = Teuchos::rcp(new CONTACT::Element(
        ele->id(), ele->owner(), ele->shape(), ele->num_node(), ele->node_ids(), true, isnurbs));

    Teuchos::RCP<Core::Elements::FaceElement> faceele =
        Teuchos::rcp_dynamic_cast<Core::Elements::FaceElement>(ele, true);
    if (faceele == Teuchos::null) FOUR_C_THROW("Cast to FaceElement failed!");
    cele->phys_type() = Mortar::Element::other;

    std::vector<Teuchos::RCP<Core::Conditions::Condition>> porocondvec;
    masterdis->get_condition("PoroCoupling", porocondvec);

    for (unsigned int i = 0; i < porocondvec.size(); ++i)
    {
      std::map<int, Teuchos::RCP<Core::Elements::Element>>::const_iterator eleitergeometry;
      for (eleitergeometry = porocondvec[i]->geometry().begin();
           eleitergeometry != porocondvec[i]->geometry().end(); ++eleitergeometry)
      {
        if (faceele->parent_element()->id() == eleitergeometry->second->id())
        {
          if (slavetype_ == 0)
            FOUR_C_THROW(
                "struct and poro slave elements on the same processor - no mixed interface "
                "supported");
          cele->phys_type() = Mortar::Element::poro;
          slavetype_ = 1;
          break;
        }
      }
    }
    if (cele->phys_type() == Mortar::Element::other)
    {
      if (slavetype_ == 1)
        FOUR_C_THROW(
            "struct and poro slave elements on the same processor - no mixed interface supported");
      cele->phys_type() = Mortar::Element::structure;
      slavetype_ = 0;
    }
    cele->set_parent_master_element(faceele->parent_element(), faceele->face_parent_number());

    if (isnurbs)
    {
      Teuchos::RCP<Core::FE::Nurbs::NurbsDiscretization> nurbsdis =
          Teuchos::rcp_dynamic_cast<Core::FE::Nurbs::NurbsDiscretization>(slavedis);

      Teuchos::RCP<Core::FE::Nurbs::Knotvector> knots = (*nurbsdis).get_knot_vector();
      std::vector<Core::LinAlg::SerialDenseVector> parentknots(dim);
      std::vector<Core::LinAlg::SerialDenseVector> mortarknots(dim - 1);

      Teuchos::RCP<Core::Elements::FaceElement> faceele =
          Teuchos::rcp_dynamic_cast<Core::Elements::FaceElement>(ele, true);
      double normalfac = 0.0;
      bool zero_size = knots->get_boundary_ele_and_parent_knots(parentknots, mortarknots, normalfac,
          faceele->parent_master_element()->id(), faceele->face_master_number());

      // store nurbs specific data to node
      cele->zero_sized() = zero_size;
      cele->knots() = mortarknots;
      cele->normal_fac() = normalfac;
    }

    interface->add_element(cele);
  }

  return;
}


/*----------------------------------------------------------------------*
 |  read mortar condition                                    Ager 02/16 |
 *----------------------------------------------------------------------*/
void Adapter::CouplingPoroMortar::create_strategy(Teuchos::RCP<Core::FE::Discretization> masterdis,
    Teuchos::RCP<Core::FE::Discretization> slavedis, Teuchos::ParameterList& input,
    int numcoupleddof)
{
  // poro lagrange strategy:

  // get problem dimension (2D or 3D) and create (Mortar::Interface)
  const int dim = Global::Problem::instance()->n_dim();

  // bools to decide which side is structural and which side is poroelastic to manage all 4
  // constellations
  // s-s, p-s, s-p, p-p
  bool poromaster = false;
  bool poroslave = false;
  bool structmaster = false;
  bool structslave = false;

  // wait for all processors to determine if they have poro or structural master or slave elements
  comm_->Barrier();
  std::vector<int> slaveTypeList(comm_->NumProc());
  std::vector<int> masterTypeList(comm_->NumProc());
  comm_->GatherAll(&slavetype_, slaveTypeList.data(), 1);
  comm_->GatherAll(&mastertype_, masterTypeList.data(), 1);
  comm_->Barrier();

  for (int i = 0; i < comm_->NumProc(); ++i)
  {
    switch (slaveTypeList[i])
    {
      case -1:
        break;
      case 1:
        if (structslave)
          FOUR_C_THROW(
              "struct and poro slave elements on the same adapter - no mixed interface supported");
        // adjust FOUR_C_THROW text, when more than one interface is supported
        poroslave = true;
        break;
      case 0:
        if (poroslave)
          FOUR_C_THROW(
              "struct and poro slave elements on the same adapter - no mixed interface supported");
        structslave = true;
        break;
      default:
        FOUR_C_THROW("this cannot happen");
        break;
    }
  }

  for (int i = 0; i < comm_->NumProc(); ++i)
  {
    switch (masterTypeList[i])
    {
      case -1:
        break;
      case 1:
        if (structmaster)
          FOUR_C_THROW(
              "struct and poro master elements on the same adapter - no mixed interface supported");
        // adjust FOUR_C_THROW text, when more than one interface is supported
        poromaster = true;
        break;
      case 0:
        if (poromaster)
          FOUR_C_THROW(
              "struct and poro master elements on the same adapter - no mixed interface supported");
        structmaster = true;
        break;
      default:
        FOUR_C_THROW("this cannot happen");
        break;
    }
  }

  const Teuchos::ParameterList& stru = Global::Problem::instance()->structural_dynamic_params();
  double theta = stru.sublist("ONESTEPTHETA").get<double>("THETA");
  // what if problem is static ? there should be an error for previous line called in a dyna_statics
  // problem and not a value of 0.5 a proper disctinction is necessary if poro meshtying is expanded
  // to other time integration strategies

  if (Core::UTILS::IntegralValue<Inpar::Solid::DynamicType>(stru, "DYNAMICTYP") ==
      Inpar::Solid::dyna_statics)
  {
    theta = 1.0;
  }
  std::vector<Teuchos::RCP<CONTACT::Interface>> interfaces;
  interfaces.push_back(interface_);
  double alphaf = 1.0 - theta;

  // build the correct data container
  Teuchos::RCP<CONTACT::AbstractStratDataContainer> data_ptr =
      Teuchos::rcp(new CONTACT::AbstractStratDataContainer());
  // create contact poro lagrange strategy for mesh tying
  porolagstrategy_ = Teuchos::rcp(new CONTACT::LagrangeStrategyPoro(data_ptr,
      masterdis->dof_row_map(), masterdis->node_row_map(), input, interfaces, dim, comm_, alphaf,
      numcoupleddof, poroslave, poromaster));

  porolagstrategy_->setup(false, true);
  porolagstrategy_->poro_mt_initialize();

  firstinit_ = true;

  return;
}


/*----------------------------------------------------------------------*
 |  complete interface (also print and parallel redist.)      Ager 02/16|
 *----------------------------------------------------------------------*/
void Adapter::CouplingPoroMortar::complete_interface(
    Teuchos::RCP<Core::FE::Discretization> masterdis, Teuchos::RCP<CONTACT::Interface>& interface)
{
  // finalize the contact interface construction
  int maxdof = masterdis->dof_row_map()->MaxAllGID();
  interface->fill_complete(true, maxdof);

  // interface->create_volume_ghosting(*masterdis);

  // store old row maps (before parallel redistribution)
  slavedofrowmap_ = Teuchos::rcp(new Epetra_Map(*interface->slave_row_dofs()));
  masterdofrowmap_ = Teuchos::rcp(new Epetra_Map(*interface->master_row_dofs()));
  slavenoderowmap_ = Teuchos::rcp(new Epetra_Map(*interface->slave_row_nodes()));

  // print parallel distribution
  interface->print_parallel_distribution();

  // TODO: is this possible?
  //**********************************************************************
  // PARALLEL REDISTRIBUTION OF INTERFACE
  //**********************************************************************
  //  if (parredist && comm_->NumProc()>1)
  //  {
  //    // redistribute optimally among all procs
  //    interface->Redistribute(1);
  //
  //    // call fill complete again
  //    interface->fill_complete();
  //
  //    // print parallel distribution again
  //    interface->print_parallel_distribution();
  //  }

  // store interface
  interface_ = interface;

  // create binary search tree
  interface_->create_search_tree();

  return;
}


/*----------------------------------------------------------------------*
 |  evaluate blockmatrices for poro meshtying             ager 10/15    |
 *----------------------------------------------------------------------*/
void Adapter::CouplingPoroMortar::evaluate_poro_mt(Teuchos::RCP<Epetra_Vector> fvel,
    Teuchos::RCP<Epetra_Vector> svel, Teuchos::RCP<Epetra_Vector> fpres,
    Teuchos::RCP<Epetra_Vector> sdisp, const Teuchos::RCP<Core::FE::Discretization> sdis,
    Teuchos::RCP<Core::LinAlg::SparseMatrix>& f, Teuchos::RCP<Core::LinAlg::SparseMatrix>& k_fs,
    Teuchos::RCP<Epetra_Vector>& frhs, Core::Adapter::Coupling& coupfs,
    Teuchos::RCP<const Epetra_Map> fdofrowmap)
{
  // safety check
  check_setup();

  // write interface values into poro contact interface data containers for integration in contact
  // integrator
  porolagstrategy_->set_state(Mortar::state_fvelocity, *fvel);
  porolagstrategy_->set_state(Mortar::state_svelocity, *svel);
  porolagstrategy_->set_state(Mortar::state_fpressure, *fpres);
  porolagstrategy_->set_state(Mortar::state_new_displacement, *sdisp);

  // store displacements of parent elements for deformation gradient determinant and its
  // linearization
  porolagstrategy_->set_parent_state("displacement", sdisp, sdis);

  interface_->initialize();
  // in the end of Evaluate coupling condition residuals and linearizations are computed in contact
  // integrator
  interface_->evaluate();

  porolagstrategy_->poro_mt_prepare_fluid_coupling();
  porolagstrategy_->poro_initialize(coupfs, fdofrowmap, firstinit_);
  if (firstinit_) firstinit_ = false;

  // do system matrix manipulations
  porolagstrategy_->evaluate_poro_no_pen_contact(k_fs, f, frhs);
  return;
}  // Adapter::CouplingNonLinMortar::EvaluatePoroMt()


/*----------------------------------------------------------------------*
 |  update poro meshtying quantities                      ager 10/15    |
 *----------------------------------------------------------------------*/
void Adapter::CouplingPoroMortar::update_poro_mt()
{
  // safety check
  check_setup();

  porolagstrategy_->poro_mt_update();
  return;
}  // Adapter::CouplingNonLinMortar::UpdatePoroMt()

/*----------------------------------------------------------------------*
 |  recover fluid coupling lagrange multiplier            ager 10/15    |
 *----------------------------------------------------------------------*/
void Adapter::CouplingPoroMortar::recover_fluid_lm_poro_mt(
    Teuchos::RCP<Epetra_Vector> disi, Teuchos::RCP<Epetra_Vector> veli)
{
  // safety check
  check_setup();

  porolagstrategy_->recover_poro_no_pen(disi, veli);
  return;
}  // Adapter::CouplingNonLinMortar::recover_fluid_lm_poro_mt

FOUR_C_NAMESPACE_CLOSE
