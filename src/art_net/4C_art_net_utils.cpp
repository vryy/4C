/*----------------------------------------------------------------------*/
/*! \file

 \brief helper functions/classes for artery problems

   \level 3

 *----------------------------------------------------------------------*/

#include "4C_art_net_utils.hpp"

#include "4C_art_net_artery.hpp"
#include "4C_art_net_explicitintegration.hpp"
#include "4C_art_net_impl_stationary.hpp"
#include "4C_fem_general_utils_createdis.hpp"
#include "4C_global_data.hpp"
#include "4C_mat_par_bundle.hpp"
#include "4C_scatra_ele.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 | create algorithm                                                      |
 *----------------------------------------------------------------------*/
Teuchos::RCP<Adapter::ArtNet> Arteries::UTILS::CreateAlgorithm(
    Inpar::ArtDyn::TimeIntegrationScheme timintscheme, Teuchos::RCP<Core::FE::Discretization> dis,
    const int linsolvernumber, const Teuchos::ParameterList& probparams,
    const Teuchos::ParameterList& artparams, Teuchos::RCP<Core::IO::DiscretizationWriter> output)
{
  // Creation of Coupled Problem algortihm.
  Teuchos::RCP<Adapter::ArtNet> algo = Teuchos::null;

  // -------------------------------------------------------------------
  // algorithm construction depending on
  // time-integration (or stationary) scheme
  // -------------------------------------------------------------------

  switch (timintscheme)
  {
    case Inpar::ArtDyn::TimeIntegrationScheme::tay_gal:
    {
      // create algorithm
      algo = Teuchos::rcp(new Arteries::ArtNetExplicitTimeInt(
          dis, linsolvernumber, probparams, artparams, *output));
      break;
    }
    case Inpar::ArtDyn::TimeIntegrationScheme::stationary:
    {
      // create algorithm
      algo = Teuchos::rcp(
          new Arteries::ArtNetImplStationary(dis, linsolvernumber, probparams, artparams, *output));
      break;
    }
    default:
      FOUR_C_THROW("Unknown time-integration scheme for artery network problem");
      break;
  }

  return algo;
}

/*----------------------------------------------------------------------*
 | exchange material pointers of both discretizations  kremheller 03/18 |
 *----------------------------------------------------------------------*/
void Arteries::UTILS::assign_material_pointers(
    const std::string& artery_disname, const std::string& scatra_disname)
{
  Global::Problem* problem = Global::Problem::instance();

  Teuchos::RCP<Core::FE::Discretization> arterydis = problem->get_dis(artery_disname);
  Teuchos::RCP<Core::FE::Discretization> scatradis = problem->get_dis(scatra_disname);

  SetMaterialPointersMatchingGrid(arterydis, scatradis);
}

/*----------------------------------------------------------------------*
 | reset Material pointers after redistribution        kremheller 03/18 |
 *----------------------------------------------------------------------*/
void Arteries::UTILS::SetMaterialPointersMatchingGrid(
    Teuchos::RCP<const Core::FE::Discretization> sourcedis,
    Teuchos::RCP<const Core::FE::Discretization> targetdis)
{
  const int numelements = targetdis->num_my_col_elements();

  for (int i = 0; i < numelements; ++i)
  {
    Core::Elements::Element* targetele = targetdis->l_col_element(i);
    const int gid = targetele->id();

    Core::Elements::Element* sourceele = sourcedis->g_element(gid);

    // for coupling we add the source material to the target element and vice versa
    targetele->add_material(sourceele->material());
    sourceele->add_material(targetele->material());
  }
}

/*----------------------------------------------------------------------*
 |                                                     kremheller 03/18 |
 *----------------------------------------------------------------------*/
bool Arteries::ArteryScatraCloneStrategy::determine_ele_type(
    Core::Elements::Element* actele, const bool ismyele, std::vector<std::string>& eletype)
{
  // clone the element
  Discret::ELEMENTS::Artery* myele = static_cast<Discret::ELEMENTS::Artery*>(actele);
  // only the pressure based artery supports this function so far
  if (myele->impl_type() == Inpar::ArtDyn::impltype_pressure_based)
  {
    // we only support transport elements here
    eletype.push_back("TRANSP");
    return true;
  }

  return false;
}


/*----------------------------------------------------------------------*
 | set the element data (protected)                    kremheller 03/18 |
 *----------------------------------------------------------------------*/
void Arteries::ArteryScatraCloneStrategy::set_element_data(
    Teuchos::RCP<Core::Elements::Element> newele, Core::Elements::Element* oldele, const int matid,
    const bool isnurbs)
{
  // We need to set material and possibly other things to complete element setup.
  // This is again really ugly as we have to extract the actual
  // element type in order to access the material property

  // note: set_material() was reimplemented by the transport element!
  Discret::ELEMENTS::Transport* trans = dynamic_cast<Discret::ELEMENTS::Transport*>(newele.get());
  if (trans != nullptr)
  {
    // set material
    trans->set_material(matid, oldele);
    // set distype as well!
    trans->set_dis_type(oldele->shape());

    // we only have one possible impltype
    Inpar::ScaTra::ImplType impltype = Inpar::ScaTra::impltype_one_d_artery;
    trans->set_impl_type(impltype);
  }
  else
  {
    FOUR_C_THROW("unsupported element type '%s'", typeid(*newele).name());
  }
  return;
}

/*----------------------------------------------------------------------*
 | check if material type is admissible                kremheller 03/18 |
 *----------------------------------------------------------------------*/
void Arteries::ArteryScatraCloneStrategy::check_material_type(const int matid)
{
  // We take the material with the ID specified by the user
  // Here we check first, whether this material is of admissible type
  Core::Materials::MaterialType mtype =
      Global::Problem::instance()->materials()->parameter_by_id(matid)->type();
  if ((mtype != Core::Materials::m_scatra) && (mtype != Core::Materials::m_matlist))
    FOUR_C_THROW("Material with ID %d is not admissible for scalar transport elements", matid);
}


/*----------------------------------------------------------------------*
 |                                                     kremheller 03/18 |
 *----------------------------------------------------------------------*/
std::map<std::string, std::string> Arteries::ArteryScatraCloneStrategy::conditions_to_copy() const
{
  return {{"TransportDirichlet", "Dirichlet"}, {"TransportPointNeumann", "PointNeumann"},
      {"Initfield", "Initfield"}, {"ArtScatraCouplConNodebased", "ArtScatraCouplConNodebased"},
      {"PoroMultiphaseScatraOxyPartPressCalcCond", "PoroMultiphaseScatraOxyPartPressCalcCond"},
      {"ArtScatraCouplConNodeToPoint", "ArtScatraCouplConNodeToPoint"}};
}

FOUR_C_NAMESPACE_CLOSE
