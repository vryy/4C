/*----------------------------------------------------------------------*/
/*! \file

 \brief helper functions/classes for artery problems

   \level 3

 *----------------------------------------------------------------------*/

#include "baci_art_net_utils.hpp"

#include "baci_art_net_artery.hpp"
#include "baci_art_net_explicitintegration.hpp"
#include "baci_art_net_impl_stationary.hpp"
#include "baci_lib_utils_createdis.hpp"
#include "baci_mat_par_bundle.hpp"
#include "baci_scatra_ele.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 | create algorithm                                                      |
 *----------------------------------------------------------------------*/
Teuchos::RCP<ADAPTER::ArtNet> ART::UTILS::CreateAlgorithm(
    INPAR::ARTDYN::TimeIntegrationScheme timintscheme, Teuchos::RCP<DRT::Discretization> dis,
    const int linsolvernumber, const Teuchos::ParameterList& probparams,
    const Teuchos::ParameterList& artparams, Teuchos::RCP<IO::DiscretizationWriter> output)
{
  // Creation of Coupled Problem algortihm.
  Teuchos::RCP<ADAPTER::ArtNet> algo = Teuchos::null;

  // -------------------------------------------------------------------
  // algorithm construction depending on
  // time-integration (or stationary) scheme
  // -------------------------------------------------------------------

  switch (timintscheme)
  {
    case INPAR::ARTDYN::TimeIntegrationScheme::tay_gal:
    {
      // create algorithm
      algo = Teuchos::rcp(
          new ART::ArtNetExplicitTimeInt(dis, linsolvernumber, probparams, artparams, *output));
      break;
    }
    case INPAR::ARTDYN::TimeIntegrationScheme::stationary:
    {
      // create algorithm
      algo = Teuchos::rcp(
          new ART::ArtNetImplStationary(dis, linsolvernumber, probparams, artparams, *output));
      break;
    }
    default:
      dserror("Unknown time-integration scheme for artery network problem");
      break;
  }

  return algo;
}

/*----------------------------------------------------------------------*
 | exchange material pointers of both discretizations  kremheller 03/18 |
 *----------------------------------------------------------------------*/
void ART::UTILS::AssignMaterialPointers(
    const std::string& artery_disname, const std::string& scatra_disname)
{
  GLOBAL::Problem* problem = GLOBAL::Problem::Instance();

  Teuchos::RCP<DRT::Discretization> arterydis = problem->GetDis(artery_disname);
  Teuchos::RCP<DRT::Discretization> scatradis = problem->GetDis(scatra_disname);

  SetMaterialPointersMatchingGrid(arterydis, scatradis);
}

/*----------------------------------------------------------------------*
 | reset Material pointers after redistribution        kremheller 03/18 |
 *----------------------------------------------------------------------*/
void ART::UTILS::SetMaterialPointersMatchingGrid(Teuchos::RCP<const DRT::Discretization> sourcedis,
    Teuchos::RCP<const DRT::Discretization> targetdis)
{
  const int numelements = targetdis->NumMyColElements();

  for (int i = 0; i < numelements; ++i)
  {
    DRT::Element* targetele = targetdis->lColElement(i);
    const int gid = targetele->Id();

    DRT::Element* sourceele = sourcedis->gElement(gid);

    // for coupling we add the source material to the target element and vice versa
    targetele->AddMaterial(sourceele->Material());
    sourceele->AddMaterial(targetele->Material());
  }
}

/*----------------------------------------------------------------------*
 |                                                     kremheller 03/18 |
 *----------------------------------------------------------------------*/
bool ART::ArteryScatraCloneStrategy::DetermineEleType(
    DRT::Element* actele, const bool ismyele, std::vector<std::string>& eletype)
{
  // clone the element
  DRT::ELEMENTS::Artery* myele = static_cast<DRT::ELEMENTS::Artery*>(actele);
  // only the pressure based artery supports this function so far
  if (myele->ImplType() == INPAR::ARTDYN::impltype_pressure_based)
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
void ART::ArteryScatraCloneStrategy::SetElementData(
    Teuchos::RCP<DRT::Element> newele, DRT::Element* oldele, const int matid, const bool isnurbs)
{
  // We need to set material and possibly other things to complete element setup.
  // This is again really ugly as we have to extract the actual
  // element type in order to access the material property

  // note: SetMaterial() was reimplemented by the transport element!
  DRT::ELEMENTS::Transport* trans = dynamic_cast<DRT::ELEMENTS::Transport*>(newele.get());
  if (trans != nullptr)
  {
    // set material
    trans->SetMaterial(matid, oldele);
    // set distype as well!
    trans->SetDisType(oldele->Shape());

    // we only have one possible impltype
    INPAR::SCATRA::ImplType impltype = INPAR::SCATRA::impltype_one_d_artery;
    trans->SetImplType(impltype);
  }
  else
  {
    dserror("unsupported element type '%s'", typeid(*newele).name());
  }
  return;
}

/*----------------------------------------------------------------------*
 | check if material type is admissible                kremheller 03/18 |
 *----------------------------------------------------------------------*/
void ART::ArteryScatraCloneStrategy::CheckMaterialType(const int matid)
{
  // We take the material with the ID specified by the user
  // Here we check first, whether this material is of admissible type
  INPAR::MAT::MaterialType mtype = GLOBAL::Problem::Instance()->Materials()->ById(matid)->Type();
  if ((mtype != INPAR::MAT::m_scatra) && (mtype != INPAR::MAT::m_matlist))
    dserror("Material with ID %d is not admissible for scalar transport elements", matid);
}


/*----------------------------------------------------------------------*
 |                                                     kremheller 03/18 |
 *----------------------------------------------------------------------*/
std::map<std::string, std::string> ART::ArteryScatraCloneStrategy::ConditionsToCopy() const
{
  return {{"TransportDirichlet", "Dirichlet"}, {"TransportPointNeumann", "PointNeumann"},
      {"Initfield", "Initfield"}, {"ArtScatraCouplConNodebased", "ArtScatraCouplConNodebased"},
      {"PoroMultiphaseScatraOxyPartPressCalcCond", "PoroMultiphaseScatraOxyPartPressCalcCond"},
      {"ArtScatraCouplConNodeToPoint", "ArtScatraCouplConNodeToPoint"}};
}

FOUR_C_NAMESPACE_CLOSE
