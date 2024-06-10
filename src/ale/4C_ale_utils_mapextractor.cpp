/*----------------------------------------------------------------------------*/
/*! \file


\brief map extractor functionalities for interface problems for ALE field

\level 1

*/
/*----------------------------------------------------------------------------*/
#include "4C_ale_utils_mapextractor.hpp"

#include "4C_fem_condition_selector.hpp"
#include "4C_fem_condition_utils.hpp"
#include "4C_fem_discretization.hpp"
#include "4C_global_data.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void ALE::UTILS::MapExtractor::Setup(const Core::FE::Discretization& dis, bool overlapping)
{
  const int ndim = Global::Problem::Instance()->NDim();
  Core::Conditions::MultiConditionSelector mcs;
  mcs.SetOverlapping(overlapping);
  mcs.AddSelector(
      Teuchos::rcp(new Core::Conditions::NDimConditionSelector(dis, "FSICoupling", 0, ndim)));
  mcs.AddSelector(
      Teuchos::rcp(new Core::Conditions::NDimConditionSelector(dis, "FREESURFCoupling", 0, ndim)));
  mcs.AddSelector(
      Teuchos::rcp(new Core::Conditions::NDimConditionSelector(dis, "StructAleCoupling", 0, ndim)));
  mcs.AddSelector(
      Teuchos::rcp(new Core::Conditions::NDimConditionSelector(dis, "AleWear", 0, ndim)));
  mcs.AddSelector(
      Teuchos::rcp(new Core::Conditions::NDimConditionSelector(dis, "BioGrCoupling", 0, ndim)));
  mcs.AddSelector(
      Teuchos::rcp(new Core::Conditions::NDimConditionSelector(dis, "ALEUPDATECoupling", 0, ndim)));
  mcs.AddSelector(
      Teuchos::rcp(new Core::Conditions::NDimConditionSelector(dis, "fpsi_coupling", 0, ndim)));
  mcs.AddSelector(
      Teuchos::rcp(new Core::Conditions::NDimConditionSelector(dis, "Mortar", 0, ndim)));
  mcs.SetupExtractor(dis, *dis.dof_row_map(), *this);
}


/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Teuchos::RCP<std::set<int>> ALE::UTILS::MapExtractor::conditioned_element_map(
    const Core::FE::Discretization& dis) const
{
  Teuchos::RCP<std::set<int>> condelements =
      Core::Conditions::conditioned_element_map(dis, "FSICoupling");
  Teuchos::RCP<std::set<int>> condelements2 =
      Core::Conditions::conditioned_element_map(dis, "FREESURFCoupling");
  Teuchos::RCP<std::set<int>> condelements3 =
      Core::Conditions::conditioned_element_map(dis, "StructAleCoupling");
  Teuchos::RCP<std::set<int>> condelements4 =
      Core::Conditions::conditioned_element_map(dis, "AleWear");
  Teuchos::RCP<std::set<int>> condelements5 =
      Core::Conditions::conditioned_element_map(dis, "BioGrCoupling");
  Teuchos::RCP<std::set<int>> condelements6 =
      Core::Conditions::conditioned_element_map(dis, "ALEUPDATECoupling");
  Teuchos::RCP<std::set<int>> condelements7 =
      Core::Conditions::conditioned_element_map(dis, "fpsi_coupling");
  Teuchos::RCP<std::set<int>> condelements8 =
      Core::Conditions::conditioned_element_map(dis, "Mortar");
  std::copy(condelements2->begin(), condelements2->end(),
      std::inserter(*condelements, condelements->begin()));
  std::copy(condelements3->begin(), condelements3->end(),
      std::inserter(*condelements, condelements->begin()));
  std::copy(condelements4->begin(), condelements4->end(),
      std::inserter(*condelements, condelements->begin()));
  std::copy(condelements5->begin(), condelements5->end(),
      std::inserter(*condelements, condelements->begin()));
  std::copy(condelements6->begin(), condelements6->end(),
      std::inserter(*condelements, condelements->begin()));
  std::copy(condelements7->begin(), condelements7->end(),
      std::inserter(*condelements, condelements->begin()));
  std::copy(condelements8->begin(), condelements8->end(),
      std::inserter(*condelements, condelements->begin()));
  return condelements;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ALE::UTILS::FsiMapExtractor::Setup(const Core::FE::Discretization& dis)
{
  const int ndim = Global::Problem::Instance()->NDim();
  Core::Conditions::MultiConditionSelector mcs;
  mcs.AddSelector(
      Teuchos::rcp(new Core::Conditions::NDimConditionSelector(dis, "FSICoupling", 0, ndim)));
  mcs.SetupExtractor(dis, *dis.dof_row_map(), *this);
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void ALE::UTILS::XFluidFluidMapExtractor::Setup(const Core::FE::Discretization& dis)
{
  const int ndim = Global::Problem::Instance()->NDim();
  Core::Conditions::MultiConditionSelector mcs;
  mcs.AddSelector(Teuchos::rcp(
      new Core::Conditions::NDimConditionSelector(dis, "FluidFluidCoupling", 0, ndim)));
  mcs.SetupExtractor(dis, *dis.dof_row_map(), *this);
}

FOUR_C_NAMESPACE_CLOSE
