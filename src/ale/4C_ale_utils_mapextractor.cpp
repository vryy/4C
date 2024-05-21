/*----------------------------------------------------------------------------*/
/*! \file


\brief map extractor functionalities for interface problems for ALE field

\level 1

*/
/*----------------------------------------------------------------------------*/
#include "4C_ale_utils_mapextractor.hpp"

#include "4C_discretization_condition_selector.hpp"
#include "4C_discretization_condition_utils.hpp"
#include "4C_global_data.hpp"
#include "4C_lib_discret.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void ALE::UTILS::MapExtractor::Setup(const DRT::Discretization& dis, bool overlapping)
{
  const int ndim = GLOBAL::Problem::Instance()->NDim();
  CORE::Conditions::MultiConditionSelector mcs;
  mcs.SetOverlapping(overlapping);
  mcs.AddSelector(
      Teuchos::rcp(new CORE::Conditions::NDimConditionSelector(dis, "FSICoupling", 0, ndim)));
  mcs.AddSelector(
      Teuchos::rcp(new CORE::Conditions::NDimConditionSelector(dis, "FREESURFCoupling", 0, ndim)));
  mcs.AddSelector(
      Teuchos::rcp(new CORE::Conditions::NDimConditionSelector(dis, "StructAleCoupling", 0, ndim)));
  mcs.AddSelector(
      Teuchos::rcp(new CORE::Conditions::NDimConditionSelector(dis, "AleWear", 0, ndim)));
  mcs.AddSelector(
      Teuchos::rcp(new CORE::Conditions::NDimConditionSelector(dis, "BioGrCoupling", 0, ndim)));
  mcs.AddSelector(
      Teuchos::rcp(new CORE::Conditions::NDimConditionSelector(dis, "ALEUPDATECoupling", 0, ndim)));
  mcs.AddSelector(
      Teuchos::rcp(new CORE::Conditions::NDimConditionSelector(dis, "FPSICoupling", 0, ndim)));
  mcs.AddSelector(
      Teuchos::rcp(new CORE::Conditions::NDimConditionSelector(dis, "Mortar", 0, ndim)));
  mcs.SetupExtractor(dis, *dis.DofRowMap(), *this);
}


/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Teuchos::RCP<std::set<int>> ALE::UTILS::MapExtractor::ConditionedElementMap(
    const DRT::Discretization& dis) const
{
  Teuchos::RCP<std::set<int>> condelements =
      CORE::Conditions::ConditionedElementMap(dis, "FSICoupling");
  Teuchos::RCP<std::set<int>> condelements2 =
      CORE::Conditions::ConditionedElementMap(dis, "FREESURFCoupling");
  Teuchos::RCP<std::set<int>> condelements3 =
      CORE::Conditions::ConditionedElementMap(dis, "StructAleCoupling");
  Teuchos::RCP<std::set<int>> condelements4 =
      CORE::Conditions::ConditionedElementMap(dis, "AleWear");
  Teuchos::RCP<std::set<int>> condelements5 =
      CORE::Conditions::ConditionedElementMap(dis, "BioGrCoupling");
  Teuchos::RCP<std::set<int>> condelements6 =
      CORE::Conditions::ConditionedElementMap(dis, "ALEUPDATECoupling");
  Teuchos::RCP<std::set<int>> condelements7 =
      CORE::Conditions::ConditionedElementMap(dis, "FPSICoupling");
  Teuchos::RCP<std::set<int>> condelements8 =
      CORE::Conditions::ConditionedElementMap(dis, "Mortar");
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
void ALE::UTILS::FsiMapExtractor::Setup(const DRT::Discretization& dis)
{
  const int ndim = GLOBAL::Problem::Instance()->NDim();
  CORE::Conditions::MultiConditionSelector mcs;
  mcs.AddSelector(
      Teuchos::rcp(new CORE::Conditions::NDimConditionSelector(dis, "FSICoupling", 0, ndim)));
  mcs.SetupExtractor(dis, *dis.DofRowMap(), *this);
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void ALE::UTILS::XFluidFluidMapExtractor::Setup(const DRT::Discretization& dis)
{
  const int ndim = GLOBAL::Problem::Instance()->NDim();
  CORE::Conditions::MultiConditionSelector mcs;
  mcs.AddSelector(Teuchos::rcp(
      new CORE::Conditions::NDimConditionSelector(dis, "FluidFluidCoupling", 0, ndim)));
  mcs.SetupExtractor(dis, *dis.DofRowMap(), *this);
}

FOUR_C_NAMESPACE_CLOSE
