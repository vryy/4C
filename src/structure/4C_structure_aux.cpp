/*----------------------------------------------------------------------*/
/*! \file
\brief structure-specific utils and auxiliary functions

\level 1

*/
/*----------------------------------------------------------------------*/



#include "4C_structure_aux.hpp"

#include "4C_fem_condition_selector.hpp"
#include "4C_fem_condition_utils.hpp"
#include "4C_fem_discretization.hpp"
#include "4C_global_data.hpp"

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*/
/* Calculate vector norm */
double STR::calculate_vector_norm(const enum Inpar::STR::VectorNorm norm,
    const Teuchos::RCP<Epetra_Vector> vect, const int numneglect)
{
  // L1 norm
  if (norm == Inpar::STR::norm_l1)
  {
    double vectnorm;
    vect->Norm1(&vectnorm);
    return vectnorm;
  }
  // L2/Euclidian norm
  else if (norm == Inpar::STR::norm_l2)
  {
    double vectnorm;
    vect->Norm2(&vectnorm);
    return vectnorm;
  }
  // RMS norm
  else if (norm == Inpar::STR::norm_rms)
  {
    double vectnorm;
    vect->Norm2(&vectnorm);
    return vectnorm / sqrt((double)(vect->GlobalLength() - numneglect));
  }
  // infinity/maximum norm
  else if (norm == Inpar::STR::norm_inf)
  {
    double vectnorm;
    vect->NormInf(&vectnorm);
    return vectnorm;
  }
  else
  {
    FOUR_C_THROW("Cannot handle vector norm");
    return 0;
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void STR::MapExtractor::Setup(
    const Discret::Discretization& dis, const Epetra_Map& fullmap, bool overlapping)
{
  const int ndim = Global::Problem::Instance()->NDim();
  Core::Conditions::MultiConditionSelector mcs;
  mcs.SetOverlapping(overlapping);
  mcs.AddSelector(
      Teuchos::rcp(new Core::Conditions::NDimConditionSelector(dis, "FSICoupling", 0, ndim)));
  mcs.AddSelector(
      Teuchos::rcp(new Core::Conditions::NDimConditionSelector(dis, "StructAleCoupling", 0, ndim)));
  mcs.AddSelector(
      Teuchos::rcp(new Core::Conditions::NDimConditionSelector(dis, "BioGrCoupling", 0, ndim)));
  mcs.AddSelector(
      Teuchos::rcp(new Core::Conditions::NDimConditionSelector(dis, "AleWear", 0, ndim)));
  mcs.AddSelector(
      Teuchos::rcp(new Core::Conditions::NDimConditionSelector(dis, "fpsi_coupling", 0, ndim)));
  mcs.AddSelector(
      Teuchos::rcp(new Core::Conditions::NDimConditionSelector(dis, "IMMERSEDCoupling", 0, ndim)));
  mcs.AddSelector(
      Teuchos::rcp(new Core::Conditions::NDimConditionSelector(dis, "ParticleWall", 0, ndim)));

  mcs.SetupExtractor(dis, fullmap, *this);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<std::set<int>> STR::MapExtractor::conditioned_element_map(
    const Discret::Discretization& dis) const
{
  Teuchos::RCP<std::set<int>> condelements =
      Core::Conditions::conditioned_element_map(dis, "FSICoupling");
  Teuchos::RCP<std::set<int>> condelements2 =
      Core::Conditions::conditioned_element_map(dis, "StructAleCoupling");
  Teuchos::RCP<std::set<int>> condelements3 =
      Core::Conditions::conditioned_element_map(dis, "BioGrCoupling");
  Teuchos::RCP<std::set<int>> condelements4 =
      Core::Conditions::conditioned_element_map(dis, "AleWear");
  Teuchos::RCP<std::set<int>> condelements5 =
      Core::Conditions::conditioned_element_map(dis, "fpsi_coupling");
  Teuchos::RCP<std::set<int>> condelements6 =
      Core::Conditions::conditioned_element_map(dis, "IMMERSEDCoupling");
  Teuchos::RCP<std::set<int>> condelements7 =
      Core::Conditions::conditioned_element_map(dis, "ParticleWall");


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
  return condelements;
}

FOUR_C_NAMESPACE_CLOSE
