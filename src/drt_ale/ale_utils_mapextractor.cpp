


#include "ale_utils_mapextractor.H"

#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_condition_selector.H"
#include "../drt_lib/drt_condition_utils.H"
#include "../drt_lib/drt_globalproblem.H"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ALE::UTILS::MapExtractor::Setup(const DRT::Discretization& dis)
{
  const int ndim = DRT::Problem::Instance()->NDim();
  DRT::UTILS::MultiConditionSelector mcs;
  mcs.AddSelector(Teuchos::rcp(new DRT::UTILS::NDimConditionSelector(dis,"FSICoupling",0,ndim)));
  mcs.AddSelector(Teuchos::rcp(new DRT::UTILS::NDimConditionSelector(dis,"FREESURFCoupling",0,ndim)));
  mcs.AddSelector(Teuchos::rcp(new DRT::UTILS::NDimConditionSelector(dis,"StructAleCoupling",0,ndim)));
  mcs.AddSelector(Teuchos::rcp(new DRT::UTILS::NDimConditionSelector(dis,"AleWear",0,ndim)));
  mcs.SetupExtractor(dis,*dis.DofRowMap(),*this);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<std::set<int> > ALE::UTILS::MapExtractor::ConditionedElementMap(const DRT::Discretization& dis) const
{
  Teuchos::RCP<std::set<int> > condelements = DRT::UTILS::ConditionedElementMap(dis,"FSICoupling");
  Teuchos::RCP<std::set<int> > condelements2 = DRT::UTILS::ConditionedElementMap(dis,"FREESURFCoupling");
  Teuchos::RCP<std::set<int> > condelements3 = DRT::UTILS::ConditionedElementMap(dis,"StructAleCoupling");
  Teuchos::RCP<std::set<int> > condelements4 = DRT::UTILS::ConditionedElementMap(dis,"AleWear");
  std::copy(condelements2->begin(),condelements2->end(),
            std::inserter(*condelements,condelements->begin()));
  std::copy(condelements3->begin(),condelements3->end(),
            std::inserter(*condelements,condelements->begin()));
  std::copy(condelements4->begin(),condelements4->end(),
            std::inserter(*condelements,condelements->begin()));
  return condelements;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ALE::UTILS::XFluidFluidMapExtractor::Setup(const DRT::Discretization& dis)
{
  const int ndim = DRT::Problem::Instance()->NDim();
  DRT::UTILS::MultiConditionSelector mcs;
  mcs.AddSelector(Teuchos::rcp(new DRT::UTILS::NDimConditionSelector(dis,"FluidFluidCoupling",0,ndim)));
  mcs.SetupExtractor(dis,*dis.DofRowMap(),*this);
}


