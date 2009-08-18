
#ifdef CCADISCRET

#include "ale_utils_mapextractor.H"

#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_condition_selector.H"
#include "../drt_lib/drt_condition_utils.H"
#include "../drt_lib/standardtypes_cpp.H"

/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ALE::UTILS::MapExtractor::Setup(const DRT::Discretization& dis)
{
  DRT::UTILS::MultiConditionSelector mcs;
  mcs.AddSelector(rcp(new DRT::UTILS::NDimConditionSelector(dis,"FSICoupling",0,genprob.ndim)));
  mcs.AddSelector(rcp(new DRT::UTILS::NDimConditionSelector(dis,"FREESURFCoupling",0,genprob.ndim)));
  mcs.SetupExtractor(dis,*dis.DofRowMap(),*this);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<std::set<int> > ALE::UTILS::MapExtractor::ConditionedElementMap(const DRT::Discretization& dis) const
{
  Teuchos::RCP<std::set<int> > condelements = DRT::UTILS::ConditionedElementMap(dis,"FSICoupling");
  Teuchos::RCP<std::set<int> > condelements2 = DRT::UTILS::ConditionedElementMap(dis,"FREESURFCoupling");
  std::copy(condelements2->begin(),condelements2->end(),
            std::inserter(*condelements,condelements->begin()));
  return condelements;
}

#endif
