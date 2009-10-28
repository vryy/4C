#ifdef CCADISCRET

#include "structure_utils_mapextractor.H"

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
void STR::UTILS::MapExtractor::Setup(const DRT::Discretization& dis, Epetra_Map fullmap)
{
  DRT::UTILS::MultiConditionSelector mcs;
  mcs.AddSelector(rcp(new DRT::UTILS::NDimConditionSelector(dis,"FSICoupling",0,genprob.ndim)));
  mcs.AddSelector(rcp(new DRT::UTILS::NDimConditionSelector(dis,"StructAleCoupling",0,genprob.ndim)));
  mcs.SetupExtractor(dis,fullmap,*this);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<std::set<int> > STR::UTILS::MapExtractor::ConditionedElementMap(const DRT::Discretization& dis) const
{
  Teuchos::RCP<std::set<int> > condelements = DRT::UTILS::ConditionedElementMap(dis,"FSICoupling");
  Teuchos::RCP<std::set<int> > condelements2 = DRT::UTILS::ConditionedElementMap(dis,"StructAleCoupling");
  std::copy(condelements2->begin(),condelements2->end(),
            std::inserter(*condelements,condelements->begin()));
  return condelements;
}

#endif
