
#ifdef CCADISCRET

#include "fluid_utils_mapextractor.H"
#include "fluid_utils.H"

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
void FLD::UTILS::MapExtractor::Setup(const DRT::Discretization& dis)
{
  DRT::UTILS::MultiConditionSelector mcs;
  mcs.AddSelector(rcp(new DRT::UTILS::NDimConditionSelector(dis,"FSICoupling",0,genprob.ndim)));
  mcs.AddSelector(rcp(new DRT::UTILS::NDimConditionSelector(dis,"FREESURFCoupling",0,genprob.ndim)));
  mcs.AddSelector(rcp(new DRT::UTILS::NDimConditionSelector(dis,"StructAleCoupling",0,genprob.ndim)));
  mcs.AddSelector(rcp(new DRT::UTILS::NDimConditionSelector(dis,"VolumetricSurfaceFlowCond",0,genprob.ndim)));
  mcs.SetupExtractor(dis,*dis.DofRowMap(),*this);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<std::set<int> > FLD::UTILS::MapExtractor::ConditionedElementMap(const DRT::Discretization& dis) const
{
  Teuchos::RCP<std::set<int> > condelements  = DRT::UTILS::ConditionedElementMap(dis,"FSICoupling");
  Teuchos::RCP<std::set<int> > condelements2 = DRT::UTILS::ConditionedElementMap(dis,"FREESURFCoupling");
  Teuchos::RCP<std::set<int> > condelements3 = DRT::UTILS::ConditionedElementMap(dis,"StructAleCoupling");
  Teuchos::RCP<std::set<int> > condelements4 = DRT::UTILS::ConditionedElementMap(dis,"VolumetricSurfaceFlowCond");
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
void FLD::UTILS::KSPMapExtractor::Setup(const DRT::Discretization& dis)
{
  DRT::UTILS::MultiConditionSelector mcs;
  mcs.AddSelector(rcp(new DRT::UTILS::NDimConditionSelector(dis,"KrylovSpaceProjection",0,genprob.ndim+1)));
  mcs.SetupExtractor(dis,*dis.DofRowMap(),*this);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<std::set<int> > FLD::UTILS::KSPMapExtractor::ConditionedElementMap(const DRT::Discretization& dis) const
{
  Teuchos::RCP<std::set<int> > condelements = DRT::UTILS::ConditionedElementMap(dis,"KrylovSpaceProjection");
  return condelements;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FLD::UTILS::VelPressExtractor::Setup(const DRT::Discretization& dis)
{
  FLD::UTILS::SetupFluidSplit( dis, genprob.ndim, *this );
}

#endif

