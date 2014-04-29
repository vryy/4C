


#include "fluid_utils_mapextractor.H"
#include "fluid_utils.H"

#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_condition_selector.H"
#include "../drt_lib/drt_condition_utils.H"
#include "../drt_lib/drt_globalproblem.H"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FLD::UTILS::MapExtractor::Setup(const DRT::Discretization& dis, bool withpressure, bool overlapping)
{
  const int ndim = DRT::Problem::Instance()->NDim();
  DRT::UTILS::MultiConditionSelector mcs;
  mcs.SetOverlapping(overlapping); //defines if maps can overlap
  mcs.AddSelector(Teuchos::rcp(new DRT::UTILS::NDimConditionSelector(dis,"FSICoupling",0,ndim+withpressure)));
  mcs.AddSelector(Teuchos::rcp(new DRT::UTILS::NDimConditionSelector(dis,"FREESURFCoupling",0,ndim+withpressure)));
  mcs.AddSelector(Teuchos::rcp(new DRT::UTILS::NDimConditionSelector(dis,"StructAleCoupling",0,ndim+withpressure)));
  mcs.AddSelector(Teuchos::rcp(new DRT::UTILS::NDimConditionSelector(dis,"VolumetricSurfaceFlowCond",0,ndim+withpressure)));
  mcs.AddSelector(Teuchos::rcp(new DRT::UTILS::NDimConditionSelector(dis,"Mortar",0,ndim+withpressure)));
  mcs.AddSelector(Teuchos::rcp(new DRT::UTILS::NDimConditionSelector(dis,"TotalTractionCorrectionCond",0,ndim+withpressure)));
  mcs.AddSelector(Teuchos::rcp(new DRT::UTILS::NDimConditionSelector(dis,"LOCALLAGRANGECoupling",0,ndim+withpressure)));
  mcs.AddSelector(Teuchos::rcp(new DRT::UTILS::NDimConditionSelector(dis,"FPSICoupling",0,ndim+withpressure)));
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
  Teuchos::RCP<std::set<int> > condelements5 = DRT::UTILS::ConditionedElementMap(dis,"Mortar");
  Teuchos::RCP<std::set<int> > condelements6 = DRT::UTILS::ConditionedElementMap(dis,"TotalTractionCorrectionCond");
  Teuchos::RCP<std::set<int> > condelements7 = DRT::UTILS::ConditionedElementMap(dis,"LOCALLAGRANGECoupling");
  Teuchos::RCP<std::set<int> > condelements8 = DRT::UTILS::ConditionedElementMap(dis,"FPSICoupling");
  std::copy(condelements2->begin(),condelements2->end(),
            std::inserter(*condelements,condelements->begin()));
  std::copy(condelements3->begin(),condelements3->end(),
            std::inserter(*condelements,condelements->begin()));
  std::copy(condelements4->begin(),condelements4->end(),
            std::inserter(*condelements,condelements->begin()));
  std::copy(condelements5->begin(),condelements5->end(),
            std::inserter(*condelements,condelements->begin()));
  std::copy(condelements7->begin(),condelements7->end(),
              std::inserter(*condelements,condelements->begin()));
  std::copy(condelements8->begin(),condelements8->end(),
            std::inserter(*condelements,condelements->begin()));
  return condelements;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FLD::UTILS::KSPMapExtractor::Setup(const DRT::Discretization& dis)
{
  DRT::UTILS::MultiConditionSelector mcs;
  mcs.AddSelector(Teuchos::rcp(new DRT::UTILS::ConditionSelector(dis,"KrylovSpaceProjection")));
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
  const int ndim = DRT::Problem::Instance()->NDim();
  FLD::UTILS::SetupFluidSplit( dis, ndim, *this );
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FLD::UTILS::FluidXFluidMapExtractor::Setup(const Epetra_Map& fullmap, Teuchos::RCP<const Epetra_Map> fluidmap, Teuchos::RCP<const Epetra_Map> xfluidmap)
{
  std::vector<Teuchos::RCP<const Epetra_Map> > maps;
  maps.push_back(fluidmap);
  maps.push_back(xfluidmap);
  MultiMapExtractor::Setup(fullmap,maps);
}

