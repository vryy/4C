/*-----------------------------------------------------------*/
/*!

\brief extracting maps of fluid discretizations

\maintainer Martin Kronbichler

\level 1

*/
/*-----------------------------------------------------------*/

#include "fluid_utils_mapextractor.H"
#include "fluid_utils.H"

#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_condition_selector.H"
#include "../drt_lib/drt_condition_utils.H"
#include "../drt_lib/drt_globalproblem.H"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FLD::UTILS::MapExtractor::Setup(
    const DRT::Discretization& dis, bool withpressure, bool overlapping, const int nds_master)
{
  const int ndim = DRT::Problem::Instance()->NDim();
  DRT::UTILS::MultiConditionSelector mcs;
  mcs.SetOverlapping(overlapping);  // defines if maps can overlap
  mcs.AddSelector(Teuchos::rcp(
      new DRT::UTILS::NDimConditionSelector(dis, "FSICoupling", 0, ndim + withpressure)));
  mcs.AddSelector(Teuchos::rcp(
      new DRT::UTILS::NDimConditionSelector(dis, "FREESURFCoupling", 0, ndim + withpressure)));
  mcs.AddSelector(Teuchos::rcp(
      new DRT::UTILS::NDimConditionSelector(dis, "StructAleCoupling", 0, ndim + withpressure)));
  mcs.AddSelector(
      Teuchos::rcp(new DRT::UTILS::NDimConditionSelector(dis, "Mortar", 0, ndim + withpressure)));
  mcs.AddSelector(Teuchos::rcp(
      new DRT::UTILS::NDimConditionSelector(dis, "ALEUPDATECoupling", 0, ndim + withpressure)));
  mcs.SetupExtractor(dis, *dis.DofRowMap(nds_master), *this);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FLD::UTILS::MapExtractor::Setup(
    Teuchos::RCP<const Epetra_Map>& additionalothermap, const FLD::UTILS::MapExtractor& extractor)
{
  // build the new othermap
  std::vector<Teuchos::RCP<const Epetra_Map>> othermaps;
  othermaps.push_back(additionalothermap);
  othermaps.push_back(extractor.OtherMap());

  if (LINALG::MultiMapExtractor::IntersectMaps(othermaps)->NumGlobalElements() > 0)
    dserror("Failed to add dofmap of foreign discretization to OtherMap. Detected overlap.");

  Teuchos::RCP<const Epetra_Map> mergedothermap = LINALG::MultiMapExtractor::MergeMaps(othermaps);

  // the vector of maps for the new map extractor consists of othermap at position 0
  // followed by the maps of conditioned DOF
  std::vector<Teuchos::RCP<const Epetra_Map>> maps;
  // append the merged other map at first position
  maps.push_back(mergedothermap);

  // append the condition maps subsequently
  for (int i = 1; i < extractor.NumMaps(); ++i) maps.push_back(extractor.Map(i));

  // merge
  Teuchos::RCP<const Epetra_Map> fullmap = LINALG::MultiMapExtractor::MergeMaps(maps);

  LINALG::MultiMapExtractor::Setup(*fullmap, maps);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<std::set<int>> FLD::UTILS::MapExtractor::ConditionedElementMap(
    const DRT::Discretization& dis) const
{
  Teuchos::RCP<std::set<int>> condelements = DRT::UTILS::ConditionedElementMap(dis, "FSICoupling");
  Teuchos::RCP<std::set<int>> condelements2 =
      DRT::UTILS::ConditionedElementMap(dis, "FREESURFCoupling");
  Teuchos::RCP<std::set<int>> condelements3 =
      DRT::UTILS::ConditionedElementMap(dis, "StructAleCoupling");
  Teuchos::RCP<std::set<int>> condelements4 = DRT::UTILS::ConditionedElementMap(dis, "Mortar");
  Teuchos::RCP<std::set<int>> condelements5 =
      DRT::UTILS::ConditionedElementMap(dis, "ALEUPDATECoupling");
  std::copy(condelements2->begin(), condelements2->end(),
      std::inserter(*condelements, condelements->begin()));
  std::copy(condelements3->begin(), condelements3->end(),
      std::inserter(*condelements, condelements->begin()));
  std::copy(condelements4->begin(), condelements4->end(),
      std::inserter(*condelements, condelements->begin()));
  std::copy(condelements5->begin(), condelements5->end(),
      std::inserter(*condelements, condelements->begin()));
  return condelements;
}

void FLD::UTILS::VolumetricFlowMapExtractor::Setup(const DRT::Discretization& dis)
{
  const int ndim = DRT::Problem::Instance()->NDim();
  DRT::UTILS::MultiConditionSelector mcs;
  mcs.SetOverlapping(true);  // defines if maps can overlap
  mcs.AddSelector(Teuchos::rcp(
      new DRT::UTILS::NDimConditionSelector(dis, "VolumetricSurfaceFlowCond", 0, ndim)));
  mcs.SetupExtractor(dis, *dis.DofRowMap(), *this);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FLD::UTILS::KSPMapExtractor::Setup(const DRT::Discretization& dis)
{
  DRT::UTILS::MultiConditionSelector mcs;
  mcs.AddSelector(Teuchos::rcp(new DRT::UTILS::ConditionSelector(dis, "KrylovSpaceProjection")));
  mcs.SetupExtractor(dis, *dis.DofRowMap(), *this);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<std::set<int>> FLD::UTILS::KSPMapExtractor::ConditionedElementMap(
    const DRT::Discretization& dis) const
{
  Teuchos::RCP<std::set<int>> condelements =
      DRT::UTILS::ConditionedElementMap(dis, "KrylovSpaceProjection");
  return condelements;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FLD::UTILS::VelPressExtractor::Setup(const DRT::Discretization& dis)
{
  const int ndim = DRT::Problem::Instance()->NDim();
  FLD::UTILS::SetupFluidSplit(dis, ndim, *this);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FLD::UTILS::FsiMapExtractor::Setup(const DRT::Discretization& dis)
{
  const int ndim = DRT::Problem::Instance()->NDim();
  DRT::UTILS::MultiConditionSelector mcs;
  mcs.AddSelector(Teuchos::rcp(new DRT::UTILS::NDimConditionSelector(dis, "FSICoupling", 0, ndim)));
  mcs.SetupExtractor(dis, *dis.DofRowMap(), *this);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FLD::UTILS::FsiMapExtractor::Setup(Teuchos::RCP<const Epetra_Map>& additionalothermap,
    const FLD::UTILS::FsiMapExtractor& extractor)
{
  // build the new othermap
  std::vector<Teuchos::RCP<const Epetra_Map>> othermaps;
  othermaps.push_back(additionalothermap);
  othermaps.push_back(extractor.OtherMap());

  if (LINALG::MultiMapExtractor::IntersectMaps(othermaps)->NumGlobalElements() > 0)
    dserror("Failed to add dofmap of foreign discretization to OtherMap. Detected overlap.");

  Teuchos::RCP<const Epetra_Map> mergedothermap = LINALG::MultiMapExtractor::MergeMaps(othermaps);

  // the vector of maps for the new map extractor consists of othermap at position 0
  // followed by the maps of conditioned DOF
  std::vector<Teuchos::RCP<const Epetra_Map>> maps;
  // append the merged other map at first position
  maps.push_back(mergedothermap);

  // append the condition maps subsequently
  for (int i = 1; i < extractor.NumMaps(); ++i) maps.push_back(extractor.Map(i));

  // merge
  Teuchos::RCP<const Epetra_Map> fullmap = LINALG::MultiMapExtractor::MergeMaps(maps);

  LINALG::MultiMapExtractor::Setup(*fullmap, maps);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FLD::UTILS::XFluidFluidMapExtractor::Setup(const Epetra_Map& fullmap,
    Teuchos::RCP<const Epetra_Map> fluidmap, Teuchos::RCP<const Epetra_Map> xfluidmap)
{
  std::vector<Teuchos::RCP<const Epetra_Map>> maps;
  maps.push_back(fluidmap);
  maps.push_back(xfluidmap);
  MultiMapExtractor::Setup(fullmap, maps);
}
