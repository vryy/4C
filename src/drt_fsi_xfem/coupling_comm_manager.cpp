/*----------------------------------------------------------------------*/
/*! \file
\brief  Coupling Communication Manager automatically creates all required coupling object to
transform matrixes, vectors, ...

\level 2

\maintainer Martin Kronbichler
            kronbichler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289 15249

*----------------------------------------------------------------------*/

#include "coupling_comm_manager.H"

#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_condition_selector.H"
#include "../drt_lib/drt_colors.H"

#include "../drt_fsi/fsi_matrixtransform.H"

#include "../linalg/linalg_mapextractor.H"
#include "../linalg/linalg_sparsematrix.H"
#include "../linalg/linalg_utils.H"

#include "../drt_adapter/adapter_coupling.H"

#include <iostream>
#include <string>


/*-----------------------------------------------------------------------------------------*
| constructor - simplified for 1 discretization to couple (public)             ager 03/2016|
*-----------------------------------------------------------------------------------------*/
XFEM::Coupling_Comm_Manager::Coupling_Comm_Manager(
    Teuchos::RCP<const DRT::Discretization> dis0, std::string cond_name, int startdim, int enddim)
    : cond_name_(cond_name), startdim_(startdim), enddim_(enddim)
{
  std::map<int, Teuchos::RCP<const DRT::Discretization>> dis_map;
  dis_map.insert(std::pair<int, Teuchos::RCP<const DRT::Discretization>>(0, dis0));

  Setup(dis_map);
}

/*-----------------------------------------------------------------------------------------*
| constructor - simplified for 2 discretizations to couple (public)            ager 03/2016|
*-----------------------------------------------------------------------------------------*/
XFEM::Coupling_Comm_Manager::Coupling_Comm_Manager(Teuchos::RCP<const DRT::Discretization> dis0,
    Teuchos::RCP<const DRT::Discretization> dis1, std::string cond_name, int startdim, int enddim)
    : cond_name_(cond_name), startdim_(startdim), enddim_(enddim)
{
  std::map<int, Teuchos::RCP<const DRT::Discretization>> dis_map;
  dis_map.insert(std::pair<int, Teuchos::RCP<const DRT::Discretization>>(0, dis0));
  dis_map.insert(std::pair<int, Teuchos::RCP<const DRT::Discretization>>(1, dis1));

  Setup(dis_map);
}

/*-----------------------------------------------------------------------------------------*
| constructor (public)                                                         ager 03/2016|
*-----------------------------------------------------------------------------------------*/
XFEM::Coupling_Comm_Manager::Coupling_Comm_Manager(
    std::map<int, Teuchos::RCP<const DRT::Discretization>> dis, std::string cond_name, int startdim,
    int enddim)
    : cond_name_(cond_name), startdim_(startdim), enddim_(enddim)
{
  Setup(dis);
}

/*-----------------------------------------------------------------------------------------*
| dstructor (public)                                                         ager 03/2016|
*-----------------------------------------------------------------------------------------*/
XFEM::Coupling_Comm_Manager::~Coupling_Comm_Manager() {}

/*--------------------------------------------------------------------------------------------------------------*
| Transfer conditioned part of Vector from Discretization A --> B with different transfer types ager
03/2016|
*--------------------------------------------------------------------------------------------------------------*/
void XFEM::Coupling_Comm_Manager::InsertVector(const int idxA,
    Teuchos::RCP<const Epetra_Vector> vecA, const int idxB, Teuchos::RCP<Epetra_Vector> vecB,
    const Coupling_Comm_Manager::transfer_type ttype, bool add, double scale)
{
  switch (ttype)
  {
    case Coupling_Comm_Manager::full_to_full:
    {
      Teuchos::RCP<LINALG::MultiMapExtractor> mmeb = GetMapExtractor(idxB);
      Teuchos::RCP<Epetra_Vector> tmpvec = Teuchos::rcp(new Epetra_Vector(*mmeb->Map(1), true));
      InsertVector(idxA, vecA, idxB, tmpvec, Coupling_Comm_Manager::full_to_partial, false, scale);
      if (!add)
        mmeb->InsertVector(tmpvec, 1, vecB);
      else
        mmeb->AddVector(tmpvec, 1, vecB);
      break;
    }
    case Coupling_Comm_Manager::full_to_partial:
    {
      InsertVector(idxA, GetMapExtractor(idxA)->ExtractVector(vecA, 1), idxB, vecB,
          Coupling_Comm_Manager::partial_to_partial, add, scale);
      break;
    }
    case Coupling_Comm_Manager::partial_to_full:
    {
      Teuchos::RCP<LINALG::MultiMapExtractor> mmeb = GetMapExtractor(idxB);
      Teuchos::RCP<Epetra_Vector> tmpvec = Teuchos::rcp(new Epetra_Vector(*mmeb->Map(1), true));
      InsertVector(
          idxA, vecA, idxB, tmpvec, Coupling_Comm_Manager::partial_to_partial, false, scale);
      if (!add)
        mmeb->InsertVector(tmpvec, 1, vecB);
      else
        mmeb->AddVector(tmpvec, 1, vecB);
      break;
    }
    case Coupling_Comm_Manager::partial_to_partial:
    {
      if (vecB == Teuchos::null)
        dserror("Coupling_Comm_Manager::InsertVector: vecB is Teuchos::null!");
      if (idxA < idxB)  // this Coupling Object is directly stored
      {
        if (!add)
          *vecB = *GetCoupling(idxA, idxB)->MasterToSlave(vecA);
        else
          vecB->Update(scale, *GetCoupling(idxA, idxB)->MasterToSlave(vecA), 1.0);
      }
      else if (idxA > idxB)  // just the inverse Coupling Object is stored
      {
        if (!add)
          *vecB = *GetCoupling(idxB, idxA)->SlaveToMaster(vecA);
        else
          vecB->Update(scale, *GetCoupling(idxB, idxA)->SlaveToMaster(vecA), 1.0);
      }
      else
      {
        if (!add)
        {
          // *vecB = *vecA; //we don't do any transformation!
          vecB->Update(1.0, *vecA, 0.0);
        }
        else
          vecB->Update(scale, *vecA, 1.0);
      }
      if (!add && scale != 1.0) vecB->Scale(scale);
      break;
    }
    case Coupling_Comm_Manager::partial_to_global:
    {
      Teuchos::RCP<LINALG::MultiMapExtractor> mme = GetFullMapExtractor();
      Teuchos::RCP<Epetra_Vector> fullvec = Teuchos::rcp(new Epetra_Vector(*mme->Map(idxB), true));
      InsertVector(idxA, vecA, idxB, fullvec, Coupling_Comm_Manager::partial_to_full, false, scale);
      if (!add)
        mme->InsertVector(fullvec, idxB, vecB);
      else
        mme->AddVector(fullvec, idxB, vecB);
      break;
    }
    case Coupling_Comm_Manager::full_to_global:
    {
      Teuchos::RCP<LINALG::MultiMapExtractor> mme = GetFullMapExtractor();
      Teuchos::RCP<Epetra_Vector> fullvec = Teuchos::rcp(new Epetra_Vector(*mme->Map(idxB), true));
      InsertVector(idxA, vecA, idxB, fullvec, Coupling_Comm_Manager::full_to_full, false, scale);
      if (!add)
        mme->InsertVector(fullvec, idxB, vecB);
      else
        mme->AddVector(fullvec, idxB, vecB);
      break;
    }
    default:
      dserror("Coupling_Comm_Manager::InsertVector: Transfer Type not implemented!");
  }
  return;
}

/*----------------------------------------------------------------------------------------------*
| //! Insert a Matrix A into Matrix B (choose type of transfer, add or scaling)     ager 03/2016|
*----------------------------------------------------------------------------------------------*/
bool XFEM::Coupling_Comm_Manager::InsertMatrix(int transform_id, int idxA,
    const LINALG::SparseMatrix& matA, int idxB, LINALG::SparseMatrix& matB,
    const Coupling_Comm_Manager::matrix_transfer_type mttype, double scale, bool exactmatch,
    bool addmatrix)
{
  switch (mttype)
  {
    case Coupling_Comm_Manager::col:
    {
      return GetTransform(transform_id)
          ->
          operator()(matA, matA.RangeMap(), matA.DomainMap(), scale, NULL,
              GetCouplingConverter(idxA, idxB).getRawPtr(), matB, exactmatch, addmatrix);
      break;
    }
    case Coupling_Comm_Manager::row:
    {
      return GetTransform(transform_id)
          ->
          operator()(matA, matA.RangeMap(), matA.DomainMap(), scale,
              GetCouplingConverter(idxA, idxB).getRawPtr(), NULL, matB, true, addmatrix);
      break;
    }
    case Coupling_Comm_Manager::row_and_col:
    {
      return GetTransform(transform_id)
          ->
          operator()(matA, matA.RangeMap(), matA.DomainMap(), scale,
              GetCouplingConverter(idxA, idxB).getRawPtr(),
              GetCouplingConverter(idxA, idxB).getRawPtr(), matB, exactmatch, addmatrix);
      break;
    }
    default:
      dserror("Coupling_Comm_Manager::InsertMatrix: Matrix Transfer Type not implemented!");
  }
  return true;
}

/*-----------------------------------------------------------------------------------------*
| Setup Coupling_Comm_Manager                                                  ager 06/2016|
*-----------------------------------------------------------------------------------------*/
void XFEM::Coupling_Comm_Manager::Setup(std::map<int, Teuchos::RCP<const DRT::Discretization>> dis)
{
  if (cond_name_ != "")  // Setup for Communication on Condition
  {
    SetupMultiMapExtractors(dis);
    SetupCouplings(dis);
  }
  else  // Setup for Communication on full Discretization
  {
    SetupFullCouplings(dis);      // First Couple full discretizations (from startdim to enddim)
    SetupFullMapExtractors(dis);  // Setup Extractors between Couplingpart of discretization and
                                  // full discretization (e.g. pres <==> vel&pres)
  }
  SetupFullExtractor(dis);
  return;
}

/*-----------------------------------------------------------------------------------------*
| Setup MultiMapExtractors for all fields                                      ager 03/2016|
*-----------------------------------------------------------------------------------------*/
void XFEM::Coupling_Comm_Manager::SetupMultiMapExtractors(
    std::map<int, Teuchos::RCP<const DRT::Discretization>> dis)
{
  for (std::map<int, Teuchos::RCP<const DRT::Discretization>>::iterator dit = dis.begin();
       dit != dis.end(); ++dit)
  {
    Teuchos::RCP<DRT::UTILS::MultiConditionSelector> mcs =
        Teuchos::rcp(new DRT::UTILS::MultiConditionSelector());
    mme_[dit->first] = Teuchos::rcp(new LINALG::MultiMapExtractor());
    mcs->AddSelector(Teuchos::rcp(
        new DRT::UTILS::NDimConditionSelector(*dit->second, cond_name_, startdim_, enddim_)));
    mcs->SetupExtractor(*dit->second, *dit->second->DofRowMap(), *mme_[dit->first]);
  }
}

/*-----------------------------------------------------------------------------------------*
| Setup MultiMapExtractors for all fields                                      ager 03/2016|
*-----------------------------------------------------------------------------------------*/
void XFEM::Coupling_Comm_Manager::SetupFullMapExtractors(
    std::map<int, Teuchos::RCP<const DRT::Discretization>> dis)
{
  if (dis.size() < 2)
    dserror(
        "SetupFullMapExtractors: Just extract from discretization to coupled map! (e.g. Fluid vel "
        "<==> Fluid vel&pres)");

  for (std::map<int, Teuchos::RCP<const DRT::Discretization>>::iterator dit = dis.begin();
       dit != dis.end(); ++dit)
  {
    Teuchos::RCP<LINALG::MapExtractor> me = Teuchos::rcp(new LINALG::MapExtractor());
    if ((uint)dit->first < dis.size() - 1)
    {
      Teuchos::RCP<ADAPTER::Coupling> coup = GetCoupling(dit->first, dit->first + 1);
      me->Setup(*dit->second->DofRowMap(), coup->MasterDofMap(),
          LINALG::SplitMap(*dit->second->DofRowMap(), *coup->MasterDofMap()));
    }
    else
    {
      Teuchos::RCP<ADAPTER::Coupling> coup = GetCoupling(dit->first - 1, dit->first);
      me->Setup(*dit->second->DofRowMap(), coup->SlaveDofMap(),
          LINALG::SplitMap(*dit->second->DofRowMap(), *coup->SlaveDofMap()));
    }
    mme_[dit->first] = me;
  }
  return;
}

/*------------------------------------------------------------------------------------------------*
| Setup Couplings between different discretizations with the same condition on in  -- ager 03/2016|
*------------------------------------------------------------------------------------------------*/
void XFEM::Coupling_Comm_Manager::SetupCouplings(
    std::map<int, Teuchos::RCP<const DRT::Discretization>> dis)
{
  if (dis.size() < 2) return;

  for (std::map<int, Teuchos::RCP<LINALG::MultiMapExtractor>>::iterator mmealpha = mme_.begin();
       mmealpha != mme_.end(); ++mmealpha)
  {
    for (std::map<int, Teuchos::RCP<LINALG::MultiMapExtractor>>::iterator mmebeta = mme_.begin();
         mmebeta != mme_.end(); ++mmebeta)
    {
      if ((*mmealpha).first >= (*mmebeta).first)
        continue;  // we  create couplings just for idxa < idxb

      std::pair<int, int> key = std::pair<int, int>((*mmealpha).first, (*mmebeta).first);

      coup_[key] = Teuchos::rcp(new ADAPTER::Coupling());

      std::map<int, Teuchos::RCP<const DRT::Discretization>>::iterator alphadis =
          dis.find((*mmealpha).first);
      if (alphadis == dis.end())
        dserror("Couldn't find discretization for key %d", (*mmealpha).first);
      std::map<int, Teuchos::RCP<const DRT::Discretization>>::iterator betadis =
          dis.find((*mmebeta).first);
      if (betadis == dis.end())
        dserror("Couldn't find discretization for key %d", (*mmebeta).first);

      coup_[key]->SetupConditionCoupling(*(*alphadis).second, (*mmealpha).second->Map(1),
          *(*betadis).second, (*mmebeta).second->Map(1), cond_name_, enddim_ - startdim_, true);
    }
  }
  return;
}

/*------------------------------------------------------------------------------------------------*
| Setup Couplings between different discretizations full coupling                  -- ager 07/2016|
*------------------------------------------------------------------------------------------------*/
void XFEM::Coupling_Comm_Manager::SetupFullCouplings(
    std::map<int, Teuchos::RCP<const DRT::Discretization>> dis)
{
  if (dis.size() < 2) return;

  for (uint idx_a = 0; idx_a < dis.size(); ++idx_a)
  {
    for (uint idx_b = 0; idx_b < dis.size(); ++idx_b)
    {
      if (idx_a >= idx_b) continue;  // we  create couplings just for idxa < idxb

      std::pair<int, int> key = std::pair<int, int>(idx_a, idx_b);

      coup_[key] = Teuchos::rcp(new ADAPTER::Coupling());

      std::map<int, Teuchos::RCP<const DRT::Discretization>>::iterator alphadis = dis.find(idx_a);
      if (alphadis == dis.end()) dserror("Couldn't find discretization for key %d", idx_a);
      std::map<int, Teuchos::RCP<const DRT::Discretization>>::iterator betadis = dis.find(idx_b);
      if (betadis == dis.end()) dserror("Couldn't find discretization for key %d", idx_b);

      coup_[key]->SetupCoupling(*(*alphadis).second, *(*betadis).second,
          *(*alphadis).second->NodeRowMap(), *(*betadis).second->NodeRowMap(), enddim_ - startdim_,
          false);
    }
  }
  return;
}

/*------------------------------------------------------------------------------------------------*
| Setup Full Extractor of all involved discretisations                                ager 03/2016|
*------------------------------------------------------------------------------------------------*/
void XFEM::Coupling_Comm_Manager::SetupFullExtractor(
    std::map<int, Teuchos::RCP<const DRT::Discretization>> dis)
{
  if (dis.size() < 2) return;

  std::vector<Teuchos::RCP<const Epetra_Map>> maps;
  for (std::map<int, Teuchos::RCP<const DRT::Discretization>>::iterator dit = dis.begin();
       dit != dis.end(); ++dit)
  {
    maps.push_back(Teuchos::rcp(new Epetra_Map(*(*dit).second->DofRowMap())));
  }

  Teuchos::RCP<Epetra_Map> fullmap = LINALG::MultiMapExtractor::MergeMaps(maps);
  fullextractor_ = Teuchos::rcp(new LINALG::MultiMapExtractor(*fullmap, maps));
}

/*------------------------------------------------------------------------------------------------*
| Get Coupling Converter between Discret A and B                                      ager 03/2016|
*------------------------------------------------------------------------------------------------*/
Teuchos::RCP<ADAPTER::CouplingConverter> XFEM::Coupling_Comm_Manager::GetCouplingConverter(
    int idxA, int idxB)
{
  if (idxA < idxB)
  {
    return Teuchos::rcp(new ADAPTER::CouplingMasterConverter(*GetCoupling(idxA, idxB)));
  }
  else if (idxA > idxB)
  {
    return Teuchos::rcp(new ADAPTER::CouplingSlaveConverter(*GetCoupling(idxB, idxA)));
  }
  else
  {
    dserror(
        "Coupling_Comm_Manager::GetCouplingConverter: Coupling Converter from %d to %d makes not "
        "really sense, does it?",
        idxA, idxB);
  }
  return Teuchos::null;  // 2 make Compiler happy :-)
}

/*------------------------------------------------------------------------------------------------*
| Get Coupling Object between Discret A and B                                      -- ager 03/2016|
*------------------------------------------------------------------------------------------------*/
Teuchos::RCP<ADAPTER::Coupling> XFEM::Coupling_Comm_Manager::GetCoupling(int idxA, int idxB)
{
  if (idxA < idxB)
  {
    std::pair<int, int> key = std::pair<int, int>(idxA, idxB);
    std::map<std::pair<int, int>, Teuchos::RCP<ADAPTER::Coupling>>::iterator cit = coup_.find(key);
    if (cit != coup_.end())
    {
      return cit->second;
    }
    else
    {
      dserror("Coupling_Comm_Manager::GetCoupling: Couldn't find Coupling for key (%d,%d)",
          key.first, key.second);
    }
  }
  else
  {
    dserror(
        "Coupling_Comm_Manager::GetCoupling: Just Coupling Objects for idxA < idxB are stored!");
  }
  return Teuchos::null;
}

/*------------------------------------------------------------------------------------------------*
| Get Transform Object between Discret A and B                                        ager 03/2016|
*------------------------------------------------------------------------------------------------*/
Teuchos::RCP<FSI::UTILS::MatrixLogicalSplitAndTransform> XFEM::Coupling_Comm_Manager::GetTransform(
    int transform_id)
{
  std::map<int, Teuchos::RCP<FSI::UTILS::MatrixLogicalSplitAndTransform>>::iterator tit =
      transform_.find(transform_id);
  if (tit != transform_.end() && transform_id != -1)
  {
    return tit->second;
  }
  else
  {
    transform_[transform_id] = Teuchos::rcp(new FSI::UTILS::MatrixLogicalSplitAndTransform());
    return transform_[transform_id];
  }
  return Teuchos::null;
}

/*------------------------------------------------------------------------------------------------*
| Get Map Extractor Object for idx                                                    ager 03/2016|
*------------------------------------------------------------------------------------------------*/
Teuchos::RCP<LINALG::MultiMapExtractor> XFEM::Coupling_Comm_Manager::GetMapExtractor(int idx)
{
  std::map<int, Teuchos::RCP<LINALG::MultiMapExtractor>>::iterator mit = mme_.find(idx);
  if (mit != mme_.end())
  {
    return mit->second;
  }
  else
  {
    dserror(
        "Coupling_Comm_Manager::GetMapExtractor: Couldn't find Map Extractor for key (%d)", idx);
  }
  return Teuchos::null;  // to guarantee that the complier feels comfortable in his skin...
}

/*---------------------------------------------------------------------------------------------------*
| Coupling_Comm_Manager Debug Output                                                     ager
03/2016|
*---------------------------------------------------------------------------------------------------*/
void XFEM::Coupling_Comm_Manager::DebugOut(
    std::string str1, std::string str2, std::string str3, std::string str4)
{
#ifdef COUP_MANAGER_DEBUG_OUT
  std::cout << GRAY_LIGHT << str1 << str2 << str3 << str4 << END_COLOR << std::flush;
#endif
  return;
}
