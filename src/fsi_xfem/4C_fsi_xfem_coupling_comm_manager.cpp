/*----------------------------------------------------------------------*/
/*! \file
\brief  Coupling Communication Manager automatically creates all required coupling object to
transform matrixes, vectors, ...

\level 2


*----------------------------------------------------------------------*/

#include "4C_fsi_xfem_coupling_comm_manager.hpp"

#include "4C_coupling_adapter.hpp"
#include "4C_coupling_adapter_converter.hpp"
#include "4C_discretization_condition_selector.hpp"
#include "4C_lib_discret.hpp"
#include "4C_linalg_mapextractor.hpp"
#include "4C_linalg_matrixtransform.hpp"
#include "4C_linalg_sparsematrix.hpp"
#include "4C_linalg_utils_sparse_algebra_manipulation.hpp"

#include <iostream>
#include <string>

FOUR_C_NAMESPACE_OPEN


/*-----------------------------------------------------------------------------------------*
| constructor - simplified for 1 discretization to couple (public)             ager 03/2016|
*-----------------------------------------------------------------------------------------*/
XFEM::CouplingCommManager::CouplingCommManager(
    Teuchos::RCP<const DRT::Discretization> dis0, std::string cond_name, int startdim, int enddim)
    : cond_name_(cond_name), startdim_(startdim), enddim_(enddim)
{
  std::map<int, Teuchos::RCP<const DRT::Discretization>> dis_map;
  dis_map.insert(std::pair<int, Teuchos::RCP<const DRT::Discretization>>(0, dis0));

  setup(dis_map);
}

/*-----------------------------------------------------------------------------------------*
| constructor - simplified for 2 discretizations to couple (public)            ager 03/2016|
*-----------------------------------------------------------------------------------------*/
XFEM::CouplingCommManager::CouplingCommManager(Teuchos::RCP<const DRT::Discretization> dis0,
    Teuchos::RCP<const DRT::Discretization> dis1, std::string cond_name, int startdim, int enddim)
    : cond_name_(cond_name), startdim_(startdim), enddim_(enddim)
{
  std::map<int, Teuchos::RCP<const DRT::Discretization>> dis_map;
  dis_map.insert(std::pair<int, Teuchos::RCP<const DRT::Discretization>>(0, dis0));
  dis_map.insert(std::pair<int, Teuchos::RCP<const DRT::Discretization>>(1, dis1));

  setup(dis_map);
}

/*-----------------------------------------------------------------------------------------*
| constructor (public)                                                         ager 03/2016|
*-----------------------------------------------------------------------------------------*/
XFEM::CouplingCommManager::CouplingCommManager(
    std::map<int, Teuchos::RCP<const DRT::Discretization>> dis, std::string cond_name, int startdim,
    int enddim)
    : cond_name_(cond_name), startdim_(startdim), enddim_(enddim)
{
  setup(dis);
}


/*--------------------------------------------------------------------------------------------------------------*
| Transfer conditioned part of Vector from Discretization A --> B with different transfer types ager
03/2016|
*--------------------------------------------------------------------------------------------------------------*/
void XFEM::CouplingCommManager::InsertVector(const int idxA, Teuchos::RCP<const Epetra_Vector> vecA,
    const int idxB, Teuchos::RCP<Epetra_Vector> vecB, const CouplingCommManager::TransferType ttype,
    bool add, double scale)
{
  switch (ttype)
  {
    case CouplingCommManager::full_to_full:
    {
      Teuchos::RCP<CORE::LINALG::MultiMapExtractor> mmeb = GetMapExtractor(idxB);
      Teuchos::RCP<Epetra_Vector> tmpvec = Teuchos::rcp(new Epetra_Vector(*mmeb->Map(1), true));
      InsertVector(idxA, vecA, idxB, tmpvec, CouplingCommManager::full_to_partial, false, scale);
      if (!add)
        mmeb->InsertVector(tmpvec, 1, vecB);
      else
        mmeb->AddVector(tmpvec, 1, vecB);
      break;
    }
    case CouplingCommManager::full_to_partial:
    {
      InsertVector(idxA, GetMapExtractor(idxA)->ExtractVector(vecA, 1), idxB, vecB,
          CouplingCommManager::partial_to_partial, add, scale);
      break;
    }
    case CouplingCommManager::partial_to_full:
    {
      Teuchos::RCP<CORE::LINALG::MultiMapExtractor> mmeb = GetMapExtractor(idxB);
      Teuchos::RCP<Epetra_Vector> tmpvec = Teuchos::rcp(new Epetra_Vector(*mmeb->Map(1), true));
      InsertVector(idxA, vecA, idxB, tmpvec, CouplingCommManager::partial_to_partial, false, scale);
      if (!add)
        mmeb->InsertVector(tmpvec, 1, vecB);
      else
        mmeb->AddVector(tmpvec, 1, vecB);
      break;
    }
    case CouplingCommManager::partial_to_partial:
    {
      if (vecB == Teuchos::null)
        FOUR_C_THROW("Coupling_Comm_Manager::InsertVector: vecB is Teuchos::null!");
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
    case CouplingCommManager::partial_to_global:
    {
      Teuchos::RCP<CORE::LINALG::MultiMapExtractor> mme = GetFullMapExtractor();
      Teuchos::RCP<Epetra_Vector> fullvec = Teuchos::rcp(new Epetra_Vector(*mme->Map(idxB), true));
      InsertVector(idxA, vecA, idxB, fullvec, CouplingCommManager::partial_to_full, false, scale);
      if (!add)
        mme->InsertVector(fullvec, idxB, vecB);
      else
        mme->AddVector(fullvec, idxB, vecB);
      break;
    }
    case CouplingCommManager::full_to_global:
    {
      Teuchos::RCP<CORE::LINALG::MultiMapExtractor> mme = GetFullMapExtractor();
      Teuchos::RCP<Epetra_Vector> fullvec = Teuchos::rcp(new Epetra_Vector(*mme->Map(idxB), true));
      InsertVector(idxA, vecA, idxB, fullvec, CouplingCommManager::full_to_full, false, scale);
      if (!add)
        mme->InsertVector(fullvec, idxB, vecB);
      else
        mme->AddVector(fullvec, idxB, vecB);
      break;
    }
    default:
      FOUR_C_THROW("Coupling_Comm_Manager::InsertVector: Transfer Type not implemented!");
  }
  return;
}

/*----------------------------------------------------------------------------------------------*
| //! Insert a Matrix A into Matrix B (choose type of transfer, add or scaling)     ager 03/2016|
*----------------------------------------------------------------------------------------------*/
bool XFEM::CouplingCommManager::InsertMatrix(int transform_id, int idxA,
    const CORE::LINALG::SparseMatrix& matA, int idxB, CORE::LINALG::SparseMatrix& matB,
    const CouplingCommManager::MatrixTransferType mttype, double scale, bool exactmatch,
    bool addmatrix)
{
  switch (mttype)
  {
    case CouplingCommManager::col:
    {
      return GetTransform(transform_id)
          ->
          operator()(matA, matA.RangeMap(), matA.DomainMap(), scale, nullptr,
              get_coupling_converter(idxA, idxB).getRawPtr(), matB, exactmatch, addmatrix);
      break;
    }
    case CouplingCommManager::row:
    {
      return GetTransform(transform_id)
          ->
          operator()(matA, matA.RangeMap(), matA.DomainMap(), scale,
              get_coupling_converter(idxA, idxB).getRawPtr(), nullptr, matB, true, addmatrix);
      break;
    }
    case CouplingCommManager::row_and_col:
    {
      return GetTransform(transform_id)
          ->
          operator()(matA, matA.RangeMap(), matA.DomainMap(), scale,
              get_coupling_converter(idxA, idxB).getRawPtr(),
              get_coupling_converter(idxA, idxB).getRawPtr(), matB, exactmatch, addmatrix);
      break;
    }
    default:
      FOUR_C_THROW("Coupling_Comm_Manager::InsertMatrix: Matrix Transfer Type not implemented!");
  }
  return true;
}

/*-----------------------------------------------------------------------------------------*
| Setup Coupling_Comm_Manager                                                  ager 06/2016|
*-----------------------------------------------------------------------------------------*/
void XFEM::CouplingCommManager::setup(std::map<int, Teuchos::RCP<const DRT::Discretization>> dis)
{
  if (cond_name_ != "")  // Setup for Communication on Condition
  {
    setup_multi_map_extractors(dis);
    setup_couplings(dis);
  }
  else  // Setup for Communication on full Discretization
  {
    setup_full_couplings(dis);       // First Couple full discretizations (from startdim to enddim)
    setup_full_map_extractors(dis);  // Setup Extractors between Couplingpart of discretization and
                                     // full discretization (e.g. pres <==> vel&pres)
  }
  setup_full_extractor(dis);
  return;
}

/*-----------------------------------------------------------------------------------------*
| Setup MultiMapExtractors for all fields                                      ager 03/2016|
*-----------------------------------------------------------------------------------------*/
void XFEM::CouplingCommManager::setup_multi_map_extractors(
    std::map<int, Teuchos::RCP<const DRT::Discretization>> dis)
{
  for (std::map<int, Teuchos::RCP<const DRT::Discretization>>::iterator dit = dis.begin();
       dit != dis.end(); ++dit)
  {
    Teuchos::RCP<CORE::Conditions::MultiConditionSelector> mcs =
        Teuchos::rcp(new CORE::Conditions::MultiConditionSelector());
    mme_[dit->first] = Teuchos::rcp(new CORE::LINALG::MultiMapExtractor());
    mcs->AddSelector(Teuchos::rcp(
        new CORE::Conditions::NDimConditionSelector(*dit->second, cond_name_, startdim_, enddim_)));
    mcs->SetupExtractor(*dit->second, *dit->second->dof_row_map(), *mme_[dit->first]);
  }
}

/*-----------------------------------------------------------------------------------------*
| Setup MultiMapExtractors for all fields                                      ager 03/2016|
*-----------------------------------------------------------------------------------------*/
void XFEM::CouplingCommManager::setup_full_map_extractors(
    std::map<int, Teuchos::RCP<const DRT::Discretization>> dis)
{
  if (dis.size() < 2)
    FOUR_C_THROW(
        "setup_full_map_extractors: Just extract from discretization to coupled map! (e.g. Fluid "
        "vel "
        "<==> Fluid vel&pres)");

  for (std::map<int, Teuchos::RCP<const DRT::Discretization>>::iterator dit = dis.begin();
       dit != dis.end(); ++dit)
  {
    Teuchos::RCP<CORE::LINALG::MapExtractor> me = Teuchos::rcp(new CORE::LINALG::MapExtractor());
    if (static_cast<std::size_t>(dit->first) < dis.size() - 1)
    {
      Teuchos::RCP<CORE::ADAPTER::Coupling> coup = GetCoupling(dit->first, dit->first + 1);
      me->Setup(*dit->second->dof_row_map(), coup->MasterDofMap(),
          CORE::LINALG::SplitMap(*dit->second->dof_row_map(), *coup->MasterDofMap()));
    }
    else
    {
      Teuchos::RCP<CORE::ADAPTER::Coupling> coup = GetCoupling(dit->first - 1, dit->first);
      me->Setup(*dit->second->dof_row_map(), coup->SlaveDofMap(),
          CORE::LINALG::SplitMap(*dit->second->dof_row_map(), *coup->SlaveDofMap()));
    }
    mme_[dit->first] = me;
  }
  return;
}

/*------------------------------------------------------------------------------------------------*
| Setup Couplings between different discretizations with the same condition on in  -- ager 03/2016|
*------------------------------------------------------------------------------------------------*/
void XFEM::CouplingCommManager::setup_couplings(
    std::map<int, Teuchos::RCP<const DRT::Discretization>> dis)
{
  if (dis.size() < 2) return;

  for (std::map<int, Teuchos::RCP<CORE::LINALG::MultiMapExtractor>>::iterator mmealpha =
           mme_.begin();
       mmealpha != mme_.end(); ++mmealpha)
  {
    for (std::map<int, Teuchos::RCP<CORE::LINALG::MultiMapExtractor>>::iterator mmebeta =
             mme_.begin();
         mmebeta != mme_.end(); ++mmebeta)
    {
      if ((*mmealpha).first >= (*mmebeta).first)
        continue;  // we  create couplings just for idxa < idxb

      std::pair<int, int> key = std::pair<int, int>((*mmealpha).first, (*mmebeta).first);

      coup_[key] = Teuchos::rcp(new CORE::ADAPTER::Coupling());

      std::map<int, Teuchos::RCP<const DRT::Discretization>>::iterator alphadis =
          dis.find((*mmealpha).first);
      if (alphadis == dis.end())
        FOUR_C_THROW("Couldn't find discretization for key %d", (*mmealpha).first);
      std::map<int, Teuchos::RCP<const DRT::Discretization>>::iterator betadis =
          dis.find((*mmebeta).first);
      if (betadis == dis.end())
        FOUR_C_THROW("Couldn't find discretization for key %d", (*mmebeta).first);

      coup_[key]->setup_condition_coupling(*(*alphadis).second, (*mmealpha).second->Map(1),
          *(*betadis).second, (*mmebeta).second->Map(1), cond_name_, enddim_ - startdim_, true);
    }
  }
  return;
}

/*------------------------------------------------------------------------------------------------*
| Setup Couplings between different discretizations full coupling                  -- ager 07/2016|
*------------------------------------------------------------------------------------------------*/
void XFEM::CouplingCommManager::setup_full_couplings(
    std::map<int, Teuchos::RCP<const DRT::Discretization>> dis)
{
  if (dis.size() < 2) return;

  for (std::size_t idx_a = 0; idx_a < dis.size(); ++idx_a)
  {
    for (std::size_t idx_b = 0; idx_b < dis.size(); ++idx_b)
    {
      if (idx_a >= idx_b) continue;  // we  create couplings just for idxa < idxb

      std::pair<int, int> key = std::pair<int, int>(idx_a, idx_b);

      coup_[key] = Teuchos::rcp(new CORE::ADAPTER::Coupling());

      std::map<int, Teuchos::RCP<const DRT::Discretization>>::iterator alphadis = dis.find(idx_a);
      if (alphadis == dis.end()) FOUR_C_THROW("Couldn't find discretization for key %d", idx_a);
      std::map<int, Teuchos::RCP<const DRT::Discretization>>::iterator betadis = dis.find(idx_b);
      if (betadis == dis.end()) FOUR_C_THROW("Couldn't find discretization for key %d", idx_b);

      coup_[key]->setup_coupling(*(*alphadis).second, *(*betadis).second,
          *(*alphadis).second->NodeRowMap(), *(*betadis).second->NodeRowMap(), enddim_ - startdim_,
          false);
    }
  }
  return;
}

/*------------------------------------------------------------------------------------------------*
| Setup Full Extractor of all involved discretisations                                ager 03/2016|
*------------------------------------------------------------------------------------------------*/
void XFEM::CouplingCommManager::setup_full_extractor(
    std::map<int, Teuchos::RCP<const DRT::Discretization>> dis)
{
  if (dis.size() < 2) return;

  std::vector<Teuchos::RCP<const Epetra_Map>> maps;
  for (std::map<int, Teuchos::RCP<const DRT::Discretization>>::iterator dit = dis.begin();
       dit != dis.end(); ++dit)
  {
    maps.push_back(Teuchos::rcp(new Epetra_Map(*(*dit).second->dof_row_map())));
  }

  Teuchos::RCP<Epetra_Map> fullmap = CORE::LINALG::MultiMapExtractor::MergeMaps(maps);
  fullextractor_ = Teuchos::rcp(new CORE::LINALG::MultiMapExtractor(*fullmap, maps));
}

/*------------------------------------------------------------------------------------------------*
| Get Coupling Converter between Discret A and B                                      ager 03/2016|
*------------------------------------------------------------------------------------------------*/
Teuchos::RCP<CORE::ADAPTER::CouplingConverter> XFEM::CouplingCommManager::get_coupling_converter(
    int idxA, int idxB)
{
  if (idxA < idxB)
  {
    return Teuchos::rcp(new CORE::ADAPTER::CouplingMasterConverter(*GetCoupling(idxA, idxB)));
  }
  else if (idxA > idxB)
  {
    return Teuchos::rcp(new CORE::ADAPTER::CouplingSlaveConverter(*GetCoupling(idxB, idxA)));
  }
  else
  {
    FOUR_C_THROW(
        "Coupling_Comm_Manager::get_coupling_converter: Coupling Converter from %d to %d makes not "
        "really sense, does it?",
        idxA, idxB);
  }
  return Teuchos::null;  // 2 make Compiler happy :-)
}

/*------------------------------------------------------------------------------------------------*
| Get Coupling Object between Discret A and B                                      -- ager 03/2016|
*------------------------------------------------------------------------------------------------*/
Teuchos::RCP<CORE::ADAPTER::Coupling> XFEM::CouplingCommManager::GetCoupling(int idxA, int idxB)
{
  if (idxA < idxB)
  {
    std::pair<int, int> key = std::pair<int, int>(idxA, idxB);
    std::map<std::pair<int, int>, Teuchos::RCP<CORE::ADAPTER::Coupling>>::iterator cit =
        coup_.find(key);
    if (cit != coup_.end())
    {
      return cit->second;
    }
    else
    {
      FOUR_C_THROW("Coupling_Comm_Manager::GetCoupling: Couldn't find Coupling for key (%d,%d)",
          key.first, key.second);
    }
  }
  else
  {
    FOUR_C_THROW(
        "Coupling_Comm_Manager::GetCoupling: Just Coupling Objects for idxA < idxB are stored!");
  }
  return Teuchos::null;
}

/*------------------------------------------------------------------------------------------------*
| Get Transform Object between Discret A and B                                        ager 03/2016|
*------------------------------------------------------------------------------------------------*/
Teuchos::RCP<CORE::LINALG::MatrixLogicalSplitAndTransform> XFEM::CouplingCommManager::GetTransform(
    int transform_id)
{
  std::map<int, Teuchos::RCP<CORE::LINALG::MatrixLogicalSplitAndTransform>>::iterator tit =
      transform_.find(transform_id);
  if (tit != transform_.end() && transform_id != -1)
  {
    return tit->second;
  }
  else
  {
    transform_[transform_id] = Teuchos::rcp(new CORE::LINALG::MatrixLogicalSplitAndTransform());
    return transform_[transform_id];
  }
  return Teuchos::null;
}

/*------------------------------------------------------------------------------------------------*
| Get Map Extractor Object for idx                                                    ager 03/2016|
*------------------------------------------------------------------------------------------------*/
Teuchos::RCP<CORE::LINALG::MultiMapExtractor> XFEM::CouplingCommManager::GetMapExtractor(int idx)
{
  std::map<int, Teuchos::RCP<CORE::LINALG::MultiMapExtractor>>::iterator mit = mme_.find(idx);
  if (mit != mme_.end())
  {
    return mit->second;
  }
  else
  {
    FOUR_C_THROW(
        "Coupling_Comm_Manager::GetMapExtractor: Couldn't find Map Extractor for key (%d)", idx);
  }
  return Teuchos::null;  // to guarantee that the complier feels comfortable in his skin...
}

/*---------------------------------------------------------------------------------------------------*
| Coupling_Comm_Manager Debug Output                                                     ager
03/2016|
*---------------------------------------------------------------------------------------------------*/
void XFEM::CouplingCommManager::DebugOut(
    std::string str1, std::string str2, std::string str3, std::string str4)
{
#ifdef COUP_MANAGER_DEBUG_OUT
  std::cout << str1 << str2 << str3 << str4 << std::flush;
#endif
  return;
}

FOUR_C_NAMESPACE_CLOSE
