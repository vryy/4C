/*----------------------------------------------------------------------------*/
/**
\file xstr_multi_discretization_wrapper.cpp

\brief Implementation of a multi discretization wrapper for eXtended structural
       problems

\maintainer Michael Hiermeier

\date Oct 6, 2016

\level 3

*/
/*----------------------------------------------------------------------------*/

#include "xstr_multi_discretization_wrapper.H"
#include "xstr_multi_io.H"

#include "../drt_lib/drt_node.H"
#include "../linalg/linalg_utils.H"
#include "../linalg/linalg_blocksparsematrix.H"

#include "../drt_lib/drt_utils_discret.H"

#include <boost/bind.hpp>
#include <Epetra_SerialDenseVector.h>

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
XSTR::MultiDiscretizationWrapper::MultiDiscretizationWrapper()
    : isinit_(false),
      issetup_(false),
      isfilled_(false),
      max_num_reserved_dofs_per_node_(-1)
{
  // Intentionally left blank.
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XSTR::MultiDiscretizationWrapper::Init(
    const std::string& name,
    const Teuchos::RCP<DRT::DiscretizationInterface>& discret_ptr,
    const int & max_num_reserved_dofs_per_node)
{
  XDisVec discret_vec(1,discret_ptr);
  Init(name,discret_vec,max_num_reserved_dofs_per_node);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XSTR::MultiDiscretizationWrapper::Init(
    const std::string& name,
    const XDisVec& discret_vec,
    const int & max_num_reserved_dofs_per_node)
{
  issetup_ = false;

  name_ = name;

  std::vector<Teuchos::RCP<DRT::DiscretizationInterface> >::const_iterator cit_vec;
  for (cit_vec=discret_vec.begin();cit_vec!=discret_vec.end();++cit_vec)
    AddDiscretization(*cit_vec);

  max_num_reserved_dofs_per_node_ = max_num_reserved_dofs_per_node;

  isinit_ = true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XSTR::MultiDiscretizationWrapper::Setup()
{
  CheckInit();
  isfilled_ = false;

  cXDisVec discret_vec(0);
  // copy the second value of the map into a vector
  std::transform( discret_map_.begin(), discret_map_.end(),
      std::back_inserter(discret_vec),
      boost::bind(&XDisMap::value_type::second,_1) );

  // build the multi field map extractor
  map_extractor_ = Teuchos::rcp(new XFEM::MultiFieldMapExtractor());
  map_extractor_->Init(discret_vec,MaxNumReservedDofsPerNode());

  // setup map extractor and build class internal maps
  FillComplete(false,false,false,true,true);

  /* Create a new discretization writer wrapper, if we added a new
   * discretization or anything else happened, which made it necessary to
   * call Setup(). */
  writers_wrapper_ = Teuchos::rcp<XSTR::IO::MultiDiscretizationWriter>(new
      XSTR::IO::MultiDiscretizationWriter(Teuchos::rcp(this,false)));

  issetup_ = true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XSTR::MultiDiscretizationWrapper::AddDiscretization(
    const Teuchos::RCP<DRT::DiscretizationInterface>& discret_ptr,
    const bool & check_insertion)
{
  // if you add a new discretization, Setup() has to be called afterwards
  issetup_ = false;

  const enum XFEM::FieldName fieldname =
      XFEM::String2FieldName(discret_ptr->Name());

  if (check_insertion)
  {
    std::pair<XDisMap::iterator,bool> ret =
        discret_map_.insert(XDisPair(fieldname,discret_ptr));

    if (not ret.second)
      dserror("The discretization \"%s\" was already inserted!",
          discret_ptr->Name().c_str());
  }
  else
    discret_map_[fieldname] = discret_ptr;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XSTR::MultiDiscretizationWrapper::ExtractVector(
    const Epetra_MultiVector & full,
    enum XFEM::FieldName       field,
    Epetra_MultiVector &       partial,
    enum ::IO::VectorType      vt) const
{
  switch (vt)
  {
    case ::IO::dofvector:
    {
      map_extractor_->ExtractVector(full,field,partial, XFEM::map_dofs);
      break;
    }
    case ::IO::nodevector:
    {
      map_extractor_->ExtractVector(full,field,partial, XFEM::map_nodes);
      break;
    }
    case ::IO::elementvector:
    {
      map_extractor_->ExtractElementVector(full,field,partial);
      break;
    }
    default:
    {
      dserror("Unknown ::IO::VectorType! (enum = %d)",vt);
      exit(EXIT_FAILURE);
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XSTR::MultiDiscretizationWrapper::InsertVector(
    const Epetra_MultiVector & partial,
    enum XFEM::FieldName       field,
    Epetra_MultiVector &       full,
    enum ::IO::VectorType      vt) const
{
  switch (vt)
  {
    case ::IO::dofvector:
    {
      map_extractor_->InsertVector(partial,field,full, XFEM::map_dofs);
      break;
    }
    case ::IO::nodevector:
    {
      map_extractor_->InsertVector(partial,field,full, XFEM::map_nodes);
      break;
    }
    case ::IO::elementvector:
    {
      map_extractor_->InsertElementVector(partial,field,full);
      break;
    }
    default:
    {
      dserror("Unknown ::IO::VectorType! (enum = %d)",vt);
      exit(EXIT_FAILURE);
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XSTR::MultiDiscretizationWrapper::AddVector(
    const Epetra_MultiVector & partial,
    enum XFEM::FieldName       field,
    Epetra_MultiVector &       full,
    double                     scale,
    enum ::IO::VectorType      vt) const
{
  switch (vt)
  {
    case ::IO::dofvector:
    {
      map_extractor_->AddVector(partial,field,full,scale, XFEM::map_dofs);
      break;
    }
    case ::IO::nodevector:
    {
      map_extractor_->AddVector(partial,field,full,scale, XFEM::map_nodes);
      break;
    }
    case ::IO::elementvector:
    {
      map_extractor_->AddElementVector(partial,field,full,scale);
      break;
    }
    default:
    {
      dserror("Unknown ::IO::VectorType! (enum = %d)",vt);
      exit(EXIT_FAILURE);
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XSTR::MultiDiscretizationWrapper::Reset()
{
  merged_node_.clear();
  merged_node_row_map_ = Teuchos::null;
  merged_node_row_ptr_.clear();

  merged_element_.clear();
  merged_ele_row_map_ = Teuchos::null;
  merged_element_row_ptr_.clear();

  isfilled_ = false;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<IO::DiscretizationWriter> XSTR::MultiDiscretizationWrapper::
    Writer() const
{
  CheckInitSetup();
  return writers_wrapper_;
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool XSTR::MultiDiscretizationWrapper::IsXFEM( enum XFEM::FieldName field )
const
{
 return map_extractor_->IsXFemDis(field);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XSTR::MultiDiscretizationWrapper::Print(std::ostream& os) const
{
  XDisMap::const_iterator cit;
  for (cit=discret_map_.begin();cit!=discret_map_.end();++cit)
  {
    if (Comm().MyPID() == 0)
    {
      os << "\n\n::::::::::::::::::::::::::::::::::::::::::::::::::::::::::\n";
      os << "            " << XFEM::FieldName2String(cit->first) <<
            "-DISCRETIZATION";
      os << "\n::::::::::::::::::::::::::::::::::::::::::::::::::::::::::\n\n";
    }
    cit->second->Print(os);
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
int XSTR::MultiDiscretizationWrapper::FillComplete(
    bool assigndegreesoffreedom,
    bool initelements,
    bool doboundaryconditions,
    bool buildsystemmaps,
    bool setupmapextractor)
{
  CheckInit();
  Reset();

  // --------------------------------------------------------------------------
  // first call FillComplete() on each single discretization
  // --------------------------------------------------------------------------
  int ret = 0;
  XDisMap::const_iterator cit;
  for (cit=discret_map_.begin();cit!=discret_map_.end();++cit)
  {
    ret = cit->second->FillComplete(assigndegreesoffreedom,initelements,
        doboundaryconditions);
    if (ret!=0)
      dserror("FillComplete failed on discretization %s with the "
          "error code %d!",cit->second->Name().c_str(),ret);
  }

  if (setupmapextractor)
    map_extractor_->Setup();

  if (buildsystemmaps)
  {
    // --------------------------------------------------------------------------
    // Build merged nodal information
    // --------------------------------------------------------------------------
    BuildMergedNodes();
    BuildMergedNodeRowMap();
    BuildMergedNodeRowPtr();

    // --------------------------------------------------------------------------
    // Build merged element information
    // --------------------------------------------------------------------------
    BuildMergedElements();
    BuildMergedElementRowMap();
    BuildMergedElementRowPtr();
  }

  isfilled_ = true;

  return ret;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XSTR::MultiDiscretizationWrapper::BuildMergedNodes()
{
  // ---------------------------------------------------------------------------
  // (1) First we add the nodes of the auxiliary interface discretization.
  //     Note that we do not access the master node row map's of each single
  //     wrapped discretization, because these maps are overlapping!
  // ---------------------------------------------------------------------------
  int * ngids = map_extractor_->INodeRowMap()->MyGlobalElements();
  int my_num_nodes = map_extractor_->INodeRowMap()->NumMyElements();
  for (int nlid=0; nlid<my_num_nodes; ++nlid)
  {
    std::pair<std::map<int,Teuchos::RCP<DRT::Node> >::iterator, bool> check =
        merged_node_.insert(std::pair<int,Teuchos::RCP<DRT::Node> >(
            ngids[nlid],Teuchos::rcp<DRT::Node>(
                map_extractor_->gINode(ngids[nlid]),false)));

    if (not check.second)
      dserror("The interface row node with GID %d could not be inserted, because "
          "the node GID is not unique!", ngids[nlid]);
  }

  // ---------------------------------------------------------------------------
  // (2) In a second attempt we add the remaining non-interface nodes of all
  //     wrapped discretizations
  // ---------------------------------------------------------------------------
  XDisMap::const_iterator cit;
  for (cit=discret_map_.begin();cit!=discret_map_.end();++cit)
  {
    ngids = map_extractor_->NodeRowMap(cit->first,
        XFEM::MULTIFIELD::block_non_interface)->MyGlobalElements();
    my_num_nodes = map_extractor_->NodeRowMap(cit->first,
        XFEM::MULTIFIELD::block_non_interface)->NumMyElements();
    for (int nlid=0;nlid<my_num_nodes;++nlid)
    {
      std::pair<std::map<int,Teuchos::RCP<DRT::Node> >::iterator, bool> check =
              merged_node_.insert(std::pair<int,Teuchos::RCP<DRT::Node> >(
                  ngids[nlid],Teuchos::rcp<DRT::Node>(
                      cit->second->gNode(ngids[nlid]),false)));

      if (not check.second)
        dserror("The non-interface row node with GID %d could not be inserted, "
            "because the node GID is not unique!", ngids[nlid]);
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XSTR::MultiDiscretizationWrapper::BuildMergedNodeRowMap()
{
  // --------------------------------------------------------------------------
  // (1) get the node row map of the auxiliary interface discretization first
  // --------------------------------------------------------------------------
  const Epetra_Map& inoderowmap = *(map_extractor_->INodeRowMap());
  merged_node_row_map_ = Teuchos::rcp(new Epetra_Map(inoderowmap));
  // --------------------------------------------------------------------------
  // (2) merge the interface row map with the non-interface node row maps
  //     of the wrapped discretizations
  // --------------------------------------------------------------------------
  XDisMap::const_iterator cit;
  for (cit=discret_map_.begin();cit!=discret_map_.end();++cit)
  {
    merged_node_row_map_ = LINALG::MergeMap(*merged_node_row_map_,
        *map_extractor_->NodeRowMap(cit->first,
            XFEM::MULTIFIELD::block_non_interface),false);
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XSTR::MultiDiscretizationWrapper::BuildMergedNodeRowPtr()
{
  int* ngids = merged_node_row_map_->MyGlobalElements();
  int my_num_nodes = merged_node_row_map_->NumMyElements();
  merged_node_row_ptr_.resize(my_num_nodes,NULL);

  for (int nlid=0;nlid<my_num_nodes;++nlid)
  {
    merged_node_row_ptr_[nlid] = merged_node_[ngids[nlid]].get();
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XSTR::MultiDiscretizationWrapper::BuildMergedElements()
{
  // add all row elements of the wrapped discretizations
  XDisMap::const_iterator cit;
  for (cit=discret_map_.begin();cit!=discret_map_.end();++cit)
  {
    int num_my_row_eles = cit->second->ElementRowMap()->NumMyElements();
    int* egids = cit->second->ElementRowMap()->MyGlobalElements();
    for (int lid=0;lid<num_my_row_eles;++lid)
    {
      std::pair<std::map<int,Teuchos::RCP<DRT::Element> >::iterator, bool> check =
          merged_element_.insert(std::pair<int,Teuchos::RCP<DRT::Element> >(
              egids[lid],Teuchos::rcp<DRT::Element>(
                  cit->second->gElement(egids[lid]),false)));

      if (not check.second)
        dserror("The row element with GID %d could not be inserted, because the"
            "element GID is not unique!", egids[lid]);
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XSTR::MultiDiscretizationWrapper::BuildMergedElementRowMap()
{
  std::vector<int> my_ele_gids;
  my_ele_gids.reserve(merged_element_.size());

  std::map<int,Teuchos::RCP<DRT::Element> >::const_iterator cit;
  for (cit=merged_element_.begin();cit!=merged_element_.end();++cit)
  {
    // sanity check
    if (cit->second->Owner()==Comm().MyPID())
      my_ele_gids.push_back(cit->first);
    else
      dserror("The merged_element map is supposed to contain only "
          "row elements (MultiDiscretizationWrapper). Something went "
          "wrong at element GID = %d!",cit->first);
  }
  // create a new element row map
  merged_ele_row_map_ = Teuchos::rcp(new Epetra_Map(-1,
      static_cast<int>(my_ele_gids.size()),&my_ele_gids[0], 0, Comm()));
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XSTR::MultiDiscretizationWrapper::BuildMergedElementRowPtr()
{
  int* ele_gids = merged_ele_row_map_->MyGlobalElements();
  int my_num_ele = merged_ele_row_map_->NumMyElements();
  merged_element_row_ptr_.resize(my_num_ele,NULL);

  for (int elid=0; elid<my_num_ele; ++elid)
  {
    merged_element_row_ptr_[elid] = merged_element_.at(ele_gids[elid]).get();
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const Epetra_Map* XSTR::MultiDiscretizationWrapper::NodeRowMap() const
{
  if (not Filled())
    dserror("FillComplete() must be called before call to NodeRowMap()");

  return merged_node_row_map_.get();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const Epetra_Map* XSTR::MultiDiscretizationWrapper::ElementRowMap() const
{
  if (not Filled())
    dserror("FillComplete() must be called before call to ElementRowMap()");

  return merged_ele_row_map_.get();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
int XSTR::MultiDiscretizationWrapper::NumMyRowNodes() const
{
  if (not Filled())
    dserror("FillComplete() must be called before call to NumMyRowNodes()");

  return merged_node_row_map_->NumMyElements();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool XSTR::MultiDiscretizationWrapper::HaveGlobalNode(int gid) const
{
  std::map<int,Teuchos::RCP<DRT::Node> >::const_iterator curr =
      merged_node_.find(gid);
  return (curr != merged_node_.end());
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
DRT::Node* XSTR::MultiDiscretizationWrapper::gNode(int gid) const
{
  std::map<int, Teuchos::RCP<DRT::Node> >::const_iterator cit =
      merged_node_.find(gid);

  if (cit == merged_node_.end())
    dserror("No direct access to the node with GID %d. The merged_node_ map "
        "is supposed to contain only row elements (MultiDiscretizationWrapper).",
        gid);

  return cit->second.get();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
DRT::Node* XSTR::MultiDiscretizationWrapper::lRowNode(int lid) const
{
  if (not Filled())
    dserror("FillComplete() must be called before call to lRowNode()");

  return merged_node_row_ptr_.at(lid);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
int XSTR::MultiDiscretizationWrapper::NumMyRowElements() const
{
  int mynumall = 0;
  XDisMap::const_iterator cit;
  for (cit=discret_map_.begin();cit!=discret_map_.end();++cit)
  {
    int mynum = cit->second->NumMyRowElements();
    if (mynum==-1)
      dserror("The discretization %s returned a number of my row elements of"
          " -1!",cit->second->Name().c_str());
    mynumall += mynum;
  }
  return mynumall;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
DRT::Element* XSTR::MultiDiscretizationWrapper::lRowElement(int lid) const
{
  if (not Filled())
    dserror("FillComplete() must be called before call to lRowElement()");

  return merged_element_row_ptr_[lid];
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
DRT::Element* XSTR::MultiDiscretizationWrapper::gElement(int gid) const
{
  std::map<int, Teuchos::RCP<DRT::Element> >::const_iterator cit =
      merged_element_.find(gid);

  if (cit == merged_element_.end())
    dserror("No direct access to the element with GID %d. The merged_element_ "
        "map is supposed to contain only row elements "
        "(MultiDiscretizationWrapper).", gid);

  return cit->second.get();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const Epetra_Map* XSTR::MultiDiscretizationWrapper::DofRowMap(unsigned nds)
const
{
  CheckInit();

  if (nds>0)
    dserror("The MultiDiscretizationWrapper provides only the combined \n"
            "DofRowMap, necessary to create state vectors and system \n"
            "matrices. To access the single DoF-sets, please access the \n"
            "wrapped discretization objects directly.");

  return map_extractor_->FullMap().get();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool XSTR::MultiDiscretizationWrapper::HaveDofs() const
{
  bool have_dofs = true;
  XDisMap::const_iterator cit;
  for (cit=discret_map_.begin();cit!=discret_map_.end();++cit)
  {
    if (not cit->second->HaveDofs())
      have_dofs = false;
  }
  return have_dofs;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
int XSTR::MultiDiscretizationWrapper::NumDofSets() const
{
  // the wrapper supports only one single dof set
  return 1;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
int XSTR::MultiDiscretizationWrapper::NumDof(unsigned dofset_id,
    const DRT::Node* node) const
{
  if (dofset_id>0)
    dserror("The MultiDiscretizationWrapper provides only the combined \n"
            "DofRowMap, necessary to create state vectors and system \n"
            "matrices. To access the single DoF-sets, please access the \n"
            "wrapped discretization objects directly.");

  // first check if the node is a considered part of the wrapped discretization
  const int ngid = node->Id();
  if (not NodeRowMap()->MyGID(ngid))
    dserror("The given node (GID: %d) is no part of the MultiDiscretizationWrapper. \n"
            "Maybe you have to call the routine directly on one of the wrapped \n"
            "discretizations to get what you want (i.e. ghosted nodes)?!", ngid);

  int numdof = 0;
  // interface nodes
  if (map_extractor_->INodeRowMap()->MyGID(ngid))
  {
    /* Necessary because the nodes were cloned during the creation of the
     * interface discretization. */
    DRT::Node * inode = map_extractor_->gINode(ngid);
    // access the interface dofset
    numdof = map_extractor_->INumDof(inode);
  }
  // non-interface nodes
  else
  {
    XDisMap::const_iterator cit;
    for (cit=discret_map_.begin();cit!=discret_map_.end();++cit)
    {
      if (cit->second->NodeRowMap()->MyGID(ngid))
      {
        numdof = cit->second->NumDof(0,node);
        break;
      }
    }
  }

  return numdof;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
int XSTR::MultiDiscretizationWrapper::NumStandardDof(const unsigned & dofset_id,
    const DRT::Node* node) const
{
  if (dofset_id>0)
    dserror("The MultiDiscretizationWrapper provides only the combined \n"
            "DofRowMap, necessary to create state vectors and system \n"
            "matrices. To access the single DoF-sets, please access the \n"
            "wrapped discretization objects directly.");

  // first check if the node is a considered part of the wrapped discretization
  const int ngid = node->Id();
  if (not NodeRowMap()->MyGID(ngid))
    dserror("The given node (GID: %d) is no part of the MultiDiscretizationWrapper. \n"
            "Maybe you have to call the routine directly on one of the wrapped \n"
            "discretizations to get what you want (i.e. ghosted nodes)?!", ngid);

  int numdof = 0;
  // interface nodes
  if (map_extractor_->INodeRowMap()->MyGID(ngid))
  {
    // access the interface dofset
    numdof = map_extractor_->INumStandardDof();
  }
  // non-interface nodes
  else
  {
    XDisMap::const_iterator cit;
    for (cit=discret_map_.begin();cit!=discret_map_.end();++cit)
    {
      if (cit->second->NodeRowMap()->MyGID(ngid))
      {
        numdof = cit->second->NumStandardDof(dofset_id,node);
        break;
      }
    }
  }
  return numdof;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
int XSTR::MultiDiscretizationWrapper::Dof(
    unsigned dofset_id, const DRT::Node* node, const int ldof) const
{
  if (dofset_id>0)
    dserror("The MultiDiscretizationWrapper provides only the combined \n"
            "DofRowMap, necessary to create state vectors and system \n"
            "matrices. To access the single DoF-sets, please access the \n"
            "wrapped discretization objects directly.");

  // first check if the node is a considered part of the wrapped discretization
  const int ngid = node->Id();
  if (not NodeRowMap()->MyGID(ngid))
    dserror("The given node (GID: %d) is no part of the MultiDiscretizationWrapper. \n"
            "Maybe you have to call the routine directly on one of the wrapped \n"
            "discretizations to get what you want (i.e. ghosted nodes)?!", ngid);

  // global dof id
  int gdof = 0;
  // interface nodes
  if (map_extractor_->INodeRowMap()->MyGID(ngid))
  {
    /* Necessary because the nodes were cloned during the creation of the
     * interface discretization. */
    DRT::Node * inode = map_extractor_->gINode(ngid);
    // access the interface dofset
    gdof = map_extractor_->IDof(inode,ldof);
  }
  // non-interface nodes
  else
  {
    XDisMap::const_iterator cit;
    for (cit=discret_map_.begin();cit!=discret_map_.end();++cit)
    {
      if (cit->second->NodeRowMap()->MyGID(ngid))
      {
        gdof = cit->second->Dof(dofset_id,node,ldof);
        break;
      }
    }
  }

  return gdof;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
std::vector<int> XSTR::MultiDiscretizationWrapper::Dof(
    unsigned nds, const DRT::Node* node) const
{
  std::vector<int> dofs(0);
  Dof(nds,node,dofs);
  return dofs;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XSTR::MultiDiscretizationWrapper::Dof(
    unsigned dofset_id, const DRT::Node* node, std::vector<int>& lm) const
{
  if (dofset_id>0)
    dserror("The MultiDiscretizationWrapper provides only the combined \n"
            "DofRowMap, necessary to create state vectors and system \n"
            "matrices. To access the single DoF-sets, please access the \n"
            "wrapped discretization objects directly.");

  // first check if the node is a considered part of the wrapped discretization
  const int ngid = node->Id();
  if (not NodeRowMap()->MyGID(ngid))
    dserror("The given node (GID: %d) is no part of the MultiDiscretizationWrapper. \n"
            "Maybe you have to call the routine directly on one of the wrapped \n"
            "discretizations to get what you want (i.e. ghosted nodes)?!", ngid);

  // interface nodes
  if (map_extractor_->INodeRowMap()->MyGID(ngid))
  {
    /* Necessary because the nodes were cloned during the creation of the
     * interface discretization. */
    DRT::Node * inode = map_extractor_->gINode(ngid);
    // access the interface dofset
    map_extractor_->IDof(inode,lm);
  }
  // non-interface nodes
  else
  {
    XDisMap::const_iterator cit;
    for (cit=discret_map_.begin();cit!=discret_map_.end();++cit)
    {
      if (cit->second->NodeRowMap()->MyGID(ngid))
      {
        cit->second->Dof(dofset_id,node,lm);
        break;
      }
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XSTR::MultiDiscretizationWrapper::Dof(std::vector<int>& dof,
    const DRT::Node* node ,unsigned dofset_id, unsigned nodaldofset_id,
    const DRT::Element* element) const
{
  if (dofset_id>0)
    dserror("The MultiDiscretizationWrapper provides only the combined \n"
            "DofRowMap, necessary to create state vectors and system \n"
            "matrices. To access the single DoF-sets, please access the \n"
            "wrapped discretization objects directly.");

  // first check if the node is a considered part of the wrapped discretization
  const int ngid = node->Id();
  if ( not NodeRowMap()->MyGID(ngid) )
    dserror("The given node is no part of the MultiDiscretizationWrapper. Maybe \n"
            "you have to call the routine directly on one of the wrapped \n"
            "discretizations to get what you want (i.e. ghosted nodes)?!");

  // interface nodes
  if (map_extractor_->INodeRowMap()->MyGID(ngid))
  {
    /* Necessary because the nodes were cloned during the creation of the
     * interface discretization. */
    DRT::Node * inode = map_extractor_->gINode(ngid);

    // access the interface dofset
    map_extractor_->IDof(dof,inode,nodaldofset_id,element);
  }
  // non-interface nodes
  else
  {
    XDisMap::const_iterator cit;
    for (cit=discret_map_.begin();cit!=discret_map_.end();++cit)
    {
      if (cit->second->NodeRowMap()->MyGID(ngid))
      {
        cit->second->Dof(dof,node,dofset_id,nodaldofset_id,element);
        break;
      }
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XSTR::MultiDiscretizationWrapper::GetCondition(const std::string& name,
    std::vector<DRT::Condition*>& out) const
{
  std::vector<Teuchos::RCP<DRT::Condition> > rcp_out;
  GetCondition(name,rcp_out);
  out.resize(rcp_out.size(),NULL);
  for (unsigned i=0;i<rcp_out.size();++i)
    out[i] = rcp_out[i].get();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XSTR::MultiDiscretizationWrapper::GetCondition(const std::string& name,
    std::vector<Teuchos::RCP<DRT::Condition> >& out) const
{
  XDisMap::const_iterator cit;
  std::vector<Teuchos::RCP<DRT::Condition> > partial_out;
  out.clear();
  for ( cit=discret_map_.begin(); cit!=discret_map_.end(); ++cit )
  {
    cit->second->GetCondition( name, partial_out );

    unsigned start_id = out.size();

    out.resize( out.size() + partial_out.size(), Teuchos::null );

    unsigned count = 0;
    for (unsigned i=start_id;i<out.size();++i)
    {
      out[i] = partial_out[count++];
    }

    if ( count != partial_out.size() )
      dserror( "Mismatch in partial number of conditions found!\n"
          "[count != partial_out.size() <--> %d != %d]",
          count,partial_out.size() );
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
DRT::Condition* XSTR::MultiDiscretizationWrapper::GetCondition(
    const std::string& name) const
{
  XDisMap::const_iterator cit;
  DRT::Condition* cond = NULL;
  for (cit=discret_map_.begin();cit!=discret_map_.end();++cit)
  {
    cond = cit->second->GetCondition(name);

    // skip discretizations which don't hold this condition
    if (cond == NULL)
      continue;

    // skip discretizations with empty conditions
    if (cond->Geometry().size()==0 and cond->Nodes()->size()==0)
      continue;

    break;
  }

  return cond;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XSTR::MultiDiscretizationWrapper::SetState(unsigned nds,
    const std::string& name, Teuchos::RCP<const Epetra_Vector> full_state)
{
  if (nds>0)
    dserror("Currently we support only one dofset in this SetState method,"
        "due to the way how the ExtractVector method is working.");

  XDisMap::const_iterator cit;
  for (cit=discret_map_.begin();cit!=discret_map_.end();++cit)
  {
    Teuchos::RCP<Epetra_Vector> partial =
        map_extractor_->ExtractVector(*full_state,cit->first);
    cit->second->SetState(nds,name,partial.getConst());
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> XSTR::MultiDiscretizationWrapper::GetState(
    unsigned nds, const std::string& name) const
{
  if (nds>0)
    dserror("Currently we support only one dofset in this GetState method,"
        "due to the way how the InsertVector method is working.");

  Teuchos::RCP<Epetra_Vector> merged_state =
      Teuchos::RCP<Epetra_Vector>(new Epetra_Vector(*DofRowMap()));

  XDisMap::const_iterator cit;
  for (cit=discret_map_.begin(); cit!=discret_map_.end(); ++cit)
  {
    Teuchos::RCP<const Epetra_Vector> partial = cit->second->GetState(nds,name);
    dserror("Check if we have to sum the interface DoF's up??? Currently we "
        "would overwrite them.");
    map_extractor_->InsertVector(partial,cit->first,merged_state);
  }
  return merged_state;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XSTR::MultiDiscretizationWrapper::ClearState(bool clearalldofsets)
{
  XDisMap::const_iterator cit;
  for (cit=discret_map_.begin(); cit!=discret_map_.end(); ++cit)
  {
    cit->second->ClearState(clearalldofsets);
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const Teuchos::RCP<LINALG::SparseOperator>& XSTR::MultiDiscretizationWrapper::
    XFieldSysMat( enum XFEM::FieldName field,
    unsigned mat_id )
{
   XSysMatMap::iterator mat_it = partial_system_matrices_.find(field);
   // ------------------------------------------------------------------------
   // insert new matrix, if there is none for the given field
   // ------------------------------------------------------------------------
   if (mat_it==partial_system_matrices_.end())
   {
     std::pair<XSysMatMap::iterator,bool> check = partial_system_matrices_.insert(
         XSysMatPair(field,std::vector<Teuchos::RCP<LINALG::SparseOperator> >(
             (mat_id+1),Teuchos::null)));

     // this should not/cannot happen
     if (not check.second)
       dserror("The new pair couldn't be inserted!");
     // set the iterator
     mat_it = check.first;
   }
   // ------------------------------------------------------------------------
   // check the size and resize if necessary
   // ------------------------------------------------------------------------
   else
   {
     if ( mat_it->second.size() < (mat_id+1) )
       mat_it->second.resize((mat_id+1),Teuchos::null);
   }

   // ------------------------------------------------------------------------
   // get the partial block system matrix
   // ------------------------------------------------------------------------
   // get the partial system matrix pointer
   Teuchos::RCP<LINALG::SparseOperator>& sys_mat = mat_it->second.at(mat_id);
   // create a new partial block system matrix (allocate memory)
   if (sys_mat.is_null())
   {
     sys_mat = Teuchos::rcp<LINALG::SparseOperator>(
         new LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>(
             map_extractor_->SlDofMapExtractor(field),
             map_extractor_->SlDofMapExtractor(field), 81,true,true));
   }
   // Zeroize an already existing partial system matrix
   else
     sys_mat->Zero();

   return sys_mat;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XSTR::MultiDiscretizationWrapper::EvaluateNeumann(
    Teuchos::ParameterList& params,
    Epetra_Vector&          full_systemvector,
    LINALG::SparseOperator* full_systemmatrix)
{
  // temporal partial systemvector
  Teuchos::RCP<Epetra_Vector> psys_vec = Teuchos::null;

  XDisMap::const_iterator cit;
  for (cit=discret_map_.begin();cit!=discret_map_.end();++cit)
  {
    // ------------------------------------------------------------------------
    // evaluate the Neumann condition on each wrapped field discretization
    // ------------------------------------------------------------------------
    psys_vec = Teuchos::rcp(new Epetra_Vector(
        *cit->second->DofRowMap(),true));

    if (full_systemmatrix)
      cit->second->EvaluateNeumann(params,psys_vec,XFieldSysMat(cit->first,0));
    else
      cit->second->EvaluateNeumann(params,psys_vec);

    // ------------------------------------------------------------------------
    // Assemble the partial system vectors into the global full system vector
    // ------------------------------------------------------------------------
    map_extractor_->AddVector(*psys_vec,cit->first,
        full_systemvector,1.0);
    // reset the partial system vector
    psys_vec = Teuchos::null;

    // ------------------------------------------------------------------------
    // Assemble the partial matrices into the global full system matrix
    // ------------------------------------------------------------------------
    if (full_systemmatrix)
      map_extractor_->AddMatrix(*XFieldSysMat(cit->first,0),cit->first,
          *full_systemmatrix,1.0);
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XSTR::MultiDiscretizationWrapper::Evaluate(
    Teuchos::ParameterList&     params,
    Teuchos::RCP<LINALG::SparseOperator> full_systemmatrix1,
    Teuchos::RCP<LINALG::SparseOperator> full_systemmatrix2,
    Teuchos::RCP<Epetra_Vector> full_systemvector1,
    Teuchos::RCP<Epetra_Vector> full_systemvector2,
    Teuchos::RCP<Epetra_Vector> full_systemvector3)
{
  CheckInitSetup();
  XDisMap::const_iterator cit;
  // put input into arrays, such that a loop treatment becomes possible
  Teuchos::RCP<Epetra_Vector> full_systemvectors[3] = {full_systemvector1,
      full_systemvector2,full_systemvector3};
  Teuchos::RCP<LINALG::SparseOperator> full_systemmatrices[2] =
      {full_systemmatrix1, full_systemmatrix2};
  // initialize temporal partial vector pointers
  Teuchos::RCP<Epetra_Vector> psys_vec[3] = {Teuchos::null,Teuchos::null,
      Teuchos::null};
  // initialize temporal partial matrix pointers
  Teuchos::RCP<LINALG::SparseOperator> psys_mat[2] = {Teuchos::null,
      Teuchos::null};

  for (cit=discret_map_.begin();cit!=discret_map_.end();++cit)
  {
    // ------------------------------------------------------------------------
    // create new temporal partial vectors (if necessary)
    // ------------------------------------------------------------------------
    for (unsigned i=0; i<3;++i)
    {
      if (not full_systemvectors[i].is_null())
        psys_vec[i] = Teuchos::rcp<Epetra_Vector>(new Epetra_Vector(
            *cit->second->DofRowMap(0),true));
    }

    // ------------------------------------------------------------------------
    // get the partial system matrices
    // ------------------------------------------------------------------------
    if (not full_systemmatrix1.is_null())
      psys_mat[0] = XFieldSysMat(cit->first,0);
    if (not full_systemmatrix2.is_null())
      psys_mat[1] = XFieldSysMat(cit->first,1);

    // ------------------------------------------------------------------------
    // evaluate elements on each wrapped field discretization
    // ------------------------------------------------------------------------
    cit->second->Evaluate(params,psys_mat[0],psys_mat[1],psys_vec[0],
        psys_vec[1],psys_vec[2]);

    // ------------------------------------------------------------------------
    // Assemble the partial system vectors into the global full system vector
    // ------------------------------------------------------------------------
    for (unsigned i=0; i<3;++i)
    {
      if (not full_systemvectors[i].is_null())
      {
        map_extractor_->AddVector(*psys_vec[i],cit->first,
            *full_systemvectors[i],1.0);
        // reset the temporal partial system vector
        psys_vec[i] = Teuchos::null;
      }
    }

    // ------------------------------------------------------------------------
    // Assemble the partial matrices into the global full system matrices
    // ------------------------------------------------------------------------
    for (unsigned i=0; i<2; ++i)
    {
      if (not full_systemmatrices[i].is_null())
      {
        const Teuchos::RCP<LINALG::SparseOperator> & psysmat_ptr = psys_mat[i];
        if (not psysmat_ptr->Filled())
          psysmat_ptr->Complete();
        map_extractor_->AddMatrix(*psysmat_ptr,cit->first,
            *full_systemmatrices[i],1.0);
      }
    }
  } // end loop over discretization map
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XSTR::MultiDiscretizationWrapper::EvaluateDirichlet(
    Teuchos::ParameterList& params,
    Teuchos::RCP<Epetra_Vector> systemvector,
    Teuchos::RCP<Epetra_Vector> systemvectord,
    Teuchos::RCP<Epetra_Vector> systemvectordd,
    Teuchos::RCP<Epetra_Vector> toggle,
    Teuchos::RCP<LINALG::MapExtractor> dbcmapextractor)
{
  DRT::UTILS::EvaluateDirichlet(*this,params,systemvector,systemvectord,
      systemvectordd,toggle,dbcmapextractor);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XSTR::MultiDiscretizationWrapper::EvaluateCondition(
    Teuchos::ParameterList& params,
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix1,
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix2,
    Teuchos::RCP<Epetra_Vector> systemvector1,
    Teuchos::RCP<Epetra_Vector> systemvector2,
    Teuchos::RCP<Epetra_Vector> systemvector3,
    const std::string& condstring,
    const int condid)
{
  dserror("Currently unsupported!");
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XSTR::MultiDiscretizationWrapper::EvaluateScalars(
    Teuchos::ParameterList& params,
    Teuchos::RCP<Epetra_SerialDenseVector> scalars)
{
  XDisMap::const_iterator cit;
  scalars->Scale(0.0);
  Teuchos::RCP<Epetra_SerialDenseVector> pscalars =
      Teuchos::rcp(new Epetra_SerialDenseVector(scalars->Length()));
  for (cit=discret_map_.begin();cit!=discret_map_.end();++cit)
  {
    // evaluate the scalars on each discretization
    cit->second->EvaluateScalars(params,pscalars);
    *scalars += *pscalars;
    pscalars->Scale(0.0);
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XSTR::MultiDiscretizationWrapper::ComputeNullSpaceIfNecessary(
    Teuchos::ParameterList& solveparams,
    bool recompute)
{
  if (Comm().MyPID() == 0)
  {
    std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
    std::cout << "XSTR::MultiDiscretizationWrapper::ComputeNullSpaceIfNecessary \n"
        "                   Currently unsupported \n";
    std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XSTR::MultiDiscretizationWrapper::CreateNewMap(
    const std::vector<int> & my_entries,
    Teuchos::RCP<Epetra_Map> & new_map ) const
{
  const int * my_entries_ptr = NULL;
  if ( my_entries.size() > 0 )
    my_entries_ptr = & my_entries[0];

  new_map = Teuchos::null;
  new_map = Teuchos::rcp<Epetra_Map>( new Epetra_Map( -1,
      static_cast<int>( my_entries.size() ), my_entries_ptr, 0, Comm() ) );
}
