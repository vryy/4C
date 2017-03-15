/*----------------------------------------------------------------------*/
/*!
\file xfem_condition_manager.cpp

\brief manages the different types of mesh and level-set based coupling conditions and thereby builds the bridge between the
xfluid class and the cut-library

\level 2

<pre>
\maintainer Benedikt Schott
            schott@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>
*/
/*----------------------------------------------------------------------*/

#include "xfem_condition_manager.H"
#include "xfem_utils.H"
#include "Epetra_IntVector.h"

#include "../linalg/linalg_utils.H"

#include "../drt_cut/cut_volumecell.H"

#include "../drt_io/io.H"
#include "../drt_io/io_pstream.H"


//constructor
XFEM::ConditionManager::ConditionManager(
    const std::map<std::string, int> &                dofset_coupling_map, ///< ???
    Teuchos::RCP<DRT::Discretization> &               bg_dis,            ///< background discretization
    std::vector<Teuchos::RCP<DRT::Discretization> > & meshcoupl_dis,     ///< mesh coupling discretizations
    std::vector<Teuchos::RCP<DRT::Discretization> > & levelsetcoupl_dis, ///< levelset coupling discretizations
    const double                                      time,              ///< time
    const int                                         step               ///< time step
  ) :
  dofset_coupling_map_(dofset_coupling_map),
  bg_dis_(bg_dis),
  levelset_gid_(-1),
  time_(time),
  step_(step),
  is_levelset_uptodate_(false),
  ele_lsc_coup_idx_col_(Teuchos::null),
  bg_phinp_(Teuchos::null),
  isinit_(false),
  issetup_(false)
{
  // create Levelset Coupling objects
  {
    std::vector<std::string> conditions_to_check;

    // NOTE: CHANGING THE ORDER HERE CAN INFLUENCE BOOLEAN COMBINATIONS WHEN CREATING UNIQUE LS-FIELD
    conditions_to_check.push_back("XFEMLevelsetNeumann"); // Neumann before Dirichlet, otherwise artificial error at outflow due to rough approximation by level-set
    conditions_to_check.push_back("XFEMLevelsetWeakDirichlet");
    conditions_to_check.push_back("XFEMLevelsetNavierSlip");
    conditions_to_check.push_back("XFEMLevelsetTwophase");
    conditions_to_check.push_back("XFEMLevelsetCombustion");

    CreateCouplings(levelsetcoupl_dis, conditions_to_check, false);
  }

  // create Mesh Coupling objects
  {
    std::vector<std::string> conditions_to_check;
    conditions_to_check.push_back("XFEMSurfNeumann");
    conditions_to_check.push_back("XFEMSurfWeakDirichlet");
    conditions_to_check.push_back("XFEMSurfFSIPart");
    conditions_to_check.push_back("XFEMSurfFSIMono");
    conditions_to_check.push_back("XFEMSurfFPIMono");
    conditions_to_check.push_back("XFEMSurfFluidFluid");
    conditions_to_check.push_back("XFEMSurfNavierSlip");

    CreateCouplings(meshcoupl_dis, conditions_to_check, true);
  }

  SetDofSetCouplingMap( dofset_coupling_map );

}


void XFEM::ConditionManager::CreateCouplings(
    std::vector<Teuchos::RCP<DRT::Discretization> > & coupl_dis,   ///< coupling discretizations
    const std::vector<std::string> & conditions_to_check,          ///< conditions for which coupling objects shall be created
    bool create_mesh_coupling                                      ///< create mesh coupling or level-set coupling object
)
{
  // check if a coupling discretization has relevant conditioned nodes
  // create new coupling object for each type of condition and each coupling discretization
  for(size_t c_idx=0; c_idx<coupl_dis.size(); c_idx++) // loop all specified mesh coupling discretizations
  {
    if(coupl_dis[c_idx] == Teuchos::null) continue;

    std::vector<std::string> names;
    coupl_dis[c_idx]->GetConditionNames( names );

    for(size_t c=0; c<conditions_to_check.size(); c++)
    {
      if(std::find(names.begin(), names.end(), conditions_to_check[c]) == names.end())
        continue;

      // get all conditions of this type, if several conditions with different coupling ids
      std::set<int> coupling_ids;
      GetCouplingIds(*(coupl_dis[c_idx]), conditions_to_check[c], coupling_ids);

      // create new coupling object for each composite
      for(std::set<int>::iterator cid=coupling_ids.begin();
          cid != coupling_ids.end(); ++cid)
      {
        if(create_mesh_coupling)
          CreateNewMeshCoupling(conditions_to_check[c], coupl_dis[c_idx], *cid);
        else
          CreateNewLevelSetCoupling(conditions_to_check[c], coupl_dis[c_idx], *cid );
      }
    }
  }
}

void XFEM::ConditionManager::GetCouplingIds(
    const DRT::Discretization & cond_dis,
    const std::string &         condition_name,
    std::set<int> &             coupling_ids
)
{
  // get all conditions of this type, if several conditions with different coupling ids
  // create an own coupling object for each coupling id
  std::vector<DRT::Condition*> conditions;
  cond_dis.GetCondition( condition_name, conditions);

  // CompositeByCouplingId
  for(size_t s=0; s<conditions.size(); ++s)
  {
    DRT::Condition* cond = conditions[s];
    const int couplingID = cond->GetInt("label");

    coupling_ids.insert(couplingID);
  }
}

void XFEM::ConditionManager::SetDofSetCouplingMap(
    const std::map<std::string, int> & dofset_coupling_map
)
{
  for(int m=0; m<NumMeshCoupling(); m++)
  {
    mesh_coupl_[m]->SetDofSetCouplingMap(dofset_coupling_map_);
  }

  for(int l=0; l<NumLevelSetCoupling(); l++)
  {
    levelset_coupl_[l]->SetDofSetCouplingMap(dofset_coupling_map_);
  }
}



void XFEM::ConditionManager::Status()
{
  int myrank = bg_dis_->Comm().MyPID();

  // -------------------------------------------------------------------
  //                       output to screen
  // -------------------------------------------------------------------
  if (myrank==0)
  {

    printf("   +--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+\n");
    printf("   +----------------------------------------------------XFEM::ConditionManager - Created Coupling objects---------------------------------------------------------------------------------+\n");
    printf("   +----------+-----------+-----------------------------+---------+-----------------------------+-----------------------------+-----------------------------+-----------------------------+\n");
    printf("   | COUP-IDX | START-SID |       CONDITION-TYPE        | COUP-ID |          CUTTER-DIS         | created from CONDITION-DIS  |        COUPLING-DIS         |     AVERAGING-STRATEGY      |\n");

    if(HasMeshCoupling())
    {
      printf("   +----------+-----------+-----------------------------+---------+-----------------------------+-----------------------------+-----------------------------+-----------------------------+\n");
      printf("   |Mesh Coupling Objects                                                                                                                                                                 |\n");
    }

    // loop all mesh coupling objects
    for(int mc=0; mc<(int)mesh_coupl_.size(); mc++)
    {
      mesh_coupl_[mc]->Status(mc, mesh_coupl_start_gid_[mc]);
    }

    if(HasLevelSetCoupling())
    {
      printf("   +----------+-----------+-----------------------------+---------+-----------------------------+-----------------------------+-----------------------------+-----------------------------+\n");
      printf("   |Levelset Coupling Objects                                                                                                                                                             |\n");
    }

    // loop all levelset coupling objects
    for(int lsc=0; lsc<(int)levelset_coupl_.size(); lsc++)
    {
      levelset_coupl_[lsc]->Status(lsc, levelset_gid_);
    }

    printf("   +----------+-----------+-----------------------------+---------+-----------------------------+-----------------------------+-----------------------------+-----------------------------+\n");
    printf("   +--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+\n");
  }
}


void XFEM::ConditionManager::IncrementTimeAndStep(
    const double dt)
{
  step_ += 1;
  time_ += dt;

  // loop all mesh coupling objects
  for(int mc=0; mc<(int)mesh_coupl_.size(); mc++)
  {
    mesh_coupl_[mc]->IncrementTimeAndStep(dt);
  }

  // loop all levelset coupling objects
  for(int lsc=0; lsc<(int)levelset_coupl_.size(); lsc++)
  {
    levelset_coupl_[lsc]->IncrementTimeAndStep(dt);
  }

}


void XFEM::ConditionManager::SetTimeAndStep(
    const double time,
    const int step)
{
  time_ = time;
  step_ = step;

  // loop all mesh coupling objects
  for(int mc=0; mc<(int)mesh_coupl_.size(); mc++)
  {
    mesh_coupl_[mc]->SetTimeAndStep(time, step);
  }

  // loop all levelset coupling objects
  for(int lsc=0; lsc<(int)levelset_coupl_.size(); lsc++)
  {
    levelset_coupl_[lsc]->SetTimeAndStep(time, step);
  }
}



void XFEM::ConditionManager::CreateNewLevelSetCoupling(
    const std::string& cond_name,
    Teuchos::RCP<DRT::Discretization> cond_dis,
    const int coupling_id
)
{
  AddLevelSetCoupling( cond_name, cond_dis, coupling_id );
}


void XFEM::ConditionManager::CreateNewMeshCoupling(
    const std::string& cond_name,
    Teuchos::RCP<DRT::Discretization> cond_dis,         ///< discretization from which the cutter discretization can be derived
    const int coupling_id
)
{
  AddMeshCoupling( cond_name, cond_dis, coupling_id );
}

void XFEM::ConditionManager::Init()
{
  issetup_ = false;

  //--------------------------------------------------------
  // loop all mesh coupling objects
  for(int mc=0; mc<(int)mesh_coupl_.size(); mc++)
  {
    mesh_coupl_[mc]->Init();
  }

  //--------------------------------------------------------
  // loop all levelset coupling objects
  for(int lc=0; lc<(int)levelset_coupl_.size(); lc++)
  {
    levelset_coupl_[lc]->Init();
  }

  isinit_=true;
}


void XFEM::ConditionManager::Setup()
{
  CheckInit();

  // do setup

  //--------------------------------------------------------
  // loop all mesh coupling objects
  for(int mc=0; mc<(int)mesh_coupl_.size(); mc++)
  {
    mesh_coupl_[mc]->Setup();
  }

  //--------------------------------------------------------
  // loop all levelset coupling objects
  for(int lc=0; lc<(int)levelset_coupl_.size(); lc++)
  {
    levelset_coupl_[lc]->Setup();
  }

  Create();

  issetup_ = true;
}



void XFEM::ConditionManager::Create()
{
  numglobal_coupling_sides = 0;
  mesh_coupl_start_gid_.reserve(mesh_coupl_.size());
  levelset_gid_ = -1;

  // set global side Ids for all Mesh coupling discretizations and level-set sides

  //--------------------------------------------------------
  // loop all mesh coupling objects
  for(int mc=0; mc<(int)mesh_coupl_.size(); mc++)
  {
    Teuchos::RCP<DRT::Discretization> mc_cutdis = mesh_coupl_[mc]->GetCutterDis();
    if(mc_cutdis == Teuchos::null) dserror("cutter dis is Teuchos::null");

    // set current number of global coupling sides as start index for global id this coupling object
    mesh_coupl_start_gid_[mc] = numglobal_coupling_sides;

    // increase total number of sides with number of global side elements of this mesh coupling object
    numglobal_coupling_sides += mc_cutdis->NumGlobalElements();
  }

  //--------------------------------------------------------
  // combine the level-set values

  if(levelset_coupl_.size() > 0)
  {
    // add one global levelset side used in the cut library
    levelset_gid_ = numglobal_coupling_sides;
    numglobal_coupling_sides+=1;

    bg_phinp_= LINALG::CreateVector(*bg_dis_->NodeRowMap(), true);

    // information about the coupling condition for level-sets is obtained via the background element
    // for which we store the index of the level-set coupling object
    // we allow for multiple level-set coupling objects however only for one level-set side

    ele_lsc_coup_idx_col_ = Teuchos::rcp(new Epetra_IntVector(*bg_dis_->ElementColMap(),true));
  }

  //--------------------------------------------------------
  // print status of conditionManager to screen
  Status();
}


void XFEM::ConditionManager::SetLevelSetField( const double time )
{
  if(levelset_coupl_.size() != 1)
    dserror("level-set field is not unique, which level-set field to be set?");

  is_levelset_uptodate_ = false;

  // update all level-set fields
  // loop all levelset coupling objects
  for(int lsc=0; lsc<(int)levelset_coupl_.size(); lsc++)
  {
    levelset_coupl_[lsc]->SetLevelSetField(time);
  }
}


void XFEM::ConditionManager::WriteAccess_GeometricQuantities(
    Teuchos::RCP<Epetra_Vector> &      scalaraf,
    Teuchos::RCP<Epetra_MultiVector> & smoothed_gradphiaf,
    Teuchos::RCP<Epetra_Vector> &      curvatureaf
   )
{
  // TOOD: when using two-phase in combination with other levelset, how to access to the right coupling twophase coupling object?
  // TODO: safety check, that there is a unique two-phase coupling object!?
  if(levelset_coupl_.size() != 1)
    dserror("level-set field is not unique, which level-set field to be set by given vector?");

  is_levelset_uptodate_ = false;

  // update the unique level-set field
  Teuchos::rcp_dynamic_cast<XFEM::LevelSetCouplingTwoPhase>(levelset_coupl_[0], true)->WriteAccess_GeometricQuantities(scalaraf, smoothed_gradphiaf, curvatureaf);
}


Teuchos::RCP<const Epetra_Vector> XFEM::ConditionManager::GetLevelSetFieldCol()
{
  if(levelset_coupl_.size()==0) return Teuchos::null;

  // export nodal level-set values to node column map
  Teuchos::RCP<Epetra_Vector> bg_phinp_col = Teuchos::rcp(new Epetra_Vector(*bg_dis_->NodeColMap()));
  LINALG::Export(*GetLevelSetField(),*bg_phinp_col);

  return bg_phinp_col;
}

void XFEM::ConditionManager::UpdateLevelSetField()
{

  //-------------------------------------------------------------------------------------------------
  // Boolean operations like \cap \cup \complementary \ ... are used to combine level-set functions
  //-------------------------------------------------------------------------------------------------
  // NOTE:
  // * we proceed for the coupling objects as they have been created (WDBC, NEUMANN, TWOPHASE, COMBUSTION...)
  // * within one TYPE of conditions (e.g. WDBC) we combine single fields sorted by their coupling ID
  // * the single groups are combined via MAX (\cap) operations,
  //   such that negative value are within fluid and positive values are non-fluid, or the second fluid phase!!!

  // assume same maps between background fluid dis and the cutterdis (scatra dis)

  // note: not only level-set values have to updated, but also which coupling condition is active in which background element
  // -> 1st: store for each node from which level-set coupling the dominating level set values stems from
  // -> 2nd: based on nodal information, we decide which coupling condition has to be evaluated on an element for which the conditions are not unique

  Teuchos::RCP<Epetra_IntVector> node_lsc_coup_idx     = Teuchos::rcp(new Epetra_IntVector(*bg_dis_->NodeRowMap(), true));
  Teuchos::RCP<Epetra_IntVector> node_lsc_coup_idx_col = Teuchos::rcp(new Epetra_IntVector(*bg_dis_->NodeColMap(), true));
  Teuchos::RCP<Epetra_IntVector> ele_lsc_coup_idx      = Teuchos::rcp(new Epetra_IntVector(*bg_dis_->ElementRowMap(), true));

  // for each row node, the dominating levelset coupling index
  for(int lsc=0; lsc<NumLevelSetCoupling(); ++lsc)
  {
    Teuchos::RCP<LevelSetCoupling> & coupling = levelset_coupl_[lsc];

    // get boolean combination w.r.t previously updated combination
    CouplingBase::LevelSetBooleanType ls_boolean_type = coupling->GetBooleanCombination();

    if(lsc == 0)   // initialize with the first coupling!
    {
      if(coupling->GetBooleanCombination() != CouplingBase::ls_none)
        dserror("the first Boundary-Condition level-set coupling (WDBC or NEUMANN) should always use BOOLEANTYPE=none ! Check your boolean operations");

      Teuchos::RCP<Epetra_Vector> tmp = coupling->GetLevelSetFieldAsNodeRowVector();
      const int err = bg_phinp_->Update(1.0, *tmp, 0.0);
      if(err) dserror("update did not work - vectors based on wrong maps?");
    }
    else  // apply boolean combinations for the further level-set fields
    {
      if(ls_boolean_type == CouplingBase::ls_none)
        dserror("there is a level-set coupling for which you did not specify the the BOOLEANTYPE! Check your boolean operations");

      //need to live here
      Teuchos::RCP<Epetra_Vector> tmp = coupling->GetLevelSetFieldAsNodeRowVector();
      CombineLevelSetField(bg_phinp_, tmp, lsc, node_lsc_coup_idx, ls_boolean_type);
    }

    if(coupling->ApplyComplementaryOperator())
      BuildComplementaryLevelSet(bg_phinp_);
  }

  // export to column vector
  LINALG::Export(*node_lsc_coup_idx, *node_lsc_coup_idx_col);

  // set the levelset coupling index for all row elements
  const Epetra_Map* elerowmap  = bg_dis_->ElementRowMap();
  const Epetra_Map* nodecolmap = bg_dis_->NodeColMap();

  // loop all row elements on the processor
  for(int leleid=0;leleid<bg_dis_->NumMyRowElements();++leleid)
  {
    const int gid = elerowmap->GID(leleid);
    DRT::Element* ele = bg_dis_->gElement(gid);
    const int numnode = ele->NumNode();
    const int * nodeids = ele->NodeIds();

    std::set<int> lsc_coupling_indices;

    for(int n=0; n< numnode; ++n)
    {
      int nlid = nodecolmap->LID(nodeids[n]);
      lsc_coupling_indices.insert((*node_lsc_coup_idx_col)[nlid]);
    }

//    if(lsc_coupling_indices.size() > 1)
//    {
//      for(std::set<int>::iterator i= lsc_coupling_indices.begin();
//          i!= lsc_coupling_indices.end();
//          ++i)
//        std::cout << "for element: " << ele->Id() << " following lsc-indices: " << *i << std::endl;
//    }

    // take the one with the lowest coupling index!
    (*ele_lsc_coup_idx)[leleid] = *(lsc_coupling_indices.begin());
  }

  LINALG::Export(*ele_lsc_coup_idx, *ele_lsc_coup_idx_col_);

  is_levelset_uptodate_ = true;
}

void XFEM::ConditionManager::CombineLevelSetField(
    Teuchos::RCP<Epetra_Vector> & vec1,
    Teuchos::RCP<Epetra_Vector> & vec2,
    const int                     lsc_index_2,
    Teuchos::RCP<Epetra_IntVector> &        node_lsc_coup_idx,
    XFEM::CouplingBase::LevelSetBooleanType ls_boolean_type
)
{
  switch(ls_boolean_type)
  {
  case XFEM::CouplingBase::ls_cut:
    SetMaximum(vec1, vec2, lsc_index_2, node_lsc_coup_idx);
    break;
  case XFEM::CouplingBase::ls_union:
    SetMinimum(vec1, vec2, lsc_index_2, node_lsc_coup_idx);
    break;
  case XFEM::CouplingBase::ls_difference:
    SetDifference(vec1, vec2, lsc_index_2, node_lsc_coup_idx);
    break;
  case XFEM::CouplingBase::ls_sym_difference:
    SetSymmetricDifference(vec1, vec2, lsc_index_2, node_lsc_coup_idx);
    break;
  default:
    dserror("unsupported type of boolean operation between two level-sets");
    break;
  }

}


void XFEM::ConditionManager::CheckForEqualMaps(
    const Teuchos::RCP<Epetra_Vector> & vec1,
    const Teuchos::RCP<Epetra_Vector> & vec2
)
{
  if(not vec1->Map().PointSameAs(vec2->Map()))
    dserror("maps do not match!");
}


void XFEM::ConditionManager::SetMinimum(
    Teuchos::RCP<Epetra_Vector> & vec1,
    Teuchos::RCP<Epetra_Vector> & vec2,
    const int                     lsc_index_2,
    Teuchos::RCP<Epetra_IntVector> &        node_lsc_coup_idx
)
{
  int err=-1;

  CheckForEqualMaps(vec1, vec2);

  // loop all nodes on the processor
  for(int lnodeid=0;lnodeid<bg_dis_->NumMyRowNodes();lnodeid++)
  {
    double val1 = (*vec1)[lnodeid];
    double val2 = (*vec2)[lnodeid];


    // std::min(val1, val2);
    int arg = -1;
    double final_val = XFEM::argmin(val1, val2, arg);

    if(arg == 2) (*node_lsc_coup_idx)[lnodeid] = lsc_index_2; // else keep the old lsc coupling

    // now copy the values
    err = vec1->ReplaceMyValue(lnodeid,0,final_val);
    if (err != 0) dserror("error while inserting value into phinp_");
  }
}


void XFEM::ConditionManager::SetMaximum(
    Teuchos::RCP<Epetra_Vector> & vec1,
    Teuchos::RCP<Epetra_Vector> & vec2,
    const int                     lsc_index_2,
    Teuchos::RCP<Epetra_IntVector> &        node_lsc_coup_idx
)
{
  int err=-1;

  CheckForEqualMaps(vec1, vec2);

  // loop all nodes on the processor
  for(int lnodeid=0;lnodeid<bg_dis_->NumMyRowNodes();lnodeid++)
  {
    double val1 = (*vec1)[lnodeid];
    double val2 = (*vec2)[lnodeid];

    // std::max(val1, val2);
    int arg = -1;
    double final_val = XFEM::argmax(val1, val2, arg);

    if(arg == 2) (*node_lsc_coup_idx)[lnodeid] = lsc_index_2; // else keep the old lsc coupling

    // now copy the values
    err = vec1->ReplaceMyValue(lnodeid,0,final_val);
    if (err != 0) dserror("error while inserting value into phinp_");
  }
}


void XFEM::ConditionManager::SetDifference(
    Teuchos::RCP<Epetra_Vector> & vec1,
    Teuchos::RCP<Epetra_Vector> & vec2,
    const int                     lsc_index_2,
    Teuchos::RCP<Epetra_IntVector> &        node_lsc_coup_idx
)
{
  int err=-1;

  CheckForEqualMaps(vec1, vec2);

  // loop all nodes on the processor
  for(int lnodeid=0;lnodeid<bg_dis_->NumMyRowNodes();lnodeid++)
  {
    double val1 = (*vec1)[lnodeid];
    double val2 = (*vec2)[lnodeid];

    // std::max(val1, -val2);
    int arg = -1;
    double final_val = XFEM::argmax(val1, -val2, arg);

    if(arg == 2) (*node_lsc_coup_idx)[lnodeid] = lsc_index_2; // else keep the old lsc coupling

    // now copy the values
    err = vec1->ReplaceMyValue(lnodeid,0,final_val);
    if (err != 0) dserror("error while inserting value into phinp_");
  }
}

void XFEM::ConditionManager::SetSymmetricDifference(
    Teuchos::RCP<Epetra_Vector> & vec1,
    Teuchos::RCP<Epetra_Vector> & vec2,
    const int                     lsc_index_2,
    Teuchos::RCP<Epetra_IntVector> &        node_lsc_coup_idx
)
{
  int err=-1;

  CheckForEqualMaps(vec1, vec2);

  // loop all nodes on the processor
  for(int lnodeid=0;lnodeid<bg_dis_->NumMyRowNodes();lnodeid++)
  {
    double val1 = (*vec1)[lnodeid];
    double val2 = (*vec2)[lnodeid];

    int arg_tmp1 = -1;
    int arg_tmp2 = -1;
    double val_tmp1 = XFEM::argmax( val1, -val2, arg_tmp1);
    double val_tmp2 = XFEM::argmax(-val1,  val2, arg_tmp2);

    int arg_tmp3 = -1;
    double final_val = XFEM::argmin(val_tmp1, val_tmp2, arg_tmp3);

    if(arg_tmp3 == 2)
      if(arg_tmp2 == 2)
        (*node_lsc_coup_idx)[lnodeid] = lsc_index_2; // else keep the old lsc coupling

    if(arg_tmp3 == 1)
      if(arg_tmp1 == 2)
        (*node_lsc_coup_idx)[lnodeid] = lsc_index_2; // else keep the old lsc coupling


    // now copy the values
    err = vec1->ReplaceMyValue(lnodeid,0,final_val);
    if (err != 0) dserror("error while inserting value into phinp_");
  }
}


void XFEM::ConditionManager::BuildComplementaryLevelSet( Teuchos::RCP<Epetra_Vector> & vec1 )
{
  vec1->Scale(-1.0);
}


void XFEM::ConditionManager::ClearState()
{
  // loop all mesh coupling objects
  for(int mc=0; mc<(int)mesh_coupl_.size(); mc++)
  {
    mesh_coupl_[mc]->ClearState();
  }
}

void XFEM::ConditionManager::SetState()
{
  // loop all mesh coupling objects
  for(int mc=0; mc<(int)mesh_coupl_.size(); mc++)
  {
    mesh_coupl_[mc]->SetState();
  }
}

void XFEM::ConditionManager::SetStateDisplacement()
{
  // loop all mesh coupling objects
  for(int mc=0; mc<(int)mesh_coupl_.size(); mc++)
  {
    mesh_coupl_[mc]->SetStateDisplacement();
  }
}

void XFEM::ConditionManager::UpdateStateVectors()
{
  // loop all mesh coupling objects
  for(int mc=0; mc<(int)mesh_coupl_.size(); mc++)
  {
    mesh_coupl_[mc]->UpdateStateVectors();
  }
}

void XFEM::ConditionManager::CompleteStateVectors()
{
  // loop all mesh coupling objects
  for(int mc=0; mc<(int)mesh_coupl_.size(); mc++)
  {
    mesh_coupl_[mc]->CompleteStateVectors();
  }
}

void XFEM::ConditionManager::ZeroStateVectors_FSI()
{
  // loop all mesh coupling objects
  for(int mc=0; mc<(int)mesh_coupl_.size(); mc++)
  {
    mesh_coupl_[mc]->ZeroStateVectors_FSI();
  }
}

void XFEM::ConditionManager::ExportGeometricQuantities()
{
  // loop all mesh coupling objects
  for(int lsc=0; lsc<(int)levelset_coupl_.size(); lsc++)
  {
    levelset_coupl_[lsc]->ExportGeometricQuantities();
  }
}

void XFEM::ConditionManager::GmshOutput(
    const std::string & filename_base,
    const int step,
    const int gmsh_step_diff,
    const bool gmsh_debug_out_screen
)
{
  // loop all mesh coupling objects
  for(int mc=0; mc<(int)mesh_coupl_.size(); mc++)
  {
    mesh_coupl_[mc]->GmshOutput(filename_base, step, gmsh_step_diff, gmsh_debug_out_screen);
  }
  // loop all level set coupling objects
  for(int lc=0; lc<(int)levelset_coupl_.size(); lc++)
  {
    levelset_coupl_[lc]->GmshOutput(filename_base, step, gmsh_step_diff, gmsh_debug_out_screen);
  }
}

void XFEM::ConditionManager::GmshOutputDiscretization(
  std::ostream& gmshfilecontent
)
{
  // loop all mesh coupling objects
  for(int mc=0; mc<(int)mesh_coupl_.size(); mc++)
  {
    mesh_coupl_[mc]->GmshOutputDiscretization(gmshfilecontent);
  }
}

void XFEM::ConditionManager::Output(
    const int step,
    const double time,
    const bool write_restart_data
)
{
  // loop all mesh coupling objects
  for(int mc=0; mc<(int)mesh_coupl_.size(); mc++)
  {
    mesh_coupl_[mc]->Output(step, time, write_restart_data);
  }

  // output for combined levelset field
  // no restart as bg_phinp can be rebuild from single level-set fields
  if(levelset_coupl_.size() > 0)
  {
    Teuchos::RCP<IO::DiscretizationWriter> output = bg_dis_->Writer();
    output->WriteVector("fluid_levelset_boundary", bg_phinp_);
  }

  // loop all level set coupling objects
  for(int lc=0; lc<(int)levelset_coupl_.size(); lc++)
  {
    levelset_coupl_[lc]->Output(step, time, write_restart_data, lc);
  }
}

void XFEM::ConditionManager::LiftDrag(
    const int step,
    const double time
)
{
  // loop all mesh coupling objects
  for(int mc=0; mc<(int)mesh_coupl_.size(); mc++)
  {
    mesh_coupl_[mc]->LiftDrag(step, time);
  }
}

void XFEM::ConditionManager::ReadRestart(
    const int step
)
{
  // loop all mesh coupling objects
  for(int mc=0; mc<(int)mesh_coupl_.size(); mc++)
  {
    mesh_coupl_[mc]->ReadRestart(step);
  }

  // loop all levelset coupling objects
  for(int lsc=0; lsc<(int)levelset_coupl_.size(); lsc++)
  {
    levelset_coupl_[lsc]->ReadRestart(step, lsc);
  }
}


void XFEM::ConditionManager::PrepareSolve()
{
  // loop all mesh coupling objects
  for(int mc=0; mc<(int)mesh_coupl_.size(); mc++)
  {
    mesh_coupl_[mc]->PrepareSolve();
  }

  // loop all levelset coupling objects
  for(int lsc=0; lsc<(int)levelset_coupl_.size(); lsc++)
  {
    levelset_coupl_[lsc]->PrepareSolve();
  }

  is_levelset_uptodate_ = false;
}

bool XFEM::ConditionManager::HasMovingInterface()
{
  // loop all mesh coupling objects
  for(int mc=0; mc<(int)mesh_coupl_.size(); mc++)
  {
    if(mesh_coupl_[mc]->HasMovingInterface()) return true;
  }

  // loop all levelset coupling objects
  for(int lsc=0; lsc<(int)levelset_coupl_.size(); lsc++)
  {
    if(levelset_coupl_[lsc]->HasMovingInterface()) return true;
  }

  return false;
}

bool XFEM::ConditionManager::HasAveragingStrategy(
  INPAR::XFEM::AveragingStrategy strategy)
{
  if(HasLevelSetCoupling())
  {
    for (size_t il = 0; il < levelset_coupl_.size(); ++il)
    {
      if (levelset_coupl_[il]->GetAveragingStrategy() == strategy)
        return true;
    }
  }

  if(HasMeshCoupling())
  {
    for (size_t im = 0; im < mesh_coupl_.size(); ++im)
    {
      if (mesh_coupl_[im]->GetAveragingStrategy() == strategy)
        return true;
    }
  }

  return false;
}

void XFEM::ConditionManager::GetVolumeCellMaterial(
  DRT::Element* actele,
  Teuchos::RCP<MAT::Material> & mat,
  const GEO::CUT::VolumeCell* vc
)
{
  XFEM::UTILS::GetVolumeCellMaterial(actele,mat,vc->Position());
}

void XFEM::ConditionManager::GetInterfaceMasterMaterial(
  DRT::Element* actele,
  Teuchos::RCP<MAT::Material> & mat,
  const GEO::CUT::VolumeCell* vc
)
{
  XFEM::UTILS::GetVolumeCellMaterial(actele,mat,vc->Position());
}

void XFEM::ConditionManager::GetInterfaceSlaveMaterial(
  DRT::Element* actele,
  Teuchos::RCP<MAT::Material> & mat,
  int coup_sid
)
{
  if(IsMeshCoupling(coup_sid))
  {
    int mc = GetMeshCouplingIndex(coup_sid);
    mesh_coupl_[mc]->GetInterfaceSlaveMaterial(actele,mat);
  }
  else if(IsLevelSetCoupling(coup_sid))
  {
    int lc = GetLevelSetCouplingIndex(actele->Id());
    levelset_coupl_[lc]->GetInterfaceSlaveMaterial(actele,mat);
  }
  else
    dserror("The coupling-side id: %d does not correspond to a mesh or levelset coupling object.", coup_sid);
}

//Get Boundary Cell Clone Information <clone_coup_idx, clone_coup_sid>
std::vector<std::pair<int,int> > XFEM::ConditionManager::GetBCCloneInformation(const int coup_sid, const int back_eid,int coup_idx)
{
  std::vector<std::pair<int,int> > BCCloneInformationvector;
  if (coup_idx == -1)
    coup_idx = GetCouplingIndex(coup_sid,back_eid);
  Teuchos::RCP<CouplingBase> coupling = GetCouplingByIdx(coup_idx);
  if (coupling != Teuchos::null)
  {
    //if there are other reasons than mcfpi ... feel free to add your case just here
    Teuchos::RCP<MeshCouplingFPI> mcfpicoupling = Teuchos::rcp_dynamic_cast<MeshCouplingFPI>(coupling);
    if (mcfpicoupling != Teuchos::null)
    {
      if (mcfpicoupling->CutGeometry()) //if this is the ps_ps_block which was loaded into the CUT --> clone from this
      {
        BCCloneInformationvector.push_back(std::pair<int,int>(coup_idx+1,coup_sid + mesh_coupl_start_gid_[coup_idx+1]-mesh_coupl_start_gid_[coup_idx]));
        BCCloneInformationvector.push_back(std::pair<int,int>(coup_idx+2,coup_sid + mesh_coupl_start_gid_[coup_idx+2]-mesh_coupl_start_gid_[coup_idx]));
        BCCloneInformationvector.push_back(std::pair<int,int>(coup_idx+3,coup_sid + mesh_coupl_start_gid_[coup_idx+3]-mesh_coupl_start_gid_[coup_idx]));
        return BCCloneInformationvector;
      }
      else
        dserror("GetBCCloneInformation: Try to clone from FPI MC != PoroStructure?");
    }
    else
      return BCCloneInformationvector;
  }
  else
    dserror("GetBCCloneInformation: Coupling is empty!");
  return BCCloneInformationvector;
}

DRT::Element* XFEM::ConditionManager::GetCouplingElement(
    const int coup_sid, ///< the overall global coupling side id
    DRT::Element * ele
)
{
  if (IsMeshCoupling(coup_sid))
  {
    // get the mesh coupling object index
    const int mc_idx =  GetMeshCouplingIndex(coup_sid);

    // compute the side id w.r.t the cutter discretization the side belongs to
    const int cutterdis_sid = GetCutterDisEleId(coup_sid, mc_idx);

    // get the boundary discretization, the side belongs to
    return mesh_coupl_[mc_idx]->GetCouplingElement(cutterdis_sid);
  }
  else if(IsLevelSetCoupling(coup_sid))
  {
    // get the level-set coupling object index for given background element
    const int lsc_idx = GetLevelSetCouplingIndex(coup_sid);

    // coupling of element with the element itself!
    const int coupldis_eid = ele->Id();

    return levelset_coupl_[lsc_idx]->GetCouplingElement(coupldis_eid);
  }
  else dserror("there is no valid mesh-/levelset-coupling condition object for side: %i", coup_sid);


  return NULL;
}

void XFEM::ConditionManager::GetCouplingEleLocationVector(
    const int coup_sid,
    std::vector<int> & patchlm)
{
  int mc = GetMeshCouplingIndex(coup_sid);
  int sid = GetCutterDisEleId(coup_sid,mc);

  mesh_coupl_[mc]->GetCouplingEleLocationVector(sid, patchlm);
  return;
}

//Get the average weights from the coupling objects
//comment: as soon as we start doing mesh and levelset coupling with overlapping interfaces, we
// need to provide the position of the volumecells down here as well to choose the correct material.
// atm the position this is hardcoded in the coupling objects
// (all assume GEO::CUT::Point::outside, exept TwoPhaseFlow Master outside, Slave inside)
void XFEM::ConditionManager::GetAverageWeights(
    const int coup_sid,                        ///< the overall global coupling side id
    DRT::Element * xfele,                      ///< xfluid ele
    double& kappa_m,                           ///< Weight parameter (parameter +/master side)
    double& kappa_s,                           ///< Weight parameter (parameter -/slave  side)
    bool   & non_xfluid_coupling
    )
{
  DRT::Element* coup_ele  = GetCouplingElement(coup_sid,xfele);
  const int coup_idx = GetCouplingIndex(coup_sid, xfele->Id());

  GetCouplingByIdx(coup_idx)->GetAverageWeights(xfele,coup_ele,kappa_m,kappa_s,non_xfluid_coupling);

  return;
}

/*--------------------------------------------------------------------------------
 * compute viscous part of Nitsche's penalty term scaling for Nitsche's method
 *--------------------------------------------------------------------------------*/
void XFEM::ConditionManager::Get_ViscPenalty_Stabfac(
    const int coup_sid,                                  ///< the overall global coupling side id
    DRT::Element * xfele,                                ///< xfluid ele
    const double& kappa_m,                               ///< Weight parameter (parameter +/master side)
    const double& kappa_s,                               ///< Weight parameter (parameter -/slave  side)
    const double& inv_h_k,                               ///< the inverse characteristic element length h_k
    const DRT::ELEMENTS::FluidEleParameterXFEM* params,  ///< parameterlist which specifies interface configuration
    double& NIT_visc_stab_fac                            ///< viscous part of Nitsche's penalty term
    )
{
  DRT::Element* coup_ele  = GetCouplingElement(coup_sid,xfele);
  const int coup_idx = GetCouplingIndex(coup_sid, xfele->Id());

  GetCouplingByIdx(coup_idx)->Get_ViscPenalty_Stabfac(xfele,coup_ele,kappa_m,kappa_s,inv_h_k,params,NIT_visc_stab_fac);

  return;
}

/*--------------------------------------------------------------------------------
* get the estimation of the penalty scaling in Nitsche's method from the trace inequality for a specific coupling side
 *--------------------------------------------------------------------------------*/
double XFEM::ConditionManager::Get_TraceEstimate_MaxEigenvalue(
const int coup_sid                                  ///< the overall global coupling side id
)
{
  // get the mesh coupling object index
  const int mc_idx =  GetMeshCouplingIndex(coup_sid);
  // compute the side id w.r.t the cutter discretization the side belongs to
  const int cutterdis_sid = GetCutterDisEleId(coup_sid, mc_idx);
  // get the boundary discretization, the side belongs to
  Teuchos::RCP<MeshVolCoupling> mvolcoupling = Teuchos::rcp_dynamic_cast<MeshVolCoupling>(mesh_coupl_[mc_idx]);
  if (mvolcoupling == Teuchos::null) dserror("Cast to MeshVolCoupling failed!");

  return mvolcoupling->Get_EstimateNitscheTraceMaxEigenvalue(mvolcoupling->GetSide(cutterdis_sid));
}
