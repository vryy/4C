/*!
\file xfem_condition_manager.cpp

\brief manages the different types of mesh and level-set based coupling conditions and thereby builds the bridge between the
xfluid class and the cut-library

<pre>
Maintainer: Benedikt Schott
            schott@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>
*/

#include "xfem_condition_manager.H"
#include "xfem_utils.H"

#include "../linalg/linalg_utils.H"

#include "../drt_cut/cut_volumecell.H"

//constructor
XFEM::ConditionManager::ConditionManager(
      Teuchos::RCP<DRT::Discretization> &               bg_dis,          ///< background discretization
      std::vector<Teuchos::RCP<DRT::Discretization> > & meshcoupl_dis,   ///< mesh coupling discretizations
      const double                                      time,            ///< time
      const int                                         step             ///< time step
  ) :
  bg_dis_(bg_dis),
  levelset_gid_(-1),
  time_(time),
  step_(step),
  is_levelset_uptodate_(false),
  bg_phinp_(Teuchos::null)
{

  // create Levelset Coupling objects
  {
    std::vector<std::string> conditions_to_check;

    conditions_to_check.push_back("XFEMLevelsetWeakDirichlet");
    conditions_to_check.push_back("XFEMLevelsetNeumann");
    conditions_to_check.push_back("XFEMLevelsetTwophase");
    conditions_to_check.push_back("XFEMLevelsetCombustion");

    std::vector<std::string> names;
    bg_dis_->GetConditionNames( names );

    // check if background discretization has relevant conditioned nodes
    // create new coupling object for each type of condition
    for(size_t c=0; c<conditions_to_check.size(); c++)
    {
      if(std::find(names.begin(), names.end(), conditions_to_check[c]) != names.end())
        CreateNewLevelSetCoupling(conditions_to_check[c]);

    }
  }

  // create Mesh Coupling objects
  {
    std::vector<std::string> conditions_to_check;
    conditions_to_check.push_back("XFEMSurfFSIPart");
    conditions_to_check.push_back("XFEMSurfFSIMono");
    conditions_to_check.push_back("XFEMSurfCrackFSIPart");
    conditions_to_check.push_back("XFEMSurfFluidFluid");
    conditions_to_check.push_back("XFEMSurfWeakDirichlet");
    conditions_to_check.push_back("XFEMSurfNeumann");


    // check if a coupling discretization has relevant conditioned nodes
    // create new coupling object for each type of condition and each coupling discretization
    for(size_t mc_idx=0; mc_idx<meshcoupl_dis.size(); mc_idx++) // loop all specified mesh coupling discretizations
    {
      if(meshcoupl_dis[mc_idx] == Teuchos::null) continue;

      std::vector<std::string> names;
      meshcoupl_dis[mc_idx]->GetConditionNames( names );

      for(size_t c=0; c<conditions_to_check.size(); c++)
      {
        if(std::find(names.begin(), names.end(), conditions_to_check[c]) != names.end())
          CreateNewMeshCoupling(conditions_to_check[c], meshcoupl_dis[mc_idx]);
      }
    }
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

    printf("   +----------------------------------------------------------------------------------------------------------------------------------------------------------------------------+\n");
    printf("   +----------------------------------------------------XFEM::ConditionManager - Created Coupling objects-----------------------------------------------------------------------+\n");
    printf("   +----------+-----------+-----------------------------+-----------------------------+-----------------------------+-----------------------------+-----------------------------+\n");
    printf("   | COUP-IDX | START-SID |       CONDITION-TYPE        |          CUTTER-DIS         | created from CONDITION-DIS  |        COUPLING-DIS         |     AVERAGING-STRATEGY      |\n");

    if(HasMeshCoupling())
    {
      printf("   +----------+-----------+-----------------------------+-----------------------------+-----------------------------+-----------------------------+-----------------------------+\n");
      printf("   |Mesh Coupling Objects                                                                                                                                                       |\n");
    }

    // loop all mesh coupling objects
    for(int mc=0; mc<(int)mesh_coupl_.size(); mc++)
    {
      mesh_coupl_[mc]->Status(mc, mesh_coupl_start_gid_[mc]);
    }

    if(HasLevelSetCoupling())
    {
      printf("   +----------+-----------+-----------------------------+-----------------------------+-----------------------------+-----------------------------+-----------------------------+\n");
      printf("   |Levelset Coupling Objects                                                                                                                                                   |\n");
    }

    // loop all levelset coupling objects
    for(int lsc=0; lsc<(int)levelset_coupl_.size(); lsc++)
    {
      levelset_coupl_[lsc]->Status(lsc, levelset_gid_);
    }

    printf("   +----------+-----------+-----------------------------+-----------------------------+-----------------------------+-----------------------------+-----------------------------+\n");
    printf("   +----------------------------------------------------------------------------------------------------------------------------------------------------------------------------+\n");
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



void XFEM::ConditionManager::CreateNewLevelSetCoupling(const std::string& cond_name)
{
  AddLevelSetCoupling( cond_name );
}


void XFEM::ConditionManager::CreateNewMeshCoupling(
    const std::string& cond_name,
    Teuchos::RCP<DRT::Discretization> cond_dis         ///< discretization from which the cutter discretization can be derived
)
{
  AddMeshCoupling( cond_name, cond_dis );
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

  // TODO: unify the level-set coupling objects to one unique level-set field
  //--------------------------------------------------------
  // combine the level-set values

  if(levelset_coupl_.size() > 0)
  {

    // add one global levelset side used in the cut library
    levelset_gid_ = numglobal_coupling_sides;
    numglobal_coupling_sides+=1;

    bg_phinp_= LINALG::CreateVector(*bg_dis_->NodeRowMap(), true);

    //TODO: note: information about the coupling condition for level-sets is obtained via the background element
    // for which we store the index of the level-set coupling object
    // we allow for multiple level-set coupling objects however only for one level-set side
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

  // update the unique level-set field
  levelset_coupl_[0]->SetLevelSetField(time);
}


void XFEM::ConditionManager::SetLevelSetField(
   Teuchos::RCP<const Epetra_Vector> scalaraf,
   Teuchos::RCP<const Epetra_Vector> curvatureaf,
   Teuchos::RCP<Epetra_MultiVector>  smoothed_gradphiaf,
   Teuchos::RCP<DRT::Discretization> scatradis
   )
{
  if(levelset_coupl_.size() != 1)
    dserror("level-set field is not unique, which level-set field to be set by given vector?");

  is_levelset_uptodate_ = false;

  // update the unique level-set field
  Teuchos::rcp_dynamic_cast<XFEM::LevelSetCouplingTwoPhase>(levelset_coupl_[0])->SetLevelSetField(scalaraf, curvatureaf, smoothed_gradphiaf, scatradis);
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

  if(levelset_coupl_.size() != 1)
      dserror("level-set field is not unique, update of the global bg_phinp by more than one phinp not implemented yet");

  // assume same maps between background fluid dis and the cutterdis (scatra dis)
  bg_phinp_->Update(1.0, *(levelset_coupl_[0]->GetLevelSetField()), 0.0);

  is_levelset_uptodate_ = true;
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

  // loop all level set coupling objects
  for(int lc=0; lc<(int)levelset_coupl_.size(); lc++)
  {
    levelset_coupl_[lc]->Output(step, time, write_restart_data);
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
    levelset_coupl_[lsc]->ReadRestart(step);
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
  Teuchos::RCP<MAT::Material> & mat
)
{
  XFEM::UTILS::GetVolumeCellMaterial(actele,mat);
}

void XFEM::ConditionManager::GetInterfaceSlaveMaterial(
  DRT::Element* actele,
  Teuchos::RCP<MAT::Material> & mat,
  int coup_sid
)
{
  int mc = GetMeshCouplingIndex(coup_sid);
  mesh_coupl_[mc]->GetInterfaceSlaveMaterial(actele,mat);
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
  else
  {
    if (GetCouplingCondition(coup_sid, ele->Id()).first == INPAR::XFEM::CouplingCond_LEVELSET_TWOPHASE)
      return ele;
  }

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
