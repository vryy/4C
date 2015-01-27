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



#include <Teuchos_TimeMonitor.hpp>

#include "xfem_condition_manager.H"

#include "../drt_inpar/inpar_xfem.H"

#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_globalproblem.H"

#include "../drt_lib/drt_condition_utils.H"
#include "../drt_lib/drt_discret_faces.H"
#include "../drt_lib/drt_discret_xfem.H"
#include "../drt_lib/drt_dofset_transparent_independent.H"
#include "../drt_lib/drt_globalproblem.H"

#include "../linalg/linalg_utils.H"

#include "../drt_inpar/inpar_xfem.H"
#include "../drt_inpar/inpar_fluid.H"

#include "../drt_crack/crackUtils.H"


void XFEM::CouplingBase::SetElementConditions()
{
  // number of column cutter boundary elements
  int nummycolele = cutter_dis_->NumMyColElements();

  cutterele_conds_.clear();
  cutterele_conds_.reserve(nummycolele);

  // initialize the vector invalid coupling-condition type "NONE"
  EleCoupCond init_pair = EleCoupCond(INPAR::XFEM::CouplingCond_NONE,NULL);
  for(int lid=0; lid<nummycolele; lid++) cutterele_conds_.push_back(init_pair);

  //-----------------------------------------------------------------------------------
  // loop all column cutting elements on this processor
  for(int lid=0; lid<nummycolele; lid++)
  {
    DRT::Element* cutele = cutter_dis_->lColElement(lid);

    // loop all possible XFEM-coupling conditions
    for(size_t cond=0; cond < conditions_to_copy_.size(); cond++)
    {
      // get all conditions with given condition name
      std::vector<DRT::Condition*> mycond;
      DRT::UTILS::FindElementConditions(cutele, conditions_to_copy_[cond], mycond);

      // safety checks
      if (mycond.size()>1)
      {
//        dserror("%i conditions of the same name for element %i! %s coupling-condition not unique!", mycond.size(), cutele->Id(), conditions_to_copy_[cond].c_str());
      }
      else if(mycond.size() == 0)
      {
        continue; // try the next condition type
      }
      else
      {
//        std::cout << "unique condition found!" << std::endl;
      }

      INPAR::XFEM::EleCouplingCondType cond_type = CondType_stringToEnum(conditions_to_copy_[cond]);

      // non-coupling condition found (e.g. FSI coupling)
      if(cond_type == INPAR::XFEM::CouplingCond_NONE) continue;

      // non-unique conditions for one cutter element
      if( cutterele_conds_[lid].first != INPAR::XFEM::CouplingCond_NONE )
      {
        dserror("There are two different condition types for the same cutter dis element with id %i: 1st %i, 2nd %i. Make the XFEM coupling conditions unique!",
            cutele->Id(), cutterele_conds_[lid].first, cond_type);
      }

      // store the unique condition pointer to the cutting element
      cutterele_conds_[lid] = EleCoupCond(cond_type, mycond[0]);
    }
  }

  //-----------------------------------------------------------------------------------
  // check if all column cutter elements have a valid condition type
  // loop all column cutting elements on this processor
  for(int lid=0; lid<nummycolele; lid++)
  {
    if(cutterele_conds_[lid].first == INPAR::XFEM::CouplingCond_NONE)
      dserror("cutter element with local id %i has no valid coupling-condition", lid);
  }

}

void XFEM::CouplingBase::SetCouplingStrategy()
{
  const INPAR::XFEM::EleCouplingCondType cond_type = CondType_stringToEnum(cond_name_);

  switch(cond_type)
  {
  case INPAR::XFEM::CouplingCond_SURF_FSI_MONO:
  {
    coupling_strategy_ = INPAR::XFEM::Xfluid_Sided_Coupling; //TODO: rename to xfluid (only one for coupling and wdbc/Neumann)
    break;
  }
  case INPAR::XFEM::CouplingCond_SURF_FLUIDFLUID:
  {
    // ask the first cutter element
    const int lid=0;
    const int val = cutterele_conds_[lid].second->GetInt("COUPSTRATEGY");
    coupling_strategy_ = static_cast<INPAR::XFEM::CouplingStrategy>(val);
    break;
  }
  case INPAR::XFEM::CouplingCond_LEVELSET_TWOPHASE:
  case INPAR::XFEM::CouplingCond_LEVELSET_COMBUSTION:
  {
    coupling_strategy_ = INPAR::XFEM::Harmonic;
    break;
  }
  case INPAR::XFEM::CouplingCond_SURF_FSI_PART:
  case INPAR::XFEM::CouplingCond_SURF_CRACK_FSI_PART:
  case INPAR::XFEM::CouplingCond_SURF_WEAK_DIRICHLET: // set this to Teuchos::null when the values are read from the function instead of the ivelnp vector
  case INPAR::XFEM::CouplingCond_SURF_NEUMANN:
  case INPAR::XFEM::CouplingCond_LEVELSET_WEAK_DIRICHLET:
  case INPAR::XFEM::CouplingCond_LEVELSET_NEUMANN:
  {
    coupling_strategy_ = INPAR::XFEM::Xfluid_Sided_weak_DBC; //TODO: rename to xfluid
    break;
  }
  default: dserror("which is the coupling discretization for this type of coupling %i?", cond_type); break;
  }
}


void XFEM::CouplingBase::SetCouplingDiscretization()
{
  const INPAR::XFEM::EleCouplingCondType cond_type = CondType_stringToEnum(cond_name_);

  switch(cond_type)
  {
  case INPAR::XFEM::CouplingCond_SURF_FSI_MONO:
  {
    coupl_dis_ = cutter_dis_; break;
    break;
  }
  case INPAR::XFEM::CouplingCond_SURF_FLUIDFLUID:
  {
    // depending on the weighting strategy
    if(coupling_strategy_==INPAR::XFEM::Xfluid_Sided_Coupling)
    {
      coupl_dis_ = cutter_dis_;
    }
    else if(coupling_strategy_==INPAR::XFEM::Embedded_Sided_Coupling or
        coupling_strategy_==INPAR::XFEM::Two_Sided_Coupling )
    {
      coupl_dis_ = cond_dis_;
    }
    else dserror("invalid coupling strategy for fluid-fluid application");
    break;
  }
  case INPAR::XFEM::CouplingCond_LEVELSET_TWOPHASE:
  case INPAR::XFEM::CouplingCond_LEVELSET_COMBUSTION:
  {
    coupl_dis_ = bg_dis_;
    break;
  }
  case INPAR::XFEM::CouplingCond_SURF_FSI_PART:
  case INPAR::XFEM::CouplingCond_SURF_CRACK_FSI_PART:
  case INPAR::XFEM::CouplingCond_SURF_WEAK_DIRICHLET: // set this to Teuchos::null when the values are read from the function instead of the ivelnp vector
  case INPAR::XFEM::CouplingCond_SURF_NEUMANN:
  {
    coupl_dis_ = cutter_dis_;
    break;
  }
  case INPAR::XFEM::CouplingCond_LEVELSET_WEAK_DIRICHLET:
  case INPAR::XFEM::CouplingCond_LEVELSET_NEUMANN:
  {
    coupl_dis_ = Teuchos::null;
    break;
  }
  default: dserror("which is the coupling discretization for this type of coupling %i?", cond_type); break;
  }
}



XFEM::MeshCoupling::MeshCoupling(
    Teuchos::RCP<DRT::Discretization>&  bg_dis,   ///< background discretization
    const std::string &                 cond_name,///< name of the condition, by which the derived cutter discretization is identified
    Teuchos::RCP<DRT::Discretization>&  cond_dis  ///< discretization from which the cutter discretization is derived
) : CouplingBase(bg_dis, cond_name, cond_dis)
{

  // set list of conditions that will be copied to the new cutter discretization
  SetConditionsToCopy();

  // create a cutter discretization from conditioned nodes of the given coupling discretization
  CreateCutterDisFromCondition();

  std::cout << *cutter_dis_ << std::endl;

  // set unique element conditions
  SetElementConditions();

  // set the coupling strategy
  SetCouplingStrategy();

  // set coupling discretization
  SetCouplingDiscretization();

  // initialize state vectors based on cutter discretization
  InitStateVectors();

}


void XFEM::MeshCoupling::SetConditionsToCopy()
{
  // fill list of conditions that will be copied to the new cutter discretization
  conditions_to_copy_.push_back(cond_name_);

  // additional conditions required for the new boundary conditions
  conditions_to_copy_.push_back("FSICoupling");  // for partitioned and monolithic XFSI
}


/*--------------------------------------------------------------------------*
 | Create the cutter discretization                                       |
 *--------------------------------------------------------------------------*/
void XFEM::MeshCoupling::CreateCutterDisFromCondition()
{
  // create name string for new cutter discretization (e.g, "boundary_of_struct" or "boundary_of_fluid")
  std::string cutterdis_name ("boundary_of_");
  cutterdis_name += cond_dis_->Name();


  //--------------------------------
  // create the new cutter discretization form the conditioned coupling discretization
  cutter_dis_ = DRT::UTILS::CreateDiscretizationFromCondition(
      cond_dis_,                ///< discretization with condition
      cond_name_,               ///< name of the condition, by which the derived discretization is identified
      cutterdis_name,           ///< name of the new discretization
      GetBELEName(cond_dis_),   ///< name/type of the elements to be created
      conditions_to_copy_       ///< list of conditions that will be copied to the new discretization
  );
  //--------------------------------


  if (cutter_dis_->NumGlobalNodes() == 0)
  {
    dserror("Empty cutter discretization detected. No coupling can be performed...");
  }

  // for parallel jobs we have to call TransparentDofSet with additional flag true
  bool parallel = cond_dis_->Comm().NumProc() > 1;
  Teuchos::RCP<DRT::DofSet> newdofset = Teuchos::rcp(new
      DRT::TransparentIndependentDofSet(cond_dis_,parallel));

  cutter_dis_->ReplaceDofSet(newdofset); //do not call this with true!!
  cutter_dis_->FillComplete();

  // create node and element distribution with elements and nodes ghosted on all processors
  const Epetra_Map noderowmap = *cutter_dis_->NodeRowMap();
  const Epetra_Map elemrowmap = *cutter_dis_->ElementRowMap();

  // put all boundary nodes and elements onto all processors
  const Epetra_Map nodecolmap = *LINALG::AllreduceEMap(noderowmap);
  const Epetra_Map elemcolmap = *LINALG::AllreduceEMap(elemrowmap);

  // redistribute nodes and elements to column (ghost) map
  cutter_dis_->ExportColumnNodes(nodecolmap);
  cutter_dis_->ExportColumnElements(elemcolmap);

  cutter_dis_->FillComplete();
}


void XFEM::MeshCoupling::InitStateVectors()
{
  const Epetra_Map* cutterdofrowmap = cutter_dis_->DofRowMap();

  ivelnp_  = LINALG::CreateVector(*cutterdofrowmap,true);
  iveln_   = LINALG::CreateVector(*cutterdofrowmap,true);
  ivelnm_  = LINALG::CreateVector(*cutterdofrowmap,true);

  idispnp_ = LINALG::CreateVector(*cutterdofrowmap,true);
  idispn_  = LINALG::CreateVector(*cutterdofrowmap,true);

  itrueresidual_ = LINALG::CreateVector(*cutterdofrowmap,true);
}

Teuchos::RCP<const Epetra_Vector> XFEM::MeshCoupling::GetCutterDispCol()
{
  // export cut-discretization mesh displacements
  Teuchos::RCP<Epetra_Vector> idispcol = LINALG::CreateVector(*cutter_dis_->DofColMap(),true);
  LINALG::Export(*idispnp_,*idispcol);

  return idispcol;
}


//! constructor
XFEM::MeshCouplingFSICrack::MeshCouplingFSICrack(
    Teuchos::RCP<DRT::Discretization>&  bg_dis,   ///< background discretization
    const std::string &                 cond_name,///< name of the condition, by which the derived cutter discretization is identified
    Teuchos::RCP<DRT::Discretization>&  cond_dis  ///< discretization from which cutter discretization can be derived
) : MeshCoupling(bg_dis,cond_name,cond_dis)
{
  InitCrackInitiationsPoints();
}


void XFEM::MeshCouplingFSICrack::SetCutterDis(Teuchos::RCP<DRT::Discretization> cutter_dis_new)
{
  cutter_dis_ = cutter_dis_new;

  // update the Coupling object

  // set unique element conditions
  SetElementConditions();

  // set the coupling strategy
  SetCouplingStrategy();

  // set coupling discretization
  SetCouplingDiscretization();
}


void XFEM::MeshCouplingFSICrack::InitCrackInitiationsPoints()
{
  tip_nodes_.clear();

  DRT::Condition* crackpts = cond_dis_->GetCondition( "CrackInitiationPoints" );

  const std::vector<int>* crackpt_nodes = const_cast<std::vector<int>* >(crackpts->Nodes());


  for(std::vector<int>::const_iterator inod=crackpt_nodes->begin(); inod!=crackpt_nodes->end();inod++ )
  {
    const int nodid = *inod;
    LINALG::Matrix<3, 1> xnod( true );

    tip_nodes_[nodid] = xnod;
  }

  if( tip_nodes_.size() == 0 )
    dserror("crack initiation points unspecified\n");

/*---------------------- POSSIBILITY 2 --- adding crack tip elements ----------------------------*/
/*
{
DRT::Condition* crackpts = soliddis_->GetCondition( "CrackInitiationPoints" );

const std::vector<int>* tipnodes = const_cast<std::vector<int>* >(crackpts->Nodes());

if( tipnodes->size() == 0 )
  dserror("crack initiation points unspecified\n");

addCrackTipElements( tipnodes );
}*/

/*  Teuchos::RCP<DRT::DofSet> newdofset1 = Teuchos::rcp(new DRT::TransparentIndependentDofSet(soliddis_,true));

boundarydis_->ReplaceDofSet(newdofset1);//do not call this with true!!
boundarydis_->FillComplete();*/
}


void XFEM::MeshCouplingFSICrack::UpdateBoundaryValuesAfterCrack(
    const std::map<int,int>& oldnewIds
)
{
  DRT::CRACK::UTILS::UpdateThisEpetraVectorCrack( cutter_dis_, ivelnp_, oldnewIds );
  DRT::CRACK::UTILS::UpdateThisEpetraVectorCrack( cutter_dis_, iveln_,  oldnewIds );
  DRT::CRACK::UTILS::UpdateThisEpetraVectorCrack( cutter_dis_, ivelnm_, oldnewIds );

  DRT::CRACK::UTILS::UpdateThisEpetraVectorCrack( cutter_dis_, idispnp_, oldnewIds );
  DRT::CRACK::UTILS::UpdateThisEpetraVectorCrack( cutter_dis_, idispn_,  oldnewIds );

  //itrueresidual_ = LINALG::CreateVector(*boundarydis_->DofRowMap(),true);
  DRT::CRACK::UTILS::UpdateThisEpetraVectorCrack( cutter_dis_, itrueresidual_, oldnewIds );

  //TODO: I guess the following lines are unnecessary (Sudhakar)
  {
    //iforcenp_ = LINALG::CreateVector(*boundarydis_->DofRowMap(),true);
    //LINALG::Export( *itrueresidual_, *iforcenp_ );

  }

  //TODO: Check whether the output in case of crack-FSI work properly (Sudhakar)
  //boundarydis_->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(boundarydis_)));
  //boundary_output_ = boundarydis_->Writer();
}


XFEM::LevelSetCoupling::LevelSetCoupling(
    Teuchos::RCP<DRT::Discretization>&  bg_dis,   ///< background discretization
    const std::string &                 cond_name ///< name of the condition, by which the derived cutter discretization is identified
) : CouplingBase(bg_dis, cond_name, bg_dis)
{
  /// level-set field is given w.r.t background mesh
  /// NOTE: more generally it would be possible cutterdis != bg_dis for the single LevelSetCoupling,
  /// however, the unique bg_phinp vector stored in the ConditionManager has to be given w.r.t bgdis
  cutter_dis_ = bg_dis;

  SetConditionsToCopy();

  SetElementConditions();

  // set the coupling strategy
  SetCouplingStrategy();

  // set coupling discretization
  SetCouplingDiscretization();

  // create node-based vector with level-set values
  phinp_ = Teuchos::rcp(new Epetra_Vector(*cutter_dis_->NodeRowMap()));

}


void XFEM::LevelSetCoupling::SetConditionsToCopy()
{
  // set only the unique given condition name
  conditions_to_copy_.push_back(cond_name_);
}


/*----------------------------------------------------------------------*
 | ... |
 *----------------------------------------------------------------------*/
void XFEM::LevelSetCoupling::SetLevelSetField(
   Teuchos::RCP<const Epetra_Vector> scalaraf,
   Teuchos::RCP<DRT::Discretization> scatradis
   )
{
  // initializations
  int err(0);
  double value(0.0);
  std::vector<int> nodedofs;


  // loop all nodes on the processor
  for(int lnodeid=0;lnodeid<cutter_dis_->NumMyRowNodes();lnodeid++)
  {
    // get the processor's local scatra node
    DRT::Node* lscatranode = scatradis->lRowNode(lnodeid);

    // find out the global dof id of the last(!) dof at the scatra node
    const int numscatradof = scatradis->NumDof(0,lscatranode);
    const int globalscatradofid = scatradis->Dof(0,lscatranode,numscatradof-1);
    const int localscatradofid = scalaraf->Map().LID(globalscatradofid);
    if (localscatradofid < 0)
      dserror("localdofid not found in map for given globaldofid");

    // now copy the values
    value = (*scalaraf)[localscatradofid];
    err = phinp_->ReplaceMyValue(lnodeid,0,value);
    if (err != 0) dserror("error while inserting value into phinp_");

  }

  return;
}


/*----------------------------------------------------------------------*
 | ... |
 *----------------------------------------------------------------------*/
void XFEM::LevelSetCoupling::SetLevelSetField(const double time)
{
  // initializations
  int err(0);
  double value(0.0);
  std::vector<int> nodedofs;

  // get the function from the first element
  const int lid=0;
  DRT::Condition* cond = cutterele_conds_[lid].second;
  const int func_no = cond->GetInt("levelsetfieldno");

  // loop all nodes on the processor
  for(int lnodeid=0;lnodeid<cutter_dis_->NumMyRowNodes();lnodeid++)
  {
    // get the processor's local scatra node
    DRT::Node* lnode = cutter_dis_->lRowNode(lnodeid);

    // get value
    value=DRT::Problem::Instance()->Funct(func_no-1).Evaluate(0,lnode->X(),time,NULL);

    // now copy the values
    err = phinp_->ReplaceMyValue(lnodeid,0,value);
    if (err != 0) dserror("error while inserting value into phinp_");
  }

  return;
}



//constructor
XFEM::ConditionManager::ConditionManager(
      Teuchos::RCP<DRT::Discretization> &               bg_dis,          ///< background discretization
      std::vector<Teuchos::RCP<DRT::Discretization> > & meshcoupl_dis    ///< mesh coupling discretizations
  ) :
  bg_dis_(bg_dis),
  levelset_gid_(-1),
  time_(0.0),
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

void XFEM::ConditionManager::CreateNewLevelSetCoupling(const std::string& cond_name)
{ std::cout << "create new LevelSetCoupling: " << cond_name << std::endl;
  AddLevelSetCoupling( cond_name );
}


void XFEM::ConditionManager::CreateNewMeshCoupling(
    const std::string& cond_name,
    Teuchos::RCP<DRT::Discretization> cond_dis         ///< discretization from which the cutter discretization can be derived
)
{std::cout << "create new MeshCoupling: " << cond_name << std::endl;
  AddMeshCoupling( cond_name, cond_dis );
}


void XFEM::ConditionManager::Create(const double time)
{
  time_ = time;
  numglobal_coupling_sides = 0;
  mesh_coupl_start_gid_.reserve(mesh_coupl_.size());
  levelset_gid_ = -1;

  // set global side Ids for all Mesh coupling discretizations and level-set sides

  // loop all mesh coupling objects
  for(int mc=0; mc<(int)mesh_coupl_.size(); mc++)
  {
    Teuchos::RCP<DRT::Discretization> mc_cutdis = mesh_coupl_[mc]->GetCutterDis();
    if(mc_cutdis == Teuchos::null) dserror("cutter dis is Teuchos::null");

    // set current number of global coupling sides as start index for global id this coupling object
    mesh_coupl_start_gid_[mc] = numglobal_coupling_sides;

    std::cout << "mesh coupling object " << mc << " starts with global side index " << mesh_coupl_start_gid_[mc] << std::endl;

    // increase total number of sides with number of global side elements of this mesh coupling object
    numglobal_coupling_sides += mc_cutdis->NumGlobalElements();
  }

  // TODO: unify the level-set coupling objects to one unique level-set field
  // combine the level-set values

  if(levelset_coupl_.size() > 0)
  {

    // add one global levelset side used in the cut library
    levelset_gid_ = numglobal_coupling_sides;
    numglobal_coupling_sides+=1;

    std::cout << "levelset coupling object " << " has global side id index " << levelset_gid_ << std::endl;

    bg_phinp_= LINALG::CreateVector(*bg_dis_->NodeRowMap(), true);

    //TODO: note: information about the coupling condition for level-sets is obtained via the background element
    // for which we store the index of the level-set coupling object
    // we allow for multiple level-set coupling objects however only for one level-set side
  }



  //TODO: Status-routine to print to screen

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
   Teuchos::RCP<DRT::Discretization> scatradis
   )
{
  if(levelset_coupl_.size() != 1)
    dserror("level-set field is not unique, which level-set field to be set by given vector?");

  is_levelset_uptodate_ = false;

  // update the unique level-set field
  levelset_coupl_[0]->SetLevelSetField(scalaraf, scatradis);
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
