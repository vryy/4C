/*----------------------------------------------------------------------*/
/*! \file

\brief manages the different types of mesh and level-set based coupling conditions and thereby
builds the bridge between the xfluid class and the cut-library

\level 2

\maintainer  Christoph Ager
             ager@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289-15249
*/
/*----------------------------------------------------------------------*/


#ifndef XFEM_CONDITION_MANAGER_H_
#define XFEM_CONDITION_MANAGER_H_

// Might want to forward declare the levelset and mesh coupling and put these functions into
// the xfem_condition_manager.cpp file in a later stage
#include "xfem_coupling_levelset.H"
#include "xfem_coupling_mesh.H"
#include "xfem_coupling_mesh_coupled_levelset.H"
#include "xfem_coupling_fpi_mesh.H"

#include <Epetra_IntVector.h>

namespace GEO
{
  namespace CUT
  {
    class VolumeCell;
  }
}  // namespace GEO

namespace DRT
{
  namespace ELEMENTS
  {
    // finally this parameter list should go and all interface relevant parameters should be stored
    // in the condition mangager or coupling objects
    class FluidEleParameterXFEM;
  }  // namespace ELEMENTS
}  // namespace DRT

namespace XFEM
{
  template <typename _Tp>
  inline const _Tp& argmin(const _Tp& __a, const _Tp& __b, int& arg)
  {
    // return __comp(__b, __a) ? __b : __a;
    if (__b < __a)
    {
      arg = 2;
      return __b;
    }
    arg = 1;
    return __a;
  }

  template <typename _Tp>
  inline const _Tp& argmax(const _Tp& __a, const _Tp& __b, int& arg)
  {
    // return __comp(__a, __b) ? __b : __a;
    if (__a < __b)
    {
      arg = 2;
      return __b;
    }
    arg = 1;
    return __a;
  }

  /*!
  \brief Manages the conditions for the xfluid (i.e. levelset/mesh cut and what BC are applied at
  these)
   */
  class ConditionManager
  {
   public:
    //! constructor
    explicit ConditionManager(const std::map<std::string, int>& dofset_coupling_map,  ///< ???
        Teuchos::RCP<DRT::Discretization>& bg_dis,  ///< background discretization
        std::vector<Teuchos::RCP<DRT::Discretization>>&
            meshcoupl_dis,  ///< mesh coupling discretizations
        std::vector<Teuchos::RCP<DRT::Discretization>>&
            levelsetcoupl_dis,  ///< levelset coupling discretizations
        const double time,      ///< time
        const int step          ///< time step
    );

    //! destructor
    ~ConditionManager() { return; }

    // TODO: we need set private and public !!! this however causes moving functions in the header
    // file!

    void SetTimeAndStep(const double time, const int step);

    void GetCouplingIds(const DRT::Discretization& cond_dis, const std::string& condition_name,
        std::set<int>& coupling_ids);

    void SetDofSetCouplingMap(const std::map<std::string, int>& dofset_coupling_map);

    void Status();

    void IncrementTimeAndStep(const double dt);

    void CreateNewLevelSetCoupling(const std::string& cond_name,
        Teuchos::RCP<DRT::Discretization>
            cond_dis,  ///< discretization from which the cutter discretization can be derived
        const int coupling_id);

    void CreateCouplings(
        std::vector<Teuchos::RCP<DRT::Discretization>>& coupl_dis,  ///< coupling discretizations
        const std::vector<std::string>&
            conditions_to_check,   ///< conditions for which coupling objects shall be created
        bool create_mesh_coupling  ///< create mesh coupling or level-set coupling object
    );

    void CreateNewMeshCoupling(const std::string& cond_name,
        Teuchos::RCP<DRT::Discretization>
            cond_dis,  ///< discretization from which the cutter discretization can be derived
        const int coupling_id);

    /// create a new mesh-coupling object based on the given coupling discretization
    void AddMeshCoupling(const std::string& cond_name, Teuchos::RCP<DRT::Discretization> cond_dis,
        const int coupling_id)
    {
      if (CondType_stringToEnum(cond_name) == INPAR::XFEM::CouplingCond_SURF_FSI_PART or
          CondType_stringToEnum(cond_name) == INPAR::XFEM::CouplingCond_SURF_FSI_MONO)
      {
        mesh_coupl_.push_back(Teuchos::rcp(
            new MeshCouplingFSI(bg_dis_, cond_name, cond_dis, coupling_id, time_, step_)));
      }
      else if (CondType_stringToEnum(cond_name) == INPAR::XFEM::CouplingCond_SURF_FPI_MONO)
      {
        mesh_coupl_.push_back(Teuchos::rcp(new MeshCouplingFPI(
            bg_dis_, cond_name, cond_dis, coupling_id, time_, step_, MeshCouplingFPI::ps_ps)));
        mesh_coupl_.push_back(Teuchos::rcp(new MeshCouplingFPI(
            bg_dis_, cond_name, cond_dis, coupling_id, time_, step_, MeshCouplingFPI::ps_pf)));
        mesh_coupl_.push_back(Teuchos::rcp(new MeshCouplingFPI(
            bg_dis_, cond_name, cond_dis, coupling_id, time_, step_, MeshCouplingFPI::pf_ps)));
        mesh_coupl_.push_back(Teuchos::rcp(new MeshCouplingFPI(
            bg_dis_, cond_name, cond_dis, coupling_id, time_, step_, MeshCouplingFPI::pf_pf)));
      }
      else if (CondType_stringToEnum(cond_name) == INPAR::XFEM::CouplingCond_SURF_WEAK_DIRICHLET)
      {
        mesh_coupl_.push_back(Teuchos::rcp(new MeshCouplingWeakDirichlet(
            bg_dis_, cond_name, cond_dis, coupling_id, time_, step_, (bg_dis_ == cond_dis))));
      }
      else if (CondType_stringToEnum(cond_name) == INPAR::XFEM::CouplingCond_SURF_NEUMANN)
      {
        mesh_coupl_.push_back(Teuchos::rcp(new MeshCouplingNeumann(
            bg_dis_, cond_name, cond_dis, coupling_id, time_, step_, (bg_dis_ == cond_dis))));
      }
      else if (CondType_stringToEnum(cond_name) == INPAR::XFEM::CouplingCond_SURF_NAVIER_SLIP)
      {
        mesh_coupl_.push_back(Teuchos::rcp(new MeshCouplingNavierSlip(
            bg_dis_, cond_name, cond_dis, coupling_id, time_, step_, (bg_dis_ == cond_dis))));
      }
      else if (CondType_stringToEnum(cond_name) ==
               INPAR::XFEM::CouplingCond_SURF_NAVIER_SLIP_TWOPHASE)
      {
        mesh_coupl_.push_back(Teuchos::rcp(new MeshCouplingNavierSlipTwoPhase(
            bg_dis_, cond_name, cond_dis, coupling_id, time_, step_, (bg_dis_ == cond_dis))));
      }
      else if (CondType_stringToEnum(cond_name) == INPAR::XFEM::CouplingCond_SURF_FLUIDFLUID)
      {
        mesh_coupl_.push_back(Teuchos::rcp(
            new MeshCouplingFluidFluid(bg_dis_, cond_name, cond_dis, coupling_id, time_, step_)));
      }
      else
      {
        mesh_coupl_.push_back(Teuchos::rcp(new MeshCoupling(
            bg_dis_, cond_name, cond_dis, coupling_id, time_, step_, "", (bg_dis_ == cond_dis))));
      }
    }

    /// add a new level-set-coupling object based on the given coupling discretization
    void AddLevelSetCoupling(const std::string& cond_name,
        Teuchos::RCP<DRT::Discretization>
            cond_dis,  ///< discretization from which the cutter discretization can be derived
        const int coupling_id)
    {
      if (CondType_stringToEnum(cond_name) == INPAR::XFEM::CouplingCond_LEVELSET_WEAK_DIRICHLET)
      {
        levelset_coupl_.push_back(Teuchos::rcp(new LevelSetCouplingWeakDirichlet(
            bg_dis_, cond_name, cond_dis, coupling_id, time_, step_)));
      }
      else if (CondType_stringToEnum(cond_name) == INPAR::XFEM::CouplingCond_LEVELSET_NEUMANN)
      {
        levelset_coupl_.push_back(Teuchos::rcp(
            new LevelSetCouplingNeumann(bg_dis_, cond_name, cond_dis, coupling_id, time_, step_)));
      }
      else if (CondType_stringToEnum(cond_name) == INPAR::XFEM::CouplingCond_LEVELSET_NAVIER_SLIP)
      {
        levelset_coupl_.push_back(Teuchos::rcp(new LevelSetCouplingNavierSlip(
            bg_dis_, cond_name, cond_dis, coupling_id, time_, step_)));
      }
      else if (CondType_stringToEnum(cond_name) == INPAR::XFEM::CouplingCond_LEVELSET_TWOPHASE)
      {
        levelset_coupl_.push_back(Teuchos::rcp(
            new LevelSetCouplingTwoPhase(bg_dis_, cond_name, cond_dis, coupling_id, time_, step_)));
      }
      else if (CondType_stringToEnum(cond_name) == INPAR::XFEM::CouplingCond_LEVELSET_COMBUSTION)
      {
        levelset_coupl_.push_back(Teuchos::rcp(new LevelSetCouplingCombustion(
            bg_dis_, cond_name, cond_dis, coupling_id, time_, step_)));
      }
      else
      {
        levelset_coupl_.push_back(Teuchos::rcp(
            new LevelSetCoupling(bg_dis_, cond_name, cond_dis, coupling_id, time_, step_)));
      }
    }



    /// Getters

    /// get cutter discretization the coupling side belongs to
    Teuchos::RCP<DRT::Discretization> GetCutterDis(
        const int coup_sid  ///< the global id of the coupling side
    )
    {
      if (IsLevelSetCoupling(coup_sid)) return Teuchos::null;

      // get the mesh coupling object index
      int mc_idx = GetMeshCouplingIndex(coup_sid);

      return mesh_coupl_[mc_idx]->GetCutterDis();
    }

    /// get cutter discretization the coupling side belongs to
    Teuchos::RCP<DRT::Discretization> GetCouplingDis(
        const int coup_sid  ///< the global id of the coupling side
    )
    {
      if (IsLevelSetCoupling(coup_sid)) return bg_dis_;

      // get the mesh coupling object index
      int mc_idx = GetMeshCouplingIndex(coup_sid);

      return mesh_coupl_[mc_idx]->GetCouplingDis();
    }

    Teuchos::RCP<MeshCoupling> GetMeshCoupling(const int mc_idx)
    {
      if (mc_idx < (int)mesh_coupl_.size() and mc_idx >= 0) return mesh_coupl_[mc_idx];

      return Teuchos::null;
    }

    Teuchos::RCP<LevelSetCoupling> GetLevelSetCoupling(const int ls_idx)
    {
      if (ls_idx < (int)levelset_coupl_.size() and ls_idx >= 0) return levelset_coupl_[ls_idx];

      return Teuchos::null;
    }

    Teuchos::RCP<CouplingBase> GetCoupling(const std::string& name)
    {
      Teuchos::RCP<CouplingBase> coupling = Teuchos::null;

      coupling = GetMeshCoupling(name);

      if (coupling != Teuchos::null) return coupling;

      coupling = GetLevelSetCoupling(name);

      if (coupling != Teuchos::null) return coupling;

      return Teuchos::null;
    }


    Teuchos::RCP<MeshCoupling> GetMeshCoupling(const std::string& name)
    {
      for (int m = 0; m < NumMeshCoupling(); m++)
      {
        if (mesh_coupl_[m]->GetName() == name) return mesh_coupl_[m];
      }
      return Teuchos::null;
    }

    int GetCouplingIndex(const std::string& name)
    {
      int coup_idx = -1;

      for (int m = 0; m < NumMeshCoupling(); m++)
      {
        if (mesh_coupl_[m]->GetName() == name) return m;
      }

      for (int l = 0; l < NumLevelSetCoupling(); l++)
      {
        if (levelset_coupl_[l]->GetName() == name) return (NumMeshCoupling() + l);
      }

      return coup_idx;
    }


    int GetMeshCouplingIndex(const std::string& name)
    {
      for (int m = 0; m < NumMeshCoupling(); m++)
      {
        if (mesh_coupl_[m]->GetName() == name) return m;
      }
      return -1;
    }

    Teuchos::RCP<LevelSetCoupling> GetLevelSetCoupling(const std::string& name)
    {
      for (int l = 0; l < NumLevelSetCoupling(); l++)
      {
        if (levelset_coupl_[l]->GetName() == name) return levelset_coupl_[l];
      }
      return Teuchos::null;
    }


    INPAR::XFEM::AveragingStrategy GetAveragingStrategy(
        const int coup_sid,  ///< the global id of the coupling side
        const int back_eid   ///< the global element id of the background mesh
    )
    {
      if (IsLevelSetCoupling(coup_sid))
      {
        // TODO: currently only one level-set field supported
        // get the level-set coupling object index for given background element
        const int lsc_idx = GetLevelSetCouplingIndex(back_eid);

        return levelset_coupl_[lsc_idx]->GetAveragingStrategy();
      }
      else if (IsMeshCoupling(coup_sid))
      {
        // get the mesh coupling object index
        const int mc_idx = GetMeshCouplingIndex(coup_sid);

        return mesh_coupl_[mc_idx]->GetAveragingStrategy();
      }
      else
        dserror(
            "there is no valid mesh-/levelset-coupling condition object for side: %i", coup_sid);

      return INPAR::XFEM::invalid;
    }

    /// ...
    int GetMeshCouplingIndex(const int coup_sid)
    {
      // safety checks
      if (coup_sid < 0)
      {
        //      dserror("invalid negative coupling side id %i", coup_sid);
        return -1;
      }
      if (levelset_gid_ >= 0 and coup_sid > levelset_gid_)
      {
        //      dserror("invalid coupling side id %i", coup_sid);
        return -1;
      }

      if (IsLevelSetCoupling(coup_sid))
      {
        //      dserror("level-set side does not have a cutter discretization. Why do you call
        //      this?");
        return -1;
      }


      const int num_mesh_coupl = mesh_coupl_.size();
      if (num_mesh_coupl == 0)
      {
        //      dserror("no mesh coupling objects available?!");
        return -1;
      }

      for (int idx = (int)mesh_coupl_.size() - 1; idx >= 0; idx--)  // inverse loop
      {
        if (coup_sid >= mesh_coupl_start_gid_[idx])
        {
          return idx;
        }
      }

      dserror("no valid mesh coupling index found!");
      return -1;
    }


    /// ...
    int GetLevelSetCouplingIndex(const int back_eid  ///< global element id
    )
    {
      // find out which level-set index is the active one for the given background element

      const Epetra_Map* elecolmap = bg_dis_->ElementColMap();
      const int lid = elecolmap->LID(back_eid);

      const int lsc_idx = (*ele_lsc_coup_idx_col_)[lid];

      return lsc_idx;
    }

    int GetCouplingIndex(const int coup_sid, const int back_eid)
    {
      int coup_idx = -1;

      if (IsLevelSetCoupling(coup_sid))
        coup_idx = NumMeshCoupling() + GetLevelSetCouplingIndex(back_eid);
      else
        coup_idx = GetMeshCouplingIndex(coup_sid);

      return coup_idx;
    }

    // Get Boundary Cell Clone Information <clone_coup_idx, clone_coup_sid>
    std::vector<std::pair<int, int>> GetBCCloneInformation(
        const int coup_sid, const int back_eid, int coup_idx = -1);

    int GetLevelSetCouplingGid() { return levelset_gid_; }

    /// check if the given coupling side corresponds the unique level-set side
    bool IsLevelSetCoupling(const int coupl_sid) { return coupl_sid == levelset_gid_; }

    bool IsMeshCoupling(const int coup_sid) { return GetMeshCouplingIndex(coup_sid) != -1; }

    bool HasLevelSetCoupling() { return levelset_coupl_.size() > 0; }

    bool HasMeshCoupling() { return mesh_coupl_.size() > 0; }

    int NumCoupling() { return (NumMeshCoupling() + NumLevelSetCoupling()); }

    Teuchos::RCP<CouplingBase> GetCouplingByIdx(const int coup_idx)
    {
      if (coup_idx >= NumMeshCoupling())
        return GetLevelSetCoupling(coup_idx - NumMeshCoupling());
      else if (coup_idx >= 0)
        return GetMeshCoupling(coup_idx);
      else
        return Teuchos::null;

      return Teuchos::null;
    }


    int NumMeshCoupling() { return mesh_coupl_.size(); }

    int NumLevelSetCoupling() { return levelset_coupl_.size(); }


    bool IsLevelSetCondition(const int coup_idx)
    {
      if (coup_idx >= NumMeshCoupling()) return true;

      return false;
    }


    bool IsMeshCondition(const int coup_idx)
    {
      if (coup_idx >= 0 and !IsLevelSetCondition(coup_idx)) return true;

      return false;
    }


    /// get the side element of the respective boundary discretization
    DRT::Element* GetSide(const int coup_sid  ///< the overall global coupling side id
    )
    {
      // get the mesh coupling object index
      const int mc_idx = GetMeshCouplingIndex(coup_sid);

      // compute the side id w.r.t the cutter discretization the side belongs to
      const int cutterdis_sid = GetCutterDisEleId(coup_sid, mc_idx);

      // get the boundary discretization, the side belongs to
      return mesh_coupl_[mc_idx]->GetSide(cutterdis_sid);
    }

    /// get the coupling element (the side for xfluid-sided averaging) for a given global coupl.
    /// side id
    DRT::Element* GetCouplingElement(const int coup_sid,  ///< the overall global coupling side id
        DRT::Element* ele);

    //! get the element from the conditioned dis for a local coupling side element id
    DRT::Element* GetCondElement(
        const int coup_sid  ///< global side element id w.r.t cutter discretization
    )
    {
      if (!IsMeshCoupling(coup_sid)) dserror("No cond. element available for non-mesh coupling!");

      // get the mesh coupling object index
      const int mc_idx = GetMeshCouplingIndex(coup_sid);

      // a map between cond. elements and side ids of the cutter dis is only available
      // for fluidfluid conditions; otherwise this is a bad request
      Teuchos::RCP<MeshCouplingFluidFluid> mc_xff =
          Teuchos::rcp_dynamic_cast<MeshCouplingFluidFluid>(mesh_coupl_[mc_idx]);
      if (mc_xff == Teuchos::null)
        dserror("Can't access cond dis elements for a given side id in non-xff cases!");
      const int cutterdis_sid = GetCutterDisEleId(coup_sid, mc_idx);
      return mc_xff->GetCondElement(cutterdis_sid);
    }

    // the cutwizard should add elements via the manager!!!

    // TODO: TransformID routines (localToGlobal, GlobalToLocal

    // get the side id w.r.t. the cutter discretization
    int GetCutterDisEleId(const int coup_sid, const int mc_idx)
    {
      return coup_sid - mesh_coupl_start_gid_[mc_idx];
    }

    // get the global coupling side id for a given mesh coupling and local side-id w.r.t. cutter
    // discretization
    int GetGlobalEleId(const int cutterdis_sid, const int mc_idx)
    {
      return cutterdis_sid + mesh_coupl_start_gid_[mc_idx];
    }

    // get the global coupling side id for a given mesh coupling and local side-id w.r.t. cutter
    // discretization
    int GetMeshCouplingStartGID(const int mc_idx) { return mesh_coupl_start_gid_[mc_idx]; }

    EleCoupCond GetCouplingCondition(const int coup_sid,  ///< the global id of the coupling side
        const int back_eid  ///< the global element id of the background mesh
    )
    {
      if (IsLevelSetCoupling(coup_sid))
      {
        // TODO: currently only one level-set field supported
        // get the level-set coupling object index for given background element
        const int lsc_idx = GetLevelSetCouplingIndex(back_eid);

        return levelset_coupl_[lsc_idx]->GetCouplingCondition(back_eid);
      }
      else if (IsMeshCoupling(coup_sid))
      {
        // get the mesh coupling object index
        const int mc_idx = GetMeshCouplingIndex(coup_sid);

        // compute the side id w.r.t the cutter discretization the side belongs to
        const int cutterdis_sid = GetCutterDisEleId(coup_sid, mc_idx);

        return mesh_coupl_[mc_idx]->GetCouplingCondition(cutterdis_sid);
      }
      else
        dserror(
            "there is no valid mesh-/levelset-coupling condition object for side: %i", coup_sid);

      return EleCoupCond(INPAR::XFEM::CouplingCond_NONE, NULL);
    }


    bool IsCoupling(const int coup_sid,  ///< the global id of the coupling side
        const int back_eid               ///< the global element id of the background mesh
    )
    {
      const EleCoupCond& coup_cond = GetCouplingCondition(coup_sid, back_eid);

      return IsCouplingCondition(coup_cond.first);
    }

    /// have coupling matrices to be evaluated or not?
    bool IsCouplingCondition(const std::string& cond_name)
    {
      return IsCouplingCondition(CondType_stringToEnum(cond_name));
    }

    /// have coupling matrices to be evaluated or not?
    bool IsCouplingCondition(const INPAR::XFEM::EleCouplingCondType& cond_type)
    {
      switch (cond_type)
      {
        case INPAR::XFEM::CouplingCond_SURF_FSI_MONO:
        case INPAR::XFEM::CouplingCond_SURF_FPI_MONO:
        case INPAR::XFEM::CouplingCond_SURF_FLUIDFLUID:
        case INPAR::XFEM::CouplingCond_LEVELSET_TWOPHASE:
        case INPAR::XFEM::CouplingCond_LEVELSET_COMBUSTION:
        {
          return true;
          break;
        }
        case INPAR::XFEM::CouplingCond_SURF_FSI_PART:
        case INPAR::XFEM::CouplingCond_SURF_WEAK_DIRICHLET:
        case INPAR::XFEM::CouplingCond_SURF_NEUMANN:
        case INPAR::XFEM::CouplingCond_SURF_NAVIER_SLIP:
        case INPAR::XFEM::CouplingCond_SURF_NAVIER_SLIP_TWOPHASE:
        case INPAR::XFEM::CouplingCond_LEVELSET_WEAK_DIRICHLET:
        case INPAR::XFEM::CouplingCond_LEVELSET_NEUMANN:
        case INPAR::XFEM::CouplingCond_LEVELSET_NAVIER_SLIP:
        {
          return false;
          break;
        }
        default:
          dserror("coupling condition type not known %i", cond_type);
          break;
      }

      return false;
    }

    void SetLevelSetField(const double time);

    void SetPointerBetweenCouplingObjects();

    void WriteAccess_GeometricQuantities(Teuchos::RCP<Epetra_Vector>& scalaraf,
        Teuchos::RCP<Epetra_MultiVector>& smoothed_gradphiaf,
        Teuchos::RCP<Epetra_Vector>& curvatureaf);

    void ExportGeometricQuantities();

    Teuchos::RCP<Epetra_Vector>& GetLevelSetField()
    {
      if (!is_levelset_uptodate_)
        UpdateLevelSetField();  // update the unique level-set field based on the background
                                // discretization

      return bg_phinp_;
    }

    Teuchos::RCP<const Epetra_Vector> GetLevelSetFieldCol();

    void ClearState();

    void SetState();

    void SetStateDisplacement();

    /// update interface field state vectors
    void UpdateStateVectors();

    void CompleteStateVectors();

    void ZeroStateVectors_FSI();

    void GmshOutput(const std::string& filename_base, const int step, const int gmsh_step_diff,
        const bool gmsh_debug_out_screen);

    void GmshOutputDiscretization(std::ostream& gmshfilecontent);

    void Output(const int step, const double time, const bool write_restart_data);

    /// compute lift and drag values by integrating the true residuals
    void LiftDrag(const int step, const double time);

    void ReadRestart(const int step);

    void PrepareSolve();

    bool HasMovingInterface();

    bool HasAveragingStrategy(INPAR::XFEM::AveragingStrategy strategy);

    void GetCouplingEleLocationVector(const int coup_sid, std::vector<int>& patchlm);

    /// Get the average weights from the coupling objects
    void GetAverageWeights(const int coup_sid,  ///< the overall global coupling side id
        DRT::Element* xfele,                    ///< xfluid ele
        double& kappa_m,                        ///< Weight parameter (parameter +/master side)
        double& kappa_s,                        ///< Weight parameter (parameter -/slave  side)
        bool& non_xfluid_coupling);

    /// compute viscous part of Nitsche's penalty term scaling for Nitsche's method
    void Get_ViscPenalty_Stabfac(const int coup_sid,  ///< the overall global coupling side id
        DRT::Element* xfele,                          ///< xfluid ele
        const double& kappa_m,  ///< Weight parameter (parameter +/master side)
        const double& kappa_s,  ///< Weight parameter (parameter -/slave  side)
        const double& inv_h_k,  ///< the inverse characteristic element length h_k
        const DRT::ELEMENTS::FluidEleParameterXFEM*
            params,                 ///< parameterlist which specifies interface configuration
        double& NIT_visc_stab_fac,  ///< viscous part of Nitsche's penalty term
        double&
            NIT_visc_stab_fac_tang  ///< viscous part of Nitsche's penalty term in tang direction
    );

    /// get the estimation of the penalty scaling in Nitsche's method from the trace inequality for
    /// a specific coupling side
    double Get_TraceEstimate_MaxEigenvalue(
        const int coup_sid  ///< the overall global coupling side id
    );

    /// set material pointer for volume
    void GetVolumeCellMaterial(
        DRT::Element* actele, Teuchos::RCP<MAT::Material>& mat, const GEO::CUT::VolumeCell* vc);

    /// set material pointer for volume cell for (coupling) master side
    void GetInterfaceMasterMaterial(
        DRT::Element* actele, Teuchos::RCP<MAT::Material>& mat, const GEO::CUT::VolumeCell* vc);

    /// set material pointer for coupling slave side
    void GetInterfaceSlaveMaterial(
        DRT::Element* actele, Teuchos::RCP<MAT::Material>& mat, int coup_sid);

    /// Initialize Fluid Intersection/Cut State
    void InitializeFluidState(Teuchos::RCP<GEO::CutWizard> cutwizard,
        Teuchos::RCP<DRT::Discretization> fluiddis,
        Teuchos::RCP<XFEM::ConditionManager> condition_manager,
        Teuchos::RCP<Teuchos::ParameterList> fluidparams);

   public:
    //! initialized the coupling object
    void Init();

    //! setup the coupling object
    void Setup();

    /// get the indicator state
    inline const bool& IsInit() const { return isinit_; };

    /// get the indicator state
    inline const bool& IsSetup() const { return issetup_; };

    /// Check if Init() and Setup() have been called, yet.
    inline void CheckInitSetup() const
    {
      if (!IsInit() or !IsSetup()) dserror("Call Init() and Setup() first!");
    }

    /// Check if Init() has been called
    inline void CheckInit() const
    {
      if (not IsInit()) dserror("Call Init() first!");
    }

   private:
    // build the whole object which then can be used
    void Create();

    ///
    void UpdateLevelSetField();

    /// combine two levelset fields via boolean type set operations and set result into vec1
    void CombineLevelSetField(Teuchos::RCP<Epetra_Vector>& vec1, Teuchos::RCP<Epetra_Vector>& vec2,
        const int lsc_index_2, Teuchos::RCP<Epetra_IntVector>& node_lsc_coup_idx,
        XFEM::CouplingBase::LevelSetBooleanType ls_boolean_type);

    /// check if the vector maps are equal
    void CheckForEqualMaps(
        const Teuchos::RCP<Epetra_Vector>& vec1, const Teuchos::RCP<Epetra_Vector>& vec2);

    /// combine two levelset fields via boolean type "union" set operation and put result into vec1
    void SetMinimum(Teuchos::RCP<Epetra_Vector>& vec1, Teuchos::RCP<Epetra_Vector>& vec2,
        const int lsc_index_2, Teuchos::RCP<Epetra_IntVector>& node_lsc_coup_idx);

    /// combine two levelset fields via boolean type "cut" set operation and put result into vec1
    void SetMaximum(Teuchos::RCP<Epetra_Vector>& vec1, Teuchos::RCP<Epetra_Vector>& vec2,
        const int lsc_index_2, Teuchos::RCP<Epetra_IntVector>& node_lsc_coup_idx);

    /// combine two levelset fields via boolean type "difference" set operation and put result into
    /// vec1
    void SetDifference(Teuchos::RCP<Epetra_Vector>& vec1, Teuchos::RCP<Epetra_Vector>& vec2,
        const int lsc_index_2, Teuchos::RCP<Epetra_IntVector>& node_lsc_coup_idx);

    /// combine two levelset fields via boolean type "sym_difference" set operation and put result
    /// into vec1
    void SetSymmetricDifference(Teuchos::RCP<Epetra_Vector>& vec1,
        Teuchos::RCP<Epetra_Vector>& vec2, const int lsc_index_2,
        Teuchos::RCP<Epetra_IntVector>& node_lsc_coup_idx);

    void BuildComplementaryLevelSet(Teuchos::RCP<Epetra_Vector>& vec1);

    ///<
    std::map<std::string, int> dofset_coupling_map_;

    ///< background discretiaztion w.r.t for which the couling manager is constructed
    Teuchos::RCP<DRT::Discretization> bg_dis_;

    /// mesh coupling objects
    std::vector<Teuchos::RCP<MeshCoupling>> mesh_coupl_;

    /// level-set coupling objects
    std::vector<Teuchos::RCP<LevelSetCoupling>> levelset_coupl_;

    /// starting index for element side global id
    std::vector<int> mesh_coupl_start_gid_;

    /// index for the unique level-set-side global id
    int levelset_gid_;

    /// global number of mesh and level-set coupling sides over all processors
    int numglobal_coupling_sides;

    /// time
    double time_;

    /// time step
    int step_;

    //! @name state vectors based on background discretization

    //! background-dis state vectors for levelset applications
    bool is_levelset_uptodate_;
    Teuchos::RCP<Epetra_IntVector> ele_lsc_coup_idx_col_;
    Teuchos::RCP<Epetra_Vector> bg_phinp_;
    //@}

    bool isinit_;  //! is conditionmanager initialized

    bool issetup_;  //! is conditionmanager set up
  };

}  // namespace XFEM

#endif /*XFEM_CONDITION_MANAGER_H_*/
