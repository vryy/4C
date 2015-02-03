/*!-----------------------------------------------------------------------------------------------*
\file xfem_edgestab.H

\brief provides the xfem fluid and ghost penalty stabilization based on EOS/CIP (edge-oriented, continuous interior penalty) scheme

<pre>
Maintainer: Benedikt Schott
            schott@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>
 *------------------------------------------------------------------------------------------------*/


#ifndef XFEM_CONDITION_MANAGER_H_
#define XFEM_CONDITION_MANAGER_H_


#include <vector>
#include <Teuchos_RCP.hpp>
//
#include "../drt_inpar/inpar_xfem.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_globalproblem.H"


#include <Epetra_Vector.h>


namespace DRT
{
class AssembleStrategy;
class Condition;
class Discretization;
class DiscretizationXFEM;
class DiscretizationFaces;
class Element;

namespace ELEMENTS
{
class Fluid;
class FluidIntFace;
}
}

namespace GEO
{
class CutWizard;

namespace CUT
{
class SideHandle;
}
}

namespace LINALG
{
class SparseMatrix;
}


namespace XFEM
{


typedef std::pair<INPAR::XFEM::EleCouplingCondType, DRT::Condition*> EleCoupCond;


static INPAR::XFEM::EleCouplingCondType CondType_stringToEnum(const std::string& condname)
{
  if     (condname == "XFEMSurfFSIPart")            return INPAR::XFEM::CouplingCond_SURF_FSI_PART;
  else if(condname == "XFEMSurfFSIMono")            return INPAR::XFEM::CouplingCond_SURF_FSI_MONO;
  else if(condname == "XFEMSurfCrackFSIPart")       return INPAR::XFEM::CouplingCond_SURF_CRACK_FSI_PART;
  else if(condname == "XFEMSurfFluidFluid")         return INPAR::XFEM::CouplingCond_SURF_FLUIDFLUID;
  else if(condname == "XFEMLevelsetWeakDirichlet")  return INPAR::XFEM::CouplingCond_LEVELSET_WEAK_DIRICHLET;
  else if(condname == "XFEMLevelsetNeumann")        return INPAR::XFEM::CouplingCond_LEVELSET_NEUMANN;
  else if(condname == "XFEMLevelsetTwophase")       return INPAR::XFEM::CouplingCond_LEVELSET_TWOPHASE;
  else if(condname == "XFEMLevelsetCombustion")     return INPAR::XFEM::CouplingCond_LEVELSET_COMBUSTION;
  else if(condname == "XFEMSurfWeakDirichlet")      return INPAR::XFEM::CouplingCond_SURF_WEAK_DIRICHLET;
  else if(condname == "XFEMSurfNeumann")            return INPAR::XFEM::CouplingCond_SURF_NEUMANN;
  //else dserror("condition type not supported: %s", condname.c_str());

  return INPAR::XFEM::CouplingCond_NONE;
}



class CouplingBase
{

public:


  //! constructor
  explicit CouplingBase(
      Teuchos::RCP<DRT::Discretization>&  bg_dis,    ///< background discretization
      const std::string &                 cond_name, ///< name of the condition, by which the derived cutter discretization is identified
      Teuchos::RCP<DRT::Discretization>&  cond_dis,  ///< full discretization from which the cutter discretization is derived
      const double                        time,      ///< time
      const int                           step       ///< time step
  ) :
  bg_dis_(bg_dis),
  cond_name_(cond_name),
  cond_dis_(cond_dis),
  cutter_dis_(Teuchos::null),
  coupl_dis_(Teuchos::null),
  coupling_strategy_(INPAR::XFEM::invalid),
  myrank_(bg_dis_->Comm().MyPID()),
  dt_(-1.0),
  time_(time),
  step_(step)
  {
  }

  //! destructor
  ~CouplingBase()
  {
    return;
  }

  void SetTimeAndStep(
      const double time,
      const int step
  )
  {
    time_ = time;
    step_ = step;
  }

  void IncrementTimeAndStep(
      const double dt
  )
  {
    dt_    = dt;
    time_ += dt;
    step_ += 1;
  }



  const EleCoupCond & GetCouplingCondition(
      const int gid ///< global element element id w.r.t cutter discretization (bgele->Id for LevelsetCoupling cut and side-Id for MeshCoupling)
  )
  {
    int lid = cutter_dis_->ElementColMap()->LID(gid);
    return cutterele_conds_[lid];
  }

  const std::string & GetName(){return cond_name_;}

  Teuchos::RCP<DRT::Discretization> GetCutterDis(){return cutter_dis_;}
  Teuchos::RCP<DRT::Discretization> GetCouplingDis(){return coupl_dis_;}
  Teuchos::RCP<DRT::Discretization> GetCondDis(){return cond_dis_;}

  const INPAR::XFEM::CouplingStrategy GetCouplingStrategy(){return coupling_strategy_;}

  virtual void PrepareSolve(){};

  virtual bool HasMovingInterface() = 0;

protected:



  void SetElementConditions();

  void SetCouplingStrategy();

  void SetCouplingDiscretization();

  ///< background discretization
  Teuchos::RCP<DRT::Discretization> bg_dis_;

  ///------------------------
  // CUTTER-DISCRETIZATION specific member
  ///------------------------

  ///< name of the condition, by which the derived cutter discretization is identified
  const std::string cond_name_;

  ///< discretization from which the cutter discretization is derived
  Teuchos::RCP<DRT::Discretization> cond_dis_;

  ///< discretization w.r.t which the interface is described and w.r.t which the state vectors describing the interface position are defined (bgdis for LevelSetCoupling and boundary dis for MeshCoupling)
  Teuchos::RCP<DRT::Discretization> cutter_dis_;

  ///< pairs of condition type and pointer to DRT::Condition for all column elements of the cutter discretization (bgdis for LevelSetCoupling and boundary dis for MeshCoupling)
  std::vector<EleCoupCond> cutterele_conds_;

  std::vector<std::string> conditions_to_copy_;   ///< list of conditions that will be copied to the new discretization and used to set for each cutter element


  ///------------------------
  // Coupling-DISCRETIZATION specific member
  ///------------------------

  ///< discretization with which the background discretization is coupled (structural dis, fluid dis, poro dis, scatra dis, boundary dis),
  ///< Teuchos::null in case that no coupling terms but only boundary terms are evaluated
  Teuchos::RCP<DRT::Discretization> coupl_dis_;

  ///< coupling strategy, type of weighting
  INPAR::XFEM::CouplingStrategy coupling_strategy_;



  int myrank_;

  double dt_; ///< current time step size

  double time_;

  int step_;

};


/*!
\brief
 */
class MeshCoupling : public CouplingBase
{
public:

  //! constructor
  explicit MeshCoupling(
      Teuchos::RCP<DRT::Discretization>&  bg_dis,   ///< background discretization
      const std::string &                 cond_name,///< name of the condition, by which the derived cutter discretization is identified
      Teuchos::RCP<DRT::Discretization>&  cond_dis, ///< discretization from which cutter discretization can be derived
      const double                        time,      ///< time
      const int                           step       ///< time step
  );

  //! destructor
  virtual ~MeshCoupling()
  {
    return;
  }

  /// get the side element of the respective boundary discretization
  DRT::Element* GetSide(
      const int sid ///< global side element id w.r.t cutter discretization
  )
  {
    return cutter_dis_->gElement(sid);
  }

  Teuchos::RCP<const Epetra_Vector> GetCutterDispCol();

  // finalized interface state vectors
  virtual void CompleteStateVectors(){};

  // zero interface state vectors for FSI
  virtual void ZeroStateVectors_FSI(){};

  /// set state vectors for cutter discretization
  virtual void SetState();

  /// set displacement state vectors for cutter discretization
  virtual void SetStateDisplacement();

  /// update interface field state vectors
  virtual void UpdateStateVectors();

  virtual void GmshOutput(
      const std::string & filename_base,
      const int step,
      const int gmsh_step_diff,
      const bool gmsh_debug_out_screen
  ){};

  virtual void GmshOutputDiscretization(
    std::ostream& gmshfilecontent
  );

  virtual void Output(
      const int step,
      const double time,
      const bool write_restart_data
  ){};

  virtual void LiftDrag(
      const int step,
      const double time
  ) const {};


  virtual void ReadRestart(
      const int step
  ){};

  virtual bool HasMovingInterface(){ return true; }


  Teuchos::RCP<Epetra_Vector> IVelnp(){ return ivelnp_; }
  Teuchos::RCP<Epetra_Vector> IVeln() { return iveln_;  }
  Teuchos::RCP<Epetra_Vector> IVelnm(){ return ivelnm_; }

  Teuchos::RCP<Epetra_Vector> IDispnp(){ return idispnp_; }
  Teuchos::RCP<Epetra_Vector> IDispn() { return idispn_;  }

private:


  void SetConditionsToCopy();

  /// return the name of the boundary element type to be created from the coupling discretization
  std::string GetBELEName(Teuchos::RCP<DRT::Discretization> dis)
  {
    const int ndim = DRT::Problem::Instance()->NDim();

    if(ndim != 3 ) dserror("only three-dimensional XFEM applications supported!");

    if(dis->NumMyColNodes() == 0) dserror("discretization does not have col nodes?!");
    const int numdof = dis->NumDof(dis->lColNode(0));

    switch(numdof)
    {
    case 3: {return "BELE3_3"; break;}   // use always 3 dofs on the structural surface, in partitioned and monolithic algorithms as well
    case 4: {return "BELE3_4"; break;}   // use always 4 dofs on the fluid surface
    default: dserror("coupling dis currently has %i dofs", numdof); break;
    }

    return "";
  }

  //! create cutting discretization from condition
  void CreateCutterDisFromCondition();



protected:

  void InitStateVectors();

  //------------------------------- vectors -----------------------------
  //TODO: these vectors are not required for Neumann and WDBC conditions, derive class
  //! @name cutter-dis state vectors
  Teuchos::RCP<Epetra_Vector>    ivelnp_;
  Teuchos::RCP<Epetra_Vector>    iveln_;
  Teuchos::RCP<Epetra_Vector>    ivelnm_;

  Teuchos::RCP<Epetra_Vector>    idispnp_;
  Teuchos::RCP<Epetra_Vector>    idispn_;
  //@}
};




/*!
\brief
 */
class MeshCouplingBC : public MeshCoupling
{
public:

  //! constructor
  explicit MeshCouplingBC(
      Teuchos::RCP<DRT::Discretization>&  bg_dis,   ///< background discretization
      const std::string &                 cond_name,///< name of the condition, by which the derived cutter discretization is identified
      Teuchos::RCP<DRT::Discretization>&  cond_dis,  ///< discretization from which cutter discretization can be derived
      const double                        time,      ///< time
      const int                           step       ///< time step
  );

  //! destructor
  virtual ~MeshCouplingBC()
  {
    return;
  }


private:

  void EvaluateInterfaceDisplacement(
      std::vector<double>& final_values,
      DRT::Node* node,
      DRT::Condition* cond,
      const double time
  );

  void EvaluateInterfaceVelocity(
      std::vector<double>& final_values,
      DRT::Node* node,
      DRT::Condition* cond,
      const double time,
      const double dt
  );

  void ComputeInterfaceVelocityFromDisplacement(
      std::vector<double>& final_values,
      DRT::Node* node,
      const double dt,
      const std::string* evaltype
  );

  void EvaluateFunction(
      std::vector<double>& final_values,
      const double* x,
      DRT::Condition* cond,
      const double time
  );

  void EvaluateImplementation(
      std::vector<double>& final_values,
      const double* x,
      DRT::Condition* cond,
      const double time,
      const std::string& function_name
  );

protected:

  virtual void SetInterfaceDisplacement();

  virtual void SetInterfaceVelocity();

  virtual void EvaluateCondition(
      Teuchos::RCP<Epetra_Vector> ivec,
      const std::string& condname,
      const double time,
      const double dt = 0.0
  );

  virtual void PrepareSolve();

  virtual bool HasMovingInterface();


};






/*!
\brief
 */
class MeshCouplingFSI : public MeshCoupling
{
public:

  //! constructor
  explicit MeshCouplingFSI(
      Teuchos::RCP<DRT::Discretization>&  bg_dis,   ///< background discretization
      const std::string &                 cond_name,///< name of the condition, by which the derived cutter discretization is identified
      Teuchos::RCP<DRT::Discretization>&  cond_dis,  ///< discretization from which cutter discretization can be derived
      const double                        time,      ///< time
      const int                           step       ///< time step
  );

  //! destructor
  virtual ~MeshCouplingFSI()
  {
    return;
  }

  // finalize the interface true residual vector
  virtual void CompleteStateVectors();

  virtual void ZeroStateVectors_FSI();

  virtual void GmshOutput(
      const std::string & filename_base,
      const int step,
      const int gmsh_step_diff,
      const bool gmsh_debug_out_screen
  );

  virtual void GmshOutputDiscretization(
    std::ostream& gmshfilecontent
  );

  virtual void LiftDrag(
      const int step,
      const double time
  ) const;

  virtual void ReadRestart(
      const int step
  );

  // interface foces
  Teuchos::RCP<Epetra_Vector> ITrueResidual(){ return itrueresidual_; }

  // for assembly of fluid interface forces
  Teuchos::RCP<Epetra_Vector> IForcecol(){ return iforcecol_; }


private:

  void PrepareCutterOutput();

  void Output(
      const int step,
      const double time,
      const bool write_restart_data
  );

protected:

  void InitStateVectors_FSI();

  virtual bool HasMovingInterface(){ return true; }


  //------------------------------- vectors -----------------------------
  //! @name cutter-dis state vectors
  Teuchos::RCP<Epetra_Vector>    itrueresidual_; //! interface forces acting on the structural surface (= -iforcenp)
  Teuchos::RCP<Epetra_Vector>    iforcecol_;     //! interface forces acting on the fluid surface (column vector assembly)
  //@}


  //---------------------------------output----------------------------

  //! @name output discretization writers
  Teuchos::RCP<IO::DiscretizationWriter> cutter_output_;

  bool firstoutputofrun_;
  //@}

};



/// class which allows to manipulate the boundary dis and to deal with crack specific staff
class MeshCouplingFSICrack : public MeshCouplingFSI
{
public:

  //! constructor
  explicit MeshCouplingFSICrack(
      Teuchos::RCP<DRT::Discretization>&  bg_dis,   ///< background discretization
      const std::string &                 cond_name,///< name of the condition, by which the derived cutter discretization is identified
      Teuchos::RCP<DRT::Discretization>&  cond_dis,  ///< discretization from which cutter discretization can be derived
      const double                        time,      ///< time
      const int                           step       ///< time step
  );

  //! destructor
  virtual ~MeshCouplingFSICrack()
  {
    return;
  }

  /// set new cutter discretization
  void SetCutterDis(Teuchos::RCP<DRT::Discretization> cutter_dis_new);

  void InitCrackInitiationsPoints();

  //! set the nodes representing crack tip in FSI with crack structure simulations
  void SetCrackTipNodes( std::map<int, LINALG::Matrix<3,1> > & tip_nodes )
  {
    tip_nodes_.clear();
    if( tip_nodes.size() > 0 )
      tip_nodes_.insert( tip_nodes.begin(), tip_nodes.end() );
  }

  //! get the nodes representing crack tip in FSI with crack structure simulations
  std::map<int, LINALG::Matrix<3,1> > & GetCrackTipNodes(){ return tip_nodes_; };

  void UpdateBoundaryValuesAfterCrack( const std::map<int,int>& oldnewIds );


private:

  std::map<int, LINALG::Matrix<3,1> > tip_nodes_;

  Teuchos::RCP<DRT::Discretization> cutterdis_init_output_;

};

/*!
\brief specialized class for coupling with an embedded fluid discretization
  \author kruse
  \date 01/15
 */
class MeshCouplingFluidFluid : public MeshCoupling
{
public:

  //! constructor
  explicit MeshCouplingFluidFluid(
      Teuchos::RCP<DRT::Discretization>&  bg_dis,   ///< background discretization
      const std::string &                 cond_name,///< name of the condition, by which the derived cutter discretization is identified
      Teuchos::RCP<DRT::Discretization>&  cond_dis,  ///< discretization from which cutter discretization can be derived
      const double                        time,      ///< time
      const int                           step       ///< time step
  );

  //! estimate the penalty scaling in Nitsche's method from the trace inequality
  //! via solving a local eigenvalue problem
  void EstimateNitscheTraceMaxEigenvalue(const Teuchos::RCP<const Epetra_Vector> & dispnp) const;

  //! get auxiliary coupling discretization (embedded elements with nodes in the cutting surface discretization)
  Teuchos::RCP<DRT::Discretization> GetAuxiliaryDiscretization()
  {
    return aux_coup_dis_;
  }

  //! Access to map from cutting surface element id to embedded element id
  int GetEmbeddedElementId(int cutting_id)
  {
    std::map<int,int>::iterator it = cutting_emb_gid_map_.find(cutting_id);

    if (it == cutting_emb_gid_map_.end())
      dserror("Could not find parent element ID for side %d.", cutting_id);

    return it->second;
  }

  //! Access to map from cutting surface element id to coupling element's face lid
  int GetFaceLidOfEmbeddedElement(int face_id)
  {
    std::map<int,int>::iterator it = cutting_emb_face_lid_map_.find(face_id);

    if (it == cutting_emb_face_lid_map_.end())
      dserror("No local id of the face-element % d w.r.t embedded element stored!", face_id);

    return it->second;
  }

  //! ghost interface-contributing embedded elements (required for error calculation in case
  //! of xfluid-sided coupling)
  void RedistributeForErrorCalculation();

private:

  //!
  virtual bool HasMovingInterface()
  {
    std::cout << "WARNING: adapt HasMovingInterface in MeshCouplingFluidFluid!" << std::endl;
    return true;
  }

  //! ghost embedded elements, that contribute to the cutting interface discretization on all procs
  void RedistributeEmbeddedDiscretization();

  //! build an auxiliary discretization out of the elements, that contribute to the cutting discretization
  void CreateAuxiliaryDiscretization();

  //! create map from cutting surface element id to embedded element id
  void CreateCuttingToEmbeddedElementMap();

  //! auxiliary discretization, based on the elements of cond_dis, that contribute to
  //! the elements of cutter_dis_ with at least one edge
  Teuchos::RCP<DRT::Discretization> aux_coup_dis_;

  //! map from cutting surface element id to embedded element id
  std::map<int,int> cutting_emb_gid_map_;

  //! map from cutting surface element id to embedded element's face lid
  std::map<int,int> cutting_emb_face_lid_map_;

};

/*!
\brief
 */
class LevelSetCoupling : public CouplingBase
{
public:

  //! constructor
  explicit LevelSetCoupling(
      Teuchos::RCP<DRT::Discretization>&  bg_dis,   ///< background discretization
      const std::string &                 cond_name, ///< name of the condition, by which the derived cutter discretization is identified
      const double                        time,      ///< time
      const int                           step       ///< time step
  );


  //! destructor
  virtual ~LevelSetCoupling()
  {
    return;
  }


  /// set levelset field by function
  void SetLevelSetField(const double time);

  /// set levelset field by given vector
  void SetLevelSetField(
      Teuchos::RCP<const Epetra_Vector> scalaraf,
      Teuchos::RCP<DRT::Discretization> scatradis
  );

  Teuchos::RCP<Epetra_Vector> GetLevelSetField(){ return phinp_; }

  virtual void ReadRestart(
      const int step
  );

  virtual bool HasMovingInterface(){ return true; }

private:

  void SetConditionsToCopy();

  //! @name state vectors

  //! scatra-dis (cutterdis) state vectors for levelset applications, prepares nonmatching discretizations between scatra and fluid
  Teuchos::RCP<Epetra_Vector>    phinp_;
  //@}

};





/*!
\brief provides the xfem fluid and ghost penalty stabilization based on EOS/CIP (edge-oriented, continuous interior penalty) scheme
 */
class ConditionManager
{
public:

  //! constructor
  explicit ConditionManager(
      Teuchos::RCP<DRT::Discretization> &               bg_dis,          ///< background discretization
      std::vector<Teuchos::RCP<DRT::Discretization> > & meshcoupl_dis,   ///< mesh coupling discretizations
      const double                                      time,            ///< time
      const int                                         step             ///< time step
  );

  //! destructor
  ~ConditionManager()
  {
    return;
  }

  void SetTimeAndStep(
      const double time,
      const int step
  );

  void IncrementTimeAndStep(
      const double dt);

  void CreateNewLevelSetCoupling(const std::string& cond_name);

  void CreateNewMeshCoupling(
      const std::string& cond_name,
      Teuchos::RCP<DRT::Discretization> cond_dis         ///< discretization from which the cutter discretization can be derived
  );

  /// create a new mesh-coupling object based on the given coupling discretization
  void AddMeshCoupling(
      const std::string &                 cond_name,
      Teuchos::RCP<DRT::Discretization>   cond_dis
  )
  {
    const int myrank = cond_dis->Comm().MyPID();
    if(CondType_stringToEnum(cond_name) == INPAR::XFEM::CouplingCond_SURF_CRACK_FSI_PART)
    {
      if (myrank == 0) std::cout << " create new MeshCouplingFSICrack" << std::endl;
      mesh_coupl_.push_back(Teuchos::rcp(new MeshCouplingFSICrack(bg_dis_, cond_name, cond_dis, time_, step_)));
    }
    else if(CondType_stringToEnum(cond_name) == INPAR::XFEM::CouplingCond_SURF_FSI_PART or
        CondType_stringToEnum(cond_name) == INPAR::XFEM::CouplingCond_SURF_FSI_MONO)
    {
      if (myrank == 0) std::cout << " create new MeshCouplingFSI" << std::endl;
      mesh_coupl_.push_back(Teuchos::rcp(new MeshCouplingFSI(bg_dis_, cond_name, cond_dis, time_, step_)));
    }
    else if(CondType_stringToEnum(cond_name) == INPAR::XFEM::CouplingCond_SURF_WEAK_DIRICHLET or
        CondType_stringToEnum(cond_name) == INPAR::XFEM::CouplingCond_SURF_NEUMANN)
    {
      if (myrank == 0) std::cout << " create new MeshCouplingBC" << std::endl;
      mesh_coupl_.push_back(Teuchos::rcp(new MeshCouplingBC(bg_dis_, cond_name, cond_dis, time_, step_)));
    }
    else if(CondType_stringToEnum(cond_name) == INPAR::XFEM::CouplingCond_SURF_FLUIDFLUID)
    {
      if (myrank == 0) std::cout << " create new MeshCoupling" << std::endl;
      mesh_coupl_.push_back(Teuchos::rcp(new MeshCouplingFluidFluid(bg_dis_, cond_name, cond_dis, time_, step_)));
    }
    else
    {
      if (myrank == 0) std::cout << " create new MeshCoupling" << std::endl;
      mesh_coupl_.push_back(Teuchos::rcp(new MeshCoupling(bg_dis_, cond_name, cond_dis, time_, step_)));
    }
  }

  /// add a new level-set-coupling object based on the given coupling discretization
  void AddLevelSetCoupling(
      const std::string &                 cond_name
  )
  {
    levelset_coupl_.push_back(Teuchos::rcp(new LevelSetCoupling(bg_dis_, cond_name, time_, step_)));
  }



  // build the whole object which then can be used
  void Create();

  //void Update(const double time); // update the meshcoupling objects (interface displ), update the level-set couplings


  /// Getters

  /// get cutter discretization the coupling side belongs to
  Teuchos::RCP<DRT::Discretization> GetCutterDis(
      const int coup_sid ///< the global id of the coupling side
  )
  {
    if(IsLevelSetCoupling(coup_sid)) return Teuchos::null;

    // get the mesh coupling object index
    int mc_idx =  GetMeshCouplingIndex(coup_sid);

    return mesh_coupl_[mc_idx]->GetCutterDis();
  }

  /// get cutter discretization the coupling side belongs to
  Teuchos::RCP<DRT::Discretization> GetCouplingDis(
      const int coup_sid ///< the global id of the coupling side
  )
  {
    if(IsLevelSetCoupling(coup_sid)) return bg_dis_;

    // get the mesh coupling object index
    int mc_idx =  GetMeshCouplingIndex(coup_sid);

    return mesh_coupl_[mc_idx]->GetCouplingDis();
  }

  Teuchos::RCP<MeshCoupling> GetMeshCoupling( const int mc_idx)
      {
    if(mc_idx < (int)mesh_coupl_.size() and mc_idx >= 0)
      return mesh_coupl_[mc_idx];

    return Teuchos::null;
      }

  Teuchos::RCP<LevelSetCoupling> GetLevelSetCoupling( const int ls_idx)
      {
    if(ls_idx < (int)levelset_coupl_.size() and ls_idx >= 0)
      return levelset_coupl_[ls_idx];

    return Teuchos::null;
      }

  Teuchos::RCP<CouplingBase> GetCoupling(const std::string & name)
      {
    Teuchos::RCP<CouplingBase> coupling = Teuchos::null;

    coupling=GetMeshCoupling(name);

    if(coupling!=Teuchos::null) return coupling;

    coupling=GetLevelSetCoupling(name);

    if(coupling!=Teuchos::null) return coupling;

    return Teuchos::null;
      }


  Teuchos::RCP<MeshCoupling> GetMeshCoupling(const std::string & name)
      {
    for(int m=0; m<NumMeshCoupling(); m++)
    {
      if(mesh_coupl_[m]->GetName() == name)
        return mesh_coupl_[m];
    }
    return Teuchos::null;
      }

  int GetMeshCouplingIndex(const std::string & name)
  {
    for(int m=0; m<NumMeshCoupling(); m++)
    {
      if(mesh_coupl_[m]->GetName() == name)
        return m;
    }
    return -1;
  }

  Teuchos::RCP<LevelSetCoupling> GetLevelSetCoupling(const std::string & name)
      {
    for(int l=0; l<NumLevelSetCoupling(); l++)
    {
      if(levelset_coupl_[l]->GetName() == name)
        return levelset_coupl_[l];
    }
    return Teuchos::null;
      }


  const INPAR::XFEM::CouplingStrategy GetCouplingStrategy(
      const int coup_sid, ///< the global id of the coupling side
      const int back_eid  ///< the global element id of the background mesh
  )
  {
    if(IsLevelSetCoupling(coup_sid))
    {
      //TODO: currently only one level-set field supported
      // get the level-set coupling object index for given background element
      const int lsc_idx = GetLevelSetCouplingIndex(back_eid);

      return levelset_coupl_[lsc_idx]->GetCouplingStrategy();
    }
    else if(IsMeshCoupling(coup_sid))
    {
      // get the mesh coupling object index
      const int mc_idx =  GetMeshCouplingIndex(coup_sid);

      return mesh_coupl_[mc_idx]->GetCouplingStrategy();
    }
    else dserror("there is no valid mesh-/levelset-coupling condition object for side: %i", coup_sid);

    return INPAR::XFEM::invalid;
  }

  /// ...
  int GetMeshCouplingIndex(
      const int coup_sid
  )
  {
    // safety checks
    if(coup_sid < 0)
    {
      //      dserror("invalid negative coupling side id %i", coup_sid);
      return -1;
    }
    if(levelset_gid_ >= 0 and coup_sid > levelset_gid_ )
    {
      //      dserror("invalid coupling side id %i", coup_sid);
      return -1;
    }

    if(IsLevelSetCoupling(coup_sid))
    {
      //      dserror("level-set side does not have a cutter discretization. Why do you call this?");
      return -1;
    }


    const int num_mesh_coupl = mesh_coupl_.size();
    if(num_mesh_coupl == 0)
    {
      //      dserror("no mesh coupling objects available?!");
      return -1;
    }



    int mc_idx = 0; // mesh coupling index to be set, init to 0 as coup_sid >= 0

    for(int idx=0; idx < (int)mesh_coupl_.size(); idx++)
    {
      if(coup_sid < mesh_coupl_start_gid_[idx]) // for the first run this is guaranteed to be false
      {
        mc_idx = idx-1;  // take the last mc-index
      }
    }

    if(mc_idx == -1 or mc_idx >= num_mesh_coupl)
    {
      dserror("invalid mesh coupling index %i", mc_idx);
      return -1;
    }

    return mc_idx;
  }


  /// ...
  int GetLevelSetCouplingIndex(
      const int back_eid
  )
  {
    // TODO: find out which level-set index is the active one for the given background element
    int lsc_idx = 0;

    return lsc_idx;
  }

  int GetCouplingIndex(const int coup_sid, const int back_eid)
  {
    int coup_idx = -1;

    if(IsLevelSetCoupling(coup_sid))
      coup_idx = NumMeshCoupling() + GetLevelSetCouplingIndex(back_eid);
    else
      coup_idx = GetMeshCouplingIndex(coup_sid);

    return coup_idx;
  }


  int GetLevelSetCouplingGid(){return levelset_gid_;}

  /// check if the given coupling side corresponds the unique level-set side
  bool IsLevelSetCoupling(const int coupl_sid){return coupl_sid == levelset_gid_;}

  bool IsMeshCoupling(const int coup_sid){return GetMeshCouplingIndex(coup_sid)!=-1;}

  bool HasLevelSetCoupling(){return levelset_coupl_.size() > 0;}

  bool HasMeshCoupling(){return mesh_coupl_.size() > 0;}

  int NumCoupling(){return (NumMeshCoupling() + NumLevelSetCoupling());}

  Teuchos::RCP<CouplingBase> GetCouplingByIdx(const int coup_idx)
      {
    if(coup_idx >= NumMeshCoupling())
      return GetLevelSetCoupling(coup_idx-NumMeshCoupling());
    else if(coup_idx >= 0)
      return GetMeshCoupling(coup_idx);
    else return Teuchos::null;

    return Teuchos::null;
      }



  int NumMeshCoupling(){return mesh_coupl_.size();}

  int NumLevelSetCoupling(){return levelset_coupl_.size();}


  bool IsLevelSetCondition(const int coup_idx)
  {
    if(coup_idx >= NumMeshCoupling()) return true;

    return false;
  }


  bool IsMeshCondition(const int coup_idx)
  {
    if(coup_idx >= 0 and !IsLevelSetCondition(coup_idx)) return true;

    return false;
  }


  /// get the side element of the respective boundary discretization
  DRT::Element* GetSide(
      const int coup_sid ///< the overall global coupling side id
  )
  {
    // get the mesh coupling object index
    const int mc_idx =  GetMeshCouplingIndex(coup_sid);

    // compute the side id w.r.t the cutter discretization the side belongs to
    const int cutterdis_sid = GetCutterDisEleId(coup_sid, mc_idx);

    // get the boundary discretization, the side belongs to
    return mesh_coupl_[mc_idx]->GetSide(cutterdis_sid);
  }

  // the cutwizard should add elements via the manager!!!

  // TODO: TransformID routines (localToGlobal, GlobalToLocal

  // get the side id w.r.t. the cutter discretization
  const int GetCutterDisEleId(const int coup_sid, const int mc_idx)
  {
    return coup_sid-mesh_coupl_start_gid_[mc_idx];
  }

  // get the global coupling side id for a given mesh coupling and local side-id w.r.t. cutter discretization
  const int GetGlobalEleId(const int cutterdis_sid, const int mc_idx)
  {
    return cutterdis_sid+mesh_coupl_start_gid_[mc_idx];
  }

  // get the global coupling side id for a given mesh coupling and local side-id w.r.t. cutter discretization
  const int GetMeshCouplingStartGID(const int mc_idx)
  {
    return mesh_coupl_start_gid_[mc_idx];
  }

  const EleCoupCond GetCouplingCondition(
      const int coup_sid, ///< the global id of the coupling side
      const int back_eid  ///< the global element id of the background mesh
  )
  {
    if(IsLevelSetCoupling(coup_sid))
    {
      //TODO: currently only one level-set field supported
      // get the level-set coupling object index for given background element
      const int lsc_idx = GetLevelSetCouplingIndex(back_eid);

      return levelset_coupl_[lsc_idx]->GetCouplingCondition(back_eid);
    }
    else if(IsMeshCoupling(coup_sid))
    {
      // get the mesh coupling object index
      const int mc_idx =  GetMeshCouplingIndex(coup_sid);

      // compute the side id w.r.t the cutter discretization the side belongs to
      const int cutterdis_sid = GetCutterDisEleId(coup_sid, mc_idx);

      return mesh_coupl_[mc_idx]->GetCouplingCondition(cutterdis_sid);
    }
    else dserror("there is no valid mesh-/levelset-coupling condition object for side: %i", coup_sid);

    return EleCoupCond(INPAR::XFEM::CouplingCond_NONE,NULL);
  }


  bool IsCoupling(
      const int coup_sid, ///< the global id of the coupling side
      const int back_eid  ///< the global element id of the background mesh
  )
  {
    const EleCoupCond & coup_cond = GetCouplingCondition(coup_sid, back_eid);

    return IsCouplingCondition(coup_cond.first);
  }

  /// have coupling matrices to be evaluated or not?
  bool IsCouplingCondition(
      const std::string & cond_name)
  {
    return IsCouplingCondition(CondType_stringToEnum(cond_name));
  }

  /// have coupling matrices to be evaluated or not?
  bool IsCouplingCondition(
      const INPAR::XFEM::EleCouplingCondType& cond_type)
  {
    switch(cond_type)
    {
    //    case INPAR::XFEM::CouplingCond_XFEMCOUPLING:
    //      return false;
    //      break;
    case INPAR::XFEM::CouplingCond_SURF_FSI_MONO:
    case INPAR::XFEM::CouplingCond_SURF_FLUIDFLUID:
    case INPAR::XFEM::CouplingCond_LEVELSET_TWOPHASE:
    case INPAR::XFEM::CouplingCond_LEVELSET_COMBUSTION:
    {
      return true;
      break;
    }
    case INPAR::XFEM::CouplingCond_SURF_FSI_PART:
    case INPAR::XFEM::CouplingCond_SURF_CRACK_FSI_PART:
    case INPAR::XFEM::CouplingCond_SURF_WEAK_DIRICHLET:
    case INPAR::XFEM::CouplingCond_SURF_NEUMANN:
    case INPAR::XFEM::CouplingCond_LEVELSET_WEAK_DIRICHLET:
    case INPAR::XFEM::CouplingCond_LEVELSET_NEUMANN:
    {
      return false;
      break;
    }
    default: dserror("coupling condition type not known %i", cond_type); break;
    }

    return false;
  }

  void SetLevelSetField( const double time );

  void SetLevelSetField(
      Teuchos::RCP<const Epetra_Vector> scalaraf,
      Teuchos::RCP<DRT::Discretization> scatradis
  );

  Teuchos::RCP<Epetra_Vector> GetLevelSetField()
  {
    if(!is_levelset_uptodate_) UpdateLevelSetField(); // update the unique level-set field based on the background discretization

    return bg_phinp_;
  }

  Teuchos::RCP<const Epetra_Vector> GetLevelSetFieldCol();

  void SetState();

  void SetStateDisplacement();

  /// update interface field state vectors
  void UpdateStateVectors();

  void CompleteStateVectors();

  void ZeroStateVectors_FSI();

  void GmshOutput(
    const std::string & filename_base,
    const int step,
    const int gmsh_step_diff,
    const bool gmsh_debug_out_screen
  );

  void GmshOutputDiscretization(
    std::ostream& gmshfilecontent
  );

  void Output(
      const int step,
      const double time,
      const bool write_restart_data
  );

  /// compute lift and drag values by integrating the true residuals
  void LiftDrag(
      const int step,
      const double time
  );

  void ReadRestart(
      const int step
  );

  void PrepareSolve();

  bool HasMovingInterface();

private:

  void UpdateLevelSetField();

  ///< background discretiaztion w.r.t for which the couling manager is constructed
  Teuchos::RCP<DRT::Discretization> bg_dis_;

  /// mesh coupling objects
  std::vector<Teuchos::RCP<MeshCoupling> > mesh_coupl_;

  /// level-set coupling objects
  std::vector<Teuchos::RCP<LevelSetCoupling> > levelset_coupl_;

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
  Teuchos::RCP<Epetra_Vector>    bg_phinp_;
  //@}

};

} // namespace XFEM

#endif /*XFEM_CONDITION_MANAGER_H_*/
