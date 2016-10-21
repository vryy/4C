/*!-----------------------------------------------------------------------------------------------*
\file xfem_coupling_mesh.H

\brief manages the different types of mesh based coupling conditions and thereby builds the bridge between the
xfluid class and the cut-library

\level 2

<pre>
\maintainer Benedikt Schott
            schott@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>
 *------------------------------------------------------------------------------------------------*/


#ifndef XFEM_COUPLING_MESH_H_
#define XFEM_COUPLING_MESH_H_

#include "xfem_coupling_base.H"

namespace XFEM
{

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
      const int                           coupling_id,///< id of composite of coupling conditions
      const double                        time,      ///< time
      const int                           step       ///< time step
  );

  //! destructor
  virtual ~MeshCoupling()
  {
    return;
  }

  //! get the coupling element (equal to the side for xfluid-sided, mesh-based coupling)
  virtual DRT::Element* GetCouplingElement(
      const int sid///< global side element id w.r.t cutter discretization
  )
  {
    return coupl_dis_->gElement(sid);
  }

  /// get the side element of the respective boundary discretization
  DRT::Element* GetSide(
      const int sid ///< global side element id w.r.t cutter discretization
  )
  {
    return cutter_dis_->gElement(sid);
  }

  Teuchos::RCP<const Epetra_Vector> GetCutterDispCol();

  /// fill lm vector for coupling element
  virtual void GetCouplingEleLocationVector(
      const int sid,
      std::vector<int> & patchlm);

  /// set material pointer for coupling slave side
  virtual void GetInterfaceSlaveMaterial(
    DRT::Element* actele,
    Teuchos::RCP<MAT::Material> & mat
  )
  {
    mat = Teuchos::null;
  }

  // finalized interface state vectors
  virtual void CompleteStateVectors(){};

  // zero interface state vectors for FSI
  virtual void ZeroStateVectors_FSI(){};

  /// clear state vectors
  virtual void ClearState();

  /// set state vectors for cutter discretization
  virtual void SetState();

  /// set displacement state vectors for cutter discretization
  virtual void SetStateDisplacement();

  /// update interface field state vectors
  virtual void UpdateStateVectors();

  /// update last iteration interface displacements
  virtual void UpdateDisplacementIterationVectors();

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
  );

  virtual void PrepareCutterOutput();

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
    const int numdof = dis->NumDof(0,dis->lColNode(0));

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

  Teuchos::RCP<Epetra_Vector>    idispnp_;  ///< current displacements at t^n+1
  Teuchos::RCP<Epetra_Vector>    idispn_;   ///< last displacements at t^n
  Teuchos::RCP<Epetra_Vector>    idispnpi_; ///< displacements of last Newton increment at t^n+1 (for monolithic approaches)


  //@}

  //! @name output discretization writers
  Teuchos::RCP<IO::DiscretizationWriter> cutter_output_;
  bool firstoutputofrun_;

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
      const int                           coupling_id,///< id of composite of coupling conditions
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

  virtual bool HasMovingInterface();


};

/*!
\brief
 */
class MeshCouplingWeakDirichlet : public MeshCouplingBC
{
public:

  //! constructor
  explicit MeshCouplingWeakDirichlet(
      Teuchos::RCP<DRT::Discretization>&  bg_dis,   ///< background discretization
      const std::string &                 cond_name,///< name of the condition, by which the derived cutter discretization is identified
      Teuchos::RCP<DRT::Discretization>&  cond_dis,  ///< discretization from which cutter discretization can be derived
      const int                           coupling_id,///< id of composite of coupling conditions
      const double                        time,      ///< time
      const int                           step       ///< time step
  );

  //! destructor
  virtual ~MeshCouplingWeakDirichlet()
  {
    return;
  }

public:

  void EvaluateCouplingConditions(
      LINALG::Matrix<3,1>& ivel,
      LINALG::Matrix<3,1>& itraction,
      const LINALG::Matrix<3,1>& x,
      const DRT::Condition* cond
  );

  void EvaluateCouplingConditionsOldState(
      LINALG::Matrix<3,1>& ivel,
      LINALG::Matrix<3,1>& itraction,
      const LINALG::Matrix<3,1>& x,
      const DRT::Condition* cond
  );

  virtual void PrepareSolve();

protected:

  //! Initializes configurationmap
  void InitConfigurationMap();

  //! Updates configurationmap for specific Gausspoint
  virtual void UpdateConfigurationMap_GP( double& kappa_m, double& visc_m, double& visc_s,
      double& visc_stab, double& full_stab, const LINALG::Matrix<3,1>& x, const DRT::Condition* cond);

};

/*!
\brief
 */
class MeshCouplingNeumann : public MeshCouplingBC
{
public:

  //! constructor
  explicit MeshCouplingNeumann(
      Teuchos::RCP<DRT::Discretization>&  bg_dis,   ///< background discretization
      const std::string &                 cond_name,///< name of the condition, by which the derived cutter discretization is identified
      Teuchos::RCP<DRT::Discretization>&  cond_dis,  ///< discretization from which cutter discretization can be derived
      const int                           coupling_id,///< id of composite of coupling conditions
      const double                        time,      ///< time
      const int                           step       ///< time step
  ) : MeshCouplingBC(bg_dis,cond_name,cond_dis, coupling_id, time, step)
  {}

  //! destructor
  virtual ~MeshCouplingNeumann()
  {
    return;
  }

public:

  void EvaluateCouplingConditions(
      LINALG::Matrix<3,1>& ivel,
      LINALG::Matrix<3,1>& itraction,
      const LINALG::Matrix<3,1>& x,
      const DRT::Condition* cond
  );

  void EvaluateCouplingConditionsOldState(
      LINALG::Matrix<3,1>& ivel,
      LINALG::Matrix<3,1>& itraction,
      const LINALG::Matrix<3,1>& x,
      const DRT::Condition* cond
  );

  virtual void PrepareSolve();

};


/*!
\brief
 */
class MeshCouplingNavierSlip : public MeshCouplingBC
{
public:

  //! constructor
  explicit MeshCouplingNavierSlip(
      Teuchos::RCP<DRT::Discretization>&  bg_dis,   ///< background discretization
      const std::string &                 cond_name,///< name of the condition, by which the derived cutter discretization is identified
      Teuchos::RCP<DRT::Discretization>&  cond_dis,  ///< discretization from which cutter discretization can be derived
      const int                           coupling_id,///< id of composite of coupling conditions
      const double                        time,      ///< time
      const int                           step       ///< time step
  );

  //! destructor
  virtual ~MeshCouplingNavierSlip()
  {
    return;
  }

public:

  void EvaluateCouplingConditions(
      LINALG::Matrix<3,1>& ivel,
      LINALG::Matrix<3,1>& itraction,
      LINALG::Matrix<3,3>& proj_matrix,
      const LINALG::Matrix<3,1>& x,
      const LINALG::Matrix<3,1>& normal,
      const DRT::Condition* cond,
      const bool & eval_dirich_at_gp,
      double& kappa_m,                         ///< fluid sided weighting
      double& visc_m,                          ///< fluid sided weighting
      double& visc_s                           ///< slave sided dynamic viscosity
  );

  /// this function has to be reviewed for usage with OST new.
  void EvaluateCouplingConditionsOldState(
      LINALG::Matrix<3,1>& ivel,
      LINALG::Matrix<3,1>& itraction,
      const LINALG::Matrix<3,1>& x,
      const DRT::Condition* cond
  );

  /// get the slip coefficient for this coupling
  void GetSlipCoefficient(
        double& slipcoeff,
        const LINALG::Matrix<3,1>& x,
        const DRT::Condition* cond
  );

  virtual void PrepareSolve();

protected:

  void SetConditionSpecificParameters(
      const std::string &                 cond_name ///< name of the condition, by which the derived cutter discretization is identified
  );

  void GetConditionByRobinId(
      const std::vector<DRT::Condition*> & mycond,
      const int coupling_id,
      std::vector<DRT::Condition*> & mynewcond
  );

  void CreateRobinIdMap(
      const std::vector< DRT::Condition* > & conditions_NS,
      const std::vector< DRT::Condition* > & conditions_robin,
      const std::string &                    robin_id_name,
      std::map< int ,DRT::Condition* > &     conditionsmap_robin
      );

  //! Initializes configurationmap
  void InitConfigurationMap();

  //! Updates configurationmap for specific Gausspoint
  virtual void UpdateConfigurationMap_GP( double& kappa_m, double& visc_m, double& visc_s,
      double& visc_stab, double& full_stab, const LINALG::Matrix<3,1>& x, const DRT::Condition* cond);

private:

  // Condition map. Corelating robin_id from Navier Slip condition and
  //                 Robin Dirichlet/Neumann input sections
  //        robin_id    Robin_cond
  std::map<   int   , DRT::Condition* >      conditionsmap_robin_dirch_;
  std::map<   int   , DRT::Condition* >      conditionsmap_robin_neumann_;

  //       Coupling Surface (E 1)          slip length    is slip length constant?
  std::map<        int       ,   std::pair<    double  ,          bool       > > sliplength_map_;
  //       Coupling Surface (E 1)  Force only tangential veloctiy to surface
  std::map<         int          ,     bool     > force_tangvel_map_;

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
      const int                           coupling_id,///< id of composite of coupling conditions
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

  virtual void GetSlipCoefficient(
      double& slipcoeff,
      const LINALG::Matrix<3,1>& x,
      const DRT::Condition* cond
  );

  // interface foces
  Teuchos::RCP<Epetra_Vector> ITrueResidual(){ return itrueresidual_; }

  // for assembly of fluid interface forces
  Teuchos::RCP<Epetra_Vector> IForcecol(){ return iforcecol_; }


private:

  void Output(
      const int step,
      const double time,
      const bool write_restart_data
  );

protected:

  void InitStateVectors_FSI();

  virtual bool HasMovingInterface(){ return true; }

  void SetConditionSpecificParameters(
      const std::string &                 cond_name ///< name of the condition, by which the derived cutter discretization is identified
  );

  //! Initializes configurationmap
  void InitConfigurationMap();

  //! Updates configurationmap for specific Gausspoint
  virtual void UpdateConfigurationMap_GP( double& kappa_m, double& visc_m, double& visc_s,
      double& visc_stab, double& full_stab, const LINALG::Matrix<3,1>& x, const DRT::Condition* cond);

  //------------------------------- vectors -----------------------------
  //! @name cutter-dis state vectors
  Teuchos::RCP<Epetra_Vector>    itrueresidual_; //! interface forces acting on the structural surface (= -iforcenp)
  Teuchos::RCP<Epetra_Vector>    iforcecol_;     //! interface forces acting on the fluid surface (column vector assembly)
  //@}

  //---------------------------------parameters--------------------------------
  //! @name parameters
  //       Coupling Surface (E 1)          slip length    is slip length constant?
  std::map<        int       ,   std::pair<    double  ,          bool       > > sliplength_map_;
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
      const int                           coupling_id,///< id of composite of coupling conditions
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
      const int                           coupling_id,///< id of composite of coupling conditions
      const double                        time,      ///< time
      const int                           step       ///< time step
  );

  void GetCouplingEleLocationVector(
    const int sid,
    std::vector<int> & patchlm);

  //! estimate the penalty scaling in Nitsche's method from the trace inequality
  //! via solving a local eigenvalue problem
  void EstimateNitscheTraceMaxEigenvalue(const Teuchos::RCP<const Epetra_Vector> & dispnp) const;

  //! get auxiliary coupling discretization (embedded elements with nodes in the cutting surface discretization)
  Teuchos::RCP<DRT::Discretization> GetAuxiliaryDiscretization()
  {
    return aux_coup_dis_;
  }

  //! get the coupling element for a local coupling side element id
  DRT::Element* GetCouplingElement(
    const int sid ///< global side element id w.r.t cutter discretization
  )
  {
    if (averaging_strategy_ == INPAR::XFEM::Xfluid_Sided)
    {
      return MeshCoupling::GetCouplingElement(sid);
    }

    return coupl_dis_->gElement(GetEmbeddedElementId(sid));
  }

  //! get the element from the conditioned dis for a local coupling side element id
  DRT::Element* GetCondElement(
    const int sid ///< global side element id w.r.t cutter discretization
  )
  {
    return cond_dis_->gElement(GetEmbeddedElementId(sid));
  }

  /// set material pointer for coupling slave side
  void GetInterfaceSlaveMaterial(
    DRT::Element* actele,
    Teuchos::RCP<MAT::Material> & mat
  );

  /// set the fluid-fluid interface fix to avoid a cut
  void SetInterfaceFixed()
  {
    // TODO:XFF-class calls this, when used in an FSI algorithm (fixed ALE)
    moving_interface_ = false;
  }

  /// free the fluid-fluid interface
  void SetInterfaceFree()
  {
    moving_interface_ = true;
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

  //! determine whether interface is fixed
  bool HasMovingInterface()
  {
    return moving_interface_;
  }

protected:

  //! Initializes configurationmap
  void InitConfigurationMap();

  //! Updates configurationmap for specific Gausspoint
  virtual void UpdateConfigurationMap_GP( double& kappa_m, double& visc_m, double& visc_s,
      double& visc_stab, double& full_stab, const LINALG::Matrix<3,1>& x, const DRT::Condition* cond);

private:

  //! Access to map from cutting surface element id to embedded element id
  int GetEmbeddedElementId(int sid)
  {
    std::map<int,int>::iterator it = cutting_emb_gid_map_.find(sid);

    if (it == cutting_emb_gid_map_.end())
      dserror("Could not find parent element ID for side %d.", sid);

    return it->second;
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

  //! whether the embedded fluid interface is moving
  bool moving_interface_;
};

} // namespace XFEM

#endif /*XFEM_COUPLING_MESH_H_*/
