/*----------------------------------------------------------------------*/
/*!
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
*/
/*----------------------------------------------------------------------*/


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
      const int                           step,      ///< time step
      const std::string &                 suffix = ""///< suffix for cutterdisname
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

protected:

  /*!
   Return a (smoothed -- soon)/non-smoothed tangiential projection of the mesh surface.
   */
  template <DRT::Element::DiscretizationType DISTYPE, class T1, class M3>
  void EvalProjectionMatrix(
      T1&       projection_matrix,         ///< Projection matrix
      M3&       normal                     ///< surface normal of cut element
  )
  {
    //Properties of a projection matrix:
    //-------------------------------------------------------------------------
    // 1) P is singular (i.e. not of full rank, no inverse exists).
    // 2) P*P = P
    // 3) P^T = P
    // 4) a*P*a \geq 0 \forall a
    //-------------------------------------------------------------------------

    // number space dimensions for element
    const size_t nsd = projection_matrix.Rows();

    #ifdef DEBUG
        if (projection_matrix.Rows() != projection_matrix.Cols() || projection_matrix.Rows() != normal.Rows())
          dserror("EvalProjectionMatrix: Rows and Cols of projection_matrix and normal vector do not fit!");
    #endif


      //create projection matrix
   // if(projtosurf_== INPAR::XFEM::Proj_normal)
    {
      for(unsigned int i =0; i<nsd; ++i)
      {
        for(unsigned int j =0; j<nsd; ++j)
        {
          projection_matrix(i,j) = (i==j) - normal(i,0) * normal(j,0);
        }
      }
    }

    return;
  }

private:


  void SetConditionsToCopy();

  //! create cutting discretization from condition
  virtual void CreateCutterDisFromCondition(std::string suffix);

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
\brief Mesh Coupling Class which handles all the communication between interface discretization and volume discretization.
This class can be used for all mesh coupling objects as base class if the interface requires volume discretization information
(e.g. FluidFluid, SlaveSided FSI, FPI,...) - (orig. from FluidFluid(Kruse))
\author ager
\date 10/16
 */
class MeshVolCoupling : public MeshCoupling
{
public:

  //! constructor
  explicit MeshVolCoupling(
      Teuchos::RCP<DRT::Discretization>&  bg_dis,                      ///< background discretization
      const std::string &                 cond_name,                   ///< name of the condition, by which the derived cutter discretization is identified
      Teuchos::RCP<DRT::Discretization>&  cond_dis,                    ///< discretization from which cutter discretization can be derived
      const int                           coupling_id,                 ///< id of composite of coupling conditions
      const double                        time,                        ///< time
      const int                           step,                        ///< time step
      const std::string &                 suffix = ""                  ///< suffix for cutterdisname
  );

  //! destructor
  virtual ~MeshVolCoupling()
  {    return;  }

  //! Initialize Volume Coupling
  void Init_VolCoupling();

  void GetCouplingEleLocationVector(
    const int sid,
    std::vector<int> & patchlm);

  //! get the coupling element for a local coupling side element id
  DRT::Element* GetCouplingElement(
    const int sid ///< global side element id w.r.t cutter discretization
  )
  {
    if (GetAveragingStrategy() == INPAR::XFEM::Xfluid_Sided)
    {
      return MeshCoupling::GetCouplingElement(sid);
    }
    //else
    DRT::FaceElement* fele = dynamic_cast<DRT::FaceElement*>(cutter_dis_->gElement(sid));
    if (!fele) dserror("Cast to FaceElement failed!");
    return fele->ParentElement();
  }

  //! get the element from the conditioned dis for a local coupling side element id
  DRT::Element* GetCondElement(
    const int sid ///< global side element id w.r.t cutter discretization
  )
  {
    DRT::FaceElement* fele = dynamic_cast<DRT::FaceElement*>(cutter_dis_->gElement(sid));
    if (!fele) dserror("Cast to FaceElement failed!");
    return fele->ParentElement();
  }

  //! get auxiliary coupling discretization (embedded elements with nodes in the cutting surface discretization)
  Teuchos::RCP<DRT::Discretization> GetAuxiliaryDiscretization()
  {
    dsassert(init_volcoupling_,"MeshVolCoupling::GetAuxiliaryDiscretization: Volume Coupling not initialized!");
    return aux_coup_dis_;
  }

  //! reset all evaluated trace estimates, next time the are required will be calculated again!
  void ResetEvaluatedTraceEstimates()
  { ele_to_max_eigenvalue_->clear(); }

  //! get the estimation of the penalty scaling in Nitsche's method from the trace inequality for a specific face element
  //! via solving a local eigenvalue problem
  double Get_EstimateNitscheTraceMaxEigenvalue(DRT::Element* ele);

protected:

  //! ghost embedded elements, that contribute to the cutting interface discretization on all procs
  void RedistributeEmbeddedDiscretization();

  //! estimate the penalty scaling in Nitsche's method from the trace inequality for a specific face element
  //! via solving a local eigenvalue problem
  virtual void EstimateNitscheTraceMaxEigenvalue(DRT::Element* ele)
  { dserror("EstimateNitscheTraceMaxEigenvalue not implemented for your coupling object!"); }

  //! build an auxiliary discretization out of the elements, that contribute to the cutting discretization
  void CreateAuxiliaryDiscretization();

  //! map stores max eigenvalues of trace estimate of the elements
  Teuchos::RCP<std::map<int,double> > ele_to_max_eigenvalue_;

private:

  //! auxiliary discretization, based on the elements of cond_dis, that contribute to
  //! the elements of cutter_dis_ with at least one edge
  Teuchos::RCP<DRT::Discretization> aux_coup_dis_;

  //! bool to indicate if volume coupling is initialized
  bool init_volcoupling_;

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
class MeshCouplingFSI : public MeshVolCoupling
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


/*!
\brief specialized class for coupling with an embedded fluid discretization
  \author kruse
  \date 01/15
 */
class MeshCouplingFluidFluid : public MeshVolCoupling
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

  //! ghost interface-contributing embedded elements (required for error calculation in case
  //! of xfluid-sided coupling)
  void RedistributeForErrorCalculation();

  //! determine whether interface is fixed
  bool HasMovingInterface()
  {
    return moving_interface_;
  }

  /// get viscosity of the slave fluid
  void GetViscositySlave(
      DRT::Element * coup_ele,                   ///< xfluid ele
      double& visc_s                             ///< viscosity slavesided
  );

  /// get scaling of the master side for penalty (viscosity, E-modulus for solids)
  virtual void GetPenaltyScalingSlave(
      DRT::Element * coup_ele,                      ///< xfluid ele
      double& penscaling_s)                         ///< penalty scaling slavesided
  {
    GetViscositySlave(coup_ele,penscaling_s);
  }

protected:

  //! estimate the penalty scaling in Nitsche's method from the trace inequality for a specific face element
  //! via solving a local eigenvalue problem
  virtual void EstimateNitscheTraceMaxEigenvalue(DRT::Element* ele);

  //! Initializes configurationmap
  void InitConfigurationMap();

  //! Updates configurationmap for specific Gausspoint
  virtual void UpdateConfigurationMap_GP( double& kappa_m, double& visc_m, double& visc_s,
      double& visc_stab, double& full_stab, const LINALG::Matrix<3,1>& x, const DRT::Condition* cond);

private:

  //! whether the embedded fluid interface is moving
  bool moving_interface_;
};

} // namespace XFEM

#endif /*XFEM_COUPLING_MESH_H_*/
