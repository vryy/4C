/*----------------------------------------------------------------------*/
/*! \file

\brief manages the different types of mesh based coupling conditions and thereby builds the bridge
between the xfluid class and the cut-library

\level 2

*/
/*----------------------------------------------------------------------*/


#ifndef FOUR_C_XFEM_COUPLING_MESH_HPP
#define FOUR_C_XFEM_COUPLING_MESH_HPP

#include "4C_config.hpp"

#include "4C_xfem_coupling_base.hpp"

FOUR_C_NAMESPACE_OPEN

namespace CORE::GEO
{
  class CutWizard;
}

namespace XFEM
{
  class ConditionManager;
  class XFluidContactComm;
  /*!
  \brief
   */
  class MeshCoupling : public CouplingBase
  {
   public:
    //! constructor
    explicit MeshCoupling(Teuchos::RCP<DRT::Discretization>& bg_dis,  ///< background discretization
        const std::string& cond_name,  ///< name of the condition, by which the derived cutter
                                       ///< discretization is identified
        Teuchos::RCP<DRT::Discretization>&
            cond_dis,           ///< discretization from which cutter discretization can be derived
        const int coupling_id,  ///< id of composite of coupling conditions
        const double time,      ///< time
        const int step,         ///< time step
        const std::string& suffix = "",  ///< suffix for cutterdisname
        bool marked_geometry = false);

    //! get the coupling element (equal to the side for xfluid-sided, mesh-based coupling)
    CORE::Elements::Element* GetCouplingElement(
        const int sid  ///< global side element id w.r.t cutter discretization
        ) override
    {
      return coupl_dis_->gElement(sid);
    }

    /// get the side element of the respective boundary discretization
    CORE::Elements::Element* get_side(
        const int sid  ///< global side element id w.r.t cutter discretization
    )
    {
      return cutter_dis_->gElement(sid);
    }

    Teuchos::RCP<const Epetra_Vector> GetCutterDispCol();

    /// fill lm vector for coupling element
    virtual void get_coupling_ele_location_vector(const int sid, std::vector<int>& patchlm);

    /// set material pointer for coupling slave side
    void get_interface_slave_material(
        CORE::Elements::Element* actele, Teuchos::RCP<CORE::MAT::Material>& mat) override
    {
      mat = Teuchos::null;
    }

    // finalized interface state vectors
    virtual void complete_state_vectors(){};

    // zero interface state vectors for FSI
    virtual void zero_state_vectors_fsi(){};

    /// clear state vectors
    virtual void ClearState();

    /// set state vectors for cutter discretization
    virtual void set_state();

    /// set displacement state vectors for cutter discretization
    virtual void set_state_displacement();

    /// update interface field state vectors
    virtual void UpdateStateVectors();

    /// update last iteration interface displacements
    virtual void update_displacement_iteration_vectors();

    virtual void gmsh_output_discretization(std::ostream& gmshfilecontent);

    virtual void Output(const int step, const double time, const bool write_restart_data);

    void prepare_cutter_output() override;

    virtual void LiftDrag(const int step, const double time) const {};


    virtual void read_restart(const int step){};

    bool HasMovingInterface() override { return true; }

    bool CutGeometry() override { return !mark_geometry_; }
    virtual bool IsMarkedGeometry() { return mark_geometry_; }

    //! do not cut, but only mark part of boundary loaded into cut
    virtual void SetMarkedGeometry(bool markgeometry) { mark_geometry_ = markgeometry; }

    Teuchos::RCP<Epetra_Vector> IVelnp() { return ivelnp_; }
    Teuchos::RCP<Epetra_Vector> IVeln() { return iveln_; }
    Teuchos::RCP<Epetra_Vector> IVelnm() { return ivelnm_; }

    Teuchos::RCP<Epetra_Vector> IDispnp() { return idispnp_; }
    Teuchos::RCP<Epetra_Vector> IDispn() { return idispn_; }

    Teuchos::RCP<Epetra_Vector> IDispnpi() { return idispnpi_; }

    /// Get background fluid mesh h scaling
    virtual double Get_h() { return h_scaling_; }

   protected:
    /*!
     Return a (smoothed -- soon)/non-smoothed tangiential projection of the mesh surface.
     */
    template <CORE::FE::CellType DISTYPE, class T1, class M3>
    void eval_projection_matrix(T1& projection_matrix,  ///< Projection matrix
        M3& normal                                      ///< surface normal of cut element
    )
    {
      // Properties of a projection matrix:
      //-------------------------------------------------------------------------
      // 1) P is singular (i.e. not of full rank, no inverse exists).
      // 2) P*P = P
      // 3) P^T = P
      // 4) a*P*a \geq 0 \forall a
      //-------------------------------------------------------------------------

      // number space dimensions for element
      // const size_t nsd = projection_matrix.Rows();

#ifdef FOUR_C_ENABLE_ASSERTIONS
      if (projection_matrix.numRows() != projection_matrix.numCols() ||
          projection_matrix.numRows() != normal.numRows())
        FOUR_C_THROW(
            "eval_projection_matrix: Rows and Cols of projection_matrix and normal vector do not "
            "fit!");
#endif

      // create projection matrix
      setup_projection_matrix(projection_matrix, normal);

      return;
    }

   private:
    //! create cutting discretization from condition
    virtual void create_cutter_dis_from_condition(std::string suffix);

   protected:
    void set_cutter_discretization() override;

    void set_conditions_to_copy() override;

    void init_state_vectors() override;

    //------------------------------- vectors -----------------------------
    // TODO: these vectors are not required for Neumann and WDBC conditions, derive class
    //! @name cutter-dis state vectors
    Teuchos::RCP<Epetra_Vector> ivelnp_;
    Teuchos::RCP<Epetra_Vector> iveln_;
    Teuchos::RCP<Epetra_Vector> ivelnm_;

    Teuchos::RCP<Epetra_Vector> idispnp_;  ///< current displacements at t^n+1
    Teuchos::RCP<Epetra_Vector> idispn_;   ///< last displacements at t^n
    Teuchos::RCP<Epetra_Vector>
        idispnpi_;  ///< displacements of last Newton increment at t^n+1 (for monolithic approaches)

    bool mark_geometry_;

    //! Background fluid mesh h scaling
    double h_scaling_;
    //@}

    //! @name output discretization writers
    Teuchos::RCP<IO::DiscretizationWriter> cutter_output_;
    bool firstoutputofrun_;

    //@}

    std::string suffix_;
  };

  /*!
  \brief Mesh Coupling Class which handles all the communication between interface discretization
  and volume discretization. This class can be used for all mesh coupling objects as base class if
  the interface requires volume discretization information (e.g. FluidFluid, SlaveSided FSI,
  FPI,...) - (orig. from FluidFluid(Kruse)) \author ager \date 10/16
   */
  class MeshVolCoupling : public MeshCoupling
  {
   public:
    //! constructor
    explicit MeshVolCoupling(
        Teuchos::RCP<DRT::Discretization>& bg_dis,  ///< background discretization
        const std::string& cond_name,  ///< name of the condition, by which the derived cutter
                                       ///< discretization is identified
        Teuchos::RCP<DRT::Discretization>&
            cond_dis,           ///< discretization from which cutter discretization can be derived
        const int coupling_id,  ///< id of composite of coupling conditions
        const double time,      ///< time
        const int step,         ///< time step
        const std::string& suffix = ""  ///< suffix for cutterdisname
    );

    void Init() override;

    //! Initialize Volume Coupling
    void Init_VolCoupling();

    void get_coupling_ele_location_vector(const int sid, std::vector<int>& patchlm) override;

    //! get the coupling element for a local coupling side element id
    CORE::Elements::Element* GetCouplingElement(
        const int sid  ///< global side element id w.r.t cutter discretization
        ) override
    {
      if (get_averaging_strategy() == INPAR::XFEM::Xfluid_Sided)
      {
        return MeshCoupling::GetCouplingElement(sid);
      }
      // else
      CORE::Elements::FaceElement* fele =
          dynamic_cast<CORE::Elements::FaceElement*>(cutter_dis_->gElement(sid));
      if (!fele) FOUR_C_THROW("Cast to FaceElement failed!");
      fele->set_parent_master_element(
          coupl_dis_->gElement(fele->ParentElementId()), fele->FaceParentNumber());
      return fele->parent_element();
    }

    //! get the element from the conditioned dis for a local coupling side element id
    CORE::Elements::Element* GetCondElement(
        const int sid  ///< global side element id w.r.t cutter discretization
    )
    {
      CORE::Elements::FaceElement* fele =
          dynamic_cast<CORE::Elements::FaceElement*>(cutter_dis_->gElement(sid));
      if (!fele) FOUR_C_THROW("Cast to FaceElement failed!");
      return fele->parent_element();
    }

    //! get auxiliary coupling discretization (embedded elements with nodes in the cutting surface
    //! discretization)
    Teuchos::RCP<DRT::Discretization> get_auxiliary_discretization()
    {
      FOUR_C_ASSERT(init_volcoupling_,
          "MeshVolCoupling::get_auxiliary_discretization: Volume Coupling not initialized!");
      return aux_coup_dis_;
    }

    //! reset all evaluated trace estimates, next time the are required will be calculated again!
    void reset_evaluated_trace_estimates();

    //! get the estimation of the penalty scaling in Nitsche's method from the trace inequality for
    //! a specific face element via solving a local eigenvalue problem
    double get_estimate_nitsche_trace_max_eigenvalue(CORE::Elements::Element* ele);

   protected:
    //! ghost embedded elements, that contribute to the cutting interface discretization on all
    //! procs
    void redistribute_embedded_discretization();

    //! estimate the penalty scaling in Nitsche's method from the trace inequality for a specific
    //! face element via solving a local eigenvalue problem
    virtual void estimate_nitsche_trace_max_eigenvalue(CORE::Elements::Element* ele)
    {
      FOUR_C_THROW(
          "estimate_nitsche_trace_max_eigenvalue not implemented for your coupling object!");
    }

    //! build an auxiliary discretization out of the elements, that contribute to the cutting
    //! discretization
    void create_auxiliary_discretization();

    //! map stores max eigenvalues of trace estimate of the elements
    Teuchos::RCP<std::map<int, double>> ele_to_max_eigenvalue_;

   private:
    //! auxiliary discretization, based on the elements of cond_dis, that contribute to
    //! the elements of cutter_dis_ with at least one edge
    Teuchos::RCP<DRT::Discretization> aux_coup_dis_;

    //! bool to indicate if volume coupling is initialized
    bool init_volcoupling_;

    //! when should the local eigenvalue problem be updated
    INPAR::XFEM::TraceEstimateEigenvalueUpdate trace_estimate_eigenvalue_update_;

    //! last reset of local eigenvalue problem
    int reset_step_;
  };

  /*!
  \brief
   */
  class MeshCouplingBC : public MeshCoupling
  {
   public:
    //! constructor
    explicit MeshCouplingBC(
        Teuchos::RCP<DRT::Discretization>& bg_dis,  ///< background discretization
        const std::string& cond_name,  ///< name of the condition, by which the derived cutter
                                       ///< discretization is identified
        Teuchos::RCP<DRT::Discretization>&
            cond_dis,           ///< discretization from which cutter discretization can be derived
        const int coupling_id,  ///< id of composite of coupling conditions
        const double time,      ///< time
        const int step,         ///< time step
        bool marked_geometry = false  ///< is this a marked geometry mesh boundary
    );


   private:
    void evaluate_interface_displacement(std::vector<double>& final_values, DRT::Node* node,
        CORE::Conditions::Condition* cond, const double time);

    void evaluate_interface_velocity(std::vector<double>& final_values, DRT::Node* node,
        CORE::Conditions::Condition* cond, const double time, const double dt);

    void compute_interface_velocity_from_displacement(std::vector<double>& final_values,
        DRT::Node* node, const double dt, const std::string* evaltype);

    void evaluate_implementation(std::vector<double>& final_values, const double* x,
        CORE::Conditions::Condition* cond, const double time, const std::string& function_name);

   protected:
    void do_condition_specific_setup() override;

    virtual void set_interface_displacement();

    virtual void set_interface_velocity();

    virtual void evaluate_condition(Teuchos::RCP<Epetra_Vector> ivec, const std::string& condname,
        const double time, const double dt = 0.0);

    bool HasMovingInterface() override;
  };

  /*!
  \brief
   */
  class MeshCouplingWeakDirichlet : public MeshCouplingBC
  {
   public:
    //! constructor
    explicit MeshCouplingWeakDirichlet(
        Teuchos::RCP<DRT::Discretization>& bg_dis,  ///< background discretization
        const std::string& cond_name,  ///< name of the condition, by which the derived cutter
                                       ///< discretization is identified
        Teuchos::RCP<DRT::Discretization>&
            cond_dis,           ///< discretization from which cutter discretization can be derived
        const int coupling_id,  ///< id of composite of coupling conditions
        const double time,      ///< time
        const int step,         ///< time step
        bool marked_geometry = false  ///< is this a marked geometry mesh boundary
    );

   public:
    void evaluate_coupling_conditions(CORE::LINALG::Matrix<3, 1>& ivel,
        CORE::LINALG::Matrix<3, 1>& itraction, const CORE::LINALG::Matrix<3, 1>& x,
        const CORE::Conditions::Condition* cond) override;

    void evaluate_coupling_conditions_old_state(CORE::LINALG::Matrix<3, 1>& ivel,
        CORE::LINALG::Matrix<3, 1>& itraction, const CORE::LINALG::Matrix<3, 1>& x,
        const CORE::Conditions::Condition* cond) override;

    void PrepareSolve() override;

   protected:
    void do_condition_specific_setup() override;

    //! Initializes configurationmap
    void setup_configuration_map() override;

    //! Updates configurationmap for specific Gausspoint
    void update_configuration_map_gp(double& kappa_m,  //< fluid sided weighting
        double& visc_m,                                //< master sided dynamic viscosity
        double& visc_s,                                //< slave sided dynamic viscosity
        double& density_m,                             //< master sided density
        double& visc_stab_tang,                        //< viscous tangential NIT Penalty scaling
        double& full_stab,                             //< full NIT Penalty scaling
        const CORE::LINALG::Matrix<3, 1>& x,           //< Position x in global coordinates
        const CORE::Conditions::Condition* cond,       //< Condition
        CORE::Elements::Element* ele,                  //< Element
        CORE::Elements::Element* bele,                 //< Boundary Element
        double* funct,  //< local shape function for Gauss Point (from fluid element)
        double* derxy,  //< local derivatives of shape function for Gauss Point (from fluid element)
        CORE::LINALG::Matrix<3, 1>& rst_slave,  //< local coord of gp on slave boundary element
        CORE::LINALG::Matrix<3, 1>& normal,     //< normal at gp
        CORE::LINALG::Matrix<3, 1>& vel_m,      //< master velocity at gp
        double* fulltraction  //< precomputed fsi traction (sigmaF n + gamma relvel)
        ) override;
  };

  /*!
  \brief
   */
  class MeshCouplingNeumann : public MeshCouplingBC
  {
   public:
    //! constructor
    explicit MeshCouplingNeumann(
        Teuchos::RCP<DRT::Discretization>& bg_dis,  ///< background discretization
        const std::string& cond_name,  ///< name of the condition, by which the derived cutter
                                       ///< discretization is identified
        Teuchos::RCP<DRT::Discretization>&
            cond_dis,           ///< discretization from which cutter discretization can be derived
        const int coupling_id,  ///< id of composite of coupling conditions
        const double time,      ///< time
        const int step,         ///< time step
        bool marked_geometry = false  ///< is this a marked geometry mesh boundary
        )
        : MeshCouplingBC(bg_dis, cond_name, cond_dis, coupling_id, time, step, marked_geometry),
          inflow_stab_(false)
    {
    }

   public:
    //! Evaluate Neumann traction 3 components
    void evaluate_coupling_conditions(CORE::LINALG::Matrix<3, 1>& ivel,
        CORE::LINALG::Matrix<3, 1>& itraction, const CORE::LINALG::Matrix<3, 1>& x,
        const CORE::Conditions::Condition* cond) override;

    //! Evaluate Neumann traction 6 components
    void evaluate_coupling_conditions(CORE::LINALG::Matrix<3, 1>& ivel,
        CORE::LINALG::Matrix<6, 1>& itraction, const CORE::LINALG::Matrix<3, 1>& x,
        const CORE::Conditions::Condition* cond) override;

    void evaluate_coupling_conditions_old_state(CORE::LINALG::Matrix<3, 1>& ivel,
        CORE::LINALG::Matrix<3, 1>& itraction, const CORE::LINALG::Matrix<3, 1>& x,
        const CORE::Conditions::Condition* cond) override;

    void PrepareSolve() override;

   protected:
    //! Do condition specific setup
    void do_condition_specific_setup() override;

    //! Initializes configurationmap
    void setup_configuration_map() override;

    //! Updates configurationmap for specific Gausspoint
    void update_configuration_map_gp(double& kappa_m,  //< fluid sided weighting
        double& visc_m,                                //< master sided dynamic viscosity
        double& visc_s,                                //< slave sided dynamic viscosity
        double& density_m,                             //< master sided density
        double& visc_stab_tang,                        //< viscous tangential NIT Penalty scaling
        double& full_stab,                             //< full NIT Penalty scaling
        const CORE::LINALG::Matrix<3, 1>& x,           //< Position x in global coordinates
        const CORE::Conditions::Condition* cond,       //< Condition
        CORE::Elements::Element* ele,                  //< Element
        CORE::Elements::Element* bele,                 //< Boundary Element
        double* funct,  //< local shape function for Gauss Point (from fluid element)
        double* derxy,  //< local derivatives of shape function for Gauss Point (from fluid element)
        CORE::LINALG::Matrix<3, 1>& rst_slave,  //< local coord of gp on slave boundary element
        CORE::LINALG::Matrix<3, 1>& normal,     //< normal at gp
        CORE::LINALG::Matrix<3, 1>& vel_m,      //< master velocity at gp
        double* fulltraction  //< precomputed fsi traction (sigmaF n + gamma relvel)
        ) override;

   private:
    //! Flag for inflow stabilization
    bool inflow_stab_;
  };


  /*!
  \brief
   */
  class MeshCouplingNavierSlip : public MeshCouplingBC
  {
   public:
    //! constructor
    explicit MeshCouplingNavierSlip(
        Teuchos::RCP<DRT::Discretization>& bg_dis,  ///< background discretization
        const std::string& cond_name,  ///< name of the condition, by which the derived cutter
                                       ///< discretization is identified
        Teuchos::RCP<DRT::Discretization>&
            cond_dis,           ///< discretization from which cutter discretization can be derived
        const int coupling_id,  ///< id of composite of coupling conditions
        const double time,      ///< time
        const int step,         ///< time step
        bool marked_geometry = false  ///< is this a marked geometry mesh boundary
    );

   public:
    void evaluate_coupling_conditions(CORE::LINALG::Matrix<3, 1>& ivel,
        CORE::LINALG::Matrix<3, 1>& itraction, CORE::LINALG::Matrix<3, 3>& proj_matrix,
        const CORE::LINALG::Matrix<3, 1>& x, const CORE::LINALG::Matrix<3, 1>& normal,
        const CORE::Conditions::Condition* cond, const bool& eval_dirich_at_gp,
        double& kappa_m,  ///< fluid sided weighting
        double& visc_m,   ///< fluid sided weighting
        double& visc_s    ///< slave sided dynamic viscosity
    );

    /// this function has to be reviewed for usage with OST new.
    void evaluate_coupling_conditions_old_state(CORE::LINALG::Matrix<3, 1>& ivel,
        CORE::LINALG::Matrix<3, 1>& itraction, const CORE::LINALG::Matrix<3, 1>& x,
        const CORE::Conditions::Condition* cond) override;

    /// get the slip coefficient for this coupling
    void GetSlipCoefficient(double& slipcoeff, const CORE::LINALG::Matrix<3, 1>& x,
        const CORE::Conditions::Condition* cond) override;

    void PrepareSolve() override;

   protected:
    void set_condition_specific_parameters() override;

    void do_condition_specific_setup() override;

    void set_interface_velocity() override;

    void get_condition_by_robin_id(const std::vector<CORE::Conditions::Condition*>& mycond,
        const int coupling_id, std::vector<CORE::Conditions::Condition*>& mynewcond);

    void create_robin_id_map(const std::vector<CORE::Conditions::Condition*>& conditions_NS,
        const std::vector<CORE::Conditions::Condition*>& conditions_robin,
        const std::string& robin_id_name,
        std::map<int, CORE::Conditions::Condition*>& conditionsmap_robin);

    //! Initializes configurationmap
    void setup_configuration_map() override;

    //! Updates configurationmap for specific Gausspoint
    void update_configuration_map_gp(double& kappa_m,  //< fluid sided weighting
        double& visc_m,                                //< master sided dynamic viscosity
        double& visc_s,                                //< slave sided dynamic viscosity
        double& density_m,                             //< master sided density
        double& visc_stab_tang,                        //< viscous tangential NIT Penalty scaling
        double& full_stab,                             //< full NIT Penalty scaling
        const CORE::LINALG::Matrix<3, 1>& x,           //< Position x in global coordinates
        const CORE::Conditions::Condition* cond,       //< Condition
        CORE::Elements::Element* ele,                  //< Element
        CORE::Elements::Element* bele,                 //< Boundary Element
        double* funct,  //< local shape function for Gauss Point (from fluid element)
        double* derxy,  //< local derivatives of shape function for Gauss Point (from fluid element)
        CORE::LINALG::Matrix<3, 1>& rst_slave,  //< local coord of gp on slave boundary element
        CORE::LINALG::Matrix<3, 1>& normal,     //< normal at gp
        CORE::LINALG::Matrix<3, 1>& vel_m,      //< master velocity at gp
        double* fulltraction  //< precomputed fsi traction (sigmaF n + gamma relvel)
        ) override;

   protected:
    // Condition map. Corelating robin_id from Navier Slip condition and
    //                 Robin Dirichlet/Neumann input sections
    //        robin_id    Robin_cond
    std::map<int, CORE::Conditions::Condition*> conditionsmap_robin_dirch_;
    std::map<int, CORE::Conditions::Condition*> conditionsmap_robin_neumann_;

    //       Coupling Surface (E 1)          slip length    is slip length constant?
    std::map<int, std::pair<double, bool>> sliplength_map_;
    //       Coupling Surface (E 1)  Force only tangential veloctiy to surface
    std::map<int, bool> force_tangvel_map_;
  };



  /*!
  \brief
   */
  class MeshCouplingFSI : public MeshVolCoupling
  {
   public:
    //! constructor
    explicit MeshCouplingFSI(
        Teuchos::RCP<DRT::Discretization>& bg_dis,  ///< background discretization
        const std::string& cond_name,  ///< name of the condition, by which the derived cutter
                                       ///< discretization is identified
        Teuchos::RCP<DRT::Discretization>&
            cond_dis,           ///< discretization from which cutter discretization can be derived
        const int coupling_id,  ///< id of composite of coupling conditions
        const double time,      ///< time
        const int step          ///< time step
    );

    // finalize the interface true residual vector
    void complete_state_vectors() override;

    void zero_state_vectors_fsi() override;

    void GmshOutput(const std::string& filename_base, const int step, const int gmsh_step_diff,
        const bool gmsh_debug_out_screen) override;

    void gmsh_output_discretization(std::ostream& gmshfilecontent) override;

    void LiftDrag(const int step, const double time) const override;

    void read_restart(const int step) override;

    void GetSlipCoefficient(double& slipcoeff, const CORE::LINALG::Matrix<3, 1>& x,
        const CORE::Conditions::Condition* cond) override;

    // interface foces
    Teuchos::RCP<Epetra_Vector> ITrueResidual() { return itrueresidual_; }

    // for assembly of fluid interface forces
    Teuchos::RCP<Epetra_Vector> IForcecol() { return iforcecol_; }

    // evaluate structural cauchy stress and linearization in case we don't have xfluid sided
    // weighting
    void evaluate_structural_cauchy_stress(CORE::Elements::Element* coupl_ele,
        CORE::LINALG::Matrix<3, 1>& rst_slave, std::vector<double>& eledisp,
        const CORE::LINALG::Matrix<3, 1>& normal,
        std::vector<CORE::LINALG::SerialDenseMatrix>& solid_stress);

    void SetTimeFac(double timefac) { timefac_ = timefac; }

    double GetTimeFac() { return timefac_; }

    void get_stress_tangent_slave(CORE::Elements::Element* coup_ele,  ///< solid ele
        double& e_s);                                                 ///< stress tangent slavesided

    /// get scaling of the master side for penalty (viscosity, E-modulus for solids)
    void get_penalty_scaling_slave(CORE::Elements::Element* coup_ele,  ///< xfluid ele
        double& penscaling_s) override  ///< penalty scaling slavesided
    {
      get_stress_tangent_slave(coup_ele, penscaling_s);
    }

    void Output(const int step, const double time, const bool write_restart_data) override;

    /// Assign communicator to contact to mesh coupling object
    void Assign_Contact_Comm(Teuchos::RCP<XFEM::XFluidContactComm> xf_c_comm)
    {
      xf_c_comm_ = xf_c_comm;
    }

    /// Get communicator to contact
    Teuchos::RCP<XFEM::XFluidContactComm> Get_Contact_Comm()
    {
      if (xf_c_comm_ == Teuchos::null)
        FOUR_C_THROW("Get_Contact_Comm: Xfluid_Contact_Communicator not assigned!");
      return xf_c_comm_;
    }

    /// Prepare solve
    void PrepareSolve() override;

    /// Get the corresponding FSI interface law
    virtual INPAR::XFEM::InterfaceLaw GetInterfaceLaw() { return interfacelaw_; }

    /// Register this side on this proc
    void RegisterSideProc(int sid);

    /// Initialize Fluid State
    bool initialize_fluid_state(Teuchos::RCP<CORE::GEO::CutWizard> cutwizard,
        Teuchos::RCP<DRT::Discretization> fluiddis,
        Teuchos::RCP<XFEM::ConditionManager> condition_manager,
        Teuchos::RCP<Teuchos::ParameterList> fluidparams);

   protected:
    //! estimate the penalty scaling in Nitsche's method from the trace inequality for a specific
    //! face element via solving a local eigenvalue problem
    void estimate_nitsche_trace_max_eigenvalue(CORE::Elements::Element* ele) override;

    void init_state_vectors() override;

    bool HasMovingInterface() override { return true; }

    void set_condition_specific_parameters() override;

    //! Initializes configurationmap
    void setup_configuration_map() override;

    //! Updates configurationmap for specific Gausspoint
    void update_configuration_map_gp(double& kappa_m,  //< fluid sided weighting
        double& visc_m,                                //< master sided dynamic viscosity
        double& visc_s,                                //< slave sided dynamic viscosity
        double& density_m,                             //< master sided density
        double& visc_stab_tang,                        //< viscous tangential NIT Penalty scaling
        double& full_stab,                             //< full NIT Penalty scaling
        const CORE::LINALG::Matrix<3, 1>& x,           //< Position x in global coordinates
        const CORE::Conditions::Condition* cond,       //< Condition
        CORE::Elements::Element* ele,                  //< Element
        CORE::Elements::Element* bele,                 //< Boundary Element
        double* funct,  //< local shape function for Gauss Point (from fluid element)
        double* derxy,  //< local derivatives of shape function for Gauss Point (from fluid element)
        CORE::LINALG::Matrix<3, 1>& rst_slave,  //< local coord of gp on slave boundary element
        CORE::LINALG::Matrix<3, 1>& normal,     //< normal at gp
        CORE::LINALG::Matrix<3, 1>& vel_m,      //< master velocity at gp
        double* fulltraction  //< precomputed fsi traction (sigmaF n + gamma relvel)
        ) override;

    //! Updates configurationmap for specific Gausspoint for FSI with contact
    void update_configuration_map_gp_contact(double& kappa_m,  //< fluid sided weighting
        double& visc_m,                                        //< master sided dynamic viscosity
        double& visc_s,                                        //< slave sided dynamic viscosity
        double& density_m,                                     //< master sided density
        double& visc_stab_tang,                   //< viscous tangential NIT Penalty scaling
        double& full_stab,                        //< full NIT Penalty scaling
        const CORE::LINALG::Matrix<3, 1>& x,      //< Position x in global coordinates
        const CORE::Conditions::Condition* cond,  //< Condition
        CORE::Elements::Element* ele,             //< Element
        CORE::Elements::Element* bele,            //< Boundary Element
        double* funct,  //< local shape function for Gauss Point (from fluid element)
        double* derxy,  //< local derivatives of shape function for Gauss Point (from fluid element)
        CORE::LINALG::Matrix<3, 1>& rst_slave,  //< local coord of gp on slave boundary element
        CORE::LINALG::Matrix<3, 1>& normal,     //< normal at gp
        CORE::LINALG::Matrix<3, 1>& vel_m,      //< master velocity at gp
        double* fulltraction  //< precomputed fsi traction (sigmaF n + gamma relvel)
    );

    //------------------------------- vectors -----------------------------
    //! @name cutter-dis state vectors
    Teuchos::RCP<Epetra_Vector>
        itrueresidual_;  //! interface forces acting on the structural surface (= -iforcenp)
    Teuchos::RCP<Epetra_Vector>
        iforcecol_;  //! interface forces acting on the fluid surface (column vector assembly)
    //@}

    //---------------------------------parameters--------------------------------
    //! @name parameters
    //       Coupling Surface (E 1)          slip length    is slip length constant?
    std::map<int, std::pair<double, bool>> sliplength_map_;

    //! theta*timestep
    double timefac_;

    //! applied interface law
    INPAR::XFEM::InterfaceLaw interfacelaw_;

    //! Xfluid Contact Communicator
    Teuchos::RCP<XFEM::XFluidContactComm> xf_c_comm_;
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
        Teuchos::RCP<DRT::Discretization>& bg_dis,  ///< background discretization
        const std::string& cond_name,  ///< name of the condition, by which the derived cutter
                                       ///< discretization is identified
        Teuchos::RCP<DRT::Discretization>&
            cond_dis,           ///< discretization from which cutter discretization can be derived
        const int coupling_id,  ///< id of composite of coupling conditions
        const double time,      ///< time
        const int step          ///< time step
    );

    /// set material pointer for coupling slave side
    void get_interface_slave_material(
        CORE::Elements::Element* actele, Teuchos::RCP<CORE::MAT::Material>& mat) override;

    /// set the fluid-fluid interface fix to avoid a cut
    void SetInterfaceFixed()
    {
      // TODO:XFF-class calls this, when used in an FSI algorithm (fixed ALE)
      moving_interface_ = false;
    }

    /// free the fluid-fluid interface
    void SetInterfaceFree() { moving_interface_ = true; }

    //! ghost interface-contributing embedded elements (required for error calculation in case
    //! of xfluid-sided coupling)
    void redistribute_for_error_calculation();

    //! determine whether interface is fixed
    bool HasMovingInterface() override { return moving_interface_; }

    /// get viscosity of the slave fluid
    void GetViscositySlave(CORE::Elements::Element* coup_ele,  ///< xfluid ele
        double& visc_s                                         ///< viscosity slavesided
    );

    /// get scaling of the master side for penalty (viscosity, E-modulus for solids)
    void get_penalty_scaling_slave(CORE::Elements::Element* coup_ele,  ///< xfluid ele
        double& penscaling_s) override  ///< penalty scaling slavesided
    {
      GetViscositySlave(coup_ele, penscaling_s);
    }


    void read_restart(const int step) override;

    void Output(const int step, const double time, const bool write_restart_data) override;

   protected:
    //! estimate the penalty scaling in Nitsche's method from the trace inequality for a specific
    //! face element via solving a local eigenvalue problem
    void estimate_nitsche_trace_max_eigenvalue(CORE::Elements::Element* ele) override;

    //! Initializes configurationmap
    void setup_configuration_map() override;

    //! Updates configurationmap for specific Gausspoint
    void update_configuration_map_gp(double& kappa_m,  //< fluid sided weighting
        double& visc_m,                                //< master sided dynamic viscosity
        double& visc_s,                                //< slave sided dynamic viscosity
        double& density_m,                             //< master sided density
        double& visc_stab_tang,                        //< viscous tangential NIT Penalty scaling
        double& full_stab,                             //< full NIT Penalty scaling
        const CORE::LINALG::Matrix<3, 1>& x,           //< Position x in global coordinates
        const CORE::Conditions::Condition* cond,       //< Condition
        CORE::Elements::Element* ele,                  //< Element
        CORE::Elements::Element* bele,                 //< Boundary Element
        double* funct,  //< local shape function for Gauss Point (from fluid element)
        double* derxy,  //< local derivatives of shape function for Gauss Point (from fluid element)
        CORE::LINALG::Matrix<3, 1>& rst_slave,  //< local coord of gp on slave boundary element
        CORE::LINALG::Matrix<3, 1>& normal,     //< normal at gp
        CORE::LINALG::Matrix<3, 1>& vel_m,      //< master velocity at gp
        double* fulltraction  //< precomputed fsi traction (sigmaF n + gamma relvel)
        ) override;

   private:
    //! whether the embedded fluid interface is moving
    bool moving_interface_;
  };

}  // namespace XFEM

FOUR_C_NAMESPACE_CLOSE

#endif
