/*----------------------------------------------------------------------*/
/*! \file

\brief manages the different types of level-set based coupling conditions and thereby builds the
bridge between the xfluid class and the cut-library

\level 2

*/
/*----------------------------------------------------------------------*/


#ifndef BACI_XFEM_COUPLING_LEVELSET_HPP
#define BACI_XFEM_COUPLING_LEVELSET_HPP

#include "baci_config.hpp"

#include "baci_cut_point.hpp"
#include "baci_discretization_fem_general_utils_local_connectivity_matrices.hpp"
#include "baci_xfem_coupling_base.hpp"
#include "baci_xfem_utils.hpp"

BACI_NAMESPACE_OPEN

namespace CORE::GEO
{
  class CutWizard;
}

namespace XFEM
{
  /*!
  \brief
   */
  class LevelSetCoupling : public CouplingBase
  {
   public:
    //! constructor
    explicit LevelSetCoupling(
        Teuchos::RCP<DRT::Discretization>& bg_dis,  ///< background discretization
        const std::string& cond_name,  ///< name of the condition, by which the derived cutter
                                       ///< discretization is identified
        Teuchos::RCP<DRT::Discretization>&
            cond_dis,  ///< full discretization from which the cutter discretization is derived
        const int coupling_id,  ///< id of composite of coupling conditions
        const double time,      ///< time
        const int step          ///< time step
    );

    void SetCouplingDofsets() override;

    bool HaveMatchingNodes(const Teuchos::RCP<DRT::Discretization>& dis_A,
        const Teuchos::RCP<DRT::Discretization>& dis_B);

    void MapCutterToBgVector(const Teuchos::RCP<DRT::Discretization>& source_dis,
        const Teuchos::RCP<Epetra_Vector>& source_vec_dofbased, const int source_nds,
        const Teuchos::RCP<DRT::Discretization>& target_dis,
        const Teuchos::RCP<Epetra_Vector>& target_vec_dofbased, const int target_nds);

    // TODO: sort the functions...

    void SetCutterDiscretization() override;

    void SetConditionSpecificParameters() override{};

    void PrepareCutterOutput() override;

    void DoConditionSpecificSetup() override;


    /// set levelset field by function, return if interface moved compared to last time step
    bool SetLevelSetField(const double time);

    /// initialize level set based state vectors
    void InitStateVectors() override;

    virtual void InitStateVectors_Bg();

    virtual void InitStateVectors_Cutter();

    /// set level-boolean type
    virtual void SetLevelSetBooleanType();

    virtual bool ApplyComplementaryOperator();

    virtual void Output(
        const int step, const double time, const bool write_restart_data, const int lsc_idx = 0);

    void GmshOutput(const std::string& filename_base, const int step, const int gmsh_step_diff,
        const bool gmsh_debug_out_screen) override;

    Teuchos::RCP<Epetra_Vector> GetLevelSetFieldAsNodeRowVector();

    virtual void ReadRestart(const int step, const int lsc_idx = 0);

    bool HasMovingInterface() override { return true; }

    void GetInterfaceSlaveMaterial(DRT::Element* actele, Teuchos::RCP<MAT::Material>& mat) override
    {
      mat = Teuchos::null;
    }

    XFEM::CouplingBase::LevelSetBooleanType GetBooleanCombination() { return ls_boolean_type_; }

    //! export row vectors storing geometric quantities to col vectors
    virtual void ExportGeometricQuantities(){};


   private:
    void SetConditionsToCopy() override;

    /// set level-set field implemented in this routine
    double FunctImplementation(const int func_no, const double* coords, const double t);

   protected:
    //! Output specific
    Teuchos::RCP<IO::DiscretizationWriter> bg_output_;

    //! @name fluid discretization related state vectors

    //! fluid-dis (bgdis) state vectors for levelset applications
    Teuchos::RCP<Epetra_Vector> phinp_;


    //@}

    //! @name scatra discretization related state vectors

    //! scatra-dis (cutterdis) state vectors for levelset applications, prepares nonmatching
    //! discretizations between scatra and fluid
    Teuchos::RCP<Epetra_Vector> cutter_phinp_;
    Teuchos::RCP<Epetra_Vector> cutter_phinp_col_;

    //! The nodal curvature and smoothed gradient of the levelset field. (Stored w.r.t to the
    //! scatra-dis = cutter-dis)
    Teuchos::RCP<Epetra_Vector> curvaturenp_node_;
    Teuchos::RCP<Epetra_MultiVector> gradphinp_smoothed_node_;
    // Teuchos::RCP<Epetra_MultiVector>   gradphi2np_smoothed_node_;

    //! and column versions
    Teuchos::RCP<Epetra_Vector> curvaturenp_node_col_;
    Teuchos::RCP<Epetra_MultiVector> gradphinp_smoothed_node_col_;

    //! boolean operation type on level-set for current ls-field and previous combination of
    //! level-set fields
    XFEM::CouplingBase::LevelSetBooleanType ls_boolean_type_;


    // Specify way of creating the projection matrix
    INPAR::XFEM::ProjToSurface projtosurf_;

    //@}

    int bg_nds_phi_;      ///<
    int cutter_nds_phi_;  ///<

    double normal_orientation_;  ///< correction factor between normal of phi-gradient and normal in
                                 ///< xfluid

    bool have_nodematching_dis_;  ///< are bgdis and cutterdis node-matching?
  };

  class LevelSetCouplingBC : public LevelSetCoupling
  {
   public:
    //! constructor
    explicit LevelSetCouplingBC(
        Teuchos::RCP<DRT::Discretization>& bg_dis,  ///< background discretization
        const std::string& cond_name,  ///< name of the condition, by which the derived cutter
                                       ///< discretization is identified
        Teuchos::RCP<DRT::Discretization>&
            cond_dis,  ///< full discretization from which the cutter discretization is derived
        const int coupling_id,  ///< id of composite of coupling conditions
        const double time,      ///< time
        const int step          ///< time step
    );


    void PrepareSolve() override;

    bool HasMovingInterface() override;

   protected:
    bool has_interface_moved_;  ///< did interface move compared to the last time step?
  };


  /*!
  \brief
   */
  class LevelSetCouplingWeakDirichlet : public LevelSetCouplingBC
  {
   public:
    //! constructor
    explicit LevelSetCouplingWeakDirichlet(
        Teuchos::RCP<DRT::Discretization>& bg_dis,  ///< background discretization
        const std::string& cond_name,  ///< name of the condition, by which the derived cutter
                                       ///< discretization is identified
        Teuchos::RCP<DRT::Discretization>&
            cond_dis,  ///< full discretization from which the cutter discretization is derived
        const int coupling_id,  ///< id of composite of coupling conditions
        const double time,      ///< time
        const int step          ///< time step
        )
        : LevelSetCouplingBC(bg_dis, cond_name, cond_dis, coupling_id, time, step)
    {
    }

   public:
    void EvaluateCouplingConditions(CORE::LINALG::Matrix<3, 1>& ivel,
        CORE::LINALG::Matrix<3, 1>& itraction, const CORE::LINALG::Matrix<3, 1>& x,
        const DRT::Condition* cond) override;

    void EvaluateCouplingConditionsOldState(CORE::LINALG::Matrix<3, 1>& ivel,
        CORE::LINALG::Matrix<3, 1>& itraction, const CORE::LINALG::Matrix<3, 1>& x,
        const DRT::Condition* cond) override;

   protected:
    //! Initializes configurationmap
    void SetupConfigurationMap() override;

    //! Updates configurationmap for specific Gausspoint
    void UpdateConfigurationMap_GP(double& kappa_m,  //< fluid sided weighting
        double& visc_m,                              //< master sided dynamic viscosity
        double& visc_s,                              //< slave sided dynamic viscosity
        double& density_m,                           //< master sided density
        double& visc_stab_tang,                      //< viscous tangential NIT Penalty scaling
        double& full_stab,                           //< full NIT Penalty scaling
        const CORE::LINALG::Matrix<3, 1>& x,         //< Position x in global coordinates
        const DRT::Condition* cond,                  //< Condition
        DRT::Element* ele,                           //< Element
        DRT::Element* bele,                          //< Boundary Element
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
  class LevelSetCouplingNeumann : public LevelSetCouplingBC
  {
   public:
    //! constructor
    explicit LevelSetCouplingNeumann(
        Teuchos::RCP<DRT::Discretization>& bg_dis,  ///< background discretization
        const std::string& cond_name,  ///< name of the condition, by which the derived cutter
                                       ///< discretization is identified
        Teuchos::RCP<DRT::Discretization>&
            cond_dis,  ///< full discretization from which the cutter discretization is derived
        const int coupling_id,  ///< id of composite of coupling conditions
        const double time,      ///< time
        const int step          ///< time step
        )
        : LevelSetCouplingBC(bg_dis, cond_name, cond_dis, coupling_id, time, step),
          inflow_stab_(false)
    {
    }

   public:
    //! Evaluate Neumann traction 3 components
    void EvaluateCouplingConditions(CORE::LINALG::Matrix<3, 1>& ivel,
        CORE::LINALG::Matrix<3, 1>& itraction, const CORE::LINALG::Matrix<3, 1>& x,
        const DRT::Condition* cond) override;

    //! Evaluate Neumann traction 6 components
    void EvaluateCouplingConditions(CORE::LINALG::Matrix<3, 1>& ivel,
        CORE::LINALG::Matrix<6, 1>& itraction, const CORE::LINALG::Matrix<3, 1>& x,
        const DRT::Condition* cond) override;

    void EvaluateCouplingConditionsOldState(CORE::LINALG::Matrix<3, 1>& ivel,
        CORE::LINALG::Matrix<3, 1>& itraction, const CORE::LINALG::Matrix<3, 1>& x,
        const DRT::Condition* cond) override;

   protected:
    //! Do condition specific setup
    void DoConditionSpecificSetup() override;

    //! Initializes configurationmap
    void SetupConfigurationMap() override;

    //! Updates configurationmap for specific Gausspoint
    void UpdateConfigurationMap_GP(double& kappa_m,  //< fluid sided weighting
        double& visc_m,                              //< master sided dynamic viscosity
        double& visc_s,                              //< slave sided dynamic viscosity
        double& density_m,                           //< master sided density
        double& visc_stab_tang,                      //< viscous tangential NIT Penalty scaling
        double& full_stab,                           //< full NIT Penalty scaling
        const CORE::LINALG::Matrix<3, 1>& x,         //< Position x in global coordinates
        const DRT::Condition* cond,                  //< Condition
        DRT::Element* ele,                           //< Element
        DRT::Element* bele,                          //< Boundary Element
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
  class LevelSetCouplingNavierSlip : public LevelSetCouplingBC
  {
   public:
    //! constructor
    explicit LevelSetCouplingNavierSlip(
        Teuchos::RCP<DRT::Discretization>& bg_dis,  ///< background discretization
        const std::string& cond_name,  ///< name of the condition, by which the derived cutter
                                       ///< discretization is identified
        Teuchos::RCP<DRT::Discretization>&
            cond_dis,  ///< full discretization from which the cutter discretization is derived
        const int coupling_id,  ///< id of composite of coupling conditions
        const double time,      ///< time
        const int step          ///< time step
    );

   public:
    void EvaluateCouplingConditions(CORE::LINALG::Matrix<3, 1>& ivel,
        CORE::LINALG::Matrix<3, 1>& itraction, const CORE::LINALG::Matrix<3, 1>& x,
        const DRT::Condition* cond) override;

    void EvaluateCouplingConditionsOldState(CORE::LINALG::Matrix<3, 1>& ivel,
        CORE::LINALG::Matrix<3, 1>& itraction, const CORE::LINALG::Matrix<3, 1>& x,
        const DRT::Condition* cond) override;

    void GetSlipCoefficient(double& slipcoeff, const CORE::LINALG::Matrix<3, 1>& x,
        const DRT::Condition* cond) override;


    /*!
     Return prescribed velocities and traction vectors for a GNBC boundary condition.
     Also returns the projection matrix (to the plane of the surface) needed for the GNBC condition.
     */
    template <CORE::FE::CellType DISTYPE, class V1, class V2, class X1, class T1, class M1,
        class M2, class M3>
    void EvaluateCouplingConditions(V1& ivel,  ///< prescribed velocity at interface
        V2& itraction,                         ///< prescribed traction at interface
        X1& x,                                 ///< coordinates of gauss point
        const DRT::Condition* cond,            ///< condition prescribed to this surface
        T1& projection_matrix,  ///< Laplace-Beltrami matrix for surface tension calculations
        int eid,                ///< element ID
        M1& funct,              ///< local shape function for Gauss Point (from fluid element)
        M2& derxy,   ///< local derivatives of shape function for Gauss Point (from fluid element)
        M3& normal,  ///< surface normal of cut element
        double& kappa_m,  ///< fluid sided weighting
        double& visc_m,   ///< fluid sided weighting
        double& visc_s    ///< slave sided dynamic viscosity
    )
    {
      EvalProjectionMatrix<DISTYPE>(projection_matrix, eid, funct, derxy, normal);
      EvaluateCouplingConditions(ivel, itraction, x, cond);

      if (has_neumann_jump_)
      {
        // This is maybe not the most efficient implementation as we evaluate dynvisc as well as the
        // sliplenght twice evaluate interface traction (given by Neumann condition) Add this to the
        // veljump!
        double sliplength = 0.0;
        GetSlipCoefficient(sliplength, x, cond);

        if (sliplength < 0.0) dserror("The slip length can not be negative.");

        if (sliplength != 0.0)
        {
          double sl_visc_fac = sliplength / (kappa_m * visc_m + (1.0 - kappa_m) * visc_s);
          V2 tmp_itraction(true);
          tmp_itraction.MultiplyTN(projection_matrix, itraction);
          // Project this into tangential direction!!!
          ivel.Update(sl_visc_fac, tmp_itraction, 1.0);
        }
        itraction.Clear();
      }

      /*Here one could do a projection of ivel to only point in the tangential direction.
        This would enforce that no spurious velocities occur in the normal direction (i.e.
        no-penetration always enforced). However, these will occur in an XFSI. Thus a solution has
        to be found which can handle this the best.
      */

      if (forcetangvel_)
      {
        // We project in the normal direction
        CORE::LINALG::Matrix<3, 1> tmp_ivel(true);
        tmp_ivel.MultiplyTN(
            projection_matrix, ivel);  // apply Projection matrix from the right. (u_0 * P^t)
        ivel.Update(1.0, tmp_ivel, 0.0);
      }
    };

    /*!
     Return a smoothed/non-smoothed tangiential projection of the level set surface.
     */
    template <CORE::FE::CellType DISTYPE, class T1, class M1, class M2, class M3>
    void EvalProjectionMatrix(T1& projection_matrix,  ///< Projection matrix
        int eid,                                      ///< element ID
        M1& funct,  ///< local shape function for Gauss Point (from fluid element)
        M2& derxy,  ///< local derivatives of shape function for Gauss Point (from fluid element)
        M3& normal  ///< surface normal of cut element
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
      const size_t nsd = CORE::FE::dim<DISTYPE>;

      // number of nodes of element
      const size_t nen = CORE::FE::num_nodes<DISTYPE>;

      // Should this be provided as well by the input?
      DRT::Element* actele = cutter_dis_->gElement(eid);

      //   Non-smoothed projection matrix
      CORE::LINALG::Matrix<nsd, 1> gradphi;
      if (projtosurf_ == INPAR::XFEM::Proj_normal)
      {
        gradphi = normal;
      }
      else if (projtosurf_ == INPAR::XFEM::Proj_smoothed)
      {
        // smoothed normal at cutter element nodes, the Gaussian point lies in
        CORE::LINALG::SerialDenseMatrix esmoothedgradphi_test(nsd, nen);
        CORE::LINALG::Matrix<nsd, nen> esmoothedgradphi(esmoothedgradphi_test, View);
        XFEM::UTILS::ExtractQuantityAtElement(esmoothedgradphi_test, actele,
            gradphinp_smoothed_node_col_, cutter_dis_, cutter_nds_phi_, nsd_);

        // Gradients @ GaussPoints
        gradphi.Multiply(esmoothedgradphi, funct);
      }
      else if (projtosurf_ == INPAR::XFEM::Proj_normal_phi)
      {
        CORE::LINALG::SerialDenseMatrix ephi_test(nen, 1);
        CORE::LINALG::Matrix<nen, 1> ephi(ephi_test, View);
        XFEM::UTILS::ExtractQuantityAtElement(
            ephi_test, actele, cutter_phinp_col_, cutter_dis_, cutter_nds_phi_, 1);

        // Gradients @ GaussPoints
        gradphi.Multiply(derxy, ephi);
      }
      else if (projtosurf_ == INPAR::XFEM::Proj_normal_smoothed_comb)
      {
        // smoothed normal at cutter element nodes, the Gaussian point lies in
        CORE::LINALG::SerialDenseMatrix esmoothedgradphi_test(nsd, nen);
        CORE::LINALG::Matrix<nsd, nen> esmoothedgradphi(esmoothedgradphi_test, View);
        XFEM::UTILS::ExtractQuantityAtElement(esmoothedgradphi_test, actele,
            gradphinp_smoothed_node_col_, cutter_dis_, cutter_nds_phi_, nsd_);

        // Gradients @ GaussPoints
        gradphi.Multiply(esmoothedgradphi, funct);

        const double normgradphi = gradphi.Norm2();
        if (normgradphi > 1e-9)  // 1e-9 is set to create a reasonable scaling.
          gradphi.Scale(1.0 / normgradphi);
        else
          gradphi.putScalar(0.0);  // This to catch the cases when gradphi \approx 0

        // normal_comb = alpha_n * normal + (1-alpha)*gradphi
        CORE::LINALG::Matrix<nsd, 1> normal_comb(true);
        double alpha_n = 0.3;
        normal_comb.Update(alpha_n, normal, -(1.0 - alpha_n), gradphi);

        gradphi = normal_comb;
      }
      else
      {
        dserror("This option for a projection matrix %d does not exist. \n", projtosurf_);
      }

      // Normalize the smoothed gradient
      const double normgradphi = gradphi.Norm2();
      if (normgradphi > 1e-9)  // 1e-9 is set to create a reasonable scaling.
        gradphi.Scale(1.0 / normgradphi);
      else
        gradphi.putScalar(0.0);  // This to catch the cases when gradphi \approx 0

      SetupProjectionMatrix(projection_matrix, gradphi);

      return;
    }

   protected:
    void SetElementConditions() override;

    void SetElementSpecificConditions(std::vector<DRT::Condition*>& cutterele_cond,
        const std::string& cond_name, const int& robin_id);

    void SetConditionSpecificParameters() override;

    void GetConditionByRobinId(const std::vector<DRT::Condition*>& mycond, const int coupling_id,
        std::vector<DRT::Condition*>& mynewcond);

    //! Initializes configurationmap
    void SetupConfigurationMap() override;

    //! Updates configurationmap for specific Gausspoint
    void UpdateConfigurationMap_GP(double& kappa_m,  //< fluid sided weighting
        double& visc_m,                              //< master sided dynamic viscosity
        double& visc_s,                              //< slave sided dynamic viscosity
        double& density_m,                           //< master sided density
        double& visc_stab_tang,                      //< viscous tangential NIT Penalty scaling
        double& full_stab,                           //< full NIT Penalty scaling
        const CORE::LINALG::Matrix<3, 1>& x,         //< Position x in global coordinates
        const DRT::Condition* cond,                  //< Condition
        DRT::Element* ele,                           //< Element
        DRT::Element* bele,                          //< Boundary Element
        double* funct,  //< local shape function for Gauss Point (from fluid element)
        double* derxy,  //< local derivatives of shape function for Gauss Point (from fluid element)
        CORE::LINALG::Matrix<3, 1>& rst_slave,  //< local coord of gp on slave boundary element
        CORE::LINALG::Matrix<3, 1>& normal,     //< normal at gp
        CORE::LINALG::Matrix<3, 1>& vel_m,      //< master velocity at gp
        double* fulltraction  //< precomputed fsi traction (sigmaF n + gamma relvel)
        ) override;

   private:
    bool forcetangvel_;
    bool is_constant_sliplength_;
    bool has_neumann_jump_;
    double sliplength_;

    // ID given the Robin Dirichlet/-Neumann conditions (need to match the one given in the "MAIN
    // condition")
    int robin_dirichlet_id_;
    int robin_neumann_id_;

    // Get the condition for the dirichlet and neumann condition associated with the Robin-condition
    std::vector<DRT::Condition*> cutterele_cond_robin_dirichlet_;
    std::vector<DRT::Condition*> cutterele_cond_robin_neumann_;


  };  // End LevelSetCouplingNavierSlip



  /// set levelset field by given vector
  void WriteAccess_GeometricQuantities(Teuchos::RCP<Epetra_Vector>& scalaraf,
      Teuchos::RCP<Epetra_MultiVector>& smoothed_gradphiaf,
      Teuchos::RCP<Epetra_Vector>& curvatureaf);


  /// set material pointer for coupling slave side
  void GetInterfaceSlaveMaterial(DRT::Element* actele, Teuchos::RCP<MAT::Material>& mat);

  template <CORE::FE::CellType DISTYPE, class M1, class M2>
  void EvaluateCurvature(double& icurvature,  ///< curvature to be computed
      int eid,                                ///< element ID
      M1& funct,  ///< local shape function for Gauss Point (from fluid element)
      M2& derxy   ///< local derivatives of shape function for Gauss Point (from fluid element)
  )
  {
    // number space dimensions for element
    const size_t nsd = CORE::FE::dim<DISTYPE>;

    // number of nodes of element
    const size_t nen = CORE::FE::num_nodes<DISTYPE>;

    // smoothed normal at cutter element nodes, the Gaussian point lies in
    CORE::LINALG::SerialDenseMatrix esmoothedgradphi(nsd, nen);
    CORE::LINALG::SerialDenseMatrix esmoothedcurvature(nen, 1);

    CORE::LINALG::Matrix<nsd, nen> esmoothedgradphi_T(esmoothedgradphi, View);
    CORE::LINALG::Matrix<nen, 1> esmoothedcurvature_T(esmoothedcurvature, View);

    return;
  }

  template <CORE::FE::CellType DISTYPE, class M1, class M2>
  void GetPhiAtGP(double& phi_gp,  ///< phi at gausspoint
      int eid,                     ///< element ID
      M1& funct,                   ///< local shape function for Gauss Point (from fluid element)
      M2& derxy  ///< local derivatives of shape function for Gauss Point (from fluid element)
  )
  {
    // number space dimensions for element
    //    const size_t nsd = CORE::FE::dim<DISTYPE>;

    // number of nodes of element
    const size_t nen = CORE::FE::num_nodes<DISTYPE>;

    // smoothed normal at cutter element nodes, the Gaussian point lies in
    CORE::LINALG::SerialDenseMatrix ephinp(nen, 1);
    CORE::LINALG::Matrix<nen, 1> ephinp_T(ephinp, View);

    phi_gp = funct.Dot(ephinp_T);
  }

  template <CORE::FE::CellType DISTYPE, class M1, class M2, class M3, class M4>
  double InterpolateCurvature(const M1& funct, const M2& derxy, const M3& esmoothedgradphi_T,
      const M4& esmoothedcurvature_T, const int numnode)
  {
    // number space dimensions for element
    const size_t nsd = CORE::FE::dim<DISTYPE>;

    double icurvature = 0.0;
    return icurvature;
  }


  template <CORE::FE::CellType DISTYPE, class T1, class M1, class M2>
  void EvalProjectionMatrix(
      T1& itraction_jump_matrix,  ///< Laplace-Beltrami matrix for surface tension calculations
      int eid,                    ///< element ID
      M1& funct,                  ///< local shape function for Gauss Point (from fluid element)
      M2& normal                  ///< surface normal of cut element
  )
  {
    // number space dimensions for element
    const size_t nsd = CORE::FE::dim<DISTYPE>;

    CORE::LINALG::Matrix<nsd, nsd> p_matrix(false);
    CORE::LINALG::Matrix<nsd, nsd> p_smoothed_matrix(false);

    //  -----------------------------------------------------------
    //
    // HERE DEPENDING ON HOW THE PROJECTION SHOULD BE CALCULATED
    // different ways of the projection matrix should be tested.
    //
    // Might, only want to use matrix_mixed_smoothed option as this is the most promising!
    //   Check for stabilized Laplace-Beltrami option!
    //    Burman, Erik and Hansbo, Peter and Larson, Mats G
    //    A stabilized cut finite element method for partial differential equations on surfaces:
    //    The Laplace--Beltrami operator Computer Methods in Applied Mechanics and Engineering
    //       2015
    //
    //  -----------------------------------------------------------

    SetupConcatenatedProjectionMatrix(itraction_jump_matrix, p_matrix, p_smoothed_matrix);

    return;
  }


};  // namespace XFEM

BACI_NAMESPACE_CLOSE

#endif  // XFEM_COUPLING_LEVELSET_H
