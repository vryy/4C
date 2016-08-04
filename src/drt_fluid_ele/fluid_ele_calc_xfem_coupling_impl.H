/*----------------------------------------------------------------------*/
/*!
\file fluid_ele_calc_xfem_coupling_impl.H

\brief Classes for interface coupling in the XFEM

\level 2

<pre>
\maintainer Benedikt Schott
            schott@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15240
</pre>
*/
/*----------------------------------------------------------------------*/

#ifndef FLUID_ELE_CALC_XFEM_COUPLING_IMPL_H_
#define FLUID_ELE_CALC_XFEM_COUPLING_IMPL_H_

#include "../drt_inpar/inpar_xfem.H"

#include "fluid_ele_calc_xfem_coupling_impl.H"

//   qnuP - option SHOULD BE ON!
//     projects the given velocity into normal direction in case there
//     is a smoothed projection matrix given.
#define PROJECT_VEL_FOR_PRESSURE_ADJOINT

//  For comparison to Urquizas paper with his slip length implementation
//   Only working for Navier-Slip (i.e. itraction_jump_ = 0)!
//   Can be modified to work for this case as well.
//#define ENFORCE_URQUIZA_GNBC

namespace DRT
{
  namespace ELEMENTS
  {
    class FluidEleParameterXFEM;
    namespace XFLUID
    {
      //! class for concrete coupling slave element
      //! this can be an arbitrary 2D/3D element and can be associated with a structure (monolithic XFSI),
      //! fluid (XFF, XFFSI, partitioned XFSI, XWDBC) or a xfluid-element with another active dofset (two-phase flow)
      template<DRT::Element::DiscretizationType distype, DRT::Element::DiscretizationType slave_distype, unsigned int slave_numdof>
      class SlaveElementRepresentation : virtual public SlaveElementInterface<distype>
      {

      public:

        /// number of nodes per master (xfem-fluid) element
        static const unsigned nen_ = DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;
        /// number of spatial dimensions of the master element (xfem-fluid)
        static const unsigned nsd_ = DRT::UTILS::DisTypeToDim<distype>::dim;
        /// number of spatial dimensions of the slave side
        static const unsigned slave_nsd_ = DRT::UTILS::DisTypeToDim<slave_distype>::dim;
        /// number of slave element's nodes
        static const unsigned slave_nen_ = DRT::UTILS::DisTypeToNumNodePerEle<slave_distype>::numNodePerElement;

        //! ctor
        SlaveElementRepresentation(Epetra_SerialDenseMatrix &  slave_xyze) :
          slave_xyze_(slave_xyze.A(), true){
          SlaveElementInterface<distype>::DefineStateNames(slave_distype,disp_statename_,vel_statename_,veln_statename_); };

        //! add coupling slave element's displacements and set current slave element node coordinates
        void AddSlaveEleDisp(
          const DRT::Discretization &  slavedis,     ///< coupling slave discretization
          const std::vector<int>&      lm            ///< local map
          );

        //! set slave element's interface velocity & pressure for current time step
        void SetSlaveState(
          const DRT::Discretization &  slavedis,     ///< coupling slave discretization
          const std::vector<int>&      lm            ///< local map
          );

        //! set slave element's interface velocity & pressure for previous time step
        void SetSlaveStaten(
          const DRT::Discretization &  slavedis,     ///< coupling slave discretization
          const std::vector<int>&      lm            ///< local map
          );

        //! @name Accessors
        //@{

        //! extract interface velocity at current time step
        void GetInterfaceVelnp(
          LINALG::Matrix<nsd_,1>& ivelint  ///< interface velocity at coupling slave side
        ) const;

        //! get interface pressure at current time step
        void GetInterfacePresnp(
          double & ipres ///< interface pressure at coupling slave side
        ) const;

        //! get interface velocity gradient at current time step
        void GetInterfaceVelGradnp(
          LINALG::Matrix<nsd_,nsd_>& velgradint  ///< interface velocity gradients at coupling slave side
        ) const;

        //! extract interface velocity at previous time step
        void GetInterfaceVeln(
          LINALG::Matrix<nsd_,1>& ivelintn  ///< interface velocity at coupling slave side
        ) const;

        //! get interface pressure at previous time step
        void GetInterfacePresn(
          double & ipresn ///< interface pressure at coupling slave side
        ) const;

        //! get interface velocity gradient at previous time step
        void GetInterfaceVelGradn(
          LINALG::Matrix<nsd_,nsd_>& velgradintn  ///< interface velocity gradients at coupling slave side
        ) const;

        //! get slave elements nodal shape functions
        void GetSlaveFunct(
          LINALG::Matrix<slave_nen_,1> & slave_funct ///< coupling slave shape functions
          ) const;

        //! set state for interface velocity jump
        void SetInterfaceJumpStatenp(
            const DRT::Discretization &  cutterdis,      ///< cutter discretization
            const std::string            state,          ///< state
            const std::vector<int>&      lm              ///< local map
        );

        //! set state for interface velocity jump for previous time step
        void SetInterfaceJumpStaten(
          const DRT::Discretization &  cutterdis,      ///< cutter discretization
          const std::string            state,          ///< state
          const std::vector<int>&      lm              ///< local map
        );

        //! get interface velocity jump at Gaussian point
        void GetInterfaceJumpVelnp(
            LINALG::Matrix<nsd_,1> & ivelint_jump ///< cutter element interface velocity jump or prescribed DBC at Gaussian point
        ) const;

        //! get interface velocity jump for previous time step at Gaussian point
        void GetInterfaceJumpVeln(
            LINALG::Matrix<nsd_,1> & ivelintn_jump ///< cutter element interface velocity jump or prescribed DBC at Gaussian point
        ) const;

        //@}

        //!  evaluate shape function, derivatives and transformation w.r.t coupling slave element at gaussian point
        void Evaluate( LINALG::Matrix<nsd_,1> & xside );

        //! compute coupling slave element's characteristic element length (for 3D-elements)
        void ComputeCharElementLength( double & h_k );

        //! compute interface force
        void ComputeInterfaceForce(
          Epetra_SerialDenseVector &   iforce,     ///< interface force vector
          LINALG::Matrix<nsd_,1> &     traction,   ///< traction vector at gaussian point
          const double &               fac         ///< integration factor
        );

        //! project gaussian point from linearized interface in normal direction onto corresponding side
        void ProjectOnSide(
          LINALG::Matrix<nsd_,1> & x_gp_lin,     ///< global coordinates of gaussian point w.r.t linearized interface
          LINALG::Matrix<nsd_,1> & x_side,       ///< projected gaussian point on side
          LINALG::Matrix<nsd_,1> & xi_side       ///< local coordinates of projected gaussian point w.r.t side
        );

        //! evaluate shape functions and derivatives at element center
        //! @return the integration factor
        double EvalShapeFuncAndDerivsAtEleCenter();

      protected:

        //! default constructor
        SlaveElementRepresentation()
        {SlaveElementInterface<distype>::DefineStateNames(slave_distype,disp_statename_,vel_statename_,veln_statename_);};

        //! evaluate shape functions and derivatives at integration point
        double EvalShapeFuncAndDerivsAtIntPoint(
          const double* gpcoord,  //< actual integration point (coords)
          double gpweight         //< actual integration point (weight)
          );

        //! @name accessors for derived classes
        //@{

        //! get spatial derivatives of slave elements nodal shape functions
        void GetSlaveFunctDeriv(LINALG::Matrix<nsd_,slave_nen_>& slave_derxy) const;

        //@}

      private:

        LINALG::Matrix<nsd_,slave_nen_> slave_xyze_;         ///< coupling slave element's node coordinates
        LINALG::Matrix<slave_nen_,1>    slave_funct_;        ///< coupling slave element's shape functions
        LINALG::Matrix<nsd_,slave_nen_> slave_derxy_;        ///< coupling slave element's local shape function derivatives
        LINALG::Matrix<nsd_,slave_nen_> slave_deriv_;        ///< coupling slave element's global shape function derivatives
        LINALG::Matrix<nsd_,slave_nen_> slave_vel_;          ///< coupling slave element's velocity at current step
        LINALG::Matrix<slave_nen_,1>    slave_pres_;         ///< coupling slave element's pressure at current step
        LINALG::Matrix<nsd_,nsd_>       slave_vderxy_;       ///< coupling slave element's velocity derivatives at current step
        LINALG::Matrix<nsd_,slave_nen_> slave_disp_;         ///< coupling slave element's displacements at current step
        LINALG::Matrix<nsd_,slave_nen_> slave_veln_;         ///< coupling slave element's velocity at previous step
        LINALG::Matrix<slave_nen_,1>    slave_presn_;        ///< coupling slave element's pressure at previous step
        LINALG::Matrix<nsd_,nsd_>       slave_vderxyn_;      ///< coupling slave element's velocity derivatives at previous step

        std::string disp_statename_; ///< name of current displacement state (for access from discretization)
        std::string vel_statename_;  ///< name of current velocity state (for access from discretization)
        std::string veln_statename_;   ///< name of previous velocity state (for access from discretization)

        //TODO: shift this vector to the same class as ProjectOnSide which is based on the cutter discretization
        LINALG::Matrix<nsd_, slave_nen_> interface_velnp_jump_; ///< cutter element's prescribed velocity jump height or prescribed DBC values
        LINALG::Matrix<nsd_, slave_nen_> interface_veln_jump_; ///< cutter element's prescribed velocity jump height or prescribed DBC values

        LINALG::Matrix<slave_nen_,1> proj_funct_;      ///< shape functions for project on side
        LINALG::Matrix<2,slave_nen_> proj_deriv_;      ///< derivatives dr, ds for project on side
        LINALG::Matrix<3,slave_nen_> proj_deriv2_;     ///< 2nd derivatives drdr, dsds, drds for project on side

        LINALG::Matrix<3,1> proj_x_;             ///< global coordinates
        LINALG::Matrix<3,2> proj_derxy_;         ///< global xyz derivatives
        LINALG::Matrix<3,3> proj_derxy2_;        ///< global xyz 2nd derivatives

        LINALG::Matrix<3,1> proj_residuum_;      ///<  residuum of the newton iteration
        LINALG::Matrix<3,3> proj_sysmat_;        ///<  matrix for the newton system
        LINALG::Matrix<3,1> proj_incr_;          ///<  increment of the newton system

        LINALG::Matrix<3,1> proj_sol_;           ///< sol carries xi_1, xi_2, d (distance)

        // get vector products
        LINALG::Matrix<3,1> proj_dx_drdr_times_dx_ds_;
        LINALG::Matrix<3,1> proj_dx_dr_times_dx_drds_;
        LINALG::Matrix<3,1> proj_dx_drds_times_dx_ds_;
        LINALG::Matrix<3,1> proj_dx_dr_times_dx_dsds_;
        LINALG::Matrix<3,1> proj_dx_dr_times_dx_ds_;



      };

      /*!
       * specialized interface class for XFluid weak Dirichlet problems with a interface given by a level-set field
       * (we then don't couple with a concrete slave element!)
       */
      //!
      template<DRT::Element::DiscretizationType distype, unsigned int slave_numdof>
      class SlaveElementRepresentation<distype,DRT::Element::dis_none,slave_numdof>
      {
      public:

        /// number of nodes per master (xfem-fluid) element
        static const unsigned nen_ = DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;
        /// number of spatial dimensions of the master element (xfem-fluid)
        static const unsigned nsd_ = DRT::UTILS::DisTypeToDim<distype>::dim;
        /// number of spatial dimensions of the slave side
        static const unsigned slave_nsd_ = DRT::UTILS::DisTypeToDim<distype>::dim;
        /// just for compatibility...
        static const unsigned slave_nen_ = nen_;

        //! get slave elements nodal shape functions - if the interface is given as a level-set
        //! field, the request is unfulfilled
        void GetSlaveFunct(
          LINALG::Matrix<slave_nen_,1> & slave_funct ///< coupling slave shape functions
          ) const
        {
          dserror("There is no concrete slave element available.");
        };

      protected:

        //! default ctor
        SlaveElementRepresentation(){};

        //! ctor
        SlaveElementRepresentation(Epetra_SerialDenseMatrix &  slave_xyze){};

        //! get nodal shape function derivatives
        void GetSlaveFunctDeriv(LINALG::Matrix<nsd_,slave_nen_>& slave_derxy) const
        {
          dserror("There is no concrete slave element available.");
        };

      };

      //! concrete evaluation class for interface coupling using Nitsche's method
      template<DRT::Element::DiscretizationType distype, DRT::Element::DiscretizationType slave_distype, unsigned int slave_numdof>
      class NitscheCoupling : public NitscheInterface<distype>, public SlaveElementRepresentation<distype,slave_distype,slave_numdof>
      {

      public:

        /// number of nodes per master (xfem-fluid) element
        static const unsigned nen_ = SlaveElementRepresentation<distype,slave_distype,slave_numdof>::nen_;
        /// number of spatial dimensions
        static const unsigned nsd_ = SlaveElementInterface<distype>::nsd_;
        /// number of nodal dof for master element (coupling master is always a fluid element!)
        static const unsigned master_numdof_ =  nsd_ + 1;
        /// number of slave element's nodes
        static const unsigned slave_nen_ = SlaveElementRepresentation<distype,slave_distype,slave_numdof>::slave_nen_;

        //! ctor for one-sided (xfluid weak dirichlet) problems (interface defined by level-set field)
        NitscheCoupling(
            Epetra_SerialDenseMatrix &  C_umum,          ///< C_umum coupling matrix
            Epetra_SerialDenseMatrix &  rhC_um,          ///< C_um coupling rhs
            const DRT::ELEMENTS::FluidEleParameterXFEM& fldparaxfem ///< specific XFEM based fluid parameters
        );

        //! ctor for one-sided (xfluid weak dirichlet) problems (interface defined by mesh)
        NitscheCoupling(
          Epetra_SerialDenseMatrix &  slave_xyze,                  ///< global node coordinates of slave element
          Epetra_SerialDenseMatrix &  C_umum,                      ///< C_umum coupling matrix
          Epetra_SerialDenseMatrix &  rhC_um,                      ///< C_um coupling rhs
          const DRT::ELEMENTS::FluidEleParameterXFEM& fldparaxfem  ///< specific XFEM based fluid parameters
        );

        //! ctor for two-sided problems
        NitscheCoupling(
          Epetra_SerialDenseMatrix &  slave_xyze,      ///< global node coordinates of slave element
          Epetra_SerialDenseMatrix &  C_umum,          ///< C_umum coupling matrix
          Epetra_SerialDenseMatrix &  C_usum,          ///< C_usum coupling matrix
          Epetra_SerialDenseMatrix &  C_umus,          ///< C_umus coupling matrix
          Epetra_SerialDenseMatrix &  C_usus,          ///< C_usus coupling matrix
          Epetra_SerialDenseMatrix &  rhC_um,          ///< C_um coupling rhs
          Epetra_SerialDenseMatrix &  rhC_us,          ///< C_us coupling rhs
          const DRT::ELEMENTS::FluidEleParameterXFEM& fldparaxfem ///< specific XFEM based fluid parameters
        );

        //! add contributions from convective stabilization
        //! this method is applied in a non-Nitsche context (e.g. MHVS) by
        //! employing shape functions and velocities from another slave element coupling object
        void ApplyConvStabTerms(
          const Teuchos::RCP<SlaveElementInterface<distype> > & slave_ele,              ///< associated slave element coupling object
          const LINALG::Matrix<nen_,1> &                        funct_m,                ///< master shape functions
          const LINALG::Matrix<nsd_,1> &                        velint_m,               ///< vector of slave shape functions
          const LINALG::Matrix<nsd_,1> &                        normal,                 ///< normal vector n^b
          const double &                                        density_m,              ///< fluid density (master)
          const double &                                        NIT_stab_fac_conv,      ///< full Nitsche's penalty term scaling (viscous+convective part)
          const double &                                        timefacfac,             ///< theta*dt
          const LINALG::Matrix<nsd_,1> &                        ivelint_jump,           ///< prescribed interface velocity, Dirichlet values or jump height for coupled problems
          const INPAR::XFEM::EleCouplingCondType &              cond_type              ///< condition type
        );

        //! build coupling matrices and assemble terms for Nitsche's (NIT) method
        void NIT_evaluateCoupling(
          const LINALG::Matrix<nsd_,1> &    normal,                 ///< outward pointing normal (defined by the coupling partner, that determines the interface traction)
          const double &                    timefacfac,             ///< theta*dt*fac
          const double &                    pres_timefacfac,        ///< pressure time scaling * fac
          const double &                    visceff_m,              ///< viscosity in coupling master fluid
          const double &                    visceff_s,              ///< viscosity in coupling slave fluid
          const double &                    kappa_m,                ///< mortaring weight for coupling master
          const double &                    kappa_s,                ///< mortaring weight for coupling slave
          const double &                    density_m,              ///< fluid density (master)
          const LINALG::Matrix<nen_,1> &    funct_m,                ///< coupling master shape functions
          const LINALG::Matrix<nsd_,nen_> & derxy_m,                ///< spatial derivatives of coupling master shape functions
          const LINALG::Matrix<nsd_,nsd_> & vderxy_m,               ///< coupling master spatial velocity derivatives
          const double &                    pres_m,                 ///< coupling master pressure
          const LINALG::Matrix<nsd_,1> &    velint_m,               ///< coupling master interface velocity
          const LINALG::Matrix<nsd_,1> &    ivelint_jump,           ///< prescribed interface velocity, Dirichlet values or jump height for coupled problems
          const LINALG::Matrix<nsd_,1> &    itraction_jump,         ///< prescribed interface traction, jump height for coupled problems
          const LINALG::Matrix<nsd_,nsd_>&  itraction_jump_matrix,  ///< prescribed projection matrix for laplace-beltrami problems
          const bool                        is_traction_jump,        ///< is it a normal traction jump or calculated throw laplace-beltrami
          std::map<INPAR::XFEM::CoupTerm, std::pair<bool,double> >& configmap ///< Interface Terms configuration map
        );

        //! add rhs contributions from old time step in Nitsche's (NIT) method
        void NIT_evaluateCouplingOldState(
          const LINALG::Matrix<nsd_,1> &           normal,                  ///< outward pointing normal (defined by the coupling partner, that determines the interface traction)
          const double &                           timefacfacn,             ///< dt*(1-theta)*fac
          bool                                     isImplPressure,          ///< flag for implicit pressure treatment
          const double &                           visceff_m,               ///< viscosity in coupling master fluid
          const double &                           visceff_s,               ///< viscosity in coupling slave fluid
          const double &                           kappa_m,                 ///< mortaring weight for coupling master
          const double &                           kappa_s,                 ///< mortaring weight for coupling slave
          const double &                           density_m,               ///< fluid density (master) USED IN XFF
          const LINALG::Matrix<nen_,1> &           funct_m,                 ///< coupling master shape functions
          const LINALG::Matrix<nsd_,nen_> &        derxy_m,                 ///< spatial derivatives of coupling master shape functions
          const LINALG::Matrix<nsd_,nsd_> &        vderxyn_m,               ///< coupling master spatial velocity derivatives
          const double &                           presn_m,                 ///< coupling master pressure
          const LINALG::Matrix<nsd_,1> &           velintn_m,               ///< coupling master interface velocity
          const LINALG::Matrix<nsd_,1> &           ivelintn_jump,           ///< prescribed interface velocity, Dirichlet values or jump height for coupled problems
          const LINALG::Matrix<nsd_,1> &           itractionn_jump,         ///< prescribed interface traction, jump height for coupled problems
          std::map<INPAR::XFEM::CoupTerm, std::pair<bool,double> >& configmap ///< Interface Terms configuration map
        );

      private:

        //! evaluate traction-consistency term for Nitsche's method
        void NIT_Traction_Consistency_Term(
          const LINALG::Matrix<nen_,1> &            funct_m_timefacfac_ks,        ///< funct * timefacfac *kappa_s
          const LINALG::Matrix<slave_nen_,1> &      funct_s_timefacfac_km,        ///< funct_s * timefacfac *kappa_m
          const LINALG::Matrix<nsd_,1> &            itraction_jump                ///< prescribed interface traction, jump height for coupled problems
        );

        //! evaluate traction-consistency term for Nitsche's method (for integration by parts aproach)
       void NIT_Projected_Traction_Consistency_Term(
         const LINALG::Matrix<nsd_,nen_> &        derxy_m_timefacfac_ks,       ///< master shape function derivatives * timefacfac * kappa_s
         const LINALG::Matrix<nsd_,slave_nen_> &  derxy_s_timefacfac_km,       ///< slave shape function derivatives * timefacfac * kappa_m
         const LINALG::Matrix<nsd_,nsd_> &        itraction_jump_matrix        ///< prescribed projection matrix
       );

        //! evaluate pressure-consistency term for Nitsche's method (scaled with master side's weighting factor)
        void NIT_p_Consistency_MasterTerms(
          const double &                            pres_m,                     ///< master pressure
          const LINALG::Matrix<nen_,1> &            funct_m,                    ///< funct
          const LINALG::Matrix<nsd_,1> &            normal_timefacfac,          ///< normal vector * timefacfac
          const std::pair<bool,double>& m_row,                                  ///< scaling for master row
          const std::pair<bool,double>& s_row,                                  ///< scaling for slave row
          const std::pair<bool,double>& m_col,                                  ///< scaling for master col
          bool only_rhs = false                                                 ///< evaluat only rhs
        );

        //! evaluate pressure-consistency term for Nitsche's method (scaled with slave side's weighting factor)
        void NIT_p_Consistency_SlaveTerms(
          const double &                            pres_s,                        ///< slave pressure
          const LINALG::Matrix<nen_,1> &            funct_m,                       ///< funct
          const LINALG::Matrix<nsd_,1> &            normal_timefacfac_ks,          ///< normal vector * timefacfac
          const std::pair<bool,double>& m_row,                                     ///< scaling for master row
          const std::pair<bool,double>& s_row,                                     ///< scaling for slave row
          const std::pair<bool,double>& s_col,                                     ///< scaling for slave col
          bool only_rhs = false                                                    ///< evaluat only rhs
        );

        //! evaluate pressure-adjoint-consistency term for Nitsche's method (scaled with master side's weighting factor)
        void NIT_p_AdjointConsistency_MasterTerms(
          const LINALG::Matrix<nen_,1> &            funct_m,                          ///< funct
          const LINALG::Matrix<nsd_,1> &            normal_timefacfac,                ///< normal vector * timefacfac
          const double &                            velint_diff_normal_timefacfac,    ///< (velint_m - velint_s) * normal * timefacfac
          const std::pair<bool,double>& m_row,                                        ///< scaling for master row
          const std::pair<bool,double>& m_col,                                        ///< scaling for master col
          const std::pair<bool,double>& s_col,                                        ///< scaling for slave row
          bool only_rhs = false                                                       ///< evaluat only rhs
        );

        //! evaluate pressure-adjoint-consistency term for Nitsche's method (scaled with slave side's weighting factor)
        void NIT_p_AdjointConsistency_SlaveTerms(
            const LINALG::Matrix<nsd_,1> &          normal_timefacfac,              ///< normal vector * timefacfac
            const double &                          velint_diff_normal_timefacfac,  ///< (velint_m - velint_s) * normal * timefacfac
            const std::pair<bool,double>& s_row,                                    ///< scaling for slave row
            const std::pair<bool,double>& m_col,                                    ///< scaling for master col
            const std::pair<bool,double>& s_col,                                    ///< scaling for slave col
            bool only_rhs = false                                                   ///< evaluat only rhs
        );

        //! evaluate viscous-consistency term for Nitsche's method (scaled with master side's weighting factor)
        void NIT_visc_Consistency_MasterTerms(
          const LINALG::Matrix<nsd_,nen_> &         derxy_m,                      ///< master deriv
          const LINALG::Matrix<nen_,1> &            funct_m,                      ///< funct_m
          const std::pair<bool,double>& m_row,                                    ///< scaling for master row
          const std::pair<bool,double>& s_row,                                    ///< scaling for slave row
          const std::pair<bool,double>& m_col,                                    ///< scaling for master col
          bool only_rhs = false                                                   ///< evaluat only rhs
        );

        //! evaluate viscous-consistency term for Nitsche's method (scaled with master side's weighting factor)
        void NIT_visc_Consistency_MasterTerms_Projected(
          const LINALG::Matrix<nsd_,nen_> &         derxy_m,                      ///< master deriv
          const LINALG::Matrix<nen_,1> &            funct_m,                      ///< funct_m
          const LINALG::Matrix<nsd_,nsd_>&          proj_matrix                   ///< pürojection matrix
        );

        //! evaluate viscous-consistency term for Nitsche's method (scaled with slave side's weighting factor)
        void NIT_visc_Consistency_SlaveTerms(
          const LINALG::Matrix<nsd_,slave_nen_> &   derxy_s,                      ///< slave shape function derivatives
          const LINALG::Matrix<nen_,1> &            funct_m,                      ///< funct_m
          const std::pair<bool,double>& m_row,                                    ///< scaling for master row
          const std::pair<bool,double>& s_row,                                    ///< scaling for slave row
          const std::pair<bool,double>& s_col,                                    ///< scaling for slave col
          bool only_rhs = false                                                   ///< evaluat only rhs
        );

        //! evaluate viscous-adjoint-consistency term for Nitsche's method (scaled with master side's weighting factor)
        void NIT_visc_AdjointConsistency_MasterTerms(
            const LINALG::Matrix<nen_,1>&                 funct_m,                      ///< funct * timefacfac
            const LINALG::Matrix<nsd_,nen_> &             derxy_m,                      ///< spatial derivatives of coupling master shape functions
            const LINALG::Matrix<nsd_,1>&                 normal,                       ///< normal-vector
            const double &                                viscm_fac,                    ///< scaling factor
            const std::pair<bool,double>& m_row,                                        ///< scaling for master row
            const std::pair<bool,double>& m_col,                                        ///< scaling for master col
            const std::pair<bool,double>& s_col,                                        ///< scaling for slave col
            bool only_rhs = false                                                       ///< evaluat only rhs
        );

        //! evaluate viscous-adjoint-consistency term for Nitsche's method (scaled with master side's weighting factor)
        void NIT_visc_AdjointConsistency_MasterTerms_Projected(
            const LINALG::Matrix<nsd_,nen_> &         derxy_m_viscm_timefacfac_km,  ///< master shape function derivatives * timefacfac * 2 * mu_m * kappa_m
            const LINALG::Matrix<nen_,1> &            funct_m,                      ///< funct_m
            const LINALG::Matrix<nsd_,1> &            normal,                        ///< normal vector
            const std::pair<bool,double>& m_row, ///< scaling for master row
            const std::pair<bool,double>& m_col, ///< scaling for master col
            const std::pair<bool,double>& s_col ///< scaling for slave col
        );

        //! evaluate traction-traction term
        void NIT_visc_Neumann_AdjointConsistency_MasterTerms_Projected(
          const LINALG::Matrix<nsd_,nen_> &        derxy_m_viscm_timefacfac_km,       ///< master shape function derivatives * timefacfac * 2 * mu_m * kappa_m
          const LINALG::Matrix<nsd_,nen_> &        derxy_m,                           ///< master deriv
          const LINALG::Matrix<nsd_,nsd_> &        vderxy_m,                          ///< coupling master spatial velocity derivatives
          const LINALG::Matrix<nen_,1> &           funct_m,                           ///< embedded element funct *mu*timefacfac
          const LINALG::Matrix<nsd_,1> &           normal,                            ///< normal vector
          const std::pair<bool,double>& m_row, ///< scaling for master row
          const std::pair<bool,double>& mstr_col ///< scaling for master col
          );

        //! evaluate viscous-adjoint-consistency term for Nitsche's method (scaled with slave side's weighting factor)
        void NIT_visc_AdjointConsistency_SlaveTerms(
          const LINALG::Matrix<nen_,1> &            funct_m,                      ///< funct_m
          const LINALG::Matrix<nsd_,slave_nen_> &   derxy_s_viscs_timefacfac_ks,  ///< master shape function derivatives * timefacfac * 2 * mu_m * kappa_m
          const LINALG::Matrix<nsd_,1> &            normal,                       ///< normal vector
          const std::pair<bool,double>& s_row,                                    ///< scaling for slave row
          const std::pair<bool,double>& m_col,                                    ///< scaling for master col
          const std::pair<bool,double>& s_col,                                    ///< scaling for slave col
          bool only_rhs = false                                                   ///< evaluat only rhs
        );

        //!  evaluate Nitsche's penalty term
        void NIT_Stab_Penalty(
            const LINALG::Matrix<nen_,1>&          funct_m,         ///< funct
            const double &                         timefacfac,      ///< time integration factor
            const std::pair<bool,double>& m_row,                    ///< scaling for master row
            const std::pair<bool,double>& s_row,                    ///< scaling for slave row
            const std::pair<bool,double>& m_col,                    ///< scaling for master col
            const std::pair<bool,double>& s_col,                    ///< scaling for slave col
            bool only_rhs = false                                   ///< evaluat only rhs
        );

        //!  evaluate Nitsche's penalty term
        void NIT_Stab_Penalty_Projected(
            const LINALG::Matrix<nen_,1>&          funct_m,                      ///< funct
            const LINALG::Matrix<nsd_,nsd_>&       projection_matrix,            ///< projection_matrix
            const LINALG::Matrix<nsd_,1>&          velint_diff_proj_matrix,      ///< velocity difference projected
            const double &                         timefacfac,                    ///< time integration factor
            const std::pair<bool,double>& m_row, ///< scaling for master row
            const std::pair<bool,double>& s_row, ///< scaling for slave row
            const std::pair<bool,double>& m_col, ///< scaling for master col
            const std::pair<bool,double>& s_col ///< scaling for slave col
        );

        //! evaluate viscous penalty for the Neumann terms of the Robin-condition
        void NIT_visc_Neumann_Penalty_MasterTerms_Projected(
          const LINALG::Matrix<nsd_,nen_> &     derxy_m,                      ///< master deriv
          const LINALG::Matrix<nen_,1> &        funct_m,                      ///< funct_m
          const double &                        km_viscm,                     ///< traction jump
          const std::pair<bool,double>& m_row, ///< scaling for penalty master row
          const std::pair<bool,double>& s_row, ///< scaling for penalty slave row
          const std::pair<bool,double>& mstr_col ///< scaling for penalty master stress col
        );

        //! add stabilizing terms due to cross-interface convective mass transport (fluid-fluid only)
        void NIT_Stab_Inflow_AveragedTerm(
          const LINALG::Matrix<nen_,1>&            funct_m,         ///< funct
          const LINALG::Matrix<nsd_,1>&            velint_m,        ///< master velocity
          const LINALG::Matrix<nsd_,1>&            normal,          ///< normal vector n^m
          const double &                           density,         ///< fluid density
          const double &                           timefacfac,      ///< timefac * fac
          bool only_rhs = false                                     ///< evaluat only rhs
        );

        //! Do Nitsche consistency and adjoint consistency terms with projection
        void Do_NIT_visc_Adjoint_and_Neumann_MasterTerms_Projected(
            const LINALG::Matrix<nen_,1>&                 funct_m,                      ///< funct * timefacfac
            const LINALG::Matrix<nsd_,nen_> &             derxy_m,                      ///< spatial derivatives of coupling master shape functions
            const LINALG::Matrix<nsd_,nsd_> &             vderxy_m,                     ///< coupling master spatial velocity derivatives
            const LINALG::Matrix<nsd_,nsd_>&              projection_matrix,            ///< projection_matrix
            const LINALG::Matrix<nsd_,1>&                 velint_diff_proj_matrix,      ///< velocity difference projected
            const LINALG::Matrix<nsd_,1>&                 normal,                       ///< normal-vector
            const double &                                km_viscm_fac,                  ///< scaling factor
            const LINALG::Matrix<nsd_,1>&                 itraction_jump,
            const std::pair<bool,double>& m_row, ///< scaling for master row
            const std::pair<bool,double>& m_col, ///< scaling for master col
            const std::pair<bool,double>& s_col, ///< scaling for slave col
            const std::pair<bool,double>& mstr_col, ///< scaling for master stress col
            //these last three parameters will be removed later as the terms are just like projected consistency terms
            const std::pair<bool,double>& m_pen_row, ///< scaling for penalty master row
            const std::pair<bool,double>& s_pen_row, ///< scaling for penalty slave row
            const std::pair<bool,double>& mstr_pen_col ///< scaling for penalty master stress col
        );


      private:

        void NIT_Create_Standard_Projection_Matrices(
            const LINALG::Matrix<nsd_,1>&          normal   ///< normal vector
        );

        // Get the stabilization parameters for the specific problem.
        void GetStabilizationParameters(
            const double &                                NIT_full_stab_fac,            ///< full Nitsche stab fac
            const double &                                NIT_visc_stab_fac,            ///< viscous Nitsche stab fac
            double &                                      stabnit,                      ///< stabilization factor NIT_Penalty
            double &                                      stabepsnit,                   ///< stabilization factor NIT_Penalty Neumann
            double &                                      stabadj,                      ///< stabilization factor Adjoint
            double &                                      stabepsadj,                   ///< stabilization factor Adjoint Neumann
            const bool &                                  sliplength_not_zero           ///< bool for now, add a inpar for more options?
        );

        //! @name useful constants for DOF-index numbering
        //@{
        static const unsigned Velx = 0;
        static const unsigned Vely = 1;
        static const unsigned Velz = 2;
        static const unsigned Pres = 3;
        //@}

        /// get global master row/col-index of element coupling matrix for a given node index and dof-index
        static unsigned mIndex(unsigned inod, unsigned idof)
        {
//          if (idof >= master_numdof_)
//            dserror("Coupling master element has only %d dof!", master_numdof_);
//          if (inod >= nen_)
//            dserror("Coupling master element has only %d nodes!", nen_);
          return inod*master_numdof_ + idof;
        }

        /// get global slave row/col-index of element coupling matrix for a given node index and dof-index
        static unsigned sIndex(unsigned inod, unsigned idof)
        {
//          if (idof >= slave_numdof)
//            dserror("Coupling slave element has only %d dof!", slave_numdof);
//          if (inod >= slave_nen_)
//            dserror("Coupling slave element has only %d nodes!", slave_nen_);
          return inod*slave_numdof + idof;
        }

        //! @name get global master row/col-index of element coupling matrix for a given node index
        //@{
        static unsigned mVelx(unsigned inod){ return inod*master_numdof_ + Velx; }
        static unsigned mVely(unsigned inod){ return inod*master_numdof_ + Vely; }
        static unsigned mVelz(unsigned inod){ return inod*master_numdof_ + Velz; }
        static unsigned mPres(unsigned inod){ return inod*master_numdof_ + Pres; }
        //@}

        /// @name get global slave row/col-index of element coupling matrix for a given node index
        //@{
        static unsigned sVelx(unsigned inod){ return inod*slave_numdof + Velx; }
        static unsigned sVely(unsigned inod){ return inod*slave_numdof + Vely; }
        static unsigned sVelz(unsigned inod){ return inod*slave_numdof + Velz; }
        static unsigned sPres(unsigned inod){ return inod*slave_numdof + Pres; }
        //@}

        /// specific XFEM based fluid parameters
        const DRT::ELEMENTS::FluidEleParameterXFEM& fldparaxfem_;

        LINALG::Matrix<master_numdof_*nen_,master_numdof_*nen_>         C_umum_;          ///< coupling matrix C_umum
        LINALG::Matrix<slave_numdof*slave_nen_,master_numdof_*nen_>     C_usum_;          ///< coupling matrix C_usum
        LINALG::Matrix<master_numdof_*nen_,slave_numdof*slave_nen_>     C_umus_;          ///< coupling matrix C_umus
        LINALG::Matrix<slave_numdof*slave_nen_,slave_numdof*slave_nen_> C_usus_;          ///< coupling matrix C_usus
        LINALG::Matrix<master_numdof_*nen_,1>                           rhC_um_;          ///< coupling rhs rhC_um
        LINALG::Matrix<slave_numdof*slave_nen_,1>                       rhC_us_;          ///< coupling rhs rhC_us

        /// scaling of Nitsche's adjoint viscous term
        const double adj_visc_scale_;

        const bool eval_coupling_; ///< do we have to evaluate coupling terms?

        LINALG::Matrix<nsd_,1> velint_s_; ///< velocity at integration point on slave side

        LINALG::Matrix<slave_nen_,1> funct_s_;

        LINALG::Matrix<slave_nen_,1> funct_s_timefacfac_km_;

        LINALG::Matrix<nen_,1> funct_m_timefacfac_ks_;

        LINALG::Matrix<nen_,nen_> funct_m_m_dyad_;

        LINALG::Matrix<slave_nen_,nen_> funct_s_m_dyad_;

        LINALG::Matrix<slave_nen_,slave_nen_> funct_s_s_dyad_;

        LINALG::Matrix<nsd_,nen_> derxy_m_viscm_timefacfac_; // dN^(nen)/dx_i * mu_m * kappa_m

        LINALG::Matrix<nsd_,1> normal_pres_timefacfac_;

        LINALG::Matrix<nsd_,1> normal_pres_timefacfac_km_;

        LINALG::Matrix<nsd_,1> normal_pres_timefacfac_ks_;

        LINALG::Matrix<nsd_,1> half_normal_;

        LINALG::Matrix<nsd_,1> half_normal_viscm_timefacfac_km_;

        LINALG::Matrix<nsd_,1> half_normal_viscs_timefacfac_ks_;

        LINALG::Matrix<nen_,1> half_normal_deriv_m_viscm_timefacfac_km_;

        LINALG::Matrix<slave_nen_,1> half_normal_deriv_s_viscs_timefacfac_ks_;

        LINALG::Matrix<nen_,1> normal_deriv_m_viscm_km_;  // dN^(nen)/dx_i * n * mu_m * kappa_m

        LINALG::Matrix<slave_nen_,1> normal_deriv_s_viscs_ks_;

        LINALG::Matrix<nsd_,1> vderxy_m_normal_;

        LINALG::Matrix<nsd_,1> vderxy_s_normal_;

        LINALG::Matrix<nsd_,1> vderxy_m_normal_transposed_viscm_timefacfac_km_;

        LINALG::Matrix<nsd_,1> vderxy_s_normal_transposed_viscs_timefacfac_ks_;

        //Velocity difference between slave and master
        LINALG::Matrix<nsd_,1> velint_diff_;
        //Velocity difference between slave and master for normal projection
        LINALG::Matrix<nsd_,1> velint_diff_normal_;
        //Velocity difference between slave and master for tangential projection
        LINALG::Matrix<nsd_,1> velint_diff_tangential_;

        //Only exists in NIT_Stab_Penalty_MasterTerms and
        //               NIT_Stab_Penalty
        //   i.e. probably should define, normal and tangential components
        LINALG::Matrix<nsd_,1> velint_diff_timefacfac_stabfac_;

        //Needed members for GNBC
        //  Projection matrices:
        LINALG::Matrix<nsd_,nsd_> proj_tangential_;
        LINALG::Matrix<nsd_,nsd_> proj_normal_;
        LINALG::Matrix<nsd_,nsd_> proj_matrix_;

        //Projected
        LINALG::Matrix<nsd_,1> vderxy_x_normal_transposed_viscx_timefacfac_kx_pmatrix_;

        //Help variables:
        LINALG::Matrix<nsd_,nsd_> velint_proj_norm_diff_dyad_normal_, velint_proj_norm_diff_dyad_normal_symm_;
        LINALG::Matrix<nsd_,1>    vderxy_m_normal_tang_, vderxy_m_normal_transposed_;

        //  Stabilization parameters:
        double dyn_visc_;

        //  projected velocity components:
        LINALG::Matrix<nsd_,1> velint_diff_proj_normal_;
        LINALG::Matrix<nsd_,1> velint_diff_proj_tangential_;
        LINALG::Matrix<nsd_,1> velint_diff_proj_matrix_;

        //  projected traction jump components:
        LINALG::Matrix<nsd_,1> itraction_jump_proj_matrix_;

        //ConsistencyNeumann help variable
        LINALG::Matrix<nsd_,nen_> proj_matrix_derxy_m_; // (beta) * 2.0 * mu_m * timefacefac * km * p_1(IX,j)

        //AdjointNeumann help variables
        LINALG::Matrix<nen_,1>    normal_deriv_m_;     // 2.0 * half_normal(k) * derxy_m(k,ix)
        LINALG::Matrix<nen_,nen_> derxy_m_P_derxy_m_;  // 2.0 * derxy_m(j,IC) P^t_{jk} * derxy_m(k,IR) * mu_m * timefacfac * km

        double velint_diff_normal_pres_timefacfac_;
        double velint_diff_pres_timefacfac_;
      };

      /// concrete class for interface coupling using mixed/hybrid stress-based Lagrange multipliers method
      template<DRT::Element::DiscretizationType distype, DRT::Element::DiscretizationType slave_distype, unsigned int slave_numdof>
      class HybridLMCoupling : public HybridLMInterface<distype>, public SlaveElementRepresentation<distype,slave_distype,slave_numdof>
      {

      public:

        /// number of nodes per master element
        static const unsigned nen_ = SlaveElementRepresentation<distype,slave_distype,slave_numdof>::nen_;
        /// number of spatial dimensions
        static const unsigned nsd_ = SlaveElementInterface<distype>::nsd_;
        /// number of nodal dof for master element (always a fluid element)
        static const unsigned master_numdof_ =  nsd_ + 1;
        /// number of slave element's nodes
        static const unsigned slave_nen_ = SlaveElementRepresentation<distype,slave_distype,slave_numdof>::slave_nen_;
        /// number of independent stress-dof
        static const unsigned numstressdof_ = HybridLMInterface<distype>::numstressdof_;

        //! ctor for one-sided (xfluid weak dirichlet) problems (interface defined by level-set field)
        HybridLMCoupling(
            bool isViscAdjointSymmetric = true
        ) :
          SlaveElementRepresentation<distype,slave_distype,slave_numdof>(),
          adj_visc_scale_(isViscAdjointSymmetric ? 1.0 : -1.0)
        {
        }

        //! ctor for xfluid weak dirichlet problem
        HybridLMCoupling(
          Epetra_SerialDenseMatrix &  slave_xyze, ///< global node coordinates of slave element
          bool isViscAdjointSymmetric = true      ///< symmetric or skew-symmetric formulation of Nitsche's adjoint viscous term
        ) :
          SlaveElementRepresentation<distype,slave_distype,slave_numdof>(slave_xyze),
          adj_visc_scale_(isViscAdjointSymmetric ? 1.0 : -1.0)
        {
        }

        //! ctor for fluid-fluid
        HybridLMCoupling(
          Epetra_SerialDenseMatrix &  slave_xyz,       ///< global node coordinates of slave element
          Epetra_SerialDenseMatrix &  C_usum,          ///< C_usum coupling matrix
          Epetra_SerialDenseMatrix &  C_umus,          ///< C_umus coupling matrix
          Epetra_SerialDenseMatrix &  rhC_us,          ///< C_us coupling rhs
          Epetra_SerialDenseMatrix &  G_s_us,          ///< \f$G_{u^s \sigma}\f$ coupling matrix
          Epetra_SerialDenseMatrix &  G_us_s,          ///< \f$G_{\sigma u^s}\f$ coupling matrix
          bool isViscAdjointSymmetric = false          ///< symmetric or skew-symmetric formulation of Nitsche's adjoint viscous term
        ) :
          SlaveElementRepresentation<distype,slave_distype,slave_numdof>(slave_xyz),
          C_usum_(C_usum.A(), true),
          C_umus_(C_umus.A(), true),
          rhC_us_(rhC_us.A(), true),
          G_sus_(G_s_us.A(), true),
          G_uss_(G_us_s.A(), true),
          adj_visc_scale_(isViscAdjointSymmetric ? 1.0 : -1.0)
        {
        }

        //! evaluate interface coupling matrices for mixed/hybrid Cauchy stress-based (MHCS) coupling
        void MHCS_buildCouplingMatrices(
          const LINALG::Matrix<nsd_,1> &                                normal,     ///< normal vector
          const double &                                                fac,        ///< integration factor
          const LINALG::Matrix<nen_,1> &                                funct,      ///< shape function
          LINALG::BlockMatrix<LINALG::Matrix<nen_,1>,numstressdof_,1> & rhs_s,      ///< block rhs vector \f$ rhs_{\sigma} \f$
          const LINALG::Matrix<nsd_,1> &                                ivelint_jump,  ///< prescribed interface velocity or interface jump height
          const LINALG::Matrix<nsd_,1> &                                itraction_jump ///< prescribed interface traction or interface jump height
        );

        //! evaluate interface matrices for mixed/hybrid viscous stress-based (MHVS) coupling
        void MHVS_buildCouplingMatrices(
          const LINALG::Matrix<nsd_,1> &                                normal,     ///< normal vector
          const double &                                                fac,        ///< integration factor
          const LINALG::Matrix<nen_,1> &                                funct,      ///< background element shape functions
          LINALG::BlockMatrix<LINALG::Matrix<nen_,1>,numstressdof_,1> & rhs_s,      ///< block rhs vector \f$ rhs_{\sigma}\f$
          const double &                                                press,      ///< background element pressure
          LINALG::Matrix<nen_,1> &                                      rhs_pmus,   ///< part of block rhs vector \f$rhs_p\f$ including interface velocity terms
          const LINALG::Matrix<nsd_,1> &                                ivelint_jump,  ///< prescribed interface velocity or interface jump height
          const LINALG::Matrix<nsd_,1> &                                itraction_jump ///< prescribed interface traction or interface jump height
        );

        //! apply the standard consistency traction interface jump term
        void MH_Traction_Consistency_Term(
            const LINALG::Matrix<slave_nen_,1> &      funct_s_timefacfac_km,        ///< funct_s * timefacfac *kappa_m
            const LINALG::Matrix<nsd_,1> &            itraction_jump                ///< prescribed interface traction, jump height for coupled problems
        );

        //! build the final coupling matrices for mixed/hybrid Cauchy or viscous stress-based coupling (MHCS or MHVS)
        void HybridLM_buildFinalCouplingMatrices(
          LINALG::BlockMatrix<LINALG::Matrix<nen_,nen_>,numstressdof_,numstressdof_> &  BinvK_ss,    ///< block inverse \f$ K^{-1}_{\sigma\sigma} \f$
          LINALG::BlockMatrix<LINALG::Matrix<nen_,nen_>,master_numdof_,numstressdof_> & BKumsInvKss, ///< block matrix \f$ K_{u\sigma} \cdot K^{-1}_{\sigma\sigma} \f$
          LINALG::BlockMatrix<LINALG::Matrix<nen_,nen_>,numstressdof_,master_numdof_> & BK_sum,      ///< block matrix \f$ K_{\sigma u} \f$
          LINALG::BlockMatrix<LINALG::Matrix<nen_,1>, numstressdof_,1> &                rhs_s        ///< block rhs vector \f$ rhs_{\sigma}\f$
        );

      protected:

        //! @name useful constants for DOF-index numbering
        //@{
        static const unsigned Velx = 0;
        static const unsigned Vely = 1;
        static const unsigned Velz = 2;
        static const unsigned Pres = 3;
        //@}

        /// get stress dof-index
        static unsigned stressIndex(unsigned xi, unsigned xj)
        {
          if ( xi > 2 || xj > 2)
            dserror("Invalid index combination (%d,%d) for stress tensor!", xi, xj);

          return (xi * xj) ? xi + xj + 1 : xi + xj;
        }

        /// get global master row/col-index of element coupling matrix for a given node index and dof-index
        static unsigned mIndex(unsigned inod, unsigned idof)
        {
//          if (idof >= master_numdof_)
//            dserror("Coupling master element has only %d dof!", master_numdof_);
//          if (inod >= nen_)
//            dserror("Coupling master element has only %d nodes!", nen_);
          return inod*master_numdof_ + idof;
        }

        /// get global slave row/col-index of element coupling matrix for a given node index and dof-index
        static unsigned sIndex(unsigned inod,unsigned idof)
        {
//          if (idof >= slave_numdof)
//            dserror("Coupling slave element has only %d dof!", slave_numdof);
//          if (inod >= slave_nen_)
//            dserror("Coupling slave element has only %d nodes!", slave_nen_);
          return inod*slave_numdof + idof;
        }

        LINALG::Matrix<slave_numdof*slave_nen_,master_numdof_*nen_>              C_usum_;       ///< coupling matrix C_usum
        LINALG::Matrix<master_numdof_*nen_,slave_numdof*slave_nen_>              C_umus_;       ///< coupling matrix C_umus
        LINALG::Matrix<slave_numdof*slave_nen_,1>                                rhC_us_;       ///< coupling rhs rhC_us
        LINALG::Matrix<numstressdof_*nen_,slave_numdof*slave_nen_>               G_sus_;        ///< G_sus coupling matrix
        LINALG::Matrix<slave_numdof*slave_nen_,numstressdof_*nen_>               G_uss_;        ///< G_uss coupling matrix
        LINALG::BlockMatrix<LINALG::Matrix<nen_,slave_nen_>,numstressdof_,nsd_>  BG_sus_;       ///< block G_sus coupling matrix
        LINALG::BlockMatrix<LINALG::Matrix<slave_nen_,nen_>,nsd_,numstressdof_>  BG_uss_;       ///< block G_uss coupling matrix
        LINALG::BlockMatrix<LINALG::Matrix<nen_,slave_nen_>,1,nsd_>              BG_pmus_;      ///< block G_pmus coupling matrix
        LINALG::BlockMatrix<LINALG::Matrix<slave_nen_,nen_>,nsd_,1>              BG_uspm_;      ///< block G_uspm coupling matrix

        /// scaling of Nitsche's adjoint viscous term
        const double adj_visc_scale_;
      };
    }
  }
}

#endif /* FLUID_ELE_CALC_XFEM_COUPLING_IMPL_H_ */
