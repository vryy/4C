/*----------------------------------------------------------------------*/
/*!
\file fluid_ele_calc_xfem_coupling_impl.H

\brief Classes for interface coupling in the XFEM

<pre>
Maintainer: Raffaela Kruse /Benedikt Schott
            kruse@lnm.mw.tum.de
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

namespace DRT
{
  namespace ELEMENTS
  {
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
          LINALG::Matrix<nsd_,1> & x_gp_lin,       ///< global coordinates of gaussian point w.r.t linearized interface
          LINALG::Matrix<nsd_,1> & x_side,         ///< projected gaussian point on side
          LINALG::Matrix<nsd_-1,1> & xi_side       ///< local coordinates of projected gaussian point w.r.t side
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
        LINALG::Matrix<nsd_,slave_nen_> slave_veln_;          ///< coupling slave element's velocity at previous step
        LINALG::Matrix<slave_nen_,1>    slave_presn_;         ///< coupling slave element's pressure at previous step
        LINALG::Matrix<nsd_,nsd_>       slave_vderxyn_;       ///< coupling slave element's velocity derivatives at previous step

        std::string disp_statename_; ///< name of current displacement state (for access from discretization)
        std::string vel_statename_;  ///< name of current velocity state (for access from discretization)
        std::string veln_statename_;   ///< name of previous velocity state (for access from discretization)

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
            bool isViscAdjointSymmetric = true
        ) :
          SlaveElementRepresentation<distype,slave_distype,slave_numdof>(),
          C_umum_(C_umum.A(), true),
          rhC_um_(rhC_um.A(), true),
          applicationType_(SlaveElementInterface<distype>::XFluidWDBC),
          adj_visc_scale_(isViscAdjointSymmetric ? 1.0 : -1.0)
        {
        }

        //! ctor for one-sided (xfluid weak dirichlet) problems (interface defined by mesh)
        NitscheCoupling(
          Epetra_SerialDenseMatrix &  slave_xyze, ///< global node coordinates of slave element
          Epetra_SerialDenseMatrix &  C_umum,          ///< C_umum coupling matrix
          Epetra_SerialDenseMatrix &  rhC_um,          ///< C_um coupling rhs
          bool isViscAdjointSymmetric = true
        ) :
          SlaveElementRepresentation<distype,slave_distype,slave_numdof>(slave_xyze),
          C_umum_(C_umum.A(), true),
          rhC_um_(rhC_um.A(), true),
          applicationType_(SlaveElementInterface<distype>::XFluidWDBC),
          adj_visc_scale_(isViscAdjointSymmetric ? 1.0 : -1.0)

        {
        }

        //! ctor for two-sided problems
        NitscheCoupling(
          Epetra_SerialDenseMatrix &  slave_xyze,      ///< global node coordinates of slave element
          Epetra_SerialDenseMatrix &  C_umum,          ///< C_umum coupling matrix
          Epetra_SerialDenseMatrix &  C_usum,          ///< C_usum coupling matrix
          Epetra_SerialDenseMatrix &  C_umus,          ///< C_umus coupling matrix
          Epetra_SerialDenseMatrix &  C_usus,          ///< C_usus coupling matrix
          Epetra_SerialDenseMatrix &  rhC_um,          ///< C_um coupling rhs
          Epetra_SerialDenseMatrix &  rhC_us,          ///< C_us coupling rhs
          typename SlaveElementInterface<distype>::ApplicationType applicationType = SlaveElementInterface<distype>::XFluidFluid, ///< application type
          bool isViscAdjointSymmetric = true           ///< symmetric or skew-symmetric formulation of Nitsche's adjoint viscous term
        ) :
          SlaveElementRepresentation<distype,slave_distype,slave_numdof>(slave_xyze),
          C_umum_(C_umum.A(), true),
          C_usum_(C_usum.A(), true),
          C_umus_(C_umus.A(), true),
          C_usus_(C_usus.A(), true),
          rhC_um_(rhC_um.A(), true),
          rhC_us_(rhC_us.A(), true),
          applicationType_(applicationType),
          adj_visc_scale_(isViscAdjointSymmetric ? 1.0 : -1.0)
        {
        }

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
          INPAR::XFEM::XFF_ConvStabScaling                      xff_conv_stab = INPAR::XFEM::XFF_ConvStabScaling_none         ///< type of convective stabilization in XFF-problems
        );

        //! build coupling matrices and assemble terms for Nitsche's (NIT) method
        void NIT_evaluateCoupling(
          const LINALG::Matrix<nsd_,1> &    normal,                 ///< outward pointing normal (defined by the coupling partner, that determines the interface traction)
          const double &                    timefacfac,             ///< theta*dt*fac
          const double &                    visceff_m,              ///< viscosity in coupling master fluid
          const double &                    visceff_s,              ///< viscosity in coupling slave fluid
          const double &                    kappa_m,                ///< mortaring weight for coupling master
          const double &                    kappa_s,                ///< mortaring weight for coupling slave
          const double &                    density_m,              ///< fluid density (master)
          const double &                    NIT_full_stab_fac,      ///< full Nitsche's penalty term scaling (viscous+convective part)
          const LINALG::Matrix<nen_,1> &    funct_m,                ///< coupling master shape functions
          const LINALG::Matrix<nsd_,nen_> & derxy_m,                ///< spatial derivatives of coupling master shape functions
          const LINALG::Matrix<nsd_,nsd_> & vderxy_m,               ///< coupling master spatial velocity derivatives
          const double &                    pres_m,                 ///< coupling master pressure
          const LINALG::Matrix<nsd_,1> &    velint_m,               ///< coupling master interface velocity
          INPAR::XFEM::XFF_ConvStabScaling  xff_conv_stab = INPAR::XFEM::XFF_ConvStabScaling_none ///< type of convective stabilization in XFF-problems
        );

        //! add rhs contributions from old time step in Nitsche's (NIT) method
        void NIT_evaluateCouplingOldState(
          const LINALG::Matrix<nsd_,1> &    normal,                  ///< outward pointing normal (defined by the coupling partner, that determines the interface traction)
          const double &                    timefacfacn,             ///< dt*(1-theta)*fac
          const double &                    visceff_m,               ///< viscosity in coupling master fluid
          const double &                    visceff_s,               ///< viscosity in coupling slave fluid
          const double &                    kappa_m,                 ///< mortaring weight for coupling master
          const double &                    kappa_s,                 ///< mortaring weight for coupling slave
          const LINALG::Matrix<nen_,1> &    funct_m,                 ///< coupling master shape functions
          const LINALG::Matrix<nsd_,nen_> & derxy_m,                 ///< spatial derivatives of coupling master shape functions
          const LINALG::Matrix<nsd_,nsd_> & vderxyn_m,               ///< coupling master spatial velocity derivatives
          const double &                    presn_m,                 ///< coupling master pressure
          const LINALG::Matrix<nsd_,1> &    velintn_m                ///< coupling master interface velocity
        );

      private:

        //! evaluate pressure-consistency term for Nitsche's method (scaled with master side's weighting factor)
        void NIT_p_Consistency_MasterTerms(
          const double &                            pres_m,                       ///< master pressure
          const LINALG::Matrix<nen_,1> &            funct_m_timefacfac_km,        ///< funct * timefacfac *kappa_m
          const LINALG::Matrix<slave_nen_,1> &      funct_s_timefacfac_km,        ///< funct_s * timefacfac *kappa_m
          const LINALG::Matrix<nsd_,1> &            normal,                       ///< normal vector
          const LINALG::Matrix<slave_nen_,nen_> &   funct_s_m_dyad_timefacfac_km, ///< (funct_s^T * funct) * timefacfac *kappa_m
          const LINALG::Matrix<nen_,nen_> &         funct_m_m_dyad_timefacfac_km  ///< (funct^T * funct) * timefacfac *kappa_m
        );

        //! evaluate pressure-consistency term for Nitsche's method (scaled with slave side's weighting factor)
        void NIT_p_Consistency_SlaveTerms(
          const double &                                pres_s,                       ///< slave pressure
          const LINALG::Matrix<nen_,1> &                funct_m_timefacfac_ks,        ///< funct * timefacfac *kappa_m
          const LINALG::Matrix<slave_nen_,1> &          funct_s_timefacfac_ks,        ///< funct_s * timefacfac *kappa_m
          const LINALG::Matrix<slave_nen_,nen_> &       funct_s_m_dyad_timefacfac_ks, ///< (funct_s^T * funct) * timefacfac *kappa_m
          const LINALG::Matrix<slave_nen_,slave_nen_> & funct_s_s_dyad_timefacfac_ks, ///< (funct^T * funct) * timefacfac *kappa_m
          const LINALG::Matrix<nsd_,1> &                normal                        ///< normal vector
        );

        //! evaluate pressure-adjoint-consistency term for Nitsche's method (scaled with master side's weighting factor)
        void NIT_p_AdjointConsistency_MasterTerms(
          const double &                          velint_normal_m,              ///< velocity in normal direction
          const double &                          velint_normal_s,              ///< interface velocity in normal direction
          const LINALG::Matrix<nen_,1> &          funct_m_timefacfac_km,        ///< funct * timefacfac *kappa_m
          const LINALG::Matrix<slave_nen_,nen_> & funct_s_m_dyad_timefacfac_km, ///< (funct^T * funct_s) * timefacfac *kappa_m
          const LINALG::Matrix<nen_,nen_> &       funct_m_m_dyad_timefacfac_km, ///< (funct^T * funct) * timefacfac *kappa_m
          const LINALG::Matrix<nsd_,1> &          normal                        ///< normal vector
        );

        //! evaluate pressure-adjoint-consistency term for Nitsche's method (scaled with slave side's weighting factor)
        void NIT_p_AdjointConsistency_SlaveTerms(
          const double &                                velint_normal_m,              ///< velocity in normal direction
          const double &                                velint_normal_s,              ///< interface velocity in normal direction
          const LINALG::Matrix<slave_nen_,1> &          funct_s_timefacfac_ks,        ///< funct_s * timefacfac * kappas
          const LINALG::Matrix<slave_nen_,nen_> &       funct_s_m_dyad_timefacfac_ks, ///< (funct_s^T * funct) * timefacfac *kappas
          const LINALG::Matrix<slave_nen_,slave_nen_> & funct_s_s_dyad_timefacfac_ks, ///< (funct_s^T * funct_s) * timefacfac *kappas
          const LINALG::Matrix<nsd_,1> &                normal                        ///< normal vector
        );

        //! evaluate viscous-consistency term for Nitsche's method (scaled with master side's weighting factor)
        void NIT_visc_Consistency_MasterTerms(
          const LINALG::Matrix<nsd_,nen_> &        derxy_m,                      ///< master deriv
          const LINALG::Matrix<nsd_,nsd_> &        vderxy_m,                     ///< master velocity gradient
          const LINALG::Matrix<nen_,1> &           funct_m_viscm_timefacfac_km,  ///< funct_m*mu_m*timefacfac
          const LINALG::Matrix<slave_nen_,1> &     funct_s_viscm_timefacfac_km,  ///< funct_s*mu_m*timefacfac
          const LINALG::Matrix<nsd_,1> &           normal                        ///< normal vector
        );

        //! evaluate viscous-consistency term for Nitsche's method (scaled with slave side's weighting factor)
        void NIT_visc_Consistency_SlaveTerms(
          const LINALG::Matrix<nsd_,slave_nen_> &  derxy_s,                      ///< slave shape function derivatives
          const LINALG::Matrix<nsd_,nsd_> &        vderxy_s,                     ///< slave velocity gradient
          const LINALG::Matrix<nen_,1> &           funct_m_viscs_timefacfac_ks,  ///< funct_m*mu_m*timefacfac
          const LINALG::Matrix<slave_nen_,1> &     funct_s_viscs_timefacfac_ks,  ///< funct_s*mu_m*timefacfac
          const LINALG::Matrix<nsd_,1> &           normal                        ///< normal vector
        );

        //! evaluate viscous-adjoint-consistency term for Nitsche's method (scaled with master side's weighting factor)
        void NIT_visc_AdjointConsistency_MasterTerms(
          const LINALG::Matrix<nsd_,1>&            velint_m,                          ///< velocity
          const LINALG::Matrix<nsd_,1>&            velint_s,                          ///< interface velocity
          const LINALG::Matrix<nsd_,nen_> &        derxy_m_viscm_timefacfac_km,       ///< master shape function derivatives * timefacfac * 2 * mu_m * kappa_m
          const LINALG::Matrix<nen_,1> &           funct_m,                           ///< embedded element funct *mu*timefacfac
          const LINALG::Matrix<slave_nen_,1> &     funct_s,                           ///< embedded element funct *mu*timefacfac
          const LINALG::Matrix<nsd_,1> &           normal                             ///< normal vector
        );

        //! evaluate viscous-adjoint-consistency term for Nitsche's method (scaled with slave side's weighting factor)
        void NIT_visc_AdjointConsistency_SlaveTerms(
          const LINALG::Matrix<nsd_,1>&            velint_m,                          ///< velocity
          const LINALG::Matrix<nsd_,1>&            velint_s,                          ///< interface velocity
          const LINALG::Matrix<nen_,1> &           funct_m,                           ///< embedded element funct *mu*timefacfac
          const LINALG::Matrix<slave_nen_,1> &     funct_s,                           ///< embedded element funct *mu*timefacfac
          const LINALG::Matrix<nsd_,slave_nen_> &  derxy_s_viscs_timefacfac_ks,       ///< master shape function derivatives * timefacfac * 2 * mu_m * kappa_m
          const LINALG::Matrix<nsd_,1> &           normal
        );

        //! evaluate rhs of pressure-consistency term for Nitsche's method (scaled with master side's weighting factor)
        void NIT_p_Consistency_MasterRHS(
          const double &                            pres_m,                       ///< master pressure
          const LINALG::Matrix<nen_,1> &            funct_m_timefacfac_km,        ///< funct * timefacfac *kappa_m
          const LINALG::Matrix<slave_nen_,1> &      funct_s_timefacfac_km,        ///< funct_s * timefacfac *kappa_m
          const LINALG::Matrix<nsd_,1> &            normal                        ///< normal vector
        );

        //! evaluate rhs of pressure-consistency term for Nitsche's method (scaled with slave side's weighting factor)
        void NIT_p_Consistency_SlaveRHS(
          const double &                                pres_s,                       ///< slave pressure
          const LINALG::Matrix<nen_,1> &                funct_m_timefacfac_ks,        ///< funct * timefacfac *kappa_m
          const LINALG::Matrix<slave_nen_,1> &          funct_s_timefacfac_ks,        ///< funct_s * timefacfac *kappa_m
          const LINALG::Matrix<nsd_,1> &                normal                        ///< normal vector
        );

        //! evaluate rhs of pressure-adjoint-consistency term for Nitsche's method (scaled with master side's weighting factor)
        void NIT_p_AdjointConsistency_MasterRHS(
          const double &                          velint_normal_m,              ///< velocity in normal direction
          const double &                          velint_normal_s,              ///< interface velocity in normal direction
          const LINALG::Matrix<nen_,1> &          funct_m_timefacfac_km,        ///< funct * timefacfac *kappa_m
          const LINALG::Matrix<nsd_,1> &          normal                        ///< normal vector
        );

        //! evaluate rhs of pressure-adjoint-consistency term for Nitsche's method (scaled with slave side's weighting factor)
        void NIT_p_AdjointConsistency_SlaveRHS(
          const double &                                velint_normal_m,              ///< velocity in normal direction
          const double &                                velint_normal_s,              ///< interface velocity in normal direction
          const LINALG::Matrix<slave_nen_,1> &          funct_s_timefacfac_ks         ///< funct_s * timefacfac * kappas
        );

        //! evaluate rhs of viscous-consistency term for Nitsche's method (scaled with master side's weighting factor)
        void NIT_visc_Consistency_MasterRHS(
          const LINALG::Matrix<nsd_,nsd_> &        vderxy_m,                     ///< master velocity gradient
          const LINALG::Matrix<nen_,1> &           funct_m_viscm_timefacfac_km,  ///< funct_m*mu_m*timefacfac
          const LINALG::Matrix<slave_nen_,1> &     funct_s_viscm_timefacfac_km,  ///< funct_s*mu_m*timefacfac
          const LINALG::Matrix<nsd_,1> &           normal                        ///< normal vector
        );

        //! evaluate rhs of viscous-consistency term for Nitsche's method (scaled with slave side's weighting factor)
        void NIT_visc_Consistency_SlaveRHS(
          const LINALG::Matrix<nsd_,nsd_> &        vderxy_s,                     ///< slave velocity gradient
          const LINALG::Matrix<nen_,1> &           funct_m_viscs_timefacfac_ks,  ///< funct_m*mu_m*timefacfac
          const LINALG::Matrix<slave_nen_,1> &     funct_s_viscs_timefacfac_ks,  ///< funct_s*mu_m*timefacfac
          const LINALG::Matrix<nsd_,1> &           normal                        ///< normal vector
        );

        //! evaluate rhs of viscous-adjoint-consistency term for Nitsche's method (scaled with master side's weighting factor)
        void NIT_visc_AdjointConsistency_MasterRHS(
          const LINALG::Matrix<nsd_,1>&            velint_m,                          ///< velocity
          const LINALG::Matrix<nsd_,1>&            velint_s,                          ///< interface velocity
          const LINALG::Matrix<nsd_,nen_> &        derxy_m_viscm_timefacfac_km,       ///< master shape function derivatives * timefacfac * 2 * mu_m * kappa_m
          const LINALG::Matrix<nsd_,1> &           normal                             ///< normal vector
        );

        //! evaluate rhs of viscous-adjoint-consistency term for Nitsche's method (scaled with slave side's weighting factor)
        void NIT_visc_AdjointConsistency_SlaveRHS(
          const LINALG::Matrix<nsd_,1>&            velint_m,                          ///< velocity
          const LINALG::Matrix<nsd_,1>&            velint_s,                          ///< interface velocity
          const LINALG::Matrix<nsd_,slave_nen_> &  derxy_s_viscs_timefacfac_ks,       ///< master shape function derivatives * timefacfac * 2 * mu_m * kappa_m
          const LINALG::Matrix<nsd_,1> &           normal
        );

        //!  evaluate Nitsche's penalty term
        void NIT_Stab_Penalty(
          const LINALG::Matrix<nsd_,1>&                 velint_m,                     ///< velocity at integration point
          const LINALG::Matrix<nsd_,1>&                 velint_s,                     ///< interface velocity at integration point
          const LINALG::Matrix<nen_,1>&                 funct_m_timefacfac,           ///< funct * timefacfac
          const LINALG::Matrix<nen_,nen_>&              funct_m_m_dyad_timefacfac,    ///< (funct^T * funct) * timefacfac
          const LINALG::Matrix<slave_nen_,1>&           funct_s_timefacfac,           ///< funct_s^T * timefacfac
          const LINALG::Matrix<slave_nen_,nen_>&        funct_s_m_dyad_timefacfac,    ///< (funct_s^T * funct) * timefacfac
          const LINALG::Matrix<slave_nen_,slave_nen_>&  funct_s_s_dyad_timefacfac,    ///< (funct_s^T * funct_s) * timefacfac
          const double &                                stabfac                       ///< stabilization factor
        );

        //!  evaluate only the master side's part of Nitsche's penalty term
        void NIT_Stab_Penalty_MasterTerms(
          const LINALG::Matrix<nsd_,1>&                 velint_m,                     ///< velocity at integration point
          const LINALG::Matrix<nsd_,1>&                 velint_s,                     ///< interface velocity at integration point
          const LINALG::Matrix<nen_,1>&                 funct_m_timefacfac,           ///< funct * timefacfac
          const LINALG::Matrix<nen_,nen_>&              funct_m_m_dyad_timefacfac,    ///< (funct^T * funct) * timefacfac
          const double &                                stabfac
        );

        //! add stabilizing terms due to cross-interface convective mass transport (fluid-fluid only)
        void NIT_Stab_Inflow_AveragedTerm(
          const LINALG::Matrix<nsd_,1>&                 velint_m,                     ///< velocity at integration point
          const LINALG::Matrix<nsd_,1>&                 velint_s,                     ///< interface velocity at integration point
          const LINALG::Matrix<nen_,1>&                 funct_m_timefacfac,           ///< funct * timefacfac
          const LINALG::Matrix<nen_,nen_>&              funct_m_m_dyad_timefacfac,    ///< (funct^T * funct) * timefacfac
          const LINALG::Matrix<slave_nen_,1>&           funct_s_timefacfac,           ///< funct_s^T * timefacfac
          const LINALG::Matrix<slave_nen_,nen_>&        funct_s_m_dyad_timefacfac,    ///< (funct_s^T * funct) * timefacfac
          const LINALG::Matrix<slave_nen_,slave_nen_> & funct_s_s_dyad_timefacfac,    ///< (funct_s^T * funct_s) * timefacfac
          const LINALG::Matrix<nsd_,1>&                 normal,                       ///< normal vector n^m
          const double &                                density                       ///< fluid density
        );

      private:

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

        LINALG::Matrix<master_numdof_*nen_,master_numdof_*nen_>         C_umum_;          ///< coupling matrix C_umum
        LINALG::Matrix<slave_numdof*slave_nen_,master_numdof_*nen_>     C_usum_;          ///< coupling matrix C_usum
        LINALG::Matrix<master_numdof_*nen_,slave_numdof*slave_nen_>     C_umus_;          ///< coupling matrix C_umus
        LINALG::Matrix<slave_numdof*slave_nen_,slave_numdof*slave_nen_> C_usus_;          ///< coupling matrix C_usus
        LINALG::Matrix<master_numdof_*nen_,1>                           rhC_um_;          ///< coupling rhs rhC_um
        LINALG::Matrix<slave_numdof*slave_nen_,1>                       rhC_us_;          ///< coupling rhs rhC_us

        typename SlaveElementInterface<distype>::ApplicationType        applicationType_; ///< application context

        /// scaling of Nitsche's adjoint viscous term
        const double adj_visc_scale_;

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

        //! ctor for xfluid weak dirichlet problem
        HybridLMCoupling(
          Epetra_SerialDenseMatrix &  slave_xyze, ///< global node coordinates of slave element
          bool isViscAdjointSymmetric = true      ///< symmetric or skew-symmetric formulation of Nitsche's adjoint viscous term
        ) :
          SlaveElementRepresentation<distype,slave_distype,slave_numdof>(slave_xyze),
          applicationType_(SlaveElementInterface<distype>::XFluidWDBC),
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
          typename SlaveElementInterface<distype>::ApplicationType applicationType = SlaveElementInterface<distype>::XFluidFluid,
          bool isViscAdjointSymmetric = false          ///< symmetric or skew-symmetric formulation of Nitsche's adjoint viscous term
        ) :
          SlaveElementRepresentation<distype,slave_distype,slave_numdof>(slave_xyz),
          C_usum_(C_usum.A(), true),
          C_umus_(C_umus.A(), true),
          rhC_us_(rhC_us.A(), true),
          G_sus_(G_s_us.A(), true),
          G_uss_(G_us_s.A(), true),
          applicationType_(applicationType),
          adj_visc_scale_(isViscAdjointSymmetric ? 1.0 : -1.0)
        {
        }

        //! evaluate interface coupling matrices for mixed/hybrid Cauchy stress-based (MHCS) coupling
        void MHCS_buildCouplingMatrices(
          const LINALG::Matrix<nsd_,1> &                                normal,     ///< normal vector
          const double &                                                fac,        ///< integration factor
          const LINALG::Matrix<nen_,1> &                                funct,      ///< shape function
          LINALG::BlockMatrix<LINALG::Matrix<nen_,1>,numstressdof_,1> & rhs_s       ///< block rhs vector \f$ rhs_{\sigma} \f$
        );

        //! evaluate interface matrices for mixed/hybrid viscous stress-based (MHVS) coupling
        void MHVS_buildCouplingMatrices(
          const LINALG::Matrix<nsd_,1> &                                normal,     ///< normal vector
          const double &                                                fac,        ///< integration factor
          const LINALG::Matrix<nen_,1> &                                funct,      ///< background element shape functions
          LINALG::BlockMatrix<LINALG::Matrix<nen_,1>,numstressdof_,1> & rhs_s,      ///< block rhs vector \f$ rhs_{\sigma}\f$
          const double &                                                press,      ///< background element pressure
          LINALG::Matrix<nen_,1> &                                      rhs_pmus    ///< part of block rhs vector \f$rhs_p\f$ including interface velocity terms
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

        typename SlaveElementInterface<distype>::ApplicationType applicationType_;              ///< application context

        /// scaling of Nitsche's adjoint viscous term
        const double adj_visc_scale_;
      };
    }
  }
}

#endif /* FLUID_ELE_CALC_XFEM_COUPLING_IMPL_H_ */
