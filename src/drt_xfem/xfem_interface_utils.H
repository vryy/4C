/*----------------------------------------------------------------------*/
/*!
\brief Basic routings to evaluate the terms for Nitsche Interface

\level 2

\maintainer  Christoph Ager
             ager@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289-15249

*/
/*----------------------------------------------------------------------*/

#ifndef XFEM_INTERFACE_UTILS_H_
#define XFEM_INTERFACE_UTILS_H_

// forward declaration of enums not possible
#include "../drt_inpar/inpar_xfem.H"
#include "../drt_lib/drt_element.H"
#include "../drt_lib/drt_element_integration_select.H"
#include "../drt_cut/cut_utils.H"

// avoid forward declaration for templated class
#include "../drt_fluid_ele/fluid_ele_calc_xfem_coupling.H"

// forward declarations
namespace DRT
{
  namespace UTILS
  {
    class GaussIntegration;
  }
}  // namespace DRT

namespace GEO
{
  namespace CUT
  {
    class BoundaryCell;
  }
}  // namespace GEO

namespace XFEM
{
  class ConditionManager;
  namespace UTILS
  {
    //! @name GetAverageWeights
    /*!
    \brief Get the std - average weights kappa_m and kappa_s for the Nitsche calculations
     */
    void GetStdAverageWeights(
        const INPAR::XFEM::AveragingStrategy averaging_strategy, double& kappa_m);

    //! @name NIT_getTraceEstimateConstant
    /*!
    \brief get the constant which satisfies the trace inequality depending on the spatial dimension
    and polynomial order of the element
     */
    double NIT_getTraceEstimateConstant(
        const DRT::Element::DiscretizationType ele_distype, const bool is_pseudo_2D);


    //! @name NIT_Compute_ViscPenalty_Stabfac
    /*!
    \brief compute viscous part of Nitsche's penalty term scaling for Nitsche's method
     */
    void NIT_Compute_ViscPenalty_Stabfac(
        const DRT::Element::DiscretizationType
            ele_distype,  ///< the discretization type of the element w.r.t which the stabilization
                          ///< factor is computed
        const double&
            penscaling,  ///< material dependent penalty scaling (e.g. visceff) divided by h
        const double& NIT_stabscaling,  ///< basic nit penalty stab scaling
        const bool& is_pseudo_2D,       ///< is pseudo 2d
        const INPAR::XFEM::ViscStab_TraceEstimate&
            visc_stab_trace_estimate,  ///< how to estimate the scaling from the trace inequality
        double& NIT_visc_stab_fac      ///< viscous part of Nitsche's penalty term
    );

    //! @name GetNavierSlipStabilizationParameters
    /*!
    \brief Get NavierSlip Stabilization Parameters for tangential direction
     */
    void GetNavierSlipStabilizationParameters(
        const double& NIT_visc_stab_fac,  ///< viscous Nitsche stab fac
        double& dynvisc,                  ///< average dynamic viscosity
        double& sliplength,               ///< sliplength
        double& stabnit,                  ///< stabilization factor NIT_Penalty
        double& stabadj                   ///< stabilization factor Adjoint
    );

    //! compute transformation factor for surface integration, normal, local and global gp
    //! coordinates
    void ComputeSurfaceTransformation(double& drs,  ///< surface transformation factor
        LINALG::Matrix<3, 1>& x_gp_lin,             ///< global coordiantes of gaussian point
        LINALG::Matrix<3, 1>& normal,               ///< normal vector on boundary cell
        GEO::CUT::BoundaryCell* bc,                 ///< boundary cell
        const LINALG::Matrix<2, 1>& eta  ///< local coordinates of gaussian point w.r.t boundarycell
    );

    //! pre-compute the measure of the element's intersecting surface
    double ComputeMeasCutSurf(const std::map<int, std::vector<DRT::UTILS::GaussIntegration>>&
                                  bintpoints,  ///< boundary cell integration points
        const std::map<int, std::vector<GEO::CUT::BoundaryCell*>>& bcells  ///< boundary cells
    );

    //! compute the measure of the elements surface with given local id
    double ComputeMeasFace(DRT::Element* ele,  ///< fluid element
        Epetra_SerialDenseMatrix& ele_xyze,    ///< element coordinates
        const int local_face_id,               ///< the local id of the face w.r.t the fluid element
        const int nsd                          ///< number of space dimensions
    );

    //! compute volume-equivalent diameter
    inline double ComputeVolEqDiameter(double vol)
    {
      // get element length for tau_Mp/tau_C: volume-equival. diameter/sqrt(3)
      const double hk = std::pow((6. * vol / M_PI), (1.0 / 3.0)) / sqrt(3.0);

      return hk;
    }

    //! evaluate element volume
    template <DRT::Element::DiscretizationType distype>
    double EvalElementVolume(
        LINALG::Matrix<3, DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement> xyze,
        LINALG::Matrix<DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement, 1>*
            nurbs_weights = NULL,
        std::vector<Epetra_SerialDenseVector>* nurbs_knots = NULL);

    //! compute characteristic element length h_k
    template <DRT::Element::DiscretizationType distype>
    double ComputeCharEleLength(DRT::Element* ele,                 ///< fluid element
        Epetra_SerialDenseMatrix& ele_xyze,                        ///< element coordinates
        const Teuchos::RCP<XFEM::ConditionManager>& cond_manager,  ///< XFEM condition manager
        const GEO::CUT::plain_volumecell_set& vcSet,  ///< volumecell sets for volume integration
        const std::map<int, std::vector<GEO::CUT::BoundaryCell*>>&
            bcells,  ///< bcells for boundary cell integration
        const std::map<int, std::vector<DRT::UTILS::GaussIntegration>>&
            bintpoints,  ///< integration points for boundary cell integration
        const INPAR::XFEM::ViscStab_hk visc_stab_hk,  ///< h definition
        Teuchos::RCP<DRT::ELEMENTS::XFLUID::SlaveElementInterface<distype>> emb =
            Teuchos::null,         ///< pointer to the embedded coupling implementation
        DRT::Element* face = NULL  ///< side element in 3D
    );

    //! compute full scaling of Nitsche's penalty term (xfluid-fluid)
    void NIT_Compute_FullPenalty_Stabfac(
        double& NIT_full_stab_fac,           ///< to be filled: full Nitsche's penalty term scaling
                                             ///< (viscous+convective part)
        const LINALG::Matrix<3, 1>& normal,  ///< interface-normal vector
        const double h_k,                    ///< characteristic element length
        const double kappa1,                 ///< Weight parameter (parameter +/master side)
        const double kappa2,                 ///< Weight parameter (parameter -/slave  side)
        const LINALG::Matrix<3, 1>& velint_m,  ///< Master side velocity at gauss-point
        const LINALG::Matrix<3, 1>& velint_s,  ///< Slave side velocity at gauss-point
        const double NIT_visc_stab_fac,        ///< Nitsche's viscous scaling part of penalty term
        const double timefac,                  ///< timefac
        const bool isstationary,               ///< isstationary
        const double densaf_master,            ///< master density
        const double densaf_slave,             ///< slave density
        INPAR::XFEM::MassConservationScaling
            MassConservationScaling,  ///< kind of mass conservation scaling
        INPAR::XFEM::MassConservationCombination
            MassConservationCombination,  ///< kind of mass conservation combination
        const double NITStabScaling,      ///< scaling of nit stab fac
        INPAR::XFEM::ConvStabScaling
            ConvStabScaling,  ///< which convective stab. scaling of inflow stab
        INPAR::XFEM::XFF_ConvStabScaling
            XFF_ConvStabScaling,            ///< which convective stab. scaling on XFF interface
        const bool IsConservative = false,  ///< conservative formulation of navier stokes
        bool error_calc = false  ///< when called in error calculation, don't add the inflow terms
    );

  }  // namespace UTILS
}  // namespace XFEM

#endif /* XFEM_INTERFACE_UTILS_H_ */
