/*----------------------------------------------------------------------*/
/*!
\file xfem_interface_utils.H
\brief Basic routings to evaluate the terms for Nitsche Interface

\level 2

<pre>
\maintainer  Ager Christoph
             ager@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289-15249
</pre>

*/
/*----------------------------------------------------------------------*/

#ifndef XFEM_INTERFACE_UTILS_H_
#define XFEM_INTERFACE_UTILS_H_

 //forward declaration of enums not possible
#include "../drt_inpar/inpar_xfem.H"
#include "../drt_lib/drt_element.H"

namespace XFEM
{
  namespace UTILS
  {
    //! @name GetAverageWeights
    /*!
    \brief Get the std - average weights kappa_m and kappa_s for the Nitsche calculations
     */
    void GetStdAverageWeights(const INPAR::XFEM::AveragingStrategy averaging_strategy,
                              double & kappa_m);

    //! @name NIT_getTraceEstimateConstant
    /*!
    \brief get the constant which satisfies the trace inequality depending on the spatial dimension and polynomial order of the element
     */
    double NIT_getTraceEstimateConstant(const DRT::Element::DiscretizationType ele_distype,
                                        const bool is_pseudo_2D);


    //! @name NIT_Compute_ViscPenalty_Stabfac
    /*!
    \brief compute viscous part of Nitsche's penalty term scaling for Nitsche's method
     */
    void NIT_Compute_ViscPenalty_Stabfac(
        const DRT::Element::DiscretizationType ele_distype,                   ///< the discretization type of the element w.r.t which the stabilization factor is computed
        const double& penscaling,                                             ///< material dependent penalty scaling (e.g. visceff) divided by h
        const double& NIT_stabscaling,                                        ///< basic nit penalty stab scaling
        const bool& is_pseudo_2D,                                             ///< is pseudo 2d
        const INPAR::XFEM::ViscStab_TraceEstimate& visc_stab_trace_estimate,  ///< how to estimate the scaling from the trace inequality
        double& NIT_visc_stab_fac                                             ///< viscous part of Nitsche's penalty term
    );

    //! @name GetNavierSlipStabilizationParameters
    /*!
    \brief Get NavierSlip Stabilization Parameters for tangential direction
     */
    void GetNavierSlipStabilizationParameters(
      const double &                                NIT_full_stab_fac,            ///< full Nitsche stab fac
      const double &                                NIT_visc_stab_fac,            ///< viscous Nitsche stab fac
      double &                                      dynvisc,                      ///< average dynamic viscosity
      double &                                      sliplength,                   ///< sliplength
      double &                                      stabnit,                      ///< stabilization factor NIT_Penalty
      double &                                      stabadj                       ///< stabilization factor Adjoint
    );
  }
}

#endif /* XFEM_INTERFACE_UTILS_H_ */
