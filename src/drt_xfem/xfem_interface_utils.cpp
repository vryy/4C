/*!----------------------------------------------------------------------
\file xfem_interface_utils.cpp
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

#include "xfem_interface_utils.H"

/*----------------------------------------------------------------------*
 * Get the std - average weights kappa_m and kappa_s for the Nitsche calculations
 *----------------------------------------------------------------------*/
void XFEM::UTILS::GetStdAverageWeights(const INPAR::XFEM::AveragingStrategy averaging_strategy,
                                       double & kappa_m)
{
  switch(averaging_strategy)
  {
    case INPAR::XFEM::Xfluid_Sided:
    {
      kappa_m = 1.0;
      break;
    }
    case INPAR::XFEM::Embedded_Sided:
    {
      kappa_m = 0.0;
      break;
    }
    case INPAR::XFEM::Mean:
    {
      kappa_m = 0.5;
      break;
    }
    default:
      dserror("XFEM::UTILS::GetAverageWeights: AveragingStrategy not implemented here!");
  }

  return;
}

/*--------------------------------------------------------------------------------
 * compute viscous part of Nitsche's penalty term scaling for Nitsche's method
 *--------------------------------------------------------------------------------*/
void XFEM::UTILS::NIT_Compute_ViscPenalty_Stabfac(
    const DRT::Element::DiscretizationType ele_distype,                   ///< the discretization type of the element w.r.t which the stabilization factor is computed
    const double& penscaling,                                             ///< material dependent penalty scaling (e.g. visceff) divided by h
    const double& NIT_stabscaling,                                        ///< basic nit penalty stab scaling
    const bool& is_pseudo_2D,                                             ///< is pseudo 2d
    const INPAR::XFEM::ViscStab_TraceEstimate& visc_stab_trace_estimate,  ///< how to estimate the scaling from the trace inequality
    double& NIT_visc_stab_fac                                             ///< viscous part of Nitsche's penalty term
)
{
  //------------------------------
  // scaling factors for Nitsche's standard stabilization term
  /*
        viscous_Nitsche-part
      /                                        \        /                           i        \
     |  NIT_visc_stab_fac *  [ v ] , [ Du ]     | =  - |   NIT_visc_stab_fac * [ v ], [ u ]   |
      \                                        /        \                                    /


      with

            NIT_visc_stab_fac =   gamma * mu * C^2

      where
           - gamma is a dimensionless user-defined parameter
           - mu is the viscosity
           - C is an estimate of the following trace inequality with  C^2 ~ 1/h
             and C depends on element type, shape and polynomial degree


                    trace inequality to be satisfied for Nitsche's method

                    || grad v_h ||         <=   C  *  || grad v_h ||
                                  Gamma_K                           K

           - C^2 \approx lambda can be estimated as the maximal eigenvalue lambda of the generalized eigenvalue problem (A*x = lambda * B*x)

                    < grad w_h,  grad v_h >          =  lambda *  ( grad w_h, grad v_h)
                                           Gamma_K                                     K

           - alternatively we can estimate C^2 as C^2 \approx \rho_safety * CT / h_K (see generalized hp-framework by Erik Burman 2006)

                      with CT depending on the dimension d and the polynomial degree p

                               CT = p(p+1)/2 * (2+1/p)^d ( for tensor-product elements (hexahedral,quadrilateral...)
                               CT = (p+1)*(p+d)/d        ( for simplices (lines, triangles, tetrahedral...)

                      and  1/h_K \approx meas(Gamma_K)/meas(K)
                      and rho_safety a safety factor depending on the element-type which is useful when the surface cuts the element in an inclined situation or when the surface
                          contains small acute angles or corners
  */

  //------------------------------
  // 1 // to be filled viscous part of Nitsche's penalty term scaling for Nitsche's method
  // 2 // scale the viscous part of the penalty scaling with the effective viscosity of both sides
  // 3 // scale the viscous part of the penalty scaling with the dimensionless user defined Nitsche-parameter gamma
  NIT_visc_stab_fac = penscaling*NIT_stabscaling;

  // compute the final viscous scaling part of Nitsche's penalty term
  if(visc_stab_trace_estimate == INPAR::XFEM::ViscStab_TraceEstimate_CT_div_by_hk)
  {
    // get an estimate of the hp-depending constant C_T satisfying the trace inequality w.r.t the corresponding element (compare different weightings)
    double C_T = NIT_getTraceEstimateConstant(ele_distype,is_pseudo_2D);

    #if(0)
        // get a safety factor dependent on the the element shape for cut elements, since the surface can cut the
        // background element in a bad way when the surface includes sharp corners or small acute angles
        //TODO: still to be set
        double rho_safety = 1.0;
        NIT_visc_stab_fac *= rho_safety;
    #endif

    // build the final viscous scaling
    NIT_visc_stab_fac *= C_T;
  }
  else if(visc_stab_trace_estimate != INPAR::XFEM::ViscStab_TraceEstimate_eigenvalue)
  {
    dserror("unknown trace-inequality-estimate type for viscous part of Nitsche's penalty term");
  }
  return;
}

/*--------------------------------------------------------------------------------
 * get the constant which satisfies the trace inequality depending on the spatial dimension and polynomial order of the element
 *--------------------------------------------------------------------------------*/
double XFEM::UTILS::NIT_getTraceEstimateConstant(
    const DRT::Element::DiscretizationType ele_distype,
    bool is_pseudo_2D
)
{
  /*
  => return
       CT = p(p+1)/2 * (2+1/p)^d  OR an estimate obtained form solving an EVP ( for tensor-product elements (hexahedral,quadrilateral...)
       CT = (p+1)*(p+d)/d                                                     ( for simplices (lines, triangles, tetrahedral...)
     constant depending on element-type/polynomial order, see below

  - C is an estimate of the following trace inequality with  C^2 ~ 1/h
     and C depends on element type, shape and polynomial degree


            trace inequality to be satisfied for Nitsche's method

            || grad v_h ||         <=   C  *  || grad v_h ||
                          Gamma_K                           K

   - we can estimate C^2 as C^2 \approx \rho_safety * CT / h_K (see generalized hp-framework by Erik Burman 2006)

              with CT depending on the dimension d and the polynomial degree p of grad(v_h)

              - For tensor-product elements (hexahedral,quadrilateral,...)
                   CT = p(p+1)/2 * (2+1/p)^d
                      -> REMARK: the theoretical estimate of the constant CT for hexahedral and quadrilateral elements seems to overestimate the actual constant
                      -> therefore we take the approximate values from solving a local eigenvalue problem...)
              - for simplices (lines, triangles, tetrahedral...)
                   CT = (p+1)*(p+d)/d        ( this estimate is sharp and leads to the same constants as solving an EVP)

              and  1/h_K \approx meas(Gamma_K)/meas(K)
              and rho_safety a safety factor depending on the element-type which is useful when the surface cuts the element in an inclined situation or when the surface
                  contains small acute angles or corners

   - Remark regarding the dimension of the problem (pseudo-2D simulations):
     -> important for CT is the degree of the element-wise polynomial which is clear for pure 2D and 3D problems
     -> for pseudo-2D with one element-layer of 3D element in z-direction, but strong Dirichlet conditions which eliminate the z-contributions
        in the polynomial, we have to adapt the constant from d=3 to d=2 although a 3D element is used
        (use the "is_pseudo_2D"-flag)
        => this ensures same results of a 2D simulation and a pseudo-2D simulation using the 3D implementation in combination with
           strong Dirichlet conditions in z-direction
     -> an important property is that the elongation of the element in z-direction does not play a role when surface-element-measure based
        hk-definitions are used since then the relation meas(Gamma_K)/meas(K) remains constant
     -> using non-strong Dirichlet conditions for the enforcement of u_z = 0, CT for 3D is required which then leads
        to differences between a 2D implementation and pseudo-2D simulation using the 3D implementation
     -> IMPORTANT: to forget to set the "is_pseudo_2D"-flag to "true" for pseudo-2D examples leads to an overestimate of CT,
        and so does not result in stability problems, it might lead to a slightly worse convergence of the linear solver,
        however, usually it has not a to worse effect.
*/

  DRT::Element::DiscretizationType trace_inequality_distype;

  if(!is_pseudo_2D)
  {
    // the standard case for real 2D or 3D simulations
    trace_inequality_distype = ele_distype;
  }
  else
  {
    // modification for pseudo 2D simulations
    switch (ele_distype)
    {
    case DRT::Element::hex8:
      trace_inequality_distype = DRT::Element::quad4; break;    // hex8 -> quad4 reduction
    case DRT::Element::hex20:
      trace_inequality_distype = DRT::Element::quad8; break;    // hex20-> quad8 reduction
    case DRT::Element::hex27:
      trace_inequality_distype = DRT::Element::quad9; break;    // hex27-> quad9 reduction
    case DRT::Element::wedge15:
      trace_inequality_distype = DRT::Element::tri6; break;     // wedge15 -> tri6 reduction (tri6 elements in 2D plane)
    case DRT::Element::wedge6:
      trace_inequality_distype = DRT::Element::tri3; break;     // wedge6 -> tri3 reduction (tri3 elements in 2D plane)
    default:
    {
      dserror("not a valid pseudo 2D element-type - what to do?");
      trace_inequality_distype = DRT::Element::dis_none;
      break;
    }
    };
  }

//  return 2.0;

  // switch over the distype which determines the right polynomial degree and dimension
  switch (trace_inequality_distype)
  {
  // for triangluar/tetradhedral elements:
  //    -> the grad-operator reduces the polynomial degree p(grad(v_h)) = p(v_h)-1
  // for quadrilateral/hexahedral and wedge/pyramid elements:
  //    -> the grad-operator does not reduce the polynomial degree due to the mixed polynomials p(grad(v_h)) = p(v_h)

  /*
  // CT = p(p+1)/2 * (2+1/p)^d
  // 3D hexahedral elements
  case DRT::Element::hex8:      return 27.0;                 break;  /// d=3, p(v_h) = 1 -> p(grad(v_h)) = 1   =>  CT=27.0
  case DRT::Element::hex20:     return 46.875;               break;  /// d=3, p(v_h) = 2 -> p(grad(v_h)) = 2   =>  CT=46.875 (375/8)
  case DRT::Element::hex27:     return 46.875;               break;  /// d=3, p(v_h) = 2 -> p(grad(v_h)) = 2   =>  CT=46.875 (375/8)
  // 2D quadrilateral elements
  case DRT::Element::quad4:     return 9.0;                  break;  /// d=2, p(v_h) = 1 -> p(grad(v_h)) = 1   =>  CT=9.0
  case DRT::Element::quad8:     return 9.0;                  break;  /// d=2, p(v_h) = 1 -> p(grad(v_h)) = 1   =>  CT=9.0
  case DRT::Element::quad9:     return 18.75;                break;  /// d=2, p(v_h) = 2 -> p(grad(v_h)) = 2   =>  CT=18.75  (75/4)
  */
  // REMARK: the theroretical estimate of the constant CT for hexahedral and quadrilateral elements seems to overestimate the actual constant
  // -> therefore we take the approximate values from solving a local eigenvalue problem
  // estimates from solving the eigenvalue problem on regular elements
  case DRT::Element::hex8:      return 1.59307;              break;  /// d=3, p(v_h) = 1 -> p(grad(v_h)) = 1   =>  CT= from eigenvalue problem
  case DRT::Element::hex20:     return 4.10462;              break;  /// d=3, p(v_h) = 2 -> p(grad(v_h)) = 2   =>  CT= from eigenvalue problem
  case DRT::Element::hex27:     return 4.27784;              break;  /// d=3, p(v_h) = 2 -> p(grad(v_h)) = 2   =>  CT= from eigenvalue problem
  // 2D quadrilateral elements
  case DRT::Element::quad4:     return 1.43426;              break;  /// d=2, p(v_h) = 1 -> p(grad(v_h)) = 1   =>  CT= from eigenvalue problem
  case DRT::Element::quad8:     return 4.06462;              break;  /// d=2, p(v_h) = 1 -> p(grad(v_h)) = 1   =>  CT= from eigenvalue problem
  case DRT::Element::quad9:     return 4.19708;              break;  /// d=2, p(v_h) = 2 -> p(grad(v_h)) = 2   =>  CT= from eigenvalue problem
  //----------------------------------------------------------------
  // CT = (p+1)*(p+d)/d (this estimate leads to the same results as solving a local eigenvalue problem)
  // 3D tetrahedral elements
  case DRT::Element::tet4:      return 1.0;                  break;  /// d=3, p(v_h) = 1 -> p(grad(v_h)) = 0   =>  CT=1.0
  case DRT::Element::tet10:     return 2.6666666666666666;   break;  /// d=3, p(v_h) = 2 -> p(grad(v_h)) = 1   =>  CT=2.6666666666666666 (8/3)
  // 2D triangular elements
  case DRT::Element::tri3:      return 1.0;                  break;  /// d=2, p(v_h) = 1 -> p(grad(v_h)) = 0   =>  CT=1.0
  case DRT::Element::tri6:      return 3.0;                  break;  /// d=2, p(v_h) = 2 -> p(grad(v_h)) = 1   =>  CT=3.0
  // 1D line elements
  case DRT::Element::line2:     return 1.0;                  break;  /// d=1, p(v_h) = 1 -> p(grad(v_h)) = 0   =>  CT=1.0
  case DRT::Element::line3:     return 4.0;                  break;  /// d=1, p(v_h) = 1 -> p(grad(v_h)) = 0   =>  CT=4.0
  //----------------------------------------------------------------
  // 3D wedge/pyramid elements, the current estimates are taken from the maximum of hex and tet elements, the correct value has to switch between the faces!
  case DRT::Element::pyramid5:  std::cout << "WARNING: calibrate this value!" << std::endl; return 1.59307;              break;  /// d=3, p(v_h) = 1 -> p(grad(v_h)) = 1   =>  CT taken from hex8
  case DRT::Element::wedge6:    std::cout << "WARNING: calibrate this value!" << std::endl; return 1.59307;              break;  /// d=3, p(v_h) = 1 -> p(grad(v_h)) = 1   =>  CT taken from hex8
  case DRT::Element::wedge15:   std::cout << "WARNING: calibrate this value!" << std::endl; return 4.10462;              break;  /// d=3, p(v_h) = 2 -> p(grad(v_h)) = 2   =>  CT taken from hex20
  default:
    dserror("constant for trace inequality not specified for this element type yet: % i", trace_inequality_distype); break;
  };

  return 0.0;
}

/*----------------------------------------------------------------------*
 * Get Stabilization Parameters (for Penalty and Adjoint)
 * for Navier Slip Formulation
 *----------------------------------------------------------------------*/
void XFEM::UTILS::GetNavierSlipStabilizationParameters(
  const double &                                NIT_visc_stab_fac,            ///< viscous Nitsche stab fac
  double &                                      dynvisc,                      ///< average dynamic viscosity
  double &                                      sliplength,                   ///< sliplength
  double &                                      stabnit,                      ///< stabilization factor NIT_Penalty
  double &                                      stabadj                       ///< stabilization factor Adjoint
)
{
  // Create stabilization parameters needed for the Robin case
  //  NIT_visc_stab_fac     =   gamma' * mu * C^2
  //  NIT_visc_stab_fac_inv =   ( (1/gamma)*h_E )

  double NIT_visc_stab_fac_inv;

  if(NIT_visc_stab_fac<=0.0)
    NIT_visc_stab_fac_inv = 1e15; //If Nitsche parameter is 0
  else
    NIT_visc_stab_fac_inv = 1.0/(NIT_visc_stab_fac/dynvisc);

  //  NIT_robin_denominator = [ mu/(epislon + gamma*h_E) ]
  double NIT_robin_denominator_no_mu = 1.0/(sliplength+NIT_visc_stab_fac_inv);
  //NIT_robin_denominator_ = dyn_visc_*NIT_robin_denominator_no_mu;

  // Nitsche penalty term stabilization in tangential direction
  // ----------------------------------------------------------------------
  // stabnit =  [ { mu }/(epsilon + gamma*h_E) ) + extra_terms]
  stabnit    = dynvisc*NIT_robin_denominator_no_mu; //NIT_robin_denominator_;

  //  stabepsnit = [ epsilon / (epsilon + gamma*h_E) + extra_terms ];
 // stabepsnit = NIT_robin_denominator_no_mu*sliplength;

  // ----------------------------------------------------------------------

  // Adjoint-terms stabilization in tangential direction
  // ----------------------------------------------------------------------
  //  stabadj = [ gamma*h_E  /(epsilon + gamma*h_E)  ]
  stabadj    = NIT_robin_denominator_no_mu*NIT_visc_stab_fac_inv;

  #ifdef ENFORCE_URQUIZA_GNBC
  // In tangential direction stabilization:
  stabnit    = (dynvisc/sliplength);  //Appears in NIT-penalty stabilization
  stabadj    = 0.0;
  stabepsnit = 0.0;
  stabepsadj = 0.0;
  #endif

  return;
}
