/*----------------------------------------------------------------------*/
/*! \file
\brief Basic routings to evaluate the terms for Nitsche Interface

\level 2

\maintainer  Christoph Ager
             ager@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289-15249

*/
/*----------------------------------------------------------------------*/

#include "xfem_interface_utils.H"
#include "xfem_condition_manager.H"

#include "../drt_cut/cut_boundarycell.H"
#include "../drt_cut/cut_volumecell.H"

#include "../drt_fem_general/drt_utils_local_connectivity_matrices.H"
#include "../drt_fem_general/drt_utils_integration.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_fem_general/drt_utils_nurbs_shapefunctions.H"
//

/*----------------------------------------------------------------------*
 * Get the std - average weights kappa_m and kappa_s for the Nitsche calculations
 *----------------------------------------------------------------------*/
void XFEM::UTILS::GetStdAverageWeights(
    const INPAR::XFEM::AveragingStrategy averaging_strategy, double &kappa_m)
{
  switch (averaging_strategy)
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
    const DRT::Element::DiscretizationType
        ele_distype,  ///< the discretization type of the element w.r.t which the stabilization
                      ///< factor is computed
    const double &penscaling,  ///< material dependent penalty scaling (e.g. visceff) divided by h
    const double &NIT_stabscaling,  ///< basic nit penalty stab scaling
    const bool &is_pseudo_2D,       ///< is pseudo 2d
    const INPAR::XFEM::ViscStab_TraceEstimate
        &visc_stab_trace_estimate,  ///< how to estimate the scaling from the trace inequality
    double &NIT_visc_stab_fac       ///< viscous part of Nitsche's penalty term
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

           - C^2 \approx lambda can be estimated as the maximal eigenvalue lambda of the generalized
     eigenvalue problem (A*x = lambda * B*x)

                    < grad w_h,  grad v_h >          =  lambda *  ( grad w_h, grad v_h)
                                           Gamma_K                                     K

           - alternatively we can estimate C^2 as C^2 \approx \rho_safety * CT / h_K (see
     generalized hp-framework by Erik Burman 2006)

                      with CT depending on the dimension d and the polynomial degree p

                               CT = p(p+1)/2 * (2+1/p)^d ( for tensor-product elements
     (hexahedral,quadrilateral...) CT = (p+1)*(p+d)/d        ( for simplices (lines, triangles,
     tetrahedral...)

                      and  1/h_K \approx meas(Gamma_K)/meas(K)
                      and rho_safety a safety factor depending on the element-type which is useful
     when the surface cuts the element in an inclined situation or when the surface contains small
     acute angles or corners
  */

  //------------------------------
  // 1 // to be filled viscous part of Nitsche's penalty term scaling for Nitsche's method
  // 2 // scale the viscous part of the penalty scaling with the effective viscosity of both sides
  // 3 // scale the viscous part of the penalty scaling with the dimensionless user defined
  // Nitsche-parameter gamma
  NIT_visc_stab_fac = penscaling * NIT_stabscaling;

  // compute the final viscous scaling part of Nitsche's penalty term
  if (visc_stab_trace_estimate == INPAR::XFEM::ViscStab_TraceEstimate_CT_div_by_hk)
  {
    // get an estimate of the hp-depending constant C_T satisfying the trace inequality w.r.t the
    // corresponding element (compare different weightings)
    double C_T = NIT_getTraceEstimateConstant(ele_distype, is_pseudo_2D);

#if (0)
    // get a safety factor dependent on the the element shape for cut elements, since the surface
    // can cut the background element in a bad way when the surface includes sharp corners or small
    // acute angles
    // TODO: still to be set
    double rho_safety = 1.0;
    NIT_visc_stab_fac *= rho_safety;
#endif

    // build the final viscous scaling
    NIT_visc_stab_fac *= C_T;
  }
  else if (visc_stab_trace_estimate != INPAR::XFEM::ViscStab_TraceEstimate_eigenvalue)
  {
    dserror("unknown trace-inequality-estimate type for viscous part of Nitsche's penalty term");
  }
  return;
}

/*--------------------------------------------------------------------------------
 * get the constant which satisfies the trace inequality depending on the spatial dimension and
 *polynomial order of the element
 *--------------------------------------------------------------------------------*/
double XFEM::UTILS::NIT_getTraceEstimateConstant(
    const DRT::Element::DiscretizationType ele_distype, bool is_pseudo_2D)
{
  /*
  => return
       CT = p(p+1)/2 * (2+1/p)^d  OR an estimate obtained form solving an EVP ( for tensor-product
  elements (hexahedral,quadrilateral...) CT = (p+1)*(p+d)/d ( for simplices (lines, triangles,
  tetrahedral...) constant depending on element-type/polynomial order, see below

  - C is an estimate of the following trace inequality with  C^2 ~ 1/h
     and C depends on element type, shape and polynomial degree


            trace inequality to be satisfied for Nitsche's method

            || grad v_h ||         <=   C  *  || grad v_h ||
                          Gamma_K                           K

   - we can estimate C^2 as C^2 \approx \rho_safety * CT / h_K (see generalized hp-framework by Erik
  Burman 2006)

              with CT depending on the dimension d and the polynomial degree p of grad(v_h)

              - For tensor-product elements (hexahedral,quadrilateral,...)
                   CT = p(p+1)/2 * (2+1/p)^d
                      -> REMARK: the theoretical estimate of the constant CT for hexahedral and
  quadrilateral elements seems to overestimate the actual constant
                      -> therefore we take the approximate values from solving a local eigenvalue
  problem...)
              - for simplices (lines, triangles, tetrahedral...)
                   CT = (p+1)*(p+d)/d        ( this estimate is sharp and leads to the same
  constants as solving an EVP)

              and  1/h_K \approx meas(Gamma_K)/meas(K)
              and rho_safety a safety factor depending on the element-type which is useful when the
  surface cuts the element in an inclined situation or when the surface contains small acute angles
  or corners

   - Remark regarding the dimension of the problem (pseudo-2D simulations):
     -> important for CT is the degree of the element-wise polynomial which is clear for pure 2D and
  3D problems
     -> for pseudo-2D with one element-layer of 3D element in z-direction, but strong Dirichlet
  conditions which eliminate the z-contributions in the polynomial, we have to adapt the constant
  from d=3 to d=2 although a 3D element is used (use the "is_pseudo_2D"-flag)
        => this ensures same results of a 2D simulation and a pseudo-2D simulation using the 3D
  implementation in combination with strong Dirichlet conditions in z-direction
     -> an important property is that the elongation of the element in z-direction does not play a
  role when surface-element-measure based hk-definitions are used since then the relation
  meas(Gamma_K)/meas(K) remains constant
     -> using non-strong Dirichlet conditions for the enforcement of u_z = 0, CT for 3D is required
  which then leads to differences between a 2D implementation and pseudo-2D simulation using the 3D
  implementation
     -> IMPORTANT: to forget to set the "is_pseudo_2D"-flag to "true" for pseudo-2D examples leads
  to an overestimate of CT, and so does not result in stability problems, it might lead to a
  slightly worse convergence of the linear solver, however, usually it has not a to worse effect.
*/

  DRT::Element::DiscretizationType trace_inequality_distype;

  if (!is_pseudo_2D)
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
        trace_inequality_distype = DRT::Element::quad4;
        break;  // hex8 -> quad4 reduction
      case DRT::Element::hex20:
        trace_inequality_distype = DRT::Element::quad8;
        break;  // hex20-> quad8 reduction
      case DRT::Element::hex27:
        trace_inequality_distype = DRT::Element::quad9;
        break;  // hex27-> quad9 reduction
      case DRT::Element::wedge15:
        trace_inequality_distype = DRT::Element::tri6;
        break;  // wedge15 -> tri6 reduction (tri6 elements in 2D plane)
      case DRT::Element::wedge6:
        trace_inequality_distype = DRT::Element::tri3;
        break;  // wedge6 -> tri3 reduction (tri3 elements in 2D plane)
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
    //    -> the grad-operator does not reduce the polynomial degree due to the mixed polynomials
    //    p(grad(v_h)) = p(v_h)

    /*
    // CT = p(p+1)/2 * (2+1/p)^d
    // 3D hexahedral elements
    case DRT::Element::hex8:      return 27.0;                 break;  /// d=3, p(v_h) = 1 ->
    p(grad(v_h)) = 1   =>  CT=27.0 case DRT::Element::hex20:     return 46.875;               break;
    /// d=3, p(v_h) = 2 -> p(grad(v_h)) = 2   =>  CT=46.875 (375/8) case DRT::Element::hex27:
    return 46.875;               break;  /// d=3, p(v_h) = 2 -> p(grad(v_h)) = 2   =>  CT=46.875
    (375/8)
    // 2D quadrilateral elements
    case DRT::Element::quad4:     return 9.0;                  break;  /// d=2, p(v_h) = 1 ->
    p(grad(v_h)) = 1   =>  CT=9.0 case DRT::Element::quad8:     return 9.0;                  break;
    /// d=2, p(v_h) = 1 -> p(grad(v_h)) = 1   =>  CT=9.0 case DRT::Element::quad9:     return 18.75;
    break;  /// d=2, p(v_h) = 2 -> p(grad(v_h)) = 2   =>  CT=18.75  (75/4)
    */
    // REMARK: the theroretical estimate of the constant CT for hexahedral and quadrilateral
    // elements seems to overestimate the actual constant
    // -> therefore we take the approximate values from solving a local eigenvalue problem
    // estimates from solving the eigenvalue problem on regular elements
    case DRT::Element::hex8:
      return 1.59307;
      break;  /// d=3, p(v_h) = 1 -> p(grad(v_h)) = 1   =>  CT= from eigenvalue problem
    case DRT::Element::hex20:
      return 4.10462;
      break;  /// d=3, p(v_h) = 2 -> p(grad(v_h)) = 2   =>  CT= from eigenvalue problem
    case DRT::Element::hex27:
      return 4.27784;
      break;  /// d=3, p(v_h) = 2 -> p(grad(v_h)) = 2   =>  CT= from eigenvalue problem
    // 2D quadrilateral elements
    case DRT::Element::quad4:
      return 1.43426;
      break;  /// d=2, p(v_h) = 1 -> p(grad(v_h)) = 1   =>  CT= from eigenvalue problem
    case DRT::Element::quad8:
      return 4.06462;
      break;  /// d=2, p(v_h) = 1 -> p(grad(v_h)) = 1   =>  CT= from eigenvalue problem
    case DRT::Element::quad9:
      return 4.19708;
      break;  /// d=2, p(v_h) = 2 -> p(grad(v_h)) = 2   =>  CT= from eigenvalue problem
    //----------------------------------------------------------------
    // CT = (p+1)*(p+d)/d (this estimate leads to the same results as solving a local eigenvalue
    // problem) 3D tetrahedral elements
    case DRT::Element::tet4:
      return 1.0;
      break;  /// d=3, p(v_h) = 1 -> p(grad(v_h)) = 0   =>  CT=1.0
    case DRT::Element::tet10:
      return 2.6666666666666666;
      break;  /// d=3, p(v_h) = 2 -> p(grad(v_h)) = 1   =>  CT=2.6666666666666666 (8/3)
    // 2D triangular elements
    case DRT::Element::tri3:
      return 1.0;
      break;  /// d=2, p(v_h) = 1 -> p(grad(v_h)) = 0   =>  CT=1.0
    case DRT::Element::tri6:
      return 3.0;
      break;  /// d=2, p(v_h) = 2 -> p(grad(v_h)) = 1   =>  CT=3.0
    // 1D line elements
    case DRT::Element::line2:
      return 1.0;
      break;  /// d=1, p(v_h) = 1 -> p(grad(v_h)) = 0   =>  CT=1.0
    case DRT::Element::line3:
      return 4.0;
      break;  /// d=1, p(v_h) = 1 -> p(grad(v_h)) = 0   =>  CT=4.0
    //----------------------------------------------------------------
    // 3D wedge/pyramid elements, the current estimates are taken from the maximum of hex and tet
    // elements, the correct value has to switch between the faces!
    case DRT::Element::pyramid5:
      std::cout << "WARNING: calibrate this value!" << std::endl;
      return 1.59307;
      break;  /// d=3, p(v_h) = 1 -> p(grad(v_h)) = 1   =>  CT taken from hex8
    case DRT::Element::wedge6:
      std::cout << "WARNING: calibrate this value!" << std::endl;
      return 1.59307;
      break;  /// d=3, p(v_h) = 1 -> p(grad(v_h)) = 1   =>  CT taken from hex8
    case DRT::Element::wedge15:
      std::cout << "WARNING: calibrate this value!" << std::endl;
      return 4.10462;
      break;  /// d=3, p(v_h) = 2 -> p(grad(v_h)) = 2   =>  CT taken from hex20
    default:
      dserror("constant for trace inequality not specified for this element type yet: % i",
          trace_inequality_distype);
      break;
  };

  return 0.0;
}

/*----------------------------------------------------------------------*
 * Get Stabilization Parameters (for Penalty and Adjoint)
 * for Navier Slip Formulation
 *----------------------------------------------------------------------*/
void XFEM::UTILS::GetNavierSlipStabilizationParameters(
    const double &NIT_visc_stab_fac,  ///< viscous Nitsche stab fac
    double &dynvisc,                  ///< average dynamic viscosity
    double &sliplength,               ///< sliplength
    double &stabnit,                  ///< stabilization factor NIT_Penalty
    double &stabadj                   ///< stabilization factor Adjoint
)
{
  // Create stabilization parameters needed for the Robin case
  //  NIT_visc_stab_fac     =   gamma' * mu * C^2
  //  NIT_visc_stab_fac_inv =   ( (1/gamma)*h_E )

  double NIT_visc_stab_fac_inv;

  if (NIT_visc_stab_fac <= 0.0)
    NIT_visc_stab_fac_inv = 1e15;  // If Nitsche parameter is 0
  else
    NIT_visc_stab_fac_inv = 1.0 / (NIT_visc_stab_fac / dynvisc);

  //  NIT_robin_denominator = [ mu/(epislon + gamma*h_E) ]
  double NIT_robin_denominator_no_mu = 1.0 / (sliplength + NIT_visc_stab_fac_inv);
  // NIT_robin_denominator_ = dyn_visc_*NIT_robin_denominator_no_mu;

  // Nitsche penalty term stabilization in tangential direction
  // ----------------------------------------------------------------------
  // stabnit =  [ { mu }/(epsilon + gamma*h_E) ) + extra_terms]
  stabnit = dynvisc * NIT_robin_denominator_no_mu;  // NIT_robin_denominator_;

  //  stabepsnit = [ epsilon / (epsilon + gamma*h_E) + extra_terms ];
  // stabepsnit = NIT_robin_denominator_no_mu*sliplength;

  // ----------------------------------------------------------------------

  // Adjoint-terms stabilization in tangential direction
  // ----------------------------------------------------------------------
  //  stabadj = [ gamma*h_E  /(epsilon + gamma*h_E)  ]
  stabadj = NIT_robin_denominator_no_mu * NIT_visc_stab_fac_inv;

#ifdef ENFORCE_URQUIZA_GNBC
  // In tangential direction stabilization:
  stabnit = (dynvisc / sliplength);  // Appears in NIT-penalty stabilization
  stabadj = 0.0;
  stabepsnit = 0.0;
  stabepsadj = 0.0;
#endif

  return;
}

/*--------------------------------------------------------------------------------
 * compute transformation factor for surface integration, normal, local and global gp coordinates
 *--------------------------------------------------------------------------------*/
void XFEM::UTILS::ComputeSurfaceTransformation(double &drs,  ///< surface transformation factor
    LINALG::Matrix<3, 1> &x_gp_lin,  ///< global coordiantes of gaussian point
    LINALG::Matrix<3, 1> &normal,    ///< normal vector on boundary cell
    GEO::CUT::BoundaryCell *bc,      ///< boundary cell
    const LINALG::Matrix<2, 1> &eta  ///< local coordinates of gaussian point w.r.t boundarycell
)
{
  normal.Clear();

  // get normal vector on linearized boundary cell, x-coordinates of gaussian point and surface
  // transformation factor
  switch (bc->Shape())
  {
    case DRT::Element::tri3:
    {
      bc->Transform<DRT::Element::tri3>(eta, x_gp_lin, normal, drs);
      break;
    }
    case DRT::Element::quad4:
    {
      bc->Transform<DRT::Element::quad4>(eta, x_gp_lin, normal, drs);
      break;
    }
    default:
      throw std::runtime_error("unsupported integration cell type");
  }

  return;
}

/*--------------------------------------------------------------------------------
 * pre-compute the measure of all side's surface cutting the element
 *--------------------------------------------------------------------------------*/
double XFEM::UTILS::ComputeMeasCutSurf(
    const std::map<int, std::vector<DRT::UTILS::GaussIntegration>>
        &bintpoints,  ///< boundary cell integration points
    const std::map<int, std::vector<GEO::CUT::BoundaryCell *>> &bcells  ///< boundary cells
)
{
  double surf = 0.0;

  //--------------------------------------------
  // loop intersecting sides
  // map of side-element id and Gauss points
  for (std::map<int, std::vector<DRT::UTILS::GaussIntegration>>::const_iterator i =
           bintpoints.begin();
       i != bintpoints.end(); ++i)
  {
    int sid = i->first;
    const std::vector<DRT::UTILS::GaussIntegration> &cutintpoints = i->second;

    // get side's boundary cells
    std::map<int, std::vector<GEO::CUT::BoundaryCell *>>::const_iterator j = bcells.find(sid);
    if (j == bcells.end()) dserror("missing boundary cell");

    const std::vector<GEO::CUT::BoundaryCell *> &bcs = j->second;
    if (bcs.size() != cutintpoints.size()) dserror("boundary cell integration rules mismatch");

    //--------------------------------------------
    // loop boundary cells w.r.t current cut side
    //--------------------------------------------
    for (std::vector<DRT::UTILS::GaussIntegration>::const_iterator i = cutintpoints.begin();
         i != cutintpoints.end(); ++i)
    {
      const DRT::UTILS::GaussIntegration &gi = *i;
      GEO::CUT::BoundaryCell *bc =
          bcs[i - cutintpoints.begin()];  // get the corresponding boundary cell

      //--------------------------------------------
      // loop gausspoints w.r.t current boundary cell
      //--------------------------------------------
      for (DRT::UTILS::GaussIntegration::iterator iquad = gi.begin(); iquad != gi.end(); ++iquad)
      {
        double drs =
            0.0;  // transformation factor between reference cell and linearized boundary cell

        const LINALG::Matrix<2, 1> eta(iquad.Point());  // xi-coordinates with respect to side

        LINALG::Matrix<3, 1> normal(true);

        LINALG::Matrix<3, 1> x_gp_lin(true);  // gp in xyz-system on linearized interface

        // compute transformation factor, normal vector and global Gauss point coordiantes
        if (bc->Shape() != DRT::Element::dis_none)  // Tessellation approach
        {
          XFEM::UTILS::ComputeSurfaceTransformation(drs, x_gp_lin, normal, bc, eta);
        }
        else  // MomentFitting approach
        {
          drs = 1.0;
          normal = bc->GetNormalVector();
          const double *gpcord = iquad.Point();
          for (int idim = 0; idim < 3; ++idim)
          {
            x_gp_lin(idim, 0) = gpcord[idim];
          }
        }

        const double surf_fac = drs * iquad.Weight();

        surf += surf_fac;

      }  // loop gausspoints w.r.t current boundary cell
    }    // loop boundary cells
  }      // loop intersecting sides

  return surf;
}

/*--------------------------------------------------------------------------------
 * compute the measure of the elements surface with given local id
 *--------------------------------------------------------------------------------*/
double XFEM::UTILS::ComputeMeasFace(DRT::Element *ele,  ///< fluid element
    Epetra_SerialDenseMatrix &ele_xyze,                 ///< element coordinates
    const int local_face_id,  ///< the local id of the face w.r.t the fluid element
    const int nsd             ///< number of space dimensions
)
{
  // get the shape of the face
  DRT::Element::DiscretizationType face_shape =
      DRT::UTILS::getEleFaceShapeType(ele->Shape(), local_face_id);

  // get the current node coordinates, extract them from the element's node coordinates
  const int numnode_face = DRT::UTILS::getNumberOfElementNodes(face_shape);
  Epetra_SerialDenseMatrix xyze_face(nsd, numnode_face);

  // map for numbering of nodes of the surfaces
  std::vector<std::vector<int>> map = DRT::UTILS::getEleNodeNumberingFaces(ele->Shape());

  // extract the surface's node coordinates from the element's nodes coordinates
  for (int n = 0; n < numnode_face; ++n)
  {
    const int node_lid = map[local_face_id][n];
    for (int idim = 0; idim < nsd; ++idim) xyze_face(idim, n) = ele_xyze(idim, node_lid);
  }

  // the metric tensor and the area of an infintesimal surface element
  Epetra_SerialDenseMatrix metrictensor(nsd - 1, nsd - 1);
  double drs = 0.0;

  if (nsd != 3)
    dserror("don't call this function for non-3D examples, adapt the following for 2D!");

  DRT::UTILS::GaussRule2D gaussrule = DRT::UTILS::intrule2D_undefined;
  switch (face_shape)
  {
    case DRT::Element::quad4:
    case DRT::Element::quad8:
    case DRT::Element::quad9:
      gaussrule = DRT::UTILS::intrule_quad_1point;
      break;
    case DRT::Element::tri3:
    case DRT::Element::tri6:
      gaussrule = DRT::UTILS::intrule_tri_1point;
      break;
    default:
      dserror("shape type unknown!\n");
      break;
  }

  double meas_face = 0.0;

  /*----------------------------------------------------------------------*
    |               start loop over integration points                     |
   *----------------------------------------------------------------------*/
  const DRT::UTILS::IntegrationPoints2D intpoints(gaussrule);
  for (int gpid = 0; gpid < intpoints.nquad; ++gpid)
  {
    const double e0 = intpoints.qxg[gpid][0];
    const double e1 = intpoints.qxg[gpid][1];

    Epetra_SerialDenseMatrix deriv(nsd - 1, numnode_face);

    // get shape functions and derivatives in the plane of the element
    DRT::UTILS::shape_function_2D_deriv1(deriv, e0, e1, face_shape);

    // compute measure tensor for surface element and the infinitesimal
    // area element drs for the integration
    DRT::UTILS::ComputeMetricTensorForSurface(xyze_face, deriv, metrictensor, &drs);

    meas_face += intpoints.qwgt[gpid] * drs;
  }

  return meas_face;
}

/*----------------------------------------------------------------------*
 | evaluate element volume                                              |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
double XFEM::UTILS::EvalElementVolume(
    LINALG::Matrix<3, DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement> xyze,
    LINALG::Matrix<DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement, 1>
        *nurbs_weights,
    std::vector<Epetra_SerialDenseVector> *nurbs_knots)
{
  static const int nsd = 3;
  static const int nen = DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;

  // use one-point Gauss rule
  DRT::UTILS::IntPointsAndWeights<nsd> intpoints_stab(
      DRT::ELEMENTS::DisTypeToStabGaussRule<distype>::rule);

  const double *gpcoord = intpoints_stab.IP().qxg[0];   // actual integration point (coords)
  const double gpweight = intpoints_stab.IP().qwgt[0];  // actual integration point (weight)

  LINALG::Matrix<nsd, 1> xsi(gpcoord, true);
  static LINALG::Matrix<nen, 1> funct;
  static LINALG::Matrix<nsd, nen> deriv;
  static LINALG::Matrix<nsd, nsd> xjm;
  static LINALG::Matrix<nsd, nsd> xji;

  switch (distype)
  {
    case DRT::Element::hex8:
    case DRT::Element::hex20:
    case DRT::Element::hex27:
    case DRT::Element::tet4:
    case DRT::Element::tet10:
    case DRT::Element::wedge6:
    case DRT::Element::wedge15:
    {
      // shape functions and their first derivatives
      DRT::UTILS::shape_function<distype>(xsi, funct);
      DRT::UTILS::shape_function_deriv1<distype>(xsi, deriv);
      break;
    }
    case DRT::Element::nurbs8:
    case DRT::Element::nurbs27:
    {
      if (nurbs_weights == NULL || nurbs_knots == NULL)
        dserror("For Nurbs elements, weights and knots are required!");
      DRT::NURBS::UTILS::nurbs_get_funct_deriv(
          funct, deriv, xsi, *nurbs_knots, *nurbs_weights, distype);
      break;
    }
    default:
      dserror("Distype not handled yet!");
  }

  //

  // get Jacobian matrix and determinant
  // actually compute its transpose....
  /*
    +-            -+ T      +-            -+
    | dx   dx   dx |        | dx   dy   dz |
    | --   --   -- |        | --   --   -- |
    | dr   ds   dt |        | dr   dr   dr |
    |              |        |              |
    | dy   dy   dy |        | dx   dy   dz |
    | --   --   -- |   =    | --   --   -- |
    | dr   ds   dt |        | ds   ds   ds |
    |              |        |              |
    | dz   dz   dz |        | dx   dy   dz |
    | --   --   -- |        | --   --   -- |
    | dr   ds   dt |        | dt   dt   dt |
    +-            -+        +-            -+
  */
  xjm.MultiplyNT(deriv, xyze);
  double det = xji.Invert(xjm);

  if (det < 1E-16) dserror("GLOBAL ELEMENT ZERO OR NEGATIVE JACOBIAN DETERMINANT: %f", det);

  // compute integration factor
  return gpweight * det;
}

/*--------------------------------------------------------------------------------
 * compute characteristic element length
 *--------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
double XFEM::UTILS::ComputeCharEleLength(DRT::Element *ele,    ///< fluid element
    Epetra_SerialDenseMatrix &ele_xyze,                        ///< element coordinates
    const Teuchos::RCP<XFEM::ConditionManager> &cond_manager,  ///< XFEM condition manager
    const GEO::CUT::plain_volumecell_set &vcSet,  ///< volumecell sets for volume integration
    const std::map<int, std::vector<GEO::CUT::BoundaryCell *>>
        &bcells,  ///< bcells for boundary cell integration
    const std::map<int, std::vector<DRT::UTILS::GaussIntegration>>
        &bintpoints,  ///< integration points for boundary cell integration
    const INPAR::XFEM::ViscStab_hk visc_stab_hk,  ///< h definition
    Teuchos::RCP<DRT::ELEMENTS::XFLUID::SlaveElementInterface<distype>>
        emb,            ///< pointer to the embedded coupling implementation
    DRT::Element *face  ///< side element in 3D
)
{
  static const int nsd = 3;
  // TEUCHOS_FUNC_TIME_MONITOR("FluidEleCalcXFEM::ComputeCharEleLength");

  LINALG::Matrix<3, DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement> xyze(
      ele_xyze, true);
  const int coup_sid = bintpoints.begin()->first;
  const INPAR::XFEM::AveragingStrategy averaging_strategy =
      cond_manager->GetAveragingStrategy(coup_sid, ele->Id());
  if (emb == Teuchos::null and averaging_strategy == INPAR::XFEM::Embedded_Sided)
    dserror("no coupling interface available, however Embedded_Sided coupling is activated!");

  // characteristic element length to be computed
  double h_k = 0.0;

  // measure of the face (surface in 3D, line in 2D) or measure of the cut-face
  double meas_surf = 0.0;

  // measure of the element volume or measure of the physical cut part of the element volume
  double meas_vol = 0.0;

  switch (visc_stab_hk)
  {
    //---------------------------------------------------
    // volume-equivalent diameter
    //---------------------------------------------------
    case INPAR::XFEM::ViscStab_hk_vol_equivalent:
    {
      // evaluate shape functions and derivatives at element center
      if (averaging_strategy == INPAR::XFEM::Embedded_Sided)
      {
        // evaluate shape functions and derivatives at element center w.r.t embedded element
        meas_vol = emb->EvalElementVolume();
      }
      else
      {
        // meas_vol =
        // XFEM::UTILS::EvalElementVolume<distype>(my::xyze_,&(my::weights_),&(my::myknots_));
        meas_vol = XFEM::UTILS::EvalElementVolume<distype>(xyze);
      }

      // compute h_k as volume-equivalent diameter and directly return the value
      return h_k = XFEM::UTILS::ComputeVolEqDiameter(meas_vol);

      break;
    }
    //---------------------------------------------------
    // compute h_k as physical/cut volume divided by physical partial/cut surface measure
    // ( used to estimate the cut-dependent inverse estimate on cut elements, not useful for sliver
    // and/or dotted cut situations)
    //---------------------------------------------------
    case INPAR::XFEM::ViscStab_hk_cut_vol_div_by_cut_surf:
    {
      if (averaging_strategy == INPAR::XFEM::Embedded_Sided)
        dserror("ViscStab_hk_cut_vol_div_by_cut_surf not reasonable for Embedded_Sided_Coupling!");

      // compute the cut surface measure
      meas_surf = XFEM::UTILS::ComputeMeasCutSurf(bintpoints, bcells);

      if (fabs(meas_surf) < 1.e-8) dserror("Element contribution to interface has zero size.");

      // compute the cut volume measure
      for (GEO::CUT::plain_volumecell_set::const_iterator i = vcSet.begin(); i != vcSet.end(); ++i)
      {
        GEO::CUT::VolumeCell *vc = *i;
        meas_vol += vc->Volume();
      }

      if (meas_vol < 0.0)
        dserror(
            " measure of cut partial volume is smaller than 0.0: %f Attention with increasing "
            "Nitsche-Parameter!!!",
            meas_vol);

      break;
    }
    //---------------------------------------------------
    // full element volume divided by physical partial/cut surface measure ( used to estimate the
    // cut-dependent inverse estimate on cut elements, however, avoids problems with sliver cuts,
    // not useful for dotted cuts)
    //---------------------------------------------------
    case INPAR::XFEM::ViscStab_hk_ele_vol_div_by_cut_surf:
    {
      if (averaging_strategy == INPAR::XFEM::Embedded_Sided)
        dserror("ViscStab_hk_ele_vol_div_by_cut_surf not reasonable for Embedded_Sided_Coupling!");

      // compute the cut surface measure
      meas_surf = XFEM::UTILS::ComputeMeasCutSurf(bintpoints, bcells);

      // compute the full element volume measure
      meas_vol = XFEM::UTILS::EvalElementVolume<distype>(xyze);

      break;
    }
    //---------------------------------------------------
    // full element volume divided by surface measure ( used for uncut situations, standard weak
    // Dirichlet boundary/coupling conditions)
    //---------------------------------------------------
    case INPAR::XFEM::ViscStab_hk_ele_vol_div_by_ele_surf:
    {
      if (averaging_strategy != INPAR::XFEM::Embedded_Sided)
        dserror("ViscStab_hk_ele_vol_div_by_ele_surf just reasonable for Embedded_Sided_Coupling!");

      DRT::FaceElement *fele = dynamic_cast<DRT::FaceElement *>(face);
      if (!fele) dserror("Cast to FaceElement failed!");

      //---------------------------------------------------
      // compute the uncut element's surface measure
      meas_surf = XFEM::UTILS::ComputeMeasFace(ele, ele_xyze, fele->FaceParentNumber(), nsd);

      // compute the full element volume measure
      meas_vol = emb->EvalElementVolume();

      break;
    }
    //---------------------------------------------------
    case INPAR::XFEM::ViscStab_hk_ele_vol_div_by_max_ele_surf:
      //---------------------------------------------------
      {
        if (averaging_strategy == INPAR::XFEM::Embedded_Sided)
          dserror(
              "ViscStab_hk_ele_vol_div_by_max_ele_surf not reasonable for "
              "Embedded_Sided_Coupling!");

        // compute the uncut element's surface measure
        const int numfaces = DRT::UTILS::getNumberOfElementFaces(ele->Shape());

        // loop all surfaces
        for (int lid = 0; lid < numfaces; ++lid)
        {
          meas_surf = std::max(meas_surf, XFEM::UTILS::ComputeMeasFace(ele, ele_xyze, lid, nsd));
        }

        // compute the full element volume measure
        meas_vol = XFEM::UTILS::EvalElementVolume<distype>(xyze);

        break;
      }
    default:
      dserror("unknown type of characteristic element length");
      break;
  }

  //--------------------------------------
  // compute the final element length if fraction-based computation and not returned yet
  h_k = meas_vol / meas_surf;
  //--------------------------------------

  // check plausibility
  if (h_k < 1e-14)
    dserror("the characteristic element length is zero or smaller, it has not been set properly!");

  return h_k;
}

/*--------------------------------------------------------------------------------
 *    compute stabilization factor for the Nitsche's penalty term
 *--------------------------------------------------------------------------------*/
void XFEM::UTILS::NIT_Compute_FullPenalty_Stabfac(
    double &NIT_full_stab_fac,             ///< to be filled: full Nitsche's penalty term scaling
                                           ///< (viscous+convective part)
    const LINALG::Matrix<3, 1> &normal,    ///< interface-normal vector
    const double h_k,                      ///< characteristic element length
    const double kappa_m,                  ///< Weight parameter (parameter +/master side)
    const double kappa_s,                  ///< Weight parameter (parameter -/slave  side)
    const LINALG::Matrix<3, 1> &velint_m,  ///< Master side velocity at gauss-point
    const LINALG::Matrix<3, 1> &velint_s,  ///< Slave side velocity at gauss-point
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
        XFF_ConvStabScaling,    ///< which convective stab. scaling on XFF interface
    const bool IsConservative,  ///< conservative formulation of navier stokes
    bool error_calc             ///< when called in error calculation, don't add the inflow terms
)
{
  // TEUCHOS_FUNC_TIME_MONITOR("XFEM::UTILS::::NIT_Compute_FullPenalty_Stabfac");

  //------------------------------------------------------------------------------
  // compute the full Nitsche parameter

  /*
   * Depending on the flow regime, the factor alpha of Nitsches penalty term
   * (\alpha * [v],[u]) can take various forms.
   * Based on INPAR::XFEM::MassConservationCombination, we choose:
   *
   *                       (1)           (2)          (3)
   *                    /  \mu    \rho             h * \rho         \
   *  NIT :=  \gamma * |    --  +  -- * |u|_inf  + ----------------- |
   *                    \   h      6               12 * \theta * dt /
   *
   *       OR:
   *                    /  \mu    \rho             h * \rho         \
   *  NIT :=  max      |    --  ;  -- * |u|_inf  ; ----------------- | *\gamma
   *                    \   h      6               12 * \theta * dt /
   *
   *          (1) NIT_visc_stab_fac = \gamma * \mu/h
   *
   *          (2) convective contribution
   *
   *          (3) transient contribution
   *
   * see Schott and Rasthofer, 'A face-oriented stabilized Nitsche-type extended variational
   * multiscale method for incompressible two-phase flow', Int. J. Numer. Meth. Engng, 2014
   *
   *
   * If INPAR::XFEM::MassConservationScaling_only_visc is set, we choose only (1),
   * no matter if the combination option max or sum is active!
   *
   */


  // (1)
  NIT_full_stab_fac = NIT_visc_stab_fac;

  if (MassConservationScaling == INPAR::XFEM::MassConservationScaling_full)
  {
    // TODO: Raffaela: which velocity has to be evaluated for these terms in ALE? the convective
    // velocity or the velint?
    double velnorminf_m = velint_m.NormInf();  // relative convective velocity
    double velnorminf_s = velint_s.NormInf();

    // take the maximum of viscous & convective contribution or the sum?
    if (MassConservationCombination == INPAR::XFEM::MassConservationCombination_max)
    {
      NIT_full_stab_fac =
          std::max(NIT_full_stab_fac, NITStabScaling *
                                          (kappa_m * densaf_master * fabs(velnorminf_m) +
                                              kappa_s * densaf_slave * fabs(velnorminf_s)) /
                                          6.0);
      if (!isstationary)
        NIT_full_stab_fac =
            std::max(NITStabScaling * h_k * (kappa_m * densaf_master + kappa_s * densaf_slave) /
                         (12.0 * timefac),
                NIT_full_stab_fac);
    }
    else  // the sum
    {
      // (2)
      NIT_full_stab_fac += NITStabScaling *
                           (kappa_m * densaf_master * fabs(velnorminf_m) +
                               kappa_s * densaf_slave * fabs(velnorminf_s)) /
                           6.0;  // THIS ONE NEEDS CHANGING!

      // (3)
      if (!isstationary)
        NIT_full_stab_fac += NITStabScaling * h_k *
                             (kappa_m * densaf_master + kappa_s * densaf_slave) / (12.0 * timefac);
    }
  }
  else if (MassConservationScaling != INPAR::XFEM::MassConservationScaling_only_visc)
    dserror("Unknown scaling choice in calculation of Nitsche's penalty parameter");

  if (IsConservative and (XFF_ConvStabScaling != INPAR::XFEM::XFF_ConvStabScaling_none or
                             ConvStabScaling != INPAR::XFEM::ConvStabScaling_none))
  {
    dserror(
        "convective stabilization is not available for conservative form of Navier-Stokes, but "
        "possible to implement!");
  }

  //----------------------------------------------------------------------------------------------
  // add inflow terms to ensure coercivity at inflow boundaries in the convective limit

  if ((ConvStabScaling == INPAR::XFEM::ConvStabScaling_none &&
          XFF_ConvStabScaling == INPAR::XFEM::XFF_ConvStabScaling_none) ||
      XFF_ConvStabScaling == INPAR::XFEM::XFF_ConvStabScaling_only_averaged || error_calc)
    return;

  const double veln_normal = velint_m.Dot(normal);

  double NIT_inflow_stab = 0.0;

  if (XFF_ConvStabScaling == INPAR::XFEM::XFF_ConvStabScaling_upwinding)
  {
    NIT_inflow_stab = fabs(veln_normal) * 0.5;
  }
  else
  {
    if (ConvStabScaling == INPAR::XFEM::ConvStabScaling_abs_inflow)
    {
      //      | u*n |
      NIT_inflow_stab = fabs(veln_normal);
    }
    else if (ConvStabScaling == INPAR::XFEM::ConvStabScaling_inflow)
    {
      //      ( -u*n ) if (u*n)<0 (inflow), conv_stabfac >= 0
      NIT_inflow_stab = std::max(0.0, -veln_normal);
    }
    else
      dserror("No valid INPAR::XFEM::ConvStabScaling for xfluid/xfsi problems");
  }

  NIT_inflow_stab *= densaf_master;  // my::densaf_;

  // Todo (kruse): it is planned to add the inflow contributions independent from the max. option!
  // This version is only kept to shift the adaption of test results to a single commit.
  if (MassConservationCombination == INPAR::XFEM::MassConservationCombination_max)
  {
    NIT_full_stab_fac = std::max(NIT_full_stab_fac, NIT_inflow_stab);
  }
  else if (MassConservationCombination == INPAR::XFEM::MassConservationCombination_sum)
  {
    NIT_full_stab_fac += NIT_inflow_stab;
  }
  else
    dserror("Unknown combination type in calculation of Nitsche's penalty parameter.");

  return;
}

template double XFEM::UTILS::EvalElementVolume<DRT::Element::hex8>(
    LINALG::Matrix<3, DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::hex8>::numNodePerElement>,
    LINALG::Matrix<DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::hex8>::numNodePerElement, 1> *,
    std::vector<Epetra_SerialDenseVector> *);
template double XFEM::UTILS::EvalElementVolume<DRT::Element::hex20>(
    LINALG::Matrix<3, DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::hex20>::numNodePerElement>,
    LINALG::Matrix<DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::hex20>::numNodePerElement, 1> *,
    std::vector<Epetra_SerialDenseVector> *);
template double XFEM::UTILS::EvalElementVolume<DRT::Element::hex27>(
    LINALG::Matrix<3, DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::hex27>::numNodePerElement>,
    LINALG::Matrix<DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::hex27>::numNodePerElement, 1> *,
    std::vector<Epetra_SerialDenseVector> *);
template double XFEM::UTILS::EvalElementVolume<DRT::Element::tet4>(
    LINALG::Matrix<3, DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::tet4>::numNodePerElement>,
    LINALG::Matrix<DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::tet4>::numNodePerElement, 1> *,
    std::vector<Epetra_SerialDenseVector> *);
template double XFEM::UTILS::EvalElementVolume<DRT::Element::tet10>(
    LINALG::Matrix<3, DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::tet10>::numNodePerElement>,
    LINALG::Matrix<DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::tet10>::numNodePerElement, 1> *,
    std::vector<Epetra_SerialDenseVector> *);
template double XFEM::UTILS::EvalElementVolume<DRT::Element::wedge6>(
    LINALG::Matrix<3, DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::wedge6>::numNodePerElement>,
    LINALG::Matrix<DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::wedge6>::numNodePerElement, 1>
        *,
    std::vector<Epetra_SerialDenseVector> *);
template double XFEM::UTILS::EvalElementVolume<DRT::Element::wedge15>(
    LINALG::Matrix<3, DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::wedge15>::numNodePerElement>,
    LINALG::Matrix<DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::wedge15>::numNodePerElement, 1>
        *,
    std::vector<Epetra_SerialDenseVector> *);

template double XFEM::UTILS::ComputeCharEleLength<DRT::Element::hex8>(DRT::Element *,
    Epetra_SerialDenseMatrix &, const Teuchos::RCP<XFEM::ConditionManager> &,
    const GEO::CUT::plain_volumecell_set &,
    const std::map<int, std::vector<GEO::CUT::BoundaryCell *>> &,
    const std::map<int, std::vector<DRT::UTILS::GaussIntegration>> &,
    const INPAR::XFEM::ViscStab_hk,
    Teuchos::RCP<DRT::ELEMENTS::XFLUID::SlaveElementInterface<DRT::Element::hex8>>, DRT::Element *);
template double XFEM::UTILS::ComputeCharEleLength<DRT::Element::hex20>(DRT::Element *,
    Epetra_SerialDenseMatrix &, const Teuchos::RCP<XFEM::ConditionManager> &,
    const GEO::CUT::plain_volumecell_set &,
    const std::map<int, std::vector<GEO::CUT::BoundaryCell *>> &,
    const std::map<int, std::vector<DRT::UTILS::GaussIntegration>> &,
    const INPAR::XFEM::ViscStab_hk,
    Teuchos::RCP<DRT::ELEMENTS::XFLUID::SlaveElementInterface<DRT::Element::hex20>>,
    DRT::Element *);
template double XFEM::UTILS::ComputeCharEleLength<DRT::Element::hex27>(DRT::Element *,
    Epetra_SerialDenseMatrix &, const Teuchos::RCP<XFEM::ConditionManager> &,
    const GEO::CUT::plain_volumecell_set &,
    const std::map<int, std::vector<GEO::CUT::BoundaryCell *>> &,
    const std::map<int, std::vector<DRT::UTILS::GaussIntegration>> &,
    const INPAR::XFEM::ViscStab_hk,
    Teuchos::RCP<DRT::ELEMENTS::XFLUID::SlaveElementInterface<DRT::Element::hex27>>,
    DRT::Element *);
template double XFEM::UTILS::ComputeCharEleLength<DRT::Element::tet4>(DRT::Element *,
    Epetra_SerialDenseMatrix &, const Teuchos::RCP<XFEM::ConditionManager> &,
    const GEO::CUT::plain_volumecell_set &,
    const std::map<int, std::vector<GEO::CUT::BoundaryCell *>> &,
    const std::map<int, std::vector<DRT::UTILS::GaussIntegration>> &,
    const INPAR::XFEM::ViscStab_hk,
    Teuchos::RCP<DRT::ELEMENTS::XFLUID::SlaveElementInterface<DRT::Element::tet4>>, DRT::Element *);
template double XFEM::UTILS::ComputeCharEleLength<DRT::Element::tet10>(DRT::Element *,
    Epetra_SerialDenseMatrix &, const Teuchos::RCP<XFEM::ConditionManager> &,
    const GEO::CUT::plain_volumecell_set &,
    const std::map<int, std::vector<GEO::CUT::BoundaryCell *>> &,
    const std::map<int, std::vector<DRT::UTILS::GaussIntegration>> &,
    const INPAR::XFEM::ViscStab_hk,
    Teuchos::RCP<DRT::ELEMENTS::XFLUID::SlaveElementInterface<DRT::Element::tet10>>,
    DRT::Element *);
template double XFEM::UTILS::ComputeCharEleLength<DRT::Element::wedge6>(DRT::Element *,
    Epetra_SerialDenseMatrix &, const Teuchos::RCP<XFEM::ConditionManager> &,
    const GEO::CUT::plain_volumecell_set &,
    const std::map<int, std::vector<GEO::CUT::BoundaryCell *>> &,
    const std::map<int, std::vector<DRT::UTILS::GaussIntegration>> &,
    const INPAR::XFEM::ViscStab_hk,
    Teuchos::RCP<DRT::ELEMENTS::XFLUID::SlaveElementInterface<DRT::Element::wedge6>>,
    DRT::Element *);
template double XFEM::UTILS::ComputeCharEleLength<DRT::Element::wedge15>(DRT::Element *,
    Epetra_SerialDenseMatrix &, const Teuchos::RCP<XFEM::ConditionManager> &,
    const GEO::CUT::plain_volumecell_set &,
    const std::map<int, std::vector<GEO::CUT::BoundaryCell *>> &,
    const std::map<int, std::vector<DRT::UTILS::GaussIntegration>> &,
    const INPAR::XFEM::ViscStab_hk,
    Teuchos::RCP<DRT::ELEMENTS::XFLUID::SlaveElementInterface<DRT::Element::wedge15>>,
    DRT::Element *);
