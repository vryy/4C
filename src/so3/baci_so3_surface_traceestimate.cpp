/*----------------------------------------------------------------------*/
/*! \file

\brief evaluation for nitsche trace inequality estimate


\level 3
*----------------------------------------------------------------------*/

#include "baci_discretization_fem_general_utils_boundary_integration.H"
#include "baci_discretization_fem_general_utils_fem_shapefunctions.H"
#include "baci_discretization_fem_general_utils_gausspoints.H"
#include "baci_discretization_fem_general_utils_nurbs_shapefunctions.H"
#include "baci_global_data.H"
#include "baci_lib_element_integration_select.H"
#include "baci_linalg_utils_densematrix_determinant.H"
#include "baci_linalg_utils_densematrix_eigen.H"
#include "baci_mat_fourieriso.H"
#include "baci_mat_service.H"
#include "baci_mat_so3_material.H"
#include "baci_nurbs_discret.H"
#include "baci_so3_surface.H"

BACI_NAMESPACE_OPEN


/*----------------------------------------------------------------------*
 |                                                           seitz 11/16|
 *----------------------------------------------------------------------*/
double DRT::ELEMENTS::StructuralSurface::EstimateNitscheTraceMaxEigenvalueCombined(
    std::vector<double>& parent_disp)
{
  switch (ParentElement()->Shape())
  {
    case CORE::FE::CellType::hex8:
      if (Shape() == CORE::FE::CellType::quad4)
        return EstimateNitscheTraceMaxEigenvalueCombined<CORE::FE::CellType::hex8,
            CORE::FE::CellType::quad4>(parent_disp);
      else
        dserror("how can an hex8 element have a surface that is not quad4 ???");
      break;
    case CORE::FE::CellType::hex27:
      return EstimateNitscheTraceMaxEigenvalueCombined<CORE::FE::CellType::hex27,
          CORE::FE::CellType::quad9>(parent_disp);
      break;
    case CORE::FE::CellType::tet4:
      return EstimateNitscheTraceMaxEigenvalueCombined<CORE::FE::CellType::tet4,
          CORE::FE::CellType::tri3>(parent_disp);
      break;
    case CORE::FE::CellType::nurbs27:
      return EstimateNitscheTraceMaxEigenvalueCombined<CORE::FE::CellType::nurbs27,
          CORE::FE::CellType::nurbs9>(parent_disp);
      break;
    default:
      dserror("parent shape not implemented");
  }

  return 0;
}


template <CORE::FE::CellType dt_vol, CORE::FE::CellType dt_surf>
double DRT::ELEMENTS::StructuralSurface::EstimateNitscheTraceMaxEigenvalueCombined(
    std::vector<double>& parent_disp)
{
  const int dim = CORE::FE::dim<dt_vol>;
  const int num_dof = CORE::FE::num_nodes<dt_vol> * CORE::FE::dim<dt_vol>;
  const int dim_image = CORE::FE::num_nodes<dt_vol> * CORE::FE::dim<dt_vol> -
                        CORE::FE::dim<dt_vol> * (CORE::FE::dim<dt_vol> + 1) / 2;

  CORE::LINALG::Matrix<CORE::FE::num_nodes<dt_vol>, 3> xrefe;
  CORE::LINALG::Matrix<CORE::FE::num_nodes<dt_vol>, 3> xcurr;

  for (int i = 0; i < ParentElement()->NumNode(); ++i)
    for (int d = 0; d < dim; ++d)
    {
      xrefe(i, d) = ParentElement()->Nodes()[i]->X()[d];
      xcurr(i, d) = xrefe(i, d) + parent_disp[i * dim + d];
    }

  CORE::LINALG::Matrix<num_dof, num_dof> vol, surf;

  TraceEstimateVolMatrix<dt_vol>(xrefe, xcurr, vol);
  TraceEstimateSurfMatrix<dt_vol, dt_surf>(xrefe, xcurr, surf);

  CORE::LINALG::Matrix<num_dof, dim_image> proj, tmp;
  SubspaceProjector<dt_vol>(xcurr, proj);

  CORE::LINALG::Matrix<dim_image, dim_image> vol_red, surf_red;

  tmp.Multiply(vol, proj);
  vol_red.MultiplyTN(proj, tmp);
  tmp.Multiply(surf, proj);
  surf_red.MultiplyTN(proj, tmp);

  CORE::LINALG::SerialDenseMatrix vol_red_sd(
      Teuchos::View, vol_red.A(), dim_image, dim_image, dim_image);
  CORE::LINALG::SerialDenseMatrix surf_red_sd(
      Teuchos::View, surf_red.A(), dim_image, dim_image, dim_image);

  return CORE::LINALG::GeneralizedEigen(surf_red_sd, vol_red_sd);
}

template <CORE::FE::CellType dt_vol>
void DRT::ELEMENTS::StructuralSurface::TraceEstimateVolMatrix(
    const CORE::LINALG::Matrix<CORE::FE::num_nodes<dt_vol>, 3>& xrefe,
    const CORE::LINALG::Matrix<CORE::FE::num_nodes<dt_vol>, 3>& xcurr,
    CORE::LINALG::Matrix<CORE::FE::num_nodes<dt_vol> * 3, CORE::FE::num_nodes<dt_vol> * 3>& vol)
{
  const int dim = CORE::FE::dim<dt_vol>;

  double jac;
  CORE::LINALG::Matrix<3, 3> defgrd;
  CORE::LINALG::Matrix<3, 3> rcg;
  CORE::LINALG::Matrix<6, 1> glstrain;
  CORE::LINALG::Matrix<6, CORE::FE::num_nodes<dt_vol> * 3> bop;
  CORE::LINALG::Matrix<CORE::FE::num_nodes<dt_vol> * 3, 6> bc;
  CORE::LINALG::Matrix<dim, CORE::FE::num_nodes<dt_vol>> N_XYZ;

  CORE::FE::IntPointsAndWeights<dim> ip(DRT::ELEMENTS::DisTypeToOptGaussRule<dt_vol>::rule);

  for (int gp = 0; gp < ip.IP().nquad; ++gp)
  {
    const CORE::LINALG::Matrix<3, 1> xi(ip.IP().qxg[gp], false);
    Strains<dt_vol>(xrefe, xcurr, xi, jac, defgrd, glstrain, rcg, bop, N_XYZ);

    CORE::LINALG::Matrix<6, 6> cmat(true);
    CORE::LINALG::Matrix<6, 1> stress(true);
    Teuchos::ParameterList params;
    Teuchos::rcp_dynamic_cast<MAT::So3Material>(ParentElement()->Material())
        ->Evaluate(&defgrd, &glstrain, params, &stress, &cmat, gp, ParentElement()->Id());
    bc.MultiplyTN(bop, cmat);
    vol.Multiply(ip.IP().qwgt[gp] * jac, bc, bop, 1.);
  }

  return;
}


template <CORE::FE::CellType dt_vol, CORE::FE::CellType dt_surf>
void DRT::ELEMENTS::StructuralSurface::TraceEstimateSurfMatrix(
    const CORE::LINALG::Matrix<CORE::FE::num_nodes<dt_vol>, 3>& xrefe,
    const CORE::LINALG::Matrix<CORE::FE::num_nodes<dt_vol>, 3>& xcurr,
    CORE::LINALG::Matrix<CORE::FE::num_nodes<dt_vol> * 3, CORE::FE::num_nodes<dt_vol> * 3>& surf)
{
  const int dim = CORE::FE::dim<dt_vol>;

  CORE::LINALG::Matrix<6, 6> id4;
  for (int i = 0; i < 3; ++i) id4(i, i) = 1.;
  for (int i = 3; i < 6; ++i) id4(i, i) = 2.;

  CORE::LINALG::SerialDenseMatrix xrefe_surf(CORE::FE::num_nodes<dt_surf>, dim);
  MaterialConfiguration(xrefe_surf);

  std::vector<double> n(3);
  CORE::LINALG::Matrix<3, 1> n_v(n.data(), true);
  CORE::LINALG::Matrix<3, 3> nn;
  double detA;
  double jac;
  CORE::LINALG::Matrix<3, 3> defgrd;
  CORE::LINALG::Matrix<3, 3> rcg;
  CORE::LINALG::Matrix<6, 1> glstrain;
  CORE::LINALG::Matrix<6, CORE::FE::num_nodes<dt_vol> * 3> bop;
  CORE::LINALG::Matrix<CORE::FE::num_nodes<dt_vol> * 3, 6> bc;
  CORE::LINALG::Matrix<dim, CORE::FE::num_nodes<dt_vol>> N_XYZ;

  CORE::FE::IntPointsAndWeights<dim - 1> ip(DRT::ELEMENTS::DisTypeToOptGaussRule<dt_surf>::rule);
  CORE::LINALG::SerialDenseMatrix deriv_surf(2, CORE::FE::num_nodes<dt_surf>);

  for (int gp = 0; gp < ip.IP().nquad; ++gp)
  {
    CORE::FE::CollectedGaussPoints intpoints =
        CORE::FE::CollectedGaussPoints(1);  // reserve just for 1 entry ...
    intpoints.Append(ip.IP().qxg[gp][0], ip.IP().qxg[gp][1], 0.0, ip.IP().qwgt[gp]);

    // get coordinates of gauss point w.r.t. local parent coordinate system
    CORE::LINALG::SerialDenseMatrix pqxg(1, 3);
    CORE::LINALG::Matrix<3, 3> derivtrafo;

    CORE::FE::BoundaryGPToParentGP<3>(
        pqxg, derivtrafo, intpoints, ParentElement()->Shape(), Shape(), FaceParentNumber());

    CORE::LINALG::Matrix<3, 1> xi;
    for (int i = 0; i < 3; ++i) xi(i) = pqxg(0, i);
    Strains<dt_vol>(xrefe, xcurr, xi, jac, defgrd, glstrain, rcg, bop, N_XYZ);

    CORE::LINALG::Matrix<6, 6> cmat(true);
    CORE::LINALG::Matrix<6, 1> stress(true);
    Teuchos::ParameterList params;
    Teuchos::rcp_dynamic_cast<MAT::So3Material>(ParentElement()->Material())
        ->Evaluate(&defgrd, &glstrain, params, &stress, &cmat, gp, ParentElement()->Id());

    double normalfac = 1.;
    if (Shape() == CORE::FE::CellType::nurbs9)
    {
      std::vector<CORE::LINALG::SerialDenseVector> parentknots(dim);
      std::vector<CORE::LINALG::SerialDenseVector> boundaryknots(dim - 1);
      dynamic_cast<DRT::NURBS::NurbsDiscretization*>(
          DRT::Problem::Instance()->GetDis("structure").get())
          ->GetKnotVector()
          ->GetBoundaryEleAndParentKnots(
              parentknots, boundaryknots, normalfac, ParentElement()->Id(), FaceParentNumber());

      CORE::LINALG::Matrix<CORE::FE::num_nodes<dt_surf>, 1> weights, shapefcn;
      for (int i = 0; i < CORE::FE::num_nodes<dt_surf>; ++i)
        weights(i) = dynamic_cast<DRT::NURBS::ControlPoint*>(Nodes()[i])->W();

      CORE::LINALG::Matrix<2, 1> xi_surf;
      xi_surf(0) = ip.IP().qxg[gp][0];
      xi_surf(1) = ip.IP().qxg[gp][1];
      CORE::FE::NURBS::nurbs_get_2D_funct_deriv(
          shapefcn, deriv_surf, xi_surf, boundaryknots, weights, dt_surf);
    }
    else
      CORE::FE::shape_function_2D_deriv1(
          deriv_surf, ip.IP().qxg[gp][0], ip.IP().qxg[gp][1], Shape());

    SurfaceIntegration(detA, n, xrefe_surf, deriv_surf);
    n_v.Scale(normalfac);
    n_v.Scale(1. / n_v.Norm2());
    nn.MultiplyNT(n_v, n_v);

    CORE::LINALG::Matrix<6, 6> cn;
    MAT::AddSymmetricHolzapfelProduct(cn, rcg, nn, .25);

    CORE::LINALG::Matrix<6, 6> tmp1, tmp2;
    tmp1.Multiply(cmat, id4);
    tmp2.Multiply(tmp1, cn);
    tmp1.Multiply(tmp2, id4);
    tmp2.Multiply(tmp1, cmat);

    CORE::LINALG::Matrix<CORE::FE::num_nodes<dt_vol> * 3, 6> tmp3;
    tmp3.MultiplyTN(bop, tmp2);

    surf.Multiply(detA * ip.IP().qwgt[gp], tmp3, bop, 1.);
  }

  return;
}

template <CORE::FE::CellType dt_vol>
void DRT::ELEMENTS::StructuralSurface::Strains(
    const CORE::LINALG::Matrix<CORE::FE::num_nodes<dt_vol>, 3>& xrefe,
    const CORE::LINALG::Matrix<CORE::FE::num_nodes<dt_vol>, 3>& xcurr,
    const CORE::LINALG::Matrix<3, 1>& xi, double& jac, CORE::LINALG::Matrix<3, 3>& defgrd,
    CORE::LINALG::Matrix<6, 1>& glstrain, CORE::LINALG::Matrix<3, 3>& rcg,
    CORE::LINALG::Matrix<6, CORE::FE::num_nodes<dt_vol> * 3>& bop,
    CORE::LINALG::Matrix<3, CORE::FE::num_nodes<dt_vol>>& N_XYZ)
{
  const int dim = CORE::FE::dim<dt_vol>;
  const int num_node = CORE::FE::num_nodes<dt_vol>;
  CORE::LINALG::Matrix<dim, num_node> deriv;

  if (dt_vol == CORE::FE::CellType::nurbs27)
  {
    std::vector<CORE::LINALG::SerialDenseVector> knots;
    dynamic_cast<DRT::NURBS::NurbsDiscretization*>(
        DRT::Problem::Instance()->GetDis("structure").get())
        ->GetKnotVector()
        ->GetEleKnots(knots, ParentElementId());

    CORE::LINALG::Matrix<CORE::FE::num_nodes<dt_vol>, 1> weights, shapefcn;

    for (int i = 0; i < CORE::FE::num_nodes<dt_vol>; ++i)
      weights(i) = dynamic_cast<DRT::NURBS::ControlPoint*>(ParentElement()->Nodes()[i])->W();

    CORE::FE::NURBS::nurbs_get_3D_funct_deriv(shapefcn, deriv, xi, knots, weights, dt_vol);
  }
  else
    CORE::FE::shape_function_deriv1<dt_vol>(xi, deriv);

  CORE::LINALG::Matrix<dim, dim> invJ;
  invJ.Multiply(deriv, xrefe);
  jac = invJ.Invert();
  N_XYZ.Multiply(invJ, deriv);
  defgrd.MultiplyTT(xcurr, N_XYZ);

  rcg.MultiplyTN(defgrd, defgrd);
  glstrain(0) = 0.5 * (rcg(0, 0) - 1.0);
  glstrain(1) = 0.5 * (rcg(1, 1) - 1.0);
  glstrain(2) = 0.5 * (rcg(2, 2) - 1.0);
  glstrain(3) = rcg(0, 1);
  glstrain(4) = rcg(1, 2);
  glstrain(5) = rcg(2, 0);

  for (int i = 0; i < num_node; ++i)
  {
    bop(0, dim * i + 0) = defgrd(0, 0) * N_XYZ(0, i);
    bop(0, dim * i + 1) = defgrd(1, 0) * N_XYZ(0, i);
    bop(0, dim * i + 2) = defgrd(2, 0) * N_XYZ(0, i);
    bop(1, dim * i + 0) = defgrd(0, 1) * N_XYZ(1, i);
    bop(1, dim * i + 1) = defgrd(1, 1) * N_XYZ(1, i);
    bop(1, dim * i + 2) = defgrd(2, 1) * N_XYZ(1, i);
    bop(2, dim * i + 0) = defgrd(0, 2) * N_XYZ(2, i);
    bop(2, dim * i + 1) = defgrd(1, 2) * N_XYZ(2, i);
    bop(2, dim * i + 2) = defgrd(2, 2) * N_XYZ(2, i);
    /* ~~~ */
    bop(3, dim * i + 0) = defgrd(0, 0) * N_XYZ(1, i) + defgrd(0, 1) * N_XYZ(0, i);
    bop(3, dim * i + 1) = defgrd(1, 0) * N_XYZ(1, i) + defgrd(1, 1) * N_XYZ(0, i);
    bop(3, dim * i + 2) = defgrd(2, 0) * N_XYZ(1, i) + defgrd(2, 1) * N_XYZ(0, i);
    bop(4, dim * i + 0) = defgrd(0, 1) * N_XYZ(2, i) + defgrd(0, 2) * N_XYZ(1, i);
    bop(4, dim * i + 1) = defgrd(1, 1) * N_XYZ(2, i) + defgrd(1, 2) * N_XYZ(1, i);
    bop(4, dim * i + 2) = defgrd(2, 1) * N_XYZ(2, i) + defgrd(2, 2) * N_XYZ(1, i);
    bop(5, dim * i + 0) = defgrd(0, 2) * N_XYZ(0, i) + defgrd(0, 0) * N_XYZ(2, i);
    bop(5, dim * i + 1) = defgrd(1, 2) * N_XYZ(0, i) + defgrd(1, 0) * N_XYZ(2, i);
    bop(5, dim * i + 2) = defgrd(2, 2) * N_XYZ(0, i) + defgrd(2, 0) * N_XYZ(2, i);
  }

  return;
}


template <CORE::FE::CellType dt_vol>
void DRT::ELEMENTS::StructuralSurface::SubspaceProjector(
    const CORE::LINALG::Matrix<CORE::FE::num_nodes<dt_vol>, 3>& xcurr,
    CORE::LINALG::Matrix<CORE::FE::num_nodes<dt_vol> * CORE::FE::dim<dt_vol>,
        CORE::FE::num_nodes<dt_vol> * CORE::FE::dim<dt_vol> -
            CORE::FE::dim<dt_vol>*(CORE::FE::dim<dt_vol> + 1) / 2>& proj)
{
  const int dim = CORE::FE::dim<dt_vol>;
  const int num_node = CORE::FE::num_nodes<dt_vol>;
  if (dim != 3) dserror("this should be 3D");

  CORE::LINALG::Matrix<3, 1> c;
  for (int r = 0; r < (int)xcurr.numRows(); ++r)
    for (int d = 0; d < (int)xcurr.numCols(); ++d) c(d) += xcurr(r, d);
  c.Scale(1. / xcurr.numRows());

  CORE::LINALG::Matrix<dim, 1> r[3];
  for (int i = 0; i < 3; ++i) r[i](i) = 1.;

  // basis, where the first six entries are the rigid body modes and the
  // remaining are constructed to be orthogonal to the rigid body modes
  CORE::LINALG::Matrix<dim * num_node, 1> basis[dim * num_node];

  // rigid body translations
  for (int i = 0; i < dim; ++i)
    for (int j = 0; j < num_node; ++j) basis[i](j * dim + i) = 1.;

  // rigid body rotations
  for (int i = 0; i < dim; ++i)
    for (int j = 0; j < num_node; ++j)
    {
      CORE::LINALG::Matrix<3, 1> x;
      for (int d = 0; d < 3; ++d) x(d) = xcurr(j, d);
      x.Update(-1., c, 1.);
      CORE::LINALG::Matrix<3, 1> cross;
      cross.CrossProduct(r[i], x);
      for (int k = 0; k < 3; ++k) basis[i + 3](j * 3 + k) = cross(k);
    }
  for (int i = 0; i < 6; ++i) basis[i].Scale(1. / basis[i].Norm2());

  // build the remaining basis vectors by generalized cross products
  for (int i = 6; i < dim * num_node; ++i)
  {
    double sign = +1.;
    int off = 0;
    bool new_basis_found = false;
    for (off = 0; (off < dim * num_node - i) && !new_basis_found; ++off)
    {
      for (int j = 0; j < i + 1; ++j)
      {
        CORE::LINALG::SerialDenseMatrix det(i, i, true);
        for (int c = 0; c < i; ++c)
        {
          for (int k = 0; k < j; ++k) det(k, c) = basis[c](k + off);
          for (int k = j; k < i; ++k) det(k, c) = basis[c](k + 1 + off);
        }
        basis[i](j + off) = CORE::LINALG::DeterminantLU(det) * sign;
        sign *= -1.;
      }
      if (basis[i].Norm2() > 1.e-6)
      {
        basis[i].Scale(1. / basis[i].Norm2());
        new_basis_found = true;
      }
    }
    if (!new_basis_found) dserror("no new basis vector found");
  }

  // at this point basis should already contain an ONB.
  // due to cut-off errors we do another sweep of Gram-Schmidt
  for (int i = 0; i < dim * num_node; ++i)
  {
    const CORE::LINALG::Matrix<dim * num_node, 1> tmp(basis[i]);
    for (int j = 0; j < i; ++j) basis[i].Update(-tmp.Dot(basis[j]), basis[j], 1.);

    basis[i].Scale(1. / basis[i].Norm2());
  }

  // hand out the projection matrix, i.e. the ONB not containing rigid body modes
  for (int i = 0; i < dim * num_node; ++i)
    for (int j = 6; j < dim * num_node; ++j) proj(i, j - 6) = basis[j](i);
}



/*----------------------------------------------------------------------*
 |                                                           seitz 11/16|
 *----------------------------------------------------------------------*/
double DRT::ELEMENTS::StructuralSurface::EstimateNitscheTraceMaxEigenvalueTSI(
    std::vector<double>& parent_disp)
{
  switch (ParentElement()->Shape())
  {
    case CORE::FE::CellType::hex8:
      if (Shape() == CORE::FE::CellType::quad4)
        return EstimateNitscheTraceMaxEigenvalueTSI<CORE::FE::CellType::hex8,
            CORE::FE::CellType::quad4>(parent_disp);
      else
        dserror("how can an hex8 element have a surface that is not quad4 ???");
      break;
    case CORE::FE::CellType::hex27:
      return EstimateNitscheTraceMaxEigenvalueTSI<CORE::FE::CellType::hex27,
          CORE::FE::CellType::quad9>(parent_disp);
    case CORE::FE::CellType::tet4:
      return EstimateNitscheTraceMaxEigenvalueTSI<CORE::FE::CellType::tet4,
          CORE::FE::CellType::tri3>(parent_disp);
    case CORE::FE::CellType::nurbs27:
      return EstimateNitscheTraceMaxEigenvalueTSI<CORE::FE::CellType::nurbs27,
          CORE::FE::CellType::nurbs9>(parent_disp);
    default:
      dserror("parent shape not implemented");
  }

  return 0;
}

template <CORE::FE::CellType dt_vol, CORE::FE::CellType dt_surf>
double DRT::ELEMENTS::StructuralSurface::EstimateNitscheTraceMaxEigenvalueTSI(
    std::vector<double>& parent_disp)
{
  const int dim = CORE::FE::dim<dt_vol>;
  const int num_dof = CORE::FE::num_nodes<dt_vol>;
  const int dim_image = CORE::FE::num_nodes<dt_vol> - 1;

  CORE::LINALG::Matrix<CORE::FE::num_nodes<dt_vol>, 3> xrefe;
  CORE::LINALG::Matrix<CORE::FE::num_nodes<dt_vol>, 3> xcurr;

  for (int i = 0; i < ParentElement()->NumNode(); ++i)
    for (int d = 0; d < dim; ++d)
    {
      xrefe(i, d) = ParentElement()->Nodes()[i]->X()[d];
      xcurr(i, d) = xrefe(i, d) + parent_disp[i * dim + d];
    }

  CORE::LINALG::Matrix<num_dof, num_dof> vol, surf;

  TraceEstimateVolMatrixTSI<dt_vol>(xrefe, xcurr, vol);
  TraceEstimateSurfMatrixTSI<dt_vol, dt_surf>(xrefe, xcurr, surf);


  CORE::LINALG::Matrix<num_dof, dim_image> proj, tmp;
  SubspaceProjectorScalar<dt_vol>(proj);

  CORE::LINALG::Matrix<dim_image, dim_image> vol_red, surf_red;

  tmp.Multiply(vol, proj);
  vol_red.MultiplyTN(proj, tmp);
  tmp.Multiply(surf, proj);
  surf_red.MultiplyTN(proj, tmp);

  CORE::LINALG::SerialDenseMatrix vol_red_sd(
      Teuchos::View, vol_red.A(), dim_image, dim_image, dim_image);
  CORE::LINALG::SerialDenseMatrix surf_red_sd(
      Teuchos::View, surf_red.A(), dim_image, dim_image, dim_image);

  return CORE::LINALG::GeneralizedEigen(surf_red_sd, vol_red_sd);
}

template <CORE::FE::CellType dt_vol>
void DRT::ELEMENTS::StructuralSurface::TraceEstimateVolMatrixTSI(
    const CORE::LINALG::Matrix<CORE::FE::num_nodes<dt_vol>, 3>& xrefe,
    const CORE::LINALG::Matrix<CORE::FE::num_nodes<dt_vol>, 3>& xcurr,
    CORE::LINALG::Matrix<CORE::FE::num_nodes<dt_vol>, CORE::FE::num_nodes<dt_vol>>& vol)
{
  const int dim = CORE::FE::dim<dt_vol>;
  const int num_node = CORE::FE::num_nodes<dt_vol>;

  double jac;
  CORE::LINALG::Matrix<3, 3> defgrd;
  CORE::LINALG::Matrix<3, 3> rcg;
  CORE::LINALG::Matrix<6, 1> glstrain;
  CORE::LINALG::Matrix<6, CORE::FE::num_nodes<dt_vol> * 3> bop;
  CORE::LINALG::Matrix<CORE::FE::num_nodes<dt_vol> * 3, 6> bc;
  CORE::LINALG::Matrix<dim, num_node> N_XYZ, iC_N_XYZ;

  CORE::FE::IntPointsAndWeights<dim> ip(DRT::ELEMENTS::DisTypeToOptGaussRule<dt_vol>::rule);

  if (ParentElement()->NumMaterial() < 2) dserror("where's my second material");
  Teuchos::RCP<MAT::FourierIso> mat_thr =
      Teuchos::rcp_dynamic_cast<MAT::FourierIso>(ParentElement()->Material(1), true);
  const double k0 = mat_thr->Conductivity();

  for (int gp = 0; gp < ip.IP().nquad; ++gp)
  {
    const CORE::LINALG::Matrix<3, 1> xi(ip.IP().qxg[gp], false);
    Strains<dt_vol>(xrefe, xcurr, xi, jac, defgrd, glstrain, rcg, bop, N_XYZ);

    CORE::LINALG::Matrix<3, 3> iC;
    iC.MultiplyTN(defgrd, defgrd);
    iC.Invert();

    iC_N_XYZ.Multiply(iC, N_XYZ);
    iC_N_XYZ.Scale(k0);

    vol.MultiplyTN(ip.IP().qwgt[gp] * jac, N_XYZ, iC_N_XYZ, 1.);
  }
}


template <CORE::FE::CellType dt_vol, CORE::FE::CellType dt_surf>
void DRT::ELEMENTS::StructuralSurface::TraceEstimateSurfMatrixTSI(
    const CORE::LINALG::Matrix<CORE::FE::num_nodes<dt_vol>, 3>& xrefe,
    const CORE::LINALG::Matrix<CORE::FE::num_nodes<dt_vol>, 3>& xcurr,
    CORE::LINALG::Matrix<CORE::FE::num_nodes<dt_vol>, CORE::FE::num_nodes<dt_vol>>& surf)
{
  const int dim = CORE::FE::dim<dt_vol>;
  const int num_node = CORE::FE::num_nodes<dt_vol>;

  double jac;
  CORE::LINALG::Matrix<3, 3> defgrd;
  CORE::LINALG::Matrix<3, 3> rcg;
  CORE::LINALG::Matrix<6, 1> glstrain;
  CORE::LINALG::Matrix<6, CORE::FE::num_nodes<dt_vol> * 3> bop;
  CORE::LINALG::Matrix<CORE::FE::num_nodes<dt_vol> * 3, 6> bc;
  CORE::LINALG::Matrix<dim, num_node> N_XYZ;
  CORE::LINALG::Matrix<1, num_node> iCn_N_XYZ;

  CORE::LINALG::SerialDenseMatrix xrefe_surf(CORE::FE::num_nodes<dt_surf>, dim);
  MaterialConfiguration(xrefe_surf);

  std::vector<double> n(3);
  CORE::LINALG::Matrix<3, 1> n_v(n.data(), true), iCn;
  double detA;

  CORE::FE::IntPointsAndWeights<dim - 1> ip(DRT::ELEMENTS::DisTypeToOptGaussRule<dt_surf>::rule);
  CORE::LINALG::SerialDenseMatrix deriv_surf(2, CORE::FE::num_nodes<dt_surf>);

  if (ParentElement()->NumMaterial() < 2) dserror("where's my second material");
  Teuchos::RCP<MAT::FourierIso> mat_thr =
      Teuchos::rcp_dynamic_cast<MAT::FourierIso>(ParentElement()->Material(1), true);
  const double k0 = mat_thr->Conductivity();

  for (int gp = 0; gp < ip.IP().nquad; ++gp)
  {
    CORE::FE::shape_function_2D_deriv1(deriv_surf, ip.IP().qxg[gp][0], ip.IP().qxg[gp][1], Shape());
    SurfaceIntegration(detA, n, xrefe_surf, deriv_surf);
    n_v.Scale(1. / n_v.Norm2());

    CORE::FE::CollectedGaussPoints intpoints =
        CORE::FE::CollectedGaussPoints(1);  // reserve just for 1 entry ...
    intpoints.Append(ip.IP().qxg[gp][0], ip.IP().qxg[gp][1], 0.0, ip.IP().qwgt[gp]);

    // get coordinates of gauss point w.r.t. local parent coordinate system
    CORE::LINALG::SerialDenseMatrix pqxg(1, 3);
    CORE::LINALG::Matrix<3, 3> derivtrafo;

    CORE::FE::BoundaryGPToParentGP<3>(
        pqxg, derivtrafo, intpoints, ParentElement()->Shape(), Shape(), FaceParentNumber());

    CORE::LINALG::Matrix<3, 1> xi;
    for (int i = 0; i < 3; ++i) xi(i) = pqxg(0, i);

    Strains<dt_vol>(xrefe, xcurr, xi, jac, defgrd, glstrain, rcg, bop, N_XYZ);

    CORE::LINALG::Matrix<3, 3> iC;
    iC.MultiplyTN(defgrd, defgrd);
    iC.Invert();
    iCn.Multiply(iC, n_v);

    iCn_N_XYZ.MultiplyTN(iCn, N_XYZ);
    iCn_N_XYZ.Scale(k0);

    surf.MultiplyTN(detA * ip.IP().qwgt[gp], iCn_N_XYZ, iCn_N_XYZ, 1.);
  }
}



template <CORE::FE::CellType dt_vol>
void DRT::ELEMENTS::StructuralSurface::SubspaceProjectorScalar(
    CORE::LINALG::Matrix<CORE::FE::num_nodes<dt_vol>, CORE::FE::num_nodes<dt_vol> - 1>& proj)
{
  const int num_node = CORE::FE::num_nodes<dt_vol>;
  CORE::LINALG::Matrix<num_node, 1> basis[num_node];

  for (int i = 0; i < num_node; ++i) basis[0](i) = 1.;

  for (int i = 1; i < num_node; ++i)
  {
    double sign = +1.;
    for (int j = 0; j < i + 1; ++j)
    {
      CORE::LINALG::SerialDenseMatrix det(i, i, true);
      for (int c = 0; c < i; ++c)
      {
        for (int k = 0; k < j; ++k) det(k, c) = basis[c](k);
        for (int k = j; k < i; ++k) det(k, c) = basis[c](k + 1);
      }
      basis[i](j) = CORE::LINALG::DeterminantLU(det) * sign;
      sign *= -1.;
    }
    basis[i].Scale(1. / basis[i].Norm2());
  }

  // hand out the projection matrix, i.e. the ONB not containing rigid body modes
  for (int i = 0; i < num_node; ++i)
    for (int j = 1; j < num_node; ++j) proj(i, j - 1) = basis[j](i);
}

BACI_NAMESPACE_CLOSE
