/*----------------------------------------------------------------------*/
/*! \file

\brief evaluation for nitsche trace inequality estimate


\level 3
*----------------------------------------------------------------------*/

#include "4C_fem_general_element_integration_select.hpp"
#include "4C_fem_general_utils_boundary_integration.hpp"
#include "4C_fem_general_utils_fem_shapefunctions.hpp"
#include "4C_fem_general_utils_gausspoints.hpp"
#include "4C_fem_general_utils_nurbs_shapefunctions.hpp"
#include "4C_fem_nurbs_discretization.hpp"
#include "4C_global_data.hpp"
#include "4C_linalg_utils_densematrix_determinant.hpp"
#include "4C_linalg_utils_densematrix_eigen.hpp"
#include "4C_mat_fourieriso.hpp"
#include "4C_mat_service.hpp"
#include "4C_mat_so3_material.hpp"
#include "4C_so3_surface.hpp"

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*
 |                                                           seitz 11/16|
 *----------------------------------------------------------------------*/
double Discret::ELEMENTS::StructuralSurface::estimate_nitsche_trace_max_eigenvalue_combined(
    std::vector<double>& parent_disp)
{
  switch (parent_element()->Shape())
  {
    case Core::FE::CellType::hex8:
      if (Shape() == Core::FE::CellType::quad4)
        return estimate_nitsche_trace_max_eigenvalue_combined<Core::FE::CellType::hex8,
            Core::FE::CellType::quad4>(parent_disp);
      else
        FOUR_C_THROW("how can an hex8 element have a surface that is not quad4 ???");
      break;
    case Core::FE::CellType::hex27:
      return estimate_nitsche_trace_max_eigenvalue_combined<Core::FE::CellType::hex27,
          Core::FE::CellType::quad9>(parent_disp);
      break;
    case Core::FE::CellType::tet4:
      return estimate_nitsche_trace_max_eigenvalue_combined<Core::FE::CellType::tet4,
          Core::FE::CellType::tri3>(parent_disp);
      break;
    case Core::FE::CellType::nurbs27:
      return estimate_nitsche_trace_max_eigenvalue_combined<Core::FE::CellType::nurbs27,
          Core::FE::CellType::nurbs9>(parent_disp);
      break;
    default:
      FOUR_C_THROW("parent shape not implemented");
  }

  return 0;
}


template <Core::FE::CellType dt_vol, Core::FE::CellType dt_surf>
double Discret::ELEMENTS::StructuralSurface::estimate_nitsche_trace_max_eigenvalue_combined(
    std::vector<double>& parent_disp)
{
  const int dim = Core::FE::dim<dt_vol>;
  const int num_dof = Core::FE::num_nodes<dt_vol> * Core::FE::dim<dt_vol>;
  const int dim_image = Core::FE::num_nodes<dt_vol> * Core::FE::dim<dt_vol> -
                        Core::FE::dim<dt_vol> * (Core::FE::dim<dt_vol> + 1) / 2;

  Core::LinAlg::Matrix<Core::FE::num_nodes<dt_vol>, 3> xrefe;
  Core::LinAlg::Matrix<Core::FE::num_nodes<dt_vol>, 3> xcurr;

  for (int i = 0; i < parent_element()->num_node(); ++i)
    for (int d = 0; d < dim; ++d)
    {
      xrefe(i, d) = parent_element()->Nodes()[i]->X()[d];
      xcurr(i, d) = xrefe(i, d) + parent_disp[i * dim + d];
    }

  Core::LinAlg::Matrix<num_dof, num_dof> vol, surf;

  trace_estimate_vol_matrix<dt_vol>(xrefe, xcurr, vol);
  trace_estimate_surf_matrix<dt_vol, dt_surf>(xrefe, xcurr, surf);

  Core::LinAlg::Matrix<num_dof, dim_image> proj, tmp;
  subspace_projector<dt_vol>(xcurr, proj);

  Core::LinAlg::Matrix<dim_image, dim_image> vol_red, surf_red;

  tmp.Multiply(vol, proj);
  vol_red.MultiplyTN(proj, tmp);
  tmp.Multiply(surf, proj);
  surf_red.MultiplyTN(proj, tmp);

  Core::LinAlg::SerialDenseMatrix vol_red_sd(
      Teuchos::View, vol_red.A(), dim_image, dim_image, dim_image);
  Core::LinAlg::SerialDenseMatrix surf_red_sd(
      Teuchos::View, surf_red.A(), dim_image, dim_image, dim_image);

  return Core::LinAlg::GeneralizedEigen(surf_red_sd, vol_red_sd);
}

template <Core::FE::CellType dt_vol>
void Discret::ELEMENTS::StructuralSurface::trace_estimate_vol_matrix(
    const Core::LinAlg::Matrix<Core::FE::num_nodes<dt_vol>, 3>& xrefe,
    const Core::LinAlg::Matrix<Core::FE::num_nodes<dt_vol>, 3>& xcurr,
    Core::LinAlg::Matrix<Core::FE::num_nodes<dt_vol> * 3, Core::FE::num_nodes<dt_vol> * 3>& vol)
{
  const int dim = Core::FE::dim<dt_vol>;

  double jac;
  Core::LinAlg::Matrix<3, 3> defgrd;
  Core::LinAlg::Matrix<3, 3> rcg;
  Core::LinAlg::Matrix<6, 1> glstrain;
  Core::LinAlg::Matrix<6, Core::FE::num_nodes<dt_vol> * 3> bop;
  Core::LinAlg::Matrix<Core::FE::num_nodes<dt_vol> * 3, 6> bc;
  Core::LinAlg::Matrix<dim, Core::FE::num_nodes<dt_vol>> N_XYZ;

  Core::FE::IntPointsAndWeights<dim> ip(Discret::ELEMENTS::DisTypeToOptGaussRule<dt_vol>::rule);

  for (int gp = 0; gp < ip.IP().nquad; ++gp)
  {
    const Core::LinAlg::Matrix<3, 1> xi(ip.IP().qxg[gp], false);
    strains<dt_vol>(xrefe, xcurr, xi, jac, defgrd, glstrain, rcg, bop, N_XYZ);

    Core::LinAlg::Matrix<6, 6> cmat(true);
    Core::LinAlg::Matrix<6, 1> stress(true);
    Teuchos::ParameterList params;
    Teuchos::rcp_dynamic_cast<Mat::So3Material>(parent_element()->Material())
        ->Evaluate(&defgrd, &glstrain, params, &stress, &cmat, gp, parent_element()->Id());
    bc.MultiplyTN(bop, cmat);
    vol.Multiply(ip.IP().qwgt[gp] * jac, bc, bop, 1.);
  }

  return;
}


template <Core::FE::CellType dt_vol, Core::FE::CellType dt_surf>
void Discret::ELEMENTS::StructuralSurface::trace_estimate_surf_matrix(
    const Core::LinAlg::Matrix<Core::FE::num_nodes<dt_vol>, 3>& xrefe,
    const Core::LinAlg::Matrix<Core::FE::num_nodes<dt_vol>, 3>& xcurr,
    Core::LinAlg::Matrix<Core::FE::num_nodes<dt_vol> * 3, Core::FE::num_nodes<dt_vol> * 3>& surf)
{
  const int dim = Core::FE::dim<dt_vol>;

  Core::LinAlg::Matrix<6, 6> id4;
  for (int i = 0; i < 3; ++i) id4(i, i) = 1.;
  for (int i = 3; i < 6; ++i) id4(i, i) = 2.;

  Core::LinAlg::SerialDenseMatrix xrefe_surf(Core::FE::num_nodes<dt_surf>, dim);
  material_configuration(xrefe_surf);

  std::vector<double> n(3);
  Core::LinAlg::Matrix<3, 1> n_v(n.data(), true);
  Core::LinAlg::Matrix<3, 3> nn;
  double detA;
  double jac;
  Core::LinAlg::Matrix<3, 3> defgrd;
  Core::LinAlg::Matrix<3, 3> rcg;
  Core::LinAlg::Matrix<6, 1> glstrain;
  Core::LinAlg::Matrix<6, Core::FE::num_nodes<dt_vol> * 3> bop;
  Core::LinAlg::Matrix<Core::FE::num_nodes<dt_vol> * 3, 6> bc;
  Core::LinAlg::Matrix<dim, Core::FE::num_nodes<dt_vol>> N_XYZ;

  Core::FE::IntPointsAndWeights<dim - 1> ip(
      Discret::ELEMENTS::DisTypeToOptGaussRule<dt_surf>::rule);
  Core::LinAlg::SerialDenseMatrix deriv_surf(2, Core::FE::num_nodes<dt_surf>);

  for (int gp = 0; gp < ip.IP().nquad; ++gp)
  {
    Core::FE::CollectedGaussPoints intpoints =
        Core::FE::CollectedGaussPoints(1);  // reserve just for 1 entry ...
    intpoints.Append(ip.IP().qxg[gp][0], ip.IP().qxg[gp][1], 0.0, ip.IP().qwgt[gp]);

    // get coordinates of gauss point w.r.t. local parent coordinate system
    Core::LinAlg::SerialDenseMatrix pqxg(1, 3);
    Core::LinAlg::Matrix<3, 3> derivtrafo;

    Core::FE::BoundaryGPToParentGP<3>(
        pqxg, derivtrafo, intpoints, parent_element()->Shape(), Shape(), FaceParentNumber());

    Core::LinAlg::Matrix<3, 1> xi;
    for (int i = 0; i < 3; ++i) xi(i) = pqxg(0, i);
    strains<dt_vol>(xrefe, xcurr, xi, jac, defgrd, glstrain, rcg, bop, N_XYZ);

    Core::LinAlg::Matrix<6, 6> cmat(true);
    Core::LinAlg::Matrix<6, 1> stress(true);
    Teuchos::ParameterList params;
    Teuchos::rcp_dynamic_cast<Mat::So3Material>(parent_element()->Material())
        ->Evaluate(&defgrd, &glstrain, params, &stress, &cmat, gp, parent_element()->Id());

    double normalfac = 1.;
    if (Shape() == Core::FE::CellType::nurbs9)
    {
      std::vector<Core::LinAlg::SerialDenseVector> parentknots(dim);
      std::vector<Core::LinAlg::SerialDenseVector> boundaryknots(dim - 1);
      dynamic_cast<Discret::Nurbs::NurbsDiscretization*>(
          Global::Problem::Instance()->GetDis("structure").get())
          ->GetKnotVector()
          ->get_boundary_ele_and_parent_knots(
              parentknots, boundaryknots, normalfac, parent_element()->Id(), FaceParentNumber());

      Core::LinAlg::Matrix<Core::FE::num_nodes<dt_surf>, 1> weights, shapefcn;
      for (int i = 0; i < Core::FE::num_nodes<dt_surf>; ++i)
        weights(i) = dynamic_cast<Discret::Nurbs::ControlPoint*>(Nodes()[i])->W();

      Core::LinAlg::Matrix<2, 1> xi_surf;
      xi_surf(0) = ip.IP().qxg[gp][0];
      xi_surf(1) = ip.IP().qxg[gp][1];
      Core::FE::Nurbs::nurbs_get_2D_funct_deriv(
          shapefcn, deriv_surf, xi_surf, boundaryknots, weights, dt_surf);
    }
    else
      Core::FE::shape_function_2D_deriv1(
          deriv_surf, ip.IP().qxg[gp][0], ip.IP().qxg[gp][1], Shape());

    surface_integration(detA, n, xrefe_surf, deriv_surf);
    n_v.Scale(normalfac);
    n_v.Scale(1. / n_v.Norm2());
    nn.MultiplyNT(n_v, n_v);

    Core::LinAlg::Matrix<6, 6> cn;
    Mat::AddSymmetricHolzapfelProduct(cn, rcg, nn, .25);

    Core::LinAlg::Matrix<6, 6> tmp1, tmp2;
    tmp1.Multiply(cmat, id4);
    tmp2.Multiply(tmp1, cn);
    tmp1.Multiply(tmp2, id4);
    tmp2.Multiply(tmp1, cmat);

    Core::LinAlg::Matrix<Core::FE::num_nodes<dt_vol> * 3, 6> tmp3;
    tmp3.MultiplyTN(bop, tmp2);

    surf.Multiply(detA * ip.IP().qwgt[gp], tmp3, bop, 1.);
  }

  return;
}

template <Core::FE::CellType dt_vol>
void Discret::ELEMENTS::StructuralSurface::strains(
    const Core::LinAlg::Matrix<Core::FE::num_nodes<dt_vol>, 3>& xrefe,
    const Core::LinAlg::Matrix<Core::FE::num_nodes<dt_vol>, 3>& xcurr,
    const Core::LinAlg::Matrix<3, 1>& xi, double& jac, Core::LinAlg::Matrix<3, 3>& defgrd,
    Core::LinAlg::Matrix<6, 1>& glstrain, Core::LinAlg::Matrix<3, 3>& rcg,
    Core::LinAlg::Matrix<6, Core::FE::num_nodes<dt_vol> * 3>& bop,
    Core::LinAlg::Matrix<3, Core::FE::num_nodes<dt_vol>>& N_XYZ)
{
  const int dim = Core::FE::dim<dt_vol>;
  const int num_node = Core::FE::num_nodes<dt_vol>;
  Core::LinAlg::Matrix<dim, num_node> deriv;

  if (dt_vol == Core::FE::CellType::nurbs27)
  {
    std::vector<Core::LinAlg::SerialDenseVector> knots;
    dynamic_cast<Discret::Nurbs::NurbsDiscretization*>(
        Global::Problem::Instance()->GetDis("structure").get())
        ->GetKnotVector()
        ->GetEleKnots(knots, ParentElementId());

    Core::LinAlg::Matrix<Core::FE::num_nodes<dt_vol>, 1> weights, shapefcn;

    for (int i = 0; i < Core::FE::num_nodes<dt_vol>; ++i)
      weights(i) = dynamic_cast<Discret::Nurbs::ControlPoint*>(parent_element()->Nodes()[i])->W();

    Core::FE::Nurbs::nurbs_get_3D_funct_deriv(shapefcn, deriv, xi, knots, weights, dt_vol);
  }
  else
    Core::FE::shape_function_deriv1<dt_vol>(xi, deriv);

  Core::LinAlg::Matrix<dim, dim> invJ;
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


template <Core::FE::CellType dt_vol>
void Discret::ELEMENTS::StructuralSurface::subspace_projector(
    const Core::LinAlg::Matrix<Core::FE::num_nodes<dt_vol>, 3>& xcurr,
    Core::LinAlg::Matrix<Core::FE::num_nodes<dt_vol> * Core::FE::dim<dt_vol>,
        Core::FE::num_nodes<dt_vol> * Core::FE::dim<dt_vol> -
            Core::FE::dim<dt_vol>*(Core::FE::dim<dt_vol> + 1) / 2>& proj)
{
  const int dim = Core::FE::dim<dt_vol>;
  const int num_node = Core::FE::num_nodes<dt_vol>;
  if (dim != 3) FOUR_C_THROW("this should be 3D");

  Core::LinAlg::Matrix<3, 1> c;
  for (int r = 0; r < (int)xcurr.numRows(); ++r)
    for (int d = 0; d < (int)xcurr.numCols(); ++d) c(d) += xcurr(r, d);
  c.Scale(1. / xcurr.numRows());

  Core::LinAlg::Matrix<dim, 1> r[3];
  for (int i = 0; i < 3; ++i) r[i](i) = 1.;

  // basis, where the first six entries are the rigid body modes and the
  // remaining are constructed to be orthogonal to the rigid body modes
  Core::LinAlg::Matrix<dim * num_node, 1> basis[dim * num_node];

  // rigid body translations
  for (int i = 0; i < dim; ++i)
    for (int j = 0; j < num_node; ++j) basis[i](j * dim + i) = 1.;

  // rigid body rotations
  for (int i = 0; i < dim; ++i)
    for (int j = 0; j < num_node; ++j)
    {
      Core::LinAlg::Matrix<3, 1> x;
      for (int d = 0; d < 3; ++d) x(d) = xcurr(j, d);
      x.Update(-1., c, 1.);
      Core::LinAlg::Matrix<3, 1> cross;
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
        Core::LinAlg::SerialDenseMatrix det(i, i, true);
        for (int c = 0; c < i; ++c)
        {
          for (int k = 0; k < j; ++k) det(k, c) = basis[c](k + off);
          for (int k = j; k < i; ++k) det(k, c) = basis[c](k + 1 + off);
        }
        basis[i](j + off) = Core::LinAlg::DeterminantLU(det) * sign;
        sign *= -1.;
      }
      if (basis[i].Norm2() > 1.e-6)
      {
        basis[i].Scale(1. / basis[i].Norm2());
        new_basis_found = true;
      }
    }
    if (!new_basis_found) FOUR_C_THROW("no new basis vector found");
  }

  // at this point basis should already contain an ONB.
  // due to cut-off errors we do another sweep of Gram-Schmidt
  for (int i = 0; i < dim * num_node; ++i)
  {
    const Core::LinAlg::Matrix<dim * num_node, 1> tmp(basis[i]);
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
double Discret::ELEMENTS::StructuralSurface::estimate_nitsche_trace_max_eigenvalue_tsi(
    std::vector<double>& parent_disp)
{
  switch (parent_element()->Shape())
  {
    case Core::FE::CellType::hex8:
      if (Shape() == Core::FE::CellType::quad4)
        return estimate_nitsche_trace_max_eigenvalue_tsi<Core::FE::CellType::hex8,
            Core::FE::CellType::quad4>(parent_disp);
      else
        FOUR_C_THROW("how can an hex8 element have a surface that is not quad4 ???");
      break;
    case Core::FE::CellType::hex27:
      return estimate_nitsche_trace_max_eigenvalue_tsi<Core::FE::CellType::hex27,
          Core::FE::CellType::quad9>(parent_disp);
    case Core::FE::CellType::tet4:
      return estimate_nitsche_trace_max_eigenvalue_tsi<Core::FE::CellType::tet4,
          Core::FE::CellType::tri3>(parent_disp);
    case Core::FE::CellType::nurbs27:
      return estimate_nitsche_trace_max_eigenvalue_tsi<Core::FE::CellType::nurbs27,
          Core::FE::CellType::nurbs9>(parent_disp);
    default:
      FOUR_C_THROW("parent shape not implemented");
  }

  return 0;
}

template <Core::FE::CellType dt_vol, Core::FE::CellType dt_surf>
double Discret::ELEMENTS::StructuralSurface::estimate_nitsche_trace_max_eigenvalue_tsi(
    std::vector<double>& parent_disp)
{
  const int dim = Core::FE::dim<dt_vol>;
  const int num_dof = Core::FE::num_nodes<dt_vol>;
  const int dim_image = Core::FE::num_nodes<dt_vol> - 1;

  Core::LinAlg::Matrix<Core::FE::num_nodes<dt_vol>, 3> xrefe;
  Core::LinAlg::Matrix<Core::FE::num_nodes<dt_vol>, 3> xcurr;

  for (int i = 0; i < parent_element()->num_node(); ++i)
    for (int d = 0; d < dim; ++d)
    {
      xrefe(i, d) = parent_element()->Nodes()[i]->X()[d];
      xcurr(i, d) = xrefe(i, d) + parent_disp[i * dim + d];
    }

  Core::LinAlg::Matrix<num_dof, num_dof> vol, surf;

  trace_estimate_vol_matrix_tsi<dt_vol>(xrefe, xcurr, vol);
  trace_estimate_surf_matrix_tsi<dt_vol, dt_surf>(xrefe, xcurr, surf);


  Core::LinAlg::Matrix<num_dof, dim_image> proj, tmp;
  subspace_projector_scalar<dt_vol>(proj);

  Core::LinAlg::Matrix<dim_image, dim_image> vol_red, surf_red;

  tmp.Multiply(vol, proj);
  vol_red.MultiplyTN(proj, tmp);
  tmp.Multiply(surf, proj);
  surf_red.MultiplyTN(proj, tmp);

  Core::LinAlg::SerialDenseMatrix vol_red_sd(
      Teuchos::View, vol_red.A(), dim_image, dim_image, dim_image);
  Core::LinAlg::SerialDenseMatrix surf_red_sd(
      Teuchos::View, surf_red.A(), dim_image, dim_image, dim_image);

  return Core::LinAlg::GeneralizedEigen(surf_red_sd, vol_red_sd);
}

template <Core::FE::CellType dt_vol>
void Discret::ELEMENTS::StructuralSurface::trace_estimate_vol_matrix_tsi(
    const Core::LinAlg::Matrix<Core::FE::num_nodes<dt_vol>, 3>& xrefe,
    const Core::LinAlg::Matrix<Core::FE::num_nodes<dt_vol>, 3>& xcurr,
    Core::LinAlg::Matrix<Core::FE::num_nodes<dt_vol>, Core::FE::num_nodes<dt_vol>>& vol)
{
  const int dim = Core::FE::dim<dt_vol>;
  const int num_node = Core::FE::num_nodes<dt_vol>;

  double jac;
  Core::LinAlg::Matrix<3, 3> defgrd;
  Core::LinAlg::Matrix<3, 3> rcg;
  Core::LinAlg::Matrix<6, 1> glstrain;
  Core::LinAlg::Matrix<6, Core::FE::num_nodes<dt_vol> * 3> bop;
  Core::LinAlg::Matrix<Core::FE::num_nodes<dt_vol> * 3, 6> bc;
  Core::LinAlg::Matrix<dim, num_node> N_XYZ, iC_N_XYZ;

  Core::FE::IntPointsAndWeights<dim> ip(Discret::ELEMENTS::DisTypeToOptGaussRule<dt_vol>::rule);

  if (parent_element()->NumMaterial() < 2) FOUR_C_THROW("where's my second material");
  Teuchos::RCP<Mat::FourierIso> mat_thr =
      Teuchos::rcp_dynamic_cast<Mat::FourierIso>(parent_element()->Material(1), true);
  const double k0 = mat_thr->Conductivity();

  for (int gp = 0; gp < ip.IP().nquad; ++gp)
  {
    const Core::LinAlg::Matrix<3, 1> xi(ip.IP().qxg[gp], false);
    strains<dt_vol>(xrefe, xcurr, xi, jac, defgrd, glstrain, rcg, bop, N_XYZ);

    Core::LinAlg::Matrix<3, 3> iC;
    iC.MultiplyTN(defgrd, defgrd);
    iC.Invert();

    iC_N_XYZ.Multiply(iC, N_XYZ);
    iC_N_XYZ.Scale(k0);

    vol.MultiplyTN(ip.IP().qwgt[gp] * jac, N_XYZ, iC_N_XYZ, 1.);
  }
}


template <Core::FE::CellType dt_vol, Core::FE::CellType dt_surf>
void Discret::ELEMENTS::StructuralSurface::trace_estimate_surf_matrix_tsi(
    const Core::LinAlg::Matrix<Core::FE::num_nodes<dt_vol>, 3>& xrefe,
    const Core::LinAlg::Matrix<Core::FE::num_nodes<dt_vol>, 3>& xcurr,
    Core::LinAlg::Matrix<Core::FE::num_nodes<dt_vol>, Core::FE::num_nodes<dt_vol>>& surf)
{
  const int dim = Core::FE::dim<dt_vol>;
  const int num_node = Core::FE::num_nodes<dt_vol>;

  double jac;
  Core::LinAlg::Matrix<3, 3> defgrd;
  Core::LinAlg::Matrix<3, 3> rcg;
  Core::LinAlg::Matrix<6, 1> glstrain;
  Core::LinAlg::Matrix<6, Core::FE::num_nodes<dt_vol> * 3> bop;
  Core::LinAlg::Matrix<Core::FE::num_nodes<dt_vol> * 3, 6> bc;
  Core::LinAlg::Matrix<dim, num_node> N_XYZ;
  Core::LinAlg::Matrix<1, num_node> iCn_N_XYZ;

  Core::LinAlg::SerialDenseMatrix xrefe_surf(Core::FE::num_nodes<dt_surf>, dim);
  material_configuration(xrefe_surf);

  std::vector<double> n(3);
  Core::LinAlg::Matrix<3, 1> n_v(n.data(), true), iCn;
  double detA;

  Core::FE::IntPointsAndWeights<dim - 1> ip(
      Discret::ELEMENTS::DisTypeToOptGaussRule<dt_surf>::rule);
  Core::LinAlg::SerialDenseMatrix deriv_surf(2, Core::FE::num_nodes<dt_surf>);

  if (parent_element()->NumMaterial() < 2) FOUR_C_THROW("where's my second material");
  Teuchos::RCP<Mat::FourierIso> mat_thr =
      Teuchos::rcp_dynamic_cast<Mat::FourierIso>(parent_element()->Material(1), true);
  const double k0 = mat_thr->Conductivity();

  for (int gp = 0; gp < ip.IP().nquad; ++gp)
  {
    Core::FE::shape_function_2D_deriv1(deriv_surf, ip.IP().qxg[gp][0], ip.IP().qxg[gp][1], Shape());
    surface_integration(detA, n, xrefe_surf, deriv_surf);
    n_v.Scale(1. / n_v.Norm2());

    Core::FE::CollectedGaussPoints intpoints =
        Core::FE::CollectedGaussPoints(1);  // reserve just for 1 entry ...
    intpoints.Append(ip.IP().qxg[gp][0], ip.IP().qxg[gp][1], 0.0, ip.IP().qwgt[gp]);

    // get coordinates of gauss point w.r.t. local parent coordinate system
    Core::LinAlg::SerialDenseMatrix pqxg(1, 3);
    Core::LinAlg::Matrix<3, 3> derivtrafo;

    Core::FE::BoundaryGPToParentGP<3>(
        pqxg, derivtrafo, intpoints, parent_element()->Shape(), Shape(), FaceParentNumber());

    Core::LinAlg::Matrix<3, 1> xi;
    for (int i = 0; i < 3; ++i) xi(i) = pqxg(0, i);

    strains<dt_vol>(xrefe, xcurr, xi, jac, defgrd, glstrain, rcg, bop, N_XYZ);

    Core::LinAlg::Matrix<3, 3> iC;
    iC.MultiplyTN(defgrd, defgrd);
    iC.Invert();
    iCn.Multiply(iC, n_v);

    iCn_N_XYZ.MultiplyTN(iCn, N_XYZ);
    iCn_N_XYZ.Scale(k0);

    surf.MultiplyTN(detA * ip.IP().qwgt[gp], iCn_N_XYZ, iCn_N_XYZ, 1.);
  }
}



template <Core::FE::CellType dt_vol>
void Discret::ELEMENTS::StructuralSurface::subspace_projector_scalar(
    Core::LinAlg::Matrix<Core::FE::num_nodes<dt_vol>, Core::FE::num_nodes<dt_vol> - 1>& proj)
{
  const int num_node = Core::FE::num_nodes<dt_vol>;
  Core::LinAlg::Matrix<num_node, 1> basis[num_node];

  for (int i = 0; i < num_node; ++i) basis[0](i) = 1.;

  for (int i = 1; i < num_node; ++i)
  {
    double sign = +1.;
    for (int j = 0; j < i + 1; ++j)
    {
      Core::LinAlg::SerialDenseMatrix det(i, i, true);
      for (int c = 0; c < i; ++c)
      {
        for (int k = 0; k < j; ++k) det(k, c) = basis[c](k);
        for (int k = j; k < i; ++k) det(k, c) = basis[c](k + 1);
      }
      basis[i](j) = Core::LinAlg::DeterminantLU(det) * sign;
      sign *= -1.;
    }
    basis[i].Scale(1. / basis[i].Norm2());
  }

  // hand out the projection matrix, i.e. the ONB not containing rigid body modes
  for (int i = 0; i < num_node; ++i)
    for (int j = 1; j < num_node; ++j) proj(i, j - 1) = basis[j](i);
}

FOUR_C_NAMESPACE_CLOSE
