/*----------------------------------------------------------------------*/
/*! \file

\brief little helpers for solid elements

\level 1

*----------------------------------------------------------------------*/

#ifndef SOLID_UTILS_H
#define SOLID_UTILS_H

#include "utils_integration.H"
#include "element.H"
#include "linalg_utils_densematrix_eigen.H"

namespace STR::UTILS
{
  template <DRT::UTILS::GaussRule3D rule>
  struct IntRuleToNquad
  {
  };
  template <>
  struct IntRuleToNquad<DRT::UTILS::GaussRule3D::hex_8point>
  {
    static const int ngp = 8;
  };
  template <>
  struct IntRuleToNquad<DRT::UTILS::GaussRule3D::hex_27point>
  {
    static const int ngp = 27;
  };
  template <>
  struct IntRuleToNquad<DRT::UTILS::GaussRule3D::hex_18point>
  {
    static const int ngp = 18;
  };
  template <>
  struct IntRuleToNquad<DRT::UTILS::GaussRule3D::pyramid_1point>
  {
    static const int ngp = 1;
  };
  template <>
  struct IntRuleToNquad<DRT::UTILS::GaussRule3D::pyramid_8point>
  {
    static const int ngp = 8;
  };
  template <>
  struct IntRuleToNquad<DRT::UTILS::GaussRule3D::tet_1point>
  {
    static const int ngp = 1;
  };
  template <>
  struct IntRuleToNquad<DRT::UTILS::GaussRule3D::tet_4point>
  {
    static const int ngp = 4;
  };
  template <>
  struct IntRuleToNquad<DRT::UTILS::GaussRule3D::tet_11point>
  {
    static const int ngp = 11;
  };


  int DisTypeToNgpOptGaussRule(DRT::Element::DiscretizationType distype);

  void Pk2ToCauchy(const LINALG::Matrix<6, 1>& pk2, const LINALG::Matrix<3, 3>& defgrd,
      LINALG::Matrix<6, 1>& cauchy);

  LINALG::Matrix<6, 1> GreenLagrangeToEulerAlmansi(
      const LINALG::Matrix<6, 1>& gl, const LINALG::Matrix<3, 3>& defgrd);

  template <unsigned dim>
  void CalcConsistentDefgrd(const LINALG::Matrix<dim, dim>& defgrd_disp,
      const LINALG::Matrix<dim*(dim + 1) / 2, 1>& glstrain_mod,
      LINALG::Matrix<dim, dim>& defgrd_mod)
  {
    LINALG::Matrix<dim, dim> R;       // rotation tensor
    LINALG::Matrix<dim, dim> U_mod;   // modified right stretch tensor
    LINALG::Matrix<dim, dim> U_disp;  // displacement-based right stretch tensor
    LINALG::Matrix<dim, dim> EW;      // temporarily store eigenvalues
    LINALG::Matrix<dim, dim> tmp;     // temporary matrix for matrix matrix matrix products
    LINALG::Matrix<dim, dim> tmp2;    // temporary matrix for matrix matrix matrix products

    // ******************************************************************
    // calculate modified right stretch tensor
    // ******************************************************************
    if (dim != 3) dserror("stop: this currently only works for 3D");
    for (unsigned i = 0; i < dim; i++) U_mod(i, i) = 2. * glstrain_mod(i) + 1.;
    U_mod(0, 1) = glstrain_mod(dim);
    U_mod(1, 0) = glstrain_mod(dim);
    U_mod(1, 2) = glstrain_mod(4);
    U_mod(2, 1) = glstrain_mod(4);
    U_mod(0, 2) = glstrain_mod(5);
    U_mod(2, 0) = glstrain_mod(5);

    LINALG::SYEV(U_mod, EW, U_mod);
    for (unsigned i = 0; i < dim; ++i) EW(i, i) = sqrt(EW(i, i));
    tmp.Multiply(U_mod, EW);
    tmp2.MultiplyNT(tmp, U_mod);
    U_mod.Update(tmp2);

    // ******************************************************************
    // calculate displacement-based right stretch tensor
    // ******************************************************************
    U_disp.MultiplyTN(defgrd_disp, defgrd_disp);

    LINALG::SYEV(U_disp, EW, U_disp);
    for (unsigned i = 0; i < dim; ++i) EW(i, i) = sqrt(EW(i, i));
    tmp.Multiply(U_disp, EW);
    tmp2.MultiplyNT(tmp, U_disp);
    U_disp.Update(tmp2);

    // ******************************************************************
    // compose consistent deformation gradient
    // ******************************************************************
    U_disp.Invert();
    R.Multiply(defgrd_disp, U_disp);
    defgrd_mod.Multiply(R, U_mod);
  }

}  // namespace STR::UTILS

#endif  // SOLID_UTILS_H
