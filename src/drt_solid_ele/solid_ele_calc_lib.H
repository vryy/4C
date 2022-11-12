/*----------------------------------------------------------------------*/
/*! \file

\brief A library of free functions for a default solid element
\level 1

*----------------------------------------------------------------------*/

#ifndef SOLID_ELE_CALC_LIB_H
#define SOLID_ELE_CALC_LIB_H

#include <Epetra_Vector.h>
#include <Teuchos_ParameterList.hpp>
#include <numeric>
#include "drt_element.H"
#include "drt_utils_local_connectivity_matrices.H"
#include "solid_ele.H"
#include "drt_utils.H"
#include "so3_material.H"
#include "voigt_notation.H"


namespace
{
  template <DRT::Element::DiscretizationType distype>
  inline static constexpr int nen = DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;

  template <DRT::Element::DiscretizationType distype>
  inline static constexpr int nsd = DRT::UTILS::DisTypeToDim<distype>::dim;

  template <DRT::Element::DiscretizationType distype>
  inline static constexpr int numstr = nsd<distype>*(nsd<distype> + 1) / 2;

  template <DRT::Element::DiscretizationType distype>
  inline static constexpr int numdofperelement = nen<distype>* nsd<distype>;
}  // namespace

namespace DRT::ELEMENTS
{
  template <DRT::Element::DiscretizationType distype>
  struct NodalCoordinates
  {
    LINALG::Matrix<nen<distype>, nsd<distype>> reference_;
    LINALG::Matrix<nen<distype>, nsd<distype>> current_;
  };

  template <DRT::Element::DiscretizationType distype>
  NodalCoordinates<distype> EvaluateNodalCoordinates(const DRT::ELEMENTS::Solid& ele,
      const DRT::Discretization& discretization, const std::vector<int>& lm)
  {
    const Epetra_Vector& displacements = *discretization.GetState("displacement");
    std::vector<double> mydisp(lm.size());
    DRT::UTILS::ExtractMyValues(displacements, mydisp, lm);

    NodalCoordinates<distype> coordinates;
    for (int i = 0; i < nen<distype>; ++i)
    {
      for (int d = 0; d < nsd<distype>; ++d)
      {
        coordinates.reference_(i, d) = ele.Nodes()[i]->X()[d];
        coordinates.current_(i, d) = coordinates.reference_(i, d) + mydisp[i * nsd<distype> + d];
      }
    }

    return coordinates;
  }

  /*!
   * @brief Evaluates the parameter coordinate of the Gauss point according the the Gauss rule
   *
   * @tparam distype : Discretization type
   * @param intpoints (in) : Gauss integration points
   * @param gp (in) : id of the Gauss point
   * @return LINALG::Matrix<nsd<distype>, 1> : Coordinates of the Gauss Point in the parameter space
   */
  template <DRT::Element::DiscretizationType distype>
  LINALG::Matrix<nsd<distype>, 1> EvaluateParameterCoordinate(
      const DRT::UTILS::GaussIntegration& intpoints, const int gp)
  {
    LINALG::Matrix<nsd<distype>, 1> xi;
    for (int d = 0; d < nsd<distype>; ++d) xi(d) = intpoints.Point(gp)[d];

    return xi;
  }

  template <DRT::Element::DiscretizationType distype>
  struct ShapeFunctionsAndDerivatives
  {
    LINALG::Matrix<nen<distype>, 1> shapefunctions_;
    LINALG::Matrix<nsd<distype>, nen<distype>> derivatives_;
  };

  /*!
   * @brief Evaluates the shape functions and their derivatives at the specified point in the
   * parameter space
   *
   * @tparam distype : Discretizationt type
   * @param xi (in) : Coordinate in the parameter space
   * @return ShapeFunctionsAndDerivatives<distype> : An object holding the shape functions and the
   * first derivatives evaluated at the respective point in the parameter space
   */
  template <DRT::Element::DiscretizationType distype>
  ShapeFunctionsAndDerivatives<distype> EvaluateShapeFunctionsAndDerivs(
      const LINALG::Matrix<nsd<distype>, 1>& xi)
  {
    ShapeFunctionsAndDerivatives<distype> shapefcns;
    DRT::UTILS::shape_function<distype>(xi, shapefcns.shapefunctions_);
    DRT::UTILS::shape_function_deriv1<distype>(xi, shapefcns.derivatives_);

    return shapefcns;
  }

  template <DRT::Element::DiscretizationType distype>
  struct JacobianMapping
  {
    double determinant_;
    LINALG::Matrix<nsd<distype>, nsd<distype>> inverse_jacobian_;
    LINALG::Matrix<nsd<distype>, nen<distype>> n_xyz_;
    double integration_factor_;
  };

  /*!
   * @brief Evaluates the jacobian mapping of the element
   *
   * @tparam distype : Discretization type
   * @param shapefcns (in) : Shape functions and derivatives evaluated at the respective point in
   * the parameter space
   * @param nodal_coordniates (in) : Reference and current coordinates of the nodes of the element
   * @param intpoints (in) : Gauss integration rule
   * @param gp (in) : Id of the Gauss point
   * @return JacobianMapping<distype> : An object holding quantities of the jacobian mapping
   * (inverse Jacobian, determinant, derivatives of the shape functions w.r.t. XYZ, integration
   * factor)
   */
  template <DRT::Element::DiscretizationType distype>
  JacobianMapping<distype> EvaluateJacobianMapping(
      const ShapeFunctionsAndDerivatives<distype>& shapefcns,
      const NodalCoordinates<distype>& nodal_coordniates,
      const DRT::UTILS::GaussIntegration& intpoints, const int gp)
  {
    JacobianMapping<distype> jacobian;

    jacobian.inverse_jacobian_.Multiply(shapefcns.derivatives_, nodal_coordniates.reference_);
    jacobian.determinant_ = jacobian.inverse_jacobian_.Invert();
    jacobian.integration_factor_ = jacobian.determinant_ * intpoints.Weight(gp);
    jacobian.n_xyz_.Multiply(jacobian.inverse_jacobian_, shapefcns.derivatives_);

    return jacobian;
  }

  template <DRT::Element::DiscretizationType distype>
  double EvaluateJacobianDeterminant(const ShapeFunctionsAndDerivatives<distype>& shapefcns,
      const NodalCoordinates<distype>& nodal_coordniates)
  {
    LINALG::Matrix<nsd<distype>, nsd<distype>> jacobian;
    jacobian.Multiply(shapefcns.derivatives_, nodal_coordniates.reference_);

    return jacobian.Determinant();
  }

  template <DRT::Element::DiscretizationType distype>
  struct Strains
  {
    LINALG::Matrix<nsd<distype>, nsd<distype>> defgrd_;
    LINALG::Matrix<nsd<distype>, nsd<distype>> rcg_;
    LINALG::Matrix<numstr<distype>, 1> gl_strain_;
  };

  /*!
   * @brief Evaluates strain measures
   *
   * @tparam distype
   * @param nodal_coordinates (in) : Reference and current coordinates of the nodes of the element
   * @param jacobian_mapping (in) : Quantities of the jacobian mapping
   * @return Strains<distype> (out) : Strain measures of the element (deformation gradient, right
   * Cauchy-Green deformation tensor, Green-Lagrange strain tensor)
   */
  template <DRT::Element::DiscretizationType distype, std::enable_if_t<nsd<distype> == 3, int> = 0>
  Strains<distype> EvaluateStrains(const NodalCoordinates<distype>& nodal_coordinates,
      const JacobianMapping<distype>& jacobian_mapping)
  {
    Strains<distype> strains;

    strains.defgrd_.MultiplyTT(nodal_coordinates.current_, jacobian_mapping.n_xyz_);
    strains.rcg_.MultiplyTN(strains.defgrd_, strains.defgrd_);

    // GL strain vector glstrain={E11,E22,E33,2*E12,2*E23,2*E31}
    strains.gl_strain_(0) = 0.5 * (strains.rcg_(0, 0) - 1.0);
    strains.gl_strain_(1) = 0.5 * (strains.rcg_(1, 1) - 1.0);
    strains.gl_strain_(2) = 0.5 * (strains.rcg_(2, 2) - 1.0);
    strains.gl_strain_(3) = strains.rcg_(0, 1);
    strains.gl_strain_(4) = strains.rcg_(1, 2);
    strains.gl_strain_(5) = strains.rcg_(2, 0);

    return strains;
  }


  /*!
   * @brief Evaluates the strain gradient (B-Operator) of the specified element
   *
   * @tparam distype
   * @param jacobian_mapping (in) : Quantities of the jacobian mapping
   * @param strains (in) : Strain measures of the element
   * @return LINALG::Matrix<numstr<distype>, nsd<distype> * nen<distype>> : B-Operator
   */
  template <DRT::Element::DiscretizationType distype, std::enable_if_t<nsd<distype> == 3, int> = 0>
  LINALG::Matrix<numstr<distype>, nsd<distype> * nen<distype>> EvaluateStrainGradient(
      const JacobianMapping<distype>& jacobian_mapping, const Strains<distype>& strains)
  {
    // B-operator
    LINALG::Matrix<numstr<distype>, nsd<distype> * nen<distype>> Bop;
    for (int i = 0; i < nen<distype>; ++i)
    {
      for (int d = 0; d < nsd<distype>; ++d)
        for (int e = 0; e < nsd<distype>; ++e)
          Bop(d, nsd<distype> * i + e) = strains.defgrd_(e, d) * jacobian_mapping.n_xyz_(d, i);

      // TODO: make this work for 2D
      Bop(3, nsd<distype> * i + 0) = strains.defgrd_(0, 0) * jacobian_mapping.n_xyz_(1, i) +
                                     strains.defgrd_(0, 1) * jacobian_mapping.n_xyz_(0, i);
      Bop(3, nsd<distype> * i + 1) = strains.defgrd_(1, 0) * jacobian_mapping.n_xyz_(1, i) +
                                     strains.defgrd_(1, 1) * jacobian_mapping.n_xyz_(0, i);
      Bop(3, nsd<distype> * i + 2) = strains.defgrd_(2, 0) * jacobian_mapping.n_xyz_(1, i) +
                                     strains.defgrd_(2, 1) * jacobian_mapping.n_xyz_(0, i);
      Bop(4, nsd<distype> * i + 0) = strains.defgrd_(0, 1) * jacobian_mapping.n_xyz_(2, i) +
                                     strains.defgrd_(0, 2) * jacobian_mapping.n_xyz_(1, i);
      Bop(4, nsd<distype> * i + 1) = strains.defgrd_(1, 1) * jacobian_mapping.n_xyz_(2, i) +
                                     strains.defgrd_(1, 2) * jacobian_mapping.n_xyz_(1, i);
      Bop(4, nsd<distype> * i + 2) = strains.defgrd_(2, 1) * jacobian_mapping.n_xyz_(2, i) +
                                     strains.defgrd_(2, 2) * jacobian_mapping.n_xyz_(1, i);
      Bop(5, nsd<distype> * i + 0) = strains.defgrd_(0, 2) * jacobian_mapping.n_xyz_(0, i) +
                                     strains.defgrd_(0, 0) * jacobian_mapping.n_xyz_(2, i);
      Bop(5, nsd<distype> * i + 1) = strains.defgrd_(1, 2) * jacobian_mapping.n_xyz_(0, i) +
                                     strains.defgrd_(1, 0) * jacobian_mapping.n_xyz_(2, i);
      Bop(5, nsd<distype> * i + 2) = strains.defgrd_(2, 2) * jacobian_mapping.n_xyz_(0, i) +
                                     strains.defgrd_(2, 0) * jacobian_mapping.n_xyz_(2, i);
    }

    return Bop;
  }

  template <DRT::Element::DiscretizationType distype>
  struct Stress
  {
    LINALG::Matrix<numstr<distype>, 1> pk2_;
    LINALG::Matrix<numstr<distype>, numstr<distype>> cmat_;
  };

  /*!
   * @brief Evaluates the material stress (2. Piola Kirchhoff stress tensor and the linearization
   * w.r.t Green-Lagrange strain)
   *
   * @tparam distype
   * @param material (in) : Reference to the material
   * @param strains (in) : Strain measures of the element
   * @param params (in) : List of additional parameter to pass quantities from the time integrator
   * to the material
   * @param gp (in) : Gauss point
   * @param eleGID (in) : Global element id
   * @return Stress<distype> : Object holding the 2. Piola Kirchhoff stress tensor and the
   * linearization w.r.t. Green Lagrange strain tensor
   */
  template <DRT::Element::DiscretizationType distype>
  Stress<distype> EvaluateMaterialStress(MAT::So3Material& material,
      const Strains<distype>& strains, Teuchos::ParameterList& params, int gp, int eleGID)
  {
    Stress<distype> stress;

    material.Evaluate(
        &strains.defgrd_, &strains.gl_strain_, params, &stress.pk2_, &stress.cmat_, gp, eleGID);
    return stress;
  }

  /*!
   * @brief Adds the internal force vector contribution of one Gauss point
   *
   * @tparam distype
   * @param Bop (in) : Strain gradient (B-Operator)
   * @param stress (in) : Stress measures
   * @param integration_fac (in) : Integration factor (Gauss point weight times the determinant of
   * the jacobian)
   * @param force_vector (in/out) : Force vector where the local contribution is added to
   */
  template <DRT::Element::DiscretizationType distype>
  void AddInternalForceVector(
      const LINALG::Matrix<numstr<distype>, nsd<distype> * nen<distype>>& Bop,
      const Stress<distype> stress, const double integration_fac,
      LINALG::Matrix<nsd<distype> * nen<distype>, 1>& force_vector)
  {
    force_vector.MultiplyTN(integration_fac, Bop, stress.pk2_, 1.);
  }

  /*!
   * @brief Add elastic stiffness matrix contribution of one Gauss point
   *
   * @tparam distype
   * @param Bop (in) : Strain gradient (B-Operator)
   * @param stress (in) : Stress measures
   * @param integration_fac (in) : Integration factor (Gauss point weight times the determinant of
   * the jacobian)
   * @param stiffness_matrix (in/out) : stiffness matrix where the local contribution is added to
   */
  template <DRT::Element::DiscretizationType distype>
  void AddElasticStiffnessMatrix(
      const LINALG::Matrix<numstr<distype>, nsd<distype> * nen<distype>>& Bop,
      const Stress<distype> stress, const double integration_fac,
      LINALG::Matrix<nsd<distype> * nen<distype>, nsd<distype> * nen<distype>>& stiffness_matrix)
  {
    LINALG::Matrix<numstr<distype>, nen<distype> * nsd<distype>> cb;
    cb.Multiply(stress.cmat_, Bop);
    stiffness_matrix.MultiplyTN(integration_fac, Bop, cb, 1.0);
  }

  /*!
   * @brief Add geometric stiffness matrix contribution of one Gauss point
   *
   * @tparam distype
   * @param jacobian_mapping (in) : Jacobian mapping at the respective Gauss point
   * @param stress (in) : Stress measures
   * @param stiffness_matrix (in/out) : stiffness matrix where the local contribution is added to
   */
  template <DRT::Element::DiscretizationType distype>
  void AddGeometricStiffnessMatrix(const JacobianMapping<distype>& jacobian_mapping,
      const Stress<distype> stress,
      LINALG::Matrix<nsd<distype> * nen<distype>, nsd<distype> * nen<distype>>& stiffness_matrix)
  {
    std::array<double, 3> SmB_L;  // intermediate Sm.B_L
    // kgeo += (B_L^T . sigma . B_L) * detJ * w(gp)  with B_L = Ni,Xj see NiliFEM-Skript
    for (int inod = 0; inod < nen<distype>; ++inod)
    {
      SmB_L[0] = stress.pk2_(0) * jacobian_mapping.n_xyz_(0, inod) +
                 stress.pk2_(3) * jacobian_mapping.n_xyz_(1, inod) +
                 stress.pk2_(5) * jacobian_mapping.n_xyz_(2, inod);
      SmB_L[1] = stress.pk2_(3) * jacobian_mapping.n_xyz_(0, inod) +
                 stress.pk2_(1) * jacobian_mapping.n_xyz_(1, inod) +
                 stress.pk2_(4) * jacobian_mapping.n_xyz_(2, inod);
      SmB_L[2] = stress.pk2_(5) * jacobian_mapping.n_xyz_(0, inod) +
                 stress.pk2_(4) * jacobian_mapping.n_xyz_(1, inod) +
                 stress.pk2_(2) * jacobian_mapping.n_xyz_(2, inod);

      for (int jnod = 0; jnod < nen<distype>; ++jnod)
      {
        double bopstrbop = 0.0;  // intermediate value
        for (int idim = 0; idim < nsd<distype>; ++idim)
          bopstrbop += jacobian_mapping.n_xyz_(idim, jnod) * SmB_L[idim];

        for (int d = 0; d < nsd<distype>; ++d)
          stiffness_matrix(nsd<distype> * inod + d, nsd<distype> * jnod + d) +=
              jacobian_mapping.integration_factor_ * bopstrbop;
      }
    }
  }

  /*!
   * @brief Add mass matrix contribution of one Gauss point
   *
   * @tparam distype
   * @param shapefunctions (in) : Shape functions and derivatives evaluated at the respective point
   * in the parameter space
   * @param integration_factor (in) : Integration factor (Gauss point weight times the determinant
   * of the jacobian)
   * @param density (in) : density at the Gauss point
   * @param mass (in/out) : mass matrix where the local contribution is added to
   */
  template <DRT::Element::DiscretizationType distype>
  void AddMassMatrix(const ShapeFunctionsAndDerivatives<distype>& shapefunctions,
      const double integration_factor, const double density,
      LINALG::Matrix<nsd<distype> * nen<distype>, nsd<distype> * nen<distype>>& mass)
  {
    for (int inod = 0; inod < nen<distype>; ++inod)
    {
      const double ifactor = shapefunctions.shapefunctions_(inod) * integration_factor * density;
      for (int jnod = 0; jnod < nen<distype>; ++jnod)
      {
        const double massfactor =
            shapefunctions.shapefunctions_(jnod) * ifactor;  // intermediate factor
        for (int d = 0; d < nsd<distype>; ++d)
          mass(nsd<distype> * inod + d, nsd<distype> * jnod + d) += massfactor;
      }
    }
  }

  template <unsigned numstr>
  void AssembleVectorToMatrixRow(
      LINALG::Matrix<numstr, 1> vector, Epetra_SerialDenseMatrix& data, int row)
  {
    for (unsigned i = 0; i < numstr; ++i) data(row, i) = vector(i);
  }

  template <DRT::Element::DiscretizationType distype>
  void AssembleStrainTypeToMatrixRow(const Strains<distype>& strains,
      INPAR::STR::StrainType strain_type, Epetra_SerialDenseMatrix& data, int row)
  {
    switch (strain_type)
    {
      case INPAR::STR::strain_gl:
      {
        LINALG::Matrix<numstr<distype>, 1> gl_strain_stress_like;
        ::UTILS::VOIGT::Strains::ToStressLike(strains.gl_strain_, gl_strain_stress_like);
        AssembleVectorToMatrixRow(gl_strain_stress_like, data, row);
        return;
      }
      case INPAR::STR::strain_ea:
      {
        const LINALG::Matrix<numstr<distype>, 1> ea =
            STR::UTILS::GreenLagrangeToEulerAlmansi(strains.gl_strain_, strains.defgrd_);
        LINALG::Matrix<numstr<distype>, 1> ea_stress_like;
        ::UTILS::VOIGT::Strains::ToStressLike(ea, ea_stress_like);
        AssembleVectorToMatrixRow(ea_stress_like, data, row);
        return;
      }
      case INPAR::STR::strain_none:
        return;
      default:
        dserror("strain type not supported");
        break;
    }
  }


  template <DRT::Element::DiscretizationType distype>
  void AssembleStressTypeToMatrixRow(const Strains<distype>& strains, const Stress<distype> stress,
      INPAR::STR::StressType stress_type, Epetra_SerialDenseMatrix& data, int row)
  {
    switch (stress_type)
    {
      case INPAR::STR::stress_2pk:
      {
        AssembleVectorToMatrixRow(stress.pk2_, data, row);
        return;
      }
      case INPAR::STR::stress_cauchy:
      {
        LINALG::Matrix<numstr<distype>, 1> cauchy;
        STR::UTILS::Pk2ToCauchy(stress.pk2_, strains.defgrd_, cauchy);
        AssembleVectorToMatrixRow(cauchy, data, row);
        return;
      }
      case INPAR::STR::stress_none:

        return;
      default:
        dserror("stress type not supported");
        break;
    }
  }

  void Serialize(const Epetra_SerialDenseMatrix& matrix, std::vector<char>& data)
  {
    DRT::PackBuffer packBuffer;
    DRT::ParObject::AddtoPack(packBuffer, matrix);
    packBuffer.StartPacking();
    DRT::ParObject::AddtoPack(packBuffer, matrix);
    std::copy(packBuffer().begin(), packBuffer().end(), std::back_inserter(data));
  }
}  // namespace DRT::ELEMENTS

#endif  // SOLID_ELE_CALC_LIB_H