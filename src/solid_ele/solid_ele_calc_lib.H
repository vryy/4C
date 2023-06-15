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
#include "lib_element.H"
#include "discretization_fem_general_utils_gausspoints.H"
#include "discretization_fem_general_utils_local_connectivity_matrices.H"
#include "solid_ele.H"
#include "lib_utils.H"
#include "mat_so3_material.H"
#include "lib_voigt_notation.H"


namespace DRT::ELEMENTS::DETAIL
{
  template <DRT::Element::DiscretizationType distype>
  inline static constexpr int nen =
      CORE::DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;

  template <DRT::Element::DiscretizationType distype>
  inline static constexpr int nsd = CORE::DRT::UTILS::DisTypeToDim<distype>::dim;

  template <DRT::Element::DiscretizationType distype>
  inline static constexpr int numstr = nsd<distype>*(nsd<distype> + 1) / 2;

  template <DRT::Element::DiscretizationType distype>
  inline static constexpr int numdofperelement = nen<distype>* nsd<distype>;
}  // namespace DRT::ELEMENTS::DETAIL

namespace DRT::ELEMENTS
{
  template <DRT::Element::DiscretizationType distype>
  constexpr auto GetGaussRuleMassMatrix()
  {
    return DRT::ELEMENTS::DisTypeToOptGaussRule<distype>::rule;
  }

  template <DRT::Element::DiscretizationType distype>
  constexpr auto GetGaussRuleStiffnessMatrix()
  {
    return DRT::ELEMENTS::DisTypeToOptGaussRule<distype>::rule;
  }

  template <>
  constexpr auto GetGaussRuleStiffnessMatrix<DRT::Element::DiscretizationType::tet10>()
  {
    return CORE::DRT::UTILS::GaussRule3D::tet_4point;
  }

  template <>
  constexpr auto GetGaussRuleStiffnessMatrix<DRT::Element::DiscretizationType::tet4>()
  {
    return CORE::DRT::UTILS::GaussRule3D::tet_1point;
  }

  // TODO: What about the stiffness matrix of tri-elements?
  template <DRT::Element::DiscretizationType distype, typename GaussRuleType>
  CORE::DRT::UTILS::GaussIntegration CreateGaussIntegration(GaussRuleType rule)
  {
    constexpr int nsd = CORE::DRT::UTILS::DisTypeToDim<distype>::dim;

    // setup default integration
    CORE::DRT::UTILS::IntPointsAndWeights<nsd> intpoints(rule);

    // format as DRT::UTILS::GaussIntegration
    Teuchos::RCP<CORE::DRT::UTILS::CollectedGaussPoints> gp =
        Teuchos::rcp(new CORE::DRT::UTILS::CollectedGaussPoints);

    std::array<double, 3> xi = {0., 0., 0.};
    for (int i = 0; i < intpoints.IP().nquad; ++i)
    {
      for (int d = 0; d < nsd; ++d) xi[d] = intpoints.IP().qxg[i][d];
      gp->Append(xi[0], xi[1], xi[2], intpoints.IP().qwgt[i]);
    }

    return CORE::DRT::UTILS::GaussIntegration(gp);
  }

  template <DRT::Element::DiscretizationType distype>
  struct NodalCoordinates
  {
    LINALG::Matrix<DETAIL::nen<distype>, DETAIL::nsd<distype>> reference_;
    LINALG::Matrix<DETAIL::nen<distype>, DETAIL::nsd<distype>> current_;
  };

  template <DRT::Element::DiscretizationType distype>
  NodalCoordinates<distype> EvaluateNodalCoordinates(const DRT::Element& ele,
      const DRT::Discretization& discretization, const std::vector<int>& lm)
  {
    const Epetra_Vector& displacements = *discretization.GetState("displacement");
    std::vector<double> mydisp(lm.size());
    DRT::UTILS::ExtractMyValues(displacements, mydisp, lm);

    NodalCoordinates<distype> coordinates;
    for (int i = 0; i < DETAIL::nen<distype>; ++i)
    {
      for (int d = 0; d < DETAIL::nsd<distype>; ++d)
      {
        coordinates.reference_(i, d) = ele.Nodes()[i]->X()[d];
        coordinates.current_(i, d) =
            coordinates.reference_(i, d) + mydisp[i * DETAIL::nsd<distype> + d];
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
  LINALG::Matrix<DETAIL::nsd<distype>, 1> EvaluateParameterCoordinate(
      const CORE::DRT::UTILS::GaussIntegration& intpoints, const int gp)
  {
    LINALG::Matrix<DETAIL::nsd<distype>, 1> xi;
    for (int d = 0; d < DETAIL::nsd<distype>; ++d) xi(d) = intpoints.Point(gp)[d];

    return xi;
  }

  template <DRT::Element::DiscretizationType distype>
  struct ShapeFunctionsAndDerivatives
  {
    LINALG::Matrix<DETAIL::nen<distype>, 1> shapefunctions_;
    LINALG::Matrix<DETAIL::nsd<distype>, DETAIL::nen<distype>> derivatives_;
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
      const LINALG::Matrix<DETAIL::nsd<distype>, 1>& xi)
  {
    ShapeFunctionsAndDerivatives<distype> shapefcns;
    CORE::DRT::UTILS::shape_function<distype>(xi, shapefcns.shapefunctions_);
    CORE::DRT::UTILS::shape_function_deriv1<distype>(xi, shapefcns.derivatives_);

    return shapefcns;
  }

  template <DRT::Element::DiscretizationType distype>
  struct JacobianMapping
  {
    double determinant_;
    LINALG::Matrix<DETAIL::nsd<distype>, DETAIL::nsd<distype>> inverse_jacobian_;
    LINALG::Matrix<DETAIL::nsd<distype>, DETAIL::nen<distype>> n_xyz_;
  };

  /*!
   * @brief Evaluates the jacobian mapping of the element
   *
   * @tparam distype : Discretization type
   * @param shapefcns (in) : Shape functions and derivatives evaluated at the respective point in
   * the parameter space
   * @param nodal_coordniates (in) : Reference and current coordinates of the nodes of the element
   * @param gp (in) : Id of the Gauss point
   * @return JacobianMapping<distype> : An object holding quantities of the jacobian mapping
   * (inverse Jacobian, determinant, derivatives of the shape functions w.r.t. XYZ, integration
   * factor)
   */
  template <DRT::Element::DiscretizationType distype>
  JacobianMapping<distype> EvaluateJacobianMapping(
      const ShapeFunctionsAndDerivatives<distype>& shapefcns,
      const NodalCoordinates<distype>& nodal_coordniates)
  {
    JacobianMapping<distype> jacobian;

    jacobian.inverse_jacobian_.Multiply(shapefcns.derivatives_, nodal_coordniates.reference_);
    jacobian.determinant_ = jacobian.inverse_jacobian_.Invert();
    jacobian.n_xyz_.Multiply(jacobian.inverse_jacobian_, shapefcns.derivatives_);

    return jacobian;
  }

  template <DRT::Element::DiscretizationType distype>
  double EvaluateJacobianDeterminant(const ShapeFunctionsAndDerivatives<distype>& shapefcns,
      const NodalCoordinates<distype>& nodal_coordniates)
  {
    LINALG::Matrix<DETAIL::nsd<distype>, DETAIL::nsd<distype>> jacobian;
    jacobian.Multiply(shapefcns.derivatives_, nodal_coordniates.reference_);

    return jacobian.Determinant();
  }

  template <DRT::Element::DiscretizationType distype>
  struct Strains
  {
    LINALG::Matrix<DETAIL::nsd<distype>, DETAIL::nsd<distype>> defgrd_;
    LINALG::Matrix<DETAIL::nsd<distype>, DETAIL::nsd<distype>> rcg_;
    LINALG::Matrix<DETAIL::numstr<distype>, 1> gl_strain_;
  };

  /*!
   * @brief Evaluates Green Lagrange strain from right cauchy green tensor
   *
   * GL strain vector glstrain={E11,E22,E33,2*E12,2*E23,2*E31}
   * @tparam distype
   * @param strains (in/out) : Strain measures
   * Cauchy-Green deformation tensor, Green-Lagrange strain tensor)
   */
  template <DRT::Element::DiscretizationType distype,
      std::enable_if_t<DETAIL::nsd<distype> == 3, int> = 0>
  void EvaluateGreenLagrangeStrain(Strains<distype>& strains)
  {
    strains.gl_strain_(0) = 0.5 * (strains.rcg_(0, 0) - 1.0);
    strains.gl_strain_(1) = 0.5 * (strains.rcg_(1, 1) - 1.0);
    strains.gl_strain_(2) = 0.5 * (strains.rcg_(2, 2) - 1.0);
    strains.gl_strain_(3) = strains.rcg_(0, 1);
    strains.gl_strain_(4) = strains.rcg_(1, 2);
    strains.gl_strain_(5) = strains.rcg_(2, 0);
  }

  /*!
   * @brief Evaluates strain measures
   *
   * @tparam distype
   * @param nodal_coordinates (in) : Reference and current coordinates of the nodes of the element
   * @param jacobian_mapping (in) : Quantities of the jacobian mapping
   * @param scale_defgrd (in) : Optional scalar factor for deformation gradient
   * @return Strains<distype> : Strain measures of the element (deformation gradient, right
   * Cauchy-Green deformation tensor, Green-Lagrange strain tensor)
   */
  template <DRT::Element::DiscretizationType distype,
      std::enable_if_t<DETAIL::nsd<distype> == 3, int> = 0>
  Strains<distype> EvaluateStrains(const NodalCoordinates<distype>& nodal_coordinates,
      const JacobianMapping<distype>& jacobian_mapping, double scale_defgrd = 1.0)
  {
    Strains<distype> strains;

    strains.defgrd_.MultiplyTT(scale_defgrd, nodal_coordinates.current_, jacobian_mapping.n_xyz_);
    strains.rcg_.MultiplyTN(strains.defgrd_, strains.defgrd_);
    EvaluateGreenLagrangeStrain(strains);

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
  template <DRT::Element::DiscretizationType distype,
      std::enable_if_t<DETAIL::nsd<distype> == 3, int> = 0>
  LINALG::Matrix<DETAIL::numstr<distype>, DETAIL::nsd<distype> * DETAIL::nen<distype>>
  EvaluateStrainGradient(
      const JacobianMapping<distype>& jacobian_mapping, const Strains<distype>& strains)
  {
    // B-operator
    LINALG::Matrix<DETAIL::numstr<distype>, DETAIL::nsd<distype> * DETAIL::nen<distype>> Bop;
    for (int i = 0; i < DETAIL::nen<distype>; ++i)
    {
      for (int d = 0; d < DETAIL::nsd<distype>; ++d)
        for (int e = 0; e < DETAIL::nsd<distype>; ++e)
          Bop(d, DETAIL::nsd<distype> * i + e) =
              strains.defgrd_(e, d) * jacobian_mapping.n_xyz_(d, i);

      // TODO: make this work for 2D
      Bop(3, DETAIL::nsd<distype> * i + 0) = strains.defgrd_(0, 0) * jacobian_mapping.n_xyz_(1, i) +
                                             strains.defgrd_(0, 1) * jacobian_mapping.n_xyz_(0, i);
      Bop(3, DETAIL::nsd<distype> * i + 1) = strains.defgrd_(1, 0) * jacobian_mapping.n_xyz_(1, i) +
                                             strains.defgrd_(1, 1) * jacobian_mapping.n_xyz_(0, i);
      Bop(3, DETAIL::nsd<distype> * i + 2) = strains.defgrd_(2, 0) * jacobian_mapping.n_xyz_(1, i) +
                                             strains.defgrd_(2, 1) * jacobian_mapping.n_xyz_(0, i);
      Bop(4, DETAIL::nsd<distype> * i + 0) = strains.defgrd_(0, 1) * jacobian_mapping.n_xyz_(2, i) +
                                             strains.defgrd_(0, 2) * jacobian_mapping.n_xyz_(1, i);
      Bop(4, DETAIL::nsd<distype> * i + 1) = strains.defgrd_(1, 1) * jacobian_mapping.n_xyz_(2, i) +
                                             strains.defgrd_(1, 2) * jacobian_mapping.n_xyz_(1, i);
      Bop(4, DETAIL::nsd<distype> * i + 2) = strains.defgrd_(2, 1) * jacobian_mapping.n_xyz_(2, i) +
                                             strains.defgrd_(2, 2) * jacobian_mapping.n_xyz_(1, i);
      Bop(5, DETAIL::nsd<distype> * i + 0) = strains.defgrd_(0, 2) * jacobian_mapping.n_xyz_(0, i) +
                                             strains.defgrd_(0, 0) * jacobian_mapping.n_xyz_(2, i);
      Bop(5, DETAIL::nsd<distype> * i + 1) = strains.defgrd_(1, 2) * jacobian_mapping.n_xyz_(0, i) +
                                             strains.defgrd_(1, 0) * jacobian_mapping.n_xyz_(2, i);
      Bop(5, DETAIL::nsd<distype> * i + 2) = strains.defgrd_(2, 2) * jacobian_mapping.n_xyz_(0, i) +
                                             strains.defgrd_(2, 0) * jacobian_mapping.n_xyz_(2, i);
    }

    return Bop;
  }

  template <DRT::Element::DiscretizationType distype>
  struct Stress
  {
    LINALG::Matrix<DETAIL::numstr<distype>, 1> pk2_;
    LINALG::Matrix<DETAIL::numstr<distype>, DETAIL::numstr<distype>> cmat_;
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
      const LINALG::Matrix<DETAIL::numstr<distype>, DETAIL::nsd<distype> * DETAIL::nen<distype>>&
          Bop,
      const Stress<distype> stress, const double integration_fac,
      LINALG::Matrix<DETAIL::nsd<distype> * DETAIL::nen<distype>, 1>& force_vector)
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
      const LINALG::Matrix<DETAIL::numstr<distype>, DETAIL::nsd<distype> * DETAIL::nen<distype>>&
          Bop,
      const Stress<distype> stress, const double integration_fac,
      LINALG::Matrix<DETAIL::nsd<distype> * DETAIL::nen<distype>,
          DETAIL::nsd<distype> * DETAIL::nen<distype>>& stiffness_matrix)
  {
    LINALG::Matrix<DETAIL::numstr<distype>, DETAIL::nen<distype> * DETAIL::nsd<distype>> cb;
    cb.Multiply(stress.cmat_, Bop);
    stiffness_matrix.MultiplyTN(integration_fac, Bop, cb, 1.0);
  }

  /*!
   * @brief Add geometric stiffness matrix contribution of one Gauss point
   *
   * @tparam distype
   * @param B_L (in) : B_L operator, i.e. derivatives of the shape functions w.r.t. XYZ
   *                   at the respective Gauss point
   * @param stress (in) : Stress measures
   * @param integration_fac (in) : Integration factor (Gauss point weight times the determinant of
   * the jacobian)
   * @param stiffness_matrix (in/out) : stiffness matrix where the local contribution is added to
   */
  template <DRT::Element::DiscretizationType distype>
  void AddGeometricStiffnessMatrix(
      const LINALG::Matrix<DETAIL::nsd<distype>, DETAIL::nen<distype>>& B_L,
      const Stress<distype> stress, const double integration_fac,
      LINALG::Matrix<DETAIL::nsd<distype> * DETAIL::nen<distype>,
          DETAIL::nsd<distype> * DETAIL::nen<distype>>& stiffness_matrix)
  {
    std::array<double, 3> SmB_L;  // intermediate Sm.B_L
    // kgeo += (B_L^T . sigma . B_L) * detJ * w(gp)  with B_L = Ni,Xj see NiliFEM-Skript
    for (int inod = 0; inod < DETAIL::nen<distype>; ++inod)
    {
      SmB_L[0] = stress.pk2_(0) * B_L(0, inod) + stress.pk2_(3) * B_L(1, inod) +
                 stress.pk2_(5) * B_L(2, inod);
      SmB_L[1] = stress.pk2_(3) * B_L(0, inod) + stress.pk2_(1) * B_L(1, inod) +
                 stress.pk2_(4) * B_L(2, inod);
      SmB_L[2] = stress.pk2_(5) * B_L(0, inod) + stress.pk2_(4) * B_L(1, inod) +
                 stress.pk2_(2) * B_L(2, inod);

      for (int jnod = 0; jnod < DETAIL::nen<distype>; ++jnod)
      {
        double bopstrbop = 0.0;  // intermediate value
        for (int idim = 0; idim < DETAIL::nsd<distype>; ++idim)
          bopstrbop += B_L(idim, jnod) * SmB_L[idim];

        for (int d = 0; d < DETAIL::nsd<distype>; ++d)
          stiffness_matrix(DETAIL::nsd<distype> * inod + d, DETAIL::nsd<distype> * jnod + d) +=
              integration_fac * bopstrbop;
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
      LINALG::Matrix<DETAIL::nsd<distype> * DETAIL::nen<distype>,
          DETAIL::nsd<distype> * DETAIL::nen<distype>>& mass)
  {
    for (int inod = 0; inod < DETAIL::nen<distype>; ++inod)
    {
      const double ifactor = shapefunctions.shapefunctions_(inod) * integration_factor * density;
      for (int jnod = 0; jnod < DETAIL::nen<distype>; ++jnod)
      {
        const double massfactor =
            shapefunctions.shapefunctions_(jnod) * ifactor;  // intermediate factor
        for (int d = 0; d < DETAIL::nsd<distype>; ++d)
          mass(DETAIL::nsd<distype> * inod + d, DETAIL::nsd<distype> * jnod + d) += massfactor;
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
        LINALG::Matrix<DETAIL::numstr<distype>, 1> gl_strain_stress_like;
        ::UTILS::VOIGT::Strains::ToStressLike(strains.gl_strain_, gl_strain_stress_like);
        AssembleVectorToMatrixRow(gl_strain_stress_like, data, row);
        return;
      }
      case INPAR::STR::strain_ea:
      {
        const LINALG::Matrix<DETAIL::numstr<distype>, 1> ea =
            STR::UTILS::GreenLagrangeToEulerAlmansi(strains.gl_strain_, strains.defgrd_);
        LINALG::Matrix<DETAIL::numstr<distype>, 1> ea_stress_like;
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
        LINALG::Matrix<DETAIL::numstr<distype>, 1> cauchy;
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

  inline void Serialize(const Epetra_SerialDenseMatrix& matrix, std::vector<char>& data)
  {
    DRT::PackBuffer packBuffer;
    DRT::ParObject::AddtoPack(packBuffer, matrix);
    packBuffer.StartPacking();
    DRT::ParObject::AddtoPack(packBuffer, matrix);
    std::copy(packBuffer().begin(), packBuffer().end(), std::back_inserter(data));
  }

  /*!
   * @brief Iterates over the jacobian mapping at each integration point of @p integration.
   *
   * @tparam distype The discretization type known at compile time.
   * @param nodal_coordinates (in) : The nodal coordinates of the element
   * @param integration (in) : The integration rule to be used.
   * @param gp_evaluator (in) : A callable object (e.g. lambda-function) with signature void(const
   * LINALG::Matrix<DETAIL::nsd<distype>, 1>& xi, const ShapeFunctionsAndDerivatives<distype>&
   * shape_functions, const JacobianMapping<distype>& jacobian_mapping, double integration_factor,
   * int gp) that will be called for each integration point.
   */
  template <DRT::Element::DiscretizationType distype, typename GaussPointEvaluator>
  inline void IterateJacobianMappingAtGaussPoints(
      const NodalCoordinates<distype>& nodal_coordinates,
      const CORE::DRT::UTILS::GaussIntegration& integration, GaussPointEvaluator gp_evaluator)
  {
    for (int gp = 0; gp < integration.NumPoints(); ++gp)
    {
      const LINALG::Matrix<DETAIL::nsd<distype>, 1> xi =
          EvaluateParameterCoordinate<distype>(integration, gp);

      const ShapeFunctionsAndDerivatives<distype> shape_functions =
          EvaluateShapeFunctionsAndDerivs<distype>(xi);

      const JacobianMapping<distype> jacobian_mapping =
          EvaluateJacobianMapping(shape_functions, nodal_coordinates);

      const double integration_factor = jacobian_mapping.determinant_ * integration.Weight(gp);

      gp_evaluator(xi, shape_functions, jacobian_mapping, integration_factor, gp);
    }
  }
}  // namespace DRT::ELEMENTS

#endif  // SOLID_ELE_CALC_LIB_H