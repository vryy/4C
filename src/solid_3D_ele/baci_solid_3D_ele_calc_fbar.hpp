/*! \file

\brief A displacement based solid element formulation with FBAR element technology

\level 1
*/

#ifndef BACI_SOLID_3D_ELE_CALC_FBAR_HPP
#define BACI_SOLID_3D_ELE_CALC_FBAR_HPP

#include "baci_config.hpp"

#include "baci_discretization_fem_general_cell_type_traits.hpp"
#include "baci_lib_element.hpp"
#include "baci_solid_3D_ele_calc.hpp"
#include "baci_solid_3D_ele_calc_lib.hpp"

BACI_NAMESPACE_OPEN

namespace DRT::ELEMENTS
{
  namespace DETAILS
  { /*!
     * @brief Evaluate the fbar factor \f[ \frac{\mathbf{F}_{\mathrm{centroid}}}{\mathbf{F}}^{1/3}
     * \f]
     *
     * @param defgrd_centroid (in) : Deformation gradient evaluated at the element centroid
     * @param defgrd_gp (in) : Deformation gradient evaluated at the Gauss point
     * @return double : Fbar factor
     */
    inline double EvaluateFbarFactor(const double& defgrd_centroid, const double& defgrd_gp)
    {
      const double fbar_factor = std::pow(defgrd_centroid / defgrd_gp, 1.0 / 3.0);
      return fbar_factor;
    }

    /*!
     * @brief Evaluates the H-Operator used in F-bar of the specified element
     *
     * @tparam celltype : Cell type
     * @param jacobian_mapping (in) : Quantities of the jacobian mapping evaluated at the Gauss
     * point
     * @param jacobian_mapping_centroid (in) : Quantities of the jacobian mapping evaluated at the
     * element centroid
     * @param spatial_material_mapping (in) :An object holding quantities of the spatial material
     * mapping (deformation_gradient, inverse_deformation_gradient,
     * determinant_deformation_gradient) evaluated at the Gauss point
     * @param spatial_material_mapping_centroid (in) : An object holding quantities of the spatial
     * material mapping (deformation_gradient, inverse_deformation_gradient,
     * determinant_deformation_gradient) evaluated at the element centroid
     * @return CORE::LINALG::Matrix<num_dof_per_ele, 1> : H-Operator
     */
    template <CORE::FE::CellType celltype, std::enable_if_t<CORE::FE::dim<celltype> == 3, int> = 0>
    inline CORE::LINALG::Matrix<CORE::FE::dim<celltype> * CORE::FE::num_nodes<celltype>, 1>
    EvaluateFbarHOperator(const DRT::ELEMENTS::JacobianMapping<celltype>& jacobian_mapping,
        const DRT::ELEMENTS::JacobianMapping<celltype>& jacobian_mapping_centroid,
        const DRT::ELEMENTS::SpatialMaterialMapping<celltype> spatial_material_mapping,
        const DRT::ELEMENTS::SpatialMaterialMapping<celltype> spatial_material_mapping_centroid)
    {
      // inverse deformation gradient at centroid
      CORE::LINALG::Matrix<CORE::FE::dim<celltype>, CORE::FE::dim<celltype>> invdefgrd_centroid;
      invdefgrd_centroid.Invert(spatial_material_mapping_centroid.deformation_gradient_);

      // inverse deformation gradient at gp
      CORE::LINALG::Matrix<CORE::FE::dim<celltype>, CORE::FE::dim<celltype>> invdefgrd;
      invdefgrd.Invert(spatial_material_mapping.deformation_gradient_);

      CORE::LINALG::Matrix<CORE::FE::dim<celltype> * CORE::FE::num_nodes<celltype>, 1> Hop(true);
      for (int idof = 0; idof < CORE::FE::dim<celltype> * CORE::FE::num_nodes<celltype>; idof++)
      {
        for (int idim = 0; idim < CORE::FE::dim<celltype>; idim++)
        {
          Hop(idof) += invdefgrd_centroid(idim, idof % CORE::FE::dim<celltype>) *
                       jacobian_mapping_centroid.N_XYZ_(idim, idof / CORE::FE::dim<celltype>);
          Hop(idof) -= invdefgrd(idim, idof % CORE::FE::dim<celltype>) *
                       jacobian_mapping.N_XYZ_(idim, idof / CORE::FE::dim<celltype>);
        }
      }

      return Hop;
    }

    /*!
     * @brief Add fbar stiffness matrix contribution of one Gauss point
     *
     * @tparam celltype : Cell type
     * @param Bop (in) : Strain gradient (B-Operator)
     * @param Hop (in) : H-Operator
     * @param f_bar_factor (in) : f_bar_factor
     * @param integration_fac (in) : Integration factor (Gauss point weight times the determinant of
     * the jacobian)
     * @param cauchyGreen (in) : An object holding the right Cauchy-Green deformation tensor and
     * its inverse
     * @param stress_bar (in) : Deviatoric part of stress measures
     * @param stiffness_matrix (in/out) : stiffness matrix where the local contribution is added to
     */
    template <CORE::FE::CellType celltype>
    inline void AddFbarStiffnessMatrix(
        const CORE::LINALG::Matrix<num_str<celltype>,
            CORE::FE::dim<celltype> * CORE::FE::num_nodes<celltype>>& Bop,
        const CORE::LINALG::Matrix<CORE::FE::dim<celltype> * CORE::FE::num_nodes<celltype>, 1>& Hop,
        const double f_bar_factor, const double integration_fac,
        const CORE::LINALG::Matrix<CORE::FE::dim<celltype>, CORE::FE::dim<celltype>> cauchyGreen,
        const DRT::ELEMENTS::Stress<celltype> stress_bar,
        CORE::LINALG::Matrix<CORE::FE::dim<celltype> * CORE::FE::num_nodes<celltype>,
            CORE::FE::dim<celltype> * CORE::FE::num_nodes<celltype>>& stiffness_matrix)
    {
      constexpr int num_dof_per_ele = CORE::FE::dim<celltype> * CORE::FE::num_nodes<celltype>;

      CORE::LINALG::Matrix<num_str<celltype>, 1> rcg_bar_voigt;
      CORE::LINALG::VOIGT::Strains::MatrixToVector(cauchyGreen, rcg_bar_voigt);

      CORE::LINALG::Matrix<num_str<celltype>, 1> ccg;
      ccg.MultiplyNN(stress_bar.cmat_, rcg_bar_voigt);

      // auxiliary integrated stress_bar
      CORE::LINALG::Matrix<num_dof_per_ele, 1> bopccg(false);
      bopccg.MultiplyTN(integration_fac * f_bar_factor / 3.0, Bop, ccg);

      CORE::LINALG::Matrix<num_dof_per_ele, 1> bops(false);
      bops.MultiplyTN(-integration_fac / f_bar_factor / 3.0, Bop, stress_bar.pk2_);

      for (int idof = 0; idof < num_dof_per_ele; idof++)
      {
        for (int jdof = 0; jdof < num_dof_per_ele; jdof++)
        {
          stiffness_matrix(idof, jdof) += Hop(jdof) * (bops(idof, 0) + bopccg(idof, 0));
        }
      }
    }
  }  // namespace DETAILS

  template <CORE::FE::CellType celltype>
  struct FBarPreparationData
  {
    /// jacobian mapping evaluated at element centroid
    JacobianMapping<celltype> jacobian_mapping_centroid;

    /// deformation gradient at element centroid
    SpatialMaterialMapping<celltype> spatial_material_mapping_centroid;
  };

  struct FBarHistoryData
  {
    // no history data needed
  };

  template <CORE::FE::CellType celltype>
  struct FBarLinearizationContainer
  {
    CORE::LINALG::Matrix<DETAILS::num_str<celltype>,
        CORE::FE::num_nodes<celltype> * CORE::FE::dim<celltype>>
        Bop{};

    CORE::LINALG::Matrix<CORE::FE::num_nodes<celltype> * CORE::FE::dim<celltype>, 1> Hop{};

    CORE::LINALG::Matrix<DETAIL::num_dim<celltype>, DETAIL::num_dim<celltype>> cauchygreen{};

    double fbar_factor = 1.0;
  };

  /*!
   * @brief A displacement based solid element formulation with FBAR element technology
   *
   * @tparam celltype
   */
  template <CORE::FE::CellType celltype>
  struct FBarFormulation
  {
    static FBarPreparationData<celltype> Prepare(const DRT::Element& ele,
        const ElementNodes<celltype>& nodal_coordinates, FBarHistoryData& history_data)
    {
      const JacobianMapping<celltype> jacobian_mapping_centroid =
          EvaluateJacobianMappingCentroid(nodal_coordinates);

      return {jacobian_mapping_centroid,
          EvaluateSpatialMaterialMapping(jacobian_mapping_centroid, nodal_coordinates)};
    }

    template <typename Evaluator>
    static auto Evaluate(const DRT::Element& ele, const ElementNodes<celltype>& nodal_coordinates,
        const CORE::LINALG::Matrix<DETAIL::num_dim<celltype>, 1>& xi,
        const ShapeFunctionsAndDerivatives<celltype>& shape_functions,
        const JacobianMapping<celltype>& jacobian_mapping,
        const FBarPreparationData<celltype>& preparation_data, FBarHistoryData& history_data,
        Evaluator evaluator)
    {
      const SpatialMaterialMapping<celltype> spatial_material_mapping =
          EvaluateSpatialMaterialMapping(jacobian_mapping, nodal_coordinates);

      // factor (detF0/detF)^1/3
      const double fbar_factor = DETAILS::EvaluateFbarFactor(
          preparation_data.spatial_material_mapping_centroid.determinant_deformation_gradient_,
          spatial_material_mapping.determinant_deformation_gradient_);

      const FBarLinearizationContainer<celltype> linearization = std::invoke(
          [&]()
          {
            FBarLinearizationContainer<celltype> linearization{};
            linearization.Bop = EvaluateStrainGradient(jacobian_mapping, spatial_material_mapping);

            linearization.Hop = DETAILS::EvaluateFbarHOperator(jacobian_mapping,
                preparation_data.jacobian_mapping_centroid, spatial_material_mapping,
                preparation_data.spatial_material_mapping_centroid);

            linearization.fbar_factor = fbar_factor;

            linearization.cauchygreen = EvaluateCauchyGreen(spatial_material_mapping);

            return linearization;
          });

      // deformation gradient F_bar and resulting strains: F_bar = (detF_0/detF)^1/3 F
      const SpatialMaterialMapping<celltype> spatial_material_mapping_bar =
          EvaluateSpatialMaterialMapping(jacobian_mapping, nodal_coordinates, fbar_factor);

      const CORE::LINALG::Matrix<CORE::FE::dim<celltype>, CORE::FE::dim<celltype>> cauchygreen_bar =
          EvaluateCauchyGreen(spatial_material_mapping_bar);

      CORE::LINALG::Matrix<DETAIL::num_str<celltype>, 1> gl_strain_bar =
          EvaluateGreenLagrangeStrain(cauchygreen_bar);

      return evaluator(
          spatial_material_mapping_bar.deformation_gradient_, gl_strain_bar, linearization);
    }

    static inline SolidFormulationLinearization<celltype> EvaluateFullLinearization(
        const DRT::Element& ele, const ElementNodes<celltype>& nodal_coordinates,
        const CORE::LINALG::Matrix<DETAIL::num_dim<celltype>, 1>& xi,
        const ShapeFunctionsAndDerivatives<celltype>& shape_functions,
        const JacobianMapping<celltype>& jacobian_mapping,
        const CORE::LINALG::Matrix<DETAIL::num_dim<celltype>, DETAIL::num_dim<celltype>>&
            deformation_gradient,
        const FBarPreparationData<celltype>& preparation_data, FBarHistoryData& history_data)
    {
      dserror(
          "The full linearization is not yet implemented for the displacement based formulation "
          "with fbar.");
    }

    static CORE::LINALG::Matrix<DETAILS::num_str<celltype>,
        CORE::FE::num_nodes<celltype> * CORE::FE::dim<celltype>>
    GetLinearBOperator(const FBarLinearizationContainer<celltype>& linearization)
    {
      return linearization.Bop;
    }

    static void AddInternalForceVector(const FBarLinearizationContainer<celltype>& linearization,
        const Stress<celltype>& stress, const double integration_factor,
        const FBarPreparationData<celltype>& preparation_data, FBarHistoryData& history_data,
        CORE::LINALG::Matrix<CORE::FE::num_nodes<celltype> * CORE::FE::dim<celltype>, 1>&
            force_vector)
    {
      DRT::ELEMENTS::AddInternalForceVector(
          linearization.Bop, stress, integration_factor / linearization.fbar_factor, force_vector);
    }

    static void AddStiffnessMatrix(const FBarLinearizationContainer<celltype>& linearization,
        const JacobianMapping<celltype>& jacobian_mapping, const Stress<celltype>& stress,
        const double integration_factor, const FBarPreparationData<celltype>& preparation_data,
        FBarHistoryData& history_data,
        CORE::LINALG::Matrix<CORE::FE::num_nodes<celltype> * CORE::FE::dim<celltype>,
            CORE::FE::num_nodes<celltype> * CORE::FE::dim<celltype>>& stiffness_matrix)
    {
      DRT::ELEMENTS::AddElasticStiffnessMatrix(linearization.Bop, stress,
          integration_factor * linearization.fbar_factor, stiffness_matrix);
      DRT::ELEMENTS::AddGeometricStiffnessMatrix(jacobian_mapping.N_XYZ_, stress,
          integration_factor / linearization.fbar_factor, stiffness_matrix);

      // additional stiffness matrix needed for fbar method
      DETAILS::AddFbarStiffnessMatrix(linearization.Bop, linearization.Hop,
          linearization.fbar_factor, integration_factor, linearization.cauchygreen, stress,
          stiffness_matrix);
    }

    static void Pack(const FBarHistoryData& history_data, CORE::COMM::PackBuffer& data)
    {
      // nothing to pack
    }

    static void Unpack(std::vector<char>::size_type& position, const std::vector<char>& data,
        FBarHistoryData& history_data)
    {
      // nothing to unpack
    }
  };

  template <CORE::FE::CellType celltype>
  using FBarSolidIntegrator = SolidEleCalc<celltype, FBarFormulation<celltype>,
      FBarPreparationData<celltype>, FBarHistoryData>;


}  // namespace DRT::ELEMENTS

BACI_NAMESPACE_CLOSE
#endif  // BACI_SOLID_3D_ELE_CALC_FBAR_HPP
