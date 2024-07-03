/*! \file

\brief A library of free functions for a default solid element

\level 1
*/

#ifndef FOUR_C_SOLID_3D_ELE_CALC_LIB_IO_HPP
#define FOUR_C_SOLID_3D_ELE_CALC_LIB_IO_HPP

#include "4C_config.hpp"

#include "4C_fem_general_cell_type.hpp"
#include "4C_fem_general_cell_type_traits.hpp"
#include "4C_fem_general_utils_gauss_point_extrapolation.hpp"
#include "4C_inpar_structure.hpp"
#include "4C_linalg_fixedsizematrix.hpp"
#include "4C_linalg_fixedsizematrix_voigt_notation.hpp"
#include "4C_so3_element_service.hpp"
#include "4C_solid_3D_ele_calc_lib.hpp"
#include "4C_solid_3D_ele_utils.hpp"
#include "4C_structure_new_gauss_point_data_output_manager.hpp"

#include <Teuchos_ParameterList.hpp>

FOUR_C_NAMESPACE_OPEN

namespace Discret::ELEMENTS
{
  namespace Details
  {
    template <Core::FE::CellType celltype>
    inline static constexpr int num_str = Core::FE::dim<celltype>*(Core::FE::dim<celltype> + 1) / 2;

    /*!
     * @brief Assemble a vector into a matrix row
     *
     * @tparam num_str
     * @param vector (in) : Vector to be assembled into matrix
     * @param data (in/out) : Matrix the vector is assembled into
     * @param row (in) : Matrix row
     */
    template <unsigned num_str>
    void assemble_vector_to_matrix_row(Core::LinAlg::Matrix<num_str, 1> vector,
        Core::LinAlg::SerialDenseMatrix& data, const int row)
    {
      for (unsigned i = 0; i < num_str; ++i) data(row, static_cast<int>(i)) = vector(i);
    }
  }  // namespace Details

  template <typename T>
  inline std::vector<char>& get_stress_data(const T& ele, const Teuchos::ParameterList& params)
  {
    if (ele.IsParamsInterface())
    {
      return *ele.params_interface().stress_data_ptr();
    }
    else
    {
      return *params.get<Teuchos::RCP<std::vector<char>>>("stress");
    }
  }

  template <typename T>
  inline std::vector<char>& get_strain_data(const T& ele, const Teuchos::ParameterList& params)
  {
    if (ele.IsParamsInterface())
    {
      return *ele.params_interface().strain_data_ptr();
    }
    else
    {
      return *params.get<Teuchos::RCP<std::vector<char>>>("strain");
    }
  }

  template <typename T>
  inline Inpar::Solid::StressType get_io_stress_type(
      const T& ele, const Teuchos::ParameterList& params)
  {
    if (ele.IsParamsInterface())
    {
      return ele.params_interface().get_stress_output_type();
    }
    else
    {
      return Core::UTILS::GetAsEnum<Inpar::Solid::StressType>(params, "iostress");
    }
  }

  template <typename T>
  inline Inpar::Solid::StrainType get_io_strain_type(
      const T& ele, const Teuchos::ParameterList& params)
  {
    if (ele.IsParamsInterface())
    {
      return ele.params_interface().get_strain_output_type();
    }
    else
    {
      return Core::UTILS::GetAsEnum<Inpar::Solid::StrainType>(params, "iostrain");
    }
  }

  /*!
   * @brief Convert Green-Lagrange strains to the desired strain type and assemble to a given matrix
   * row in stress-like Voigt notation
   *
   * @tparam celltype : Cell type
   * @param gl_strain (in) : Green-Lagrange strain
   * @param defgrd (in) : Deformation gradient
   * @param strain_type (in) : Strain type, i.e., Green-Lagrange or Euler-Almansi
   * @param data (in/out) : Matrix the strains are assembled into
   * @param row (in) : Matrix row
   */
  template <Core::FE::CellType celltype>
  void assemble_strain_type_to_matrix_row(
      const Core::LinAlg::Matrix<Details::num_str<celltype>, 1>& gl_strain,
      const Core::LinAlg::Matrix<Core::FE::dim<celltype>, Core::FE::dim<celltype>>& defgrd,
      const Inpar::Solid::StrainType strain_type, Core::LinAlg::SerialDenseMatrix& data,
      const int row)
  {
    switch (strain_type)
    {
      case Inpar::Solid::strain_gl:
      {
        Core::LinAlg::Matrix<Details::num_str<celltype>, 1> gl_strain_stress_like;
        Core::LinAlg::Voigt::Strains::to_stress_like(gl_strain, gl_strain_stress_like);
        Details::assemble_vector_to_matrix_row(gl_strain_stress_like, data, row);
        return;
      }
      case Inpar::Solid::strain_ea:
      {
        const Core::LinAlg::Matrix<Details::num_str<celltype>, 1> ea =
            Solid::UTILS::green_lagrange_to_euler_almansi(gl_strain, defgrd);
        Core::LinAlg::Matrix<Details::num_str<celltype>, 1> ea_stress_like;
        Core::LinAlg::Voigt::Strains::to_stress_like(ea, ea_stress_like);
        Details::assemble_vector_to_matrix_row(ea_stress_like, data, row);
        return;
      }
      case Inpar::Solid::strain_log:
      {
        const Core::LinAlg::Matrix<Details::num_str<celltype>, 1> log_strain =
            Solid::UTILS::green_lagrange_to_log_strain(gl_strain);
        Core::LinAlg::Matrix<Details::num_str<celltype>, 1> log_strain_stress_like;
        Core::LinAlg::Voigt::Strains::to_stress_like(log_strain, log_strain_stress_like);
        Details::assemble_vector_to_matrix_row(log_strain_stress_like, data, row);
        return;
      }
      case Inpar::Solid::strain_none:
        return;
      default:
        FOUR_C_THROW("strain type not supported");
        break;
    }
  }

  /*!
   * @brief Convert 2nd Piola-Kirchhoff stresses to the desired stress type and assemble to a given
   * matrix row in stress-like Voigt notation
   *
   * @tparam celltype : Cell type
   * @param defgrd (in) : Deformation gradient
   * @param stress (in) : 2nd Piola-Kirchhoff stress
   * @param stress_type (in) : Stress type, i.e., 2nd Piola-Kirchhoff or Cauchy
   * @param data (in/out) : Matrix the stresses are assembled into
   * @param row (in) : Matrix row
   */
  template <Core::FE::CellType celltype>
  void assemble_stress_type_to_matrix_row(
      const Core::LinAlg::Matrix<Core::FE::dim<celltype>, Core::FE::dim<celltype>>& defgrd,
      const Stress<celltype>& stress, const Inpar::Solid::StressType stress_type,
      Core::LinAlg::SerialDenseMatrix& data, const int row)
  {
    switch (stress_type)
    {
      case Inpar::Solid::stress_2pk:
      {
        Details::assemble_vector_to_matrix_row(stress.pk2_, data, row);
        return;
      }
      case Inpar::Solid::stress_cauchy:
      {
        Core::LinAlg::Matrix<DETAIL::num_str<celltype>, 1> cauchy;
        Solid::UTILS::pk2_to_cauchy(stress.pk2_, defgrd, cauchy);
        Details::assemble_vector_to_matrix_row(cauchy, data, row);
        return;
      }
      case Inpar::Solid::stress_none:

        return;
      default:
        FOUR_C_THROW("stress type not supported");
        break;
    }
  }

  /*!
   * @brief Serialize a matrix by conversion to a vector representation

   * @param matrix (in) : Matrix
   * @param serialized_matrix (in/out) : Serialized matrix
   */
  inline void serialize(
      const Core::LinAlg::SerialDenseMatrix& matrix, std::vector<char>& serialized_matrix)
  {
    Core::Communication::PackBuffer packBuffer;
    Core::Communication::ParObject::add_to_pack(packBuffer, matrix);
    std::copy(packBuffer().begin(), packBuffer().end(), std::back_inserter(serialized_matrix));
  }

  /*!
   * @brief Asks the material for the Gauss Point output quantities and adds the information to
   * the Gauss point output data manager
   *
   * @param num_gp (in) : Number of Gauss Points of the element
   * @param solid_material (in) : Solid material of the element
   * @param gp_data_output_manager (in/out) : Gauss point data output manager
   *                                          (only for new structure time integration)
   */
  inline void ask_and_add_quantities_to_gauss_point_data_output(const int num_gp,
      const Mat::So3Material& solid_material,
      Solid::MODELEVALUATOR::GaussPointDataOutputManager& gp_data_output_manager)
  {
    // Save number of Gauss Points of the element for gauss point data output
    gp_data_output_manager.add_element_number_of_gauss_points(num_gp);

    // holder for output quantity names and their size
    std::unordered_map<std::string, int> quantities_map{};

    // Ask material for the output quantity names and sizes
    solid_material.register_output_data_names(quantities_map);

    // Add quantities to the Gauss point output data manager (if they do not already exist)
    gp_data_output_manager.merge_quantities(quantities_map);
  }

  /*!
   * @brief Collect Gauss Point output data from material and assemble/interpolate depending on
   * output type to element center, Gauss Points, or nodes
   *
   * @tparam celltype : Cell type
   * @param stiffness_matrix_integration (in) : Container holding the integration points
   * @param solid_material (in) : Solid material of the element
   * @param ele (in) : Reference to the element
   * @param gp_data_output_manager (in/out) : Gauss point data output manager
   *                                          (only for new structure time integration)
   */
  template <Core::FE::CellType celltype>
  inline void collect_and_assemble_gauss_point_data_output(
      const Core::FE::GaussIntegration& stiffness_matrix_integration,
      const Mat::So3Material& solid_material, const Core::Elements::Element& ele,
      Solid::MODELEVALUATOR::GaussPointDataOutputManager& gp_data_output_manager)
  {
    // Collection and assembly of gauss point data
    for (const auto& quantity : gp_data_output_manager.get_quantities())
    {
      const std::string& quantity_name = quantity.first;
      const int quantity_size = quantity.second;

      // Step 1: Collect the data for each Gauss point for the material
      Core::LinAlg::SerialDenseMatrix gp_data(
          stiffness_matrix_integration.NumPoints(), quantity_size, true);
      bool data_available = solid_material.EvaluateOutputData(quantity_name, gp_data);

      // Step 2: Assemble data based on output type (elecenter, postprocessed to nodes, Gauss
      // point)
      if (data_available)
      {
        switch (gp_data_output_manager.get_output_type())
        {
          case Inpar::Solid::GaussPointDataOutputType::element_center:
          {
            // compute average of the quantities
            Teuchos::RCP<Epetra_MultiVector> global_data =
                gp_data_output_manager.get_element_center_data().at(quantity_name);
            Core::FE::AssembleAveragedElementValues(*global_data, gp_data, ele);
            break;
          }
          case Inpar::Solid::GaussPointDataOutputType::nodes:
          {
            Teuchos::RCP<Epetra_MultiVector> global_data =
                gp_data_output_manager.get_nodal_data().at(quantity_name);

            Epetra_IntVector& global_nodal_element_count =
                *gp_data_output_manager.get_nodal_data_count().at(quantity_name);

            Core::FE::ExtrapolateGPQuantityToNodesAndAssemble<celltype>(
                ele, gp_data, *global_data, false, stiffness_matrix_integration);
            Discret::ELEMENTS::AssembleNodalElementCount(global_nodal_element_count, ele);
            break;
          }
          case Inpar::Solid::GaussPointDataOutputType::gauss_points:
          {
            std::vector<Teuchos::RCP<Epetra_MultiVector>>& global_data =
                gp_data_output_manager.get_gauss_point_data().at(quantity_name);
            Discret::ELEMENTS::AssembleGaussPointValues(global_data, gp_data, ele);
            break;
          }
          case Inpar::Solid::GaussPointDataOutputType::none:
            FOUR_C_THROW(
                "You specified a Gauss point data output type of none, so you should not end up "
                "here.");
          default:
            FOUR_C_THROW("Unknown Gauss point data output type.");
        }
      }
    }
  }
}  // namespace Discret::ELEMENTS
FOUR_C_NAMESPACE_CLOSE

#endif