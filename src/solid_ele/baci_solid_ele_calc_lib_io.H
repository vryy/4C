/*! \file

\brief A library of free functions for a default solid element

\level 1
*/

#ifndef BACI_SOLID_ELE_CALC_LIB_IO_H
#define BACI_SOLID_ELE_CALC_LIB_IO_H

#include "baci_config.H"

#include "baci_discretization_fem_general_cell_type.H"
#include "baci_discretization_fem_general_cell_type_traits.H"
#include "baci_discretization_fem_general_utils_gauss_point_extrapolation.H"
#include "baci_inpar_structure.H"
#include "baci_linalg_fixedsizematrix.H"
#include "baci_linalg_fixedsizematrix_voigt_notation.H"
#include "baci_so3_element_service.H"
#include "baci_solid_ele_calc_lib.H"
#include "baci_solid_ele_utils.H"
#include "baci_structure_new_gauss_point_data_output_manager.H"

#include <Teuchos_ParameterList.hpp>

BACI_NAMESPACE_OPEN

namespace DRT::ELEMENTS
{
  namespace DETAILS
  {
    template <CORE::FE::CellType celltype>
    inline static constexpr int num_str = CORE::FE::dim<celltype>*(CORE::FE::dim<celltype> + 1) / 2;

    /*!
     * @brief Assemble a vector into a matrix row
     *
     * @tparam num_str
     * @param vector (in) : Vector to be assembled into matrix
     * @param data (in/out) : Matrix the vector is assembled into
     * @param row (in) : Matrix row
     */
    template <unsigned num_str>
    void AssembleVectorToMatrixRow(CORE::LINALG::Matrix<num_str, 1> vector,
        CORE::LINALG::SerialDenseMatrix& data, const int row)
    {
      for (unsigned i = 0; i < num_str; ++i) data(row, static_cast<int>(i)) = vector(i);
    }
  }  // namespace DETAILS

  template <typename T>
  inline std::vector<char>& GetStressData(const T& ele, const Teuchos::ParameterList& params)
  {
    if (ele.IsParamsInterface())
    {
      return *ele.ParamsInterface().StressDataPtr();
    }
    else
    {
      return *params.get<Teuchos::RCP<std::vector<char>>>("stress");
    }
  }

  template <typename T>
  inline std::vector<char>& GetStrainData(const T& ele, const Teuchos::ParameterList& params)
  {
    if (ele.IsParamsInterface())
    {
      return *ele.ParamsInterface().StrainDataPtr();
    }
    else
    {
      return *params.get<Teuchos::RCP<std::vector<char>>>("strain");
    }
  }

  template <typename T>
  inline INPAR::STR::StressType GetIOStressType(const T& ele, const Teuchos::ParameterList& params)
  {
    if (ele.IsParamsInterface())
    {
      return ele.ParamsInterface().GetStressOutputType();
    }
    else
    {
      return DRT::INPUT::get<INPAR::STR::StressType>(params, "iostress");
    }
  }

  template <typename T>
  inline INPAR::STR::StrainType GetIOStrainType(const T& ele, const Teuchos::ParameterList& params)
  {
    if (ele.IsParamsInterface())
    {
      return ele.ParamsInterface().GetStrainOutputType();
    }
    else
    {
      return DRT::INPUT::get<INPAR::STR::StrainType>(params, "iostrain");
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
  template <CORE::FE::CellType celltype>
  void AssembleStrainTypeToMatrixRow(
      const CORE::LINALG::Matrix<DETAILS::num_str<celltype>, 1>& gl_strain,
      const CORE::LINALG::Matrix<CORE::FE::dim<celltype>, CORE::FE::dim<celltype>>& defgrd,
      const INPAR::STR::StrainType strain_type, CORE::LINALG::SerialDenseMatrix& data,
      const int row)
  {
    switch (strain_type)
    {
      case INPAR::STR::strain_gl:
      {
        CORE::LINALG::Matrix<DETAILS::num_str<celltype>, 1> gl_strain_stress_like;
        CORE::LINALG::VOIGT::Strains::ToStressLike(gl_strain, gl_strain_stress_like);
        DETAILS::AssembleVectorToMatrixRow(gl_strain_stress_like, data, row);
        return;
      }
      case INPAR::STR::strain_ea:
      {
        const CORE::LINALG::Matrix<DETAILS::num_str<celltype>, 1> ea =
            STR::UTILS::GreenLagrangeToEulerAlmansi(gl_strain, defgrd);
        CORE::LINALG::Matrix<DETAILS::num_str<celltype>, 1> ea_stress_like;
        CORE::LINALG::VOIGT::Strains::ToStressLike(ea, ea_stress_like);
        DETAILS::AssembleVectorToMatrixRow(ea_stress_like, data, row);
        return;
      }
      case INPAR::STR::strain_none:
        return;
      default:
        dserror("strain type not supported");
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
  template <CORE::FE::CellType celltype>
  void AssembleStressTypeToMatrixRow(
      const CORE::LINALG::Matrix<CORE::FE::dim<celltype>, CORE::FE::dim<celltype>>& defgrd,
      const Stress<celltype>& stress, const INPAR::STR::StressType stress_type,
      CORE::LINALG::SerialDenseMatrix& data, const int row)
  {
    switch (stress_type)
    {
      case INPAR::STR::stress_2pk:
      {
        DETAILS::AssembleVectorToMatrixRow(stress.pk2_, data, row);
        return;
      }
      case INPAR::STR::stress_cauchy:
      {
        CORE::LINALG::Matrix<DETAIL::num_str<celltype>, 1> cauchy;
        STR::UTILS::Pk2ToCauchy(stress.pk2_, defgrd, cauchy);
        DETAILS::AssembleVectorToMatrixRow(cauchy, data, row);
        return;
      }
      case INPAR::STR::stress_none:

        return;
      default:
        dserror("stress type not supported");
        break;
    }
  }

  /*!
   * @brief Serialize a matrix by conversion to a vector representation

   * @param matrix (in) : Matrix
   * @param serialized_matrix (in/out) : Serialized matrix
   */
  inline void Serialize(
      const CORE::LINALG::SerialDenseMatrix& matrix, std::vector<char>& serialized_matrix)
  {
    CORE::COMM::PackBuffer packBuffer;
    CORE::COMM::ParObject::AddtoPack(packBuffer, matrix);
    packBuffer.StartPacking();
    CORE::COMM::ParObject::AddtoPack(packBuffer, matrix);
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
  inline void AskAndAddQuantitiesToGaussPointDataOutput(const int num_gp,
      const MAT::So3Material& solid_material,
      STR::MODELEVALUATOR::GaussPointDataOutputManager& gp_data_output_manager)
  {
    // Save number of Gauss Points of the element for gauss point data output
    gp_data_output_manager.AddElementNumberOfGaussPoints(num_gp);

    // holder for output quantity names and their size
    std::unordered_map<std::string, int> quantities_map{};

    // Ask material for the output quantity names and sizes
    solid_material.RegisterOutputDataNames(quantities_map);

    // Add quantities to the Gauss point output data manager (if they do not already exist)
    gp_data_output_manager.MergeQuantities(quantities_map);
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
  template <CORE::FE::CellType celltype>
  inline void CollectAndAssembleGaussPointDataOutput(
      const CORE::FE::GaussIntegration& stiffness_matrix_integration,
      const MAT::So3Material& solid_material, const DRT::Element& ele,
      STR::MODELEVALUATOR::GaussPointDataOutputManager& gp_data_output_manager)
  {
    // Collection and assembly of gauss point data
    for (const auto& quantity : gp_data_output_manager.GetQuantities())
    {
      const std::string& quantity_name = quantity.first;
      const int quantity_size = quantity.second;

      // Step 1: Collect the data for each Gauss point for the material
      CORE::LINALG::SerialDenseMatrix gp_data(
          stiffness_matrix_integration.NumPoints(), quantity_size, true);
      bool data_available = solid_material.EvaluateOutputData(quantity_name, gp_data);

      // Step 2: Assemble data based on output type (elecenter, postprocessed to nodes, Gauss
      // point)
      if (data_available)
      {
        switch (gp_data_output_manager.GetOutputType())
        {
          case INPAR::STR::GaussPointDataOutputType::element_center:
          {
            // compute average of the quantities
            Teuchos::RCP<Epetra_MultiVector> global_data =
                gp_data_output_manager.GetElementCenterData().at(quantity_name);
            CORE::FE::AssembleAveragedElementValues(*global_data, gp_data, ele);
            break;
          }
          case INPAR::STR::GaussPointDataOutputType::nodes:
          {
            Teuchos::RCP<Epetra_MultiVector> global_data =
                gp_data_output_manager.GetNodalData().at(quantity_name);

            Epetra_IntVector& global_nodal_element_count =
                *gp_data_output_manager.GetNodalDataCount().at(quantity_name);

            CORE::FE::ExtrapolateGPQuantityToNodesAndAssemble<celltype>(
                ele, gp_data, *global_data, false, stiffness_matrix_integration);
            DRT::ELEMENTS::AssembleNodalElementCount(global_nodal_element_count, ele);
            break;
          }
          case INPAR::STR::GaussPointDataOutputType::gauss_points:
          {
            std::vector<Teuchos::RCP<Epetra_MultiVector>>& global_data =
                gp_data_output_manager.GetGaussPointData().at(quantity_name);
            DRT::ELEMENTS::AssembleGaussPointValues(global_data, gp_data, ele);
            break;
          }
          case INPAR::STR::GaussPointDataOutputType::none:
            dserror(
                "You specified a Gauss point data output type of none, so you should not end up "
                "here.");
          default:
            dserror("Unknown Gauss point data output type.");
        }
      }
    }
  }
}  // namespace DRT::ELEMENTS
BACI_NAMESPACE_CLOSE

#endif  // SOLID_ELE_CALC_LIB_IO_H