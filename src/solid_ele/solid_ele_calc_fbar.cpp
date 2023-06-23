/*----------------------------------------------------------------------*/
/*! \file

\brief main file containing routines for calculation of solid element
       simple displacement based
\level 1
*/
/*----------------------------------------------------------------------*/

#include "solid_ele_calc_lib.H"
#include "solid_ele_calc_fbar.H"
#include <Teuchos_ParameterList.hpp>
#include <memory>
#include <optional>
#include "utils_exceptions.H"
#include "discretization_fem_general_utils_integration.H"
#include "lib_utils.H"
#include "lib_voigt_notation.H"
#include "solid_ele.H"
#include "lib_discret.H"
#include "mat_so3_material.H"
#include "solid_ele_utils.H"
#include "discretization_fem_general_utils_local_connectivity_matrices.H"
#include "fiber_node.H"
#include "fiber_utils.H"
#include "fiber_nodal_fiber_holder.H"
#include "discretization_fem_general_utils_gauss_point_postprocess.H"
#include "discretization_fem_general_utils_gauss_point_extrapolation.H"

#include "structure_new_gauss_point_data_output_manager.H"
#include "so3_element_service.H"
#include "utils_singleton_owner.H"

namespace
{
  /*!
   * @brief Evaluate the fbar factor \f[ \frac{\mathbf{F}_{\mathrm{centroid}}}{\mathbf{F}}^{1/3} \f]
   *
   * @param defgrd_centroid (in) : Deformation gradient evaluated at the element centroid
   * @param defgrd_gp (in) : Deformation gradient evaluated at the Gauss point
   * @return double : Fbar factor
   */
  double EvaluateFbarFactor(const double& defgrd_centroid, const double& defgrd_gp)
  {
    const double fbar_factor = std::pow(defgrd_centroid / defgrd_gp, 1.0 / 3.0);
    return fbar_factor;
  }

  /*!
   * @brief Evaluates the parameter coordinate of at the element centroid
   *
   * Currently this returns a zero-filled matrix, i.e. xi_i = 0. It is valid for HEXes but not
   * necessarily for other element types.
   *
   * @tparam distype : Discretization type
   * @return LINALG::Matrix<DETAIL::nsd<distype>, 1> : Coordinates of the centroid in the parameter
   * space
   */
  template <DRT::Element::DiscretizationType distype,
      std::enable_if_t<DRT::is_hex_v<distype>, int> = 0>
  LINALG::Matrix<DRT::ELEMENTS::DETAIL::nsd<distype>, 1> EvaluateParameterCoordinateCentroid()
  {
    LINALG::Matrix<DRT::ELEMENTS::DETAIL::nsd<distype>, 1> xi;
    for (int d = 0; d < DRT::ELEMENTS::DETAIL::nsd<distype>; ++d) xi(d) = 0;

    return xi;
  }

  /*!
   * @brief Evaluate the strain measures at the element centroid
   *
   * @tparam distype : Discretization type
   * @param nodal_coordinates (in) : Reference and current coordinates of the nodes of the element
   * @return double : Determinant of the deformation gradient at the centroid
   */
  template <DRT::Element::DiscretizationType distype,
      std::enable_if_t<DRT::ELEMENTS::DETAIL::nsd<distype> == 3, int> = 0>
  double EvaluateDeformationGradientDeterminantCentroid(
      DRT::ELEMENTS::NodalCoordinates<distype> nodal_coordinates)
  {
    // set coordinates in parameter space at centroid as zero -> xi = [0; 0; 0]
    LINALG::Matrix<DRT::ELEMENTS::DETAIL::nsd<distype>, 1> xi_centroid =
        EvaluateParameterCoordinateCentroid<distype>();

    // shape functions and derivatives evaluated at element centroid
    const DRT::ELEMENTS::ShapeFunctionsAndDerivatives<distype> shape_functions_centroid =
        DRT::ELEMENTS::EvaluateShapeFunctionsAndDerivs<distype>(xi_centroid);

    // jacobian mapping evaluated at centroid
    const DRT::ELEMENTS::JacobianMapping<distype> jacobian_mapping_centroid =
        DRT::ELEMENTS::EvaluateJacobianMapping(shape_functions_centroid, nodal_coordinates);

    // deformation gradient and strains at centroid of element
    const DRT::ELEMENTS::Strains<distype> strains_centroid =
        DRT::ELEMENTS::EvaluateStrains<distype>(nodal_coordinates, jacobian_mapping_centroid);

    return strains_centroid.defgrd_.Determinant();
  }

  /*!
   * @brief Evaluates the H-Operator used in F-bar of the specified element
   *
   * @tparam distype : Discretization type
   * @param jacobian_mapping (in) : Quantities of the jacobian mapping evaluated at the Gauss point
   * @param jacobian_mapping_centroid (in) : Quantities of the jacobian mapping evaluated at the
   * element centroid
   * @param strains (in) : Strain measures evaluated at the Gauss point
   * @param strains_centroid (in) : Strain measures evaluated at the element centroid
   * @return LINALG::Matrix<DETAIL:numdofperelement<distype>, 1> : H-Operator
   */
  template <DRT::Element::DiscretizationType distype,
      std::enable_if_t<DRT::ELEMENTS::DETAIL::nsd<distype> == 3, int> = 0>
  const LINALG::Matrix<DRT::ELEMENTS::DETAIL::numdofperelement<distype>, 1> EvaluateFbarHOperator(
      const DRT::ELEMENTS::JacobianMapping<distype>& jacobian_mapping,
      const DRT::ELEMENTS::JacobianMapping<distype>& jacobian_mapping_centroid,
      const DRT::ELEMENTS::Strains<distype> strains,
      const DRT::ELEMENTS::Strains<distype> strains_centroid)
  {
    // inverse deformation gradient at centroid
    LINALG::Matrix<DRT::ELEMENTS::DETAIL::nsd<distype>, DRT::ELEMENTS::DETAIL::nsd<distype>>
        invdefgrd_centroid;
    invdefgrd_centroid.Invert(strains_centroid.defgrd_);

    // inverse deformation gradient at gp
    LINALG::Matrix<DRT::ELEMENTS::DETAIL::nsd<distype>, DRT::ELEMENTS::DETAIL::nsd<distype>>
        invdefgrd;
    invdefgrd.Invert(strains.defgrd_);

    LINALG::Matrix<DRT::ELEMENTS::DETAIL::numdofperelement<distype>, 1> Hop(true);
    for (int idof = 0; idof < DRT::ELEMENTS::DETAIL::numdofperelement<distype>; idof++)
    {
      for (int idim = 0; idim < DRT::ELEMENTS::DETAIL::nsd<distype>; idim++)
      {
        Hop(idof) +=
            invdefgrd_centroid(idim, idof % DRT::ELEMENTS::DETAIL::nsd<distype>) *
            jacobian_mapping_centroid.n_xyz_(idim, idof / DRT::ELEMENTS::DETAIL::nsd<distype>);
        Hop(idof) -= invdefgrd(idim, idof % DRT::ELEMENTS::DETAIL::nsd<distype>) *
                     jacobian_mapping.n_xyz_(idim, idof / DRT::ELEMENTS::DETAIL::nsd<distype>);
      }
    }

    return Hop;
  }

  /*!
   * @brief Add fbar stiffness matrix contribution of one Gauss point
   *
   * @tparam distype : Discretization type
   * @param Bop (in) : Strain gradient (B-Operator)
   * @param Hop (in) : H-Operator
   * @param f_bar_factor (in) : f_bar_factor
   * @param integration_fac (in) : Integration factor (Gauss point weight times the determinant of
   * the jacobian)
   * @param strains (in) : Strain measures
   * @param stress_bar (in) : Deviatoric part of stress measures
   * @param stiffness_matrix (in/out) : stiffness matrix where the local contribution is added to
   */
  template <DRT::Element::DiscretizationType distype>
  void AddFbarStiffnessMatrix(const LINALG::Matrix<DRT::ELEMENTS::DETAIL::numstr<distype>,
                                  DRT::ELEMENTS::DETAIL::numdofperelement<distype>>& Bop,
      const LINALG::Matrix<DRT::ELEMENTS::DETAIL::numdofperelement<distype>, 1>& Hop,
      const double f_bar_factor, const double integration_fac,
      const DRT::ELEMENTS::Strains<distype> strains,
      const DRT::ELEMENTS::Stress<distype> stress_bar,
      LINALG::Matrix<DRT::ELEMENTS::DETAIL::numdofperelement<distype>,
          DRT::ELEMENTS::DETAIL::numdofperelement<distype>>& stiffness_matrix)
  {
    LINALG::Matrix<DRT::ELEMENTS::DETAIL::numstr<distype>, 1> rcg_bar_voigt;
    ::UTILS::VOIGT::Strains::MatrixToVector(strains.rcg_, rcg_bar_voigt);

    LINALG::Matrix<DRT::ELEMENTS::DETAIL::numstr<distype>, 1> ccg;
    ccg.MultiplyNN(stress_bar.cmat_, rcg_bar_voigt);

    // auxiliary integrated stress_bar
    LINALG::Matrix<DRT::ELEMENTS::DETAIL::numdofperelement<distype>, 1> bopccg(false);
    bopccg.MultiplyTN(integration_fac * f_bar_factor / 3.0, Bop, ccg);

    LINALG::Matrix<DRT::ELEMENTS::DETAIL::numdofperelement<distype>, 1> bops(false);
    bops.MultiplyTN(-integration_fac / f_bar_factor / 3.0, Bop, stress_bar.pk2_);

    for (int idof = 0; idof < DRT::ELEMENTS::DETAIL::numdofperelement<distype>; idof++)
    {
      for (int jdof = 0; jdof < DRT::ELEMENTS::DETAIL::numdofperelement<distype>; jdof++)
      {
        stiffness_matrix(idof, jdof) += Hop(jdof) * (bops(idof, 0) + bopccg(idof, 0));
      }
    }
  }

  template <DRT::Element::DiscretizationType distype>
  DRT::ELEMENTS::Strains<distype> EvaluateStrainsBar(
      const DRT::ELEMENTS::NodalCoordinates<distype>& nodal_coordinates,
      const DRT::ELEMENTS::JacobianMapping<distype>& jacobian_mapping, double detF_centroid)
  {
    const DRT::ELEMENTS::Strains<distype> strains =
        DRT::ELEMENTS::EvaluateStrains<distype>(nodal_coordinates, jacobian_mapping);

    // factor (detF0/detF)^1/3
    const double fbar_factor = EvaluateFbarFactor(detF_centroid, strains.defgrd_.Determinant());

    // deformation gradient F_bar and resulting strains: F_bar = (detF_0/detF)^1/3 F
    return DRT::ELEMENTS::EvaluateStrains<distype>(
        nodal_coordinates, jacobian_mapping, fbar_factor);
  }
}  // namespace

template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::SolidEleCalcFbar<distype>* DRT::ELEMENTS::SolidEleCalcFbar<distype>::Instance(
    CORE::UTILS::SingletonAction action)
{
  static auto singleton_owner = CORE::UTILS::MakeSingletonOwner(
      []()
      {
        return std::unique_ptr<DRT::ELEMENTS::SolidEleCalcFbar<distype>>(
            new DRT::ELEMENTS::SolidEleCalcFbar<distype>());
      });
  return singleton_owner.Instance(action);
}

template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::SolidEleCalcFbar<distype>::SolidEleCalcFbar()
    : DRT::ELEMENTS::SolidEleCalcInterface::SolidEleCalcInterface(),
      stiffness_matrix_integration_(
          CreateGaussIntegration<distype>(GetGaussRuleStiffnessMatrix<distype>())),
      mass_matrix_integration_(CreateGaussIntegration<distype>(GetGaussRuleMassMatrix<distype>()))
{
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::SolidEleCalcFbar<distype>::EvaluateNonlinearForceStiffnessMass(
    const DRT::Element& ele, MAT::So3Material& solid_material,
    const DRT::Discretization& discretization, const std::vector<int>& lm,
    Teuchos::ParameterList& params, Epetra_SerialDenseVector* force_vector,
    Epetra_SerialDenseMatrix* stiffness_matrix, Epetra_SerialDenseMatrix* mass_matrix)
{
  // Create views to SerialDenseMatrices
  std::optional<LINALG::Matrix<numdofperelement_, numdofperelement_>> stiff = {};
  std::optional<LINALG::Matrix<numdofperelement_, numdofperelement_>> mass = {};
  std::optional<LINALG::Matrix<numdofperelement_, 1>> force = {};
  if (stiffness_matrix != nullptr) stiff.emplace(*stiffness_matrix, true);
  if (mass_matrix != nullptr) mass.emplace(*mass_matrix, true);
  if (force_vector != nullptr) force.emplace(*force_vector, true);

  // get current nodal coordinates
  const NodalCoordinates<distype> nodal_coordinates =
      EvaluateNodalCoordinates<distype>(ele, discretization, lm);

  // TODO: This is a quite unsafe check, whether the same integrations are used
  bool equal_integration_mass_stiffness =
      mass_matrix_integration_.NumPoints() == stiffness_matrix_integration_.NumPoints();

  // initialize element density
  double mean_density = 0.0;

  // set coordinates in parameter space at centroid as zero -> xi = [0; 0; 0]
  LINALG::Matrix<nsd_, 1> xi_centroid = EvaluateParameterCoordinateCentroid<distype>();

  // shape functions and derivatives evaluated at element centroid
  const ShapeFunctionsAndDerivatives<distype> shape_functions_centroid =
      EvaluateShapeFunctionsAndDerivs<distype>(xi_centroid);

  // jacobian mapping evaluated at element centroid
  const JacobianMapping<distype> jacobian_mapping_centroid =
      EvaluateJacobianMapping(shape_functions_centroid, nodal_coordinates);

  // deformation gradient and strains at element centroid
  const Strains<distype> strains_centroid =
      EvaluateStrains<distype>(nodal_coordinates, jacobian_mapping_centroid);

  IterateJacobianMappingAtGaussPoints<distype>(nodal_coordinates, stiffness_matrix_integration_,
      [&](const LINALG::Matrix<DETAIL::nsd<distype>, 1>& xi,
          const ShapeFunctionsAndDerivatives<distype>& shape_functions,
          const JacobianMapping<distype>& jacobian_mapping, double integration_factor, int gp)
      {
        const Strains<distype> strains =
            EvaluateStrains<distype>(nodal_coordinates, jacobian_mapping);

        LINALG::Matrix<numstr_, numdofperelement_> Bop =
            EvaluateStrainGradient(jacobian_mapping, strains);

        // factor (detF0/detF)^1/3
        const double fbar_factor = EvaluateFbarFactor(
            strains_centroid.defgrd_.Determinant(), strains.defgrd_.Determinant());

        const LINALG::Matrix<numdofperelement_, 1> Hop = EvaluateFbarHOperator(
            jacobian_mapping, jacobian_mapping_centroid, strains, strains_centroid);

        // deformation gradient F_bar and resulting strains: F_bar = (detF_0/detF)^1/3 F
        const Strains<distype> strains_bar =
            EvaluateStrains<distype>(nodal_coordinates, jacobian_mapping, fbar_factor);

        // stress stress_bar evaluated using strains_bar
        const Stress<distype> stress_bar =
            EvaluateMaterialStress(solid_material, strains_bar, params, gp, ele.Id());

        // evaluate internal force vector using only the deviatoric component
        if (force.has_value())
          AddInternalForceVector(Bop, stress_bar, integration_factor / fbar_factor, *force);

        if (stiff.has_value())
        {
          // evaluate elastic stiffness matrix using only the deviatoric component Fbar
          AddElasticStiffnessMatrix(Bop, stress_bar, integration_factor * fbar_factor, *stiff);

          // evaluate geometric stiffness matrix using only the deviatoric component Fbar
          AddGeometricStiffnessMatrix(
              jacobian_mapping.n_xyz_, stress_bar, integration_factor / fbar_factor, *stiff);

          // additional stiffness matrix needed for fbar method
          AddFbarStiffnessMatrix(
              Bop, Hop, fbar_factor, integration_factor, strains, stress_bar, *stiff);
        }

        if (mass.has_value())
        {
          if (equal_integration_mass_stiffness)
          {
            AddMassMatrix(shape_functions, integration_factor, solid_material.Density(gp), *mass);
          }
          else
          {
            mean_density += solid_material.Density(gp) / stiffness_matrix_integration_.NumPoints();
          }
        }
      });

  if (mass.has_value() && !equal_integration_mass_stiffness)
  {
    // integrate mass matrix
    dsassert(mean_density > 0, "It looks like the density is 0.0");
    IterateJacobianMappingAtGaussPoints<distype>(nodal_coordinates, mass_matrix_integration_,
        [&](const LINALG::Matrix<DETAIL::nsd<distype>, 1>& xi,
            const ShapeFunctionsAndDerivatives<distype>& shape_functions,
            const JacobianMapping<distype>& jacobian_mapping, double integration_factor, int gp)
        { AddMassMatrix(shape_functions, integration_factor, mean_density, *mass); });
  }
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::SolidEleCalcFbar<distype>::Recover(const DRT::Element& ele,
    const DRT::Discretization& discretization, const std::vector<int>& lm,
    Teuchos::ParameterList& params)
{
  // nothing to do for a standard element
  // except...
  // TODO: We need to recover history information of materials!
  // which was also not implemented in the old elements
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::SolidEleCalcFbar<distype>::Update(const DRT::Element& ele,
    MAT::So3Material& solid_material, const DRT::Discretization& discretization,
    const std::vector<int>& lm, Teuchos::ParameterList& params)
{
  const NodalCoordinates<distype> nodal_coordinates =
      EvaluateNodalCoordinates<distype>(ele, discretization, lm);

  // deformation gradient and strains at centroid of element
  auto detF_centroid = EvaluateDeformationGradientDeterminantCentroid<distype>(nodal_coordinates);

  // Loop over all Gauss points
  IterateJacobianMappingAtGaussPoints<distype>(nodal_coordinates, stiffness_matrix_integration_,
      [&](const LINALG::Matrix<DETAIL::nsd<distype>, 1>& xi,
          const ShapeFunctionsAndDerivatives<distype>& shape_functions,
          const JacobianMapping<distype>& jacobian_mapping, double integration_factor, int gp)
      {
        const Strains<distype> strains_bar =
            EvaluateStrainsBar(nodal_coordinates, jacobian_mapping, detF_centroid);

        solid_material.Update(strains_bar.defgrd_, gp, params, ele.Id());
      });

  solid_material.Update();
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::SolidEleCalcFbar<distype>::CalculateStress(const DRT::Element& ele,
    MAT::So3Material& solid_material, const StressIO& stressIO, const StrainIO& strainIO,
    const DRT::Discretization& discretization, const std::vector<int>& lm,
    Teuchos::ParameterList& params)
{
  // TODO: If we get rid of post_*, we don't need this here anymore. We could directly use
  // InitializeGaussPointDataOutput and EvaluateGaussPointDataOutput and write the stresses there.
  if (discretization.Comm().MyPID() != ele.Owner()) return;

  std::vector<char>& serialized_stress_data = stressIO.mutable_data;
  std::vector<char>& serialized_strain_data = strainIO.mutable_data;
  Epetra_SerialDenseMatrix stress_data(stiffness_matrix_integration_.NumPoints(), numstr_);
  Epetra_SerialDenseMatrix strain_data(stiffness_matrix_integration_.NumPoints(), numstr_);

  const NodalCoordinates<distype> nodal_coordinates =
      EvaluateNodalCoordinates<distype>(ele, discretization, lm);

  // deformation gradient and strains at centroid of element
  auto detF_centroid = EvaluateDeformationGradientDeterminantCentroid<distype>(nodal_coordinates);

  // Loop over all Gauss points
  IterateJacobianMappingAtGaussPoints<distype>(nodal_coordinates, stiffness_matrix_integration_,
      [&](const LINALG::Matrix<DETAIL::nsd<distype>, 1>& xi,
          const ShapeFunctionsAndDerivatives<distype>& shape_functions,
          const JacobianMapping<distype>& jacobian_mapping, double integration_factor, int gp)
      {
        const Strains<distype> strains_bar =
            EvaluateStrainsBar(nodal_coordinates, jacobian_mapping, detF_centroid);

        const Stress<distype> stress_bar =
            EvaluateMaterialStress(solid_material, strains_bar, params, gp, ele.Id());

        AssembleStrainTypeToMatrixRow(strains_bar, strainIO.type, strain_data, gp);
        AssembleStressTypeToMatrixRow(strains_bar, stress_bar, stressIO.type, stress_data, gp);
      });

  Serialize(stress_data, serialized_stress_data);
  Serialize(strain_data, serialized_strain_data);
}

template <DRT::Element::DiscretizationType distype>
double DRT::ELEMENTS::SolidEleCalcFbar<distype>::CalculateInternalEnergy(const DRT::Element& ele,
    MAT::So3Material& solid_material, const DRT::Discretization& discretization,
    const std::vector<int>& lm, Teuchos::ParameterList& params)
{
  double intenergy = 0.0;
  const NodalCoordinates<distype> nodal_coordinates =
      EvaluateNodalCoordinates<distype>(ele, discretization, lm);

  // deformation gradient and strains at centroid of element
  auto detF_centroid = EvaluateDeformationGradientDeterminantCentroid<distype>(nodal_coordinates);

  IterateJacobianMappingAtGaussPoints<distype>(nodal_coordinates, stiffness_matrix_integration_,
      [&](const LINALG::Matrix<DETAIL::nsd<distype>, 1>& xi,
          const ShapeFunctionsAndDerivatives<distype>& shape_functions,
          const JacobianMapping<distype>& jacobian_mapping, double integration_factor, int gp)
      {
        const Strains<distype> strains_bar =
            EvaluateStrainsBar(nodal_coordinates, jacobian_mapping, detF_centroid);

        double psi = 0.0;
        solid_material.StrainEnergy(strains_bar.gl_strain_, psi, gp, ele.Id());

        intenergy += psi * integration_factor;
      });

  return intenergy;
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::SolidEleCalcFbar<distype>::Setup(
    MAT::So3Material& solid_material, DRT::INPUT::LineDefinition* linedef)
{
  solid_material.Setup(stiffness_matrix_integration_.NumPoints(), linedef);
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::SolidEleCalcFbar<distype>::MaterialPostSetup(
    const DRT::Element& ele, MAT::So3Material& solid_material)
{
  Teuchos::ParameterList params{};
  if (DRT::FIBER::UTILS::HaveNodalFibers<distype>(ele.Nodes()))
  {
    // This element has fiber nodes.
    // Interpolate fibers to the Gauss points and pass them to the material

    // Get shape functions
    const static std::vector<LINALG::Matrix<nen_, 1>> shapefcts = std::invoke(
        [&]
        {
          std::vector<LINALG::Matrix<nen_, 1>> shapefcns(stiffness_matrix_integration_.NumPoints());
          for (int gp = 0; gp < stiffness_matrix_integration_.NumPoints(); ++gp)
          {
            LINALG::Matrix<nsd_, 1> xi(stiffness_matrix_integration_.Point(gp), true);
            CORE::DRT::UTILS::shape_function<distype>(xi, shapefcns[gp]);
          }
          return shapefcns;
        });

    // add fibers to the ParameterList
    DRT::FIBER::NodalFiberHolder fiberHolder;

    // Do the interpolation
    DRT::FIBER::UTILS::ProjectFibersToGaussPoints<distype>(ele.Nodes(), shapefcts, fiberHolder);

    params.set("fiberholder", fiberHolder);
  }

  // Call PostSetup of material
  solid_material.PostSetup(params, ele.Id());
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::SolidEleCalcFbar<distype>::InitializeGaussPointDataOutput(
    const DRT::Element& ele, const MAT::So3Material& solid_material,
    STR::MODELEVALUATOR::GaussPointDataOutputManager& gp_data_output_manager) const
{
  dsassert(ele.IsParamsInterface(),
      "This action type should only be called from the new time integration framework!");

  // Save number of Gauss of the element for gauss point data output
  gp_data_output_manager.AddElementNumberOfGaussPoints(stiffness_matrix_integration_.NumPoints());

  // holder for output quantity names and their size
  std::unordered_map<std::string, int> quantities_map{};

  // Ask material for the output quantity names and sizes
  solid_material.RegisterVtkOutputDataNames(quantities_map);

  // Add quantities to the Gauss point output data manager (if they do not already exist)
  gp_data_output_manager.MergeQuantities(quantities_map);
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::SolidEleCalcFbar<distype>::EvaluateGaussPointDataOutput(const DRT::Element& ele,
    const MAT::So3Material& solid_material,
    STR::MODELEVALUATOR::GaussPointDataOutputManager& gp_data_output_manager) const
{
  dsassert(ele.IsParamsInterface(),
      "This action type should only be called from the new time integration framework!");

  // Collection and assembly of gauss point data
  for (const auto& quantity : gp_data_output_manager.GetQuantities())
  {
    const std::string& quantity_name = quantity.first;
    const int quantity_size = quantity.second;

    // Step 1: Collect the data for each Gauss point for the material
    LINALG::SerialDenseMatrix gp_data(
        stiffness_matrix_integration_.NumPoints(), quantity_size, true);
    bool data_available = solid_material.EvaluateVtkOutputData(quantity_name, gp_data);

    // Step 3: Assemble data based on output type (elecenter, postprocessed to nodes, Gauss
    // point)
    if (data_available)
    {
      switch (gp_data_output_manager.GetOutputType())
      {
        case INPAR::STR::GaussPointDataOutputType::element_center:
        {
          // compute average of the quantities
          Teuchos::RCP<Epetra_MultiVector> global_data =
              gp_data_output_manager.GetMutableElementCenterData().at(quantity_name);
          CORE::DRT::ELEMENTS::AssembleAveragedElementValues(*global_data, gp_data, ele);
          break;
        }
        case INPAR::STR::GaussPointDataOutputType::nodes:
        {
          Teuchos::RCP<Epetra_MultiVector> global_data =
              gp_data_output_manager.GetMutableNodalData().at(quantity_name);

          Epetra_IntVector& global_nodal_element_count =
              *gp_data_output_manager.GetMutableNodalDataCount().at(quantity_name);

          CORE::DRT::UTILS::ExtrapolateGPQuantityToNodesAndAssemble<distype>(
              ele, gp_data, *global_data, false, stiffness_matrix_integration_);
          DRT::ELEMENTS::AssembleNodalElementCount(global_nodal_element_count, ele);
          break;
        }
        case INPAR::STR::GaussPointDataOutputType::gauss_points:
        {
          std::vector<Teuchos::RCP<Epetra_MultiVector>>& global_data =
              gp_data_output_manager.GetMutableGaussPointData().at(quantity_name);
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

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::SolidEleCalcFbar<distype>::ResetAll(
    const DRT::Element& ele, MAT::So3Material& solid_material)
{
  solid_material.ResetAll(stiffness_matrix_integration_.NumPoints());
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::SolidEleCalcFbar<distype>::ResetToLastConverged(
    const DRT::Element& ele, MAT::So3Material& solid_material)
{
  solid_material.ResetStep();
}

// template classes
template class DRT::ELEMENTS::SolidEleCalcFbar<DRT::Element::hex8>;