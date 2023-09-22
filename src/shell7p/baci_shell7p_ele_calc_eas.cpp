/*! \file

\brief Implementation of routines for calculation of shell element with EAS element technology

\level 3
*/

#include "baci_shell7p_ele_calc_eas.H"

#include "baci_lib_discret.H"
#include "baci_lib_voigt_notation.H"
#include "baci_linalg_serialdensematrix.H"
#include "baci_linalg_serialdensevector.H"
#include "baci_linalg_utils_densematrix_inverse.H"
#include "baci_mat_so3_material.H"
#include "baci_shell7p_ele.H"
#include "baci_shell7p_ele_calc_eas_utils.H"
#include "baci_utils_exceptions.H"

#include <Teuchos_ParameterList.hpp>

namespace
{

  /*!
   * @brief Evaluates the enhanced strains scalar increment alpha
   *
   * @tparam distype : Discretization type
   * @param old_eas_data (in) : EAS iteration data
   * @param neas (in) : Number of EAS parameter
   * @param residual (in) : Residual displacement increment
   * @param alpha_inc (out) : Enhanced strains scalar increment
   */
  template <DRT::Element::DiscretizationType distype>
  void EvaluateAlphaIncrement(DRT::ELEMENTS::ShellEASIterationData& old_eas_data, const int& neas,
      const std::vector<double>& residual, CORE::LINALG::SerialDenseMatrix& delta_alpha)
  {
    // we need the (residual) displacement at the previous step
    CORE::LINALG::SerialDenseVector disp_inc(
        DRT::ELEMENTS::SHELL::DETAIL::numdofperelement<distype>);
    for (int i = 0; i < DRT::ELEMENTS::SHELL::DETAIL::numdofperelement<distype>; ++i)
      disp_inc(i) = residual[i];

    CORE::LINALG::SerialDenseMatrix eashelp(neas, 1);
    // make multiplication eashelp = - old L^T * disp_incr[kstep]
    eashelp.multiply(
        Teuchos::NO_TRANS, Teuchos::NO_TRANS, 1.0, old_eas_data.transL_, disp_inc, 0.0);
    // add old RTilde to eashelp
    eashelp += old_eas_data.RTilde_;
    // make multiplication alpha_inc = - old invDTilde * eashelp
    delta_alpha.multiply(
        Teuchos::NO_TRANS, Teuchos::NO_TRANS, -1.0, old_eas_data.invDTilde_, eashelp, 0.0);
  }  // end of EvaluateAlphaIncrement

  /*!
   * @brief Integrate the EAS data
   *
   *  Function needs to be called for every gaussian point
   *
   * @tparam distype : Discretization type
   * @param stress_enh (in) : An object holding the enhanced stress resultants
   * @param M (in) : EAS shapefunctions
   * @param Bop (in) : B-operator matrix
   * @param eas_data (in/out) : An object holding the EAS data
   * @param integration_factor (in) : Integration factor
   */
  template <DRT::Element::DiscretizationType distype>
  void IntegrateEAS(const DRT::ELEMENTS::SHELL::StressEnhanced& stress_enh,
      const CORE::LINALG::SerialDenseMatrix& M, const CORE::LINALG::SerialDenseMatrix& Bop,
      DRT::ELEMENTS::ShellEASIterationData& eas_data, const double& integration_factor,
      const int& neas)
  {
    // integrate D_Tilde += M^T * D * M  * detJ * w(gp)
    // IMPORTANT: here we save D_Tilde in invDTilde_, since after the loop over all Gaussian points,
    // we invert the matrix. At this point, this is still D_Tilde and NOT invD_Tilde
    CORE::LINALG::SerialDenseMatrix MTDmat(
        neas, DRT::ELEMENTS::SHELL::DETAIL::num_internal_variables);
    MTDmat.multiply(Teuchos::TRANS, Teuchos::NO_TRANS, 1.0, M, stress_enh.dmat_, 0.0);

    eas_data.invDTilde_.multiply(
        Teuchos::NO_TRANS, Teuchos::NO_TRANS, integration_factor, MTDmat, M, 1.);

    //  integrate transL (L^T) += M^T * D * B * detJ * w(gp)
    eas_data.transL_.multiply(
        Teuchos::NO_TRANS, Teuchos::NO_TRANS, integration_factor, MTDmat, Bop, 1.);

    //  integrate Rtilde (R_Tilde) : Rtilde  += M^T * stress_r * detJ * w(gp)
    eas_data.RTilde_.multiply(
        Teuchos::TRANS, Teuchos::NO_TRANS, integration_factor, M, stress_enh.stress_, 1.);
  }
}  // namespace


template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::Shell7pEleCalcEas<distype>::Shell7pEleCalcEas()
    : DRT::ELEMENTS::Shell7pEleCalcInterface::Shell7pEleCalcInterface(),
      intpoints_midsurface_(
          SHELL::CreateGaussIntegrationPoints<distype>(SHELL::GetGaussRule<distype>()))
{
  old_step_length_ = 0.0;
  cur_thickness_.resize(intpoints_midsurface_.NumPoints(), shell_data_.thickness);
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Shell7pEleCalcEas<distype>::Setup(DRT::Element& ele,
    MAT::So3Material& solid_material, DRT::INPUT::LineDefinition* linedef,
    const STR::ELEMENTS::ShellLockingTypes& locking_types,
    const STR::ELEMENTS::ShellData& shell_data)
{
  shell_data_ = shell_data;
  cur_thickness_.resize(intpoints_midsurface_.NumPoints(), shell_data_.thickness);
  locking_types_ = locking_types;

  // init sizes of EAS data for integration
  eas_iteration_data_.alpha_.shape(locking_types.total, 1);
  eas_iteration_data_.RTilde_.shape(locking_types_.total, 1);
  eas_iteration_data_.invDTilde_.shape(locking_types_.total, locking_types_.total);
  eas_iteration_data_.transL_.shape(locking_types_.total, SHELL::DETAIL::numdofperelement<distype>);

  //  set up of materials with GP data (e.g., history variables)
  solid_material.Setup(intpoints_midsurface_.NumPoints(), linedef);
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Shell7pEleCalcEas<distype>::Pack(DRT::PackBuffer& data) const
{
  DRT::ELEMENTS::Shell7p::AddtoPack(data, shell_data_.sdc);
  DRT::ELEMENTS::Shell7p::AddtoPack(data, shell_data_.thickness);
  DRT::ELEMENTS::Shell7p::AddtoPack(data, shell_data_.num_ans);

  DRT::ELEMENTS::Shell7p::AddtoPack(data, eas_iteration_data_.alpha_);
  DRT::ELEMENTS::Shell7p::AddtoPack(data, eas_iteration_data_.RTilde_);
  DRT::ELEMENTS::Shell7p::AddtoPack(data, eas_iteration_data_.invDTilde_);
  DRT::ELEMENTS::Shell7p::AddtoPack(data, eas_iteration_data_.transL_);

  // number of total EAS parameters
  DRT::ELEMENTS::Shell7p::AddtoPack(data, locking_types_.membrane);
  DRT::ELEMENTS::Shell7p::AddtoPack(data, locking_types_.bending);
  DRT::ELEMENTS::Shell7p::AddtoPack(data, locking_types_.thickness);
  DRT::ELEMENTS::Shell7p::AddtoPack(data, locking_types_.transverse_shear_strain_const);
  DRT::ELEMENTS::Shell7p::AddtoPack(data, locking_types_.transverse_shear_strain_lin);
  DRT::ELEMENTS::Shell7p::AddtoPack(data, locking_types_.total);

  DRT::ELEMENTS::Shell7p::AddtoPack(data, old_step_length_);
  DRT::ELEMENTS::Shell7p::AddtoPack(data, cur_thickness_);
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Shell7pEleCalcEas<distype>::Unpack(
    std::vector<char>::size_type& position, const std::vector<char>& data)
{
  DRT::ParObject::ExtractfromPack(position, data, shell_data_.sdc);
  DRT::ParObject::ExtractfromPack(position, data, shell_data_.thickness);
  DRT::ParObject::ExtractfromPack(position, data, shell_data_.num_ans);

  DRT::ParObject::ExtractfromPack(position, data, eas_iteration_data_.alpha_);
  DRT::ParObject::ExtractfromPack(position, data, eas_iteration_data_.RTilde_);
  DRT::ParObject::ExtractfromPack(position, data, eas_iteration_data_.invDTilde_);
  DRT::ParObject::ExtractfromPack(position, data, eas_iteration_data_.transL_);
  // number of total EAS parameters
  DRT::ParObject::ExtractfromPack(position, data, locking_types_.membrane);
  DRT::ParObject::ExtractfromPack(position, data, locking_types_.bending);
  DRT::ParObject::ExtractfromPack(position, data, locking_types_.thickness);
  DRT::ParObject::ExtractfromPack(position, data, locking_types_.transverse_shear_strain_const);
  DRT::ParObject::ExtractfromPack(position, data, locking_types_.transverse_shear_strain_lin);
  DRT::ParObject::ExtractfromPack(position, data, locking_types_.total);

  DRT::ParObject::ExtractfromPack(position, data, old_step_length_);
  DRT::ParObject::ExtractfromPack(position, data, cur_thickness_);
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Shell7pEleCalcEas<distype>::MaterialPostSetup(
    DRT::Element& ele, MAT::So3Material& solid_material)
{
  // element/nodal wise defined data
  Teuchos::ParameterList params{};
  // Call PostSetup of material
  solid_material.PostSetup(params, ele.Id());
}


template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Shell7pEleCalcEas<distype>::ResetToLastConverged(
    DRT::Element& ele, MAT::So3Material& solid_material)
{
  solid_material.ResetStep();
}

template <DRT::Element::DiscretizationType distype>
double DRT::ELEMENTS::Shell7pEleCalcEas<distype>::CalculateInternalEnergy(DRT::Element& ele,
    MAT::So3Material& solid_material, const DRT::Discretization& discretization,
    const CORE::LINALG::SerialDenseMatrix& nodal_directors, const std::vector<int>& dof_index_array,
    Teuchos::ParameterList& params)
{
  double intenergy = 0.0;

  Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
  Teuchos::RCP<const Epetra_Vector> res = discretization.GetState("residual displacement");
  std::vector<double> displacement(dof_index_array.size());
  DRT::UTILS::ExtractMyValues(*disp, displacement, dof_index_array);
  std::vector<double> residual(dof_index_array.size());
  DRT::UTILS::ExtractMyValues(*res, residual, dof_index_array);

  // init scale factor for scaled director approach (SDC)
  const double condfac = shell_data_.sdc;

  // get nodal coordinates
  SHELL::NodalCoordinates<distype> nodal_coordinates = SHELL::EvaluateNodalCoordinates<distype>(
      ele.Nodes(), displacement, shell_data_.thickness, nodal_directors, condfac);

  // init gauss point in thickness direction that will be modified via SDC
  double zeta = 0.0;

  // Enhanced Assumed Strain (EAS) Technology: declare, initialize, set up, and alpha history
  CORE::LINALG::SerialDenseMatrix delta_alpha(locking_types_.total, 1);

  // EAS Update of alphas: the current alphas are (re-)evaluated out of DTilde and L^T of previous
  // step to avoid additional element call
  if (not ele.IsParamsInterface())
  {
    // compute the EAS increment delta_alpha
    EvaluateAlphaIncrement<distype>(
        eas_iteration_data_, locking_types_.total, residual, delta_alpha);
    // update alpha += 1.0 * delta_alpha
    eas_iteration_data_.alpha_ += delta_alpha;
  }

  // Assumed Natural Strains (ANS) Technology to remedy transverse shear strain locking
  std::vector<double> shape_functions_ans(true);

  // for a_13 and a_23 each
  const int total_ansq = 2 * shell_data_.num_ans;
  std::vector<SHELL::ShapefunctionsAndDerivatives<distype>> shapefunctions_collocation(total_ansq);
  std::vector<SHELL::BasisVectorsAndMetrics<distype>> metrics_collocation_reference(total_ansq);
  std::vector<SHELL::BasisVectorsAndMetrics<distype>> metrics_collocation_current(total_ansq);

  if (shell_data_.num_ans > 0)
  {
    SHELL::SetupANS(shapefunctions_collocation, metrics_collocation_reference,
        metrics_collocation_current, nodal_coordinates, total_ansq);
  }

  // metric of element at centroid point (for EAS)
  const std::array<double, 2> centroid_point = {0.0, 0.0};
  SHELL::ShapefunctionsAndDerivatives<distype> shapefunctions_centroid =
      SHELL::EvaluateShapefunctionsAndDerivs<distype>(centroid_point);
  SHELL::BasisVectorsAndMetrics<distype> metrics_centroid_reference;
  SHELL::BasisVectorsAndMetrics<distype> metrics_centroid_current;

  SHELL::EvaluateMetrics(shapefunctions_centroid, metrics_centroid_reference,
      metrics_centroid_current, nodal_coordinates, 0.0);

  // init metric tensor and basis vectors of element mid-surface
  SHELL::BasisVectorsAndMetrics<distype> a_reference;
  SHELL::BasisVectorsAndMetrics<distype> a_current;

  // init metric tensor and basis vectors of element shell body
  SHELL::BasisVectorsAndMetrics<distype> g_reference;
  SHELL::BasisVectorsAndMetrics<distype> g_current;

  DRT::ELEMENTS::SHELL::Strains strains;
  // init enhanced strain for shell
  CORE::LINALG::SerialDenseVector strain_enh(SHELL::DETAIL::num_internal_variables);

  // init EAS shape function matrix
  CORE::LINALG::SerialDenseMatrix M(SHELL::DETAIL::num_internal_variables, locking_types_.total);

  SHELL::ForEachGaussPoint<distype>(nodal_coordinates, intpoints_midsurface_,
      [&](const std::array<double, 2>& xi_gp,
          const SHELL::ShapefunctionsAndDerivatives<distype>& shape_functions,
          SHELL::BasisVectorsAndMetrics<distype>& a_current,
          SHELL::BasisVectorsAndMetrics<distype>& a_reference, double gpweight, double da, int gp)
      {
        double integration_factor = gpweight * da;

        //  make shape functions for incompatible strains
        M = SHELL::EAS::EvaluateEasShapeFunctions(
            xi_gp, locking_types_, a_reference, metrics_centroid_reference);
        SHELL::EAS::EvaluateEasStrains(strain_enh, eas_iteration_data_.alpha_, M);

        // integration loop in thickness direction, here we prescribe 2 integration points
        for (int gpt = 0; gpt < intpoints_thickness_.NumPoints(); ++gpt)
        {
          zeta = intpoints_thickness_.qxg[gpt][0] / condfac;

          SHELL::EvaluateMetrics(shape_functions, g_reference, g_current, nodal_coordinates, zeta);

          // modify the current kovariant metric tensor to neglect the quadratic terms in thickness
          // directions
          if (shell_data_.num_ans > 0)
          {
            // modify the current kovariant metric tensor due to transverse shear strain ANS
            SHELL::ModifyKovariantMetricsAns(g_reference, g_current, a_reference, a_current, zeta,
                shape_functions_ans, metrics_collocation_reference, metrics_collocation_current,
                shell_data_.num_ans);
          }
          else
          {
            SHELL::ModifyKovariantMetrics(g_reference, g_current, a_reference, a_current, zeta);
          }

          // change to current metrics due to eas
          SHELL::EAS::UpdateCurrentMetricsEAS(g_current, strain_enh, zeta);

          // evaluate Green-Lagrange strains and deformation gradient in cartesian coordinate system
          strains = SHELL::EvaluateStrains(g_reference, g_current);

          // call material for evaluation of strain energy function
          double psi = 0.0;

          solid_material.StrainEnergy(strains.gl_strain_, psi, gp, ele.Id());

          double thickness = 0.0;
          for (int i = 0; i < SHELL::DETAIL::num_node<distype>; ++i)
            thickness += thickness * shape_functions.shapefunctions_(i);

          intenergy += psi * integration_factor * 0.5 * thickness;
        }
      });

  return intenergy;
}


template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Shell7pEleCalcEas<distype>::CalculateStressesStrains(DRT::Element& ele,
    MAT::So3Material& solid_material, const ShellStressIO& stressIO, const ShellStrainIO& strainIO,
    const DRT::Discretization& discretization,
    const CORE::LINALG::SerialDenseMatrix& nodal_directors, const std::vector<int>& dof_index_array,
    Teuchos::ParameterList& params)
{
  if (discretization.Comm().MyPID() != ele.Owner()) return;

  std::vector<char>& serialized_stress_data = stressIO.mutable_data;
  std::vector<char>& serialized_strain_data = strainIO.mutable_data;
  CORE::LINALG::SerialDenseMatrix stress_data(
      intpoints_midsurface_.NumPoints(), MAT::NUM_STRESS_3D);
  CORE::LINALG::SerialDenseMatrix strain_data(
      intpoints_midsurface_.NumPoints(), MAT::NUM_STRESS_3D);

  Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
  Teuchos::RCP<const Epetra_Vector> res = discretization.GetState("residual displacement");
  std::vector<double> displacement(dof_index_array.size());
  DRT::UTILS::ExtractMyValues(*disp, displacement, dof_index_array);
  std::vector<double> residual(dof_index_array.size());
  DRT::UTILS::ExtractMyValues(*res, residual, dof_index_array);

  // init gauss point in thickness direction that will be modified via SDC
  double zeta = 0.0;

  // init scale factor for scaled director approach (SDC)
  const double condfac = shell_data_.sdc;

  // get nodal coordinates
  SHELL::NodalCoordinates<distype> nodal_coordinates = SHELL::EvaluateNodalCoordinates<distype>(
      ele.Nodes(), displacement, shell_data_.thickness, nodal_directors, condfac);

  // Enhanced Assumed Strain (EAS) Technology: declare, initialize, set up, and alpha history

  // EAS Update of alphas: the current alphas are (re-)evaluated out of DTilde and L^T of previous
  // step to avoid additional element call
  CORE::LINALG::SerialDenseMatrix delta_alpha(locking_types_.total, 1);

  if (not ele.IsParamsInterface())
  {
    // compute the EAS increment delta_alpha
    EvaluateAlphaIncrement<distype>(
        eas_iteration_data_, locking_types_.total, residual, delta_alpha);
    // update alpha += 1.0 * delta_alpha
    eas_iteration_data_.alpha_ += delta_alpha;
  }

  // Assumed Natural Strains (ANS) Technology to remedy transverse shear strain locking
  std::vector<double> shape_functions_ans(true);
  // for a_13 and a_23 each
  const int total_ansq = 2 * shell_data_.num_ans;
  std::vector<SHELL::ShapefunctionsAndDerivatives<distype>> shapefunctions_collocation(total_ansq);
  std::vector<SHELL::BasisVectorsAndMetrics<distype>> metrics_collocation_reference(total_ansq);
  std::vector<SHELL::BasisVectorsAndMetrics<distype>> metrics_collocation_current(total_ansq);

  if (shell_data_.num_ans > 0)
  {
    SHELL::SetupANS(shapefunctions_collocation, metrics_collocation_reference,
        metrics_collocation_current, nodal_coordinates, total_ansq);
  }

  // metric of element at centroid point (for EAS)
  const std::array<double, 2> centroid_point = {0.0, 0.0};
  SHELL::ShapefunctionsAndDerivatives<distype> shapefunctions_centroid =
      SHELL::EvaluateShapefunctionsAndDerivs<distype>(centroid_point);
  SHELL::BasisVectorsAndMetrics<distype> metrics_centroid_reference;
  SHELL::BasisVectorsAndMetrics<distype> metrics_centroid_current;

  SHELL::EvaluateMetrics(shapefunctions_centroid, metrics_centroid_reference,
      metrics_centroid_current, nodal_coordinates, 0.0);

  // init metric tensor and basis vectors of element mid-surface
  SHELL::BasisVectorsAndMetrics<distype> a_reference;
  SHELL::BasisVectorsAndMetrics<distype> a_current;

  // init metric tensor and basis vectors of element shell body
  SHELL::BasisVectorsAndMetrics<distype> g_reference;
  SHELL::BasisVectorsAndMetrics<distype> g_current;

  // init enhanced strains for shell
  CORE::LINALG::SerialDenseVector strain_enh(SHELL::DETAIL::num_internal_variables);
  SHELL::StressEnhanced stress_enh;

  // init EAS shape function matrix
  CORE::LINALG::SerialDenseMatrix M(SHELL::DETAIL::num_internal_variables, locking_types_.total);

  SHELL::ForEachGaussPoint<distype>(nodal_coordinates, intpoints_midsurface_,
      [&](const std::array<double, 2>& xi_gp,
          const SHELL::ShapefunctionsAndDerivatives<distype>& shape_functions,
          SHELL::BasisVectorsAndMetrics<distype>& a_current,
          SHELL::BasisVectorsAndMetrics<distype>& a_reference, double gpweight, double da, int gp)
      {
        //  evaluate shape functions for incompatible strains
        M = SHELL::EAS::EvaluateEasShapeFunctions(
            xi_gp, locking_types_, a_reference, metrics_centroid_reference);

        SHELL::EAS::EvaluateEasStrains(strain_enh, eas_iteration_data_.alpha_, M);

        // integration loop in thickness direction, here we prescribe 2 integration points to avoid
        // nonlinear poisson stiffening
        for (int gpt = 0; gpt < intpoints_thickness_.NumPoints(); ++gpt)
        {
          zeta = intpoints_thickness_.qxg[gpt][0] / condfac;
          SHELL::EvaluateMetrics(shape_functions, g_reference, g_current, nodal_coordinates, zeta);

          // modify the current kovariant metric tensor to neglect the quadratic terms in thickness
          // directions
          if (shell_data_.num_ans > 0)
          {
            // modify the current kovariant metric tensor due to transverse shear strain ANS
            SHELL::ModifyKovariantMetricsAns(g_reference, g_current, a_reference, a_current, zeta,
                shape_functions_ans, metrics_collocation_reference, metrics_collocation_current,
                shell_data_.num_ans);
          }
          else
          {
            SHELL::ModifyKovariantMetrics(g_reference, g_current, a_reference, a_current, zeta);
          }

          // change to current metrics due to eas
          SHELL::EAS::UpdateCurrentMetricsEAS(g_current, strain_enh, zeta);

          // evaluate Green-Lagrange strains and deformationgradient in cartesian coordinate system
          auto strains = SHELL::EvaluateStrains(g_reference, g_current);

          // update the deformation gradient
          CORE::LINALG::Matrix<SHELL::DETAIL::num_dim, SHELL::DETAIL::num_dim> defgrd_enh(false);
          SHELL::CalcConsistentDefgrd<SHELL::DETAIL::num_dim>(
              strains.defgrd_, strains.gl_strain_, defgrd_enh);
          strains.defgrd_ = defgrd_enh;

          // evaluate stress in local cartesian system
          auto stress = SHELL::EvaluateMaterialStressCartesianSystem<SHELL::DETAIL::num_dim>(
              solid_material, strains, params, gp, ele.Id());
          SHELL::AssembleStrainTypeToMatrixRow<distype>(
              strains, strainIO.type, strain_data, gp, 0.5);
          SHELL::AssembleStressTypeToMatrixRow<distype>(
              strains, stress, stressIO.type, stress_data, gp, 0.5);
        }
      });
  SHELL::Serialize(stress_data, serialized_stress_data);
  SHELL::Serialize(strain_data, serialized_strain_data);
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Shell7pEleCalcEas<distype>::EvaluateNonlinearForceStiffnessMass(
    DRT::Element& ele, MAT::So3Material& solid_material, const DRT::Discretization& discretization,
    const CORE::LINALG::SerialDenseMatrix& nodal_directors, const std::vector<int>& dof_index_array,
    Teuchos::ParameterList& params, CORE::LINALG::SerialDenseVector* force_vector,
    CORE::LINALG::SerialDenseMatrix* stiffness_matrix, CORE::LINALG::SerialDenseMatrix* mass_matrix)
{
  Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
  Teuchos::RCP<const Epetra_Vector> res = discretization.GetState("residual displacement");
  std::vector<double> displacement(dof_index_array.size());
  DRT::UTILS::ExtractMyValues(*disp, displacement, dof_index_array);
  std::vector<double> residual(dof_index_array.size());
  DRT::UTILS::ExtractMyValues(*res, residual, dof_index_array);

  // init gauss point in thickness direction that will be modified via SDC
  double zeta = 0.0;

  // init scale factor for scaled director approach (SDC)
  const double condfac = shell_data_.sdc;

  // get nodal coordinates
  SHELL::NodalCoordinates<distype> nodal_coordinates = SHELL::EvaluateNodalCoordinates<distype>(
      ele.Nodes(), displacement, shell_data_.thickness, nodal_directors, condfac);

  // Enhanced Assumed Strain (EAS) Technology: declare, initialize, set up, and alpha history

  // EAS Update of alphas: the current alphas are (re-)evaluated out of DTilde and L^T of previous
  // step to avoid additional element call
  CORE::LINALG::SerialDenseMatrix delta_alpha(locking_types_.total, 1);

  if (not ele.IsParamsInterface())
  {
    // compute the EAS increment delta_alpha
    EvaluateAlphaIncrement<distype>(
        eas_iteration_data_, locking_types_.total, residual, delta_alpha);
    // update alpha += 1.0 * delta_alpha
    eas_iteration_data_.alpha_ += delta_alpha;
  }

  // clear EAS data for integration
  eas_iteration_data_.RTilde_.shape(locking_types_.total, 1);
  eas_iteration_data_.invDTilde_.shape(locking_types_.total, locking_types_.total);
  eas_iteration_data_.transL_.shape(locking_types_.total, SHELL::DETAIL::numdofperelement<distype>);

  // Assumed Natural Strains (ANS) Technology to remedy transverse shear strain locking
  std::vector<double> shape_functions_ans(true);
  // for a_13 and a_23 each
  const int total_ansq = 2 * shell_data_.num_ans;
  std::vector<SHELL::ShapefunctionsAndDerivatives<distype>> shapefunctions_collocation(total_ansq);
  std::vector<SHELL::BasisVectorsAndMetrics<distype>> metrics_collocation_reference(total_ansq);
  std::vector<SHELL::BasisVectorsAndMetrics<distype>> metrics_collocation_current(total_ansq);

  if (shell_data_.num_ans > 0)
  {
    SHELL::SetupANS(shapefunctions_collocation, metrics_collocation_reference,
        metrics_collocation_current, nodal_coordinates, total_ansq);
  }

  // metric of element at centroid point (for EAS)
  const std::array<double, 2> centroid_point = {0.0, 0.0};
  SHELL::ShapefunctionsAndDerivatives<distype> shapefunctions_centroid =
      SHELL::EvaluateShapefunctionsAndDerivs<distype>(centroid_point);
  SHELL::BasisVectorsAndMetrics<distype> metrics_centroid_reference;
  SHELL::BasisVectorsAndMetrics<distype> metrics_centroid_current;

  SHELL::EvaluateMetrics(shapefunctions_centroid, metrics_centroid_reference,
      metrics_centroid_current, nodal_coordinates, 0.0);

  // init metric tensor and basis vectors of element mid-surface
  SHELL::BasisVectorsAndMetrics<distype> a_reference;
  SHELL::BasisVectorsAndMetrics<distype> a_current;

  // init metric tensor and basis vectors of element shell body
  SHELL::BasisVectorsAndMetrics<distype> g_reference;
  SHELL::BasisVectorsAndMetrics<distype> g_current;

  // init enhanced strain for shell
  constexpr auto num_internal_variables = SHELL::DETAIL::num_internal_variables;
  CORE::LINALG::SerialDenseVector strain_enh(num_internal_variables);
  SHELL::StressEnhanced stress_enh;

  // init EAS shape function matrix
  CORE::LINALG::SerialDenseMatrix M(num_internal_variables, locking_types_.total);

  SHELL::ForEachGaussPoint<distype>(nodal_coordinates, intpoints_midsurface_,
      [&](const std::array<double, 2>& xi_gp,
          const SHELL::ShapefunctionsAndDerivatives<distype>& shape_functions,
          SHELL::BasisVectorsAndMetrics<distype>& a_current,
          SHELL::BasisVectorsAndMetrics<distype>& a_reference, double gpweight, double da, int gp)
      {
        double integration_factor = gpweight * da;

        // update current thickness at gauss point
        cur_thickness_[gp] = SHELL::UpdateGaussPointThickness<distype>(
            nodal_coordinates.a3_curr_, shape_functions.shapefunctions_);

        // reset mid-surface material tensor and stress resultants to zero
        stress_enh.dmat_.shape(num_internal_variables, num_internal_variables);
        stress_enh.stress_.size(num_internal_variables);

        // init mass matrix variables
        SHELL::MassMatrixVariables mass_matrix_variables;

        //  evaluate shape functions for incompatible strains
        M = SHELL::EAS::EvaluateEasShapeFunctions(
            xi_gp, locking_types_, a_reference, metrics_centroid_reference);
        SHELL::EAS::EvaluateEasStrains(strain_enh, eas_iteration_data_.alpha_, M);

        // calculate B-operator for compatible strains (displacement)
        CORE::LINALG::SerialDenseMatrix Bop = SHELL::CalcBOperator<distype>(
            a_current.kovariant_, a_current.partial_derivative_, shape_functions);

        // modifications due to ANS with B-bar method (Hughes (1980))
        if (shell_data_.num_ans > 0)
        {
          shape_functions_ans = SHELL::GetShapefunctionsForAns<distype>(xi_gp);
          SHELL::ModifyBOperatorAns(Bop, shape_functions_ans, shapefunctions_collocation,
              metrics_collocation_current, shell_data_.num_ans);
        }

        // integration loop in thickness direction, here we prescribe 2 integration points to avoid
        // nonlinear poisson stiffening
        for (int gpt = 0; gpt < intpoints_thickness_.NumPoints(); ++gpt)
        {
          zeta = intpoints_thickness_.qxg[gpt][0] / condfac;
          double factor = intpoints_thickness_.qwgt[gpt];

          SHELL::EvaluateMetrics(shape_functions, g_reference, g_current, nodal_coordinates, zeta);

          // evaluate metric tensor at gp in shell body
          if (shell_data_.num_ans > 0)
          {
            // modify the current kovariant metric tensor due to transverse shear strain ANS
            SHELL::ModifyKovariantMetricsAns(g_reference, g_current, a_reference, a_current, zeta,
                shape_functions_ans, metrics_collocation_reference, metrics_collocation_current,
                shell_data_.num_ans);
          }
          else
          {
            SHELL::ModifyKovariantMetrics(g_reference, g_current, a_reference, a_current, zeta);
          }

          // calc shell shifter and put it in the integration factor
          factor *= (1.0 / condfac) * (g_reference.detJ_ / da);

          // change to current metrics due to EAS
          SHELL::EAS::UpdateCurrentMetricsEAS(g_current, strain_enh, zeta);

          // evaluate Green-Lagrange strains and deformation gradient in cartesian coordinate system
          auto strains = SHELL::EvaluateStrains(g_reference, g_current);

          // update the deformation gradient (if needed?)
          if (solid_material.NeedsDefgrd())
          {
            CORE::LINALG::Matrix<SHELL::DETAIL::num_dim, SHELL::DETAIL::num_dim> defgrd_enh(false);
            SHELL::CalcConsistentDefgrd<SHELL::DETAIL::num_dim>(
                strains.defgrd_, strains.gl_strain_, defgrd_enh);
            strains.defgrd_ = defgrd_enh;
          }

          auto stress = SHELL::EvaluateMaterialStressCartesianSystem<SHELL::DETAIL::num_dim>(
              solid_material, strains, params, gp, ele.Id());
          SHELL::MapMaterialStressToCurvilinearSystem(stress, g_reference);
          SHELL::ThicknessIntegration<distype>(stress_enh, stress, factor, zeta);
          // thickness integration of mass matrix variables
          if (mass_matrix != nullptr)
          {
            double tmp_integration_factor = intpoints_thickness_.qwgt[gpt] * g_reference.detJ_;
            mass_matrix_variables.factor_v_ += tmp_integration_factor;
            mass_matrix_variables.factor_w_ += tmp_integration_factor *
                                               intpoints_thickness_.qxg[gpt][0] *
                                               intpoints_thickness_.qxg[gpt][0];
            mass_matrix_variables.factor_vw_ +=
                tmp_integration_factor * intpoints_thickness_.qxg[gpt][0];
          }
        }

        // integration of EAS matrices
        IntegrateEAS<distype>(
            stress_enh, M, Bop, eas_iteration_data_, integration_factor, locking_types_.total);

        // add stiffness matrix
        if (stiffness_matrix != nullptr)
        {
          // elastic stiffness matrix Ke
          SHELL::AddElasticStiffnessMatrix<distype>(
              Bop, stress_enh.dmat_, integration_factor, *stiffness_matrix);
          // geometric stiffness matrix Kg
          SHELL::AddGeometricStiffnessMatrix(shapefunctions_collocation, shape_functions_ans,
              shape_functions, stress_enh.stress_, shell_data_.num_ans, integration_factor,
              *stiffness_matrix);
        }
        // add internal force vector
        if (force_vector != nullptr)
        {
          SHELL::AddInternalForceVector<distype>(
              Bop, stress_enh.stress_, integration_factor, *force_vector);
        }
        // add internal mass_matrix
        if (mass_matrix != nullptr)
        {
          double density = solid_material.Density(gp);
          mass_matrix_variables.factor_v_ *= gpweight * density;
          mass_matrix_variables.factor_w_ *= gpweight * density;
          mass_matrix_variables.factor_vw_ *= gpweight * density;
          SHELL::AddMassMatrix(
              shape_functions, mass_matrix_variables, shell_data_.thickness, *mass_matrix);
        }
      });

  // compute inverse of DTilde = invDTilde
  CORE::LINALG::SymmetricInverse(eas_iteration_data_.invDTilde_, locking_types_.total);

  // compute  L * DTilde^-1  which is later needed for force and stiffness update
  CORE::LINALG::SerialDenseMatrix LinvDTilde(
      SHELL::DETAIL::numdofperelement<distype>, locking_types_.total);
  LinvDTilde.multiply(Teuchos::TRANS, Teuchos::NO_TRANS, 1., eas_iteration_data_.transL_,
      eas_iteration_data_.invDTilde_, 0.);
  if (stiffness_matrix != nullptr)
  {
    SHELL::EAS::AddEASStiffnessMatrix(LinvDTilde, eas_iteration_data_.transL_, *stiffness_matrix);
  }

  if (force_vector != nullptr)
  {
    SHELL::EAS::AddEASInternalForce(LinvDTilde, eas_iteration_data_.RTilde_, *force_vector);
  }

  if (stiffness_matrix != nullptr)
  {
    // make stiffness matrix absolute symmetric
    for (int i = 0; i < SHELL::DETAIL::numdofperelement<distype>; ++i)
    {
      for (int j = i + 1; j < SHELL::DETAIL::numdofperelement<distype>; ++j)
      {
        const double average = 0.5 * ((*stiffness_matrix)(i, j) + (*stiffness_matrix)(j, i));
        (*stiffness_matrix)(i, j) = average;
        (*stiffness_matrix)(j, i) = average;
      }
    }
  }
}



template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Shell7pEleCalcEas<distype>::Recover(DRT::Element& ele,
    const DRT::Discretization& discretization, const std::vector<int>& dof_index_array,
    Teuchos::ParameterList& params, STR::ELEMENTS::ParamsInterface& interface_ptr)
{
  Teuchos::RCP<const Epetra_Vector> res = discretization.GetState("residual displacement");
  if (res == Teuchos::null) dserror("Cannot get residual displacement state vector");
  std::vector<double> residual(dof_index_array.size());
  DRT::UTILS::ExtractMyValues(*res, residual, dof_index_array);

  // get access to the interface parameters
  double step_length = interface_ptr.GetStepLength();

  // access general EAS history stuff stored in element
  CORE::LINALG::SerialDenseMatrix delta_alpha(locking_types_.total, 1);

  // if it is a default step, we have to recover the condensed solution vectors
  if (interface_ptr.IsDefaultStep())
  {
    // first, store the eas state of the previous accepted Newton step
    interface_ptr.SumIntoMyPreviousSolNorm(NOX::NLN::StatusTest::quantity_eas, locking_types_.total,
        eas_iteration_data_.alpha_[0], ele.Owner());

    // compute the EAS increment delta_alpha
    EvaluateAlphaIncrement<distype>(
        eas_iteration_data_, locking_types_.total, residual, delta_alpha);
    // update alpha += step_length * delta_alpha
    delta_alpha.scale(step_length);
    eas_iteration_data_.alpha_ += delta_alpha;
  }
  // if it is no default step, we can correct the update and the current eas state without the
  // need for any matrix-vector products.
  else
  {
    // The first step has to be a default step!
    if (old_step_length_ < 0.0) dserror("The old step length was not defined!");
    // if this is no full step, we have to adjust the length of the enhanced assumed strain
    // incremental step.
    // undo the previous step:
    //            alpha_new = alpha_old - old_step * delta_alpha
    // and update the solution variable with the new step length:
    //           alpha_new = alpha_new + new_step * delta_alpha
    delta_alpha.scale(step_length - old_step_length_);
    eas_iteration_data_.alpha_ += delta_alpha;
  }

  // Check if delta alpha is tested and if yes, calculate the element
  // contribution to the norm
  interface_ptr.SumIntoMyUpdateNorm(NOX::NLN::StatusTest::quantity_eas, locking_types_.total,
      delta_alpha[0], eas_iteration_data_.alpha_[0], step_length, ele.Owner());

  // save the old step length
  old_step_length_ = step_length;
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Shell7pEleCalcEas<distype>::Update(DRT::Element& ele,
    MAT::So3Material& solid_material, const DRT::Discretization& discretization,
    const CORE::LINALG::SerialDenseMatrix& nodal_directors, const std::vector<int>& dof_index_array,
    Teuchos::ParameterList& params)
{
  Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
  if (disp == Teuchos::null) dserror("Cannot get state vectors 'displacement' ");
  std::vector<double> displacement(dof_index_array.size());
  DRT::UTILS::ExtractMyValues(*disp, displacement, dof_index_array);

  // No need to update alpha here. Update is called to copy states from t_{n+1} to
  // t_{n} after the time step and output Hence, there are no more Newton iterations that would
  // require an update of alpha

  // calculate and update inelastic deformation gradient if needed
  if (solid_material.UsesExtendedUpdate())
  {
    // init scale factor for scaled director approach (SDC)
    const double condfac = shell_data_.sdc;

    // init gauss point in thickness direction that will be modified via SDC
    double zeta = 0.0;

    // get nodal coordinates
    SHELL::NodalCoordinates<distype> nodal_coordinates = SHELL::EvaluateNodalCoordinates<distype>(
        ele.Nodes(), displacement, shell_data_.thickness, nodal_directors, condfac);

    // metric of element centroid point (for EAS)
    const std::array<double, 2> centroid_point = {0.0, 0.0};
    SHELL::ShapefunctionsAndDerivatives<distype> shapefunctions_centroid =
        SHELL::EvaluateShapefunctionsAndDerivs<distype>(centroid_point);
    SHELL::BasisVectorsAndMetrics<distype> metrics_centroid_reference;
    SHELL::BasisVectorsAndMetrics<distype> metrics_centroid_current;

    SHELL::EvaluateMetrics(shapefunctions_centroid, metrics_centroid_reference,
        metrics_centroid_current, nodal_coordinates, 0.0);

    // Assumed Natural Strains (ANS) Technology to remedy transverse shear strain locking
    std::vector<double> shape_functions_ans(true);
    // for a_13 and a_23 each
    const int total_ansq = 2 * shell_data_.num_ans;
    std::vector<SHELL::ShapefunctionsAndDerivatives<distype>> shapefunctions_collocation(
        total_ansq);
    std::vector<SHELL::BasisVectorsAndMetrics<distype>> metrics_collocation_reference(total_ansq);
    std::vector<SHELL::BasisVectorsAndMetrics<distype>> metrics_collocation_current(total_ansq);

    if (shell_data_.num_ans > 0)
    {
      SHELL::SetupANS(shapefunctions_collocation, metrics_collocation_reference,
          metrics_collocation_current, nodal_coordinates, total_ansq);
    }

    // init metric tensor and basis vectors of mid-surface
    SHELL::BasisVectorsAndMetrics<distype> a_reference;
    SHELL::BasisVectorsAndMetrics<distype> a_current;

    // init metric tensor and basis vectors of shell body
    SHELL::BasisVectorsAndMetrics<distype> g_reference;
    SHELL::BasisVectorsAndMetrics<distype> g_current;

    // enhanced strain for shell
    constexpr auto num_internal_variables = SHELL::DETAIL::num_internal_variables;
    CORE::LINALG::SerialDenseVector strain_enh(num_internal_variables);

    // init EAS shape function matrix
    CORE::LINALG::SerialDenseMatrix M(num_internal_variables, num_internal_variables);

    SHELL::ForEachGaussPoint<distype>(nodal_coordinates, intpoints_midsurface_,
        [&](const std::array<double, 2>& xi_gp,
            const SHELL::ShapefunctionsAndDerivatives<distype>& shape_functions,
            SHELL::BasisVectorsAndMetrics<distype>& a_current,
            SHELL::BasisVectorsAndMetrics<distype>& a_reference, double gpweight, double da, int gp)
        {
          //  make shape functions for incompatible strains  M
          M = SHELL::EAS::EvaluateEasShapeFunctions(
              xi_gp, locking_types_, a_reference, metrics_centroid_reference);
          SHELL::EAS::EvaluateEasStrains(strain_enh, eas_iteration_data_.alpha_, M);

          // integration loop in thickness direction, here we prescribe 2 integration points
          for (int gpt = 0; gpt < intpoints_thickness_.NumPoints(); ++gpt)
          {
            zeta = intpoints_thickness_.qxg[gpt][0] / condfac;

            SHELL::EvaluateMetrics(
                shape_functions, g_reference, g_current, nodal_coordinates, zeta);

            // modify the current kovariant metric tensor to neglect the quadratic terms in
            // thickness directions
            if (shell_data_.num_ans > 0)
            {
              // modify the current kovariant metric tensor due to transverse shear strain ANS
              SHELL::ModifyKovariantMetricsAns(g_reference, g_current, a_reference, a_current, zeta,
                  shape_functions_ans, metrics_collocation_reference, metrics_collocation_current,
                  shell_data_.num_ans);
            }
            else
            {
              SHELL::ModifyKovariantMetrics(g_reference, g_current, a_reference, a_current, zeta);
            }

            SHELL::EAS::UpdateCurrentMetricsEAS(g_current, strain_enh, zeta);

            auto strains = EvaluateStrains(g_reference, g_current);

            // calculate deformation gradient consistent with modified GL strain tensor
            if (solid_material.NeedsDefgrd())
            {
              // update the deformation gradient (if needed)
              CORE::LINALG::Matrix<DRT::ELEMENTS::SHELL::DETAIL::num_dim,
                  DRT::ELEMENTS::SHELL::DETAIL::num_dim>
                  defgrd_enh(false);
              SHELL::CalcConsistentDefgrd<SHELL::DETAIL::num_dim>(
                  strains.defgrd_, strains.gl_strain_, defgrd_enh);
              strains.defgrd_ = defgrd_enh;
            }
            solid_material.Update(strains.defgrd_, gp, params, ele.Id());
          }
        });
  }
  solid_material.Update();
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Shell7pEleCalcEas<distype>::VisData(
    const std::string& name, std::vector<double>& data)
{
  if (name == "thickness")
  {
    if (data.size() != 1) dserror("size mismatch");
    for (auto& thickness_data : cur_thickness_)
    {
      data[0] += thickness_data;
    }
    data[0] = data[0] / intpoints_midsurface_.NumPoints();
  }

}  // VisData()

// template classes
template class DRT::ELEMENTS::Shell7pEleCalcEas<DRT::Element::quad4>;
template class DRT::ELEMENTS::Shell7pEleCalcEas<DRT::Element::quad8>;
template class DRT::ELEMENTS::Shell7pEleCalcEas<DRT::Element::quad9>;
template class DRT::ELEMENTS::Shell7pEleCalcEas<DRT::Element::tri3>;
template class DRT::ELEMENTS::Shell7pEleCalcEas<DRT::Element::tri6>;
