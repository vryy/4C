/*----------------------------------------------------------------------*/
/*! \file
\brief unit testing functionality for stvenantkirchhoff material
\level 2

*/
/*----------------------------------------------------------------------*/

#include "gtest/gtest.h"

#include "matpar_material.H"
#include "stvenantkirchhoff.H"

#include "unittests_assertions.h"

namespace
{
  class StVenantKirchhoffTest : public ::testing::Test
  {
   protected:
    void SetUp() override
    {
      // initialize container for material parameters
      const Teuchos::RCP<MAT::PAR::Material> container = Teuchos::rcp(new MAT::PAR::Material());

      // add material parameters to container
      container->Add("YOUNG", young_);
      container->Add("NUE", nu_);
      container->Add("DENS", rho_);
      container->Add("THEXPANS", thermal_);

      // initialize parameter class for StVenantKirchhoff material with container
      parameters_stvenantkirchhoff_ = Teuchos::rcp(new MAT::PAR::StVenantKirchhoff(container));

      // initialize stvenantkirchhoff material with parameter class
      stvenantkirchhoff_ =
          Teuchos::rcp(new MAT::StVenantKirchhoff(parameters_stvenantkirchhoff_.get()));
    }

    //! material parameters
    const double young_ = 210.;
    const double nu_ = 0.3;
    const double rho_ = 1.0;      // dummy value (needed for construction)
    const double thermal_ = 1.0;  // dummy value (needed for construction)
    Teuchos::RCP<MAT::PAR::StVenantKirchhoff> parameters_stvenantkirchhoff_;

    //! material class
    Teuchos::RCP<MAT::StVenantKirchhoff> stvenantkirchhoff_;

    //! Test Green-Lagrange Strain
    std::array<double, 6> input_glstrain_ = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0};

    //! calculate reference results for stress
    const double ref_stress_normal_ =
        (young_ / ((1.0 + nu_) * (1.0 - (2.0 * nu_)))) * ((1.0 - nu_) + nu_ + nu_);

    const double ref_stress_shear_ =
        (young_ / ((1.0 + nu_) * (1.0 - (2.0 * nu_)))) * ((1.0 - (2.0 * nu_)) / 2.0);

    const std::array<double, 6> ref_stress_ = {ref_stress_normal_, ref_stress_normal_,
        ref_stress_normal_, ref_stress_shear_, ref_stress_shear_, ref_stress_shear_};
  };

  TEST_F(StVenantKirchhoffTest, TestEvaluateEpetraSerialDenseMatrix)
  {
    // Input strain
    const Epetra_SerialDenseVector input_glstrain(Copy, input_glstrain_.data(), 6);

    // Resulting material stiffness matrix
    Teuchos::RCP<Epetra_SerialDenseMatrix> result_cmat =
        Teuchos::rcp(new Epetra_SerialDenseMatrix(6, 6));

    // Resulting stress
    Teuchos::RCP<Epetra_SerialDenseVector> result_stress =
        Teuchos::rcp(new Epetra_SerialDenseVector(6));

    // Call evaluate function with test strain
    stvenantkirchhoff_->Evaluate(&input_glstrain, result_cmat.get(), result_stress.get());

    // Test member function results using reference stress values
    BACI_EXPECT_NEAR(result_stress->Values(), ref_stress_.data(), 1.0e-4);
  }

  TEST_F(StVenantKirchhoffTest, TestEvaluateLinalgMatrix)
  {
    // Resulting stress
    LINALG::Matrix<6, 1> result_stress(true);

    // Resulting material stiffness matrix
    LINALG::Matrix<6, 6> result_cmat(true);

    // Input deformation gradient, which is not used here
    LINALG::Matrix<3, 3> defgrad(true);

    // ParameterList, also not used here
    Teuchos::ParameterList paras;

    // Input strain
    const LINALG::Matrix<6, 1> input_strain(input_glstrain_.data(), false);

    // Reference stress
    const LINALG::Matrix<6, 1> ref_stress(ref_stress_.data(), false);

    // Call evaluate function with test strain
    stvenantkirchhoff_->Evaluate(
        &defgrad, &input_strain, paras, &result_stress, &result_cmat, 0, 0);

    // Test member function results using reference stress values
    BACI_EXPECT_NEAR(result_stress, ref_stress, 1.0e-4);
  }

  TEST_F(StVenantKirchhoffTest, TestStrainEnergy)
  {
    // define reference result for strain energy
    const double ref_strain_energy = 908.6538;

    // Input strain
    const LINALG::Matrix<6, 1> test_glstrain(input_glstrain_.data(), false);

    // result strain energy
    double result_psi;

    int eleGID = 1;

    // Call evaluate function with test strain
    stvenantkirchhoff_->StrainEnergy(test_glstrain, result_psi, 0, eleGID);

    // test result with respect to reference result
    EXPECT_NEAR(result_psi, ref_strain_energy, 1.0e-4);
  }
}  // namespace