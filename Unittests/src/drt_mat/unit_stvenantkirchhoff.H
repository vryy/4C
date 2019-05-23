/*----------------------------------------------------------------------*/
/*!
\file unit_stvenantkirchhoff.H

\brief unit testing functionality for stvenantkirchhoff material

\level 2

\maintainer Jan Schnabel
            jan-eike.schnabel@hzg.de
*/
/*----------------------------------------------------------------------*/

#ifndef UNITTESTS_SRC_DRT_MAT_UNIT_STVENANTKIRCHHOFF_H_
#define UNITTESTS_SRC_DRT_MAT_UNIT_STVENANTKIRCHHOFF_H_

#include <cxxtest/TestSuite.h>

#include "src/drt_mat/stvenantkirchhoff.H"

#include <Epetra_SerialDenseVector.h>
#include <Epetra_SerialDenseMatrix.h>

#include <vector>

// declaration
namespace MAT
{
  class StVenantKirchhoff_TestSuite;
}

// class implementation
class MAT::StVenantKirchhoff_TestSuite : public CxxTest::TestSuite
{
 public:
  //! set up unit test for StVenantKirchhoff material
  void setUp()
  {
    // initialize container for material parameters
    const Teuchos::RCP<MAT::PAR::Material> container = Teuchos::rcp(new MAT::PAR::Material());

    // add material parameters to container
    container->Add("YOUNG", 210.0);
    container->Add("NUE", 0.3);
    container->Add("DENS", 8.7);
    container->Add("THEXPANS", 11.0e-6);

    // initialize parameter class for StVenantKirchhoff material with container
    parameters_stvenantkirchhoff_ = Teuchos::rcp(new MAT::PAR::StVenantKirchhoff(container));

    // initialize stvenantkirchhoff material with parameter class
    stvenantkirchhoff_ =
        Teuchos::rcp(new MAT::StVenantKirchhoff(parameters_stvenantkirchhoff_.get()));

    return;
  }

  //! tear down unit testing for StVenantKirchhoff material
  void tearDown() { return; }

  //! evaluate function unit testing for StVenantKirchhoff material
  void test_Evaluate()
  {
    // calculate reference results for stress
    const double Ref_Stress_calc1 =
        (210.0 / ((1.0 + 0.3) * (1.0 - (2.0 * 0.3)))) * ((1.0 - 0.3) + 0.3 + 0.3);
    const double Ref_Stress_calc2 =
        (210.0 / ((1.0 + 0.3) * (1.0 - (2.0 * 0.3)))) * ((1.0 - (2.0 * 0.3)) / 2.0);

    const double Ref_Stress[6] = {Ref_Stress_calc1, Ref_Stress_calc1, Ref_Stress_calc1,
        Ref_Stress_calc2, Ref_Stress_calc2, Ref_Stress_calc2};

    // define tolerances for test
    const double tolerances_Evaluate_Test[] = {1.e-4, 1.e-4, 1.e-4, 1.e-4, 1.e-4, 1.e-4};

    // Input strain
    Teuchos::RCP<const Epetra_SerialDenseVector> Test_glstrain =
        Teuchos::rcp(new const Epetra_SerialDenseVector(Copy, Input_glstrain, 6));

    // Resulting material stiffness matrix
    Teuchos::RCP<Epetra_SerialDenseMatrix> Result_cmat =
        Teuchos::rcp(new Epetra_SerialDenseMatrix(6, 6));

    // Resulting stress
    Teuchos::RCP<Epetra_SerialDenseVector> Result_stress =
        Teuchos::rcp(new Epetra_SerialDenseVector(6));

    // Call evaluate function with test strain
    stvenantkirchhoff_->Evaluate(Test_glstrain.get(), Result_cmat.get(), Result_stress.get());

    // test member function results using reference stress values
    for (unsigned i = 0; i < 6; ++i)
      TS_ASSERT_DELTA((*Result_stress)[i], Ref_Stress[i], tolerances_Evaluate_Test[i]);

    return;
  }
  //! strain energy unit testing for StVenantKirchhoff material
  void test_StrainEnergy()
  {
    // define reference result for strain energy
    const double Ref_StrainEnergy = 908.6538;

    // define tolerance for result check
    const double tolerance_StrainEnergy_Test = 1.e-4;

    // Input strain
    const LINALG::Matrix<6, 1> Test_glstrain2(Input_glstrain, false);

    // result strain energy
    double Result_psi;

    int eleGID = 1;

    // Call evaluate function with test strain
    stvenantkirchhoff_->StrainEnergy(Test_glstrain2, Result_psi, eleGID);

    // test result with respect to reference result
    TS_ASSERT_DELTA(Result_psi, Ref_StrainEnergy, tolerance_StrainEnergy_Test);

    return;
  }

 private:
  //! material parameters
  Teuchos::RCP<MAT::PAR::StVenantKirchhoff> parameters_stvenantkirchhoff_;

  //! material class
  Teuchos::RCP<MAT::StVenantKirchhoff> stvenantkirchhoff_;

  //! Test Green-Lagrange Strain
  double Input_glstrain[6] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
};

#endif
