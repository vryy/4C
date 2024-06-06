/*----------------------------------------------------------------------*/
/*! \file

\brief unit testing functionality for the mirco contact constitutivelaw with pressure-based Green
function

\level 2

*/
/*----------------------------------------------------------------------*/

#include <gtest/gtest.h>

#include "4C_contact_constitutivelaw_mirco_contactconstitutivelaw.hpp"
#include "4C_contact_rough_node.hpp"
#include "4C_global_data.hpp"
#include "4C_mat_material_factory.hpp"
#include "4C_mat_par_bundle.hpp"
#include "4C_utils_function.hpp"

#ifdef FOUR_C_WITH_MIRCO

namespace
{
  using namespace FourC;

  class MircoConstitutiveLawPressureTest : public ::testing::Test
  {
   public:
    MircoConstitutiveLawPressureTest()
    {
      const int problemid(0);
      Global::Problem& problem = (*Global::Problem::Instance());
      problem.Materials()->SetReadFromProblem(problemid);

      Teuchos::RCP<Core::UTILS::SymbolicFunctionOfSpaceTime<1>> FFUNCT1 =
          Teuchos::rcp(new Core::UTILS::SymbolicFunctionOfSpaceTime<1>({"0.7"}, {}));
      Teuchos::RCP<Core::UTILS::SymbolicFunctionOfSpaceTime<1>> FFUNCT2 =
          Teuchos::rcp(new Core::UTILS::SymbolicFunctionOfSpaceTime<1>({"20.0"}, {}));

      Teuchos::RCP<Core::UTILS::FunctionOfSpaceTime> FUNCT1 = FFUNCT1;
      Teuchos::RCP<Core::UTILS::FunctionOfSpaceTime> FUNCT2 = FFUNCT2;

      Core::UTILS::FunctionManager functionmanager_;
      functionmanager_.SetFunctions<std::any>({FUNCT1, FUNCT2});
      problem.SetFunctionManager(std::move(functionmanager_));

      // set up material to be added to problem instance
      Core::IO::InputParameterContainer mat_stvenant;
      mat_stvenant.Add("YOUNG", 1.0);
      mat_stvenant.Add("NUE", 0.3);
      mat_stvenant.Add("DENS", 1.0);

      problem.Materials()->insert(
          1, Mat::make_parameter(1, Core::Materials::MaterialType::m_stvenant, mat_stvenant));

      // initialize container for material parameters
      const Teuchos::RCP<CONTACT::CONSTITUTIVELAW::Container> container =
          Teuchos::rcp(new CONTACT::CONSTITUTIVELAW::Container(
              1, Inpar::CONTACT::ConstitutiveLawType::colaw_mirco, "Mirco Constitutivelaw"));

      // add parameters to container
      container->Add("FirstMatID", 1);
      container->Add("SecondMatID", 1);
      container->Add("LateralLength", 1000.0);
      container->Add("Resolution", 6);
      container->Add("PressureGreenFunFlag", true);
      container->Add("InitialTopologyStdDeviationFunct", 2);
      container->Add("HurstExponentFunct", 1);
      container->Add("RandomTopologyFlag", true);
      container->Add("RandomSeedFlag", false);
      container->Add("RandomGeneratorSeed", 95);
      container->Add("Tolerance", 0.01);
      container->Add("MaxIteration", 100);
      container->Add("WarmStartingFlag", true);
      container->Add("Offset", 2.0);
      container->Add("FiniteDifferenceFraction", 0.001);
      container->Add("ActiveGapTolerance", 1e-6);
      container->Add("TopologyFilePath", std::string("sup6.dat"));

      const Teuchos::RCP<CONTACT::CONSTITUTIVELAW::ConstitutiveLaw> mircococonstlaw =
          CONTACT::CONSTITUTIVELAW::ConstitutiveLaw::Factory(container);
      coconstlaw_ = mircococonstlaw;

      std::vector<double> x(3, 0.0);
      std::vector<int> dofs(3);

      double hurstexponentfunct = container->Get<int>("HurstExponentFunct");
      double initialtopologystddevfunct = container->Get<int>("InitialTopologyStdDeviationFunct");
      int resolution = container->Get<int>("Resolution");
      bool randomtopologyflag = container->Get<bool>("RandomTopologyFlag");
      bool randomseedflag = container->Get<bool>("RandomSeedFlag");
      int randomgeneratorseed = container->Get<int>("RandomGeneratorSeed");

      cnode = Teuchos::rcp(new CONTACT::RoughNode(1, x, 1, dofs, true, true, hurstexponentfunct,
          initialtopologystddevfunct, resolution, randomtopologyflag, randomseedflag,
          randomgeneratorseed));
    }
    Teuchos::RCP<CONTACT::CONSTITUTIVELAW::ConstitutiveLaw> coconstlaw_;

    Teuchos::RCP<CONTACT::Node> cnode;
  };

  //! test member function Evaluate
  TEST_F(MircoConstitutiveLawPressureTest, TestEvaluate)
  {
    // gap < 0
    EXPECT_ANY_THROW(coconstlaw_->Evaluate(1.0, cnode.get()));
    // 0< gap < offset
    EXPECT_ANY_THROW(coconstlaw_->Evaluate(-0.25, cnode.get()));
    // offset < gap
    EXPECT_NEAR(coconstlaw_->Evaluate(-12.0, cnode.get()), -0.0005861475487657709, 1.e-10);
  }

  //! test member function EvaluateDeriv
  TEST_F(MircoConstitutiveLawPressureTest, TestEvaluateDeriv)
  {
    EXPECT_NEAR(coconstlaw_->EvaluateDeriv(-12.0, cnode.get()), 1.56329102801896e-04, 1.e-10);
    EXPECT_ANY_THROW(coconstlaw_->EvaluateDeriv(-0.25, cnode.get()));
  }
}  // namespace

#endif
