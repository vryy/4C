/*----------------------------------------------------------------------*/
/*! \file

\brief unit testing functionality for the mirco contact constitutivelaw with point-force-based Green
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

  class MircoConstitutiveLawForceTest : public ::testing::Test
  {
   public:
    MircoConstitutiveLawForceTest()
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
          1, Core::UTILS::LazyPtr<Core::Mat::PAR::Parameter>(
                 Mat::make_parameter(1, Core::Materials::MaterialType::m_stvenant, mat_stvenant)));

      // initialize container for material parameters
      const Teuchos::RCP<CONTACT::CONSTITUTIVELAW::Container> container =
          Teuchos::rcp(new CONTACT::CONSTITUTIVELAW::Container(
              1, Inpar::CONTACT::ConstitutiveLawType::colaw_mirco, "Mirco Constitutivelaw"));

      // add parameters to container
      container->Add("FirstMatID", 1);
      container->Add("SecondMatID", 1);
      container->Add("LateralLength", 1000.0);
      container->Add("Resolution", 6);
      container->Add("PressureGreenFunFlag", false);
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

      double hurstexponentfunct = container->get<int>("HurstExponentFunct");
      double initialtopologystddevfunct = container->get<int>("InitialTopologyStdDeviationFunct");
      int resolution = container->get<int>("Resolution");
      bool randomtopologyflag = container->get<bool>("RandomTopologyFlag");
      bool randomseedflag = container->get<bool>("RandomSeedFlag");
      int randomgeneratorseed = container->get<int>("RandomGeneratorSeed");

      cnode = Teuchos::rcp(new CONTACT::RoughNode(1, x, 1, dofs, true, true, hurstexponentfunct,
          initialtopologystddevfunct, resolution, randomtopologyflag, randomseedflag,
          randomgeneratorseed));
    }
    Teuchos::RCP<CONTACT::CONSTITUTIVELAW::ConstitutiveLaw> coconstlaw_;

    Teuchos::RCP<CONTACT::Node> cnode;
  };  // namespace


  //! test member function Evaluate
  TEST_F(MircoConstitutiveLawForceTest, TestEvaluate)
  {
    // gap < 0
    EXPECT_ANY_THROW(coconstlaw_->evaluate(1.0, cnode.get()));
    // 0< gap < offset
    EXPECT_ANY_THROW(coconstlaw_->evaluate(-0.25, cnode.get()));
    // offset < gap
    EXPECT_NEAR(coconstlaw_->evaluate(-12.0, cnode.get()), -0.0005942101230076477, 1.e-10);
  }

  //! test member function EvaluateDeriv
  TEST_F(MircoConstitutiveLawForceTest, TestEvaluateDeriv)
  {
    EXPECT_NEAR(coconstlaw_->EvaluateDeriv(-12.0, cnode.get()), 1.34284789678326e-04, 1.e-10);
    EXPECT_ANY_THROW(coconstlaw_->EvaluateDeriv(-0.25, cnode.get()));
  }
}  // namespace

#endif
