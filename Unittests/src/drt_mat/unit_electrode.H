/*----------------------------------------------------------------------*/
/*!
\file unit_electrode.H

\brief unit testing functionality for electrode material

\level 2

<pre>
\maintainer Rui Fang
            fang@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/
            089-289-15251
</pre>
*/
/*----------------------------------------------------------------------*/
#ifndef UNIT_ELECTRODE_H
#define UNIT_ELECTRODE_H

#include <cxxtest/TestSuite.h>

#include "src/drt_inpar/inpar_elch.H"

#include "src/drt_mat/electrode.H"
#include "src/drt_mat/matpar_material.H"

// forward declaration
namespace MAT
{
  class Electrode_TestSuite;
}

// class implementation
class MAT::Electrode_TestSuite : public CxxTest::TestSuite
{
  public:

    //! set up unit testing for electrode material
    void setUp()
    {
      // initialize container for material parameters
      const Teuchos::RCP<MAT::PAR::Material> container = Teuchos::rcp(new MAT::PAR::Material());

      // add dummy parameters to container
      container->Add("DIFFCOEF",0);
      container->Add("DIFF_PARA_NUM",0);
      container->Add("DIFF_PARA",std::vector<double>(0,0.));
      container->Add("COND",0);
      container->Add("COND_PARA_NUM",0);
      container->Add("COND_PARA",std::vector<double>(0,0.));

      // choose semi-empirical Redlich-Kister expansion to model half-cell open-circuit equilibrium potential
      container->Add("OCP_MODEL","Redlich-Kister");

      // add anode parameters to container according to Goldin et al., Electrochimica Acta 64 (2012) 118-129
      container->Add("C_MAX",16.1);
      container->Add("OCP_PARA_NUM",16);
      std::vector<double> ocp_para(16,0.);
      ocp_para[0] = 1.1652e4;
      ocp_para[1] = -3.268e3;
      ocp_para[2] = 3.955e3;
      ocp_para[3] = -4.573e3;
      ocp_para[4] = 6.147e3;
      ocp_para[5] = -3.339e3;
      ocp_para[6] = 1.117e4;
      ocp_para[7] = 2.997e2;
      ocp_para[8] = -4.866e4;
      ocp_para[9] = 1.362e2;
      ocp_para[10] = 1.373e5;
      ocp_para[11] = -2.129e4;
      ocp_para[12] = -1.722e5;
      ocp_para[13] = 3.956e4;
      ocp_para[14] = 9.302e4;
      ocp_para[15] = -3.280e4;
      container->Add("OCP_PARA",ocp_para);

      // initialize parameter class for anode material
      parameters_anode_ = Teuchos::rcp(new MAT::PAR::Electrode(container));

      // initialize anode material
      anode_ = Teuchos::rcp(new MAT::Electrode(parameters_anode_.get()));

      // define sample concentration values for anode material
      // cf. Goldin et al., Electrochimica Acta 64 (2012) 118-129
      concentrations_anode_.resize(3,0.);
      concentrations_anode_[0] = 2.029;   // anode concentration at 0% state of charge
      concentrations_anode_[1] = 8.05;    // anode concentration at 50% intercalation fraction (constitutes singularity in original form of Redlich-Kister expansion)
      concentrations_anode_[2] = 10.88;   // anode concentration at 100% state of charge

      // add cathode parameters to container according to Goldin et al., Electrochimica Acta 64 (2012) 118-129
      container->Add("C_MAX",23.9);
      container->Add("OCP_PARA_NUM",21);
      ocp_para.resize(21,0.);
      ocp_para[0] = 3.954616e5;
      ocp_para[1] = -7.676e4;
      ocp_para[2] = 3.799e4;
      ocp_para[3] = -2.873e4;
      ocp_para[4] = 1.169e4;
      ocp_para[5] = 1.451e4;
      ocp_para[6] = -8.938e4;
      ocp_para[7] = 1.671e5;
      ocp_para[8] = -7.236e4;
      ocp_para[9] = -1.746e5;
      ocp_para[10] = -4.067e5;
      ocp_para[11] = 9.534e5;
      ocp_para[12] = 5.897e5;
      ocp_para[13] = -7.455e5;
      ocp_para[14] = -1.102e6;
      ocp_para[15] = -2.927e5;
      ocp_para[16] = 7.214e5;
      ocp_para[17] = 9.029e5;
      ocp_para[18] = -1.599e5;
      ocp_para[19] = 6.658e5;
      ocp_para[20] = -1.084e6;
      container->Add("OCP_PARA",ocp_para);

      // initialize parameter class for cathode material
      parameters_cathode_ = Teuchos::rcp(new MAT::PAR::Electrode(container));

      // initialize cathode material
      cathode_ = Teuchos::rcp(new MAT::Electrode(parameters_cathode_.get()));

      // define sample concentration values for cathode material
      // cf. Goldin et al., Electrochimica Acta 64 (2012) 118-129
      concentrations_cathode_.resize(3,0.);
      concentrations_cathode_[0] = 10.56;   // cathode concentration at 100% state of charge
      concentrations_cathode_[1] = 11.95;   // cathode concentration at 50% intercalation fraction (constitutes singularity in original form of Redlich-Kister expansion)
      concentrations_cathode_[2] = 22.37;   // cathode concentration at 0% state of charge

      // calculate factor F/(RT) for temperature 298.0 K
      frt_ = INPAR::ELCH::faraday_const/(INPAR::ELCH::gas_const*298.);

      return;
    }

    //! tear down unit testing for electrode material
    void tearDown()
    {
      // free memory by resetting pointers and clearing vectors
      anode_ = Teuchos::null;
      cathode_ = Teuchos::null;
      parameters_anode_ = Teuchos::null;
      parameters_cathode_ = Teuchos::null;
      concentrations_anode_.clear();
      concentrations_cathode_.clear();
      frt_ = 0.;

      return;
    }

    //! test member function ComputeOpenCircuitPotential of electrode material
    void test_ComputeOpenCircuitPotential()
    {
      // define results and tolerances for anode material
      const double results_anode[] = {1.78149102529067321354e-01,1.00269118707846299765e-01,8.08307883078566852264e-02};
      const double tolerances_anode[] = {1.e-17,1.e-17,1.e-18};

      // test member function using sample concentration values for anode material
      for(unsigned i=0; i<concentrations_anode_.size(); ++i)
        TS_ASSERT_DELTA(anode_->ComputeOpenCircuitPotential(concentrations_anode_[i],INPAR::ELCH::faraday_const,frt_),results_anode[i],tolerances_anode[i]);

      // define results and tolerances for cathode material
      const double results_cathode [] = {3.97322289583641552468e+00,3.90180104449214848472e+00,3.55923291247763096123e+00};
      const double tolerances_cathode[] = {1.e-16,1.e-16,1.e-16};

      // test member function using sample concentration values for cathode material
      for(unsigned i=0; i<concentrations_cathode_.size(); ++i)
        TS_ASSERT_DELTA(cathode_->ComputeOpenCircuitPotential(concentrations_cathode_[i],INPAR::ELCH::faraday_const,frt_),results_cathode[i],tolerances_cathode[i]);

      return;
    }

    //! test member function ComputeFirstDerivOpenCircuitPotential of electrode material
    void test_ComputeFirstDerivOpenCircuitPotential()
    {
      // define results and tolerances for anode material
      const double results_anode[] = {-4.67265915814377499893e-02,-4.69987444027383099998e-03,-8.94809228199507609991e-03};
      const double tolerances_anode[] = {1.e-18,1.e-19,1.e-19};

      // test member function using sample concentration values for anode material
      for(unsigned i=0; i<concentrations_anode_.size(); ++i)
        TS_ASSERT_DELTA(anode_->ComputeFirstDerivOpenCircuitPotential(concentrations_anode_[i],INPAR::ELCH::faraday_const,frt_),results_anode[i],tolerances_anode[i]);

      // define results and tolerances for cathode material
      const double results_cathode [] = {-5.87374545823974442138e-02,-4.59544093325152322449e-02,-2.40034452721541308184e-02};
      const double tolerances_cathode[] = {1.e-18,1.e-18,1.e-18};

      // test member function using sample concentration values for cathode material
      for(unsigned i=0; i<concentrations_cathode_.size(); ++i)
        TS_ASSERT_DELTA(cathode_->ComputeFirstDerivOpenCircuitPotential(concentrations_cathode_[i],INPAR::ELCH::faraday_const,frt_),results_cathode[i],tolerances_cathode[i]);

      return;
    }

    //! test member function ComputeSecondDerivOpenCircuitPotential of electrode material
    void test_ComputeSecondDerivOpenCircuitPotential()
    {
      // define results and tolerances for anode material
      const double results_anode[] = {3.04927774052812655681e-02,-1.05174077799887171619e-03,-7.35165984480485820186e-04};
      const double tolerances_anode[] = {1.e-18,1.e-19,1.e-20};

      // test member function using sample concentration values for anode material
      for(unsigned i=0; i<concentrations_anode_.size(); ++i)
        TS_ASSERT_DELTA(anode_->ComputeSecondDerivOpenCircuitPotential(concentrations_anode_[i],INPAR::ELCH::faraday_const,frt_),results_anode[i],tolerances_anode[i]);

      // define results and tolerances for cathode material
      const double results_cathode [] = {1.41256487214221713306e-02,5.72637600860917590773e-03,-7.56499159049272597299e-03};
      const double tolerances_cathode[] = {1.e-18,1.e-19,1.e-19};

      // test member function using sample concentration values for cathode material
      for(unsigned i=0; i<concentrations_cathode_.size(); ++i)
        TS_ASSERT_DELTA(cathode_->ComputeSecondDerivOpenCircuitPotential(concentrations_cathode_[i],INPAR::ELCH::faraday_const,frt_),results_cathode[i],tolerances_cathode[i]);

      return;
    }

  private:

    //! anode material
    Teuchos::RCP<const MAT::Electrode> anode_;

    //! cathode material
    Teuchos::RCP<const MAT::Electrode> cathode_;

    //! parameters for anode material
    Teuchos::RCP<MAT::PAR::Electrode> parameters_anode_;

    //! parameters for cathode material
    Teuchos::RCP<MAT::PAR::Electrode> parameters_cathode_;

    //! sample concentration values for anode material
    std::vector<double> concentrations_anode_;

    //! sample concentration values for cathode material
    std::vector<double> concentrations_cathode_;

    //! factor F/(RT)
    double frt_;
};
#endif
