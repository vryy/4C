// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_MAT_VISCOELAST_GENERALIZEDMAXWELL_HPP
#define FOUR_C_MAT_VISCOELAST_GENERALIZEDMAXWELL_HPP

#include "4C_config.hpp"

#include "4C_linalg_fixedsizematrix.hpp"
#include "4C_mat_elast_summand.hpp"
#include "4C_material_parameter_base.hpp"

#include <functional>
#include <vector>

FOUR_C_NAMESPACE_OPEN

namespace Mat
{
  namespace Elastic
  {
    namespace PAR
    {
      /*!
       * @brief material parameters for viscous contribution to a viscoelastic branch of a
       * generalized Maxwell model
       *
       * <h3>Input line</h3>
       *  MAT 1 VISCO_GeneralizedMaxwell NUMBRANCH 3 MATIDS 4 5 6 SOLVE
       *  ExponentialTimeDiscretization
       *  MAT 1 VISCO_GeneralizedMaxwell NUMBRANCH 3 MATIDS 4 5 6 SOLVE OneStepTheta
       */
      class GeneralizedMaxwell : public Core::Mat::PAR::Parameter
      {
       public:
        /// standard constructor
        GeneralizedMaxwell(const Core::Mat::PAR::Parameter::Data& matdata);

        /// @name material parameters
        //@{

        /// material parameters
        int numbranch_;
        const std::vector<int> matids_;
        std::string solve_;
        //@}

        /// create material instance of matching type with my parameters

        std::shared_ptr<Core::Mat::Material> create_material() override { return nullptr; };
      };  // class GeneralizedMaxwell


      /*!
       * @brief material parameters for viscous contribution to a viscoelastic branch of a
       * generalized Maxwell model
       *
       * <h3>Input line</h3>
       * MAT 1 VISCO_GeneralizedMaxwellBranch TAU 1.5 MATID 2
       */
      class ViscoBranch : public Core::Mat::PAR::Parameter
      {
       public:
        /// standard constructor
        ViscoBranch(const Core::Mat::PAR::Parameter::Data& matdata);

        /// @name material parameters
        //@{

        /// material parameters

        double tau_;
        int matid_;


        //@}

        /// Override this method and throw error, as the material should be created in within the
        /// Factory method of the elastic summand
        std::shared_ptr<Core::Mat::Material> create_material() override
        {
          FOUR_C_THROW(
              "Cannot create a material from this method, as it should be created in "
              "Mat::Elastic::Summand::Factory.");
          return nullptr;
        };
      };  // class ViscoBranch

    }  // namespace PAR



    class GeneralizedMaxwell : public Summand
    {
     public:
      /// constructor with given material parameters
      GeneralizedMaxwell(Mat::Elastic::PAR::GeneralizedMaxwell* params);

      /// @name Access material constants
      //@{

      /// material type
      Core::Materials::MaterialType material_type() const override
      {
        return Core::Materials::mes_generalizedmaxwell;
      }
      //@}

      /// Read material parameters
      void read_material_parameters(int& numbranch,  ///< number of viscoelastic branches
          const std::vector<int>*& matids,           ///< material IDs of the viscoelastic branches
          std::string& solve  /// variant of the solution of the evolution integral
          ) override;

      /// @name Access methods
      //@{
      const std::vector<std::vector<std::shared_ptr<Mat::Elastic::Summand>>>& get_branchespotsum()
          const
      {
        return branchespotsum_;
      }

      const std::vector<double>& get_branchtaus() const { return branchtau_; }
      //@}

      /// Indicator for formulation
      void specify_formulation(
          bool& isoprinc,     ///< global indicator for isotropic principal formulation
          bool& isomod,       ///< global indicator for isotropic split formulation
          bool& anisoprinc,   ///< global indicator for anisotropic principal formulation
          bool& anisomod,     ///< global indicator for anisotropic split formulation
          bool& viscogeneral  ///< general indicator, if one viscoelastic formulation is used
          ) override
      {
        viscogeneral = true;
        return;
      };

      /// Indicator for the chosen viscoelastic formulations
      void specify_visco_formulation(
          bool& visco_iso_rate,  ///< global indicator for isotropic rate-dependent visco response
          bool& visco_generalized_maxwell,  ///< global indicator for generalized Maxwell model
          bool& visco_fsls                  ///< global indicator for FSLS model
          ) override
      {
        visco_generalized_maxwell = true;
        return;
      };


     private:
      /// my material parameters
      Mat::Elastic::PAR::GeneralizedMaxwell* params_;

     protected:
      /// summands of the GeneralizedMaxwell material or each branch
      std::vector<std::vector<std::shared_ptr<Mat::Elastic::Summand>>> branchespotsum_;
      /// branch relaxation parameters
      std::vector<double> branchtau_;
      /// summands in one particular branch
      std::vector<std::shared_ptr<Mat::Elastic::Summand>> internalpotsum_;
    };

    class ViscoBranch : public Summand
    {
     public:
      /// constructor with given material parameters
      ViscoBranch(Mat::Elastic::PAR::ViscoBranch* params);

      /// @name Access material constants
      //@{

      /// material type
      Core::Materials::MaterialType material_type() const override
      {
        return Core::Materials::mes_viscobranch;
      }

      //@}

      /// Read material parameters
      void read_material_parameters(double& tau,  ///< branch relaxation time
          int& matid                              ///< material ID of branch elasticity
          ) override;

      /// Indicator for formulation
      void specify_formulation(
          bool& isoprinc,     ///< global indicator for isotropic principal formulation
          bool& isomod,       ///< global indicator for isotropic split formulation
          bool& anisoprinc,   ///< global indicator for anisotropic principal formulation
          bool& anisomod,     ///< global indicator for anisotropic split formulation
          bool& viscogeneral  ///< general indicator, if one viscoelastic formulation is used
          ) override
      {
        viscogeneral = true;
        return;
      };


     private:
      /// my material parameters
      Mat::Elastic::PAR::ViscoBranch* params_;

    };  // class ViscoBranch


  }  // namespace Elastic


  namespace ViscoElast
  {
    namespace Kernels
    {
      enum class GeneralizedMaxwellSolveKind
      {
        one_step_theta,
        exponential_time_discretization
      };

      using StressVector = Core::LinAlg::Matrix<6, 1>;
      using TangentMatrix = Core::LinAlg::Matrix<6, 6>;
      using PointHistory = std::vector<StressVector>;

      struct GeneralizedMaxwellKernelInput
      {
        int visco_mat_id = -1;
        int gp = -1;
        int ele_gid = -1;
        double dt = 0.0;
        double one_step_theta = 0.5;
        GeneralizedMaxwellSolveKind solve_kind =
            GeneralizedMaxwellSolveKind::exponential_time_discretization;
        const PointHistory* previous_branch_elastic_stress = nullptr;
        const PointHistory* previous_branch_stress = nullptr;
      };

      using BranchResponseEvaluator = std::function<void(
          int branch_index, StressVector& branch_elastic_stress, TangentMatrix& branch_cmat)>;

      void evaluate_generalized_maxwell_kernel(StressVector& q_total, TangentMatrix& cmatq_total,
          PointHistory& current_branch_elastic_stress, PointHistory& current_branch_stress,
          const std::vector<double>& branch_taus, const GeneralizedMaxwellKernelInput& input,
          const BranchResponseEvaluator& evaluate_branch_response);
    }  // namespace Kernels
  }  // namespace ViscoElast
}  // namespace Mat

FOUR_C_NAMESPACE_CLOSE

#endif
