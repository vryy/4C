// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_MAT_VISCOELASTHYPER_HPP
#define FOUR_C_MAT_VISCOELASTHYPER_HPP



#include "4C_config.hpp"

#include "4C_comm_parobjectfactory.hpp"
#include "4C_mat_elasthyper.hpp"
#include "4C_mat_material_factory.hpp"
#include "4C_mat_so3_material.hpp"
#include "4C_mat_viscoelast_contribution.hpp"
#include "4C_mat_viscoelast_state.hpp"
#include "4C_mat_viscoelast_summand.hpp"
#include "4C_material_parameter_base.hpp"

#include <array>
#include <memory>
#include <optional>
#include <string>
#include <vector>

FOUR_C_NAMESPACE_OPEN

// forward declaration due to avoid header definition
namespace Mat
{
  namespace Elastic
  {
    class Summand;
  }

  namespace ViscoElast
  {
    class Summand;
  }

  // forward declaration
  class ViscoElastHyper;

  namespace PAR
  {
    /*----------------------------------------------------------------------*/
    /// Collection of viscohyperelastic materials
    ///
    /// Storage map of hyperelastic summands.


    class ViscoElastHyper : public Mat::PAR::ElastHyper
    //    class ViscoElastHyper : public Core::Mat::PAR::Parameter
    {
      friend class Mat::ViscoElastHyper;

     private:
      struct SummandSplit
      {
        int numelast = 0;
        std::vector<int> elast_matids;
        int numvisco = 0;
        std::vector<int> visco_matids;
        bool uses_legacy_matids = true;
      };

      static SummandSplit parse_summand_split(const Core::Mat::PAR::Parameter::Data& matdata);
      ViscoElastHyper(
          const Core::Mat::PAR::Parameter::Data& matdata, const SummandSplit& summand_split);

     public:
      /// standard constructor
      ///
      /// This constructor recursively calls the constructors of the
      /// parameter sets of the hyperelastic summands.
      ViscoElastHyper(const Core::Mat::PAR::Parameter::Data& matdata);

      /// create material instance of matching type with my parameters
      std::shared_ptr<Core::Mat::Material> create_material() override;

      /// number of elastic summands in split-mode input
      const int numelast_;

      /// list of elastic summand material IDs in split-mode input
      const std::vector<int> elast_matids_;

      /// number of visco summands in split-mode input
      const int numvisco_;

      /// list of visco summand material IDs in split-mode input
      const std::vector<int> visco_matids_;

      /// true if legacy NUMMAT/MATIDS input was used and auto-partitioned
      const bool uses_legacy_matids_;

      //@}

    };  // class ViscoElastHyper

  }  // namespace PAR

  class ViscoElastHyperType : public Core::Communication::ParObjectType
  {
   public:
    std::string name() const override { return "ViscoElastHyperType"; }

    static ViscoElastHyperType& instance() { return instance_; };

    Core::Communication::ParObject* create(Core::Communication::UnpackBuffer& buffer) override;

   private:
    static ViscoElastHyperType instance_;
  };


  /*----------------------------------------------------------------------*/
  /// Collection of hyperelastic materials
  ///
  /// This collection offers to additively compose a stress response
  /// based on summands defined separately.  This is possible, because
  /// we deal with hyperelastic materials, which are composed
  /// of (Helmholtz free energy density) potentials.  Effectively, we want
  ///\f[
  ///  \Psi(\boldsymbol{C}) = \sum_i \Psi_i(\boldsymbol{C})
  ///\f]
  /// in which the individual \f$\Psi_i\f$ is implemented as #Mat::Elastic::Summand.
  ///
  /// Quite often the right Cauchy-Green 2-tensor \f$\boldsymbol{C}\f$
  /// is replaced by its various invariant forms as argument.
  ///
  /// The task of ElastHyper is the evaluation of the
  /// potential energies and their derivatives to obtain the actual
  /// stress response and the elasticity tensor. The storage is located
  /// at the associated member #params_.
  ///
  /// <h3>References</h3>
  /// <ul>
  /// <li> [1] GA Holzapfel, "Nonlinear solid mechanics", Wiley, 2000.
  /// </ul>
  ///

  class Material;

  class ViscoElastHyper : public Mat::ElastHyper
  {
   public:
    /// construct empty material object
    ViscoElastHyper();

    /// construct the material object given material parameters
    explicit ViscoElastHyper(Mat::PAR::ViscoElastHyper* params);

    ///@name Packing and Unpacking
    //@{

    /// \brief Return unique ParObject id
    ///
    /// every class implementing ParObject needs a unique id defined at the
    /// top of parobject.H (this file) and should return it in this method.
    int unique_par_object_id() const override
    {
      return ViscoElastHyperType::instance().unique_par_object_id();
    }

    /// \brief Pack this class so it can be communicated
    ///
    /// Resizes the vector data and stores all information of a class in it.
    /// The first information to be stored in data has to be the
    /// unique parobject id delivered by unique_par_object_id() which will then
    /// identify the exact class on the receiving processor.
    ///
    /// \param data (in/out): char vector to store class information
    void pack(Core::Communication::PackBuffer& data) const override;

    /// \brief Unpack data from a char vector into this class
    ///
    /// The vector data contains all information to rebuild the
    /// exact copy of an instance of a class on a different processor.
    /// The first entry in data has to be an integer which is the unique
    /// parobject id defined at the top of this file and delivered by
    /// unique_par_object_id().
    ///
    /// \param data (in) : vector storing all data to be unpacked into this
    ///                    instance.
    void unpack(Core::Communication::UnpackBuffer& buffer) override;

    //@}

    /// material type
    Core::Materials::MaterialType material_type() const override
    {
      return Core::Materials::m_viscoelasthyper;
    }

    /// check if element kinematics and material kinematics are compatible
    void valid_kinematics(Inpar::Solid::KinemType kinem) override
    {
      FOUR_C_ASSERT_ALWAYS(kinem == Inpar::Solid::KinemType::nonlinearTotLag,
          "element and material kinematics are not compatible");
    }

    /// return copy of this material object
    std::shared_ptr<Core::Mat::Material> clone() const override
    {
      return std::make_shared<ViscoElastHyper>(*this);
    }

    /// Check if history variables are already initialized
    virtual bool initialized() const { return state_.initialized(); }

    /// hyperelastic stress response plus elasticity tensor
    void evaluate(const Core::LinAlg::Tensor<double, 3, 3>* defgrad,
        const Core::LinAlg::SymmetricTensor<double, 3, 3>& glstrain,
        const Teuchos::ParameterList& params, const EvaluationContext<3>& context,
        Core::LinAlg::SymmetricTensor<double, 3, 3>& stress,
        Core::LinAlg::SymmetricTensor<double, 3, 3, 3, 3>& cmat, int gp,
        int eleGID) override;  ///< Constitutive matrix

    /// setup material description
    void setup(int numgp, const Discret::Elements::Fibers& fibers,
        const std::optional<Discret::Elements::CoordinateSystem>& coord_system) override;

    /// update history variables
    void update() override;

   protected:
    /// @name Flags to specify the viscous formulations
    //@{
    bool isovisco_;                   ///< global indicator for isotropic split viscous formulation
    bool visco_generalized_maxwell_;  ///< global indicator for viscous contribution of branches
                                      ///< according to the generalized Maxwell model
    bool visco_fsls_;  ///< global indicator for viscous contribution according the FSLS model
                       //@}

   private:
    struct EvaluateWorkspace
    {
      EvaluateWorkspace();

      Core::LinAlg::Matrix<6, 1> glstrain_mat;
      Core::LinAlg::Matrix<6, 1> c_stress;
      Core::LinAlg::Matrix<6, 1> i_c_stress;
      Core::LinAlg::Matrix<6, 1> c_strain;
      Core::LinAlg::Matrix<6, 1> mod_c_strain;
      Core::LinAlg::Matrix<6, 1> id2;
      Core::LinAlg::Matrix<6, 6> id4;
      Core::LinAlg::Matrix<6, 6> id4sharp;
      Core::LinAlg::Matrix<3, 1> prinv;
      Core::LinAlg::Matrix<3, 1> modinv;
      Core::LinAlg::Matrix<7, 1> rateinv;
      Core::LinAlg::Matrix<7, 1> modrateinv;
      Core::LinAlg::Matrix<3, 1> dPI;
      Core::LinAlg::Matrix<6, 1> ddPII;
      Core::LinAlg::Matrix<6, 1> scgrate;
      Core::LinAlg::Matrix<6, 1> modrcgrate;
      Core::LinAlg::Matrix<8, 1> mu;
      Core::LinAlg::Matrix<8, 1> modmu;
      Core::LinAlg::Matrix<33, 1> xi;
      Core::LinAlg::Matrix<33, 1> modxi;
      Core::LinAlg::SymmetricTensor<double, 3, 3> c{};
      Core::LinAlg::SymmetricTensor<double, 3, 3> i_c{};
      double dt = 0.0;
    };

    using ViscoModelKind = ViscoElast::ViscoModelKind;
    using ActiveModelSequence = std::vector<ViscoModelKind>;

    [[nodiscard]] static constexpr std::array<ViscoModelKind, 3> visco_model_registry()
    {
      return {ViscoModelKind::iso_rate, ViscoModelKind::generalized_maxwell, ViscoModelKind::fsls};
    }

    [[nodiscard]] static bool is_visco_material_type(Core::Materials::MaterialType material_type);
    [[nodiscard]] const Mat::PAR::ViscoElastHyper* visco_params() const;
    [[nodiscard]] int visco_mat_id(unsigned index) const;
    void rebuild_summand_sets();
    void rebuild_effective_summand_properties();

    [[nodiscard]] bool is_model_flag_enabled(ViscoModelKind model_kind) const;
    [[nodiscard]] bool is_model_active(ViscoModelKind model_kind) const;
    void rebuild_active_model_sequence();
    void rebuild_contributions();
    void setup_contributions(int gp, int eleGID);

    [[nodiscard]] ViscoElastState::ActiveModels active_models_from_flags() const;
    [[nodiscard]] ViscoElastState::ActiveModels active_models_from_sequence() const;
    void ensure_model_activation_consistency(const char* context) const;

    [[nodiscard]] ViscoElastState::ActiveModels active_models() const;
    [[nodiscard]] std::size_t history_entry_count_for_setup(ViscoModelKind model_kind) const;
    [[nodiscard]] unsigned int history_capacity_for_update(ViscoModelKind model_kind) const;
    [[nodiscard]] double read_visco_time_step_size(
        const EvaluationContext<3>& context, int gp, int eleGID) const;

    void prepare_evaluate_kinematics(const Core::LinAlg::SymmetricTensor<double, 3, 3>& glstrain,
        const EvaluationContext<3>& context, EvaluateWorkspace& workspace, int gp,
        int eleGID) const;
    void initialize_elastic_response(const EvaluateWorkspace& workspace,
        Core::LinAlg::Matrix<6, 1>& stress_view, Core::LinAlg::Matrix<6, 6>& cmat_view,
        Core::LinAlg::SymmetricTensor<double, 3, 3>& stress,
        Core::LinAlg::SymmetricTensor<double, 3, 3, 3, 3>& cmat) const;
    void add_visco_contributions_in_sequence(const Teuchos::ParameterList& params,
        EvaluateWorkspace& workspace, Core::LinAlg::Matrix<6, 1>& stress_view,
        Core::LinAlg::Matrix<6, 6>& cmat_view, int gp, int eleGID);
    void add_post_elastic_composition_hooks(const Teuchos::ParameterList& params,
        const EvaluationContext<3>& context, const EvaluateWorkspace& workspace,
        Core::LinAlg::SymmetricTensor<double, 3, 3>& stress,
        Core::LinAlg::SymmetricTensor<double, 3, 3, 3, 3>& cmat, int gp, int eleGID);

    ActiveModelSequence active_model_sequence_;
    std::vector<std::shared_ptr<ViscoElast::Contribution>> contributions_;
    std::vector<std::shared_ptr<Mat::Elastic::Summand>> elast_potsum_;
    std::vector<std::shared_ptr<ViscoElast::Summand>> visco_potsum_;
    SummandProperties elast_summand_properties_;
    SummandProperties effective_summand_properties_;
    ViscoElastState state_;  ///< unified viscoelastic history state
  };  // class ViscoElastHyper

}  // namespace Mat

FOUR_C_NAMESPACE_CLOSE

#endif
