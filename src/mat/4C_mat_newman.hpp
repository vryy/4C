// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_MAT_NEWMAN_HPP
#define FOUR_C_MAT_NEWMAN_HPP

#include "4C_config.hpp"

#include "4C_comm_parobjectfactory.hpp"
#include "4C_mat_elchsinglemat.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Mat
{
  namespace PAR
  {
    /*----------------------------------------------------------------------*/
    /// material parameters for convection-diffusion
    class Newman : public ElchSingleMat
    {
     public:
      /// standard constructor
      explicit Newman(const Core::Mat::PAR::Parameter::Data& matdata);


      //! gets the function defining the concentration scaling of the transference number from the
      //! global problem for the first call and returns it for all later calls
      [[nodiscard]] const Core::Utils::FunctionOfTime&
      transference_number_concentration_scaling_funct();

      [[nodiscard]] bool has_transference_number_concentration_scaling() const
      {
        return transference_number_concentration_scaling_funct_num_.has_value();
      }

      //! gets the function defining the concentration scaling of the thermodynamic factor from the
      //! global problem for the first call and returns it for all later calls
      [[nodiscard]] const Core::Utils::FunctionOfTime&
      thermodynamic_factor_concentration_scaling_funct();

      [[nodiscard]] bool has_thermodynamic_factor_concentration_scaling() const
      {
        return thermodynamic_factor_concentration_scaling_funct_num_.has_value();
      }

      /// @name material parameters
      ///@{
      /// valence (= charge number)
      const double valence_;

      //! transference number
      const double transference_number_;

      //! function number defining the optional concentration scaling of the transference number
      const std::optional<int> transference_number_concentration_scaling_funct_num_;
      //! function defining the optional concentration scaling of the transference number
      std::optional<std::reference_wrapper<const Core::Utils::FunctionOfTime>>
          transference_number_concentration_scaling_funct_{std::nullopt};


      //! thermodynamic factor
      const double thermodynamic_factor_;

      //! function number defining the optional concentration scaling of the thermodynamic factor
      const std::optional<int> thermodynamic_factor_concentration_scaling_funct_num_;
      //! function defining the optional concentration scaling of the thermodynamic factor
      std::optional<std::reference_wrapper<const Core::Utils::FunctionOfTime>>
          thermodynamic_factor_concentration_scaling_funct_{std::nullopt};
      ///@}

      /// create material instance of matching type with my parameters
      std::shared_ptr<Core::Mat::Material> create_material() override;
    };  // class Newman

  }  // namespace PAR

  class NewmanType : public Core::Communication::ParObjectType
  {
   public:
    std::string name() const override { return "NewmanType"; }

    static NewmanType& instance() { return instance_; };

    Core::Communication::ParObject* create(Core::Communication::UnpackBuffer& buffer) override;

   private:
    static NewmanType instance_;
  };

  /*----------------------------------------------------------------------*/
  /// Wrapper for the material properties of an ion species in an electrolyte solution
  class Newman : public ElchSingleMat
  {
   public:
    Newman() = default;
    /// construct the material object given material parameters
    explicit Newman(Mat::PAR::Newman* params);

    //! @name Packing and Unpacking

    /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of parobject.H (this file) and should return it in this method.
    */
    int unique_par_object_id() const override
    {
      return NewmanType::instance().unique_par_object_id();
    }

    /*!
      \brief Pack this class so it can be communicated

      Resizes the vector data and stores all information of a class in it.
      The first information to be stored in data has to be the
      unique parobject id delivered by unique_par_object_id() which will then
      identify the exact class on the receiving processor.

      \param data (in/out): char vector to store class information
    */
    void pack(Core::Communication::PackBuffer& data) const override;

    /*!
      \brief Unpack data from a char vector into this class

      The vector data contains all information to rebuild the
      exact copy of an instance of a class on a different processor.
      The first entry in data has to be an integer which is the unique
      parobject id defined at the top of this file and delivered by
      unique_par_object_id().

      \param data (in) : vector storing all data to be unpacked into this
      instance.
    */
    void unpack(Core::Communication::UnpackBuffer& buffer) override;

    //@}

    /// material type
    Core::Materials::MaterialType material_type() const override
    {
      return Core::Materials::m_newman;
    }

    /// return copy of this material object
    std::shared_ptr<Core::Mat::Material> clone() const override
    {
      return std::make_shared<Newman>(*this);
    }

    /// valence (= charge number)
    double valence() const { return params_->valence_; }

    /// computation of the transference number
    [[nodiscard]] double compute_transference_number(double concentration) const;
    /// computation of the first derivative of the transference number w.r.t. concentration
    [[nodiscard]] double compute_concentration_derivative_of_transference_number(
        double concentration) const;

    /// computation of the thermodynamic factor
    [[nodiscard]] double compute_thermodynamic_factor(double concentration) const;
    /// computation of the first derivative of the thermodynamic factor w.r.t. concentration
    [[nodiscard]] double compute_concentration_derivative_of_thermodynamic_factor(
        double concentration) const;

    [[nodiscard]] Core::Mat::PAR::Parameter* parameter() const override { return params_; }

   private:
    /// my material parameters
    Mat::PAR::Newman* params_;
  };
}  // namespace Mat

FOUR_C_NAMESPACE_CLOSE

#endif
