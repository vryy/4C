/*----------------------------------------------------------------------*/
/*! \file
 \brief material defining the degree of freedom of phase within a multiphase porous fluid.

   \level 3

 *----------------------------------------------------------------------*/

#ifndef FOUR_C_MAT_FLUIDPORO_SINGLEPHASEDOF_HPP
#define FOUR_C_MAT_FLUIDPORO_SINGLEPHASEDOF_HPP


#include "4C_config.hpp"

#include "4C_inpar_material.hpp"
#include "4C_linalg_serialdensematrix.hpp"
#include "4C_material_parameter_base.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Mat
{
  namespace PAR
  {
    // forward declaration
    class FluidPoroPhaseLaw;

    //! interface class for generic phase degree of freedom
    class FluidPoroPhaseDof : public Core::Mat::PAR::Parameter
    {
     public:
      /// standard constructor
      explicit FluidPoroPhaseDof(const Core::Mat::PAR::Parameter::Data& matdata);

      //! build the phase dof
      static FluidPoroPhaseDof* create_phase_dof(int phasedofId);

      /// create material instance of matching type with my parameters
      Teuchos::RCP<Core::Mat::Material> create_material() override { return Teuchos::null; };

      /// initialize
      virtual void initialize() = 0;

      /// return type of phase law
      virtual Core::Materials::MaterialType poro_phase_law_type() const = 0;

      /// mark dofs associated with this phase in a given row (=numphase) in a matrix
      virtual void fill_do_f_matrix(
          Core::LinAlg::SerialDenseMatrix& dofmat, int numphase) const = 0;

      /// evaluate saturation of the phase
      virtual double evaluate_saturation(int phasenum, const std::vector<double>& state,
          const std::vector<double>& pressure) const = 0;

      /// evaluate the generalized(!) pressure of this phase
      virtual double evaluate_gen_pressure(
          int phasenum, const std::vector<double>& state) const = 0;

      //! evaluate derivative of saturation with respect to pressure
      virtual double evaluate_deriv_of_saturation_wrt_pressure(
          int phasenum, int doftoderive, const std::vector<double>& pressure) const = 0;

      //! evaluate 2nd derivative of saturation with respect to pressure
      virtual double evaluate_second_deriv_of_saturation_wrt_pressure(int phasenum,
          int firstdoftoderive, int seconddoftoderive,
          const std::vector<double>& pressure) const = 0;

      //! evaluate derivative of degree of freedom with respect to pressure
      virtual double evaluate_deriv_of_dof_wrt_pressure(
          int phasenum, int doftoderive, const std::vector<double>& state) const = 0;
    };

    //! pressure degree of freedom of a single poro fluid phase
    class FluidPoroPhaseDofPressure : public FluidPoroPhaseDof
    {
     public:
      /// standard constructor
      explicit FluidPoroPhaseDofPressure(const Core::Mat::PAR::Parameter::Data& matdata);

      /// initialize
      void initialize() override;

      /// return type of phase law
      Core::Materials::MaterialType poro_phase_law_type() const override;

      /// mark dofs associated with this phase in a given row (=numphase) in a matrix
      void fill_do_f_matrix(Core::LinAlg::SerialDenseMatrix& dofmat, int numphase) const override;

      /// evaluate saturation of the phase
      double evaluate_saturation(int phasenum, const std::vector<double>& state,
          const std::vector<double>& pressure) const override;

      /// evaluate the generalized(!) pressure of this phase
      double evaluate_gen_pressure(int phasenum, const std::vector<double>& state) const override;

      //! evaluate derivative of saturation with respect to pressure
      double evaluate_deriv_of_saturation_wrt_pressure(
          int phasenum, int doftoderive, const std::vector<double>& pressure) const override;

      //! evaluate 2nd derivative of saturation with respect to pressure
      double evaluate_second_deriv_of_saturation_wrt_pressure(int phasenum, int firstdoftoderive,
          int seconddoftoderive, const std::vector<double>& pressure) const override;

      //! evaluate derivative of degree of freedom with respect to pressure
      double evaluate_deriv_of_dof_wrt_pressure(
          int phasenum, int doftoderive, const std::vector<double>& state) const override;

     protected:
      /// ID of pressure-saturation law
      const int phaselawId_;
      /// implementation of pressure-saturation law
      FluidPoroPhaseLaw* phaselaw_;
    };

    //! differential pressure degree of freedom of a single phase of porous multiphase fluid
    class FluidPoroPhaseDofDiffPressure : public FluidPoroPhaseDof
    {
     public:
      /// standard constructor
      explicit FluidPoroPhaseDofDiffPressure(const Core::Mat::PAR::Parameter::Data& matdata);

      /// initialize
      void initialize() override;

      /// return type of phase law
      Core::Materials::MaterialType poro_phase_law_type() const override;

      /// mark dofs associated with this phase in a given row (=numphase) in a matrix
      void fill_do_f_matrix(Core::LinAlg::SerialDenseMatrix& dofmat, int numphase) const override;

      /// evaluate saturation of the phase
      double evaluate_saturation(int phasenum, const std::vector<double>& state,
          const std::vector<double>& pressure) const override;

      /// evaluate the generalized(!) pressure of this phase
      double evaluate_gen_pressure(int phasenum, const std::vector<double>& state) const override;

      //! evaluate derivative of saturation with respect to pressure
      double evaluate_deriv_of_saturation_wrt_pressure(
          int phasenum, int doftoderive, const std::vector<double>& pressure) const override;

      //! evaluate 2nd derivative of saturation with respect to pressure
      double evaluate_second_deriv_of_saturation_wrt_pressure(int phasenum, int firstdoftoderive,
          int seconddoftoderive, const std::vector<double>& pressure) const override;

      //! evaluate derivative of degree of freedom with respect to pressure
      double evaluate_deriv_of_dof_wrt_pressure(
          int phasenum, int doftoderive, const std::vector<double>& state) const override;

     protected:
      ///  pressure-coefficients defining the differential pressure
      const std::vector<int> diffpresCoeffs_;
      /// ID of pressure-saturation law
      const int phaselawId_;
      /// implementation of pressure-saturation law
      FluidPoroPhaseLaw* phaselaw_;
    };

    //! saturation degree of freedom of a single phase within a multiphase porous fluid
    class FluidPoroPhaseDofSaturation : public FluidPoroPhaseDof
    {
     public:
      /// standard constructor
      explicit FluidPoroPhaseDofSaturation(const Core::Mat::PAR::Parameter::Data& matdata);

      /// initialize
      void initialize() override;

      /// return type of phase law
      Core::Materials::MaterialType poro_phase_law_type() const override;

      /// mark dofs associated with this phase in a given row (=numphase) in a matrix
      void fill_do_f_matrix(Core::LinAlg::SerialDenseMatrix& dofmat, int numphase) const override;

      /// evaluate saturation of the phase
      double evaluate_saturation(int phasenum, const std::vector<double>& state,
          const std::vector<double>& pressure) const override;

      /// evaluate the generalized(!) pressure of this phase
      double evaluate_gen_pressure(int phasenum, const std::vector<double>& state) const override;

      //! evaluate derivative of saturation with respect to pressure
      double evaluate_deriv_of_saturation_wrt_pressure(
          int phasenum, int doftoderive, const std::vector<double>& pressure) const override;

      //! evaluate 2nd derivative of saturation with respect to pressure
      double evaluate_second_deriv_of_saturation_wrt_pressure(int phasenum, int firstdoftoderive,
          int seconddoftoderive, const std::vector<double>& pressure) const override;

      //! evaluate derivative of degree of freedom with respect to pressure
      double evaluate_deriv_of_dof_wrt_pressure(
          int phasenum, int doftoderive, const std::vector<double>& pressure) const override;

     protected:
      /// ID of pressure-saturation law
      const int phaselawId_;
      /// implementation of pressure-saturation law
      FluidPoroPhaseLaw* phaselaw_;
    };
  }  // namespace PAR
}  // namespace Mat

FOUR_C_NAMESPACE_CLOSE

#endif
