/*----------------------------------------------------------------------*/
/*! \file
 \brief material defining the degree of freedom of phase within a multiphase porous fluid.

   \level 3

 *----------------------------------------------------------------------*/

#ifndef BACI_MAT_FLUIDPORO_SINGLEPHASEDOF_H
#define BACI_MAT_FLUIDPORO_SINGLEPHASEDOF_H


#include "baci_config.H"

#include "baci_inpar_material.H"
#include "baci_linalg_serialdensematrix.H"
#include "baci_mat_par_parameter.H"

BACI_NAMESPACE_OPEN

namespace MAT
{
  namespace PAR
  {
    // forward declaration
    class FluidPoroPhaseLaw;

    //! interface class for generic phase degree of freedom
    class FluidPoroPhaseDof : public Parameter
    {
     public:
      /// standard constructor
      explicit FluidPoroPhaseDof(Teuchos::RCP<MAT::PAR::Material> matdata);

      //! build the phase dof
      static FluidPoroPhaseDof* CreatePhaseDof(int phasedofId);

      /// create material instance of matching type with my parameters
      Teuchos::RCP<MAT::Material> CreateMaterial() override { return Teuchos::null; };

      /// initialize
      virtual void Initialize() = 0;

      /// return type of phase law
      virtual INPAR::MAT::MaterialType PoroPhaseLawType() const = 0;

      /// mark dofs associated with this phase in a given row (=numphase) in a matrix
      virtual void FillDoFMatrix(CORE::LINALG::SerialDenseMatrix& dofmat, int numphase) const = 0;

      /// evaluate saturation of the phase
      virtual double EvaluateSaturation(int phasenum, const std::vector<double>& state,
          const std::vector<double>& pressure) const = 0;

      /// evaluate the generalized(!) pressure of this phase
      virtual double EvaluateGenPressure(int phasenum, const std::vector<double>& state) const = 0;

      //! evaluate derivative of saturation with respect to pressure
      virtual double EvaluateDerivOfSaturationWrtPressure(
          int phasenum, int doftoderive, const std::vector<double>& pressure) const = 0;

      //! evaluate 2nd derivative of saturation with respect to pressure
      virtual double EvaluateSecondDerivOfSaturationWrtPressure(int phasenum, int firstdoftoderive,
          int seconddoftoderive, const std::vector<double>& pressure) const = 0;

      //! evaluate derivative of degree of freedom with respect to pressure
      virtual double EvaluateDerivOfDofWrtPressure(
          int phasenum, int doftoderive, const std::vector<double>& state) const = 0;
    };

    //! pressure degree of freedom of a single poro fluid phase
    class FluidPoroPhaseDofPressure : public FluidPoroPhaseDof
    {
     public:
      /// standard constructor
      explicit FluidPoroPhaseDofPressure(Teuchos::RCP<MAT::PAR::Material> matdata);

      /// initialize
      void Initialize() override;

      /// return type of phase law
      INPAR::MAT::MaterialType PoroPhaseLawType() const override;

      /// mark dofs associated with this phase in a given row (=numphase) in a matrix
      void FillDoFMatrix(CORE::LINALG::SerialDenseMatrix& dofmat, int numphase) const override;

      /// evaluate saturation of the phase
      double EvaluateSaturation(int phasenum, const std::vector<double>& state,
          const std::vector<double>& pressure) const override;

      /// evaluate the generalized(!) pressure of this phase
      double EvaluateGenPressure(int phasenum, const std::vector<double>& state) const override;

      //! evaluate derivative of saturation with respect to pressure
      double EvaluateDerivOfSaturationWrtPressure(
          int phasenum, int doftoderive, const std::vector<double>& pressure) const override;

      //! evaluate 2nd derivative of saturation with respect to pressure
      double EvaluateSecondDerivOfSaturationWrtPressure(int phasenum, int firstdoftoderive,
          int seconddoftoderive, const std::vector<double>& pressure) const override;

      //! evaluate derivative of degree of freedom with respect to pressure
      double EvaluateDerivOfDofWrtPressure(
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
      explicit FluidPoroPhaseDofDiffPressure(Teuchos::RCP<MAT::PAR::Material> matdata);

      /// initialize
      void Initialize() override;

      /// return type of phase law
      INPAR::MAT::MaterialType PoroPhaseLawType() const override;

      /// mark dofs associated with this phase in a given row (=numphase) in a matrix
      void FillDoFMatrix(CORE::LINALG::SerialDenseMatrix& dofmat, int numphase) const override;

      /// evaluate saturation of the phase
      double EvaluateSaturation(int phasenum, const std::vector<double>& state,
          const std::vector<double>& pressure) const override;

      /// evaluate the generalized(!) pressure of this phase
      double EvaluateGenPressure(int phasenum, const std::vector<double>& state) const override;

      //! evaluate derivative of saturation with respect to pressure
      double EvaluateDerivOfSaturationWrtPressure(
          int phasenum, int doftoderive, const std::vector<double>& pressure) const override;

      //! evaluate 2nd derivative of saturation with respect to pressure
      double EvaluateSecondDerivOfSaturationWrtPressure(int phasenum, int firstdoftoderive,
          int seconddoftoderive, const std::vector<double>& pressure) const override;

      //! evaluate derivative of degree of freedom with respect to pressure
      double EvaluateDerivOfDofWrtPressure(
          int phasenum, int doftoderive, const std::vector<double>& state) const override;

     protected:
      ///  pressure-coefficients defining the differential pressure
      const std::vector<int>* diffpresCoeffs_;
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
      explicit FluidPoroPhaseDofSaturation(Teuchos::RCP<MAT::PAR::Material> matdata);

      /// initialize
      void Initialize() override;

      /// return type of phase law
      INPAR::MAT::MaterialType PoroPhaseLawType() const override;

      /// mark dofs associated with this phase in a given row (=numphase) in a matrix
      void FillDoFMatrix(CORE::LINALG::SerialDenseMatrix& dofmat, int numphase) const override;

      /// evaluate saturation of the phase
      double EvaluateSaturation(int phasenum, const std::vector<double>& state,
          const std::vector<double>& pressure) const override;

      /// evaluate the generalized(!) pressure of this phase
      double EvaluateGenPressure(int phasenum, const std::vector<double>& state) const override;

      //! evaluate derivative of saturation with respect to pressure
      double EvaluateDerivOfSaturationWrtPressure(
          int phasenum, int doftoderive, const std::vector<double>& pressure) const override;

      //! evaluate 2nd derivative of saturation with respect to pressure
      double EvaluateSecondDerivOfSaturationWrtPressure(int phasenum, int firstdoftoderive,
          int seconddoftoderive, const std::vector<double>& pressure) const override;

      //! evaluate derivative of degree of freedom with respect to pressure
      double EvaluateDerivOfDofWrtPressure(
          int phasenum, int doftoderive, const std::vector<double>& pressure) const override;

     protected:
      /// ID of pressure-saturation law
      const int phaselawId_;
      /// implementation of pressure-saturation law
      FluidPoroPhaseLaw* phaselaw_;
    };
  }  // namespace PAR
}  // namespace MAT

BACI_NAMESPACE_CLOSE

#endif  // MAT_FLUIDPORO_SINGLEPHASEDOF_H
