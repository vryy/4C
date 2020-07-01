/*----------------------------------------------------------------------*/
/*! \file
 \brief material defining the degree of freedom of phase within a multiphase porous fluid.

   \level 3

 *----------------------------------------------------------------------*/

#ifndef SRC_DRT_MAT_FLUIDPORO_SINGLEPHASEDOF_H_
#define SRC_DRT_MAT_FLUIDPORO_SINGLEPHASEDOF_H_


#include "matpar_parameter.H"
#include "../drt_inpar/inpar_material.H"

#include <Epetra_SerialDenseMatrix.h>

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

      /// destructor
      virtual ~FluidPoroPhaseDof() { ; }

      //! build the phase dof
      static FluidPoroPhaseDof* CreatePhaseDof(int phasedofId);

      /// create material instance of matching type with my parameters
      virtual Teuchos::RCP<MAT::Material> CreateMaterial() { return Teuchos::null; };

      /// initialize
      virtual void Initialize() = 0;

      /// return type of phase law
      virtual INPAR::MAT::MaterialType PoroPhaseLawType() const = 0;

      /// mark dofs associated with this phase in a given row (=numphase) in a matrix
      virtual void FillDoFMatrix(Epetra_SerialDenseMatrix& dofmat, int numphase) const = 0;

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

      /// destructor
      virtual ~FluidPoroPhaseDofPressure() { ; }

      /// initialize
      virtual void Initialize();

      /// return type of phase law
      INPAR::MAT::MaterialType PoroPhaseLawType() const;

      /// mark dofs associated with this phase in a given row (=numphase) in a matrix
      virtual void FillDoFMatrix(Epetra_SerialDenseMatrix& dofmat, int numphase) const;

      /// evaluate saturation of the phase
      virtual double EvaluateSaturation(int phasenum, const std::vector<double>& state,
          const std::vector<double>& pressure) const;

      /// evaluate the generalized(!) pressure of this phase
      virtual double EvaluateGenPressure(int phasenum, const std::vector<double>& state) const;

      //! evaluate derivative of saturation with respect to pressure
      virtual double EvaluateDerivOfSaturationWrtPressure(
          int phasenum, int doftoderive, const std::vector<double>& pressure) const;

      //! evaluate 2nd derivative of saturation with respect to pressure
      virtual double EvaluateSecondDerivOfSaturationWrtPressure(int phasenum, int firstdoftoderive,
          int seconddoftoderive, const std::vector<double>& pressure) const;

      //! evaluate derivative of degree of freedom with respect to pressure
      virtual double EvaluateDerivOfDofWrtPressure(
          int phasenum, int doftoderive, const std::vector<double>& state) const;

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

      /// destructor
      virtual ~FluidPoroPhaseDofDiffPressure() { ; }

      /// initialize
      virtual void Initialize();

      /// return type of phase law
      INPAR::MAT::MaterialType PoroPhaseLawType() const;

      /// mark dofs associated with this phase in a given row (=numphase) in a matrix
      virtual void FillDoFMatrix(Epetra_SerialDenseMatrix& dofmat, int numphase) const;

      /// evaluate saturation of the phase
      virtual double EvaluateSaturation(int phasenum, const std::vector<double>& state,
          const std::vector<double>& pressure) const;

      /// evaluate the generalized(!) pressure of this phase
      virtual double EvaluateGenPressure(int phasenum, const std::vector<double>& state) const;

      //! evaluate derivative of saturation with respect to pressure
      virtual double EvaluateDerivOfSaturationWrtPressure(
          int phasenum, int doftoderive, const std::vector<double>& pressure) const;

      //! evaluate 2nd derivative of saturation with respect to pressure
      virtual double EvaluateSecondDerivOfSaturationWrtPressure(int phasenum, int firstdoftoderive,
          int seconddoftoderive, const std::vector<double>& pressure) const;

      //! evaluate derivative of degree of freedom with respect to pressure
      virtual double EvaluateDerivOfDofWrtPressure(
          int phasenum, int doftoderive, const std::vector<double>& state) const;

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

      /// destructor
      virtual ~FluidPoroPhaseDofSaturation() { ; }

      /// initialize
      virtual void Initialize();

      /// return type of phase law
      INPAR::MAT::MaterialType PoroPhaseLawType() const;

      /// mark dofs associated with this phase in a given row (=numphase) in a matrix
      virtual void FillDoFMatrix(Epetra_SerialDenseMatrix& dofmat, int numphase) const;

      /// evaluate saturation of the phase
      virtual double EvaluateSaturation(int phasenum, const std::vector<double>& state,
          const std::vector<double>& pressure) const;

      /// evaluate the generalized(!) pressure of this phase
      virtual double EvaluateGenPressure(int phasenum, const std::vector<double>& state) const;

      //! evaluate derivative of saturation with respect to pressure
      virtual double EvaluateDerivOfSaturationWrtPressure(
          int phasenum, int doftoderive, const std::vector<double>& pressure) const;

      //! evaluate 2nd derivative of saturation with respect to pressure
      virtual double EvaluateSecondDerivOfSaturationWrtPressure(int phasenum, int firstdoftoderive,
          int seconddoftoderive, const std::vector<double>& pressure) const;

      //! evaluate derivative of degree of freedom with respect to pressure
      virtual double EvaluateDerivOfDofWrtPressure(
          int phasenum, int doftoderive, const std::vector<double>& state) const;

     protected:
      /// ID of pressure-saturation law
      const int phaselawId_;
      /// implementation of pressure-saturation law
      FluidPoroPhaseLaw* phaselaw_;
    };
  }  // namespace PAR
}  // namespace MAT

#endif /* SRC_DRT_MAT_FLUIDPORO_SINGLEPHASEDOF_H_ */
