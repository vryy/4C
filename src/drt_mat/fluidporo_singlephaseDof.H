/*----------------------------------------------------------------------*/
/*!
 \file fluidporo_singlephaseDof.H

 \brief material defining the degree of freedom of phase within a multiphase porous fluid.

   \level 3

   \maintainer  Anh-Tu Vuong
                vuong@lnm.mw.tum.de
                http://www.lnm.mw.tum.de
                089 - 289-15251
 *----------------------------------------------------------------------*/

#ifndef SRC_DRT_MAT_FLUIDPORO_SINGLEPHASEDOF_H_
#define SRC_DRT_MAT_FLUIDPORO_SINGLEPHASEDOF_H_


#include "matpar_parameter.H"

#include <Epetra_SerialDenseMatrix.h>

namespace MAT
{
namespace PAR
{
  //forward declaration
  class FluidPoroPhaseLaw;

  //! interface class for generic phase degree of freedom
  class FluidPoroPhaseDof: public Parameter
  {

  public:

    /// standard constructor
    explicit FluidPoroPhaseDof(Teuchos::RCP<MAT::PAR::Material> matdata);

    /// destructor
    virtual ~FluidPoroPhaseDof()
    {
      ;
    }

    /// create material instance of matching type with my parameters
    virtual Teuchos::RCP<MAT::Material> CreateMaterial() {return Teuchos::null;};

    /// mark dofs associated with this phase in a given row (=numphase) in a matrix
    virtual void FillDoFMatrix(
        Epetra_SerialDenseMatrix& dofmat,
        int numphase) const =0;

    /// evaluate saturation of the phase
    virtual double EvaluateSaturation(
        int phasenum,
        const std::vector<double>& state,
        const std::vector<double>& pressure) const =0;

    /// evaluate the generalized(!) pressure of this phase
    virtual double EvaluateGenPressure(
            int phasenum,
            const std::vector<double>& state) const =0;

    //! evaluate derivative of saturation with respect to pressure
    virtual double EvaluateDerivOfSaturationWrtPressure(
        int phasenum,
        int doftoderive,
        const std::vector<double>& state) const =0;

    //! evaluate derivative of degree of freedom with respect to pressure
    virtual double EvaluateDerivOfDofWrtPressure(
            int phasenum,
            int doftoderive,
            const std::vector<double>& state) const =0;

  };

  //! pressure degree of freedom of a single poro fluid phase
  class FluidPoroPhaseDofPressure: public FluidPoroPhaseDof
  {

  public:

    /// standard constructor
    explicit FluidPoroPhaseDofPressure(Teuchos::RCP<MAT::PAR::Material> matdata);

    /// destructor
    virtual ~FluidPoroPhaseDofPressure()
    {
      ;
    }

    /// mark dofs associated with this phase in a given row (=numphase) in a matrix
    virtual void FillDoFMatrix(
        Epetra_SerialDenseMatrix& dofmat,
        int numphase) const ;

    /// evaluate saturation of the phase
    virtual double EvaluateSaturation(
        int phasenum,
        const std::vector<double>& state,
        const std::vector<double>& pressure) const;

    /// evaluate the generalized(!) pressure of this phase
    virtual double EvaluateGenPressure(
            int phasenum,
            const std::vector<double>& state) const;

    //! evaluate derivative of saturation with respect to pressure
    virtual double EvaluateDerivOfSaturationWrtPressure(
        int phasenum,
        int doftoderive,
        const std::vector<double>& state) const;

    //! evaluate derivative of degree of freedom with respect to pressure
    virtual double EvaluateDerivOfDofWrtPressure(
            int phasenum,
            int doftoderive,
            const std::vector<double>& state) const;

  protected:
    /// ID of pressure-saturation law
    const int phaselawId_;
    /// implementation of pressure-saturation law
    FluidPoroPhaseLaw* phaselaw_;
  };

  //! sumed pressure degree of freedom (the last phase) within a porous fluid
  class FluidPoroPhaseDofPressureSum: public FluidPoroPhaseDof
  {

  public:

    /// standard constructor
    explicit FluidPoroPhaseDofPressureSum(Teuchos::RCP<MAT::PAR::Material> matdata);

    /// destructor
    virtual ~FluidPoroPhaseDofPressureSum()
    {
      ;
    }

    /// mark dofs associated with this phase in a given row (=numphase) in a matrix
    virtual void FillDoFMatrix(
        Epetra_SerialDenseMatrix& dofmat,
        int numphase) const ;

    /// evaluate saturation of the phase
    virtual double EvaluateSaturation(
        int phasenum,
        const std::vector<double>& state,
        const std::vector<double>& pressure) const;

    /// evaluate the generalized(!) pressure of this phase
    virtual double EvaluateGenPressure(
            int phasenum,
            const std::vector<double>& state) const;

    //! evaluate derivative of saturation with respect to pressure
    virtual double EvaluateDerivOfSaturationWrtPressure(
        int phasenum,
        int doftoderive,
        const std::vector<double>& state) const;

    //! evaluate derivative of degree of freedom with respect to pressure
    virtual double EvaluateDerivOfDofWrtPressure(
            int phasenum,
            int doftoderive,
            const std::vector<double>& state) const;

  };

  //! differential pressure degree of freedom of a single phase of porous multiphase fluid
  class FluidPoroPhaseDofDiffPressure: public FluidPoroPhaseDof
  {

  public:

    /// standard constructor
    explicit FluidPoroPhaseDofDiffPressure(Teuchos::RCP<MAT::PAR::Material> matdata);

    /// destructor
    virtual ~FluidPoroPhaseDofDiffPressure()
    {
      ;
    }

    /// mark dofs associated with this phase in a given row (=numphase) in a matrix
    virtual void FillDoFMatrix(
        Epetra_SerialDenseMatrix& dofmat,
        int numphase) const ;

    /// evaluate saturation of the phase
    virtual double EvaluateSaturation(
        int phasenum,
        const std::vector<double>& state,
        const std::vector<double>& pressure) const;

    /// evaluate the generalized(!) pressure of this phase
    virtual double EvaluateGenPressure(
            int phasenum,
            const std::vector<double>& state) const;

    //! evaluate derivative of saturation with respect to pressure
    virtual double EvaluateDerivOfSaturationWrtPressure(
        int phasenum,
        int doftoderive,
        const std::vector<double>& state) const;

    //! evaluate derivative of degree of freedom with respect to pressure
    virtual double EvaluateDerivOfDofWrtPressure(
            int phasenum,
            int doftoderive,
            const std::vector<double>& state) const;

  protected:
    ///  pressure-coefficients defining the differential pressure
    const std::vector<int>* diffpresCoeffs_;
    /// ID of pressure-saturation law
    const int phaselawId_;
    /// implementation of pressure-saturation law
    FluidPoroPhaseLaw* phaselaw_;
  };

  //! saturation degree of freedom of a single phase within a multiphase porous fluid
  class FluidPoroPhaseDofSaturation: public FluidPoroPhaseDof
  {

  public:

    /// standard constructor
    explicit FluidPoroPhaseDofSaturation(Teuchos::RCP<MAT::PAR::Material> matdata);

    /// destructor
    virtual ~FluidPoroPhaseDofSaturation()
    {
      ;
    }

    /// mark dofs associated with this phase in a given row (=numphase) in a matrix
    virtual void FillDoFMatrix(
        Epetra_SerialDenseMatrix& dofmat,
        int numphase) const ;

    /// evaluate saturation of the phase
    virtual double EvaluateSaturation(
        int phasenum,
        const std::vector<double>& state,
        const std::vector<double>& pressure) const;

    /// evaluate the generalized(!) pressure of this phase
    virtual double EvaluateGenPressure(
            int phasenum,
            const std::vector<double>& state) const;

    //! evaluate derivative of saturation with respect to pressure
    virtual double EvaluateDerivOfSaturationWrtPressure(
        int phasenum,
        int doftoderive,
        const std::vector<double>& state) const;

    //! evaluate derivative of degree of freedom with respect to pressure
    virtual double EvaluateDerivOfDofWrtPressure(
            int phasenum,
            int doftoderive,
            const std::vector<double>& state) const;

  protected:
    /// ID of pressure-saturation law
    const int phaselawId_;
    /// implementation of pressure-saturation law
    FluidPoroPhaseLaw* phaselaw_;
  };
}
}

#endif /* SRC_DRT_MAT_FLUIDPORO_SINGLEPHASEDOF_H_ */
