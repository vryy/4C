/*----------------------------------------------------------------------*/
/*!
 \file fluidporo_singlephaselaw.H

 \brief a material defining the pressure-saturation relationship of
        fluid phase within a multiphase porous fluid

 \level 3

 \maintainer   Anh-Tu Vuong
               vuong@lnm.mw.tum.de
               http://www.lnm.mw.tum.de
               089 - 289-15251
 *----------------------------------------------------------------------*/

#ifndef SRC_DRT_MAT_FLUIDPORO_SINGLEPHASELAW_H_
#define SRC_DRT_MAT_FLUIDPORO_SINGLEPHASELAW_H_

#include "matpar_parameter.H"

namespace DRT
{
namespace UTILS
{
  class VariableExprFunction;
}
}

namespace MAT
{
namespace PAR
{
  //! interface class for generic phase (pressure-saturation) law
  class FluidPoroPhaseLaw: public Parameter
  {

  public:

    /// standard constructor
    explicit FluidPoroPhaseLaw(Teuchos::RCP<MAT::PAR::Material> matdata);

    /// destructor
    virtual ~FluidPoroPhaseLaw()
    {
      ;
    }

    /// create material instance of matching type with my parameters
    virtual Teuchos::RCP<MAT::Material> CreateMaterial() {return Teuchos::null;};

    /// initialize
    virtual void Initialize()=0;

    /// @name material parameters
    //@{
    //! number of degrees of freedom
    const int numdof_;
    //! coefficients defining the pressures involved int the pressure-saturation law
    const std::vector<int>* presids_;
    //@}

    //! evaluate saturation of phase
    virtual double EvaluateSaturation(
        const std::vector<double>& pressure) const = 0;

    //! evaluate derivative of saturation w.r.t. pressure
    virtual double EvaluateDerivOfSaturationWrtPressure(
        int doftoderive,
        const std::vector<double>& state) const = 0;

    //! evaluate generalized pressure related to phase law
    virtual double EvaluateGenPressure(double saturation) const = 0;

    //! evaluate derivative of pressure w.r.t. saturation
    virtual double EvaluateDerivOfPressureWrtSaturation(int doftoderive,
        double saturation) const = 0;
  };

  /*----------------------------------------------------------------------*/
  //! a linear pressure-saturation relation. Only for testing, no physical meaning.
  class FluidPoroPhaseLawLinear: public FluidPoroPhaseLaw
  {
  public:

    /// standard constructor
    explicit FluidPoroPhaseLawLinear(Teuchos::RCP<MAT::PAR::Material> matdata);

    /// destructor
    virtual ~FluidPoroPhaseLawLinear()
    {
      ;
    }

    /// initialize
    virtual void Initialize(){return;};

    /// @name material parameters
    //@{
    //! 'relative tension' parameter (proportionality coefficient)
    const double reltensions_;
    //! saturation value at zero pressure
    const double  sat0_;
    //@}

    //! evaluate saturation of phase
    virtual double EvaluateSaturation(
        const std::vector<double>& pressure) const;

    //! evaluate generalized pressure related to phase law
    virtual double EvaluateGenPressure(double saturation) const;

    //! evaluate derivative of saturation w.r.t. pressure
    virtual double EvaluateDerivOfSaturationWrtPressure(
        int doftoderive,
        const std::vector<double>& state) const;

    //! evaluate derivative of pressure w.r.t. saturtion
    virtual double EvaluateDerivOfPressureWrtSaturation(int doftoderive,
        double saturation) const;

  };
  // class FluidPoroPhaseLawLinear

  /*----------------------------------------------------------------------*/
  //! tangent pressure-saturation relationship
  class FluidPoroPhaseLawTangent: public FluidPoroPhaseLaw
  {
  public:

    /// standard constructor
    explicit FluidPoroPhaseLawTangent(Teuchos::RCP<MAT::PAR::Material> matdata);

    /// destructor
    virtual ~FluidPoroPhaseLawTangent()
    {
      ;
    }

    /// initialize
    virtual void Initialize(){return;};

    /// @name material parameters
    //@{
    //! relative tension coefficient
    const double reltensions_;
    //! exponent within phase law
    const double exp_;
    //! saturation value at zero pressure
    const double  sat0_;
    //@}

    //! evaluate saturation of phase
    virtual double EvaluateSaturation(const std::vector<double>& pressure) const;

    //! evaluate derivative of saturation w.r.t. pressure
    virtual double EvaluateDerivOfSaturationWrtPressure(
        int doftoderive,
        const std::vector<double>& state) const;

    //! evaluate generalized pressure related to phase law
    virtual double EvaluateGenPressure(double saturation) const;

    //! evaluate derivative of pressure w.r.t. saturation
    virtual double EvaluateDerivOfPressureWrtSaturation(
        int doftoderive,
        double saturation) const;

  };
  // class FluidPoroPhaseLawTangent

  /*----------------------------------------------------------------------*/
  //! pressure-saturation relationship defined by function
  class FluidPoroPhaseLawByFunction: public FluidPoroPhaseLaw
  {
  public:

    /// standard constructor
    explicit FluidPoroPhaseLawByFunction(Teuchos::RCP<MAT::PAR::Material> matdata);

    /// destructor
    virtual ~FluidPoroPhaseLawByFunction()
    {
      ;
    }

    /// initialize
    virtual void Initialize();

    /// @name material parameters
    //@{
    //! relative tension coefficient
    const int functionID_saturation_;
    //! exponent within phase law
    const int functionID_pressure_;
    //@}

    //! evaluate saturation of phase
    virtual double EvaluateSaturation(const std::vector<double>& pressure) const;

    //! evaluate derivative of saturation w.r.t. pressure
    virtual double EvaluateDerivOfSaturationWrtPressure(
        int doftoderive,
        const std::vector<double>& state) const;

    //! evaluate generalized pressure related to phase law
    virtual double EvaluateGenPressure(double saturation) const;

    //! evaluate derivative of pressure w.r.t. saturation
    virtual double EvaluateDerivOfPressureWrtSaturation(
        int doftoderive,
        double saturation) const;

  private:
    inline DRT::UTILS::VariableExprFunction& Function(int functnum) const;
  };
  // class FluidPoroPhaseLawTangent

}
}

#endif /* SRC_DRT_MAT_FLUIDPORO_SINGLEPHASELAW_H_ */
