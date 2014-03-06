/*----------------------------------------------------------------------*/
/*!
 \file growth_law.H

 \brief

 <pre>
   Maintainer: Anh-Tu Vuong
               vuong@lnm.mw.tum.de
               http://www.lnm.mw.tum.de
               089 - 289-15251
 </pre>
 *----------------------------------------------------------------------*/


#ifndef GROWTH_LAW_H_
#define GROWTH_LAW_H_

#include <Teuchos_RCP.hpp>
#include "matpar_parameter.H"
#include "material.H"
#include "../linalg/linalg_fixedsizematrix.H"


namespace MAT
{
  class GrowthLaw
  {

  public:

    /// construct empty material object
    GrowthLaw(){};

    /// destructor
    virtual ~GrowthLaw() { ; }

    /// Evaluate growth function
    virtual void EvaluateGrowthFunction
    (
        double & growthfunc,
        double traceM,
        double theta
    ) =0;

    /// Evaluate derivative of growth function w.r.t. growth factor
    virtual void EvaluateGrowthFunctionDerivTheta
    (
        double & dgrowthfunctheta,
        double traceM,
        double theta,
        const LINALG::Matrix<NUM_STRESS_3D, 1>& Cdach,
        const LINALG::Matrix<NUM_STRESS_3D, NUM_STRESS_3D>& cmatelastic
    ) =0;

    /// Evaluate derivative of growth function w.r.t. Green-Lagrange-Strain
    virtual void EvaluateGrowthFunctionDerivC
    (
        LINALG::Matrix<NUM_STRESS_3D, 1>& dgrowthfuncdC,
        double traceM,
        double theta,
        const LINALG::Matrix<NUM_STRESS_3D, 1>& C,
        const LINALG::Matrix<NUM_STRESS_3D, 1>& S,
        const LINALG::Matrix<NUM_STRESS_3D, NUM_STRESS_3D>& cmat
    ) =0;

  };

  namespace PAR
  {
    /*----------------------------------------------------------------------*/
    /// material parameters
    class GrowthLawLinear : public Parameter
    {
    public:

      /// standard constructor
      GrowthLawLinear(
        Teuchos::RCP<MAT::PAR::Material> matdata
        );

      /// destructor
      virtual ~GrowthLawLinear() { ; }

      /// @name material parameters
      //@{
      /// growth law parameter kthetaplus
      const double kthetaplus_;
      /// growth law parameter mthetaplus
      const double mthetaplus_;
      /// growth law parameter kthetaminus
      const double kthetaminus_;
      /// growth law parameter mthetaminus
      const double mthetaminus_;
      /// homeostatic value for mandelstress
      const double hommandel_;
      //@}

      /// create material instance of matching type with my parameters
      virtual Teuchos::RCP<MAT::Material> CreateMaterial();

      /// create growth law instance of matching type with my parameters
      virtual Teuchos::RCP<MAT::GrowthLaw> CreateGrowthLaw();

    };  // class Growth

  }  // namespace PAR

  class GrowthLawLinear: public GrowthLaw
  {

  public:

    /// construct empty material object
    GrowthLawLinear();
    explicit GrowthLawLinear(MAT::PAR::GrowthLawLinear* params);

    /// destructor
    virtual ~GrowthLawLinear() { ; }

    /// Evaluate growth function
    virtual void EvaluateGrowthFunction
    (
        double & growthfunc,
        double traceM,
        double theta
    ) ;

    /// Evaluate derivative of growth function w.r.t. growth factor
    virtual void EvaluateGrowthFunctionDerivTheta
    (
        double & dgrowthfunctheta,
        double traceM,
        double theta,
        const LINALG::Matrix<NUM_STRESS_3D, 1>& Cdach,
        const LINALG::Matrix<NUM_STRESS_3D, NUM_STRESS_3D>& cmatelastic
    ) ;

    /// Evaluate derivative of growth function w.r.t. Green-Lagrange-Strain
    virtual void EvaluateGrowthFunctionDerivC
    (
        LINALG::Matrix<NUM_STRESS_3D, 1>& dgrowthfuncdC,
        double traceM,
        double theta,
        const LINALG::Matrix<NUM_STRESS_3D, 1>& C,
        const LINALG::Matrix<NUM_STRESS_3D, 1>& S,
        const LINALG::Matrix<NUM_STRESS_3D, NUM_STRESS_3D>& cmat
    ) ;

  private:
    MAT::PAR::GrowthLawLinear* params_;

  };

  namespace PAR
  {
    /*----------------------------------------------------------------------*/
    /// material parameters
    class GrowthLawExp : public Parameter
    {
    public:

      /// standard constructor
      GrowthLawExp(
        Teuchos::RCP<MAT::PAR::Material> matdata
        );

      /// destructor
      virtual ~GrowthLawExp() { ; }

      /// @name material parameters
      //@{
      /// reference value for trace of mandel stress
      const double mandel_;
      //@}

      /// create material instance of matching type with my parameters
      virtual Teuchos::RCP<MAT::Material> CreateMaterial();

      /// create growth law instance of matching type with my parameters
      virtual Teuchos::RCP<MAT::GrowthLaw> CreateGrowthLaw();

    };  // class Growth

  }  // namespace PAR

  class GrowthLawExp: public GrowthLaw
  {

  public:

    /// construct empty material object
    GrowthLawExp();
    explicit GrowthLawExp(MAT::PAR::GrowthLawExp* params);

    /// destructor
    virtual ~GrowthLawExp() { ; }

    /// Evaluate growth function
    virtual void EvaluateGrowthFunction
    (
        double & growthfunc,
        double traceM,
        double theta
    ) ;

    /// Evaluate derivative of growth function w.r.t. growth factor
    virtual void EvaluateGrowthFunctionDerivTheta
    (
        double & dgrowthfunctheta,
        double traceM,
        double theta,
        const LINALG::Matrix<NUM_STRESS_3D, 1>& Cdach,
        const LINALG::Matrix<NUM_STRESS_3D, NUM_STRESS_3D>& cmatelastic
    ) ;

    /// Evaluate derivative of growth function w.r.t. Green-Lagrange-Strain
    virtual void EvaluateGrowthFunctionDerivC
    (
        LINALG::Matrix<NUM_STRESS_3D, 1>& dgrowthfuncdC,
        double traceM,
        double theta,
        const LINALG::Matrix<NUM_STRESS_3D, 1>& C,
        const LINALG::Matrix<NUM_STRESS_3D, 1>& S,
        const LINALG::Matrix<NUM_STRESS_3D, NUM_STRESS_3D>& cmat
    ) ;

  private:
    MAT::PAR::GrowthLawExp* params_;

  };

} // namespace MAT


#endif /* GROWTH_LAW_H_ */
