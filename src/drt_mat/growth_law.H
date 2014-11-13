/*----------------------------------------------------------------------*/
/*!
 \file growth_law.H

 \brief

This file contains routines needed for the calculation of the volumetric growth
parameter theta. The returned results can be either theta itself or its derivative

 <pre>
   Maintainer: Moritz Thon
               thon@lnm.mw.tum.de
               http://www.mhpc.mw.tum.de
               089 - 289-10364
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

//    /// construct empty material object
//    GrowthLaw(){};

    /// construct empty material object
    GrowthLaw();
    explicit GrowthLaw(MAT::PAR::Parameter* params);

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

    /// Calculate Theta from the given mean concentrations
    virtual double CalculateTheta
    (
        Teuchos::RCP<std::vector<double> > concentrations
    ) =0;

    /// material type
    virtual INPAR::MAT::MaterialType MaterialType() const = 0;

    MAT::PAR::Parameter* Parameter() {return params_;};

  private:

    MAT::PAR::Parameter* params_;
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
      /// maximal growth stretch
      const double thetaplus_;
      /// growth law parameter kthetaplus
      const double kthetaplus_;
      /// growth law parameter mthetaplus
      const double mthetaplus_;
      /// minimal growth stretch
      const double thetaminus_;
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

    /// Calculate Theta from the given mean concentrations
    virtual double CalculateTheta
    (
        Teuchos::RCP<std::vector<double> > concentrations
    ) {dserror("In case of GrowthLawLinear one doesn't have to compute this!"); return -1.0; };

    /// material type
    virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::m_growth_linear; };

    MAT::PAR::GrowthLawLinear* Parameter() {return dynamic_cast<MAT::PAR::GrowthLawLinear*>(MAT::GrowthLaw::Parameter()); };
  private:
    //MAT::PAR::GrowthLawLinear* params_;

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

    /// Calculate Theta from the given mean concentrations
    virtual double CalculateTheta
    (
        Teuchos::RCP<std::vector<double> > concentrations
    ) {dserror("In case of GrowthLawExp one doesn't have to compute this!"); return -1.0; };

    /// material type
    virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::m_growth_exponential; };

    MAT::PAR::GrowthLawExp* Parameter() {return dynamic_cast<MAT::PAR::GrowthLawExp*>(MAT::GrowthLaw::Parameter()); };
  private:
    //MAT::PAR::GrowthLawExp* params_;

  };


  namespace PAR
  {
    /*----------------------------------------------------------------------*/
    /// material parameters
    class GrowthLawAC : public Parameter
    {
    public:

      /// standard constructor
      GrowthLawAC(
        Teuchos::RCP<MAT::PAR::Material> matdata
        );

      /// destructor
      virtual ~GrowthLawAC() { ; }

      /// @name material parameters
      //@{
      /// maximal growth stretch
      const int Sc1_;
      /// growth law parameter kthetaplus
      const double alpha_;
      /// growth law parameter mthetaplus
      const int Sc2_;
      /// minimal growth stretch
      const double beta_;
      //@}

      /// create material instance of matching type with my parameters
      virtual Teuchos::RCP<MAT::Material> CreateMaterial();

      /// create growth law instance of matching type with my parameters
      virtual Teuchos::RCP<MAT::GrowthLaw> CreateGrowthLaw();

    };  // class Growth

  }  // namespace PAR

class GrowthLawAC: public GrowthLaw
  {

  public:

    /// construct empty material object
    GrowthLawAC();
    explicit GrowthLawAC(MAT::PAR::GrowthLawAC* params);

    /// destructor
    virtual ~GrowthLawAC() { ; }

    /// Evaluate growth function
    virtual void EvaluateGrowthFunction
    (
        double & growthfunc,
        double traceM,
        double theta
    ) {dserror("In case of GrowthLawAC one doesn't have to compute this!");};

    /// Evaluate derivative of growth function w.r.t. growth factor
    virtual void EvaluateGrowthFunctionDerivTheta
    (
        double & dgrowthfunctheta,
        double traceM,
        double theta,
        const LINALG::Matrix<NUM_STRESS_3D, 1>& Cdach,
        const LINALG::Matrix<NUM_STRESS_3D, NUM_STRESS_3D>& cmatelastic
    ) {dserror("In case of GrowthLawAC one doesn't have to compute this!");};

    /// Evaluate derivative of growth function w.r.t. Green-Lagrange-Strain
    virtual void EvaluateGrowthFunctionDerivC
    (
        LINALG::Matrix<NUM_STRESS_3D, 1>& dgrowthfuncdC,
        double traceM,
        double theta,
        const LINALG::Matrix<NUM_STRESS_3D, 1>& C,
        const LINALG::Matrix<NUM_STRESS_3D, 1>& S,
        const LINALG::Matrix<NUM_STRESS_3D, NUM_STRESS_3D>& cmat
    ) {dserror("In case of GrowthLawAC one doesn't have to compute this!");};

    /// Calculate Theta from the given mean concentrations
    virtual double CalculateTheta
    (
        Teuchos::RCP<std::vector<double> > concentrations //pointer to a vector containing element wise mean concentrations
    );

    /// material type
    virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::m_growth_ac; };

    MAT::PAR::GrowthLawAC* Parameter() {return dynamic_cast<MAT::PAR::GrowthLawAC*>(MAT::GrowthLaw::Parameter()); };

  private:
    //MAT::PAR::GrowthLawAC* params_;

  };

} // namespace MAT


#endif /* GROWTH_LAW_H_ */
