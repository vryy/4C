/*----------------------------------------------------------------------*/
/*!
 \file scatra_reaction_coupling.H

 \brief helper class encapsulating the reaction terms and its derivatives

   \level 3

   \maintainer  Anh-Tu Vuong
                vuong@lnm.mw.tum.de
                http://www.lnm.mw.tum.de
                089 - 289-15251
 *----------------------------------------------------------------------*/

#ifndef SRC_DRT_MAT_SCATRA_REACTION_COUPLING_H_
#define SRC_DRT_MAT_SCATRA_REACTION_COUPLING_H_

#include "matpar_parameter.H"
#include "material.H"
#include "../drt_lib/drt_parobjectfactory.H"

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
namespace REACTIONCOUPLING
{
  //! interface class for generic reaction coupling
  class ReactionBase
  {

  public:

    /// standard constructor
    ReactionBase():isinit_(false){};

    /// destructor
    virtual ~ReactionBase() { ; }

    ///helper for calculating advanced reaction terms
    virtual double CalcReaBodyForceTerm(
        const int k,                         //!< current scalar id
        int numscal,                         //!< number of scalars
        const std::vector<double>& phinp,    //!< scalar values at t_(n+1)
        const std::vector<double>& couprole, //!< coupling role vector
        const double scale_reac              //!< scaling factor for reaction term (= reaction coefficient * stoichometry)
        ) =0;

    ///helper for calculating advanced reaction term derivatives
    virtual void CalcReaBodyForceDeriv(
        const int k,                         //!< current scalar id
        int numscal,                         //!< number of scalars
        std::vector<double>& derivs,         //!< vector with derivatives (to be filled)
        const std::vector<double>& phinp,    //!< scalar values at t_(n+1)
        const std::vector<double>& couprole, //!< coupling role vector
        const std::vector<double>& reacstart,//!< reaction start vector
        const double scale_reac              //!< scaling factor for reaction term (= reaction coefficient * stoichometry)
        ) =0;

    /// initialization (to be called by derived classes)
    virtual void Initialize(
        int numscal,                         //!< number of scalars
        const std::vector<double>& couprole  //!< coupling role vector
        ) {isinit_=true;};

    /// check for initialization
    bool IsInit() const {return isinit_;};

  private:
    // initialization flag
    bool isinit_;
  };

  /*----------------------------------------------------------------------*
   * **********************************************************************
   *----------------------------------------------------------------------*/

  //! simple multiplicative reaction coupling
  class SimpleMultiplicative : public ReactionBase
  {

  public:

    /// standard constructor
    SimpleMultiplicative(){};

    /// destructor
    virtual ~SimpleMultiplicative() { ; }

    /// initialization
    void Initialize(
        int numscal,                         //!< number of scalars
        const std::vector<double>& couprole  //!< coupling role vector
        ) ;

    ///helper for calculating advanced reaction terms
    double CalcReaBodyForceTerm(
        int k,                               //!< current scalar id
        int numscal,                         //!< number of scalars
        const std::vector<double>& phinp,    //!< scalar values at t_(n+1)
        const std::vector<double>& couprole, //!< coupling role vector
        double scale_reac                    //!< scaling factor for reaction term (= reaction coefficient * stoichometry)
        ) ;

    ///helper for calculating advanced reaction term derivatives
    void CalcReaBodyForceDeriv(
        int k,                               //!< current scalar id
        int numscal,                         //!< number of scalars
        std::vector<double>& derivs,         //!< vector with derivatives (to be filled)
        const std::vector<double>& phinp,    //!< scalar values at t_(n+1)
        const std::vector<double>& couprole, //!< coupling role vector
        const std::vector<double>& reacstart,//!< reaction start vector
        double scale_reac                    //!< scaling factor for reaction term (= reaction coefficient * stoichometry)
        ) ;

  };

  /*----------------------------------------------------------------------*
   * **********************************************************************
   *----------------------------------------------------------------------*/

  //! power multiplicative reaction coupling
  class PowerMultiplicative : public ReactionBase
  {

  public:

    /// standard constructor
    PowerMultiplicative(){};

    /// destructor
    virtual ~PowerMultiplicative() { ; }

    /// initialization
    void Initialize(
        int numscal,                         //!< number of scalars
        const std::vector<double>& couprole  //!< coupling role vector
        ) ;

    ///helper for calculating advanced reaction terms
    double CalcReaBodyForceTerm(
        const int k,                         //!< current scalar id
        int numscal,                         //!< number of scalars
        const std::vector<double>& phinp,    //!< scalar values at t_(n+1)
        const std::vector<double>& couprole, //!< coupling role vector
        const double scale_reac              //!< scaling factor for reaction term (= reaction coefficient * stoichometry)
        ) ;

    ///helper for calculating advanced reaction term derivatives
    void CalcReaBodyForceDeriv(
        const int k,                         //!< current scalar id
        int numscal,                         //!< number of scalars
        std::vector<double>& derivs,         //!< vector with derivatives (to be filled)
        const std::vector<double>& phinp,    //!< scalar values at t_(n+1)
        const std::vector<double>& couprole, //!< coupling role vector
        const std::vector<double>& reacstart,//!< reaction start vector
        const double scale_reac              //!< scaling factor for reaction term (= reaction coefficient * stoichometry)
        ) ;

  };

  /*----------------------------------------------------------------------*
   * **********************************************************************
   *----------------------------------------------------------------------*/

  //! simple multiplicative reaction coupling
  class Constant : public ReactionBase
  {

  public:

    /// standard constructor
    Constant(){};

    /// destructor
    virtual ~Constant() { ; }

    /// initialization
    void Initialize(
        int numscal,                         //!< number of scalars
        const std::vector<double>& couprole  //!< coupling role vector
        ) ;

    ///helper for calculating advanced reaction terms
    double CalcReaBodyForceTerm(
        const int k,                         //!< current scalar id
        int numscal,                         //!< number of scalars
        const std::vector<double>& phinp,    //!< scalar values at t_(n+1)
        const std::vector<double>& couprole, //!< coupling role vector
        const double scale_reac              //!< scaling factor for reaction term (= reaction coefficient * stoichometry)
        ) ;

    ///helper for calculating advanced reaction term derivatives
    void CalcReaBodyForceDeriv(
        const int k,                         //!< current scalar id
        int numscal,                         //!< number of scalars
        std::vector<double>& derivs,         //!< vector with derivatives (to be filled)
        const std::vector<double>& phinp,    //!< scalar values at t_(n+1)
        const std::vector<double>& couprole, //!< coupling role vector
        const std::vector<double>& reacstart,//!< reaction start vector
        const double scale_reac              //!< scaling factor for reaction term (= reaction coefficient * stoichometry)
        ) ;

  };

  /*----------------------------------------------------------------------*
   * **********************************************************************
   *----------------------------------------------------------------------*/

  //! simple multiplicative reaction coupling
  class MichaelisMenten : public ReactionBase
  {

  public:

    /// standard constructor
    MichaelisMenten(){};

    /// destructor
    virtual ~MichaelisMenten() { ; }

    /// initialization
    void Initialize(
        int numscal,                         //!< number of scalars
        const std::vector<double>& couprole  //!< coupling role vector
        ) ;

    ///helper for calculating advanced reaction terms
    double CalcReaBodyForceTerm(
        const int k,                         //!< current scalar id
        int numscal,                         //!< number of scalars
        const std::vector<double>& phinp,    //!< scalar values at t_(n+1)
        const std::vector<double>& couprole, //!< coupling role vector
        const double scale_reac              //!< scaling factor for reaction term (= reaction coefficient * stoichometry)
        ) ;

    ///helper for calculating advanced reaction term derivatives
    void CalcReaBodyForceDeriv(
        const int k,                         //!< current scalar id
        int numscal,                         //!< number of scalars
        std::vector<double>& derivs,         //!< vector with derivatives (to be filled)
        const std::vector<double>& phinp,    //!< scalar values at t_(n+1)
        const std::vector<double>& couprole, //!< coupling role vector
        const std::vector<double>& reacstart,//!< reaction start vector
        const double scale_reac              //!< scaling factor for reaction term (= reaction coefficient * stoichometry)
        ) ;

  };

  /*----------------------------------------------------------------------*
   * **********************************************************************
   *----------------------------------------------------------------------*/

  //! simple multiplicative reaction coupling
  class ByFunction : public ReactionBase
  {

  public:

    /// standard constructor
    ByFunction(){};

    /// destructor
    virtual ~ByFunction() { ; }

    /// initialization
    void Initialize(
        int numscal,                         //!< number of scalars
        const std::vector<double>& couprole  //!< coupling role vector
        ) ;

    ///helper for calculating advanced reaction terms
    double CalcReaBodyForceTerm(
        const int k,                         //!< current scalar id
        int numscal,                         //!< number of scalars
        const std::vector<double>& phinp,    //!< scalar values at t_(n+1)
        const std::vector<double>& couprole, //!< coupling role vector
        const double scale_reac              //!< scaling factor for reaction term (= reaction coefficient * stoichometry)
        ) ;

    ///helper for calculating advanced reaction term derivatives
    void CalcReaBodyForceDeriv(
        const int k,                         //!< current scalar id
        int numscal,                         //!< number of scalars
        std::vector<double>& derivs,         //!< vector with derivatives (to be filled)
        const std::vector<double>& phinp,    //!< scalar values at t_(n+1)
        const std::vector<double>& couprole, //!< coupling role vector
        const std::vector<double>& reacstart,//!< reaction start vector
        const double scale_reac              //!< scaling factor for reaction term (= reaction coefficient * stoichometry)
        ) ;

  private:
    /// helper for evaluation by function
    void BuildPhiVectorForFunction(
        const std::vector<double>& phinp_org, //!< scalar values at t_(n+1)
        int numscal                           //!< number of scalars
        );

    /// helper method for accessing variable expression function
    inline DRT::UTILS::VariableExprFunction& Function(int functnum) const;

    /// variable vector for function evaluation
    std::vector<std::pair<std::string,double> > variables_;

  };


} // namespace REACTIONCOUPLING

} // namespace PAR
} // namespace MAT




#endif /* SRC_DRT_MAT_SCATRA_REACTION_COUPLING_H_ */
