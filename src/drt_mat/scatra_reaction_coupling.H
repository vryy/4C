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
#include "scatra_reaction_mat.H"
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
  class ReactionInterface
  {

  public:

    /// factory method
    static Teuchos::RCP<ReactionInterface> CreateReaction(
        MAT::PAR::reaction_coupling couplingtype,  //!< coupling type definig reaction
        bool isreacstart,                          //!< flag for reaction start feature
        const std::vector<double>& reacstart       //!< reaction start vector
        );

    /// standard constructor
    ReactionInterface(){};

    /// destructor
    virtual ~ReactionInterface() { ; }

    /// initialization (to be called by derived classes)
    virtual void Initialize(
        int numscal,                         //!< number of scalars
        const std::vector<double>& couprole  //!< coupling role vector
        ) =0;

    /// check for initialization
    virtual bool IsInit() const =0;

    ///helper for calculating advanced reaction terms
    virtual double CalcReaBodyForceTerm(
        int k,                               //!< current scalar id
        int numscal,                         //!< number of scalars
        const std::vector<double>& phinp,    //!< scalar values at t_(n+1)
        const std::vector<double>& couprole, //!< coupling role vector
        double scale_reac,                   //!< scaling factor for reaction term (= reaction coefficient * stoichometry)
        double scale_phi                     //!< scaling factor for scalar values (used for reference concentrations)
        ) =0;

    ///helper for calculating advanced reaction term derivatives
    virtual void CalcReaBodyForceDeriv(
        int k,                               //!< current scalar id
        int numscal,                         //!< number of scalars
        std::vector<double>& derivs,         //!< vector with derivatives (to be filled)
        const std::vector<double>& phinp,    //!< scalar values at t_(n+1)
        const std::vector<double>& couprole, //!< coupling role vector
        double scale_reac,                   //!< scaling factor for reaction term (= reaction coefficient * stoichometry)
        double scale_phi                     //!< scaling factor for scalar values (used for reference concentrations)
        ) =0;

    ///helper for calculating advanced reaction terms
    virtual double CalcReaBodyForceTerm(
        int k,                                                           //!< current scalar id
        int numscal,                                                     //!< number of scalars
        const std::vector<double>& phinp,                                //!< scalar values at t_(n+1)
        const std::vector<std::pair<std::string,double> >& constants,    //!< vector containing values which are independent of the scalars
        const std::vector<double>& couprole,                             //!< coupling role vector
        double scale_reac,                                               //!< scaling factor for reaction term (= reaction coefficient * stoichometry)
        double scale_phi                                                 //!< scaling factor for scalar values (used for reference concentrations)
        ) =0;

    ///helper for calculating advanced reaction term derivatives
    virtual void CalcReaBodyForceDeriv(
        int k,                                                           //!< current scalar id
        int numscal,                                                     //!< number of scalars
        std::vector<double>& derivs,                                     //!< vector with derivatives (to be filled)
        const std::vector<double>& phinp,                                //!< scalar values at t_(n+1)
        const std::vector<std::pair<std::string,double> >& constants,    //!< vector containing values which are independent of the scalars
        const std::vector<double>& couprole,                             //!< coupling role vector
        double scale_reac,                                               //!< scaling factor for reaction term (= reaction coefficient * stoichometry)
        double scale_phi                                                 //!< scaling factor for scalar values (used for reference concentrations)
        ) =0;
  };

  /*----------------------------------------------------------------------*
   * **********************************************************************
   *----------------------------------------------------------------------*/

  //! reaction start feature for reaction coupling
  //! it wraps another reaction class
  class ReacStart : public ReactionInterface
  {

  public:

    /// standard constructor
    ReacStart(Teuchos::RCP<ReactionInterface> reaction,const std::vector<double>& reacstart) :
        reaction_(reaction),
        reacstart_(reacstart)
    {};

    /// destructor
    virtual ~ReacStart() { ; }

    /// initialization
    void Initialize(
        int numscal,                         //!< number of scalars
        const std::vector<double>& couprole  //!< coupling role vector
        ) ;

    /// check for initialization
    bool IsInit() const {return reaction_->IsInit();};

    ///helper for calculating advanced reaction terms
    virtual double CalcReaBodyForceTerm(
        int k,                               //!< current scalar id
        int numscal,                         //!< number of scalars
        const std::vector<double>& phinp,    //!< scalar values at t_(n+1)
        const std::vector<double>& couprole, //!< coupling role vector
        double scale_reac,                   //!< scaling factor for reaction term (= reaction coefficient * stoichometry)
        double scale_phi                     //!< scaling factor for scalar values (used for reference concentrations)
        );

    ///helper for calculating advanced reaction term derivatives
    virtual void CalcReaBodyForceDeriv(
        int k,                               //!< current scalar id
        int numscal,                         //!< number of scalars
        std::vector<double>& derivs,         //!< vector with derivatives (to be filled)
        const std::vector<double>& phinp,    //!< scalar values at t_(n+1)
        const std::vector<double>& couprole, //!< coupling role vector
        double scale_reac,                   //!< scaling factor for reaction term (= reaction coefficient * stoichometry)
        double scale_phi                     //!< scaling factor for scalar values (used for reference concentrations)
        );

    ///helper for calculating advanced reaction terms
    virtual double CalcReaBodyForceTerm(
        int k,                                                           //!< current scalar id
        int numscal,                                                     //!< number of scalars
        const std::vector<double>& phinp,                                //!< scalar values at t_(n+1)
        const std::vector<std::pair<std::string,double> >& constants,    //!< vector containing values which are independent of the scalars
        const std::vector<double>& couprole,                             //!< coupling role vector
        double scale_reac,                                               //!< scaling factor for reaction term (= reaction coefficient * stoichometry)
        double scale_phi                                                 //!< scaling factor for scalar values (used for reference concentrations)
        );

    ///helper for calculating advanced reaction term derivatives
    virtual void CalcReaBodyForceDeriv(
        int k,                                                           //!< current scalar id
        int numscal,                                                     //!< number of scalars
        std::vector<double>& derivs,                                     //!< vector with derivatives (to be filled)
        const std::vector<double>& phinp,                                //!< scalar values at t_(n+1)
        const std::vector<std::pair<std::string,double> >& constants,    //!< vector containing values which are independent of the scalars
        const std::vector<double>& couprole,                             //!< coupling role vector
        double scale_reac,                                               //!< scaling factor for reaction term (= reaction coefficient * stoichometry)
        double scale_phi                                                 //!< scaling factor for scalar values (used for reference concentrations)
        );

  private:

    /// helper for reaction start feature
    virtual std::vector<double> ModifyPhi(
        const std::vector<double>& phinp
        );

    /// actual reaction
    Teuchos::RCP<ReactionInterface> reaction_;
    /// reacstart vector
    const std::vector<double>& reacstart_;
  };

  /*----------------------------------------------------------------------*
   * **********************************************************************
   *----------------------------------------------------------------------*/

  //! wrapper class for reaction coupling with potential phi scaling
  //! it wraps another reaction class
  class ReactionWithPhiScaling : public ReactionInterface
  {

  public:

    /// standard constructor
    ReactionWithPhiScaling(Teuchos::RCP<ReactionInterface> reaction):
      reaction_(reaction)
    {};

    /// destructor
    virtual ~ReactionWithPhiScaling() { ; }

    /// initialization (to be called by derived classes)
    virtual void Initialize(
        int numscal,                         //!< number of scalars
        const std::vector<double>& couprole  //!< coupling role vector
        ) {reaction_->Initialize(numscal,couprole);};

    /// check for initialization
    bool IsInit() const {return reaction_->IsInit();};

    ///helper for calculating advanced reaction terms
    virtual double CalcReaBodyForceTerm(
        int k,                               //!< current scalar id
        int numscal,                         //!< number of scalars
        const std::vector<double>& phinp,    //!< scalar values at t_(n+1)
        const std::vector<double>& couprole, //!< coupling role vector
        double scale_reac,                   //!< scaling factor for reaction term (= reaction coefficient * stoichometry)
        double scale_phi                     //!< scaling factor for scalar values (used for reference concentrations)
        );

    ///helper for calculating advanced reaction term derivatives
    virtual void CalcReaBodyForceDeriv(
        int k,                               //!< current scalar id
        int numscal,                         //!< number of scalars
        std::vector<double>& derivs,         //!< vector with derivatives (to be filled)
        const std::vector<double>& phinp,    //!< scalar values at t_(n+1)
        const std::vector<double>& couprole, //!< coupling role vector
        double scale_reac,                   //!< scaling factor for reaction term (= reaction coefficient * stoichometry)
        double scale_phi                     //!< scaling factor for scalar values (used for reference concentrations)
        ) ;

    ///helper for calculating advanced reaction terms
    virtual double CalcReaBodyForceTerm(
        int k,                                                           //!< current scalar id
        int numscal,                                                     //!< number of scalars
        const std::vector<double>& phinp,                                //!< scalar values at t_(n+1)
        const std::vector<std::pair<std::string,double> >& constants,    //!< vector containing values which are independent of the scalars
        const std::vector<double>& couprole,                             //!< coupling role vector
        double scale_reac,                                               //!< scaling factor for reaction term (= reaction coefficient * stoichometry)
        double scale_phi                                                 //!< scaling factor for scalar values (used for reference concentrations)
        );

    ///helper for calculating advanced reaction term derivatives
    virtual void CalcReaBodyForceDeriv(
        int k,                                                           //!< current scalar id
        int numscal,                                                     //!< number of scalars
        std::vector<double>& derivs,                                     //!< vector with derivatives (to be filled)
        const std::vector<double>& phinp,                                //!< scalar values at t_(n+1)
        const std::vector<std::pair<std::string,double> >& constants,    //!< vector containing values which are independent of the scalars
        const std::vector<double>& couprole,                             //!< coupling role vector
        double scale_reac,                                               //!< scaling factor for reaction term (= reaction coefficient * stoichometry)
        double scale_phi                                                 //!< scaling factor for scalar values (used for reference concentrations)
        );

  protected:

    /// actual reaction
    Teuchos::RCP<ReactionInterface> reaction_;

    /// helper for scaling
    virtual std::vector<double> ModifyPhi(
        const std::vector<double>& phinp,
        double scale_phi
        );

  };

  /*----------------------------------------------------------------------*
   * **********************************************************************
   *----------------------------------------------------------------------*/

  //! base class for reaction coupling kinetics
  class ReactionBase : public ReactionInterface
  {

  public:

    /// standard constructor
    ReactionBase():isinit_(false){};

    /// destructor
    virtual ~ReactionBase() { ; }

    /// initialization (to be called by derived classes)
    virtual void Initialize(
        int numscal,                         //!< number of scalars
        const std::vector<double>& couprole  //!< coupling role vector
        ) {isinit_=true;};

    /// check for initialization
    bool IsInit() const {return isinit_;};

    ///helper for calculating advanced reaction terms
    virtual double CalcReaBodyForceTerm(
        int k,                               //!< current scalar id
        int numscal,                         //!< number of scalars
        const std::vector<double>& phinp,    //!< scalar values at t_(n+1)
        const std::vector<double>& couprole, //!< coupling role vector
        double scale_reac,                   //!< scaling factor for reaction term (= reaction coefficient * stoichometry)
        double scale_phi                     //!< scaling factor for scalar values (used for reference concentrations)
        ) ;

    ///helper for calculating advanced reaction term derivatives
    virtual void CalcReaBodyForceDeriv(
        int k,                               //!< current scalar id
        int numscal,                         //!< number of scalars
        std::vector<double>& derivs,         //!< vector with derivatives (to be filled)
        const std::vector<double>& phinp,    //!< scalar values at t_(n+1)
        const std::vector<double>& couprole, //!< coupling role vector
        double scale_reac,                   //!< scaling factor for reaction term (= reaction coefficient * stoichometry)
        double scale_phi                     //!< scaling factor for scalar values (used for reference concentrations)
        ) ;

    ///helper for calculating advanced reaction terms
    virtual double CalcReaBodyForceTerm(
        int k,                                                           //!< current scalar id
        int numscal,                                                     //!< number of scalars
        const std::vector<double>& phinp,                                //!< scalar values at t_(n+1)
        const std::vector<std::pair<std::string,double> >& constants,    //!< vector containing values which are independent of the scalars
        const std::vector<double>& couprole,                             //!< coupling role vector
        double scale_reac,                                               //!< scaling factor for reaction term (= reaction coefficient * stoichometry)
        double scale_phi                                                 //!< scaling factor for scalar values (used for reference concentrations)
        )
    {
      //default: ignore constants
      return CalcReaBodyForceTerm(k,numscal,phinp,couprole,scale_reac,scale_phi);
    };

    ///helper for calculating advanced reaction term derivatives
    virtual void CalcReaBodyForceDeriv(
        int k,                                                           //!< current scalar id
        int numscal,                                                     //!< number of scalars
        std::vector<double>& derivs,                                     //!< vector with derivatives (to be filled)
        const std::vector<double>& phinp,                                //!< scalar values at t_(n+1)
        const std::vector<std::pair<std::string,double> >& constants,    //!< vector containing values which are independent of the scalars
        const std::vector<double>& couprole,                             //!< coupling role vector
        double scale_reac,                                               //!< scaling factor for reaction term (= reaction coefficient * stoichometry)
        double scale_phi                                                 //!< scaling factor for scalar values (used for reference concentrations)
        )
    {
      //default: ignore constants
      CalcReaBodyForceDeriv(k,numscal,derivs,phinp,couprole,scale_reac,scale_phi);
      return;
    };

  protected:

    ///helper for calculating advanced reaction terms
    virtual double CalcReaBodyForceTerm(
        int k,                               //!< current scalar id
        int numscal,                         //!< number of scalars
        const std::vector<double>& phinp,    //!< scalar values at t_(n+1)
        const std::vector<double>& couprole, //!< coupling role vector
        double scale_reac                    //!< scaling factor for reaction term (= reaction coefficient * stoichometry)
        ) =0;

    ///helper for calculating advanced reaction term derivatives
    virtual void CalcReaBodyForceDeriv(
        int k,                               //!< current scalar id
        int numscal,                         //!< number of scalars
        std::vector<double>& derivs,         //!< vector with derivatives (to be filled)
        const std::vector<double>& phinp,    //!< scalar values at t_(n+1)
        const std::vector<double>& couprole, //!< coupling role vector
        double scale_reac                    //!< scaling factor for reaction term (= reaction coefficient * stoichometry)
        ) =0;

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

  private:

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
        double scale_reac                    //!< scaling factor for reaction term (= reaction coefficient * stoichometry)
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

  private:

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
        const int k,                         //!< current scalar id
        int numscal,                         //!< number of scalars
        std::vector<double>& derivs,         //!< vector with derivatives (to be filled)
        const std::vector<double>& phinp,    //!< scalar values at t_(n+1)
        const std::vector<double>& couprole, //!< coupling role vector
        double scale_reac                    //!< scaling factor for reaction term (= reaction coefficient * stoichometry)
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

  private:

    ///helper for calculating advanced reaction terms
    double CalcReaBodyForceTerm(
        const int k,                         //!< current scalar id
        int numscal,                         //!< number of scalars
        const std::vector<double>& phinp,    //!< scalar values at t_(n+1)
        const std::vector<double>& couprole, //!< coupling role vector
        double scale_reac                    //!< scaling factor for reaction term (= reaction coefficient * stoichometry)
        ) ;

    ///helper for calculating advanced reaction term derivatives
    void CalcReaBodyForceDeriv(
        const int k,                         //!< current scalar id
        int numscal,                         //!< number of scalars
        std::vector<double>& derivs,         //!< vector with derivatives (to be filled)
        const std::vector<double>& phinp,    //!< scalar values at t_(n+1)
        const std::vector<double>& couprole, //!< coupling role vector
        double scale_reac                    //!< scaling factor for reaction term (= reaction coefficient * stoichometry)
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
    virtual double CalcReaBodyForceTerm(
        int k,                                                           //!< current scalar id
        int numscal,                                                     //!< number of scalars
        const std::vector<double>& phinp,                                //!< scalar values at t_(n+1)
        const std::vector<std::pair<std::string,double> >& constants,    //!< vector containing values which are independent of the scalars
        const std::vector<double>& couprole,                             //!< coupling role vector
        double scale_reac,                                               //!< scaling factor for reaction term (= reaction coefficient * stoichometry)
        double scale_phi                                                 //!< scaling factor for scalar values (used for reference concentrations)
        );

    ///helper for calculating advanced reaction term derivatives
    virtual void CalcReaBodyForceDeriv(
        int k,                                                           //!< current scalar id
        int numscal,                                                     //!< number of scalars
        std::vector<double>& derivs,                                     //!< vector with derivatives (to be filled)
        const std::vector<double>& phinp,                                //!< scalar values at t_(n+1)
        const std::vector<std::pair<std::string,double> >& constants,    //!< vector containing values which are independent of the scalars
        const std::vector<double>& couprole,                             //!< coupling role vector
        double scale_reac,                                               //!< scaling factor for reaction term (= reaction coefficient * stoichometry)
        double scale_phi                                                 //!< scaling factor for scalar values (used for reference concentrations)
        );

  private:

    ///helper for calculating advanced reaction terms
    double CalcReaBodyForceTerm(
        const int k,                         //!< current scalar id
        int numscal,                         //!< number of scalars
        const std::vector<double>& phinp,    //!< scalar values at t_(n+1)
        const std::vector<double>& couprole, //!< coupling role vector
        double scale_reac                    //!< scaling factor for reaction term (= reaction coefficient * stoichometry)
        ) ;

    ///helper for calculating advanced reaction term derivatives
    void CalcReaBodyForceDeriv(
        const int k,                         //!< current scalar id
        int numscal,                         //!< number of scalars
        std::vector<double>& derivs,         //!< vector with derivatives (to be filled)
        const std::vector<double>& phinp,    //!< scalar values at t_(n+1)
        const std::vector<double>& couprole, //!< coupling role vector
        double scale_reac                    //!< scaling factor for reaction term (= reaction coefficient * stoichometry)
        ) ;

    ///helper for calculating advanced reaction terms
    double CalcReaBodyForceTerm(
        const int k,                                                     //!< current scalar id
        int numscal,                                                     //!< number of scalars
        const std::vector<double>& phinp,                                //!< scalar values at t_(n+1)
        const std::vector<std::pair<std::string,double> >& constants,    //!< vector containing values which are independent of the scalars
        const std::vector<double>& couprole,                             //!< coupling role vector
        double scale_reac                                                //!< scaling factor for reaction term (= reaction coefficient * stoichometry)
        ) ;

    ///helper for calculating advanced reaction term derivatives
    void CalcReaBodyForceDeriv(
        const int k,                                                     //!< current scalar id
        int numscal,                                                     //!< number of scalars
        std::vector<double>& derivs,                                     //!< vector with derivatives (to be filled)
        const std::vector<double>& phinp,                                //!< scalar values at t_(n+1)
        const std::vector<std::pair<std::string,double> >& constants,    //!< vector containing values which are independent of the scalars
        const std::vector<double>& couprole,                             //!< coupling role vector
        double scale_reac                                                //!< scaling factor for reaction term (= reaction coefficient * stoichometry)
        ) ;

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
