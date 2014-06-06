/*----------------------------------------------------------------------*/
/*!
 \file scatra_ele_calc_advanced_reaction.H

 \brief main file containing routines for calculation of scatra element with advanced reaction terms


 <pre>
   Maintainer: Moritz Thon
               thon@mhpc.mw.tum.de
               http://www.lnm.mw.tum.de
 </pre>
 *----------------------------------------------------------------------*/


#ifndef SCATRA_ELE_CALC_ADVANCED_REACTION_H_
#define SCATRA_ELE_CALC_ADVANCED_REACTION_H_

#include "scatra_ele_calc.H"

namespace DRT
{
namespace ELEMENTS
{

enum reaction_coupling
  {
    simple_multiplicative
  };

template<DRT::Element::DiscretizationType distype>
class ScaTraEleCalcAdvReac : public virtual ScaTraEleCalc<distype>
{
protected:
    /// (private) protected constructor, since we are a Singleton.
  ScaTraEleCalcAdvReac(const int numdofpernode,const int numscal);

private:

  typedef ScaTraEleCalc<distype> my;

public:
  virtual ~ScaTraEleCalcAdvReac() {}

  /// Singleton access method
  static ScaTraEleCalcAdvReac<distype> * Instance(
    const int numdofpernode,
    const int numscal,
    bool create=true );

  /// called upon destruction
  virtual void Done();

protected:

  //! evaluate material
  virtual void Materials(
    const Teuchos::RCP<const MAT::Material> material, //!< pointer to current material
    const int                               k,        //!< id of current scalar
    double&                                 densn,    //!< density at t_(n)
    double&                                 densnp,   //!< density at t_(n+1) or t_(n+alpha_F)
    double&                                 densam,   //!< density at t_(n+alpha_M)
    Teuchos::RCP<ScaTraEleDiffManager>      diffmanager,  //!< diffusion manager handling diffusivity / diffusivities (in case of systems) or (thermal conductivity/specific heat) in case of loma
    Teuchos::RCP<ScaTraEleReaManager>       reamanager,   //!< reaction manager
    double&                                 visc,      //!< fluid viscosity
    const int                               iquad=-1   //!< id of current gauss point (default = -1)
  );


  //! material ScaTra
  virtual void MatScaTra(
    const Teuchos::RCP<const MAT::Material> material, //!< pointer to current material
    const int                               k,        //!< id of current scalar
    double&                                 densn,    //!< density at t_(n)
    double&                                 densnp,   //!< density at t_(n+1) or t_(n+alpha_F)
    double&                                 densam,   //!< density at t_(n+alpha_M)
    Teuchos::RCP<ScaTraEleDiffManager>      diffmanager,  //!< diffusion manager handling diffusivity / diffusivities (in case of systems) or (thermal conductivity/specific heat) in case of loma
    Teuchos::RCP<ScaTraEleReaManager>       reamanager,   //!< reaction manager
    double&                                 visc,      //!< fluid viscosity
    const int                               iquad=-1   //!< id of current gauss point (default = -1)
  );


  //! material BioFilm
  virtual void MatBioFilm(
      const Teuchos::RCP<const MAT::Material> material, //!< pointer to current material
      const int                               k,        //!< id of current scalar
      double&                                 densn,    //!< density at t_(n)
      double&                                 densnp,   //!< density at t_(n+1) or t_(n+alpha_F)
      double&                                 densam,   //!< density at t_(n+alpha_M)
      Teuchos::RCP<ScaTraEleDiffManager>      diffmanager,  //!< diffusion manager handling diffusivity / diffusivities (in case of systems) or (thermal conductivity/specific heat) in case of loma
      Teuchos::RCP<ScaTraEleReaManager>       reamanager,   //!< reaction manager
      double&                                 visc,      //!< fluid viscosity
      const int                               iquad=-1   //!< id of current gauss point (default = -1)
  );


  //! material ScatraGrowthScd
  virtual void MatGrowthScd(
      const Teuchos::RCP<const MAT::Material> material, //!< pointer to current material
      const int                               k,        //!< id of current scalar
      double&                                 densn,    //!< density at t_(n)
      double&                                 densnp,   //!< density at t_(n+1) or t_(n+alpha_F)
      double&                                 densam,   //!< density at t_(n+alpha_M)
      Teuchos::RCP<ScaTraEleDiffManager>      diffmanager,  //!< diffusion manager handling diffusivity / diffusivities (in case of systems) or (thermal conductivity/specific heat) in case of loma
      Teuchos::RCP<ScaTraEleReaManager>       reamanager,   //!< reaction manager
      double&                                 visc,      //!< fluid viscosity
      const int                               iquad=-1   //!< id of current gauss point (default = -1)
  );


  virtual void GetRhsInt(
    double&      rhsint,  //!< rhs containing bodyforce at Gauss point
    const double densnp,  //!< density at t_(n+1)
    const int    k        //!< index of current scalar
  );


  //! check if homogenous scatra couplung condition exists and if, read it
  virtual bool IsCoupledAndRead(
      const DRT::Discretization&            ScaTraDiscretization //!<discretisation of the ScaTra field
    );


  //!Calculate K(c)
  virtual double CalcReaCoeff(
      const int                        k                       //!< id of current scalar
      );


  virtual double CalcReaCoeffFac(
      const std::vector<int>                    stoich,                  //!<stoichometrie of current condition
      const DRT::ELEMENTS::reaction_coupling    couplingtype,            //!<type of coupling the stoichometry coefficients
      const int                                 k                       //!< id of current scalar
  );


  //! calculate \frac{partial}{\partial c} K(c)
  virtual double CalcReaCoeffDerivMatrix(
      const int                 k,                  //!< id of current scalar
      const int                 j                   //!< concentration to be derived to
      );


  virtual double CalcReaCoeffDerivFac(
      const std::vector<int>                  stoich,                  //!<stoichometrie of current condition
      const DRT::ELEMENTS::reaction_coupling  couplingtype,            //!<type of coupling the stoichometry coefficients
      const int                               toderive,                //!<concentration to be derived to
      const int                               k                       //!< id of current scalar
  );


  //! calculate f(c)
  virtual double CalcReaBodyForceTerm(
      const int                                 k                       //!< id of current scalar
);


  virtual double CalcReaBodyForceTermFac(
      const std::vector<int>                      stoich,                 //!<stoichometrie of current condition
      const DRT::ELEMENTS::reaction_coupling      couplingtype            //!<type of coupling the stoichometry coefficients
  );


  //! calculate \frac{partial}{\partial c} f(c)
  virtual double CalcReaBodyForceDerivMatrix(
      const int                 k,                  //!< id of current scalar
      const int                 j                   //!< concentration to be derived to
);


  virtual double CalcReaBodyForceDerivFac(
      const std::vector<int>                    stoich,                  //!<stoichometrie of current condition
      const DRT::ELEMENTS::reaction_coupling    couplingtype,            //!<type of coupling the stoichometry coefficients
      const int                                 toderive                 //!<concentration to be derived to
  );


  //! Calculate (K(c)-reacstart(c)./c)_{+}
  virtual void ReacStartForReaCoeff(
  const int                               k,            //!< id of current scalar
  const int                               condnum,       //!< id of current condition
  double&                                 value         //!< current reaction value
  );

  //! Calculate \frac{\partial}{\partial_c} (K(c)-reacstart(c)./c)_{+}
  void ReacStartForReaCoeffDeriv(
    const int                               k,              //!< id of current scalar
    const int                               toderive,       //!<concentration to be derived to
    const int                               condnum,        //!< id of current condition
    double&                                 value,          //!< current reaction value
    const std::vector<int>                  stoich,         //!<stoichometrie of current condition
    const DRT::ELEMENTS::reaction_coupling  couplingtype    //!<type of coupling the stoichometry coefficients
    );


  //! Calculate (f(c)-reacstart(c))_{+}
  virtual void ReacStartForReaBF(
  const int                               k,            //!< id of current scalar
  const int                               condnum,       //!< id of current condition
  double&                                 value         //!< current reaction value
  );

  //! Calculate \frac{\partial}{\partial_c} (f(c)-reacstart(c))_{+}
  void ReacStartForReaBFDeriv(
    const int                               k,              //!< id of current scalar
    const int                               toderive,       //!<concentration to be derived to
    const int                               condnum,        //!< id of current condition
    double&                                 value,          //!< current reaction value
    const std::vector<int>                  stoich,         //!<stoichometrie of current condition
    const DRT::ELEMENTS::reaction_coupling  couplingtype    //!<type of coupling the stoichometry coefficients
    );


  virtual int Evaluate(DRT::ELEMENTS::Transport*     ele,
                         Teuchos::ParameterList&       params,
                         DRT::Discretization &         discretization,
                         const std::vector<int> &      lm,
                         Epetra_SerialDenseMatrix&     elemat1_epetra,
                         Epetra_SerialDenseMatrix&     elemat2_epetra,
                         Epetra_SerialDenseVector&     elevec1_epetra,
                         Epetra_SerialDenseVector&     elevec2_epetra,
                         Epetra_SerialDenseVector&     elevec3_epetra);


  //! calculation of reactive element matrix
  virtual void CalcMatReact(
    Epetra_SerialDenseMatrix&          emat,        //!< element matrix to be filled
    const int                          k,           //!< index of current scalar
    const double                       timefacfac,  //!< domain-integration factor times time-integration factor
    const double                       timetaufac,  //!< domain-integration factor times time-integration factor times tau
    const double                       taufac,      //!< domain-integration factor times tau
    const double                       densnp,      //!< density at time_(n+1)
    const double                       phinp,       //!< scalar at time_(n+1)
    Teuchos::RCP<ScaTraEleReaManager>  reamanager,  //!< reaction manager
    const LINALG::Matrix<my::nen_,1>&      conv,        //!< convective operator
    const LINALG::Matrix<my::nen_,1>&      sgconv,      //!< subgrid-scale convective operator
    const LINALG::Matrix<my::nen_,1>&      diff         //!< laplace term
  );

  void SetAdvancedReactionTerms(Teuchos::RCP<ScaTraEleReaManager>       reamanager,
                                const int                               k,
                                const double                            scale=1.0);

  bool iscoupled_;                           /// Flag for homogenous scatra coupling volume condition
  bool isinit_;                              /// Flag if is checked for coupling condition or Material is allready done
  int numcond_;                              /// number of homogenous scatra coupling conditions
  std::vector<std::vector<int> > stoich_;    /// vector containing the stoichometrie
  std::vector<double> reaconst_;             ///constant by which the reaction coefficient is multiplied
  std::vector<DRT::ELEMENTS::reaction_coupling> couplingtype_;///how to calculate the different couplings
  std::vector<double> reacstart_;             ///constant by which the reaction starts


}; //end ScaTraEleCalcAdvReac



/// Scatra reaction manager for Advanced_Reaction
/*!
  This class keeps all Advanced_Reaction-specific reaction parameter needed for the evaluation of an element. The
  ScaTraEleReaManagerAdvReac is derived from the standard ScaTraEleReaManager.
*/
class ScaTraEleReaManagerAdvReac : public ScaTraEleReaManager
  {
  public:

    ScaTraEleReaManagerAdvReac(int numscal): ScaTraEleReaManager(numscal),
      reacoeff_(numscal,0.0),
      reacoeffderivmatrix_(numscal,std::vector<double>(numscal,0.0 )), // size of matrix + initialized to zero
      reabodyforce_(numscal,0.0), // size of vector + initialized to zero
      reabodyforcederiv_(numscal,std::vector<double>(numscal,0.0 )),
      include_me_(false)
      {return;}

    virtual ~ScaTraEleReaManagerAdvReac(){};

    //! @name set routines

    //! Set the reaction coefficient
    virtual void SetReaCoeff(const double reacoeff, const int k)
    {
      reacoeff_[k] = reacoeff;
      if( reacoeff!=0.0 )
        include_me_=true;

      return;
    }

    virtual void SetReaCoeffDerivMatrix(const double reacoeffderivmatrix, const int k, const int j)
    {
      (reacoeffderivmatrix_[k])[j] = reacoeffderivmatrix;
      return;
    }

    //! Set the derivative of reaction coefficient
    virtual void SetReaBodyForce(const double reabodyforce, const int k)
    {
      reabodyforce_[k] = reabodyforce;
      return;
    }

    virtual void SetReaBodyForceDerivMatrix(const double reabodyforcederiv, const int k, const int j)
    {
      (reabodyforcederiv_[k])[j] = reabodyforcederiv;
      return;
    }

    //! @name access routines

    //! Return the reaction coefficient
    virtual double GetReaCoeff(const int k){return reacoeff_[k];}

    //! Return the reaction coefficient
    virtual double GetReaCoeffDerivMatrix(const int k , const int j){return (reacoeffderivmatrix_[k])[j];}

    //! Return the reaction coefficient
    virtual double GetReaBodyForce(const int k){return reabodyforce_[k];}

    //! Return the reaction coefficient
    virtual double GetReaBodyForceDerivMatrix(const int k, const int j){return (reabodyforcederiv_[k])[j];}

    //! return flag: reaction activated
    virtual bool Active(){return include_me_;}

  private:

    //! scalar reaction coefficient
    std::vector<double> reacoeff_;

    //! reaction derivative for coupled problems
    std::vector<std::vector<double> > reacoeffderivmatrix_;

    //! scalar reaction coefficient
    std::vector<double> reabodyforce_;

    //! scalar reaction coefficient
    std::vector<std::vector<double> > reabodyforcederiv_;

    //! flag for reaction
    bool include_me_;
  };


}

}


#endif /* SCATRA_ELE_CALC_ADVANCED_REACTION_H_ */
