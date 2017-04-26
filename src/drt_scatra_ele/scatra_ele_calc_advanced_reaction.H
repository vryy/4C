/*----------------------------------------------------------------------*/
/*!
 \file scatra_ele_calc_advanced_reaction.H

 \brief main file containing routines for calculation of scatra element with advanced reaction terms


 <pre>
  \level 2

   \maintainer Moritz Thon
               thon@mhpc.mw.tum.de
               http://www.lnm.mw.tum.de
               089 - 289-10364
 </pre>
 *----------------------------------------------------------------------*/


#ifndef SCATRA_ELE_CALC_ADVANCED_REACTION_H_
#define SCATRA_ELE_CALC_ADVANCED_REACTION_H_

#include "scatra_ele_calc.H"

// forward declarations
namespace MAT
{
  class MatListReactions;
}

namespace DRT
{
namespace ELEMENTS
{
// forward declaration
class ScaTraEleReaManagerAdvReac;

/*
 This class calls all 'advanced reaction terms' calculations and applies them correctly. Thereby, no assumption on
 the shape of the (potentially nonlinear) reaction term f(c_1,...,c_n) have to be made. The actual calculations
 of the reaction term f(c_1,...,c_n) as well as all its linearizations \partial_c f(c) are done within the material
 MAT_matlist_reactions and MAT_scatra_reaction
 */
template<DRT::Element::DiscretizationType distype,int probdim=DRT::UTILS::DisTypeToDim<distype>::dim>
class ScaTraEleCalcAdvReac : public virtual ScaTraEleCalc<distype,probdim>
{
protected:
    /// (private) protected constructor, since we are a Singleton.
  ScaTraEleCalcAdvReac(const int numdofpernode,const int numscal,const std::string& disname);

private:

  typedef ScaTraEleCalc<distype,probdim> my;

public:
  virtual ~ScaTraEleCalcAdvReac(){return;};

  /// Singleton access method
  static ScaTraEleCalcAdvReac<distype,probdim> * Instance(
    const int numdofpernode,
    const int numscal,
    const std::string& disname,
    const ScaTraEleCalcAdvReac* delete_me = NULL  //!< creation/destruction indication
    );

  /// called upon destruction
  virtual void Done();

protected:

  //! set internal variables
  virtual void SetInternalVariablesForMatAndRHS();

  //! get the material parameters
  virtual void GetMaterialParams(
    const DRT::Element*  ele,       //!< the element we are dealing with
    std::vector<double>& densn,     //!< density at t_(n)
    std::vector<double>& densnp,    //!< density at t_(n+1) or t_(n+alpha_F)
    std::vector<double>& densam,    //!< density at t_(n+alpha_M)
    double&              visc,       //!< fluid viscosity
    const int            iquad=-1    //!< id of current gauss point (default = -1)
  );


  //! evaluate material
  virtual void Materials(
    const Teuchos::RCP<const MAT::Material> material, //!< pointer to current material
    const int                               k,        //!< id of current scalar
    double&                                 densn,    //!< density at t_(n)
    double&                                 densnp,   //!< density at t_(n+1) or t_(n+alpha_F)
    double&                                 densam,   //!< density at t_(n+alpha_M)
    double&                                 visc,      //!< fluid viscosity
    const int                               iquad=-1   //!< id of current gauss point (default = -1)
  );


  //! Get right hand side including reaction bodyforce term
  virtual void GetRhsInt(
    double&      rhsint,  //!< rhs containing bodyforce at Gauss point
    const double densnp,  //!< density at t_(n+1)
    const int    k        //!< index of current scalar
  );


  //! calculation of reactive element matrix
  virtual void CalcMatReact(
    Epetra_SerialDenseMatrix&          emat,        //!< element matrix to be filled
    const int                          k,           //!< index of current scalar
    const double                       timefacfac,  //!< domain-integration factor times time-integration factor
    const double                       timetaufac,  //!< domain-integration factor times time-integration factor times tau
    const double                       taufac,      //!< domain-integration factor times tau
    const double                       densnp,      //!< density at time_(n+1)
    const LINALG::Matrix<my::nen_,1>&      sgconv,  //!< subgrid-scale convective operator
    const LINALG::Matrix<my::nen_,1>&      diff     //!< laplace term
  );


  //! Set advanced reaction terms and derivatives
  virtual void SetAdvancedReactionTerms(
      const int                                 k,           //!< index of current scalar
      const Teuchos::RCP<MAT::MatListReactions> matreaclist, //!< index of current scalar
      const double* gpcoord                                  //!< current Gauss-point coordinates
  );

  //! evaluate shape functions and their derivatives at element center
  virtual double EvalShapeFuncAndDerivsAtEleCenter();

  //! array for shape function at element center
  LINALG::Matrix<my::nen_,1> funct_elementcenter_;

  //! get current Gauss-point coordinates
  virtual const double* GetGpCoord() const { return gpcoord_;}

  //! get reaction manager for advanced reaction
  Teuchos::RCP<ScaTraEleReaManagerAdvReac> ReaManager(){return Teuchos::rcp_static_cast<ScaTraEleReaManagerAdvReac>(my::reamanager_);};

private:

  //! current Gauss-point coordinates
  double gpcoord_[my::nsd_];

}; //end ScaTraEleCalcAdvReac


/// Scatra reaction manager for Advanced_Reaction
/*!
  This class keeps all advanced reaction terms specific sutffs needed for the evaluation of an element. The
  ScaTraEleReaManagerAdvReac is derived from the standard ScaTraEleReaManager.
*/
class ScaTraEleReaManagerAdvReac : public ScaTraEleReaManager
  {
  public:

    ScaTraEleReaManagerAdvReac(int numscal): ScaTraEleReaManager(numscal),
      reabodyforce_(numscal,0.0), // size of vector + initialized to zero
      reabodyforcederiv_(numscal,std::vector<double>(numscal,0.0 ))
      {return;}

    virtual ~ScaTraEleReaManagerAdvReac(){};

    //! @name set routines

    //! Clear everything and resize to length numscal
    void Clear(int numscal)
    {
      //clear base class
      ScaTraEleReaManager::Clear(numscal);
      //clear
      reabodyforce_.resize(0);
      reabodyforcederiv_.resize(0);
      //resize
      reabodyforce_.resize(numscal,0.0);
      reabodyforcederiv_.resize(numscal,std::vector<double>(numscal,0.0 ));
      return;
    }

    //! Add to the body force due to reaction
    void AddToReaBodyForce(const double reabodyforce, const int k)
    {
      reabodyforce_[k] += reabodyforce;
      if( reabodyforce!=0.0 )
        include_me_=true;

      return;
    }

    //! Return one line of the jacobian of the reaction vector
    std::vector<double>& GetReaBodyForceDerivVector(const int k){return reabodyforcederiv_[k];}

    //! Add to the derivative of the body force due to reaction
    void AddToReaBodyForceDerivMatrix(const double reabodyforcederiv, const int k, const int j)
    {
      (reabodyforcederiv_[k])[j] += reabodyforcederiv;
      return;
    }

    //@}

    //! @name access routines

    //! Return the reaction coefficient
    double GetReaBodyForce(const int k) const {return reabodyforce_[k];}

    //! Return the reaction coefficient
    double GetReaBodyForceDerivMatrix(const int k, const int j) const {return (reabodyforcederiv_[k])[j];}

    //! Return the stabilization coefficient
    virtual double GetStabilizationCoeff(const int k, const double phinp_k) const
    {
      double stabboeff = ScaTraEleReaManager::GetStabilizationCoeff(k,phinp_k);

      if (phinp_k>1.0e-10)
        stabboeff += fabs(reabodyforce_[k]/phinp_k);

      return stabboeff;
    }

    //@}

  private:
    //! @name protected variables

    //! scalar reaction coefficient
    std::vector<double> reabodyforce_;

    //! scalar reaction coefficient
    std::vector<std::vector<double> > reabodyforcederiv_;

    //@}
  };


}

}


#endif /* SCATRA_ELE_CALC_ADVANCED_REACTION_H_ */
