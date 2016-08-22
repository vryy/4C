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
#include "../drt_mat/scatra_reaction_mat.H"

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

  //! get the material parameters
  virtual void GetMaterialParams(
    const DRT::Element* ele,       //!< the element we are dealing with
    double&             densn,     //!< density at t_(n)
    double&             densnp,    //!< density at t_(n+1) or t_(n+alpha_F)
    double&             densam,    //!< density at t_(n+alpha_M)
    double&             visc,       //!< fluid viscosity
    const int           iquad=-1    //!< id of current gauss point (default = -1)
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


  //! material BioFilm
  void MatBioFilm(
      const Teuchos::RCP<const MAT::Material> material, //!< pointer to current material
      const int                               k,        //!< id of current scalar
      double&                                 densn,    //!< density at t_(n)
      double&                                 densnp,   //!< density at t_(n+1) or t_(n+alpha_F)
      double&                                 densam,   //!< density at t_(n+alpha_M)
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


  //! Set reac. body force, reaction coefficient and derivatives
  virtual void SetAdvancedReactionTerms(
      const int                                 k,          //!< index of current scalar
      const Teuchos::RCP<MAT::MatListReactions> matreaclist //!< index of current scalar
  );

  //! evaluate shape functions and their derivatives at element center
  virtual double EvalShapeFuncAndDerivsAtEleCenter();

  LINALG::Matrix<my::nen_,1> funct_elementcenter_;            //! array for shape function at element center

  //! get reaction manager for advanced reaction
  Teuchos::RCP<ScaTraEleReaManagerAdvReac> ReaManager(){return Teuchos::rcp_static_cast<ScaTraEleReaManagerAdvReac>(my::reamanager_);};
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
      reabodyforce_(numscal,0.0), // size of vector + initialized to zero
      reabodyforcederiv_(numscal,std::vector<double>(numscal,0.0 ))
      {return;}

    virtual ~ScaTraEleReaManagerAdvReac(){};

    //! @name set routines

    //! Clear everything and resize to length numscal
    virtual void Clear(int numscal)
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

    //! Set the reaction coefficient
    virtual void SetReaCoeff(const double reacoeff, const int k)
    {
      reacoeff_[k] = reacoeff;
      if( reacoeff!=0.0 )
        include_me_=true;

      return;
    }

    //! Add to the derivative of reaction coefficient
    virtual void AddToReaCoeffDerivMatrix(const double reacoeffderivmatrix, const int k, const int j)
    {
      (reacoeffderivmatrix_[k])[j] += reacoeffderivmatrix;
      return;
    }

    //! Add to the body force due to reaction
    virtual void AddToReaBodyForce(const double reabodyforce, const int k)
    {
      reabodyforce_[k] += reabodyforce;
      if( reabodyforce!=0.0 )
        include_me_=true;

      return;
    }

    //! Add to the derivative of the body force due to reaction
    virtual void AddToReaBodyForceDerivMatrix(const double reabodyforcederiv, const int k, const int j)
    {
      (reabodyforcederiv_[k])[j] += reabodyforcederiv;
      return;
    }

    //@}

    //! @name access routines

    //! Return the reaction coefficient
    virtual double GetReaBodyForce(const int k){return reabodyforce_[k];}

    //! Return the reaction coefficient
    virtual double GetReaBodyForceDerivMatrix(const int k, const int j){return (reabodyforcederiv_[k])[j];}

    //@}

  protected:

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
