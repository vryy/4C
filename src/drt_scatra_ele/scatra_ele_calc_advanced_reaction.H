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


/*!
\brief Implementation of advanced deterministic reaction rate equations


<h3>Homogeneous Scatra Coupling</h3>

Note for the implementation of the homogenous scatra coupling:
assume the following reaction: 1*A + 2*B  --> 3*C with reaction coefficient 4.0

If we assume the reaction is depending on the product of all
reactants (this corresponds to couplingtype "simple_multiplicative"),
the corresponding equations are: \partial_t A = -(4*1*B)*A  (negative since reactant)
                                 \partial_t B = -(4*2*A)*B  (negative since reactant)
                                 \partial_t C = + 4*3*A*B   (positive since product)

This equation is achieved in BACI via the MAT_scatra_reaction material:
----------------------------------------------------------MATERIALS
MAT 1 MAT_matlist_reactions LOCAL No NUMMAT 3 MATIDS 2 4 5 NUMREAC 1 REACIDS 3 END //collect Concentrations
MAT 2 MAT_scatra DIFFUSIVITY 0.0
MAT 4 MAT_scatra DIFFUSIVITY 0.0
MAT 5 MAT_scatra DIFFUSIVITY 0.0
MAT 3 MAT_scatra_reaction NUMSCAL 3 STOICH -1 -2 3 REACOEFF 4.0 COUPLING simple_multiplicative

Implementation is of form: \partial_t c_i + K_i(c)*c_i = f_i(c), were f_i(c) is supposed not to linear depend on c_i
hence we have to calculate and set K(c)=(4*B;8*A;0) and f(c)=(0;0;12*A*B) and corresponding derivatives.


<h3>Special Reaction Kinetics</h3>

Note that for the implementation of the michaelis-menten scatra coupling feature:
reactant A promotes C and reactant B influences only until certain limit

MAT 1 MAT_matlist_reactions LOCAL No NUMMAT 3 MATIDS 2 4 5 NUMREAC 1 REACIDS 3 END //collect Concentrations
MAT 2 MAT_scatra DIFFUSIVITY 0.0
MAT 4 MAT_scatra DIFFUSIVITY 0.0
MAT 5 MAT_scatra DIFFUSIVITY 0.0
MAT 3 MAT_scatra_reaction NUMSCAL 3 STOICH -1 -1 1 REACOEFF 5.0 COUPLING michaelis_menten ROLE -1 3 0

The corresponding equations are
            \partial_t A = -5*A*(B/(3+B))
            \partial_t B = -5*A*(B/(3+B))
            \partial_t C =  5*A*(B/(3+B))

Implementation is of form: \partial_t c_i + K_i(c)*c_i = f_i(c), were f_i(c) is supposed not to linear depend on c_i
hence we have to calculate and set K(c)=(5*(B/(3+B));5*A*(1/(3+B));0) and f(c)=(0;0;5*A*(B/(3+B))) and corresponding derivatives.

Thereby ROLE does describe of how to build the reaction term (negative value -1.3: simply multiply by A,
positive value 3.2: multiply by B/(B+3.2) ) and STOICH does describe on which scalar the reaction term should be applied.
Here another example:

MAT 3 MAT_scatra_reaction NUMSCAL 4 STOICH +2 -1 0 0  REACOEFF 2.1 COUPLING michaelis_menten ROLE -1 0 1.2 3.0

The corresponding equations are
           \partial_t A =  2*2.1*A*C/(C+1.2)*D/(D+3.0)
           \partial_t B = -1*2.1*A*C/(C+1.2)*D/(D+3.0)
           \partial_t C =  0
           \partial_t D =  0

Note for the implementation of the constant scatra coupling feature:
Product A is constantly produced

MAT 1 MAT_matlist_reactions LOCAL No NUMMAT 1 MATIDS 2 NUMREAC 1 REACIDS 3 END //collect Concentrations
MAT 2 MAT_scatra DIFFUSIVITY 0.0
MAT 3 MAT_scatra_reaction NUMSCAL 1 STOICH 2 REACOEFF 5.0 COUPLING constant

The corresponding equation is:
            \partial_t A = 5*2

Implementation is of form: \partial_t c_i + K_i(c)*c_i = f_i(c), were f_i(c) is supposed not to depend linearly on c_i
hence we have to calculate and set K(c)=(0) and f(c)=(5*2) and zero derivatives.

<h3>Power Law Kinetics</h3>

Note for the implementation of the power law feature:
assume the following reaction: 1*A + 2*B  --> 3*C with reaction coefficient 4.0

If we assume the reaction is depending on the product of all reactants via a power law
(this corresponds to couplingtype "power_multiplicative"),

The corresponding equations are: \partial_t A = -(4*1*B^2)*A^3  (negative since reactant)
                                 \partial_t B = -(4*2*A^3)*B^2  (negative since reactant)
                                 \partial_t C = + 4*3*A^3*B^2   (positive since product)

This equation is in BACI achieved by the MAT_scatra_reaction material:
----------------------------------------------------------MATERIALS
MAT 1 MAT_matlist_reactions LOCAL No NUMMAT 3 MATIDS 2 4 5 NUMREAC 1 REACIDS 3 END //collect Concentrations
MAT 2 MAT_scatra DIFFUSIVITY 0.0
MAT 4 MAT_scatra DIFFUSIVITY 0.0
MAT 5 MAT_scatra DIFFUSIVITY 0.0
MAT 3 MAT_scatra_reaction NUMSCAL 3 STOICH -1 -2 3 REACOEFF 4.0 COUPLING power_multiplicative ROLE 3 2 0

Implementation is of form: \partial_t c_i + K_i(c)*c_i = f_i(c), were f_i(c) is supposed not to linear depend on c_i
hence we have to calculate and set K(c)=(4*B;8*A;0) and f(c)=(0;0;12*A*B) and corresponding derivatives.

<h3>Kinetics Defined by Function</h3>

Note for the implementation of the reaction-by-function feature:
assume the following reaction: 1*A + 2*B  --> 3*C with reaction coefficient 4.0

The reaction is defined by a function from the input file
(this corresponds to couplingtype "by_function").

It can reproduce all other implemented reaction types (and more). Also, linearizations are computed automatically
However, keep in mind, that it might be little slower ...

For instance if you want to model mechaelis-menten type kinematics of the form

           \partial_t A =  2*2.1*A*C/(C+1.2)*D/(D+3.0)
           \partial_t B = -1*2.1*A*C/(C+1.2)*D/(D+3.0)
           \partial_t C =  0
           \partial_t D =  0

This equation is in BACI achieved by the MAT_scatra_reaction material combined with the reaction by function feature as:
----------------------------------------------------------MATERIALS
MAT 1 MAT_matlist_reactions LOCAL No NUMMAT 3 MATIDS 2 4 5 NUMREAC 1 REACIDS 3 END //collect Concentrations
MAT 2 MAT_scatra DIFFUSIVITY 0.0
MAT 4 MAT_scatra DIFFUSIVITY 0.0
MAT 5 MAT_scatra DIFFUSIVITY 0.0
MAT 3 MAT_scatra_reaction NUMSCAL 4 STOICH +2 -1 0 0  REACOEFF 2.1 COUPLING by_function ROLE 1 1 0 0
-------------------------------------------------------------FUNCT1
FUNCT 1 COMPONENT 0 VAREXPR 0 0 0 FUNCTION phi1*phi3/(phi3+1.2)*phi4/(phi4+3.0)

Note: You have to use the VAREXPR function in order to be able to read the 'phix' variables.
      The function is defined within the role list

<h3>REACSTART Feature</h3>

Note for the implementation of the reacstart feature:
Assume concentration A is reproducing with reaction coefficient 1.0 and if the concentration
exceeds some threshold 2.0 if starts to react A->3*B with reacion coefficient 4.0.

The corresponding equations are:
            \partial_t A = -(-1.0)*A - 4.0*(A - 2.0)_{+} (first termn postive, since equivalent as reactant with negative reaction coefficient)
            \partial_t B = 3.0*4.0 (A - 2.0)_{+}   (positive since product)

This equation is achieved in BACI via the boundary condition:
MAT 1 MAT_matlist_reactions LOCAL No NUMMAT 1 MATIDS 2 NUMREAC 2 REACIDS 3 4 END //collect Concentrations
MAT 2 MAT_scatra DIFFUSIVITY 0.0
MAT 3 MAT_scatra_reaction NUMSCAL 2 STOICH -1 0 REACOEFF -1.0 COUPLING simple_multiplicative
MAT 4 MAT_scatra_reaction NUMSCAL 2 STOICH -1 3 REACOEFF 4.0 COUPLING simple_multiplicative REACSTART 2.0

Implementation is of form: \partial_t c_i + K_i(c)*c_i = f_i(c), were f_i(c) is supposed not to linear depend on c_i
hence we have to calculate and set K(c)=(-A + 4*(A-2)_{+};0) and f(c)=(0;12*(A-2)_{+}) and corresponding derivatives.


\author thon (thon@mhpc.mw.tum.de)
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
