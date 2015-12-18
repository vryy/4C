/*--------------------------------------------------------------------------*/
/*!
\file lubrication_ele_calc.H

\brief main file containing routines for calculation of lubrication element

<pre>
Maintainer: Andy Wirtz
            wirtz@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089-289-15270
</pre>
*/
/*--------------------------------------------------------------------------*/
#ifndef SRC_DRT_LUBRICATION_ELE_LUBRICATION_ELE_CALC_H_
#define SRC_DRT_LUBRICATION_ELE_LUBRICATION_ELE_CALC_H_


#include "../drt_lib/drt_element.H"

#include "../drt_fem_general/drt_utils_local_connectivity_matrices.H"
#include "../drt_fem_general/drt_utils_integration.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_fem_general/drt_utils_boundary_integration.H"
#include "../drt_lubrication_ele/lubrication_ele_interface.H"

#include "../drt_lubrication_ele/lubrication_ele_action.H"




namespace DRT
{
namespace ELEMENTS
{

  // forward declaration
  class LubricationEleParameter;

  class LubricationEleViscManager;
  template<int NSD,int NEN> class LubricationEleInternalVariableManager;

  /// Lubrication element implementation
  /*!
    This internal class keeps all the working arrays needed to
    calculate the Lubrication element. Additionally, the method Sysmat()
    provides a clean and fast element implementation.

    <h3>Purpose</h3>

    The idea is to separate the element maintenance (class Lubrication) from the
    mathematical contents (this class). There are different
    implementations of the Lubrication element, this is just one such
    implementation.

    The Lubrication element will allocate exactly one object of this class for all
    Lubrication elements with the same number of nodes in the mesh. This
    allows us to use exactly matching working arrays (and keep them
    around.)

    The code is meant to be as clean as possible. This is the only way
    to keep it fast. The number of working arrays has to be reduced to
    a minimum so that the element fits into the cache. (There might be
    room for improvements.)

    <h3>Usability</h3>

    The calculations are done by the Evaluate() method. There are two
    version. The virtual method that is inherited from LubricationEleInterface
    (and called from Lubrication) and the non-virtual one that does the actual
    work. The non-virtual Evaluate() method must be callable without an actual
    Lubrication object.
  */

  template<DRT::Element::DiscretizationType distype,int probdim=DRT::UTILS::DisTypeToDim<distype>::dim>
  class LubricationEleCalc : public LubricationEleInterface
  {
  protected:

    /// (private) protected constructor, since we are a Singleton.
    /// this constructor is called from a derived class
    /// -> therefore, it has to be protected instead of private
    LubricationEleCalc(const std::string& disname);

  public:

    virtual ~LubricationEleCalc() {};

    //! singleton access method
    static LubricationEleCalc<distype,probdim> * Instance(
        const std::string& disname,
        const LubricationEleCalc* delete_me = NULL
        );

    /// called upon destruction
    virtual void Done();

    /// In this class we do not define a static LubricationEle...* Instance
    /// since only derived child classes are free to be allocated!!

    /// Setup element evaluation
    virtual int SetupCalc(
        DRT::Element*               ele,
        DRT::Discretization&        discretization
        );

    /// Evaluate the element
     /*!
       Generic virtual interface function. Called via base pointer.
      */
     virtual int Evaluate(DRT::Element*                 ele,
                          Teuchos::ParameterList&       params,
                          DRT::Discretization&          discretization,
                          DRT::Element::LocationArray&  la,
                          Epetra_SerialDenseMatrix&     elemat1_epetra,
                          Epetra_SerialDenseMatrix&     elemat2_epetra,
                          Epetra_SerialDenseVector&     elevec1_epetra,
                          Epetra_SerialDenseVector&     elevec2_epetra,
                          Epetra_SerialDenseVector&     elevec3_epetra);

     //! evaluate action
     virtual int EvaluateAction(
         DRT::Element*                 ele,
         Teuchos::ParameterList&       params,
         DRT::Discretization&          discretization,
         const LUBRICATION::Action&    action,
         DRT::Element::LocationArray&  la,
         Epetra_SerialDenseMatrix&     elemat1_epetra,
         Epetra_SerialDenseMatrix&     elemat2_epetra,
         Epetra_SerialDenseVector&     elevec1_epetra,
         Epetra_SerialDenseVector&     elevec2_epetra,
         Epetra_SerialDenseVector&     elevec3_epetra
         );

     //! evaluate service routine
     int EvaluateService(
         DRT::Element*                 ele,
         Teuchos::ParameterList&       params,
         DRT::Discretization &         discretization,
         DRT::Element::LocationArray&  la,
         Epetra_SerialDenseMatrix&     elemat1_epetra,
         Epetra_SerialDenseMatrix&     elemat2_epetra,
         Epetra_SerialDenseVector&     elevec1_epetra,
         Epetra_SerialDenseVector&     elevec2_epetra,
         Epetra_SerialDenseVector&     elevec3_epetra
         );

     /*========================================================================*/
     //! @name static member variables
     /*========================================================================*/

     //! number of element nodes (nomenclature: T. Hughes, The finite element method)
     static const int nen_ = DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;

     //! number of space dimensions
     static const int nsd_ = probdim;
     //! space dimension of Lubrication element (only for flat domains nsd_ele_ = nsd_)
     static const int nsd_ele_ = DRT::UTILS::DisTypeToDim<distype>::dim;

  protected:

    /*========================================================================*/
    //! @name general framework
    /*========================================================================*/

    //! extract element based or nodal values
    //  return extracted values of prenp
    virtual void ExtractElementAndNodeValues(
      DRT::Element*                 ele,
      Teuchos::ParameterList&       params,
      DRT::Discretization&          discretization,
      DRT::Element::LocationArray&  la
    );

    //! calculate matrix and rhs. Here the whole thing is hidden.
    virtual void Sysmat(
      DRT::Element*                         ele,       //!< the element we are dealing with
      Epetra_SerialDenseMatrix&             emat,      //!< element matrix to calculate
      Epetra_SerialDenseVector&             erhs       //!< element rhs to calculate
    );

    //! get the body force
    virtual void BodyForce(
      const DRT::Element*  ele //!< the element we are dealing with
    );

    //! read element coordinates
    virtual void ReadElementCoordinates(
        const DRT::Element*     ele
    );

    //! evaluate shape functions and their derivatives at current integration point
    double EvalShapeFuncAndDerivsAtIntPoint(
      const DRT::UTILS::IntPointsAndWeights<nsd_ele_>& intpoints, //!< integration points
      const int                                        iquad      //!< id of current Gauss point
    );

    //! evaluate shape functions and their derivatives at current integration point
    double EvalShapeFuncAndDerivsInParameterSpace();

    //! set internal variables
    virtual void SetInternalVariablesForMatAndRHS();

    /*========================================================================*/
    //! @name routines for additional element evaluations (called from EvaluateAction)
    /*========================================================================*/

    //! calculate error of numerical solution with respect to analytical solution
    void CalErrorComparedToAnalytSolution(
      const DRT::Element*                   ele,        //!< the element we are dealing with
      Teuchos::ParameterList&               params,     //!< parameter list
      Epetra_SerialDenseVector&             errors      //!< vector containing L2-error norm
    );

    //! calculate pressure(s) and domain integral
    virtual void CalculatePressures(
      const DRT::Element*             ele,       //!< the element we are dealing with
      Epetra_SerialDenseVector&       pressures,   //!< pressure to be computed
      const bool                      inverting  //!< bool indicating inversion
    );


    /*========================================================================*/
    //! @name material and related functions
    /*========================================================================*/

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
      double&                                 densn,    //!< density at t_(n)
      double&                                 densnp,   //!< density at t_(n+1) or t_(n+alpha_F)
      double&                                 densam,   //!< density at t_(n+alpha_M)
      double&                                 visc,      //!< fluid viscosity
      const int                               iquad=-1   //!< id of current gauss point (default = -1)
    );

    //! material Lubrication
    virtual void MatLubrication(
      const Teuchos::RCP<const MAT::Material> material, //!< pointer to current material
      double&                                 densn,    //!< density at t_(n)
      double&                                 densnp,   //!< density at t_(n+1) or t_(n+alpha_F)
      double&                                 densam,   //!< density at t_(n+alpha_M)
      double&                                 visc,      //!< fluid viscosity
      const int                               iquad=-1   //!< id of current gauss point (default = -1)
    );
    /*========================================================================*/
    //! @name methods for evaluation of individual terms
    /*========================================================================*/

    //! calculate linearization of the Laplacian (weak form)
    void GetLaplacianWeakForm(
      double&                          val,   //!< value of linearization of weak laplacian
      const int                        vi,    //!< node index for the weighting function
      const int                        ui     //!< node index for the shape function
    )
    {
      val = 0.0;
      for (int j = 0; j<nsd_; j++)
      {
        val += derxy_(j, vi)*derxy_(j, ui);
      }
      return;
    };

    //! calculate the Laplacian (weak form)
    void GetLaplacianWeakFormRHS(
      double& val,                              //!< value of weak laplacian
      const LINALG::Matrix<nsd_,1>&    gradpre, //!< pressure gradient at gauss point
      const int vi                              //!< node index for the weighting function
    )
    {
      val = 0.0;
      for (int j = 0; j<nsd_; j++)
      {
        val += derxy_(j,vi)*gradpre(j);
      }
      return;
    };

    //! compute rhs containing bodyforce
    virtual void GetRhsInt(
      double&      rhsint,  //!< rhs containing bodyforce at integration point
      const double densnp  //!< density at t_(n+1)
    );

    //! calculation of diffusive element matrix
    virtual void CalcMatDiff(
      Epetra_SerialDenseMatrix&     emat,        //!< element matrix to be filled
      const double                  timefacfac   //!< domain-integration factor times time-integration factor
    );

    //! standard Galerkin transient, old part of rhs and source term
    void CalcRHSHistAndSource(
      Epetra_SerialDenseVector&     erhs,        //!< element vector to be filled
      const double                  fac,         //!< domain-integration factor
      const double                  rhsint       //!< rhs at Gauss point
    );

    //! standard Galerkin diffusive term on right hand side
    virtual void CalcRHSDiff(
      Epetra_SerialDenseVector&     erhs,        //!< element vector to be filled
      const double                  rhsfac       //!< time-integration factor for rhs times domain-integration factor
    );

    /*========================================================================*/
    //! @name parameter lists
    /*========================================================================*/

    //! pointer to general lubrication parameter class
    DRT::ELEMENTS::LubricationEleParameter* lubricationpara_;

    /*========================================================================*/
    //! @name pressure degrees of freedom and related
    /*========================================================================*/

    //! state variables at t_(n+1) or t_(n+alpha_F)
    LINALG::Matrix<nen_,1> eprenp_;

    /*========================================================================*/
    //! @name Galerkin approximation and related
    /*========================================================================*/

    //! coordinates of current integration point in reference coordinates
    LINALG::Matrix<nsd_ele_,1> xsi_;
    //! node coordinates
    LINALG::Matrix<nsd_,nen_> xyze_;
    //! array for shape functions
    LINALG::Matrix<nen_,1> funct_;
    //! array for shape function derivatives w.r.t r,s,t
    LINALG::Matrix<nsd_,nen_> deriv_;
    //! global derivatives of shape functions w.r.t x,y,z
    LINALG::Matrix<nsd_,nen_> derxy_;

    //! transposed jacobian "dx/ds"
    LINALG::Matrix<nsd_,nsd_> xjm_;
    //! inverse of transposed jacobian "ds/dx"
    LINALG::Matrix<nsd_,nsd_> xij_;

    //! bodyforce in element nodes
    LINALG::Matrix<nen_,1> bodyforce_;

    /*========================================================================*/
    //! @name manager classes for efficient application to various problems
    /*========================================================================*/

    //! manager for viscosity
    Teuchos::RCP<LubricationEleViscManager> viscmanager_;

    //! variable manager for Gauss point values
    Teuchos::RCP<LubricationEleInternalVariableManager<nsd_,nen_> > lubricationvarmanager_;

    /*========================================================================*/
    //! @name can be very useful
    /*========================================================================*/

    //! global element id
    int eid_;
    //! current element
    DRT::Element* ele_;


  };

  /// LubricationEleInternalVariableManager implementation
  /*!
    This class manages all internal variables needed for the evaluation of an element.

    All formulation-specific internal variables are stored and managed by a class derived from this class.
  */
  template<int NSD, int NEN>
  class LubricationEleInternalVariableManager
  {
    public:

    LubricationEleInternalVariableManager():
      prenp_(0.0),
      gradpre_(true)
      {return;}

    virtual ~LubricationEleInternalVariableManager(){};

    // compute and set internal variables
    void SetInternalVariables(
        LINALG::Matrix<NEN,1>& funct,  //! array for shape functions
        LINALG::Matrix<NSD,NEN>& derxy, //! global derivatives of shape functions w.r.t x,y,z
        LINALG::Matrix<NEN,1>& eprenp  //! pressure at t_(n+1) or t_(n+alpha_F)
    )
    {
      // calculate pressure at t_(n+1) or t_(n+alpha_F)
      prenp_ = funct.Dot(eprenp);
      // spatial gradient of current pressure value
      gradpre_.Multiply(derxy,eprenp);

      return;
    };

  /*========================================================================*/
  //! @name return methods for internal variables
  /*========================================================================*/

    //! return pressure values at t_(n+1) or t_(n+alpha_F)
    virtual const double& Prenp() const {return prenp_;};
    //! return spatial gradient of all pressure values
    virtual const LINALG::Matrix<NSD,1>& GradPre() const {return gradpre_;};

  protected:

  /*========================================================================*/
  //! @name internal variables evaluated at element center or Gauss point
  /*========================================================================*/

    //! pressure
    double prenp_;
    //! spatial gradient of current pressure value
    LINALG::Matrix<NSD,1> gradpre_;

  };

  /// Lubrication diffusion manager
  /*!
      This is a basic class to handle diffusion. It exclusively contains
      the isotropic diffusion coefficient. For anisotropic diffusion or
      more advanced diffusion laws, e.g., nonlinear ones, a derived class
      has to be constructed in the problem-dependent subclass for element
      evaluation.
  */
  class LubricationEleViscManager
  {
  public:

    LubricationEleViscManager():
      visc_(0.0)
      {return;}

    virtual ~LubricationEleViscManager(){};

    //! Set the isotropic diffusion coefficient
    virtual void SetIsotropicVisc(const double visc)
    {
      if (visc < 0.0) dserror("negative (physical) viscosity: %f",0,visc);

      visc_ = visc;
      return;
    }

    //! Return the stored isotropic diffusion coefficients
    virtual double GetIsotropicDiff(){return visc_;}

  protected:

    //! lubricant viscosity
    double visc_;

  };

}
}

#endif /* SRC_DRT_LUBRICATION_ELE_LUBRICATION_ELE_CALC_H_ */
