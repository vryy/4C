/*----------------------------------------------------------------------*/
/*!
 \file porofluidmultiphase_ele_calc.H

 \brief implementation of the evaluation routines of the porofluidmultiphase element

   \level 3

   \maintainer  Anh-Tu Vuong
                vuong@lnm.mw.tum.de
                http://www.lnm.mw.tum.de
                089 - 289-15251
 *----------------------------------------------------------------------*/

#ifndef SRC_DRT_POROFLUIDMULTIPHASE_ELE_FLUIDPOROMULTIPHASE_ELE_CALC_H_
#define SRC_DRT_POROFLUIDMULTIPHASE_ELE_FLUIDPOROMULTIPHASE_ELE_CALC_H_

#include "porofluidmultiphase_ele_interface.H"

#include "porofluidmultiphase_ele_calc_utils.H"
#include "porofluidmultiphase_ele_action.H"

#include "../drt_fem_general/drt_utils_local_connectivity_matrices.H"
#include "../drt_fem_general/drt_utils_integration.H"

namespace MAT
{
  class StructPoro;
}

namespace DRT
{
namespace ELEMENTS
{
// forward declarations
class PoroFluidMultiPhaseEleParameter;

namespace POROFLUIDMANAGER
{
class PhaseManagerInterface;
template<int,int> class VariableManagerInterface;
}

  template<DRT::Element::DiscretizationType distype>
  class PoroFluidMultiPhaseEleCalc : public PoroFluidMultiPhaseEleInterface
  {
  protected:

    /// protected constructor, since we are a Singleton.
    /// this constructor is called from a derived class
    /// -> therefore, it has to be protected instead of private
    PoroFluidMultiPhaseEleCalc(const int numdofpernode, const std::string& disname);

  public:

    /// virtual destructor
    virtual ~PoroFluidMultiPhaseEleCalc() {};

    /// called upon destruction
    virtual void Done();

    //! Singleton access method
    static PoroFluidMultiPhaseEleCalc<distype> * Instance(
        const int numdofpernode,
        const std::string& disname,
        const PoroFluidMultiPhaseEleCalc* delete_me=NULL);

    /*========================================================================*/
    //! @name static member variables
    /*========================================================================*/

    //! number of element nodes (nomenclature: T. Hughes, The finite element method)
    static const int nen_ = DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;

    //! number of space dimensions
    static const int nsd_ = DRT::UTILS::DisTypeToDim<distype>::dim;

    //! number of components necessary to store second derivatives
    // 1 component  for nsd=1:  (N,xx)
    // 3 components for nsd=2:  (N,xx ; N,yy ; N,xy)
    // 6 components for nsd=3:  (N,xx ; N,yy ; N,zz ; N,xy ; N,xz ; N,yz)
    static const int numderiv2_ = DRT::UTILS::DisTypeToNumDeriv2<distype>::numderiv2;

    //! element-type specific flag if second derivatives are needed
    static const bool use2ndderiv_ = POROFLUIDMULTIPHASE::ELEUTILS::Use2ndDerivs<distype>::use;

    /// In this class we do not define a static ScaTraEle...* Instance
    /// since only derived child classes are free to be allocated!!

    /// Setup element evaluation
    virtual int SetupCalc(
        DRT::Element*               ele,
        DRT::Discretization&        discretization,
        const POROFLUIDMULTIPHASE::Action& action
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

//    //! element-type specific flag if second derivatives are needed
//    static const bool use2ndderiv_ = POROFLUIDMULTIPHASE::Use2ndDerivs<distype>::use;
//
//    //! number of components necessary to store second derivatives
//    // 1 component  for nsd=1:  (N,xx)
//    // 3 components for nsd=2:  (N,xx ; N,yy ; N,xy)
//    // 6 components for nsd=3:  (N,xx ; N,yy ; N,zz ; N,xy ; N,xz ; N,yz)
//    static const int numderiv2_ = DRT::UTILS::DisTypeToNumDeriv2<distype>::numderiv2;

  protected:
    /*========================================================================*/
    //! @name general framework
    /*========================================================================*/

    //! evaluate action
    virtual int EvaluateAction(
        DRT::Element*                           ele,
        Teuchos::ParameterList&                 params,
        DRT::Discretization&                    discretization,
        const POROFLUIDMULTIPHASE::Action&      action,
        DRT::Element::LocationArray&            la,
        Epetra_SerialDenseMatrix&               elemat1_epetra,
        Epetra_SerialDenseMatrix&               elemat2_epetra,
        Epetra_SerialDenseVector&               elevec1_epetra,
        Epetra_SerialDenseVector&               elevec2_epetra,
        Epetra_SerialDenseVector&               elevec3_epetra
        );

    //! extract element based or nodal values
    //  return extracted values of phinp
    virtual void ExtractElementAndNodeValues(
      DRT::Element*                 ele,
      Teuchos::ParameterList&       params,
      DRT::Discretization&          discretization,
      DRT::Element::LocationArray&  la
    );

    /// access structure material of corresponding (poro) element
    void GetStructMaterial(DRT::Element*                 ele);

    /// Setup element evaluation
    virtual void PrepareGaussPointLoop(
        DRT::Element*                         ele        ///< the element whose matrix is calculated
        );

    virtual void GaussPointLoop(
        const DRT::UTILS::IntPointsAndWeights<nsd_>& intpoints,  ///< integration points
        DRT::Element*                         ele,       //!< the element we are dealing with
        Epetra_SerialDenseMatrix&             emat,      //!< element matrix to calculate
        Epetra_SerialDenseVector&             erhs       //!< element rhs to calculate
        );

    //! calculate matrix and rhs. Here the whole thing is hidden.
    virtual void Sysmat(
      DRT::Element*                         ele,       //!< the element we are dealing with
      Epetra_SerialDenseMatrix&             emat,      //!< element matrix to calculate
      Epetra_SerialDenseVector&             erhs       //!< element rhs to calculate
    );

    //! evaluate shape functions and their derivatives at current integration point
    double EvalShapeFuncAndDerivsAtIntPoint(
      const DRT::UTILS::IntPointsAndWeights<nsd_>& intpoints, //!< integration points
      const int                                    iquad      //!< id of current Gauss point
    );

    //! evaluate shape functions and their derivatives at current integration point
    double EvalShapeFuncAndDerivsInParameterSpace();

    //! calculate divergence of vector field (e.g., velocity)
    void GetDivergence(
      double&                          vdiv,
      const LINALG::Matrix<nsd_,nen_>& evel
    );

//    //! compute rhs containing bodyforce
//    virtual void GetRhsInt(
//      double&      rhsint,  //!< rhs containing bodyforce at integration point
//      const double densnp,  //!< density at t_(n+1)
//      const int    k        //!< index of current scalar
//    );

    //! calculation of convective element matrix in convective form
    void CalcMatConv(
      Epetra_SerialDenseMatrix&     emat,        //!< element matrix to be filled
      const int                     k,           //!< index of current scalar
      const POROFLUIDMANAGER::PhaseManagerInterface&  phasemanager,//!< phase manager
      const double                  porosity,    //!< porosity at gauss point
      const double                  timefacfac,  //!< domain-integration factor times time-integration factor
      const LINALG::Matrix<nen_,1>& conv         //!< convective term
    );


    //! calculation of convective element matrix: add conservative contributions
    virtual void CalcMatConvCons(
      Epetra_SerialDenseMatrix&     emat,           //!< element matrix to be filled
      const int                     k,              //!< index of current scalar
      const POROFLUIDMANAGER::PhaseManagerInterface&  phasemanager,   //!< phase manager
      const double                  timefacfac,     //!< domain-integration factor times time-integration factor
      const double                  vdiv            //!< velocity divergence
    );

    //! calculation of diffusive element matrix
    virtual void CalcMatDiff(
      Epetra_SerialDenseMatrix&         emat,         //!< element matrix to be filled
      const int                         curphase,     //!< index of current phase
      const int                         phasetoadd,   //!< index of current phase
      const POROFLUIDMANAGER::PhaseManagerInterface&      phasemanager, //!< phase manager
      const MAT::Material&              material,     //!< material of the element
      const double                      timefacfac,   //!< domain-integration factor times time-integration factor
      const LINALG::Matrix<nsd_,nsd_>&  difftensor    //!< diffusion tensor
    );

    //! calculation of reactive element matrix
    virtual void CalcMatReac(
      Epetra_SerialDenseMatrix&          emat,               //!< element matrix to be filled
      const int                          curphase,           //!< index of current phase
      const int                          phasetoadd,         //!< index of current phase
      const POROFLUIDMANAGER::PhaseManagerInterface&       phasemanager,       //!< phase manager
      const MAT::Material&               material,           //!< material of the element
      const double                       timefacfac,         //!< domain-integration factor times time-integration factor
      bool                               scalewithsaturation //!< flag wether to scale with saturation
    );

    //! calculation of mass element matrix
    void CalcMatMassPressure(
      Epetra_SerialDenseMatrix&     emat,          //!< element matrix to be filled
      const int                     curphase,      //!< index of current phase
      const int                     phasetoadd,    //!< index of current phase
      const POROFLUIDMANAGER::PhaseManagerInterface&  phasemanager,  //!< phase manager
      const MAT::Material&          material,      //!< material of the element
      const double                  fac,           //!< domain-integration factor
      const double                  timefacfac,    //!< domain-integration factor times time-integration factor
      const double                  porosity,      //!< porosity at gauss point
      const double                  hist,          //!< history of time integartion
      const std::vector<double>&    phinp,         //!< state phi at t_{n+1}
      const std::vector<double>&    phidtnp        //!< time derivative state phi at t_{n+1}
    );

    //! calculation of mass element matrix
    void CalcMatMassSolidPressure(
      Epetra_SerialDenseMatrix&     emat,                //!< element matrix to be filled
      const int                     curphase,            //!< index of current scalar
      const POROFLUIDMANAGER::PhaseManagerInterface&  phasemanager,        //!< phase manager
      const MAT::Material&          material,            //!< material of the element
      const double                  fac,                 //!< domain-integration factor
      const double                  timefacfac,          //!< domain-integration factor times time-integration factor
      const double                  porosity,            //!< porosity at gauss point
      const double                  hist,                //!< history of time integartion
      const std::vector<double>&    phinp,               //!< state phi at t_{n+1}
      const std::vector<double>&    phidtnp,             //!< time derivative state phi at t_{n+1}
      bool                          scalewithsaturation  //!< flag for saturation scaling
    );

    //! calculation of mass element matrix
    void CalcMatMassSaturation(
      Epetra_SerialDenseMatrix&     emat,         //!< element matrix to be filled
      const int                     curphase,     //!< index of current scalar
      const POROFLUIDMANAGER::PhaseManagerInterface&  phasemanager, //!< phase manager
      const MAT::Material&          material,     //!< material of the element
      const double                  fac,          //!< domain-integration factor
      const double                  porosity,     //!< porosity at gauss point
      const std::vector<double>&    phinp         //!< state phi at t_{n+1}
    );

    //! calculation of linearized mass rhs vector
    void CalcRHSMassPressure(
      Epetra_SerialDenseVector&     erhs,         //!< element vector to be filled
      const int                     curphase,     //!< index of current phase
      const int                     phasetoadd,   //!< index of current phase
      const POROFLUIDMANAGER::PhaseManagerInterface&  phasemanager, //!< phase manager
      const MAT::Material&          material,     //!< material of the element
      const double                  rhsfac,       //!< time-integration factor for rhs times domain-integration factor
      const double                  fac,          //!< domain-integration factor
      const double                  hist,         //!< history of time integartion
      const double                  porosity,     //!< porosity at gauss point
      const std::vector<double>&    phinp,        //!< state phi at t_{n+1}
      const std::vector<double>&    phidtnp       //!< time derivative state phi at t_{n+1}
    );

    //! calculation of linearized mass rhs vector
    void CalcRHSMassSolidPressure(
      Epetra_SerialDenseVector&     erhs,                 //!< element vector to be filled
      const int                     curphase,             //!< index of current phase
      const int                     phasetoadd,           //!< index of current phase
      const POROFLUIDMANAGER::PhaseManagerInterface&  phasemanager,         //!< phase manager
      const MAT::Material&          material,             //!< material of the element
      const double                  rhsfac,               //!< time-integration factor for rhs times domain-integration factor
      const double                  fac,                  //!< domain-integration factor
      const double                  hist,                 //!< history of time integartion
      const double                  porosity,             //!< porosity at gauss point
      const std::vector<double>&    phinp,                //!< state phi at t_{n+1}
      const std::vector<double>&    phidtnp,              //!< time derivative state phi at t_{n+1}
      bool                          scalewithsaturation   //!< flag for saturation scaling
    );

    //! calculation of linearized mass rhs vector
    void CalcRHSMassSaturation(
      Epetra_SerialDenseVector&     erhs,         //!< element vector to be filled
      const int                     curphase,     //!< index of current phase
      const int                     phasetoadd,   //!< index of current phase
      const POROFLUIDMANAGER::PhaseManagerInterface&  phasemanager, //!< phase manager
      const MAT::Material&          material,     //!< material of the element
      const double                  rhsfac,       //!< time-integration factor for rhs times domain-integration factor
      const double                  fac,          //!< domain-integration factor
      const double                  hist,         //!< history of time integartion
      const double                  porosity,     //!< porosity at gauss point
      const std::vector<double>&    phinp,        //!< state phi at t_{n+1}
      const std::vector<double>&    phidtnp       //!< time derivative state phi at t_{n+1}
    );

    //! standard Galerkin transient, old part of rhs and source term
    void CalcRHSHistAndSource(
      Epetra_SerialDenseVector&     erhs,        //!< element vector to be filled
      const int                     curphase,    //!< index of current phase
      const double                  fac,         //!< domain-integration factor
      const double                  rhsint       //!< rhs at Gauss point
    );

    //! standard Galerkin convective term on right hand side
    void CalcRHSConv(
      Epetra_SerialDenseVector&     erhs,         //!< element vector to be filled
      const int                     curphase,     //!< index of current phase
      const POROFLUIDMANAGER::PhaseManagerInterface&  phasemanager, //!< phase manager
      const double                  porosity,     //!< porosity at gauss point
      const double                  rhsfac,       //!< time-integration factor for rhs times domain-integration factor
      const std::vector<double>&    conv_phi      //!< convective term
    );

    //! standard Galerkin convective term (conservative contribution) on right hand side
    void CalcRHSConvCons(
      Epetra_SerialDenseVector&     erhs,                 //!< element vector to be filled
      const int                     curphase,             //!< index of current phase
      const POROFLUIDMANAGER::PhaseManagerInterface&  phasemanager,         //!< phase manager
      const double                  rhsfac,               //!< time-integration factor for rhs times domain-integration factor
      const double                  vdiv,                 //!< convective term
      bool                          scalewithsaturation   //!< flag wether to scale with saturation
    );

    //! standard Galerkin diffusive term on right hand side
    virtual void CalcRHSDiff(
      Epetra_SerialDenseVector&                   erhs,         //!< element vector to be filled
      const int                                   curphase,     //!< index of current phase
      const int                                   phasetoadd,   //!< index of current phase
      const POROFLUIDMANAGER::PhaseManagerInterface&                phasemanager, //!< phase manager
      const double                                rhsfac,       //!< time-integration factor for rhs times domain-integration factor
      const std::vector<LINALG::Matrix<nsd_,1> >& gradphi,      //!< gradient of scalar
      const LINALG::Matrix<nsd_,nsd_>&            difftensor    //!< pre factor of diffusive term
    );

    //! reaction terms on right hand side
    virtual void CalcRHSReac(
      Epetra_SerialDenseVector&             erhs,               //!< element vector to be filled
      const int                             curphase,           //!< index of current phase
      const int                             phasetoadd,         //!< index of current phase
      const POROFLUIDMANAGER::PhaseManagerInterface&          phasemanager,       //!< phase manager
      const MAT::Material&                  material,           //!< material of the element
      const double                          rhsfac,             //!< time-integration factor for rhs times domain-integration factor
      bool                                  scalewithsaturation //!< flag wether to scale with saturation
    );

    void ComputeDiffTensor(
        const POROFLUIDMANAGER::PhaseManagerInterface&  phasemanager, //!< phase manager
        const int                     phase,        //!< index of current phase
        const MAT::Material&          material,     //!< material of the element
        const std::vector<double>&    phinp,        //!< state phi at t_{n+1}
        LINALG::Matrix<nsd_,nsd_>&    difftensor    //!< diffusion tensor
        ) const;

    /*========================================================================*/
    //! @name routines for additional element evaluations (called from EvaluateAction)
    /*========================================================================*/

    //! calculate mass matrix and rhs for determining initial time derivative
    virtual void CalcPressureAndSaturation(
      DRT::Element*                 ele,              //!< current element
      Epetra_SerialDenseVector&     pressure,         //!< pressure vector
      Epetra_SerialDenseVector&     saturation,       //!< saturation vector
      Epetra_SerialDenseVector&     counter,          //!< counter
      Teuchos::ParameterList&       params,           //!< parameter list
      DRT::Discretization&          discretization,   //!< discretization
      DRT::Element::LocationArray&  la                //!< location array
    );

    virtual void CalcSolidPressure(
      DRT::Element*                 ele,              //!< current element
      Epetra_SerialDenseVector&     solidpressure,    //!< pressure vector
      Epetra_SerialDenseVector&     counter,          //!< counter
      Teuchos::ParameterList&       params,           //!< parameter list
      DRT::Discretization&          discretization,   //!< discretization
      DRT::Element::LocationArray&  la                //!< location array
    );

    virtual void ReconFlux(
      DRT::Element*                 ele,              //!< current element
      Epetra_SerialDenseMatrix&     elemat1,          //!< linearizations
      Epetra_SerialDenseMatrix&     elemat2,          //!< rhs
      Teuchos::ParameterList&       params,           //!< parameter list
      DRT::Discretization&          discretization,   //!< discretization
      DRT::Element::LocationArray&  la                //!< location array
    );

    /*========================================================================*/
    //! @name material and related functions
    /*========================================================================*/

//    //! get the material parameters
//    virtual void GetMaterialParams(
//      const DRT::Element* ele,       //!< the element we are dealing with
//      double&             densn,     //!< density at t_(n)
//      double&             densnp,    //!< density at t_(n+1) or t_(n+alpha_F)
//      double&             densam,    //!< density at t_(n+alpha_M)
//      double&             visc,      //!< fluid viscosity
//      const int           iquad      //!< id of current gauss point (default = -1)
//    );

  private:

    /*========================================================================*/
    //! @name can be very useful
    /*========================================================================*/

    //! element
    DRT::Element* ele_;

    /*========================================================================*/
    //! @name dofs and nodes
    /*========================================================================*/

    //! number of dof per node (= number of phases)
    const int numdofpernode_;

    /*========================================================================*/
    //! @name parameter lists
    /*========================================================================*/

    //! pointer to general scalar transport parameter class
    DRT::ELEMENTS::PoroFluidMultiPhaseEleParameter* para_;

    /*========================================================================*/
    //! @name Galerkin approximation and related
    /*========================================================================*/

    //! coordinates of current integration point in reference coordinates
    LINALG::Matrix<nsd_,1> xsi_;
    //! initial node coordinates
    LINALG::Matrix<nsd_,nen_> xyze0_;
    //! current node coordinates
    LINALG::Matrix<nsd_,nen_> xyze_;
    //! array for shape functions
    LINALG::Matrix<nen_,1> funct_;
    //! array for shape function derivatives w.r.t r,s,t
    LINALG::Matrix<nsd_,nen_> deriv_;
    //! array for second derivatives of shape function w.r.t r,s,t
    LINALG::Matrix<numderiv2_,nen_> deriv2_;
    //! global derivatives of shape functions w.r.t x,y,z
    LINALG::Matrix<nsd_,nen_> derxy_;
    //! global second derivatives of shape functions w.r.t x,y,z
    LINALG::Matrix<numderiv2_,nen_> derxy2_;

    //! transposed jacobian "dx/ds"
    LINALG::Matrix<nsd_,nsd_> xjm_;
    //! inverse of transposed jacobian "ds/dx"
    LINALG::Matrix<nsd_,nsd_> xij_;
    //! determinant of jacobian "dx/ds"
    double det_;
    //! determinant of deformation gradient "dx/dX"
    double J_;

    /*========================================================================*/
    //! @name scalar degrees of freedom and related
    /*========================================================================*/
    //! manager class for variables
    Teuchos::RCP< POROFLUIDMANAGER::VariableManagerInterface<nsd_,nen_> > variablemanager_;

    //! corresponding poro structure material
    Teuchos::RCP< MAT::StructPoro > structmat_;

    //! manager class for handling phases and corresponding DOFs
    Teuchos::RCP< POROFLUIDMANAGER::PhaseManagerInterface > phasemanager_;

  };

}
}

#endif /* SRC_DRT_POROFLUIDMULTIPHASE_ELE_FLUIDPOROMULTIPHASE_ELE_CALC_H_ */
