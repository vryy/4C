/*!------------------------------------------------------------------------------------------------*
\file topopt_optimizer_ele_impl.H

\brief 

<pre>
Maintainer: Martin Winklmaier
            winklmaier@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>
 *------------------------------------------------------------------------------------------------*/


#ifndef TOPOPT_OPTIMIZER_ELE_IMPL_H_
#define TOPOPT_OPTIMIZER_ELE_IMPL_H_


#include "../drt_fem_general/drt_utils_gausspoints.H"


namespace DRT
{
namespace ELEMENTS
{
// forward declarations
class TopOptParam;


//! Interface base class for ScaTraImpl
/*!
  This class exists to provide a common interface for all template
  versions of ScaTraImpl. The only function
  this class actually defines is Impl, which returns a pointer to
  the appropriate version of ScaTraImpl.
  Deriving from DRT::SingletonDestruction provides proper registration
  of singleton at DRT::Problem instance and - as a result -
  correct memory handling at the end of a BACI run.
  (correctly delete allocated ScaTraImplInterface instances)
 */
class TopOptImplInterface : public DRT::SingletonDestruction
{
public:
  //! Empty constructor
  TopOptImplInterface() {}

  //! Empty destructor
  virtual ~TopOptImplInterface() {}

  //! Evaluate the objective value of element
  /*!
    This class does not provide a definition for this function, it
    must be defined in TopOptImpl.
   */
  virtual int EvaluateObjective(
      DRT::Element*              ele,
      ParameterList&             params,
      DRT::Discretization&       discretization,
      vector<int>&               lm
  ) = 0;

  //! Evaluate the gradient of the objective of element
  /*!
    This class does not provide a definition for this function, it
    must be defined in TopOptImpl.
   */
  virtual int EvaluateGradient(
      DRT::Element*              ele,
      ParameterList&             params,
      DRT::Discretization&       discretization,
      vector<int>&               lm,
      Epetra_SerialDenseVector&  elevec1_epetra
  ) = 0;

  //! return appropriate internal implementation class (acts as a simple factory)
  static TopOptImplInterface* Impl(const DRT::Element* ele);

}; // class TopOptImplInterface



//! Internal Scalar transport element implementation
/*!
  This internal class keeps all the working arrays needed to
  calculate the transport element. Additionally the method Sysmat()
  provides a clean and fast element implementation.

  <h3>Purpose</h3>

  The idea is to separate the element maintenance (class Transport)
  from the mathematical contents (this class). Of course there are
  different implementations of the Transport element, this is just one
  such implementation.

  The Transport element will allocate exactly one object of this class
  for all transport elements with the same number of nodes in the mesh.
  This allows us to use exactly matching working arrays (and keep them
  around.)

  The code is meant to be as clean as possible. This is the only way
  to keep it fast. The number of working arrays has to be reduced to
  a minimum so that the element fits into the cache. (There might be
  room for improvements.)

  <h3>History</h3>

  The implementation here is the standard convection-diffusion element
  capable of dealing with systems of transported scalars.

  Right now we do not read any stabilization parameters from the
  input file but have a fixed version.

  \author gjb
  \date 08/08
 */
template<DRT::Element::DiscretizationType distype>
class TopOptImpl: public TopOptImplInterface
{
  friend class TopOptParam;

public:
  //! Constructor
  TopOptImpl();

  //! Singleton access method
  static TopOptImpl<distype> * Instance(bool create=true );

  //! called upon destruction
  virtual void Done();

  //! Evaluate objective value
  virtual int EvaluateObjective(
      DRT::Element*              ele,            //!< the element we are dealing with
      ParameterList&             params,         //!< the parameter list
      DRT::Discretization&       discretization, //!< the discretization the element is part of
      vector<int>&               lm              //!< location matrix
  );

  //! Evaluate objective value
  virtual int EvaluateObjective(
      DRT::Element*                 ele,            //!< the element we are dealing with
      ParameterList&                params,         //!< the parameter list
      DRT::Discretization&          discretization, //!< the discretization the element is part of
      vector<int>&                  lm,             //!< location matrix
      DRT::UTILS::GaussIntegration& intpoints       ///< gauss points
  );

  //! Evaluate gradient of objective
  virtual int EvaluateGradient(
      DRT::Element*              ele,            //!< the element we are dealing with
      ParameterList&             params,         //!< the parameter list
      DRT::Discretization&       discretization, //!< the discretization the element is part of
      vector<int>&               lm,             //!< location matrix
      Epetra_SerialDenseVector&  elevec1_epetra
  );

  //! Evaluate gradient of objective
  virtual int EvaluateGradient(
      DRT::Element*                 ele,            //!< the element we are dealing with
      ParameterList&                params,         //!< the parameter list
      DRT::Discretization&          discretization, //!< the discretization the element is part of
      vector<int>&                  lm,             //!< location matrix
      Epetra_SerialDenseVector&     elevec1_epetra,
      DRT::UTILS::GaussIntegration& intpoints       ///< gauss points
  );

  /*========================================================================*/
  //! @name static member variables
  /*========================================================================*/

  //! number of element nodes (nomenclature: T. Hughes, The finite element method)
  static const int nen_ = DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;

  //! number of space dimensions
  static const int nsd_ = DRT::UTILS::DisTypeToDim<distype>::dim;

private:

  /*========================================================================*/
  //! @name general framework
  /*========================================================================*/

//  //! Calculate matrix and rhs. Here the whole thing is hidden.
//  virtual void Sysmat(
//      DRT::Element*                         ele,       //!< the element we are dealing with
//      Epetra_SerialDenseMatrix&             emat,      //!< element matrix to calculate
//      Epetra_SerialDenseVector&             erhs,      //!< element rhs to calculate
//      Epetra_SerialDenseVector&             subgrdiff, //!< subgrid-diff.-scaling vector
//      const double                          time,      //!< current simulation time
//      const double                          dt,        //!< current time-step length
//      const double                          timefac,   //!< time discretization factor
//      const double                          alphaF,    //!< factor for generalized-alpha time integration
//      const enum INPAR::SCATRA::AssgdType   whichassgd,//!< all-scale subgrid-diffusivity definition
//      const enum INPAR::SCATRA::FSSUGRDIFF  whichfssgd,//!< fine-scale subgrid-diffusivity definition
//      const bool                            assgd,     //!< all-scale subgrid-diff. flag
//      const bool                            fssgd,     //!< fine-scale subgrid-diff. flag
//      const double                          Cs,        //!< Smagorinsky constant
//      const double                          tpn,       //!< turbulent Prandtl number
//      const double                          Csgs_sgvel,//!< parameter of multifractal subgrid-scales (velocity)
//      const double                          alpha,     //!< grid-filter to test-filter ratio
//      const bool                            calc_N,    //!< flag to activate calculation of N
//      const double                          N_vel,     //!< value for N if not calculated
//      const enum INPAR::FLUID::RefVelocity  refvel,    //!< reference velocity
//      const enum INPAR::FLUID::RefLength    reflength, //!< reference length
//      const double                          c_nu,      //!< scaling for Re
//      const bool                            nwl,       //!< flag to activate near-wall limit
//      const double                          Csgs_sgphi,//!< parameter of multifractal subgrid-scales (phi)
//      const double                          c_diff,    //!< scaling for Re*Pr
//      const bool                            BD_gp,     //!< evaluation of model coefficient at gp
//      const double                          frt,       //!< factor F/RT needed for ELCH calculations
//      const enum INPAR::SCATRA::ScaTraType  scatratype //!< type of scalar transport problem
//  );
//
//  //! get the material parameters
//  virtual void GetMaterialParams(
//      const DRT::Element*                   ele,       //!< the element we are dealing with
//      const enum INPAR::SCATRA::ScaTraType  scatratype //!< type of scalar transport problem
//  );
//
//  //! integral of shape functions over the element
//  void IntegrateShapeFunctions(
//      const DRT::Element*                ele,    //!< the current element
//      Epetra_SerialDenseVector&          elevec1,//!< result vector (to be assembled)
//      const Epetra_IntSerialDenseVector& dofids  //!< for which dof we need to integrate?
//  );
//
//  //! evaluate shape functions and their derivatives at current integration point
//  virtual double EvalShapeFuncAndDerivsAtIntPoint(
//      const DRT::UTILS::IntPointsAndWeights<nsd_>& intpoints, //!< integration points
//      const int                                    iquad,     //!< id of current Gauss point
//      const int                                    eleid      //!< the element id
//  );
//
//  //! calculate error compared to analytical solutions
//  virtual void CalErrorComparedToAnalytSolution(
//      const DRT::Element*          ele,                 //!< the element
//      const enum INPAR::SCATRA::ScaTraType  scatratype, //!< type of scalar transport problem
//      ParameterList&               params,              //!< parameter list
//      Epetra_SerialDenseVector&    errors               //!< computed error contributions
//  );
//
//  //! calculate characteristic element length
//  virtual double CalcCharEleLength(
//      const double  vol,     //!< element volume
//      const double  vel_norm //!?
//  );

  /// extract element data from global vector
  void ExtractValuesFromGlobalVector(
      Teuchos::RCP<DRT::Discretization> discretization, ///< discretization
      const vector<int>&           lm,                  ///<
      LINALG::Matrix<nsd_,nen_> *  matrixtofill,        ///< vector field
      LINALG::Matrix<nen_,1> *     vectortofill,        ///< scalar field
      RCP<Epetra_Vector>&          globalvector         ///< global vector
  ) const;


  //! @name flow data

  /// fluid velocity at gauss point
  std::map<int,LINALG::Matrix<nsd_,1> > fluidvel_;
  /// adjoint velocity at gauss point
  std::map<int,LINALG::Matrix<nsd_,1> > adjointvel_;
  //@}

  //! @name gauss point data

  // gauss points
  DRT::UTILS::GaussIntegration intpoints_;
  //! coordinates of current integration point in reference coordinates
  LINALG::Matrix<nsd_,1> xsi_;
  //! Jacobian determinant
  double det_;
  //! integration factor
  double fac_;
  //@}


  //! @name physical coefficients

  //! physical viscosity
  double visc_;
  //! reaction coefficient
  double reacoeff_;
  //! density
  double dens_;
  //@}

  //! Flag to (de)activate higher order elements
  //! elements with only mixed second order derivatives are not counted as higher order elements
  //! (see definition of higher order elements in fluid3_ele_impl_utils.cpp)
  bool is_higher_order_ele_;

  // general optimization parameters
  Teuchos::RCP<TopOptParam> optiparams_;

}; // end class DRT::ELEMENTS::ScaTraImpl



/// Interface base class for ScaTraImpl
/*!
  This class exists to provide a common interface for all template
  versions of ScaTraImpl. The only function
  this class actually defines is Impl, which returns a pointer to
  the appropriate version of ScaTraImpl.
 */
class TopOptBoundaryImplInterface : public DRT::SingletonDestruction
{

public:
  /// Empty constructor
  TopOptBoundaryImplInterface() {}
  /// Empty destructor
  virtual ~TopOptBoundaryImplInterface() {}

  //! Evaluate the objective value of element
  /*!
    This class does not provide a definition for this function, it
    must be defined in TopOptBoundaryImpl.
   */
  virtual int EvaluateBoundaryObjective(
      DRT::Element*              ele,
      ParameterList&             params,
      DRT::Discretization&       discretization,
      vector<int>&               lm
  ) = 0;

  //! Evaluate the gradient of the objective of element
  /*!
    This class does not provide a definition for this function, it
    must be defined in TopOptBoundaryImpl.
   */
  virtual int EvaluateBoundaryGradient(
      DRT::Element*              ele,
      ParameterList&             params,
      DRT::Discretization&       discretization,
      vector<int>&               lm,
      Epetra_SerialDenseVector&  elevec1_epetra
  ) = 0;

  /// Internal implementation class for scalar transport elements
  static TopOptBoundaryImplInterface* Impl(const DRT::Element* ele);
};



/// Internal Scalar transport element implementation
/*!
    This internal class keeps all the working arrays needed to
    calculate the transport element. Additionally the method Sysmat()
    provides a clean and fast element implementation.

    <h3>Purpose</h3>

    The idea is to separate the element maintenance (class Transport)
    from the mathematical contents (this class). Of course there are
    different implementations of the Transport element, this is just one
    such implementation.

    The Transport element will allocate exactly one object of this class
    for all transport elements with the same number of nodes in the mesh.
    This allows us to use exactly matching working arrays (and keep them
    around.)

    The code is meant to be as clean as possible. This is the only way
    to keep it fast. The number of working arrays has to be reduced to
    a minimum so that the element fits into the cache. (There might be
    room for improvements.)

    <h3>History</h3>

    The implementation here is the standard convection-diffusion element
    capable of dealing with systems of transported scalars.

    Right now we do not read any stabilization parameters from the
    input file but have a fixed version.

    \author gjb
    \date 08/08
 */
template<DRT::Element::DiscretizationType distype>
class TopOptBoundaryImpl: public TopOptBoundaryImplInterface
{
  friend class TopOptParam;

public:
  /// Constructor
  TopOptBoundaryImpl();

  /// Singleton access method
  static TopOptBoundaryImpl<distype> * Instance(bool create = true);

  //! called upon destruction
  virtual void Done();

  //! number of element nodes (nomenclature: T. Hughes, The finite element method)
  static const int nen_ = DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;

  //! number of boundary(!) space dimensions
  static const int nsd_ = DRT::UTILS::DisTypeToDim<distype>::dim;

  //! Evaluate objective value
  virtual int EvaluateBoundaryObjective(
      DRT::Element*              ele,            //!< the element we are dealing with
      ParameterList&             params,         //!< the parameter list
      DRT::Discretization&       discretization, //!< the discretization the element is part of
      vector<int>&               lm              //!< location matrix
  );

  //! Evaluate gradient of objective
  virtual int EvaluateBoundaryGradient(
      DRT::Element*              ele,            //!< the element we are dealing with
      ParameterList&             params,         //!< the parameter list
      DRT::Discretization&       discretization, //!< the discretization the element is part of
      vector<int>&               lm,             //!< location matrix
      Epetra_SerialDenseVector&  elevec1_epetra
  );

private:

//  //! evaluate shape functions and derivatives at int. point
//  double EvalShapeFuncAndIntFac(
//      const DRT::UTILS::IntPointsAndWeights<nsd_>& intpoints,  ///< integration points
//      const int                                    iquad,      ///< id of current Gauss point
//      const int                                    eleid,      ///< the element id
//      LINALG::Matrix<1 + nsd_,1>*  normalvec = NULL  ///< normal vector at Gauss point(optional)
//  );
//
//  //! integral of shape functions over boundary surface
//  void IntegrateShapeFunctions(
//      const DRT::Element*       ele,    ///< the actual boundary element
//      ParameterList&            params, ///< the parameter list
//      Epetra_SerialDenseVector& elevec1,///< result vector (to be assembled)
//      const bool                addarea ///< flag for area calculation
//  );

  /// extract element data from global vector
  void ExtractValuesFromGlobalVector(
      Teuchos::RCP<DRT::Discretization> discretization, ///< discretization
      const vector<int>&           lm,                  ///<
      LINALG::Matrix<nsd_,nen_> *  matrixtofill,        ///< vector field
      LINALG::Matrix<nen_,1> *     vectortofill,        ///< scalar field
      RCP<Epetra_Vector>&          globalvector         ///< global vector
  ) const;


  //! @name gauss point data

  // gauss points
  DRT::UTILS::GaussIntegration intpoints_;
  //! coordinates of current integration point in reference coordinates
  LINALG::Matrix<nsd_,1> xsi_;
  //! Jacobian determinant
  double det_;
  //! integration factor
  double fac_;
  //@}


  //! @name physical coefficients

  //! physical viscosity
  double visc_;
  //! reaction coefficient
  double reacoeff_;
  //! density
  double dens_;
  //@}

  //! Flag to (de)activate higher order elements
  //! elements with only mixed second order derivatives are not counted as higher order elements
  //! (see definition of higher order elements in fluid3_ele_impl_utils.cpp)
  bool is_higher_order_ele_;

  // general optimization parameters
  Teuchos::RCP<TopOptParam> optiparams_;

};

}
}

#endif /* TOPOPT_OPTIMIZER_ELE_IMPL_H_ */
