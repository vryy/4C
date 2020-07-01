/*---------------------------------------------------------------------*/
/*! \file

\brief element routines of the topology optimization element

\level 3


*/
/*---------------------------------------------------------------------*/

#ifndef TOPOPT_OPTIMIZER_ELE_IMPL_H_
#define TOPOPT_OPTIMIZER_ELE_IMPL_H_


#include "../drt_fem_general/drt_utils_gausspoints.H"
#include "../drt_fem_general/drt_utils_integration.H"


namespace DRT
{
  namespace ELEMENTS
  {
    // forward declarations
    class TopOptParam;

    /*
    \brief Interface base class for TopOptImpl

      This class exists to provide a common interface for all template
      versions of TopOptImpl. The only function
      this class actually defines is Impl, which returns a pointer to
      the appropriate version of ScaTraImpl.
      Deriving from DRT::SingletonDestruction provides proper registration
      of singleton at DRT::Problem instance and - as a result -
      correct memory handling at the end of a BACI run.
      (correctly delete allocated TopOptImplInterface instances)

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
      virtual int EvaluateValues(DRT::Element* ele, Teuchos::ParameterList& params,
          DRT::Discretization& optidis, Teuchos::RCP<MAT::Material> mat) = 0;

      //! Evaluate the gradient of the objective of element
      /*!
        This class does not provide a definition for this function, it
        must be defined in TopOptImpl.
       */
      virtual int EvaluateGradients(DRT::Element* ele, Teuchos::ParameterList& params,
          DRT::Discretization& optidis, Teuchos::RCP<MAT::Material> mat) = 0;

      //! return appropriate internal implementation class (acts as a simple factory)
      static TopOptImplInterface* Impl(const DRT::Element* ele);

    };  // class TopOptImplInterface



    /*
    \brief Internal topology optimization element implementation

      This internal class keeps all the working arrays needed to
      calculate the topopt element.

      <h3>Purpose</h3>

      The idea is to separate the element maintenance (class TopOpt)
      from the mathematical contents (this class).

      The topopt element will allocate exactly one object of this class
      for all topopt elements with the same number of nodes in the mesh.
      This allows us to use exactly matching working arrays (and keep them
      around.)

      The code is meant to be as clean as possible. This is the only way
      to keep it fast. The number of working arrays has to be reduced to
      a minimum so that the element fits into the cache. (There might be
      room for improvements.)

      <h3>History</h3>

      \author winklmaier
     */
    template <DRT::Element::DiscretizationType distype>
    class TopOptImpl : public TopOptImplInterface
    {
      friend class TopOptParam;

     public:
      //! Constructor
      TopOptImpl();

      //! Singleton access method
      static TopOptImpl<distype>* Instance(bool create = true);

      //! called upon destruction
      virtual void Done();

      //! Evaluate objective value
      virtual int EvaluateValues(DRT::Element* ele,  //!< the element we are dealing with
          Teuchos::ParameterList& params,            //!< the parameter list
          DRT::Discretization& optidis,              //!< the discretization the element is part of
          Teuchos::RCP<MAT::Material> mat);

      //! Evaluate objective value
      virtual int EvaluateValues(DRT::Element* ele,  //!< the element we are dealing with
          Teuchos::ParameterList& params,            //!< the parameter list
          DRT::Discretization& optidis,              //!< the discretization the element is part of
          Teuchos::RCP<MAT::Material> mat,
          DRT::UTILS::GaussIntegration& intpoints  ///< gauss points
      );

      //! Evaluate gradient of objective
      virtual int EvaluateGradients(DRT::Element* ele,  //!< the element we are dealing with
          Teuchos::ParameterList& params,               //!< the parameter list
          DRT::Discretization& optidis,  //!< the discretization the element is part of
          Teuchos::RCP<MAT::Material> mat);

      //! Evaluate gradient of objective
      virtual int EvaluateGradients(DRT::Element* ele,  //!< the element we are dealing with
          Teuchos::ParameterList& params,               //!< the parameter list
          DRT::Discretization& optidis,  //!< the discretization the element is part of
          Teuchos::RCP<MAT::Material> mat,
          DRT::UTILS::GaussIntegration& intpoints  ///< gauss points
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

      /// calculate objective value - main computation function
      virtual void Values(const int eid, std::map<int, LINALG::Matrix<nsd_, nen_>>& efluidvel,
          LINALG::Matrix<nen_, 1>& edens, double& objective, Epetra_SerialDenseVector& constraints,
          Teuchos::RCP<MAT::Material> mat, DRT::UTILS::GaussIntegration& intpoints);

      /// Calculate gradient of objective - main computation function
      virtual void Gradients(const int eid, std::map<int, LINALG::Matrix<nsd_, nen_>>& efluidvels,
          std::map<int, LINALG::Matrix<nsd_, nen_>>& eadjointvels,
          std::map<int, LINALG::Matrix<nen_, 1>>& eadjointpress, LINALG::Matrix<nen_, 1>& edens,
          Epetra_SerialDenseVector& obj_deriv, Epetra_SerialDenseVector& constr_deriv,
          Teuchos::RCP<MAT::Material> mat, DRT::UTILS::GaussIntegration& intpoints);

      //! evaluate shape functions and their derivatives at current integration point
      virtual void EvalShapeFuncAndDerivsAtIntPoint(
          DRT::UTILS::GaussIntegration::iterator& iquad,  ///< actual integration point
          const int eleid                                 //!< the element id
      );

      /// evaluate porosity at integration point
      void EvalPorosityAtIntPoint(LINALG::Matrix<nen_, 1> edens  /// nodal densities
      );

      /// calculate stabilization parameter
      void CalcStabParameter(const double vol  /// volume
      );

      //! calculate characteristic element length
      void CalcCharEleLength(const double vol,  ///< volume
          const double fluidvel_norm,           ///< norm of velocity vector
          double& strle,                        ///< streamlength (for tau_Mu)
          double& hk) const;                    ///< length for tau_Mp

      /// extract element data from global vector
      void ExtractValuesFromGlobalVector(DRT::Discretization& optidis,  ///< discretization
          const std::vector<int>& lm,                                   ///<
          LINALG::Matrix<nsd_, nen_>* matrixtofill,                     ///< vector field
          LINALG::Matrix<nen_, 1>* vectortofill,                        ///< scalar field
          Teuchos::RCP<Epetra_Vector>& globalvector                     ///< global vector
          ) const;


      //! @name coordinate data

      //! node coordinates
      LINALG::Matrix<nsd_, nen_> xyze_;
      //@}


      //! @name shape functions

      //! array for shape functions
      LINALG::Matrix<nen_, 1> funct_;
      //! array for shape function derivatives w.r.t r,s,t
      LINALG::Matrix<nsd_, nen_> deriv_;
      //! global derivative of shape functions w.r.t x,y,z
      LINALG::Matrix<nsd_, nen_> derxy_;
      //@}


      //! @name flow data

      /// adjoint pressure of one time step
      LINALG::Matrix<nen_, 1> eadjointpres_;

      /// fluid velocity of one time step at gauss point
      LINALG::Matrix<nsd_, 1> fluidvelint_;
      /// adjoint velocity of one time step at gauss point
      LINALG::Matrix<nsd_, 1> adjointvelint_;

      /// first derivative of fluid velocity of one time step at gauss point w.r.t. x,y,z
      LINALG::Matrix<nsd_, nsd_> fluidvelxy_;
      /// first derivative of adjoint velocity of one time step at gauss point w.r.t. x,y,z
      LINALG::Matrix<nsd_, nsd_> adjointvelxy_;
      /// first derivative of adjoint pressure of one time step at gauss point w.r.t. x,y,z
      LINALG::Matrix<nsd_, 1> adjointpresxy_;
      /// convective adjoint velocity of one time step
      LINALG::Matrix<nsd_, 1> adjointconvvel_;

      /// fluid velocity of last time step at gauss point
      LINALG::Matrix<nsd_, 1> fluidvelint_old_;

      //! the stabilisation parameters -> it is a (3,1) vector for 2D and 3D
      LINALG::Matrix<3, 1> tau_;  //@}

      /// @name optimization data

      /// porosity at integration point
      double poroint_;
      /// first derivative of porosity at gauss point w.r.t. density function
      double poroderdens_;
      //@}

      //! @name gauss point data

      // gauss points
      DRT::UTILS::GaussIntegration intpoints_;
      //! coordinates of current integration point in reference coordinates
      LINALG::Matrix<nsd_, 1> xsi_;
      /// jacobian dx/ds
      LINALG::Matrix<nsd_, nsd_> xjm_;
      /// invers transposed jacobian ds/dx
      LINALG::Matrix<nsd_, nsd_> xji_;

      //! integration factor
      double fac_;
      //@}


      //! Flag to (de)activate higher order elements
      //! elements with only mixed second order derivatives are not counted as higher order elements
      //! (see definition of higher order elements in fluid3_ele_impl_utils.cpp)
      bool is_higher_order_ele_;

      // general optimization parameters
      Teuchos::RCP<TopOptParam> optiparams_;

    };  // end class DRT::ELEMENTS::ScaTraImpl



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
      virtual int EvaluateBoundaryValues(DRT::Element* ele, Teuchos::ParameterList& params,
          DRT::Discretization& optidis, Teuchos::RCP<MAT::Material> mat, std::vector<int>& lm) = 0;

      //! Evaluate the gradient of the objective of element
      /*!
        This class does not provide a definition for this function, it
        must be defined in TopOptBoundaryImpl.
       */
      virtual int EvaluateBoundaryGradients(DRT::Element* ele, Teuchos::ParameterList& params,
          DRT::Discretization& optidis, Teuchos::RCP<MAT::Material> mat, std::vector<int>& lm,
          Epetra_SerialDenseVector& elevec1_epetra) = 0;

      /// Internal implementation class for scalar transport elements
      static TopOptBoundaryImplInterface* Impl(const DRT::Element* ele);
    };



    /*
    \brief Internal topology optimization element implementation

      This internal class keeps all the working arrays needed to
      calculate the boundary terms of the topopt element.

      <h3>Purpose</h3>

      The idea is to separate the element maintenance (class TopOpt)
      from the mathematical contents (this class).

      The topopt element will allocate exactly one object of this class
      for all topopt elements with the same number of nodes in the mesh.
      This allows us to use exactly matching working arrays (and keep them
      around.)

      The code is meant to be as clean as possible. This is the only way
      to keep it fast. The number of working arrays has to be reduced to
      a minimum so that the element fits into the cache. (There might be
      room for improvements.)

      <h3>History</h3>

      \author winklmaier
     */
    template <DRT::Element::DiscretizationType distype>
    class TopOptBoundaryImpl : public TopOptBoundaryImplInterface
    {
      friend class TopOptParam;

     public:
      /// Constructor
      TopOptBoundaryImpl();

      /// Singleton access method
      static TopOptBoundaryImpl<distype>* Instance(bool create = true);

      //! called upon destruction
      virtual void Done();

      //! number of element nodes (nomenclature: T. Hughes, The finite element method)
      static const int nen_ = DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;

      //! number of boundary(!) space dimensions
      static const int nsd_ = DRT::UTILS::DisTypeToDim<distype>::dim;

      //! Evaluate objective value
      virtual int EvaluateBoundaryValues(DRT::Element* ele,  //!< the element we are dealing with
          Teuchos::ParameterList& params,                    //!< the parameter list
          DRT::Discretization& optidis,  //!< the discretization the element is part of
          Teuchos::RCP<MAT::Material> mat,
          std::vector<int>& lm  //!< location matrix
      );

      //! Evaluate gradient of objective
      virtual int EvaluateBoundaryGradients(DRT::Element* ele,  //!< the element we are dealing with
          Teuchos::ParameterList& params,                       //!< the parameter list
          DRT::Discretization& optidis,  //!< the discretization the element is part of
          Teuchos::RCP<MAT::Material> mat,
          std::vector<int>& lm,  //!< location matrix
          Epetra_SerialDenseVector& elevec1_epetra);

     private:
      //  //! evaluate shape functions and derivatives at int. point
      //  double EvalShapeFuncAndIntFac(
      //      const DRT::UTILS::IntPointsAndWeights<nsd_>& intpoints,  ///< integration points
      //      const int                                    iquad,      ///< id of current Gauss
      //      point const int                                    eleid,      ///< the element id
      //      LINALG::Matrix<1 + nsd_,1>*  normalvec = NULL  ///< normal vector at Gauss
      //      point(optional)
      //  );
      //
      //  //! integral of shape functions over boundary surface
      //  void IntegrateShapeFunctions(
      //      const DRT::Element*       ele,    ///< the actual boundary element
      //      Teuchos::ParameterList&   params, ///< the parameter list
      //      Epetra_SerialDenseVector& elevec1,///< result vector (to be assembled)
      //      const bool                addarea ///< flag for area calculation
      //  );

      /// extract element data from global vector
      void ExtractValuesFromGlobalVector(DRT::Discretization& optidis,  ///< discretization
          const std::vector<int>& lm,                                   ///<
          LINALG::Matrix<nsd_, nen_>* matrixtofill,                     ///< vector field
          LINALG::Matrix<nen_, 1>* vectortofill,                        ///< scalar field
          Teuchos::RCP<Epetra_Vector>& globalvector                     ///< global vector
          ) const;


      //! @name gauss point data

      // gauss points
      DRT::UTILS::GaussIntegration intpoints_;
      //! coordinates of current integration point in reference coordinates
      LINALG::Matrix<nsd_, 1> xsi_;
      //! integration factor
      double fac_;
      //@}

      // general optimization parameters
      Teuchos::RCP<TopOptParam> optiparams_;
    };

  }  // namespace ELEMENTS
}  // namespace DRT

#endif /* TOPOPT_OPTIMIZER_ELE_IMPL_H_ */
