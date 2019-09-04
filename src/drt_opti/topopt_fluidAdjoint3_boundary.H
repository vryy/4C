/*---------------------------------------------------------------------*/
/*! \file

\brief boundary element implementation of fluid adjoint equations for topology optimization

\maintainer Martin Kronbichler

\level 3

*/
/*---------------------------------------------------------------------*/

#ifndef TOPOPT_FLUIDADJOINT3_BOUNDARY_H_
#define TOPOPT_FLUIDADJOINT3_BOUNDARY_H_


#include "../drt_fem_general/drt_utils_gausspoints.H"


namespace DRT
{
  class Discretization;
  namespace ELEMENTS
  {
    class FluidBoundary;
    class FluidAdjoint3ImplParameter;



    class FluidAdjoint3BoundaryImplInterface : public DRT::SingletonDestruction
    {
     public:
      /// Empty constructor
      FluidAdjoint3BoundaryImplInterface() {}

      /// Evaluate a Neumann boundary condition
      /*!
        This class does not provide a definition for this function, it
        must be defined in FluidBoundaryImpl.
       */
      virtual int EvaluateNeumann(DRT::ELEMENTS::FluidBoundary* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, const std::vector<int>& lm,
          Epetra_SerialDenseVector& elevec) = 0;

      /// Internal implementation class for FluidAdjointBoundary elements
      static FluidAdjoint3BoundaryImplInterface* Impl(const DRT::Element* ele);
    };

    /// Internal FluidAdjoint3Boundary element implementation
    /*!
      This internal class keeps all the working arrays needed to
      calculate the FluidAdjoint3Boundary element.

      The FluidAdjoint3Boundary element will allocate exactly one object of this class
      for all FluidAdjoint3Boundary elements with the same number of nodes in the mesh.
      This allows us to use exactly matching working arrays (and keep them
      around.)

      \author winklmaier
      \date 03/12
    */
    template <DRT::Element::DiscretizationType distype>
    class FluidAdjoint3BoundaryImpl : public FluidAdjoint3BoundaryImplInterface
    {
      friend class FluidEleParameter;

     public:
      /// Singleton access method
      static FluidAdjoint3BoundaryImpl<distype>* Instance(bool create = true);

      /// called upon destruction
      virtual void Done();

      /// Constructor
      FluidAdjoint3BoundaryImpl();

      //! number of element nodes
      static const int nen_ = DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;

      //! number of space dimensions of the FluidBoundary element
      static const int bdrynsd_ = DRT::UTILS::DisTypeToDim<distype>::dim;

      //! number of space dimensions of the parent element
      static const int nsd_ = bdrynsd_ + 1;

      //! number of dof's per node
      static const int numdofpernode_ = nsd_ + 1;

      //! Evaluate a Neumann boundary condition
      virtual int EvaluateNeumann(DRT::ELEMENTS::FluidBoundary* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, const std::vector<int>& lm,
          Epetra_SerialDenseVector& elevec);

     private:
      //! evaluate shape functions/derivatives, metric tensor and
      //! optionally, the infinitesimal area element and the normalized normal vector the at
      //! integration point
      void EvalShapeFuncAtBouIntPoint(
          DRT::UTILS::GaussIntegration::iterator& iquad,  // actual integration point
          const int eleid                                 // element ID
      );

      //! get material parameters
      void GetMaterialParams(
          Teuchos::RCP<const MAT::Material> material  ///< reference pointer to material
      );

      /*!
       * \brief fill elment matrix and vectors with the global values
       */
      void ExtractValuesFromGlobalVector(
          const DRT::Discretization& discretization,  ///< discretization
          const std::vector<int>& lm,                 ///<
          LINALG::Matrix<nsd_, nen_>* matrixtofill,   ///< vector field
          LINALG::Matrix<nen_, 1>* vectortofill,      ///< scalar field
          const std::string state                     ///< state of the global vector
          ) const;

      //! pointer to parameter list
      Teuchos::RCP<DRT::ELEMENTS::FluidAdjoint3ImplParameter> fldAdPara_;

      // gauss points
      DRT::UTILS::GaussIntegration intpoints_;

      //! node coordinates for boundary element
      LINALG::Matrix<nsd_, nen_> xyze_;
      //! coordinates of current integration point in reference coordinates
      LINALG::Matrix<bdrynsd_, 1> xsi_;
      //! array for shape functions for boundary element
      LINALG::Matrix<nen_, 1> funct_;
      //! array for shape function derivatives for boundary element
      LINALG::Matrix<bdrynsd_, nen_> deriv_;

      //! normal vector pointing out of the domain
      LINALG::Matrix<nsd_, 1> unitnormal_;

      //! velocity vector at integration point
      LINALG::Matrix<nsd_, 1> velint_;
      //! velocity vector at integration point at old time step
      LINALG::Matrix<nsd_, 1> velint_old_;

      //! fluid velocity vector at integration point
      LINALG::Matrix<nsd_, 1> fluidvelint_;
      //! fluid velocity vector at integration point
      LINALG::Matrix<nsd_, 1> fluidvelint_old_;

      //! pressure at integration point
      double pressint_;
      //! pressure at integration point
      double pressint_old_;

      //! infinitesimal area element drs
      double drs_;
      //! integration factor
      double fac_;
      //! physical viscosity
      double visc_;
      //! density at t_(n+alpha_F) or t_(n+1)
      double dens_;

    };  // end class fluid3boundaryImpl

  }  // namespace ELEMENTS
}  // namespace DRT

#endif /* TOPOPT_FLUIDADJOINT3_BOUNDARY_H_ */
