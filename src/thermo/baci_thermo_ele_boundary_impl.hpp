/*----------------------------------------------------------------------*/
/*! \file

\brief internal implementation of thermo boundary elements

\level 1

*/

/*----------------------------------------------------------------------*
 | definitions                                               dano 09/09 |
 *----------------------------------------------------------------------*/
#ifndef BACI_THERMO_ELE_BOUNDARY_IMPL_HPP
#define BACI_THERMO_ELE_BOUNDARY_IMPL_HPP


/*----------------------------------------------------------------------*
 | headers                                                   dano 09/09 |
 *----------------------------------------------------------------------*/
#include "baci_config.hpp"

#include "baci_thermo_ele_impl_utils.hpp"
#include "baci_thermo_element.hpp"

BACI_NAMESPACE_OPEN


/*----------------------------------------------------------------------*
 |                                                           dano 09/09 |
 *----------------------------------------------------------------------*/
namespace DRT
{
  namespace ELEMENTS
  {
    /// Interface base class for TemperImpl
    //!
    //! This class exists to provide a common interface for all template
    //! versions of TemperImpl. The only function
    //! this class actually defines is Impl, which returns a pointer to
    //! the appropriate version of TemperImpl.
    class TemperBoundaryImplInterface
    {
     public:
      //! Empty constructor
      TemperBoundaryImplInterface() {}
      //! Empty destructor
      virtual ~TemperBoundaryImplInterface() = default;
      //! Evaluate the element in the case of multiple dofsets
      //!
      //! This class does not provide a definition for this function, it
      //! must be defined in TemperBoundaryImpl.
      virtual int Evaluate(DRT::ELEMENTS::ThermoBoundary* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, DRT::Element::LocationArray& la,
          CORE::LINALG::SerialDenseMatrix& elemat1_epetra,
          CORE::LINALG::SerialDenseMatrix& elemat2_epetra,
          CORE::LINALG::SerialDenseVector& elevec1_epetra,
          CORE::LINALG::SerialDenseVector& elevec2_epetra,
          CORE::LINALG::SerialDenseVector& elevec3_epetra) = 0;

      //! Evaluate a Neumann boundary condition
      //!
      //! This class does not provide a definition for this function, it
      //! must be defined in TemperBoundaryImpl.
      virtual int EvaluateNeumann(DRT::Element* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, DRT::Condition& condition, std::vector<int>& lm,
          CORE::LINALG::SerialDenseVector& elevec1_epetra) = 0;

      //! Internal implementation class for thermo elements
      static TemperBoundaryImplInterface* Impl(DRT::Element* ele);

    };  // TemperBoundaryImplInterface

    //! Internal thermo element implementation
    //!
    //! This internal class keeps all the working arrays needed to
    //! calculate the Thermo element. Additionally the method Sysmat()
    //! provides a clean and fast element implementation.
    //!
    //!  <h3>Purpose</h3>
    //!
    //! The idea is to separate the element maintenance (class Thermo)
    //! from the mathematical contents (this class). Of course there are
    //! different implementations of the Thermo element, this is just one
    //! such implementation.
    //!
    //! The Thermo element will allocate exactly one object of this class
    //! for all transport elements with the same number of nodes in the mesh.
    //! This allows us to use exactly matching working arrays (and keep them
    //! around.)
    //!
    //! The code is meant to be as clean as possible. This is the only way
    //! to keep it fast. The number of working arrays has to be reduced to
    //! a minimum so that the element fits into the cache. (There might be
    //! room for improvements.)
    //!
    //! <h3>History</h3>
    //!
    //! The implementation here is the standard convection-diffusion element
    //! capable of dealing with systems of transported scalars.
    //!
    //! Right now we do not read any stabilization parameters from the
    //! input file but have a fixed version.
    //!
    //! \author gjb
    //! \date 08/08
    template <CORE::FE::CellType distype>
    class TemperBoundaryImpl : public TemperBoundaryImplInterface
    {
     public:
      /// Constructor
      TemperBoundaryImpl(int numdofpernode);

      //! number of nodes
      static constexpr int nen_ = CORE::FE::num_nodes<distype>;

      //! number of space dimensions
      static constexpr int nsd_ = CORE::FE::dim<distype>;

      //! number of Gauss points
      static constexpr int nquad_ = THR::DisTypeToNumGaussPoints<distype>::nquad;


      //! Evaluate (la required in case of multiple dofsets)
      int Evaluate(DRT::ELEMENTS::ThermoBoundary* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, DRT::Element::LocationArray& la,
          CORE::LINALG::SerialDenseMatrix& elemat1_epetra,
          CORE::LINALG::SerialDenseMatrix& elemat2_epetra,
          CORE::LINALG::SerialDenseVector& elevec1_epetra,
          CORE::LINALG::SerialDenseVector& elevec2_epetra,
          CORE::LINALG::SerialDenseVector& elevec3_epetra) override;

      //! Evaluate a Neumann boundary condition
      int EvaluateNeumann(DRT::Element* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, DRT::Condition& condition, std::vector<int>& lm,
          CORE::LINALG::SerialDenseVector& elevec1_epetra) override;

     private:
      //! prepare the evaluation of NURBS shape functions
      virtual void PrepareNurbsEval(DRT::Element* ele,  //!< the element whose matrix is calculated
          DRT::Discretization& discretization           //!< current discretisation
      );

      //! evaluate shape functions and derivatives at int. point
      void EvalShapeFuncAndIntFac(
          const CORE::FE::IntPointsAndWeights<nsd_>& intpoints,  //!< integration points
          const int& iquad,                                      //!< id of current Gauss point
          const int& eleid                                       //!< the element id
      );

      //! integral of shape functions over boundary surface
      void IntegrateShapeFunctions(const DRT::Element* ele,  //!< the actual boundary element
          Teuchos::ParameterList& params,                    //!< the parameter list
          CORE::LINALG::SerialDenseVector& elevec1,          //!< result vector (to be assembled)
          const bool addarea                                 //!< flag for area calculation
      );

      //! Compute a constant normal vector for a boundary element
      void GetConstNormal(
          CORE::LINALG::Matrix<nsd_ + 1, 1>& normal,        //!< the constant normal vector
          const CORE::LINALG::Matrix<nsd_ + 1, nen_>& xyze  //!< element node coordinates
      );

      // evaluate a Neumann boundary condition
      // this method evaluates normal and detA at gaussian point
      // deriv (in)  : derivatives of shape functions
      void SurfaceIntegration(double& detA,  //!< area at GP
          CORE::LINALG::Matrix<nsd_ + 1, 1>&
              normal,  //!< (o) normal at gaussian point, length is detA!
          const CORE::LINALG::Matrix<nen_, nsd_ + 1>& xcurr  //!< (i)current coordinates of nodes
                                                             //!< nodal coords in either material
                                                             //!< or spatial frame
      );

      //! evaluate a heat convection boundary condition
      //!
      //! This is a non-linear heat load type, due to its dependence on the
      //! current heat state (follower load-like).
      //! Using the current temperature solution requires the linearisation.
      void CalculateConvectionFintCond(const DRT::Element* ele,  //!< the actual boundary element
          CORE::LINALG::Matrix<nen_, nen_>* econd,               // view only!
                                                                 //!< tangent of the thermal problem
          CORE::LINALG::Matrix<nen_, 1>* efext,                  // view only!
                                                                 //!< heat flux to be applied
          const double coeff,                                    //!< heat transfer coefficient
          const double surtemp,                                  //!< surrounding temperature
          const std::string tempstate                            //!< desired temperature state
      );

      //! evaluate a geometrically nonlinear heat convection boundary condition
      //!
      //! This is a non-linear heat load type, due to its dependence on the
      //! current heat AND the current deformation state (follower load-like).
      //!
      //! Using the current temperature solution T_{n+1} requires the linearisation.
      //! In addition: using the current displacements u_{n+1} requires linearisation
      void CalculateNlnConvectionFintCond(DRT::Element* ele,  //!< the actual boundary element
          std::vector<double>& disp,                          //!< current displacements d_{n+1}
          CORE::LINALG::Matrix<nen_, nen_>* econd,            // view only!
                                                    //!< tangent of the thermal problem k_TT
          CORE::LINALG::Matrix<nen_, (nsd_ + 1) * nen_>* etangcoupl,  // view only!
                                                                      //!< coupling tangent k_Td
          CORE::LINALG::Matrix<nen_, 1>* efext,                       // view only!
                                                                      //!< heat flux to be applied
          const double coeff,                                         //!< heat transfer coefficient
          const double surtemp,                                       //!< surrounding temperature
          const std::string tempstate                                 //!< desired temperature state
      );

      //! actual values of temperatures
      CORE::LINALG::Matrix<nen_, 1> etemp_;

      //! number of dof per node
      const int numdofpernode_;
      //! node coordinates
      CORE::LINALG::Matrix<nsd_ + 1, nen_> xyze_;
      //! coordinates of current integration point in reference coordinates
      CORE::LINALG::Matrix<nsd_, 1> xsi_;
      //! array for shape functions
      CORE::LINALG::Matrix<nen_, 1> funct_;
      //! array for shape function derivatives w.r.t. r,s,t
      CORE::LINALG::Matrix<nsd_, nen_> deriv_;
      //! global derivatives of shape functions w.r.t. x,y,z
      CORE::LINALG::Matrix<nsd_, nen_> derxy_;
      //! unit normal vector at integration point
      CORE::LINALG::Matrix<nsd_ + 1, 1> normal_;
      //! metric tensor at integration point
      CORE::LINALG::Matrix<nsd_, nsd_> metrictensor_;
      //! integration factor for current GP: fac = GaussWeight * drs
      double fac_;

      //! nurbs specific: element knots
      std::vector<CORE::LINALG::SerialDenseVector> myknots_;
      //! nurbs specific: control point weights
      CORE::LINALG::Matrix<nen_, 1> weights_;
      // normal fac
      double normalfac_;

    };  // TemperBoundaryImpl

  }  // namespace ELEMENTS

}  // namespace DRT


/*----------------------------------------------------------------------*/
BACI_NAMESPACE_CLOSE

#endif  // THERMO_ELE_BOUNDARY_IMPL_H
