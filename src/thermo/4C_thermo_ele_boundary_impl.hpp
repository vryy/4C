/*----------------------------------------------------------------------*/
/*! \file

\brief internal implementation of thermo boundary elements

\level 1

*/

/*----------------------------------------------------------------------*
 | definitions                                               dano 09/09 |
 *----------------------------------------------------------------------*/
#ifndef FOUR_C_THERMO_ELE_BOUNDARY_IMPL_HPP
#define FOUR_C_THERMO_ELE_BOUNDARY_IMPL_HPP


/*----------------------------------------------------------------------*
 | headers                                                   dano 09/09 |
 *----------------------------------------------------------------------*/
#include "4C_config.hpp"

#include "4C_thermo_ele_impl_utils.hpp"
#include "4C_thermo_element.hpp"

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*
 |                                                           dano 09/09 |
 *----------------------------------------------------------------------*/
namespace Discret
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
      virtual int Evaluate(const Discret::ELEMENTS::ThermoBoundary* ele,
          Teuchos::ParameterList& params, const Discret::Discretization& discretization,
          const Core::Elements::Element::LocationArray& la,
          Core::LinAlg::SerialDenseMatrix& elemat1_epetra,
          Core::LinAlg::SerialDenseMatrix& elemat2_epetra,
          Core::LinAlg::SerialDenseVector& elevec1_epetra,
          Core::LinAlg::SerialDenseVector& elevec2_epetra,
          Core::LinAlg::SerialDenseVector& elevec3_epetra) = 0;

      //! Evaluate a Neumann boundary condition
      //!
      //! This class does not provide a definition for this function, it
      //! must be defined in TemperBoundaryImpl.
      virtual int evaluate_neumann(const Core::Elements::Element* ele,
          Teuchos::ParameterList& params, const Discret::Discretization& discretization,
          const Core::Conditions::Condition& condition, const std::vector<int>& lm,
          Core::LinAlg::SerialDenseVector& elevec1_epetra) = 0;

      //! Internal implementation class for thermo elements
      static TemperBoundaryImplInterface* Impl(const Core::Elements::Element* ele);

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
    template <Core::FE::CellType distype>
    class TemperBoundaryImpl : public TemperBoundaryImplInterface
    {
     public:
      /// Constructor
      TemperBoundaryImpl(int numdofpernode);

      //! number of nodes
      static constexpr int nen_ = Core::FE::num_nodes<distype>;

      //! number of space dimensions
      static constexpr int nsd_ = Core::FE::dim<distype>;

      //! number of Gauss points
      static constexpr int nquad_ = THR::DisTypeToNumGaussPoints<distype>::nquad;


      //! Evaluate (la required in case of multiple dofsets)
      int Evaluate(const Discret::ELEMENTS::ThermoBoundary* ele, Teuchos::ParameterList& params,
          const Discret::Discretization& discretization,
          const Core::Elements::Element::LocationArray& la,
          Core::LinAlg::SerialDenseMatrix& elemat1_epetra,
          Core::LinAlg::SerialDenseMatrix& elemat2_epetra,
          Core::LinAlg::SerialDenseVector& elevec1_epetra,
          Core::LinAlg::SerialDenseVector& elevec2_epetra,
          Core::LinAlg::SerialDenseVector& elevec3_epetra) override;

      //! Evaluate a Neumann boundary condition
      int evaluate_neumann(const Core::Elements::Element* ele, Teuchos::ParameterList& params,
          const Discret::Discretization& discretization,
          const Core::Conditions::Condition& condition, const std::vector<int>& lm,
          Core::LinAlg::SerialDenseVector& elevec1_epetra) override;

     private:
      //! prepare the evaluation of NURBS shape functions
      virtual void prepare_nurbs_eval(
          const Core::Elements::Element* ele,            //!< the element whose matrix is calculated
          const Discret::Discretization& discretization  //!< current discretisation
      );

      //! evaluate shape functions and derivatives at int. point
      void eval_shape_func_and_int_fac(
          const Core::FE::IntPointsAndWeights<nsd_>& intpoints,  //!< integration points
          const int& iquad,                                      //!< id of current Gauss point
          const int& eleid                                       //!< the element id
      );

      //! integral of shape functions over boundary surface
      void integrate_shape_functions(
          const Core::Elements::Element* ele,        //!< the actual boundary element
          Teuchos::ParameterList& params,            //!< the parameter list
          Core::LinAlg::SerialDenseVector& elevec1,  //!< result vector (to be assembled)
          const bool addarea                         //!< flag for area calculation
      );

      //! Compute a constant normal vector for a boundary element
      void get_const_normal(
          Core::LinAlg::Matrix<nsd_ + 1, 1>& normal,        //!< the constant normal vector
          const Core::LinAlg::Matrix<nsd_ + 1, nen_>& xyze  //!< element node coordinates
      ) const;

      // evaluate a Neumann boundary condition
      // this method evaluates normal and detA at gaussian point
      // deriv (in)  : derivatives of shape functions
      void surface_integration(double& detA,  //!< area at GP
          Core::LinAlg::Matrix<nsd_ + 1, 1>&
              normal,  //!< (o) normal at gaussian point, length is detA!
          const Core::LinAlg::Matrix<nen_, nsd_ + 1>& xcurr  //!< (i)current coordinates of nodes
                                                             //!< nodal coords in either material
                                                             //!< or spatial frame
      );

      //! evaluate a heat convection boundary condition
      //!
      //! This is a non-linear heat load type, due to its dependence on the
      //! current heat state (follower load-like).
      //! Using the current temperature solution requires the linearisation.
      void calculate_convection_fint_cond(
          const Core::Elements::Element* ele,       //!< the actual boundary element
          Core::LinAlg::Matrix<nen_, nen_>* econd,  // view only!
                                                    //!< tangent of the thermal problem
          Core::LinAlg::Matrix<nen_, 1>* efext,     // view only!
                                                    //!< heat flux to be applied
          const double coeff,                       //!< heat transfer coefficient
          const double surtemp,                     //!< surrounding temperature
          const std::string& tempstate              //!< desired temperature state
      );

      //! evaluate a geometrically nonlinear heat convection boundary condition
      //!
      //! This is a non-linear heat load type, due to its dependence on the
      //! current heat AND the current deformation state (follower load-like).
      //!
      //! Using the current temperature solution T_{n+1} requires the linearisation.
      //! In addition: using the current displacements u_{n+1} requires linearisation
      void calculate_nln_convection_fint_cond(
          const Core::Elements::Element* ele,       //!< the actual boundary element
          const std::vector<double>& disp,          //!< current displacements d_{n+1}
          Core::LinAlg::Matrix<nen_, nen_>* econd,  // view only!
                                                    //!< tangent of the thermal problem k_TT
          Core::LinAlg::Matrix<nen_, (nsd_ + 1) * nen_>* etangcoupl,  // view only!
                                                                      //!< coupling tangent k_Td
          Core::LinAlg::Matrix<nen_, 1>* efext,                       // view only!
                                                                      //!< heat flux to be applied
          const double coeff,                                         //!< heat transfer coefficient
          const double surtemp,                                       //!< surrounding temperature
          const std::string& tempstate                                //!< desired temperature state
      );

      //! actual values of temperatures
      Core::LinAlg::Matrix<nen_, 1> etemp_;

      //! number of dof per node
      const int numdofpernode_;
      //! node coordinates
      Core::LinAlg::Matrix<nsd_ + 1, nen_> xyze_;
      //! coordinates of current integration point in reference coordinates
      Core::LinAlg::Matrix<nsd_, 1> xsi_;
      //! array for shape functions
      Core::LinAlg::Matrix<nen_, 1> funct_;
      //! array for shape function derivatives w.r.t. r,s,t
      Core::LinAlg::Matrix<nsd_, nen_> deriv_;
      //! global derivatives of shape functions w.r.t. x,y,z
      Core::LinAlg::Matrix<nsd_, nen_> derxy_;
      //! unit normal vector at integration point
      Core::LinAlg::Matrix<nsd_ + 1, 1> normal_;
      //! metric tensor at integration point
      Core::LinAlg::Matrix<nsd_, nsd_> metrictensor_;
      //! integration factor for current GP: fac = GaussWeight * drs
      double fac_;

      //! nurbs specific: element knots
      std::vector<Core::LinAlg::SerialDenseVector> myknots_;
      //! nurbs specific: control point weights
      Core::LinAlg::Matrix<nen_, 1> weights_;
      // normal fac
      double normalfac_;

    };  // TemperBoundaryImpl

  }  // namespace ELEMENTS

}  // namespace Discret


/*----------------------------------------------------------------------*/
FOUR_C_NAMESPACE_CLOSE

#endif
