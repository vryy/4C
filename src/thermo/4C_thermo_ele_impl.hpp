/*----------------------------------------------------------------------*/
/*! \file

\brief Internal implementation of thermo elements

\level 1


*/

/*----------------------------------------------------------------------*
 | definitions                                                gjb 01/08 |
 *----------------------------------------------------------------------*/
#ifndef FOUR_C_THERMO_ELE_IMPL_HPP
#define FOUR_C_THERMO_ELE_IMPL_HPP

/*----------------------------------------------------------------------*
 | headers                                                    gjb 01/08 |
 *----------------------------------------------------------------------*/
#include "4C_config.hpp"

#include "4C_discretization_fem_general_utils_local_connectivity_matrices.hpp"
#include "4C_inpar_thermo.hpp"
#include "4C_thermo_ele_impl_utils.hpp"
#include "4C_utils_singleton_owner.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 |                                                            gjb 01/08 |
 *----------------------------------------------------------------------*/
namespace DRT
{
  namespace ELEMENTS
  {
    //! Interface base class for TemperImpl
    //!
    //!  This class exists to provide a common interface for all template
    //!  versions of TemperImpl. The only function
    //!  this class actually defines is Impl, which returns a pointer to
    //!  the appropriate version of TemperImpl.
    class TemperImplInterface
    {
     public:
      //! Empty constructor
      TemperImplInterface() = default;

      //! Virtual destructor.
      virtual ~TemperImplInterface() = default;

      //! Evaluate the element
      //!  This class does not provide a definition for this function, it
      //!  must be defined in TemperImpl.
      virtual int Evaluate(const DRT::Element* ele,         //!< current element
          Teuchos::ParameterList& params,                   //!< parameter list, containing e.g., dt
          const DRT::Discretization& discretization,        //!< current discretisation
          const DRT::Element::LocationArray& la,            //!< location array
          CORE::LINALG::SerialDenseMatrix& elemat1_epetra,  //!< conductivity matrix
          CORE::LINALG::SerialDenseMatrix& elemat2_epetra,  //!< capacity matrix
          CORE::LINALG::SerialDenseVector&
              elevec1_epetra,  //!< internal force, view on heat flux in x-direction
          CORE::LINALG::SerialDenseVector&
              elevec2_epetra,  //!< external force, view on heat flux in y-direction
          CORE::LINALG::SerialDenseVector& elevec3_epetra  //!< view on heat flux in z-direction
          ) = 0;

      //! Evaluate the element
      //!
      //!  This class does not provide a definition for this function, it
      //!  must be defined in TemperImpl.
      virtual int EvaluateNeumann(const DRT::Element* ele,  //!< current element
          const Teuchos::ParameterList& params,             //!< parameter list
          const DRT::Discretization& discretization,        //!< current discretisation
          const std::vector<int>&
              lm,  //!< location vector, EvalNeumann is called only on own discretisation
          CORE::LINALG::SerialDenseVector& elevec1_epetra,  //!< view on external force vector
          CORE::LINALG::SerialDenseMatrix* elemat1_epetra   //!< matrix is not needed
          ) = 0;

      //! Internal implementation class for thermo elements
      static TemperImplInterface* Impl(const DRT::Element* ele);

    };  // class TemperImplInterface


    //! Internal Thermo element implementation
    //!
    //!  This internal class keeps all the working arrays needed to
    //!  calculate the thermo element. Additionally the method Sysmat()
    //!  provides a clean and fast element implementation.
    //!
    //!  <h3>Purpose</h3>
    //!
    //!  The idea is to separate the element maintenance (class Thermo)
    //!  from the mathematical contents (this class). Of course there are
    //!  different implementations of the Thermo element, this is just one
    //!  such implementation.
    //!
    //!  The Thermo element will allocate exactly one object of this class
    //!  for all thermo elements with the same number of nodes in the mesh.
    //!  This allows us to use exactly matching working arrays (and keep them
    //!  around.)
    //!
    //!  The code is meant to be as clean as possible. This is the only way
    //!  to keep it fast. The number of working arrays has to be reduced to
    //!  a minimum so that the element fits into the cache. (There might be
    //!  room for improvements.)
    //!
    //!  <h3>History</h3>
    //!
    //!  \author dano
    //!  \date 09/09
    //!
    template <CORE::FE::CellType distype>
    class TemperImpl : public TemperImplInterface
    {
     public:
      //! Constructor
      TemperImpl();

      //! Singleton access method
      static TemperImpl<distype>* Instance(
          CORE::UTILS::SingletonAction action = CORE::UTILS::SingletonAction::create);

      //! number of nodes
      static constexpr int nen_ = CORE::FE::num_nodes<distype>;

      //! number of space dimensions
      static constexpr int nsd_ = CORE::FE::dim<distype>;

      //! number of dof per node
      static constexpr int numdofpernode_ = 1;

      //! number of Gauss points
      static constexpr int nquad_ = THR::DisTypeToNumGaussPoints<distype>::nquad;

      //! Evaluate for multiple dofsets
      int Evaluate(const DRT::Element* ele,                 //!< current element
          Teuchos::ParameterList& params,                   //!< parameter list, containing e.g., dt
          const DRT::Discretization& discretization,        //!< current discretisation
          const DRT::Element::LocationArray& la,            //!< location array
          CORE::LINALG::SerialDenseMatrix& elemat1_epetra,  //!< conductivity matrix
          CORE::LINALG::SerialDenseMatrix& elemat2_epetra,  //!< capacity matrix
          CORE::LINALG::SerialDenseVector&
              elevec1_epetra,  //!< internal force, view on heat flux in x-direction
          CORE::LINALG::SerialDenseVector&
              elevec2_epetra,  //!< external force, view on heat flux in y-direction
          CORE::LINALG::SerialDenseVector& elevec3_epetra  //!< view on heat flux in z-direction
          ) override;

      //! Evaluate the element
      int EvaluateNeumann(const DRT::Element* ele,    //!< current element
          const Teuchos::ParameterList& params,       //!< parameter list
          const DRT::Discretization& discretization,  //!< current discretisation
          const std::vector<int>&
              lm,  //!< location vector, EvalNeumann is called only on own discretisation
          CORE::LINALG::SerialDenseVector& elevec1_epetra,  //!< view on external force vector
          CORE::LINALG::SerialDenseMatrix* elemat1_epetra   //!< matrix is not needed
          ) override;

     private:
      //! @name general thermal functions
      //! @{

      //! calculate complete internal force, tangent matrix k_TT and capacity matrix
      //!
      //! builds quantitites from linear/nonlinear and thermo/TSI specific routines
      void evaluate_tang_capa_fint(const Element* ele, const double time,
          const Discretization& discretization, const Element::LocationArray& la,
          CORE::LINALG::Matrix<nen_ * numdofpernode_, nen_ * numdofpernode_>* etang,
          CORE::LINALG::Matrix<nen_ * numdofpernode_, nen_ * numdofpernode_>* ecapa,
          CORE::LINALG::Matrix<nen_ * numdofpernode_, nen_ * numdofpernode_>* ecapalin,
          CORE::LINALG::Matrix<nen_ * numdofpernode_, 1>* efint, Teuchos::ParameterList& params);

      /*!
       * evaluate complete coupled tangent matrix k_Td
       * @param ele the element whose matrix is calculated
       * @param discretization discretization containing this element
       * @param la LocationArray of this element inside discretization
       * @param etangcoupl matrix k_Td to be filled
       * @param params ParameterList of options
       */
      void EvaluateCoupledTang(const Element* ele, const Discretization& discretization,
          const Element::LocationArray& la,
          CORE::LINALG::Matrix<nen_ * numdofpernode_, nen_ * nsd_ * numdofpernode_>* etangcoupl,
          Teuchos::ParameterList& params);

      /*!
       * evaluate external body loads
       * @param ele the element whose matrix is calculated
       * @param time time for function evaluation
       * @param efext external force vector
       */
      void EvaluateFext(const Element* ele, const double time,
          CORE::LINALG::Matrix<nen_ * numdofpernode_, 1>& efext);

      //! Calculate element force vectors and a few matrices
      void linear_thermo_contribution(
          const DRT::Element* ele,  //!< the element whose matrix is calculated
          const double time,        //!< current time
          CORE::LINALG::Matrix<nen_ * numdofpernode_, nen_ * numdofpernode_>*
              econd,  //!< conductivity matrix
          CORE::LINALG::Matrix<nen_ * numdofpernode_, nen_ * numdofpernode_>*
              ecapa,  //!< capacity matrix
          CORE::LINALG::Matrix<nen_ * numdofpernode_, nen_ * numdofpernode_>*
              ecapalin,  //!< linearization contribution of capacity
          CORE::LINALG::Matrix<nen_ * numdofpernode_, 1>* efint  //!< internal force
      );


      //! @}

      //! @name geometrically linear TSI
      //! @{

      //! Calculate element vectors (internal/external) and a few matrices
      //! considering current displacement solution
      void linear_disp_contribution(const DRT::Element* ele, const double time,
          const std::vector<double>& disp, const std::vector<double>& vel,
          CORE::LINALG::Matrix<nen_ * numdofpernode_, nen_ * numdofpernode_>* econd,
          CORE::LINALG::Matrix<nen_ * numdofpernode_, 1>* efint,
          const Teuchos::ParameterList& params);

      //! calculate thermal-mechanical system matrix term needed in monolithic TSI
      void LinearCoupledTang(const DRT::Element* ele,  //!< the element whose matrix is calculated
          const std::vector<double>& disp,             //!< current displacements
          const std::vector<double>& vel,              //!< current velocities
          CORE::LINALG::Matrix<nen_ * numdofpernode_, nsd_ * nen_ * numdofpernode_>*
              etangcoupl,  //!< k_Tu matrix
          const Teuchos::ParameterList& params);

      //! @}

      //! @name linear, small strain thermoplasticity solved with TSI
      //! @{

      //! calculate internal dissipation arising when a thermo-elasto-plastic
      //! material is used
      //! Clausius-Duhem inequality is no longer = 0, but > 0:
      //! mechanical energy dissipates as heat
      void linear_dissipation_fint(
          const DRT::Element* ele,  //!< the element whose matrix is calculated
          CORE::LINALG::Matrix<nen_ * numdofpernode_, 1>* efint,  //!< internal force
          Teuchos::ParameterList& params);

      //! calculate terms of dissipation for thermo-mechanical
      //! system matrix k_Td used in case of plastic material
      void linear_dissipation_coupled_tang(
          const DRT::Element* ele,  // the element whose matrix is calculated
          CORE::LINALG::Matrix<nen_ * numdofpernode_, nsd_ * nen_ * numdofpernode_>*
              etangcoupl,  // k_Td
          Teuchos::ParameterList& params);

      //! @}

      //! @name geometrically nonlinear TSI analysis
      //! @{

      //! calculate element vectors (internal/external) and a few matrices
      //! considering current displacement solution
      //! --> all terms are coupled to the displacements/velocities
      void nonlinear_thermo_disp_contribution(
          const DRT::Element* ele,          //!< the element whose matrix is calculated
          const double time,                //!< current time
          const std::vector<double>& disp,  //!< current displacements
          const std::vector<double>& vel,   //!< current velocities
          CORE::LINALG::Matrix<nen_ * numdofpernode_, nen_ * numdofpernode_>*
              econd,  //!< conductivity matrix
          CORE::LINALG::Matrix<nen_ * numdofpernode_, nen_ * numdofpernode_>*
              ecapa,  //!< capacity matrix
          CORE::LINALG::Matrix<nen_ * numdofpernode_, nen_ * numdofpernode_>*
              ecapalin,  //!< partial linearization dC/dT of capacity matrix
          CORE::LINALG::Matrix<nen_ * numdofpernode_, 1>* efint,  //!< internal force
          Teuchos::ParameterList& params);

      //! calculate thermal-mechanical system matrix k_Td needed in monolithic TSI
      void nonlinear_coupled_tang(
          const DRT::Element* ele,          //!< current element whose terms are calculated
          const std::vector<double>& disp,  //!< current displacements
          const std::vector<double>& vel,   //!< current velocities
          CORE::LINALG::Matrix<nen_ * numdofpernode_, nsd_ * nen_ * numdofpernode_>*
              etangcoupl,                 //!< k_Tu matrix
          Teuchos::ParameterList& params  //!< parameter list, containing e.g., dt,theta
      );

      //! build nonlinear B-operator
      void CalculateBop(
          CORE::LINALG::Matrix<6, nsd_ * nen_ * numdofpernode_>* bop,  //!< nonlinear B-operator
          const CORE::LINALG::Matrix<nsd_, nsd_>* defgrd,              //!< deformation gradient
          const CORE::LINALG::Matrix<nsd_, nen_>* N_XYZ                //!< gradient-operator
      ) const;

      //! build linearisation of Jacobian w.r.t. d: dJ_dd
      void calculate_linearisation_of_jacobian(
          CORE::LINALG::Matrix<1, nsd_ * nen_ * numdofpernode_>& dJ_dd,  //!<  [out] dJ_dd
          const double J,                                                //!< Jacobian
          const CORE::LINALG::Matrix<nsd_, nen_>& N_XYZ,  //!< linear gradient of shape functions
          const CORE::LINALG::Matrix<nsd_, nsd_>& defgrd_inv  //!< inverse of F
      ) const;

      //! build derivatives of right Cauchy-Green deformation tensor C
      //! build the inverse of C^{-1} and the time derivative C'
      void calculate_cauchy_greens(
          CORE::LINALG::Matrix<6, 1>& Cratevct,            //!< right Cauchy-Green rate vector
          CORE::LINALG::Matrix<6, 1>& Cinvvct,             //!< inverse of right Cauchy-Green vector
          CORE::LINALG::Matrix<nsd_, nsd_>& Cinv,          //!< inverse right Cauchy-Green tensor
          const CORE::LINALG::Matrix<nsd_, nsd_>* defgrd,  //!< deformation gradient tensor
          const CORE::LINALG::Matrix<nsd_, nsd_>* defgrdrate,  //!< velocity gradient tensor
          const CORE::LINALG::Matrix<nsd_, nsd_>*
              invdefgrd  //!< inverse deformation gradient tensor
      ) const;

      /// @}

      //! @name finite strain thermoplasticity solved with TSI
      //! @{

      //! calculate internal dissipation arising when a thermo-elasto-plastic
      //! material is used within geometrically nonlinear analysis
      //! Clausius-Duhem inequality is no longer = 0, but > 0:
      //! mechanical energy dissipates as heat
      void nonlinear_dissipation_fint_tang(
          const DRT::Element* ele,          //!< the element whose matrix is calculated
          const std::vector<double>& disp,  //!< current displacements
          CORE::LINALG::Matrix<nen_ * numdofpernode_, nen_ * numdofpernode_>*
              econd,                                              //!< conductivity matrix
          CORE::LINALG::Matrix<nen_ * numdofpernode_, 1>* efint,  //!< internal force
          Teuchos::ParameterList& params);

      //! calculate terms of dissipation for thermo-mechanical system matrix k_Td
      //! used in case of plastic material within geometrically nonlinear analysis
      void nonlinear_dissipation_coupled_tang(
          const DRT::Element* ele,          //!< the element whose matrix is calculated
          const std::vector<double>& disp,  //!< current displacements
          const std::vector<double>& vel,   //!< current velocities
          CORE::LINALG::Matrix<nen_ * numdofpernode_, nsd_ * nen_ * numdofpernode_>*
              etangcoupl,  //!< k_Td
          Teuchos::ParameterList& params);

      /// @}

      //! get the body force
      virtual void Radiation(const DRT::Element* ele,  //!< current element we are dealing with
          const double time                            //!< current times
      );

      //! build linear B-operator
      void CalculateBoplin(
          CORE::LINALG::Matrix<6, nsd_ * nen_ * numdofpernode_>* boplin,  //!< linear B-operator
          const CORE::LINALG::Matrix<nsd_, nen_>* N_XYZ                   //!< gradient-operator
      ) const;

      //! get corresponding structural material
      Teuchos::RCP<CORE::MAT::Material> GetSTRMaterial(
          const DRT::Element* ele  //!< the element whose matrix is calculated
      ) const;

      //! calculate reactive term
      void calculate_reactive_term(
          const CORE::LINALG::Matrix<6, 1>* ctemp,     //!< temperature-dependent material tangent
          const CORE::LINALG::Matrix<6, 1>* strainvel  //!< strain rate
      ) const;

      //! determine heat flux and conductivity tensor
      //! based on material law
      virtual void Materialize(const DRT::Element* ele,  //!< current element
          int gp                                         //!< current GP id
      );

      //! evaluate shape functions and their derivatives at current integration point
      virtual void eval_shape_func_and_derivs_at_int_point(
          const CORE::FE::IntPointsAndWeights<nsd_>& intpoints,  //!< integration points
          int iquad,                                             //!< id of current Gauss point
          int eleid                                              //!< the element id
      );

      //! compute heatflux and temperature gradient in linear case
      void linear_heatflux_tempgrad(const Element* ele,   //!< the current element
          CORE::LINALG::Matrix<nquad_, nsd_>* eheatflux,  //!< [out] heat fluxes at Gauss points
          CORE::LINALG::Matrix<nquad_, nsd_>*
              etempgrad  //!< [out] temperature gradients at Gauss points
      );

      //! compute heatflux and temperature gradient in nonlinear case
      void nonlinear_heatflux_tempgrad(const Element* ele,  //!< the current element
          const std::vector<double>& disp,                  //!< element displacements
          const std::vector<double>& vel,                   //!< element velocities
          CORE::LINALG::Matrix<nquad_, nsd_>* eheatflux,    //!< [out] heat fluxes at Gauss points
          CORE::LINALG::Matrix<nquad_, nsd_>*
              etempgrad,                  //!< [out] temperature gradients at Gauss points
          Teuchos::ParameterList& params  //!< parameters containing type of flux and grad
      );

      //! calculate lumped capacity matrix in case of explicit time integration
      void CalculateLumpMatrix(
          CORE::LINALG::Matrix<nen_ * numdofpernode_, nen_ * numdofpernode_>* ecapa) const;

      //! calculate characteristic element length
      double calculate_char_ele_length() const;

      //! Compute the error compared to an analytical solution from input file
      virtual void ComputeError(
          const DRT::Element* ele,  //!< current element whose terms are calculated
          CORE::LINALG::Matrix<nen_ * numdofpernode_, 1>& elevec1,  //!< element vectorr
          Teuchos::ParameterList& params  //!< parameter list, containing analytical solution
      );

      //! Compute nodal position and velocity
      inline void initial_and_current_nodal_position_velocity(const Element* ele,
          const std::vector<double>& disp, const std::vector<double>& vel,
          CORE::LINALG::Matrix<nen_, nsd_>& xcurr, CORE::LINALG::Matrix<nen_, nsd_>& xcurrrate);

      //! prepare the evaluation of NURBS shape functions
      virtual void PrepareNurbsEval(
          const DRT::Element* ele,                   //!< the element whose matrix is calculated
          const DRT::Discretization& discretization  //!< current discretisation
      );

      //! integral of shape functions over the element
      void integrate_shape_functions(const DRT::Element* ele,  //!< current element
          CORE::LINALG::SerialDenseVector& elevec1,            //!< result vector (to be assembled)
          const CORE::LINALG::IntSerialDenseVector& dofids  //!< for which dof we need to integrate?
      );

      //! extrapolate from Gauss points to nodes, needed for postprocessing
      void extrapolate_from_gauss_points_to_nodes(const DRT::Element* ele,  //!< the actual element
          const CORE::LINALG::Matrix<nquad_, nsd_>& gpheatflux,  //!< heat flux at each Gauss Point
          CORE::LINALG::Matrix<nen_ * numdofpernode_, 1>&
              efluxx,  //!< element heat flux in x-direction
          CORE::LINALG::Matrix<nen_ * numdofpernode_, 1>&
              efluxy,  //!< element heat flux in y-direction
          CORE::LINALG::Matrix<nen_ * numdofpernode_, 1>&
              efluxz  //!< element heat flux in z-direction
      );

      //! extract displacement and velocity vector from discretization
      void ExtractDispVel(const Discretization& discretization, const Element::LocationArray& la,
          std::vector<double>& mydisp, std::vector<double>& myvel) const;

      //! copy matrix contents into character vector
      void copy_matrix_into_char_vector(
          std::vector<char>& data, const CORE::LINALG::Matrix<nquad_, nsd_>& stuff) const;

      //! FDcheck of conductivity matrix on element level
      void fd_check_coupl_nln_fint_cond_capa(
          const DRT::Element* ele,          //!< the element whose matrix is calculated
          const double time,                //!< current time
          const std::vector<double>& disp,  //!< current displacements
          const std::vector<double>& vel,   //!< current velocities
          CORE::LINALG::Matrix<nen_ * numdofpernode_, nen_ * numdofpernode_>*
              etang,                                              //!< tangent conductivity matrix
          CORE::LINALG::Matrix<nen_ * numdofpernode_, 1>* efint,  //!< internal force);)
          Teuchos::ParameterList& params) const;

      //! FDcheck of linearized capacity matrix on element level
      void FDCheckCapalin(const DRT::Element* ele,  //!< the element whose matrix is calculated
          const double time,                        //!< current time
          const std::vector<double>& disp,          //!< current displacements
          const std::vector<double>& vel,           //!< current velocities
          CORE::LINALG::Matrix<nen_ * numdofpernode_, nen_ * numdofpernode_>*
              ecapa,  //!< capacity matrix
          CORE::LINALG::Matrix<nen_ * numdofpernode_, nen_ * numdofpernode_>*
              ecapalin,  //!< linearization term from capacity matrix
          Teuchos::ParameterList& params) const;

      //! actual values of temperatures T_{n+1}
      CORE::LINALG::Matrix<nen_, 1> etempn_;
      //! temperatures in last time step T_{n}
      CORE::LINALG::Matrix<nen_, 1> etemp_;

      //! node reference coordinates
      CORE::LINALG::Matrix<nsd_, nen_> xyze_;
      //! radiation in element nodes
      CORE::LINALG::Matrix<numdofpernode_, 1> radiation_;
      //! coordinates of current integration point in reference coordinates
      CORE::LINALG::Matrix<nsd_, 1> xsi_;
      //! array for shape functions
      CORE::LINALG::Matrix<nen_, 1> funct_;
      //! array for shape function derivatives w.r.t r,s,t
      CORE::LINALG::Matrix<nsd_, nen_> deriv_;
      //! transposed jacobian "dx/ds"
      CORE::LINALG::Matrix<nsd_, nsd_> xjm_;
      //! inverse of transposed jacobian "ds/dx"
      CORE::LINALG::Matrix<nsd_, nsd_> xij_;
      //! global derivatives of shape functions w.r.t x,y,z
      CORE::LINALG::Matrix<nsd_, nen_> derxy_;
      //! integration factor for current GP: fac = GaussWeight * det(J)
      double fac_;
      //! (global) gradient of temperature at integration point
      CORE::LINALG::Matrix<nsd_, 1> gradtemp_;
      //! (global) gradient of temperature at integration point
      CORE::LINALG::Matrix<nsd_, 1> heatflux_;
      //! (global) conductivity 2-tensor
      CORE::LINALG::Matrix<nsd_, nsd_> cmat_;
      //! (global) derivative of conductivity 2-tensor w.r.t. T
      CORE::LINALG::Matrix<nsd_, nsd_> dercmat_;
      //! capacity density
      double capacoeff_;
      //! derivative of capacity w.r.t. T
      double dercapa_;

      //! @name material related stuff
      //! @{

      //! flag plastic material is used
      bool plasticmat_;

      //! nurbs specific: element knots
      std::vector<CORE::LINALG::SerialDenseVector> myknots_;
      //! nurbs specific: control point weights
      CORE::LINALG::Matrix<nen_, 1> weights_;

      /// @}

    };  // class TemperImpl

  }  // namespace ELEMENTS

}  // namespace DRT


/*----------------------------------------------------------------------*/
FOUR_C_NAMESPACE_CLOSE

#endif
