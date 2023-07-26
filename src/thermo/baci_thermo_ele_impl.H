/*----------------------------------------------------------------------*/
/*! \file

\brief Internal implementation of thermo elements

\level 1


*/

/*----------------------------------------------------------------------*
 | definitions                                                gjb 01/08 |
 *----------------------------------------------------------------------*/
#ifndef BACI_THERMO_ELE_IMPL_H
#define BACI_THERMO_ELE_IMPL_H

/*----------------------------------------------------------------------*
 | headers                                                    gjb 01/08 |
 *----------------------------------------------------------------------*/
#include "baci_thermo_ele_impl_utils.H"
#include "baci_inpar_thermo.H"

#include "baci_utils_singleton_owner.H"
#include "baci_discretization_fem_general_utils_local_connectivity_matrices.H"

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
      virtual int Evaluate(DRT::Element* ele,               //!< current element
          Teuchos::ParameterList& params,                   //!< parameter list, containing e.g., dt
          DRT::Discretization& discretization,              //!< current discretisation
          DRT::Element::LocationArray& la,                  //!< location array
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
      virtual int EvaluateNeumann(DRT::Element* ele,  //!< current element
          Teuchos::ParameterList& params,             //!< parameter list
          DRT::Discretization& discretization,        //!< current discretisation
          std::vector<int>&
              lm,  //!< location vector, EvalNeumann is called only on own discretisation
          CORE::LINALG::SerialDenseVector& elevec1_epetra,  //!< view on external force vector
          CORE::LINALG::SerialDenseMatrix* elemat1_epetra   //!< matrix is not needed
          ) = 0;

      //! Internal implementation class for thermo elements
      static TemperImplInterface* Impl(DRT::Element* ele);

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
    template <DRT::Element::DiscretizationType distype>
    class TemperImpl : public TemperImplInterface
    {
     public:
      //! Constructor
      TemperImpl();

      //! Singleton access method
      static TemperImpl<distype>* Instance(
          CORE::UTILS::SingletonAction action = CORE::UTILS::SingletonAction::create);

      //! number of nodes
      static const int nen_ = CORE::DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;

      //! number of space dimensions
      static const int nsd_ = CORE::DRT::UTILS::DisTypeToDim<distype>::dim;

      //! number of dof per node
      static const int numdofpernode_ = 1;

      //! number of Gauss points
      static const int nquad_ = THR::DisTypeToNumGaussPoints<distype>::nquad;

      //! Evaluate for multiple dofsets
      int Evaluate(DRT::Element* ele,                       //!< current element
          Teuchos::ParameterList& params,                   //!< parameter list, containing e.g., dt
          DRT::Discretization& discretization,              //!< current discretisation
          DRT::Element::LocationArray& la,                  //!< location array
          CORE::LINALG::SerialDenseMatrix& elemat1_epetra,  //!< conductivity matrix
          CORE::LINALG::SerialDenseMatrix& elemat2_epetra,  //!< capacity matrix
          CORE::LINALG::SerialDenseVector&
              elevec1_epetra,  //!< internal force, view on heat flux in x-direction
          CORE::LINALG::SerialDenseVector&
              elevec2_epetra,  //!< external force, view on heat flux in y-direction
          CORE::LINALG::SerialDenseVector& elevec3_epetra  //!< view on heat flux in z-direction
          ) override;

      //! Evaluate the element
      int EvaluateNeumann(DRT::Element* ele,    //!< current element
          Teuchos::ParameterList& params,       //!< parameter list
          DRT::Discretization& discretization,  //!< current discretisation
          std::vector<int>&
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
      void EvaluateTangCapaFint(Element* ele, const double& time, Discretization& discretization,
          Element::LocationArray& la,
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
      void EvaluateCoupledTang(Element* ele, const Discretization& discretization,
          Element::LocationArray& la,
          CORE::LINALG::Matrix<nen_ * numdofpernode_, nen_ * nsd_ * numdofpernode_>* etangcoupl,
          Teuchos::ParameterList& params);

      /*!
       * evaluate external body loads
       * @param ele the element whose matrix is calculated
       * @param time time for function evaluation
       * @param efext external force vector
       */
      void EvaluateFext(
          Element* ele, const double& time, CORE::LINALG::Matrix<nen_ * numdofpernode_, 1>& efext);

      //! Calculate element force vectors and a few matrices
      void LinearThermoContribution(DRT::Element* ele,  //!< the element whose matrix is calculated
          const double& time,                           //!< current time
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
      void LinearDispContribution(DRT::Element* ele, const double& time, std::vector<double>& disp,
          std::vector<double>& vel,
          CORE::LINALG::Matrix<nen_ * numdofpernode_, nen_ * numdofpernode_>* econd,
          CORE::LINALG::Matrix<nen_ * numdofpernode_, 1>* efint, Teuchos::ParameterList& params);

      //! calculate thermal-mechanical system matrix term needed in monolithic TSI
      void LinearCoupledTang(DRT::Element* ele,  //!< the element whose matrix is calculated
          std::vector<double>& disp,             //!< current displacements
          std::vector<double>& vel,              //!< current velocities
          CORE::LINALG::Matrix<nen_ * numdofpernode_, nsd_ * nen_ * numdofpernode_>*
              etangcoupl,  //!< k_Tu matrix
          Teuchos::ParameterList& params);

      //! @}

      //! @name linear, small strain thermoplasticity solved with TSI
      //! @{

      //! calculate internal dissipation arising when a thermo-elasto-plastic
      //! material is used
      //! Clausius-Duhem inequality is no longer = 0, but > 0:
      //! mechanical energy dissipates as heat
      void LinearDissipationFint(DRT::Element* ele,  //!< the element whose matrix is calculated
          CORE::LINALG::Matrix<nen_ * numdofpernode_, 1>* efint,  //!< internal force
          Teuchos::ParameterList& params);

      //! calculate terms of dissipation for thermo-mechanical
      //! system matrix k_Td used in case of plastic material
      void LinearDissipationCoupledTang(
          DRT::Element* ele,  // the element whose matrix is calculated
          CORE::LINALG::Matrix<nen_ * numdofpernode_, nsd_ * nen_ * numdofpernode_>*
              etangcoupl,  // k_Td
          Teuchos::ParameterList& params);

      //! @}

      //! @name geometrically nonlinear TSI analysis
      //! @{

      //! calculate element vectors (internal/external) and a few matrices
      //! considering current displacement solution
      //! --> all terms are coupled to the displacements/velocities
      void NonlinearThermoDispContribution(
          DRT::Element* ele,          //!< the element whose matrix is calculated
          const double& time,         //!< current time
          std::vector<double>& disp,  //!< current displacements
          std::vector<double>& vel,   //!< current velocities
          CORE::LINALG::Matrix<nen_ * numdofpernode_, nen_ * numdofpernode_>*
              econd,  //!< conductivity matrix
          CORE::LINALG::Matrix<nen_ * numdofpernode_, nen_ * numdofpernode_>*
              ecapa,  //!< capacity matrix
          CORE::LINALG::Matrix<nen_ * numdofpernode_, nen_ * numdofpernode_>*
              ecapalin,  //!< partial linearization dC/dT of capacity matrix
          CORE::LINALG::Matrix<nen_ * numdofpernode_, 1>* efint,  //!< internal force
          Teuchos::ParameterList& params);

      //! calculate thermal-mechanical system matrix k_Td needed in monolithic TSI
      void NonlinearCoupledTang(DRT::Element* ele,  //!< current element whose terms are calculated
          std::vector<double>& disp,                //!< current displacements
          std::vector<double>& vel,                 //!< current velocities
          CORE::LINALG::Matrix<nen_ * numdofpernode_, nsd_ * nen_ * numdofpernode_>*
              etangcoupl,                 //!< k_Tu matrix
          Teuchos::ParameterList& params  //!< parameter list, containing e.g., dt,theta
      );

      //! build nonlinear B-operator
      void CalculateBop(
          CORE::LINALG::Matrix<6, nsd_ * nen_ * numdofpernode_>* bop,  //!< nonlinear B-operator
          CORE::LINALG::Matrix<nsd_, nsd_>* defgrd,                    //!< deformation gradient
          CORE::LINALG::Matrix<nsd_, nen_>* N_XYZ                      //!< gradient-operator
      );

      //! build linearisation of Jacobian w.r.t. d: dJ_dd
      void CalculateLinearisationOfJacobian(
          CORE::LINALG::Matrix<1, nsd_ * nen_ * numdofpernode_>& dJ_dd,  //!<  [out] dJ_dd
          const double& J,                                               //!< Jacobian
          const CORE::LINALG::Matrix<nsd_, nen_>& N_XYZ,  //!< linear gradient of shape functions
          const CORE::LINALG::Matrix<nsd_, nsd_>& defgrd_inv  //!< inverse of F
      );

      //! build derivatives of right Cauchy-Green deformation tensor C
      //! build the inverse of C^{-1} and the time derivative C'
      void CalculateCauchyGreens(
          CORE::LINALG::Matrix<6, 1>& Cratevct,          //!< right Cauchy-Green rate vector
          CORE::LINALG::Matrix<6, 1>& invCvct,           //!< inverse of right Cauchy-Green vector
          CORE::LINALG::Matrix<nsd_, nsd_>& invC,        //!< inverse right Cauchy-Green tensor
          CORE::LINALG::Matrix<nsd_, nsd_>* defgrd,      //!< deformation gradient tensor
          CORE::LINALG::Matrix<nsd_, nsd_>* defgrdrate,  //!< velocity gradient tensor
          CORE::LINALG::Matrix<nsd_, nsd_>* invdefgrd    //!< inverse deformation gradient tensor
      );

      /// @}

      //! @name finite strain thermoplasticity solved with TSI
      //! @{

      //! calculate internal dissipation arising when a thermo-elasto-plastic
      //! material is used within geometrically nonlinear analysis
      //! Clausius-Duhem inequality is no longer = 0, but > 0:
      //! mechanical energy dissipates as heat
      void NonlinearDissipationFintTang(
          DRT::Element* ele,          //!< the element whose matrix is calculated
          std::vector<double>& disp,  //!< current displacements
          CORE::LINALG::Matrix<nen_ * numdofpernode_, nen_ * numdofpernode_>*
              econd,                                              //!< conductivity matrix
          CORE::LINALG::Matrix<nen_ * numdofpernode_, 1>* efint,  //!< internal force
          Teuchos::ParameterList& params);

      //! calculate terms of dissipation for thermo-mechanical system matrix k_Td
      //! used in case of plastic material within geometrically nonlinear analysis
      void NonlinearDissipationCoupledTang(
          DRT::Element* ele,          //!< the element whose matrix is calculated
          std::vector<double>& disp,  //!< current displacements
          std::vector<double>& vel,   //!< current velocities
          CORE::LINALG::Matrix<nen_ * numdofpernode_, nsd_ * nen_ * numdofpernode_>*
              etangcoupl,  //!< k_Td
          Teuchos::ParameterList& params);

      /// @}

      //! get the body force
      virtual void Radiation(DRT::Element* ele,  //!< current element we are dealing with
          double time                            //!< current times
      );

      //! build linear B-operator
      void CalculateBoplin(
          CORE::LINALG::Matrix<6, nsd_ * nen_ * numdofpernode_>* boplin,  //!< linear B-operator
          CORE::LINALG::Matrix<nsd_, nen_>* N_XYZ                         //!< gradient-operator
      );

      //! get corresponding structural material
      Teuchos::RCP<MAT::Material> GetSTRMaterial(
          DRT::Element* ele  //!< the element whose matrix is calculated
      );

      //! calculate reactive term
      void CalculateReactiveTerm(
          CORE::LINALG::Matrix<6, 1>* ctemp,     //!< temperature-dependent material tangent
          CORE::LINALG::Matrix<6, 1>* strainvel  //!< strain rate
      );

      //! determine heat flux and conductivity tensor
      //! based on material law
      virtual void Materialize(const DRT::Element* ele,  //!< current element
          int gp                                         //!< current GP id
      );

      //! evaluate shape functions and their derivatives at current integration point
      virtual void EvalShapeFuncAndDerivsAtIntPoint(
          const CORE::DRT::UTILS::IntPointsAndWeights<nsd_>& intpoints,  //!< integration points
          int iquad,  //!< id of current Gauss point
          int eleid   //!< the element id
      );

      //! compute heatflux and temperature gradient in linear case
      void LinearHeatfluxTempgrad(Element* ele,           //!< the current element
          CORE::LINALG::Matrix<nquad_, nsd_>* eheatflux,  //!< [out] heat fluxes at Gauss points
          CORE::LINALG::Matrix<nquad_, nsd_>*
              etempgrad  //!< [out] temperature gradients at Gauss points
      );

      //! compute heatflux and temperature gradient in nonlinear case
      void NonlinearHeatfluxTempgrad(Element* ele,        //!< the current element
          std::vector<double>& disp,                      //!< element displacements
          std::vector<double>& vel,                       //!< element velocities
          CORE::LINALG::Matrix<nquad_, nsd_>* eheatflux,  //!< [out] heat fluxes at Gauss points
          CORE::LINALG::Matrix<nquad_, nsd_>*
              etempgrad,                  //!< [out] temperature gradients at Gauss points
          Teuchos::ParameterList& params  //!< parameters containing type of flux and grad
      );

      //! calculate lumped capacity matrix in case of explicit time integration
      void CalculateLumpMatrix(
          CORE::LINALG::Matrix<nen_ * numdofpernode_, nen_ * numdofpernode_>* ecapa);

      //! calculate characteristic element length
      double CalculateCharEleLength();

      //! Compute the error compared to an analytical solution from input file
      virtual void ComputeError(DRT::Element* ele,  //!< current element whose terms are calculated
          CORE::LINALG::Matrix<nen_ * numdofpernode_, 1>& evector,  //!< element vectorr
          Teuchos::ParameterList& params  //!< parameter list, containing analytical solution
      );

      //! Compute nodal position and velocity
      inline void InitialAndCurrentNodalPositionVelocity(const Element* ele,
          const std::vector<double>& disp, const std::vector<double>& vel,
          CORE::LINALG::Matrix<nen_, nsd_>& xcurr, CORE::LINALG::Matrix<nen_, nsd_>& xcurrrate);

      //! prepare the evaluation of NURBS shape functions
      virtual void PrepareNurbsEval(DRT::Element* ele,  //!< the element whose matrix is calculated
          DRT::Discretization& discretization           //!< current discretisation
      );

      //! integral of shape functions over the element
      void IntegrateShapeFunctions(const DRT::Element* ele,  //!< current element
          CORE::LINALG::SerialDenseVector& elevec1,          //!< result vector (to be assembled)
          const CORE::LINALG::IntSerialDenseVector& dofids  //!< for which dof we need to integrate?
      );

      //! extrapolate from Gauss points to nodes, needed for postprocessing
      void ExtrapolateFromGaussPointsToNodes(DRT::Element* ele,  //!< the actual element
          const CORE::LINALG::Matrix<nquad_, nsd_>& gpheatflux,  //!< heat flux at each Gauss Point
          CORE::LINALG::Matrix<nen_ * numdofpernode_, 1>&
              efluxx,  //!< element heat flux in x-direction
          CORE::LINALG::Matrix<nen_ * numdofpernode_, 1>&
              efluxy,  //!< element heat flux in y-direction
          CORE::LINALG::Matrix<nen_ * numdofpernode_, 1>&
              efluxz  //!< element heat flux in z-direction
      );

      //! extract displacement and velocity vector from discretization
      void ExtractDispVel(const Discretization& discretization, Element::LocationArray& la,
          std::vector<double>& mydisp, std::vector<double>& myvel) const;

      //! copy matrix contents into character vector
      void CopyMatrixIntoCharVector(
          std::vector<char>& data, CORE::LINALG::Matrix<nquad_, nsd_>& stuff);

      //! FDcheck of conductivity matrix on element level
      void FDCheckCouplNlnFintCondCapa(
          DRT::Element* ele,          //!< the element whose matrix is calculated
          const double& time,         //!< current time
          std::vector<double>& disp,  //!< current displacements
          std::vector<double>& vel,   //!< current velocities
          CORE::LINALG::Matrix<nen_ * numdofpernode_, nen_ * numdofpernode_>*
              etang,                                              //!< tangent conductivity matrix
          CORE::LINALG::Matrix<nen_ * numdofpernode_, 1>* efint,  //!< internal force);)
          Teuchos::ParameterList& params);

      //! FDcheck of linearized capacity matrix on element level
      void FDCheckCapalin(DRT::Element* ele,  //!< the element whose matrix is calculated
          const double& time,                 //!< current time
          std::vector<double>& disp,          //!< current displacements
          std::vector<double>& vel,           //!< current velocities
          CORE::LINALG::Matrix<nen_ * numdofpernode_, nen_ * numdofpernode_>*
              ecapa,  //!< capacity matrix
          CORE::LINALG::Matrix<nen_ * numdofpernode_, nen_ * numdofpernode_>*
              ecapalin,  //!< linearization term from capacity matrix
          Teuchos::ParameterList& params);

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
#endif  // THERMO_ELE_IMPL_H
