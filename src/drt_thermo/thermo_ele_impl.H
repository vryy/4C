/*----------------------------------------------------------------------*/
/*! \file

\brief Internal implementation of thermo elements

\level 1

\maintainer Sebastian Proell

*/

/*----------------------------------------------------------------------*
 | definitions                                                gjb 01/08 |
 *----------------------------------------------------------------------*/
#ifdef D_THERMO

#ifndef THERMO_ELE_IMPL_H
#define THERMO_ELE_IMPL_H

/*----------------------------------------------------------------------*
 | headers                                                    gjb 01/08 |
 *----------------------------------------------------------------------*/
#include "thermo_ele_impl_utils.H"
#include "../drt_inpar/inpar_thermo.H"

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
    class TemperImplInterface : SingletonDestruction
    {
     public:
      //! Empty constructor
      TemperImplInterface() = default;
      //! Empty destructor
      ~TemperImplInterface() override = default;

      //! Evaluate the element
      //!  This class does not provide a definition for this function, it
      //!  must be defined in TemperImpl.
      virtual int Evaluate(DRT::Element* ele,        //!< current element
          Teuchos::ParameterList& params,            //!< parameter list, containing e.g., dt
          DRT::Discretization& discretization,       //!< current discretisation
          DRT::Element::LocationArray& la,           //!< location array
          Epetra_SerialDenseMatrix& elemat1_epetra,  //!< conductivity matrix
          Epetra_SerialDenseMatrix& elemat2_epetra,  //!< capacity matrix
          Epetra_SerialDenseVector&
              elevec1_epetra,  //!< internal force, view on heat flux in x-direction
          Epetra_SerialDenseVector&
              elevec2_epetra,  //!< external force, view on heat flux in y-direction
          Epetra_SerialDenseVector& elevec3_epetra  //!< view on heat flux in z-direction
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
          Epetra_SerialDenseVector& elevec1_epetra,  //!< view on external force vector
          Epetra_SerialDenseMatrix* elemat1_epetra   //!< matrix is not needed
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
      static TemperImpl<distype>* Instance(bool create = true);

      //! called upon destruction
      void Done() override;

      //! number of nodes
      static const int nen_ = DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;

      //! number of space dimensions
      static const int nsd_ = DRT::UTILS::DisTypeToDim<distype>::dim;

      //! number of dof per node
      static const int numdofpernode_ = 1;

      //! number of Gauss points
      static const int nquad_ = THR::DisTypeToNumGaussPoints<distype>::nquad;

      // new method for volume coupling
      //! Evaluate for multiple dofsets
      int Evaluate(DRT::Element* ele,                //!< current element
          Teuchos::ParameterList& params,            //!< parameter list, containing e.g., dt
          DRT::Discretization& discretization,       //!< current discretisation
          DRT::Element::LocationArray& la,           //!< location array
          Epetra_SerialDenseMatrix& elemat1_epetra,  //!< conductivity matrix
          Epetra_SerialDenseMatrix& elemat2_epetra,  //!< capacity matrix
          Epetra_SerialDenseVector&
              elevec1_epetra,  //!< internal force, view on heat flux in x-direction
          Epetra_SerialDenseVector&
              elevec2_epetra,  //!< external force, view on heat flux in y-direction
          Epetra_SerialDenseVector& elevec3_epetra  //!< view on heat flux in z-direction
          ) override;

      //! Evaluate the element
      int EvaluateNeumann(DRT::Element* ele,    //!< current element
          Teuchos::ParameterList& params,       //!< parameter list
          DRT::Discretization& discretization,  //!< current discretisation
          std::vector<int>&
              lm,  //!< location vector, EvalNeumann is called only on own discretisation
          Epetra_SerialDenseVector& elevec1_epetra,  //!< view on external force vector
          Epetra_SerialDenseMatrix* elemat1_epetra   //!< matrix is not needed
          ) override;

     private:
      //! @name general thermal functions

      //! calculate complete internal force, tangent matrix k_TT and capacity matrix
      //!
      //! builds quantitites from linear/nonlinear and thermo/TSI specific routines
      void EvaluateFintTang(Element* ele, const double& time, Discretization& discretization,
          Element::LocationArray& la,
          LINALG::Matrix<nen_ * numdofpernode_, nen_ * numdofpernode_>* etang,
          LINALG::Matrix<nen_ * numdofpernode_, nen_ * numdofpernode_>* ecapa,
          LINALG::Matrix<nen_ * numdofpernode_, 1>* efint, Teuchos::ParameterList& params);

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
          LINALG::Matrix<nen_ * numdofpernode_, nen_ * nsd_ * numdofpernode_>* etangcoupl,
          Teuchos::ParameterList& params);

      /*!
       * evaluate external body loads
       * @param ele the element whose matrix is calculated
       * @param time time for function evaluation
       * @param efext external force vector
       */
      void EvaluateFext(
          Element* ele, const double& time, LINALG::Matrix<nen_ * numdofpernode_, 1>& efext);

      //! Calculate element force vectors and a few matrices
      void LinearThermoContribution(DRT::Element* ele,  //!< the element whose matrix is calculated
          const double& time,                           //!< current time
          LINALG::Matrix<nen_ * numdofpernode_, nen_ * numdofpernode_>*
              econd,  //!< conductivity matrix
          LINALG::Matrix<nen_ * numdofpernode_, nen_ * numdofpernode_>* ecapa,  //!< capacity matrix
          LINALG::Matrix<nen_ * numdofpernode_, nen_ * numdofpernode_>*
              ecapalin,  //!< linearization contribution of capacity
          LINALG::Matrix<nen_ * numdofpernode_, 1>* efint  //!< internal force
      );

      //! get the body force
      virtual void Radiation(DRT::Element* ele,  //!< current element we are dealing with
          double time                            //!< current times
      );

      //! post-Newton iteration update of temperatures and source terms to account for latent heat
      //! of phase change
      virtual void CalculatePhaseChangeIncrement(
          DRT::Element* ele,  //!< the element whose matrix is calculated
          LINALG::Matrix<nen_ * numdofpernode_, 1>*
              esourceinc,  //!< artificial increment in source term
          LINALG::Matrix<nen_ * numdofpernode_, 1>* etempi, Teuchos::ParameterList& params);

      virtual void CalculateTotalLatentHeat(DRT::Element* ele, Teuchos::ParameterList& params);

      //! determine heat flux and conductivity tensor
      //! based on material law
      virtual void Materialize(const DRT::Element* ele,  //!< current element
          int gp                                         //!< current GP id
      );

      //! evaluate shape functions and their derivatives at current integration point
      virtual void EvalShapeFuncAndDerivsAtIntPoint(
          const DRT::UTILS::IntPointsAndWeights<nsd_>& intpoints,  //!< integration points
          int iquad,                                               //!< id of current Gauss point
          int eleid                                                //!< the element id
      );

      //! Compute nodal position and velocity
      inline void InitialAndCurrentNodalPositionVelocity(const Element* ele,
          const std::vector<double>& disp, const std::vector<double>& vel,
          LINALG::Matrix<nen_, nsd_>& xcurr, LINALG::Matrix<nen_, nsd_>& xcurrrate);

      //! prepare the evaluation of NURBS shape functions
      virtual void PrepareNurbsEval(DRT::Element* ele,  //!< the element whose matrix is calculated
          DRT::Discretization& discretization           //!< current discretisation
      );

      //! integral of shape functions over the element
      void IntegrateShapeFunctions(const DRT::Element* ele,  //!< current element
          Epetra_SerialDenseVector& elevec1,                 //!< result vector (to be assembled)
          const Epetra_IntSerialDenseVector& dofids  //!< for which dof we need to integrate?
      );

      //! extrapolate from Gauss points to nodes, needed for postprocessing
      void ExtrapolateFromGaussPointsToNodes(DRT::Element* ele,  //!< the actual element
          const LINALG::Matrix<nquad_, nsd_>& gpheatflux,        //!< heat flux at each Gauss Point
          LINALG::Matrix<nen_ * numdofpernode_, 1>& efluxx,  //!< element heat flux in x-direction
          LINALG::Matrix<nen_ * numdofpernode_, 1>& efluxy,  //!< element heat flux in y-direction
          LINALG::Matrix<nen_ * numdofpernode_, 1>& efluxz   //!< element heat flux in z-direction
      );

      //! compute heatflux and temperature gradient in linear case
      void LinearHeatfluxTempgrad(Element* ele,     //!< the current element
          LINALG::Matrix<nquad_, nsd_>* eheatflux,  //!< [out] heat fluxes at Gauss points
          LINALG::Matrix<nquad_, nsd_>* etempgrad   //!< [out] temperature gradients at Gauss points
      );

      //! compute heatflux and temperature gradient in nonlinear case
      void NonlinearHeatfluxTempgrad(Element* ele,  //!< the current element
          std::vector<double>& disp,                //!< element displacements
          std::vector<double>& vel,                 //!< element velocities
          LINALG::Matrix<nquad_, nsd_>* eheatflux,  //!< [out] heat fluxes at Gauss points
          LINALG::Matrix<nquad_, nsd_>* etempgrad,  //!< [out] temperature gradients at Gauss points
          Teuchos::ParameterList& params            //!< parameters containing type of flux and grad
      );

      //! calculate lumped capacity matrix in case of explicit time integration
      void CalculateLumpMatrix(LINALG::Matrix<nen_ * numdofpernode_, nen_ * numdofpernode_>* ecapa);

      //! calculate characteristic element length
      double CalculateCharEleLength();

      /// Evaluate element ERROR
      /*!
          general function to compute the error (analytical solution from input file)
       */
      virtual void ComputeError(DRT::Element* ele,  //!< current element whose terms are calculated
          LINALG::Matrix<nen_ * numdofpernode_, 1>& evector,  //!< element vectorr
          Teuchos::ParameterList& params  //!< parameter list, containing e.g., dt,theta
      );
      //! @}

      //! @name geometrically linear TSI
      //! @{

      //! Calculate element vectors (internal/external) and a few matrices
      //! considering current displacement solution
      void LinearDispContribution(DRT::Element* ele, const double& time, std::vector<double>& disp,
          std::vector<double>& vel,
          LINALG::Matrix<nen_ * numdofpernode_, nen_ * numdofpernode_>* econd,
          LINALG::Matrix<nen_ * numdofpernode_, 1>* efint, Teuchos::ParameterList& params);

      //! calculate thermal-mechanical system matrix term needed in monolithic TSI
      void LinearCoupledTang(DRT::Element* ele,  //!< the element whose matrix is calculated
          std::vector<double>& disp,             //!< current displacements
          std::vector<double>& vel,              //!< current velocities
          LINALG::Matrix<nen_ * numdofpernode_, nsd_ * nen_ * numdofpernode_>*
              etangcoupl,  //!< k_Tu matrix
          Teuchos::ParameterList& params);

      //! build linear B-operator
      void CalculateBoplin(
          LINALG::Matrix<6, nsd_ * nen_ * numdofpernode_>* boplin,  //!< linear B-operator
          LINALG::Matrix<nsd_, nen_>* N_XYZ                         //!< gradient-operator
      );

      //! get corresponding structural material
      Teuchos::RCP<MAT::Material> GetSTRMaterial(
          DRT::Element* ele  //!< the element whose matrix is calculated
      );

      //! calculate reactive term
      void CalculateReactiveTerm(
          LINALG::Matrix<6, 1>* ctemp,     //!< temperature-dependent material tangent
          LINALG::Matrix<6, 1>* strainvel  //!< strain rate
      );

      //! @}

      //! @name linear, small strain thermoplasticity solved with TSI
      //! @{

      //! calculate internal dissipation arising when a thermo-elasto-plastic
      //! material is used
      //! Clausius-Duhem inequality is no longer = 0, but > 0:
      //! mechanical energy dissipates as heat
      void LinearDissipationFint(DRT::Element* ele,  //!< the element whose matrix is calculated
          LINALG::Matrix<nen_ * numdofpernode_, 1>* efint,  //!< internal force
          Teuchos::ParameterList& params);

      //! calculate terms of dissipation for thermo-mechanical
      //! system matrix k_Td used in case of plastic material
      void LinearDissipationCoupledTang(
          DRT::Element* ele,  // the element whose matrix is calculated
          LINALG::Matrix<nen_ * numdofpernode_, nsd_ * nen_ * numdofpernode_>* etangcoupl,  // k_Td
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
          LINALG::Matrix<nen_ * numdofpernode_, nen_ * numdofpernode_>*
              econd,  //!< conductivity matrix
          LINALG::Matrix<nen_ * numdofpernode_, nen_ * numdofpernode_>* ecapa,  //!< capacity matrix
          LINALG::Matrix<nen_ * numdofpernode_, nen_ * numdofpernode_>*
              ecapalin,  //!< partial linearization dC/dT of capacity matrix
          LINALG::Matrix<nen_ * numdofpernode_, 1>* efint,  //!< internal force
          Teuchos::ParameterList& params);

      //! calculate thermal-mechanical system matrix k_Td needed in monolithic TSI
      void NonlinearCoupledTang(DRT::Element* ele,  //!< current element whose terms are calculated
          std::vector<double>& disp,                //!< current displacements
          std::vector<double>& vel,                 //!< current velocities
          LINALG::Matrix<nen_ * numdofpernode_, nsd_ * nen_ * numdofpernode_>*
              etangcoupl,                 //!< k_Tu matrix
          Teuchos::ParameterList& params  //!< parameter list, containing e.g., dt,theta
      );

      //! build nonlinear B-operator
      void CalculateBop(
          LINALG::Matrix<6, nsd_ * nen_ * numdofpernode_>* bop,  //!< nonlinear B-operator
          LINALG::Matrix<nsd_, nsd_>* defgrd,                    //!< deformation gradient
          LINALG::Matrix<nsd_, nen_>* N_XYZ                      //!< gradient-operator
      );

      //! build linearisation of Jacobian w.r.t. d: dJ_dd
      void CalculateLinearisationOfJacobian(
          LINALG::Matrix<1, nsd_ * nen_ * numdofpernode_>& dJ_dd,  //  (o) dJ_dd
          const double& J,                                         // Jacobian
          const LINALG::Matrix<nsd_, nen_>& N_XYZ,      // linear gradient of shape functions
          const LINALG::Matrix<nsd_, nsd_>& defgrd_inv  // inverse of F
      );

      //! build derivatives of right Cauchy-Green deformation tensor C
      //! build the inverse of C^{-1} and the time derivative C'
      void CalculateCauchyGreens(
          LINALG::Matrix<6, 1>& Cratevct,          //!< right Cauchy-Green rate vector
          LINALG::Matrix<6, 1>& invCvct,           //!< inverse of right Cauchy-Green vector
          LINALG::Matrix<nsd_, nsd_>& invC,        //!< inverse right Cauchy-Green tensor
          LINALG::Matrix<nsd_, nsd_>* defgrd,      //!< deformation gradient tensor
          LINALG::Matrix<nsd_, nsd_>* defgrdrate,  //!< velocity gradient tensor
          LINALG::Matrix<nsd_, nsd_>* invdefgrd    //!< inverse deformation gradient tensor
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
          LINALG::Matrix<nen_ * numdofpernode_, nen_ * numdofpernode_>*
              econd,                                        //!< conductivity matrix
          LINALG::Matrix<nen_ * numdofpernode_, 1>* efint,  //!< internal force
          Teuchos::ParameterList& params);

      //! calculate terms of dissipation for thermo-mechanical system matrix k_Td
      //! used in case of plastic material within geometrically nonlinear analysis
      void NonlinearDissipationCoupledTang(
          DRT::Element* ele,          //!< the element whose matrix is calculated
          std::vector<double>& disp,  //!< current displacements
          std::vector<double>& vel,   //!< current velocities
          LINALG::Matrix<nen_ * numdofpernode_, nsd_ * nen_ * numdofpernode_>* etangcoupl,  // k_Td
          Teuchos::ParameterList& params);

      /// @}

      //! extract displacement and velocity vector from discretization
      void ExtractDispVel(const Discretization& discretization, Element::LocationArray& la,
          std::vector<double>& mydisp, std::vector<double>& myvel) const;

      // extract the node based phase information
      void ExtractPhaseInformation(DRT::Element* ele,         ///< current element
          LINALG::Matrix<nen_ * numdofpernode_, 1>& ephase);  ///< vector to fill, view on elevec1

      //! copy matrix contents into character vector
      void CopyMatrixIntoCharVector(std::vector<char>& data, LINALG::Matrix<nquad_, nsd_>& stuff);

      // FDcheck of conductivity matrix on element level
      void FDCheckCouplNlnFintCondCapa(
          DRT::Element* ele,          //!< the element whose matrix is calculated
          const double& time,         //!< current time
          std::vector<double>& disp,  //!< current displacements
          std::vector<double>& vel,   //!< current velocities
          LINALG::Matrix<nen_ * numdofpernode_, nen_ * numdofpernode_>*
              etang,                                        //!< tangent conductivity matrix
          LINALG::Matrix<nen_ * numdofpernode_, 1>* efint,  //!< internal force);)
          Teuchos::ParameterList& params);

      // FDcheck of linearized capacity matrix on element level
      void FDCheckCapalin(DRT::Element* ele,  //!< the element whose matrix is calculated
          const double& time,                 //!< current time
          std::vector<double>& disp,          //!< current displacements
          std::vector<double>& vel,           //!< current velocities
          LINALG::Matrix<nen_ * numdofpernode_, nen_ * numdofpernode_>* ecapa,  //!< capacity matrix
          LINALG::Matrix<nen_ * numdofpernode_, nen_ * numdofpernode_>*
              ecapalin,  //!< linearization term from capacity matrix
          Teuchos::ParameterList& params);

      //! actual values of temperatures T_{n+1}
      LINALG::Matrix<nen_, 1> etempn_;
      //! temperatures in last time step T_{n}
      LINALG::Matrix<nen_, 1> etemp_;
      //! capacity matrix
      LINALG::Matrix<nen_, nen_> ecapa_;

      //! node reference coordinates
      LINALG::Matrix<nsd_, nen_> xyze_;
      //! radiation in element nodes
      LINALG::Matrix<numdofpernode_, 1> radiation_;
      //! coordinates of current integration point in reference coordinates
      LINALG::Matrix<nsd_, 1> xsi_;
      //! array for shape functions
      LINALG::Matrix<nen_, 1> funct_;
      //! array for shape function derivatives w.r.t r,s,t
      LINALG::Matrix<nsd_, nen_> deriv_;
      //! transposed jacobian "dx/ds"
      LINALG::Matrix<nsd_, nsd_> xjm_;
      //! inverse of transposed jacobian "ds/dx"
      LINALG::Matrix<nsd_, nsd_> xij_;
      //! global derivatives of shape functions w.r.t x,y,z
      LINALG::Matrix<nsd_, nen_> derxy_;
      //! integration factor for current GP: fac = GaussWeight * det(J)
      double fac_;
      //! (global) gradient of temperature at integration point
      LINALG::Matrix<nsd_, 1> gradtemp_;
      //! (global) gradient of temperature at integration point
      LINALG::Matrix<nsd_, 1> heatflux_;
      //! (global) conductivity 2-tensor
      LINALG::Matrix<nsd_, nsd_> cmat_;
      //! (global) derivative of conductivity 2-tensor w.r.t. T
      LINALG::Matrix<nsd_, nsd_> dercmat_;
      //! capacity density
      double capacoeff_;
      //! derivative of capacity w.r.t. T
      double dercapa_;

      //! @name material related stuff
      //! @{

      //! flag plastic material is used
      bool plasticmat_;

      //! nurbs specific: element knots
      std::vector<Epetra_SerialDenseVector> myknots_;
      //! nurbs specific: control point weights
      LINALG::Matrix<nen_, 1> weights_;

      /// @}

    };  // class TemperImpl

  }  // namespace ELEMENTS

}  // namespace DRT


/*----------------------------------------------------------------------*/
#endif  // THERMO_ELE_IMPL_H
#endif  // D_THERMO
