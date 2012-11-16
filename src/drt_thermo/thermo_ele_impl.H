/*----------------------------------------------------------------------*/
/*!
\file thermo_ele_impl.H

\brief Internal implementation of thermo elements

<pre>
Maintainer: Caroline Danowski
            danowski@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15253
</pre>
*/

/*----------------------------------------------------------------------*
 |  definitions                                                gjb 01/08|
 *----------------------------------------------------------------------*/
#ifdef D_THERMO

#ifndef THERMO_ELE_IMPL_H
#define THERMO_ELE_IMPL_H

/*----------------------------------------------------------------------*
 |  headers                                                    gjb 01/08|
 *----------------------------------------------------------------------*/
#include "thermo_ele_impl_utils.H"
#include "../drt_fem_general/drt_utils_local_connectivity_matrices.H"
#include <Epetra_IntSerialDenseVector.h>

/*----------------------------------------------------------------------*
 |                                                             gjb 01/08|
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
    TemperImplInterface() {}
    //! Empty destructor
    virtual ~TemperImplInterface() {}

    //! Evaluate the element
    //!  This class does not provide a definition for this function, it
    //!  must be defined in TemperImpl.
    virtual int Evaluate(
      DRT::Element* ele,  //!< current element
      Teuchos::ParameterList& params,  //!< parameter list, containing e.g., dt
      DRT::Discretization& discretization,  //!< current discretisation
      DRT::Element::LocationArray& la,  //!< location array
      Epetra_SerialDenseMatrix& elemat1_epetra,  //!< conductivity matrix
      Epetra_SerialDenseMatrix& elemat2_epetra,  //!< capacity matrix
      Epetra_SerialDenseVector& elevec1_epetra,  //!< internal force, view on heat flux in x-direction
      Epetra_SerialDenseVector& elevec2_epetra,  //!< external force, view on heat flux in y-direction
      Epetra_SerialDenseVector& elevec3_epetra  //!< view on heat flux in z-direction
      ) = 0;

    //! Evaluate the element
    //!
    //!  This class does not provide a definition for this function, it
    //!  must be defined in TemperImpl.
    virtual int EvaluateNeumann(
      DRT::Element* ele,  //!< current element
      Teuchos::ParameterList& params,  //!< parameter list
      DRT::Discretization&  discretization,  //!< current discretisation
      std::vector<int>& lm,  //!< location vector, EvalNeumann is called only on own discretisation
      Epetra_SerialDenseVector& elevec1_epetra,  //!< view on external force vector
      Epetra_SerialDenseMatrix* elemat1_epetra  //!< matrix is not needed
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
  template<DRT::Element::DiscretizationType distype>
  class TemperImpl: public TemperImplInterface
  {
  public:

    //! Constructor
    TemperImpl();

    //! Singleton access method
    static TemperImpl<distype> * Instance(bool create=true);

    //! called upon destruction
    virtual void Done();

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
    virtual int Evaluate(
      DRT::Element* ele,  //!< current element
      Teuchos::ParameterList& params,  //!< parameter list, containing e.g., dt
      DRT::Discretization& discretization,  //!< current discretisation
      DRT::Element::LocationArray& la,  //!< location array
      Epetra_SerialDenseMatrix& elemat1_epetra,  //!< conductivity matrix
      Epetra_SerialDenseMatrix& elemat2_epetra,  //!< capacity matrix
      Epetra_SerialDenseVector& elevec1_epetra,  //!< internal force, view on heat flux in x-direction
      Epetra_SerialDenseVector& elevec2_epetra,  //!< external force, view on heat flux in y-direction
      Epetra_SerialDenseVector& elevec3_epetra  //!< view on heat flux in z-direction
      );

    //! Evaluate the element
    virtual int EvaluateNeumann(
      DRT::Element* ele,  //!< current element
      Teuchos::ParameterList& params,  //!< parameter list
      DRT::Discretization&  discretization,  //!< current discretisation
      std::vector<int>& lm,  //!< location vector, EvalNeumann is called only on own discretisation
      Epetra_SerialDenseVector& elevec1_epetra,  //!< view on external force vector
      Epetra_SerialDenseMatrix* elemat1_epetra  //!< matrix is not needed
      );

  private:

    //! @name general thermal functions

    //! Calculate element vectors (internal/external) and a few matrices
    void CalculateFintCondCapa(
      DRT::Element* ele,  //!< the element whose matrix is calculated
      const double& time,  //!< current time
      LINALG::Matrix<nen_*numdofpernode_,nen_*numdofpernode_>* econd,  //!< conductivity matrix
      LINALG::Matrix<nen_*numdofpernode_,nen_*numdofpernode_>* ecapa,  //!< capacity matrix
      LINALG::Matrix<nen_*numdofpernode_,1>* efint,  //!< internal force
      LINALG::Matrix<nen_*numdofpernode_,1>* efext,  //!< external force
      LINALG::Matrix<nquad_,nsd_>* eheatflux,  //!< heat fluxes at Gauss points
      LINALG::Matrix<nquad_,nsd_>* etempgrad  //!< temperature gradients at Gauss points
      );

    //! get the body force
    virtual void Radiation(
      const DRT::Element* ele,  //!< current element we are dealing with
      const double time  //!< current times
      );

    //! determine heat flux and conductivity tensor
    //! based on material law
    virtual void Materialize(
      const DRT::Element* ele  //!< current element
      );

    //! evaluate shape functions and their derivatives at current integration point
    virtual void EvalShapeFuncAndDerivsAtIntPoint(
      const DRT::UTILS::IntPointsAndWeights<nsd_>& intpoints,  //!< integration points
      const int iquad,  //!< id of current Gauss point
      const int eleid  //!< the element id
      );

    //! integral of shape functions over the element
    void IntegrateShapeFunctions(
      const DRT::Element* ele,  //!< current element
      Epetra_SerialDenseVector& elevec1,  //!< result vector (to be assembled)
      const Epetra_IntSerialDenseVector& dofids  //!< for which dof we need to integrate?
      );

    //! extrapolate from Gauss points to nodes, needed for postprocessing
    void ExtrapolateFromGaussPointsToNodes(
      DRT::Element* ele,  //!< the actual element
      const LINALG::Matrix<nquad_,nsd_>& gpheatflux,  //!< heat flux at each Gauss Point
      LINALG::Matrix<nen_*numdofpernode_,1>& efluxx,  //!< element heat flux in x-direction
      LINALG::Matrix<nen_*numdofpernode_,1>& efluxy,  //!< element heat flux in y-direction
      LINALG::Matrix<nen_*numdofpernode_,1>& efluxz  //!< element heat flux in z-direction
      );

    //! calculate lumped capacity matrix in case of explicit time integration
    void CalculateLumpMatrix(LINALG::Matrix<nen_*numdofpernode_,nen_*numdofpernode_>* ecapa);

    //! calculate characteristic element length
    double CalculateCharEleLength();

    //@}

    //! @name TSI related stuff

    //! Calculate element vectors (internal/external) and a few matrices
    //! considering current displacement solution
    void CalculateCouplFintCondCapa(
      DRT::Element* ele,  //!< the element whose matrix is calculated
      const double& time,  //!< current time
      std::vector<double>& disp,  //!< current displacements
      std::vector<double>& vel,  //!< current velocities
      const double& stepsize,  //!< time step size
      LINALG::Matrix<nen_*numdofpernode_,nen_*numdofpernode_>* econd,  //!< conductivity matrix
      LINALG::Matrix<nen_*numdofpernode_,nen_*numdofpernode_>* ecapa,  //!< capacity matrix
      LINALG::Matrix<nen_*numdofpernode_,1>* efint,  //!< internal force
      LINALG::Matrix<nen_*numdofpernode_,1>* efext,  //!< external force
      LINALG::Matrix<nquad_,nsd_>* eheatflux,  //!< heat fluxes at Gauss points
      LINALG::Matrix<nquad_,nsd_>* etempgrad  //!< temperature gradients at Gauss points
      );

    //! calculate internal dissipation arising when a thermo-elasto-plastic
    //! material is used
    //! Clausius-Duhem inequality is no longer = 0, but > 0:
    //! mechanical energy dissipates as heat
    void CalculateInternalDissipation(
      DRT::Element* ele,  //!< the element whose matrix is calculated
      std::vector<double>& vel,  //!< current velocities
      const double& stepsize,  //!< step size Delta t
      LINALG::Matrix<nen_*numdofpernode_,nen_*numdofpernode_>* etang,  //!< conductivity matrix
      LINALG::Matrix<nen_*numdofpernode_,1>* efint  //!< internal force
      );

    //! calculate thermal-mechanical system matrix term needed in monolithic TSI
    void CalculateCouplCond(
      DRT::Element* ele,  //!< the element whose matrix is calculated
      std::vector<double>& disp,  //!< current displacements
      std::vector<double>& vel,  //!< current velocities
      LINALG::Matrix<nen_*numdofpernode_,nsd_*nen_*numdofpernode_>* econdcoupl  //!< k_Tu matrix
      );

    //@}

    //! actual values of temperatures
    LINALG::Matrix<nen_,1> etemp_;
    //! capacity matrix
    LINALG::Matrix<nen_,nen_> ecapa_;

    //! node coordinates
    LINALG::Matrix<nsd_,nen_> xyze_;
    //! radiation in element nodes
    LINALG::Matrix<numdofpernode_,1> radiation_;
    //! coordinates of current integration point in reference coordinates
    LINALG::Matrix<nsd_,1> xsi_;
    //! array for shape functions
    LINALG::Matrix<nen_,1> funct_;
    //! array for shape function derivatives w.r.t r,s,t
    LINALG::Matrix<nsd_,nen_> deriv_;
    //! transposed jacobian "dx/ds"
    LINALG::Matrix<nsd_,nsd_> xjm_;
    //! inverse of transposed jacobian "ds/dx"
    LINALG::Matrix<nsd_,nsd_> xij_;
    //! global derivatives of shape functions w.r.t x,y,z
    LINALG::Matrix<nsd_,nen_> derxy_;
    //! global temperature rate derivatives in gausspoint w.r.t x,y,z
    LINALG::Matrix<nsd_,nsd_> rderxy_;
    //! integration factor for current GP: fac = GaussWeight * det(J)
    double fac_;
    //! (global) gradient of temperature at integration point
    LINALG::Matrix<nsd_,1> gradtemp_;
    //! (global) gradient of temperature at integration point
    LINALG::Matrix<nsd_,1> heatflux_;
    //! (global) conductivity 2-tensor
    LINALG::Matrix<nsd_,nsd_> cmat_;
    //! capacity density
    double capacoeff_;

    //! @name material related stuff

    //! material tangent of the structure
    LINALG::Matrix<6,1> ctemp_;

    //! flag plastic material is used
    bool plasticmat_;

    //@}

  };  // class TemperImpl

}  // namespace ELEMENTS

}  // namespace DRT


/*----------------------------------------------------------------------*/
#endif // THERMO_ELE_IMPL_H
#endif // D_THERMO
