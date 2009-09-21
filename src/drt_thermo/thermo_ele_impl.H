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
#ifdef CCADISCRET
#ifdef D_THERMO

#ifndef THERMO_ELE_IMPL_H
#define THERMO_ELE_IMPL_H

/*----------------------------------------------------------------------*
 |  headers                                                    gjb 01/08|
 *----------------------------------------------------------------------*/
#include "thermo_ele_impl_utils.H"
#include "../drt_fem_general/drt_utils_local_connectivity_matrices.H"
#include <Epetra_IntSerialDenseVector.h>


namespace DRT
  {
  namespace ELEMENTS
  {

  /// Interface base class for TemperImpl
  /*!
    This class exists to provide a common interface for all template
    versions of TemperImpl. The only function
    this class actually defines is Impl, which returns a pointer to
    the appropriate version of TemperImpl.
   */
  class TemperImplInterface
  {

  public:
    /// Empty constructor
    TemperImplInterface() {}
    /// Empty destructor
    virtual ~TemperImplInterface() {}
    /// Evaluate the element
    /*!
      This class does not provide a definition for this function, it
      must be defined in TemperImpl.
     */
    virtual int Evaluate(DRT::Element*              ele,
                         Teuchos::ParameterList&    params,
                         DRT::Discretization&       discretization,
                         std::vector<int>&          lm,
                         Epetra_SerialDenseMatrix&  elemat1_epetra,
                         Epetra_SerialDenseMatrix&  elemat2_epetra,
                         Epetra_SerialDenseVector&  elevec1_epetra,
                         Epetra_SerialDenseVector&  elevec2_epetra,
                         Epetra_SerialDenseVector&  elevec3_epetra
    ) = 0;

    /// Evaluate the element
    /*!
      This class does not provide a definition for this function, it
      must be defined in TemperImpl.
     */
    virtual int EvaluateNeumann(DRT::Element*              ele,
                                Teuchos::ParameterList&    params,
                                DRT::Discretization&       discretization,
                                std::vector<int>&          lm,
                                Epetra_SerialDenseVector&  elevec1_epetra,
                                Epetra_SerialDenseMatrix*  elemat1_epetra
      ) = 0;

    /// Internal implementation class for thermo elements
    static TemperImplInterface* Impl(DRT::Element* ele);

  };

    /// Internal Thermo element implementation
    /*!
      This internal class keeps all the working arrays needed to
      calculate the thermo element. Additionally the method Sysmat()
      provides a clean and fast element implementation.

      <h3>Purpose</h3>

      The idea is to separate the element maintenance (class Thermo)
      from the mathematical contents (this class). Of course there are
      different implementations of the Thermo element, this is just one
      such implementation.

      The Thermo element will allocate exactly one object of this class
      for all thermo elements with the same number of nodes in the mesh.
      This allows us to use exactly matching working arrays (and keep them
      around.)

      The code is meant to be as clean as possible. This is the only way
      to keep it fast. The number of working arrays has to be reduced to
      a minimum so that the element fits into the cache. (There might be
      room for improvements.)

      <h3>History</h3>

      \author dano
      \date 09/09
    */
  template<DRT::Element::DiscretizationType distype>
  class TemperImpl: public TemperImplInterface
  {
  public:
    //! Constructor
    TemperImpl(int numdofpernode);

    //! number of nodes
    static const int iel = DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;

    //! number of space dimensions
    static const int nsd_ = DRT::UTILS::DisTypeToDim<distype>::dim;

    //! number of dof per node
    static const int numdofpernode_ = 1;

    //! Evaluate
    virtual int Evaluate(DRT::Element*              ele,
                         Teuchos::ParameterList&    params,
                         DRT::Discretization&       discretization,
                         std::vector<int>&          lm,
                         Epetra_SerialDenseMatrix&  elemat1_epetra,
                         Epetra_SerialDenseMatrix&  elemat2_epetra,
                         Epetra_SerialDenseVector&  elevec1_epetra,
                         Epetra_SerialDenseVector&  elevec2_epetra,
                         Epetra_SerialDenseVector&  elevec3_epetra
                         );

    //! Evaluate the element
    virtual int EvaluateNeumann(DRT::Element*              ele,
                                Teuchos::ParameterList&    params,
                                DRT::Discretization&       discretization,
                                std::vector<int>&          lm,
                                Epetra_SerialDenseVector&  elevec1_epetra,
                                Epetra_SerialDenseMatrix*  elemat1_epetra
                                );

     private:

     //! Calculate element vectors (internal/external) and a few matrices
     void CalculateFintTangCapa(
         DRT::Element* ele,  ///< the element whose matrix is calculated
         const double& time,  ///< current time
         LINALG::Matrix<iel*numdofpernode_,iel*numdofpernode_>* ekond,  ///< conductivity matrix
         LINALG::Matrix<iel*numdofpernode_,iel*numdofpernode_>* ecapa,  ///< capacity matrix
         LINALG::Matrix<iel*numdofpernode_,1>* efint,  ///< internal force
         LINALG::Matrix<iel*numdofpernode_,1>* efext  ///< external force
         );

      /*
      //! calculate mass matrix and rhs for determining initial time derivative
      virtual void InitialTimeDerivative(
          DRT::Element*                         ele,    ///< the element
          Epetra_SerialDenseMatrix&             emat,   ///< element matrix to calculate
          Epetra_SerialDenseVector&             erhs   ///< element rhs to calculate
      );
*/

      //! get the body force
      virtual void Radiation(
          const DRT::Element*    ele,  ///< the element we are dealing with
          const double time            ///< current times
      );

      //! determine heat flux and conductivity tensor
      //! based on material law
      virtual void Materialize(
          const DRT::Element* ele                ///< the element
      );

      //! evaluate shape functions and their derivatives at current integration point
      virtual void EvalShapeFuncAndDerivsAtIntPoint(
          const DRT::UTILS::IntPointsAndWeights<nsd_>& intpoints, ///< integration points
          const int                                    iquad,     ///< id of current Gauss point
          const int                                    eleid      ///< the element id
      );

      /*
      //! calculate flux vector for actual solution
      virtual void CalculateHeatFlux(
          LINALG::Matrix<3,iel>&          flux,
          const DRT::Element*             ele,
          const vector<double>&           etempnp,
          const Epetra_SerialDenseVector& erate,
      );
      */

     //! integral of shape functions over the element
     void IntegrateShapeFunctions(
         const DRT::Element*                ele,      ///< the actual element
         Epetra_SerialDenseVector&          elevec1,  ///< result vector (to be assembled)
         const Epetra_IntSerialDenseVector& dofids    ///< for which dof we need to integrate?
     );

      //! actual values of temperatures
      LINALG::Matrix<iel,1> etemp_;

      //! node coordinates
      LINALG::Matrix<nsd_,iel> xyze_;
      //! radiation in element nodes
      LINALG::Matrix<numdofpernode_,1> radiation_;
      //! coordinates of current integration point in reference coordinates
      LINALG::Matrix<nsd_,1> xsi_;
      //! array for shape functions
      LINALG::Matrix<iel,1> funct_;
      //! array for shape function derivatives w.r.t r,s,t
      LINALG::Matrix<nsd_,iel> deriv_;
      //! transposed jacobian "dx/ds"
      LINALG::Matrix<nsd_,nsd_> xjm_;
      //! inverse of transposed jacobian "ds/dx"
      LINALG::Matrix<nsd_,nsd_> xij_;
      //! global derivatives of shape functions w.r.t x,y,z
      LINALG::Matrix<nsd_,iel> derxy_;
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

    }; //TemperImpl

  } // ELEMENTS

} // DRT


#endif // THERMO_ELE_IMPL_H
#endif // D_THERMO
#endif // CCADISCRET
