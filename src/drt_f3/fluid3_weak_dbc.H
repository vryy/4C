/*----------------------------------------------------------------------*/
/*!
\file fluid3_genalpha_resVMM.H

\brief Internal implementation of Fluid3 element with a generalisbed alpha
       time integration.

<pre>
Maintainer: Peter Gamnitzer
            gamnitzer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>
*/
/*----------------------------------------------------------------------*/
#ifdef D_FLUID3
#ifdef CCADISCRET

#ifndef FLUID3_WEAKDBC_H
#define FLUID3_WEAKDBC_H

#include "fluid3.H"

#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_timecurve.H"
#include "../drt_lib/drt_function.H"

#include "../drt_fem_general/drt_utils_local_connectivity_matrices.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_fem_general/drt_utils_integration.H"
#include "../drt_fem_general/drt_utils_boundary_integration.H"

#include "../drt_mat/newtonianfluid.H"
#include "../drt_mat/carreauyasuda.H"
#include "../drt_mat/modpowerlaw.H"

namespace DRT
{
  namespace ELEMENTS
  {

    //-----------------------------------------------------------------
    ///  an interface class for the weak Dirichlet boundary condition
    ///     for every surface element distype this class allocates
    ///       one instance of the weak Dirichlet implementation
    //-----------------------------------------------------------------
   
    class Fluid3SurfaceWeakDBCInterface
    {
    public:
      /// Empty constructor
      Fluid3SurfaceWeakDBCInterface() {}
      /// Empty destructor
      virtual ~Fluid3SurfaceWeakDBCInterface() {}

      /// Evaluate the surface elements weak Dirichlet boundary 
      /// conditions
      /*!
        This class does not provide a definition for this function, it
        is defined in the implementation.
      */
      virtual int EvaluateWeakDBC(
        Fluid3Surface*             surfele       ,
        ParameterList&             params        ,
        DRT::Discretization&       discretization,
        vector<int>&               lm            ,
        Epetra_SerialDenseMatrix&  elemat_epetra ,
        Epetra_SerialDenseVector&  elevec_epetra ) = 0;

      
      /// Allocate one static instance of the internal implementation 
      /// class for weak dirichlet condition and return pointer to it
      static Fluid3SurfaceWeakDBCInterface* Impl(DRT::ELEMENTS::Fluid3Surface* f3surf);
      
    };

    //-----------------------------------------------------------------
    ///     Internal Fluid3Surface weak Dirichlet implementation
    //-----------------------------------------------------------------
    template<DRT::Element::DiscretizationType distype ,
             DRT::Element::DiscretizationType pdistype>
    class Fluid3SurfaceWeakDBC: public Fluid3SurfaceWeakDBCInterface
    {
    public:
      /// Constructor with number of nodes
      Fluid3SurfaceWeakDBC();

      /// empty destructor
      ~Fluid3SurfaceWeakDBC(){return;};

      /// number of nodes
      static const int iel  = DRT::UTILS::DisTypeToNumNodePerEle<distype >::numNodePerElement;

      /// number of parentnodes
      static const int piel = DRT::UTILS::DisTypeToNumNodePerEle<pdistype>::numNodePerElement;

      /// Evaluate implementation 
      virtual int EvaluateWeakDBC(
        Fluid3Surface*             surfele       ,
        ParameterList&             params        ,
        DRT::Discretization&       discretization,
        vector<int>&               lm            ,
        Epetra_SerialDenseMatrix&  elemat_epetra ,
        Epetra_SerialDenseVector&  elevec_epetra );

    private:

      // nodal arrays
      // ------------
      //! node coordinates of parent element
      LINALG::Matrix<3,piel> pxyze_;
      //! array of nodal velocities, intermediate time level
      LINALG::Matrix<3,piel> pevelaf_;
      //! array of nodal velocities, new time level
      LINALG::Matrix<3,piel> pevelnp_;
      //! array of nodal pressure, new time level
      LINALG::Matrix<piel,1> peprenp_;
      //! node coordinates of boundary element
      LINALG::Matrix<3, iel> xyze_;



      // shape functions and derivatives, mapping from reference element to actual geometry
      // ----------------------------------------------------------------------------------
      //! transpose of the jacobian matrix of the mapping (r,s,t)->(x,y,z)
      LINALG::Matrix<3,3>    pxjm_;
      //! its inverse
      LINALG::Matrix<3,3>    pxji_;
      //! vector of shape functions, parent element
      LINALG::Matrix<piel,1> pfunct_;
      //! vector of shape function derivatives in reference coordinate system, parent element
      LINALG::Matrix<3,piel> pderiv_;
      //! vector of shape function derivatives in global coordinate system
      LINALG::Matrix<3,piel> pderxy_;

      //! vector of shape functions, boundary element
      LINALG::Matrix<iel ,1> funct_;
      //! vector of shape function derivatives in reference coordinate system, boundary element
      LINALG::Matrix<2, iel> deriv_;
      //! vector for outward normal, boundary element
      LINALG::Matrix<3,   1> n_;
      //! derivatives of surface in all reference directions
      LINALG::Matrix<2,3>    dxyzdrs_;
      //! the metric tensor
      LINALG::Matrix<2,2>    metrictensor_;
      //! the area of an infintesimal surface element
      double                 drs_;



      // values of non geometrical quantities in gausspoints
      // ---------------------------------------------------

      //! velocity in gausspoint, time n+af
      LINALG::Matrix<3,1>    velintaf_;
      //! velocity in gausspoint, time n+1
      LINALG::Matrix<3,1>    velintnp_;
      //! velocity derivatives in gausspoint, time n+af
      LINALG::Matrix<3,3>    vderxyaf_;
      //! pressure in gausspoint, time n+1
      double                 prenp_;

    };
  }
}

#endif

#endif
#endif
