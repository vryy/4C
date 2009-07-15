/*----------------------------------------------------------------------*/
/*!
\file fluid2_weak_dbc.H

\brief weak Dirichlet boundary conditions for fluid problems.

Literature:

    Weak imposition of Dirichlet boundary conditions in fluid mechanics
    Y. Bazilevs, T.J.R. Hughes
    Computers & Fluids 36 (2007) 12–26


    Weak Dirichlet Boundary Conditions for Wall-Bounded Turbulent Flows
    Y. Bazilevs, C. Michler, V.M. Calo, and T.J.R. Hughes
    Computer Methods in Applied Mechanics and Engineering
    Volume 196, Issues 49-52, 1 November 2007, Pages 4853-4862

<pre>
Maintainer: Peter Gamnitzer
            gamnitzer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>
*/
/*----------------------------------------------------------------------*/
#ifdef D_FLUID2
#ifdef CCADISCRET

#ifndef FLUID2_WEAKDBC_H
#define FLUID2_WEAKDBC_H

#include "fluid2.H"

#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_timecurve.H"
#include "../drt_lib/drt_function.H"

#include "../drt_nurbs_discret/drt_nurbs_discret.H"
#include "../drt_fem_general/drt_utils_nurbs_shapefunctions.H"

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
    //-----------------------------------------------------------------
    //
    //                        INTERFACE CLASS
    //
    //-----------------------------------------------------------------
    //-----------------------------------------------------------------

    /*-----------------------------------------------------------------

    \brief an interface class for the weak Dirichlet boundary condition
           for every surface element distype this class allocates one
           instance of the weak Dirichlet implementation

    \author gammi (gamnitzer@lnm.mw.tum.de), 11/09

    -----------------------------------------------------------------*/
    class Fluid2LineWeakDBCInterface
    {
    public:
      //! Empty constructor
      Fluid2LineWeakDBCInterface() {}
      //! Empty destructor
      virtual ~Fluid2LineWeakDBCInterface() {}

      /*!
        \brief Evaluate the line elements weak Dirichlet boundary
               conditions

        This class does not provide a definition for this function, it
        is defined in the implementation.

        \author gammi (gamnitzer@lnm.mw.tum.de), 11/09

        \param lineele             (in): line element in geometry
                                         associated with weak dbc
        \param params              (in): parameters (time integration
                                         parent lms, condition)
        \param discretization      (in): the discretization
        \param lm                  (in): the location vector of surfele
        \param elemat_epetra   (in/out): element matrix (due to wdbc)
        \param elevec_epetra   (in/out): element force vector (due to wdbc)

      */
      virtual int EvaluateWeakDBC(
        Fluid2Line*                lineele       ,
        ParameterList&             params        ,
        DRT::Discretization&       discretization,
        vector<int>&               lm            ,
        Epetra_SerialDenseMatrix&  elemat_epetra ,
        Epetra_SerialDenseVector&  elevec_epetra ) = 0;

      /*!
        \brief Allocate one static instance of the internal
               implementation class for weak dirichlet condition and
               return pointer to it

        \author gammi (gamnitzer@lnm.mw.tum.de), 11/09

        \param f2line (in): fluid line element

      */
      static Fluid2LineWeakDBCInterface* Impl(DRT::ELEMENTS::Fluid2Line* f2line);

    };


    //-----------------------------------------------------------------
    //-----------------------------------------------------------------
    //
    //                        IMPLEMENTATION
    //
    //-----------------------------------------------------------------
    //-----------------------------------------------------------------

    /*!
    \brief Internal Fluid2Line weak Dirichlet implementation
           (line and parent element specific)
    */

    //-----------------------------------------------------------------
    ///
    //-----------------------------------------------------------------
    template<DRT::Element::DiscretizationType distype ,
             DRT::Element::DiscretizationType pdistype>
    class Fluid2LineWeakDBC: public Fluid2LineWeakDBCInterface
    {
    public:
      /// Constructor with number of nodes
      Fluid2LineWeakDBC();

      /// empty destructor
      ~Fluid2LineWeakDBC(){return;};

      /// number of nodes
      static const int iel  = DRT::UTILS::DisTypeToNumNodePerEle<distype >::numNodePerElement;

      /// number of parentnodes
      static const int piel = DRT::UTILS::DisTypeToNumNodePerEle<pdistype>::numNodePerElement;


      /*!
        \brief Evaluate weakly imposed Dirichlet conditions

        This method calculates the contributions to rhs and matrix of
        weak Dirichlet boundary conditions for a generalized alpha
        system.

        Literature:
        Weak imposition of Dirichlet boundary conditions in fluid mechanics
        Y. Bazilevs, T.J.R. Hughes
        Computers & Fluids 36 (2007) 12–26

      */
      virtual int EvaluateWeakDBC(
        Fluid2Line*                lineele       ,
        ParameterList&             params        ,
        DRT::Discretization&       discretization,
        vector<int>&               lm            ,
        Epetra_SerialDenseMatrix&  elemat_epetra ,
        Epetra_SerialDenseVector&  elevec_epetra );

    private:

      // nodal arrays
      // ------------
      //! node coordinates of parent element
      LINALG::Matrix<2,piel> pxyze_;
      //! array of nodal velocities, intermediate time level
      LINALG::Matrix<2,piel> pevelaf_;
      //! array of nodal velocities, new time level
      LINALG::Matrix<2,piel> pevelnp_;
      //! array of nodal pressure, new time level
      LINALG::Matrix<piel,1> peprenp_;
      //! array of nodal grid displacements, parent element, new time level
      LINALG::Matrix<2,piel> pedispnp_;
      //! array of nodal grid displacements, surface element, new time level
      LINALG::Matrix<2,piel> edispnp_;
      //! node coordinates of boundary element
      LINALG::Matrix<2, iel> xyze_;



      // shape functions and derivatives, mapping from reference element to actual geometry
      // ----------------------------------------------------------------------------------
      //! transpose of the jacobian matrix of the mapping (r,s,t)->(x,y,z)
      LINALG::Matrix<2,2>    pxjm_;
      //! its inverse
      LINALG::Matrix<2,2>    pxji_;
      //! vector of shape functions, parent element
      LINALG::Matrix<piel,1> pfunct_;
      //! vector of shape function derivatives in reference coordinate system, parent element
      LINALG::Matrix<2,piel> pderiv_;
      //! vector of shape function derivatives in global coordinate system
      LINALG::Matrix<2,piel> pderxy_;


      //! vector of shape functions, boundary element
      LINALG::Matrix<iel ,1> funct_;
      //! vector of shape function derivatives in reference coordinate system, boundary element
      LINALG::Matrix<iel, 1> deriv_;
      //! vector for outward normal, boundary element
      LINALG::Matrix<2,   1> n_;
      //! derivatives of surface in reference direction
      LINALG::Matrix<2,   1> dxydr_;
      //! the length of an infintesimal line element
      double                 dr_;



      // values of non geometrical quantities in gausspoints
      // ---------------------------------------------------

      //! velocity in gausspoint, time n+af
      LINALG::Matrix<2,1>    velintaf_;
      //! velocity in gausspoint, time n+1
      LINALG::Matrix<2,1>    velintnp_;
      //! velocity derivatives in gausspoint, time n+af
      LINALG::Matrix<2,2>    vderxyaf_;
      //! pressure in gausspoint, time n+1
      double                 prenp_;

    };
  }
}
#endif

#endif
#endif
