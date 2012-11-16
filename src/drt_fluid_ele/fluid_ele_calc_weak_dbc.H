/*----------------------------------------------------------------------*/
/*!
\file fluid_ele_calc_weak_dbc.H

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
Maintainer: Volker Gravemeier
            gravemeier@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15245
</pre>
*/
/*----------------------------------------------------------------------*/

#ifndef FLUID_ELE_CALC_WEAK_DBC_H
#define FLUID_ELE_CALC_WEAK_DBC_H


#include "fluid_ele_parameter.H"
#include "../drt_fem_general/drt_utils_integration.H"
#include "../drt_fem_general/drt_utils_local_connectivity_matrices.H"

;
namespace DRT
{
  class Discretization;

  namespace ELEMENTS
  {
    class FluidBoundary;


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
    class FluidBoundaryWeakDBCInterface
    {
    public:
      //! Empty constructor
      FluidBoundaryWeakDBCInterface() {}
      //! Empty destructor
      virtual ~FluidBoundaryWeakDBCInterface() {}

      /*!
        \brief Evaluate the surface elements weak Dirichlet boundary
               conditions

        This class does not provide a definition for this function, it
        is defined in the implementation.

        \author gammi (gamnitzer@lnm.mw.tum.de), 11/09

        \param surfele             (in): surface element in geometry
                                         associated with weak dbc
        \param params              (in): parameters (time integration
                                         parent lms, condition)
        \param discretization      (in): the discretization
        \param lm                  (in): the location vector of surfele
        \param elemat_epetra   (in/out): element matrix (due to wdbc)
        \param elevec_epetra   (in/out): element force vector (due to wdbc)

      */
      virtual int EvaluateWeakDBC(
        FluidBoundary*             surfele       ,
        Teuchos::ParameterList&    params        ,
        DRT::Discretization&       discretization,
        std::vector<int>&          lm            ,
        Epetra_SerialDenseMatrix&  elemat_epetra ,
        Epetra_SerialDenseVector&  elevec_epetra ) = 0;

      /*!
        \brief Allocate one static instance of the internal
               implementation class for weak dirichlet condition and
               return pointer to it

        \author gammi (gamnitzer@lnm.mw.tum.de), 11/09

        \param f3surf (in): fluid surface element

      */
      static FluidBoundaryWeakDBCInterface* Impl(DRT::ELEMENTS::FluidBoundary* f3surf);

    };


    //-----------------------------------------------------------------
    //-----------------------------------------------------------------
    //
    //                        IMPLEMENTATION
    //
    //-----------------------------------------------------------------
    //-----------------------------------------------------------------

    /*!
    \brief Internal FluidSurface weak Dirichlet implementation
           (surface and parent element specific)
    */

    //-----------------------------------------------------------------
    ///
    //-----------------------------------------------------------------
    template<DRT::Element::DiscretizationType distype ,
             DRT::Element::DiscretizationType pdistype>
    class FluidSurfaceWeakDBC: public FluidBoundaryWeakDBCInterface
    {
    public:
      /// Constructor with number of nodes
      FluidSurfaceWeakDBC();

      /// empty destructor
      ~FluidSurfaceWeakDBC(){return;};

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
        FluidBoundary*             surfele       ,
        Teuchos::ParameterList&    params        ,
        DRT::Discretization&       discretization,
        std::vector<int>&          lm            ,
        Epetra_SerialDenseMatrix&  elemat_epetra ,
        Epetra_SerialDenseVector&  elevec_epetra );

    private:

      //! pointer to parameter list
      Teuchos::RCP<DRT::ELEMENTS::FluidEleParameter> fldpara_;

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
      //! array of nodal grid displacements, parent element, new time level
      LINALG::Matrix<3,piel> pedispnp_;
      //! array of nodal grid displacements, surface element, new time level
      LINALG::Matrix<3,piel> edispnp_;
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

    /*!
     \brief Internal Fluid2Line weak Dirichlet implementation
            (line and parent element specific)
     */

     //-----------------------------------------------------------------
     ///
     //-----------------------------------------------------------------
     template<DRT::Element::DiscretizationType distype ,
              DRT::Element::DiscretizationType pdistype>
     class FluidLineWeakDBC: public FluidBoundaryWeakDBCInterface
     {
     public:
       /// empty constructor
       FluidLineWeakDBC();

       /// empty destructor
       ~FluidLineWeakDBC(){return;};

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
         FluidBoundary*                lineele       ,
         Teuchos::ParameterList&    params        ,
         DRT::Discretization&       discretization,
         std::vector<int>&          lm            ,
         Epetra_SerialDenseMatrix&  elemat_epetra ,
         Epetra_SerialDenseVector&  elevec_epetra );

     private:

        //! pointer to parameter list
        Teuchos::RCP<DRT::ELEMENTS::FluidEleParameter> fldpara_;

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


    //-----------------------------------------------------------------
    //-----------------------------------------------------------------
    //
    //                    SPALDINGS LAW OF THE WALL
    //
    //-----------------------------------------------------------------
    //-----------------------------------------------------------------

    /*!
      \brief a class providing a method for the dynamic computation of
             the penalty parameter tauB for the weak Dirichlet boundary
             condition. It is computed such that the resulting traction
             corresponds to Spaldings formulation of the law of the
             wall

      Literature:
      Weak Dirichlet Boundary Conditions for Wall-Bounded Turbulent Flows
      Y. Bazilevs, C. Michler, V.M. Calo, and T.J.R. Hughes
      Computer Methods in Applied Mechanics and Engineering
      Volume 196, Issues 49-52, 1 November 2007, Pages 4853-4862

      and references therein

    */
    class FluidSurfaceWeakDBCSpaldingsLaw
    {
    public:
      /*! \brief constructor

      \param chi_in (in) empirical paramter of Spaldings law, usually 0.4
      \param B_in   (in) empirical paramter of Spaldings law, usually 5.5

      */
      FluidSurfaceWeakDBCSpaldingsLaw(const double chi_in,
                                       const double B_in  );

      /// Empty destructor
      virtual ~FluidSurfaceWeakDBCSpaldingsLaw() {}

      /*! \brief dynamic computation of the penalty paramter using Spaldings
                 law


      \param  tau_B (in/out) start/final value for penalty parameter tau_B
      \param  normu (in)     norm velocity in this point on surface
      \param  h     (in)     cell size in wall normal direction
      \param  Cb    (in)     constant originating from weak Dirichlet
                             framework
      \param  visc  (in)     kinematic viscosity

      */
      void ComputeTauBUsingSpaldingsLaw(
        double&       tau_B,
        const double& normu,
        const double& h    ,
        const double& Cb   ,
        const double& visc );

    private:
      //! first parameter of the Spalding model
      double chi_;

      //! second parameter of the Spalding model
      double B_;

      /*! \brief evaluate the residual of Spaldings law of the wall



      <pre>

                               +      / +\
                        res = y  - f | u  |
                                      \  /

                                  /                                                  \
                                 |                           /     +\ 2    /     +\ 3 |
                                 |       +                  | chi*u  |    | chi*u  |  |
    +     / +\     +    -chi*B   |  chi*u               +    \      /      \      /   |
   y = f | u  | = u  + e       * | e       - 1.0 - chi*u  - ----------- - ----------- |
          \  /                   |                              2.0           6.0     |
                                 |                                                    |
                                  \                                                  /


      </pre>

      \param  y     (in)     thickness of modeled boundary layer
      \param  nu    (in)     kinematic viscosity
      \param  tau_B (in)     value for penalty parameter tau_B
      \param  normu (in)     norm velocity in this point on surface

      \return residual of Spaldings equation
      */
      double SpaldingResidual(
        const double y     ,
        const double nu    ,
        const double tau_B ,
        const double normu
        );

      /*! \brief evaluate the residuals Jacobian of Spaldings law of the wall

      \param  y     (in)     thickness of modeled boundary layer
      \param  nu    (in)     kinematic viscosity
      \param  tau_B (in)     value for penalty parameter tau_B
      \param  normu (in)     norm velocity in this point on surface

      \return derivative of residual of Spaldings equation w.r.t. tau_B

      */
      double JacobianSpaldingResidual(
        const double y    ,
        const double nu   ,
        const double tau_B,
        const double normu);

      /*! \brief compute dimensionless velocity u+

      <pre>

                                               +-------------+
                     ||  n+af ||              /  ||  n+af||
           +         || u     ||             /   || u    ||
          u  = ---------------------- =     /   -------------
                 +------------------+    \ /        tau
                / tau  * ||  n+af ||      v            B
               v     B   || u     ||

      </pre>

      \param  tau_B (in)     value for penalty parameter tau_B
      \param  normu (in)     norm velocity in this point on surface

      \return u+

       */
      double Uplus(const double normu,
                   const double tau_B);

      /*! \brief compute dimensionless thickness of modeled layer y+

      <pre>
                   +-------------+              +-------------+
                  / tau  * ||u||          h *  / tau  * ||u||
           +     v     B                   b  v     B
          y  =  ------------------- * y = ---------------------- =
                        nu                      C  * nu
                                                 b
                  +-------------+
                 / tau  * ||u||
                v     B
             = -------------------
                     tau
                        B,0

      </pre>

      \param  normu (in)     norm velocity in this point on surface
      \param  tau_B (in)     value for penalty parameter tau_B
      \param  visc  (in)     kinematic viscosity
      \param  y     (in)     thickness of modeled boundary layer

      \return y+
       */
      double Yplus(const double normu,
                   const double tau_B,
                   const double visc ,
                   const double y    );

    };

  }
}

#endif
