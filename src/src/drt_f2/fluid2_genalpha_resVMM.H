/*----------------------------------------------------------------------*/
/*!
\file fluid2_genalpha_resVMM.H

\brief Internal implementation of Fluid2 element with a generalised alpha
       time integration.

       This element is designed for the solution of the Navier-Stokes
       equations using a residual based stabilised method. The
       stabilisation terms are derived in a variational multiscale sense.

       Subscales are either treated as quasi-static or time dependent.

       Both versions provide a conservative and an advective form. ALE
       is fully supported.
Internal implementation of Fluid2 element with a generalised alpha
       time integration.

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

#ifndef FLUID2_GENALPHA_H
#define FLUID2_GENALPHA_H

#include "fluid2.H"
#include "fluid2_nurbs.H"
#include "../drt_fem_general/drt_utils_nurbs_shapefunctions.H"
#include "../drt_lib/drt_utils.H"

namespace DRT
{
  namespace ELEMENTS
  {
    /// Interface base class for Fluid2GenalphaResVMM
    /*!
      This class exists to provide a common interface for all template
      versions of Fluid2GenalphaResVMM. The only function
      this class actually defines is Impl, which returns a pointer to
      the appropriate version of Fluid2GenalphaResVMM.
    */
    class Fluid2GenalphaResVMMInterface
    {
    public:
      /// Empty constructor
      Fluid2GenalphaResVMMInterface() {}
      /// Empty destructor
      virtual ~Fluid2GenalphaResVMMInterface() {}
      /// Evaluate the element
      /*!
        This class does not provide a definition for this function, it
        must be defined in Fluid3GenalphaResVMM.
      */
      virtual int Evaluate(Fluid2*                    ele,
                           ParameterList&             params,
                           DRT::Discretization&       discretization,
                           vector<int>&               lm,
                           Epetra_SerialDenseMatrix&  elemat1_epetra,
                           Epetra_SerialDenseMatrix&  elemat2_epetra,
                           Epetra_SerialDenseVector&  elevec1_epetra,
                           Epetra_SerialDenseVector&  elevec2_epetra,
                           Epetra_SerialDenseVector&  elevec3_epetra,
                           RefCountPtr<MAT::Material> mat) = 0;

      /// Internal implementation class for fluid element
      static Fluid2GenalphaResVMMInterface* Impl(DRT::ELEMENTS::Fluid2* f2);

    };


    /*----------------------------------------------------------------------*/
    /*----------------------------------------------------------------------*/
    /*----------------------------------------------------------------------*/
    /*----------------------------------------------------------------------*/


    /// Internal Fluid2 implementation (genalpha version)
    /*!

    This internal class keeps all the working arrays needed to
    calculate the Fluid2 element (genalpha implementation).
    Additionally the method Sysmat() provides this element
    implementation.

    <h3>Purpose</h3>

    The idea is to separate the element maintenance (class Fluid2)
    from the mathematical contents (this class). The differences
    to the default fluid implementation is the time integration
    scheme and the residual based stabilisation. This affects the
    factors of some of the contributions and some additional
    element contributions arise.

    The Fluid2 will allocate exactly one object of this class for all
    fluid elements with the same number of nodes in the mesh. This
    allows us to use exactly matching working arrays (and keep them
    around.)

    The number of working arrays has to be reduced to
    a minimum so that the element fits into the cache.
    (There might be room for improvements.)

    <h3>Description</h3>

    The implementation here is based on a generalised alpha time
    discretisation with an implicit treatment of the continuity
    equation and pressure, assuring that every solution is divergence
    free on output.

    There are multiple stabilisation options included in this
    implementation. For the standard residual based method these are:
    <pre>

    o Newton/Quasi-Newton (on/off)
      -------------------
      Reactive terms are/are not linearised

    o SUPG Stabilisation
      ------------------

    o PSPG Stabilisation
      ------------------

    o Viscous Stabilisation
      ---------------------
      Of the (A)GLS type

    o Cross stress (ResVMM)
      ---------------------
      Fix point treatment, no linearisation

    o Reynolds stress (ResVMM)
      ------------------------
      Fix point treatment, no linearisation

    </pre>
    For the choice of the stabilisation parameter see the inline
    documentation.

    For time dependent subscale stabilisation we have an additional
    option:
    <pre>

    o Stabilisation using subscale accelerations
      ------------------------------------------
      Equivalent to include a reaction like expression into the
      stablilisation operator --- should be included if any remarkable
      change on the matrix should occur in comparison to quasistatic
      subscales formulation.

      It means to stop to neglect the time derivative of the subgrid
      scales in the coarse scale equation. Obviously, this one is used
      only in the case of time-dependent subgrid scales. We allow to
      switch it off for debugging purposes.

    </pre>

    The stabilisation flags are read from the input file and passed
    through the evaluate interface.

    \auther gammi
    \date 05/09
  */

  template<DRT::Element::DiscretizationType distype>
  class Fluid2GenalphaResVMM: public Fluid2GenalphaResVMMInterface
  {
  public:
    /// Constructor with number of nodes
    Fluid2GenalphaResVMM();

    /// empty destructor
    ~Fluid2GenalphaResVMM(){return;};

    //! number of nodes
    static const int iel = DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;


    //! @name Control routine calling the specific stabilised (non)conservative form

    /*!
      Evaluate

      The evaluate function for the generalized-alpha fluid case.

      Get statevetors, extract values required by the element using the
      location vector.

      Set parameters for time integration and nonlinear iteration  (i.e.
      read control parameters from the parameterlist)

      Do setup for stabilisation and so on.

      Get knots for the nurbs implementation.

      In the end, call the specific Sysmat implementation, i.e. assemble
      element residual (and matrix).

      \param ele            (i)   the element to evaluate
      \param params         (i/o) list of all parameters --- Smagorinsky model
                                  stores some stuff in there to return it to
                                  the time integration level
      \param discretization (i)   the discretisation containing knot vectors etc
      \param lm             (i)   location vector to acces element dofs in
                                  global vectors
      \param elemat1_epetra (i/o) element matrix
      \param elemat2_epetra (-)   dummy
      \param elevec1_epetra (i/o) element residual (negative)
      \param elevec2_epetra (-)   dummy
      \param elevec3_epetra (-)   dummy
      \param mat            (i)   interface class for materials in elements

    */
      virtual int Evaluate(Fluid2*                    ele,
                           ParameterList&             params,
                           DRT::Discretization&       discretization,
                           vector<int>&               lm,
                           Epetra_SerialDenseMatrix&  elemat1_epetra,
                           Epetra_SerialDenseMatrix&  elemat2_epetra,
                           Epetra_SerialDenseVector&  elevec1_epetra,
                           Epetra_SerialDenseVector&  elevec2_epetra,
                           Epetra_SerialDenseVector&  elevec3_epetra,
                           RefCountPtr<MAT::Material> mat);

      //@}

      //! @name routines for calculation of element matrix and rhs

      /*!
        Calculate matrix and rhs in advective form. Here the whole thing is hidden
        for the quasistatic approach.

        \param ele                (i) the element those matrix is calculated
    	\param myknots            (i) the elements knot vector for nurbs elements
        \param elemat             (o) element matrix to calculate
        \param elevec             (o) element rhs to calculate
        \param edispnp            (i) mesh displacement (ALE formulation)
        \param egridvaf           (i) intermedoate grid velocity (ALE formulation)
        \param evelnp             (i) nodal velocities at n+1
        \param eprenp             (i) nodal pressure, time n+1
        \param eaccam             (i) nodal accelerations at n+alphaM
        \param evelaf             (i) nodal velocities at n+alphaF
        \param fsevelaf           (i) nodal fine-scale velocity, time n+af
        \param material           (i) fluid material
        \param alphaM             (i) evaluate accelerations at this intermediate level
        \param alphaF             (i) evaluate velocities at this intermediate level
        \param gamma              (i) one step theta based parameter
        \param dt                 (i) time step size
        \param time               (i) current simulation time
        \param newton             (i) flag for use of full linearisation of reactive terms
        \param higher_order_ele   (i) include second derivatives in residual for stabilization
        \param inertia            (i) flag for use of subscale accelerations for stabilisation (only for tds==true)
        \param pspg               (i) flag for use of nabla q testfunctions
        \param supg               (i) flag for use of (u o nabla) v testfunctions
        \param vstab              (i) flag for use of 2 * nu * grad o eps(v) testfunctions
        \param cross              (i) flag to add (u,(u_sub o nabla) v) to right hand side
        \param reynolds           (i) flag to add (u_sub,(u_sub o nabla) v) to right hand side
        \param whichtau           (i) flag choosing the tau version used for stabilization
        \param visceff            (o) effective viscosity computed for non Newtonian fluid etc.
        \param compute_elemat     (i) if false, element matrices are not computed

      */
      void Sysmat_adv_qs(
        Fluid2*                                ele,
        std::vector<Epetra_SerialDenseVector>& myknots,
        LINALG::Matrix<3*iel,3*iel>&           elemat,
        LINALG::Matrix<3*iel,    1>&           elevec,
        const LINALG::Matrix<2,iel>&           edispnp,
        const LINALG::Matrix<2,iel>&           egridvaf,
        const LINALG::Matrix<2,iel>&           evelnp,
        const LINALG::Matrix<iel,1>&           eprenp,
        const LINALG::Matrix<2,iel>&           eaccam,
        const LINALG::Matrix<2,iel>&           evelaf,
        Teuchos::RCP<const MAT::Material>      material,
        const double                           alphaM,
        const double                           alphaF,
        const double                           gamma,
        const double                           dt,
        const double                           time,
        const enum Fluid2::LinearisationAction newton,
        const bool                             higher_order_ele,
        const enum Fluid2::StabilisationAction inertia,
        const enum Fluid2::StabilisationAction pspg,
        const enum Fluid2::StabilisationAction supg,
        const enum Fluid2::StabilisationAction vstab,
        const enum Fluid2::StabilisationAction cstab,
        const enum Fluid2::StabilisationAction cross,
        const enum Fluid2::StabilisationAction reynolds,
        const enum Fluid2::TauType             whichtau,
        double&                                visceff,
        const bool                             compute_elemat
        );

    /*!
      Calculate matrix and rhs in advective form using time-dependent subgrid scales

        \param ele                (i) the element those matrix is calculated
    	\param myknots            (i) the elements knot vector for nurbs elements
        \param elemat             (o) element matrix to calculate
        \param elevec             (o) element rhs to calculate
        \param edispnp            (i) mesh displacement (ALE formulation)
        \param egridvaf           (i) intermedoate grid velocity (ALE formulation)
        \param evelnp             (i) nodal velocities at n+1
        \param eprenp             (i) nodal pressure, time n+1
        \param eaccam             (i) nodal accelerations at n+alphaM
        \param evelaf             (i) nodal velocities at n+alphaF
        \param material           (i) fluid material
        \param alphaM             (i) evaluate accelerations at this intermediate level
        \param alphaF             (i) evaluate velocities at this intermediate level
        \param gamma              (i) one step theta based parameter
        \param dt                 (i) time step size
        \param time               (i) current simulation time
        \param newton             (i) flag for use of full linearisation of reactive terms
        \param higher_order_ele   (i) include second derivatives in residual for stabilization
        \param inertia            (i) flag for use of subscale accelerations for stabilisation (only for tds==true)
        \param pspg               (i) flag for use of nabla q testfunctions
        \param supg               (i) flag for use of (u o nabla) v testfunctions
        \param vstab              (i) flag for use of 2 * nu * grad o eps(v) testfunctions
        \param cross              (i) flag to add (u,(u_sub o nabla) v) to right hand side
        \param reynolds           (i) flag to add (u_sub,(u_sub o nabla) v) to right hand side
        \param whichtau           (i) flag choosing the tau version used for stabilization
        \param visceff            (o) effective viscosity computed for non Newtonian fluid etc.
        \param compute_elemat     (i) if false, element matrices are not computed

      */
      void Sysmat_adv_td(
        Fluid2*                                ele,
        std::vector<Epetra_SerialDenseVector>& myknots,
        LINALG::Matrix<3*iel,3*iel>&           elemat,
        LINALG::Matrix<3*iel,    1>&           elevec,
        const LINALG::Matrix<2,iel>&           edispnp,
        const LINALG::Matrix<2,iel>&           egridvaf,
        const LINALG::Matrix<2,iel>&           evelnp,
        const LINALG::Matrix<iel,1>&           eprenp,
        const LINALG::Matrix<2,iel>&           eaccam,
        const LINALG::Matrix<2,iel>&           evelaf,
        Teuchos::RCP<const MAT::Material>      material,
        const double                           alphaM,
        const double                           alphaF,
        const double                           gamma,
        const double                           dt,
        const double                           time,
        const enum Fluid2::LinearisationAction newton,
        const bool                             higher_order_ele,
        const enum Fluid2::StabilisationAction inertia,
        const enum Fluid2::StabilisationAction pspg,
        const enum Fluid2::StabilisationAction supg,
        const enum Fluid2::StabilisationAction vstab,
        const enum Fluid2::StabilisationAction cstab,
        const enum Fluid2::StabilisationAction cross,
        const enum Fluid2::StabilisationAction reynolds,
        const enum Fluid2::TauType             whichtau,
        double&                                visceff,
        const bool                             compute_elemat
        );

      /*!
        Calculate matrix and rhs. Here the whole thing is hidden.
        Conservative, quasistatic version

        \param ele                (i) the element those matrix is calculated
    	\param myknots            (i) the elements knot vector for nurbs elements
        \param elemat             (o) element matrix to calculate
        \param elevec             (o) element rhs to calculate
        \param edispnp            (i) mesh displacement (ALE formulation)
        \param egridvaf           (i) intermedoate grid velocity (ALE formulation)
        \param evelnp             (i) nodal velocities at n+1
        \param eprenp             (i) nodal pressure, time n+1
        \param eaccam             (i) nodal accelerations at n+alphaM
        \param evelaf             (i) nodal velocities at n+alphaF
        \param material           (i) fluid material
        \param alphaM             (i) evaluate accelerations at this intermediate level
        \param alphaF             (i) evaluate velocities at this intermediate level
        \param gamma              (i) one step theta based parameter
        \param dt                 (i) time step size
        \param time               (i) current simulation time
        \param newton             (i) flag for use of full linearisation of reactive terms
        \param higher_order_ele   (i) include second derivatives in residual for stabilization
        \param pspg               (i) flag for use of nabla q testfunctions
        \param supg               (i) flag for use of (u o nabla) v testfunctions
        \param vstab              (i) flag for use of 2 * nu * grad o eps(v) testfunctions
        \param cross              (i) flag to add (u,(u_sub o nabla) v) to right hand side
        \param reynolds           (i) flag to add (u_sub,(u_sub o nabla) v) to right hand side
        \param whichtau           (i) flag choosing the tau version used for stabilization
        \param visceff            (o) effective viscosity computed for non Newtonian fluid etc.
        \param compute_elemat     (i) if false, element matrices are not computed

      */
      void Sysmat_cons_qs(
        Fluid2*                                ele,
        std::vector<Epetra_SerialDenseVector>& myknots,
        LINALG::Matrix<3*iel,3*iel>&           elemat,
        LINALG::Matrix<3*iel,    1>&           elevec,
        const LINALG::Matrix<2,iel>&           edispnp,
        const LINALG::Matrix<2,iel>&           egridvaf,
        const LINALG::Matrix<2,iel>&           evelnp,
        const LINALG::Matrix<iel,1>&           eprenp,
        const LINALG::Matrix<2,iel>&           eaccam,
        const LINALG::Matrix<2,iel>&           evelaf,
        Teuchos::RCP<const MAT::Material>      material,
        const double                           alphaM,
        const double                           alphaF,
        const double                           gamma,
        const double                           dt,
        const double                           time,
        const enum Fluid2::LinearisationAction newton,
        const bool                             higher_order_ele,
        const enum Fluid2::StabilisationAction pspg,
        const enum Fluid2::StabilisationAction supg,
        const enum Fluid2::StabilisationAction vstab,
        const enum Fluid2::StabilisationAction cstab,
        const enum Fluid2::StabilisationAction cross,
        const enum Fluid2::StabilisationAction reynolds,
        const enum Fluid2::TauType             whichtau,
        double&                                visceff,
        const bool                             compute_elemat
        );

      /*!
        Calculate matrix and rhs. Here the whole thing is hidden.
        Conservative, time-dependent version

        \param ele                (i) the element those matrix is calculated
    	\param myknots            (i) the elements knot vector for nurbs elements
        \param elemat             (o) element matrix to calculate
        \param elevec             (o) element rhs to calculate
        \param edispnp            (i) mesh displacement (ALE formulation)
        \param egridvaf           (i) intermedoate grid velocity (ALE formulation)
        \param evelnp             (i) nodal velocities at n+1
        \param eprenp             (i) nodal pressure, time n+1
        \param eaccam             (i) nodal accelerations at n+alphaM
        \param evelaf             (i) nodal velocities at n+alphaF
        \param material           (i) fluid material
        \param alphaM             (i) evaluate accelerations at this intermediate level
        \param alphaF             (i) evaluate velocities at this intermediate level
        \param gamma              (i) one step theta based parameter
        \param dt                 (i) time step size
        \param time               (i) current simulation time
        \param newton             (i) flag for use of full linearisation of reactive terms
        \param higher_order_ele   (i) include second derivatives in residual for stabilization
        \param inertia            (i) flag for use of subscale accelerations for stabilisation (only for tds==true)
        \param pspg               (i) flag for use of nabla q testfunctions
        \param supg               (i) flag for use of (u o nabla) v testfunctions
        \param vstab              (i) flag for use of 2 * nu * grad o eps(v) testfunctions
        \param cross              (i) flag to add (u,(u_sub o nabla) v) to right hand side
        \param reynolds           (i) flag to add (u_sub,(u_sub o nabla) v) to right hand side
        \param whichtau           (i) flag choosing the tau version used for stabilization
        \param visceff            (o) effective viscosity computed for non Newtonian fluid etc.
        \param compute_elemat     (i) if false, element matrices are not computed

      */
      void Sysmat_cons_td(
        Fluid2*                                ele,
        std::vector<Epetra_SerialDenseVector>& myknots,
        LINALG::Matrix<3*iel,3*iel>&           elemat,
        LINALG::Matrix<3*iel,    1>&           elevec,
        const LINALG::Matrix<2,iel>&           edispnp,
        const LINALG::Matrix<2,iel>&           egridvaf,
        const LINALG::Matrix<2,iel>&           evelnp,
        const LINALG::Matrix<iel,1>&           eprenp,
        const LINALG::Matrix<2,iel>&           eaccam,
        const LINALG::Matrix<2,iel>&           evelaf,
        Teuchos::RCP<const MAT::Material>      material,
        const double                           alphaM,
        const double                           alphaF,
        const double                           gamma,
        const double                           dt,
        const double                           time,
        const enum Fluid2::LinearisationAction newton,
        const bool                             higher_order_ele,
        const enum Fluid2::StabilisationAction inertia,
        const enum Fluid2::StabilisationAction pspg,
        const enum Fluid2::StabilisationAction supg,
        const enum Fluid2::StabilisationAction vstab,
        const enum Fluid2::StabilisationAction cstab,
        const enum Fluid2::StabilisationAction cross,
        const enum Fluid2::StabilisationAction reynolds,
        const enum Fluid2::TauType             whichtau,
        double&                                visceff,
        const bool                             compute_elemat
        );

      //@}

    private:

      //! @name service methods

      ///  extract velocities, pressure and accelerations etc from the global distributed vectors

      //
      //  \param is_ale         (i) Control parameter: get mesh/displacement
      //  \param discretization (i) discretization containing the global vectors
      //  \param lm             (i) location vector for element dofs in global vector
      //  \param eprenp         (o) nodal pressure values, t_{n+1}
      //  \param evelnp         (o) nodal velocities, t_{n+1}
      //  \param evelaf         (o) nodal velocities, t_{n+alpha_F}
      //  \param eaccam         (o) nodal accelerations, t_{n+alpha_M}
      //  \param edispnp        (o) nodal mesh displacements, t_{n+1} (for ALE)
      //  \param egridvelaf     (o) nodal mesh velocities, t_{n+alpha_F} (for ALE)

      void ExtractValuesFromGlobalVectors(
        const bool                        is_ale        ,
        const DRT::Discretization&        discretization,
        const vector<int>&                lm            ,
        LINALG::Matrix<iel,1>&            eprenp        ,
        LINALG::Matrix<2,iel>&            evelnp        ,
        LINALG::Matrix<2,iel>&            evelaf        ,
        LINALG::Matrix<2,iel>&            eaccam        ,
        LINALG::Matrix<2,iel>&            edispnp       ,
        LINALG::Matrix<2,iel>&            egridvelaf
        );


      ///  Get the nodal values of the bodyforce. The values are stored in edeadnf
      //
      //  \param ele         (i) element pointer for access to the neumann condition
      //  \param time        (i) time when the condition is evaluated
      void GetNodalBodyForce(const Fluid2* ele ,
                             const double  time);

      /*!  Compute the rate of strain of a (fine-scale) velocity

      //  \param evel       (i) nodal velocity values
      //  \param derxy      (i) shape function derivatives
      //  \param velderxy   (o) velocity derivatives
      //
      //  \return computed rate of strain

      */
      double GetStrainRate(
        const LINALG::Matrix<2,iel>& evel     ,
        const LINALG::Matrix<2,iel>& derxy    ,
        LINALG::Matrix<2,2>        & velderxy
        )
        {
          double rateofstrain=0;

          // get velocity derivatives at integration point
          //
          //              +-----  dN (x)
          //   dvel (x)    \        k
          //   -------- =   +     ------ * vel
          //      dx       /        dx        k
          //        j     +-----      j
          //              node k
          //
          // j : direction of derivative x/y
          //
          for(int rr=0;rr<2;++rr)
          {
            for(int mm=0;mm<2;++mm)
            {
              velderxy(rr,mm)=derxy(mm,0)*evel(rr,0);
              for(int nn=1;nn<iel;++nn)
              {
                velderxy(rr,mm)+=derxy(mm,nn)*evel(rr,nn);
              }
            }
          }

          // compute (resolved) rate of strain
          //
          //          +-                                 -+ 1
          //          |          /   \           /   \    | -
          //          | 2 * eps | vel |   * eps | vel |   | 2
          //          |          \   / ij        \   / ij |
          //          +-                                 -+
          //
          LINALG::Matrix<2,2> two_epsilon;

          for(int rr=0;rr<2;++rr)
          {
            for(int mm=0;mm<2;++mm)
            {
              two_epsilon(rr,mm) = velderxy(rr,mm) + velderxy(mm,rr);
            }
          }

          for(int rr=0;rr<2;rr++)
          {
            for(int mm=0;mm<2;mm++)
            {
              rateofstrain += two_epsilon(rr,mm)*two_epsilon(rr,mm);
            }
          }
          // sqrt(two_epsilon(rr,mm)*two_epsilon(rr,mm)/4.0*2.0)

          return(sqrt(rateofstrain/2.0));
        }

      /*!
        \brief Get all global shape functions, first and eventually second
               derivatives in a gausspoint

    	\param ele       (i) pointer to element for debug io
    	\param iquad     (i) number of gausspoint
    	\param intpoints (i) array with gausspoints and weights
    	\param myknots   (i) knotvector for NURBS elements
    	\param hoel      (i) flag for the evaluation of second derivatives

        \return integration factor (weight)

       */
      double ShapeFunctionsFirstAndSecondDerivatives(
        const Fluid2*                                 ele      ,
        const int                                  &  iquad    ,
        const DRT::UTILS::IntegrationPoints2D      &  intpoints,
        const std::vector<Epetra_SerialDenseVector>&  myknots  ,
        const bool                                    hoel
        );


      /*!
        \brief calculates viscosity at a given point (element center or Gaussian point)

    	\param material		 	(i) material
        \param visc        		(o) kinematic material viscosity
    	\param rateofshear		(i) shear rate
       */
      void CalVisc(Teuchos::RCP<const MAT::Material> material,
                   double&                           visc,
                   const double &                    rateofshear);

      /*!
        \brief calculation of stabilisation parameter (element center or Gaussian point)

    	\param whichtau   (i) select tau definition
    	\param tds        (i) tau definitions differ for time dependent subscales
        \param gamma      (i) time integration dependent parameter
        \param dt         (i) time step size
        \param hk         (i) an element size required for all stabpars except the Bazilevs parameter
    	\param mk         (i) element dependent constant
        \param visceff    (i) effective viscosity
       */
      void CalcTau(const enum Fluid2::TauType             whichtau ,
                   const enum Fluid2::StabilisationAction tds      ,
                   const double &                         gamma    ,
                   const double &                         dt       ,
                   const double &                         hk       ,
                   const double &                         mk       ,
                   const double &                         visceff  );


      /*!
        \brief calculates all quantities which are defined at the element center
               or for the whole element
                        o element geometry (xyze_ etc)
                        o element volume vol_
                        o element size hk, constant mk from inverse estimate
                        o dead load
                        o viscosity, effective viscosity


     	\param ele                      (i) pointer to element
        \param edispnp                  (i) nodal mesh displacements, t_{n+1} (for ALE)
        \param evelaf                   (i) nodal velocities, t_{n+alpha_F}
    	\param myknots                  (i) knotvector for NURBS elements
        \param timealphaF               (i) alphaF*t_{n+1}+(1-alphaF)*t_{n}
        \param hk                       (o) element length for stabilisation parameter
        \param mk                       (o) element type constant mk for tau and the fssgv_artificial approach
    	\param material		 	(i) material
        \param visc                     (o) kinematicmaterial viscosity
        \param visceff        		(o) effective material viscosity

       */
      void SetElementData
      (
        Fluid2*                                      ele            ,
        const LINALG::Matrix<2,iel>                & edispnp        ,
        const LINALG::Matrix<2,iel>                & evelaf         ,
        const std::vector<Epetra_SerialDenseVector>& myknots        ,
        const double                               & timealphaF     ,
        double                                     & hk             ,
        double                                     & mk             ,
        Teuchos::RCP<const MAT::Material>            material       ,
        double                                     & visc           ,
        double                                     & visceff        );


      /*!
        \brief Interpolates standard quantities to gausspoint, including

        o accintam_
        o velintaf_
        o velintnp_
        o bodyforce_
        o prenp_
        o pderxynp_
        o vderxyaf_
        o vderxynp_
        o divunp_
        o u_G_af
        o aleconvintaf_
        o conv_af_old_
        o resM_
        o conv_c_af_
        o viscs2_ (if hoel)
        o viscaf_old_ (if hoel)


     	\param ele               (i) pointer to element
        \param egridvelaf        (i) nodal mesh velocities, t_{n+alpha_F} (for ALE)
        \param evelnp            (i) nodal velocities, t_{n+1}
        \param eprenp            (i) nodal pressure values, t_{n+1}
        \param eaccam            (i) nodal accelerations, t_{n+alpha_M}
        \param evelaf            (i) nodal velocities, t_{n+alpha_F}
        \param visceff        	 (i) effective material viscosity
        \param higher_order_ele  (i) include second derivatives in residual for stabilization

       */
      void InterpolateToGausspoint(
        Fluid2*                                     ele             ,
        const LINALG::Matrix<2,iel>               & egridvaf        ,
        const LINALG::Matrix<2,iel>               & evelnp          ,
        const LINALG::Matrix<iel,1>               & eprenp          ,
        const LINALG::Matrix<2,iel>               & eaccam          ,
        const LINALG::Matrix<2,iel>               & evelaf          ,
        const double                              & visceff         ,
        const bool                                  higher_order_ele);

      //@}

      //! constant body force
      bool   constant_bodyforce_;

      // nodal arrays
      // ------------
      //! node coordinates
      LINALG::Matrix<2,iel> xyze_;
      //! weights for nurbs elements
      LINALG::Matrix<iel,1> weights_;
      //! value of body force in the nodes of the element
      LINALG::Matrix<2,iel> edeadaf_;

      // shape functions and derivatives, mapping from reference element to actual geometry
      // ----------------------------------------------------------------------------------
      //! vector of shape functions
      LINALG::Matrix<iel,1> funct_;
      //! vector of shape function derivatives in reference coordinate system
      LINALG::Matrix<2,iel> deriv_;
      //! vector of second shape function derivatives in reference coordinate system
      LINALG::Matrix<3,iel> deriv2_;
      //! vector of shape function derivatives in global coordinate system
      LINALG::Matrix<2,iel> derxy_;
      //! vector of second shape function derivatives in global coordinate system
      LINALG::Matrix<3,iel> derxy2_;
      //! derivatives needed for linearisation of viscous term including 2nd derivatives
      /*
                       corresponds to a part of the operator

                                  grad o epsilon

                       on the discret level
      */
      LINALG::Matrix<2,iel> viscs2_;
      //! transpose of the jacobian matrix of the mapping (r,s)->(x,y)
      LINALG::Matrix<2,2>   xjm_;
      //! its inverse
      LINALG::Matrix<2,2>   xji_;
      //! second derivatives of the mapping (r,s)->(x,y)
      LINALG::Matrix<3,2>   xder2_;

      // values of non geometrical quantities in gausspoints
      // ---------------------------------------------------
      //! pressure in gausspoint, time n+1
      double                prenp_;
      //! acceleration in gausspoint, time n+am
      LINALG::Matrix<2,1>   accintam_;
      //! velocity in gausspoint, time n+1
      LINALG::Matrix<2,1>   velintnp_;
      //! velocity in gausspoint, time n+af
      LINALG::Matrix<2,1>   velintaf_;
      // this copy will be used to store the normed velocity at time n+af
      LINALG::Matrix<2,1>   normed_velintaf_;
      //! fine-scale velocity in gausspoint, time n+af
      LINALG::Matrix<2,1>   fsvelintaf_;
      //! grid velocity at time level n+af
      LINALG::Matrix<2,1>   u_G_af_;
      //! ale convective velocity at intermediate time level n+af
      //! aleconvint_af_ = velintaf_ - gridvelintaf_
      LINALG::Matrix<2,1>   aleconvintaf_;
      //! pressure gradient in gausspoint, time n+1
      LINALG::Matrix<2,1>   pderxynp_;
      //! velocity derivatives in gausspoint, time n+1
      LINALG::Matrix<2,2>   vderxynp_;
      //! velocity derivatives in gausspoint, time n+af
      LINALG::Matrix<2,2>   vderxyaf_;
      //! value of body force interpolated to the gauss point, time n+af
      LINALG::Matrix<2,1>   bodyforceaf_;
      //! linearisation of convection, convective part, time n+af
      LINALG::Matrix<iel,1> conv_c_af_;
      //! linearisation of convection test function, convective part, time n+af
      //! may differ from conv_c_af_ if reynolds stress is used
      LINALG::Matrix<iel,1> conv_c_plus_svel_af_;
      //! convective part introduced by grid velocity, time n+af
      LINALG::Matrix<iel,1> conv_u_G_af_;

      // element data
      // ------------
      //! the stabilisation parameters tauM,tauMp and tauC
      LINALG::Matrix<3,1>   tau_;
      //! matrix of subscale velocities in gausspoints of this element, time n+alphaF
      LINALG::Matrix<2,1>   svelaf_;
      //! Convective term, most recent iteration value
      LINALG::Matrix<2,1>   convaf_old_;
      //! Convectivion by mesh motion, most recent iteration value
      LINALG::Matrix<2,1>   convu_G_af_old_;
      /*
                         /  n+af        \   n+af
                        |  u     o grad  | u
                         \  (i)         /   (i)
      */
      //! Subgrid velocity convective term, most recent iteration value
      LINALG::Matrix<2,1>   convsubaf_old_;
      /* Most recent value for subgrid velocity convective term

                        /~n+af         \   n+af
                       | u      o nabla | u
                        \   (i)        /   (i)
      */
      //! Viscous term, most recent iteration value
      LINALG::Matrix<2,1>   viscaf_old_;
      /*
                                        /  n+af \
                       nabla o  epsilon | u     |
                                        \  (i)  /
      */
      //! Current residual of the momentum equation
      LINALG::Matrix<2,1>   resM_;

      //! This is the convection operator using resM_ as the convecting velocity
      LINALG::Matrix<iel,1> conv_resM_;
      /*

                     /               \
                    | resM    o nabla |
                     \    (i)        /
      */
      //! This is the convection operator using the subscale velocity as the convecting velocity
      LINALG::Matrix<iel,1> conv_subaf_;
      /*
                    /~n+af         \
                   | u      o nabla |
                    \   (i)        /
      */

      //! Current residual of the continuity equation, time t_{n+1}
      double                divunp_;


      //! the element area, used to define element length
      double                area_;


      // temporary working arrays
      // ------------------------
      //! array used for the computation of jacobian
      //! determinant and second derivatives
      LINALG::Matrix<3,3>   bm_;

      //! An epetra wrapper for LAPACK library functions
      Epetra_LAPACK         solver_;

  };

  }
}

#endif

#endif
#endif
