/*----------------------------------------------------------------------*/
/*!
\file fluid2_genalpha_resVMM.H

\brief Internal implementation of Fluid2 element with a generalised alpha
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
#include <blitz/array.h>

namespace DRT
{
  namespace ELEMENTS
  {
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
      
    </pre>
        
    The stabilisation flags are read from the input file and passed
    through the evaluate interface.

    \auther gammi
    \date 02/08
  */
    class Fluid2GenalphaResVMM
    {
    public:
      /// Constructor with number of nodes
      explicit Fluid2GenalphaResVMM(int iel); 

      /// calculating element matrix and rhs
      /*
        Calculate matrix and rhs. Here the whole thing is hidden.

        \param ele             (i) the element those matrix is calculated
        \param myknots         (i) the elements knotvector segments for nurbs shape
                                   function evaluation
        \param estif           (o) element matrix to calculate
        \param eforce          (o) element rhs to calculate
        \param evelnp          (i) nodal velocities at n+1
        \param eprenp          (i) nodal pressure
        \param eaccam          (i) nodal accelerations at n+alphaM
        \param evelaf          (i) nodal velocities at n+alphaF
        \param material        (i) fluid material
        \param alphaM          (i) evaluate accelerations at this intermediate level
        \param alphaF          (i) evaluate velocities at this intermediate level
        \param gamma           (i) one step theta based parameter
        \param dt              (i) time step size
        \param time            (i) current simulation time
        \param newton          (i) flag for use of full linearisation of reactive terms
        \param tds             (i) flag for use of time dependent subscale stabilisation
        \param inertia         (i) flag for use of subscale accelerations for stabilisation (only for tds==true)
        \param pspg            (i) flag for use of nabla q testfunctions 
        \param supg            (i) flag for use of (u o nabla) v testfunctions 
        \param vstab           (i) flag for use of 2 * nu * grad o eps(v) testfunctions
        \param cross           (i) flag to add (u,(u_sub o nabla) v) to right hand side
        \param reynolds        (i) flag to add (u_sub,(u_sub o nabla) v) to right hand side
        \param compute_elemat  (i) if false, element matrices are not computed
        
      */
      void Sysmat(Fluid2*                                               ele,
                  std::vector<Epetra_SerialDenseVector>&                myknots,
                  Epetra_SerialDenseMatrix&                             elemat,
                  Epetra_SerialDenseVector&                             elevec,
                  const Epetra_SerialDenseMatrix&                       edispnp,
                  const Epetra_SerialDenseMatrix&                       egridvaf,
                  const Epetra_SerialDenseMatrix&                       evelnp,
                  const Epetra_SerialDenseVector&                       eprenp,
                  const Epetra_SerialDenseMatrix&                       eaccam,
                  const Epetra_SerialDenseMatrix&                       evelaf,
                  Teuchos::RCP<const MAT::Material>                     material,
                  const double                                          alphaM,
                  const double                                          alphaF,
                  const double                                          gamma,
                  const double                                          dt,
                  const double                                          time,
                  const bool                                            newton,
                  const enum Fluid2::StabilisationAction                tds,
                  const enum Fluid2::StabilisationAction                inertia,
                  const enum Fluid2::StabilisationAction                pspg,
                  const enum Fluid2::StabilisationAction                supg,
                  const enum Fluid2::StabilisationAction                vstab,
                  const enum Fluid2::StabilisationAction                cstab,
                  const enum Fluid2::StabilisationAction                cross,
                  const enum Fluid2::StabilisationAction                reynolds,
                  const bool                                            compute_elemat
        );



      void CalcRes(Fluid2*                                               ele,
                   const Epetra_SerialDenseMatrix&                       evelnp,
                   const Epetra_SerialDenseVector&                       eprenp,
                   const Epetra_SerialDenseMatrix&                       eaccam,
                   const Epetra_SerialDenseMatrix&                       evelaf,
                   Teuchos::RCP<const MAT::Material>                     material,
                   const double                                          alphaM,
                   const double                                          alphaF,
                   const double                                          gamma,
                   const double                                          dt,
                   const double                                          time,
                   const enum Fluid2::StabilisationAction                tds,
                   Epetra_SerialDenseVector&                             mean_res,
                   Epetra_SerialDenseVector&                             mean_sacc,
                   Epetra_SerialDenseVector&                             mean_res_sq,
                   Epetra_SerialDenseVector&                             mean_sacc_sq);

     
    private:
      ///  Get the nodal values of the bodyforce. The values are stored in edeadnf
      //
      //  param ele         (i) element pointer for access to the neumann condition
      //  param time        (i) time when the condition is evaluated
      void GetNodalBodyForce(Fluid2*      ele,
                             const double time);
      
      /*!
        \brief calculates viscosity at a given point (element center or Gaussian point)

        \param material                 (i) material
        \param visc                     (o) kinematic material viscosity
        */
      void CalVisc(     Teuchos::RCP<const MAT::Material>       material,
                        double&                                 visc);
      //! number of nodes --- determines array sizes
      int iel_;                             

      // nodal arrays
      // ------------
      //! node coordinates
      Epetra_SerialDenseMatrix xyze_;
      //! weights for nurbs elements
      Epetra_SerialDenseVector weights_;
      //! value of body force in the nodes of the element
      Epetra_SerialDenseMatrix edeadaf_;      

      // shape functions and derivatives, mapping from reference element to actual geometry
      // ----------------------------------------------------------------------------------
      //! vector of shape functions
      Epetra_SerialDenseVector funct_;
      //! vector of shape function derivatives in reference coordinate system
      Epetra_SerialDenseMatrix deriv_;
      //! vector of second shape function derivatives in reference coordinate system
      Epetra_SerialDenseMatrix deriv2_;
      //! vector of shape function derivatives in global coordinate system
      Epetra_SerialDenseMatrix derxy_;
      //! vector of second shape function derivatives in global coordinate system
      Epetra_SerialDenseMatrix derxy2_;
      //! derivatives needed for linearisation of viscous term including 2nd derivatives
      blitz::Array<double,3> viscs2_;
      /*  
                       corresponds to the operator
                     
                                 - grad o epsilon 
                       
                       on the discret level
      */
      //! transpose of the jacobian matrix of the mapping (r,s,t)->(x,y,z)
      Epetra_SerialDenseMatrix xjm_;
      //! its inverse
      Epetra_SerialDenseMatrix xji_;
      //! second derivatives of the mapping (r,s,t)->(x,y,z)
      Epetra_SerialDenseMatrix xder2_;       

      // values of non geometrical quantities in gausspoints
      // ---------------------------------------------------
      //! pressure in gausspoint, time n+1
      double                   prenp_;        
      //! acceleration in gausspoint, time n+am
      Epetra_SerialDenseVector accintam_;
      //! velocity in gausspoint, time n+1
      Epetra_SerialDenseVector velintnp_;
      //! velocity in gausspoint, time n+af
      Epetra_SerialDenseVector velintaf_;
      //! grid velocity at intermediate time level n+af
      Epetra_SerialDenseVector ugrid_af_;
      //! pressure gradient in gausspoint, time n+1
      Epetra_SerialDenseVector pderxynp_;
      //! velocity derivatives in gausspoint, time n+1
      Epetra_SerialDenseMatrix vderxynp_;
      //! velocity derivatives in gausspoint, time n+af
      Epetra_SerialDenseMatrix vderxyaf_;
      //! second derivatives of velocity in gausspoint, time n+af
      Epetra_SerialDenseMatrix vderxy2af_;
      //! value of body force interpolated to the gauss point, time n+af
      Epetra_SerialDenseVector bodyforceaf_;
      //! linearisation of convection, convective part, time n+af
      Epetra_SerialDenseVector conv_c_af_; 
      //! linearisation of convection, reactive part, time n+af
      blitz::Array<double,3> conv_r_af_;
      //! convective part introduced by grid velocity, time n+af
      Epetra_SerialDenseVector conv_g_af_; 

      // element data
      // ------------
      //! the stabilisation parameters tauM,tauMp and tauC
      Epetra_SerialDenseVector tau_;
      //! matrix of subscale velocities in gausspoints of this element, time n+alphaF
      Epetra_SerialDenseVector svelaf_;
      //! Convective term, most recent iteration value
      Epetra_SerialDenseVector convaf_old_;   
      /*            
                         /  n+af        \   n+af
                        |  u     o grad  | u
                         \  (i)         /   (i)
      */
      //! Subgrid velocity convective term, most recent iteration value 
      Epetra_SerialDenseVector convsubaf_old_;
      /* Most recent value for subgrid velocity convective term
        
                        /~n+af         \   n+af
                       | u      o nabla | u
                        \   (i)        /   (i)
      */
      //! Viscous term, most recent iteration value
      Epetra_SerialDenseVector viscaf_old_;
      /*                          
                                        /  n+af \
                       nabla o  epsilon | u     |
                                        \  (i)  /
      */
      //! Current residual of the momentum equation
      Epetra_SerialDenseVector resM_;
      //! This is the convection operator using resM_ as the convecting velocity
      Epetra_SerialDenseVector conv_resM_;
      /*  

                     /               \
                    | resM    o nabla |
                     \    (i)        /
      */
      //! This is the convection operator using the subscale velocity as the convecting velocity
      Epetra_SerialDenseVector  conv_subaf_; 
      /*  
                    /~n+af         \
                   | u      o nabla |
                    \   (i)        /
      */

      //! the element area, used to define element length
      double                    area_;
    };
  }
}

#endif

#endif
#endif
