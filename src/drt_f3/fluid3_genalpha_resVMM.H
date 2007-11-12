/*----------------------------------------------------------------------*/
/*!
\file fluid3_genalpha_resVMM.H

\brief Internal implementation of Fluid3 element with a generalised alpha
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

#ifndef FLUID3_GENALPHA_H
#define FLUID3_GENALPHA_H

#include "fluid3.H"
#include <blitz/array.h>

namespace DRT
{
  namespace Elements
  {
  /// Internal Fluid3 implementation (genalpha version)
  /*!
    This internal class keeps all the working arrays needed to
    calculate the Fluid3 element (genalpha implementation).
    Additionally the method Sysmat() provides this element
    implementation.

    <h3>Purpose</h3>

    The idea is to separate the element maintenance (class Fluid3)
    from the mathematical contents (this class). The differences
    to the default fluid implementation is the time integration
    scheme and the residual based stabilisation. This affects the
    factors of some of the contributions and some additional
    element contributions arise.

    The Fluid3 will allocate exactly one object of this class for all
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

    In addition to the stabilisation options, the method allows the use of
    classical turbulence models
    <pre>
    
    o Standard Smagorinsky
      --------------------
      With constant Smagorinsky coefficient

    o Smagorinsky with Van Driest damping
      -----------------------------------
      Exponential damping of mixing lenght close to the wall      
      
    </pre>
        
    The stabilisation flags as well as turbulence model options are have been read
    from the input file and passed through the evaluate interface.

    There is no ALE-ability of the element up to now.
    
    \auther gammi
    \date 11/07
  */
    class Fluid3GenalphaResVMM
    {
    public:
      /// Constructor with number of nodes
      explicit Fluid3GenalphaResVMM(int iel); 

      /// calculating element matrix and rhs
      /*
        Calculate matrix and rhs. Here the whole thing is hidden.

        \param ele             (i) the element those matrix is calculated
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
        \param linearisation   (i) flag for use of full linearisation of reactive terms
        \param tds             (i) flag for use of time dependent subscale stabilisation
        \param inertia         (i) flag for use of subscale accelerations for stabilisation (only for tds==true)
        \param pspg            (i) flag for use of nabla q testfunctions 
        \param supg            (i) flag for use of (u o nabla) v testfunctions 
        \param agls            (i) flag for use of 2 * nu * grad o eps(v) testfunctions
        \param cross           (i) flag to add (u,(u_sub o nabla) v) to right hand side
        \param reynolds        (i) flag to add (u_sub,(u_sub o nabla) v) to right hand side
        \param turb_mod_action (i) action flag for the application of a turbulence model
        \param Cs              (i) value for Smagorinsky constant (zero if unused)
        \param l_tau           (i) viscous length scale for van Driest damping (zero if unused)
        \param compute_elemat  (i) if false, element matrices are not computed
        
      */
      void Sysmat(Fluid3*                                               ele,
                  Epetra_SerialDenseMatrix&                             elemat,
                  Epetra_SerialDenseVector&                             elevec,
                  const blitz::Array<double,2>&                         evelnp,
                  const blitz::Array<double,1>&                         eprenp,
                  const blitz::Array<double,2>&                         eaccam,
                  const blitz::Array<double,2>&                         evelaf,
                  const struct _MATERIAL*                               material,
                  const double                                          alphaM,
                  const double                                          alphaF,
                  const double                                          gamma,
                  const double                                          dt,
                  const double                                          time,
                  const bool                                            newton,
                  const enum Fluid3::StabilisationAction                tds,
                  const enum Fluid3::StabilisationAction                inertia,
                  const enum Fluid3::StabilisationAction                pspg,
                  const enum Fluid3::StabilisationAction                supg,
                  const enum Fluid3::StabilisationAction                agls,
                  const enum Fluid3::StabilisationAction                cstab,
                  const enum Fluid3::StabilisationAction                cross,
                  const enum Fluid3::StabilisationAction                reynolds,
                  const enum Fluid3::TurbModelAction                    turb_mod_action,
                  const double                                          Cs,
                  const double                                          l_tau,
                  const bool                                            compute_elemat);
    private:
      ///  Get the nodal values of the bodyforce. The values are stored in edeadnf
      //
      //  param ele         (i) element pointer for access to the neumann condition
      //  param time        (i) time when the condition is evaluated
      void GetNodalBodyForce(Fluid3*      ele,
                             const double time);
      //! number of nodes --- determines array sizes
      int iel_;                             

      // nodal arrays
      // ------------
      //! node coordinates
      blitz::Array<double,2> xyze_;
      //! value of body force in the nodes of the element
      blitz::Array<double,2> edeadaf_;      

      // shape functions and derivatives, mapping from reference element to actual geometry
      // ----------------------------------------------------------------------------------
      //! vector of shape functions
      blitz::Array<double,1> funct_;
      //! vector of shape function derivatives in reference coordinate system
      blitz::Array<double,2> deriv_;
      //! vector of second shape function derivatives in reference coordinate system
      blitz::Array<double,2> deriv2_;
      //! vector of shape function derivatives in global coordinate system
      blitz::Array<double,2> derxy_;
      //! vector of second shape function derivatives in global coordinate system
      blitz::Array<double,2> derxy2_;
      //! derivatives needed for linearisation of viscous term including 2nd derivatives
      blitz::Array<double,3> viscs2_;
      /*  
                       corresponds to the operator
                     
                                 - grad o epsilon 
                       
                       on the discret level
      */
      //! transpose of the jacobian matrix of the mapping (r,s,t)->(x,y,z)
      blitz::Array<double,2> xjm_;
      //! its inverse
      blitz::Array<double,2> xji_;
      //! second derivatives of the mapping (r,s,t)->(x,y,z)
      blitz::Array<double,2> xder2_;       

      // values of non geometrical quantities in gausspoints
      // ---------------------------------------------------
      //! pressure in gausspoint, time n+1
      double                 prenp_;        
      //! acceleration in gausspoint, time n+am
      blitz::Array<double,1> accintam_;
      //! velocity in gausspoint, time n+1
      blitz::Array<double,1> velintnp_;
      //! velocity in gausspoint, time n+af
      blitz::Array<double,1> velintaf_;
      //! pressure gradient in gausspoint, time n+1
      blitz::Array<double,1> pderxynp_;
      //! velocity derivatives in gausspoint, time n+1
      blitz::Array<double,2> vderxynp_;
      //! velocity derivatives in gausspoint, time n+af
      blitz::Array<double,2> vderxyaf_;
      //! second derivatives of velocity in gausspoint, time n+af
      blitz::Array<double,2> vderxy2af_;
      //! value of body force interpolated to the gauss point, time n+af
      blitz::Array<double,1> bodyforceaf_;
      //! linearisation of convection, convective part, time n+af
      blitz::Array<double,1> conv_c_af_; 
      //! linearisation of convection, reactive part, time n+af
      blitz::Array<double,3> conv_r_af_;

      // element data
      // ------------
      //! the stabilisation parameters tauM,tauMp and tauC
      blitz::Array<double,1> tau_;
             //! matrix of subscale velocities in gausspoints of this element, time n+alphaF
      blitz::Array<double,1> svelaf_;
      //! Convective term, most recent iteration value
      blitz::Array<double,1> convaf_old_;   
      /*            
                         /  n+af        \   n+af
                        |  u     o grad  | u
                         \  (i)         /   (i)
      */
      //! Subgrid velocity convective term, most recent iteration value 
      blitz::Array<double,1> convsubaf_old_;
      /* Most recent value for subgrid velocity convective term
        
                        /~n+af         \   n+af
                       | u      o nabla | u
                        \   (i)        /   (i)
      */
      //! Viscous term, most recent iteration value
      blitz::Array<double,1> viscaf_old_;
      /*                          
                                        /  n+af \
                       nabla o  epsilon | u     |
                                        \  (i)  /
      */
      //! Current residual of the momentum equation
      blitz::Array<double,1> resM_;
      //! This is the convection operator using resM_ as the convecting velocity
      blitz::Array<double,1> conv_resM_;
      /*  

                     /               \
                    | resM    o nabla |
                     \    (i)        /
      */
      //! This is the convection operator using the subscale velocity as the convecting velocity
      blitz::Array<double,1>  conv_subaf_; 
      /*  
                    /~n+af         \
                   | u      o nabla |
                    \   (i)        /
      */
      //! the element volume, used to define element length
      double                  vol_;
      //! the effective viscosity (material+turbulent viscosity)
      double                  visceff_;
    };
  }
}

#endif

#endif
#endif
