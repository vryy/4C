
/*----------------------------------------------------------------------*/
/*!
\file fluid_genalpha_resVMM.H
FineSubgridVisc
\brief Internal implementation of Fluid element with a generalised alpha
       time integration. The major difference to the standard
       one-step-theta/BDF2 implementation is that the pressure and the
       continuity equation are treated fully implicitly.

       This element is designed for the solution of the Navier-Stokes
       equations using a residual based stabilised method. The
       stabilisation terms are derived in a variational multiscale sense.

       Subscales are either treated as quasi-static or time dependent.

       Both versions provide a conservative and an advective form. ALE
       is fully supported.

       There are options for Smagorinsky-type turbulence models.

<pre>
Maintainer: Peter Gamnitzer
            gamnitzer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>
*/
/*----------------------------------------------------------------------*/

#ifndef FLUID_GENALPHA_H
#define FLUID_GENALPHA_H


#include "../drt_inpar/inpar_fluid.H"
#include "../drt_lib/drt_element.H"
#include "../drt_fem_general/drt_utils_integration.H"
#include "../drt_fem_general/drt_utils_local_connectivity_matrices.H"


class ParameterList;
namespace DRT
{
  class Discretization;

  namespace ELEMENTS
  {
    class Fluid;

    /// Interface base class for FluidGenalphaResVMM

    /*!
      This class exists to provide a common interface for all template
      versions of FluidGenalphaResVMM. The only function
      this class actually defines is Impl, which returns a pointer to
      the appropriate version of FluidGenalphaResVMM.
    */
    class FluidGenalphaResVMMInterface
    {
    public:
      /// Empty constructor
      FluidGenalphaResVMMInterface() {}
      /// Empty destructor
      virtual ~FluidGenalphaResVMMInterface() {}
      /// Evaluate the element
      /*!
        This class does not provide a definition for this function, it
        must be defined in FluidGenalphaResVMM.
      */
      virtual int Evaluate(Fluid*                    ele,
                           ParameterList&             params,
                           DRT::Discretization&       discretization,
                           vector<int>&               lm,
                           Epetra_SerialDenseMatrix&  elemat1_epetra,
                           Epetra_SerialDenseMatrix&  elemat2_epetra,
                           Epetra_SerialDenseVector&  elevec1_epetra,
                           Epetra_SerialDenseVector&  elevec2_epetra,
                           Epetra_SerialDenseVector&  elevec3_epetra,
                           Teuchos::RCP<MAT::Material> mat) = 0;

      /// Calculate averages for residual, dissipation rates etc for
      /// elements in a turbulent channel flow. The computed values are
      /// are added to a vector associated with an element layer
      /*!
        This class does not provide a definition for this function, it
        must be defined in FluidGenalphaResVMM.
      */
      virtual int CalcResAvgs(Fluid*                    ele,
                              ParameterList&             params,
                              DRT::Discretization&       discretization,
                              vector<int>&               lm,
                              Teuchos::RCP<MAT::Material> mat) = 0;

      /// Internal implementation class for fluid element
      static FluidGenalphaResVMMInterface* Impl(DRT::ELEMENTS::Fluid* f3);

    };


    /*----------------------------------------------------------------------*/
    /*----------------------------------------------------------------------*/
    /*----------------------------------------------------------------------*/
    /*----------------------------------------------------------------------*/


    /// Internal Fluid implementation (genalpha version)
    /*!
      This internal class keeps all the working arrays needed to
      calculate the Fluid element (genalpha implementation).
      Additionally the methods Sysmat_...() provide this elements
      implementations.

      <h3>Purpose</h3>

      The idea is to separate the element maintenance (class Fluid)
      from the mathematical contents (this class). The differences
      to the default fluid implementation is the time integration
      scheme and the residual based stabilisation. This affects the
      factors of some of the contributions and some additional
      element contributions arise.

      The Fluid will allocate exactly one object of this class for all
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

      o Newton/Quasi-Newton/minimal (on/off)
        ------------------------------------
        Reactive terms are/are not linearised, minimal does not
        include second derivatives in tangent matrix

      o SUPG Stabilisation
        ------------------

      o PSPG Stabilisation
        ------------------

      o Viscous Stabilisation
        ---------------------
        Of the (A)GLS type

      o Cross stress (ResVMM)
        ---------------------
        Optional: Fix point treatment, no linearisation

      o Reynolds stress (ResVMM)
        ------------------------
        Optional: Fix point treatment, no linearisation

      </pre>

      Available versions of the stabilisation parameters are
      <pre>

      o tauM: Codina
        tauC: Codina

        (only for time-dependent subgrid scales)

      o tauM: Barrenechea, Franca , Valentin
        tauC: Wall

      o tauM: Barrenechea, Franca , Valentin
        tauC: Codina

        (only for quasistatic subgrid scales)

      o tauM: Bazilevs
        tauC: Bazilevs

      </pre>
      For further details on the choice of the stabilisation parameter
      see the inline documentation. Note that the formulations differ
      for time-dependent and quasistatic subgrid scales.

      For time dependent subscale stabilisation we have an additional
      option:
      <pre>

      o Stabilisation using subscale accelerations
        ------------------------------------------
        Stop to neglect the time derivative of the subgrid scales in the
        coarse scale equation. Obviously, this one is used only in the
        case of time-dependent subgrid scales. We allow to switch it off
        for debugging purposes.

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

      o Dynamic Smagorinsky
        -------------------
        Based on Germano's model, see for example. Note that Cs is averaged over
        planes in the channel flow and clipping of negative values of Cs is
        applied (and necessary to keep the algorithm stable)

        M. Germano, U. Piomelli, P. Moin, W.H. Cabot:
        A dynamic subgrid-scale eddy viscosity model
        (Phys. Fluids 1991)

        or

        D.K. Lilly:
        A proposed modification of the Germano subgrid-scale closure method
        (Phys. Fluids 1992)

      </pre>

      The stabilisation flags as well as turbulence model options have been read
      from the input file and are passed through the evaluate interface.

      \author gammi
      \date 04/08
    */

    template<DRT::Element::DiscretizationType distype>
    class FluidGenalphaResVMM: public FluidGenalphaResVMMInterface
    {
    public:
      /// Constructor with number of nodes
      FluidGenalphaResVMM(int numdofpernode);

      /// empty destructor
      ~FluidGenalphaResVMM(){return;};

      //! number of nodes
      static const int iel = DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;

      //! number of space dimensions
      static const int nsd_ = DRT::UTILS::DisTypeToDim<distype>::dim;

      //! @name Control routine calling the specific stabilised (non)conservative form

      /*!
        Evaluate

        The evaluate function for the generalized-alpha fluid case.

        Get statevetors, extract values required by the element using the
        location vector.

        Set parameters for time integration and nonlinear iteration  (i.e.
        read control parameters from the parameterlist)

        Do setup for stabilisation, turbulence models and so on.

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
        \param elevec1_epetra (i/o) (negative) element residual
        \param elevec2_epetra (-)   dummy
        \param elevec3_epetra (-)   dummy
        \param mat            (i)   interface class for materials in elements

      */
      virtual int Evaluate(Fluid*                    ele,
                           ParameterList&             params,
                           DRT::Discretization&       discretization,
                           vector<int>&               lm,
                           Epetra_SerialDenseMatrix&  elemat1_epetra,
                           Epetra_SerialDenseMatrix&  elemat2_epetra,
                           Epetra_SerialDenseVector&  elevec1_epetra,
                           Epetra_SerialDenseVector&  elevec2_epetra,
                           Epetra_SerialDenseVector&  elevec3_epetra,
                           Teuchos::RCP<MAT::Material> mat);

      //@}

      //! @name Data collection for turbulent channel flow computations

      /*!
        \brief evaluate residuals resM and resC on the element and return
               averaged values together with averaged subgrid-scale quantities.

        The calculation of averages for residual, dissipation rates etc for
        elements in a turbulent channel flow is based on a generalized-alpha
        time integration.

        The computed values are added to a vector associated with an
        element layer and communicated with the time integration via the
        parameterlist.

        The following things are done:

        o Get statevetors, extract values required by the element using the
          location vector.

        o Set parameters for time integration (i.e. read control parameters
          from the parameterlist)

        o Do setup for stabilization, turbulence models and so on.

        o Get knots for the nurbs implementation.

        o evaluate the gausspoint quantities, average over the element

        o add values to layer vectors

        Note:

        In large parts, this method is a copy of Sysmat... (concerning
        the computations of the residuals). Its only purpose is to
        allow an output of time averaged residuals and subscale
        accelerations to study their influence


        \param ele            (i)   the element to evaluate
        \param params         (i/o) list of all parameters; returns averages to
                                    time integration
        \param discretization (i)   the discretisation containing knot vectors etc
        \param lm             (i)   location vector to acces element dofs in
                                    global vectors
        \param mat            (i)   interface class for materials in elements
        \param actmat         (i)   actual material storing density

      */

      virtual int CalcResAvgs(Fluid*                    ele,
                              ParameterList&             params,
                              DRT::Discretization&       discretization,
                              vector<int>&               lm,
                              Teuchos::RCP<MAT::Material> mat);
      //@}

      //! @name routines for calculation of element matrix and rhs

      /*!
        Calculate matrix and rhs in advective form. Here the whole thing is hidden
        for the quasistatic approach.

        \param ele                (i) the element those matrix is calculated
    	\param myknots            (i) the elements knot vector for nurbs elements
        \param elemat             (o) element matrix to calculate
        \param elemat             (o) linearisations wrt ale motion (for monolithic fsi)
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
        \param fssgv              (i) flag for type of fine-scale subgrid viscosity
        \param tds                (i) flag for use of time dependent subscale stabilisation
        \param inertia            (i) flag for use of subscale accelerations for stabilisation (only for tds==true)
        \param pspg               (i) flag for use of nabla q testfunctions
        \param supg               (i) flag for use of (u o nabla) v testfunctions
        \param vstab              (i) flag for use of 2 * nu * grad o eps(v) testfunctions
        \param cross              (i) flag to add (u,(u_sub o nabla) v) to right hand side
        \param reynolds           (i) flag to add (u_sub,(u_sub o nabla) v) to right hand side
        \param whichtau           (i) flag choosing the tau version used for stabilization
        \param turb_mod_action    (i) action flag for the application of a turbulence model
        \param Cs               (i/o) value for Smagorinsky constant (zero if unused)
                                      in case of the dynamic model, it is only equal to the
                                      Smagorinsky constant on output.
        \param Cs_delta_sq        (i) Used for the dynamic Smagorinsky model. Corresponds to (Cs*hk)^2
        \param visceff            (o) effective viscosity computed by dynamic Smagorinsky model,
                                      non Newtonian fluid etc.
        \param l_tau              (i) viscous length scale for van Driest damping (zero if unused)
        \param compute_elemat     (i) if false, element matrices are not computed

      */
      void Sysmat_adv_qs(Fluid*                                    ele,
                         std::vector<Epetra_SerialDenseVector>&     myknots,
                         LINALG::Matrix<(nsd_+1)*iel,(nsd_+1)*iel>& elemat,
                         LINALG::Matrix<(nsd_+1)*iel,(nsd_+1)*iel>& meshmat,
                         LINALG::Matrix<(nsd_+1)*iel,    1>&        elevec,
                         const LINALG::Matrix<nsd_,iel>&            edispnp,
                         const LINALG::Matrix<nsd_,iel>&            egridvaf,
                         const LINALG::Matrix<nsd_,iel>&            evelnp,
                         const LINALG::Matrix<iel,1>&               eprenp,
                         const LINALG::Matrix<nsd_,iel>&            eaccam,
                         const LINALG::Matrix<nsd_,iel>&            evelaf,
                         const LINALG::Matrix<nsd_,iel>&            fsevelaf,
                         Teuchos::RCP<const MAT::Material>          material,
                         const double                               alphaM,
                         const double                               alphaF,
                         const double                               gamma,
                         const double                               dt,
                         const double                               time,
                         const enum INPAR::FLUID::LinearisationAction     newton,
                         const bool                                 higher_order_ele,
                         const enum INPAR::FLUID::FineSubgridVisc         fssgv,
                         const enum INPAR::FLUID::Transient         inertia,
                         const enum INPAR::FLUID::PSPG              pspg,
                         const enum INPAR::FLUID::SUPG              supg,
                         const enum INPAR::FLUID::VStab             vstab,
                         const enum INPAR::FLUID::CStab             cstab,
                         const enum INPAR::FLUID::CrossStress       cross,
                         const enum INPAR::FLUID::ReynoldsStress    reynolds,
                         const enum INPAR::FLUID::TauType_genalpha                 whichtau,
                         const enum INPAR::FLUID::TurbModelAction         turb_mod_action,
                         double&                                    Cs,
                         double&                                    Cs_delta_sq,
                         double&                                    visceff,
                         const double                               l_tau,
                         const bool                                 compute_elemat
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
        \param fsevelaf           (i) nodal fine-scale velocity, time n+af
        \param material           (i) fluid material
        \param alphaM             (i) evaluate accelerations at this intermediate level
        \param alphaF             (i) evaluate velocities at this intermediate level
        \param gamma              (i) one step theta based parameter
        \param dt                 (i) time step size
        \param time               (i) current simulation time
        \param newton             (i) flag for use of full linearisation of reactive terms
        \param higher_order_ele   (i) include second derivatives in residual for stabilization
        \param fssgv              (i) flag for type of fine-scale subgrid viscosity
        \param tds                (i) flag for use of time dependent subscale stabilisation
        \param inertia            (i) flag for use of subscale accelerations for stabilisation (only for tds==true)
        \param pspg               (i) flag for use of nabla q testfunctions
        \param supg               (i) flag for use of (u o nabla) v testfunctions
        \param vstab              (i) flag for use of 2 * nu * grad o eps(v) testfunctions
        \param cross              (i) flag to add (u,(u_sub o nabla) v) to right hand side
        \param reynolds           (i) flag to add (u_sub,(u_sub o nabla) v) to right hand side
        \param whichtau           (i) flag choosing the tau version used for stabilization
        \param turb_mod_action    (i) action flag for the application of a turbulence model
        \param Cs               (i/o) value for Smagorinsky constant (zero if unused)
                                      in case of the dynamic model, it is only equal to the
                                      Smagorinsky constant on output.
        \param Cs_delta_sq        (i) Used for the dynamic Smagorinsky model. Corresponds to (Cs*hk)^2
        \param visceff            (o) effective viscosity computed by dynamic Smagorinsky model,
                                      non Newtonian fluid etc.
        \param l_tau              (i) viscous length scale for van Driest damping (zero if unused)
        \param compute_elemat     (i) if false, element matrices are not computed

      */
      void Sysmat_adv_td(Fluid*                                    ele,
                         std::vector<Epetra_SerialDenseVector>&     myknots,
                         LINALG::Matrix<(nsd_+1)*iel,(nsd_+1)*iel>& elemat,
                         LINALG::Matrix<(nsd_+1)*iel,    1>&        elevec,
                         const LINALG::Matrix<nsd_,iel>&            edispnp,
                         const LINALG::Matrix<nsd_,iel>&            egridvaf,
                         const LINALG::Matrix<nsd_,iel>&            evelnp,
                         const LINALG::Matrix<iel,1>&               eprenp,
                         const LINALG::Matrix<nsd_,iel>&            eaccam,
                         const LINALG::Matrix<nsd_,iel>&            evelaf,
                         const LINALG::Matrix<nsd_,iel>&            fsevelaf,
                         Teuchos::RCP<const MAT::Material>          material,
                         const double                               alphaM,
                         const double                               alphaF,
                         const double                               gamma,
                         const double                               dt,
                         const double                               time,
                         const enum INPAR::FLUID::LinearisationAction     newton,
                         const bool                                 higher_order_ele,
                         const enum INPAR::FLUID::FineSubgridVisc         fssgv,
                         const enum INPAR::FLUID::Transient         inertia,
                         const enum INPAR::FLUID::PSPG              pspg,
                         const enum INPAR::FLUID::SUPG              supg,
                         const enum INPAR::FLUID::VStab             vstab,
                         const enum INPAR::FLUID::CStab             cstab,
                         const enum INPAR::FLUID::CrossStress       cross,
                         const enum INPAR::FLUID::ReynoldsStress    reynolds,
                         const enum INPAR::FLUID::TauType_genalpha                 whichtau,
                         const enum INPAR::FLUID::TurbModelAction         turb_mod_action,
                         double&                                    Cs,
                         double&                                    Cs_delta_sq,
                         double&                                    visceff,
                         const double                               l_tau,
                         const bool                                 compute_elemat
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
        \param fsevelaf           (i) nodal fine-scale velocity, time n+af
        \param material           (i) fluid material
        \param alphaM             (i) evaluate accelerations at this intermediate level
        \param alphaF             (i) evaluate velocities at this intermediate level
        \param gamma              (i) one step theta based parameter
        \param dt                 (i) time step size
        \param time               (i) current simulation time
        \param newton             (i) flag for use of full linearisation of reactive terms
        \param higher_order_ele   (i) include second derivatives in residual for stabilization
        \param fssgv              (i) flag for type of fine-scale subgrid viscosity
        \param pspg               (i) flag for use of nabla q testfunctions
        \param supg               (i) flag for use of (u o nabla) v testfunctions
        \param vstab              (i) flag for use of 2 * nu * grad o eps(v) testfunctions
        \param cross              (i) flag to add (u,(u_sub o nabla) v) to right hand side
        \param reynolds           (i) flag to add (u_sub,(u_sub o nabla) v) to right hand side
        \param whichtau           (i) flag choosing the tau version used for stabilization
        \param turb_mod_action    (i) action flag for the application of a turbulence model
        \param Cs               (i/o) value for Smagorinsky constant (zero if unused)
                                      in case of the dynamic model, it is only equal to the
                                      Smagorinsky constant on output.
        \param Cs_delta_sq        (i) Used for the dynamic Smagorinsky model. Corresponds to (Cs*hk)^2
        \param visceff            (o) effective viscosity computed by dynamic Smagorinsky model,
                                      non Newtonian fluid etc.
        \param l_tau              (i) viscous length scale for van Driest damping (zero if unused)
        \param compute_elemat     (i) if false, element matrices are not computed

      */
      void Sysmat_cons_qs(Fluid*                                    ele,
                          std::vector<Epetra_SerialDenseVector>&     myknots,
                          LINALG::Matrix<(nsd_+1)*iel,(nsd_+1)*iel>& elemat,
                          LINALG::Matrix<(nsd_+1)*iel,    1>&        elevec,
                          const LINALG::Matrix<nsd_,iel>&            edispnp,
                          const LINALG::Matrix<nsd_,iel>&            egridvaf,
                          const LINALG::Matrix<nsd_,iel>&            evelnp,
                          const LINALG::Matrix<iel,1>&               eprenp,
                          const LINALG::Matrix<nsd_,iel>&            eaccam,
                          const LINALG::Matrix<nsd_,iel>&            evelaf,
                          const LINALG::Matrix<nsd_,iel>&            fsevelaf,
                          Teuchos::RCP<const MAT::Material>          material,
                          const double                               alphaM,
                          const double                               alphaF,
                          const double                               gamma,
                          const double                               dt,
                          const double                               time,
                          const enum INPAR::FLUID::LinearisationAction     newton,
                          const bool                                 higher_order_ele,
                          const enum INPAR::FLUID::FineSubgridVisc         fssgv,
                          const enum INPAR::FLUID::PSPG              pspg,
                          const enum INPAR::FLUID::SUPG              supg,
                          const enum INPAR::FLUID::VStab             vstab,
                          const enum INPAR::FLUID::CStab             cstab,
                          const enum INPAR::FLUID::CrossStress       cross,
                          const enum INPAR::FLUID::ReynoldsStress    reynolds,
                          const enum INPAR::FLUID::TauType_genalpha                 whichtau,
                          const enum INPAR::FLUID::TurbModelAction         turb_mod_action,
                          double&                                    Cs,
                          double&                                    Cs_delta_sq,
                          double&                                    visceff,
                          const double                               l_tau,
                          const bool                                 compute_elemat
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
        \param fsevelaf           (i) nodal fine-scale velocity, time n+af
        \param material           (i) fluid material
        \param alphaM             (i) evaluate accelerations at this intermediate level
        \param alphaF             (i) evaluate velocities at this intermediate level
        \param gamma              (i) one step theta based parameter
        \param dt                 (i) time step size
        \param time               (i) current simulation time
        \param newton             (i) flag for use of full linearisation of reactive terms
        \param higher_order_ele   (i) include second derivatives in residual for stabilization
        \param fssgv              (i) flag for type of fine-scale subgrid viscosity
        \param inertia            (i) flag for use of subscale accelerations for stabilisation (only for tds==true)
        \param pspg               (i) flag for use of nabla q testfunctions
        \param supg               (i) flag for use of (u o nabla) v testfunctions
        \param vstab              (i) flag for use of 2 * nu * grad o eps(v) testfunctions
        \param cross              (i) flag to add (u,(u_sub o nabla) v) to right hand side
        \param reynolds           (i) flag to add (u_sub,(u_sub o nabla) v) to right hand side
        \param whichtau           (i) flag choosing the tau version used for stabilization
        \param turb_mod_action    (i) action flag for the application of a turbulence model
        \param Cs               (i/o) value for Smagorinsky constant (zero if unused)
                                      in case of the dynamic model, it is only equal to the
                                      Smagorinsky constant on output.
        \param Cs_delta_sq        (i) Used for the dynamic Smagorinsky model. Corresponds to (Cs*hk)^2
        \param visceff            (o) effective viscosity computed by dynamic Smagorinsky model,
                                      non Newtonian fluid etc.
        \param l_tau              (i) viscous length scale for van Driest damping (zero if unused)
        \param compute_elemat     (i) if false, element matrices are not computed

      */
      void Sysmat_cons_td(Fluid*                                    ele,
                          std::vector<Epetra_SerialDenseVector>&     myknots,
                          LINALG::Matrix<(nsd_+1)*iel,(nsd_+1)*iel>& elemat,
                          LINALG::Matrix<(nsd_+1)*iel,    1>&        elevec,
                          const LINALG::Matrix<nsd_,iel>&            edispnp,
                          const LINALG::Matrix<nsd_,iel>&            egridvaf,
                          const LINALG::Matrix<nsd_,iel>&            evelnp,
                          const LINALG::Matrix<iel,1>&               eprenp,
                          const LINALG::Matrix<nsd_,iel>&            eaccam,
                          const LINALG::Matrix<nsd_,iel>&            evelaf,
                          const LINALG::Matrix<nsd_,iel>&            fsevelaf,
                          Teuchos::RCP<const MAT::Material>          material,
                          const double                               alphaM,
                          const double                               alphaF,
                          const double                               gamma,
                          const double                               dt,
                          const double                               time,
                          const enum INPAR::FLUID::LinearisationAction     newton,
                          const bool                                 higher_order_ele,
                          const enum INPAR::FLUID::FineSubgridVisc         fssgv,
                          const enum INPAR::FLUID::Transient         inertia,
                          const enum INPAR::FLUID::PSPG              pspg,
                          const enum INPAR::FLUID::SUPG              supg,
                          const enum INPAR::FLUID::VStab             vstab,
                          const enum INPAR::FLUID::CStab             cstab,
                          const enum INPAR::FLUID::CrossStress       cross,
                          const enum INPAR::FLUID::ReynoldsStress    reynolds,
                          const enum INPAR::FLUID::TauType_genalpha                 whichtau,
                          const enum INPAR::FLUID::TurbModelAction         turb_mod_action,
                          double&                                    Cs,
                          double&                                    Cs_delta_sq,
                          double&                                    visceff,
                          const double                               l_tau,
                          const bool                                 compute_elemat
        );


      void Sysmat_adv_qs_2D(
        Fluid*                                ele,
        std::vector<Epetra_SerialDenseVector>& myknots,
        LINALG::Matrix<(nsd_+1)*iel,(nsd_+1)*iel>&           elemat,
        LINALG::Matrix<(nsd_+1)*iel,    1>&           elevec,
        const LINALG::Matrix<nsd_,iel>&           edispnp,
        const LINALG::Matrix<nsd_,iel>&           egridvaf,
        const LINALG::Matrix<nsd_,iel>&           evelnp,
        const LINALG::Matrix<iel,1>&           eprenp,
        const LINALG::Matrix<nsd_,iel>&           eaccam,
        const LINALG::Matrix<nsd_,iel>&           evelaf,
        Teuchos::RCP<const MAT::Material>      material,
        const double                           alphaM,
        const double                           alphaF,
        const double                           gamma,
        const double                           dt,
        const double                           time,
        const enum INPAR::FLUID::LinearisationAction newton,
        const bool                             higher_order_ele,
        const enum INPAR::FLUID::Transient         inertia,
        const enum INPAR::FLUID::PSPG              pspg,
        const enum INPAR::FLUID::SUPG              supg,
        const enum INPAR::FLUID::VStab             vstab,
        const enum INPAR::FLUID::CStab             cstab,
        const enum INPAR::FLUID::CrossStress       cross,
        const enum INPAR::FLUID::ReynoldsStress    reynolds,
        const enum INPAR::FLUID::TauType_genalpha             whichtau,
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
      void Sysmat_adv_td_2D(
        Fluid*                                ele,
        std::vector<Epetra_SerialDenseVector>& myknots,
        LINALG::Matrix<(nsd_+1)*iel,(nsd_+1)*iel>&           elemat,
        LINALG::Matrix<(nsd_+1)*iel,    1>&           elevec,
        const LINALG::Matrix<nsd_,iel>&           edispnp,
        const LINALG::Matrix<nsd_,iel>&           egridvaf,
        const LINALG::Matrix<nsd_,iel>&           evelnp,
        const LINALG::Matrix<iel,1>&           eprenp,
        const LINALG::Matrix<nsd_,iel>&           eaccam,
        const LINALG::Matrix<nsd_,iel>&           evelaf,
        Teuchos::RCP<const MAT::Material>      material,
        const double                           alphaM,
        const double                           alphaF,
        const double                           gamma,
        const double                           dt,
        const double                           time,
        const enum INPAR::FLUID::LinearisationAction newton,
        const bool                             higher_order_ele,
        const enum INPAR::FLUID::Transient         inertia,
        const enum INPAR::FLUID::PSPG              pspg,
        const enum INPAR::FLUID::SUPG              supg,
        const enum INPAR::FLUID::VStab             vstab,
        const enum INPAR::FLUID::CStab             cstab,
        const enum INPAR::FLUID::CrossStress       cross,
        const enum INPAR::FLUID::ReynoldsStress    reynolds,
        const enum INPAR::FLUID::TauType_genalpha             whichtau,
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
      void Sysmat_cons_qs_2D(
        Fluid*                                ele,
        std::vector<Epetra_SerialDenseVector>& myknots,
        LINALG::Matrix<(nsd_+1)*iel,(nsd_+1)*iel>&           elemat,
        LINALG::Matrix<(nsd_+1)*iel,    1>&           elevec,
        const LINALG::Matrix<nsd_,iel>&           edispnp,
        const LINALG::Matrix<nsd_,iel>&           egridvaf,
        const LINALG::Matrix<nsd_,iel>&           evelnp,
        const LINALG::Matrix<iel,1>&           eprenp,
        const LINALG::Matrix<nsd_,iel>&           eaccam,
        const LINALG::Matrix<nsd_,iel>&           evelaf,
        Teuchos::RCP<const MAT::Material>      material,
        const double                           alphaM,
        const double                           alphaF,
        const double                           gamma,
        const double                           dt,
        const double                           time,
        const enum INPAR::FLUID::LinearisationAction newton,
        const bool                             higher_order_ele,
        const enum INPAR::FLUID::PSPG              pspg,
        const enum INPAR::FLUID::SUPG              supg,
        const enum INPAR::FLUID::VStab             vstab,
        const enum INPAR::FLUID::CStab             cstab,
        const enum INPAR::FLUID::CrossStress       cross,
        const enum INPAR::FLUID::ReynoldsStress    reynolds,
        const enum INPAR::FLUID::TauType_genalpha             whichtau,
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
      void Sysmat_cons_td_2D(
        Fluid*                                ele,
        std::vector<Epetra_SerialDenseVector>& myknots,
        LINALG::Matrix<(nsd_+1)*iel,(nsd_+1)*iel>&           elemat,
        LINALG::Matrix<(nsd_+1)*iel,    1>&           elevec,
        const LINALG::Matrix<nsd_,iel>&           edispnp,
        const LINALG::Matrix<nsd_,iel>&           egridvaf,
        const LINALG::Matrix<nsd_,iel>&           evelnp,
        const LINALG::Matrix<iel,1>&           eprenp,
        const LINALG::Matrix<nsd_,iel>&           eaccam,
        const LINALG::Matrix<nsd_,iel>&           evelaf,
        Teuchos::RCP<const MAT::Material>      material,
        const double                           alphaM,
        const double                           alphaF,
        const double                           gamma,
        const double                           dt,
        const double                           time,
        const enum INPAR::FLUID::LinearisationAction newton,
        const bool                             higher_order_ele,
        const enum INPAR::FLUID::Transient         inertia,
        const enum INPAR::FLUID::PSPG              pspg,
        const enum INPAR::FLUID::SUPG              supg,
        const enum INPAR::FLUID::VStab             vstab,
        const enum INPAR::FLUID::CStab             cstab,
        const enum INPAR::FLUID::CrossStress       cross,
        const enum INPAR::FLUID::ReynoldsStress    reynolds,
        const enum INPAR::FLUID::TauType_genalpha             whichtau,
        double&                                visceff,
        const bool                             compute_elemat
        );


      //@}

    private:

      //! @name service methods

      ///  Set parameters for turbulence models

      //
      //  \param ele              (i) the element those matrix is calculated
      //  \param turbmodelparams  (i) parameterlist for turbulence models
      //  \param fssgv            (i) Control parameter: use a fine-scale subgrid viscosity
      //  \param turb_mod_action  (o) enum specifying the turbulence approach used
      //  \param Cs               (o) Smagorinsky constant, read from parameterlist
      //  \param Cs_delta_sq      (o) Computed for dynamic Smagorinsky model
      //  \param l_tau            (o) viscous length scale of problem (set for van Driest
      //                              damping)
      //  \param nlayer           (o) number of the element layer in a turbulent channel
      //                              this one is used for a dynamic Smagorinsky model
      //                              with an average of MijLij etc. over channel planes
      //                              In addition we use it for the output of averaged
      //                              values of visceff, Cs_delta_sq, Cs

      //! number of components necessary to store second derivatives
      // 1 component  for nsd=1:  (N,xx)
      // 3 components for nsd=2:  (N,xx ; N,yy ; N,xy)
      // 6 components for nsd=3:  (N,xx ; N,yy ; N,zz ; N,xy ; N,xz ; N,yz)
      static const int numderiv2_ = DRT::UTILS::DisTypeToNumDeriv2<distype>::numderiv2;

      void SetParametersForTurbulenceModel(
        const Fluid*                       ele            ,
        ParameterList                     & turbmodelparams,
        ParameterList                     & sgvsicparams,
        const INPAR::FLUID::FineSubgridVisc     & fssgv          ,
        INPAR::FLUID::TurbModelAction           & turb_mod_action,
        double                            & Cs             ,
        double                            & Cs_delta_sq    ,
        double                            & l_tau          ,
        int                               & nlayer
        );

      ///  extract velocities, pressure and accelerations etc from the global distributed vectors

      //
      //  \param fssgv          (i) Control parameter: use a fine-scale subgrid viscosity
      //  \param is_ale         (i) Control parameter: get mesh/displacement
      //  \param discretization (i) discretization containing the global vectors
      //  \param lm             (i) location vector for element dofs in global vector
      //  \param eprenp         (o) nodal pressure values, t_{n+1}
      //  \param evelnp         (o) nodal velocities, t_{n+1}
      //  \param evelaf         (o) nodal velocities, t_{n+alpha_F}
      //  \param eaccam         (o) nodal accelerations, t_{n+alpha_M}
      //  \param edispnp        (o) nodal mesh displacements, t_{n+1} (for ALE)
      //  \param egridvelaf     (o) nodal mesh velocities, t_{n+alpha_F} (for ALE)
      //  \param fsevelaf       (o) fine-scale velocity

      void ExtractValuesFromGlobalVectors(
        const INPAR::FLUID::FineSubgridVisc fssgv         ,
        const bool                    is_ale        ,
        const DRT::Discretization&    discretization,
        const vector<int>&            lm,
        LINALG::Matrix<iel,1>&        eprenp        ,
        LINALG::Matrix<nsd_,iel>&     evelnp        ,
        LINALG::Matrix<nsd_,iel>&     evelaf        ,
        LINALG::Matrix<nsd_,iel>&     eaccam        ,
        LINALG::Matrix<nsd_,iel>&     edispnp       ,
        LINALG::Matrix<nsd_,iel>&     egridvelaf    ,
        LINALG::Matrix<nsd_,iel>&     fsevelaf
        );


      ///  Get the nodal values of the bodyforce. The values are stored in edeadnf
      //
      //  \param ele         (i) element pointer for access to the neumann condition
      //  \param time        (i) time when the condition is evaluated
      void GetNodalBodyForce(const Fluid* ele ,
                             const double  time);

      /*!  Compute the rate of strain of a (fine-scale) velocity

      //  \param evel       (i) nodal velocity values
      //  \param derxy      (i) shape function derivatives
      //  \param velderxy   (o) velocity derivatives
      //
      //  \return computed rate of strain

      */
      double GetStrainRate(
        const LINALG::Matrix<nsd_,iel>& evel     ,
        const LINALG::Matrix<nsd_,iel>& derxy    ,
        LINALG::Matrix<nsd_,nsd_>     & velderxy
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
          // j : direction of derivative x/y/z
          //
          for(int rr=0;rr<nsd_;++rr)
          {
            for(int mm=0;mm<nsd_;++mm)
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
          LINALG::Matrix<nsd_,nsd_> two_epsilon;

          for(int rr=0;rr<nsd_;++rr)
          {
            for(int mm=0;mm<nsd_;++mm)
            {
              two_epsilon(rr,mm) = velderxy(rr,mm) + velderxy(mm,rr);
            }
          }

          for(int rr=0;rr<nsd_;rr++)
          {
            for(int mm=0;mm<nsd_;mm++)
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
        const Fluid*                                 ele      ,
        const int                                  &  iquad    ,
        const DRT::UTILS::IntPointsAndWeights<nsd_>&  intpoints,
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
      void CalcTau(const enum INPAR::FLUID::TauType_genalpha             whichtau ,
                   const enum INPAR::FLUID::SubscalesTD           tds      ,
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
        \param fsevelaf                 (i) fine-scale velocity
    	\param myknots                  (i) knotvector for NURBS elements
        \param timealphaF               (i) alphaF*t_{n+1}+(1-alphaF)*t_{n}
        \param hk                       (o) element length for stabilisation parameter
        \param mk                       (o) element type constant mk for tau and the fssgv_artificial approach
    	\param material		 	(i) material
        \param visc                     (o) kinematicmaterial viscosity
        \param fssgv                    (i) Control parameter: use a fine-scale subgrid viscosity
        \param turb_mod_action          (i) enum specifying the turbulence approach used
        \param l_tau                    (i) viscous length scale of problem (set for van Driest damping)
        \param Cs                     (i/o) Smagorinsky constant, read from parameterlist for constant coefficient models,
                                            returned (computed value)for dynamic Smagorinsky model and van Driest damping
        \param Cs_delta_sq              (o) Computed for dynamic Smagorinsky model
        \param visceff        		(o) effective material viscosity

       */
      void SetElementData
      (
        Fluid*                                      ele            ,
        const LINALG::Matrix<nsd_,iel>             & edispnp        ,
        const LINALG::Matrix<nsd_,iel>             & evelaf         ,
        const LINALG::Matrix<nsd_,iel>             & fsevelaf       ,
        const std::vector<Epetra_SerialDenseVector>& myknots        ,
        const double                               & timealphaF     ,
        double                                     & hk             ,
        double                                     & mk             ,
        Teuchos::RCP<const MAT::Material>            material       ,
        double                                     & visc           ,
        const enum INPAR::FLUID::FineSubgridVisc           fssgv          ,
        const enum INPAR::FLUID::TurbModelAction           turb_mod_action,
        const double                                 l_tau          ,
        double                                     & Cs             ,
        double                                     & Cs_delta_sq    ,
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
        o fsvderxyaf_ (if fssgv)
        o viscs2_ (if hoel)
        o viscaf_old_ (if hoel)


     	\param ele               (i) pointer to element
        \param egridvelaf        (i) nodal mesh velocities, t_{n+alpha_F} (for ALE)
        \param evelnp            (i) nodal velocities, t_{n+1}
        \param eprenp            (i) nodal pressure values, t_{n+1}
        \param eaccam            (i) nodal accelerations, t_{n+alpha_M}
        \param evelaf            (i) nodal velocities, t_{n+alpha_F}
        \param fsevelaf          (i) fine-scale velocity
        \param visceff        	 (i) effective material viscosity
        \param fssgv             (i) Control parameter: use a fine-scale subgrid viscosity
        \param higher_order_ele  (i) include second derivatives in residual for stabilization

       */
      void InterpolateToGausspoint(
        Fluid*                                     ele             ,
        const LINALG::Matrix<nsd_,iel>            & egridvaf        ,
        const LINALG::Matrix<nsd_,iel>            & evelnp          ,
        const LINALG::Matrix<iel,1>               & eprenp          ,
        const LINALG::Matrix<nsd_,iel>            & eaccam          ,
        const LINALG::Matrix<nsd_,iel>            & evelaf          ,
        const LINALG::Matrix<nsd_,iel>            & fsevelaf        ,
        const double                              & visceff         ,
        const enum INPAR::FLUID::FineSubgridVisc          fssgv           ,
        const bool                                  higher_order_ele);


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

      void ExtractValuesFromGlobalVectors_2D(
        const bool                        is_ale        ,
        const DRT::Discretization&        discretization,
        const vector<int>&                lm            ,
        LINALG::Matrix<iel,1>&            eprenp        ,
        LINALG::Matrix<nsd_,iel>&            evelnp        ,
        LINALG::Matrix<nsd_,iel>&            evelaf        ,
        LINALG::Matrix<nsd_,iel>&            eaccam        ,
        LINALG::Matrix<nsd_,iel>&            edispnp       ,
        LINALG::Matrix<nsd_,iel>&            egridvelaf
        );


      ///  Get the nodal values of the bodyforce. The values are stored in edeadnf
      //
      //  \param ele         (i) element pointer for access to the neumann condition
      //  \param time        (i) time when the condition is evaluated
      void GetNodalBodyForce_2D(const Fluid* ele ,
                             const double  time);

      /*!  Compute the rate of strain of a (fine-scale) velocity

      //  \param evel       (i) nodal velocity values
      //  \param derxy      (i) shape function derivatives
      //  \param velderxy   (o) velocity derivatives
      //
      //  \return computed rate of strain

      */
      double GetStrainRate_2D(
        const LINALG::Matrix<nsd_,iel>& evel     ,
        const LINALG::Matrix<nsd_,iel>& derxy    ,
        LINALG::Matrix<nsd_,nsd_>        & velderxy
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
      double ShapeFunctionsFirstAndSecondDerivatives_2D(
        const Fluid*                                 ele      ,
        const int                                  &  iquad    ,
        const DRT::UTILS::IntPointsAndWeights<nsd_>&  intpoints,
        const std::vector<Epetra_SerialDenseVector>&  myknots  ,
        const bool                                    hoel
        );


      /*!
        \brief calculates viscosity at a given point (element center or Gaussian point)

      \param material     (i) material
        \param visc           (o) kinematic material viscosity
      \param rateofshear    (i) shear rate
       */
      void CalVisc_2D(Teuchos::RCP<const MAT::Material> material,
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
      void CalcTau_2D(const enum INPAR::FLUID::TauType_genalpha             whichtau ,
                   const enum INPAR::FLUID::SubscalesTD           tds      ,
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
      \param material     (i) material
        \param visc                     (o) kinematicmaterial viscosity
        \param visceff            (o) effective material viscosity

       */
      void SetElementData_2D
      (
        Fluid*                                      ele            ,
        const LINALG::Matrix<nsd_,iel>                & edispnp        ,
        const LINALG::Matrix<nsd_,iel>                & evelaf         ,
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
        \param visceff           (i) effective material viscosity
        \param higher_order_ele  (i) include second derivatives in residual for stabilization

       */
      void InterpolateToGausspoint_2D(
        Fluid*                                     ele             ,
        const LINALG::Matrix<nsd_,iel>               & egridvaf        ,
        const LINALG::Matrix<nsd_,iel>               & evelnp          ,
        const LINALG::Matrix<iel,1>               & eprenp          ,
        const LINALG::Matrix<nsd_,iel>               & eaccam          ,
        const LINALG::Matrix<nsd_,iel>               & evelaf          ,
        const double                              & visceff         ,
        const bool                                  higher_order_ele);


      //@}

      const int numdofpernode_;

      //! fine-scale subgrid viscosity
      double vart_;

      //! constant body force
      bool   constant_bodyforce_;

      // nodal arrays
      // ------------
      //! node coordinates
      LINALG::Matrix<nsd_,iel> xyze_;
      //! weights for nurbs elements
      LINALG::Matrix<iel,1> weights_;
      //! value of body force in the nodes of the element
      LINALG::Matrix<nsd_,iel> edeadaf_;

      // shape functions and derivatives, mapping from reference element to actual geometry
      // ----------------------------------------------------------------------------------
      //! vector of shape functions
      LINALG::Matrix<iel,1> funct_;
      //! vector of shape function derivatives in reference coordinate system
      LINALG::Matrix<nsd_,iel> deriv_;
      //! vector of second shape function derivatives in reference coordinate system
      LINALG::Matrix<numderiv2_,iel> deriv2_;
      //! vector of shape function derivatives in global coordinate system
      LINALG::Matrix<nsd_,iel> derxy_;
      //! vector of second shape function derivatives in global coordinate system
      LINALG::Matrix<numderiv2_,iel> derxy2_;
      //! derivatives needed for linearisation of viscous term including 2nd derivatives
      /*
                       corresponds to a part of the operator

                                  grad o epsilon

                       on the discret level
      */
      LINALG::Matrix<nsd_,iel> viscs2_;
      //! transpose of the jacobian matrix of the mapping (r,s,t)->(x,y,z)
      LINALG::Matrix<nsd_,nsd_>   xjm_;
      //! its inverse
      LINALG::Matrix<nsd_,nsd_>   xji_;
      //! second derivatives of the mapping (r,s,t)->(x,y,z)
      LINALG::Matrix<numderiv2_,nsd_>   xder2_;

      // values of non geometrical quantities in gausspoints
      // ---------------------------------------------------
      //! pressure in gausspoint, time n+1
      double                   prenp_;
      //! acceleration in gausspoint, time n+am
      LINALG::Matrix<nsd_,1>   accintam_;
      //! velocity in gausspoint, time n+1
      LINALG::Matrix<nsd_,1>   velintnp_;
      //! velocity in gausspoint, time n+af
      LINALG::Matrix<nsd_,1>   velintaf_;
      // this copy will be used to store the normed velocity at time n+af
      LINALG::Matrix<nsd_,1>   normed_velintaf_;
      //! fine-scale velocity in gausspoint, time n+af
      LINALG::Matrix<nsd_,1>   fsvelintaf_;
      //! grid velocity at time level n+af
      LINALG::Matrix<nsd_,1>   u_G_af_;
      //! ale convective velocity at intermediate time level n+af
      //! aleconvint_af_ = velintaf_ - gridvelintaf_
      LINALG::Matrix<nsd_,1>   aleconvintaf_;
      //! pressure gradient in gausspoint, time n+1
      LINALG::Matrix<nsd_,1>   pderxynp_;
      //! velocity derivatives in gausspoint, time n+1
      LINALG::Matrix<nsd_,nsd_>   vderxynp_;
      //! velocity derivatives in gausspoint, time n+af
      LINALG::Matrix<nsd_,nsd_>   vderxyaf_;
      //! fine-scale velocity derivatives in gausspoint, time n+af
      LINALG::Matrix<nsd_,nsd_>   fsvderxyaf_;
      //! value of body force interpolated to the gauss point, time n+af
      LINALG::Matrix<nsd_,1>   bodyforceaf_;
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
      LINALG::Matrix<nsd_,1>   svelaf_;
      //! matrix of subscale velocities in gausspoints of this element, time n+1
      LINALG::Matrix<nsd_,1>   svelnp_;
      //! Convective term, most recent iteration value
      LINALG::Matrix<nsd_,1>   convaf_old_;
      //! Convectivion by mesh motion, most recent iteration value
      LINALG::Matrix<nsd_,1>   convu_G_af_old_;
      //! Convective part of scale-similarity term, most recent iteration value
      LINALG::Matrix<nsd_,1>   convaf_s_;
      /*
                         /  n+af        \   n+af
                        |  u     o grad  | u
                         \  (i)         /   (i)
      */
      //! Subgrid velocity convective term, most recent iteration value
      LINALG::Matrix<nsd_,1>   convsubaf_old_;
      /* Most recent value for subgrid velocity convective term

                        /~n+af         \   n+af
                       | u      o nabla | u
                        \   (i)        /   (i)
      */
      //! Viscous term, most recent iteration value
      LINALG::Matrix<nsd_,1>   viscaf_old_;
      /*
                                        /  n+af \
                       nabla o  epsilon | u     |
                                        \  (i)  /
      */
      //! Current residual of the momentum equation
      LINALG::Matrix<nsd_,1>   resM_;

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

      //! the element volume, used to define element length
      double                vol_;

      //! Jacobian determinant
      double                det_;

      //! velocity derivatives w.r.t. r,s,t.
      LINALG::Matrix<nsd_,nsd_> vderiv_;

      // temporary working arrays
      // ------------------------
      //! array used for the computation of jacobian
      //! determinant and second derivatives
      LINALG::Matrix<numderiv2_,numderiv2_>   bm_;

      //! An epetra wrapper for LAPACK library functions
      Epetra_LAPACK         solver_;

    };
  }
}

#endif

