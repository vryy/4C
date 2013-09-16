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
Maintainers: Ursula Rasthofer & Volker Gravemeier
             {rasthofer,vgravem}@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289-15236/-245
</pre>
*/
/*----------------------------------------------------------------------*/

#ifndef FLUID_ELE_BOUNDARY_PARENT_CALC_H
#define FLUID_ELE_BOUNDARY_PARENT_CALC_H


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

    /// Interface base class for FluidBoundaryParent
    /*!
      This class exists to provide a common interface for all template
      versions of FluidBoundaryParent. The only function
      this class actually defines is Impl, which returns a pointer to
      the appropriate version of FluidBoundaryImpl.
    */
    class FluidBoundaryParentInterface : DRT::SingletonDestruction
    {
    public:

      //! empty constructor
      FluidBoundaryParentInterface() {}

      //! empty destructor
      virtual ~FluidBoundaryParentInterface() {}

      virtual void FlowDepPressureBC(
        DRT::ELEMENTS::FluidBoundary*  ele1,
        Teuchos::ParameterList&        params,
        DRT::Discretization&           discretization,
        std::vector<int>&              lm,
        Epetra_SerialDenseMatrix&      elemat1,
        Epetra_SerialDenseVector&      elevec1) = 0;

      virtual void EvaluateWeakDBC(
        DRT::ELEMENTS::FluidBoundary*  ele1,
        Teuchos::ParameterList&        params,
        DRT::Discretization&           discretization,
        std::vector<int>&              lm,
        Epetra_SerialDenseMatrix&      elemat1,
        Epetra_SerialDenseVector&      elevec1) = 0;

      virtual void EvaluateNitschePar(
        DRT::ELEMENTS::FluidBoundary*  ele1,
        Teuchos::ParameterList&        params,
        DRT::Discretization&           discretization,
        std::vector<int>&              lm,
        Epetra_SerialDenseMatrix&      elemat1,
        Epetra_SerialDenseMatrix&      elemat2) = 0;

      virtual void MixHybDirichlet(
        DRT::ELEMENTS::FluidBoundary*  ele1,
        Teuchos::ParameterList&        params,
        DRT::Discretization&           discretization,
        std::vector<int>&              lm,
        Epetra_SerialDenseMatrix&      elemat1,
        Epetra_SerialDenseVector&      elevec1) = 0;

      /// Internal implementation class for FluidBoundaryParent elements
      static FluidBoundaryParentInterface* Impl(DRT::ELEMENTS::FluidBoundary* ele);

    protected:
      static std::map<int,std::map<int,FluidBoundaryParentInterface*> * > instances_;

    };


  /// Internal FluidBoundaryParent element implementation
  /*!
    This internal class keeps all the working arrays needed to
    calculate the FluidBoundaryParent element. Additionally the method Sysmat()
    provides a clean and fast element implementation.

  */
    template<DRT::Element::DiscretizationType distype>
    class FluidBoundaryParent: public FluidBoundaryParentInterface
    {
      friend class FluidEleParameter;
      typedef FluidBoundaryParentInterface my;

    public:

      //! Singleton access method
      static FluidBoundaryParentInterface* Instance(bool create=true, int num=0);

      /// called upon destruction
      virtual void Done();

      /// Constructor with number of nodes
      FluidBoundaryParent(int num);

      void  FlowDepPressureBC(
        DRT::ELEMENTS::FluidBoundary*  surfele,
        Teuchos::ParameterList&        params,
        DRT::Discretization&           discretization,
        std::vector<int>&              lm,
        Epetra_SerialDenseMatrix&      elemat,
        Epetra_SerialDenseVector&      elevec);

      void EvaluateWeakDBC(
        DRT::ELEMENTS::FluidBoundary*  surfele,
        Teuchos::ParameterList&        params,
        DRT::Discretization&           discretization,
        std::vector<int>&              lm,
        Epetra_SerialDenseMatrix&      elemat,
        Epetra_SerialDenseVector&      elevec);

      void EvaluateNitschePar(
        DRT::ELEMENTS::FluidBoundary*  surfele,
        Teuchos::ParameterList&        params,
        DRT::Discretization&           discretization,
        std::vector<int>&              lm,
        Epetra_SerialDenseMatrix&      elemat1,
        Epetra_SerialDenseMatrix&      elemat2);

      void  MixHybDirichlet(
        DRT::ELEMENTS::FluidBoundary*  surfele,
        Teuchos::ParameterList&        params,
        DRT::Discretization&           discretization,
        std::vector<int>&              lm,
        Epetra_SerialDenseMatrix&      elemat,
        Epetra_SerialDenseVector&      elevec);


    private:

      template <DRT::Element::DiscretizationType bdistype,
                DRT::Element::DiscretizationType pdistype>
      void FlowDepPressureBC(
        DRT::ELEMENTS::FluidBoundary*  ele,
        Teuchos::ParameterList&        params,
        DRT::Discretization&           discretization,
        std::vector<int>&              lm,
        Epetra_SerialDenseMatrix&      elemat,
        Epetra_SerialDenseVector&      elevec);

      template <DRT::Element::DiscretizationType bdistype,
                DRT::Element::DiscretizationType pdistype>
      void EvaluateWeakDBC(
        DRT::ELEMENTS::FluidBoundary*  ele,
        Teuchos::ParameterList&        params,
        DRT::Discretization&           discretization,
        std::vector<int>&              lm,
        Epetra_SerialDenseMatrix&      elemat,
        Epetra_SerialDenseVector&      elevec);

      template <DRT::Element::DiscretizationType bdistype,
                DRT::Element::DiscretizationType pdistype>
      void EvaluateNitschePar(
        DRT::ELEMENTS::FluidBoundary*  ele,
        Teuchos::ParameterList&        params,
        DRT::Discretization&           discretization,
        std::vector<int>&              lm,
        Epetra_SerialDenseMatrix&      elemat1,
        Epetra_SerialDenseMatrix&      elemat2);

      template <DRT::Element::DiscretizationType bdistype,
                DRT::Element::DiscretizationType pdistype>
      void MixHybDirichlet(
        DRT::ELEMENTS::FluidBoundary*  ele,
        Teuchos::ParameterList&        params,
        DRT::Discretization&           discretization,
        std::vector<int>&              lm,
        Epetra_SerialDenseMatrix&      elemat,
        Epetra_SerialDenseVector&      elevec);

      void GetDensityAndViscosity(Teuchos::RCP<const MAT::Material> material,
                                  const double                      rateofstrain);


      //! pointer to parameter list
      Teuchos::RCP<DRT::ELEMENTS::FluidEleParameter> fldpara_;

      //! infinitesimal area element drs
      double drs_;
      //! integration factor
      double fac_;
      //! physical viscosity
      double visc_;
      //! density at t_(n+alpha_F) or t_(n+1)
      double densaf_;

      //-----------------------------------------------------------------
      //-----------------------------------------------------------------
      //
      //                    SPALDINGS LAW OF THE WALL
      //
      //-----------------------------------------------------------------
      //-----------------------------------------------------------------
      //-----------------------------------------------------------------
      //  evaluate the residual of Spaldings law of the wall
      //                                           (private) gammi 11/09
      //-----------------------------------------------------------------
      double SpaldingResidual(const double y,
                              const double visc,
                              const double tau_B,
                              const double normu)
      {
        // get dimensionless velocity
        const double up = sqrt(normu/tau_B);

        // constants
        const double chi=0.4;
        const double B  =5.5;

        //      +
        // get y , a dimensionless boundary layer thickness
        const double yp = y*sqrt(tau_B*normu)/visc;

        return(yp-(up+exp(-chi*B)*(exp(chi*up)-1.0-chi*up*(1+chi*up/2.0*(1.0+chi*up/3.0)))));
      }

      //-----------------------------------------------------------------
      //  evaluate the residual of Spaldings law of the wall
      //                                           (private) gammi 11/09
      //-----------------------------------------------------------------
      double SpaldingResidual_utau(const double y,
                                   const double visc,
                                   const double utau,
                                   const double normu)
      {
        // get dimensionless velocity
        if (abs(utau)<1.0E-14) dserror("utau is zero!");
        const double up = normu/utau;

        // constants
        const double chi=0.4;
        const double B  =5.5;

        //      +
        // get y , a dimensionless boundary layer thickness
        if (visc<1.0E-14) dserror("visc is zero or negative!");
        const double yp = y*utau/visc;

        return(yp-(up+exp(-chi*B)*(exp(chi*up)-1.0-chi*up*(1+chi*up/2.0*(1.0+chi*up/3.0)))));
      }

      //-----------------------------------------------------------------
      //     evaluate the residual of Spaldings law of the wall
      //                                           (private) gammi 11/09
      //-----------------------------------------------------------------
      double JacobianSpaldingResidual(const double y,
                                      const double visc,
                                      const double tau_B,
                                      const double normu)
      {
        // constants
        const double chi=0.4;
        const double B  =5.5;

        // get dimensionless velocity
        const double up = sqrt(normu/tau_B);

        // compute the derivative of the Spalding residual w.r.t. tau_B
        double drdtauB = y/(2.0*visc*sqrt(tau_B))*sqrt(normu);

        drdtauB += (1+chi*exp(-chi*B)*(exp(chi*up)-1.0-chi*up*(1.0+0.5*chi*up)))
             *0.5*sqrt(normu)/(sqrt(tau_B)*sqrt(tau_B)*sqrt(tau_B));

        return(drdtauB);
      }

      //-----------------------------------------------------------------
      //     evaluate the residual of Spaldings law of the wall
      //                                           (private) gammi 11/09
      //-----------------------------------------------------------------
      double JacobianSpaldingResidual_utau(const double y,
                                           const double visc,
                                           const double utau,
                                           const double normu)
      {
        // get dimensionless velocity
        const double up = normu/utau;

        // constants
        const double chi=0.4;
        const double B  =5.5;

        //                    +
        // get derivative of y , a dimensionless boundary layer thickness
        const double dyplus_dutau=y/visc;

        //                                   +
        // get derivative of function f wrt u
        const double df_duplus=1+exp(-chi*B)*(chi*exp(chi*up)-chi-chi*chi*up-chi*chi*chi*up*up/2.0);

        //                    +
        // get derivative of u  wrt u_tau
        const double duplus_dutau=-normu/(utau*utau);

        return(dyplus_dutau-df_duplus*duplus_dutau);
      }

      //-----------------------------------------------------------------
      //     evaluate the residual of Spaldings law of the wall
      //                                           (private) gammi 11/09
      //-----------------------------------------------------------------
      double JacobianSpaldingResidual_u(const double y,
                                        const double visc,
                                        const double utau,
                                        const double normu)
      {
        // get dimensionless velocity
        const double up = normu/utau;

        // constants
        const double chi=0.4;
        const double B  =5.5;

        //                    +
        // get derivative of y , a dimensionless boundary layer thickness
        const double dyplus_du=0.0;

        //                                   +
        // get derivative of function f wrt u
        const double df_duplus=1+exp(-chi*B)*(chi*exp(chi*up)-chi-chi*chi*up-chi*chi*chi*up*up/2.0);

        //                    +
        // get derivative of u  wrt u
        const double duplus_du=1.0/utau;

        return(dyplus_du-df_duplus*duplus_du);
      }

      //-----------------------------------------------------------------
      //     evaluate the residual of Spaldings law of the wall
      //                                           (private) gammi 11/09
      //-----------------------------------------------------------------
      double JacobianSpaldingResidual_uplus(const double y,
                                            const double visc,
                                            const double utau,
                                            const double normu)
      {
        // get dimensionless velocity
        const double up = normu/utau;

        // constants
        const double chi=0.4;
        const double B  =5.5;

        //                    +
        // get derivative of y , a dimensionless boundary layer thickness
        const double dyplus_duplus=0.0;
        /*
        // Spaldings law of the wall
        //                     /                                                  \
        //                    |                           /     +\ 2    /     +\ 3 |
        //                    |       +                  | chi*u  |    | chi*u  |  |
        //  +   +    -chi*B   |  chi*u               +    \      /      \      /   |
        // y - u  - e       * | e       - 1.0 - chi*u  - ----------- - ----------- | = 0
        //                    |                              2.0           6.0     |
        //                    |                                                    |
        //                     \                                                  /
        */

        //                                   +
        // get derivative of function f wrt u
        const double df_duplus=1+exp(-chi*B)*(chi*exp(chi*up)-chi-chi*chi*up-chi*chi*chi*up*up/2.0);

        return(dyplus_duplus-df_duplus);
      }
    };  // end class fluid3boundaryImpl
  }
}

#endif
