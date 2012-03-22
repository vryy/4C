/*----------------------------------------------------------------------*/
/*!
\file fluid_ele_calc.cpp

\brief Internal implementation of Fluid3 element

<pre>
Maintainer: Shadan Shahmiri /Benedikt Schott
            shahmiri@lnm.mw.tum.de
            schott@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15240
</pre>
*/
/*----------------------------------------------------------------------*/
#include <Teuchos_TimeMonitor.hpp>

#include <fstream>

#include "fluid_ele_calc_xfem.H"
#include "fluid_ele.H"
#include "fluid_ele_parameter.H"
// TODO: Benedikt: Einbinden der hk über ele_utils.H, nicht fluid_ele_calc_stabilization?!?!
#include "fluid_ele_calc_stabilization.H"
//#include "fluid_ele_utils.H"

#include "../drt_bele3/bele3.H"
#include "../drt_bele3/bele3_4.H"

#include "../drt_cut/cut_boundarycell.H"
#include "../drt_cut/cut_position.H"

#include "../drt_geometry/position_array.H"

#include "../drt_inpar/inpar_xfem.H"

#include "../drt_lib/drt_elementtype.H"
#include "../drt_lib/drt_element.H"

#include "../linalg/linalg_fixedsizeblockmatrix.H"
#include "../linalg/linalg_utils.H"


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::FluidEleCalcXFEM<distype> * DRT::ELEMENTS::FluidEleCalcXFEM<distype>::Instance( bool create )
{
  static FluidEleCalcXFEM<distype> * instance;
  if ( create )
  {
    if ( instance==NULL )
    {
      instance = new FluidEleCalcXFEM<distype>();
    }
  }
  else
  {
    if ( instance!=NULL )
      delete instance;
    instance = NULL;
  }
  return instance;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidEleCalcXFEM<distype>::Done()
{
  // delete this pointer! Afterwards we have to go! But since this is a
  // cleanup call, we can do it this way.
  Instance( false );
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::FluidEleCalcXFEM<distype>::FluidEleCalcXFEM()
  : DRT::ELEMENTS::FluidEleCalc<distype>::FluidEleCalc()
{

}


// TODO: Viel Spass euch beiden, Shadan und Benedikt!!!!!
namespace DRT
{
  namespace ELEMENTS
  {
    namespace XFLUID
    {

      template<DRT::Element::DiscretizationType distype>
      class SideInterface
      {
      public:


        static const int nen_ = DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;
        static const int nsd_ = DRT::UTILS::DisTypeToDim<distype>::dim;

        static Teuchos::RCP<SideInterface<distype> > Impl(DRT::Element * side,
                                                          Epetra_SerialDenseMatrix & C_uiu,
                                                          Epetra_SerialDenseMatrix & C_uui,
                                                          Epetra_SerialDenseMatrix & rhC_ui,
                                                          Epetra_SerialDenseMatrix & Gsui,
                                                          Epetra_SerialDenseMatrix & Guis,
                                                          Epetra_SerialDenseMatrix & side_xyze
                                                          );

        static Teuchos::RCP<SideInterface<distype> > Impl(DRT::Element * side,
                                                          Epetra_SerialDenseMatrix & C_uiu,
                                                          Epetra_SerialDenseMatrix & C_uui,
                                                          Epetra_SerialDenseMatrix & rhC_ui,
                                                          Epetra_SerialDenseMatrix & C_uiui,
                                                          Epetra_SerialDenseMatrix & side_xyze
                                                          );

        static Teuchos::RCP<SideInterface<distype> > Impl(DRT::Element * side,
                                                          Epetra_SerialDenseMatrix & side_xyze
                                                          );

        virtual ~SideInterface() {}

        virtual void Evaluate(const LINALG::Matrix<2,1> & eta,
                              LINALG::Matrix<3,1> & x,
                              LINALG::Matrix<3,1> & normal,
                              double & drs) = 0;


        virtual void ProjectOnSide(LINALG::Matrix<3,1> & x_gp_lin,
                                   LINALG::Matrix<3,1> & x_side,
                                   LINALG::Matrix<2,1> & xi_side) = 0;


        virtual void eivel(const DRT::Discretization &  cutdis,
                           const std::string            state,
                           const vector<int>&           lm) = 0;


        virtual void addeidisp(const DRT::Discretization &  cutdis,
                               const std::string            state,
                               const vector<int>&           lm,
                               Epetra_SerialDenseMatrix  &  side_xyze) = 0;


        virtual void buildInterfaceForce( const Teuchos::RCP<Epetra_Vector> &  iforcecol,
                                          const DRT::Discretization &          cutdis,
                                          const vector<int>&                   lm,
                                          LINALG::Matrix<nsd_,1> &             traction,
                                          double &                             fac ) = 0;


        virtual void buildCouplingMatrices(LINALG::Matrix<3,1> & normal,
                                           const double          fac,
                                           LINALG::Matrix<nen_,1>  & funct,
                                           LINALG::BlockMatrix<LINALG::Matrix<nen_,  1>,6,1> &  rhs
                                           ) = 0;

        virtual void buildFinalCouplingMatrices(LINALG::BlockMatrix<LINALG::Matrix<nen_,nen_>,6,6> &  invK_ss,
                                                LINALG::BlockMatrix<LINALG::Matrix<nen_,nen_>,4,6> &  K_iK,
                                                LINALG::BlockMatrix<LINALG::Matrix<nen_,nen_>,6,4> &  K_su,
                                                LINALG::BlockMatrix<LINALG::Matrix<nen_,  1>,6,1> &  rhs
                                                ) = 0;

        virtual void buildCouplingMatricesNitsche(   Epetra_SerialDenseMatrix &    C_uu_,          // standard bg-bg-matrix
                                                     Epetra_SerialDenseVector &    rhs_Cu_,        // standard bg-rhs
                                                     bool &                        coupling,       // assemble coupling terms (yes/no)
                                                     bool &                        bg_mortaring,   // yes: background-sided mortaring, no: coupling between two meshes (mixed mortaring)
                                                     LINALG::Matrix<nsd_,1> &      normal,         // normal vector
                                                     const double                  timefacfac,     // theta*dt
                                                     const double                  visceff_1,      // viscosity in background fluid
                                                     const double                  visceff_2,      // viscosity in embedded fluid
                                                     double &                      kappa1,         // mortaring weighting
                                                     double &                      kappa2,         // mortaring weighting
                                                     double &                      stabfac,        // Nitsche non-dimensionless stabilization factor
                                                     double &                      stabfac_conv,   // Nitsche convective non-dimensionless stabilization factor
                                                     LINALG::Matrix<nen_,1> &      funct_,         // bg shape functions
                                                     LINALG::Matrix<nsd_,nen_> &   derxy_,         // bg deriv
                                                     LINALG::Matrix<nsd_,nsd_> &   vderxy_,        // bg deriv^n
                                                     double &                      press,          // bg p^n
                                                     LINALG::Matrix<nsd_,1> &      velint,          // bg u^n
                                                     LINALG::Matrix<nsd_,1> &      ivelint_WDBC_JUMP // Dirichlet velocity vector or prescribed jump vector
                                                  ) = 0;

        virtual void get_vel_WeakDBC (LINALG::Matrix<nsd_,1> & ivelint) = 0;

      };

      template<DRT::Element::DiscretizationType distype,
               DRT::Element::DiscretizationType side_distype,
               const int numdof>

      class SideImpl : public SideInterface<distype>
      {
      public:

        //! nen_: number of element nodes (P. Hughes: The Finite Element Method)
        static const int side_nen_ = DRT::UTILS::DisTypeToNumNodePerEle<side_distype>::numNodePerElement;
        static const int nen_ = DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;
        static const int nsd_ = DRT::UTILS::DisTypeToDim<distype>::dim;

        // for stress-based fluid-fluid coupling
        SideImpl(DRT::Element * side,
                 Epetra_SerialDenseMatrix & C_uiu,
                 Epetra_SerialDenseMatrix & C_uui,
                 Epetra_SerialDenseMatrix & rhC_ui,
                 Epetra_SerialDenseMatrix & Gsui,
                 Epetra_SerialDenseMatrix & Guis,
                 Epetra_SerialDenseMatrix & side_xyze
                 )
          : C_uiu_(C_uiu.A(),true),
            C_uui_(C_uui.A(),true),
            rhC_ui_(rhC_ui.A(),true),
            K_sui_(Gsui.A(),true),
            K_uis_(Guis.A(),true),
            xyze_(side_xyze.A(),true)
        {
        }

        // for Nitsche-based fluid-fluid coupling
        SideImpl(DRT::Element * side,
                Epetra_SerialDenseMatrix & C_uiu,
                Epetra_SerialDenseMatrix & C_uui,
                Epetra_SerialDenseMatrix & rhC_ui,
                Epetra_SerialDenseMatrix & C_uiui,
                Epetra_SerialDenseMatrix & side_xyze
                 )
          : C_uiu_(C_uiu.A(),true),
            C_uui_(C_uui.A(),true),
            rhC_ui_(rhC_ui.A(),true),
            C_uiui_(C_uiui.A(),true),
            xyze_(side_xyze.A(),true)
        {
        }

        // without any coupling
        SideImpl(DRT::Element * side,
                 Epetra_SerialDenseMatrix & side_xyze
                 )
          : xyze_(side_xyze.A(),true)
        {
        }

        virtual void Evaluate(const LINALG::Matrix<nsd_-1,1>  & eta,
                              LINALG::Matrix<nsd_,1>          & x,
                              LINALG::Matrix<nsd_,1>          & normal,
                              double                          & drs
                              )
        {
          LINALG::Matrix<nsd_-1,2> metrictensor;
          DRT::UTILS::shape_function_2D( side_funct_, eta( 0 ), eta( 1 ), side_distype );
          DRT::UTILS::shape_function_2D_deriv1( side_deriv_, eta( 0 ), eta( 1 ), side_distype );
          DRT::UTILS::ComputeMetricTensorForBoundaryEle<side_distype>( xyze_,side_deriv_, metrictensor, drs, &normal );
          x.Multiply( xyze_,side_funct_ );

        }


        virtual void ProjectOnSide( LINALG::Matrix<3,1> & x_gp_lin,
                                    LINALG::Matrix<3,1> & x_side,
                                    LINALG::Matrix<2,1> & xi_side
                                  )
        {
          TEUCHOS_FUNC_TIME_MONITOR( "FLD::XFluid::XFluidState::ProjectOnSide" );
          // Initialization
          LINALG::Matrix<side_nen_,1> funct(true);          // shape functions
          LINALG::Matrix<2,side_nen_> deriv(true);          // derivatives dr, ds
          LINALG::Matrix<3,side_nen_> deriv2(true);          // 2nd derivatives drdr, dsds, drds


          LINALG::Matrix<3,1> x(true);

          LINALG::Matrix<3,2> derxy (true);
          LINALG::Matrix<3,1> dx_dr (true);
          LINALG::Matrix<3,1> dx_ds (true);

          LINALG::Matrix<3,3> derxy2 (true);
          LINALG::Matrix<3,1> dx_drdr (true);
          LINALG::Matrix<3,1> dx_dsds (true);
          LINALG::Matrix<3,1> dx_drds (true);

          LINALG::Matrix<3,1> residuum(true);             // residuum of the newton iteration
          LINALG::Matrix<3,3> sysmat(true);               // matrix for the newton system
          LINALG::Matrix<3,1> incr(true);                 // increment of the newton system

          LINALG::Matrix<3,1> sol(true); // sol carries xi_1, xi_2, d (distance)

          if(side_distype == DRT::Element::tri3 or
              side_distype == DRT::Element::tri6)
          {
            sol(0) = 0.333333333333333;
            sol(1) = 0.333333333333333;
          }
          else if( side_distype == DRT::Element::quad4 or
              side_distype == DRT::Element::quad8 or
              side_distype == DRT::Element::quad9)
          {
            sol(0) = 0.0;
            sol(1) = 0.0;
          }
          else
          {
            dserror("define start side xi-coordinates for unsupported cell type");
          }

          const double absTolIncr = 1.0e-9;   // rel tolerance for the local coordinates increment
          const double relTolRes  = 1.0e-9;   // rel tolerance for the whole residual
          const double absTOLdist = 1.0e-9;   // abs tolerance for distance

          int iter=0;
          const int maxiter = 7;

          bool converged = false;

          while(iter < maxiter && !converged)
          {

            iter++;


            // get current values
            DRT::UTILS::shape_function_2D( funct, sol( 0 ), sol( 1 ), side_distype );
            DRT::UTILS::shape_function_2D_deriv1( deriv, sol( 0 ), sol( 1 ), side_distype );
            DRT::UTILS::shape_function_2D_deriv2( deriv2, sol( 0 ), sol( 1 ), side_distype );

            x.Multiply(xyze_, funct);

            derxy.MultiplyNT(xyze_, deriv);

            derxy2.MultiplyNT(xyze_, deriv2);

            // set dx_dr and dx_ds
            for (int i=0; i< 3; i++)
            {
              dx_dr(i) = derxy(i,0);
              dx_ds(i) = derxy(i,1);

              dx_drdr(i) = derxy2(i,0);
              dx_dsds(i) = derxy2(i,1);
              dx_drds(i) = derxy2(i,2);
            }

            // get vector products
            LINALG::Matrix<3,1> dx_drdr_times_dx_ds(true);
            LINALG::Matrix<3,1> dx_dr_times_dx_drds(true);
            LINALG::Matrix<3,1> dx_drds_times_dx_ds(true);
            LINALG::Matrix<3,1> dx_dr_times_dx_dsds(true);
            LINALG::Matrix<3,1> dx_dr_times_dx_ds(true);

            dx_drdr_times_dx_ds(0) = dx_drdr(1)*dx_ds(2)-dx_ds(1)*dx_drdr(2);
            dx_drdr_times_dx_ds(1) = dx_drdr(2)*dx_ds(0)-dx_ds(2)*dx_drdr(0);
            dx_drdr_times_dx_ds(2) = dx_drdr(0)*dx_ds(1)-dx_ds(0)*dx_drdr(1);

            dx_dr_times_dx_drds(0) = dx_dr(1)*dx_drds(2)-dx_drds(1)*dx_dr(2);
            dx_dr_times_dx_drds(1) = dx_dr(2)*dx_drds(0)-dx_drds(2)*dx_dr(0);
            dx_dr_times_dx_drds(2) = dx_dr(0)*dx_drds(1)-dx_drds(0)*dx_dr(1);

            dx_drds_times_dx_ds(0) = dx_drds(1)*dx_ds(2)-dx_ds(1)*dx_drds(2);
            dx_drds_times_dx_ds(1) = dx_drds(2)*dx_ds(0)-dx_ds(2)*dx_drds(0);
            dx_drds_times_dx_ds(2) = dx_drds(0)*dx_ds(1)-dx_ds(0)*dx_drds(1);

            dx_dr_times_dx_dsds(0) = dx_dr(1)*dx_dsds(2)-dx_dsds(1)*dx_dr(2);
            dx_dr_times_dx_dsds(1) = dx_dr(2)*dx_dsds(0)-dx_dsds(2)*dx_dr(0);
            dx_dr_times_dx_dsds(2) = dx_dr(0)*dx_dsds(1)-dx_dsds(0)*dx_dr(1);

            dx_dr_times_dx_ds(0) = dx_dr(1)*dx_ds(2)-dx_ds(1)*dx_dr(2);
            dx_dr_times_dx_ds(1) = dx_dr(2)*dx_ds(0)-dx_ds(2)*dx_dr(0);
            dx_dr_times_dx_ds(2) = dx_dr(0)*dx_ds(1)-dx_ds(0)*dx_dr(1);

            // define sysmat
            for(int i=0; i< 3; i++)
            {
              // d/dr
              sysmat(i,0) = dx_dr(i) - sol(2) * (dx_drdr_times_dx_ds(i) + dx_dr_times_dx_drds(i));

              // d/ds
              sysmat(i,1) = dx_ds(i) - sol(2) * (dx_drds_times_dx_ds(i) + dx_dr_times_dx_dsds(i));

              // d/d(dist)
              sysmat(i,2) = - dx_dr_times_dx_ds(i);


              // residual
              residuum(i) = x(i) - sol(2) * dx_dr_times_dx_ds(i) - x_gp_lin(i);

            }



            sysmat.Invert();

            //solve Newton iteration
            incr.Clear();
            incr.Multiply(-1.0,sysmat,residuum); // incr = -Systemmatrix^-1 * residuum

            // update solution
            sol.Update(1.0, incr, 1.0);

            if ( (incr.Norm2()/sol.Norm2() < absTolIncr) && (residuum.Norm2()/sol.Norm2() < relTolRes) )
            {
              converged = true;
            }

            // check ° relative criterion for local coordinates (between [-1,1]^2)
            //       ° absolute criterion for distance (-> 0)
            //       ° relative criterion for whole residuum
            if(    //sqrt(incr(0)*incr(0)+incr(1)*incr(1))/sqrt(sol(0)*sol(0)+sol(1)*sol(1)) <  relTolIncr
                sqrt(incr(0)*incr(0)+incr(1)*incr(1)) <  absTolIncr
                && incr(2) < absTOLdist
                && residuum.Norm2()/sol.Norm2() < relTolRes)
            {
              converged = true;
            }

          }

          if(!converged)
          {
            cout.precision(15);

            cout << "increment criterion loc coord "
                //<< sqrt(incr(0)*incr(0)+incr(1)*incr(1))/sqrt(sol(0)*sol(0)+sol(1)*sol(1))
                << sqrt(incr(0)*incr(0)+incr(1)*incr(1))
                << " \tabsTOL: " << absTolIncr
                << endl;
            cout << "absolute criterion for distance "
                << incr(2)
                << " \tabsTOL: " << absTOLdist
                << endl;
            cout << "relative criterion whole residuum "
                << residuum.Norm2()/sol.Norm2()
                << " \trelTOL: " << relTolRes
                << endl;


            cout << "sysmat.Invert" << sysmat << endl;
            cout << "sol-norm " << sol.Norm2() << endl;
            cout << "sol " << sol << endl;
            cout << "x_gp_lin" << x_gp_lin << endl;
            cout << "side " << xyze_ << endl;

            dserror( "newton scheme in ProjectOnSide not converged! " );
          }



          // evaluate shape function at solution
          DRT::UTILS::shape_function_2D( side_funct_, sol( 0 ), sol( 1 ), side_distype );

          // get projected gauss point
          x_side.Multiply(xyze_, side_funct_);

          xi_side(0) = sol(0);
          xi_side(1) = sol(1);

        }


        virtual void eivel(const DRT::Discretization &  cutdis,
                           const std::string            state,
                           const vector<int>&           lm)
        {
          // get state of the global vector
          Teuchos::RCP<const Epetra_Vector> matrix_state = cutdis.GetState(state);
          if(matrix_state == null)
            dserror("Cannot get state vector %s", state.c_str());

          // extract local values of the global vectors
          std::vector<double> mymatrix(lm.size());
          DRT::UTILS::ExtractMyValues(*matrix_state,mymatrix,lm);

          for (int inode=0; inode<side_nen_; ++inode)  // number of nodes
          {
            for(int idim=0; idim<3; ++idim) // number of dimensions
            {
              (eivel_)(idim,inode) = mymatrix[idim+(inode*numdof)];  // state vector includes velocity and pressure
            }
            if(numdof == 4) (eipres_)(inode,0) = mymatrix[3+(inode*numdof)];
          }
        }

        virtual void addeidisp(const DRT::Discretization &  cutdis,
                               const std::string            state,
                               const vector<int>&           lm,
                               Epetra_SerialDenseMatrix  &  side_xyze)
        {
          // get state of the global vector
          Teuchos::RCP<const Epetra_Vector> matrix_state = cutdis.GetState(state);
          if(matrix_state == null)
            dserror("Cannot get state vector %s", state.c_str());

          // extract local values of the global vectors
          std::vector<double> mymatrix(lm.size());
          DRT::UTILS::ExtractMyValues(*matrix_state,mymatrix,lm);

          for (int inode=0; inode<side_nen_; ++inode)  // number of nodes
          {
            for(int idim=0; idim<3; ++idim) // number of dimensions
            {
              (eidisp_)(idim,inode) = mymatrix[idim+(inode*numdof)]; // attention! disp state vector has 3+1 dofs for displacement (the same as for (u,p))
            }
          }

          // add the displacement of the interface
          for (int inode = 0; inode < side_nen_; ++inode)
          {
            xyze_(0,inode) += eidisp_(0, inode);
            xyze_(1,inode) += eidisp_(1, inode);
            xyze_(2,inode) += eidisp_(2, inode);
          }

        }

        virtual void buildInterfaceForce( const Teuchos::RCP<Epetra_Vector> &   iforcecol,
                                          const DRT::Discretization &           cutdis,
                                          const vector<int>&                    lm,
                                          LINALG::Matrix<nsd_,1> &              traction,
                                          double &                              fac )
        {

          if(numdof != nsd_) dserror(" pay attention in buildInterfaceForce: numdof != nsd_");

          const Epetra_Map* dofcolmap = cutdis.DofColMap();

          if((int) lm.size() != side_nen_*numdof) dserror("mismatch between number of side nodes and lm.size()");

          for (int inode = 0; inode < side_nen_; ++inode)
          {

            for(int idim=0; idim<3; ++idim )
            {
              int gdof = lm[idim+(inode*numdof)];

              // f^i = ( N^i, t ) = ( N^i, (-pI+2mu*eps(u))*n )
              (*iforcecol)[dofcolmap->LID(gdof)] += side_funct_(inode) * traction(idim) * fac;
            }

          }

        }

        virtual void buildCouplingMatrices(LINALG::Matrix<nsd_,1>        & normal,
                                            const double                fac,
                                            LINALG::Matrix<nen_,1>     & funct,
                                            LINALG::BlockMatrix<LINALG::Matrix<nen_,  1>,6,1> &  rhs
                                            )
        {
          LINALG::Matrix<nen_,side_nen_> bKi_ss;
          LINALG::Matrix<side_nen_,nen_> bKiT_ss;

          // interface velocity vector in gausspoint
          LINALG::Matrix<nsd_,1> ivelint;

          const unsigned Sigmaxx = 0;
          const unsigned Sigmaxy = 1;
          const unsigned Sigmaxz = 2;
          const unsigned Sigmayx = 1;
          const unsigned Sigmayy = 3;
          const unsigned Sigmayz = 4;
          const unsigned Sigmazx = 2;
          const unsigned Sigmazy = 4;
          const unsigned Sigmazz = 5;

          const unsigned Velxi = 0;
          const unsigned Velyi = 1;
          const unsigned Velzi = 2;

          // get velocity at integration point
          // (values at n+alpha_F for generalized-alpha scheme, n+1 otherwise)
          ivelint.Multiply(eivel_,side_funct_);

          for(int i=0;i<nen_;++i)
          {
            for(int j=0;j<side_nen_;++j)
            {
              bKi_ss(i,j) = funct(i)*side_funct_(j);
            }
          }

          // G_sui (coupling only for sigma and u (not for sigma and p)
          BK_sui_( Sigmaxx, Velxi )->Update( fac*normal(0), bKi_ss, 1.0 );
          BK_sui_( Sigmaxy, Velxi )->Update( fac*normal(1), bKi_ss, 1.0 );
          BK_sui_( Sigmaxz, Velxi )->Update( fac*normal(2), bKi_ss, 1.0 );
          BK_sui_( Sigmayx, Velyi )->Update( fac*normal(0), bKi_ss, 1.0 );
          BK_sui_( Sigmayy, Velyi )->Update( fac*normal(1), bKi_ss, 1.0 );
          BK_sui_( Sigmayz, Velyi )->Update( fac*normal(2), bKi_ss, 1.0 );
          BK_sui_( Sigmazx, Velzi )->Update( fac*normal(0), bKi_ss, 1.0 );
          BK_sui_( Sigmazy, Velzi )->Update( fac*normal(1), bKi_ss, 1.0 );
          BK_sui_( Sigmazz, Velzi )->Update( fac*normal(2), bKi_ss, 1.0 );


          rhs( Sigmaxx, 0 )->Update( -fac*normal(0)*ivelint(0), funct, 1.0 );
          rhs( Sigmaxy, 0 )->Update( -fac*normal(1)*ivelint(0), funct, 1.0 );
          rhs( Sigmaxz, 0 )->Update( -fac*normal(2)*ivelint(0), funct, 1.0 );
          rhs( Sigmayx, 0 )->Update( -fac*normal(0)*ivelint(1), funct, 1.0 );
          rhs( Sigmayy, 0 )->Update( -fac*normal(1)*ivelint(1), funct, 1.0 );
          rhs( Sigmayz, 0 )->Update( -fac*normal(2)*ivelint(1), funct, 1.0 );
          rhs( Sigmazx, 0 )->Update( -fac*normal(0)*ivelint(2), funct, 1.0 );
          rhs( Sigmazy, 0 )->Update( -fac*normal(1)*ivelint(2), funct, 1.0 );
          rhs( Sigmazz, 0 )->Update( -fac*normal(2)*ivelint(2), funct, 1.0 );


          bKiT_ss.UpdateT(bKi_ss);

          // G_uis
          BK_uis_( Velxi, Sigmaxx )->Update( fac*normal(0), bKiT_ss, 1.0 );
          BK_uis_( Velxi, Sigmaxy )->Update( fac*normal(1), bKiT_ss, 1.0 );
          BK_uis_( Velxi, Sigmaxz )->Update( fac*normal(2), bKiT_ss, 1.0 );
          BK_uis_( Velyi, Sigmayx )->Update( fac*normal(0), bKiT_ss, 1.0 );
          BK_uis_( Velyi, Sigmayy )->Update( fac*normal(1), bKiT_ss, 1.0 );
          BK_uis_( Velyi, Sigmayz )->Update( fac*normal(2), bKiT_ss, 1.0 );
          BK_uis_( Velzi, Sigmazx )->Update( fac*normal(0), bKiT_ss, 1.0 );
          BK_uis_( Velzi, Sigmazy )->Update( fac*normal(1), bKiT_ss, 1.0 );
          BK_uis_( Velzi, Sigmazz )->Update( fac*normal(2), bKiT_ss, 1.0 );

        }

        virtual void buildFinalCouplingMatrices(LINALG::BlockMatrix<LINALG::Matrix<nen_,nen_>,6,6> &  BinvK_ss,
                                                LINALG::BlockMatrix<LINALG::Matrix<nen_,nen_>,4,6> &  K_iK,
                                                LINALG::BlockMatrix<LINALG::Matrix<nen_,nen_>,6,4> &  K_su,
                                                LINALG::BlockMatrix<LINALG::Matrix<nen_,  1>,6,1> &  rhs
                                               )
        {
          LINALG::BlockMatrix<LINALG::Matrix<side_nen_,nen_>,3,6>      BKi_iK;   //G_uis*Inv(K_ss)
          LINALG::BlockMatrix<LINALG::Matrix<nen_,side_nen_>,4,3>      BCuui;    //K_iK*G_sui
          LINALG::BlockMatrix<LINALG::Matrix<side_nen_,nen_>,3,4>      BCuiu;    //Ki_iK*K_su
          LINALG::BlockMatrix<LINALG::Matrix<side_nen_,side_nen_>,3,3> BCuiui;   //Ki_iK*G_sui
          LINALG::BlockMatrix<LINALG::Matrix<side_nen_, 1>,3,1>       extrhsi;  //G_uis*Inv(K_ss)*rhs


          BCuui  .Multiply( K_iK, BK_sui_ );
          BKi_iK  .Multiply( BK_uis_, BinvK_ss );
          BCuiu  .Multiply( BKi_iK, K_su );
          BCuiui .Multiply( BKi_iK, BK_sui_ );
          extrhsi.Multiply( BKi_iK, rhs );



          // build C_uiu
           for ( unsigned icb=0; icb<4; ++icb ) // u1,u2,u3,p
           {
             for ( unsigned irb=0; irb<3; ++irb ) // no pressure coupling required (one-sided mortaring)
             {
               if ( BCuiu.IsUsed( irb, icb ) )
               {
                 LINALG::Matrix<side_nen_,nen_> & local_BCuiu = *BCuiu( irb, icb );
                 for ( int ic=0; ic<nen_; ++ic ) // ele-nodes
                 {
                   int c = ( 4 )*ic + icb;
                   for ( int ir=0; ir<side_nen_; ++ir )
                   {
                     int r = ( numdof )*ir + irb;
                     C_uiu_( r, c ) -= local_BCuiu( ir, ic );
                   }
                 }
               }
             }
           }

           // irb = 3, pressure component
           if(numdof == 4)
           {
           int irb = 3;
           for ( unsigned icb=0; icb<4; ++icb )
           {
               for ( int ic=0; ic<nen_; ++ic )
               {
                 int c = ( numdof )*ic + icb;
                 for ( int ir=0; ir<side_nen_; ++ir )
                 {
                     int r = ( numdof )*ir + irb;
                     C_uiu_( r , c) = 0.0;
                 }
               }
           }
           }

           // build C_uui
           for ( unsigned icb=0; icb<3; ++icb )  // no pressure coupling required (one-sided mortaring)
           {
             for ( unsigned irb=0; irb<4; ++irb )
             {
               if ( BCuui.IsUsed( irb, icb ) )
               {
                 LINALG::Matrix<nen_,side_nen_> & local_BCuui = *BCuui( irb, icb );
                 for ( int ic=0; ic<side_nen_; ++ic )
                 {
                   int c = ( numdof )*ic + icb;
                   for ( int ir=0; ir<nen_; ++ir )
                   {
                     int r = ( 4 )*ir + irb;
                     C_uui_( r, c ) -= local_BCuui( ir, ic );
                   }
                 }
               }
             }
           }

           if(numdof == 4)
           {
           // icb=3
           int icb=3;
           for ( unsigned irb=0; irb<4; ++irb )
           {
        	   for ( int ic=0; ic<side_nen_; ++ic )
        	   {
        	      int c = ( numdof )*ic + icb;
        	      for ( int ir=0; ir<nen_; ++ir )
        	      {
        	         int r = ( 4 )*ir + irb;
        	         C_uui_( r, c ) = 0.0;
        	      }
        	   }
           }
           }


           // build C_ui (right hand side)
           for ( unsigned irb=0; irb<3; ++irb )  // no pressure coupling required (one-sided mortaring)
           {
             if ( extrhsi.IsUsed( irb, 0 ) )
             {
               LINALG::Matrix<side_nen_,1> & local_extrhsi = *extrhsi( irb, 0 );
               for ( int ir=0; ir<side_nen_; ++ir )
               {
                 unsigned r = numdof*ir + irb;
                 rhC_ui_( r, 0 ) -= local_extrhsi( ir, 0 );
               }
             }
           }

           if(numdof == 4)
           {
           int irb = 3;
           for ( int ir=0; ir<side_nen_; ++ir )
           {
               unsigned r = numdof*ir + irb;
               rhC_ui_( r,0) = 0.0;
           }
           }

           // build K_sui_
           // icb und irb: Richtungen
            for ( unsigned icb=0; icb<3; ++icb )
            {
              for ( unsigned irb=0; irb<6; ++irb )
              {
                if ( BK_sui_.IsUsed( irb, icb ) )
                {
                  LINALG::Matrix<nen_,side_nen_> & local_BK_sui = *BK_sui_( irb, icb );
                  for ( int ic=0; ic<side_nen_; ++ic )
                  {
                    int c = ( numdof )*ic + icb;
                    for ( int ir=0; ir<nen_; ++ir )
                    {
                      int r = ( 6 )*ir + irb;
                      K_sui_( r, c ) = local_BK_sui( ir, ic );
                    }
                  }
                }
              }
            }

            // build K_uis_
            for ( unsigned icb=0; icb<6; ++icb )
            {
              for ( unsigned irb=0; irb<3; ++irb )
              {
                if ( BK_uis_.IsUsed( irb, icb ) )
                {
                  LINALG::Matrix<side_nen_,nen_> & local_BK_uis = *BK_uis_( irb, icb );
                  for ( int ic=0; ic<nen_; ++ic )
                  {
                    int c = ( 6 )*ic + icb;
                    for ( int ir=0; ir<side_nen_; ++ir )
                    {
                      int r = ( numdof )*ir + irb;
                      K_uis_( r, c ) = local_BK_uis( ir, ic );
                    }
                  }
                }
              }
           }



        }



        virtual void buildCouplingMatricesNitsche(  Epetra_SerialDenseMatrix &    C_uu_,          // standard bg-bg-matrix
                                                    Epetra_SerialDenseVector &    rhs_Cu_,        // standard bg-rhs
                                                    bool &                        coupling,       // assemble coupling terms (yes/no)
                                                    bool &                        bg_mortaring,   // yes: background-sided mortaring, no: coupling between two meshes (mixed mortaring)
                                                    LINALG::Matrix<nsd_,1> &      normal,         // normal vector
                                                    const double                  timefacfac,     // theta*dt
                                                    const double                  visceff_1,      // viscosity in background fluid
                                                    const double                  visceff_2,      // viscosity in embedded fluid
                                                    double &                      kappa1,         // mortaring weighting
                                                    double &                      kappa2,         // mortaring weighting
                                                    double &                      stabfac,        // Nitsche non-dimensionless stabilization factor
                                                    double &                      stabfac_conv,   // Nitsche convective non-dimensionless stabilization factor
                                                    LINALG::Matrix<nen_,1> &      funct_,         // bg shape functions
                                                    LINALG::Matrix<nsd_,nen_> &   derxy_,         // bg deriv
                                                    LINALG::Matrix<nsd_,nsd_> &   vderxy_,        // bg deriv^n
                                                    double &                      press,          // bg p^n
                                                    LINALG::Matrix<nsd_,1> &      velint,         // bg u^n
                                                    LINALG::Matrix<nsd_,1> &      ivelint_WDBC_JUMP // Dirichlet velocity vector or prescribed jump vector
                                               )
        {

          const unsigned Velx = 0;
          const unsigned Vely = 1;
          const unsigned Velz = 2;

          //--------------------------------------------

          // define the coupling between two not matching grids
          // for fluidfluidcoupling
          // domain Omega^1 := Xfluid
          // domain Omega^2 := Alefluid( or monolithic: structure) ( not available for non-coupling (Dirichlet) )

          // [| v |] := v1 - v2
          //  { v }  := kappa1 * v1 + kappa2 * v2 = kappa1 * v1 (for Dirichlet coupling k1=1.0, k2 = 0.0)
          //  < v >  := kappa2 * v1 + kappa1 * v2 = kappa1 * v2 (for Dirichlet coupling k1=1.0, k2 = 0.0)
//
          double k1mu1_fac = 2.0 * timefacfac * kappa1 * visceff_1;
//          //const double k2mu2_fac = 2.0 * timefacfac * kappa2 * visceff_2;


          //--------------------------------------------
          // get fields at interface


          // get velocity at integration point
          // (values at n+alpha_F for generalized-alpha scheme, n+1 otherwise)
          // interface velocity vector in gausspoint
          LINALG::Matrix<nsd_,1> ivelint;
          ivelint.Multiply(eivel_,side_funct_);


          // funct_ * timefac * fac * funct_ (dyadic product)
          LINALG::Matrix<nen_,1> funct_timefacfac(true);
          funct_timefacfac.Update(timefacfac,funct_,0.0);

          LINALG::Matrix<side_nen_,1> side_funct_timefacfac(true);
          side_funct_timefacfac.Update(timefacfac,side_funct_,0.0);


                     /*                  \       /          i      \
                  + |  [ v ],   {Dp}*n    | = - | [ v ], { p }* n   |
                     \                   /       \                */

          //-----------------------------------------------
          //    + (v1, k1 *(Dp1)*n)
          //-----------------------------------------------
          LINALG::Matrix<nen_,nen_> funct_dyad_timefacfac(true);
          LINALG::Matrix<nen_,nen_> funct_dyad_k1_timefacfac(true);
          funct_dyad_timefacfac.MultiplyNT(funct_timefacfac, funct_);
          funct_dyad_k1_timefacfac.Update(kappa1,funct_dyad_timefacfac,0.0);

          for(int ir = 0; ir<nen_; ir++)
          {
            int idVelx = ir*(nsd_+1) + 0;
            int idVely = ir*(nsd_+1) + 1;
            int idVelz = ir*(nsd_+1) + 2;

            // (v,Dp*n)
            for(int ic =0; ic<nen_; ic++)
            {
              int iPres = ic*(nsd_+1)+3;

              C_uu_(idVelx, iPres) += funct_dyad_k1_timefacfac(ir,ic)*normal(Velx);
              C_uu_(idVely, iPres) += funct_dyad_k1_timefacfac(ir,ic)*normal(Vely);
              C_uu_(idVelz, iPres) += funct_dyad_k1_timefacfac(ir,ic)*normal(Velz);
            }

            // -(v,p*n)
            double funct_k1_timefacfac_press = funct_timefacfac(ir)*press*kappa1;
            rhs_Cu_(idVelx,0) -= funct_k1_timefacfac_press*normal(Velx);
            rhs_Cu_(idVely,0) -= funct_k1_timefacfac_press*normal(Vely);
            rhs_Cu_(idVelz,0) -= funct_k1_timefacfac_press*normal(Velz);
          }


          LINALG::Matrix<side_nen_,nen_> side_funct_dyad_timefacfac(true);
          LINALG::Matrix<side_nen_,nen_> side_funct_dyad_k1_timefacfac(true);
          side_funct_dyad_timefacfac.MultiplyNT(side_funct_timefacfac, funct_);
          side_funct_dyad_k1_timefacfac.Update(kappa1, side_funct_dyad_timefacfac, 0.0);

          if(coupling)
          {
          //-----------------------------------------------
          //    - (v2, k1 *(Dp1)*n)
          //-----------------------------------------------
          for(int ir = 0; ir<side_nen_; ir++)
          {
            int idVelx = ir*(nsd_+1) + 0;
            int idVely = ir*(nsd_+1) + 1;
            int idVelz = ir*(nsd_+1) + 2;

            // (v,Dp*n)
            for(int ic =0; ic<nen_; ic++)
            {
              int iPres = ic*(nsd_+1)+3;

              C_uiu_(idVelx, iPres) -= side_funct_dyad_k1_timefacfac(ir,ic)*normal(Velx);
              C_uiu_(idVely, iPres) -= side_funct_dyad_k1_timefacfac(ir,ic)*normal(Vely);
              C_uiu_(idVelz, iPres) -= side_funct_dyad_k1_timefacfac(ir,ic)*normal(Velz);
            }

            // -(v,p*n)
            double side_funct_k1_timefacfac_press = side_funct_timefacfac(ir)*press*kappa1;
            rhC_ui_(idVelx,0) += side_funct_k1_timefacfac_press*normal(Velx);
            rhC_ui_(idVely,0) += side_funct_k1_timefacfac_press*normal(Vely);
            rhC_ui_(idVelz,0) += side_funct_k1_timefacfac_press*normal(Velz);
          }
          }// end coupling



              /*                   \     /              i   \
           - |  { q }*n ,[ Du ]     | = |  { q }*n  ,[ u ]   |
              \                    /     \                 */

          // -1.0 antisymmetric
          // +1.0 symmetric
//          double alpha_p = -1.0;
          double alpha_p = -1.0;

         //-----------------------------------------------
         //    - (q1*n, k1 *(Du1))
         //-----------------------------------------------
         for(int ir = 0; ir<nen_; ir++)
         {
          int idPres = ir*(nsd_+1) + 3;

          // -(q*n,Du)
          for(int ic =0; ic<nen_; ic++)
          {
            int iVelx = ic*(nsd_+1)+0;
            int iVely = ic*(nsd_+1)+1;
            int iVelz = ic*(nsd_+1)+2;

            C_uu_(idPres, iVelx) += alpha_p*funct_dyad_k1_timefacfac(ir,ic)*normal(Velx);
            C_uu_(idPres, iVely) += alpha_p*funct_dyad_k1_timefacfac(ir,ic)*normal(Vely);
            C_uu_(idPres, iVelz) += alpha_p*funct_dyad_k1_timefacfac(ir,ic)*normal(Velz);
          }

          // (q*n,u)
          double velint_normal = velint.Dot(normal);
          rhs_Cu_(idPres,0) -= alpha_p*funct_timefacfac(ir)*kappa1*velint_normal;

          if(!coupling) // weak Dirichlet case
          {
            // -(q*n,u_DBC)
            double ivelint_WDBC_JUMP_normal = ivelint_WDBC_JUMP.Dot(normal);
            rhs_Cu_(idPres,0) += alpha_p*funct_timefacfac(ir)*kappa1*ivelint_WDBC_JUMP_normal;
          }
         }


         if(coupling)
         {
         //-----------------------------------------------
         //    + (q1*n, k1 *(Du2))
         //-----------------------------------------------
         for(int ir = 0; ir<nen_; ir++)
         {
          int idPres = ir*(nsd_+1) + 3;

          // -(q*n,Du)
          for(int ic =0; ic<side_nen_; ic++)
          {
            int iVelx = ic*(nsd_+1)+0;
            int iVely = ic*(nsd_+1)+1;
            int iVelz = ic*(nsd_+1)+2;

            C_uui_(idPres, iVelx) -= alpha_p*side_funct_dyad_k1_timefacfac(ic,ir)*normal(Velx);
            C_uui_(idPres, iVely) -= alpha_p*side_funct_dyad_k1_timefacfac(ic,ir)*normal(Vely);
            C_uui_(idPres, iVelz) -= alpha_p*side_funct_dyad_k1_timefacfac(ic,ir)*normal(Velz);
          }

          // -(q*n,u)
          double ivelint_normal = ivelint.Dot(normal);
          rhs_Cu_(idPres,0) += alpha_p*funct_timefacfac(ir)*kappa1*ivelint_normal;

         }
         }// end coupling






          /*                           \       /                   i      \
       - |  [ v ],  { 2mu eps(u) }*n    | = + | [ v ],  { 2mu eps(u ) }*n  |
          \                            /       \                         */

          //-----------------------------------------------
          //    - (v1, (2*k1*mu1) *eps(Du1)*n)
          //-----------------------------------------------


          LINALG::Matrix<nen_,1> e_funct_visc1_timefacfac(true);
          e_funct_visc1_timefacfac.Update(k1mu1_fac, funct_, 0.0);

          //            LINALG::Matrix<side_nen_,1> s_funct_visc_timefacfac(true);
          //            s_funct_visc_timefacfac.Update(k2mu2_fac, side_funct_, 0.0);

          for(int ir = 0; ir<nen_; ir++)
          {
            int idVelx = ir*(nsd_+1) + 0;
            int idVely = ir*(nsd_+1) + 1;
            int idVelz = ir*(nsd_+1) + 2;
//            int idPres = ir*(nsd_+1) + 3;


            for(int ic =0; ic<nen_; ic++)
            {
              int iVelx = ic*(nsd_+1)+0;
              int iVely = ic*(nsd_+1)+1;
              int iVelz = ic*(nsd_+1)+2;
//              int iPres = ic*(nsd_+1)+3;


              // - (v1, (2*k1*mu1) *eps(Du1)*n)

              //(x,x)
              C_uu_(idVelx, iVelx) -= e_funct_visc1_timefacfac(ir)*(         normal(Velx)*derxy_(Velx,ic)
                                                                     + 0.5 * normal(Vely)*derxy_(Vely,ic)
                                                                     + 0.5 * normal(Velz)*derxy_(Velz,ic)  );
              //(x,y)
              C_uu_(idVelx, iVely) -= e_funct_visc1_timefacfac(ir)*    0.5 * normal(Vely)*derxy_(Velx,ic);
              //(x,z)
              C_uu_(idVelx, iVelz) -= e_funct_visc1_timefacfac(ir)*    0.5 * normal(Velz)*derxy_(Velx,ic);

              //(y,x)
              C_uu_(idVely, iVelx) -= e_funct_visc1_timefacfac(ir)*    0.5 * normal(Velx)*derxy_(Vely,ic);
              //(y,y)
              C_uu_(idVely, iVely) -= e_funct_visc1_timefacfac(ir)*(   0.5 * normal(Velx)*derxy_(Velx,ic)
                                                                     +       normal(Vely)*derxy_(Vely,ic)
                                                                     + 0.5 * normal(Velz)*derxy_(Velz,ic)  );
              //(y,z)
              C_uu_(idVely, iVelz) -= e_funct_visc1_timefacfac(ir)*    0.5 * normal(Velz)*derxy_(Vely,ic);

              //(z,x)
              C_uu_(idVelz, iVelx) -= e_funct_visc1_timefacfac(ir)*    0.5 * normal(Velx)*derxy_(Velz,ic);
              //(z,y)
              C_uu_(idVelz, iVely) -= e_funct_visc1_timefacfac(ir)*    0.5 * normal(Vely)*derxy_(Velz,ic);
              //(z,z)
              C_uu_(idVelz, iVelz) -= e_funct_visc1_timefacfac(ir)*(   0.5 * normal(Velx)*derxy_(Velx,ic)
                                                                     + 0.5 * normal(Vely)*derxy_(Vely,ic)
                                                                     +       normal(Velz)*derxy_(Velz,ic)  );
            }

            // - (v1, (2*k1*mu1) *eps(Du1)*n)
            rhs_Cu_(idVelx,0) += e_funct_visc1_timefacfac(ir)*(            vderxy_(Velx,Velx)                      *normal(Velx)
                                                                 + 0.5 * ( vderxy_(Velx,Vely) + vderxy_(Vely,Velx))*normal(Vely)
                                                                 + 0.5 * ( vderxy_(Velx,Velz) + vderxy_(Velz,Velx))*normal(Velz)  );
            rhs_Cu_(idVely,0) += e_funct_visc1_timefacfac(ir)*(    0.5 * ( vderxy_(Vely,Velx) + vderxy_(Velx,Vely))*normal(Velx)
                                                                 +         vderxy_(Vely,Vely)                      *normal(Vely)
                                                                 + 0.5 * ( vderxy_(Vely,Velz) + vderxy_(Velz,Vely))*normal(Velz)  );
            rhs_Cu_(idVelz,0) += e_funct_visc1_timefacfac(ir)*(    0.5 * ( vderxy_(Velz,Velx) + vderxy_(Velx,Velz))*normal(Velx)
                                                                 + 0.5 * ( vderxy_(Velz,Vely) + vderxy_(Vely,Velz))*normal(Vely)
                                                                 +         vderxy_(Velz,Velz)                      *normal(Velz)  );

          }


          LINALG::Matrix<side_nen_,1> s_funct_visc1_timefacfac(true);
          s_funct_visc1_timefacfac.Update(k1mu1_fac, side_funct_, 0.0);
          if(coupling)
          {
          //-----------------------------------------------
          //    + (v2, (2*k1*mu1) *eps(Du1)*n)
          //-----------------------------------------------
          for(int ir = 0; ir<side_nen_; ir++)
          {
            int idVelx = ir*(nsd_+1) + 0;
            int idVely = ir*(nsd_+1) + 1;
            int idVelz = ir*(nsd_+1) + 2;
            //int idPres = ir*(nsd_+1) + 3;


            for(int ic =0; ic<nen_; ic++)
            {
              int iVelx = ic*(nsd_+1)+0;
              int iVely = ic*(nsd_+1)+1;
              int iVelz = ic*(nsd_+1)+2;
              //int iPres = ic*(nsd_+1)+3;


              // + (v2, (2*k1*mu1) *eps(Du1)*n)

              //(x,x)
              C_uiu_(idVelx, iVelx) += s_funct_visc1_timefacfac(ir)*(         normal(Velx)*derxy_(Velx,ic)
                                                                      + 0.5 * normal(Vely)*derxy_(Vely,ic)
                                                                      + 0.5 * normal(Velz)*derxy_(Velz,ic)  );
              //(x,y)
              C_uiu_(idVelx, iVely) += s_funct_visc1_timefacfac(ir)*    0.5 * normal(Vely)*derxy_(Velx,ic);
              //(x,z)
              C_uiu_(idVelx, iVelz) += s_funct_visc1_timefacfac(ir)*    0.5 * normal(Velz)*derxy_(Velx,ic);

              //(y,x)
              C_uiu_(idVely, iVelx) += s_funct_visc1_timefacfac(ir)*    0.5 * normal(Velx)*derxy_(Vely,ic);
              //(y,y)
              C_uiu_(idVely, iVely) += s_funct_visc1_timefacfac(ir)*(   0.5 * normal(Velx)*derxy_(Velx,ic)
                                                                      +       normal(Vely)*derxy_(Vely,ic)
                                                                      + 0.5 * normal(Velz)*derxy_(Velz,ic)  );
              //(y,z)
              C_uiu_(idVely, iVelz) += s_funct_visc1_timefacfac(ir)*    0.5 * normal(Velz)*derxy_(Vely,ic);

              //(z,x)
              C_uiu_(idVelz, iVelx) += s_funct_visc1_timefacfac(ir)*    0.5 * normal(Velx)*derxy_(Velz,ic);
              //(z,y)
              C_uiu_(idVelz, iVely) += s_funct_visc1_timefacfac(ir)*    0.5 * normal(Vely)*derxy_(Velz,ic);
              //(z,z)
              C_uiu_(idVelz, iVelz) += s_funct_visc1_timefacfac(ir)*(   0.5 * normal(Velx)*derxy_(Velx,ic)
                                                                      + 0.5 * normal(Vely)*derxy_(Vely,ic)
                                                                      +       normal(Velz)*derxy_(Velz,ic)  );
            }

            // - (v2, (2*k1*mu1) *eps(Du1)*n)
            rhC_ui_(idVelx,0) -= s_funct_visc1_timefacfac(ir)*(            vderxy_(Velx,Velx)                      *normal(Velx)
                                                                 + 0.5 * ( vderxy_(Velx,Vely) + vderxy_(Vely,Velx))*normal(Vely)
                                                                 + 0.5 * ( vderxy_(Velx,Velz) + vderxy_(Velz,Velx))*normal(Velz)  );
            rhC_ui_(idVely,0) -= s_funct_visc1_timefacfac(ir)*(    0.5 * ( vderxy_(Vely,Velx) + vderxy_(Velx,Vely))*normal(Velx)
                                                                 +         vderxy_(Vely,Vely)                      *normal(Vely)
                                                                 + 0.5 * ( vderxy_(Vely,Velz) + vderxy_(Velz,Vely))*normal(Velz)  );
            rhC_ui_(idVelz,0) -= s_funct_visc1_timefacfac(ir)*(    0.5 * ( vderxy_(Velz,Velx) + vderxy_(Velx,Velz))*normal(Velx)
                                                                 + 0.5 * ( vderxy_(Velz,Vely) + vderxy_(Vely,Velz))*normal(Vely)
                                                                 +         vderxy_(Velz,Velz)                      *normal(Velz)  );

          }
          }// end coupling




            /*                                \       /                             i   \
         - |  alpha* { 2mu*eps(v) }*n , [ Du ] |  =  |  alpha* { 2mu eps(v) }*n ,[ u ]   |
            \                                 /       \                                */
         // antisymmetric formulation (see Burman, Fernandez 2009)

          // +1.0 symmetric
          // -1.0 antisymmetric
//          double alpha = +1.0;
          double alpha = +1.0;


         //-----------------------------------------------
         //    - ((2*k1*mu1) *eps(v1)*n , u1)
         //-----------------------------------------------
         for(int ir = 0; ir<nen_; ir++)
         {
         int idVelx = ir*(nsd_+1) + 0;
         int idVely = ir*(nsd_+1) + 1;
         int idVelz = ir*(nsd_+1) + 2;

         // -(2mu*eps(v)*n, Du)
         for(int ic =0; ic<nen_; ic++)
         {
           int iVelx = ic*(nsd_+1)+0;
           int iVely = ic*(nsd_+1)+1;
           int iVelz = ic*(nsd_+1)+2;

           //(x,x)
           C_uu_(idVelx, iVelx) -= alpha*e_funct_visc1_timefacfac(ic)*(         normal(Velx)*derxy_(Velx,ir)
                                                                        + 0.5 * normal(Vely)*derxy_(Vely,ir)
                                                                        + 0.5 * normal(Velz)*derxy_(Velz,ir)  );
           //(y,x)
           C_uu_(idVely, iVelx) -= alpha*e_funct_visc1_timefacfac(ic)*    0.5 * normal(Vely)*derxy_(Velx,ir);
           //(z,x)
           C_uu_(idVelz, iVelx) -= alpha*e_funct_visc1_timefacfac(ic)*    0.5 * normal(Velz)*derxy_(Velx,ir);

           //(x,y)
           C_uu_(idVelx, iVely) -= alpha*e_funct_visc1_timefacfac(ic)*    0.5 * normal(Velx)*derxy_(Vely,ir);
           //(y,y)
           C_uu_(idVely, iVely) -= alpha*e_funct_visc1_timefacfac(ic)*(   0.5 * normal(Velx)*derxy_(Velx,ir)
                                                                        +       normal(Vely)*derxy_(Vely,ir)
                                                                        + 0.5 * normal(Velz)*derxy_(Velz,ir)  );
           //(z,y)
           C_uu_(idVelz, iVely) -= alpha*e_funct_visc1_timefacfac(ic)*    0.5 * normal(Velz)*derxy_(Vely,ir);

           //(x,z)
           C_uu_(idVelx, iVelz) -= alpha*e_funct_visc1_timefacfac(ic)*    0.5 * normal(Velx)*derxy_(Velz,ir);
           //(y,z)
           C_uu_(idVely, iVelz) -= alpha*e_funct_visc1_timefacfac(ic)*    0.5 * normal(Vely)*derxy_(Velz,ir);
           //(z,z)
           C_uu_(idVelz, iVelz) -= alpha*e_funct_visc1_timefacfac(ic)*(   0.5 * normal(Velx)*derxy_(Velx,ir)
                                                                        + 0.5 * normal(Vely)*derxy_(Vely,ir)
                                                                        +       normal(Velz)*derxy_(Velz,ir)  );
         }
         //  (2mu1*eps(v1)*n, u1)
         double timefacfac_visc = alpha*timefacfac*2.0*visceff_1;
         rhs_Cu_(idVelx,0) += timefacfac_visc* (     derxy_(Velx,ir) *       normal(Velx) * velint(Velx)
                                                   + derxy_(Vely,ir) * 0.5* (normal(Vely) * velint(Velx) + normal(Velx)*velint(Vely))
                                                   + derxy_(Velz,ir) * 0.5* (normal(Velz) * velint(Velx) + normal(Velx)*velint(Velz)));

         rhs_Cu_(idVely,0) += timefacfac_visc* (     derxy_(Velx,ir) * 0.5* (normal(Vely) * velint(Velx) + normal(Velx)*velint(Vely))
                                                   + derxy_(Vely,ir) *       normal(Vely) * velint(Vely)
                                                   + derxy_(Velz,ir) * 0.5* (normal(Velz) * velint(Vely) + normal(Vely)*velint(Velz)));

         rhs_Cu_(idVelz,0) += timefacfac_visc* (     derxy_(Velx,ir) * 0.5* (normal(Velx) * velint(Velz) + normal(Velz)*velint(Velx))
                                                   + derxy_(Vely,ir) * 0.5* (normal(Vely) * velint(Velz) + normal(Velz)*velint(Vely))
                                                   + derxy_(Velz,ir) *       normal(Velz) * velint(Velz));

         if(!coupling) // weak Dirichlet case
         {
           // -(2mu*eps(v)*n, u_DBC)
           rhs_Cu_(idVelx,0) -= timefacfac_visc* (  derxy_(Velx,ir) *       normal(Velx) * ivelint_WDBC_JUMP(Velx)
                                                  + derxy_(Vely,ir) * 0.5* (normal(Vely) * ivelint_WDBC_JUMP(Velx) + normal(Velx)*ivelint_WDBC_JUMP(Vely))
                                                  + derxy_(Velz,ir) * 0.5* (normal(Velz) * ivelint_WDBC_JUMP(Velx) + normal(Velx)*ivelint_WDBC_JUMP(Velz)));

           rhs_Cu_(idVely,0) -= timefacfac_visc* (  derxy_(Velx,ir) * 0.5* (normal(Vely) * ivelint_WDBC_JUMP(Velx) + normal(Velx)*ivelint_WDBC_JUMP(Vely))
                                                  + derxy_(Vely,ir) *       normal(Vely) * ivelint_WDBC_JUMP(Vely)
                                                  + derxy_(Velz,ir) * 0.5* (normal(Velz) * ivelint_WDBC_JUMP(Vely) + normal(Vely)*ivelint_WDBC_JUMP(Velz)));

           rhs_Cu_(idVelz,0) -= timefacfac_visc* (  derxy_(Velx,ir) * 0.5* (normal(Velx) * ivelint_WDBC_JUMP(Velz) + normal(Velz)*ivelint_WDBC_JUMP(Velx))
                                                  + derxy_(Vely,ir) * 0.5* (normal(Vely) * ivelint_WDBC_JUMP(Velz) + normal(Velz)*ivelint_WDBC_JUMP(Vely))
                                                  + derxy_(Velz,ir) *       normal(Velz) * ivelint_WDBC_JUMP(Velz));

         }
         }

         if(coupling)
         {
         //-----------------------------------------------
         //    + ((2*k1*mu1) *eps(v1)*n , u2)
         //-----------------------------------------------
         for(int ir = 0; ir<nen_; ir++)
         {
         int idVelx = ir*(nsd_+1) + 0;
         int idVely = ir*(nsd_+1) + 1;
         int idVelz = ir*(nsd_+1) + 2;

         // -(2mu*eps(v1)*n, Du2)
         for(int ic =0; ic<side_nen_; ic++)
         {
           int iVelx = ic*(nsd_+1)+0;
           int iVely = ic*(nsd_+1)+1;
           int iVelz = ic*(nsd_+1)+2;

           //(x,x)
           C_uui_(idVelx, iVelx) += alpha*s_funct_visc1_timefacfac(ic)*(         normal(Velx)*derxy_(Velx,ir)
                                                                         + 0.5 * normal(Vely)*derxy_(Vely,ir)
                                                                         + 0.5 * normal(Velz)*derxy_(Velz,ir)  );
           //(y,x)
           C_uui_(idVely, iVelx) += alpha*s_funct_visc1_timefacfac(ic)*    0.5 * normal(Vely)*derxy_(Velx,ir);
           //(z,x)
           C_uui_(idVelz, iVelx) += alpha*s_funct_visc1_timefacfac(ic)*    0.5 * normal(Velz)*derxy_(Velx,ir);

           //(x,y)
           C_uui_(idVelx, iVely) += alpha*s_funct_visc1_timefacfac(ic)*    0.5 * normal(Velx)*derxy_(Vely,ir);
           //(y,y)
           C_uui_(idVely, iVely) += alpha*s_funct_visc1_timefacfac(ic)*(   0.5 * normal(Velx)*derxy_(Velx,ir)
                                                                         +       normal(Vely)*derxy_(Vely,ir)
                                                                         + 0.5 * normal(Velz)*derxy_(Velz,ir)  );
           //(z,y)
           C_uui_(idVelz, iVely) += alpha*s_funct_visc1_timefacfac(ic)*    0.5 * normal(Velz)*derxy_(Vely,ir);

           //(x,z)
           C_uui_(idVelx, iVelz) += alpha*s_funct_visc1_timefacfac(ic)*    0.5 * normal(Velx)*derxy_(Velz,ir);
           //(y,z)
           C_uui_(idVely, iVelz) += alpha*s_funct_visc1_timefacfac(ic)*    0.5 * normal(Vely)*derxy_(Velz,ir);
           //(z,z)
           C_uui_(idVelz, iVelz) += alpha*s_funct_visc1_timefacfac(ic)*(   0.5 * normal(Velx)*derxy_(Velx,ir)
                                                                         + 0.5 * normal(Vely)*derxy_(Vely,ir)
                                                                         +       normal(Velz)*derxy_(Velz,ir)  );
         }
         //  (2mu1*eps(v1)*n, u2)
         double timefacfac_visc = alpha*timefacfac*2.0*visceff_1;
         rhs_Cu_(idVelx,0) -= timefacfac_visc* (     derxy_(Velx,ir) *       normal(Velx) * ivelint(Velx)
                                                   + derxy_(Vely,ir) * 0.5* (normal(Vely) * ivelint(Velx) + normal(Velx)*ivelint(Vely))
                                                   + derxy_(Velz,ir) * 0.5* (normal(Velz) * ivelint(Velx) + normal(Velx)*ivelint(Velz)));

         rhs_Cu_(idVely,0) -= timefacfac_visc* (     derxy_(Velx,ir) * 0.5* (normal(Vely) * ivelint(Velx) + normal(Velx)*ivelint(Vely))
                                                   + derxy_(Vely,ir) *       normal(Vely) * ivelint(Vely)
                                                   + derxy_(Velz,ir) * 0.5* (normal(Velz) * ivelint(Vely) + normal(Vely)*ivelint(Velz)));

         rhs_Cu_(idVelz,0) -= timefacfac_visc* (     derxy_(Velx,ir) * 0.5* (normal(Velx) * ivelint(Velz) + normal(Velz)*ivelint(Velx))
                                                   + derxy_(Vely,ir) * 0.5* (normal(Vely) * ivelint(Velz) + normal(Velz)*ivelint(Vely))
                                                   + derxy_(Velz,ir) *       normal(Velz) * ivelint(Velz));

//         // -(2mu*eps(v)*n, u_DBC)
//         elevec1_epetra(idVelx,0) -= timefacfac_visc* (  derxy_(Velx,ir) *       normal(Velx) * velint_WDBC(Velx)
//                                                   + derxy_(Vely,ir) * 0.5* (normal(Vely) * velint_WDBC(Velx) + normal(Velx)*velint_WDBC(Vely))
//                                                   + derxy_(Velz,ir) * 0.5* (normal(Velz) * velint_WDBC(Velx) + normal(Velx)*velint_WDBC(Velz)));
//
//         elevec1_epetra(idVely,0) -= timefacfac_visc* (  derxy_(Velx,ir) * 0.5* (normal(Vely) * velint_WDBC(Velx) + normal(Velx)*velint_WDBC(Vely))
//                                                   + derxy_(Vely,ir) *       normal(Vely) * velint_WDBC(Vely)
//                                                   + derxy_(Velz,ir) * 0.5* (normal(Velz) * velint_WDBC(Vely) + normal(Vely)*velint_WDBC(Velz)));
//
//         elevec1_epetra(idVelz,0) -= timefacfac_visc* (  derxy_(Velx,ir) * 0.5* (normal(Velx) * velint_WDBC(Velz) + normal(Velz)*velint_WDBC(Velx))
//                                                   + derxy_(Vely,ir) * 0.5* (normal(Vely) * velint_WDBC(Velz) + normal(Velz)*velint_WDBC(Vely))
//                                                   + derxy_(Velz,ir) *       normal(Velz) * velint_WDBC(Velz));
         }
         }// end coupling




            /*                                  \        /                           i   \
           |  gamma*mu/h_K *  [ v ] , [ Du ]     | =  - |   gamma*mu/h_K * [ v ], [ u ]   |
            \                                   /        \                              */


         // + gamma*mu/h_K (v1, u1))
         for(int ir=0; ir<nen_; ir++)
         {
           int idVelx = ir*(nsd_+1) + 0;
           int idVely = ir*(nsd_+1) + 1;
           int idVelz = ir*(nsd_+1) + 2;

           for(int ic=0; ic<nen_; ic++)
           {
             int iVelx = ic*(nsd_+1)+0;
             int iVely = ic*(nsd_+1)+1;
             int iVelz = ic*(nsd_+1)+2;

             C_uu_(idVelx, iVelx) += funct_dyad_timefacfac(ir,ic)*stabfac;
             C_uu_(idVely, iVely) += funct_dyad_timefacfac(ir,ic)*stabfac;
             C_uu_(idVelz, iVelz) += funct_dyad_timefacfac(ir,ic)*stabfac;
           }

           // -(stab * v, u)
           rhs_Cu_(idVelx,0) -= funct_timefacfac(ir)*stabfac*velint(Velx);
           rhs_Cu_(idVely,0) -= funct_timefacfac(ir)*stabfac*velint(Vely);
           rhs_Cu_(idVelz,0) -= funct_timefacfac(ir)*stabfac*velint(Velz);

           if(!coupling) // weak Dirichlet case
           {
             // +(stab * v, u_DBC)
             rhs_Cu_(idVelx,0) += funct_timefacfac(ir)*stabfac*ivelint_WDBC_JUMP(Velx);
             rhs_Cu_(idVely,0) += funct_timefacfac(ir)*stabfac*ivelint_WDBC_JUMP(Vely);
             rhs_Cu_(idVelz,0) += funct_timefacfac(ir)*stabfac*ivelint_WDBC_JUMP(Velz);
           }

         }

         if(coupling)
         {
         // - gamma*mu/h_K (v1, u2))
         for(int ir=0; ir<nen_; ir++)
         {
           int idVelx = ir*(nsd_+1) + 0;
           int idVely = ir*(nsd_+1) + 1;
           int idVelz = ir*(nsd_+1) + 2;

           for(int ic=0; ic<side_nen_; ic++)
           {
             int iVelx = ic*(nsd_+1)+0;
             int iVely = ic*(nsd_+1)+1;
             int iVelz = ic*(nsd_+1)+2;

             C_uui_(idVelx, iVelx) -= side_funct_dyad_timefacfac(ic,ir)*stabfac;
             C_uui_(idVely, iVely) -= side_funct_dyad_timefacfac(ic,ir)*stabfac;
             C_uui_(idVelz, iVelz) -= side_funct_dyad_timefacfac(ic,ir)*stabfac;
           }

           // -(stab * v, u)
           rhs_Cu_(idVelx,0) += funct_timefacfac(ir)*stabfac*ivelint(Velx);
           rhs_Cu_(idVely,0) += funct_timefacfac(ir)*stabfac*ivelint(Vely);
           rhs_Cu_(idVelz,0) += funct_timefacfac(ir)*stabfac*ivelint(Velz);


         }


         // - gamma*mu/h_K (v2, u1))
         for(int ir=0; ir<side_nen_; ir++)
         {
           int idVelx = ir*(nsd_+1) + 0;
           int idVely = ir*(nsd_+1) + 1;
           int idVelz = ir*(nsd_+1) + 2;

           for(int ic=0; ic<nen_; ic++)
           {
             int iVelx = ic*(nsd_+1)+0;
             int iVely = ic*(nsd_+1)+1;
             int iVelz = ic*(nsd_+1)+2;

             C_uiu_(idVelx, iVelx) -= side_funct_dyad_timefacfac(ir,ic)*stabfac;
             C_uiu_(idVely, iVely) -= side_funct_dyad_timefacfac(ir,ic)*stabfac;
             C_uiu_(idVelz, iVelz) -= side_funct_dyad_timefacfac(ir,ic)*stabfac;
           }

           // +(stab * v2, u1)
           rhC_ui_(idVelx,0) += side_funct_timefacfac(ir)*stabfac*velint(Velx);
           rhC_ui_(idVely,0) += side_funct_timefacfac(ir)*stabfac*velint(Vely);
           rhC_ui_(idVelz,0) += side_funct_timefacfac(ir)*stabfac*velint(Velz);

         }

         LINALG::Matrix<side_nen_,side_nen_> side_side_dyad_timefacfac(true);
         side_side_dyad_timefacfac.MultiplyNT(side_funct_timefacfac, side_funct_);

         // + gamma*mu/h_K (v2, u2))
         for(int ir=0; ir<side_nen_; ir++)
         {
           int idVelx = ir*(nsd_+1) + 0;
           int idVely = ir*(nsd_+1) + 1;
           int idVelz = ir*(nsd_+1) + 2;

           for(int ic=0; ic<side_nen_; ic++)
           {
             int iVelx = ic*(nsd_+1)+0;
             int iVely = ic*(nsd_+1)+1;
             int iVelz = ic*(nsd_+1)+2;

             C_uiui_(idVelx, iVelx) += side_side_dyad_timefacfac(ir,ic)*stabfac;
             C_uiui_(idVely, iVely) += side_side_dyad_timefacfac(ir,ic)*stabfac;
             C_uiui_(idVelz, iVelz) += side_side_dyad_timefacfac(ir,ic)*stabfac;
           }

           // -(stab * v2, u2)
           rhC_ui_(idVelx,0) -= side_funct_timefacfac(ir)*stabfac*ivelint(Velx);
           rhC_ui_(idVely,0) -= side_funct_timefacfac(ir)*stabfac*ivelint(Vely);
           rhC_ui_(idVelz,0) -= side_funct_timefacfac(ir)*stabfac*ivelint(Velz);

         }
         } // end coupling


#if 1
        // convective stabilization
         /*                           \        /                       i   _     \
          |  gamma/h_K *  v*n , Du*n     | =  - |   gamma/h_K *  v*n , (u  - u)*n   |
         \                            /        \                                */

        for(int ir=0; ir<nen_; ir++)
        {
          int idVelx = ir*(nsd_+1) + 0;
          int idVely = ir*(nsd_+1) + 1;
          int idVelz = ir*(nsd_+1) + 2;

          // (stab * v, Du)
          for(int ic=0; ic<nen_; ic++)
          {
            int iVelx = ic*(nsd_+1)+0;
            int iVely = ic*(nsd_+1)+1;
            int iVelz = ic*(nsd_+1)+2;

            C_uu_(idVelx, iVelx) += funct_dyad_timefacfac(ir,ic)*stabfac_conv*normal(Velx)*normal(Velx);
            C_uu_(idVelx, iVely) += funct_dyad_timefacfac(ir,ic)*stabfac_conv*normal(Velx)*normal(Vely);
            C_uu_(idVelx, iVelz) += funct_dyad_timefacfac(ir,ic)*stabfac_conv*normal(Velx)*normal(Velz);

            C_uu_(idVely, iVelx) += funct_dyad_timefacfac(ir,ic)*stabfac_conv*normal(Vely)*normal(Velx);
            C_uu_(idVely, iVely) += funct_dyad_timefacfac(ir,ic)*stabfac_conv*normal(Vely)*normal(Vely);
            C_uu_(idVely, iVelz) += funct_dyad_timefacfac(ir,ic)*stabfac_conv*normal(Vely)*normal(Velz);

            C_uu_(idVelz, iVelx) += funct_dyad_timefacfac(ir,ic)*stabfac_conv*normal(Velz)*normal(Velx);
            C_uu_(idVelz, iVely) += funct_dyad_timefacfac(ir,ic)*stabfac_conv*normal(Velz)*normal(Vely);
            C_uu_(idVelz, iVelz) += funct_dyad_timefacfac(ir,ic)*stabfac_conv*normal(Velz)*normal(Velz);
          }

          double velint_normal = velint.Dot(normal);

          // -(stab * v*n, u*n)
          rhs_Cu_(idVelx) -= funct_timefacfac(ir)*normal(Velx)*stabfac_conv*velint_normal;
          rhs_Cu_(idVely) -= funct_timefacfac(ir)*normal(Vely)*stabfac_conv*velint_normal;
          rhs_Cu_(idVelz) -= funct_timefacfac(ir)*normal(Velz)*stabfac_conv*velint_normal;


          double velint_WDBC_normal = ivelint_WDBC_JUMP.Dot(normal);
          // +(stab * v*n, u_DBC*n)
          rhs_Cu_(idVelx) += funct_timefacfac(ir)*normal(Velx)*stabfac_conv*velint_WDBC_normal;
          rhs_Cu_(idVely) += funct_timefacfac(ir)*normal(Vely)*stabfac_conv*velint_WDBC_normal;
          rhs_Cu_(idVelz) += funct_timefacfac(ir)*normal(Velz)*stabfac_conv*velint_WDBC_normal;
        }
#endif



        }


        // set prescribed WDBC at Gaussian point
        virtual void get_vel_WeakDBC( LINALG::Matrix<3,1> & ivelint )
        {
            ivelint.Clear();
            ivelint.Multiply(eivel_,side_funct_);
        }


        LINALG::BlockMatrix<LINALG::Matrix<nen_,side_nen_>,6,3> BK_sui_;
        LINALG::BlockMatrix<LINALG::Matrix<nen_,   1>,6,1>     rhsi_;
        LINALG::BlockMatrix<LINALG::Matrix<side_nen_,nen_>,3,6> BK_uis_;

        LINALG::Matrix<numdof*side_nen_,4*nen_>            C_uiu_;    // row: sidenode1:u1,u2,u3(,p), sidenode2:u1,u2,u3,p ... | col: elenode1:u1,u2,u3,p, elenode2:u1,u2,u3,p, ...
        LINALG::Matrix<4*nen_,numdof*side_nen_>            C_uui_;    // includes ui(,pi)
        LINALG::Matrix<numdof*side_nen_,1>                 rhC_ui_;   // includes (ui,pi)
        LINALG::Matrix<numdof*side_nen_,numdof*side_nen_>  C_uiui_;   // includes (ui,pi) only for Nitsche coupling
        LINALG::Matrix<side_nen_,1>                        side_funct_;
        LINALG::Matrix<nsd_-1,side_nen_>                   side_deriv_;
        LINALG::Matrix<6*nen_,numdof*side_nen_>            K_sui_;
        LINALG::Matrix<numdof*side_nen_,6*nen_>            K_uis_;    // G_uis
        LINALG::Matrix<3,side_nen_>                        eivel_;
        LINALG::Matrix<side_nen_,1>                        eipres_;
        LINALG::Matrix<3,side_nen_>                        eidisp_;
        LINALG::Matrix<3,side_nen_>                        xyze_;
      };

      template<DRT::Element::DiscretizationType distype>
      Teuchos::RCP<SideInterface<distype> > SideInterface<distype>::Impl(DRT::Element * side,
                                                                         Epetra_SerialDenseMatrix &  C_uiu,
                                                                         Epetra_SerialDenseMatrix &  C_uui,
                                                                         Epetra_SerialDenseMatrix &  rhC_ui,
                                                                         Epetra_SerialDenseMatrix &  Gsui,
                                                                         Epetra_SerialDenseMatrix &  Guis,
                                                                         Epetra_SerialDenseMatrix &  side_xyze
                                                                         )
      {
        SideInterface * si = NULL;

        if (side->ElementType() == DRT::ELEMENTS::Bele3Type::Instance()) // three dofs per node, for standard Dirichlet coupling
        {
          const int numdofpernode = 3;

          switch ( side->Shape() )
          {
          case DRT::Element::tri3:
          {
            typedef SideImpl<distype,DRT::Element::tri3, numdofpernode> SideImplType;
            si = new SideImplType(side,C_uiu,C_uui,rhC_ui,Gsui,Guis,side_xyze);
            break;
          }
          case DRT::Element::tri6:
          {
            typedef SideImpl<distype,DRT::Element::tri6, numdofpernode> SideImplType;
            si = new SideImplType(side,C_uiu,C_uui,rhC_ui,Gsui,Guis,side_xyze);
            break;
          }
          case DRT::Element::quad4:
          {
            typedef SideImpl<distype,DRT::Element::quad4, numdofpernode> SideImplType;
            si = new SideImplType(side,C_uiu,C_uui,rhC_ui,Gsui,Guis,side_xyze);
            break;
          }
          case DRT::Element::quad8:
          {
            typedef SideImpl<distype,DRT::Element::quad8, numdofpernode> SideImplType;
            si = new SideImplType(side,C_uiu,C_uui,rhC_ui,Gsui,Guis,side_xyze);
            break;
          }
          case DRT::Element::quad9:
          {
            typedef SideImpl<distype,DRT::Element::quad9, numdofpernode> SideImplType;
            si = new SideImplType(side,C_uiu,C_uui,rhC_ui,Gsui,Guis,side_xyze);
            break;
          }
          default:
            dserror( "unsupported side shape %d", side->Shape() );
          }
        }
        else if (side->ElementType() == DRT::ELEMENTS::Bele3_4Type::Instance()) // four dofs per node, for standard Dirichlet coupling
        {
          const int numdofpernode = 4;

          switch ( side->Shape() )
          {
          case DRT::Element::tri3:
          {
            typedef SideImpl<distype,DRT::Element::tri3, numdofpernode> SideImplType;
            si = new SideImplType(side,C_uiu,C_uui,rhC_ui,Gsui,Guis,side_xyze);
            break;
          }
          case DRT::Element::tri6:
          {
            typedef SideImpl<distype,DRT::Element::tri6, numdofpernode> SideImplType;
            si = new SideImplType(side,C_uiu,C_uui,rhC_ui,Gsui,Guis,side_xyze);
            break;
          }
          case DRT::Element::quad4:
          {
            typedef SideImpl<distype,DRT::Element::quad4, numdofpernode> SideImplType;
            si = new SideImplType(side,C_uiu,C_uui,rhC_ui,Gsui,Guis,side_xyze);
            break;
          }
          case DRT::Element::quad8:
          {
            typedef SideImpl<distype,DRT::Element::quad8, numdofpernode> SideImplType;
            si = new SideImplType(side,C_uiu,C_uui,rhC_ui,Gsui,Guis,side_xyze);
            break;
          }
          case DRT::Element::quad9:
          {
            typedef SideImpl<distype,DRT::Element::quad9, numdofpernode> SideImplType;
            si = new SideImplType(side,C_uiu,C_uui,rhC_ui,Gsui,Guis,side_xyze);
            break;
          }
          default:
            dserror( "unsupported side shape %d", side->Shape() );
          }
        }
        else dserror("unknown boundary element Type!");

        return Teuchos::rcp(si);
      }


      // for Nitsche coupling
      template<DRT::Element::DiscretizationType distype>
      Teuchos::RCP<SideInterface<distype> > SideInterface<distype>::Impl(DRT::Element * side,
                                                                         Epetra_SerialDenseMatrix &  C_uiu,
                                                                         Epetra_SerialDenseMatrix &  C_uui,
                                                                         Epetra_SerialDenseMatrix &  rhC_ui,
                                                                         Epetra_SerialDenseMatrix &  C_uiui,
                                                                         Epetra_SerialDenseMatrix &  side_xyze
                                                                         )
      {
        SideInterface * si = NULL;

        if (side->ElementType() == DRT::ELEMENTS::Bele3Type::Instance()) // three dofs per node, for coupling
        {
          const int numdofpernode = 3;

          switch ( side->Shape() )
          {
          case DRT::Element::tri3:
          {
            typedef SideImpl<distype,DRT::Element::tri3, numdofpernode> SideImplType;
            si = new SideImplType(side,C_uiu,C_uui,rhC_ui,C_uiui,side_xyze);
            break;
          }
          case DRT::Element::tri6:
          {
            typedef SideImpl<distype,DRT::Element::tri6, numdofpernode> SideImplType;
            si = new SideImplType(side,C_uiu,C_uui,rhC_ui,C_uiui,side_xyze);
            break;
          }
          case DRT::Element::quad4:
          {
            typedef SideImpl<distype,DRT::Element::quad4, numdofpernode> SideImplType;
            si = new SideImplType(side,C_uiu,C_uui,rhC_ui,C_uiui,side_xyze);
            break;
          }
          case DRT::Element::quad8:
          {
            typedef SideImpl<distype,DRT::Element::quad8, numdofpernode> SideImplType;
            si = new SideImplType(side,C_uiu,C_uui,rhC_ui,C_uiui,side_xyze);
            break;
          }
          case DRT::Element::quad9:
          {
            typedef SideImpl<distype,DRT::Element::quad9, numdofpernode> SideImplType;
            si = new SideImplType(side,C_uiu,C_uui,rhC_ui,C_uiui,side_xyze);
            break;
          }
          default:
            dserror( "unsupported side shape %d", side->Shape() );
          }
//          return Teuchos::rcp(si);
        }
        else if (side->ElementType() == DRT::ELEMENTS::Bele3_4Type::Instance()) // four dofs per node, for coupling
        {
          const int numdofpernode = 4;

          switch ( side->Shape() )
          {
          case DRT::Element::tri3:
          {
            typedef SideImpl<distype,DRT::Element::tri3, numdofpernode> SideImplType;
            si = new SideImplType(side,C_uiu,C_uui,rhC_ui,C_uiui,side_xyze);
            break;
          }
          case DRT::Element::tri6:
          {
            typedef SideImpl<distype,DRT::Element::tri6, numdofpernode> SideImplType;
            si = new SideImplType(side,C_uiu,C_uui,rhC_ui,C_uiui,side_xyze);
            break;
          }
          case DRT::Element::quad4:
          {
            typedef SideImpl<distype,DRT::Element::quad4, numdofpernode> SideImplType;
            si = new SideImplType(side,C_uiu,C_uui,rhC_ui,C_uiui,side_xyze);
            break;
          }
          case DRT::Element::quad8:
          {
            typedef SideImpl<distype,DRT::Element::quad8, numdofpernode> SideImplType;
            si = new SideImplType(side,C_uiu,C_uui,rhC_ui,C_uiui,side_xyze);
            break;
          }
          case DRT::Element::quad9:
          {
            typedef SideImpl<distype,DRT::Element::quad9, numdofpernode> SideImplType;
            si = new SideImplType(side,C_uiu,C_uui,rhC_ui,C_uiui,side_xyze);
            break;
          }
          default:
            dserror( "unsupported side shape %d", side->Shape() );
          }
//          return Teuchos::rcp(si);
        }
        else dserror("unknown boundary element Type!");

        return Teuchos::rcp(si);
      }


    template<DRT::Element::DiscretizationType distype>
    Teuchos::RCP<SideInterface<distype> > SideInterface<distype>::Impl(DRT::Element * side,
                                                                       Epetra_SerialDenseMatrix &  side_xyze
                                                                       )
    {
      SideInterface * si = NULL;

      if (side->ElementType() == DRT::ELEMENTS::Bele3Type::Instance()) // three dofs per node, for standard Dirichlet coupling
      {
        const int numdofpernode = 3;

        switch ( side->Shape() )
        {
        case DRT::Element::tri3:
        {
          typedef SideImpl<distype,DRT::Element::tri3, numdofpernode> SideImplType;
          si = new SideImplType(side,side_xyze);
          break;
        }
        case DRT::Element::tri6:
        {
          typedef SideImpl<distype,DRT::Element::tri6, numdofpernode> SideImplType;
          si = new SideImplType(side,side_xyze);
          break;
        }
        case DRT::Element::quad4:
        {
          typedef SideImpl<distype,DRT::Element::quad4, numdofpernode> SideImplType;
          si = new SideImplType(side,side_xyze);
          break;
        }
        case DRT::Element::quad8:
        {
          typedef SideImpl<distype,DRT::Element::quad8, numdofpernode> SideImplType;
          si = new SideImplType(side,side_xyze);
          break;
        }
        case DRT::Element::quad9:
        {
          typedef SideImpl<distype,DRT::Element::quad9, numdofpernode> SideImplType;
          si = new SideImplType(side,side_xyze);
          break;
        }
        default:
          dserror( "unsupported side shape %d", side->Shape() );
        }
//        return Teuchos::rcp(si);
      }
      else if (side->ElementType() == DRT::ELEMENTS::Bele3_4Type::Instance()) // three dofs per node, for standard Dirichlet coupling
      {
        const int numdofpernode = 4;

        switch ( side->Shape() )
        {
        case DRT::Element::tri3:
        {
          typedef SideImpl<distype,DRT::Element::tri3, numdofpernode> SideImplType;
          si = new SideImplType(side,side_xyze);
          break;
        }
        case DRT::Element::tri6:
        {
          typedef SideImpl<distype,DRT::Element::tri6, numdofpernode> SideImplType;
          si = new SideImplType(side,side_xyze);
          break;
        }
        case DRT::Element::quad4:
        {
          typedef SideImpl<distype,DRT::Element::quad4, numdofpernode> SideImplType;
          si = new SideImplType(side,side_xyze);
          break;
        }
        case DRT::Element::quad8:
        {
          typedef SideImpl<distype,DRT::Element::quad8, numdofpernode> SideImplType;
          si = new SideImplType(side,side_xyze);
          break;
        }
        case DRT::Element::quad9:
        {
          typedef SideImpl<distype,DRT::Element::quad9, numdofpernode> SideImplType;
          si = new SideImplType(side,side_xyze);
          break;
        }
        default:
          dserror( "unsupported side shape %d", side->Shape() );
        }
//        return Teuchos::rcp(si);
      }
      else dserror("unknown boundary element Type!");

      return Teuchos::rcp(si);
    }



    template<DRT::Element::DiscretizationType distype>
    class EmbCoupling
    {
    public:


      static const int nen_ = DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;
      static const int nsd_ = DRT::UTILS::DisTypeToDim<distype>::dim;


      static Teuchos::RCP<EmbCoupling<distype> > TwoSidedImpl(DRT::Element * emb_ele,
                                                        Epetra_SerialDenseMatrix & C_uiu,
                                                        Epetra_SerialDenseMatrix & C_uui,
                                                        Epetra_SerialDenseMatrix & rhC_ui,
                                                        Epetra_SerialDenseMatrix & C_uiui,
                                                        Epetra_SerialDenseMatrix & emb_xyze
                                                        );


      virtual ~EmbCoupling() {}

      virtual void EvaluateEmb( LINALG::Matrix<nsd_,1> & xside ) = 0;



      virtual void emb_vel(const DRT::Discretization &  cutdis,
                         const std::string            state,
                         const vector<int>&           lm) = 0;


      virtual void addembdisp(const DRT::Discretization &  cutdis,
                             const std::string            state,
                             const vector<int>&           lm,
                             Epetra_SerialDenseMatrix  &  side_xyze) = 0;

      virtual void element_length( double & hk_emb ) = 0;


      virtual void buildCouplingMatricesNitscheTwoSided(   Epetra_SerialDenseMatrix &    C_uu_,          // standard bg-bg-matrix
                                                   Epetra_SerialDenseVector &    rhs_Cu_,        // standard bg-rhs
                                                   bool &                        coupling,       // assemble coupling terms (yes/no)
                                                   bool &                        bg_mortaring,   // yes: background-sided mortaring, no: coupling between two meshes (mixed mortaring)
                                                   LINALG::Matrix<nsd_,1> &      normal,         // normal vector
                                                   const double                  timefacfac,     // theta*dt
                                                   const double                  visceff_1,      // viscosity in background fluid
                                                   const double                  visceff_2,      // viscosity in embedded fluid
                                                   double &                      kappa1,         // mortaring weighting
                                                   double &                      kappa2,         // mortaring weighting
                                                   double &                      stabfac,        // Nitsche non-dimensionless stabilization factor
                                                   double &                      stabfac_conv,   // Nitsche convective non-dimensionless stabilization factor
                                                   LINALG::Matrix<nen_,1> &      funct_,         // bg shape functions
                                                   LINALG::Matrix<nsd_,nen_> &   derxy_,         // bg deriv
                                                   LINALG::Matrix<nsd_,nsd_> &   vderxy_,        // bg deriv^n
                                                   double &                      press,          // bg p^n
                                                   LINALG::Matrix<nsd_,1> &      velint,          // bg u^n
                                                   LINALG::Matrix<nsd_,1> &      ivelint_WDBC_JUMP // Dirichlet velocity vector or prescribed jump vector
                                                ) = 0;

      virtual void get_vel_WeakDBC (LINALG::Matrix<nsd_,1> & emb_velint) = 0;

    };


    template<DRT::Element::DiscretizationType distype,
             DRT::Element::DiscretizationType emb_distype>
    class EmbImpl : public EmbCoupling<distype>
    {
    public:

      //! nen_: number of element nodes (P. Hughes: The Finite Element Method)
      static const int emb_nen_ = DRT::UTILS::DisTypeToNumNodePerEle<emb_distype>::numNodePerElement;
      static const int nen_ = DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;
      static const int nsd_ = DRT::UTILS::DisTypeToDim<distype>::dim;


      // for Nitsche-based fluid-fluid coupling
      EmbImpl(DRT::Element * emb_ele,
              Epetra_SerialDenseMatrix & C_uiu,
              Epetra_SerialDenseMatrix & C_uui,
              Epetra_SerialDenseMatrix & rhC_ui,
              Epetra_SerialDenseMatrix & C_uiui,
              Epetra_SerialDenseMatrix & emb_xyze
               )
        : C_uiu_(C_uiu.A(),true),
          C_uui_(C_uui.A(),true),
          rhC_ui_(rhC_ui.A(),true),
          C_uiui_(C_uiui.A(),true),
          emb_xyze_(emb_xyze.A(),true)
      {
      }


      virtual void EvaluateEmb( LINALG::Matrix<nsd_,1> & xside )
      {
        // find element local position of gauss point
        GEO::CUT::Position<emb_distype> pos( emb_xyze_, xside );
        pos.Compute();

        const LINALG::Matrix<3,1> & rst_emb = pos.LocalCoordinates();

        DRT::UTILS::shape_function_3D( emb_funct_, rst_emb( 0 ), rst_emb( 1 ), rst_emb( 2 ), emb_distype );
        DRT::UTILS::shape_function_3D_deriv1( emb_deriv_, rst_emb( 0 ), rst_emb( 1 ), rst_emb( 2 ), emb_distype );


        LINALG::Matrix<nsd_,nsd_> emb_xjm(true);
        LINALG::Matrix<nsd_,nsd_> emb_xji(true);

        emb_xjm.MultiplyNT(emb_deriv_,emb_xyze_);
        emb_xji.Invert(emb_xjm);

        // compute global first derivates
        emb_derxy_.Multiply(emb_xji,emb_deriv_);

        // get velocity derivatives at integration point
        // (values at n+alpha_F for generalized-alpha scheme, n+1 otherwise)
        emb_vderxy_.MultiplyNT(emb_vel_,emb_derxy_);


      }



      virtual void emb_vel(const DRT::Discretization &  embdis,
                         const std::string            state,
                         const vector<int>&           lm)
      {
        // get state of the global vector
        Teuchos::RCP<const Epetra_Vector> matrix_state = embdis.GetState(state);
        if(matrix_state == null)
          dserror("Cannot get state vector %s", state.c_str());

        // extract local values of the global vectors
        std::vector<double> mymatrix(lm.size());
        DRT::UTILS::ExtractMyValues(*matrix_state,mymatrix,lm);

        for (int inode=0; inode<emb_nen_; ++inode)  // number of nodes
        {
          for(int idim=0; idim<3; ++idim) // number of dimensions
          {
            (emb_vel_)(idim,inode) = mymatrix[idim+(inode*4)];  // state vector includes velocity and pressure
          }
          (emb_pres_)(inode,0) = mymatrix[3+(inode*4)];
        }
      }

      virtual void addembdisp(const DRT::Discretization &  embdis,
                             const std::string            state,
                             const vector<int>&           lm,
                             Epetra_SerialDenseMatrix  &  side_xyze)
      {
        // get state of the global vector
        Teuchos::RCP<const Epetra_Vector> matrix_state = embdis.GetState(state);
        if(matrix_state == null)
          dserror("Cannot get state vector %s", state.c_str());

        // extract local values of the global vectors
        std::vector<double> mymatrix(lm.size());
        DRT::UTILS::ExtractMyValues(*matrix_state,mymatrix,lm);

        for (int inode=0; inode<emb_nen_; ++inode)  // number of nodes
        {
          for(int idim=0; idim<3; ++idim) // number of dimensions
          {
            (emb_disp_)(idim,inode) = mymatrix[idim+(inode*4)]; // attention! disp state vector has 3+1 dofs for displacement (the same as for (u,p))
          }
        }

        // add the displacement of the interface
        for (int inode = 0; inode < emb_nen_; ++inode)
        {
          emb_xyze_(0,inode) += emb_disp_(0, inode);
          emb_xyze_(1,inode) += emb_disp_(1, inode);
          emb_xyze_(2,inode) += emb_disp_(2, inode);
        }

      }


      virtual void element_length( double & hk_emb )
      {
        LINALG::Matrix<1,1> dummy(true);
        hk_emb = FLD::UTILS::HK<emb_distype>(dummy, emb_xyze_);
      }


      virtual void buildCouplingMatricesNitscheTwoSided(  Epetra_SerialDenseMatrix &    C_uu_,          // standard bg-bg-matrix
                                                  Epetra_SerialDenseVector &    rhs_Cu_,        // standard bg-rhs
                                                  bool &                        coupling,       // assemble coupling terms (yes/no)
                                                  bool &                        bg_mortaring,   // yes: background-sided mortaring, no: coupling between two meshes (mixed mortaring)
                                                  LINALG::Matrix<nsd_,1> &      normal,         // normal vector
                                                  const double                  timefacfac,     // theta*dt
                                                  const double                  visceff_1,      // viscosity in background fluid
                                                  const double                  visceff_2,      // viscosity in embedded fluid
                                                  double &                      kappa1,         // mortaring weighting
                                                  double &                      kappa2,         // mortaring weighting
                                                  double &                      stabfac,        // Nitsche non-dimensionless stabilization factor
                                                  double &                      stabfac_conv,   // Nitsche convective non-dimensionless stabilization factor
                                                  LINALG::Matrix<nen_,1> &      funct_,         // bg shape functions
                                                  LINALG::Matrix<nsd_,nen_> &   derxy_,         // bg deriv
                                                  LINALG::Matrix<nsd_,nsd_> &   vderxy_,        // bg deriv^n
                                                  double &                      press,          // bg p^n
                                                  LINALG::Matrix<nsd_,1> &      velint,         // bg u^n
                                                  LINALG::Matrix<nsd_,1> &      ivelint_WDBC_JUMP // Dirichlet velocity vector or prescribed jump vector
                                             )
      {

        const unsigned Velx = 0;
        const unsigned Vely = 1;
        const unsigned Velz = 2;
        //const unsigned Pres = 3;

        //--------------------------------------------

        // define the coupling between two not matching grids
        // for fluidfluidcoupling
        // domain Omega^1 := Xfluid
        // domain Omega^2 := Alefluid( or monolithic: structure) ( not available for non-coupling (Dirichlet) )

        // [| v |] := v1 - v2
        //  { v }  := kappa1 * v1 + kappa2 * v2 = kappa1 * v1 (for Dirichlet coupling k1=1.0, k2 = 0.0)
        //  < v >  := kappa2 * v1 + kappa1 * v2 = kappa1 * v2 (for Dirichlet coupling k1=1.0, k2 = 0.0)


        double k1mu1_fac = 2.0 * timefacfac * kappa1 * visceff_1;
        double k2mu2_fac = 2.0 * timefacfac * kappa2 * visceff_2;


        //--------------------------------------------
        // get fields at interface


        // get velocity at integration point
        // (values at n+alpha_F for generalized-alpha scheme, n+1 otherwise)
        // interface velocity vector in gausspoint
        LINALG::Matrix<nsd_,1> emb_velint(true);
        emb_velint.Multiply(emb_vel_,emb_funct_);

        double emb_press = emb_funct_.Dot(emb_pres_);

        // funct_ * timefac * fac * funct_ (dyadic product)
        LINALG::Matrix<nen_,1> funct_timefacfac(true);
        funct_timefacfac.Update(timefacfac,funct_,0.0);

        LINALG::Matrix<emb_nen_,1> emb_funct_timefacfac(true);
        emb_funct_timefacfac.Update(timefacfac,emb_funct_,0.0);

        LINALG::Matrix<emb_nen_,emb_nen_> emb_emb_dyad_timefacfac(true);
        emb_emb_dyad_timefacfac.MultiplyNT(emb_funct_timefacfac, emb_funct_);

        LINALG::Matrix<emb_nen_,emb_nen_> emb_emb_dyad_k2_timefacfac(true);
        emb_emb_dyad_k2_timefacfac.Update(kappa2, emb_emb_dyad_timefacfac, 0.0);


        LINALG::Matrix<emb_nen_,nen_> emb_funct_dyad_timefacfac(true);
        LINALG::Matrix<emb_nen_,nen_> emb_funct_dyad_k1_timefacfac(true);
        emb_funct_dyad_timefacfac.MultiplyNT(emb_funct_timefacfac, funct_);
        emb_funct_dyad_k1_timefacfac.Update(kappa1, emb_funct_dyad_timefacfac, 0.0);


        LINALG::Matrix<emb_nen_,nen_> emb_funct_dyad_k2_timefacfac(true);
        emb_funct_dyad_k2_timefacfac.Update(kappa2, emb_funct_dyad_timefacfac, 0.0);


                   /*                  \       /          i      \
                + |  [ v ], + {Dp}*n    | = - | [ v ], { p }* n   |
                   \                   /       \                */

        //-----------------------------------------------
        //    + (v1, k1 *(Dp1)*n)
        //-----------------------------------------------
        LINALG::Matrix<nen_,nen_> funct_dyad_timefacfac(true);
        LINALG::Matrix<nen_,nen_> funct_dyad_k1_timefacfac(true);
        funct_dyad_timefacfac.MultiplyNT(funct_timefacfac, funct_);
        funct_dyad_k1_timefacfac.Update(kappa1,funct_dyad_timefacfac,0.0);

        for(int ir = 0; ir<nen_; ir++)
        {
          int idVelx = ir*(nsd_+1) + 0;
          int idVely = ir*(nsd_+1) + 1;
          int idVelz = ir*(nsd_+1) + 2;

          // (v,Dp*n)
          for(int ic =0; ic<nen_; ic++)
          {
            int iPres = ic*(nsd_+1)+3;

            C_uu_(idVelx, iPres) += funct_dyad_k1_timefacfac(ir,ic)*normal(Velx);
            C_uu_(idVely, iPres) += funct_dyad_k1_timefacfac(ir,ic)*normal(Vely);
            C_uu_(idVelz, iPres) += funct_dyad_k1_timefacfac(ir,ic)*normal(Velz);
          }

          // -(v,p*n)
          double funct_k1_timefacfac_press = funct_timefacfac(ir)*press*kappa1;
          rhs_Cu_(idVelx,0) -= funct_k1_timefacfac_press*normal(Velx);
          rhs_Cu_(idVely,0) -= funct_k1_timefacfac_press*normal(Vely);
          rhs_Cu_(idVelz,0) -= funct_k1_timefacfac_press*normal(Velz);
        }



        if(coupling)
        {
        //-----------------------------------------------
        //    - (v2, k1 *(Dp1)*n)
        //-----------------------------------------------
        for(int ir = 0; ir<emb_nen_; ir++)
        {
          int idVelx = ir*(nsd_+1) + 0;
          int idVely = ir*(nsd_+1) + 1;
          int idVelz = ir*(nsd_+1) + 2;

          // (v,Dp*n)
          for(int ic =0; ic<nen_; ic++)
          {
            int iPres = ic*(nsd_+1)+3;

            C_uiu_(idVelx, iPres) -= emb_funct_dyad_k1_timefacfac(ir,ic)*normal(Velx);
            C_uiu_(idVely, iPres) -= emb_funct_dyad_k1_timefacfac(ir,ic)*normal(Vely);
            C_uiu_(idVelz, iPres) -= emb_funct_dyad_k1_timefacfac(ir,ic)*normal(Velz);
          }

          // -(v,p*n)
          double emb_funct_k1_timefacfac_press = emb_funct_timefacfac(ir)*press*kappa1;
          rhC_ui_(idVelx,0) += emb_funct_k1_timefacfac_press*normal(Velx);
          rhC_ui_(idVely,0) += emb_funct_k1_timefacfac_press*normal(Vely);
          rhC_ui_(idVelz,0) += emb_funct_k1_timefacfac_press*normal(Velz);
        }


        if(!bg_mortaring)
        {
          LINALG::Matrix<emb_nen_,nen_> emb_funct_dyad_k2_timefacfac(true);
          emb_funct_dyad_k2_timefacfac.Update(kappa2, emb_funct_dyad_timefacfac, 0.0);

          //-----------------------------------------------
          //    + (v1, k2 *(Dp2)*n)
          //-----------------------------------------------

          for(int ir = 0; ir<nen_; ir++)
          {
            int idVelx = ir*(nsd_+1) + 0;
            int idVely = ir*(nsd_+1) + 1;
            int idVelz = ir*(nsd_+1) + 2;

            // (v,Dp*n)
            for(int ic =0; ic<emb_nen_; ic++)
            {
              int iPres = ic*(nsd_+1)+3;

              C_uui_(idVelx, iPres) += emb_funct_dyad_k2_timefacfac(ic,ir)*normal(Velx);
              C_uui_(idVely, iPres) += emb_funct_dyad_k2_timefacfac(ic,ir)*normal(Vely);
              C_uui_(idVelz, iPres) += emb_funct_dyad_k2_timefacfac(ic,ir)*normal(Velz);
            }

            // -(v,p*n)
            double funct_k2_timefacfac_press = funct_timefacfac(ir)*emb_press*kappa2;
            rhs_Cu_(idVelx,0) -= funct_k2_timefacfac_press*normal(Velx);
            rhs_Cu_(idVely,0) -= funct_k2_timefacfac_press*normal(Vely);
            rhs_Cu_(idVelz,0) -= funct_k2_timefacfac_press*normal(Velz);
          }




          //-----------------------------------------------
          //    - (v2, k1 *(Dp2)*n)
          //-----------------------------------------------
          for(int ir = 0; ir<emb_nen_; ir++)
          {
            int idVelx = ir*(nsd_+1) + 0;
            int idVely = ir*(nsd_+1) + 1;
            int idVelz = ir*(nsd_+1) + 2;

            // (v,Dp*n)
            for(int ic =0; ic<nen_; ic++)
            {
              int iPres = ic*(nsd_+1)+3;

              C_uiui_(idVelx, iPres) -= emb_emb_dyad_k2_timefacfac(ir,ic)*normal(Velx);
              C_uiui_(idVely, iPres) -= emb_emb_dyad_k2_timefacfac(ir,ic)*normal(Vely);
              C_uiui_(idVelz, iPres) -= emb_emb_dyad_k2_timefacfac(ir,ic)*normal(Velz);
            }

            // -(v,p*n)
            double emb_funct_k2_timefacfac_press = emb_funct_timefacfac(ir)*emb_press*kappa2;
            rhC_ui_(idVelx,0) += emb_funct_k2_timefacfac_press*normal(Velx);
            rhC_ui_(idVely,0) += emb_funct_k2_timefacfac_press*normal(Vely);
            rhC_ui_(idVelz,0) += emb_funct_k2_timefacfac_press*normal(Velz);
          }

        } // end !bgmortaring



        }// end coupling



            /*                   \     /              i   \
         - |  { q }*n ,[ Du ]     | = |  { q }*n  ,[ u ]   |
            \                    /     \                 */
        double alpha_p = +1.0; // (+1) antisymmetric, (-1) symmetric

       //-----------------------------------------------
       //    - (q1*n, k1 *(Du1))
       //-----------------------------------------------
       for(int ir = 0; ir<nen_; ir++)
       {
        int idPres = ir*(nsd_+1) + 3;

        // -(q*n,Du)
        for(int ic =0; ic<nen_; ic++)
        {
          int iVelx = ic*(nsd_+1)+0;
          int iVely = ic*(nsd_+1)+1;
          int iVelz = ic*(nsd_+1)+2;

          C_uu_(idPres, iVelx) -= alpha_p*funct_dyad_k1_timefacfac(ir,ic)*normal(Velx);
          C_uu_(idPres, iVely) -= alpha_p*funct_dyad_k1_timefacfac(ir,ic)*normal(Vely);
          C_uu_(idPres, iVelz) -= alpha_p*funct_dyad_k1_timefacfac(ir,ic)*normal(Velz);
        }

        // (q*n,u)
        double velint_normal = velint.Dot(normal);
        rhs_Cu_(idPres,0) += alpha_p*funct_timefacfac(ir)*kappa1*velint_normal;

//        if(!coupling) // weak Dirichlet case
//        {
//          // -(q*n,u_DBC)
//          double ivelint_WDBC_JUMP_normal = ivelint_WDBC_JUMP.Dot(normal);
//          rhs_Cu_(idPres,0) -= funct_timefacfac(ir)*ivelint_WDBC_JUMP_normal;
//        }
       }


       if(coupling)
       {
       //-----------------------------------------------
       //    + (q1*n, k1 *(Du2))
       //-----------------------------------------------
       for(int ir = 0; ir<nen_; ir++)
       {
        int idPres = ir*(nsd_+1) + 3;

        // -(q*n,Du)
        for(int ic =0; ic<emb_nen_; ic++)
        {
          int iVelx = ic*(nsd_+1)+0;
          int iVely = ic*(nsd_+1)+1;
          int iVelz = ic*(nsd_+1)+2;

          C_uui_(idPres, iVelx) += alpha_p*emb_funct_dyad_k1_timefacfac(ic,ir)*normal(Velx);
          C_uui_(idPres, iVely) += alpha_p*emb_funct_dyad_k1_timefacfac(ic,ir)*normal(Vely);
          C_uui_(idPres, iVelz) += alpha_p*emb_funct_dyad_k1_timefacfac(ic,ir)*normal(Velz);
        }

        // -(q*n,u)
        double emb_velint_normal = emb_velint.Dot(normal);
        rhs_Cu_(idPres,0) -= alpha_p*funct_timefacfac(ir)*kappa1*emb_velint_normal;

       }
       }// end coupling

       if(!bg_mortaring)
       {
         //-----------------------------------------------
         //    - (q2*n, k2 *(Du1))
         //-----------------------------------------------
         for(int ir = 0; ir<emb_nen_; ir++)
         {
          int idPres = ir*(nsd_+1) + 3;

          // -(q*n,Du)
          for(int ic =0; ic<nen_; ic++)
          {
            int iVelx = ic*(nsd_+1)+0;
            int iVely = ic*(nsd_+1)+1;
            int iVelz = ic*(nsd_+1)+2;

            C_uiu_(idPres, iVelx) -= alpha_p*emb_funct_dyad_k2_timefacfac(ir,ic)*normal(Velx);
            C_uiu_(idPres, iVely) -= alpha_p*emb_funct_dyad_k2_timefacfac(ir,ic)*normal(Vely);
            C_uiu_(idPres, iVelz) -= alpha_p*emb_funct_dyad_k2_timefacfac(ir,ic)*normal(Velz);
          }

          // (q*n,u)
          double velint_normal = velint.Dot(normal);
          rhC_ui_(idPres,0) += alpha_p*emb_funct_timefacfac(ir)*kappa2*velint_normal;


         }


         //-----------------------------------------------
         //    + (q2*n, k2 *(Du2))
         //-----------------------------------------------
         for(int ir = 0; ir<emb_nen_; ir++)
         {
          int idPres = ir*(nsd_+1) + 3;

          // -(q*n,Du)
          for(int ic =0; ic<emb_nen_; ic++)
          {
            int iVelx = ic*(nsd_+1)+0;
            int iVely = ic*(nsd_+1)+1;
            int iVelz = ic*(nsd_+1)+2;

            C_uiui_(idPres, iVelx) += alpha_p*emb_emb_dyad_k2_timefacfac(ic,ir)*normal(Velx);
            C_uiui_(idPres, iVely) += alpha_p*emb_emb_dyad_k2_timefacfac(ic,ir)*normal(Vely);
            C_uiui_(idPres, iVelz) += alpha_p*emb_emb_dyad_k2_timefacfac(ic,ir)*normal(Velz);
          }

          // -(q*n,u)
          double emb_velint_normal = emb_velint.Dot(normal);
          rhC_ui_(idPres,0) -= alpha_p*emb_funct_timefacfac(ir)*kappa2*emb_velint_normal;

         }

       }





        /*                           \       /                   i      \
     - |  [ v ],  { 2mu eps(u) }*n    | = + | [ v ],  { 2mu eps(u ) }*n  |
        \                            /       \                         */

        //-----------------------------------------------
        //    - (v1, (2*k1*mu1) *eps(Du1)*n)
        //-----------------------------------------------


        LINALG::Matrix<nen_,1> e_funct_visc1_timefacfac(true);
        e_funct_visc1_timefacfac.Update(k1mu1_fac, funct_, 0.0);

        LINALG::Matrix<nen_,1> e_funct_visc2_timefacfac(true);
        e_funct_visc2_timefacfac.Update(k2mu2_fac, funct_, 0.0);

        //            LINALG::Matrix<emb_nen_,1> s_funct_visc_timefacfac(true);
        //            s_funct_visc_timefacfac.Update(k2mu2_fac, emb_funct_, 0.0);

        for(int ir = 0; ir<nen_; ir++)
        {
          int idVelx = ir*(nsd_+1) + 0;
          int idVely = ir*(nsd_+1) + 1;
          int idVelz = ir*(nsd_+1) + 2;

          for(int ic =0; ic<nen_; ic++)
          {
            int iVelx = ic*(nsd_+1)+0;
            int iVely = ic*(nsd_+1)+1;
            int iVelz = ic*(nsd_+1)+2;


            // - (v1, (2*k1*mu1) *eps(Du1)*n)

            //(x,x)
            C_uu_(idVelx, iVelx) -= e_funct_visc1_timefacfac(ir)*(         normal(Velx)*derxy_(Velx,ic)
                                                                   + 0.5 * normal(Vely)*derxy_(Vely,ic)
                                                                   + 0.5 * normal(Velz)*derxy_(Velz,ic)  );
            //(x,y)
            C_uu_(idVelx, iVely) -= e_funct_visc1_timefacfac(ir)*    0.5 * normal(Vely)*derxy_(Velx,ic);
            //(x,z)
            C_uu_(idVelx, iVelz) -= e_funct_visc1_timefacfac(ir)*    0.5 * normal(Velz)*derxy_(Velx,ic);

            //(y,x)
            C_uu_(idVely, iVelx) -= e_funct_visc1_timefacfac(ir)*    0.5 * normal(Velx)*derxy_(Vely,ic);
            //(y,y)
            C_uu_(idVely, iVely) -= e_funct_visc1_timefacfac(ir)*(   0.5 * normal(Velx)*derxy_(Velx,ic)
                                                                   +       normal(Vely)*derxy_(Vely,ic)
                                                                   + 0.5 * normal(Velz)*derxy_(Velz,ic)  );
            //(y,z)
            C_uu_(idVely, iVelz) -= e_funct_visc1_timefacfac(ir)*    0.5 * normal(Velz)*derxy_(Vely,ic);

            //(z,x)
            C_uu_(idVelz, iVelx) -= e_funct_visc1_timefacfac(ir)*    0.5 * normal(Velx)*derxy_(Velz,ic);
            //(z,y)
            C_uu_(idVelz, iVely) -= e_funct_visc1_timefacfac(ir)*    0.5 * normal(Vely)*derxy_(Velz,ic);
            //(z,z)
            C_uu_(idVelz, iVelz) -= e_funct_visc1_timefacfac(ir)*(   0.5 * normal(Velx)*derxy_(Velx,ic)
                                                                   + 0.5 * normal(Vely)*derxy_(Vely,ic)
                                                                   +       normal(Velz)*derxy_(Velz,ic)  );
          }

          // - (v1, (2*k1*mu1) *eps(Du1)*n)
          rhs_Cu_(idVelx,0) += e_funct_visc1_timefacfac(ir)*(            vderxy_(Velx,Velx)                      *normal(Velx)
                                                               + 0.5 * ( vderxy_(Velx,Vely) + vderxy_(Vely,Velx))*normal(Vely)
                                                               + 0.5 * ( vderxy_(Velx,Velz) + vderxy_(Velz,Velx))*normal(Velz)  );
          rhs_Cu_(idVely,0) += e_funct_visc1_timefacfac(ir)*(    0.5 * ( vderxy_(Vely,Velx) + vderxy_(Velx,Vely))*normal(Velx)
                                                               +         vderxy_(Vely,Vely)                      *normal(Vely)
                                                               + 0.5 * ( vderxy_(Vely,Velz) + vderxy_(Velz,Vely))*normal(Velz)  );
          rhs_Cu_(idVelz,0) += e_funct_visc1_timefacfac(ir)*(    0.5 * ( vderxy_(Velz,Velx) + vderxy_(Velx,Velz))*normal(Velx)
                                                               + 0.5 * ( vderxy_(Velz,Vely) + vderxy_(Vely,Velz))*normal(Vely)
                                                               +         vderxy_(Velz,Velz)                      *normal(Velz)  );

        }


        LINALG::Matrix<emb_nen_,1> s_funct_visc1_timefacfac(true);
        s_funct_visc1_timefacfac.Update(k1mu1_fac, emb_funct_, 0.0);
        LINALG::Matrix<emb_nen_,1> s_funct_visc2_timefacfac(true);
        s_funct_visc2_timefacfac.Update(k2mu2_fac, emb_funct_, 0.0);
        if(coupling)
        {
        //-----------------------------------------------
        //    + (v2, (2*k1*mu1) *eps(Du1)*n)
        //-----------------------------------------------
        for(int ir = 0; ir<emb_nen_; ir++)
        {
          int idVelx = ir*(nsd_+1) + 0;
          int idVely = ir*(nsd_+1) + 1;
          int idVelz = ir*(nsd_+1) + 2;
          //int idPres = ir*(nsd_+1) + 3;


          for(int ic =0; ic<nen_; ic++)
          {
            int iVelx = ic*(nsd_+1)+0;
            int iVely = ic*(nsd_+1)+1;
            int iVelz = ic*(nsd_+1)+2;
            //int iPres = ic*(nsd_+1)+3;


            // + (v2, (2*k1*mu1) *eps(Du1)*n)

            //(x,x)
            C_uiu_(idVelx, iVelx) += s_funct_visc1_timefacfac(ir)*(         normal(Velx)*derxy_(Velx,ic)
                                                                    + 0.5 * normal(Vely)*derxy_(Vely,ic)
                                                                    + 0.5 * normal(Velz)*derxy_(Velz,ic)  );
            //(x,y)
            C_uiu_(idVelx, iVely) += s_funct_visc1_timefacfac(ir)*    0.5 * normal(Vely)*derxy_(Velx,ic);
            //(x,z)
            C_uiu_(idVelx, iVelz) += s_funct_visc1_timefacfac(ir)*    0.5 * normal(Velz)*derxy_(Velx,ic);

            //(y,x)
            C_uiu_(idVely, iVelx) += s_funct_visc1_timefacfac(ir)*    0.5 * normal(Velx)*derxy_(Vely,ic);
            //(y,y)
            C_uiu_(idVely, iVely) += s_funct_visc1_timefacfac(ir)*(   0.5 * normal(Velx)*derxy_(Velx,ic)
                                                                    +       normal(Vely)*derxy_(Vely,ic)
                                                                    + 0.5 * normal(Velz)*derxy_(Velz,ic)  );
            //(y,z)
            C_uiu_(idVely, iVelz) += s_funct_visc1_timefacfac(ir)*    0.5 * normal(Velz)*derxy_(Vely,ic);

            //(z,x)
            C_uiu_(idVelz, iVelx) += s_funct_visc1_timefacfac(ir)*    0.5 * normal(Velx)*derxy_(Velz,ic);
            //(z,y)
            C_uiu_(idVelz, iVely) += s_funct_visc1_timefacfac(ir)*    0.5 * normal(Vely)*derxy_(Velz,ic);
            //(z,z)
            C_uiu_(idVelz, iVelz) += s_funct_visc1_timefacfac(ir)*(   0.5 * normal(Velx)*derxy_(Velx,ic)
                                                                    + 0.5 * normal(Vely)*derxy_(Vely,ic)
                                                                    +       normal(Velz)*derxy_(Velz,ic)  );
          }

          // - (v2, (2*k1*mu1) *eps(Du1)*n)
          rhC_ui_(idVelx,0) -= s_funct_visc1_timefacfac(ir)*(            vderxy_(Velx,Velx)                      *normal(Velx)
                                                               + 0.5 * ( vderxy_(Velx,Vely) + vderxy_(Vely,Velx))*normal(Vely)
                                                               + 0.5 * ( vderxy_(Velx,Velz) + vderxy_(Velz,Velx))*normal(Velz)  );
          rhC_ui_(idVely,0) -= s_funct_visc1_timefacfac(ir)*(    0.5 * ( vderxy_(Vely,Velx) + vderxy_(Velx,Vely))*normal(Velx)
                                                               +         vderxy_(Vely,Vely)                      *normal(Vely)
                                                               + 0.5 * ( vderxy_(Vely,Velz) + vderxy_(Velz,Vely))*normal(Velz)  );
          rhC_ui_(idVelz,0) -= s_funct_visc1_timefacfac(ir)*(    0.5 * ( vderxy_(Velz,Velx) + vderxy_(Velx,Velz))*normal(Velx)
                                                               + 0.5 * ( vderxy_(Velz,Vely) + vderxy_(Vely,Velz))*normal(Vely)
                                                               +         vderxy_(Velz,Velz)                      *normal(Velz)  );

        }
        }// end coupling

        if(!bg_mortaring)
        {
          //-----------------------------------------------
          //    - (v1, (2*k2*mu2) *eps(Du2)*n)
          //-----------------------------------------------


          LINALG::Matrix<nen_,1> e_funct_visc2_timefacfac(true);
          e_funct_visc2_timefacfac.Update(k2mu2_fac, funct_, 0.0);


          for(int ir = 0; ir<nen_; ir++)
          {
            int idVelx = ir*(nsd_+1) + 0;
            int idVely = ir*(nsd_+1) + 1;
            int idVelz = ir*(nsd_+1) + 2;


            for(int ic =0; ic<emb_nen_; ic++)
            {
              int iVelx = ic*(nsd_+1)+0;
              int iVely = ic*(nsd_+1)+1;
              int iVelz = ic*(nsd_+1)+2;

              // - (v1, (2*k2*mu2) *eps(Du2)*n)

              //(x,x)
              C_uui_(idVelx, iVelx) -= e_funct_visc2_timefacfac(ir)*(         normal(Velx)*emb_derxy_(Velx,ic)
                                                                     + 0.5 * normal(Vely)*emb_derxy_(Vely,ic)
                                                                     + 0.5 * normal(Velz)*emb_derxy_(Velz,ic)  );
              //(x,y)
              C_uui_(idVelx, iVely) -= e_funct_visc2_timefacfac(ir)*    0.5 * normal(Vely)*emb_derxy_(Velx,ic);
              //(x,z)
              C_uui_(idVelx, iVelz) -= e_funct_visc2_timefacfac(ir)*    0.5 * normal(Velz)*emb_derxy_(Velx,ic);

              //(y,x)
              C_uui_(idVely, iVelx) -= e_funct_visc2_timefacfac(ir)*    0.5 * normal(Velx)*emb_derxy_(Vely,ic);
              //(y,y)
              C_uui_(idVely, iVely) -= e_funct_visc2_timefacfac(ir)*(   0.5 * normal(Velx)*emb_derxy_(Velx,ic)
                                                                     +       normal(Vely)*emb_derxy_(Vely,ic)
                                                                     + 0.5 * normal(Velz)*emb_derxy_(Velz,ic)  );
              //(y,z)
              C_uui_(idVely, iVelz) -= e_funct_visc2_timefacfac(ir)*    0.5 * normal(Velz)*emb_derxy_(Vely,ic);

              //(z,x)
              C_uui_(idVelz, iVelx) -= e_funct_visc2_timefacfac(ir)*    0.5 * normal(Velx)*emb_derxy_(Velz,ic);
              //(z,y)
              C_uui_(idVelz, iVely) -= e_funct_visc2_timefacfac(ir)*    0.5 * normal(Vely)*emb_derxy_(Velz,ic);
              //(z,z)
              C_uui_(idVelz, iVelz) -= e_funct_visc2_timefacfac(ir)*(   0.5 * normal(Velx)*emb_derxy_(Velx,ic)
                                                                     + 0.5 * normal(Vely)*emb_derxy_(Vely,ic)
                                                                     +       normal(Velz)*emb_derxy_(Velz,ic)  );
            }

            // - (v1, (2*k2*mu2) *eps(u2)*n)
            rhs_Cu_(idVelx,0) += e_funct_visc2_timefacfac(ir)*(            emb_vderxy_(Velx,Velx)                          *normal(Velx)
                                                                 + 0.5 * ( emb_vderxy_(Velx,Vely) + emb_vderxy_(Vely,Velx))*normal(Vely)
                                                                 + 0.5 * ( emb_vderxy_(Velx,Velz) + emb_vderxy_(Velz,Velx))*normal(Velz)  );
            rhs_Cu_(idVely,0) += e_funct_visc2_timefacfac(ir)*(    0.5 * ( emb_vderxy_(Vely,Velx) + emb_vderxy_(Velx,Vely))*normal(Velx)
                                                                 +         emb_vderxy_(Vely,Vely)                          *normal(Vely)
                                                                 + 0.5 * ( emb_vderxy_(Vely,Velz) + emb_vderxy_(Velz,Vely))*normal(Velz)  );
            rhs_Cu_(idVelz,0) += e_funct_visc2_timefacfac(ir)*(    0.5 * ( emb_vderxy_(Velz,Velx) + emb_vderxy_(Velx,Velz))*normal(Velx)
                                                                 + 0.5 * ( emb_vderxy_(Velz,Vely) + emb_vderxy_(Vely,Velz))*normal(Vely)
                                                                 +         emb_vderxy_(Velz,Velz)                          *normal(Velz)  );

          }


          LINALG::Matrix<emb_nen_,1> s_funct_visc2_timefacfac(true);
          s_funct_visc2_timefacfac.Update(k2mu2_fac, emb_funct_, 0.0);

          //-----------------------------------------------
          //    + (v2, (2*k2*mu2) *eps(Du2)*n)
          //-----------------------------------------------
          for(int ir = 0; ir<emb_nen_; ir++)
          {
            int idVelx = ir*(nsd_+1) + 0;
            int idVely = ir*(nsd_+1) + 1;
            int idVelz = ir*(nsd_+1) + 2;
            //int idPres = ir*(nsd_+1) + 3;


            for(int ic =0; ic<emb_nen_; ic++)
            {
              int iVelx = ic*(nsd_+1)+0;
              int iVely = ic*(nsd_+1)+1;
              int iVelz = ic*(nsd_+1)+2;
              //int iPres = ic*(nsd_+1)+3;


              // + (v2, (2*k2*mu2) *eps(Du2)*n)

              //(x,x)
              C_uiui_(idVelx, iVelx) += s_funct_visc2_timefacfac(ir)*(         normal(Velx)*emb_derxy_(Velx,ic)
                                                                       + 0.5 * normal(Vely)*emb_derxy_(Vely,ic)
                                                                       + 0.5 * normal(Velz)*emb_derxy_(Velz,ic)  );
              //(x,y)
              C_uiui_(idVelx, iVely) += s_funct_visc2_timefacfac(ir)*    0.5 * normal(Vely)*emb_derxy_(Velx,ic);
              //(x,z)
              C_uiui_(idVelx, iVelz) += s_funct_visc2_timefacfac(ir)*    0.5 * normal(Velz)*emb_derxy_(Velx,ic);

              //(y,x)
              C_uiui_(idVely, iVelx) += s_funct_visc2_timefacfac(ir)*    0.5 * normal(Velx)*emb_derxy_(Vely,ic);
              //(y,y)
              C_uiui_(idVely, iVely) += s_funct_visc2_timefacfac(ir)*(   0.5 * normal(Velx)*emb_derxy_(Velx,ic)
                                                                       +       normal(Vely)*emb_derxy_(Vely,ic)
                                                                       + 0.5 * normal(Velz)*emb_derxy_(Velz,ic)  );
              //(y,z)
              C_uiui_(idVely, iVelz) += s_funct_visc2_timefacfac(ir)*    0.5 * normal(Velz)*emb_derxy_(Vely,ic);

              //(z,x)
              C_uiui_(idVelz, iVelx) += s_funct_visc2_timefacfac(ir)*    0.5 * normal(Velx)*emb_derxy_(Velz,ic);
              //(z,y)
              C_uiui_(idVelz, iVely) += s_funct_visc2_timefacfac(ir)*    0.5 * normal(Vely)*emb_derxy_(Velz,ic);
              //(z,z)
              C_uiui_(idVelz, iVelz) += s_funct_visc2_timefacfac(ir)*(   0.5 * normal(Velx)*emb_derxy_(Velx,ic)
                                                                       + 0.5 * normal(Vely)*emb_derxy_(Vely,ic)
                                                                       +       normal(Velz)*emb_derxy_(Velz,ic)  );
            }

            // - (v2, (2*k2*mu2) *eps(Du2)*n)
            rhC_ui_(idVelx,0) -= s_funct_visc2_timefacfac(ir)*(            emb_vderxy_(Velx,Velx)                          *normal(Velx)
                                                                 + 0.5 * ( emb_vderxy_(Velx,Vely) + emb_vderxy_(Vely,Velx))*normal(Vely)
                                                                 + 0.5 * ( emb_vderxy_(Velx,Velz) + emb_vderxy_(Velz,Velx))*normal(Velz)  );
            rhC_ui_(idVely,0) -= s_funct_visc2_timefacfac(ir)*(    0.5 * ( emb_vderxy_(Vely,Velx) + emb_vderxy_(Velx,Vely))*normal(Velx)
                                                                 +         emb_vderxy_(Vely,Vely)                          *normal(Vely)
                                                                 + 0.5 * ( emb_vderxy_(Vely,Velz) + emb_vderxy_(Velz,Vely))*normal(Velz)  );
            rhC_ui_(idVelz,0) -= s_funct_visc2_timefacfac(ir)*(    0.5 * ( emb_vderxy_(Velz,Velx) + emb_vderxy_(Velx,Velz))*normal(Velx)
                                                                 + 0.5 * ( emb_vderxy_(Velz,Vely) + emb_vderxy_(Vely,Velz))*normal(Vely)
                                                                 +         emb_vderxy_(Velz,Velz)                          *normal(Velz)  );

          }

        }




          /*                                \       /                             i   \
       - |  alpha* { 2mu*eps(v) }*n , [ Du ] |  =  |  alpha* { 2mu eps(v) }*n ,[ u ]   |
          \                                 /       \                                */
       // antisymmetric formulation (see Burman, Fernandez 2009)
       double alpha = +1.0; // (+1)symmetric, (-1) unsymmetric


       //-----------------------------------------------
       //    - ((2*k1*mu1) *eps(v1)*n , u1)
       //-----------------------------------------------
       for(int ir = 0; ir<nen_; ir++)
       {
       int idVelx = ir*(nsd_+1) + 0;
       int idVely = ir*(nsd_+1) + 1;
       int idVelz = ir*(nsd_+1) + 2;

       // -(2mu*eps(v)*n, Du)
       for(int ic =0; ic<nen_; ic++)
       {
         int iVelx = ic*(nsd_+1)+0;
         int iVely = ic*(nsd_+1)+1;
         int iVelz = ic*(nsd_+1)+2;

         //(x,x)
         C_uu_(idVelx, iVelx) -= alpha*e_funct_visc1_timefacfac(ic)*(         normal(Velx)*derxy_(Velx,ir)
                                                                      + 0.5 * normal(Vely)*derxy_(Vely,ir)
                                                                      + 0.5 * normal(Velz)*derxy_(Velz,ir)  );
         //(y,x)
         C_uu_(idVely, iVelx) -= alpha*e_funct_visc1_timefacfac(ic)*    0.5 * normal(Vely)*derxy_(Velx,ir);
         //(z,x)
         C_uu_(idVelz, iVelx) -= alpha*e_funct_visc1_timefacfac(ic)*    0.5 * normal(Velz)*derxy_(Velx,ir);

         //(x,y)
         C_uu_(idVelx, iVely) -= alpha*e_funct_visc1_timefacfac(ic)*    0.5 * normal(Velx)*derxy_(Vely,ir);
         //(y,y)
         C_uu_(idVely, iVely) -= alpha*e_funct_visc1_timefacfac(ic)*(   0.5 * normal(Velx)*derxy_(Velx,ir)
                                                                      +       normal(Vely)*derxy_(Vely,ir)
                                                                      + 0.5 * normal(Velz)*derxy_(Velz,ir)  );
         //(z,y)
         C_uu_(idVelz, iVely) -= alpha*e_funct_visc1_timefacfac(ic)*    0.5 * normal(Velz)*derxy_(Vely,ir);

         //(x,z)
         C_uu_(idVelx, iVelz) -= alpha*e_funct_visc1_timefacfac(ic)*    0.5 * normal(Velx)*derxy_(Velz,ir);
         //(y,z)
         C_uu_(idVely, iVelz) -= alpha*e_funct_visc1_timefacfac(ic)*    0.5 * normal(Vely)*derxy_(Velz,ir);
         //(z,z)
         C_uu_(idVelz, iVelz) -= alpha*e_funct_visc1_timefacfac(ic)*(   0.5 * normal(Velx)*derxy_(Velx,ir)
                                                                      + 0.5 * normal(Vely)*derxy_(Vely,ir)
                                                                      +       normal(Velz)*derxy_(Velz,ir)  );
       }
       //  (2mu1k1*eps(v1)*n, u1)
       double timefacfac_visc = alpha*timefacfac*2.0*visceff_1*kappa1;
       rhs_Cu_(idVelx,0) += timefacfac_visc* (     derxy_(Velx,ir) *       normal(Velx) * velint(Velx)
                                                 + derxy_(Vely,ir) * 0.5* (normal(Vely) * velint(Velx) + normal(Velx)*velint(Vely))
                                                 + derxy_(Velz,ir) * 0.5* (normal(Velz) * velint(Velx) + normal(Velx)*velint(Velz)));

       rhs_Cu_(idVely,0) += timefacfac_visc* (     derxy_(Velx,ir) * 0.5* (normal(Vely) * velint(Velx) + normal(Velx)*velint(Vely))
                                                 + derxy_(Vely,ir) *       normal(Vely) * velint(Vely)
                                                 + derxy_(Velz,ir) * 0.5* (normal(Velz) * velint(Vely) + normal(Vely)*velint(Velz)));

       rhs_Cu_(idVelz,0) += timefacfac_visc* (     derxy_(Velx,ir) * 0.5* (normal(Velx) * velint(Velz) + normal(Velz)*velint(Velx))
                                                 + derxy_(Vely,ir) * 0.5* (normal(Vely) * velint(Velz) + normal(Velz)*velint(Vely))
                                                 + derxy_(Velz,ir) *       normal(Velz) * velint(Velz));

//       if(!coupling) // weak Dirichlet case
//       {
//         // -(2mu*eps(v)*n, u_DBC)
//         rhs_Cu_(idVelx,0) -= timefacfac_visc* (  derxy_(Velx,ir) *       normal(Velx) * ivelint_WDBC_JUMP(Velx)
//                                                + derxy_(Vely,ir) * 0.5* (normal(Vely) * ivelint_WDBC_JUMP(Velx) + normal(Velx)*ivelint_WDBC_JUMP(Vely))
//                                                + derxy_(Velz,ir) * 0.5* (normal(Velz) * ivelint_WDBC_JUMP(Velx) + normal(Velx)*ivelint_WDBC_JUMP(Velz)));
//
//         rhs_Cu_(idVely,0) -= timefacfac_visc* (  derxy_(Velx,ir) * 0.5* (normal(Vely) * ivelint_WDBC_JUMP(Velx) + normal(Velx)*ivelint_WDBC_JUMP(Vely))
//                                                + derxy_(Vely,ir) *       normal(Vely) * ivelint_WDBC_JUMP(Vely)
//                                                + derxy_(Velz,ir) * 0.5* (normal(Velz) * ivelint_WDBC_JUMP(Vely) + normal(Vely)*ivelint_WDBC_JUMP(Velz)));
//
//         rhs_Cu_(idVelz,0) -= timefacfac_visc* (  derxy_(Velx,ir) * 0.5* (normal(Velx) * ivelint_WDBC_JUMP(Velz) + normal(Velz)*ivelint_WDBC_JUMP(Velx))
//                                                + derxy_(Vely,ir) * 0.5* (normal(Vely) * ivelint_WDBC_JUMP(Velz) + normal(Velz)*ivelint_WDBC_JUMP(Vely))
//                                                + derxy_(Velz,ir) *       normal(Velz) * ivelint_WDBC_JUMP(Velz));
//
//       }
       }

       if(coupling)
       {
       //-----------------------------------------------
       //    + ((2*k1*mu1) *eps(v1)*n , u2)
       //-----------------------------------------------
       for(int ir = 0; ir<nen_; ir++)
       {
       int idVelx = ir*(nsd_+1) + 0;
       int idVely = ir*(nsd_+1) + 1;
       int idVelz = ir*(nsd_+1) + 2;

       // -(2mu*eps(v1)*n, Du2)
       for(int ic =0; ic<emb_nen_; ic++)
       {
         int iVelx = ic*(nsd_+1)+0;
         int iVely = ic*(nsd_+1)+1;
         int iVelz = ic*(nsd_+1)+2;

         //(x,x)
         C_uui_(idVelx, iVelx) += alpha*s_funct_visc1_timefacfac(ic)*(         normal(Velx)*derxy_(Velx,ir)
                                                                       + 0.5 * normal(Vely)*derxy_(Vely,ir)
                                                                       + 0.5 * normal(Velz)*derxy_(Velz,ir)  );
         //(y,x)
         C_uui_(idVely, iVelx) += alpha*s_funct_visc1_timefacfac(ic)*    0.5 * normal(Vely)*derxy_(Velx,ir);
         //(z,x)
         C_uui_(idVelz, iVelx) += alpha*s_funct_visc1_timefacfac(ic)*    0.5 * normal(Velz)*derxy_(Velx,ir);

         //(x,y)
         C_uui_(idVelx, iVely) += alpha*s_funct_visc1_timefacfac(ic)*    0.5 * normal(Velx)*derxy_(Vely,ir);
         //(y,y)
         C_uui_(idVely, iVely) += alpha*s_funct_visc1_timefacfac(ic)*(   0.5 * normal(Velx)*derxy_(Velx,ir)
                                                                       +       normal(Vely)*derxy_(Vely,ir)
                                                                       + 0.5 * normal(Velz)*derxy_(Velz,ir)  );
         //(z,y)
         C_uui_(idVelz, iVely) += alpha*s_funct_visc1_timefacfac(ic)*    0.5 * normal(Velz)*derxy_(Vely,ir);

         //(x,z)
         C_uui_(idVelx, iVelz) += alpha*s_funct_visc1_timefacfac(ic)*    0.5 * normal(Velx)*derxy_(Velz,ir);
         //(y,z)
         C_uui_(idVely, iVelz) += alpha*s_funct_visc1_timefacfac(ic)*    0.5 * normal(Vely)*derxy_(Velz,ir);
         //(z,z)
         C_uui_(idVelz, iVelz) += alpha*s_funct_visc1_timefacfac(ic)*(   0.5 * normal(Velx)*derxy_(Velx,ir)
                                                                       + 0.5 * normal(Vely)*derxy_(Vely,ir)
                                                                       +       normal(Velz)*derxy_(Velz,ir)  );
       }
       //  (2k1mu1*eps(v1)*n, u2)
       double timefacfac_visc = alpha*timefacfac*2.0*visceff_1*kappa1;
       rhs_Cu_(idVelx,0) -= timefacfac_visc* (     derxy_(Velx,ir) *       normal(Velx) * emb_velint(Velx)
                                                 + derxy_(Vely,ir) * 0.5* (normal(Vely) * emb_velint(Velx) + normal(Velx)*emb_velint(Vely))
                                                 + derxy_(Velz,ir) * 0.5* (normal(Velz) * emb_velint(Velx) + normal(Velx)*emb_velint(Velz)));

       rhs_Cu_(idVely,0) -= timefacfac_visc* (     derxy_(Velx,ir) * 0.5* (normal(Vely) * emb_velint(Velx) + normal(Velx)*emb_velint(Vely))
                                                 + derxy_(Vely,ir) *       normal(Vely) * emb_velint(Vely)
                                                 + derxy_(Velz,ir) * 0.5* (normal(Velz) * emb_velint(Vely) + normal(Vely)*emb_velint(Velz)));

       rhs_Cu_(idVelz,0) -= timefacfac_visc* (     derxy_(Velx,ir) * 0.5* (normal(Velx) * emb_velint(Velz) + normal(Velz)*emb_velint(Velx))
                                                 + derxy_(Vely,ir) * 0.5* (normal(Vely) * emb_velint(Velz) + normal(Velz)*emb_velint(Vely))
                                                 + derxy_(Velz,ir) *       normal(Velz) * emb_velint(Velz));

//         // -(2mu*eps(v)*n, u_DBC)
//         elevec1_epetra(idVelx,0) -= timefacfac_visc* (  derxy_(Velx,ir) *       normal(Velx) * velint_WDBC(Velx)
//                                                   + derxy_(Vely,ir) * 0.5* (normal(Vely) * velint_WDBC(Velx) + normal(Velx)*velint_WDBC(Vely))
//                                                   + derxy_(Velz,ir) * 0.5* (normal(Velz) * velint_WDBC(Velx) + normal(Velx)*velint_WDBC(Velz)));
//
//         elevec1_epetra(idVely,0) -= timefacfac_visc* (  derxy_(Velx,ir) * 0.5* (normal(Vely) * velint_WDBC(Velx) + normal(Velx)*velint_WDBC(Vely))
//                                                   + derxy_(Vely,ir) *       normal(Vely) * velint_WDBC(Vely)
//                                                   + derxy_(Velz,ir) * 0.5* (normal(Velz) * velint_WDBC(Vely) + normal(Vely)*velint_WDBC(Velz)));
//
//         elevec1_epetra(idVelz,0) -= timefacfac_visc* (  derxy_(Velx,ir) * 0.5* (normal(Velx) * velint_WDBC(Velz) + normal(Velz)*velint_WDBC(Velx))
//                                                   + derxy_(Vely,ir) * 0.5* (normal(Vely) * velint_WDBC(Velz) + normal(Velz)*velint_WDBC(Vely))
//                                                   + derxy_(Velz,ir) *       normal(Velz) * velint_WDBC(Velz));
       }
       }// end coupling

       if(!bg_mortaring)
       {
         //-----------------------------------------------
         //    - ((2*k2*mu2) *eps(v2)*n , u1)
         //-----------------------------------------------
         for(int ir = 0; ir<emb_nen_; ir++)
         {
         int idVelx = ir*(nsd_+1) + 0;
         int idVely = ir*(nsd_+1) + 1;
         int idVelz = ir*(nsd_+1) + 2;

         // -(2mu*eps(v2)*n, Du1)
         for(int ic =0; ic<nen_; ic++)
         {
           int iVelx = ic*(nsd_+1)+0;
           int iVely = ic*(nsd_+1)+1;
           int iVelz = ic*(nsd_+1)+2;

           //(x,x)
           C_uiu_(idVelx, iVelx) -= alpha*e_funct_visc2_timefacfac(ic)*(         normal(Velx)*emb_derxy_(Velx,ir)
                                                                         + 0.5 * normal(Vely)*emb_derxy_(Vely,ir)
                                                                         + 0.5 * normal(Velz)*emb_derxy_(Velz,ir)  );
           //(y,x)
           C_uiu_(idVely, iVelx) -= alpha*e_funct_visc2_timefacfac(ic)*    0.5 * normal(Vely)*emb_derxy_(Velx,ir);
           //(z,x)
           C_uiu_(idVelz, iVelx) -= alpha*e_funct_visc2_timefacfac(ic)*    0.5 * normal(Velz)*emb_derxy_(Velx,ir);

           //(x,y)
           C_uiu_(idVelx, iVely) -= alpha*e_funct_visc2_timefacfac(ic)*    0.5 * normal(Velx)*emb_derxy_(Vely,ir);
           //(y,y)
           C_uiu_(idVely, iVely) -= alpha*e_funct_visc2_timefacfac(ic)*(   0.5 * normal(Velx)*emb_derxy_(Velx,ir)
                                                                         +       normal(Vely)*emb_derxy_(Vely,ir)
                                                                         + 0.5 * normal(Velz)*emb_derxy_(Velz,ir)  );
           //(z,y)
           C_uiu_(idVelz, iVely) -= alpha*e_funct_visc2_timefacfac(ic)*    0.5 * normal(Velz)*emb_derxy_(Vely,ir);

           //(x,z)
           C_uiu_(idVelx, iVelz) -= alpha*e_funct_visc2_timefacfac(ic)*    0.5 * normal(Velx)*emb_derxy_(Velz,ir);
           //(y,z)
           C_uiu_(idVely, iVelz) -= alpha*e_funct_visc2_timefacfac(ic)*    0.5 * normal(Vely)*emb_derxy_(Velz,ir);
           //(z,z)
           C_uiu_(idVelz, iVelz) -= alpha*e_funct_visc2_timefacfac(ic)*(   0.5 * normal(Velx)*emb_derxy_(Velx,ir)
                                                                         + 0.5 * normal(Vely)*emb_derxy_(Vely,ir)
                                                                         +       normal(Velz)*emb_derxy_(Velz,ir)  );
         }
         //  (2k2mu2*eps(v2)*n, u1)
         double timefacfac_visc = alpha*timefacfac*2.0*visceff_2*kappa2;
         rhC_ui_(idVelx,0) += timefacfac_visc* (     emb_derxy_(Velx,ir) *       normal(Velx) * velint(Velx)
                                                   + emb_derxy_(Vely,ir) * 0.5* (normal(Vely) * velint(Velx) + normal(Velx)*velint(Vely))
                                                   + emb_derxy_(Velz,ir) * 0.5* (normal(Velz) * velint(Velx) + normal(Velx)*velint(Velz)));

         rhC_ui_(idVely,0) += timefacfac_visc* (     emb_derxy_(Velx,ir) * 0.5* (normal(Vely) * velint(Velx) + normal(Velx)*velint(Vely))
                                                   + emb_derxy_(Vely,ir) *       normal(Vely) * velint(Vely)
                                                   + emb_derxy_(Velz,ir) * 0.5* (normal(Velz) * velint(Vely) + normal(Vely)*velint(Velz)));

         rhC_ui_(idVelz,0) += timefacfac_visc* (     emb_derxy_(Velx,ir) * 0.5* (normal(Velx) * velint(Velz) + normal(Velz)*velint(Velx))
                                                   + emb_derxy_(Vely,ir) * 0.5* (normal(Vely) * velint(Velz) + normal(Velz)*velint(Vely))
                                                   + emb_derxy_(Velz,ir) *       normal(Velz) * velint(Velz));

         }

         //-----------------------------------------------
         //    + ((2*k2*mu2) *eps(v2)*n , u2)
         //-----------------------------------------------
         for(int ir = 0; ir<emb_nen_; ir++)
         {
         int idVelx = ir*(nsd_+1) + 0;
         int idVely = ir*(nsd_+1) + 1;
         int idVelz = ir*(nsd_+1) + 2;

         // -(2mu2*eps(v2)*n, Du2)
         for(int ic =0; ic<emb_nen_; ic++)
         {
           int iVelx = ic*(nsd_+1)+0;
           int iVely = ic*(nsd_+1)+1;
           int iVelz = ic*(nsd_+1)+2;

           //(x,x)
           C_uiui_(idVelx, iVelx) += alpha*s_funct_visc2_timefacfac(ic)*(         normal(Velx)*emb_derxy_(Velx,ir)
                                                                          + 0.5 * normal(Vely)*emb_derxy_(Vely,ir)
                                                                          + 0.5 * normal(Velz)*emb_derxy_(Velz,ir)  );
           //(y,x)
           C_uiui_(idVely, iVelx) += alpha*s_funct_visc2_timefacfac(ic)*    0.5 * normal(Vely)*emb_derxy_(Velx,ir);
           //(z,x)
           C_uiui_(idVelz, iVelx) += alpha*s_funct_visc2_timefacfac(ic)*    0.5 * normal(Velz)*emb_derxy_(Velx,ir);

           //(x,y)
           C_uiui_(idVelx, iVely) += alpha*s_funct_visc2_timefacfac(ic)*    0.5 * normal(Velx)*emb_derxy_(Vely,ir);
           //(y,y)
           C_uiui_(idVely, iVely) += alpha*s_funct_visc2_timefacfac(ic)*(   0.5 * normal(Velx)*emb_derxy_(Velx,ir)
                                                                          +       normal(Vely)*emb_derxy_(Vely,ir)
                                                                          + 0.5 * normal(Velz)*emb_derxy_(Velz,ir)  );
           //(z,y)
           C_uiui_(idVelz, iVely) += alpha*s_funct_visc2_timefacfac(ic)*    0.5 * normal(Velz)*emb_derxy_(Vely,ir);

           //(x,z)
           C_uiui_(idVelx, iVelz) += alpha*s_funct_visc2_timefacfac(ic)*    0.5 * normal(Velx)*emb_derxy_(Velz,ir);
           //(y,z)
           C_uiui_(idVely, iVelz) += alpha*s_funct_visc2_timefacfac(ic)*    0.5 * normal(Vely)*emb_derxy_(Velz,ir);
           //(z,z)
           C_uiui_(idVelz, iVelz) += alpha*s_funct_visc2_timefacfac(ic)*(   0.5 * normal(Velx)*emb_derxy_(Velx,ir)
                                                                          + 0.5 * normal(Vely)*emb_derxy_(Vely,ir)
                                                                          +       normal(Velz)*emb_derxy_(Velz,ir)  );
         }
         //  (2mu2*eps(v2)*n, u2)
         double timefacfac_visc = alpha*timefacfac*2.0*visceff_2*kappa2;
         rhC_ui_(idVelx,0) -= timefacfac_visc* (     emb_derxy_(Velx,ir) *       normal(Velx) * emb_velint(Velx)
                                                   + emb_derxy_(Vely,ir) * 0.5* (normal(Vely) * emb_velint(Velx) + normal(Velx)*emb_velint(Vely))
                                                   + emb_derxy_(Velz,ir) * 0.5* (normal(Velz) * emb_velint(Velx) + normal(Velx)*emb_velint(Velz)));

         rhC_ui_(idVely,0) -= timefacfac_visc* (     emb_derxy_(Velx,ir) * 0.5* (normal(Vely) * emb_velint(Velx) + normal(Velx)*emb_velint(Vely))
                                                   + emb_derxy_(Vely,ir) *       normal(Vely) * emb_velint(Vely)
                                                   + emb_derxy_(Velz,ir) * 0.5* (normal(Velz) * emb_velint(Vely) + normal(Vely)*emb_velint(Velz)));

         rhC_ui_(idVelz,0) -= timefacfac_visc* (     emb_derxy_(Velx,ir) * 0.5* (normal(Velx) * emb_velint(Velz) + normal(Velz)*emb_velint(Velx))
                                                   + emb_derxy_(Vely,ir) * 0.5* (normal(Vely) * emb_velint(Velz) + normal(Velz)*emb_velint(Vely))
                                                   + emb_derxy_(Velz,ir) *       normal(Velz) * emb_velint(Velz));

        }

       }




          /*                                  \        /                           i   \
         |  gamma*mu/h_K *  [ v ] , [ Du ]     | =  - |   gamma*mu/h_K * [ v ], [ u ]   |
          \                                   /        \                              */


       // + gamma*mu/h_K (v1, u1))
       for(int ir=0; ir<nen_; ir++)
       {
         int idVelx = ir*(nsd_+1) + 0;
         int idVely = ir*(nsd_+1) + 1;
         int idVelz = ir*(nsd_+1) + 2;

         for(int ic=0; ic<nen_; ic++)
         {
           int iVelx = ic*(nsd_+1)+0;
           int iVely = ic*(nsd_+1)+1;
           int iVelz = ic*(nsd_+1)+2;

           C_uu_(idVelx, iVelx) += funct_dyad_timefacfac(ir,ic)*stabfac;
           C_uu_(idVely, iVely) += funct_dyad_timefacfac(ir,ic)*stabfac;
           C_uu_(idVelz, iVelz) += funct_dyad_timefacfac(ir,ic)*stabfac;
         }

         // -(stab * v, u)
         rhs_Cu_(idVelx,0) -= funct_timefacfac(ir)*stabfac*velint(Velx);
         rhs_Cu_(idVely,0) -= funct_timefacfac(ir)*stabfac*velint(Vely);
         rhs_Cu_(idVelz,0) -= funct_timefacfac(ir)*stabfac*velint(Velz);

         if(!coupling) // weak Dirichlet case
         {
           // +(stab * v, u_DBC)
           rhs_Cu_(idVelx,0) += funct_timefacfac(ir)*stabfac*ivelint_WDBC_JUMP(Velx);
           rhs_Cu_(idVely,0) += funct_timefacfac(ir)*stabfac*ivelint_WDBC_JUMP(Vely);
           rhs_Cu_(idVelz,0) += funct_timefacfac(ir)*stabfac*ivelint_WDBC_JUMP(Velz);
         }

       }

       if(coupling)
       {
       // - gamma*mu/h_K (v1, u2))
       for(int ir=0; ir<nen_; ir++)
       {
         int idVelx = ir*(nsd_+1) + 0;
         int idVely = ir*(nsd_+1) + 1;
         int idVelz = ir*(nsd_+1) + 2;

         for(int ic=0; ic<emb_nen_; ic++)
         {
           int iVelx = ic*(nsd_+1)+0;
           int iVely = ic*(nsd_+1)+1;
           int iVelz = ic*(nsd_+1)+2;

           C_uui_(idVelx, iVelx) -= emb_funct_dyad_timefacfac(ic,ir)*stabfac;
           C_uui_(idVely, iVely) -= emb_funct_dyad_timefacfac(ic,ir)*stabfac;
           C_uui_(idVelz, iVelz) -= emb_funct_dyad_timefacfac(ic,ir)*stabfac;
         }

         // -(stab * v, u)
         rhs_Cu_(idVelx,0) += funct_timefacfac(ir)*stabfac*emb_velint(Velx);
         rhs_Cu_(idVely,0) += funct_timefacfac(ir)*stabfac*emb_velint(Vely);
         rhs_Cu_(idVelz,0) += funct_timefacfac(ir)*stabfac*emb_velint(Velz);


       }


       // - gamma*mu/h_K (v2, u1))
       for(int ir=0; ir<emb_nen_; ir++)
       {
         int idVelx = ir*(nsd_+1) + 0;
         int idVely = ir*(nsd_+1) + 1;
         int idVelz = ir*(nsd_+1) + 2;

         for(int ic=0; ic<nen_; ic++)
         {
           int iVelx = ic*(nsd_+1)+0;
           int iVely = ic*(nsd_+1)+1;
           int iVelz = ic*(nsd_+1)+2;

           C_uiu_(idVelx, iVelx) -= emb_funct_dyad_timefacfac(ir,ic)*stabfac;
           C_uiu_(idVely, iVely) -= emb_funct_dyad_timefacfac(ir,ic)*stabfac;
           C_uiu_(idVelz, iVelz) -= emb_funct_dyad_timefacfac(ir,ic)*stabfac;
         }

         // +(stab * v2, u1)
         rhC_ui_(idVelx,0) += emb_funct_timefacfac(ir)*stabfac*velint(Velx);
         rhC_ui_(idVely,0) += emb_funct_timefacfac(ir)*stabfac*velint(Vely);
         rhC_ui_(idVelz,0) += emb_funct_timefacfac(ir)*stabfac*velint(Velz);

       }

//       LINALG::Matrix<emb_nen_,emb_nen_> emb_emb_dyad_timefacfac(true);
//       emb_emb_dyad_timefacfac.MultiplyNT(emb_funct_timefacfac, emb_funct_);
       // + gamma*mu/h_K (v2, u2))
       for(int ir=0; ir<emb_nen_; ir++)
       {
         int idVelx = ir*(nsd_+1) + 0;
         int idVely = ir*(nsd_+1) + 1;
         int idVelz = ir*(nsd_+1) + 2;

         for(int ic=0; ic<emb_nen_; ic++)
         {
           int iVelx = ic*(nsd_+1)+0;
           int iVely = ic*(nsd_+1)+1;
           int iVelz = ic*(nsd_+1)+2;

           C_uiui_(idVelx, iVelx) += emb_emb_dyad_timefacfac(ir,ic)*stabfac;
           C_uiui_(idVely, iVely) += emb_emb_dyad_timefacfac(ir,ic)*stabfac;
           C_uiui_(idVelz, iVelz) += emb_emb_dyad_timefacfac(ir,ic)*stabfac;
         }

         // -(stab * v2, u2)
         rhC_ui_(idVelx,0) -= emb_funct_timefacfac(ir)*stabfac*emb_velint(Velx);
         rhC_ui_(idVely,0) -= emb_funct_timefacfac(ir)*stabfac*emb_velint(Vely);
         rhC_ui_(idVelz,0) -= emb_funct_timefacfac(ir)*stabfac*emb_velint(Velz);

       }
       } // end coupling



       if(fabs(stabfac_conv) > 0.0)
       {

         // convective stabilization
         /*                           \        /                       i   _     \
        |  gamma/h_K *  v*n , Du*n     | =  - |   gamma/h_K *  v*n , (u  - u)*n   |
         \                            /        \                                */


         for(int ir=0; ir<nen_; ir++)
         {
           int idVelx = ir*(nsd_+1) + 0;
           int idVely = ir*(nsd_+1) + 1;
           int idVelz = ir*(nsd_+1) + 2;

           // (stab * v1, Du1)
           for(int ic=0; ic<nen_; ic++)
           {
             int iVelx = ic*(nsd_+1)+0;
             int iVely = ic*(nsd_+1)+1;
             int iVelz = ic*(nsd_+1)+2;

             C_uu_(idVelx, iVelx) += funct_dyad_timefacfac(ir,ic)*stabfac_conv*normal(Velx)*normal(Velx);
             C_uu_(idVelx, iVely) += funct_dyad_timefacfac(ir,ic)*stabfac_conv*normal(Velx)*normal(Vely);
             C_uu_(idVelx, iVelz) += funct_dyad_timefacfac(ir,ic)*stabfac_conv*normal(Velx)*normal(Velz);

             C_uu_(idVely, iVelx) += funct_dyad_timefacfac(ir,ic)*stabfac_conv*normal(Vely)*normal(Velx);
             C_uu_(idVely, iVely) += funct_dyad_timefacfac(ir,ic)*stabfac_conv*normal(Vely)*normal(Vely);
             C_uu_(idVely, iVelz) += funct_dyad_timefacfac(ir,ic)*stabfac_conv*normal(Vely)*normal(Velz);

             C_uu_(idVelz, iVelx) += funct_dyad_timefacfac(ir,ic)*stabfac_conv*normal(Velz)*normal(Velx);
             C_uu_(idVelz, iVely) += funct_dyad_timefacfac(ir,ic)*stabfac_conv*normal(Velz)*normal(Vely);
             C_uu_(idVelz, iVelz) += funct_dyad_timefacfac(ir,ic)*stabfac_conv*normal(Velz)*normal(Velz);
           }

           double velint_normal = velint.Dot(normal);

           // -(stab * v1*n, u1*n)
           rhs_Cu_(idVelx) -= funct_timefacfac(ir)*normal(Velx)*stabfac_conv*velint_normal;
           rhs_Cu_(idVely) -= funct_timefacfac(ir)*normal(Vely)*stabfac_conv*velint_normal;
           rhs_Cu_(idVelz) -= funct_timefacfac(ir)*normal(Velz)*stabfac_conv*velint_normal;

           if (!coupling)
           {
             double velint_WDBC_normal = ivelint_WDBC_JUMP.Dot(normal);
             // +(stab * v1*n, u_DBC*n)
             rhs_Cu_(idVelx) += funct_timefacfac(ir)*normal(Velx)*stabfac_conv*velint_WDBC_normal;
             rhs_Cu_(idVely) += funct_timefacfac(ir)*normal(Vely)*stabfac_conv*velint_WDBC_normal;
             rhs_Cu_(idVelz) += funct_timefacfac(ir)*normal(Velz)*stabfac_conv*velint_WDBC_normal;
           }
         }

         if(coupling)
         {

           for(int ir=0; ir<nen_; ir++)
           {
             int idVelx = ir*(nsd_+1) + 0;
             int idVely = ir*(nsd_+1) + 1;
             int idVelz = ir*(nsd_+1) + 2;

             // (stab * v1, Du2)
             for(int ic=0; ic<emb_nen_; ic++)
             {
               int iVelx = ic*(nsd_+1)+0;
               int iVely = ic*(nsd_+1)+1;
               int iVelz = ic*(nsd_+1)+2;

               C_uui_(idVelx, iVelx) -= emb_funct_dyad_timefacfac(ic,ir)*stabfac_conv*normal(Velx)*normal(Velx);
               C_uui_(idVelx, iVely) -= emb_funct_dyad_timefacfac(ic,ir)*stabfac_conv*normal(Velx)*normal(Vely);
               C_uui_(idVelx, iVelz) -= emb_funct_dyad_timefacfac(ic,ir)*stabfac_conv*normal(Velx)*normal(Velz);

               C_uui_(idVely, iVelx) -= emb_funct_dyad_timefacfac(ic,ir)*stabfac_conv*normal(Vely)*normal(Velx);
               C_uui_(idVely, iVely) -= emb_funct_dyad_timefacfac(ic,ir)*stabfac_conv*normal(Vely)*normal(Vely);
               C_uui_(idVely, iVelz) -= emb_funct_dyad_timefacfac(ic,ir)*stabfac_conv*normal(Vely)*normal(Velz);

               C_uui_(idVelz, iVelx) -= emb_funct_dyad_timefacfac(ic,ir)*stabfac_conv*normal(Velz)*normal(Velx);
               C_uui_(idVelz, iVely) -= emb_funct_dyad_timefacfac(ic,ir)*stabfac_conv*normal(Velz)*normal(Vely);
               C_uui_(idVelz, iVelz) -= emb_funct_dyad_timefacfac(ic,ir)*stabfac_conv*normal(Velz)*normal(Velz);
             }

             double emb_velint_normal = emb_velint.Dot(normal);

             // -(stab * v1*n, u2*n)
             rhs_Cu_(idVelx) += funct_timefacfac(ir)*normal(Velx)*stabfac_conv*emb_velint_normal;
             rhs_Cu_(idVely) += funct_timefacfac(ir)*normal(Vely)*stabfac_conv*emb_velint_normal;
             rhs_Cu_(idVelz) += funct_timefacfac(ir)*normal(Velz)*stabfac_conv*emb_velint_normal;

           }

           for(int ir=0; ir<emb_nen_; ir++)
           {
             int idVelx = ir*(nsd_+1) + 0;
             int idVely = ir*(nsd_+1) + 1;
             int idVelz = ir*(nsd_+1) + 2;

             // (stab * v2, Du1)
             for(int ic=0; ic<nen_; ic++)
             {
               int iVelx = ic*(nsd_+1)+0;
               int iVely = ic*(nsd_+1)+1;
               int iVelz = ic*(nsd_+1)+2;

               C_uiu_(idVelx, iVelx) -= emb_funct_dyad_timefacfac(ir,ic)*stabfac_conv*normal(Velx)*normal(Velx);
               C_uiu_(idVelx, iVely) -= emb_funct_dyad_timefacfac(ir,ic)*stabfac_conv*normal(Velx)*normal(Vely);
               C_uiu_(idVelx, iVelz) -= emb_funct_dyad_timefacfac(ir,ic)*stabfac_conv*normal(Velx)*normal(Velz);

               C_uiu_(idVely, iVelx) -= emb_funct_dyad_timefacfac(ir,ic)*stabfac_conv*normal(Vely)*normal(Velx);
               C_uiu_(idVely, iVely) -= emb_funct_dyad_timefacfac(ir,ic)*stabfac_conv*normal(Vely)*normal(Vely);
               C_uiu_(idVely, iVelz) -= emb_funct_dyad_timefacfac(ir,ic)*stabfac_conv*normal(Vely)*normal(Velz);

               C_uiu_(idVelz, iVelx) -= emb_funct_dyad_timefacfac(ir,ic)*stabfac_conv*normal(Velz)*normal(Velx);
               C_uiu_(idVelz, iVely) -= emb_funct_dyad_timefacfac(ir,ic)*stabfac_conv*normal(Velz)*normal(Vely);
               C_uiu_(idVelz, iVelz) -= emb_funct_dyad_timefacfac(ir,ic)*stabfac_conv*normal(Velz)*normal(Velz);
             }

             double velint_normal = velint.Dot(normal);

             // -(stab * v2*n, u1*n)
             rhC_ui_(idVelx) += emb_funct_timefacfac(ir)*normal(Velx)*stabfac_conv*velint_normal;
             rhC_ui_(idVely) += emb_funct_timefacfac(ir)*normal(Vely)*stabfac_conv*velint_normal;
             rhC_ui_(idVelz) += emb_funct_timefacfac(ir)*normal(Velz)*stabfac_conv*velint_normal;

           }

         }

         for(int ir=0; ir<emb_nen_; ir++)
         {
           int idVelx = ir*(nsd_+1) + 0;
           int idVely = ir*(nsd_+1) + 1;
           int idVelz = ir*(nsd_+1) + 2;

           // (stab * v2, Du2)
           for(int ic=0; ic<emb_nen_; ic++)
           {
             int iVelx = ic*(nsd_+1)+0;
             int iVely = ic*(nsd_+1)+1;
             int iVelz = ic*(nsd_+1)+2;

             C_uiui_(idVelx, iVelx) += emb_emb_dyad_timefacfac(ir,ic)*stabfac_conv*normal(Velx)*normal(Velx);
             C_uiui_(idVelx, iVely) += emb_emb_dyad_timefacfac(ir,ic)*stabfac_conv*normal(Velx)*normal(Vely);
             C_uiui_(idVelx, iVelz) += emb_emb_dyad_timefacfac(ir,ic)*stabfac_conv*normal(Velx)*normal(Velz);

             C_uiui_(idVely, iVelx) += emb_emb_dyad_timefacfac(ir,ic)*stabfac_conv*normal(Vely)*normal(Velx);
             C_uiui_(idVely, iVely) += emb_emb_dyad_timefacfac(ir,ic)*stabfac_conv*normal(Vely)*normal(Vely);
             C_uiui_(idVely, iVelz) += emb_emb_dyad_timefacfac(ir,ic)*stabfac_conv*normal(Vely)*normal(Velz);

             C_uiui_(idVelz, iVelx) += emb_emb_dyad_timefacfac(ir,ic)*stabfac_conv*normal(Velz)*normal(Velx);
             C_uiui_(idVelz, iVely) += emb_emb_dyad_timefacfac(ir,ic)*stabfac_conv*normal(Velz)*normal(Vely);
             C_uiui_(idVelz, iVelz) += emb_emb_dyad_timefacfac(ir,ic)*stabfac_conv*normal(Velz)*normal(Velz);
           }

           double emb_velint_normal = emb_velint.Dot(normal);

           // -(stab * v1*n, u1*n)
           rhC_ui_(idVelx) -= emb_funct_timefacfac(ir)*normal(Velx)*stabfac_conv*emb_velint_normal;
           rhC_ui_(idVely) -= emb_funct_timefacfac(ir)*normal(Vely)*stabfac_conv*emb_velint_normal;
           rhC_ui_(idVelz) -= emb_funct_timefacfac(ir)*normal(Velz)*stabfac_conv*emb_velint_normal;


         }



       }


      }


      // set prescribed WDBC at Gaussian point
      virtual void get_vel_WeakDBC( LINALG::Matrix<3,1> & emb_velint )
      {
          emb_velint.Clear();
          emb_velint.Multiply(emb_vel_,emb_funct_);
      }


      LINALG::Matrix<4*emb_nen_ ,4*nen_>      C_uiu_;    // row: sidenode1:u1,u2,u3,p, sidenode2:u1,u2,u3,p ... | col: elenode1:u1,u2,u3,p, elenode2:u1,u2,u3,p, ...
      LINALG::Matrix<4*nen_,4*emb_nen_>      C_uui_;    // includes (ui,pi)
      LINALG::Matrix<4*emb_nen_,1>           rhC_ui_;   // includes (ui,pi)
      LINALG::Matrix<4*emb_nen_,4*emb_nen_>  C_uiui_;   // includes (ui,pi) only for Nitsche coupling
      LINALG::Matrix<emb_nen_,1>             emb_funct_;
      LINALG::Matrix<nsd_,emb_nen_>          emb_deriv_;
      LINALG::Matrix<nsd_,nsd_>              emb_vderxy_;
      LINALG::Matrix<nsd_,emb_nen_>          emb_derxy_;
      LINALG::Matrix<3,emb_nen_>             emb_vel_;
      LINALG::Matrix<emb_nen_,1>             emb_pres_;
      LINALG::Matrix<3,emb_nen_>             emb_disp_;
      LINALG::Matrix<3,emb_nen_>             emb_xyze_;
    };




    // for Nitsche coupling
    template<DRT::Element::DiscretizationType distype>
    Teuchos::RCP<EmbCoupling<distype> > EmbCoupling<distype>::TwoSidedImpl(DRT::Element * emb_ele,
                                                                       Epetra_SerialDenseMatrix &  C_uiu,
                                                                       Epetra_SerialDenseMatrix &  C_uui,
                                                                       Epetra_SerialDenseMatrix &  rhC_ui,
                                                                       Epetra_SerialDenseMatrix &  C_uiui,
                                                                       Epetra_SerialDenseMatrix &  emb_xyze
                                                                       )
    {

      EmbCoupling * emb = NULL;
      switch ( emb_ele->Shape() )
      {
       case DRT::Element::tet4:
       {
         typedef EmbImpl<distype,DRT::Element::tet4> EmbImplType;
         emb = new EmbImplType(emb_ele,C_uiu,C_uui,rhC_ui,C_uiui,emb_xyze);
         break;
       }
       case DRT::Element::tet10:
       {
         typedef EmbImpl<distype,DRT::Element::tet10> EmbImplType;
         emb = new EmbImplType(emb_ele,C_uiu,C_uui,rhC_ui,C_uiui,emb_xyze);
         break;
       }
       case DRT::Element::hex8:
       {
         typedef EmbImpl<distype,DRT::Element::hex8> EmbImplType;
         emb = new EmbImplType(emb_ele,C_uiu,C_uui,rhC_ui,C_uiui,emb_xyze);
         break;
       }
       case DRT::Element::hex20:
       {
         typedef EmbImpl<distype,DRT::Element::hex20> EmbImplType;
         emb = new EmbImplType(emb_ele,C_uiu,C_uui,rhC_ui,C_uiui,emb_xyze);
         break;
       }
       case DRT::Element::hex27:
       {
         typedef EmbImpl<distype,DRT::Element::hex27> EmbImplType;
         emb = new EmbImplType(emb_ele,C_uiu,C_uui,rhC_ui,C_uiui,emb_xyze);
         break;
       }
      default:
        dserror( "unsupported side shape %d", emb_ele->Shape() );
      }
      return Teuchos::rcp(emb);
    }

    } // end namespace Xfluid


//  }
//}



/*--------------------------------------------------------------------------------
 *--------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidEleCalcXFEM<distype>::ElementXfemInterface(
  DRT::ELEMENTS::Fluid3 * ele,
  DRT::Discretization & dis,
  const std::vector<int> & lm,
  const DRT::UTILS::GaussIntegration & intpoints,
  DRT::Discretization & cutdis,
  const std::map<int, std::vector<GEO::CUT::BoundaryCell*> > & bcells,
  const std::map<int, std::vector<DRT::UTILS::GaussIntegration> > & bintpoints,
  std::map<int, std::vector<Epetra_SerialDenseMatrix> > & side_coupling,
  Teuchos::ParameterList&    params,
  Epetra_SerialDenseMatrix&  elemat1_epetra,
  Epetra_SerialDenseVector&  elevec1_epetra,
  Epetra_SerialDenseMatrix&  Cuiui
  )
{
  const Teuchos::RCP<Epetra_Vector> iforcecol = params.get<Teuchos::RCP<Epetra_Vector> >("iforcenp", Teuchos::null);


  //----------------------------------------------------------------------------
  //                         ELEMENT GEOMETRY
  //----------------------------------------------------------------------------

  // get node coordinates
  GEO::fillInitialPositionArray< distype, my::nsd_, LINALG::Matrix<my::nsd_,my::nen_> >( ele, my::xyze_ );

  LINALG::Matrix<my::nsd_,my::nen_> evelaf(true);
  LINALG::Matrix<my::nen_,1> epreaf(true);
  my::ExtractValuesFromGlobalVector(dis, lm, *my::rotsymmpbc_, &evelaf, &epreaf, "velaf");

  int eid = ele->Id();
  LINALG::Matrix<my::nen_,my::nen_> bK_ss;
  LINALG::Matrix<my::nen_,my::nen_> invbK_ss( true );
  LINALG::Matrix<my::nen_,my::nen_> half_invbK_ss;
  LINALG::Matrix<my::nen_,my::nen_> conv_x;
  LINALG::Matrix<my::nen_,my::nen_> conv_y;
  LINALG::Matrix<my::nen_,my::nen_> conv_z;

  LINALG::Matrix<my::nen_,1> dx;
  LINALG::Matrix<my::nen_,1> dy;
  LINALG::Matrix<my::nen_,1> dz;

  // get viscosity
  // check here, if we really have a fluid !!
//   Teuchos::RCP<const MAT::Material> material = ele->Material();
//   dsassert(material->MaterialType() == INPAR::MAT::m_fluid, "Material law is not of type m_fluid.");
//   const MAT::NewtonianFluid* actmat = dynamic_cast<const MAT::NewtonianFluid*>(material.get());
//   const double dens = actmat->Density();
//   // dynamic viscosity \mu
//   const double dynvisc = actmat->Viscosity() * dens;

//   const double viscfac = 1.0/(2.0*dynvisc);

  LINALG::BlockMatrix<LINALG::Matrix<my::nen_,my::nen_>,6,(my::nsd_+1)> K_su;
  LINALG::BlockMatrix<LINALG::Matrix<my::nen_,my::nen_>,(my::nsd_+1),6> K_us;
  LINALG::BlockMatrix<LINALG::Matrix<my::nen_,my::nen_>,6,6>        invK_ss;
  LINALG::BlockMatrix<LINALG::Matrix<my::nen_,   1>,6,1>        rhs;


  const unsigned Velx = 0;
  const unsigned Vely = 1;
  const unsigned Velz = 2;
  const unsigned Pres = 3;

  // unused variable
  //const unsigned Velxi = 0;
  //const unsigned Velyi = 1;
  //const unsigned Velzi = 2;

  const unsigned Sigmaxx = 0;
  const unsigned Sigmaxy = 1;
  const unsigned Sigmaxz = 2;
  const unsigned Sigmayx = 1;
  const unsigned Sigmayy = 3;
  const unsigned Sigmayz = 4;
  const unsigned Sigmazx = 2;
  const unsigned Sigmazy = 4;
  const unsigned Sigmazz = 5;

  // volume integral

  for ( DRT::UTILS::GaussIntegration::iterator iquad=intpoints.begin(); iquad!=intpoints.end(); ++iquad )
  {
    // evaluate shape functions and derivatives at integration point
    my::EvalShapeFuncAndDerivsAtIntPoint(iquad,eid);

    //----------------------------------------------------------------------
    // set time-integration factors for left- and right-hand side
    // (two right-hand-side factors: general and for residuals)
    //----------------------------------------------------------------------
    const double timefacfac = my::f3Parameter_->timefac_ * my::fac_;
#if 0
    // rhsresfac does not exist anymore
    double rhsfac           = timefacfac;
    double rhsresfac        = fac_;
    // modify integration factors for right-hand side such that they
    // are identical in case of generalized-alpha time integration:
    if (f3Parameter_->is_genalpha_)
    {
      rhsfac   /= f3Parameter_->alphaF_;
      rhsresfac = rhsfac;
    }
    else
    {
      // modify residual integration factor for right-hand side in instat. case:
      if (not f3Parameter_->is_stationary_) rhsresfac *= f3Parameter_->dt_;
    }
#endif

    //const double visceff_timefacfac = visceff_*timefacfac;

    //const double viscfac = 1.0/(2.0*dynvisc);
    const double viscfac = 1.0/(2.0*my::visceff_);

    // get velocity at integration point
    // (values at n+alpha_F for generalized-alpha scheme, n+1 otherwise)
    my::velint_.Multiply(evelaf,my::funct_);

    // get velocity derivatives at integration point
    // (values at n+alpha_F for generalized-alpha scheme, n+1 otherwise)
    my::vderxy_.MultiplyNT(evelaf,my::derxy_);

    // get pressure at integration point
    // (value at n+alpha_F for generalized-alpha scheme, n+1 otherwise)
    double press = my::funct_.Dot(epreaf);


    for ( int i=0; i<my::nen_; ++i )
    {
      dx( i ) = my::derxy_( 0, i );
      dy( i ) = my::derxy_( 1, i );
      dz( i ) = my::derxy_( 2, i );
    }

    // block - K_ss
    bK_ss.MultiplyNT( my::funct_, my::funct_ );


    conv_x.MultiplyNT( my::funct_, dx );
    conv_y.MultiplyNT( my::funct_, dy );
    conv_z.MultiplyNT( my::funct_, dz );

                       /*                     \
                    - |  virt tau , eps(Dtau)  |
                       \                     */

    invbK_ss.Update( -viscfac*timefacfac, bK_ss, 1.0 );

                   /*                 \
                  | virt tau , eps(Du) |
                   \                 */

    // K_su

    K_su( Sigmaxx, Velx )->Update( timefacfac, conv_x, 1.0 );
    K_su( Sigmaxy, Velx )->Update( timefacfac, conv_y, 1.0 );
    K_su( Sigmayx, Vely )->Update( timefacfac, conv_x, 1.0 );
    K_su( Sigmaxz, Velx )->Update( timefacfac, conv_z, 1.0 );
    K_su( Sigmazx, Velz )->Update( timefacfac, conv_x, 1.0 );
    K_su( Sigmayy, Vely )->Update( timefacfac, conv_y, 1.0 );
    K_su( Sigmayz, Vely )->Update( timefacfac, conv_z, 1.0 );
    K_su( Sigmazy, Velz )->Update( timefacfac, conv_y, 1.0 );
    K_su( Sigmazz, Velz )->Update( timefacfac, conv_z, 1.0 );

    // r_su

    rhs( Sigmaxx, 0 )->Update( - timefacfac* my::vderxy_(0, 0)                 , my::funct_, 1.0 );
    rhs( Sigmaxy, 0 )->Update( - timefacfac*(my::vderxy_(0, 1) + my::vderxy_(1, 0)), my::funct_, 1.0 );
    rhs( Sigmaxz, 0 )->Update( - timefacfac*(my::vderxy_(0, 2) + my::vderxy_(2, 0)), my::funct_, 1.0 );
    rhs( Sigmayy, 0 )->Update( - timefacfac* my::vderxy_(1, 1)                 , my::funct_, 1.0 );
    rhs( Sigmayz, 0 )->Update( - timefacfac*(my::vderxy_(1, 2) + my::vderxy_(2, 1)), my::funct_, 1.0 );
    rhs( Sigmazz, 0 )->Update( - timefacfac* my::vderxy_(2, 2)                 , my::funct_, 1.0 );

    // stressbar-pressure coupling
    /*
                     /                    \
                    |                      |
                  - | tr(virt tau^e) , p I |
                    |                      |
                     \                    /
    */

    // K_sp

    K_su( Sigmaxx, Pres )->Update( -viscfac*timefacfac, bK_ss, 1.0 );
    K_su( Sigmayy, Pres )->Update( -viscfac*timefacfac, bK_ss, 1.0 );
    K_su( Sigmazz, Pres )->Update( -viscfac*timefacfac, bK_ss, 1.0 );

    // r_sp

    rhs( Sigmaxx, 0 )->Update( viscfac*timefacfac*press, my::funct_, 1.0 );
    rhs( Sigmayy, 0 )->Update( viscfac*timefacfac*press, my::funct_, 1.0 );
    rhs( Sigmazz, 0 )->Update( viscfac*timefacfac*press, my::funct_, 1.0 );
  }

  // integrate surface

  DRT::Element::LocationArray cutla( 1 );

  LINALG::Matrix<3,1> normal;
  LINALG::Matrix<3,1> x_side;

  bool fluidfluidcoupling = false;

  std::map<int, Teuchos::RCP<XFLUID::SideInterface<distype> > > side_impl;
  Teuchos::RCP<XFLUID::SideInterface<distype> > si;

  // find all the intersecting elements of actele
  std::set<int> begids;
  for (std::map<int,  std::vector<GEO::CUT::BoundaryCell*> >::const_iterator bc=bcells.begin();
       bc!=bcells.end(); ++bc )
  {
    int sid = bc->first;
    begids.insert(sid);
  }

  // map of boundary element gids, to coupling matrices Gsui and Guis (Cuiui = - Guis*Kss^-1*Gsui)
  std::map<int, std::vector<Epetra_SerialDenseMatrix> > Cuiui_coupling;

  // lm vector of all intersecting elements (boundary elements which intersect the current element)
  std::vector<int> patchelementslmv;
  std::vector<int> patchelementslmowner;
  for (std::set<int>::const_iterator bgid=begids.begin(); bgid!=begids.end(); ++bgid)
  {
    DRT::Element * side = cutdis.gElement(*bgid); // for each boundary element there is one corresponding side
    vector<int> patchlm;
    vector<int> patchlmowner;
    vector<int> patchlmstride;
    side->LocationVector(cutdis, patchlm, patchlmowner, patchlmstride);

    patchelementslmv.reserve( patchelementslmv.size() + patchlm.size());
    patchelementslmv.insert(patchelementslmv.end(), patchlm.begin(), patchlm.end());

    patchelementslmowner.reserve( patchelementslmowner.size() + patchlmowner.size());
    patchelementslmowner.insert( patchelementslmowner.end(), patchlmowner.begin(), patchlmowner.end());

    // get coupling matrices for the current side (boundary element)
    std::vector<Epetra_SerialDenseMatrix> & Cuiui_matrices = Cuiui_coupling[*bgid];

    Cuiui_matrices.resize(2);
    Cuiui_matrices[0].Reshape(my::nen_*6,patchlm.size()); //Gsui (coupling between background elements sigma and current side!)
    Cuiui_matrices[1].Reshape(patchlm.size(),my::nen_*6); //Guis
  }


  // coupling between domain and all! sides (boundary elements) that cut the element
  Epetra_SerialDenseMatrix Gsui(my::nen_*6,patchelementslmv.size());
  Epetra_SerialDenseMatrix Guis(patchelementslmv.size(),my::nen_*6);
  Epetra_SerialDenseMatrix InvKss(my::nen_*6,my::nen_*6);
  Epetra_SerialDenseMatrix GuisInvKss(patchelementslmv.size(),my::nen_*6);

  // map of side-element id and Gauss points
  for ( std::map<int, std::vector<DRT::UTILS::GaussIntegration> >::const_iterator i=bintpoints.begin();
        i!=bintpoints.end();
        ++i )
  {
    int sid = i->first;
    const std::vector<DRT::UTILS::GaussIntegration> & cutintpoints = i->second;

    std::map<int, std::vector<GEO::CUT::BoundaryCell*> >::const_iterator j = bcells.find( sid );
    if ( j==bcells.end() )
      dserror( "missing boundary cell" );

    const std::vector<GEO::CUT::BoundaryCell*> & bcs = j->second;
    if ( bcs.size()!=cutintpoints.size() )
      dserror( "boundary cell integration rules mismatch" );

    DRT::Element * side = cutdis.gElement( sid );
    side->LocationVector(cutdis,cutla,false);

    const int numnodes = side->NumNode();
    DRT::Node ** nodes = side->Nodes();
    Epetra_SerialDenseMatrix side_xyze( 3, numnodes );
    for ( int i=0; i<numnodes; ++i )
    {
      const double * x = nodes[i]->X();
      std::copy( x, x+3, &side_xyze( 0, i ) );
    }

    std::map<int,std::vector<Epetra_SerialDenseMatrix> >::iterator c = side_coupling.find( sid );

    std::vector<Epetra_SerialDenseMatrix> & side_matrices = c->second;

    if ( side_matrices.size()==3 )
      fluidfluidcoupling = true;


    if(fluidfluidcoupling)
    {
    	// coupling matrices between background element and one! side
        Epetra_SerialDenseMatrix & C_uiu  = side_matrices[0];
        Epetra_SerialDenseMatrix & C_uui  = side_matrices[1];
        Epetra_SerialDenseMatrix & rhC_ui = side_matrices[2];

        // coupling matrices between one side and itself via the element Kss
        std::map<int,std::vector<Epetra_SerialDenseMatrix> >::iterator c2 = Cuiui_coupling.find( sid );
        std::vector<Epetra_SerialDenseMatrix> & Cuiui_matrices = c2->second;
        Epetra_SerialDenseMatrix & eleGsui = Cuiui_matrices[0];
        Epetra_SerialDenseMatrix & eleGuis = Cuiui_matrices[1];
        Epetra_SerialDenseMatrix  eleGuisKssInv;

        si = XFLUID::SideInterface<distype>::Impl(side,C_uiu,C_uui,rhC_ui,eleGsui,eleGuis,side_xyze);
    }
    else
    {
        si = XFLUID::SideInterface<distype>::Impl(side,side_xyze);
    }

    side_impl[sid] = si;

    // get velocity at integration point of boundary dis

    si->eivel(cutdis,"ivelnp",cutla[0].lm_);
    si->addeidisp(cutdis,"idispnp",cutla[0].lm_,side_xyze);



    // loop gausspoints
    for ( std::vector<DRT::UTILS::GaussIntegration>::const_iterator i=cutintpoints.begin();
          i!=cutintpoints.end();
          ++i )
    {
      const DRT::UTILS::GaussIntegration & gi = *i;
      GEO::CUT::BoundaryCell * bc = bcs[i - cutintpoints.begin()]; // get the corresponding boundary cell

      //gi.Print();


      for ( DRT::UTILS::GaussIntegration::iterator iquad=gi.begin(); iquad!=gi.end(); ++iquad )
      {
        double drs = 0.0; // transformation factor between reference cell and linearized boundary cell
        LINALG::Matrix<3,1> x_gp_lin(true); // gp in xyz-system on linearized interface
        if(bc->Shape()==DRT::Element::tri3 || bc->Shape()==DRT::Element::quad4)
        {
#ifdef BOUNDARYCELL_TRANSFORMATION_OLD
          const LINALG::Matrix<2,1> eta( iquad.Point() ); // xi-coordinates with respect to side

          si->Evaluate(eta,x_side,normal,drs);

          const double fac = drs * iquad.Weight() * f3Parameter_->timefac_;

          // find element local position of gauss point at interface
          GEO::CUT::Position<distype> pos( xyze_, x_side );
          pos.Compute();
          const LINALG::Matrix<3,1> & rst = pos.LocalCoordinates();
  #else
          const LINALG::Matrix<2,1> eta( iquad.Point() ); // eta-coordinates with respect to cell

          normal.Clear();

          // get normal vector on linearized boundary cell, x-coordinates of gaussian point and surface transformation factor
          switch ( bc->Shape() )
          {
          case DRT::Element::tri3:
          {
              bc->Transform<DRT::Element::tri3>(eta, x_gp_lin, normal, drs);
            break;
          }
          case DRT::Element::quad4:
          {
              bc->Transform<DRT::Element::quad4>(eta, x_gp_lin, normal, drs);
            break;
          }
          default:
            throw std::runtime_error( "unsupported integration cell type" );
          }
        }
        else if(bc->Shape()==DRT::Element::dis_none)
        {
          drs = 1.0;
          normal = bc->GetNormalVector();
          const double* gpcord = iquad.Point();
          for (int idim=0;idim<3;idim++)
          {
             x_gp_lin(idim,0) = gpcord[idim];
          }
        }

        const double fac = drs * iquad.Weight() * my::f3Parameter_->timefac_;

        // find element local position of gauss point
        GEO::CUT::Position<distype> pos( my::xyze_, x_gp_lin );
        pos.Compute();
        const LINALG::Matrix<3,1> & rst = pos.LocalCoordinates();

        // project gaussian point from linearized interface to warped side (get local side coordinates)
        LINALG::Matrix<2,1> xi_side(true);
        si->ProjectOnSide(x_gp_lin, x_side, xi_side);

#endif


        // evaluate shape functions
        DRT::UTILS::shape_function<distype>( rst, my::funct_ );

        // evaluate shape functions and derivatives at integration point
        //EvalShapeFuncAndDerivsAtIntPoint(iquad,eid);

        // get velocity at integration point
        // (values at n+alpha_F for generalized-alpha scheme, n+1 otherwise)
        my::velint_.Multiply(evelaf,my::funct_);

        // get velocity derivatives at integration point
        // (values at n+alpha_F for generalized-alpha scheme, n+1 otherwise)
        //vderxy_.MultiplyNT(evelaf,derxy_);

        // get pressure at integration point
        // (value at n+alpha_F for generalized-alpha scheme, n+1 otherwise)
        //double press = funct_.Dot(epreaf);

        bK_ss.MultiplyNT( my::funct_, my::funct_ );

               /*                      \
            - |  (virt tau) * n^f , Du  |
               \                      */

        // G_su

        K_su( Sigmaxx, Velx )->Update( -fac*normal(0), bK_ss, 1.0 );
        K_su( Sigmaxy, Velx )->Update( -fac*normal(1), bK_ss, 1.0 );
        K_su( Sigmaxz, Velx )->Update( -fac*normal(2), bK_ss, 1.0 );
        K_su( Sigmayx, Vely )->Update( -fac*normal(0), bK_ss, 1.0 );
        K_su( Sigmayy, Vely )->Update( -fac*normal(1), bK_ss, 1.0 );
        K_su( Sigmayz, Vely )->Update( -fac*normal(2), bK_ss, 1.0 );
        K_su( Sigmazx, Velz )->Update( -fac*normal(0), bK_ss, 1.0 );
        K_su( Sigmazy, Velz )->Update( -fac*normal(1), bK_ss, 1.0 );
        K_su( Sigmazz, Velz )->Update( -fac*normal(2), bK_ss, 1.0 );

        rhs( Sigmaxx, 0 )->Update( fac*normal(0)*my::velint_(0), my::funct_, 1.0 );
        rhs( Sigmaxy, 0 )->Update( fac*normal(1)*my::velint_(0), my::funct_, 1.0 );
        rhs( Sigmaxz, 0 )->Update( fac*normal(2)*my::velint_(0), my::funct_, 1.0 );
        rhs( Sigmayx, 0 )->Update( fac*normal(0)*my::velint_(1), my::funct_, 1.0 );
        rhs( Sigmayy, 0 )->Update( fac*normal(1)*my::velint_(1), my::funct_, 1.0 );
        rhs( Sigmayz, 0 )->Update( fac*normal(2)*my::velint_(1), my::funct_, 1.0 );
        rhs( Sigmazx, 0 )->Update( fac*normal(0)*my::velint_(2), my::funct_, 1.0 );
        rhs( Sigmazy, 0 )->Update( fac*normal(1)*my::velint_(2), my::funct_, 1.0 );
        rhs( Sigmazz, 0 )->Update( fac*normal(2)*my::velint_(2), my::funct_, 1.0 );

if(!fluidfluidcoupling)
{
        /*                   _  \
       |  (virt tau) * n^f , u   |
        \                      */


        LINALG::Matrix<my::nsd_,1> velint_WDBC(true);
        si->get_vel_WeakDBC(velint_WDBC);

        rhs( Sigmaxx, 0 )->Update( -fac*normal(0)*velint_WDBC(0), my::funct_, 1.0 );
        rhs( Sigmaxy, 0 )->Update( -fac*normal(1)*velint_WDBC(0), my::funct_, 1.0 );
        rhs( Sigmaxz, 0 )->Update( -fac*normal(2)*velint_WDBC(0), my::funct_, 1.0 );
        rhs( Sigmayx, 0 )->Update( -fac*normal(0)*velint_WDBC(1), my::funct_, 1.0 );
        rhs( Sigmayy, 0 )->Update( -fac*normal(1)*velint_WDBC(1), my::funct_, 1.0 );
        rhs( Sigmayz, 0 )->Update( -fac*normal(2)*velint_WDBC(1), my::funct_, 1.0 );
        rhs( Sigmazx, 0 )->Update( -fac*normal(0)*velint_WDBC(2), my::funct_, 1.0 );
        rhs( Sigmazy, 0 )->Update( -fac*normal(1)*velint_WDBC(2), my::funct_, 1.0 );
        rhs( Sigmazz, 0 )->Update( -fac*normal(2)*velint_WDBC(2), my::funct_, 1.0 );
}

               /*               \
            - |  v , Dtau * n^f  |
               \               */

        // G_us
        K_us( Velx, Sigmaxx )->Update( -fac*normal(0), bK_ss, 1.0 );
        K_us( Velx, Sigmaxy )->Update( -fac*normal(1), bK_ss, 1.0 );
        K_us( Velx, Sigmaxz )->Update( -fac*normal(2), bK_ss, 1.0 );
        K_us( Vely, Sigmayx )->Update( -fac*normal(0), bK_ss, 1.0 );
        K_us( Vely, Sigmayy )->Update( -fac*normal(1), bK_ss, 1.0 );
        K_us( Vely, Sigmayz )->Update( -fac*normal(2), bK_ss, 1.0 );
        K_us( Velz, Sigmazx )->Update( -fac*normal(0), bK_ss, 1.0 );
        K_us( Velz, Sigmazy )->Update( -fac*normal(1), bK_ss, 1.0 );
        K_us( Velz, Sigmazz )->Update( -fac*normal(2), bK_ss, 1.0 );


        // evaluate the derivatives of shape functions
        DRT::UTILS::shape_function_deriv1<distype>(rst,my::deriv_);
        my::xjm_.MultiplyNT(my::deriv_,my::xyze_);
        my::det_ = my::xji_.Invert(my::xjm_);

        // compute global first derivates
        my::derxy_.Multiply(my::xji_,my::deriv_);


        // get velocity derivatives at integration point
        // (values at n+alpha_F for generalized-alpha scheme, n+1 otherwise)
        my::vderxy_.MultiplyNT(evelaf,my::derxy_);

        // get pressure at integration point
        // (value at n+alpha_F for generalized-alpha scheme, n+1 otherwise)
        double press = my::funct_.Dot(epreaf);


        // calculate interface forces
        if(!fluidfluidcoupling)
        {

          // compute the stresses at the current Gaussian point for computing the interface force
          LINALG::Matrix<my::nsd_,my::nsd_> eps(true);
          for(int i=0; i<my::nsd_; i++)
          {
            for(int j=0; j<my::nsd_; j++)
            {
              eps(i,j) = 0.5 * (my::vderxy_(i,j) + my::vderxy_(j,i));
            }
          }


          // t = ( -pI + 2mu eps(u) )*n^f
          LINALG::Matrix<my::nsd_,1> traction (true);
          traction.Multiply(eps, normal);
          traction.Scale(2.0*my::visceff_);

          // add the pressure part
          traction.Update( -press, normal, 1.0);

          // we need normal vector on fluid
          traction.Scale(-1.0);

          double surf_fac = drs*iquad.Weight();

          si->buildInterfaceForce(iforcecol, cutdis, cutla[0].lm_, traction, surf_fac );


        }

//
//        if(stress_with_l2_proj == false)
//        {
//
          /*                  \       /          i      \
       - |  [ v ], - {Dp}*n    | = - | [ v ], { p }* n   |
          \                   /       \                */
//
//          //-----------------------------------------------
//          //    + (v1, k1 *(Dp1)*n)
//          //-----------------------------------------------
//
//          LINALG::Matrix<my::nen_,1> my::funct_timefacfac(true);
//          funct_timefacfac.Update(fac,my::funct_,0.0);
//
//          LINALG::Matrix<my::nen_,my::nen_> funct_dyad_timefacfac(true);
//          LINALG::Matrix<my::nen_,my::nen_> funct_dyad_k1_timefacfac(true);
//          funct_dyad_timefacfac.MultiplyNT(funct_timefacfac, funct_);
//
//          for(int ir = 0; ir<my::nen_; ir++)
//          {
//            int idVelx = ir*(my::nsd_+1) + 0;
//            int idVely = ir*(my::nsd_+1) + 1;
//            int idVelz = ir*(my::nsd_+1) + 2;
//
//            // (v,Dp*n)
//            for(int ic =0; ic<my::nen_; ic++)
//            {
//              int iPres = ic*(my::nsd_+1)+3;
//
//              elemat1_epetra(idVelx, iPres) += funct_dyad_timefacfac(ir,ic)*normal(Velx);
//              elemat1_epetra(idVely, iPres) += funct_dyad_timefacfac(ir,ic)*normal(Vely);
//              elemat1_epetra(idVelz, iPres) += funct_dyad_timefacfac(ir,ic)*normal(Velz);
//            }
//
//            // -(v,p*n)
//            double funct_timefacfac_press = funct_timefacfac(ir)*press;
//            elevec1_epetra(idVelx,0) -= funct_timefacfac_press*normal(Velx);
//            elevec1_epetra(idVely,0) -= funct_timefacfac_press*normal(Vely);
//            elevec1_epetra(idVelz,0) -= funct_timefacfac_press*normal(Velz);
//          }
//
//
//
          /*                           \       /                   i      \
       - |  [ v ],  { 2mu eps(u) }*n    | = + | [ v ],  { 2mu eps(u ) }*n  |
          \                            /       \                         */
//
//          //-----------------------------------------------
//          //    - (v1, (2*k1*mu1) *eps(Du1)*n)
//          //-----------------------------------------------
//
//
//          LINALG::Matrix<my::nen_,1> e_funct_visc1_timefacfac(true);
//          e_funct_visc1_timefacfac.Update(2.0 * fac * visceff_, funct_, 0.0);
//
//          //            LINALG::Matrix<side_my::nen_,1> s_funct_visc_timefacfac(true);
//          //            s_funct_visc_timefacfac.Update(k2mu2_fac, side_funct_, 0.0);
//
//          for(int ir = 0; ir<my::nen_; ir++)
//          {
//            int idVelx = ir*(my::nsd_+1) + 0;
//            int idVely = ir*(my::nsd_+1) + 1;
//            int idVelz = ir*(my::nsd_+1) + 2;
//
//
//            for(int ic =0; ic<my::nen_; ic++)
//            {
//              int iVelx = ic*(my::nsd_+1)+0;
//              int iVely = ic*(my::nsd_+1)+1;
//              int iVelz = ic*(my::nsd_+1)+2;
//
//
//              // - (v1, (2*k1*mu1) *eps(Du1)*n)
//
//              //(x,x)
//              elemat1_epetra(idVelx, iVelx) -= e_funct_visc1_timefacfac(ir)*(         normal(Velx)*derxy_(Velx,ic)
//                                                                              + 0.5 * normal(Vely)*derxy_(Vely,ic)
//                                                                              + 0.5 * normal(Velz)*derxy_(Velz,ic)  );
//              //(x,y)
//              elemat1_epetra(idVelx, iVely) -= e_funct_visc1_timefacfac(ir)*    0.5 * normal(Vely)*derxy_(Velx,ic);
//              //(x,z)
//              elemat1_epetra(idVelx, iVelz) -= e_funct_visc1_timefacfac(ir)*    0.5 * normal(Velz)*derxy_(Velx,ic);
//
//              //(y,x)
//              elemat1_epetra(idVely, iVelx) -= e_funct_visc1_timefacfac(ir)*    0.5 * normal(Velx)*derxy_(Vely,ic);
//              //(y,y)
//              elemat1_epetra(idVely, iVely) -= e_funct_visc1_timefacfac(ir)*(   0.5 * normal(Velx)*derxy_(Velx,ic)
//                                                                              +       normal(Vely)*derxy_(Vely,ic)
//                                                                              + 0.5 * normal(Velz)*derxy_(Velz,ic)  );
//              //(y,z)
//              elemat1_epetra(idVely, iVelz) -= e_funct_visc1_timefacfac(ir)*    0.5 * normal(Velz)*derxy_(Vely,ic);
//
//              //(z,x)
//              elemat1_epetra(idVelz, iVelx) -= e_funct_visc1_timefacfac(ir)*    0.5 * normal(Velx)*derxy_(Velz,ic);
//              //(z,y)
//              elemat1_epetra(idVelz, iVely) -= e_funct_visc1_timefacfac(ir)*    0.5 * normal(Vely)*derxy_(Velz,ic);
//              //(z,z)
//              elemat1_epetra(idVelz, iVelz) -= e_funct_visc1_timefacfac(ir)*(   0.5 * normal(Velx)*derxy_(Velx,ic)
//                                                                              + 0.5 * normal(Vely)*derxy_(Vely,ic)
//                                                                              +       normal(Velz)*derxy_(Velz,ic)  );
//            }
//
//            // - (v1, (2*k1*mu1) *eps(Du1)*n)
//            elevec1_epetra(idVelx) += e_funct_visc1_timefacfac(ir)*(            vderxy_(Velx,Velx)                      *normal(Velx)
//                                                                      + 0.5 * ( vderxy_(Velx,Vely) + vderxy_(Vely,Velx))*normal(Vely)
//                                                                      + 0.5 * ( vderxy_(Velx,Velz) + vderxy_(Velz,Velx))*normal(Velz)  );
//            elevec1_epetra(idVely) += e_funct_visc1_timefacfac(ir)*(    0.5 * ( vderxy_(Vely,Velx) + vderxy_(Velx,Vely))*normal(Velx)
//                                                                      +         vderxy_(Vely,Vely)                      *normal(Vely)
//                                                                      + 0.5 * ( vderxy_(Vely,Velz) + vderxy_(Velz,Vely))*normal(Velz)  );
//            elevec1_epetra(idVelz) += e_funct_visc1_timefacfac(ir)*(    0.5 * ( vderxy_(Velz,Velx) + vderxy_(Velx,Velz))*normal(Velx)
//                                                                      + 0.5 * ( vderxy_(Velz,Vely) + vderxy_(Vely,Velz))*normal(Vely)
//                                                                      +         vderxy_(Velz,Velz)                      *normal(Velz)  );
//
//          }




//        }




#if 0
              /*                      \
             |  (virt tau) * n^f , Dui |
              \                      */

        // G_si

        patchassembler.template Matrix<Sigmaxx,Velxiface>(*(couplmats.Gsui_uncond), shp_tau.d0, fac*normal(0), shp_iface_d0);
        patchassembler.template Matrix<Sigmaxy,Velxiface>(*(couplmats.Gsui_uncond), shp_tau.d0, fac*normal(1), shp_iface_d0);
        patchassembler.template Matrix<Sigmaxz,Velxiface>(*(couplmats.Gsui_uncond), shp_tau.d0, fac*normal(2), shp_iface_d0);
        patchassembler.template Matrix<Sigmayx,Velyiface>(*(couplmats.Gsui_uncond), shp_tau.d0, fac*normal(0), shp_iface_d0);
        patchassembler.template Matrix<Sigmayy,Velyiface>(*(couplmats.Gsui_uncond), shp_tau.d0, fac*normal(1), shp_iface_d0);
        patchassembler.template Matrix<Sigmayz,Velyiface>(*(couplmats.Gsui_uncond), shp_tau.d0, fac*normal(2), shp_iface_d0);
        patchassembler.template Matrix<Sigmazx,Velziface>(*(couplmats.Gsui_uncond), shp_tau.d0, fac*normal(0), shp_iface_d0);
        patchassembler.template Matrix<Sigmazy,Velziface>(*(couplmats.Gsui_uncond), shp_tau.d0, fac*normal(1), shp_iface_d0);
        patchassembler.template Matrix<Sigmazz,Velziface>(*(couplmats.Gsui_uncond), shp_tau.d0, fac*normal(2), shp_iface_d0);

        assembler.template Vector<Sigmaxx>(shp_tau.d0, -fac*normal(0)*interface_gpvelnp(0));
        assembler.template Vector<Sigmaxy>(shp_tau.d0, -fac*normal(1)*interface_gpvelnp(0));
        assembler.template Vector<Sigmaxz>(shp_tau.d0, -fac*normal(2)*interface_gpvelnp(0));
        assembler.template Vector<Sigmayx>(shp_tau.d0, -fac*normal(0)*interface_gpvelnp(1));
        assembler.template Vector<Sigmayy>(shp_tau.d0, -fac*normal(1)*interface_gpvelnp(1));
        assembler.template Vector<Sigmayz>(shp_tau.d0, -fac*normal(2)*interface_gpvelnp(1));
        assembler.template Vector<Sigmazx>(shp_tau.d0, -fac*normal(0)*interface_gpvelnp(2));
        assembler.template Vector<Sigmazy>(shp_tau.d0, -fac*normal(1)*interface_gpvelnp(2));
        assembler.template Vector<Sigmazz>(shp_tau.d0, -fac*normal(2)*interface_gpvelnp(2));

        if (monolithic_FSI)
        {
          const double facvelx = monolithic_FSI ? fac*(gpvelnp(0) - interface_gpvelnp(0)) : 0.0;
          const double facvely = monolithic_FSI ? fac*(gpvelnp(1) - interface_gpvelnp(1)) : 0.0;
          const double facvelz = monolithic_FSI ? fac*(gpvelnp(2) - interface_gpvelnp(2)) : 0.0;
          patchassembler.template Matrix<Sigmaxx,Dispxiface>(*couplmats.GNsdi_uncond, shp_tau.d0, -facvelx, normalderiv.dnxdx);
          patchassembler.template Matrix<Sigmaxy,Dispxiface>(*couplmats.GNsdi_uncond, shp_tau.d0, -facvelx, normalderiv.dnydx);
          patchassembler.template Matrix<Sigmaxz,Dispxiface>(*couplmats.GNsdi_uncond, shp_tau.d0, -facvelx, normalderiv.dnzdx);
          patchassembler.template Matrix<Sigmayx,Dispyiface>(*couplmats.GNsdi_uncond, shp_tau.d0, -facvely, normalderiv.dnxdy);
          patchassembler.template Matrix<Sigmayy,Dispyiface>(*couplmats.GNsdi_uncond, shp_tau.d0, -facvely, normalderiv.dnydy);
          patchassembler.template Matrix<Sigmayz,Dispyiface>(*couplmats.GNsdi_uncond, shp_tau.d0, -facvely, normalderiv.dnzdy);
          patchassembler.template Matrix<Sigmazx,Dispziface>(*couplmats.GNsdi_uncond, shp_tau.d0, -facvelz, normalderiv.dnxdz);
          patchassembler.template Matrix<Sigmazy,Dispziface>(*couplmats.GNsdi_uncond, shp_tau.d0, -facvelz, normalderiv.dnydz);
          patchassembler.template Matrix<Sigmazz,Dispziface>(*couplmats.GNsdi_uncond, shp_tau.d0, -facvelz, normalderiv.dnzdz);
        }
#endif

        if ( fluidfluidcoupling )
          si->buildCouplingMatrices(normal,fac,my::funct_,rhs);

//         if (monolithic_FSI)
//         {
//           const double nfac1 = monolithic_FSI ? fac : 0.0;
//           patchassembler.template Matrix<Velx,Dispxiface>(*couplmats.GNudi_uncond, shp.d0, -nfac1, tau_times_nderiv.xx);
//           patchassembler.template Matrix<Velx,Dispyiface>(*couplmats.GNudi_uncond, shp.d0, -nfac1, tau_times_nderiv.xy);
//           patchassembler.template Matrix<Velx,Dispziface>(*couplmats.GNudi_uncond, shp.d0, -nfac1, tau_times_nderiv.xz);
//           patchassembler.template Matrix<Vely,Dispxiface>(*couplmats.GNudi_uncond, shp.d0, -nfac1, tau_times_nderiv.yx);
//           patchassembler.template Matrix<Vely,Dispyiface>(*couplmats.GNudi_uncond, shp.d0, -nfac1, tau_times_nderiv.yy);
//           patchassembler.template Matrix<Vely,Dispziface>(*couplmats.GNudi_uncond, shp.d0, -nfac1, tau_times_nderiv.yz);
//           patchassembler.template Matrix<Velz,Dispxiface>(*couplmats.GNudi_uncond, shp.d0, -nfac1, tau_times_nderiv.zx);
//           patchassembler.template Matrix<Velz,Dispyiface>(*couplmats.GNudi_uncond, shp.d0, -nfac1, tau_times_nderiv.zy);
//           patchassembler.template Matrix<Velz,Dispziface>(*couplmats.GNudi_uncond, shp.d0, -nfac1, tau_times_nderiv.zz);
//         }

//         LINALG::Matrix<nsd,1> disctau_times_nf;
//         disctau_times_nf.Multiply(tau,normal);
//         //cout << "sigmaijnj : " << disctau_times_n << endl;
//         assembler.template Vector<Velx>(shp.d0, fac*disctau_times_nf(0));
//         assembler.template Vector<Vely>(shp.d0, fac*disctau_times_nf(1));
//         assembler.template Vector<Velz>(shp.d0, fac*disctau_times_nf(2));

        // integrate the force boundary

        // here the interface force is integrated
        // this is done using test
        // shape functions of the boundary mesh
        // hence, we can't use the local assembler here
//        for (size_t inode = 0; inode < numnode_boundary; ++inode)
//        {
//          force_boundary(0,inode) += funct_boundary(inode) * -(disctau_times_nf(0) * fac);
//          force_boundary(1,inode) += funct_boundary(inode) * -(disctau_times_nf(1) * fac);
//          force_boundary(2,inode) += funct_boundary(inode) * -(disctau_times_nf(2) * fac);
//        }

#if 0
              /*                  \
             |  v^i , Dtau * n^f   |
              \                  */


        patchassembler.template Matrix<Velxiface,Sigmaxx>(*(couplmats.Guis_uncond), shp_iface_d0, fac*normal(0), shp_tau.d0);
        patchassembler.template Matrix<Velxiface,Sigmaxy>(*(couplmats.Guis_uncond), shp_iface_d0, fac*normal(1), shp_tau.d0);
        patchassembler.template Matrix<Velxiface,Sigmaxz>(*(couplmats.Guis_uncond), shp_iface_d0, fac*normal(2), shp_tau.d0);
        patchassembler.template Matrix<Velyiface,Sigmayx>(*(couplmats.Guis_uncond), shp_iface_d0, fac*normal(0), shp_tau.d0);
        patchassembler.template Matrix<Velyiface,Sigmayy>(*(couplmats.Guis_uncond), shp_iface_d0, fac*normal(1), shp_tau.d0);
        patchassembler.template Matrix<Velyiface,Sigmayz>(*(couplmats.Guis_uncond), shp_iface_d0, fac*normal(2), shp_tau.d0);
        patchassembler.template Matrix<Velziface,Sigmazx>(*(couplmats.Guis_uncond), shp_iface_d0, fac*normal(0), shp_tau.d0);
        patchassembler.template Matrix<Velziface,Sigmazy>(*(couplmats.Guis_uncond), shp_iface_d0, fac*normal(1), shp_tau.d0);
        patchassembler.template Matrix<Velziface,Sigmazz>(*(couplmats.Guis_uncond), shp_iface_d0, fac*normal(2), shp_tau.d0);

        if (monolithic_FSI)
        {
          const double nfac1 = monolithic_FSI ? fac : 0.0;
          patchassembler.template Matrix<Dispxiface,Dispxiface>(*couplmats.GNdidi_uncond, shp_iface_d0, -nfac1, tau_times_nderiv.xx);
          patchassembler.template Matrix<Dispxiface,Dispyiface>(*couplmats.GNdidi_uncond, shp_iface_d0, -nfac1, tau_times_nderiv.xy);
          patchassembler.template Matrix<Dispxiface,Dispziface>(*couplmats.GNdidi_uncond, shp_iface_d0, -nfac1, tau_times_nderiv.xz);
          patchassembler.template Matrix<Dispyiface,Dispxiface>(*couplmats.GNdidi_uncond, shp_iface_d0, -nfac1, tau_times_nderiv.yx);
          patchassembler.template Matrix<Dispyiface,Dispyiface>(*couplmats.GNdidi_uncond, shp_iface_d0, -nfac1, tau_times_nderiv.yy);
          patchassembler.template Matrix<Dispyiface,Dispziface>(*couplmats.GNdidi_uncond, shp_iface_d0, -nfac1, tau_times_nderiv.yz);
          patchassembler.template Matrix<Dispziface,Dispxiface>(*couplmats.GNdidi_uncond, shp_iface_d0, -nfac1, tau_times_nderiv.zx);
          patchassembler.template Matrix<Dispziface,Dispyiface>(*couplmats.GNdidi_uncond, shp_iface_d0, -nfac1, tau_times_nderiv.zy);
          patchassembler.template Matrix<Dispziface,Dispziface>(*couplmats.GNdidi_uncond, shp_iface_d0, -nfac1, tau_times_nderiv.zz);
        }

        patchassembler.template Vector<Velxiface>(*(couplmats.rhsui_uncond), shp_iface_d0, -fac*disctau_times_nf(0));
        patchassembler.template Vector<Velyiface>(*(couplmats.rhsui_uncond), shp_iface_d0, -fac*disctau_times_nf(1));
        patchassembler.template Vector<Velziface>(*(couplmats.rhsui_uncond), shp_iface_d0, -fac*disctau_times_nf(2));

        // here the interface force is integrated
        // this is done using test
        // shape functions of the boundary mesh
        // hence, we can't use the local assembler here
        for (size_t inode = 0; inode < numnode_boundary; ++inode)
        {
          force_boundary(0,inode) += funct_boundary(inode) * -(disctau_times_nf(0) * fac);
          force_boundary(1,inode) += funct_boundary(inode) * -(disctau_times_nf(1) * fac);
          force_boundary(2,inode) += funct_boundary(inode) * -(disctau_times_nf(2) * fac);
        }
#endif

#if 0
        // TODO: timefac not used here?!?!?
  double timefacfac = fac;

  // funct_ * timefac * fac
  LINALG::Matrix<my::nen_,1> funct_timefacfac(true);
  funct_timefacfac.Update(timefacfac,funct_,0.0);

  // funct_ * timefac * fac * funct_ (dyadic product)
  LINALG::Matrix<my::nen_,my::nen_> funct_dyad_timefacfac(true);
  funct_dyad_timefacfac.MultiplyNT(funct_timefacfac,funct_);

			  // convective stabilization
				 /*                           \        /                       i   _     \
			    |  gamma/h_K *  v*n , Du*n     | =  - |   gamma/h_K *  v*n , (u  - u)*n   |
				 \                            /        \                                */

			  const double gamma_conv = 100.0;
			  const double h_K = 1.0/20.0;

			  const double stab_fac_conv = gamma_conv/h_K;

			  for(int ir=0; ir<my::nen_; ir++)
			  {
					int idVelx = ir*(my::nsd_+1) + 0;
					int idVely = ir*(my::nsd_+1) + 1;
					int idVelz = ir*(my::nsd_+1) + 2;

					// (stab * v, Du)
					for(int ic=0; ic<my::nen_; ic++)
					{
						int iVelx = ic*(my::nsd_+1)+0;
						int iVely = ic*(my::nsd_+1)+1;
						int iVelz = ic*(my::nsd_+1)+2;

						elemat1_epetra(idVelx, iVelx) += funct_dyad_timefacfac(ir,ic)*stab_fac_conv*normal(Velx)*normal(Velx);
						elemat1_epetra(idVelx, iVely) += funct_dyad_timefacfac(ir,ic)*stab_fac_conv*normal(Velx)*normal(Vely);
						elemat1_epetra(idVelx, iVelz) += funct_dyad_timefacfac(ir,ic)*stab_fac_conv*normal(Velx)*normal(Velz);

						elemat1_epetra(idVely, iVelx) += funct_dyad_timefacfac(ir,ic)*stab_fac_conv*normal(Vely)*normal(Velx);
						elemat1_epetra(idVely, iVely) += funct_dyad_timefacfac(ir,ic)*stab_fac_conv*normal(Vely)*normal(Vely);
						elemat1_epetra(idVely, iVelz) += funct_dyad_timefacfac(ir,ic)*stab_fac_conv*normal(Vely)*normal(Velz);

						elemat1_epetra(idVelz, iVelx) += funct_dyad_timefacfac(ir,ic)*stab_fac_conv*normal(Velz)*normal(Velx);
						elemat1_epetra(idVelz, iVely) += funct_dyad_timefacfac(ir,ic)*stab_fac_conv*normal(Velz)*normal(Vely);
						elemat1_epetra(idVelz, iVelz) += funct_dyad_timefacfac(ir,ic)*stab_fac_conv*normal(Velz)*normal(Velz);
					}

					double velint_normal = my::velint_.Dot(normal);

					// -(stab * v*n, u*n)
					elevec1_epetra(idVelx) -= funct_timefacfac(ir)*normal(Velx)*stab_fac_conv*velint_normal;
					elevec1_epetra(idVely) -= funct_timefacfac(ir)*normal(Vely)*stab_fac_conv*velint_normal;
					elevec1_epetra(idVelz) -= funct_timefacfac(ir)*normal(Velz)*stab_fac_conv*velint_normal;


//					double velint_WDBC_normal = velint_WDBC.Dot(normal);
//					// +(stab * v*n, u_DBC*n)
//					elevec1_epetra(idVelx) += funct_timefacfac(ir)*normal(Velx)*stab_fac_conv*velint_WDBC_normal;
//					elevec1_epetra(idVely) += funct_timefacfac(ir)*normal(Vely)*stab_fac_conv*velint_WDBC_normal;
//					elevec1_epetra(idVelz) += funct_timefacfac(ir)*normal(Velz)*stab_fac_conv*velint_WDBC_normal;
			  }
#endif

      }
    }
  }

  // construct views
  LINALG::Matrix<(my::nsd_+1)*my::nen_,(my::nsd_+1)*my::nen_> elemat1(elemat1_epetra,true);
  LINALG::Matrix<(my::nsd_+1)*my::nen_,            1> elevec1(elevec1_epetra,true);

#if 0

  // DEBUG

   LINALG::Matrix<my::nen_,my::nen_> two_invbK_ss( true );
   two_invbK_ss.Update( 2, invbK_ss, 0.0 );
   LINALG::BlockMatrix<LINALG::Matrix<my::nen_,my::nen_>,6,6> K_ss;

   K_ss.AddView( Sigmaxx, Sigmaxx,     invbK_ss );
   K_ss.AddView( Sigmaxy, Sigmaxy, two_invbK_ss );
   K_ss.AddView( Sigmaxz, Sigmaxz, two_invbK_ss );
   K_ss.AddView( Sigmayy, Sigmayy,     invbK_ss );
   K_ss.AddView( Sigmayz, Sigmayz, two_invbK_ss );
   K_ss.AddView( Sigmazz, Sigmazz,     invbK_ss );

   LINALG::Matrix<4*my::nen_, 6*my::nen_> real_K_us( true );
   LINALG::Matrix<6*my::nen_, 4*my::nen_> real_K_su( true );
   LINALG::Matrix<6*my::nen_, 6*my::nen_> real_K_ss( true );

   K_us.template AssembleTo<4*my::nen_, 6*my::nen_>( real_K_us, 1. );
   K_su.template AssembleTo<6*my::nen_, 4*my::nen_>( real_K_su, 1. );
   K_ss.template AssembleTo<6*my::nen_, 6*my::nen_>( real_K_ss, 1. );

   std::cout << real_K_us << "*\n" << real_K_su << "*\n" << real_K_ss << "***\n";

   LINALG::FixedSizeSerialDenseSolver<6*my::nen_,6*my::nen_> solver;
   solver.SetMatrix( real_K_ss );
   solver.Invert();

   LINALG::Matrix<(my::nsd_+1)*my::nen_,6*my::nen_> K_iK( true );
   LINALG::Matrix<(my::nsd_+1)*my::nen_,(my::nsd_+1)*my::nen_> extK( true );
   //LINALG::Matrix<(my::nsd_+1)*my::nen_,   1> extrhs;

   K_iK  .Multiply( real_K_us, real_K_ss );
   extK  .Multiply( K_iK, real_K_su );
   //extrhs.Multiply( K_iK, rhs );

   elemat1.Update( -1, extK, 1 );
   //elevec1.Update( -1, extrhs, 1 );

#else

  // invert block mass matrix
  LINALG::FixedSizeSerialDenseSolver<my::nen_,my::nen_> solver;
  solver.SetMatrix( invbK_ss );
  solver.Invert();

  half_invbK_ss.Update( 0.5, invbK_ss, 0.0 );

  invK_ss.AddView( Sigmaxx, Sigmaxx,      invbK_ss );
  invK_ss.AddView( Sigmaxy, Sigmaxy, half_invbK_ss );
  invK_ss.AddView( Sigmaxz, Sigmaxz, half_invbK_ss );
  invK_ss.AddView( Sigmayy, Sigmayy,      invbK_ss );
  invK_ss.AddView( Sigmayz, Sigmayz, half_invbK_ss );
  invK_ss.AddView( Sigmazz, Sigmazz,      invbK_ss );

  LINALG::BlockMatrix<LINALG::Matrix<my::nen_,my::nen_>,(my::nsd_+1),6>        K_iK;
  LINALG::BlockMatrix<LINALG::Matrix<my::nen_,my::nen_>,(my::nsd_+1),(my::nsd_+1)> extK;
  LINALG::BlockMatrix<LINALG::Matrix<my::nen_,   1>,(my::nsd_+1),1>        extrhs;

  K_iK  .Multiply( K_us, invK_ss );
  extK  .Multiply( K_iK, K_su );
  extrhs.Multiply( K_iK, rhs );

  for ( unsigned icb=0; icb<my::nsd_+1; ++icb )
  {
    for ( unsigned irb=0; irb<my::nsd_+1; ++irb )
    {
      if ( extK.IsUsed( irb, icb ) )
      {
        LINALG::Matrix<my::nen_,my::nen_> & local_extK = *extK( irb, icb );
        for ( int ic=0; ic<my::nen_; ++ic )
        {
          unsigned c = ( my::nsd_+1 )*ic + icb;
          for ( int ir=0; ir<my::nen_; ++ir )
          {
            unsigned r = ( my::nsd_+1 )*ir + irb;
            elemat1( r, c ) -= local_extK( ir, ic );
          }
        }
      }
    }
  }

  for ( unsigned irb=0; irb<my::nsd_+1; ++irb )
  {
    if ( extrhs.IsUsed( irb, 0 ) )
    {
      LINALG::Matrix<my::nen_,1> & local_extrhs = *extrhs( irb, 0 );
      for ( int ir=0; ir<my::nen_; ++ir )
      {
        unsigned r = ( my::nsd_+1 )*ir + irb;
        elevec1( r, 0 ) -= local_extrhs( ir, 0 );
      }
    }
  }

  if ( fluidfluidcoupling )
  {
    // build fluid-fluid matrices
    for (typename std::map<int, Teuchos::RCP<XFLUID::SideInterface<distype> > >::iterator i=side_impl.begin();  i!=side_impl.end(); ++i)
    {
      XFLUID::SideInterface<distype> * si = &*i->second;
      si->buildFinalCouplingMatrices(invK_ss,K_iK,K_su,rhs);
    }

    // build InvKss
    for ( unsigned icb=0; icb<6; ++icb )
    {
      for ( unsigned irb=0; irb<6; ++irb )
      {
        if ( invK_ss.IsUsed( irb, icb ) )
        {
          LINALG::Matrix<my::nen_,my::nen_> & local_invK_ss = *invK_ss( irb, icb );
          for ( int ic=0; ic<my::nen_; ++ic )
          {
            int c = ( 6 )*ic + icb;
            for ( int ir=0; ir<my::nen_; ++ir )
            {
              int r = ( 6 )*ir + irb;
              InvKss( r, c ) -= local_invK_ss( ir, ic );
            }
          }
        }
      }
    }

   // build Gsui and Guis
    int ipatchsizesbefore = 0;
    for (std::map<int, std::vector<Epetra_SerialDenseMatrix> >::const_iterator m=Cuiui_coupling.begin();
         m!=Cuiui_coupling.end(); ++m)
    {
      int bid = m->first;
      std::vector<Epetra_SerialDenseMatrix> & Cuiui_mats = Cuiui_coupling[bid];
      // assemble Gsui
      for ( int icb=0; icb<Cuiui_mats[0].N(); ++icb ) // Cuiui includes only ui,ui coupling, not (ui,p) ...
      {
        for ( unsigned irb=0; irb<6*my::nen_; ++irb )
        {
          Gsui(irb,icb+ipatchsizesbefore) = Cuiui_mats[0](irb,icb);
        }
      }

      // assemble Guis
      for ( unsigned icb=0; icb<6*my::nen_; ++icb )
      {
        for ( int irb=0; irb<Cuiui_mats[1].M(); ++irb )
        {
          Guis(irb+ipatchsizesbefore,icb) = Cuiui_mats[1](irb,icb);
        }
      }

      ipatchsizesbefore += Cuiui_mats[0].N();

    }


    GuisInvKss.Multiply('N','N',1.0,Guis,InvKss,1.0);
    Cuiui.Multiply('N','N',1.0,GuisInvKss,Gsui,1.0);

  }

#endif

//   std::cout << elemat1;
//   std::cout << elevec1;

//   elemat1_epetra.Print( std::cout );
}

/*--------------------------------------------------------------------------------
 *--------------------------------------------------------------------------------*/
template <>
void DRT::ELEMENTS::FluidEleCalcXFEM<DRT::Element::tri3>::ElementXfemInterface(
  DRT::ELEMENTS::Fluid3 * ele,
  DRT::Discretization & dis,
  const std::vector<int> & lm,
  const DRT::UTILS::GaussIntegration & intpoints,
  DRT::Discretization & cutdis,
  const std::map<int, std::vector<GEO::CUT::BoundaryCell*> > & bcells,
  const std::map<int, std::vector<DRT::UTILS::GaussIntegration> > & bintpoints,
  std::map<int, std::vector<Epetra_SerialDenseMatrix> > & side_coupling,
  Teuchos::ParameterList&    params,
  Epetra_SerialDenseMatrix&  elemat1_epetra,
  Epetra_SerialDenseVector&  elevec1_epetra,
  Epetra_SerialDenseMatrix&  Cuiui
  )
{
  dserror( "distype not supported" );
}

/*--------------------------------------------------------------------------------
 *--------------------------------------------------------------------------------*/
template <>
void DRT::ELEMENTS::FluidEleCalcXFEM<DRT::Element::tri6>::ElementXfemInterface(
  DRT::ELEMENTS::Fluid3 * ele,
  DRT::Discretization & dis,
  const std::vector<int> & lm,
  const DRT::UTILS::GaussIntegration & intpoints,
  DRT::Discretization & cutdis,
  const std::map<int, std::vector<GEO::CUT::BoundaryCell*> > & bcells,
  const std::map<int, std::vector<DRT::UTILS::GaussIntegration> > & bintpoints,
  std::map<int, std::vector<Epetra_SerialDenseMatrix> > & side_coupling,
  Teuchos::ParameterList&    params,
  Epetra_SerialDenseMatrix&  elemat1_epetra,
  Epetra_SerialDenseVector&  elevec1_epetra,
  Epetra_SerialDenseMatrix&  Cuiui
  )
{
  dserror( "distype not supported" );
}

/*--------------------------------------------------------------------------------
 *--------------------------------------------------------------------------------*/
template <>
void DRT::ELEMENTS::FluidEleCalcXFEM<DRT::Element::quad4>::ElementXfemInterface(
  DRT::ELEMENTS::Fluid3 * ele,
  DRT::Discretization & dis,
  const std::vector<int> & lm,
  const DRT::UTILS::GaussIntegration & intpoints,
  DRT::Discretization & cutdis,
  const std::map<int, std::vector<GEO::CUT::BoundaryCell*> > & bcells,
  const std::map<int, std::vector<DRT::UTILS::GaussIntegration> > & bintpoints,
  std::map<int, std::vector<Epetra_SerialDenseMatrix> > & side_coupling,
  Teuchos::ParameterList&    params,
  Epetra_SerialDenseMatrix&  elemat1_epetra,
  Epetra_SerialDenseVector&  elevec1_epetra,
  Epetra_SerialDenseMatrix&  Cuiui
  )
{
  dserror( "distype not supported" );
}

/*--------------------------------------------------------------------------------
 *--------------------------------------------------------------------------------*/
template <>
void DRT::ELEMENTS::FluidEleCalcXFEM<DRT::Element::quad8>::ElementXfemInterface(
  DRT::ELEMENTS::Fluid3 * ele,
  DRT::Discretization & dis,
  const std::vector<int> & lm,
  const DRT::UTILS::GaussIntegration & intpoints,
  DRT::Discretization & cutdis,
  const std::map<int, std::vector<GEO::CUT::BoundaryCell*> > & bcells,
  const std::map<int, std::vector<DRT::UTILS::GaussIntegration> > & bintpoints,
  std::map<int, std::vector<Epetra_SerialDenseMatrix> > & side_coupling,
  Teuchos::ParameterList&    params,
  Epetra_SerialDenseMatrix&  elemat1_epetra,
  Epetra_SerialDenseVector&  elevec1_epetra,
  Epetra_SerialDenseMatrix&  Cuiui
  )
{
  dserror( "distype not supported" );
}

/*--------------------------------------------------------------------------------
 *--------------------------------------------------------------------------------*/
template <>
void DRT::ELEMENTS::FluidEleCalcXFEM<DRT::Element::quad9>::ElementXfemInterface(
  DRT::ELEMENTS::Fluid3 * ele,
  DRT::Discretization & dis,
  const std::vector<int> & lm,
  const DRT::UTILS::GaussIntegration & intpoints,
  DRT::Discretization & cutdis,
  const std::map<int, std::vector<GEO::CUT::BoundaryCell*> > & bcells,
  const std::map<int, std::vector<DRT::UTILS::GaussIntegration> > & bintpoints,
  std::map<int, std::vector<Epetra_SerialDenseMatrix> > & side_coupling,
  Teuchos::ParameterList&    params,
  Epetra_SerialDenseMatrix&  elemat1_epetra,
  Epetra_SerialDenseVector&  elevec1_epetra,
  Epetra_SerialDenseMatrix&  Cuiui
  )
{
  dserror( "distype not supported" );
}

/*--------------------------------------------------------------------------------
 *--------------------------------------------------------------------------------*/
template <>
void DRT::ELEMENTS::FluidEleCalcXFEM<DRT::Element::nurbs9>::ElementXfemInterface(
  DRT::ELEMENTS::Fluid3 * ele,
  DRT::Discretization & dis,
  const std::vector<int> & lm,
  const DRT::UTILS::GaussIntegration & intpoints,
  DRT::Discretization & cutdis,
  const std::map<int, std::vector<GEO::CUT::BoundaryCell*> > & bcells,
  const std::map<int, std::vector<DRT::UTILS::GaussIntegration> > & bintpoints,
  std::map<int, std::vector<Epetra_SerialDenseMatrix> > & side_coupling,
  Teuchos::ParameterList&    params,
  Epetra_SerialDenseMatrix&  elemat1_epetra,
  Epetra_SerialDenseVector&  elevec1_epetra,
  Epetra_SerialDenseMatrix&  Cuiui
  )
{
  dserror( "distype not supported" );
}



/*--------------------------------------------------------------------------------
 *--------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidEleCalcXFEM<distype>::ElementXfemInterfaceNitsche(
  DRT::ELEMENTS::Fluid3 * ele,
  DRT::Discretization & dis,
  const std::vector<int> & lm,
  const DRT::UTILS::GaussIntegration & intpoints,
  DRT::Discretization & cutdis,
  const std::map<int, std::vector<GEO::CUT::BoundaryCell*> > & bcells,
  const std::map<int, std::vector<DRT::UTILS::GaussIntegration> > & bintpoints,
  std::map<int, std::vector<Epetra_SerialDenseMatrix> > & side_coupling,
  Teuchos::ParameterList&    params,
  Epetra_SerialDenseMatrix&  elemat1_epetra,
  Epetra_SerialDenseVector&  elevec1_epetra,
  Epetra_SerialDenseMatrix&  Cuiui
  )
{

  const Teuchos::RCP<Epetra_Vector> iforcecol = params.get<Teuchos::RCP<Epetra_Vector> >("iforcenp", Teuchos::null);

  double nitsche_stab      = params.get<double>("nitsche_stab");
  double nitsche_stab_conv = params.get<double>("nitsche_stab_conv");

  // volume integral, just for the partial volume
  double meas_partial_volume = 0.0;

  int eid = ele->Id();

  for ( DRT::UTILS::GaussIntegration::iterator iquad=intpoints.begin(); iquad!=intpoints.end(); ++iquad )
  {
    // evaluate shape functions and derivatives at integration point
    my::EvalShapeFuncAndDerivsAtIntPoint(iquad,eid);

    //----------------------------------------------------------------------
    // set time-integration factors for left- and right-hand side
    // (two right-hand-side factors: general and for residuals)
    //----------------------------------------------------------------------
    meas_partial_volume += my::fac_;
  }


  //----------------------------------------------------------------------------
  //                         ELEMENT GEOMETRY
  //----------------------------------------------------------------------------

  // get node coordinates
  GEO::fillInitialPositionArray< distype, my::nsd_, LINALG::Matrix<my::nsd_,my::nen_> >( ele, my::xyze_ );

  LINALG::Matrix<my::nsd_,my::nen_> evelaf(true);
  LINALG::Matrix<my::nen_,1> epreaf(true);
  ExtractValuesFromGlobalVector(dis, lm, *my::rotsymmpbc_, &evelaf, &epreaf, "velaf");

  // integrate surface

  DRT::Element::LocationArray cutla( 1 );

  LINALG::Matrix<3,1> normal;
  LINALG::Matrix<3,1> x_side;

  bool fluidfluidcoupling = false;

  std::map<int, Teuchos::RCP<XFLUID::SideInterface<distype> > > side_impl;
  Teuchos::RCP<XFLUID::SideInterface<distype> > si;

  // find all the intersecting elements of actele
  std::set<int> begids;
  for (std::map<int,  std::vector<GEO::CUT::BoundaryCell*> >::const_iterator bc=bcells.begin();
      bc!=bcells.end(); ++bc )
  {
    int sid = bc->first;
    begids.insert(sid);
  }

  // map of boundary element gids, to coupling matrices Gsui and Guis (Cuiui = - Guis*Kss^-1*Gsui)
  std::map<int, std::vector<Epetra_SerialDenseMatrix> > Cuiui_coupling;

  // lm vector of all intersecting elements (boundary elements which intersect the current element)
  std::vector<int> patchelementslmv;
  std::vector<int> patchelementslmowner;
  for (std::set<int>::const_iterator bgid=begids.begin(); bgid!=begids.end(); ++bgid)
  {
    DRT::Element * side = cutdis.gElement(*bgid); // for each boundary element there is one corresponding side
    vector<int> patchlm;
    vector<int> patchlmowner;
    vector<int> patchlmstride;
    side->LocationVector(cutdis, patchlm, patchlmowner, patchlmstride);

    patchelementslmv.reserve( patchelementslmv.size() + patchlm.size());
    patchelementslmv.insert(patchelementslmv.end(), patchlm.begin(), patchlm.end());

    patchelementslmowner.reserve( patchelementslmowner.size() + patchlmowner.size());
    patchelementslmowner.insert( patchelementslmowner.end(), patchlmowner.begin(), patchlmowner.end());

    // get coupling matrices for the current side (boundary element)
    std::vector<Epetra_SerialDenseMatrix> & Cuiui_matrices = Cuiui_coupling[*bgid];

    Cuiui_matrices.resize(1);
    Cuiui_matrices[0].Reshape(patchlm.size(),patchlm.size()); //Cuiui (coupling between background elements sigma and current side!)
  }


   // map of side-element id and Guass points
   for ( std::map<int, std::vector<DRT::UTILS::GaussIntegration> >::const_iterator i=bintpoints.begin();
       i!=bintpoints.end();
       ++i )
   {
     int sid = i->first;
     const std::vector<DRT::UTILS::GaussIntegration> & cutintpoints = i->second;

     std::map<int, std::vector<GEO::CUT::BoundaryCell*> >::const_iterator j = bcells.find( sid );
     if ( j==bcells.end() )
       dserror( "missing boundary cell" );

     const std::vector<GEO::CUT::BoundaryCell*> & bcs = j->second;
     if ( bcs.size()!=cutintpoints.size() )
       dserror( "boundary cell integration rules mismatch" );

     DRT::Element * side = cutdis.gElement( sid );
     side->LocationVector(cutdis,cutla,false);

     const int numnodes = side->NumNode();
     DRT::Node ** nodes = side->Nodes();
     Epetra_SerialDenseMatrix side_xyze( 3, numnodes );
     for ( int i=0; i<numnodes; ++i )
     {
       const double * x = nodes[i]->X();
       std::copy( x, x+3, &side_xyze( 0, i ) );
     }

     std::map<int,std::vector<Epetra_SerialDenseMatrix> >::iterator c = side_coupling.find( sid );

     std::vector<Epetra_SerialDenseMatrix> & side_matrices = c->second;

     if ( side_matrices.size()==3 )
       fluidfluidcoupling = true;


     if(fluidfluidcoupling)
     {
       // coupling matrices between background element and one! side
       Epetra_SerialDenseMatrix & C_uiu  = side_matrices[0];
       Epetra_SerialDenseMatrix & C_uui  = side_matrices[1];
       Epetra_SerialDenseMatrix & rhC_ui = side_matrices[2];

       // coupling matrices between one side and itself via the element Kss
       std::map<int,std::vector<Epetra_SerialDenseMatrix> >::iterator c2 = Cuiui_coupling.find( sid );
       std::vector<Epetra_SerialDenseMatrix> & Cuiui_matrices = c2->second;
       Epetra_SerialDenseMatrix & eleCuiui = Cuiui_matrices[0];

       si = XFLUID::SideInterface<distype>::Impl(side,C_uiu,C_uui,rhC_ui,eleCuiui,side_xyze);
     }
     else
     {
       si = XFLUID::SideInterface<distype>::Impl(side,side_xyze);
     }

     side_impl[sid] = si;

     // get velocity at integration point of boundary dis

     si->eivel(cutdis,"ivelnp",cutla[0].lm_);
     si->addeidisp(cutdis,"idispnp",cutla[0].lm_,side_xyze);


     double meas_surface = 0.0;

     // pre-evaluate for element-size
     // loop gausspoints
     for ( std::vector<DRT::UTILS::GaussIntegration>::const_iterator i=cutintpoints.begin();
         i!=cutintpoints.end();
         ++i )
     {
       const DRT::UTILS::GaussIntegration & gi = *i;
       GEO::CUT::BoundaryCell * bc = bcs[i - cutintpoints.begin()]; // get the corresponding boundary cell

       //gi.Print();

       // TODO: do this transformation not twice!!!
       for ( DRT::UTILS::GaussIntegration::iterator iquad=gi.begin(); iquad!=gi.end(); ++iquad )
       {
         double drs = 0.0; // transformation factor between reference cell and linearized boundary cell
         LINALG::Matrix<3,1> x_gp_lin(true); // gp in xyz-system on linearized interface
         if(bc->Shape()==DRT::Element::tri3 || bc->Shape()==DRT::Element::quad4)
         {
#ifdef BOUNDARYCELL_TRANSFORMATION_OLD
           const LINALG::Matrix<2,1> eta( iquad.Point() ); // xi-coordinates with respect to side

           double drs = 0;

           si->Evaluate(eta,x_side,normal,drs);

           const double fac = drs*iquad.Weight();

           // find element local position of gauss point at interface
           GEO::CUT::Position<distype> pos( my::xyze_, x_side );
           pos.Compute();
           const LINALG::Matrix<3,1> & rst = pos.LocalCoordinates();

  #else
           const LINALG::Matrix<2,1> eta( iquad.Point() ); // eta-coordinates with respect to cell

           normal.Clear();

           // get normal vector on linearized boundary cell, x-coordinates of gaussian point and surface transformation factor
           switch ( bc->Shape() )
           {
           case DRT::Element::tri3:
           {
             bc->Transform<DRT::Element::tri3>(eta, x_gp_lin, normal, drs);
             break;
           }
           case DRT::Element::quad4:
           {
             bc->Transform<DRT::Element::quad4>(eta, x_gp_lin, normal, drs);
             break;
           }
           default:
             throw std::runtime_error( "unsupported integration cell type" );
           }
         }
         else if(bc->Shape()==DRT::Element::dis_none)
         {
           drs = 1.0;
           normal = bc->GetNormalVector();
           const double* gpcord = iquad.Point();
           for (int idim=0;idim<3;idim++)
           {
              x_gp_lin(idim,0) = gpcord[idim];
           }
         }
#endif


         meas_surface += drs*iquad.Weight();
       }
     }



      //-----------------------------------------------------------------------------------

      double stabfac = 0.0;         // Nitsche stabilization factor
      double stabfac_conv = 0.0;    // Nitsche convective stabilization factor


      // define stabilization parameters and mortaring weights

      double kappa1 = 1.0;      // Xfluid-sided mortaring

      if(meas_partial_volume < 0.0) dserror(" measure of cut partial volume is smaller than 0.0: %f Attention with increasing Nitsche-Parameter!!!", meas_partial_volume);




      if( kappa1 > 1.0 || kappa1 < 0.0) dserror("Nitsche weights for inverse estimate kappa1 lies not in [0,1]: %d", kappa1);

      double kappa2 = 1.0-kappa1;

      // loop gausspoints
      for ( std::vector<DRT::UTILS::GaussIntegration>::const_iterator i=cutintpoints.begin();
            i!=cutintpoints.end();
            ++i )
      {
        const DRT::UTILS::GaussIntegration & gi = *i;
        GEO::CUT::BoundaryCell * bc = bcs[i - cutintpoints.begin()]; // get the corresponding boundary cell

        //gi.Print();

        for ( DRT::UTILS::GaussIntegration::iterator iquad=gi.begin(); iquad!=gi.end(); ++iquad )
        {
          double drs = 0; // transformation factor between reference cell and linearized boundary cell
          LINALG::Matrix<3,1> x_gp_lin(true); // gp in xyz-system on linearized interface
          if(bc->Shape()==DRT::Element::tri3 || bc->Shape()==DRT::Element::quad4)
          {
#ifdef BOUNDARYCELL_TRANSFORMATION_OLD
            const LINALG::Matrix<2,1> eta( iquad.Point() ); // xi-coordinates with respect to side

            double drs = 0;

            si->Evaluate(eta,x_side,normal,drs);

            const double fac = drs*iquad.Weight();

            // find element local position of gauss point at interface
            GEO::CUT::Position<distype> pos( my::xyze_, x_side );
            pos.Compute();
            const LINALG::Matrix<3,1> & rst = pos.LocalCoordinates();

    #else
            const LINALG::Matrix<2,1> eta( iquad.Point() ); // eta-coordinates with respect to cell

            normal.Clear();

            // get normal vector on linearized boundary cell, x-coordinates of gaussian point and surface transformation factor
            switch ( bc->Shape() )
            {
            case DRT::Element::tri3:
            {
                bc->Transform<DRT::Element::tri3>(eta, x_gp_lin, normal, drs);
              break;
            }
            case DRT::Element::quad4:
            {
                bc->Transform<DRT::Element::quad4>(eta, x_gp_lin, normal, drs);
              break;
            }
            default:
              throw std::runtime_error( "unsupported integration cell type" );
            }
          }
          else if(bc->Shape()==DRT::Element::dis_none)
          {
            drs = 1.0;
            normal = bc->GetNormalVector();
            const double* gpcord = iquad.Point();
            for (int idim=0;idim<3;idim++)
            {
               x_gp_lin(idim,0) = gpcord[idim];
            }
          }


        const double fac = drs*iquad.Weight();

        // find element local position of gauss point
        GEO::CUT::Position<distype> pos( my::xyze_, x_gp_lin );
        pos.Compute();
        const LINALG::Matrix<3,1> & rst = pos.LocalCoordinates();


        // project gaussian point from linearized interface to warped side (get local side coordinates)
        LINALG::Matrix<2,1> xi_side(true);
        si->ProjectOnSide(x_gp_lin, x_side, xi_side);

#endif



        // evaluate shape functions
        DRT::UTILS::shape_function<distype>( rst, my::funct_ );

        // evaluate shape functions and derivatives at integration point

        DRT::UTILS::shape_function_deriv1<distype>(rst,my::deriv_);
        my::xjm_.MultiplyNT(my::deriv_,my::xyze_);
        my::det_ = my::xji_.Invert(my::xjm_);


        // compute global first derivates
        my::derxy_.Multiply(my::xji_,my::deriv_);


        // get velocity at integration point
        // (values at n+alpha_F for generalized-alpha scheme, n+1 otherwise)
        my::velint_.Multiply(evelaf,my::funct_);


        // get velocity derivatives at integration point
        // (values at n+alpha_F for generalized-alpha scheme, n+1 otherwise)
        my::vderxy_.MultiplyNT(evelaf,my::derxy_);

        // get pressure at integration point
        // (value at n+alpha_F for generalized-alpha scheme, n+1 otherwise)
        double press = my::funct_.Dot(epreaf);

        const double timefacfac = my::f3Parameter_->timefac_ * fac;

        double h_k= 0.4/28.0;

        //      stabfac      = nitsche_stab * visceff_ * meas_surface / meas_partial_volume;
        //      stabfac_conv      = nitsche_stab_conv * meas_surface / meas_partial_volume;
        stabfac = nitsche_stab * my::visceff_ / h_k;
        stabfac_conv = nitsche_stab_conv / h_k;

        //=================================================================================
        // definition in Burman 2007
        // Interior penalty variational multiscale method for the incompressible Navier-Stokes equation:
        // Monitoring artificial dissipation
        /*
        //      viscous_Nitsche-part, convective inflow part
        //                |                 |
        //                    mu                                  /              \
        //  max( gamma_Nit * ----  , | u_h * n | )       *       |  u_h - u, v_h  |
        //                    h_k                                 \              /
        */

        stabfac = max(nitsche_stab*my::visceff_/h_k, fabs(my::velint_.Dot(normal)));
        stabfac_conv = nitsche_stab_conv * max(1.0, max(fabs(my::velint_.Dot(normal)) , my::visceff_ / h_k) );
        //=================================================================================


      if(fluidfluidcoupling)
      {
        bool bg_mortaring = true; // one-sided background fluid mortaring (kappa1=1, kappa2=0)

        //zero velocity jump for fluidfluidcoupling
        LINALG::Matrix<my::nsd_,1> ivelint_WDBC_JUMP(true);


        si->buildCouplingMatricesNitsche( elemat1_epetra,          // standard bg-bg-matrix
                                          elevec1_epetra,          // standard bg-rhs
                                          fluidfluidcoupling,      // assemble coupling terms (yes/no)
                                          bg_mortaring,            // yes: background-sided mortaring, no: coupling between two meshes (mixed mortaring)
                                          normal,                  // normal vector
                                          timefacfac,              // theta*dt
                                          my::visceff_,                // viscosity in background fluid
                                          my::visceff_,                // viscosity in embedded fluid
                                          kappa1,                  // mortaring weighting
                                          kappa2,                  // mortaring weighting
                                          stabfac,                 // Nitsche non-dimensionless stabilization factor
                                          stabfac_conv,            // Nitsche convective non-dimensionless stabilization factor
                                          my::funct_,                  // bg shape functions
                                          my::derxy_,                  // bg deriv
                                          my::vderxy_,                 // bg deriv^n
                                          press,                   // bg p^n
                                          my::velint_,                 // bg u^n
                                          ivelint_WDBC_JUMP         // Dirichlet velocity vector or prescribed jump vector
                                          );



      }
      if(!fluidfluidcoupling)
      {
        // case for one-sided weak Dirichlet
        bool bg_mortaring = true; // one-sided background fluid mortaring (kappa1=1, kappa2=0)

        // prescribed velocity vector at weak Dirichlet boundary
        LINALG::Matrix<my::nsd_,1> ivelint_WDBC_JUMP(true);
        si->get_vel_WeakDBC(ivelint_WDBC_JUMP);

        si->buildCouplingMatricesNitsche( elemat1_epetra,          // standard bg-bg-matrix
                                          elevec1_epetra,          // standard bg-rhs
                                          fluidfluidcoupling,      // assemble coupling terms (yes/no)
                                          bg_mortaring,            // yes: background-sided mortaring, no: coupling between two meshes (mixed mortaring)
                                          normal,                  // normal vector
                                          timefacfac,              // theta*dt
                                          my::visceff_,                // viscosity in background fluid
                                          0.0,                     // viscosity in embedded fluid
                                          kappa1,                  // mortaring weighting
                                          kappa2,                  // mortaring weighting
                                          stabfac,                 // Nitsche non-dimensionless stabilization factor
                                          stabfac_conv,            // Nitsche convective non-dimensionless stabilization factor
                                          my::funct_,                  // bg shape functions
                                          my::derxy_,                  // bg deriv
                                          my::vderxy_,                 // bg deriv^n
                                          press,                   // bg p^n
                                          my::velint_,                  // bg u^n
                                          ivelint_WDBC_JUMP         // Dirichlet velocity vector or prescribed jump vector
                                          );
      }

      // calculate interface forces
      if(!fluidfluidcoupling)
      {


        // get pressure at integration point
        // (value at n+alpha_F for generalized-alpha scheme, n+1 otherwise)
        double press = my::funct_.Dot(epreaf);


        // compute the stresses at the current Gaussian point for computing the interface force
        LINALG::Matrix<my::nsd_,my::nsd_> eps(true);
        for(int i=0; i<my::nsd_; i++)
        {
          for(int j=0; j<my::nsd_; j++)
          {
            eps(i,j) = 0.5 * (my::vderxy_(i,j) + my::vderxy_(j,i));
          }
        }


        // t = ( -pI + 2mu eps(u) )*n^f
        LINALG::Matrix<my::nsd_,1> traction (true);
        traction.Multiply(eps, normal);
        traction.Scale(2.0*my::visceff_);

        // add the pressure part
        traction.Update( -press, normal, 1.0);

        // we need normal vector on fluid
        traction.Scale(-1.0);

        double surf_fac = drs*iquad.Weight();

        si->buildInterfaceForce(iforcecol, cutdis, cutla[0].lm_, traction, surf_fac );

      }

      }
    }
  }

  if ( fluidfluidcoupling )
  {
    // build Gsui and Guis
    int ipatchsizesbefore = 0;
    for (std::map<int, std::vector<Epetra_SerialDenseMatrix> >::const_iterator m=Cuiui_coupling.begin();
        m!=Cuiui_coupling.end(); ++m)
    {

      int bid = m->first;
      std::vector<Epetra_SerialDenseMatrix> & Cuiui_mats = Cuiui_coupling[bid];

      //Cuiui matrices in Cuiui_mats[0]

      // assemble Cuiui
      for ( int ic=0; ic<Cuiui_mats[0].N(); ++ic) // Cuiui includes only ui,ui coupling, not (ui,p) ...
      {
        for ( int ir=0; ir<Cuiui_mats[0].M(); ++ir )
        {
          Cuiui(ir+ipatchsizesbefore,ic+ipatchsizesbefore) = Cuiui_mats[0](ir,ic);
        }
      }

      ipatchsizesbefore += Cuiui_mats[0].N();

    }
  }
}

/*--------------------------------------------------------------------------------
 *--------------------------------------------------------------------------------*/
template <>
void DRT::ELEMENTS::FluidEleCalcXFEM<DRT::Element::tri3>::ElementXfemInterfaceNitsche(
  DRT::ELEMENTS::Fluid3 * ele,
  DRT::Discretization & dis,
  const std::vector<int> & lm,
  const DRT::UTILS::GaussIntegration & intpoints,
  DRT::Discretization & cutdis,
  const std::map<int, std::vector<GEO::CUT::BoundaryCell*> > & bcells,
  const std::map<int, std::vector<DRT::UTILS::GaussIntegration> > & bintpoints,
  std::map<int, std::vector<Epetra_SerialDenseMatrix> > & side_coupling,
  Teuchos::ParameterList&    params,
  Epetra_SerialDenseMatrix&  elemat1_epetra,
  Epetra_SerialDenseVector&  elevec1_epetra,
  Epetra_SerialDenseMatrix&  Cuiui
  )
{
  dserror( "distype not supported" );
}

/*--------------------------------------------------------------------------------
 *--------------------------------------------------------------------------------*/
template <>
void DRT::ELEMENTS::FluidEleCalcXFEM<DRT::Element::tri6>::ElementXfemInterfaceNitsche(
  DRT::ELEMENTS::Fluid3 * ele,
  DRT::Discretization & dis,
  const std::vector<int> & lm,
  const DRT::UTILS::GaussIntegration & intpoints,
  DRT::Discretization & cutdis,
  const std::map<int, std::vector<GEO::CUT::BoundaryCell*> > & bcells,
  const std::map<int, std::vector<DRT::UTILS::GaussIntegration> > & bintpoints,
  std::map<int, std::vector<Epetra_SerialDenseMatrix> > & side_coupling,
  Teuchos::ParameterList&    params,
  Epetra_SerialDenseMatrix&  elemat1_epetra,
  Epetra_SerialDenseVector&  elevec1_epetra,
  Epetra_SerialDenseMatrix&  Cuiui
  )
{
  dserror( "distype not supported" );
}

/*--------------------------------------------------------------------------------
 *--------------------------------------------------------------------------------*/
template <>
void DRT::ELEMENTS::FluidEleCalcXFEM<DRT::Element::quad4>::ElementXfemInterfaceNitsche(
  DRT::ELEMENTS::Fluid3 * ele,
  DRT::Discretization & dis,
  const std::vector<int> & lm,
  const DRT::UTILS::GaussIntegration & intpoints,
  DRT::Discretization & cutdis,
  const std::map<int, std::vector<GEO::CUT::BoundaryCell*> > & bcells,
  const std::map<int, std::vector<DRT::UTILS::GaussIntegration> > & bintpoints,
  std::map<int, std::vector<Epetra_SerialDenseMatrix> > & side_coupling,
  Teuchos::ParameterList&    params,
  Epetra_SerialDenseMatrix&  elemat1_epetra,
  Epetra_SerialDenseVector&  elevec1_epetra,
  Epetra_SerialDenseMatrix&  Cuiui
  )
{
  dserror( "distype not supported" );
}

/*--------------------------------------------------------------------------------
 *--------------------------------------------------------------------------------*/
template <>
void DRT::ELEMENTS::FluidEleCalcXFEM<DRT::Element::quad8>::ElementXfemInterfaceNitsche(
  DRT::ELEMENTS::Fluid3 * ele,
  DRT::Discretization & dis,
  const std::vector<int> & lm,
  const DRT::UTILS::GaussIntegration & intpoints,
  DRT::Discretization & cutdis,
  const std::map<int, std::vector<GEO::CUT::BoundaryCell*> > & bcells,
  const std::map<int, std::vector<DRT::UTILS::GaussIntegration> > & bintpoints,
  std::map<int, std::vector<Epetra_SerialDenseMatrix> > & side_coupling,
  Teuchos::ParameterList&    params,
  Epetra_SerialDenseMatrix&  elemat1_epetra,
  Epetra_SerialDenseVector&  elevec1_epetra,
  Epetra_SerialDenseMatrix&  Cuiui
  )
{
  dserror( "distype not supported" );
}

/*--------------------------------------------------------------------------------
 *--------------------------------------------------------------------------------*/
template <>
void DRT::ELEMENTS::FluidEleCalcXFEM<DRT::Element::quad9>::ElementXfemInterfaceNitsche(
  DRT::ELEMENTS::Fluid3 * ele,
  DRT::Discretization & dis,
  const std::vector<int> & lm,
  const DRT::UTILS::GaussIntegration & intpoints,
  DRT::Discretization & cutdis,
  const std::map<int, std::vector<GEO::CUT::BoundaryCell*> > & bcells,
  const std::map<int, std::vector<DRT::UTILS::GaussIntegration> > & bintpoints,
  std::map<int, std::vector<Epetra_SerialDenseMatrix> > & side_coupling,
  Teuchos::ParameterList&    params,
  Epetra_SerialDenseMatrix&  elemat1_epetra,
  Epetra_SerialDenseVector&  elevec1_epetra,
  Epetra_SerialDenseMatrix&  Cuiui
  )
{
  dserror( "distype not supported" );
}

/*--------------------------------------------------------------------------------
 *--------------------------------------------------------------------------------*/
template <>
void DRT::ELEMENTS::FluidEleCalcXFEM<DRT::Element::nurbs9>::ElementXfemInterfaceNitsche(
  DRT::ELEMENTS::Fluid3 * ele,
  DRT::Discretization & dis,
  const std::vector<int> & lm,
  const DRT::UTILS::GaussIntegration & intpoints,
  DRT::Discretization & cutdis,
  const std::map<int, std::vector<GEO::CUT::BoundaryCell*> > & bcells,
  const std::map<int, std::vector<DRT::UTILS::GaussIntegration> > & bintpoints,
  std::map<int, std::vector<Epetra_SerialDenseMatrix> > & side_coupling,
  Teuchos::ParameterList&    params,
  Epetra_SerialDenseMatrix&  elemat1_epetra,
  Epetra_SerialDenseVector&  elevec1_epetra,
  Epetra_SerialDenseMatrix&  Cuiui
  )
{
  dserror( "distype not supported" );
}



/*--------------------------------------------------------------------------------
 *--------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidEleCalcXFEM<distype>::ElementXfemInterfaceNitscheTwoSided(
  DRT::ELEMENTS::Fluid3 * ele,
  DRT::Discretization & dis,
  const std::vector<int> & lm,
  const DRT::UTILS::GaussIntegration & intpoints,
  DRT::Discretization & cutdis,
  const std::map<int, std::vector<GEO::CUT::BoundaryCell*> > & bcells,
  const std::map<int, std::vector<DRT::UTILS::GaussIntegration> > & bintpoints,
  std::map<int, std::vector<Epetra_SerialDenseMatrix> > & side_coupling,
  Teuchos::ParameterList&    params,
  DRT::Discretization &      alediscret,
  map<int,int> &             boundary_emb_gid_map,
  Epetra_SerialDenseMatrix&  elemat1_epetra,
  Epetra_SerialDenseVector&  elevec1_epetra,
  Epetra_SerialDenseMatrix&  Cuiui
  )
{

  INPAR::XFEM::CouplingStrategy coupling_strategy = params.get<INPAR::XFEM::CouplingStrategy>("coupling_strategy");

  double nitsche_stab      = params.get<double>("nitsche_stab");
  double nitsche_stab_conv = params.get<double>("nitsche_stab_conv");

  // volume integral, just for the partial volume
  double meas_partial_volume = 0.0;

  int eid = ele->Id();

  // for two-sided mortaring the stabilization parameter depends on interface/volume fraction
  if(coupling_strategy == INPAR::XFEM::Two_Sided_Mortaring)
  {
    for ( DRT::UTILS::GaussIntegration::iterator iquad=intpoints.begin(); iquad!=intpoints.end(); ++iquad )
    {
      // evaluate shape functions and derivatives at integration point
      my::EvalShapeFuncAndDerivsAtIntPoint(iquad,eid);

      meas_partial_volume += my::fac_;
    }
  }


  //----------------------------------------------------------------------------
  //                         ELEMENT GEOMETRY
  //----------------------------------------------------------------------------

  // get node coordinates
  GEO::fillInitialPositionArray< distype, my::nsd_, LINALG::Matrix<my::nsd_,my::nen_> >( ele, my::xyze_ );

  LINALG::Matrix<my::nsd_,my::nen_> evelaf(true);
  LINALG::Matrix<my::nen_,1> epreaf(true);
  ExtractValuesFromGlobalVector(dis, lm, *my::rotsymmpbc_, &evelaf, &epreaf, "velaf");


  // numbering for velocity components and pressure
  //const unsigned Velx = 0;
  //const unsigned Vely = 1;
  //const unsigned Velz = 2;
  //const unsigned Pres = 3;


  // integrate surface

   DRT::Element::LocationArray alela( 1 );
   DRT::Element::LocationArray cutla( 1 );

   LINALG::Matrix<3,1> normal;
   LINALG::Matrix<3,1> x_side;

   bool fluidfluidcoupling = false;

   std::map<int, Teuchos::RCP<XFLUID::EmbCoupling<distype> > > emb_impl;
   std::map<int, Teuchos::RCP<XFLUID::SideInterface<distype> > > side_impl;
   Teuchos::RCP<XFLUID::SideInterface<distype> > si;
   Teuchos::RCP<XFLUID::EmbCoupling<distype> > emb;

   // find all the intersecting elements of actele
   std::set<int> begids;
   for (std::map<int,  std::vector<GEO::CUT::BoundaryCell*> >::const_iterator bc=bcells.begin();
        bc!=bcells.end(); ++bc )
   {
     int sid = bc->first;
     begids.insert(sid);
   }

   // map of boundary element gids, to coupling matrices Cuiui
   std::map<int, std::vector<Epetra_SerialDenseMatrix> > Cuiui_coupling;

   // lm vector of all intersecting elements (boundary elements which intersect the current element)
   std::vector<int> patchelementslmv;
   std::vector<int> patchelementslmowner;
   for (std::set<int>::const_iterator bgid=begids.begin(); bgid!=begids.end(); ++bgid)
   {
     int emb_gid = boundary_emb_gid_map.find(*bgid)->second;
     DRT::Element * emb_ele = alediscret.gElement(emb_gid);

     vector<int> patchlm;
     vector<int> patchlmowner;
     vector<int> patchlmstride;
     emb_ele->LocationVector(alediscret, patchlm, patchlmowner, patchlmstride);

     patchelementslmv.reserve( patchelementslmv.size() + patchlm.size());
     patchelementslmv.insert(patchelementslmv.end(), patchlm.begin(), patchlm.end());

     patchelementslmowner.reserve( patchelementslmowner.size() + patchlmowner.size());
     patchelementslmowner.insert( patchelementslmowner.end(), patchlmowner.begin(), patchlmowner.end());

     // get coupling matrices for the current side (boundary element)
     std::vector<Epetra_SerialDenseMatrix> & Cuiui_matrices = Cuiui_coupling[*bgid];

     Cuiui_matrices.resize(1);
     Cuiui_matrices[0].Reshape(patchlm.size(),patchlm.size()); //Cuiui (coupling between background elements sigma and current side!)
   }


    // map of side-element id and Guass points
    for ( std::map<int, std::vector<DRT::UTILS::GaussIntegration> >::const_iterator i=bintpoints.begin();
          i!=bintpoints.end();
          ++i )
    {
      int sid = i->first;
      const std::vector<DRT::UTILS::GaussIntegration> & cutintpoints = i->second;

      std::map<int, std::vector<GEO::CUT::BoundaryCell*> >::const_iterator j = bcells.find( sid );
      if ( j==bcells.end() )
        dserror( "missing boundary cell" );

      const std::vector<GEO::CUT::BoundaryCell*> & bcs = j->second;
      if ( bcs.size()!=cutintpoints.size() )
        dserror( "boundary cell integration rules mismatch" );

      DRT::Element * side = cutdis.gElement( sid );
      side->LocationVector(cutdis,cutla,false);
      DRT::Element * emb_ele = alediscret.gElement( boundary_emb_gid_map.find(sid)->second );
      emb_ele->LocationVector(alediscret,alela,false);

      const int emb_numnodes = emb_ele->NumNode();
      DRT::Node ** emb_nodes = emb_ele->Nodes();
      Epetra_SerialDenseMatrix emb_xyze( 3, emb_numnodes );
      for ( int i=0; i<emb_numnodes; ++i )
      {
        const double * x = emb_nodes[i]->X();
        std::copy( x, x+3, &emb_xyze( 0, i ) );
      }

      const int numnodes = side->NumNode();
      DRT::Node ** nodes = side->Nodes();
      Epetra_SerialDenseMatrix side_xyze( 3, numnodes );
      for ( int i=0; i<numnodes; ++i )
      {
        const double * x = nodes[i]->X();
        std::copy( x, x+3, &side_xyze( 0, i ) );
      }

      std::map<int,std::vector<Epetra_SerialDenseMatrix> >::iterator c = side_coupling.find( sid );

      std::vector<Epetra_SerialDenseMatrix> & side_matrices = c->second;

      if ( side_matrices.size()==3 )
        fluidfluidcoupling = true;


      if(fluidfluidcoupling)
      {
        // coupling matrices between background element and one! side
          Epetra_SerialDenseMatrix & C_uiu  = side_matrices[0];
          Epetra_SerialDenseMatrix & C_uui  = side_matrices[1];
          Epetra_SerialDenseMatrix & rhC_ui = side_matrices[2];

          // coupling matrices between one side and itself via the element Kss
          std::map<int,std::vector<Epetra_SerialDenseMatrix> >::iterator c2 = Cuiui_coupling.find( sid );
          std::vector<Epetra_SerialDenseMatrix> & Cuiui_matrices = c2->second;
          Epetra_SerialDenseMatrix & eleCuiui = Cuiui_matrices[0];

          emb = XFLUID::EmbCoupling<distype>::TwoSidedImpl(emb_ele,C_uiu,C_uui,rhC_ui,eleCuiui,emb_xyze);
          si  = XFLUID::SideInterface<distype>::Impl(side,side_xyze);
      }
      else
      {
          dserror("InterfaceNitscheTwoSided should not be called for non-fluidfluidcoupling!");
      }

      emb_impl[sid] = emb;
      side_impl[sid] = si;

      // get velocity at integration point of boundary dis

      emb->emb_vel(alediscret,"velaf",alela[0].lm_);
      emb->addembdisp(alediscret,"dispnp",alela[0].lm_,emb_xyze);

//      si->eivel(cutdis,"ivelnp",cutla[0].lm_);
      si->addeidisp(cutdis,"idispnp",cutla[0].lm_,side_xyze);


      double meas_surface = 0.0;
      //TODO: do this not twice ?!
      if(coupling_strategy == INPAR::XFEM::Two_Sided_Mortaring)
      {

        // pre-evaluate for element-size
        // loop gausspoints
        for ( std::vector<DRT::UTILS::GaussIntegration>::const_iterator i=cutintpoints.begin();
            i!=cutintpoints.end();
            ++i )
        {
          const DRT::UTILS::GaussIntegration & gi = *i;
          GEO::CUT::BoundaryCell * bc = bcs[i - cutintpoints.begin()]; // get the corresponding boundary cell

          //gi.Print();

          // TODO: do this transformation not twice!!!
          for ( DRT::UTILS::GaussIntegration::iterator iquad=gi.begin(); iquad!=gi.end(); ++iquad )
          {
            double drs = 0.0; // transformation factor between reference cell and linearized boundary cell
            LINALG::Matrix<3,1> x_gp_lin(true); // gp in xyz-system on linearized interface
            if(bc->Shape()==DRT::Element::tri3 || bc->Shape()==DRT::Element::quad4)
            {
#ifdef BOUNDARYCELL_TRANSFORMATION_OLD
              const LINALG::Matrix<2,1> eta( iquad.Point() ); // xi-coordinates with respect to side

              double drs = 0;

              si->Evaluate(eta,x_side,normal,drs);

              const double fac = drs*iquad.Weight();

              // find element local position of gauss point at interface
              GEO::CUT::Position<distype> pos( my::xyze_, x_side );
              pos.Compute();
              const LINALG::Matrix<3,1> & rst = pos.LocalCoordinates();

  #else
              const LINALG::Matrix<2,1> eta( iquad.Point() ); // eta-coordinates with respect to cell

              normal.Clear();

              // get normal vector on linearized boundary cell, x-coordinates of gaussian point and surface transformation factor
              switch ( bc->Shape() )
              {
              case DRT::Element::tri3:
              {
                bc->Transform<DRT::Element::tri3>(eta, x_gp_lin, normal, drs);
                break;
              }
              case DRT::Element::quad4:
              {
                bc->Transform<DRT::Element::quad4>(eta, x_gp_lin, normal, drs);
                break;
              }
              default:
                throw std::runtime_error( "unsupported integration cell type" );
              }
            }
            else if(bc->Shape()==DRT::Element::dis_none)
            {
              drs = 1.0;
              normal = bc->GetNormalVector();
              const double* gpcord = iquad.Point();
              for (int idim=0;idim<3;idim++)
              {
                 x_gp_lin(idim,0) = gpcord[idim];
              }
            }
#endif


            meas_surface += drs*iquad.Weight();
          }
        }
      }





      //-----------------------------------------------------------------------------------

      // loop gausspoints
      for ( std::vector<DRT::UTILS::GaussIntegration>::const_iterator i=cutintpoints.begin();
            i!=cutintpoints.end();
            ++i )
      {
        const DRT::UTILS::GaussIntegration & gi = *i;
        GEO::CUT::BoundaryCell * bc = bcs[i - cutintpoints.begin()]; // get the corresponding boundary cell

        //gi.Print();

        for ( DRT::UTILS::GaussIntegration::iterator iquad=gi.begin(); iquad!=gi.end(); ++iquad )
        {
          double drs = 0.0; // transformation factor between reference cell and linearized boundary cell
          LINALG::Matrix<3,1> x_gp_lin(true); // gp in xyz-system on linearized interface
          if(bc->Shape()==DRT::Element::tri3 || bc->Shape()==DRT::Element::quad4)
          {
#ifdef BOUNDARYCELL_TRANSFORMATION_OLD
            const LINALG::Matrix<2,1> eta( iquad.Point() ); // xi-coordinates with respect to side

            double drs = 0;

            si->Evaluate(eta,x_side,normal,drs);

            const double fac = drs*iquad.Weight();

            // find element local position of gauss point at interface
            GEO::CUT::Position<distype> pos( my::xyze_, x_side );
            pos.Compute();
            const LINALG::Matrix<3,1> & rst = pos.LocalCoordinates();

  #else
            const LINALG::Matrix<2,1> eta( iquad.Point() ); // eta-coordinates with respect to cell

            normal.Clear();

            // get normal vector on linearized boundary cell, x-coordinates of gaussian point and surface transformation factor
            switch ( bc->Shape() )
            {
            case DRT::Element::tri3:
            {
              bc->Transform<DRT::Element::tri3>(eta, x_gp_lin, normal, drs);
              break;
            }
            case DRT::Element::quad4:
            {
              bc->Transform<DRT::Element::quad4>(eta, x_gp_lin, normal, drs);
              break;
            }
            default:
              throw std::runtime_error( "unsupported integration cell type" );
            }
          }
          else if(bc->Shape()==DRT::Element::dis_none)
          {
            drs = 1.0;
            normal = bc->GetNormalVector();
            const double* gpcord = iquad.Point();
            for (int idim=0;idim<3;idim++)
            {
               x_gp_lin(idim,0) = gpcord[idim];
            }
          }


          const double fac = drs*iquad.Weight();

          // find element local position of gauss point
          GEO::CUT::Position<distype> pos( my::xyze_, x_gp_lin );
          pos.Compute();
          const LINALG::Matrix<3,1> & rst = pos.LocalCoordinates();


          // project gaussian point from linearized interface to warped side (get local side coordinates)
          LINALG::Matrix<2,1> xi_side(true);
          si->ProjectOnSide(x_gp_lin, x_side, xi_side);

#endif


          // evaluate embedded element shape functions
          emb->EvaluateEmb( x_side );


          // evaluate shape functions
          DRT::UTILS::shape_function<distype>( rst, my::funct_ );

          // evaluate shape functions and derivatives at integration point
          //          EvalShapeFuncAndDerivsAtIntPoint(iquad,eid);

          DRT::UTILS::shape_function_deriv1<distype>(rst,my::deriv_);
          my::xjm_.MultiplyNT(my::deriv_,my::xyze_);
          my::det_ = my::xji_.Invert(my::xjm_);


          // compute global first derivates
          my::derxy_.Multiply(my::xji_,my::deriv_);


          // get velocity at integration point
          // (values at n+alpha_F for generalized-alpha scheme, n+1 otherwise)
          my::velint_.Multiply(evelaf,my::funct_);


          // get velocity derivatives at integration point
          // (values at n+alpha_F for generalized-alpha scheme, n+1 otherwise)
          my::vderxy_.MultiplyNT(evelaf,my::derxy_);

          // get pressure at integration point
          // (value at n+alpha_F for generalized-alpha scheme, n+1 otherwise)
          double press = my::funct_.Dot(epreaf);

          const double timefacfac = my::f3Parameter_->timefac_ * fac;


          //-----------------------------------------------------------------------------------

          double stabfac = 0.0;         // Nitsche stabilization factor
          double stabfac_conv = 0.0;    // Nitsche convective stabilization factor


          // define stabilization parameters and mortaring weights
          double visceff_max = my::visceff_;

          double kappa1 = 1.0;      // Xfluid-sided mortaring

          if(coupling_strategy == INPAR::XFEM::Xfluid_Sided_Mortaring)
          {
            if(meas_partial_volume < 1e-008) dserror(" measure of cut partial volume is smaller than 1d-008: %d Attention with increasing Nitsche-Parameter!!!", meas_partial_volume);

            stabfac = nitsche_stab * visceff_max * meas_surface / meas_partial_volume;
            stabfac_conv = nitsche_stab_conv * meas_surface / meas_partial_volume;

            kappa1 = 1.0;
          }
          else if(coupling_strategy == INPAR::XFEM::Embedded_Sided_Mortaring)
          {
            // get element diameter
            double hk_emb = 0.0;
            emb->element_length(hk_emb);

            if(hk_emb < 1e-006) dserror("element length is smaller than 1e-006");
            stabfac = nitsche_stab * visceff_max /hk_emb;
            stabfac_conv = nitsche_stab_conv * max(1.0, max(fabs(my::velint_.Dot(normal)) , visceff_max /hk_emb) );

            kappa1 = 0.0;
          }
          else if(coupling_strategy == INPAR::XFEM::Two_Sided_Mortaring)
          {
            if(meas_partial_volume < 1e-008) dserror(" measure of cut partial volume is smaller than 1d-008: %d Attention with increasing Nitsche-Parameter!!!", meas_partial_volume);
            stabfac = nitsche_stab * visceff_max * meas_surface / meas_partial_volume;
            stabfac_conv = nitsche_stab_conv * meas_surface / meas_partial_volume;

            kappa1 = 0.5;
          }
          else dserror("coupling strategy not known");

          if( kappa1 > 1.0 || kappa1 < 0.0) dserror("Nitsche weights for inverse estimate kappa1 lies not in [0,1]: %d", kappa1);

          double kappa2 = 1.0-kappa1;







          if(fluidfluidcoupling)
          {
            bool bg_mortaring = false; // one-sided background fluid mortaring (kappa1=1, kappa2=0)

            if(coupling_strategy == INPAR::XFEM::Xfluid_Sided_Mortaring) bg_mortaring = true;

            //zero velocity jump for fluidfluidcoupling
            LINALG::Matrix<my::nsd_,1> ivelint_WDBC_JUMP(true);


            emb->buildCouplingMatricesNitscheTwoSided( elemat1_epetra,          // standard bg-bg-matrix
                elevec1_epetra,          // standard bg-rhs
                fluidfluidcoupling,      // assemble coupling terms (yes/no)
                bg_mortaring,            // yes: background-sided mortaring, no: coupling between two meshes (mixed mortaring)
                normal,                  // normal vector
                timefacfac,              // theta*dt
                my::visceff_,                // viscosity in background fluid
                my::visceff_,                // viscosity in embedded fluid
                kappa1,                  // mortaring weighting
                kappa2,                  // mortaring weighting
                stabfac,                 // Nitsche non-dimensionless stabilization factor
                stabfac_conv,            // Nitsche convective non-dimensionless stabilization factor
                my::funct_,                  // bg shape functions
                my::derxy_,                  // bg deriv
                my::vderxy_,                 // bg deriv^n
                press,                   // bg p^n
                my::velint_,                 // bg u^n
                ivelint_WDBC_JUMP         // Dirichlet velocity vector or prescribed jump vector
            );


          }
          if(!fluidfluidcoupling)
          {
            dserror(" no two sided mortaring for non-fluidfluidcoupling");
          }



        }
      }
    }




  if ( fluidfluidcoupling )
  {


    // build Gsui and Guis
    int ipatchsizesbefore = 0;
    for (std::map<int, std::vector<Epetra_SerialDenseMatrix> >::const_iterator m=Cuiui_coupling.begin();
        m!=Cuiui_coupling.end(); ++m)
    {

      int bid = m->first;
      std::vector<Epetra_SerialDenseMatrix> & Cuiui_mats = Cuiui_coupling[bid];

      //Cuiui matrices in Cuiui_mats[0]

      // assemble Cuiui
      for ( int ic=0; ic<Cuiui_mats[0].N(); ++ic) // Cuiui includes only ui,ui coupling, not (ui,p) ...
      {
        for ( int ir=0; ir<Cuiui_mats[0].M(); ++ir )
        {
          Cuiui(ir+ipatchsizesbefore,ic+ipatchsizesbefore) = Cuiui_mats[0](ir,ic);
        }
      }

      ipatchsizesbefore += Cuiui_mats[0].N();

    }
  }
}

/*--------------------------------------------------------------------------------
 *--------------------------------------------------------------------------------*/
template <>
void DRT::ELEMENTS::FluidEleCalcXFEM<DRT::Element::tri3>::ElementXfemInterfaceNitscheTwoSided(
  DRT::ELEMENTS::Fluid3 * ele,
  DRT::Discretization & dis,
  const std::vector<int> & lm,
  const DRT::UTILS::GaussIntegration & intpoints,
  DRT::Discretization & cutdis,
  const std::map<int, std::vector<GEO::CUT::BoundaryCell*> > & bcells,
  const std::map<int, std::vector<DRT::UTILS::GaussIntegration> > & bintpoints,
  std::map<int, std::vector<Epetra_SerialDenseMatrix> > & side_coupling,
  Teuchos::ParameterList&    params,
  DRT::Discretization &  alediscret,
  map<int,int> & boundary_emb_gid_map,
  Epetra_SerialDenseMatrix&  elemat1_epetra,
  Epetra_SerialDenseVector&  elevec1_epetra,
  Epetra_SerialDenseMatrix&  Cuiui
  )
{
  dserror( "distype not supported" );
}

/*--------------------------------------------------------------------------------
 *--------------------------------------------------------------------------------*/
template <>
void DRT::ELEMENTS::FluidEleCalcXFEM<DRT::Element::tri6>::ElementXfemInterfaceNitscheTwoSided(
  DRT::ELEMENTS::Fluid3 * ele,
  DRT::Discretization & dis,
  const std::vector<int> & lm,
  const DRT::UTILS::GaussIntegration & intpoints,
  DRT::Discretization & cutdis,
  const std::map<int, std::vector<GEO::CUT::BoundaryCell*> > & bcells,
  const std::map<int, std::vector<DRT::UTILS::GaussIntegration> > & bintpoints,
  std::map<int, std::vector<Epetra_SerialDenseMatrix> > & side_coupling,
  Teuchos::ParameterList&    params,
  DRT::Discretization &  alediscret,
  map<int,int> & boundary_emb_gid_map,
  Epetra_SerialDenseMatrix&  elemat1_epetra,
  Epetra_SerialDenseVector&  elevec1_epetra,
  Epetra_SerialDenseMatrix&  Cuiui
  )
{
  dserror( "distype not supported" );
}

/*--------------------------------------------------------------------------------
 *--------------------------------------------------------------------------------*/
template <>
void DRT::ELEMENTS::FluidEleCalcXFEM<DRT::Element::quad4>::ElementXfemInterfaceNitscheTwoSided(
  DRT::ELEMENTS::Fluid3 * ele,
  DRT::Discretization & dis,
  const std::vector<int> & lm,
  const DRT::UTILS::GaussIntegration & intpoints,
  DRT::Discretization & cutdis,
  const std::map<int, std::vector<GEO::CUT::BoundaryCell*> > & bcells,
  const std::map<int, std::vector<DRT::UTILS::GaussIntegration> > & bintpoints,
  std::map<int, std::vector<Epetra_SerialDenseMatrix> > & side_coupling,
  Teuchos::ParameterList&    params,
  DRT::Discretization &  alediscret,
  map<int,int> & boundary_emb_gid_map,
  Epetra_SerialDenseMatrix&  elemat1_epetra,
  Epetra_SerialDenseVector&  elevec1_epetra,
  Epetra_SerialDenseMatrix&  Cuiui
  )
{
  dserror( "distype not supported" );
}

/*--------------------------------------------------------------------------------
 *--------------------------------------------------------------------------------*/
template <>
void DRT::ELEMENTS::FluidEleCalcXFEM<DRT::Element::quad8>::ElementXfemInterfaceNitscheTwoSided(
  DRT::ELEMENTS::Fluid3 * ele,
  DRT::Discretization & dis,
  const std::vector<int> & lm,
  const DRT::UTILS::GaussIntegration & intpoints,
  DRT::Discretization & cutdis,
  const std::map<int, std::vector<GEO::CUT::BoundaryCell*> > & bcells,
  const std::map<int, std::vector<DRT::UTILS::GaussIntegration> > & bintpoints,
  std::map<int, std::vector<Epetra_SerialDenseMatrix> > & side_coupling,
  Teuchos::ParameterList&    params,
  DRT::Discretization &  alediscret,
  map<int,int> & boundary_emb_gid_map,
  Epetra_SerialDenseMatrix&  elemat1_epetra,
  Epetra_SerialDenseVector&  elevec1_epetra,
  Epetra_SerialDenseMatrix&  Cuiui
  )
{
  dserror( "distype not supported" );
}

/*--------------------------------------------------------------------------------
 *--------------------------------------------------------------------------------*/
template <>
void DRT::ELEMENTS::FluidEleCalcXFEM<DRT::Element::quad9>::ElementXfemInterfaceNitscheTwoSided(
  DRT::ELEMENTS::Fluid3 * ele,
  DRT::Discretization & dis,
  const std::vector<int> & lm,
  const DRT::UTILS::GaussIntegration & intpoints,
  DRT::Discretization & cutdis,
  const std::map<int, std::vector<GEO::CUT::BoundaryCell*> > & bcells,
  const std::map<int, std::vector<DRT::UTILS::GaussIntegration> > & bintpoints,
  std::map<int, std::vector<Epetra_SerialDenseMatrix> > & side_coupling,
  Teuchos::ParameterList&    params,
  DRT::Discretization &  alediscret,
  map<int,int> & boundary_emb_gid_map,
  Epetra_SerialDenseMatrix&  elemat1_epetra,
  Epetra_SerialDenseVector&  elevec1_epetra,
  Epetra_SerialDenseMatrix&  Cuiui
  )
{
  dserror( "distype not supported" );
}

/*--------------------------------------------------------------------------------
 *--------------------------------------------------------------------------------*/
template <>
void DRT::ELEMENTS::FluidEleCalcXFEM<DRT::Element::nurbs9>::ElementXfemInterfaceNitscheTwoSided(
  DRT::ELEMENTS::Fluid3 * ele,
  DRT::Discretization & dis,
  const std::vector<int> & lm,
  const DRT::UTILS::GaussIntegration & intpoints,
  DRT::Discretization & cutdis,
  const std::map<int, std::vector<GEO::CUT::BoundaryCell*> > & bcells,
  const std::map<int, std::vector<DRT::UTILS::GaussIntegration> > & bintpoints,
  std::map<int, std::vector<Epetra_SerialDenseMatrix> > & side_coupling,
  Teuchos::ParameterList&    params,
  DRT::Discretization &  alediscret,
  map<int,int> & boundary_emb_gid_map,
  Epetra_SerialDenseMatrix&  elemat1_epetra,
  Epetra_SerialDenseVector&  elevec1_epetra,
  Epetra_SerialDenseMatrix&  Cuiui
  )
{
  dserror( "distype not supported" );
}


//
///*--------------------------------------------------------------------------------
// *--------------------------------------------------------------------------------*/
//template <DRT::Element::DiscretizationType distype>
//void FluidEleCalcXFEM<distype>::ElementXfemInterfaceNeumann(
//  DRT::ELEMENTS::Fluid3 * ele,
//  DRT::Discretization & dis,
//  const std::vector<int> & lm,
//  const DRT::UTILS::GaussIntegration & intpoints,
//  DRT::Discretization & cutdis,
//  const std::map<int, std::vector<GEO::CUT::BoundaryCell*> > & bcells,
//  const std::map<int, std::vector<DRT::UTILS::GaussIntegration> > & bintpoints,
//  std::map<int, std::vector<Epetra_SerialDenseMatrix> > & side_coupling,
//  Teuchos::ParameterList&    params,
//  Epetra_SerialDenseMatrix&  elemat1_epetra,
//  Epetra_SerialDenseVector&  elevec1_epetra,
//  Epetra_SerialDenseMatrix&  Cuiui
//  )
//{
//	  //----------------------------------------------------------------------------
//	  //                         ELEMENT GEOMETRY
//	  //----------------------------------------------------------------------------
//
//	  // get node coordinates
//	  GEO::fillInitialPositionArray< distype, my::nsd_, LINALG::Matrix<my::nsd_,my::nen_> >( ele, my::xyze_ );
//
//	  LINALG::Matrix<my::nsd_,my::nen_> evelaf(true);
//	  LINALG::Matrix<my::nen_,1> epreaf(true);
//	  ExtractValuesFromGlobalVector(dis, lm, *my::rotsymmpbc_, &evelaf, &epreaf, "velaf");
//
//
//	  // integrate surface
//
//	  DRT::Element::LocationArray cutla( 1 );
//
//	  LINALG::Matrix<3,1> normal;
//	  LINALG::Matrix<3,1> x_side;
//
//	  bool fluidfluidcoupling = false;
//
//	  std::map<int, Teuchos::RCP<XFLUID::SideInterface<distype> > > side_impl;
//	  Teuchos::RCP<XFLUID::SideInterface<distype> > si;
//
//	  // find all the intersecting elements of actele
//	  std::set<int> begids;
//	  for (std::map<int,  std::vector<GEO::CUT::BoundaryCell*> >::const_iterator bc=bcells.begin();
//	       bc!=bcells.end(); ++bc )
//	  {
//	    int sid = bc->first;
//	    begids.insert(sid);
//	  }
//
//	  std::map<int, std::vector<Epetra_SerialDenseMatrix> > Cuiui_coupling;
//
//	  // lm vector of all intersecting elements
//	  std::vector<int> patchelementslmv;
//	  std::vector<int> patchelementslmowner;
//	  for (std::set<int>::const_iterator bgid=begids.begin(); bgid!=begids.end(); ++bgid)
//	  {
//	    DRT::Element * side = cutdis.gElement(*bgid);
//	    vector<int> patchlm;
//	    vector<int> patchlmowner;
//	    vector<int> patchlmstride;
//	    side->LocationVector(cutdis, patchlm, patchlmowner, patchlmstride);
//
//	    patchelementslmv.reserve( patchelementslmv.size() + patchlm.size());
//	    patchelementslmv.insert(patchelementslmv.end(), patchlm.begin(), patchlm.end());
//
//	    patchelementslmowner.reserve( patchelementslmowner.size() + patchlmowner.size());
//	    patchelementslmowner.insert( patchelementslmowner.end(), patchlmowner.begin(), patchlmowner.end());
//
//	    std::vector<Epetra_SerialDenseMatrix> & Cuiui_matrices = Cuiui_coupling[*bgid];
//
//	    Cuiui_matrices.resize(2);
//	    Cuiui_matrices[0].Reshape(my::nen_*6,patchlm.size()); //Gsui
//	    Cuiui_matrices[1].Reshape(patchlm.size(),my::nen_*6); //Guis
//	  }
//
//	  Epetra_SerialDenseMatrix Gsui(my::nen_*6,patchelementslmv.size());
//	  Epetra_SerialDenseMatrix Guis(patchelementslmv.size(),my::nen_*6);
//	  Epetra_SerialDenseMatrix InvKss(my::nen_*6,my::nen_*6);
//	  Epetra_SerialDenseMatrix GuisInvKss(patchelementslmv.size(),my::nen_*6);
//
//	  // map of side-element id and Guass points
//	  for ( std::map<int, std::vector<DRT::UTILS::GaussIntegration> >::const_iterator i=bintpoints.begin();
//	        i!=bintpoints.end();
//	        ++i )
//	  {
//	    int sid = i->first;
//	    const std::vector<DRT::UTILS::GaussIntegration> & cutintpoints = i->second;
//
//	    std::map<int, std::vector<GEO::CUT::BoundaryCell*> >::const_iterator j = bcells.find( sid );
//	    if ( j==bcells.end() )
//	      dserror( "missing boundary cell" );
//
//	    const std::vector<GEO::CUT::BoundaryCell*> & bcs = j->second;
//	    if ( bcs.size()!=cutintpoints.size() )
//	      dserror( "boundary cell integration rules mismatch" );
//
//	    DRT::Element * side = cutdis.gElement( sid );
//	    side->LocationVector(cutdis,cutla,false);
//
//	    const int numnodes = side->NumNode();
//	    DRT::Node ** nodes = side->Nodes();
//	    Epetra_SerialDenseMatrix side_xyze( 3, numnodes );
//	    for ( int i=0; i<numnodes; ++i )
//	    {
//	      const double * x = nodes[i]->X();
//	      std::copy( x, x+3, &side_xyze( 0, i ) );
//	    }
//
//	    std::map<int,std::vector<Epetra_SerialDenseMatrix> >::iterator c = side_coupling.find( sid );
//
//	    std::vector<Epetra_SerialDenseMatrix> & side_matrices = c->second;
//
//	    if ( side_matrices.size()==3 )
//	      fluidfluidcoupling = true;
//
//
//	    if(fluidfluidcoupling)
//	    {
//	        Epetra_SerialDenseMatrix & C_uiu  = side_matrices[0];
//	        Epetra_SerialDenseMatrix & C_uui  = side_matrices[1];
//	        Epetra_SerialDenseMatrix & rhC_ui = side_matrices[2];
//
//	        std::map<int,std::vector<Epetra_SerialDenseMatrix> >::iterator c2 = Cuiui_coupling.find( sid );
//	        std::vector<Epetra_SerialDenseMatrix> & Cuiui_matrices = c2->second;
//	        Epetra_SerialDenseMatrix & eleGsui = Cuiui_matrices[0];
//	        Epetra_SerialDenseMatrix & eleGuis = Cuiui_matrices[1];
//	        Epetra_SerialDenseMatrix  eleGuisKssInv;
//
//	        si = XFLUID::SideInterface<distype>::Impl(side,C_uiu,C_uui,rhC_ui,eleGsui,eleGuis,side_xyze);
//	    }
//	    else
//	    {
//	        si = XFLUID::SideInterface<distype>::Impl(side,side_xyze);
//	    }
//
//	    side_impl[sid] = si;
//
//	    // get velocity at integration point of boundary dis
//	    si->eivel(cutdis,"ivelnp",cutla[0].lm_);
//	    si->addeidisp(cutdis,"idispnp",cutla[0].lm_,side_xyze);
//
//	    // loop gausspoints
//	    for ( std::vector<DRT::UTILS::GaussIntegration>::const_iterator i=cutintpoints.begin();
//	          i!=cutintpoints.end();
//	          ++i )
//	    {
//	      const DRT::UTILS::GaussIntegration & gi = *i;
//	      //const GEO::CUT::BoundaryCell & bc = *bcs[i - cutintpoints.begin()];
//	      GEO::CUT::BoundaryCell * bc = bcs[i - cutintpoints.begin()]; // get the corresponding boundary cell
//
//	      //gi.Print();
//
//	      for ( DRT::UTILS::GaussIntegration::iterator iquad=gi.begin(); iquad!=gi.end(); ++iquad )
//	      {
//#ifdef BOUNDARYCELL_TRANSFORMATION_OLD
//        const LINALG::Matrix<2,1> eta( iquad.Point() ); // xi-coordinates with respect to side
//
//        double drs = 0;
//
//        si->Evaluate(eta,x_side,normal,drs);
//
//        const double fac = drs*iquad.Weight();
//
//        // find element local position of gauss point at interface
//        GEO::CUT::Position<distype> pos( my::xyze_, x_side );
//        pos.Compute();
//        const LINALG::Matrix<3,1> & rst = pos.LocalCoordinates();
//
//#else
//        const LINALG::Matrix<2,1> eta( iquad.Point() ); // eta-coordinates with respect to cell
//
//        double drs = 0; // transformation factor between reference cell and linearized boundary cell
//
//        LINALG::Matrix<3,1> x_gp_lin(true); // gp in xyz-system on linearized interface
//        normal.Clear();
//
//        // get normal vector on linearized boundary cell, x-coordinates of gaussian point and surface transformation factor
//        switch ( bc->Shape() )
//        {
//        case DRT::Element::tri3:
//        {
//            bc->Transform<DRT::Element::tri3>(eta, x_gp_lin, normal, drs);
//          break;
//        }
//        case DRT::Element::quad4:
//        {
//            bc->Transform<DRT::Element::quad4>(eta, x_gp_lin, normal, drs);
//          break;
//        }
//        default:
//          throw std::runtime_error( "unsupported integration cell type" );
//        }
//
//
//        const double fac = drs*iquad.Weight();
//
//        // find element local position of gauss point
//        GEO::CUT::Position<distype> pos( my::xyze_, x_gp_lin );
//        pos.Compute();
//        const LINALG::Matrix<3,1> & rst = pos.LocalCoordinates();
//
//
//        // project gaussian point from linearized interface to warped side (get local side coordinates)
//        LINALG::Matrix<2,1> xi_side(true);
//        si->ProjectOnSide(x_gp_lin, x_side, xi_side);
//
//#endif
//
//	        // evaluate shape functions
//	        DRT::UTILS::shape_function<distype>( rst, funct_ );
//
//
////	        // get velocity at integration point
////	        // (values at n+alpha_F for generalized-alpha scheme, n+1 otherwise)
////	        velint_.Multiply(evelaf,funct_);
////
////
//////	        get pressure at integration point
//////	        (value at n+alpha_F for generalized-alpha scheme, n+1 otherwise)
////	        double press = funct_.Dot(epreaf);
//
////	        const double timefacfac = my::f3Parameter_->timefac_ * fac;
//
//
//	        // pure Neumann boundary condition
//	        // evaluate  Neumann boundary condition
//        LINALG::Matrix<my::nsd_,1> h_N(true);
//
//        double shear_xx =   0.0;
//        double shear_xy =  30.0;
//        double shear_xz =  60.0;
//
//        double shear_yx =   0.0;
//        double shear_yy =   0.0;
//        double shear_yz =   0.0;
//
//        double shear_zx =   0.0;
//        double shear_zy =   0.0;
//        double shear_zz =   0.0;
//
//        double pres_N= 5.0;
//
//        // attention normal has the wrong sign!!! (is the normal vector normal to the "side")
//        h_N(0) = -pres_N*normal(0) + 2.0*my::visceff_*(      shear_xx          *normal(0)
//                                                   +0.5*(shear_xy+shear_yx)*normal(1)
//                                                   +0.5*(shear_xz+shear_zx)*normal(2)  );
//        h_N(1) = -pres_N*normal(1) + 2.0*my::visceff_*( 0.5*(shear_yx+shear_xy)*normal(0)
//                                                   +     shear_yy          *normal(1)
//                                                   +0.5*(shear_yz+shear_zy)*normal(2)  );
//        h_N(2) = -pres_N*normal(2) + 2.0*my::visceff_*( 0.5*(shear_zx+shear_xz)*normal(0)
//                                                   +0.5*(shear_zy+shear_yz)*normal(1)
//                                                   +     shear_zz          *normal(2)  );
//
//
//        for(int r=0; r<my::nen_; r++)
//        {
//        	int rind = r*(my::nsd_+1);
//        	int dVelx=rind+0;
//        	int dVely=rind+1;
//        	int dVelz=rind+2;
//
//        	elevec1_epetra(dVelx,0) += funct_(r)*fac*h_N(0);
//        	elevec1_epetra(dVely,0) += funct_(r)*fac*h_N(1);
//        	elevec1_epetra(dVelz,0) += funct_(r)*fac*h_N(2);
//
//        }
//
//
//
//	      }
//	    }
//	  }
//
//
//}
//
///*--------------------------------------------------------------------------------
// *--------------------------------------------------------------------------------*/
//template <>
//void FluidEleCalcXFEM<DRT::Element::tri3>::ElementXfemInterfaceNeumann(
//  DRT::ELEMENTS::Fluid3 * ele,
//  DRT::Discretization & dis,
//  const std::vector<int> & lm,
//  const DRT::UTILS::GaussIntegration & intpoints,
//  DRT::Discretization & cutdis,
//  const std::map<int, std::vector<GEO::CUT::BoundaryCell*> > & bcells,
//  const std::map<int, std::vector<DRT::UTILS::GaussIntegration> > & bintpoints,
//  std::map<int, std::vector<Epetra_SerialDenseMatrix> > & side_coupling,
//  //std::map<int, std::vector<RCP<Epetra_SerialDenseMatrix> > > side_coupling,
//  Teuchos::ParameterList&    params,
//  Epetra_SerialDenseMatrix&  elemat1_epetra,
//  Epetra_SerialDenseVector&  elevec1_epetra,
//  Epetra_SerialDenseMatrix&  Cuiui
//  )
//{
//  dserror( "distype not supported" );
//}
//
///*--------------------------------------------------------------------------------
// *--------------------------------------------------------------------------------*/
//template <>
//void FluidEleCalcXFEM<DRT::Element::tri6>::ElementXfemInterfaceNeumann(
//  DRT::ELEMENTS::Fluid3 * ele,
//  DRT::Discretization & dis,
//  const std::vector<int> & lm,
//  const DRT::UTILS::GaussIntegration & intpoints,
//  DRT::Discretization & cutdis,
//  const std::map<int, std::vector<GEO::CUT::BoundaryCell*> > & bcells,
//  const std::map<int, std::vector<DRT::UTILS::GaussIntegration> > & bintpoints,
//  std::map<int, std::vector<Epetra_SerialDenseMatrix> > & side_coupling,
////  std::map<int, std::vector<RCP<Epetra_SerialDenseMatrix> > > side_coupling,
//  Teuchos::ParameterList&    params,
//  Epetra_SerialDenseMatrix&  elemat1_epetra,
//  Epetra_SerialDenseVector&  elevec1_epetra,
//  Epetra_SerialDenseMatrix&  Cuiui
//  )
//{
//  dserror( "distype not supported" );
//}
//
///*--------------------------------------------------------------------------------
// *--------------------------------------------------------------------------------*/
//template <>
//void FluidEleCalcXFEM<DRT::Element::quad4>::ElementXfemInterfaceNeumann(
//  DRT::ELEMENTS::Fluid3 * ele,
//  DRT::Discretization & dis,
//  const std::vector<int> & lm,
//  const DRT::UTILS::GaussIntegration & intpoints,
//  DRT::Discretization & cutdis,
//  const std::map<int, std::vector<GEO::CUT::BoundaryCell*> > & bcells,
//  const std::map<int, std::vector<DRT::UTILS::GaussIntegration> > & bintpoints,
//  std::map<int, std::vector<Epetra_SerialDenseMatrix> > & side_coupling,
////  std::map<int, std::vector<RCP<Epetra_SerialDenseMatrix> > > side_coupling,
//  Teuchos::ParameterList&    params,
//  Epetra_SerialDenseMatrix&  elemat1_epetra,
//  Epetra_SerialDenseVector&  elevec1_epetra,
//  Epetra_SerialDenseMatrix&  Cuiui
//  )
//{
//  dserror( "distype not supported" );
//}
//
///*--------------------------------------------------------------------------------
// *--------------------------------------------------------------------------------*/
//template <>
//void FluidEleCalcXFEM<DRT::Element::quad8>::ElementXfemInterfaceNeumann(
//  DRT::ELEMENTS::Fluid3 * ele,
//  DRT::Discretization & dis,
//  const std::vector<int> & lm,
//  const DRT::UTILS::GaussIntegration & intpoints,
//  DRT::Discretization & cutdis,
//  const std::map<int, std::vector<GEO::CUT::BoundaryCell*> > & bcells,
//  const std::map<int, std::vector<DRT::UTILS::GaussIntegration> > & bintpoints,
//  std::map<int, std::vector<Epetra_SerialDenseMatrix> > & side_coupling,
//  //std::map<int, std::vector<RCP<Epetra_SerialDenseMatrix> > > side_coupling,
//  Teuchos::ParameterList&    params,
//  Epetra_SerialDenseMatrix&  elemat1_epetra,
//  Epetra_SerialDenseVector&  elevec1_epetra,
//  Epetra_SerialDenseMatrix&  Cuiui
//  )
//{
//  dserror( "distype not supported" );
//}
//
///*--------------------------------------------------------------------------------
// *--------------------------------------------------------------------------------*/
//template <>
//void FluidEleCalcXFEM<DRT::Element::quad9>::ElementXfemInterfaceNeumann(
//  DRT::ELEMENTS::Fluid3 * ele,
//  DRT::Discretization & dis,
//  const std::vector<int> & lm,
//  const DRT::UTILS::GaussIntegration & intpoints,
//  DRT::Discretization & cutdis,
//  const std::map<int, std::vector<GEO::CUT::BoundaryCell*> > & bcells,
//  const std::map<int, std::vector<DRT::UTILS::GaussIntegration> > & bintpoints,
//  std::map<int, std::vector<Epetra_SerialDenseMatrix> > & side_coupling,
//  //std::map<int, std::vector<RCP<Epetra_SerialDenseMatrix> > > side_coupling,
//  Teuchos::ParameterList&    params,
//  Epetra_SerialDenseMatrix&  elemat1_epetra,
//  Epetra_SerialDenseVector&  elevec1_epetra,
//  Epetra_SerialDenseMatrix&  Cuiui
//  )
//{
//  dserror( "distype not supported" );
//}
//
///*--------------------------------------------------------------------------------
// *--------------------------------------------------------------------------------*/
//template <>
//void FluidEleCalcXFEM<DRT::Element::nurbs9>::ElementXfemInterfaceNeumann(
//  DRT::ELEMENTS::Fluid3 * ele,
//  DRT::Discretization & dis,
//  const std::vector<int> & lm,
//  const DRT::UTILS::GaussIntegration & intpoints,
//  DRT::Discretization & cutdis,
//  const std::map<int, std::vector<GEO::CUT::BoundaryCell*> > & bcells,
//  const std::map<int, std::vector<DRT::UTILS::GaussIntegration> > & bintpoints,
//  std::map<int, std::vector<Epetra_SerialDenseMatrix> > & side_coupling,
////  std::map<int, std::vector<RCP<Epetra_SerialDenseMatrix> > > side_coupling,
//  Teuchos::ParameterList&    params,
//  Epetra_SerialDenseMatrix&  elemat1_epetra,
//  Epetra_SerialDenseVector&  elevec1_epetra,
//  Epetra_SerialDenseMatrix&  Cuiui
//  )
//{
//  dserror( "distype not supported" );
//}

  }
}


/*--------------------------------------------------------------------------------
 *--------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidEleCalcXFEM<distype>::CalculateContinuityXFEM(
  DRT::ELEMENTS::Fluid3 * ele,
  DRT::Discretization & dis,
  const std::vector<int> & lm,
  Epetra_SerialDenseVector&  elevec1_epetra,
  const DRT::UTILS::GaussIntegration & intpoints
  )
{
  LINALG::Matrix<(my::nsd_+1)*my::nen_,1> elevec1(elevec1_epetra,true);
  int eid = ele->Id();

  //------------------------------------------------------------------------
  //  start loop over integration points
  //------------------------------------------------------------------------
  for ( DRT::UTILS::GaussIntegration::const_iterator iquad=intpoints.begin(); iquad!=intpoints.end(); ++iquad )
  {
    // evaluate shape functions and derivatives at integration point
    my::EvalShapeFuncAndDerivsAtIntPoint(iquad,eid);

    for(int ui=0; ui<my::nen_; ++ui)
    {
      for (int idim = 0; idim <my::nsd_; ++idim)
      {
        const int fui = (my::nsd_+1)*ui;
        /* continuity term */
        /*
             /           \
            |             |
            | nabla o Du  |
            |             |
             \           /
        */

        elevec1(fui+idim) += my::fac_*my::derxy_(idim,ui);
      }
    }

  }
  return;
}

/*--------------------------------------------------------------------------------
 *--------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidEleCalcXFEM<distype>::CalculateContinuityXFEM(
  DRT::ELEMENTS::Fluid3 * ele,
  DRT::Discretization & dis,
  const std::vector<int> & lm,
  Epetra_SerialDenseVector&  elevec1_epetra
  )
{

  CalculateContinuityXFEM(ele,
                          dis,
                          lm,
                          elevec1_epetra,
                          my::intpoints_);
}


// Ursula is responsible for this comment!
template class DRT::ELEMENTS::FluidEleCalcXFEM<DRT::Element::hex8>;
template class DRT::ELEMENTS::FluidEleCalcXFEM<DRT::Element::hex20>;
template class DRT::ELEMENTS::FluidEleCalcXFEM<DRT::Element::hex27>;
template class DRT::ELEMENTS::FluidEleCalcXFEM<DRT::Element::tet4>;
template class DRT::ELEMENTS::FluidEleCalcXFEM<DRT::Element::tet10>;
template class DRT::ELEMENTS::FluidEleCalcXFEM<DRT::Element::wedge6>;
template class DRT::ELEMENTS::FluidEleCalcXFEM<DRT::Element::pyramid5>;
template class DRT::ELEMENTS::FluidEleCalcXFEM<DRT::Element::quad4>;
template class DRT::ELEMENTS::FluidEleCalcXFEM<DRT::Element::quad8>;
template class DRT::ELEMENTS::FluidEleCalcXFEM<DRT::Element::quad9>;
template class DRT::ELEMENTS::FluidEleCalcXFEM<DRT::Element::tri3>;
template class DRT::ELEMENTS::FluidEleCalcXFEM<DRT::Element::tri6>;
template class DRT::ELEMENTS::FluidEleCalcXFEM<DRT::Element::nurbs9>;
template class DRT::ELEMENTS::FluidEleCalcXFEM<DRT::Element::nurbs27>;


