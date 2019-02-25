/*----------------------------------------------------------------------*/
/*!
\file fluid_ele_calc_xfem_coupling_impl.cpp

\brief Classes for interface coupling in the XFEM

\level 2

<pre>
\maintainer  Ager Christoph
             ager@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289-15249
</pre>
*/
/*----------------------------------------------------------------------*/

#include "fluid_ele_calc_xfem_coupling.H"
#include "fluid_ele_calc_xfem_coupling_impl.H"
#include "fluid_ele_parameter_xfem.H"

#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_element_integration_select.H"

#include "../linalg/linalg_utils.H"

#include "../drt_cut/cut_boundarycell.H"
#include "../drt_cut/cut_position.H"

#include "../drt_xfem/xfem_interface_utils.H"

#include <Teuchos_TimeMonitor.hpp>

namespace DRT
{
  namespace ELEMENTS
  {
    namespace XFLUID
    {
      /*----------------------------------------------------------------------*
       *----------------------------------------------------------------------*/
      template <DRT::Element::DiscretizationType distype,
          DRT::Element::DiscretizationType slave_distype, unsigned int slave_numdof>
      void SlaveElementRepresentation<distype, slave_distype, slave_numdof>::AddSlaveEleDisp(
          const DRT::Discretization& slavedis,  ///< coupling slave discretization
          const std::vector<int>& lm            ///< local map
      )
      {
        std::vector<double> mymatrix(lm.size());
        AddSlaveEleDisp(slavedis, lm, mymatrix);
        return;
      }

      /*----------------------------------------------------------------------*
       *----------------------------------------------------------------------*/
      template <DRT::Element::DiscretizationType distype,
          DRT::Element::DiscretizationType slave_distype, unsigned int slave_numdof>
      void SlaveElementRepresentation<distype, slave_distype, slave_numdof>::AddSlaveEleDisp(
          const DRT::Discretization& slavedis,  ///< coupling slave discretization
          const std::vector<int>& lm,           ///< local map
          std::vector<double>& mymatrix         ///< slave element displacement vector
      )
      {
        // leave, if displacements are not set
        if (!slavedis.HasState(disp_statename_)) return;
        // get state of the global vector
        Teuchos::RCP<const Epetra_Vector> matrix_state = slavedis.GetState(disp_statename_);
        if (matrix_state == Teuchos::null)
          dserror("Cannot get state vector %s", disp_statename_.c_str());

        // extract local values of the global vector
        DRT::UTILS::ExtractMyValues(*matrix_state, mymatrix, lm);

        for (unsigned inode = 0; inode < slave_nen_; ++inode)  // number of nodes
        {
          for (unsigned idim = 0; idim < nsd_; ++idim)  // number of dimensions
          {
            (slave_disp_)(idim, inode) =
                mymatrix[idim + (inode * slave_numdof)];  // attention! for fluidfluid disp state
                                                          // vector has 3+1 dofs for displacement
                                                          // (the same as for (u,p))
          }
        }

        // add the displacement of the interface
        for (unsigned inode = 0; inode < slave_nen_; ++inode)
        {
          slave_xyze_(0, inode) += slave_disp_(0, inode);
          slave_xyze_(1, inode) += slave_disp_(1, inode);
          slave_xyze_(2, inode) += slave_disp_(2, inode);
        }

        return;
      }

      /*----------------------------------------------------------------------*
       *----------------------------------------------------------------------*/
      template <DRT::Element::DiscretizationType distype,
          DRT::Element::DiscretizationType slave_distype, unsigned int slave_numdof>
      void SlaveElementRepresentation<distype, slave_distype, slave_numdof>::SetSlaveState(
          const DRT::Discretization& slavedis,  ///< coupling slave discretization
          const std::vector<int>& lm            ///< local map
      )
      {
        // get state of the global vector
        Teuchos::RCP<const Epetra_Vector> matrix_state = slavedis.GetState(vel_statename_);
        if (matrix_state == Teuchos::null)
          dserror("Cannot get state vector %s", vel_statename_.c_str());

        // extract local values of the global vectors
        std::vector<double> mymatrix(lm.size());
        DRT::UTILS::ExtractMyValues(*matrix_state, mymatrix, lm);

        for (unsigned inode = 0; inode < slave_nen_; ++inode)  // number of nodes
        {
          for (unsigned idim = 0; idim < nsd_; ++idim)  // number of dimensions
          {
            (slave_vel_)(idim, inode) =
                mymatrix[idim +
                         (inode * slave_numdof)];  // state vector includes velocity and pressure
          }
          if (slave_numdof == nsd_ + 1)
            (slave_pres_)(inode, 0) = mymatrix[nsd_ + (inode * slave_numdof)];
        }

        return;
      }

      /*----------------------------------------------------------------------*
       *----------------------------------------------------------------------*/
      template <DRT::Element::DiscretizationType distype,
          DRT::Element::DiscretizationType slave_distype, unsigned int slave_numdof>
      void SlaveElementRepresentation<distype, slave_distype, slave_numdof>::SetSlaveStaten(
          const DRT::Discretization& slavedis,  ///< coupling slave discretization
          const std::vector<int>& lm            ///< local map
      )
      {
        // get state of the global vector
        Teuchos::RCP<const Epetra_Vector> matrix_state = slavedis.GetState(veln_statename_);
        if (matrix_state == Teuchos::null)
          dserror("Cannot get state vector %s", veln_statename_.c_str());

        // extract local values of the global vectors
        std::vector<double> mymatrix(lm.size());
        DRT::UTILS::ExtractMyValues(*matrix_state, mymatrix, lm);

        for (unsigned inode = 0; inode < slave_nen_; ++inode)  // number of nodes
        {
          for (unsigned idim = 0; idim < nsd_; ++idim)  // number of dimensions
          {
            (slave_veln_)(idim, inode) =
                mymatrix[idim +
                         (inode * slave_numdof)];  // state vector includes velocity and pressure
          }
          if (slave_numdof == nsd_ + 1)
            (slave_presn_)(inode, 0) = mymatrix[nsd_ + (inode * slave_numdof)];
        }

        return;
      }

      /*----------------------------------------------------------------------*
       *----------------------------------------------------------------------*/
      template <DRT::Element::DiscretizationType distype,
          DRT::Element::DiscretizationType slave_distype, unsigned int slave_numdof>
      void SlaveElementRepresentation<distype, slave_distype, slave_numdof>::GetInterfaceVelnp(
          LINALG::Matrix<nsd_, 1>& ivelint  ///< interface velocity at coupling slave side
          ) const
      {
        ivelint.Multiply(slave_vel_, slave_funct_);
      }

      /*----------------------------------------------------------------------*
       *----------------------------------------------------------------------*/
      template <DRT::Element::DiscretizationType distype,
          DRT::Element::DiscretizationType slave_distype, unsigned int slave_numdof>
      void SlaveElementRepresentation<distype, slave_distype, slave_numdof>::GetInterfaceVeln(
          LINALG::Matrix<nsd_, 1>& ivelintn  ///< interface velocity at coupling slave side
          ) const
      {
        ivelintn.Multiply(slave_veln_, slave_funct_);
      }

      /*----------------------------------------------------------------------*
       *----------------------------------------------------------------------*/
      template <DRT::Element::DiscretizationType distype,
          DRT::Element::DiscretizationType slave_distype, unsigned int slave_numdof>
      void SlaveElementRepresentation<distype, slave_distype, slave_numdof>::GetInterfacePresnp(
          double& ipres  ///< interface pressure at coupling slave side
          ) const
      {
        // pressure at current gauss-point
        ipres = slave_funct_.Dot(slave_pres_);
      }

      /*----------------------------------------------------------------------*
       *----------------------------------------------------------------------*/
      template <DRT::Element::DiscretizationType distype,
          DRT::Element::DiscretizationType slave_distype, unsigned int slave_numdof>
      void SlaveElementRepresentation<distype, slave_distype, slave_numdof>::GetInterfacePresn(
          double& ipresn  ///< interface pressure at coupling slave side
          ) const
      {
        // pressure at current gauss-point
        ipresn = slave_funct_.Dot(slave_presn_);
      }

      /*----------------------------------------------------------------------*
       *----------------------------------------------------------------------*/
      template <DRT::Element::DiscretizationType distype,
          DRT::Element::DiscretizationType slave_distype, unsigned int slave_numdof>
      void SlaveElementRepresentation<distype, slave_distype, slave_numdof>::GetInterfaceVelGradnp(
          LINALG::Matrix<nsd_, nsd_>&
              velgradint  ///< interface velocity gradients at coupling slave side
          ) const
      {
        velgradint = slave_vderxy_;
      }

      /*----------------------------------------------------------------------*
       *----------------------------------------------------------------------*/
      template <DRT::Element::DiscretizationType distype,
          DRT::Element::DiscretizationType slave_distype, unsigned int slave_numdof>
      void SlaveElementRepresentation<distype, slave_distype, slave_numdof>::GetInterfaceVelGradn(
          LINALG::Matrix<nsd_, nsd_>&
              velgradintn  ///< interface velocity gradients at coupling slave side
          ) const
      {
        velgradintn = slave_vderxyn_;
      }

      /*----------------------------------------------------------------------*
       *----------------------------------------------------------------------*/
      template <DRT::Element::DiscretizationType distype,
          DRT::Element::DiscretizationType slave_distype, unsigned int slave_numdof>
      void SlaveElementRepresentation<distype, slave_distype, slave_numdof>::GetSlaveFunct(
          LINALG::Matrix<slave_nen_, 1>& slave_funct  ///< coupling slave shape functions
          ) const
      {
        slave_funct = slave_funct_;
      }

      /*----------------------------------------------------------------------*
       *----------------------------------------------------------------------*/
      template <DRT::Element::DiscretizationType distype,
          DRT::Element::DiscretizationType slave_distype, unsigned int slave_numdof>
      void
      SlaveElementRepresentation<distype, slave_distype, slave_numdof>::SetInterfaceJumpStatenp(
          const DRT::Discretization& cutterdis,  ///< cutter discretization
          const std::string state,               ///< state
          const std::vector<int>& lm             ///< local map
      )
      {
        // get state of the global vector
        Teuchos::RCP<const Epetra_Vector> matrix_state = cutterdis.GetState(state);
        if (matrix_state == Teuchos::null) dserror("Cannot get state vector %s", state.c_str());

        // extract local values of the global vectors
        std::vector<double> mymatrix(lm.size());
        DRT::UTILS::ExtractMyValues(*matrix_state, mymatrix, lm);

        for (unsigned inode = 0; inode < slave_nen_; ++inode)  // number of nodes
        {
          for (unsigned idim = 0; idim < nsd_; ++idim)  // number of dimensions
          {
            (interface_velnp_jump_)(idim, inode) =
                mymatrix[idim +
                         (inode * slave_numdof)];  // state vector includes velocity and pressure
            // no pressure jump required
          }
        }

        return;
      }

      /*----------------------------------------------------------------------*
       *----------------------------------------------------------------------*/
      template <DRT::Element::DiscretizationType distype,
          DRT::Element::DiscretizationType slave_distype, unsigned int slave_numdof>
      void SlaveElementRepresentation<distype, slave_distype, slave_numdof>::SetInterfaceJumpStaten(
          const DRT::Discretization& cutterdis,  ///< cutter discretization
          const std::string state,               ///< state
          const std::vector<int>& lm             ///< local map
      )
      {
        // get state of the global vector
        Teuchos::RCP<const Epetra_Vector> matrix_state = cutterdis.GetState(state);
        if (matrix_state == Teuchos::null) dserror("Cannot get state vector %s", state.c_str());

        // extract local values of the global vectors
        std::vector<double> mymatrix(lm.size());
        DRT::UTILS::ExtractMyValues(*matrix_state, mymatrix, lm);

        for (unsigned inode = 0; inode < slave_nen_; ++inode)  // number of nodes
        {
          for (unsigned idim = 0; idim < nsd_; ++idim)  // number of dimensions
          {
            (interface_veln_jump_)(idim, inode) =
                mymatrix[idim +
                         (inode * slave_numdof)];  // state vector includes velocity and pressure
            // no pressure jump required
          }
        }

        return;
      }

      /*----------------------------------------------------------------------*
       *----------------------------------------------------------------------*/
      template <DRT::Element::DiscretizationType distype,
          DRT::Element::DiscretizationType slave_distype, unsigned int slave_numdof>
      void SlaveElementRepresentation<distype, slave_distype, slave_numdof>::GetInterfaceJumpVelnp(
          LINALG::Matrix<nsd_, 1>& ivelint_jump  ///< cutter element interface velocity jump or
                                                 ///< prescribed DBC at Gaussian point
          ) const
      {
        ivelint_jump.Multiply(interface_velnp_jump_, slave_funct_);
      }

      /*----------------------------------------------------------------------*
       *----------------------------------------------------------------------*/
      template <DRT::Element::DiscretizationType distype,
          DRT::Element::DiscretizationType slave_distype, unsigned int slave_numdof>
      void SlaveElementRepresentation<distype, slave_distype, slave_numdof>::GetInterfaceJumpVeln(
          LINALG::Matrix<nsd_, 1>& ivelintn_jump  ///< cutter element interface velocity jump or
                                                  ///< prescribed DBC at Gaussian point
          ) const
      {
        ivelintn_jump.Multiply(interface_veln_jump_, slave_funct_);
      }


      /*----------------------------------------------------------------------*
       *----------------------------------------------------------------------*/
      template <DRT::Element::DiscretizationType distype,
          DRT::Element::DiscretizationType slave_distype, unsigned int slave_numdof>
      void SlaveElementRepresentation<distype, slave_distype, slave_numdof>::Evaluate(
          LINALG::Matrix<nsd_, 1>& xslave)
      {
        LINALG::Matrix<3, 1> rst_slave(true);
        Evaluate(xslave, rst_slave);
      }

      /*----------------------------------------------------------------------*
       *----------------------------------------------------------------------*/
      template <DRT::Element::DiscretizationType distype,
          DRT::Element::DiscretizationType slave_distype, unsigned int slave_numdof>
      void SlaveElementRepresentation<distype, slave_distype, slave_numdof>::Evaluate(
          LINALG::Matrix<nsd_, 1>& xslave, LINALG::Matrix<nsd_, 1>& rst_slave)
      {
        // coupling with a 2D element
        if (slave_nsd_ == nsd_ - 1)
        {
          // evaluate shape function at solution
          DRT::UTILS::shape_function_2D(slave_funct_, xslave(0), xslave(1), slave_distype);
          rst_slave(0) = xslave(0);
          rst_slave(1) = xslave(1);
          //    dserror("You called 3D evaluation routine when coupling with a 2D element.");
          return;
        }

        // find element local position of gauss point
        Teuchos::RCP<GEO::CUT::Position> pos =
            GEO::CUT::PositionFactory::BuildPosition<nsd_, slave_distype>(slave_xyze_, xslave);
        pos->Compute();

        if (slave_nsd_ == nsd_)
        {
          pos->LocalCoordinates(rst_slave);
          DRT::UTILS::shape_function_3D(
              slave_funct_, rst_slave(0), rst_slave(1), rst_slave(2), slave_distype);
          DRT::UTILS::shape_function_3D_deriv1(
              slave_deriv_, rst_slave(0), rst_slave(1), rst_slave(2), slave_distype);
        }
        else
          dserror("Unsupported dimension clash!");


        LINALG::Matrix<nsd_, nsd_> slave_xjm(true);
        LINALG::Matrix<nsd_, nsd_> slave_xji(true);

        slave_xjm.MultiplyNT(slave_deriv_, slave_xyze_);
        slave_xji.Invert(slave_xjm);

        // compute global first derivates
        slave_derxy_.Multiply(slave_xji, slave_deriv_);

        // get velocity derivatives at integration point
        // (values at n+alpha_F for generalized-alpha scheme, n+1 otherwise)
        slave_vderxy_.MultiplyNT(slave_vel_, slave_derxy_);
        // previous time step
        slave_vderxyn_.MultiplyNT(slave_veln_, slave_derxy_);

        return;
      }

      /*----------------------------------------------------------------------*
       *----------------------------------------------------------------------*/
      template <DRT::Element::DiscretizationType distype,
          DRT::Element::DiscretizationType slave_distype, unsigned int slave_numdof>
      void SlaveElementRepresentation<distype, slave_distype, slave_numdof>::ComputeInterfaceForce(
          Epetra_SerialDenseVector& iforce,   ///< interface force vector
          LINALG::Matrix<nsd_, 1>& traction,  ///< traction vector at gaussian point
          const double& fac                   ///< integration factor
      )
      {
        for (unsigned inode = 0; inode < slave_nen_; ++inode)
        {
          double tmp = slave_funct_(inode) * fac;

          for (unsigned idim = 0; idim < nsd_; ++idim)
          {
            // f^i = ( N^i, t ) = ( N^i, (-pI+2mu*eps(u))*n )
            iforce[idim + inode * nsd_] += tmp * traction(idim);
          }
        }

        return;
      }

      /*----------------------------------------------------------------------*
       *----------------------------------------------------------------------*/
      template <DRT::Element::DiscretizationType distype,
          DRT::Element::DiscretizationType slave_distype, unsigned int slave_numdof>
      void SlaveElementRepresentation<distype, slave_distype, slave_numdof>::ProjectOnSide(
          LINALG::Matrix<nsd_, 1>&
              x_gp_lin,  ///< global coordinates of gaussian point w.r.t linearized interface
          LINALG::Matrix<nsd_, 1>& x_side,  ///< projected gaussian point on side
          LINALG::Matrix<nsd_, 1>&
              xi_side  ///< local coordinates of projected gaussian point w.r.t side
      )
      {
        //  TEUCHOS_FUNC_TIME_MONITOR( "FLD::XFluid::XFluidState::ProjectOnSide" );

#if (0)
        // REMARK: the current cut_kernel implementation is a factor of 3 slower than this
        // implementation

        GEO::CUT::Position2d<slave_distype> pos(slave_xyze_, x_gp_lin);
        pos.Compute(true);
        pos.LocalCoordinates(xi_side);

        // evaluate shape function at solution
        DRT::UTILS::shape_function_2D(slave_funct_, xi_side(0), xi_side(1), slave_distype);

        // get projected gauss point
        x_side.Multiply(slave_xyze_, slave_funct_);
#endif

        // check, if called on a 3D-element
#ifdef DEBUG
        if (slave_nsd_ == nsd_)
        {
          dserror(
              "You can't project onto a 3D coupling slave element directly. You need an associated "
              "boundary element!");
        }
#endif

        TEUCHOS_FUNC_TIME_MONITOR("FLD::XFluid::XFluidState::ProjectOnSide");


        if (slave_distype == DRT::Element::tri3 or slave_distype == DRT::Element::tri6)
        {
          proj_sol_(0) = 0.333333333333333;
          proj_sol_(1) = 0.333333333333333;
        }
        else if (slave_distype == DRT::Element::quad4 or slave_distype == DRT::Element::quad8 or
                 slave_distype == DRT::Element::quad9)
        {
          proj_sol_(0) = 0.0;
          proj_sol_(1) = 0.0;
        }
        else
        {
          dserror("Define start side xi-coordinates for unsupported cell type");
        }
        proj_sol_(2) = 0.0;

        // we only use absolute tolerances, since we compute local coordinates
        const double absTolIncr = 1.0e-9;  // abs tolerance for the local coordinates increment
        const double absTolRes = 1.0e-9;   // abs tolerance for the whole residual
        const double absTOLdist = 1.0e-9;  // abs tolerance for distance

        unsigned iter = 0;
        const unsigned maxiter = 7;

        bool converged = false;

        while (iter < maxiter && !converged)
        {
          ++iter;

          // get current values
          DRT::UTILS::shape_function_2D(proj_funct_, proj_sol_(0), proj_sol_(1), slave_distype);
          DRT::UTILS::shape_function_2D_deriv1(
              proj_deriv_, proj_sol_(0), proj_sol_(1), slave_distype);
          DRT::UTILS::shape_function_2D_deriv2(
              proj_deriv2_, proj_sol_(0), proj_sol_(1), slave_distype);

          proj_x_.Multiply(slave_xyze_, proj_funct_);
          proj_derxy_.MultiplyNT(slave_xyze_, proj_deriv_);
          proj_derxy2_.MultiplyNT(slave_xyze_, proj_deriv2_);

          proj_dx_drdr_times_dx_ds_(0) =
              proj_derxy2_(1, 0) * proj_derxy_(2, 1) - proj_derxy_(1, 1) * proj_derxy2_(2, 0);
          proj_dx_drdr_times_dx_ds_(1) =
              proj_derxy2_(2, 0) * proj_derxy_(0, 1) - proj_derxy_(2, 1) * proj_derxy2_(0, 0);
          proj_dx_drdr_times_dx_ds_(2) =
              proj_derxy2_(0, 0) * proj_derxy_(1, 1) - proj_derxy_(0, 1) * proj_derxy2_(1, 0);

          proj_dx_dr_times_dx_drds_(0) =
              proj_derxy_(1, 0) * proj_derxy2_(2, 2) - proj_derxy2_(1, 2) * proj_derxy_(2, 0);
          proj_dx_dr_times_dx_drds_(1) =
              proj_derxy_(2, 0) * proj_derxy2_(0, 2) - proj_derxy2_(2, 2) * proj_derxy_(0, 0);
          proj_dx_dr_times_dx_drds_(2) =
              proj_derxy_(0, 0) * proj_derxy2_(1, 2) - proj_derxy2_(0, 2) * proj_derxy_(1, 0);

          proj_dx_drds_times_dx_ds_(0) =
              proj_derxy2_(1, 2) * proj_derxy_(2, 1) - proj_derxy_(1, 1) * proj_derxy2_(2, 2);
          proj_dx_drds_times_dx_ds_(1) =
              proj_derxy2_(2, 2) * proj_derxy_(0, 1) - proj_derxy_(2, 1) * proj_derxy2_(0, 2);
          proj_dx_drds_times_dx_ds_(2) =
              proj_derxy2_(0, 2) * proj_derxy_(1, 1) - proj_derxy_(0, 1) * proj_derxy2_(1, 2);

          proj_dx_dr_times_dx_dsds_(0) =
              proj_derxy_(1, 0) * proj_derxy2_(2, 1) - proj_derxy2_(1, 1) * proj_derxy_(2, 0);
          proj_dx_dr_times_dx_dsds_(1) =
              proj_derxy_(2, 0) * proj_derxy2_(0, 1) - proj_derxy2_(2, 1) * proj_derxy_(0, 0);
          proj_dx_dr_times_dx_dsds_(2) =
              proj_derxy_(0, 0) * proj_derxy2_(1, 1) - proj_derxy2_(0, 1) * proj_derxy_(1, 0);

          proj_dx_dr_times_dx_ds_(0) =
              proj_derxy_(1, 0) * proj_derxy_(2, 1) - proj_derxy_(1, 1) * proj_derxy_(2, 0);
          proj_dx_dr_times_dx_ds_(1) =
              proj_derxy_(2, 0) * proj_derxy_(0, 1) - proj_derxy_(2, 1) * proj_derxy_(0, 0);
          proj_dx_dr_times_dx_ds_(2) =
              proj_derxy_(0, 0) * proj_derxy_(1, 1) - proj_derxy_(0, 1) * proj_derxy_(1, 0);

          // define sysmat
          for (unsigned i = 0; i < 3; ++i)
          {
            // d/dr
            proj_sysmat_(i, 0) =
                proj_derxy_(i, 0) -
                proj_sol_(2) * (proj_dx_drdr_times_dx_ds_(i) + proj_dx_dr_times_dx_drds_(i));

            // d/ds
            proj_sysmat_(i, 1) =
                proj_derxy_(i, 1) -
                proj_sol_(2) * (proj_dx_drds_times_dx_ds_(i) + proj_dx_dr_times_dx_dsds_(i));

            // d/d(dist)
            proj_sysmat_(i, 2) = -proj_dx_dr_times_dx_ds_(i);


            // residual
            proj_residuum_(i) =
                proj_x_(i) - proj_sol_(2) * proj_dx_dr_times_dx_ds_(i) - x_gp_lin(i);
          }

          proj_sysmat_.Invert();

          // solve Newton iteration
          proj_incr_.Multiply(
              -1.0, proj_sysmat_, proj_residuum_);  // incr = -Systemmatrix^-1 * residuum

          // update solution
          proj_sol_.Update(1.0, proj_incr_, 1.0);

          // check ° absolute criterion for local coordinates (between [-1,1]^2)
          //       ° absolute criterion for distance (-> 0)
          //       ° absolute criterion for whole residuum
          if (proj_incr_(2) < absTOLdist &&
              sqrt(proj_incr_(0) * proj_incr_(0) + proj_incr_(1) * proj_incr_(1)) < absTolIncr &&
              proj_residuum_.Norm2() < absTolRes)
          {
            converged = true;
          }
        }

        if (!converged)
        {
          std::cout.precision(15);

          std::cout << "increment criterion loc coord "
                    //<< sqrt(incr(0)*incr(0)+incr(1)*incr(1))/sqrt(sol(0)*sol(0)+sol(1)*sol(1))
                    << sqrt(proj_incr_(0) * proj_incr_(0) + proj_incr_(1) * proj_incr_(1))
                    << " \tabsTOL: " << absTolIncr << std::endl;
          std::cout << "absolute criterion for distance " << proj_incr_(2)
                    << " \tabsTOL: " << absTOLdist << std::endl;
          std::cout << "relative criterion whole residuum " << proj_residuum_.Norm2()
                    << " \tabsTOL: " << absTolRes << std::endl;


          std::cout << "sysmat.Invert" << proj_sysmat_ << std::endl;
          std::cout << "sol-norm " << proj_sol_.Norm2() << std::endl;
          std::cout << "sol " << proj_sol_ << std::endl;
          std::cout << "x_gp_lin" << x_gp_lin << std::endl;
          std::cout << "side " << slave_xyze_ << std::endl;

          dserror("Newton scheme in ProjectOnSide not converged! ");
        }

        // evaluate shape function at solution
        DRT::UTILS::shape_function_2D(slave_funct_, proj_sol_(0), proj_sol_(1), slave_distype);

        // get projected gauss point
        x_side.Multiply(slave_xyze_, slave_funct_);

        // set local coordinates w.r.t side
        xi_side(0) = proj_sol_(0);
        xi_side(1) = proj_sol_(1);
        xi_side(2) = 0.0;  // actually 2D coordinates
      }

      /*----------------------------------------------------------------------*
       *----------------------------------------------------------------------*/
      template <DRT::Element::DiscretizationType distype,
          DRT::Element::DiscretizationType slave_distype, unsigned int slave_numdof>
      double SlaveElementRepresentation<distype, slave_distype, slave_numdof>::EvalElementVolume()
      {
        switch (DRT::UTILS::DisTypeToDim<slave_distype>::dim)
        {
          case 3:
          {
            return XFEM::UTILS::EvalElementVolume<slave_distype>(slave_xyze_);
            break;
          }
          default:
          {
            dserror("Element volume for non 3D element type?");
            return 0.0;
          }
        }
      }

      /*----------------------------------------------------------------------*
       *----------------------------------------------------------------------*/
      template <DRT::Element::DiscretizationType distype,
          DRT::Element::DiscretizationType slave_distype, unsigned int slave_numdof>
      void SlaveElementRepresentation<distype, slave_distype, slave_numdof>::GetSlaveFunctDeriv(
          LINALG::Matrix<nsd_, slave_nen_>& slave_derxy) const
      {
        slave_derxy = slave_derxy_;
      }

      /*----------------------------------------------------------------------*
       *----------------------------------------------------------------------*/
      template <DRT::Element::DiscretizationType distype,
          DRT::Element::DiscretizationType slave_distype, unsigned int slave_numdof>
      NitscheCoupling<distype, slave_distype, slave_numdof>::NitscheCoupling(
          Epetra_SerialDenseMatrix& C_umum,  ///< C_umum coupling matrix
          Epetra_SerialDenseMatrix& rhC_um,  ///< C_um coupling rhs
          const DRT::ELEMENTS::FluidEleParameterXFEM&
              fldparaxfem  ///< specific XFEM based fluid parameters
          )
          : SlaveElementRepresentation<distype, slave_distype, slave_numdof>(),
            fldparaxfem_(fldparaxfem),
            C_umum_(C_umum.A(), true),
            rhC_um_(rhC_um.A(), true),
            adj_visc_scale_(fldparaxfem_.GetViscousAdjointScaling()),
            eval_coupling_(false)
      {
      }

      /*----------------------------------------------------------------------*
       *----------------------------------------------------------------------*/
      template <DRT::Element::DiscretizationType distype,
          DRT::Element::DiscretizationType slave_distype, unsigned int slave_numdof>
      NitscheCoupling<distype, slave_distype, slave_numdof>::NitscheCoupling(
          Epetra_SerialDenseMatrix& slave_xyze,  ///< global node coordinates of slave element
          Epetra_SerialDenseMatrix& C_umum,      ///< C_umum coupling matrix
          Epetra_SerialDenseMatrix& rhC_um,      ///< C_um coupling rhs
          const DRT::ELEMENTS::FluidEleParameterXFEM&
              fldparaxfem  ///< specific XFEM based fluid parameters
          )
          : SlaveElementRepresentation<distype, slave_distype, slave_numdof>(slave_xyze),
            fldparaxfem_(fldparaxfem),
            C_umum_(C_umum.A(), true),
            rhC_um_(rhC_um.A(), true),
            adj_visc_scale_(fldparaxfem_.GetViscousAdjointScaling()),
            eval_coupling_(false)
      {
      }

      /*----------------------------------------------------------------------*
       *----------------------------------------------------------------------*/
      template <DRT::Element::DiscretizationType distype,
          DRT::Element::DiscretizationType slave_distype, unsigned int slave_numdof>
      NitscheCoupling<distype, slave_distype, slave_numdof>::NitscheCoupling(
          Epetra_SerialDenseMatrix& slave_xyze,  ///< global node coordinates of slave element
          Epetra_SerialDenseMatrix& C_umum,      ///< C_umum coupling matrix
          Epetra_SerialDenseMatrix& C_usum,      ///< C_usum coupling matrix
          Epetra_SerialDenseMatrix& C_umus,      ///< C_umus coupling matrix
          Epetra_SerialDenseMatrix& C_usus,      ///< C_usus coupling matrix
          Epetra_SerialDenseMatrix& rhC_um,      ///< C_um coupling rhs
          Epetra_SerialDenseMatrix& rhC_us,      ///< C_us coupling rhs
          const DRT::ELEMENTS::FluidEleParameterXFEM&
              fldparaxfem  ///< specific XFEM based fluid parameters
          )
          : SlaveElementRepresentation<distype, slave_distype, slave_numdof>(slave_xyze),
            fldparaxfem_(fldparaxfem),
            C_umum_(C_umum.A(), true),
            C_usum_(C_usum.A(), true),
            C_umus_(C_umus.A(), true),
            C_usus_(C_usus.A(), true),
            rhC_um_(rhC_um.A(), true),
            rhC_us_(rhC_us.A(), true),
            adj_visc_scale_(fldparaxfem_.GetViscousAdjointScaling()),
            eval_coupling_(true)
      {
      }

      /*----------------------------------------------------------------------*
       *----------------------------------------------------------------------*/
      template <DRT::Element::DiscretizationType distype,
          DRT::Element::DiscretizationType slave_distype, unsigned int slave_numdof>
      void NitscheCoupling<distype, slave_distype, slave_numdof>::ApplyConvStabTerms(
          const Teuchos::RCP<SlaveElementInterface<distype>>&
              slave_ele,                            ///< associated slave element coupling object
          const LINALG::Matrix<nen_, 1>& funct_m,   ///< master shape functions
          const LINALG::Matrix<nsd_, 1>& velint_m,  ///< vector of slave shape functions
          const LINALG::Matrix<nsd_, 1>& normal,    ///< normal vector n^b
          const double& density_m,                  ///< fluid density (master)
          const double&
              NIT_stab_fac_conv,  ///< full Nitsche's penalty term scaling (viscous+convective part)
          const double& timefacfac,  ///< theta*dt
          const LINALG::Matrix<nsd_, 1>&
              ivelint_jump,  ///< prescribed interface velocity, Dirichlet values or jump height for
                             ///< coupled problems
          const INPAR::XFEM::EleCouplingCondType& cond_type  ///< condition type
      )
      {
        if (cond_type == INPAR::XFEM::CouplingCond_SURF_FLUIDFLUID &&
            fldparaxfem_.XffConvStabScaling() == INPAR::XFEM::XFF_ConvStabScaling_none)
          dserror("Cannot apply convective stabilization terms for XFF_ConvStabScaling_none!");

        // funct_m * timefac * fac * funct_m  * kappa_m (dyadic product)
        funct_m_m_dyad_.MultiplyNT(funct_m, funct_m);

        // velint_s
        velint_s_.Clear();

        if (eval_coupling_) slave_ele->GetInterfaceVelnp(velint_s_);

        // add the prescribed interface velocity for weak Dirichlet boundary conditions or the jump
        // height for coupled problems
        velint_s_.Update(1.0, ivelint_jump, 1.0);

        velint_diff_.Update(1.0, velint_m, -1.0, velint_s_, 0.0);


        // REMARK:
        // the (additional) convective stabilization is included in NIT_full_stab_fac
        // (in case of mixed/hybrid LM approaches, we don't compute the penalty term explicitly -
        // it 'evolves'); in that case we therefore don't chose the maximum, but add the penalty
        // term scaled with conv_stab_fac to the viscous counterpart; this happens by calling
        // NIT_Stab_Penalty

        switch (cond_type)
        {
          case INPAR::XFEM::CouplingCond_LEVELSET_WEAK_DIRICHLET:
          case INPAR::XFEM::CouplingCond_SURF_WEAK_DIRICHLET:
          case INPAR::XFEM::CouplingCond_SURF_FSI_PART:
          {
            NIT_Stab_Penalty(
                funct_m, timefacfac, std::pair<bool, double>(true, NIT_stab_fac_conv),  // F_Pen_Row
                std::pair<bool, double>(false, 0.0),                                    // X_Pen_Row
                std::pair<bool, double>(true, 1.0),                                     // F_Pen_Col
                std::pair<bool, double>(false, 0.0)                                     // X_Pen_Col
            );
            break;
          }
          case INPAR::XFEM::CouplingCond_SURF_FLUIDFLUID:
          {
            // funct_s
            Teuchos::RCP<SlaveElementRepresentation<distype, slave_distype, slave_numdof>> ser =
                Teuchos::rcp_dynamic_cast<
                    SlaveElementRepresentation<distype, slave_distype, slave_numdof>>(slave_ele);
            if (ser == Teuchos::null)
              dserror("Failed to cast slave_ele to SlaveElementRepresentation!");
            LINALG::Matrix<slave_nen_, 1> funct_s;
            ser->GetSlaveFunct(funct_s);

            // funct_s * timefac * fac * funct_s * kappa_s (dyadic product)
            funct_s_s_dyad_.MultiplyNT(funct_s_, funct_s);

            funct_s_m_dyad_.MultiplyNT(funct_s_, funct_m);

            if (fldparaxfem_.XffConvStabScaling() == INPAR::XFEM::XFF_ConvStabScaling_upwinding)
            {
              NIT_Stab_Penalty(funct_m, timefacfac,
                  std::pair<bool, double>(true, NIT_stab_fac_conv),  // F_Pen_Row
                  std::pair<bool, double>(true, NIT_stab_fac_conv),  // X_Pen_Row
                  std::pair<bool, double>(true, 1.0),                // F_Pen_Col
                  std::pair<bool, double>(true, 1.0)                 // X_Pen_Col
              );
            }

            // prevent instabilities due to convective mass transport across the fluid-fluid
            // interface
            if (fldparaxfem_.XffConvStabScaling() == INPAR::XFEM::XFF_ConvStabScaling_upwinding ||
                fldparaxfem_.XffConvStabScaling() == INPAR::XFEM::XFF_ConvStabScaling_only_averaged)
            {
              NIT_Stab_Inflow_AveragedTerm(funct_m, velint_m, normal, density_m, timefacfac);
            }
            break;
          }
          case INPAR::XFEM::CouplingCond_SURF_FSI_MONO:
          {
            dserror("Convective stabilization in monolithic XFSI is not yet available!");
            break;
          }
          default:
            dserror(
                "Unsupported coupling condition type. Cannot apply convective stabilization "
                "terms.");
            break;
        }
      }

      /*----------------------------------------------------------------------*
       *----------------------------------------------------------------------*/
      template <DRT::Element::DiscretizationType distype,
          DRT::Element::DiscretizationType slave_distype, unsigned int slave_numdof>
      void NitscheCoupling<distype, slave_distype, slave_numdof>::NIT_evaluateCoupling(
          const LINALG::Matrix<nsd_, 1>&
              normal,  ///< outward pointing normal (defined by the coupling partner, that
                       ///< determines the interface traction)
          const double& timefacfac,                ///< theta*dt*fac
          const double& pres_timefacfac,           ///< scaling for pressure part
          const double& visceff_m,                 ///< viscosity in coupling master fluid
          const double& visceff_s,                 ///< viscosity in coupling slave fluid
          const double& density_m,                 ///< fluid density (master) USED IN XFF
          const LINALG::Matrix<nen_, 1>& funct_m,  ///< coupling master shape functions
          const LINALG::Matrix<nsd_, nen_>&
              derxy_m,  ///< spatial derivatives of coupling master shape functions
          const LINALG::Matrix<nsd_, nsd_>&
              vderxy_m,          ///< coupling master spatial velocity derivatives
          const double& pres_m,  ///< coupling master pressure
          const LINALG::Matrix<nsd_, 1>& velint_m,  ///< coupling master interface velocity
          const LINALG::Matrix<nsd_, 1>&
              ivelint_jump,  ///< prescribed interface velocity, Dirichlet values or jump height for
                             ///< coupled problems
          const LINALG::Matrix<nsd_, 1>&
              itraction_jump,  ///< prescribed interface traction, jump height for coupled problems
          const LINALG::Matrix<nsd_, nsd_>& proj_tangential,  ///< tangential projection matrix
          const LINALG::Matrix<nsd_, nsd_>&
              LB_proj_matrix,  ///< prescribed projection matrix for laplace-beltrami problems
          const std::vector<Epetra_SerialDenseMatrix>&
              solid_stress,  ///< structural cauchy stress and linearization
          std::map<INPAR::XFEM::CoupTerm, std::pair<bool, double>>&
              configmap  ///< Interface Terms configuration map
      )
      {
        TEUCHOS_FUNC_TIME_MONITOR("FLD::NIT_evaluateCoupling");

        //--------------------------------------------

        // define the coupling between two not matching grids
        // for fluidfluidcoupling
        // domain Omega^m := Coupling master (XFluid)
        // domain Omega^s := Alefluid( or monolithic: structure) ( not available for non-coupling
        // (Dirichlet) )

        // [| v |] := vm - vs
        //  { v }  := kappa_m * vm + kappas * vs = kappam * vm (for Dirichlet coupling km=1.0, ks =
        //  0.0) < v >  := kappa_s * vm + kappam * vs = kappas * vm (for Dirichlet coupling km=1.0,
        //  ks = 0.0)
        //--------------------------------------------

        // Create projection matrices
        //--------------------------------------------
        proj_tangential_ = proj_tangential;
        proj_normal_.Scale(0.0);
        for (unsigned i = 0; i < nsd_; i++) proj_normal_(i, i) = 1.0;
        proj_normal_.Update(-1.0, proj_tangential_, 1.0);

        half_normal_.Update(0.5, normal, 0.0);
        normal_pres_timefacfac_.Update(pres_timefacfac, normal, 0.0);

        // get velocity at integration point
        // (values at n+alpha_F for generalized-alpha scheme, n+1 otherwise)
        // interface velocity vector in gausspoint
        velint_s_.Clear();

        if (configmap.at(INPAR::XFEM::X_Adj_Col).first ||
            configmap.at(INPAR::XFEM::X_Pen_Col).first ||
            configmap.at(INPAR::XFEM::X_Adj_n_Col).first ||
            configmap.at(INPAR::XFEM::X_Pen_n_Col).first ||
            configmap.at(INPAR::XFEM::X_Adj_t_Col).first ||
            configmap.at(INPAR::XFEM::X_Pen_t_Col).first)
          this->GetInterfaceVelnp(velint_s_);

        // Calc full veldiff
        if (configmap.at(INPAR::XFEM::F_Adj_Row).first ||
            configmap.at(INPAR::XFEM::XF_Adj_Row).first ||
            configmap.at(INPAR::XFEM::XS_Adj_Row).first ||
            configmap.at(INPAR::XFEM::F_Pen_Row).first ||
            configmap.at(INPAR::XFEM::X_Pen_Row).first)
        {
          velint_diff_.Update(configmap.at(INPAR::XFEM::F_Adj_Col).second, velint_m,
              -configmap.at(INPAR::XFEM::X_Adj_Col).second, velint_s_, 0.0);
          // add the prescribed interface velocity for weak Dirichlet boundary conditions or the
          // jump height for coupled problems
          velint_diff_.Update(-1.0, ivelint_jump, 1.0);

#ifdef PROJECT_VEL_FOR_PRESSURE_ADJOINT
          LINALG::Matrix<nsd_, 1> tmp_pval;
          tmp_pval.Multiply(proj_normal_, normal_pres_timefacfac_);
          // Project the velocity jump [|u|] in the pressure term with the projection matrix.
          //  Useful if smoothed normals are used (performs better for rotating cylinder case).
          velint_diff_pres_timefacfac_ = velint_diff_.Dot(tmp_pval);
#else
          velint_diff_pres_timefacfac_ = velint_diff_.Dot(normal_pres_timefacfac_);
#endif
        }

        // Calc normal-veldiff
        if (configmap.at(INPAR::XFEM::F_Adj_n_Row).first ||
            configmap.at(INPAR::XFEM::XF_Adj_n_Row).first ||
            configmap.at(INPAR::XFEM::XS_Adj_n_Row).first ||
            configmap.at(INPAR::XFEM::F_Pen_n_Row).first ||
            configmap.at(INPAR::XFEM::X_Pen_n_Row).first)
        {
          // velint_diff_proj_normal_ = (u^m_k - u^s_k - u^{jump}_k) P^n_{kj}
          // (([|u|]-u_0)*P^n) Apply from right for consistency
          velint_diff_normal_.Update(configmap.at(INPAR::XFEM::F_Pen_n_Col).second, velint_m,
              -configmap.at(INPAR::XFEM::X_Pen_n_Col).second, velint_s_, 0.0);
          // add the prescribed interface velocity for weak Dirichlet boundary conditions or the
          // jump height for coupled problems
          velint_diff_normal_.Update(-1.0, ivelint_jump, 1.0);
          velint_diff_proj_normal_.MultiplyTN(proj_normal_, velint_diff_normal_);

#ifdef PROJECT_VEL_FOR_PRESSURE_ADJOINT
          LINALG::Matrix<nsd_, 1> tmp_pval;
          tmp_pval.Multiply(proj_normal_, normal_pres_timefacfac_);
          // Project the velocity jump [|u|] in the pressure term with the projection matrix.
          //  Useful if smoothed normals are used (performs better for rotating cylinder case).
          velint_diff_normal_pres_timefacfac_ = velint_diff_normal_.Dot(tmp_pval);
#else
          velint_diff_normal_pres_timefacfac_ = velint_diff_normal_.Dot(normal_pres_timefacfac_);
#endif
        }

        // Calc tangential-veldiff
        if (configmap.at(INPAR::XFEM::F_Adj_t_Row).first ||
            configmap.at(INPAR::XFEM::XF_Adj_t_Row).first ||
            configmap.at(INPAR::XFEM::XS_Adj_t_Row).first ||
            configmap.at(INPAR::XFEM::F_Pen_t_Row).first ||
            configmap.at(INPAR::XFEM::X_Pen_t_Row).first)
        {
          // velint_diff_proj_tangential_ = (u^m_k - u^s_k - u^{jump}_k) P^t_{kj}
          // (([|u|]-u_0)*P^t) Apply from right for consistency
          velint_diff_tangential_.Update(configmap.at(INPAR::XFEM::F_Pen_t_Col).second, velint_m,
              -configmap.at(INPAR::XFEM::X_Pen_t_Col).second, velint_s_, 0.0);
          // add the prescribed interface velocity for weak Dirichlet boundary conditions or the
          // jump height for coupled problems
          velint_diff_tangential_.Update(-1.0, ivelint_jump, 1.0);
          velint_diff_proj_tangential_.MultiplyTN(proj_tangential_, velint_diff_tangential_);
        }

        // funct_s * timefac * fac
        funct_s_.Clear();
        if (slave_distype != DRT::Element::dis_none) this->GetSlaveFunct(funct_s_);

        // funct_m * timefac * fac * funct_m  * kappa_m (dyadic product)
        funct_m_m_dyad_.MultiplyNT(funct_m, funct_m);

        // funct_s * timefac * fac * funct_s * kappa_s (dyadic product)
        funct_s_s_dyad_.MultiplyNT(funct_s_, funct_s_);

        // funct_s * timefac * fac * funct_m (dyadic product)
        funct_s_m_dyad_.MultiplyNT(funct_s_, funct_m);

        //-----------------------------------------------------------------
        // viscous stability term
        // REMARK: this term includes also inflow coercivity in case of XFSI
        // with modified stabfac (see NIT_ComputeStabfac)

        if (configmap.at(INPAR::XFEM::F_Pen_n_Row).first ||
            configmap.at(INPAR::XFEM::X_Pen_n_Row).first)
        {
          // Normal Terms!
          NIT_Stab_Penalty_Projected(funct_m, proj_normal_, velint_diff_proj_normal_, timefacfac,
              configmap.at(INPAR::XFEM::F_Pen_n_Row), configmap.at(INPAR::XFEM::X_Pen_n_Row),
              configmap.at(INPAR::XFEM::F_Pen_n_Col), configmap.at(INPAR::XFEM::X_Pen_n_Col));
        }

        if (configmap.at(INPAR::XFEM::F_Pen_t_Row).first ||
            configmap.at(INPAR::XFEM::X_Pen_t_Row).first)
        {
          // Tangential Terms!
          NIT_Stab_Penalty_Projected(funct_m, proj_tangential_, velint_diff_proj_tangential_,
              timefacfac, configmap.at(INPAR::XFEM::F_Pen_t_Row),
              configmap.at(INPAR::XFEM::X_Pen_t_Row), configmap.at(INPAR::XFEM::F_Pen_t_Col),
              configmap.at(INPAR::XFEM::X_Pen_t_Col));
        }

        if (configmap.at(INPAR::XFEM::F_Pen_Row).first ||
            configmap.at(INPAR::XFEM::X_Pen_Row).first)
        {
          NIT_Stab_Penalty(funct_m, timefacfac, configmap.at(INPAR::XFEM::F_Pen_Row),
              configmap.at(INPAR::XFEM::X_Pen_Row), configmap.at(INPAR::XFEM::F_Pen_Col),
              configmap.at(INPAR::XFEM::X_Pen_Col));

          if (configmap.at(INPAR::XFEM::F_Pen_Row_linF1).first)
          {
            if (!configmap.at(INPAR::XFEM::F_Pen_Row_linF1).first ||
                !configmap.at(INPAR::XFEM::F_Pen_Row_linF2).first ||
                !configmap.at(INPAR::XFEM::F_Pen_Row_linF3).first)
              dserror("Linearization for Penalty Term not set for all Components!");

            NIT_Stab_Penalty_lin(funct_m, timefacfac, configmap.at(INPAR::XFEM::F_Pen_Row),
                configmap.at(INPAR::XFEM::F_Pen_Row_linF1),
                configmap.at(INPAR::XFEM::F_Pen_Row_linF2),
                configmap.at(INPAR::XFEM::F_Pen_Row_linF3));
          }
        }

        // add averaged term (TODO: For XFF? How does this work for non-master coupled? @Benedikt?)
        // Todo: is not handled by configmap yet as it has the shape of a penalty term and therefore
        // will be evaluated there at the end!
        if (fldparaxfem_.XffConvStabScaling() == INPAR::XFEM::XFF_ConvStabScaling_upwinding ||
            fldparaxfem_.XffConvStabScaling() == INPAR::XFEM::XFF_ConvStabScaling_only_averaged)
        {
          NIT_Stab_Inflow_AveragedTerm(funct_m, velint_m, normal, density_m, timefacfac);
        }
        //-------------------------------------------- Nitsche-Stab penalty added
        //-----------------------------------------------------------------

        // evaluate the terms, that contribute to the background fluid
        // system - standard Dirichlet case/pure xfluid-sided case
        // AND
        // system - two-sided or xfluid-sided:

        // 2 * mu_m * timefac * fac
        const double km_viscm_fac = 2.0 * timefacfac * visceff_m;
        half_normal_viscm_timefacfac_km_.Update(km_viscm_fac, half_normal_, 0.0);

        half_normal_deriv_m_viscm_timefacfac_km_.MultiplyTN(
            derxy_m, half_normal_);  // 0.5*normal(k)*derxy_m(k,ic)
        half_normal_deriv_m_viscm_timefacfac_km_.Scale(km_viscm_fac);

        // 0.5* (\nabla u + (\nabla u)^T) * normal
        vderxy_m_normal_.Multiply(vderxy_m, half_normal_);
        vderxy_m_normal_transposed_viscm_timefacfac_km_.MultiplyTN(vderxy_m, half_normal_);
        vderxy_m_normal_transposed_viscm_timefacfac_km_.Update(1.0, vderxy_m_normal_, 1.0);
        vderxy_m_normal_transposed_viscm_timefacfac_km_.Scale(km_viscm_fac);

        //-----------------------------------------------------------------
        // pressure consistency term
        if (configmap.at(INPAR::XFEM::F_Con_Col).first)
        {
          NIT_p_Consistency_MasterTerms(pres_m, funct_m, normal_pres_timefacfac_,
              configmap.at(INPAR::XFEM::F_Con_Row), configmap.at(INPAR::XFEM::X_Con_Row),
              configmap.at(INPAR::XFEM::F_Con_Col));
        }

        if (configmap.at(INPAR::XFEM::F_Con_n_Col)
                .first)  //(COMMENT: evaluating this seperatly seems to be more efficient for our
                         // cases)
        {
          NIT_p_Consistency_MasterTerms(pres_m, funct_m, normal_pres_timefacfac_,
              configmap.at(INPAR::XFEM::F_Con_n_Row), configmap.at(INPAR::XFEM::X_Con_n_Row),
              configmap.at(INPAR::XFEM::F_Con_n_Col));
        }

        //-----------------------------------------------------------------
        // viscous consistency term
        if (configmap.at(INPAR::XFEM::F_Con_Col).first)
        {
#ifndef ENFORCE_URQUIZA_GNBC
          // Comment: Here vderxy_m_normal_transposed_viscm_timefacfac_km_ is used!
          NIT_visc_Consistency_MasterTerms(derxy_m, funct_m, configmap.at(INPAR::XFEM::F_Con_Row),
              configmap.at(INPAR::XFEM::X_Con_Row), configmap.at(INPAR::XFEM::F_Con_Col));

#else  // Todo: @Magnus, what to do with this?
          NIT_visc_Consistency_MasterTerms_Projected(derxy_m, funct_m, proj_normal_, km_viscm_fac,
              1.0,   // not to chang anything in ENFORCE_URQUIZA_GNBC case
              1.0,   // not to chang anything in ENFORCE_URQUIZA_GNBC case
              1.0);  // not to chang anything in ENFORCE_URQUIZA_GNBC case
#endif
        }

        if (configmap.at(INPAR::XFEM::F_Con_n_Col).first)
        {
          NIT_visc_Consistency_MasterTerms_Projected(derxy_m, funct_m, proj_normal_, km_viscm_fac,
              configmap.at(INPAR::XFEM::F_Con_n_Row), configmap.at(INPAR::XFEM::X_Con_n_Row),
              configmap.at(INPAR::XFEM::F_Con_n_Col));
        }

        if (configmap.at(INPAR::XFEM::F_Con_t_Col).first)
        {
          NIT_visc_Consistency_MasterTerms_Projected(derxy_m, funct_m, proj_tangential_,
              km_viscm_fac, configmap.at(INPAR::XFEM::F_Con_t_Row),
              configmap.at(INPAR::XFEM::X_Con_t_Row), configmap.at(INPAR::XFEM::F_Con_t_Col));
        }

        //-----------------------------------------------------------------
        // pressure adjoint consistency term
        if (configmap.at(INPAR::XFEM::F_Adj_Row).first)
        {
          //-----------------------------------------------------------------
          // +++ qnuP option added! +++
          NIT_p_AdjointConsistency_MasterTerms(funct_m, normal_pres_timefacfac_,
              velint_diff_pres_timefacfac_, configmap.at(INPAR::XFEM::F_Adj_Row),
              configmap.at(INPAR::XFEM::F_Adj_Col), configmap.at(INPAR::XFEM::X_Adj_Col));
        }

        if (configmap.at(INPAR::XFEM::F_Adj_n_Row)
                .first)  //(COMMENT: evaluating this seperatly seems to be more efficient for our
                         // cases)
        {
          //-----------------------------------------------------------------
          // +++ qnuP option added! +++
          NIT_p_AdjointConsistency_MasterTerms(funct_m, normal_pres_timefacfac_,
              velint_diff_normal_pres_timefacfac_, configmap.at(INPAR::XFEM::F_Adj_n_Row),
              configmap.at(INPAR::XFEM::F_Adj_n_Col), configmap.at(INPAR::XFEM::X_Adj_n_Col));
        }

        //-----------------------------------------------------------------
        // viscous adjoint consistency term (and for NavierSlip Penalty Term ([v],{sigma}))
        // Normal Terms!

        if (configmap.at(INPAR::XFEM::F_Adj_n_Row).first)
        {
          Do_NIT_visc_Adjoint_and_Neumann_MasterTerms_Projected(funct_m,  ///< funct * timefacfac
              derxy_m,                   ///< spatial derivatives of coupling master shape functions
              vderxy_m,                  ///< coupling master spatial velocity derivatives
              proj_normal_,              ///< projection_matrix
              velint_diff_proj_normal_,  ///< velocity difference projected
              normal,                    ///< normal-vector
              km_viscm_fac, configmap.at(INPAR::XFEM::F_Adj_n_Row),
              configmap.at(INPAR::XFEM::F_Adj_n_Col), configmap.at(INPAR::XFEM::X_Adj_n_Col),
              configmap.at(INPAR::XFEM::FStr_Adj_n_Col));
        }
        if (configmap.at(INPAR::XFEM::FStr_Adj_n_Col).first)
          dserror("(NOT SUPPORTED FOR NORMAL DIR! Check Coercivity!)");

        // Tangential Terms!
        if (configmap.at(INPAR::XFEM::F_Adj_t_Row).first)
        {
          Do_NIT_visc_Adjoint_and_Neumann_MasterTerms_Projected(funct_m,  ///< funct * timefacfac
              derxy_m,           ///< spatial derivatives of coupling master shape functions
              vderxy_m,          ///< coupling master spatial velocity derivatives
              proj_tangential_,  ///< projection_matrix
              velint_diff_proj_tangential_,  ///< velocity difference projected
              normal,                        ///< normal-vector
              km_viscm_fac, configmap.at(INPAR::XFEM::F_Adj_t_Row),
              configmap.at(INPAR::XFEM::F_Adj_t_Col), configmap.at(INPAR::XFEM::X_Adj_t_Col),
              configmap.at(INPAR::XFEM::FStr_Adj_t_Col));
        }

        if (configmap.at(INPAR::XFEM::F_Adj_Row).first)
        {
          NIT_visc_AdjointConsistency_MasterTerms(funct_m,  ///< funct * timefacfac
              derxy_m,       ///< spatial derivatives of coupling master shape functions
              normal,        ///< normal-vector
              km_viscm_fac,  ///< scaling factor
              configmap.at(INPAR::XFEM::F_Adj_Row), configmap.at(INPAR::XFEM::F_Adj_Col),
              configmap.at(INPAR::XFEM::X_Adj_Col));

          if (configmap.at(INPAR::XFEM::FStr_Adj_Col).first)
            dserror("visc Adjoint Stress Term without projection not implemented - feel free!");
        }

        if ((configmap.at(INPAR::XFEM::XF_Con_Col).first ||
                configmap.at(INPAR::XFEM::XF_Con_n_Col).first ||
                configmap.at(INPAR::XFEM::XF_Con_t_Col).first ||
                configmap.at(INPAR::XFEM::XF_Adj_Row).first ||
                configmap.at(INPAR::XFEM::XF_Adj_n_Row).first ||
                configmap.at(INPAR::XFEM::XF_Adj_t_Row).first))
        {
          // TODO: @Christoph:
          //--------------------------------------------------------------------------------
          // This part needs to be adapted if a Robin-condition needs to be applied not only
          //  xfluid_sided (i.e. kappa^m =/= 1.0).
          // Should be more or less analogue to the above implementation.
          //--------------------------------------------------------------------------------

          //-----------------------------------------------------------------
          // the following quantities are only required for two-sided coupling
          // kappa_s > 0.0

          //-----------------------------------------------------------------
          // pressure consistency term

          double pres_s = 0.0;
          // must use this-pointer because of two-stage lookup!
          this->GetInterfacePresnp(pres_s);

          if (configmap.at(INPAR::XFEM::XF_Con_Col).first)
          {
            NIT_p_Consistency_SlaveTerms(pres_s, funct_m, normal_pres_timefacfac_,
                configmap.at(INPAR::XFEM::F_Con_Row), configmap.at(INPAR::XFEM::X_Con_Row),
                configmap.at(INPAR::XFEM::XF_Con_Col));
          }

          if (configmap.at(INPAR::XFEM::XF_Con_n_Col)
                  .first)  //(COMMENT: evaluating this seperatly seems to be more efficient for our
                           // cases)
          {
            NIT_p_Consistency_SlaveTerms(pres_s, funct_m, normal_pres_timefacfac_,
                configmap.at(INPAR::XFEM::F_Con_n_Row), configmap.at(INPAR::XFEM::X_Con_n_Row),
                configmap.at(INPAR::XFEM::XF_Con_n_Col));
          }

          //-----------------------------------------------------------------
          // pressure adjoint consistency term
          // HAS PROJECTION FOR VELOCITY IMPLEMENTED!!!
          if (configmap.at(INPAR::XFEM::XF_Adj_Row).first)
          {
            NIT_p_AdjointConsistency_SlaveTerms(normal_pres_timefacfac_,
                velint_diff_pres_timefacfac_, configmap.at(INPAR::XFEM::XF_Adj_Row),
                configmap.at(INPAR::XFEM::F_Adj_Col), configmap.at(INPAR::XFEM::X_Adj_Col));
          }
          if (configmap.at(INPAR::XFEM::XF_Adj_n_Row)
                  .first)  //(COMMENT: evaluating this seperatly seems to be more efficient for our
                           // cases)
          {
            NIT_p_AdjointConsistency_SlaveTerms(normal_pres_timefacfac_,
                velint_diff_normal_pres_timefacfac_, configmap.at(INPAR::XFEM::XF_Adj_n_Row),
                configmap.at(INPAR::XFEM::F_Adj_n_Col), configmap.at(INPAR::XFEM::X_Adj_n_Col));
          }

          //-----------------------------------------------------------------
          // viscous consistency term

          // Shape function derivatives for slave side
          LINALG::Matrix<nsd_, slave_nen_> derxy_s;
          this->GetSlaveFunctDeriv(derxy_s);

          // Spatial velocity gradient for slave side
          LINALG::Matrix<nsd_, nsd_> vderxy_s;
          this->GetInterfaceVelGradnp(vderxy_s);

          // 2 * mu_s * kappa_s * timefac * fac
          const double ks_viscs_fac = 2.0 * visceff_s * timefacfac;
          half_normal_viscs_timefacfac_ks_.Update(ks_viscs_fac, half_normal_, 0.0);
          half_normal_deriv_s_viscs_timefacfac_ks_.MultiplyTN(
              derxy_s, half_normal_);  // half_normal(k)*derxy_s(k,ic);
          half_normal_deriv_s_viscs_timefacfac_ks_.Scale(ks_viscs_fac);
          vderxy_s_normal_.Multiply(vderxy_s, half_normal_);
          vderxy_s_normal_transposed_viscs_timefacfac_ks_.MultiplyTN(vderxy_s, half_normal_);
          vderxy_s_normal_transposed_viscs_timefacfac_ks_.Update(1.0, vderxy_s_normal_, 1.0);
          vderxy_s_normal_transposed_viscs_timefacfac_ks_.Scale(ks_viscs_fac);

          if (configmap.at(INPAR::XFEM::XF_Con_Col).first)
          {
            NIT_visc_Consistency_SlaveTerms(derxy_s, funct_m, configmap.at(INPAR::XFEM::F_Con_Row),
                configmap.at(INPAR::XFEM::X_Con_Row), configmap.at(INPAR::XFEM::XF_Con_Col));
          }
          if (configmap.at(INPAR::XFEM::XF_Con_n_Col).first ||
              configmap.at(INPAR::XFEM::XF_Con_t_Col).first)
            dserror("Want to implement projected slave consistency?");

          //-----------------------------------------------------------------
          // viscous adjoint consistency term

          LINALG::Matrix<nsd_, slave_nen_> derxy_s_viscs_timefacfac_ks(derxy_s);
          derxy_s_viscs_timefacfac_ks.Scale(adj_visc_scale_ * ks_viscs_fac);

          // TODO: Needs added Projection. (If deemed necessary!)
          if (configmap.at(INPAR::XFEM::XF_Adj_Row).first)
          {
            NIT_visc_AdjointConsistency_SlaveTerms(funct_m, derxy_s_viscs_timefacfac_ks, normal,
                configmap.at(INPAR::XFEM::XF_Adj_Row), configmap.at(INPAR::XFEM::F_Adj_Col),
                configmap.at(INPAR::XFEM::X_Adj_Col));
          }
          if (configmap.at(INPAR::XFEM::XF_Adj_n_Row).first ||
              configmap.at(INPAR::XFEM::XF_Adj_t_Row).first)
            dserror("Want to  implement projected slave adjoint consistency?");

          //-----------------------------------------------------------------
          // standard consistency traction jump term
          // Only needed for XTPF
          if (configmap.at(INPAR::XFEM::F_TJ_Rhs).first ||
              configmap.at(INPAR::XFEM::X_TJ_Rhs).first)
          {
            // funct_s * timefac * fac * kappa_m
            funct_s_timefacfac_km_.Update(
                configmap.at(INPAR::XFEM::X_TJ_Rhs).second * timefacfac, funct_s_, 0.0);

            // funct_m * timefac * fac * kappa_s
            funct_m_timefacfac_ks_.Update(
                configmap.at(INPAR::XFEM::F_TJ_Rhs).second * timefacfac, funct_m, 0.0);

            NIT_Traction_Consistency_Term(
                funct_m_timefacfac_ks_, funct_s_timefacfac_km_, itraction_jump);
          }

          //-----------------------------------------------------------------
          // projection matrix approach (Laplace-Beltrami)
          if (configmap.at(INPAR::XFEM::F_LB_Rhs).first ||
              configmap.at(INPAR::XFEM::X_LB_Rhs).first)
          {
            LINALG::Matrix<nsd_, slave_nen_> derxy_s_timefacfac_km(derxy_s);
            derxy_s_timefacfac_km.Scale(configmap.at(INPAR::XFEM::X_LB_Rhs).second * timefacfac);

            LINALG::Matrix<nsd_, nen_> derxy_m_timefacfac_ks(derxy_m);
            derxy_m_timefacfac_ks.Scale(configmap.at(INPAR::XFEM::F_LB_Rhs).second * timefacfac);

            NIT_Projected_Traction_Consistency_Term(
                derxy_m_timefacfac_ks, derxy_s_timefacfac_km, LB_proj_matrix);
          }
          //-------------------------------------------- Traction-Jump added (XTPF)
        }

        // Structural Stress Terms (e.g. non xfluid sided FSI)
        if (configmap.at(INPAR::XFEM::XS_Con_Col).first ||
            configmap.at(INPAR::XFEM::XS_Con_n_Col).first ||
            configmap.at(INPAR::XFEM::XS_Con_t_Col).first ||
            configmap.at(INPAR::XFEM::XS_Adj_Row).first ||
            configmap.at(INPAR::XFEM::XS_Adj_n_Row).first ||
            configmap.at(INPAR::XFEM::XS_Adj_t_Row).first)
        {
          traction_ = LINALG::Matrix<nsd_, 1>(solid_stress[0].A(), true);
          dtraction_vel_ = LINALG::Matrix<nsd_ * slave_nen_, nsd_>(solid_stress[1].A(), true);

          d2traction_vel_[0] =
              LINALG::Matrix<nsd_ * slave_nen_, nsd_ * slave_nen_>(solid_stress[2].A(), true);
          d2traction_vel_[1] =
              LINALG::Matrix<nsd_ * slave_nen_, nsd_ * slave_nen_>(solid_stress[3].A(), true);
          d2traction_vel_[2] =
              LINALG::Matrix<nsd_ * slave_nen_, nsd_ * slave_nen_>(solid_stress[4].A(), true);

          if (configmap.at(INPAR::XFEM::XS_Con_Col).first)
          {
            NIT_solid_Consistency_SlaveTerms(funct_m, timefacfac,
                configmap.at(INPAR::XFEM::F_Con_Row), configmap.at(INPAR::XFEM::X_Con_Row),
                configmap.at(INPAR::XFEM::XS_Con_Col));
          }

          if (configmap.at(INPAR::XFEM::XS_Con_n_Col).first)
          {
            NIT_solid_Consistency_SlaveTerms_Projected(funct_m, proj_normal_, timefacfac,
                configmap.at(INPAR::XFEM::F_Con_n_Row), configmap.at(INPAR::XFEM::X_Con_n_Row),
                configmap.at(INPAR::XFEM::XS_Con_n_Col));
          }

          if (configmap.at(INPAR::XFEM::XS_Con_t_Col).first)
          {
            NIT_solid_Consistency_SlaveTerms_Projected(funct_m, proj_tangential_, timefacfac,
                configmap.at(INPAR::XFEM::F_Con_t_Row), configmap.at(INPAR::XFEM::X_Con_t_Row),
                configmap.at(INPAR::XFEM::XS_Con_t_Col));
          }

          if (configmap.at(INPAR::XFEM::XS_Adj_Row).first)
          {
            NIT_solid_AdjointConsistency_SlaveTerms(funct_m, timefacfac, velint_diff_,
                dtraction_vel_, configmap.at(INPAR::XFEM::XS_Adj_Row),
                configmap.at(INPAR::XFEM::F_Adj_Col), configmap.at(INPAR::XFEM::X_Adj_Col));
          }

          if (configmap.at(INPAR::XFEM::XS_Adj_n_Row).first)
          {
            NIT_solid_AdjointConsistency_SlaveTerms_Projected(funct_m, timefacfac, proj_normal_,
                velint_diff_proj_normal_, dtraction_vel_, configmap.at(INPAR::XFEM::XS_Adj_n_Row),
                configmap.at(INPAR::XFEM::F_Adj_n_Col), configmap.at(INPAR::XFEM::X_Adj_n_Col));
          }

          if (configmap.at(INPAR::XFEM::XS_Adj_t_Row).first)
          {
            NIT_solid_AdjointConsistency_SlaveTerms_Projected(funct_m, timefacfac, proj_tangential_,
                velint_diff_proj_tangential_, dtraction_vel_,
                configmap.at(INPAR::XFEM::XS_Adj_t_Row), configmap.at(INPAR::XFEM::F_Adj_t_Col),
                configmap.at(INPAR::XFEM::X_Adj_t_Col));
          }
        }

        return;
      }

      /*----------------------------------------------------------------------*
       *----------------------------------------------------------------------*/
      template <DRT::Element::DiscretizationType distype,
          DRT::Element::DiscretizationType slave_distype, unsigned int slave_numdof>
      void NitscheCoupling<distype, slave_distype, slave_numdof>::NIT_solid_Consistency_SlaveTerms(
          const LINALG::Matrix<nen_, 1>& funct_m,  ///< funct_m
          const double& timefacfac,                ///< theta*dt*fac
          const std::pair<bool, double>& m_row,    ///< scaling for master row
          const std::pair<bool, double>& s_row,    ///< scaling for slave row
          const std::pair<bool, double>& s_col,    ///< scaling for slave col
          bool only_rhs)
      {
        const double facms = m_row.second * s_col.second;
        const double facss = s_row.second * s_col.second;

        for (unsigned ir = 0; ir < nen_; ++ir)
        {
          const double tmp_val = funct_m(ir) * facms * timefacfac;
          for (unsigned ivel = 0; ivel < nsd_; ++ivel)
          {
            // rhs
            rhC_um_(mIndex(ir, ivel), 0) += tmp_val * traction_(ivel);
          }
        }

        for (unsigned ir = 0; ir < slave_nen_; ++ir)
        {
          const double tmp_val = funct_s_(ir) * facss * timefacfac;
          for (unsigned ivel = 0; ivel < nsd_; ++ivel)
          {
            // rhs
            rhC_us_(sIndex(ir, ivel), 0) -= tmp_val * traction_(ivel);
          }
        }

        if (only_rhs) return;

        for (unsigned ic = 0; ic < slave_nen_; ++ic)
        {
          for (unsigned jvel = 0; jvel < nsd_; ++jvel)
          {
            for (unsigned ivel = 0; ivel < nsd_; ++ivel)
            {
              const unsigned col = sIndex(ic, jvel);
              for (unsigned ir = 0; ir < nen_; ++ir)
              {
                //-----------------------------------------------
                //    - (vm,
                //-----------------------------------------------
                C_umus_(mIndex(ir, ivel), col) -=
                    funct_m(ir) * dtraction_vel_(col, ivel) * facms * timefacfac;
              }

              for (unsigned ir = 0; ir < slave_nen_; ++ir)
              {
                //-----------------------------------------------
                //    + (vs,
                //-----------------------------------------------
                // diagonal block
                C_usus_(sIndex(ir, ivel), col) +=
                    funct_s_(ir) * dtraction_vel_(col, ivel) * facss * timefacfac;
              }
            }
          }
        }

        return;
      }

      /*----------------------------------------------------------------------*
       *----------------------------------------------------------------------*/
      template <DRT::Element::DiscretizationType distype,
          DRT::Element::DiscretizationType slave_distype, unsigned int slave_numdof>
      void NitscheCoupling<distype, slave_distype,
          slave_numdof>::NIT_solid_Consistency_SlaveTerms_Projected(const LINALG::Matrix<nen_, 1>&
                                                                        funct_m,  ///< funct_m
          const LINALG::Matrix<nsd_, nsd_>& proj_matrix,  ///< projection matrix
          const double& timefacfac,                       ///< theta*dt*fac
          const std::pair<bool, double>& m_row,           ///< scaling for master row
          const std::pair<bool, double>& s_row,           ///< scaling for slave row
          const std::pair<bool, double>& s_col,           ///< scaling for slave col
          bool only_rhs)
      {
        static LINALG::Matrix<nsd_, 1> proj_traction;
        proj_traction.MultiplyTN(proj_matrix, traction_);

        const double facms = m_row.second * s_col.second;
        const double facss = s_row.second * s_col.second;

        for (unsigned ir = 0; ir < nen_; ++ir)
        {
          const double tmp_val = funct_m(ir) * facms * timefacfac;
          for (unsigned ivel = 0; ivel < nsd_; ++ivel)
          {
            // rhs
            rhC_um_(mIndex(ir, ivel), 0) += tmp_val * proj_traction(ivel);
          }
        }

        for (unsigned ir = 0; ir < slave_nen_; ++ir)
        {
          const double tmp_val = funct_s_(ir) * facss * timefacfac;
          for (unsigned ivel = 0; ivel < nsd_; ++ivel)
          {
            // rhs
            rhC_us_(sIndex(ir, ivel), 0) -= tmp_val * proj_traction(ivel);
          }
        }

        if (only_rhs) return;

        static LINALG::Matrix<nsd_ * slave_nen_, nsd_> proj_dtraction_vel(true);
        proj_dtraction_vel.Clear();
        for (unsigned col = 0; col < nsd_ * slave_nen_; ++col)
        {
          for (unsigned j = 0; j < nsd_; ++j)
          {
            for (unsigned i = 0; i < nsd_; ++i)
              proj_dtraction_vel(col, j) += dtraction_vel_(col, i) * proj_matrix(i, j);
          }
        }

        for (unsigned ic = 0; ic < slave_nen_; ++ic)
        {
          for (unsigned jvel = 0; jvel < nsd_; ++jvel)
          {
            for (unsigned ivel = 0; ivel < nsd_; ++ivel)
            {
              const unsigned col = sIndex(ic, jvel);
              for (unsigned ir = 0; ir < nen_; ++ir)
              {
                //-----------------------------------------------
                //    - (vm,
                //-----------------------------------------------
                C_umus_(mIndex(ir, ivel), col) -=
                    funct_m(ir) * proj_dtraction_vel(col, ivel) * facms * timefacfac;
              }

              for (unsigned ir = 0; ir < slave_nen_; ++ir)
              {
                //-----------------------------------------------
                //    + (vs,
                //-----------------------------------------------
                // diagonal block
                C_usus_(sIndex(ir, ivel), col) +=
                    funct_s_(ir) * proj_dtraction_vel(col, ivel) * facss * timefacfac;
              }
            }
          }
        }

        return;
      }

      /*----------------------------------------------------------------------*
       *----------------------------------------------------------------------*/
      template <DRT::Element::DiscretizationType distype,
          DRT::Element::DiscretizationType slave_distype, unsigned int slave_numdof>
      void NitscheCoupling<distype, slave_distype,
          slave_numdof>::NIT_solid_AdjointConsistency_SlaveTerms(const LINALG::Matrix<nen_, 1>&
                                                                     funct_m,  ///< funct_m
          const double& timefacfac,                                            ///< theta*dt*fac
          const LINALG::Matrix<nsd_, 1>& velint_diff,  ///< (velint_m - velint_s)
          const LINALG::Matrix<nsd_ * slave_nen_, nsd_>&
              dtraction_vel,                     ///< derivative of solid traction w.r.t. velocities
          const std::pair<bool, double>& s_row,  ///< scaling for slave row
          const std::pair<bool, double>& m_col,  ///< scaling for master col
          const std::pair<bool, double>& s_col,  ///< scaling for slave col
          bool only_rhs)
      {
        //
        // RHS: dv<d(sigma)/dv|u*n,uF-uS>
        // Lin: dv<d(sigma)/dv|u*n>duF-dv<d(sigma)/dv|u*n>duS+dv<d(sigma)/dv|du/dus*n,uF-uS>duS

        const double facs = s_row.second * timefacfac * adj_visc_scale_;
        for (unsigned ir = 0; ir < nen_; ++ir)
        {
          for (unsigned jvel = 0; jvel < nsd_; ++jvel)
          {
            for (unsigned ivel = 0; ivel < nsd_; ++ivel)
            {
              const unsigned row = sIndex(ir, ivel);
              // rhs
              rhC_us_(row, 0) += dtraction_vel(row, jvel) * velint_diff(jvel, 0) * facs;
            }
          }
        }

        if (only_rhs) return;

        const double facsm = s_row.second * m_col.second * timefacfac * adj_visc_scale_;
        const double facss = s_row.second * s_col.second * timefacfac * adj_visc_scale_;

        for (unsigned ir = 0; ir < slave_nen_; ++ir)
        {
          for (unsigned jvel = 0; jvel < nsd_; ++jvel)
          {
            for (unsigned ivel = 0; ivel < nsd_; ++ivel)
            {
              const unsigned row = sIndex(ir, ivel);
              for (unsigned ic = 0; ic < nen_; ++ic)
              {
                const unsigned col = mIndex(ic, jvel);
                C_usum_(row, col) -= funct_m(ic) * dtraction_vel(row, jvel) * facsm;
              }

              for (unsigned ic = 0; ic < slave_nen_; ++ic)
              {
                const unsigned col = sIndex(ic, jvel);
                C_usus_(row, col) += funct_s_(ic) * dtraction_vel(row, jvel) * facss;
                for (unsigned k = 0; k < nsd_; ++k)
                {
                  C_usus_(row, col) -= d2traction_vel_[k](row, col) * velint_diff(k, 0) * facs;
                }
              }
            }
          }
        }

        return;
      }

      /*----------------------------------------------------------------------*
       *----------------------------------------------------------------------*/
      template <DRT::Element::DiscretizationType distype,
          DRT::Element::DiscretizationType slave_distype, unsigned int slave_numdof>
      void NitscheCoupling<distype, slave_distype, slave_numdof>::
          NIT_solid_AdjointConsistency_SlaveTerms_Projected(
              const LINALG::Matrix<nen_, 1>& funct_m,           ///< funct_m
              const double& timefacfac,                         ///< theta*dt*fac
              const LINALG::Matrix<nsd_, nsd_>& proj_matrix,    ///< projection matrix
              const LINALG::Matrix<nsd_, 1>& proj_velint_diff,  ///< P^T*(velint_m - velint_s)
              const LINALG::Matrix<nsd_ * slave_nen_, nsd_>&
                  dtraction_vel,  ///< derivative of solid traction w.r.t. velocities
              const std::pair<bool, double>& s_row,  ///< scaling for slave row
              const std::pair<bool, double>& m_col,  ///< scaling for master col
              const std::pair<bool, double>& s_col,  ///< scaling for slave col
              bool only_rhs)
      {
        static LINALG::Matrix<nsd_ * slave_nen_, nsd_> proj_dtraction_vel(true);
        if (!only_rhs)
        {
          proj_dtraction_vel.Clear();
          for (unsigned col = 0; col < nsd_ * slave_nen_; ++col)
          {
            for (unsigned j = 0; j < nsd_; ++j)
            {
              for (unsigned i = 0; i < nsd_; ++i)
                proj_dtraction_vel(col, j) += dtraction_vel(col, i) * proj_matrix(i, j);
            }
          }
        }

        // Call Std Version with the projected quantities!
        NIT_solid_AdjointConsistency_SlaveTerms(funct_m, timefacfac, proj_velint_diff,
            proj_dtraction_vel, s_row, m_col, s_col, only_rhs);

        return;
      }

      /*----------------------------------------------------------------------*
       *----------------------------------------------------------------------*/
      template <DRT::Element::DiscretizationType distype,
          DRT::Element::DiscretizationType slave_distype, unsigned int slave_numdof>
      void NitscheCoupling<distype, slave_distype, slave_numdof>::NIT_evaluateCouplingOldState(
          const LINALG::Matrix<nsd_, 1>&
              normal,  ///< outward pointing normal (defined by the coupling partner, that
                       ///< determines the interface traction)
          const double& timefacfac,                ///< dt*(1-theta)*fac
          bool isImplPressure,                     ///< flag for implicit pressure treatment
          const double& visceff_m,                 ///< viscosity in coupling master fluid
          const double& visceff_s,                 ///< viscosity in coupling slave fluid
          const double& density_m,                 ///< fluid density (master) USED IN XFF
          const LINALG::Matrix<nen_, 1>& funct_m,  ///< coupling master shape functions
          const LINALG::Matrix<nsd_, nen_>&
              derxy_m,  ///< spatial derivatives of coupling master shape functions
          const LINALG::Matrix<nsd_, nsd_>&
              vderxy_m,          ///< coupling master spatial velocity derivatives
          const double& pres_m,  ///< coupling master pressure
          const LINALG::Matrix<nsd_, 1>& velint_m,  ///< coupling master interface velocity
          const LINALG::Matrix<nsd_, 1>&
              ivelint_jump,  ///< prescribed interface velocity, Dirichlet values or jump height for
                             ///< coupled problems
          const LINALG::Matrix<nsd_, nsd_>& proj_tangential,  ///< tangential projection matrix
          const LINALG::Matrix<nsd_, 1>&
              itraction_jump,  ///< prescribed interface traction, jump height for coupled problems
          std::map<INPAR::XFEM::CoupTerm, std::pair<bool, double>>&
              configmap  ///< Interface Terms configuration map
      )
      {
        //--------------------------------------------

        // define the coupling between two not matching grids
        // for fluidfluidcoupling
        // domain Omega^m := Coupling master (XFluid)
        // domain Omega^s := Alefluid( or monolithic: structure) ( not available for non-coupling
        // (Dirichlet) )

        // [| v |] := vm - vs
        //  { v }  := kappa_m * vm + kappas * vs = kappam * vm (for Dirichlet coupling km=1.0, ks =
        //  0.0) < v >  := kappa_s * vm + kappam * vs = kappas * vm (for Dirichlet coupling km=1.0,
        //  ks = 0.0)

        //--------------------------------------------

        // TODO: @XFEM-Team Add possibility to use new One Step Theta with Robin Boundary Condition.

        // Create projection matrices
        //--------------------------------------------
        proj_tangential_ = proj_tangential;
        proj_normal_.Scale(0.0);
        for (unsigned i = 0; i < nsd_; i++) proj_normal_(i, i) = 1.0;
        proj_normal_.Update(-1.0, proj_tangential_, 1.0);

        half_normal_.Update(0.5, normal, 0.0);
        normal_pres_timefacfac_.Update(timefacfac, normal, 0.0);

        // get velocity at integration point
        // (values at n)
        // interface velocity vector in gausspoint
        velint_s_.Clear();

        if (configmap.at(INPAR::XFEM::X_Adj_Col).first ||
            configmap.at(INPAR::XFEM::X_Pen_Col).first ||
            configmap.at(INPAR::XFEM::X_Adj_n_Col).first ||
            configmap.at(INPAR::XFEM::X_Pen_n_Col).first ||
            configmap.at(INPAR::XFEM::X_Adj_t_Col).first ||
            configmap.at(INPAR::XFEM::X_Pen_t_Col).first)
          this->GetInterfaceVeln(velint_s_);

        // Calc full veldiff
        if (configmap.at(INPAR::XFEM::F_Adj_Row).first ||
            configmap.at(INPAR::XFEM::XF_Adj_Row).first ||
            configmap.at(INPAR::XFEM::XS_Adj_Row).first ||
            configmap.at(INPAR::XFEM::F_Pen_Row).first ||
            configmap.at(INPAR::XFEM::X_Pen_Row).first)
        {
          velint_diff_.Update(configmap.at(INPAR::XFEM::F_Adj_Col).second, velint_m,
              -configmap.at(INPAR::XFEM::X_Adj_Col).second, velint_s_, 0.0);
          // add the prescribed interface velocity for weak Dirichlet boundary conditions or the
          // jump height for coupled problems
          velint_diff_.Update(-1.0, ivelint_jump, 1.0);

          //    #ifdef PROJECT_VEL_FOR_PRESSURE_ADJOINT //Todo: commented this not to change results
          //    with this commit
          //      LINALG::Matrix<nsd_,1> tmp_pval;
          //      tmp_pval.Multiply(proj_normal_,normal_pres_timefacfac_);
          //      //Project the velocity jump [|u|] in the pressure term with the projection matrix.
          //      //  Useful if smoothed normals are used (performs better for rotating cylinder
          //      case). velint_diff_pres_timefacfac_ = velint_diff_.Dot(tmp_pval);
          //    #else
          velint_diff_pres_timefacfac_ = velint_diff_.Dot(normal_pres_timefacfac_);
          //    #endif
        }

        // Calc normal-veldiff
        if (configmap.at(INPAR::XFEM::F_Adj_n_Row).first ||
            configmap.at(INPAR::XFEM::XF_Adj_n_Row).first ||
            configmap.at(INPAR::XFEM::XS_Adj_n_Row).first ||
            configmap.at(INPAR::XFEM::F_Pen_n_Row).first ||
            configmap.at(INPAR::XFEM::X_Pen_n_Row).first)
        {
          // velint_diff_proj_normal_ = (u^m_k - u^s_k - u^{jump}_k) P^n_{kj}
          // (([|u|]-u_0)*P^n) Apply from right for consistency
          velint_diff_normal_.Update(configmap.at(INPAR::XFEM::F_Adj_n_Col).second, velint_m,
              -configmap.at(INPAR::XFEM::X_Adj_n_Col).second, velint_s_, 0.0);
          // add the prescribed interface velocity for weak Dirichlet boundary conditions or the
          // jump height for coupled problems
          velint_diff_normal_.Update(-1.0, ivelint_jump, 1.0);
          velint_diff_proj_normal_.MultiplyTN(proj_normal_, velint_diff_normal_);

          //    #ifdef PROJECT_VEL_FOR_PRESSURE_ADJOINT //Todo: commented this not to change results
          //    with this commit
          //      LINALG::Matrix<nsd_,1> tmp_pval;
          //      tmp_pval.Multiply(proj_normal_,normal_pres_timefacfac_);
          //      //Project the velocity jump [|u|] in the pressure term with the projection matrix.
          //      //  Useful if smoothed normals are used (performs better for rotating cylinder
          //      case). velint_diff_normal_pres_timefacfac_ = velint_diff_normal_.Dot(tmp_pval);
          //    #else
          velint_diff_normal_pres_timefacfac_ = velint_diff_normal_.Dot(normal_pres_timefacfac_);
          //    #endif
        }

        // Calc tangential-veldiff
        if (configmap.at(INPAR::XFEM::F_Adj_t_Row).first ||
            configmap.at(INPAR::XFEM::XF_Adj_t_Row).first ||
            configmap.at(INPAR::XFEM::F_Pen_t_Row).first ||
            configmap.at(INPAR::XFEM::X_Pen_t_Row).first)
        {
          // velint_diff_proj_tangential_ = (u^m_k - u^s_k - u^{jump}_k) P^t_{kj}
          // (([|u|]-u_0)*P^t) Apply from right for consistency
          velint_diff_tangential_.Update(configmap.at(INPAR::XFEM::F_Adj_t_Col).second, velint_m,
              -configmap.at(INPAR::XFEM::X_Adj_t_Col).second, velint_s_, 0.0);
          // add the prescribed interface velocity for weak Dirichlet boundary conditions or the
          // jump height for coupled problems
          velint_diff_tangential_.Update(-1.0, ivelint_jump, 1.0);
          velint_diff_proj_tangential_.MultiplyTN(proj_tangential_, velint_diff_tangential_);
        }

        // funct_s * timefac * fac
        funct_s_.Clear();
        if (slave_distype != DRT::Element::dis_none) this->GetSlaveFunct(funct_s_);

        // funct_m * funct_m (dyadic product)
        funct_m_m_dyad_.MultiplyNT(funct_m, funct_m);

        // funct_s * funct_s (dyadic product)
        funct_s_s_dyad_.MultiplyNT(funct_s_, funct_s_);

        // funct_s * funct_m (dyadic product)
        funct_s_m_dyad_.MultiplyNT(funct_s_, funct_m);

        // penalty term
        if (fldparaxfem_.InterfaceTermsPreviousState() == INPAR::XFEM::PreviousState_full)
        {
          if (configmap.at(INPAR::XFEM::F_Pen_Row).first ||
              configmap.at(INPAR::XFEM::X_Pen_Row).first)
          {
            NIT_Stab_Penalty(funct_m, timefacfac, configmap.at(INPAR::XFEM::F_Pen_Row),
                configmap.at(INPAR::XFEM::X_Pen_Row), configmap.at(INPAR::XFEM::F_Pen_Col),
                configmap.at(INPAR::XFEM::X_Pen_Col), true);
          }

          // add averaged term
          if (fldparaxfem_.XffConvStabScaling() == INPAR::XFEM::XFF_ConvStabScaling_upwinding ||
              fldparaxfem_.XffConvStabScaling() == INPAR::XFEM::XFF_ConvStabScaling_only_averaged)
          {
            NIT_Stab_Inflow_AveragedTerm(funct_m, velint_m, normal, density_m, timefacfac, true);
          }
        }

        //-----------------------------------------------------------------
        // evaluate the terms, that contribute to the background fluid
        // system - standard Dirichlet case/pure xfluid-sided case
        // AND
        // system - two-sided or xfluid-sided:

        // 2 * mu_m * kappa_m * timefac * fac
        const double km_viscm_fac = 2.0 * timefacfac * visceff_m;
        half_normal_viscm_timefacfac_km_.Update(km_viscm_fac, half_normal_, 0.0);

        // 0.5* (\nabla u + (\nabla u)^T) * normal
        vderxy_m_normal_.Multiply(vderxy_m, half_normal_);
        vderxy_m_normal_transposed_viscm_timefacfac_km_.MultiplyTN(vderxy_m, half_normal_);
        vderxy_m_normal_transposed_viscm_timefacfac_km_.Update(1.0, vderxy_m_normal_, 1.0);
        vderxy_m_normal_transposed_viscm_timefacfac_km_.Scale(km_viscm_fac);

        // evaluate the terms, that contribute to the background fluid
        // system - two-sided or xfluid-sided:
        //-----------------------------------------------------------------
        // pressure consistency term
        if (not isImplPressure)
        {
          //-----------------------------------------------------------------
          // pressure consistency term
          if (configmap.at(INPAR::XFEM::F_Con_Col).first)
          {
            NIT_p_Consistency_MasterTerms(pres_m, funct_m, normal_pres_timefacfac_,
                configmap.at(INPAR::XFEM::F_Con_Row), configmap.at(INPAR::XFEM::X_Con_Row),
                configmap.at(INPAR::XFEM::F_Con_Col), true);
          }

          if (configmap.at(INPAR::XFEM::F_Con_n_Col)
                  .first)  //(COMMENT: evaluating this seperatly seems to be more efficient for our
                           // cases)
          {
            NIT_p_Consistency_MasterTerms(pres_m, funct_m, normal_pres_timefacfac_,
                configmap.at(INPAR::XFEM::F_Con_n_Row), configmap.at(INPAR::XFEM::X_Con_n_Row),
                configmap.at(INPAR::XFEM::F_Con_n_Col), true);
          }
        }

        //-----------------------------------------------------------------
        // viscous consistency term
        if (configmap.at(INPAR::XFEM::F_Con_Col).first)
        {
#ifndef ENFORCE_URQUIZA_GNBC
          const LINALG::Matrix<nsd_, nen_>
              dummy;  // as for the evaluation of the rhs this parameter is not used!
          // Comment: Here vderxy_m_normal_transposed_viscm_timefacfac_km_ is used!
          NIT_visc_Consistency_MasterTerms(dummy, funct_m, configmap.at(INPAR::XFEM::F_Con_Row),
              configmap.at(INPAR::XFEM::X_Con_Row), configmap.at(INPAR::XFEM::F_Con_Col), true);
#else
          dserror(
              "ENFORCE_URQUIZA_GNBC for NIT_visc_Consistency_MasterRHS?");  //@Magnus: What
                                                                            // should we do here?
#endif
        }
        if (configmap.at(INPAR::XFEM::F_Con_n_Col).first) dserror("F_Con_n_Col will come soon");
        if (configmap.at(INPAR::XFEM::F_Con_t_Col).first) dserror("F_Con_t_Col will come soon");

        if (fldparaxfem_.InterfaceTermsPreviousState() == INPAR::XFEM::PreviousState_full)
        {
          if (not isImplPressure)
          {
            //-----------------------------------------------------------------
            // pressure adjoint consistency term
            if (configmap.at(INPAR::XFEM::F_Adj_Row).first)
            {
              //-----------------------------------------------------------------
              // +++ qnuP option added! +++
              NIT_p_AdjointConsistency_MasterTerms(funct_m, normal_pres_timefacfac_,
                  velint_diff_pres_timefacfac_, configmap.at(INPAR::XFEM::F_Adj_Row),
                  configmap.at(INPAR::XFEM::F_Adj_Col), configmap.at(INPAR::XFEM::X_Adj_Col), true);
            }

            if (configmap.at(INPAR::XFEM::F_Adj_n_Row)
                    .first)  //(COMMENT: evaluating this seperatly seems to be more efficient for
                             // our cases)
            {
              //-----------------------------------------------------------------
              // +++ qnuP option added! +++
              NIT_p_AdjointConsistency_MasterTerms(funct_m, normal_pres_timefacfac_,
                  velint_diff_normal_pres_timefacfac_, configmap.at(INPAR::XFEM::F_Adj_n_Row),
                  configmap.at(INPAR::XFEM::F_Adj_n_Col), configmap.at(INPAR::XFEM::X_Adj_n_Col),
                  true);
            }
          }

          // Normal Terms!
          if (configmap.at(INPAR::XFEM::F_Adj_n_Row).first)
            dserror("Implement normal Adjoint Consistency term RHS for NEW OST !");
          if (configmap.at(INPAR::XFEM::FStr_Adj_n_Col).first)
            dserror("(NOT SUPPORTED FOR NORMAL DIR! Check Coercivity!)");
          if (configmap.at(INPAR::XFEM::F_Adj_t_Row).first)
            dserror("Implement tangential Adjoint Consistency term RHS for NEW OST !");

          //-----------------------------------------------------------------
          // viscous adjoint consistency term
          if (configmap.at(INPAR::XFEM::F_Adj_Row).first)
          {
            const LINALG::Matrix<nsd_, nen_>
                dummy;  // as for the evaluation of the rhs this parameter is not used!
            NIT_visc_AdjointConsistency_MasterTerms(funct_m,  ///< funct * timefacfac
                dummy,         ///< spatial derivatives of coupling master shape functions
                normal,        ///< normal-vector
                km_viscm_fac,  ///< scaling factor
                configmap.at(INPAR::XFEM::F_Adj_Row), configmap.at(INPAR::XFEM::F_Adj_Col),
                configmap.at(INPAR::XFEM::X_Adj_Col));
          }
        }

        //-----------------------------------------------------------------
        // the following quantities are only required for two-sided coupling
        // kappa_s > 0.0
        if ((configmap.at(INPAR::XFEM::XF_Con_Col).first ||
                configmap.at(INPAR::XFEM::XF_Con_n_Col).first ||
                configmap.at(INPAR::XFEM::XF_Con_t_Col).first ||
                configmap.at(INPAR::XFEM::XF_Adj_Row).first ||
                configmap.at(INPAR::XFEM::XF_Adj_n_Row).first ||
                configmap.at(INPAR::XFEM::XF_Adj_t_Row).first))
        {
          //-----------------------------------------------------------------
          // pressure consistency term
          if (configmap.at(INPAR::XFEM::XF_Con_Col).first ||
              configmap.at(INPAR::XFEM::XF_Con_n_Col).first)
          {
            if (not isImplPressure)
            {
              double presn_s = 0.0;
              // must use this-pointer because of two-stage lookup!
              this->GetInterfacePresn(presn_s);

              if (configmap.at(INPAR::XFEM::XF_Con_Col).first)
              {
                NIT_p_Consistency_SlaveTerms(presn_s, funct_m, normal_pres_timefacfac_,
                    configmap.at(INPAR::XFEM::F_Con_Row), configmap.at(INPAR::XFEM::X_Con_Row),
                    configmap.at(INPAR::XFEM::XF_Con_Col), true);
              }

              if (configmap.at(INPAR::XFEM::XF_Con_n_Col)
                      .first)  //(COMMENT: evaluating this seperatly seems to be more efficient for
                               // our cases)
              {
                NIT_p_Consistency_SlaveTerms(presn_s, funct_m, normal_pres_timefacfac_,
                    configmap.at(INPAR::XFEM::F_Con_n_Row), configmap.at(INPAR::XFEM::X_Con_n_Row),
                    configmap.at(INPAR::XFEM::XF_Con_n_Col), true);
              }
            }
          }

          //-----------------------------------------------------------------
          // viscous consistency term

          // Spatial velocity gradient for slave side
          LINALG::Matrix<nsd_, nsd_> vderxyn_s;
          this->GetInterfaceVelGradn(vderxyn_s);

          // 2 * mu_s * kappa_s * timefac * fac
          double ks_viscs_fac = 2.0 * visceff_s * timefacfac;

          vderxy_s_normal_.Multiply(vderxyn_s, half_normal_);
          vderxy_s_normal_transposed_viscs_timefacfac_ks_.MultiplyTN(vderxyn_s, half_normal_);
          vderxy_s_normal_transposed_viscs_timefacfac_ks_.Update(1.0, vderxy_s_normal_, 1.0);
          vderxy_s_normal_transposed_viscs_timefacfac_ks_.Scale(ks_viscs_fac);

          if (configmap.at(INPAR::XFEM::XF_Con_Col).first)
          {
            const LINALG::Matrix<nsd_, slave_nen_>
                dummy;  // as for the evaluation of the rhs this parameter is not used!
            NIT_visc_Consistency_SlaveTerms(dummy, funct_m, configmap.at(INPAR::XFEM::F_Con_Row),
                configmap.at(INPAR::XFEM::X_Con_Row), configmap.at(INPAR::XFEM::XF_Con_Col), true);
          }
          if (configmap.at(INPAR::XFEM::XF_Con_n_Col).first ||
              configmap.at(INPAR::XFEM::XF_Con_t_Col).first)
            dserror("Want to implement projected slave consistency?");

          // consistency terms evaluated
          if (fldparaxfem_.InterfaceTermsPreviousState() == INPAR::XFEM::PreviousState_full)
          {
            if (not isImplPressure)
            {
              //-----------------------------------------------------------------
              // pressure adjoint consistency term
              // HAS PROJECTION FOR VELOCITY IMPLEMENTED!!!
              if (configmap.at(INPAR::XFEM::XF_Adj_Row).first)
              {
                NIT_p_AdjointConsistency_SlaveTerms(normal_pres_timefacfac_,
                    velint_diff_pres_timefacfac_, configmap.at(INPAR::XFEM::XF_Adj_Row),
                    configmap.at(INPAR::XFEM::F_Adj_Col), configmap.at(INPAR::XFEM::X_Adj_Col),
                    true);
              }
              if (configmap.at(INPAR::XFEM::XF_Adj_n_Row)
                      .first)  //(COMMENT: evaluating this seperatly seems to be more efficient for
                               // our cases)
              {
                NIT_p_AdjointConsistency_SlaveTerms(normal_pres_timefacfac_,
                    velint_diff_normal_pres_timefacfac_, configmap.at(INPAR::XFEM::XF_Adj_n_Row),
                    configmap.at(INPAR::XFEM::F_Adj_n_Col), configmap.at(INPAR::XFEM::X_Adj_n_Col),
                    true);
              }
            }

            //-----------------------------------------------------------------
            // viscous adjoint consistency term
            // Shape function derivatives for slave side
            LINALG::Matrix<nsd_, slave_nen_> derxy_s;
            this->GetSlaveFunctDeriv(derxy_s);

            LINALG::Matrix<nsd_, slave_nen_> derxy_s_viscs_timefacfac_ks(derxy_s);
            derxy_s_viscs_timefacfac_ks.Scale(adj_visc_scale_ * ks_viscs_fac);

            // TODO: Needs added Projection. (If deemed necessary!)
            if (configmap.at(INPAR::XFEM::XF_Adj_Row).first)
            {
              NIT_visc_AdjointConsistency_SlaveTerms(funct_m, derxy_s_viscs_timefacfac_ks, normal,
                  configmap.at(INPAR::XFEM::XF_Adj_Row), configmap.at(INPAR::XFEM::F_Adj_Col),
                  configmap.at(INPAR::XFEM::X_Adj_Col), true);
            }
            if (configmap.at(INPAR::XFEM::XF_Adj_n_Row).first ||
                configmap.at(INPAR::XFEM::XF_Adj_t_Row).first)
              dserror("Want to  implement projected slave adjoint consistency?");
          }
        }

        //-----------------------------------------------------------------
        // standard consistency traction jump term
        // Only needed for XTPF
        if (configmap.at(INPAR::XFEM::F_TJ_Rhs).first || configmap.at(INPAR::XFEM::X_TJ_Rhs).first)
        {
          // funct_s * timefac * fac * kappa_m
          funct_s_timefacfac_km_.Update(
              configmap.at(INPAR::XFEM::F_TJ_Rhs).second * timefacfac, funct_s_, 0.0);

          // funct_m * timefac * fac * kappa_s
          funct_m_timefacfac_ks_.Update(
              configmap.at(INPAR::XFEM::X_TJ_Rhs).second * timefacfac, funct_m, 0.0);

          NIT_Traction_Consistency_Term(
              funct_m_timefacfac_ks_, funct_s_timefacfac_km_, itraction_jump);
        }

        //-----------------------------------------------------------------
        // projection matrix approach (Laplace-Beltrami)
        if (configmap.at(INPAR::XFEM::F_LB_Rhs).first || configmap.at(INPAR::XFEM::X_LB_Rhs).first)
        {
          dserror(
              "Check if we need the (Laplace-Beltrami) for the old timestep, "
              "then you should not forget to add the LB_proj_matrix as member to this function?");
          //    LINALG::Matrix<nsd_,slave_nen_> derxy_s_timefacfac_km(derxy_s);
          //    derxy_s_timefacfac_km.Scale(configmap.at(INPAR::XFEM::F_LB_Rhs).second*timefacfac);
          //
          //    LINALG::Matrix<nsd_,nen_> derxy_m_timefacfac_ks(derxy_m);
          //    derxy_m_timefacfac_ks.Scale(configmap.at(INPAR::XFEM::X_LB_Rhs).second*timefacfac);
          //
          //    NIT_Projected_Traction_Consistency_Term(
          //    derxy_m_timefacfac_ks,
          //    derxy_s_timefacfac_km,
          //    LB_proj_matrix);
        }

        return;
      }

      /*----------------------------------------------------------------------*
       *----------------------------------------------------------------------*/
      template <DRT::Element::DiscretizationType distype,
          DRT::Element::DiscretizationType slave_distype, unsigned int slave_numdof>
      void NitscheCoupling<distype, slave_distype, slave_numdof>::NIT_Traction_Consistency_Term(
          const LINALG::Matrix<nen_, 1>& funct_m_timefacfac_ks,  ///< funct * timefacfac *kappa_s
          const LINALG::Matrix<slave_nen_, 1>&
              funct_s_timefacfac_km,  ///< funct_s * timefacfac *kappa_m
          const LINALG::Matrix<nsd_, 1>&
              itraction_jump  ///< prescribed interface traction, jump height for coupled problems
      )
      {
        /*            \
     - |  < v >,   t   |   with t = [sigma * n]
        \             /     */

        // Two-Phase Flow:
        //
        //     t_{n+1}          [| sigma*n |] =      gamma * curv * n
        //                                                             with curv =
        //                                                             div(grad(phi)/||grad(phi)||)
        //
        //      t_{n}           [| sigma*n |] = [|  -pI + \mu*[\nabla u + (\nabla u)^T]  |] * n
        //

        // Combustion:        TO BE IMPLEMENTED.
        //

        // All else:            [| sigma*n |] = 0


        // loop over velocity components
        for (unsigned ivel = 0; ivel < nsd_; ++ivel)
        {
          //-----------------------------------------------
          //    - (vm, ks * t)
          //-----------------------------------------------
          for (unsigned ir = 0; ir < nen_; ++ir)
          {
            const double funct_m_ks_timefacfac_traction =
                funct_m_timefacfac_ks(ir) * itraction_jump(ivel);

            const unsigned row = mIndex(ir, ivel);
            rhC_um_(row, 0) += funct_m_ks_timefacfac_traction;
          }

          //-----------------------------------------------
          //    + (vs, km * t)
          //-----------------------------------------------
          for (unsigned ir = 0; ir < slave_nen_; ++ir)
          {
            const double funct_s_km_timefacfac_traction =
                funct_s_timefacfac_km(ir) * itraction_jump(ivel);

            const unsigned row = sIndex(ir, ivel);
            rhC_us_(row, 0) += funct_s_km_timefacfac_traction;
          }
        }  // end loop over velocity components
      }

      /*----------------------------------------------------------------------*
       *----------------------------------------------------------------------*/
      template <DRT::Element::DiscretizationType distype,
          DRT::Element::DiscretizationType slave_distype, unsigned int slave_numdof>
      void NitscheCoupling<distype, slave_distype, slave_numdof>::
          NIT_Projected_Traction_Consistency_Term(
              const LINALG::Matrix<nsd_, nen_>&
                  derxy_m_timefacfac_ks,  ///< master shape function derivatives * timefacfac *
                                          ///< kappa_s
              const LINALG::Matrix<nsd_, slave_nen_>&
                  derxy_s_timefacfac_km,  ///< slave shape function derivatives * timefacfac *
                                          ///< kappa_m
              const LINALG::Matrix<nsd_, nsd_>&
                  itraction_jump_matrix  ///< prescribed projection matrix
          )
      {
        /*                   \
     - |  < \nabla v > : P    |   with P = (I - n (x) n )
        \                    /     */

        // Two-Phase Flow:
        //
        //     t_{n+1}          ( < \nabla v > : P )
        //                                          P can be calculated in different ways.
        //                                          P_smooth*P_cut is best approach so far.
        //
        //      t_{n}           [| sigma*n |] = [|  -pI + \mu*[\nabla u + (\nabla u)^T]  |] * n


        // Two-Phase Flow, Laplace Beltrami approach:
        // loop over velocity components
        for (unsigned ivel = 0; ivel < nsd_; ++ivel)
        {
          //-----------------------------------------------
          //    - (\nabla vm, ks * P)
          //-----------------------------------------------
          for (unsigned ir = 0; ir < nen_; ++ir)
          {
            double derxy_m_ks_timefacfac_sum = 0.0;
            // Sum over space dimensions
            for (unsigned idum = 0; idum < nsd_; ++idum)
            {
              derxy_m_ks_timefacfac_sum +=
                  derxy_m_timefacfac_ks(idum, ir) * itraction_jump_matrix(idum, ivel);
            }

            const unsigned row = mIndex(ir, ivel);
            rhC_um_(row, 0) -= derxy_m_ks_timefacfac_sum;
          }

          //-----------------------------------------------
          //    + (\nabla vs, km * P)
          //-----------------------------------------------
          for (unsigned ir = 0; ir < slave_nen_; ++ir)
          {
            double derxy_s_km_timefacfac_sum = 0.0;
            // Sum over space dimensions
            for (unsigned idum = 0; idum < nsd_; ++idum)
            {
              derxy_s_km_timefacfac_sum +=
                  derxy_s_timefacfac_km(idum, ir) * itraction_jump_matrix(idum, ivel);
            }

            const unsigned row = sIndex(ir, ivel);
            rhC_us_(row, 0) -= derxy_s_km_timefacfac_sum;
          }
        }  // end loop over velocity components
      }



      /*----------------------------------------------------------------------*
       *----------------------------------------------------------------------*/
      template <DRT::Element::DiscretizationType distype,
          DRT::Element::DiscretizationType slave_distype, unsigned int slave_numdof>
      void NitscheCoupling<distype, slave_distype, slave_numdof>::NIT_p_Consistency_MasterTerms(
          const double& pres_m,                              ///< master pressure
          const LINALG::Matrix<nen_, 1>& funct_m,            ///< funct
          const LINALG::Matrix<nsd_, 1>& normal_timefacfac,  ///< normal vector * timefacfac
          const std::pair<bool, double>& m_row,              ///< scaling for master row
          const std::pair<bool, double>& s_row,              ///< scaling for slave row
          const std::pair<bool, double>& m_col,              ///< scaling for master col
          bool only_rhs)
      {
        TEUCHOS_FUNC_TIME_MONITOR("FLD::NIT_p_Consistency_MasterTerms");


        /*                   \       /          i      \
      + |  [ v ],   {Dp}*n    | = - | [ v ], { p }* n   |
        \                    /       \                */

        //-----------------------------------------------
        //    + (vm, km *(Dpm)*n)
        //-----------------------------------------------
        const double facmm = m_row.second * m_col.second;

        for (unsigned ir = 0; ir < nen_; ++ir)
        {
          const double funct_m_pres = funct_m(ir) * pres_m * facmm;

          // loop over velocity components
          for (unsigned ivel = 0; ivel < nsd_; ++ivel)
          {
            // -(v,p*n)
            rhC_um_(mIndex(ir, ivel), 0) -= funct_m_pres * normal_timefacfac(ivel);
          }
        }

        double facsm = 0.0;
        if (s_row.first)
        {
          facsm = s_row.second * m_col.second;
          for (unsigned ir = 0; ir < slave_nen_; ++ir)
          {
            const double funct_s_pres = funct_s_(ir) * pres_m * facsm;

            for (unsigned ivel = 0; ivel < nsd_; ++ivel)
            {
              // -(v,p*n)
              rhC_us_(sIndex(ir, ivel), 0) += funct_s_pres * normal_timefacfac(ivel);
            }
          }  // end loop over velocity components
        }

        if (only_rhs) return;

        for (unsigned ic = 0; ic < nen_; ++ic)
        {
          const unsigned col = mPres(ic);

          for (unsigned ir = 0; ir < nen_; ++ir)
          {
            const double tmp = funct_m_m_dyad_(ir, ic) * facmm;

            // loop over velocity components
            for (unsigned ivel = 0; ivel < nsd_; ++ivel)
            {
              // (v,Dp*n)
              C_umum_(mIndex(ir, ivel), col) += tmp * normal_timefacfac(ivel);
            }
          }
        }

        if (s_row.first)
        {
          for (unsigned ic = 0; ic < nen_; ++ic)
          {
            const unsigned col = mPres(ic);

            for (unsigned ir = 0; ir < slave_nen_; ++ir)
            {
              const double tmp = funct_s_m_dyad_(ir, ic) * facsm;

              for (unsigned ivel = 0; ivel < nsd_; ++ivel)
              {
                //-----------------------------------------------
                //    - (vs, km *(Dpm)*n)
                //-----------------------------------------------

                // (v,Dp*n)
                C_usum_(sIndex(ir, ivel), col) -= tmp * normal_timefacfac(ivel);
              }
            }
          }
        }

        return;
      }

      /*----------------------------------------------------------------------*
       *----------------------------------------------------------------------*/
      template <DRT::Element::DiscretizationType distype,
          DRT::Element::DiscretizationType slave_distype, unsigned int slave_numdof>
      void NitscheCoupling<distype, slave_distype, slave_numdof>::NIT_p_Consistency_SlaveTerms(
          const double& pres_s,                              ///< slave pressure
          const LINALG::Matrix<nen_, 1>& funct_m,            ///< funct
          const LINALG::Matrix<nsd_, 1>& normal_timefacfac,  ///< normal vector * timefacfac
          const std::pair<bool, double>& m_row,              ///< scaling for master row
          const std::pair<bool, double>& s_row,              ///< scaling for slave row
          const std::pair<bool, double>& s_col,              ///< scaling for slave col
          bool only_rhs)
      {
        const double facms = m_row.second * s_col.second;
        const double facss = s_row.second * s_col.second;

        for (unsigned ir = 0; ir < nen_; ++ir)
        {
          const double tmp = funct_m(ir) * pres_s * facms;
          // loop over velocity components
          for (unsigned ivel = 0; ivel < nsd_; ++ivel)
          {
            // -(vm, ks * ps*n)
            rhC_um_(mIndex(ir, ivel), 0) -= tmp * normal_timefacfac(ivel);
          }
        }

        for (unsigned ir = 0; ir < slave_nen_; ++ir)
        {
          const double tmp = funct_s_(ir) * pres_s * facss;
          for (unsigned ivel = 0; ivel < nsd_; ++ivel)
          {
            // +(vs,ks * ps*n)
            rhC_us_(sIndex(ir, ivel), 0) += tmp * normal_timefacfac(ivel);
          }
        }  // end loop over velocity components

        if (only_rhs) return;

        for (unsigned ic = 0; ic < slave_nen_; ++ic)
        {
          //-----------------------------------------------
          //    + (vm, ks *(Dps)*n)
          //-----------------------------------------------
          unsigned col = sPres(ic);

          for (unsigned ir = 0; ir < nen_; ++ir)
          {
            const double tmp = funct_s_m_dyad_(ic, ir) * facms;
            for (unsigned ivel = 0; ivel < nsd_; ++ivel)
            {
              // (vm, ks * Dps*n)
              C_umus_(mIndex(ir, ivel), col) += tmp * normal_timefacfac(ivel);
            }
          }

          //-----------------------------------------------
          //    - (vs, ks *(Dps)*n)
          //-----------------------------------------------
          for (unsigned ir = 0; ir < slave_nen_; ++ir)
          {
            const double tmp = funct_s_s_dyad_(ir, ic) * facss;
            for (unsigned ivel = 0; ivel < nsd_; ++ivel)
            {
              C_usus_(sIndex(ir, ivel), col) -= tmp * normal_timefacfac(ivel);
            }
          }
        }

        return;
      }

      /*----------------------------------------------------------------------*
       *----------------------------------------------------------------------*/
      template <DRT::Element::DiscretizationType distype,
          DRT::Element::DiscretizationType slave_distype, unsigned int slave_numdof>
      void
      NitscheCoupling<distype, slave_distype, slave_numdof>::NIT_p_AdjointConsistency_MasterTerms(
          const LINALG::Matrix<nen_, 1>& funct_m,            ///< funct
          const LINALG::Matrix<nsd_, 1>& normal_timefacfac,  ///< normal vector * timefacfac
          const double&
              velint_diff_normal_timefacfac,     ///< (velint_m - velint_s) * normal * timefacfac
          const std::pair<bool, double>& m_row,  ///< scaling for master row
          const std::pair<bool, double>& m_col,  ///< scaling for master col
          const std::pair<bool, double>& s_col,  ///< scaling for slave col
          bool only_rhs)
      {
        TEUCHOS_FUNC_TIME_MONITOR("FLD::NIT_p_AdjointConsistency_MasterTerms");

        // 1) No-split no qunP option:
        /*                   \     /              i   \
      - |  { q }*n ,[ Du ]     | = |  { q }*n  ,[ u ]   |
        \                    /     \                 */

        // 2) qunP option:
        /*                       \     /              i      \
      - |  { q }*n ,[ Du ] P^n    | = |  { q }*n  ,[ u ] P^n  |
        \                        /     \                    */

        // REMARK:
        // the sign of the pressure adjoint consistency term is opposite to the sign of the pressure
        // consistency term (interface), as a non-symmetric pressure formulation is chosen in the
        // standard fluid the sign of the standard volumetric pressure consistency term is opposite
        // to the (chosen) sign of the pressure-weighted continuity residual; think about the
        // Schur-complement for the Stokes problem: S_pp = A_pp - A_pu A_uu^-1 A_up
        // (--> A_pu == -A_up^T; sgn(A_pp) == sgn(- A_pu A_uu^-1 Aup), where A_pp comes from
        // pressure-stabilizing terms) a symmetric adjoint pressure consistency term would also
        // affect the sign of the pressure stabilizing terms for Stokes' problem, this sign choice
        // leads to a symmetric, positive definite Schur-complement matrix S (v, p*n)--> A_up;
        // -(q,u*n)--> -A_up^T; S_pp = A_pp + A_up^T A_uu A_up


        const double velint_diff_normal_timefacfac_km =
            velint_diff_normal_timefacfac * m_row.second;
        for (unsigned ir = 0; ir < nen_; ++ir)
        {
          // (qm*n, km * um)
          // -(qm*n,km * u_DBC) for weak DBC or
          // -(qm*n,km * us)
          rhC_um_(mPres(ir), 0) += funct_m(ir) * velint_diff_normal_timefacfac_km;
          //    rhC_um_(mPres(ir),0) -= funct_m(ir)*velint_diff_normal_timefacfac; //TESTING!
        }

        if (only_rhs) return;

#ifdef PROJECT_VEL_FOR_PRESSURE_ADJOINT
        LINALG::Matrix<nsd_, 1> proj_norm_timefacfac;
        proj_norm_timefacfac.Multiply(proj_normal_, normal_timefacfac);
#endif
        //-----------------------------------------------
        //    - (qm*n, km *(Dum))
        //-----------------------------------------------
        const double facmm = m_row.second * m_col.second;
        for (unsigned ic = 0; ic < nen_; ++ic)
        {
          for (unsigned ir = 0; ir < nen_; ++ir)
          {
            const unsigned row = mPres(ir);

            const double tmp = funct_m_m_dyad_(ir, ic) * facmm;

            for (unsigned ivel = 0; ivel < nsd_; ++ivel)
            {
              // - (qm*n, km *(Dum))
#ifndef PROJECT_VEL_FOR_PRESSURE_ADJOINT
              C_umum_(row, mIndex(ic, ivel)) -= tmp * normal_timefacfac(ivel);
              //        C_umum_(row, mIndex(ic,ivel)) += tmp*normal_timefacfac_km(ivel); //TESTING!
#else
              C_umum_(row, mIndex(ic, ivel)) -= tmp * proj_norm_timefacfac(ivel);
#endif
            }
          }
        }

        if (s_col.first)
        {
          const double facms = m_row.second * s_col.second;
          //-----------------------------------------------
          //    + (qm*n, km *(Dus))
          //-----------------------------------------------
          for (unsigned ic = 0; ic < slave_nen_; ++ic)
          {
            for (unsigned ir = 0; ir < nen_; ++ir)
            {
              const unsigned row = mPres(ir);
              const double tmp = funct_s_m_dyad_(ic, ir) * facms;

              for (unsigned ivel = 0; ivel < nsd_; ++ivel)
              {
                // -(qm*n, km * Dus)
#ifndef PROJECT_VEL_FOR_PRESSURE_ADJOINT
                C_umus_(row, sIndex(ic, ivel)) += tmp * normal_timefacfac(ivel);
#else
                C_umus_(row, sIndex(ic, ivel)) += tmp * proj_norm_timefacfac(ivel);
#endif
              }
            }
          }
        }

        return;
      }

      /*----------------------------------------------------------------------*
       *----------------------------------------------------------------------*/
      template <DRT::Element::DiscretizationType distype,
          DRT::Element::DiscretizationType slave_distype, unsigned int slave_numdof>
      void
      NitscheCoupling<distype, slave_distype, slave_numdof>::NIT_p_AdjointConsistency_SlaveTerms(
          const LINALG::Matrix<nsd_, 1>& normal_timefacfac,  ///< normal vector * timefacfac
          const double& velint_diff_normal_timefacfac,  ///< (velint_m - velint_s) * n * timefacfac
          const std::pair<bool, double>& s_row,         ///< scaling for slave row
          const std::pair<bool, double>& m_col,         ///< scaling for master col
          const std::pair<bool, double>& s_col,         ///< scaling for slave col
          bool only_rhs)
      {
        // 1) No-split no qunP option:
        /*                   \     /              i   \
      - |  { q }*n ,[ Du ]     | = |  { q }*n  ,[ u ]   |
        \                    /     \                 */

        // 2) qunP option:
        /*                       \     /              i      \
      - |  { q }*n ,[ Du ] P^n    | = |  { q }*n  ,[ u ] P^n  |
        \                        /     \                    */

#ifdef PROJECT_VEL_FOR_PRESSURE_ADJOINT
        LINALG::Matrix<nsd_, 1> proj_norm_timefacfac;
        proj_norm_timefacfac.Multiply(proj_normal_, normal_timefacfac);
#endif

        const double velint_diff_normal_timefacfac_ks =
            velint_diff_normal_timefacfac * s_row.second;
        for (unsigned ir = 0; ir < slave_nen_; ++ir)
        {
          // (qs*n,ks* um)
          rhC_us_(sPres(ir), 0) += funct_s_(ir) * velint_diff_normal_timefacfac_ks;
        }

        if (only_rhs) return;

        //-----------------------------------------------
        //    - (qs*n, ks *(Dum))
        //-----------------------------------------------
        const double facsm = s_row.second * m_col.second;
        for (unsigned ic = 0; ic < nen_; ++ic)
        {
          for (unsigned ir = 0; ir < slave_nen_; ++ir)
          {
            const unsigned row = sPres(ir);

            const double tmp = funct_s_m_dyad_(ir, ic) * facsm;

            for (unsigned ivel = 0; ivel < nsd_; ++ivel)
            {
              // -(qs*n, ks* Dum)
#ifndef PROJECT_VEL_FOR_PRESSURE_ADJOINT
              C_usum_(row, mIndex(ic, ivel)) -= tmp * normal_timefacfac(ivel);
              //        C_usum_(row,mIndex(ic,ivel)) += tmp*normal_timefacfac_ks(ivel); //TESTING!
#else
              C_usum_(row, mIndex(ic, ivel)) -= tmp * proj_norm_timefacfac(ivel);
#endif
            }
          }
        }

        //-----------------------------------------------
        //    + (qs*n, ks *(Dus))
        //-----------------------------------------------
        const double facss = s_row.second * s_col.second;
        for (unsigned ic = 0; ic < slave_nen_; ++ic)
        {
          for (unsigned ir = 0; ir < slave_nen_; ++ir)
          {
            const unsigned row = sPres(ir);

            const double tmp = funct_s_s_dyad_(ir, ic) * facss;

            for (unsigned ivel = 0; ivel < nsd_; ++ivel)
            {
              // +(qs*n, ks* Dus)
#ifndef PROJECT_VEL_FOR_PRESSURE_ADJOINT
              C_usus_(row, sIndex(ic, ivel)) += tmp * normal_timefacfac(ivel);
#else
              C_usus_(row, sIndex(ic, ivel)) += tmp * proj_norm_timefacfac(ivel);
#endif
            }
          }
        }
        return;
      }

      template <DRT::Element::DiscretizationType distype,
          DRT::Element::DiscretizationType slave_distype, unsigned int slave_numdof>
      void NitscheCoupling<distype, slave_distype, slave_numdof>::NIT_visc_Consistency_MasterTerms(
          const LINALG::Matrix<nsd_, nen_>& derxy_m,  ///< master deriv
          const LINALG::Matrix<nen_, 1>& funct_m,     ///< funct_m
          const std::pair<bool, double>& m_row,       ///< scaling for master row
          const std::pair<bool, double>& s_row,       ///< scaling for slave row
          const std::pair<bool, double>& m_col,       ///< scaling for master col
          bool only_rhs)
      {
        // viscous consistency term

        /*                           \       /                   i      \
      - |  [ v ],  { 2mu eps(u) }*n    | = + | [ v ],  { 2mu eps(u ) }*n  |
        \                            /       \                         */


        // here we use a non-optimal order to assemble the values into C_umum
        // however for this term we have to save operations
        const double facmm = m_row.second * m_col.second;
        const double facsm = s_row.second * m_col.second;

        for (unsigned ir = 0; ir < nen_; ++ir)
        {
          const double tmp_val = funct_m(ir) * facmm;

          for (unsigned ivel = 0; ivel < nsd_; ++ivel)
          {
            //-----------------------------------------------
            //    - (vm, (2*km*mum) *eps(Dum)*n)
            //-----------------------------------------------
            rhC_um_(mIndex(ir, ivel), 0) +=
                tmp_val * vderxy_m_normal_transposed_viscm_timefacfac_km_(ivel);
          }
        }


        if (s_row.first)
        {
          for (unsigned ir = 0; ir < slave_nen_; ++ir)
          {
            const double tmp_val = funct_s_(ir) * facsm;

            for (unsigned ivel = 0; ivel < nsd_; ++ivel)
            {
              //-----------------------------------------------
              //    + (vs, (2*km*mum) *eps(Dum)*n)
              //-----------------------------------------------

              // rhs
              rhC_us_(sIndex(ir, ivel), 0) -=
                  tmp_val * vderxy_m_normal_transposed_viscm_timefacfac_km_(ivel);
            }
          }
        }

        if (only_rhs) return;

        for (unsigned ic = 0; ic < nen_; ++ic)
        {
          const double normal_deriv_tmp = half_normal_deriv_m_viscm_timefacfac_km_(ic);

          for (unsigned ivel = 0; ivel < nsd_; ++ivel)
          {
            const double tmp_derxy_m = derxy_m(ivel, ic);
            for (unsigned jvel = 0; jvel < nsd_; ++jvel)
            {
              const unsigned col = mIndex(ic, jvel);

              double tmp = half_normal_viscm_timefacfac_km_(jvel) * tmp_derxy_m;
              if (ivel == jvel) tmp += normal_deriv_tmp;

              const double tmpm = tmp * facmm;
              for (unsigned ir = 0; ir < nen_; ++ir)
              {
                C_umum_(mIndex(ir, ivel), col) -= funct_m(ir) * tmpm;
              }

              if (s_row.first)
              {
                const double tmps = tmp * facsm;
                for (unsigned ir = 0; ir < slave_nen_; ++ir)
                {
                  C_usum_(sIndex(ir, ivel), col) += funct_s_(ir) * tmps;
                }
              }
            }
          }
        }

        return;
      }


      template <DRT::Element::DiscretizationType distype,
          DRT::Element::DiscretizationType slave_distype, unsigned int slave_numdof>
      void NitscheCoupling<distype, slave_distype, slave_numdof>::
          NIT_visc_Consistency_MasterTerms_Projected(
              const LINALG::Matrix<nsd_, nen_>& derxy_m,      ///< master deriv
              const LINALG::Matrix<nen_, 1>& funct_m,         ///< funct_m
              const LINALG::Matrix<nsd_, nsd_>& proj_matrix,  ///< projection matrix
              const double& km_viscm_fac,                     ///< scaling factor
              const std::pair<bool, double>& m_row,           ///< scaling for master row
              const std::pair<bool, double>& s_row,           ///< scaling for slave row
              const std::pair<bool, double>& m_col            ///< scaling for master col
          )
      {
        // 1) No-split WDBC option:
        //----------------------------------------------------------------------------------------
        /*                            \       /                   i      \
      - |  [ v ],  { 2mu eps(u) }*n    | = + | [ v ],  { 2mu eps(u ) }*n  |
        \                             /       \                         */
        //----------------------------------------------------------------------------------------

        // 2) (Normal - Tangential split):
        //----------------------------------------------------------------------------------------
        /*                                         \
      - |   { 2mu*eps(v) }*n  ,  [Du] P_n           |  =
        \                                         */

        /*                               i                   \
      + |  alpha* { 2mu*eps(v) }*n  , [ u ]  P_n              |
        \                                                   */

        //----------------------------------------------------------------------------------------

        // 2.0 * timefacfac * visceff_m * k_m * [\nabla N^(IX)]_k P^t_{kj}
        proj_matrix_derxy_m_.MultiplyTN(proj_matrix, derxy_m);  // Apply from right for consistency
        proj_matrix_derxy_m_.Scale(km_viscm_fac);

        // P_norm * {2.0 * timefacfac * visceff_m * 0.5 * (\nabla u + (\nabla u)^T)}
        vderxy_x_normal_transposed_viscx_timefacfac_kx_pmatrix_.MultiplyTN(
            proj_matrix, vderxy_m_normal_transposed_viscm_timefacfac_km_);

        // here we use a non-optimal order to assemble the values into C_umum
        // however for this term we have to save operations
        const double facmm = m_row.second * m_col.second;
        const double facsm = s_row.second * m_col.second;
        for (unsigned ic = 0; ic < nen_; ++ic)
        {
          // half_normal_deriv_m_viscm_timefacfac_km_ = 2.0 * timefacfac *
          // visceff_m*(0.5*normal(k)*derxy_m(k,ic))
          const double normal_deriv_tmp = half_normal_deriv_m_viscm_timefacfac_km_(ic);

          for (unsigned ivel = 0; ivel < nsd_; ++ivel)
          {
            const double tmp_derxy_m = derxy_m(ivel, ic);
            for (unsigned jvel = 0; jvel < nsd_; ++jvel)
            {
              const unsigned col = mIndex(ic, jvel);

              double tmp = half_normal_viscm_timefacfac_km_(jvel) * tmp_derxy_m;
              if (ivel == jvel) tmp += normal_deriv_tmp;

              for (unsigned ir = 0; ir < nen_; ++ir)
              {
                // C_umum_(mIndex(ir,ivel), col) -= funct_m(ir) * tmp;
                C_umum_(mIndex(ir, ivel), col) -=
                    funct_m(ir) * facmm *
                    (proj_matrix(jvel, ivel) * normal_deriv_tmp +
                        proj_matrix_derxy_m_(ivel, ic) * half_normal_(jvel));
              }


              if (!s_row.first) continue;


              for (unsigned ir = 0; ir < slave_nen_; ++ir)
              {
                C_usum_(sIndex(ir, ivel), col) +=
                    funct_s_(ir) * facsm *
                    (proj_matrix(jvel, ivel) * normal_deriv_tmp +
                        proj_matrix_derxy_m_(ivel, ic) * half_normal_(jvel));
              }
            }
          }
        }


        for (unsigned ir = 0; ir < nen_; ++ir)
        {
          const double tmp_val = facmm * funct_m(ir);

          for (unsigned ivel = 0; ivel < nsd_; ++ivel)
          {
            //-----------------------------------------------
            //    - (vm, (2*km*mum) *eps(Dum)*n)
            //-----------------------------------------------
            rhC_um_(mIndex(ir, ivel), 0) +=
                tmp_val * vderxy_x_normal_transposed_viscx_timefacfac_kx_pmatrix_(ivel);
          }
        }


        if (!s_row.first) return;


        for (unsigned ir = 0; ir < slave_nen_; ++ir)
        {
          const double tmp_val = facsm * funct_s_(ir);

          for (unsigned ivel = 0; ivel < nsd_; ++ivel)
          {
            //-----------------------------------------------
            //    + (vs, (2*km*mum) *eps(Dum)*n)
            //-----------------------------------------------

            // rhs
            rhC_us_(sIndex(ir, ivel), 0) -=
                tmp_val * vderxy_x_normal_transposed_viscx_timefacfac_kx_pmatrix_(ivel);
          }
        }
        return;
      }


      /*----------------------------------------------------------------------*
       *----------------------------------------------------------------------*/
      template <DRT::Element::DiscretizationType distype,
          DRT::Element::DiscretizationType slave_distype, unsigned int slave_numdof>
      void NitscheCoupling<distype, slave_distype, slave_numdof>::NIT_visc_Consistency_SlaveTerms(
          const LINALG::Matrix<nsd_, slave_nen_>& derxy_s,  ///< slave shape function derivatives
          const LINALG::Matrix<nen_, 1>& funct_m,           ///< funct_m
          const std::pair<bool, double>& m_row,             ///< scaling for master row
          const std::pair<bool, double>& s_row,             ///< scaling for slave row
          const std::pair<bool, double>& s_col,             ///< scaling for slave col
          bool only_rhs)
      {
        // diagonal block (i,i): +/-2*ks*mus * ...
        /*       nsd_
         *       *---*
         *       \    dN                     dN
         *  N  *  *   --  * 0.5 * n_j + N  * --  * n_i * 0.5
         *       /    dxj                    dxi
         *       *---*
         *       j = 1
         */

        // off-diagonal block (i,j) : +/-2*ks*mus * ...
        /*       dN
         *  N  * -- * n_j * 0.5
         *       dxi
         */

        const double facms = m_row.second * s_col.second;
        const double facss = s_row.second * s_col.second;

        for (unsigned ir = 0; ir < nen_; ++ir)
        {
          const double tmp_val = funct_m(ir) * facms;
          for (unsigned ivel = 0; ivel < nsd_; ++ivel)
          {
            // rhs
            rhC_um_(mIndex(ir, ivel), 0) +=
                tmp_val * vderxy_s_normal_transposed_viscs_timefacfac_ks_(ivel);
          }
        }

        for (unsigned ir = 0; ir < slave_nen_; ++ir)
        {
          const double tmp_val = funct_s_(ir) * facss;
          for (unsigned ivel = 0; ivel < nsd_; ++ivel)
          {
            // rhs
            rhC_us_(sIndex(ir, ivel), 0) -=
                tmp_val * vderxy_s_normal_transposed_viscs_timefacfac_ks_(ivel);
          }
        }

        if (only_rhs) return;

        for (unsigned ic = 0; ic < slave_nen_; ++ic)
        {
          const double normal_deriv_tmp = half_normal_deriv_s_viscs_timefacfac_ks_(ic);

          for (unsigned ivel = 0; ivel < nsd_; ++ivel)
          {
            const double tmp_derxy_s = derxy_s(ivel, ic);

            for (unsigned jvel = 0; jvel < nsd_; ++jvel)
            {
              const unsigned col = sIndex(ic, jvel);

              double tmp = half_normal_viscs_timefacfac_ks_(jvel) * tmp_derxy_s;

              if (ivel == jvel) tmp += normal_deriv_tmp;

              const double tmpm = tmp * facms;
              for (unsigned ir = 0; ir < nen_; ++ir)
              {
                //-----------------------------------------------
                //    - (vm, (2*ks*mus) *eps(Dus)*n)
                //-----------------------------------------------
                C_umus_(mIndex(ir, ivel), col) -= funct_m(ir) * tmpm;
              }

              const double tmps = tmp * facss;
              for (unsigned ir = 0; ir < slave_nen_; ++ir)
              {
                //-----------------------------------------------
                //    + (vs, (2*ks*mus) *eps(Dus)*n)
                //-----------------------------------------------
                // diagonal block
                C_usus_(sIndex(ir, ivel), col) += funct_s_(ir) * tmps;
              }
            }
          }
        }
        return;
      }

      template <DRT::Element::DiscretizationType distype,
          DRT::Element::DiscretizationType slave_distype, unsigned int slave_numdof>
      void NitscheCoupling<distype, slave_distype, slave_numdof>::
          NIT_visc_AdjointConsistency_MasterTerms(
              const LINALG::Matrix<nen_, 1>& funct_m,  ///< funct * timefacfac
              const LINALG::Matrix<nsd_, nen_>&
                  derxy_m,  ///< spatial derivatives of coupling master shape functions
              const LINALG::Matrix<nsd_, 1>& normal,  ///< normal-vector
              const double& viscm_fac,                ///< scaling factor
              const std::pair<bool, double>& m_row,   ///< scaling for master row
              const std::pair<bool, double>& m_col,   ///< scaling for master col
              const std::pair<bool, double>& s_col,   ///< scaling for slave col
              bool only_rhs)
      {
        /*                                \       /                             i   \
      - |  alpha* { 2mu*eps(v) }*n , [ Du ] |  =  |  alpha* { 2mu eps(v) }*n ,[ u ]   |
        \                                 /       \                                */
        // (see Burman, Fernandez 2009)
        // +1.0 symmetric
        // -1.0 antisymmetric

        //-----------------------------------------------------------------
        // viscous adjoint consistency term

        double tmp_fac = adj_visc_scale_ * viscm_fac;
        derxy_m_viscm_timefacfac_.Update(tmp_fac, derxy_m);  // 2 * mu_m * timefacfac *
                                                             // derxy_m(k,ic)

        static LINALG::Matrix<nsd_, nsd_> velint_diff_dyad_normal, velint_diff_dyad_normal_symm;
        velint_diff_dyad_normal.MultiplyNT(velint_diff_, normal);

        for (unsigned jvel = 0; jvel < nsd_; ++jvel)
        {
          for (unsigned ivel = 0; ivel < nsd_; ++ivel)
          {
            velint_diff_dyad_normal_symm(ivel, jvel) =
                velint_diff_dyad_normal(ivel, jvel) + velint_diff_dyad_normal(jvel, ivel);
          }
        }

        const double facm = m_row.second * 0.5;
        for (unsigned ir = 0; ir < nen_; ++ir)
        {
          for (unsigned jvel = 0; jvel < nsd_; ++jvel)
          {
            const double derxy_m_viscm_timefacfac_km_half_tmp =
                derxy_m_viscm_timefacfac_(jvel, ir) * facm;

            for (unsigned ivel = 0; ivel < nsd_; ++ivel)
            {
              // rhs
              rhC_um_(mIndex(ir, ivel), 0) +=
                  derxy_m_viscm_timefacfac_km_half_tmp * velint_diff_dyad_normal_symm(ivel, jvel);
            }
          }
        }

        if (only_rhs) return;

        const double facmm = m_row.second * m_col.second;
        const double facms = m_row.second * s_col.second;

        normal_deriv_m_viscm_km_.MultiplyTN(
            derxy_m_viscm_timefacfac_, half_normal_);  // half_normal(k)*derxy_m(k,ic)*viscm*km

        // here we use a non-optimal order to assemble the values into C_umum
        // however for this term we have to save operations
        for (unsigned ir = 0; ir < nen_; ++ir)
        {
          const double normal_deriv_tmp = normal_deriv_m_viscm_km_(ir);

          for (unsigned jvel = 0; jvel < nsd_; ++jvel)
          {
            const double tmp_derxy_m = derxy_m_viscm_timefacfac_(jvel, ir);
            for (unsigned ivel = 0; ivel < nsd_; ++ivel)
            {
              const unsigned row = mIndex(ir, ivel);

              double tmp = half_normal_(ivel) * tmp_derxy_m;
              if (ivel == jvel) tmp += normal_deriv_tmp;

              const double tmpm = tmp * facmm;
              for (unsigned ic = 0; ic < nen_; ++ic)
              {
                C_umum_(row, mIndex(ic, jvel)) -= funct_m(ic) * tmpm;
              }


              if (s_col.first)
              {
                const double tmps = tmp * facms;
                for (unsigned ic = 0; ic < slave_nen_; ++ic)
                {
                  C_umus_(row, sIndex(ic, jvel)) += funct_s_(ic) * tmps;
                }
              }
            }
          }
        }
        return;
      }

      /*----------------------------------------------------------------------*
       *----------------------------------------------------------------------*/
      template <DRT::Element::DiscretizationType distype,
          DRT::Element::DiscretizationType slave_distype, unsigned int slave_numdof>
      void NitscheCoupling<distype, slave_distype, slave_numdof>::
          NIT_visc_AdjointConsistency_MasterTerms_Projected(
              const LINALG::Matrix<nsd_, nen_>&
                  derxy_m_viscm_timefacfac_km,  ///< master shape function derivatives * timefacfac
                                                ///< * 2 * mu_m * kappa_m
              const LINALG::Matrix<nen_, 1>& funct_m,  ///< embedded element funct *mu*timefacfac
              const LINALG::Matrix<nsd_, 1>& normal,   ///< normal vector
              const std::pair<bool, double>& m_row,    ///< scaling for master row
              const std::pair<bool, double>& m_col,    ///< scaling for master col
              const std::pair<bool, double>& s_col     ///< scaling for slave col
          )
      {
        // 1) No-split WDBC option:
        //----------------------------------------------------------------------------------------
        /*                                 \       /                             i   \
      - |  alpha* { 2mu*eps(v) }*n , [ Du ] |  =  |  alpha* { 2mu eps(v) }*n ,[ u ]   |
        \                                  /       \                                */
        //----------------------------------------------------------------------------------------

        // 2) (Normal - Tangential split):
        //----------------------------------------------------------------------------------------
        /*                                                       \
      - |  alpha* { 2mu*eps(v) }*n  ,  [Du] *  stab_fac *  P      |  =
        \                                                       */

        /*                               i                       \
      + |  alpha* { 2mu*eps(v) }*n  , [ u ]  *  stab_fac *  P     |
        \                                                       */

        //----------------------------------------------------------------------------------------

        // (see Burman, Fernandez 2009)
        // alpha =  +1.0 symmetric
        //          -1.0 antisymmetric

        // timefacfac                   = theta*dt*fac
        // derxy_m_viscm_timefacfac_km  = alpha * 2 * mu_m * timefacfac * derxy_m(k,IX)

        // normal_deriv_m_viscm_km_    = alpha * half_normal(k)* 2 * mu_m * timefacfac *
        // derxy_m(k,IX)
        //                             = alpha * mu_m * timefacfac * c(IX)

        // proj_matrix_derxy_m_        = alpha * 2 * mu_m * timefacfac * derxy_m_(k,ir) * P_{jk}
        //                             = alpha * 2 * mu_m * timefacfac * p^t_1(ir,j)

        //  // proj_tang_derxy_m = 2.0 * P^t_{jk} * derxy_m(k,IX) * mu_m * timefacfac * km
        //  //                   = 2.0 * mu_m * timefacefac * km * p_1(IX,j)

        // here we use a non-optimal order to assemble the values into C_umum
        // however for this term we have to save operations
        const double facmm = m_row.second * m_col.second;
        const double facms = m_row.second * s_col.second;
        for (unsigned ir = 0; ir < nen_; ++ir)
        {
          // alpha * mu_m * timefacfac * \sum_k dN^(ir)/dx_k * n_k
          const double normal_deriv_tmp = normal_deriv_m_viscm_km_(ir);

          for (unsigned jvel = 0; jvel < nsd_; ++jvel)
          {
            for (unsigned ivel = 0; ivel < nsd_; ++ivel)
            {
              const unsigned row = mIndex(ir, ivel);

              const double tmpm = facmm * (proj_matrix_(jvel, ivel) * normal_deriv_tmp +
                                              proj_matrix_derxy_m_(jvel, ir) * half_normal_(ivel));
              for (unsigned ic = 0; ic < nen_; ++ic)
              {
                C_umum_(row, mIndex(ic, jvel)) -= funct_m(ic) * tmpm;
              }

              if (s_col.first)
              {
                const double tmps =
                    facms * (proj_matrix_(jvel, ivel) * normal_deriv_tmp +
                                proj_matrix_derxy_m_(jvel, ir) * half_normal_(ivel));
                for (unsigned ic = 0; ic < slave_nen_; ++ic)
                {
                  C_umus_(row, sIndex(ic, jvel)) += funct_s_(ic) * tmps;
                }
              }
            }
          }
        }

        // Can this be made more effective?
        // static LINALG::Matrix<nsd_,nsd_> velint_proj_norm_diff_dyad_normal,
        // velint_proj_norm_diff_dyad_normal_symm;
        // velint_diff_proj_normal_ = (u^m_k - u^s_k) P^n_{kj} * n
        velint_proj_norm_diff_dyad_normal_.MultiplyNT(velint_diff_proj_matrix_, normal);

        for (unsigned jvel = 0; jvel < nsd_; ++jvel)
        {
          for (unsigned ivel = 0; ivel < nsd_; ++ivel)
          {
            velint_proj_norm_diff_dyad_normal_symm_(ivel, jvel) =
                velint_proj_norm_diff_dyad_normal_(ivel, jvel) +
                velint_proj_norm_diff_dyad_normal_(jvel, ivel);
          }
        }

        const double facm = m_row.second * 0.5;
        for (unsigned ir = 0; ir < nen_; ++ir)
        {
          for (unsigned jvel = 0; jvel < nsd_; ++jvel)
          {
            const double derxy_m_viscm_timefacfac_km_half_tmp =
                derxy_m_viscm_timefacfac_km(jvel, ir) * facm;

            for (unsigned ivel = 0; ivel < nsd_; ++ivel)
            {
              // rhs
              rhC_um_(mIndex(ir, ivel), 0) += derxy_m_viscm_timefacfac_km_half_tmp *
                                              velint_proj_norm_diff_dyad_normal_symm_(ivel, jvel);
            }
          }
        }

        return;
      }



      /*----------------------------------------------------------------------*
       *----------------------------------------------------------------------*/
      template <DRT::Element::DiscretizationType distype,
          DRT::Element::DiscretizationType slave_distype, unsigned int slave_numdof>
      void
      NitscheCoupling<distype, slave_distype, slave_numdof>::NIT_visc_AdjointConsistency_SlaveTerms(
          const LINALG::Matrix<nen_, 1>& funct_m,  ///< embedded element funct *mu*timefacfac
          const LINALG::Matrix<nsd_, slave_nen_>&
              derxy_s_viscs_timefacfac_ks,  ///< master shape function derivatives * timefacfac * 2
                                            ///< * mu_m * kappa_m
          const LINALG::Matrix<nsd_, 1>& normal,
          const std::pair<bool, double>& s_row,  ///< scaling for slave row
          const std::pair<bool, double>& m_col,  ///< scaling for master col
          const std::pair<bool, double>& s_col,  ///< scaling for slave col
          bool only_rhs)
      {
        /*                                \       /                             i   \
      - |  alpha* { 2mu*eps(v) }*n , [ Du ] |  =  |  alpha* { 2mu eps(v) }*n ,[ u ]   |
        \                                 /       \                                */
        // (see Burman, Fernandez 2009)
        // +1.0 symmetric
        // -1.0 antisymmetric

        // diagonal block (i,i): +/-2*km*mum * alpha * ...
        /*
         *       nsd_
         *       *---*
         *       \    dN                    dN
         *        *   --  * 0.5 * n_j *N +  --  * n_i * 0.5 *N
         *       /    dxj                   dxi
         *       *---*
         *       j = 1
         */

        // off-diagonal block (i,j) : +/-2*km*mum * alpha * ...
        /*   dN
         *   -- * n_i * 0.5 * N
         *   dxj
         */

        static LINALG::Matrix<nsd_, nsd_> velint_diff_dyad_normal, velint_diff_dyad_normal_symm;
        velint_diff_dyad_normal.MultiplyNT(velint_diff_, normal);

        for (unsigned jvel = 0; jvel < nsd_; ++jvel)
        {
          for (unsigned ivel = 0; ivel < nsd_; ++ivel)
          {
            velint_diff_dyad_normal_symm(ivel, jvel) =
                velint_diff_dyad_normal(ivel, jvel) + velint_diff_dyad_normal(jvel, ivel);
          }
        }

        const double facs = s_row.second * 0.5;
        for (unsigned ir = 0; ir < nen_; ++ir)
        {
          for (unsigned jvel = 0; jvel < nsd_; ++jvel)
          {
            const double derxy_s_viscs_timefacfac_ks_half_tmp =
                derxy_s_viscs_timefacfac_ks(jvel, ir) * facs;

            for (unsigned ivel = 0; ivel < nsd_; ++ivel)
            {
              // rhs
              rhC_us_(sIndex(ir, ivel), 0) +=
                  derxy_s_viscs_timefacfac_ks_half_tmp * velint_diff_dyad_normal_symm(ivel, jvel);
            }
          }
        }

        if (only_rhs) return;

        const double facsm = s_row.second * m_col.second;
        const double facss = s_row.second * s_col.second;
        normal_deriv_s_viscs_ks_.MultiplyTN(
            derxy_s_viscs_timefacfac_ks, half_normal_);  // half_normal(k)*derxy_m(k,ic)*viscm*km

        for (unsigned ir = 0; ir < slave_nen_; ++ir)
        {
          const double normal_deriv_tmp = normal_deriv_s_viscs_ks_(ir);

          for (unsigned jvel = 0; jvel < nsd_; ++jvel)
          {
            const double tmp_derxy_s = derxy_s_viscs_timefacfac_ks(jvel, ir);
            for (unsigned ivel = 0; ivel < nsd_; ++ivel)
            {
              const unsigned row = sIndex(ir, ivel);

              double tmp = half_normal_(ivel) * tmp_derxy_s;
              if (ivel == jvel) tmp += normal_deriv_tmp;

              const double tmpm = tmp * facsm;
              for (unsigned ic = 0; ic < nen_; ++ic)
              {
                C_usum_(row, mIndex(ic, jvel)) -= funct_m(ic) * tmpm;
              }

              if (s_col.first)
              {
                const double tmps = tmp * facss;
                for (unsigned ic = 0; ic < slave_nen_; ++ic)
                {
                  C_usus_(row, sIndex(ic, jvel)) += funct_s_(ic) * tmps;
                }
              }
            }
          }
        }
      }


      /*----------------------------------------------------------------------*
       *----------------------------------------------------------------------*/
      template <DRT::Element::DiscretizationType distype,
          DRT::Element::DiscretizationType slave_distype, unsigned int slave_numdof>
      void NitscheCoupling<distype, slave_distype, slave_numdof>::
          NIT_visc_Neumann_AdjointConsistency_MasterTerms_Projected(
              const LINALG::Matrix<nsd_, nen_>&
                  derxy_m_viscm_timefacfac_km,  ///< master shape function derivatives * timefacfac
                                                ///< * 2 * mu_m * kappa_m
              const LINALG::Matrix<nsd_, nen_>& derxy_m,  ///< master deriv
              const LINALG::Matrix<nsd_, nsd_>&
                  vderxy_m,  ///< coupling master spatial velocity derivatives
              const LINALG::Matrix<nen_, 1>& funct_m,  ///< embedded element funct *mu*timefacfac
              const LINALG::Matrix<nsd_, 1>& normal,   ///< normal vector
              const std::pair<bool, double>& m_row,    ///< scaling for master row
              const std::pair<bool, double>& mstr_col  ///< scaling for master col
          )
      {
        // 1) No-split WDBC option:
        /*           \
        |  v  ,  0   |
        \           */

        // 2) (Normal - Tangential split):
        /*                                                                                                      \
      - |  alpha* (\epsilon * \gamma * h_E)/(epsilon + \gamma *h_E)  { 2*eps(v) }*n ,{ 2mu eps(Du)
      }*n P_t       | =
        \ */

        /*                                                                               i \
      + |  alpha* (\epsilon * \gamma * h_E)/(epsilon + \gamma *h_E)  { 2*eps(v) }*n ,{ 2mu eps(u)
      }*n P_t       |
        \ */

        /*                                                                                \
      - |  alpha* { 2mu*eps(v) }*n  ,   g  ( epsilon*gamma*h_E/(gamma*h_E+epsilon) * P_t)  |
        \                                                                                */

        // (see Burman, Fernandez 2009)
        // +1.0 symmetric
        // -1.0 antisymmetric

        const double facmm = m_row.second * mstr_col.second;

        //  //normal_deriv_m_viscm_km_ = 2.0 * alpha * half_normal(k) * mu_m * timefacfac * km *
        //  derxy_m(k,ic)
        //  //                         = mu_m * alpha * timefacfac *km * c(ix)
        //  normal_deriv_m_viscm_km_.MultiplyTN(derxy_m_viscm_timefacfac_km, half_normal_);

        //  LINALG::Matrix<nen_,1> normal_deriv_m;
        normal_deriv_m_.MultiplyTN(derxy_m, half_normal_);
        normal_deriv_m_.Scale(2.0);  // 2.0 * half_normal(k) * derxy_m(k,ix) =c(ix)

        //  LINALG::Matrix<nsd_,nen_> proj_tang_derxy_m_1;
        //  // proj_matrix_derxy_m_ = alpha * 2.0 * P^t_{jk} * derxy_m(k,IX) * mu_m * timefacfac *
        //  km
        //  //                      = 2.0 * mu_m * timefacefac * km * p_1(IX,j)
        //  // IF P^t_{jk} IS SYMMETRIC: p_1(IX,j) = p_2(IX,j)
        //  proj_tang_derxy_m_1.Multiply(proj_matrix_,derxy_m_viscm_timefacfac_km);

        //  LINALG::Matrix<nsd_,nen_> proj_tang_derxy_m_2;
        //  // proj_tang_derxy_m = 2.0 * P^t_{jk} * derxy_m(k,IX) * mu_m * timefacfac * km
        //  //                   = 2.0 * mu_m * timefacefac * km * p_2(IX,j)
        //  // IF P^t_{jk} IS SYMMETRIC: p_1(IX,j) = p_2(IX,j)
        //  proj_tang_derxy_m_2.MultiplyTN(proj_matrix_,derxy_m_viscm_timefacfac_km);

        //  LINALG::Matrix<nen_,nen_> derxy_m_P_derxy_m;
        // derxy_m_P_derxy_m = 2.0 * derxy_m(j,IC) P^t_{jk} * derxy_m(k,IR) * mu_m * timefacfac * km
        //                   = 2.0 * C(IC,IR) * mu_m * timefacfac * km
        //  derxy_m_P_derxy_m.MultiplyTN(derxy_m,proj_tang_derxy_m_1);
        derxy_m_P_derxy_m_.MultiplyTN(derxy_m, proj_matrix_derxy_m_);

        // here we use a non-optimal order to assemble the values into C_umum
        // however for this term we have to save operations
        for (unsigned ir = 0; ir < nen_; ++ir)
        {
          const double normal_deriv_tmp = normal_deriv_m_viscm_km_(
              ir);  // alpha * mu_m * timefacfac *km * \sum_k dN^(ir)/dx_k * n_k

          for (unsigned jvel = 0; jvel < nsd_; ++jvel)
          {
            for (unsigned ivel = 0; ivel < nsd_; ++ivel)
            {
              const unsigned row = mIndex(ir, ivel);

              for (unsigned ic = 0; ic < nen_; ++ic)
              {
                C_umum_(row, mIndex(ic, jvel)) -=
                    facmm *
                    (normal_deriv_m_(ic) *
                            (proj_matrix_(jvel, ivel) * normal_deriv_tmp +
                                proj_matrix_derxy_m_(jvel, ir) * half_normal_(ivel)) +
                        normal(ivel) * half_normal_(jvel) * derxy_m_P_derxy_m_(ic, ir) +
                        normal_deriv_m_(ir) * proj_matrix_derxy_m_(ivel, ic) * half_normal_(jvel));
              }
            }
          }
        }

        // 2.0 * timefacfac * visceff_m * 0.5* (\nabla u + (\nabla u)^T) * normal
        // vderxy_m_normal_transposed_viscm_timefacfac_km_

        vderxy_m_normal_tang_.Multiply(vderxy_m, normal);
        vderxy_m_normal_transposed_.MultiplyTN(vderxy_m, normal);
        vderxy_m_normal_transposed_.Update(
            1.0, vderxy_m_normal_tang_, 1.0);  // (\nabla u + (\nabla u)^T) * normal

        // (\nabla u + (\nabla u)^T) * normal * P^t
        vderxy_m_normal_tang_.MultiplyTN(proj_matrix_, vderxy_m_normal_transposed_);

        // 2.0 * derxy_m(k,IX) * mu_m * timefacfac * km ( (\nabla u + (\nabla u)^T) * normal * P^t
        // )_k
        static LINALG::Matrix<nen_, 1> tmp_rhs;
        tmp_rhs.MultiplyTN(derxy_m_viscm_timefacfac_km, vderxy_m_normal_tang_);

        for (unsigned ir = 0; ir < nen_; ++ir)
        {
          // alpha * mu_m * timefacfac *km * \sum_k dN^(ir)/dx_k * n_k
          const double normal_deriv_tmp = normal_deriv_m_viscm_km_(ir);

          for (unsigned ivel = 0; ivel < nsd_; ++ivel)
          {
            rhC_um_(mIndex(ir, ivel), 0) +=
                facmm *
                (normal_deriv_tmp * vderxy_m_normal_tang_(ivel) + tmp_rhs(ir) * half_normal_(ivel));
          }
        }
      }

      /*----------------------------------------------------------------------*
       *----------------------------------------------------------------------*/
      template <DRT::Element::DiscretizationType distype,
          DRT::Element::DiscretizationType slave_distype, unsigned int slave_numdof>
      void NitscheCoupling<distype, slave_distype, slave_numdof>::NIT_Stab_Penalty(
          const LINALG::Matrix<nen_, 1>& funct_m,  ///< funct
          const double& timefacfac,                ///< time integration factor
          const std::pair<bool, double>& m_row,    ///< scaling for master row
          const std::pair<bool, double>& s_row,    ///< scaling for slave row
          const std::pair<bool, double>& m_col,    ///< scaling for master col
          const std::pair<bool, double>& s_col,    ///< scaling for slave col
          bool only_rhs)
      {
        TEUCHOS_FUNC_TIME_MONITOR("FLD::NIT_Stab_Penalty");

        // viscous stability term

        // combined viscous and inflow stabilization for one-sided problems (XFSI)
        // gamma_combined = max(alpha*mu/hk, |u*n| )
        /*                      _        \        /                     i   _   \
       |  gamma_combined *  v , u - u     | =  - |   gamma/h_K *  v , (u  - u)   |
        \                                /        \                            */

        // just viscous stabilization for two-sided problems (XFF, XFFSI)
        /*                                  \        /                           i   \
       |  gamma*mu/h_K *  [ v ] , [ Du ]     | =  - |   gamma*mu/h_K * [ v ], [ u ]   |
        \                                   /        \                              */

        // + gamma*mu/h_K (vm, um))

        const double stabfac_timefacfac_m = timefacfac * m_row.second;
        velint_diff_timefacfac_stabfac_.Update(stabfac_timefacfac_m, velint_diff_, 0.0);

        for (unsigned ir = 0; ir < nen_; ++ir)
        {
          const double tmp_val = funct_m(ir);

          for (unsigned ivel = 0; ivel < nsd_; ++ivel)
          {
            // +(stab * vm, u_DBC) (weak dirichlet case) or from
            // +(stab * vm, u_s)
            rhC_um_(mIndex(ir, ivel), 0) -= tmp_val * velint_diff_timefacfac_stabfac_(ivel);
          }
        }

        if (s_row.first)
        {
          const double stabfac_timefacfac_s = timefacfac * s_row.second;
          velint_diff_timefacfac_stabfac_.Update(stabfac_timefacfac_s, velint_diff_, 0.0);

          for (unsigned ir = 0; ir < slave_nen_; ++ir)
          {
            const double tmp_val = funct_s_(ir);

            for (unsigned ivel = 0; ivel < nsd_; ++ivel)
            {
              // +(stab * vs, um)
              // -(stab * vs, us)
              rhC_us_(sIndex(ir, ivel), 0) += tmp_val * velint_diff_timefacfac_stabfac_(ivel);
            }
          }
        }

        if (only_rhs) return;

        const double stabfac_timefacfac_mm = timefacfac * m_row.second * m_col.second;

        for (unsigned ic = 0; ic < nen_; ++ic)
        {
          for (unsigned ir = 0; ir < nen_; ++ir)
          {
            const double tmp_val = funct_m_m_dyad_(ir, ic) * stabfac_timefacfac_mm;

            for (unsigned ivel = 0; ivel < nsd_; ++ivel)
            {
              C_umum_(mIndex(ir, ivel), mIndex(ic, ivel)) += tmp_val;
            }
          }
        }

        if (s_col.first)
        {
          // - gamma*mu/h_K (vm, us))
          // - gamma*mu/h_K (vs, um))

          const double stabfac_timefacfac_ms = timefacfac * m_row.second * s_col.second;

          for (unsigned ic = 0; ic < slave_nen_; ++ic)
          {
            for (unsigned ir = 0; ir < nen_; ++ir)
            {
              const double tmp_val = funct_s_m_dyad_(ic, ir) * stabfac_timefacfac_ms;

              for (unsigned ivel = 0; ivel < nsd_; ++ivel)
              {
                C_umus_(mIndex(ir, ivel), sIndex(ic, ivel)) -= tmp_val;
              }
            }
          }
        }

        if (s_row.first && s_col.first)
        {
          const double stabfac_timefacfac_ss = timefacfac * s_row.second * s_col.second;

          for (unsigned ic = 0; ic < slave_nen_; ++ic)
          {
            // + gamma*mu/h_K (vs, us))
            for (unsigned ir = 0; ir < slave_nen_; ++ir)
            {
              const double tmp_val = funct_s_s_dyad_(ir, ic) * stabfac_timefacfac_ss;

              for (unsigned ivel = 0; ivel < nsd_; ++ivel)
              {
                C_usus_(sIndex(ir, ivel), sIndex(ic, ivel)) += tmp_val;
              }
            }
          }
        }

        if (s_row.first)
        {
          const double stabfac_timefacfac_sm = timefacfac * s_row.second * m_col.second;

          for (unsigned ic = 0; ic < nen_; ++ic)
          {
            for (unsigned ir = 0; ir < slave_nen_; ++ir)
            {
              const double tmp_val = funct_s_m_dyad_(ir, ic) * stabfac_timefacfac_sm;

              for (unsigned ivel = 0; ivel < nsd_; ++ivel)
              {
                C_usum_(sIndex(ir, ivel), mIndex(ic, ivel)) -= tmp_val;
              }
            }
          }
        }

        return;
      }

      /*----------------------------------------------------------------------*
       *----------------------------------------------------------------------*/
      template <DRT::Element::DiscretizationType distype,
          DRT::Element::DiscretizationType slave_distype, unsigned int slave_numdof>
      void NitscheCoupling<distype, slave_distype, slave_numdof>::NIT_Stab_Penalty_lin(
          const LINALG::Matrix<nen_, 1>& funct_m,  ///< funct
          const double& timefacfac,                ///< time integration factor
          const std::pair<bool, double>& m_row,    ///< scaling for master row
          const std::pair<bool, double>&
              m_row_linm1,  ///< linearization of scaling for master row w.r.t. master comp. one
          const std::pair<bool, double>&
              m_row_linm2,  ///< linearization of scaling for master row w.r.t. master comp. two
          const std::pair<bool, double>&
              m_row_linm3,  ///< linearization of scaling for master row w.r.t. master comp. three
          bool only_rhs)
      {
        TEUCHOS_FUNC_TIME_MONITOR("FLD::NIT_Stab_Penalty_linearization");

        if (only_rhs) return;

        const double stabfac_timefacfac_m = timefacfac;
        velint_diff_timefacfac_stabfac_.Update(stabfac_timefacfac_m, velint_diff_, 0.0);

        for (unsigned ir = 0; ir < nen_; ++ir)
        {
          const double tmp_val = funct_m(ir);

          for (unsigned ivel = 0; ivel < nsd_; ++ivel)
          {
            for (unsigned ic = 0; ic < nen_; ++ic)
            {
              C_umum_(mIndex(ir, ivel), mIndex(ic, 0)) += tmp_val *
                                                          velint_diff_timefacfac_stabfac_(ivel) *
                                                          funct_m(ic) * m_row_linm1.second;
              C_umum_(mIndex(ir, ivel), mIndex(ic, 1)) += tmp_val *
                                                          velint_diff_timefacfac_stabfac_(ivel) *
                                                          funct_m(ic) * m_row_linm2.second;
              C_umum_(mIndex(ir, ivel), mIndex(ic, 2)) += tmp_val *
                                                          velint_diff_timefacfac_stabfac_(ivel) *
                                                          funct_m(ic) * m_row_linm3.second;
            }
          }
        }
        return;
      }


      /*----------------------------------------------------------------------*
       *----------------------------------------------------------------------*/
      template <DRT::Element::DiscretizationType distype,
          DRT::Element::DiscretizationType slave_distype, unsigned int slave_numdof>
      void NitscheCoupling<distype, slave_distype, slave_numdof>::NIT_Stab_Penalty_Projected(
          const LINALG::Matrix<nen_, 1>& funct_m,               ///< funct
          const LINALG::Matrix<nsd_, nsd_>& projection_matrix,  ///< projection_matrix
          const LINALG::Matrix<nsd_, 1>&
              velint_diff_proj_matrix,           ///< velocity difference projected
          const double& timefacfac,              ///< time integration factor
          const std::pair<bool, double>& m_row,  ///< scaling for master row
          const std::pair<bool, double>& s_row,  ///< scaling for slave row
          const std::pair<bool, double>& m_col,  ///< scaling for master col
          const std::pair<bool, double>& s_col   ///< scaling for slave col
      )
      {
        TEUCHOS_FUNC_TIME_MONITOR("FLD::NIT_Stab_Penalty");

        // viscous stability term

        // combined viscous and inflow stabilization for one-sided problems (XFSI)
        // gamma_combined = max(alpha*mu/hk, |u*n| )
        /*                      _        \        /                     i   _   \
       |  gamma_combined *  v , u - u     | =  - |   gamma/h_K *  v , (u  - u)   |
        \                                /        \                            */

        // just viscous stabilization for two-sided problems (XFF, XFFSI)
        /*                                  \        /                           i   \
       |  gamma*mu/h_K *  [ v ] , [ Du ]     | =  - |   gamma*mu/h_K * [ v ], [ u ]   |
        \                                   /        \                              */

        // + gamma*mu/h_K (vm, um))


        // 2) (Normal - Tangential split):
        /*                                                                                       \
      + |   [v]  ,  [ Du ] ( gamma_comb_n P_n + {mu} * gamma*h_E/(gamma*h_E+epsilon) * P_t )     | =
        \                                                                                       */

        /*             i                                                                         \
      - |   [v]  ,  [ u ] ( gamma_comb_n P_n + {mu} * gamma*h_E/(gamma*h_E+epsilon) * P_t )      | =
        \                                                                                       */

        velint_diff_proj_matrix_ = velint_diff_proj_matrix;
        proj_matrix_ = projection_matrix;

        const double stabfac_timefacfac_mm = timefacfac * m_row.second * m_col.second;

        for (unsigned ic = 0; ic < nen_; ++ic)
        {
          for (unsigned ir = 0; ir < nen_; ++ir)
          {
            const double stab_funct_m_m_dyad_iric = funct_m_m_dyad_(ir, ic) * stabfac_timefacfac_mm;
            for (unsigned ivel = 0; ivel < nsd_; ++ivel)
            {
              const unsigned col = mIndex(ic, ivel);
              //        C_umum_(mIndex(ir,ivel), mIndex(ic,ivel)) += tmp_val;
              for (unsigned jvel = 0; jvel < nsd_; ++jvel)
              {
                C_umum_(mIndex(ir, jvel), col) +=
                    stab_funct_m_m_dyad_iric * proj_matrix_(ivel, jvel);
              }
            }
          }
        }

        const double stabfac_timefacfac_m = timefacfac * m_row.second;
        velint_diff_timefacfac_stabfac_.Update(stabfac_timefacfac_m, velint_diff_proj_matrix_, 0.0);
        for (unsigned ir = 0; ir < nen_; ++ir)
        {
          const double tmp_val = funct_m(ir);

          for (unsigned ivel = 0; ivel < nsd_; ++ivel)
          {
            // +(stab * vm, u_DBC) (weak dirichlet case) or from
            // +(stab * vm, u_s)
            rhC_um_(mIndex(ir, ivel), 0) -= tmp_val * (velint_diff_timefacfac_stabfac_(ivel));
          }
        }


        if (s_col.first)
        {
          // - gamma*mu/h_K (vm, us))
          // - gamma*mu/h_K (vs, um))

          const double stabfac_timefacfac_ms = timefacfac * m_row.second * s_col.second;

          for (unsigned ic = 0; ic < slave_nen_; ++ic)
          {
            for (unsigned ir = 0; ir < nen_; ++ir)
            {
              const double tmp_val = funct_s_m_dyad_(ic, ir) * stabfac_timefacfac_ms;

              for (unsigned ivel = 0; ivel < nsd_; ++ivel)
              {
                for (unsigned jvel = 0; jvel < nsd_; ++jvel)
                {
                  C_umus_(mIndex(ir, jvel), sIndex(ic, ivel)) -= tmp_val * proj_matrix_(ivel, jvel);
                }
              }
            }
          }
        }

        if (s_row.first && s_col.first)
        {
          const double stabfac_timefacfac_ss = timefacfac * s_row.second * s_col.second;

          for (unsigned ic = 0; ic < slave_nen_; ++ic)
          {
            // + gamma*mu/h_K (vs, us))
            for (unsigned ir = 0; ir < slave_nen_; ++ir)
            {
              const double tmp_val = funct_s_s_dyad_(ir, ic) * stabfac_timefacfac_ss;

              for (unsigned ivel = 0; ivel < nsd_; ++ivel)
              {
                for (unsigned jvel = 0; jvel < nsd_; ++jvel)
                {
                  C_usus_(sIndex(ir, jvel), sIndex(ic, ivel)) += tmp_val * proj_matrix_(ivel, jvel);
                }
              }
            }
          }
        }

        if (s_row.first)
        {
          const double stabfac_timefacfac_sm = timefacfac * s_row.second * m_col.second;

          for (unsigned ic = 0; ic < nen_; ++ic)
          {
            for (unsigned ir = 0; ir < slave_nen_; ++ir)
            {
              const double tmp_val = funct_s_m_dyad_(ir, ic) * stabfac_timefacfac_sm;

              for (unsigned ivel = 0; ivel < nsd_; ++ivel)
              {
                for (unsigned jvel = 0; jvel < nsd_; ++jvel)
                {
                  C_usum_(sIndex(ir, jvel), mIndex(ic, ivel)) -= tmp_val * proj_matrix_(ivel, jvel);
                }
              }
            }
          }

          const double stabfac_timefacfac_s = timefacfac * s_row.second;
          velint_diff_timefacfac_stabfac_.Update(
              stabfac_timefacfac_s, velint_diff_proj_matrix_, 0.0);

          for (unsigned ir = 0; ir < slave_nen_; ++ir)
          {
            const double tmp_val = funct_s_(ir);

            for (unsigned ivel = 0; ivel < nsd_; ++ivel)
            {
              // +(stab * vs, um)
              // -(stab * vs, us)
              rhC_us_(sIndex(ir, ivel), 0) += tmp_val * velint_diff_timefacfac_stabfac_(ivel);
            }
          }
        }

        return;
      }



      /*----------------------------------------------------------------------*
       * add averaged term to balance instabilities due to convective
       * mass transport across the fluid-fluid interface
       *----------------------------------------------------------------------*/
      template <DRT::Element::DiscretizationType distype,
          DRT::Element::DiscretizationType slave_distype, unsigned int slave_numdof>
      void NitscheCoupling<distype, slave_distype, slave_numdof>::NIT_Stab_Inflow_AveragedTerm(
          const LINALG::Matrix<nen_, 1>& funct_m,   ///< funct
          const LINALG::Matrix<nsd_, 1>& velint_m,  ///< master velocity
          const LINALG::Matrix<nsd_, 1>& normal,    ///< normal vector n^m
          const double& density,                    ///< fluid density
          const double& timefacfac,                 ///< timefac * fac
          bool only_rhs)
      {
        //
        /*                                        \
       -|  [rho * (beta * n)] *  { v }_m , [   u ] |
        \ ----stab_avg-----                       / */

        // { v }_m = 0.5* (v^b + v^e) leads to the scaling with 0.5;
        // beta: convective velocity, currently beta=u^b_Gamma;
        // n:= n^b
        const double stabfac_avg_scaled = 0.5 * velint_m.Dot(normal) * density * timefacfac;

        for (unsigned ivel = 0; ivel < nsd_; ++ivel)
        {
          for (unsigned ir = 0; ir < nen_; ++ir)
          {
            const unsigned mrow = mIndex(ir, ivel);

            const double tmp = funct_m(ir) * stabfac_avg_scaled;
            rhC_um_(mrow, 0) += tmp * velint_diff_(ivel);
          }

          for (unsigned ir = 0; ir < slave_nen_; ++ir)
          {
            const unsigned srow = sIndex(ir, ivel);

            const double tmp = funct_s_(ir) * stabfac_avg_scaled;
            rhC_us_(srow, 0) += tmp * velint_diff_(ivel);
          }
        }

        if (only_rhs) return;

        for (unsigned ivel = 0; ivel < nsd_; ++ivel)
        {
          //  [rho * (beta * n^b)] (0.5*vb,ub)
          for (unsigned ir = 0; ir < nen_; ++ir)
          {
            const unsigned mrow = mIndex(ir, ivel);

            for (unsigned ic = 0; ic < nen_; ++ic)
            {
              C_umum_(mrow, mIndex(ic, ivel)) -= funct_m_m_dyad_(ir, ic) * stabfac_avg_scaled;
            }

            //  -[rho * (beta * n^b)] (0.5*vb,ue)
            for (unsigned ic = 0; ic < slave_nen_; ++ic)
            {
              C_umus_(mrow, sIndex(ic, ivel)) += funct_s_m_dyad_(ic, ir) * stabfac_avg_scaled;
            }
          }

          //  [rho * (beta * n^b)] (0.5*ve,ub)
          for (unsigned ir = 0; ir < slave_nen_; ++ir)
          {
            const unsigned srow = sIndex(ir, ivel);

            for (unsigned ic = 0; ic < nen_; ++ic)
            {
              C_usum_(srow, mIndex(ic, ivel)) -= funct_s_m_dyad_(ir, ic) * stabfac_avg_scaled;
            }

            //-[rho * (beta * n^b)] (0.5*ve,ue)
            for (unsigned ic = 0; ic < slave_nen_; ++ic)
            {
              C_usus_(srow, sIndex(ic, ivel)) += funct_s_s_dyad_(ir, ic) * stabfac_avg_scaled;
            }
          }
        }
        return;
      }

      /*----------------------------------------------------------------------*
       *----------------------------------------------------------------------*/
      template <DRT::Element::DiscretizationType distype,
          DRT::Element::DiscretizationType slave_distype, unsigned int slave_numdof>
      void NitscheCoupling<distype, slave_distype,
          slave_numdof>::NIT_Create_Standard_Projection_Matrices(const LINALG::Matrix<nsd_, 1>&
              normal  ///< normal vector n^b
      )
      {
        // Create the identity matrix (Probably not the fastest way...) Might make it global?
        proj_tangential_.Scale(0.0);
        for (unsigned int i = 0; i < nsd_; ++i) proj_tangential_(i, i) = 1;

        //   Non-smoothed projection matrix
        for (unsigned int i = 0; i < nsd_; ++i)
        {
          for (unsigned int j = 0; j < nsd_; ++j)
          {
            proj_tangential_(i, j) -= normal(i, 0) * normal(j, 0);
          }
        }

        proj_normal_.Scale(0.0);
        for (unsigned i = 0; i < nsd_; ++i) proj_normal_(i, i) = 1.0;
        proj_normal_.Update(-1.0, proj_tangential_, 1.0);
      }

      /*----------------------------------------------------------------------*
       *----------------------------------------------------------------------*/
      template <DRT::Element::DiscretizationType distype,
          DRT::Element::DiscretizationType slave_distype, unsigned int slave_numdof>
      void HybridLMCoupling<distype, slave_distype, slave_numdof>::MHCS_buildCouplingMatrices(
          const LINALG::Matrix<nsd_, 1>& normal,  ///< normal vector
          const double& fac,                      ///< integration factor
          const LINALG::Matrix<nen_, 1>& funct,   ///< shape function
          LINALG::BlockMatrix<LINALG::Matrix<nen_, 1>, numstressdof_, 1>&
              rhs_s,  ///< block rhs vector \f$ rhs_{\sigma} \f$
          const LINALG::Matrix<nsd_, 1>&
              ivelint_jump,  ///< prescribed interface velocity or interface jump height
          const LINALG::Matrix<nsd_, 1>&
              itraction_jump  ///< prescribed interface traction or interface jump height
      )
      {
        LINALG::Matrix<nen_, slave_nen_> bK_ms;
        LINALG::Matrix<slave_nen_, nen_> bK_sm;

        // interface velocity at gauss-point (current gauss-point in calling method)
        LINALG::Matrix<nsd_, 1> velint_s(true);
        this->GetInterfaceVelnp(velint_s);

        // add the prescribed interface velocity for weak Dirichlet boundary conditions or the jump
        // height for coupled problems
        velint_s.Update(1.0, ivelint_jump, 1.0);

        // get nodal shape function vector
        LINALG::Matrix<slave_nen_, 1> slave_funct(true);
        this->GetSlaveFunct(slave_funct);

        bK_ms.MultiplyNT(funct, slave_funct);
        bK_sm.UpdateT(bK_ms);

        for (unsigned ivel = 0; ivel < nsd_; ++ivel)
        {
          for (unsigned jvel = 0; jvel < nsd_; ++jvel)
          {
            const double tmp = fac * normal(jvel);

            const unsigned sigma = stressIndex(ivel, jvel);
            // G_sus
            BG_sus_(sigma, ivel)->Update(tmp, bK_ms, 1.0);
            rhs_s(sigma, 0)->Update(-tmp * velint_s(ivel), funct, 1.0);

            // G_uss
            BG_uss_(ivel, sigma)->Update(tmp, bK_sm, 1.0);
          }
        }

        const double km = 1.0;                                // only master-sided weighting
        LINALG::Matrix<slave_nen_, 1> funct_s_timefacfac_km;  ///< funct_s * timefacfac *kappa_m
        funct_s_timefacfac_km.Update(km, slave_funct, 0.0);

        // Traction Standard Consistency term
        MH_Traction_Consistency_Term(funct_s_timefacfac_km, itraction_jump);
      }

      /*----------------------------------------------------------------------*
       *----------------------------------------------------------------------*/
      template <DRT::Element::DiscretizationType distype,
          DRT::Element::DiscretizationType slave_distype, unsigned int slave_numdof>
      void HybridLMCoupling<distype, slave_distype, slave_numdof>::MHVS_buildCouplingMatrices(
          const LINALG::Matrix<nsd_, 1>& normal,  ///< normal vector
          const double& fac,                      ///< integration factor
          const LINALG::Matrix<nen_, 1>& funct,   ///< background element shape functions
          LINALG::BlockMatrix<LINALG::Matrix<nen_, 1>, numstressdof_, 1>&
              rhs_s,                          ///< block rhs vector \f$ rhs_{\sigma}\f$
          const double& press,                ///< background element pressure
          LINALG::Matrix<nen_, 1>& rhs_pmus,  ///< part of block rhs vector \f$rhs_p\f$ including
                                              ///< interface velocity terms
          const LINALG::Matrix<nsd_, 1>&
              ivelint_jump,  ///< prescribed interface velocity or interface jump height
          const LINALG::Matrix<nsd_, 1>&
              itraction_jump  ///< prescribed interface traction or interface jump height
      )
      {
        // interface velocity at gauss-point (current gauss-point in calling method)
        LINALG::Matrix<nsd_, 1> velint_s(true);
        this->GetInterfaceVelnp(velint_s);

        // add the prescribed interface velocity for weak Dirichlet boundary conditions or the jump
        // height for coupled problems
        velint_s.Update(1.0, ivelint_jump, 1.0);

        // block submatrices for interface coupling; s: slave side, m: master side (always
        // background here)
        LINALG::Matrix<nen_, slave_nen_> bG_ms(true);
        LINALG::Matrix<slave_nen_, nen_> bG_sm(true);

        LINALG::Matrix<slave_nen_, 1> slave_funct;
        this->GetSlaveFunct(slave_funct);

        bG_ms.MultiplyNT(funct, slave_funct);
        bG_sm.MultiplyNT(slave_funct, funct);

        for (unsigned ivel = 0; ivel < nsd_; ++ivel)
        {
          for (unsigned jvel = 0; jvel < nsd_; ++jvel)
          {
            /*
             * G_sus
             *
             * from:
             *   /              \
             *  | \tau^m n, u^s  |
             *   \              /
             */
            BG_sus_(stressIndex(ivel, jvel), ivel)->Update(fac * normal(jvel), bG_ms, 1.0);
            rhs_s(stressIndex(ivel, jvel), 0)
                ->Update(-fac * normal(jvel) * velint_s(ivel), funct, 1.0);

            /*
             *  G_uss
             *
             *  from:
             *   /                 \
             *  | v^s, \sigma^m n   |
             *   \                 /
             *
             */
            BG_uss_(ivel, stressIndex(ivel, jvel))->Update(fac * normal(jvel), bG_sm, 1.0);
          }

          // Build cross-interface pressure-velocity coupling matrices G_uip, G_pui!

          // G_pmus - from adjoint pressure consistency term
          /*
           *  /            \
           *  | q^m, u^s n  |
           *  \            /
           *
           */

          BG_pmus_(0, ivel)->Update(fac * normal(ivel), bG_ms, 1.0);


          // G_uspm - from pressure consistency term
          /*
           *  /           \
           *  | -v^s, p n  |
           *  \           /
           *
           */

          BG_uspm_(ivel, 0)->Update(-fac * normal(ivel), bG_sm, 1.0);
        }

        // add normal interface velocity to rhs vector (pressure row)
        const double svelnormal = velint_s.Dot(normal);

        // (Viscous) stress/Velocities

        // rhs_pm_us, residual from:
        /*
         *  /            \
         *  | q^m, u^s n  |
         *  \            /
         *
         */
        rhs_pmus.Update(-fac * svelnormal, funct, 1.0);

        // rhs_us_pm, residual from:
        /*
         *  /            \
         *  | -v^s, p^m n |
         *  \            /
         *
         */
        // belongs to the side and therefore contributes to rhC_ui_!
        for (unsigned ir = 0; ir < slave_nen_; ++ir)
        {
          for (unsigned ivel = 0; ivel < nsd_; ++ivel)
          {
            rhC_us_(sIndex(ir, ivel), 0) += press * fac * normal(ivel) * slave_funct(ir);
          }

          rhC_us_(sIndex(ir, Pres), 0) = 0.0;
        }

        const double km = 1.0;                                // only master-sided weighting
        LINALG::Matrix<slave_nen_, 1> funct_s_timefacfac_km;  ///< funct_s * timefacfac *kappa_m
        funct_s_timefacfac_km.Update(km, slave_funct, 0.0);

        // Traction Standard Consistency term
        MH_Traction_Consistency_Term(funct_s_timefacfac_km, itraction_jump);


        return;
      }


      /*----------------------------------------------------------------------*
       *----------------------------------------------------------------------*/
      template <DRT::Element::DiscretizationType distype,
          DRT::Element::DiscretizationType slave_distype, unsigned int slave_numdof>
      void HybridLMCoupling<distype, slave_distype, slave_numdof>::MH_Traction_Consistency_Term(
          const LINALG::Matrix<slave_nen_, 1>&
              funct_s_timefacfac_km,  ///< funct_s * timefacfac *kappa_m
          const LINALG::Matrix<nsd_, 1>&
              itraction_jump  ///< prescribed interface traction, jump height for coupled problems
      )
      {
        /*            \
     - |  < v >,   t   |   with t = [sigma * n]
        \             /     */

        // loop over velocity components
        for (unsigned ivel = 0; ivel < nsd_; ++ivel)
        {
          /*
          //-----------------------------------------------
          //    - (vm, ks * t) = 0 as ks=0
          //-----------------------------------------------
          */

          //-----------------------------------------------
          //    + (vs, km * t)
          //-----------------------------------------------
          for (unsigned ir = 0; ir < slave_nen_; ++ir)
          {
            const double funct_s_km_timefacfac_traction =
                funct_s_timefacfac_km(ir) * itraction_jump(ivel);

            const unsigned row = sIndex(ir, ivel);
            rhC_us_(row, 0) += funct_s_km_timefacfac_traction;
          }
        }  // end loop over velocity components
      }


      /*----------------------------------------------------------------------*
       *----------------------------------------------------------------------*/
      template <DRT::Element::DiscretizationType distype,
          DRT::Element::DiscretizationType slave_distype, unsigned int slave_numdof>
      void
      HybridLMCoupling<distype, slave_distype, slave_numdof>::HybridLM_buildFinalCouplingMatrices(
          LINALG::BlockMatrix<LINALG::Matrix<nen_, nen_>, numstressdof_, numstressdof_>&
              BinvK_ss,  ///< block inverse \f$ K^{-1}_{\sigma\sigma} \f$
          LINALG::BlockMatrix<LINALG::Matrix<nen_, nen_>, master_numdof_, numstressdof_>&
              BKumsInvKss,  ///< block matrix \f$ K_{u\sigma} \cdot K^{-1}_{\sigma\sigma} \f$
          LINALG::BlockMatrix<LINALG::Matrix<nen_, nen_>, numstressdof_, master_numdof_>&
              BK_sum,  ///< block matrix \f$ K_{\sigma u} \f$
          LINALG::BlockMatrix<LINALG::Matrix<nen_, 1>, numstressdof_, 1>&
              rhs_s  ///< block rhs vector \f$ rhs_{\sigma}\f$
      )
      {
        // final coupling matrices in block form
        LINALG::BlockMatrix<LINALG::Matrix<nen_, slave_nen_>, master_numdof_, nsd_> BCumus;
        LINALG::BlockMatrix<LINALG::Matrix<slave_nen_, nen_>, nsd_, master_numdof_> BCusum;

        // auxiliary matrices for intermediate calculation results of the matrix product
        LINALG::BlockMatrix<LINALG::Matrix<slave_nen_, nen_>, nsd_, numstressdof_> BGussInvKss;
        LINALG::BlockMatrix<LINALG::Matrix<slave_nen_, 1>, nsd_, 1> BGussinvKssrhs_s;

        // BKusInvKss:
        // (K_ums + G_ums) * K_ss^-1 (MHVS) or
        // G_ums * K_ss^-1 (MHCS)

        // (K_ums + G_ums) * K_ss^-1 G_sus (MHVS) or
        // G_ums * K_ss^-1 G_sus (MHCS)
        BCumus.Multiply(BKumsInvKss, BG_sus_);

        // G_uss * K_ss^-1
        BGussInvKss.Multiply(BG_uss_, BinvK_ss);

        // G_uss * K_ss^-1 * (K_sum + G_sum) (MHVS) or
        // G_uss * K_ss^-1 * (K_sum + G_sum + K_spm) (MHCS)
        BCusum.Multiply(BGussInvKss, BK_sum);

        // G_uss K_ss^-1 rhs_s
        BGussinvKssrhs_s.Multiply(BGussInvKss, rhs_s);

        // transfer the entries from Cumus,Cusum, rhCus (in case of MHCS, coupling term us-p is
        // included in Cusum!)

        // loop over slave velocity dof
        for (unsigned isvel = 0; isvel < nsd_; ++isvel)
        {
          // loop over background element dof (velocity & pressure)
          for (unsigned imdof = 0; imdof < master_numdof_; ++imdof)
          {
            // (um-us)
            if (BCumus.IsUsed(imdof, isvel))
            {
              const LINALG::Matrix<nen_, slave_nen_>& bCumus = *BCumus(imdof, isvel);
              // loop over slave element nodes
              for (unsigned isn = 0; isn < slave_nen_; ++isn)
              {
                // loop over background element nodes
                for (unsigned imn = 0; imn < nen_; ++imn)
                {
                  C_umus_(mIndex(imn, imdof), sIndex(isn, isvel)) -= bCumus(imn, isn);
                }
              }
            }  // (um-us)

            // (us-um), MHCS: (us-pm)
            if (BCusum.IsUsed(isvel, imdof))
            {
              const LINALG::Matrix<slave_nen_, nen_>& bCusum = *BCusum(isvel, imdof);
              // loop over slave element nodes
              for (unsigned isn = 0; isn < slave_nen_; ++isn)
              {
                // loop over background element nodes
                for (unsigned imn = 0; imn < nen_; ++imn)
                {
                  C_usum_(sIndex(isn, isvel), mIndex(imn, imdof)) -= bCusum(isn, imn);
                }
              }
            }  // (us-um), MHCS: (us-pm)
          }

          // add surface-based pressure coupling terms (only MHVS)

          // (us-pm)
          if (BG_uspm_.IsUsed(isvel, 0))
          {
            const LINALG::Matrix<slave_nen_, nen_>& bGuspm = *BG_uspm_(isvel, 0);
            // loop over slave element nodes
            for (unsigned isn = 0; isn < slave_nen_; ++isn)
            {
              // loop over background element nodes
              for (unsigned imn = 0; imn < nen_; ++imn)
              {
                C_usum_(sIndex(isn, isvel), mIndex(imn, Pres)) = bGuspm(isn, imn);
              }
            }
          }  // (us-pm)

          // (pm-us)
          if (BG_pmus_.IsUsed(0, isvel))
          {
            const LINALG::Matrix<nen_, slave_nen_>& bGpmus = *BG_pmus_(0, isvel);
            // loop over slave element nodes
            for (unsigned isn = 0; isn < slave_nen_; ++isn)
            {
              // loop over background element nodes
              for (unsigned imn = 0; imn < nen_; ++imn)
              {
                C_umus_(mIndex(imn, Pres), sIndex(isn, isvel)) = bGpmus(imn, isn);
              }
            }
          }  // (pm-us)

          // rhs - us
          if (BGussinvKssrhs_s.IsUsed(isvel, 0))
          {
            const LINALG::Matrix<slave_nen_, 1>& bGussinvKssrhs_s = *BGussinvKssrhs_s(isvel, 0);

            // loop over slave element nodes
            for (unsigned isn = 0; isn < slave_nen_; ++isn)
            {
              rhC_us_(sIndex(isn, isvel), 0) -= bGussinvKssrhs_s(isn, 0);
            }
          }  // rhs - us
        }    // end loop over slave velocity dof

        // finally, build G_uss & G_sus for C_usus

        // G_sus
        // loop over block rows
        for (unsigned ibr = 0; ibr < numstressdof_; ++ibr)
        {
          // loop over block columns (interface velocity)
          for (unsigned ibc = 0; ibc < nsd_; ++ibc)
          {
            // extract the stress-velocity coupling submatrix
            if (BG_sus_.IsUsed(ibr, ibc))
            {
              LINALG::Matrix<nen_, slave_nen_>& bGsus = *BG_sus_(ibr, ibc);

              // transfer the entries
              for (unsigned ir = 0; ir < nen_; ++ir)
              {
                unsigned stressrow = ibr + ir * numstressdof_;

                for (unsigned ic = 0; ic < slave_nen_; ++ic)
                {
                  unsigned slavevelcol = ibc + ic * slave_numdof;

                  G_sus_(stressrow, slavevelcol) = bGsus(ir, ic);
                }
              }
            }
          }
        }  // G_sus filled

        // fill G_uss_ from BG_uss_

        // loop over block columns (sigmaxx, sigmaxy, ...)
        for (unsigned ibc = 0; ibc < numstressdof_; ++ibc)
        {
          // loop over block rows (interface velocity)
          for (unsigned ibr = 0; ibr < nsd_; ++ibr)
          {
            if (BG_uss_.IsUsed(ibr, ibc))
            {
              LINALG::Matrix<slave_nen_, nen_>& bGuss = *BG_uss_(ibr, ibc);

              // transfer the entries
              for (unsigned ic = 0; ic < nen_; ++ic)
              {
                unsigned stresscol = ibc + ic * numstressdof_;

                for (unsigned ir = 0; ir < slave_nen_; ++ir)
                {
                  unsigned slavevelrow = ibr + ir * slave_numdof;

                  G_uss_(slavevelrow, stresscol) = bGuss(ir, ic);
                }
              }
            }
          }
        }  // G_uss filled
      }

      /*----------------------------------------------------------------------*
       *----------------------------------------------------------------------*/
      template <DRT::Element::DiscretizationType distype,
          DRT::Element::DiscretizationType slave_distype, unsigned int slave_numdof>
      void NitscheCoupling<distype, slave_distype, slave_numdof>::
          Do_NIT_visc_Adjoint_and_Neumann_MasterTerms_Projected(
              const LINALG::Matrix<nen_, 1>& funct_m,  ///< funct * timefacfac
              const LINALG::Matrix<nsd_, nen_>&
                  derxy_m,  ///< spatial derivatives of coupling master shape functions
              const LINALG::Matrix<nsd_, nsd_>&
                  vderxy_m,  ///< coupling master spatial velocity derivatives
              const LINALG::Matrix<nsd_, nsd_>& projection_matrix,  ///< projection_matrix
              const LINALG::Matrix<nsd_, 1>&
                  velint_diff_proj_matrix,             ///< velocity difference projected
              const LINALG::Matrix<nsd_, 1>& normal,   ///< normal-vector
              const double& km_viscm_fac,              ///< scaling factor
              const std::pair<bool, double>& m_row,    ///< scaling for master row
              const std::pair<bool, double>& m_col,    ///< scaling for master col
              const std::pair<bool, double>& s_col,    ///< scaling for slave col
              const std::pair<bool, double>& mstr_col  ///< scaling for master stress col
          )
      {
        velint_diff_proj_matrix_ = velint_diff_proj_matrix;
        proj_matrix_ = projection_matrix;

        // 2.0 * timefacfac * visceff_m * k_m * [\nabla N^(IX)]_k P^t_{kj}
        proj_matrix_derxy_m_.MultiplyTN(proj_matrix_, derxy_m);  // Apply from right for consistency
        proj_matrix_derxy_m_.Scale(km_viscm_fac);

        //-----------------------------------------------------------------
        // viscous adjoint consistency term

        double tmp_fac = adj_visc_scale_ * km_viscm_fac;
        derxy_m_viscm_timefacfac_.Update(tmp_fac, derxy_m);  // 2 * mu_m * timefacfac *
                                                             // derxy_m(k,ic)

        // Scale with adjoint viscous scaling {-1,+1}
        proj_matrix_derxy_m_.Scale(adj_visc_scale_);

        // Same as half_normal_deriv_m_viscm_timefacfac_km_? Might be unnecessary?
        // normal_deriv_m_viscm_km_ = alpha * half_normal(k)* 2 * km * mu_m * timefacfac *
        // derxy_m(k,IX)
        //                          = alpha * mu_m * viscfac_km * c(IX)
        normal_deriv_m_viscm_km_.MultiplyTN(derxy_m_viscm_timefacfac_, half_normal_);

        NIT_visc_AdjointConsistency_MasterTerms_Projected(
            derxy_m_viscm_timefacfac_, funct_m, normal, m_row, m_col, s_col);

#ifndef ENFORCE_URQUIZA_GNBC
        //-----------------------------------------------------------------
        // Terms needed for Neumann consistency terms

        if (mstr_col.first)
        {
          NIT_visc_Neumann_AdjointConsistency_MasterTerms_Projected(
              derxy_m_viscm_timefacfac_, derxy_m, vderxy_m, funct_m, normal, m_row, mstr_col);
        }
        //-----------------------------------------------------------------
#endif  // ENFORCE_URQUIZA_GNBC

        return;
      }


    }  // namespace XFLUID
  }    // namespace ELEMENTS
}  // namespace DRT

// pairs with numdof=3
// template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::hex8,  DRT::Element::tri3,3>;
// template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::hex8,  DRT::Element::tri6,3>;
template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::hex8, DRT::Element::quad4, 3>;
template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::hex8, DRT::Element::quad8, 3>;
template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::hex8, DRT::Element::quad9, 3>;
// template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::hex20, DRT::Element::tri3,3>;
// template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::hex20, DRT::Element::tri6,3>;
template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::hex20, DRT::Element::quad4, 3>;
template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::hex20, DRT::Element::quad8, 3>;
template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::hex20, DRT::Element::quad9, 3>;
// template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::hex27, DRT::Element::tri3,3>;
// template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::hex27, DRT::Element::tri6,3>;
template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::hex27, DRT::Element::quad4, 3>;
template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::hex27, DRT::Element::quad8, 3>;
template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::hex27, DRT::Element::quad9, 3>;
// template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::tet4,  DRT::Element::tri3,3>;
// template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::tet4,  DRT::Element::tri6,3>;
template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::tet4, DRT::Element::quad4, 3>;
template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::tet4, DRT::Element::quad8, 3>;
template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::tet4, DRT::Element::quad9, 3>;
// template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::tet10, DRT::Element::tri3,3>;
// template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::tet10, DRT::Element::tri6,3>;
template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::tet10, DRT::Element::quad4, 3>;
template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::tet10, DRT::Element::quad8, 3>;
template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::tet10, DRT::Element::quad9, 3>;
// template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::wedge6,
// DRT::Element::tri3,3>; template class
// DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::wedge6, DRT::Element::tri6,3>;
template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::wedge6, DRT::Element::quad4, 3>;
template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::wedge6, DRT::Element::quad8, 3>;
template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::wedge6, DRT::Element::quad9, 3>;
// template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::wedge15,
// DRT::Element::tri3,3>; template class
// DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::wedge15,  DRT::Element::tri6,3>;
template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::wedge15, DRT::Element::quad4,
    3>;
template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::wedge15, DRT::Element::quad8,
    3>;
template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::wedge15, DRT::Element::quad9,
    3>;

template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::hex8, DRT::Element::dis_none,
    3>;
template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::hex20, DRT::Element::dis_none,
    3>;
template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::hex27, DRT::Element::dis_none,
    3>;
template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::tet4, DRT::Element::dis_none,
    3>;
template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::tet10, DRT::Element::dis_none,
    3>;
template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::wedge6, DRT::Element::dis_none,
    3>;
template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::wedge15, DRT::Element::dis_none,
    3>;

// volume coupled with numdof = 3, FSI Slavesided, FPI
template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::hex8, DRT::Element::hex8, 3>;
template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::hex20, DRT::Element::hex8, 3>;
template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::hex27, DRT::Element::hex8, 3>;
template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::tet4, DRT::Element::hex8, 3>;
template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::tet10, DRT::Element::hex8, 3>;
template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::wedge6, DRT::Element::hex8, 3>;
template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::wedge15, DRT::Element::hex8, 3>;

// pairs with numdof=4
// template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::hex8,  DRT::Element::tri3,4>;
// template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::hex8,  DRT::Element::tri6,4>;
template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::hex8, DRT::Element::quad4, 4>;
template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::hex8, DRT::Element::quad8, 4>;
template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::hex8, DRT::Element::quad9, 4>;
// template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::hex20, DRT::Element::tri3,4>;
// template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::hex20, DRT::Element::tri6,4>;
template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::hex20, DRT::Element::quad4, 4>;
template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::hex20, DRT::Element::quad8, 4>;
template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::hex20, DRT::Element::quad9, 4>;
// template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::hex27, DRT::Element::tri3,4>;
// template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::hex27, DRT::Element::tri6,4>;
template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::hex27, DRT::Element::quad4, 4>;
template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::hex27, DRT::Element::quad8, 4>;
template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::hex27, DRT::Element::quad9, 4>;
// template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::tet4,  DRT::Element::tri3,4>;
// template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::tet4,  DRT::Element::tri6,4>;
template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::tet4, DRT::Element::quad4, 4>;
template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::tet4, DRT::Element::quad8, 4>;
template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::tet4, DRT::Element::quad9, 4>;
// template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::tet10, DRT::Element::tri3,4>;
// template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::tet10, DRT::Element::tri6,4>;
template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::tet10, DRT::Element::quad4, 4>;
template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::tet10, DRT::Element::quad8, 4>;
template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::tet10, DRT::Element::quad9, 4>;
// template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::wedge6,
// DRT::Element::tri3,4>; template class
// DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::wedge6, DRT::Element::tri6,4>;
template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::wedge6, DRT::Element::quad4, 4>;
template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::wedge6, DRT::Element::quad8, 4>;
template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::wedge6, DRT::Element::quad9, 4>;
// template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::wedge15,
// DRT::Element::tri3,4>; template class
// DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::wedge15, DRT::Element::tri6,4>;
template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::wedge15, DRT::Element::quad4,
    4>;
template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::wedge15, DRT::Element::quad8,
    4>;
template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::wedge15, DRT::Element::quad9,
    4>;
//

// template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::hex8,  DRT::Element::tet4,
// 4>; template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::hex8,
// DRT::Element::tet10,4>;
template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::hex8, DRT::Element::hex8, 4>;
template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::hex8, DRT::Element::hex20, 4>;
template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::hex8, DRT::Element::hex27, 4>;
// template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::hex8,
// DRT::Element::wedge15,4>; template class
// DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::hex20, DRT::Element::tet4,4>; template class
// DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::hex20, DRT::Element::tet10,4>;
template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::hex20, DRT::Element::hex8, 4>;
template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::hex20, DRT::Element::hex20, 4>;
template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::hex20, DRT::Element::hex27, 4>;
// template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::hex20,
// DRT::Element::wedge15,4>; template class
// DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::hex27, DRT::Element::tet4,4>; template class
// DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::hex27, DRT::Element::tet10,4>;
template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::hex27, DRT::Element::hex8, 4>;
template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::hex27, DRT::Element::hex20, 4>;
template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::hex27, DRT::Element::hex27, 4>;
// template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::hex27,
// DRT::Element::wedge15,4>; template class
// DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::tet4,  DRT::Element::tet4,4>; template class
// DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::tet4,  DRT::Element::tet10,4>;
template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::tet4, DRT::Element::hex8, 4>;
template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::tet4, DRT::Element::hex20, 4>;
template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::tet4, DRT::Element::hex27, 4>;
// template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::tet4,
// DRT::Element::wedge15,4>; template class
// DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::tet10, DRT::Element::tet4,4>; template class
// DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::tet10, DRT::Element::tet10,4>;
template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::tet10, DRT::Element::hex8, 4>;
template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::tet10, DRT::Element::hex20, 4>;
template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::tet10, DRT::Element::hex27, 4>;
// template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::tet10,
// DRT::Element::wedge15,4>; template class
// DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::wedge6, DRT::Element::tet4,4>; template
// class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::wedge6, DRT::Element::tet10,4>;
template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::wedge6, DRT::Element::hex8, 4>;
template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::wedge6, DRT::Element::hex20, 4>;
template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::wedge6, DRT::Element::hex27, 4>;
// template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::wedge6,
// DRT::Element::wedge15,4>; template class
// DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::wedge15, DRT::Element::tet4,4>; template
// class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::wedge15, DRT::Element::tet10,4>;
template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::wedge15, DRT::Element::hex8, 4>;
template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::wedge15, DRT::Element::hex20,
    4>;
template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::wedge15, DRT::Element::hex27,
    4>;
// template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::wedge15,
// DRT::Element::wedge15,4>;


template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::hex8, DRT::Element::dis_none,
    4>;
template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::hex20, DRT::Element::dis_none,
    4>;
template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::hex27, DRT::Element::dis_none,
    4>;
template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::tet4, DRT::Element::dis_none,
    4>;
template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::tet10, DRT::Element::dis_none,
    4>;
template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::wedge6, DRT::Element::dis_none,
    4>;
template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::wedge15, DRT::Element::dis_none,
    4>;


// pairs with numdof=3
// template class DRT::ELEMENTS::XFLUID::HybridLMCoupling<DRT::Element::hex8, DRT::Element::tri3,3>;
// template class DRT::ELEMENTS::XFLUID::HybridLMCoupling<DRT::Element::hex8, DRT::Element::tri6,3>;
template class DRT::ELEMENTS::XFLUID::HybridLMCoupling<DRT::Element::hex8, DRT::Element::quad4, 3>;
template class DRT::ELEMENTS::XFLUID::HybridLMCoupling<DRT::Element::hex8, DRT::Element::quad8, 3>;
template class DRT::ELEMENTS::XFLUID::HybridLMCoupling<DRT::Element::hex8, DRT::Element::quad9, 3>;
// template class DRT::ELEMENTS::XFLUID::HybridLMCoupling<DRT::Element::hex20,
// DRT::Element::tri3,3>; template class
// DRT::ELEMENTS::XFLUID::HybridLMCoupling<DRT::Element::hex20, DRT::Element::tri6,3>;
template class DRT::ELEMENTS::XFLUID::HybridLMCoupling<DRT::Element::hex20, DRT::Element::quad4, 3>;
template class DRT::ELEMENTS::XFLUID::HybridLMCoupling<DRT::Element::hex20, DRT::Element::quad8, 3>;
template class DRT::ELEMENTS::XFLUID::HybridLMCoupling<DRT::Element::hex20, DRT::Element::quad9, 3>;
// template class DRT::ELEMENTS::XFLUID::HybridLMCoupling<DRT::Element::hex27,
// DRT::Element::tri3,3>; template class
// DRT::ELEMENTS::XFLUID::HybridLMCoupling<DRT::Element::hex27, DRT::Element::tri6,3>;
template class DRT::ELEMENTS::XFLUID::HybridLMCoupling<DRT::Element::hex27, DRT::Element::quad4, 3>;
template class DRT::ELEMENTS::XFLUID::HybridLMCoupling<DRT::Element::hex27, DRT::Element::quad8, 3>;
template class DRT::ELEMENTS::XFLUID::HybridLMCoupling<DRT::Element::hex27, DRT::Element::quad9, 3>;
// template class DRT::ELEMENTS::XFLUID::HybridLMCoupling<DRT::Element::tet4, DRT::Element::tri3,3>;
// template class DRT::ELEMENTS::XFLUID::HybridLMCoupling<DRT::Element::tet4, DRT::Element::tri6,3>;
template class DRT::ELEMENTS::XFLUID::HybridLMCoupling<DRT::Element::tet4, DRT::Element::quad4, 3>;
template class DRT::ELEMENTS::XFLUID::HybridLMCoupling<DRT::Element::tet4, DRT::Element::quad8, 3>;
template class DRT::ELEMENTS::XFLUID::HybridLMCoupling<DRT::Element::tet4, DRT::Element::quad9, 3>;
// template class DRT::ELEMENTS::XFLUID::HybridLMCoupling<DRT::Element::tet10,
// DRT::Element::tri3,3>; template class
// DRT::ELEMENTS::XFLUID::HybridLMCoupling<DRT::Element::tet10, DRT::Element::tri6,3>;
template class DRT::ELEMENTS::XFLUID::HybridLMCoupling<DRT::Element::tet10, DRT::Element::quad4, 3>;
template class DRT::ELEMENTS::XFLUID::HybridLMCoupling<DRT::Element::tet10, DRT::Element::quad8, 3>;
template class DRT::ELEMENTS::XFLUID::HybridLMCoupling<DRT::Element::tet10, DRT::Element::quad9, 3>;
// template class DRT::ELEMENTS::XFLUID::HybridLMCoupling<DRT::Element::wedge6,
// DRT::Element::tri3,3>; template class
// DRT::ELEMENTS::XFLUID::HybridLMCoupling<DRT::Element::wedge6,  DRT::Element::tri6,3>;
template class DRT::ELEMENTS::XFLUID::HybridLMCoupling<DRT::Element::wedge6, DRT::Element::quad4,
    3>;
template class DRT::ELEMENTS::XFLUID::HybridLMCoupling<DRT::Element::wedge6, DRT::Element::quad8,
    3>;
template class DRT::ELEMENTS::XFLUID::HybridLMCoupling<DRT::Element::wedge6, DRT::Element::quad9,
    3>;
// template class DRT::ELEMENTS::XFLUID::HybridLMCoupling<DRT::Element::wedge15,
// DRT::Element::tri3,3>; template class
// DRT::ELEMENTS::XFLUID::HybridLMCoupling<DRT::Element::wedge15,  DRT::Element::tri6,3>;
template class DRT::ELEMENTS::XFLUID::HybridLMCoupling<DRT::Element::wedge15, DRT::Element::quad4,
    3>;
template class DRT::ELEMENTS::XFLUID::HybridLMCoupling<DRT::Element::wedge15, DRT::Element::quad8,
    3>;
template class DRT::ELEMENTS::XFLUID::HybridLMCoupling<DRT::Element::wedge15, DRT::Element::quad9,
    3>;

template class DRT::ELEMENTS::XFLUID::HybridLMCoupling<DRT::Element::hex8, DRT::Element::dis_none,
    3>;
template class DRT::ELEMENTS::XFLUID::HybridLMCoupling<DRT::Element::hex20, DRT::Element::dis_none,
    3>;
template class DRT::ELEMENTS::XFLUID::HybridLMCoupling<DRT::Element::hex27, DRT::Element::dis_none,
    3>;
template class DRT::ELEMENTS::XFLUID::HybridLMCoupling<DRT::Element::tet4, DRT::Element::dis_none,
    3>;
template class DRT::ELEMENTS::XFLUID::HybridLMCoupling<DRT::Element::tet10, DRT::Element::dis_none,
    3>;
template class DRT::ELEMENTS::XFLUID::HybridLMCoupling<DRT::Element::wedge6, DRT::Element::dis_none,
    3>;
template class DRT::ELEMENTS::XFLUID::HybridLMCoupling<DRT::Element::wedge15,
    DRT::Element::dis_none, 3>;


// pairs with numdof=4
// template class DRT::ELEMENTS::XFLUID::HybridLMCoupling<DRT::Element::hex8, DRT::Element::tri3,4>;
// template class DRT::ELEMENTS::XFLUID::HybridLMCoupling<DRT::Element::hex8, DRT::Element::tri6,4>;
template class DRT::ELEMENTS::XFLUID::HybridLMCoupling<DRT::Element::hex8, DRT::Element::quad4, 4>;
template class DRT::ELEMENTS::XFLUID::HybridLMCoupling<DRT::Element::hex8, DRT::Element::quad8, 4>;
template class DRT::ELEMENTS::XFLUID::HybridLMCoupling<DRT::Element::hex8, DRT::Element::quad9, 4>;
// template class DRT::ELEMENTS::XFLUID::HybridLMCoupling<DRT::Element::hex20,
// DRT::Element::tri3,4>; template class
// DRT::ELEMENTS::XFLUID::HybridLMCoupling<DRT::Element::hex20, DRT::Element::tri6,4>;
template class DRT::ELEMENTS::XFLUID::HybridLMCoupling<DRT::Element::hex20, DRT::Element::quad4, 4>;
template class DRT::ELEMENTS::XFLUID::HybridLMCoupling<DRT::Element::hex20, DRT::Element::quad8, 4>;
template class DRT::ELEMENTS::XFLUID::HybridLMCoupling<DRT::Element::hex20, DRT::Element::quad9, 4>;
// template class DRT::ELEMENTS::XFLUID::HybridLMCoupling<DRT::Element::hex27,
// DRT::Element::tri3,4>; template class
// DRT::ELEMENTS::XFLUID::HybridLMCoupling<DRT::Element::hex27, DRT::Element::tri6,4>;
template class DRT::ELEMENTS::XFLUID::HybridLMCoupling<DRT::Element::hex27, DRT::Element::quad4, 4>;
template class DRT::ELEMENTS::XFLUID::HybridLMCoupling<DRT::Element::hex27, DRT::Element::quad8, 4>;
template class DRT::ELEMENTS::XFLUID::HybridLMCoupling<DRT::Element::hex27, DRT::Element::quad9, 4>;
// template class DRT::ELEMENTS::XFLUID::HybridLMCoupling<DRT::Element::tet4, DRT::Element::tri3,4>;
// template class DRT::ELEMENTS::XFLUID::HybridLMCoupling<DRT::Element::tet4, DRT::Element::tri6,4>;
template class DRT::ELEMENTS::XFLUID::HybridLMCoupling<DRT::Element::tet4, DRT::Element::quad4, 4>;
template class DRT::ELEMENTS::XFLUID::HybridLMCoupling<DRT::Element::tet4, DRT::Element::quad8, 4>;
template class DRT::ELEMENTS::XFLUID::HybridLMCoupling<DRT::Element::tet4, DRT::Element::quad9, 4>;
// template class DRT::ELEMENTS::XFLUID::HybridLMCoupling<DRT::Element::tet10,
// DRT::Element::tri3,4>; template class
// DRT::ELEMENTS::XFLUID::HybridLMCoupling<DRT::Element::tet10, DRT::Element::tri6,4>;
template class DRT::ELEMENTS::XFLUID::HybridLMCoupling<DRT::Element::tet10, DRT::Element::quad4, 4>;
template class DRT::ELEMENTS::XFLUID::HybridLMCoupling<DRT::Element::tet10, DRT::Element::quad8, 4>;
template class DRT::ELEMENTS::XFLUID::HybridLMCoupling<DRT::Element::tet10, DRT::Element::quad9, 4>;
// template class DRT::ELEMENTS::XFLUID::HybridLMCoupling<DRT::Element::wedge6,
// DRT::Element::tri3,4>; template class
// DRT::ELEMENTS::XFLUID::HybridLMCoupling<DRT::Element::wedge6, DRT::Element::tri6,4>;
template class DRT::ELEMENTS::XFLUID::HybridLMCoupling<DRT::Element::wedge6, DRT::Element::quad4,
    4>;
template class DRT::ELEMENTS::XFLUID::HybridLMCoupling<DRT::Element::wedge6, DRT::Element::quad8,
    4>;
template class DRT::ELEMENTS::XFLUID::HybridLMCoupling<DRT::Element::wedge6, DRT::Element::quad9,
    4>;
// template class DRT::ELEMENTS::XFLUID::HybridLMCoupling<DRT::Element::wedge15,
// DRT::Element::tri3,4>; template class
// DRT::ELEMENTS::XFLUID::HybridLMCoupling<DRT::Element::wedge15, DRT::Element::tri6,4>;
template class DRT::ELEMENTS::XFLUID::HybridLMCoupling<DRT::Element::wedge15, DRT::Element::quad4,
    4>;
template class DRT::ELEMENTS::XFLUID::HybridLMCoupling<DRT::Element::wedge15, DRT::Element::quad8,
    4>;
template class DRT::ELEMENTS::XFLUID::HybridLMCoupling<DRT::Element::wedge15, DRT::Element::quad9,
    4>;

// pairs with numdof=3
// template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::hex8,
// DRT::Element::tri3,3>; template class
// DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::hex8,  DRT::Element::tri6,3>;
template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::hex8,
    DRT::Element::quad4, 3>;
template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::hex8,
    DRT::Element::quad8, 3>;
template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::hex8,
    DRT::Element::quad9, 3>;
// template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::hex20,
// DRT::Element::tri3,3>; template class
// DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::hex20,  DRT::Element::tri6,3>;
template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::hex20,
    DRT::Element::quad4, 3>;
template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::hex20,
    DRT::Element::quad8, 3>;
template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::hex20,
    DRT::Element::quad9, 3>;
// template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::hex27,
// DRT::Element::tri3,3>; template class
// DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::hex27,  DRT::Element::tri6,3>;
template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::hex27,
    DRT::Element::quad4, 3>;
template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::hex27,
    DRT::Element::quad8, 3>;
template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::hex27,
    DRT::Element::quad9, 3>;
// template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::tet4,
// DRT::Element::tri3,3>; template class
// DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::tet4,  DRT::Element::tri6,3>;
template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::tet4,
    DRT::Element::quad4, 3>;
template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::tet4,
    DRT::Element::quad8, 3>;
template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::tet4,
    DRT::Element::quad9, 3>;
// template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::tet10,
// DRT::Element::tri3,3>; template class
// DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::tet10,  DRT::Element::tri6,3>;
template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::tet10,
    DRT::Element::quad4, 3>;
template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::tet10,
    DRT::Element::quad8, 3>;
template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::tet10,
    DRT::Element::quad9, 3>;
// template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::wedge6,
// DRT::Element::tri3,3>; template class
// DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::wedge6,  DRT::Element::tri6,3>;
template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::wedge6,
    DRT::Element::quad4, 3>;
template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::wedge6,
    DRT::Element::quad8, 3>;
template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::wedge6,
    DRT::Element::quad9, 3>;
// template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::wedge15,
// DRT::Element::tri3,3>; template class
// DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::wedge15,  DRT::Element::tri6,3>;
template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::wedge15,
    DRT::Element::quad4, 3>;
template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::wedge15,
    DRT::Element::quad8, 3>;
template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::wedge15,
    DRT::Element::quad9, 3>;

// volume coupled with numdof = 3, FSI Slavesided, FPI
template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::hex8,
    DRT::Element::hex8, 3>;
template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::hex20,
    DRT::Element::hex8, 3>;
template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::hex27,
    DRT::Element::hex8, 3>;
template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::tet4,
    DRT::Element::hex8, 3>;
template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::tet10,
    DRT::Element::hex8, 3>;
template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::wedge6,
    DRT::Element::hex8, 3>;
template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::wedge15,
    DRT::Element::hex8, 3>;

// pairs with numdof=4
// template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::hex8,
// DRT::Element::tri3,4>; template class
// DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::hex8,  DRT::Element::tri6,4>;
template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::hex8,
    DRT::Element::quad4, 4>;
template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::hex8,
    DRT::Element::quad8, 4>;
template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::hex8,
    DRT::Element::quad9, 4>;
// template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::hex20,
// DRT::Element::tri3,4>; template class
// DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::hex20, DRT::Element::tri6,4>;
template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::hex20,
    DRT::Element::quad4, 4>;
template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::hex20,
    DRT::Element::quad8, 4>;
template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::hex20,
    DRT::Element::quad9, 4>;
// template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::hex27,
// DRT::Element::tri3,4>; template class
// DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::hex27, DRT::Element::tri6,4>;
template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::hex27,
    DRT::Element::quad4, 4>;
template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::hex27,
    DRT::Element::quad8, 4>;
template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::hex27,
    DRT::Element::quad9, 4>;
// template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::tet4,
// DRT::Element::tri3,4>; template class
// DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::tet4,  DRT::Element::tri6,4>;
template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::tet4,
    DRT::Element::quad4, 4>;
template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::tet4,
    DRT::Element::quad8, 4>;
template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::tet4,
    DRT::Element::quad9, 4>;
// template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::tet10,
// DRT::Element::tri3,4>; template class
// DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::tet10, DRT::Element::tri6,4>;
template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::tet10,
    DRT::Element::quad4, 4>;
template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::tet10,
    DRT::Element::quad8, 4>;
template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::tet10,
    DRT::Element::quad9, 4>;
// template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::wedge6,
// DRT::Element::tri3,4>; template class
// DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::wedge6, DRT::Element::tri6,4>;
template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::wedge6,
    DRT::Element::quad4, 4>;
template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::wedge6,
    DRT::Element::quad8, 4>;
template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::wedge6,
    DRT::Element::quad9, 4>;
// template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::wedge15,
// DRT::Element::tri3,4>; template class
// DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::wedge15, DRT::Element::tri6,4>;
template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::wedge15,
    DRT::Element::quad4, 4>;
template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::wedge15,
    DRT::Element::quad8, 4>;
template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::wedge15,
    DRT::Element::quad9, 4>;

// template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::hex8,
// DRT::Element::tet4, 4>; template class
// DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::hex8,  DRT::Element::tet10,4>;
template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::hex8,
    DRT::Element::hex8, 4>;
template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::hex8,
    DRT::Element::hex20, 4>;
template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::hex8,
    DRT::Element::hex27, 4>;
// template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::hex20,
// DRT::Element::tet4,4>; template class
// DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::hex20, DRT::Element::tet10,4>;
template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::hex20,
    DRT::Element::hex8, 4>;
template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::hex20,
    DRT::Element::hex20, 4>;
template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::hex20,
    DRT::Element::hex27, 4>;
// template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::hex27,
// DRT::Element::tet4,4>; template class
// DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::hex27, DRT::Element::tet10,4>;
template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::hex27,
    DRT::Element::hex8, 4>;
template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::hex27,
    DRT::Element::hex20, 4>;
template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::hex27,
    DRT::Element::hex27, 4>;
// template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::tet4,
// DRT::Element::tet4,4>; template class
// DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::tet4,  DRT::Element::tet10,4>;
template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::tet4,
    DRT::Element::hex8, 4>;
template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::tet4,
    DRT::Element::hex20, 4>;
template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::tet4,
    DRT::Element::hex27, 4>;
// template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::tet10,
// DRT::Element::tet4,4>; template class
// DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::tet10, DRT::Element::tet10,4>;
template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::tet10,
    DRT::Element::hex8, 4>;
template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::tet10,
    DRT::Element::hex20, 4>;
template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::tet10,
    DRT::Element::hex27, 4>;
// template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::wedge6,
// DRT::Element::tet4,4>; template class
// DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::wedge6, DRT::Element::tet10,4>;
template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::wedge6,
    DRT::Element::hex8, 4>;
template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::wedge6,
    DRT::Element::hex20, 4>;
template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::wedge6,
    DRT::Element::hex27, 4>;
// template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::wedge15,
// DRT::Element::tet4,4>; template class
// DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::wedge15, DRT::Element::tet10,4>;
template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::wedge15,
    DRT::Element::hex8, 4>;
template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::wedge15,
    DRT::Element::hex20, 4>;
template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::wedge15,
    DRT::Element::hex27, 4>;
