/*----------------------------------------------------------------------*/
/*! \file

\brief Template classes for interface coupling in the XFEM with slave element representation

\level 2


*/
/*----------------------------------------------------------------------*/

#include <Teuchos_TimeMonitor.hpp>

#include "../drt_cut/cut_position.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_xfem/xfem_interface_utils.H"
#include "fluid_ele_calc_xfem_coupling_impl.H"

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


    }  // namespace XFLUID
  }    // namespace ELEMENTS
}  // namespace DRT


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
