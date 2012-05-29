/*----------------------------------------------------------------------*/
/*!
\file fluid_ele_calc_intfaces_stab.H

\brief edge-oriented/continuous interior penalty stabilization for fluid (especially xfluid) problems.

Literature:

    Edge stabilization for the incompressible Navier-Stokes equations: a continuous interior penalty finite element method
    E.Burman, M.A.Fernandez and P.Hansbo (2006)


    Finite element methods with symmetric stabilization for the transient convection-diffusion-reaction equation
    E.Burman, M.A.Fernandez
    Comput. Methods Appl. Mech. Engrg. 198 (2009) 2508-2519


<pre>
Maintainer: Benedikt Schott
            schott@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>
*/
/*----------------------------------------------------------------------*/

#ifndef FLUID_ELE_CALC_INTFACES_STAB_H
#define FLUID_ELE_CALC_INTFACES_STAB_H

#include "fluid_ele.H"
#include "fluid_ele_parameter.H"

#include "../drt_lib/drt_utils.H"

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

    \brief an interface class for edge-oriented/continuous interior penalty stabilization
           for every surface element distype this class allocates one
           instance of the edge-oriented stabilization implementation

    \author schott (schott@lnm.mw.tum.de), 02/12

    -----------------------------------------------------------------*/
    class FluidIntFaceStab : DRT::SingletonDestruction
    {
    public:

      //! Empty constructor
      FluidIntFaceStab() {}
      //! Empty destructor
      virtual ~FluidIntFaceStab() {}

      /*!
        \brief Evaluate the surface elements edge-oriented stabilization and ghost penalty

        This class does not provide a definition for this function, it
        is defined in the implementation.

        \author schott (schott@lnm.mw.tum.de), 02/12

      */
      virtual int EvaluateEdgeBasedStabilization(
          DRT::ELEMENTS::FluidIntFace*       intface,              ///< internal face element
          Teuchos::ParameterList&            params,               ///< parameter list
          DRT::Discretization&               discretization,       ///< discretization
          vector<int>&                       patchlm,              ///< patch local map
          vector<int>&                       lm_masterToPatch,     ///< local map between master dofs and patchlm
          vector<int>&                       lm_slaveToPatch,      ///< local map between slave dofs and patchlm
          vector<int>&                       lm_faceToPatch,       ///< local map between face dofs and patchlm
          std::vector<int>&                  lm_masterNodeToPatch, ///< local map between master nodes and nodes in patch
          std::vector<int>&                  lm_slaveNodeToPatch,  ///< local map between slave nodes and nodes in patch
          vector<Epetra_SerialDenseMatrix>&  elemat_blocks,        ///< element matrix blocks
          vector<Epetra_SerialDenseVector>&  elevec_blocks         ///< element vector blocks
        ) = 0;


      /*!
        \brief Allocate one static instance of the internal
               implementation class for edge oriented stabilization and
               return pointer to it

        \author schott (schott@lnm.mw.tum.de), 02/12

        \param surfele (in):   fluid internal surface element

      */
      static FluidIntFaceStab* Impl(DRT::ELEMENTS::FluidIntFace* surfele);

    };


    //-----------------------------------------------------------------
    //-----------------------------------------------------------------
    //
    //                        IMPLEMENTATION
    //
    //-----------------------------------------------------------------
    //-----------------------------------------------------------------

    /*!
    \brief Internal FluidInternalSurfaceStab EOS-implementation
           (surface and parent/neighbor element specific)
    */

    //-----------------------------------------------------------------
    ///
    //-----------------------------------------------------------------
    template<DRT::Element::DiscretizationType distype ,
             DRT::Element::DiscretizationType pdistype,
             DRT::Element::DiscretizationType ndistype>
    class FluidInternalSurfaceStab: public FluidIntFaceStab
    {
    public:

      /// Singleton access method
      static FluidInternalSurfaceStab<distype,pdistype,ndistype> * Instance(bool create=true);

      /// called upon destruction
      virtual void Done();

      /// Constructor with number of nodes
      FluidInternalSurfaceStab();

      /// empty destructor
      ~FluidInternalSurfaceStab(){return;};

      /// number of space dimensions of the FluidIntFace element
      static const int facensd_ = DRT::UTILS::DisTypeToDim<distype>::dim;

      /// number of space dimensions of the parent element
      static const int nsd_ = facensd_+1;

      /// number of dof's per node
      static const int numdofpernode_ = nsd_ + 1;

      /// number of nodes
      static const int iel  = DRT::UTILS::DisTypeToNumNodePerEle<distype >::numNodePerElement;

      /// number of parentnodes
      static const int piel = DRT::UTILS::DisTypeToNumNodePerEle<pdistype>::numNodePerElement;

      /// number of parentnodes of neighbor element
      static const int niel = DRT::UTILS::DisTypeToNumNodePerEle<ndistype>::numNodePerElement;



      /*!
        \brief Evaluate EOS stabilization

        This method calculates the contributions to rhs and matrix of
        EOS and ghost penalty stabilization for a generalized alpha
        system.

        Literature:
        Burman Fernandez and Hansbo (2006-2009)

      */
      virtual int EvaluateEdgeBasedStabilization(
          DRT::ELEMENTS::FluidIntFace*       intface,              ///< internal face element
          Teuchos::ParameterList&            params,               ///< parameter list
          DRT::Discretization&               discretization,       ///< discretization
          vector<int>&                       patchlm,              ///< patch local map
          vector<int>&                       lm_masterToPatch,     ///< local map between master dofs and patchlm
          vector<int>&                       lm_slaveToPatch,      ///< local map between slave dofs and patchlm
          vector<int>&                       lm_faceToPatch,       ///< local map between face dofs and patchlm
          std::vector<int>&                  lm_masterNodeToPatch, ///< local map between master nodes and nodes in patch
          std::vector<int>&                  lm_slaveNodeToPatch,  ///< local map between slave nodes and nodes in patch
          vector<Epetra_SerialDenseMatrix>&  elemat_blocks,        ///< element matrix blocks
          vector<Epetra_SerialDenseVector>&  elevec_blocks         ///< element vector blocks
      );

    private:

      //! reassemble matrix block from master-slave pairs to patch-node block for field (row, col)
      void ReassembleMATBlock(
          int                                               row_block,     ///< row block
          int                                               col_block,     ///< column block
          Epetra_SerialDenseMatrix&                         mat_block,     ///< matrix block
          LINALG::Matrix<4*piel, 4*piel>&                   elematrix_mm,  ///< element matrix master-master block
          LINALG::Matrix<4*piel, 4*niel>&                   elematrix_ms,  ///< element matrix master-slave block
          LINALG::Matrix<4*niel, 4*piel>&                   elematrix_sm,  ///< element matrix slave-master block
          LINALG::Matrix<4*niel, 4*niel>&                   elematrix_ss,  ///< element matrix slave-slave block
          std::vector<int>&                                 lm_masterNodeToPatch, ///< local map between master nodes and nodes in patch
          std::vector<int>&                                 lm_slaveNodeToPatch   ///< local map between slave nodes and nodes in patch
      );

      //! reassemble rhs block from master/slave rhs to patch-node block for field (row)
      void ReassembleRHSBlock(
          int                            row_block,            ///< row block
          Epetra_SerialDenseVector&      rhs_block,            ///< rhs block
          LINALG::Matrix<4*piel, 1>&     elevector_m,          ///< element vector master block
          LINALG::Matrix<4*niel, 1>&     elevector_s,          ///< element vector slave block
          std::vector<int>&              lm_masterNodeToPatch, ///< local map between master nodes and nodes in patch
          std::vector<int>&              lm_slaveNodeToPatch   ///< local map between slave nodes and nodes in patch
      );

      //! fill element vectors with extracted data
      void GetElementData(
          FluidIntFace*              surfele,          ///< surface FluidIntFace element
          Fluid*                    master_ele,       ///< master parent element
          Fluid*                    slave_ele,        ///< slave  parent element
          double &                   kinvisc,          ///< patch kinematic viscosity
          double &                   dens,             ///< patch density
          vector<double>&            mypvelaf,         ///< master velaf
          vector<double>&            mypvelnp,         ///< master velnp
          vector<double>&            mypedispnp,       ///< master dispnp
          vector<double>&            myedispnp,        ///< surfele dispnp
          vector<double>&            mynvelaf,         ///< slave velaf
          vector<double>&            mynvelnp,         ///< slave velnp
          vector<double>&            mynedispnp        ///< slave dispnp
      );

      //! evaluate shape functions and derivatives at integr. point
      double EvalShapeFuncAndDerivsAtIntPoint(
          DRT::UTILS::IntegrationPoints2D&  intpoints,     ///< reference to 2D integration points
          int                               iquad,         ///< actual integration point
          int                               master_eid,    ///< master parent element
          int                               slave_eid,     ///< slave parent element
          bool                              use2ndderiv = false
      );

      //! Provide ghost penalty assembly for xfluid case
      void GhostPenalty(
          LINALG::Matrix<4*piel, 4*piel>&                   elematrix_mm,  ///< element matrix master-master block
          LINALG::Matrix<4*piel, 4*niel>&                   elematrix_ms,  ///< element matrix master-slave block
          LINALG::Matrix<4*niel, 4*piel>&                   elematrix_sm,  ///< element matrix slave-master block
          LINALG::Matrix<4*niel, 4*niel>&                   elematrix_ss,  ///< element matrix slave-slave block
          LINALG::Matrix<4*piel, 1>&                        elevector_m,   ///< element vector master block
          LINALG::Matrix<4*niel, 1>&                        elevector_s,   ///< element vector slave block
          const double &                                    timefacfac,    ///< (time factor) x (integration factor)
          double &                                          tau_grad,      ///< penalty parameter for ghost penalty terms
          bool &                                            ghost_penalty, ///< boolian if ghost penalty terms has to be assembled
          bool &                                            use2ndderiv
      );

      //! Provide pressure (EOS) stabilization assembly for fluid
      void pressureEOS(
          LINALG::Matrix<4*piel, 4*piel>&                   elematrix_mm,  ///< element matrix master-master block
          LINALG::Matrix<4*piel, 4*niel>&                   elematrix_ms,  ///< element matrix master-slave block
          LINALG::Matrix<4*niel, 4*piel>&                   elematrix_sm,  ///< element matrix slave-master block
          LINALG::Matrix<4*niel, 4*niel>&                   elematrix_ss,  ///< element matrix slave-slave block
          LINALG::Matrix<4*piel, 1>&                        elevector_m,   ///< element vector master block
          LINALG::Matrix<4*niel, 1>&                        elevector_s,   ///< element vector slave block
          const double &                                    timefacfacpre, ///< (time factor pressure) x (integration factor)
          double &                                          tau_p          ///< penalty parameter for pressure stabilization terms
      );

      //! Provide divergence and streamline (EOS) stabilization assembly for fluid
      void div_streamline_EOS(
          LINALG::Matrix<4*piel, 4*piel>&                   elematrix_mm,  ///< element matrix master-master block
          LINALG::Matrix<4*piel, 4*niel>&                   elematrix_ms,  ///< element matrix master-slave block
          LINALG::Matrix<4*niel, 4*piel>&                   elematrix_sm,  ///< element matrix slave-master block
          LINALG::Matrix<4*niel, 4*niel>&                   elematrix_ss,  ///< element matrix slave-slave block
          LINALG::Matrix<4*piel, 1>&                        elevector_m,   ///< element vector master block
          LINALG::Matrix<4*niel, 1>&                        elevector_s,   ///< element vector slave block
          const double &                                    timefacfac,          ///< (time factor ) x (integration factor)
          double &                                          tau_div_streamline,  ///< combined penalty parameter for divergence and streamline stabilization terms
          LINALG::Matrix<3,3>&                              vderxyaf_diff        ///< velocity derivatives (neighbor-parent) element
      );


      //! Provide pressure (EOS) stabilization assembly for fluid (gradients in normal direction)
      void pressureEOSnormal(
          LINALG::Matrix<4*piel, 4*piel>&                   elematrix_mm,  ///< element matrix master-master block
          LINALG::Matrix<4*piel, 4*niel>&                   elematrix_ms,  ///< element matrix master-slave block
          LINALG::Matrix<4*niel, 4*piel>&                   elematrix_sm,  ///< element matrix slave-master block
          LINALG::Matrix<4*niel, 4*niel>&                   elematrix_ss,  ///< element matrix slave-slave block
          LINALG::Matrix<4*piel, 1>&                        elevector_m,   ///< element vector master block
          LINALG::Matrix<4*niel, 1>&                        elevector_s,   ///< element vector slave block
          const double &                                    timefacfacpre, ///< (time factor pressure) x (integration factor)
          double &                                          tau_p          ///< penalty parameter for pressure stabilization terms
      );


      //! Provide divergence (EOS) stabilization assembly for fluid
      void div_EOS(
          LINALG::Matrix<4*piel, 4*piel>&                   elematrix_mm,  ///< element matrix master-master block
          LINALG::Matrix<4*piel, 4*niel>&                   elematrix_ms,  ///< element matrix master-slave block
          LINALG::Matrix<4*niel, 4*piel>&                   elematrix_sm,  ///< element matrix slave-master block
          LINALG::Matrix<4*niel, 4*niel>&                   elematrix_ss,  ///< element matrix slave-slave block
          LINALG::Matrix<4*piel, 1>&                        elevector_m,   ///< element vector master block
          LINALG::Matrix<4*niel, 1>&                        elevector_s,   ///< element vector slave block
          const double &                                    timefacfac,    ///< (time factor ) x (integration factor)
          double &                                          tau_div,       ///< penalty parameter for divergence stabilization terms
          LINALG::Matrix<3,3>&                              vderxyaf_diff  ///< velocity derivatives (neighbor-parent) element
      );


      //! Provide streamline (EOS) stabilization assembly for fluid
      void streamline_EOS(
          LINALG::Matrix<4*piel, 4*piel>&                   elematrix_mm,  ///< element matrix master-master block
          LINALG::Matrix<4*piel, 4*niel>&                   elematrix_ms,  ///< element matrix master-slave block
          LINALG::Matrix<4*niel, 4*piel>&                   elematrix_sm,  ///< element matrix slave-master block
          LINALG::Matrix<4*niel, 4*niel>&                   elematrix_ss,  ///< element matrix slave-slave block
          LINALG::Matrix<4*piel, 1>&                        elevector_m,   ///< element vector master block
          LINALG::Matrix<4*niel, 1>&                        elevector_s,   ///< element vector slave block
          const double &                                    timefacfac,    ///< (time factor ) x (integration factor)
          double &                                          tau_streamline,///< combined penalty parameter for streamline stabilization terms
          LINALG::Matrix<3,3>&                              vderxyaf_diff, ///< velocity derivatives (neighbor-parent) element
          LINALG::Matrix<3,1>&                              conv_diff      ///< convective velocity (neighbor-parent) element
      );


      //! compute h_k w.r.t master and slave element
      void compute_patch_hk ( double &   patch_hk,   ///< patch h_k
                              Fluid*    master,     ///< master fluid element
                              Fluid*    slave       ///< slave fluid element
      );

      //! compute surface diameter w.r.t to parent master element or parent slave element (flag master=true/false)
      template <int numnode>
      void diameter(     bool                master,           ///< master or slave parent element
                         std::vector<int> &  connectivity,     ///< connectivity vector for parent element
                         double &            h_e               ///< element length w.r.t parent element
      )
      {
        LINALG::Matrix<3,numnode> xyz_surf(true);

        if(connectivity.size() != numnode) dserror("wrong number of nodes for parent's surface");

        if(master == true) // is parent element the master element?
        {
          // extract node coords
          for(int i=0;i<(int)numnode;++i)
          {
            int col = connectivity[i];
            xyz_surf(0,i)=pxyze_(0, col);
            xyz_surf(1,i)=pxyze_(1, col);
            xyz_surf(2,i)=pxyze_(2, col);
          }
        }
        else // parent element is the slave element
        {
          // extract node coords
          for(int i=0;i<(int)numnode;++i)
          {
            int col = connectivity[i];
            xyz_surf(0,i)=nxyze_(0, col);
            xyz_surf(1,i)=nxyze_(1, col);
            xyz_surf(2,i)=nxyze_(2, col);
          }
        }


        if(numnode == 4) // quad4 surface
        {
          double diam1 = 0.0;
          double diam2 = 0.0;
          double line0 = 0.0;
          double line1 = 0.0;
          double line2 = 0.0;
          double line3 = 0.0;

          for(int i=0; i<3; i++)
          {

            // diagonals
            // line nodes (0,2)
            diam1 += (xyz_surf(i,2) - xyz_surf(i,0))*(xyz_surf(i,2) - xyz_surf(i,0));

            // line nodes (1,3)
            diam2 += (xyz_surf(i,3) - xyz_surf(i,1))*(xyz_surf(i,3) - xyz_surf(i,1));

            // lines
            line0 += (xyz_surf(i,1) - xyz_surf(i,0))*(xyz_surf(i,1) - xyz_surf(i,0));
            line1 += (xyz_surf(i,2) - xyz_surf(i,1))*(xyz_surf(i,2) - xyz_surf(i,1));
            line2 += (xyz_surf(i,3) - xyz_surf(i,2))*(xyz_surf(i,3) - xyz_surf(i,2));
            line3 += (xyz_surf(i,0) - xyz_surf(i,3))*(xyz_surf(i,0) - xyz_surf(i,3));
          }
          diam1 = sqrt(diam1);
          diam2 = sqrt(diam2);
          line0 = sqrt(line0);
          line1 = sqrt(line1);
          line2 = sqrt(line2);
          line3 = sqrt(line3);


          h_e = max(diam1, max(diam2, max(line0, max(line1, max(line2, line3)))));
        }
        else if(numnode == 3) // tri3 surface
        {
          double line0 = 0.0;
          double line1 = 0.0;
          double line2 = 0.0;

          for(int i=0; i<3; i++)
          {
            // no diagonals

            // lines
            line0 += (xyz_surf(i,1) - xyz_surf(i,0))*(xyz_surf(i,1) - xyz_surf(i,0));
            line1 += (xyz_surf(i,2) - xyz_surf(i,1))*(xyz_surf(i,2) - xyz_surf(i,1));
            line2 += (xyz_surf(i,0) - xyz_surf(i,3))*(xyz_surf(i,0) - xyz_surf(i,3));
          }

          line0 = sqrt(line0);
          line1 = sqrt(line1);
          line2 = sqrt(line2);


          h_e = max(line0, max(line1, line2));
        }
        else dserror("unknown number of nodes");

        if(h_e <= 0.0) dserror("negative or zero diameter for current face!");

        return;
       }



      //! pointer to parameter list
      DRT::ELEMENTS::FluidEleParameter* fldpara_;

      // nodal arrays
      // ------------
      //! node coordinates of parent (master) element
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

      //! node coordinates of neighbor (slave) element
      LINALG::Matrix<3,niel> nxyze_;
      //! array of nodal velocities, intermediate time level
      LINALG::Matrix<3,niel> nevelaf_;
      //! array of nodal velocities, new time level
      LINALG::Matrix<3,niel> nevelnp_;
      //! array of nodal pressure, new time level
      LINALG::Matrix<niel,1> neprenp_;
      //! array of nodal grid displacements, parent element, new time level
      LINALG::Matrix<3,niel> nedispnp_;

      //! node coordinates of boundary element
      LINALG::Matrix<3, iel> xyze_;

      //! linearisation of convection, convective part for parent element
      LINALG::Matrix<piel,1> p_conv_c;
      //! linearisation of convection, convective part for neighbor element
      LINALG::Matrix<niel,1> n_conv_c;


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
      //! vector of shape function (2nd) derivatives in reference coordinate system, parent element
      LINALG::Matrix<6,piel> pderiv2_;
      //! vector of shape function (2nd) derivatives in global coordinate system
      LINALG::Matrix<6,piel> pderxy2_;

      // ----------------------------------------------------------------------------------
      //! transpose of the jacobian matrix of the mapping (r,s,t)->(x,y,z)
      LINALG::Matrix<3,3>    nxjm_;
      //! its inverse
      LINALG::Matrix<3,3>    nxji_;
      //! vector of shape functions, parent element
      LINALG::Matrix<niel,1> nfunct_;
      //! vector of shape function derivatives in reference coordinate system, parent element
      LINALG::Matrix<3,niel> nderiv_;
      //! vector of shape function derivatives in global coordinate system
      LINALG::Matrix<3,niel> nderxy_;
      //! vector of shape function (2nd) derivatives in reference coordinate system, parent element
      LINALG::Matrix<6,niel> nderiv2_;
      //! vector of shape function (2nd) derivatives in global coordinate system
      LINALG::Matrix<6,niel> nderxy2_;


      //! vector of shape functions, boundary element
      LINALG::Matrix<iel ,1> funct_;
      //! vector of shape function derivatives in reference coordinate system, boundary element
      LINALG::Matrix<2, iel> deriv_;
      //! normal vector on surface element (outward pointing from parent element)
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
      //! parent velocity derivatives in gausspoint, time n+af
      LINALG::Matrix<3,3>    pvderxyaf_;
      //! neighbor velocity derivatives in gausspoint, time n+af
      LINALG::Matrix<3,3>    nvderxyaf_;
      //! parent velocity (2nd) derivatives in gausspoint, time n+af
      LINALG::Matrix<3,6>    pvderxy2af_;
      //! neighbor velocity (2nd) derivatives in gausspoint, time n+af
      LINALG::Matrix<3,6>    nvderxy2af_;
      //! pressure in gausspoint, time n+1
      double                 prenp_;
      //! parent pressure derivatives in gausspoint, time n+1
      LINALG::Matrix<3,1>    pprederxy_;
      //! neighbor velocity derivatives in gausspoint, time n+1
      LINALG::Matrix<3,1>    nprederxy_;



    };




    //-----------------------------------------------------------------
    //-----------------------------------------------------------------
    //
    //                    EDGE_BASED_STABILIZATION
    //
    //-----------------------------------------------------------------
    //-----------------------------------------------------------------

    /*!
      \brief ...

      Literature: ...

      and references therein

    */
    class FluidEdgeBasedStab
    {
    public:


      /*! \brief constructor
      */
      FluidEdgeBasedStab( double patch_hk );

      /// Empty destructor
      virtual ~FluidEdgeBasedStab() {}


      /*!
       \brief computation of the EOS-stabilization parameter and ghost penalty parameter
      */
      void ComputeStabilizationParams(
          double&       tau_grad,
          double&       tau_u,
          double&       tau_div,
          double&       tau_p,
          double&       tau_u_lin,
          double&       tau_div_lin,
          double&       tau_p_lin,
          double&       kinvisc,
          double&       density,
          double&       max_pevelnp,
          const double&       timefac,
          const double&       gamma_ghost_penalty);


    private:

      //! longest edge in parent element
      double p_hk_;



    };


  } // namespace ELEMENTS
} // namespace DRT

#endif
