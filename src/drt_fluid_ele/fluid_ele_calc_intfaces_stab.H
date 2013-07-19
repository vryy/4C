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
#include "../drt_fem_general/drt_utils_gausspoints.H"

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
          DRT::ELEMENTS::FluidIntFace*       intface,                   ///< internal face element
          DRT::ELEMENTS::FluidEleParameter&  fldpara_,
          Teuchos::ParameterList&            params,                    ///< parameter list
          DRT::Discretization&               discretization,            ///< discretization
          std::vector<int>&                  patchlm,                   ///< patch local map
          std::vector<int>&                  lm_masterToPatch,          ///< local map between master dofs and patchlm
          std::vector<int>&                  lm_slaveToPatch,           ///< local map between slave dofs and patchlm
          std::vector<int>&                  lm_faceToPatch,            ///< local map between face dofs and patchlm
          std::vector<int>&                  lm_masterNodeToPatch,      ///< local map between master nodes and nodes in patch
          std::vector<int>&                  lm_slaveNodeToPatch,       ///< local map between slave nodes and nodes in patch
          std::vector<Epetra_SerialDenseMatrix>&  elemat_blocks,        ///< element matrix blocks
          std::vector<Epetra_SerialDenseVector>&  elevec_blocks         ///< element vector blocks
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

      /// number of second order derivatives for master element
      static const int numderiv2_p = DRT::UTILS::DisTypeToNumDeriv2<pdistype>::numderiv2;

      /// number of second order derivatives for slave element
      static const int numderiv2_n = DRT::UTILS::DisTypeToNumDeriv2<ndistype>::numderiv2;



      /*!
        \brief Evaluate EOS stabilization

        This method calculates the contributions to rhs and matrix of
        EOS and ghost penalty stabilization for a generalized alpha
        system.

        Literature:
        Burman Fernandez and Hansbo (2006-2009)

      */
      virtual int EvaluateEdgeBasedStabilization(
          DRT::ELEMENTS::FluidIntFace*       intface,                   ///< internal face element
          DRT::ELEMENTS::FluidEleParameter&  fldpara,                   ///< fluid parameter
          Teuchos::ParameterList&            params,                    ///< parameter list
          DRT::Discretization&               discretization,            ///< discretization
          std::vector<int>&                  patchlm,                   ///< patch local map
          std::vector<int>&                  lm_masterToPatch,          ///< local map between master dofs and patchlm
          std::vector<int>&                  lm_slaveToPatch,           ///< local map between slave dofs and patchlm
          std::vector<int>&                  lm_faceToPatch,            ///< local map between face dofs and patchlm
          std::vector<int>&                  lm_masterNodeToPatch,      ///< local map between master nodes and nodes in patch
          std::vector<int>&                  lm_slaveNodeToPatch,       ///< local map between slave nodes and nodes in patch
          std::vector<Epetra_SerialDenseMatrix>&  elemat_blocks,        ///< element matrix blocks
          std::vector<Epetra_SerialDenseVector>&  elevec_blocks         ///< element vector blocks
      );

    private:

      //! reassemble matrix block from master-slave pairs to patch-node block for field (row, col)
      void ReassembleMATBlock(
          int                                                          row_block,            ///< row block
          int                                                          col_block,            ///< column block
          Epetra_SerialDenseMatrix&                                    mat_block,            ///< matrix block
          LINALG::Matrix<numdofpernode_*piel, numdofpernode_*piel>&    elematrix_mm,         ///< element matrix master-master block
          LINALG::Matrix<numdofpernode_*piel, numdofpernode_*niel>&    elematrix_ms,         ///< element matrix master-slave block
          LINALG::Matrix<numdofpernode_*niel, numdofpernode_*piel>&    elematrix_sm,         ///< element matrix slave-master block
          LINALG::Matrix<numdofpernode_*niel, numdofpernode_*niel>&    elematrix_ss,         ///< element matrix slave-slave block
          std::vector<int>&                                            lm_masterNodeToPatch, ///< local map between master nodes and nodes in patch
          std::vector<int>&                                            lm_slaveNodeToPatch   ///< local map between slave nodes and nodes in patch
      );

      //! reassemble rhs block from master/slave rhs to patch-node block for field (row)
      void ReassembleRHSBlock(
          int                                         row_block,            ///< row block
          Epetra_SerialDenseVector&                   rhs_block,            ///< rhs block
          LINALG::Matrix<numdofpernode_*piel, 1>&     elevector_m,          ///< element vector master block
          LINALG::Matrix<numdofpernode_*niel, 1>&     elevector_s,          ///< element vector slave block
          std::vector<int>&                           lm_masterNodeToPatch, ///< local map between master nodes and nodes in patch
          std::vector<int>&                           lm_slaveNodeToPatch   ///< local map between slave nodes and nodes in patch
      );

      //! fill element vectors with extracted data
      void GetElementData(
          FluidIntFace*              surfele,          ///< surface FluidIntFace element
          Fluid*                     master_ele,       ///< master parent element
          Fluid*                     slave_ele,        ///< slave  parent element
          double &                   kinvisc,          ///< patch kinematic viscosity
          double &                   dens,             ///< patch density
          std::vector<double>&       mypvelaf,         ///< master velaf
          std::vector<double>&       mypvelnp,         ///< master velnp
          std::vector<double>&       mypedispnp,       ///< master dispnp
          std::vector<double>&       myedispnp,        ///< surfele dispnp
          std::vector<double>&       mynvelaf,         ///< slave velaf
          std::vector<double>&       mynvelnp,         ///< slave velnp
          std::vector<double>&       mynedispnp        ///< slave dispnp
      );

      double EvalShapeFuncAndDerivsAtIntPoint(
          const double                             wquad,                ///< Gaussian weight
          const LINALG::Matrix<facensd_,1> &       xi_gp,                ///< local coordinates of gaussian point w.r.t the master's face
          const LINALG::Matrix<nsd_,1> &           p_xi_gp,              ///< local coordinates of gaussian point w.r.t master element
          const LINALG::Matrix<nsd_,1> &           n_xi_gp,              ///< local coordinates of gaussian point w.r.t slave element
          int                                      master_eid,           ///< master parent element
          int                                      slave_eid,            ///< slave parent element
          bool                                     use2ndderiv = false   ///< flag to use 2nd order derivatives
      );

      //! evaluate shape functions and derivatives at integr. point
      double EvalShapeFuncAndDerivsAtIntPoint(
          DRT::UTILS::GaussIntegration::iterator & iquad,                 ///< actual integration point
          int                                      master_eid,            ///< master parent element
          int                                      slave_eid,             ///< slave parent element
          bool                                     use2ndderiv = false    ///< flag to use 2nd order derivatives
      );

      //! Provid pressure and viscous u (EOS) ghost penalty stabilization for 2nd order derivatives
      void GhostPenalty2nd(
                  LINALG::Matrix<numdofpernode_*piel, numdofpernode_*piel>&     elematrix_mm,  ///< element matrix master-master block
                  LINALG::Matrix<numdofpernode_*piel, numdofpernode_*niel>&     elematrix_ms,  ///< element matrix master-slave block
                  LINALG::Matrix<numdofpernode_*niel, numdofpernode_*piel>&     elematrix_sm,  ///< element matrix slave-master block
                  LINALG::Matrix<numdofpernode_*niel, numdofpernode_*niel>&     elematrix_ss,  ///< element matrix slave-slave block
                  LINALG::Matrix<numdofpernode_*piel, 1>&                       elevector_m,   ///< element vector master block
                  LINALG::Matrix<numdofpernode_*niel, 1>&                       elevector_s,   ///< element vector slave block
                  const double &                                                timefacfac,
                  const double &                                                timefacfacpre,
                  double &                                                      tau_u_2nd,
                  double &                                                      tau_p_2nd
                  );

      //! Provide pressure (EOS) stabilization assembly for fluid
      void pressureEOS(
          LINALG::Matrix<numdofpernode_*piel, numdofpernode_*piel>&      elematrix_mm,  ///< element matrix master-master block
          LINALG::Matrix<numdofpernode_*piel, numdofpernode_*niel>&      elematrix_ms,  ///< element matrix master-slave block
          LINALG::Matrix<numdofpernode_*niel, numdofpernode_*piel>&      elematrix_sm,  ///< element matrix slave-master block
          LINALG::Matrix<numdofpernode_*niel, numdofpernode_*niel>&      elematrix_ss,  ///< element matrix slave-slave block
          LINALG::Matrix<numdofpernode_*piel, 1>&                        elevector_m,   ///< element vector master block
          LINALG::Matrix<numdofpernode_*niel, 1>&                        elevector_s,   ///< element vector slave block
          const double &                                                 timefacfacpre, ///< (time factor pressure) x (integration factor)
          const double &                                                 timefacfacrhs, ///< (time factor rhs) x (integration factor)
          double &                                                       tau_p          ///< pressure stabilization parameter
      );

      //! Provide divergence and streamline (EOS) stabilization assembly for fluid
      void div_streamline_EOS(
          LINALG::Matrix<numdofpernode_*piel, numdofpernode_*piel>&      elematrix_mm,        ///< element matrix master-master block
          LINALG::Matrix<numdofpernode_*piel, numdofpernode_*niel>&      elematrix_ms,        ///< element matrix master-slave block
          LINALG::Matrix<numdofpernode_*niel, numdofpernode_*piel>&      elematrix_sm,        ///< element matrix slave-master block
          LINALG::Matrix<numdofpernode_*niel, numdofpernode_*niel>&      elematrix_ss,        ///< element matrix slave-slave block
          LINALG::Matrix<numdofpernode_*piel, 1>&                        elevector_m,         ///< element vector master block
          LINALG::Matrix<numdofpernode_*niel, 1>&                        elevector_s,         ///< element vector slave block
          const double &                                                 timefacfac,          ///< timefac x integration factor
          const double &                                                 timefacfacrhs,       ///< (time factor rhs) x (integration factor)
          double &                                                       tau_div_streamline,  ///< streamline stabilization parameter
          LINALG::Matrix<nsd_,nsd_>&                                     vderxyaf_diff        ///< difference of velocity gradients
      );

      //! Provide divergence (EOS) stabilization assembly for fluid
      void div_EOS(
          LINALG::Matrix<numdofpernode_*piel, numdofpernode_*piel>&      elematrix_mm,  ///< element matrix master-master block
          LINALG::Matrix<numdofpernode_*piel, numdofpernode_*niel>&      elematrix_ms,  ///< element matrix master-slave block
          LINALG::Matrix<numdofpernode_*niel, numdofpernode_*piel>&      elematrix_sm,  ///< element matrix slave-master block
          LINALG::Matrix<numdofpernode_*niel, numdofpernode_*niel>&      elematrix_ss,  ///< element matrix slave-slave block
          LINALG::Matrix<numdofpernode_*piel, 1>&                        elevector_m,   ///< element vector master block
          LINALG::Matrix<numdofpernode_*niel, 1>&                        elevector_s,   ///< element vector slave block
          const double &                                                 timefacfacpre, ///< (time factor pressure) x (integration factor)
          const double &                                                 timefacfacrhs, ///< (time factor rhs) x (integration factor)
          double &                                                       tau_div        ///< combined penalty parameter for divergence stabilization terms
      );


      //! compute h_k w.r.t master and slave element
      void compute_patch_hk ( double &   patch_hk,   ///< patch h_k
                              Fluid*    master,     ///< master fluid element
                              Fluid*    slave,       ///< slave fluid element
                              DRT::ELEMENTS::FluidIntFace*       intface,
                              const INPAR::FLUID::EOS_ElementLength&    eos_element_length
      );

      //! compute surface diameter w.r.t to parent master element or parent slave element (flag master=true/false)
      template <int numnode>
      void diameter2D(     bool                master,           ///< master or slave parent element
                           std::vector<int> &  connectivity,     ///< connectivity vector for parent element
                           double &            h_e               ///< element length w.r.t parent element
      )
      {
        LINALG::Matrix<nsd_,numnode> xyz_surf(true);

        if(connectivity.size() != numnode) dserror("wrong number of nodes for parent's surface");

        if(master == true) // is parent element the master element?
        {
          // extract node coords
          for(int i=0;i<(int)numnode;++i)
          {
            int col = connectivity[i];

            for(int isd=0; isd<nsd_; isd++)
              xyz_surf(isd,i)=pxyze_(isd, col);

          }
        }
        else // parent element is the slave element
        {
          // extract node coords
          for(int i=0;i<(int)numnode;++i)
          {
            int col = connectivity[i];

            for(int isd=0; isd<nsd_; isd++)
              xyz_surf(isd,i)=nxyze_(isd, col);

          }
        }


        if(numnode == 4 or numnode == 8) // quad4 surface and take the corner points also for quad8 surfaces
        {
          double diam1 = 0.0;
          double diam2 = 0.0;
          double line0 = 0.0;
          double line1 = 0.0;
          double line2 = 0.0;
          double line3 = 0.0;

          for(int i=0; i<nsd_; i++)
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


          h_e = std::max(diam1, std::max(diam2, std::max(line0, std::max(line1, std::max(line2, line3)))));
        }
        else if(numnode == 3) // tri3 surface
        {
          double line0 = 0.0;
          double line1 = 0.0;
          double line2 = 0.0;

          for(int i=0; i<nsd_; i++)
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


          h_e = std::max(line0, std::max(line1, line2));
        }
        else dserror("unknown number of nodes");

        if(h_e <= 0.0) dserror("negative or zero diameter for current face!");

        return;
       }

      //! compute surface diameter w.r.t to parent master element or parent slave element (flag master=true/false)
      template <int numnode>
      void diameter1D(     bool                master,           ///< master or slave parent element
                           std::vector<int> &  connectivity,     ///< connectivity vector for parent element
                           double &            h_e               ///< element length w.r.t parent element
      )
      {
        LINALG::Matrix<nsd_,numnode> xyz_surf(true);

        if(connectivity.size() != numnode) dserror("wrong number of nodes for parent's surface");

        if(master == true) // is parent element the master element?
        {
          // extract node coords
          for(int i=0;i<(int)numnode;++i)
          {
            int col = connectivity[i];

            for(int isd=0; isd<nsd_; isd++)
              xyz_surf(isd,i)=pxyze_(isd, col);

          }
        }
        else // parent element is the slave element
        {
          // extract node coords
          for(int i=0;i<(int)numnode;++i)
          {
            int col = connectivity[i];

            for(int isd=0; isd<nsd_; isd++)
              xyz_surf(isd,i)=nxyze_(isd, col);

          }
        }


        if(numnode == 2) // line 2 face
        {
          double diam1 = 0.0;

          for(int i=0; i<nsd_; i++)
          {

            // diagonals
            // line nodes (0,1) -> (0,1)
            diam1 += (xyz_surf(i,1) - xyz_surf(i,0))*(xyz_surf(i,1) - xyz_surf(i,0));

          }
          diam1 = sqrt(diam1);
          h_e = diam1;
        }
        else if(numnode == 3) // line3 face
        {
          double diam1 = 0.0;

          for(int i=0; i<nsd_; i++)
          {

            // diagonals
            // line nodes (0,2,1) -> (0,1)
            diam1 += (xyz_surf(i,1) - xyz_surf(i,0))*(xyz_surf(i,1) - xyz_surf(i,0));

          }
          diam1 = sqrt(diam1);
          h_e = diam1;
        }
        else dserror("unknown number of nodes");

        if(h_e <= 0.0) dserror("negative or zero diameter for current face!");

        return;
       }

      template <int numnode>
      void distance2D( bool                              master,           ///< master or slave parent element
    		           std::map< int,std::vector<int> >  p_lines_nodes,
                       double &                          h_e               ///< element length w.r.t parent element
      )
      {
    	  LINALG::Matrix<nsd_ ,iel*2> xyze_distance_lines;

    	  if (master == true) // is parent element the master element?
    	  {
    		  //find the distance
    		  int count = 0;
    		  for(std::map<int,std::vector<int> >::iterator iter = p_lines_nodes.begin(); iter!= p_lines_nodes.end();
    				  iter++)
    		  {
    			  std::vector<int> nodesofline = iter->second;
    			  for(int isd=0; isd<nsd_; isd++)
    			  {
    				  xyze_distance_lines(isd,count) = pxyze_(isd, nodesofline.at(0));
    				  xyze_distance_lines(isd,count+1) = pxyze_(isd, nodesofline.at(1));
    			  }
    			  count=count+2;
    		  }
    	  }
    	  else
    	  {
    		  //find the distance
    		  int count = 0;
    		  for(std::map<int,std::vector<int> >::iterator iter = p_lines_nodes.begin(); iter!= p_lines_nodes.end();
    				  iter++)
    		  {
    			  std::vector<int> nodesofline = iter->second;
    			  for(int isd=0; isd<nsd_; isd++)
    			  {
    				  xyze_distance_lines(isd,count) = nxyze_(isd, nodesofline.at(0));
    				  xyze_distance_lines(isd,count+1) = nxyze_(isd, nodesofline.at(1));
    			  }
    			  count=count+2;
    		  }

       	  }
          if (numnode == 4 or numnode == 8)
          {
        	  double line0 = 0.0;
        	  double line1 = 0.0;
        	  double line2 = 0.0;
        	  double line3 = 0.0;

        	  for(int i=0; i<nsd_; i++)
        	  {
        		  line0 += (xyze_distance_lines(i,0) - xyze_distance_lines(i,1))*(xyze_distance_lines(i,0) - xyze_distance_lines(i,1));
        		  line1 += (xyze_distance_lines(i,2) - xyze_distance_lines(i,3))*(xyze_distance_lines(i,2) - xyze_distance_lines(i,3));
        		  line2 += (xyze_distance_lines(i,4) - xyze_distance_lines(i,5))*(xyze_distance_lines(i,4) - xyze_distance_lines(i,5));
        		  line3 += (xyze_distance_lines(i,6) - xyze_distance_lines(i,7))*(xyze_distance_lines(i,6) - xyze_distance_lines(i,7));
        	  }

        	  line0 = sqrt(line0);
        	  line1 = sqrt(line1);
        	  line2 = sqrt(line2);
        	  line3 = sqrt(line3);

        	  h_e = std::max(line0, std::max(line1, std::max(line2, line3)));
          }
          if (numnode == 3)
          {
        	  double line0 = 0.0;

        	  for(int i=0; i<nsd_; i++)
        	  {
        		  line0 += (xyze_distance_lines(i,0) - xyze_distance_lines(i,1))*(xyze_distance_lines(i,0) - xyze_distance_lines(i,1));
        	  }
        	  line0 = sqrt(line0);

        	  h_e = line0;
          }
      }

      template <int numnode>
      void distance1D( bool                              master,           ///< master or slave parent element
    		           std::map< int,std::vector<int> >  p_lines_nodes,
                       double &                          h_e               ///< element length w.r.t parent element
      )
      {
    	  dserror("The implementation of distance1D is not tested 2D elements!");
    	  LINALG::Matrix<nsd_ ,iel*2> xyze_distance_lines;

    	  if (master == true) // is parent element the master element?
    	  {
    		  //find the distance
    		  int count = 0;
    		  for(std::map<int,std::vector<int> >::iterator iter = p_lines_nodes.begin(); iter!= p_lines_nodes.end();
    				  iter++)
    		  {
    			  std::vector<int> nodesofline = iter->second;
    			  for(int isd=0; isd<nsd_; isd++)
    			  {
    				  xyze_distance_lines(isd,count) = pxyze_(isd, nodesofline.at(0));
    				  xyze_distance_lines(isd,count+1) = pxyze_(isd, nodesofline.at(1));
    			  }
    			  count=count+2;
    		  }
    	  }
    	  else
    	  {
    		  //find the distance
    		  int count = 0;
    		  for(std::map<int,std::vector<int> >::iterator iter = p_lines_nodes.begin(); iter!= p_lines_nodes.end();
    				  iter++)
    		  {
    			  std::vector<int> nodesofline = iter->second;
    			  for(int isd=0; isd<nsd_; isd++)
    			  {
    				  xyze_distance_lines(isd,count) = nxyze_(isd, nodesofline.at(0));
    				  xyze_distance_lines(isd,count+1) = nxyze_(isd, nodesofline.at(1));
    			  }
    			  count=count+2;
    		  }

       	  }
          if (numnode == 4) // quad4 parent element
          {
        	  double line0 = 0.0;
        	  double line1 = 0.0;

        	  for(int i=0; i<nsd_; i++)
        	  {
        		  line0 += (xyze_distance_lines(i,0) - xyze_distance_lines(i,1))*(xyze_distance_lines(i,0) - xyze_distance_lines(i,1));
        		  line1 += (xyze_distance_lines(i,2) - xyze_distance_lines(i,3))*(xyze_distance_lines(i,2) - xyze_distance_lines(i,3));
        	  }

        	  line0 = sqrt(line0);
        	  line1 = sqrt(line1);
        	  h_e = std::max(line0, line1);
          }
          if (numnode == 3) //tri3 parent element
          {
        	  double line0 = 0.0;

        	  for(int i=0; i<nsd_; i++)
        	  {
        		  line0 += (xyze_distance_lines(i,0) - xyze_distance_lines(i,1))*(xyze_distance_lines(i,0) - xyze_distance_lines(i,1));
        	  }
        	  line0 = sqrt(line0);

        	  h_e = line0;
          }
      }
      //! find the surface at the opposite side
      template <int numnode>
      void FindOppositeSurface2D(int surfacelocalid,int & surf_neighbor)
      {
	  if ((numnode == 4) or (numnode == 8))
	  {
	    switch (surfacelocalid)
	    {
            case 0:
    	        surf_neighbor = 5;
    	    break;
            case 1:
             	surf_neighbor = 3;
    	    break;
            case 2:
    	        surf_neighbor = 4;
    	    break;
            case 3:
    	        surf_neighbor = 1;
    	    break;
            case 4:
    	        surf_neighbor = 2;
    	    break;
            case 5:
    	        surf_neighbor = 0;
    	    break;
            default:
    	        dserror("wrong number!!");
    	    break;
           }
	 }
	 else if (numnode == 3) 
	 {
	    switch (surfacelocalid)
	    {
	    	case 0:
	    		surf_neighbor = 3;
	    		break;
	    	case 1:
	    		surf_neighbor = 3;
	    		break;
            case 2:
            	surf_neighbor = 1;
            	break;
            case 3:
            	surf_neighbor = 0;
            	break;
            default:
            	dserror("wrong number!!");
            break;
	    } 
	 }
     }

     //! find the surface at the opposite side
     template <int numnodeparentele>
     void FindOppositeSurface1D(int surfacelocalid,int & surf_neighbor)
     {
    	 dserror("The implementation of FindOppositeSurface1D is not tested 2D elements!");
    	 if ((numnodeparentele == 4) or (numnodeparentele == 8))
    	  {
    	    switch (surfacelocalid)
    	    {
                case 0:
        	        surf_neighbor = 2;
        	    break;
                case 1:
                 	surf_neighbor = 3;
        	    break;
                case 2:
        	        surf_neighbor = 0;
        	    break;
                case 3:
        	        surf_neighbor = 1;
        	    break;
                default:
        	        dserror("wrong number!!");
        	    break;
               }
    	  }
    	  if ((numnodeparentele == 3) or (numnodeparentele == 6))
    	  {
    	    switch (surfacelocalid)
    	    {
                case 0:
        	        surf_neighbor = 1;
        	    break;
                case 1:
                 	surf_neighbor = 0;
        	    break;
                case 2:
        	        surf_neighbor = 0;
        	    break;
                default:
        	        dserror("wrong number!!");
        	    break;
            }
    	  }
      }

      //! find the lines which connect this surface to the surface at the other side
      template <int numnode>
      void FindConnectingLines2D(std::vector< std::vector<int> > connectivity_line_surf,
    		                     int                             masterlocalid,
    		                     int                             slavelocalid,
    		                     std::set<int>                &  p_lines_m,
    		                     std::set<int>                &  p_lines_s,
    		                     int                             p_surf_neighbor_m,
    		                     int                             p_surf_neighbor_s)
      {
    	  if ((numnode == 4) or (numnode == 8))
          {
    		    for (size_t i_line=0; i_line< connectivity_line_surf.size(); i_line++)
    		    {
    		    	//find all adjacent surfaces to this line
    		 	    std::vector<int> line_surfs = connectivity_line_surf.at(i_line);

    		 	    int countsurf_m = 0;
    		 	    int countsurf_s = 0;
    		 	    for (size_t i_surf=0; i_surf<line_surfs.size(); i_surf++)
    		 	    {
    		 		    if ((line_surfs.at(i_surf) != masterlocalid) and (line_surfs.at(i_surf) != p_surf_neighbor_m))
    		 		    {
    		 			    countsurf_m++;
    		 			    // insert the line, where the two surfaces does not belong to it
    		 			    if (countsurf_m == 2)
    		 			      p_lines_m.insert(i_line);

    		 		   }
    		 		   if ((line_surfs.at(i_surf) != slavelocalid) and (line_surfs.at(i_surf) != p_surf_neighbor_s))
    		 		   {
    		 			   countsurf_s++;
    		 			   // insert the line, where the two surfaces does not belong to it
    		 			   if (countsurf_s == 2)
    		 			     p_lines_s.insert(i_line);

    		 		   }
    		 	   }
    		    }
          }
    	  else if (numnode ==3)
    	  {
    		  dserror("The implementation of FindConnectingLines is not tested for Tet4 elements!");

    		  for (size_t i_line=0; i_line< connectivity_line_surf.size(); i_line++)
    		  {
    			  //find all adjacent surfaces to this line
    			  std::vector<int> line_surfs = connectivity_line_surf.at(i_line);

    			  int countsurf_m = 0;
    			  int countsurf_s = 0;
    			  for (size_t i_surf=0; i_surf<line_surfs.size(); i_surf++)
    			  {
    				  // take the line if it belongs to both surfaces (original surface and neighbor)
    				  if ((line_surfs.at(i_surf) == masterlocalid) and (line_surfs.at(i_surf) == p_surf_neighbor_m))
    				  {
    					  countsurf_m++;
    					  // insert the line, where the two surfaces belong to it
    					  if (countsurf_m == 2)
    						  p_lines_m.insert(i_line);
    				  }
    				  if ((line_surfs.at(i_surf) == slavelocalid) and (line_surfs.at(i_surf) == p_surf_neighbor_s))
    				  {
    					  countsurf_s++;
    					  // insert the line, where the two surfaces belong to it
    					  if (countsurf_s == 2)
    						  p_lines_s.insert(i_line);
    				  }
    			  }
    		  }
    	  }
      }
      //! find the lines which connect this line (masterlocalid) to the neighboring line (slavelocalid) at the other side
       template <int numnode>
       void FindConnectingLines1D(//std::vector< std::vector<int> > connectivity_line_surf,
     		                     int                              masterlocalid,
     		                     int                              slavelocalid,
     		                     std::set<int>                 &  p_lines_m,
     		                     std::set<int>                 &  p_lines_s,
     		                     int                              p_surf_neighbor_m,
     		                     int                              p_surf_neighbor_s)
       {
    	  dserror("The implementation of FindConnectingLines is not 2D elements!");
     	  if (numnode == 4)
           {
     		  if (((masterlocalid == 3) and (p_surf_neighbor_m == 1)) or
     				  ((masterlocalid == 1) and (p_surf_neighbor_m == 3)))
     		  {
     			p_lines_m.insert(6);
     			p_lines_m.insert(4);
     		  }
     		  else if (((slavelocalid == 0) and (p_surf_neighbor_s == 2)) or
     				 ((slavelocalid == 2) and (p_surf_neighbor_s == 0)))
     		  {
       			p_lines_s.insert(3);
       			p_lines_s.insert(1);
     		  }
           }
     	  if (numnode == 3)
           {
     		  if (((masterlocalid == 0) and (p_surf_neighbor_m == 1)) or
     				  ((masterlocalid == 1) and (p_surf_neighbor_m == 0)))
     		  {
     			p_lines_m.insert(2);
     		  }
     		  else if (((slavelocalid == 0) and (p_surf_neighbor_s == 2)) or
     				 ((slavelocalid == 2) and (p_surf_neighbor_s == 0)))
     		  {
       			p_lines_s.insert(1);
     		  }
           }
       }

//      //! pointer to parameter list
//      Teuchos::RCP<DRT::ELEMENTS::FluidEleParameter> fldpara_;

      // nodal arrays
      // ------------
      //! node coordinates of parent (master) element
      LINALG::Matrix<nsd_,piel> pxyze_;
      //! array of nodal velocities, intermediate time level
      LINALG::Matrix<nsd_,piel> pevelaf_;
      //! array of nodal velocities, new time level
      LINALG::Matrix<nsd_,piel> pevelnp_;
      //! array of nodal pressure, new time level
      LINALG::Matrix<piel,1> peprenp_;
      //! array of nodal grid displacements, parent element, new time level
      LINALG::Matrix<nsd_,piel> pedispnp_;
      //! array of nodal grid displacements, surface element, new time level
      LINALG::Matrix<nsd_,piel> edispnp_;

      //! node coordinates of neighbor (slave) element
      LINALG::Matrix<nsd_,niel> nxyze_;
      //! array of nodal velocities, intermediate time level
      LINALG::Matrix<nsd_,niel> nevelaf_;
      //! array of nodal velocities, new time level
      LINALG::Matrix<nsd_,niel> nevelnp_;
      //! array of nodal pressure, new time level
      LINALG::Matrix<niel,1> neprenp_;
      //! array of nodal grid displacements, parent element, new time level
      LINALG::Matrix<nsd_,niel> nedispnp_;

      //! node coordinates of boundary element
      LINALG::Matrix<nsd_, iel> xyze_;

      //! linearisation of convection, convective part for parent element
      LINALG::Matrix<piel,1> p_conv_c;
      //! linearisation of convection, convective part for neighbor element
      LINALG::Matrix<niel,1> n_conv_c;


      // shape functions and derivatives, mapping from reference element to actual geometry
      // ----------------------------------------------------------------------------------
      //! transpose of the jacobian matrix of the mapping (r,s,t)->(x,y,z)
      LINALG::Matrix<nsd_,nsd_> pxjm_;
      //! its inverse
      LINALG::Matrix<nsd_,nsd_> pxji_;
      //! vector of shape functions, parent element
      LINALG::Matrix<piel,1> pfunct_;
      //! vector of shape function derivatives in reference coordinate system, parent element
      LINALG::Matrix<nsd_,piel> pderiv_;
      //! vector of shape function derivatives in global coordinate system
      LINALG::Matrix<nsd_,piel> pderxy_;
      //! vector of shape function (2nd) derivatives in reference coordinate system, parent element
      LINALG::Matrix<numderiv2_p,piel> pderiv2_;
      //! vector of shape function (2nd) derivatives in global coordinate system
      LINALG::Matrix<numderiv2_p,piel> pderxy2_;

      // ----------------------------------------------------------------------------------
      //! transpose of the jacobian matrix of the mapping (r,s,t)->(x,y,z)
      LINALG::Matrix<nsd_,nsd_> nxjm_;
      //! its inverse
      LINALG::Matrix<nsd_,nsd_> nxji_;
      //! vector of shape functions, parent element
      LINALG::Matrix<niel,1> nfunct_;
      //! vector of shape function derivatives in reference coordinate system, parent element
      LINALG::Matrix<nsd_,niel> nderiv_;
      //! vector of shape function derivatives in global coordinate system
      LINALG::Matrix<nsd_,niel> nderxy_;
      //! vector of shape function (2nd) derivatives in reference coordinate system, parent element
      LINALG::Matrix<numderiv2_n,niel> nderiv2_;
      //! vector of shape function (2nd) derivatives in global coordinate system
      LINALG::Matrix<numderiv2_n,niel> nderxy2_;


      //! vector of shape functions, boundary element
      LINALG::Matrix<iel ,1> funct_;
      //! vector of shape function derivatives in reference coordinate system, boundary element
      LINALG::Matrix<facensd_, iel> deriv_;
      //! normal vector on surface element (outward pointing from parent element)
      LINALG::Matrix<nsd_, 1> n_;
      //! derivatives of surface in all reference directions
      LINALG::Matrix<facensd_,nsd_> dxyzdrs_;
      //! the metric tensor
      LINALG::Matrix<facensd_,facensd_> metrictensor_;
      //! the area of an infintesimal surface element
      double drs_;



      // values of non geometrical quantities in gausspoints
      // ---------------------------------------------------

      //! local coordinates w.r.t face
      LINALG::Matrix<facensd_,1> xsi_;
      //! velocity in gausspoint, time n+af
      LINALG::Matrix<nsd_,1> velintaf_;
      //! velocity in gausspoint, time n+1
      LINALG::Matrix<nsd_,1> velintnp_;
      //! parent velocity derivatives in gausspoint, time n+af
      LINALG::Matrix<nsd_,nsd_> pvderxyaf_;
      //! neighbor velocity derivatives in gausspoint, time n+af
      LINALG::Matrix<nsd_,nsd_> nvderxyaf_;
      //! parent velocity derivatives in gausspoint, time n+1
      LINALG::Matrix<nsd_,nsd_> pvderxynp_;
      //! neighbor velocity derivatives in gausspoint, time n+1
      LINALG::Matrix<nsd_,nsd_> nvderxynp_;
      //! parent velocity (2nd) derivatives in gausspoint, time n+af
      LINALG::Matrix<nsd_,numderiv2_p> pvderxy2af_;
      //! neighbor velocity (2nd) derivatives in gausspoint, time n+af
      LINALG::Matrix<nsd_,numderiv2_n> nvderxy2af_;

      //! parent velocity (2nd) derivatives in gausspoint, time n+af
      LINALG::Matrix<1,numderiv2_p> ppderxy2af_;
      //! neighbor velocity (2nd) derivatives in gausspoint, time n+af
      LINALG::Matrix<1,numderiv2_n> npderxy2af_;

      //! pressure in gausspoint, time n+1
      double prenp_;
      //! parent pressure derivatives in gausspoint, time n+1
      LINALG::Matrix<nsd_,1> pprederxy_;
      //! neighbor velocity derivatives in gausspoint, time n+1
      LINALG::Matrix<nsd_,1> nprederxy_;

    protected:

      DRT::UTILS::GaussIntegration intpoints_;

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
          const INPAR::FLUID::EOS_TauType tautype,
          double&       tau_grad,
          double&       tau_u,
          double&       tau_div,
          double&       tau_p,
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
