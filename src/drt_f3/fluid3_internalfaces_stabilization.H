/*----------------------------------------------------------------------*/
/*!
\file fluid3_internalfaces_stabilization.H

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
#ifdef D_FLUID3
#ifdef CCADISCRET

#ifndef FLUID3_INTERNALFACES_STABILIZATION_H
#define FLUID3_INTERNALFACES_STABILIZATION_H

#include "fluid3.H"
#include "../drt_f3_impl/fluid3_impl_parameter.H"

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
    class Fluid3InternalFacesStabilization
    {
    public:
      //! Empty constructor
      Fluid3InternalFacesStabilization() {}
      //! Empty destructor
      virtual ~Fluid3InternalFacesStabilization() {}

      /*!
        \brief Evaluate the surface elements edge-oriented stabilization and ghost penalty

        This class does not provide a definition for this function, it
        is defined in the implementation.

        \author schott (schott@lnm.mw.tum.de), 02/12

      */
      virtual int EvaluateEdgeBasedStabilization(
        Fluid3Boundary*             surfele,            ///< surface Fluid3Boundary element)
        Fluid3*                     ele1,               ///< parent Fluid3 element)
        Fluid3*                     ele2,               ///< neighbor Fluid3 element)
        ParameterList&              params,             ///< parameter list )
        DRT::Discretization&        discretization,     ///< reference to discretization)
        vector<int>&                lm,                 ///< local map for surface element)
        vector<int>&                lm1,                ///< local map for parent element)
        vector<int>&                lm2,                ///< local map for neighbor element)
        Epetra_SerialDenseMatrix&   elemat_epetra,      ///< element matrix combined for parent and neighbor element)
        Epetra_SerialDenseVector&   elevec_epetra       ///< element vector combined for parent and neighbor element)
        ) = 0;


      /*!
        \brief Allocate one static instance of the internal
               implementation class for edge oriented stabilization and
               return pointer to it

        \author schott (schott@lnm.mw.tum.de), 02/12

        \param f3surf (in):   fluid surface element
        \param neighbor (in): fluid3 neighbor element

      */
      static Fluid3InternalFacesStabilization* Impl(DRT::ELEMENTS::Fluid3Boundary* f3surf, DRT::ELEMENTS::Fluid3* neighbor);

    };


    //-----------------------------------------------------------------
    //-----------------------------------------------------------------
    //
    //                        IMPLEMENTATION
    //
    //-----------------------------------------------------------------
    //-----------------------------------------------------------------

    /*!
    \brief Internal Fluid3Surface EOS-implementation
           (surface and parent/neighbor element specific)
    */

    //-----------------------------------------------------------------
    ///
    //-----------------------------------------------------------------
    template<DRT::Element::DiscretizationType distype ,
             DRT::Element::DiscretizationType pdistype,
             DRT::Element::DiscretizationType ndistype>
    class Fluid3InternalSurfaceStabilization: public Fluid3InternalFacesStabilization
    {
    public:
      /// Constructor with number of nodes
      Fluid3InternalSurfaceStabilization();

      /// empty destructor
      ~Fluid3InternalSurfaceStabilization(){return;};

      /// number of nodes
      static const int iel  = DRT::UTILS::DisTypeToNumNodePerEle<distype >::numNodePerElement;

      /// number of parentnodes
      static const int piel = DRT::UTILS::DisTypeToNumNodePerEle<pdistype>::numNodePerElement;

      /// number of parentnodes of neighbor element
      static const int niel = DRT::UTILS::DisTypeToNumNodePerEle<pdistype>::numNodePerElement;


      /*!
        \brief Evaluate EOS stabilization

        This method calculates the contributions to rhs and matrix of
        EOS and ghost penalty stabilization for a generalized alpha
        system.

        Literature:
        Burman Fernandez and Hansbo (2006-2009)

      */
      virtual int EvaluateEdgeBasedStabilization(
        Fluid3Boundary*            surfele,            ///< surface Fluid3Boundary element)
        Fluid3*                    pele,               ///< parent Fluid3 element)
        Fluid3*                    nele,               ///< neighbor Fluid3 element)
        ParameterList&             params,             ///< parameter list )
        DRT::Discretization&       discretization,     ///< reference to discretization)
        vector<int>&               lm,                 ///< local map for surface element)
        vector<int>&               plm,                ///< local map for parent element)
        vector<int>&               nlm,                ///< local map for neighbor element)
        Epetra_SerialDenseMatrix&  elemat_epetra,      ///< element matrix combined for parent and neighbor element)
        Epetra_SerialDenseVector&  elevec_epetra       ///< element vector combined for parent and neighbor element)
      );

    private:

      //! Provide ghost penalty assembly for xfluid case
      void GhostPenalty(
          LINALG::Matrix<4*piel + 4*niel,4*piel + 4*niel>&  elemat,        ///< element matrix combined for parent and neighbor element)
          LINALG::Matrix<4*piel + 4*niel,     1>&           elevec,        ///< element vector combined for parent and neighbor element)
          int                                               noffset,       ///< offset for neighbor element in element matrix and vector
          const double &                                    timefacfac,    ///< (time factor) x (integration factor)
          double &                                          tau_grad,      ///< penalty parameter for ghost penalty terms
          bool &                                            ghost_penalty  ///< boolian if ghost penalty terms has to be assembled
      );

      //! Provide pressure (EOS) stabilization assembly for fluid
      void pressureEOS(
          LINALG::Matrix<4*piel + 4*niel,4*piel + 4*niel>&  elemat,        ///< element matrix combined for parent and neighbor element)
          LINALG::Matrix<4*piel + 4*niel,     1>&           elevec,        ///< element vector combined for parent and neighbor element)
          int                                               noffset,       ///< offset for neighbor element in element matrix and vector
          const double &                                    timefacfacpre, ///< (time factor pressure) x (integration factor)
          double &                                          tau_p          ///< penalty parameter for pressure stabilization terms
      );

      //! Provide divergence and streamline (EOS) stabilization assembly for fluid
      void div_streamline_EOS(
          LINALG::Matrix<4*piel + 4*niel,4*piel + 4*niel>&  elemat,              ///< element matrix combined for parent and neighbor element)
          LINALG::Matrix<4*piel + 4*niel,     1>&           elevec,              ///< element vector combined for parent and neighbor element)
          int                                               noffset,             ///< offset for neighbor element in element matrix and vector
          const double &                                    timefacfac,          ///< (time factor ) x (integration factor)
          double &                                          tau_div_streamline,  ///< combined penalty parameter for divergence and streamline stabilization terms
          LINALG::Matrix<3,3>&                              vderxyaf_diff        ///< velocity derivatives (neighbor-parent) element
      );


      //! Provide pressure (EOS) stabilization assembly for fluid (gradients in normal direction)
      void pressureEOSnormal(
          LINALG::Matrix<4*piel + 4*niel,4*piel + 4*niel>&  elemat,        ///< element matrix combined for parent and neighbor element)
          LINALG::Matrix<4*piel + 4*niel,     1>&           elevec,        ///< element vector combined for parent and neighbor element)
          int                                               noffset,       ///< offset for neighbor element in element matrix and vector
          const double &                                    timefacfacpre, ///< (time factor pressure) x (integration factor)
          double &                                          tau_p          ///< penalty parameter for pressure stabilization terms
      );


      //! Provide divergence (EOS) stabilization assembly for fluid
      void div_EOS(
          LINALG::Matrix<4*piel + 4*niel,4*piel + 4*niel>&  elemat,              ///< element matrix combined for parent and neighbor element)
          LINALG::Matrix<4*piel + 4*niel,     1>&           elevec,              ///< element vector combined for parent and neighbor element)
          int                                               noffset,             ///< offset for neighbor element in element matrix and vector
          const double &                                    timefacfac,          ///< (time factor ) x (integration factor)
          double &                                          tau_div,             ///< penalty parameter for divergence stabilization terms
          LINALG::Matrix<3,3>&                              vderxyaf_diff        ///< velocity derivatives (neighbor-parent) element
      );


      //! Provide streamline (EOS) stabilization assembly for fluid
      void streamline_EOS(
          LINALG::Matrix<4*piel + 4*niel,4*piel + 4*niel>&  elemat,              ///< element matrix combined for parent and neighbor element)
          LINALG::Matrix<4*piel + 4*niel,     1>&           elevec,              ///< element vector combined for parent and neighbor element)
          int                                               noffset,             ///< offset for neighbor element in element matrix and vector
          const double &                                    timefacfac,          ///< (time factor ) x (integration factor)
          double &                                          tau_streamline,      ///< combined penalty parameter for streamline stabilization terms
          LINALG::Matrix<3,3>&                              vderxyaf_diff,       ///< velocity derivatives (neighbor-parent) element
          LINALG::Matrix<3,1>&                              conv_diff            ///< convective velocity (neighbor-parent) element
      );


      //! pointer to parameter list
      DRT::ELEMENTS::Fluid3ImplParameter* f3Parameter_;

      // nodal arrays
      // ------------
      //! node coordinates of parent element
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
    class Fluid3EdgeBasedStabilization
    {
    public:


      /*! \brief constructor
      */
      Fluid3EdgeBasedStabilization();

      /// Empty destructor
      virtual ~Fluid3EdgeBasedStabilization() {}




      template <DRT::Element::DiscretizationType surf_distype>
       void diameter(    RCP<DRT::Element>                     actsurf,
                         Fluid3Boundary*                       surfele,
                         double &                              h_e,
                         vector<int> &                         lm,
                         DRT::Discretization&                  discretization)
       {
         const size_t numnode = DRT::UTILS::DisTypeToNumNodePerEle<surf_distype>::numNodePerElement;

         vector<double> mysurfdispnp ((lm ).size());
         if (surfele->ParentElement()->IsAle())
         {
           // mesh displacements, new time step, n+1
           RefCountPtr<const Epetra_Vector> dispnp = discretization.GetState("dispnp");
           if (dispnp==null)
           {
             dserror("Cannot get state vector 'dispnp'");
           }

           DRT::UTILS::ExtractMyValues(*dispnp,mysurfdispnp ,lm  );
         }


         LINALG::Matrix<3,numnode> xyz_surf(true);
         LINALG::Matrix<3,numnode> psurfdispnp (true);

         // extract node coords
         for(int i=0;i<(int)numnode;++i)
         {
           xyz_surf(0,i)=actsurf->Nodes()[i]->X()[0];
           xyz_surf(1,i)=actsurf->Nodes()[i]->X()[1];
           xyz_surf(2,i)=actsurf->Nodes()[i]->X()[2];
         }

         if (surfele->ParentElement()->IsAle())
         {
           for (int i=0;i<(int)numnode;++i)
           {
             const int fi=4*i;

             psurfdispnp(0,i) = mysurfdispnp[  fi];
             psurfdispnp(1,i) = mysurfdispnp[1+fi];
             psurfdispnp(2,i) = mysurfdispnp[2+fi];
           }
         }


         if (surfele->ParentElement()->IsAle())
         {
           for (size_t i=0;i<numnode;++i)
           {
             xyz_surf(0,i) += psurfdispnp(0,i);
             xyz_surf(1,i) += psurfdispnp(1,i);
             xyz_surf(2,i) += psurfdispnp(2,i);
           }
         }

         double diam1 = 0.0;
         double diam2 = 0.0;

         if(numnode == 4)
         {
           for(int i=0; i<3; i++)
           {
             // line nodes (0,2)
             diam1 += (xyz_surf(i,2) - xyz_surf(i,0))*(xyz_surf(i,2) - xyz_surf(i,0));

             // line nodes (1,3)
             diam2 += (xyz_surf(i,3) - xyz_surf(i,1))*(xyz_surf(i,3) - xyz_surf(i,1));
           }
           diam1 = sqrt(diam1);
           diam2 = sqrt(diam2);
         }
         else dserror("wrong number of nodes");

         h_e = max(diam1, diam2);

         if(h_e <= 0.0) dserror("negative or zero diameter for current face!");

         return;
       }




//      template <DRT::Element::DiscretizationType line_distype>
//       void line_length( RCP<DRT::Element>          actline,
//                         Fluid3Boundary*            surfele,
//                         double &                   length,
//                         vector<int> &              lm,
//                         DRT::Discretization&       discretization)
//       {
//         const size_t numnode = DRT::UTILS::DisTypeToNumNodePerEle<line_distype>::numNodePerElement;
//
//         vector<double> mylinedispnp ((lm ).size());
//         if (surfele->ParentElement()->IsAle())
//         {
//           // mesh displacements, new time step, n+1
//           RefCountPtr<const Epetra_Vector> dispnp = discretization.GetState("dispnp");
//           if (dispnp==null)
//           {
//             dserror("Cannot get state vector 'dispnp'");
//           }
//
//           DRT::UTILS::ExtractMyValues(*dispnp,mylinedispnp ,lm  );
//         }
//
//
//         LINALG::Matrix<3,numnode> xyz_line (true);
//         LINALG::Matrix<3,numnode> plinedispnp (true);
//
//         // extract node coords
//         for(int i=0;i<(int)numnode;++i)
//         {
//           xyz_line(0,i)=actline->Nodes()[i]->X()[0];
//           xyz_line(1,i)=actline->Nodes()[i]->X()[1];
//           xyz_line(2,i)=actline->Nodes()[i]->X()[2];
//         }
//
//         if (surfele->ParentElement()->IsAle())
//         {
//           for (int i=0;i<(int)numnode;++i)
//           {
//             const int fi=4*i;
//
//             plinedispnp(0,i) = mylinedispnp[  fi];
//             plinedispnp(1,i) = mylinedispnp[1+fi];
//             plinedispnp(2,i) = mylinedispnp[2+fi];
//           }
//         }
//
//
//         if (surfele->ParentElement()->IsAle())
//         {
//           for (int i=0;i<numnode;++i)
//           {
//             xyz_line(0,i) += plinedispnp(0,i);
//             xyz_line(1,i) += plinedispnp(1,i);
//             xyz_line(2,i) += plinedispnp(2,i);
//           }
//         }
//
//         if(numnode == 2)
//         {
//           for(int i=0; i<3; i++)
//           {
//             length += (xyz_line(i,1) - xyz_line(i,0))*(xyz_line(i,1) - xyz_line(i,0));
//           }
//           length = sqrt(length);
//         }
//         else if(numnode == 3)
//         {
//           for(int i=0; i<3; i++)
//           {
//             length += (xyz_line(i,2) - xyz_line(i,0))*(xyz_line(i,2) - xyz_line(i,0));
//           }
//           length = sqrt(length);
//         }
//         else dserror("this number of points of lines is not 2 or 3");
//
//         return;
//       }



      /*!
       \brief computation of longest diameter of all faces in the parent element
      */
      void h_k(
          DRT::Discretization &      discretization,
          Fluid3Boundary*            surfele);



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
          const double&       timefac);


    private:

      //! longest edge in parent element
      double p_hk_;



    };



#if(0)
    /*!
     \brief Internal Fluid2Line weak Dirichlet implementation
            (line and parent element specific)
     */

     //-----------------------------------------------------------------
     ///
     //-----------------------------------------------------------------
     template<DRT::Element::DiscretizationType distype ,
              DRT::Element::DiscretizationType pdistype>
     class Fluid3LineWeakDBC: public Fluid3BoundaryWeakDBCInterface
     {
     public:
       /// empty constructor
       Fluid3LineWeakDBC();

       /// empty destructor
       ~Fluid3LineWeakDBC(){return;};

       /// number of nodes
       static const int iel  = DRT::UTILS::DisTypeToNumNodePerEle<distype >::numNodePerElement;

       /// number of parentnodes
       static const int piel = DRT::UTILS::DisTypeToNumNodePerEle<pdistype>::numNodePerElement;


       /*!
         \brief Evaluate weakly imposed Dirichlet conditions

         This method calculates the contributions to rhs and matrix of
         weak Dirichlet boundary conditions for a generalized alpha
         system.

         Literature:
         Weak imposition of Dirichlet boundary conditions in fluid mechanics
         Y. Bazilevs, T.J.R. Hughes
         Computers & Fluids 36 (2007) 12–26

       */
       virtual int EvaluateWeakDBC(
         Fluid3Boundary*                lineele       ,
         ParameterList&             params        ,
         DRT::Discretization&       discretization,
         vector<int>&               lm            ,
         Epetra_SerialDenseMatrix&  elemat_epetra ,
         Epetra_SerialDenseVector&  elevec_epetra );

     private:

        //! pointer to parameter list
        DRT::ELEMENTS::Fluid3ImplParameter* f3Parameter_;

       // nodal arrays
       // ------------
       //! node coordinates of parent element
       LINALG::Matrix<2,piel> pxyze_;
       //! array of nodal velocities, intermediate time level
       LINALG::Matrix<2,piel> pevelaf_;
       //! array of nodal velocities, new time level
       LINALG::Matrix<2,piel> pevelnp_;
       //! array of nodal pressure, new time level
       LINALG::Matrix<piel,1> peprenp_;
       //! array of nodal grid displacements, parent element, new time level
       LINALG::Matrix<2,piel> pedispnp_;
       //! array of nodal grid displacements, surface element, new time level
       LINALG::Matrix<2,piel> edispnp_;
       //! node coordinates of boundary element
       LINALG::Matrix<2, iel> xyze_;



       // shape functions and derivatives, mapping from reference element to actual geometry
       // ----------------------------------------------------------------------------------
       //! transpose of the jacobian matrix of the mapping (r,s,t)->(x,y,z)
       LINALG::Matrix<2,2>    pxjm_;
       //! its inverse
       LINALG::Matrix<2,2>    pxji_;
       //! vector of shape functions, parent element
       LINALG::Matrix<piel,1> pfunct_;
       //! vector of shape function derivatives in reference coordinate system, parent element
       LINALG::Matrix<2,piel> pderiv_;
       //! vector of shape function derivatives in global coordinate system
       LINALG::Matrix<2,piel> pderxy_;


       //! vector of shape functions, boundary element
       LINALG::Matrix<iel ,1> funct_;
       //! vector of shape function derivatives in reference coordinate system, boundary element
       LINALG::Matrix<iel, 1> deriv_;
       //! vector for outward normal, boundary element
       LINALG::Matrix<2,   1> n_;
       //! derivatives of surface in reference direction
       LINALG::Matrix<2,   1> dxydr_;
       //! the length of an infintesimal line element
       double                 dr_;



       // values of non geometrical quantities in gausspoints
       // ---------------------------------------------------

       //! velocity in gausspoint, time n+af
       LINALG::Matrix<2,1>    velintaf_;
       //! velocity in gausspoint, time n+1
       LINALG::Matrix<2,1>    velintnp_;
       //! velocity derivatives in gausspoint, time n+af
       LINALG::Matrix<2,2>    vderxyaf_;
       //! pressure in gausspoint, time n+1
       double                 prenp_;

     };
#endif

#if(0)
    //-----------------------------------------------------------------
    //-----------------------------------------------------------------
    //
    //                    SPALDINGS LAW OF THE WALL
    //
    //-----------------------------------------------------------------
    //-----------------------------------------------------------------

    /*!
      \brief a class providing a method for the dynamic computation of
             the penalty parameter tauB for the weak Dirichlet boundary
             condition. It is computed such that the resulting traction
             corresponds to Spaldings formulation of the law of the
             wall

      Literature:
      Weak Dirichlet Boundary Conditions for Wall-Bounded Turbulent Flows
      Y. Bazilevs, C. Michler, V.M. Calo, and T.J.R. Hughes
      Computer Methods in Applied Mechanics and Engineering
      Volume 196, Issues 49-52, 1 November 2007, Pages 4853-4862

      and references therein

    */
    class Fluid3SurfaceWeakDBCSpaldingsLaw
    {
    public:
      /*! \brief constructor

      \param chi_in (in) empirical paramter of Spaldings law, usually 0.4
      \param B_in   (in) empirical paramter of Spaldings law, usually 5.5

      */
      Fluid3SurfaceWeakDBCSpaldingsLaw(const double chi_in,
                                       const double B_in  );

      /// Empty destructor
      virtual ~Fluid3SurfaceWeakDBCSpaldingsLaw() {}

      /*! \brief dynamic computation of the penalty paramter using Spaldings
                 law


      \param  tau_B (in/out) start/final value for penalty parameter tau_B
      \param  normu (in)     norm velocity in this point on surface
      \param  h     (in)     cell size in wall normal direction
      \param  Cb    (in)     constant originating from weak Dirichlet
                             framework
      \param  visc  (in)     kinematic viscosity

      */
      void ComputeTauBUsingSpaldingsLaw(
        double&       tau_B,
        const double& normu,
        const double& h    ,
        const double& Cb   ,
        const double& visc );

    private:
      //! first parameter of the Spalding model
      double chi_;

      //! second parameter of the Spalding model
      double B_;

      /*! \brief evaluate the residual of Spaldings law of the wall



      <pre>

                               +      / +\
                        res = y  - f | u  |
                                      \  /

                                  /                                                  \
                                 |                           /     +\ 2    /     +\ 3 |
                                 |       +                  | chi*u  |    | chi*u  |  |
    +     / +\     +    -chi*B   |  chi*u               +    \      /      \      /   |
   y = f | u  | = u  + e       * | e       - 1.0 - chi*u  - ----------- - ----------- |
          \  /                   |                              2.0           6.0     |
                                 |                                                    |
                                  \                                                  /


      </pre>

      \param  y     (in)     thickness of modeled boundary layer
      \param  nu    (in)     kinematic viscosity
      \param  tau_B (in)     value for penalty parameter tau_B
      \param  normu (in)     norm velocity in this point on surface

      \return residual of Spaldings equation
      */
      double SpaldingResidual(
        const double y     ,
        const double nu    ,
        const double tau_B ,
        const double normu
        );

      /*! \brief evaluate the residuals Jacobian of Spaldings law of the wall

      \param  y     (in)     thickness of modeled boundary layer
      \param  nu    (in)     kinematic viscosity
      \param  tau_B (in)     value for penalty parameter tau_B
      \param  normu (in)     norm velocity in this point on surface

      \return derivative of residual of Spaldings equation w.r.t. tau_B

      */
      double JacobianSpaldingResidual(
        const double y    ,
        const double nu   ,
        const double tau_B,
        const double normu);

      /*! \brief compute dimensionless velocity u+

      <pre>

                                               +-------------+
                     ||  n+af ||              /  ||  n+af||
           +         || u     ||             /   || u    ||
          u  = ---------------------- =     /   -------------
                 +------------------+    \ /        tau
                / tau  * ||  n+af ||      v            B
               v     B   || u     ||

      </pre>

      \param  tau_B (in)     value for penalty parameter tau_B
      \param  normu (in)     norm velocity in this point on surface

      \return u+

       */
      double Uplus(const double normu,
                   const double tau_B);

      /*! \brief compute dimensionless thickness of modeled layer y+

      <pre>
                   +-------------+              +-------------+
                  / tau  * ||u||          h *  / tau  * ||u||
           +     v     B                   b  v     B
          y  =  ------------------- * y = ---------------------- =
                        nu                      C  * nu
                                                 b
                  +-------------+
                 / tau  * ||u||
                v     B
             = -------------------
                     tau
                        B,0

      </pre>

      \param  normu (in)     norm velocity in this point on surface
      \param  tau_B (in)     value for penalty parameter tau_B
      \param  visc  (in)     kinematic viscosity
      \param  y     (in)     thickness of modeled boundary layer

      \return y+
       */
      double Yplus(const double normu,
                   const double tau_B,
                   const double visc ,
                   const double y    );

    };
#endif

  }
}

#endif

#endif
#endif
