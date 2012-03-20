
#include "fluid_ele_calc_intfaces_stab.H"

#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"

#include "../drt_lib/drt_globalproblem.H"

#include "../drt_cut/cut_position.H"

#include "../drt_mat/newtonianfluid.H"
#include "../drt_mat/carreauyasuda.H"
#include "../drt_mat/modpowerlaw.H"


//-----------------------------------------------------------------
//-----------------------------------------------------------------
//
//                        INTERFACE CLASS
//
//-----------------------------------------------------------------
//-----------------------------------------------------------------

//-----------------------------------------------------------------
//   Allocate one static instance of the internal implementation
//   class for weak dirichlet condition and return pointer to it
//-----------------------------------------------------------------
DRT::ELEMENTS::Fluid3InternalFacesStabilization* DRT::ELEMENTS::Fluid3InternalFacesStabilization::Impl(
  DRT::ELEMENTS::Fluid3Boundary* f3bdry,
  DRT::ELEMENTS::Fluid3*         f3neighbor
  )
{
  switch (f3bdry->Shape())
  {
  // 3D:
  case DRT::Element::quad4:
  {
    static Fluid3InternalSurfaceStabilization<DRT::Element::quad4,DRT::Element::hex8,DRT::Element::hex8>* fsurfq4;

    if(f3bdry->ParentElement()->Shape()==DRT::Element::hex8
        && f3neighbor->Shape()== DRT::Element::hex8)
    {
      if (fsurfq4==NULL)
        fsurfq4 = new Fluid3InternalSurfaceStabilization<DRT::Element::quad4,DRT::Element::hex8,DRT::Element::hex8>();
    }
    else
    {
      dserror("expected combination quad4/hex8/hex8 for surface/parent/neighbor pair");
    }
    return fsurfq4;
  }
  // 2D:
  case DRT::Element::line2:
  {
    dserror("Edgebased stabilization not implemented for 2D elements!");
    break;
//    static Fluid3InternalLineStabilization<DRT::Element::line2,DRT::Element::quad4>* fline2;
//
//    if(f3bdry->ParentElement()->Shape()==DRT::Element::quad4)
//    {
//      if (fline2==NULL)
//      {
//        fline2 = new Fluid3InternalLineStabilization<DRT::Element::line2,DRT::Element::quad4>();
//      }
//    }
//    else
//    {
//      dserror("expected combination line2/quad4 for line/parent pair");
//    }
//
//    return fline2;
  }
  default:
    dserror("shape %d (%d nodes) not supported by internalfaces stabilization", f3bdry->Shape(), f3bdry->NumNode());
  }

  return NULL;
}


//-----------------------------------------------------------------
//-----------------------------------------------------------------
//
//                        IMPLEMENTATION
//
//-----------------------------------------------------------------
//-----------------------------------------------------------------

//-----------------------------------------------------------------
//                       empty constructor
//-----------------------------------------------------------------
template <DRT::Element::DiscretizationType distype,
          DRT::Element::DiscretizationType pdistype,
          DRT::Element::DiscretizationType ndistype>
DRT::ELEMENTS::Fluid3InternalSurfaceStabilization<distype,pdistype,ndistype>::Fluid3InternalSurfaceStabilization()
{
  // pointer to class Fluid3ImplParameter (access to the general parameter)
  f3Parameter_ = DRT::ELEMENTS::FluidEleParameter::Instance();

  return;
}

//-----------------------------------------------------------------
//             evaluate implementation for internal surface stabilization
//-----------------------------------------------------------------
template <DRT::Element::DiscretizationType distype,
          DRT::Element::DiscretizationType pdistype,
          DRT::Element::DiscretizationType ndistype>
int DRT::ELEMENTS::Fluid3InternalSurfaceStabilization<distype,pdistype, ndistype>::EvaluateEdgeBasedStabilization(
  Fluid3Boundary*            surfele       ,
  Fluid3*                    pele,
  Fluid3*                    nele,
  ParameterList&             params        ,
  DRT::Discretization&       discretization,
  vector<int>&               lm,
  vector<int>&               plm            ,
  vector<int>&               nlm            ,
  Epetra_SerialDenseMatrix&  elemat_epetra ,
  Epetra_SerialDenseVector&  elevec_epetra )
{

    bool ghost_penalty   = params.get<bool>("ghost_penalty");
    bool edge_based_stab = params.get<bool>("edge_based_stab");

    if( (!ghost_penalty) and (!edge_based_stab)) dserror("do not call EvaluateEdgeBasedStabilization if no stab is required!");

//    const double timefac      = f3Parameter_->timefac_;     // timefac_ = theta_*dt_;
//    const double timefacpre   = f3Parameter_->timefacpre_;  // special factor for pressure terms in genalpha time integration
//    const double timefacrhs   = f3Parameter_->timefacrhs_;  // factor for rhs (for OST: also theta_*dt_), modified just for genalpha time integration

    // modified time factors

    // for the streamline and divergence stabilization a full matrix pattern is applied
    // fully implicit integration of j_stream(u_h,v_h)
    // Literature: E.Burman, M.A.Fernandez 2009
    // "Finite element methods with symmetric stabilization for the transient convection-diffusion-reaction equation"

    double timefac    = 0.0;
    double timefacpre = 0.0;

    // full matrix pattern (implicit) for streamline and div-stab
    if(not f3Parameter_->is_stationary_)
    {
      timefac    = f3Parameter_->dt_; // set theta = 1.0
      timefacpre = f3Parameter_->dt_; // set theta = 1.0
    }
    else
    {
      timefac    = 1.0;
      timefacpre = 1.0;
    }

    // --------------------------------------------------
    // Reshape element matrices and vectors and init to zero, construct views
    const int eledim = 4*piel + 4*niel;

    elemat_epetra.Shape(eledim,eledim);
    elevec_epetra.Size (eledim);
    //
    LINALG::Matrix<eledim,eledim> elemat(elemat_epetra.A(),true);
    LINALG::Matrix<eledim,     1> elevec(elevec_epetra.A(),true);


    // --------------------------------------------------
    // extract velocities from global distributed vectors

    // velocities (intermediate time step, n+alpha_F)
    RefCountPtr<const Epetra_Vector> velaf = discretization.GetState("velaf");
    if (velaf==null)
      dserror("Cannot get state vector 'velaf'");

    vector<double> mypvelaf((plm).size());
    DRT::UTILS::ExtractMyValues(*velaf,mypvelaf,plm);
    vector<double> mynvelaf((nlm).size());
    DRT::UTILS::ExtractMyValues(*velaf,mynvelaf,nlm);

    // velocities n+1 for parent element and neighbor element
    vector<double> mypvelnp((plm).size());
    vector<double> mynvelnp((nlm).size());

    if((f3Parameter_->timealgo_==INPAR::FLUID::timeint_gen_alpha) or
        (f3Parameter_->timealgo_==INPAR::FLUID::timeint_npgenalpha))
    {
      // velocities (intermediate time step, n+1)
      RefCountPtr<const Epetra_Vector> velnp = discretization.GetState("velnp");
      if (velnp==null)
        dserror("Cannot get state vector 'velnp'");

      DRT::UTILS::ExtractMyValues(*velnp,mypvelnp,plm);
      DRT::UTILS::ExtractMyValues(*velnp,mynvelnp,nlm);
    }
    // mypvelnp = mypvelaf
    else
    {
      DRT::UTILS::ExtractMyValues(*velaf,mypvelnp,plm);
      DRT::UTILS::ExtractMyValues(*velaf,mynvelnp,nlm);
    }


    vector<double> myedispnp ((lm ).size());
    vector<double> mypedispnp((plm).size());
    vector<double> mynedispnp((nlm).size());
    if (surfele->ParentElement()->IsAle())
    {
      // mesh displacements, new time step, n+1
      RefCountPtr<const Epetra_Vector> dispnp
      =
          discretization.GetState("dispnp");
      if (dispnp==null)
      {
        dserror("Cannot get state vector 'dispnp'");
      }

      DRT::UTILS::ExtractMyValues(*dispnp,myedispnp ,lm  );
      DRT::UTILS::ExtractMyValues(*dispnp,mypedispnp,plm);
      DRT::UTILS::ExtractMyValues(*dispnp,mynedispnp,nlm);
    }



    //--------------------------------------------------
    //                GET PARENT DATA
    //--------------------------------------------------

    // extract intermediate velocities
    for(int i=0;i<piel;++i)
    {
      const int fi=4*i;

      pevelaf_(0,i) = mypvelaf[  fi];
      pevelaf_(1,i) = mypvelaf[1+fi];
      pevelaf_(2,i) = mypvelaf[2+fi];
    }

    // extract current velocities and pressure
    for(int i=0;i<piel;++i)
    {
      const int fi=4*i;

      pevelnp_(0,i) = mypvelnp[  fi];
      pevelnp_(1,i) = mypvelnp[1+fi];
      pevelnp_(2,i) = mypvelnp[2+fi];

      peprenp_(  i) = mypvelnp[3+fi];
    }

    if (surfele->ParentElement()->IsAle())
    {
      for (int i=0;i<piel;++i)
      {
        const int fi=4*i;

        pedispnp_(0,i) = mypedispnp[  fi];
        pedispnp_(1,i) = mypedispnp[1+fi];
        pedispnp_(2,i) = mypedispnp[2+fi];
      }

      for (int i=0;i<iel;++i)
      {
        const int fi=4*i;

        edispnp_(0,i) = myedispnp[  fi];
        edispnp_(1,i) = myedispnp[1+fi];
        edispnp_(2,i) = myedispnp[2+fi];
      }
    }

    // extract node coords
    for(int i=0;i<piel;++i)
    {
      pxyze_(0,i)=surfele->ParentElement()->Nodes()[i]->X()[0];
      pxyze_(1,i)=surfele->ParentElement()->Nodes()[i]->X()[1];
      pxyze_(2,i)=surfele->ParentElement()->Nodes()[i]->X()[2];
    }

    if (surfele->ParentElement()->IsAle())
    {
      for (int i=0;i<piel;++i)
      {
        pxyze_(0,i) += pedispnp_(0,i);
        pxyze_(1,i) += pedispnp_(1,i);
        pxyze_(2,i) += pedispnp_(2,i);
      }
    }



    //--------------------------------------------------
     //                GET NEIGHBOR DATA
     //--------------------------------------------------

     // extract intermediate velocities
     for(int i=0;i<niel;++i)
     {
       const int fi=4*i;

       nevelaf_(0,i) = mynvelaf[  fi];
       nevelaf_(1,i) = mynvelaf[1+fi];
       nevelaf_(2,i) = mynvelaf[2+fi];
     }

     // extract current velocities and pressure
     for(int i=0;i<niel;++i)
     {
       const int fi=4*i;

       nevelnp_(0,i) = mynvelnp[  fi];
       nevelnp_(1,i) = mynvelnp[1+fi];
       nevelnp_(2,i) = mynvelnp[2+fi];

       neprenp_(  i) = mynvelnp[3+fi];
     }

     if (nele->IsAle())
     {
       for (int i=0;i<niel;++i)
       {
         const int fi=4*i;

         nedispnp_(0,i) = mynedispnp[  fi];
         nedispnp_(1,i) = mynedispnp[1+fi];
         nedispnp_(2,i) = mynedispnp[2+fi];
       }

     }

     // extract node coords
     for(int i=0;i<niel;++i)
     {
       nxyze_(0,i)=nele->Nodes()[i]->X()[0];
       nxyze_(1,i)=nele->Nodes()[i]->X()[1];
       nxyze_(2,i)=nele->Nodes()[i]->X()[2];
     }

     if (nele->IsAle())
     {
       for (int i=0;i<niel;++i)
       {
         nxyze_(0,i) += nedispnp_(0,i);
         nxyze_(1,i) += nedispnp_(1,i);
         nxyze_(2,i) += nedispnp_(2,i);
       }
     }


     //--------------------------------------------------
     // get material of volume element this surface belongs to
     RCP<MAT::Material> pmat = surfele->ParentElement()->Material();
     RCP<MAT::Material> nmat = nele->Material();

     if( pmat->MaterialType()    != INPAR::MAT::m_carreauyasuda
         && pmat->MaterialType() != INPAR::MAT::m_modpowerlaw
         && pmat->MaterialType() != INPAR::MAT::m_fluid)
       dserror("Material law for parent element is not a fluid");

     if( nmat->MaterialType()    != INPAR::MAT::m_carreauyasuda
         && nmat->MaterialType() != INPAR::MAT::m_modpowerlaw
         && nmat->MaterialType() != INPAR::MAT::m_fluid)
       dserror("Material law for neighbor element is not a fluid");

     // get parent viscosity
     double pkinvisc = 0.0;
     double pdens = 0.0;

     if(pmat->MaterialType() == INPAR::MAT::m_fluid)
     {
       const MAT::NewtonianFluid* actmat = static_cast<const MAT::NewtonianFluid*>(pmat.get());
       // we need the kinematic viscosity (nu ~ m^2/s) here
       pkinvisc = actmat->Viscosity()/actmat->Density();
       pdens = actmat->Density();
       if (actmat->Density() != 1.0)
         dserror("density 1.0 expected: the density need to be included in the linearization terms");
     }
     else
     {
       dserror("up to now I expect a constant viscosity for edge stabilization\n");
     }

     // get neighbor viscosity
     double nkinvisc = 0.0;
     double ndens = 0.0;
     if(nmat->MaterialType() == INPAR::MAT::m_fluid)
     {
       const MAT::NewtonianFluid* actmat = static_cast<const MAT::NewtonianFluid*>(nmat.get());
       // we need the kinematic viscosity here
       nkinvisc = actmat->Viscosity()/actmat->Density();
       ndens = actmat->Density();
       if (actmat->Density() != 1.0)
         dserror("density 1.0 expected: the density need to be included in the linearization terms");
     }
     else
     {
       dserror("up to now I expect a constant viscosity for edge stabilization\n");
     }

     double kinvisc = 0.0;
     double dens = 0.0;

     if(nkinvisc == pkinvisc) kinvisc = pkinvisc;
     else dserror("parent and neighbor element do not have the same viscosity!");

     if(ndens == pdens) dens = pdens;
     else dserror("parent and neighbor element do not have the same density!");

     //--------------------------------------------------
     //          GET BOUNDARY ELEMENT DATA
     //--------------------------------------------------

     // local surface id
//     int psurfaceid = surfele->SurfaceNumber();


     // extract node coords
     for(int i=0;i<iel;++i)
     {
       xyze_(0,i)=surfele->Nodes()[i]->X()[0];
       xyze_(1,i)=surfele->Nodes()[i]->X()[1];
       xyze_(2,i)=surfele->Nodes()[i]->X()[2];
     }

     if (surfele->ParentElement()->IsAle())
     {
       for (int i=0;i<iel;++i)
       {
         xyze_(0,i) += edispnp_(0,i);
         xyze_(1,i) += edispnp_(1,i);
         xyze_(2,i) += edispnp_(2,i);
       }
     }

     //--------------------------------------------------
     // get gausspoints to integrate over boundary element

     // get gauss rule
     DRT::UTILS::GaussRule2D gaussrule=DRT::UTILS::intrule2D_undefined;
     switch (distype)
     {
     case DRT::Element::quad4:
     {
       gaussrule = DRT::UTILS::intrule_quad_4point;
       break;
     }
     case DRT::Element::nurbs9:
     {
       gaussrule = DRT::UTILS::intrule_quad_9point;
       break;
     }
     default:
       dserror("invalid discretization type for Fluid3InternalSurfaceStabilization");
     }

     // gaussian points on surface
     const DRT::UTILS::IntegrationPoints2D intpoints(gaussrule);


     //--------------------------------------------------
     // the gausspoints above have to be mapped to the
     // parent element to be able to evaluate one sided
     // derivatives on the boundary
     //
     // in addition, get information on the orientation of the
     // outward normal

//     Epetra_SerialDenseMatrix pqxg(intpoints.nquad,3);

//     SurfaceGPToParentGP(pqxg     ,
//         intpoints,
//         pdistype ,
//         distype  ,
//         psurfaceid);



     //--------------------------------------------------
     // the gausspoints above have to be mapped to the
     // neighbor element to be able to evaluate one sided
     // derivatives on the boundary
     //
     // in addition, get information on the orientation of the
     // outward normal

//     Epetra_SerialDenseMatrix nqxg(intpoints.nquad,3);
//
//     SurfaceGPToParentGP(nqxg     ,
//         intpoints,
//         ndistype ,
//         distype  ,
//         nsurfaceid);

     //get physical coordinates of Gaussian point



//     double max_vel_L2_norm = params.get<double>("max_vel_L2_norm");

     //------------------------------------------------------------------
     //                       INTEGRATION LOOP
     //------------------------------------------------------------------
     for (int iquad=0;iquad<intpoints.nquad;++iquad)
     {
       // gaussian weight
       const double wquad = intpoints.qwgt[iquad];

       // gaussian point in boundary elements local coordinates
       const double xi    = intpoints.qxg [iquad][0];
       const double eta   = intpoints.qxg [iquad][1];

//       // gaussian point in parent elements local coordinates
//       const double pr     = pqxg(iquad,0);
//       const double ps     = pqxg(iquad,1);
//       const double pt     = pqxg(iquad,2);



       if(!(distype == DRT::Element::nurbs9))
       {
         // ------------------------------------------------
         // shape function derivs of boundary element at gausspoint
         DRT::UTILS::shape_function_2D       (funct_,xi,eta,distype);
         DRT::UTILS::shape_function_2D_deriv1(deriv_,xi,eta,distype);
       }
       else
       {
         dserror("not implemented for nurbs");
       }

       LINALG::Matrix<3,1> x_gp(true);
       x_gp.Multiply(xyze_, funct_);

       // compute local coordinates with respect to neighbor element
       GEO::CUT::Position<ndistype> n_pos( nxyze_, x_gp );
       n_pos.Compute();
       const LINALG::Matrix<3,1> & nqxg = n_pos.LocalCoordinates();


       // gaussian point in neighbor elements local coordinates
       const double nr     = nqxg(0);
       const double ns     = nqxg(1);
       const double nt     = nqxg(2);

       // compute local coordinates with respect to neighbor element
       GEO::CUT::Position<pdistype> p_pos( pxyze_, x_gp );
       p_pos.Compute();
       const LINALG::Matrix<3,1> & pqxg = p_pos.LocalCoordinates();


       // gaussian point in neighbor elements local coordinates
       const double pr     = pqxg(0);
       const double ps     = pqxg(1);
       const double pt     = pqxg(2);

       // ------------------------------------------------
       // compute measure tensor for surface element and the infinitesimal
       // area element drs for the integration

       /*
           |                                              0 1 2
           |                                             +-+-+-+
           |       0 1 2              0...iel-1          | | | | 0
           |      +-+-+-+             +-+-+-+-+          +-+-+-+
           |      | | | | 1           | | | | | 0        | | | | .
           |      +-+-+-+       =     +-+-+-+-+       *  +-+-+-+ .
           |      | | | | 2           | | | | | 1        | | | | .
           |      +-+-+-+             +-+-+-+-+          +-+-+-+
           |                                             | | | | iel-1
           |                                             +-+-+-+
           |
           |       dxyzdrs             deriv              xyze^T
           |
           |
           |                                 +-            -+
           |                                 | dx   dy   dz |
           |                                 | --   --   -- |
           |                                 | dr   dr   dr |
           |     yields           dxyzdrs =  |              |
           |                                 | dx   dy   dz |
           |                                 | --   --   -- |
           |                                 | ds   ds   ds |
           |                                 +-            -+
           |
        */
       dxyzdrs_.MultiplyNT(deriv_,xyze_);
       /*
           |
           |      +-           -+    +-            -+   +-            -+ T
           |      |             |    | dx   dy   dz |   | dx   dy   dz |
           |      |  g11   g12  |    | --   --   -- |   | --   --   -- |
           |      |             |    | dr   dr   dr |   | dr   dr   dr |
           |      |             |  = |              | * |              |
           |      |             |    | dx   dy   dz |   | dx   dy   dz |
           |      |  g21   g22  |    | --   --   -- |   | --   --   -- |
           |      |             |    | ds   ds   ds |   | ds   ds   ds |
           |      +-           -+    +-            -+   +-            -+
           |
           | the calculation of g21 is redundant since g21=g12
        */
       metrictensor_.MultiplyNT(dxyzdrs_,dxyzdrs_);

       /*
                               +--------------+
                              /               |
                sqrtdetg =   /  g11*g22-g12^2
                           \/
        */

       drs_= sqrt(metrictensor_(0,0)*metrictensor_(1,1)
           -
           metrictensor_(0,1)*metrictensor_(1,0));


       // total integration factor
       const double fac = drs_*wquad;

       // ------------------------------------------------
       // compute normal
       if(distype!=DRT::Element::nurbs9)
       {
         double length = 0.0;
         n_(0) = (xyze_(1,1)-xyze_(1,0))*(xyze_(2,2)-xyze_(2,0))
           -
           (xyze_(2,1)-xyze_(2,0))*(xyze_(1,2)-xyze_(1,0));
         n_(1) = (xyze_(2,1)-xyze_(2,0))*(xyze_(0,2)-xyze_(0,0))
           -
           (xyze_(0,1)-xyze_(0,0))*(xyze_(2,2)-xyze_(2,0));
         n_(2) = (xyze_(0,1)-xyze_(0,0))*(xyze_(1,2)-xyze_(1,0))
           -
           (xyze_(1,1)-xyze_(1,0))*(xyze_(0,2)-xyze_(0,0));

         length = n_.Norm2();

         for(int i=0;i<3;++i)
         {
           n_(i)/=length;
         }
       }
       else
       {
         dserror("not implemented for nurbs");
       }


       // determine coordinates of current Gauss point
       LINALG::Matrix<3,1> coordgp;
       for(int i=0;i<3;++i)
       {
         coordgp(i)= 0.0;
       }
       for (int i=0;i<iel;++i)
       {
         for(int j=0;j<3;++j)
         {
           coordgp(j)+=xyze_(j,i)*funct_(i);
         }
       }


       // ------------------------------------------------
       // shape functions and derivs of corresponding parent at gausspoint
       if(!(pdistype == DRT::Element::nurbs27))
       {
         DRT::UTILS::shape_function_3D       (pfunct_,pr,ps,pt,pdistype);
         DRT::UTILS::shape_function_3D_deriv1(pderiv_,pr,ps,pt,pdistype);
       }
       else
       {
         dserror("not implemented for nurbs");
       }

       // ------------------------------------------------
       // shape functions and derivs of corresponding parent at gausspoint
       if(!(ndistype == DRT::Element::nurbs27))
       {
         DRT::UTILS::shape_function_3D       (nfunct_,nr,ns,nt,ndistype);
         DRT::UTILS::shape_function_3D_deriv1(nderiv_,nr,ns,nt,ndistype);
       }
       else
       {
         dserror("not implemented for nurbs");
       }

       //-----------------------------------------------------
       //
       //                       +-------------+
       //                      / /  T       \ |
       //           h = 2 * \ / |  n * G * n |
       //            b       +   \          /
       //

       // get Jacobian matrix and determinant
       pxjm_=0;

       for(int i=0;i<piel;++i)
       {
         for(int rr=0;rr<3;++rr)
         {
           for(int mm=0;mm<3;++mm)
           {
             pxjm_(rr,mm)+=pderiv_(rr,i)*pxyze_(mm,i);
           }
         }
       }

       const double pdet =
           pxjm_(0,0)*pxjm_(1,1)*pxjm_(2,2)+
           pxjm_(0,1)*pxjm_(1,2)*pxjm_(2,0)+
           pxjm_(0,2)*pxjm_(1,0)*pxjm_(2,1)-
           pxjm_(0,2)*pxjm_(1,1)*pxjm_(2,0)-
           pxjm_(0,0)*pxjm_(1,2)*pxjm_(2,1)-
           pxjm_(0,1)*pxjm_(1,0)*pxjm_(2,2);

       // check for degenerated elements
       if (pdet < 0.0)
       {
         dserror("GLOBAL ELEMENT NO.%i\nNEGATIVE JACOBIAN DETERMINANT: %f",
             surfele->ParentElement()->Id(),
             pdet);
       }

       //-----------------------------------------------------
       //
       //                       +-------------+
       //                      / /  T       \ |
       //           h = 2 * \ / |  n * G * n |
       //            b       +   \          /
       //

       // get Jacobian matrix and determinant
       nxjm_=0;

       for(int i=0;i<niel;++i)
       {
         for(int rr=0;rr<3;++rr)
         {
           for(int mm=0;mm<3;++mm)
           {
             nxjm_(rr,mm)+=nderiv_(rr,i)*nxyze_(mm,i);
           }
         }
       }

       const double ndet =
           nxjm_(0,0)*nxjm_(1,1)*nxjm_(2,2)+
           nxjm_(0,1)*nxjm_(1,2)*nxjm_(2,0)+
           nxjm_(0,2)*nxjm_(1,0)*nxjm_(2,1)-
           nxjm_(0,2)*nxjm_(1,1)*nxjm_(2,0)-
           nxjm_(0,0)*nxjm_(1,2)*nxjm_(2,1)-
           nxjm_(0,1)*nxjm_(1,0)*nxjm_(2,2);

       // check for degenerated elements
       if (ndet < 0.0)
       {
         dserror("GLOBAL ELEMENT NO.%i\nNEGATIVE JACOBIAN DETERMINANT: %f",
             nele->Id(),
             ndet);
       }


       //-----------------------------------------------------
       //
       //             compute global first derivates
       //
       /*
         Use the Jacobian and the known derivatives in element coordinate
         directions on the right hand side to compute the derivatives in
         global coordinate directions

               +-                 -+     +-    -+      +-    -+
               |  dx    dy    dz   |     | dN_k |      | dN_k |
               |  --    --    --   |     | ---- |      | ---- |
               |  dr    dr    dr   |     |  dx  |      |  dr  |
               |                   |     |      |      |      |
               |  dx    dy    dz   |     | dN_k |      | dN_k |
               |  --    --    --   |  *  | ---- |   =  | ---- | for all k
               |  ds    ds    ds   |     |  dy  |      |  ds  |
               |                   |     |      |      |      |
               |  dx    dy    dz   |     | dN_k |      | dN_k |
               |  --    --    --   |     | ---- |      | ---- |
               |  dt    dt    dt   |     |  dz  |      |  dt  |
               +-                 -+     +-    -+      +-    -+

        */

       // inverse of jacobian (transposed)
       /*
               +-                 -+     +-                 -+ -1
               |  dr    ds    dt   |     |  dx    dy    dz   |
               |  --    --    --   |     |  --    --    --   |
               |  dx    dx    dx   |     |  dr    dr    dr   |
               |                   |     |                   |
               |  dr    ds    dt   |     |  dx    dy    dz   |
               |  --    --    --   |  =  |  --    --    --   |
               |  dy    dy    dy   |     |  ds    ds    ds   |
               |                   |     |                   |
               |  dr    ds    dt   |     |  dx    dy    dz   |
               |  --    --    --   |     |  --    --    --   |
               |  dz    dz    dz   |     |  dt    dt    dt   |
               +-                 -+     +-                 -+

        */
       pxji_(0,0) = (  pxjm_(1,1)*pxjm_(2,2) - pxjm_(2,1)*pxjm_(1,2))/pdet;
       pxji_(1,0) = (- pxjm_(1,0)*pxjm_(2,2) + pxjm_(2,0)*pxjm_(1,2))/pdet;
       pxji_(2,0) = (  pxjm_(1,0)*pxjm_(2,1) - pxjm_(2,0)*pxjm_(1,1))/pdet;
       pxji_(0,1) = (- pxjm_(0,1)*pxjm_(2,2) + pxjm_(2,1)*pxjm_(0,2))/pdet;
       pxji_(1,1) = (  pxjm_(0,0)*pxjm_(2,2) - pxjm_(2,0)*pxjm_(0,2))/pdet;
       pxji_(2,1) = (- pxjm_(0,0)*pxjm_(2,1) + pxjm_(2,0)*pxjm_(0,1))/pdet;
       pxji_(0,2) = (  pxjm_(0,1)*pxjm_(1,2) - pxjm_(1,1)*pxjm_(0,2))/pdet;
       pxji_(1,2) = (- pxjm_(0,0)*pxjm_(1,2) + pxjm_(1,0)*pxjm_(0,2))/pdet;
       pxji_(2,2) = (  pxjm_(0,0)*pxjm_(1,1) - pxjm_(1,0)*pxjm_(0,1))/pdet;

       nxji_(0,0) = (  nxjm_(1,1)*nxjm_(2,2) - nxjm_(2,1)*nxjm_(1,2))/ndet;
       nxji_(1,0) = (- nxjm_(1,0)*nxjm_(2,2) + nxjm_(2,0)*nxjm_(1,2))/ndet;
       nxji_(2,0) = (  nxjm_(1,0)*nxjm_(2,1) - nxjm_(2,0)*nxjm_(1,1))/ndet;
       nxji_(0,1) = (- nxjm_(0,1)*nxjm_(2,2) + nxjm_(2,1)*nxjm_(0,2))/ndet;
       nxji_(1,1) = (  nxjm_(0,0)*nxjm_(2,2) - nxjm_(2,0)*nxjm_(0,2))/ndet;
       nxji_(2,1) = (- nxjm_(0,0)*nxjm_(2,1) + nxjm_(2,0)*nxjm_(0,1))/ndet;
       nxji_(0,2) = (  nxjm_(0,1)*nxjm_(1,2) - nxjm_(1,1)*nxjm_(0,2))/ndet;
       nxji_(1,2) = (- nxjm_(0,0)*nxjm_(1,2) + nxjm_(1,0)*nxjm_(0,2))/ndet;
       nxji_(2,2) = (  nxjm_(0,0)*nxjm_(1,1) - nxjm_(1,0)*nxjm_(0,1))/ndet;



       //-----------------------------------------------------
       // compute global derivates at integration point
       //
       //   dN    +-----  dN (xi)    dxi
       //     i    \        i           k
       //   --- =   +     ------- * -----
       //   dx     /        dxi      dx
       //     j   +-----       k       j
       //         node k
       //
       // j : direction of derivative x/y/z
       //
       for(int nn=0;nn<piel;++nn)
       {
         for(int rr=0;rr<3;++rr)
         {
           pderxy_(rr,nn)=pxji_(rr,0)*pderiv_(0,nn);

           for(int mm=1;mm<3;++mm)
           {
             pderxy_(rr,nn)+=pxji_(rr,mm)*pderiv_(mm,nn);
           }
         }
       }

       for(int nn=0;nn<niel;++nn)
       {
         for(int rr=0;rr<3;++rr)
         {
           nderxy_(rr,nn)=nxji_(rr,0)*nderiv_(0,nn);

           for(int mm=1;mm<3;++mm)
           {
             nderxy_(rr,nn)+=nxji_(rr,mm)*nderiv_(mm,nn);
           }
         }
       }

       //-----------------------------------------------------
       // get velocities (n+1,i) at integration point
       //
       //                +-----
       //       n+1       \                  n+1
       //    vel   (x) =   +      N (x) * vel
       //                 /        j         j
       //                +-----
       //                node j
       //
       for(int rr=0;rr<3;++rr)
       {
         velintnp_(rr)=pfunct_(0)*pevelnp_(rr,0);
         for(int nn=1;nn<piel;++nn)
         {
           velintnp_(rr)+=pfunct_(nn)*pevelnp_(rr,nn);
         }
       }

       //-----------------------------------------------------
       // get pressure (n+1,i) at integration point
       //
       //                +-----
       //       n+1       \                  n+1
       //    pre   (x) =   +      N (x) * pre
       //                 /        i         i
       //                +-----
       //                node i
       //

       // pressure and velocity are continuous
       prenp_=pfunct_(0)*peprenp_(0);
       for(int nn=1;nn<piel;++nn)
       {
         prenp_+=pfunct_(nn)*peprenp_(nn);
       }

       // pressure derivatives at integration point
       // parent element
       pprederxy_.MultiplyNN(pderxy_, peprenp_);

       // neighbor element
       nprederxy_.MultiplyNN(nderxy_, neprenp_);


       //-----------------------------------------------------
       // get velocities (n+alpha_F,i) at integration point
       //
       //                 +-----
       //       n+af       \                  n+af
       //    vel    (x) =   +      N (x) * vel
       //                  /        j         j
       //                 +-----
       //                 node j
       //
       for(int rr=0;rr<3;++rr)
       {
         velintaf_(rr)=pfunct_(0)*pevelaf_(rr,0);
         for(int nn=1;nn<piel;++nn)
         {
           velintaf_(rr)+=pfunct_(nn)*pevelaf_(rr,nn);
         }
       }



       //-----------------------------------------------------
       // get velocity (n+alpha_F,i) derivatives at integration point
       //
       //       n+af      +-----  dN (x)
       //   dvel    (x)    \        k         n+af
       //   ----------- =   +     ------ * vel
       //       dx         /        dx        k
       //         j       +-----      j
       //                 node k
       //
       // j : direction of derivative x/y/z
       //
       for(int rr=0;rr<3;++rr)
       {
         for(int mm=0;mm<3;++mm)
         {
           // parent element
           pvderxyaf_(rr,mm)=pderxy_(mm,0)*pevelaf_(rr,0);
           for(int nn=1;nn<piel;++nn)
           {
             pvderxyaf_(rr,mm)+=pderxy_(mm,nn)*pevelaf_(rr,nn);
           }

           // neighbor element
           nvderxyaf_(rr,mm)=nderxy_(mm,0)*nevelaf_(rr,0);
           for(int nn=1;nn<niel;++nn)
           {
             nvderxyaf_(rr,mm)+=nderxy_(mm,nn)*nevelaf_(rr,nn);
           }

         }
       }

       //-----------------------------------------------------
       // get convective velocities at integration point
       p_conv_c.MultiplyTN(pderxy_, velintaf_);
       n_conv_c.MultiplyTN(nderxy_, velintaf_);

       // get convective term at old increment for rhs
       LINALG::Matrix<3,1> p_conv_old(true);
       p_conv_old.Multiply(pvderxyaf_,velintaf_);

       LINALG::Matrix<3,1> n_conv_old(true);
       n_conv_old.Multiply(nvderxyaf_,velintaf_);

       LINALG::Matrix<3,3> vderxyaf_diff(true);
       vderxyaf_diff.Update(1.0, nvderxyaf_, -1.0, pvderxyaf_);

//       LINALG::Matrix<3,3> derxyaf_diff(true);
//       derxyaf_diff.Update(1.0, nvderxyaf_, -1.0, pvderxyaf_);

       LINALG::Matrix<3,1> conv_diff(true);
       conv_diff.Multiply(vderxyaf_diff,velintaf_);



       int noffset = 4*piel; // offset where neighbor dofs are placed in elemat and elevec


       double max_pevelnp = 0.0;
       // get the L_inf-norm of the parent's element velocity for stabilization
       for(int r=0; r<3; r++)
       {
         for(int c=0; c<piel; c++)
         {
           if( fabs(pevelnp_(r,c)) > max_pevelnp ) max_pevelnp = fabs(pevelnp_(r,c));
         }
       }

       double max_vel_L2_norm = max_pevelnp;

//       double max_vel_L2_norm = velintaf_.Norm2();



       // additional values for stabilization factors
       Fluid3EdgeBasedStabilization EdgeBasedStabilization;

       // compute the longest element length for the parent element
       EdgeBasedStabilization.h_k(discretization, surfele);



       double tau_u   = 0.0;
       double tau_div = 0.0;
       double tau_p   = 0.0;

       double tau_grad = 0.0;

       double tau_u_lin   = 0.0;
       double tau_p_lin   = 0.0;
       double tau_div_lin = 0.0;


       double timefacfac     = fac*timefac;
       double timefacfac_pre = fac*timefacpre;


       EdgeBasedStabilization.ComputeStabilizationParams(tau_grad, tau_u, tau_div, tau_p, tau_u_lin, tau_div_lin, tau_p_lin, kinvisc, dens, max_vel_L2_norm, timefac);


       // assemble ghost penalty terms for xfluid application
       GhostPenalty( elemat,
                     elevec,
                     noffset,
                     timefacfac,
                     tau_grad,
                     ghost_penalty);


#if(1)

       // EOS stabilization terms for whole Reynolds number regime
       if(edge_based_stab)
       {

         // assemble pressure (EOS) stabilization terms for fluid
         pressureEOS(  elemat,
                       elevec,
                       noffset,
                       timefacfac_pre,
                       tau_p);


         // assemble combined divergence and streamline(EOS) stabilization terms for fluid
         double normal_vel_lin_space = fabs(velintaf_.Dot(n_));
         double div_streamline_tau = tau_div + tau_u * normal_vel_lin_space;

         div_streamline_EOS(  elemat,
                              elevec,
                              noffset,
                              timefacfac,
                              div_streamline_tau,
                              vderxyaf_diff);

       }

#else
       // first version of Burman's EOS stabilization
       if(edge_based_stab)
       {
         // assemble pressure (EOS) stabilization terms for fluid (gradients in normal direction)
         pressureEOSnormal(     elemat,
                                elevec,
                                noffset,
                                timefacfac_pre,
                                tau_p);


         // assemble divergence (EOS) stabilization terms for fluid
         div_EOS(  elemat,
                   elevec,
                   noffset,
                   timefacfac,
                   tau_div,
                   vderxyaf_diff);


         // assemble streamline (EOS) stabilization terms for fluid
         streamline_EOS(      elemat,
                              elevec,
                              noffset,
                              timefacfac,
                              tau_u,
                              vderxyaf_diff,
                              conv_diff);

       }
#endif

#if(0)

       bool u_lin    = false; //true; //true;
       bool div_lin  = false; //true;
       bool pres_lin = false; //true;



       if(ghost_penalty)
       {


         // get grad(u)*n
         LINALG::Matrix<piel,1> p_normal_deriv(true);
         p_normal_deriv.MultiplyTN(pderxy_,n_);

         LINALG::Matrix<niel,1> n_normal_deriv(true);
         n_normal_deriv.MultiplyTN(nderxy_,n_);

         // additional stability of gradients
         // parent column
         for (int ui=0; ui<piel; ++ui)
         {
  //         const double v=timefacfac_densaf*conv_c_(ui);

           // v_parent * u_parent
           //parent row
           for (int vi=0; vi<piel; ++vi)
           {
             for (int ijdim = 0; ijdim <3; ++ijdim) // combined components of u and v
             {
               elemat(vi*4+ijdim  ,ui*4+ijdim) += tau_grad*fac*timefac*p_normal_deriv(vi)*p_normal_deriv(ui);
             }
           }

           for (int vi=0; vi<niel; ++vi)
           {
             for (int ijdim = 0; ijdim <3; ++ijdim) // combined components of u and v
             {
               elemat(noffset + vi*4+ijdim  ,ui*4+ijdim) -= tau_grad*fac*timefac*n_normal_deriv(vi)*p_normal_deriv(ui);
             }
           }
         }

         for (int ui=0; ui<niel; ++ui)
         {
  //         const double v=timefacfac_densaf*conv_c_(ui);

           // v_parent * u_parent
           //parent row
           for (int vi=0; vi<piel; ++vi)
           {
             for (int ijdim = 0; ijdim <3; ++ijdim) // combined components of u and v
             {
               elemat(vi*4+ijdim  ,noffset + ui*4+ijdim) -= tau_grad*fac*timefac*p_normal_deriv(vi)*n_normal_deriv(ui);
             }
           }

           for (int vi=0; vi<niel; ++vi)
           {
             for (int ijdim = 0; ijdim <3; ++ijdim) // combined components of u and v
             {
               elemat(noffset + vi*4+ijdim  ,noffset + ui*4+ijdim) += tau_grad*fac*timefac*n_normal_deriv(vi)*n_normal_deriv(ui);
             }
           }
         }

         LINALG::Matrix<3,1> p_grad_u_n(true);
         p_grad_u_n.Multiply(pvderxyaf_,n_);
         LINALG::Matrix<3,1> n_grad_u_n(true);
         n_grad_u_n.Multiply(nvderxyaf_,n_);

         // v_parent (u_neighbor-u_parent)
         for (int vi=0; vi<piel; ++vi)
         {
           for(int idim = 0; idim <3; ++idim)
           {
  //           elevec(vi*4+idim,0) += tau_u*fac*timefacrhs * (p_conv_c(vi)* (n_conv_old(idim) - p_conv_old(idim)));
             elevec(vi*4+idim,0) += tau_grad*fac*timefac * p_normal_deriv(vi)*(n_grad_u_n(idim)-p_grad_u_n(idim));
           }
         }

         // v_neighbor (u_neighbor-u_parent)
         for (int vi=0; vi<niel; ++vi)
         {
           for(int idim = 0; idim <3; ++idim)
           {
  //           elevec(noffset + vi*4+idim,0) -= tau_u*fac*timefacrhs * (n_conv_c(vi)* (n_conv_old(idim) - p_conv_old(idim)));
             elevec(noffset + vi*4+idim,0) -= tau_grad*fac*timefac * n_normal_deriv(vi)*(n_grad_u_n(idim)-p_grad_u_n(idim));

           }
         }


       }// end if ghost_penalty


       if(edge_based_stab)
       {



         // IMPLEMENTATION OF INTERIOR PENALTY MULTISCALE METHOD for the incompressible Navier-Stokes equation (Burman 2007)

         // pressure part

         //--------------------------------------------------
         // edge stabilization: pressure
         /*
                  //
                  //
                  //             /                             \
                  //            |                               |
                  //  + tau_p * |  |[ grad q ]| , |[ grad p ]|  |
                  //            |                               |
                  //             \                             / surface
                  //
          */

         // grad(p_neighbor) - grad(p_parent)
         LINALG::Matrix<3,1> prederxy_jump(true);
         prederxy_jump.Update(1.0, nprederxy_, -1.0, pprederxy_);

         LINALG::Matrix<piel,piel> pderiv_dyad_pderiv(true);
         pderiv_dyad_pderiv.MultiplyTN(pderxy_, pderxy_);

         LINALG::Matrix<piel,niel> pderiv_dyad_nderiv(true);
         pderiv_dyad_nderiv.MultiplyTN(pderxy_, nderxy_);

         LINALG::Matrix<niel,niel> nderiv_dyad_nderiv(true);
         nderiv_dyad_nderiv.MultiplyTN(nderxy_, nderxy_);


         for (int ui=0; ui<piel; ++ui)
         {
           // q_parent * p_parent
           for (int vi=0; vi<piel; ++vi)
           {
             elemat(vi*4+3  ,ui*4+3) += tau_p*fac*timefacpre*pderiv_dyad_pderiv(vi,ui);
           }

           // q_neighbor * p_parent
           for (int vi=0; vi<niel; ++vi)
           {
             elemat(noffset + vi*4+3,ui*4+3) -= tau_p*fac*timefacpre*pderiv_dyad_nderiv(ui,vi);
           }
         }

         for (int ui=0; ui<niel; ++ui)
         {
           // q_parent * p_neighbor
           for (int vi=0; vi<piel; ++vi)
           {
             elemat(vi*4+3, noffset + ui*4+3) -= tau_p*fac*timefacpre*pderiv_dyad_nderiv(vi,ui);
           }

           // q_neighbor * p_neighbor

           for (int vi=0; vi<niel; ++vi)
           {
             elemat(noffset + vi*4+3, noffset + ui*4+3) += tau_p*fac*timefacpre*nderiv_dyad_nderiv(vi,ui);
           }
         }


         // q_parent (p_neighbor-p_parent)
         for (int vi=0; vi<piel; ++vi)
         {
           elevec(vi*4+3,0) += tau_p*fac*timefacpre * (pderxy_(0,vi)*prederxy_jump(0) +
                                                       pderxy_(1,vi)*prederxy_jump(1) +
                                                       pderxy_(2,vi)*prederxy_jump(2)   );
         }

         // -q_neighbor (p_neighbor-p_parent)
         for (int vi=0; vi<niel; ++vi)
         {
           elevec(noffset + vi*4+3,0) -= tau_p*fac*timefacpre * (nderxy_(0,vi)*prederxy_jump(0) +
                                                                 nderxy_(1,vi)*prederxy_jump(1) +
                                                                 nderxy_(2,vi)*prederxy_jump(2)   );
         }


         // divergence part

         //--------------------------------------------------
         // edge stabilization: divergence (
         /*
                  //
                  //
                  //               /                               \
                  //              |                                 |
                  //  + tau_div * |   |[ grad Du ]| : |[ grad v ]|  |
                  //              |                                 |
                  //               \                               / surface
                  //
          */

         // parent column
         for (int ui=0; ui<piel; ++ui)
         {
           // v_parent * u_parent
           //parent row
           for (int vi=0; vi<piel; ++vi)
           {
             for (int idim = 0; idim <3; ++idim) // combined components of u and v
             {
               for(int jdim = 0; jdim<3; ++jdim) // derivative components
               {
                 elemat(vi*4+idim  ,ui*4+idim) += tau_div*fac*timefac*pderxy_(jdim,vi)*pderxy_(jdim,ui);
               }
             }
           }

           // v_neighbor * u_parent
           //parent row
           for (int vi=0; vi<niel; ++vi)
           {
             for (int idim = 0; idim <3; ++idim) // combined components of u and v
             {
               for(int jdim = 0; jdim<3; ++jdim) // derivative components
               {
                 elemat(noffset + vi*4+idim  ,ui*4+idim) -= tau_div*fac*timefac*nderxy_(jdim,vi)*pderxy_(jdim,ui);
               }
             }
           }
         }

         for (int ui=0; ui<niel; ++ui)
         {
           // v_parent * u_parent
           //parent row
           for (int vi=0; vi<piel; ++vi)
           {
             for (int idim = 0; idim <3; ++idim) // combined components of u and v
             {
               for(int jdim = 0; jdim<3; ++jdim) // derivative components
               {
                 elemat(vi*4+idim  ,noffset + ui*4+idim) -= tau_div*fac*timefac*pderxy_(jdim,vi)*nderxy_(jdim,ui);
               }
             }
           }

           // v_neighbor * u_parent
           //parent row
           for (int vi=0; vi<niel; ++vi)
           {
             for (int idim = 0; idim <3; ++idim) // combined components of u and v
             {
               for(int jdim = 0; jdim<3; ++jdim) // derivative components
               {
                 elemat(noffset + vi*4+idim  , noffset + ui*4+idim) += tau_div*fac*timefac*nderxy_(jdim,vi)*nderxy_(jdim,ui);
               }
             }
           }
         }


         // v_parent (u_neighbor-u_parent)
         for (int vi=0; vi<piel; ++vi)
         {
           for(int idim = 0; idim <3; ++idim)
           {
             for(int jdim = 0; jdim<3;++jdim)
             {
               elevec(vi*4+idim,0) += tau_div*fac*timefac * pderxy_(jdim, vi)*vderxyaf_diff(idim, jdim);
             }
           }
         }

         // v_neighbor (u_neighbor-u_parent)
         for (int vi=0; vi<niel; ++vi)
         {
           for(int idim = 0; idim <3; ++idim)
           {
             for(int jdim = 0; jdim<3;++jdim)
             {
               elevec(noffset + vi*4+idim,0) -= tau_div*fac*timefac * nderxy_(jdim, vi)*vderxyaf_diff(idim, jdim);
             }
           }
         }






         // velocity part

         //--------------------------------------------------
         // edge stabilization: velocity
         /*
                 //
                 //
                 //             /                                              \
                 //            |   lin   i                 i                    |
                 //  + tau_u * | | P  ( u )*n | * |[ grad Du ]| : |[ grad v ]|  |
                 //            |                                                |
                 //             \                                              / surface
                 //
          */

         // projected velocity in linear space
         // TODO: this has to be done for higher order elements!!!
         if(pele->Shape() != DRT::Element::hex8) dserror("implement the velocity projection on linear space for non-hex8 elements");

         double normal_vel_lin_space = fabs(velintaf_.Dot(n_));
         //squared value
//         normal_vel_lin_space = normal_vel_lin_space*normal_vel_lin_space;

         // parent column
         for (int ui=0; ui<piel; ++ui)
         {
           // v_parent * u_parent
           //parent row
           for (int vi=0; vi<piel; ++vi)
           {
             for (int idim = 0; idim <3; ++idim) // combined components of u and v
             {
               for(int jdim = 0; jdim<3; ++jdim) // derivative components
               {
                 elemat(vi*4+idim  ,ui*4+idim) += tau_u*fac*normal_vel_lin_space*timefac*pderxy_(jdim,vi)*pderxy_(jdim,ui);
               }
             }
           }

           // v_neighbor * u_parent
           //parent row
           for (int vi=0; vi<niel; ++vi)
           {
             for (int idim = 0; idim <3; ++idim) // combined components of u and v
             {
               for(int jdim = 0; jdim<3; ++jdim) // derivative components
               {
                 elemat(noffset + vi*4+idim  ,ui*4+idim) -= tau_u*fac*normal_vel_lin_space*timefac*nderxy_(jdim,vi)*pderxy_(jdim,ui);
               }
             }
           }
         }

         for (int ui=0; ui<niel; ++ui)
         {
           // v_parent * u_parent
           //parent row
           for (int vi=0; vi<piel; ++vi)
           {
             for (int idim = 0; idim <3; ++idim) // combined components of u and v
             {
               for(int jdim = 0; jdim<3; ++jdim) // derivative components
               {
                 elemat(vi*4+idim  ,noffset + ui*4+idim) -= tau_u*fac*normal_vel_lin_space*timefac*pderxy_(jdim,vi)*nderxy_(jdim,ui);
               }
             }
           }

           // v_neighbor * u_parent
           //parent row
           for (int vi=0; vi<niel; ++vi)
           {
             for (int idim = 0; idim <3; ++idim) // combined components of u and v
             {
               for(int jdim = 0; jdim<3; ++jdim) // derivative components
               {
                 elemat(noffset + vi*4+idim  , noffset + ui*4+idim) += tau_u*fac*normal_vel_lin_space*timefac*nderxy_(jdim,vi)*nderxy_(jdim,ui);
               }
             }
           }
         }


         // v_parent (u_neighbor-u_parent)
         for (int vi=0; vi<piel; ++vi)
         {
           for(int idim = 0; idim <3; ++idim)
           {
             for(int jdim = 0; jdim<3;++jdim)
             {
               elevec(vi*4+idim,0) += tau_u*fac*normal_vel_lin_space*timefac * pderxy_(jdim, vi)*vderxyaf_diff(idim, jdim);
             }
           }
         }

         // v_neighbor (u_neighbor-u_parent)
         for (int vi=0; vi<niel; ++vi)
         {
           for(int idim = 0; idim <3; ++idim)
           {
             for(int jdim = 0; jdim<3;++jdim)
             {
               elevec(noffset + vi*4+idim,0) -= tau_u*fac*normal_vel_lin_space*timefac * nderxy_(jdim, vi)*vderxyaf_diff(idim, jdim);
             }
           }
         }












         //=========================================================================================================




                // edge stabilization: velocity (Oseen-part, Linearization part I)
                /*
                         //
                         //
                         //             /                                      \
                         //            |   i                   i                |
                         //  + tau_u * |  u * |[ grad Du ]| , u * |[ grad v ]|  |
                         //            |                                        |
                         //             \                                      / surface
                         //
                 */

                // parent column
                for (int ui=0; ui<piel; ++ui)
                {
         //         const double v=timefacfac_densaf*conv_c_(ui);

                  // v_parent * u_parent
                  //parent row
                  for (int vi=0; vi<piel; ++vi)
                  {
                    for (int ijdim = 0; ijdim <3; ++ijdim) // combined components of u and v
                    {
                      elemat(vi*4+ijdim  ,ui*4+ijdim) += tau_u*fac*timefac*p_conv_c(vi)*p_conv_c(ui);
                    }
                  }

                  // v_neighbor * u_parent
                  //parent row
                  for (int vi=0; vi<niel; ++vi)
                  {
                    for (int ijdim = 0; ijdim <3; ++ijdim) // combined components of u and v
                    {
                      elemat(noffset + vi*4+ijdim  ,ui*4+ijdim) -= tau_u*fac*timefac*n_conv_c(vi)*p_conv_c(ui);
                    }
                  }
                }

                for (int ui=0; ui<niel; ++ui)
                {
         //         const double v=timefacfac_densaf*conv_c_(ui);

                  // v_parent * u_neighbor
                  //parent row
                  for (int vi=0; vi<piel; ++vi)
                  {
                    for (int ijdim = 0; ijdim <3; ++ijdim) // combined components of u and v
                    {
                      elemat(vi*4+ijdim  ,noffset + ui*4+ijdim) -= tau_u*fac*timefac*p_conv_c(vi)*n_conv_c(ui);
                    }
                  }

                  // v_neighbor * u_neighbor
                  //parent row
                  for (int vi=0; vi<niel; ++vi)
                  {
                    for (int ijdim = 0; ijdim <3; ++ijdim) // combined components of u and v
                    {
                      elemat(noffset + vi*4+ijdim  ,noffset + ui*4+ijdim) += tau_u*fac*timefac*n_conv_c(vi)*n_conv_c(ui);
                    }
                  }
                }


                // v_parent (u_neighbor-u_parent)
                for (int vi=0; vi<piel; ++vi)
                {
                  for(int idim = 0; idim <3; ++idim)
                  {
         //           elevec(vi*4+idim,0) += tau_u*fac*timefacrhs * (p_conv_c(vi)* (n_conv_old(idim) - p_conv_old(idim)));
                    elevec(vi*4+idim,0) += tau_u*fac*timefacrhs * (p_conv_c(vi)* (conv_diff(idim)));
                  }
                }

                // v_neighbor (u_neighbor-u_parent)
                for (int vi=0; vi<niel; ++vi)
                {
                  for(int idim = 0; idim <3; ++idim)
                  {
         //           elevec(noffset + vi*4+idim,0) -= tau_u*fac*timefacrhs * (n_conv_c(vi)* (n_conv_old(idim) - p_conv_old(idim)));
                    elevec(noffset + vi*4+idim,0) -= tau_u*fac*timefacrhs * (n_conv_c(vi)* (conv_diff(idim)));

                  }
                }



                // edge stabilization: velocity ( Linearization part II)
                /*
                         //
                         //
                         //             /                                      \
                         //            |                i      i                |
                         //  + tau_u * |  Du * |[ grad u ]| , u * |[ grad v ]|  |
                         //            |                                        |
                         //             \                                      / surface
                         //
                 */

                // |[ grad u ]| = derxyaf_diff

                // just parent column because of Du := Du_parent
                for (int ui=0; ui<piel; ++ui)
                {

                  for(int idim=0; idim<3; ++idim) // dimensions of Du
                  {
                    // v_parent * u_parent
                    //parent row
                    for (int vi=0; vi<piel; ++vi)
                    {
                      for (int jdim = 0; jdim <3; ++jdim) // combined components of u and v
                      {
                        elemat(vi*4+jdim  ,ui*4+idim) -= tau_u*fac*timefac*  p_conv_c(vi) * vderxyaf_diff(jdim,idim) * pfunct_(ui);
                      }
                    }

                    // v_neighbor * u_parent
                    //neighbor row
                    for (int vi=0; vi<niel; ++vi)
                    {
                      for (int jdim = 0; jdim <3; ++jdim) // combined components of u and v
                      {
                        elemat(noffset+ vi*4+jdim  ,ui*4+idim) += tau_u*fac*timefac*  n_conv_c(vi) * vderxyaf_diff(jdim,idim) * pfunct_(ui);
                      }
                    }
                  }

                }



                // edge stabilization: velocity ( Linearization part III)
                /*
                         //
                         //
                         //             /                                       \
                         //            |    i           i                        |
                         //  + tau_u * |   u * |[ grad u ]| , Du * |[ grad v ]|  |
                         //            |                                         |
                         //             \                                       / surface
                         //
                 */

                // just parent column because of Du := Du_parent
                for (int ui=0; ui<piel; ++ui)
                {

                  for(int idim=0; idim<3; ++idim) // dimensions of Du
                  {
                    // v_parent * u_parent
                    //parent row
                    for (int vi=0; vi<piel; ++vi)
                    {
                      for (int jdim = 0; jdim <3; ++jdim) // combined components of u and v
                      {
                        elemat(vi*4+jdim  ,ui*4+idim) -= tau_u*fac*timefac*  conv_diff(jdim) * pderxy_(idim,vi) * pfunct_(ui);
                      }
                    }

                    // v_neighbor * u_parent
                    //neighbor row
                    for (int vi=0; vi<niel; ++vi)
                    {
                      for (int jdim = 0; jdim <3; ++jdim) // combined components of u and v
                      {
                        elemat(noffset + vi*4+jdim  ,ui*4+idim) += tau_u*fac*timefac*  conv_diff(jdim) * nderxy_(idim,vi) * pfunct_(ui);
                      }
                    }
                  }

                }


                // edge stabilization: velocity ( Linearization of stabilization parameter)
                /*
                         //
                         //
                 */
                if(u_lin)
                {
                  // just parent column because of Du := Du_parent
                  for (int ui=0; ui<piel; ++ui)
                  {

                    for(int idim=0; idim<3; ++idim) // dimensions of Du
                    {
                      // v_parent * u_parent
                      //parent row
                      for (int vi=0; vi<piel; ++vi)
                      {
                        for (int jdim = 0; jdim <3; ++jdim) // combined components of u and v
                        {
                          elemat(vi*4+jdim  ,ui*4+idim) += tau_u_lin*fac*timefac*  conv_diff(jdim) * pderxy_(jdim,vi) * velintaf_(idim) * pfunct_(ui);
                        }
                      }

                      // v_neighbor * u_parent
                      //neighbor row
                      for (int vi=0; vi<niel; ++vi)
                      {
                        for (int jdim = 0; jdim <3; ++jdim) // combined components of u and v
                        {
                          elemat(noffset + vi*4+jdim  ,ui*4+idim) -= tau_u_lin*fac*timefac*  conv_diff(jdim) * nderxy_(jdim,vi) * velintaf_(idim) * pfunct_(ui);
                        }
                      }
                    }

                  }
                }



                // edge stabilization: divergence
                /*
                         //
                         //
                         //               /                               \
                         //              |                                 |
                         //  + tau_div * |   |[ div(u) ]| ,  |[ div(v) ]|  |
                         //              |                                 |
                         //               \                               / surface
                         //
                 */

                // ... not yet implemented
                // parent column
                for (int ui=0; ui<piel; ++ui)
                {
                  // v_parent * u_parent
                  for (int idim = 0; idim <3; ++idim) // components of u
                  {
                    //parent row
                    for (int vi=0; vi<piel; ++vi)
                    {
                      for (int jdim = 0; jdim <3; ++jdim) // components of v
                      {
                        elemat(vi*4+jdim  ,ui*4+idim) += tau_div*fac*timefac*pderxy_(jdim,vi)*pderxy_(idim,ui);
                      }
                    }
                  }

                  // v_neighbor * u_parent
                  for (int idim = 0; idim <3; ++idim) // components of u
                  {
                    //neighbor row
                    for (int vi=0; vi<niel; ++vi)
                    {
                      for (int jdim = 0; jdim <3; ++jdim) // components of v
                      {
                        elemat(noffset + vi*4+jdim  ,ui*4+idim) -= tau_div*fac*timefac*nderxy_(jdim,vi)*pderxy_(idim,ui);
                      }
                    }
                  }
                }

                // neighbor column
                for (int ui=0; ui<niel; ++ui)
                {
                  // v_parent * u_neighbor
                  for (int idim = 0; idim <3; ++idim) // components of u
                  {
                    //parent row
                    for (int vi=0; vi<piel; ++vi)
                    {
                      for (int jdim = 0; jdim <3; ++jdim) // components of v
                      {
                        elemat(vi*4+jdim  , noffset + ui*4+idim) -= tau_div*fac*timefac*pderxy_(jdim,vi)*nderxy_(idim,ui);
                      }
                    }
                  }

                  // v_neighbor * u_neighbor
                  for (int idim = 0; idim <3; ++idim) // components of u
                  {
                    //neighbor row
                    for (int vi=0; vi<niel; ++vi)
                    {
                      for (int jdim = 0; jdim <3; ++jdim) // components of v
                      {
                        elemat(noffset + vi*4+jdim  ,noffset + ui*4+idim) += tau_div*fac*timefac*nderxy_(jdim,vi)*nderxy_(idim,ui);
                      }
                    }
                  }
                }

                // parent divergence
                double p_div = pvderxyaf_(0,0) + pvderxyaf_(1,1) + pvderxyaf_(2,2);

                // neighbor divergence
                double n_div = nvderxyaf_(0,0) + nvderxyaf_(1,1) + nvderxyaf_(2,2);

                // v_parent (u_neighbor-u_parent)
                for (int vi=0; vi<piel; ++vi)
                {
                  for (int idim = 0; idim <3; ++idim) // components of v
                  {
                    elevec(vi*4+idim,0) += tau_div*fac*timefacrhs * pderxy_(idim,vi)*(n_div - p_div);
                  }
                }

                // -v_neighbor (u_neighbor-u_parent)
                for (int vi=0; vi<niel; ++vi)
                {
                  for (int idim = 0; idim <3; ++idim) // components of v
                  {
                    elevec(noffset + vi*4+idim,0) -= tau_div*fac*timefacrhs * nderxy_(idim,vi)*(n_div - p_div);
                  }
                }


                // divergence stabilization parameter linearization

                if(div_lin)
                {
                  for(int ui=0; ui<piel; ++ui)
                  {
                    for(int idim=0; idim<3; ++idim)
                    {
                      // v_parent * u
                      for (int vi=0; vi<piel; ++vi)
                      {
                        for (int jdim = 0; jdim <3; ++jdim) // components of v
                        {
                          elemat(vi*4+jdim,ui*4+idim) -= tau_div_lin*fac*timefac * pderxy_(jdim,vi)*(n_div - p_div)*velintaf_(idim)*pfunct_(ui);
                        }
                      }
                      // v_neighbor * u
                      for (int vi=0; vi<niel; ++vi)
                      {
                        for (int jdim = 0; jdim <3; ++jdim) // components of v
                        {
                          elemat(noffset + vi*4+jdim,ui*4+idim) += tau_div_lin*fac*timefac * nderxy_(jdim,vi)*(n_div - p_div)*velintaf_(idim)*pfunct_(ui);
                        }
                      }
                    }
                  }

                }


                //--------------------------------------------------
                // edge stabilization: pressure
                /*
                         //
                         //
                         //             /                             \
                         //            |                               |
                         //  + tau_p * |  |[ grad q ]| , |[ grad p ]|  |
                         //            |                               |
                         //             \                             / surface
                         //
                 */

                LINALG::Matrix<piel,piel> pderiv_dyad_pderiv(true);
                pderiv_dyad_pderiv.MultiplyTN(pderxy_, pderxy_);

                LINALG::Matrix<piel,niel> pderiv_dyad_nderiv(true);
                pderiv_dyad_nderiv.MultiplyTN(pderxy_, nderxy_);

                LINALG::Matrix<niel,niel> nderiv_dyad_nderiv(true);
                nderiv_dyad_nderiv.MultiplyTN(nderxy_, nderxy_);

                for (int ui=0; ui<piel; ++ui)
                {
                  for (int vi=0; vi<piel; ++vi)
                  {
                    for(int idim=0; idim<3; idim++)
                    {
                      for(int jdim=0; jdim<3; jdim++)
                      {
                        elemat(vi*4+3  ,ui*4+3) += tau_p*fac*timefacpre*pderxy_(idim,vi)*n_(idim)*n_(jdim)*pderxy_(idim,ui);

                      }
                    }
                  }

                  for (int vi=0; vi<niel; ++vi)
                  {
                    for(int idim=0; idim<3; idim++)
                    {
//                      elemat(noffset+vi*4+3  ,ui*4+3) -= tau_p*fac*timefacpre*nderxy_(idim,vi)*pderxy_(idim,ui);
                      for(int jdim=0; jdim<3; jdim++)
                      {
                        elemat(noffset+vi*4+3  ,ui*4+3) -= tau_p*fac*timefacpre*nderxy_(idim,vi)*n_(idim)*n_(jdim)*pderxy_(idim,ui);
                      }
                    }
                  }
                }
                for (int ui=0; ui<niel; ++ui)
                {
                  for (int vi=0; vi<piel; ++vi)
                  {
                    for(int idim=0; idim<3; idim++)
                    {
                      for(int jdim=0; jdim<3; jdim++)
                      {
                        elemat(vi*4+3  ,noffset+ui*4+3) -= tau_p*fac*timefacpre*pderxy_(idim,vi)*n_(idim)*n_(jdim)*nderxy_(idim,ui);
                      }
                    }
                  }

                  for (int vi=0; vi<niel; ++vi)
                  {
                    for(int idim=0; idim<3; idim++)
                    {
                      for(int jdim=0; jdim<3; jdim++)
                      {
                        elemat(noffset+vi*4+3  ,noffset+ui*4+3) += tau_p*fac*timefacpre*nderxy_(idim,vi)*n_(idim)*n_(jdim)*nderxy_(idim,ui);
                      }
                    }
                  }
                }



                    //                {
                for (int ui=0; ui<piel; ++ui)
                {
                  // q_parent * p_parent
                  for (int vi=0; vi<piel; ++vi)
                  {
                    elemat(vi*4+3  ,ui*4+3) += tau_p*fac*timefacpre*pderiv_dyad_pderiv(vi,ui);
//                    elemat(vi*4+3  ,ui*4+3) += 1.0;
                  }

                  // q_neighbor * p_parent
                  for (int vi=0; vi<niel; ++vi)
                  {
                    elemat(noffset + vi*4+3,ui*4+3) -= tau_p*fac*timefacpre*pderiv_dyad_nderiv(ui,vi);
//                    elemat(noffset + vi*4+3,ui*4+3) -= 1.0;
                  }
                }

                for (int ui=0; ui<niel; ++ui)
                {
                  // q_parent * p_neighbor
                  for (int vi=0; vi<piel; ++vi)
                  {
                    elemat(vi*4+3, noffset + ui*4+3) -= tau_p*fac*timefacpre*pderiv_dyad_nderiv(vi,ui);
//                    elemat(vi*4+3, noffset + ui*4+3) -= -1.0;
                  }

                  // q_neighbor * p_neighbor

                  for (int vi=0; vi<niel; ++vi)
                  {
                    elemat(noffset + vi*4+3, noffset + ui*4+3) += tau_p*fac*timefacpre*nderiv_dyad_nderiv(vi,ui);
//                    elemat(noffset + vi*4+3, noffset + ui*4+3) += -1.0;
                  }
                }

                // grad(p_neighbor) - grad(p_parent)
                LINALG::Matrix<3,1> prederxy_jump(true);
                prederxy_jump.Update(1.0, nprederxy_, -1.0, pprederxy_);

                double prederxy_jump_normal = prederxy_jump.Dot(n_);

                // q_parent (p_neighbor-p_parent)
                for (int vi=0; vi<piel; ++vi)
                {
                  elevec(vi*4+3,0) += tau_p*fac*timefacrhs * (pderxy_(0,vi)*n_(0) +
                                                              pderxy_(1,vi)*n_(1) +
                                                              pderxy_(2,vi)*n_(2)   )*prederxy_jump_normal;
                }

                // -q_neighbor (p_neighbor-p_parent)
                for (int vi=0; vi<niel; ++vi)
                {
                  elevec(noffset + vi*4+3,0) -= tau_p*fac*timefacrhs * (nderxy_(0,vi)*n_(0) +
                                                                        nderxy_(1,vi)*n_(1) +
                                                                        nderxy_(2,vi)*n_(2)   )*prederxy_jump_normal;
                }



                // q_parent (p_neighbor-p_parent)
                for (int vi=0; vi<piel; ++vi)
                {
                  elevec(vi*4+3,0) += tau_p*fac*timefacrhs * (pderxy_(0,vi)*prederxy_jump(0) +
                                                              pderxy_(1,vi)*prederxy_jump(1) +
                                                              pderxy_(2,vi)*prederxy_jump(2)   );
                }

                // -q_neighbor (p_neighbor-p_parent)
                for (int vi=0; vi<niel; ++vi)
                {
                  elevec(noffset + vi*4+3,0) -= tau_p*fac*timefacrhs * (nderxy_(0,vi)*prederxy_jump(0) +
                                                                        nderxy_(1,vi)*prederxy_jump(1) +
                                                                        nderxy_(2,vi)*prederxy_jump(2)   );
                }



                if(pres_lin)
                {
                // linearization of pressure stabilization parameter
                for (int ui=0; ui<piel; ++ui)
                {
                  for(int idim=0; idim<3; ++idim)
                  {
                    // q_parent * u_parent
                    for (int vi=0; vi<piel; ++vi)
                    {
                      for(int jdim=0; jdim<3; ++jdim)
                      {
                        elemat(vi*4+3, ui*4+idim) += tau_p_lin*fac*timefacpre*pderxy_(jdim,vi)*prederxy_jump(jdim)*velintaf_(idim)*pfunct_(ui);
                      }

                    }

                    // q_neighbor * u_parent
                    for (int vi=0; vi<niel; ++vi)
                    {
                      for(int jdim=0; jdim<3; ++jdim)
                      {
                        elemat(noffset + vi*4+3, ui*4+idim) -= tau_p_lin*fac*timefacpre*nderxy_(jdim,vi)*prederxy_jump(jdim)*velintaf_(idim)*pfunct_(ui);
                      }

                    }
                  }

                }

                }// end if pres_lin



       } // end if edge_based_stab
#endif






     } // end gaussloop

//     elemat.Print(cout);


  return 0;
}



template <DRT::Element::DiscretizationType distype,
          DRT::Element::DiscretizationType pdistype,
          DRT::Element::DiscretizationType ndistype>
void DRT::ELEMENTS::Fluid3InternalSurfaceStabilization<distype,pdistype, ndistype>::GhostPenalty(
            LINALG::Matrix<4*piel + 4*niel,4*piel + 4*niel>&  elemat,
            LINALG::Matrix<4*piel + 4*niel,     1>&           elevec,
            int                                               noffset,
            const double &                                    timefacfac,
            double &                                          tau_grad,
            bool &                                            ghost_penalty)
{

  if(ghost_penalty)
  {
    double tau_timefacfac = tau_grad * timefacfac;

    // get grad(u)*n
    LINALG::Matrix<piel,1> p_normal_deriv(true);
    p_normal_deriv.MultiplyTN(pderxy_,n_);

    LINALG::Matrix<niel,1> n_normal_deriv(true);
    n_normal_deriv.MultiplyTN(nderxy_,n_);

    // additional stability of gradients
    // parent column
    for (int ui=0; ui<piel; ++ui)
    {
      for (int ijdim = 0; ijdim <3; ++ijdim) // combined components of u and v
      {
        // v_parent * u_parent
        //parent row
        for (int vi=0; vi<piel; ++vi)
        {
          elemat(vi*4+ijdim  ,ui*4+ijdim) += tau_timefacfac*p_normal_deriv(vi)*p_normal_deriv(ui);
        }

        // neighbor row
        for (int vi=0; vi<niel; ++vi)
        {
          elemat(noffset + vi*4+ijdim  ,ui*4+ijdim) -= tau_timefacfac*n_normal_deriv(vi)*p_normal_deriv(ui);
        }
      }
    }

    for (int ui=0; ui<niel; ++ui)
    {
      for (int ijdim = 0; ijdim <3; ++ijdim) // combined components of u and v
      {
        // v_parent * u_parent
        //parent row
        for (int vi=0; vi<piel; ++vi)
        {
          elemat(vi*4+ijdim  ,noffset + ui*4+ijdim) -= tau_timefacfac*p_normal_deriv(vi)*n_normal_deriv(ui);
        }

        //neighbor row
        for (int vi=0; vi<niel; ++vi)
        {
          elemat(noffset + vi*4+ijdim  ,noffset + ui*4+ijdim) += tau_timefacfac*n_normal_deriv(vi)*n_normal_deriv(ui);
        }
      }

    }

    LINALG::Matrix<3,1> p_grad_u_n(true);
    p_grad_u_n.Multiply(pvderxyaf_,n_);
    LINALG::Matrix<3,1> n_grad_u_n(true);
    n_grad_u_n.Multiply(nvderxyaf_,n_);


    for(int idim = 0; idim <3; ++idim)
    {

      double diff_grad_u_n = tau_timefacfac * (n_grad_u_n(idim)-p_grad_u_n(idim));

      // v_parent (u_neighbor-u_parent)
      for (int vi=0; vi<piel; ++vi)
      {
        elevec(vi*4+idim,0) +=  p_normal_deriv(vi)*diff_grad_u_n;
      }

      // v_neighbor (u_neighbor-u_parent)
      for (int vi=0; vi<niel; ++vi)
      {
        elevec(noffset + vi*4+idim,0) -=  n_normal_deriv(vi)*diff_grad_u_n;
      }
    }



  }// end if ghost_penalty

  return;
}


template <DRT::Element::DiscretizationType distype,
          DRT::Element::DiscretizationType pdistype,
          DRT::Element::DiscretizationType ndistype>
void DRT::ELEMENTS::Fluid3InternalSurfaceStabilization<distype,pdistype, ndistype>::pressureEOS(
            LINALG::Matrix<4*piel + 4*niel,4*piel + 4*niel>&  elemat,
            LINALG::Matrix<4*piel + 4*niel,     1>&           elevec,
            int                                               noffset,
            const double &                                    timefacfacpre,
            double &                                          tau_p)
{

  // pressure part

  //--------------------------------------------------
  // edge stabilization: pressure
  /*
           //
           //
           //             /                             \
           //            |                               |
           //  + tau_p * |  |[ grad q ]| , |[ grad p ]|  |
           //            |                               |
           //             \                             / surface
           //
   */

  // grad(p_neighbor) - grad(p_parent)
  LINALG::Matrix<3,1> prederxy_jump(true);
  prederxy_jump.Update(1.0, nprederxy_, -1.0, pprederxy_);

  LINALG::Matrix<piel,piel> pderiv_dyad_pderiv(true);
  pderiv_dyad_pderiv.MultiplyTN(pderxy_, pderxy_);

  LINALG::Matrix<piel,niel> pderiv_dyad_nderiv(true);
  pderiv_dyad_nderiv.MultiplyTN(pderxy_, nderxy_);

  LINALG::Matrix<niel,niel> nderiv_dyad_nderiv(true);
  nderiv_dyad_nderiv.MultiplyTN(nderxy_, nderxy_);

  double tau_timefacfacpre = tau_p*timefacfacpre;

  for (int ui=0; ui<piel; ++ui)
  {
    // q_parent * p_parent
    for (int vi=0; vi<piel; ++vi)
    {
      elemat(vi*4+3  ,ui*4+3) += tau_timefacfacpre*pderiv_dyad_pderiv(vi,ui);
    }

    // q_neighbor * p_parent
    for (int vi=0; vi<niel; ++vi)
    {
      elemat(noffset + vi*4+3,ui*4+3) -= tau_timefacfacpre*pderiv_dyad_nderiv(ui,vi);
    }
  }

  for (int ui=0; ui<niel; ++ui)
  {
    // q_parent * p_neighbor
    for (int vi=0; vi<piel; ++vi)
    {
      elemat(vi*4+3, noffset + ui*4+3) -= tau_timefacfacpre*pderiv_dyad_nderiv(vi,ui);
    }

    // q_neighbor * p_neighbor

    for (int vi=0; vi<niel; ++vi)
    {
      elemat(noffset + vi*4+3, noffset + ui*4+3) += tau_timefacfacpre*nderiv_dyad_nderiv(vi,ui);
    }
  }


  // q_parent (p_neighbor-p_parent)
  for (int vi=0; vi<piel; ++vi)
  {
    elevec(vi*4+3,0) += tau_timefacfacpre * (pderxy_(0,vi)*prederxy_jump(0) +
                                             pderxy_(1,vi)*prederxy_jump(1) +
                                             pderxy_(2,vi)*prederxy_jump(2)   );
  }

  // -q_neighbor (p_neighbor-p_parent)
  for (int vi=0; vi<niel; ++vi)
  {
    elevec(noffset + vi*4+3,0) -= tau_timefacfacpre * (nderxy_(0,vi)*prederxy_jump(0) +
                                                       nderxy_(1,vi)*prederxy_jump(1) +
                                                       nderxy_(2,vi)*prederxy_jump(2)   );
  }


  return;
}


template <DRT::Element::DiscretizationType distype,
          DRT::Element::DiscretizationType pdistype,
          DRT::Element::DiscretizationType ndistype>
void DRT::ELEMENTS::Fluid3InternalSurfaceStabilization<distype,pdistype, ndistype>::div_streamline_EOS(
            LINALG::Matrix<4*piel + 4*niel,4*piel + 4*niel>&  elemat,
            LINALG::Matrix<4*piel + 4*niel,     1>&           elevec,
            int                                               noffset,
            const double &                                    timefacfac,
            double &                                          tau_div_streamline,
            LINALG::Matrix<3,3>&                              vderxyaf_diff)
{

  // tau_div_streamline = tau_div + tau_u * | P  ( u )*n |   combined divergence and streamline parameter

  //--------------------------------------------------
  // edge stabilization: divergence
  /*
           //
           //
           //               /                               \
           //              |                                 |
           //  + tau_div * |   |[ grad Du ]| : |[ grad v ]|  |
           //              |                                 |
           //               \                               / surface
           //
   */

  //--------------------------------------------------
  // edge stabilization: velocity
  /*
           //
           //
           //             /                                              \
           //            |   lin   i                 i                    |
           //  + tau_u * | | P  ( u )*n | * |[ grad Du ]| : |[ grad v ]|  |
           //            |                                                |
           //             \                                              / surface
           //
  */

  double tau_timefacfac = tau_div_streamline * timefacfac;


  for (int idim = 0; idim <3; ++idim) // combined components of u and v
  {
    for(int jdim = 0; jdim<3; ++jdim) // derivative components
    {
      // parent column
      for (int ui=0; ui<piel; ++ui)
      {
        // v_parent * u_parent
        //parent row
        for (int vi=0; vi<piel; ++vi)
        {
          elemat(vi*4+idim  ,ui*4+idim) += tau_timefacfac*pderxy_(jdim,vi)*pderxy_(jdim,ui);
        }

        // v_neighbor * u_parent
        //parent row
        for (int vi=0; vi<niel; ++vi)
        {
          elemat(noffset + vi*4+idim  ,ui*4+idim) -= tau_timefacfac*nderxy_(jdim,vi)*pderxy_(jdim,ui);
        }
      }

      for (int ui=0; ui<niel; ++ui)
      {
        // v_parent * u_parent
        //parent row
        for (int vi=0; vi<piel; ++vi)
        {
          elemat(vi*4+idim  ,noffset + ui*4+idim) -= tau_timefacfac*pderxy_(jdim,vi)*nderxy_(jdim,ui);
        }

        // v_neighbor * u_parent
        //parent row
        for (int vi=0; vi<niel; ++vi)
        {
          elemat(noffset + vi*4+idim  , noffset + ui*4+idim) += tau_timefacfac*nderxy_(jdim,vi)*nderxy_(jdim,ui);
        }
      }


      // v_parent (u_neighbor-u_parent)
      for (int vi=0; vi<piel; ++vi)
      {
        elevec(vi*4+idim,0) += tau_timefacfac * pderxy_(jdim, vi)*vderxyaf_diff(idim, jdim);
      }

      // v_neighbor (u_neighbor-u_parent)
      for (int vi=0; vi<niel; ++vi)
      {
        elevec(noffset + vi*4+idim,0) -= tau_timefacfac * nderxy_(jdim, vi)*vderxyaf_diff(idim, jdim);
      }
    }
  }


//  // parent column
//  for (int ui=0; ui<piel; ++ui)
//  {
//    // v_parent * u_parent
//    //parent row
//    for (int vi=0; vi<piel; ++vi)
//    {
//      for (int idim = 0; idim <3; ++idim) // combined components of u and v
//      {
//        for(int jdim = 0; jdim<3; ++jdim) // derivative components
//        {
//          elemat(vi*4+idim  ,ui*4+idim) += tau_timefacfac*pderxy_(jdim,vi)*pderxy_(jdim,ui);
//        }
//      }
//    }
//
//    // v_neighbor * u_parent
//    //parent row
//    for (int vi=0; vi<niel; ++vi)
//    {
//      for (int idim = 0; idim <3; ++idim) // combined components of u and v
//      {
//        for(int jdim = 0; jdim<3; ++jdim) // derivative components
//        {
//          elemat(noffset + vi*4+idim  ,ui*4+idim) -= tau_timefacfac*nderxy_(jdim,vi)*pderxy_(jdim,ui);
//        }
//      }
//    }
//  }
//
//  for (int ui=0; ui<niel; ++ui)
//  {
//    // v_parent * u_parent
//    //parent row
//    for (int vi=0; vi<piel; ++vi)
//    {
//      for (int idim = 0; idim <3; ++idim) // combined components of u and v
//      {
//        for(int jdim = 0; jdim<3; ++jdim) // derivative components
//        {
//          elemat(vi*4+idim  ,noffset + ui*4+idim) -= tau_timefacfac*pderxy_(jdim,vi)*nderxy_(jdim,ui);
//        }
//      }
//    }
//
//    // v_neighbor * u_parent
//    //parent row
//    for (int vi=0; vi<niel; ++vi)
//    {
//      for (int idim = 0; idim <3; ++idim) // combined components of u and v
//      {
//        for(int jdim = 0; jdim<3; ++jdim) // derivative components
//        {
//          elemat(noffset + vi*4+idim  , noffset + ui*4+idim) += tau_timefacfac*nderxy_(jdim,vi)*nderxy_(jdim,ui);
//        }
//      }
//    }
//  }
//
//
//  // v_parent (u_neighbor-u_parent)
//  for (int vi=0; vi<piel; ++vi)
//  {
//    for(int idim = 0; idim <3; ++idim)
//    {
//      for(int jdim = 0; jdim<3;++jdim)
//      {
//        elevec(vi*4+idim,0) += tau_timefacfac * pderxy_(jdim, vi)*vderxyaf_diff(idim, jdim);
//      }
//    }
//  }
//
//  // v_neighbor (u_neighbor-u_parent)
//  for (int vi=0; vi<niel; ++vi)
//  {
//    for(int idim = 0; idim <3; ++idim)
//    {
//      for(int jdim = 0; jdim<3;++jdim)
//      {
//        elevec(noffset + vi*4+idim,0) -= tau_timefacfac * nderxy_(jdim, vi)*vderxyaf_diff(idim, jdim);
//      }
//    }
//  }


  return;
}


template <DRT::Element::DiscretizationType distype,
          DRT::Element::DiscretizationType pdistype,
          DRT::Element::DiscretizationType ndistype>
void DRT::ELEMENTS::Fluid3InternalSurfaceStabilization<distype,pdistype, ndistype>::pressureEOSnormal(
            LINALG::Matrix<4*piel + 4*niel,4*piel + 4*niel>&  elemat,
            LINALG::Matrix<4*piel + 4*niel,     1>&           elevec,
            int                                               noffset,
            const double &                                    timefacfacpre,
            double &                                          tau_p)
{

  //--------------------------------------------------
  // edge stabilization: pressure in normal direction
  /*
          //
          //
          //             /                                 \
          //            |                                   |
          //  + tau_p * |  |[ grad q ]|*n , |[ grad p ]|*n  |
          //            |                                   |
          //             \                                 / surface
          //
   */

  LINALG::Matrix<piel,piel> pderiv_dyad_pderiv(true);
  pderiv_dyad_pderiv.MultiplyTN(pderxy_, pderxy_);

  LINALG::Matrix<piel,niel> pderiv_dyad_nderiv(true);
  pderiv_dyad_nderiv.MultiplyTN(pderxy_, nderxy_);

  LINALG::Matrix<niel,niel> nderiv_dyad_nderiv(true);
  nderiv_dyad_nderiv.MultiplyTN(nderxy_, nderxy_);

  for(int idim=0; idim<3; idim++)
  {
    for(int jdim=0; jdim<3; jdim++)
    {
      for (int ui=0; ui<piel; ++ui)
      {
        for (int vi=0; vi<piel; ++vi)
        {
          elemat(vi*4+3  ,ui*4+3) += timefacfacpre*pderxy_(idim,vi)*n_(idim)*n_(jdim)*pderxy_(idim,ui);
        }

        for (int vi=0; vi<niel; ++vi)
        {
          elemat(noffset+vi*4+3  ,ui*4+3) -= timefacfacpre*nderxy_(idim,vi)*n_(idim)*n_(jdim)*pderxy_(idim,ui);
        }
      }
      for (int ui=0; ui<niel; ++ui)
      {
        for (int vi=0; vi<piel; ++vi)
        {
          elemat(vi*4+3  ,noffset+ui*4+3) -= timefacfacpre*pderxy_(idim,vi)*n_(idim)*n_(jdim)*nderxy_(idim,ui);
        }

        for (int vi=0; vi<niel; ++vi)
        {
          elemat(noffset+vi*4+3  ,noffset+ui*4+3) += timefacfacpre*nderxy_(idim,vi)*n_(idim)*n_(jdim)*nderxy_(idim,ui);
        }
      }
    }
  }



  // grad(p_neighbor) - grad(p_parent)
  LINALG::Matrix<3,1> prederxy_jump(true);
  prederxy_jump.Update(1.0, nprederxy_, -1.0, pprederxy_);

  double prederxy_jump_normal = prederxy_jump.Dot(n_);

  // q_parent (p_neighbor-p_parent)
  for (int vi=0; vi<piel; ++vi)
  {
    elevec(vi*4+3,0) += timefacfacpre * (pderxy_(0,vi)*n_(0) +
                                         pderxy_(1,vi)*n_(1) +
                                         pderxy_(2,vi)*n_(2)   )*prederxy_jump_normal;
  }

  // -q_neighbor (p_neighbor-p_parent)
  for (int vi=0; vi<niel; ++vi)
  {
    elevec(noffset + vi*4+3,0) -= timefacfacpre * (nderxy_(0,vi)*n_(0) +
                                                   nderxy_(1,vi)*n_(1) +
                                                   nderxy_(2,vi)*n_(2)   )*prederxy_jump_normal;
  }



  return;
}




template <DRT::Element::DiscretizationType distype,
          DRT::Element::DiscretizationType pdistype,
          DRT::Element::DiscretizationType ndistype>
void DRT::ELEMENTS::Fluid3InternalSurfaceStabilization<distype,pdistype, ndistype>::div_EOS(
            LINALG::Matrix<4*piel + 4*niel,4*piel + 4*niel>&  elemat,
            LINALG::Matrix<4*piel + 4*niel,     1>&           elevec,
            int                                               noffset,
            const double &                                    timefacfac,
            double &                                          tau_div,
            LINALG::Matrix<3,3>&                              vderxyaf_diff)
{


  double tau_timefacfac = tau_div * timefacfac;


  // edge stabilization: divergence
  /*
           //
           //
           //               /                               \
           //              |                                 |
           //  + tau_div * |   |[ div(u) ]| ,  |[ div(v) ]|  |
           //              |                                 |
           //               \                               / surface
           //
   */

  // v_parent * u_parent
  for (int idim = 0; idim <3; ++idim) // components of u
  {
    for (int jdim = 0; jdim <3; ++jdim) // components of v
    {
      // parent column
      for (int ui=0; ui<piel; ++ui)
      {
        //parent row
        for (int vi=0; vi<piel; ++vi)
        {
          elemat(vi*4+jdim  ,ui*4+idim) += tau_timefacfac*pderxy_(jdim,vi)*pderxy_(idim,ui);

        }

        //neighbor row
        for (int vi=0; vi<niel; ++vi)
        {
          elemat(noffset + vi*4+jdim  ,ui*4+idim) -= tau_timefacfac*nderxy_(jdim,vi)*pderxy_(idim,ui);
        }
      }

      // neighbor column
      for (int ui=0; ui<niel; ++ui)
      {
        //parent row
        for (int vi=0; vi<piel; ++vi)
        {
          elemat(vi*4+jdim  , noffset + ui*4+idim) -= tau_timefacfac*pderxy_(jdim,vi)*nderxy_(idim,ui);
        }

        //neighbor row
        for (int vi=0; vi<niel; ++vi)
        {
          elemat(noffset + vi*4+jdim  ,noffset + ui*4+idim) += tau_timefacfac*nderxy_(jdim,vi)*nderxy_(idim,ui);
        }
      }
    }
  }




  // parent divergence
  double p_div = pvderxyaf_(0,0) + pvderxyaf_(1,1) + pvderxyaf_(2,2);

  // neighbor divergence
  double n_div = nvderxyaf_(0,0) + nvderxyaf_(1,1) + nvderxyaf_(2,2);

  for (int idim = 0; idim <3; ++idim) // components of v
  {
    // v_parent (u_neighbor-u_parent)
    for (int vi=0; vi<piel; ++vi)
    {
      elevec(vi*4+idim,0) += tau_timefacfac * pderxy_(idim,vi)*(n_div - p_div);
    }

    // -v_neighbor (u_neighbor-u_parent)
    for (int vi=0; vi<niel; ++vi)
    {
      elevec(noffset + vi*4+idim,0) -= tau_timefacfac * nderxy_(idim,vi)*(n_div - p_div);
    }
  }



  return;
}


template <DRT::Element::DiscretizationType distype,
          DRT::Element::DiscretizationType pdistype,
          DRT::Element::DiscretizationType ndistype>
void DRT::ELEMENTS::Fluid3InternalSurfaceStabilization<distype,pdistype, ndistype>::streamline_EOS(
            LINALG::Matrix<4*piel + 4*niel,4*piel + 4*niel>&  elemat,
            LINALG::Matrix<4*piel + 4*niel,     1>&           elevec,
            int                                               noffset,
            const double &                                    timefacfac,
            double &                                          tau_streamline,
            LINALG::Matrix<3,3>&                              vderxyaf_diff,
            LINALG::Matrix<3,1>&                              conv_diff)
{
  // newton flag, if convective stabilization shall be linearized
  // TODO: add input parameter
  bool newton = false;

  double tau_timefacfac = tau_streamline * timefacfac;


  // edge stabilization: velocity (Oseen-part, Linearization part I)
  /*
         //
         //
         //             /                                      \
         //            |   i                   i                |
         //  + tau_u * |  u * |[ grad Du ]| , u * |[ grad v ]|  |
         //            |                                        |
         //             \                                      / surface
         //
   */

  for (int ijdim = 0; ijdim <3; ++ijdim) // combined components of u and v
  {
    // parent column
    for (int ui=0; ui<piel; ++ui)
    {

      // v_parent * u_parent
      //parent row
      for (int vi=0; vi<piel; ++vi)
      {
          elemat(vi*4+ijdim  ,ui*4+ijdim) += tau_timefacfac*p_conv_c(vi)*p_conv_c(ui);
      }

      // v_neighbor * u_parent
      //parent row
      for (int vi=0; vi<niel; ++vi)
      {
          elemat(noffset + vi*4+ijdim  ,ui*4+ijdim) -= tau_timefacfac*n_conv_c(vi)*p_conv_c(ui);
      }
    }

    for (int ui=0; ui<niel; ++ui)
    {
      // v_parent * u_neighbor
      //parent row
      for (int vi=0; vi<piel; ++vi)
      {
          elemat(vi*4+ijdim  ,noffset + ui*4+ijdim) -= tau_timefacfac*p_conv_c(vi)*n_conv_c(ui);
      }

      // v_neighbor * u_neighbor
      //parent row
      for (int vi=0; vi<niel; ++vi)
      {
          elemat(noffset + vi*4+ijdim  ,noffset + ui*4+ijdim) += tau_timefacfac*n_conv_c(vi)*n_conv_c(ui);
      }
    }


    // v_parent (u_neighbor-u_parent)
    for (int vi=0; vi<piel; ++vi)
    {
        elevec(vi*4+ijdim,0) += tau_timefacfac * (p_conv_c(vi)* (conv_diff(ijdim)));
    }

    // v_neighbor (u_neighbor-u_parent)
    for (int vi=0; vi<niel; ++vi)
    {
        elevec(noffset + vi*4+ijdim,0) -= tau_timefacfac * (n_conv_c(vi)* (conv_diff(ijdim)));
    }
  }



  if(newton)
  {

    // edge stabilization: velocity ( Linearization part II)
    /*
       //
       //
       //             /                                      \
       //            |                i      i                |
       //  + tau_u * |  Du * |[ grad u ]| , u * |[ grad v ]|  |
       //            |                                        |
       //             \                                      / surface
       //
     */

    // |[ grad u ]| = derxyaf_diff


    for(int idim=0; idim<3; ++idim) // dimensions of Du
    {
      for (int jdim = 0; jdim <3; ++jdim) // combined components of u and v
      {
        // just parent column because of Du := Du_parent
        for (int ui=0; ui<piel; ++ui)
        {
          // v_parent * u_parent
          //parent row
          for (int vi=0; vi<piel; ++vi)
          {
            elemat(vi*4+jdim  ,ui*4+idim) -= tau_timefacfac*  p_conv_c(vi) * vderxyaf_diff(jdim,idim) * pfunct_(ui);
          }

          // v_neighbor * u_parent
          //neighbor row
          for (int vi=0; vi<niel; ++vi)
          {
            elemat(noffset+ vi*4+jdim  ,ui*4+idim) += tau_timefacfac*  n_conv_c(vi) * vderxyaf_diff(jdim,idim) * pfunct_(ui);
          }
        }

      }
    }





    // edge stabilization: velocity ( Linearization part III)
    /*
        //
        //
        //             /                                       \
        //            |    i           i                        |
        //  + tau_u * |   u * |[ grad u ]| , Du * |[ grad v ]|  |
        //            |                                         |
        //             \                                       / surface
        //
     */

    for(int idim=0; idim<3; ++idim) // dimensions of Du
    {
      for (int jdim = 0; jdim <3; ++jdim) // combined components of u and v
      {
        // just parent column because of Du := Du_parent
        for (int ui=0; ui<piel; ++ui)
        {
          // v_parent * u_parent
          //parent row
          for (int vi=0; vi<piel; ++vi)
          {
            elemat(vi*4+jdim  ,ui*4+idim) -= tau_timefacfac * conv_diff(jdim) * pderxy_(idim,vi) * pfunct_(ui);
          }

          // v_neighbor * u_parent
          //neighbor row
          for (int vi=0; vi<niel; ++vi)
          {
            elemat(noffset + vi*4+jdim  ,ui*4+idim) += tau_timefacfac * conv_diff(jdim) * nderxy_(idim,vi) * pfunct_(ui);
          }
        }

      }
    }

  }// end if newton




//                  // edge stabilization: velocity ( Linearization of stabilization parameter)
//                  /*
//                           //
//                           //
//                   */
//                  if(u_lin)
//                  {
//                    // just parent column because of Du := Du_parent
//                    for (int ui=0; ui<piel; ++ui)
//                    {
//
//                      for(int idim=0; idim<3; ++idim) // dimensions of Du
//                      {
//                        // v_parent * u_parent
//                        //parent row
//                        for (int vi=0; vi<piel; ++vi)
//                        {
//                          for (int jdim = 0; jdim <3; ++jdim) // combined components of u and v
//                          {
//                            elemat(vi*4+jdim  ,ui*4+idim) += tau_u_lin*fac*timefac*  conv_diff(jdim) * pderxy_(jdim,vi) * velintaf_(idim) * pfunct_(ui);
//                          }
//                        }
//
//                        // v_neighbor * u_parent
//                        //neighbor row
//                        for (int vi=0; vi<niel; ++vi)
//                        {
//                          for (int jdim = 0; jdim <3; ++jdim) // combined components of u and v
//                          {
//                            elemat(noffset + vi*4+jdim  ,ui*4+idim) -= tau_u_lin*fac*timefac*  conv_diff(jdim) * nderxy_(jdim,vi) * velintaf_(idim) * pfunct_(ui);
//                          }
//                        }
//                      }
//
//                    }
//                  }

  return;
}





//-----------------------------------------------------------------
//                          constructor
//                                            (public) schott 02/12
//-----------------------------------------------------------------
DRT::ELEMENTS::Fluid3EdgeBasedStabilization::Fluid3EdgeBasedStabilization()
{
  return;
}






/*---------------------------------------------------------------------------------*
 |  calculation of longest face's diameter over faces of parent-ele   schott 01/12 |
 *---------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Fluid3EdgeBasedStabilization::h_k(
    DRT::Discretization &      discretization,
    Fluid3Boundary*            surfele
    )
{
    // get h_k = max h_e over all faces e, whereas h_e is the diameter of face e

      vector<RCP<DRT::Element> > surfaces = surfele->ParentElement()->Surfaces();

      p_hk_ = 0.0;

      for (unsigned int surf=0; surf<surfaces.size(); ++surf)
      {

        RCP< DRT::Element > actsurf = surfaces[surf];

        vector<int> lm;
        vector<int> lmowner;
        vector<int> lmstride;
        actsurf->LocationVector(discretization,lm,lmowner,lmstride);

        double h_e = 0.0;

        switch (actsurf->Shape())
        {
        case DRT::Element::quad4:
          diameter<DRT::Element::quad4>(actsurf, surfele, h_e, lm, discretization );
        break;
//        case DRT::Element::line3:
//          line_length<DRT::Element::line3>(actline, surfele, act_length, lm, discretization );
//        break;
        default:
          dserror("line type in h_k not templated yet");
        };

        // take the longest surface diameter
        p_hk_ = max(p_hk_, h_e);

      }


  return;
}



//------------------------------------------------------------------------------------------------
// computation of edge-oriented/ghost-penalty stabilization parameter                 schott 12/02
//------------------------------------------------------------------------------------------------
void DRT::ELEMENTS::Fluid3EdgeBasedStabilization::ComputeStabilizationParams(
  double&       tau_grad,
  double&       tau_u,
  double&       tau_div,
  double&       tau_p,
  double&       tau_u_lin,
  double&       tau_div_lin,
  double&       tau_p_lin,
  double&       kinvisc, // kinematic viscosity (nu = mu/rho ~ m^2/2)
  double&       density,
  double&       max_vel_L2_norm,
  const double&       timefac)
{


  // dimensionless factors
  double gamma_grad = 0.0; // scaling factor for ghost penalty stabilization

  double gamma_u   = 0.0; //scaling factor for streamline stabilization
  double gamma_div = 0.0; //scaling factor for divergence stabilization
  double gamma_p   = 0.0; //scaling factor for pressure stabilization


  bool tau_lin = false; // boolian for linearization of tau parameter


  //--------------------------------------------------------------------------------------------------------------
  //                       edge-oriented stabilization (EOS), continuous interior penalty (CIP)
  //--------------------------------------------------------------------------------------------------------------

  //TODO: get the tau-definition via input parameters
//  INPAR::FLUID::TauType_EOS_CIP tautype = INPAR::FLUID::tau_EOS_franca_barrenechea_valentin_wall;
//  INPAR::FLUID::TauType_EOS_CIP tautype = INPAR::FLUID::tau_EOS_braack_burman_2007;
//  INPAR::FLUID::TauType_EOS_CIP tautype = INPAR::FLUID::tau_EOS_burman_fernandez_hansbo_2006;
  INPAR::FLUID::TauType_EOS_CIP tautype = INPAR::FLUID::tau_EOS_burman_fernandez;


  switch(tautype)
  {
  case INPAR::FLUID::tau_EOS_burman_fernandez_hansbo_2006:
  {
    // E.Burman, M.A.Fernandez and P.Hansbo 2006
    // "Edge stabilization for the incompressible Navier-Stokes equations: a continuous interior penalty finite element method"
    //
    // velocity:
    //                              h_K^2 * rho
    //             gamma_u *  -----------------------
    //                           || u ||_0,inf_,K
    //
    // divergence:
    //
    //             gamma_div * h_K^2 * || u ||_0,inf_,K * rho
    //
    // pressure:
    //
    //                              h_K^2                                          || u ||_0,inf_,K   *   h_K
    //  gamma_p * min(1,Re_K) * --------------------------       with    Re_K =  --------------------------------
    //                           || u ||_0,inf_,K  * rho                                      nu
    //

    // element Reynold's number Re_K
    double Re_K = max_vel_L2_norm*p_hk_/ kinvisc;

    // streamline/velocity:
    if(Re_K < 1.0)
    {
      tau_u   = gamma_u * p_hk_*p_hk_*p_hk_ / kinvisc * density;
      tau_u_lin = 0.0;
    }
    else
    {
      tau_u   = gamma_u * p_hk_*p_hk_ / max_vel_L2_norm * density;
      if(tau_lin) tau_u_lin = gamma_u * p_hk_*p_hk_ * density;
    }

    // divergence:
    if(max_vel_L2_norm > 1.0e-14)
    {
      tau_div = gamma_div  * p_hk_* p_hk_ * max_vel_L2_norm * density;
      if(tau_lin) tau_div_lin = gamma_div  * p_hk_* p_hk_ / max_vel_L2_norm * density;
    }
    else
    {
      tau_div = 0.0;
      tau_div_lin = 0.0;
    }

    // pressure stabilization
    // switch between low and high Reynolds numbers
    if(Re_K < 1.0)
    {
      tau_p   = gamma_p    * p_hk_* p_hk_* p_hk_/ (kinvisc * density);
      if(tau_lin) tau_p_lin = 0.0;
    }
    else
    {
      tau_p   = gamma_p    *        p_hk_* p_hk_/ (max_vel_L2_norm * density);
      if(tau_lin) tau_p_lin = gamma_p    *        p_hk_* p_hk_/ (max_vel_L2_norm*max_vel_L2_norm*max_vel_L2_norm * density);
    }
  }
  break;
  case INPAR::FLUID::tau_EOS_burman_fernandez:
  {
    // E.Burman, M.A.Fernandez 2009
    // "Finite element methods with symmetric stabilization for the transient convection-diffusion-reaction equation"
    //
    // velocity:
    //                      1                                         1.0
    //             gamma_u --- *  h_E^2 * rho        with gamma_u = -------
    //                      2                                        100.0
    //
    //--------------------------------------------------------------------------------------------------------------------------
    //
    // E.Burman 2007
    // "Interior penalty variational multiscale method for the incompressible Navier-Stokes equation: Monitoring artificial dissipation"
    //
    // velocity:
    //                                                                  1.0
    //             gamma_u *  h_E^2 * rho            with gamma_u   = -------
    //                                                                 100.0
    // divergence:
    //
    //           gamma_div *  h_E^2 * rho            with gamma_div = 0.05*gamma_u
    //
    // pressure:
    //                                                                1.0
    //             gamma_p *  h_E^2 * rho            with gamma_p = -------
    //                                                               100.0
    //--------------------------------------------------------------------------------------------------------------------------
    // E.Burman, P.Hansbo 2006
    // "Edge stabilization for the generalized Stokes problem: A continuous interior penalty method"
    //
    // pressure:
    //                      1                                             1.0
    //             gamma_p --- *  h_E^(s+1) * rho        with gamma_u = -------
    //                      2                                            100.0
    //
    //                                                   with s=2 (nu>=h) viscous case
    //                                                   with s=1 (nu<h)  non-viscous case
    //
    // divergence:
    //                      1                                             1.0
    //           gamma_div --- *  h_E^(s+1) * rho        with gamma_u = -------
    //                      2                                            100.0
    //
    //                                                   with s=2 (nu>=h) viscous case
    //                                                   with s=1 (nu<h)  non-viscous case
    //
    //                                              1
    //           nu-weighting: gamma*(h_E^2) * -----------  (smoothing between h_E^3/nu and h_E^2 )
    //                                          (1+ nu/h)
    //--------------------------------------------------------------------------------------------------------------------------
// TODO: scaling with density!!!

    //-----------------------------------------------
    // pressure
    bool nu_weighting = true;

//    gamma_p = 0.5 / 100.0;
    gamma_p = 1.0 / 100.0;

    //scaling with h^2 (non-viscous case)
    tau_p = gamma_p * p_hk_*p_hk_;

    // nu-weighting
    if(nu_weighting) // viscous -> non-viscous case
    {
      tau_p /= (1.0 + (kinvisc/p_hk_));
    }
    else //viscous case
    {
      if(kinvisc >= p_hk_) tau_p /= (kinvisc/p_hk_);
    }

    //-----------------------------------------------
    // streamline
    tau_u = tau_p;

    //-----------------------------------------------
    // divergence
    tau_div= 0.05*tau_p;

  }
  break;
  case INPAR::FLUID::tau_EOS_braack_burman_2007:
  {
    dserror("Braack_Burman_2007 tau-def not implemented yet");
  }
  break;
  case INPAR::FLUID::tau_EOS_franca_barrenechea_valentin_wall:
  {
    // stationary definition of stabilization parameters

    // Barrenechea/Valentin, Franca/Valentin
    double mk = 1.0/3.0;
    double Re_K = mk*max_vel_L2_norm*p_hk_/ (2.0*kinvisc);

    double xi = max(1.0, Re_K);

    gamma_p = 1.0/30.0;

    //multibody
    //1.0, 1.0/8.0 not possible, 1.0/10.0 okay, 1.0/50.0, 1.0/100.0 unstable for multibody

    //cylinder2D fine
    // 1.0/10.0, 1.0/20.0 not okay, 1.0/30.0 okay, 1.0/50.0 okay

    tau_p = gamma_p * p_hk_*p_hk_*p_hk_*mk/(4.0*density*kinvisc*xi);

    tau_p = gamma_p * p_hk_*p_hk_*p_hk_*mk/(4.0*density*kinvisc);
  //  tau_p = 1.0/40.0    *        p_hk_* p_hk_/max(1.0, max_vel_L2_norm);

    tau_div = 0.0;

    gamma_u = 1.0/40.0; //gamma_p;
  //  tau_u   = gamma_u * p_hk_*p_hk_*p_hk_*mk/(4.0*density*kinvisc*xi);
    tau_u   = gamma_u * p_hk_*p_hk_/(2.0*density);


    gamma_div = gamma_p*1.0/10.0;//1.0/10.0;
    tau_div = gamma_div  * p_hk_* p_hk_ * max_vel_L2_norm * density; // * min(1.0, Re_K);
  }
  break;
  default: dserror("unknown definition for tau\n %i  ", tautype);
  }


  //--------------------------------------------------------------------------------------------------------------
  //                                               ghost penalty
  //--------------------------------------------------------------------------------------------------------------

  gamma_grad = 1.0;
//  tau_grad = gamma_grad*kinvisc * density * p_hk_;
  tau_grad = gamma_grad*0.00001;



  return;
}


