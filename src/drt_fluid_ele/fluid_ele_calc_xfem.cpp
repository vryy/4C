/*----------------------------------------------------------------------*/
/*!
\file fluid_ele_calc_xfem.cpp

\brief Internal implementation of XFluid element interface coupling

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

#include "../drt_cut/cut_boundarycell.H"
#include "../drt_cut/cut_position.H"

#include "../drt_geometry/position_array.H"

#include "../drt_inpar/inpar_xfem.H"

#include "../linalg/linalg_utils.H"

#include "fluid_ele.H"
#include "fluid_ele_parameter.H"
#include "fluid_ele_calc_xfem_coupling.H"
#include "fluid_ele_calc_xfem.H"

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


namespace DRT
{
namespace ELEMENTS
{


/*--------------------------------------------------------------------------------
 *--------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void FluidEleCalcXFEM<distype>::ElementXfemInterface(
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

  std::map<int, Teuchos::RCP<DRT::ELEMENTS::XFLUID::SideInterface<distype> > > side_impl;
  Teuchos::RCP<DRT::ELEMENTS::XFLUID::SideInterface<distype> > si;

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

        si = DRT::ELEMENTS::XFLUID::SideInterface<distype>::Impl(side,C_uiu,C_uui,rhC_ui,eleGsui,eleGuis,side_xyze);
    }
    else
    {
        si = DRT::ELEMENTS::XFLUID::SideInterface<distype>::Impl(side,side_xyze);
    }

//    Teuchos::RCP<DRT::ELEMENTS::XFLUID::CouplingInterface> cp = DRT::ELEMENTS::XFLUID::CouplingInterface::Impl(ele,side,side_xyze);

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
    for (typename std::map<int, Teuchos::RCP<DRT::ELEMENTS::XFLUID::SideInterface<distype> > >::iterator i=side_impl.begin();  i!=side_impl.end(); ++i)
    {
      DRT::ELEMENTS::XFLUID::SideInterface<distype> * si = &*i->second;
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
template <DRT::Element::DiscretizationType distype>
void FluidEleCalcXFEM<distype>::ElementXfemInterfaceNitsche(
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

  std::map<int, Teuchos::RCP<DRT::ELEMENTS::XFLUID::SideInterface<distype> > > side_impl;
  Teuchos::RCP<DRT::ELEMENTS::XFLUID::SideInterface<distype> > si;

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

       si = DRT::ELEMENTS::XFLUID::SideInterface<distype>::Impl(side,C_uiu,C_uui,rhC_ui,eleCuiui,side_xyze);
     }
     else
     {
       si = DRT::ELEMENTS::XFLUID::SideInterface<distype>::Impl(side,side_xyze);
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
template <DRT::Element::DiscretizationType distype>
void FluidEleCalcXFEM<distype>::ElementXfemInterfaceNitscheTwoSided(
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

   std::map<int, Teuchos::RCP<DRT::ELEMENTS::XFLUID::EmbCoupling<distype> > > emb_impl;
   std::map<int, Teuchos::RCP<DRT::ELEMENTS::XFLUID::SideInterface<distype> > > side_impl;
   Teuchos::RCP<DRT::ELEMENTS::XFLUID::SideInterface<distype> > si;
   Teuchos::RCP<DRT::ELEMENTS::XFLUID::EmbCoupling<distype> > emb;

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

          emb = DRT::ELEMENTS::XFLUID::EmbCoupling<distype>::TwoSidedImpl(emb_ele,C_uiu,C_uui,rhC_ui,eleCuiui,emb_xyze);
          si  = DRT::ELEMENTS::XFLUID::SideInterface<distype>::Impl(side,side_xyze);
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
template <DRT::Element::DiscretizationType distype>
void FluidEleCalcXFEM<distype>::CalculateContinuityXFEM(
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
void FluidEleCalcXFEM<distype>::CalculateContinuityXFEM(
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




/*--------------------------------------------------------------------------------
 *--------------------------------------------------------------------------------*/
template <>
void FluidEleCalcXFEM<DRT::Element::tri3>::ElementXfemInterface(
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
void FluidEleCalcXFEM<DRT::Element::tri6>::ElementXfemInterface(
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
void FluidEleCalcXFEM<DRT::Element::quad4>::ElementXfemInterface(
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
void FluidEleCalcXFEM<DRT::Element::quad8>::ElementXfemInterface(
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
void FluidEleCalcXFEM<DRT::Element::quad9>::ElementXfemInterface(
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
void FluidEleCalcXFEM<DRT::Element::nurbs9>::ElementXfemInterface(
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
void FluidEleCalcXFEM<DRT::Element::tri3>::ElementXfemInterfaceNitsche(
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
void FluidEleCalcXFEM<DRT::Element::tri6>::ElementXfemInterfaceNitsche(
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
void FluidEleCalcXFEM<DRT::Element::quad4>::ElementXfemInterfaceNitsche(
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
void FluidEleCalcXFEM<DRT::Element::quad8>::ElementXfemInterfaceNitsche(
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
void FluidEleCalcXFEM<DRT::Element::quad9>::ElementXfemInterfaceNitsche(
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
void FluidEleCalcXFEM<DRT::Element::nurbs9>::ElementXfemInterfaceNitsche(
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
void FluidEleCalcXFEM<DRT::Element::tri3>::ElementXfemInterfaceNitscheTwoSided(
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
void FluidEleCalcXFEM<DRT::Element::tri6>::ElementXfemInterfaceNitscheTwoSided(
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
void FluidEleCalcXFEM<DRT::Element::quad4>::ElementXfemInterfaceNitscheTwoSided(
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
void FluidEleCalcXFEM<DRT::Element::quad8>::ElementXfemInterfaceNitscheTwoSided(
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
void FluidEleCalcXFEM<DRT::Element::quad9>::ElementXfemInterfaceNitscheTwoSided(
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
void FluidEleCalcXFEM<DRT::Element::nurbs9>::ElementXfemInterfaceNitscheTwoSided(
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


  } // end namespace ELEMENTS
} // end namespace DRT







// Ursula is responsible for this comment!
// all these templates objects have to be created when all element types are created in instance
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


