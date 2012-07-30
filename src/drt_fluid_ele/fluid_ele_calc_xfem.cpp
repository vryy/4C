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
#include "../drt_cut/cut_volumecell.H"

#include "../drt_geometry/position_array.H"

#include "../linalg/linalg_utils.H"

#include "fluid_ele_calc_stabilization.H"
#include "fluid_ele.H"
#include "fluid_ele_parameter.H"
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
 * add mixed/stress/hybrid (MSH) interface condition to element matrix and rhs
 *--------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void FluidEleCalcXFEM<distype>::ElementXfemInterfaceMSH(
    DRT::ELEMENTS::Fluid *                                              ele,               ///< fluid element
    DRT::Discretization &                                               dis,               ///< background discretization
    const std::vector<int> &                                            lm,                ///< element local map
    const std::vector<DRT::UTILS::GaussIntegration> &                   intpoints,         ///< background element integration points
    DRT::Discretization &                                               cutdis,            ///< cut discretization
    const std::map<int, std::vector<GEO::CUT::BoundaryCell*> > &        bcells,            ///< boundary cells
    const std::map<int, std::vector<DRT::UTILS::GaussIntegration> > &   bintpoints,        ///< boundary integration points
    std::map<int, std::vector<Epetra_SerialDenseMatrix> > &             side_coupling,     ///< side coupling matrices
    Teuchos::ParameterList&                                             params,            ///< parameter list
    Epetra_SerialDenseMatrix&                                           elemat1_epetra,    ///< element matrix
    Epetra_SerialDenseVector&                                           elevec1_epetra,    ///< element vector
    Epetra_SerialDenseMatrix&                                           Cuiui,             ///< ui-ui coupling matrix
    std::string&                                                        VCellGaussPts,     ///< Method of volumecell gauss point generation
    const GEO::CUT::plain_volumecell_set&                               cells              ///< Volumecells in the present set
  )
{

  const Teuchos::RCP<Epetra_Vector> iforcecol = params.get<Teuchos::RCP<Epetra_Vector> >("iforcenp", Teuchos::null);
  bool assemble_iforce = false;
  if(iforcecol != Teuchos::null) assemble_iforce = true;

  double nitsche_stab_conv = params.get<double>("conv_stab_fac");
  INPAR::XFEM::ConvStabScaling conv_stab_scaling = params.get<INPAR::XFEM::ConvStabScaling>("conv_stab_scaling");


  //----------------------------------------------------------------------------
  //                         ELEMENT GEOMETRY
  //----------------------------------------------------------------------------

  // get node coordinates
  GEO::fillInitialPositionArray< distype, my::nsd_, LINALG::Matrix<my::nsd_,my::nen_> >( ele, my::xyze_ );

  // get element-wise velocity/pressure field
  LINALG::Matrix<my::nsd_,my::nen_> evelaf(true);
  LINALG::Matrix<my::nen_,1> epreaf(true);
  my::ExtractValuesFromGlobalVector(dis, lm, *my::rotsymmpbc_, &evelaf, &epreaf, "velaf");

  int eid = ele->Id();

  //--------------------------------------------

  LINALG::Matrix<my::nen_,my::nen_> bK_ss;                // (N * N^T)
  LINALG::Matrix<my::nen_,my::nen_> invbK_ss( true );     // (N * N^T)^(-1)
  LINALG::Matrix<my::nen_,my::nen_> half_invbK_ss;        // 1/2 * (N * N^T)^(-1)

  // block matrices for couplings between stress-components and (ux,uy,uz,p)-components
  LINALG::BlockMatrix<LINALG::Matrix<my::nen_,my::nen_>,6,(my::nsd_+1)> K_su;
  LINALG::BlockMatrix<LINALG::Matrix<my::nen_,my::nen_>,(my::nsd_+1),6> K_us;
  LINALG::BlockMatrix<LINALG::Matrix<my::nen_,my::nen_>,6,6>            invK_ss;
  LINALG::BlockMatrix<LINALG::Matrix<my::nen_,   1>,6,1>                rhs;

  //--------------------------------------------

  const unsigned Velx = 0;
  const unsigned Vely = 1;
  const unsigned Velz = 2;
  //const unsigned Pres = 3;

  const unsigned Sigmaxx = 0;
  const unsigned Sigmaxy = 1;
  const unsigned Sigmaxz = 2;
  const unsigned Sigmayx = 1;
  const unsigned Sigmayy = 3;
  const unsigned Sigmayz = 4;
  const unsigned Sigmazx = 2;
  const unsigned Sigmazy = 4;
  const unsigned Sigmazz = 5;


  INPAR::XFEM::MSH_L2_Proj msh_l2_proj =  params.get<INPAR::XFEM::MSH_L2_Proj>("msh_l2_proj");

  MSH_Build_K_Matrices(msh_l2_proj, intpoints, VCellGaussPts, cells, eid, evelaf, epreaf, bK_ss, invbK_ss, K_su,rhs);


  //----------------------------------------------------------------------------
  //            surface integral --- build G_su, G_us matrices and rhs
  //----------------------------------------------------------------------------

  DRT::Element::LocationArray cutla( 1 );

  LINALG::Matrix<3,1> normal;
  LINALG::Matrix<3,1> x_side;

  bool fluidfluidcoupling = false;

  // side coupling implementation between background element and each cut side (map<sid, side_impl)
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

  // map of boundary element gids and coupling matrices, [0]: Gsui, [1]: Guis
  std::map<int, std::vector<Epetra_SerialDenseMatrix> > Cuiui_coupling;

  // lm vector of all intersecting boundary elements that intersect the current background element
  std::vector<int> patchelementslmv;
  std::vector<int> patchelementslmowner;

  // create location vectors for intersecting boundary elements and reshape coupling matrices
  PatchLocationVector(begids,cutdis,patchelementslmv,patchelementslmowner, Cuiui_coupling, "MixedStressHybrid");


  // coupling between domain and all sides (boundary elements) that cut the element
  Epetra_SerialDenseMatrix Gsui(my::nen_*6,patchelementslmv.size());
  Epetra_SerialDenseMatrix Guis(patchelementslmv.size(),my::nen_*6);
  Epetra_SerialDenseMatrix InvKss(my::nen_*6,my::nen_*6);
  Epetra_SerialDenseMatrix GuisInvKss(patchelementslmv.size(),my::nen_*6);


  // evaluate shape function derivatives
  bool eval_deriv = false;

  if(!fluidfluidcoupling) eval_deriv = true; // evaluate derivatives to evaluate traction vector



  //--------------------------------------------
  // loop intersecting sides
  //--------------------------------------------
  // map of side-element id and Gauss points
  for ( std::map<int, std::vector<DRT::UTILS::GaussIntegration> >::const_iterator i=bintpoints.begin();
        i!=bintpoints.end();
        ++i )
  {
    int sid = i->first;
    const std::vector<DRT::UTILS::GaussIntegration> & cutintpoints = i->second;

    // get side's boundary cells
    std::map<int, std::vector<GEO::CUT::BoundaryCell*> >::const_iterator j = bcells.find( sid );
    if ( j==bcells.end() )
      dserror( "missing boundary cell" );

    const std::vector<GEO::CUT::BoundaryCell*> & bcs = j->second;
    if ( bcs.size()!=cutintpoints.size() )
      dserror( "boundary cell integration rules mismatch" );

    // side and location vector
    DRT::Element * side = cutdis.gElement( sid );
    side->LocationVector(cutdis,cutla,false);

    // side geometry
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

    // create side impl
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

      si = DRT::ELEMENTS::XFLUID::SideInterface<distype>::Impl(side,C_uiu,C_uui,rhC_ui,eleGsui,eleGuis,side_xyze);
    }
    else
    {
      si = DRT::ELEMENTS::XFLUID::SideInterface<distype>::Impl(side,side_xyze);
    }


    side_impl[sid] = si;

    // get velocity at integration point of boundary dis
    si->eivel(cutdis,"ivelnp",cutla[0].lm_);

    // set displacement of side
    si->addeidisp(cutdis,"idispnp",cutla[0].lm_);

    // define interface force vector w.r.t side
    Epetra_SerialDenseVector iforce;
    iforce.Size(cutla[0].lm_.size());


    //--------------------------------------------
    // loop boundary cells w.r.t current cut side
    //--------------------------------------------
    for ( std::vector<DRT::UTILS::GaussIntegration>::const_iterator i=cutintpoints.begin();
          i!=cutintpoints.end();
          ++i )
    {
      const DRT::UTILS::GaussIntegration & gi = *i;
      GEO::CUT::BoundaryCell * bc = bcs[i - cutintpoints.begin()]; // get the corresponding boundary cell

      //--------------------------------------------
      // loop gausspoints w.r.t current boundary cell
      //--------------------------------------------
      for ( DRT::UTILS::GaussIntegration::iterator iquad=gi.begin(); iquad!=gi.end(); ++iquad )
      {
        double drs = 0.0; // transformation factor between reference cell and linearized boundary cell

        const LINALG::Matrix<2,1> eta( iquad.Point() ); // xi-coordinates with respect to side

        LINALG::Matrix<3,1> rst(true); // local coordinates w.r.t background element


#ifdef BOUNDARYCELL_TRANSFORMATION_OLD

        si->Evaluate(eta,x_side,normal,drs);

        // find element local position of gauss point at interface
        GEO::CUT::Position<distype> pos( my::xyze_, x_side );
        pos.Compute();
        rst = pos.LocalCoordinates();

#else
        LINALG::Matrix<3,1> x_gp_lin(true); // gp in xyz-system on linearized interface

        // compute transformation factor, normal vector and global Gauss point coordiantes
        if(bc->Shape() != DRT::Element::dis_none) // Tessellation approach
        {
          ComputeSurfaceTransformation(drs, x_gp_lin, normal, bc, eta);
        }
        else // MomentFitting approach
        {
          drs = 1.0;
          normal = bc->GetNormalVector();
          const double* gpcord = iquad.Point();
          for (int idim=0;idim<3;idim++)
          {
            x_gp_lin(idim,0) = gpcord[idim];
          }
        }

        // find element local position of gauss point
        GEO::CUT::Position<distype> pos( my::xyze_, x_gp_lin );
        pos.Compute();
        rst = pos.LocalCoordinates();

        // project gaussian point from linearized interface to warped side (get/set local side coordinates in SideImpl)
        LINALG::Matrix<2,1> xi_side(true);
        si->ProjectOnSide(x_gp_lin, x_side, xi_side);
 #endif

        const double surf_fac = drs*iquad.Weight();

        const double fac = surf_fac * my::fldpara_->TimeFac();


        //--------------------------------------------

        // evaluate shape functions (and derivatives)

        if(eval_deriv)
        {
          EvalFuncAndDeriv( rst );
        }
        else
        {
          DRT::UTILS::shape_function<distype>( rst, my::funct_ );
        }


        // get velocity at integration point
        // (values at n+alpha_F for generalized-alpha scheme, n+1 otherwise)
        my::velint_.Multiply(evelaf,my::funct_);


        //--------------------------------------------

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


        //--------------------------------------------
        // evaluate weak dirichlet boundary term or build coupling matrices

        if( fluidfluidcoupling )
        {
          si->MSH_buildCouplingMatrices(normal,fac,my::funct_,rhs);
        }
        else
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

#if(1)

        //--------------------------------------------
        // compute stabilization factors

        double stabfac_visc = 0.0;
        double stabfac_conv = 0.0;

        NIT_ComputeStabfac(fluidfluidcoupling,         // if fluidfluidcoupling
                           stabfac_visc,               // stabfac 1 for standard Nitsche term to set
                           stabfac_conv,               // stabfac 2 for additional stabilization term to set
                           0.0,                        // viscous Nitsche prefactor
                           conv_stab_scaling,          // type of scaling for convective stabilization term
                           nitsche_stab_conv,          // stabilization factor for additional stabilization
                           my::velint_.Dot(normal)     // velocity in normal direction
                           );

        //--------------------------------------------
        // evaluate additional inflow/convective stabilization terms


        if(fluidfluidcoupling)
        {

          bool bg_mortaring = true; // one-sided background fluid mortaring (kappa1=1, kappa2=0)

          //zero velocity jump for fluidfluidcoupling
          LINALG::Matrix<my::nsd_,1> ivelint_WDBC_JUMP(true);


          si->MSH_Stab_InflowCoercivity(
              elemat1_epetra,          // standard bg-bg-matrix
              elevec1_epetra,          // standard bg-rhs
              fluidfluidcoupling,      // assemble coupling terms (yes/no)
              bg_mortaring,            // yes: background-sided mortaring, no: coupling between two meshes (mixed mortaring)
              normal,                  // normal vector
              fac,                     // theta*dt
              my::visceff_,            // viscosity in background fluid
              0.0,                     // viscosity in embedded fluid
              0.5,                     // mortaring weighting
              0.5,                     // mortaring weighting
              stabfac_conv,            // Nitsche convective non-dimensionless stabilization factor
              my::funct_,              // bg shape functions
              my::derxy_,              // bg deriv
              my::vderxy_,             // bg deriv^n
              my::velint_,             // bg u^n
              ivelint_WDBC_JUMP        // Dirichlet velocity vector or prescribed jump vector
          );


        }
        if(!fluidfluidcoupling)
        {
          // case for one-sided weak Dirichlet
          bool bg_mortaring = true; // one-sided background fluid mortaring (kappa1=1, kappa2=0)

          // prescribed velocity vector at weak Dirichlet boundary
          LINALG::Matrix<my::nsd_,1> ivelint_WDBC_JUMP(true);
          si->get_vel_WeakDBC(ivelint_WDBC_JUMP);


          si->MSH_Stab_InflowCoercivity(
              elemat1_epetra,          // standard bg-bg-matrix
              elevec1_epetra,          // standard bg-rhs
              fluidfluidcoupling,      // assemble coupling terms (yes/no)
              bg_mortaring,            // yes: background-sided mortaring, no: coupling between two meshes (mixed mortaring)
              normal,                  // normal vector
              fac,                     // theta*dt
              my::visceff_,            // viscosity in background fluid
              0.0,                     // viscosity in embedded fluid
              1.0,                     // mortaring weighting
              0.0,                     // mortaring weighting
              stabfac_visc,            // Nitsche convective non-dimensionless stabilization factor
              my::funct_,              // bg shape functions
              my::derxy_,              // bg deriv
              my::vderxy_,             // bg deriv^n
              my::velint_,             // bg u^n
              ivelint_WDBC_JUMP        // Dirichlet velocity vector or prescribed jump vector
          );
        }
#endif

        //--------------------------------------------
        // calculate interface forces for XFSI
        if(assemble_iforce)
        {

          // get velocity derivatives at integration point
          // (values at n+alpha_F for generalized-alpha scheme, n+1 otherwise)
          my::vderxy_.MultiplyNT(evelaf,my::derxy_);

          // get pressure at integration point
          // (value at n+alpha_F for generalized-alpha scheme, n+1 otherwise)
          double press = my::funct_.Dot(epreaf);

          //-------------------------------
          LINALG::Matrix<my::nsd_,1> traction(true);

          buildTractionVector( traction, press, normal );

          si->InterfaceForce(iforce, traction, surf_fac );
        } // buildInterfaceForce


      } // end loop gauss points of boundary cell
    } // end loop boundary cells of side

    if(assemble_iforce) AssembleInterfaceForce(iforcecol, cutdis, cutla[0].lm_, iforce);

  } // end loop cut sides


  // construct views
  LINALG::Matrix<(my::nsd_+1)*my::nen_,(my::nsd_+1)*my::nen_> elemat1(elemat1_epetra,true);
  LINALG::Matrix<(my::nsd_+1)*my::nen_,            1> elevec1(elevec1_epetra,true);



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

  LINALG::BlockMatrix<LINALG::Matrix<my::nen_,my::nen_>,(my::nsd_+1),6>            K_iK;
  LINALG::BlockMatrix<LINALG::Matrix<my::nen_,my::nen_>,(my::nsd_+1),(my::nsd_+1)> extK;
  LINALG::BlockMatrix<LINALG::Matrix<my::nen_,   1>,(my::nsd_+1),1>                extrhs;

  K_iK  .Multiply( K_us, invK_ss );
  extK  .Multiply( K_iK, K_su );
  extrhs.Multiply( K_iK, rhs );


  //--------------------------------------------
  // fill element matrix
  //--------------------------------------------
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

  //--------------------------------------------
  // fill element vector
  //--------------------------------------------
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


  //----------------------------------------------------------------------------
  //           build (G_sui + K_sui), (G_uis + K_uis) coupling matrices
  //----------------------------------------------------------------------------

  if ( fluidfluidcoupling )
  {
    // build fluid-fluid matrices
    for (typename std::map<int, Teuchos::RCP<DRT::ELEMENTS::XFLUID::SideInterface<distype> > >::iterator i=side_impl.begin();  i!=side_impl.end(); ++i)
    {
      DRT::ELEMENTS::XFLUID::SideInterface<distype> * si = &*i->second;
      si->MSH_buildFinalCouplingMatrices(invK_ss,K_iK,K_su,rhs);
    }
    //--------------------------------------------
    // build InvKss ( K_ss^(-1) )
    //--------------------------------------------
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

    //--------------------------------------------
    // add G_sui to K_sui => (G_sui+K_sui) and
    // add G_uis to K_uis => (G_uis+K_uis)
    //--------------------------------------------
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

    //----------------------------------------------------------------------------
    //  build condensed   C_uiui=[ (G_sui + K_sui) * K_ss^(-1) * (G_uis + K_uis) ]
    //----------------------------------------------------------------------------

    GuisInvKss.Multiply('N','N',1.0,Guis,InvKss,1.0);
    Cuiui.Multiply('N','N',1.0,GuisInvKss,Gsui,1.0);

  } // end build coupling matrices



//   std::cout << elemat1;
//   std::cout << elevec1;

//   elemat1_epetra.Print( std::cout );

  return;

} // ElementXfemInterface



/*--------------------------------------------------------------------------------
 * build volume-based K matrices for coupling stress fields
 *-------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void FluidEleCalcXFEM<distype>::MSH_Build_K_Matrices(
    INPAR::XFEM::MSH_L2_Proj                                               msh_l2_proj,   ///< full or partial l2 projection for MSH method
    const std::vector<DRT::UTILS::GaussIntegration> &                      intpoints,     ///< background element integration points
    std::string &                                                          VCellGaussPts, ///< volumecell gaussian points method
    const GEO::CUT::plain_volumecell_set&                                  cells,         ///< volumecells in the present set
    int                                                                    eid,           ///< element ID
    LINALG::Matrix<my::nsd_,my::nen_>&                                     evelaf,        ///< element velocity
    LINALG::Matrix<my::nen_,1>&                                            epreaf,        ///< element pressure
    LINALG::Matrix<my::nen_,my::nen_>&                                     bK_ss,         ///< block K_ss matrix
    LINALG::Matrix<my::nen_,my::nen_>&                                     invbK_ss,      ///< inverse of block K_ss matrix
    LINALG::BlockMatrix<LINALG::Matrix<my::nen_,my::nen_>,6,(my::nsd_+1)>& K_su,          ///< K_su matrix
    LINALG::BlockMatrix<LINALG::Matrix<my::nen_,   1>,6,1>&                rhs            ///< rhs vector
    )
{

  //----------------------------------------------------------------------------
  //            volume integral --- build K_su, K_sp matrices and rhs
  //----------------------------------------------------------------------------
  if( msh_l2_proj == INPAR::XFEM::MSH_L2_Proj_full ) // projection over the full background element
  {
    for ( DRT::UTILS::GaussIntegration::iterator iquad=my::intpoints_.begin(); iquad!=my::intpoints_.end(); ++iquad )
    {
      // evaluate shape functions and derivatives at integration point
      my::EvalShapeFuncAndDerivsAtIntPoint(iquad,eid);
      MSH_EvaluateMatrices(evelaf,epreaf,bK_ss,invbK_ss,K_su,rhs);
    }
  }
  else //projection over the volumecells (partial projection)
  {
    if( VCellGaussPts!="DirectDivergence" ) // standard case for tessellation and momentfitting
    {
      for( std::vector<DRT::UTILS::GaussIntegration>::const_iterator i=intpoints.begin();i!=intpoints.end();++i )
      {
        const DRT::UTILS::GaussIntegration intcell = *i;
        for ( DRT::UTILS::GaussIntegration::iterator iquad=intcell.begin(); iquad!=intcell.end(); ++iquad )
        {
          // evaluate shape functions and derivatives at integration point
          my::EvalShapeFuncAndDerivsAtIntPoint(iquad,eid);
          MSH_EvaluateMatrices(evelaf,epreaf,bK_ss,invbK_ss,K_su,rhs);
        }
      }
    }
    else  // DirectDivergence method
    {
      for( std::vector<DRT::UTILS::GaussIntegration>::const_iterator i=intpoints.begin();i!=intpoints.end();++i )
      {
        const DRT::UTILS::GaussIntegration intcell = *i;
        GEO::CUT::VolumeCell * vc = cells[i-intpoints.begin()];
        //----------------------------------------------------------------------
        //integration over the main gauss points to get the required integral
        //----------------------------------------------------------------------
        int mainPtno = 0;
        for ( DRT::UTILS::GaussIntegration::iterator iquad=intcell.begin(); iquad!=intcell.end(); ++iquad )
        {
          LINALG::Matrix<my::nen_,my::nen_> invbK_ssTemp( true );
          LINALG::BlockMatrix<LINALG::Matrix<my::nen_,my::nen_>,6,(my::nsd_+1)> K_suTemp;
          LINALG::BlockMatrix<LINALG::Matrix<my::nen_,   1>,6,1>                rhsTemp;

          // get internal Gaussian rule for every main Gauss point
          DRT::UTILS::GaussIntegration gint = vc->GetInternalRule( mainPtno );
          mainPtno++;

          //----------------------------------------------------------------------
          //integration over the internal gauss points - to get modified integrand
          //----------------------------------------------------------------------
          for ( DRT::UTILS::GaussIntegration::iterator quadint=gint.begin(); quadint!=gint.end(); ++quadint )
          {
            my::EvalShapeFuncAndDerivsAtIntPoint( quadint, eid );
            MSH_EvaluateMatrices( evelaf, epreaf, bK_ss, invbK_ssTemp, K_suTemp, rhsTemp );
          }

          my::EvalShapeFuncAndDerivsAtIntPoint( iquad, eid );
          bK_ss.MultiplyNT( my::funct_, my::funct_ );

          invbK_ss.Update( iquad.Weight(), invbK_ssTemp, 1.0 );

          const unsigned Velx = 0;
          const unsigned Vely = 1;
          const unsigned Velz = 2;
          const unsigned Pres = 3;

          const unsigned Sigmaxx = 0;
          const unsigned Sigmaxy = 1;
          const unsigned Sigmaxz = 2;
          const unsigned Sigmayx = 1;
          const unsigned Sigmayy = 3;
          const unsigned Sigmayz = 4;
          const unsigned Sigmazx = 2;
          const unsigned Sigmazy = 4;
          const unsigned Sigmazz = 5;

          K_su( Sigmaxx, Velx )->Update( iquad.Weight(), *K_suTemp( Sigmaxx, Velx ), 1.0 );
          K_su( Sigmaxy, Velx )->Update( iquad.Weight(), *K_suTemp( Sigmaxy, Velx ), 1.0 );
          K_su( Sigmayx, Vely )->Update( iquad.Weight(), *K_suTemp( Sigmayx, Vely ), 1.0 );
          K_su( Sigmaxz, Velx )->Update( iquad.Weight(), *K_suTemp( Sigmaxz, Velx ), 1.0 );
          K_su( Sigmazx, Velz )->Update( iquad.Weight(), *K_suTemp( Sigmazx, Velz ), 1.0 );
          K_su( Sigmayy, Vely )->Update( iquad.Weight(), *K_suTemp( Sigmayy, Vely ), 1.0 );
          K_su( Sigmayz, Vely )->Update( iquad.Weight(), *K_suTemp( Sigmayz, Vely ), 1.0 );
          K_su( Sigmazy, Velz )->Update( iquad.Weight(), *K_suTemp( Sigmazy, Velz ), 1.0 );
          K_su( Sigmazz, Velz )->Update( iquad.Weight(), *K_suTemp( Sigmazz, Velz ), 1.0 );

          rhs( Sigmaxx, 0 )->Update( iquad.Weight(), *rhsTemp( Sigmaxx, 0 ), 1.0 );
          rhs( Sigmaxy, 0 )->Update( iquad.Weight(), *rhsTemp( Sigmaxy, 0 ), 1.0 );
          rhs( Sigmaxz, 0 )->Update( iquad.Weight(), *rhsTemp( Sigmaxz, 0 ), 1.0 );
          rhs( Sigmayy, 0 )->Update( iquad.Weight(), *rhsTemp( Sigmayy, 0 ), 1.0 );
          rhs( Sigmayz, 0 )->Update( iquad.Weight(), *rhsTemp( Sigmayz, 0 ), 1.0 );
          rhs( Sigmazz, 0 )->Update( iquad.Weight(), *rhsTemp( Sigmazz, 0 ), 1.0 );

          K_su( Sigmaxx, Pres )->Update( iquad.Weight(), *K_suTemp( Sigmaxx, Pres ), 1.0 );
          K_su( Sigmayy, Pres )->Update( iquad.Weight(), *K_suTemp( Sigmayy, Pres ), 1.0 );
          K_su( Sigmazz, Pres )->Update( iquad.Weight(), *K_suTemp( Sigmazz, Pres ), 1.0 );
        }
      }
    }
  }

  return;

} // Build_K_Matrices_MSH


/*--------------------------------------------------------------------------------
 * evaluate volume-based K matrix terms at Gaussian point
 *-------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void FluidEleCalcXFEM<distype>::MSH_EvaluateMatrices(
    LINALG::Matrix<my::nsd_,my::nen_>&    evelaf,                                ///< element velocity
    LINALG::Matrix<my::nen_,1>&           epreaf,                                ///< element pressure
    LINALG::Matrix<my::nen_,my::nen_>&    bK_ss,                                 ///< block K_ss matrix
    LINALG::Matrix<my::nen_,my::nen_>&    invbK_ss,                              ///< inverse of block K_ss matrix
    LINALG::BlockMatrix<LINALG::Matrix<my::nen_,my::nen_>,6,(my::nsd_+1)>& K_su, ///< K_su matrix
    LINALG::BlockMatrix<LINALG::Matrix<my::nen_,   1>,6,1>&                rhs   ///< rhs vector
    )
{
  const unsigned Velx = 0;
  const unsigned Vely = 1;
  const unsigned Velz = 2;
  const unsigned Pres = 3;

  const unsigned Sigmaxx = 0;
  const unsigned Sigmaxy = 1;
  const unsigned Sigmaxz = 2;
  const unsigned Sigmayx = 1;
  const unsigned Sigmayy = 3;
  const unsigned Sigmayz = 4;
  const unsigned Sigmazx = 2;
  const unsigned Sigmazy = 4;
  const unsigned Sigmazz = 5;

  LINALG::Matrix<my::nen_,1> dx;
  LINALG::Matrix<my::nen_,1> dy;
  LINALG::Matrix<my::nen_,1> dz;

  LINALG::Matrix<my::nen_,my::nen_> conv_x;
  LINALG::Matrix<my::nen_,my::nen_> conv_y;
  LINALG::Matrix<my::nen_,my::nen_> conv_z;

  //----------------------------------------------------------------------
  // set time-integration factors for left- and right-hand side
  // (two right-hand-side factors: general and for residuals)
  //----------------------------------------------------------------------

  // TODO: check if this parameter should be fac*1.0*dt (full implicit stabilization)
  const double timefacfac = my::fldpara_->TimeFac() * my::fac_;

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


  //--------------------------------------------

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

   rhs( Sigmaxx, 0 )->Update( - timefacfac* my::vderxy_(0, 0)                     , my::funct_, 1.0 );
   rhs( Sigmaxy, 0 )->Update( - timefacfac*(my::vderxy_(0, 1) + my::vderxy_(1, 0)), my::funct_, 1.0 );
   rhs( Sigmaxz, 0 )->Update( - timefacfac*(my::vderxy_(0, 2) + my::vderxy_(2, 0)), my::funct_, 1.0 );
   rhs( Sigmayy, 0 )->Update( - timefacfac* my::vderxy_(1, 1)                     , my::funct_, 1.0 );
   rhs( Sigmayz, 0 )->Update( - timefacfac*(my::vderxy_(1, 2) + my::vderxy_(2, 1)), my::funct_, 1.0 );
   rhs( Sigmazz, 0 )->Update( - timefacfac* my::vderxy_(2, 2)                     , my::funct_, 1.0 );

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

  return;

} //EvaluateMatricesMSH

/*--------------------------------------------------------------------------------
 * add Nitsche (NIT) interface condition to element matrix and rhs
 *--------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void FluidEleCalcXFEM<distype>::ElementXfemInterfaceNIT(
    DRT::ELEMENTS::Fluid *                                              ele,               ///< fluid element
    DRT::Discretization &                                               dis,               ///< background discretization
    const std::vector<int> &                                            lm,                ///< element local map
    DRT::Discretization &                                               cutdis,            ///< cut discretization
    const std::map<int, std::vector<GEO::CUT::BoundaryCell*> > &        bcells,            ///< boundary cells
    const std::map<int, std::vector<DRT::UTILS::GaussIntegration> > &   bintpoints,        ///< boundary integration points
    std::map<int, std::vector<Epetra_SerialDenseMatrix> > &             side_coupling,     ///< side coupling matrices
    Teuchos::ParameterList&                                             params,            ///< parameter list
    Epetra_SerialDenseMatrix&                                           elemat1_epetra,    ///< element matrix
    Epetra_SerialDenseVector&                                           elevec1_epetra,    ///< element vector
    Epetra_SerialDenseMatrix&                                           Cuiui,             ///< ui-ui coupling matrix
    const GEO::CUT::plain_volumecell_set&                               vcSet							 ///< volumecell sets in this element
  )
{

  bool compute_meas_surf = false;
  bool compute_meas_vol  = false;


  const Teuchos::RCP<Epetra_Vector> iforcecol = params.get<Teuchos::RCP<Epetra_Vector> >("iforcenp", Teuchos::null);

  bool assemble_iforce = false;
  if(iforcecol != Teuchos::null) assemble_iforce = true;

  double nitsche_stab      = params.get<double>("visc_stab_fac");
  double nitsche_stab_conv = params.get<double>("conv_stab_fac");
  INPAR::XFEM::ViscStabScaling visc_stab_scaling = params.get<INPAR::XFEM::ViscStabScaling>("visc_stab_scaling");
  INPAR::XFEM::ConvStabScaling conv_stab_scaling = params.get<INPAR::XFEM::ConvStabScaling>("conv_stab_scaling");
  INPAR::XFEM::ViscStab_hk visc_stab_hk = params.get<INPAR::XFEM::ViscStab_hk>("visc_stab_hk");

  if(visc_stab_hk == INPAR::XFEM::ViscStab_hk_vol_div_by_surf)
  {
    compute_meas_surf = true;
    compute_meas_vol  = true;
  }

  //----------------------------------------------------------------------------
  //     pre-evaluate fluid volume measure and intersecting surface measure
  //----------------------------------------------------------------------------
  double meas_partial_volume = 0.0;
  double meas_surf = 0.0;

  if(compute_meas_vol)
  {
    meas_partial_volume = 0.0;
    for( GEO::CUT::plain_volumecell_set::const_iterator i=vcSet.begin();i!=vcSet.end();i++ )
    {
      GEO::CUT::VolumeCell* vc = *i;
      meas_partial_volume += vc->Volume();
    }
  }

  if(compute_meas_surf)
  {
    meas_surf = ComputeMeasSurf(bintpoints, bcells);
  }

  //----------------------------------------------------------------------------
  //                         ELEMENT GEOMETRY
  //----------------------------------------------------------------------------

  // get node coordinates
  GEO::fillInitialPositionArray< distype, my::nsd_, LINALG::Matrix<my::nsd_,my::nen_> >( ele, my::xyze_ );

  // get element-wise velocity/pressure field
  LINALG::Matrix<my::nsd_,my::nen_> evelaf(true);
  LINALG::Matrix<my::nen_,1> epreaf(true);
  my::ExtractValuesFromGlobalVector(dis, lm, *my::rotsymmpbc_, &evelaf, &epreaf, "velaf");


  //----------------------------------------------------------------------------
  //      surface integral --- build Cuiui, Cuui, Cuiu and Cuu matrix and rhs
  //----------------------------------------------------------------------------

  DRT::Element::LocationArray cutla( 1 );

  LINALG::Matrix<3,1> normal;
  LINALG::Matrix<3,1> x_side;

  bool fluidfluidcoupling = false;

  // side coupling implementation between background element and each cut side (map<sid, side_impl)
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

  // map of boundary element gids and coupling matrices, [0]: Cuiui matrix
  std::map<int, std::vector<Epetra_SerialDenseMatrix> > Cuiui_coupling;

  // lm vector of all intersecting boundary elements that intersect the current background element
  std::vector<int> patchelementslmv;
  std::vector<int> patchelementslmowner;

  // create location vectors for intersecting boundary elements and reshape coupling matrices
  PatchLocationVector(begids,cutdis,patchelementslmv,patchelementslmowner, Cuiui_coupling, "Nitsche");


  //-----------------------------------------------------------------------------------
  //         evaluate element length, stabilization factors and average weights
  //-----------------------------------------------------------------------------------

  // element length
  double h_k = 0.0;

  if(visc_stab_hk == INPAR::XFEM::ViscStab_hk_vol_equivalent)
  {
    LINALG::Matrix<1,1> dummy(true);
    h_k = FLD::UTILS::HK<distype>(dummy, my::xyze_);
  }
  else if(visc_stab_hk == INPAR::XFEM::ViscStab_hk_vol_div_by_surf)
  {
    h_k = meas_partial_volume / meas_surf;
  }
  else if(visc_stab_hk == INPAR::XFEM::ViscStab_hk_longest_ele_length)
  {
    dserror("longest element length for hk not supported yet");
  }
  else dserror("unknown type of characteristic element length");

  if( h_k <= 0.0 ) dserror("element length is <= 0.0");


  //------------------------------
  // scaling factors for Nitsche's standard stabilization term
  double NIT_stab_fac         = 1.0;

  if(visc_stab_scaling == INPAR::XFEM::ViscStabScaling_visc_div_by_hk)
  {
    NIT_stab_fac = nitsche_stab * my::visceff_ / h_k;
  }
  else if(visc_stab_scaling == INPAR::XFEM::ViscStabScaling_inv_hk)
  {
    NIT_stab_fac = nitsche_stab / h_k;
  }
  else if(visc_stab_scaling == INPAR::XFEM::ViscStabScaling_const)
  {
    NIT_stab_fac = nitsche_stab;
  }
  else dserror("unknown scaling for viscous stabilization term");



  //------------------------------
  // define average weights
  double kappa1 = 1.0;      // Xfluid-sided mortaring

  if(meas_partial_volume < 0.0) dserror(" measure of cut partial volume is smaller than 0.0: %f Attention with increasing Nitsche-Parameter!!!", meas_partial_volume);

  if( kappa1 > 1.0 || kappa1 < 0.0) dserror("Nitsche weights for inverse estimate kappa1 lies not in [0,1]: %d", kappa1);

  double kappa2 = 1.0-kappa1;


  
  // evaluate shape function derivatives
  bool eval_deriv = true;


  //--------------------------------------------
  // loop intersecting sides
  //--------------------------------------------
  // map of side-element id and Gauss points
  for ( std::map<int, std::vector<DRT::UTILS::GaussIntegration> >::const_iterator i=bintpoints.begin();
        i!=bintpoints.end();
        ++i )
  {
    int sid = i->first;
    const std::vector<DRT::UTILS::GaussIntegration> & cutintpoints = i->second;

    // get side's boundary cells
    std::map<int, std::vector<GEO::CUT::BoundaryCell*> >::const_iterator j = bcells.find( sid );
    if ( j==bcells.end() )
      dserror( "missing boundary cell" );

    const std::vector<GEO::CUT::BoundaryCell*> & bcs = j->second;
    if ( bcs.size()!=cutintpoints.size() )
      dserror( "boundary cell integration rules mismatch" );

    // side and location vector
    DRT::Element * side = cutdis.gElement( sid );
    side->LocationVector(cutdis,cutla,false);

    // side geometry
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

    // create side impl
    if(fluidfluidcoupling)
    {
      // coupling matrices between background element and one! side
      Epetra_SerialDenseMatrix & C_uiu  = side_matrices[0];
      Epetra_SerialDenseMatrix & C_uui  = side_matrices[1];
      Epetra_SerialDenseMatrix & rhC_ui = side_matrices[2];

      // coupling matrices between one side and itself
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

    // set displacement of side
    si->addeidisp(cutdis,"idispnp",cutla[0].lm_);

    // define interface force vector w.r.t side
    Epetra_SerialDenseVector iforce;
    iforce.Size(cutla[0].lm_.size());


    //--------------------------------------------
    // loop boundary cells w.r.t current cut side
    //--------------------------------------------
    for ( std::vector<DRT::UTILS::GaussIntegration>::const_iterator i=cutintpoints.begin();
          i!=cutintpoints.end();
          ++i )
    {
      const DRT::UTILS::GaussIntegration & gi = *i;
      GEO::CUT::BoundaryCell * bc = bcs[i - cutintpoints.begin()]; // get the corresponding boundary cell

      //--------------------------------------------
      // loop gausspoints w.r.t current boundary cell
      //--------------------------------------------
      for ( DRT::UTILS::GaussIntegration::iterator iquad=gi.begin(); iquad!=gi.end(); ++iquad )
      {
        double drs = 0.0; // transformation factor between reference cell and linearized boundary cell

        const LINALG::Matrix<2,1> eta( iquad.Point() ); // xi-coordinates with respect to side

        LINALG::Matrix<3,1> rst(true); // local coordinates w.r.t background element

#ifdef BOUNDARYCELL_TRANSFORMATION_OLD

        si->Evaluate(eta,x_side,normal,drs);

        // find element local position of gauss point at interface
        GEO::CUT::Position<distype> pos( my::xyze_, x_side );
        pos.Compute();
        rst = pos.LocalCoordinates();

#else
        LINALG::Matrix<3,1> x_gp_lin(true); // gp in xyz-system on linearized interface

        // compute transformation factor, normal vector and global Gauss point coordiantes
        if(bc->Shape() != DRT::Element::dis_none) // Tessellation approach
        {
          ComputeSurfaceTransformation(drs, x_gp_lin, normal, bc, eta);
        }
        else // MomentFitting approach
        {
          drs = 1.0;
          normal = bc->GetNormalVector();
          const double* gpcord = iquad.Point();
          for (int idim=0;idim<3;idim++)
          {
            x_gp_lin(idim,0) = gpcord[idim];
          }
        }

        // find element local position of gauss point
        GEO::CUT::Position<distype> pos( my::xyze_, x_gp_lin );
        pos.Compute();
        rst = pos.LocalCoordinates();

        // project gaussian point from linearized interface to warped side (get/set local side coordinates in SideImpl)
        LINALG::Matrix<2,1> xi_side(true);
        si->ProjectOnSide(x_gp_lin, x_side, xi_side);
 #endif

        const double surf_fac = drs*iquad.Weight();

        const double timefacfac = surf_fac * my::fldpara_->TimeFac();


        //--------------------------------------------

        // evaluate shape functions (and derivatives)

        if(eval_deriv)
        {
          EvalFuncAndDeriv( rst );
        }
        else
        {
          DRT::UTILS::shape_function<distype>( rst, my::funct_ );
        }


        // get velocity at integration point
        // (values at n+alpha_F for generalized-alpha scheme, n+1 otherwise)
        my::velint_.Multiply(evelaf,my::funct_);

        // get velocity derivatives at integration point
        // (values at n+alpha_F for generalized-alpha scheme, n+1 otherwise)
        my::vderxy_.MultiplyNT(evelaf,my::derxy_);

        // get pressure at integration point
        // (value at n+alpha_F for generalized-alpha scheme, n+1 otherwise)
        double press = my::funct_.Dot(epreaf);

        //--------------------------------------------
        // compute stabilization factors

        double stabfac_visc = 0.0;
        double stabfac_conv = 0.0;

        NIT_ComputeStabfac(fluidfluidcoupling,         // if fluidfluidcoupling
                           stabfac_visc,               // stabfac 1 for standard Nitsche term to set
                           stabfac_conv,               // stabfac 2 for additional stabilization term to set
                           NIT_stab_fac,               // viscous Nitsche prefactor
                           conv_stab_scaling,          // type of scaling for convective stabilization term
                           nitsche_stab_conv,          // stabilization factor for additional stabilization
                           my::velint_.Dot(normal)     // velocity in normal direction
                           );


        //--------------------------------------------

        if(fluidfluidcoupling)
        {
          bool bg_mortaring = true; // one-sided background fluid mortaring (kappa1=1, kappa2=0)

          //zero velocity jump for fluidfluidcoupling
          LINALG::Matrix<my::nsd_,1> ivelint_WDBC_JUMP(true);


          si->NIT_buildCouplingMatrices(
              elemat1_epetra,          // standard bg-bg-matrix
              elevec1_epetra,          // standard bg-rhs
              fluidfluidcoupling,      // assemble coupling terms (yes/no)
              bg_mortaring,            // yes: background-sided mortaring, no: coupling between two meshes (mixed mortaring)
              normal,                  // normal vector
              timefacfac,              // theta*dt
              my::visceff_,            // viscosity in background fluid
              my::visceff_,            // viscosity in embedded fluid
              kappa1,                  // mortaring weighting
              kappa2,                  // mortaring weighting
              stabfac_visc,            // Nitsche non-dimensionless stabilization factor
              stabfac_conv,            // Nitsche convective non-dimensionless stabilization factor
              my::funct_,              // bg shape functions
              my::derxy_,              // bg deriv
              my::vderxy_,             // bg deriv^n
              press,                   // bg p^n
              my::velint_,             // bg u^n
              ivelint_WDBC_JUMP        // Dirichlet velocity vector or prescribed jump vector
          );



        }
        if(!fluidfluidcoupling)
        {
          // case for one-sided weak Dirichlet
          bool bg_mortaring = true; // one-sided background fluid mortaring (kappa1=1, kappa2=0)

          // prescribed velocity vector at weak Dirichlet boundary
          LINALG::Matrix<my::nsd_,1> ivelint_WDBC_JUMP(true);
          si->get_vel_WeakDBC(ivelint_WDBC_JUMP);

          si->NIT_buildCouplingMatrices(
              elemat1_epetra,          // standard bg-bg-matrix
              elevec1_epetra,          // standard bg-rhs
              fluidfluidcoupling,      // assemble coupling terms (yes/no)
              bg_mortaring,            // yes: background-sided mortaring, no: coupling between two meshes (mixed mortaring)
              normal,                  // normal vector
              timefacfac,              // theta*dt
              my::visceff_,            // viscosity in background fluid
              0.0,                     // viscosity in embedded fluid
              kappa1,                  // mortaring weighting
              kappa2,                  // mortaring weighting
              stabfac_visc,            // Nitsche non-dimensionless stabilization factor
              stabfac_conv,            // Nitsche convective non-dimensionless stabilization factor
              my::funct_,              // bg shape functions
              my::derxy_,              // bg deriv
              my::vderxy_,             // bg deriv^n
              press,                   // bg p^n
              my::velint_,             // bg u^n
              ivelint_WDBC_JUMP        // Dirichlet velocity vector or prescribed jump vector
          );
        }


        //--------------------------------------------
        // calculate interface forces for XFSI
        if(assemble_iforce)
        {

          // get pressure at integration point
          // (value at n+alpha_F for generalized-alpha scheme, n+1 otherwise)
          double press = my::funct_.Dot(epreaf);

          //-------------------------------
          LINALG::Matrix<my::nsd_,1> traction(true);

          buildTractionVector( traction, press, normal );

          si->InterfaceForce(iforce, traction, surf_fac );

        } // buildInterfaceForce

      } // end loop gauss points of boundary cell
    } // end loop boundary cells of side

    if(assemble_iforce) AssembleInterfaceForce(iforcecol, cutdis, cutla[0].lm_, iforce);

  } // end loop cut sides


  //----------------------------------------------------------------------------
  // build Cuiui coupling matrix (includes patch of Cuiui matrices for all sides)
  //----------------------------------------------------------------------------

  if ( fluidfluidcoupling )
    NIT_BuildPatchCuiui(Cuiui, Cuiui_coupling);

  return;
}



/*--------------------------------------------------------------------------------
 *--------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void FluidEleCalcXFEM<distype>::ElementXfemInterfaceNIT2(
    DRT::ELEMENTS::Fluid *                                              ele,
    DRT::Discretization &                                               dis,
    const std::vector<int> &                                            lm,
    DRT::Discretization &                                               cutdis,
    const std::map<int, std::vector<GEO::CUT::BoundaryCell*> > &        bcells,
    const std::map<int, std::vector<DRT::UTILS::GaussIntegration> > &   bintpoints,
    std::map<int, std::vector<Epetra_SerialDenseMatrix> > &             side_coupling,
    Teuchos::ParameterList&                                             params,
    DRT::Discretization &                                               alediscret,
    map<int,int> &                                                      boundary_emb_gid_map,
    Epetra_SerialDenseMatrix&                                           elemat1_epetra,
    Epetra_SerialDenseVector&                                           elevec1_epetra,
    Epetra_SerialDenseMatrix&                                           Cuiui,
    const GEO::CUT::plain_volumecell_set &                              vcSet
  )
{

  bool compute_meas_surf = false;
  bool compute_meas_vol  = false;

  INPAR::XFEM::CouplingStrategy coupling_strategy = params.get<INPAR::XFEM::CouplingStrategy>("coupling_strategy");
  if(coupling_strategy == INPAR::XFEM::Two_Sided_Mortaring) compute_meas_surf = true;

  double nitsche_stab      = params.get<double>("visc_stab_fac");
  double nitsche_stab_conv = params.get<double>("conv_stab_fac");
//  INPAR::XFEM::ViscStabScaling visc_stab_scaling = params.get<INPAR::XFEM::ViscStabScaling>("visc_stab_scaling");
//  INPAR::XFEM::ConvStabScaling conv_stab_scaling = params.get<INPAR::XFEM::ConvStabScaling>("conv_stab_scaling");
  INPAR::XFEM::ViscStab_hk visc_stab_hk = params.get<INPAR::XFEM::ViscStab_hk>("visc_stab_hk");

  if(visc_stab_hk == INPAR::XFEM::ViscStab_hk_vol_div_by_surf)
  {
    compute_meas_surf = true;
    compute_meas_vol  = true;
  }

  //----------------------------------------------------------------------------
  //     pre-evaluate fluid volume measure and intersecting surface measure
  //----------------------------------------------------------------------------
  double meas_partial_volume = 0.0;
  double meas_surf = 0.0;

  // for two-sided mortaring the stabilization parameter depends on interface/volume fraction
  if(coupling_strategy == INPAR::XFEM::Two_Sided_Mortaring)
  {
    if(compute_meas_vol)
    {
      meas_partial_volume = 0.0;
      for( GEO::CUT::plain_volumecell_set::const_iterator i=vcSet.begin();i!=vcSet.end();i++ )
      {
        GEO::CUT::VolumeCell* vc = *i;
        meas_partial_volume += vc->Volume();
      }
    }
  }

  if(compute_meas_surf)
  {
    meas_surf = ComputeMeasSurf(bintpoints, bcells);
  }

  //----------------------------------------------------------------------------
  //                         ELEMENT GEOMETRY
  //----------------------------------------------------------------------------

  // get node coordinates
  GEO::fillInitialPositionArray< distype, my::nsd_, LINALG::Matrix<my::nsd_,my::nen_> >( ele, my::xyze_ );

  // get element-wise velocity/pressure field
  LINALG::Matrix<my::nsd_,my::nen_> evelaf(true);
  LINALG::Matrix<my::nen_,1> epreaf(true);
  my::ExtractValuesFromGlobalVector(dis, lm, *my::rotsymmpbc_, &evelaf, &epreaf, "velaf");


  //----------------------------------------------------------------------------
  //      surface integral --- build Cuiui, Cuui, Cuiu and Cuu matrix and rhs
  //----------------------------------------------------------------------------

  DRT::Element::LocationArray alela( 1 );
  DRT::Element::LocationArray cutla( 1 );

  LINALG::Matrix<3,1> normal;
  LINALG::Matrix<3,1> x_side;

  bool fluidfluidcoupling = false;

  // embedded element coupling implementation between background element and each cutting embedded element (map<sid, emb_impl)
  std::map<int, Teuchos::RCP<DRT::ELEMENTS::XFLUID::EmbCoupling<distype> > > emb_impl;
  // side coupling implementation between background element and each cut side (map<sid, side_impl)
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

  // create location vectors for intersecting embedded elements and reshape coupling matrices
  PatchLocationVector(begids,alediscret,patchelementslmv,patchelementslmowner,Cuiui_coupling,boundary_emb_gid_map,"Nitsche");

  //-----------------------------------------------------------------------------------
  //         evaluate element length, stabilization factors and average weights
  //-----------------------------------------------------------------------------------

  // element length
  double h_k = 0.0;

  if(coupling_strategy != INPAR::XFEM::Embedded_Sided_Mortaring)
  {
    if(visc_stab_hk == INPAR::XFEM::ViscStab_hk_vol_equivalent)
    {
      LINALG::Matrix<1,1> dummy(true);
      h_k = FLD::UTILS::HK<distype>(dummy, my::xyze_);
    }
    else if(visc_stab_hk == INPAR::XFEM::ViscStab_hk_vol_div_by_surf)
    {
      h_k = meas_partial_volume / meas_surf;
    }
    else if(visc_stab_hk == INPAR::XFEM::ViscStab_hk_longest_ele_length)
    {
      dserror("longest element length for hk not supported yet");
    }
    else dserror("unknown type of characteristic element length");

    if( h_k <= 0.0 ) dserror("element length is <= 0.0");
  }


  double kappa1 = 1.0;      // Xfluid-sided mortaring

  // scaling factors for stabilization terms
  double visceff_max = my::visceff_;

  double stabfac_scaling      = 0.0;  // scaling factor for standard Nitsche stabilization
  double stabfac_conv_scaling = 0.0;  // scaling factor for convecitve Nitsche stabilization


  if(coupling_strategy == INPAR::XFEM::Xfluid_Sided_Mortaring)
  {

    stabfac_scaling      = nitsche_stab * visceff_max / h_k;  // scaling factor for standard Nitsche stabilization
    stabfac_conv_scaling = nitsche_stab_conv;            // scaling factor for convecitve Nitsche stabilization

    if( h_k <= 0.0 ) dserror("element length is <= 0.0");

    kappa1 = 1.0;
  }
  else if(coupling_strategy == INPAR::XFEM::Embedded_Sided_Mortaring)
  {
    kappa1 = 0.0;
  }
  else if(coupling_strategy == INPAR::XFEM::Two_Sided_Mortaring)
  {
    //if(meas_partial_volume < 1e-008) dserror(" measure of cut partial volume is smaller than 1d-008: %d Attention with increasing Nitsche-Parameter!!!", meas_partial_volume);

    if( h_k <= 0.0 ) dserror("element length is <= 0.0");

    stabfac_scaling      = nitsche_stab * visceff_max / h_k;  // scaling factor for standard Nitsche stabilization
    stabfac_conv_scaling = nitsche_stab_conv;            // scaling factor for convecitve Nitsche stabilization

    kappa1 = 0.5;
  }
  else dserror("coupling strategy not known");





  //------------------------------
  // define average weights
  if(meas_partial_volume < 0.0) dserror(" measure of cut partial volume is smaller than 0.0: %f Attention with increasing Nitsche-Parameter!!!", meas_partial_volume);

  if( kappa1 > 1.0 || kappa1 < 0.0) dserror("Nitsche weights for inverse estimate kappa1 lies not in [0,1]: %d", kappa1);

  double kappa2 = 1.0-kappa1;




  // evaluate shape function derivatives
  bool eval_deriv = true;

  //----------------------------------------------------------------------------
  //      surface integral --- loop sides
  //----------------------------------------------------------------------------
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

    // set displacement of embedded element
    emb->addembdisp(alediscret,"dispnp",alela[0].lm_);

    // set displacement of side
    si->addeidisp(cutdis,"idispnp",cutla[0].lm_);


    // set the embedded element length dependent on side in case of Emb Embedded_Sided_Mortaring
    if(coupling_strategy == INPAR::XFEM::Embedded_Sided_Mortaring)
    {
      emb->element_length(h_k);

      if(h_k < 1e-006) dserror("element length is smaller than 1e-006");

      stabfac_scaling      = nitsche_stab * visceff_max / h_k;  // scaling factor for standard Nitsche stabilization
      stabfac_conv_scaling = nitsche_stab_conv;            // scaling factor for convecitve Nitsche stabilization
    }




    //--------------------------------------------
    // loop boundary cells w.r.t current cut side
    //--------------------------------------------
    for ( std::vector<DRT::UTILS::GaussIntegration>::const_iterator i=cutintpoints.begin();
          i!=cutintpoints.end();
          ++i )
    {
      const DRT::UTILS::GaussIntegration & gi = *i;
      GEO::CUT::BoundaryCell * bc = bcs[i - cutintpoints.begin()]; // get the corresponding boundary cell

      //--------------------------------------------
      // loop gausspoints w.r.t current boundary cell
      //--------------------------------------------
      for ( DRT::UTILS::GaussIntegration::iterator iquad=gi.begin(); iquad!=gi.end(); ++iquad )
      {
        double drs = 0.0; // transformation factor between reference cell and linearized boundary cell

        const LINALG::Matrix<2,1> eta( iquad.Point() ); // xi-coordinates with respect to side

        LINALG::Matrix<3,1> rst(true); // local coordinates w.r.t background element

#ifdef BOUNDARYCELL_TRANSFORMATION_OLD

        si->Evaluate(eta,x_side,normal,drs);

        // find element local position of gauss point at interface
        GEO::CUT::Position<distype> pos( my::xyze_, x_side );
        pos.Compute();
        rst = pos.LocalCoordinates();

#else
        LINALG::Matrix<3,1> x_gp_lin(true); // gp in xyz-system on linearized interface

        // compute transformation factor, normal vector and global Gauss point coordiantes
        if(bc->Shape() != DRT::Element::dis_none) // Tessellation approach
        {
          ComputeSurfaceTransformation(drs, x_gp_lin, normal, bc, eta);
        }
        else // MomentFitting approach
        {
          drs = 1.0;
          normal = bc->GetNormalVector();
          const double* gpcord = iquad.Point();
          for (int idim=0;idim<3;idim++)
          {
            x_gp_lin(idim,0) = gpcord[idim];
          }
        }

        // find element local position of gauss point
        GEO::CUT::Position<distype> pos( my::xyze_, x_gp_lin );
        pos.Compute();
        rst = pos.LocalCoordinates();

        // project gaussian point from linearized interface to warped side (get/set local side coordinates in SideImpl)
        LINALG::Matrix<2,1> xi_side(true);
        si->ProjectOnSide(x_gp_lin, x_side, xi_side);
 #endif

        const double surf_fac = drs*iquad.Weight();

        const double timefacfac = surf_fac * my::fldpara_->TimeFac();



        // evaluate embedded element shape functions
        if(visc_stab_hk == INPAR::XFEM::ViscStab_hk_vol_equivalent)
        {
          emb->EvaluateEmb( x_side );
        }
        else dserror("choose vol_equivalent characteristic element length for embedded sided mortaring");

        //--------------------------------------------

        // evaluate shape functions (and derivatives)

        if(eval_deriv)
        {
          EvalFuncAndDeriv( rst );
        }
        else
        {
          DRT::UTILS::shape_function<distype>( rst, my::funct_ );
        }



        // get velocity at integration point
        // (values at n+alpha_F for generalized-alpha scheme, n+1 otherwise)
        my::velint_.Multiply(evelaf,my::funct_);

        // get velocity derivatives at integration point
        // (values at n+alpha_F for generalized-alpha scheme, n+1 otherwise)
        my::vderxy_.MultiplyNT(evelaf,my::derxy_);

        // get pressure at integration point
        // (value at n+alpha_F for generalized-alpha scheme, n+1 otherwise)
        double press = my::funct_.Dot(epreaf);

        //--------------------------------------------
        // compute stabilization factors

        double stabfac      = 0.0;
        double stabfac_conv = 0.0;

//        NIT_ComputeStabfac(stabfac, stabfac_conv, stabfac_scaling, stabfac_conv_scaling, my::velint_.Dot(normal));

        stabfac = stabfac_scaling;
        stabfac_conv = stabfac_conv_scaling;

        //--------------------------------------------

        if(fluidfluidcoupling)
        {
          bool bg_mortaring = false; // one-sided background fluid mortaring (kappa1=1, kappa2=0)

          if(coupling_strategy == INPAR::XFEM::Xfluid_Sided_Mortaring) bg_mortaring = true;

          //zero velocity jump for fluidfluidcoupling
          LINALG::Matrix<my::nsd_,1> ivelint_WDBC_JUMP(true);


          emb->NIT2_buildCouplingMatrices(
              elemat1_epetra,          // standard bg-bg-matrix
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
        else  dserror(" no two sided mortaring for non-fluidfluidcoupling");




      } // end loop gauss points of boundary cell
    } // end loop boundary cells of side
  } // end loop cut sides


  //----------------------------------------------------------------------------
  // build Cuiui coupling matrix (includes patch of Cuiui matrices for all sides)
  //----------------------------------------------------------------------------

  if ( fluidfluidcoupling )
    NIT_BuildPatchCuiui(Cuiui, Cuiui_coupling);

  return;
}


/*--------------------------------------------------------------------------------
 * build the patch coupling matrix Cuiui containing Cuiui for all cutting sides
 *--------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void FluidEleCalcXFEM<distype>::NIT_BuildPatchCuiui(
    Epetra_SerialDenseMatrix &                              Cuiui,            ///< ui-ui patch coupling matrix containing Cuiui for all cutting sides
    std::map<int, std::vector<Epetra_SerialDenseMatrix> >&  Cuiui_coupling    ///< Cuiui matrices for all cutting sides
)
{

  // build patch-Cuiui matrix
  int ipatchsizesbefore = 0;
  for (std::map<int, std::vector<Epetra_SerialDenseMatrix> >::const_iterator m=Cuiui_coupling.begin();
      m!=Cuiui_coupling.end(); ++m)
  {

    int bid = m->first;
    std::vector<Epetra_SerialDenseMatrix> & Cuiui_mats = Cuiui_coupling[bid];

    // Cuiui matrices in Cuiui_mats[0]

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

  return;
}

/*--------------------------------------------------------------------------------
 *    compute stabilization factor for Nitsche's method
 *--------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void FluidEleCalcXFEM<distype>::NIT_ComputeStabfac(
    bool                            fluidfluidcoupling,
    double &                        stabfac,              ///< stabilization factor 1
    double &                        stabfac_conv,         ///< stabilization factor 2
    const double                    NIT_stab_fac,         ///< stabilization factor for Nitsche term
    INPAR::XFEM::ConvStabScaling    conv_stab_scaling,
    const double                    nitsche_stab_conv,
    const double                    veln_normal
)
{
  // additional stabilization for convective stabilization
  double conv_stabfac = 0.0;

  if(conv_stab_scaling == INPAR::XFEM::ConvStabScaling_abs_normal_vel)
  {
    //      | u*n |
    conv_stabfac =  fabs(veln_normal);
  }
  else if(conv_stab_scaling == INPAR::XFEM::ConvStabScaling_inflow)
  {
    //      ( -u*n ) if (u*n)<0 (inflow)
    conv_stabfac = max(0.0,-veln_normal);
  }
  else if(conv_stab_scaling == INPAR::XFEM::ConvStabScaling_const)
  {
    //      const = conv_stab_fac
    conv_stabfac = 1.0;
  }
  else if(conv_stab_scaling == INPAR::XFEM::ConvStabScaling_none)
  {
    //      const = conv_stab_fac
    conv_stabfac = 0.0;
  }
  else dserror("unknown scaling for viscous stabilization term");

  if(!fluidfluidcoupling)
  {
    //=================================================================================
    // definition in Burman 2007
    // Interior penalty variational multiscale method for the incompressible Navier-Stokes equation:
    // Monitoring artificial dissipation
    /*
    //      viscous_Nitsche-part, convective inflow part
    //
    //                    mu                                               /       _      \
    //  max( gamma_Nit * ----  , gamma_conv * | u_h * n | )       *       |  u_h - u, v_h  |
    //                    h_k                                              \              /
    */

    // final stabilization factors
    stabfac      = max(NIT_stab_fac, nitsche_stab_conv*conv_stabfac);
    stabfac_conv = 0.0;
    //=================================================================================
  }
  else //fluidfluidcoupling
  {

    /*
    //      viscous_Nitsche-part
    //   /                                   \        /                           i   \
    //  |  gamma*mu/h_K *  [ v ] , [ Du ]     | =  - |   gamma*mu/h_K * [ v ], [ u ]   |
    //   \                                   /        \                               /
    */

    stabfac      = NIT_stab_fac;

    /*
    //    /                                           \        /                        i               \
    //   |  gamma_conv * | u * n | *  { v } , [ Du ]   | =  - |  gamma_conv*  | u * n | * { v }, [ u ]   |
    //    \                                           /        \                                        /
    //
    */

    if(nitsche_stab_conv != 0.0) dserror("no special inflow stabilization implemented for xff");
    stabfac_conv = nitsche_stab_conv*conv_stabfac;

  }

  return;
}

/*--------------------------------------------------------------------------------
 * pre-compute the measure of all side's surface cutting the element
 *--------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
double FluidEleCalcXFEM<distype>::ComputeMeasSurf(
    const std::map<int, std::vector<DRT::UTILS::GaussIntegration> > &   bintpoints,        ///< boundary cell integration points
    const std::map<int, std::vector<GEO::CUT::BoundaryCell*> > &        bcells             ///< boundary cells
    )
{
  double surf = 0.0;

  //--------------------------------------------
  // loop intersecting sides
  // map of side-element id and Gauss points
  for ( std::map<int, std::vector<DRT::UTILS::GaussIntegration> >::const_iterator i=bintpoints.begin();
        i!=bintpoints.end();
        ++i )
  {
    int sid = i->first;
    const std::vector<DRT::UTILS::GaussIntegration> & cutintpoints = i->second;

    // get side's boundary cells
    std::map<int, std::vector<GEO::CUT::BoundaryCell*> >::const_iterator j = bcells.find( sid );
    if ( j==bcells.end() )
      dserror( "missing boundary cell" );

    const std::vector<GEO::CUT::BoundaryCell*> & bcs = j->second;
    if ( bcs.size()!=cutintpoints.size() )
      dserror( "boundary cell integration rules mismatch" );

    //--------------------------------------------
    // loop boundary cells w.r.t current cut side
    //--------------------------------------------
    for ( std::vector<DRT::UTILS::GaussIntegration>::const_iterator i=cutintpoints.begin();
          i!=cutintpoints.end();
          ++i )
    {
      const DRT::UTILS::GaussIntegration & gi = *i;
      GEO::CUT::BoundaryCell * bc = bcs[i - cutintpoints.begin()]; // get the corresponding boundary cell

      //--------------------------------------------
      // loop gausspoints w.r.t current boundary cell
      //--------------------------------------------
      for ( DRT::UTILS::GaussIntegration::iterator iquad=gi.begin(); iquad!=gi.end(); ++iquad )
      {
        double drs = 0.0; // transformation factor between reference cell and linearized boundary cell

        const LINALG::Matrix<2,1> eta( iquad.Point() ); // xi-coordinates with respect to side

        LINALG::Matrix<3,1> normal(true);

#ifdef BOUNDARYCELL_TRANSFORMATION_OLD
        dserror("at the moment not available -> fix it");
//        si->Evaluate(eta,x_side,normal,drs);

#else
        LINALG::Matrix<3,1> x_gp_lin(true); // gp in xyz-system on linearized interface

        // compute transformation factor, normal vector and global Gauss point coordiantes
        if(bc->Shape() != DRT::Element::dis_none) // Tessellation approach
        {
          ComputeSurfaceTransformation(drs, x_gp_lin, normal, bc, eta);
        }
        else // MomentFitting approach
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

        const double surf_fac = drs*iquad.Weight();


        surf += surf_fac;

      } //loop gausspoints w.r.t current boundary cell
    } // loop boundary cells
  } // loop intersecting sides

  return surf;
}


/*--------------------------------------------------------------------------------
 * create location vector w.r.t patch of intersecting boundary elements and reshape coupling matrices
 *--------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void FluidEleCalcXFEM<distype>::PatchLocationVector(
    std::set<int> &                                         begids,                  ///< ids of intersecting boundary elements
    DRT::Discretization &                                   cutdis,                  ///< cut discretization
    std::vector<int> &                                      patchelementslmv,        ///< lm vector for patch of boundary elements
    std::vector<int> &                                      patchelementslmowner,    ///< lmowner vector for patch of boundary elements
    std::map<int, std::vector<Epetra_SerialDenseMatrix> > & Cuiui_coupling,          ///< coupling matrices
    string                                                  coupl_method             ///< coupling method
)
{
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

    if(coupl_method == "MixedStressHybrid")
    {
      Cuiui_matrices.resize(2);
      Cuiui_matrices[0].Reshape(my::nen_*6,patchlm.size()); //Gsui (coupling between background elements sigma and current side!)
      Cuiui_matrices[1].Reshape(patchlm.size(),my::nen_*6); //Guis
    }
    else if(coupl_method == "Nitsche")
    {
      Cuiui_matrices.resize(1);
      Cuiui_matrices[0].Reshape(patchlm.size(),patchlm.size()); //Cuiui
    }
    else dserror("not supported coupling method");

  }

  return;
}


/*--------------------------------------------------------------------------------
 * create location vector w.r.t patch of intersecting boundary elements and
 * reshape coupling matrices (for embedded coupling)
 *--------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void FluidEleCalcXFEM<distype>::PatchLocationVector(
    std::set<int> &                                         begids,                  ///< ids of intersecting boundary elements
    DRT::Discretization &                                   alediscret,                  ///< cut discretization
    std::vector<int> &                                      patchelementslmv,        ///< lm vector for patch of boundary elements
    std::vector<int> &                                      patchelementslmowner,    ///< lmowner vector for patch of boundary elements
    std::map<int, std::vector<Epetra_SerialDenseMatrix> > & Cuiui_coupling,          ///< coupling matrices
    map<int,int> &                                          boundary_emb_gid_map,    ///< map between boundary sid and corresponding embedded element id
    string                                                  coupl_method             ///< coupling method
)
{

  for (std::set<int>::const_iterator bgid=begids.begin(); bgid!=begids.end(); ++bgid)
  {
    DRT::Element * emb_ele = alediscret.gElement(boundary_emb_gid_map.find(*bgid)->second);

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

    if(coupl_method == "MixedStressHybrid")
    {
      dserror("embedded coupling for MixedStressHybrid not available!");
    }
    else if(coupl_method == "Nitsche")
    {
      Cuiui_matrices.resize(1);
      Cuiui_matrices[0].Reshape(patchlm.size(),patchlm.size()); //Cuiui
    }
    else dserror("not supported coupling method");

  }

  return;
}




/*--------------------------------------------------------------------------------
 * compute transformation factor for surface integration, normal, local and global gp coordinates
 *--------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void FluidEleCalcXFEM<distype>::ComputeSurfaceTransformation(
    double &                    drs,         ///< surface transformation factor
    LINALG::Matrix<3,1> &       x_gp_lin,    ///< global coordiantes of gaussian point
    LINALG::Matrix<3,1> &       normal,      ///< normal vector on boundary cell
    GEO::CUT::BoundaryCell *    bc,          ///< boundary cell
    const LINALG::Matrix<2,1> & eta          ///< local coordinates of gaussian point w.r.t boundarycell
)
{

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

  return;
}



/*----------------------------------------------------------------------*
 | evaluate shape functions and derivatives at given local coordinates  |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void FluidEleCalcXFEM<distype>::EvalFuncAndDeriv( LINALG::Matrix<3,1> &  rst )
{

  // evaluate shape functions
  DRT::UTILS::shape_function<distype>( rst, my::funct_ );

  // evaluate the derivatives of shape functions
  DRT::UTILS::shape_function_deriv1<distype>(rst,my::deriv_);
  my::xjm_.MultiplyNT(my::deriv_,my::xyze_);
  my::det_ = my::xji_.Invert(my::xjm_);

  // compute global first derivates
  my::derxy_.Multiply(my::xji_,my::deriv_);


  return;
}


/*----------------------------------------------------------------------*
 | build traction vector                                                |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void FluidEleCalcXFEM<distype>::buildTractionVector(
    LINALG::Matrix<my::nsd_,1> &  traction,   ///< traction vector
    double &                      press,      ///< pressure at gaussian point
    LINALG::Matrix<my::nsd_,1> &  normal      ///< normal vector
)
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
  //-------------------------------

  // t = ( -pI + 2mu eps(u) )*n^f
  traction.Clear();
  traction.Multiply(eps, normal);
  traction.Scale(2.0*my::visceff_);

  // add the pressure part
  traction.Update( -press, normal, 1.0);

  // we need normal vector on fluid
  traction.Scale(-1.0);

  return;
}

/*----------------------------------------------------------------------*
 | assemble side's interface force                                      |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void FluidEleCalcXFEM<distype>::AssembleInterfaceForce(
    Teuchos::RCP<Epetra_Vector>            iforcecol, ///< interface force column vector
    DRT::Discretization &                  cutdis,    ///< cut discretization
    std::vector<int> &                     lm,        ///< local dof map
    Epetra_SerialDenseVector &             iforce     ///< interface force vector
)
{

  const Epetra_Map* dofcolmap = cutdis.DofColMap();

  for (int idof = 0; idof < (int)(lm.size()); ++idof)
  {
      int gdof = lm[idof];

      // f^i = ( N^i, t ) = ( N^i, (-pI+2mu*eps(u))*n )
      (*iforcecol)[dofcolmap->LID(gdof)] += iforce[idof];
  }

  return;
}

///*--------------------------------------------------------------------------------
// * create location vector w.r.t patch of intersecting boundary elements and reshape coupling matrices
// *--------------------------------------------------------------------------------*/
//template <DRT::Element::DiscretizationType distype>
//Teuchos::RCP<DRT::ELEMENTS::XFLUID::SideInterface<distype> > FluidEleCalcXFEM<distype>::CreateSideImpl(
//    bool                                                    fluidfluidcoupling,
//    DRT::Element *                                          side,
//    Epetra_SerialDenseMatrix &                              side_xyze,
//    std::vector<Epetra_SerialDenseMatrix> &                 side_matrices,
//    std::map<int, std::vector<Epetra_SerialDenseMatrix> >   Cuiui_coupling
//    )
//{
//  Teuchos::RCP<DRT::ELEMENTS::XFLUID::SideInterface<distype> > si = Teuchos::null;
//
//  if(fluidfluidcoupling)
//  {
//    // coupling matrices between background element and one! side
//    Epetra_SerialDenseMatrix & C_uiu  = side_matrices[0];
//    Epetra_SerialDenseMatrix & C_uui  = side_matrices[1];
//    Epetra_SerialDenseMatrix & rhC_ui = side_matrices[2];
//
//    // coupling matrices between one side and itself via the element Kss
//    std::map<int,std::vector<Epetra_SerialDenseMatrix> >::iterator c2 = Cuiui_coupling.find( side->Id() );
//    std::vector<Epetra_SerialDenseMatrix> & Cuiui_matrices = c2->second;
//    Epetra_SerialDenseMatrix & eleGsui = Cuiui_matrices[0];
//    Epetra_SerialDenseMatrix & eleGuis = Cuiui_matrices[1];
//    Epetra_SerialDenseMatrix  eleGuisKssInv;
//
//    si = DRT::ELEMENTS::XFLUID::SideInterface<distype>::Impl(side,C_uiu,C_uui,rhC_ui,eleGsui,eleGuis,side_xyze);
//  }
//  else
//  {
//    si = DRT::ELEMENTS::XFLUID::SideInterface<distype>::Impl(side,side_xyze);
//  }
//
//  return si;
//}

/*--------------------------------------------------------------------------------
 *--------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void FluidEleCalcXFEM<distype>::CalculateContinuityXFEM(
    DRT::ELEMENTS::Fluid *              ele,            ///< fluid element
    DRT::Discretization &                dis,            ///< discretization
    const std::vector<int> &             lm,             ///< local map
    Epetra_SerialDenseVector&            elevec1_epetra, ///< element vector
    const DRT::UTILS::GaussIntegration & intpoints       ///< integration points
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
    DRT::ELEMENTS::Fluid *     ele,                ///< fluid element
    DRT::Discretization &       dis,                ///< discretization
    const std::vector<int> &    lm,                 ///< local map
    Epetra_SerialDenseVector&   elevec1_epetra      ///< element vector
  )
{

  CalculateContinuityXFEM(ele,
                          dis,
                          lm,
                          elevec1_epetra,
                          my::intpoints_);
}


  } // end namespace ELEMENTS
} // end namespace DRT




// Ursula is responsible for this comment!
template class DRT::ELEMENTS::FluidEleCalcXFEM<DRT::Element::hex8>;
template class DRT::ELEMENTS::FluidEleCalcXFEM<DRT::Element::hex20>;
//template class DRT::ELEMENTS::FluidEleCalcXFEM<DRT::Element::hex27>;
//template class DRT::ELEMENTS::FluidEleCalcXFEM<DRT::Element::tet4>;
//template class DRT::ELEMENTS::FluidEleCalcXFEM<DRT::Element::tet10>;
//template class DRT::ELEMENTS::FluidEleCalcXFEM<DRT::Element::wedge6>;
//template class DRT::ELEMENTS::FluidEleCalcXFEM<DRT::Element::pyramid5>;


