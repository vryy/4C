#ifdef STKADAPTIVE

#include "fluid_errorestimate.H"

#include "../stk_lib/stk_iterator.H"
#include "../stk_lib/stk_discret.H"
#include "../stk_refine/stk_mesh.H"

#include "../drt_fem_general/drt_utils_integration.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"

#include "../drt_lib/standardtypes_cpp.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_element.H"

#include "../drt_mat/matpar_parameter.H"
#include "../drt_mat/material.H"
#include "../drt_mat/newtonianfluid.H"

// #include "../linalg/linalg_utils.H"

namespace STK
{
  namespace FLD
  {
    // there should be a better place for this
    DRT::Element::DiscretizationType ElementShape( const CellTopologyData & topology );
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::Element::DiscretizationType STK::FLD::ElementShape( const CellTopologyData & topology )
{
  switch ( topology.key )
  {
  case shards::Hexahedron<20>::key :         return DRT::Element::hex20;
  case shards::Hexahedron<27>::key :         return DRT::Element::hex27;
  case shards::Hexahedron<8>::key :          return DRT::Element::hex8;
  case shards::Line<2>::key :                return DRT::Element::line2;
  case shards::Line<3>::key :                return DRT::Element::line3;
  case shards::Pyramid< 5>::key :            return DRT::Element::pyramid5;
  case shards::Quadrilateral<4>::key :       return DRT::Element::quad4;
  case shards::Quadrilateral<8>::key :       return DRT::Element::quad8;
  case shards::Quadrilateral<9>::key :       return DRT::Element::quad9;
  case shards::Tetrahedron< 4>::key :        return DRT::Element::tet4;
  case shards::Tetrahedron<10>::key :        return DRT::Element::tet10;
  case shards::Triangle<3>::key :            return DRT::Element::tri3;
  case shards::Triangle<6>::key :            return DRT::Element::tri6;
  case shards::Wedge< 6>::key :              return DRT::Element::wedge6;
  case shards::Wedge<15>::key :              return DRT::Element::wedge15;

  case shards::Node::key :
  case shards::Pyramid<13>::key :
  case shards::Pyramid<14>::key :
  case shards::ShellLine<2>::key :
  case shards::ShellLine<3>::key :
  case shards::ShellQuadrilateral<4>::key :
  case shards::ShellQuadrilateral<8>::key :
  case shards::ShellQuadrilateral<9>::key :
  case shards::ShellTriangle<3>::key :
  case shards::ShellTriangle<6>::key :
  case shards::Wedge<18>::key :

  default:
    dserror( "unknown topology" );
  }
  return DRT::Element::dis_none;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
STK::FLD::FluidJumpIntegrator::FluidJumpIntegrator( STK::Discretization & dis )
  : dis_( dis )
{
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STK::FLD::FluidJumpIntegrator::Integrate( stk::mesh::VectorField & coords,
                                               stk::mesh::VectorField & velnp,
                                               stk::mesh::ScalarField & pressure,
                                               stk::mesh::ScalarField & error,
                                               stk::mesh::ScalarField & volume,
                                               stk::mesh::VectorField & resvel,
                                               stk::mesh::ScalarField & respres )
{
  stk::mesh::BulkData & bulk = dis_.GetMesh().BulkData();
  stk::mesh::Part & active = dis_.GetMesh().ActivePart();
  //stk::mesh::Part & owned  = dis_.GetMesh().OwnedPart();

  // we need to
//   Teuchos::RCP<Epetra_Vector> colresidual = Teuchos::rcp( new Epetra_Vector( dis_.DofColMap() ) );
//   LINALG::Export( *residual, *colresidual );

  // iterate all elements
  for ( STK::ElementIterator i( bulk, active );
        not i.done();
        ++i )
  {
    FluidEvaluateHelper e( dis_, i.element(), i.topology(), coords, velnp, pressure, &resvel, &respres );

    // iterate all sides
    for ( UniqueSideIterator j( i ); not j.done(); ++j )
    {
      unsigned side_number = j.side_number();

      FluidSideEvaluateHelper side( e, side_number );

      if ( j.is_boundary_side() )
      {
      }
      else if ( j.is_inner_side() )
      {
        FluidEvaluateHelper oe( dis_, j.other_element(), j.other_topology(), coords, velnp, pressure );
        unsigned other_side_number = j.other_side_number();

        side.Integrate( error, e, side_number, oe, other_side_number, false );
      }
      else if ( j.is_small_hanging_side() )
      {
        FluidEvaluateHelper oe( dis_, j.other_element(), j.other_topology(), coords, velnp, pressure );
        unsigned other_side_number = j.other_side_number();

        side.Integrate( error, e, side_number, oe, other_side_number, true );
      }
      else if ( j.is_big_hanging_side() )
      {
        dserror( "should never happen with UniqueSideIterator" );
      }
      else
      {
        dserror( "side type undecided" );
      }
    }

    // Final element integration. Might override side contributions if necessary.
    e.Integrate( error, volume );
  }

  // Summation is done. Calcuate square root.

  for ( STK::ElementIterator i( bulk, active );
        not i.done();
        ++i )
  {
    stk::mesh::Entity & e = i.element();
    double & err = *reinterpret_cast<double*>( stk::mesh::field_data( error , e ) );
    err = sqrt( err );
  }

  // should communicate field values to aura
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
STK::FLD::FluidEvaluateHelper::FluidEvaluateHelper( STK::Discretization & dis,
                                                    stk::mesh::Entity & e,
                                                    const CellTopologyData & topology,
                                                    stk::mesh::VectorField & coords,
                                                    stk::mesh::VectorField & velnp,
                                                    stk::mesh::ScalarField & pressure,
                                                    stk::mesh::VectorField * resvel,
                                                    stk::mesh::ScalarField * respres )
  : e_( e ),
    topology_( topology ),
    nu_( 0.0 )
{
  stk::mesh::PairIterRelation nodes = e.relations( stk::mesh::Node );

  ExtractMyValues( xyz_,   coords,   nodes, genprob.ndim );
  ExtractMyValues( vel_,   velnp,    nodes, genprob.ndim );
  ExtractMyValues( press_, pressure, nodes );

  if ( resvel!=NULL and respres!=NULL )
  {
    ExtractMyValues( resvel_,   *resvel,  nodes, genprob.ndim );
    ExtractMyValues( respress_, *respres, nodes );
  }

  // create material on-demand
  MAT::PAR::Parameter * matpat = dis.MaterialParameter( e.bucket() );
  if ( matpat==NULL )
    dserror( "no material parameters" );
  Teuchos::RCP<MAT::Material> mat = matpat->CreateMaterial();

  if ( mat->MaterialType() != INPAR::MAT::m_fluid )
    dserror( "newtonian fluid material expected but got type %d", mat->MaterialType() );

  double mu = static_cast<MAT::NewtonianFluid*>( &*mat )->Viscosity();
  double dens = static_cast<MAT::NewtonianFluid*>( &*mat )->Density();
  nu_ = mu/dens;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STK::FLD::FluidEvaluateHelper::Integrate( stk::mesh::ScalarField & error, stk::mesh::ScalarField & volume )
{
  const int iel = NumNode();

  const int dim = genprob.ndim;

  if ( dim==2 )
  {
    Epetra_SerialDenseVector func( iel );
    Epetra_SerialDenseMatrix dxdxi( dim, dim );
    Epetra_SerialDenseMatrix xji( dim, dim );
    Epetra_SerialDenseMatrix deriv( dim, iel );
    Epetra_SerialDenseMatrix derxy( dim, iel );
    double normedRtau = 0;
    double normedRdiv = 0;
    Epetra_SerialDenseMatrix grad_u( dim, dim );

    const DRT::Element::DiscretizationType distype = ElementShape( topology() );

    double & err = *reinterpret_cast<double*>( stk::mesh::field_data( error , element() ) );
    double & vol = *reinterpret_cast<double*>( stk::mesh::field_data( volume, element() ) );

    DRT::UTILS::GaussRule2D gaussrule;

    switch ( topology().key )
    {
    case shards::Triangle<3>::key:
      gaussrule = DRT::UTILS::intrule_tri_37point;
      break;
    case shards::Quadrilateral<4>::key:
      gaussrule = DRT::UTILS::intrule_quad_4point;
      break;
    default:
      dserror( "unknown topology for gaussrule initialization" );
    }

    const DRT::UTILS::IntegrationPoints2D gaussquadrature( gaussrule );

    // integral of the L2 norm
    for ( int intpoint=0; intpoint<gaussquadrature.nquad; ++intpoint )
    {
      DRT::UTILS::shape_function_2D(func,
                                    gaussquadrature.qxg[intpoint][0],
                                    gaussquadrature.qxg[intpoint][1],
                                    distype);
      DRT::UTILS::shape_function_2D_deriv1(deriv,
                                           gaussquadrature.qxg[intpoint][0],
                                           gaussquadrature.qxg[intpoint][1],
                                           distype);

      dxdxi.Multiply( 'N', 'T', 1.0, xyz(), deriv, 0.0 );
      double det = dxdxi( 0, 0 )*dxdxi( 1, 1 ) - dxdxi( 0, 1 )*dxdxi( 1, 0 );

      xji(0,0) =  dxdxi(1,1)/det;
      xji(1,0) = -dxdxi(1,0)/det;
      xji(0,1) = -dxdxi(0,1)/det;
      xji(1,1) =  dxdxi(0,0)/det;

      derxy.Multiply( 'N', 'N', 1.0, xji, deriv, 0.0 );

      grad_u.Multiply( 'N', 'T', 1.0, vel(), derxy, 0.0 );

      double Rtau_0 = 0;
      double Rtau_1 = 0;
      double Rtau_2 = 0;
      for ( int i=0; i<iel; ++i )
      {
        double f = func( i );
        Rtau_0 += f * resvel_( 0, i );
        Rtau_1 += f * resvel_( 1, i );

        // add pressure residual here?
        Rtau_2 += f * respress_( i );
      }

      double fac = det * gaussquadrature.qwgt[intpoint];
      normedRtau += ( Rtau_0*Rtau_0 + Rtau_1*Rtau_1 + Rtau_2*Rtau_2 ) * fac;

      double v1 = grad_u(0,0)+grad_u(1,1);
      normedRdiv += v1 * v1 * fac;

      //double v2 = -0.5 * ( grad_u(0,1)-grad_u(1,0) );
      //w += v2 * v2 * fac;
    }

    // R1 = alpha*||Rtau||²
    // R2 = || Rdiv ||² = || div u ||²
    // R4 = R_vorticity  = || curl u ||² or normed_w = || curl u ||

    double alpha = Area();

    double R1 = alpha * normedRtau;

    double R2 = normedRdiv;

    err += R1 + R2;
    vol  = alpha;
  }
  else
    dserror( "unsupported dimension" );
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STK::FLD::FluidEvaluateHelper::LocalCoordinates( Epetra_SerialDenseMatrix & rst, unsigned side_number )
{
  unsigned dim = topology().dimension;
  const CellTopologyData_Subcell & side = topology().subcell[dim-1][side_number] ;
  const CellTopologyData * side_topology = side.topology;

  rst.Reshape( side_topology->node_count, dim );

  switch ( topology().key )
  {
  case shards::Triangle<3>::key :
  case shards::Triangle<6>::key :
  {
    for ( unsigned j=0; j<side_topology->node_count; ++j )
    {
      int elenodenumber = side.node[j];
      for ( unsigned i=0; i<dim; ++i )
        rst( j, i ) = DRT::UTILS::eleNodeNumbering_tri6_nodes_reference[elenodenumber][i];
    }
    break;
  }

  case shards::Quadrilateral<4>::key :
  case shards::Quadrilateral<8>::key :
  case shards::Quadrilateral<9>::key :
  {
    for ( unsigned j=0; j<side_topology->node_count; ++j )
    {
      int elenodenumber = side.node[j];
      for ( unsigned i=0; i<dim; ++i )
        rst( j, i ) = DRT::UTILS::eleNodeNumbering_quad9_nodes_reference[elenodenumber][i];
    }
    break;
  }

  case shards::Hexahedron<8>::key :
  case shards::Hexahedron<20>::key :
  case shards::Hexahedron<27>::key :
  {
    for ( unsigned j=0; j<side_topology->node_count; ++j )
    {
      int elenodenumber = side.node[j];
      for ( unsigned i=0; i<dim; ++i )
        rst( j, i ) = DRT::UTILS::eleNodeNumbering_hex27_nodes_reference[elenodenumber][i];
    }
    break;
  }

  case shards::Tetrahedron< 4>::key :
  case shards::Tetrahedron<10>::key :
  {
    for ( unsigned j=0; j<side_topology->node_count; ++j )
    {
      int elenodenumber = side.node[j];
      for ( unsigned i=0; i<dim; ++i )
        rst( j, i ) = DRT::UTILS::eleNodeNumbering_tet10_nodes_reference[elenodenumber][i];
    }
    break;
  }

  case shards::Pyramid< 5>::key :

  case shards::Pyramid<13>::key :

  case shards::Pyramid<14>::key :

  case shards::Wedge< 6>::key :

  case shards::Wedge<15>::key :

  case shards::Wedge<18>::key :

  case shards::Line<2>::key :

  case shards::Node::key :

  case shards::Line<3>::key :

  case shards::ShellLine<2>::key :

  case shards::ShellLine<3>::key :

  case shards::ShellTriangle<3>::key :

  case shards::ShellTriangle<6>::key :

  case shards::ShellQuadrilateral<4>::key :

  case shards::ShellQuadrilateral<8>::key :

  case shards::ShellQuadrilateral<9>::key :

  default:
    dserror( "unknown element topology" );
  }
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STK::FLD::FluidEvaluateHelper::ExtractMyValues( Epetra_SerialDenseMatrix & emat,
                                                     stk::mesh::FieldBase & f,
                                                     const stk::mesh::PairIterRelation & nodes,
                                                     int dim )
{
  unsigned numnodes = nodes.size();
  int err = emat.Reshape( dim, numnodes );
  if ( err!=0 )
    dserror( "Reshape failed: err=%d", err );

  for ( unsigned i=0; i < numnodes; ++i )
  {
    stk::mesh::Entity & n = *nodes[i].entity();
    double * c = reinterpret_cast<double*>( stk::mesh::field_data( f , n ) );
    for ( int d=0; d<dim; ++d )
    {
      emat( d, i ) = c[d];
    }
  }
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STK::FLD::FluidEvaluateHelper::ExtractMyValues( Epetra_SerialDenseVector & evec,
                                                     stk::mesh::FieldBase & f,
                                                     const stk::mesh::PairIterRelation & nodes )
{
  unsigned numnodes = nodes.size();
  int err = evec.Resize( numnodes );
  if ( err!=0 )
    dserror( "Reshape failed: err=%d", err );

  for ( unsigned i=0; i < numnodes; ++i )
  {
    stk::mesh::Entity & n = *nodes[i].entity();
    double * c = reinterpret_cast<double*>( stk::mesh::field_data( f , n ) );
    evec( i ) = *c;
  }
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
double STK::FLD::FluidEvaluateHelper::Area()
{
  // quadratic elements are treated as linear shaped elements

  switch ( topology().key )
  {
  case shards::Triangle<3>::key:
  case shards::Triangle<6>::key:
  {
    double x1 = xyz_(0,0);
    double y1 = xyz_(1,0);

    double x2 = xyz_(0,1);
    double y2 = xyz_(1,1);

    double x3 = xyz_(0,2);
    double y3 = xyz_(1,2);

    return x1*y2/2 + x2*y3/2 + x3*y1/2 - x1*y3/2 - x2*y1/2 - x3*y2/2;
  }
  case shards::Quadrilateral<4>::key:
  case shards::Quadrilateral<8>::key:
  case shards::Quadrilateral<9>::key:
  {
    double x1 = xyz_(0,0);
    double y1 = xyz_(1,0);

    double x2 = xyz_(0,1);
    double y2 = xyz_(1,1);

    double x3 = xyz_(0,2);
    double y3 = xyz_(1,2);

    double x4 = xyz_(0,3);
    double y4 = xyz_(1,3);

    return x1*y2/2 + x2*y3/2 + x3*y4/2 + x4*y1/2 - x1*y4/2 - x2*y1/2 - x3*y2/2 - x4*y3/2;
  }
  default:
    dserror( "unsupported topology" );
  }
  return 0.;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
STK::FLD::FluidSideEvaluateHelper::FluidSideEvaluateHelper( STK::FLD::FluidEvaluateHelper & e,
                                                            unsigned side_number )
{
  unsigned dim = e.topology().dimension - 1;
  const CellTopologyData_Subcell & side = e.topology().subcell[dim][side_number] ;
  const CellTopologyData * side_topology = side.topology;

  ExtractMyValues( e.xyz(),   xyz_,   side, side_topology );
  ExtractMyValues( e.vel(),   vel_,   side, side_topology );
  ExtractMyValues( e.press(), press_, side, side_topology );
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STK::FLD::FluidSideEvaluateHelper::Integrate( stk::mesh::ScalarField & error,
                                                   STK::FLD::FluidEvaluateHelper & le,
                                                   unsigned l_side_number,
                                                   STK::FLD::FluidEvaluateHelper & re,
                                                   unsigned r_side_number,
                                                   bool hanging )
{
  const int inode = NumNode();
  const int l_inode = le.NumNode();
  const int r_inode = re.NumNode();

  const int dim = genprob.ndim;

  const CellTopologyData_Subcell & l_side = le.topology().subcell[dim-1][l_side_number] ;
  const CellTopologyData_Subcell & r_side = re.topology().subcell[dim-1][r_side_number] ;

  const CellTopologyData * side_topology = l_side.topology;

  DRT::Element::DiscretizationType sideshape = ElementShape( *side_topology );
  DRT::Element::DiscretizationType l_sideshape = ElementShape( le.topology() );
  DRT::Element::DiscretizationType r_sideshape = ElementShape( re.topology() );

  Epetra_SerialDenseMatrix l_rst;
  Epetra_SerialDenseMatrix r_rst;

  le.LocalCoordinates( l_rst, l_side_number );
  re.LocalCoordinates( r_rst, r_side_number );

  if ( hanging )
  {
    // We are looking at a hanging node side. The right side element is the
    // large one. We need to set r_rst to the local coordiantes of the small
    // side's nodes.

    // Since the orientation of the side element is arbitrary, we cannot just
    // assign local coordiantes here. We need to calculate coordiantes from
    // known ones.

    stk::mesh::PairIterRelation l_nodes = le.element().relations( stk::mesh::Node );
    stk::mesh::PairIterRelation r_nodes = re.element().relations( stk::mesh::Node );

    switch ( side_topology->key )
    {
    case shards::Line<2>::key :
    {
      double dr = r_rst( 1, 0 )-r_rst( 0, 0 );
      double ds = r_rst( 1, 1 )-r_rst( 0, 1 );

      stk::mesh::Entity & l_n0 = * l_nodes[l_side.node[0]].entity();
      stk::mesh::Entity & l_n1 = * l_nodes[l_side.node[1]].entity();
      stk::mesh::Entity & r_n0 = * r_nodes[r_side.node[0]].entity();
      stk::mesh::Entity & r_n1 = * r_nodes[r_side.node[1]].entity();

      if ( l_n0.key() == r_n0.key() or
           l_n1.key() == r_n0.key() )
      {
        r_rst( 1, 0 ) = r_rst( 0, 0 ) + .5*dr;
        r_rst( 1, 1 ) = r_rst( 0, 1 ) + .5*ds;
      }
      else if ( l_n0.key() == r_n1.key() or
                l_n1.key() == r_n1.key() )
      {
        r_rst( 0, 0 ) = r_rst( 1, 0 ) - .5*dr;
        r_rst( 0, 1 ) = r_rst( 1, 1 ) - .5*ds;
      }
      else
        dserror( "line2 nodes do not match" );
      break;
    }
    default:
      dserror( "topology not supported" );
    }
  }

  if ( dim==2 )
  {
    DRT::UTILS::GaussRule1D gaussrule;

    switch ( side_topology->key )
    {
    case shards::Line<2>::key:
      gaussrule = DRT::UTILS::intrule_line_1point;
      break;
    case shards::Line<3>::key:
      gaussrule = DRT::UTILS::intrule_line_2point;
      dserror( "actually, the normal vector cannot be calculated here" );
      break;
    default:
      dserror( "unknown topology for gaussrule initialization" );
    }

    // for curved surfaces we would have to think anew at each gauss point
    Epetra_SerialDenseVector n_vector( dim );
    n_vector( 0 ) = (  1 )*( xyz_( 1, 1 ) - xyz_( 1, 0 ) );
    n_vector( 1 ) = ( -1 )*( xyz_( 0, 1 ) - xyz_( 0, 0 ) );
    double length = n_vector.Norm2();
    n_vector.Scale( 1. / length );

    const DRT::UTILS::IntegrationPoints1D gaussquadrature( gaussrule );

    Epetra_SerialDenseVector sidefunc( inode );
    Epetra_SerialDenseMatrix sidederiv( dim-1, inode );

    Epetra_SerialDenseVector l_func( l_inode );
    Epetra_SerialDenseVector r_func( r_inode );
    Epetra_SerialDenseMatrix l_deriv( dim, l_inode );
    Epetra_SerialDenseMatrix r_deriv( dim, r_inode );

    Epetra_SerialDenseVector l_xi( dim );
    Epetra_SerialDenseVector r_xi( dim );

    Epetra_SerialDenseMatrix l_dxdxi( dim, dim );
    Epetra_SerialDenseMatrix r_dxdxi( dim, dim );

    Epetra_SerialDenseMatrix l_xji( dim, dim );
    Epetra_SerialDenseMatrix r_xji( dim, dim );

    Epetra_SerialDenseVector l_dxideta( dim );
    Epetra_SerialDenseVector r_dxideta( dim );
    Epetra_SerialDenseVector l_dxdeta( dim );
    Epetra_SerialDenseVector r_dxdeta( dim );
    Epetra_SerialDenseMatrix l_derxy( dim, l_inode );
    Epetra_SerialDenseMatrix r_derxy( dim, r_inode );

    double l_det;
    double r_det;

    Epetra_SerialDenseMatrix l_grad_u( dim, dim );
    Epetra_SerialDenseMatrix r_grad_u( dim, dim );

    double l_p;
    double r_p;

    Epetra_SerialDenseVector left( dim );
    Epetra_SerialDenseVector right( dim );

    Epetra_SerialDenseVector jump( dim );

    double & l_err = *reinterpret_cast<double*>( stk::mesh::field_data( error , le.element() ) );
    double & r_err = *reinterpret_cast<double*>( stk::mesh::field_data( error , re.element() ) );

    // integrate over the side
    for ( int intpoint=0; intpoint<gaussquadrature.nquad; ++intpoint )
    {
      // coordiantes of the current (line-)integration points, eta
      const double eta1 = gaussquadrature.qxg[intpoint][0];

      // side-shape functions and their derivatives
      DRT::UTILS::shape_function_1D( sidefunc, eta1, sideshape );
      DRT::UTILS::shape_function_1D_deriv1( sidederiv, eta1, sideshape );

      // mapping the gauss points coordinates from 1D eta to 2D xi |
      // xi(dim) = rst(inode,dim)^T . sidefunc(inode)
      l_rst.Multiply( true, sidefunc, l_xi );
      r_rst.Multiply( true, sidefunc, r_xi );

      DRT::UTILS::shape_function_2D( l_func, l_xi(0), l_xi(1), l_sideshape );
      DRT::UTILS::shape_function_2D( r_func, r_xi(0), r_xi(1), r_sideshape );

      DRT::UTILS::shape_function_2D_deriv1( l_deriv, l_xi(0), l_xi(1), l_sideshape );
      DRT::UTILS::shape_function_2D_deriv1( r_deriv, r_xi(0), r_xi(1), r_sideshape );

      // mapping from global to 2D | dxdxi(dim,dim) = xyz(dim,leftinode) . deriv(dim,leftinode)^T
      l_dxdxi.Multiply( 'N', 'T', 1.0, le.xyz(), l_deriv, 0.0 );
      r_dxdxi.Multiply( 'N', 'T', 1.0, re.xyz(), r_deriv, 0.0 );

      // from 2D to 1D | dxideta(dim) = rst(inode,dim)^T . sidederiv(1,inode)^T
      l_dxideta.Multiply( 'T', 'T', 1.0, l_rst, sidederiv, 0.0 );
      r_dxideta.Multiply( 'T', 'T', 1.0, r_rst, sidederiv, 0.0 );

      // finally from global to 1D | dxdeta(dim) = dxdxi(dim,dim) . dxideta(dim)
      l_dxdeta.Multiply( 'N', 'N', 1.0, l_dxdxi, l_dxideta, 0.0 );
      r_dxdeta.Multiply( 'N', 'N', 1.0, r_dxdxi, r_dxideta, 0.0 );

      // jacobian in 2D
      //fluid2->f2_jaco( le.xyz(), l_deriv, l_xjacob, &l_det, l_inode );
      //fluid2->f2_jaco( re.xyz(), r_deriv, r_xjacob, &r_det, r_inode );

      //l_xjacob.Multiply( 'N', 'T', 1.0, le.xyz(), l_deriv, 0.0 );
      //r_xjacob.Multiply( 'N', 'T', 1.0, re.xyz(), r_deriv, 0.0 );

      l_det = l_dxdxi( 0, 0 )*l_dxdxi( 1, 1 ) - l_dxdxi( 0, 1 )*l_dxdxi( 1, 0 );
      r_det = r_dxdxi( 0, 0 )*r_dxdxi( 1, 1 ) - r_dxdxi( 0, 1 )*r_dxdxi( 1, 0 );

      // dN / dx , derivatives of shape funct.s by global coord.s
      //fluid2->f2_gder( l_derxy, l_deriv, l_xjacob, l_detjacob, l_inode ); //global derivatives
      //fluid2->f2_gder( r_derxy, r_deriv, r_xjacob, r_detjacob, r_inode );

      l_xji(0,0) =  l_dxdxi(1,1)/l_det;
      l_xji(1,0) = -l_dxdxi(1,0)/l_det;
      l_xji(0,1) = -l_dxdxi(0,1)/l_det;
      l_xji(1,1) =  l_dxdxi(0,0)/l_det;

      r_xji(0,0) =  r_dxdxi(1,1)/r_det;
      r_xji(1,0) = -r_dxdxi(1,0)/r_det;
      r_xji(0,1) = -r_dxdxi(0,1)/r_det;
      r_xji(1,1) =  r_dxdxi(0,0)/r_det;

      l_derxy.Multiply( 'N', 'N', 1.0, l_xji, l_deriv, 0.0 );
      r_derxy.Multiply( 'N', 'N', 1.0, r_xji, r_deriv, 0.0 );

      // calculating the grad u

      l_grad_u.Multiply( 'N', 'T', 1.0, le.vel(), l_derxy, 0.0 );
      r_grad_u.Multiply( 'N', 'T', 1.0, re.vel(), r_derxy, 0.0 );

      // due to continuous pressure interpolation, l_p==r_p, but for the sake
      // of generality
      l_p = le.press().Dot( l_func );
      r_p = re.press().Dot( r_func );

      // let's do the || grad u . n - p . n || for left and right element (on the side !)
      left = n_vector;
      left.Multiply( 'N', 'N', le.nu(), l_grad_u, n_vector, -l_p );
      left.Scale( l_dxdeta.Norm2() );

      right = n_vector;
      right.Multiply( 'N', 'N', re.nu(), r_grad_u, n_vector, -r_p );
      right.Scale( (-1)*r_dxdeta.Norm2() );

      jump = left;
      jump += right;

      double norm = jump.Norm2();
      double normedReps = norm*norm*gaussquadrature.qwgt[intpoint];

      // we sum the squares of our error guess here

      l_err += 0.5*length*normedReps;
      r_err += 0.5*length*normedReps;
    }
  }
  else if ( dim==3 )
  {
    dserror( "not implemented" );
  }
  else
  {
    dserror( "unsupported problem dimension" );
  }
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STK::FLD::FluidSideEvaluateHelper::ExtractMyValues( const Epetra_SerialDenseMatrix & emat,
                                                         Epetra_SerialDenseMatrix & sidemat,
                                                         const CellTopologyData_Subcell & side,
                                                         const CellTopologyData * side_topology )
{
  int dim = emat.M();
  sidemat.Reshape( dim, side_topology->node_count );
  for ( unsigned j=0; j<side_topology->node_count; ++j )
  {
    unsigned col = side.node[j];
    for ( int i=0; i<dim; ++i )
    {
      sidemat( i, j ) = emat( i, col );
    }
  }
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STK::FLD::FluidSideEvaluateHelper::ExtractMyValues( const Epetra_SerialDenseVector & evec,
                                                         Epetra_SerialDenseVector & sidevec,
                                                         const CellTopologyData_Subcell & side,
                                                         const CellTopologyData * side_topology )
{
  sidevec.Resize( side_topology->node_count );
  for ( unsigned j=0; j<side_topology->node_count; ++j )
  {
    unsigned col = side.node[j];
    sidevec( j ) = evec( col );
  }
}


#endif
