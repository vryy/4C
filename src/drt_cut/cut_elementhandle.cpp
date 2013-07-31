
#include <Teuchos_TimeMonitor.hpp>

#include "cut_integrationcell.H"
#include "cut_meshintersection.H"
#include "cut_position.H"
#include "cut_volumecell.H"

#include "../drt_inpar/inpar_xfem.H"

#include<fstream>


/*-------------------------------------------------------------------------------------------------------------*
          Project the integration rule available in the local coordinates of the integationcells to the
          local coordinates of volumecells
 *-------------------------------------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
Teuchos::RCP<DRT::UTILS::GaussPoints> GEO::CUT::ElementHandle::CreateProjected( GEO::CUT::IntegrationCell * ic )
{
  const unsigned nen = DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;

  LINALG::Matrix<3, nen> xie;

  const std::vector<GEO::CUT::Point*> & cpoints = ic->Points();
  if ( cpoints.size() != nen )
    throw std::runtime_error( "non-matching number of points" );

  for ( unsigned i=0; i<nen; ++i )
  {
    GEO::CUT::Point * p = cpoints[i];
    const LINALG::Matrix<3,1> & xi = LocalCoordinates( p );
    std::copy( xi.A(), xi.A()+3, &xie( 0, i ) );
  }

  Teuchos::RCP<DRT::UTILS::GaussPoints> gp =
    DRT::UTILS::GaussIntegration::CreateProjected<distype>( xie, ic->CubatureDegree( Shape() ) );
  return gp;
}

/*----------------------------------------------------------------------------------------------------------*
                  Collect the Gauss points of all the volumecells belong to this element
                  in such a way that gaussian rule for every vol.cell can be separated
*-----------------------------------------------------------------------------------------------------------*/
void GEO::CUT::ElementHandle::VolumeCellGaussPoints( plain_volumecell_set & cells,
                                                     std::vector<DRT::UTILS::GaussIntegration> & intpoints,
                                                     std::string gausstype )
{
  /*cells.clear(); //check whether any problems with gmsh output
  GetVolumeCells( cells );*/

  intpoints.clear();
  intpoints.reserve( cells.size() );


  if(gausstype == "Tessellation")
  {
        for ( plain_volumecell_set::iterator i=cells.begin(); i!=cells.end(); ++i )
        {
            GEO::CUT::VolumeCell * vc = *i;

            Teuchos::RCP<DRT::UTILS::GaussPointsComposite> gpc =
                Teuchos::rcp( new DRT::UTILS::GaussPointsComposite( 0 ) );

            const plain_integrationcell_set & cells = vc->IntegrationCells();
            for ( plain_integrationcell_set::const_iterator i=cells.begin(); i!=cells.end(); ++i )
            {
              GEO::CUT::IntegrationCell * ic = *i;
              switch ( ic->Shape() )
              {
              case DRT::Element::hex8:
              {
                Teuchos::RCP<DRT::UTILS::GaussPoints> gp = CreateProjected<DRT::Element::hex8>( ic );
                gpc->Append( gp );
                break;
              }
              case DRT::Element::tet4:
              {
                Teuchos::RCP<DRT::UTILS::GaussPoints> gp = CreateProjected<DRT::Element::tet4>( ic );
                gpc->Append( gp );
                break;
              }
              case DRT::Element::wedge6:
              {
                Teuchos::RCP<DRT::UTILS::GaussPoints> gp = CreateProjected<DRT::Element::wedge6>( ic );
                gpc->Append( gp );
                break;
              }
              case DRT::Element::pyramid5:
              {
                Teuchos::RCP<DRT::UTILS::GaussPoints> gp = CreateProjected<DRT::Element::pyramid5>( ic );
                gpc->Append( gp );
                break;
              }
              default:
                throw std::runtime_error( "unsupported integration cell type" );
              }
          }

        intpoints.push_back( DRT::UTILS::GaussIntegration( gpc ) );
        }
  }

  else if( gausstype == "MomentFitting" || gausstype=="DirectDivergence" )
  {
       for(plain_volumecell_set::iterator i=cells.begin(); i!=cells.end(); ++i)
       {
               GEO::CUT::VolumeCell *vc = *i;
               Teuchos::RCP<DRT::UTILS::GaussPoints> gp = vc->GetGaussRule();
               Teuchos::RCP<DRT::UTILS::GaussPointsComposite> gpc =
                        Teuchos::rcp( new DRT::UTILS::GaussPointsComposite( 0 ) );
               gpc->Append(gp);
               intpoints.push_back( DRT::UTILS::GaussIntegration( gpc ) );
       }
  }


  /*if(IsCut())
  {
	  static int eeno=0;
	  eeno++;
	  if(eeno==1 || eeno==2 || eeno==3)
	  {
  for(std::vector<DRT::UTILS::GaussIntegration>::iterator i=intpoints.begin();i!=intpoints.end();i++)
  {
	  DRT::UTILS::GaussIntegration ga = *i;
	  static int sideno = 0;
          sideno++;
	  std::string filename="wrong";
   	  std::ofstream file;

          std::stringstream out;
          out <<"gauss"<<sideno<<".dat";
          filename = out.str();
          file.open(filename.c_str());
	  for ( DRT::UTILS::GaussIntegration::const_iterator iquad=ga.begin(); iquad!=ga.end(); ++iquad )
	  {
		  const double* gpp = iquad.Point();
		  file<<gpp[0]<<"\t"<<gpp[1]<<"\t"<<gpp[2]<<"\t";
		  file<<iquad.Weight()<<std::endl;
	  }
	  file.close();
  }
  }
  }*/

#if 0

  // assume normal integration of element
  DRT::UTILS::GaussIntegration element_intpoints( Shape() );

  // see if we benefit from a "negative" integartion
  int numgps = element_intpoints.NumPoints();
  for ( std::vector<DRT::UTILS::GaussIntegration>::iterator i=intpoints.begin(); i!=intpoints.end(); ++i )
  {
    const DRT::UTILS::GaussIntegration & gi = *i;
    numgps += gi.NumPoints();
  }

  for ( std::vector<DRT::UTILS::GaussIntegration>::iterator i=intpoints.begin(); i!=intpoints.end(); ++i )
  {
    DRT::UTILS::GaussIntegration & gi = *i;
    if ( gi.NumPoints() > numgps / 2 )
    {
      // Here we have the expensive bastard.
      // Exchange gauss points. Use all other points instead.

      Teuchos::RCP<DRT::UTILS::GaussPointsComposite> gpc =
        Teuchos::rcp( new DRT::UTILS::GaussPointsComposite( intpoints.size() ) );

      gpc->Append( element_intpoints.Points() );

      for ( std::vector<DRT::UTILS::GaussIntegration>::iterator j=intpoints.begin(); j!=intpoints.end(); ++j )
      {
        const DRT::UTILS::GaussIntegration & gi = *j;
        if ( i!=j )
        {
          gpc->Append( gi.Points() );
        }
      }

      gi.SetPoints( gpc );

      break;
    }
  }
#endif

//   if ( not include_inner )
//   {
//     int count = 0;
//     for ( plain_volumecell_set::iterator i=cells.begin(); i!=cells.end(); ++i )
//     {
//       GEO::CUT::VolumeCell * vc = *i;
//       if ( vc->Position()==Point::inside )
//       {
//         intpoints[count].Clear();
//       }
//       count += 1;
//     }
//   }
}


void GEO::CUT::ElementHandle::GetVolumeCellsDofSets ( std::vector<plain_volumecell_set > & cellsets ,
                                                      std::vector< std::vector< int > >  & nds_sets )
{
	bool include_inner = false;

  const std::vector<plain_volumecell_set> & ele_vc_sets_inside = GetVcSetsInside();
  const std::vector<plain_volumecell_set> & ele_vc_sets_outside = GetVcSetsOutside();

  std::vector<std::vector<int> > & nodaldofset_vc_sets_inside = GetNodalDofSet_VcSets_Inside();
  std::vector<std::vector<int> > & nodaldofset_vc_sets_outside = GetNodalDofSet_VcSets_Outside();

	if(include_inner)
	{
		std::copy( ele_vc_sets_inside.begin(), ele_vc_sets_inside.end(), std::inserter(cellsets,cellsets.end()) );
		std::copy( nodaldofset_vc_sets_inside.begin(), nodaldofset_vc_sets_inside.end(), std::inserter(nds_sets,nds_sets.end()) );
	}

	std::copy( ele_vc_sets_outside.begin(), ele_vc_sets_outside.end(), std::inserter(cellsets,cellsets.end()) );
	std::copy( nodaldofset_vc_sets_outside.begin(), nodaldofset_vc_sets_outside.end(), std::inserter(nds_sets,nds_sets.end()) );


}

/*---------------------------------------------------------------------------------------------------------------------------*
                        Collect the Gauss points of all the volumecells belong to this element
                        The integration rule over all the volumecells are connected.
                        Cannot identify the integration rules separately for each cell
 *---------------------------------------------------------------------------------------------------------------------------*/
Teuchos::RCP<DRT::UTILS::GaussPointsComposite> GEO::CUT::ElementHandle::GaussPointsConnected( plain_volumecell_set & cells,
                                                                                              std::string gausstype )
{

    Teuchos::RCP<DRT::UTILS::GaussPointsComposite> gpc =
      Teuchos::rcp( new DRT::UTILS::GaussPointsComposite( 0 ) );

  for ( plain_volumecell_set::iterator i=cells.begin(); i!=cells.end(); ++i )
  {
    GEO::CUT::VolumeCell * vc = *i;

    const plain_integrationcell_set & cells = vc->IntegrationCells();


    if(gausstype == "Tessellation")
    {
        for ( plain_integrationcell_set::const_iterator i=cells.begin(); i!=cells.end(); ++i )
        {
          GEO::CUT::IntegrationCell * ic = *i;
          switch ( ic->Shape() )
          {
          case DRT::Element::hex8:
          {
            Teuchos::RCP<DRT::UTILS::GaussPoints> gp = CreateProjected<DRT::Element::hex8>( ic );
            gpc->Append( gp );
            break;
          }
          case DRT::Element::tet4:
          {
            Teuchos::RCP<DRT::UTILS::GaussPoints> gp = CreateProjected<DRT::Element::tet4>( ic );
            gpc->Append( gp );
            break;
          }
          case DRT::Element::wedge6:
          {
            Teuchos::RCP<DRT::UTILS::GaussPoints> gp = CreateProjected<DRT::Element::wedge6>( ic );
            gpc->Append( gp );
            break;
          }
          case DRT::Element::pyramid5:
          {
            Teuchos::RCP<DRT::UTILS::GaussPoints> gp = CreateProjected<DRT::Element::pyramid5>( ic );
            gpc->Append( gp );
            break;
          }
          default:
            throw std::runtime_error( "unsupported integration cell type" );
          }
        }
    }
    else if( gausstype == "MomentFitting" || gausstype=="DirectDivergence" )
    {
        Teuchos::RCP<DRT::UTILS::GaussPoints> gp = vc->GetGaussRule();
        gpc->Append(gp);
    }
  }

  return gpc;
}


#if 0
template <DRT::Element::DiscretizationType distype>
Teuchos::RCP<DRT::UTILS::GaussPoints> GEO::CUT::ElementHandle::TransformVolumeCellGaussPoints( GEO::CUT::IntegrationCell * ic )
{

  const int nsd = DRT::UTILS::DisTypeToDim<distype>::dim;
  const int nen = DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;

  int degree = ic->CubatureDegree( Shape() );

  Teuchos::RCP<DRT::UTILS::GaussPoints> gp = DRT::UTILS::GaussPointCache::Instance().Create( distype, degree );
  Teuchos::RCP<DRT::UTILS::CollectedGaussPoints> cgp = Teuchos::rcp( new DRT::UTILS::CollectedGaussPoints( gp->NumPoints() ) );

  LINALG::Matrix<nen,1> funct;
  LINALG::Matrix<nsd,nen> deriv;

  LINALG::Matrix<nsd,nsd> xjm;
//  LINALG::Matrix<nsd,1> xi;

  DRT::UTILS::GaussIntegration intpoints( gp );
  for ( DRT::UTILS::GaussIntegration::iterator iquad=intpoints.begin();
        iquad!=intpoints.end();
        ++iquad )
  {
    LINALG::Matrix<nsd,1> eta( iquad.Point() );

    // cell shape functions and their first derivatives
    DRT::UTILS::shape_function<distype>(eta,funct);
    DRT::UTILS::shape_function_deriv1<distype>(eta,deriv);

    // local coordinates of gauss point
//    xi.Multiply( xie, funct );

    // get transposed of the jacobian matrix d x / d \xi
    // xjm(i,j) = deriv(i,k)*xyze(j,k)
//    const Epetra_SerialDenseMatrix & xyz_cell = ic->Coordinates();


      LINALG::Matrix<3, nen> xyz (true);

      const std::vector<GEO::CUT::Point*> & cpoints = ic->Points();
      if ( (int)cpoints.size() != nen )
        throw std::runtime_error( "non-matching number of points" );

      for ( int i=0; i<nen; ++i )
      {
        GEO::CUT::Point * p = cpoints[i];
        LINALG::Matrix<3,1> xyz_p (true);
        p->Coordinates(xyz_p.A());

//        std::cout << "coordinates of p: " << xyz_p << std::endl;

        for(int r=0; r<3; r++)
        {
        	xyz(r,i) = xyz_p(r);
        }
      }

//      std::cout << "xyz " << xyz << std::endl;

//    LINALG::Matrix<nsd,nen> xyz( xyz_cell );

//    LINALG::Matrix<nsd,nen> xyz_cell( ic->Coordinates() );

    xjm.MultiplyNT( deriv, xyz );
//    xjm.MultiplyNT( deriv, xie );

    // global coordinates of Gaussian point
    LINALG::Matrix<3,1> xyz_gp (true);
    xyz_gp.Multiply( xyz , funct );

    LINALG::Matrix<3,1> xi (true); // (true);
    LocalCoordinates( xyz_gp, xi );

//    std::cout << "xi: " << xi << std::endl;

    double det = xjm.Determinant();
if (det < 0.0 ) dserror( "det is negative");

double det_ele =0.0;
// additional transformation for element
{

    //! node coordinates
    LINALG::Matrix<3,20> xyze_;
    //! array for shape functions
    LINALG::Matrix<20,1> funct_;
    //! array for shape function derivatives w.r.t r,s,t
    LINALG::Matrix<3,20> deriv_;
    //! array for second derivatives of shape function w.r.t r,s,t
    LINALG::Matrix<6,20> deriv2_;
    //! transposed jacobian "dx/ds"
    LINALG::Matrix<3,3> xjm_;
    //! inverse of transposed jacobian "ds/dx"
    LINALG::Matrix<3,3> xji_;
    //! global velocity derivatives in gausspoint w.r.t x,y,z

    std::vector<Node* > nodes = Nodes();
    if(nodes.size() != 20) dserror("not 20 nodes!!!");

    for(int i=0; i< (int)nodes.size(); i++)
    {
    	LINALG::Matrix<3,1> n_coord;
    	nodes[i]->Coordinates(n_coord.A());

    	for(int r=0; r< 3; r++) xyze_(r,i) = n_coord(r);
    }

  // shape functions and their first derivatives
  DRT::UTILS::shape_function<DRT::Element::hex20>(xi,funct_);
  DRT::UTILS::shape_function_deriv1<DRT::Element::hex20>(xi,deriv_);

  if (true)
  {
    // get the second derivatives of standard element at current GP
    DRT::UTILS::shape_function_deriv2<DRT::Element::hex20>(xi,deriv2_);
  }

  xjm_.MultiplyNT(deriv_,xyze_);
  det_ele = xji_.Invert(xjm_);

//  if (det_ < 1E-16)
//    dserror("GLOBAL ELEMENT NO.%i\nZERO OR NEGATIVE JACOBIAN DETERMINANT: %f", eleid, det_);
//
//  // compute integration factor
//  fac_ = iquad.Weight()*det_;
//
//  // compute global first derivates
//  derxy_.Multiply(xji_,deriv_);

}




    cgp->Append( xi , iquad.Weight()*det / det_ele);
  }
  return cgp;



}



void GEO::CUT::ElementHandle::VolumeCellGaussPointsLin( plain_volumecell_set & cells,
                                                     std::vector<DRT::UTILS::GaussIntegration> & intpoints )
{
  GetVolumeCells( cells );

  intpoints.clear();
  intpoints.reserve( cells.size() );

  for ( plain_volumecell_set::iterator i=cells.begin(); i!=cells.end(); ++i )
  {
    GEO::CUT::VolumeCell * vc = *i;


//    if(vc->Position() == GEO::CUT::Point::outside) std::cout << "outside vc" << std::endl;
//    else std::cout << "inside vc" << std::endl;

    Teuchos::RCP<DRT::UTILS::GaussPointsComposite> gpc =
      Teuchos::rcp( new DRT::UTILS::GaussPointsComposite( 0 ) );

    const plain_integrationcell_set & cells = vc->IntegrationCells();
    for ( plain_integrationcell_set::const_iterator i=cells.begin(); i!=cells.end(); ++i )
    {
      GEO::CUT::IntegrationCell * ic = *i;
      switch ( ic->Shape() )
      {
      case DRT::Element::hex8:
      {
        Teuchos::RCP<DRT::UTILS::GaussPoints> gp = TransformVolumeCellGaussPoints<DRT::Element::hex8>( ic );
        gpc->Append( gp );
        break;
      }
      case DRT::Element::tet4:
      {
        Teuchos::RCP<DRT::UTILS::GaussPoints> gp = TransformVolumeCellGaussPoints<DRT::Element::tet4>( ic );
        gpc->Append( gp );
        break;
      }
      case DRT::Element::wedge6:
      {
        Teuchos::RCP<DRT::UTILS::GaussPoints> gp = TransformVolumeCellGaussPoints<DRT::Element::wedge6>( ic );
        gpc->Append( gp );
        break;
      }
      case DRT::Element::pyramid5:
      {
        Teuchos::RCP<DRT::UTILS::GaussPoints> gp = TransformVolumeCellGaussPoints<DRT::Element::pyramid5>( ic );
        gpc->Append( gp );
        break;
      }
      default:
        throw std::runtime_error( "unsupported integration cell type" );
      }
    }

    intpoints.push_back( DRT::UTILS::GaussIntegration( gpc ) );
  }


}

#endif

void GEO::CUT::ElementHandle::BoundaryCellGaussPoints( const std::map<int, std::vector<GEO::CUT::BoundaryCell*> > & bcells,
                                                       std::map<int, std::vector<DRT::UTILS::GaussIntegration> > & intpoints )
{
  for ( std::map<int, std::vector<GEO::CUT::BoundaryCell*> >::const_iterator i=bcells.begin(); i!=bcells.end(); ++i )
  {
    int sid = i->first;
    const std::vector<GEO::CUT::BoundaryCell*> & cells = i->second;
    std::vector<DRT::UTILS::GaussIntegration> & cell_points = intpoints[sid];

    cell_points.clear();
    cell_points.reserve( cells.size() );

    for ( std::vector<GEO::CUT::BoundaryCell*>::const_iterator i=cells.begin(); i!=cells.end(); ++i )
    {
      GEO::CUT::BoundaryCell * bc = *i;

      switch ( bc->Shape() )
      {
      case DRT::Element::tri3:
      {
        cell_points.push_back( DRT::UTILS::GaussIntegration( DRT::Element::tri3,
                                                             bc->CubatureDegree( DRT::Element::quad9 ) ) );
        break;
      }
      case DRT::Element::quad4:
      {
        cell_points.push_back( DRT::UTILS::GaussIntegration( DRT::Element::quad4,
                                                             bc->CubatureDegree( DRT::Element::quad9 ) ) );
        break;
      }
      default:
        throw std::runtime_error( "unsupported integration cell type" );
      }
    }
  }
}

void GEO::CUT::ElementHandle::BoundaryCellGaussPoints( MeshIntersection & mesh,
                                                       int mi,
                                                       const std::map<int, std::vector<GEO::CUT::BoundaryCell*> > & bcells,
                                                       std::map<int, std::vector<DRT::UTILS::GaussIntegration> > & intpoints )
{
  for ( std::map<int, std::vector<GEO::CUT::BoundaryCell*> >::const_iterator i=bcells.begin(); i!=bcells.end(); ++i )
  {
    int sid = i->first;
    const std::vector<GEO::CUT::BoundaryCell*> & cells = i->second;
    std::vector<DRT::UTILS::GaussIntegration> & cell_points = intpoints[sid];

    SideHandle * side = mesh.GetCutSide( sid, mi );
    if ( side==NULL )
    {
      throw std::runtime_error( "no such side" );
    }

    cell_points.clear();
    cell_points.reserve( cells.size() );

    for ( std::vector<GEO::CUT::BoundaryCell*>::const_iterator i=cells.begin(); i!=cells.end(); ++i )
    {
      GEO::CUT::BoundaryCell * bc = *i;

      switch ( bc->Shape() )
      {
      case DRT::Element::tri3:
      {
        Teuchos::RCP<DRT::UTILS::GaussPoints> gp = side->CreateProjected<DRT::Element::tri3>( bc );
        cell_points.push_back( DRT::UTILS::GaussIntegration( gp ) );
        break;
      }
      case DRT::Element::quad4:
      {
        Teuchos::RCP<DRT::UTILS::GaussPoints> gp = side->CreateProjected<DRT::Element::quad4>( bc );
        cell_points.push_back( DRT::UTILS::GaussIntegration( gp ) );
        break;
      }
      default:
        throw std::runtime_error( "unsupported integration cell type" );
      }
    }
  }
}


void GEO::CUT::ElementHandle::BoundaryCellGaussPointsLin( MeshIntersection & mesh,
                                                       int mi,
                                                       const std::map<int, std::vector<GEO::CUT::BoundaryCell*> > & bcells,
                                                       std::map<int, std::vector<DRT::UTILS::GaussIntegration> > & intpoints )
{
  for ( std::map<int, std::vector<GEO::CUT::BoundaryCell*> >::const_iterator i=bcells.begin(); i!=bcells.end(); ++i )
  {
    int sid = i->first;
    const std::vector<GEO::CUT::BoundaryCell*> & cells = i->second;
    std::vector<DRT::UTILS::GaussIntegration> & cell_points = intpoints[sid];

    SideHandle * side = mesh.GetCutSide( sid, mi );
    if ( side==NULL )
    {
      throw std::runtime_error( "no such side" );
    }

    cell_points.clear();
    cell_points.reserve( cells.size() );

    for ( std::vector<GEO::CUT::BoundaryCell*>::const_iterator i=cells.begin(); i!=cells.end(); ++i )
    {
      GEO::CUT::BoundaryCell * bc = *i;

      // Create (unmodified) gauss points for integration cell with requested
      // polynomial order. This is supposed to be fast, since there is a cache.
      cell_points.push_back(bc->gaussRule());
    }
  }

/*  static int eeno=0;
  eeno++;
  if(eeno<5)
  {
    std::string filename="wrong";
    std::ofstream file;

    std::stringstream out;
    out <<"boundary"<<eeno<<".dat";
    filename = out.str();
    file.open(filename.c_str());
    for ( std::map<int, std::vector<GEO::CUT::BoundaryCell*> >::const_iterator i=bcells.begin(); i!=bcells.end(); ++i )
    {
      int sid = i->first;
      std::vector<DRT::UTILS::GaussIntegration> ga1 = intpoints[sid];

      const std::vector<GEO::CUT::BoundaryCell*> cells = i->second;
      int m=0;
      for ( std::vector<GEO::CUT::BoundaryCell*>::const_iterator i=cells.begin(); i!=cells.end(); ++i )
      {
        GEO::CUT::BoundaryCell * bc = *i;
        file<<"boundary cell type = "<<bc->Shape()<<"\n";
        DRT::UTILS::GaussIntegration ga = ga1[m];
        for ( DRT::UTILS::GaussIntegration::const_iterator iquad=ga.begin(); iquad!=ga.end(); ++iquad )
        {
          const double* gpp = iquad.Point();
          file<<gpp[0]<<"\t"<<gpp[1]<<"\t"<<gpp[2]<<"\t";
          file<<iquad.Weight()<<std::endl;
        }
        m++;
      }

    }
    file.close();
  }*/
}

void GEO::CUT::ElementHandle::BoundaryCellGaussPointsLevelset( LevelSetIntersection & mesh,
                                                               const std::map<int, std::vector<GEO::CUT::BoundaryCell*> > & bcells,
                                                               std::map<int, std::vector<DRT::UTILS::GaussIntegration> > & intpoints )
{
#if 1
  dserror( "todo" );
#else
  for ( std::map<int, std::vector<GEO::CUT::BoundaryCell*> >::const_iterator i=bcells.begin(); i!=bcells.end(); ++i )
  {
    int sid = i->first;
    const std::vector<GEO::CUT::BoundaryCell*> & cells = i->second;
    std::vector<DRT::UTILS::GaussIntegration> & cell_points = intpoints[sid];

    SideHandle * side = mesh.GetCutSide( sid, mi );
    if ( side==NULL )
    {
      throw std::runtime_error( "no such side" );
    }

    cell_points.clear();
    cell_points.reserve( cells.size() );

    for ( std::vector<GEO::CUT::BoundaryCell*>::const_iterator i=cells.begin(); i!=cells.end(); ++i )
    {
      GEO::CUT::BoundaryCell * bc = *i;

      switch ( bc->Shape() )
      {
      case DRT::Element::tri3:
      {
        Teuchos::RCP<DRT::UTILS::GaussPoints> gp = side->CreateProjected<DRT::Element::tri3>( bc );
        cell_points.push_back( DRT::UTILS::GaussIntegration( gp ) );
        break;
      }
      case DRT::Element::quad4:
      {
        Teuchos::RCP<DRT::UTILS::GaussPoints> gp = side->CreateProjected<DRT::Element::quad4>( bc );
        cell_points.push_back( DRT::UTILS::GaussIntegration( gp ) );
        break;
      }
      default:
        throw std::runtime_error( "unsupported integration cell type" );
      }
    }
  }
#endif
}

void GEO::CUT::LinearElementHandle::GetVolumeCells( plain_volumecell_set & cells )
{
  const plain_volumecell_set & cs = element_->VolumeCells();
  std::copy( cs.begin(), cs.end(), std::inserter( cells, cells.begin() ) );
}


void GEO::CUT::LinearElementHandle::GetCellSets_DofSets_GaussPoints ( std::vector<plain_volumecell_set > & cell_sets ,
                                                                      std::vector< std::vector< int > >  & nds_sets,
                                                                      std::vector<std::vector< DRT::UTILS::GaussIntegration > > & intpoints_sets,
                                                                      std::string gausstype)
{
    TEUCHOS_FUNC_TIME_MONITOR( "FLD::XFluid::XFluidState::Get_CellSets_nds_GaussPoints" );

    GetVolumeCellsDofSets( cell_sets, nds_sets );

    intpoints_sets.clear();
    //intpoints_sets.reserve( cell_sets.size() );

    for(std::vector<plain_volumecell_set>::iterator i= cell_sets.begin(); i!=cell_sets.end(); i++)
    {
        plain_volumecell_set & cells = *i;

        if(cells.size() != 1) dserror("number of volumecells in set not equal 1, this should not be for linear elements!");

        std::vector<DRT::UTILS::GaussIntegration> gaussCellsets;
        VolumeCellGaussPoints( cells, gaussCellsets, gausstype );

        intpoints_sets.push_back( gaussCellsets );
    }
}


void GEO::CUT::LinearElementHandle::CollectVolumeCells ( bool include_inner, plain_volumecell_set & cells_inside, plain_volumecell_set & cells_outside )
{

	    const plain_volumecell_set & ecells = element_->VolumeCells();

	    // sort for inside and outside volume cells
	    for (plain_volumecell_set::const_iterator i = ecells.begin(); i!=ecells.end(); i++)
	    {
	    	if( (*i)->Position() == GEO::CUT::Point::outside )
	    	{
	    		cells_outside.insert( *i );
	    	}
	    	else // inside vc
	    	{
	    		if( include_inner)
	    		{
	    			cells_inside.insert( *i );
	    		}
	    	}
	    }

}


void GEO::CUT::LinearElementHandle::VolumeCellSets ( bool include_inner, std::vector<plain_volumecell_set> & ele_vc_sets_inside, std::vector<plain_volumecell_set> & ele_vc_sets_outside)
{

	if(!cells_set_)
	{
	    const plain_volumecell_set & ecells = element_->VolumeCells();

	    // sort for inside and outside volume cells
	    for (plain_volumecell_set::const_iterator i = ecells.begin(); i!=ecells.end(); i++)
	    {
	    	if( (*i)->Position() == GEO::CUT::Point::outside )
	    	{
	    		plain_volumecell_set s;
	    		s.insert(*i);
	    		vc_sets_outside_.push_back(s);
	    	}
	    	else // inside vc
	    	{
	    		if( include_inner)
	    		{
	    			plain_volumecell_set s;
	    			s.insert(*i);
	    			vc_sets_inside_.push_back(s);
	    		}
	    	}
	    }

	    cells_set_ = true;
	}



	if(include_inner)
	{
	std::copy( vc_sets_inside_.begin(), vc_sets_inside_.end(), std::inserter(ele_vc_sets_inside, ele_vc_sets_inside.begin()));
	}

	std::copy( vc_sets_outside_.begin(), vc_sets_outside_.end(), std::inserter(ele_vc_sets_outside, ele_vc_sets_outside.begin()));

}

//void GEO::CUT::LinearElementHandle::GetNodalDofSets ( bool include_inner, std::vector<std::vector<int> > & nodaldofset_vc_sets_inside, std::vector<std::vector<int> > & nodaldofset_vc_sets_outside)
//{
//	if(include_inner)
//	{
//		nodaldofset_vc_sets_inside = nodaldofset_vc_sets_inside_;
//	}
//	nodaldofset_vc_sets_outside = nodaldofset_vc_sets_outside_;
//}

void GEO::CUT::QuadraticElementHandle::BuildCellSets ( plain_volumecell_set & cells_to_connect, std::vector<plain_volumecell_set> & connected_sets)
{
	plain_volumecell_set done;

	  for ( plain_volumecell_set::const_iterator i=cells_to_connect.begin();
	        i!=cells_to_connect.end();
	        ++i )
	  {
	    VolumeCell * cell = *i;
	    if ( done.count( cell )==0 ) // cell currently not-done
	    {
	      plain_volumecell_set connected;
	      // REMARK: here use the version without! elements check:
	      // here we build cell sets within one global element with vcs of subelements
	      // maybe the vcs of one subelement are not connected within one subelement, but within one global element,
	      // therefore more than one vc of one subelements may be connected.
//	      cell->Neighbors( NULL, cells_to_connect, done, connected, elements );
	      cell->Neighbors( NULL, cells_to_connect, done, connected);

	      if ( connected.size()>0 )
	      {
        	connected_sets.push_back( connected );
	        std::copy( connected.begin(), connected.end(), std::inserter( done, done.begin() ) );
	      }
	    }
	  }
}



bool GEO::CUT::QuadraticElementHandle::IsCut()
{
  for ( std::vector<Element*>::iterator i=subelements_.begin(); i!=subelements_.end(); ++i )
  {
    Element * e = *i;
    if ( e->IsCut() )
    {
      return true;
    }
  }
  return false;
}

void GEO::CUT::QuadraticElementHandle::GetVolumeCells( plain_volumecell_set & cells )
{
  for ( std::vector<Element*>::iterator i=subelements_.begin(); i!=subelements_.end(); ++i )
  {
    Element * e = *i;
    const plain_volumecell_set & cs = e->VolumeCells();
    std::copy( cs.begin(), cs.end(), std::inserter( cells, cells.begin() ) );
  }
}

void GEO::CUT::QuadraticElementHandle::GetIntegrationCells( plain_integrationcell_set & cells )
{
  for ( std::vector<Element*>::iterator i=subelements_.begin(); i!=subelements_.end(); ++i )
  {
    Element * e = *i;
    e->GetIntegrationCells( cells );
  }
}

void GEO::CUT::QuadraticElementHandle::GetBoundaryCells( plain_boundarycell_set & bcells )
{
  for ( std::vector<Element*>::iterator i=subelements_.begin(); i!=subelements_.end(); ++i )
  {
    Element * e = *i;
    e->GetBoundaryCells( bcells );
  }
}

void GEO::CUT::QuadraticElementHandle::VolumeCells( plain_volumecell_set & cells ) const
{
  for ( std::vector<Element*>::const_iterator i=subelements_.begin(); i!=subelements_.end(); ++i )
  {
    Element * e = *i;
    const plain_volumecell_set & ecells = e->VolumeCells();
    std::copy( ecells.begin(), ecells.end(), std::inserter( cells, cells.begin() ) );
  }
}


void GEO::CUT::QuadraticElementHandle::GetCellSets_DofSets_GaussPoints ( std::vector<plain_volumecell_set > & cell_sets ,
                                                               std::vector< std::vector< int > >  & nds_sets,
                                                               std::vector<std::vector< DRT::UTILS::GaussIntegration > > & intpoints_sets,
                                                               std::string gausstype)
{
    GetVolumeCellsDofSets( cell_sets, nds_sets );

    intpoints_sets.clear();
    intpoints_sets.reserve( cell_sets.size() );

    for(std::vector<plain_volumecell_set>::iterator i= cell_sets.begin(); i!=cell_sets.end(); i++)
    {
        plain_volumecell_set & cells = *i;

        std::vector<DRT::UTILS::GaussIntegration> gaussCellsets;
        VolumeCellGaussPoints( cells, gaussCellsets, gausstype );

        intpoints_sets.push_back( gaussCellsets );
    }

}



void GEO::CUT::QuadraticElementHandle::CollectVolumeCells ( bool include_inner, plain_volumecell_set & cells_inside, plain_volumecell_set & cells_outside )
{

	  for ( std::vector<Element*>::const_iterator i=subelements_.begin(); i!=subelements_.end(); ++i )
	  {
	    Element * e = *i;
	    const plain_volumecell_set & ecells = e->VolumeCells();

	    // sort for inside and outside volume cells
	    for (plain_volumecell_set::const_iterator i = ecells.begin(); i!=ecells.end(); i++)
	    {
	    	if( (*i)->Position() == GEO::CUT::Point::outside )
	    	{
	    		cells_outside.insert( *i );
	    	}
	    	else // inside vc
	    	{
	    		if( include_inner)
	    		{
	    			cells_inside.insert( *i );
	    		}
	    	}
	    }

	  }
}


void GEO::CUT::QuadraticElementHandle::VolumeCellSets ( bool include_inner, std::vector<plain_volumecell_set> & ele_vc_sets_inside, std::vector<plain_volumecell_set> & ele_vc_set_ouside)
{

	// connect volumecells of subelements
	ConnectVolumeCells( include_inner );

	// get inside and outside connected volumecellsets for the current element
	GetConnectedVolumeCellSets( include_inner, ele_vc_sets_inside, ele_vc_set_ouside );

}




void GEO::CUT::QuadraticElementHandle::GetConnectedVolumeCellSets ( bool include_inner, std::vector<plain_volumecell_set> & connected_vc_sets_inside, std::vector<plain_volumecell_set> & connected_vc_sets_outside)
{

	if(include_inner)
	{
		std::copy(connected_vc_sets_inside_.begin(), connected_vc_sets_inside_.end(), std::inserter(connected_vc_sets_inside, connected_vc_sets_inside.end()));
	}

	std::copy(connected_vc_sets_outside_.begin(), connected_vc_sets_outside_.end(), std::inserter(connected_vc_sets_outside, connected_vc_sets_outside.end()));
}



void GEO::CUT::QuadraticElementHandle::ConnectVolumeCells ( bool include_inner )
{
	// find the connection between volumecells of all subelements for the current element (hex8, hex20, tet4, tet10 etc.)
	// remark: this function determines not the connection outside this element

	// get the volumecells of all subelements stored in two plainvolume_sets (inside and outside)
	if( !cells_connected_ )
	{

	plain_volumecell_set e_vcs_inside;
	plain_volumecell_set e_vcs_outside;

	CollectVolumeCells(include_inner, e_vcs_inside, e_vcs_outside);

    if( include_inner )
    {
    	BuildCellSets(e_vcs_inside , connected_vc_sets_inside_);
    }
    BuildCellSets(e_vcs_outside , connected_vc_sets_outside_);

    cells_connected_ = true;

	}

}

GEO::CUT::Hex20ElementHandle::Hex20ElementHandle( Mesh & mesh, int eid, const std::vector<int> & nodes )
  : QuadraticElementHandle()
{
  subelements_.reserve( 8 );

  const CellTopologyData * top_data = shards::getCellTopologyData< shards::Hexahedron<8> >();

  // create middle nodes

  LINALG::Matrix<3,1> xyz;
  LINALG::Matrix<1,1> lsv;

  plain_int_set node_nids;

  LINALG::Matrix<3,8> side_xyze;
  LINALG::Matrix<1,8> side_lsvs;
  LINALG::Matrix<8,1> side_funct;
  std::vector<Node*> side_nodes( 8 );

  std::vector<Node*> center_nodes( 6 );

  for ( int localsideid = 0; localsideid < 6; ++localsideid )
  {
    node_nids.clear();
    for ( int i=0; i<8; ++i )
    {
      int localnodeid = DRT::UTILS::eleNodeNumbering_hex27_surfaces[localsideid][i];
      Node * n = mesh.GetNode( nodes[localnodeid], static_cast<double*>( NULL ) );
      side_nodes[i] = n;
      node_nids.insert( nodes[localnodeid] );
      n->Coordinates( &side_xyze( 0, i ) );
      side_lsvs( i ) = n->LSV();
    }

    DRT::UTILS::shape_function_2D( side_funct, 0, 0, DRT::Element::quad8 );
    xyz.Multiply( side_xyze, side_funct );
    lsv.Multiply( side_lsvs, side_funct );

    center_nodes[localsideid] = mesh.GetNode( node_nids, xyz.A(), lsv( 0 ) );
  }

  Node* node20 = center_nodes[0];
  int node20_id = node20->Id();

  Node* node21 = center_nodes[1];
  int node21_id = node21->Id();

  Node* node22 = center_nodes[2];
  int node22_id = node22->Id();

  Node* node23 = center_nodes[3];
  int node23_id = node23->Id();

  Node* node24 = center_nodes[4];
  int node24_id = node24->Id();

  Node* node25 = center_nodes[5];
  int node25_id = node25->Id();

  LINALG::Matrix<3,20> xyze;
  LINALG::Matrix<1,20> lsvs;
  nodes_.reserve( 20 );
  for ( int i=0; i<20; ++i )
  {
    Node * n = mesh.GetNode( nodes[i], static_cast<double*>( NULL ) );
    nodes_.push_back( n );
    n->Coordinates( &xyze( 0,i ) );
    lsvs( i ) = n->LSV();
  }

  LINALG::Matrix<20,1> funct;
  DRT::UTILS::shape_function_3D( funct, 0, 0, 0, DRT::Element::hex20 );

  xyz.Multiply( xyze, funct );
  lsv.Multiply( lsvs, funct );
  node_nids.clear();
  std::copy( nodes.begin(), nodes.end(), std::inserter( node_nids, node_nids.begin() ) );
  Node* node26 = mesh.GetNode( node_nids, xyz.A(), lsv( 0 ) );
  int node26_id = node26->Id();


  std::vector<int> nids( 8 );

  nids[0] = nodes[ 0];
  nids[1] = nodes[ 8];
  nids[2] = node20_id;
  nids[3] = nodes[11];
  nids[4] = nodes[12];
  nids[5] = node21_id;
  nids[6] = node26_id;
  nids[7] = node24_id;
  subelements_.push_back( mesh.GetElement( -1, nids, *top_data ) );

  nids[0] = nodes[ 8];
  nids[1] = nodes[ 1];
  nids[2] = nodes[ 9];
  nids[3] = node20_id;
  nids[4] = node21_id;
  nids[5] = nodes[13];
  nids[6] = node22_id;
  nids[7] = node26_id;
  subelements_.push_back( mesh.GetElement( -1, nids, *top_data ) );

  nids[0] = node20_id;
  nids[1] = nodes[ 9];
  nids[2] = nodes[ 2];
  nids[3] = nodes[10];
  nids[4] = node26_id;
  nids[5] = node22_id;
  nids[6] = nodes[14];
  nids[7] = node23_id;
  subelements_.push_back( mesh.GetElement( -1, nids, *top_data ) );

  nids[0] = nodes[11];
  nids[1] = node20_id;
  nids[2] = nodes[10];
  nids[3] = nodes[ 3];
  nids[4] = node24_id;
  nids[5] = node26_id;
  nids[6] = node23_id;
  nids[7] = nodes[15];
  subelements_.push_back( mesh.GetElement( -1, nids, *top_data ) );

  /////////////////////////////////////////////////////////////////

  nids[0] = nodes[12];
  nids[1] = node21_id;
  nids[2] = node26_id;
  nids[3] = node24_id;
  nids[4] = nodes[ 4];
  nids[5] = nodes[16];
  nids[6] = node25_id;
  nids[7] = nodes[19];
  subelements_.push_back( mesh.GetElement( -1, nids, *top_data ) );

  nids[0] = node21_id;
  nids[1] = nodes[13];
  nids[2] = node22_id;
  nids[3] = node26_id;
  nids[4] = nodes[16];
  nids[5] = nodes[ 5];
  nids[6] = nodes[17];
  nids[7] = node25_id;
  subelements_.push_back( mesh.GetElement( -1, nids, *top_data ) );

  nids[0] = node26_id;
  nids[1] = node22_id;
  nids[2] = nodes[14];
  nids[3] = node23_id;
  nids[4] = node25_id;
  nids[5] = nodes[17];
  nids[6] = nodes[ 6];
  nids[7] = nodes[18];
  subelements_.push_back( mesh.GetElement( -1, nids, *top_data ) );

  nids[0] = node24_id;
  nids[1] = node26_id;
  nids[2] = node23_id;
  nids[3] = nodes[15];
  nids[4] = nodes[19];
  nids[5] = node25_id;
  nids[6] = nodes[18];
  nids[7] = nodes[ 7];
  subelements_.push_back( mesh.GetElement( -1, nids, *top_data ) );
}

GEO::CUT::Hex27ElementHandle::Hex27ElementHandle( Mesh & mesh, int eid, const std::vector<int> & nodes )
  : QuadraticElementHandle()
{
  subelements_.reserve( 8 );

  nodes_.reserve( 27 );
  for ( int i=0; i<27; ++i )
  {
    Node * n = mesh.GetNode( nodes[i], static_cast<double*>( NULL ) );
    nodes_.push_back( n );
  }

  const CellTopologyData * top_data = shards::getCellTopologyData< shards::Hexahedron<8> >();

  std::vector<int> nids( 8 );

  nids[0] = nodes[ 0];
  nids[1] = nodes[ 8];
  nids[2] = nodes[20];
  nids[3] = nodes[11];
  nids[4] = nodes[12];
  nids[5] = nodes[21];
  nids[6] = nodes[26];
  nids[7] = nodes[24];
  subelements_.push_back( mesh.GetElement( -1, nids, *top_data ) );

  nids[0] = nodes[ 8];
  nids[1] = nodes[ 1];
  nids[2] = nodes[ 9];
  nids[3] = nodes[20];
  nids[4] = nodes[21];
  nids[5] = nodes[13];
  nids[6] = nodes[22];
  nids[7] = nodes[26];
  subelements_.push_back( mesh.GetElement( -1, nids, *top_data ) );

  nids[0] = nodes[20];
  nids[1] = nodes[ 9];
  nids[2] = nodes[ 2];
  nids[3] = nodes[10];
  nids[4] = nodes[26];
  nids[5] = nodes[22];
  nids[6] = nodes[14];
  nids[7] = nodes[23];
  subelements_.push_back( mesh.GetElement( -1, nids, *top_data ) );

  nids[0] = nodes[11];
  nids[1] = nodes[20];
  nids[2] = nodes[10];
  nids[3] = nodes[ 3];
  nids[4] = nodes[24];
  nids[5] = nodes[26];
  nids[6] = nodes[23];
  nids[7] = nodes[15];
  subelements_.push_back( mesh.GetElement( -1, nids, *top_data ) );

  /////////////////////////////////////////////////////////////////

  nids[0] = nodes[12];
  nids[1] = nodes[21];
  nids[2] = nodes[26];
  nids[3] = nodes[24];
  nids[4] = nodes[ 4];
  nids[5] = nodes[16];
  nids[6] = nodes[25];
  nids[7] = nodes[19];
  subelements_.push_back( mesh.GetElement( -1, nids, *top_data ) );

  nids[0] = nodes[21];
  nids[1] = nodes[13];
  nids[2] = nodes[22];
  nids[3] = nodes[26];
  nids[4] = nodes[16];
  nids[5] = nodes[ 5];
  nids[6] = nodes[17];
  nids[7] = nodes[25];
  subelements_.push_back( mesh.GetElement( -1, nids, *top_data ) );

  nids[0] = nodes[26];
  nids[1] = nodes[22];
  nids[2] = nodes[14];
  nids[3] = nodes[23];
  nids[4] = nodes[25];
  nids[5] = nodes[17];
  nids[6] = nodes[ 6];
  nids[7] = nodes[18];
  subelements_.push_back( mesh.GetElement( -1, nids, *top_data ) );

  nids[0] = nodes[24];
  nids[1] = nodes[26];
  nids[2] = nodes[23];
  nids[3] = nodes[15];
  nids[4] = nodes[19];
  nids[5] = nodes[25];
  nids[6] = nodes[18];
  nids[7] = nodes[ 7];
  subelements_.push_back( mesh.GetElement( -1, nids, *top_data ) );
}

GEO::CUT::Tet10ElementHandle::Tet10ElementHandle( Mesh & mesh, int eid, const std::vector<int> & nids )
  : QuadraticElementHandle()
{
  subelements_.reserve( 8 );

  nodes_.reserve( 10 );
  for ( int i=0; i<10; ++i )
  {
    Node * n = mesh.GetNode( nids[i], static_cast<double*>( NULL ) );
    nodes_.push_back( n );
  }

  const CellTopologyData * top_data = shards::getCellTopologyData< shards::Tetrahedron<4> >();

  std::vector<int> subnids( 4 );

  subnids[0] = nids[ 0];
  subnids[1] = nids[ 4];
  subnids[2] = nids[ 6];
  subnids[3] = nids[ 7];
  subelements_.push_back( mesh.GetElement( -1, subnids, *top_data ) );

  subnids[0] = nids[ 4];
  subnids[1] = nids[ 1];
  subnids[2] = nids[ 5];
  subnids[3] = nids[ 8];
  subelements_.push_back( mesh.GetElement( -1, subnids, *top_data ) );

  subnids[0] = nids[ 6];
  subnids[1] = nids[ 5];
  subnids[2] = nids[ 2];
  subnids[3] = nids[ 9];
  subelements_.push_back( mesh.GetElement( -1, subnids, *top_data ) );

  subnids[0] = nids[ 7];
  subnids[1] = nids[ 8];
  subnids[2] = nids[ 9];
  subnids[3] = nids[ 3];
  subelements_.push_back( mesh.GetElement( -1, subnids, *top_data ) );

  /////////////////////////////////////////////////////////////////

  subnids[0] = nids[ 4];
  subnids[1] = nids[ 5];
  subnids[2] = nids[ 6];
  subnids[3] = nids[ 8];
  subelements_.push_back( mesh.GetElement( -1, subnids, *top_data ) );

  subnids[0] = nids[ 6];
  subnids[1] = nids[ 9];
  subnids[2] = nids[ 7];
  subnids[3] = nids[ 8];
  subelements_.push_back( mesh.GetElement( -1, subnids, *top_data ) );

  subnids[0] = nids[ 4];
  subnids[1] = nids[ 6];
  subnids[2] = nids[ 7];
  subnids[3] = nids[ 8];
  subelements_.push_back( mesh.GetElement( -1, subnids, *top_data ) );

  subnids[0] = nids[ 9];
  subnids[1] = nids[ 6];
  subnids[2] = nids[ 5];
  subnids[3] = nids[ 8];
  subelements_.push_back( mesh.GetElement( -1, subnids, *top_data ) );
}

void GEO::CUT::Hex20ElementHandle::LocalCoordinates( const LINALG::Matrix<3,1> & xyz, LINALG::Matrix<3,1> & rst )
{
  LINALG::Matrix<3, 20> xyze;

  for ( int i=0; i<20; ++i )
  {
    Node * n = nodes_[i];
    n->Coordinates( &xyze( 0, i ) );
  }

  Position<DRT::Element::hex20> pos( xyze, xyz );
  bool success = pos.Compute();
  if ( not success )
  {
	  dserror("local coordinates for hex20 element could not be determined");
  }
  rst = pos.LocalCoordinates();
}

void GEO::CUT::Hex27ElementHandle::LocalCoordinates( const LINALG::Matrix<3,1> & xyz, LINALG::Matrix<3,1> & rst )
{
  LINALG::Matrix<3, 27> xyze;

  for ( int i=0; i<27; ++i )
  {
    Node * n = nodes_[i];
    n->Coordinates( &xyze( 0, i ) );
  }

  Position<DRT::Element::hex27> pos( xyze, xyz );
  bool success = pos.Compute();
  if ( not success )
  {
  }
  rst = pos.LocalCoordinates();
}

void GEO::CUT::Tet10ElementHandle::LocalCoordinates( const LINALG::Matrix<3,1> & xyz, LINALG::Matrix<3,1> & rst )
{
  LINALG::Matrix<3, 10> xyze;

  for ( int i=0; i<10; ++i )
  {
    Node * n = nodes_[i];
    n->Coordinates( &xyze( 0, i ) );
  }

  Position<DRT::Element::tet10> pos( xyze, xyz );
  bool success = pos.Compute();
  if ( not success )
  {
  }
  rst = pos.LocalCoordinates();
}

