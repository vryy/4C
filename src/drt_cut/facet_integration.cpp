/*---------------------------------------------------------------------*/
/*!
\file facet_integration.cpp

\brief Integrates base functions over the facet for both volumecell facets and for boundarycells
equations

\level 2

<pre>
\maintainer Christoph Ager
            ager@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>

*----------------------------------------------------------------------*/
#include "line_integration.H"
#include "boundarycell_integration.H"

#include "cut_boundarycell.H"
#include "cut_triangulateFacet.H"
#include "cut_kernel.H"
#include "cut_side.H"

#include "cut_options.H"

#ifdef DEBUGCUTLIBRARY
#include "cut_output.H"
#endif

#include "../drt_lib/drt_colors.H"

#include <Teuchos_TimeMonitor.hpp>


/*-----------------------------------------------------------------------------------------------------*
      compute the equation of the plane Ax+By+Cz=D with the local coordinates of corner points
*------------------------------------------------------------------------------------------------------*/
std::vector<double> GEO::CUT::FacetIntegration::equation_plane(const std::vector<std::vector<double> > & cornersLocal)
{
  //TODO: use references for return!!!
#if 1  //Newell's method of determining equation of plane
  std::vector<double> eqn_plane = KERNEL::EqnPlaneOfPolygon( cornersLocal );
#endif


#if 0 // old method of deleting the inlin points etc etc..
  std::vector<double> eqn_plane(4);
  int cornSize = cornersLocal.size();

  // construct a temporary point list to find the concave points of the facet in local coords
  std::vector<GEO::CUT::Point*> ptlist(cornersLocal.size());
  GEO::CUT::Options options;
  GEO::CUT::Mesh mesh(options);

  for( int i=0;i<cornSize;i++ )
  {
    std::vector<double> coordPt = cornersLocal[i];
    ptlist[i] = mesh.NewPoint( &coordPt[0], NULL, NULL );
  }

  CUT::KERNEL::DeleteInlinePts( ptlist );

  // For some geometries, cut produces a facet with all points
  // lying on a line. all coefficients are zero in such cases
  if( ptlist.size()==0 )
  {
#if DEBUGCUTLIBRARY
    std::cout<<"WARNING:::cut algorithm produced a facet with all points on a line\n";
#endif
    for( unsigned i=0;i<4;i++ )
      eqn_plane[i] = 0.0;
  }
  else
  {
    eqn_plane = GEO::CUT::KERNEL::EqnPlanePolygon( ptlist );
  }
#endif

  return eqn_plane;
}

/*---------------------------------------------------------------------------------------------------------------------*
  A facet is clockwise ordered meaning that the normal vector of the facet is acting on the wrong direction.
  This routine checks if the facet is ordered clockwise                                                        sudhakar 05/15

  We have two cases:

  Case 1 : the facet is formed from a cut side
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      We compute the normal vectors of the parent side and the facet. If they are not acting on the opposite
      directions, then the facet is clockwise ordered.
      REMEMBER: This idea cannot be used in the case 2 below because of the use of Shards element topology, wherein
      all the surfaces are not ensured to have correct normal

  Case 2 : the facet is formed from the element side
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      Construct a reference vector from pointing from the element centre to the geometric centre of the facet.
      This vector must be in the correct normal direction. Take the dot product of reference vector and the
      normal vector of facet. If the dot product < 0, then the facet is clockwise ordered
*----------------------------------------------------------------------------------------------------------------------*/
void GEO::CUT::FacetIntegration::IsClockwise( const std::vector<double> & eqn_plane,
                                              const std::vector<std::vector<double> > & cornersLocal )
{
  orderingComputed_ = true;
  clockwise_ = false;

#if 1
  bool iscut = face1_->OnCutSide();

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Case 1: Facet is formed from the cut side
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if(iscut)
  {
    Side* parent = face1_->ParentSide();

    double dotProduct = 0.0;

    if(not face1_->BelongsToLevelSetSide())
    {
      const std::vector<Node*> &par_nodes = parent->Nodes();
      std::vector<std::vector<double> > corners(par_nodes.size());
      int mm=0;
      for(std::vector<Node*>::const_iterator i=par_nodes.begin();i!=par_nodes.end();i++)
      {
        Node *nod = *i;
        double x1[3];
        nod->Coordinates(x1);

        std::vector<double> pt_local(3);
#ifdef LOCAL
        LINALG::Matrix<3,1> glo,loc;

        for(int nodno=0;nodno<3;nodno++)
          glo(nodno,0) = x1[nodno];
        elem1_->LocalCoordinates(glo,loc);

        pt_local[0] = loc(0,0);
        pt_local[1] = loc(1,0);
        pt_local[2] = loc(2,0);
#else
        pt_local[0] = x1[0];
        pt_local[1] = x1[1];
        pt_local[2] = x1[2];
#endif

        corners[mm] = pt_local;
        mm++;
      }

      std::vector<double> eqn_par = equation_plane(corners);

#ifdef LOCAL
      dotProduct = eqn_plane[0]*eqn_par[0];
#else
     dotProduct = eqn_plane[0]*eqn_par[0] +
         eqn_plane[1]*eqn_par[1] +
         eqn_plane[2]*eqn_par[2];
#endif
    }
    else //If the facet is on a cut-side which also is a level set side.
    {
      std::vector<double> phi_deriv1;

      LINALG::Matrix<3,1> coord;

      //First entry is midpoint of triangulation!!
      //Third entry could possibly be a concave poin (take care of in Triangulation of Facet).
      // Thus-> Choose second entry, which should be "safe".
      coord(0,0)= cornersLocal[1][0];
      coord(1,0)= cornersLocal[1][1];
      coord(2,0)= cornersLocal[1][2];
#ifdef LOCAL
      phi_deriv1 = elem1_->GetLevelSetGradientAtLocalCoordsInLocalCoords(coord);
#else
      phi_deriv1 = elem1_->GetLevelSetGradient(coord, false);
#endif

      dotProduct = eqn_plane[0]*phi_deriv1[0] +
                   eqn_plane[1]*phi_deriv1[1] +
                   eqn_plane[2]*phi_deriv1[2];

#ifdef DIRECTDIV_EXTENDED_DEBUG_OUTPUT
      std::cout << "cornersLocalFacetForGrad: " << cornersLocal[1][0] << ", " << cornersLocal[1][1] << ", " << cornersLocal[1][2] << std::endl;
      std::cout << "phi_deriv1: " << phi_deriv1[0] << ", " << phi_deriv1[1] << ", " << phi_deriv1[2] << std::endl;
      std::cout << "eqn_plane_IC: " << eqn_plane[0] << ", " << eqn_plane[1] << ", " << eqn_plane[2] << std::endl;

      std::cout << "dotProduct: " << dotProduct << std::endl;
      std::cout << "position_ : " << position_ << std::endl;
#endif
    }

    // We use the fact that the normal to any cut side is always pointing towards fluid domain
    if( position_ == CUT::Point::inside)
    {
      if( dotProduct < 0.0 )
        clockwise_ = true;
    }
    else if (position_ == CUT::Point::outside)
    {
      if( dotProduct > 0.0 )
        clockwise_ = true;
    }
    else
    {
      dserror("VC position not defined!");
    }
  }

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Case 2: Facet is formed from the element  side
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  else
  {
    //-----
    // STEP 1: Get centroid of the parent element
    //-----
    LINALG::Matrix<3,1> elecen;
#ifdef LOCAL
    switch( elem1_->Shape() )
    {
    case DRT::Element::hex8:
    {
      elecen = DRT::UTILS::getLocalCenterPosition<3>(DRT::Element::hex8);
      break;
    }
    case DRT::Element::tet4:
    {
      elecen = DRT::UTILS::getLocalCenterPosition<3>(DRT::Element::tet4);
      break;
    }
    default:
    {
      dserror("Add other elements not only here but in complete direct divergence procedure");
      break;
    }
    }
#else
    elem1_->ElementCenter( elecen );
#endif

    //-----
    // STEP 2: Get geometric centre point of the facet
    // For concave facets, this is not the actual centre, but it does not matter
    //-----
    LINALG::Matrix<3,1> facecen;
    unsigned npts = cornersLocal.size();
    for( std::vector<std::vector<double> >::const_iterator fit = cornersLocal.begin(); fit != cornersLocal.end(); fit++ )
    {
      const std::vector<double> & fpt = *fit;
      for( unsigned dim = 0; dim < 3; dim++ )
        facecen(dim,0) += fpt[dim];
    }

    for( unsigned dim = 0; dim < 3; dim++ )
      facecen(dim,0) = facecen(dim,0) / npts;

    //-----
    // STEP 3: Construct a unit vector that points FROM element centre TO the facet centre
    // This reference vector is in the correct normal direction
    //-----
    LINALG::Matrix<3,1> ref_vec;
    ref_vec.Update( 1.0, facecen, -1.0, elecen );
    double l2_ref_vec = ref_vec.Norm2();
    ref_vec.Scale( 1.0/l2_ref_vec );

    //-----
    // STEP 4: Take dot product with the normal of facet
    // If both are in the opposite direction, then the facet nodes are arranged clockwise
    //-----
    LINALG::Matrix<3,1> norm_fac;
    for( unsigned dim = 0; dim < 3; dim++ )
      norm_fac(dim,0) = eqn_plane[dim];

    double dotProduct = ref_vec.Dot( norm_fac );
    if( dotProduct < 0.0 )
      clockwise_ = 1;
  }
  //std::cout<<"clockwise = "<<clockwise_<<"\t"<<"is cut side = "<<iscut<<"\n";
  return;

  // Old implementation of checking the normal direction
  // This involves the  use of element topology data extensively
  // TODO: Must go away after extensive testing of the above procedure
#else
  {
    std::string ordering;

    if(eqn_plane[0]>0.0)
      ordering = "acw";
    else
      ordering = "cw";

    clockwise_ = 0;
    Side* parent = face1_->ParentSide();
    const std::vector<Node*> &par_nodes = parent->Nodes();
    std::vector<std::vector<double> > corners(par_nodes.size());
    int mm=0;

    for(std::vector<Node*>::const_iterator i=par_nodes.begin();i!=par_nodes.end();i++)
    {
      Node *nod = *i;
      double x1[3];
      nod->Coordinates(x1);
      LINALG::Matrix<3,1> glo,loc;
      std::vector<double> pt_local(3);
      for(int nodno=0;nodno<3;nodno++)
        glo(nodno,0) = x1[nodno];

       elem1_->LocalCoordinates(glo,loc);

       pt_local[0] = loc(0,0);
       pt_local[1] = loc(1,0);
       pt_local[2] = loc(2,0);

       corners[mm] = pt_local;
       mm++;
      }

      std::vector<double> eqn_par = equation_plane(corners);

      // the cut sides are always convex
      // just comparing the sign of coordinates of eqn of plane will do
      // We use the fact that the normal to any cut side is always pointing towards fluid domain
      bool iscut = face1_->OnCutSide();
      if(iscut)
      {
        if( position_ == GEO::CUT::Point::inside )
        {
          if(eqn_plane[0]*eqn_par[0]<0.0)
          {
            clockwise_ = 1;
          }
        }
        else
        {
          if(eqn_plane[0]*eqn_par[0]>0.0)
          {
            clockwise_ = 1;
          }
        }
      }
      else
      {
        const std::vector<Side*> &ele_sides = elem1_->Sides();

        int parentSideno = 0;
        for(std::vector<Side*>::const_iterator i=ele_sides.begin();i!=ele_sides.end();i++)
        {
          Side*sss = *i;
          if(sss==parent)
          {
            break;
          }
          parentSideno++;

        }

  //should check whether this is sufficient or do we need to find the number of side in the element and
  //use their orientation to get clockwise ordering
  //                std::cout<<"parent side no = "<<parentSideno<<"\n";
  //      std::cout<<"parentSideno = "<<parentSideno<<"corner = "<<corners[0][0]<<"\n";
  //ParentSideno=1 is x=1 face and 3 is x=-1 face
      // After an element is mapped in local coordinate system, it has
      // non-zero x-normal component only for checked "parentSideno"
      // (See BaciReport for ordering of sides and nodes in Baci)
      switch ( elem1_->Shape() )
      {
        case DRT::Element::hex8:
        {
          if(parentSideno==0 and ordering=="cw")
            clockwise_ = 1;
          else if(parentSideno==1 and ordering=="cw")
            clockwise_ = 1;
          else if(parentSideno==2 and ordering=="acw")
            clockwise_ = 1;
          else if(parentSideno==3 and ordering=="acw")
            clockwise_ = 1;
          else if(parentSideno==4 and ordering=="acw")
            clockwise_ = 1;
          else if(parentSideno==5 and ordering=="cw")
            clockwise_ = 1;
          break;
        }
        case DRT::Element::tet4:
        {
          if(parentSideno==1 && ordering=="cw")
           clockwise_ = 1;
         if(parentSideno==2 && ordering=="acw")
           clockwise_ = 1;
          break;
        }
        case DRT::Element::wedge6:
        {
          if(parentSideno==0 && ordering=="cw")
           clockwise_ = 1;
         if(parentSideno==1 && ordering=="acw")
           clockwise_ = 1;
          break;
        }
        case DRT::Element::pyramid5:
        {
          if(parentSideno==1 && ordering=="cw")
            clockwise_ = 1;
          if(parentSideno==3 && ordering=="acw")
            clockwise_ = 1;
          break;
        }
        default:
          dserror( "unsupported integration cell type ( cell type = %s )",
              DRT::DistypeToString( elem1->Shape() ).c_str() );
          exit( EXIT_FAILURE );
      }
    }
    //std::cout<<"clockwise = "<<clockwise_<<"\t"<<"is cut side = "<<iscut<<"\n";
  }
#endif

}

/*-----------------------------------------------------------------------------------------------*
          Returns true if vertices of facets are ordered clockwise              sudhakar 07/12
*------------------------------------------------------------------------------------------------*/
bool GEO::CUT::FacetIntegration::IsClockwiseOrdering()
{
  if( orderingComputed_ )
    return clockwise_;

  std::vector<std::vector<double> > cornersLocal;
  face1_->CornerPointsLocal(elem1_,cornersLocal);
  if( eqn_plane_.size()==0)
  {
    eqn_plane_ = equation_plane(cornersLocal);
  }

  this->IsClockwise( eqn_plane_, cornersLocal );
  return clockwise_;
}

/*-----------------------------------------------------------------------------------------------*
                            computes x=f(y,z) from the plane equation
                  equation of this form is used to replace x in the line integral
*------------------------------------------------------------------------------------------------*/
std::vector<double> GEO::CUT::FacetIntegration::compute_alpha( std::vector<double> & eqn_plane,
                                                               GEO::CUT::ProjectionDirection intType )
{
  std::vector<double> alfa(3);
  double a = eqn_plane[0];
  double b = eqn_plane[1];
  double c = eqn_plane[2];
  double d = eqn_plane[3];


  if(intType==GEO::CUT::proj_x)
  {
    alfa[0] = d/a;
    alfa[1] = -1.0*b/a;
    alfa[2] = -1.0*c/a;
  }
  else if(intType==GEO::CUT::proj_y)
  {
    alfa[0] = d/b;
    alfa[1] = -1.0*c/b;
    alfa[2] = -1.0*a/b;
  }
  else if(intType==GEO::CUT::proj_z)
  {
    alfa[0] = d/c;
    alfa[1] = -1.0*a/c;
    alfa[2] = -1.0*b/c;
  }
  else
  {
    dserror("The facet integration type undefined");
    exit(1);
  }
  return alfa;
}

/*---------------------------------------------------------------------------------*
      return the absolute normal of the facet in a particular direction
*----------------------------------------------------------------------------------*/
double GEO::CUT::FacetIntegration::getNormal(GEO::CUT::ProjectionDirection intType)
{
  double normalScale=0.0;
  for(unsigned i=0;i<3;i++)
    normalScale += eqn_plane_[i]*eqn_plane_[i];
  normalScale = sqrt(normalScale);

  if(intType==GEO::CUT::proj_x)
    return eqn_plane_[0]/normalScale;
  else if(intType==GEO::CUT::proj_y)
    return eqn_plane_[1]/normalScale;
  else if(intType==GEO::CUT::proj_z)
    return eqn_plane_[2]/normalScale;
  else
  {
    dserror("The normal direction is unspecified");
    return 0.0;
  }
}

/*-------------------------------------------------------------------------------------*
                      perform integration over the facet
*--------------------------------------------------------------------------------------*/
double GEO::CUT::FacetIntegration::integrate_facet()
{
#ifdef LOCAL
  std::vector<std::vector<double> > cornersLocal;
  face1_->CornerPointsLocal(elem1_,cornersLocal,true);

  // Special case when using moment fitting for boundary integration
  if( global_ )
  {
    cornersLocal.clear();
    cornersLocal = face1_->CornerPointsGlobal(elem1_,true);
  }

#else
  std::vector<std::vector<double> > cornersLocal = face1_->CornerPointsGlobal(elem1_,true);
#endif

    eqn_plane_ = equation_plane(cornersLocal);

// the face is in the x-y or in y-z plane which gives zero facet integral
    if(fabs(eqn_plane_[0])<TOL_EQN_PLANE  && bcellInt_==false)
      return 0.0;
//x=0 plane which also do not contribute to facet integral
    if(fabs(eqn_plane_[1])<TOL_EQN_PLANE && fabs(eqn_plane_[2])<TOL_EQN_PLANE  &&
       fabs(eqn_plane_[3])<TOL_EQN_PLANE && bcellInt_==false)
      return 0.0;

    if(bcellInt_==true) // the integral value of boundarycell will not change w.r.t the ordering of vertices
      clockwise_ = false;
    else
      IsClockwise(eqn_plane_,cornersLocal);

    //integrating over each line of the facet
    double facet_integ = 0.0;
    if(bcellInt_==false)
    {
      std::vector<double> alpha;
      alpha = compute_alpha(eqn_plane_,GEO::CUT::proj_x);
      for(std::vector<std::vector<double> >::const_iterator k=cornersLocal.begin();k!=cornersLocal.end();k++)
      {
        const std::vector<double> coords1 = *k;
        std::vector<double> coords2;
        //for the last line the end point is the first point of the facet
        if(k!=(cornersLocal.end()-1))
          coords2 = *(k+1);
        else
          coords2= *(cornersLocal.begin());

  //first index decides the x or y coordinate, second index decides the start point or end point
        LINALG::Matrix<2,2> coordLine;

  //The facet is projected over y-z plane and then the integration is performed
  //so only y- and z-coordinates are passed to make the lines
  //[0]-- indicates y and [1] indicate z because we are now in y-z plane
        coordLine(0,0) = coords1[1];
        coordLine(1,0) = coords1[2];
        coordLine(0,1) = coords2[1];
        coordLine(1,1) = coords2[2];

        LineIntegration line1(coordLine,inte_num_,alpha,bcellInt_);
        facet_integ += line1.integrate_line();
      }
    }
    else
    {
//to reduce the truncation error introduced during the projection of plane,
//the plane, along which the normal component is maximum, is chosen
      GEO::CUT::ProjectionDirection projType;
      if(fabs(eqn_plane_[0])<1e-8)
      {
        if(fabs(eqn_plane_[1])<1e-8)
          projType = GEO::CUT::proj_z;
        else if(fabs(eqn_plane_[2])<1e-8)
          projType = GEO::CUT::proj_y;
        else
        {
          if(fabs(eqn_plane_[1])>fabs(eqn_plane_[2]))
            projType = GEO::CUT::proj_y;
          else
            projType = GEO::CUT::proj_z;
        }
      }
      else if(fabs(eqn_plane_[1])<1e-8)
      {
        if(fabs(eqn_plane_[2])<1e-8)
          projType = GEO::CUT::proj_x;
        else
        {
          if(fabs(eqn_plane_[0])>fabs(eqn_plane_[2]))
            projType = GEO::CUT::proj_x;
          else
            projType = GEO::CUT::proj_z;
        }
      }
      else if(fabs(eqn_plane_[2])<1e-8)
      {
        if(fabs(eqn_plane_[1])>fabs(eqn_plane_[0]))
          projType = GEO::CUT::proj_y;
        else
          projType = GEO::CUT::proj_x;
      }
      else
      {
        if(fabs(eqn_plane_[0])>=fabs(eqn_plane_[1]) && fabs(eqn_plane_[0])>=fabs(eqn_plane_[2]))
          projType = GEO::CUT::proj_x;
        else if(fabs(eqn_plane_[1])>=fabs(eqn_plane_[2]) && fabs(eqn_plane_[1])>=fabs(eqn_plane_[0]))
          projType = GEO::CUT::proj_y;
        else
          projType = GEO::CUT::proj_z;
      }

      if(projType!=GEO::CUT::proj_x and projType!=GEO::CUT::proj_y and projType!=GEO::CUT::proj_z)
      {
        dserror("projection plane is not defined");
        exit(1);
      }

      BoundaryFacetIntegration(cornersLocal,facet_integ,projType);

    }

//this condition results in negative normal for all the lines in the line integral
    if(clockwise_ && bcellInt_==false)
    {
      facet_integ = -1.0*facet_integ;
      for(int i=0;i<4;i++)
        eqn_plane_[i] = -1.0*eqn_plane_[i];
    }

    return facet_integ;
}

/*-----------------------------------------------------------------------------------------------*
                            Performs integration over the boundarycell
*------------------------------------------------------------------------------------------------*/
void GEO::CUT::FacetIntegration::BoundaryFacetIntegration( const std::vector<std::vector<double> > & cornersLocal,
                                                           double &facet_integ,
                                                           GEO::CUT::ProjectionDirection intType )
{
  std::vector<double> alpha;
  double abs_normal=0.0;

  for(std::vector<std::vector<double> >::const_iterator k=cornersLocal.begin();k!=cornersLocal.end();k++)
  {
    const std::vector<double> coords1 = *k;
    std::vector<double> coords2;

    //for the last line the end point is the first point of the facet
    if(k!=(cornersLocal.end()-1))
      coords2 = *(k+1);
    else
      coords2= *(cornersLocal.begin());

    //first index decides the x or y coordinate, second index decides the start point or end point
    LINALG::Matrix<2,2> coordLine;
    if(intType==GEO::CUT::proj_x)
    {
      if(k==cornersLocal.begin())
      {
        alpha = compute_alpha(eqn_plane_,GEO::CUT::proj_x);
        abs_normal = getNormal(intType);
      }
  //The facet is projected over y-z plane and then the integration is performed
  //so only y- and z-coordinates are passed to make the lines
  //(0,i)-- indicates y and (1,i) indicate z because we are now in y-z plane
      coordLine(0,0) = coords1[1];
      coordLine(1,0) = coords1[2];
      coordLine(0,1) = coords2[1];
      coordLine(1,1) = coords2[2];

      LineIntegration line1(coordLine,inte_num_,alpha,bcellInt_);
      line1.set_integ_type(GEO::CUT::proj_x);
      facet_integ += line1.integrate_line();
    }
    else if(intType==GEO::CUT::proj_y)
    {
      if(k==cornersLocal.begin())
      {
        alpha = compute_alpha(eqn_plane_,GEO::CUT::proj_y);
        abs_normal = getNormal(intType);
      }
  //The facet is projected over y-z plane and then the integration is performed
  //so only y- and z-coordinates are passed to make the lines
  //(0,i)-- indicates y and (1,i) indicate z because we are now in y-z plane
      coordLine(0,0) = coords1[2];
      coordLine(1,0) = coords1[0];
      coordLine(0,1) = coords2[2];
      coordLine(1,1) = coords2[0];
      LineIntegration line1(coordLine,inte_num_,alpha,bcellInt_);
      line1.set_integ_type(GEO::CUT::proj_y);
      facet_integ += line1.integrate_line();
    }

    else if(intType==GEO::CUT::proj_z)
    {
      if(k==cornersLocal.begin())
      {
        alpha = compute_alpha(eqn_plane_,GEO::CUT::proj_z);
        abs_normal = getNormal(intType);
      }
  //The facet is projected over y-z plane and then the integration is performed
  //so only y- and z-coordinates are passed to make the lines
  //(0,i)-- indicates y and (1,i) indicate z because we are now in y-z plane
      coordLine(0,0) = coords1[0];
      coordLine(1,0) = coords1[1];
      coordLine(0,1) = coords2[0];
      coordLine(1,1) = coords2[1];
      LineIntegration line1(coordLine,inte_num_,alpha,bcellInt_);
      line1.set_integ_type(GEO::CUT::proj_z);
      facet_integ += line1.integrate_line();
    }
    else
    {
      dserror("The facet integration type not supported");
      exit(1);
    }
  }

  facet_integ = facet_integ/abs_normal;
}

/*-------------------------------------------------------------------------------------------------*
      Generate integration rule for the facet if the divergence theorem is used      Sudhakar 03/12
      directly to generate Gauss integration rule for the facet
*--------------------------------------------------------------------------------------------------*/
void GEO::CUT::FacetIntegration::DivergenceIntegrationRule( Mesh &mesh,
                                                            Teuchos::RCP<DRT::UTILS::CollectedGaussPoints> & cgp )
{
  TEUCHOS_FUNC_TIME_MONITOR( "GEO::CUT::FacetIntegration::DivergenceIntegrationRule" );

  std::list<Teuchos::RCP<BoundaryCell> > divCells;

  //the last two parameters has no influence when called from the first parameter is set to true
  GenerateDivergenceCells( true, mesh, divCells );

  double normalX = getNormal(GEO::CUT::proj_x);//make sure eqn of plane is available before calling this

  if(clockwise_) //because if ordering is clockwise the contribution of this facet must be subtracted
    normalX = -1.0*normalX;

#ifdef DIRECTDIV_EXTENDED_DEBUG_OUTPUT
  static int mm = 0;
  std::stringstream str;
  str << "divCell" << mm << ".pos";
  std::ofstream file( str.str().c_str() );
  file << "View \"" << "FacetBCellInfo" << "\" {\n";
#endif
  for ( std::list<Teuchos::RCP<BoundaryCell> >::iterator i=divCells.begin();
                i!=divCells.end();
                ++i )
  {
    BoundaryCell * bcell = &**i;

#ifdef DIRECTDIV_EXTENDED_DEBUG_OUTPUT
    if(face1_->IsFacetSplit())
    {
      std::vector<std::vector<Point*> > facetSplit = face1_->GetSplitCells();
      for(std::vector<std::vector<Point*> >::iterator k = facetSplit.begin();
                                             k!=facetSplit.end();
                                             k++)
      {
        std::cout << "Split cell output." << std::endl;
        GEO::CUT::OUTPUT::GmshTriSideDump(file,*k);
      }
    }
    else
      GEO::CUT::OUTPUT::GmshTriSideDump(file,bcell->Points());

    LINALG::Matrix<3,1> midpt;
    bcell->ElementCenter(midpt);
    GEO::CUT::OUTPUT::GmshVector(file,midpt,GEO::CUT::OUTPUT::GetEqOfPlane(bcell->Points()), true);
#endif


#if 0//DEBUGCUTLIBRARY //write separate file for each bcell along with the distribution of Gauss points
    static int facetno = 0;
    facetno++;
    std::stringstream str;
    str << "facetid" << facetno << ".pos";
    std::ofstream file( str.str().c_str() );

    const std::vector<Point*> ptl = bcell->Points();
    int mm=0;
    for( unsigned ii=0;ii<ptl.size();ii++ )
    {
      Point* pt1 = ptl[ii];
      LINALG::Matrix<3,1> global1,local1;
      pt1->Coordinates(global1.A());
      elem1_->LocalCoordinates( global1, local1 );
      //file<<"Point("<<mm<<")="<<"{"<<local1(0,0)<<","<<local1(1,0)<<","<<local1(2,0)<<"};\n";
      file<<"Point("<<mm<<")="<<"{"<<global1(0,0)<<","<<global1(1,0)<<","<<global1(2,0)<<"};\n";
      mm++;
    }
    for( unsigned ii=0;ii<ptl.size();ii++ )
    {
      file<<"Line("<<ii<<")="<<"{"<<ii<<","<<(ii+1)%ptl.size()<<"};\n";
    }
#endif

    DRT::UTILS::GaussIntegration gi_temp = DRT::UTILS::GaussIntegration( bcell->Shape(), DIRECTDIV_GAUSSRULE );

    if( bcell->Area() < REF_AREA_BCELL )
      continue;

    for ( DRT::UTILS::GaussIntegration::iterator iquad=gi_temp.begin(); iquad!=gi_temp.end(); ++iquad )
    {
      double drs = 0.0;
      LINALG::Matrix<3,1> x_gp_loc(true), normal(true);
      const LINALG::Matrix<2,1> eta( iquad.Point() );

      switch ( bcell->Shape() )
      {
        case DRT::Element::tri3:
        {
#ifdef LOCAL
          bcell->TransformLocalCoords<DRT::Element::tri3>(elem1_,eta, x_gp_loc, normal, drs, true);
#else
          bcell->Transform<DRT::Element::tri3>( eta, x_gp_loc, normal, drs );
#endif
          break;
        }
        case DRT::Element::quad4:
        {
#ifdef LOCAL
          bcell->TransformLocalCoords<DRT::Element::quad4>(elem1_,eta, x_gp_loc, normal, drs, true);
#else
          bcell->Transform<DRT::Element::quad4>( eta, x_gp_loc, normal, drs );
#endif
          break;
        }
        default:
          throw std::runtime_error( "unsupported integration cell type" );
      }
      double wei = iquad.Weight()*drs*normalX;

      cgp->Append( x_gp_loc, wei );

#if 0//DEBUGCUTLIBRARY //write separate file for each bcell along with the distribution of Gauss points (contd...)
      file<<"Point("<<mm<<")="<<"{"<<x_gp_loc(0,0)<<","<<x_gp_loc(1,0)<<","<<x_gp_loc(2,0)<<"};\n";
      mm++;
#endif
    }
  }
#ifdef DIRECTDIV_EXTENDED_DEBUG_OUTPUT
  file << "};\n";
  mm++;
#endif
}

/*-----------------------------------------------------------------------------------------------------*
      Split the facet into auxillary divergence cells which is either Tri or Quad         Sudhakar 03/12
*------------------------------------------------------------------------------------------------------*/
void GEO::CUT::FacetIntegration::GenerateDivergenceCells( bool divergenceRule, //if called to generate direct divergence rule
                                                          Mesh &mesh,
                                                          std::list<Teuchos::RCP<BoundaryCell> >& divCells)
{
#ifdef LOCAL
  std::vector<std::vector<double> > cornersLocal;
  face1_->CornerPointsLocal(elem1_,cornersLocal);
#else
  std::vector<std::vector<double> > cornersLocal = face1_->CornerPointsGlobal(elem1_);
#endif

  eqn_plane_ = equation_plane(cornersLocal);

// the face is in the x-y or in y-z plane which will not be considered when divergence theorem is applied
  if(divergenceRule)
  {
    if(fabs(eqn_plane_[0])<TOL_EQN_PLANE)
    {
      return;
    }
  }

  if( !divergenceRule && !face1_->OnCutSide() )
    return;

  IsClockwise(eqn_plane_,cornersLocal);

  std::vector<Point*> corners = face1_->CornerPoints();

  if(divergenceRule)
  {
    //-----------------------------------------------------------------
    // only facets with 3 corners are considered as the special case and a Tri3 is created directly.
    // we cannot directly create a Quad4 by checking the number of corners, because it can be a
    // concave facet with 4 corners for which Gaussian rules are not available
    //-----------------------------------------------------------------
    if(corners.size()==3)
    {
      TemporaryTri3(corners, divCells);
    }
    else //split the aribtrary noded facet into cells
    {
      std::string splitMethod;

      std::vector<std::vector<Point*> > split;

      // if the facet is warped, do centre point triangulation --> reduced error (??)
      if( not face1_->IsPlanar( mesh, face1_->CornerPoints() ) )
      {
        if(!face1_->IsTriangulated())
          face1_->DoTriangulation( mesh, corners );
        split  = face1_->Triangulation();
      }
      // the facet is not warped
      else
      {
#if 1 // split facet
        if( !face1_->IsFacetSplit() )
          face1_->SplitFacet( corners );
        split = face1_->GetSplitCells();
        splitMethod = "split";
#endif

#if 0 // triangulate facet
        //std::cout<<"!!! WARNING !!! Facets are triangulated instead of getting splitted ---> more Gauss points\n";
        TriangulateFacet tf( corners );

        std::vector<int> ptconc;
        tf.EarClipping( ptconc, true );

        split = tf.GetSplitCells();
        splitMethod = "triangulation";
#endif
      }

      for ( std::vector<std::vector<Point*> >::const_iterator j=split.begin();
                                                              j!=split.end(); ++j )
      {
        const std::vector<Point*> & tri = *j;
        if(tri.size()==3)
          TemporaryTri3(tri, divCells);
        else if(tri.size()==4) // split algorithm always gives convex quad
          TemporaryQuad4(tri, divCells);
        else
        {
          std::cout<<"number of sides = "<<tri.size();
          dserror("Splitting created neither tri3 or quad4");
        }
      }

#ifdef DEBUGCUTLIBRARY // check the area of facet computed from splitting and triangulation
  if(splitMethod=="split")
  {
    DebugAreaCheck( divCells, splitMethod, mesh );
  }
#endif
    }
  }

}

/*-------------------------------------------------------------------------------------------*
                            temporarily create a tri3 cell
                    this is temporary because this is not stored for the volumecell
*--------------------------------------------------------------------------------------------*/
void GEO::CUT::FacetIntegration::TemporaryTri3( const std::vector<Point*>& corners,
                                                std::list<Teuchos::RCP<BoundaryCell> >& divCells )
{
  Epetra_SerialDenseMatrix xyz( 3, 3 );
  for ( int i=0; i<3; ++i )
    corners[i]->Coordinates( &xyz( 0, i ) );
  Tri3BoundaryCell * bc = new Tri3BoundaryCell( xyz, face1_, corners );
  divCells.push_back( Teuchos::rcp( bc ) );
}

/*-------------------------------------------------------------------------------------------*
                            temporarily create a quad4 cell
                    this is temporary because this is not stored for the volumecell
*--------------------------------------------------------------------------------------------*/
void GEO::CUT::FacetIntegration::TemporaryQuad4( const std::vector<Point*>& corners,
                                                 std::list<Teuchos::RCP<BoundaryCell> >& divCells )
{
  Epetra_SerialDenseMatrix xyz( 3, 4 );
  for ( int i=0; i<4; ++i )
    corners[i]->Coordinates( &xyz( 0, i ) );
  Quad4BoundaryCell * bc = new Quad4BoundaryCell( xyz, face1_, corners );
  divCells.push_back( Teuchos::rcp( bc ) );
}

/*----------------------------------------------------------------------------------------------*
          Check whether the area of the facet computed by the splitting and
          triangular procedure are the same.                                       sudhakar 05/12
          This is a check for the adopted splitting procedures
*-----------------------------------------------------------------------------------------------*/
void GEO::CUT::FacetIntegration::DebugAreaCheck( std::list<Teuchos::RCP<BoundaryCell> > & divCells,
                                                 std::string alreadyDone,
                                                 Mesh &mesh )
{
  std::cout << RED << "FacetIntegration::DebugAreaCheck --- This Test is broken for facets with wholes inside!!!" << END_COLOR << std::endl;
  //finally we should just have one triangulation/splitting proceducre which handles automatically holes, not very useful if we have to treat
  //this on every step where we need triangulation ... Todo ...error is commented out because of that!!!

  double area1 = 0.0;
  for ( std::list<Teuchos::RCP<BoundaryCell> >::iterator i=divCells.begin();
                i!=divCells.end();
                ++i )
  {
    BoundaryCell * bcell = &**i;
    area1 += bcell->Area();
  }

  std::vector<Point*> corners = face1_->CornerPoints();

  std::vector<std::vector<GEO::CUT::Point*> > split1;
  if( alreadyDone=="split" )
  {
    TriangulateFacet tf( corners );
    std::vector<int> ptconc;

    tf.EarClipping( ptconc, true, true );
    split1 = tf.GetSplitCells();
  }
  else
  {
    TriangulateFacet tf( corners );
    tf.SplitFacet();
    split1 = tf.GetSplitCells();
  }

  std::list<Teuchos::RCP<BoundaryCell> > tempCells;

  for ( std::vector<std::vector<Point*> >::const_iterator j=split1.begin();
                                                          j!=split1.end(); ++j )
  {
    std::vector<Point*> tri = *j;
    if(tri.size()==3)
      TemporaryTri3(tri, tempCells);
    else if(tri.size()==4)
      TemporaryQuad4(tri, tempCells);
    else
    {
      std::cout<<"number of sides = "<<tri.size();
      dserror("Splitting created neither tri3 or quad4");
    }
  }

  double area2 = 0.0;
  for ( std::list<Teuchos::RCP<BoundaryCell> >::iterator i=tempCells.begin();
          i!=tempCells.end();
          ++i )
  {
    BoundaryCell * bcell = &**i;
    area2 += bcell->Area();
  }

  //std::cout.precision(15);
  if( fabs(area1-area2)>1e-10 )
  {
    std::cout<<"The coordinates of the facet\n";
    for( unsigned i=0;i<corners.size();i++ )
    {
      Point* pt = corners[i];
      double coo[3];
      pt->Coordinates(coo);
      std::cout<<std::setprecision(20)<<coo[0]<<"\t"<<coo[1]<<"\t"<<coo[2]<<"\n";
    }

    std::cout<<"cells produced by Method 1\n";
    for ( std::list<Teuchos::RCP<BoundaryCell> >::iterator i=divCells.begin();
              i!=divCells.end();
              ++i )
    {
      BoundaryCell * bcell = &**i;
      const std::vector<Point*> pts = bcell->Points();
      std::cout<<"cells\n";
      for( unsigned j=0;j<pts.size();j++ )
      {
        Point *pt1 = pts[j];
        double coo[3];
        pt1->Coordinates(coo);
        std::cout<<coo[0]<<"\t"<<coo[1]<<"\t"<<coo[2]<<"\n";
      }
    }

    std::cout<<"cells produced by Method 2\n";
    for ( std::list<Teuchos::RCP<BoundaryCell> >::iterator i=tempCells.begin();
              i!=tempCells.end();
              ++i )
    {
      BoundaryCell * bcell = &**i;
      const std::vector<Point*> pts = bcell->Points();
      std::cout<<"cells\n";
      for( unsigned j=0;j<pts.size();j++ )
      {
        Point *pt1 = pts[j];
        double coo[3];
        pt1->Coordinates(coo);
        std::cout<<coo[0]<<"\t"<<coo[1]<<"\t"<<coo[2]<<"\n";
      }
    }

    std::cout<<"Area1 = "<<area1<<"\t"<<"Area2 = "<<area2<<"\n";
    std::cout<<RED<<"!!!WARNING!!! area predicted by splitting and triangulation are not the same\n"<<END_COLOR;
    //throw std::runtime_error( "error: area predicted by splitting and triangulation are not the same" );
    //dserror( "error: area predicted by splitting and triangulation are not the same" );


    /***************************************************************************************/
    // change splitting cells into triangulated cells
    /*divCells.clear();

    for ( std::vector<std::vector<Point*> >::const_iterator j=split1.begin();
                                                            j!=split1.end(); ++j )
    {
      std::vector<Point*> tri = *j;
      if(tri.size()==3)
        TemporaryTri3(tri, divCells);
      else if(tri.size()==4) // split algorithm always gives convex quad
        TemporaryQuad4(tri, divCells);
      else
      {
        std::cout<<"number of sides = "<<tri.size();
        dserror("Splitting created neither tri3 or quad4");
      }
    }*/
    /***************************************************************************************/
  }
}


/*-------------------------------------------------------------------------------------------------*
      Generate integration rule for the facet if the divergence theorem is used      Sudhakar 03/12
      directly to generate Gauss integration rule for the facet
*--------------------------------------------------------------------------------------------------*/
void GEO::CUT::FacetIntegration::DivergenceIntegrationRuleNew( Mesh &mesh,
                                                            Teuchos::RCP<DRT::UTILS::CollectedGaussPoints> & cgp )
{
  TEUCHOS_FUNC_TIME_MONITOR( "GEO::CUT::FacetIntegration::DivergenceIntegrationRule" );

  std::list<Teuchos::RCP<BoundaryCell> > divCells;

  //the last two parameters has no influence when called from the first parameter is set to true
  std::vector<std::vector<double> >        eqn_plane_divCell;

  //If the facet is not planar it will be triangulated in DirectDivergence::ListFacets().
  // Might want to split the facet for the case it is a planar quad -> less divCells.

#ifndef TRIANGULATE_ALL_FACETS_FOR_DIVERGENCECELLS
   bool triangulate_and_levelset = ( face1_->CornerPoints().size()>3 and face1_->BelongsToLevelSetSide() );
#else
   bool triangulate_and_levelset = ( face1_->CornerPoints().size()>3 );
#endif

//  if(face1_->CornerPoints().size()>3 and face1_->BelongsToLevelSetSide())
  if(triangulate_and_levelset)
    face1_->DoTriangulation(mesh,face1_->CornerPoints());

#ifndef TRIANGULATE_ALL_FACETS_FOR_DIVERGENCECELLS
  triangulate_and_levelset = ( (face1_->IsTriangulated() or face1_->IsFacetSplit()) and face1_->BelongsToLevelSetSide());
#else
  triangulate_and_levelset = (face1_->IsTriangulated() or face1_->IsFacetSplit());
#endif

//  if((face1_->IsTriangulated() or face1_->IsFacetSplit()) and face1_->BelongsToLevelSetSide())
  if(triangulate_and_levelset)
  {
    std::vector<std::vector<Point*> > facet_triang;
    if(face1_->IsTriangulated())
      facet_triang = face1_->Triangulation();
    else
      facet_triang  = face1_->GetSplitCells();
#ifdef DIRECTDIV_EXTENDED_DEBUG_OUTPUT
    std::cout << "TRIANGULATED OR SPLIT FACET!" << std::endl;
#endif

    unsigned counterDivCells =0;
    for ( std::vector<std::vector<Point*> >::const_iterator j=facet_triang.begin();
        j!=facet_triang.end(); ++j )
    {
      GenerateDivergenceCellsNew(true,mesh,divCells,*j);
      while(counterDivCells < divCells.size())
      {
        counterDivCells++;
        eqn_plane_divCell.push_back(eqn_plane_);
      }
    }
  }
  else
  {
    //OLD IMPLEMENTATION!
    //Maybe triangulate for cut side?
#ifdef DIRECTDIV_EXTENDED_DEBUG_OUTPUT
    std::cout << "NOT TRIANGULATED FACET." << std::endl;
#endif

    GenerateDivergenceCellsNew(true,mesh,divCells,face1_->CornerPoints());
    //GenerateDivergenceCells(true,mesh,divCells);
    for(unsigned i=0; i<divCells.size(); i++)
    {
      eqn_plane_divCell.push_back(eqn_plane_);
    }
  }

  //SAFETY-CHECK
#ifdef DEBUGCUTLIBRARY
  if(eqn_plane_divCell.size() != divCells.size())
  {
    std::cout << "FINAL!: " << std::endl;
    std::vector<std::vector<double> > cornersLocal;
    face1_->CornerPointsLocal(elem1_,cornersLocal);
    for( std::vector<std::vector<double> >::const_iterator j=cornersLocal.begin();
                                           j!=cornersLocal.end(); ++j)
    {
      std::cout << "cornerLocalFacet: " << (*j)[0] << ", " << (*j)[1] << ", " << (*j)[2] << std::endl;
    }
    std::cout << "eqn_plane_divCell.size(): " << eqn_plane_divCell.size() << std::endl;
    std::cout << "divCells.size():  " << divCells.size() << std::endl;
    //dserror("Something wrong with divCell and clockwise assignment");
    throw std::runtime_error("Something wrong with divCell and clockwise assignment.");
  }
#endif

  double normalX;

#ifdef DIRECTDIV_EXTENDED_DEBUG_OUTPUT
  static int mm = 0;
  std::stringstream str;
  str << "divCell" << mm << ".pos";
  std::ofstream file( str.str().c_str() );
  file << "View \"" << "FacetBCellInfo" << "\" {\n";
#endif
  int zz=0;
  for ( std::list<Teuchos::RCP<BoundaryCell> >::iterator i=divCells.begin();
                i!=divCells.end();
                ++i )
  {
    BoundaryCell * bcell = &**i;

    //Get equation of plane for divergence Cell.
    std::vector<double> eqn_plane_bcell = eqn_plane_divCell[zz];

    double normalScale=0.0;
    for(unsigned i=0;i<3;i++)
      normalScale += eqn_plane_bcell[i]*eqn_plane_bcell[i];
    normalScale = sqrt(normalScale);
    normalX = eqn_plane_bcell[0]/normalScale;

#ifdef DIRECTDIV_EXTENDED_DEBUG_OUTPUT
    if(face1_->IsFacetSplit())
    {
      std::vector<std::vector<Point*> > facetSplit = face1_->GetSplitCells();
      for(std::vector<std::vector<Point*> >::iterator k = facetSplit.begin();
                                             k!=facetSplit.end();
                                             k++)
      {
        LINALG::Matrix<3,1> midpt(true);

        std::vector<std::vector<double> > corners_split;
        for(std::vector<Point*>::iterator l = (*k).begin();
                                                     l!=(*k).end();
                                                     l++)
        {
          const double* coords = (*l)->X();
          midpt(0,0) += coords[0];
          midpt(1,0) += coords[1];
          midpt(2,0) += coords[2];
        }

        midpt.Scale(1.0/(*k).size());

//        bcell->ElementCenter(midpt);
        GEO::CUT::OUTPUT::GmshTriSideDump(file,*k);
        GEO::CUT::OUTPUT::GmshVector(file,midpt,GEO::CUT::OUTPUT::GetEqOfPlane(*k), true);
      }
    }
    else
    {
      LINALG::Matrix<3,1> midpt(true);
      bcell->ElementCenter(midpt);
      GEO::CUT::OUTPUT::GmshTriSideDump(file,bcell->Points());
      GEO::CUT::OUTPUT::GmshVector(file,midpt,GEO::CUT::OUTPUT::GetEqOfPlane(bcell->Points()), true);
    }
#endif


#if 0//DEBUGCUTLIBRARY //write separate file for each bcell along with the distribution of Gauss points
    static int facetno = 0;
    facetno++;
    std::stringstream str;
    str << "facetid" << facetno << ".pos";
    std::ofstream file( str.str().c_str() );

    const std::vector<Point*> ptl = bcell->Points();
    int mm=0;
    for( unsigned ii=0;ii<ptl.size();ii++ )
    {
      Point* pt1 = ptl[ii];
      LINALG::Matrix<3,1> global1,local1;
      pt1->Coordinates(global1.A());
      elem1_->LocalCoordinates( global1, local1 );
      //file<<"Point("<<mm<<")="<<"{"<<local1(0,0)<<","<<local1(1,0)<<","<<local1(2,0)<<"};\n";
      file<<"Point("<<mm<<")="<<"{"<<global1(0,0)<<","<<global1(1,0)<<","<<global1(2,0)<<"};\n";
      mm++;
    }
    for( unsigned ii=0;ii<ptl.size();ii++ )
    {
      file<<"Line("<<ii<<")="<<"{"<<ii<<","<<(ii+1)%ptl.size()<<"};\n";
    }
#endif

    DRT::UTILS::GaussIntegration gi_temp = DRT::UTILS::GaussIntegration( bcell->Shape(), DIRECTDIV_GAUSSRULE );

    if( bcell->Area() < REF_AREA_BCELL )
      continue;

    for ( DRT::UTILS::GaussIntegration::iterator iquad=gi_temp.begin(); iquad!=gi_temp.end(); ++iquad )
    {
      double drs = 0.0;
      LINALG::Matrix<3,1> x_gp_loc(true), normal(true);
      const LINALG::Matrix<2,1> eta( iquad.Point() );

      switch ( bcell->Shape() )
      {
        case DRT::Element::tri3:
        {
#ifdef LOCAL
          bcell->TransformLocalCoords<DRT::Element::tri3>(elem1_,eta, x_gp_loc, normal, drs, true);
#else
          bcell->Transform<DRT::Element::tri3>( eta, x_gp_loc, normal, drs );
#endif
          break;
        }
        case DRT::Element::quad4:
        {
#ifdef LOCAL
          bcell->TransformLocalCoords<DRT::Element::quad4>(elem1_,eta, x_gp_loc, normal, drs, true);
#else
          bcell->Transform<DRT::Element::quad4>( eta, x_gp_loc, normal, drs );
#endif
          break;
        }
        default:
          dserror( "unsupported integration cell type ( cell type = %s )",
              DRT::DistypeToString( bcell->Shape() ).c_str() );
          exit( EXIT_FAILURE );
      }
      double wei = iquad.Weight()*drs*normalX;

      cgp->Append( x_gp_loc, wei );

#if 0//DEBUGCUTLIBRARY //write separate file for each bcell along with the distribution of Gauss points (contd...)
      file<<"Point("<<mm<<")="<<"{"<<x_gp_loc(0,0)<<","<<x_gp_loc(1,0)<<","<<x_gp_loc(2,0)<<"};\n";
      mm++;
#endif
    }
    zz++; //Iterator of divCells.
  }
#ifdef DIRECTDIV_EXTENDED_DEBUG_OUTPUT
  file << "};\n";
  mm++;
#endif
}

void GEO::CUT::FacetIntegration::GenerateDivergenceCellsNew(bool       divergenceRule,
                             Mesh                                      &mesh,
                             std::list<Teuchos::RCP<BoundaryCell> >    & divCells,
                             const std::vector<Point*>&                cornersGlobal)
{
//First convert format...
#ifdef LOCAL
  std::vector<std::vector<double> > cornersLocal;//(cornersGlobal.size(),3);
  std::vector<std::vector<double> > cornersGlobalVec;//(cornersGlobal.size(),3);

  for(std::vector<Point*> ::const_iterator j=cornersGlobal.begin();
      j!=cornersGlobal.end(); ++j )
  {

    LINALG::Matrix<3,1> cornersLocMatrix;
    LINALG::Matrix<3,1> cornersGloMatrix((*j)->X());

    elem1_->LocalCoordinates(cornersGloMatrix,cornersLocMatrix);

    std::vector<double> cornerLocal(3);
    cornerLocal[0] = cornersLocMatrix(0,0);
    cornerLocal[1] = cornersLocMatrix(1,0);
    cornerLocal[2] = cornersLocMatrix(2,0);
    cornersLocal.push_back(cornerLocal);

//    std::cout << "cornerLocal: " << cornerLocal[0] << ", " << cornerLocal[1] << ", " << cornerLocal[2] << std::endl;

    std::vector<double> cornerGlobal(3);
    cornerGlobal[0] = cornersGloMatrix(0,0);
    cornerGlobal[1] = cornersGloMatrix(1,0);
    cornerGlobal[2] = cornersGloMatrix(2,0);
    cornersGlobalVec.push_back(cornerGlobal);
  }
#else
  std::vector<std::vector<double> > cornersGlobalVec;//(cornersGlobal.size(),3);

  for(std::vector<Point*> ::const_iterator j=cornersGlobal.begin();
      j!=cornersGlobal.end(); ++j )
  {
    LINALG::Matrix<3,1> cornersGloMatrix((*j)->X());

    std::vector<double> cornerGlobal(3);
    cornerGlobal[0] = cornersGloMatrix(0,0);
    cornerGlobal[1] = cornersGloMatrix(1,0);
    cornerGlobal[2] = cornersGloMatrix(2,0);
    cornersGlobalVec.push_back(cornerGlobal);
  }
  std::vector<std::vector<double> > cornersLocal = cornersGlobalVec; //face1_->CornerPointsGlobal(elem1_);
#endif

  //Equation of plane for triangulated/split facet or full facet
  eqn_plane_ = equation_plane(cornersLocal);

// the face is in the x-y or in y-z plane which will not be considered when divergence theorem is applied
  if(divergenceRule)
  {
    if(fabs(eqn_plane_[0])<TOL_EQN_PLANE/10.0)
    {
#ifdef DIRECTDIV_EXTENDED_DEBUG_OUTPUT
      for(std::vector<std::vector<double> > ::const_iterator j=cornersLocal.begin();
          j!=cornersLocal.end(); ++j )
      {
        std::cout << "cornerLocal: " << (*j)[0] << ", " << (*j)[1] << ", " << (*j)[2] << std::endl;
      }
      std::cout << "eqn_plane_: " << eqn_plane_[0] << ", " << eqn_plane_[1] << ", " << eqn_plane_[2] << ", " << eqn_plane_[3] << std::endl;
      std::cout << "RETURN BECAUSE eqn_plane_ too small in x-dir." << std::endl;
#endif
      return;
    }
  }

  if( !divergenceRule && !face1_->OnCutSide() )
    return;

  IsClockwise(eqn_plane_,cornersLocal);

  std::vector<Point*> corners = cornersGlobal;

  if(clockwise_)
  {
    //CHANGE SIGN For equation of plane.
    eqn_plane_[0] = -eqn_plane_[0];
    eqn_plane_[1] = -eqn_plane_[1];
    eqn_plane_[2] = -eqn_plane_[2];
    eqn_plane_[3] = -eqn_plane_[3];
  }

  if(divergenceRule)
  {
    //-----------------------------------------------------------------
    // only facets with 3 corners are considered as the special case and a Tri3 is created directly.
    // we cannot directly create a Quad4 by checking the number of corners, because it can be a
    // concave facet with 4 corners for which Gaussian rules are not available
    //-----------------------------------------------------------------
    if(corners.size()==3)
    {
      TemporaryTri3(corners, divCells);
    }
    else //split the aribtrary noded facet into cells
    {
      std::string splitMethod;

      std::vector<std::vector<Point*> > split;

      // if the facet is warped, do centre point triangulation --> reduced error (??)
      if( not face1_->IsPlanar( mesh, face1_->CornerPoints() ) )
      {
        if(!face1_->IsTriangulated())
          face1_->DoTriangulation( mesh, corners );
        split  = face1_->Triangulation();
      }
      // the facet is not warped
      else
      {
#if 1 // split facet
        if( !face1_->IsFacetSplit() )
          face1_->SplitFacet( corners );
        split = face1_->GetSplitCells();
        splitMethod = "split";
#endif

#if 0 // triangulate facet
        //std::cout<<"!!! WARNING !!! Facets are triangulated instead of getting splitted ---> more Gauss points\n";
        TriangulateFacet tf( corners );

        std::vector<int> ptconc;
        tf.EarClipping( ptconc, true );

        split = tf.GetSplitCells();
        splitMethod = "triangulation";
#endif
      }

      for ( std::vector<std::vector<Point*> >::const_iterator j=split.begin();
                                                              j!=split.end(); ++j )
      {
        const std::vector<Point*> & tri = *j;
        if(tri.size()==3)
          TemporaryTri3(tri, divCells);
        else if(tri.size()==4) // split algorithm always gives convex quad
          TemporaryQuad4(tri, divCells);
        else
        {
          std::cout<<"number of sides = "<<tri.size();
          dserror("Splitting created neither tri3 or quad4");
        }
      }

#ifdef DEBUGCUTLIBRARY // check the area of facet computed from splitting and triangulation
  if(splitMethod=="split")
  {
    DebugAreaCheck( divCells, splitMethod, mesh );
  }
#endif
    }
  }
}
