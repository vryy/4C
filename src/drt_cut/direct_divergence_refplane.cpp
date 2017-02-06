/*!-----------------------------------------------------------------------------------------------*
\file direct_divergence_refplane.cpp

\brief Construct reference plane for direct divergence method when used in global
coordinate system

\level 2

<pre>
\maintainer Magnus Winter
            winter@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>
 *------------------------------------------------------------------------------------------------*/
#include "direct_divergence_refplane.H"
#include "cut_kernel.H"
#include "cut_side.H"

/*-----------------------------------------------------------------------------------*
 * Perform all the operations related to computing reference plane           sudhakar 06/15
 *-----------------------------------------------------------------------------------*/
std::vector<double> GEO::CUT::DirectDivergenceGlobalRefplane::GetReferencePlane()
{



  if( elem1_->Shape() != DRT::Element::hex8 )
  {
    std::cout << "Element Shape: " << elem1_->Shape() << std::endl;
    throw std::runtime_error("Currently can handle only hexagonal family\n");
    //dserror("Currently can handle only hexagonal family\n");
  }

  std::vector<double> RefPlaneEqn( 4, 0.0 );

  bool comp_ref_plane = false;

  //---
  // First estimate -- Compute reference plane based on the diagonal information
  //---
  comp_ref_plane = DiagonalBasedRef( RefPlaneEqn );
  if( comp_ref_plane )
  {
    return RefPlaneEqn;
  }
  refPtsGmsh_.clear();

  //---
  // Second estimate -- Compute reference plane based on Side information
  //---
  comp_ref_plane = SideBasedRef( RefPlaneEqn );
  if( comp_ref_plane )
  {
    return RefPlaneEqn;
  }

  if( not comp_ref_plane )
    dserror( "Proper reference plane not found" );

  return RefPlaneEqn;
}

/*-------------------------------------------------------------------------------------------------------*
 * Computation of reference plane based on the diagonal of background element                    sudhakar 06/15
 * This considers all 6 diagonals of background Hex element, and choose the one
 * that has maximum normal component in x-direction
 *-------------------------------------------------------------------------------------------------------*/
bool GEO::CUT::DirectDivergenceGlobalRefplane::DiagonalBasedRef( std::vector<double>& RefPlaneEqn )
{
  // Any method should involve the following two steps

  //---
  // STEP 1: Estimate equation of reference plane and store the corresponding points for gmsh output.
  // Here we take all possible 6 diagonals of the Hex element, and choose the diagonal which has
  // maximum normal component in x-direction as the reference plane
  //---
  std::vector<Point*> ptslist = elem1_->Points();
  std::vector<std::vector<Point*> >diagonals;

  std::vector<Point*> diag;
  diag.push_back( ptslist[0] );
  diag.push_back( ptslist[1] );
  diag.push_back( ptslist[6] );
  diag.push_back( ptslist[7] );
  diagonals.push_back( diag );

  diag.clear();
  diag.push_back( ptslist[2] );
  diag.push_back( ptslist[3] );
  diag.push_back( ptslist[4] );
  diag.push_back( ptslist[5] );
  diagonals.push_back( diag );

  diag.clear();
  diag.push_back( ptslist[5] );
  diag.push_back( ptslist[6] );
  diag.push_back( ptslist[3] );
  diag.push_back( ptslist[0] );
  diagonals.push_back( diag );

  diag.clear();
  diag.push_back( ptslist[4] );
  diag.push_back( ptslist[7] );
  diag.push_back( ptslist[2] );
  diag.push_back( ptslist[1] );
  diagonals.push_back( diag );

  diag.clear();
  diag.push_back( ptslist[0] );
  diag.push_back( ptslist[4] );
  diag.push_back( ptslist[6] );
  diag.push_back( ptslist[2] );
  diagonals.push_back( diag );

  diag.clear();
  diag.push_back( ptslist[5] );
  diag.push_back( ptslist[1] );
  diag.push_back( ptslist[3] );
  diag.push_back( ptslist[7] );
  diagonals.push_back( diag );

  double xnormal = 0.0;
  bool found_refplane = false;
  for( std::vector<std::vector<Point*> >::iterator itd = diagonals.begin();
                                                   itd != diagonals.end(); itd++ )
  {
    std::vector<Point*> ptl = *itd;
    std::vector<double> RefPlaneTemp = KERNEL::EqnPlaneOfPolygon( ptl );

    scaleEquationOfPlane( RefPlaneTemp );

    if( fabs(RefPlaneTemp[0]) < REF_PLANE_DIRDIV )
      continue;

    //---
    // STEP 2: Project all the corner points of the Hex element onto the reference plane.
    // If all these projected points are within the background element, then we consider this as a possible reference plane!
    if  (isAllProjectedCornersInsideEle( RefPlaneTemp ))
    {
      double fac = sqrt( pow(RefPlaneTemp[0],2)+pow(RefPlaneTemp[1],2)+pow(RefPlaneTemp[2],2) );
      //---
      //STEP 3: Take the reference plane with biggest component in x-direction of the normal vector!
      double xn = fabs(RefPlaneTemp[0]) / fac;
      if( xn > xnormal )
      {
        xnormal = xn;
        RefPlaneEqn = RefPlaneTemp;
        refPtsGmsh_ = ptl;
        found_refplane = true;
      }
    }
  }

  //---
  //Basically using a diagonal reference plane should be enought, otherwise we have to look into that again!
  if (!found_refplane) dserror("Couldn't find a diagonal reference plane with all gausspoints inside!!!");
  return found_refplane;
}

/*-------------------------------------------------------------------------------------------------------*
 * Computation of reference plane based on the sides of background element                    sudhakar 06/15
 * Sort all the sides based on n_x (the one has more n_x gets on the top)
 * Iterate through all the sides to get the correct reference plane
 *-------------------------------------------------------------------------------------------------------*/
bool GEO::CUT::DirectDivergenceGlobalRefplane::SideBasedRef( std::vector<double>& RefPlaneEqn )
{
  const std::vector<Side*> & allsides = elem1_->Sides();

  //---
  // STEP 1: Estimate equation of reference plane and store the corresponding points for gmsh output.
  // First get all the sides of the element and compute the equation of plane for each side
  // Store them in a data structure which stores the sides based on n_x in descending order
  //---
  std::multimap<double,std::pair<std::vector<double>,std::vector<Point*> >,compareClass> side_data;
  for( std::vector<Side*>::const_iterator it = allsides.begin(); it != allsides.end(); it++ )
  {
    const Side* s = *it;
    const std::vector<Node*> nds = s->Nodes();

    std::vector<Point*> ptside;
    for( std::vector<Node*>::const_iterator itn = nds.begin(); itn != nds.end(); itn++ )
      ptside.push_back( (*itn)->point() );

    std::vector<double> RefPlaneTemp = KERNEL::EqnPlaneOfPolygon( ptside );
    scaleEquationOfPlane( RefPlaneTemp );
    if( fabs(RefPlaneTemp[0]) < REF_PLANE_DIRDIV )
      continue;

    side_data.insert( std::make_pair( RefPlaneTemp[0], std::make_pair( RefPlaneTemp, ptside )  ) );
  }

  //---
  // STEP 2: Project all the corner points of the Hex element onto the reference plane.
  // If all these projected points are within the background element, then we have the right choice
  // of reference plane
  //---
  for( std::multimap<double,std::pair<std::vector<double>,std::vector<Point*> > >::iterator it = side_data.begin();
                                                                                           it != side_data.end(); it++ )
  {
    RefPlaneEqn = it->second.first;
    bool isWithin = isAllProjectedCornersInsideEle( RefPlaneEqn );
    if( isWithin )
    {
      refPtsGmsh_ = it->second.second;
      return true;
    }
  }

  return false;
}

/*---------------------------------------------------------------------------------------------------------------*
 * In order to check whether the chosen reference plane is a correct choice,
 * we project all the corner points of the element onto this reference plane                              sudhakar 06/15
 * If all these projected points are within the element, then we got the correct ref plane
 *---------------------------------------------------------------------------------------------------------------*/
bool GEO::CUT::DirectDivergenceGlobalRefplane::isAllProjectedCornersInsideEle( std::vector<double>& RefPlaneEqn )
{
  std::vector<Point*> corners = elem1_->Points();

  for( std::vector<Point*>::iterator it = corners.begin(); it != corners.end(); it++ )
  {
    Point * pt = *it;

    LINALG::Matrix<3,1> coo;
    pt->Coordinates( &coo(0,0) );

    LINALG::Matrix<3,1> xyz_proj( coo ), rst_proj;

    // x-coordinate of corner pt projected in the reference plane
    // y- and z-coordinate remain the same
    xyz_proj(0,0) = (RefPlaneEqn[3]-RefPlaneEqn[1]*coo(1,0)-RefPlaneEqn[2]*coo(2,0))/RefPlaneEqn[0];

    // get the local coordinates of the projected point
    elem1_->LocalCoordinates( xyz_proj, rst_proj );

    // Check whether the local coordinate of the projected point is within the specified limits
    if( std::abs(rst_proj(0,0)) > 1.0+1e-8 or std::abs(rst_proj(1,0)) > 1.0+1e-8 or std::abs(rst_proj(2,0)) > 1.0+1e-8 or
       std::isnan(rst_proj(0,0)) or std::isnan(rst_proj(1,0)) or std::isnan(rst_proj(2,0)))
      return false;
  }
  return true;
}

/*----------------------------------------------------------------------------------------------------------*
 * Scale the equation of plane to enable comparison of normals                                  sudhakar 07/15
 *----------------------------------------------------------------------------------------------------------*/
void GEO::CUT::DirectDivergenceGlobalRefplane::scaleEquationOfPlane( std::vector<double>& RefPlaneEqn )
{
  double scale = sqrt(pow( RefPlaneEqn[0], 2.0 ) + pow( RefPlaneEqn[1], 2.0 ) + pow( RefPlaneEqn[2], 2.0 ));
  for( unsigned i=0; i<4; i++ )
    RefPlaneEqn[i] /= scale;
}
