/*!-----------------------------------------------------------------------------------------------*
\file direct_divergence.cpp

\brief Generate main Gauss points when using "DirectDivergence" approach.
equations
 *------------------------------------------------------------------------------------------------*/

#include "direct_divergence.H"
#include "facet_integration.H"
#include "volume_integration.H"

#include "cut_boundingbox.H"

/*-------------------------------------------------------------------------------------------------------------------*
  Create integration points on the facets of the volumecell by triangulating the facets
  A reference facet is identified on which integration weights are set to zero                        Sudhakar 04/12
*--------------------------------------------------------------------------------------------------------------------*/
Teuchos::RCP<DRT::UTILS::GaussPoints> GEO::CUT::DirectDivergence::VCIntegrationRule( std::vector<double> &RefPlaneEqn )
{
  std::vector<plain_facet_set::const_iterator> facetIterator; //iterators of facets which need to be considered for integration rule
  plain_facet_set::const_iterator IteratorRefFacet;           //iterator for the reference facet
  isRef_ = false;                                             //whether ref plane is falling on facet?

  // get integration facets and reference plane
  ListFacets( facetIterator, RefPlaneEqn, IteratorRefFacet, isRef_ );

  if( isRef_ )
    refFacet_ = *IteratorRefFacet;

  if( facetIterator.size()==0 )
    dserror( "x-component normal is zero on all the facets? It should not be." );

  Teuchos::RCP<DRT::UTILS::CollectedGaussPoints> cgp = Teuchos::rcp( new DRT::UTILS::CollectedGaussPoints(0) );

  for( unsigned i=0;i<facetIterator.size();i++ )
  {
    plain_facet_set::const_iterator iter = facetIterator[i];
    Facet * fe = *iter;
    FacetIntegration faee1(fe,elem1_,position_,false,false);

    faee1.DivergenceIntegrationRule( mesh_, cgp );
  }

  DRT::UTILS::GaussIntegration gi(cgp);

#if 0 //integrate specified functions using the Gaussian rule generated -- used in postprocessing
  IntegrateSpecificFuntions( gi, RefPlaneEqn );  //integrate specific functions
#endif

  return cgp;
}

/*-------------------------------------------------------------------------------------------------------*
        Identify the list of facets which need to be triangulated (called as integration facets)
        Get the reference facet that will be used to get internal Gauss rule                  (sudhakar 04/12)
        As far as possible, the reference facet is set on a cut side to reduce the no of Gauss pts
        Reference facet is selected so that all internal points are inside the volumecell
*--------------------------------------------------------------------------------------------------------*/
void GEO::CUT::DirectDivergence::ListFacets( std::vector<plain_facet_set::const_iterator>& facetIterator,
                                             std::vector<double>& RefPlaneEqn,
                                             plain_facet_set::const_iterator& IteratorRefFacet,
                                             bool & IsRefFacet )
{
  const plain_facet_set & facete = volcell_->Facets();

  bool RefOnCutSide=false;
  std::vector<std::vector<double> > eqnAllFacets(facete.size());

  // store the iterators of  all warped facets
  std::vector<plain_facet_set::const_iterator> warpFac;

  // check whether all facets of this vc are oriented in a plane
  // if not then some sides are warped
  // we need to generate quadrature rules in global coordinates
  bool inGlobal = false;
  for(plain_facet_set::const_iterator i=facete.begin();i!=facete.end();i++)
  {
    Facet *fe = *i;
    std::vector<Point*> corn = fe->CornerPoints();
    bool isPlanar = fe->IsPlanar( mesh_, corn );

    if ( isPlanar == false )
    {
      inGlobal = true;
      warpFac.push_back(i);
      std::cout<<"encountered a WARPED side\n";
    }
  }

  for(plain_facet_set::const_iterator i=facete.begin();i!=facete.end();i++)
  {
    Facet *fe = *i;

    std::vector<std::vector<double> > cornersLocal = fe->CornerPointsLocal(elem1_);
    FacetIntegration faee1(fe,elem1_,position_,false,false);

    std::vector<double> RefPlaneTemp = faee1.equation_plane(cornersLocal);
    eqnAllFacets[i-facete.begin()] = RefPlaneTemp;

    // consider only facet whose x-direction normal componenet is non-zero
    if( fabs(RefPlaneTemp[0])>TOL_EQN_PLANE )
    {

      if( warpFac.size() > 0 ) // if there are warped facets that are not yet processed
      {
         if( i == warpFac[0] )
         {
           // reference plane cant be defined over warped facet since the facet itself is not in a plane
           facetIterator.push_back(i);
           warpFac.erase( warpFac.begin() );
           continue;
         }
      }

      // store the non-cut facet as reference
      if( isRef_==false && !fe->OnCutSide() && fabs(RefPlaneTemp[0])>REF_PLANE_DIRDIV )
      {
        RefPlaneEqn = RefPlaneTemp;
        IteratorRefFacet = i;
        isRef_ = true;
        continue;
      }
      // as far as possible, take cut side as reference
      // because it it possible to delet one or more from facetList
      if( RefOnCutSide==false && fe->OnCutSide() && fabs(RefPlaneTemp[0])>REF_PLANE_DIRDIV )
      {
        bool addRef = false;

        // if any cut side is on the plane x=C definitely this can be a reference facet
        // this produces all internal points within background element
        if( fabs(RefPlaneTemp[1])<TOL_EQN_PLANE && fabs(RefPlaneTemp[2])<TOL_EQN_PLANE )
        {
          addRef = true;
        }

        // for any other cut facet, make sure when it is extended over the whole
        // volumecell, it falls within the background element --> inside internal points
        else
        {
          double x1=0.0,x2=0.0,x3=0.0,x4=0.0;
          double y1=-1.0,y2=1.0,z1=-1.0,z2=1.0;
          x1 = (RefPlaneTemp[3]-RefPlaneTemp[2]*z1-RefPlaneTemp[1]*y1)/RefPlaneTemp[0];
          x2 = (RefPlaneTemp[3]-RefPlaneTemp[2]*z1-RefPlaneTemp[1]*y2)/RefPlaneTemp[0];
          x3 = (RefPlaneTemp[3]-RefPlaneTemp[2]*z2-RefPlaneTemp[1]*y1)/RefPlaneTemp[0];
          x4 = (RefPlaneTemp[3]-RefPlaneTemp[2]*z2-RefPlaneTemp[1]*y2)/RefPlaneTemp[0];

          //TODO: This is specific to hex8 background element. Extend this for other elements also
          if( fabs(x1) < (1.0+1e-8) && fabs(x2) < (1.0+1e-8)  &&
              fabs(x3) < (1.0+1e-8) && fabs(x4) < (1.0+1e-8) )
          {
            addRef = true;
          }
        }

        if( addRef )
        {
          if( isRef_ ) // push already existing reference facet from non-cut side
            facetIterator.push_back(IteratorRefFacet);
          RefPlaneEqn = RefPlaneTemp;
          IteratorRefFacet = i;
          isRef_ = true;
          RefOnCutSide = true;
          continue;
        }
      }
      facetIterator.push_back(i); // if not a referece side, push this to facet consideration
    }
  }

  // if no reference side found --> no element side and no cut side with x=C
  // this means we create a plane x=0 and assumes this as reference plane
  if( isRef_==false )
  {
    RefPlaneEqn[0] = 1.0;
    for( unsigned i=1;i<4;i++ )
      RefPlaneEqn[i] = 0.0;
  }

  // if a1x+a2y+a3z=a4 is the equation of reference plane and
  //    b1x+b2y+b3z=b4 is equation of the considered facet
  // if (a1/a4==b1/b4 && a2/a4==b2/b4 && a3/a4==b3/b4 ) then both reference and the
  //   considered facet are in the same plane, so delete this facet
  if( RefOnCutSide )
  {
    for( unsigned i=0;i<facetIterator.size();i++ )
    {
      plain_facet_set::const_iterator iter = facetIterator[i];

      double facetx = eqnAllFacets[iter-facete.begin()][0];
      double facety = eqnAllFacets[iter-facete.begin()][1];
      double facetz = eqnAllFacets[iter-facete.begin()][2];
      double facetRhs = eqnAllFacets[iter-facete.begin()][3];

      if( fabs(RefPlaneEqn[3])>TOL_EQN_PLANE && fabs(facetRhs)>TOL_EQN_PLANE ) // planes for which ax+by+cz=d
      {
        if( fabs(RefPlaneEqn[0]/RefPlaneEqn[3]-facetx/facetRhs)<TOL_EQN_PLANE &&
            fabs(RefPlaneEqn[1]/RefPlaneEqn[3]-facety/facetRhs)<TOL_EQN_PLANE &&
            fabs(RefPlaneEqn[2]/RefPlaneEqn[3]-facetz/facetRhs)<TOL_EQN_PLANE )
        {
          facetIterator.erase( facetIterator.begin()+i );
          i--;
        }
      }
      else                            // planes for which ax+by+cz=0
      {
        if( fabs(RefPlaneEqn[1]/RefPlaneEqn[0]-facety/facetx)<TOL_EQN_PLANE &&
            fabs(RefPlaneEqn[2]/RefPlaneEqn[0]-facetz/facetx)<TOL_EQN_PLANE &&
            fabs(RefPlaneEqn[3]/RefPlaneEqn[0]-facetRhs/facetx)<TOL_EQN_PLANE )
        {
          facetIterator.erase( facetIterator.begin()+i );
          i--;
        }
      }
    }
  }

}

/*--------------------------------------------------------------------------------------------------------------*
                   Geometry of volumecell and main Gauss pts for visualization                     sudhakar 04/12
*---------------------------------------------------------------------------------------------------------------*/
void GEO::CUT::DirectDivergence::DivengenceCellsGMSH( const DRT::UTILS::GaussIntegration & gpv,
                                                      const std::vector<DRT::UTILS::GaussIntegration> intGRule )
{
  static int sideno = 0;
  sideno++;

  std::stringstream str;
  str << "divergenceCells" << sideno << ".pos";
  std::ofstream file( str.str().c_str() );

  volcell_->DumpGmsh( file );

  //-----------------------------------------------------------------------
  // write main Gauss points
  file<<"Geometry.PointSize=8.0;\n";      // Increase the point size
  file<<"View \"Main points \" {\n";
  for ( DRT::UTILS::GaussIntegration::iterator iquad=gpv.begin(); iquad!=gpv.end(); ++iquad )
  {
    const LINALG::Matrix<3,1> etaFacet( iquad.Point() );
    file<<"SP("<<etaFacet(0,0)<<","<<etaFacet(1,0)<<","<<etaFacet(2,0)<<","<<"1"<<"){0.0};"<<std::endl;
  }
  file<<"};\n";
  file<<"View[PostProcessing.NbViews-1].ColorTable = { {0,100,0} };\n"; // Changing color to red
  file<<"View[PostProcessing.NbViews-1].Light=0;\n";  // Disable the lighting
  file<<"View[PostProcessing.NbViews-1].ShowScale=0;\n";  // Disable legend
  file<<"View[PostProcessing.NbViews-1].PointSize = 6.0;"; // increase point size

  //-----------------------------------------------------------------------
  // write internal Gauss points
  file<<"Geometry.PointSize=8.0;\n";      // Increase the point size
  file<<"View \"Internal points \" {\n";
  int nu = 0;
  for ( DRT::UTILS::GaussIntegration::iterator iquad=gpv.begin(); iquad!=gpv.end(); ++iquad )
  {
    DRT::UTILS::GaussIntegration gi = intGRule[nu];
    for ( DRT::UTILS::GaussIntegration::iterator iqu=gi.begin(); iqu!=gi.end(); ++iqu )
    {
      const LINALG::Matrix<3,1> eta( iqu.Point() );
      file<<"SP("<<eta(0,0)<<","<<eta(1,0)<<","<<eta(2,0)<<","<<"1"<<"){0.0};"<<std::endl;
    }
    nu++;
  }
  file<<"};\n";
  file<<"View[PostProcessing.NbViews-1].ColorTable = { {0,0,0} };\n"; // Changing color to red
  file<<"View[PostProcessing.NbViews-1].Light=0;\n";  // Disable the lighting
  file<<"View[PostProcessing.NbViews-1].ShowScale=0;\n";  // Disable legend
  file<<"View[PostProcessing.NbViews-1].PointSize = 4.0;"; // increase point size

  //---------------------------------------------------------------------------------------
  // write the coordinates of the reference facet, and change its color for better visualization
  file<<"View \"Ref Facet \" {\n";
  //ref plane is on a facet
  if( isRef_ )
  {
    std::vector<std::vector<double> > corners = refFacet_->CornerPointsLocal(elem1_);
    for( unsigned i=0;i<corners.size();i++ )
    {
      const std::vector<double> coords1 = corners[i];
      const std::vector<double> coords2 = corners[(i+1)%corners.size()];
      file<<"SL("<<coords1[0]<<","<<coords1[1]<<","<<coords1[2]<<","<<
          coords2[0]<<","<<coords2[1]<<","<<coords2[2]<<")"<<"{0,0};\n";
    }
  }
  else // ref plane is x=0
  {
    file<<"SL(0,-1,-1,0,1,-1)"<<"{0,0};\n";
    file<<"SL(0,1,-1,0,1,1)"<<"{0,0};\n";
    file<<"SL(0,1,1,0,-1,1)"<<"{0,0};\n";
    file<<"SL(0,-1,1,0,-1,-1)"<<"{0,0};\n";
  }
  file<<"};\n";
  file<<"View[PostProcessing.NbViews-1].ColorTable = { {255,0,0} };\n"; // Changing color to red
  file<<"View[PostProcessing.NbViews-1].Light=0;\n";  // Disable the lighting
  file<<"View[PostProcessing.NbViews-1].ShowScale=0;\n";  // Disable legend
  file<<"View[PostProcessing.NbViews-1].LineWidth = 4.0;"; // increase line width*/
}

/*--------------------------------------------------------------------------------------------------------------*
     Compute the volume of the considered cell by integrating 1 using the Gauss rule obtained.     sudhakar 04/12
     Then the volume in local coordinates is converted to global coordinate value
*---------------------------------------------------------------------------------------------------------------*/
void GEO::CUT::DirectDivergence::DebugVolume( const DRT::UTILS::GaussIntegration & gpv,
                                              const std::vector<double> &RefPlaneEqn,
                                              const std::vector<DRT::UTILS::GaussIntegration> intGRule )
{

  int numint=0;
  double TotalInteg=0.0;

  for ( DRT::UTILS::GaussIntegration::iterator iquad=gpv.begin(); iquad!=gpv.end(); ++iquad )
  {
    const LINALG::Matrix<3,1> etaFacet( iquad.Point() );
    const double weiFacet = iquad.Weight();

    double integVal = 0.0;
    DRT::UTILS::GaussIntegration gi = intGRule[numint];

    for ( DRT::UTILS::GaussIntegration::iterator iqu=gi.begin(); iqu!=gi.end(); ++iqu )
    {
      double weight = iqu.Weight();
      integVal += 1.0*weight; //Integration of 1.0 since volume is computed
    }
    TotalInteg += integVal*weiFacet;
    numint++;
  }

  //set the volume of this volumecell
  //the volume from local coordinates is converted in terms of global coordinates
  double volGlobal=0.0;
  switch ( elem1_->Shape() )
  {
    case DRT::Element::hex8:
    {
      volGlobal = elem1_->ScalarFromLocalToGlobal<DRT::Element::hex8>(TotalInteg,"LocalToGlobal");
      break;
    }
    case DRT::Element::tet4:
    {
      volGlobal = elem1_->ScalarFromLocalToGlobal<DRT::Element::tet4>(TotalInteg,"LocalToGlobal");
      break;
    }
    case DRT::Element::wedge6:
    {
      volGlobal = elem1_->ScalarFromLocalToGlobal<DRT::Element::wedge6>(TotalInteg,"LocalToGlobal");
      break;
    }
    case DRT::Element::pyramid5:
    {
      volGlobal = elem1_->ScalarFromLocalToGlobal<DRT::Element::pyramid5>(TotalInteg,"LocalToGlobal");
      break;
    }
    default:
      throw std::runtime_error( "unsupported integration cell type" );
  }

  if( volGlobal<0.0 || TotalInteg<0.0 )
  {
    std::cout<<"volume in local coordinates = "<<TotalInteg<<"\t volume in global coordinates = "<<volGlobal<<"\n";
    dserror("negative volume predicted by the DirectDivergence integration rule");
  }

#ifdef DEBUGCUTLIBRARY //check the volume with the moment fitting and check the values
  VolumeIntegration vi( volcell_, elem1_, volcell_->Position(), 1);
  Epetra_SerialDenseVector volMom = vi.compute_rhs_moment();
  volMom(0) = volcell_->Volume();

  if( fabs(volGlobal-volMom(0))>1e-6 )
  {
    std::cout<<"comparison of volume prediction\n";
    std::cout<<std::setprecision(15)<<volGlobal<<"\t"<<volMom(0)<<"\n";
    dserror("volume prediction is wrong");
  }
#endif

  volcell_->SetVolume(volGlobal);
}

/*--------------------------------------------------------------------------------------------------------------*
         Integrate given polynomials using the gaussian rule generated using directDivergence.   sudhakar 04/12
         Can be used for post-processing
*---------------------------------------------------------------------------------------------------------------*/
void GEO::CUT::DirectDivergence::IntegrateSpecificFuntions( const DRT::UTILS::GaussIntegration & gpv,
                                                            const std::vector<double> &RefPlaneEqn,
                                                            const std::vector<DRT::UTILS::GaussIntegration> intGRule )
{

  int numint=0;
  double TotalInteg=0.0;

  for ( DRT::UTILS::GaussIntegration::iterator iquad=gpv.begin(); iquad!=gpv.end(); ++iquad )
  {
    const LINALG::Matrix<3,1> etaFacet( iquad.Point() );
    const double weiFacet = iquad.Weight();

    double integVal = 0.0;
    DRT::UTILS::GaussIntegration gi = intGRule[numint];

    for ( DRT::UTILS::GaussIntegration::iterator iqu=gi.begin(); iqu!=gi.end(); ++iqu )
    {
      double weight = iqu.Weight();

      const LINALG::Matrix<3,1> eta( iqu.Point() );
      double xx = eta(0,0);
      double yy = eta(1,0);
      double zz = eta(2,0);

      integVal += (pow(xx,6)+xx*pow(yy,4)*zz+xx*xx*yy*yy*zz*zz+pow(zz,6))*weight; //Integration of 1.0 since volume is computed
    }
    TotalInteg += integVal*weiFacet;
    numint++;
  }

  std::cout<<std::setprecision(20)<<"the integral = "<<TotalInteg<<"\n";

}
