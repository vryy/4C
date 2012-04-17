#include "direct_divergence.H"
#include "facet_integration.H"

/*-------------------------------------------------------------------------------------------------------------------*
  Create integration points on the facets of the volumecell by triangulating the facets
  A reference facet is identified on which integration weights are set to zero                        Sudhakar 04/12
  This ref. facet is used to construct the modified integrand in Fluid integration part(???)
*--------------------------------------------------------------------------------------------------------------------*/
Teuchos::RCP<DRT::UTILS::GaussPoints> GEO::CUT::DirectDivergence::VCIntegrationRule( std::vector<double> &RefPlaneEqn )
{
  std::vector<plain_facet_set::const_iterator> facetIterator; //iterators of facets which need to be considered for integration rule
  plain_facet_set::const_iterator IteratorRefFacet;           //iterator for the reference facets

  ListFacets( facetIterator, RefPlaneEqn, IteratorRefFacet );

  std::cout<<"number of facets for integration = "<<facetIterator.size()<<"\n";
  if( facetIterator.size()==0 )
    dserror( "x-component normal is zero on all the facets? It should not be." );
//  dserror("break tempo");

  Teuchos::RCP<DRT::UTILS::CollectedGaussPoints> cgp = Teuchos::rcp( new DRT::UTILS::CollectedGaussPoints(0) );
  /*Teuchos::RCP<DRT::UTILS::GaussPointsComposite> gpc =
        Teuchos::rcp( new DRT::UTILS::GaussPointsComposite( 0 ) );*/

  for( unsigned i=0;i<facetIterator.size();i++ )
  {
    plain_facet_set::const_iterator iter = facetIterator[i];
    Facet * fe = *iter;
    FacetIntegration faee1(fe,elem1_,position_,false,false);

    faee1.DivergenceIntegrationRule( mesh_, cgp );
  }

  /******************************************************************************/
 /* std::cout<<"the corner coordinates of the reference facet\n";
  Facet * fe = *IteratorRefFacet;
  std::vector<std::vector<double> > cornersLocal = fe->CornerPointsLocal(elem1_);
  for( std::vector<std::vector<double> >::iterator m=cornersLocal.begin();m!=cornersLocal.end();m++ )
  {
    std::vector<double> coo = *m;
    std::cout<<coo[0]<<"\t"<<coo[1]<<"\t"<<coo[2]<<"\n";
  }

  for( unsigned i=0;i<facetIterator.size();i++ )
  {
    plain_facet_set::const_iterator iter = facetIterator[i];
    Facet * fe = *iter;
    std::vector<std::vector<double> > cornersLocal = fe->CornerPointsLocal(elem1_);
    std::cout<<"\nfacet "<<i+1<<"\n";
    for( std::vector<std::vector<double> >::iterator m=cornersLocal.begin();m!=cornersLocal.end();m++ )
    {
      std::vector<double> coo = *m;
      std::cout<<coo[0]<<"\t"<<coo[1]<<"\t"<<coo[2]<<"\n";
    }
  }*/
  /******************************************************************************/

#ifdef DEBUGCUTLIBRARY
  DivengenceCellsGMSH( IteratorRefFacet, facetIterator );
#endif

  DRT::UTILS::GaussIntegration gi(cgp);
  std::cout<<"number of Gauss points = "<<gi.NumPoints()<<"\n";
  ComputeVolume( gi, RefPlaneEqn );               //compute volume of the cell by integrating 1.0
#if 0 //integrate specified functions using the Gaussian rule generated -- used in postprocessing
  IntegrateSpecificFuntions( gi, RefPlaneEqn );  //integrate specific functions
#endif
  Teuchos::RCP<DRT::UTILS::GaussPoints> cp1;
  return cp1;
}

/*-------------------------------------------------------------------------------------------------------*
                  Identify the list of facets which need to be triangulated
                  Get the reference facet that will be used in xfluid part                  sudhakar 04/12
                  As far as possible, the reference facet is set on a cut side
*--------------------------------------------------------------------------------------------------------*/
void GEO::CUT::DirectDivergence::ListFacets( std::vector<plain_facet_set::const_iterator>& facetIterator,
                                             std::vector<double>& RefPlaneEqn,
                                             plain_facet_set::const_iterator& IteratorRefFacet )
{
#if 1
  const plain_facet_set & facete = volcell_->Facets();

  bool IsRefFacet = false,RefOnCut=false;
  std::vector<std::vector<double> > eqnAllFacets(facete.size());

  for(plain_facet_set::const_iterator i=facete.begin();i!=facete.end();i++)
  {
    Facet *fe = *i;
    std::vector<std::vector<double> > cornersLocal = fe->CornerPointsLocal(elem1_);

    FacetIntegration faee1(fe,elem1_,position_,false,false);

    std::vector<double> RefPlaneTemp = faee1.equation_plane(cornersLocal);
    eqnAllFacets[i-facete.begin()] = RefPlaneTemp;

    if( fabs(RefPlaneTemp[0])>1e-10 )
    {
 //     std::cout<<"another equation = ";
 //     std::cout<<RefPlaneTemp[0]<<"\t"<<RefPlaneTemp[1]<<"\t"<<RefPlaneTemp[2]<<"\t"<<RefPlaneTemp[3]<<"\n";
      if( IsRefFacet==false )
      {
        RefPlaneEqn = RefPlaneTemp;
        IteratorRefFacet = i;
        IsRefFacet = true;
        if( fe->OnCutSide() )
          RefOnCut = true;
      }
      else
      {
        if( fe->OnCutSide() && RefOnCut == false )
        {
          facetIterator.push_back(IteratorRefFacet); //when the ref facet is changed, the old ref facet is included for integration
          RefPlaneEqn = RefPlaneTemp;
          IteratorRefFacet = i;
          RefOnCut = true;
        }
        else
          facetIterator.push_back(i);
      }
    }
  }

/*  std::cout<<"the equation of refere plane ="<<RefPlaneEqn[0]<<"\t"<<RefPlaneEqn[1]<<"\t"<<RefPlaneEqn[2]<<"\t"
      <<RefPlaneEqn[3]<<"\n";
  std::cout<<"number of facets before erasing = "<<facetIterator.size()<<"\n";*/

  // if a1x+a2y+a3z=a4 is the equation of reference plane and
  //    b1x+b2y+b3z=b4 is equation of the considered facet
  // if (a1/a4==b1/b4 && a2/a4==b2/b4 && a3/a4==b3/b4 ) then both reference and the
  //   considered facet are in the same plane, so delete this facet
  for( unsigned i=0;i<facetIterator.size();i++ )
  {
    plain_facet_set::const_iterator iter = facetIterator[i];

    double facetx = eqnAllFacets[iter-facete.begin()][0];
    double facety = eqnAllFacets[iter-facete.begin()][1];
    double facetz = eqnAllFacets[iter-facete.begin()][2];
    double facetRhs = eqnAllFacets[iter-facete.begin()][3];
    if( fabs(RefPlaneEqn[3])>1e-10 ) //check whether it is x=0 plane
    {
      if( fabs(RefPlaneEqn[0]/RefPlaneEqn[3]-facetx/facetRhs)<1e-8 &&
          fabs(RefPlaneEqn[1]/RefPlaneEqn[3]-facety/facetRhs)<1e-8 &&
          fabs(RefPlaneEqn[2]/RefPlaneEqn[3]-facetz/facetRhs)<1e-8 )
      {
        facetIterator.erase( facetIterator.begin()+i );
        i--;
      }
    }
    else
    {
      if( fabs(facetRhs)<1e-10             &&
          fabs(RefPlaneEqn[0]-facetx)<1e-8 &&
          fabs(RefPlaneEqn[1]-facety)<1e-8 &&
          fabs(RefPlaneEqn[2]-facetz)<1e-8 )
      {
        facetIterator.erase( facetIterator.begin()+i );
        i--;
      }
    }
  }
//  std::cout<<"number of facets after erasing = "<<facetIterator.size()<<"\n";
#endif
}

void GEO::CUT::DirectDivergence::DivengenceCellsGMSH( plain_facet_set::const_iterator& IteratorRefFacet,
                                                      std::vector<plain_facet_set::const_iterator>& facetIterator )
{
  std::string filename="side";
  std::ofstream file;

  static int sideno = 0;
  sideno++;
  std::stringstream out;
  out <<"divergenceCells"<<sideno<<".pos";
  filename = out.str();
  file.open(filename.c_str());

  volcell_->DumpGmsh( file );

  /*const plain_facet_set & facete = volcell_->Facets();

  static int pointno=1,point_begin,point_end,lineno=1,line_begin,line_end,surf=1,surf_begin, surf_end;
  surf_begin = surf;
  for(plain_facet_set::const_iterator i=facete.begin();i!=facete.end();i++)
  {
     Facet *fe = *i;
     point_begin = pointno;
     const std::vector<Point*> corners = fe->CornerPoints();
     for( std::vector<Point*>::const_iterator k=corners.begin();k!=corners.end();k++ )
     {
       Point *pt = *k;
       double coords[3];
       pt->Coordinates(coords);
       file<<"Point("<<pointno<<")={"<<coords[0]<<","<<coords[1]<<","<<coords[2]<<","<<"1"<<"};"<<std::endl;
       pointno++;

     }
     point_end = pointno;
     line_begin = lineno;
     for(int i=point_begin;i!=point_end;i++)
     {
       if(i!=point_end-1)
         file<<"Line("<<lineno<<")={"<<i<<","<<i+1<<"};"<<std::endl;
       else
         file<<"Line("<<lineno<<")={"<<i<<","<<point_begin<<"};"<<std::endl;
       lineno++;
     }
     line_end = lineno;
     file<<"Line Loop("<<surf<<")={";
     for(int i=line_begin;i!=line_end;i++)
     {
       file<<i;
       if(i!=line_end-1)
         file<<",";
     }
     file<<"};\n";
     file<<"Plane Surface("<<surf<<")={"<<surf<<"};"<<std::endl;

     surf++;
  }
  surf_end = surf;
  file<<"Surface Loop("<<surf<<")={";
  for(int i=surf_begin;i<surf;i++)
  {
    file<<i;
    if(i!=(surf-1))
      file<<",";
  }
  file<<"};\n";
  file<<"Volume("<<surf<<")={"<<surf<<"};"<<std::endl;*/
}

/*--------------------------------------------------------------------------------------------------------------*
     Compute the volume of the considered cell by integrating 1 using the Gauss rule obtained.     sudhakar 04/12
     Then the volume in local coordinates is converted to global coordinate value
*---------------------------------------------------------------------------------------------------------------*/
void GEO::CUT::DirectDivergence::ComputeVolume( const DRT::UTILS::GaussIntegration & gpv,
                                                const std::vector<double> &RefPlaneEqn    )
{
  double TotalInteg=0.0;
  for ( DRT::UTILS::GaussIntegration::iterator iquad=gpv.begin(); iquad!=gpv.end(); ++iquad )
  {
    const LINALG::Matrix<3,1> etaFacet( iquad.Point() );
    const double weiFacet = iquad.Weight();

    double integVal = 0.0;
    DRT::UTILS::GaussIntegration gi( DRT::Element::line2, 5 );
    for ( DRT::UTILS::GaussIntegration::iterator iqu=gi.begin(); iqu!=gi.end(); ++iqu )
    {
      const LINALG::Matrix<1,1> eta( iqu.Point() );
      double weight = iqu.Weight();
      double xbegin = (RefPlaneEqn[3]-RefPlaneEqn[1]*etaFacet(1,0)-RefPlaneEqn[2]*etaFacet(2,0))/RefPlaneEqn[0];
      double jac = fabs(xbegin-etaFacet(0,0))*0.5;

      int fac = 1;
      if(xbegin>etaFacet(0,0))
        fac = -1;

      integVal += 1.0*weight*jac*fac; //Integration of 1.0 since volume is computed
    }
    TotalInteg += integVal*weiFacet;
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
  volcell_->SetVolume(volGlobal);

  std::cout<<"volume = "<<volGlobal<<"\n";
}

/*--------------------------------------------------------------------------------------------------------------*
         Integrate given polynomials using the gaussian rule generated using directDivergence.   sudhakar 04/12
         Can be used for post-processing
*---------------------------------------------------------------------------------------------------------------*/
void GEO::CUT::DirectDivergence::IntegrateSpecificFuntions( const DRT::UTILS::GaussIntegration & gpv,
                                                            const std::vector<double> &RefPlaneEqn    )
{
  double TotalInteg=0.0;
  for ( DRT::UTILS::GaussIntegration::iterator iquad=gpv.begin(); iquad!=gpv.end(); ++iquad )
  {
    const LINALG::Matrix<3,1> etaFacet( iquad.Point() );
    const double weiFacet = iquad.Weight();

    double integVal = 0.0;
    DRT::UTILS::GaussIntegration gi( DRT::Element::line2, 5 );
    for ( DRT::UTILS::GaussIntegration::iterator iqu=gi.begin(); iqu!=gi.end(); ++iqu )
    {
      const LINALG::Matrix<1,1> eta( iqu.Point() );
      double weight = iqu.Weight();
      double xbegin = (RefPlaneEqn[3]-RefPlaneEqn[1]*etaFacet(1,0)-RefPlaneEqn[2]*etaFacet(2,0))/RefPlaneEqn[0];
      double jac = fabs(xbegin-etaFacet(0,0))*0.5;

      double xmid = 0.5*(xbegin+etaFacet(0,0));

      double intGausspt = (xmid-xbegin)*eta(0,0)+xmid;

      int fac = 1;
      if( xbegin>etaFacet(0,0) )
        fac = -1;

      //integVal += intGausspt*intGausspt*weight*jac; //integration of x2
      integVal += etaFacet(2,0)*weight*jac*fac; //integration of z
    }
    TotalInteg += integVal*weiFacet;
  }
  //std::cout<<"integration of z = "<<TotalInteg<<"\n";
}
