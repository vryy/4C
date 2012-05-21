#include "direct_divergence.H"
#include "facet_integration.H"
#include "volume_integration.H"

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

//  std::cout<<"number of facets for integration = "<<facetIterator.size()<<"\n";
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

  /******************************************************************************/
  /*std::cout<<"the corner coordinates of the reference facet\n";
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

  DRT::UTILS::GaussIntegration gi(cgp);

#if 0 //integrate specified functions using the Gaussian rule generated -- used in postprocessing
  ComputeVolume( gi, RefPlaneEqn );               //compute volume of the cell by integrating 1.0
  IntegrateSpecificFuntions( gi, RefPlaneEqn );  //integrate specific functions
#endif

#ifdef DEBUGCUTLIBRARY
  DivengenceCellsGMSH( IteratorRefFacet, facetIterator, gi );
#endif
  return cgp;
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
}

void GEO::CUT::DirectDivergence::DivengenceCellsGMSH( plain_facet_set::const_iterator& IteratorRefFacet,
                                                      std::vector<plain_facet_set::const_iterator>& facetIterator,
                                                      const DRT::UTILS::GaussIntegration & gpv )
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

  int nu = 10001;
  for ( DRT::UTILS::GaussIntegration::iterator iquad=gpv.begin(); iquad!=gpv.end(); ++iquad )
  {
    const LINALG::Matrix<3,1> etaFacet( iquad.Point() );
    file<<"Point("<<nu<<")={"<<etaFacet(0,0)<<","<<etaFacet(1,0)<<","<<etaFacet(2,0)<<","<<"1"<<"};"<<std::endl;
    nu++;
  }
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

#ifdef DEBUGCUTLIBRARY //check the volume with the moment fitting and check the values
  VolumeIntegration vi( volcell_, elem1_, volcell_->Position(), 1);
  Epetra_SerialDenseVector volMom = vi.compute_rhs_moment();
  volMom(0) = volcell_->Volume();
  std::cout<<"comparison of volume prediction\n";
  std::cout<<volGlobal<<"\t"<<volMom(0)<<"\n";
/*  if( fabs(TotalInteg-volMom(0))>1e-12 )
    dserror("volume prediction is wrong");*/
#endif

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

      weight = weight*jac;
      if( xbegin>etaFacet(0,0) )
        weight = -1*weight;

      //integVal += intGausspt*intGausspt*weight*jac; //integration of x2
      integVal += intGausspt*etaFacet(2,0)*weight; //integration of z
    }
    TotalInteg += integVal*weiFacet;
  }
  //std::cout<<"integration of z = "<<TotalInteg<<"\n";
}
