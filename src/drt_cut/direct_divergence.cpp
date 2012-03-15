#include "direct_divergence.H"
#include "facet_integration.H"

//Create integration points on the facets of the volumecell
Teuchos::RCP<DRT::UTILS::GaussPoints> GEO::CUT::DirectDivergence::VCIntegrationRule(std::vector<double> &RefPlaneEqn)
{
  //const plain_facet_set & facete = volcell_->Facets();

  std::vector<plain_facet_set::const_iterator> facetIterator;
  plain_facet_set::const_iterator IteratorRefFacet;

  //Try to call the new function instead, and if necessary delete the unused variables

  /*bool IsRefFacet = false;
  for(plain_facet_set::const_iterator i=facete.begin();i!=facete.end();i++)
  {
    Facet *fe = *i;
    std::vector<std::vector<double> > cornersLocal = fe->CornerPointsLocal(elem1_);
    FacetIntegration faee1(fe,elem1_,position_,false,false);
    std::vector<double> RefPlaneTemp = faee1.equation_plane(cornersLocal);

    if( fabs(RefPlaneTemp[0])>0.0000001 )
    {
      if( IsRefFacet==false )
      {
        RefPlaneEqn = RefPlaneTemp;
        IteratorRefFacet = i;
        IsRefFacet = true;
      }
      else
      {
        if( fe->OnCutSide() )
        {
          facetIterator.push_back(IteratorRefFacet);
          RefPlaneEqn = RefPlaneTemp;
          IteratorRefFacet = i;
        }
        else
          facetIterator.push_back(i);
      }
    }
    else
      continue;
  }*/

  for(unsigned i=0;i<facetIterator.size();i++)
  {
    plain_facet_set::const_iterator iter = facetIterator[i];
    Facet * fe = *iter;
    FacetIntegration faee1(fe,elem1_,position_,false,false);

    faee1.DivergenceIntegrationRule(mesh_);

    std::vector<std::vector<double> > cornersLocal = fe->CornerPointsLocal(elem1_);
    /*std::cout<<"The corner points in local coordinates \n";
    for(unsigned j=0;j<cornersLocal.size();j++)
    {
      std::vector<double> a = cornersLocal[j];
      std::cout<<a[0]<<"\t"<<a[1]<<"\t"<<a[2]<<"\n";
    }*/
  }

  /*std::cout<<"the equation of the reference plane is "<<RefPlaneEqn[0]<<"x+"<<
      RefPlaneEqn[1]<<"y+"<<RefPlaneEqn[2]<<"z = "<<RefPlaneEqn[3]<<"\n";*/
  dserror("one volumecell done");

  Teuchos::RCP<DRT::UTILS::GaussPoints> gp1;
  return gp1;
}

//Identify the list of facets which need to be triangulated, and also get the reference facet that will be used in xfluid part
void GEO::CUT::DirectDivergence::ListFacets( std::vector<plain_facet_set::const_iterator>& facetIterator,
                                             std::vector<double>& RefPlaneEqn )
{
  const plain_facet_set & facete = volcell_->Facets();

  plain_facet_set::const_iterator IteratorRefFacet;
  bool IsRefFacet = false;
  for(plain_facet_set::const_iterator i=facete.begin();i!=facete.end();i++)
  {
    Facet *fe = *i;
    std::vector<std::vector<double> > cornersLocal = fe->CornerPointsLocal(elem1_);
    FacetIntegration faee1(fe,elem1_,position_,false,false);

    std::vector<double> RefPlaneTemp = faee1.equation_plane(cornersLocal);

    if( fabs(RefPlaneTemp[0])>0.0000001 )
    {
      if( IsRefFacet==false )
      {
        RefPlaneEqn = RefPlaneTemp;
        IteratorRefFacet = i;
        IsRefFacet = true;
      }
      else
      {
        if( fe->OnCutSide() )
        {
          facetIterator.push_back(IteratorRefFacet);
          RefPlaneEqn = RefPlaneTemp;
          IteratorRefFacet = i;
        }
        else
          facetIterator.push_back(i);
      }
    }
    else
      continue;
  }
}
