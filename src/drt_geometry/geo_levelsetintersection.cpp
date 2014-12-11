/*!-----------------------------------------------------------------------------------------------*
\file geo_levelsetintersection.cpp

\brief class that providesthe wizard for the functionality for a mesh cut based on a level set field

<pre>
Maintainer: Benedikt Schott
            schott@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>
 *------------------------------------------------------------------------------------------------*/
#include <Teuchos_TimeMonitor.hpp>

#include "../drt_lib/drt_discret.H"
#include "../drt_cut/cut_levelsetintersection.H"

#include "geo_levelsetintersection.H"

/// Constructor
GEO::CutWizardLevelSet::CutWizardLevelSet( const DRT::Discretization & dis )
:   CutWizard( dis )
{
  levelsetintersection_ = Teuchos::rcp( new GEO::CUT::LevelSetIntersection( myrank_ ) );
  mesh_ = levelsetintersection_;

}

/*------------------------------------------------------------------------------------------------*
 * check if element is cut and add element to cut libraries if this is the case.     winter 08/14 *
 * myphinp needs to on the node map for the element. This has to be done before the call.         *
 *------------------------------------------------------------------------------------------------*/
void GEO::CutWizardLevelSet::AddElement(DRT::Element * ele, std::vector<double> myphinp, bool include_inner)
{
  const int numnode = ele->NumNode();
  const DRT::Node * const * nodes = ele->Nodes();
  const int * nodeids = ele->NodeIds();

  Epetra_SerialDenseMatrix xyze( 3, numnode );

  for ( int i=0; i < numnode; ++i )
  {
    const DRT::Node & node = *nodes[i];
    std::copy( node.X(), node.X()+3, &xyze( 0, i ) );
  }

  std::vector<int> nids( nodeids, nodeids+numnode );
  //If include_inner == false then add elements with negative level set values to discretization.
  levelsetintersection_->AddElement( ele->Id(), nids, xyze, ele->Shape(), &myphinp[0], !include_inner );
}



void GEO::CutWizardLevelSet::FindNodePositions()
{
  GEO::CUT::Mesh & m = mesh_->NormalMesh();
  m.FindLSNodePositions();
}

/*------------------------------------------------------------------------------------------------*
 *cut routine for standard non-parallel framework (only for cuttest)                 schott 03/12 *
 *------------------------------------------------------------------------------------------------*/
void GEO::CutWizardLevelSet::Cut(
    bool include_inner,
    INPAR::CUT::VCellGaussPts VCellgausstype,  //!< Gauss point generation method for Volumecell
    INPAR::CUT::BCellGaussPts BCellgausstype,  //!< Gauss point generation method for Boundarycell
    bool screenoutput
)
{
  levelsetintersection_->Cut( include_inner, screenoutput );
}
