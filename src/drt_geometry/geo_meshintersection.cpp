/*!-----------------------------------------------------------------------------------------------*
\file geo_meshintersection.cpp

\brief class that provides to set up a mesh cut based on surface meshes

<pre>
Maintainer: Benedikt Schott
            schott@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>
 *------------------------------------------------------------------------------------------------*/
#include <Teuchos_TimeMonitor.hpp>

#include "../drt_lib/drt_discret.H"
#include "../drt_cut/cut_meshintersection.H"

#include "../drt_geometry/geo_intersection.H"

#include "geo_meshintersection.H"

/// Constructor
GEO::CutWizardMesh::CutWizardMesh( DRT::Discretization & dis, int numcutmesh )
:   CutWizard( dis )
{
  meshintersection_ = Teuchos::rcp( new GEO::CUT::MeshIntersection( numcutmesh, myrank_ ) );
  mesh_ = meshintersection_;

}

void GEO::CutWizardMesh::AddCutSide( int mi, DRT::Element * ele, const Epetra_SerialDenseMatrix & xyze )
{
  const int numnode = ele->NumNode();
  const int * nodeids = ele->NodeIds();

  std::vector<int> nids( nodeids, nodeids+numnode );
  meshintersection_->AddCutSide( ele->Id(), nids, xyze, ele->Shape(), mi );
}

void GEO::CutWizardMesh::AddElement( DRT::Element * ele, const Epetra_SerialDenseMatrix & xyze)
{
  const int numnode = ele->NumNode();
  const int * nodeids = ele->NodeIds();

  std::vector<int> nids( nodeids, nodeids+numnode );
  meshintersection_->AddElement( ele->Id(), nids, xyze, ele->Shape() );
}

GEO::CUT::SideHandle * GEO::CutWizardMesh::GetCutSide( int sid, int mi )
{
  return meshintersection_->GetCutSide( sid, mi );
}

/*------------------------------------------------------------------------------------------------*
 * build the bounding volume tree for the collision detection in the context of the selfcut       *
 *                                                                                    wirtz 09/14 *
 *------------------------------------------------------------------------------------------------*/
void GEO::CutWizardMesh::BuildBVTree()
{
  meshintersection_->BuildBVTree();
}

/*------------------------------------------------------------------------------------------------*
 * build the static search tree for the collision detection                           wirtz 08/14 *
 *------------------------------------------------------------------------------------------------*/
void GEO::CutWizardMesh::BuildStaticSearchTree()
{
  meshintersection_->BuildStaticSearchTree();
}

/*------------------------------------------------------------------------------------------------*
 * find positions                                                                    schott 11/14 *
 *------------------------------------------------------------------------------------------------*/
void GEO::CutWizardMesh::FindNodePositions()
{
  GEO::CUT::Mesh & m = mesh_->NormalMesh();
  m.FindNodePositions();
}

/*------------------------------------------------------------------------------------------------*
 *cut routine for standard non-parallel framework (only for cuttest)                 schott 03/12 *
 *------------------------------------------------------------------------------------------------*/
void GEO::CutWizardMesh::Cut(
    bool include_inner,
    INPAR::CUT::VCellGaussPts VCellgausstype,  //!< Gauss point generation method for Volumecell
    INPAR::CUT::BCellGaussPts BCellgausstype,  //!< Gauss point generation method for Boundarycell
    bool screenoutput
)
{
  meshintersection_->Cut( include_inner, VCellgausstype, BCellgausstype,screenoutput );
}
