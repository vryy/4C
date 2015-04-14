/*!-----------------------------------------------------------------------------------------------*
\file cut_output.cpp

\brief Handles file writing of all cut related stuff

<pre>
Maintainer: Sudhakar
            sudhakar@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15257
</pre>
 *------------------------------------------------------------------------------------------------*/

#include<iosfwd>
#include<vector>

#include "cut_output.H"
#include "cut_element.H"
#include "cut_line.H"
#include "cut_edge.H"

/*--------------------------------------------------------------------------------------*
 * Write GMSH output of given element                                           sudhakar 03/14
 *--------------------------------------------------------------------------------------*/
void GEO::CUT::OUTPUT::GmshElementDump( std::ofstream & file, Element * ele )
{
  const std::vector<Node*> & nodes = ele->Nodes();
  char elementtype;
  switch ( nodes.size() )
  {
  case 8:
    elementtype = 'H';
    break;
  case 4:
    elementtype = 'S';
    break;
  case 6:
    elementtype = 'I';
    break;
  default:
    throw std::runtime_error( "unknown element type in GmshElementDump" );
  }
  GmshElementDump( file, nodes, elementtype );
}

/*--------------------------------------------------------------------------------------*
 * Write GMSH output of given element                                           sudhakar 03/14
 *--------------------------------------------------------------------------------------*/
void GEO::CUT::OUTPUT::GmshElementDump( std::ofstream & file, const std::vector<GEO::CUT::Node*> & nodes, char elementtype )
{
  file << "S" << elementtype
       << "(";
  for ( std::vector<GEO::CUT::Node*>::const_iterator i=nodes.begin(); i!=nodes.end(); ++i )
  {
    GEO::CUT::Node * n = *i;
    double x[3];
    n->Coordinates( x );
    if ( i!=nodes.begin() )
      file << ",";
    file << x[0] << "," << x[1] << "," << x[2];
  }
  file << "){";
  for ( std::vector<GEO::CUT::Node*>::const_iterator i=nodes.begin(); i!=nodes.end(); ++i )
  {
    GEO::CUT::Node * n = *i;
    GEO::CUT::Point * p = n->point();
    if ( i!=nodes.begin() )
      file << ",";
    file << p->Position();
    //file << n->DofSets().size();
  }
  file << "};\n";
}

/*--------------------------------------------------------------------------------------*
 * Write GMSH output of given side                                              sudhakar 03/14
 *--------------------------------------------------------------------------------------*/
void GEO::CUT::OUTPUT::GmshSideDump( std::ofstream & file, const Side* s )
{
  const std::vector<Node*> & nodes = s->Nodes();
  char elementtype;
  switch ( nodes.size() )
  {
  case 3:
    elementtype = 'T';
    break;
  case 4:
    elementtype = 'Q';
    break;
  default:
    throw std::runtime_error( "unknown element type in GmshSideDump" );
  }
  GmshElementDump( file, nodes, elementtype );
}

/*--------------------------------------------------------------------------------------*
 * Write GMSH output of element along with all its cut sides                sudhakar 03/14
 *--------------------------------------------------------------------------------------*/
void GEO::CUT::OUTPUT::GmshCompleteCutElement( std::ofstream & file, Element * ele )
{
  // write details of background element
  file << "View \"" << "Element" << "\" {\n";
  GmshElementDump( file, ele );
  file<<"};";

  // write details of cut sides
  file << "View \"" << "Cut sides" << "\" {\n";
  const plain_side_set & cutsides = ele->CutSides();
  for( plain_side_set::const_iterator its = cutsides.begin(); its != cutsides.end(); its++ )
  {
    const Side * s = *its;
    GmshSideDump( file, s );
  }
  file<<"};";
}

/*--------------------------------------------------------------------------------------*
 * Write GMSH output of given line                                           ager 04/15
 *--------------------------------------------------------------------------------------*/
void GEO::CUT::OUTPUT::GmshLineDump( std::ofstream & file, GEO::CUT::Line*  line)
{
     file << "SL (";
     Point* p1 = line->BeginPoint();
     Point* p2 = line->EndPoint();
     file << p1->X()[0] << "," << p1->X()[1]
         << "," << p1->X()[2] << ",";
     file << p2->X()[0] << "," << p2->X()[1]
         << "," << p2->X()[2];
     file << "){";
     file << p1->Id()<< ",";
     file << p2->Id();
     file << "};\n";
}

/*--------------------------------------------------------------------------------------*
 * Write GMSH output of given edge                                           ager 04/15
 *--------------------------------------------------------------------------------------*/
void GEO::CUT::OUTPUT::GmshEdgeDump( std::ofstream & file, GEO::CUT::Edge*  edge)
{
     file << "SL (";
     Point* p1 = edge->BeginNode()->point();
     Point* p2 = edge->EndNode()->point();
     file << p1->X()[0] << "," << p1->X()[1]
         << "," << p1->X()[2] << ",";
     file << p2->X()[0] << "," << p2->X()[1]
         << "," << p2->X()[2];
     file << "){";
     file << edge->BeginNode()->Id()<< ",";
     file << edge->EndNode()->Id();
     file << "};\n";
}

/*--------------------------------------------------------------------------------------*
 * Write GMSH output of given node                                           ager 04/15
 *--------------------------------------------------------------------------------------*/
void GEO::CUT::OUTPUT::GmshNodeDump( std::ofstream & file, GEO::CUT::Node*  node)
{
     file << "SP (";
     LINALG::Matrix<3, 1> nodecoordinates;
     node->Coordinates(nodecoordinates.A());
     file << nodecoordinates(0, 0) << "," << nodecoordinates(1, 0)
         << "," << nodecoordinates(2, 0);
     file << "){";
     file << node->Id();
     file << "};\n";
}

/*--------------------------------------------------------------------------------------*
 * Write GMSH output of given point                                           ager 04/15
 *--------------------------------------------------------------------------------------*/
void GEO::CUT::OUTPUT::GmshPointDump( std::ofstream & file, GEO::CUT::Point*  point, int idx)
{
     file << "SP (";
     LINALG::Matrix<3, 1> pointcoordinates;
     point->Coordinates(pointcoordinates.A());
     file << pointcoordinates(0, 0) << "," << pointcoordinates(1, 0)
         << "," << pointcoordinates(2, 0);
     file << "){";
     file << idx;
     file << "};\n";
}


/*--------------------------------------------------------------------------------------*
 * Write GMSH output of given point                                           ager 04/15
 *--------------------------------------------------------------------------------------*/
void GEO::CUT::OUTPUT::GmshPointDump( std::ofstream & file, GEO::CUT::Point*  point)
{
  GmshPointDump(file,point, point->Position());
}

/*-------------------------------------------------------------------------------*
 * Write cuttest for this element!                                     ager 04/15
 *-------------------------------------------------------------------------------*/
void GEO::CUT::OUTPUT::GmshElementCutTest( std::ofstream & file, GEO::CUT::Element* ele)
{
  std::cout << "Write Cut Test for Element " << ele->Id() << " ... " << std::flush;
  file << "// This test was automatically generated by CUT::OUTPUT::GmshElementCutTest(), " << "\n";
  file << "// as the cut crashed for this configuration!" << "\n";
  file << "" << "\n";
  file << "#include <iostream>" << "\n";
  file << "#include <map>" << "\n";
  file << "#include <string>" << "\n";
  file << "#include <vector>" << "\n";
  file << "" << "\n";
  file << "#include \"cut_test_utils.H\"" << "\n";
  file << "" << "\n";
  file << "#include \"../../src/drt_cut/cut_side.H\"" << "\n";
  file << "#include \"../../src/drt_cut/cut_meshintersection.H\"" << "\n";
  file << "#include \"../../src/drt_cut/cut_tetmeshintersection.H\"" << "\n";
  file << "#include \"../../src/drt_cut/cut_options.H\"" << "\n";
  file << "#include \"../../src/drt_cut/cut_volumecell.H\"" << "\n";
  file << "" << "\n";
  file << "#include \"../../src/drt_fem_general/drt_utils_local_connectivity_matrices.H\"" << "\n";
  file << "" << "\n";
  file << "void test_bacigenerated_" << ele->Id() << "()" << "\n";
  file << "{" << "\n";
  file << "  GEO::CUT::MeshIntersection intersection;" << "\n";
  file << "  std::vector<int> nids;" << "\n";
  file << "" << "\n";
  file << "  int sidecount = 0;" << "\n";
  //add sides
  const plain_side_set & cutsides = ele->CutSides();
  for (plain_side_set::const_iterator i = cutsides.begin(); i != cutsides.end();++i)
  {
    file << "  {" << "\n";
    file << "    Epetra_SerialDenseMatrix tri3_xyze( 3, 3 );" << "\n";
    file << "" << "\n";
    Side * s = *i;
    const std::vector<Node*> & side_nodes = s->Nodes();
    int nodelid = -1;
    file << "    nids.clear();" << "\n";
    for (std::vector<Node*>::const_iterator j = side_nodes.begin(); j != side_nodes.end(); ++j)
    {
      nodelid++;
      Node * n = *j;
      for (int dim = 0; dim < 3; ++dim)
      {
        file << "    tri3_xyze(" << dim << "," << nodelid << ") = " << n->point()->X()[dim] << ";" << "\n";
      }
      file << "    nids.push_back( " << n->Id() << " );" << "\n";
    }
    file << "    intersection.AddCutSide( ++sidecount, nids, tri3_xyze, DRT::Element::tri3 );" << "\n";
    file << "  }" << "\n";
  }

  //add background element
  file << "  {" << "\n";
  file << "  Epetra_SerialDenseMatrix hex" << ele->Nodes().size() <<"_xyze( 3, " << ele->Nodes().size() << " );" << "\n";
  file << "" << "\n";
  file << "    nids.clear();" << "\n";
  for (uint i = 0; i < ele->Nodes().size(); ++i)
  {
    for (uint dim = 0; dim < 3; ++dim)
    {
      file << "  hex8_xyze(" << dim << "," << i << ") = " << ele->Nodes()[i]->point()->X()[dim] << ";" << "\n";
    }
    file << "  nids.push_back( " << ele->Nodes()[i]->Id() << " );" << "\n";
  }
  file << "" << "\n";
  file << "  intersection.AddElement( 1, nids, hex8_xyze, DRT::Element::hex8 );" << "\n";
  file << "" << "\n";
  file << "  intersection.Status();" << "\n";
  file << "" << "\n";
  file << "  intersection.Cut( true, INPAR::CUT::VCellGaussPts_Tessellation );" << "\n";
  file << "  }" << "\n";
  file << "" << "\n";

  //Compare Integration
  file << "  std::vector<double> tessVol,momFitVol,dirDivVol;" << "\n";
  file << "" << "\n";
  file << "  GEO::CUT::Mesh mesh = intersection.NormalMesh();" << "\n";
  file << "  const std::list<Teuchos::RCP<GEO::CUT::VolumeCell> > & other_cells = mesh.VolumeCells();" << "\n";
  file << "  for ( std::list<Teuchos::RCP<GEO::CUT::VolumeCell> >::const_iterator i=other_cells.begin();" << "\n";
  file << "        i!=other_cells.end();" << "\n";
  file << "        ++i )" << "\n";
  file << "  {" << "\n";
  file << "    GEO::CUT::VolumeCell * vc = &**i;" << "\n";
  file << "    tessVol.push_back(vc->Volume());" << "\n";
  file << "  }" << "\n";
  file << "" << "\n";
  file << "  intersection.Status();" << "\n";
  file << "" << "\n";
  file << "  for ( std::list<Teuchos::RCP<GEO::CUT::VolumeCell> >::const_iterator i=other_cells.begin();" << "\n";
  file << "        i!=other_cells.end();" << "\n";
  file << "        ++i )" << "\n";
  file << "  {" << "\n";
  file << "    GEO::CUT::VolumeCell * vc = &**i;" << "\n";
  file << "    vc->MomentFitGaussWeights(vc->ParentElement(),mesh,true,INPAR::CUT::BCellGaussPts_Tessellation);" << "\n";
  file << "    momFitVol.push_back(vc->Volume());" << "\n";
  file << "  }" << "\n";
  file << "" << "\n";
  file << "  for ( std::list<Teuchos::RCP<GEO::CUT::VolumeCell> >::const_iterator i=other_cells.begin();" << "\n";
  file << "           i!=other_cells.end();" << "\n";
  file << "           ++i )" << "\n";
  file << "   {" << "\n";
  file << "     GEO::CUT::VolumeCell * vc = &**i;" << "\n";
  file << "     vc->DirectDivergenceGaussRule(vc->ParentElement(),mesh,true,INPAR::CUT::BCellGaussPts_DirectDivergence);" << "\n";
  file << "     dirDivVol.push_back(vc->Volume());" << "\n";
  file << "   }" << "\n";
  file << "" << "\n";
  file << "  std::cout<<\"the volumes predicted by\\n tessellation \\t MomentFitting \\t DirectDivergence\\n\";" << "\n";
  file << "  for(unsigned i=0;i<tessVol.size();i++)" << "\n";
  file << "  {" << "\n";
  file << "    std::cout<<tessVol[i]<<\"\\t\"<<momFitVol[i]<<\"\\t\"<<dirDivVol[i]<<\"\\n\";" << "\n";
  file << "    if( fabs(tessVol[i]-momFitVol[i])>1e-9 || fabs(dirDivVol[i]-momFitVol[i])>1e-9 )" << "\n";
  file << "      dserror(\"volume predicted by either one of the method is wrong\");" << "\n";
  file << "  }" << "\n";
  file << "}" << "\n";
  std::cout << "done " << std::endl;
}



