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
