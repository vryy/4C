/*!-----------------------------------------------------------------------------------------------*
\file cut_output.H

\brief Handles file writing of all cut related stuff

<pre>
Maintainer: Sudhakar
            sudhakar@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15257
</pre>
 *------------------------------------------------------------------------------------------------*/

#ifndef CUTOUTPUT_H_
#define CUTOUTPUT_H_

#include "cut_point.H"

namespace GEO {
namespace CUT {

class Element;
class Side;
class Node;
class Point;
class Line;
class Edge;

namespace OUTPUT {


/*!
\brief Write output of background element geometry
 */
void GmshElementDump( std::ofstream & file, Element * ele );

/*!
\brief Write output of background element geometry
 */
void GmshElementDump( std::ofstream & file, const std::vector<Node*> & nodes, char elementtype );

/*!
\brief Write output of a side
 */
void GmshSideDump( std::ofstream & file, const Side* s );

/*!
\brief Write output of a Triside
 */
void GmshTriSideDump( std::ofstream & file, std::vector<Point*> points );

/*!
\brief Write output of the background element and all the cut sides corresponding to this element
 */
void GmshCompleteCutElement( std::ofstream & file, Element * ele );

/*!
\brief Write output of a line
 */
void GmshLineDump( std::ofstream & file, GEO::CUT::Line*  line);

/*!
\brief Write output of a edge
 */
void GmshEdgeDump( std::ofstream & file, GEO::CUT::Edge*  edge);

/*!
\brief Write output of a node
 */
void GmshNodeDump( std::ofstream & file, GEO::CUT::Node*  node);

/*!
\brief Write output of a point with special idx
 */
void GmshPointDump( std::ofstream & file, GEO::CUT::Point*  point,int idx);

/*!
\brief Write output of a point
 */
void GmshPointDump( std::ofstream & file, GEO::CUT::Point*  point);

/*!
\brief Write cuttest for this element!
*/
void GmshElementCutTest( std::ofstream & file, GEO::CUT::Element* ele, bool haslevelsetside=false);

/*!
 \brief Generate filename for gmsh output with specific ending
 */
std::string GenerateGmshOutputFilename(const std::string& filename_tail);


}

} /* namespace CUT */
} /* namespace GEO */
#endif /* CUTOUTPUT_H_ */
