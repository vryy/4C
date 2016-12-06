/*---------------------------------------------------------------------*/
/*!
\file cut_output.H

\brief Handles file writing of all cut related stuff (gmsh)

\level 2

<pre>
\maintainer Christoph Ager
            ager@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15249
</pre>

*----------------------------------------------------------------------*/

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
class Cycle;

namespace OUTPUT {

/*!
\brief Write output of background element geometry
 */
void GmshElementDump( std::ofstream & file, Element * ele , bool to_local = false);

/*!
\brief Write output of background element geometry
 */
void GmshElementDump( std::ofstream & file, const std::vector<Node*> & nodes, char elementtype , bool to_local = false, Element* ele = NULL);

/*!
\brief Write output of a side
 */
void GmshSideDump( std::ofstream & file, const Side* s, bool to_local = false, Element* ele = NULL );

/*!
\brief Write output of a Triside
 */
void GmshTriSideDump( std::ofstream & file, std::vector<Point*> points , bool to_local = false, Element* ele = NULL);

/*!
\brief Write output of a facet
 */
void GmshFacetDump( std::ofstream & file, Facet* facet, const std::string& visualizationtype = "sides" , bool print_all = false, bool to_local = false, Element* ele = NULL);

/*!
\brief Write output of a volumecell
 */
void GmshVolumecellDump( std::ofstream & file, VolumeCell* VC, const std::string& visualizationtype = "sides" , bool print_all = false, bool to_local = false, Element* ele = NULL);

/*!
\brief Write output of a cylce
 */
void GmshCycleDump( std::ofstream & file, Cycle* cycle, const std::string& visualizationtype = "sides" , bool to_local = false, Element* ele = NULL);

/*!
\brief Write output of the background element and all the cut sides corresponding to this element
 */
void GmshCompleteCutElement( std::ofstream & file, Element * ele,bool to_local = false);

/*!
\brief Write output of a line
 */
void GmshLineDump( std::ofstream & file, GEO::CUT::Line*  line, bool to_local = false, Element* ele = NULL);

void GmshLineDump( std::ofstream & file, GEO::CUT::Point* p1, GEO::CUT::Point* p2, bool to_local = false, Element* ele = NULL);

void GmshLineDump( std::ofstream & file, GEO::CUT::Point* p1, GEO::CUT::Point* p2, int idx1, int idx2, bool to_local = false, Element* ele = NULL);

/*!
\brief Write output of a edge
 */
void GmshEdgeDump( std::ofstream & file, GEO::CUT::Edge*  edge, bool to_local = false, Element* ele = NULL);

/*!
\brief Write output of a node
 */
void GmshNodeDump( std::ofstream & file, GEO::CUT::Node*  node, bool to_local = false, Element* ele = NULL);

/*!
\brief Write output of a point with special idx
 */
void GmshPointDump( std::ofstream & file, GEO::CUT::Point*  point,int idx, bool to_local = false, Element* ele = NULL);

/*!
\brief Write output of a point with point position as idx
 */
void GmshPointDump( std::ofstream & file, GEO::CUT::Point*  point, bool to_local = false, Element* ele = NULL);

/*!
\brief Write level set value on cut surface (information taken from facet).
 */
void GmshLevelSetValueDump( std::ofstream & file, Element * ele, bool dumpnodevalues=false, bool to_local = false );

/*!
\brief Write level set gradient on cut surface (information taken from facet).
 */
void GmshLevelSetGradientDump( std::ofstream & file, Element * ele , bool to_local = false);

/*!
\brief Write level set value on cut surface (information taken from facet).
 */
void GmshLevelSetValueZeroSurfaceDump( std::ofstream & file, Element * ele, bool to_local = false);

/*!
 * Write Level Set Gradient Orientation of Boundary-Cell Normal and LevelSet
 */
void GmshLevelSetOrientationDump( std::ofstream & file, Element * ele , bool to_local = false);

/*!
\brief Write Eqn of plane normal for all facets (used for DirectDivergence).
*/
void GmshEqnPlaneNormalDump(std::ofstream & file, Element * ele, bool normalize = false , bool to_local = false);
void GmshEqnPlaneNormalDump(std::ofstream & file, Facet * facet, bool normalize = false, bool to_local = false, Element* ele = NULL);
/*!
\brief Get equation of plane as implemented in DirectDivergence routine.
 */
std::vector<double> GetEqOfPlane(std::vector<Point*> pts);

/*!
\brief Simplify output of for normal output options.
 */
void GmshScalar(std::ofstream & file, LINALG::Matrix<3,1> coord, double scalar, bool to_local = false, Element* ele = NULL);
void GmshVector(std::ofstream & file, LINALG::Matrix<3,1> coord, std::vector<double> vector, bool normalize=false, bool to_local = false, Element* ele = NULL);

/*!
\brief Write cuttest for this element!
*/
void GmshElementCutTest( std::ofstream & file, GEO::CUT::Element* ele, bool haslevelsetside=false);

/*!
 \brief Generate filename for gmsh output with specific ending
 */
std::string GenerateGmshOutputFilename(const std::string& filename_tail);

/*!
 \brief Write new Section in Gmsh file (eventually end section from before...)
 */
void GmshNewSection(std::ofstream & file, const std::string & section, bool first_endsection = false);

/*!
 \brief End Section in Gmsh file
 */
void GmshEndSection(std::ofstream & file, bool close_file = false);


//void GmshObjectDump2(std::ofstream & file, GEO::CUT::Facet* f)
//{
//  return;
//};//GmshFacetDump(file,f);}

/*!
 \brief Writes the whole container into a gmsh file!
 */
//template<class T>
//void GmshPrintContainer(std::ofstream & file, const std::string & section, std::vector<T> container)
//{
//  GmshNewSection(file, section);
//  for (typename std::vector<T>::iterator t = container.begin(); t!=container.end(); ++t)
//    //GmshObjectDump( file, (T)(*t));
//  GmshEndSection(file);
//}
/*!
 \brief Write Coordinates in Gmsh file (for internal use)
 //to_local ... transform to local coordinates of the ele?
 */
void GmshWriteCoords(std::ofstream & file, std::vector<double> coord, bool to_local = false, Element* ele = NULL);

void GmshWriteCoords( std::ofstream & file, LINALG::Matrix<3,1> coord, bool to_local = false, Element* ele = NULL);

void GmshWriteCoords(std::ofstream & file, Node* node, bool to_local = false, Element* ele = NULL);

void GmshWriteCoords(std::ofstream & file, Point* point, bool to_local = false, Element* ele = NULL);
}

} /* namespace CUT */
} /* namespace GEO */
#endif /* CUTOUTPUT_H_ */
