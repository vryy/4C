/*---------------------------------------------------------------------*/
/*! \file

\brief Handles file writing of all cut related stuff (gmsh)

\level 2

\maintainer Martin Kronbichler
            kronbichler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15249

*----------------------------------------------------------------------*/

#ifndef CUTOUTPUT_H_
#define CUTOUTPUT_H_

#include "cut_point.H"
#include "cut_boundarycell.H"

namespace GEO
{
  namespace CUT
  {
    class Element;
    class Side;
    class Node;
    class Point;
    class Line;
    class Edge;
    class Cycle;

    namespace OUTPUT
    {
      char GmshElementType(DRT::Element::DiscretizationType shape);

      /** \brief Write a Gmsh output file, which contains only the given volume cells
       *
       *  The volume cells are numbered and stored individually.
       *
       *  \author hiermeier \date 12/16 */
      void GmshVolumeCellsOnly(const plain_volumecell_set& vcells);

      /** \brief Write a Gmsh output file, which contains only the given facets
       *
       *  The facets are numbered and stored individually.
       *
       *  \author hiermeier \date 12/16 */
      void GmshFacetsOnly(
          const plain_facet_set& facets, Element* ele, const std::string& file_affix = "");
      inline void GmshFacetsOnly(
          const plain_facet_set& facets, Element* ele, const int& file_affix = -100)
      {
        std::ostringstream affix;
        if (file_affix != -100) affix << file_affix;
        GmshFacetsOnly(facets, ele, affix.str());
      }

      /** \brief Write a Gmsh output file, which contains only the given volume cells
       *
       *  The volume cells are numbered and stored individually.
       *
       *  \author hiermeier \date 12/16 */
      void GmshEdgesOnly(const plain_edge_set& edges);

      /** \brief Write Gmsh output for a cell of given discretization type
       *
       *  \param file    (out) : file to write into
       *  \param shape    (in) : discretization type of the cell
       *  \param xyze     (in) : nodal coordinates
       *  \param position (in) : position of the cell ( optional )
       *  \param value    (in) : value of the cell ( optional ). If a value is given,
       *                         it will overrule the position input.
       *
       *  \author hiermeier \date 01/17 */
      void GmshCellDump(std::ofstream& file, DRT::Element::DiscretizationType shape,
          const Epetra_SerialDenseMatrix& xyze, const Point::PointPosition* position = NULL,
          const int* value = NULL);

      /*!
      \brief Write output of background element geometry
       */
      void GmshElementDump(std::ofstream& file, Element* ele, bool to_local = false);

      /*!
      \brief Write output of background element geometry
       */
      void GmshElementDump(std::ofstream& file, const std::vector<Node*>& nodes, char elementtype,
          bool to_local = false, Element* ele = NULL);

      /*!
      \brief Write output of a side
       */
      void GmshSideDump(
          std::ofstream& file, const Side* s, bool to_local = false, Element* ele = NULL);

      /*!
      \brief Write output of a side with given name
       */
      void GmshSideDump(std::ofstream& file, const Side* s, const std::string& sname,
          bool to_local = false, Element* ele = NULL);

      /*!
      \brief Write output of a Triside
       */
      void GmshTriSideDump(std::ofstream& file, std::vector<Point*> points, bool to_local = false,
          Element* ele = NULL);

      /*!
      \brief Write output of a facet
       */
      void GmshFacetDump(std::ofstream& file, Facet* facet,
          const std::string& visualizationtype = "sides", bool print_all = false,
          bool to_local = false, Element* ele = NULL);

      /*!
      \brief Write output of a volumecell
       */
      void GmshVolumecellDump(std::ofstream& file, VolumeCell* VC,
          const std::string& visualizationtype = "sides", bool print_all = false,
          bool to_local = false, Element* ele = NULL);

      /*!
      \brief Write output of a cylce
       */
      void GmshCycleDump(std::ofstream& file, Cycle* cycle,
          const std::string& visualizationtype = "sides", bool to_local = false,
          Element* ele = NULL);

      /*!
      \brief Write output of the background element and all the cut sides corresponding to this
      element
       */
      void GmshCompleteCutElement(std::ofstream& file, Element* ele, bool to_local = false);

      /*!
      \brief Write output of a coord as point with idx
       */
      void GmshCoordDump(std::ofstream& file, LINALG::Matrix<3, 1> coord, double idx,
          bool to_local = false, Element* ele = NULL);


      /*!
      \brief Write output of a line
       */
      void GmshLineDump(std::ofstream& file, GEO::CUT::Point* p1, GEO::CUT::Point* p2, int idx1,
          int idx2, bool to_local = false, Element* ele = NULL);

      inline void GmshLineDump(std::ofstream& file, GEO::CUT::Point* p1, GEO::CUT::Point* p2)
      {
        GmshLineDump(file, p1, p2, p1->Id(), p2->Id(), false, NULL);
      }

      inline void GmshLineDump(
          std::ofstream& file, GEO::CUT::Point* p1, GEO::CUT::Point* p2, bool to_local)
      {
        GmshLineDump(file, p1, p2, p1->Id(), p2->Id(), to_local, NULL);
      }

      inline void GmshLineDump(std::ofstream& file, GEO::CUT::Point* p1, GEO::CUT::Point* p2,
          bool to_local, Element* ele)
      {
        GmshLineDump(file, p1, p2, p1->Id(), p2->Id(), to_local, ele);
      }

      void GmshLineDump(
          std::ofstream& file, GEO::CUT::Line* line, bool to_local = false, Element* ele = NULL);
      /*!
      \brief Write output of a edge
       */
      void GmshEdgeDump(
          std::ofstream& file, GEO::CUT::Edge* edge, bool to_local = false, Element* ele = NULL);

      /*!
      \brief Write output of a edge
       */
      void GmshEdgeDump(std::ofstream& file, GEO::CUT::Edge* edge, const std::string& ename,
          bool to_local = false, Element* ele = NULL);

      /*!
      \brief Write output of a node
       */
      void GmshNodeDump(
          std::ofstream& file, GEO::CUT::Node* node, bool to_local = false, Element* ele = NULL);

      /*!
      \brief Write output of a point with special idx
       */
      void GmshPointDump(std::ofstream& file, GEO::CUT::Point* point, int idx,
          bool to_local = false, Element* ele = NULL);

      /*!
      \brief Write output of a point with special idx and special name
       */
      void GmshPointDump(std::ofstream& file, GEO::CUT::Point* point, int idx,
          const std::string& pname, bool to_local, Element* ele);

      /*!
      \brief Write output of a point with point position as idx
       */
      void GmshPointDump(
          std::ofstream& file, GEO::CUT::Point* point, bool to_local = false, Element* ele = NULL);

      /*!
      \brief Write level set value on cut surface (information taken from facet).
       */
      void GmshLevelSetValueDump(
          std::ofstream& file, Element* ele, bool dumpnodevalues = false, bool to_local = false);

      /*!
      \brief Write level set gradient on cut surface (information taken from facet).
       */
      void GmshLevelSetGradientDump(std::ofstream& file, Element* ele, bool to_local = false);

      /*!
      \brief Write level set value on cut surface (information taken from facet).
       */
      void GmshLevelSetValueZeroSurfaceDump(
          std::ofstream& file, Element* ele, bool to_local = false);

      /*!
       * Write Level Set Gradient Orientation of Boundary-Cell Normal and LevelSet
       */
      void GmshLevelSetOrientationDump(std::ofstream& file, Element* ele, bool to_local = false);

      /*!
      \brief Write Eqn of plane normal for all facets (used for DirectDivergence).
      */
      void GmshEqnPlaneNormalDump(
          std::ofstream& file, Element* ele, bool normalize = false, bool to_local = false);
      void GmshEqnPlaneNormalDump(std::ofstream& file, Facet* facet, bool normalize = false,
          bool to_local = false, Element* ele = NULL);
      /*!
      \brief Get equation of plane as implemented in DirectDivergence routine.
       */
      std::vector<double> GetEqOfPlane(std::vector<Point*> pts);

      /*!
      \brief Simplify output of for normal output options.
       */
      void GmshScalar(std::ofstream& file, LINALG::Matrix<3, 1> coord, double scalar,
          bool to_local = false, Element* ele = NULL);
      void GmshVector(std::ofstream& file, LINALG::Matrix<3, 1> coord, std::vector<double> vector,
          bool normalize = false, bool to_local = false, Element* ele = NULL);

      /*!
      \brief Write cuttest for this element!
      */
      void GmshElementCutTest(
          std::ofstream& file, GEO::CUT::Element* ele, bool haslevelsetside = false);

      /*!
       \brief Generate filename for gmsh output with specific ending
       */
      std::string GenerateGmshOutputFilename(const std::string& filename_tail);

      /*!
       \brief Write new Section in Gmsh file (eventually end section from before...)
       */
      void GmshNewSection(
          std::ofstream& file, const std::string& section, bool first_endsection = false);

      /*!
       \brief End Section in Gmsh file
       */
      void GmshEndSection(std::ofstream& file, bool close_file = false);

      /// generate combination of output of this particular edge and intersection, useful for
      /// debugging cut libraries
      void GmshCutPairDump(
          std::ofstream& file, Side* side, Edge* edge, int id, const std::string& suffix);

      void GmshCutPairDump(std::ofstream& file, const std::pair<Side*, Edge*>& pair, int id,
          const std::string& suffix);

      // void GmshObjectDump2(std::ofstream & file, GEO::CUT::Facet* f)
      //{
      //  return;
      //};//GmshFacetDump(file,f);}

      /*!
       \brief Writes the whole container into a gmsh file!
       */
      // template<class T>
      // void GmshPrintContainer(std::ofstream & file, const std::string & section, std::vector<T>
      // container)
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
      void GmshWriteCoords(std::ofstream& file, std::vector<double> coord, bool to_local = false,
          Element* ele = NULL);

      void GmshWriteCoords(std::ofstream& file, LINALG::Matrix<3, 1> coord, bool to_local = false,
          Element* ele = NULL);

      void GmshWriteCoords(
          std::ofstream& file, Node* node, bool to_local = false, Element* ele = NULL);

      void GmshWriteCoords(
          std::ofstream& file, Point* point, bool to_local = false, Element* ele = NULL);

      // Generate debug output for various critical cases
      void DebugDump_ThreePointsOnEdge(
          Side* first, Side* second, Edge* e, Point* p, const PointSet& cut);

      void DebugDump_MoreThanTwoIntersectionPoints(
          Edge* edge, Side* other, const std::vector<Point*>& point_stack);

      void DebugDump_MultipleCutPointsSpecial(Side* first, Side* second, const PointSet& cut,
          const PointSet& collected_points, const point_line_set& new_lines);

      /*!
      \brief Dumps object into gmsh file!
      */
      template <class T>
      void GmshObjectDump(std::ofstream& file, T* obj, bool to_local = false, Element* ele = NULL)
      {
        dserror("GmshObjectDump: no specific implementation for your Object type!");
      }

      template <>
      inline void GmshObjectDump<Point>(
          std::ofstream& file, Point* obj, bool to_local, Element* ele)
      {
        GmshPointDump(file, obj, to_local, ele);
      }

      template <>
      inline void GmshObjectDump<Node>(std::ofstream& file, Node* obj, bool to_local, Element* ele)
      {
        GmshNodeDump(file, obj, to_local, ele);
      }

      template <>
      inline void GmshObjectDump<Element>(
          std::ofstream& file, Element* obj, bool to_local, Element* ele)
      {
        GmshElementDump(file, obj, to_local);
      }

      template <>
      inline void GmshObjectDump<Edge>(std::ofstream& file, Edge* obj, bool to_local, Element* ele)
      {
        GmshEdgeDump(file, obj, to_local, ele);
      }

      template <>
      inline void GmshObjectDump<Side>(std::ofstream& file, Side* obj, bool to_local, Element* ele)
      {
        GmshSideDump(file, obj, to_local, ele);
      }

      template <>
      inline void GmshObjectDump<Line>(std::ofstream& file, Line* obj, bool to_local, Element* ele)
      {
        GmshLineDump(file, obj, to_local, ele);
      }

      template <>
      inline void GmshObjectDump<Facet>(
          std::ofstream& file, Facet* obj, bool to_local, Element* ele)
      {
        GmshFacetDump(file, obj, "sides", true, to_local, ele);
      }

      template <>
      inline void GmshObjectDump<VolumeCell>(
          std::ofstream& file, VolumeCell* obj, bool to_local, Element* ele)
      {
        GmshVolumecellDump(file, obj, "sides", true, to_local, ele);
      }

      template <>
      inline void GmshObjectDump<Cycle>(
          std::ofstream& file, Cycle* obj, bool to_local, Element* ele)
      {
        GmshCycleDump(file, obj, "lines", to_local, ele);
      }

      template <>
      inline void GmshObjectDump<BoundaryCell>(
          std::ofstream& file, BoundaryCell* obj, bool to_local, Element* ele)
      {
        GmshCellDump(file, obj->Shape(), obj->Coordinates());
      }

      /*!
       \brief Writes the geom. object  into a section gmsh file!
       */
      // for std::vector<T*>
      template <class T>
      void GmshWriteSection(std::ofstream& file, const std::string& section,
          std::vector<T*> container, bool close_file = false, bool to_local = false,
          Element* ele = NULL)
      {
        if (section != "") GmshNewSection(file, section);
        for (typename std::vector<T*>::iterator t = container.begin(); t != container.end(); ++t)
          GmshObjectDump<T>(file, (*t), to_local, ele);
        if (section != "") GmshEndSection(file, close_file);
      }

      // for std::vector<Teuchos::RCP<T> >
      template <class T>
      void GmshWriteSection(std::ofstream& file, const std::string& section,
          std::vector<Teuchos::RCP<T>> container, bool close_file = false, bool to_local = false,
          Element* ele = NULL)
      {
        if (section != "") GmshNewSection(file, section);
        for (typename std::vector<Teuchos::RCP<T>>::iterator t = container.begin();
             t != container.end(); ++t)
          GmshObjectDump<T>(file, &(*(*t)), to_local, ele);
        if (section != "") GmshEndSection(file, close_file);
      }

      // for std::map<int, Teuchos::RCP<T> >
      template <class T>
      void GmshWriteSection(std::ofstream& file, const std::string& section,
          std::map<int, Teuchos::RCP<T>> container, bool close_file = false, bool to_local = false,
          Element* ele = NULL)
      {
        if (section != "") GmshNewSection(file, section);
        for (typename std::map<int, Teuchos::RCP<T>>::iterator t = container.begin();
             t != container.end(); ++t)
          GmshObjectDump<T>(file, &(*(t->second)), to_local, ele);
        if (section != "") GmshEndSection(file, close_file);
      }

      // for std::map<plain_int_set, Teuchos::RCP<T> >
      template <class T>
      void GmshWriteSection(std::ofstream& file, const std::string& section,
          std::map<plain_int_set, Teuchos::RCP<T>> container, bool close_file = false,
          bool to_local = false, Element* ele = NULL)
      {
        if (section != "") GmshNewSection(file, section);
        for (typename std::map<plain_int_set, Teuchos::RCP<T>>::iterator t = container.begin();
             t != container.end(); ++t)
        {
          GmshObjectDump<T>(file, &(*(t->second)), to_local, ele);
        }
        if (section != "") GmshEndSection(file, close_file);
      }

      // for sorted_vector<T*>
      template <class T>
      void GmshWriteSection(std::ofstream& file, const std::string& section,
          sorted_vector<T*> container, bool close_file = false, bool to_local = false,
          Element* ele = NULL)
      {
        if (section != "") GmshNewSection(file, section);
        for (typename sorted_vector<T*>::iterator t = container.begin(); t != container.end(); ++t)
          GmshObjectDump<T>(file, (*t), to_local, ele);
        if (section != "") GmshEndSection(file, close_file);
      }

      // for sorted_vector<T*,true,PointPidLess>
      template <class T>
      void GmshWriteSection(std::ofstream& file, const std::string& section,
          sorted_vector<T*, true, PointPidLess> container, bool close_file = false,
          bool to_local = false, Element* ele = NULL)
      {
        if (section != "") GmshNewSection(file, section);
        for (typename sorted_vector<T*, true, PointPidLess>::iterator t = container.begin();
             t != container.end(); ++t)
          GmshObjectDump<T>(file, (*t), to_local, ele);
        if (section != "") GmshEndSection(file, close_file);
      }

      // for std::set<T*>
      template <class T>
      void GmshWriteSection(std::ofstream& file, const std::string& section, std::set<T*> container,
          bool close_file = false, bool to_local = false, Element* ele = NULL)
      {
        if (section != "") GmshNewSection(file, section);
        for (typename std::set<T*>::iterator t = container.begin(); t != container.end(); ++t)
          GmshObjectDump<T>(file, (*t), to_local, ele);
        if (section != "") GmshEndSection(file, close_file);
      }

      // for T*
      template <class T>
      void GmshWriteSection(std::ofstream& file, const std::string& section, T* obj,
          bool close_file = false, bool to_local = false, Element* ele = NULL)
      {
        if (section != "") GmshNewSection(file, section);
        GmshObjectDump<T>(file, obj, to_local, ele);
        if (section != "") GmshEndSection(file, close_file);
      }

    }  // namespace OUTPUT

  } /* namespace CUT */
} /* namespace GEO */
#endif /* CUTOUTPUT_H_ */
