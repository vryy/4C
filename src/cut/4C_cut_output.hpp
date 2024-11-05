// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_CUT_OUTPUT_HPP
#define FOUR_C_CUT_OUTPUT_HPP

#include "4C_config.hpp"

#include "4C_cut_boundarycell.hpp"
#include "4C_cut_point.hpp"

FOUR_C_NAMESPACE_OPEN


namespace Cut
{
  class Element;
  class Side;
  class Node;
  class Point;
  class Line;
  class Edge;
  class Cycle;

  namespace Output
  {
    char gmsh_element_type(Core::FE::CellType shape);

    /** \brief Write a Gmsh output file, which contains only the given volume cells
     *
     *  The volume cells are numbered and stored individually.
     *
     *  \author hiermeier \date 12/16 */
    void gmsh_volume_cells_only(const plain_volumecell_set& vcells);

    /** \brief Write a Gmsh output file, which contains only the given facets
     *
     *  The facets are numbered and stored individually.
     *
     *  \author hiermeier \date 12/16 */
    void gmsh_facets_only(
        const plain_facet_set& facets, Element* ele, const std::string& file_affix = "");
    inline void gmsh_facets_only(
        const plain_facet_set& facets, Element* ele, const int& file_affix = -100)
    {
      std::ostringstream affix;
      if (file_affix != -100) affix << file_affix;
      gmsh_facets_only(facets, ele, affix.str());
    }

    /** \brief Write a Gmsh output file, which contains only the given volume cells
     *
     *  The volume cells are numbered and stored individually.
     *
     *  \author hiermeier \date 12/16 */
    void gmsh_edges_only(const plain_edge_set& edges);

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
    void gmsh_cell_dump(std::ofstream& file, Core::FE::CellType shape,
        const Core::LinAlg::SerialDenseMatrix& xyze, const Point::PointPosition* position = nullptr,
        const int* value = nullptr);

    /*!
    \brief Write output of background element geometry
     */
    void gmsh_element_dump(std::ofstream& file, Element* ele, bool to_local = false);

    /*!
    \brief Write output of background element geometry
     */
    void gmsh_element_dump(std::ofstream& file, const std::vector<Node*>& nodes, char elementtype,
        bool to_local = false, Element* ele = nullptr);

    /*!
    \brief Write output of a side
     */
    void gmsh_side_dump(
        std::ofstream& file, const Side* s, bool to_local = false, Element* ele = nullptr);

    /*!
    \brief Write output of a side with given name
     */
    void gmsh_side_dump(std::ofstream& file, const Side* s, const std::string& sname,
        bool to_local = false, Element* ele = nullptr);

    /*!
    \brief Write output of a Triside
     */
    void gmsh_tri_side_dump(std::ofstream& file, std::vector<Point*> points, bool to_local = false,
        Element* ele = nullptr);

    /*!
    \brief Write output of a facet
     */
    void gmsh_facet_dump(std::ofstream& file, Facet* facet,
        const std::string& visualizationtype = "sides", bool print_all = false,
        bool to_local = false, Element* ele = nullptr);

    /*!
    \brief Write output of a volumecell
     */
    void gmsh_volumecell_dump(std::ofstream& file, VolumeCell* VC,
        const std::string& visualizationtype = "sides", bool print_all = false,
        bool to_local = false, Element* ele = nullptr);

    /*!
    \brief Write output of a cylce
     */
    void gmsh_cycle_dump(std::ofstream& file, Cycle* cycle,
        const std::string& visualizationtype = "sides", bool to_local = false,
        Element* ele = nullptr);

    /*!
    \brief Write output of the background element and all the cut sides corresponding to this
    element
     */
    void gmsh_complete_cut_element(std::ofstream& file, Element* ele, bool to_local = false);

    /*!
    \brief Write output of a coord as point with idx
     */
    void gmsh_coord_dump(std::ofstream& file, Core::LinAlg::Matrix<3, 1> coord, double idx,
        bool to_local = false, Element* ele = nullptr);


    /*!
    \brief Write output of a line
     */
    void gmsh_line_dump(std::ofstream& file, Cut::Point* p1, Cut::Point* p2, int idx1, int idx2,
        bool to_local = false, Element* ele = nullptr);

    inline void gmsh_line_dump(std::ofstream& file, Cut::Point* p1, Cut::Point* p2)
    {
      gmsh_line_dump(file, p1, p2, p1->id(), p2->id(), false, nullptr);
    }

    inline void gmsh_line_dump(std::ofstream& file, Cut::Point* p1, Cut::Point* p2, bool to_local)
    {
      gmsh_line_dump(file, p1, p2, p1->id(), p2->id(), to_local, nullptr);
    }

    inline void gmsh_line_dump(
        std::ofstream& file, Cut::Point* p1, Cut::Point* p2, bool to_local, Element* ele)
    {
      gmsh_line_dump(file, p1, p2, p1->id(), p2->id(), to_local, ele);
    }

    void gmsh_line_dump(
        std::ofstream& file, Cut::Line* line, bool to_local = false, Element* ele = nullptr);
    /*!
    \brief Write output of a edge
     */
    void gmsh_edge_dump(
        std::ofstream& file, Cut::Edge* edge, bool to_local = false, Element* ele = nullptr);

    /*!
    \brief Write output of a edge
     */
    void gmsh_edge_dump(std::ofstream& file, Cut::Edge* edge, const std::string& ename,
        bool to_local = false, Element* ele = nullptr);

    /*!
    \brief Write output of a node
     */
    void gmsh_node_dump(
        std::ofstream& file, Cut::Node* node, bool to_local = false, Element* ele = nullptr);

    /*!
    \brief Write output of a point with special idx
     */
    void gmsh_point_dump(std::ofstream& file, Cut::Point* point, int idx, bool to_local = false,
        Element* ele = nullptr);

    /*!
    \brief Write output of a point with special idx and special name
     */
    void gmsh_point_dump(std::ofstream& file, Cut::Point* point, int idx, const std::string& pname,
        bool to_local, Element* ele);

    /*!
    \brief Write output of a point with point position as idx
     */
    void gmsh_point_dump(
        std::ofstream& file, Cut::Point* point, bool to_local = false, Element* ele = nullptr);

    /*!
    \brief Write level set value on cut surface (information taken from facet).
     */
    void gmsh_level_set_value_dump(
        std::ofstream& file, Element* ele, bool dumpnodevalues = false, bool to_local = false);

    /*!
    \brief Write level set gradient on cut surface (information taken from facet).
     */
    void gmsh_level_set_gradient_dump(std::ofstream& file, Element* ele, bool to_local = false);

    /*!
    \brief Write level set value on cut surface (information taken from facet).
     */
    void gmsh_level_set_value_zero_surface_dump(
        std::ofstream& file, Element* ele, bool to_local = false);

    /*!
     * Write Level Set Gradient Orientation of Boundary-Cell Normal and LevelSet
     */
    void gmsh_level_set_orientation_dump(std::ofstream& file, Element* ele, bool to_local = false);

    /*!
    \brief Write Eqn of plane normal for all facets (used for DirectDivergence).
    */
    void gmsh_eqn_plane_normal_dump(
        std::ofstream& file, Element* ele, bool normalize = false, bool to_local = false);
    void gmsh_eqn_plane_normal_dump(std::ofstream& file, Facet* facet, bool normalize = false,
        bool to_local = false, Element* ele = nullptr);
    /*!
    \brief Get equation of plane as implemented in DirectDivergence routine.
     */
    std::vector<double> get_eq_of_plane(std::vector<Point*> pts);

    /*!
    \brief Simplify output of for normal output options.
     */
    void gmsh_scalar(std::ofstream& file, Core::LinAlg::Matrix<3, 1> coord, double scalar,
        bool to_local = false, Element* ele = nullptr);
    void gmsh_vector(std::ofstream& file, Core::LinAlg::Matrix<3, 1> coord,
        std::vector<double> vector, bool normalize = false, bool to_local = false,
        Element* ele = nullptr);

    /*!
    \brief Write cuttest for this element!
    */
    void gmsh_element_cut_test(
        std::ofstream& file, Cut::Element* ele, bool haslevelsetside = false);

    /*!
     \brief Generate filename for gmsh output with specific ending
     */
    std::string generate_gmsh_output_filename(const std::string& filename_tail);

    /*!
     \brief Write new Section in Gmsh file (eventually end section from before...)
     */
    void gmsh_new_section(
        std::ofstream& file, const std::string& section, bool first_endsection = false);

    /*!
     \brief End Section in Gmsh file
     */
    void gmsh_end_section(std::ofstream& file, bool close_file = false);

    /// generate combination of output of this particular edge and intersection, useful for
    /// debugging cut libraries
    void gmsh_cut_pair_dump(
        std::ofstream& file, Side* side, Edge* edge, int id, const std::string& suffix);

    void gmsh_cut_pair_dump(std::ofstream& file, const std::pair<Side*, Edge*>& pair, int id,
        const std::string& suffix);

    /*!
     \brief Write Coordinates in Gmsh file (for internal use)
     //to_local ... transform to local coordinates of the ele?
     */
    void gmsh_write_coords(std::ofstream& file, std::vector<double> coord, bool to_local = false,
        Element* ele = nullptr);

    void gmsh_write_coords(std::ofstream& file, Core::LinAlg::Matrix<3, 1> coord,
        bool to_local = false, Element* ele = nullptr);

    void gmsh_write_coords(
        std::ofstream& file, Node* node, bool to_local = false, Element* ele = nullptr);

    void gmsh_write_coords(
        std::ofstream& file, Point* point, bool to_local = false, Element* ele = nullptr);

    // Generate debug output for various critical cases
    void debug_dump_three_points_on_edge(
        Side* first, Side* second, Edge* e, Point* p, const PointSet& cut);

    void debug_dump_more_than_two_intersection_points(
        Edge* edge, Side* other, const std::vector<Point*>& point_stack);

    void debug_dump_multiple_cut_points_special(Side* first, Side* second, const PointSet& cut,
        const PointSet& collected_points, const point_line_set& new_lines);

    /*!
    \brief Dumps object into gmsh file!
    */
    template <class T>
    void gmsh_object_dump(
        std::ofstream& file, T* obj, bool to_local = false, Element* ele = nullptr)
    {
      FOUR_C_THROW("gmsh_object_dump: no specific implementation for your Object type!");
    }

    template <>
    inline void gmsh_object_dump<Point>(
        std::ofstream& file, Point* obj, bool to_local, Element* ele)
    {
      gmsh_point_dump(file, obj, to_local, ele);
    }

    template <>
    inline void gmsh_object_dump<Node>(std::ofstream& file, Node* obj, bool to_local, Element* ele)
    {
      gmsh_node_dump(file, obj, to_local, ele);
    }

    template <>
    inline void gmsh_object_dump<Element>(
        std::ofstream& file, Element* obj, bool to_local, Element* ele)
    {
      gmsh_element_dump(file, obj, to_local);
    }

    template <>
    inline void gmsh_object_dump<Edge>(std::ofstream& file, Edge* obj, bool to_local, Element* ele)
    {
      gmsh_edge_dump(file, obj, to_local, ele);
    }

    template <>
    inline void gmsh_object_dump<Side>(std::ofstream& file, Side* obj, bool to_local, Element* ele)
    {
      gmsh_side_dump(file, obj, to_local, ele);
    }

    template <>
    inline void gmsh_object_dump<Line>(std::ofstream& file, Line* obj, bool to_local, Element* ele)
    {
      gmsh_line_dump(file, obj, to_local, ele);
    }

    template <>
    inline void gmsh_object_dump<Facet>(
        std::ofstream& file, Facet* obj, bool to_local, Element* ele)
    {
      gmsh_facet_dump(file, obj, "sides", true, to_local, ele);
    }

    template <>
    inline void gmsh_object_dump<VolumeCell>(
        std::ofstream& file, VolumeCell* obj, bool to_local, Element* ele)
    {
      gmsh_volumecell_dump(file, obj, "sides", true, to_local, ele);
    }

    template <>
    inline void gmsh_object_dump<Cycle>(
        std::ofstream& file, Cycle* obj, bool to_local, Element* ele)
    {
      gmsh_cycle_dump(file, obj, "lines", to_local, ele);
    }

    template <>
    inline void gmsh_object_dump<BoundaryCell>(
        std::ofstream& file, BoundaryCell* obj, bool to_local, Element* ele)
    {
      gmsh_cell_dump(file, obj->shape(), obj->coordinates());
    }

    /*!
     \brief Writes the geom. object  into a section gmsh file!
     */
    // for std::vector<T*>
    template <class T>
    void gmsh_write_section(std::ofstream& file, const std::string& section,
        std::vector<T*> container, bool close_file = false, bool to_local = false,
        Element* ele = nullptr)
    {
      if (section != "") gmsh_new_section(file, section);
      for (typename std::vector<T*>::iterator t = container.begin(); t != container.end(); ++t)
        gmsh_object_dump<T>(file, (*t), to_local, ele);
      if (section != "") gmsh_end_section(file, close_file);
    }

    // for std::vector<std::shared_ptr<T> >
    template <class T>
    void gmsh_write_section(std::ofstream& file, const std::string& section,
        std::vector<std::shared_ptr<T>> container, bool close_file = false, bool to_local = false,
        Element* ele = nullptr)
    {
      if (section != "") gmsh_new_section(file, section);
      for (typename std::vector<std::shared_ptr<T>>::iterator t = container.begin();
           t != container.end(); ++t)
        gmsh_object_dump<T>(file, &(*(*t)), to_local, ele);
      if (section != "") gmsh_end_section(file, close_file);
    }

    // for std::map<int, std::shared_ptr<T> >
    template <class T>
    void gmsh_write_section(std::ofstream& file, const std::string& section,
        std::map<int, std::shared_ptr<T>> container, bool close_file = false, bool to_local = false,
        Element* ele = nullptr)
    {
      if (section != "") gmsh_new_section(file, section);
      for (typename std::map<int, std::shared_ptr<T>>::iterator t = container.begin();
           t != container.end(); ++t)
        gmsh_object_dump<T>(file, &(*(t->second)), to_local, ele);
      if (section != "") gmsh_end_section(file, close_file);
    }

    // for std::map<plain_int_set, std::shared_ptr<T> >
    template <class T>
    void gmsh_write_section(std::ofstream& file, const std::string& section,
        std::map<plain_int_set, std::shared_ptr<T>> container, bool close_file = false,
        bool to_local = false, Element* ele = nullptr)
    {
      if (section != "") gmsh_new_section(file, section);
      for (typename std::map<plain_int_set, std::shared_ptr<T>>::iterator t = container.begin();
           t != container.end(); ++t)
      {
        gmsh_object_dump<T>(file, &(*(t->second)), to_local, ele);
      }
      if (section != "") gmsh_end_section(file, close_file);
    }

    // for sorted_vector<T*>
    template <class T>
    void gmsh_write_section(std::ofstream& file, const std::string& section,
        SortedVector<T*> container, bool close_file = false, bool to_local = false,
        Element* ele = nullptr)
    {
      if (section != "") gmsh_new_section(file, section);
      for (typename SortedVector<T*>::iterator t = container.begin(); t != container.end(); ++t)
        gmsh_object_dump<T>(file, (*t), to_local, ele);
      if (section != "") gmsh_end_section(file, close_file);
    }

    // for sorted_vector<T*,true,PointPidLess>
    template <class T>
    void gmsh_write_section(std::ofstream& file, const std::string& section,
        SortedVector<T*, true, PointPidLess> container, bool close_file = false,
        bool to_local = false, Element* ele = nullptr)
    {
      if (section != "") gmsh_new_section(file, section);
      for (typename SortedVector<T*, true, PointPidLess>::iterator t = container.begin();
           t != container.end(); ++t)
        gmsh_object_dump<T>(file, (*t), to_local, ele);
      if (section != "") gmsh_end_section(file, close_file);
    }

    // for std::set<T*>
    template <class T>
    void gmsh_write_section(std::ofstream& file, const std::string& section, std::set<T*> container,
        bool close_file = false, bool to_local = false, Element* ele = nullptr)
    {
      if (section != "") gmsh_new_section(file, section);
      for (typename std::set<T*>::iterator t = container.begin(); t != container.end(); ++t)
        gmsh_object_dump<T>(file, (*t), to_local, ele);
      if (section != "") gmsh_end_section(file, close_file);
    }

    // for T*
    template <class T>
    void gmsh_write_section(std::ofstream& file, const std::string& section, T* obj,
        bool close_file = false, bool to_local = false, Element* ele = nullptr)
    {
      if (section != "") gmsh_new_section(file, section);
      gmsh_object_dump<T>(file, obj, to_local, ele);
      if (section != "") gmsh_end_section(file, close_file);
    }

  }  // namespace Output

} /* namespace Cut */


FOUR_C_NAMESPACE_CLOSE

#endif
