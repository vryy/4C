/*---------------------------------------------------------------------*/
/*! \file

\brief cut line

\level 3


*----------------------------------------------------------------------*/

#ifndef FOUR_C_CUT_LINE_HPP
#define FOUR_C_CUT_LINE_HPP

#include "4C_config.hpp"

#include "4C_cut_point.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Core::Geo
{
  namespace Cut
  {
    /*!
    \begin Line between two points. These lines result from cuts and there are no cut points on a
    line.
     */
    class Line
    {
     public:
      Line(Point* p1, Point* p2, Side* cut_side1, Side* cut_side2, Element* cut_element);

      void add_side(Side* cut_side);

      void add_element(Element* cut_element);

      bool is_cut(Side* s1, Side* s2)
      {
        return cut_sides_.count(s1) > 0 and cut_sides_.count(s2) > 0;
      }

      bool is_cut(Element* element) { return cut_elements_.count(element) > 0; }

      bool is_cut(Side* side)
      {
        return cut_sides_.count(side) > 0;
        //     return ( cut_sides_.count( side ) > 0 and
        //              BeginPoint()->IsCut( side ) and
        //              EndPoint()->IsCut( side ) );
      }

      bool is_internal_cut(Side* side);

      bool on_edge(Edge* edge) { return p1_->is_cut(edge) and p2_->is_cut(edge); }

      Point* other_point(Point* point)
      {
        if (p1_ == point)
          return p2_;
        else if (p2_ == point)
          return p1_;
        else
          FOUR_C_THROW("foreign point provided");
      }

      Point* begin_point() { return p1_; }
      const Point* begin_point() const { return p1_; }

      Point* end_point() { return p2_; }
      const Point* end_point() const { return p2_; }

      bool between(Point* p1, Point* p2)
      {
        return ((p1 == p1_ and p2 == p2_) or (p1 == p2_ and p2 == p1_));
      }

      /*! \brief Print the coordinates of the points on screen */
      void print()
      {
        p1_->print();
        std::cout << "--";
        p2_->print();
        std::cout << "\n";
      }

      void plot(std::ofstream& f) const
      {
        f << "# line\n";
        p1_->plot(f);
        p2_->plot(f);
        f << "\n\n";
      }

      void intersection(plain_side_set& sides)
      {
        plain_side_set intersection;
        std::set_intersection(cut_sides_.begin(), cut_sides_.end(), sides.begin(), sides.end(),
            std::inserter(intersection, intersection.begin()));
        std::swap(sides, intersection);
      }

      const plain_side_set& cut_sides() const { return cut_sides_; }

      /// Replace point p by point p_new in the line
      void replace(Point* p, Point* p_new)
      {
        auto& replace_p = (p == p1_) ? p1_ : p2_;
        replace_p = p_new;
      }

     private:
      Point* p1_;
      Point* p2_;

      plain_side_set cut_sides_;

      plain_element_set cut_elements_;
    };

  }  // namespace Cut
}  // namespace Core::Geo

FOUR_C_NAMESPACE_CLOSE

#endif
