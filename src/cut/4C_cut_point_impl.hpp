/*---------------------------------------------------------------------*/
/*! \file

\brief tbd

\level 3


*----------------------------------------------------------------------*/

#ifndef FOUR_C_CUT_POINT_IMPL_HPP
#define FOUR_C_CUT_POINT_IMPL_HPP

#include "4C_config.hpp"

#include "4C_cut_line.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Core::Geo
{
  namespace Cut
  {
    namespace Impl
    {
      class LineBetweenFilter
      {
       public:
        LineBetweenFilter(Point* me, Point* other) : me_(me), other_(other) {}

        bool operator()(Line* line) { return line->Between(me_, other_); }

       private:
        Point* me_;
        Point* other_;
      };

      class LineHasSideFilter
      {
       public:
        explicit LineHasSideFilter(Side* side) : side_(side) {}

        /// true if the line is cut by the side but not on any side's edges
        bool operator()(Line* line) { return line->IsInternalCut(side_); }

       private:
        Side* side_;
      };

      class NextLineOnElementCutFilter
      {
       public:
        NextLineOnElementCutFilter(Line* line, Side* side, Element* element)
            : line_(line), side_(side), element_(element)
        {
        }

        bool operator()(Line* line)
        {
          return line != line_ and line->IsCut(side_) and
                 (element_ == nullptr or line->IsCut(element_));
        }

       private:
        Line* line_;
        Side* side_;
        Element* element_;
      };

    }  // namespace Impl

    template <class Filter>
    Line* Point::find(Filter& filter, bool unique)
    {
      Line* line_found = nullptr;
      for (plain_line_set::iterator i = lines_.begin(); i != lines_.end(); ++i)
      {
        Line* line = *i;
        if (filter(line))
        {
          if (line_found == nullptr)
          {
            line_found = line;
            if (not unique)
            {
              break;
            }
          }
          else
          {
            FOUR_C_THROW("not unique");
          }
        }
      }
      return line_found;
    }

    inline Line* Point::CommonLine(Point* other)
    {
      Impl::LineBetweenFilter filter(this, other);
      return find(filter, true);
    }

    inline Line* Point::CutLine(Side* side, bool unique)
    {
      Impl::LineHasSideFilter filter(side);
      return find(filter, unique);
    }

    inline Line* Point::CutLine(Line* line, Side* side, Element* element)
    {
      Impl::NextLineOnElementCutFilter filter(line, side, element);
      return find(filter, true);
    }

  }  // namespace Cut
}  // namespace Core::Geo

FOUR_C_NAMESPACE_CLOSE

#endif
