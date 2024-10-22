// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_UTILS_STD_CXX20_RANGES_HPP
#define FOUR_C_UTILS_STD_CXX20_RANGES_HPP

#include "4C_config.hpp"

#include <iterator>

FOUR_C_NAMESPACE_OPEN

namespace std_20  // NOLINT
{
  namespace ranges::views  // NOLINT
  {
    namespace views = ranges::views;
    namespace Internal
    {
      template <typename Iterator>
      class IteratorRange
      {
       public:
        using iterator = Iterator;

        IteratorRange(Iterator begin, Iterator end) : begin_(begin), end_(end) {}

        [[nodiscard]] Iterator begin() const { return begin_; }  // NOLINT
        [[nodiscard]] Iterator end() const { return end_; }      // NOLINT

       private:
        Iterator begin_;
        Iterator end_;
      };

      template <typename Range, typename Predicate>
      class FilterIterator
      {
       public:
        using iterator = typename Range::iterator;

        FilterIterator(iterator current, iterator end, Predicate pred)
            : current_(current), end_(end), pred_(pred)
        {
          advance_to_next_valid();
        }

        FilterIterator& operator++()
        {
          ++current_;
          advance_to_next_valid();
          return *this;
        }

        auto operator*() const { return *current_; }
        auto operator->() const { return current_; }

        bool operator==(const FilterIterator& other) const { return current_ == other.current_; }

        bool operator!=(const FilterIterator& other) const { return current_ != other.current_; }

       private:
        void advance_to_next_valid()
        {
          while (current_ != end_ && !pred_(*current_))
          {
            ++current_;
          }
        }

        iterator current_;
        iterator end_;
        Predicate pred_;
      };

      template <typename Range, typename Predicate>
      class FilterRange
      {
       public:
        using iterator = FilterIterator<Range, Predicate>;

        FilterRange(Range range, Predicate pred) : range_(range), pred_(pred) {}

        [[nodiscard]] iterator begin() const
        {
          return iterator(range_.begin(), range_.end(), pred_);
        }

        [[nodiscard]] iterator end() const { return iterator(range_.end(), range_.end(), pred_); }

       private:
        Range range_;
        Predicate pred_;
      };
    }  // namespace Internal

    /**
     * \brief Returns a view that includes all elements of the given @p container.
     *
     */
    template <typename Container>
    auto all(Container& container)
    {
      return Internal::IteratorRange(std::begin(container), std::end(container));
    }

    /**
     * \brief Returns a view that includes only elements of @p range satisfying the given @p
     * predicate.
     */
    template <typename Range, typename Predicate>
    auto filter(Range range, Predicate predicate)
    {
      return Internal::FilterRange<Range, Predicate>(range, predicate);
    }
  }  // namespace ranges::views
}  // namespace std_20

FOUR_C_NAMESPACE_CLOSE

#endif
