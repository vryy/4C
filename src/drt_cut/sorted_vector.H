/*---------------------------------------------------------------------*/
/*!
\file sorted_vector.H

\brief  sorted vector that emulates a set (and is supposed to be much more efficient)

\level 3

<pre>
\maintainer Christoph Ager
            ager@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15249
</pre>

*----------------------------------------------------------------------*/

/* STL-conforming "sorted vector" container
 *
 * (C) 2002 Martin Holzherr (holzherr@infobrain.com). All rights reserved.
 *
 * Permission is granted to use, distribute and modify this code provided that:
 *   - this copyright notice appears,
 *
 * The author welcomes any suggestions on the code or reportings of actual
 * use of the code. Please send your comments to holzherr@infobrain.com.
 *
 * The author makes NO WARRANTY or representation, either express or implied,
 * with respect to this code, its quality, accuracy, merchantability, or
 * fitness for a particular purpose.  This software is provided "AS IS", and
 * you, its user, assume the entire risk as to its quality and accuracy.
 *
 * Created:         November 19th, 2002
 * Last modified:   November 27th, 2002
                        (changed namespace from std to codeproject;
                        uses template member functions for MSCVER>=1300)

 */

#ifndef SORTED_VECTOR_H
#define SORTED_VECTOR_H

#include <algorithm>
#include <vector>
#include <utility>
#include <functional>
#include <stdexcept>

#include <boost/operators.hpp>

namespace GEO
{
  namespace CUT
  {
    /// sorted vector that emulates a set (and is supposed to be much more efficient)
    template <class K, bool bNoDuplicates = true, class Pr = std::less<K>,
        class A = std::allocator<K>>
    class sorted_vector
        : public boost::less_than_comparable<sorted_vector<K, bNoDuplicates, Pr, A>>,
          boost::equality_comparable<sorted_vector<K, bNoDuplicates, Pr, A>>
    {
     public:
      typedef sorted_vector<K, bNoDuplicates, Pr, A> Myt_;
      typedef std::vector<K, A> Cont;
      typedef typename Cont::allocator_type allocator_type;
      typedef typename Cont::size_type size_type;
      typedef typename Cont::difference_type difference_type;
      typedef typename Cont::reference reference;
      typedef typename Cont::const_reference const_reference;
      typedef typename Cont::value_type value_type;
      typedef K key_type;
      typedef typename Cont::iterator iterator;
      typedef typename Cont::const_iterator const_iterator;
      typedef Pr key_compare;
      typedef Pr value_compare;

      typedef typename Cont::const_reverse_iterator const_reverse_iterator;
      typedef typename Cont::reverse_iterator reverse_iterator;

      typedef std::pair<iterator, iterator> Pairii_;
      typedef std::pair<const_iterator, const_iterator> Paircc_;
      typedef std::pair<iterator, bool> Pairib_;

      explicit sorted_vector(const Pr& pred = Pr(), const A& al = A()) : less_(pred), vec_(al) {}

      template <class It>
      sorted_vector(It first, It beyond, const Pr& pred = Pr(), const A& al = A())
          : less_(pred), vec_(first, beyond, al)
      {
        stable_sort();
      }

      sorted_vector(const Myt_& x) : less_(x.less_), vec_(x.vec_) {}

      ~sorted_vector() {}

      Myt_& operator=(const Myt_& x)
      {
        vec_.operator=(x.vec_);
        less_ = x.less_;
        return *this;
      }

      Myt_& operator=(const Cont& x)
      {
        vec_.operator=(x);
        sort();
        return *this;
      }

      void reserve(size_type n) { vec_.reserve(n); }
      iterator begin() { return vec_.begin(); }
      const_iterator begin() const { return vec_.begin(); }
      iterator end() { return vec_.end(); }
      const_iterator end() const { return vec_.end(); }
      reverse_iterator rbegin() { return vec_.rbegin(); }
      const_reverse_iterator rbegin() const { return vec_.rbegin(); }

      reverse_iterator rend() { return vec_.rend(); }
      const_reverse_iterator rend() const { return vec_.rend(); }

      size_type size() const { return vec_.size(); }
      size_type max_size() const { return vec_.max_size(); }
      bool empty() const { return vec_.empty(); }
      A get_allocator() const { return vec_.get_allocator(); }
      const_reference at(size_type p) const { return vec_.at(p); }
      reference at(size_type p) { return vec_.at(p); }
      const_reference operator[](size_type p) const { return vec_.operator[](p); }

      reference operator[](size_type p) { return vec_.operator[](p); }
      reference front() { return vec_.front(); }
      const_reference front() const { return vec_.front(); }
      reference back() { return vec_.back(); }
      const_reference back() const { return vec_.back(); }
      void pop_back() { vec_.pop_back(); }

      void assign(const_iterator first, const_iterator beyond) { vec_.assign(first, beyond); }

      void assign(size_type n, const K& x = K()) { vec_.assign(n, x); }

      /// insert members
      Pairib_ insert(const value_type& x)
      {
        if (bNoDuplicates)
        {
          iterator p = lower_bound(x);
          if (p == end() or less_(x, *p))
          {
            return Pairib_(InsertImpl_(p, x), true);
          }
          else
          {
            return Pairib_(p, false);
          }
        }
        else
        {
          iterator p = upper_bound(x);
          return Pairib_(InsertImpl_(p, x), true);
        }
      }

      iterator insert(iterator it, const value_type& x)  // it is the hint
      {
#if 0
      if ( it != end() )
      {
        if ( bNoDuplicates )
        {
          if ( less_( *it, x ) )
          {
            if ( ( it+1 ) == end() or Gt_( *(it+1), x ) )
            {
              // use hint
              return InsertImpl_(it+1,x);
            }
            else if ( Geq_( *(it+1), x ) )
            {
              return end();
            }
          }
        }
        else
        {
          if ( Leq_(*it,x) and ( (it+1)==end() or Geq_( *(it+1),x ) ) )
          {
            return InsertImpl_(it+1,x);
          }
        }
      }
#endif
        return insert(x).first;
      }

      template <class It>
      void insert(It first, It beyond)
      {
        size_type n = std::distance(first, beyond);
        reserve(size() + n);
        for (; first != beyond; ++first)
        {
          insert(*first);
        }
      }

      iterator ierase(iterator p) { return vec_.erase(p); }

      iterator erase(iterator first, iterator beyond) { return vec_.erase(first, beyond); }

      size_type erase(const K& key)
      {
        Pairii_ begEnd = equal_range(key);
        size_type n = std::distance(begEnd.first, begEnd.second);
        erase(begEnd.first, begEnd.second);
        return n;
      }

      void clear() { vec_.clear(); }

      void swap(Myt_& x)
      {
        vec_.swap(x.vec_);
        std::swap(less_, x.less_);
      }

      friend void swap(Myt_& x, Myt_& Y_) { x.swap(Y_); }

      key_compare key_comp() const { return less_; }
      value_compare value_comp() const { return (key_comp()); }

      iterator find(const K& k)
      {
        iterator p = lower_bound(k);
        return (p == end() or less_(k, *p)) ? end() : p;
      }

      const_iterator find(const K& k) const
      {
        const_iterator p = lower_bound(k);
        return (p == end() or less_(k, *p)) ? end() : p;
      }

      size_type count(const K& k) const
      {
        Paircc_ Ans_ = equal_range(k);
        size_type n = std::distance(Ans_.first, Ans_.second);
        return n;
      }

      iterator lower_bound(const K& k) { return std::lower_bound(begin(), end(), k, less_); }

      const_iterator lower_bound(const K& k) const
      {
        return std::lower_bound(begin(), end(), k, less_);
      }

      iterator upper_bound(const K& k) { return std::upper_bound(begin(), end(), k, less_); }

      const_iterator upper_bound(const K& k) const
      {
        return std::upper_bound(begin(), end(), k, less_);
      }

      Pairii_ equal_range(const K& k) { return std::equal_range(begin(), end(), k, less_); }

      Paircc_ equal_range(const K& k) const { return std::equal_range(begin(), end(), k, less_); }

      /// functions for use with direct std::vector-access
      Cont& get_container() { return vec_; }

      /// restore sorted order after low level access
      void sort()
      {
        std::sort(vec_.begin(), vec_.end(), less_);
        if (bNoDuplicates)
        {
          vec_.erase(Unique_(), vec_.end());
        }
      }

      /// restore sorted order after low level access
      void stable_sort()
      {
        std::stable_sort(vec_.begin(), vec_.end(), less_);
        if (bNoDuplicates)
        {
          erase(Unique_(), end());
        }
      }

      bool operator==(const sorted_vector<K, bNoDuplicates, Pr, A>& other) const
      {
        return Eq(other);
      }

      bool operator<(const sorted_vector<K, bNoDuplicates, Pr, A>& other) const
      {
        return Lt(other);
      }

     protected:
      iterator Unique_()
      {
        iterator front_ = vec_.begin(), out_ = vec_.end(), end_ = vec_.end();
        bool bCopy_ = false;
        for (iterator prev_; (prev_ = front_) != end_ and ++front_ != end_;)
        {
          if (less_(*prev_, *front_))
          {
            if (bCopy_)
            {
              *out_ = *front_;
              ++out_;
            }
          }
          else
          {
            if (not bCopy_)
            {
              out_ = front_;
              bCopy_ = true;
            }
          }
        }
        return out_;
      }

      bool Eq(const Myt_& x) const
      {
        return (size() == x.size() and std::equal(begin(), end(), x.begin()));
      }

      bool Lt(const Myt_& x) const
      {
        return (std::lexicographical_compare(begin(), end(), x.begin(), x.end()));
      }

      iterator InsertImpl_(iterator p, const value_type& x)
      {
        iterator i = vec_.insert(p, x);
#ifdef DEBUGCUTLIBRARY
        std::set<K, key_compare> testset(less_);
        testset.insert(begin(), end());
        if (testset.size() != vec_.size() or not std::equal(begin(), end(), testset.begin()))
          throw std::runtime_error("set invariance failed");
#endif
        return i;
      }

      bool Leq_(const K& ty0, const K& ty1) { return not less_(ty1, ty0); }

      bool Geq_(const K& ty0, const K& ty1) { return not less_(ty0, ty1); }

      bool Lt_(const K& ty0, const K& ty1) { return less_(ty0, ty1); }

      bool Gt_(const K& ty0, const K& ty1) { return less_(ty1, ty0); }

      key_compare less_;
      Cont vec_;
    };

  }  // namespace CUT
}  // namespace GEO

#endif
