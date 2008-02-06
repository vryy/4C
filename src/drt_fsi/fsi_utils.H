
#ifdef CCADISCRET

#ifndef FSI_UTILS_H
#define FSI_UTILS_H

#include <set>
#include <string>
#include <functional>
#include <iostream>

#include <NOX_Epetra_Interface_Required.H>
#include <NOX_Abstract_Vector.H>
#include <NOX_Abstract_Group.H>
#include <Teuchos_RefCountPtr.hpp>

#include "../drt_lib/drt_discret.H"


namespace FSI
{
  /// FSI::UTILS: Random stuff that might be helpful when dealing with FSI problems
  namespace UTILS
  {

    /// return the first slot of a pair
    /*!
      To be used with stl algorithms.

      This should be part of stl but is not (yet). So we define our
      own version.
    */
    template <typename pair_type>
    struct select1st: public std::unary_function<const pair_type &,
                                                 const typename pair_type::first_type &>
    {
      const typename pair_type::first_type &operator()(const pair_type &v) const
        {
          return v.first;
        }
    };

    /// return the second slot of a pair
    /*!
      To be used with stl algorithms.

      This should be part of stl but is not (yet). So we define our
      own version.
    */
    template <typename pair_type>
    struct select2nd: public std::unary_function<const pair_type &,
                                                 const typename pair_type::second_type &>
    {
      const typename pair_type::second_type &operator()(const pair_type &v) const
        {
          return v.second;
        }
    };

    /// std unary function version of Epetra_Map::MyGID()
    struct MyGID : public std::unary_function<int, bool>
    {
      const Epetra_Map* emap_;
      MyGID(const Epetra_Map* emap) : emap_(emap) {}
      bool operator()(int gid) const
        {
          return emap_->MyGID(gid);
        }
    };

    /// write interface jacobian to file
    void DumpJacobian(NOX::Epetra::Interface::Required& interface,
                      double alpha,
                      double beta,
                      Teuchos::RefCountPtr<Epetra_Vector> soln,
                      std::string filename);

  }
}

#endif
#endif
