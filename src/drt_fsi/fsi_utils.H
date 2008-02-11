
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
#include "../drt_lib/linalg_mapextractor.H"


namespace FSI
{
  /// FSI::UTILS: Random stuff that might be helpful when dealing with FSI problems
  namespace UTILS
  {

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

    ///
    void SetupInterfaceExtractor(const DRT::Discretization& dis,
                                 std::string condname,
                                 LINALG::MapExtractor& extractor);
  }
}

#endif
#endif
