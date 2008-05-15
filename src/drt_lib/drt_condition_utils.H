/*----------------------------------------------------------------------*/
/*!
\file drt_condition_utils.H

\brief

<pre>
Maintainer: Axel Gerstenberger
            gerstenberger@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>
*/
/*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#ifndef DRT_CONDITION_UTILS_H
#define DRT_CONDITION_UTILS_H

#include "drt_discret.H"
#include "linalg_mapextractor.H"

namespace DRT
{
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
    
    /// collect all local nodes and elements in a condition, e.g. an interface
    void FindInterfaceObjects(
      const DRT::Discretization& dis,
      map<int, DRT::Node*>& nodes,
      map< int, RefCountPtr<DRT::Element> >& elements,
      const string&              condname
      );
    
    /// create ale discretization parallel to the fluid one
    void CreateAleDiscretization();
    
    //! create discretization from a given condition
    /*
     * This file is the right place to have such functionionality, but it needs to be generalized
     * A. Gerstenberger 05/08
     */
    Teuchos::RCP<DRT::Discretization> CreateDiscretizationFromCondition(
      Teuchos::RCP<DRT::Discretization> sourcedis,    ///< discretization with condition
      const std::string&                condname,     ///< name of the condition, by which the derived discretization is identified
      const std::string&                discret_name, ///< name of the new condition
      const std::string&                element_name,  ///< name/type of the elements to be created
      const std::vector<std::string>&   conditions_to_copy  ///< list of conditions that will be copied to the new discretization
      );
    
    //! Given a nodal row and coloumn map, this routines creates the corresponding element maps and call FillComplete()
    void RedistributeWithNewNodalDistribution(
        DRT::Discretization&     dis,
        const Epetra_Map&        noderowmap,
        const Epetra_Map&        nodecolmap
        );
  }
}

#endif
#endif
