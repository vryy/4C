/*----------------------------------------------------------------------*/
/*!
\file drt_condition_utils.H

\brief

<pre>
-------------------------------------------------------------------------
                 BACI finite element library subsystem
            Copyright (2008) Technical University of Munich

Under terms of contract T004.008.000 there is a non-exclusive license for use
of this work by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library is proprietary software. It must not be published, distributed,
copied or altered in any form or any media without written permission
of the copyright holder. It may be used under terms and conditions of the
above mentioned license by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library may solemnly used in conjunction with the BACI contact library
for purposes described in the above mentioned contract.

This library contains and makes use of software copyrighted by Sandia Corporation
and distributed under LGPL licence. Licensing does not apply to this or any
other third party software used here.

Questions? Contact Dr. Michael W. Gee (gee@lnm.mw.tum.de)
                   or
                   Prof. Dr. Wolfgang A. Wall (wall@lnm.mw.tum.de)

http://www.lnm.mw.tum.de

-------------------------------------------------------------------------
</pre>

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
      const string& condname
      );

    /// collect all local nodes and elements in a condition including ghosts
    void FindInterfaceObjects(const DRT::Discretization& dis,
                              map<int, RCP<DRT::Element> >& elements,
                              const string& condname);

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
