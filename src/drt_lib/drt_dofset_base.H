/*!----------------------------------------------------------------------
\file drt_dofset_base.H

\brief A set of degrees of freedom

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
Maintainer: Thomas Kloeppel
            kloeppel@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
</pre>

*----------------------------------------------------------------------*/
#ifndef DOFSETBASE_H
#define DOFSETBASE_H

#include <list>

#include <Teuchos_RCP.hpp>
#include <Epetra_Map.h>

namespace DRT
{

class Discretization;


/*!
\brief Base class set of degrees of freedom

\note This is a pure virtual class all classes managing sets of degrees of freedom
should inherit from. This base class only manages the static list, all DofSets are
written into.

\author tk
*/
class DofSetBase
{
  friend class DofSetProxy;
public:

  /*!
  \brief Standard Constructor

  */
  DofSetBase();

  /*!
  \brief Destructor

  */
  virtual ~DofSetBase();


  //! @name Access methods

  ///Get Number of Global Elements of degree of freedom row map
  virtual int NumGlobalElements() const;

  /// Get maximum GID of degree of freedom row map
  virtual int MaxAllGID() const;

  /// Get minimum GID of degree of freedom row map
  virtual int MinAllGID() const;

  /// Get degree of freedom row map
  virtual const Epetra_Map* DofRowMap() const;

  /// Get degree of freedom column map
  virtual const Epetra_Map* DofColMap() const;

  /// are the dof maps already initialized?
  virtual bool Initialized() const;

  /// Get Max of all GID assigned in the DofSets in front of current one in the list #static_dofsets_
  virtual int MaxGIDinList(const Epetra_Comm& comm) const;

  //@}

  //! @name Construction

  /// Add Dof Set to list #static_dofsets_
  virtual void AddDofSettoList();

  /// Replace a Dof Set in list #static_dofsets_ with this
  virtual void ReplaceInStaticDofsets(Teuchos::RCP<DofSetBase> olddofset);

  /// reset all internal variables
  virtual void Reset()=0;

  //@}

  /// Print the dofsets in the static_dofsets_ list
  virtual void PrintAllDofsets(const Epetra_Comm& comm) const;

protected:

  /// unique row map of degrees of freedom (node, face, and element dofs))
  Teuchos::RCP<Epetra_Map>                       dofrowmap_;

  /// unique column map of degrees of freedom (node, face, and element dofs)
  Teuchos::RCP<Epetra_Map>                       dofcolmap_;

private:

  /// store dofset in static list, if derived class chooses so using AddDofSettoList()
  /*!
    This is hack to get unique dof numbers on all dof sets. With this in place
    we can combine the maps from any dof sets to form block systems and the
    like.
   */
  static std::list<DofSetBase*>                 static_dofsets_;

}; // class DofSetBase
} // namespace DRT

#endif  // #ifndef DofSetBase_H
