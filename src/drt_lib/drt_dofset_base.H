/*!----------------------------------------------------------------------
\file drt_dofset_base.H

\brief Declaration A set of degrees of freedom

<pre>
\brief Declaration
\level 0
\maintainer Martin Kronbichler
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
\brief Declaration Base class set of degrees of freedom

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
  \brief Declaration Standard Constructor

  */
  DofSetBase();

  /*!
  \brief Declaration Destructor

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
