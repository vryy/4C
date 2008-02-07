#ifdef CCADISCRET
#ifndef DRT_UTILS_MAPEXTRACTOR_H
#define DRT_UTILS_MAPEXTRACTOR_H

#include <map>
#include <set>
#include <string>
#include <vector>
#include <algorithm>

#include <Teuchos_RCP.hpp>
#include <Epetra_Import.h>
#include <Epetra_Map.h>
#include <Epetra_Vector.h>


namespace LINALG
{

  /// Split a row map into a set of partial maps and establish the communication pattern back and forth
  /*!

    A general purpose class that contains a nonoverlapping full map and a set
    of partial maps. The sum of all partial maps is equals the full map. There
    is no overlap, neither within the partial maps nor between them.

    Communication from full vectors to partial vectors is supported.

    \note The MultiMapExtractor does not do the actual splitting. Thus no
    assumption on the items covered by the maps is made. The actual splitting
    has to be performed by the user.

    \author u.kue
    \date 02/08
   */
  class MultiMapExtractor
  {
  public:

    /// create an uninitialized (empty) extractor
    MultiMapExtractor();

    /// create an extractor from fullmap to the given set of maps
    MultiMapExtractor(const Epetra_Map& fullmap, const std::vector<Teuchos::RCP<Epetra_Map> >& maps);

    /// setup of an empty extractor
    /*!
      \warning The fullmap has to be nonoverlapping. The list of maps has to
      be nonoverlapping as well and its summ has to equal the fullmap.
     */
    void Setup(const Epetra_Map& fullmap, const std::vector<Teuchos::RCP<Epetra_Map> >& maps);

    /** \name Maps */
    //@{

    /// number of partial maps
    int NumMaps() const { return maps_.size(); }

    /// get the map
    const Teuchos::RCP<Epetra_Map>& Map(int i) const { return maps_[i]; }

    //@}

    /** \name Extract from full vector */
    //@{

    /// extract a partial vector from a full vector
    Teuchos::RCP<Epetra_Vector> ExtractVector(Teuchos::RCP<const Epetra_Vector> full, int i) const;

    /// extract a partial vector from a full vector
    void ExtractVector(Teuchos::RCP<const Epetra_Vector> full, int i, Teuchos::RCP<Epetra_Vector> cond) const;

    //@}

    /** \name Insert from full dof vector */
    //@{

    /// Put a partial vector into a full vector
    Teuchos::RCP<Epetra_Vector> InsertVector(Teuchos::RCP<const Epetra_Vector> cond, int i) const;

    /// Put a partial vector into a full vector
    void InsertVector(Teuchos::RCP<const Epetra_Vector> cond, int i, Teuchos::RCP<Epetra_Vector> full) const;

    //@}

  private:

    /// the full row map
    Teuchos::RCP<Epetra_Map> fullmap_;

    /// the list of nonoverlapping partial row maps that sums up to the full map
    std::vector<Teuchos::RCP<Epetra_Map> > maps_;

    /// communication between condition dof map and full row dof map
    std::vector<Teuchos::RCP<Epetra_Import> > importer_;
  };


  /// Split a dof row map in two and establish the communication pattern between those maps
  /*!

  Special convenience version of MultiMapExtractor that knows exactly two
  partial maps.

  Examples of such splits include the velocity -- pressure split of the dof
  row map of a fluid problem or the interface -- interior split in FSI
  problems. Many more examples are possible. This is the class to use each
  time a submap needs to be managed.

  \note We work on row maps. The maps we deal with are meant to be
  nonoverlapping.

  At the core there are the CondMap(), the map of all selected dofs, and
  OtherMap(), the map of all remaining dofs. This duality also exists in
  the extraction methods ExtractCondVector() and ExtractOtherVector(), that
  extract a subvector from a full one, and the insertion methods
  InsertCondVector() and InsertOtherVector(), that copy a subvector into a
  full vector. These extractions and insertions are termed communications,
  because internally an Epetra_Import class is used, even though there is no
  communication required once the Epetra_Import object is created.

  \author u.kue
  \date 01/08
  */
  class MapExtractor : public MultiMapExtractor
  {
  public:

    MapExtractor();

    MapExtractor(const Epetra_Map& fullmap, Teuchos::RCP<Epetra_Map> condmap, Teuchos::RCP<Epetra_Map> othermap);

    /** \name Setup */
    //@{

    /// setup from known maps
    void Setup(const Epetra_Map& fullmap, Teuchos::RCP<Epetra_Map> condmap, Teuchos::RCP<Epetra_Map> othermap);

    //@}

    /** \name Maps */
    //@{

    /// get the condition dof map
    const Teuchos::RCP<Epetra_Map>& CondMap() const { return Map(1); }

    /// get the other dof map
    const Teuchos::RCP<Epetra_Map>& OtherMap() const { return Map(0); }

    //@}

    /** \name Extract from full dof vector */
    //@{

    /// extract a condition dof vector from a full row (or column) vector
    Teuchos::RCP<Epetra_Vector> ExtractCondVector(Teuchos::RCP<const Epetra_Vector> full) const;

    /// extract an other dof vector from a full row (or column) vector
    Teuchos::RCP<Epetra_Vector> ExtractOtherVector(Teuchos::RCP<const Epetra_Vector> full) const;

    /// extract a condition dof vector from a full row (or column) vector
    void ExtractCondVector(Teuchos::RCP<const Epetra_Vector> full, Teuchos::RCP<Epetra_Vector> cond) const;

    /// extract an other dof vector from a full row (or column) vector
    void ExtractOtherVector(Teuchos::RCP<const Epetra_Vector> full, Teuchos::RCP<Epetra_Vector> other) const;

    //@}

    /** \name Insert from full dof vector */
    //@{

    /// Put a condition dof vector into a full row vector
    Teuchos::RCP<Epetra_Vector> InsertCondVector(Teuchos::RCP<const Epetra_Vector> cond) const;

    /// Put an other dof vector into a full row vector
    Teuchos::RCP<Epetra_Vector> InsertOtherVector(Teuchos::RCP<const Epetra_Vector> other) const;

    /// Put a condition dof vector into a full row vector
    void InsertCondVector(Teuchos::RCP<const Epetra_Vector> cond, Teuchos::RCP<Epetra_Vector> full) const;

    /// Put an other dof vector into a full row vector
    void InsertOtherVector(Teuchos::RCP<const Epetra_Vector> other, Teuchos::RCP<Epetra_Vector> full) const;

    //@}

  private:

  };

}


#endif
#endif
