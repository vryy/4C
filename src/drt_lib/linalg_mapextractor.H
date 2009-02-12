/*----------------------------------------------------------------------*/
/*!
\file linalg_mapextractor.H

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
Maintainer: Ulrich Kuettler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15238
</pre>
*/
/*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifndef LINALG_MAPEXTRACTOR_H
#define LINALG_MAPEXTRACTOR_H

#include <map>
#include <set>
#include <string>
#include <vector>
#include <algorithm>

#include <Teuchos_RCP.hpp>
#include <Epetra_Import.h>
#include <Epetra_Map.h>
#include <Epetra_Vector.h>

/*!
\brief LINALG: namespace of the ccarat linear algebra module

*/
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
    MultiMapExtractor(const Epetra_Map& fullmap, const std::vector<Teuchos::RCP<const Epetra_Map> >& maps);

    /// setup of an empty extractor
    /*!
      \warning The fullmap has to be nonoverlapping. The list of maps has to
      be nonoverlapping as well and its summ has to equal the fullmap.
     */
    void Setup(const Epetra_Map& fullmap, const std::vector<Teuchos::RCP<const Epetra_Map> >& maps);

    /// merge set of unique maps
    /*!
      \warning There must be no overlap in these maps.
     */
    static Teuchos::RCP<Epetra_Map> MergeMaps(const std::vector<Teuchos::RCP<const Epetra_Map> >& maps);

    /// intersect set of unique maps
    /*!
      \warning There must be no overlap in these maps.
     */
    static Teuchos::RCP<Epetra_Map> IntersectMaps(const std::vector<Teuchos::RCP<const Epetra_Map> >& maps);

    /** \name Maps */
    //@{

    /// number of partial maps
    int NumMaps() const { return maps_.size(); }

    /// get the map
    const Teuchos::RCP<const Epetra_Map>& Map(int i) const { return maps_[i]; }

    /// the full map
    const Teuchos::RCP<const Epetra_Map>& FullMap() const { return fullmap_; }

    //@}

    /** \name Vector creation */
    //@{

    /// create vector to map i
    Teuchos::RCP<Epetra_Vector> Vector(int i) const
    { return Teuchos::rcp(new Epetra_Vector(*Map(i))); }

    /// create multi vector to map i
    Teuchos::RCP<Epetra_MultiVector> Vector(int i, int numvec) const
    { return Teuchos::rcp(new Epetra_MultiVector(*Map(i),numvec)); }

    //@

    /** \name Extract from full vector */
    //@{

    /// extract a partial vector from a full vector
    /*!
      \param full vector on the full map
      \param block block number of vector to extract
     */
    Teuchos::RCP<Epetra_Vector> ExtractVector(const Epetra_Vector& full, int block) const;

    /// extract a partial vector from a full vector
    /*!
      \param full vector on the full map
      \param block block number of vector to extract
     */
    Teuchos::RCP<Epetra_MultiVector> ExtractVector(const Epetra_MultiVector& full, int block) const;

    /// extract a partial vector from a full vector
    /*!
      \param full vector on the full map
      \param block block number of vector to extract
     */
    Teuchos::RCP<Epetra_Vector> ExtractVector(Teuchos::RCP<Epetra_Vector> full, int block) const
    { return ExtractVector(*full,block); }

    /// extract a partial vector from a full vector
    /*!
      \param full vector on the full map
      \param block block number of vector to extract
     */
    Teuchos::RCP<Epetra_MultiVector> ExtractVector(Teuchos::RCP<Epetra_MultiVector> full, int block) const
    { return ExtractVector(*full,block); }

    /// extract a partial vector from a full vector
    /*!
      \param full vector on the full map
      \param block block number of vector to extract
     */
    Teuchos::RCP<Epetra_Vector> ExtractVector(Teuchos::RCP<const Epetra_Vector> full, int block) const
    { return ExtractVector(*full,block); }

    /// extract a partial vector from a full vector
    /*!
      \param full vector on the full map
      \param block block number of vector to extract
     */
    Teuchos::RCP<Epetra_MultiVector> ExtractVector(Teuchos::RCP<const Epetra_MultiVector> full, int block) const
    { return ExtractVector(*full,block); }

    /// extract a partial vector from a full vector
    /*!
      \param full vector on the full map
      \param block block number of vector to extract
      \param partial partial vector to fill
     */
    void ExtractVector(const Epetra_MultiVector& full, int block, Epetra_MultiVector& partial) const;

    /// extract a partial vector from a full vector
    /*!
      \param full vector on the full map
      \param block block number of vector to extract
      \param partial partial vector to fill
     */
    void ExtractVector(Teuchos::RCP<const Epetra_Vector> full, int block, Teuchos::RCP<Epetra_Vector> partial) const
    { ExtractVector(*full,block,*partial); }

    //@}

    /** \name Insert from full dof vector */
    //@{

    /// Put a partial vector into a full vector
    /*!
      \param partial vector to copy into full vector
      \param block block number of partial vector
     */
    Teuchos::RCP<Epetra_Vector> InsertVector(const Epetra_Vector& partial, int block) const;

    /// Put a partial vector into a full vector
    /*!
      \param partial vector to copy into full vector
      \param block block number of partial vector
     */
    Teuchos::RCP<Epetra_MultiVector> InsertVector(const Epetra_MultiVector& partial, int block) const;

    /// Put a partial vector into a full vector
    /*!
      \param partial vector to copy into full vector
      \param block block number of partial vector
     */
    Teuchos::RCP<Epetra_Vector> InsertVector(Teuchos::RCP<const Epetra_Vector> partial, int block) const
    { return InsertVector(*partial,block); }

    /// Put a partial vector into a full vector
    /*!
      \param partial vector to copy into full vector
      \param block block number of partial vector
     */
    Teuchos::RCP<Epetra_MultiVector> InsertVector(Teuchos::RCP<const Epetra_MultiVector> partial, int block) const
    { return InsertVector(*partial,block); }

    /// Put a partial vector into a full vector
    /*!
      \param partial vector to copy into full vector
      \param block block number of partial vector
     */
    Teuchos::RCP<Epetra_Vector> InsertVector(Teuchos::RCP<Epetra_Vector> partial, int block) const
    { return InsertVector(*partial,block); }

    /// Put a partial vector into a full vector
    /*!
      \param partial vector to copy into full vector
      \param block block number of partial vector
     */
    Teuchos::RCP<Epetra_MultiVector> InsertVector(Teuchos::RCP<Epetra_MultiVector> partial, int block) const
    { return InsertVector(*partial,block); }

    /// Put a partial vector into a full vector
    /*!
      \param partial vector to copy into full vector
      \param block block number of partial vector
      \param full vector to copy into
     */
    void InsertVector(const Epetra_MultiVector& partial, int block, Epetra_MultiVector& full) const;

    /// Put a partial vector into a full vector
    /*!
      \param partial vector to copy into full vector
      \param block block number of partial vector
      \param full vector to copy into
     */
    void InsertVector(Teuchos::RCP<const Epetra_Vector> partial, int block, Teuchos::RCP<Epetra_Vector> full) const
    { InsertVector(*partial,block,*full); }

    //@}

    /** \name Add from full dof vector */
    //@{

    /// Put a partial vector into a full vector
    /*!
      \param partial vector to copy into full vector
      \param block block number of partial vector
      \param full vector to copy into
      \param scale scaling factor for partial vector
     */
    void AddVector(const Epetra_MultiVector& partial, int block, Epetra_MultiVector& full, double scale=1.0) const;

    /// Put a partial vector into a full vector
    /*!
      \param partial vector to copy into full vector
      \param block block number of partial vector
      \param full vector to copy into
      \param scale scaling factor for partial vector
     */
    void AddVector(Teuchos::RCP<const Epetra_Vector> partial, int block, Teuchos::RCP<Epetra_Vector> full, double scale=1.0) const
    { AddVector(*partial,block,*full,scale); }

    //@}

  private:

    /// the full row map
    Teuchos::RCP<const Epetra_Map> fullmap_;

    /// the list of nonoverlapping partial row maps that sums up to the full map
    std::vector<Teuchos::RCP<const Epetra_Map> > maps_;

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

  \note The two partial maps (cond and other) are stored in the parent member variable maps_,
  where other has index 0 and cond has index 1.

  \author u.kue
  \date 01/08
  */
  class MapExtractor : public MultiMapExtractor
  {
  public:

    MapExtractor();

    MapExtractor(const Epetra_Map& fullmap, Teuchos::RCP<const Epetra_Map> condmap, Teuchos::RCP<const Epetra_Map> othermap);

    /// constructor builds non-overlapping othermap which is complementary to condmap with respect to fullmap
    /// \author bborn
    /// \date 10/08
    MapExtractor(const Epetra_Map& fullmap, Teuchos::RCP<const Epetra_Map> condmap);

    /** \name Setup */
    //@{

    /// setup from known maps
    void Setup(const Epetra_Map& fullmap, Teuchos::RCP<const Epetra_Map> condmap, Teuchos::RCP<const Epetra_Map> othermap);

    //@}

    /** \name Maps */
    //@{

    /// get the condition dof map
    const Teuchos::RCP<const Epetra_Map>& CondMap() const { return Map(1); }

    /// get the other dof map
    const Teuchos::RCP<const Epetra_Map>& OtherMap() const { return Map(0); }

    //@}

    /// tell if this extractor has any meaning
    /*!
      We tend to setup extractors even if there is no condition at all. In
      such a case we will end up with an empty condition map and there is no
      point to the extractor all together.
     */
    bool Relevant() const { return CondMap()->NumGlobalElements()!=0; }

    /** \name Extract from full dof vector */
    //@{

    /// extract a condition dof vector from a full row (or column) vector
    Teuchos::RCP<Epetra_Vector> ExtractCondVector(const Epetra_Vector& full) const
    { return MultiMapExtractor::ExtractVector(full,1); }

    /// extract an other dof vector from a full row (or column) vector
    Teuchos::RCP<Epetra_Vector> ExtractOtherVector(const Epetra_Vector& full) const
    { return MultiMapExtractor::ExtractVector(full,0); }

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

    /** \name Add from full dof vector */
    /// These are convenience methods that create a temporary vector internally.

    /// Put a condition dof vector into a full row vector
    void AddCondVector(Teuchos::RCP<const Epetra_Vector> cond, Teuchos::RCP<Epetra_Vector> full) const;

    /// Put a scaled condition dof vector into a full row vector
    void AddCondVector(double scale, Teuchos::RCP<const Epetra_Vector> cond, Teuchos::RCP<Epetra_Vector> full) const;

    /// Put an other dof vector into a full row vector
    void AddOtherVector(Teuchos::RCP<const Epetra_Vector> other, Teuchos::RCP<Epetra_Vector> full) const;

    /// Put a scaled other dof vector into a full row vector
    void AddOtherVector(double scale,Teuchos::RCP<const Epetra_Vector> other, Teuchos::RCP<Epetra_Vector> full) const;

    //@}

  private:

  };

}


#endif
#endif
