/*----------------------------------------------------------------------*/
/*!
\file linalg_mapextractor.H

\brief Declaration of map extractor class

<pre>
\level 0
\maintainer Martin Kronbichler
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>
*/
/*----------------------------------------------------------------------*/
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

namespace LINALG
{
  /// Split a row map into a set of partial maps and establish the communication pattern back and
  /// forth
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

    /// destructor
    virtual ~MultiMapExtractor(){};

    /// create an extractor from fullmap to the given set of maps
    MultiMapExtractor(
        const Epetra_Map& fullmap, const std::vector<Teuchos::RCP<const Epetra_Map>>& maps);

    /// setup of an empty extractor
    /*!
      \warning The fullmap has to be nonoverlapping. The list of maps has to
      be nonoverlapping as well and its sum has to equal the fullmap.
     */
    void Setup(const Epetra_Map& fullmap, const std::vector<Teuchos::RCP<const Epetra_Map>>& maps);

    /// debug helper
    /*!
      loop all maps in the list of nonoverlapping partial row maps unequal
      Teuchos::null and check whether they have a valid DataPtr() and are
      unique

      \note hidden calls to Redistribute may render maps in maps_ obsolete.
      This function is intendend to simplify debugging for these cases.
    */
    void CheckForValidMapExtractor() const;

    /// merge set of unique maps
    /*!
      \warning There must be no overlap in these maps.
               The order of the GIDs is destroyed
     */
    static Teuchos::RCP<Epetra_Map> MergeMaps(
        const std::vector<Teuchos::RCP<const Epetra_Map>>& maps);

    /// merge set of unique maps
    /*!
      \warning There must be no overlap in these maps.
    */
    static Teuchos::RCP<Epetra_Map> MergeMapsKeepOrder(
        const std::vector<Teuchos::RCP<const Epetra_Map>>& maps);

    /// intersect set of unique maps
    /*!
      \warning There must be no overlap in these maps.
     */
    static Teuchos::RCP<Epetra_Map> IntersectMaps(
        const std::vector<Teuchos::RCP<const Epetra_Map>>& maps);

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
    {
      return Teuchos::rcp(new Epetra_Vector(*Map(i)));
    }

    /// create multi vector to map i
    Teuchos::RCP<Epetra_MultiVector> Vector(int i, int numvec) const
    {
      return Teuchos::rcp(new Epetra_MultiVector(*Map(i), numvec));
    }

    //@

    /** \name Extract from full vector */
    //@{

    /// extract a partial vector from a full vector
    /*!
      \param full vector on the full map
      \param block number of vector to extract
     */
    Teuchos::RCP<Epetra_Vector> ExtractVector(const Epetra_Vector& full, int block) const;

    /// extract a partial vector from a full vector
    /*!
      \param full vector on the full map
      \param block number of vector to extract
     */
    Teuchos::RCP<Epetra_MultiVector> ExtractVector(const Epetra_MultiVector& full, int block) const;

    /// extract a partial vector from a full vector
    /*!
      \param full vector on the full map
      \param block number of vector to extract
     */
    Teuchos::RCP<Epetra_Vector> ExtractVector(Teuchos::RCP<Epetra_Vector> full, int block) const
    {
      return ExtractVector(*full, block);
    }

    /// extract a partial vector from a full vector
    /*!
      \param full vector on the full map
      \param block number of vector to extract
     */
    Teuchos::RCP<Epetra_MultiVector> ExtractVector(
        Teuchos::RCP<Epetra_MultiVector> full, int block) const
    {
      return ExtractVector(*full, block);
    }

    /// extract a partial vector from a full vector
    /*!
      \param full vector on the full map
      \param block number of vector to extract
     */
    Teuchos::RCP<Epetra_Vector> ExtractVector(
        Teuchos::RCP<const Epetra_Vector> full, int block) const
    {
      return ExtractVector(*full, block);
    }

    /// extract a partial vector from a full vector
    /*!
      \param full vector on the full map
      \param block number of vector to extract
     */
    Teuchos::RCP<Epetra_MultiVector> ExtractVector(
        Teuchos::RCP<const Epetra_MultiVector> full, int block) const
    {
      return ExtractVector(*full, block);
    }

    /// extract a partial vector from a full vector
    /*!
      \param full vector on the full map
      \param block number of vector to extract
      \param partial vector to fill
     */
    virtual void ExtractVector(
        const Epetra_MultiVector& full, int block, Epetra_MultiVector& partial) const;

    /// extract a partial vector from a full vector
    /*!
      \param full vector on the full map
      \param block number of vector to extract
      \param partial vector to fill
     */
    void ExtractVector(Teuchos::RCP<const Epetra_Vector> full, int block,
        Teuchos::RCP<Epetra_Vector> partial) const
    {
      ExtractVector(*full, block, *partial);
    }

    //@}

    /** \name Insert from full dof vector */
    //@{

    /// Put a partial vector into a full vector
    /*!
      \param partial vector to copy into full vector
      \param block number of partial vector
     */
    Teuchos::RCP<Epetra_Vector> InsertVector(const Epetra_Vector& partial, int block) const;

    /// Put a partial vector into a full vector
    /*!
      \param partial vector to copy into full vector
      \param block number of partial vector
     */
    Teuchos::RCP<Epetra_MultiVector> InsertVector(
        const Epetra_MultiVector& partial, int block) const;

    /// Put a partial vector into a full vector
    /*!
      \param partial vector to copy into full vector
      \param block number of partial vector
     */
    Teuchos::RCP<Epetra_Vector> InsertVector(
        Teuchos::RCP<const Epetra_Vector> partial, int block) const
    {
      return InsertVector(*partial, block);
    }

    /// Put a partial vector into a full vector
    /*!
      \param partial vector to copy into full vector
      \param block number of partial vector
     */
    Teuchos::RCP<Epetra_MultiVector> InsertVector(
        Teuchos::RCP<const Epetra_MultiVector> partial, int block) const
    {
      return InsertVector(*partial, block);
    }

    /// Put a partial vector into a full vector
    /*!
      \param partial vector to copy into full vector
      \param block number of partial vector
     */
    Teuchos::RCP<Epetra_Vector> InsertVector(Teuchos::RCP<Epetra_Vector> partial, int block) const
    {
      return InsertVector(*partial, block);
    }

    /// Put a partial vector into a full vector
    /*!
      \param partial vector to copy into full vector
      \param block number of partial vector
     */
    Teuchos::RCP<Epetra_MultiVector> InsertVector(
        Teuchos::RCP<Epetra_MultiVector> partial, int block) const
    {
      return InsertVector(*partial, block);
    }

    /// Put a partial vector into a full vector
    /*!
      \param partial vector to copy into full vector
      \param block number of partial vector
      \param full vector to copy into
     */
    virtual void InsertVector(
        const Epetra_MultiVector& partial, int block, Epetra_MultiVector& full) const;

    /// Put a partial vector into a full vector
    /*!
      \param partial vector to copy into full vector
      \param block number of partial vector
      \param full vector to copy into
     */
    void InsertVector(Teuchos::RCP<const Epetra_Vector> partial, int block,
        Teuchos::RCP<Epetra_Vector> full) const
    {
      InsertVector(*partial, block, *full);
    }

    //@}

    /** \name Add from full dof vector */
    //@{

    /// Put a partial vector into a full vector
    /*!
      \param partial vector to copy into full vector
      \param block number of partial vector
      \param full vector to copy into
      \param scale scaling factor for partial vector
     */
    virtual void AddVector(const Epetra_MultiVector& partial, int block, Epetra_MultiVector& full,
        double scale = 1.0) const;

    /// Put a partial vector into a full vector
    /*!
      \param partial vector to copy into full vector
      \param block number of partial vector
      \param full vector to copy into
      \param scale scaling factor for partial vector
     */
    void AddVector(Teuchos::RCP<const Epetra_Vector> partial, int block,
        Teuchos::RCP<Epetra_Vector> full, double scale = 1.0) const
    {
      AddVector(*partial, block, *full, scale);
    }

    //@}

    /// PutScalar to one block only
    void PutScalar(Epetra_Vector& full, int block, double scalar) const;

    /// L2-norm of one block only
    double Norm2(const Epetra_Vector& full, int block) const;

    /// Scale one block only
    void Scale(Epetra_Vector& full, int block, double scalar) const;

    /// Scale one block only
    void Scale(Epetra_MultiVector& full, int block, double scalar) const;

   protected:
    /// the full row map
    Teuchos::RCP<const Epetra_Map> fullmap_;

    /// the list of nonoverlapping partial row maps that sums up to the full map
    std::vector<Teuchos::RCP<const Epetra_Map>> maps_;

    /// communication between condition dof map and full row dof map
    std::vector<Teuchos::RCP<Epetra_Import>> importer_;
  };


/// Add all kinds of support methods to derived classes of MultiMapExtractor.
#define MAP_EXTRACTOR_VECTOR_METHODS(name, pos)                                                    \
  Teuchos::RCP<Epetra_Vector> Extract##name##Vector(const Epetra_Vector& full) const               \
  {                                                                                                \
    return MultiMapExtractor::ExtractVector(full, pos);                                            \
  }                                                                                                \
                                                                                                   \
  Teuchos::RCP<Epetra_Vector> Extract##name##Vector(Teuchos::RCP<const Epetra_Vector> full) const  \
  {                                                                                                \
    return MultiMapExtractor::ExtractVector(full, pos);                                            \
  }                                                                                                \
                                                                                                   \
  void Extract##name##Vector(                                                                      \
      Teuchos::RCP<const Epetra_Vector> full, Teuchos::RCP<Epetra_Vector> cond) const              \
  {                                                                                                \
    ExtractVector(full, pos, cond);                                                                \
  }                                                                                                \
                                                                                                   \
  Teuchos::RCP<Epetra_Vector> Insert##name##Vector(Teuchos::RCP<const Epetra_Vector> cond) const   \
  {                                                                                                \
    return InsertVector(cond, pos);                                                                \
  }                                                                                                \
                                                                                                   \
  void Insert##name##Vector(                                                                       \
      Teuchos::RCP<const Epetra_Vector> cond, Teuchos::RCP<Epetra_Vector> full) const              \
  {                                                                                                \
    InsertVector(cond, pos, full);                                                                 \
  }                                                                                                \
                                                                                                   \
  void Add##name##Vector(Teuchos::RCP<const Epetra_Vector> cond, Teuchos::RCP<Epetra_Vector> full) \
      const                                                                                        \
  {                                                                                                \
    AddVector(cond, pos, full);                                                                    \
  }                                                                                                \
                                                                                                   \
  void Add##name##Vector(double scale, Teuchos::RCP<const Epetra_Vector> cond,                     \
      Teuchos::RCP<Epetra_Vector> full) const                                                      \
  {                                                                                                \
    AddVector(cond, pos, full, scale);                                                             \
  }                                                                                                \
                                                                                                   \
  const Teuchos::RCP<const Epetra_Map>& name##Map() const { return Map(pos); }                     \
                                                                                                   \
  bool name##Relevant() const { return name##Map()->NumGlobalElements() != 0; }                    \
                                                                                                   \
  void name##PutScalar(Epetra_Vector& full, double scalar) const { PutScalar(full, pos, scalar); } \
                                                                                                   \
  double name##Norm2(const Epetra_Vector& full) const { return Norm2(full, pos); }


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
    /** \brief empty constructor
     *
     *  You have to call a Setup() routine of your choice. */
    MapExtractor();

    /** \brief  constructor
     *
     *  Calls Setup() from known maps */
    MapExtractor(const Epetra_Map& fullmap, Teuchos::RCP<const Epetra_Map> condmap,
        Teuchos::RCP<const Epetra_Map> othermap);

    /** \brief constructor
     *
     *  Calls Setup() to create non-overlapping othermap/condmap which is complementary
     *  to condmap/othermap with respect to fullmap depending on boolean 'iscondmap'  */
    MapExtractor(const Epetra_Map& fullmap,         //< full map
        Teuchos::RCP<const Epetra_Map> partialmap,  //< partial map, ie condition or other map
        bool iscondmap = true                       //< true, if partialmap is condition map
    );

    /** \name Setup */
    //@{

    /// setup from known maps
    void Setup(const Epetra_Map& fullmap, const Teuchos::RCP<const Epetra_Map>& condmap,
        const Teuchos::RCP<const Epetra_Map>& othermap);

    /// setup creates non-overlapping othermap/condmap which is complementary to condmap/othermap
    /// with respect to fullmap depending on boolean 'iscondmap'
    /// \author bborn
    /// \date 10/08
    void Setup(const Epetra_Map& fullmap, const Teuchos::RCP<const Epetra_Map>& partialmap,
        bool iscondmap = true);

    //@}

    MAP_EXTRACTOR_VECTOR_METHODS(Cond, 1)
    MAP_EXTRACTOR_VECTOR_METHODS(Other, 0)

   private:
  };

}  // namespace LINALG


#endif
