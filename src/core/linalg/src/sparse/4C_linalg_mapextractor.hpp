/*----------------------------------------------------------------------*/
/*! \file

\brief Declaration of map extractor class

\level 0
*/
/*----------------------------------------------------------------------*/
#ifndef FOUR_C_LINALG_MAPEXTRACTOR_HPP
#define FOUR_C_LINALG_MAPEXTRACTOR_HPP

#include "4C_config.hpp"

#include <Epetra_Import.h>
#include <Epetra_Map.h>
#include <Epetra_Vector.h>
#include <Teuchos_RCP.hpp>

#include <algorithm>
#include <map>
#include <set>
#include <string>
#include <vector>

FOUR_C_NAMESPACE_OPEN

namespace Core::LinAlg
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
    virtual ~MultiMapExtractor() = default;

    /// create an extractor from fullmap to the given set of maps
    MultiMapExtractor(
        const Epetra_Map& fullmap, const std::vector<Teuchos::RCP<const Epetra_Map>>& maps);

    /// setup of an empty extractor
    /*!
      \warning The fullmap has to be nonoverlapping. The list of maps has to
      be nonoverlapping as well and its sum has to equal the fullmap.
     */
    void setup(const Epetra_Map& fullmap, const std::vector<Teuchos::RCP<const Epetra_Map>>& maps);

    /// debug helper
    /*!
      loop all maps in the list of nonoverlapping partial row maps unequal
      Teuchos::null and check whether they have a valid DataPtr() and are
      unique

      \note hidden calls to redistribute() may render maps in maps_ obsolete.
      This function is intendend to simplify debugging for these cases.
    */
    void check_for_valid_map_extractor() const;

    /// merge set of unique maps
    /*!
      \warning There must be no overlap in these maps.
               The order of the GIDs is destroyed
     */
    static Teuchos::RCP<Epetra_Map> merge_maps(
        const std::vector<Teuchos::RCP<const Epetra_Map>>& maps);

    /// merge set of unique maps
    /*!
      \warning There must be no overlap in these maps.
    */
    static Teuchos::RCP<Epetra_Map> merge_maps_keep_order(
        const std::vector<Teuchos::RCP<const Epetra_Map>>& maps);

    /// intersect set of unique maps
    /*!
      \warning There must be no overlap in these maps.
     */
    static Teuchos::RCP<Epetra_Map> intersect_maps(
        const std::vector<Teuchos::RCP<const Epetra_Map>>& maps);

    /** \name Maps */
    //@{

    /// number of partial maps
    int num_maps() const { return maps_.size(); }

    /// get the map
    const Teuchos::RCP<const Epetra_Map>& Map(int i) const { return maps_[i]; }

    /// the full map
    const Teuchos::RCP<const Epetra_Map>& full_map() const { return fullmap_; }

    //@}

    /** \name Vector creation */
    //@{

    /// create vector to map i
    Teuchos::RCP<Epetra_Vector> vector(int i) const
    {
      return Teuchos::rcp(new Epetra_Vector(*Map(i)));
    }

    /// create multi vector to map i
    Teuchos::RCP<Epetra_MultiVector> vector(int i, int numvec) const
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
    Teuchos::RCP<Epetra_Vector> extract_vector(const Epetra_Vector& full, int block) const;

    /// extract a partial vector from a full vector
    /*!
      \param full vector on the full map
      \param block number of vector to extract
     */
    Teuchos::RCP<Epetra_MultiVector> extract_vector(
        const Epetra_MultiVector& full, int block) const;

    /// extract a partial vector from a full vector
    /*!
      \param full vector on the full map
      \param block number of vector to extract
     */
    Teuchos::RCP<Epetra_Vector> extract_vector(Teuchos::RCP<Epetra_Vector> full, int block) const
    {
      return extract_vector(*full, block);
    }

    /// extract a partial vector from a full vector
    /*!
      \param full vector on the full map
      \param block number of vector to extract
     */
    Teuchos::RCP<Epetra_MultiVector> extract_vector(
        Teuchos::RCP<Epetra_MultiVector> full, int block) const
    {
      return extract_vector(*full, block);
    }

    /// extract a partial vector from a full vector
    /*!
      \param full vector on the full map
      \param block number of vector to extract
     */
    Teuchos::RCP<Epetra_Vector> extract_vector(
        Teuchos::RCP<const Epetra_Vector> full, int block) const
    {
      return extract_vector(*full, block);
    }

    /// extract a partial vector from a full vector
    /*!
      \param full vector on the full map
      \param block number of vector to extract
     */
    Teuchos::RCP<Epetra_MultiVector> extract_vector(
        Teuchos::RCP<const Epetra_MultiVector> full, int block) const
    {
      return extract_vector(*full, block);
    }

    /// extract a partial vector from a full vector
    /*!
      \param full vector on the full map
      \param block number of vector to extract
      \param partial vector to fill
     */
    virtual void extract_vector(
        const Epetra_MultiVector& full, int block, Epetra_MultiVector& partial) const;

    /// extract a partial vector from a full vector
    /*!
      \param full vector on the full map
      \param block number of vector to extract
      \param partial vector to fill
     */
    void extract_vector(Teuchos::RCP<const Epetra_Vector> full, int block,
        Teuchos::RCP<Epetra_Vector> partial) const
    {
      extract_vector(*full, block, *partial);
    }

    //@}

    /** \name Insert from full dof vector */
    //@{

    /// Put a partial vector into a full vector
    /*!
      \param partial vector to copy into full vector
      \param block number of partial vector
     */
    Teuchos::RCP<Epetra_Vector> insert_vector(const Epetra_Vector& partial, int block) const;

    /// Put a partial vector into a full vector
    /*!
      \param partial vector to copy into full vector
      \param block number of partial vector
     */
    Teuchos::RCP<Epetra_MultiVector> insert_vector(
        const Epetra_MultiVector& partial, int block) const;

    /// Put a partial vector into a full vector
    /*!
      \param partial vector to copy into full vector
      \param block number of partial vector
     */
    Teuchos::RCP<Epetra_Vector> insert_vector(
        Teuchos::RCP<const Epetra_Vector> partial, int block) const
    {
      return insert_vector(*partial, block);
    }

    /// Put a partial vector into a full vector
    /*!
      \param partial vector to copy into full vector
      \param block number of partial vector
     */
    Teuchos::RCP<Epetra_MultiVector> insert_vector(
        Teuchos::RCP<const Epetra_MultiVector> partial, int block) const
    {
      return insert_vector(*partial, block);
    }

    /// Put a partial vector into a full vector
    /*!
      \param partial vector to copy into full vector
      \param block number of partial vector
     */
    Teuchos::RCP<Epetra_Vector> insert_vector(Teuchos::RCP<Epetra_Vector> partial, int block) const
    {
      return insert_vector(*partial, block);
    }

    /// Put a partial vector into a full vector
    /*!
      \param partial vector to copy into full vector
      \param block number of partial vector
     */
    Teuchos::RCP<Epetra_MultiVector> insert_vector(
        Teuchos::RCP<Epetra_MultiVector> partial, int block) const
    {
      return insert_vector(*partial, block);
    }

    /// Put a partial vector into a full vector
    /*!
      \param partial vector to copy into full vector
      \param block number of partial vector
      \param full vector to copy into
     */
    virtual void insert_vector(
        const Epetra_MultiVector& partial, int block, Epetra_MultiVector& full) const;

    /// Put a partial vector into a full vector
    /*!
      \param partial vector to copy into full vector
      \param block number of partial vector
      \param full vector to copy into
     */
    void insert_vector(Teuchos::RCP<const Epetra_Vector> partial, int block,
        Teuchos::RCP<Epetra_Vector> full) const
    {
      insert_vector(*partial, block, *full);
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
    virtual void add_vector(const Epetra_MultiVector& partial, int block, Epetra_MultiVector& full,
        double scale = 1.0) const;

    /// Put a partial vector into a full vector
    /*!
      \param partial vector to copy into full vector
      \param block number of partial vector
      \param full vector to copy into
      \param scale scaling factor for partial vector
     */
    void add_vector(Teuchos::RCP<const Epetra_Vector> partial, int block,
        Teuchos::RCP<Epetra_Vector> full, double scale = 1.0) const
    {
      add_vector(*partial, block, *full, scale);
    }

    //@}

    /// PutScalar to one block only
    void put_scalar(Epetra_Vector& full, int block, double scalar) const;

    /// L2-norm of one block only
    double norm2(const Epetra_Vector& full, int block) const;

    /// Scale one block only
    void scale(Epetra_Vector& full, int block, double scalar) const;

    /// Scale one block only
    void scale(Epetra_MultiVector& full, int block, double scalar) const;

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
  Teuchos::RCP<Epetra_Vector> extract_##name##_vector(const Epetra_Vector& full) const             \
  {                                                                                                \
    return MultiMapExtractor::extract_vector(full, pos);                                           \
  }                                                                                                \
                                                                                                   \
  Teuchos::RCP<Epetra_Vector> extract_##name##_vector(Teuchos::RCP<const Epetra_Vector> full)      \
      const                                                                                        \
  {                                                                                                \
    return MultiMapExtractor::extract_vector(full, pos);                                           \
  }                                                                                                \
                                                                                                   \
  void extract_##name##_vector(                                                                    \
      Teuchos::RCP<const Epetra_Vector> full, Teuchos::RCP<Epetra_Vector> cond) const              \
  {                                                                                                \
    extract_vector(full, pos, cond);                                                               \
  }                                                                                                \
                                                                                                   \
  Teuchos::RCP<Epetra_Vector> insert_##name##_vector(Teuchos::RCP<const Epetra_Vector> cond) const \
  {                                                                                                \
    return insert_vector(cond, pos);                                                               \
  }                                                                                                \
                                                                                                   \
  void insert_##name##_vector(                                                                     \
      Teuchos::RCP<const Epetra_Vector> cond, Teuchos::RCP<Epetra_Vector> full) const              \
  {                                                                                                \
    insert_vector(cond, pos, full);                                                                \
  }                                                                                                \
                                                                                                   \
  void add_##name##_vector(                                                                        \
      Teuchos::RCP<const Epetra_Vector> cond, Teuchos::RCP<Epetra_Vector> full) const              \
  {                                                                                                \
    add_vector(cond, pos, full);                                                                   \
  }                                                                                                \
                                                                                                   \
  void add_##name##_vector(double scale, Teuchos::RCP<const Epetra_Vector> cond,                   \
      Teuchos::RCP<Epetra_Vector> full) const                                                      \
  {                                                                                                \
    add_vector(cond, pos, full, scale);                                                            \
  }                                                                                                \
                                                                                                   \
  const Teuchos::RCP<const Epetra_Map>& name##_map() const { return Map(pos); }                    \
                                                                                                   \
  bool name##_relevant() const { return name##_map()->NumGlobalElements() != 0; }                  \
                                                                                                   \
  void name##_put_scalar(Epetra_Vector& full, double scalar) const                                 \
  {                                                                                                \
    put_scalar(full, pos, scalar);                                                                 \
  }                                                                                                \
                                                                                                   \
  double name##_norm2(const Epetra_Vector& full) const { return norm2(full, pos); }


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

  At the core there are the cond_map(), the map of all selected dofs, and
  other_map(), the map of all remaining dofs. This duality also exists in
  the extraction methods extract_cond_vector() and extract_other_vector(), that
  extract a subvector from a full one, and the insertion methods
  insert_cond_vector() and insert_other_vector(), that copy a subvector into a
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
     *  You have to call a setup() routine of your choice. */
    MapExtractor();

    /** \brief  constructor
     *
     *  Calls setup() from known maps */
    MapExtractor(const Epetra_Map& fullmap, Teuchos::RCP<const Epetra_Map> condmap,
        Teuchos::RCP<const Epetra_Map> othermap);

    /** \brief constructor
     *
     *  Calls setup() to create non-overlapping othermap/condmap which is complementary
     *  to condmap/othermap with respect to fullmap depending on boolean 'iscondmap'  */
    MapExtractor(const Epetra_Map& fullmap,         //< full map
        Teuchos::RCP<const Epetra_Map> partialmap,  //< partial map, ie condition or other map
        bool iscondmap = true                       //< true, if partialmap is condition map
    );

    /** \name Setup */
    //@{

    /// setup from known maps
    void setup(const Epetra_Map& fullmap, const Teuchos::RCP<const Epetra_Map>& condmap,
        const Teuchos::RCP<const Epetra_Map>& othermap);

    /// setup creates non-overlapping othermap/condmap which is complementary to condmap/othermap
    /// with respect to fullmap depending on boolean 'iscondmap'
    /// \author bborn
    /// \date 10/08
    void setup(const Epetra_Map& fullmap, const Teuchos::RCP<const Epetra_Map>& partialmap,
        bool iscondmap = true);

    //@}

    MAP_EXTRACTOR_VECTOR_METHODS(cond, 1)
    MAP_EXTRACTOR_VECTOR_METHODS(other, 0)

   private:
  };

}  // namespace Core::LinAlg

FOUR_C_NAMESPACE_CLOSE

#endif
