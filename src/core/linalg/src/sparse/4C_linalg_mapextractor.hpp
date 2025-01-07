// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_LINALG_MAPEXTRACTOR_HPP
#define FOUR_C_LINALG_MAPEXTRACTOR_HPP

#include "4C_config.hpp"

#include "4C_linalg_vector.hpp"

#include <Epetra_Import.h>
#include <Epetra_Map.h>

#include <algorithm>
#include <map>
#include <memory>
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
        const Epetra_Map& fullmap, const std::vector<std::shared_ptr<const Epetra_Map>>& maps);

    /// setup of an empty extractor
    /*!
      \warning The fullmap has to be nonoverlapping. The list of maps has to
      be nonoverlapping as well and its sum has to equal the fullmap.
     */
    void setup(
        const Epetra_Map& fullmap, const std::vector<std::shared_ptr<const Epetra_Map>>& maps);

    /// debug helper
    /*!
      loop all maps in the list of nonoverlapping partial row maps unequal
      nullptr and check whether they have a valid DataPtr() and are
      unique

      \note hidden calls to redistribute() may render maps in maps_ obsolete.
      This function is intended to simplify debugging for these cases.
    */
    void check_for_valid_map_extractor() const;

    /// merge set of unique maps
    /*!
      \warning There must be no overlap in these maps.
               The order of the GIDs is destroyed
     */
    static std::shared_ptr<Epetra_Map> merge_maps(
        const std::vector<std::shared_ptr<const Epetra_Map>>& maps);

    /// merge set of unique maps
    /*!
      \warning There must be no overlap in these maps.
    */
    static std::shared_ptr<Epetra_Map> merge_maps_keep_order(
        const std::vector<std::shared_ptr<const Epetra_Map>>& maps);

    /// intersect set of unique maps
    /*!
      \warning There must be no overlap in these maps.
     */
    static std::shared_ptr<Epetra_Map> intersect_maps(
        const std::vector<std::shared_ptr<const Epetra_Map>>& maps);

    /** \name Maps */
    //@{

    /// number of partial maps
    int num_maps() const { return maps_.size(); }

    /// get the map
    const std::shared_ptr<const Epetra_Map>& Map(int i) const { return maps_[i]; }

    /// the full map
    const std::shared_ptr<const Epetra_Map>& full_map() const { return fullmap_; }

    //@}

    /** \name Vector creation */
    //@{

    /// create vector to map i
    std::shared_ptr<Core::LinAlg::Vector<double>> vector(int i) const
    {
      return std::make_shared<Core::LinAlg::Vector<double>>(*Map(i));
    }

    /// create multi vector to map i
    std::shared_ptr<Core::LinAlg::MultiVector<double>> vector(int i, int numvec) const
    {
      return std::make_shared<Core::LinAlg::MultiVector<double>>(*Map(i), numvec);
    }

    //@

    /** \name Extract from full vector */
    //@{

    /// extract a partial vector from a full vector
    /*!
      \param full vector on the full map
      \param block number of vector to extract
     */
    std::shared_ptr<Core::LinAlg::Vector<double>> extract_vector(
        const Core::LinAlg::Vector<double>& full, int block) const;

    /// extract a partial vector from a full vector
    /*!
      \param full vector on the full map
      \param block number of vector to extract
     */
    std::shared_ptr<Core::LinAlg::MultiVector<double>> extract_vector(
        const Core::LinAlg::MultiVector<double>& full, int block) const;

    /// extract a partial vector from a full vector
    /*!
      \param full vector on the full map
      \param block number of vector to extract
      \param partial vector to fill
     */
    virtual void extract_vector(const Core::LinAlg::MultiVector<double>& full, int block,
        Core::LinAlg::MultiVector<double>& partial) const;


    //@}

    /** \name Insert from full dof vector */
    //@{

    /// Put a partial vector into a full vector
    /*!
      \param partial vector to copy into full vector
      \param block number of partial vector
     */
    std::shared_ptr<Core::LinAlg::Vector<double>> insert_vector(
        const Core::LinAlg::Vector<double>& partial, int block) const;

    /// Put a partial vector into a full vector
    /*!
      \param partial vector to copy into full vector
      \param block number of partial vector
     */
    std::shared_ptr<Core::LinAlg::MultiVector<double>> insert_vector(
        const Core::LinAlg::MultiVector<double>& partial, int block) const;

    /// Put a partial vector into a full vector
    /*!
      \param partial vector to copy into full vector
      \param block number of partial vector
      \param full vector to copy into
     */
    virtual void insert_vector(const Core::LinAlg::MultiVector<double>& partial, int block,
        Core::LinAlg::MultiVector<double>& full) const;

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
    virtual void add_vector(const Core::LinAlg::MultiVector<double>& partial, int block,
        Core::LinAlg::MultiVector<double>& full, double scale = 1.0) const;

    //@}

    /// PutScalar to one block only
    void put_scalar(Core::LinAlg::Vector<double>& full, int block, double scalar) const;

    /// L2-norm of one block only
    double norm2(const Core::LinAlg::Vector<double>& full, int block) const;

    /// Scale one block only
    void scale(Core::LinAlg::Vector<double>& full, int block, double scalar) const;

    /// Scale one block only
    void scale(Core::LinAlg::MultiVector<double>& full, int block, double scalar) const;

   protected:
    /// the full row map
    std::shared_ptr<const Epetra_Map> fullmap_;

    /// the list of nonoverlapping partial row maps that sums up to the full map
    std::vector<std::shared_ptr<const Epetra_Map>> maps_;

    /// communication between condition dof map and full row dof map
    std::vector<std::shared_ptr<Epetra_Import>> importer_;
  };


/// Add all kinds of support methods to derived classes of MultiMapExtractor.
#define MAP_EXTRACTOR_VECTOR_METHODS(name, pos)                                           \
  std::shared_ptr<Core::LinAlg::Vector<double>> extract_##name##_vector(                  \
      const Core::LinAlg::Vector<double>& full) const                                     \
  {                                                                                       \
    return MultiMapExtractor::extract_vector(full, pos);                                  \
  }                                                                                       \
                                                                                          \
  void extract_##name##_vector(                                                           \
      const Core::LinAlg::Vector<double>& full, Core::LinAlg::Vector<double>& cond) const \
  {                                                                                       \
    extract_vector(full, pos, cond);                                                      \
  }                                                                                       \
                                                                                          \
  std::shared_ptr<Core::LinAlg::Vector<double>> insert_##name##_vector(                   \
      const Core::LinAlg::Vector<double>& cond) const                                     \
  {                                                                                       \
    return insert_vector(cond, pos);                                                      \
  }                                                                                       \
                                                                                          \
  void insert_##name##_vector(                                                            \
      const Core::LinAlg::Vector<double>& cond, Core::LinAlg::Vector<double>& full) const \
  {                                                                                       \
    insert_vector(cond, pos, full);                                                       \
  }                                                                                       \
                                                                                          \
  void add_##name##_vector(                                                               \
      const Core::LinAlg::Vector<double>& cond, Core::LinAlg::Vector<double>& full) const \
  {                                                                                       \
    add_vector(cond, pos, full);                                                          \
  }                                                                                       \
                                                                                          \
  void add_##name##_vector(double scale, const Core::LinAlg::Vector<double>& cond,        \
      Core::LinAlg::Vector<double>& full) const                                           \
  {                                                                                       \
    add_vector(cond, pos, full, scale);                                                   \
  }                                                                                       \
                                                                                          \
  const std::shared_ptr<const Epetra_Map>& name##_map() const { return Map(pos); }        \
                                                                                          \
  bool name##_relevant() const { return name##_map()->NumGlobalElements() != 0; }         \
                                                                                          \
  void name##_put_scalar(Core::LinAlg::Vector<double>& full, double scalar) const         \
  {                                                                                       \
    put_scalar(full, pos, scalar);                                                        \
  }                                                                                       \
                                                                                          \
  double name##_norm2(const Core::LinAlg::Vector<double>& full) const { return norm2(full, pos); }


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
    MapExtractor(const Epetra_Map& fullmap, std::shared_ptr<const Epetra_Map> condmap,
        std::shared_ptr<const Epetra_Map> othermap);

    /** \brief constructor
     *
     *  Calls setup() to create non-overlapping othermap/condmap which is complementary
     *  to condmap/othermap with respect to fullmap depending on boolean 'iscondmap'  */
    MapExtractor(const Epetra_Map& fullmap,            //< full map
        std::shared_ptr<const Epetra_Map> partialmap,  //< partial map, ie condition or other map
        bool iscondmap = true                          //< true, if partialmap is condition map
    );

    /** \name Setup */
    //@{

    /// setup from known maps
    void setup(const Epetra_Map& fullmap, const std::shared_ptr<const Epetra_Map>& condmap,
        const std::shared_ptr<const Epetra_Map>& othermap);

    /// setup creates non-overlapping othermap/condmap which is complementary to condmap/othermap
    /// with respect to fullmap depending on boolean 'iscondmap'
    /// \author bborn
    /// \date 10/08
    void setup(const Epetra_Map& fullmap, const std::shared_ptr<const Epetra_Map>& partialmap,
        bool iscondmap = true);

    //@}

    MAP_EXTRACTOR_VECTOR_METHODS(cond, 1)
    MAP_EXTRACTOR_VECTOR_METHODS(other, 0)

   private:
  };

}  // namespace Core::LinAlg

FOUR_C_NAMESPACE_CLOSE

#endif
