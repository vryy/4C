// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_STRUCTURE_AUX_HPP
#define FOUR_C_STRUCTURE_AUX_HPP

#include "4C_config.hpp"

#include "4C_linalg_mapextractor.hpp"
#include "4C_structure_new_input.hpp"

FOUR_C_NAMESPACE_OPEN

// forward declarations
namespace Core::FE
{
  class Discretization;
}  // namespace Core::FE

namespace Solid
{
  /// Determine norm of force residual
  double calculate_vector_norm(const Solid::VectorNorm norm,  ///< type of norm to use
      const Core::LinAlg::Vector<double>& vect,               ///< the vector of interest
      const int numneglect =
          0  ///< number of DOFs that have to be neglected for possible length scaling
  );

  /// specific MultiMapExtractor to handle the structure field
  class MapExtractor : public Core::LinAlg::MultiMapExtractor
  {
   public:
    enum
    {
      cond_other = 0,
      cond_fsi = 1,
      cond_lung_asi = 2,
      cond_bio_gr = 3,
      cond_ale_wear = 4,
      cond_fpsi = 5,
      cond_immersed = 6,
      cond_pasi = 7
    };

    /// setup the whole thing
    void setup(const Core::FE::Discretization& dis, const Core::LinAlg::Map& fullmap,
        bool overlapping = false);

    /// get all element gids those nodes are touched by any condition
    std::shared_ptr<std::set<int>> conditioned_element_map(
        const Core::FE::Discretization& dis) const;

    std::shared_ptr<Core::LinAlg::Vector<double>> extract_other_vector(
        const Core::LinAlg::Vector<double>& full) const
    {
      return MultiMapExtractor::extract_vector(full, cond_other);
    }
    void extract_other_vector(
        const Core::LinAlg::Vector<double>& full, Core::LinAlg::Vector<double>& cond) const
    {
      extract_vector(full, cond_other, cond);
    }
    std::shared_ptr<Core::LinAlg::Vector<double>> insert_other_vector(
        const Core::LinAlg::Vector<double>& cond) const
    {
      return insert_vector(cond, cond_other);
    }
    void insert_other_vector(
        const Core::LinAlg::Vector<double>& cond, Core::LinAlg::Vector<double>& full) const
    {
      insert_vector(cond, cond_other, full);
    }
    void add_other_vector(
        const Core::LinAlg::Vector<double>& cond, Core::LinAlg::Vector<double>& full) const
    {
      add_vector(cond, cond_other, full);
    }
    void add_other_vector(double scale, const Core::LinAlg::Vector<double>& cond,
        Core::LinAlg::Vector<double>& full) const
    {
      add_vector(cond, cond_other, full, scale);
    }
    const std::shared_ptr<const Core::LinAlg::Map>& other_map() const { return map(cond_other); }
    bool other_relevant() const { return other_map()->num_global_elements() != 0; }
    void other_put_scalar(Core::LinAlg::Vector<double>& full, double scalar) const
    {
      put_scalar(full, cond_other, scalar);
    }
    std::shared_ptr<Core::LinAlg::Vector<double>> extract_fsi_cond_vector(
        const Core::LinAlg::Vector<double>& full) const
    {
      return MultiMapExtractor::extract_vector(full, cond_fsi);
    }
    void extract_fsi_cond_vector(
        const Core::LinAlg::Vector<double>& full, Core::LinAlg::Vector<double>& cond) const
    {
      extract_vector(full, cond_fsi, cond);
    }
    std::shared_ptr<Core::LinAlg::Vector<double>> insert_fsi_cond_vector(
        const Core::LinAlg::Vector<double>& cond) const
    {
      return insert_vector(cond, cond_fsi);
    }
    void insert_fsi_cond_vector(
        const Core::LinAlg::Vector<double>& cond, Core::LinAlg::Vector<double>& full) const
    {
      insert_vector(cond, cond_fsi, full);
    }
    void add_fsi_cond_vector(
        const Core::LinAlg::Vector<double>& cond, Core::LinAlg::Vector<double>& full) const
    {
      add_vector(cond, cond_fsi, full);
    }
    void add_fsi_cond_vector(double scale, const Core::LinAlg::Vector<double>& cond,
        Core::LinAlg::Vector<double>& full) const
    {
      add_vector(cond, cond_fsi, full, scale);
    }
    const std::shared_ptr<const Core::LinAlg::Map>& fsi_cond_map() const { return map(cond_fsi); }
    bool fsi_cond_relevant() const { return fsi_cond_map()->num_global_elements() != 0; }
    void fsi_cond_put_scalar(Core::LinAlg::Vector<double>& full, double scalar) const
    {
      put_scalar(full, cond_fsi, scalar);
    }
    std::shared_ptr<Core::LinAlg::Vector<double>> extract_lung_asi_cond_vector(
        const Core::LinAlg::Vector<double>& full) const
    {
      return MultiMapExtractor::extract_vector(full, cond_lung_asi);
    }
    void extract_lung_asi_cond_vector(
        const Core::LinAlg::Vector<double>& full, Core::LinAlg::Vector<double>& cond) const
    {
      extract_vector(full, cond_lung_asi, cond);
    }
    std::shared_ptr<Core::LinAlg::Vector<double>> insert_lung_asi_cond_vector(
        const Core::LinAlg::Vector<double>& cond) const
    {
      return insert_vector(cond, cond_lung_asi);
    }
    void insert_lung_asi_cond_vector(
        const Core::LinAlg::Vector<double>& cond, Core::LinAlg::Vector<double>& full) const
    {
      insert_vector(cond, cond_lung_asi, full);
    }
    void add_lung_asi_cond_vector(
        const Core::LinAlg::Vector<double>& cond, Core::LinAlg::Vector<double>& full) const
    {
      add_vector(cond, cond_lung_asi, full);
    }
    void add_lung_asi_cond_vector(double scale, const Core::LinAlg::Vector<double>& cond,
        Core::LinAlg::Vector<double>& full) const
    {
      add_vector(cond, cond_lung_asi, full, scale);
    }
    const std::shared_ptr<const Core::LinAlg::Map>& lung_asi_cond_map() const
    {
      return map(cond_lung_asi);
    }
    bool lung_asi_cond_relevant() const { return lung_asi_cond_map()->num_global_elements() != 0; }
    void lung_asi_cond_put_scalar(Core::LinAlg::Vector<double>& full, double scalar) const
    {
      put_scalar(full, cond_lung_asi, scalar);
    }
    std::shared_ptr<Core::LinAlg::Vector<double>> extract_ale_wear_cond_vector(
        const Core::LinAlg::Vector<double>& full) const
    {
      return MultiMapExtractor::extract_vector(full, cond_ale_wear);
    }
    void extract_ale_wear_cond_vector(
        const Core::LinAlg::Vector<double>& full, Core::LinAlg::Vector<double>& cond) const
    {
      extract_vector(full, cond_ale_wear, cond);
    }
    std::shared_ptr<Core::LinAlg::Vector<double>> insert_ale_wear_cond_vector(
        const Core::LinAlg::Vector<double>& cond) const
    {
      return insert_vector(cond, cond_ale_wear);
    }
    void insert_ale_wear_cond_vector(
        const Core::LinAlg::Vector<double>& cond, Core::LinAlg::Vector<double>& full) const
    {
      insert_vector(cond, cond_ale_wear, full);
    }
    void add_ale_wear_cond_vector(
        const Core::LinAlg::Vector<double>& cond, Core::LinAlg::Vector<double>& full) const
    {
      add_vector(cond, cond_ale_wear, full);
    }
    void add_ale_wear_cond_vector(double scale, const Core::LinAlg::Vector<double>& cond,
        Core::LinAlg::Vector<double>& full) const
    {
      add_vector(cond, cond_ale_wear, full, scale);
    }
    const std::shared_ptr<const Core::LinAlg::Map>& ale_wear_cond_map() const
    {
      return map(cond_ale_wear);
    }
    bool ale_wear_cond_relevant() const { return ale_wear_cond_map()->num_global_elements() != 0; }
    void ale_wear_cond_put_scalar(Core::LinAlg::Vector<double>& full, double scalar) const
    {
      put_scalar(full, cond_ale_wear, scalar);
    }
    std::shared_ptr<Core::LinAlg::Vector<double>> extract_fpsi_cond_vector(
        const Core::LinAlg::Vector<double>& full) const
    {
      return MultiMapExtractor::extract_vector(full, cond_fpsi);
    }
    void extract_fpsi_cond_vector(
        const Core::LinAlg::Vector<double>& full, Core::LinAlg::Vector<double>& cond) const
    {
      extract_vector(full, cond_fpsi, cond);
    }
    std::shared_ptr<Core::LinAlg::Vector<double>> insert_fpsi_cond_vector(
        const Core::LinAlg::Vector<double>& cond) const
    {
      return insert_vector(cond, cond_fpsi);
    }
    void insert_fpsi_cond_vector(
        const Core::LinAlg::Vector<double>& cond, Core::LinAlg::Vector<double>& full) const
    {
      insert_vector(cond, cond_fpsi, full);
    }
    void add_fpsi_cond_vector(
        const Core::LinAlg::Vector<double>& cond, Core::LinAlg::Vector<double>& full) const
    {
      add_vector(cond, cond_fpsi, full);
    }
    void add_fpsi_cond_vector(double scale, const Core::LinAlg::Vector<double>& cond,
        Core::LinAlg::Vector<double>& full) const
    {
      add_vector(cond, cond_fpsi, full, scale);
    }
    const std::shared_ptr<const Core::LinAlg::Map>& fpsi_cond_map() const { return map(cond_fpsi); }
    bool fpsi_cond_relevant() const { return fpsi_cond_map()->num_global_elements() != 0; }
    void fpsi_cond_put_scalar(Core::LinAlg::Vector<double>& full, double scalar) const
    {
      put_scalar(full, cond_fpsi, scalar);
    }
    std::shared_ptr<Core::LinAlg::Vector<double>> extract_immersed_cond_vector(
        const Core::LinAlg::Vector<double>& full) const
    {
      return MultiMapExtractor::extract_vector(full, cond_immersed);
    }
    void extract_immersed_cond_vector(
        const Core::LinAlg::Vector<double>& full, Core::LinAlg::Vector<double>& cond) const
    {
      extract_vector(full, cond_immersed, cond);
    }
    std::shared_ptr<Core::LinAlg::Vector<double>> insert_immersed_cond_vector(
        const Core::LinAlg::Vector<double>& cond) const
    {
      return insert_vector(cond, cond_immersed);
    }
    void insert_immersed_cond_vector(
        const Core::LinAlg::Vector<double>& cond, Core::LinAlg::Vector<double>& full) const
    {
      insert_vector(cond, cond_immersed, full);
    }
    void add_immersed_cond_vector(
        const Core::LinAlg::Vector<double>& cond, Core::LinAlg::Vector<double>& full) const
    {
      add_vector(cond, cond_immersed, full);
    }
    void add_immersed_cond_vector(double scale, const Core::LinAlg::Vector<double>& cond,
        Core::LinAlg::Vector<double>& full) const
    {
      add_vector(cond, cond_immersed, full, scale);
    }
    const std::shared_ptr<const Core::LinAlg::Map>& immersed_cond_map() const
    {
      return map(cond_immersed);
    }
    bool immersed_cond_relevant() const { return immersed_cond_map()->num_global_elements() != 0; }
    void immersed_cond_put_scalar(Core::LinAlg::Vector<double>& full, double scalar) const
    {
      put_scalar(full, cond_immersed, scalar);
    }
    std::shared_ptr<Core::LinAlg::Vector<double>> extract_pasi_cond_vector(
        const Core::LinAlg::Vector<double>& full) const
    {
      return MultiMapExtractor::extract_vector(full, cond_pasi);
    }
    void extract_pasi_cond_vector(
        const Core::LinAlg::Vector<double>& full, Core::LinAlg::Vector<double>& cond) const
    {
      extract_vector(full, cond_pasi, cond);
    }
    std::shared_ptr<Core::LinAlg::Vector<double>> insert_pasi_cond_vector(
        const Core::LinAlg::Vector<double>& cond) const
    {
      return insert_vector(cond, cond_pasi);
    }
    void insert_pasi_cond_vector(
        const Core::LinAlg::Vector<double>& cond, Core::LinAlg::Vector<double>& full) const
    {
      insert_vector(cond, cond_pasi, full);
    }
    void add_pasi_cond_vector(
        const Core::LinAlg::Vector<double>& cond, Core::LinAlg::Vector<double>& full) const
    {
      add_vector(cond, cond_pasi, full);
    }
    void add_pasi_cond_vector(double scale, const Core::LinAlg::Vector<double>& cond,
        Core::LinAlg::Vector<double>& full) const
    {
      add_vector(cond, cond_pasi, full, scale);
    }
    const std::shared_ptr<const Core::LinAlg::Map>& pasi_cond_map() const { return map(cond_pasi); }
    bool pasi_cond_relevant() const { return pasi_cond_map()->num_global_elements() != 0; }
    void pasi_cond_put_scalar(Core::LinAlg::Vector<double>& full, double scalar) const
    {
      put_scalar(full, cond_pasi, scalar);
    }
  };
}  // namespace Solid

FOUR_C_NAMESPACE_CLOSE

#endif
