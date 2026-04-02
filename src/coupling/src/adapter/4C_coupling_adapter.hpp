// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_COUPLING_ADAPTER_HPP
#define FOUR_C_COUPLING_ADAPTER_HPP

#include "4C_config.hpp"

#include "4C_coupling_adapter_base.hpp"
#include "4C_linalg_fevector.hpp"
#include "4C_linalg_map.hpp"

#include <map>
#include <memory>

FOUR_C_NAMESPACE_OPEN

namespace Core::LinAlg
{
  class SparseMatrix;

}  // namespace Core::LinAlg

namespace Core::FE
{
  class Discretization;
}  // namespace Core::FE

namespace Coupling::Adapter
{
  /*! \class Coupling
   *  \brief Management of coupling between two (matching) sets of nodes from
   *  two discretizations
   *
   *  Interface coupled problems with matching meshes need to transfer
   *  results between their interfaces. That is values which belong to
   *  the dofs of one side have to be accessed by the other side as
   *  well. In a parallel setting, of course, these sides do not in
   *  general reside on the same processors.
   *
   *  This class handles the transfer of dof values at the interface
   *  (Core::LinAlg::Vectors build on the interface dof map of either side) to
   *  the other side. To distinguish both sides lets speak of target and
   *  source, even though no side really dominates the other. On the
   *  contrary we provide perfect symmetry once the setup is done.
   *
   *  The idea is simple: We have a target dof map that describes the
   *  distribution of the target's interface dofs. And we have a source
   *  dof map that describes the distribution of the source's interface
   *  dofs. Both maps, however, live in their respective
   *  communicator. We cannot transfer data between them. Furthermore
   *  the dof gids are totally independent of each other.
   *
   *  So we build a mapping between the nodes from both meshes during
   *  the setup phase. This establishes the connection. From that
   *  connection we deduce (in a somewhat painful way) a permuted
   *  target dof map and a permuted source dof map. These permuted
   *  maps are bound to have the same layout as the normal maps from the
   *  other side. So we can exchange dof values between fields by simply
   *  copying from a normal Core::LinAlg::Vector<double> to the permuted
   *  Core::LinAlg::Vector<double> from the other side without actually looking at the
   *  respective maps. Afterwards the communication happens within one
   *  field in the usual fashion. So the transfer functions
   *  target_to_source() and source_to_target() are quite simple. The hard
   *  work happens (once) during setup.
   *

   */
  class Coupling : public CouplingBase
  {
   public:
    /// empty constructor
    Coupling();

    /** \name Setup */
    //@{

    /// setup coupling of nodes marked with condition
    /*!
      Setup the whole thing. Find matching nodes via octtree and build
      appropriate dof maps.

      first variant will couple the first numdof DOFs of both discretizations

      second variant accepts two vectors to couple the DOFs specified in these
      vectors, e.g. masterdofs = [0 3 5] and slavedofs = [1 2 4] will couple
      DOF 0 (target) with DOF 1 source, DOF 3 (target) with DOF 2 source, etc.

     */
    void setup_condition_coupling(const Core::FE::Discretization& masterdis,
        std::shared_ptr<const Core::LinAlg::Map> mastercondmap,
        const Core::FE::Discretization& slavedis,
        std::shared_ptr<const Core::LinAlg::Map> slavecondmap, const std::string& condname,
        const int numdof, bool matchall = true, const int nds_master = 0, const int nds_slave = 0);

    void setup_condition_coupling(const Core::FE::Discretization& masterdis,
        std::shared_ptr<const Core::LinAlg::Map> mastercondmap,
        const Core::FE::Discretization& slavedis,
        std::shared_ptr<const Core::LinAlg::Map> slavecondmap, const std::string& condname,
        const std::vector<int>& masterdofs, const std::vector<int>& slavedofs, bool matchall = true,
        const int nds_master = 0, const int nds_slave = 0);

    /*! \brief Setup coupling of given nodes
     *
     *  Setup the whole thing. Find matching nodes via octtree and build
     *  appropriate dof maps.
     *
     *  \note All source nodes are required to find a match. In contrast
     *  target nodes do not need to have a match. So it is legal to hand
     *  in more target nodes that source nodes as long as there is one
     *  target node for each source node.
     *
     *  We need some way to guess the tolerance for the octree. It must not be
     *  too small, otherwise we won't find matching nodes. Too large a tolerance
     *  will not hurt that much. It just means we will have to test more nodes.
     *  (But it could hurt, if the geometry is very small. It should be
     *   around the size of the smallest element)
     *
     *  \param masterdis   (i) target side mesh
     *  \param slavedis    (i) source side mesh
     *  \param masternodes (i) gids of nodes on target side to be coupled
     *  \param slavenodes  (i) gids of nodes on source side to be coupled
     *  \param masterdofs  (i) vector with DOF ids of target side to be coupled
     *  \param slavedofs   (i) vector with DOF ids of side side to be coupled
     *  \param matchall    (i) flag indicating matching source and target nodes
     *  \param tolerance   (i) tolerance for octree for node matching
     */
    void setup_coupling(const Core::FE::Discretization& masterdis,
        const Core::FE::Discretization& slavedis, const std::vector<int>& masternodes,
        const std::vector<int>& slavenodes, const std::vector<int>& masterdofs,
        const std::vector<int>& slavedofs, const bool matchall = true,
        const double tolerance = 1.e-3, const int nds_master = 0, const int nds_slave = 0);

    /*! \brief Setup coupling of given nodes
     *
     *  only difference to the variant above is that the first numdof
     *  DOFs on both source and target side will be coupled
     *
     *  \param numdof      (i) number of dofs per node to be coupled
     */
    void setup_coupling(const Core::FE::Discretization& masterdis,
        const Core::FE::Discretization& slavedis, const std::vector<int>& masternodes,
        const std::vector<int>& slavenodes, const int numdof, const bool matchall = true,
        const double tolerance = 1.e-3, const int nds_master = 0, const int nds_slave = 0);

    /*! \brief Setup coupling of given nodes
     *
     *  Setup the whole thing. Find matching nodes via octtree and build
     *  appropriate dof maps.
     *
     *  \note All source nodes are required to find a match. In contrast
     *  target nodes do not need to have a match. So it is legal to hand
     *  in more target nodes that source nodes as long as there is one
     *  target node for each source node.
     *
     *  We need some way to guess the tolerance for the octree. It must not be
     *  too small, otherwise we won't find matching nodes. Too large a tolerance
     *  will not hurt that much. It just means we will have to test more nodes.
     *
     *  \param masterdis   (i) target side mesh
     *  \param slavedis    (i) source side mesh
     *  \param masternodes (i) gids of nodes on target side to be coupled
     *  \param slavenodes  (i) gids of nodes on source side to be coupled
     *  \param numdof      (i) number of dofs per node to be coupled
     *  \param matchall    (i) flag indicating matching source and target nodes
     *  \param tolerance   (i) tolerance for octree for node matching
     */
    void setup_coupling(const Core::FE::Discretization& masterdis,
        const Core::FE::Discretization& slavedis, const Core::LinAlg::Map& masternodes,
        const Core::LinAlg::Map& slavenodes, const int numdof, const bool matchall = true,
        const double tolerance = 1.e-3, const int nds_master = 0, const int nds_slave = 0);

    /*! \brief Setup coupling of given nodes and node maps
     *
     *  Setup the whole thing. Build appropriate dof maps (for given, matching nodes on both
     * sides).
     *
     *  \note All target and source nodes are required to match.
     *  The permuted source node map is given as input, thus no search tree to find matching nodes
     * is used. Only the corresponding dof maps are build here (+ some checks). This setup method
     * can be used, if one discretization was cloned from the others discretization, and thus
     * matching node GIDs are assured (for Fluid-ALE or poro, for instance). In this case, this
     * setup gives the same results as the standard version above, but without the need to perform
     * the matching node search.
     *
     *  \param masterdis         (i) target side mesh
     *  \param slavedis          (i) source side mesh
     *  \param masternodemap     (i) gids of nodes on target side to be coupled
     *  \param slavenodemap      (i) gids of nodes on source side to be coupled
     *  \param permslavenodemap  (i) permuted node map of the source side
     *  \param numdof            (i) number of dofs per node to be coupled
     */
    void setup_coupling(const Core::FE::Discretization& masterdis,
        const Core::FE::Discretization& slavedis, const Core::LinAlg::Map& masternodemap,
        const Core::LinAlg::Map& slavenodemap, const Core::LinAlg::Map& permslavenodemap,
        const int numdof);

    /*! \brief Setup coupling of given target and source discretization
     *
     *  Setup the whole thing for matching dofs on both sides
     *
     *  \note All target and source dofs are required to match.
     *  This setup method can be used, if one discretization was cloned from the other
     * discretization, and thus matching dof GIDs are assured. All dofs of all nodes are coupled
     * this way. In this case, this setup gives the same results as the standard version above,
     * but neither a matching node search nor building dof maps is necessary.
     *
     *  \param masterdis         (i) target side mesh
     *  \param slavedis          (i) source side mesh
     */
    void setup_coupling(
        const Core::FE::Discretization& masterdis, const Core::FE::Discretization& slavedis);

    /*! \brief Setup coupling of given nodes
     *
     *  Setup the whole thing. Find matching nodes via octtree and build
     *  appropriate dof maps.
     *
     *  \note All source nodes are required to find a match. In contrast
     *  target nodes do not need to have a match.
     *
     *  Nodes are clustered, such that coupling is performed between clusters
     *  (e.g. masternodes_vec.at(0) <-> slavenodes_vec.at(0), masternodes_vec.at(1) <->
     * slavenodes_vec.at(1), ...)
     *
     *  \param masterdis       (i) target side mesh
     *  \param slavedis        (i) source side mesh
     *  \param masternodes_vec (i) gids of nodes on target side to be coupled
     *  \param slavenodes_vec  (i) gids of nodes on source side to be coupled
     *  \param numdof          (i) number of dofs per node to be coupled
     *  \param matchall        (i) flag indicating matching source and target nodes
     *  \param tolerance       (i) tolerance for octree for node matching
     *  \param nds_master      (i) dofset ID of target side
     *  \param nds_slave       (i) dofset ID of source side
     */
    void setup_coupling(const Core::FE::Discretization& masterdis,
        const Core::FE::Discretization& slavedis,
        const std::vector<std::vector<int>>& masternodes_vec,
        const std::vector<std::vector<int>>& slavenodes_vec, const int numdof,
        const bool matchall = true, const double tolerance = 1.0e-3, const int nds_master = 0,
        const int nds_slave = 0);

    /*! Setup coupling based on dof maps. This can be used if no search is required anymore and
     * all dof maps are known already
     *
     * \param slavedofmap       (i) source side dof map
     * \param permslavedofmap   (i) permuted dofs on source side to match dofs on target side
     * \param masterdofmap      (i) target side dof map
     * \param permmasterdofmap  (i) permuted dofs on target side to match dofs on source side
     */
    void setup_coupling(std::shared_ptr<const Core::LinAlg::Map> slavedofmap,
        std::shared_ptr<const Core::LinAlg::Map> permslavedofmap,
        std::shared_ptr<const Core::LinAlg::Map> masterdofmap,
        std::shared_ptr<const Core::LinAlg::Map> permmasterdofmap);
    //@}

    /** \name Conversion between target and source */
    //@{
    /// There are different versions to satisfy all needs. The basic
    /// idea is the same for all of them.

    /// transfer a dof vector from target to source
    std::shared_ptr<Core::LinAlg::Vector<double>> target_to_source(
        const Core::LinAlg::Vector<double>& mv  ///< target vector (to be transferred)
    ) const override;

    /// transfer a dof vector from source to target
    std::shared_ptr<Core::LinAlg::Vector<double>> source_to_target(
        const Core::LinAlg::Vector<double>& sv  ///< source vector (to be transferred)
    ) const override;

    /// transfer a dof vector from target to source
    std::shared_ptr<Core::LinAlg::FEVector<double>> target_to_source(
        const Core::LinAlg::FEVector<double>& mv  ///< target vector (to be transferred)
    ) const;

    /// transfer a dof vector from source to target
    std::shared_ptr<Core::LinAlg::FEVector<double>> source_to_target(
        const Core::LinAlg::FEVector<double>& sv  ///< source vector (to be transferred)
    ) const;

    /// transfer a dof vector from target to source
    std::shared_ptr<Core::LinAlg::MultiVector<double>> target_to_source(
        const Core::LinAlg::MultiVector<double>& mv  ///< target vector (to be transferred)
    ) const override;

    /// transfer a dof vector from source to target
    std::shared_ptr<Core::LinAlg::MultiVector<double>> source_to_target(
        const Core::LinAlg::MultiVector<double>& sv  ///< source vector (to be transferred)
    ) const override;

    /// transfer a dof vector from target to source
    void target_to_source(
        const Core::LinAlg::MultiVector<double>& mv,  ///< target vector (to be transferred)
        Core::LinAlg::MultiVector<double>& sv         ///< source vector (containing result)
    ) const override;

    /// transfer a dof vector from source to target
    void source_to_target(
        const Core::LinAlg::MultiVector<double>& sv,  ///< source vector (to be transferred)
        Core::LinAlg::MultiVector<double>& mv         ///< target vector (containing result)
    ) const override;

    /// transfer a dof vector from target to source
    void target_to_source(
        const Core::LinAlg::Vector<int>& mv,  ///< target vector (to be transferred)
        Core::LinAlg::Vector<int>& sv         ///< source vector (containing result)
    ) const;

    /// transfer a dof vector from source to target
    void source_to_target(
        const Core::LinAlg::Vector<int>& sv,  ///< source vector (to be transferred)
        Core::LinAlg::Vector<int>& mv         ///< target vector (containing result)
    ) const;

    //@}

    //! @name Access to coupled maps
    //@{

    /// the interface dof map of the target side
    std::shared_ptr<const Core::LinAlg::Map> target_dof_map() const override
    {
      return masterdofmap_;
    }

    /// the interface dof map of the source side
    std::shared_ptr<const Core::LinAlg::Map> source_dof_map() const override
    {
      return slavedofmap_;
    }

    /// the permuted interface dof map of the target side
    std::shared_ptr<const Core::LinAlg::Map> perm_master_dof_map() const
    {
      return permmasterdofmap_;
    }

    /// the permuted interface dof map of the source side
    std::shared_ptr<const Core::LinAlg::Map> perm_source_dof_map() const
    {
      return permslavedofmap_;
    }

    //@}

    //! @name Matrix transform
    //@{

    /// fill rowmap with target -> source pairs
    void fill_master_to_slave_map(std::map<int, int>& rowmap) const;

    /// fill rowmap with source -> target pairs
    void fill_slave_to_master_map(std::map<int, int>& rowmap) const;

    /// fill partial mastermap with gid of partial slavemap
    std::shared_ptr<Core::LinAlg::Map> slave_to_master_map(Core::LinAlg::Map& source);

    /// fill partial slavemap with gid of partial mastermap
    std::shared_ptr<Core::LinAlg::Map> target_to_source_map(Core::LinAlg::Map& target);

    /// redistribute crsmatrix from target row map to permuted target row map
    std::shared_ptr<Core::LinAlg::SparseMatrix> master_to_perm_master(
        const Core::LinAlg::SparseMatrix& sm) const;

    /// redistribute crsmatrix from source row map to permuted source row map
    std::shared_ptr<Core::LinAlg::SparseMatrix> slave_to_perm_slave(
        const Core::LinAlg::SparseMatrix& sm) const;

    //@}

    /// \name Lagrangian coupling helpers

    /// create coupling matrices for Lagrangian coupling conditions
    void setup_coupling_matrices(const Core::LinAlg::Map& shiftedmastermap,
        const Core::LinAlg::Map& masterdomainmap, const Core::LinAlg::Map& slavedomainmap);

    std::shared_ptr<Core::LinAlg::SparseMatrix> master_to_master_mat() const { return matmm_; }
    std::shared_ptr<Core::LinAlg::SparseMatrix> slave_to_master_mat() const { return matsm_; }
    std::shared_ptr<Core::LinAlg::SparseMatrix> master_to_master_mat_trans() const
    {
      return matmm_trans_;
    }
    std::shared_ptr<Core::LinAlg::SparseMatrix> slave_to_master_mat_trans() const
    {
      return matsm_trans_;
    }

    //@}

   protected:
    virtual void build_dof_maps(const Core::FE::Discretization& masterdis,
        const Core::FE::Discretization& slavedis,
        const std::shared_ptr<const Core::LinAlg::Map>& masternodemap,
        const std::shared_ptr<const Core::LinAlg::Map>& slavenodemap,
        const std::shared_ptr<const Core::LinAlg::Map>& permmasternodemap,
        const std::shared_ptr<const Core::LinAlg::Map>& permslavenodemap,
        const std::vector<int>& masterdofs, const std::vector<int>& slavedofs,
        const int nds_master = 0, const int nds_slave = 0);

    /*! \brief Helper function for creating vector from numdof
     *
     *  will return just an ascending vector [0, 1, ..., numdof-1]
     *  if called with numdof > 0, otherwise will return [-1],
     *  which is used inside xfield_field_coupling to determine
     *  if base class is called
     */
    virtual std::vector<int> build_dof_vector_from_num_dof(const int numdof);

   private:
    /*! \brief Do the actual matching of the target <-> source pairs
     *
     *  Here the octree is used. Afterwards all not paired target nodes
     *  are removed.
     *
     *  \param masterdis   (i) target side mesh
     *  \param slavedis    (i) source side mesh
     *  \param masternodes (i/o) all target node gids. on output those that have a match
     *  \param permslavenodes (o) source node gids permuted to match target node gids
     *  \param slavenodes (i) source node gids
     *  \param matchall (i) bool indicating match of all source and target nodes
     *  \param tolerance   (i) tolerance for octree for node matching
     */
    void match_nodes(const Core::FE::Discretization& masterdis,
        const Core::FE::Discretization& slavedis, std::vector<int>& masternodes,
        std::vector<int>& permslavenodes, const std::vector<int>& slavenodes, const bool matchall,
        const double tolerance);

    /// build source to target permutation and dof all maps
    void finish_coupling(const Core::FE::Discretization& masterdis,
        const Core::FE::Discretization& slavedis, std::shared_ptr<Core::LinAlg::Map> masternodemap,
        std::shared_ptr<Core::LinAlg::Map> slavenodemap,
        std::shared_ptr<Core::LinAlg::Map> permslavenodemap, const std::vector<int>& masterdofs,
        const std::vector<int>& slavedofs, const int nds_master = 0, const int nds_slave = 0);

    /// build dof maps from node maps
    /*!
      \note It is assumed that the first numdof dofs of each
      node as specified in vector coupled_dofs are of interest.
     */
    void build_dof_maps(const Core::FE::Discretization& dis, const Core::LinAlg::Map& nodemap,
        const Core::LinAlg::Map& permnodemap, std::shared_ptr<const Core::LinAlg::Map>& dofmap,
        std::shared_ptr<const Core::LinAlg::Map>& permdofmap,
        std::shared_ptr<Core::LinAlg::Export>& exporter, const std::vector<int>& coupled_dofs,
        const int nds = 0) const;

   protected:
    /// @name accessors to the private class members for derived classes
    /// @{

    /// access the interface DoF map of the target side
    std::shared_ptr<const Core::LinAlg::Map>& ma_dof_map_ptr();
    const Core::LinAlg::Map& ma_dof_map() const;

    /// access the permuted interface DoF map of the target side
    std::shared_ptr<const Core::LinAlg::Map>& permuted_ma_dof_map_ptr();
    const Core::LinAlg::Map& permuted_ma_dof_map() const;

    /// access the interface DoF map of the source side
    std::shared_ptr<const Core::LinAlg::Map>& sl_dof_map_ptr();
    const Core::LinAlg::Map& sl_dof_map() const;

    /// access the permuted interface DoF map of the source side
    std::shared_ptr<const Core::LinAlg::Map>& permuted_sl_dof_map_ptr();
    const Core::LinAlg::Map& permuted_sl_dof_map() const;

    /// access the permuted target dof map to target dof map exporter
    std::shared_ptr<Core::LinAlg::Export>& ma_exporter_ptr();
    const Core::LinAlg::Export& ma_exporter() const;

    /// access the permuted source dof map to source dof map exporter
    std::shared_ptr<Core::LinAlg::Export>& sl_exporter_ptr();
    const Core::LinAlg::Export& sl_exporter() const;

    /// @}


   private:
    /*! @name Fundamental dof maps
     *
     *  We keep the target and source dof map as well as permuted versions
     *  that match the respective other side.
     */
    //@{

    //! the interface dof map of the target side
    std::shared_ptr<const Core::LinAlg::Map> masterdofmap_;

    //! the permuted interface dof map of the target side
    std::shared_ptr<const Core::LinAlg::Map> permmasterdofmap_;

    //! the interface dof map of the source side
    std::shared_ptr<const Core::LinAlg::Map> slavedofmap_;

    //! the permuted interface dof map of the source side
    std::shared_ptr<const Core::LinAlg::Map> permslavedofmap_;

    //@}

    //! @name Communication object
    //@{

    //! permuted target dof map to target dof map exporter
    std::shared_ptr<Core::LinAlg::Export> masterexport_;

    //! permuted source dof map to source dof map exporter
    std::shared_ptr<Core::LinAlg::Export> slaveexport_;

    //@}

    //! @name coupling matrices for Lagrangian multiplier coupling
    //@{

    std::shared_ptr<Core::LinAlg::SparseMatrix> matmm_;
    std::shared_ptr<Core::LinAlg::SparseMatrix> matsm_;
    std::shared_ptr<Core::LinAlg::SparseMatrix> matmm_trans_;
    std::shared_ptr<Core::LinAlg::SparseMatrix> matsm_trans_;

    //@}
  };
}  // namespace Coupling::Adapter

FOUR_C_NAMESPACE_CLOSE

#endif
