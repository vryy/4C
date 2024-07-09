/*---------------------------------------------------------------------*/
/*! \file

\brief A set of degrees of freedom special for contact

\level 1


*/
/*---------------------------------------------------------------------*/

#ifndef FOUR_C_FEM_DOFSET_TRANSPARENT_HPP
#define FOUR_C_FEM_DOFSET_TRANSPARENT_HPP


#include "4C_config.hpp"

#include "4C_comm_exporter.hpp"
#include "4C_fem_discretization.hpp"
#include "4C_fem_dofset.hpp"

FOUR_C_NAMESPACE_OPEN


namespace Core::DOFSets
{
  /// Alias dofset that shares dof numbers with another dofset
  /*!
  A special set of degrees of freedom, implemented in order to assign the same degrees of freedom to
  nodes belonging to two discretizations. This way two discretizations can assemble into the same
  position of the system matrix. As internal variable it holds a source discretization
  (Constructor). If such a nodeset is assigned to a sub-discretization, its dofs are assigned
  according to the dofs of the source.

  */
  class TransparentDofSet : public virtual Core::DOFSets::DofSet
  {
   public:
    /*!
    \brief Standard Constructor
    */
    explicit TransparentDofSet(
        Teuchos::RCP<Core::FE::Discretization> sourcedis, bool parallel = false);



    /// create a copy of this object
    Teuchos::RCP<DofSet> clone() override { return Teuchos::rcp(new TransparentDofSet(*this)); }

    /// Assign dof numbers to all elements and nodes of the discretization.
    int assign_degrees_of_freedom(
        const Core::FE::Discretization& dis, const unsigned dspos, const int start) override;

    /// Assign dof numbers for new discretization using dof numbering from source discretization.
    void transfer_degrees_of_freedom(const Core::FE::Discretization& sourcedis,  ///< source discret
        const Core::FE::Discretization&
            newdis,      ///< discretization that gets dof numbering from source discret
        const int start  ///< offset for dof numbering (obsolete)
    );


    /// Assign dof numbers for new discretization using dof numbering from source discretization.
    /// for this version, newdis is allowed to be distributed completely different; the
    /// communication  of the dofs is done internally.
    void parallel_transfer_degrees_of_freedom(
        const Core::FE::Discretization& sourcedis,  ///< source discret
        const Core::FE::Discretization&
            newdis,      ///< discretization that gets dof numbering from source discret
        const int start  ///< offset for dof numbering (obsolete)
    );

    /// helper for parallel_transfer_degrees_of_freedom; unpack the received block to
    /// generate the current map node gid -> its dofs
    void unpack_local_source_dofs(
        std::map<int, std::vector<int>>& gid_to_dofs, std::vector<char>& rblock);

    /// helper for parallel_transfer_degrees_of_freedom; pack the current map
    /// node gid -> its dofs into a send block
    void pack_local_source_dofs(
        std::map<int, std::vector<int>>& gid_to_dofs, Core::Communication::PackBuffer& sblock);

    /// helper for parallel_transfer_degrees_of_freedom; add processor local information
    /// to the map unpack the received block to the current map node gid -> its dofs
    void set_source_dofs_available_on_this_proc(std::map<int, std::vector<int>>& gid_to_dofs);

    /// helper for parallel_transfer_degrees_of_freedom, an MPI send call
    void send_block(int numproc, int myrank, std::vector<char>& sblock,
        Core::Communication::Exporter& exporter, MPI_Request& request);

    /// helper for parallel_transfer_degrees_of_freedom, an MPI receive call
    void receive_block(int numproc, int myrank, std::vector<char>& rblock,
        Core::Communication::Exporter& exporter, MPI_Request& request);

   protected:
    Teuchos::RCP<Core::FE::Discretization> sourcedis_;  ///< source discretization

    bool parallel_;  ///< call parallel_transfer_degrees_of_freedom instead of
                     ///< transfer_degrees_of_freedom

  };  // class TransparentDofSet
}  // namespace Core::DOFSets

FOUR_C_NAMESPACE_CLOSE

#endif
