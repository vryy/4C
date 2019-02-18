/*!
\file beam_to_solid_mortar_manager.H

\brief Manage the creation of additional DOFs for mortar couplings between beams and solids.

\level 3
\maintainer Ivo Steinbrecher
*/
// End doxygen header.


#ifndef BEAM_TO_SOLID_MORTAR_MANAGER_H_
#define BEAM_TO_SOLID_MORTAR_MANAGER_H_


#include <Teuchos_RCP.hpp>

#include "../drt_inpar/inpar_beaminteraction.H"


// Forward declarations.
class Epetra_Map;
class Epetra_MultiVector;
namespace DRT
{
  class Discretization;
}  // namespace DRT
namespace STR
{
  namespace TIMINT
  {
    class BaseDataGlobalState;
  }
}  // namespace STR
namespace BEAMINTERACTION
{
  class BeamContactParams;
  class BeamContactPair;
}  // namespace BEAMINTERACTION


namespace BEAMINTERACTION
{
  /**
   * \brief In beam to solid interactions with mortar contact discretization, we need to create a
   * map with the Lagrange multiplier DOFs (in contrast to solid meshtying / mortar we do not create
   * a own discretization for the interface).
   *
   * The created DOF can be split into two groups:
   *   - Lagrange multiplier DOFs on  nodes that are physical nodes of the system. The do not need
   *     to have the same number of nodal values as the physical node or even the same dimension
   *     (although in most cases the Lagrange multiplier have 3 components for each nodal value).
   *   - Lagrange multiplier DOFs on elements. For example when we have a two noded beam element and
   *     we want a quadratic interpolation of the Lagrange multipliers, we 'give' the element
   *     additional DOFs that represent the values at the middle node.
   *
   * By defining the Lagrange multipliers like described above, each additional DOF can be
   * identified by either the global id of the physical node it is defined on or by the global id of
   * the element it is defined on.
   *
   * This class manages the connection between the created nodes and the global node / element DOFs.
   * For the created maps a offset can be chosen, so the new DOFs fit into a global saddle-point
   * system.
   */
  class BeamToSolidMortarManager
  {
   public:
    /**
     * \brief Standard Constructor
     */
    BeamToSolidMortarManager(
        const INPAR::BEAMINTERACTION::BeamToSolidVolumeMortarShapefunctions& mortar_shape_function);

    /**
     * \brief This method builds the global maps for the global node / element IDs to the Lagrange
     * multiplier DOFs.
     *
     * Some nodes / elements in the discretization need additional Lagrange multiplier DOFs. We need
     * to be able to know which pair refers to which Lagrange multipliers. In this setup routine, a
     * Epetra multi vector is created, that maps all centerline nodes and beam elements, to a
     * Lagrange multiplier DOF.
     *
     * @param discret Pointer to the discretization.
     */
    void Setup(const Teuchos::RCP<DRT::Discretization> discret);

    /**
     * \brief This method builds the local maps from the global multi vector created in Setup.
     *
     * Since some nodes of this pair, that have Lagrange multipliers may not be on this processor,
     * we need to get the node ID to Lagrange multiplier ID form the processor that holds the
     * node. All relevant global node / element to global Lagrange multiplier maps for the given
     * contact pairs are stored in a standard maps in this object. The keys in those maps are the
     * global node / element id and the value is a vector with the corresponding Lagrange multiplier
     * gids. By doing so we only have to communicate between the ranks once per timestep (to be more
     * precise: only once for each set of contact pairs. If they do not change between timesteps and
     * do not switch rank, we can keep the created maps).
     *
     * @param discret Pointer to the discretization.
     * @param contact_pairs All contact pairs on this processor.
     */
    void SetLocalMaps(const Teuchos::RCP<DRT::Discretization> discret,
        const std::vector<Teuchos::RCP<BEAMINTERACTION::BeamContactPair>> contact_pairs);

    /**
     * \brief Get the global IDs of all Lagrange multipliers for the contact pair.
     * @param contact_pair (in) pointer to contact pair.
     * @param lambda_row (out) Standard vector with the global IDs of the Lagrange multipliers for
     * this pair.
     */
    void LocationVector(const Teuchos::RCP<BEAMINTERACTION::BeamContactPair> contact_pair,
        std::vector<int>& lambda_row) const;

   protected:
    /**
     * \brief Throw an error if setup was not called on the object prior to this function call.
     */
    inline void CheckSetup() const
    {
      if (!is_setup_) dserror("Setup not called on BeamToSolidMortarManager!");
    }

   private:
    //! Flag if setup was called.
    bool is_setup_;

    //! Base for the index used for the Lagrange multiplier DOFs.
    int index_base_;

    //! Number of Lagrange multiplier DOFs on a node.
    unsigned int n_lambda_node_;

    //! Number of Lagrange multiplier DOFs on an element.
    unsigned int n_lambda_element_;

    //! Rowmap of the additional Lagrange multiplier DOFs.
    Teuchos::RCP<Epetra_Map> lambda_dof_rowmap_;

    //! Multivector that connects the global node IDs with the Lagrange multiplier DOF IDs.
    Teuchos::RCP<Epetra_MultiVector> node_gid_to_lambda_gid_;

    //! Multivector that connects the global element IDs with the Lagrange multiplier DOF IDs.
    Teuchos::RCP<Epetra_MultiVector> element_gid_to_lambda_gid_;

    //! Standard map from global node ids to global Lagrange multiplier ids, for all
    //! nodes used on this rank.
    std::map<int, std::vector<int>> node_gid_to_lambda_gid_map_;

    //! Standard map from global element ids to global Lagrange multiplier ids, for all elements
    //! used on this rank.
    std::map<int, std::vector<int>> element_gid_to_lambda_gid_map_;
  };
}  // namespace BEAMINTERACTION

#endif
