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
    BeamToSolidMortarManager(Teuchos::RCP<BEAMINTERACTION::BeamContactParams> beam_contact_params);

    /**
     * \brief This method builds the maps for the global node / element ids to the Lagrange
     * multiplier DOFs.
     *
     * Some nodes / elements in the discretization need additional Lagrange multiplier DOFs. We need
     * to be able to know which pair refers to which Lagrange multipliers. This is done with a
     * standard map for nodes / elements, where the key is the global node / element id and the
     * value is a vector with the corresponding Lagrange multiplier gids.
     *
     * @param discret Pointer to the discretization.
     * @param contact_pairs All contact pairs on this processor.
     */
    void GenerateMaps(const Teuchos::RCP<DRT::Discretization> discret,
        const std::vector<Teuchos::RCP<BEAMINTERACTION::BeamContactPair>> contact_pairs);

   private:
    //! Base for the index used for the Lagrange multiplier DOFs.
    int index_base_;

    //! Pointer to the beam contact parameters.
    Teuchos::RCP<BEAMINTERACTION::BeamContactParams> beam_contact_params_;

    //! Standard map from global node ids to global Lagrange multiplier ids, for all nodes used on
    //! this rank.
    std::map<int, std::vector<int>> node_gid_to_lambda_gid_map_;

    //! Standard map from global element ids to global Lagrange multiplier ids, for all elements
    //! used on this rank.
    std::map<int, std::vector<int>> element_gid_to_lambda_gid_map_;
  };
}  // namespace BEAMINTERACTION

#endif
