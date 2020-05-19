/*----------------------------------------------------------------------*/
/*! \file

\brief Manage the creation of additional DOFs for mortar couplings between beams and fluid elements.

\level 3

\maintainer Nora Hagmeyer
*/
/*----------------------------------------------------------------------*/


#ifndef BEAM_TO_FLUID_MORTAR_MANAGER_H
#define BEAM_TO_FLUID_MORTAR_MANAGER_H


#include <Teuchos_RCP.hpp>

#include "../drt_inpar/inpar_beaminteraction.H"


// Forward declarations.
class Epetra_Map;
class Epetra_MultiVector;
class Epetra_Vector;
class Epetra_FEVector;

namespace FBI
{
  class BeamToFluidMeshtyingParams;
}
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

namespace LINALG
{
  class SparseMatrix;
}  // namespace LINALG

namespace BEAMINTERACTION
{
  /**
   * \brief Manage Lagrange mulitplier DOFs for BeamToFluid mortar coupling
   *
   * In beam to fluid interactions with mortar contact discretization, we need to create a
   * map with the Lagrange multiplier DOFs (in contrast to solid meshtying / mortar we do not create
   * a own discretization for the interface).
   *
   * The created DOF can be split into two groups:
   *   - Lagrange multiplier DOFs on  nodes that are physical nodes of the system. They do not need
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
   * The start value for the Lagrange multiplier global IDs can be explicitly given. This is usually
   * the number of solid DOFs + beam DOFs + Lagrange multipliers from other beam-to-fluid couplings
   * preceding this mortar manager in the model.
   * The Lagrange multiplier DOFs are then numbered the following way, and used in \ref
   * lambda_dof_rowmap_.
   *   - Lagrange multiplier DOFs on nodes of processor 0
   *   - Lagrange multiplier DOFs on elements of processor 0
   *   - Lagrange multiplier DOFs on nodes of processor 1
   *   - Lagrange multiplier DOFs on elements of processor 1
   *   - ...
   *
   * This class manages the connection between the created nodes and the global node / element DOFs.
   * For the created maps a offset can be chosen, so the new DOFs fit into a global saddle-point
   * system.
   */
  class BeamToFluidMortarManager
  {
   public:
    /**
     * \brief Standard Constructor
     *
     * @param[in] discretization1 Pointer to the structure discretization.
     * @param[in] discretization2 Pointer to the fluid discretization.
     * @param[in] params Parameters for the beam contact.
     * @param[in] start_value_lambda_gid Start value for the Lagrange multiplier global IDs.
     */
    BeamToFluidMortarManager(Teuchos::RCP<const DRT::Discretization> discretization1,
        Teuchos::RCP<const DRT::Discretization> discretization2,
        Teuchos::RCP<const FBI::BeamToFluidMeshtyingParams> params, int start_value_lambda_gid);

    /**
     * \brief This method builds the global maps for the global node / element IDs to the Lagrange
     * multiplier DOFs.
     *
     * Some nodes / elements in the discretization need additional Lagrange multiplier DOFs. We need
     * to be able to know which pair refers to which Lagrange multipliers. In this setup routine, a
     * Epetra multi vector is created, that maps all centerline nodes and beam elements, to a
     * Lagrange multiplier DOF.
     *
     */
    void Setup();

    /**
     * \brief Calculate the maps for the beam and fluid dofs. The calculated maps are used in
     * Complete of the mortar matrices.
     */
    void SetGlobalMaps();

    /**
     * \brief This method builds the local maps from the global multi vector created in Setup. The
     * global mortar matrices are also created.
     *
     * Since some nodes of this pair, that have Lagrange multipliers, may not be on this processor,
     * we need to get the node ID to Lagrange multiplier ID form the processor that holds the
     * node. All relevant global node / element to global Lagrange multiplier maps for the given
     * contact pairs are stored in a standard maps in this object. The keys in those maps are the
     * global node / element id and the value is a vector with the corresponding Lagrange multiplier
     * gids. By doing so we only have to communicate between the ranks once per timestep (to be more
     * precise: only once for each set of contact pairs. If they do not change between timesteps and
     * do not switch rank, we can keep the created maps).
     *
     * @param contact_pairs All contact pairs on this processor.
     */
    void SetLocalMaps(
        const std::vector<Teuchos::RCP<BEAMINTERACTION::BeamContactPair>>& contact_pairs);

    /**
     * \brief Get the global IDs of all Lagrange multipliers for the contact pair.
     * @param contact_pair (in) pointer to contact pair.
     * @param lambda_row (out) Standard vector with the global IDs of the Lagrange multipliers for
     * this pair.
     */
    void LocationVector(const Teuchos::RCP<const BEAMINTERACTION::BeamContactPair>& contact_pair,
        std::vector<int>& lambda_row) const;

    /**
     * \brief Evaluate D and M on all pairs and assemble them into the global matrices.
     * @param[in] contact_pairs Vector with all beam contact pairs in the model evaluator.
     */
    void EvaluateGlobalDM(
        const std::vector<Teuchos::RCP<BEAMINTERACTION::BeamContactPair>>& contact_pairs);

    /**
     * \brief Add the mortar penalty contributions to the global force vector and stiffness matrix.
     * @param[in] beam_vel Global beam velocity vector
     * @param[in] fluid_vel Global fluid velocity vector
     * @param[out] kbf Global stiffness matrix relating the fluid dofs to the structure residual
     * @param[out] kfb Global stiffness matrix relating the structure dofs to the fluid residual
     * @param[out] kbb Global stiffness matrix relating the structure dofs to the structure residual
     * @param[out] kff Global stiffness matrix relating the fluid dofs to the fluid residual
     * @param[out] fluid_force Global force vector acting on the fluid
     * @param[out] beam_force Global force vector acting on the beam
     */
    void AddGlobalForceStiffnessContributions(Teuchos::RCP<Epetra_FEVector> fluid_force,
        Teuchos::RCP<Epetra_FEVector> beam_force, Teuchos::RCP<LINALG::SparseMatrix> kbb,
        Teuchos::RCP<LINALG::SparseMatrix> kbf, Teuchos::RCP<LINALG::SparseMatrix> kff,
        Teuchos::RCP<LINALG::SparseMatrix> kfb, Teuchos::RCP<const Epetra_Vector> beam_vel,
        Teuchos::RCP<const Epetra_Vector> fluid_vel) const;

    /**
     * \brief Get the global vector of Lagrange multipliers.
     * @param[in] vel Global velocity vector.
     * @return Global vector of Lagrange multipliers.
     */
    Teuchos::RCP<Epetra_Vector> GetGlobalLambda(Teuchos::RCP<const Epetra_Vector> vel) const;

    /**
     * \brief Get the global vector of Lagrange multipliers, with the maps being the colum maps of
     * the Lagrange GID. on the ranks where they are used.
     * @param disp (in) Global displacement vector.
     * @return Global vector of Lagrange multipliers.
     */
    Teuchos::RCP<Epetra_Vector> GetGlobalLambdaCol(Teuchos::RCP<const Epetra_Vector> vel) const;

   protected:
    /**
     * \brief Throw an error if setup was not called on the object prior to this function call.
     */
    inline void CheckSetup() const
    {
      if (!is_setup_) dserror("Setup not called on BeamToSolidMortarManager!");
    }

    /**
     * \brief Throw an error if the global maps were not build.
     */
    inline void CheckGlobalMaps() const
    {
      if (!is_global_maps_build_) dserror("Global maps are not build in BeamToSolidMortarManager!");
    }

    /**
     * \brief Throw an error if the local maps were not build.
     */
    inline void CheckLocalMaps() const
    {
      if (!is_local_maps_build_) dserror("Local maps are not build in BeamToSolidMortarManager!");
    }

    /**
     * \brief Invert the scaling vector \ref global_kappa_ vector with accounting for non active
     * Lagrange multipliers.
     *
     * @return Inverted global_kappa_ vector.
     */
    Teuchos::RCP<Epetra_Vector> InvertKappa() const;

   private:
    //! Flag if setup was called.
    bool is_setup_;

    //! Flag if local maps were build.
    bool is_local_maps_build_;

    //! Flag if global maps were build.
    bool is_global_maps_build_;

    //! The start value for the Lagrange multiplier global IDs.
    int start_value_lambda_gid_;

    //! Number of Lagrange multiplier DOFs on a node.
    unsigned int n_lambda_node_;

    //! Number of Lagrange multiplier DOFs on an element.
    unsigned int n_lambda_element_;

    //! structure discretization
    Teuchos::RCP<const DRT::Discretization> discretization_structure_;

    //! fluid discretization
    Teuchos::RCP<const DRT::Discretization> discretization_fluid_;

    //! Pointer to the beam contact parameters.
    Teuchos::RCP<const FBI::BeamToFluidMeshtyingParams> beam_contact_parameters_ptr_;

    //! Row map of the additional Lagrange multiplier DOFs.
    Teuchos::RCP<Epetra_Map> lambda_dof_rowmap_;

    //! Column map of the additional Lagrange multiplier DOFs.
    Teuchos::RCP<Epetra_Map> lambda_dof_colmap_;

    //! Row map of the beam DOFs.
    Teuchos::RCP<Epetra_Map> beam_dof_rowmap_;

    //! Row map of the fluid DOFs.
    Teuchos::RCP<Epetra_Map> fluid_dof_rowmap_;

    //! Multivector that connects the global node IDs with the Lagrange multiplier DOF IDs.
    //! The global row ID of the multi vector is the global ID of the node that a Lagrange
    //! multiplier is defined on. The columns hold the corresponding global IDs of the Lagrange
    //! multipliers.
    Teuchos::RCP<Epetra_MultiVector> node_gid_to_lambda_gid_;

    //! Multivector that connects the global element IDs with the Lagrange multiplier DOF IDs.
    //! The global row ID of the multi vector is the global ID of the element that a Lagrange
    //! multiplier is defined on. The columns hold the corresponding global IDs of the Lagrange
    //! multipliers.
    Teuchos::RCP<Epetra_MultiVector> element_gid_to_lambda_gid_;

    //! Standard map from global node ids to global Lagrange multiplier ids, for all
    //! nodes used on this rank.
    std::map<int, std::vector<int>> node_gid_to_lambda_gid_map_;

    //! Standard map from global element ids to global Lagrange multiplier ids, for all elements
    //! used on this rank.
    std::map<int, std::vector<int>> element_gid_to_lambda_gid_map_;

    //! Global \f$D\f$ matrix.
    Teuchos::RCP<LINALG::SparseMatrix> global_D_;

    //! Global \f$M\f$ matrix.
    Teuchos::RCP<LINALG::SparseMatrix> global_M_;

    //! Global \f$\kappa\f$ vector. This vector is used to scale the mortar matrices. See Yang et
    //! al: Two dimensional mortar contact methods for large deformation frictional sliding (eq.
    //! 37).
    //! With this scaling correct units and pass patch tests are achieved (in the penalty case).
    Teuchos::RCP<Epetra_FEVector> global_kappa_;

    //! This vector keeps tack of all Lagrange multipliers that are active. This is needed when the
    //! kappa vector is inverted and some entries are zero, because no active contributions act on
    //! that Lagrange multiplier.
    Teuchos::RCP<Epetra_FEVector> global_active_lambda_;
  };
}  // namespace BEAMINTERACTION

#endif
