/*----------------------------------------------------------------------*/
/*! \file

\brief Abstract class to be overloaded by different constraint enforcement techniques for fluid-beam
interaction.

\level 2

\maintainer Nora Hagmeyer
*----------------------------------------------------------------------*/
#ifndef CONSTRAINTENFORCER_FBI_H
#define CONSTRAINTENFORCER_FBI_H

#include "../drt_lib/drt_dserror.H"

#include <Teuchos_RCP.hpp>
#include <Epetra_Vector.h>
#include <vector>

namespace ADAPTER
{
  class FSIStructureWrapper;
  class FluidMovingBoundary;

}  // namespace ADAPTER

namespace DRT
{
  class Discretization;
  class Element;
}  // namespace DRT
namespace FBI
{
  class FBIGeometryCoupler;
}
namespace LINALG
{
  class SparseMatrix;
}
namespace ADAPTER
{
  class ConstraintEnforcerFactory;
  class FBIConstraintBridge;

  /**
   * \brief Abstract class to be overloaded by different constraint enforcement techniques for
   * fluid-beam interaction
   *
   * Depending on the constraint enforcement technique used to couple embedded meshes, e.g. through
   * the penalty method, Lagrange multiplier method, Nitsche method, ect., very different
   * information have to be passed to the participating fields.
   * This class is designed to decouple the decision of which information to pass from the actual
   * (partitioned) algorithm.
   *
   * The interface to the outside world is just a Setup and Evaluate routine, as well as the
   * FluidToStruct() and StructToFluid() routines, which only return one vector, as customary for
   * the FSI::DirichletNeumann algorithm. Everything else, like contributions to the stiffness
   * matrix of the fluid, it provides by calling internal functions. To make this possible, the
   * constraint enforcer need information on the two field.
   *
   */
  class FBIConstraintenforcer
  {
    friend ConstraintEnforcerFactory;

   public:
    /// empty destructor
    virtual ~FBIConstraintenforcer(){};

    /**
     * \brief Sets up the constraint enforcer
     *
     *\param[in] structure wrapper for the structure solver
     *\param[in] fluid moving boundary wrapper for the fluid solver
     */
    virtual void Setup(Teuchos::RCP<ADAPTER::FSIStructureWrapper> structure,
        Teuchos::RCP<ADAPTER::FluidMovingBoundary> fluid);

    /**
     * \brief Computes the coupling matrices
     *
     * This is where the magic happens. The stiffness contributions are integrated using information
     * of the structure elements, the fluid elements and their position relative to each other
     */

    virtual void Evaluate();

    /**
     * \brief Abstractly, we do everything we have to, to introduce the coupling condition into the
     * structure field.
     *
     * Depending on the constraint enforcement strategy, either only an interface force is returned
     * (Mortar-Lagrangemultiplier partitioned, linearized penalty force partitioned), or a force
     * vector as well as a stiffness matrix with additional
     * information is returned (monolithic formulation, full penalty partitioned).
     *
     *
     * \returns f residuum vector
     */

    virtual Teuchos::RCP<Epetra_Vector> FluidToStructure();

    /**
     * \brief Abstractly, we do everything we have to, to introduce the coupling condition into the
     * slave field.
     *
     * Depending on the constraint enforcement strategy, either only an interface force is returned
     * (Mortar-Lagrangemultiplier partitioned, linearized penalty force partitioned), or force
     * vector with additional contributions as well as a stiffness matrix with additional
     * information is returned (monolithic formulation, full penalty partitioned, weak Dirichlet).
     *
     *
     * \returns v velocity
     */

    virtual Teuchos::RCP<Epetra_Vector> StructureToFluid();

    /// Interface to do preparations to solve the fluid
    virtual void PrepareFluidSolve(){};

    /// Get function for the structure field strcuture_
    virtual Teuchos::RCP<const ADAPTER::FSIStructureWrapper> GetStructure() const final
    {
      return structure_;
    };
    /// Get function for the fluid field fluid_
    virtual Teuchos::RCP<const ADAPTER::FluidMovingBoundary> GetFluid() const final
    {
      return fluid_;
    };

    /// Get function for the bridge object bridge_
    virtual const Teuchos::RCP<ADAPTER::FBIConstraintBridge> GetBridge() const final
    {
      return bridge_;
    };

    /// Get function for the structure and the fluid discretization in the vector discretizations_
    virtual Teuchos::RCP<const std::vector<Teuchos::RCP<DRT::Discretization>>> GetDiscretizations()
        const final
    {
      return discretizations_;
    }

    virtual void PrintViolation();

   protected:
    /** \brief You will have to use the ADAPTER::ConstraintEnforcerFactory
     *
     * \param[in] bridge an object managing the pair contributins
     * \param[in] geometrycoupler an object managing the search, parallel communication, ect.
     */
    FBIConstraintenforcer(Teuchos::RCP<ADAPTER::FBIConstraintBridge> bridge,
        Teuchos::RCP<FBI::FBIGeometryCoupler> geometrycoupler);

    /**
     * \brief Creates all possible interaction pairs
     *
     * \param[in] pairids a map relating all beam element ids to a set of fluid
     * elements ids which they potentially cut
     */
    virtual void CreatePairs(Teuchos::RCP<std::map<int, std::vector<int>>> pairids);

    /**
     * \brief Extracts current element dofs that are needed for the computations on pair level
     *
     *\param[in] elements elements belonging to the pair
     *\param[out] beam_dofvec current positions and velocities of the beam element
     *\param[out] fluid_dofvec current positions and velocities of the fluid element
     */
    virtual void ExtractCurrentElementDofs(std::vector<DRT::Element const*> elements,
        Teuchos::RCP<std::vector<double>>& beam_dofvec,
        Teuchos::RCP<std::vector<double>>& fluid_dofvec) const;

    /**
     * \brief Computes the contributions to the stiffness matrix of the fluid field.
     *
     * This has to be implemented differently depending on the concrete constraint enforcement
     * strategy.
     *
     * \param[inout] k system/stiffness matrix
     */
    virtual Teuchos::RCP<LINALG::SparseMatrix> AssembleFluidStiffness() const
    {
      dserror("Not yet implemented! This has to be overloaded by a derived class.\n");
      return Teuchos::null;
    };

    /**
     * \brief Computes the contributions to the stiffness matrix of the structure field.
     *
     * This has to be implemented differently depending on the concrete constraint enforcement
     * strategy.
     *
     * \returns k system/stiffness matrix
     */
    virtual Teuchos::RCP<LINALG::SparseMatrix> AssembleStructureStiffness() const
    {
      dserror("Not yet implemented! This has to be overloaded by a derived class.\n");
      return Teuchos::null;
    };

    /**
     * \brief Computes the contributions to the rhs of the structure field.
     *
     * This has to be implemented differently depending on the concrete constraint enforcement
     * strategy.
     *
     * \returns f force vector acting on structure
     */
    virtual Teuchos::RCP<Epetra_Vector> AssembleStructureForce() const
    {
      dserror("Not yet implemented! This has to be overloaded by a derived class.\n");
      return Teuchos::null;
    };

    /**
     * \brief Computes the contributions to the residuum of the fluid field.
     *
     * This has to be implemented differently depending on the concrete constraint enforcement
     * strategy.
     *
     * \returns f force vector acting on fluid
     */
    virtual Teuchos::RCP<Epetra_Vector> AssembleFluidForce() const
    {
      dserror("Not yet implemented! This has to be overloaded by a derived class.\n");
      return Teuchos::null;
    };

   private:
    FBIConstraintenforcer() = delete;

    /// underlying fluid of the FSI problem
    Teuchos::RCP<ADAPTER::FluidMovingBoundary> fluid_;

    /// underlying structure of the FSI problem
    Teuchos::RCP<ADAPTER::FSIStructureWrapper> structure_;

    /// Vector containing both (fluid and structure) field discretizations
    Teuchos::RCP<std::vector<Teuchos::RCP<DRT::Discretization>>> discretizations_;

    /**
     * \brief Object bridging the gap between the specific implementation of the constraint
     * enforcement technique and the specific implementation of the meshtying discretization
     * approach
     */
    Teuchos::RCP<ADAPTER::FBIConstraintBridge> bridge_;

    /**
     * \brief Object handling geometric operations like the search of embedded pairs as well as the
     * parallel communication
     */
    Teuchos::RCP<FBI::FBIGeometryCoupler> geometrycoupler_;

    /// Displacement of the structural column nodes on the current proc
    Teuchos::RCP<const Epetra_Vector> column_structure_displacement_;
    /// Velocity of the structural column nodes on the current proc
    Teuchos::RCP<const Epetra_Vector> column_structure_velocity_;
    /// Velocity of the fluid column nodes on the current proc
    Teuchos::RCP<const Epetra_Vector> column_fluid_velocity_;
  };
}  // namespace ADAPTER
#endif
