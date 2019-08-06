/*----------------------------------------------------------------------*/
/*! \file
\file constraintenforcer_fbi.H

\brief Abstract class to be overloaded by different constraintenforcement techniques for fluid-beam
interaction.

\level 2

\maintainer Nora Hagmeyer
*----------------------------------------------------------------------*/
#ifndef CONSTRAINTENFORCER_FBI_H
#define CONSTRAINTENFORCER_FBI_H

#include <Teuchos_RCP.hpp>
#include "../linalg/linalg_utils.H"

namespace ADAPTER
{
  class FSIStructureWrapper;
  class FluidMovingBoundary;

}  // namespace ADAPTER
namespace GEO
{
  class SearchTree;
}
namespace ADAPTER
{
  class ConstraintEnforcerFactory;
  class FBIConstraintBridge;
  class FBIConstraintenforcer
  {
    friend ConstraintEnforcerFactory;

   public:
    /// empty destructor
    virtual ~FBIConstraintenforcer(){};

    /**
     * \brief Sets up the search tree.
     *
     */
    virtual void Setup(Teuchos::RCP<ADAPTER::FSIStructureWrapper> structure,
        Teuchos::RCP<ADAPTER::FluidMovingBoundary> fluid);

    /**
     * \brief Computes the coupling matrices
     *
     * This is where the magic happens. The stiffness contributions are integrated using information
     * of the beam elements, the fluid elements and their position relative to each other
     */

    virtual void Evaluate();

    /**
     * \brief Abstractly, we do everything we have to, to introduce the coupling condition into the
     * master field.
     *
     * Depending on the constraint enforcement strategy, either only an interface force is returned
     * (Mortar-Lagrangemultiplier partitioned, linearized penalty force partitioned), or force
     * vector with additional contributions as well as a stiffness matrix with additional
     * information is returned (monolithic formulation, full penalty partitioned, weak Dirichlet).
     *
     * What exactly has to be done, has to be implemented in the derived classes.
     *
     * \returns f residuum vector
     */

    virtual void ExtendBeamGhosting();

    virtual void PreparePairCreation(std::map<int, std::vector<int>>& pairids);

    virtual void CreatePairs(std::map<int, std::vector<int>> pairids);

    virtual Teuchos::RCP<Epetra_Vector> SlaveToMaster();

    /**
     * \brief Abstractly, we do everything we have to, to introduce the coupling condition into the
     * slave field.
     *
     * Depending on the constraint enforcement strategy, either only an interface force is returned
     * (Mortar-Lagrangemultiplier partitioned, linearized penalty force partitioned), or force
     * vector with additional contributions as well as a stiffness matrix with additional
     * information is returned (monolithic formulation, full penalty partitioned, weak Dirichlet).
     *
     * What exactly has to be done, has to be implemented in the derived classes.
     *
     * \returns f residuum vector
     */
    virtual Teuchos::RCP<Epetra_Vector> MasterToSlave();

    /// Interface to do preparations to solve the fluid
    virtual void PrepareFluidSolve(){};

    virtual Teuchos::RCP<const ADAPTER::FSIStructureWrapper> GetStructure() const final
    {
      return structure_;
    };
    virtual Teuchos::RCP<const ADAPTER::FluidMovingBoundary> GetFluid() const final
    {
      return fluid_;
    };
    virtual Teuchos::RCP<const GEO::SearchTree> GetSearchtree() const final { return searchtree_; };
    virtual const Teuchos::RCP<ADAPTER::FBIConstraintBridge> GetBridge() const final
    {
      return bridge_;
    };
    virtual Teuchos::RCP<const std::map<int, LINALG::Matrix<3, 1>>> GetFluidPositions() const final
    {
      return fluidpositions_;
    };
    virtual Teuchos::RCP<const std::map<int, LINALG::Matrix<3, 1>>> GetBeamPositions() const final
    {
      return beampositions_;
    };

    virtual Teuchos::RCP<const std::vector<Teuchos::RCP<DRT::Discretization>>> GetDiscretizations()
        const final
    {
      return discretizations_;
    }

   protected:
    /// Standard constructor
    FBIConstraintenforcer(Teuchos::RCP<ADAPTER::FBIConstraintBridge> bridge);

    // Do we need set-routines here?

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
     * \brief Computes the contributions to the stiffness matrix of the slave field.
     *
     * This has to be implemented differently depending on the concrete constraint enforcement
     * strategy.
     *
     * \param[inout] k system/stiffness matrix
     */
    virtual Teuchos::RCP<LINALG::SparseMatrix> AssembleMasterStiffness()
        const  // todo Is this really necessary?
    {
      dserror("Not yet implemented! This has to be overloaded by a derived class.\n");
      return Teuchos::null;
    };

    /**
     * \brief Computes the contributions to the stiffness matrix of the master field.
     *
     * This has to be implemented differently depending on the concrete constraint enforcement
     * strategy.
     *
     * \returns k system/stiffness matrix
     */
    virtual Teuchos::RCP<LINALG::SparseMatrix> AssembleSlaveStiffness()
        const  // todo Is this really necessary?
    {
      dserror("Not yet implemented! This has to be overloaded by a derived class.\n");
      return Teuchos::null;
    };

    /**
     * \brief Computes the contributions to the residuum of the master field.
     *
     * This has to be implemented differently depending on the concrete constraint enforcement
     * strategy.
     *
     * \returns f force/rhs matrix
     */
    virtual Teuchos::RCP<Epetra_Vector> AssembleMasterForce()
        const  // todo Is this really necessary?
    {
      dserror("Not yet implemented! This has to be overloaded by a derived class.\n");
      return Teuchos::null;
    };

    /**
     * \brief Computes the contributions to the residuum of the slave field.
     *
     * This has to be implemented differently depending on the concrete constraint enforcement
     * strategy.
     *
     * \returns f force/rhs matrix
     */
    virtual Teuchos::RCP<Epetra_Vector> AssembleSlaveForce()
        const  // todo Is this really necessary?
    {
      dserror("Not yet implemented! This has to be overloaded by a derived class.\n");
      return Teuchos::null;
    };

    /**
     * \brief Computes the reference nodal positions needed for the search
     *
     */
    virtual void ComputeFixedPositions(Teuchos::RCP<DRT::Discretization> dis,
        Teuchos::RCP<std::map<int, LINALG::Matrix<3, 1>>> positions) const;

    /**
     * \brief Computes the reference current positions needed for the search
     */
    virtual void ComputeCurrentPositions(Teuchos::RCP<DRT::Discretization> dis,
        Teuchos::RCP<std::map<int, LINALG::Matrix<3, 1>>> positions,
        Teuchos::RCP<const Epetra_Vector> disp) const;

   private:
    /// underlying fluid of the FSI problem
    Teuchos::RCP<ADAPTER::FluidMovingBoundary> fluid_;

    /// underlying structure of the FSI problem
    Teuchos::RCP<ADAPTER::FSIStructureWrapper> structure_;

    /// 3D seach tree for embedded discretization
    Teuchos::RCP<GEO::SearchTree> searchtree_;

    /// Map storing the nodal positions of the fluid for the search
    Teuchos::RCP<std::map<int, LINALG::Matrix<3, 1>>> fluidpositions_;

    /// Map storing the centerline positions of the beam for the search
    Teuchos::RCP<std::map<int, LINALG::Matrix<3, 1>>> beampositions_;

    /// Vector containing both field discretizations
    Teuchos::RCP<std::vector<Teuchos::RCP<DRT::Discretization>>> discretizations_;

    /**
     * \brief Object bridging the gap between the specific implementation of the constraint
     * enforcement technique and the specific implementation of the meshtying discretization
     * approach
     *
     * The ConstraintBridge connects the constraint enforcement technique with the
     * discretization approach. It does all the work. The correct matrices are computed and stored
     * in here. Depending on the constraint enforcement strategy, different constraint bridges will
     * be created in a separate factory This class serves as an interface to keep the constraint
     * enforcement technique (penalty, lagrange, ...) separate from the discretization appraoch
     * (GPTS, Mortar, ...)
     */
    Teuchos::RCP<ADAPTER::FBIConstraintBridge> bridge_;

    /// Flag determining the fluid stabilization type
    bool edgebased_fluidstabilization;

    Teuchos::RCP<const Epetra_Vector> column_structure_displacement_;
    Teuchos::RCP<const Epetra_Vector> column_structure_velocity_;
    Teuchos::RCP<const Epetra_Vector> column_fluid_velocity_;
  };
}  // namespace ADAPTER
#endif
