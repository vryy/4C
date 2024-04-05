/*--------------------------------------------------------------------------*/
/*! \file

\brief base class for all elastohydrodynamic lubrication (lubrication structure interaction)
algorithms

\level 3

*/
/*--------------------------------------------------------------------------*/
#ifndef FOUR_C_EHL_BASE_HPP
#define FOUR_C_EHL_BASE_HPP


#include "baci_config.hpp"

#include "baci_adapter_algorithmbase.hpp"
#include "baci_coupling_adapter.hpp"
#include "baci_inpar_ehl.hpp"
#include "baci_linalg_utils_sparse_algebra_math.hpp"

#include <Epetra_Vector.h>

BACI_NAMESPACE_OPEN

// forward declarations
namespace ADAPTER
{
  class LubricationBaseAlgorithm;
  class Structure;
  class Coupling;
  class CouplingEhlMortar;
}  // namespace ADAPTER

namespace CORE::LINALG
{
  class MapExtractor;
}

namespace EHL
{
  class Base : public ADAPTER::AlgorithmBase
  {
   public:
    /// create using a Epetra_Comm
    explicit Base(const Epetra_Comm& comm, const Teuchos::ParameterList& globaltimeparams,
        const Teuchos::ParameterList& lubricationparams, const Teuchos::ParameterList& structparams,
        const std::string struct_disname,
        const std::string lubrication_disname);  // Problem builder

    /// setup
    virtual void SetupSystem() = 0;

    /// timeloop of coupled problem
    virtual void Timeloop() = 0;

    /// test results (if necessary)
    void TestResults(const Epetra_Comm& comm);

    /// read restart
    void ReadRestart(int restart) override;

    //! access to structural field
    const Teuchos::RCP<ADAPTER::Structure>& StructureField() { return structure_; }

    /// set structure solution on lubrication field
    void SetStructSolution(Teuchos::RCP<const Epetra_Vector> disp);

    /// set lubrication solution on structure field
    void SetLubricationSolution(Teuchos::RCP<const Epetra_Vector> pressure);

    /// evaluate fluid forces on structure
    Teuchos::RCP<Epetra_Vector> EvaluateFluidForce(Teuchos::RCP<const Epetra_Vector> pressure);

   protected:
    void AddPressureForce(
        Teuchos::RCP<Epetra_Vector> slaveiforce, Teuchos::RCP<Epetra_Vector> masteriforce);

    void AddPoiseuilleForce(
        Teuchos::RCP<Epetra_Vector> slaveiforce, Teuchos::RCP<Epetra_Vector> masteriforce);

    void AddCouetteForce(
        Teuchos::RCP<Epetra_Vector> slaveiforce, Teuchos::RCP<Epetra_Vector> masteriforce);

    /// underlying structure of the EHL problem
    Teuchos::RCP<ADAPTER::Structure> structure_;

    /// underlying lubrication problem of the EHL problem
    Teuchos::RCP<ADAPTER::LubricationBaseAlgorithm> lubrication_;

    //! Type of coupling strategy between the two fields of the EHL problems
    const INPAR::EHL::FieldCoupling fieldcoupling_;

    //! adapter for coupling the nodes of the lubrication field with the nodes from the master side
    //! of the structure
    Teuchos::RCP<ADAPTER::CouplingEhlMortar> mortaradapter_;

    //! Interface traction vector in the slave str dof map
    Teuchos::RCP<Epetra_Vector> stritraction_D_;
    Teuchos::RCP<Epetra_Vector> stritraction_M_;

    //! Transformation matrix for lubrication pre dof map <-> lubrication disp dof map
    Teuchos::RCP<CORE::LINALG::SparseMatrix> lubrimaptransform_;

    //! Mapextractors for dealing with interface vectors of the structure field
    Teuchos::RCP<CORE::LINALG::MapExtractor> slaverowmapextr_;
    Teuchos::RCP<CORE::LINALG::MapExtractor> masterrowmapextr_;
    Teuchos::RCP<CORE::LINALG::MapExtractor> mergedrowmapextr_;

    //! Transformation matrix for slave side node map <-> slave side disp dof map
    Teuchos::RCP<CORE::LINALG::SparseMatrix> slavemaptransform_;

    //! several adapters to transform maps
    Teuchos::RCP<CORE::ADAPTER::Coupling> ada_strDisp_to_lubDisp_;
    Teuchos::RCP<CORE::ADAPTER::Coupling> ada_strDisp_to_lubPres_;
    Teuchos::RCP<CORE::ADAPTER::Coupling> ada_lubPres_to_lubDisp_;

    //! height old vector to calculate the time derivative of height (Squeeze term)
    Teuchos::RCP<const Epetra_Vector> heightold_;

    //! use of a dry contact model
    bool dry_contact_;

    /// setup adapters for EHL on boundary
    virtual void SetupFieldCoupling(
        const std::string struct_disname, const std::string lubrication_disname);

    //! take current results for converged and save for next time step
    void Update() override;

    //! write output
    virtual void Output(bool forced_writerestart = false);
    Teuchos::RCP<Epetra_Vector> inf_gap_toggle_lub_;

    /// velocity calculation given the displacements
    Teuchos::RCP<Epetra_Vector> CalcVelocity(Teuchos::RCP<const Epetra_Vector> dispnp);

   private:
    /// setup discretizations and dofsets
    void SetupDiscretizations(const Epetra_Comm& comm, const std::string struct_disname,
        const std::string lubrication_disname);

    /// set structure mesh displacement on lubrication field
    void SetMeshDisp(Teuchos::RCP<const Epetra_Vector> disp);

    /// set average tangential interface velocity (ie structural velocities
    /// this is invariant w.r.t. rigid body rotations
    void SetAverageVelocityField();

    /// set relative tangential interface velocity (ie structural velocities
    /// this is invariant w.r.t. rigid body rotations
    void SetRelativeVelocityField();

    /// set film height of the lubrication field
    void SetHeightField();

    /// set Time derivative of height (squeeze term)
    void SetHeightDot();

    /// Create DBC map for unprojectable nodes
    void SetupUnprojectableDBC();
  };
}  // namespace EHL


BACI_NAMESPACE_CLOSE

#endif  // EHL_BASE_H
