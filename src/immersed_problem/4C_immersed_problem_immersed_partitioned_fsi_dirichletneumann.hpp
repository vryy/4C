/*----------------------------------------------------------------------*/
/*! \file

\brief partitioned immersed fsi algorithm for neumann-(dirichlet-neumann) like coupling

\level 2


*----------------------------------------------------------------------*/
#ifndef FOUR_C_IMMERSED_PROBLEM_IMMERSED_PARTITIONED_FSI_DIRICHLETNEUMANN_HPP
#define FOUR_C_IMMERSED_PROBLEM_IMMERSED_PARTITIONED_FSI_DIRICHLETNEUMANN_HPP

#include "4C_config.hpp"

#include "4C_fluid_ele_action.hpp"
#include "4C_immersed_problem_fsi_partitioned_immersed.hpp"
#include "4C_immersed_problem_immersed_base.hpp"

namespace Teuchos
{
  class Time;
}

FOUR_C_NAMESPACE_OPEN

namespace FLD
{
  class FluidImplicitTimeInt;
}

namespace CORE::GEO
{
  class SearchTree;
}

namespace ADAPTER
{
  class FSIStructureWrapperImmersed;
}

namespace IMMERSED
{
  class ImmersedPartitionedFSIDirichletNeumann : public ImmersedBase,
                                                 public FSI::PartitionedImmersed
  {
   public:
    /// setup partitioned ImmersedFSI algorithm
    explicit ImmersedPartitionedFSIDirichletNeumann(const Epetra_Comm& comm);

    /*! \brief Initialize this object

    \warning none
    \return int
    \date 08/16
    \author rauch  */
    int Init(const Teuchos::ParameterList& params) override;

    /*! \brief Setup all class internal objects and members

    \warning none
    \return void
    \date 08/16
    \author rauch  */
    void Setup() override;

    /// initialize search tree for structure discretization
    void setup_structural_discretization();

    /// read restart data
    void ReadRestart(int step) override;

   protected:
    /*!
    \brief composed FSI operator

    \note Relaxation for displacement coupled scheme is just experimental. By default displacement
    coupling + relaxation invokes a FOUR_C_THROW. NOX relaxes the vector x which is supposed to be
    the vector that is handed over to the other field. Since x in our case is the boundary
    displacement (in displacement coupled scheme) of the immersed structure but the vector applied
    to the fluid field is the volume Dirichlet stored in the vector fluid_artifical_velocity_,
    relaxation of x has no effect. That is why we get the factory from the fsi partitioned base
    class to extract the relaxation parameter from there. Now we are able to relax our artificial
    velocity in this class. This does not seem to work properly. For force coupled immersed fsi,
    everything works fine.

    */
    void FSIOp(const Epetra_Vector& x, Epetra_Vector& F, const FillType fillFlag) override;

    /// interface fluid operator
    Teuchos::RCP<Epetra_Vector> FluidOp(
        Teuchos::RCP<Epetra_Vector> fluid_artificial_velocity, const FillType fillFlag) override;

    /// interface structural operator
    Teuchos::RCP<Epetra_Vector> StructOp(
        Teuchos::RCP<Epetra_Vector> struct_bdry_traction, const FillType fillFlag) override;
    /// initial guess
    Teuchos::RCP<Epetra_Vector> InitialGuess() override;

    /// get immersed nodes and determine their dofs
    void build_immersed_dirich_map(Teuchos::RCP<DRT::Discretization> dis,
        Teuchos::RCP<Epetra_Map>& dirichmap,
        const Teuchos::RCP<const Epetra_Map>& dirichmap_original);

    /// add immersed dirichlet values from immersed dis to systemvector of background dis
    void do_immersed_dirichlet_cond(Teuchos::RCP<Epetra_Vector> statevector,
        Teuchos::RCP<Epetra_Vector> dirichvals, Teuchos::RCP<Epetra_Map> dbcmap);

    /// set state necessary state vectors
    virtual void SetStatesFluidOP();

    /// set state necessary state vectors
    virtual void set_states_velocity_correction();

    /// set state necessary state vectors
    virtual void SetStatesStructOP();

    /// call the nonlinear fluid solver
    virtual void SolveFluid();

    /// call the nonlinear structural solver
    virtual void SolveStruct();

    /// determine elements cut by the boundary
    void PrepareFluidOp();

    /// call to special extraction method
    virtual Teuchos::RCP<Epetra_Vector> extract_interface_dispnp();

    /// call to special application of interface forces
    virtual void apply_interface_forces(Teuchos::RCP<Epetra_Vector> full_traction_vec);

    /// call to special routine that adds dirichlet values to fluid field
    virtual void AddDirichCond();

    /// call to special routine that removes dirichlet values from fluid field
    virtual void RemoveDirichCond();

    /// calc the fsi residual
    virtual int CalcResidual(Epetra_Vector& F, const Teuchos::RCP<Epetra_Vector> newstate,
        const Teuchos::RCP<Epetra_Vector> oldstate);

    /*!
    \brief calc the current fluid tractions integrated over structural surface

    Calc new tractions on immersed boundary. After leaving this method the vector
    struct_bdry_traction_ contains the current tractions on the immersed boundary.

    */
    virtual void calc_fluid_tractions_on_structure();

    /*!
    \brief calc the current artificial velocity by projection from structure velocity onto fluid

    Calc new velocity in artificial fluid domain. After leaving this method the vector
    fluid_artificial_velocity_ contains the current velocity projected from the structure onto the
    fluid.

     \note Artificial velocity is only projected if the member artificial_velocity_isvalid_ is=false
    on entry. if artificial_velocity_isvalid_=true the current artificial velocity stored in
    fluid_artificial_velocity_ is returned. This is important when restart is requested since then
    the we need to perform the projection in the first restart step. This is done in \ref
    InitialGuess(). In case of no restart, we can just return the velocity since it had recently
    been projected after the last structural solve. This work would be done twice, else.
    */
    virtual Teuchos::RCP<Epetra_Vector> calc_artificial_velocity();

    /// apply given vector as Dirichlet to artificial fluid domain
    virtual void apply_immersed_dirichlet(Teuchos::RCP<Epetra_Vector> artificial_velocity);

    /// improve quality of solution near the interface
    virtual void correct_interface_velocity();

    /// reset immersed information in fluid dis
    /// e.g. isimmersed_, isboundaryimmersed_
    virtual void reset_immersed_information();


    //! @name Various global forces
    //@{
    Teuchos::RCP<Epetra_Vector> struct_bdry_traction_;  //!< bdry traction rhs on struct FSI surface
    Teuchos::RCP<Epetra_Vector> fluid_artificial_velocity_;  //!< background velocity interpolated
                                                             //!< from immersed dis (current state)
    Teuchos::RCP<Epetra_Map> dbcmap_immersed_;              //!< dirichlet bc map of immersed values
    Teuchos::RCP<CORE::GEO::SearchTree> fluid_SearchTree_;  //!< search tree for fluid domain
    Teuchos::RCP<CORE::GEO::SearchTree>
        structure_SearchTree_;  //!< search tree for structure domain
    std::map<int, std::set<int>>
        curr_subset_of_fluiddis_;  //!< fluid elements to evaluate the dirichlet interpolation
    std::map<int, CORE::LINALG::Matrix<3, 1>>
        currpositions_struct_;  //!< map of vectors for search tree containing current structural
                                //!< positions
    std::map<int, CORE::LINALG::Matrix<3, 1>>
        currpositions_fluid_;  //!< map of vectors for search tree containing current fluid
                               //!< positions
    //@}

    //! @name  basic information for parallelism
    int myrank_;
    int numproc_;
    //@}

    /// pointer to global problem
    GLOBAL::Problem* globalproblem_;

    //! @name various bools and switches
    //@{
    bool displacementcoupling_;  //!< true if increment of displacement is chosen for convergence
                                 //!< check of partitioned scheme.
    bool multibodysimulation_;   //!< true if more than one body is present. Each body needs to be
                                 //!< labeled with "ImmersedSearchbox" condition.
    bool output_evry_nlniter_;   //!< true if output after every iteration of the outer fsi loop.
                                 //!< Useful for debugging.
    bool is_relaxation_;         //!< true if relaxation is used
    int correct_boundary_velocities_;  //!< true if interface accuracy correction is requested
    int degree_gp_fluid_bound_;  //!< determines number of gp in fluid elements cut by boundary. gps
                                 //!< covered by structure are "compressible" gps.
    //@}

    //! @name various validation switches
    //@{
    bool artificial_velocity_isvalid_;  //!< this flag is true if the current vel. in artificial
                                        //!< fluid domain is consistent with current state of
                                        //!< structure
    bool boundary_traction_isvalid_;    //!< this flag is true if the current bdry. traction on the
                                        //!< structure is consistent with current state of fluid
    bool immersed_info_isvalid_;  //!< true if struct. configuration and fluid configuration are
                                  //!< consistently reflected in immersed infos, like e.g.
                                  //!< isimmersed_
    //@}

    //! @name pointer to discretizations
    //@{
    Teuchos::RCP<DRT::Discretization> fluiddis_;
    Teuchos::RCP<DRT::Discretization> structdis_;
    //@}

    /// pointer to immersed structure adapter
    Teuchos::RCP<ADAPTER::FSIStructureWrapperImmersed> immersedstructure_;

  };  // class ImmersedFSI
}  // namespace IMMERSED

FOUR_C_NAMESPACE_CLOSE

#endif
