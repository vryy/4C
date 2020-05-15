/*----------------------------------------------------------------------*/
/*! \file

\brief Data container holding all beam to fluid volume meshtying input parameters.

\level 1

\maintainer Nora Hagmeyer
*/


#ifndef BEAM_TO_FLUID_MESHTYING_PARAMS_H
#define BEAM_TO_FLUID_MESHTYING_PARAMS_H

#include "../drt_beaminteraction/beam_contact_params.H"
#include "../drt_fem_general/drt_utils_integration.H"

namespace INPAR
{
  namespace FBI
  {
    enum class BeamToFluidConstraintEnforcement;
    enum class BeamToFluidDiscretization;
    enum class BeamToFluidMeshtingMortarShapefunctions;
  }  // namespace FBI
}  // namespace INPAR
namespace FBI
{
  /**
   * \brief class containing all relevant parameters for beam to fluid meshtying
   */
  class BeamToFluidMeshtyingVtkOutputParams;

  /**
   * \brief Class for beam to fluid meshtying parameters.
   */
  class BeamToFluidMeshtyingParams : public BEAMINTERACTION::BeamContactParams
  {
   public:
    /**
     * \brief Constructor.
     */
    BeamToFluidMeshtyingParams();

    /**
     * \brief empty Destructor
     */
    virtual ~BeamToFluidMeshtyingParams(){};

    /**
     * \brief Initialize with the stuff coming from input file
     */
    void Init();

    /**
     * \brief Setup
     */
    void Setup();

    /// Sets the flag to compute only force contributions from the beam
    void SetWeakDirichletFlag() { calcfluidweakdirichletforce_ = true; }

    /// Sets the flag to compute force contributions from beam and fluid
    void UnsetWeakDirichletFlag() { calcfluidweakdirichletforce_ = false; }

    /// Returns \ref calcfluidweakdirichletforce_
    bool GetWeakDirichletFlag() { return calcfluidweakdirichletforce_; }

    /**
     * \brief Returns the isinit_ flag.
     */
    inline const bool& IsInit() const { return isinit_; };

    /**
     * \brief Returns the issetup_ flag.
     */
    inline const bool& IsSetup() const { return issetup_; };

    /**
     * \brief Checks the init and setup status.
     */
    inline void CheckInitSetup() const
    {
      if (!IsInit() or !IsSetup()) dserror("Call Init() and Setup() first!");
    }

    /**
     * \brief Checks the init status.
     */
    inline void CheckInit() const
    {
      if (!IsInit()) dserror("Init() has not been called, yet!");
    }

    /**
     * \brief Returns the contact discretization method.
     */
    inline INPAR::FBI::BeamToFluidConstraintEnforcement GetConstraintEnforcement() const
    {
      return constraint_enforcement_;
    }

    /**
     * \brief Returns constraints enforcement strategy.
     */
    inline INPAR::FBI::BeamToFluidDiscretization GetContactDiscretization() const
    {
      return meshtying_discretization_;
    }

    /**
     * \brief Returns the penalty parameter.
     * \returns penalty parameter.
     */
    inline double GetPenaltyParameter() const { return penalty_parameter_; }

    /**
     * \brief Returns the Gauss rule.
     * \returns gauss rule.
     */
    inline DRT::UTILS::GaussRule1D GetGaussRule() const { return gauss_rule_; }

    /**
     * \brief Returns a pointer to the vtk output parameters.
     * @return Pointer to vtk output parameters.
     */
    const Teuchos::RCP<const FBI::BeamToFluidMeshtyingVtkOutputParams> GetVtkOuputParamsPtr() const
    {
      return output_params_;
    }

    /**
     * \brief Returns the shape function for the mortar Lagrange-multiplicators.
     */
    inline INPAR::FBI::BeamToFluidMeshtingMortarShapefunctions GetMortarShapeFunctionType() const
    {
      return mortar_shape_function_;
    }


   private:
    /// Flag if Init was called.
    bool isinit_;

    /// Flag if Setup was called.
    bool issetup_;

    /// Enforcement strategy for constraints.
    INPAR::FBI::BeamToFluidConstraintEnforcement constraint_enforcement_;

    /// Discretization used for the contact.
    INPAR::FBI::BeamToFluidDiscretization meshtying_discretization_;

    /// Penalty parameter.
    double penalty_parameter_;

    /// Gauss rule to be used.
    DRT::UTILS::GaussRule1D gauss_rule_;

    /**
     * \brief Flag to keep track if the RHS contribution for the weak Dirichlet enforcement of the
     * kinematic continuity condition is requested
     *
     * In the case of a DirichletNeumann algorithm the right-hand side contribution C_fs*v_s leads
     * to the unnecessary costly creation of a non-local sparse matrix and matrix-vector product.
     * Instead, the contributions of the fluid "force" introduced by the beam DOFs can be calculated
     * directly on pair level
     */
    bool calcfluidweakdirichletforce_;

    /// VTK output params. For now I see no reason to overload this
    Teuchos::RCP<FBI::BeamToFluidMeshtyingVtkOutputParams> output_params_;

    //! Shape function for the mortar Lagrange-multiplicators
    INPAR::FBI::BeamToFluidMeshtingMortarShapefunctions mortar_shape_function_;
  };

}  // namespace FBI

#endif
