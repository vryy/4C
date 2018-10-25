/*!-----------------------------------------------------------------------------------------------*
\file scatra_timint_cardiac_monodomain.H

\brief scatra time integration for cardiac monodomain

\level 2

<pre>
\maintainer Lasse Jagschies
            lasse.jagschies@tum.de
            http://www.lnm.mw.tum.de
            089 - 289-10365
</pre>
 *------------------------------------------------------------------------------------------------*/

#ifndef SCATRA_TIMINT_CARDIAC_MONODOMAIN_H
#define SCATRA_TIMINT_CARDIAC_MONODOMAIN_H

#include "../drt_scatra/scatra_timint_implicit.H"


/*==========================================================================*/
// forward declarations
/*==========================================================================*/


namespace SCATRA
{
  class TimIntCardiacMonodomain : public virtual ScaTraTimIntImpl
  {
   public:
    /*========================================================================*/
    //! @name Constructors and destructors and related methods
    /*========================================================================*/

    //! Standard Constructor
    TimIntCardiacMonodomain(Teuchos::RCP<DRT::Discretization> dis,
        Teuchos::RCP<LINALG::Solver> solver, Teuchos::RCP<Teuchos::ParameterList> params,
        Teuchos::RCP<Teuchos::ParameterList> sctratimintparams,
        Teuchos::RCP<Teuchos::ParameterList> extraparams,
        Teuchos::RCP<IO::DiscretizationWriter> output);

    //! virtual destructor
    virtual ~TimIntCardiacMonodomain() { return; };

    //! setup algorithm
    virtual void Setup();

    //! time update of time-dependent materials
    void ElementMaterialTimeUpdate();

    void OutputState();

    //! Set ep-specific parameters
    void SetElementSpecificScaTraParameters(Teuchos::ParameterList& eleparams) const;

    //! return pointer to old scatra dofrowmap for coupled problems
    virtual Teuchos::RCP<const Epetra_Map> DofRowMapScatra() { return dofmap_; };

    /*========================================================================*/
    //! @name electrophysiology
    /*========================================================================*/

    //! activation_time at times n+1
    Teuchos::RCP<Epetra_Vector> activation_time_np_;

    //! activation threshold for postprocessing
    double activation_threshold_;

    //! maximum expected number of material internal state variables
    int nb_max_mat_int_state_vars_;

    //! material internal state at times n+1
    Teuchos::RCP<Epetra_MultiVector> material_internal_state_np_;

    //! one component of the material internal state at times n+1 (for separated postprocessing)
    Teuchos::RCP<Epetra_Vector> material_internal_state_np_component_;

    //! maximum expected number of material ionic currents
    int nb_max_mat_ionic_currents_;

    //! material ionic currents at times n+1
    Teuchos::RCP<Epetra_MultiVector> material_ionic_currents_np_;

    //! one component of the material ionic currents at times n+1 (for separated postprocessing)
    Teuchos::RCP<Epetra_Vector> material_ionic_currents_np_component_;

    //! parameter list
    const Teuchos::RCP<Teuchos::ParameterList> ep_params_;

    //! row map of dofs (one dof for each node)
    Teuchos::RCP<Epetra_Map> dofmap_;
  };

};  // namespace SCATRA

#endif
