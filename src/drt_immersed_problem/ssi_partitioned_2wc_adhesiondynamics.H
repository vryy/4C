/*----------------------------------------------------------------------*/
/*! \file

\brief specialization of ssi2wc, including adhesion dynamics

\level 2

\maintainer Jonas Eichinger
 *----------------------------------------------------------------------*/


#ifndef SSI_PARTITIONED2_WC_ADHESIONDYNAMICS_H_
#define SSI_PARTITIONED2_WC_ADHESIONDYNAMICS_H_


#include "../drt_ssi/ssi_partitioned_2wc.H"
#include "immersed_field_exchange_manager.H"


namespace SSI
{
  //! specialized class to deal with partitioned 2WC SSI during adhesion dynamics in cell migration.
  class SSI_Part2WC_ADHESIONDYNAMICS : public SSI_Part2WC
  {
   public:
    //! constructor (ctor)
    SSI_Part2WC_ADHESIONDYNAMICS(
        const Epetra_Comm& comm, const Teuchos::ParameterList& globaltimeparams);

    //! virtual destructor
    virtual ~SSI_Part2WC_ADHESIONDYNAMICS(){};

    /*! \brief Initialize this object

    Hand in all objects/parameters/etc. from outside.
    Construct and manipulate internal objects.

    \return int
    \date 08/16
    \author rauch  */
    virtual int Init(const Epetra_Comm& comm, const Teuchos::ParameterList& globaltimeparams,
        const Teuchos::ParameterList& scatraparams, const Teuchos::ParameterList& structparams,
        const std::string struct_disname, const std::string scatra_disname, bool isAle);

    /*!
    \brief Setup this object

     Initializes members and performs problem specific setup.

    \return void
    \date 08/16
    \author rauch  */
    virtual void Setup();

    //! update time step and print to screen
    virtual void UpdateAndOutput();


   protected:
    //! update bond traction vector
    void EvaluateBondTraction();

    /*! \brief build map containing row dofs of conditioned nodes

      \param condition  (in)  : condition to be considered
      \param dis        (in)  : discretization wearing the condition
      \param conddofmap (out) : dofrowmap for conditioned nodes
     */
    void BuildConditionDofRowMap(const DRT::Condition* condition,
        const Teuchos::RCP<const DRT::Discretization>& dis, Teuchos::RCP<Epetra_Map>& condnodemap);

    /*! \brief build map containing col dofs of conditioned nodes

      \param condition  (in)  : condition to be considered
      \param dis        (in)  : discretization wearing the condition
      \param conddofmap (out) : dofcolmap for conditioned nodes
     */
    void BuildConditionDofColMap(const DRT::Condition* condition,
        const Teuchos::RCP<const DRT::Discretization>& dis, Teuchos::RCP<Epetra_Map>& condnodemap);

    //! update the current states in every iteration
    //! states are set to the last solutions obtained
    virtual void IterUpdateStates();


   protected:
    //! solve field 1
    virtual void Operator1() { DoScatraStep(); };

    //! solve field 2
    virtual void Operator2() { DoStructStep(); };

    //! pre operator called before first field operator
    virtual void PreOperator1(){};

    //! pre operator called before second field operator
    virtual void PreOperator2();

    //! post operator called after first field operator
    virtual void PostOperator1();

    //! post operator called after second field operator
    virtual void PostOperator2();


   private:
    //! @name Various global vectors
    //@{
    Teuchos::RCP<Epetra_Vector> cell_adhesion_forces_;  //!< cell adhesion forces
    Teuchos::RCP<Epetra_Vector>
        surface_traction_;  //!< cell surface traction at adhesions row layout
    Teuchos::RCP<Epetra_Vector>
        surface_traction_col_;  //!< cell surface traction at adhesions col layout
    Teuchos::RCP<Epetra_Vector>
        cell_adhesion_fixpoints_;  //!< adhesion fixpoints (the concrete spatial coordinates to
                                   //!< which the integrins adhere)
    //@}

    //! @name Various global vectors
    //@{
    Teuchos::RCP<Epetra_Map>
        conditiondofrowmap_;  //!< map containing the dofs of conditioned adhesion nodes row layout
    Teuchos::RCP<Epetra_Map>
        conditiondofcolmap_;  //!< map containing the dofs of conditioned adhesion nodes col layout
    //@}

    //! manager for data exchange between fields
    DRT::ImmersedFieldExchangeManager* exchange_manager_;

  };  // class SSI_Part2WC_ADHESIONDYNAMICS


}  // namespace SSI

#endif /* SSI_PARTITIONED2_WC_ADHESIONDYNAMICS_H_ */
