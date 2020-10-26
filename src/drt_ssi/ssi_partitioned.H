/*----------------------------------------------------------------------*/
/*! \file
 \brief base class for partitioned scalar structure interaction

 \level 2


 *------------------------------------------------------------------------------------------------*/

#ifndef SSI_PARTITIONED_H_
#define SSI_PARTITIONED_H_

#include "ssi_base.H"

namespace SSI
{
  class SSIPart : public SSIBase
  {
   public:
    /// setup SSI algorithm
    SSIPart(const Epetra_Comm& comm, const Teuchos::ParameterList& globaltimeparams);

    /// virtual destructor
    ~SSIPart() override = default;
    ;

    /*!
    \brief Initialize this object

     Initializes members and performs problem specific setup.

    \note Must only be called after parallel (re-)distribution of discretizations is finished !
          Otherwise, vectors may have wrong maps.

    \warning none
    \return void
    \date 08/16
    \author rauch
    */
    int Init(const Epetra_Comm& comm, const Teuchos::ParameterList& globaltimeparams,
        const Teuchos::ParameterList& scatraparams, const Teuchos::ParameterList& structparams,
        const std::string struct_disname, const std::string scatra_disname, bool isAle) override;

    /*!
    \brief Setup vectors and pointers

     Initializes members.

    \note Must only be called after parallel (re-)distribution of discretizations is finished !
          Otherwise, vectors may have wrong maps.

    \warning none
    \return void
    \date 08/16
    \author rauch
    */
    void Setup() override;

    /// setup of single fields (if needed)
    void SetupSystem() override;

    /// time loop of coupled problem
    void Timeloop() override = 0;

   protected:
    /// prepare time step of single fields
    virtual void PrepareTimeStep(bool printheader = true) = 0;

    /// do one time step (one way coupled) or one inner iteration loop step (two way coupled),
    /// depending on coupling algorithm
    virtual void DoStructStep() = 0;

    /// do one time step (one way coupled) or one inner iteration loop step (two way coupled),
    /// depending on coupling algorithm
    virtual void DoScatraStep() = 0;

   private:
    //! set up structural model evaluator for scalar-structure interaction
    void SetupModelEvaluator() const override;
  };

}  // namespace SSI

#endif /* SSI_PARTITIONED_H_ */
