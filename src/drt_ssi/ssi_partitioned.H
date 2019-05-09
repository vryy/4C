/*!------------------------------------------------------------------------------------------------*
 \brief base class for partitioned scalar structure interaction

 \level 2

 \maintainer Christoph Schmidt

 *------------------------------------------------------------------------------------------------*/

#ifndef SSI_PARTITIONED_H_
#define SSI_PARTITIONED_H_

#include "ssi_base.H"

namespace SSI
{
  class SSI_Part : public SSI_Base
  {
   public:
    /// setup SSI algorithm
    SSI_Part(const Epetra_Comm& comm, const Teuchos::ParameterList& globaltimeparams);

    /// virtual destructor
    virtual ~SSI_Part(){};

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
    virtual int Init(const Epetra_Comm& comm, const Teuchos::ParameterList& globaltimeparams,
        const Teuchos::ParameterList& scatraparams, const Teuchos::ParameterList& structparams,
        const std::string struct_disname, const std::string scatra_disname, bool isAle);

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
    virtual void Setup();

    /// setup of single fields (if needed)
    virtual void SetupSystem();

    /// time loop of coupled problem
    virtual void Timeloop() = 0;

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
    void SetupModelEvaluator() const;
  };

}  // namespace SSI

#endif /* SSI_PARTITIONED_H_ */
