/*----------------------------------------------------------------------*/
/*! \file
 \brief base class for partitioned scalar structure interaction

 \level 2


 *------------------------------------------------------------------------------------------------*/

#ifndef FOUR_C_SSI_PARTITIONED_HPP
#define FOUR_C_SSI_PARTITIONED_HPP

#include "baci_config.hpp"

#include "baci_ssi_base.hpp"

FOUR_C_NAMESPACE_OPEN

namespace SSI
{
  class SSIPart : public SSIBase
  {
   public:
    /// setup SSI algorithm
    SSIPart(const Epetra_Comm& comm, const Teuchos::ParameterList& globaltimeparams);

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
    void Init(const Epetra_Comm& comm, const Teuchos::ParameterList& globaltimeparams,
        const Teuchos::ParameterList& scatraparams, const Teuchos::ParameterList& structparams,
        const std::string& struct_disname, const std::string& scatra_disname, bool isAle) override;

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
    void SetupModelEvaluator() override;
  };

}  // namespace SSI

FOUR_C_NAMESPACE_CLOSE

#endif
