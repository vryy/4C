/*----------------------------------------------------------------------*/
/*! \file
 \brief one way coupled partitioned scalar structure interaction

 \level 2


 *------------------------------------------------------------------------------------------------*/

#ifndef SSI_PARTITIONED_1WC_H_
#define SSI_PARTITIONED_1WC_H_

#include "ssi_partitioned.H"

namespace SSI
{
  class SSIPart1WC : public SSIPart
  {
   public:
    //! constructor
    SSIPart1WC(const Epetra_Comm& comm, const Teuchos::ParameterList& globaltimeparams);

    //! virtual destructor
    ~SSIPart1WC() override = default;

    /*!
    \brief Setup this object

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

    //! actual time loop (implemented by derived class)
    void Timeloop() override = 0;

   protected:
    //! prepare time step for single fields
    void PrepareTimeStep(bool printheader = true) override = 0;

    //! perform one time step of structure field
    void DoStructStep() override;

    //! perform one time step of scatra field
    void DoScatraStep() override;

    // Flag for reading scatra result from restart files
    bool isscatrafromfile_;
  };

  class SSIPart1WCSolidToScatra : public SSIPart1WC
  {
   public:
    //! constructor
    SSIPart1WCSolidToScatra(
        const Epetra_Comm& comm, const Teuchos::ParameterList& globaltimeparams);

    //! virtual destructor
    ~SSIPart1WCSolidToScatra() override = default;

    /*!
    \brief Setup this object

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

    //! actual time loop
    void Timeloop() override;

    //! prepare time step for single fields
    void PrepareTimeStep(bool printheader = true) override;
  };

  class SSIPart1WCScatraToSolid : public SSIPart1WC
  {
   public:
    //! constructor
    SSIPart1WCScatraToSolid(
        const Epetra_Comm& comm, const Teuchos::ParameterList& globaltimeparams);

    //! virtual destructor
    ~SSIPart1WCScatraToSolid() override = default;

    /*!
    \brief Setup this object

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

    //! actual time loop
    void Timeloop() override;

    //! prepare time step for single fields
    void PrepareTimeStep(bool printheader = true) override;
  };

}  // namespace SSI

#endif /* SSI_PARTITIONED_1WC_H_ */
