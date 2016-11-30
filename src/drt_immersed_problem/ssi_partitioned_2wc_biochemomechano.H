/*!----------------------------------------------------------------------
\file ssi_partitioned_2wc_biochemomechano.H

\brief specialization of ssi2wc for biochemo-mechano coupled active cell material

\level 3

\maintainer  Andreas Rauch
             rauch@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 -15240

*----------------------------------------------------------------------*/
#ifndef SSI_PARTITIONED2_WC_BIOCHEMOMECHANO_H_
#define SSI_PARTITIONED2_WC_BIOCHEMOMECHANO_H_

#include "../drt_ssi/ssi_partitioned_2wc.H"

namespace POROELAST
{
  //class Monolithic;
  class PoroScatraBase;
}

namespace ADAPTER
{
  class FSIStructureWrapper;
}

namespace DRT
{
  class Element;
  class ImmersedFieldExchangeManager;
}

namespace SSI
{
  //! specialized class to deal with partitioned 2WC SSI during cell contraction in cell migration.
  class SSI_Part2WC_BIOCHEMOMECHANO : public SSI_Part2WC
  {

  public:

    //! constructor
    SSI_Part2WC_BIOCHEMOMECHANO(const Epetra_Comm& comm,
        const Teuchos::ParameterList& globaltimeparams);

    //!virtual destructor
    virtual ~SSI_Part2WC_BIOCHEMOMECHANO(){};

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
    virtual int Init(
        const Epetra_Comm& comm,
        const Teuchos::ParameterList& globaltimeparams,
        const Teuchos::ParameterList& scatraparams,
        const Teuchos::ParameterList& structparams,
        const std::string struct_disname,
        const std::string scatra_disname,
        bool isAle);

    //! special init entry method for this class
    bool Init(
        const Epetra_Comm& comm,
        const Teuchos::ParameterList& params,
        const Teuchos::ParameterList& globaltimeparams,
        const Teuchos::ParameterList& scatraparams,
        const Teuchos::ParameterList& structparams,
        const std::string struct_disname,
        const std::string scatra_disname);

    //! setup this class
    virtual void Setup();


  protected:

    //! update Phinp
    void UpdateScalars();

    //! perform iteration step of cell structure field
    virtual void DoStructStep();

    //! perform iteration step of biochemo scatra field
    virtual void DoScatraStep();

    //! prepare time step for single fields
    virtual void PrepareTimeStep(bool printheader = true);

    //! update time step and print to screen
    virtual void UpdateAndOutput();

    //! check convergence (absolute OR relative increments)
    virtual bool ConvergenceCheck(int itnum);

  private:

    //! pointer to poroelast subproblem (poro-scatra interaction)
    Teuchos::RCP<POROELAST::PoroScatraBase> poroscatra_subproblem_;

    //! communicator
    int myrank_;

    //! underlying structure of the SSI problem
    Teuchos::RCP< ::ADAPTER::FSIStructureWrapper>    specialized_structure_;

    //! pointer to ImmersedFieldExchangeManager singleton object
    DRT::ImmersedFieldExchangeManager* exchange_manager_ ;

    //! vector containing rate of change in actin concentration
    Teuchos::RCP<Epetra_MultiVector> ratesactin_;

    //! vector containing rate of change in 'ROCK' concentration
    Teuchos::RCP<Epetra_MultiVector> rates_;

    //! vector containing concentrations
    Teuchos::RCP<Epetra_MultiVector> phi_;


  }; // class
} // namespace SSI

#endif /* SSI_PARTITIONED2_WC_BIOCHEMOMECHANO_H_ */
