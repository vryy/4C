/*----------------------------------------------------------------------*/
/*! \file

 \brief  partitioned one way coupled poroelasticity scalar transport interaction algorithms

\level 2

 *----------------------------------------------------------------------*/


#ifndef BACI_POROELAST_SCATRA_PART_1WC_H
#define BACI_POROELAST_SCATRA_PART_1WC_H

#include "baci_config.H"

#include "baci_poroelast_scatra_part.H"

BACI_NAMESPACE_OPEN

namespace POROELASTSCATRA
{
  class PoroScatraPart1WC : public PoroScatraPart
  {
   public:
    explicit PoroScatraPart1WC(const Epetra_Comm& comm, const Teuchos::ParameterList& timeparams)
        : PoroScatraPart(comm, timeparams){};

    //! solve one time step of porous media problem
    void DoPoroStep() override;
    //! solve one time step of scalar transport problem
    void DoScatraStep() override;

    //! prepare output
    void PrepareOutput() override;

    //! update time step
    void Update() override;

    //! write output print to screen
    void Output() override;
  };

  class PoroScatraPart1WCPoroToScatra : public PoroScatraPart1WC
  {
   public:
    //! constructor
    explicit PoroScatraPart1WCPoroToScatra(
        const Epetra_Comm& comm, const Teuchos::ParameterList& timeparams);

    //! actual time loop
    void Timeloop() override;

    //! increment time and step and print header
    void PrepareTimeStep(bool printheader = true) override;

    //! perform iteration loop between fields
    void Solve() override;

    //! read and set fields needed for restart
    void ReadRestart(int restart) override;
  };

  class PoroScatraPart1WCScatraToPoro : public PoroScatraPart1WC
  {
   public:
    //! constructor
    explicit PoroScatraPart1WCScatraToPoro(
        const Epetra_Comm& comm, const Teuchos::ParameterList& timeparams);

    //! actual time loop
    void Timeloop() override;

    //! increment time and step and print header
    void PrepareTimeStep(bool printheader = true) override;

    //! perform iteration loop between fields
    void Solve() override;

    //! read and set fields needed for restart
    void ReadRestart(int restart) override;
  };
}  // namespace POROELASTSCATRA


BACI_NAMESPACE_CLOSE

#endif  // POROELAST_SCATRA_PART_1WC_H
