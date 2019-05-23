/*----------------------------------------------------------------------*/
/*!

 \brief  partitioned one way coupled poroelasticity scalar transport interaction algorithms

\level 2

\maintainer Christoph Ager
            ager@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289 15249
 *----------------------------------------------------------------------*/


#ifndef PORO_SCATRA_PART_1WC_H_
#define PORO_SCATRA_PART_1WC_H_

#include "poro_scatra_part.H"

namespace POROELAST
{
  class PoroScatraPart1WC : public PoroScatraPart
  {
   public:
    explicit PoroScatraPart1WC(const Epetra_Comm& comm, const Teuchos::ParameterList& timeparams)
        : PoroScatraPart(comm, timeparams){};

    //! solve one time step of porous media problem
    void DoPoroStep();
    //! solve one time step of scalar transport problem
    void DoScatraStep();

    //! prepare output
    void PrepareOutput();

    //! update time step
    void Update();

    //! write output print to screen
    void Output();
  };

  class PoroScatraPart1WCPoroToScatra : public PoroScatraPart1WC
  {
   public:
    //! constructor
    explicit PoroScatraPart1WCPoroToScatra(
        const Epetra_Comm& comm, const Teuchos::ParameterList& timeparams);

    //! virtual destructor
    virtual ~PoroScatraPart1WCPoroToScatra(){};

    //! actual time loop
    void Timeloop();

    //! increment time and step and print header
    void PrepareTimeStep(bool printheader = true);

    //! perform iteration loop between fields
    void Solve();

    //! read and set fields needed for restart
    void ReadRestart(int restart);
  };

  class PoroScatraPart1WCScatraToPoro : public PoroScatraPart1WC
  {
   public:
    //! constructor
    explicit PoroScatraPart1WCScatraToPoro(
        const Epetra_Comm& comm, const Teuchos::ParameterList& timeparams);

    //! virtual destructor
    virtual ~PoroScatraPart1WCScatraToPoro(){};

    //! actual time loop
    void Timeloop();

    //! increment time and step and print header
    void PrepareTimeStep(bool printheader = true);

    //! perform iteration loop between fields
    void Solve();

    //! read and set fields needed for restart
    void ReadRestart(int restart);
  };
}  // namespace POROELAST


#endif /* PORO_SCATRA_PART_1WC_H_ */
