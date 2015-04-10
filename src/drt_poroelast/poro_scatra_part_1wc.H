/*----------------------------------------------------------------------*/
/*!
 \file poro_scatra_part_1wc.H

 \brief  partitioned one way coupled poroelasticity scalar transport interaction algorithms

 <pre>
   Maintainer: Anh-Tu Vuong
               vuong@lnm.mw.tum.de
               http://www.lnm.mw.tum.de
               089 - 289-15251
 </pre>
 *----------------------------------------------------------------------*/


#ifndef PORO_SCATRA_PART_1WC_H_
#define PORO_SCATRA_PART_1WC_H_

#include "poro_scatra_part.H"

namespace POROELAST
{
  class PORO_SCATRA_Part_1WC : public PORO_SCATRA_Part
  {
  public:
    explicit PORO_SCATRA_Part_1WC(const Epetra_Comm& comm,
        const Teuchos::ParameterList& timeparams):
        PORO_SCATRA_Part(comm,timeparams){};

    //! solve one time step of porous media problem
    void DoPoroStep();
    //! solve one time step of scalar transport problem
    void DoScatraStep();

    //! increment time and step and print header
    void PrepareTimeStep();

    //! prepare output
    void PrepareOutput();

    //! update time step
    void Update();

    //! write output print to screen
    void Output();

  };

  class PORO_SCATRA_Part_1WC_PoroToScatra : public PORO_SCATRA_Part_1WC
  {
  public:
    //! constructor
    explicit PORO_SCATRA_Part_1WC_PoroToScatra(const Epetra_Comm& comm,
        const Teuchos::ParameterList& timeparams);

    //! virtual destructor
    virtual ~PORO_SCATRA_Part_1WC_PoroToScatra(){};

    //! actual time loop
    void Timeloop();

    //! perform iteration loop between fields
    void Solve();

    //! read and set fields needed for restart
    void ReadRestart(int restart);

  };

  class PORO_SCATRA_Part_1WC_ScatraToPoro : public PORO_SCATRA_Part_1WC
  {
  public:
    //! constructor
    explicit PORO_SCATRA_Part_1WC_ScatraToPoro(const Epetra_Comm& comm,
        const Teuchos::ParameterList& timeparams);

    //! virtual destructor
    virtual ~PORO_SCATRA_Part_1WC_ScatraToPoro(){};

    //! actual time loop
    void Timeloop();

    //! perform iteration loop between fields
    void Solve();

    //! read and set fields needed for restart
    void ReadRestart(int restart);

  };
}


#endif /* PORO_SCATRA_PART_1WC_H_ */
