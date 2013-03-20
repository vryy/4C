/*----------------------------------------------------------------------*/
/*!
 \file poro_scatra_part_2wc.H

 \brief  partitioned two way coupled poroelasticity scalar transport interaction algorithms

 <pre>
   Maintainer: Anh-Tu Vuong
               vuong@lnm.mw.tum.de
               http://www.lnm.mw.tum.de
               089 - 289-15264
 </pre>
 *----------------------------------------------------------------------*/


#ifndef PORO_SCATRA_PART_2WC_H_
#define PORO_SCATRA_PART_2WC_H_

#include "poro_scatra_part.H"

namespace POROELAST
{
  class PORO_SCATRA_Part_2WC : public PORO_SCATRA_Part
  {
  public:

    //!explicit constructor
    explicit PORO_SCATRA_Part_2WC(const Epetra_Comm& comm,
        const Teuchos::ParameterList& timeparams);

    //!virtual destructor
    virtual ~PORO_SCATRA_Part_2WC(){};

    //!full time loop
    virtual void Timeloop();

  protected:
    //! prepare time step for single fields
    void PrepareTimeStep();

    //! perform iteration step of structure field
    void DoPoroStep();

    //! perform iteration step of scatra field
    void DoScatraStep();

    //! perform iteration loop between fields
    void OuterLoop();

    //! update time step and print to screen
    void UpdateAndOutput();

    //! convergence check of outer loop
    bool ConvergenceCheck(int itnum);

    //! scalar increment of the outer loop
    Teuchos::RCP<Epetra_Vector> scaincnp_;
    //! structure increment of the outer loop
    Teuchos::RCP<Epetra_Vector> structincnp_;
    //! fluid increment of the outer loop
    Teuchos::RCP<Epetra_Vector> fluidincnp_;

    //! maximum iteration steps
    int itmax_;
    //! convergence tolerance
    double ittol_;
  };

}


#endif /* PORO_SCATRA_PART_2WC_H_ */
