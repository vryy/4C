/*!------------------------------------------------------------------------------------------------*
 \file ssi_partitioned_2wc.H

 \brief two way coupled partitioned scalar structure interaction

 <pre>
   Maintainer: Anh-Tu Vuong
               vuong@lnm.mw.tum.de
               http://www.lnm.mw.tum.de
               089 - 289-15264
 </pre>
 *------------------------------------------------------------------------------------------------*/

#ifndef SSI_PARTITIONED2_WC_H_
#define SSI_PARTITIONED2_WC_H_

#include "ssi_partitioned.H"

namespace SSI
{
  class SSI_Part2WC : public SSI_Part
  {
  public:

    //!explicit constructor
    explicit SSI_Part2WC(const Epetra_Comm& comm,
        const Teuchos::ParameterList& globaltimeparams,
        const Teuchos::ParameterList& scatraparams,
        const Teuchos::ParameterList& structparams);

    //!virtual destructor
    virtual ~SSI_Part2WC(){};

    //!full time loop
    virtual void Timeloop();

  protected:
    //! prepare time step for single fields
    void PrepareTimeStep();

    //! perform iteration step of structure field
    void DoStructStep();

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
    //! displacement increment of the outer loop
    Teuchos::RCP<Epetra_Vector> dispincnp_;

    //! maximum iteration steps
    int itmax_;
    //! convergence tolerance
    double ittol_;
  };

}

#endif /* SSI_PARTITIONED2_WC_H_ */
