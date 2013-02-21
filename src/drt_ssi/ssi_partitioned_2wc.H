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

    explicit SSI_Part2WC(const Epetra_Comm& comm,
        const Teuchos::ParameterList& timeparams);

    virtual ~SSI_Part2WC(){};

    //full time loop
    virtual void Timeloop();

  protected:
    void PrepareTimeStep();

    void DoStructStep();

    void DoScatraStep();

    void OuterLoop();

    void UpdateAndOutput();

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
