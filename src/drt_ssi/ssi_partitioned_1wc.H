/*!------------------------------------------------------------------------------------------------*
 \file ssi_partitioned_1wc.H

 \brief one way coupled partitioned scalar structure interaction

 <pre>
   Maintainer: Anh-Tu Vuong
               vuong@lnm.mw.tum.de
               http://www.lnm.mw.tum.de
               089 - 289-15264
 </pre>
 *------------------------------------------------------------------------------------------------*/

#ifndef SSI_PARTITIONED_1WC_H_
#define SSI_PARTITIONED_1WC_H_

#include "ssi_partitioned.H"

namespace SSI
{
  class SSI_Part1WC : public SSI_Part
  {
  public:

    explicit SSI_Part1WC(const Epetra_Comm& comm,
        const Teuchos::ParameterList& timeparams);

    virtual ~SSI_Part1WC(){};

    virtual void Timeloop()=0;

  protected:
    void PrepareTimeStep();

    void DoStructStep();

    void DoScatraStep();
  };

  class SSI_Part1WC_SolidToScatra : public SSI_Part1WC
  {
  public:

    explicit SSI_Part1WC_SolidToScatra(const Epetra_Comm& comm,
        const Teuchos::ParameterList& timeparams);

    virtual ~SSI_Part1WC_SolidToScatra(){};

    virtual void Timeloop();

  };

  class SSI_Part1WC_ScatraToSolid : public SSI_Part1WC
  {
  public:

    explicit SSI_Part1WC_ScatraToSolid(const Epetra_Comm& comm,
        const Teuchos::ParameterList& timeparams);

    virtual ~SSI_Part1WC_ScatraToSolid(){};

    virtual void Timeloop();

  };

}

#endif /* SSI_PARTITIONED_1WC_H_ */
