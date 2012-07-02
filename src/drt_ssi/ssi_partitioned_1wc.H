/*
 * ssi_partitioned_1wc.H
 *
 *  Created on: Jun 29, 2012
 *      Author: vuong
 */

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

    virtual void Timeloop();

  protected:
    void PrepareTimeStep();

    void DoStructStep();

    void DoScatraStep();
  };

}

#endif /* SSI_PARTITIONED_1WC_H_ */
