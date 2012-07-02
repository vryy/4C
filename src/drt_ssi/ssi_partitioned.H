/*
 * ssi_partitioned.H
 *
 *  Created on: Jun 28, 2012
 *      Author: vuong
 */

#ifndef SSI_PARTITIONED_H_
#define SSI_PARTITIONED_H_

#include "ssi_base.H"
#include "../drt_adapter/ad_str_wrapper.H"
#include "../drt_adapter/adapter_scatra_base_algorithm.H"

namespace SSI {



class SSI_Part : public SSI_Base
{
public:

  /// setup SSI algorithm
  explicit SSI_Part(const Epetra_Comm& comm,
      const Teuchos::ParameterList& timeparams);// Problem builder

  /// virtual destructor
  virtual ~SSI_Part(){};

  /// setup
  virtual void SetupSystem();

  /// timeloop of coupled problem
  virtual void Timeloop() = 0;

protected:

  virtual void PrepareTimeStep() = 0 ;

  //one time step or one inner iteration loop step, depending on coupling algorithm
  virtual void DoStructStep() = 0 ;

  //one time step or one inner iteration loop step, depending on coupling algorithm
  virtual void DoScatraStep() = 0 ;

  virtual void SetStructSolution() ;

  virtual void SetMeshDisp() ;

  virtual void SetVelocityFields() ;

private:

};

}

#endif /* SSI_PARTITIONED_H_ */
