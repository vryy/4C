/*!------------------------------------------------------------------------------------------------*
 \file ssi_partitioned.H

 \brief base class for partitioned scalar structure interaction

 <pre>
   Maintainer: Anh-Tu Vuong
               vuong@lnm.mw.tum.de
               http://www.lnm.mw.tum.de
               089 - 289-15264
 </pre>
 *------------------------------------------------------------------------------------------------*/

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

  virtual void SetScatraSolution() ;

  virtual void SetMeshDisp() ;

  virtual void SetVelocityFields() ;

private:

};

}

#endif /* SSI_PARTITIONED_H_ */
