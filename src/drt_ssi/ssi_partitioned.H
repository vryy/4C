/*!------------------------------------------------------------------------------------------------*
 \file ssi_partitioned.H

 \brief base class for partitioned scalar structure interaction

 \level 1

 \maintainer Anh-Tu Vuong
             vuong@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289-15264

 *------------------------------------------------------------------------------------------------*/

#ifndef SSI_PARTITIONED_H_
#define SSI_PARTITIONED_H_

#include "ssi_base.H"

namespace SSI {

class SSI_Part : public SSI_Base
{
public:

  /// setup SSI algorithm
  SSI_Part(const Epetra_Comm& comm,
      const Teuchos::ParameterList& globaltimeparams);

  /// virtual destructor
  virtual ~SSI_Part(){};

  /*!
  \brief Setup this object

   Initializes members and performs problem specific setup.

  \note Must only be called after parallel (re-)distribution of discretizations is finished !
        Otherwise, vectors may have wrong maps.

  \warning none
  \return void
  \date 08/16
  \author rauch
  */
  virtual int Init(const Epetra_Comm& comm,
      const Teuchos::ParameterList& globaltimeparams,
      const Teuchos::ParameterList& scatraparams,
      const Teuchos::ParameterList& structparams,
      const std::string struct_disname,
      const std::string scatra_disname);

  /*!
  \brief Init vectors and pointers

   Initializes members.

  \note Must only be called after parallel (re-)distribution of discretizations is finished !
        Otherwise, vectors may have wrong maps.

  \warning none
  \return void
  \date 08/16
  \author rauch
  */
  virtual void Setup();

  /// setup of single fields (if needed)
  virtual void SetupSystem();

  /// time loop of coupled problem
  virtual void Timeloop() = 0;

protected:

  ///prepare time step of single fields
  virtual void PrepareTimeStep(bool printheader = true) = 0 ;

  ///do one time step (one way coupled) or one inner iteration loop step (two way coupled), depending on coupling algorithm
  virtual void DoStructStep() = 0 ;

  ///do one time step (one way coupled) or one inner iteration loop step (two way coupled), depending on coupling algorithm
  virtual void DoScatraStep() = 0 ;

};

}

#endif /* SSI_PARTITIONED_H_ */
