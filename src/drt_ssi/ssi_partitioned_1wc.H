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
    //! constructor
    explicit SSI_Part1WC(const Epetra_Comm& comm,
        const Teuchos::ParameterList& globaltimeparams,
        const Teuchos::ParameterList& scatraparams,
        const Teuchos::ParameterList& structparams,
        const std::string struct_disname,
        const std::string scatra_disname);

    //! virtual destructor
    virtual ~SSI_Part1WC(){};

    //! actual time loop (implemented by derived class)
    virtual void Timeloop()=0;

  protected:
    //! prepare time step for single fields
    void PrepareTimeStep()=0;

    //! perform one time step of structure field
    void DoStructStep();

    //! perform one time step of scatra field
    void DoScatraStep();

    // Flag for reading scatra result from restart files
     bool isscatrafromfile_;

  };

  class SSI_Part1WC_SolidToScatra : public SSI_Part1WC
  {
  public:
    //! constructor
    explicit SSI_Part1WC_SolidToScatra(const Epetra_Comm& comm,
        const Teuchos::ParameterList& globaltimeparams,
        const Teuchos::ParameterList& scatraparams,
        const Teuchos::ParameterList& structparams,
        const std::string struct_disname,
        const std::string scatra_disname);

    //! virtual destructor
    virtual ~SSI_Part1WC_SolidToScatra(){};

    //! actual time loop
    void Timeloop();

    //! prepare time step for single fields
    void PrepareTimeStep();

  };

  class SSI_Part1WC_ScatraToSolid : public SSI_Part1WC
  {
  public:
    //! constructor
    explicit SSI_Part1WC_ScatraToSolid(const Epetra_Comm& comm,
        const Teuchos::ParameterList& globaltimeparams,
        const Teuchos::ParameterList& scatraparams,
        const Teuchos::ParameterList& structparams,
        const std::string struct_disname,
        const std::string scatra_disname) ;

    //! virtual destructor
    virtual ~SSI_Part1WC_ScatraToSolid(){};

    //! actual time loop
    void Timeloop();

    //! prepare time step for single fields
    void PrepareTimeStep();
  };

}

#endif /* SSI_PARTITIONED_1WC_H_ */
