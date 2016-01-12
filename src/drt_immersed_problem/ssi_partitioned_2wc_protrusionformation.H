/*!----------------------------------------------------------------------
\file ssi_partitioned_2wc_protrusionformation.H

\brief specialization of ssi2wc, which includes "structale"-surface growth

<pre>
Maintainers: Andreas Rauch
             rauch@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 -15240
</pre>
*----------------------------------------------------------------------*/
#ifndef SSI_PARTITIONED2_WC_PROTRUSIONFORMATION_H_
#define SSI_PARTITIONED2_WC__PROTRUSIONFORMATION_H_

#include "../drt_ssi/ssi_partitioned_2wc.H"
#include "../drt_inpar/inpar_ssi.H"

namespace SSI
{
  //! specialized class to deal with partitioned 2WC SSI during protrusion formation in cell migration.
  class SSI_Part2WC_PROTRUSIONFORMATION : public SSI_Part2WC
  {
  public:

    //!explicit constructor
    explicit SSI_Part2WC_PROTRUSIONFORMATION(const Epetra_Comm& comm,
        const Teuchos::ParameterList& globaltimeparams,
        const Teuchos::ParameterList& scatraparams,
        const Teuchos::ParameterList& structparams,
        const std::string struct_disname,
        const std::string scatra_disname);

    //!virtual destructor
    virtual ~SSI_Part2WC_PROTRUSIONFORMATION(){};

  }; // class
} // namespace SSI

#endif /* SSI_PARTITIONED2_WC__PROTRUSIONFORMATION_H_ */
