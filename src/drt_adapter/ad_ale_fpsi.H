/*----------------------------------------------------------------------------*/
/*! \file

 \brief FPSI wrapper for the ALE time integration

 \level 2

\maintainer  Johannes Kremheller
 */
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/

#ifndef AD_ALE_FPSI_H_
#define AD_ALE_FPSI_H_


/*----------------------------------------------------------------------------*/
/* header inclusions */
#include "ad_ale_wrapper.H"

/*----------------------------------------------------------------------------*/
/* class definitions */
namespace ADAPTER
{
  class AleFpsiWrapper : public AleWrapper
  {
   public:
    //! @name Construction / Destruction
    //@{

    //! constructor
    explicit AleFpsiWrapper(Teuchos::RCP<Ale> ale);

    //! Destructor
    virtual ~AleFpsiWrapper() {}

    //! specialized method to apply displacements to fpsi interface
    void ApplyInterfaceDisplacements(Teuchos::RCP<const Epetra_Vector> idisp);

    //! specialized method to apply displacements to fsi interface
    void ApplyFSIInterfaceDisplacements(Teuchos::RCP<const Epetra_Vector> idisp);

    //! communicate object at the interface
    Teuchos::RCP<const ALE::UTILS::MapExtractor> Interface() const;

    //@}

   private:
    //! interface map extractor
    Teuchos::RCP<ALE::UTILS::MapExtractor> interface_;


  };  // class AleFpsiWrapper
}  // namespace ADAPTER

#endif /* AD_ALE_FPSI_H_ */
