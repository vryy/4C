/*----------------------------------------------------------------------------*/
/*!

 \brief FSI Wrapper for the ALE time integration

 \level 1

 \maintainer Matthias Mayr
 */
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/

#ifndef AD_ALE_FSI_H_
#define AD_ALE_FSI_H_


/*----------------------------------------------------------------------------*/
/* header inclusions */
#include "ad_ale_wrapper.H"

/*----------------------------------------------------------------------------*/
/* forward declarations */
namespace ALE
{
  namespace UTILS
  {
    class MapExtractor;
  }
}  // namespace ALE

/*----------------------------------------------------------------------------*/
/* class definitions */
namespace ADAPTER
{
  /*! \brief ALE Wrapper for FSI Problems
   *
   *  Provide FSI specific ALE functionalities here by overloading the respective
   *  routines from ADAPTER::AleWrapper
   *
   *  \sa ADAPTER::Ale, ADAPTER::AleWrapper
   *
   *  \author mayr.mt \date 10/2014
   */
  class AleFsiWrapper : public AleWrapper
  {
   public:
    //! @name Construction / Destruction
    //@{

    //! constructor
    explicit AleFsiWrapper(Teuchos::RCP<Ale> ale);

    //! Destructor
    virtual ~AleFsiWrapper() {}

    //@}

    //! communicate object at the interface
    Teuchos::RCP<const ALE::UTILS::MapExtractor> Interface() const;

    //! apply interface displacements
    void ApplyInterfaceDisplacements(Teuchos::RCP<const Epetra_Vector> idisp)
    {
      interface_->InsertFSICondVector(idisp, WriteAccessDispnp());
    }

    //! get Dirichlet map extractor
    virtual const Teuchos::RCP<const LINALG::MapExtractor> GetDBCMapExtractor()
    {
      return AleWrapper::GetDBCMapExtractor();
    }

   private:
    Teuchos::RCP<ALE::UTILS::MapExtractor> interface_;

  };  // class AleFsiWrapper
}  // namespace ADAPTER

#endif /* AD_ALE_FSI_H_ */
