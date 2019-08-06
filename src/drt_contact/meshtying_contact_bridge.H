/*!----------------------------------------------------------------------
\level 2

\maintainer Matthias Mayr

\brief Bridge for accessing meshtying & contact from STR time integration
*/
/*-----------------------------------------------------------------------*/
#ifndef MESHTYING_CONTACT_BRIDGE_H_
#define MESHTYING_CONTACT_BRIDGE_H_

#include "Teuchos_RCP.hpp"
#include "Epetra_Vector.h"
#include "../drt_inpar/inpar_contact.H"
#include "Epetra_Comm.h"
#include "../drt_lib/drt_dserror.H"

// forward declarations
namespace DRT
{
  class Discretization;
  class Condition;
}  // namespace DRT

namespace IO
{
  class DiscretizationWriter;
  class DiscretizationReader;
}  // namespace IO

namespace LINALG
{
  class MapExtractor;
}  // namespace LINALG

namespace MORTAR
{
  class ManagerBase;
  class StrategyBase;
}  // namespace MORTAR

namespace CONTACT
{
  /*!
  \brief Bridge to enable unified access to contact and meshtying managers

  This bridge wraps contact and meshtying managers, such that the structure time integration does
  not have to distinguish between contact and meshtying operations, but has a single interface to
  both of them. The distinction between contact and meshtying operations is hidden in here.
  */
  class MeshtyingContactBridge
  {
   public:
    /*!
    \brief Constructor

    @param dis Structure discretization
    @param mtcond List of meshtying conditions as given in input file
    @param ccond List of contact conditions as given in input file
    @param timeIntegrationMidPoint Generalized mid-point of time integration scheme
    */
    MeshtyingContactBridge(DRT::Discretization& dis, std::vector<DRT::Condition*>& mtcond,
        std::vector<DRT::Condition*>& ccond, double timeIntegrationMidPoint);

    /*!
    \brief Destructor

    */
    virtual ~MeshtyingContactBridge() = default;

    //! @name Access methods

    /*!
    \brief Get Epetra communicator

    */
    const Epetra_Comm& Comm() const;

    /*!
    \brief Get contact manager

    */
    Teuchos::RCP<MORTAR::ManagerBase> ContactManager() const;

    /*!
    \brief Get meshtying manager

    */
    Teuchos::RCP<MORTAR::ManagerBase> MtManager() const;

    /*!
    \brief Get strategy of mt/contact problem

    */
    MORTAR::StrategyBase& GetStrategy();

    /*!
    \brief return bool indicating if contact is defined

    */
    bool HaveContact() { return (cman_ != Teuchos::null); }

    /*!
    \brief return bool indicating if meshtying is defined

    */
    bool HaveMeshtying() { return (mtman_ != Teuchos::null); }

    /*!
    \brief Write results for visualization for meshtying/contact problems

    This routine does some postprocessing (e.g. computing interface tractions) and then writes
    results to disk through the structure discretization's output writer \c output.

    \param[in] output Output writer of structure discretization to write results to disk
    */
    void PostprocessQuantities(Teuchos::RCP<IO::DiscretizationWriter>& output);

    /*!
    \brief Write results for visualization separately for each meshtying/contact interface

    Call each interface, such that each interface can handle its own output of results.

    \param[in] outputParams Parameter list with stuff required by interfaces to write output
    */
    void PostprocessQuantitiesPerInterface(Teuchos::RCP<Teuchos::ParameterList> outputParams);

    /*!
    \brief read restart

    */
    void ReadRestart(IO::DiscretizationReader& reader, Teuchos::RCP<Epetra_Vector> dis,
        Teuchos::RCP<Epetra_Vector> zero);
    /*!
    \brief recover lagr. mult. for contact/meshtying and slave displ for mesht.

    */
    void Recover(Teuchos::RCP<Epetra_Vector> disi);

    /*!
    \brief set state vector

    */
    void SetState(Teuchos::RCP<Epetra_Vector> zeros);

    /*!
    \brief store dirichlet status

    */
    void StoreDirichletStatus(Teuchos::RCP<LINALG::MapExtractor> dbcmaps);

    /*!
    \brief update

    */
    void Update(Teuchos::RCP<Epetra_Vector> dis);

    /*!
    \brief visualize stuff with gmsh

    */
    void VisualizeGmsh(const int istep, const int iter = -1);

    /*!
    \brief write restart

    @param[in] output Output writer to be used for writing outpu
    @oaram[in] forcedrestart Force to write restart data

    */
    void WriteRestart(Teuchos::RCP<IO::DiscretizationWriter>& output, bool forcedrestart = false);

   private:
    // don't want cctor (= operator impossible anyway for abstract class)
    MeshtyingContactBridge(const MeshtyingContactBridge& old) = delete;

    // internal variables:

    //! Contact manager
    Teuchos::RCP<MORTAR::ManagerBase> cman_;

    //! Meshtying manager
    Teuchos::RCP<MORTAR::ManagerBase> mtman_;

  };  // class MeshtyingContactBridge
}  // namespace CONTACT

#endif /* MESHTYING_CONTACT_BRIDGE_H_ */
