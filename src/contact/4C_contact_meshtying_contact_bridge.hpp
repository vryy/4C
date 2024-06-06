/*----------------------------------------------------------------------*/
/*! \file
\level 2


\brief Bridge for accessing meshtying & contact from STR time integration
*/
/*-----------------------------------------------------------------------*/
#ifndef FOUR_C_CONTACT_MESHTYING_CONTACT_BRIDGE_HPP
#define FOUR_C_CONTACT_MESHTYING_CONTACT_BRIDGE_HPP

#include "4C_config.hpp"

#include "4C_discretization_condition.hpp"
#include "4C_inpar_contact.hpp"
#include "4C_utils_exceptions.hpp"

#include <Epetra_Comm.h>
#include <Epetra_Vector.h>
#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN

// forward declarations
namespace Discret
{
  class Discretization;
}  // namespace Discret

namespace Core::IO
{
  class DiscretizationWriter;
  class DiscretizationReader;
}  // namespace Core::IO

namespace Core::LinAlg
{
  class MapExtractor;
}  // namespace Core::LinAlg

namespace Mortar
{
  class ManagerBase;
  class StrategyBase;
}  // namespace Mortar

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
    @param meshtyingConditions List of meshtying conditions as given in input file
    @param contactConditions List of contact conditions as given in input file
    @param timeIntegrationMidPoint Generalized mid-point of time integration scheme
    */
    MeshtyingContactBridge(Discret::Discretization& dis,
        std::vector<Core::Conditions::Condition*>& meshtyingConditions,
        std::vector<Core::Conditions::Condition*>& contactConditions,
        double timeIntegrationMidPoint);

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
    Teuchos::RCP<Mortar::ManagerBase> ContactManager() const;

    /*!
    \brief Get meshtying manager

    */
    Teuchos::RCP<Mortar::ManagerBase> MtManager() const;

    /*!
    \brief Get strategy of meshtying/contact problem

    */
    Mortar::StrategyBase& GetStrategy() const;

    /*!
    \brief return bool indicating if contact is defined

    */
    bool HaveContact() const { return (cman_ != Teuchos::null); }

    /*!
    \brief return bool indicating if meshtying is defined

    */
    bool HaveMeshtying() const { return (mtman_ != Teuchos::null); }

    /*!
    \brief Write results for visualization for meshtying/contact problems

    This routine does some postprocessing (e.g. computing interface tractions) and then writes
    results to disk through the structure discretization's output writer \c output.

    \param[in] output Output writer of structure discretization to write results to disk
    */
    void postprocess_quantities(Teuchos::RCP<Core::IO::DiscretizationWriter>& output);

    /*!
    \brief Write results for visualization separately for each meshtying/contact interface

    Call each interface, such that each interface can handle its own output of results.

    \param[in] outputParams Parameter list with stuff required by interfaces to write output
    */
    void postprocess_quantities_per_interface(Teuchos::RCP<Teuchos::ParameterList> outputParams);

    /*!
    \brief read restart

    */
    void read_restart(Core::IO::DiscretizationReader& reader, Teuchos::RCP<Epetra_Vector> dis,
        Teuchos::RCP<Epetra_Vector> zero);
    /*!
    \brief recover lagr. mult. for contact/meshtying and slave displ for mesht.

    */
    void Recover(Teuchos::RCP<Epetra_Vector> disi);

    /*!
    \brief set state vector

    */
    void set_state(Teuchos::RCP<Epetra_Vector> zeros);

    /*!
    \brief store dirichlet status

    */
    void store_dirichlet_status(Teuchos::RCP<Core::LinAlg::MapExtractor> dbcmaps);

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
    void write_restart(
        Teuchos::RCP<Core::IO::DiscretizationWriter>& output, bool forcedrestart = false);

   private:
    //! don't want cctor (= operator impossible anyway for abstract class)
    MeshtyingContactBridge(const MeshtyingContactBridge& old) = delete;

    //! Contact manager
    Teuchos::RCP<Mortar::ManagerBase> cman_;

    //! Meshtying manager
    Teuchos::RCP<Mortar::ManagerBase> mtman_;

  };  // class meshtying_contact_bridge
}  // namespace CONTACT

FOUR_C_NAMESPACE_CLOSE

#endif
