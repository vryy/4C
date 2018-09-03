/*!----------------------------------------------------------------------
\file immersed_partitioned_cellmigration.H
\level 2

\brief partitioned immersed cell migration algorithm

\maintainer  Andreas Rauch
             rauch@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 -15240

*----------------------------------------------------------------------*/
#ifndef IMMERSED_PARTITIONED_CELLMIGRATION_H
#define IMMERSED_PARTITIONED_CELLMIGRATION_H

#include "immersed_partitioned.H"
#include "immersed_field_exchange_manager.H"

/*----------------------------------------------------------------------*
 | forward declarations                                                  |
 *----------------------------------------------------------------------*/


/*----------------------------------------------------------------------*/

namespace IMMERSED
{
  class ImmersedPartitionedCellMigration : public ImmersedPartitioned
  {
   public:
    /// setup partitioned immersed cell migration algorithm
    explicit ImmersedPartitionedCellMigration(
        const Teuchos::ParameterList& params, const Epetra_Comm& comm);

    /*! \brief Initialize this object

    \warning none
    \return int
    \date 08/16
    \author rauch  */
    virtual int Init(const Teuchos::ParameterList& params);

    /*! \brief Setup all class internal objects and members

    \warning none
    \return void
    \date 08/16
    \author rauch  */
    virtual void Setup();

   private:
    /// pointer to global problem
    DRT::Problem* globalproblem_;

    //! @name Manager for data exchange between fields
    //@{
    DRT::ImmersedFieldExchangeManager* exchange_manager_;
    //@}

    //! @name Interaction Switches
    //@{
    bool fluid_interaction_;  //!< switch for cell-fluid-interaction
    bool ecm_interaction_;    //!< switch for cell compression by ecm
    bool adhesion_dynamics_;  //!< switch for cell-ecm interaction at adhesion sites
    //@}

  };  // class ImmersedPartitionedCellMigration
}  // namespace IMMERSED

#endif
