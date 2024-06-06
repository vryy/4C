/*---------------------------------------------------------------------------*/
/*! \file
\brief wall data state container for particle wall handler
\level 2
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
#ifndef FOUR_C_PARTICLE_WALL_DATASTATE_HPP
#define FOUR_C_PARTICLE_WALL_DATASTATE_HPP

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include "4C_config.hpp"

#include <Epetra_Vector.h>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>

#include <memory>

FOUR_C_NAMESPACE_OPEN

/*---------------------------------------------------------------------------*
 | forward declarations                                                      |
 *---------------------------------------------------------------------------*/
namespace Discret
{
  class Discretization;
}

/*---------------------------------------------------------------------------*
 | class declarations                                                        |
 *---------------------------------------------------------------------------*/
namespace PARTICLEWALL
{
  /*!
   * \brief wall data state container
   *
   * \author Sebastian Fuchs \date 05/2019
   */
  class WallDataState final
  {
   public:
    /*!
     * \brief constructor
     *
     * \author Sebastian Fuchs \date 05/2019
     *
     * \param[in] params particle simulation parameter list
     */
    explicit WallDataState(const Teuchos::ParameterList& params);

    /*!
     * \brief init wall data state container
     *
     * \author Sebastian Fuchs \date 05/2019
     *
     * \param[in] walldiscretization wall discretization
     */
    void Init(const Teuchos::RCP<Discret::Discretization> walldiscretization);

    /*!
     * \brief setup wall data state container
     *
     * \author Sebastian Fuchs \date 05/2019
     */
    void Setup();

    /*!
     * \brief check for correct maps
     *
     * \author Sebastian Fuchs \date 05/2019
     */
    void CheckForCorrectMaps();

    /*!
     * \brief update maps of state vectors
     *
     * \author Sebastian Fuchs \date 05/2019
     */
    void update_maps_of_state_vectors();

    //! @name get states (read only access)
    //! @{

    //! get wall displacements (row map based)
    inline Teuchos::RCP<const Epetra_Vector> GetDispRow() const { return disp_row_; };

    //! get wall displacements (column map based)
    inline Teuchos::RCP<const Epetra_Vector> GetDispCol() const { return disp_col_; };

    //! get wall displacements (row map based) after last transfer
    inline Teuchos::RCP<const Epetra_Vector> get_disp_row_last_transfer() const
    {
      return disp_row_last_transfer_;
    };

    //! get wall velocities (column map based)
    inline Teuchos::RCP<const Epetra_Vector> GetVelCol() const { return vel_col_; };

    //! get wall accelerations (column map based)
    inline Teuchos::RCP<const Epetra_Vector> GetAccCol() const { return acc_col_; };

    //! get wall forces (column map based)
    inline Teuchos::RCP<const Epetra_Vector> GetForceCol() const { return force_col_; };

    //! @}

    //! @name get states (read and write access)
    //! @{

    //! get wall displacements (row map based)
    inline Teuchos::RCP<Epetra_Vector> GetDispRow() { return disp_row_; };
    inline Teuchos::RCP<Epetra_Vector>& GetRefDispRow() { return disp_row_; };

    //! get wall displacements (column map based)
    inline Teuchos::RCP<Epetra_Vector> GetDispCol() { return disp_col_; };
    inline Teuchos::RCP<Epetra_Vector>& GetRefDispCol() { return disp_col_; };

    //! get wall displacements (row map based) after last transfer
    inline Teuchos::RCP<Epetra_Vector> get_disp_row_last_transfer()
    {
      return disp_row_last_transfer_;
    };

    //! get wall velocities (column map based)
    inline Teuchos::RCP<Epetra_Vector> GetVelCol() { return vel_col_; };

    //! get wall accelerations (column map based)
    inline Teuchos::RCP<Epetra_Vector> GetAccCol() { return acc_col_; };

    //! get wall forces (column map based)
    inline Teuchos::RCP<Epetra_Vector> GetForceCol() { return force_col_; };

    //! @}

   private:
    //! particle simulation parameter list
    const Teuchos::ParameterList& params_;

    //! wall discretization
    Teuchos::RCP<Discret::Discretization> walldiscretization_;

    //! current dof row map
    Teuchos::RCP<Epetra_Map> curr_dof_row_map_;

    //! @name stored states
    //! @{

    //! wall displacements (row map based)
    Teuchos::RCP<Epetra_Vector> disp_row_;

    //! wall displacements (column map based)
    Teuchos::RCP<Epetra_Vector> disp_col_;

    //! wall displacements (row map based) after last transfer
    Teuchos::RCP<Epetra_Vector> disp_row_last_transfer_;

    //! wall velocities (column map based)
    Teuchos::RCP<Epetra_Vector> vel_col_;

    //! wall accelerations (column map based)
    Teuchos::RCP<Epetra_Vector> acc_col_;

    //! wall forces (column map based)
    Teuchos::RCP<Epetra_Vector> force_col_;

    //! @}
  };

}  // namespace PARTICLEWALL

/*---------------------------------------------------------------------------*/
FOUR_C_NAMESPACE_CLOSE

#endif
