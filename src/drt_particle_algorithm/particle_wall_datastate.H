/*---------------------------------------------------------------------------*/
/*!
\file particle_wall_datastate.H

\brief wall data state container for particle wall handler

\level 3

\maintainer  Sebastian Fuchs
             fuchs@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 -15262

*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                sfuchs 05/2019 |
 *---------------------------------------------------------------------------*/
#ifndef PARTICLE_WALL_DATASTATE_H
#define PARTICLE_WALL_DATASTATE_H

/*---------------------------------------------------------------------------*
 | headers                                                    sfuchs 05/2019 |
 *---------------------------------------------------------------------------*/
#include <Epetra_Vector.h>
#include <Teuchos_RCP.hpp>

/*---------------------------------------------------------------------------*
 | forward declarations                                       sfuchs 05/2019 |
 *---------------------------------------------------------------------------*/
namespace DRT
{
  class Discretization;
}

/*---------------------------------------------------------------------------*
 | class definitions                                          sfuchs 05/2019 |
 *---------------------------------------------------------------------------*/
namespace PARTICLEALGORITHM
{
  class WallDataState
  {
   public:
    //! constructor
    explicit WallDataState();

    //! destructor
    ~WallDataState(){};

    //! init wall data state container
    void Init();

    //! setup wall data state container
    void Setup(Teuchos::RCP<const DRT::Discretization> const& walldiscretization, bool ismoving,
        bool isloaded);

    //! check for correct maps
    void CheckForCorrectMaps(Teuchos::RCP<const DRT::Discretization> const& walldiscretization);

    //! update maps of state vectors
    void UpdateMapsOfStateVectors(
        Teuchos::RCP<const DRT::Discretization> const& walldiscretization);

    //! @name get states (read only access)
    //@{

    //! get wall displacements (row map based)
    inline Teuchos::RCP<const Epetra_Vector> GetDispRow() const { return disp_row_; };

    //! get wall displacements (column map based)
    inline Teuchos::RCP<const Epetra_Vector> GetDispCol() const { return disp_col_; };

    //! get wall displacements (row map based) after last transfer
    inline Teuchos::RCP<const Epetra_Vector> GetDispRowLastTransfer() const
    {
      return disp_row_last_transfer_;
    };

    //! get wall velocities (column map based)
    inline Teuchos::RCP<const Epetra_Vector> GetVelCol() const { return vel_col_; };

    //! get wall accelerations (column map based)
    inline Teuchos::RCP<const Epetra_Vector> GetAccCol() const { return acc_col_; };

    //! get wall forces (column map based)
    inline Teuchos::RCP<const Epetra_Vector> GetForceCol() const { return force_col_; };

    //@}

    //! @name get states (read and write access)
    //@{

    //! get wall displacements (row map based)
    inline Teuchos::RCP<Epetra_Vector> GetMutableDispRow() { return disp_row_; };
    inline Teuchos::RCP<Epetra_Vector>& GetRefMutableDispRow() { return disp_row_; };

    //! get wall displacements (column map based)
    inline Teuchos::RCP<Epetra_Vector> GetMutableDispCol() { return disp_col_; };
    inline Teuchos::RCP<Epetra_Vector>& GetRefMutableDispCol() { return disp_col_; };

    //! get wall displacements (row map based) after last transfer
    inline Teuchos::RCP<Epetra_Vector> GetMutableDispRowLastTransfer()
    {
      return disp_row_last_transfer_;
    };

    //! get wall velocities (column map based)
    inline Teuchos::RCP<Epetra_Vector> GetMutableVelCol() { return vel_col_; };

    //! get wall accelerations (column map based)
    inline Teuchos::RCP<Epetra_Vector> GetMutableAccCol() { return acc_col_; };

    //! get wall forces (column map based)
    inline Teuchos::RCP<Epetra_Vector> GetMutableForceCol() { return force_col_; };

    //@}

   private:
    //! current dof row map
    Teuchos::RCP<Epetra_Map> curr_dof_row_map_;

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
  };

}  // namespace PARTICLEALGORITHM

/*---------------------------------------------------------------------------*/
#endif
