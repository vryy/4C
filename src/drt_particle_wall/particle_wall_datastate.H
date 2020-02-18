/*---------------------------------------------------------------------------*/
/*! \file
\brief wall data state container for particle wall handler

\level 2

\maintainer Sebastian Fuchs
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
#ifndef PARTICLE_WALL_DATASTATE_H
#define PARTICLE_WALL_DATASTATE_H

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include <memory>

#include <Epetra_Vector.h>
#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>

/*---------------------------------------------------------------------------*
 | forward declarations                                                      |
 *---------------------------------------------------------------------------*/
namespace IO
{
  class DiscretizationReader;
}

namespace DRT
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

    //! destructor
    ~WallDataState(){};

    /*!
     * \brief init wall data state container
     *
     * \author Sebastian Fuchs \date 05/2019
     *
     * \param[in] walldiscretization wall discretization
     */
    void Init(const Teuchos::RCP<DRT::Discretization> walldiscretization);

    /*!
     * \brief setup wall data state container
     *
     * \author Sebastian Fuchs \date 05/2019
     */
    void Setup();

    /*!
     * \brief write restart of wall data state container
     *
     * \author Sebastian Fuchs \date 09/2019
     *
     * \param[in] step restart step
     * \param[in] time restart time
     */
    void WriteRestart(const int step, const double time) const;

    /*!
     * \brief read restart of wall data state container
     *
     * \author Sebastian Fuchs \date 09/2019
     *
     * \param[in] reader discretization reader
     */
    void ReadRestart(const std::shared_ptr<IO::DiscretizationReader> reader);

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
    void UpdateMapsOfStateVectors();

    //! @name get states (read only access)
    //! @{

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

    //! @}

    //! @name get states (read and write access)
    //! @{

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

    //! @}

   private:
    //! particle simulation parameter list
    const Teuchos::ParameterList& params_;

    //! wall discretization
    Teuchos::RCP<DRT::Discretization> walldiscretization_;

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
#endif
