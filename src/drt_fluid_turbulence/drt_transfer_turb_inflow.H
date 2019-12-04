/*----------------------------------------------------------------------*/
/*! \file

\brief Methods to transfer a turbulent inflow profile a from a master
boundary to a slave boundary on a seperate domain.
The slave boundary must have an additional Dirichlet condition, the
master boundary will usually be a periodic boundary (but is not required
to)

\maintainer Martin Kronbichler

\level 2

*/
/*----------------------------------------------------------------------*/

#ifndef TRANSFERTURBINF_H
#define TRANSFERTURBINF_H
#include "Teuchos_RCP.hpp"
#include "../linalg/linalg_utils_sparse_algebra_math.H"

/*

*/

// forward declarations
namespace DRT
{
  class Discretization;
}

namespace FLD
{
  class TransferTurbulentInflowCondition
  {
   public:
    /*!
    \brief Standard Constructor

    */
    TransferTurbulentInflowCondition(
        Teuchos::RCP<DRT::Discretization> dis, Teuchos::RCP<LINALG::MapExtractor> dbcmaps);

    /*!
    \brief Destructor

    */
    virtual ~TransferTurbulentInflowCondition();

    /*!
    \brief Transfer process copying values from master boundary to slave
           boundary (slave must be of Dirichlet type, otherwise this
           operation doesn't make to mauch sense)

           Intended to be called after ApplyDirichlet, overwriting the
           dummy Dirichlet values on the slave boundary by the values
           of the last time step on the master boundary
    */
    virtual void Transfer(const Teuchos::RCP<Epetra_Vector> veln, Teuchos::RCP<Epetra_Vector> velnp,
        const double time);

   protected:
    //! there are two types of tranfer conditions. values are transferred
    //! from master to slave conditions
    enum ToggleType
    {
      none,
      master,
      slave
    };

    //! get all data from condtion
    void GetData(int& id, int& direction, ToggleType& type, const DRT::Condition*);

#ifdef PARALLEL
    //! receive a block in the round robin communication pattern
    void ReceiveBlock(std::vector<char>& rblock, DRT::Exporter& exporter, MPI_Request& request);

    //! send a block in the round robin communication pattern
    void SendBlock(std::vector<char>& sblock, DRT::Exporter& exporter, MPI_Request& request);
#endif

    //! unpack all master values contained in receive block
    void UnpackLocalMasterValues(std::vector<int>& mymasters,
        std::vector<std::vector<double>>& mymasters_vel, std::vector<char>& rblock);

    //! pack all master values into a send block
    void PackLocalMasterValues(std::vector<int>& mymasters,
        std::vector<std::vector<double>>& mymasters_vel, DRT::PackBuffer& sblock);

    //! for all values avaible on the processor, do the final setting of the value
    virtual void SetValuesAvailableOnThisProc(std::vector<int>& mymasters,
        std::vector<std::vector<double>>& mymasters_vel, Teuchos::RCP<Epetra_Vector> velnp);

    //! flag active boundary condition (may be used to switch off everything)
    bool active_;

    //! the discretisation
    Teuchos::RCP<DRT::Discretization> dis_;

    //! info on DIirchlet boundary
    Teuchos::RCP<LINALG::MapExtractor> dbcmaps_;

    //! the connectivity of the boundary condition
    std::map<int, std::vector<int>> midtosid_;

    //! time curve number
    int curve_;

    int numveldof_;
  };

  class TransferTurbulentInflowConditionXW : public TransferTurbulentInflowCondition
  {
   public:
    /*!
    \brief Standard Constructor

    */
    TransferTurbulentInflowConditionXW(
        Teuchos::RCP<DRT::Discretization> dis, Teuchos::RCP<LINALG::MapExtractor> dbcmaps);

    /*!
    \brief Destructor

    */
    virtual ~TransferTurbulentInflowConditionXW();

    /*!
    \brief Transfer process copying values from master boundary to slave
           boundary (slave must be of Dirichlet type, otherwise this
           operation doesn't make to mauch sense)

           Intended to be called after ApplyDirichlet, overwriting the
           dummy Dirichlet values on the slave boundary by the values
           of the last time step on the master boundary
    */
    void Transfer(const Teuchos::RCP<Epetra_Vector> veln, Teuchos::RCP<Epetra_Vector> velnp,
        const double time);

   private:
    //! for all values avaible on the processor, do the final setting of the value
    void SetValuesAvailableOnThisProc(std::vector<int>& mymasters,
        std::vector<std::vector<double>>& mymasters_vel, Teuchos::RCP<Epetra_Vector> velnp);
  };


  class TransferTurbulentInflowConditionNodal : public TransferTurbulentInflowCondition
  {
   public:
    /*!
    \brief Standard Constructor

    */
    TransferTurbulentInflowConditionNodal(
        Teuchos::RCP<DRT::Discretization> dis, Teuchos::RCP<LINALG::MapExtractor> dbcmaps);

    /*!
    \brief Destructor

    */
    virtual ~TransferTurbulentInflowConditionNodal();

    /*!
    \brief Transfer process copying values from master boundary to slave
           boundary (slave must be of Dirichlet type, otherwise this
           operation doesn't make to mauch sense)

           Intended to be called after ApplyDirichlet, overwriting the
           dummy Dirichlet values on the slave boundary by the values
           of the last time step on the master boundary
    */
    void Transfer(const Teuchos::RCP<Epetra_Vector> invec, Teuchos::RCP<Epetra_Vector> outvec,
        const double time);

    bool IsActive() { return active_; }

   private:
    //! for all values avaible on the processor, do the final setting of the value
    void SetValuesAvailableOnThisProc(std::vector<int>& mymasters,
        std::vector<std::vector<double>>& mymasters_vec, Teuchos::RCP<Epetra_Vector> outvec);
  };

}  // end namespace FLD

#endif  // #ifndef TRANSFERTURBINF_H
