/*----------------------------------------------------------------------*/
/*! \file

\brief Methods to transfer a turbulent inflow profile a from a master
boundary to a slave boundary on a seperate domain.
The slave boundary must have an additional Dirichlet condition, the
master boundary will usually be a periodic boundary (but is not required
to)


\level 2

*/
/*----------------------------------------------------------------------*/

#ifndef FOUR_C_FLUID_TURBULENCE_TRANSFER_TURB_INFLOW_HPP
#define FOUR_C_FLUID_TURBULENCE_TRANSFER_TURB_INFLOW_HPP

#include "4C_config.hpp"

#include "4C_comm_exporter.hpp"
#include "4C_linalg_utils_sparse_algebra_math.hpp"

#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN

namespace Core::FE
{
  class Discretization;
}  // namespace Core::FE

namespace FLD
{
  class TransferTurbulentInflowCondition
  {
   public:
    /*!
    \brief Standard Constructor

    */
    TransferTurbulentInflowCondition(Teuchos::RCP<Core::FE::Discretization> dis,
        Teuchos::RCP<Core::LinAlg::MapExtractor> dbcmaps);

    /*!
    \brief Destructor

    */
    virtual ~TransferTurbulentInflowCondition() = default;

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
    void get_data(int& id, int& direction, ToggleType& type, const Core::Conditions::Condition*);

    //! receive a block in the round robin communication pattern
    void receive_block(
        std::vector<char>& rblock, Core::Communication::Exporter& exporter, MPI_Request& request);

    //! send a block in the round robin communication pattern
    void send_block(
        std::vector<char>& sblock, Core::Communication::Exporter& exporter, MPI_Request& request);

    //! unpack all master values contained in receive block
    void unpack_local_master_values(std::vector<int>& mymasters,
        std::vector<std::vector<double>>& mymasters_vel, std::vector<char>& rblock);

    //! pack all master values into a send block
    void pack_local_master_values(std::vector<int>& mymasters,
        std::vector<std::vector<double>>& mymasters_vel, Core::Communication::PackBuffer& sblock);

    //! for all values avaible on the processor, do the final setting of the value
    virtual void set_values_available_on_this_proc(std::vector<int>& mymasters,
        std::vector<std::vector<double>>& mymasters_vel, Teuchos::RCP<Epetra_Vector> velnp);

    //! flag active boundary condition (may be used to switch off everything)
    bool active_;

    //! the discretisation
    Teuchos::RCP<Core::FE::Discretization> dis_;

    //! info on DIirchlet boundary
    Teuchos::RCP<Core::LinAlg::MapExtractor> dbcmaps_;

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
    TransferTurbulentInflowConditionXW(Teuchos::RCP<Core::FE::Discretization> dis,
        Teuchos::RCP<Core::LinAlg::MapExtractor> dbcmaps);


    /*!
    \brief Transfer process copying values from master boundary to slave
           boundary (slave must be of Dirichlet type, otherwise this
           operation doesn't make to mauch sense)

           Intended to be called after ApplyDirichlet, overwriting the
           dummy Dirichlet values on the slave boundary by the values
           of the last time step on the master boundary
    */
    void Transfer(const Teuchos::RCP<Epetra_Vector> veln, Teuchos::RCP<Epetra_Vector> velnp,
        const double time) override;

   private:
    //! for all values avaible on the processor, do the final setting of the value
    void set_values_available_on_this_proc(std::vector<int>& mymasters,
        std::vector<std::vector<double>>& mymasters_vel,
        Teuchos::RCP<Epetra_Vector> velnp) override;
  };


  class TransferTurbulentInflowConditionNodal : public TransferTurbulentInflowCondition
  {
   public:
    /*!
    \brief Standard Constructor

    */
    TransferTurbulentInflowConditionNodal(Teuchos::RCP<Core::FE::Discretization> dis,
        Teuchos::RCP<Core::LinAlg::MapExtractor> dbcmaps);


    /*!
    \brief Transfer process copying values from master boundary to slave
           boundary (slave must be of Dirichlet type, otherwise this
           operation doesn't make to mauch sense)

           Intended to be called after ApplyDirichlet, overwriting the
           dummy Dirichlet values on the slave boundary by the values
           of the last time step on the master boundary
    */
    void Transfer(const Teuchos::RCP<Epetra_Vector> invec, Teuchos::RCP<Epetra_Vector> outvec,
        const double time) override;

    bool is_active() { return active_; }

   private:
    //! for all values avaible on the processor, do the final setting of the value
    void set_values_available_on_this_proc(std::vector<int>& mymasters,
        std::vector<std::vector<double>>& mymasters_vec,
        Teuchos::RCP<Epetra_Vector> outvec) override;
  };

}  // end namespace FLD

FOUR_C_NAMESPACE_CLOSE

#endif
