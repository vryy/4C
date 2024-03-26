/*---------------------------------------------------------------------*/
/*! \file

\brief A collection of helper methods for namespace DRT

\level 0


*/
/*---------------------------------------------------------------------*/

#ifndef BACI_LIB_UTILS_HPP
#define BACI_LIB_UTILS_HPP

#include "baci_config.hpp"

#include "baci_lib_element.hpp"
#include "baci_lib_node.hpp"
#include "baci_linalg_serialdensematrix.hpp"
#include "baci_linalg_serialdensevector.hpp"

#include <signal.h>
#include <stdio.h>
#include <Teuchos_RCP.hpp>

#include <random>

BACI_NAMESPACE_OPEN

// forward declarations
namespace CORE::LINALG
{
  class SparseMatrix;
}

namespace DRT
{
  namespace UTILS
  {
    /*!
      \brief compute L2 projection of a dof based field onto a node based field in a least
      squares sense.
      WARNING: Make sure to pass down a dofrowmap appropriate for your discretization.

      \return an Epetra_MultiVector based on the discret's node row map containing numvec vectors
              with the projected state

      \author Georg Hammerl
      \date 06/14
     */
    Teuchos::RCP<Epetra_MultiVector> ComputeNodalL2Projection(
        Teuchos::RCP<DRT::Discretization> dis,  ///< underlying discretization
        const std::string& statename,           ///< name of state which will be set
        const int& numvec,                      ///< number of entries per node to project
        Teuchos::ParameterList& params,         ///< parameter list that contains the element action
        const Teuchos::ParameterList&
            solverparams);  ///< solver parameters for solving the resulting global system;

    Teuchos::RCP<Epetra_MultiVector> EvaluateAndSolveNodalL2Projection(Discretization& dis,
        const Epetra_Map& noderowmap, const std::string& statename, const int& numvec,
        Teuchos::ParameterList& params, const Teuchos::ParameterList& solverparams,
        const Epetra_Map* fullnoderowmap = nullptr,
        const std::map<int, int>* slavetomastercolnodesmap = nullptr);

    Teuchos::RCP<Epetra_MultiVector> SolveNodalL2Projection(CORE::LINALG::SparseMatrix& massmatrix,
        Epetra_MultiVector& rhs, const Epetra_Comm& comm, const int& numvec,
        const Teuchos::ParameterList& solverparams, const Epetra_Map& noderowmap,
        const Epetra_Map* fullnoderowmap = nullptr,
        const std::map<int, int>* slavetomastercolnodesmap = nullptr);

    /*!
      \brief reconstruct nodal values via superconvergent patch recovery

      \return an Epetra_MultiVector based on the discret's node row map containing numvec vectors
              with the reconstruced state

      \author Georg Hammerl
      \date 05/15
     */
    template <int dim>
    Teuchos::RCP<Epetra_MultiVector> ComputeSuperconvergentPatchRecovery(
        Teuchos::RCP<DRT::Discretization> dis,    ///< underlying discretization
        Teuchos::RCP<const Epetra_Vector> state,  ///< state vector needed on element level
        const std::string statename,              ///< name of state which will be set
        const int numvec,                         ///< number of entries per node to project
        Teuchos::ParameterList& params  ///< parameter list that contains the element action
    );


    /*!
      \brief Return Element center coordinates
    */
    std::vector<double> ElementCenterRefeCoords(const DRT::Element* const ele);

    /*!
    \brief handles restart after a certain walltime interval, step interval or on a user signal

    \author hammerl
    */
    class RestartManager
    {
     public:
      RestartManager();

      virtual ~RestartManager() = default;

      /// setup of restart manager
      void SetupRestartManager(const double restartinterval, const int restartevry);

      /// return whether it is time for a restart
      /// \param step [in] : current time step for multi-field syncronisation
      /// \param comm [in] : get access to involved procs
      bool Restart(const int step, const Epetra_Comm& comm);

      /// the signal handler that gets passed to the kernel and listens for SIGUSR1 and SIGUSR2
      static void restart_signal_handler(
          int signal_number, siginfo_t* signal_information, void* ignored);

     protected:
      /// @name wall time parameters
      //@{

      /// start time of simulation
      double startwalltime_;

      /// after this wall time interval a restart is enforced
      double restartevrytime_;

      /// check to enforce restart only once during time interval
      int restartcounter_;

      //@}

      /// store the step which was allowed to write restart
      int lastacceptedstep_;

      /// member to detect time step increment
      int lasttestedstep_;

      /// after this number of steps a restart is enforced
      int restartevrystep_;

      /// signal which was caught by the signal handler
      volatile static int signal_;
    };

    /**
     * \brief Default error handling of scanf().
     *
     * \param output (in): output provided by the call of scanf()
     * \throws dserror() occurs if the function returns without reading any
     */
    void Checkscanf(int output);
  }  // namespace UTILS
}  // namespace DRT


BACI_NAMESPACE_CLOSE

#endif  // LIB_UTILS_H
