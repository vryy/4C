/*!----------------------------------------------------------------------
\file drt_utils.H

\brief A collection of helper methods for namespace DRT

<pre>
\level 0

\maintainer Martin Kronbichler
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>

*----------------------------------------------------------------------*/
#ifndef UTILS_H
#define UTILS_H

#include <signal.h>
#include <boost/random.hpp>
#include <stdio.h>

#include <Teuchos_RCP.hpp>
#include <Epetra_SerialDenseMatrix.h>
#include <Epetra_SerialDenseVector.h>

#include "drt_node.H"
#include "drt_element.H"

#include "../drt_inpar/inpar_scatra.H"

namespace LINALG
{
  class SparseMatrix;
}

namespace DRT
{
  /// collection of useful and important but less central stuff
  namespace UTILS
  {
    /*!
    \brief Locally extract a subset of values from an Epetra_Vector

    Extracts lm.size() values from a distributed epetra vector and stores them into local.
    this is NOT a parallel method, meaning that all values to be extracted on a processor
    must be present in global on that specific processor. This usually means that global
    has to be in column map style.

    \param global (in): global distributed vector with values to be extracted
    \param local (out): vector or matrix holding values extracted from global
    \param lm     (in): vector containing global ids to be extracted. Size of lm
                        determines number of values to be extracted.
    */
    void ExtractMyValues(
        const Epetra_Vector& global, std::vector<double>& local, const std::vector<int>& lm);

    void ExtractMyValues(
        const Epetra_Vector& global, Epetra_SerialDenseVector& local, const std::vector<int>& lm);

    void ExtractMyValues(
        const Epetra_MultiVector& global, std::vector<double>& local, const std::vector<int>& lm);

    template <class matrix>
    void ExtractMyValues(
        const Epetra_Vector& global, std::vector<matrix>& local, const std::vector<int>& lm)
    {
      // safety check
      if (local[0].N() != 1 or local.size() * (unsigned)local[0].M() != lm.size())
        dserror("Received matrix vector of wrong size!");

      // loop over all nodes of current element
      for (unsigned inode = 0; inode < local[0].M(); ++inode)
      {
        // loop over all dofs of current node
        for (unsigned idof = 0; idof < local.size(); ++idof)
        {
          // extract local ID of current dof
          const int lid = global.Map().LID(lm[inode * local.size() + idof]);

          // safety check
          if (lid < 0)
            dserror("Proc %d: Cannot find gid=%d in Epetra_Vector", global.Comm().MyPID(),
                lm[inode * local.size() + idof]);

          // store current dof in local matrix vector consisting of ndof matrices of size nnode x 1,
          // where nnode denotes the number of element nodes and ndof denotes the number of degrees
          // of freedom per element node.
          local[idof](inode, 0) = global[lid];
        }
      }

      return;
    }

    template <class matrix>
    void ExtractMyValues(const Epetra_Vector& global, matrix& local, const std::vector<int>& lm)
    {
      // safety check
      if ((unsigned)(local.M() * local.N()) != lm.size()) dserror("Received matrix of wrong size!");

      // loop over all columns of local matrix
      for (unsigned icol = 0; icol < local.N(); ++icol)
      {
        // loop over all rows of local matrix
        for (unsigned irow = 0; irow < local.M(); ++irow)
        {
          // extract local ID of current dof
          const unsigned index = icol * local.M() + irow;
          const int lid = global.Map().LID(lm[index]);

          // safety check
          if (lid < 0)
            dserror(
                "Proc %d: Cannot find gid=%d in Epetra_Vector", global.Comm().MyPID(), lm[index]);

          // store current dof in local matrix, which is filled column-wise with the dofs listed in
          // the lm vector
          local(irow, icol) = global[lid];
        }
      }

      return;
    }

    /// Locally extract a subset of values from a (column)-nodemap-based Epetra_MultiVector
    /*  \author henke
     *  \date 06/09
     */
    void ExtractMyNodeBasedValues(const DRT::Element* ele,  ///< pointer to current element
        std::vector<double>& local,                         ///< local vector on element-level
        const Epetra_MultiVector& global                    ///< global (multi) vector
    );


    /// Locally extract a subset of values from a (column)-nodemap-based Epetra_MultiVector
    /*  \author g.bau
     *  \date 08/08
     */
    void ExtractMyNodeBasedValues(const DRT::Element* ele,  ///< pointer to current element
        Epetra_SerialDenseVector& local,                    ///< local vector on element-level
        const Teuchos::RCP<Epetra_MultiVector>& global,     ///< global vector
        const int nsd                                       ///< number of space dimensions
    );

    /// Locally extract a subset of values from a (column)-nodemap-based Epetra_MultiVector
    /*  \author schott
     *  \date 12/16
     */
    void ExtractMyNodeBasedValues(const DRT::Node* node,  ///< pointer to current element
        Epetra_SerialDenseVector& local,                  ///< local vector on node-level
        const Teuchos::RCP<Epetra_MultiVector>& global,   ///< global vector
        const int nsd                                     ///< number of space dimensions
    );


    /// Locally extract a subset of values from a (column)-nodemap-based Epetra_MultiVector
    /// and fill a local matrix that has implemented the (.,.) operator
    /*  \author g.bau
     *  \date 04/09
     */
    template <class M>
    void ExtractMyNodeBasedValues(const DRT::Element* ele,  ///< pointer to current element
        M& localmatrix,                                     ///< local matrix on element-level
        const Teuchos::RCP<Epetra_MultiVector>& global,     ///< global vector
        const int nsd                                       ///< number of space dimensions
    )
    {
      if (global == Teuchos::null) dserror("received a TEUCHOS::null pointer");
      if (nsd > global->NumVectors())
        dserror("Requested %d of %d available columns", nsd, global->NumVectors());
      const int iel = ele->NumNode();  // number of nodes
      if (((int)localmatrix.N()) != iel) dserror("local matrix has wrong number of columns");
      if (((int)localmatrix.M()) != nsd) dserror("local matrix has wrong number of rows");

      for (int i = 0; i < nsd; i++)
      {
        // access actual component column of multi-vector
        double* globalcolumn = (*global)[i];
        // loop over the element nodes
        for (int j = 0; j < iel; j++)
        {
          const int nodegid = (ele->Nodes()[j])->Id();
          const int lid = global->Map().LID(nodegid);
          if (lid < 0)
            dserror(
                "Proc %d: Cannot find gid=%d in Epetra_Vector", (*global).Comm().MyPID(), nodegid);
          localmatrix(i, j) = globalcolumn[lid];
        }
      }
      return;
    }

    /*!
    \brief extract local values from global node-based (multi) vector

    This function returns a column vector!

    \author henke
   */
    template <class M>
    void ExtractMyNodeBasedValues(
        const DRT::Element* ele, M& local, const Epetra_MultiVector& global)
    {
      const int numnode = ele->NumNode();
      const int numcol = global.NumVectors();
      if (((int)local.N()) != 1) dserror("local matrix must have one column");
      if (((int)local.M()) != numnode * numcol) dserror("local matrix has wrong number of rows");

      // loop over element nodes
      for (int i = 0; i < numnode; ++i)
      {
        const int nodegid = (ele->Nodes()[i])->Id();
        const int lid = global.Map().LID(nodegid);
        if (lid < 0)
          dserror("Proc %d: Cannot find gid=%d in Epetra_Vector", global.Comm().MyPID(), nodegid);

        // loop over multi vector columns (numcol=1 for Epetra_Vector)
        for (int col = 0; col < numcol; col++)
        {
          double* globalcolumn = (global)[col];
          local((col + (numcol * i)), 0) = globalcolumn[lid];
        }
      }
      return;
    }


    /*!
    \brief extract location vector of the element dofs based on the numdof of the discretization

    Instead of the element function "location vector" here the discretization is
    used to identify the number of degrees of freedom of each node

    \author henke
   */
    void DisBasedLocationVector(const DRT::Discretization& dis, const DRT::Element& ele,
        std::vector<int>& lm, const int num);

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
        const int& solvernumber,  ///< solver number for solving the resulting global system
        const enum INPAR::SCATRA::L2ProjectionSystemType& l2_proj_type =
            INPAR::SCATRA::l2_proj_system_std);

    Teuchos::RCP<Epetra_MultiVector> ComputeNodalL2Projection(Discretization& dis,
        const Epetra_Map& noderowmap, const Epetra_Map& elecolmap, const std::string& statename,
        const int& numvec, Teuchos::ParameterList& params, const int& solvernumber,
        const enum INPAR::SCATRA::L2ProjectionSystemType& l2_proj_type,
        const Epetra_Map* fullnoderowmap = NULL,
        const std::map<int, int>* slavetomastercolnodesmap = NULL,
        Epetra_Vector* const sys_mat_diagonal_ptr = NULL);

    Teuchos::RCP<Epetra_MultiVector> ComputeNodalL2Projection(Discretization& dis,
        const Epetra_Map& noderowmap, DRT::Element* const* coleles, const unsigned& numcolele,
        const std::string& statename, const int& numvec, Teuchos::ParameterList& params,
        const int& solvernumber, const enum INPAR::SCATRA::L2ProjectionSystemType& l2_proj_type,
        const Epetra_Map* fullnoderowmap = NULL,
        const std::map<int, int>* slavetomastercolnodesmap = NULL,
        Epetra_Vector* const sys_mat_diagonal_ptr = NULL);

    Teuchos::RCP<Epetra_MultiVector> SolveDiagonalNodalL2Projection(
        LINALG::SparseMatrix& massmatrix, Epetra_MultiVector& rhs, const int& numvec,
        const Epetra_Map& noderowmap, const Epetra_Map* fullnoderowmap = NULL,
        const std::map<int, int>* slavetomastercolnodesmap = NULL,
        Epetra_Vector* const sys_mat_diagonal_ptr = NULL);

    Teuchos::RCP<Epetra_MultiVector> SolveNodalL2Projection(LINALG::SparseMatrix& massmatrix,
        Epetra_MultiVector& rhs, const Epetra_Comm& comm, const int& numvec,
        const int& solvernumber, const Epetra_Map& noderowmap,
        const Epetra_Map* fullnoderowmap = NULL,
        const std::map<int, int>* slavetomastercolnodesmap = NULL);

    Teuchos::RCP<Epetra_MultiVector> PostSolveNodalL2Projection(
        const Teuchos::RCP<Epetra_MultiVector>& nodevec, const Epetra_Map& noderowmap,
        const Epetra_Map* fullnoderowmap, const std::map<int, int>* slavetomastercolnodesmap);

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
    \brief handles random numbers through boost

    \author wichmann
    */
    class Random
    {
     public:
      Random();

      virtual ~Random();

      /// get a random number
      double Uni();

      /// get a vector of random numbers of size count
      void Uni(std::vector<double>& randvec, int count);

      /// get a random number
      double Normal();

      /// get a vector of random numbers of size count
      void Normal(std::vector<double>& randvec, int count);

      /// set the random seed
      void SetRandSeed(const unsigned int seed);

      /// set the range for the uniform rng
      void SetRandRange(const double lower, const double upper);

      /// set the mean and variance for the normal rng
      void SetMeanVariance(const double mean, const double var);

     protected:
/// @name Random number generation
#if (BOOST_MAJOR_VERSION == 1) && (BOOST_MINOR_VERSION >= 47)
      /// random number engine
      boost::random::mt19937 rand_engine_;

      /// uniform random number distribution
      boost::random::uniform_real_distribution<double> uni_dist_;

      /// uniform random number generator
      boost::random::variate_generator<boost::random::mt19937&,
          boost::random::uniform_real_distribution<double>>
          uni_rand_no_;

      /// normal random number distribution
      boost::random::normal_distribution<double> norm_dist_;

      /// uniform random number generator
      boost::random::variate_generator<boost::random::mt19937&,
          boost::random::normal_distribution<double>>
          norm_rand_no_;
#else
      /********************************************************
       * The following declarations are deprecated in current *
       * Boosts, but they are all we got on FC8 and F12       *
       ********************************************************/
      /// random number engine
      boost::mt19937 rand_engine_;

      /// uniform random number distribution
      boost::uniform_real<double> uni_dist_;

      /// uniform random number generator
      boost::variate_generator<boost::mt19937&, boost::uniform_real<double>> uni_rand_no_;

      /// normal random number distribution
      boost::normal_distribution<double> norm_dist_;

      /// uniform random number generator
      boost::variate_generator<boost::mt19937&, boost::normal_distribution<double>> norm_rand_no_;
#endif
      //@}
    };

    /*!
    \brief handles restart after a certain walltime interval, step interval or on a user signal

    \author hammerl
    */
    class RestartManager
    {
     public:
      RestartManager();

      virtual ~RestartManager();

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
     *            \brief Default error handling of fgets().
     *
     *            \note In debug mode, this function prints out a warning if the
     *            End-of-File was reached without reading any character.
     *            \param output (in): output provided by the call of fgets()
     *            \param stream (in): the FILE pointer passed to fgets()
     *            \param filename (in): name of the file being processed by fgets()
     *            \throws dserror() occurs if the error indicator was set while reading
     *            \author Hagmeyer
     */
    void Checkfgets(char* output, FILE* stream, std::string filenename);
    /**
     *            \brief Default error handling of scanf().
     *
     *            \param output (in): output provided by the call of scanf()
     *            \throws dserror() occurs if the function returns without reading any
     * input \author Hagmeyer
     */
    void Checkscanf(int output);
  }  // namespace UTILS
}  // namespace DRT


#endif  // #ifndef UTILS_H
