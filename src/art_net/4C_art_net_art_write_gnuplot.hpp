/*----------------------------------------------------------------------*/
/*! \file
\brief Method to print the arteries in a way that could be displayed by
\gnuplot

\level 3


*----------------------------------------------------------------------*/

#ifndef FOUR_C_ART_NET_ART_WRITE_GNUPLOT_HPP
#define FOUR_C_ART_NET_ART_WRITE_GNUPLOT_HPP

#include "4C_config.hpp"

#include "4C_fem_discretization.hpp"
#include "4C_io.hpp"
#include "4C_linalg_utils_sparse_algebra_math.hpp"

#include <Epetra_MpiComm.h>
#include <Teuchos_RCP.hpp>

#include <iostream>

FOUR_C_NAMESPACE_OPEN

namespace Arteries
{
  namespace UTILS
  {
    //--------------------------------------------------------------------
    // Wrapper class (to be called from outside) for outputing
    //--------------------------------------------------------------------

    /*!
    \brief 1d-artery gnuplot output condition wrapper
    this class is meant to do some organisation stuff
    */
    class ArtWriteGnuplotWrapper
    {
      friend class ArtNetExplicitTimeInt;


     public:
      /*!
      \brief Standard Constructor
      */
      ArtWriteGnuplotWrapper(
          Teuchos::RCP<Core::FE::Discretization> actdis, Teuchos::ParameterList& params);

      /*!
      \brief Destructor
      */
      virtual ~ArtWriteGnuplotWrapper() = default;

      /*!
      \brief Standard write
      */
      void write(Teuchos::ParameterList& params);


     private:
      /*!
      \brief all single artery write conditions
     */
      std::map<const int, Teuchos::RCP<class ArtWriteGnuplot>> agmap_;
      std::map<const int, const std::vector<int>*> agnode_map_;


      //! 1d artery discretization
      Teuchos::RCP<Core::FE::Discretization> discret_;

    };  // class ArtWriteGnuplotWrapper



    //--------------------------------------------------------------------
    // Actual artery gnuplot output condition
    //--------------------------------------------------------------------
    /*!
    \brief 1d-artery gnuplot output condition

    */

    class ArtWriteGnuplot
    {
      friend class ArtWriteGnuplotWrapper;

     public:
      /*!
      \brief Standard Constructor
     */
      ArtWriteGnuplot(int ArteryNum);

      /*!
      \brief Empty Constructor
      */
      ArtWriteGnuplot();

      /*!
      \brief Destructor
      */
      virtual ~ArtWriteGnuplot() = default;

     protected:
      /*!
      \Solve the write the results of an artery
      */
      void write(Teuchos::RCP<Core::FE::Discretization> discret, Teuchos::ParameterList& params,
          const std::vector<int>* nodes);


     private:
      //! An epetra wrapper for Dense matrix solver
      Teuchos::RCP<std::ofstream> fout_;

      //! the Artery number
      int artery_num_;


    };  // class ArtWriteGnuplot
  }     // namespace UTILS
}  // namespace Arteries


FOUR_C_NAMESPACE_CLOSE

#endif
