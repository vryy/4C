// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_ART_NET_ART_WRITE_GNUPLOT_HPP
#define FOUR_C_ART_NET_ART_WRITE_GNUPLOT_HPP

#include "4C_config.hpp"

#include "4C_fem_discretization.hpp"
#include "4C_io.hpp"
#include "4C_linalg_utils_sparse_algebra_math.hpp"

#include <Epetra_MpiComm.h>

#include <iostream>
#include <memory>

FOUR_C_NAMESPACE_OPEN

namespace Arteries
{
  namespace Utils
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
          std::shared_ptr<Core::FE::Discretization> actdis, Teuchos::ParameterList& params);

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
      std::map<const int, std::shared_ptr<class ArtWriteGnuplot>> agmap_;
      std::map<const int, const std::vector<int>*> agnode_map_;


      //! 1d artery discretization
      std::shared_ptr<Core::FE::Discretization> discret_;

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
      void write(Core::FE::Discretization& discret, Teuchos::ParameterList& params,
          const std::vector<int>* nodes);


     private:
      //! An epetra wrapper for Dense matrix solver
      std::shared_ptr<std::ofstream> fout_;

      //! the Artery number
      int artery_num_;


    };  // class ArtWriteGnuplot
  }  // namespace Utils
}  // namespace Arteries


FOUR_C_NAMESPACE_CLOSE

#endif
