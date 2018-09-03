/*!----------------------------------------------------------------------
\file art_write_gnuplot.H
\brief Method to print the arteries in a way that could be displayed by
\gnuplot

\level 3

\maintainer Lena Yoshihara

*----------------------------------------------------------------------*/

#ifndef ART_WRITE_GNUPLOT_H_
#define ART_WRITE_GNUPLOT_H_

#ifdef PARALLEL
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif /* PARALLEL */

#include <iostream>
#include "Teuchos_RCP.hpp"
#include "../drt_lib/drt_discret.H"
#include "../linalg/linalg_utils.H"
#include "../drt_io/io.H"

namespace ART
{
  namespace UTILS
  {
    //--------------------------------------------------------------------
    // Wrapper class (to be called from outside) for outputing
    //--------------------------------------------------------------------

    /*!
    \brief 1d-artery gnuplot output condition wrapper
    this class is meant to do some organisation stuff
    \author ismail (ismail@lnm.mw.tum.de)
    */
    class ArtWriteGnuplotWrapper
    {
      friend class ArtNetExplicitTimeInt;


     public:
      /*!
      \brief Standard Constructor
      */
      ArtWriteGnuplotWrapper(
          Teuchos::RCP<DRT::Discretization> actdis, Teuchos::ParameterList& params);

      /*!
      \brief Destructor
      */
      virtual ~ArtWriteGnuplotWrapper(){};

      /*!
      \brief Standard write
      */
      void Write(Teuchos::ParameterList& params);


     private:
      /*!
      \brief all single artery write conditions
     */
      std::map<const int, Teuchos::RCP<class ArtWriteGnuplot>> agmap_;
      std::map<const int, const std::vector<int>*> agnode_map_;


      //! 1d artery discretization
      Teuchos::RCP<DRT::Discretization> discret_;

    };  // class ArtWriteGnuplotWrapper



    //--------------------------------------------------------------------
    // Actual artery gnuplot output condition
    //--------------------------------------------------------------------
    /*!
    \brief 1d-artery gnuplot output condition

    \author ismail (ismail@lnm.mw.tum.de)
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
      virtual ~ArtWriteGnuplot() {}

     protected:
      /*!
      \Solve the write the results of an artery
      */
      void Write(Teuchos::RCP<DRT::Discretization> discret, Teuchos::ParameterList& params,
          const std::vector<int>* nodes);


     private:
      //! An epetra wrapper for Dense matrix solver
      Teuchos::RCP<std::ofstream> fout_;

      //! the Artery number
      int ArteryNum_;


    };  // class ArtWriteGnuplot
  }     // namespace UTILS
}  // namespace ART


#endif  // ART_WRITE_GNUPLOT_H_
