/*!-----------------------------------------------------------------------------------------------*
\file scale_sep_gmo.H

  \brief provides scale separation operators for les of turbulent flows

<pre>
Maintainer: Ursula Rasthofer
            rasthofer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>
 *------------------------------------------------------------------------------------------------*/

/*
 * This class provides a general framework for dealing with scale separation for turbulence modeling.
 * So far, it contains a scale separation operator of the form of a geometric multigrid operator, only. For more details
 * on this filter, read the comment at the beginning of the function ConstructSepMatGeoMultigrid().
 * This class could be extended so that the complete filtering process is controlled by it. This includes
 * the inclusion of the already existing box-filter for the dynamic Smagorinsky model and the algebraic
 * multigrid operator for AVM3.
 */


#ifndef SCALE_SEP_GMO_H
#define SCALE_SEP_GMO_H

#ifdef PARALLEL
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif /* PARALLEL */

#include "../drt_inpar/inpar_fluid.H"
#include "../drt_inpar/inpar_turbulence.H"
#include "../drt_lib/drt_discret.H"
#include "../linalg/linalg_utils.H"

namespace FLD {

class LESScaleSeparation
{
public:
   /*!
  \brief Constructor

  */
  LESScaleSeparation(
  INPAR::FLUID::ScaleSeparation scale_sep, //! separation type
  Teuchos::RCP<DRT::Discretization> discret);

   /*!
  \brief Destructor

  */
  virtual ~LESScaleSeparation();

  void ApplyScaleSeparation(Teuchos::RCP<Epetra_Vector> vel,
                            Teuchos::RCP<Epetra_Vector> fsvel);


protected:

  /*!
  \brief sort criterium for double values up to a tolerance of 10-9

  This is used to create sets of doubles (e.g. coordinates)

  */
  class LineSortCriterion
  {
  public:
    bool operator() (
      const double& p1,
      const double& p2) const {
      return (p1 < p2 - 1E-9);
    }

  protected:
  private:

  };


private:
  //!
  void ConstructSepMatGeoMultigrid();

  //! functions which help
  DRT::Node* FindNode(double x,
                      double y,
                      double z);
  //! store values and ids in vector
  void StoreInVector(const std::vector<DRT::Node*>& nodes,
                     const std::vector<double>& funct,
                     std::vector<double>& weights,
                     std::vector<int>& colids);

  //! filter type
  enum INPAR::FLUID::ScaleSeparation scale_sep_;
  //! pointer to discretization
  Teuchos::RCP<DRT::Discretization> discret_;
  //! bool: true if Sep_-matrix has been build
  bool sepmat_build_;
  //! separation matrix
  Teuchos::RCP<LINALG::SparseMatrix> Sep_;
};

}

#endif

