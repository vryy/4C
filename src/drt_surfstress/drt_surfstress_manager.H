/*!----------------------------------------------------------------------
\file drt_surfstress_manager.H

<pre>
Maintainer: Lena Wiechert
            wiechert@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15303
</pre>

*----------------------------------------------------------------------*/

#ifndef SURFSTRESS_MANAGER_H
#define SURFSTRESS_MANAGER_H

#include <Teuchos_RefCountPtr.hpp>
#include <Epetra_Vector.h>

#include <Teuchos_ParameterList.hpp>
#include "../drt_lib/linalg_sparsematrix.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_io/io.H"


using namespace std;
using namespace Teuchos;


namespace UTILS
{
/*!
\brief Class controlling surface stresses due to interfacial phenomena
and containing all necessary history data
*/
class SurfStressManager
{
public:

   /*!
    \brief Standard Constructor

    \param discret (in): Discretization
  */
  SurfStressManager(DRT::Discretization& discret);

  /*!
    \brief Destructor
  */
  ~SurfStressManager() {};

  /*!
    \brief Call discretization to calculate additional contributions
    due to interfacial energy

    \param p (in): parameter list
    \param dism (in): displacements (n+1-alphaf)
    \param disn (in): displacements (n+1)
    \param fint (out): internal forces
    \param stiff (out): stiffness matrix
  */
  void EvaluateSurfStress(ParameterList& p,
                          const RefCountPtr<Epetra_Vector> dism,
                          const RefCountPtr<Epetra_Vector> disn,
                          RefCountPtr<Epetra_Vector> fint,
                          RefCountPtr<LINALG::SparseMatrix> stiff);

  /*!
    \brief update history variables
  */
  void Update();

  /*!
    \brief get history variables needed for writing restart. The column map based
    vectors used for calculations are exported to row map based ones
    needed for writing.
  */
  void GetHistory(RCP<Epetra_Vector> A_last_row,
                  RCP<Epetra_Vector> con_last_row,
                  RCP<Epetra_Vector> gamma_last_row)
  {
    LINALG::Export(*A_last_, *A_last_row);
    LINALG::Export(*con_last_, *con_last_row);
    LINALG::Export(*gamma_last_, *gamma_last_row);
  }

  /*!
    \brief get row map of surface elements
  */
  RCP<Epetra_Map> GetSurfRowmap() { return surfrowmap_; }

  /*!
    \brief set history variables after reading restart. The row map
    based vectors written in case of restart are exported to column
    based ones needed for calculations again.
  */
  void SetHistory(RCP<Epetra_Vector> A_last,
                  RCP<Epetra_Vector> con_last,
                  RCP<Epetra_Vector> gamma_last)
  {
    LINALG::Export(*A_last, *A_last_);
    LINALG::Export(*con_last, *con_last_);
    LINALG::Export(*gamma_last, *gamma_last_);
  }

  /*!
    \brief Routine to calculate additional internal forces and
    stiffnesses due to interfacial phenomena

    \param curvenum (in): number of time curve
    \param A (in): current interfacial area at (n+1-alphaf)
    \param Adiff (in): first partial derivatives of surface area at (n+1-alphaf)
    \param Adiff2 (in): second partial derivatives of surface area at (n+1-alphaf)
    \param A (in): current interfacial area at (n+1)
    \param Adiff (in): first partial derivatives of surface area at (n+1)
    \param fint (out): element internal force vector
    \param stiff (out): element stiffness matrix
    \param ID (in): ID of surface condition
    \param time (in): current time
    \param dt (in): timestep size
    \param surface_flag (in): flag for surfactant model (0) or ideal liquid (1)
    \param const_gamma (in): constant surface tension in case surface_flag=1
    \param k1xC (in): adsorption coefficient times bulk concentration
    \param k2 (in): desorption coefficient
    \param m1 (in): 1st isothermal slope
    \param m2 (in): 2nd isothermal slope
    \param gamma_0 (in): surface tension of water
    \param gamma_min (in): minimum surface stress
    \param gamma_min_eq (in): minimum equilibrium surface stress
    \param con_quot_max (in): maximum surfactant concentration
    \param con_quot_eq (in): equilibrium surfactant concentration
    \param alphaf (in): generalized-alpha parameter alphaf
    \param newstep (in): flag for new time step
  */
  void StiffnessAndInternalForces(const int curvenum,
                                  const double& A,
                                  const RCP<Epetra_SerialDenseVector> Adiff,
                                  const RCP<Epetra_SerialDenseMatrix> Adiff2,
                                  const double& Anew,
                                  const RCP<Epetra_SerialDenseVector> Adiffnew,
                                  Epetra_SerialDenseVector& fint,
                                  Epetra_SerialDenseMatrix& stiff,
                                  const int ID,
                                  const double time,
                                  const double dt,
                                  const int surface_flag,
                                  const double const_gamma,
                                  const double k1xC,
                                  const double k2,
                                  const double m1,
                                  const double m2,
                                  const double gamma_0,
                                  const double gamma_min,
                                  const double gamma_min_eq,
                                  const double con_quot_max,
                                  const double con_quot_eq,
                                  const double alphaf,
                                  const bool newstep);

  /*!
    \brief Surfactant model equations

    \param ID (in): ID of surface condition
    \param dgamma (out): current derivative of surface stress
    \param dt (in): timestep size
    \param k1xC (in): adsorption coefficient times bulk concentration
    \param k2 (in): desorption coefficient
    \param m1 (in): 1st isothermal slope
    \param m2 (in): 2nd isothermal slope
    \param gamma_0 (in): surface tension of water
    \param gamma_min (in): minimum surface stress
    \param gamma_min_eq (in): minimum equilibrium surface stress
    \param con_quot_max (in) maximum surfactant concentration
    \param alphaf (in): generalized-alpha parameter alphaf
    \param update (in): flag for new time step
  */
  void SurfactantModel(const int ID,
                       double& dgamma,
                       const double dt,
                       double k1xC,
                       double k2,
                       const double m1,
                       const double m2,
                       const double gamma_0,
                       const double gamma_min,
                       const double gamma_min_eq,
                       const double con_quot_max,
                       const double alphaf,
                       const bool newstep);

private:

  //don't want = operator and cctor
  SurfStressManager operator = (const SurfStressManager& old);
  SurfStressManager(const UTILS::SurfStressManager& old);


  //! name@ Private Attributes

  // current discretization
  DRT::Discretization& discret_;

  // row map of surfaces
  RCP<Epetra_Map> surfrowmap_;

  // current surface areas (n+1)^i
  RefCountPtr<Epetra_Vector> A_current_;

  // surface areas of former time step (n)
  RefCountPtr<Epetra_Vector> A_last_;

  // current non-dimensionalized surfactant concentrations (n+1)^i
  RefCountPtr<Epetra_Vector> con_current_;

  // non-dimensionalized surfactant concentrations of former time step (n)
  RefCountPtr<Epetra_Vector> con_last_;

  // current surface stresses (n+1)^i
  RefCountPtr<Epetra_Vector> gamma_current_;

  // surface stresses of former time step (n)
  RefCountPtr<Epetra_Vector> gamma_last_;

  //@}
}; //class
}//namespace DRT

#endif /*SURFSTRESS_MANAGER_H*/
