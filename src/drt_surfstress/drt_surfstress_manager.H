/*!----------------------------------------------------------------------
\file drt_surfstress_manager.H

\class DRT::SurfStressManager

\brief Class controlling surface stresses due to interfacial phenomena
and containing all necessary history data

<pre>
Maintainer: Lena Wiechert
            wiechert@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15303
</pre>

*----------------------------------------------------------------------*/

#ifndef SURFSTRESS_MANAGER_H
#define SURFSTRESS_MANAGER_H

#include "Teuchos_RefCountPtr.hpp"
#include "Epetra_Vector.h"

#include "Teuchos_ParameterList.hpp"
#include "../drt_lib/drt_discret.H"


using namespace std;
using namespace Teuchos;


namespace DRT
{
class SurfStressManager
{
public:

   /*!
    \brief Standard Constructor

    \param discret (in): Discretization
  */
  SurfStressManager(DRT::Discretization& discret, const int numsurf,
                    const Epetra_Map& surfmap);

  /*!
    \brief Destructor
  */
  ~SurfStressManager() {};

  /*!
    \brief Call discretization to calculate additional contributions
    due to interfacial energy

    \param p (in): parameter list
    \param fint (out): internal forces
    \param stiff (out): stiffness matrix
  */
  void EvaluateSurfStress(ParameterList& p,
                          RefCountPtr<Epetra_Vector> disp,
                          RefCountPtr<Epetra_Vector> fint,
                          RefCountPtr<Epetra_CrsMatrix> stiff);

  /*!
    \brief Routine to calculate additional internal forces and
    stiffnesses due to interfacial phenomena

    \param curvenum (in): number of time curve
    \param xs (in): current nodal positions
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
    \param con_quot_max (in) maximum surfactant concentration
    \param con_quot_eq (in) equilibrium surfactant concentration
  */
  void StiffnessAndInternalForces(const int curvenum,
                                  const Epetra_SerialDenseMatrix& xs,
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
                                  const double con_quot_eq);

  /*!
    \brief Routine to calculate current interfacial area and its
    derivatives with respect to the nodal displacements

    \param xs (in): current nodal positions
    \param A (out): current interfacial area
    \param Adiff (out): first partial derivative of interfacial area
    \param Adiff2 (out): second partial derivative of interfacial area
  */
  void soh8_surface_calc(const Epetra_SerialDenseMatrix& xs,
                         double& A,
                         Epetra_SerialDenseVector& Adiff,
                         Epetra_SerialDenseMatrix& Adiff2);

  /*!
    \brief Routine to calculate current surface stress and its
    derivative with respect to the change in interfacial area

    \param gamma (out): current surface stress
    \param dgamma (out): current derivative of surface stress
    \param con_quot (in/out): old/current interfacial surfactant concentration
    \param A_old (in): surface area of last time step
    \param A_new (in): current surface area
    \param dt (in): timestep size
    \param k1xC (in): adsorption coefficient times bulk concentration
    \param k2 (in): desorption coefficient
    \param m1 (in): 1st isothermal slope
    \param m2 (in): 2nd isothermal slope
    \param gamma_0 (in): surface tension of water
    \param gamma_min (in): minimum surface stress
    \param gamma_min_eq (in): minimum equilibrium surface stress
    \param con_quot_max (in) maximum surfactant concentration
  */
  void gamma_calc(double& gamma,
                  double& dgamma,
                  double& con_quot,
                  const double A_old,
                  const double A_new,
                  const double dt,
                  const double k1xC,
                  const double k2,
                  const double m1,
                  const double m2,
                  const double gamma_0,
                  const double gamma_min,
                  const double gamma_min_eq,
                  const double con_quot_max);

private:

  //don't want = operator and cctor
  SurfStressManager operator = (const SurfStressManager& old);
  SurfStressManager(const DRT::SurfStressManager& old);


  //! name@ Private Attributes

  // current discretization
  DRT::Discretization& discret_;

  // surface areas of former iteration
  RefCountPtr<Epetra_Vector> A_old_temp_;

  // surface areas of former time step
  RefCountPtr<Epetra_Vector> A_old_;

  // non-dimensionalized surfactant concentrations of former iteration
  RefCountPtr<Epetra_Vector> con_quot_temp_;

  // non-dimensionalized surfactant concentrations of former time step
  RefCountPtr<Epetra_Vector> con_quot_;

  // current time
  RefCountPtr<Epetra_Vector> time_;

  //@}
}; //class
}//namespace DRT

#endif /*SURFSTRESS_MANAGER_H*/
