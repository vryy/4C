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

#include <Teuchos_RefCountPtr.hpp>
#include <Epetra_Vector.h>

#include <Teuchos_ParameterList.hpp>
#include "../drt_lib/linalg_systemmatrix.H"
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
                          RefCountPtr<LINALG::SparseMatrix> stiff);

  /*!
    \brief update history variables and time if necessary
  */
  void Update(const double time, const int surfaceflag, const int ID);

  /*!
    \brief calculate interfacial area and its derivatives with respect
    to the displacements
    \param dxyzdrs (in): dXYZ / drs
    \param derivatives (in): derivatives of shape functions
    \params gpccord (in): coordinates of Gauss points
    \params gpweight (in): Gaussian weights
    \params ngp (in): number of Gauss points
    \params nnode (in): number of nodes
    \param A (out): current interfacial area
    \param Adiff (out): first partial derivatives of interfacial area
    \param Adiff2 (out): second partial derivatives of interfacial area
  */
  void SurfaceCalc(const vector<Epetra_SerialDenseMatrix>& dxyzdrs,
                   const vector<Epetra_SerialDenseMatrix>& derivatives,
                   const Epetra_SerialDenseMatrix& gpccord,
                   const Epetra_SerialDenseVector& gpweight,
                   const int ngp,
                   const int nnode,
                   double& A,
                   Epetra_SerialDenseVector& Adiff,
                   Epetra_SerialDenseMatrix& Adiff2);

  /*!
    \brief Routine to calculate additional internal forces and
    stiffnesses due to interfacial phenomena

    \param curvenum (in): number of time curve
    \param A (in): current interfacial area
    \param Adiff (in): first partial derivatives of surface area
    \param Adiff2 (in): second partial derivatives of surface area
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
                                  const double& A,
                                  const Epetra_SerialDenseVector& Adiff,
                                  const Epetra_SerialDenseMatrix& Adiff2,
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
