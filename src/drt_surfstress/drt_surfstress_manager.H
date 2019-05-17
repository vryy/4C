/*---------------------------------------------------------------------*/
/*!

\brief surface stresses due to interfacial phenomena

\maintainer Martin Kronbichler

\level 2

*/
/*---------------------------------------------------------------------*/

#ifndef SURFSTRESS_MANAGER_H
#define SURFSTRESS_MANAGER_H

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Epetra_Vector.h>
#include <Epetra_SerialDenseMatrix.h>
#include <Epetra_SerialDenseVector.h>


// forward declarations

namespace DRT
{
  class Discretization;
}

namespace IO
{
  class DiscretizationWriter;
}

namespace LINALG
{
  class SparseOperator;
}

namespace UTILS
{
  /*!
  \brief Class controlling surface stresses due to interfacial phenomena
  and containing all necessary history data

  \author wiechert (yoshihara@lnm.mw.tum.de)
  */
  class SurfStressManager
  {
   public:
    /*!
      \brief Standard Constructor

      \param discret (in): Discretization
      \param sdynparams (in): structural dynamics parameters
      \param file_prefix (in): prefix for output file name

      Note that the prefix needs to be handed in since the
      SurfStressManager can not distinguish whether it is allocated on
      the macro- or a possible micro-scale (and if so, in which
      macro-scale element and at which macro-scale GP).
    */
    SurfStressManager(Teuchos::RCP<DRT::Discretization> discret, Teuchos::ParameterList sdynparams,
        const std::string file_prefix);


    /*!
      \brief Destructor
    */
    ~SurfStressManager(){};


    /*!
      \brief Return if there are surface stress conditions
    */
    bool HaveSurfStress() const { return havesurfstress_; };

    /*!
      \brief Write binary results (current gamma, Gamma, A)
    */
    void WriteResults(const int istep, const double timen);

    /*!
      \brief Write restart (current gamma, Gamma, A)
    */
    void WriteRestart(const int istep, const double timen);

    /*!
      \brief Read restart (current gamma, Gamma, A)
    */
    void ReadRestart(const int step, const std::string file_prefix, const bool serial = false);

    /*!
      \brief Call discretization to calculate additional contributions
      due to interfacial energy

      \param p (in): parameter list
      \param disn (in): displacements (n+1)
      \param fint (out): internal forces
      \param stiff (out): stiffness matrix
    */
    void EvaluateSurfStress(Teuchos::ParameterList& p, const Teuchos::RCP<Epetra_Vector> disn,
        Teuchos::RCP<Epetra_Vector> fint, Teuchos::RCP<LINALG::SparseOperator> stiff);

    /*!
      \brief Update history variables
    */
    void Update();

    /*!
      \brief Calculate additional internal forces and
      stiffnesses due to interfacial phenomena

      \param curvenum (in): number of time curve
      \param A (in): current interfacial area at (n+1)
      \param Adiff (in): first partial derivatives of surface area at (n+1)
      \param Adiff2 (in): second partial derivatives of surface area at (n+1)
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
      \param newstep (in): flag for new time step
    */
    void StiffnessAndInternalForces(const int curvenum, const double& A,
        const Teuchos::RCP<Epetra_SerialDenseVector> Adiff,
        const Teuchos::RCP<Epetra_SerialDenseMatrix> Adiff2, Epetra_SerialDenseVector& fint,
        Epetra_SerialDenseMatrix& stiff, const int ID, const double time, const double dt,
        const int surface_flag, const double const_gamma, const double k1xC, const double k2,
        const double m1, const double m2, const double gamma_0, const double gamma_min,
        const double gamma_min_eq, const double con_quot_max, const double con_quot_eq,
        const bool newstep);

    /*!
      \brief Surfactant model equations

      \param ID (in): ID of surface condition
      \param gamma (out): current surface stresses
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
      \param update (in): flag for new time step
    */
    void SurfactantModel(const int ID, double& gamma, double& dgamma, const double dt, double k1xC,
        double k2, const double m1, const double m2, const double gamma_0, const double gamma_min,
        const double gamma_min_eq, const double con_quot_max, const bool newstep);

   private:
    // don't want = operator and cctor
    SurfStressManager operator=(const SurfStressManager& old);
    SurfStressManager(const UTILS::SurfStressManager& old);


    //! name@ Private Attributes

    bool havesurfstress_;
    double timen_;

    // current discretization
    Teuchos::RCP<DRT::Discretization> discret_;

    // current surface condition discretization
    Teuchos::RCP<DRT::Discretization> surfdiscret_;

    // output
    Teuchos::RCP<IO::DiscretizationWriter> surfoutput_;

    // row map of surfaces
    Teuchos::RCP<Epetra_Map> surfrowmap_;

    // current surface areas (n+1)^i
    Teuchos::RCP<Epetra_Vector> A_current_;

    // surface areas of former time step (n)
    Teuchos::RCP<Epetra_Vector> A_last_;

    // current non-dimensionalized surfactant concentrations (n+1)^i
    Teuchos::RCP<Epetra_Vector> con_current_;

    // non-dimensionalized surfactant concentrations of former time step (n)
    Teuchos::RCP<Epetra_Vector> con_last_;

    // current surface stresses (n+1)^i
    Teuchos::RCP<Epetra_Vector> gamma_current_;

    //@}
  };  // class
}  // namespace UTILS

#endif /*SURFSTRESS_MANAGER_H*/
