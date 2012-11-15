/*!----------------------------------------------------------------------
\file drt_potential_manager.H

\brief  Class controlling surface stresses due to potential forces
        between mesoscopic structures

<pre>
Maintainer: Ursula Mayer
            mayer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15257
</pre>

*----------------------------------------------------------------------*/

#ifndef DRT_POTENTIAL_MANAGER_H_
#define DRT_POTENTIAL_MANAGER_H_

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>

#include "../drt_potential/drt_potential_surface.H"
#include "../drt_potential/drt_potential_volume.H"
#include "../drt_geometry/searchtree.H"

using namespace std;


namespace POTENTIAL
{
/*!
\brief Class controlling surface and volume stresses due to potential forces
as well as their approximations
*/
class PotentialManager
{
public:

  /*!
    \brief Standard Constructor

    \param discret (in): Discretization
   */
  PotentialManager(
      const Teuchos::RCP<DRT::Discretization>     discretRCP,
      DRT::Discretization&                        discret
      );

  /*!
    \brief Destructor
   */
  ~PotentialManager() {};

  
  /*!
    \brief Read parameter
   */
  void ReadParameter();
  
  
  /*!
    \brief Call discretization to calculate additional contributions
    due to interfacial energy

    \param p (in): parameter list
    \param fint (out): internal forces
    \param stiff (out): stiffness matrix
   */
  void EvaluatePotential(
      ParameterList&                      p,
      RCP<Epetra_Vector>          disp,
      RCP<Epetra_Vector>          fint,
      RCP<LINALG::SparseMatrix>   stiff);

  
  /*!
    \brief TEST Call discretization to calculate additional contributions
    due to interfacial energy

    \param p (in): parameter list
    \param fint (out): internal forces
    \param stiff (out): stiffness matrix
   */
  void TestEvaluatePotential( 
      ParameterList&                    p,
      RCP<Epetra_Vector>        disp,
      RCP<Epetra_Vector>        fint,
      RCP<LINALG::SparseMatrix> stiff,
      const double                      time,
      const int                         step
      );
  
  
  void StiffnessAndInternalForcesPotential(
      const DRT::Element*             element,
      const DRT::UTILS::GaussRule3D&  gaussrule,
      ParameterList&                  eleparams,
      vector<int>&                    lm,
      Epetra_SerialDenseMatrix&       K_stiff,
      Epetra_SerialDenseVector&       F_int);
  
  

  /*!
    \brief Routine to calculate additional internal forces and
    stiffnesses due to interfacial phenomena for surface elements

    \param element    (in)  : element
    \param gaussrule  (in)  : gaussrule of element
    \param params     (in)  : first partial derivatives of surface area
    \param K_surf     (out) : element stiffness matrix
    \param fint       (out) : element internal force vector
   */
  void StiffnessAndInternalForcesPotential(
      const DRT::Element*               element,
      const DRT::UTILS::GaussRule2D&    gaussrule,
      ParameterList&                    eleparams,
      vector<int>&                      lm,
      Epetra_SerialDenseMatrix&         K_stiff,
      Epetra_SerialDenseVector&         fint
      );


  /*!
    \brief Routine to calculate additional internal forces and
    stiffnesses due to interfacial phenomena for line elements

    \param element    (in)  : element
    \param gaussrule  (in)  : gaussrule of element
    \param params     (in)  : first partial derivatives of surface area
    \param K_surf     (out) : element stiffness matrix
    \param fint       (out) : element internal force vector
   */
  void StiffnessAndInternalForcesPotential(
      const DRT::Element*               element,
      const DRT::UTILS::GaussRule1D&    gaussrule,
      ParameterList&                    eleparams,
      vector<int>&                      lm,
      Epetra_SerialDenseMatrix&         K_stiff,
      Epetra_SerialDenseVector&         fint
      );
  
  bool ComputeAnalyticalSolution();

private:

  //don't want = operator and cctor
  PotentialManager operator = (const PotentialManager& old);
  PotentialManager(const POTENTIAL::PotentialManager& old);

  //! name@ Private Attributes

  //!  current discretization
  RCP<DRT::Discretization>                      discretRCP_;

  //!  current discretization
  DRT::Discretization&                          discret_;

  //!  parameter list
  Teuchos::ParameterList                        params_;

  //! surface potential
  Teuchos::RCP<POTENTIAL::SurfacePotential>     surfacePotential_;

  //! volume potential
  Teuchos::RCP<POTENTIAL::VolumePotential>      volumePotential_;

  //! surface
  bool                                          surface_;

  //! volume
  bool                                          volume_;
  
  
  GEO::TreeType                                 treetype_;

}; //class

}//namespace UTILS

#endif /*DRT_POTENTIAL_MANAGER_H_*/
