/*!----------------------------------------------------------------------
\file drt_potential_manager.H

<pre>
Maintainer: Ursula Mayer
            mayer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15257
</pre>

*----------------------------------------------------------------------*/

#ifndef DRT_POTENTIAL_MANAGER_H_
#define DRT_POTENTIAL_MANAGER_H_

#include <Teuchos_RefCountPtr.hpp>
#include <Teuchos_ParameterList.hpp>

#include <Epetra_Vector.h>
#include <Epetra_Map.h>
#include <Epetra_Operator.h>

#include "../drt_lib/linalg_systemmatrix.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_io/io.H"
#include "../drt_xfem/xfsi_searchtree.H"


using namespace std;
using namespace Teuchos;


namespace DRT
{
/*!
\brief Class controlling surface stresses due to potential forces
and containing all necessary history data
*/
class PotentialManager
{
public:

  /*!
    \brief Standard Constructor

    \param discret (in): Discretization
   */
  PotentialManager( Teuchos::RCP<DRT::Discretization> discretRCP,
      DRT::Discretization& discret);

  /*!
    \brief Destructor
   */
  ~PotentialManager() {};

  /*!
    \brief Call discretization to calculate additional contributions
    due to interfacial energy

    \param p (in): parameter list
    \param fint (out): internal forces
    \param stiff (out): stiffness matrix
   */
  void EvaluatePotential( ParameterList& p,
      RefCountPtr<Epetra_Vector> disp,
      RefCountPtr<Epetra_Vector> fint,
      RefCountPtr<LINALG::SparseMatrix> stiff);

  /*!
    \brief update history variables
   */
  void Update();

  /*!
    \brief get history variables needed for writing restart. Note that
    the temporal vectors needs to be written since in the old ones we
    still have the data of the former step. The column map based
    vectors used for calculations are exported to row map based ones
    needed for writing.
   */
  void GetHistory(RCP<Epetra_Vector> A_old_temp_row);

  /*!
    \brief get row map of surface elements
   */
  RCP<Epetra_Map> GetSurfRowmap() { return surfrowmap_; }

  /*!
    \brief set history variables after reading restart. The row map
    based vectors written in case of restart are exported to column
    based ones needed for calculations again.
   */
  void SetHistory(RCP<Epetra_Vector> A_old)
  {
    LINALG::Export(*A_old, *A_old_);
  }

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
      Epetra_SerialDenseVector& fint,
      Epetra_SerialDenseMatrix& K_surf,
      const int ID,
      const double time,
      const double dt,
      const double label,
      const double depth,
      const double rootDist,
      const double cutOff);


private:

  //don't want = operator and cctor
  PotentialManager operator = (const PotentialManager& old);
  PotentialManager(const DRT::PotentialManager& old);

  /*!
  \brief Update displacements
   */
  void UpdateDisplacementsOfBoundaryDiscretization(
    Teuchos::RCP<Epetra_Vector> idisp_solid
    );

  /*!
  \brief  set history variables after reading restart. The row map
          based vectors written in case of restart are exported to column
          based ones needed for calculations again.
   */
  void searchElementsInCutOffRadius(    
    const vector<double>&   point,
    set< int >&             potentialElementIds, 
    const double            cutOff
    );
  
  
  //! name@ Private Attributes

  // current discretization
  RCP<DRT::Discretization> discretRCP_;

  // current discretization
  DRT::Discretization& discret_;

  // row map of surfaces
  RCP<Epetra_Map> surfrowmap_;

  // surface areas of former iteration
  RefCountPtr<Epetra_Vector> A_old_temp_;

  // surface areas of former time step
  RefCountPtr<Epetra_Vector> A_old_;

  //! XSearchTree xTree_
  Teuchos::RCP<XFEM::XSearchTree>  xTree_;

  //! current position of boundary elements
  std::map<int,vector<double> > currentboundarypositions_;

  LINALG::MapExtractor                   potboundary_;

  //! @name local copies of input parameters
  Teuchos::RCP<DRT::Discretization>      boundarydis_;
  Teuchos::RCP<Teuchos::ParameterList>   params_;
  //@}

  //! @name state vectors using the solid boundary parallel distribution (row vectors)
  Teuchos::RCP<Epetra_Vector> idisp_total_; 
  Teuchos::RCP<Epetra_Vector> idisp_onproc_;


}; //class
}//namespace DRT

#endif /*DRT_POTENTIAL_MANAGER_H_*/
