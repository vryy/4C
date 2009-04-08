/*!----------------------------------------------------------------------
\file drt_potential_manager.H

\brief Class controlling surface stresses due to potential forces
between interfaces

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
#include "Epetra_Import.h"

#include "../drt_lib/linalg_sparsematrix.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_lib/linalg_serialdensematrix.H"
#include "../drt_lib/linalg_serialdensevector.H"
#include "../drt_io/io.H"
#include "../drt_geometry/searchtree.H"
#include "../drt_geometry/intersection_service.H"
#include "../drt_fem_general/drt_utils_integration.H"
#include "../drt_surfstress/potential_dofset.H"



using namespace std;
using namespace Teuchos;


namespace UTILS
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
  PotentialManager(
      Teuchos::RCP<DRT::Discretization>     discretRCP,
      DRT::Discretization&            discret
      );

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
      ParameterList&                    params,
      vector<int>&                      lm,
      Epetra_SerialDenseMatrix&         K_surf,
      Epetra_SerialDenseVector&         fint);
  
  
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
      ParameterList&                    params,
      vector<int>&                      lm,
      Epetra_SerialDenseMatrix&         K_surf,
      Epetra_SerialDenseVector&         fint);


private:

  //don't want = operator and cctor
  PotentialManager operator = (const PotentialManager& old);
  PotentialManager(const UTILS::PotentialManager& old);


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
    const LINALG::Matrix<3,1>&        point,
    std::map<int,std::set<int> >&     potentialElementIds,
    const double                      radius
    );


  /*!
  \brief  set history variables after reading restart. The row map
          based vectors written in case of restart are exported to column
          based ones needed for calculations again.
   */
  void treeSearchElementsInCutOffRadius(
      const LINALG::Matrix<3,1>&      point,
      std::map<int,std::set<int> >&   potentialElementIds,
      const double                    radius,
      const int                       label
      );


  /*!
  \brief  evaluate potential condition based on Epetra_FECrsMatrix
   */
  void EvaluatePotentialCondition(
      ParameterList&                          params,
      RefCountPtr<LINALG::SparseMatrix>       systemmatrix1,
      RefCountPtr<LINALG::SparseMatrix>       systemmatrix2,
      RefCountPtr<Epetra_Vector>              systemvector1,
      Teuchos::RCP<Epetra_Vector>             systemvector2,
      Teuchos::RCP<Epetra_Vector>             systemvector3,
      const string&                           condstring
      );
    

  /*!
  \brief  compute internal force and stiffness matrix
          for surface elements
   */
  void computeFandK(
      const DRT::Element*              actEle,
      const DRT::UTILS::GaussRule2D&   gaussrule,
      std::map<int,std::set<int> >&    potElements,
      vector<int>&                     lm,
  	  Epetra_SerialDenseMatrix&        K_surf,
  	  Epetra_SerialDenseVector&        F_int,
  	  RefCountPtr<DRT::Condition> 		 cond,
  		const int                        label,
      const double                     curvefac
      );
  
  
  /*!
  \brief  compute internal force and stiffness matrix
          for line elements
   */
  void computeFandK(
      const DRT::Element*              actEle,
      const DRT::UTILS::GaussRule1D&   gaussrule,
      std::map<int,std::set<int> >&    potElements,
      vector<int>&                     lm,
      Epetra_SerialDenseMatrix&        K_surf,
      Epetra_SerialDenseVector&        F_int,
      RefCountPtr<DRT::Condition>      cond,
      const int                        label,
      const double                     curvefac
      );
  
  
  /*!
  \brief  select correct potential out of condition
   */
  void EvaluatePotentialfromCondition(
		RefCountPtr<DRT::Condition> 	cond,
		const LINALG::Matrix<3,1>&    x,
    const LINALG::Matrix<3,1>&    y,
    LINALG::Matrix<3,1>&          potderiv1,
    LINALG::Matrix<3,3>&          potderiv2
    );
  
  
  /*!
  \brief  compute first and second derivative of the Lennard-Jones potential
   */
  void EvaluateLennardJonesPotential(
      const double                depth,
      const double                rootDist,
      const LINALG::Matrix<3,1>&  x,
      const LINALG::Matrix<3,1>&  y,
      LINALG::Matrix<3,1>&        potderiv1,
      LINALG::Matrix<3,3>&        potderiv2
      );
  
      
  /*!
  \brief  compute first and second derivative of the Zeta potential
   */
  void EvaluateZetaPotential(
    const double                  zeta_param_1,
    const double                  zeta_param_2,
    const LINALG::Matrix<3,1>&    x,
    const LINALG::Matrix<3,1>&    y,
    LINALG::Matrix<3,1>&          potderiv1,
    LINALG::Matrix<3,3>&          potderiv2
    );    
  
  
  /*!
  \brief  computes distance related structures
   */
  void computeDistance(
      const LINALG::Matrix<3,1>&  x,
      const LINALG::Matrix<3,1>&  y,
      LINALG::Matrix<3,3>&        du_tensor_du,
      LINALG::Matrix<3,1>&        dist_vec,
      LINALG::Matrix<3,1>&        dist_unit,
      double&                     distance
      );
  
  
  /*!
  \brief  collect global column indices for stiffness matrices with non-local entries
   */
  void CollectLmcol(
      std::map<int,std::set<int> >&     potentialElementIds,
      vector<int>&                      lmcol
      );
  
  
  //! get local index
  int GetLocalIndex(
      vector<int>&    lmcol,
      int             value
      );
  
  
  //! compute element nodes in reference configuration
  void ReferenceConfiguration(
      const DRT::Element*         element,
      LINALG::SerialDenseMatrix&  X, 
      const int                   numdim
      ) const;
  
  
  //! compute element nodes in spatial configuration
  void SpatialConfiguration(
      const DRT::Element*         element,
      LINALG::SerialDenseMatrix&  x, 
      const int                   numdim
      ) const;
  
  
  //! compute integration factors for surface elements
  double ComputeFactor(
      const DRT::Element*                     element, 
      LINALG::SerialDenseVector&              funct, 
      LINALG::SerialDenseMatrix&              deriv, 
      const DRT::UTILS::IntegrationPoints2D&  intpoints,
      const int                               gp,
      LINALG::Matrix<3,1>&                    x_gp,
      const double                            curve_fac
      );
 
  
  //! compute integration factors for line elements
  double ComputeFactor(
      const DRT::Element*                     element, 
      LINALG::SerialDenseVector&              funct, 
      LINALG::SerialDenseMatrix&              deriv, 
      const DRT::UTILS::IntegrationPoints1D&  intpoints,
      const int                               gp,
      LINALG::Matrix<3,1>&                    x_gp,
      const double                            curve_fac
      );
 
  
  //! name@ Private Attributes
  
  //! problem dimension
  int                                         prob_dim_;

  //!  current discretization
  RCP<DRT::Discretization>                    discretRCP_;

  //!  current discretization
  DRT::Discretization&                        discret_;

  //! SearchTree xTree_
  Teuchos::RCP<GEO::SearchTree>               searchTree_;

  //! current position of boundary elements
  std::map<int,LINALG::Matrix<3,1> >          currentboundarypositions_;

  //! for each potential condition label there is a set of boundary elements
 std::map<int,std::set<int> >                 elementsByLabel_;

  LINALG::MapExtractor                        potboundary_;

  //! @name local copies of input parameters
  Teuchos::RCP<DRT::Discretization>           boundarydis_;
  
  Teuchos::RCP<Teuchos::ParameterList>        params_;
  //@}

  //! @name displacement vector of redundant boundary discretization
  Teuchos::RCP<Epetra_Vector>                 idisp_total_;
  
  //! @name displacement vector of boundary discretization, which is part of the local solid dis
  Teuchos::RCP<Epetra_Vector>                 idisp_onproc_;
  
  Teuchos::RCP<Epetra_Import>                 importer_;


}; //class

}//namespace UTILS

#endif /*DRT_POTENTIAL_MANAGER_H_*/
