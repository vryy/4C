/*!----------------------------------------------------------------------
\file drt_potential.H

\brief  Base class controlling surface stresses due to potential forces
        between mesoscopic structures

<pre>
Maintainer: Ursula Mayer
            mayer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15257
</pre>

*----------------------------------------------------------------------*/

#ifndef DRT_POTENTIAL_H_
#define DRT_POTENTIAL_H_

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


namespace POTENTIAL
{
/*!
\brief  Base class controlling surface and volume stresses due
        to potential forces and containing all necessary history data
*/

class Potential
{
public:

  /*!
    \brief Standard Constructor

    \param discret (in): Discretization
   */
  Potential(
      Teuchos::RCP<DRT::Discretization>     discretRCP,
      DRT::Discretization&                  discret
      );


  //! copy constructor
  Potential(const POTENTIAL::Potential& old);


  /*!
    \brief Destructor
   */
  ~Potential() {};



protected:


  /*!
  \brief  set history variables after reading restart. The row map
          based vectors written in case of restart are exported to column
          based ones needed for calculations again.
   */
  void searchElementsInCutOffRadius(
    const Teuchos::RCP<DRT::Discretization>     potentialdis,
    const std::map<int,LINALG::Matrix<3,1> >&   currentpositions,
    const LINALG::Matrix<3,1>&                  point,
    std::map<int,std::set<int> >&               potentialElementIds,
    const double                                radius
    );


  /*!
  \brief  set history variables after reading restart. The row map
          based vectors written in case of restart are exported to column
          based ones needed for calculations again.
   */
  void treeSearchElementsInCutOffRadius(
    const Teuchos::RCP<DRT::Discretization>     potentialdis,
    const std::map<int,LINALG::Matrix<3,1> >&   currentpositions,
    const LINALG::Matrix<3,1>&                  point,
    std::map<int,std::set<int> >&               potentialElementIds,
    const double                                radius,
    const int                                   label
    );
  
  
  /*!
  \brief  set history variables after reading restart. The row map
          based vectors written in case of restart are exported to column
          based ones needed for calculations again. 
   */
  void treeSearchElementsInCutOffRadius(
    const Teuchos::RCP<DRT::Discretization>   potentialdis,
    std::map<int,LINALG::Matrix<3,2> >&       elemXAABBList,
    const DRT::Element*                       element,
    std::map<int,std::set<int> >&             potentialElementIds,
    const double                              radius,
    const int                                 label
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
  \brief  select correct potential out of condition
   */
  void EvaluatePotentialfromCondition(
      RefCountPtr<DRT::Condition>   cond,
      const LINALG::Matrix<3,1>&    x,
      const LINALG::Matrix<3,1>&    y,
      LINALG::Matrix<3,1>&          potderiv1,
      LINALG::Matrix<3,3>&          potderiv2
      );
  
  
  /*!
   \brief  select correct potential out of condition
   */
  void EvaluatePotentialfromCondition(
      RefCountPtr<DRT::Condition>   cond,
      const LINALG::Matrix<2,1>&    x,
      const LINALG::Matrix<2,1>&    y,
      LINALG::Matrix<2,1>&          potderiv1,
      LINALG::Matrix<2,2>&          potderiv2
      );
  
  
  /*!
    \brief  select correct potential out of condition
     */
  void EvaluatePotentialfromCondition_Approx1(
      RefCountPtr<DRT::Condition>   cond,
      const LINALG::Matrix<3,1>&    x,
      const LINALG::Matrix<3,1>&    y,
      LINALG::Matrix<3,1>&          potderiv1,
      LINALG::Matrix<3,3>&          potderiv2
      );
  
  /*!
      \brief  select correct potential out of condition
       */
    void EvaluatePotentialfromCondition_Approx2(
        RefCountPtr<DRT::Condition>   cond,
        const LINALG::Matrix<3,1>&    x,
        const LINALG::Matrix<3,1>&    y,
        LINALG::Matrix<3,1>&          Fs,
        LINALG::Matrix<3,3>&          Fsderiv
        );


  /*!
  \brief  compute first and second derivative of the Lennard-Jones potential 
          for 3D problems
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
    \brief  compute first and second derivative of the Lennard-Jones potential 
            for 2D problems
   */
  void EvaluateLennardJonesPotential(
      const double                  depth,
      const double                  rootDist,
      const LINALG::Matrix<2,1>&    x,
      const LINALG::Matrix<2,1>&    y,
      LINALG::Matrix<2,1>&          potderiv1,
      LINALG::Matrix<2,2>&          potderiv2
      );
  
  
  /*!
     \brief  compute first and second derivative of the approximated Lennard-Jones  potential
     */
  void EvaluateLennardJonesPotential_Approx1(
      const double                  depth,
      const double                  rootDist,
      const LINALG::Matrix<3,1>&    x,
      const LINALG::Matrix<3,1>&    y,
      LINALG::Matrix<3,1>&          potderiv1,
      LINALG::Matrix<3,3>&          potderiv2
      );
  
  /*!
       \brief  compute first and second derivative of the approximated Lennard-Jones  potential
       */
    void EvaluateLennardJonesPotential_Approx2(
        const double                  depth,
        const double                  rootDist,
        const LINALG::Matrix<3,1>&    x,
        const LINALG::Matrix<3,1>&    y,
        LINALG::Matrix<3,1>&          Fs,
        LINALG::Matrix<3,3>&          Fsderiv
        );


  /*!
  \brief  compute first and second derivative of the Zeta potential
          for 3D problems
   */
  void EvaluateElectrostaticRepulsion(
    const double                  zeta_param_1,
    const double                  zeta_param_2,
    const LINALG::Matrix<3,1>&    x,
    const LINALG::Matrix<3,1>&    y,
    LINALG::Matrix<3,1>&          potderiv1,
    LINALG::Matrix<3,3>&          potderiv2
    );
  
  
  /*!
    \brief  compute first and second derivative of the Zeta potential
            for 2D problems
   */
  void EvaluateElectrostaticRepulsion(
    const double                  zeta_param_1,   //depth
    const double                  zeta_param_2,   //rootdist
    const LINALG::Matrix<2,1>&    x,
    const LINALG::Matrix<2,1>&    y,
    LINALG::Matrix<2,1>&          potderiv1,
    LINALG::Matrix<2,2>&          potderiv2
    );


  /*!
    \brief  compute first and second derivative of the Van-Der-Waals potential
            for 3D problems
   */
  void EvaluateVanDerWaals(
      const double                  lambda,
      const LINALG::Matrix<3,1>&    x,
      const LINALG::Matrix<3,1>&    y,
      LINALG::Matrix<3,1>&          potderiv1,
      LINALG::Matrix<3,3>&          potderiv2
      );
  
  
  /*!
    \brief  compute first and second derivative of the Van-Der-Waals potential
            for 3D problems (volume approximation)
   */
  void EvaluateVanDerWaals_Approx1(
      const double                  lambda,
      const LINALG::Matrix<3,1>&    x,
      const LINALG::Matrix<3,1>&    y,
      LINALG::Matrix<3,1>&          potderiv1,
      LINALG::Matrix<3,3>&          potderiv2
  );
  
  /*!
      \brief  compute first and second derivative of the Van-Der-Waals potential
              for 3D problems (volume approximation)
     */
    void EvaluateVanDerWaals_Approx2(
        const double                  lambda,
        const LINALG::Matrix<3,1>&    x,
        const LINALG::Matrix<3,1>&    y,
        LINALG::Matrix<3,1>&          Fs,
        LINALG::Matrix<3,3>&          Fsderiv
    );
  
  
  /*!
  \brief  computes distance related structures 3D
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
    \brief  computes distance related structures 2D
     */
  void computeDistance(
      const LINALG::Matrix<2,1>&  x,
      const LINALG::Matrix<2,1>&  y,
      LINALG::Matrix<2,2>&        du_tensor_du,
      LINALG::Matrix<2,1>&        dist_vec,
      LINALG::Matrix<2,1>&        dist_unit,
      double&                     distance
      );


  /*!
  \brief  collect global column indices for stiffness matrices with non-local entries
   */
  void CollectLmcol(
      const Teuchos::RCP<DRT::Discretization>     potentialdis,
      std::map<int,std::set<int> >&               potentialElementIds,
      vector<int>&                                lmcol
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
      const std::map<int,LINALG::Matrix<3,1> >&   currentpositions,
      const DRT::Element*                         element,
      LINALG::SerialDenseMatrix&                  x,
      const int                                   numdim
      ) const;


  //! name@ Protected Attributes

  //!  current discretization
  RCP<DRT::Discretization>                    discretRCP_;

  //!  current discretization
  DRT::Discretization&                        discret_;

  //! SearchTree xTree_
  Teuchos::RCP<GEO::SearchTree>               searchTree_;

  //! problem dimension
  int                                         prob_dim_;



private:

  //don't want = operator and cctor
  Potential operator = (const POTENTIAL::Potential& old);


}; //class potential

}//namespace UTILS

#endif /*DRT_POTENTIAL_H_*/

