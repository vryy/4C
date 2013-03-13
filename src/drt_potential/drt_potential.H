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

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>

#include <Epetra_Vector.h>
#include <Epetra_Map.h>
#include <Epetra_Operator.h>
#include "Epetra_Import.h"

#include "../linalg/linalg_sparsematrix.H"
#include "../drt_lib/drt_discret.H"
#include "../linalg/linalg_utils.H"
#include "../linalg/linalg_serialdensematrix.H"
#include "../linalg/linalg_serialdensevector.H"
#include "../drt_io/io.H"
#include "../drt_geometry/searchtree.H"
#include "../drt_geometry/intersection_service.H"
#include "../drt_fem_general/drt_utils_integration.H"
#include "../drt_potential/drt_potential_dofset.H"
#include "../drt_potential/drt_potential_container.H"



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
  \brief  select correct potential out of condition
   */
  void EvaluatePotentialfromCondition(
      RCP<DRT::Condition>   cond,
      const LINALG::Matrix<3,1>&    x,
      const LINALG::Matrix<3,1>&    y,
      LINALG::Matrix<3,1>&          potderiv1,
      LINALG::Matrix<3,3>&          potderiv2
      );
  
  
  /*!
   \brief  select correct potential out of condition
   */
  void EvaluatePotentialfromCondition(
      RCP<DRT::Condition>   cond,
      const LINALG::Matrix<2,1>&    x,
      const LINALG::Matrix<2,1>&    y,
      LINALG::Matrix<2,1>&          potderiv1,
      LINALG::Matrix<2,2>&          potderiv2
      );
  
  
  /*!
    \brief  select correct potential out of condition
     */
  void EvaluatePotentialfromCondition_Approx1(
      RCP<DRT::Condition>   cond,
      const LINALG::Matrix<3,1>&    x,
      const LINALG::Matrix<3,1>&    y,
      LINALG::Matrix<3,1>&          potderiv1,
      LINALG::Matrix<3,3>&          potderiv2
      );
  
  /*!
      \brief  select correct potential out of condition
       */
    void EvaluatePotentialfromCondition_Approx2(
        RCP<DRT::Condition>   cond,
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
            for 2D problems (plane strain approximation)
   */
  void EvaluateVanDerWaals(
      const double                  lambda,
      const LINALG::Matrix<2,1>&    x,
      const LINALG::Matrix<2,1>&    y,
      LINALG::Matrix<2,1>&          potderiv1,
      LINALG::Matrix<2,2>&          potderiv2
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
      std::vector<int>&                           lmcol
      );
  
  
  /*!
    \brief  collect global column indices for stiffness matrices with loacl and non-local entries
     */
  void CollectLmcol(
      const Teuchos::RCP<DRT::Discretization>                       potentialdis,
      std::map<int,std::set<int> >&                                 potentialElementIds,
      std::map<int,std::vector<PotentialElementContainer> >&        nonlocalPecs,
      std::vector<int>&                                             lmcol
      );
  
  
  void CollectLmcol(
      const Teuchos::RCP<DRT::Discretization>                       potentialdis,
      const std::map< int, std::map<int, GEO::NearestObject> >&     potentialObjectsAtGP,
      std::vector<int> &                                            lmcol
      );
  
  
  int GetElementId( 
      const Teuchos::RCP<DRT::Discretization>         potentialdis,
      const GEO::NearestObject&                       potObject
      );


  //! get local index
  int GetLocalIndex(
      std::vector<int>& lmcol,
      int               value
      );
  
  
  //! compute element nodes in reference configuration
  void ReferenceConfiguration(
      const DRT::Element*         element,
      Epetra_SerialDenseMatrix&  X,
      const int                   numdim
      ) const;


   //! compute element nodes in reference configuration
  Epetra_SerialDenseMatrix ReferenceConfiguration(
    const DRT::Element*                         element,
    const int                                   numdim
    ) const;
  
  
  //! compute element nodes in spatial configuration
  void SpatialConfiguration(
      const std::map<int,LINALG::Matrix<3,1> >&   currentpositions,
      const DRT::Element*                         element,
      Epetra_SerialDenseMatrix&                   x,
      const int                                   numdim
      ) const;


  //! compute element nodes in spatial configuration
  Epetra_SerialDenseMatrix SpatialConfiguration(
    const std::map<int,LINALG::Matrix<3,1> >&   currentpositions,
    const DRT::Element*                         element,
    const int                                   numdim
    ) const;

   
  //! get time curve factor
  double GetTimeCurveFactor(
    Teuchos::ParameterList&                        		params
	);
  
  
  //! invert elements by label
  void InvertElementsPerLabel(
      std::map<int,std::set<int> >&	elementsByLabel,
      std::map<int,int >&				labelByElements
  	  );
  
  
  ///////////////////////  Test methods ///////////////////////
  
  void computeTestVanDerWaalsSpheres(
    const Teuchos::RCP<DRT::Discretization> potentialsurfdis,
    const std::map<int,std::set<int> >&     elementsByLabel_Vol,
    const std::map<int,std::set<int> >&     elementsByLabel_Surf,
    const RCP<const Epetra_Vector>  disp,
    const RCP<Epetra_Vector>        fint,
    const double                            time,
    const int                               step,
    const double                            vdw_radius,
    const double                            n_offset
    );
  
  
  void computeTestVanDerWaalsMembranes(
    const Teuchos::RCP<DRT::Discretization> potentialsurfdis,
    const std::map<int,std::set<int> >&     elementsByLabel_Vol,
    const std::map<int,std::set<int> >&     elementsByLabel_Surf,
    const RCP<const Epetra_Vector>  disp,
    const RCP<Epetra_Vector>        fint,
    const double                            time,
    const int                               step,
    const double                            vdw_radius,
    const double                            n_offset,
    const double                            thickness
    );
  
  
  // compute local forces and center of gravity
  double computeLocalForceAndCOG(  
    const Teuchos::RCP<DRT::Discretization> potentialsurfdis,
    LINALG::Matrix<3,1>&        fpot_sphere,
    LINALG::Matrix<3,1>&        cog_sphere,
    const RCP<Epetra_Vector>    fint,
    const RCP<const Epetra_Vector>    disp,
    const std::set<int>&        elementIds,
    const std::set<int>&        surfElementIds
    );
  
  
  // compute global forces and center of gravity
  double computeGlobalForceAndCOG(
    double                vol_sphere_local, 
    LINALG::Matrix<3,1>&  fpot_sphere_local,
    LINALG::Matrix<3,1>&  cog_sphere_local,
    LINALG::Matrix<3,1>&  fpot_sphere_global,
    LINALG::Matrix<3,1>&  cog_sphere_global
    );
  
  
  //! get atomic density
  double GetAtomicDensity(
      const int                       elementId,
      const string&                   conditionName,
      const std::map<int,int >&       labelByElements
      );
  
  
  //! get physical eleCoords
  void getPhysicalEleCoords(
      Teuchos::RCP<DRT::Discretization>     dis,  
      Teuchos::RCP<const Epetra_Vector>     idisp_solid,
      const DRT::Element*                   element,
      LINALG::SerialDenseMatrix&            xyze
      );
  
  
  //! write output
  void WriteTestOutput(
    const double        distance,
    const double        force1,
    const double        force2,
    const double        force_analytical,
    const double        time,
    const int           step,
    const std::string   name
    );
  
  
  ///////////////////////  Test methods ///////////////////////


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

