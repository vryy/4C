/*----------------------------------------------------------------------*/
/*!
 * \file objective_funct_surfcurr.H
 *
<pre>
Maintainer: Sebastian Kehl
            kehl@mhpc.mw.tum.de
            089 - 289-10361
</pre>
*/
/*----------------------------------------------------------------------*/

#ifndef OBJECTIVE_FUNCT_SURFCURR_H
#define OBJECTIVE_FUNCT_SURFCURR_H

#include "objective_funct.H"

#include "Teuchos_RCP.hpp"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_Vector.h"
#include "Epetra_Map.h"

// forward declarations
namespace DRT
{
  class Discretization;
  class Condition;
}

namespace LINALG
{
  class SerialDenseMatrix;
}

namespace STR
{

namespace INVANA
{

class MatParManager;
class ObjectiveFunctSurfCurr;

/*!
  \brief An objective function based on a surface current representation

  This is a handler to have easy access to currents based on different conditions,
  i.e every source surface current condition should be evaluated against the corresponding
  target surface current condition.

  Handling of timesteps should eventually be in here to. For now there is only one target
  discretization representing the "final" configuration. So only the final time step contributes
  to the evaluation functional.

  */


class ObjectiveFunctSurfCurrRepresentation : public ObjectiveFunct
{

private:

  //! discretization of the source model
  Teuchos::RCP<DRT::Discretization> sourcedis_;

  //! discretization of the target model
  Teuchos::RCP<DRT::Discretization> targetdis_;

  //! simulated timesteps
  Teuchos::RCP<std::vector<double> > timesteps_;

  //! number of vectors in MultiVectors mdisp_ and disp
  int msteps_;

  //! vector of single surface objective functions
  std::vector<Teuchos::RCP<STR::INVANA::ObjectiveFunctSurfCurr> > currents_;


public:

  //! constructor
  ObjectiveFunctSurfCurrRepresentation(Teuchos::RCP<DRT::Discretization> discret, int steps, Teuchos::RCP<std::vector<double> > timesteps);

  //! standard destructor
  ~ObjectiveFunctSurfCurrRepresentation(){};

  //! set global scale
  void SetScale(double sigmaW);

  //! evaluate value of the objective function
  virtual void Evaluate(Teuchos::RCP<Epetra_MultiVector> disp, double& val);

  //! evaluate gradient of objective function w.r.t to displacements
  virtual void EvaluateGradient(Teuchos::RCP<Epetra_MultiVector> disp, Teuchos::RCP<Epetra_MultiVector> gradient);

};


  /*!
  \brief An objective function based on a surface current representation

  documentation will follow

  */
class ObjectiveFunctSurfCurr
{
private:

  //! discretization of the source model
  Teuchos::RCP<DRT::Discretization> sourcedis_;

  //! discretization of the target model
  Teuchos::RCP<DRT::Discretization> targetdis_;

  //! Surface condition of the source discretization
  DRT::Condition* sourcecond_;

  //! Surface conditions of the target discretization
  DRT::Condition* targetcond_;

  //! dofrowmap of the source discretization
  const Epetra_Map* dofrowmap_;

  //! distributed map of source surfaces facets global ids
  Teuchos::RCP<Epetra_Map> sourcemap_;

  //! distributed map of target surfaces facets global ids
  Teuchos::RCP<Epetra_Map> targetmap_;

  //! distributed map of source surfaces facets global ids
  Teuchos::RCP<Epetra_Map> sourcemapred_;

  //! distributed map of target surfaces facets global ids
  Teuchos::RCP<Epetra_Map> targetmapred_;

  //! data map containing target normals
  //! <facet gid, facet normal <nx, ny, nz> >
  std::map<int,Teuchos::RCP<Epetra_SerialDenseMatrix> > normal_t_;

  //! data map containing target centers
  //! <facet gid, facet center <cx, cy, cz> >
  std::map<int,Teuchos::RCP<Epetra_SerialDenseMatrix> > center_t_;

  //! data map containing source normals
  //! <facet gid, facet normal <nx, ny, nz> >
  std::map<int,Teuchos::RCP<Epetra_SerialDenseMatrix> > normal_s_;

  //! data map containing source normal derivative wrt displacements
  std::map<int,Teuchos::RCP<Epetra_SerialDenseMatrix> > derivnormal_s_;

  //! data map containing source centers
  std::map<int,Teuchos::RCP<Epetra_SerialDenseMatrix> > center_s_;

  //! data map containing each source facets dof gids
  std::map<int, std::vector<int> > facetdofmap_;

  //! squared scale of the kernel
  double sigmaW2_;

  //! scaling of the objective function?
  bool scaling_;

  //! structural integrity component of the target current, which can be preevalauted
  double structinttarget_;

  //! set up maps to hold each facets gid
  const Epetra_Map SetupConditionMap(DRT::Condition* cond, Teuchos::RCP<DRT::Discretization> dis);

  //! Compute normals and centers in a material configuration
  void ComputeNormalCenterMaterialConfig(DRT::Condition* cond,
                                         Teuchos::RCP<DRT::Discretization> discret,
                                         std::map<int,Teuchos::RCP<Epetra_SerialDenseMatrix> >* normals,
                                         std::map<int,Teuchos::RCP<Epetra_SerialDenseMatrix> >* centers);

  //! Compute normals and centers in a spatial configuration
  void ComputeNormalCenterSpatialConfig(DRT::Condition* cond,
                                        Teuchos::RCP<DRT::Discretization> discret,
                                        std::map<int,Teuchos::RCP<Epetra_SerialDenseMatrix> >* normals,
                                        std::map<int,Teuchos::RCP<Epetra_SerialDenseMatrix> >* centers,
                                        std::map<int,Teuchos::RCP<Epetra_SerialDenseMatrix> >* derivnormal,
                                        std::map<int, std::vector<int> >* facetdofmap,
                                        bool deriv=false);

  //! compute normal and centers for a tri facet (i.e there go exactly 3 nodes in here)
  void ComputeNormalCenter(LINALG::SerialDenseMatrix& x, Epetra_SerialDenseMatrix& n, Epetra_SerialDenseMatrix& c);

  //! overload to compute normal and center and derivative of the normal for a tri facet (i.e there go exactly 3 nodes in here)
  void ComputeNormalDeriv(LINALG::SerialDenseMatrix& x,Epetra_SerialDenseMatrix& dn);

  //! Compute partial summation between two maps
  double Convolute(std::map<int,Teuchos::RCP<Epetra_SerialDenseMatrix> >& n1,
                   std::map<int,Teuchos::RCP<Epetra_SerialDenseMatrix> >& c1,
                   std::map<int,Teuchos::RCP<Epetra_SerialDenseMatrix> >& n2,
                   std::map<int,Teuchos::RCP<Epetra_SerialDenseMatrix> >& c2);

  //! Compute derivative of partial summation between two maps, type 1
  void DerivComponent1(Teuchos::RCP<Epetra_Vector> lgradient,
                       std::map<int,std::vector<int> >& facetdofmap,
                       std::map<int,Teuchos::RCP<Epetra_SerialDenseMatrix> >& dn1,
                       std::map<int,Teuchos::RCP<Epetra_SerialDenseMatrix> >& c1,
                       std::map<int,Teuchos::RCP<Epetra_SerialDenseMatrix> >& n2,
                       std::map<int,Teuchos::RCP<Epetra_SerialDenseMatrix> >& c2);

  //! Compute derivative of partial summation between two maps, type 2
  void DerivComponent2(Teuchos::RCP<Epetra_Vector> lgradient,
                       std::map<int,std::vector<int> >& facetdofmap,
                       std::map<int,Teuchos::RCP<Epetra_SerialDenseMatrix> >& n1,
                       std::map<int,Teuchos::RCP<Epetra_SerialDenseMatrix> >& c1,
                       std::map<int,Teuchos::RCP<Epetra_SerialDenseMatrix> >& n2,
                       std::map<int,Teuchos::RCP<Epetra_SerialDenseMatrix> >& c2);

  //! Compute derivative of partial summation between two maps, type 2
  void DerivComponent3(Teuchos::RCP<Epetra_Vector> lgradient,
                       std::map<int,std::vector<int> >& facetdofmap,
                       std::map<int,Teuchos::RCP<Epetra_SerialDenseMatrix> >& n1,
                       std::map<int,Teuchos::RCP<Epetra_SerialDenseMatrix> >& c1,
                       std::map<int,Teuchos::RCP<Epetra_SerialDenseMatrix> >& n2,
                       std::map<int,Teuchos::RCP<Epetra_SerialDenseMatrix> >& c2);

  //! gauss kernel evaluation
  inline double kernel(Epetra_SerialDenseMatrix x, Epetra_SerialDenseMatrix y);

  //! derivativ of the gauss kernel wrt the dofs of the facet x (fid=0) or y (fid=1)
  Epetra_SerialDenseMatrix kernelderiv(Epetra_SerialDenseMatrix x, Epetra_SerialDenseMatrix y, int fid);

  //! dot product of Epetra_SerialDenseMatrices
  inline double dotp(Epetra_SerialDenseMatrix x, Epetra_SerialDenseMatrix y);

  //! debug printout: data map to screen
  void PrintDataToScreen(std::map<int,Teuchos::RCP<Epetra_SerialDenseMatrix> > data);

  //! debug printout: data map to screen
  void PrintDataToScreen(std::map<int,std::vector<int> > data);

  //! debug printout: Epetra Vector to screen with precision
  void PrintDataToScreen(Epetra_Vector& vec);

public:

  //! constructor
  ObjectiveFunctSurfCurr(Teuchos::RCP<DRT::Discretization> sourcedis,
                         Teuchos::RCP<DRT::Discretization> targetdis,
                         DRT::Condition* sourcecond,
                         DRT::Condition* targetcond);

  //! standard destructor
  ~ObjectiveFunctSurfCurr(){};

  //! get scale
  double GetScale(){return sqrt(sigmaW2_);};

  //! set scale
  void SetScale(double sigmaW){sigmaW2_=sigmaW*sigmaW;};

  //! Compute Norm on current difference in W-space
  double WSpaceNorm();

  //! Compute gradient of Norm on current difference in W-space
  void GradientWSpaceNorm(Teuchos::RCP<Epetra_MultiVector> gradient,int vind);


}; // class ObjectiveFunctSurfCurr
}  // namespace INVANA
}  // namespace STR

#endif // OBJECTIVE_FUNCT_SURFCURR_H
