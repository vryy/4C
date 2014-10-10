/*!----------------------------------------------------------------------
\file fluid_xwall.H
\brief

<pre>
Maintainer: Benjamin Krank
            krank@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
</pre>

*----------------------------------------------------------------------*/

#ifndef FLD_XWALL_H
#define FLD_XWALL_H

#include <Teuchos_RCP.hpp>
#include <Epetra_Vector.h>
#include <Teuchos_ParameterList.hpp>
#include "../drt_inpar/inpar_fluid.H"

//forward declarations
namespace DRT
{
  class Discretization;
}

namespace LINALG
{
  class SparseMatrix;
  class Solver;
}


namespace FLD
{
class FluidImplicitTimeInt;

class XWall
{

public:
  /// Standard Constructor
  XWall(Teuchos::RCP<DRT::Discretization>      dis,
      int                           nsd,
      Teuchos::RCP<Teuchos::ParameterList>& params);

  /// Destructor
  virtual ~XWall(){return;};

  // set element params for xwall EnrichmentType
  void SetXWallParams(Teuchos::ParameterList& eleparams);

  // adapt ml nullspace for aggregation for scale separation (mfs)
  void AdaptMLNullspace(const Teuchos::RCP<LINALG::Solver>&   solver);

  // smooth vector via aggregation
  void FilterResVectorForOutput(Teuchos::RCP<Epetra_Vector>   residual);

  // get output vector of enriched dofs
  Teuchos::RCP<Epetra_Vector> GetOutputVector(Teuchos::RCP<Epetra_Vector>   vel);

  // zero out enriched dofs of a row vector
  void ZeroEnrDofs(Teuchos::RCP<Epetra_Vector>   vec);

  // returns, if Properties for GenAlpha have to be updated
  void UpdateTauW(int step,Teuchos::RCP<Epetra_Vector>   trueresidual,int itnum,Teuchos::RCP<Epetra_Vector>   accn,Teuchos::RCP<Epetra_Vector>   velnp,Teuchos::RCP<Epetra_Vector>   veln, FluidImplicitTimeInt& fluid);

protected:

  // set element params for xwall EnrichmentType, distributed for xwall discretization
  void SetXWallParamsXWDis(Teuchos::ParameterList& eleparams);

  // setup XWall
  void Setup();

  // initialize some maps
  void InitXWallMaps();

  // initialize a toggle vector for the element
  void InitToggleVector();

  // initialize wall distance
  void InitWallDist();

  // setup xwall discretization
  void SetupXWallDis();

  // setup l2 projection
  void SetupL2Projection();

  // calculate wall shear stress
  void CalcTauW(int step, Teuchos::RCP<Epetra_Vector>   trueresidual, Teuchos::RCP<Epetra_Vector>   velnp,Teuchos::RCP<Epetra_Vector>   wss);

  // l2 project vectors
  void L2ProjectVector(Teuchos::RCP<Epetra_Vector>   veln,Teuchos::RCP<Epetra_Vector>   velnp,Teuchos::RCP<Epetra_Vector>   accn);

  // l2 project vectors with continuity constraint
  void L2ProjectVectorWithContinuityConstraint(Teuchos::RCP<Epetra_Vector>   veln,Teuchos::RCP<Epetra_Vector>   velnp,Teuchos::RCP<Epetra_Vector>   accn);

  // calculate parameter for stabilization parameter mk
  void CalcMK();

  // make enrichment function constant on interface resulting in a natural blending
  void BlendingViaModificationOfTauw(Teuchos::RCP<Epetra_Vector>   newtauw2);

private:

  //! discretisation
  Teuchos::RCP<DRT::Discretization> discret_;

  //! fluid params
  Teuchos::RCP<Teuchos::ParameterList> params_;

  //! the processor ID from the communicator
  int                           myrank_;

  //! map including all wall nodes, redundant map
  Teuchos::RCP<Epetra_Map> dircolnodemap_;

  //! xwall node row map
  Teuchos::RCP<Epetra_Map> xwallrownodemap_;

  //! xwall node col map
  Teuchos::RCP<Epetra_Map> xwallcolnodemap_;

  //! map including the enriched dofs, row map
  Teuchos::RCP<Epetra_Map> enrdofrowmap_;

  //! map including the unused pressure dofs, row map
  Teuchos::RCP<Epetra_Map> lagrdofrowmap_;

  //! map including all enriched dofs plus unused pressure dofs
  Teuchos::RCP<Epetra_Map> mergedmap_;

  //! wall distance, local vector
  Teuchos::RCP<Epetra_Vector> walldist_;

  //! wall distance, standard node row map
  Teuchos::RCP<Epetra_Vector> wdist_;

  //! wall distance, row map of redistributed discretization
  Teuchos::RCP<Epetra_Vector> wdistxwdis_;

  //! vector on same dofs for tauw
  Teuchos::RCP<Epetra_Vector> tauw_;

  //! vector on same dofs for tauw
  Teuchos::RCP<Epetra_Vector> tauwxwdis_;

  //! vector on same dofs for deltatauw (increment)
  Teuchos::RCP<Epetra_Vector> inctauw_;

  //! vector on same dofs for deltatauw (increment)
  Teuchos::RCP<Epetra_Vector> inctauwxwdis_;

  //! matrix projecting the wall shear stress to off-wall nodes
  Teuchos::RCP<LINALG::SparseMatrix> tauwcouplingmattrans_;

  //! toggle vector, standard node row map
  Teuchos::RCP<Epetra_Vector> xwalltoggle_;

  //! toggle vector, xwall discretization
  Teuchos::RCP<Epetra_Vector> xwalltogglexwdis_;

  //! toggle vector, local vector
  Teuchos::RCP<Epetra_Vector> xtoggleloc_;

  //! redistributed xwall discretization
  Teuchos::RCP<DRT::Discretization> xwdiscret_;

  //! mass matrix for projection
  Teuchos::RCP<LINALG::SparseMatrix> massmatrix_;

  //! solver for projection
  Teuchos::RCP<LINALG::Solver> solver_;

  //! increment of veln during projection
  Teuchos::RCP<Epetra_Vector> incveln_;

  //! increment of velnp during projection
  Teuchos::RCP<Epetra_Vector> incvelnp_;

  //! increment of accn during projection
  Teuchos::RCP<Epetra_Vector> incaccn_;

  //! veln for state of xwall discretization during projection
  Teuchos::RCP<Epetra_Vector> stateveln_;

  //! velnp for state of xwall discretization during projection
  Teuchos::RCP<Epetra_Vector> statevelnp_;

  //! accn for state of xwall discretization during projection
  Teuchos::RCP<Epetra_Vector> stateaccn_;

  //! MK for standard discretization
  Teuchos::RCP<Epetra_Vector> mkstate_;

  //! MK for xwall discretization
  Teuchos::RCP<Epetra_Vector> mkxwstate_;

  //! scale separation matrix for smoothing of tauw/ residual
  Teuchos::RCP<LINALG::SparseMatrix> SepEnr_;

  //! increment factor of tauw
  double fac_;

  //! increment of tauw
  double inctauwnorm_;

  //! constant tauw from input file
  double constant_tauw_;

  //! minimum tauw from input file
  double min_tauw_;

  //! number of gauss points in wall-parallel and normal direction
  int gp_norm_;
  int gp_par_;

  //! penalty parameter
  double penalty_param_;

  //! viscosity
  double visc_;

  //! density
  double dens_;

  //! when and how to update tauw
  enum INPAR::FLUID::XWallTauwType tauwtype_;

  //! how to calculate tauw
  enum INPAR::FLUID::XWallTauwCalcType tauwcalctype_;

  //! how to blend
  enum INPAR::FLUID::XWallBlendingType blendingtype_;

  //! projection
  bool proj_;

  //! projection with continuity constraint
  bool projconstr_;

  //! smoothing through aggregation of residual
  bool smooth_res_aggregation_;

  //! switch from gradient-based to residual-based calculation of tauw
  int switch_step_;

}; // class XWall

} // namespace FLD


#endif  // #ifndef FLD_XWALL_H
