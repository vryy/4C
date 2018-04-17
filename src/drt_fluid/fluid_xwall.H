/*!----------------------------------------------------------------------
\file fluid_xwall.H
\brief

\level 2

<pre>
\maintainer Martin Kronbichler
            kronbichler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15235/-235
</pre>

*----------------------------------------------------------------------*/

#ifndef FLD_XWALL_H
#define FLD_XWALL_H

#include <Teuchos_RCP.hpp>
#include <Epetra_Vector.h>
#include <Teuchos_ParameterList.hpp>
#include "../drt_inpar/inpar_fluid.H"
#include "../linalg/linalg_utils.H"

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

namespace IO
{
  class DiscretizationReader;
}

namespace FLD
{
class FluidImplicitTimeInt;
class TransferTurbulentInflowConditionNodal;
namespace UTILS
{
class StressManager;
}

class XWall
{

public:
  /// Standard Constructor
  XWall(Teuchos::RCP<DRT::Discretization>      dis,
      int                           nsd,
      Teuchos::RCP<Teuchos::ParameterList>& params,
      Teuchos::RCP<LINALG::MapExtractor> dbcmaps,
      Teuchos::RCP<FLD::UTILS::StressManager> wssmanager);

  /// Destructor
  virtual ~XWall(){return;};

  // set element params for xwall EnrichmentType
  virtual void SetXWallParams(Teuchos::ParameterList& eleparams);

  // adapt ml nullspace for aggregation for scale separation (mfs)
  void AdaptMLNullspace(const Teuchos::RCP<LINALG::Solver>&   solver);

  // get output vector of enriched dofs
  Teuchos::RCP<Epetra_Vector> GetOutputVector(Teuchos::RCP<Epetra_Vector>   vel);

  // returns, if Properties for GenAlpha have to be updated
  virtual void UpdateTauW(int step,Teuchos::RCP<Epetra_Vector>   trueresidual,int itnum,Teuchos::RCP<Epetra_Vector>   accn,Teuchos::RCP<Epetra_Vector>   velnp,Teuchos::RCP<Epetra_Vector>   veln);

  // returns tauw of discret_
  Teuchos::RCP<Epetra_Vector>   GetTauw();

  Teuchos::RCP<Epetra_Vector> GetTauwVector()
      {
        Teuchos::RCP<Epetra_Vector> tauw=Teuchos::rcp(new Epetra_Vector(*(discret_->NodeRowMap()),true));
        LINALG::Export(*tauw_,*tauw);
        return tauw;
      }

  //read restart including wall stresses
  void ReadRestart(IO::DiscretizationReader& reader);

  //fix residual at Dirichlet-inflow nodes such that the wss can be calculated
  Teuchos::RCP<Epetra_Vector>  FixDirichletInflow(Teuchos::RCP<Epetra_Vector>   trueresidual);

  //set current non-linear iteration number
  void SetIter(int iter){iter_=iter;return;}

protected:

  // set element params for xwall EnrichmentType, distributed for xwall discretization
  virtual void SetXWallParamsXWDis(Teuchos::ParameterList& eleparams);

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
  void CalcTauW(int step, Teuchos::RCP<Epetra_Vector>   velnp,Teuchos::RCP<Epetra_Vector>   wss);

  // l2 project vectors
  void L2ProjectVector(Teuchos::RCP<Epetra_Vector>   veln,Teuchos::RCP<Epetra_Vector>   velnp,Teuchos::RCP<Epetra_Vector>   accn);

  // calculate parameter for stabilization parameter mk
  void CalcMK();

  // for inflowchannel
  void TransferAndSaveTauw();

  // for inflowchannel
  void OverwriteTransferredValues();

  //! discretisation
  Teuchos::RCP<DRT::Discretization> discret_;

  //! fluid params
  Teuchos::RCP<Teuchos::ParameterList> params_;

  //! manager for wall shear stress
  Teuchos::RCP<FLD::UTILS::StressManager> mystressmanager_;

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

  //! vector on same dofs for old tauw
  Teuchos::RCP<Epetra_Vector> oldtauw_;

  //! vector on same dofs for old tauw
  Teuchos::RCP<Epetra_Vector> oldinctauw_;

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

  Teuchos::RCP<Epetra_Vector> restart_wss_;

  Teuchos::RCP<TransferTurbulentInflowConditionNodal> turbulent_inflow_condition_;

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
  int gp_norm_ow_;
  int gp_par_;

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

  //! smoothing through aggregation of residual
  bool smooth_res_aggregation_;

  //! mfs on fine scales
  bool fix_residual_on_inflow_;

  //! switch from gradient-based to residual-based calculation of tauw
  int switch_step_;

  int iter_;

}; // class XWall

class XWallAleFSI: public XWall
{
public:
  /// Standard Constructor
  XWallAleFSI(Teuchos::RCP<DRT::Discretization>      dis,
      int                           nsd,
      Teuchos::RCP<Teuchos::ParameterList>& params,
      Teuchos::RCP<LINALG::MapExtractor> dbcmaps,
      Teuchos::RCP<FLD::UTILS::StressManager> wssmanager,
      Teuchos::RCP<Epetra_Vector>   dispnp,
      Teuchos::RCP<Epetra_Vector>   gridv);

  /// Destructor
  virtual ~XWallAleFSI(){return;};

  void UpdateWDistWALE();

  // set element params for xwall EnrichmentType
  void SetXWallParams(Teuchos::ParameterList& eleparams);

  // returns, if Properties for GenAlpha have to be updated
  void UpdateTauW(int step,Teuchos::RCP<Epetra_Vector>   trueresidual,int itnum,Teuchos::RCP<Epetra_Vector>   accn,Teuchos::RCP<Epetra_Vector>   velnp,Teuchos::RCP<Epetra_Vector>   veln);

private:
  // set element params for xwall EnrichmentType, distributed for xwall discretization
  void SetXWallParamsXWDis(Teuchos::ParameterList& eleparams);
  Teuchos::RCP<Epetra_Vector>   mydispnp_;
  Teuchos::RCP<Epetra_Vector>   mygridv_;

  //! wall distance, row map of redistributed discretization
  Teuchos::RCP<Epetra_Vector> incwdistxwdis_;
};

} // namespace FLD


#endif  // #ifndef FLD_XWALL_H
