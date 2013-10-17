

/*!-----------------------------------------------------------------------------------------------*
\file xfluidfluid_timeInt.H

\brief provides the basic xfluidfluid time integration

<pre>
Maintainer: Shadan Shahmiri
            shahmiri@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15240
</pre>
 *------------------------------------------------------------------------------------------------*/


#ifndef XFLUIDFLUID_TIMEINTEGRATION_H_
#define XFLUIDFLUID_TIMEINTEGRATION_H_


#include "../linalg/linalg_fixedsizematrix.H"
#include "../drt_inpar/inpar_xfem.H"


class Epetra_Map;
class Epetra_Vector;

namespace DRT
{
  class Discretization;
  class Element;
  class PackBuffer;
  class Exporter;
  class Node;
}

namespace LINALG
{
  class MapExtractor;
}

namespace Teuchos
{
  class ParameterList;
}

namespace XFEM
{
  class FluidWizard;

/*!
\brief this class is the basic TIMEINT class for the projection, adaption or
       something else in XFEM-problems between consecutive time steps
 */
class XFluidFluidTimeIntegration
{
public:

  //! constructor
  explicit XFluidFluidTimeIntegration(
    const Teuchos::RCP<DRT::Discretization> bgdis,
    const Teuchos::RCP<DRT::Discretization> embdis,
    Teuchos::RCP<XFEM::FluidWizard>         wizard,
    int                            step,
    enum INPAR::XFEM::XFluidFluidTimeInt xfem_timeintapproach,
    const Teuchos::ParameterList&              params
    );

  //! destructor
  ~XFluidFluidTimeIntegration()
  {
    return;
  }


  /// Create maps of bgfluid's nodes needed for time integration
  void CreateBgNodeMaps(const Teuchos::RCP<DRT::Discretization> bgdis,
                        Teuchos::RCP<XFEM::FluidWizard>         wizard);

  /// Save the old state
  void SaveBgNodeMaps();

  /// Save the old maps of bg nodes and create new map of bg nodes
  int SaveAndCreateNewBgNodeMaps(Teuchos::RCP<DRT::Discretization> bgdis,
                                 Teuchos::RCP<XFEM::FluidWizard>   wizard);

  void CreateBgNodeMapsForRestart(Teuchos::RCP<DRT::Discretization> bgdis,
                                  Teuchos::RCP<XFEM::FluidWizard>   wizard);

  void SetNewBgStatevectorAndProjectEmbToBg(const Teuchos::RCP<DRT::Discretization>   bgdis,
		  	  	  	  	  	  	  	  	  	Teuchos::RCP<Epetra_Vector>               bgstatevn1,
                                            Teuchos::RCP<Epetra_Vector>               bgstatevnp1,
                                            Teuchos::RCP<Epetra_Vector>               embstatevn1,
                                            Teuchos::RCP<Epetra_Vector>           	  bgstatevn2,
                                            Teuchos::RCP<Epetra_Vector>           	  bgstatevnp2,
                                            Teuchos::RCP<Epetra_Vector>           	  embstatevn2,
                                            Teuchos::RCP<Epetra_Vector>               bgstatevn3,
                                            Teuchos::RCP<Epetra_Vector>               bgstatevnp3,
                                            Teuchos::RCP<Epetra_Vector>               embstatevn3,
                                            Teuchos::RCP<Epetra_Vector>               bgstatevn4,
                                            Teuchos::RCP<Epetra_Vector>               bgstatevnp4,
                                            Teuchos::RCP<Epetra_Vector>               embstatevn4,
                                            Teuchos::RCP<Epetra_Vector>               bgstatevn5,
                                            Teuchos::RCP<Epetra_Vector>               bgstatevnp5,
                                            Teuchos::RCP<Epetra_Vector>               embstatevn5,
                                            Teuchos::RCP<Epetra_Vector>               aledispn);

  /// set values to background state vector. Do not do the projection if
  /// ghost-values are available.
  void SetNewBgStatevectorKeepGhostValues(const Teuchos::RCP<DRT::Discretization>     bgdis,
		  	  	  	  	  	  	  	  	  Teuchos::RCP<Epetra_Vector>                 bgstatevn1,
		  	  	  	  	  	  	  	  	  Teuchos::RCP<Epetra_Vector>                 bgstatevnp1,
		  	  	  	  	  	  	  	  	  Teuchos::RCP<Epetra_Vector>                 embstatevn1,
		  	  	  	  	  	  	  	  	  Teuchos::RCP<Epetra_Vector>                 bgstatevn2,
		  	  	  	  	  	  	  	  	  Teuchos::RCP<Epetra_Vector>                 bgstatevnp2,
		  	  	  	  	  	  	  	  	  Teuchos::RCP<Epetra_Vector>           	  embstatevn2,
		  	  	  	  	  	  	  	  	  Teuchos::RCP<Epetra_Vector>                 bgstatevn3,
		  	  	  	  	  	  	  	  	  Teuchos::RCP<Epetra_Vector>                 bgstatevnp3,
		  	  	  	  	  	  	  	  	  Teuchos::RCP<Epetra_Vector>                 embstatevn3,
		  	  	  	  	  	  	  	  	  Teuchos::RCP<Epetra_Vector>                 bgstatevn4,
		  	  	  	  	  	  	  	  	  Teuchos::RCP<Epetra_Vector>                 bgstatevnp4,
		  	  	  	  	  	  	  	  	  Teuchos::RCP<Epetra_Vector>                 embstatevn4,
		  	  	  	  	  	  	  	  	  Teuchos::RCP<Epetra_Vector>     	          bgstatevn5,
		  	  	  	  	  	  	  	  	  Teuchos::RCP<Epetra_Vector>                 bgstatevnp5,
		  	  	  	  	  	  	  	  	  Teuchos::RCP<Epetra_Vector>                 embstatevn5,
		  	  	  	  	  	  	  	  	  Teuchos::RCP<Epetra_Vector>                 aledispn);

  /// set values to background state vector. Always project values from
  /// embedded discretization if no history values are available.
  void SetNewBgStatevectorFullProjection(const Teuchos::RCP<DRT::Discretization>      bgdis,
		  	  	  	  	  	  	  	  	 Teuchos::RCP<Epetra_Vector>                  bgstatevn1,
		  	  	  	  	  	  	  	  	 Teuchos::RCP<Epetra_Vector>                  bgstatevnp1,
		  	  	  	  	  	  	  	  	 Teuchos::RCP<Epetra_Vector>                  embstatevn1,
		  	  	  	  	  	  	  	  	 Teuchos::RCP<Epetra_Vector>                  bgstatevn2,
		  	  	  	  	  	  	  	  	 Teuchos::RCP<Epetra_Vector>                  bgstatevnp2,
		  	  	  	  	  	  	  	  	 Teuchos::RCP<Epetra_Vector>                  embstatevn2,
		  	  	  	  	  	  	  	  	 Teuchos::RCP<Epetra_Vector>                  bgstatevn3,
		  	  	  	  	  	  	  	  	 Teuchos::RCP<Epetra_Vector>                  bgstatevnp3,
		  	  	  	  	  	  	  	  	 Teuchos::RCP<Epetra_Vector>                  embstatevn3,
		  	  	  	  	  	  	  	  	 Teuchos::RCP<Epetra_Vector>                  bgstatevn4,
		  	  	  	  	  	  	  	  	 Teuchos::RCP<Epetra_Vector>                  bgstatevnp4,
		  	  	  	  	  	  	  	  	 Teuchos::RCP<Epetra_Vector>                  embstatevn4,
		  	  	  	  	  	  	  	  	 Teuchos::RCP<Epetra_Vector>                  bgstatevn5,
		  	  	  	  	  	  	  	  	 Teuchos::RCP<Epetra_Vector>                  bgstatevnp5,
		  	  	  	  	  	  	  	  	 Teuchos::RCP<Epetra_Vector>                  embstatevn5,
		  	  	  	  	  	  	  	  	 Teuchos::RCP<Epetra_Vector>                  aledispn);

  void SetNewEmbStatevector(const Teuchos::RCP<DRT::Discretization>  bgdis,
                            Teuchos::RCP<Epetra_Vector>              statevbg_n,
                            Teuchos::RCP<Epetra_Vector>              statevemb_n,
                            Teuchos::RCP<Epetra_Vector>              statevembnew_n,
                            Teuchos::RCP<Epetra_Vector>              aledispnp,
                            Teuchos::RCP<Epetra_Vector>              aledispnpoldstate);

  /// Write the values of node from bgstatevn to bgstatevnp
  void WriteValuestoBgStateVector(const Teuchos::RCP<DRT::Discretization>   bgdis,
                                  DRT::Node*                                bgnode,
                                  std::vector<int>                          gdofs_n,
                                  Teuchos::RCP<Epetra_Vector>               bgstatevnp,
                                  Teuchos::RCP<Epetra_Vector>               bgstatevn);

  ///  Check whether the unknown node lies in the patch-element.If yes fill the
  ///  interpolatedvec and return true.
  bool ComputeSpacialToElementCoordAndProject(DRT::Element*                       pele,
                                              LINALG::Matrix<3,1>&                x,
                                              LINALG::Matrix<20,1>&               interpolatedvec,
                                              Teuchos::RCP<Epetra_Vector>         embstate_n1,
                                              Teuchos::RCP<Epetra_Vector>         embstate_n2,
                                              Teuchos::RCP<Epetra_Vector>         embstate_n3,
                                              Teuchos::RCP<Epetra_Vector>         embstate_n4,
                                              Teuchos::RCP<Epetra_Vector>         embstate_n5,
                                              Teuchos::RCP<Epetra_Vector>         embeddeddisp,
                                              Teuchos::RCP<DRT::Discretization>   sourcedis);

  bool ComputeSpacialToElementCoordAndProject2(DRT::Element*                       pele,
                                               LINALG::Matrix<3,1>&                x,
                                               LINALG::Matrix<4,1>&                interpolatedvec,
                                               Teuchos::RCP<Epetra_Vector>         embstate_n,
                                               Teuchos::RCP<Epetra_Vector>         embeddeddisp,
                                               Teuchos::RCP<DRT::Discretization>   sourcedis);


  void GmshOutput(const Teuchos::RCP<DRT::Discretization>        bgdis);

  void GmshOutputForInterpolateFSI(const Teuchos::RCP<DRT::Discretization>    bgdis,
                                   Teuchos::RCP<Epetra_Vector>                aledispnp,
                                   Teuchos::RCP<Epetra_Vector>                aledispnpoldstate);

  void PatchelementForIncompressibility(const Teuchos::RCP<DRT::Discretization>  bgdis,
                                        Teuchos::RCP<XFEM::FluidWizard>          wizard_n,
                                        Teuchos::RCP<XFEM::FluidWizard>          wizard_np,
                                        Teuchos::RCP<LINALG::MapExtractor>       dbcmaps);

  /// build an incompressibility discretization
  void PrepareIncompDiscret(const Teuchos::RCP<DRT::Discretization>     bgdis,
                            Teuchos::RCP<XFEM::FluidWizard>             wizard_np);


  void EvaluateIncompressibility(const Teuchos::RCP<DRT::Discretization>  bgdis,
                                 Teuchos::RCP<XFEM::FluidWizard>          wizard);

  void SolveIncompOptProb(Teuchos::RCP<Epetra_Vector>     initialvel);


  /// Transfer data in the round robin communication pattern
  void CommunicateNodes(const Teuchos::RCP<DRT::Discretization>    bgdis,
                        std::vector<LINALG::Matrix<3,1> >        & bgnodes_coords,
                        std::vector<LINALG::Matrix<20,1> >       & interpolated_vecs,
                        std::vector<int>                         & bgnodeidwithnohistory,
                        Teuchos::RCP<Epetra_Vector>                aledispn,
                        Teuchos::RCP<Epetra_Vector>     	   embstatevn1,
                        Teuchos::RCP<Epetra_Vector>         	   bgstatevnp1,
                        Teuchos::RCP<Epetra_Vector>		   bgstatevn1,
                        Teuchos::RCP<Epetra_Vector>                embstatevn2,
                        Teuchos::RCP<Epetra_Vector>                bgstatevnp2,
                        Teuchos::RCP<Epetra_Vector>                bgstatevn2,
                        Teuchos::RCP<Epetra_Vector>                embstatevn3,
                        Teuchos::RCP<Epetra_Vector>                bgstatevnp3,
                        Teuchos::RCP<Epetra_Vector>                bgstatevn3,
                        Teuchos::RCP<Epetra_Vector>                embstatevn4,
                        Teuchos::RCP<Epetra_Vector>                bgstatevnp4,
                        Teuchos::RCP<Epetra_Vector>                bgstatevn4,
                        Teuchos::RCP<Epetra_Vector>                embstatevn5,
                        Teuchos::RCP<Epetra_Vector>                bgstatevnp5,
                        Teuchos::RCP<Epetra_Vector>                bgstatevn5);

  /// receive a block in the round robin communication pattern
  void ReceiveBlock(std::vector<char>  & rblock,
                    DRT::Exporter      & exporter,
                    MPI_Request        & request);

  /// send a block in the round robin communication pattern
  void SendBlock(std::vector<char>  & sblock,
                 DRT::Exporter      & exporter,
                 MPI_Request        & request );

  /// pack values in the round robin communication pattern
  void PackValues(std::vector<LINALG::Matrix<3,1> > & bgnodes_coords,
                  std::vector<LINALG::Matrix<20,1> >& interpolatedvec,
                  std::vector<int>                  & bgnodeidwithnohistory,
                  std::vector<int>                  & NodeDone,
                  std::vector<char>                 & sblock);


  void FindEmbeleAndInterpolatevalues(std::vector<LINALG::Matrix<3,1> > & bgnodes_coords,
		  	  	  	  	  	  	  	  std::vector<LINALG::Matrix<20,1> > & interpolated_vecs,
		  	  	  	  	  	  	  	  std::vector<int>                  & NodeDone,
		  	  	  	  	  	  	  	  Teuchos::RCP<Epetra_Vector>         aledispn,
		  	  	  	  	  	  	  	  Teuchos::RCP<Epetra_Vector>         embstatevn1,
		  	  	  	  	  	  	  	  Teuchos::RCP<Epetra_Vector>         embstatevn2,
		  	  	  	  	  	  	  	  Teuchos::RCP<Epetra_Vector>         embstatevn3,
		  	  	  	  	  	  	  	  Teuchos::RCP<Epetra_Vector>         embstatevn4,
		  	  	  	  	  	  	  	  Teuchos::RCP<Epetra_Vector>         embstatevn5);

  void SearchRadius();

  private:

  Teuchos::RCP<DRT::Discretization> embdis_; // embedded  discretization
  Teuchos::RCP<DRT::Discretization> incompdis_; // incompressibility  discretization

  /// maps needed for xfem time integration
  std::map<int, std::vector<int> > stdnoden_;
  std::map<int, std::vector<int> > stdnodenp_;

  std::map<int, std::vector<int> > enrichednoden_;
  std::map<int, std::vector<int> > enrichednodenp_;

  // map of every node to it's parent elements. This could be more
  // than one if we have more than one volume cells for the node
  std::map<int, std::vector<int> > nodetoparentele_;

  // map of the parent elements (from the map above) to their nds
  // vectos
  std::map<int, std::vector< std::vector< int > > > parenteletondsset_;

  std::set<int> projectednodeids_;

  // current processor id
  int myrank_;

  int numproc_;

  const Epetra_Map* currentbgdofmap_;
  const Epetra_Map* oldbgdofmap_;

  // if the status of nodes remains the same don't do any
  // interpolation from the embedded fluid
  bool samemaps_;

  enum INPAR::XFEM::XFluidFluidTimeInt timeintapproach_;

  Teuchos::RCP<Epetra_Vector> C_; // incompressibility vector

  // incompressibility elements and nodes
  std::set<int> incompnodeids_set_;
  std::set<int> incompelementids_set_;

  Teuchos::ParameterList        params_;

  int step_;

  bool gmsh_debug_out_;

  // min radius needed for the search tree
  double minradius_;

  // search radius factor
  double searchradius_fac_;

};

} // namespace XFEM

#endif /*TIMEINTEGRATION_H_*/

