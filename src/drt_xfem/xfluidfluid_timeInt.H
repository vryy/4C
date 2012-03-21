

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


  // Create maps of bgfluid's nodes needed for time integration
  void CreateBgNodeMaps(const Teuchos::RCP<DRT::Discretization> bgdis,
                        Teuchos::RCP<XFEM::FluidWizard>         wizard);

  // Save the old state
  void SaveBgNodeMaps();


  int SaveAndCreateNewBgNodeMaps(Teuchos::RCP<DRT::Discretization> bgdis,
                                 Teuchos::RCP<XFEM::FluidWizard>   wizard);

  void SetNewBgStatevectorAndProjectEmbToBg(const Teuchos::RCP<DRT::Discretization>        bgdis,
                                            Teuchos::RCP<Epetra_Vector>           bgstatevn,
                                            Teuchos::RCP<Epetra_Vector>           bgstatevnp,
                                            Teuchos::RCP<Epetra_Vector>           embstatevn,
                                            Teuchos::RCP<Epetra_Vector>           aledispn);

  void SetNewBgStatevectorKeepGhostValues(const Teuchos::RCP<DRT::Discretization>        bgdis,
                                          Teuchos::RCP<Epetra_Vector>           bgstatevn,
                                          Teuchos::RCP<Epetra_Vector>           bgstatevnp,
                                          Teuchos::RCP<Epetra_Vector>           embstatevn,
                                          Teuchos::RCP<Epetra_Vector>           aledispn);

  void SetNewBgStatevectorFullProjection(const Teuchos::RCP<DRT::Discretization>        bgdis,
                                         Teuchos::RCP<Epetra_Vector>           bgstatevn,
                                         Teuchos::RCP<Epetra_Vector>           bgstatevnp,
                                         Teuchos::RCP<Epetra_Vector>           embstatevn,
                                         Teuchos::RCP<Epetra_Vector>           aledispn);

  void SetNewEmbStatevector(const Teuchos::RCP<DRT::Discretization>        bgdis,
                            Teuchos::RCP<Epetra_Vector>    statevbg_n,
                            Epetra_Vector                  statevemb_n,
                            Teuchos::RCP<Epetra_Vector>    statevembnew_n,
                            Teuchos::RCP<Epetra_Vector>    aledispnp,
                            Teuchos::RCP<Epetra_Vector>    aledispnpoldstate);

  bool ComputeSpacialToElementCoordAndProject(DRT::Element*                       pele,
                                              LINALG::Matrix<3,1>&                x,
                                              LINALG::Matrix<4,1>&                interpolatedvec,
                                              Epetra_Vector                       embstate_n,
                                              Teuchos::RCP<Epetra_Vector>         embeddeddisp,
                                              Teuchos::RCP<DRT::Discretization>   sourcedis);

  bool ComputeSpacialToElementCoordAndProject2(DRT::Element*                       pele,
                                               LINALG::Matrix<3,1>&                x,
                                               LINALG::Matrix<4,1>&                interpolatedvec,
                                               Epetra_Vector                       embstate_n,
                                               Teuchos::RCP<Epetra_Vector>         embeddeddisp,
                                               Teuchos::RCP<DRT::Discretization>   sourcedis);

  void GmshOutput(const Teuchos::RCP<DRT::Discretization>        bgdis);

  void PatchelementForIncompressibility(const Teuchos::RCP<DRT::Discretization>     bgdis,
                                        XFEM::FluidWizard                  wizard_n,
                                        XFEM::FluidWizard                  wizard_np,
                                        Teuchos::RCP<LINALG::MapExtractor> dbcmaps);

  void EnforceIncompressibility(const Teuchos::RCP<DRT::Discretization>  bgdis,
                                XFEM::FluidWizard               wizard,
                                Teuchos::RCP<Epetra_Vector>     initialvel);


//  void CreatePatchBoxes(std::map<int, GEO::CUT::BoundingBox> &
//  patchboxes);

  private:

  Teuchos::RCP<DRT::Discretization> embdis_; // bg  discretization

  // maps needed for xfem time integration
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

  const Epetra_Map* currentbgdofmap_;
  const Epetra_Map* oldbgdofmap_;

  // if the status of nodes remains the same don't do any
  // intrepolation from the embedded fluid
  bool samemaps_;

  enum INPAR::XFEM::XFluidFluidTimeInt timeintapproach_;

  // elements of the incompressibility patch
  std::vector<int> incompelementidsAllproc_;

  // Incompressibilitypatch dofrowmap
  // map of all dofs of incompressibility patch
  Teuchos::RCP<const Epetra_Map>  incompressibilitydofrowmap_;
  // map of all velcoity-dofs of incompressibility patch (no pressure)
  Teuchos::RCP<const Epetra_Map>  incompressibilityveldofrowmap_;

  Teuchos::ParameterList        params_;

  int step_;
};

} // namespace XFEM

#endif /*TIMEINTEGRATION_H_*/

