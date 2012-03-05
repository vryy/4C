

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

#ifdef CCADISCRET

#ifndef XFLUIDFLUID_TIMEINTEGRATION_H_
#define XFLUIDFLUID_TIMEINTEGRATION_H_


#include "xfem_fluidwizard.H"
#include "../drt_io/io_ostream0.H"
#include "../drt_io/io.H"
#include "../drt_io/io_gmsh.H"
#include "../drt_fluid/fluid_utils.H"

namespace XFEM
{
/*!
\brief this class is the basic TIMEINT class for the projection, adaption or
       something else in XFEM-problems between consecutive time steps
 */
class XFluidFluidTimeIntegration
{
public:

  //! constructor
  explicit XFluidFluidTimeIntegration(
    const RCP<DRT::Discretization> bgdis,
    const RCP<DRT::Discretization> embdis,
    RCP<XFEM::FluidWizard>         wizard,
    int                            step,
    enum INPAR::XFEM::XFluidFluidTimeInt xfem_timeintapproach
    );

  //! destructor
  ~XFluidFluidTimeIntegration()
  {
    return;
  }


  // Create maps of bgfluid's nodes needed for time integration
  void CreateBgNodeMaps(const RCP<DRT::Discretization> bgdis,
                        RCP<XFEM::FluidWizard>         wizard);

  // Save the old state
  void SaveBgNodeMaps();

  int SaveAndCreateNewBgNodeMaps(RCP<DRT::Discretization> bgdis,
                                 RCP<XFEM::FluidWizard>   wizard);

  void SetNewBgStatevectorAndProjectEmbToBg(const RCP<DRT::Discretization>        bgdis,
                                            Teuchos::RCP<Epetra_Vector>           bgstatevn,
                                            Teuchos::RCP<Epetra_Vector>           bgstatevnp,
                                            Teuchos::RCP<Epetra_Vector>           embstatevn,
                                            Teuchos::RCP<Epetra_Vector>           aledispn);

  void SetNewBgStatevectorKeepGhostValues(const RCP<DRT::Discretization>        bgdis,
                                          Teuchos::RCP<Epetra_Vector>           bgstatevn,
                                          Teuchos::RCP<Epetra_Vector>           bgstatevnp,
                                          Teuchos::RCP<Epetra_Vector>           embstatevn,
                                          Teuchos::RCP<Epetra_Vector>           aledispn);

  void SetNewBgStatevectorFullProjection(const RCP<DRT::Discretization>        bgdis,
                                         Teuchos::RCP<Epetra_Vector>           bgstatevn,
                                         Teuchos::RCP<Epetra_Vector>           bgstatevnp,
                                         Teuchos::RCP<Epetra_Vector>           embstatevn,
                                         Teuchos::RCP<Epetra_Vector>           aledispn);

  void SetNewEmbStatevector(const RCP<DRT::Discretization>        bgdis,
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

  void GmshOutput(const RCP<DRT::Discretization>        bgdis);

//  void CreatePatchBoxes(std::map<int, GEO::CUT::BoundingBox> &
//  patchboxes);

  private:

  RCP<DRT::Discretization> embdis_; // bg  discretization

  // maps needed for xfem time integration
  std::map<int, vector<int> > stdnoden_;
  std::map<int, vector<int> > stdnodenp_;

  std::map<int, vector<int> > enrichednoden_;
  std::map<int, vector<int> > enrichednodenp_;

  // map of every node to it's parent elements. This could be more
  // than one if we have more than one volume cells for the node
  std::map<int, vector<int> > nodetoparentele_;

  // map of the parent elements (from the map above) to their nds
  // vectos
  std::map<int, std::vector< std::vector< int > > > parenteletondsset_;

  const Epetra_Map* currentbgdofmap_;
  const Epetra_Map* oldbgdofmap_;

  // if the status of nodes remains the same don't do any
  // intrepolation from the embedded fluid
  bool samemaps_;

  enum INPAR::XFEM::XFluidFluidTimeInt timeintapproach_;

  int step_;
};

} // namespace XFEM

#endif /*TIMEINTEGRATION_H_*/
#endif // CCADISCRET
