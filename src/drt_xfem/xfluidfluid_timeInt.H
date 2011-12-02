

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
#include "../drt_cut/cut_position.H"
#include "../drt_io/io_ostream0.H"
#include "../drt_io/io.H"
#include "../drt_io/io_gmsh.H"
#include "../drt_fluid/fluid_utils.H"
#include "../drt_lib/drt_condition_selector.H"

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
      XFEM::FluidWizard              wizard,
      int                            step
    );

  //! destructor
  ~XFluidFluidTimeIntegration()
  {
    return;
  }


  // Create maps of bgfluid's nodes needed for time integration
  void CreateBgNodeMaps();

  // Save the old state
  void SaveBgNodeMaps();

  void SaveAndCreateNewBgNodeMaps();

  void SetNewBgStatevectorAndProjectEmbToBg(Teuchos::RCP<Epetra_Vector>           bgstatevn,
                                            Teuchos::RCP<Epetra_Vector>           bgstatevnp,
                                            Teuchos::RCP<Epetra_Vector>           embstatevn,
                                            Teuchos::RCP<Epetra_Vector>           aledispn);

  void SetNewEmbStatevector(Teuchos::RCP<Epetra_Vector>    statevbg_n,
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

  void GmshOutput();

//  void CreatePatchBoxes(std::map<int, GEO::CUT::BoundingBox> &
//  patchboxes);

  private:

  RCP<DRT::Discretization> bgdis_; // bg  discretization
  RCP<DRT::Discretization> embdis_; // bg  discretization


  // maps needed for xfem time integration
  std::map<int, vector<int> > stdnoden_;
  std::map<int, vector<int> > stdnodenp_;

  std::map<int, vector<int> > enrichednoden_;
  std::map<int, vector<int> > enrichednodenp_;

  XFEM::FluidWizard wizard_;

  int step_;
};

} // namespace XFEM

#endif /*TIMEINTEGRATION_H_*/
#endif // CCADISCRET
