/*!

\file dyn_vreman.H

\brief Filter methods for the dynamic Smagorinsky model

Reference:

<pre>
    Subgrid viscosity:
    D. You, P. Moin:
    A dynamic global-coefficient
    subgrid-scale eddy-viscosity model for large-eddy simulation
    in complex geometries
    (Phys. Fluids 2007)

    Subgrid diffusivity:
    D. You, P. Moin:
    A dynamic global-coefficient subgrid-scale odel for large-eddy simulation
    of turbulent scalar transport in complex geometries
    (Phys. Fluids 2009)

    Remark:
    alphaij=transpose(velderxy)
    rows and colums of the jacobian alphaij are swapped in comparison
    to the jacobian velderxy used in baci

</pre>

<pre>
Maintainer: Benjamin Krank
            krank@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
</pre>

*/



#ifndef DYN_VREMAN_H
#define DYN_VREMAN_H


#include "../drt_lib/drt_discret.H"
#include "../linalg/linalg_utils.H"


#include "Epetra_Vector.h"

#include "Teuchos_RCP.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "Teuchos_ParameterList.hpp"

#include "../drt_inpar/inpar_fluid.H"
#include "../drt_inpar/inpar_scatra.H"



namespace FLD
{
class Boxfilter;

class Vreman
{

public:

  /*!
  \brief Standard Constructor (public)

  */
  Vreman(    RCP<DRT::Discretization>     actdis             ,
      RCP<std::map<int,std::vector<int> > >  pbcmapmastertoslave,
      Teuchos::ParameterList&      params);

  /*!
  \brief Destructor

  */
  virtual ~Vreman();



  void ApplyFilterForDynamicComputationOfCv(
      const Teuchos::RCP<const Epetra_Vector>             velocity,
      const Teuchos::RCP<const Epetra_Vector>             scalar,
      const double                                        thermpress,
      const Teuchos::RCP<const Epetra_Vector>             dirichtoggle);

  void ApplyFilterForDynamicComputationOfDt(
    const Teuchos::RCP<const Epetra_MultiVector>        velocity,
    const Teuchos::RCP<const Epetra_Vector>             scalar,
    const double                                        thermpress,
    const Teuchos::RCP<const Epetra_Vector>             dirichtoggle,
    Teuchos::ParameterList&                             extraparams
    );

  void AddScatra(
    RCP<DRT::Discretization>     scatradis,
    INPAR::SCATRA::ScaTraType    scatratype,
    RCP<std::map<int,std::vector<int> > >  scatra_pbcmapmastertoslave);

  void GetCv(double& Cv)
  {
    Cv=Cv_;
    return ;
  }

  double Cv_;

private:

  /// provide access to the box filter
  Teuchos::RCP<FLD::Boxfilter> Boxfilter();
  //Boxfilter
  Teuchos::RCP<FLD::Boxfilter> Boxf_;
  Teuchos::RCP<FLD::Boxfilter> Boxfsc_;



  //! the discretization
  RCP<DRT::Discretization>     discret_;
  //! connection between master and slave nodes on this proc
  RCP<std::map<int,std::vector<int> > >  pbcmapmastertoslave_;
  //! parameterlist including time params, stabilization params and turbulence sublist
  Teuchos::ParameterList&      params_ ;
  //! flag for physical type of fluid flow
  INPAR::FLUID::PhysicalType   physicaltype_;
  //scatra specific
  RCP<DRT::Discretization>     scatradiscret_;
  INPAR::SCATRA::ScaTraType    scatratype_;
  RCP<std::map<int,std::vector<int> > >  scatra_pbcmapmastertoslave_;


  double DynVremanComputeCv();

  void DynVremanComputeDt(Teuchos::ParameterList& extraparams);

  //@}

  //! @name vectors used for filtering (for dynamic Smagorinsky model)
  //        --------------------------

  //! the filtered reystress exported to column map
  RCP<Epetra_MultiVector>      col_filtered_strainrate_;
  RCP<Epetra_Vector>      col_filtered_expression_;
  RCP<Epetra_MultiVector>      col_filtered_alphaij_;
  RCP<Epetra_Vector>      col_filtered_alpha2_;
  RCP<Epetra_MultiVector>      col_filtered_phi_;
  RCP<Epetra_Vector>      col_filtered_phi2_;
  RCP<Epetra_Vector>      col_filtered_phiexpression_;
  RCP<Epetra_MultiVector>      col_filtered_alphaijsc_;
  //@}

  //! @name homogeneous flow specials
  //        -------------------------------

  //@}

}; // end class Vreman

} // end namespace FLD

#endif // not DYN_VREMAN_H

