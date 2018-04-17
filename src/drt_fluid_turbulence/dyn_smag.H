/*!

\file dyn_smag.H

\brief Filter methods for the dynamic Smagorinsky model

References are

<pre>

    M. Germano, U. Piomelli, P. Moin, W.H. Cabot:
    A dynamic subgrid-scale eddy viscosity model
    (Phys. Fluids 1991)

    or

    D.K. Lilly:
    A proposed modification of the Germano subgrid-scale closure method
    (Phys. Fluids 1992)

    or
    A.E. Tejada-Martinez
    Dynamic subgrid-scale modeling for large eddy simulation of turbulent
    flows with a stabilized finite element method
    (Phd thesis, Rensselaer Polytechnic Institute, Troy, New York)

</pre>

\level 2

<pre>
\maintainer Martin Kronbichler
            kronbichler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>

*/



#ifndef DYN_SMAG_H
#define DYN_SMAG_H


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

class DynSmagFilter
{

public:

  /*!
  \brief Standard Constructor (public)

  */
  DynSmagFilter(
    Teuchos::RCP<DRT::Discretization>     actdis,
    Teuchos::ParameterList&               params);

  /*!
  \brief Destructor

  */
  virtual ~DynSmagFilter();

  void AddScatra(
    Teuchos::RCP<DRT::Discretization>     scatradis
  );

  /*!
  \brief Perform box filter operation, compare filtered quantities
  to solution to get an estimate for Cs (using clpping), average
  over element layers in turbulent channel flows.

  This method initialises element quantities (standard case) or
  provides information for the element via the parameter list
  (in plane averaging for channel flow)

  \param solution     (in) velocity field to filter and to
                           determine Cs from
  \param dirichtoggle (in) information on dirichlet dofs to be
                           able to exclude boundary nodes from
                           filtering

  */
  void ApplyFilterForDynamicComputationOfCs(
    const Teuchos::RCP<const Epetra_Vector>             velocity,
    const Teuchos::RCP<const Epetra_Vector>             scalar,
    const double                                        thermpress,
    const Teuchos::RCP<const Epetra_Vector>             dirichtoggle
    );

  void ApplyFilterForDynamicComputationOfPrt(
    const Teuchos::RCP<const Epetra_Vector>             scalar,
    const double                                        thermpress,
    const Teuchos::RCP<const Epetra_Vector>             dirichtoggle,
    Teuchos::ParameterList&                             extraparams,
    const int                                           ndsvel
    );


  /*!
  \brief Output of averaged velocity vector for paraview IO

  \param outvec  (in/out) vector in dofrowmap-format to use for
                          output of averaged solution


  */





private:
  /*!
  \brief perform box filtering in five steps

  1) Integrate element Heaviside functions against the quantities
     which are filtered. Add the result to the nodevectors
     (we get a contribution for every node of the element)
     This is an element call!
  2) send/add values from slaves to masters
  3) zero out dirichlet nodes
  4) do normalization by division by the patchvolume
     (Heaviside function -> box filter function)
  5) Communication part: Export filtered quantities from row to
     column map

     \param velocity     (i) the velocity defining the
                             unfiltered quantities
     \param dirichtoggle (i) specifying which nodes have to be
                             set to zero

  */

  /// provide access to the box filter
  Teuchos::RCP<FLD::Boxfilter> Boxfilter();
  /*!
  \brief Compute Cs using the filtered quantities.
  This is an element call!

  For a turbulent channel flow, the averaging of the Smagorinsky
  constant is done in here.
  */
  void DynSmagComputeCs();

  /*!
  \brief Compute Prt using the filtered quantities.
  This is an element call!

  For a turbulent channel flow, the averaging of the turbulent
  Prandtl number is done in here.
  */
  void DynSmagComputePrt(Teuchos::ParameterList& extraparams,
                         int& numele_layer);


  //! @name input arguments of the constructor
  //

  //Boxfilter
  Teuchos::RCP<FLD::Boxfilter> Boxf_;
  Teuchos::RCP<FLD::Boxfilter> Boxfsc_;

  //! the discretization
  Teuchos::RCP<DRT::Discretization>     discret_;
  //! parameterlist including time params, stabilization params and turbulence sublist
  Teuchos::ParameterList&      params_ ;
  //! flag for physical type of fluid flow
  INPAR::FLUID::PhysicalType   physicaltype_;
  //@}

  //! @name control parameters
  bool                         homdir_;
  std::string                       special_flow_homdir_;
  bool                         apply_dynamic_smagorinsky_;
  bool                         calc_Ci_;
  //@}

  //! @name special scatra variables
  //! the discretization
  Teuchos::RCP<DRT::Discretization>     scatradiscret_;
  //@}

  //! @name vectors used for filtering (for dynamic Smagorinsky model)
  //        --------------------------

  //! the filtered vel exported to column map
  Teuchos::RCP<Epetra_MultiVector>      col_filtered_vel_;
  //! the filtered reystress exported to column map
  Teuchos::RCP<Epetra_MultiVector>      col_filtered_reynoldsstress_;
  //! the modeled subgrid stresses exported to column map
  Teuchos::RCP<Epetra_MultiVector>      col_filtered_modeled_subgrid_stress_;
  //! the filtered velocities times rho exported to column map
  Teuchos::RCP<Epetra_MultiVector>      col_filtered_dens_vel_;
  //! the filtered density exported to column map
  Teuchos::RCP<Epetra_Vector>           col_filtered_dens_;
  //! the filtered strainrate times rho exported to column map
  Teuchos::RCP<Epetra_Vector>           col_filtered_dens_strainrate_;
  //! the modeled fine scale velocities exported to column map
  Teuchos::RCP<Epetra_MultiVector>      col_fs_vel_;
  //! the filtered density times temperature times velocity exported to column map (scalar)
  Teuchos::RCP<Epetra_MultiVector>      col_filtered_dens_vel_temp_;
  //! the filtered density times temperature gradient times rate of strain exported to column map (scalar)
  Teuchos::RCP<Epetra_MultiVector>      col_filtered_dens_rateofstrain_temp_;
//  //! the filtered temperature gradient exported to column map (scalar)
//  Teuchos::RCP<Epetra_MultiVector>      col_filtered_gradtemp_;
  //! the filtered temperature exported to column map (scalar)
  Teuchos::RCP<Epetra_Vector>           col_filtered_temp_;
  //! the filtered density times temperature exported to column map (scalar)
  Teuchos::RCP<Epetra_Vector>           col_filtered_dens_temp_;
  //@}

  //! @name homogeneous flow specials
  //        -------------------------------

  //! the direction coordinates for the above mentioned averaging procedure
  Teuchos::RCP<std::vector<double> >         dir1coords_;
  Teuchos::RCP<std::vector<double> >         dir2coords_;
  //@}

}; // end class DynSmagFilter

} // end namespace FLD

#endif // not DYN_SMAG_H

