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

<pre>
Maintainer: Peter Gamnitzer
            gamnitzer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>

*/


#ifdef CCADISCRET
#ifndef DYN_SMAG_H
#define DYN_SMAG_H


#include "../drt_lib/drt_discret.H"
#include "../drt_lib/linalg_utils.H"


#include "Epetra_Vector.h"

#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "Teuchos_ParameterList.hpp"

using namespace std;
using namespace Teuchos;

namespace FLD 
{

class DynSmagFilter
{

public:

  /*!
  \brief Standard Constructor (public)

  */
  DynSmagFilter(
    RCP<DRT::Discretization>     actdis             ,
    RCP<map<int,vector<int> > >  pbcmapmastertoslave,
    ParameterList&               params             );

  /*!
  \brief Destructor
  
  */
  virtual ~DynSmagFilter();

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
    Teuchos::RCP<Epetra_Vector>           solution    ,
    Teuchos::RCP<Epetra_Vector>           dirichtoggle
    );

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
  void DynSmagBoxFilter(
    Teuchos::RCP<Epetra_Vector>           velocity    ,
    Teuchos::RCP<Epetra_Vector>           dirichtoggle
    );

  /*!
  \brief Compute Cs using the filtered quantities.
  This is an element call!

  For a turbulent channel flow, the averaging of the Smagorinsky 
  constant is done in here.
  */
  void DynSmagComputeCs();

  //! @name input arguments of the constructor
  //
  //! the discretization
  RCP<DRT::Discretization>     discret_;
  //! connection between master and slave nodes on this proc
  RCP<map<int,vector<int> > >  pbcmapmastertoslave_;
  //! parameterlist including time params, stabilization params and turbulence sublist
  ParameterList&               params_ ;
  //@}

  //! @name control parameters
  bool                         channel_flow_;
  bool                         apply_dynamic_smagorinsky_;
  //@}


  //! @name vectors used for filtering (for dynamic Smagorinsky model)
  //        --------------------------

  //! the box filtered velocities in nodes (3 vectors)
  RCP<Epetra_MultiVector>      filtered_vel_;
  //! the box filtered reynoldsstresses in nodes (9 vectors)
  RCP<Epetra_MultiVector>      filtered_reynoldsstress_;
  //! the modeled subgrid stress in nodes (9 vectors)
  RCP<Epetra_MultiVector>      filtered_modeled_subgrid_stress_;
  //! the filtered vel exported to column map
  RCP<Epetra_MultiVector>      col_filtered_vel_;
  //! the filtered reystress exported to column map
  RCP<Epetra_MultiVector>      col_filtered_reynoldsstress_;
  //! the modeled subgrid stresses exported to column map
  RCP<Epetra_MultiVector>      col_filtered_modeled_subgrid_stress_;
  //@}

  //! @name turbulent channel flow specials
  //        -------------------------------

  //! the plane coordinates for the above mentioned averaging procedure
  RCP<vector<double> >         planecoords_;
  //@}

}; // end class DynSmagFilter

} // end namespace FLD

#endif // not DYN_SMAG_H

#endif
