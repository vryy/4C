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
Maintainer: Ursula Rasthofer
            grasthofer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
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

class DynSmagFilter
{

public:

  /*!
  \brief Standard Constructor (public)

  */
  DynSmagFilter(
    RCP<DRT::Discretization>     actdis             ,
    RCP<map<int,std::vector<int> > >  pbcmapmastertoslave,
    Teuchos::ParameterList&      params             );

  /*!
  \brief Destructor

  */
  virtual ~DynSmagFilter();

  void AddScatra(
    RCP<DRT::Discretization>     scatradis,
    INPAR::SCATRA::ScaTraType    scatratype,
    RCP<map<int,std::vector<int> > >  scatra_pbcmapmastertoslave);

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
    Teuchos::RCP<Epetra_Vector>             velocity,
    Teuchos::RCP<Epetra_Vector>             scalar,
    const double                            thermpress,
    const Teuchos::RCP<const Epetra_Vector> dirichtoggle
    );

  void ApplyFilterForDynamicComputationOfPrt(
    Teuchos::RCP<Epetra_MultiVector>        velocity,
    Teuchos::RCP<Epetra_Vector>             scalar,
    const double                            thermpress,
    const Teuchos::RCP<const Epetra_Vector> dirichtoggle,
    Teuchos::ParameterList&                          extraparams
    );

  void ApplyFilter(
          Teuchos::RCP<Epetra_Vector>             velocity,
          Teuchos::RCP<Epetra_Vector>             scalar,
          const double                            thermpress,
    const Teuchos::RCP<const Epetra_Vector> dirichtoggle
    );


  /*!
  \brief Output of averaged velocity vector for paraview IO

  \param outvec  (in/out) vector in dofrowmap-format to use for
                          output of averaged solution


  */
  void OutputofAveragedVel(
    Teuchos::RCP<Epetra_Vector> outvec
    )
    {
      // loop all elements on this proc (including ghosted ones)
      for (int n=0;n<discret_->NumMyRowNodes();++n)
      {
        // get the node
        DRT::Node* node = discret_->lRowNode(n);

        vector<int> dofs= discret_->Dof(node);
        
        for(int d=0;d<discret_->NumDof(node)-1;++d)
        {
          double val=(*((*filtered_vel_)(d)))[n];
          
          int    id =dofs[d];

          (outvec)->ReplaceGlobalValues(1,&val,&id);
        }
      }
      return;
    }


  void OutputofFineScaleVel(
    Teuchos::RCP<Epetra_Vector> outvec
    )
    {
      // loop all elements on this proc (including ghosted ones)
      for (int n=0;n<discret_->NumMyRowNodes();++n)
      {
        // get the node
        DRT::Node* node = discret_->lRowNode(n);

        vector<int> dofs= discret_->Dof(node);

        for(int d=0;d<discret_->NumDof(node)-1;++d)
        {
          double val=(*((*fs_vel_)(d)))[n];

          int    id =dofs[d];

          (outvec)->ReplaceGlobalValues(1,&val,&id);
        }
      }
      return;
    }


  void FilteredVelComp(
    Teuchos::RCP<Epetra_Vector> outvec,
    const int                   i,
    const int                   j
    )
    {
      // loop all elements on this proc (including ghosted ones)
      for (int nid=0;nid<discret_->NumMyRowNodes();++nid)
      {
        // get the node
        DRT::Node* node = discret_->lRowNode(nid);
        int id = node->Id();

        double vali=(*((*filtered_vel_)(i-1)))[nid];
        double valj=(*((*filtered_vel_)(j-1)))[nid];
        double val = vali*valj;

        (outvec)->ReplaceGlobalValues(1,&val,&id);
      }
      return;
    }

  /*!
  \brief Output of averaged reynoldsstress vector for paraview IO (one component)

  \param outvec  (in/out) vector in dofrowmap-format to use for
                          output of averaged solution
  \param i       which component?
  \param j       which component?


  */
  void FilteredReyStrComp(
    Teuchos::RCP<Epetra_Vector> outvec,
    const int                   i,
    const int                   j
    )
    {
      // loop all elements on this proc (including ghosted ones)
      int ij = 3*(i-1)+(j-1);
      for (int nid=0;nid<discret_->NumMyRowNodes();++nid)
      {
        // get the node
        DRT::Node* node = discret_->lRowNode(nid);
        int id = node->Id();

        double val = (*((*filtered_reynoldsstress_) (ij)))[nid];

        (outvec)->ReplaceGlobalValues(1,&val,&id);
      }
      return;
    }


  void GetFilteredVelocity(
    Teuchos::RCP<Epetra_MultiVector> velocity
    )
    {
      // loop all elements on this proc (including ghosted ones)
      for (int nid=0;nid<discret_->NumMyColNodes();++nid)
      {
        for (int i =0;i<3;++i)
        {
          double val = (*((*col_filtered_vel_)(i)))[nid];
          int err = ((*velocity)(i))->ReplaceMyValues(1,&val,&nid);
          if (err!=0) dserror("dof not on proc");
        }
      }
      return;
    }


  void GetFilteredReynoldsStress(
    Teuchos::RCP<Epetra_MultiVector> reystre
    )
    {
      // loop all elements on this proc (including ghosted ones)
      for (int nid=0;nid<discret_->NumMyColNodes();++nid)
      {
        for (int i =0;i<3;++i)
        {
          for (int j =0;j<3;++j)
          {
             const int ij = 3*i+j;
             double val = (*((*col_filtered_reynoldsstress_) (ij)))[nid];
             int err = ((*reystre) (ij))->ReplaceMyValues(1,&val,&nid);
             if (err!=0) dserror("dof not on proc");
          }
        }
      }
      return;
    }


  void GetFineScaleVelocity(
    Teuchos::RCP<Epetra_MultiVector> velocity
    )
    {
      // loop all elements on this proc (including ghosted ones)
      for (int nid=0;nid<discret_->NumMyColNodes();++nid)
      {
        for (int i =0;i<3;++i)
        {
          double val = (*((*col_fs_vel_)(i)))[nid];
          int err = ((*velocity)(i))->ReplaceMyValues(1,&val,&nid);
          if (err!=0) dserror("dof not on proc");
        }
      }
      return;
    }


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

  void ApplyBoxFilter(
    const Teuchos::RCP<const Epetra_Vector> velocity,
    const Teuchos::RCP<const Epetra_Vector> scalar,
    const double                            thermpress,
    const Teuchos::RCP<const Epetra_Vector> dirichtoggle
    );

  void ApplyBoxFilterScatra(
    const Teuchos::RCP<const Epetra_MultiVector> velocity,
    const Teuchos::RCP<const Epetra_Vector> scalar,
    const double                            thermpress,
    const Teuchos::RCP<const Epetra_Vector> dirichtoggle
    );

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
  //! the discretization
  RCP<DRT::Discretization>     discret_;
  //! connection between master and slave nodes on this proc
  RCP<map<int,std::vector<int> > >  pbcmapmastertoslave_;
  //! parameterlist including time params, stabilization params and turbulence sublist
  Teuchos::ParameterList&      params_ ;
  //! flag for physical type of fluid flow
  INPAR::FLUID::PhysicalType   physicaltype_;
  //@}

  //! @name control parameters
  bool                         homdir_;
  string                       special_flow_homdir_;
  bool                         apply_dynamic_smagorinsky_;
  bool                         apply_box_filter_;
  //@}

  //! @name special scatra variables
  //! the discretization
  RCP<DRT::Discretization>     scatradiscret_;
  INPAR::SCATRA::ScaTraType    scatratype_;
  RCP<map<int,std::vector<int> > >  scatra_pbcmapmastertoslave_;
  //@}

  //! @name vectors used for filtering (for dynamic Smagorinsky model)
  //        --------------------------

  //! the box filtered velocities in nodes (3 vectors)
  RCP<Epetra_MultiVector>      filtered_vel_;
  //! the box filtered reynoldsstresses in nodes (9 vectors)
  RCP<Epetra_MultiVector>      filtered_reynoldsstress_;
  //! the modeled subgrid stress in nodes (9 vectors)
  RCP<Epetra_MultiVector>      filtered_modeled_subgrid_stress_;
  //! the filtered velocities times rho (3 vectors)
  RCP<Epetra_MultiVector>      filtered_dens_vel_;
  //! the filtered density (vector)
  RCP<Epetra_Vector>           filtered_dens_;
  //! the filtered strainrate times rho (vector)
  RCP<Epetra_Vector>           filtered_dens_strainrate_;
  //! the modeled fine scale velocities in nodes (3 vectors)
  RCP<Epetra_MultiVector>      fs_vel_;
  //! the filtered density times temperature times velocity (scalar)
  RCP<Epetra_MultiVector>      filtered_dens_vel_temp_;
  //! the filtered density times temperature gradient times rate of strain (scalar)
  RCP<Epetra_MultiVector>      filtered_dens_rateofstrain_temp_;
//  //! the filtered temperature gradient (scalar)
//  RCP<Epetra_MultiVector>      filtered_gradtemp_;
  //! the filtered temperature (scalar)
  RCP<Epetra_Vector>           filtered_temp_;
  //! the filtered density times temperature (scalar)
  RCP<Epetra_Vector>           filtered_dens_temp_;
  //! the filtered vel exported to column map
  RCP<Epetra_MultiVector>      col_filtered_vel_;
  //! the filtered reystress exported to column map
  RCP<Epetra_MultiVector>      col_filtered_reynoldsstress_;
  //! the modeled subgrid stresses exported to column map
  RCP<Epetra_MultiVector>      col_filtered_modeled_subgrid_stress_;
  //! the filtered velocities times rho exported to column map
  RCP<Epetra_MultiVector>      col_filtered_dens_vel_;
  //! the filtered density exported to column map
  RCP<Epetra_Vector>           col_filtered_dens_;
  //! the filtered strainrate times rho exported to column map
  RCP<Epetra_Vector>           col_filtered_dens_strainrate_;
  //! the modeled fine scale velocities exported to column map
  RCP<Epetra_MultiVector>      col_fs_vel_;
  //! the filtered density times temperature times velocity exported to column map (scalar)
  RCP<Epetra_MultiVector>      col_filtered_dens_vel_temp_;
  //! the filtered density times temperature gradient times rate of strain exported to column map (scalar)
  RCP<Epetra_MultiVector>      col_filtered_dens_rateofstrain_temp_;
//  //! the filtered temperature gradient exported to column map (scalar)
//  RCP<Epetra_MultiVector>      col_filtered_gradtemp_;
  //! the filtered temperature exported to column map (scalar)
  RCP<Epetra_Vector>           col_filtered_temp_;
  //! the filtered density times temperature exported to column map (scalar)
  RCP<Epetra_Vector>           col_filtered_dens_temp_;
  //@}

  //! @name homogeneous flow specials
  //        -------------------------------

  //! the direction coordinates for the above mentioned averaging procedure
  RCP<std::vector<double> >         dir1coords_;
  RCP<std::vector<double> >         dir2coords_;
  //@}

}; // end class DynSmagFilter

} // end namespace FLD

#endif // not DYN_SMAG_H

