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
#include "../linalg/linalg_utils.H"


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
    Teuchos::RCP<Epetra_Vector>             solution    ,
    const Teuchos::RCP<const Epetra_Vector> dirichtoggle
    );

  void ApplyFilter(
    Teuchos::RCP<Epetra_Vector>             solution    ,
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
          ((*velocity)(i))->ReplaceMyValues(1,&val,&nid);
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
             ((*reystre) (ij))->ReplaceMyValues(1,&val,&nid);
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
          ((*velocity)(i))->ReplaceMyValues(1,&val,&nid);
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
//  void DynSmagBoxFilter(
//    Teuchos::RCP<Epetra_Vector>             velocity    ,
//    const Teuchos::RCP<const Epetra_Vector> dirichtoggle
//    );

  void ApplyBoxFilter(
    const Teuchos::RCP<const Epetra_Vector> velocity,
    const Teuchos::RCP<const Epetra_Vector> dirichtoggle
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
  bool                         homdir_;
  string                       special_flow_homdir_;
  bool                         apply_dynamic_smagorinsky_;
  bool                         apply_box_filter_;
  //@}


  //! @name vectors used for filtering (for dynamic Smagorinsky model)
  //        --------------------------

  //! the box filtered velocities in nodes (3 vectors)
  RCP<Epetra_MultiVector>      filtered_vel_;
  //! the box filtered reynoldsstresses in nodes (9 vectors)
  RCP<Epetra_MultiVector>      filtered_reynoldsstress_;
  //! the modeled subgrid stress in nodes (9 vectors)
  RCP<Epetra_MultiVector>      filtered_modeled_subgrid_stress_;
  //! the modeled fine scale velocities in nodes (3 vectors)
  RCP<Epetra_MultiVector>      fs_vel_;
  //! the filtered vel exported to column map
  RCP<Epetra_MultiVector>      col_filtered_vel_;
  //! the filtered reystress exported to column map
  RCP<Epetra_MultiVector>      col_filtered_reynoldsstress_;
  //! the modeled subgrid stresses exported to column map
  RCP<Epetra_MultiVector>      col_filtered_modeled_subgrid_stress_;
  //! the modeled fine scale velocities exported to column map
  RCP<Epetra_MultiVector>      col_fs_vel_;
  //@}

  //! @name homogeneous flow specials
  //        -------------------------------

  //! the direction coordinates for the above mentioned averaging procedure
  RCP<vector<double> >         dir1coords_;
  RCP<vector<double> >         dir2coords_;
  //@}

}; // end class DynSmagFilter

} // end namespace FLD

#endif // not DYN_SMAG_H

#endif
