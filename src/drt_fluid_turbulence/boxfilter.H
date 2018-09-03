/*!

\file boxfilter.H

\brief Filter methods for the dynamic Smagorinsky, dynamic Vreman and any other model requiring a
box filter

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



#ifndef BOXFILTER_H
#define BOXFILTER_H


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
  class Boxfilter
  {
   public:
    /*!
    \brief Standard Constructor (public)

    */
    Boxfilter(Teuchos::RCP<DRT::Discretization> actdis, Teuchos::ParameterList& params);

    /*!
    \brief Destructor

    */
    virtual ~Boxfilter();

    void AddScatra(Teuchos::RCP<DRT::Discretization> scatradis);

    void InitializeVremanScatra(Teuchos::RCP<DRT::Discretization> scatradis);

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

    void InitializeVreman();

    void ApplyFilter(const Teuchos::RCP<const Epetra_Vector> velocity,
        const Teuchos::RCP<const Epetra_Vector> scalar, const double thermpress,
        const Teuchos::RCP<const Epetra_Vector> dirichtoggle);

    void ApplyFilterScatra(const Teuchos::RCP<const Epetra_Vector> scalar, const double thermpress,
        const Teuchos::RCP<const Epetra_Vector> dirichtoggle, const int ndsvel);


    /*!
    \brief Output of averaged velocity vector for paraview IO

    \param outvec  (in/out) vector in dofrowmap-format to use for
                            output of averaged solution


    */
    /*
    void OutputofAveragedVel(
      Teuchos::RCP<Epetra_Vector> outvec
      )
      {
        // loop all elements on this proc (including ghosted ones)
        for (int n=0;n<discret_->NumMyRowNodes();++n)
        {
          // get the node
          DRT::Node* node = discret_->lRowNode(n);

          std::vector<int> dofs= discret_->Dof(node);

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

          std::vector<int> dofs= discret_->Dof(node);

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
  */
    /*!
    \brief Output of averaged reynoldsstress vector for paraview IO (one component)

    \param outvec  (in/out) vector in dofrowmap-format to use for
                            output of averaged solution
    \param i       which component?
    \param j       which component?


    */
    /*
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

  */
    void GetFilteredVelocity(Teuchos::RCP<Epetra_MultiVector> velocity)
    {
      // loop all elements on this proc (including ghosted ones)
      for (int nid = 0; nid < discret_->NumMyColNodes(); ++nid)
      {
        for (int i = 0; i < 3; ++i)
        {
          double val = (*((*col_filtered_vel_)(i)))[nid];
          int err = ((*velocity)(i))->ReplaceMyValues(1, &val, &nid);
          if (err != 0) dserror("dof not on proc");
        }
      }
      return;
    }



    void GetFilteredDensVelocity(Teuchos::RCP<Epetra_MultiVector> densvelocity)
    {
      // loop all elements on this proc (including ghosted ones)
      for (int nid = 0; nid < discret_->NumMyColNodes(); ++nid)
      {
        for (int i = 0; i < 3; ++i)
        {
          double val = (*((*col_filtered_dens_vel_)(i)))[nid];
          int err = ((*densvelocity)(i))->ReplaceMyValues(1, &val, &nid);
          if (err != 0) dserror("dof not on proc");
        }
      }
      return;
    }

    void GetFilteredDensVelocityTemp(Teuchos::RCP<Epetra_MultiVector> densvelocity)
    {
      // loop all elements on this proc (including ghosted ones)
      for (int nid = 0; nid < discret_->NumMyColNodes(); ++nid)
      {
        for (int i = 0; i < 3; ++i)
        {
          double val = (*((*col_filtered_dens_vel_temp_)(i)))[nid];
          int err = ((*densvelocity)(i))->ReplaceMyValues(1, &val, &nid);
          if (err != 0) dserror("dof not on proc");
        }
      }
      return;
    }

    void GetFilteredDensRateofstrainTemp(Teuchos::RCP<Epetra_MultiVector> densvelocity)
    {
      // loop all elements on this proc (including ghosted ones)
      for (int nid = 0; nid < discret_->NumMyColNodes(); ++nid)
      {
        for (int i = 0; i < 3; ++i)
        {
          double val = (*((*col_filtered_dens_rateofstrain_temp_)(i)))[nid];
          int err = ((*densvelocity)(i))->ReplaceMyValues(1, &val, &nid);
          if (err != 0) dserror("dof not on proc");
        }
      }
      return;
    }

    void GetFilteredReynoldsStress(Teuchos::RCP<Epetra_MultiVector> reystr)
    {
      // loop all elements on this proc (including ghosted ones)
      for (int nid = 0; nid < discret_->NumMyColNodes(); ++nid)
      {
        for (int i = 0; i < 3; ++i)
        {
          for (int j = 0; j < 3; ++j)
          {
            const int ij = 3 * i + j;
            double val = (*((*col_filtered_reynoldsstress_)(ij)))[nid];
            int err = ((*reystr)(ij))->ReplaceMyValues(1, &val, &nid);
            if (err != 0) dserror("dof not on proc");
          }
        }
      }
      return;
    }

    void GetFilteredModeledSubgridStress(Teuchos::RCP<Epetra_MultiVector> mosubstr)
    {
      // loop all elements on this proc (including ghosted ones)
      for (int nid = 0; nid < discret_->NumMyColNodes(); ++nid)
      {
        for (int i = 0; i < 3; ++i)
        {
          for (int j = 0; j < 3; ++j)
          {
            const int ij = 3 * i + j;
            double val = (*((*col_filtered_modeled_subgrid_stress_)(ij)))[nid];
            int err = ((*mosubstr)(ij))->ReplaceMyValues(1, &val, &nid);
            if (err != 0) dserror("dof not on proc");
          }
        }
      }
      return;
    }



    void GetFilteredVremanStrainrate(Teuchos::RCP<Epetra_MultiVector> strainr)
    {
      // loop all elements on this proc (including ghosted ones)
      for (int nid = 0; nid < discret_->NumMyColNodes(); ++nid)
      {
        for (int i = 0; i < 3; ++i)
        {
          for (int j = 0; j < 3; ++j)
          {
            const int ij = 3 * i + j;
            double val = (*((*col_filtered_strainrate_)(ij)))[nid];
            int err = ((*strainr)(ij))->ReplaceMyValues(1, &val, &nid);
            if (err != 0) dserror("dof not on proc");
          }
        }
      }
      return;
    }

    void GetFilteredVremanAlphaij(Teuchos::RCP<Epetra_MultiVector> alphij)
    {
      // loop all elements on this proc (including ghosted ones)
      for (int nid = 0; nid < discret_->NumMyColNodes(); ++nid)
      {
        for (int i = 0; i < 3; ++i)
        {
          for (int j = 0; j < 3; ++j)
          {
            const int ij = 3 * i + j;
            double val = (*((*col_filtered_alphaij_)(ij)))[nid];
            int err = ((*alphij)(ij))->ReplaceMyValues(1, &val, &nid);
            if (err != 0) dserror("dof not on proc");
          }
        }
      }
      return;
    }
    void GetFilteredVremanAlphaijsc(Teuchos::RCP<Epetra_MultiVector> alphijsc)
    {
      // loop all elements on this proc (including ghosted ones)
      for (int nid = 0; nid < discret_->NumMyColNodes(); ++nid)
      {
        for (int i = 0; i < 3; ++i)
        {
          for (int j = 0; j < 3; ++j)
          {
            const int ij = 3 * i + j;
            double val = (*((*col_filtered_alphaijsc_)(ij)))[nid];
            int err = ((*alphijsc)(ij))->ReplaceMyValues(1, &val, &nid);
            if (err != 0) dserror("dof not on proc");
          }
        }
      }
      return;
    }


    void GetExpression(Teuchos::RCP<Epetra_Vector> expr)
    {
      // loop all elements on this proc (including ghosted ones)
      for (int nid = 0; nid < discret_->NumMyColNodes(); ++nid)
      {
        for (int i = 0; i < 1; ++i)
        {
          double val = (*((*col_filtered_expression_)(i)))[nid];
          int err = ((*expr)(i))->ReplaceMyValues(1, &val, &nid);
          if (err != 0) dserror("dof not on proc");
        }
      }
      return;
    }


    void GetDensity(Teuchos::RCP<Epetra_Vector> dens)
    {
      // loop all elements on this proc (including ghosted ones)
      for (int nid = 0; nid < discret_->NumMyColNodes(); ++nid)
      {
        for (int i = 0; i < 1; ++i)
        {
          double val = (*((*col_filtered_dens_)(i)))[nid];
          int err = ((*dens)(i))->ReplaceMyValues(1, &val, &nid);
          if (err != 0) dserror("dof not on proc");
        }
      }
      return;
    }

    void GetTemp(Teuchos::RCP<Epetra_Vector> dens)
    {
      // loop all elements on this proc (including ghosted ones)
      for (int nid = 0; nid < discret_->NumMyColNodes(); ++nid)
      {
        for (int i = 0; i < 1; ++i)
        {
          double val = (*((*col_filtered_temp_)(i)))[nid];
          int err = ((*dens)(i))->ReplaceMyValues(1, &val, &nid);
          if (err != 0) dserror("dof not on proc");
        }
      }
      return;
    }

    void GetDensTemp(Teuchos::RCP<Epetra_Vector> dens)
    {
      // loop all elements on this proc (including ghosted ones)
      for (int nid = 0; nid < discret_->NumMyColNodes(); ++nid)
      {
        for (int i = 0; i < 1; ++i)
        {
          double val = (*((*col_filtered_dens_temp_)(i)))[nid];
          int err = ((*dens)(i))->ReplaceMyValues(1, &val, &nid);
          if (err != 0) dserror("dof not on proc");
        }
      }
      return;
    }


    void GetDensityStrainrate(Teuchos::RCP<Epetra_Vector> densstrainr)
    {
      // loop all elements on this proc (including ghosted ones)
      for (int nid = 0; nid < discret_->NumMyColNodes(); ++nid)
      {
        for (int i = 0; i < 1; ++i)
        {
          double val = (*((*col_filtered_dens_strainrate_)(i)))[nid];
          int err = ((*densstrainr)(i))->ReplaceMyValues(1, &val, &nid);
          if (err != 0) dserror("dof not on proc");
        }
      }
      return;
    }

    void GetAlpha2(Teuchos::RCP<Epetra_Vector> alph2)
    {
      // loop all elements on this proc (including ghosted ones)
      for (int nid = 0; nid < discret_->NumMyColNodes(); ++nid)
      {
        for (int i = 0; i < 1; ++i)
        {
          double val = (*((*col_filtered_alpha2_)(i)))[nid];
          int err = ((*alph2)(i))->ReplaceMyValues(1, &val, &nid);
          if (err != 0) dserror("dof not on proc");
        }
      }
      return;
    }



    void GetFineScaleVelocity(Teuchos::RCP<Epetra_MultiVector> velocity)
    {
      // loop all elements on this proc (including ghosted ones)
      for (int nid = 0; nid < discret_->NumMyColNodes(); ++nid)
      {
        for (int i = 0; i < 3; ++i)
        {
          double val = (*((*col_fs_vel_)(i)))[nid];
          int err = ((*velocity)(i))->ReplaceMyValues(1, &val, &nid);
          if (err != 0) dserror("dof not on proc");
        }
      }
      return;
    }


    void FilteredReyStrComp(Teuchos::RCP<Epetra_Vector> outvec, const int i, const int j)
    {
      // loop all elements on this proc (including ghosted ones)
      int ij = 3 * (i - 1) + (j - 1);
      for (int nid = 0; nid < discret_->NumMyRowNodes(); ++nid)
      {
        // get the node
        DRT::Node* node = discret_->lRowNode(nid);
        int id = node->Id();

        double val = (*((*filtered_reynoldsstress_)(ij)))[nid];

        (outvec)->ReplaceGlobalValues(1, &val, &id);
      }
      return;
    }

    void FilteredVelComp(Teuchos::RCP<Epetra_Vector> outvec, const int i, const int j)
    {
      // loop all elements on this proc (including ghosted ones)
      for (int nid = 0; nid < discret_->NumMyRowNodes(); ++nid)
      {
        // get the node
        DRT::Node* node = discret_->lRowNode(nid);
        int id = node->Id();

        double vali = (*((*filtered_vel_)(i - 1)))[nid];
        double valj = (*((*filtered_vel_)(j - 1)))[nid];
        double val = vali * valj;

        (outvec)->ReplaceGlobalValues(1, &val, &id);
      }
      return;
    }

    void OutputofFineScaleVel(Teuchos::RCP<Epetra_Vector> outvec)
    {
      // loop all elements on this proc (including ghosted ones)
      for (int n = 0; n < discret_->NumMyRowNodes(); ++n)
      {
        // get the node
        DRT::Node* node = discret_->lRowNode(n);

        std::vector<int> dofs = discret_->Dof(node);

        for (int d = 0; d < discret_->NumDof(node) - 1; ++d)
        {
          double val = (*((*fs_vel_)(d)))[n];

          int id = dofs[d];

          (outvec)->ReplaceGlobalValues(1, &val, &id);
        }
      }
      return;
    }

    void OutputofAveragedVel(Teuchos::RCP<Epetra_Vector> outvec)
    {
      // loop all elements on this proc (including ghosted ones)
      for (int n = 0; n < discret_->NumMyRowNodes(); ++n)
      {
        // get the node
        DRT::Node* node = discret_->lRowNode(n);

        std::vector<int> dofs = discret_->Dof(node);

        for (int d = 0; d < discret_->NumDof(node) - 1; ++d)
        {
          double val = (*((*filtered_vel_)(d)))[n];

          int id = dofs[d];

          (outvec)->ReplaceGlobalValues(1, &val, &id);
        }
      }
      return;
    }

    void GetFilteredPhi(Teuchos::RCP<Epetra_MultiVector> phi)
    {
      // loop all elements on this proc (including ghosted ones)
      for (int nid = 0; nid < discret_->NumMyColNodes(); ++nid)
      {
        for (int i = 0; i < 3; ++i)
        {
          double val = (*((*col_filtered_phi_)(i)))[nid];
          int err = ((*phi)(i))->ReplaceMyValues(1, &val, &nid);
          if (err != 0) dserror("dof not on proc");
        }
      }
      return;
    }

    void GetFilteredPhi2(Teuchos::RCP<Epetra_Vector> phi2)
    {
      // loop all elements on this proc (including ghosted ones)
      for (int nid = 0; nid < discret_->NumMyColNodes(); ++nid)
      {
        for (int i = 0; i < 1; ++i)
        {
          double val = (*((*col_filtered_phi2_)(i)))[nid];
          int err = ((*phi2)(i))->ReplaceMyValues(1, &val, &nid);
          if (err != 0) dserror("dof not on proc");
        }
      }
      return;
    }

    void GetFilteredPhiexpression(Teuchos::RCP<Epetra_Vector> phiexpression)
    {
      // loop all elements on this proc (including ghosted ones)
      for (int nid = 0; nid < discret_->NumMyColNodes(); ++nid)
      {
        for (int i = 0; i < 1; ++i)
        {
          double val = (*((*col_filtered_phiexpression_)(i)))[nid];
          int err = ((*phiexpression)(i))->ReplaceMyValues(1, &val, &nid);
          if (err != 0) dserror("dof not on proc");
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



    void ApplyBoxFilter(const Teuchos::RCP<const Epetra_Vector> velocity,
        const Teuchos::RCP<const Epetra_Vector> scalar, const double thermpress,
        const Teuchos::RCP<const Epetra_Vector> dirichtoggle);

    void ApplyBoxFilterScatra(const Teuchos::RCP<const Epetra_Vector> scalar,
        const double thermpress, const Teuchos::RCP<const Epetra_Vector> dirichtoggle,
        const int ndsvel);



    //! @name input arguments of the constructor
    //
    //! the discretization
    Teuchos::RCP<DRT::Discretization> discret_;
    //! parameterlist including time params, stabilization params and turbulence sublist
    Teuchos::ParameterList& params_;
    //! flag for physical type of fluid flow
    INPAR::FLUID::PhysicalType physicaltype_;
    //@}

    //! @name control parameters
    bool apply_dynamic_smagorinsky_;
    bool vreman_dynamic_;
    bool apply_box_filter_;
    bool dynsmag_loma_on_;
    bool loma_;
    bool incomp_;
    // bool for field values
    bool velocity_;
    bool reynoldsstress_;
    bool modeled_subgrid_stress_;
    bool expression_;
    bool strainrate_;
    bool alphaij_;
    bool alpha2_;
    bool finescale_velocity_;
    bool densvelocity_;
    bool densstrainrate_;
    bool density_;
    bool phi_;
    bool phi2_;
    bool phiexpression_;
    bool alphaijsc_;

    //@}

    //! @name special scatra variables
    //! the discretization
    Teuchos::RCP<DRT::Discretization> scatradiscret_;
    //@}

    //! @name vectors used for filtering (for dynamic Smagorinsky model)
    //        --------------------------

    //! the box filtered velocities in nodes (3 vectors)
    Teuchos::RCP<Epetra_MultiVector> filtered_vel_;
    //! the box filtered reynoldsstresses in nodes (9 vectors)
    Teuchos::RCP<Epetra_MultiVector> filtered_reynoldsstress_;
    //! the modeled subgrid stress in nodes (9 vectors)
    Teuchos::RCP<Epetra_MultiVector> filtered_modeled_subgrid_stress_;
    //! the filtered velocities times rho (3 vectors)
    Teuchos::RCP<Epetra_MultiVector> filtered_dens_vel_;
    //! the filtered density (vector)
    Teuchos::RCP<Epetra_Vector> filtered_dens_;
    //! the filtered strainrate times rho (vector)
    Teuchos::RCP<Epetra_Vector> filtered_dens_strainrate_;
    //! the modeled fine scale velocities in nodes (3 vectors)
    Teuchos::RCP<Epetra_MultiVector> fs_vel_;
    //! Vreman model: strain rate in nodes (9 vectors)
    Teuchos::RCP<Epetra_MultiVector> filtered_strainrate_;
    //! Vreman model: expression (Vector)
    Teuchos::RCP<Epetra_Vector> filtered_expression_;
    //! Vreman model: alphaij in nodes (9 vectors)
    Teuchos::RCP<Epetra_MultiVector> filtered_alphaij_;
    //! Vreman model: alpha2 (Vector)
    Teuchos::RCP<Epetra_Vector> filtered_alpha2_;
    //! the filtered density times temperature times velocity (scalar)
    Teuchos::RCP<Epetra_MultiVector> filtered_dens_vel_temp_;
    //! the filtered density times temperature gradient times rate of strain (scalar)
    Teuchos::RCP<Epetra_MultiVector> filtered_dens_rateofstrain_temp_;
    //  //! the filtered temperature gradient (scalar)
    //  Teuchos::RCP<Epetra_MultiVector>      filtered_gradtemp_;
    //! the filtered temperature (scalar)
    Teuchos::RCP<Epetra_Vector> filtered_temp_;
    //! the filtered density times temperature (scalar)
    Teuchos::RCP<Epetra_Vector> filtered_dens_temp_;
    //! filtered phi: gradient
    Teuchos::RCP<Epetra_MultiVector> filtered_phi_;
    //! filtered grad(phi_j)*grad(phi_j)
    Teuchos::RCP<Epetra_Vector> filtered_phi2_;
    //! filtered visc*grad(phi_j)*grad(phi_j)
    Teuchos::RCP<Epetra_Vector> filtered_phiexpression_;
    //! alphaij tensor for scatra
    Teuchos::RCP<Epetra_MultiVector> filtered_alphaijsc_;
    //! the filtered vel exported to column map
    Teuchos::RCP<Epetra_MultiVector> col_filtered_vel_;
    //! the filtered reystress exported to column map
    Teuchos::RCP<Epetra_MultiVector> col_filtered_reynoldsstress_;
    //! the modeled subgrid stresses exported to column map
    Teuchos::RCP<Epetra_MultiVector> col_filtered_modeled_subgrid_stress_;
    //! the filtered velocities times rho exported to column map
    Teuchos::RCP<Epetra_MultiVector> col_filtered_dens_vel_;
    //! the filtered density exported to column map
    Teuchos::RCP<Epetra_Vector> col_filtered_dens_;
    //! the filtered strainrate times rho exported to column map
    Teuchos::RCP<Epetra_Vector> col_filtered_dens_strainrate_;
    //! the modeled fine scale velocities exported to column map
    Teuchos::RCP<Epetra_MultiVector> col_fs_vel_;
    //! Vreman model: strain rate exported to column map
    Teuchos::RCP<Epetra_MultiVector> col_filtered_strainrate_;
    //! Vreman model: expression (Vector) exported to column map
    Teuchos::RCP<Epetra_Vector> col_filtered_expression_;
    //! Vreman model: alphaij in nodes (9 vectors) exported to column map
    Teuchos::RCP<Epetra_MultiVector> col_filtered_alphaij_;
    //! Vreman model: alpha2 (Vector) exported to column map
    Teuchos::RCP<Epetra_Vector> col_filtered_alpha2_;
    //! the filtered density times temperature times velocity exported to column map (scalar)
    Teuchos::RCP<Epetra_MultiVector> col_filtered_dens_vel_temp_;
    //! the filtered density times temperature gradient times rate of strain exported to column map
    //! (scalar)
    Teuchos::RCP<Epetra_MultiVector> col_filtered_dens_rateofstrain_temp_;
    //  //! the filtered temperature gradient exported to column map (scalar)
    //  Teuchos::RCP<Epetra_MultiVector>      col_filtered_gradtemp_;
    //! the filtered temperature exported to column map (scalar)
    Teuchos::RCP<Epetra_Vector> col_filtered_temp_;
    //! the filtered density times temperature exported to column map (scalar)
    Teuchos::RCP<Epetra_Vector> col_filtered_dens_temp_;
    //! filtered phi: gradient
    Teuchos::RCP<Epetra_MultiVector> col_filtered_phi_;
    //! filtered grad(phi_j)*grad(phi_j)
    Teuchos::RCP<Epetra_Vector> col_filtered_phi2_;
    //! filtered visc*grad(phi_j)*grad(phi_j)
    Teuchos::RCP<Epetra_Vector> col_filtered_phiexpression_;
    //! alphaij tensor for scatra
    Teuchos::RCP<Epetra_MultiVector> col_filtered_alphaijsc_;
    //@}

    //! @name homogeneous flow specials
    //        -------------------------------

    //! the direction coordinates for the above mentioned averaging procedure
    Teuchos::RCP<std::vector<double>> dir1coords_;
    Teuchos::RCP<std::vector<double>> dir2coords_;
    //@}

  };  // end class Boxfilter

}  // end namespace FLD

#endif  // not BOXFILTER_H
