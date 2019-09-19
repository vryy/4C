/*----------------------------------------------------------------------*/
/*! \file

\brief manages a mesh coupling object with knowledge of a level set field

\level 3

\maintainer Christoph Ager
            ager@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
*/
/*----------------------------------------------------------------------*/


#ifndef XFEM_COUPLING_MESH_COUPLED_LEVELSET_H_
#define XFEM_COUPLING_MESH_COUPLED_LEVELSET_H_

#include "xfem_coupling_base.H"
#include "xfem_coupling_levelset.H"
#include "xfem_interface_utils.H"

namespace XFEM
{
  /*!
  \brief
   */
  class MeshCouplingNavierSlipTwoPhase : public MeshCouplingNavierSlip
  {
   public:
    //! constructor
    explicit MeshCouplingNavierSlipTwoPhase(
        Teuchos::RCP<DRT::Discretization>& bg_dis,  ///< background discretization
        const std::string& cond_name,  ///< name of the condition, by which the derived cutter
                                       ///< discretization is identified
        Teuchos::RCP<DRT::Discretization>&
            cond_dis,           ///< discretization from which cutter discretization can be derived
        const int coupling_id,  ///< id of composite of coupling conditions
        const double time,      ///< time
        const int step,         ///< time step
        bool marked_geometry = false  ///< is this a marked geometry mesh boundary
    );

    //! destructor
    virtual ~MeshCouplingNavierSlipTwoPhase() { return; }

    void SetTwoPhaseCouplingPointer(
        Teuchos::RCP<XFEM::LevelSetCouplingTwoPhase> twophase_coupling_object);

    /*!
    Return prescribed velocities and traction vectors for a GNBC boundary condition.
    Also returns the projection matrix (to the plane of the surface) needed for the GNBC condition.
    */
    template <DRT::Element::DiscretizationType DISTYPE, class V1, class V2, class X1, class T1,
        class M1, class M2, class M3>
    void EvaluateCouplingConditions(V1& ivel,  ///< prescribed velocity at interface
        V2& itraction,                         ///< prescribed traction at interface
        X1& x,                                 ///< coordinates of gauss point
        const DRT::Condition* cond,            ///< condition prescribed to this surface
        T1& proj_matrix,  ///< Laplace-Beltrami matrix for surface tension calculations
        int eid,          ///< element ID
        M1& funct,        ///< local shape function for Gauss Point (from fluid element)
        M2& derxy,   ///< local derivatives of shape function for Gauss Point (from fluid element)
        M3& normal,  ///< surface normal of cut element
        const bool& eval_dirich_at_gp,
        double& kappa_m,  ///< fluid sided weighting
        double& visc_m,   ///< fluid sided weighting
        double& visc_s    ///< slave sided dynamic viscosity
    )
    {
      SetupProjectionMatrix(proj_matrix, normal);

      // help variable
      int robin_id_dirch;

      if (eval_dirich_at_gp)
      {
        // evaluate interface velocity (given by weak Dirichlet condition)
        robin_id_dirch = cond->GetInt("robin_id_dirch");
        // Check if int is negative (signbit(x) -> x<0 true, x=>0 false)
        if (!std::signbit(robin_id_dirch))
          EvaluateDirichletFunction(
              ivel, x, conditionsmap_robin_dirch_.find(robin_id_dirch)->second, time_);

          // Safety checks
#ifdef DEBUG
        if ((conditionsmap_robin_dirch_.find(robin_id_dirch)) == conditionsmap_robin_dirch_.end())
        {
          dserror("Key was not found in this instance!! Fatal error! (conditionsmap_robin_dirch_)");
        }
#endif
      }

      // evaluate interface traction (given by Neumann condition)
      robin_id_dirch = cond->GetInt("robin_id_neumann");
      if (!std::signbit(robin_id_dirch))
      {
        // This is maybe not the most efficient implementation as we evaluate dynvisc as well as the
        // sliplenght twice (also done in UpdateConfigurationMap_GP ... as soon as this gets
        // relevant we should merge this functions)

        // evaluate interface traction (given by Neumann condition)
        // Add this to the veljump!

        double sliplength = 0.0;
        GetSliplengthAtGP<DISTYPE>(sliplength, eid, cond, funct, derxy);

        if (sliplength < 0.0) dserror("The slip length can not be negative.");

        if (sliplength != 0.0)
        {
          EvaluateNeumannFunction(
              itraction, x, conditionsmap_robin_neumann_.find(robin_id_dirch)->second, time_);

          double sl_visc_fac = sliplength / (kappa_m * visc_m + (1.0 - kappa_m) * visc_s);
          LINALG::Matrix<3, 1> tmp_itraction(true);
          tmp_itraction.MultiplyTN(proj_matrix, itraction);
          // Project this into tangential direction!!!

          ivel.Update(sl_visc_fac, tmp_itraction, 1.0);

          itraction.Clear();
        }
      }

      if (force_tangvel_map_.find(cond->Id())->second)
      {
        LINALG::Matrix<3, 1> tmp_ivel(true);
        tmp_ivel.MultiplyTN(
            proj_matrix, ivel);  // apply Projection matrix from the right. (u_0 * P^t)
        ivel.Update(1.0, tmp_ivel, 0.0);
      }

// Safety checks
#ifdef DEBUG
      if (!std::signbit(robin_id_dirch))
      {
        if ((conditionsmap_robin_neumann_.find(robin_id_dirch)) ==
            conditionsmap_robin_neumann_.end())
        {
          dserror(
              "Key was not found in this instance!! Fatal error! (conditionsmap_robin_neumann_)");
        }
      }
      std::map<int, bool>::iterator it_bool;
      if ((it_bool = force_tangvel_map_.find(cond->Id())) == force_tangvel_map_.end())
      {
        dserror("Key was not found in this instance!! Fatal error! (force_tangvel_map_)");
      }
#endif
    };


    // template <DRT::Element::DiscretizationType DISTYPE>//,class M1, class M2>
    //! Updates configurationmap for specific Gausspoint
    virtual void UpdateConfigurationMap_GP(double& kappa_m,  //< fluid sided weighting
        double& visc_m,                                      //< master sided dynamic viscosity
        double& visc_s,                                      //< slave sided dynamic viscosity
        double& density_m,                                   //< master sided density
        double& visc_stab_tang,         //< viscous tangential NIT Penalty scaling
        double& full_stab,              //< full NIT Penalty scaling
        const LINALG::Matrix<3, 1>& x,  //< Position x in global coordinates
        const DRT::Condition* cond,     //< Condition
        DRT::Element* ele,              //< Element
        DRT::Element* bele,             //< Boundary Element
        double* funct,  //< local shape function for Gauss Point (from fluid element)
        double* derxy,  //< local derivatives of shape function for Gauss Point (from fluid element)
        LINALG::Matrix<3, 1>& rst_slave,  //< local coord of gp on slave boundary element
        LINALG::Matrix<3, 1>& normal,     //< normal at gp
        LINALG::Matrix<3, 1>& vel_m,      //< master velocity at gp
        double* fulltraction              //< precomputed fsi traction (sigmaF n + gamma relvel)
    )
    {
      double dynvisc = (kappa_m * visc_m + (1.0 - kappa_m) * visc_s);
      double sliplength = 0.0;

      int eid = ele->Id();

      if (ele->Shape() == DRT::Element::hex8)
      {
        const DRT::Element::DiscretizationType shape = DRT::Element::hex8;
        //
        const size_t nsd = DRT::UTILS::DisTypeToDim<shape>::dim;
        const size_t nen = DRT::UTILS::DisTypeToNumNodePerEle<shape>::numNodePerElement;
        LINALG::Matrix<nen, 1> funct_(funct, true);
        LINALG::Matrix<nen, nsd> derxy_(derxy, true);
        GetSliplengthAtGP<DRT::Element::hex8>(sliplength, eid, cond, funct_, derxy_);
      }
      else if (ele->Shape() == DRT::Element::hex27)
      {
        const DRT::Element::DiscretizationType shape = DRT::Element::hex27;
        //
        const size_t nsd = DRT::UTILS::DisTypeToDim<shape>::dim;
        const size_t nen = DRT::UTILS::DisTypeToNumNodePerEle<shape>::numNodePerElement;
        LINALG::Matrix<nen, 1> funct_(funct, true);
        LINALG::Matrix<nen, nsd> derxy_(derxy, true);
        GetSliplengthAtGP<DRT::Element::hex27>(sliplength, eid, cond, funct_, derxy_);
      }
      else if (ele->Shape() == DRT::Element::hex20)
      {
        const DRT::Element::DiscretizationType shape = DRT::Element::hex20;
        //
        const size_t nsd = DRT::UTILS::DisTypeToDim<shape>::dim;
        const size_t nen = DRT::UTILS::DisTypeToNumNodePerEle<shape>::numNodePerElement;
        LINALG::Matrix<nen, 1> funct_(funct, true);
        LINALG::Matrix<nen, nsd> derxy_(derxy, true);
        GetSliplengthAtGP<DRT::Element::hex20>(sliplength, eid, cond, funct_, derxy_);
      }
      else
        dserror("Element not considered.");

      if (sliplength < 0.0) dserror("The slip length can not be negative.");

      // Extract correct slip length - bool pair for this condition ID.
      double& scaled_normal_pen = scalednormalpenalty_map_.find(cond->Id())->second;

      if (sliplength != 0.0)
      {
        double stabnit = 0.0;
        double stabadj = 0.0;
        XFEM::UTILS::GetNavierSlipStabilizationParameters(
            visc_stab_tang, dynvisc, sliplength, stabnit, stabadj);
        configuration_map_[INPAR::XFEM::F_Pen_t_Row].second = stabnit;
        configuration_map_[INPAR::XFEM::F_Con_t_Row] =
            std::pair<bool, double>(true, -stabnit);  //+sign for penalty!
        configuration_map_[INPAR::XFEM::F_Con_t_Col] =
            std::pair<bool, double>(true, sliplength / dynvisc);
        configuration_map_[INPAR::XFEM::F_Adj_t_Row].second = stabadj;
        configuration_map_[INPAR::XFEM::FStr_Adj_t_Col] = std::pair<bool, double>(true, sliplength);
      }
      else
      {
        configuration_map_[INPAR::XFEM::F_Pen_t_Row].second = full_stab;
        configuration_map_[INPAR::XFEM::F_Con_t_Row] = std::pair<bool, double>(false, 0.0);
        configuration_map_[INPAR::XFEM::F_Con_t_Col] = std::pair<bool, double>(false, 0.0);
        configuration_map_[INPAR::XFEM::F_Adj_t_Row].second = 1.0;
        configuration_map_[INPAR::XFEM::FStr_Adj_t_Col] = std::pair<bool, double>(false, 0.0);
      }

      // Configuration of Penalty Terms
      configuration_map_[INPAR::XFEM::F_Pen_n_Row].second =
          scaled_normal_pen * full_stab;  // full_stab <-- to keep results!
    };

   protected:
    template <DRT::Element::DiscretizationType DISTYPE, class M1, class M2>
    void GetSliplengthAtGP(double& sliplength,  ///< phi at gausspoint
        int eid,                                ///< element ID
        const DRT::Condition* cond,
        M1& funct,  ///< local shape function for Gauss Point (from fluid element)
        M2& derxy   ///< local derivatives of shape function for Gauss Point (from fluid element)
    )
    {
      double phi_at_gp = 0.0;
      if (twophase_coupling_object_ == Teuchos::null)
        dserror("twophase_coupling_object_, object is null-pointer.");

      twophase_coupling_object_->GetPhiAtGP<DISTYPE>(phi_at_gp, eid, funct, derxy);

      // Extract correct slip length - bool pair for this condition ID.
      double& slipsmear_length = slipsmear_map_.find(cond->Id())->second;
      // Extract correct slip length - bool pair for this condition ID.
      std::pair<double, bool>& tmp_pair = sliplength_map_.find(cond->Id())->second;

      // max smearlength at phi = 0;
      if (abs(phi_at_gp) <= slipsmear_length)
      {
        sliplength = tmp_pair.first;  // 1000.0;
        // double max_smearlength = 100;
        // double sqrt_max_smearlength = sqrt(max_smearlength);
        // sliplength=(0.5*(sqrt_max_smearlength +
        // sqrt_max_smearlength*std::cos(M_PI*phi_at_gp/slipsmear_length)))^2;
      }
    };

    virtual void SetConditionSpecificParameters();

   private:
    //! two-phase coupling object
    Teuchos::RCP<XFEM::LevelSetCouplingTwoPhase> twophase_coupling_object_;

    //       Coupling Surface (E 1)  Length of slip smear
    std::map<int, double> slipsmear_map_;
    std::map<int, double> scalednormalpenalty_map_;
  };

}  // namespace XFEM

#endif /*XFEM_COUPLING_MESH_COUPLED_LEVELSET_H_*/
