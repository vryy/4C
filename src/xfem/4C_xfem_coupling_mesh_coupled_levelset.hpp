/*----------------------------------------------------------------------*/
/*! \file

\brief manages a mesh coupling object with knowledge of a level set field

\level 3

*/
/*----------------------------------------------------------------------*/


#ifndef FOUR_C_XFEM_COUPLING_MESH_COUPLED_LEVELSET_HPP
#define FOUR_C_XFEM_COUPLING_MESH_COUPLED_LEVELSET_HPP

#include "4C_config.hpp"

#include "4C_xfem_coupling_base.hpp"
#include "4C_xfem_coupling_levelset.hpp"
#include "4C_xfem_coupling_mesh.hpp"
#include "4C_xfem_interface_utils.hpp"

FOUR_C_NAMESPACE_OPEN

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


    /*!
    Return prescribed velocities and traction vectors for a GNBC boundary condition.
    Also returns the projection matrix (to the plane of the surface) needed for the GNBC condition.
    */
    template <CORE::FE::CellType DISTYPE, class V1, class V2, class X1, class T1, class M1,
        class M2, class M3>
    void evaluate_coupling_conditions(V1& ivel,   ///< prescribed velocity at interface
        V2& itraction,                            ///< prescribed traction at interface
        X1& x,                                    ///< coordinates of gauss point
        const CORE::Conditions::Condition* cond,  ///< condition prescribed to this surface
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
      setup_projection_matrix(proj_matrix, normal);

      // help variable
      int robin_id_dirch;

      if (eval_dirich_at_gp)
      {
        // evaluate interface velocity (given by weak Dirichlet condition)
        robin_id_dirch = cond->parameters().Get<int>("robin_id_dirch");
        // Check if int is negative (signbit(x) -> x<0 true, x=>0 false)
        if (!std::signbit(static_cast<double>(robin_id_dirch)))
          evaluate_dirichlet_function(
              ivel, x, conditionsmap_robin_dirch_.find(robin_id_dirch)->second, time_);

          // Safety checks
#ifdef FOUR_C_ENABLE_ASSERTIONS
        if ((conditionsmap_robin_dirch_.find(robin_id_dirch)) == conditionsmap_robin_dirch_.end())
        {
          FOUR_C_THROW(
              "Key was not found in this instance!! Fatal error! (conditionsmap_robin_dirch_)");
        }
#endif
      }

      // evaluate interface traction (given by Neumann condition)
      robin_id_dirch = cond->parameters().Get<int>("robin_id_neumann");
      if (!std::signbit(static_cast<double>(robin_id_dirch)))
      {
        // This is maybe not the most efficient implementation as we evaluate dynvisc as well as the
        // sliplenght twice (also done in update_configuration_map_gp ... as soon as this gets
        // relevant we should merge this functions)

        // evaluate interface traction (given by Neumann condition)
        // Add this to the veljump!

        double sliplength = 0.0;

        if (sliplength < 0.0) FOUR_C_THROW("The slip length can not be negative.");

        if (sliplength != 0.0)
        {
          evaluate_neumann_function(
              itraction, x, conditionsmap_robin_neumann_.find(robin_id_dirch)->second, time_);

          double sl_visc_fac = sliplength / (kappa_m * visc_m + (1.0 - kappa_m) * visc_s);
          CORE::LINALG::Matrix<3, 1> tmp_itraction(true);
          tmp_itraction.MultiplyTN(proj_matrix, itraction);
          // Project this into tangential direction!!!

          ivel.Update(sl_visc_fac, tmp_itraction, 1.0);

          itraction.Clear();
        }
      }

      if (force_tangvel_map_.find(cond->Id())->second)
      {
        CORE::LINALG::Matrix<3, 1> tmp_ivel(true);
        tmp_ivel.MultiplyTN(
            proj_matrix, ivel);  // apply Projection matrix from the right. (u_0 * P^t)
        ivel.Update(1.0, tmp_ivel, 0.0);
      }

// Safety checks
#ifdef FOUR_C_ENABLE_ASSERTIONS
      if (!std::signbit(static_cast<double>(robin_id_dirch)))
      {
        if ((conditionsmap_robin_neumann_.find(robin_id_dirch)) ==
            conditionsmap_robin_neumann_.end())
        {
          FOUR_C_THROW(
              "Key was not found in this instance!! Fatal error! (conditionsmap_robin_neumann_)");
        }
      }
      std::map<int, bool>::iterator it_bool;
      if ((it_bool = force_tangvel_map_.find(cond->Id())) == force_tangvel_map_.end())
      {
        FOUR_C_THROW("Key was not found in this instance!! Fatal error! (force_tangvel_map_)");
      }
#endif
    };


    // template <CORE::FE::CellType DISTYPE>//,class M1, class M2>
    //! Updates configurationmap for specific Gausspoint
    void update_configuration_map_gp(double& kappa_m,  //< fluid sided weighting
        double& visc_m,                                //< master sided dynamic viscosity
        double& visc_s,                                //< slave sided dynamic viscosity
        double& density_m,                             //< master sided density
        double& visc_stab_tang,                        //< viscous tangential NIT Penalty scaling
        double& full_stab,                             //< full NIT Penalty scaling
        const CORE::LINALG::Matrix<3, 1>& x,           //< Position x in global coordinates
        const CORE::Conditions::Condition* cond,       //< Condition
        CORE::Elements::Element* ele,                  //< Element
        CORE::Elements::Element* bele,                 //< Boundary Element
        double* funct,  //< local shape function for Gauss Point (from fluid element)
        double* derxy,  //< local derivatives of shape function for Gauss Point (from fluid element)
        CORE::LINALG::Matrix<3, 1>& rst_slave,  //< local coord of gp on slave boundary element
        CORE::LINALG::Matrix<3, 1>& normal,     //< normal at gp
        CORE::LINALG::Matrix<3, 1>& vel_m,      //< master velocity at gp
        double* fulltraction  //< precomputed fsi traction (sigmaF n + gamma relvel)
        ) override
    {
      double dynvisc = (kappa_m * visc_m + (1.0 - kappa_m) * visc_s);
      double sliplength = 0.0;

      if (ele->Shape() == CORE::FE::CellType::hex8)
      {
        const CORE::FE::CellType shape = CORE::FE::CellType::hex8;
        //
        const size_t nsd = CORE::FE::dim<shape>;
        const size_t nen = CORE::FE::num_nodes<shape>;
        CORE::LINALG::Matrix<nen, 1> funct_(funct, true);
        CORE::LINALG::Matrix<nen, nsd> derxy_(derxy, true);
      }
      else if (ele->Shape() == CORE::FE::CellType::hex27)
      {
        const CORE::FE::CellType shape = CORE::FE::CellType::hex27;
        //
        const size_t nsd = CORE::FE::dim<shape>;
        const size_t nen = CORE::FE::num_nodes<shape>;
        CORE::LINALG::Matrix<nen, 1> funct_(funct, true);
        CORE::LINALG::Matrix<nen, nsd> derxy_(derxy, true);
      }
      else if (ele->Shape() == CORE::FE::CellType::hex20)
      {
        const CORE::FE::CellType shape = CORE::FE::CellType::hex20;
        //
        const size_t nsd = CORE::FE::dim<shape>;
        const size_t nen = CORE::FE::num_nodes<shape>;
        CORE::LINALG::Matrix<nen, 1> funct_(funct, true);
        CORE::LINALG::Matrix<nen, nsd> derxy_(derxy, true);
      }
      else
        FOUR_C_THROW("Element not considered.");

      if (sliplength < 0.0) FOUR_C_THROW("The slip length can not be negative.");

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
    };

    void set_condition_specific_parameters() override;
  };

}  // namespace XFEM

FOUR_C_NAMESPACE_CLOSE

#endif
