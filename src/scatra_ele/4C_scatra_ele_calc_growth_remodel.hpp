// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_SCATRA_ELE_CALC_GROWTH_REMODEL_HPP
#define FOUR_C_SCATRA_ELE_CALC_GROWTH_REMODEL_HPP

#include "4C_config.hpp"

#include "4C_mat_mixture.hpp"
#include "4C_mat_scatra_growth_remodel.hpp"
#include "4C_mixture_constituent_remodelfiber_ssi.hpp"
#include "4C_scatra_ele_calc.hpp"
FOUR_C_NAMESPACE_OPEN

namespace Discret
{
  namespace Elements
  {
    // forward declaration

    template <Core::FE::CellType distype, int probdim = Core::FE::dim<distype>>

    class ScaTraEleCalcGrowthRemodel : public ScaTraEleCalc<distype, probdim>
    {
     private:
      /// (private) protected constructor, since we are a Singleton.
      ScaTraEleCalcGrowthRemodel(
          const int numdofpernode, const int numscal, const std::string& disname);

      /// base class alias
      using my = ScaTraEleCalc<distype, probdim>;

      /// cached constituent pointer of the structural material
      std::vector<const Mixture::MixtureConstituentRemodelFiberSsi*> constituent_;

     public:
      /// Singleton access method
      static ScaTraEleCalcGrowthRemodel<distype, probdim>* instance(
          const int numdofpernode, const int numscal, const std::string& disname);

      int setup_calc(
          Core::Elements::Element* ele, Core::FE::Discretization& discretization) override;

     private:
      using my::nen_;
      using my::nsd_;

      void get_material_params(const Core::Elements::Element* ele, std::vector<double>& densn,
          std::vector<double>& densnp, std::vector<double>& densam, double& visc,
          const int iquad = -1) override;


      std::shared_ptr<Mat::Mixture> get_struct_material(const Core::Elements::Element* ele);

      void materials(const std::shared_ptr<const Core::Mat::Material> material, const int k,
          double& densn, double& densnp, double& densam, double& visc,
          const int iquad = -1) override;

      void mat_scatra_gr(const std::shared_ptr<const Core::Mat::Material> material, const int k,
          double& densn, double& densnp, double& densam, double& visc, const int iquad = -1);

      double get_stress_dependent_reaction_coeff(
          const Mixture::MixtureConstituentRemodelFiberSsi& constituent,
          const Mat::PAR::ScatraGrowthRemodelMat::ScalarQuantity scalar_type, const int gp);

      std::shared_ptr<ScaTraEleReaManager> rea_manager()
      {
        return std::static_pointer_cast<ScaTraEleReaManager>(my::reamanager_);
      }

      std::shared_ptr<ScaTraEleDiffManager> diff_manager()
      {
        return std::static_pointer_cast<ScaTraEleDiffManager>(my::diffmanager_);
      }
    };

  }  // namespace Elements

}  // namespace Discret


FOUR_C_NAMESPACE_CLOSE

#endif  // FOUR_C_SCATRA_ELE_CALC_GROWTH_REMODEL_HPP
