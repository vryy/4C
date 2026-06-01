// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_MAT_VISCOELAST_SUMMAND_HPP
#define FOUR_C_MAT_VISCOELAST_SUMMAND_HPP

#include "4C_config.hpp"

#include "4C_comm_parobject.hpp"
#include "4C_legacy_enum_definitions_materials.hpp"
#include "4C_linalg_fixedsizematrix.hpp"
#include "4C_solid_ele_fibers.hpp"

#include <memory>
#include <optional>
#include <string>
#include <vector>

namespace Teuchos
{
  class ParameterList;
}

FOUR_C_NAMESPACE_OPEN

namespace Mat::ViscoElast
{
  class Summand : public Core::Communication::ParObject
  {
   public:
    Summand() = default;
    ~Summand() override = default;

    int unique_par_object_id() const override;
    void pack(Core::Communication::PackBuffer& data) const override;
    void unpack(Core::Communication::UnpackBuffer& buffer) override;

    virtual Core::Materials::MaterialType material_type() const = 0;

    static std::shared_ptr<Summand> factory(int matnum);

    virtual void pack_summand(Core::Communication::PackBuffer& data) const;
    virtual void unpack_summand(Core::Communication::UnpackBuffer& buffer);

    virtual void setup(int numgp, const Discret::Elements::Fibers& fibers,
        const std::optional<Discret::Elements::CoordinateSystem>& coord_system);

    virtual void update();

    virtual void specify_formulation(
        bool& isoprinc, bool& isomod, bool& anisoprinc, bool& anisomod, bool& viscogeneral);

    virtual void add_coefficients_visco_principal(const Core::LinAlg::Matrix<3, 1>& inv,
        Core::LinAlg::Matrix<8, 1>& mu, Core::LinAlg::Matrix<33, 1>& xi,
        Core::LinAlg::Matrix<7, 1>& rateinv, const Teuchos::ParameterList& params, double dt,
        int gp, int eleGID);

    virtual void add_coefficients_visco_modified(const Core::LinAlg::Matrix<3, 1>& modinv,
        Core::LinAlg::Matrix<8, 1>& modmu, Core::LinAlg::Matrix<33, 1>& modxi,
        Core::LinAlg::Matrix<7, 1>& modrateinv, const Teuchos::ParameterList& params, double dt,
        int gp, int eleGID);

    virtual void read_material_parameters_visco(
        double& tau, double& beta, double& alpha, std::string& solve);

    virtual void read_material_parameters(
        int& numbranch, const std::vector<int>*& matids, std::string& solve);

    virtual void read_material_parameters(double& tau, int& matid);

    virtual void specify_visco_formulation(
        bool& visco_iso_rate, bool& visco_generalized_maxwell, bool& visco_fsls);
  };
}  // namespace Mat::ViscoElast

FOUR_C_NAMESPACE_CLOSE

#endif
