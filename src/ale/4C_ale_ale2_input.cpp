// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_ale_ale2.hpp"
#include "4C_io_linedefinition.hpp"
#include "4C_mat_so3_material.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
bool Discret::Elements::Ale2::read_element(const std::string& eletype, const std::string& distype,
    const Core::IO::InputParameterContainer& container)
{
  // read number of material model
  int material_id = container.get<int>("MAT");
  set_material(0, Mat::factory(material_id));

  // get gauss rule
  const Core::FE::GaussRule2D gaussrule = get_optimal_gaussrule(shape());
  const Core::FE::IntegrationPoints2D intpoints(gaussrule);
  const int numgp = intpoints.nquad;

  // get material
  Teuchos::RCP<Core::Mat::Material> mat = material();
  Teuchos::RCP<Mat::So3Material> so3mat = Teuchos::rcp_dynamic_cast<Mat::So3Material>(mat, true);

  // call material setup
  so3mat->setup(numgp, container);
  return true;
}

FOUR_C_NAMESPACE_CLOSE
