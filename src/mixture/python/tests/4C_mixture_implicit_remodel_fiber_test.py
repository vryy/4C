# This file is part of 4C multiphysics licensed under the
# GNU Lesser General Public License v3.0 or later.
#
# See the LICENSE.md file in the top-level for license information.
#
# SPDX-License-Identifier: LGPL-3.0-or-later

import py4C
import pytest


def test_implicit_remodel_fiber():
    material = py4C.mixture.ExponentialFiberMaterial(1.0, 1.0, True)
    gnr = py4C.mixture.LinearGrowthPoissonTurnover(1.0, 10.0)
    fiber = py4C.mixture.ImplicitRemodelFiber(material, gnr, 1.1)

    fiber.recompute_state(1.11, 0.1)
    fiber.update()
    fiber.recompute_state(1.12, 0.1)

    assert fiber.cauchy_stress == pytest.approx(1.626608880501962)
    assert fiber.growth_scalar == pytest.approx(1.395405107490959)
    assert fiber.dcauchy_stress_dlambda == pytest.approx(14.034275247072243)
    assert fiber.dgrowth_scalar_dlambda == pytest.approx(2.055629068748841)
    assert fiber.lambda_r == pytest.approx(0.9286126849322407)
