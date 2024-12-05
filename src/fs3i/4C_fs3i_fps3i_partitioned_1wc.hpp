// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_FS3I_FPS3I_PARTITIONED_1WC_HPP
#define FOUR_C_FS3I_FPS3I_PARTITIONED_1WC_HPP


#include "4C_config.hpp"

#include "4C_fs3i_fps3i_partitioned.hpp"

FOUR_C_NAMESPACE_OPEN


namespace FS3I
{
  class PartFpS3I1Wc : public PartFPS3I
  {
   public:
    //! constructor of one-way coupled FPS3I
    PartFpS3I1Wc(MPI_Comm comm);

    //! initialize this class
    void init() override;

    //! setup this class
    void setup() override;

    /// timeloop of coupled problem
    void timeloop() override;

    /// FPSI step
    void do_fpsi_step();

    /// Scatra step
    void do_scatra_step();

    //! routine for preparing time step
    virtual void prepare_time_step();

    //! check convergence of monolithic ScaTra problem
    virtual bool scatra_convergence_check(int itnum);
  };
}  // namespace FS3I

FOUR_C_NAMESPACE_CLOSE

#endif
