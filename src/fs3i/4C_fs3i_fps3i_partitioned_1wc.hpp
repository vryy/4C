/*----------------------------------------------------------------------*/
/*! \file
\brief H-file associated with algorithmic routines for partitioned
     solution approaches to fluid-porous-structure-scalar-scalar interaction
     (FPS3I) specifically related to one-way-coupled problem
     configurations

\level 3



*----------------------------------------------------------------------*/


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
    PartFpS3I1Wc(const Epetra_Comm& comm);

    //! initialize this class
    void Init() override;

    //! setup this class
    void Setup() override;

    /// timeloop of coupled problem
    void Timeloop() override;

    /// FPSI step
    void DoFPSIStep();

    /// Scatra step
    void DoScatraStep();

    //! routine for preparing time step
    virtual void PrepareTimeStep();

    //! check convergence of monolithic ScaTra problem
    virtual bool scatra_convergence_check(int itnum);
  };
}  // namespace FS3I

FOUR_C_NAMESPACE_CLOSE

#endif
