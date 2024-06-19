/*----------------------------------------------------------------------*/
/*! \file
\brief H-file associated with algorithmic routines for partitioned
       solution approaches to fluid-structure-scalar-scalar interaction
       (FS3I) specifically related to two-way-coupled problem
       configurations

\level 2



*----------------------------------------------------------------------*/


#ifndef FOUR_C_FS3I_PARTITIONED_2WC_HPP
#define FOUR_C_FS3I_PARTITIONED_2WC_HPP


#include "4C_config.hpp"

#include "4C_fs3i_partitioned.hpp"

FOUR_C_NAMESPACE_OPEN


namespace FS3I
{
  class PartFS3I2Wc : public PartFS3I
  {
   public:
    PartFS3I2Wc(const Epetra_Comm& comm);

    //! initialize this class
    void init() override;

    //! setup this class
    void setup() override;

    void Timeloop() override;

    void initial_calculations();

    void prepare_time_step() override;

    void outer_loop();

    // void SetFSIValuesInScaTra();

    void set_sca_tra_values_in_fsi();

    bool convergence_check(int itnum);

    bool scatra_convergence_check(int itnum) override;

    void TimeUpdateAndOutput();

   private:
    //! @name  (preliminary) maximum number of iterations and tolerance for outer iteration
    int itmax_;
    double ittol_;
    //@}

    /// flag for constant thermodynamic pressure
    std::string consthermpress_;

    /// fluid- and structure-based scalar transport problem
    Teuchos::RCP<Adapter::ScaTraBaseAlgorithm> fluidscatra_;
    Teuchos::RCP<Adapter::ScaTraBaseAlgorithm> structurescatra_;
  };
}  // namespace FS3I

FOUR_C_NAMESPACE_CLOSE

#endif
