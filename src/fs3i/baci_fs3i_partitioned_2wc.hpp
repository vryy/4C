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


#include "baci_config.hpp"

#include "baci_fs3i_partitioned.hpp"

FOUR_C_NAMESPACE_OPEN


namespace FS3I
{
  class PartFS3I_2WC : public PartFS3I
  {
   public:
    PartFS3I_2WC(const Epetra_Comm& comm);

    //! initialize this class
    void Init() override;

    //! setup this class
    void Setup() override;

    void Timeloop() override;

    void InitialCalculations();

    void PrepareTimeStep() override;

    void OuterLoop();

    // void SetFSIValuesInScaTra();

    void SetScaTraValuesInFSI();

    bool ConvergenceCheck(int itnum);

    bool ScatraConvergenceCheck(int itnum) override;

    void TimeUpdateAndOutput();

   private:
    //! @name  (preliminary) maximum number of iterations and tolerance for outer iteration
    int itmax_;
    double ittol_;
    //@}

    /// flag for constant thermodynamic pressure
    std::string consthermpress_;

    /// fluid- and structure-based scalar transport problem
    Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> fluidscatra_;
    Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> structurescatra_;
  };
}  // namespace FS3I

FOUR_C_NAMESPACE_CLOSE

#endif
