#ifdef CCADISCRET

#ifndef FS3I_PART_2WC_H
#define FS3I_PART_2WC_H

#include "fs3i_partitioned.H"

namespace FS3I
{
  class PartFS3I_2WC : public PartFS3I
  {
  public:

    PartFS3I_2WC(const Epetra_Comm& comm);

    virtual ~PartFS3I_2WC(){};

    virtual void Timeloop();

    void InitialCalculations();

    void PrepareTimeStep();

    void OuterLoop();

    //void SetFSIValuesInScaTra();

    void SetScaTraValuesInFSI();

    bool ConvergenceCheck(int itnum);

    virtual bool ScatraConvergenceCheck(int itnum);

    void TimeUpdateAndOutput();

  private:

    //! @name  (preliminary) maximum number of iterations and tolerance for outer iteration
    int     itmax_;
    double  ittol_;
    //@}

    /// flag for constant thermodynamic pressure
    string consthermpress_;

    /// fluid- and structure-based scalar transport problem
    Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> fluidscatra_;
    Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> structurescatra_;
  };
}

#endif
#endif
