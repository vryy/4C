#ifndef STR_GENALPHA_H
#define STR_GENALPHA_H

#ifdef STKADAPTIVE

#include "../stk_lib/stk_algorithm.H"
#include "../stk_refine/stk_io.H"
#include "../stk_lib/stk_linearsystem.H"

namespace LINALG
{
  class FixedSparseMatrix;
  class Solver;
}

namespace DRT
{
  class ResultTest;
}

namespace STK
{
  class Discretization;

  namespace STR
  {

    /// central fluid algorithm
    class Structure : public STK::Algorithm
    {
      friend struct StructureState;
    public:

      Structure( STK::Discretization & dis, Teuchos::RCP<LINALG::Solver> solver );

      STK::Discretization & Discretization() const { return dis_; }

      void Integrate();

      Teuchos::RCP<DRT::ResultTest> CreateFieldTest();

      virtual void declare_fields( stk::mesh::MetaData & meta );

      virtual void collect_unknowns( std::vector<stk::mesh::FieldBase*> & fields );

      virtual void notify_state_changed();

    private:

      void ConstantPredictor();
      void ConsistentPredictor();
      void FullNewton();
      void UpdateandOutput();

      STK::Discretization & dis_;

      //! @name time stepping variables
      double  time_;        ///< physical time
      int     step_;        ///< timestep
      int     stepmax_;     ///< maximal number of timesteps
      double  maxtime_;     ///< maximal physical computation time
      //@}

    };

  }
}

#endif
#endif
