/*----------------------------------------------------------------------*/
/*! \file
\brief write debug information for fsi applications
\level 2
*----------------------------------------------------------------------*/

#ifndef FOUR_C_FSI_DEBUGWRITER_HPP
#define FOUR_C_FSI_DEBUGWRITER_HPP

#include "baci_config.hpp"

#include "baci_coupling_adapter.hpp"

#include <Epetra_Vector.h>
#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN

namespace ADAPTER
{
  class Coupling;
}

namespace IO
{
  class OutputControl;
  class DiscretizationWriter;
}  // namespace IO

namespace DRT
{
  class Discretization;
}

namespace FSI
{
  class Monolithic;

  namespace UTILS
  {
    /// Helper to write unconverged results at the FSI interface
    class DebugWriter
    {
     public:
      /// create FSI debug writer on given field discretization
      explicit DebugWriter(Teuchos::RCP<DRT::Discretization> dis);

      /// announce new time step
      /*!
        Results in a new control file
       */
      void NewTimeStep(int step, std::string name = "");

      /// announce new interface iteration step
      void NewIteration();

      /// write interface vector
      /*!
        This includes an internal transfer of the interface vector to the
        interface discretization.
       */
      void WriteVector(const std::string& name, const Epetra_Vector& v);

     private:
      /// internal interface discretization
      Teuchos::RCP<DRT::Discretization> dis_;

      /// coupling of field discretization and interface discretization
      Teuchos::RCP<CORE::ADAPTER::Coupling> coup_;

      /// current control file
      Teuchos::RCP<IO::OutputControl> control_;

      /// writer to control file
      Teuchos::RCP<IO::DiscretizationWriter> writer_;

      /// internal FSI iteration count
      int itnum_;
    };



    /// Helper to write unconverged results
    class SimpleDebugWriter
    {
     public:
      SimpleDebugWriter(Teuchos::RCP<DRT::Discretization> dis, const std::string& name);

      virtual ~SimpleDebugWriter() = default;
      /// announce new time step
      /*!
        Results in a new control file
       */
      virtual void NewLinearSystem(int step, std::string name = "");

      /// announce new interface iteration step
      virtual void NewIteration();

      /// write interface vector
      /*!
        This includes an internal transfer of the interface vector to the
        interface discretization.
       */
      virtual void WriteVector(const std::string& name, Epetra_Vector& v);

     private:
      /// discretization
      Teuchos::RCP<DRT::Discretization> dis_;

      /// name of output file
      std::string name_;

      /// current control file
      Teuchos::RCP<IO::OutputControl> control_;

      /// writer to control file
      Teuchos::RCP<IO::DiscretizationWriter> writer_;

      /// internal FSI iteration count
      int itnum_;
    };


    /// FSI debug writer
    class MonolithicDebugWriter
    {
     public:
      MonolithicDebugWriter(FSI::Monolithic& algorithm);

      virtual ~MonolithicDebugWriter() = default;
      /// announce new time step
      /*!
        Results in a new control file
       */
      virtual void NewLinearSystem();

      /// announce new interface iteration step
      virtual void NewIteration();

      /// write interface vector
      /*!
        This includes an internal transfer of the interface vector to the
        interface discretization.
       */
      virtual void WriteVector(const std::string& name, const Teuchos::RCP<Epetra_Vector>& v);

     private:
      FSI::Monolithic& algorithm_;

      Teuchos::RCP<SimpleDebugWriter> struct_writer_;
      Teuchos::RCP<SimpleDebugWriter> fluid_writer_;
      Teuchos::RCP<SimpleDebugWriter> ale_writer_;

      int counter_;
    };
  }  // namespace UTILS
}  // namespace FSI

FOUR_C_NAMESPACE_CLOSE

#endif
