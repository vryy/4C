/*----------------------------------------------------------------------*/
/*! \file
\brief write debug information for fsi applications
\level 2
*----------------------------------------------------------------------*/

#ifndef FOUR_C_FSI_DEBUGWRITER_HPP
#define FOUR_C_FSI_DEBUGWRITER_HPP

#include "4C_config.hpp"

#include "4C_coupling_adapter.hpp"

#include <Epetra_Vector.h>
#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN

namespace Adapter
{
  class Coupling;
}

namespace Core::IO
{
  class OutputControl;
  class DiscretizationWriter;
}  // namespace Core::IO

namespace Core::FE
{
  class Discretization;
}  // namespace Core::FE

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
      explicit DebugWriter(Teuchos::RCP<Core::FE::Discretization> dis);

      /// announce new time step
      /*!
        Results in a new control file
       */
      void new_time_step(int step, std::string name = "");

      /// announce new interface iteration step
      void new_iteration();

      /// write interface vector
      /*!
        This includes an internal transfer of the interface vector to the
        interface discretization.
       */
      void write_vector(const std::string& name, const Epetra_Vector& v);

     private:
      /// internal interface discretization
      Teuchos::RCP<Core::FE::Discretization> dis_;

      /// coupling of field discretization and interface discretization
      Teuchos::RCP<Core::Adapter::Coupling> coup_;

      /// current control file
      Teuchos::RCP<Core::IO::OutputControl> control_;

      /// writer to control file
      Teuchos::RCP<Core::IO::DiscretizationWriter> writer_;

      /// internal FSI iteration count
      int itnum_;
    };



    /// Helper to write unconverged results
    class SimpleDebugWriter
    {
     public:
      SimpleDebugWriter(Teuchos::RCP<Core::FE::Discretization> dis, const std::string& name);

      virtual ~SimpleDebugWriter() = default;
      /// announce new time step
      /*!
        Results in a new control file
       */
      virtual void new_linear_system(int step, std::string name = "");

      /// announce new interface iteration step
      virtual void new_iteration();

      /// write interface vector
      /*!
        This includes an internal transfer of the interface vector to the
        interface discretization.
       */
      virtual void write_vector(const std::string& name, Epetra_Vector& v);

     private:
      /// discretization
      Teuchos::RCP<Core::FE::Discretization> dis_;

      /// name of output file
      std::string name_;

      /// current control file
      Teuchos::RCP<Core::IO::OutputControl> control_;

      /// writer to control file
      Teuchos::RCP<Core::IO::DiscretizationWriter> writer_;

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
      virtual void new_linear_system();

      /// announce new interface iteration step
      virtual void new_iteration();

      /// write interface vector
      /*!
        This includes an internal transfer of the interface vector to the
        interface discretization.
       */
      virtual void write_vector(const std::string& name, const Teuchos::RCP<Epetra_Vector>& v);

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
