// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_FSI_DEBUGWRITER_HPP
#define FOUR_C_FSI_DEBUGWRITER_HPP

#include "4C_config.hpp"

#include "4C_coupling_adapter.hpp"
#include "4C_linalg_vector.hpp"

#include <memory>

FOUR_C_NAMESPACE_OPEN

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

  namespace Utils
  {
    /// Helper to write unconverged results at the FSI interface
    class DebugWriter
    {
     public:
      /// create FSI debug writer on given field discretization
      explicit DebugWriter(std::shared_ptr<Core::FE::Discretization> dis);

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
      void write_vector(const std::string& name, const Core::LinAlg::Vector<double>& v);

     private:
      /// internal interface discretization
      std::shared_ptr<Core::FE::Discretization> dis_;

      /// coupling of field discretization and interface discretization
      std::shared_ptr<Coupling::Adapter::Coupling> coup_;

      /// current control file
      std::shared_ptr<Core::IO::OutputControl> control_;

      /// writer to control file
      std::shared_ptr<Core::IO::DiscretizationWriter> writer_;

      /// internal FSI iteration count
      int itnum_;
    };



    /// Helper to write unconverged results
    class SimpleDebugWriter
    {
     public:
      SimpleDebugWriter(std::shared_ptr<Core::FE::Discretization> dis, const std::string& name);

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
      virtual void write_vector(const std::string& name, Core::LinAlg::Vector<double>& v);

     private:
      /// discretization
      std::shared_ptr<Core::FE::Discretization> dis_;

      /// name of output file
      std::string name_;

      /// current control file
      std::shared_ptr<Core::IO::OutputControl> control_;

      /// writer to control file
      std::shared_ptr<Core::IO::DiscretizationWriter> writer_;

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
      virtual void write_vector(
          const std::string& name, const std::shared_ptr<Core::LinAlg::Vector<double>>& v);

     private:
      FSI::Monolithic& algorithm_;

      std::shared_ptr<SimpleDebugWriter> struct_writer_;
      std::shared_ptr<SimpleDebugWriter> fluid_writer_;
      std::shared_ptr<SimpleDebugWriter> ale_writer_;

      int counter_;
    };
  }  // namespace Utils
}  // namespace FSI

FOUR_C_NAMESPACE_CLOSE

#endif
