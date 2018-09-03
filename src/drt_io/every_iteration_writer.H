/*----------------------------------------------------------------------------*/
/*!
\file every_iteration_writer.H

\brief Write output for each Newton step during one load step in
       an extra output file.

\maintainer Michael Hiermeier

\date Aug 10, 2017

\level 1

*/
/*----------------------------------------------------------------------------*/

#ifndef SRC_DRT_IO_EVERY_ITERATION_WRITER_H_
#define SRC_DRT_IO_EVERY_ITERATION_WRITER_H_

#include "../drt_lib/drt_dserror.H"
#include <Teuchos_RCP.hpp>

namespace Teuchos
{
  class ParameterList;
}
namespace IO
{
  class DiscretizationWriter;
  class OutputControl;

  /*--------------------------------------------------------------------------*/
  /** @brief Abstract interface class for the EveryIterationWriter
   *
   *  You have to give the @class EveryIterationWriter access to your Output
   *  routines. All methods which have to be accessed are summarized in this
   *  EveryIterationWriterInterface class. Just add this class as base to your
   *  inheritance hierarchy of the class which contains these functions.
   *  See STR::TIMINT::Base for an example.
   *
   *  @author hiermeier @date 08/17 */
  class EveryIterationWriterInterface
  {
   public:
    /// empty constructor
    EveryIterationWriterInterface(){};

    /// destructor
    virtual ~EveryIterationWriterInterface(){};

    /// Prepare the output before it is accessed via OutputState.
    virtual void PrepareOutput() = 0;

    /** @brief Access the entire output state
     *
     *  @param[in] iowriter     EveryIterationWriter object.
     *  @param[in] writer_owner If TRUE the owner of each node and element
     *                          is written.
     *
     *  @author hiermeier @date 08/17 */
    virtual void OutputState(IO::DiscretizationWriter& iowriter, bool write_owner) const = 0;

    /** @brief Get the current time/load step number
     *
     *  Np stands for n+1, i.e. the first step is equal to one. This is just a
     *  definition thing, you can return any kind of meaningful unique
     *  identification number.
     *
     *  @author hiermeier @date 08/17 */
    virtual int GetStepNp() const = 0;
  };

  /*--------------------------------------------------------------------------*/
  /** @brief Write output for each Newton step during one load step in
   *  an extra output file.
   *
   *  This class is designed for debugging purposes only. The corresponding
   *  input parameters are summarized in the section "---IO/EVERY ITERATION".
   *  The implementation is meant to work as plug-in solution for all kinds of
   *  problems. Anyway, the first implemented use can be found in the structural
   *  time integration, such that it can be seen as a blue print.
   *
   *  In general, if you consider to use this class, you will have to construct
   *  an object and call Init() and Setup() in consecutive order. Pass the
   *  output-writer of your problem to the Init() function. This output
   *  writer will be copied, such that the every iteration output should generate
   *  the same output each Newton step just as the parent writer does in each
   *  time/load step.
   *
   *  Secondly, you have to give this class access to your Output routines.
   *  All methods which have to be accessed are summarized in the
   *  EveryIterationWriterInterface class. Just add the abstract interface class
   *  as base to your inheritance hierarchy of the class which contains these
   *  functions. See STR::TIMINT::Base for an example.
   *
   *  Finally, such that the output is actually written you have to call the
   *  public class methods at the correct point in your implementation.
   *  The method InitNewtonIteration() should be called before you enter your
   *  Newton loop. AddNewtonIteration is supposed to be called just in the end
   *  of each Newton iteration.
   *
   *  You can find your every iteration output in an extra directory located in
   *  your specified output folder. For example: If your output path is
   *  "../o/fancy_sim", the every iteration output is written to
   *  "../o/fancy_sim_every_iter/fancy_sim_step_<CURRENT_STEP_NUMBER>". Each
   *  step gets an extra control, result and mesh file. It is also possible to
   *  restrict the output to one single time/load step. If you perform different
   *  runs of the same simulation and the same load/time step, just use the
   *  input parameter RUN_NUMBER. If a positive number is provided the result of
   *  each run can be found in the directory
   *
   *            "../o/fancy_sim_every_iter/run_<RUN_NUMBER>/".
   *
   *  Finally, a short note related to restarted steps: Think of a scenario
   *  where you want to restart step 2 of your simulation "fancy_sim" and you
   *  want to write the corresponding Newton iteration output as well. In this
   *  case the directory path stays untouched, while the final file name is
   *  extended by the restart counter number. Therefore, in the first restart
   *  run the entire path will be
   *
   *      "../o/fancy_sim_every_iter/run_<RUN_NUMBER>/fancy_sim-1_step_2".
   *
   *  This makes sense since you will not overwrite any previous (debugging)
   *  results which were generated by a possible completely different head file.
   *
   *  For more information see the mentioned input section.
   *
   *  That's all. I hope this will help you to find some nasty bugs.
   *
   *  @author hiermeier @date 08/17 */
  class EveryIterationWriter
  {
   public:
    /// Empty default constructor
    EveryIterationWriter();

    /// destructor
    ~EveryIterationWriter(){};

    /** @brief Initialize the object
     *
     *  @param[in] parent_writer  Discretization which is currently used to
     *                            generate the load/time step output.
     *  @param[in] interface      Access to the problem specific output routines
     *                            and load/time step numbers.
     *  @param[in] params         EVERY ITERATION parameter list.
     *
     *  @author hiermeier @date 08/17 */
    void Init(const IO::DiscretizationWriter* parent_writer,
        EveryIterationWriterInterface* interface, const Teuchos::ParameterList& params);

    /// Setup the class object
    void Setup();

    /** @brief Initialize a new Newton loop
     *
     *  This routine has to be called at the beginning of a new Newton loop.
     *  It will generate all necessary control, result and mesh files.
     *
     *  @author hiermeier @date 08/17 */
    void InitNewtonIteration();

    /** @brief Add a Newton iteration to the current output
     *
     *  Call this method in the very end of each Newton step.
     *
     *  @param[in] newton_iteration  Number of the current Newton iteration.
     *
     *  @author hiermeier @date 08/17 */
    void AddNewtonIteration(const int newton_iteration);

    /** @brief Add a line search iteration to the current output
     *
     *  Call this method before the step length will be modified.
     *
     *  @param[in] newton_iteration  Number of the current Newton iteration.
     *  @param[in] linesearch_iteration  Number of the current line search iteration.
     *
     *  @author hiermeier @date 11/17 */
    void AddLineSearchIteration(const int newton_iteration, const int linesearch_iteration);

   private:
    /// Throw if Init() has not been called.
    inline void ThrowIfNotInitialized(const int line) const
    {
      if (not isinit_) dserror("LINE %d: Call Init() first!", line);
    }

    /// Throw if Setup() has not been called.
    inline void ThrowIfNotSetup(const int line) const
    {
      if (not issetup_) dserror("LINE %d: Call Setup() first!", line);
    }

    /// Returns true if the current load/time step is supposed to be written.
    bool WriteThisStep() const;

    /** Create a new run directory if the corresponding input parameter
     *  is specified. */
    std::string CreateRunDirectory(const std::string& file_dir_path) const;

    /// Create a new directory.
    void CreateDirectory(const std::string& dir_path) const;

    /// Adjust the specified steps per file.
    void AdjustStepsPerFile(IO::OutputControl& control) const;

    /// Read-only access to the parent discretization writer.
    const IO::DiscretizationWriter& ParentWriter() const
    {
      if (not parent_writer_) dserror("The parent writer has not been initialized correctly.");

      return *parent_writer_;
    }

    /// Access to the interface object.
    EveryIterationWriterInterface& Interface()
    {
      if (not interface_)
        dserror("The every iteration interface has not been initialized correctly.");

      return *interface_;
    }

    /// print current output path to the os
    void PrintPath2Screen(const std::string& path) const;

   private:
    bool isinit_;
    bool issetup_;
    bool isnewton_initialized_;
    int myrank_;
    int run_number_;
    int write_only_this_step_;
    bool write_owner_each_newton_iteration_;
    std::string base_filename_;

    const IO::DiscretizationWriter* parent_writer_;

    EveryIterationWriterInterface* interface_;

    Teuchos::RCP<IO::DiscretizationWriter> every_iter_writer_;

    constexpr static unsigned MAX_NUMBER_LINE_SEARCH_ITERATIONS_ = 100;
  };

  /// Extract the path of a full filename.
  std::string ExtractPath(const std::string& full_filename);

  std::string ExtractFileName(const std::string& full_filename);

  /// remove the restart counter from the file name
  std::string RemoveRestartCountFromFileName(const std::string& filename);

  /// Create a new directory.
  void CreateDirectory(const std::string& dir_path, const int myrank);

  /// Return number of lines in the given file
  int CountLinesInFile(const std::string& filepath);
}  // namespace IO


#endif /* SRC_DRT_IO_EVERY_ITERATION_WRITER_H_ */
