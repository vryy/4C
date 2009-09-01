/*----------------------------------------------------------------------*/
/*!
\file drt_globalproblem.H

\brief global list of problems

<pre>
Maintainer: Ulrich Kuettler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15238
</pre>
*/
/*----------------------------------------------------------------------*/

#ifdef CCADISCRET
#ifndef DRT_GLOBALPROBLEM_H
#define DRT_GLOBALPROBLEM_H

#include <vector>

#include <Epetra_Comm.h>
#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterListAcceptorDefaultBase.hpp>

#include "../drt_mat/matpar_bundle.H"
#include "drt_discret.H"
#include "../drt_nurbs_discret/drt_nurbs_discret.H"
#include "drt_inputreader.H"
#include "drt_function.H"
#include "drt_timecurve.H"
#include "drt_resulttest.H"

#include "standardtypes_cpp.H"

using namespace std;
using namespace Teuchos;

namespace IO
{
  class OutputControl;
  class InputControl;
  class ErrorFileControl;
}

namespace MAT
{
namespace PAR
{
class Bundle;
}
}


namespace DRT
{

namespace INPUT
{
class DatFileReader;
}


/// Global problem instance that keeps the discretizations
/*!

  The global problem represents the input file passed to baci. This
  class organizes the reading of a dat file (utilizing the
  DatFileReader of course). That is way, in all but the most eccentric
  cases there will be exactly one object of this class during a baci
  run. This object contains all parameters read from the input file as
  well as any material definitions and even all the discretizations.

  <h3>Input parameters</h3>

  All input parameters are known by the global problem object. These
  parameters are guaranteed to be valid (because they passed the
  validation) and are guaranteed to be there (because default values
  have been set for all parameters missing from the input file.) This
  is Teuchos::ParameterList magic, that just requires the list of
  valid parameters in the file drt_validparameters.cpp to be
  complete.

  The algorithms are meant to ask the global problem object for their
  parameters and extract them from the respective parameter list
  directly.

  <h3>Discretizations</h3>

  The global problem object knows the discretizations defined by the
  input file. In order to access a particular discretization you get
  the global problem object and ask.

  <h3>Materials</h3>

  The global problem object knows the material descriptions from the
  input file. These are not to be confused with the material classes
  the elements know and work with. The global problem object does not
  keep track of gauss point material values, all that is known here
  are the definitions from the input file.

  <h3>Purpose</h3>

  baci evolved out of ccarat in a (somewhat) smooth fashion. The old
  ccarat used global C struct variables all over the place. It is
  indeed the aim of the global problem to collect these variables in
  one place and provide a sane interface to the data. The transition
  is already pretty advanced, however there are still many places
  where global variables are used. These places will have to be
  revised some day.

  The global problem class is still evolving.

  \note It is illegal to introduce new global variables in baci. It is
  illegal to used global variables where other means to do the job are
  available. We just did not get around to remove all of them yet.

  <h3>Singleton behavior</h3>

  \warning This is a guru only section!

  The global problem behaves like a singleton, so there is always one
  instance available. But you can have more than one instance of
  Problem. In normal situations this will not be needed. So don't
  bother. Just call the static Instance() function to get the global
  instance and access your discretizations.

  In the special case that you want to read more that one input file,
  however, you will need to handle the fields from each file
  separately. This is something you cannot do with ccarat's field
  variable so in that case you will have to rely on a second Problem
  instance.

  One artefact that comes from using global variables together with
  multiple Problem objects is the notion of activating of problem
  object. This translates to setting the global variable pointers to
  this Problem object's internal variables. Normally, if there is just
  one Problem object, this is done by default. If you need more that
  one, however, you will have to activate the global problems
  yourself.

  \author u.kue
  \date 06/07
 */
class Problem : public Teuchos::ParameterListAcceptorDefaultBase
{
public:

  /// @name Instances

  /// return an instance of this class
  static RCP<Problem> Instance(int num=0);

  /// return number of problem instances
  static unsigned NumInstances() { return instances_.size(); }

  /// calculation done, clean up
  /*!
    There can be a variety of objects to a problem. Some of them might
    require proper cleanup. Make sure we always do it.
   */
  static void Done();

  //@}

  /// @name Input

  /// read global parameters
  void ReadParameter(DRT::INPUT::DatFileReader& reader);

  /// input of not mesh or time based problem data
  void InputControl();

  /// setup the discretizations
  void ReadFields(DRT::INPUT::DatFileReader& reader);

  /// input of materials
  void ReadMaterials(const DRT::INPUT::DatFileReader& reader);

  /// input of conditions
  void ReadConditions(const DRT::INPUT::DatFileReader& reader);

  /// input of knots for isogeometric analysis
  void ReadKnots(const DRT::INPUT::DatFileReader& reader);

  //@}

  /// @name General query methods
  /// Once and for all definitions

  /// give string name of my problem type
  std::string ProblemType() const;

  // give string name of my basis function type (polynomial/nurbs)
  std::string SpatialApproximation() const;

  //@}

  /// @name Control file

  /// create control file for output and read restart data if required
  void OpenControlFile(const Epetra_Comm& comm, std::string inputfile, std::string prefix);

  /// error files
  void OpenErrorFile(const Epetra_Comm& comm, std::string prefix);

  /// control file for restart read
  Teuchos::RCP<IO::InputControl> InputControlFile() { return inputcontrol_; }

  /// control file for normal output
  Teuchos::RCP<IO::OutputControl> OutputControlFile() { return outputcontrol_; }

  /// control error files
  Teuchos::RCP<IO::ErrorFileControl> ErrorFile() { return errorfilecontrol_; }

  /// write parameters read from dat file for documentation
  void WriteInputParameters();

  //@}

  /// @name Parameters read from file

  /// Set parameters from a parameter list and return with default values.
  void setParameterList(Teuchos::RCP< Teuchos::ParameterList > const &paramList);

  /// Return a const parameter list of all of the valid parameters that this->setParameterList(...) will accept.
  Teuchos::RCP<const Teuchos::ParameterList> getValidParameters() const;

  //@}

  /// @name Input parameter sections
  /// direct access to parameters from input file sections

  const Teuchos::ParameterList& ProblemTypeParams() const { return getParameterList()->sublist("PROBLEM TYP"); }
  const Teuchos::ParameterList& ProblemSizeParams() const { return getParameterList()->sublist("PROBLEM SIZE"); }
  const Teuchos::ParameterList& IOParams() const { return getParameterList()->sublist("IO"); }
  const Teuchos::ParameterList& DiscretisationParams() const { return getParameterList()->sublist("DISCRETISATION"); }
  const Teuchos::ParameterList& DesignDescriptionParams() const { return getParameterList()->sublist("DESIGN DESCRIPTION"); }
  const Teuchos::ParameterList& StructuralDynamicParams() const { return getParameterList()->sublist("STRUCTURAL DYNAMIC"); }
  const Teuchos::ParameterList& InverseAnalysisParams() const { return getParameterList()->sublist("INVERSE ANALYSIS"); }
  const Teuchos::ParameterList& StructuralContactParams() const { return getParameterList()->sublist("STRUCTURAL CONTACT"); }
  const Teuchos::ParameterList& InteractionPotentialParams() const { return getParameterList()->sublist("INTERACTION POTENTIALS"); }
  const Teuchos::ParameterList& BrownianMotionParams() const { return getParameterList()->sublist("BROWNIAN MOTION"); }
  const Teuchos::ParameterList& StatisticalMechanicsParams() const { return getParameterList()->sublist("STATISTICAL MECHANICS"); }
  const Teuchos::ParameterList& ThermalDynamicParams() const { return getParameterList()->sublist("THERMAL DYNAMIC"); }
  const Teuchos::ParameterList& FluidDynamicParams() const { return getParameterList()->sublist("FLUID DYNAMIC"); }
  const Teuchos::ParameterList& ScalarTransportDynamicParams() const { return getParameterList()->sublist("SCALAR TRANSPORT DYNAMIC"); }
  const Teuchos::ParameterList& CombustionDynamicParams() const { return getParameterList()->sublist("COMBUSTION CONTROL"); }
  const Teuchos::ParameterList& AleDynamicParams() const { return getParameterList()->sublist("ALE DYNAMIC"); }
  const Teuchos::ParameterList& FSIDynamicParams() const { return getParameterList()->sublist("FSI DYNAMIC"); }
  const Teuchos::ParameterList& XFEMGeneralParams() const { return getParameterList()->sublist("XFEM GENERAL"); }
  const Teuchos::ParameterList& LOMAControlParams() const { return getParameterList()->sublist("LOMA CONTROL"); }
  const Teuchos::ParameterList& ELCHControlParams() const { return getParameterList()->sublist("ELCH CONTROL"); }
  const Teuchos::ParameterList& ArterialDynamicParams() const { return getParameterList()->sublist("ARTERIAL DYNAMIC"); }

  const Teuchos::ParameterList& StructuralNoxParams() const { return getParameterList()->sublist("STRUCT NOX"); }

  const Teuchos::ParameterList& StructSolverParams() const { return getParameterList()->sublist("STRUCT SOLVER"); }
  const Teuchos::ParameterList& FluidSolverParams() const { return getParameterList()->sublist("FLUID SOLVER"); }
  const Teuchos::ParameterList& FluidPressureSolverParams() const { return getParameterList()->sublist("FLUID PRESSURE SOLVER"); }
  const Teuchos::ParameterList& AleSolverParams() const { return getParameterList()->sublist("ALE SOLVER"); }
  const Teuchos::ParameterList& ThermalSolverParams() const { return getParameterList()->sublist("THERMAL SOLVER"); }
  const Teuchos::ParameterList& ScalarTransportSolverParams() const { return getParameterList()->sublist("SCALAR TRANSPORT SOLVER"); }
  const Teuchos::ParameterList& ScalarTransportElectricPotentialSolverParams() const { return getParameterList()->sublist("SCALAR TRANSPORT ELECTRIC POTENTIAL SOLVER"); }
  const Teuchos::ParameterList& ArteryNetworkSolverParams() const { return getParameterList()->sublist("ARTERY NETWORK SOLVER"); }
  const Teuchos::ParameterList& ConditionNamesParams() const { return getParameterList()->sublist("CONDITION NAMES"); }

  //@}

  /// @name Discretizations

  /// find a particular discretization
  RCP<Discretization> Dis(int fieldnum, int disnum) const;

  /// tell number of known fields
  unsigned NumFields() const { return discretizations_.size(); }

  /// tell number of discretizations in a field
  unsigned NumDis(int fieldnum) const { return discretizations_[fieldnum].size(); }

  /// add a discretization to a specific field
  void AddDis(int fieldnum, RCP<Discretization> dis);

  /// set given discretization to specific field at specific position
  void SetDis(int fieldnum, int disnum, RCP<Discretization> dis);

  //@}

  /// @name Materials

  /// return pointer to materials bundled to the problem
  Teuchos::RCP<MAT::PAR::Bundle> Materials() { return materials_; }

  //@}

  /// @name Spatial Functions

  /// return a specified function
  DRT::UTILS::Function& Funct(int num);

  //@}

  /// @name Spatial Functions

  /// return a specified time curve
  DRT::UTILS::TimeCurve& Curve(int num);

  //@}

  /// @name Result Tests

  /// Do the testing
  void TestAll(const Epetra_Comm& comm) { resulttest_.TestAll(comm); }

  /// add field specific result test object
  void AddFieldTest(Teuchos::RCP<ResultTest> test) { resulttest_.AddFieldTest(test); }

  //@}

private:

  /// private default constructor to disallow creation of instances
  Problem() {}

  /// no copying
  Problem(const Problem&);

  /// no copying
  Problem& operator=(const Problem&);

  /// the single instance
  static vector<RCP<Problem> > instances_;

  /// discretizations of this problem
  vector<vector<RCP<Discretization> > > discretizations_;

  /// material bundle
  Teuchos::RCP<MAT::PAR::Bundle> materials_;

  /// basket of spatial function
  DRT::UTILS::FunctionManager functionmanager_;

  /// basket of time curves
  DRT::UTILS::TimeCurveManager timecurvemanager_;

  /// all test values we might have
  DRT::ResultTestManager resulttest_;

  /// @name File IO

  Teuchos::RCP<IO::InputControl>  inputcontrol_;
  Teuchos::RCP<IO::OutputControl> outputcontrol_;
  Teuchos::RCP<IO::ErrorFileControl> errorfilecontrol_;  ///< error file handles

  //@}
};

}

#endif
#endif
