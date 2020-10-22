/*---------------------------------------------------------------------------*/
/*! \file
\brief write visualization output for rigid bodies in vtk/vtp format at runtime
\level 1
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include "particle_rigidbody_runtime_vtp_writer.H"

#include "particle_rigidbody_datastate.H"

#include "../drt_io/io.H"
#include "../drt_io/io_control.H"
#include "../drt_io/runtime_vtp_writer.H"

#include "../drt_lib/drt_globalproblem.H"

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
PARTICLERIGIDBODY::RigidBodyRuntimeVtpWriter::RigidBodyRuntimeVtpWriter(const Epetra_Comm& comm)
    : comm_(comm), setuptime_(0.0), fieldname_("rigidbody")
{
  // empty constructor
}

void PARTICLERIGIDBODY::RigidBodyRuntimeVtpWriter::Init(
    const std::shared_ptr<PARTICLERIGIDBODY::RigidBodyDataState> rigidbodydatastate)
{
  // set rigid body data state container
  rigidbodydatastate_ = rigidbodydatastate;

  // construct the writer object
  runtime_vtpwriter_ = std::make_shared<RuntimeVtpWriter>();
}

void PARTICLERIGIDBODY::RigidBodyRuntimeVtpWriter::Setup(bool write_binary_output)
{
  // determine path of output directory
  const std::string outputfilename(DRT::Problem::Instance()->OutputControlFile()->FileName());
  size_t pos = outputfilename.find_last_of("/");
  if (pos == outputfilename.npos)
    pos = 0ul;
  else
    pos++;
  const std::string output_directory_path(outputfilename.substr(0ul, pos));

  // we need a better upper bound for total number of time steps here
  // however, this 'only' affects the number of leading zeros in the vtk file names
  const unsigned int max_number_timesteps_to_be_written = 1.0e+6;

  // initialize the writer object
  runtime_vtpwriter_->Initialize(comm_.MyPID(), comm_.NumProc(), max_number_timesteps_to_be_written,
      output_directory_path, DRT::Problem::Instance()->OutputControlFile()->FileNameOnlyPrefix(),
      fieldname_, DRT::Problem::Instance()->OutputControlFile()->RestartName(), setuptime_,
      write_binary_output);
}

void PARTICLERIGIDBODY::RigidBodyRuntimeVtpWriter::ReadRestart(
    const std::shared_ptr<IO::DiscretizationReader> reader)
{
  // get restart time
  setuptime_ = reader->ReadDouble("time");
}

void PARTICLERIGIDBODY::RigidBodyRuntimeVtpWriter::ResetTimeAndTimeStep(
    double time, unsigned int timestep)
{
  runtime_vtpwriter_->SetupForNewTimeStepAndGeometry(time, timestep, fieldname_);
}

void PARTICLERIGIDBODY::RigidBodyRuntimeVtpWriter::SetRigidBodyPositionsAndStates(
    const std::vector<int>& ownedrigidbodies)
{
  // rigid body position
  {
    // get and prepare storage for position data
    std::vector<double>& posdata = runtime_vtpwriter_->GetMutablePointCoordinateVector();
    posdata.clear();
    posdata.reserve(3 * ownedrigidbodies.size());

    // get reference to rigid body position
    const std::vector<std::vector<double>>& pos = rigidbodydatastate_->GetRefPosition();

    // copy rigid body position data
    for (int rigidbody_k : ownedrigidbodies)
      posdata.insert(posdata.end(), pos[rigidbody_k].begin(), pos[rigidbody_k].end());
  }

  // rigid body mass
  {
    // prepare rigid body mass data
    std::vector<double> massdata;
    massdata.reserve(ownedrigidbodies.size());

    // get reference to rigid body mass
    const std::vector<double>& mass = rigidbodydatastate_->GetRefMass();

    // copy rigid body mass data
    for (int rigidbody_k : ownedrigidbodies) massdata.push_back(mass[rigidbody_k]);

    // append rigid body mass data to vtp writer
    runtime_vtpwriter_->AppendVisualizationPointDataVector(massdata, 1, "mass");
  }

  // rigid body velocity
  {
    // prepare rigid body velocity data
    std::vector<double> veldata;
    veldata.reserve(3 * ownedrigidbodies.size());

    // get reference to rigid body velocity
    const std::vector<std::vector<double>>& vel = rigidbodydatastate_->GetRefVelocity();

    // copy rigid body velocity data
    for (int rigidbody_k : ownedrigidbodies)
      veldata.insert(veldata.end(), vel[rigidbody_k].begin(), vel[rigidbody_k].end());

    // append rigid body velocity data to vtp writer
    runtime_vtpwriter_->AppendVisualizationPointDataVector(veldata, 3, "velocity");
  }

  // rigid body acceleration
  {
    // prepare rigid body acceleration data
    std::vector<double> accdata;
    accdata.reserve(3 * ownedrigidbodies.size());

    // get reference to rigid body acceleration
    const std::vector<std::vector<double>>& acc = rigidbodydatastate_->GetRefAcceleration();

    // copy rigid body acceleration data
    for (int rigidbody_k : ownedrigidbodies)
      accdata.insert(accdata.end(), acc[rigidbody_k].begin(), acc[rigidbody_k].end());

    // append rigid body acceleration data to vtp writer
    runtime_vtpwriter_->AppendVisualizationPointDataVector(accdata, 3, "acceleration");
  }

  // rigid body angular velocity
  {
    // prepare rigid body angular velocity data
    std::vector<double> angveldata;
    angveldata.reserve(3 * ownedrigidbodies.size());

    // get reference to rigid body angular velocity
    const std::vector<std::vector<double>>& angvel = rigidbodydatastate_->GetRefAngularVelocity();

    // copy rigid body angular velocity data
    for (int rigidbody_k : ownedrigidbodies)
      angveldata.insert(angveldata.end(), angvel[rigidbody_k].begin(), angvel[rigidbody_k].end());

    // append rigid body angular velocity data to vtp writer
    runtime_vtpwriter_->AppendVisualizationPointDataVector(angveldata, 3, "angular velocity");
  }

  // rigid body angular acceleration
  {
    // prepare rigid body angular acceleration data
    std::vector<double> angaccdata;
    angaccdata.reserve(3 * ownedrigidbodies.size());

    // get reference to rigid body angular acceleration
    const std::vector<std::vector<double>>& angacc =
        rigidbodydatastate_->GetRefAngularAcceleration();

    // copy rigid body angular acceleration data
    for (int rigidbody_k : ownedrigidbodies)
      angaccdata.insert(angaccdata.end(), angacc[rigidbody_k].begin(), angacc[rigidbody_k].end());

    // append rigid body angular acceleration data to vtp writer
    runtime_vtpwriter_->AppendVisualizationPointDataVector(angaccdata, 3, "angular acceleration");
  }

  // rigid body force
  {
    // prepare rigid body force data
    std::vector<double> forcedata;
    forcedata.reserve(3 * ownedrigidbodies.size());

    // get reference to rigid body force
    const std::vector<std::vector<double>>& force = rigidbodydatastate_->GetRefForce();

    // copy rigid body force data
    for (int rigidbody_k : ownedrigidbodies)
      forcedata.insert(forcedata.end(), force[rigidbody_k].begin(), force[rigidbody_k].end());

    // append rigid body force data to vtp writer
    runtime_vtpwriter_->AppendVisualizationPointDataVector(forcedata, 3, "force");
  }

  // rigid body torque
  {
    // prepare rigid body torque data
    std::vector<double> torquedata;
    torquedata.reserve(3 * ownedrigidbodies.size());

    // get reference to rigid torque force
    const std::vector<std::vector<double>>& torque = rigidbodydatastate_->GetRefTorque();

    // copy rigid body torque data
    for (int rigidbody_k : ownedrigidbodies)
      torquedata.insert(torquedata.end(), torque[rigidbody_k].begin(), torque[rigidbody_k].end());

    // append rigid body torque data to vtp writer
    runtime_vtpwriter_->AppendVisualizationPointDataVector(torquedata, 3, "torque");
  }

  // rigid body global id
  {
    // prepare rigid body global id data
    std::vector<double> globaliddata;
    globaliddata.reserve(ownedrigidbodies.size());

    // copy rigid body global id data
    for (int rigidbody_k : ownedrigidbodies) globaliddata.push_back(rigidbody_k);

    // append rigid body global id data to vtp writer
    runtime_vtpwriter_->AppendVisualizationPointDataVector(globaliddata, 1, "globalid");
  }

  // rigid body owner
  {
    // set rigid body owner data
    std::vector<double> ownerdata(ownedrigidbodies.size(), comm_.MyPID());

    // append owner of rigid bodies to vtp writer
    runtime_vtpwriter_->AppendVisualizationPointDataVector(ownerdata, 1, "owner");
  }
}

void PARTICLERIGIDBODY::RigidBodyRuntimeVtpWriter::WriteFiles()
{
  runtime_vtpwriter_->WriteFiles();
}

void PARTICLERIGIDBODY::RigidBodyRuntimeVtpWriter::WriteCollectionFileOfAllWrittenFiles()
{
  runtime_vtpwriter_->WriteCollectionFileOfAllWrittenFiles(
      DRT::Problem::Instance()->OutputControlFile()->FileNameOnlyPrefix() + "-" + fieldname_);
}
