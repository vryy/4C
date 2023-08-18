/*---------------------------------------------------------------------------*/
/*! \file
\brief write visualization output for rigid bodies in vtk/vtp format at runtime
\level 1
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include "baci_particle_rigidbody_runtime_vtp_writer.H"

#include "baci_io.H"
#include "baci_io_control.H"
#include "baci_io_visualization_manager.H"
#include "baci_lib_globalproblem.H"
#include "baci_particle_rigidbody_datastate.H"

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
PARTICLERIGIDBODY::RigidBodyRuntimeVtpWriter::RigidBodyRuntimeVtpWriter(const Epetra_Comm& comm)
    : comm_(comm), setuptime_(0.0)
{
  // empty constructor
}

void PARTICLERIGIDBODY::RigidBodyRuntimeVtpWriter::Init(
    const std::shared_ptr<PARTICLERIGIDBODY::RigidBodyDataState> rigidbodydatastate)
{
  // set rigid body data state container
  rigidbodydatastate_ = rigidbodydatastate;

  // construct the writer object
  visualization_manager_ = std::make_shared<IO::VisualizationManager>(
      IO::VisualizationParametersFactory(
          DRT::Problem::Instance()->IOParams().sublist("RUNTIME VTK OUTPUT")),
      comm_, "rigidbody");
}

void PARTICLERIGIDBODY::RigidBodyRuntimeVtpWriter::ReadRestart(
    const std::shared_ptr<IO::DiscretizationReader> reader)
{
  // get restart time
  setuptime_ = reader->ReadDouble("time");
}

void PARTICLERIGIDBODY::RigidBodyRuntimeVtpWriter::SetRigidBodyPositionsAndStates(
    const std::vector<int>& ownedrigidbodies)
{
  auto& visualization_data = visualization_manager_->GetVisualizationDataMutable();

  // rigid body position
  {
    // get and prepare storage for position data
    std::vector<double>& posdata = visualization_data.GetPointCoordinatesMutable();
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
    visualization_data.SetPointDataVector<double>("mass", massdata, 1);
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
    visualization_data.SetPointDataVector<double>("velocity", veldata, 3);
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
    visualization_data.SetPointDataVector<double>("acceleration", accdata, 3);
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
    visualization_data.SetPointDataVector<double>("angular velocity", angveldata, 3);
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
    visualization_data.SetPointDataVector<double>("angular acceleration", angaccdata, 3);
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
    visualization_data.SetPointDataVector<double>("force", forcedata, 3);
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
    visualization_data.SetPointDataVector<double>("torque", torquedata, 3);
  }

  // rigid body global id
  {
    // prepare rigid body global id data
    std::vector<double> globaliddata;
    globaliddata.reserve(ownedrigidbodies.size());

    // copy rigid body global id data
    for (int rigidbody_k : ownedrigidbodies) globaliddata.push_back(rigidbody_k);

    // append rigid body global id data to vtp writer
    visualization_data.SetPointDataVector<double>("globalid", globaliddata, 1);
  }

  // rigid body owner
  {
    // set rigid body owner data
    std::vector<double> ownerdata(ownedrigidbodies.size(), comm_.MyPID());

    // append owner of rigid bodies to vtp writer
    visualization_data.SetPointDataVector<double>("owner", ownerdata, 1);
  }
}

void PARTICLERIGIDBODY::RigidBodyRuntimeVtpWriter::WriteToDisk(
    const double time, const unsigned int timestep_number)
{
  visualization_manager_->WriteToDisk(time, timestep_number);
}
