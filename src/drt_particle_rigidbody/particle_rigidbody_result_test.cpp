/*---------------------------------------------------------------------------*/
/*! \file
\brief rigid body result test for particle simulations
\level 1
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include "particle_rigidbody_result_test.H"

#include "particle_rigidbody_interface.H"
#include "particle_rigidbody_datastate.H"

#include "../drt_lib/drt_linedefinition.H"
#include "../drt_lib/drt_dserror.H"

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
PARTICLERIGIDBODY::RigidBodyResultTest::RigidBodyResultTest() : DRT::ResultTest("RIGIDBODY")
{
  // empty constructor
}

void PARTICLERIGIDBODY::RigidBodyResultTest::Init()
{
  // nothing to do
}

void PARTICLERIGIDBODY::RigidBodyResultTest::Setup(
    const std::shared_ptr<PARTICLERIGIDBODY::RigidBodyHandlerInterface> particlerigidbodyinterface)
{
  // set interface to rigid body handler
  particlerigidbodyinterface_ = particlerigidbodyinterface;
}

void PARTICLERIGIDBODY::RigidBodyResultTest::TestSpecial(
    DRT::INPUT::LineDefinition& res, int& nerr, int& test_count)
{
  // get owned rigid bodies by this processor
  const std::vector<int>& ownedrigidbodies = particlerigidbodyinterface_->GetOwnedRigidBodies();

  // extract global id of rigid body
  int globalid;
  res.ExtractInt("ID", globalid);

  // rigid body owned by this processor
  if (std::find(ownedrigidbodies.begin(), ownedrigidbodies.end(), globalid) !=
      ownedrigidbodies.end())
  {
    // get rigid body data state container
    std::shared_ptr<PARTICLERIGIDBODY::RigidBodyDataState> rigidbodydatastate =
        particlerigidbodyinterface_->GetRigidBodyDataState();

    // get result
    std::string quantity;
    res.ExtractString("QUANTITY", quantity);

    // init actual result
    double actresult = 0.0;

    // position
    if (quantity == "posx" or quantity == "posy" or quantity == "posz")
    {
      // get reference to rigid body position
      const std::vector<std::vector<double>>& pos = rigidbodydatastate->GetRefPosition();

      // get actual result
      if (quantity == "posx")
        actresult = pos[globalid][0];
      else if (quantity == "posy")
        actresult = pos[globalid][1];
      else if (quantity == "posz")
        actresult = pos[globalid][2];
    }
    // velocity
    else if (quantity == "velx" or quantity == "vely" or quantity == "velz")
    {
      // get reference to rigid body velocity
      const std::vector<std::vector<double>>& vel = rigidbodydatastate->GetRefVelocity();

      // get actual result
      if (quantity == "velx")
        actresult = vel[globalid][0];
      else if (quantity == "vely")
        actresult = vel[globalid][1];
      else if (quantity == "velz")
        actresult = vel[globalid][2];
    }
    // velocity
    else if (quantity == "angvelx" or quantity == "angvely" or quantity == "angvelz")
    {
      // get reference to rigid body angular velocity
      const std::vector<std::vector<double>>& angvel = rigidbodydatastate->GetRefAngularVelocity();

      // get actual result
      if (quantity == "angvelx")
        actresult = angvel[globalid][0];
      else if (quantity == "angvely")
        actresult = angvel[globalid][1];
      else if (quantity == "angvelz")
        actresult = angvel[globalid][2];
    }
    // mass
    else if (quantity == "mass")
    {
      // get reference to rigid body mass
      const std::vector<double>& mass = rigidbodydatastate->GetRefMass();

      // get actual result
      actresult = mass[globalid];
    }
    else
      dserror("result check failed with unknown quantity '%s'!", quantity.c_str());

    // compare values
    const int err = CompareValues(actresult, "SPECIAL", res);
    nerr += err;
    test_count++;
  }
}
