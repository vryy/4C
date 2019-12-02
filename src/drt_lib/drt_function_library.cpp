/*----------------------------------------------------------------------*/
/*! \file

\brief Collection of functions that are not problem-specific

The functions in this file are not problem-specific and may be useful for a number of applications.

\maintainer Sebastian Proell

\level 3

*/
/*----------------------------------------------------------------------*/

#include <Sacado.hpp>
#include "drt_function.H"
#include "drt_function_library.H"
#include "drt_globalproblem.H"



/*----------------------------------------------------------------------* |
  Constructor of ControlledRotation                          hahn 04 / 13 |
 *----------------------------------------------------------------------*/
DRT::UTILS::ControlledRotationFunction::ControlledRotationFunction(
    std::string fileName, std::string type, double origin_x, double origin_y, double origin_z)
    : Function(), NUMMANEUVERCELLS_(4)
{
  // Initialize condition type
  if (type == "STRUCTURE")
    physicsType_ = PhysicalField::Structure;
  else if (type == "FLUID")
    physicsType_ = PhysicalField::Fluid;
  else
    dserror(
        "When using the function CONTROLLEDROTATION, the type must be either 'STRUCTURE' or "
        "'FLUID'");

  // Initialize origin, about which the rotation is performed
  origin_(0, 0) = origin_x;
  origin_(1, 0) = origin_y;
  origin_(2, 0) = origin_z;

  // Initialize time of previous time step (at t-deltaT)
  timeOld_ = 0.0;

  // Initialize previous angular acceleration (at t-deltaT)
  omegaDotOld_B_(0, 0) = 0.0;
  omegaDotOld_B_(1, 0) = 0.0;
  omegaDotOld_B_(2, 0) = 0.0;

  // Current angular rate (at t)
  omega_B_(0, 0) = 0.0;
  omega_B_(1, 0) = 0.0;
  omega_B_(2, 0) = 0.0;

  // Initialize satellite's current attitude trafo matrix from B- to I-system (at t)
  satAtt_dcm_IB_(0, 0) = 1.0;
  satAtt_dcm_IB_(0, 1) = 0.0;
  satAtt_dcm_IB_(0, 2) = 0.0;
  satAtt_dcm_IB_(1, 0) = 0.0;
  satAtt_dcm_IB_(1, 1) = 1.0;
  satAtt_dcm_IB_(1, 2) = 0.0;
  satAtt_dcm_IB_(2, 0) = 0.0;
  satAtt_dcm_IB_(2, 1) = 0.0;
  satAtt_dcm_IB_(2, 2) = 1.0;

  // Initialize satellite's current attitude quaternion from B- to I-system (at t)
  satAtt_q_IB_(0, 0) = 0.0;
  satAtt_q_IB_(1, 0) = 0.0;
  satAtt_q_IB_(2, 0) = 0.0;
  satAtt_q_IB_(3, 0) = 1.0;

  // Initialize maneuvers
  maneuvers_.clear();

  // Read maneuver file and fill maneuvers variable
  // *****************************************************************************

  std::string line;
  std::stringstream lineStream;
  std::string cell;

  // Open file
  std::ifstream file(fileName.c_str());
  if (!file.is_open())
  {
    dserror("Unable to open file: %s", fileName.c_str());
  }

  // Loop through all lines
  while (getline(file, line))
  {
    if (!line.empty())
    {
      // Clear local variables
      lineStream.clear();
      cell.clear();

      // Obtain all numManeuverCells=4 values from current line (t, omegaDot_x_B, omegaDot_y_B,
      // omegaDot_z_B)
      lineStream << line;
      for (int i = 0; i < NUMMANEUVERCELLS_; i++)
      {
        // Obtain i-th cell from current line
        getline(lineStream, cell, ' ');

        // If empty cell, than either empty line or one cell in the line
        // missing, anyhow an error.
        if (cell.empty())
        {
          dserror("Error during reading of file: %s", fileName.c_str());
        }

        // Convert input cell from string to double
        double cellDouble = (double)strtod(cell.c_str(), NULL);

        // Add cell to maneuvers vector
        maneuvers_.push_back(cellDouble);
      }
    }
  }

  // Close file
  file.close();

  // Determine number of maneuvers
  numManeuvers_ = (int)(maneuvers_.size()) / (int)NUMMANEUVERCELLS_;

  // Output maneuver list
  printf("\n=================================================================================\n");
  printf("ControlledRotation - %s\n", type.c_str());
  printf("---------------------------------------------------------------------------------\n");
  printf("The following %d maneuvers have been loaded from file %s:\n", numManeuvers_,
      fileName.c_str());
  for (int i = 0; i < numManeuvers_; i++)
  {
    printf("Time: %e  OmegaDot_B: %e, %e, %e \n", maneuvers_[i * NUMMANEUVERCELLS_ + 0],
        maneuvers_[i * NUMMANEUVERCELLS_ + 1], maneuvers_[i * NUMMANEUVERCELLS_ + 2],
        maneuvers_[i * NUMMANEUVERCELLS_ + 3]);
  }
  printf("=================================================================================\n\n");
}

/*----------------------------------------------------------------------*
 | Evaluate ControlledRotation and return for structures the current    |
 | displacement and for fluids the current velocity of the respective   |
 | node for the given index.                                 hahn 04/13 |
 *----------------------------------------------------------------------*/
double DRT::UTILS::ControlledRotationFunction::Evaluate(const int index, const double* xp, double t)
{
  // Check, if a restart has been performed
  // *****************************************************************************
  const int step = DRT::Problem::Instance()->Restart();
  if ((step > 0) && (timeOld_ == 0.0))
  {
    dserror(
        "When using the function CONTROLLEDROTATION, the restart functionality cannot be used!");
  }

  // If new time step, apply angular acceleration (if desired) and determine
  // new attitude satAtt_dcm_IB_
  // *****************************************************************************
  // Determine time difference
  double deltaT = t - timeOld_;

  if (deltaT > 1e-12)
  {  // new time step

    // Determine current angular acceleration (at t)
    // -----------------------------------------------------------------------------
    LINALG::Matrix<3, 1> omegaDot_B;
    omegaDot_B(0, 0) = 0.0;
    omegaDot_B(1, 0) = 0.0;
    omegaDot_B(2, 0) = 0.0;

    for (int i = 0; i < numManeuvers_; i++)
    {
      if (t >= maneuvers_[i * NUMMANEUVERCELLS_ + 0])
      {
        omegaDot_B(0, 0) = maneuvers_[i * NUMMANEUVERCELLS_ + 1];
        omegaDot_B(1, 0) = maneuvers_[i * NUMMANEUVERCELLS_ + 2];
        omegaDot_B(2, 0) = maneuvers_[i * NUMMANEUVERCELLS_ + 3];
      }
    }

    // Calculate current angular rate (at t) by integration
    // of angular acceleration (trapezoidal rule):
    // omega_(t) = deltaT * (omegaDotOld + omegaDot) / 2 + omega_(t-deltaT)
    // -----------------------------------------------------------------------------
    LINALG::Matrix<3, 1> deltaOmega;
    deltaOmega.Update(omegaDotOld_B_, omegaDot_B);  // 1) deltaOmega <- omegaDotOld_ + omegaDot
    deltaOmega.Scale(deltaT / 2.0);                 // 2) deltaOmega <- deltaOmega * deltaT / 2.0
    omega_B_ += deltaOmega;                         // 3) omega_ <- omega_ + deltaOmega

    /* // Debugging output
       cout << "omegaDot: "; omegaDot_B.Print(cout); // Print omegaDot_B
       cout << "omega:    "; omega_B_.Print(cout);   // Print omega_B_
    */

    omegaDotOld_B_ = omegaDot_B;  // Set omegaDotOld_B_ for next time step

    // Calculate new attitude quaternion satAtt_q_IB_ [Wertz, p. 511f]
    // -----------------------------------------------------------------------------
    LINALG::Matrix<4, 4> mOmega;  // Skew-symmetric matrix containing angular velocity components
    mOmega(0, 0) = 0.0;
    mOmega(0, 1) = omega_B_(2, 0);
    mOmega(0, 2) = -omega_B_(1, 0);
    mOmega(0, 3) = omega_B_(0, 0);
    mOmega(1, 0) = -omega_B_(2, 0);
    mOmega(1, 1) = 0.0;
    mOmega(1, 2) = omega_B_(0, 0);
    mOmega(1, 3) = omega_B_(1, 0);
    mOmega(2, 0) = omega_B_(1, 0);
    mOmega(2, 1) = -omega_B_(0, 0);
    mOmega(2, 2) = 0.0;
    mOmega(2, 3) = omega_B_(2, 0);
    mOmega(3, 0) = -omega_B_(0, 0);
    mOmega(3, 1) = -omega_B_(1, 0);
    mOmega(3, 2) = -omega_B_(2, 0);
    mOmega(3, 3) = 0.0;

    mOmega.Scale(deltaT / 2.0);
    mOmega(0, 0) = 1.0;
    mOmega(1, 1) = 1.0;
    mOmega(2, 2) = 1.0;
    mOmega(3, 3) = 1.0;

    LINALG::Matrix<4, 1> satAtt_q_IB_TMP(true);
    satAtt_q_IB_TMP.Multiply(mOmega, satAtt_q_IB_);
    satAtt_q_IB_ = satAtt_q_IB_TMP;

    satAtt_q_IB_.Scale(1 / satAtt_q_IB_.Norm2());  // Normalize attitude quaternion

    // Create transformation matrix satAtt_dcm_IB_ [Wertz, (E-8)]
    // -----------------------------------------------------------------------------
    const double q1 = satAtt_q_IB_(0, 0);
    const double q2 = satAtt_q_IB_(1, 0);
    const double q3 = satAtt_q_IB_(2, 0);
    const double q4 = satAtt_q_IB_(3, 0);

    satAtt_dcm_IB_(0, 0) = q1 * q1 - q2 * q2 - q3 * q3 + q4 * q4;
    satAtt_dcm_IB_(0, 1) = 2.0 * (q1 * q2 - q3 * q4);
    satAtt_dcm_IB_(0, 2) = 2.0 * (q1 * q3 + q2 * q4);
    satAtt_dcm_IB_(1, 0) = 2.0 * (q1 * q2 + q3 * q4);
    satAtt_dcm_IB_(1, 1) = -q1 * q1 + q2 * q2 - q3 * q3 + q4 * q4;
    satAtt_dcm_IB_(1, 2) = 2.0 * (q2 * q3 - q1 * q4);
    satAtt_dcm_IB_(2, 0) = 2.0 * (q1 * q3 - q2 * q4);
    satAtt_dcm_IB_(2, 1) = 2.0 * (q2 * q3 + q1 * q4);
    satAtt_dcm_IB_(2, 2) = -q1 * q1 - q2 * q2 + q3 * q3 + q4 * q4;

    // Update time of last time step
    // -----------------------------------------------------------------------------
    timeOld_ = t;
  }

  // Obtain the current node position in the inertial system
  // *****************************************************************************
  // NOTE: 1) Here it is assumed that the difference between the mesh system and the
  //          body system is just given by a constant displacement named origin_.
  //          Hence the variable origin_ specifies the point, given in the mesh
  //          system, about which is being rotated.
  //       2) The inertial system used here has position and attitude of the body
  //          system at initial time. Thus it is deviating from the ECI system J2000
  //          at most by a constant displacement that is due to the satellite position
  //          at initial time and a constant rotation that is due to the satellite's
  //          attitude at initial time. Due to Galilei-invariance this shouldn't be
  //          a problem and it can easily be converted to the J2000 system.

  // Node reference position, given in the body system
  LINALG::Matrix<3, 1> nodeReferencePos_B;
  nodeReferencePos_B(0, 0) = xp[0] - origin_(0, 0);
  nodeReferencePos_B(1, 0) = xp[1] - origin_(1, 0);
  nodeReferencePos_B(2, 0) = xp[2] - origin_(2, 0);

  // Node position, given in the inertial system
  LINALG::Matrix<3, 1> nodePos_I;
  nodePos_I.Multiply(satAtt_dcm_IB_, nodeReferencePos_B);

  // Calculate and return the displacement/velocity of the node for the given index
  // *****************************************************************************

  if (t <= 0)
  {
    // Return zero displacement/translational velocity of the node for the given index
    return 0.0;
  }

  // Displacement of the node for the given index
  const double dispNodePos_I = nodePos_I(index, 0) - nodeReferencePos_B(index, 0);

  switch (physicsType_)
  {
    case Structure:
    {
      // Return the displacement of the node for the given index
      return dispNodePos_I;
    }
    case Fluid:
    {
      // Node velocity, given in the inertial system: v = omega x r
      double nodeVel_I[3];
      nodeVel_I[0] = omega_B_(1, 0) * nodePos_I(2, 0) - omega_B_(2, 0) * nodePos_I(1, 0);
      nodeVel_I[1] = omega_B_(2, 0) * nodePos_I(0, 0) - omega_B_(0, 0) * nodePos_I(2, 0);
      nodeVel_I[2] = omega_B_(0, 0) * nodePos_I(1, 0) - omega_B_(1, 0) * nodePos_I(0, 0);

      // Return the translational velocity of the node for the given index
      return (nodeVel_I[index]);
    }
    default:
    {
      dserror(
          "When using the function CONTROLLEDROTATION, the type must be either 'STRUCTURE' or "
          "'FLUID'");
      return 0.0;
    }
  }
}

/*----------------------------------------------------------------------*
 | Constructor of AccelerationProfile                        hahn 09/13 |
 *----------------------------------------------------------------------*/
DRT::UTILS::AccelerationProfileFunction::AccelerationProfileFunction(std::string fileName)
    : Function(), NUMACCELERATIONCELLS_(4)
{
  // Initialize variables
  // *****************************************************************************

  // Initialize time of previous time step (at t-deltaT)
  timeOld_ = 0.0;

  // Initialize accelerations
  accelerations_.clear();

  // Initialize current acceleration (at t)
  acc_B_(0, 0) = 0.0;
  acc_B_(1, 0) = 0.0;
  acc_B_(2, 0) = 0.0;

  // Read acceleration profile file and fill acceleration variable
  // *****************************************************************************

  std::string line;
  std::stringstream lineStream;
  std::string cell;

  // Open file
  std::ifstream file(fileName.c_str());
  if (!file.is_open())
  {
    dserror("Unable to open file: %s", fileName.c_str());
  }

  // Loop through all lines
  while (getline(file, line))
  {
    if (!line.empty())
    {
      // Clear local variables
      lineStream.clear();
      cell.clear();

      // Obtain all numAccelerationCells=4 values from current line (t, acc_x_B, acc_y_B, acc_z_B)
      lineStream << line;
      for (int i = 0; i < NUMACCELERATIONCELLS_; i++)
      {
        // Obtain i-th cell from current line
        getline(lineStream, cell, ' ');

        // If empty cell, than either empty line or one cell in the line
        // missing, anyhow an error.
        if (cell.empty())
        {
          dserror("Error during reading of file: %s", fileName.c_str());
        }

        // Convert input cell from string to double
        double cellDouble = (double)strtod(cell.c_str(), NULL);

        // Add cell to acceleration vector
        accelerations_.push_back(cellDouble);
      }
    }
  }

  // Close file
  file.close();

  // Determine number of accelerations
  numAccelerations_ = (int)(accelerations_.size()) / (int)NUMACCELERATIONCELLS_;

  // Output acceleration list
  printf("\n=================================================================================\n");
  printf("AccelerationProfile\n");
  printf("---------------------------------------\n");
  printf("The following %d acceleration rows have been loaded from file %s:\n", numAccelerations_,
      fileName.c_str());
  for (int i = 0; i < numAccelerations_; i++)
  {
    printf("Time: %e  acc_B: %e, %e, %e \n", accelerations_[i * NUMACCELERATIONCELLS_ + 0],
        accelerations_[i * NUMACCELERATIONCELLS_ + 1],
        accelerations_[i * NUMACCELERATIONCELLS_ + 2],
        accelerations_[i * NUMACCELERATIONCELLS_ + 3]);
  }
  printf("=================================================================================\n\n");
}

/*---------------------------------------------------------------------*
 | Evaluate AccelerationProfile and return the respective acceleration |
 | of the node for the given index.                         hahn 09/13 |
 *---------------------------------------------------------------------*/
double DRT::UTILS::AccelerationProfileFunction::Evaluate(
    const int index, const double* xp, double t)
{
  // Determine time difference
  double deltaT = t - timeOld_;

  // If new time step, determine current acceleration
  // *****************************************************************************
  if (deltaT > 1e-9)
  {  // new time step

    // Determine current acceleration (at t)
    // -------------------------------------------------------------------------
    acc_B_(0, 0) = 0.0;
    acc_B_(1, 0) = 0.0;
    acc_B_(2, 0) = 0.0;

    for (int i = 0; i < numAccelerations_; i++)
    {
      if (t >= accelerations_[i * NUMACCELERATIONCELLS_ + 0])
      {
        acc_B_(0, 0) = accelerations_[i * NUMACCELERATIONCELLS_ + 1];
        acc_B_(1, 0) = accelerations_[i * NUMACCELERATIONCELLS_ + 2];
        acc_B_(2, 0) = accelerations_[i * NUMACCELERATIONCELLS_ + 3];
      }
    }

    // Update time of last time step
    // -------------------------------------------------------------------------
    timeOld_ = t;
  }

  // Return the acceleration of the node for the given index
  // *****************************************************************************
  return acc_B_(index, 0);
}

DRT::UTILS::FastPolynomialFunction::FastPolynomialFunction(std::vector<double>* coefficients)
    : Function(), mypoly_(new Polynomial(*coefficients))
{
}

/*---------------------------------------------------------------------*
 | Evaluate                                                proell 01/19 |
 *---------------------------------------------------------------------*/
double DRT::UTILS::FastPolynomialFunction::Evaluate(const double argument) const
{
  return mypoly_->Evaluate(argument);
}

/*---------------------------------------------------------------------*
 | Evaluate derivative                                    proell 01/19 |
 *---------------------------------------------------------------------*/
double DRT::UTILS::FastPolynomialFunction::EvaluateDerivative(const double argument) const
{
  LINALG::Matrix<2, 1> derivs(false);
  mypoly_->Evaluate(argument, derivs);
  return derivs(1);
}

DRT::UTILS::TranslatedFunction::TranslatedFunction(
    Teuchos::RCP<Function> origin, Teuchos::RCP<Function> local)
{
  if (origin->NumberComponents() != nsd_originTranslation)
    dserror("Origin function needs to have exactly %d components but %d were given.",
        nsd_originTranslation, origin->NumberComponents());
  originFunction_ = origin;
  localFunction_ = local;
}

double DRT::UTILS::TranslatedFunction::Evaluate(const int index, const double* x, double t)
{
  if (index < 0 or index >= localFunction_->NumberComponents())
    dserror("Index must be between 0 and %d but is %d.", localFunction_->NumberComponents(), index);

  double new_coord[3];
  for (int i = 0; i < nsd_originTranslation; i++)
  {
    new_coord[i] = x[i] - originFunction_->Evaluate(i, x, t);
  }
  return localFunction_->Evaluate(index, new_coord, t);
}

std::vector<double> DRT::UTILS::TranslatedFunction::EvaluateTimeDerivative(
    const int index, const double* x, const double t, const unsigned deg)
{
  if (deg == 0) return std::vector<double>{Evaluate(index, x, t)};

  if (deg != 1) dserror("Time derivative only implemented for degree <= 1");

  if (index < 0 or index >= localFunction_->NumberComponents())
    dserror("Index must be between 0 and %d but is %d.", localFunction_->NumberComponents(), index);

  std::vector<double> result(deg + 1);
  result[0] = Evaluate(index, x, t);

  double translatedCoord[3];
  double translatedDeriv[3];
  for (int i = 0; i < nsd_originTranslation; i++)
  {
    auto evalResult = originFunction_->EvaluateTimeDerivative(i, x, t, 1);
    translatedCoord[i] = x[i] - evalResult[0];
    translatedDeriv[i] = -evalResult[1];
  }

  auto localSpatialDeriv = localFunction_->EvaluateSpatialDerivative(index, translatedCoord, t);
  auto localValues = localFunction_->EvaluateTimeDerivative(index, translatedCoord, t, 1);
  // total time derivative according to chain rule
  // dh/dt = -df/dx*dg/dt + df/dt
  result[1] = localSpatialDeriv[0] * translatedDeriv[0] +
              localSpatialDeriv[1] * translatedDeriv[1] +
              localSpatialDeriv[2] * translatedDeriv[2] + localValues[1];
  return result;
}
