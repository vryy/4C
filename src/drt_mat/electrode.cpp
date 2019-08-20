/*----------------------------------------------------------------------*/
/*! \file
\brief electrode material carrying concentration and electric potential as degrees of freedom

\level 2

\maintainer Christoph Schmidt
*/
/*----------------------------------------------------------------------*/
#include "electrode.H"

#include "../drt_io/io_control.H"

#include "../drt_lib/drt_globalproblem.H"

#include "../drt_mat/matpar_bundle.H"

#include <Epetra_SerialDenseSolver.h>

/*----------------------------------------------------------------------*
 | constructor                                               fang 08/15 |
 *----------------------------------------------------------------------*/
MAT::PAR::Electrode::Electrode(Teuchos::RCP<MAT::PAR::Material> matdata)
    : ElchSingleMat(matdata),
      cmax_(matdata->GetDouble("C_MAX")),
      ocpmodel_(StringToOCPModel(*matdata->Get<std::string>("OCP_MODEL"))),
      ocpparanum_(matdata->GetInt("OCP_PARA_NUM")),
      ocppara_(*matdata->Get<std::vector<double>>("OCP_PARA")),
      X_(0, 0.),
      b_(0, 0.),
      a_(0, 0.),
      m_(0, 0.),
      xmin_(matdata->GetDouble("X_MIN")),
      xmax_(matdata->GetDouble("X_MAX"))
{
  // safety checks
  if (cmax_ < 1.e-12)
    dserror("Saturation value c_max of intercalated Lithium concentration is too small!");
  if ((int)ocppara_.size() != ocpparanum_)
    dserror(
        "Length of coefficient vector for electrode half cell open circuit potential doesn't match "
        "prescribed number of coefficients!");
  if ((xmin_ > 1.0) or (xmax_ > 1.0))
    dserror(
        "Lower bound (X_MIN) and upper bound (X_MAX) of range of validity for ocp calculation "
        "model cannot be larger than one since X "
        "is calculated as c/c_max! If you do not want to prescribe bounds, you have to set the two "
        "variables to negative values. "
        "If you set the bounds to realistic values (i.e. [0,1]) you will get a warning printed to "
        "the screen if bounds are violated throughout the simulation time!");
  if (xmin_ > xmax_) dserror("X_MIN cannot be larger than X_MAX!");

  // additional preparations
  std::string ocpcsv(*matdata->Get<std::string>("OCP_CSV"));
  switch (ocpmodel_)
  {
    case ocp_csv:
    {
      // safety checks
      if (ocpcsv.length() == 0)
        dserror("You forgot to specify the *.csv file for the half cell open circuit potential!");
      if (ocpparanum_)
        dserror(
            "Must not specify any parameters in case half-cell open-circuit potential is to be "
            "determined via a *.csv file!");

      // parse *.csv file
      if (ocpcsv[0] != '/')
      {
        std::string ocpcsvpath = DRT::Problem::Instance()->OutputControlFile()->InputFileName();
        ocpcsvpath = ocpcsvpath.substr(0, ocpcsvpath.rfind('/') + 1);
        ocpcsv.insert(ocpcsv.begin(), ocpcsvpath.begin(), ocpcsvpath.end());
      }
      std::ifstream file(ocpcsv);
      if (!file.good()) dserror("Invalid file!");
      std::string line, value;
      std::vector<double> ocp(0, 0.);
      while (getline(file, line))
      {
        std::istringstream linestream(line);
        std::getline(linestream, value, ',');
        try
        {
          X_.push_back(std::stod(value));
        }
        catch (...)
        {
          continue;
        }

        std::getline(linestream, value);

        try
        {
          ocp.push_back(std::stod(value));
        }
        catch (...)
        {
          dserror("Invalid *.csv file!");
        }
      }
      file.close();

      // safety checks
      if (X_.size() != ocp.size()) dserror("Internal error! Vector lengths have to match!");
      if (X_.size() < 2) dserror("Need at least two data points for cubic spline interpolation!");
      for (unsigned i = 0; i < X_.size() - 1; ++i)
        if (X_[i + 1] <= X_[i]) dserror("Data points must be sorted in ascending order!");

      // build coefficient matrix and right-hand side
      const unsigned N = X_.size() - 2;
      Epetra_SerialDenseMatrix A(N, N);
      Epetra_SerialDenseVector M(N), B(N);
      for (unsigned i = 0; i < N; ++i)
      {
        const double Xm = X_[i + 1] - X_[i], Xp = X_[i + 2] - X_[i + 1], ocpm = ocp[i + 1] - ocp[i],
                     ocpp = ocp[i + 2] - ocp[i + 1];
        if (i > 0) A(i, i - 1) = Xm;
        A(i, i) = 2. * (Xm + Xp);
        if (i < N - 1) A(i, i + 1) = Xp;
        B(i) = -ocpm / Xm + ocpp / Xp;
      }

      // solve for third-order coefficients for cubic spline interpolation
      Epetra_SerialDenseSolver solver;
      solver.SetMatrix(A);
      solver.SetVectors(M, B);
      solver.FactorWithEquilibration(true);
      solver.SolveToRefinedSolution(true);
      if (solver.Factor() or solver.Solve())
        dserror("Solution of linear system of equations failed!");

      // fill coefficient vectors
      m_.resize(X_.size(), 0.);
      for (unsigned i = 1; i < m_.size() - 1; ++i) m_[i] = M(i - 1);
      b_.resize(X_.size() - 1, 0.);
      a_.resize(X_.size() - 1, 0.);
      for (unsigned i = 0; i < b_.size(); ++i)
      {
        const double Xm = X_[i + 1] - X_[i];
        b_[i] = ocp[i] - Xm * Xm * m_[i];
        a_[i] = (ocp[i + 1] - ocp[i]) / Xm - Xm * (m_[i + 1] - m_[i]);
      }

      break;
    }

    case ocp_polynomial:
    case ocp_redlichkister:
    case ocp_taralov:
    {
      // safety checks
      if (ocpparanum_ < 1)
        dserror("No parameters found for electrode half cell open circuit potential!");
      if (ocpcsv.length())
        dserror(
            "Must not specify *.csv file with data points for chosen half cell open circuit "
            "potential model!");
      if (ocpmodel_ == ocp_taralov and ocpparanum_ != 13)
        dserror(
            "Electrode half cell open circuit potential according to Taralov, Taralova, Popov, "
            "Iliev, Latz, and Zausch (2012) needs to be specified by exactly 13 coefficients!");

      break;
    }

    default:
    {
      // safety check
      dserror("Invalid model for half-cell open-circuit potential!");

      break;
    }
  }

  return;
}


/*----------------------------------------------------------------------*
 | create instance of electrode material                     fang 02/15 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<MAT::Material> MAT::PAR::Electrode::CreateMaterial()
{
  return Teuchos::rcp(new MAT::Electrode(this));
}


/*---------------------------------------------------------------------------*
 | convert string to model for half cell open circuit potential   fang 08/15 |
 *---------------------------------------------------------------------------*/
MAT::PAR::OCPModels MAT::PAR::Electrode::StringToOCPModel(const std::string& ocpmodelstring) const
{
  OCPModels ocpmodelenum(ocp_undefined);

  // Redlich-Kister expansion
  if (ocpmodelstring == "Redlich-Kister") ocpmodelenum = ocp_redlichkister;

  // empirical correlation given in Taralov, Taralova, Popov, Iliev, Latz, and Zausch (2012)
  else if (ocpmodelstring == "Taralov")
    ocpmodelenum = ocp_taralov;

  // polynomial
  else if (ocpmodelstring == "Polynomial")
    ocpmodelenum = ocp_polynomial;

  // *.csv file
  else if (ocpmodelstring == "csv")
    ocpmodelenum = ocp_csv;

  // unknown model
  else
    ocpmodelenum = ocp_undefined;

  return ocpmodelenum;
}


MAT::ElectrodeType MAT::ElectrodeType::instance_;


DRT::ParObject* MAT::ElectrodeType::Create(const std::vector<char>& data)
{
  MAT::Electrode* electrode = new MAT::Electrode();
  electrode->Unpack(data);
  return electrode;
}


/*----------------------------------------------------------------------*
 | construct empty electrode material                        fang 02/15 |
 *----------------------------------------------------------------------*/
MAT::Electrode::Electrode() : params_(NULL) { return; }


/*-----------------------------------------------------------------------------*
 | construct electrode material with specific material parameters   fang 02/15 |
 *-----------------------------------------------------------------------------*/
MAT::Electrode::Electrode(MAT::PAR::Electrode* params) : params_(params) { return; }


/*----------------------------------------------------------------------*
 | pack material for communication purposes                  fang 02/15 |
 *----------------------------------------------------------------------*/
void MAT::Electrode::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data, type);

  int matid = -1;
  if (params_ != NULL) matid = params_->Id();  // in case we are in post-process mode
  AddtoPack(data, matid);

  return;
}


/*----------------------------------------------------------------------*
 | unpack data from a char vector                            fang 02/15 |
 *----------------------------------------------------------------------*/
void MAT::Electrode::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  // extract type
  int type = 0;
  ExtractfromPack(position, data, type);
  if (type != UniqueParObjectId()) dserror("Wrong instance type data!");

  // matid and recover params_
  int matid;
  ExtractfromPack(position, data, matid);
  params_ = NULL;
  if (DRT::Problem::Instance()->Materials() != Teuchos::null)
    if (DRT::Problem::Instance()->Materials()->Num() != 0)
    {
      const int probinst = DRT::Problem::Instance()->Materials()->GetReadFromProblem();
      MAT::PAR::Parameter* mat =
          DRT::Problem::Instance(probinst)->Materials()->ParameterById(matid);
      if (mat->Type() == MaterialType())
        params_ = static_cast<MAT::PAR::Electrode*>(mat);
      else
        dserror("Type of parameter material %d does not fit to calling type %d", mat->Type(),
            MaterialType());
    }

  if (position != data.size()) dserror("Mismatch in size of data %d <-> %d", data.size(), position);

  return;
}


/*----------------------------------------------------------------------*
 | compute half cell open circuit potential                  fang 08/15 |
 *----------------------------------------------------------------------*/
double MAT::Electrode::ComputeOpenCircuitPotential(const double concentration,  //!< concentration
    const double faraday,  //!< Faraday constant
    const double frt       //!< factor F/RT
    ) const
{
  double ocp(0.);

  // intercalation fraction
  const double X = concentration / params_->cmax_;

  // print warning to screen if prescribed interval of validity for ocp calculation model is given
  // but not satisfied
  if (((X < params_->xmin_) or (X > params_->xmax_)) and !(params_->xmax_ < 0.))
  {
    std::cout << "WARNING: intercalation fraction X = c/c_max is violating prescribed bounds of "
                 "ocp calculation model. Calculated "
                 "values might therefore not be reasonable!"
              << std::endl;
    std::cout << "X: " << X << " lower bound is: " << params_->xmin_
              << "  upper bound is: " << params_->xmax_ << std::endl
              << std::endl;
  }

  // physically reasonable intercalation fraction
  if (X > 0. and X < 1.)
  {
    switch (params_->ocpmodel_)
    {
      // half cell open circuit potential obtained from cubic spline interpolation of *.csv data
      // points
      case MAT::PAR::ocp_csv:
      {
        // safety check
        if (X < params_->X_.front() or X > params_->X_.back())
          dserror("Intercalation fraction X = %lf lies outside sampling point range!", X);

        // evaluate cubic spline interpolation
        for (unsigned i = 0; i < params_->m_.size() - 1; ++i)
        {
          if (X <= params_->X_[i + 1])
          {
            const double invdX = 1. / (params_->X_[i + 1] - params_->X_[i]);
            ocp = params_->m_[i] * invdX * pow(params_->X_[i + 1] - X, 3) +
                  params_->m_[i + 1] * invdX * pow(X - params_->X_[i], 3) +
                  params_->a_[i] * (X - params_->X_[i]) + params_->b_[i];
            break;
          }
        }

        break;
      }

      // half cell open circuit potential according to Redlich-Kister expansion
      case MAT::PAR::ocp_redlichkister:
      {
        // cf. Colclasure and Kee, Electrochimica Acta 55 (2010) 8960:
        // ocppara_[0]         = DeltaG
        // ocppara_[1,2,3,...] = Redlich-Kister coefficients

        // terms not associated with any Redlich-Kister coefficient
        ocp = params_->ocppara_[0] + faraday / frt * log((1. - X) / X);

        // terms associated with first and second Redlich-Kister coefficients
        // these two terms are separated from the remaining sum and simplified thereafter to remove
        // singularities in the expansion in case X == 0.5
        ocp += params_->ocppara_[1] * (2. * X - 1.) +
               params_->ocppara_[2] * (6. * X * X - 6. * X + 1.);

        // terms associated with remaining Redlich-Kister coefficients
        for (int i = 2; i < params_->ocpparanum_ - 1; ++i)
          ocp += params_->ocppara_[i + 1] *
                 (pow(2. * X - 1., i + 1) - 2. * i * X * (1. - X) * pow(2. * X - 1., i - 1));

        // final scaling
        ocp /= faraday;

        break;
      }

      case MAT::PAR::ocp_taralov:
      {
        // cf. Taralov, Taralova, Popov, Iliev, Latz, and Zausch (2012)
        ocp = params_->ocppara_[0] +
              params_->ocppara_[1] * tanh(params_->ocppara_[2] * X + params_->ocppara_[3]) +
              params_->ocppara_[4] * exp(params_->ocppara_[5] * pow(X, 8.0)) +
              params_->ocppara_[6] * (1 / (pow((params_->ocppara_[7] - X), params_->ocppara_[8])) +
                                         params_->ocppara_[9]) +
              params_->ocppara_[10] * exp(params_->ocppara_[11] * (X + params_->ocppara_[12]));

        break;
      }

      // polynomial ocp
      case MAT::PAR::ocp_polynomial:
      {
        // add constant
        ocp = params_->ocppara_[0];

        // add higher polynomial order terms
        for (int i = 1; i < params_->ocpparanum_; ++i) ocp += params_->ocppara_[i] * pow(X, i);

        break;
      }

      default:
      {
        dserror("Model for half cell open circuit potential not recognized!");
        break;
      }
    }
  }

  // non-physical intercalation fraction
  else
    ocp = std::numeric_limits<double>::infinity();

  return ocp;
}  // MAT::Electrode::ComputeOpenCircuitPotential


/*---------------------------------------------------------------------------------------------------------*
 | compute first derivative of half cell open circuit potential with respect to concentration   fang
 08/15 |
 *---------------------------------------------------------------------------------------------------------*/
double MAT::Electrode::ComputeFirstDerivOpenCircuitPotential(
    const double concentration,  //!< concentration
    const double faraday,        //!< Faraday constant
    const double frt             //!< factor F/RT
    ) const
{
  double ocpderiv(0.);

  // intercalation fraction
  const double X = concentration / params_->cmax_;

  // physically reasonable intercalation fraction
  if (X > 0. and X < 1.)
  {
    switch (params_->ocpmodel_)
    {
      // derivative of half cell open circuit potential w.r.t. concentration, obtained from cubic
      // spline interpolation of *.csv data points
      case MAT::PAR::ocp_csv:
      {
        // safety check
        if (X < params_->X_.front() or X > params_->X_.back())
          dserror("Intercalation fraction X = %lf lies outside sampling point range!", X);

        // evaluate derivative of cubic spline interpolation w.r.t. concentration
        for (unsigned i = 0; i < params_->m_.size() - 1; ++i)
        {
          if (X <= params_->X_[i + 1])
          {
            const double invdX = 1. / (params_->X_[i + 1] - params_->X_[i]);
            ocpderiv = -3. * params_->m_[i] * invdX * pow(params_->X_[i + 1] - X, 2) +
                       3. * params_->m_[i + 1] * invdX * pow(X - params_->X_[i], 2) +
                       params_->a_[i];
            break;
          }
        }

        break;
      }

      // derivative of half cell open circuit potential w.r.t. concentration according to
      // Redlich-Kister expansion
      case MAT::PAR::ocp_redlichkister:
      {
        // cf. Colclasure and Kee, Electrochimica Acta 55 (2010) 8960:
        // ocppara_[0]         = DeltaG
        // ocppara_[1,2,3,...] = Redlich-Kister coefficients

        // term not associated with any Redlich-Kister coefficient
        ocpderiv = faraday / (2. * frt * X * (X - 1.));

        // terms associated with first, second, and third Redlich-Kister coefficients
        // these three terms are separated from the remaining sum and simplified thereafter to
        // remove singularities in the derivative of the expansion in case X == 0.5
        ocpderiv += params_->ocppara_[1] + params_->ocppara_[2] * (6. * X - 3.) +
                    params_->ocppara_[3] * (24. * X * X - 24. * X + 5.);

        // terms associated with remaining Redlich-Kister coefficients
        for (int i = 3; i < params_->ocpparanum_ - 1; ++i)
          ocpderiv += params_->ocppara_[i + 1] *
                      ((2. * i + 1.) * pow(2. * X - 1., i) +
                          2. * X * i * (X - 1.) * (i - 1.) * pow(2. * X - 1., i - 2));

        // intermediate scaling
        ocpderiv *= 2. / faraday;

        break;
      }

      // derivative of half cell open circuit potential w.r.t. concentration according to Taralov,
      // Taralova, Popov, Iliev, Latz, and Zausch (2012)
      case MAT::PAR::ocp_taralov:
      {
        ocpderiv = params_->ocppara_[1] * params_->ocppara_[2] /
                       pow(cosh(params_->ocppara_[2] * X + params_->ocppara_[3]), 2) +
                   8. * params_->ocppara_[4] * params_->ocppara_[5] *
                       exp(params_->ocppara_[5] * pow(X, 8)) * pow(X, 7) +
                   params_->ocppara_[6] * params_->ocppara_[8] /
                       pow(params_->ocppara_[7] - X, params_->ocppara_[8] + 1.) +
                   params_->ocppara_[10] * params_->ocppara_[11] *
                       exp(params_->ocppara_[11] * (X + params_->ocppara_[12]));

        break;
      }

      // derivative of polynomial half cell open circuit potential w.r.t. concentration
      case MAT::PAR::ocp_polynomial:
      {
        for (int i = 1; i < params_->ocpparanum_; ++i)
          ocpderiv += i * params_->ocppara_[i] * pow(X, i - 1);

        break;
      }

      default:
      {
        dserror("Model for half cell open circuit potential not recognized!");
        break;
      }
    }

    // final scaling
    ocpderiv /= params_->cmax_;
  }

  // non-physical intercalation fraction
  else
    ocpderiv = std::numeric_limits<double>::infinity();

  return ocpderiv;
}  // MAT::Electrode::ComputeFirstDerivOpenCircuitPotential


/*----------------------------------------------------------------------------------------------------------*
 | compute second derivative of half cell open circuit potential with respect to concentration fang
 08/15 |
 *----------------------------------------------------------------------------------------------------------*/
double MAT::Electrode::ComputeSecondDerivOpenCircuitPotential(
    const double concentration,  //!< concentration
    const double faraday,        //!< Faraday constant
    const double frt             //!< factor F/RT
    ) const
{
  double ocpderiv2(0.);

  // intercalation fraction
  const double X = concentration / params_->cmax_;

  // physically reasonable intercalation fraction
  if (X > 0. and X < 1.)
  {
    switch (params_->ocpmodel_)
    {
      // second derivative of half cell open circuit potential w.r.t. concentration, obtained from
      // cubic spline interpolation of *.csv data points
      case MAT::PAR::ocp_csv:
      {
        // safety check
        if (X < params_->X_.front() or X > params_->X_.back())
          dserror("Intercalation fraction X = %lf lies outside sampling point range!", X);

        // evaluate second derivative of cubic spline interpolation w.r.t. concentration
        for (unsigned i = 0; i < params_->m_.size() - 1; ++i)
        {
          if (X <= params_->X_[i + 1])
          {
            const double invdX = 1. / (params_->X_[i + 1] - params_->X_[i]);
            ocpderiv2 = 6. * params_->m_[i] * invdX * (params_->X_[i + 1] - X) +
                        6. * params_->m_[i + 1] * invdX * (X - params_->X_[i]);
            break;
          }
        }

        break;
      }

      // second derivative of half cell open circuit potential w.r.t. concentration according to
      // Redlich-Kister expansion
      case MAT::PAR::ocp_redlichkister:
      {
        // cf. Colclasure and Kee, Electrochimica Acta 55 (2010) 8960:
        // ocppara_[0]         = DeltaG
        // ocppara_[1,2,3,...] = Redlich-Kister coefficients

        // term not associated with any Redlich-Kister coefficient
        ocpderiv2 = -faraday * (2. * X - 1.) / (4. * frt * X * X * (X - 1.) * (X - 1.));

        // term associated with first Redlich-Kister coefficient vanishes

        // terms associated with second, third, and fourth Redlich-Kister coefficients
        // these three terms are separated from the remaining sum and simplified thereafter to
        // remove singularities in the second derivative of the expansion in case X == 0.5
        ocpderiv2 += 3. * params_->ocppara_[2] + params_->ocppara_[3] * (24. * X - 12.) +
                     params_->ocppara_[4] * (120. * X * X - 120. * X + 27.);

        // terms associated with remaining Redlich-Kister coefficients
        for (int i = 4; i < params_->ocpparanum_ - 1; ++i)
          ocpderiv2 += params_->ocppara_[i + 1] *
                       (3. * i * i * pow(2. * X - 1., i - 1) +
                           2. * i * (i - 1.) * (i - 2.) * X * (X - 1.) * pow(2. * X - 1., i - 3));

        // intermediate scaling
        ocpderiv2 *= 4. / faraday;

        break;
      }

      // second derivative of half cell open circuit potential w.r.t. concentration according to
      // Taralov, Taralova, Popov, Iliev, Latz, and Zausch (2012)
      case MAT::PAR::ocp_taralov:
      {
        ocpderiv2 = -2. * params_->ocppara_[1] * pow(params_->ocppara_[2], 2) /
                        pow(cosh(params_->ocppara_[2] * X + params_->ocppara_[3]), 2) *
                        tanh(params_->ocppara_[2] * X + params_->ocppara_[3]) +
                    8. * params_->ocppara_[4] * params_->ocppara_[5] * pow(X, 6) *
                        exp(params_->ocppara_[5] * pow(X, 8)) *
                        (7. + 8. * params_->ocppara_[5] * pow(X, 8)) +
                    params_->ocppara_[6] * params_->ocppara_[8] * (params_->ocppara_[8] + 1.) /
                        pow(params_->ocppara_[7] - X, params_->ocppara_[8] + 2.) +
                    params_->ocppara_[10] * pow(params_->ocppara_[11], 2) *
                        exp(params_->ocppara_[11] * (X + params_->ocppara_[12]));

        break;
      }

      // second derivative of polynomial half cell open circuit potential w.r.t. concentration
      case MAT::PAR::ocp_polynomial:
      {
        for (int i = 2; i < params_->ocpparanum_; ++i)
          ocpderiv2 += i * (i - 1) * params_->ocppara_[i] * pow(X, i - 2);

        break;
      }

      default:
      {
        dserror("Model for half cell open circuit potential not recognized!");
        break;
      }
    }

    // final scaling
    ocpderiv2 /= params_->cmax_ * params_->cmax_;
  }

  // non-physical intercalation fraction
  else
    ocpderiv2 = std::numeric_limits<double>::infinity();

  return ocpderiv2;
}  // MAT::Electrode::ComputeSecondDerivOpenCircuitPotential
