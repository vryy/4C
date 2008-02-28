
#ifdef CCADISCRET

#include "fsi_statustest.H"
#include "../drt_lib/drt_dserror.H"

#include <NOX_Common.H>
#include <NOX_Abstract_Vector.H>
#include <NOX_Abstract_Group.H>
#include <NOX_Solver_Generic.H>
#include <NOX_Utils.H>

#include <NOX_Epetra_Vector.H>

#include <Thyra_DefaultProductVector.hpp>



/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FSI::GenericNormF::GenericNormF(std::string name,
                                 double tolerance,
                                 ScaleType stype)
  : status_(NOX::StatusTest::Unevaluated),
    normType_(NOX::Abstract::Vector::TwoNorm),
    scaleType_(stype),
    toleranceType_(Absolute),
    specifiedTolerance_(tolerance),
    initialTolerance_(1.0),
    trueTolerance_(tolerance),
    normF_(0.0),
    name_(name)
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double FSI::GenericNormF::computeNorm(const Epetra_Vector& v)
{
  int n = v.GlobalLength();
  double norm;
  int err;

  switch (normType_)
  {
  case NOX::Abstract::Vector::TwoNorm:
    err = v.Norm2(&norm);
    if (err!=0)
      dserror("norm failed");
    if (scaleType_ == Scaled)
      norm /= sqrt(1.0 * n);
    break;

  case NOX::Abstract::Vector::OneNorm:
    err = v.Norm1(&norm);
    if (err!=0)
      dserror("norm failed");
    if (scaleType_ == Scaled)
      norm /= n;
    break;

  case NOX::Abstract::Vector::MaxNorm:
    err = v.NormInf(&norm);
    if (err!=0)
      dserror("norm failed");
    if (scaleType_ == Scaled)
      norm /= n;
    break;

  default:
    dserror("norm type confusion");
    break;
  }

  return norm;
}


#if 0
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::GenericNormF::relativeSetup(NOX::Abstract::Group& initialGuess)
{
  NOX::Abstract::Group::ReturnType rtype;
  rtype = initialGuess.computeF();
  if (rtype != NOX::Abstract::Group::Ok)
  {
    utils.err() << "NOX::StatusTest::NormF::NormF - Unable to compute F"
		<< endl;
    throw "NOX Error";
  }

  initialTolerance_ = computeNorm(initialGuess);
  trueTolerance_ = specifiedTolerance_ / initialTolerance_;
}
#endif


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
NOX::StatusTest::StatusType
FSI::GenericNormF::checkStatus(const NOX::Solver::Generic& problem,
                                NOX::StatusTest::CheckType checkType)
{
  if (checkType == NOX::StatusTest::None)
  {
    normF_ = 0.0;
    status_ = NOX::StatusTest::Unevaluated;
  }
  else
  {
    normF_ = computeNorm( problem.getSolutionGroup() );
    if ((normF_ != -1) and (normF_ < trueTolerance_))
    {
      status_ = NOX::StatusTest::Converged;
    }
    else
    {
      status_ = NOX::StatusTest::Unconverged;
    }
  }

  return status_;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
NOX::StatusTest::StatusType FSI::GenericNormF::getStatus() const
{
  return status_;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::ostream& FSI::GenericNormF::print(std::ostream& stream, int indent) const
{
  for (int j = 0; j < indent; j ++)
    stream << ' ';

  stream << status_
         << name_ << "-Norm = " << NOX::Utils::sciformat(normF_,3)
         << " < " << NOX::Utils::sciformat(trueTolerance_, 3)
         << "\n";

  for (int j = 0; j < indent; j ++)
    stream << ' ';

  stream << setw(13) << " (";

  if (scaleType_ == Scaled)
    stream << "Length-Scaled";
  else
    stream << "Unscaled";

  stream << " ";

  if (normType_ == NOX::Abstract::Vector::TwoNorm)
    stream << "Two-Norm";
  else if (normType_ == NOX::Abstract::Vector::OneNorm)
    stream << "One-Norm";
  else if (normType_ == NOX::Abstract::Vector::MaxNorm)
    stream << "Max-Norm";

  stream << ", ";

  if (toleranceType_ == Absolute)
    stream << "Absolute Tolerance";
  else
    stream << "Relative Tolerance";

  stream << ")\n";

  return stream;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double FSI::GenericNormF::getNormF() const
{
  return normF_;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double FSI::GenericNormF::getTrueTolerance() const
{
  return trueTolerance_;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double FSI::GenericNormF::getSpecifiedTolerance() const
{
  return specifiedTolerance_;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double FSI::GenericNormF::getInitialTolerance() const
{
  return initialTolerance_;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FSI::PartialNormF::PartialNormF(std::string name,
                                const Epetra_Map &fullmap,
                                const Epetra_Map &innermap,
                                double tolerance,
                                ScaleType stype)
  : GenericNormF(name,tolerance,stype),
    innermap_(innermap)
{
  extractor_ = Teuchos::rcp(new Epetra_Import(innermap_, fullmap));
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double FSI::PartialNormF::computeNorm(const NOX::Abstract::Group& grp)
{
  if (!grp.isF())
    return -1.0;

  // extract the block epetra vector

  const NOX::Abstract::Vector& abstract_f = grp.getF();
  const NOX::Epetra::Vector& f = Teuchos::dyn_cast<const NOX::Epetra::Vector>(abstract_f);

  // extract the inner vector elements we are interessted in

  Teuchos::RCP<Epetra_Vector> v = Teuchos::rcp(new Epetra_Vector(innermap_));
  int err = v->Import(f.getEpetraVector(), *extractor_, Insert);
  if (err!=0)
    dserror("import failed with err=%d", err);

  return FSI::GenericNormF::computeNorm(*v);
}


#if 0

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FSI::InterfaceNormF::InterfaceNormF(double structfac,
                                     const Epetra_Map &structblockmap,
                                     const Epetra_Map &structinterfacemap,
                                     double fluidfac,
                                     const Epetra_Map &fluidblockmap,
                                     const Epetra_Map &fluidinterfacemap,
                                     const FSI::Coupling& coupsf,
                                     double tolerance,
                                     ScaleType stype)
  : GenericNormF("FSI interface",tolerance,stype),
    structfac_(structfac),
    fluidfac_(fluidfac),
    structblockmap_(structblockmap),
    structinterfacemap_(structinterfacemap),
    fluidblockmap_(fluidblockmap),
    fluidinterfacemap_(fluidinterfacemap),
    coupsf_(coupsf)
{
  structextractor_ = Teuchos::rcp(new Epetra_Import(structinterfacemap_, structblockmap_));
  fluidextractor_ = Teuchos::rcp(new Epetra_Import(fluidinterfacemap_, fluidblockmap_));
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double FSI::InterfaceNormF::computeNorm(const NOX::Abstract::Group& grp)
{
  if (!grp.isF())
    return -1.0;

  // extract the block epetra vector

  const NOX::Abstract::Vector& abstract_res = grp.getF();
  const NOX::Thyra::Vector& thyra_res = Teuchos::dyn_cast<const NOX::Thyra::Vector>(abstract_res);

  const Thyra::DefaultProductVector<double>& res =
    Teuchos::dyn_cast<const Thyra::DefaultProductVector<double> >(*thyra_res.getThyraRCPVector());

  // extract structure interface residual

  Teuchos::RCP<const Thyra::VectorBase<double> > thyra_s = res.getVectorBlock(0);

  // This is ugly. We have to strip the RCP on const VectorBase and
  // create a new (not owning) one, because the original RCP contains
  // the original (non-const) Epetra_Vector and we cannot extract it
  // from a const RCP. :(
  Teuchos::RCP<const Epetra_Vector> epetra_s = Thyra::get_Epetra_Vector(structblockmap_,
                                                                        Teuchos::rcp(&*thyra_s,false));

  Teuchos::RCP<Epetra_Vector> sv = Teuchos::rcp(new Epetra_Vector(structinterfacemap_));
  int err = sv->Import(*epetra_s, *structextractor_, Insert);
  if (err!=0)
    dserror("import failed with err=%d", err);

  // extract fluid interface residual

  Teuchos::RCP<const Thyra::VectorBase<double> > thyra_f = res.getVectorBlock(1);

  // This is ugly. We have to strip the RCP on const VectorBase and
  // create a new (not owning) one, because the original RCP contains
  // the original (non-const) Epetra_Vector and we cannot extract it
  // from a const RCP. :(
  Teuchos::RCP<const Epetra_Vector> epetra_f = Thyra::get_Epetra_Vector(fluidblockmap_,
                                                                        Teuchos::rcp(&*thyra_f,false));

  Teuchos::RCP<Epetra_Vector> fv = Teuchos::rcp(new Epetra_Vector(fluidinterfacemap_));
  err = fv->Import(*epetra_f, *fluidextractor_, Insert);
  if (err!=0)
    dserror("import failed with err=%d", err);

  // transfer the fluid interface vector to the structural side
  fv = coupsf_.SlaveToMaster(fv);

  // add both residual with appropriate scaling
  sv->Update(fluidfac_, *fv, structfac_);

  return FSI::GenericNormF::computeNorm(*sv);
}

#endif

#endif
