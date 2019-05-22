/*----------------------------------------------------------------------------*/
/*!
\brief Output for the inverse analysis control routine

\level 3

\maintainer Sebastian Brandstaeter
*/
/*----------------------------------------------------------------------------*/
/* headers */
#include "invana_writer.H"

// baci
#include "../drt_io/io.H"
#include "../drt_lib/drt_dserror.H"

/*----------------------------------------------------------------------------*/
INVANA::InvanaWriter::InvanaWriter() : output_(Teuchos::null) { return; }

/*----------------------------------------------------------------------------*/
void INVANA::InvanaWriter::Init(Teuchos::RCP<IO::DiscretizationWriter> writer)
{
  output_ = writer;

  return;
}

/*----------------------------------------------------------------------------*/
void INVANA::InvanaWriter::WriteNewStep(const int step, const double time)
{
  output_->NewStep(step, time);
  return;
}

/*----------------------------------------------------------------------------*/
void INVANA::InvanaWriter::WriteNamedInt(const std::string name, const int aint)
{
  output_->WriteInt(name, aint);
  return;
}

/*----------------------------------------------------------------------------*/
void INVANA::InvanaWriter::WriteNamedDouble(const std::string name, const double adouble)
{
  output_->WriteDouble(name, adouble);
  return;
}

/*----------------------------------------------------------------------------*/
void INVANA::InvanaWriter::WriteNamedVector(const std::string name,
    Teuchos::RCP<const Epetra_MultiVector> vector, IO::VectorType type) const
{
  output_->WriteVector(name, vector, type);
  return;
}

/*----------------------------------------------------------------------------*/
void INVANA::InvanaWriter::WriteNamedVectors(const std::string name,
    Teuchos::RCP<const Epetra_MultiVector> vector, IO::VectorType type) const
{
  // if there will be more vectors to write use some string functionality!
  // do it when all configs are on C++11 because without its just annoying to convert
  // integers to strings (do the same in the post filter)
  // also in C++11 the vector can direclty be constructed from the range!
  const char* init[] = {"01", "02", "03", "04", "05", "05", "07", "08", "09", "10"};
  std::vector<std::string> numbers(init, init + 10);

  if ((unsigned int)vector->NumVectors() > numbers.size())
    dserror("extent the list above or make it automatic! Don't forget the post filter.");


  for (int i = 0; i < vector->NumVectors(); i++)
    output_->WriteVector(name + numbers[i], Teuchos::rcp((*vector)(i), false), type);

  return;
}


/*----------------------------------------------------------------------------*/
void INVANA::InvanaWriter::WriteMesh(int step, double time)
{
  output_->WriteMesh(step, time);

  return;
}
