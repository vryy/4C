/*----------------------------------------------------------------------------*/
/*!
\file invana_writer.cpp

\brief Output for the inverse analysis control routine

<pre>
\level 3
\maintainer Sebastian Kehl
            kehl@mhpc.mw.tum.de
            089 - 289-10361
</pre>
*/
/*----------------------------------------------------------------------------*/
/* headers */
#include "invana_writer.H"

// baci
#include "../drt_io/io.H"
#include "../drt_lib/drt_dserror.H"

/*----------------------------------------------------------------------------*/
INVANA::InvanaWriter::InvanaWriter() :
output_(Teuchos::null)
{
  return;
}

/*----------------------------------------------------------------------------*/
void INVANA::InvanaWriter::Init(Teuchos::RCP<IO::DiscretizationWriter> writer)
{
  output_=writer;

  return;
}

/*----------------------------------------------------------------------------*/
void INVANA::InvanaWriter::WriteNewStep(const int step)
{
  output_->NewStep(step,double(step));
  return;
}

void INVANA::InvanaWriter::WriteNamedInt(const std::string name, const int aint)
{
  output_->WriteInt(name,aint);
  return;
}

/*----------------------------------------------------------------------------*/
void INVANA::InvanaWriter::WriteNamedVector(const std::string name,
    Teuchos::RCP<const Epetra_MultiVector> vector) const
{
  output_->WriteVector(name,vector);
  return;
}
