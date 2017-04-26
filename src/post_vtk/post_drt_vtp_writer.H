/*----------------------------------------------------------------------*/
/*!
\file post_drt_vtp_writer.H

\brief VTP filter specialized for particle output

\level 2

\maintainer Martin Kronbichler
*-----------------------------------------------------------------------*/
#ifndef POST_DRT_VTP_H
#define POST_DRT_VTP_H


#include "post_drt_vtk_writer.H"

#include <string>
#include <map>
#include <vector>


// forward declarations
class PostField;


/*
 \brief Base class for VTP output generation

 \author hammerl
 \date 03/16
*/
class PostVtpWriter : public PostVtkWriter
{
public:
  //! constructor. Initializes the writer to a certain field.
  PostVtpWriter (PostField* field,
             const std::string &name);

  //! destructor
  virtual ~PostVtpWriter()
  {}

protected:

  //! Return the opening xml tag for this writer type
  virtual const std::string& WriterOpeningTag() const;

  //! Return the parallel opening xml tag for this writer type
  virtual const std::string& WriterPOpeningTag() const;

  //! Return a vector of parallel piece tags for each file
  virtual const std::vector<std::string>& WriterPPieceTags() const;

  //! Give every writer a chance to do preparations before writing
  virtual void WriterPrepTimestep() {};

  //! Return the parallel file suffix including the dot for this file type
  virtual const std::string& WriterPSuffix() const;

  //! Return the string of this writer type
  virtual const std::string& WriterString() const;

  //! Return the file suffix including the dot for this file type
  virtual const std::string& WriterSuffix() const;

  //! Write a single result step
  virtual void WriteDofResultStep(
      std::ofstream&                                                file,
      const Teuchos::RCP<Epetra_Vector>&                            data,
      std::map<std::string, std::vector<std::ofstream::pos_type> >& resultfilepos,
      const std::string&                                            groupname,
      const std::string&                                            name,
      const int                                                     numdf,
      const int                                                     from,
      const bool                                                    fillzeros
      );

  //! Write a single result step
  virtual void WriteNodalResultStep(
      std::ofstream&                                                file,
      const Teuchos::RCP<Epetra_MultiVector>&                       data,
      std::map<std::string, std::vector<std::ofstream::pos_type> >& resultfilepos,
      const std::string&                                            groupname,
      const std::string&                                            name,
      const int                                                     numdf
      );

  virtual void WriteElementResultStep(
      std::ofstream&                                                file,
      const Teuchos::RCP<Epetra_MultiVector>&                       data,
      std::map<std::string, std::vector<std::ofstream::pos_type> >& resultfilepos,
      const std::string&                                            groupname,
      const std::string&                                            name,
      const int                                                     numdf,
      const int                                                     from
      );

  //! write the geometry of one time step
  virtual void WriteGeo();

};

#endif
