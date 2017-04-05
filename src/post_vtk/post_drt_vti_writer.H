/*----------------------------------------------------------------------*/
/*!
\file post_drt_vti_writer.H

\brief VTI filter

\maintainer Karl-Robert Wichmann

\level 2
*/
/*----------------------------------------------------------------------*/
#ifndef POST_DRT_VTI_H
#define POST_DRT_VTI_H


#include "post_drt_vtk_writer.H"

#include <string>
#include <map>
#include <vector>

#include <Teuchos_RCP.hpp>


// forward declarations
class Epetra_Vector;
class PostField;
class PostResult;
namespace DRT
{
  class Discretization;
  class Node;
}


/*
 \brief Base class for VTU output generation

 \author kronbichler
 \date 03/14
*/
class PostVtiWriter : public PostVtkWriter
{
public:
  //! constructor. Initializes the writer to a certain field.
  PostVtiWriter (PostField* field,
             const std::string &name);

  //! destructor
  virtual ~PostVtiWriter()
  {}

protected:

  //! Return the opening xml tag for this writer type
  virtual const std::string& WriterOpeningTag() const;

  //! Return the parallel opening xml tag for this writer type
  virtual const std::string& WriterPOpeningTag() const;

  //! Return a vector of parallel piece tags for each file
  virtual const std::vector<std::string>& WriterPPieceTags() const;

  //! Give every writer a chance to do preparations before writing
  virtual void WriterPrepTimestep();

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

  //! Write a single result step
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

  //! origin of the ImageData-grid
  double origin_[3];

  //! spacing of the ImageData-grid
  double spacing_[3];

  //! global extent of the ImageData-grid (x_min x_max y_min y_max z_min z_max)
  int globalextent_[6];

  //! local extent of the ImageData-grid (x_min x_max y_min y_max z_min z_max)
  int localextent_[6];

  //! Mapping between nodeids and their position on an ImageData-grid in a (z*Ny+y)*Nx+x form
  std::map<int,int> idmapping_;

  //! Mapping between elementids and their position on an ImageData-grid in a (z*Ny+y)*Nx+x form
  std::map<int,int> eidmapping_;
};

#endif
