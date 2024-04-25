/*----------------------------------------------------------------------*/
/*! \file

\brief VTI filter


\level 2
*/
/*----------------------------------------------------------------------*/
#ifndef FOUR_C_POST_VTK_VTI_WRITER_HPP
#define FOUR_C_POST_VTK_VTI_WRITER_HPP


#include "4C_config.hpp"

#include "4C_post_vtk_writer.hpp"

#include <Teuchos_RCP.hpp>

#include <map>
#include <string>
#include <vector>

// forward declarations
class Epetra_Vector;

FOUR_C_NAMESPACE_OPEN
class PostField;
class PostResult;


namespace DRT
{
  class Discretization;
  class Node;
}  // namespace DRT


/*
 \brief Base class for VTU output generation

 \author kronbichler
 \date 03/14
*/
class PostVtiWriter : public PostVtkWriter
{
 public:
  //! constructor. Initializes the writer to a certain field.
  PostVtiWriter(PostField* field, const std::string& name);

 protected:
  //! Return the opening xml tag for this writer type
  const std::string& WriterOpeningTag() const override;

  //! Return the parallel opening xml tag for this writer type
  const std::string& WriterPOpeningTag() const override;

  //! Return a vector of parallel piece tags for each file
  const std::vector<std::string>& WriterPPieceTags() const override;

  //! Give every writer a chance to do preparations before writing
  void WriterPrepTimestep() override;

  //! Return the parallel file suffix including the dot for this file type
  const std::string& WriterPSuffix() const override;

  //! Return the string of this writer type
  const std::string& WriterString() const override;

  //! Return the file suffix including the dot for this file type
  const std::string& WriterSuffix() const override;

  //! Write a single result step
  void WriteDofResultStep(std::ofstream& file, const Teuchos::RCP<Epetra_Vector>& data,
      std::map<std::string, std::vector<std::ofstream::pos_type>>& resultfilepos,
      const std::string& groupname, const std::string& name, const int numdf, const int from,
      const bool fillzeros) override;

  //! Write a single result step
  void WriteNodalResultStep(std::ofstream& file, const Teuchos::RCP<Epetra_MultiVector>& data,
      std::map<std::string, std::vector<std::ofstream::pos_type>>& resultfilepos,
      const std::string& groupname, const std::string& name, const int numdf) override;

  //! Write a single result step
  void WriteElementResultStep(std::ofstream& file, const Teuchos::RCP<Epetra_MultiVector>& data,
      std::map<std::string, std::vector<std::ofstream::pos_type>>& resultfilepos,
      const std::string& groupname, const std::string& name, const int numdf,
      const int from) override;

  //! write the geometry of one time step
  void WriteGeo() override;

  //! origin of the ImageData-grid
  double origin_[3];

  //! spacing of the ImageData-grid
  double spacing_[3];

  //! global extent of the ImageData-grid (x_min x_max y_min y_max z_min z_max)
  int globalextent_[6];

  //! local extent of the ImageData-grid (x_min x_max y_min y_max z_min z_max)
  int localextent_[6];

  //! Mapping between nodeids and their position on an ImageData-grid in a (z*Ny+y)*Nx+x form
  std::map<int, int> idmapping_;

  //! Mapping between elementids and their position on an ImageData-grid in a (z*Ny+y)*Nx+x form
  std::map<int, int> eidmapping_;
};

FOUR_C_NAMESPACE_CLOSE

#endif
