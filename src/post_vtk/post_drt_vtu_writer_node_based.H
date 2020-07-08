/*----------------------------------------------------------------------*/
/*! \file

\brief node based VTU filter

\level 3

*/
/*----------------------------------------------------------------------*/
#ifndef POST_DRT_VTU_NODE_H
#define POST_DRT_VTU_NODE_H


#include "post_drt_vtu_writer.H"

#include <string>
#include <map>
#include <vector>

#include "../drt_lib/drt_element.H"


// forward declarations
class PostField;
class PostResult;
namespace DRT
{
  class Discretization;
  class Node;

  namespace ELEMENTS
  {
    class Beam3Base;
  }
}  // namespace DRT


/*
 \brief Base class for VTU node based output generation

 \author biehler
 \date 12/16
*/
class PostVtuWriterNode : public PostVtuWriter
{
 public:
  //! constructor. Initializes the writer to a certain field.
  PostVtuWriterNode(PostField* field, const std::string& name);

  //! destructor
  virtual ~PostVtuWriterNode() {}

 protected:
  //! Return the opening xml tag for this writer type
  virtual const std::string& WriterOpeningTag() const;

  //! Return the parallel opening xml tag for this writer type
  virtual const std::string& WriterPOpeningTag() const;

  //! Return a vector of parallel piece tags for each file
  virtual const std::vector<std::string>& WriterPPieceTags() const;

  //! Give every writer a chance to do preparations before writing
  virtual void WriterPrepTimestep(){};

  //! Return the parallel file suffix including the dot for this file type
  virtual const std::string& WriterPSuffix() const;

  //! Return the string of this writer type
  virtual const std::string& WriterString() const;

  //! Return the file suffix including the dot for this file type
  virtual const std::string& WriterSuffix() const;

  //! Write a single result step
  virtual void WriteDofResultStep(std::ofstream& file, const Teuchos::RCP<Epetra_Vector>& data,
      std::map<std::string, std::vector<std::ofstream::pos_type>>& resultfilepos,
      const std::string& groupname, const std::string& name, const int numdf, const int from,
      const bool fillzeros);

  //! Write a single result step
  virtual void WriteNodalResultStep(std::ofstream& file,
      const Teuchos::RCP<Epetra_MultiVector>& data,
      std::map<std::string, std::vector<std::ofstream::pos_type>>& resultfilepos,
      const std::string& groupname, const std::string& name, const int numdf);

  //! Write a single result step
  virtual void WriteElementResultStep(std::ofstream& file,
      const Teuchos::RCP<Epetra_MultiVector>& data,
      std::map<std::string, std::vector<std::ofstream::pos_type>>& resultfilepos,
      const std::string& groupname, const std::string& name, const int numdf, const int from);

  //! write the geometry of one time step
  virtual void WriteGeo();

  //! write the geometry of Nurbs Element
  virtual void WriteGeoNurbsEle(const DRT::Element* ele, std::vector<uint8_t>& celltypes,
      int& outNodeId, std::vector<int32_t>& celloffset, std::vector<double>& coordinates);

  //! write the geometry of beam element (special treatment due to Hermite interpolation)
  virtual void WriteGeoBeamEle(const DRT::ELEMENTS::Beam3Base* beamele,
      std::vector<uint8_t>& celltypes, int& outNodeId, std::vector<int32_t>& celloffset,
      std::vector<double>& coordinates);

  //! Write a single result step for one Nurbs Element
  virtual void WirteDofResultStepNurbsEle(const DRT::Element* ele, int ncomponents, const int numdf,
      std::vector<double>& solution, Teuchos::RCP<Epetra_Vector> ghostedData, const int from,
      const bool fillzeros);

  virtual void WriteDofResultStepBeamEle(const DRT::ELEMENTS::Beam3Base* beamele,
      const int& ncomponents, const int& numdf, std::vector<double>& solution,
      Teuchos::RCP<Epetra_Vector>& ghostedData, const int& from, const bool fillzeros);

  //! Write a single result step for one Nurbs Element
  virtual void WriteNodalResultStepNurbsEle(const DRT::Element* ele, int ncomponents,
      const int numdf, std::vector<double>& solution, Teuchos::RCP<Epetra_MultiVector> ghostedData);
};

#endif
