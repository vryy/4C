/*!
 \file post_drt_ensight_writer.H

 \brief contains basis class for the Ensight filter

 <pre>
 Maintainer: Martin Kronbichler
 kronbichler@lnm.mw.tum.de
 http://www.lnm.mw.tum.de/Members/kronbichler
 089 - 289-15235
 </pre>

 */


#ifndef POST_DRT_ENSIGHT_WRITER_H_
#define POST_DRT_ENSIGHT_WRITER_H_

#include <iostream>
#include <sstream>
#include <cstdio>
#include <vector>
#include <map>
#include <fstream>
#include <string>

#include "../drt_lib/drt_element.H"               // due to DiscretizationType
#include "../post_drt_common/post_writer_base.H"  // base class PostWriterBase

// class PostProblem;
class PostField;
class PostResult;


// forward declaration
namespace DRT
{
  class Discretization;
  class Node;

  namespace NURBS
  {
    class NurbsDiscretization;
  }
}  // namespace DRT

//! defines how 2 line2 elements are constructed from a line3
const int sublinemap[2][2] = {{0, 2}, {2, 1}};


//! defines how 4 quad4 elements are constructed from a quad9
const int subquadmap[4][4] = {{0, 4, 8, 7}, {4, 1, 5, 8}, {8, 5, 2, 6}, {7, 8, 6, 3}};

//! defines how 8 hex8 elements are constructed from a hex27
//  ;-) its symetric for some reason
const int subhexmap[8][8] = {{0, 8, 20, 11, 12, 21, 26, 24}, {8, 1, 9, 20, 21, 13, 22, 26},
    {20, 9, 2, 10, 26, 22, 14, 23}, {11, 20, 10, 3, 24, 26, 23, 15},
    {12, 21, 26, 24, 4, 16, 25, 19}, {21, 13, 22, 26, 16, 5, 17, 25},
    {26, 22, 14, 23, 25, 17, 6, 18}, {24, 26, 23, 15, 19, 25, 18, 7}};

//! defines how 4 hex8 elements are constructed from a hex16
const int subhex16map[4][8] = {
    {0, 4, 7, 8, 12, 15}, {4, 1, 5, 12, 9, 13}, {5, 2, 6, 13, 10, 14}, {7, 6, 3, 14, 14, 12}};

//! defines how 4 hex8 elements are constructed from a hex18
const int subhex18map[4][8] = {{0, 4, 8, 7, 9, 13, 17, 16}, {4, 1, 5, 8, 13, 10, 14, 17},
    {8, 5, 2, 6, 17, 14, 11, 15}, {7, 8, 6, 3, 16, 17, 15, 12}};


//! basis class for the Ensight filter
class EnsightWriter : public PostWriterBase
{
 public:
  typedef std::map<DRT::Element::DiscretizationType, int> NumElePerDisType;

  typedef std::map<DRT::Element::DiscretizationType, std::vector<int>> EleGidPerDisType;

  //! constructor, does nothing (SetField must be called before usage)
  EnsightWriter(PostField* field, const std::string& name);

  //! destructor
  virtual ~EnsightWriter() {}

  //! write the whole thing
  virtual void WriteFiles(PostFilterBase& filter);

  //! write the whole thing for changing geometries
  virtual void WriteFilesChangingGeom(PostFilterBase& filter);

 protected:
  /*!
   \brief write all time steps of a result

   Write nodal results. The results are taken from a reconstructed
   Epetra_Vector. In many cases this vector will contain just one
   variable (displacements) and thus is easy to write as a whole. At
   other times, however, there is more than one result (velocity,
   pressure) and we want to write just one part of it. So we have to
   specify which part.

   \author u.kue
   \date 03/07
   */
  void WriteResult(const std::string groupname,  ///< name of the result group in the control file
      const std::string name,                    ///< name of the result to be written
      const ResultType restype,     ///< type of the result to be written (nodal-/element-based)
      const int numdf,              ///< number of dofs per node to this result
      const int from = 0,           ///< start position of values in nodes
      const bool fillzeros = false  ///< zeros are filled to ensight file when no data is available
  );


  /*!
   \brief write all time steps of a result in one time step

   Write nodal results. The results are taken from a reconstructed
   Epetra_Vector. In many cases this vector will contain just one
   variable (displacements) and thus is easy to write as a whole. At
   other times, however, there is more than one result (velocity,
   pressure) and we want to write just one part of it. So we have to
   specify which part. Currently file continuation is not supported
   because Paraview is not able to load it due to some weird wild card
   issue.

    originally
   \author u.kue
   \date 03/07
   adapted
   \author ghamm
   \date 03/13
   */
  void WriteResultOneTimeStep(PostResult& result,  ///< result group in the control file
      const std::string groupname,                 ///< name of the result group in the control file
      const std::string name,                      ///< name of the result to be written
      const ResultType restype,  ///< type of the result to be written (nodal-/element-based)
      const int numdf,           ///< number of dofs per node to this result
      bool firststep,            ///< bool whether this is the first time step
      bool laststep,             ///< bool whether this is the last time step
      const int from = 0         ///< start position of values in nodes
  );

  /*!
   \brief write a particular variable to file

   Write results. Some variables need interaction with the post filter,
   e.g. structural stresses that do some element computations before output.
   To allow for a generic interface, the calling site needs to supply a
   class derived from SpecialFieldInterface that knows which function to call.

   \author kronbichler
   \date 04/14
   */
  virtual void WriteSpecialField(SpecialFieldInterface& special,
      PostResult& result,  ///< result group in the control file
      const ResultType restype, const std::string& groupname,
      const std::vector<std::string>& fieldnames, const std::string& outinfo);

  template <class T>
  void Write(std::ofstream& os, T i) const
  {
    // only processor 0 does the writing !!
    if (myrank_ == 0) os.write(reinterpret_cast<const char*>(&i), sizeof(T));
  }

  void Write(std::ofstream& os, const std::string s) const
  {
    // only processor 0 does the writing !!
    if (myrank_ == 0) WriteString(os, s);
  }

  void Write(std::ofstream& os, const char* s) const
  {
    // only processor 0 does the writing !!
    if (myrank_ == 0) WriteString(os, s);
  }

  void WriteString(std::ofstream& stream,  ///< filestream we are writing to
      const std::string str                ///< string to be written to file
      ) const;
  void WriteGeoFile(const std::string& geofilename);
  void WriteGeoFileOneTimeStep(std::ofstream& file,
      std::map<std::string, std::vector<std::ofstream::pos_type>>& resultfilepos,
      const std::string name);

  Teuchos::RCP<Epetra_Map> WriteCoordinates(
      std::ofstream& geofile,  ///< filestream for the geometry
      const Teuchos::RCP<DRT::Discretization>
          dis  ///< discretization where the nodal positions are take from
  );

  /*! \brief Write the coordinates for a Polynomial discretization
      The cordinates of the vizualisation points (i.e. the corner
      nodes of elements displayed in paraview) are just the node
      coordinates of the nodes in the discretization.
    */
  void WriteCoordinatesForPolynomialShapefunctions(std::ofstream& geofile,
      const Teuchos::RCP<DRT::Discretization> dis, Teuchos::RCP<Epetra_Map>& proc0map);

  /*! \brief Write the coordinates for a Nurbs discretization
    The cordinates of the vizualisation points (i.e. the corner
    nodes of elements displayed in paraview) are not the control point
    coordinates of the nodes in the discretization but the points the
    knot values are mapped to.
  */
  void WriteCoordinatesForNurbsShapefunctions(std::ofstream& geofile,
      const Teuchos::RCP<DRT::Discretization> dis, Teuchos::RCP<Epetra_Map>& proc0map);

  virtual void WriteCells(std::ofstream& geofile,  ///< filestream for the geometry
      const Teuchos::RCP<DRT::Discretization>
          dis,  ///< discretization where the nodal positions are take from
      const Teuchos::RCP<Epetra_Map>&
          proc0map  ///< current proc0 node map, created by WriteCoordinatesPar
      ) const;

  /*! \brief Write the cells for a Nurbs discretization
    quadratic nurbs split one element in knot space into
    four(2d)/eight(3d) cells. The global numbering of the
    vizualisation points (i.e. the corner points of the
    cells) is computed from the local patch numbering and
    the patch offset.                              (gammi)

    \param DRT::Element::DiscretizationType (i)          the nurbs discretisation type
    \param int                              (i)          global element id
    \param std::ofstream                    (used for o) direct print to file
    \param std::vector<int>                 (o)          remember node values for parallel IO
    \param Teuchos::RCP<DRT::Discretization> (i)          the discretisation holding knots etc
    \param Teuchos::RCP<Epetra_Map>          (i)          an allreduced nodemap

  */
  void WriteNurbsCell(const DRT::Element::DiscretizationType distype, const int gid,
      std::ofstream& geofile, std::vector<int>& nodevector,
      const Teuchos::RCP<DRT::Discretization> dis, const Teuchos::RCP<Epetra_Map>& proc0map) const;

  /*! \brief Quadratic nurbs split one nurbs27 element
    in knot space into eight(3d) cells. The global
    numbering of the vizualisation points (i.e. the corner
    points of the cells) are computed from the local patch
    numbering and the patch offset. This method appends
    the node connectivity for one hex8 cell to the vector
    of cell nodes                                 (gammi)

    \param DRT::Element::DiscretizationType (i/o)        the vector to which the node
                                                         connectivity is appended to
    \param int                              (i)          0: left      1: right (which hex to
    generate) \param int                              (i)          0: front     1: rear  (which hex
    to generate) \param int                              (i)          0: bottom    1: top   (which
    hex to generate) \param std::vector<int>                 (i)          cartesian element ids in
    patch \param int                              (i)          number of visualisation points in u
    direction \param int                              (i)          number of visualisation points in
    v direction \param int                              (i)          number of patch

    \return void

  */
  void AppendNurbsSubHex(std::vector<int>& cellnodes, const int& l, const int& m, const int& n,
      const std::vector<int>& ele_cart_id, const int& nvpu, const int& nvpv,
      const int& npatch) const
  {
    int twoid[3];
    twoid[0] = 2 * ele_cart_id[0];
    twoid[1] = 2 * ele_cart_id[1];
    twoid[2] = 2 * ele_cart_id[2];

    cellnodes.push_back((twoid[0] + l) + ((twoid[1] + m) + (twoid[2] + n) * nvpv) * nvpu);
    cellnodes.push_back((twoid[0] + 1 + l) + ((twoid[1] + m) + (twoid[2] + n) * nvpv) * nvpu);
    cellnodes.push_back((twoid[0] + 1 + l) + ((twoid[1] + 1 + m) + (twoid[2] + n) * nvpv) * nvpu);
    cellnodes.push_back((twoid[0] + l) + ((twoid[1] + 1 + m) + (twoid[2] + n) * nvpv) * nvpu);
    cellnodes.push_back((twoid[0] + l) + ((twoid[1] + m) + (twoid[2] + 1 + n) * nvpv) * nvpu);
    cellnodes.push_back((twoid[0] + 1 + l) + ((twoid[1] + m) + (twoid[2] + 1 + n) * nvpv) * nvpu);
    cellnodes.push_back(
        (twoid[0] + 1 + l) + ((twoid[1] + 1 + m) + (twoid[2] + 1 + n) * nvpv) * nvpu);
    cellnodes.push_back((twoid[0] + l) + ((twoid[1] + 1 + m) + (twoid[2] + 1 + n) * nvpv) * nvpu);

    return;
  };


  void WriteNodeConnectivityPar(std::ofstream& geofile, const Teuchos::RCP<DRT::Discretization> dis,
      const std::vector<int>& nodevector, const Teuchos::RCP<Epetra_Map> proc0map) const;
  void WriteDofResultStep(std::ofstream& file, PostResult& result,
      std::map<std::string, std::vector<std::ofstream::pos_type>>& resultfilepos,
      const std::string& groupname, const std::string& name, const int numdf, const int from,
      const bool fillzeros) const;


  /*! \brief Write the results for a NURBS discretisation
    (dof based).

    On input, result data for an ndimensional computation
    is provided (from the result file)

    This element data is communicated in such a way that
    all elements have access to their (dof-accessible) data.
    Here we seperate velocity/displacement and pressure
    output, since for velocity/displacement and pressure
    different dofs are required.

    Then, all elements are looped and function values are
    evaluated at visualisation points. This is the place
    where we need the dof data (again, different data for
    velocity/displacement and pressure output)

    The resulting vector is allreduced on proc0 and written.

                         .                              (gammi)

    \param std::ofstream                    (used for o) direct print to file
    \param int                              (i)          number of degrees of freedom
    \param Teuchos::RCP<Epetra_Vector>      (i)          the result data read from the baci output
    \param string                           (i)          name of the thing we are writing
                                                         (velocity, pressure etc.)
    \param int                              (i)          potential offset in dof numbering

  */
  void WriteDofResultStepForNurbs(std::ofstream& file, const int numdf,
      const Teuchos::RCP<Epetra_Vector> data, const std::string name, const int offset) const;

  //! perform interpolation of result data to visualization points.
  void InterpolateNurbsResultToVizPoints(Teuchos::RCP<Epetra_MultiVector> idata, const int dim,
      const int npatch, const std::vector<int>& vpoff, const std::vector<int>& ele_cart_id,
      const DRT::Element* actele, DRT::NURBS::NurbsDiscretization* nurbsdis,
      const std::vector<Epetra_SerialDenseVector>& eleknots,
      const Epetra_SerialDenseVector& weights, const int numdf,
      const std::vector<double>& my_data) const;

  void WriteNodalResultStepForNurbs(std::ofstream& file, const int numdf,
      const Teuchos::RCP<Epetra_MultiVector> data, const std::string name, const int offset) const;

  void WriteNodalResultStep(std::ofstream& file, PostResult& result,
      std::map<std::string, std::vector<std::ofstream::pos_type>>& resultfilepos,
      const std::string& groupname, const std::string& name, const int numdf);
  virtual void WriteNodalResultStep(std::ofstream& file,
      const Teuchos::RCP<Epetra_MultiVector>& data,
      std::map<std::string, std::vector<std::ofstream::pos_type>>& resultfilepos,
      const std::string& groupname, const std::string& name, const int numdf);
  void WriteElementDOFResultStep(std::ofstream& file, PostResult& result,
      std::map<std::string, std::vector<std::ofstream::pos_type>>& resultfilepos,
      const std::string& groupname, const std::string& name, const int numdf, const int from) const;
  void WriteElementResultStep(std::ofstream& file, PostResult& result,
      std::map<std::string, std::vector<std::ofstream::pos_type>>& resultfilepos,
      const std::string& groupname, const std::string& name, const int numdf, const int from);
  virtual void WriteElementResultStep(std::ofstream& file,
      const Teuchos::RCP<Epetra_MultiVector>& data,
      std::map<std::string, std::vector<std::ofstream::pos_type>>& resultfilepos,
      const std::string& groupname, const std::string& name, const int numdf, const int from);
  void WriteIndexTable(
      std::ofstream& file, const std::vector<std::ofstream::pos_type>& filepos) const;

  /*!
   * \brief create string for the VARIABLE section
   *        that corresponds to the current field
   */
  std::string GetVariableEntryForCaseFile(
      int numdf,  ///< degrees of freedom per node for this field
      unsigned int fileset, const std::string name, const std::string filename,
      const int timeset) const;

  /*!
   * \brief create string for the VARIABLE section
   *        for all fields in the variablemap
   */
  std::string GetVariableSection(std::map<std::string, std::vector<int>> filesetmap,
      std::map<std::string, int> variablenumdfmap,
      std::map<std::string, std::string> variablefilenamemap) const;

  /*!
   * \brief estimate, how many elements of each distype will be written
   * \return map between distype and number of elements to be written
   */
  NumElePerDisType GetNumElePerDisType(const Teuchos::RCP<DRT::Discretization> dis) const;

  std::string GetEnsightString(const DRT::Element::DiscretizationType distype) const;

  /*!
   * \brief if files become to big, this routine switches to a new file
   */
  void FileSwitcher(std::ofstream& file,                    ///< filestream that is switched
      bool& multiple_files,                                 ///< ???
      std::map<std::string, std::vector<int>>& filesetmap,  ///< ???
      std::map<std::string, std::vector<std::ofstream::pos_type>>& resultfilepos,  ///< ???
      const int stepsize,                                                          ///< ???
      const std::string name,                                                      ///< ???
      const std::string filename  ///< constant part of the filename
      ) const;

  int GetNumEleOutput(const DRT::Element::DiscretizationType distype, const int numele) const;

  //! number of subelements (equals 1 if no element splitting has to be done)
  int GetNumSubEle(const DRT::Element::DiscretizationType distype) const;

  /*!
   * \brief store, which elements are belonging to each present distype
   * \return map between distype and vector containing the global ids of th corresponding elements
   */
  EleGidPerDisType GetEleGidPerDisType(
      const Teuchos::RCP<DRT::Discretization> dis, NumElePerDisType numeleperdistype) const;

  //! create string for one TIME section in the case file
  std::string GetTimeSectionString(const int timeset,  ///< number of timeset to be written
      const std::vector<double>& times  ///< vector with time value for each time step
      ) const;

  //! create string for TIME section in the case file
  std::string GetTimeSectionStringFromTimesets(
      const std::map<std::string, std::vector<double>>& timesetmap) const;

  //! create string for FILE section in the case file
  std::string GetFileSectionStringFromFilesets(const std::map<std::string, std::vector<int>>&
          filesetmap  ///< filesets when using multiple huge binary files
      ) const;

  bool nodeidgiven_;  ///< indicates whether BACI global node ids are written to geometry file.
                      ///< default value: true
  bool
      writecp_;  ///< NURBS-specific: defines if control point information should be written to file

  NumElePerDisType numElePerDisType_;  ///< number of elements per element discretization type
  EleGidPerDisType
      eleGidPerDisType_;  ///< global ids of corresponding elements per element discretization type

  Teuchos::RCP<Epetra_Map> proc0map_;  ///< allreduced node row map for proc 0, empty on other procs

  Teuchos::RCP<Epetra_Map> vispointmap_;  ///< map for all visualisation points

  std::map<std::string, std::vector<int>> filesetmap_;
  std::map<std::string, std::vector<double>> timesetmap_;

  std::map<std::string, int> variablenumdfmap_;
  std::map<std::string, std::string> variablefilenamemap_;
  std::map<std::string, std::string>
      variableresulttypemap_;  ///< STL map for storing the result-type per variable name

  std::map<std::string, int> timesetnumbermap_;
  std::map<std::string, int> filesetnumbermap_;

  std::map<std::string, std::vector<std::ofstream::pos_type>> resultfilepos_;

  //! maps a distype to the corresponding Ensight cell type
  std::map<DRT::Element::DiscretizationType, std::string> distype2ensightstring_;

  //! maximum file size
  static const unsigned FILE_SIZE_LIMIT_ = 0x7fffffff;  // 2GB
  // static const unsigned FILE_SIZE_LIMIT_ = 1024*10; // 10kB ... useful for debugging ;-)
};

#endif /*POST_DRT_ENSIGHT_WRITER_H_*/
