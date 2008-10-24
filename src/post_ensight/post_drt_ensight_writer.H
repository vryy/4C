/*!
 \file post_drt_ensight_writer.H

 \brief contains basis class for the Ensight filter

 <pre>
 Maintainer: Ulrich Kuettler
 kuettler@lnm.mw.tum.de
 http://www.lnm.mw.tum.de/Members/kuettler
 089 - 289-15238
 </pre>

 */
#ifdef CCADISCRET

#ifndef POST_DRT_ENSIGHT_WRITER_H_
#define POST_DRT_ENSIGHT_WRITER_H_

#include <iostream>
#include <sstream>
#include <cstdio>
#include <vector>
#include <map>
#include <fstream>
#include <string>

#include "../post_drt_common/post_drt_common.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_nurbs_discret.H"
#include "../drt_fem_general/drt_utils_nurbs_shapefunctions.H"
#include "../drt_lib/drt_control_point.H"
#include "../drt_lib/drt_node.H"
#include "../drt_lib/drt_element.H"
#include "../drt_lib/linalg_utils.H"

using namespace std;

//! defines how 4 quad4 elements are constructed from a quad9
const int subquadmap[4][4] =
    {
        { 0, 4, 8, 7 },
        { 4, 1, 5, 8 },
        { 8, 5, 2, 6 },
        { 7, 8, 6, 3 } };

//! defines how 8 hex8 elements are constructed from a hex27
//  ;-) its symetric for some reason
const int subhexmap[8][8] =
    {
        { 0, 8, 20, 11, 12, 21, 26, 24 },
        { 8, 1, 9, 20, 21, 13, 22, 26 },
        { 20, 9, 2, 10, 26, 22, 14, 23 },
        { 11, 20, 10, 3, 24, 26, 23, 15 },
        { 12, 21, 26, 24, 4, 16, 25, 19 },
        { 21, 13, 22, 26, 16, 5, 17, 25 },
        { 26, 22, 14, 23, 25, 17, 6, 18 },
        { 24, 26, 23, 15, 19, 25, 18, 7 } };

typedef map<DRT::Element::DiscretizationType, int> NumElePerDisType;

typedef map<DRT::Element::DiscretizationType, vector<int> > EleGidPerDisType;

//! defines different result types (node-based, element-based, Gauss-Points?)
enum ResultType
{
  no_restype,          		// unknown result type
  dofbased,                     // dof-based result values
  nodebased,           		// node-based result values
  elementbased,      		// element-based result values
  elementdof,     		// element degree of freedom values
  max_restype          		// end marker. must be the last entry
};

//! basis class for the Ensight filter
class EnsightWriter
{
public:
    //! constructor
    EnsightWriter(
            PostField* field,
            const string filename);
    //! destructor
    virtual ~EnsightWriter()
    {
    }

    //! write the whole thing
    void WriteFiles();

protected:

    //! look for problem dependent result entries and write them
    virtual void WriteAllResults(PostField* field) = 0;

    /// write all element based results
    /*!
      \note This method sees only those element results that are defined in
      the first result group.
     */
    void WriteElementResults(PostField* field) { WriteAnyResults(field,"element",elementbased); }

    /// write all node based results
    /*!
      \note This method sees only those node results that are defined in
      the first result group.
     */
    void WriteNodeResults(PostField* field) { WriteAnyResults(field,"node",nodebased); }

    /// write all dof based results
    /*!
      \note This method sees only those dof results that are defined in
      the first result group.
     */
    void WriteDofResults(PostField* field) { WriteAnyResults(field,"dof",dofbased); }

    /// write all results of the given type
    /*!
      \note This method sees only those results that are defined in the first
      result group (that is the first time step).
     */
    void WriteAnyResults(PostField* field, const char* type, const ResultType restype);

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
    void WriteResult(
            const string groupname, ///< name of the result group in the control file
            const string name, ///< name of the result to be written
            const ResultType restype, ///< type of the result to be written (nodal-/element-based)
            const int numdf, ///< number of dofs per node to this result
            const int from=0 ///< start position of values in nodes
            );

    template<class T> void Write(
            ofstream& os,
            T i) const
    {
    	// only processor 0 does the writing !!
        if (myrank_==0)
    	os.write(reinterpret_cast<const char*>(&i), sizeof(T));
    }

    void Write(
            ofstream& os,
            const string s) const
    {
    	// only processor 0 does the writing !!
        if (myrank_==0)
        WriteString(os, s);
    }

    void Write(
            ofstream& os,
            const char* s) const
    {
    	// only processor 0 does the writing !!
        if (myrank_==0)
        WriteString(os, s);
    }

    void WriteString(
            ofstream& stream,  ///< filestream we are writing to
            const string str   ///< string to be written to file
            ) const;
    void WriteGeoFile(
            const string& geofilename);
    void WriteGeoFileOneTimeStep(
            ofstream& file,
            map<string, vector<ofstream::pos_type> >& resultfilepos,
            const string name);

    RefCountPtr<Epetra_Map> WriteCoordinates(
            ofstream& geofile,                         ///< filestream for the geometry
            const RefCountPtr<DRT::Discretization> dis ///< discretization where the nodal positions are take from
            ) ;

  /*! \brief Write the coordinates for a Polynomial discretization
    The ccordinates of the vizualisation points (i.e. the corner 
    nodes of elements displayed in paraview) are just the node 
    coordinates of the nodes in the discretization.
  */
    void WriteCoordinatesForPolynomialShapefunctions(
      ofstream&                              geofile ,
      const RefCountPtr<DRT::Discretization> dis     ,
      RefCountPtr<Epetra_Map>&               proc0map
      );

  /*! \brief Write the coordinates for a Nurbs discretization
    The ccordinates of the vizualisation points (i.e. the corner 
    nodes of elements displayed in paraview) are not the control point 
    coordinates of the nodes in the discretization but the points the
    knot values are mapped to.
  */
    void WriteCoordinatesForNurbsShapefunctions(
      ofstream&                              geofile ,
      const RefCountPtr<DRT::Discretization> dis     ,
      RefCountPtr<Epetra_Map>&               proc0map
      );


    void WriteCells(
            ofstream& geofile,///< filestream for the geometry
            const RefCountPtr<DRT::Discretization> dis, ///< discretization where the nodal positions are take from
            const RefCountPtr<Epetra_Map>& proc0map ///< current proc0 node map, created by WriteCoordinatesPar
            ) const;


  /*! \brief Write the cells for a Nurbs discretization
    quadratic nurbs split one element in knot space into 
    four(2d)/eight(3d) cells. The global numbering of the 
    vizualisation points (i.e. the corner points of the 
    cells) is computed from the local patch numbering and 
    the patch offset.                              (gammi)

    \param DRT::Element::DiscretizationType (i)          the nurbs discretisation type
    \param int                              (i)          global element id
    \param ofstream                         (used for o) direct print to file
    \param vector<int>                      (o)          remember node values for parallel IO
    \param RefCountPtr<DRT::Discretization> (i)          the discretisation holding knots etc
    \param RefCountPtr<Epetra_Map>          (i)          an allreduced nodemap

  */  
    void WriteNurbsCell(const DRT::Element::DiscretizationType distype   ,
                        const int                              gid       ,
                        ofstream&                              geofile   ,
                        vector<int>                            nodevector,
                        const RefCountPtr<DRT::Discretization> dis       ,
                        const RefCountPtr<Epetra_Map>&         proc0map  ) const;


    void WriteNodeConnectivityPar(
            ofstream& geofile,
            const RefCountPtr<DRT::Discretization> dis,
    	    const vector<int>& nodevector,
            const RefCountPtr<Epetra_Map> proc0map
            ) const;
    void WriteDofResultStep(
            ofstream& file,
            PostResult& result,
            map<string, vector<ofstream::pos_type> >& resultfilepos,
            const string groupname,
            const string name,
            const int numdf,
            const int from) const;


  /*! \brief Write the results for a NURBS discretisation
    (dof based).
    
    On input, result data for an ndimensional computation
    is provided (from the result file)

    This element data is communicated in such a way that
    all elements have access to their (dof-accessible) data.
    Here we seperate velocity and pressure output, since
    for velocity and pressure different dofs are required.

    Then, all elements are looped and function values are
    evaluated at visualisation points. This is the place 
    where we need the dof data (again, different data for 
    velocity and pressure output)

    The resulting vector is allreduced on proc0 and written.

                         .                              (gammi)

    \param ofstream                         (used for o) direct print to file
    \param int                              (i)          number of degrees of freedom
    \param RefCountPtr<Epetra_Vector>       (i)          the result data read from the baci output
    \param string                           (i)          name of the thing we are writing
                                                         (velocity, pressure etc.)

  */  
  void WriteDofResultStepForNurbs(
    ofstream&                        file ,
    const int                        numdf,
    const RefCountPtr<Epetra_Vector> data ,
    const string                     name
    ) const;


    void WriteNodalResultStep(
            ofstream& file,
            PostResult& result,
            map<string, vector<ofstream::pos_type> >& resultfilepos,
            const string groupname,
            const string name,
            const int numdf) const;
    void WriteNodalStressStep(
            ofstream& file,
            const RefCountPtr<Epetra_MultiVector> nodalstresses,
            map<string, vector<ofstream::pos_type> >& resultfilepos,
            const string name,
            const int numdf) const;
    void WriteElementCenterStressStep(
            ofstream& file,
            const RefCountPtr<Epetra_MultiVector> data,
            map<string, vector<ofstream::pos_type> >& resultfilepos,
            const string name,
            const int numdf) const;
    void WriteElementDOFResultStep(
            ofstream& file,
            PostResult& result,
            map<string, vector<ofstream::pos_type> >& resultfilepos,
            const string groupname,
            const string name,
            const int numdf,
            const int from) const;
    void WriteElementResultStep(
            ofstream& file,
            PostResult& result,
            map<string, vector<ofstream::pos_type> >& resultfilepos,
            const string groupname,
            const string name,
            const int numdf,
            const int from) const;
    void WriteIndexTable(
            ofstream& file,
            const vector<ofstream::pos_type>& filepos) const;
    /*!
     * \brief create string for the VARIABLE section
     *        that corresponds to the current field
     */
    string GetVariableEntryForCaseFile(
            int numdf,               ///< degrees of freedom per node for this field
            unsigned int fileset,
            const string name,
            const string filename,
            const int timeset
            ) const;
    /*!
     * \brief create string for the VARIABLE section
     *        for all fields in the variablemap
     */
    string GetVariableSection(
            map<string,vector<int> > filesetmap,
            map<string,int>          variablenumdfmap,
            map<string,string>       variablefilenamemap
            ) const;


    /*!
     * \brief estimate, how many elements of each distype will be written
     * \return map between distype and number of elements to be written
     */
    NumElePerDisType GetNumElePerDisType(
            const RefCountPtr<DRT::Discretization> dis) const;

    string GetEnsightString(
            const DRT::Element::DiscretizationType distype) const;

    /*!
     * \brief if files become to big, this routine switches to a new file
     */
    void FileSwitcher(
            ofstream& file,                                          ///< filestream that is switched
            bool& multiple_files,                                    ///< ???
            map<string,vector<int> >& filesetmap,                    ///< ???
            map<string, vector<ofstream::pos_type> >& resultfilepos, ///< ???
            const int stepsize,                                      ///< ???
            const string name,                                       ///< ???
            const string filename                                    ///< constant part of the filename
            ) const;

    int GetNumEleOutput(
            const DRT::Element::DiscretizationType distype,
            const int numele) const;

    /*!
     * \brief store, which elements are belonging to each present distype
     * \return map between distype and vector containing the global ids of th corresponding elements
     */
    EleGidPerDisType GetEleGidPerDisType(
            const RefCountPtr<DRT::Discretization> dis,
            NumElePerDisType numeleperdistype
            ) const;

    //! create string for one TIME section in the case file
    string GetTimeSectionString(
            const int timeset,            ///< number of timeset to be written
            const vector<double>& times   ///< vector with time value for each time step
            ) const;

    //! create string for TIME section in the case file
    string GetTimeSectionStringFromTimesets(
      const map<string,vector<double> >& timesetmap
      ) const;

    //! create string for FILE section in the case file
    string GetFileSectionStringFromFilesets(
            const map<string,vector<int> >& filesetmap  ///< filesets when using multiple huge binary files
            ) const;

    PostField* field_;
    string filename_;
    int myrank_; ///< global processor id
    bool nodeidgiven_; ///< indicates whether BACI global node ids are written to geometry file. default value: true

    NumElePerDisType numElePerDisType_; ///< number of elements per element discretization type
    EleGidPerDisType eleGidPerDisType_; ///< global ids of corresponding elements per element discretization type

    RCP<Epetra_Map> proc0map_; ///< allreduced node row map for proc 0, empty on other procs

    Teuchos::RCP<Epetra_Map> vispointmap_; ///< map for all visualisation points
    
    map<string,vector<int> > filesetmap_;
    map<string,vector<double> > timesetmap_;

    map<string,int>          variablenumdfmap_;
    map<string,string>       variablefilenamemap_;
    map<string,string>       variableresulttypemap_; ///< STL map for storing the result-type per variable name

    map<string,int> timesetnumbermap_;
    map<string,int> filesetnumbermap_;

    //! maps a distype to the corresponding Ensight cell type
    map<DRT::Element::DiscretizationType, string> distype2ensightstring_;

    //! maximum file size
    static const unsigned FILE_SIZE_LIMIT_ = 0x7fffffff; // 2GB
    //static const unsigned FILE_SIZE_LIMIT_ = 1024*10; // 10kB ... useful for debugging ;-)

};

#endif /*POST_DRT_ENSIGHT_WRITER_H_*/
#endif
