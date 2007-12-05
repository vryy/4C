/*!
 \file post_drt_ensight_writer.H

 \brief basis class for the Ensight filter 

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
#include "../drt_lib/drt_node.H"
#include "../drt_lib/drt_element.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_xfem/dof_management.H"

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

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
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
    virtual void WriteAllResults(
            PostField* field) = 0;

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
    
    /*!
     \brief write a formated string to already opened filestream.
            currently used for writing the ensight casefile
     \author gjb
     \date 11/07
     */
    void WriteFormatedString(
            ofstream& os,
            const string s) const
    {
    	// only processor 0 does the writing !! 
    	if (myrank_==0)
    	  {
    		// first ensure, that the file is already open 
    	    if (os.is_open()) 
    	      {os << s;}
    	    else
    	      dserror("File is not open");
    	  }
    }

    /*!
     * \brief writes int like a formated string to already open filestream
     */  
    void WriteFormatedString(
            ofstream& os,
            const int i) const
    {
    	// only processor 0 does the writing !! 
    	if (myrank_==0)
    	  {
    		// first ensure, that the file is already open 
    	    if (os.is_open()) 
    	      {os << i;}
    	    else
    	      dserror("File is not open");
    	  }
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
    void WriteCoordinates(
            ofstream& geofile,                         ///< filestream for the geometry
            const RefCountPtr<DRT::Discretization> dis ///< discretization where the nodal positions are take from
            ) const;
    RefCountPtr<Epetra_Map> WriteCoordinatesPar(
            ofstream& geofile,                         ///< filestream for the geometry
            const RefCountPtr<DRT::Discretization> dis ///< discretization where the nodal positions are take from
            ) const;
    void WriteCells(
            ofstream& geofile,
            const RefCountPtr<DRT::Discretization> dis) const;
    void WriteCellsPar(
            ofstream& geofile,///< filestream for the geometry
            const RefCountPtr<DRT::Discretization> dis, ///< discretization where the nodal positions are take from
            const RefCountPtr<Epetra_Map> proc0map ///< current proc0 node map, created by WriteCoordinatesPar
            ) const;
    void WriteNodeConnectivityPar(
            ofstream& geofile,
            const RefCountPtr<DRT::Discretization> dis,
    	    const vector<int>& nodevector,
            const RefCountPtr<Epetra_Map> proc0map
            ) const;
    void WriteResultStep(
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
            string name,
            string filename
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
    
    //! create string for one TIME section in the case file
    string GetTimeSectionString(
            const int timeset,            ///< number of timeset to be written
            const vector<double>& times   ///< vector with time value for each time step
            ) const;
    
    //! create string for FILE section in the case file
    string GetFileSectionStringFromFilesets(
            const map<string,vector<int> >& filesetmap  ///< filesets when using multiple huge binary files
            ) const;
    
    PostField* field_;
    string filename_;
    int myrank_; ///< global processor id
    bool nodeidgiven_; ///< indicates whether BACI global node ids are written to geometry file. default value: true
    
    RCP<Epetra_Map> proc0map_; ///< allreduced nodal map for proc 0

    map<string,vector<int> > filesetmap_;
    
    map<string,vector<int> > timesetmap_;
    
    map<string,int>          variablenumdfmap_;
    map<string,string>       variablefilenamemap_;


    //! maps a distype to the corresponding Ensight cell type
    map<DRT::Element::DiscretizationType, string> distype2ensightstring_;
    
    // maximum file size
    static const unsigned FILE_SIZE_LIMIT_ = 0x7fffffff; // 2GB
    //static const unsigned FILE_SIZE_LIMIT_ = 1024*10; // 10kB ... useful for debugging ;-)
    
};


#endif /*POST_DRT_ENSIGHT_WRITER_H_*/
#endif
