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

//! defines that a qudratic tet10 is temporarily written just as a tet4.
const int subtet10map[1][4] =
    {{0,1,2,3}};

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
    virtual string WriteAllResults(
            PostField* field) = 0;

    /*!
     \brief write all time steps of a result

     Write nodal results. The results are taken from a reconstructed
     Epetra_Vector. In many cases this vector will contain just one
     variable (displacements) and thus is easy to write as a whole. At
     other times, however, there is more than one result (velocity,
     pressure) and we want to write just one part of it. So we have to
     specify which part.

     Finally, after writing to the result file, a string is returned that
     describes the result for the case file VARIABLE section

     \return string with entry for VARIABLE section in case file
     \author u.kue
     \date 03/07
     */
    string WriteResult(
            const string groupname, ///< name of the result group in the control file
            const string name, ///< name of the result to be written
            const int numdf, ///< number of dofs per node to this result
            const int from=0 ///< start position of values in nodes
            );

    template<class T> void Write(
            ofstream& os,
            T i) const
    {
        os.write(reinterpret_cast<const char*>(&i), sizeof(T));
    }
    void Write(
            ofstream& os,
            const string s) const
    {
        WriteString(os, s);
    }
    void Write(
            ofstream& os,
            const char* s) const
    {
        WriteString(os, s);
    }

    void WriteString(
            ofstream& stream,  ///< filestream we are writing to
            const string str   ///< string to be written to file
            ) const;
    void WriteGeoFile(
            const string& geofilename) const;
    void WriteCoordinates(
            ofstream& geofile,                         ///< filestream for the geometry
            const RefCountPtr<DRT::Discretization> dis ///< discretization where the nodal positions are take from
            ) const;
    void WriteCells(
            ofstream& geofile,
            const RefCountPtr<DRT::Discretization> dis) const;
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
            const int numdf, ///< degrees of freedom per node for this field
            const unsigned int fileset,
            const string name,
            const string filename) const;
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
            unsigned int& fileset,                                   ///< ???
            vector<vector<int> >& filesets,                          ///< ???
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
    
    PostField* field_;
    string filename_;
    int myrank_;

    vector<vector<int> > filesets_;

    //! maps a distype to the corresponding Ensight cell type
    map<DRT::Element::DiscretizationType, string> distype2ensightstring_;
};


#endif /*POST_DRT_ENSIGHT_WRITER_H_*/
#endif
