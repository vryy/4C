/*!
 \file post_drt_ensight_single_field_writers.H

 \brief main routine of the Ensight filter

 <pre>
 Maintainer: Axel Gerstenberger
 gerstenberger@lnm.mw.tum.de
 http://www.lnm.mw.tum.de/Members/gerstenberger
 089 - 289-15236
 </pre>

 */
#ifdef CCADISCRET

#ifndef POST_DRT_ENSIGHT_SINGLE_FIELD_WRITER_H_
#define POST_DRT_ENSIGHT_SINGLE_FIELD_WRITER_H_

#include "post_drt_ensight_writer.H"
#include "../drt_xfem/physics.H"

/*
 \brief Writer for structural problems
 */
class StructureEnsightWriter : public EnsightWriter
{
public:
    StructureEnsightWriter(
            PostField* field,
            string filename) :
        EnsightWriter(field, filename)
    {
    }

protected:

    virtual void WriteAllResults(
            PostField* field);
};

/*----------------------------------------------------------------------*/
/*
 \brief Writer for fluid problems
 */
/*----------------------------------------------------------------------*/
class FluidEnsightWriter : public EnsightWriter
{
public:
    FluidEnsightWriter(
            PostField* field,
            string filename) :
        EnsightWriter(field, filename)
    {
    }

protected:

    virtual void WriteAllResults(
            PostField* field);
};

/*----------------------------------------------------------------------*/
/*
 \brief Writer for ale problems
 */
/*----------------------------------------------------------------------*/
class AleEnsightWriter : public EnsightWriter
{
public:
    AleEnsightWriter(
            PostField* field,
            string filename) :
        EnsightWriter(field, filename)
    {
    }

protected:

    virtual void WriteAllResults(
            PostField* field);
};


/*----------------------------------------------------------------------*/
/* 
 \brief Writer for convection-diffusion problems
 
 \author gjb
 \date 12/07
*/
/*----------------------------------------------------------------------*/
class ConDifEnsightWriter : public EnsightWriter
{
public:
    ConDifEnsightWriter(
            PostField* field,
            string filename) :
        EnsightWriter(field, filename)
    {
    }

protected:

    virtual void WriteAllResults(
            PostField* field);
};


/*----------------------------------------------------------------------*/
/*
 \brief Writer for fluid problems including xfem interfaces
 
 this might be handled more general in the near future, when integration cells
 are saved along the data for each timestep
 
 At the moment, we just intersect in the post filter again and conlude on the DOF 
 distribution which then should fit to the unknown vector. 
 hence, the xfem post filter is doing quite some work at this moment...
 
 \author a.ger
 \date 11/07
 
 */
/*----------------------------------------------------------------------*/
class XFluidEnsightWriter : public EnsightWriter
{
public:
    XFluidEnsightWriter(
            PostField* field,           ///< Field to be plotted
            PostField* cutterfield,     ///< intersecting field (for fsi)
            string filename) :
        EnsightWriter(field, filename),
        cutterfield_(cutterfield)
    {
    }
            
    //! write the whole thing
    void WriteFiles();
            
protected:

    virtual void WriteAllResults(
            PostField* field);
    void WriteGeoFile(
            const string& geofilename) const;
    void WriteGeoFileOneTimeStep(
            ofstream& file,
            map<string, vector<ofstream::pos_type> >& resultfilepos,
            const string name) const;
    void WriteGeoFilePart(
            ofstream& file,
            map<string, vector<ofstream::pos_type> >& resultfilepos,
            const string name,
            const RefCountPtr<XFEM::InterfaceHandle> ih ///< interfacehandle        
            ) const;
    int NumNodesPerField(
            const RefCountPtr<XFEM::InterfaceHandle> ih ///< interfacehandle
            ) const;
    void WriteCoordinates(
            ofstream& geofile,                          ///< filestream for the geometry
            const RefCountPtr<DRT::Discretization> dis, ///< discretization where the nodal positions are take from
            const RefCountPtr<XFEM::InterfaceHandle> ih ///< interfacehandle
            ) const;
    void WriteCells(
            ofstream& geofile,
            const RefCountPtr<DRT::Discretization> dis,
            const RefCountPtr<XFEM::InterfaceHandle> ih ///< interfacehandle
            ) const;
    /*!
     * \brief estimate, how many elements of each distype will be written
     * \return map between distype and number of elements to be written
     */
    NumElePerDisType GetNumElePerDisType(
            const RefCountPtr<DRT::Discretization> dis,
            const RefCountPtr<XFEM::InterfaceHandle> ih) const;
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
    void WriteResult(
            const string groupname, ///< name of the result group in the control file
            const string name,      ///< name of the result to be written
            const set<XFEM::PHYSICS::Field> fieldset
            );
    void WriteNodalResultStep(
            ofstream& file,
            PostResult& result,
            map<string, vector<ofstream::pos_type> >& resultfilepos,
            const string groupname,
            const string name,
            const set<XFEM::PHYSICS::Field> fieldset,
            const RCP<XFEM::InterfaceHandle> ih,
            const RCP<XFEM::DofManager> dofman
            ) const;
    
    PostField* cutterfield_;
};


#endif /*POST_DRT_ENSIGHT_SINGLE_FIELD_WRITER_H_*/
#endif
