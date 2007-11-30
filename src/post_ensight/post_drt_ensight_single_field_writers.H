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

#include "post_drt_ensight_writer.H"

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
    
//    /*!
//     \brief write all time steps of a result for xfem problems with varying number of unknowns per node
//
//     Write nodal results. The results are taken from a reconstructed
//     Epetra_Vector. We always assume a dofmanager that knows about the unknowns distribution
//     So we have to specify which XFEM::Physics::Field.
//
//     Finally, after writing to the result file, a string is returned that
//     describes the result for the case file VARIABLE section
//
//     \return string with entry for VARIABLE section in case file
//     \author a.ger
//     \date 11/07
//     */
//    string WriteResult(
//            const string groupname,                   ///< name of the result group in the control file
//            const string name,                        ///< name of the result to be written
//            const vector<XFEM::Physics::Field> fields ///< Fieldnames to print, could be just one (pressure) or a list (Dispx,Dispy,Dispz)
//            );
    
    
    
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
    
    PostField* cutterfield_;
};

#endif
