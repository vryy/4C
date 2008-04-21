/*----------------------------------------------------------------------*/
/*!
\file fsi_create_boundary.H

\brief

<pre>
Maintainer: Axel Gerstenberger
            gerstenberger@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>
*/
/*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#ifndef FSI_CREATE_BOUNDARY_H
#define FSI_CREATE_BOUNDARY_H

//! create discretization from a given condition
RCP<DRT::Discretization> CreateDiscretizationFromCondition(
        RCP<DRT::Discretization>  cutterdis,    ///< condition with discretization
        const string&             condname,     ///< name of the condition, by which the derived discretization is identified 
        const string&             discret_name, ///< name of the new condition
        const string&             element_name  ///< name/type of the elements to be created
        );

#endif
#endif
