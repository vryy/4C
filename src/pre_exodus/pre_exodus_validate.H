/*----------------------------------------------------------------------*/
/*!
\file pre_exodus_validate.H

\brief validate a given .dat-file

<pre>
Maintainer: Georg Bauer
            bauer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/bauer
            089 - 289-15252
</pre>

Validate a given BACI input file (after all preprocessing steps)

*/

#ifdef D_EXODUS
#ifdef CCADISCRET

#ifndef PRE_EXODUS_VALIDATE_H_
#define PRE_EXODUS_VALIDATE_H_

#include <string>
#include "pre_exodus_reader.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_inputreader.H"
#include "../drt_lib/standardtypes_cpp.H"
#include "../drt_fem_general/drt_utils_integration.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include <Teuchos_RefCountPtr.hpp>

#ifdef PARALLEL
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

using namespace std;

namespace EXODUS{
//! validate a given datfile
void ValidateInputFile(const string datfile);

//! Check Elements for positive Jacobian and otherwise 'rewind' them
void ValidateMeshElementJacobians(EXODUS::Mesh& mymesh);

//! Check for positive Jacobian for Element of distype and otherwise 'rewind' them
void ValidateElementJacobian(EXODUS::Mesh& mymesh,const DRT::Element::DiscretizationType distype, RCP<EXODUS::ElementBlock>);
//! Check Elements of distype with full gauss integration rule for positive det at all gps and return number of negative dets
int ValidateElementJacobian_fullgp(Mesh& mymesh, const DRT::Element::DiscretizationType distype, RCP<ElementBlock> eb);

//! Check one element for positive Jacobian
bool PositiveEle(const vector<int>& nodes,const EXODUS::Mesh& mymesh,const Epetra_SerialDenseMatrix& deriv);
int EleSaneSign(const vector<int>& nodes,const map<int,vector<double> >& nodecoords);

//! Rewind one Element
vector<int> RewindEle(vector<int> old_nodeids, const DRT::Element::DiscretizationType distype);

}

#endif // PRE_EXODUS_VALIDATE_H_

#endif
#endif
