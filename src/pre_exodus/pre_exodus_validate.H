/*----------------------------------------------------------------------*/
/*!
\file pre_exodus_validate.H

\brief validate a given .dat-file

\maintainer Martin Kronbichler

\level 1

Validate a given BACI input file (after all preprocessing steps)

*/

#ifndef PRE_EXODUS_VALIDATE_H_
#define PRE_EXODUS_VALIDATE_H_

#include <string>
#include <Teuchos_RCP.hpp>
#include <Epetra_SerialDenseMatrix.h>
#include "../drt_lib/drt_element.H"


namespace EXODUS{

// forward declaration
class Mesh;
class ElementBlock;

//! validate a given datfile
void ValidateInputFile(const Teuchos::RCP<Epetra_Comm> comm, const std::string datfile);

//! Check Elements for positive Jacobian and otherwise 'rewind' them
void ValidateMeshElementJacobians(EXODUS::Mesh& mymesh);

//! Check for positive Jacobian for Element of distype and otherwise 'rewind' them
void ValidateElementJacobian(EXODUS::Mesh& mymesh,const DRT::Element::DiscretizationType distype, Teuchos::RCP<EXODUS::ElementBlock>);

//! Check Elements of distype with full gauss integration rule for positive det at all gps and return number of negative dets
int ValidateElementJacobian_fullgp(Mesh& mymesh, const DRT::Element::DiscretizationType distype, Teuchos::RCP<ElementBlock> eb);

//! Check one element for positive Jacobi determinant
bool PositiveEle(const int& eleid, const std::vector<int>& nodes,const EXODUS::Mesh& mymesh,const Epetra_SerialDenseMatrix& deriv);
int EleSaneSign(const std::vector<int>& nodes,const std::map<int,std::vector<double> >& nodecoords);

//! Rewind one Element
std::vector<int> RewindEle(std::vector<int> old_nodeids, const DRT::Element::DiscretizationType distype);

}

#endif // PRE_EXODUS_VALIDATE_H_
