/*----------------------------------------------------------------------*/
/*!
\file drt_elementtype.H

\brief Type definitions for elements

<pre>
-------------------------------------------------------------------------
                 BACI finite element library subsystem
            Copyright (2008) Technical University of Munich

Under terms of contract T004.008.000 there is a non-exclusive license for use
of this work by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library is proprietary software. It must not be published, distributed,
copied or altered in any form or any media without written permission
of the copyright holder. It may be used under terms and conditions of the
above mentioned license by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library may solemnly used in conjunction with the BACI contact library
for purposes described in the above mentioned contract.

This library contains and makes use of software copyrighted by Sandia Corporation
and distributed under LGPL licence. Licensing does not apply to this or any
other third party software used here.

Questions? Contact Dr. Michael W. Gee (gee@lnm.mw.tum.de)
                   or
                   Prof. Dr. Wolfgang A. Wall (wall@lnm.mw.tum.de)

http://www.lnm.mw.tum.de

-------------------------------------------------------------------------
</pre>

<pre>
Maintainer: Martin Kronbichler
            kronbichler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>
*/
/*----------------------------------------------------------------------*/

#ifndef DRT_ELEMENTTPYE_H
#define DRT_ELEMENTTPYE_H

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Epetra_Vector.h>

#include "drt_parobjectfactory.H"
#include "drt_element.H"

namespace LINALG
{
  class SparseOperator;
}

namespace DRT
{

  class Element;
  class Discretization;

/// Subclass of ParObjectType that adds element type specific methods
/*!
  Element types need to be initialized. Furthermore, there is a PreEvaluate
  method and the ability to read elements from dat files. And finally the
  element specific setup of null spaces for multi grid preconditioning is
  here, too.

  \note There are boundary elements that do not need all of this
  functionality.

  \author u.kue
  \date 06/10
 */
class ElementType : public ParObjectType
{
protected:

  // only derived classes might create an instance
  ElementType();

public:

  /// setup the dat file input line definitions for this type of element
  virtual void SetupElementDefinition( std::map<std::string,std::map<std::string,DRT::INPUT::LineDefinition> > & definitions ) {}

  /// create an element from a dat file specifier
  virtual Teuchos::RCP<Element> Create( const std::string eletype,
                                        const std::string eledistype,
                                        const int id,
                                        const int owner )
  { return Teuchos::null; }

  /// create an empty element
  virtual Teuchos::RCP<DRT::Element> Create( const int id, const int owner ) = 0;

  /// initialize the element type
  virtual int Initialize(DRT::Discretization& dis);

  /// preevaluation
  virtual inline void PreEvaluate(DRT::Discretization& dis,
                                  Teuchos::ParameterList& p,
                                  Teuchos::RCP<LINALG::SparseOperator> systemmatrix1,
                                  Teuchos::RCP<LINALG::SparseOperator> systemmatrix2,
                                  Teuchos::RCP<Epetra_Vector>          systemvector1,
                                  Teuchos::RCP<Epetra_Vector>          systemvector2,
                                  Teuchos::RCP<Epetra_Vector>          systemvector3)
  { return; }

  /// nodal block information to create a null space description
  virtual void NodalBlockInformation( DRT::Element * dwele, int & numdf, int & dimns, int & nv, int & np ) = 0;

  /// do the null space computation
  virtual void ComputeNullSpace( DRT::Discretization & dis, std::vector<double> & ns, const double * x0, int numdf, int dimns ) = 0;

};

}

#endif
