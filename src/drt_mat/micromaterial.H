/*----------------------------------------------------------------------*/
/*!
\file micromaterial.H

\brief class for handling of micro-macro transitions

<pre>
Maintainer: Lena Wiechert
            wiechert@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15303
</pre>
*/
/*----------------------------------------------------------------------*/

#ifndef MICROMATERIAL_H
#define MICROMATERIAL_H

#ifdef D_SOH8
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE

#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <list>
#include <string>
#include <algorithm>

#include <Epetra_Comm.h>
#include <Teuchos_RefCountPtr.hpp>

#include "../drt_lib/drt_discret.H"
#include "material.H"

using namespace std;
using namespace Teuchos;

class MicroMaterial : public Material
{
public:

  /// construct an instance of MicroMaterial for a given Gauss point
  /// and microscale discretization
  MicroMaterial(int gp, int microdis_num);

  /// destructor
  ~MicroMaterial();


  /// test routine for calculating stresses, constitutive matrix and density in
  /// case of St Venant Kirchhoff material
  void CalcStressStiffDens (Epetra_SerialDenseVector* stress,
                            Epetra_SerialDenseMatrix* cmat,
                            double* density,
                            const Epetra_SerialDenseVector* glstrain);

private:

  /// my Gauss point
  int gp_;

  /// my micro discretization
  RefCountPtr<DRT::Discretization> microdis_;

  /// my vector of displacements
  RefCountPtr<Epetra_Vector> disp_;

};

#endif
#endif
#endif
#endif
