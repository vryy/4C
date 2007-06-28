/*!----------------------------------------------------------------------
\file so_material.cpp
\brief Everything concerning structural material selection for solid elements

<pre>
Maintainer: Moritz Frenzel
            frenzel@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15240
</pre>

*----------------------------------------------------------------------*/
#ifdef D_SOH8
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE

// This is just here to get the c++ mpi header, otherwise it would
// use the c version included inside standardtypes.h
#ifdef PARALLEL
#include "mpi.h"
#endif
#include "so_hex8.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_exporter.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_lib/linalg_serialdensematrix.H"
#include "../drt_lib/linalg_serialdensevector.H"
#include "Epetra_SerialDenseSolver.h"

#include "../drt_mat/micromaterial.H"
#include "../drt_mat/stvenantkirchhoff.H"


//extern "C"
//{
//#include "../headers/standardtypes.h"
//}
using namespace std; // cout etc.
using namespace LINALG; // our linear algebra

/*----------------------------------------------------------------------*
 | material laws for So_hex8                                   maf 04/07|
 *----------------------------------------------------------------------*/
void DRT::Elements::So_hex8::soh8_mat_sel(
      Epetra_SerialDenseVector* stress,
      Epetra_SerialDenseMatrix* cmat,
      double* density,
      const Epetra_SerialDenseVector* glstrain,
      const Epetra_SerialDenseMatrix* defgrd,
      int gp)
{
  RefCountPtr<MAT::Material> mat = Material();
  switch (mat->MaterialType())
  {
    case m_stvenant: /*------------------ st.venant-kirchhoff-material */
    {
      MAT::StVenantKirchhoff* stvk = static_cast <MAT::StVenantKirchhoff*>(mat.get());
      
      stvk->Evaluate(glstrain,cmat,stress);
      
      *density = stvk->Density();
      
      break;
    }

#if 0
    case m_struct_multiscale: /*------------------- multiscale approach */
    {
      // Here macro-micro transition (localization) will take place

      int microdis_num = material->m.struct_multiscale->microdis;

      if (gp > static_cast<int>(mat_.size())-1)
      {
        mat_.resize(gp+1);
        mat_[gp] = rcp(new MAT::MicroMaterial(gp, microdis_num));
      }

      MAT::MicroMaterial* micromat =
        dynamic_cast<MAT::MicroMaterial*>(mat_[gp].get());

      if (micromat == NULL)
        dserror("Wrong type of derived material class");

      micromat->CalcStressStiffDens(stress, cmat, density, glstrain);
      break;
    }
#endif

    default:
      dserror("Illegal type %d of material for element solid3 hex8", mat->MaterialType());
      break;
  }

  /*--------------------------------------------------------------------*/
  return;
}  // of soh8_mat_sel

#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_SOH8
