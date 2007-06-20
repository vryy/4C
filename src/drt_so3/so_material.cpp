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


extern "C"
{
#include "../headers/standardtypes.h"
}
#include "../drt_lib/dstrc.H"
using namespace std; // cout etc.
using namespace LINALG; // our linear algebra

/*----------------------------------------------------------------------*
 | material laws for So_hex8                                   maf 04/07|
 *----------------------------------------------------------------------*/
void DRT::Elements::So_hex8::soh8_mat_sel(
      struct _MATERIAL* material,
      Epetra_SerialDenseVector* stress,
      Epetra_SerialDenseMatrix* cmat,
      double* density,
      const Epetra_SerialDenseVector* glstrain,
      const Epetra_SerialDenseMatrix* defgrd,
      int gp)
{
  DSTraceHelper dst("So_hex8::soh8_mat_sel");

  switch (material->mattyp)
  {
    case m_stvenant: /*------------------ st.venant-kirchhoff-material */
    {
      // get material parameters
      double Emod = material->m.stvenant->youngs;    // Young's modulus (modulus of elasticity)
      double nu = material->m.stvenant->possionratio;// Poisson's ratio (Querdehnzahl)
      (*density) = material->m.stvenant->density;    // density, returned to evaluate mass matrix

      /*--------------------------------------------------------------------*/
      /* isotropic elasticity tensor C in matrix notion */
      /*                       [ 1-nu     nu     nu |          0    0    0 ]
       *                       [        1-nu     nu |          0    0    0 ]
       *           E           [               1-nu |          0    0    0 ]
       *   C = --------------- [ ~~~~   ~~~~   ~~~~   ~~~~~~~~~~  ~~~  ~~~ ]
       *       (1+nu)*(1-2*nu) [                    | (1-2*nu)/2    0    0 ]
       *                       [                    |      (1-2*nu)/2    0 ]
       *                       [ symmetric          |           (1-2*nu)/2 ]
       */
      double mfac = Emod/((1.0+nu)*(1.0-2.0*nu));  /* factor */
      /* write non-zero components */
      (*cmat)(0,0) = mfac*(1.0-nu);
      (*cmat)(0,1) = mfac*nu;
      (*cmat)(0,2) = mfac*nu;
      (*cmat)(1,0) = mfac*nu;
      (*cmat)(1,1) = mfac*(1.0-nu);
      (*cmat)(1,2) = mfac*nu;
      (*cmat)(2,0) = mfac*nu;
      (*cmat)(2,1) = mfac*nu;
      (*cmat)(2,2) = mfac*(1.0-nu);
      /* ~~~ */
      (*cmat)(3,3) = mfac*0.5*(1.0-2.0*nu);
      (*cmat)(4,4) = mfac*0.5*(1.0-2.0*nu);
      (*cmat)(5,5) = mfac*0.5*(1.0-2.0*nu);

      // evaluate stresses
      (*cmat).Multiply('N',(*glstrain),(*stress));   // sigma = C . epsilon
      break;
    }

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

    default:
      dserror("Illegal type of material for element solid3 hex8");
      break;
  }

  /*--------------------------------------------------------------------*/
  return;
}  // of soh8_mat_sel

#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_SOH8
