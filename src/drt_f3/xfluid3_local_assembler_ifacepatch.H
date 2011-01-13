/*----------------------------------------------------------------------*/
/*!
\file xfluid3_local_assembler_ifacepatch.H

\brief routines to assemble element stiffness matrix and vector

<pre>
Maintainer: Axel Gerstenberger
            gerstenberger@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>
*/
/*----------------------------------------------------------------------*/

#ifdef CCADISCRET

#ifndef XFLUID3_LOCAL_ASSEMBLER_IFACEPATCH_H
#define XFLUID3_LOCAL_ASSEMBLER_IFACEPATCH_H

#include "../drt_xfem/physics.H"
#include "../drt_xfem/dof_management_element.H"
#include <Epetra_SerialDenseMatrix.h>
#include <Epetra_SerialDenseVector.h>
#include "../linalg/linalg_fixedsizematrix.H"

namespace XFLUID
{

const static std::size_t IerrorCheckValue = 77777777;

//! for non-xfem elements, figure out the order of unknowns per node in the array

template<XFEM::PHYSICS::Field field>
struct IfaceFieldNum {static const std::size_t pos = IerrorCheckValue;};
template<> struct IfaceFieldNum<XFEM::PHYSICS::Velxiface> {static const std::size_t pos = 0;};
template<> struct IfaceFieldNum<XFEM::PHYSICS::Velyiface> {static const std::size_t pos = 1;};
template<> struct IfaceFieldNum<XFEM::PHYSICS::Velziface> {static const std::size_t pos = 2;};
template<> struct IfaceFieldNum<XFEM::PHYSICS::Dispxiface> {static const std::size_t pos = 0;};
template<> struct IfaceFieldNum<XFEM::PHYSICS::Dispyiface> {static const std::size_t pos = 1;};
template<> struct IfaceFieldNum<XFEM::PHYSICS::Dispziface> {static const std::size_t pos = 2;};

//! for the viscous stress, replace shear stress components with symmetric part
template <XFEM::PHYSICS::Field FIELD>
struct SymmetryManglePatch {static const XFEM::PHYSICS::Field name = FIELD;};
template<> struct SymmetryManglePatch<XFEM::PHYSICS::Sigmayx> {static const XFEM::PHYSICS::Field name = XFEM::PHYSICS::Sigmaxy;};
template<> struct SymmetryManglePatch<XFEM::PHYSICS::Sigmazx> {static const XFEM::PHYSICS::Field name = XFEM::PHYSICS::Sigmaxz;};
template<> struct SymmetryManglePatch<XFEM::PHYSICS::Sigmazy> {static const XFEM::PHYSICS::Field name = XFEM::PHYSICS::Sigmayz;};
template<> struct SymmetryManglePatch<XFEM::PHYSICS::Tauyx>   {static const XFEM::PHYSICS::Field name = XFEM::PHYSICS::Tauxy;};
template<> struct SymmetryManglePatch<XFEM::PHYSICS::Tauzx>   {static const XFEM::PHYSICS::Field name = XFEM::PHYSICS::Tauxz;};
template<> struct SymmetryManglePatch<XFEM::PHYSICS::Tauzy>   {static const XFEM::PHYSICS::Field name = XFEM::PHYSICS::Tauyz;};


//! class that assembles one after another integral of the weak form into the stiffness matrix
template <DRT::Element::DiscretizationType DISTYPE,
          std::size_t NUMDOF>
class IFacePatchLocalAssembler
{
public:

    explicit IFacePatchLocalAssembler(
            const XFEM::ElementDofManager&    dofman
            ) :
        dofman_(dofman)
    {
        return;
    }

    ~IFacePatchLocalAssembler()
    {
        return;
    }


public:
    //! assemble one expression for the element stiffness matrix for a given test and trial field
    template <XFEM::PHYSICS::Field testfield,
              XFEM::PHYSICS::Field trialfield,
              class V1,
              class V2>
    void Matrix(
            Epetra_SerialDenseMatrix& emat,
            const V1&           testshape,
            const double&       fac,
            const V2&           trialshape
            ) const
    {

      if (testfield == XFEM::PHYSICS::Velxiface or
          testfield == XFEM::PHYSICS::Velyiface or
          testfield == XFEM::PHYSICS::Velziface or
          testfield == XFEM::PHYSICS::Dispxiface or
          testfield == XFEM::PHYSICS::Dispyiface or
          testfield == XFEM::PHYSICS::Dispziface )
      {
        if (trialfield == XFEM::PHYSICS::Velx or
            trialfield == XFEM::PHYSICS::Vely or
            trialfield == XFEM::PHYSICS::Velz or
            trialfield == XFEM::PHYSICS::Sigmaxx or
            trialfield == XFEM::PHYSICS::Sigmaxy or
            trialfield == XFEM::PHYSICS::Sigmaxz or
            trialfield == XFEM::PHYSICS::Sigmayx or
            trialfield == XFEM::PHYSICS::Sigmayy or
            trialfield == XFEM::PHYSICS::Sigmayz or
            trialfield == XFEM::PHYSICS::Sigmazx or
            trialfield == XFEM::PHYSICS::Sigmazy or
            trialfield == XFEM::PHYSICS::Sigmazz )
        {
          Matrix_Guis_uncond<testfield,trialfield>(emat, testshape, fac, trialshape);
        }
        else if (trialfield == XFEM::PHYSICS::Velxiface or
                 trialfield == XFEM::PHYSICS::Velyiface or
                 trialfield == XFEM::PHYSICS::Velziface or
                 trialfield == XFEM::PHYSICS::Dispxiface or
                 trialfield == XFEM::PHYSICS::Dispyiface or
                 trialfield == XFEM::PHYSICS::Dispziface )
        {
          Matrix_Guiui_uncond<testfield,trialfield>(emat, testshape, fac, trialshape);
        }
        else
        {
          dserror(" Wrong trial position! ");
        }
      }
      else if (testfield == XFEM::PHYSICS::Velx or
               testfield == XFEM::PHYSICS::Vely or
               testfield == XFEM::PHYSICS::Velz or
               testfield == XFEM::PHYSICS::Sigmaxx or
               testfield == XFEM::PHYSICS::Sigmaxy or
               testfield == XFEM::PHYSICS::Sigmaxz or
               testfield == XFEM::PHYSICS::Sigmayx or
               testfield == XFEM::PHYSICS::Sigmayy or
               testfield == XFEM::PHYSICS::Sigmayz or
               testfield == XFEM::PHYSICS::Sigmazx or
               testfield == XFEM::PHYSICS::Sigmazy or
               testfield == XFEM::PHYSICS::Sigmazz )
      {
        if (trialfield == XFEM::PHYSICS::Velxiface or
            trialfield == XFEM::PHYSICS::Velyiface or
            trialfield == XFEM::PHYSICS::Velziface or
            trialfield == XFEM::PHYSICS::Dispxiface or
            trialfield == XFEM::PHYSICS::Dispyiface or
            trialfield == XFEM::PHYSICS::Dispziface )
        {
          Matrix_Gsui_uncond<testfield,trialfield>(emat, testshape, fac, trialshape);
        }
        else
        {
          dserror(" Wrong trial position! ");
        }
      }
      else
        dserror(" Wrong test position! ");
    }

    template <XFEM::PHYSICS::Field testfield,
              XFEM::PHYSICS::Field trialfield,
              class V1,
              class V2>
    void Matrix_Guis_uncond(
            Epetra_SerialDenseMatrix& emat,
            const V1&           testshape,
            const double&       fac,
            const V2&           trialshape
            ) const
    {

      const std::size_t numdof = 3;
      const std::vector<int>& trialdof(dofman_.LocalDofPosPerField<SymmetryMangle<trialfield>::name>());

      const std::size_t numnodetest = emat.M()/numdof;
      const std::size_t numparamtrial = trialdof.size();

      for (std::size_t ui=0; ui<numparamtrial; ++ui)
      {
        const std::size_t trialpos = trialdof[ui];
        const double multiply = fac * trialshape(ui);
        for (std::size_t vi=0; vi<numnodetest; ++vi)
        {
          const std::size_t testpos = numdof*vi + IfaceFieldNum<testfield>::pos;
          emat(testpos, trialpos) += testshape(vi)*multiply;
        }
      }
    }

    template <XFEM::PHYSICS::Field testfield,
                  XFEM::PHYSICS::Field trialfield,
                  class V1,
                  class V2>
    void Matrix_Gsui_uncond(
            Epetra_SerialDenseMatrix& emat,
            const V1&           testshape,
            const double&       fac,
            const V2&           trialshape
            ) const
    {

      const std::size_t numdof = 3;

      const std::vector<int>& testdof (dofman_.LocalDofPosPerField<SymmetryMangle<testfield>::name>());

      const std::size_t numparamtest = testdof.size();
      const std::size_t numnodetrial = emat.N()/numdof;
      for (std::size_t uinode=0; uinode<numnodetrial; ++uinode)
      {
        const std::size_t trialpos = numdof*uinode + IfaceFieldNum<trialfield>::pos;
        const double multiply = fac * trialshape(uinode);
        for (std::size_t vi=0; vi<numparamtest; ++vi)
        {
          const std::size_t testpos = testdof[vi];
          emat(testpos, trialpos) += testshape(vi)*multiply;
        }
      }
    }

    template <XFEM::PHYSICS::Field testfield,
              XFEM::PHYSICS::Field trialfield,
              class V1,
              class V2>
    void Matrix_Guiui_uncond(
        Epetra_SerialDenseMatrix& emat,
        const V1&           testshape,
        const double&       fac,
        const V2&           trialshape
    ) const
    {
      const std::size_t numdof = 3;
      const std::size_t numnodetest = emat.M()/numdof;
      const std::size_t numnodetrial = emat.N()/numdof;
      for (std::size_t uinode=0; uinode<numnodetrial; ++uinode)
      {
        const std::size_t trialpos = numdof*uinode + IfaceFieldNum<trialfield>::pos;
        const double multiply = fac * trialshape(uinode);
        for (std::size_t vinode=0; vinode<numnodetest; ++vinode)
        {
          const std::size_t testpos = numdof*vinode + IfaceFieldNum<testfield>::pos;
          emat(testpos, trialpos) += testshape(vinode)*multiply;
        }
      }
    }

    //! assemble one expression for the element rhs for a given test field
    template <XFEM::PHYSICS::Field testfield,
              class V1>
    void Vector(
            Epetra_SerialDenseVector& evec,
            const V1&           testshape,
            const double&       fac
            ) const
    {
      if (testfield == XFEM::PHYSICS::Velxiface or
          testfield == XFEM::PHYSICS::Velyiface or
          testfield == XFEM::PHYSICS::Velziface)
      {
        Vector_rhsui<testfield>(evec, testshape, fac);
      }
      else
        dserror(" Wrong position! ");
    }

    template <XFEM::PHYSICS::Field testfield,
                  class V1>
    void Vector_rhsui(
        Epetra_SerialDenseVector& evec,
        const V1&           testshape,
        const double&       fac
        ) const
    {
        const std::size_t numdof = 3;

        const std::size_t numnodetest = evec.Length()/numdof;
        for (std::size_t vi=0; vi<numnodetest; ++vi)
        {
            const std::size_t testpos = numdof*vi + IfaceFieldNum<testfield>::pos;
            evec(testpos) += fac*testshape(vi);
        }

    }

private:
    // don't want = operator and cctor
    IFacePatchLocalAssembler operator = (const IFacePatchLocalAssembler& old);
    IFacePatchLocalAssembler(const IFacePatchLocalAssembler& old);

    const XFEM::ElementDofManager& dofman_;
};

}

#endif
#endif
