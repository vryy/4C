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
#include "../drt_lib/linalg_fixedsizematrix.H"


const static std::size_t IerrorCheckValue = 77777777;

//! for non-xfem elements, figure out the order of unknowns per node in the array
template<XFEM::PHYSICS::Field field>
struct FieldNumPatch {static const std::size_t pos = IerrorCheckValue;};
template<> struct FieldNumPatch<XFEM::PHYSICS::Dispx> {static const std::size_t pos = 0;};
template<> struct FieldNumPatch<XFEM::PHYSICS::Dispy> {static const std::size_t pos = 1;};
template<> struct FieldNumPatch<XFEM::PHYSICS::Dispz> {static const std::size_t pos = 2;};

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

      if (testfield == XFEM::PHYSICS::Dispx or
          testfield == XFEM::PHYSICS::Dispy or
          testfield == XFEM::PHYSICS::Dispz )
      {
        if (trialfield == XFEM::PHYSICS::Dispx or
            trialfield == XFEM::PHYSICS::Dispy or
            trialfield == XFEM::PHYSICS::Dispz )
        {
          Matrix_dd<testfield,trialfield>(emat, testshape, fac, trialshape);
        }
        else
        {
          Matrix_ds<testfield,trialfield>(emat, testshape, fac, trialshape);
        }
      }
      else
      {
        Matrix_sd<testfield,trialfield>(emat, testshape, fac, trialshape);
      }
    }

    //! assemble one expression for the element stiffness matrix for a given test and trial field
    template <XFEM::PHYSICS::Field testfield,
              XFEM::PHYSICS::Field trialfield,
              class V1,
              class V2>
    void Matrix_ds(
            Epetra_SerialDenseMatrix& emat,
            const V1&           testshape,
            const double&       fac,
            const V2&           trialshape
            ) const
    {

      const int numdof = 3;
//cout << emat << endl;
//        const std::vector<int>& testdof (dofman_.LocalDofPosPerField<SymmetryMangle<testfield>::name>());
        const std::vector<int>& trialdof(dofman_.LocalDofPosPerField<SymmetryMangle<trialfield>::name>());

        const std::size_t numparamtest = emat.M()/numdof;
//        cout << numparamtest << endl;
        const std::size_t numparamtrial = trialdof.size();


        for (std::size_t ui=0; ui<numparamtrial; ++ui)
        {
            const std::size_t trialpos = trialdof[ui];
//            cout << trialpos << endl;
            const double multiply = fac * trialshape(ui);
//            cout << fac << endl;
//            cout<< trialshape(ui) << endl;
            for (std::size_t vi=0; vi<numparamtest; ++vi)
            {
                const std::size_t testpos = numdof*vi + FieldNumPatch<testfield>::pos;
//                cout << testpos << ", " << trialpos << endl;
                emat(testpos, trialpos) += testshape(vi)*multiply ;
            }
        }
    }

    //! assemble one expression for the element stiffness matrix for a given test and trial field
    template <XFEM::PHYSICS::Field testfield,
              XFEM::PHYSICS::Field trialfield,
              class V1,
              class V2>
    void Matrix_sd(
            Epetra_SerialDenseMatrix& emat,
            const V1&           testshape,
            const double&       fac,
            const V2&           trialshape
            ) const
    {

        const int numdof = 3;

        const std::vector<int>& testdof (dofman_.LocalDofPosPerField<SymmetryMangle<testfield>::name>());
//        const std::vector<int>& trialdof(dofman_.LocalDofPosPerField<SymmetryMangle<trialfield>::name>());

        const std::size_t numparamtest  = testdof.size();
        const std::size_t numparamtrial = emat.N()/numdof;


        for (std::size_t ui=0; ui<numparamtrial; ++ui)
        {
            const std::size_t trialpos = numdof*ui + FieldNumPatch<testfield>::pos;
            const double multiply = fac * trialshape(ui);
            for (std::size_t vi=0; vi<numparamtest; ++vi)
            {
                const std::size_t testpos = testdof[vi];
                emat(testpos, trialpos) += testshape(vi)*multiply ;
            }
        }
    }

    //! assemble one expression for the element stiffness matrix for a given test and trial field
    template <XFEM::PHYSICS::Field testfield,
              XFEM::PHYSICS::Field trialfield,
              class V1,
              class V2>
    void Matrix_dd(
            Epetra_SerialDenseMatrix& emat,
            const V1&           testshape,
            const double&       fac,
            const V2&           trialshape
            ) const
    {
        const int numdof = 3;

        const std::size_t numparamtest = emat.M()/numdof;
        const std::size_t numparamtrial = emat.N()/numdof;
        dsassert(numparamtest==numparamtrial, "bug. matrix not quadratic!");

        for (std::size_t ui=0; ui<numparamtrial; ++ui)
        {
            const std::size_t trialpos = numdof*ui + FieldNumPatch<testfield>::pos;
            const double multiply = fac * trialshape(ui);
            for (std::size_t vi=0; vi<numparamtest; ++vi)
            {
                const std::size_t testpos = numdof*vi + FieldNumPatch<testfield>::pos;
                emat(testpos, trialpos) += testshape(vi)*multiply ;
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
        const int numdof = 3;

        const std::size_t numparamtest = evec.Length()/numdof;
        for (std::size_t vi=0; vi<numparamtest; ++vi)
        {
            const std::size_t testpos = numdof*vi + FieldNumPatch<testfield>::pos;
            evec(testpos) += fac*testshape(vi);
        }

    }

private:
    // don't want = operator and cctor
    IFacePatchLocalAssembler operator = (const IFacePatchLocalAssembler& old);
    IFacePatchLocalAssembler(const IFacePatchLocalAssembler& old);

    const XFEM::ElementDofManager& dofman_;
//    Epetra_SerialDenseMatrix& estif_;
//    Epetra_SerialDenseVector& eforce_;
};


#endif
#endif
