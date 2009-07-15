/*----------------------------------------------------------------------*/
/*!
\file xfluid3_local_assembler.H

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

#ifndef XFLUID3_LOCAL_ASSEMBLER_H
#define XFLUID3_LOCAL_ASSEMBLER_H

#include "../drt_xfem/physics.H"
#include "../drt_xfem/dof_management_element.H"
#include <Epetra_SerialDenseMatrix.h>
#include <Epetra_SerialDenseVector.h>
#include "../drt_lib/linalg_fixedsizematrix.H"


const static std::size_t errorCheckValue = 77777777;

//! for non-xfem elements, figure out the order of unknowns per node in the array
template<XFEM::PHYSICS::Field field>
struct FieldNum {static const std::size_t pos = errorCheckValue;};
template<> struct FieldNum<XFEM::PHYSICS::Temp> {static const std::size_t pos = 0;};

//! for the viscous stress, replace shear stress components with symmetric part
template <XFEM::PHYSICS::Field FIELD>
struct SymmetryMangle {static const XFEM::PHYSICS::Field name = FIELD;};

////! which shape functions are equal anyway
//template <XFEM::PHYSICS::Field FIELD>
//struct ShapeMangle;
//template<> struct ShapeMangle<XFEM::PHYSICS::Velx> {static const XFEM::PHYSICS::Field name = XFEM::PHYSICS::Velx;};
//template<> struct ShapeMangle<XFEM::PHYSICS::Vely> {static const XFEM::PHYSICS::Field name = XFEM::PHYSICS::Velx;};
//template<> struct ShapeMangle<XFEM::PHYSICS::Velz> {static const XFEM::PHYSICS::Field name = XFEM::PHYSICS::Velx;};
//template<> struct ShapeMangle<XFEM::PHYSICS::Pres> {static const XFEM::PHYSICS::Field name = XFEM::PHYSICS::Velx;};
//template<> struct ShapeMangle<XFEM::PHYSICS::DiscPres> {static const XFEM::PHYSICS::Field name = XFEM::PHYSICS::DiscPres;};
//template<> struct ShapeMangle<XFEM::PHYSICS::Sigmaxx> {static const XFEM::PHYSICS::Field name = XFEM::PHYSICS::Sigmaxx;};
//template<> struct ShapeMangle<XFEM::PHYSICS::Sigmaxy> {static const XFEM::PHYSICS::Field name = XFEM::PHYSICS::Sigmaxx;};
//template<> struct ShapeMangle<XFEM::PHYSICS::Sigmaxz> {static const XFEM::PHYSICS::Field name = XFEM::PHYSICS::Sigmaxx;};
//template<> struct ShapeMangle<XFEM::PHYSICS::Sigmayx> {static const XFEM::PHYSICS::Field name = XFEM::PHYSICS::Sigmaxx;};
//template<> struct ShapeMangle<XFEM::PHYSICS::Sigmayy> {static const XFEM::PHYSICS::Field name = XFEM::PHYSICS::Sigmaxx;};
//template<> struct ShapeMangle<XFEM::PHYSICS::Sigmayz> {static const XFEM::PHYSICS::Field name = XFEM::PHYSICS::Sigmaxx;};
//template<> struct ShapeMangle<XFEM::PHYSICS::Sigmazx> {static const XFEM::PHYSICS::Field name = XFEM::PHYSICS::Sigmaxx;};
//template<> struct ShapeMangle<XFEM::PHYSICS::Sigmazy> {static const XFEM::PHYSICS::Field name = XFEM::PHYSICS::Sigmaxx;};
//template<> struct ShapeMangle<XFEM::PHYSICS::Sigmazz> {static const XFEM::PHYSICS::Field name = XFEM::PHYSICS::Sigmaxx;};
//template<> struct ShapeMangle<XFEM::PHYSICS::Tauxx>   {static const XFEM::PHYSICS::Field name = XFEM::PHYSICS::Tauxx;};
//template<> struct ShapeMangle<XFEM::PHYSICS::Tauxy>   {static const XFEM::PHYSICS::Field name = XFEM::PHYSICS::Tauxx;};
//template<> struct ShapeMangle<XFEM::PHYSICS::Tauxz>   {static const XFEM::PHYSICS::Field name = XFEM::PHYSICS::Tauxx;};
//template<> struct ShapeMangle<XFEM::PHYSICS::Tauyx>   {static const XFEM::PHYSICS::Field name = XFEM::PHYSICS::Tauxx;};
//template<> struct ShapeMangle<XFEM::PHYSICS::Tauyy>   {static const XFEM::PHYSICS::Field name = XFEM::PHYSICS::Tauxx;};
//template<> struct ShapeMangle<XFEM::PHYSICS::Tauyz>   {static const XFEM::PHYSICS::Field name = XFEM::PHYSICS::Tauxx;};
//template<> struct ShapeMangle<XFEM::PHYSICS::Tauzx>   {static const XFEM::PHYSICS::Field name = XFEM::PHYSICS::Tauxx;};
//template<> struct ShapeMangle<XFEM::PHYSICS::Tauzy>   {static const XFEM::PHYSICS::Field name = XFEM::PHYSICS::Tauxx;};
//template<> struct ShapeMangle<XFEM::PHYSICS::Tauzz>   {static const XFEM::PHYSICS::Field name = XFEM::PHYSICS::Tauxx;};

//! class that assembles one after another integral of the weak form into the stiffness matrix
template <DRT::Element::DiscretizationType DISTYPE,
          XFEM::AssemblyType ASSEMBLYTYPE,
          std::size_t NUMDOF>
class LocalAssembler
{
public:

    explicit LocalAssembler(
            const XFEM::ElementDofManager&    dofman,
            Epetra_SerialDenseMatrix& estif,
            Epetra_SerialDenseVector& eforce
            ) :
        dofman_(dofman),
        estif_(estif),
        eforce_(eforce)
    {
        return;
    }

    ~LocalAssembler()
    {
        return;
    }

private:
    template<std::size_t numparamtest,
    std::size_t numparamtrial,
             class V1,
             class V2>
    void IntegrateMatrix(
        const V1&                 testshape,
        const std::vector<int>&   testdof,
        const double&             fac,
        const V2&                 trialshape,
        const std::vector<int>&   trialdof
        ) const
    {
      vector<int>::const_iterator trialpos = trialdof.begin();
      for (std::size_t ui=0; ui!=numparamtrial; ++ui)
      {
        const double multiply = fac * trialshape(ui);
        vector<int>::const_iterator testpos = testdof.begin();
        for (std::size_t vi=0; vi!=numparamtest; ++vi)
        {
//          const std::size_t testpos = testdof[vi];
          estif_(*testpos, *trialpos) += testshape(vi)*multiply ;
          ++testpos;
        }
        ++trialpos;
      }
    }

    template<std::size_t numparamtest,
             class V1>
    void IntegrateVector(
        const V1&                 testshape,
        const std::vector<int>&   testdof,
        const double&             fac
        ) const
    {
      vector<int>::const_iterator testpos = testdof.begin();
      for (std::size_t vi=0; vi!=numparamtest; ++vi)
      {
          eforce_(*testpos) += fac*testshape(vi);
          ++testpos;
      }
    }

public:
    //! assemble one expression for the element stiffness matrix for a given test and trial field
    template <XFEM::PHYSICS::Field testfield,
              XFEM::PHYSICS::Field trialfield,
              class V1,
              class V2>
    void Matrix(
            const V1&           testshape,
            const double&       fac,
            const V2&           trialshape
            ) const
    {

        const std::vector<int>& testdof (dofman_.LocalDofPosPerField<SymmetryMangle<testfield>::name>());
        const std::vector<int>& trialdof(dofman_.LocalDofPosPerField<SymmetryMangle<trialfield>::name>());

        const std::size_t numparamtest  = testdof.size();
        const std::size_t numparamtrial = trialdof.size();

        bool found_template = true;
        switch (numparamtest)
        {
        case 20:
        {
          switch (numparamtrial) {
          case 20: this->template IntegrateMatrix<20,20>(testshape,testdof,fac,trialshape,trialdof); break;
          default: found_template = false;
          }
          break;
        }
        case 8:
        {
          switch (numparamtrial) {
          case 8: this->template IntegrateMatrix<8,8>(testshape,testdof,fac,trialshape,trialdof); break;
          default: found_template = false;
          }
          break;
        }
        default:
          found_template = false;
        }

        if (not found_template)
        {
//          cout << numparamtest << "   " << numparamtrial << endl;
//          for (std::size_t ui=0; ui<numparamtrial; ++ui)
//          {
//              const std::size_t trialpos = trialdof[ui];
//              const double multiply = fac * trialshape(ui);
//              for (std::size_t vi=0; vi<numparamtest; ++vi)
//              {
//                  const std::size_t testpos = testdof[vi];
//                  estif_(testpos, trialpos) += testshape(vi)*multiply ;
//              }
//          }
        std::size_t ui = 0;
        for (std::vector<int>::const_iterator trialpos = trialdof.begin(); trialpos != trialdof.end(); ++trialpos)
        {
            const double multiply = fac * trialshape(ui);
            std::size_t vi = 0;
            for (std::vector<int>::const_iterator testpos = testdof.begin(); testpos != testdof.end(); ++testpos)
            {
                estif_(*testpos, *trialpos) += testshape(vi)*multiply;
                ++vi;
            }
            ++ui;
        }
        }
    }

    //! assemble one expression for the element rhs for a given test field
    template <XFEM::PHYSICS::Field testfield,
              class V1>
    void Vector(
            const V1&           testshape,
            const double&       fac
            ) const
    {
        const std::vector<int>& testdof(dofman_.LocalDofPosPerField<SymmetryMangle<testfield>::name>());

        const std::size_t numparamtest  = testdof.size();

        bool found_template = true;
        switch (numparamtest)
        {
        case 20:
        {
          this->template IntegrateVector<20>(testshape,testdof,fac); break;
        }
        case 8:
        {
          this->template IntegrateVector<8>(testshape,testdof,fac); break;
        }
        default:
          found_template = false;
        }
//
//        for (std::size_t vi=0; vi<numparamtest; ++vi)
//        {
//            const int testpos = testdof[vi];
//            eforce_(testpos) += fac*testshape(vi);
//        }
        if (not found_template)
        {
        std::size_t vi = 0;
        for (std::vector<int>::const_iterator testpos = testdof.begin(); testpos != testdof.end(); ++testpos)
        {
            eforce_(*testpos) += fac*testshape(vi);
            vi++;
        }
        }
    }

private:
    // don't want = operator and cctor
    LocalAssembler operator = (const LocalAssembler& old);
    LocalAssembler(const LocalAssembler& old);

    const XFEM::ElementDofManager& dofman_;
    Epetra_SerialDenseMatrix& estif_;
    Epetra_SerialDenseVector& eforce_;
};

/*!
 * class that assembles one after another integral of the weak form into the stiffness matrix
 *
 * specialization for standard assembly
 */
template <DRT::Element::DiscretizationType DISTYPE, int NUMDOF>
class LocalAssembler<DISTYPE, XFEM::standard_assembly, NUMDOF>
{
public:

    explicit LocalAssembler(
            const XFEM::ElementDofManager&    ,
            Epetra_SerialDenseMatrix&                         estif,
            Epetra_SerialDenseVector&                         eforce
            ) :
        estif_(LINALG::Matrix<NUMDOF*DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement,NUMDOF*DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement>(estif,true)),
        eforce_(LINALG::Matrix<NUMDOF*DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement,1>(eforce,true))
    {
        return;
    }

    explicit LocalAssembler(
            Epetra_SerialDenseMatrix&                         estif,
            Epetra_SerialDenseVector&                         eforce
            ) :
              estif_(LINALG::Matrix<NUMDOF*DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement,NUMDOF*DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement>(estif,true)),
              eforce_(LINALG::Matrix<NUMDOF*DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement,1>(eforce,true))
    {
        return;
    }

    ~LocalAssembler()
    {
        return;
    }

    //! assemble one expression for the element stiffness matrix for a given test and trial field
    template <XFEM::PHYSICS::Field testfield,
              XFEM::PHYSICS::Field trialfield,
              class V1,
              class V2>
    void Matrix(
            const V1&           testshape,
            const double&       fac,
            const V2&           trialshape
            )
    {
#ifdef DEBUG
      if (FieldNum<trialfield>::pos == errorCheckValue or FieldNum<testfield>::pos == errorCheckValue)
        dserror("no stress allowed here");
#endif
      const std::size_t numparam = DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;
        for (std::size_t ui=0; ui!=numparam; ++ui)
        {
            const std::size_t trialpos = NUMDOF*ui + FieldNum<trialfield>::pos;
            const double multiply = fac*trialshape(ui);
            for (std::size_t vi=0; vi!=numparam; ++vi)
            {
                const std::size_t testpos = NUMDOF*vi + FieldNum<testfield>::pos;
                estif_(testpos, trialpos) += testshape(vi)*multiply ;
            }
        }
    }

    //! assemble one expression for the element rhs for a given test field
    template <XFEM::PHYSICS::Field testfield,
              class V1>
    void Vector(
            const V1&           testshape,
            const double&       fac
            )
    {
#ifdef DEBUG
      if (FieldNum<testfield>::pos == errorCheckValue)
        dserror("no stress allowed here");
#endif
        const std::size_t numparam = DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;
        for (std::size_t vi=0; vi<numparam; ++vi)
        {
            const std::size_t testpos = NUMDOF*vi + FieldNum<testfield>::pos;
            eforce_(testpos) += fac*testshape(vi);
        }
    }

private:
    // don't want = operator and cctor
    LocalAssembler operator = (const LocalAssembler& old);
    LocalAssembler(const LocalAssembler& old);

    LINALG::Matrix<NUMDOF*DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement,NUMDOF*DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement> estif_;
    LINALG::Matrix<NUMDOF*DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement,1> eforce_;
};


#endif

#endif
