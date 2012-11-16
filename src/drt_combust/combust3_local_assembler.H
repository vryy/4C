/*----------------------------------------------------------------------*/
/*!
\file combust3_local_assembler.H

\brief routines to assemble element stiffness matrix and vector

<pre>
Maintainer: Florian Henke
            henke@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15265
</pre>
*/
/*----------------------------------------------------------------------*/

#ifndef COMBUST3_LOCAL_ASSEMBLER_H
#define COMBUST3_LOCAL_ASSEMBLER_H


#include "combust_defines.H"
#include "../drt_xfem/xfem_utils.H"
#include "../drt_xfem/physics.H"
#include "../drt_xfem/dof_management_element.H"
#include "../drt_xfem/xfem_element_utils.H"


namespace COMBUST
{

const static std::size_t errorCheckValue = 77777777;

//! for non-xfem elements, figure out the order of unknowns per node in the array
template<XFEM::PHYSICS::Field field>
struct FieldNum {static const std::size_t pos = errorCheckValue;};
template<> struct FieldNum<XFEM::PHYSICS::Velx> {static const std::size_t pos = 0;};
template<> struct FieldNum<XFEM::PHYSICS::Vely> {static const std::size_t pos = 1;};
template<> struct FieldNum<XFEM::PHYSICS::Velz> {static const std::size_t pos = 2;};
template<> struct FieldNum<XFEM::PHYSICS::Pres> {static const std::size_t pos = 3;};

//! for the viscous stress, replace shear stress components with symmetric part
template <XFEM::PHYSICS::Field FIELD>
struct SymmetryMangle {static const XFEM::PHYSICS::Field name = FIELD;};
template<> struct SymmetryMangle<XFEM::PHYSICS::Sigmayx> {static const XFEM::PHYSICS::Field name = XFEM::PHYSICS::Sigmaxy;};
template<> struct SymmetryMangle<XFEM::PHYSICS::Sigmazx> {static const XFEM::PHYSICS::Field name = XFEM::PHYSICS::Sigmaxz;};
template<> struct SymmetryMangle<XFEM::PHYSICS::Sigmazy> {static const XFEM::PHYSICS::Field name = XFEM::PHYSICS::Sigmayz;};
template<> struct SymmetryMangle<XFEM::PHYSICS::Epsilonyx>   {static const XFEM::PHYSICS::Field name = XFEM::PHYSICS::Epsilonxy;};
template<> struct SymmetryMangle<XFEM::PHYSICS::Epsilonzx>   {static const XFEM::PHYSICS::Field name = XFEM::PHYSICS::Epsilonxz;};
template<> struct SymmetryMangle<XFEM::PHYSICS::Epsilonzy>   {static const XFEM::PHYSICS::Field name = XFEM::PHYSICS::Epsilonyz;};


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
      std::vector<int>::const_iterator trialpos = trialdof.begin();
      for (std::size_t ui=0; ui!=numparamtrial; ++ui)
      {
        const double multiply = fac * trialshape(ui);
        std::vector<int>::const_iterator testpos = testdof.begin();
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
      std::vector<int>::const_iterator testpos = testdof.begin();
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
        //------------------------------
        // short cut for diagonal blocks
        //------------------------------
        if (numparamtest == numparamtrial) // diagonal blocks
        {
          switch (numparamtest)
          {
            case  8: this->template IntegrateMatrix< 8, 8>(testshape,testdof,fac,trialshape,trialdof); break;
            case 20: this->template IntegrateMatrix<20,20>(testshape,testdof,fac,trialshape,trialdof); break;
            case 27: this->template IntegrateMatrix<27,27>(testshape,testdof,fac,trialshape,trialdof); break;
            default: found_template = false;
          }
        }
        else
        {
          found_template = false;
        }

        //----------------------
        // this the regular path
        //----------------------
        if (not found_template)
        {
          // this loop does the same thing as the one below, but it is slower
          //cout << numparamtest << "   " << numparamtrial << endl;
          //for (std::size_t ui=0; ui<numparamtrial; ++ui)
          //{
          //  const std::size_t trialpos = trialdof[ui];
          //  const double multiply = fac * trialshape(ui);
          //  for (std::size_t vi=0; vi<numparamtest; ++vi)
          //  {
          //    const std::size_t testpos = testdof[vi];
          //    estif_(testpos, trialpos) += testshape(vi)*multiply ;
          //  }
          //}
        // this loop does the same thing as the one above (commented), but it is faster
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

    //! assemble one expression for the element stiffness matrix for a given test and trial field
    template <XFEM::PHYSICS::Field testfield,
              XFEM::PHYSICS::Field trialfield,
              size_t shpVecSize>
    void MatrixSN(
        const LINALG::Matrix<shpVecSize,1>& testshape,
        const double&                       fac,
        const XFEM::ApproxFuncNormalPair<DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement>& trialshape
    ) const
    {
#ifndef COMBUST_NORMAL_ENRICHMENT
      dserror("wrong function!");
#endif
      Matrix<testfield,trialfield>         (testshape, fac, trialshape.s);
      Matrix<testfield,XFEM::PHYSICS::Veln>(testshape, fac, trialshape.n);
    }

    //! assemble one expression for the element stiffness matrix for a given test and trial field
    template <XFEM::PHYSICS::Field testfield,
              XFEM::PHYSICS::Field trialfield,
              size_t shpVecSize>
    void MatrixNS(
        const XFEM::ApproxFuncNormalPair<DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement>& testshape,
        const double&                       fac,
        const LINALG::Matrix<shpVecSize,1>& trialshape
    ) const
    {
#ifndef COMBUST_NORMAL_ENRICHMENT
      dserror("wrong function!");
#endif
      Matrix<testfield          ,trialfield> (testshape.s, fac, trialshape);
      Matrix<XFEM::PHYSICS::Veln,trialfield> (testshape.n, fac, trialshape);
    }

    //! assemble one expression for the element stiffness matrix for a given test and trial field
    template <XFEM::PHYSICS::Field testfield,
              XFEM::PHYSICS::Field trialfield>
    void MatrixNN(
        const XFEM::ApproxFuncNormalPair<DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement>& testshape,
        const double& fac,
        const XFEM::ApproxFuncNormalPair<DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement>& trialshape
    ) const
    {
#ifndef COMBUST_NORMAL_ENRICHMENT
      dserror("wrong function!");
#endif
      //if(numparamtest != numparamtrial) dserror("so nicht");

      //cout << physVarToString(testfield) << endl;

      //for (int i=0;i<numparamtest;i++)
      //{
      //  cout << "testparam " << i << " " << testdof[i] << endl;
      //}
      //cout << physVarToString(trialfield) << endl;
      //for (int i=0;i<numparamtrial;i++)
      //{
      //  cout << "trialparam " << i << " " << trialdof[i] << endl;
      //}

      Matrix<testfield          ,trialfield>         (testshape.s, fac, trialshape.s);
      Matrix<testfield          ,XFEM::PHYSICS::Veln>(testshape.s, fac, trialshape.n);
      Matrix<XFEM::PHYSICS::Veln,trialfield>         (testshape.n, fac, trialshape.s);
      Matrix<XFEM::PHYSICS::Veln,XFEM::PHYSICS::Veln>(testshape.n, fac, trialshape.n);
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
          case 20: this->template IntegrateVector<20>(testshape,testdof,fac); break;
          case  8: this->template IntegrateVector< 8>(testshape,testdof,fac); break;
          default: found_template = false;
        }

        //for (std::size_t vi=0; vi<numparamtest; ++vi)
        //{
        //  const int testpos = testdof[vi];
        //  eforce_(testpos) += fac*testshape(vi);
        //}
        // this loop does the same thing as the one above (commented), but it is faster
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

    //! assemble one expression for the element rhs for a given test field
    template <XFEM::PHYSICS::Field testfield>
    void VectorN(
        const XFEM::ApproxFuncNormalPair<DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement>& testshape,
        const double& fac
    ) const
    {
#ifndef COMBUST_NORMAL_ENRICHMENT
      dserror("wrong function!");
#endif

      Vector<testfield>          (testshape.s, fac);
      Vector<XFEM::PHYSICS::Veln>(testshape.n, fac);
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
template <DRT::Element::DiscretizationType DISTYPE, std::size_t NUMDOF>
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

    //! assemble one expression for the element stiffness matrix for a given test and trial field
    template <XFEM::PHYSICS::Field testfield,
              XFEM::PHYSICS::Field trialfield,
              size_t shpVecSize>
    void MatrixNS(
        const XFEM::ApproxFuncNormalPair<DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement>& testshape,
        const double& fac,
        const LINALG::Matrix<shpVecSize,1>& trialshape
    )
    {
#ifndef COMBUST_NORMAL_ENRICHMENT
      dserror("wrong function!");
#endif
      Matrix<testfield,trialfield>(testshape.s,fac,trialshape);
    }

    //! assemble one expression for the element stiffness matrix for a given test and trial field
    template <XFEM::PHYSICS::Field testfield,
              XFEM::PHYSICS::Field trialfield,
              size_t shpVecSize>
    void MatrixSN(
        const LINALG::Matrix<shpVecSize,1>& testshape,
        const double& fac,
        const XFEM::ApproxFuncNormalPair<DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement>& trialshape
    )
    {
#ifndef COMBUST_NORMAL_ENRICHMENT
      dserror("wrong function!");
#endif
      Matrix<testfield,trialfield>(testshape,fac,trialshape.s);
    }

    //! assemble one expression for the element stiffness matrix for a given test and trial field
    template <XFEM::PHYSICS::Field testfield,
              XFEM::PHYSICS::Field trialfield>
    void MatrixNN(
        const XFEM::ApproxFuncNormalPair<DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement>& testshape,
        const double& fac,
        const XFEM::ApproxFuncNormalPair<DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement>& trialshape
    )
    {
#ifndef COMBUST_NORMAL_ENRICHMENT
      dserror("wrong function!");
#endif
      Matrix<testfield,trialfield>(testshape.s,fac,trialshape.s);
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

    //! assemble one expression for the element rhs for a given test field
    template <XFEM::PHYSICS::Field testfield>
    void VectorN(
            const XFEM::ApproxFuncNormalPair<DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement>& testshape,
            const double& fac
            )
    {
#ifndef COMBUST_NORMAL_ENRICHMENT
      dserror("wrong function!");
#endif
      Vector<testfield>(testshape.s,fac);
    }

private:
    // don't want = operator and cctor
    LocalAssembler operator = (const LocalAssembler& old);
    LocalAssembler(const LocalAssembler& old);

    LINALG::Matrix<NUMDOF*DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement,NUMDOF*DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement> estif_;
    LINALG::Matrix<NUMDOF*DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement,1> eforce_;
};

}

#endif
