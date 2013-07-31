/*----------------------------------------------------------------------*/
/*!
\file combust3_local_assembler.H

\brief routines to assemble element stiffness matrix and vector

<pre>
Maintainer: Ursula Rasthofer
            rasthofer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>
*/
/*----------------------------------------------------------------------*/

#ifndef COMBUST3_LOCAL_ASSEMBLER_H
#define COMBUST3_LOCAL_ASSEMBLER_H


#include "combust_defines.H"
#include "../drt_inpar/inpar_fluid.H"
#include "../drt_xfem/xfem_utils.H"
#include "../drt_xfem/physics.H"
#include "../drt_xfem/dof_management_element.H"
#include "../drt_xfem/xfem_element_utils.H"


namespace COMBUST
{

namespace UTILS
{
  enum patch_side
  {
    master,
    slave
  };
}

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
          //std::cout << numparamtest << "   " << numparamtrial << std::endl;
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

      //std::cout << physVarToString(testfield) << std::endl;

      //for (int i=0;i<numparamtest;i++)
      //{
      //  std::cout << "testparam " << i << " " << testdof[i] << std::endl;
      //}
      //std::cout << physVarToString(trialfield) << std::endl;
      //for (int i=0;i<numparamtrial;i++)
      //{
      //  std::cout << "trialparam " << i << " " << trialdof[i] << std::endl;
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


//! class that assembles one after another integral of the weak form into the stiffness matrix
//! patch-integration version
template <DRT::Element::DiscretizationType M_DISTYPE,
          DRT::Element::DiscretizationType S_DISTYPE,
          XFEM::AssemblyType M_ASSEMBLYTYPE,
          XFEM::AssemblyType S_ASSEMBLYTYPE,
          std::size_t NUMDOF>
class LocalPatchAssembler
{
public:

    explicit LocalPatchAssembler(
            const XFEM::ElementDofManager&   m_dofman,
            const XFEM::ElementDofManager&   s_dofman,
            std::vector<Epetra_SerialDenseMatrix>& estif,
            std::vector<Epetra_SerialDenseVector>& eforce,
            const INPAR::FLUID::EOS_GP_Pattern& pattern,
            const std::vector<int>& lm_masterNodeToPatch,
            const std::vector<int>& lm_slaveNodeToPatch
            ) :
        master_dofman_(m_dofman),
        slave_dofman_(s_dofman),
        estif_block_(estif),
        eforce_block_(eforce),
        pattern_(pattern),
        lm_masterNodeToPatch_(lm_masterNodeToPatch),
        lm_slaveNodeToPatch_(lm_slaveNodeToPatch)
    {
      // all values are automatically initialized to zero
      elematrix_mm_ = Teuchos::rcp(new Epetra_SerialDenseMatrix);
      elematrix_ms_ = Teuchos::rcp(new Epetra_SerialDenseMatrix);
      elematrix_sm_ = Teuchos::rcp(new Epetra_SerialDenseMatrix);
      elematrix_ss_ = Teuchos::rcp(new Epetra_SerialDenseMatrix);

      elevector_m_ = Teuchos::rcp(new Epetra_SerialDenseVector);
      elevector_s_ = Teuchos::rcp(new Epetra_SerialDenseVector);
      return;
    }

    ~LocalPatchAssembler()
    {
      return;
    }

private:

    template<std::size_t numparamtest,
             std::size_t numparamtrial,
             COMBUST::UTILS::patch_side testside,
             COMBUST::UTILS::patch_side trialside,
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
          switch (testside)
          {
           case COMBUST::UTILS::master:
           {
             switch (trialside)
             {
               case COMBUST::UTILS::master:
               {
                 (*elematrix_mm_)(*testpos, *trialpos) += testshape(vi)*multiply ;
                 break;
               }
               case COMBUST::UTILS::slave:
               {
                 (*elematrix_ms_)(*testpos, *trialpos) += testshape(vi)*multiply ;
                 break;
               }
               default:
               {
                 dserror("Unknown patch side");
                 break;
               }
             }
             break;
           }
           case COMBUST::UTILS::slave:
           {
             switch (trialside)
             {
               case COMBUST::UTILS::master:
               {
                 (*elematrix_sm_)(*testpos, *trialpos) += testshape(vi)*multiply ;
                 break;
               }
               case COMBUST::UTILS::slave:
               {
                 (*elematrix_ss_)(*testpos, *trialpos) += testshape(vi)*multiply ;
                 break;
               }
               default:
               {
                 dserror("Unknown patch side");
                 break;
               }
             }
             break;
           }
           default:
           {
             dserror("Unknown patch side");
             break;
           }
          }
          ++testpos;
        }
        ++trialpos;
      }
    }

    template<std::size_t numparamtest,
             COMBUST::UTILS::patch_side testside,
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
        switch (testside)
        {
          case COMBUST::UTILS::master:
          {
            (*elevector_m_)(*testpos) += fac*testshape(vi);
            break;
          }
          case COMBUST::UTILS::slave:
          {
            (*elevector_s_)(*testpos) += fac*testshape(vi);
            break;
          }
          default:
          {
            dserror("Unknown patch side");
            break;
          }
        }
        ++testpos;
      }
    }

public:

    //! assemble one expression for the element stiffness matrix for a given test and trial field
    template <XFEM::PHYSICS::Field testfield,
              XFEM::PHYSICS::Field trialfield,
              COMBUST::UTILS::patch_side testside,
              COMBUST::UTILS::patch_side trialside,
              class V1,
              class V2>
    void Matrix(
            const V1&           testshape,
            const double&       fac,
            const V2&           trialshape
            ) const
    {
        std::vector<int> testdof;
        std::vector<int> trialdof;

        switch (testside)
        {
         case COMBUST::UTILS::master:
         {
           testdof = master_dofman_.LocalDofPosPerField<SymmetryMangle<testfield>::name>();
           break;
         }
         case COMBUST::UTILS::slave:
         {
           testdof = slave_dofman_.LocalDofPosPerField<SymmetryMangle<testfield>::name>();
           break;
         }
         default:
         {
           dserror("Unknown patch side!");
           break;
         }
        }
        switch (trialside)
        {
         case COMBUST::UTILS::master:
         {
           trialdof = master_dofman_.LocalDofPosPerField<SymmetryMangle<trialfield>::name>();
           break;
         }
         case COMBUST::UTILS::slave:
         {
           trialdof = slave_dofman_.LocalDofPosPerField<SymmetryMangle<trialfield>::name>();
           break;
         }
         default:
         {
           dserror("Unknown patch side!");
           break;
         }
        }

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
            case  8: this->template IntegrateMatrix< 8, 8,testside,trialside>(testshape,testdof,fac,trialshape,trialdof); break;
            case 20: this->template IntegrateMatrix<20,20,testside,trialside>(testshape,testdof,fac,trialshape,trialdof); break;
            case 27: this->template IntegrateMatrix<27,27,testside,trialside>(testshape,testdof,fac,trialshape,trialdof); break;
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
                switch (testside)
                          {
                           case COMBUST::UTILS::master:
                           {
                             switch (trialside)
                             {
                               case COMBUST::UTILS::master:
                               {
                                 (*elematrix_mm_)(*testpos, *trialpos) += testshape(vi)*multiply ;
                                 break;
                               }
                               case COMBUST::UTILS::slave:
                               {
                                 (*elematrix_ms_)(*testpos, *trialpos) += testshape(vi)*multiply ;
                                 break;
                               }
                               default:
                               {
                                 dserror("Unknown patch side");
                                 break;
                               }
                             }
                             break;
                           }
                           case COMBUST::UTILS::slave:
                           {
                             switch (trialside)
                             {
                               case COMBUST::UTILS::master:
                               {
                                 (*elematrix_sm_)(*testpos, *trialpos) += testshape(vi)*multiply ;
                                 break;
                               }
                               case COMBUST::UTILS::slave:
                               {
                                 (*elematrix_ss_)(*testpos, *trialpos) += testshape(vi)*multiply ;
                                 break;
                               }
                               default:
                               {
                                 dserror("Unknown patch side");
                                 break;
                               }
                             }
                             break;
                           }
                           default:
                           {
                             dserror("Unknown patch side");
                             break;
                           }
                          }
                ++vi;
            }
            ++ui;
        }
      }

        return;
    }

    //! assemble one expression for the element rhs for a given test field
    template <XFEM::PHYSICS::Field testfield,
              COMBUST::UTILS::patch_side testside,
              class V1>
    void Vector(
            const V1&           testshape,
            const double&       fac
            ) const
    {
        std::vector<int> testdof;
        switch (testside)
        {
          case COMBUST::UTILS::master:
          {
            testdof = master_dofman_.LocalDofPosPerField<SymmetryMangle<testfield>::name>();
            break;
          }
          case COMBUST::UTILS::slave:
          {
            testdof = slave_dofman_.LocalDofPosPerField<SymmetryMangle<testfield>::name>();
            break;
          }
          default:
          {
            dserror("Unknown patch side");
            break;
          }
        }

        const std::size_t numparamtest  = testdof.size();

        bool found_template = true;
        switch (numparamtest)
        {
          case 20: this->template IntegrateVector<20,testside>(testshape,testdof,fac); break;
          case  8: this->template IntegrateVector< 8,testside>(testshape,testdof,fac); break;
          default: found_template = false;
        }

        if (not found_template)
        {
        std::size_t vi = 0;
        for (std::vector<int>::const_iterator testpos = testdof.begin(); testpos != testdof.end(); ++testpos)
        {
            switch (testside)
            {
              case COMBUST::UTILS::master:
              {
                (*elevector_m_)(*testpos) += fac*testshape(vi);
                break;
              }
              case COMBUST::UTILS::slave:
              {
                (*elevector_s_)(*testpos) += fac*testshape(vi);
                break;
              }
              default:
              {
                dserror("Unknown patch side");
                break;
              }
            }
            vi++;
        }
        }
        return;
    }

    void PrintMat()
    {
      dserror("Not yet implemented");
      return;
    }

    void PrintRhs()
    {
      dserror("Not yet implemented");
      return;
    }

    // assembles master, slave and coupled matrices into the patch matrix blocks
    void Finish()
    {
      // assemble all blocks
      dserror("Not yet implemented");

      return;
    }

private:

    // don't want = operator and cctor
    LocalPatchAssembler operator = (const LocalPatchAssembler& old);
    LocalPatchAssembler(const LocalPatchAssembler& old);

    const XFEM::ElementDofManager& master_dofman_;
    const XFEM::ElementDofManager& slave_dofman_;
    std::vector<Epetra_SerialDenseMatrix>& estif_block_;
    std::vector<Epetra_SerialDenseVector>& eforce_block_;

    const INPAR::FLUID::EOS_GP_Pattern& pattern_;

    const std::vector<int>& lm_masterNodeToPatch_; ///< local map between master nodes and nodes in patch
    const std::vector<int>& lm_slaveNodeToPatch_;  ///< local map between slave nodes and nodes in patch

    Teuchos::RCP<Epetra_SerialDenseMatrix>  elematrix_mm_;         // element matrix master-master block
    Teuchos::RCP<Epetra_SerialDenseMatrix>  elematrix_ms_;         // element matrix master-slave block
    Teuchos::RCP<Epetra_SerialDenseMatrix>  elematrix_sm_;         // element matrix slave-master block
    Teuchos::RCP<Epetra_SerialDenseMatrix>  elematrix_ss_;         // element matrix slave-slave block

    Teuchos::RCP<Epetra_SerialDenseVector>  elevector_m_;          // element vector master block
    Teuchos::RCP<Epetra_SerialDenseVector>  elevector_s_;          // element vector slave block

};


//! class that assembles one after another integral of the weak form into the stiffness matrix
//! patch-integration version
//! specialization for standard assembly
template <DRT::Element::DiscretizationType M_DISTYPE,
          DRT::Element::DiscretizationType S_DISTYPE,
          std::size_t NUMDOF>
class LocalPatchAssembler<M_DISTYPE,S_DISTYPE,XFEM::standard_assembly,XFEM::standard_assembly,NUMDOF>
{
public:

    explicit LocalPatchAssembler(
            const XFEM::ElementDofManager&   ,
            const XFEM::ElementDofManager&   ,
            std::vector<Epetra_SerialDenseMatrix>& estif,
            std::vector<Epetra_SerialDenseVector>& eforce,
            const INPAR::FLUID::EOS_GP_Pattern& pattern,
            const std::vector<int>& lm_masterNodeToPatch,
            const std::vector<int>& lm_slaveNodeToPatch
            ) :
        estif_block_(estif),
        eforce_block_(eforce),
        pattern_(pattern),
        lm_masterNodeToPatch_(lm_masterNodeToPatch),
        lm_slaveNodeToPatch_(lm_slaveNodeToPatch),
        elematrix_mm_(true),
        elematrix_ms_(true),
        elematrix_sm_(true),
        elematrix_ss_(true),
        elevector_m_(true),
        elevector_s_(true)
    {
      return;
    }

    explicit LocalPatchAssembler(
            std::vector<Epetra_SerialDenseMatrix>& estif,
            std::vector<Epetra_SerialDenseVector>& eforce,
            const INPAR::FLUID::EOS_GP_Pattern& pattern,
            const std::vector<int>& lm_masterNodeToPatch,
            const std::vector<int>& lm_slaveNodeToPatch
            ) :
        estif_block_(estif),
        eforce_block_(eforce),
        pattern_(pattern),
        lm_masterNodeToPatch_(lm_masterNodeToPatch),
        lm_slaveNodeToPatch_(lm_slaveNodeToPatch),
        elematrix_mm_(true),
        elematrix_ms_(true),
        elematrix_sm_(true),
        elematrix_ss_(true),
        elevector_m_(true),
        elevector_s_(true)

    {
      return;
    }

    ~LocalPatchAssembler()
    {
      return;
    }

public:

    //! assemble one expression for the element stiffness matrix for a given test and trial field
    template <XFEM::PHYSICS::Field testfield,
              XFEM::PHYSICS::Field trialfield,
              COMBUST::UTILS::patch_side testside,
              COMBUST::UTILS::patch_side trialside,
              class V1,
              class V2>
    void Matrix(
            const V1&           testshape,
            const double&       fac,
            const V2&           trialshape
            )
    {
      std::size_t numparam_test = 0;
      switch (testside)
      {
        case COMBUST::UTILS::master:
        {
          numparam_test = DRT::UTILS::DisTypeToNumNodePerEle<M_DISTYPE>::numNodePerElement;
          break;
        }
        case COMBUST::UTILS::slave:
        {
          numparam_test = DRT::UTILS::DisTypeToNumNodePerEle<S_DISTYPE>::numNodePerElement;
          break;
        }
        default:
        {
          dserror("Unknown patch side");
          break;
        }
      }
      std::size_t numparam_trial = 0;
      switch (trialside)
      {
        case COMBUST::UTILS::master:
        {
          numparam_trial = DRT::UTILS::DisTypeToNumNodePerEle<M_DISTYPE>::numNodePerElement;
          break;
        }
        case COMBUST::UTILS::slave:
        {
          numparam_trial = DRT::UTILS::DisTypeToNumNodePerEle<S_DISTYPE>::numNodePerElement;
          break;
        }
        default:
        {
          dserror("Unknown patch side");
          break;
        }
      }

        for (std::size_t ui=0; ui!=numparam_trial; ++ui)
        {
            const std::size_t trialpos = NUMDOF*ui + FieldNum<trialfield>::pos;
            const double multiply = fac*trialshape(ui);
            for (std::size_t vi=0; vi!=numparam_test; ++vi)
            {
                const std::size_t testpos = NUMDOF*vi + FieldNum<testfield>::pos;
                switch (trialside)
                {
                  case COMBUST::UTILS::master:
                  {
                    switch (testside)
                    {
                      case COMBUST::UTILS::master:
                      {
                        elematrix_mm_(testpos, trialpos) += testshape(vi)*multiply ;
                        break;
                      }
                      case COMBUST::UTILS::slave:
                      {
                        elematrix_sm_(testpos, trialpos) += testshape(vi)*multiply ;
                        break;
                      }
                      default:
                      {
                        dserror("Unknown patch side");
                        break;
                      }
                    }
                    break;
                  }
                  case COMBUST::UTILS::slave:
                  {
                    switch (testside)
                    {
                      case COMBUST::UTILS::master:
                      {
                        elematrix_ms_(testpos, trialpos) += testshape(vi)*multiply ;
                        break;
                      }
                      case COMBUST::UTILS::slave:
                      {
                        elematrix_ss_(testpos, trialpos) += testshape(vi)*multiply ;
                        break;
                      }
                      default:
                      {
                        dserror("Unknown patch side");
                        break;
                      }
                    }
                    break;
                  }
                  default:
                  {
                    dserror("Unknown patch side");
                    break;
                  }
                }
            }
        }
    }

    //! assemble one expression for the element rhs for a given test field
    template <XFEM::PHYSICS::Field testfield,
              COMBUST::UTILS::patch_side testside,
              class V1>
    void Vector(
            const V1&           testshape,
            const double&       fac
            )
    {
        std::size_t numparam = 0;
        switch (testside)
        {
          case COMBUST::UTILS::master:
          {
            numparam = DRT::UTILS::DisTypeToNumNodePerEle<M_DISTYPE>::numNodePerElement;
            break;
          }
          case COMBUST::UTILS::slave:
          {
            numparam = DRT::UTILS::DisTypeToNumNodePerEle<S_DISTYPE>::numNodePerElement;
            break;
          }
          default:
          {
            dserror("Unknown patch side");
            break;
          }
        }

        for (std::size_t vi=0; vi<numparam; ++vi)
        {
          const std::size_t testpos = NUMDOF*vi + FieldNum<testfield>::pos;
          switch (testside)
          {
            case COMBUST::UTILS::master:
            {
              elevector_m_(testpos) += fac*testshape(vi);
              break;
            }
            case COMBUST::UTILS::slave:
            {
              elevector_s_(testpos) += fac*testshape(vi);
              break;
            }
            default:
            {
              dserror("Unknown patch side");
              break;
              }
            }
        }
    }

    void PrintMat()
    {
        std::cout << "elematrix_mm " << elematrix_mm_;
        std::cout << "elematrix_ms " << elematrix_ms_;
        std::cout << "elematrix_sm " << elematrix_sm_;
        std::cout << "elematrix_ss " << elematrix_ss_;
        std::cout<< "Hey" << std::endl;
     for (std::size_t i=0; i<estif_block_.size(); i++)
         estif_block_[i].Print(std::cout);
      return;
    }

    void PrintRhs()
    {
        std::cout << "elevector_m_ " << elevector_m_;
        std::cout << "elevector_s_ " << elevector_s_;
        std::cout<< "Hey" << std::endl;
        for (std::size_t i=0; i<eforce_block_.size(); i++)
            eforce_block_[i].Print(std::cout);
      return;
    }

    // assembles master, slave and coupled matrices into the patch matrix blocks
    void Finish()
    {
      const std::size_t numparam_master = DRT::UTILS::DisTypeToNumNodePerEle<M_DISTYPE>::numNodePerElement;
      const std::size_t numparam_slave = DRT::UTILS::DisTypeToNumNodePerEle<S_DISTYPE>::numNodePerElement;

      // assemble all blocks, i.e. mat and rhs
      for(std::size_t idof = 0; idof < NUMDOF; idof++)
      {
        //-----------------------------------
        // rhs
        //-----------------------------------
        // master row
        for (std::size_t vi=0; vi<numparam_master; ++vi)
        {
          int ridx = vi*NUMDOF+idof;
          int rpatch =lm_masterNodeToPatch_[vi];

          (eforce_block_[idof])(rpatch) += elevector_m_(ridx);
//          for (std::size_t i=0; i<estif_block_.size(); i++)
//              eforce_block_[i].Print(std::cout);
        }
        // slave row
        for (std::size_t vi=0; vi<numparam_slave; ++vi)
        {
          int ridx = vi*NUMDOF+idof;
          int rpatch = lm_slaveNodeToPatch_[vi];

          (eforce_block_[idof])(rpatch) += elevector_s_(ridx);
//          for (std::size_t i=0; i<estif_block_.size(); i++)
//                        eforce_block_[i].Print(std::cout);
        }

        //-----------------------------------
        // mat
        //-----------------------------------
        for(std::size_t jdof=0; jdof < NUMDOF; jdof++)
        {
          ///////////////////////////////
          // test
          bool assemble = false;
          int block = 100;
          if (pattern_ == INPAR::FLUID::EOS_GP_Pattern_uvwp)
          {
            if (idof == jdof)
            {
              assemble = true;
              block = idof;
            }
          }
          else if (pattern_ == INPAR::FLUID::EOS_GP_Pattern_up)
          {
            if (idof < (NUMDOF-1) and jdof < (NUMDOF-1))
            {
              assemble = true;
              block = idof*(NUMDOF-1)+jdof;
            }
            else
            {
              if (idof == (NUMDOF-1) and jdof == (NUMDOF-1))//(idof == jdof)
              {
                assemble = true;
                block = (NUMDOF-1)*(NUMDOF-1);
              }
            }
          }
          else
            dserror("Unknown pattern!");
          /////////////////////////////////
//          std::cout << "assembly " << assemble << std::endl;
//          std::cout << "idof " << idof << std::endl;
//          std::cout << "jdof " << jdof << std::endl;

          if (assemble){
//              std::cout << "block " << block << std::endl;
          // master row
          for (std::size_t vi=0; vi<numparam_master; ++vi)
          {
            int ridx = vi*NUMDOF+idof;
            int rpatch =lm_masterNodeToPatch_[vi];

            //master col
            for (std::size_t ui=0; ui<numparam_master; ++ui)
            {
              int cidx = ui*NUMDOF+jdof;
              int cpatch =lm_masterNodeToPatch_[ui];

              (estif_block_[block])(rpatch,cpatch) += elematrix_mm_(ridx ,cidx); //idof*NUMDOF+jdof
            }
          }
          // slave row
          for (std::size_t vi=0; vi<numparam_slave; ++vi)
          {
            int ridx = vi*NUMDOF+idof;
            int rpatch = lm_slaveNodeToPatch_[vi];

            //master col
            for (std::size_t ui=0; ui<numparam_master; ++ui)
            {
              int cidx = ui*NUMDOF+jdof;
              int cpatch = lm_masterNodeToPatch_[ui];

              (estif_block_[block])(rpatch,cpatch) += elematrix_sm_(ridx ,cidx);
            }
          }
          // master row
          for (std::size_t vi=0; vi<numparam_master; ++vi)
          {
            int ridx = vi*NUMDOF+idof;
            int rpatch = lm_masterNodeToPatch_[vi];

            // slave col
            for (std::size_t ui=0; ui<numparam_slave; ++ui)
            {
              int cidx = ui*NUMDOF+jdof;
              int cpatch = lm_slaveNodeToPatch_[ui];

              (estif_block_[block])(rpatch,cpatch) += elematrix_ms_(ridx ,cidx);
            }
          }
          // slave row
          for (std::size_t vi=0; vi<numparam_slave; ++vi)
          {
            int ridx = vi*NUMDOF+idof;
            int rpatch = lm_slaveNodeToPatch_[vi];

            // slave col
            for (std::size_t ui=0; ui<numparam_slave; ++ui)
            {
              int cidx = ui*NUMDOF+jdof;
              int cpatch = lm_slaveNodeToPatch_[ui];

              (estif_block_[block])(rpatch,cpatch) += elematrix_ss_(ridx ,cidx);
            }
          }
        }
        }
      }

      return;
    }

private:

    // don't want = operator and cctor
    LocalPatchAssembler operator = (const LocalPatchAssembler& old);
    LocalPatchAssembler(const LocalPatchAssembler& old);

    std::vector<Epetra_SerialDenseMatrix>& estif_block_;
    std::vector<Epetra_SerialDenseVector>& eforce_block_;

    const INPAR::FLUID::EOS_GP_Pattern& pattern_;

    const std::vector<int>& lm_masterNodeToPatch_; ///< local map between master nodes and nodes in patch
    const std::vector<int>& lm_slaveNodeToPatch_;  ///< local map between slave nodes and nodes in patch

    LINALG::Matrix<NUMDOF*DRT::UTILS::DisTypeToNumNodePerEle<M_DISTYPE>::numNodePerElement,NUMDOF*DRT::UTILS::DisTypeToNumNodePerEle<M_DISTYPE>::numNodePerElement> elematrix_mm_;  // element matrix master-master block
    LINALG::Matrix<NUMDOF*DRT::UTILS::DisTypeToNumNodePerEle<M_DISTYPE>::numNodePerElement,NUMDOF*DRT::UTILS::DisTypeToNumNodePerEle<S_DISTYPE>::numNodePerElement> elematrix_ms_;  // element matrix master-slave block
    LINALG::Matrix<NUMDOF*DRT::UTILS::DisTypeToNumNodePerEle<S_DISTYPE>::numNodePerElement,NUMDOF*DRT::UTILS::DisTypeToNumNodePerEle<M_DISTYPE>::numNodePerElement> elematrix_sm_;  // element matrix slave-master block
    LINALG::Matrix<NUMDOF*DRT::UTILS::DisTypeToNumNodePerEle<S_DISTYPE>::numNodePerElement,NUMDOF*DRT::UTILS::DisTypeToNumNodePerEle<S_DISTYPE>::numNodePerElement> elematrix_ss_;  // element matrix slave-slave block

    LINALG::Matrix<NUMDOF*DRT::UTILS::DisTypeToNumNodePerEle<M_DISTYPE>::numNodePerElement,1> elevector_m_; // element vector master block
    LINALG::Matrix<NUMDOF*DRT::UTILS::DisTypeToNumNodePerEle<S_DISTYPE>::numNodePerElement,1> elevector_s_; // element vector slave block

};


}
#endif
