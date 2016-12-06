/*----------------------------------------------------------------------*/
/*!
\file combust3_local_assembler.H

\brief routines to assemble element stiffness matrix and vector

\level 2

<pre>
\maintainer Benedikt Schott
            schott@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>
*/
/*----------------------------------------------------------------------*/

#ifndef COMBUST3_LOCAL_ASSEMBLER_H
#define COMBUST3_LOCAL_ASSEMBLER_H


#include "combust_defines.H"
#include "../drt_inpar/inpar_fluid.H"
#include "../drt_xfem/xfem_utils_combust.H"
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

    void PrintMat()
    {
      estif_.Print(std::cout);
      return;
    }

    void PrintRHS()
    {
      eforce_.Print(std::cout);
      return;
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

    void PrintMat()
    {
      std::cout << estif_ << std::endl;
      return;
    }

    void PrintRHS()
    {
      std::cout << eforce_ << std::endl;
      return;
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
            std::map<std::pair<XFEM::PHYSICS::Field,XFEM::PHYSICS::Field>,Epetra_SerialDenseMatrix>& estif,
            std::map<XFEM::PHYSICS::Field,Epetra_SerialDenseVector>& eforce,
            const INPAR::FLUID::EOS_GP_Pattern& pattern,
            const std::map<XFEM::PHYSICS::Field,std::vector<int> >& lm_masterDofPerFieldToPatch,
            const std::map<XFEM::PHYSICS::Field,std::vector<int> >& lm_slaveDofPerFieldToPatch
            ) :
        master_dofman_(m_dofman),
        slave_dofman_(s_dofman),
        estif_block_(estif),
        eforce_block_(eforce),
        pattern_(pattern),
        lm_masterDofPerFieldToPatch_(lm_masterDofPerFieldToPatch),
        lm_slaveDofPerFieldToPatch_(lm_slaveDofPerFieldToPatch)
    {
      const size_t m_numdofs = master_dofman_.NumNodeDof();
      const size_t s_numdofs = slave_dofman_.NumNodeDof();

      // all values are automatically initialized to zero
      elematrix_mm_ = Teuchos::rcp(new Epetra_SerialDenseMatrix(m_numdofs,m_numdofs));
      elematrix_ms_ = Teuchos::rcp(new Epetra_SerialDenseMatrix(m_numdofs,s_numdofs));
      elematrix_sm_ = Teuchos::rcp(new Epetra_SerialDenseMatrix(s_numdofs,m_numdofs));
      elematrix_ss_ = Teuchos::rcp(new Epetra_SerialDenseMatrix(s_numdofs,s_numdofs));

      elevector_m_ = Teuchos::rcp(new Epetra_SerialDenseVector(m_numdofs));
      elevector_s_ = Teuchos::rcp(new Epetra_SerialDenseVector(s_numdofs));

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

        return;
    }

    void PrintMat()
    {
      elematrix_mm_->Print(std::cout);
      elematrix_ms_->Print(std::cout);
      elematrix_sm_->Print(std::cout);
      elematrix_ss_->Print(std::cout);

      std::map<std::pair<XFEM::PHYSICS::Field,XFEM::PHYSICS::Field>,Epetra_SerialDenseMatrix>::iterator it;
      for (it=estif_block_.begin(); it!= estif_block_.end(); it++)
      {
        std::cout << XFEM::PHYSICS::physVarToString((it->first).first) << " " << XFEM::PHYSICS::physVarToString((it->first).second) << std::endl;
        it->second.Print(std::cout);
      }

      return;
    }

    void PrintRhs()
    {
      elevector_m_->Print(std::cout);
      elevector_s_->Print(std::cout);

      std::map<XFEM::PHYSICS::Field,Epetra_SerialDenseVector>::iterator it;
      for (it=eforce_block_.begin(); it!= eforce_block_.end(); it++)
      {
        std::cout << XFEM::PHYSICS::physVarToString(it->first) << std::endl;
        it->second.Print(std::cout);
      }

      return;
    }

    // assembles master, slave and coupled matrices into the patch matrix blocks
    void Finish()
    {
      // assemble all blocks

      int numfieldinpatch = master_dofman_.NumFields();
      if (numfieldinpatch != slave_dofman_.NumFields()) dserror("Same number of fields expected!");

      const std::vector<XFEM::PHYSICS::Field> fields = master_dofman_.GetFields();

      // define vector of blocks
      std::vector<std::pair<XFEM::PHYSICS::Field,XFEM::PHYSICS::Field> > block_vec;

      if (pattern_ == INPAR::FLUID::EOS_GP_Pattern_uvwp)
      {
        // 3D: 4 blocks =  u-u block, v-v block, w-w block and p-p block
        for (std::size_t ifield = 0; ifield < fields.size(); ifield++)
          block_vec.push_back(std::make_pair(fields[ifield],fields[ifield]));
      }
      else if(pattern_ == INPAR::FLUID::EOS_GP_Pattern_up)
      {
        // 3D: 10 blocks = 3x3 u-u blocks + 1x1 p-p block
        for (std::size_t ifield = 0; ifield < fields.size(); ifield++)
        {
          for (std::size_t jfield = 0; jfield < fields.size(); jfield++)
          {
            bool add = false;
            if (fields[ifield] == fields[jfield])
            {
              add = true;
            }
            else
            {
              if (fields[ifield] != XFEM::PHYSICS::Pres and fields[jfield] != XFEM::PHYSICS::Pres)
                add = true;
            }

            if (add)
             block_vec.push_back(std::make_pair(fields[ifield],fields[jfield]));
          }
        }
      }
      else if(pattern_ == INPAR::FLUID::EOS_GP_Pattern_full)
      {
        // 3D: 16 blocks = 4x4 uvwp blocks
        dserror("Full pattern does not make sense currently");
        // for the time being no terms that require all blocks are implemented
      }
      else dserror("unknown matrix pattern");

      // assemble all blocks, i.e. mat and rhs

      //-----------------------------------
      // rhs
      //-----------------------------------
      for(std::size_t ifield = 0; ifield < fields.size(); ifield++) // das sollten die Felder sein
      {
        // master row
        // get dof position for current field
        std::vector<int> m_testdof = master_dofman_.LocalDofPosPerField(fields[ifield]);

        std::map<XFEM::PHYSICS::Field,std::vector<int> >::const_iterator m_patchdofloc = lm_masterDofPerFieldToPatch_.find(fields[ifield]);
        if (m_patchdofloc == lm_masterDofPerFieldToPatch_.end()) dserror("Field not found");
        std::vector<int> m_patchdof = m_patchdofloc->second;
        if (m_patchdof.size() != m_testdof.size()) dserror("Mismatch in number of dofs!");

        int vi = 0;
        for (std::vector<int>::const_iterator testpos = m_testdof.begin(); testpos != m_testdof.end(); ++testpos)
        {
          int rpatch = m_patchdof[vi]; //(lm_masterDofPerFieldToPatch_[XFEM::PHYSICS::Velx])[vi];
          (eforce_block_[fields[ifield]])(rpatch) += (*elevector_m_)(*testpos);
          vi++;
        }

        // slave row
        // get dof position for current field
//        std::vector<int> s_testdof = slave_dofman_.LocalDofPosPerField<SymmetryMangle<fields[ifield]>::name>();
        std::vector<int> s_testdof = slave_dofman_.LocalDofPosPerField(fields[ifield]);

        std::map<XFEM::PHYSICS::Field,std::vector<int> >::const_iterator s_patchdofloc = lm_slaveDofPerFieldToPatch_.find(fields[ifield]);
        if (s_patchdofloc == lm_slaveDofPerFieldToPatch_.end()) dserror("Field not found");
        std::vector<int> s_patchdof = s_patchdofloc->second;
        if (s_patchdof.size() != s_testdof.size()) dserror("Mismatch in number of dofs!");

        vi = 0;
        for (std::vector<int>::const_iterator testpos = s_testdof.begin(); testpos != s_testdof.end(); ++testpos)
        {
          int rpatch = s_patchdof[vi]; //(lm_slaveDofPerFieldToPatch_[fields[ifield]])[vi];
          (eforce_block_[fields[ifield]])(rpatch) += (*elevector_s_)(*testpos);
          vi++;
        }
      }

      //-----------------------------------
      // mat
      //-----------------------------------
      for(std::size_t ib = 0; ib < block_vec.size(); ib++)
      {
        std::vector<int> m_testdof = master_dofman_.LocalDofPosPerField(block_vec[ib].first);
        std::vector<int> m_trialdof = master_dofman_.LocalDofPosPerField(block_vec[ib].second);

        std::map<XFEM::PHYSICS::Field,std::vector<int> >::const_iterator m_patchdofloc_test = lm_masterDofPerFieldToPatch_.find(block_vec[ib].first);
        if (m_patchdofloc_test == lm_masterDofPerFieldToPatch_.end()) dserror("Field not found");
        std::vector<int> m_patchdof_test = m_patchdofloc_test->second;
        if (m_patchdof_test.size() != m_testdof.size()) dserror("Mismatch in number of dofs!");

        std::map<XFEM::PHYSICS::Field,std::vector<int> >::const_iterator m_patchdofloc_trial = lm_masterDofPerFieldToPatch_.find(block_vec[ib].second);
        if (m_patchdofloc_trial == lm_masterDofPerFieldToPatch_.end()) dserror("Field not found");
        std::vector<int> m_patchdof_trial = m_patchdofloc_trial->second;
        if (m_patchdof_trial.size() != m_trialdof.size()) dserror("Mismatch in number of dofs!");


        std::vector<int> s_testdof = slave_dofman_.LocalDofPosPerField(block_vec[ib].first);
        std::vector<int> s_trialdof = slave_dofman_.LocalDofPosPerField(block_vec[ib].second);

        std::map<XFEM::PHYSICS::Field,std::vector<int> >::const_iterator s_patchdofloc_test = lm_slaveDofPerFieldToPatch_.find(block_vec[ib].first);
        if (s_patchdofloc_test == lm_slaveDofPerFieldToPatch_.end()) dserror("Field not found");
        std::vector<int> s_patchdof_test = s_patchdofloc_test->second;
        if (s_patchdof_test.size() != s_testdof.size()) dserror("Mismatch in number of dofs!");

        std::map<XFEM::PHYSICS::Field,std::vector<int> >::const_iterator s_patchdofloc_trial = lm_slaveDofPerFieldToPatch_.find(block_vec[ib].second);
        if (s_patchdofloc_trial == lm_slaveDofPerFieldToPatch_.end()) dserror("Field not found");
        std::vector<int> s_patchdof_trial = s_patchdofloc_trial->second;
        if (s_patchdof_trial.size() != s_trialdof.size()) dserror("Mismatch in number of dofs!");

        int vi = 0;
        for (std::vector<int>::const_iterator testpos = m_testdof.begin(); testpos != m_testdof.end(); ++testpos)
        {
          int rpatch = m_patchdof_test[vi];

          int ui = 0;
          for (std::vector<int>::const_iterator trialpos = m_trialdof.begin(); trialpos != m_trialdof.end(); ++trialpos)
          {
            int cpatch = m_patchdof_trial[ui];
            (estif_block_[block_vec[ib]])(rpatch,cpatch) += (*elematrix_mm_)(*testpos ,*trialpos);
            ui++;
          }

          ui = 0;
          for (std::vector<int>::const_iterator trialpos = s_trialdof.begin(); trialpos != s_trialdof.end(); ++trialpos)
          {
            int cpatch = s_patchdof_trial[ui];
            (estif_block_[block_vec[ib]])(rpatch,cpatch) += (*elematrix_ms_)(*testpos ,*trialpos);
            ui++;
          }

          vi++;
        }

        vi = 0;
        for (std::vector<int>::const_iterator testpos = s_testdof.begin(); testpos != s_testdof.end(); ++testpos)
        {
          int rpatch = s_patchdof_test[vi];

          int ui = 0;
          for (std::vector<int>::const_iterator trialpos = s_trialdof.begin(); trialpos != s_trialdof.end(); ++trialpos)
          {
            int cpatch = s_patchdof_trial[ui];
            (estif_block_[block_vec[ib]])(rpatch,cpatch) += (*elematrix_ss_)(*testpos ,*trialpos);
            ui++;
          }

          ui = 0;
          for (std::vector<int>::const_iterator trialpos = m_trialdof.begin(); trialpos != m_trialdof.end(); ++trialpos)
          {
            int cpatch = m_patchdof_trial[ui];
            (estif_block_[block_vec[ib]])(rpatch,cpatch) += (*elematrix_sm_)(*testpos ,*trialpos);
            ui++;
          }

          vi++;
        }

      }

      return;
    }

private:

    // don't want = operator and cctor
    LocalPatchAssembler operator = (const LocalPatchAssembler& old);
    LocalPatchAssembler(const LocalPatchAssembler& old);

    const XFEM::ElementDofManager& master_dofman_;
    const XFEM::ElementDofManager& slave_dofman_;
    std::map<std::pair<XFEM::PHYSICS::Field,XFEM::PHYSICS::Field>,Epetra_SerialDenseMatrix>& estif_block_;
    std::map<XFEM::PHYSICS::Field,Epetra_SerialDenseVector>& eforce_block_;

    const INPAR::FLUID::EOS_GP_Pattern& pattern_;

    const std::map<XFEM::PHYSICS::Field,std::vector<int> >& lm_masterDofPerFieldToPatch_; ///< local map between master nodes and nodes in patch
    const std::map<XFEM::PHYSICS::Field,std::vector<int> >& lm_slaveDofPerFieldToPatch_;  ///< local map between slave nodes and nodes in patch

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
            const XFEM::ElementDofManager&   m_dofman,
            const XFEM::ElementDofManager&   s_dofman,
            std::map<std::pair<XFEM::PHYSICS::Field,XFEM::PHYSICS::Field>,Epetra_SerialDenseMatrix>& estif,
            std::map<XFEM::PHYSICS::Field,Epetra_SerialDenseVector>& eforce,
            const INPAR::FLUID::EOS_GP_Pattern& pattern,
            const std::map<XFEM::PHYSICS::Field,std::vector<int> >& lm_masterDofPerFieldToPatch,
            const std::map<XFEM::PHYSICS::Field,std::vector<int> >& lm_slaveDofPerFieldToPatch
            ) :
        estif_block_(estif),
        eforce_block_(eforce),
        pattern_(pattern),
        lm_masterDofPerFieldToPatch_(lm_masterDofPerFieldToPatch),
        lm_slaveDofPerFieldToPatch_(lm_slaveDofPerFieldToPatch),
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
            std::map<std::pair<XFEM::PHYSICS::Field,XFEM::PHYSICS::Field>,Epetra_SerialDenseMatrix>& estif,
            std::map<XFEM::PHYSICS::Field,Epetra_SerialDenseVector>& eforce,
            const INPAR::FLUID::EOS_GP_Pattern& pattern,
            const std::map<XFEM::PHYSICS::Field,std::vector<int> >& lm_masterDofPerFieldToPatch,
            const std::map<XFEM::PHYSICS::Field,std::vector<int> >& lm_slaveDofPerFieldToPatch
            ) :
        estif_block_(estif),
        eforce_block_(eforce),
        pattern_(pattern),
        lm_masterDofPerFieldToPatch_(lm_masterDofPerFieldToPatch),
        lm_slaveDofPerFieldToPatch_(lm_slaveDofPerFieldToPatch),
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

      std::map<std::pair<XFEM::PHYSICS::Field,XFEM::PHYSICS::Field>,Epetra_SerialDenseMatrix>::iterator it;
      for (it=estif_block_.begin(); it!= estif_block_.end(); it++)
      {
        std::cout << XFEM::PHYSICS::physVarToString((it->first).first) << " " << XFEM::PHYSICS::physVarToString((it->first).second) << std::endl;
        it->second.Print(std::cout);
      }
      return;
    }

    void PrintRhs()
    {
      std::cout << "elevector_m_ " << elevector_m_;
      std::cout << "elevector_s_ " << elevector_s_;

      std::map<XFEM::PHYSICS::Field,Epetra_SerialDenseVector>::iterator it;
      for (it=eforce_block_.begin(); it!= eforce_block_.end(); it++)
      {
        std::cout << XFEM::PHYSICS::physVarToString(it->first) << std::endl;
         it->second.Print(std::cout);
      }
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
        XFEM::PHYSICS::Field currentfield = XFEM::PHYSICS::Velx;
        if (idof == 1) currentfield = XFEM::PHYSICS::Vely;
        if (idof == 2) currentfield = XFEM::PHYSICS::Velz;
        if (idof == 3) currentfield = XFEM::PHYSICS::Pres;

        std::map<XFEM::PHYSICS::Field,std::vector<int> >::const_iterator m_patchdofloc = lm_masterDofPerFieldToPatch_.find(currentfield);
        if (m_patchdofloc == lm_masterDofPerFieldToPatch_.end()) dserror("Field not found");
        std::vector<int> m_patchdof = m_patchdofloc->second;
//        if (m_patchdof.size() != m_testdof.size()) dserror("Mismatch in number of dofs!");

        std::map<XFEM::PHYSICS::Field,std::vector<int> >::const_iterator s_patchdofloc = lm_slaveDofPerFieldToPatch_.find(currentfield);
        if (s_patchdofloc == lm_slaveDofPerFieldToPatch_.end()) dserror("Field not found");
        std::vector<int> s_patchdof = s_patchdofloc->second;
//        if (s_patchdof.size() != s_testdof.size()) dserror("Mismatch in number of dofs!");
        //-----------------------------------
        // rhs
        //-----------------------------------
        // master row
        for (std::size_t vi=0; vi<numparam_master; ++vi)
        {
          int ridx = vi*NUMDOF+idof;
          int rpatch = m_patchdof[vi];//lm_masterDofPerFieldToPatch_[vi];

          (eforce_block_[currentfield])(rpatch) += elevector_m_(ridx);

        }
        // slave row
        for (std::size_t vi=0; vi<numparam_slave; ++vi)
        {
          int ridx = vi*NUMDOF+idof;
          int rpatch = s_patchdof[vi];//lm_slaveDofPerFieldToPatch_[vi];

          (eforce_block_[currentfield])(rpatch) += elevector_s_(ridx);

        }

        //-----------------------------------
        // mat
        //-----------------------------------
        for(std::size_t jdof=0; jdof < NUMDOF; jdof++)
        {
          XFEM::PHYSICS::Field scurrentfield = XFEM::PHYSICS::Velx;
          if (jdof == 1) scurrentfield = XFEM::PHYSICS::Vely;
          if (jdof == 2) scurrentfield = XFEM::PHYSICS::Velz;
          if (jdof == 3) scurrentfield = XFEM::PHYSICS::Pres;

          bool assemble = false;

          if (pattern_ == INPAR::FLUID::EOS_GP_Pattern_uvwp)
          {
            if (idof == jdof)
            {
              assemble = true;
            }
          }
          else if (pattern_ == INPAR::FLUID::EOS_GP_Pattern_up)
          {
            if (idof < (NUMDOF-1) and jdof < (NUMDOF-1))
            {
              assemble = true;
            }
            else
            {
              if (idof == (NUMDOF-1) and jdof == (NUMDOF-1))//(idof == jdof)
              {
                assemble = true;
              }
            }
          }
          else
            dserror("Unknown pattern!");


          if (assemble){

          // master row
          for (std::size_t vi=0; vi<numparam_master; ++vi)
          {
            int ridx = vi*NUMDOF+idof;
            int rpatch = m_patchdof[vi];//lm_masterDofPerFieldToPatch_[vi];

            //master col
            for (std::size_t ui=0; ui<numparam_master; ++ui)
            {
              int cidx = ui*NUMDOF+jdof;
              int cpatch = m_patchdof[ui];//lm_masterDofPerFieldToPatch_[ui];

              (estif_block_[std::make_pair(currentfield,scurrentfield)])(rpatch,cpatch) += elematrix_mm_(ridx ,cidx); //idof*NUMDOF+jdof
            }
          }
          // slave row
          for (std::size_t vi=0; vi<numparam_slave; ++vi)
          {
            int ridx = vi*NUMDOF+idof;
            int rpatch = s_patchdof[vi]; //lm_slaveDofPerFieldToPatch_[vi];

            //master col
            for (std::size_t ui=0; ui<numparam_master; ++ui)
            {
              int cidx = ui*NUMDOF+jdof;
              int cpatch = m_patchdof[ui];//lm_masterDofPerFieldToPatch_[ui];

              (estif_block_[std::make_pair(currentfield,scurrentfield)])(rpatch,cpatch) += elematrix_sm_(ridx ,cidx);
            }
          }
          // master row
          for (std::size_t vi=0; vi<numparam_master; ++vi)
          {
            int ridx = vi*NUMDOF+idof;
            int rpatch = m_patchdof[vi];//lm_masterDofPerFieldToPatch_[vi];

            // slave col
            for (std::size_t ui=0; ui<numparam_slave; ++ui)
            {
              int cidx = ui*NUMDOF+jdof;
              int cpatch = s_patchdof[ui];//lm_slaveDofPerFieldToPatch_[ui];

              (estif_block_[std::make_pair(currentfield,scurrentfield)])(rpatch,cpatch) += elematrix_ms_(ridx ,cidx);
            }
          }
          // slave row
          for (std::size_t vi=0; vi<numparam_slave; ++vi)
          {
            int ridx = vi*NUMDOF+idof;
            int rpatch = s_patchdof[vi]; //lm_slaveDofPerFieldToPatch_[vi];

            // slave col
            for (std::size_t ui=0; ui<numparam_slave; ++ui)
            {
              int cidx = ui*NUMDOF+jdof;
              int cpatch = s_patchdof[ui];//lm_slaveDofPerFieldToPatch_[ui];

              (estif_block_[std::make_pair(currentfield,scurrentfield)])(rpatch,cpatch) += elematrix_ss_(ridx ,cidx);
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

    std::map<std::pair<XFEM::PHYSICS::Field,XFEM::PHYSICS::Field>,Epetra_SerialDenseMatrix>& estif_block_;
    std::map<XFEM::PHYSICS::Field,Epetra_SerialDenseVector>& eforce_block_;

    const INPAR::FLUID::EOS_GP_Pattern& pattern_;

    const std::map<XFEM::PHYSICS::Field,std::vector<int> >& lm_masterDofPerFieldToPatch_; ///< local map between master nodes and nodes in patch
    const std::map<XFEM::PHYSICS::Field,std::vector<int> >& lm_slaveDofPerFieldToPatch_;  ///< local map between slave nodes and nodes in patch

    LINALG::Matrix<NUMDOF*DRT::UTILS::DisTypeToNumNodePerEle<M_DISTYPE>::numNodePerElement,NUMDOF*DRT::UTILS::DisTypeToNumNodePerEle<M_DISTYPE>::numNodePerElement> elematrix_mm_;  // element matrix master-master block
    LINALG::Matrix<NUMDOF*DRT::UTILS::DisTypeToNumNodePerEle<M_DISTYPE>::numNodePerElement,NUMDOF*DRT::UTILS::DisTypeToNumNodePerEle<S_DISTYPE>::numNodePerElement> elematrix_ms_;  // element matrix master-slave block
    LINALG::Matrix<NUMDOF*DRT::UTILS::DisTypeToNumNodePerEle<S_DISTYPE>::numNodePerElement,NUMDOF*DRT::UTILS::DisTypeToNumNodePerEle<M_DISTYPE>::numNodePerElement> elematrix_sm_;  // element matrix slave-master block
    LINALG::Matrix<NUMDOF*DRT::UTILS::DisTypeToNumNodePerEle<S_DISTYPE>::numNodePerElement,NUMDOF*DRT::UTILS::DisTypeToNumNodePerEle<S_DISTYPE>::numNodePerElement> elematrix_ss_;  // element matrix slave-slave block

    LINALG::Matrix<NUMDOF*DRT::UTILS::DisTypeToNumNodePerEle<M_DISTYPE>::numNodePerElement,1> elevector_m_; // element vector master block
    LINALG::Matrix<NUMDOF*DRT::UTILS::DisTypeToNumNodePerEle<S_DISTYPE>::numNodePerElement,1> elevector_s_; // element vector slave block

};

}

#endif
