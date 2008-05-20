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

#include <blitz/array.h>
#include "../drt_xfem/xfem.H"
#include "../drt_lib/drt_utils.H"


enum PhysDerivDir
{
  D0 = 0,
  Dx = 0,
  Dy = 1,
  Dz = 2,
  DxDx = 0,
  DyDy = 1,
  DzDz = 2,
  DxDy = 3,
  DxDz = 4,
  DyDz = 5
};
//! indicates direction
enum LocalDerivDir
{
  Dr = 0,
  Ds = 1,
  Dt = 2,
  DrDr = 0,
  DsDs = 1,
  DtDt = 2,
  DrDs = 3,
  DrDt = 4,
  DsDt = 5
};


//! for non-xfem elements, figure out the order of unknowns per node in the array
template<XFEM::PHYSICS::Field field>
struct FieldNum {static const int pos = -1;};
template<> struct FieldNum<XFEM::PHYSICS::Velx> {static const int pos = 0;};
template<> struct FieldNum<XFEM::PHYSICS::Vely> {static const int pos = 1;};
template<> struct FieldNum<XFEM::PHYSICS::Velz> {static const int pos = 2;};
template<> struct FieldNum<XFEM::PHYSICS::Pres> {static const int pos = 3;};

//! for the viscous stress, replace shear stress components with symmetric part
template <XFEM::PHYSICS::Field FIELD>
struct MangleField {static const XFEM::PHYSICS::Field name = FIELD;};
template<> struct MangleField<XFEM::PHYSICS::Sigmayx> {static const XFEM::PHYSICS::Field name = XFEM::PHYSICS::Sigmaxy;};
template<> struct MangleField<XFEM::PHYSICS::Sigmazx> {static const XFEM::PHYSICS::Field name = XFEM::PHYSICS::Sigmaxz;};
template<> struct MangleField<XFEM::PHYSICS::Sigmazy> {static const XFEM::PHYSICS::Field name = XFEM::PHYSICS::Sigmayz;};
template<> struct MangleField<XFEM::PHYSICS::Tauyx>   {static const XFEM::PHYSICS::Field name = XFEM::PHYSICS::Tauxy;};
template<> struct MangleField<XFEM::PHYSICS::Tauzx>   {static const XFEM::PHYSICS::Field name = XFEM::PHYSICS::Tauxz;};
template<> struct MangleField<XFEM::PHYSICS::Tauzy>   {static const XFEM::PHYSICS::Field name = XFEM::PHYSICS::Tauyz;};

//! class that assembles one after another integral of the weak form into the stiffness matrix
template <DRT::Element::DiscretizationType DISTYPE, 
          XFEM::AssemblyType ASSEMBLYTYPE>
class LocalAssembler
{
public:
    
    explicit LocalAssembler(
            const XFEM::ElementDofManager&    dofman,
            BlitzMat& estif,
            BlitzVec& eforce
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
    
    
    //! assemble one expression for the element stiffness matrix for a given test and trial field
    template <XFEM::PHYSICS::Field testfield,
              XFEM::PHYSICS::Field trialfield>
    void Matrix(
            const BlitzVec&     testshape,
            const double&       fac,
            const BlitzVec&     trialshape
            ) const
    {
        const std::vector<int>& testdof  = dofman_.LocalDofPosPerField(MangleField<testfield>::name);
        const std::vector<int>& trialdof = dofman_.LocalDofPosPerField(MangleField<trialfield>::name);
        
//        const int numparamtest  = testdof.size();
//        const int numparamtrial = trialdof.size();
//        
//        for (int ui=0; ui<numparamtrial; ++ui)
//        {
//            const int trialpos = trialdof[ui];
//
//            for (int vi=0; vi<numparamtest; ++vi)
//            {
//                const int testpos = testdof[vi];
//                estif_(testpos, trialpos) += fac*testshape(vi)*trialshape(ui) ;
//            }
//        }
        int ui = 0;
        for (std::vector<int>::const_iterator trialpos = trialdof.begin(); trialpos != trialdof.end(); ++trialpos)
        {
            int vi = 0;
            for (std::vector<int>::const_iterator testpos = testdof.begin(); testpos != testdof.end(); ++testpos)
            {
                estif_(*testpos, *trialpos) += fac*testshape(vi)*trialshape(ui) ;
                vi++;
            }
            ui++;
        }
    }
    
    //! assemble one expression for the element rhs for a given test field
    template <XFEM::PHYSICS::Field testfield>
    void Vector(
            const BlitzVec&     testshape,
            const double&       fac
            ) const
    {
        const std::vector<int>& testdof = dofman_.LocalDofPosPerField(MangleField<testfield>::name);
        
//        const int numparamtest  = testdof.size();
//        
//        for (int vi=0; vi<numparamtest; ++vi)
//        {
//            const int testpos = testdof[vi];
//            eforce_(testpos) += fac*testshape(vi);
//        }
        int vi = 0;
        for (std::vector<int>::const_iterator testpos = testdof.begin(); testpos != testdof.end(); ++testpos)
        {
            eforce_(*testpos) += fac*testshape(vi);
            vi++;
        }
    }
    
private:
    // don't want = operator and cctor
    LocalAssembler operator = (const LocalAssembler& old);
    LocalAssembler(const LocalAssembler& old);
    
    const XFEM::ElementDofManager& dofman_;
    BlitzMat& estif_;
    BlitzVec& eforce_;
};

template <DRT::Element::DiscretizationType DISTYPE>
class LocalAssembler<DISTYPE, XFEM::standard_assembly>
{
public:
    
    explicit LocalAssembler(
            const XFEM::ElementDofManager&    dofman,
            BlitzMat&                         estif,
            BlitzVec&                         eforce
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
    
    //! assemble one expression for the element stiffness matrix for a given test and trial field
    template <XFEM::PHYSICS::Field testfield,
              XFEM::PHYSICS::Field trialfield>
    void Matrix(
            const BlitzVec&     testshape,
            const double&       fac,
            const BlitzVec&     trialshape
            ) const
    {
        const int numparam = DRT::UTILS::_switchDisType<DISTYPE>::numNodePerElement;
        for (int ui=0; ui<numparam; ++ui)
        {
            const int trialpos = 4*ui + FieldNum<trialfield>::pos;
            const double multiply = fac*trialshape(ui);
            for (int vi=0; vi<numparam; ++vi)
            {
                const int testpos = 4*vi + FieldNum<testfield>::pos;
                estif_(testpos, trialpos) += testshape(vi)*multiply ;
            }
        }
    }
    
    //! assemble one expression for the element rhs for a given test field
    template <XFEM::PHYSICS::Field testfield>
    void Vector(
            const BlitzVec&     testshape,
            const double&       fac
            ) const
    {
        const int numparam = DRT::UTILS::_switchDisType<DISTYPE>::numNodePerElement;
        for (int vi=0; vi<numparam; ++vi)
        {
            const int testpos = 4*vi + FieldNum<testfield>::pos;
            eforce_(testpos) += fac*testshape(vi);
        }
    }
    
private:
    // don't want = operator and cctor
    LocalAssembler operator = (const LocalAssembler& old);
    LocalAssembler(const LocalAssembler& old);
    
    const XFEM::ElementDofManager& dofman_;
    BlitzMat& estif_;
    BlitzVec& eforce_;
};

  
#endif

#endif
