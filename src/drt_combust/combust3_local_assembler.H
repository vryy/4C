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

#ifdef CCADISCRET

#ifndef COMBUST3_LOCAL_ASSEMBLER_H
#define COMBUST3_LOCAL_ASSEMBLER_H

#include "../drt_fem_general/drt_utils_local_connectivity_matrices.H"
#include "../drt_xfem/physics.H"
#include "../drt_xfem/dof_management_element.H"
#include <Epetra_SerialDenseMatrix.h>
#include <Epetra_SerialDenseVector.h>

enum DerivDir
{
  D0 = 0,
  Dx = 1,
  Dy = 2,
  Dz = 3,
  DxDx = 4,
  DyDy = 5,
  DzDz = 6,
  DxDy = 7,
  DxDz = 8,
  DyDz = 9,
  DyDx = 7,
  DzDx = 8,
  DzDy = 9
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
struct SymmetryMangle {static const XFEM::PHYSICS::Field name = FIELD;};
template<> struct SymmetryMangle<XFEM::PHYSICS::Sigmayx> {static const XFEM::PHYSICS::Field name = XFEM::PHYSICS::Sigmaxy;};
template<> struct SymmetryMangle<XFEM::PHYSICS::Sigmazx> {static const XFEM::PHYSICS::Field name = XFEM::PHYSICS::Sigmaxz;};
template<> struct SymmetryMangle<XFEM::PHYSICS::Sigmazy> {static const XFEM::PHYSICS::Field name = XFEM::PHYSICS::Sigmayz;};
template<> struct SymmetryMangle<XFEM::PHYSICS::Tauyx>   {static const XFEM::PHYSICS::Field name = XFEM::PHYSICS::Tauxy;};
template<> struct SymmetryMangle<XFEM::PHYSICS::Tauzx>   {static const XFEM::PHYSICS::Field name = XFEM::PHYSICS::Tauxz;};
template<> struct SymmetryMangle<XFEM::PHYSICS::Tauzy>   {static const XFEM::PHYSICS::Field name = XFEM::PHYSICS::Tauyz;};

//! which shape functions are equal anyway
template <XFEM::PHYSICS::Field FIELD>
struct ShapeMangle;
template<> struct ShapeMangle<XFEM::PHYSICS::Velx> {static const XFEM::PHYSICS::Field name = XFEM::PHYSICS::Velx;};
template<> struct ShapeMangle<XFEM::PHYSICS::Vely> {static const XFEM::PHYSICS::Field name = XFEM::PHYSICS::Velx;};
template<> struct ShapeMangle<XFEM::PHYSICS::Velz> {static const XFEM::PHYSICS::Field name = XFEM::PHYSICS::Velx;};
template<> struct ShapeMangle<XFEM::PHYSICS::Pres> {static const XFEM::PHYSICS::Field name = XFEM::PHYSICS::Velx;};
template<> struct ShapeMangle<XFEM::PHYSICS::DiscPres> {static const XFEM::PHYSICS::Field name = XFEM::PHYSICS::DiscPres;};
template<> struct ShapeMangle<XFEM::PHYSICS::Sigmaxx> {static const XFEM::PHYSICS::Field name = XFEM::PHYSICS::Sigmaxx;};
template<> struct ShapeMangle<XFEM::PHYSICS::Sigmaxy> {static const XFEM::PHYSICS::Field name = XFEM::PHYSICS::Sigmaxx;};
template<> struct ShapeMangle<XFEM::PHYSICS::Sigmaxz> {static const XFEM::PHYSICS::Field name = XFEM::PHYSICS::Sigmaxx;};
template<> struct ShapeMangle<XFEM::PHYSICS::Sigmayx> {static const XFEM::PHYSICS::Field name = XFEM::PHYSICS::Sigmaxx;};
template<> struct ShapeMangle<XFEM::PHYSICS::Sigmayy> {static const XFEM::PHYSICS::Field name = XFEM::PHYSICS::Sigmaxx;};
template<> struct ShapeMangle<XFEM::PHYSICS::Sigmayz> {static const XFEM::PHYSICS::Field name = XFEM::PHYSICS::Sigmaxx;};
template<> struct ShapeMangle<XFEM::PHYSICS::Sigmazx> {static const XFEM::PHYSICS::Field name = XFEM::PHYSICS::Sigmaxx;};
template<> struct ShapeMangle<XFEM::PHYSICS::Sigmazy> {static const XFEM::PHYSICS::Field name = XFEM::PHYSICS::Sigmaxx;};
template<> struct ShapeMangle<XFEM::PHYSICS::Sigmazz> {static const XFEM::PHYSICS::Field name = XFEM::PHYSICS::Sigmaxx;};
template<> struct ShapeMangle<XFEM::PHYSICS::Tauxx>   {static const XFEM::PHYSICS::Field name = XFEM::PHYSICS::Tauxx;};
template<> struct ShapeMangle<XFEM::PHYSICS::Tauxy>   {static const XFEM::PHYSICS::Field name = XFEM::PHYSICS::Tauxx;};
template<> struct ShapeMangle<XFEM::PHYSICS::Tauxz>   {static const XFEM::PHYSICS::Field name = XFEM::PHYSICS::Tauxx;};
template<> struct ShapeMangle<XFEM::PHYSICS::Tauyx>   {static const XFEM::PHYSICS::Field name = XFEM::PHYSICS::Tauxx;};
template<> struct ShapeMangle<XFEM::PHYSICS::Tauyy>   {static const XFEM::PHYSICS::Field name = XFEM::PHYSICS::Tauxx;};
template<> struct ShapeMangle<XFEM::PHYSICS::Tauyz>   {static const XFEM::PHYSICS::Field name = XFEM::PHYSICS::Tauxx;};
template<> struct ShapeMangle<XFEM::PHYSICS::Tauzx>   {static const XFEM::PHYSICS::Field name = XFEM::PHYSICS::Tauxx;};
template<> struct ShapeMangle<XFEM::PHYSICS::Tauzy>   {static const XFEM::PHYSICS::Field name = XFEM::PHYSICS::Tauxx;};
template<> struct ShapeMangle<XFEM::PHYSICS::Tauzz>   {static const XFEM::PHYSICS::Field name = XFEM::PHYSICS::Tauxx;};

template <XFEM::PHYSICS::Field FIELD>
struct ShapeSelector {static BlitzMat shp;};
//template<> struct ShapeSelector<XFEM::PHYSICS::Velx> {static const XFEM::PHYSICS::Field name = XFEM::PHYSICS::Sigmaxy;};
//template<> struct ShapeSelector<XFEM::PHYSICS::Vely> {static const XFEM::PHYSICS::Field name = XFEM::PHYSICS::Sigmaxz;};
//template<> struct ShapeSelector<XFEM::PHYSICS::Velz> {static const XFEM::PHYSICS::Field name = XFEM::PHYSICS::Sigmayz;};
//template<> struct ShapeSelector<XFEM::PHYSICS::Pres> {static const XFEM::PHYSICS::Field name = XFEM::PHYSICS::Tauxy;};
//template<> struct ShapeSelector<XFEM::PHYSICS::Tauzx>   {static const XFEM::PHYSICS::Field name = XFEM::PHYSICS::Tauxz;};
//template<> struct ShapeSelector<XFEM::PHYSICS::Tauzy>   {static const XFEM::PHYSICS::Field name = XFEM::PHYSICS::Tauyz;};

//! class that assembles one after another integral of the weak form into the stiffness matrix
template <DRT::Element::DiscretizationType DISTYPE, 
          XFEM::AssemblyType ASSEMBLYTYPE>
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
        const std::vector<int>& testdof  = dofman_.LocalDofPosPerField(SymmetryMangle<testfield>::name);
        const std::vector<int>& trialdof = dofman_.LocalDofPosPerField(SymmetryMangle<trialfield>::name);
        
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
    
    //! assemble one expression for the element stiffness matrix for a given test and trial field
    template <XFEM::PHYSICS::Field testfield,
              DerivDir testDerivId,
              XFEM::PHYSICS::Field trialfield,
              DerivDir trialDerivId>
    void Matrix(
            const double&       fac
            ) const
    {
        const std::vector<int>& testdof  = dofman_.LocalDofPosPerField<SymmetryMangle<testfield>::name>();
        const std::vector<int>& trialdof = dofman_.LocalDofPosPerField<SymmetryMangle<trialfield>::name>();
        
        int ui = 0;
        for (std::vector<int>::const_iterator trialpos = trialdof.begin(); trialpos != trialdof.end(); ++trialpos)
        {
            int vi = 0;
            for (std::vector<int>::const_iterator testpos = testdof.begin(); testpos != testdof.end(); ++testpos)
            {
                estif_(*testpos, *trialpos) += fac*ShapeSelector<ShapeMangle<SymmetryMangle<testfield>::name>::name >::shp(testDerivId,vi)*ShapeSelector<trialfield>::shp(trialDerivId,ui) ;
                vi++;
            }
            ui++;
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
        const std::vector<int>& testdof = dofman_.LocalDofPosPerField<SymmetryMangle<testfield>::name>();
        
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
    Epetra_SerialDenseMatrix& estif_;
    Epetra_SerialDenseVector& eforce_;
};

/*!
 * class that assembles one after another integral of the weak form into the stiffness matrix
 * 
 * specialization for standard assembly
 */
template <DRT::Element::DiscretizationType DISTYPE>
class LocalAssembler<DISTYPE, XFEM::standard_assembly>
{
public:
    
    explicit LocalAssembler(
            const XFEM::ElementDofManager&    ,
            Epetra_SerialDenseMatrix&                         estif,
            Epetra_SerialDenseVector&                         eforce
            ) :
        estif_(estif),
        eforce_(eforce)
    {
        return;
    }
            
    explicit LocalAssembler(
            Epetra_SerialDenseMatrix&                         estif,
            Epetra_SerialDenseVector&                         eforce
            ) :
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
              XFEM::PHYSICS::Field trialfield,
              class V1,
              class V2>
    void Matrix(
            const V1&           testshape,
            const double&       fac,
            const V2&           trialshape
            ) const
    {
        const int numparam = DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;
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
    template <XFEM::PHYSICS::Field testfield,
              class V1>
    void Vector(
            const V1&           testshape,
            const double&       fac
            ) const
    {
        const int numparam = DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;
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
    
    Epetra_SerialDenseMatrix& estif_;
    Epetra_SerialDenseVector& eforce_;
};

  
#endif

#endif
