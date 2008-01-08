/*----------------------------------------------------------------------*/
/*!
\file xfluid3_stationary.H

\brief Internal implementation of Fluid3 element -- stationary formulation

<pre>
Maintainer: Axel Gerstenberger
            gerstenberger@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>
*/
/*----------------------------------------------------------------------*/

#ifdef D_FLUID3
#ifdef CCADISCRET

#ifndef XFLUID3_INTEGRATORY_H
#define XFLUID3_INTEGRATORY_H

#include "xfluid3.H"
#include <blitz/array.h>
#include "../drt_mat/newtonianfluid.H"
#include "../drt_lib/drt_timecurve.H"
#include "../drt_xfem/xfem.H"
#include "../drt_xfem/intersection_service.H"
#include "../drt_lib/drt_utils.H"
#include "fluid3_stabilization.H"

#include <Epetra_SerialDenseSolver.h>

using namespace std;
using namespace XFEM::PHYSICS;


typedef blitz::Array<double,1> BlitzVec; ///< a short-cut for a one-dimensional BLITZ vector
typedef blitz::Array<double,2> BlitzMat; ///< a short-cut for a two-dimensional BLITZ matrix


//! for non-xfem elements, figure out the position in the array
template <XFEM::PHYSICS::Field field> int getFieldNumber()
{
    dserror("this field does not belong to the 3d fluid");
    return -1;
}; 
template <> int getFieldNumber<XFEM::PHYSICS::Velx>() {return 0;};
template <> int getFieldNumber<XFEM::PHYSICS::Vely>() {return 1;};
template <> int getFieldNumber<XFEM::PHYSICS::Velz>() {return 2;};
template <> int getFieldNumber<XFEM::PHYSICS::Pres>() {return 3;};




template <DRT::Element::DiscretizationType DISTYPE, 
          XFEM::AssemblyType ASSEMBLYTYPE>
class LocalAssembler
{
public:
    //! assemble one expression for the element stiffness matrix for a given test and trial field
    template <XFEM::PHYSICS::Field testfield,
              XFEM::PHYSICS::Field trialfield>
    void Matrix(
            const XFEM::ElementDofManager&    dofman,
            blitz::Array<double,2>&           estif,
            const blitz::Array<double,1>&     testshape,
            const double                      fac,
            const blitz::Array<double,1>&     trialshape
            ) const
    {
        const std::vector<int> testdof  = dofman.LocalDofPosPerField(testfield);
        const std::vector<int> trialdof = dofman.LocalDofPosPerField(trialfield);
        
        const int numparamtest  = testdof.size();
        const int numparamtrial = trialdof.size();

        for (int ui=0; ui<numparamtrial; ++ui)
        {
            const int trialpos = trialdof[ui];

            for (int vi=0; vi<numparamtest; ++vi)
            {
                const int testpos = testdof[vi];

                estif(testpos, trialpos) += fac*testshape(vi)*trialshape(ui) ;
            }
        }
    };
    
    //! assemble one expression for the element rhs for a given test field
    template <XFEM::PHYSICS::Field testfield>
    void Vector(
            const XFEM::ElementDofManager&    dofman,
            blitz::Array<double,1>&           eforce,
            const blitz::Array<double,1>&     testshape,
            const double                      fac
            ) const
    {
        const std::vector<int> testdof = dofman.LocalDofPosPerField(testfield);
        const int numparamtest = testdof.size();
        for (int vi=0; vi<numparamtest; ++vi)
        {
            const int testpos = testdof[vi];
            eforce(testpos) += fac*testshape(vi);
        }
    };
};

template <DRT::Element::DiscretizationType DISTYPE>
class LocalAssembler<DISTYPE, XFEM::standard_assembly>
{
public:
    //! assemble one expression for the element stiffness matrix for a given test and trial field
    template <XFEM::PHYSICS::Field testfield,
              XFEM::PHYSICS::Field trialfield>
    void Matrix(
            const XFEM::ElementDofManager&    dofman,
            blitz::Array<double,2>&           estif,
            const blitz::Array<double,1>&     testshape,
            const double                      fac,
            const blitz::Array<double,1>&     trialshape
            ) const
    {
        const int numparam = DRT::Utils::getNumberOfElementNodes<DISTYPE>();
        for (int ui=0; ui<numparam; ++ui)
        {
            const int trialpos = 4*ui + getFieldNumber<trialfield>();
            for (int vi=0; vi<numparam; ++vi)
            {
                const int testpos = 4*vi + getFieldNumber<testfield>();
                estif(testpos, trialpos) += fac*testshape(vi)*trialshape(ui) ;
            }
        }
    };
    
    //! assemble one expression for the element rhs for a given test field
    template <XFEM::PHYSICS::Field testfield>
    void Vector(
            const XFEM::ElementDofManager&    dofman,
            blitz::Array<double,1>&           eforce,
            const blitz::Array<double,1>&     testshape,
            const double                      fac
            ) const
    {
        const int numparam = DRT::Utils::getNumberOfElementNodes<DISTYPE>();
        for (int vi=0; vi<numparam; ++vi)
        {
            const int testpos = 4*vi + getFieldNumber<testfield>();
            eforce(testpos) += fac*testshape(vi);
        }
    };
};


//void MatrixFast(
//        const XFEM::ElementDofManager&    dofman,
//        blitz::Array<double,2>&           estif,
//        const XFEM::PHYSICS::Field&       testfield,
//        const XFEM::PHYSICS::Field&       trialfield,
//        const blitz::Array<double,2>&     shapeprod
//        )
//{
//    const vector<int> testdof  = dofman.LocalDofPosPerField(testfield);
//    const vector<int> trialdof = dofman.LocalDofPosPerField(trialfield);
//    
//    const int numparamtest  = testdof.size();
//    const int numparamtrial = trialdof.size();
//
//    for (int ui=0; ui<numparamtrial; ++ui)
//    {
//        const int trialpos = trialdof[ui];
//
//        for (int vi=0; vi<numparamtest; ++vi)
//        {
//            const int testpos = testdof[vi];
//
//            estif(testpos, trialpos) += shapeprod(vi, ui) ;
//        }
//    }
//}

//void integrateVectorBlock(
//        const XFEM::ElementDofManager&    dofman,
//        blitz::Array<double,1>&           eforce,
//        const XFEM::PHYSICS::Field&       testfield,
//        const blitz::Array<double,1>&     testshape
//        )
//{
//
//    const vector<int> testdof  = dofman.LocalDofPosPerField(testfield);
//    
//    const int numparamtest = testdof.size();
//
//    for (int vi=0; vi<numparamtest; ++vi)
//    {
//        const int testpos = testdof[vi];
//
//        eforce(testpos) += testshape(vi);
//    }
//}



  
#endif

#endif
#endif
