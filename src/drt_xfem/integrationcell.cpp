/*!
\file integrationcell.cpp

\brief integration cell

<pre>
Maintainer: Axel Gerstenberger
            gerstenberger@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>
*/

#ifdef D_XFEM
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE

#include "integrationcell.H"
#include <string>
#include <sstream>

using namespace std;
using namespace XFEM;

#define MFOREACH(TYPE,VAL,VALS) for( TYPE::iterator VAL = VALS.begin(); VAL != VALS.end(); ++VAL )
#define MCONST_FOREACH(TYPE,VAL,VALS) for( TYPE::const_iterator VAL = VALS.begin(); VAL != VALS.end(); ++VAL )
#define MPFOREACH(TYPE,VAL,VALS) for( TYPE::const_iterator VAL = VALS->begin(); VAL != VALS->end(); ++VAL )

//
//  ctor
//
Integrationcell::Integrationcell(	const int id, 
									const vector< vector<double> > coordinates) :
id_(id), 
coordinates_(coordinates)
{
    return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor                                                mwgee 11/06|
 *----------------------------------------------------------------------*/
Integrationcell::Integrationcell(const Integrationcell& old) : 
id_(old.id_), 
coordinates_(old.coordinates_)
{
    return;   
}
     
//
//  get coordinates
//
vector< vector<double> >  Integrationcell::GetCoord() const
{
    return coordinates_;   
}

string Integrationcell::Print() const
{
    stringstream s;
    s << "IntCell" << id_ << "\n";
    MCONST_FOREACH(vector< vector<double> >, coordinate, coordinates_)
    {
        s << "[";
        MPFOREACH(vector<double>, val, coordinate)
        {
            s << *val << " ";
        };
        s << "]\n";
    };
    return s.str();
}

#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_XFEM


