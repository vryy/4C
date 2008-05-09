/*----------------------------------------------------------------------*/
/*!
\file pre_exodus_writedat.H

\brief pre_exodus .dat-file writer 

<pre>
Maintainer: Moritz
            frenzel@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/frenzel
            089 - 289-15240
</pre>

Here is everything related with writing a dat-file
*/
/*----------------------------------------------------------------------*/
#ifdef D_EXODUS
#ifndef PRE_EXODUS_WRITEDAT_H_
#define PRE_EXODUS_WRITEDAT_H_

#include <string>
#include <iostream>
#include <vector>
#include <Teuchos_RefCountPtr.hpp>
#include "../drt_lib/drt_dserror.H"
#include "pre_exodus_reader.H"
#include "pre_exodus_readbc.H"

using namespace std;
using namespace Teuchos;

namespace EXODUS
{
//! write datfile
int WriteDatFile(const string datfile, const EXODUS::Mesh mymesh,
    const string headfile, const vector<map<int,EXODUS::bc_entity> > bc_specs);

//! datfile Intro
string WriteDatIntro(const EXODUS::Mesh& mymesh);

//! headfile content into datfile
string WriteDatHead(const string& headfile);

//! design description
string WriteDatDesign(const vector<EXODUS::bc_entity>& baciconds);

//! conditions
string WriteDatConditions(const vector<EXODUS::bc_entity>& baciconds);

//! DesignNode - Node Topology
string WriteDatDesignTopology(const vector<EXODUS::bc_entity>& baciconds, const EXODUS::Mesh& mymesh);

//! return set of nodes in bc_entity
const set<int> GetNsFromBCEntity(const EXODUS::bc_entity& e, const EXODUS::Mesh& m);

//! Nodes into datfile
string WriteDatNodes(const EXODUS::Mesh& mymesh);

//! Elements into datfile
string WriteDatEles(const vector<map<int,EXODUS::bc_entity> >& bc_specs, const EXODUS::Mesh& mymesh);

//! Elements from eblock
string DatEles(const EXODUS::ElementBlock& eb, const EXODUS::bc_entity& acte, int& struele);

//! as e.g. also ElementBlocks can hold 'BACI-Conditions' we find everything what will be a 'BACI-Condition"
vector<EXODUS::bc_entity> FindBACIConditions(const vector<map<int,EXODUS::bc_entity> >& bc_specs);

} // of namespace EXODUS


#endif /*PRE_EXODUS_WRITEDAT_H_*/
#endif /*EXODUS*/

