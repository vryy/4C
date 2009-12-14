/*!----------------------------------------------------------------------
\file drt_utils.H
\brief A collection of helper methods for namespace DRT

<pre>
-------------------------------------------------------------------------
                 BACI finite element library subsystem
            Copyright (2008) Technical University of Munich

Under terms of contract T004.008.000 there is a non-exclusive license for use
of this work by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library is proprietary software. It must not be published, distributed,
copied or altered in any form or any media without written permission
of the copyright holder. It may be used under terms and conditions of the
above mentioned license by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library may solemnly used in conjunction with the BACI contact library
for purposes described in the above mentioned contract.

This library contains and makes use of software copyrighted by Sandia Corporation
and distributed under LGPL licence. Licensing does not apply to this or any
other third party software used here.

Questions? Contact Dr. Michael W. Gee (gee@lnm.mw.tum.de)
                   or
                   Prof. Dr. Wolfgang A. Wall (wall@lnm.mw.tum.de)

http://www.lnm.mw.tum.de

-------------------------------------------------------------------------
</pre>
<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifndef UTILS_H
#define UTILS_H


#include "Epetra_Export.h"

#include "drt_discret.H"
#include "../drt_fem_general/drt_utils_local_connectivity_matrices.H"

using namespace std;
using namespace Teuchos;




namespace DRT
{
/// collection of useful and important but less central stuff
namespace UTILS
{

  /*!
  \brief Create an instance of a ParObject depending on the type stored in data

  An instance of ParObject is allocated and returned. The type of instance
  depends on the first entry in data and is an int (defined at the top of
  ParObject.H)

  \param data (in): A char vector used for communication or io

  \warning All instances of ParObject have to store the parobject id
           at the vey beginning of the char vector

  \return Allocates and returns the correct instance of ParObject where data is
          already unpacked in the instance. The calling method is responsible for
          freeing this instance!
  */
  DRT::ParObject* Factory(const vector<char>& data);

  /*!
  \brief Create an instance of a finite element depending on the type of element

  A DRT::Element derived class is allocated and returned. The type of element
  allocated depends on the input parameter eletype.

  \param eletype (in): A string containing the type of element
  \param distype (in): A string containing the distype of the element
  \param id      (in): id of the new element to be created
  \param owner   (in): owner of the new element

  */
  RefCountPtr<DRT::Element> Factory(const string eletype,
                                    const string distype,
                                    const int id        ,
                                    const int owner     );

  /*!
  \brief Partition a graph with metis

  Take an input graph with any row map and a vector of point weights with the same row
  map and redistribute the graph using the partitioning library metis.
  The output graph has a row map according to the new distribution.

  \param graph : (in) any graph
  \param weights: (int) a vector of weight factors for each row of graph.
                  Map of weights must be pointwise equal to graph.RowMap().

  \return graph with new distribution
  */
  RefCountPtr<Epetra_CrsGraph> PartGraphUsingMetis(const Epetra_CrsGraph& graph,
                                                   const Epetra_Vector& weights);


#ifdef PARMETIS

  /// redistribute a discretization using ParMetis
  /*!
    Actually, only the \p rownodes and \p colnodes maps are created.
   */
  void PartUsingParMetis(RCP<DRT::Discretization> dis,
                         RCP<Epetra_Map> roweles,
                         RCP<Epetra_Map>& rownodes,
                         RCP<Epetra_Map>& colnodes,
                         vector<int>& nids,
                         int nblock,
                         RCP<Epetra_Comm> comm,
                         Epetra_Time& time,
                         bool outflag);
#else

  /// redistribute a discretization using metis
  /*!
    Actually, only the \p rownodes and \p colnodes maps are created.
   */
  void PartUsingMetis(RCP<Epetra_Map>& rownodes,
                      RCP<Epetra_Map>& colnodes,
                      list<vector<int> >& elementnodes,
                      RCP<Epetra_Comm> comm);

#endif


  /*!
  \brief Locally extract a subset of values from an Epetra_Vector

  Extracts lm.size() values from a distributed epetra vector and stores them into local.
  this is NOT a parallel method, meaning that all values to be extracted on a processor
  must be present in global on that specific processor. This usually means that global
  has to be in column map style.

  \param global (in): global distributed vector with values to be extracted
  \param local (out): vector holding values extracted from global
  \param lm     (in): vector containing global ids to be extracted. Size of lm
                      determines number of values to be extracted.
  */
  void ExtractMyValues(const Epetra_Vector& global, vector<double>& local, const vector<int>& lm);

  void ExtractMyValues(const Epetra_Vector& global, Epetra_SerialDenseVector& local, const vector<int>& lm);

  /// Locally extract a subset of values from a (column)-nodemap-based Epetra_MultiVector
  /*  \author henke
   *  \date 06/09
   */
  void ExtractMyNodeBasedValues
  (
    const DRT::Element* ele,            ///< pointer to current element
    std::vector<double>& local,         ///< local vector on element-level
    const Epetra_MultiVector& global    ///< global (multi) vector
    );


  /// Locally extract a subset of values from a (column)-nodemap-based Epetra_MultiVector
  /*  \author g.bau
   *  \date 08/08
   */
  void ExtractMyNodeBasedValues
  (
    const DRT::Element* ele,                ///< pointer to current element
    Epetra_SerialDenseVector& local,        ///< local vector on element-level
    const RCP<Epetra_MultiVector>& global,  ///< global vector
    const int nsd                           ///< number of space dimensions
    );


  /// Locally extract a subset of values from a (column)-nodemap-based Epetra_MultiVector
  /// and fill a local matrix that has implemented the (.,.) operator
  /*  \author g.bau
   *  \date 04/09
   */
  template <class M>
  void ExtractMyNodeBasedValues(
      const DRT::Element* ele,              ///< pointer to current element
      M& localmatrix,                       ///< local matrix on element-level
      const RCP<Epetra_MultiVector>& global,///< global vector
      const int nsd                         ///< number of space dimensions
      )
  {
    if (global==null) dserror("received a TEUCHOS::null pointer");
    if (nsd > global->NumVectors())
      dserror("Requested %d of %d available columns", nsd,global->NumVectors());
    const int iel = ele->NumNode(); // number of nodes
    if (((int) localmatrix.N())!=iel) dserror("local matrix has wrong number of columns");
    if (((int) localmatrix.M())!=nsd) dserror("local matrix has wrong number of rows");

    for (int i=0; i<nsd; i++)
    {
      // access actual component column of multi-vector
      double* globalcolumn = (*global)[i];
      // loop over the element nodes
      for (int j=0;j<iel;j++)
      {
        const int nodegid = (ele->Nodes()[j])->Id();
        const int lid = global->Map().LID(nodegid);
        if (lid<0) dserror("Proc %d: Cannot find gid=%d in Epetra_Vector",(*global).Comm().MyPID(),nodegid);
        localmatrix(i,j) = globalcolumn[lid];
      }
    }
    return;
  }

  /*!
  \brief extract local values from global node-based (multi) vector

  This function returns a column vector!

  \author henke
 */
  template <class M>
  void ExtractMyNodeBasedValues(
    const DRT::Element* ele,
    M& local,
    const Epetra_MultiVector& global)
{
  const int numnode = ele->NumNode();
  const int numcol = global.NumVectors();
  if (((int) local.N())!=1) dserror("local matrix must have one column");
  if (((int) local.M())!=numnode) dserror("local matrix has wrong number of rows");

  // loop over element nodes
  for (int i=0; i<numnode; ++i)
  {
    const int nodegid = (ele->Nodes()[i])->Id();
    const int lid = global.Map().LID(nodegid);
    if (lid<0) dserror("Proc %d: Cannot find gid=%d in Epetra_Vector",global.Comm().MyPID(),nodegid);

    // loop over multi vector columns (numcol=1 for Epetra_Vector)
    for (int col=0; col<numcol; col++)
    {
      double* globalcolumn = (global)[col];
      local((col+(numcol*i)),0) = globalcolumn[lid];
    }
  }
  return;
}

  //! flag, whether surfaces or lines have to be created in the ElementBoundaryFactory
  enum BoundaryBuildType
  {
    buildSurfaces,  ///< build surfaces
    buildLines      ///< build lines
  };

  /*!
   * create new instances of volume / surface / line elements for a given parent element
   *
   * template version
   *
   * \tparam BoundaryEle class name of desired volume/surface/line element, e.g. Fluid3Surface
   * \tparam ParentEle   class name of parent element, e.g. Fluid3
   *
   * this templated function creates all volume / surface / line elements for a given element
   * and fills handed over vectors with RCPs and raw pointers of these boundary elements
   * This method is a very powerful helper function for implementing the necessary
   * Surfaces() and Lines() methods for every element class (especially in 3D)
   * using the BACI conventions for element connectivity
   *
   * \return boundaryeles   vector filled with RCPs of allocated boundary elements
   *
   * \author gjb
   * \date 05/08
   */
  template <class BoundaryEle, class ParentEle>
  vector<RCP<DRT::Element> > ElementBoundaryFactory(
      const BoundaryBuildType  buildtype,  ///< flag, whether volumes, surfaces or lines have to be created
      ParentEle*               ele         ///< pointer on the parent element
      )
  {
    // do we have to build volume, surface or line elements?
    // get node connectivity for specific distype of parent element
    unsigned int nele=0;
    const DRT::Element::DiscretizationType distype = ele->Shape();
    vector< vector<int> > connectivity;
    switch (buildtype)
    {
      case buildSurfaces:
      {
        nele = ele->NumSurface();
        connectivity = DRT::UTILS::getEleNodeNumberingSurfaces(distype);
        break;
      }
      case buildLines:
      {
        nele = ele->NumLine();
        connectivity = DRT::UTILS::getEleNodeNumberingLines(distype);
        break;
      }
    }
    // create vectors that will contain the volume, surface or line elements
    vector<RCP<DRT::Element> > boundaryeles(nele);

    // does DRT::UTILS convention match your implementation of NumSurface() or NumLine()?
    if (nele != connectivity.size()) dserror("number of surfaces or lines does not match!");

    // now, build the new surface/line elements
    for (unsigned int iele = 0; iele < nele; iele++)
    {
      // allocate node vectors
      unsigned int nnode = connectivity[iele].size(); // this number changes for pyramids or wedges
      vector<int> nodeids(nnode);
      vector<DRT::Node*> nodes(nnode);

      // get connectivity infos
      for (unsigned int inode=0;inode<nnode;inode++)
      {
        nodeids[inode] = ele->NodeIds()[connectivity[iele][inode]];
        nodes  [inode] = ele->Nodes  ()[connectivity[iele][inode]];
      }

      // allocate a new boundary element
      boundaryeles[iele] = rcp(new BoundaryEle(iele          ,
                                               ele->Owner()  ,
                                               nodeids.size(),
                                               &nodeids[0]   ,
                                               &nodes[0]     ,
                                               ele           ,
                                               iele          ));

    }

    return boundaryeles;
  }

  //! Given a nodal row and coloumn map, this routines creates the corresponding element maps and call FillComplete()
  void RedistributeWithNewNodalDistribution(
    DRT::Discretization&     dis,
    const Epetra_Map&        noderowmap,
    const Epetra_Map&        nodecolmap
    );

  /// print information about parallel mesh distribution
  void PrintParallelDistribution(const DRT::Discretization& dis);

  /*!
  \brief Get a column vector made of a row vector.

  Using this method, a reference to a vector is returned.
  If the vector is supplied in dis.DofColMap() the vector itself will be returned.
  If the vector is NOT supplied in dis.DofColMap(), but in dis.DofRowMap(),
  a new vector with column map is allocated and the supplied vector is exported to it and returned.

  Note that the very same functionality is used in Discretization::SetState()!

  \param name (in): discretization
  \param state (in): vector of some data
  */
  Teuchos::RCP<const Epetra_Vector> GetColVersionOfRowVector(
      const Teuchos::RCP<const DRT::Discretization> dis,
      const Teuchos::RCP<const Epetra_Vector> state);

} // namespace UTILS
} // namespace DRT


#endif  // #ifndef UTILS_H
#endif  // #ifdef CCADISCRET




