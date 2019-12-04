/*----------------------------------------------------------------------*/
/*! \file

\brief A collection of communication methods for namespace LINALG

\level 0
\maintainer Martin Kronbichler
*/
/*----------------------------------------------------------------------*/

#ifndef SRC_LINALG_LINALG_UTILS_DENSEMATRIX_COMMUNICATION_H_
#define SRC_LINALG_LINALG_UTILS_DENSEMATRIX_COMMUNICATION_H_

#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_CrsGraph.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Export.h"
#include "Epetra_Import.h"
#include "Teuchos_RCP.hpp"
#include "../drt_lib/drt_exporter.H"
#include "linalg_blocksparsematrix.H"

namespace LINALG
{
  // forward declaration
  // class Solver;

  /*!
   \brief Gather information of type vector<T> on a subset of processors

   This template gathers information provided in sdata on a subset of
   processors tprocs, where the length of the array tprocs is ntargetprocs.
   The redistributed data is returned in rdata which has appropiate size
   on output (size of rdata on input is arbitrary). ntargetprocs can be
   one to reduce data to one proc, it also can be equal to the total number
   of processors to make sdata redundant on all procs.

   \note Functionality of this method is equal to that of Epetra_Comm::GatherAll
   except for that the Epetra version demands the data to be of constant
   size over all procs which this method does not require!

   \param sdata (in) : Information to be gathered on tprocs.
   Length of sdata can be different on every proc.
   \param rdata (out): Information from sdata gathered on a subset of procs.
   size of rdata on input is arbitrary, it is exact on output.
   \param ntargetprocs (in): length of tprocs
   \param tprocs (in): vector of procs ids the information in sdata shall be
   gathered on.
   \param comm (in):   communicator to be used.


   */
  template <typename T>
  void Gather(std::vector<T>& sdata, std::vector<T>& rdata, const int ntargetprocs,
      const int* tprocs, const Epetra_Comm& comm)
  {
    const int myrank = comm.MyPID();
    const int numproc = comm.NumProc();
    if (numproc == 1)
    {
      rdata = sdata;
      return;  // nothing to do in serial
    }
    // build a map of data
    std::map<int, std::vector<T>> datamap;
    datamap[myrank] = sdata;
    // build a source map
    Epetra_Map source(numproc, 1, &myrank, 0, comm);
    // build a target map which is redundant on all target procs and zero everywhere else
    bool iamtarget = false;
    for (int i = 0; i < ntargetprocs; ++i)
      if (tprocs[i] == myrank)
      {
        iamtarget = true;
        break;
      }
    std::vector<int> targetvec(0);
    if (iamtarget)
    {
      targetvec.resize(numproc);
      for (int i = 0; i < numproc; ++i) targetvec[i] = i;
    }
    const int tnummyelements = (int)targetvec.size();
    Epetra_Map target(-1, tnummyelements, &targetvec[0], 0, comm);
    // build an exporter and export data
    DRT::Exporter exporter(source, target, comm);
    exporter.Export(datamap);
    // put data from map in rdata
    rdata.clear();
    int count = 0;
    typename std::map<int, std::vector<T>>::const_iterator curr;
    for (curr = datamap.begin(); curr != datamap.end(); ++curr)
    {
      const std::vector<T>& current = curr->second;
      const int size = (int)current.size();
      rdata.resize((int)rdata.size() + size);
      for (int i = 0; i < size; ++i) rdata[count + i] = current[i];
      count += size;
    }
    return;
  }

  /*!
   \brief Gather information of type set<T> on a subset of processors

   This template gathers information provided in sdata on a subset of
   processors tprocs, where the length of the array tprocs is ntargetprocs.
   The redistributed data is returned in rdata which has appropiate size
   on output (size of rdata on input is arbitrary). ntargetprocs can be
   one to reduce data to one proc, it also can be equal to the total number
   of processors to make sdata redundant on all procs.

   \note Functionality of this method is equal to that of Epetra_Comm::GatherAll
   except for that the Epetra version demands the data to be of constant
   size over all procs which this method does not require!

   \param sdata (in) : Information to be gathered on tprocs.
   Length of sdata can be different on every proc.
   \param rdata (out): Information from sdata gathered on a subset of procs.
   size of rdata on input is arbitrary, it is exact on output.
   \param ntargetprocs (in): length of tprocs
   \param tprocs (in): vector of procs ids the information in sdata shall be
   gathered on.
   \param comm (in):   communicator to be used.


   */
  template <typename T>
  void Gather(std::set<T>& sdata, std::set<T>& rdata, const int ntargetprocs, const int* tprocs,
      const Epetra_Comm& comm)
  {
    const int myrank = comm.MyPID();
    const int numproc = comm.NumProc();
    if (numproc == 1)
    {
      rdata = sdata;
      return;  // nothing to do in serial
    }
    // build a map of data
    std::map<int, std::set<T>> datamap;
    datamap[myrank] = sdata;
    // build a source map
    Epetra_Map source(numproc, 1, &myrank, 0, comm);
    // build a target map which is redundant on all target procs and zero everywhere else
    bool iamtarget = false;
    for (int i = 0; i < ntargetprocs; ++i)
      if (tprocs[i] == myrank)
      {
        iamtarget = true;
        break;
      }
    std::vector<int> targetvec(0);
    if (iamtarget)
    {
      targetvec.resize(numproc);
      for (int i = 0; i < numproc; ++i) targetvec[i] = i;
    }
    const int tnummyelements = (int)targetvec.size();
    Epetra_Map target(-1, tnummyelements, &targetvec[0], 0, comm);
    // build an exporter and export data
    DRT::Exporter exporter(source, target, comm);
    exporter.Export(datamap);
    // put data from map in rdata
    rdata.clear();
    typename std::map<int, std::set<T>>::const_iterator curr;
    for (curr = datamap.begin(); curr != datamap.end(); ++curr)
    {
      const std::set<T>& current = curr->second;
      typename std::set<T>::const_iterator setiter;
      for (setiter = current.begin(); setiter != current.end(); ++setiter)
      {
        rdata.insert(*setiter);
      }
    }
    return;
  }

  /*!
   \brief Gather information of type std::map<int, std::set<T> > on a subset of processors

   This template gathers information provided in sdata on a subset of
   processors tprocs, where the length of the array tprocs is ntargetprocs.
   The redistributed data is returned in rdata which has appropiate size
   on output (size of rdata on input is arbitrary). ntargetprocs can be
   one to reduce data to one proc, it also can be equal to the total number
   of processors to make sdata redundant on all procs.

   \note Functionality of this method is equal to that of Epetra_Comm::GatherAll
   except for that the Epetra version demands the data to be of constant
   size over all procs which this method does not require!

   \param sdata (in) : Information to be gathered on tprocs.
   Length of sdata can be different on every proc.
   \param rdata (out): Information from sdata gathered on a subset of procs.
   size of rdata on input is arbitrary, it is exact on output.
   \param ntargetprocs (in): length of tprocs
   \param tprocs (in): vector of procs ids the information in sdata shall be
   gathered on.
   \param comm (in):   communicator to be used.


   */
  template <typename T>
  void Gather(std::map<int, std::set<T>>& sdata, std::map<int, std::set<T>>& rdata,
      const int ntargetprocs, const int* tprocs, const Epetra_Comm& comm)
  {
    const int myrank = comm.MyPID();
    const int numproc = comm.NumProc();
    if (numproc == 1)
    {
      rdata = sdata;
      return;  // nothing to do in serial
    }
    // build a map of data
    std::map<int, std::map<int, std::set<T>>> datamap;
    datamap[myrank] = sdata;
    // build a source map
    Epetra_Map source(numproc, 1, &myrank, 0, comm);
    // build a target map which is redundant on all target procs and zero everywhere else
    bool iamtarget = false;
    for (int i = 0; i < ntargetprocs; ++i)
      if (tprocs[i] == myrank)
      {
        iamtarget = true;
        break;
      }
    std::vector<int> targetvec(0);
    if (iamtarget)
    {
      targetvec.resize(numproc);
      for (int i = 0; i < numproc; ++i) targetvec[i] = i;
    }
    const int tnummyelements = (int)targetvec.size();
    Epetra_Map target(-1, tnummyelements, &targetvec[0], 0, comm);
    // build an exporter and export data
    DRT::Exporter exporter(source, target, comm);
    exporter.Export(datamap);
    // put data from map in rdata
    rdata.clear();
    typename std::map<int, std::map<int, std::set<T>>>::const_iterator curr;
    for (curr = datamap.begin(); curr != datamap.end(); ++curr)
    {
      const std::map<int, std::set<T>>& innercurr = curr->second;
      typename std::map<int, std::set<T>>::const_iterator inneriter;
      for (inneriter = innercurr.begin(); inneriter != innercurr.end(); ++inneriter)
      {
        rdata[inneriter->first].insert(inneriter->second.begin(), inneriter->second.end());
      }
    }
    return;
  }

  /*!
   \brief Gather information of type std::map<int, std::vector<T> > on a subset of processors

   This template gathers information provided in sdata on a subset of
   processors tprocs, where the length of the array tprocs is ntargetprocs.
   The redistributed data is returned in rdata which has appropiate size
   on output (size of rdata on input is arbitrary). ntargetprocs can be
   one to reduce data to one proc, it also can be equal to the total number
   of processors to make sdata redundant on all procs.

   \note Functionality of this method is equal to that of Epetra_Comm::GatherAll
   except for that the Epetra version demands the data to be of constant
   size over all procs which this method does not require!

   \param sdata (in) : Information to be gathered on tprocs.
   Length of sdata can be different on every proc.
   \param rdata (out): Information from sdata gathered on a subset of procs.
   size of rdata on input is arbitrary, it is exact on output.
   \param ntargetprocs (in): length of tprocs
   \param tprocs (in): vector of procs ids the information in sdata shall be
   gathered on.
   \param comm (in):   communicator to be used.


   */
  template <typename T>
  void Gather(std::map<int, std::vector<T>>& sdata, std::map<int, std::vector<T>>& rdata,
      const int ntargetprocs, const int* tprocs, const Epetra_Comm& comm)
  {
    const int myrank = comm.MyPID();
    const int numproc = comm.NumProc();
    if (numproc == 1)
    {
      rdata = sdata;
      return;  // nothing to do in serial
    }
    // build a map of data
    std::map<int, std::map<int, std::vector<T>>> datamap;
    datamap[myrank] = sdata;
    // build a source map
    Epetra_Map source(numproc, 1, &myrank, 0, comm);
    // build a target map which is redundant on all target procs and zero everywhere else
    bool iamtarget = false;
    for (int i = 0; i < ntargetprocs; ++i)
      if (tprocs[i] == myrank)
      {
        iamtarget = true;
        break;
      }
    std::vector<int> targetvec(0);
    if (iamtarget)
    {
      targetvec.resize(numproc);
      for (int i = 0; i < numproc; ++i) targetvec[i] = i;
    }
    const int tnummyelements = (int)targetvec.size();
    Epetra_Map target(-1, tnummyelements, &targetvec[0], 0, comm);
    // build an exporter and export data
    DRT::Exporter exporter(source, target, comm);
    exporter.Export(datamap);
    // put data from map in rdata
    rdata.clear();
    typename std::map<int, std::map<int, std::vector<T>>>::const_iterator curr1;
    typename std::map<int, std::vector<T>>::const_iterator curr2;
    typename std::vector<T>::const_iterator curr3;
    for (curr1 = datamap.begin(); curr1 != datamap.end(); ++curr1)
    {
      const std::map<int, std::vector<T>>& data = curr1->second;
      for (curr2 = data.begin(); curr2 != data.end(); ++curr2)
      {
        const std::vector<T>& vectordata = curr2->second;
        for (curr3 = vectordata.begin(); curr3 != vectordata.end(); ++curr3)
        {
          rdata[curr2->first].push_back(*curr3);
        }
      }
    }
    return;
  }

  /*!
   \brief Gather information of type map<T,U> on a subset of processors

   This template gathers information provided in sdata on a subset of
   processors tprocs, where the length of the array tprocs is ntargetprocs.
   The redistributed data is returned in rdata which has appropiate size
   on output (size of rdata on input is arbitrary). ntargetprocs can be
   one to reduce data to one proc, it also can be equal to the total number
   of processors to make sdata redundant on all procs.

   \param sdata (in) : Information to be gathered on tprocs.
   Length of sdata can be different on every proc.
   \param rdata (out): Information from sdata gathered on a subset of procs.
   size of rdata on input is arbitrary, it is exact on output.
   \param ntargetprocs (in): length of tprocs
   \param tprocs (in): vector of procs ids the information in sdata shall be
   gathered on.
   \param comm (in):   communicator to be used.

   */

  template <typename T, typename U>
  void Gather(std::map<T, U>& sdata, std::map<T, U>& rdata, const int ntargetprocs,
      const int* tprocs, const Epetra_Comm& comm)
  {
    const int myrank = comm.MyPID();
    const int numproc = comm.NumProc();
    if (numproc == 1)
    {
      rdata = sdata;
      return;  // nothing to do in serial
    }

    // build a map of data
    std::map<int, std::map<T, U>> datamap;
    datamap[myrank] = sdata;

    // build a source map
    Epetra_Map source(numproc, 1, &myrank, 0, comm);
    // build a target map which is redundant on all target procs and zero everywhere else
    bool iamtarget = false;
    for (int i = 0; i < ntargetprocs; ++i)
      if (tprocs[i] == myrank)
      {
        iamtarget = true;
        break;
      }
    std::vector<int> targetvec(0);
    if (iamtarget)
    {
      targetvec.resize(numproc);
      for (int i = 0; i < numproc; ++i) targetvec[i] = i;
    }
    const int tnummyelements = (int)targetvec.size();
    Epetra_Map target(-1, tnummyelements, &targetvec[0], 0, comm);
    // build an exporter and export data
    DRT::Exporter exporter(source, target, comm);
    exporter.Export(datamap);
    // put data from map in rdata
    rdata.clear();
    typename std::map<int, std::map<T, U>>::const_iterator curr;
    for (curr = datamap.begin(); curr != datamap.end(); ++curr)
    {
      const std::map<T, U>& current = curr->second;
      typename std::map<T, U>::const_iterator mapiter;
      for (mapiter = current.begin(); mapiter != current.end(); ++mapiter)
      {
        rdata.insert(std::make_pair(mapiter->first, mapiter->second));
      }
    }
    return;
  }

  /*!
   \brief Gather information of type set<T> from all processors

   This template gathers information provided in data on all processors.
   The redistributed data is returned in data which has appropriate size
   on output.

   \note Functionality of this method is equal to that of Epetra_Comm::GatherAll
   except for that the Epetra version demands the data to be of constant
   size over all procs which this method does not require!

   \param data (in/out) : Information to be gathered.
   Length of data can be different on every proc.
   \param comm (in):   communicator to be used.

   */
  template <typename T>
  void GatherAll(std::set<T>& data, const Epetra_Comm& comm)
  {
    // ntargetprocs is equal to the total number of processors to make data redundant on all procs
    const int numprocs = comm.NumProc();
    int allproc[numprocs];
    for (int i = 0; i < numprocs; ++i) allproc[i] = i;

    Gather<T>(data, data, numprocs, allproc, comm);
    return;
  }

  // nagler 07/2012
  /*!
   \brief Gather information of type map<T,U> from all processors

   This template gathers information provided in data on all processors.
   The redistributed data is returned in data which has appropriate size
   on output.

   \param data (in/out) : Information to be gathered.
   Length of data can be different on every proc.
   \param comm (in):   communicator to be used.

   */

  template <typename T, typename U>
  void GatherAll(std::map<T, U>& data, const Epetra_Comm& comm)
  {
    const int numprocs = comm.NumProc();
    int allproc[numprocs];
    for (int i = 0; i < numprocs; ++i) allproc[i] = i;

    Gather<T, U>(data, data, numprocs, allproc, comm);

    return;
  }

  //                                                              sudhakar 02/2014
  /*!
   \brief Gather information of type map<int,vector<T> > from all processors

   This template gathers information provided in data on all processors.
   The redistributed data is returned in data which has appropriate size
   on output.

   \param data (in/out) : Information to be gathered.
   Length of data can be different on every proc.
   \param comm (in):   communicator to be used.

   */

  template <typename T>
  void GatherAll(std::map<int, std::vector<T>>& data, const Epetra_Comm& comm)
  {
    const int numprocs = comm.NumProc();
    int allproc[numprocs];
    for (int i = 0; i < numprocs; ++i) allproc[i] = i;

    Gather<T>(data, data, numprocs, allproc, comm);

    return;
  }

  /*!
   \brief Gather information of type vector<T> from all processors

   This template gathers information provided in data on all processors.
   The redistributed data is returned in data which has appropriate size
   on output.

   \note Functionality of this method is equal to that of Epetra_Comm::GatherAll
   except for that the Epetra version demands the data to be of constant
   size over all procs which this method does not require!

   \param data (in/out) : Information to be gathered.
   Length of data can be different on every proc.
   \param comm (in):   communicator to be used.

   */
  template <typename T>
  void GatherAll(std::vector<T>& data, const Epetra_Comm& comm)
  {
    // ntargetprocs is equal to the total number of processors to make data redundant on all procs
    const int numprocs = comm.NumProc();
    int allproc[numprocs];
    for (int i = 0; i < numprocs; ++i) allproc[i] = i;

    Gather<T>(data, data, numprocs, allproc, comm);
    return;
  }

  /*!
   \brief Create an allreduced vector of gids from the given Epetra_Map

   We have nodes and elements arbitrary global ids. On rare occasions, however,
   we need to allreduce a particular map to one or more processors.
   This is a building block for such occasions. We allreduce the gids
   of the given Epetra_Map into a vector ordered by processor number.

   \note You are not supposed to use redundant vectors in normal
   situations. If you happen to need this method you are probably
   about to do something illegal.

   \param rredundant (o) redundant vector of global ids
   \param emap (i) unique distributed Epetra_Map

   \author u.kue
   \date 05/07
   */
  void AllreduceEMap(std::vector<int>& rredundant, const Epetra_Map& emap);

  /// Create an allreduced gid to index map from the given Epetra_Map
  /*!
   We have nodes and elements with unique but otherwise arbitrary
   global ids. But unfortunately we need an allreduced vector of dof
   numbers during the dof assignment phase. In order to use such a
   vector we need to map from global ids to vector indexes. Here we
   provide that map.

   \note You are not supposed to use redundant vectors in normal
   situations. If you happen to need this method you are probably
   about to do something illegal.

   \param idxmap (o) map from global ids to (redundant) vector indexes
   \param emap (i) unique distributed Epetra_Map

   \author u.kue
   \date 05/07
   */
  void AllreduceEMap(std::map<int, int>& idxmap, const Epetra_Map& emap);

  /*!
   \brief Create an allreduced gid to index map from the given Epetra_Map
          on a distinct processor, all other procs create empty maps instead.

   This method is currently used within the parallel post_drt_ensight
   filter in order to import all values stored in a distributed Epetra_Vector
   to processor 0 for writing them into file.

   \note see also documentation for the usual AllreduceEMap methods

   \param emap (i) any distributed Epetra_Map
   \param pid (i)  processor id where you want to have the allreduced map
   exclusively
   \author gjb
   \date 11/07
   */
  Teuchos::RCP<Epetra_Map> AllreduceEMap(const Epetra_Map& emap, const int pid);

  /*!
   \brief Create an allreduced Epetra_Map from the given Epetra_Map
          and give it to all processors.

   This method is currently used within the constraint management, since
   current values of constraint values and langrange multipliers are distributed
   uniquely for computation. At some places we need the full information of these
   values on every processor, so this method has to be used.

   \note You are not supposed to use redundant vectors in normal
   situations. If you happen to need this method you are probably
   about to do something illegal.

   \param emap (i) any distributed Epetra_Map

   \author tk
   \date 04/08
   */
  Teuchos::RCP<Epetra_Map> AllreduceEMap(const Epetra_Map& emap);

  /*!
   \brief Create an allreduced Epetra_Map from the given Epetra_Map
          and give it to all processors.

   Here, we have a overlapping source map and still want to have a fully
   redundant map on all processors without duplicated entries.

   \author u.kue
   \date 08/09
   */
  Teuchos::RCP<Epetra_Map> AllreduceOverlappingEMap(const Epetra_Map& emap);

  /*!
   \brief Create an allreduced Epetra_Map from the given Epetra_Map
          on a distinct processor, all other procs create empty maps instead.

   \param emap (i) any distributed overlapping Epetra_Map
   \param pid (i)  processor id where you want to have the allreduced and sorted map exclusively

   \author ghamm
   \date 10/14
   */
  Teuchos::RCP<Epetra_Map> AllreduceOverlappingEMap(const Epetra_Map& emap, const int pid);

  /*!
   \brief Find position of my map elements in a consecutive vector

   The idea is to put the entries of a given map into a redundant
   vector, ordered by processor number. The map is assumed to be
   nonoverlapping. Here we figure out the index of our first entry in
   that vector.

   \note You are not supposed to use redundant vectors in normal
   situations. If you happen to need this method you are probably
   about to do something illegal.

   \param nummyelements (i) number of elements on this proc
   \param comm (i) communicator

   \return vector position of first entry on each processor

   \author u.kue
   \date 05/07
   */
  int FindMyPos(int nummyelements, const Epetra_Comm& comm);

  /// create an allreduced sorted copy of the source vectors
  void AllreduceVector(
      const std::vector<int>& src, std::vector<int>& dest, const Epetra_Comm& comm);

  /*!
   \brief Communication between all pairs of processes, with distinct data for each.

   Sends a different vector<int> to each processes. The size of each vector may
   be different, zero-length vectors are allowed.
   Communication is implemented with the MPI function MPI_Alltoallv.

   \param comm (i) communicator
   \param send (i) vector of length comm.NumProc(), j-th element to be send to j-th processor.
   \param recv (o) vector of length comm.NumProc(), j-th element received from j-th processor.

   \author h.kue
   \date 09/07
   */
  void AllToAllCommunication(const Epetra_Comm& comm, const std::vector<std::vector<int>>& send,
      std::vector<std::vector<int>>& recv);

  /*!
   \brief Return the first slot of a pair

   To be used with stl algorithms.

   This should be part of stl but is not. So we define our own version.
   */
  template <typename pair_type>
  struct select1st
      : public std::unary_function<const pair_type&, const typename pair_type::first_type&>
  {
    const typename pair_type::first_type& operator()(const pair_type& v) const { return v.first; }
  };

}  // namespace LINALG

#endif /* SRC_LINALG_LINALG_UTILS_DENSEMATRIX_COMMUNICATION_H_ */
