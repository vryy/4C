/*!----------------------------------------------------------------------
\file drt_exporter.H

\brief Management of explicit MPI communications

<pre>
\level 0

\maintainer Martin Kronbichler
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>

*----------------------------------------------------------------------*/
#ifndef EXPORTER_H
#define EXPORTER_H

#include <vector>
#include <map>
#include <set>
#include <Epetra_Map.h>
#include <Epetra_Comm.h>
#include <Epetra_MpiComm.h>
#include <Teuchos_RCP.hpp>

#include "drt_node.H"
#include "drt_utils_factory.H"
#include "drt_parobject.H"


namespace DRT
{
  /*!
  \brief A class to manage explicit mpi communications

  The discretization management module uses this class to do most of its
  communication. It is used to redistribute grids and to do point-to-point
  communications. It is therefore the only place on DRT where explicit calls to MPI
  methods are done.<br>
  It has strong capabilities in gathering and scattering information in a collective AND
  an individual way. Whenever you need explicit communication, check this class first before
  implementing your own mpi stuff.

  \author gee (gee@lnm.mw.tum.de)
  */
  class Exporter
  {
    class ExporterHelper;

   public:
    /*!
    \brief Standard Constructor

    this ctor constructs an exporter with no maps. It can than be used to do
    point-to-point communication only, map based exportes are not possible!

    \param comm    (in): Communicator that shall be used in exports
    */
    Exporter(const Epetra_Comm& comm);

    /*!
    \brief Standard Constructor

    \param frommap (in): The source map data shall be exported from
    \param tomap   (in): The target map data shall be exported to
    \param comm    (in): Communicator that shall be used in exports
    */
    Exporter(const Epetra_Map& frommap, const Epetra_Map& tomap, const Epetra_Comm& comm);

    /*!
    \brief Copy Constructor

    */
    Exporter(const DRT::Exporter& old);

    /*!
    \brief Destructor

    */
    virtual ~Exporter();


    //! @name Acess methods

    /*!
    \brief Get communicator
    */
    inline const Epetra_Comm& Comm() const { return comm_; }

    /*!
    \brief Get source map
    */
    inline const Epetra_Map& SourceMap() const { return frommap_; }

    /*!
    \brief Get target map
    */
    inline const Epetra_Map& TargetMap() const { return tomap_; }

    //@}

    //! @name Communication methods

    /*!
    \brief Communicate a map of objects that implement ParObject

    This method takes a map of objects and redistributes them according to
    the send and receive plans. It is implicitly assumed, that the key in
    the map of objects pointwise matches SourceMap(). It is also assumed
    (and tested), that type T implements the ParObject class.

    \param parobjects (in/out): A map of classes T that implement the
                                class ParObject. On input, the map
                                has a distribution matching SourceMap().
                                On output, the map has a distribution of
                                TargetMap().
    */
    template <typename T>
    void Export(std::map<int, Teuchos::RCP<T>>& parobjects);

    /*!
    \brief Communicate a map of vectors of some basic data type T

    This method takes a map of vectors and redistributes them according to
    the send and receive plans. It is implicitly assumed, that the key in
    the map of vectors pointwise matches SourceMap().

    \note T can be int, double or char. The method will not compile will other
          then these basic data types (and will give a rater kryptic error message)

    \param data (in/out): A map of vectors<T>. On input, the map
                          has a distribution matching SourceMap().
                          On output, the map has a distribution of
                          TargetMap().
    */
    template <typename T>
    void Export(std::map<int, std::vector<T>>& data);

    /*!
    \brief Communicate a map of sets of some basic data type T (currently only int)

    This method takes a map of sets and redistributes them according to
    the send and receive plans. It is implicitly assumed, that the key in
    the map of sets pointwise matches SourceMap().

    \note T can be int. The method will not compile will other
          then these basic data types (and will give a rater kryptic error message)

    \param data (in/out): A map of sets<T>. On input, the map
                          has a distribution matching SourceMap().
                          On output, the map has a distribution of
                          TargetMap().
    */
    template <typename T>
    void Export(std::map<int, std::set<T>>& data);

    //                                                             nagler 07/2012
    /*!
    \brief Communicate a map of sets of some basic data type T and U

    This method takes a map of sets and redistributes them according to
    the send and receive plans. It is implicitly assumed, that the key in
    the map of sets pointwise matches SourceMap().

    \note T and U can have arbitrarily types

    \param data (in/out): A map of maps<T,U>. On input, the map
                          has a distribution matching SourceMap().
                          On output, the map has a distribution of
                          TargetMap().
    */
    template <typename T, typename U>
    void Export(std::map<int, std::map<T, U>>& data);

    /*!
    \brief Communicate a map of int values

    This method takes a map of ints and redistributes them according to
    the send and receive plans. It is implicitly assumed, that the key in
    the map of objects pointwise matches SourceMap().

    \param parobjects (in/out): A map of ints. On input, the map
                                has a distribution matching SourceMap().
                                On output, the map has a distribution of
                                TargetMap().
    */
    void Export(std::map<int, int>& data);

    /*!
    \brief Communicate a map of double values

    This method takes a map of doubles and redistributes them according to
    the send and receive plans. It is implicitly assumed, that the key in
    the map of objects pointwise matches SourceMap().

    \param parobjects (in/out): A map of doubles. On input, the map
                                has a distribution matching SourceMap().
                                On output, the map has a distribution of
                                TargetMap().
    */
    void Export(std::map<int, double>& data);

    /*!
    \brief Communicate a map of serial dense matrices

    This method takes a map of serial dense matrices and redistributes them
    according to the send and receive plans. It is implicitly assumed, that the
    key in the map of objects pointwise matches SourceMap().

    \param parobjects (in/out): A map of serial dense matrices. On input, the map
                                has a distribution matching SourceMap().
                                On output, the map has a distribution of
                                TargetMap().
    */
    void Export(std::map<int, Teuchos::RCP<Epetra_SerialDenseMatrix>>& data);

    /*!
    \brief Send data from one processor to another (nonblocking)

    The method will send an array of chars in a nonblocking way meaning that this method
    will return immediately on the calling processor - even if communcation has not
    finished yet. The char array must not be altered or destroyed as long as the
    communication might still be in progress. This can be tested for using DRT::Exporter::Wait
    and the request handle returned.
    The receiving processor should call DRT::Exporter::ReceiveAny to receive the message.
    Note that messages from one explicit proc to another explicit proc are non-overtaking meaning
    they will arrive in the order they have been sent.

    \note This is an individual call.

    \param frompid (in)  : sending processors' pid
    \param topid (in)    : target processors' pid
    \param data (in)     : ptr to data to be send
    \param dsize (in)    : size of data (no. of chars)
    \param tag (in)      : tag to be used with message
    \param request (out) : mpi request handle to be used for testing completion of the
                           communication. data may not be altered or destroyed before
                           communcation finalized! One can use DRT::Exporter::Wait for this.

    \note This is an individual call
    */
    void ISend(const int frompid, const int topid, const char* data, const int dsize, const int tag,
        MPI_Request& request);

    /*!
    \brief Send data from one processor to another (nonblocking)

    The method will send an array of ints in a nonblocking way meaning that this method
    will return immediately on the calling processor - even if communcation has not
    finished yet. The int array must not be altered or destroyed as long as the
    communication might still be in progress. This can be tested for using DRT::Exporter::Wait
    and the request handle returned.
    The receiving processor should call DRT::Exporter::ReceiveAny to receive the message.
    Note that messages from one explicit proc to another explicit proc are non-overtaking meaning
    they will arrive in the order they have been sent.

    \note This is an individual call.

    \param frompid (in)  : sending processors' pid
    \param topid (in)    : target processors' pid
    \param data (in)     : ptr to data to be send
    \param dsize (in)    : size of data (no. of integers)
    \param tag (in)      : tag to be used with message
    \param request (out) : mpi request handle to be used for testing completion of the
                           communication. data may not be altered or destroyed before
                           communcation finalized! One can use DRT::Exporter::Wait for this.

    \note This is an individual call
    */
    void ISend(const int frompid, const int topid, const int* data, const int dsize, const int tag,
        MPI_Request& request);

    /*!
    \brief Send data from one processor to another (nonblocking)

    The method will send an array of doubles in a nonblocking way meaning that this method
    will return immediately on the calling processor - even if communcation has not
    finished yet. The double array must not be altered or destroyed as long as the
    communication might still be in progress. This can be tested for using DRT::Exporter::Wait
    and the request handle returned.
    The receiving processor should call DRT::Exporter::ReceiveAny to receive the message.
    Note that messages from one explicit proc to another explicit proc are non-overtaking meaning
    they will arrive in the order they have been sent.

    \note This is an individual call.

    \param frompid (in)  : sending processors' pid
    \param topid (in)    : target processors' pid
    \param data (in)     : ptr to data to be send
    \param dsize (in)    : size of data (no. of doubles)
    \param tag (in)      : tag to be used with message
    \param request (out) : mpi request handle to be used for testing completion of the
                           communication. data may not be altered or destroyed before
                           communcation finalized! One can use DRT::Exporter::Wait for this.

    \note This is an individual call
    */
    void ISend(const int frompid, const int topid, const double* data, const int dsize,
        const int tag, MPI_Request& request);

    /*!
    \brief Receive anything joker (blocking)

    This method receives an MPI_CHAR string message from any source proc with
    any message tag of any length. It simply takes the first message that's
    coming in no matter from which sender of with which tag.
    recvbuff is resized to fit received message.
    the method is blocking for the calling (receiving) proc but for none of the
    other processors.
    It is used together with ISend and Wait to do nonblocking chaotic
    point to point communication.

    \note This is an individual call.

    \warning There is absolutely no guarantee about the order messages are
             received with this method except for one: Messages from the SAME sender
             to the SAME receiver will not overtake each other (which is not a really strong
             statement).

    \param source (output): source the message came from
    \param tag (output): message tag of message received
    \param recvbuff (output): buffer containing received data
    \param length (output): length of message upon receive
    */
    void ReceiveAny(int& source, int& tag, std::vector<char>& recvbuff, int& length);
    void Receive(const int source, const int tag, std::vector<char>& recvbuff, int& length);

    /*!
    \brief Receive anything joker (blocking)

    This method receives an MPI_INT message from any source proc with
    any message tag of any length. It simply takes the first message that's
    coming in no matter from which sender of with which tag.
    recvbuff is resized to fit received message.
    the method is blocking for the calling (receiving) proc but for none of the
    other processors.
    It is used together with ISend and Wait to do nonblocking chaotic
    point to point communication.

    \note This is an individual call.

    \warning There is absolutely no guarantee about the order messages are
             received with this method except for one: Messages from the SAME sender
             to the SAME receiver will not overtake each other (which is not a really strong
             statement).

    \param source (output): source the message came from
    \param tag (output): message tag of message received
    \param recvbuff (output): buffer containing received data
    \param length (output): length of message upon receive
    */
    void ReceiveAny(int& source, int& tag, std::vector<int>& recvbuff, int& length);
    void Receive(const int source, const int tag, std::vector<int>& recvbuff, int& length);

    /*!
    \brief Receive anything joker (blocking)

    This method receives an MPI_DOUBLE message from any source proc with
    any message tag of any length. It simply takes the first message that's
    coming in no matter from which sender of with which tag.
    recvbuff is resized to fit received message.
    the method is blocking for the calling (receiving) proc but for none of the
    other processors.
    It is used together with ISend and Wait to do nonblocking chaotic
    point to point communication.

    \note This is an individual call.

    \warning There is absolutely no guarantee about the order messages are
             received with this method except for one: Messages from the SAME sender
             to the SAME receiver will not overtake each other (which is not a really strong
             statement).

    \param source (output): source the message came from
    \param tag (output): message tag of message received
    \param recvbuff (output): buffer containing received data
    \param length (output): length of message upon receive
    */
    void ReceiveAny(int& source, int& tag, std::vector<double>& recvbuff, int& length);
    void Receive(const int source, const int tag, std::vector<double>& recvbuff, int& length);

    /*!
    \brief wait for nonblocking send to finish

    The method is used together with Isend and ReceiveAny to guarantee finalization
    of a communication. It is an individual call done by the sending processor to guarantee
    that message was taken from the sendbuffer before destroying the sendbuffer.
    This method is blocking and will return one communication associated with request has
    left the sender.

    \param request (in): mpi request handle

    */
    void Wait(MPI_Request& request)
    {
      MPI_Status status;
      MPI_Wait(&request, &status);
      return;
    }


    /*!
    \brief performs an allreduce operation on all processors
         and sends the result to all processors

    \param sendbuff (input): buffer containing data that has to be sent
    \param recvbuff (output): buffer containing received data
    \param mpi_op   (input): MPI operation
    */
    void Allreduce(std::vector<int>& sendbuff, std::vector<int>& recvbuff, MPI_Op mpi_op);

    //@}

   private:
    /*!
    \brief Do initialization of the exporter
    */
    void ConstructExporter();

    /*!
    \brief Get PID
    */
    inline int MyPID() const { return myrank_; }
    /*!
    \brief Get no. of processors
    */
    inline int NumProc() const { return numproc_; }

    /*!
    \brief Get sendplan_
    */
    inline std::vector<std::set<int>>& SendPlan() { return sendplan_; }

#if 0
  /*!
  \brief Get recvplan_
  */
  inline vector<std::set<int> > & RecvPlan() { return recvplan_; }
#endif

    /*!
    \brief generic export algorithm that delegates the specific pack/unpack to a helper
     */
    void GenericExport(ExporterHelper& helper);

   private:
    //! dummy map in case of empty exporter
    Epetra_Map dummymap_;
    //! source layout
    const Epetra_Map& frommap_;
    //! target map
    const Epetra_Map& tomap_;
    //! communicator
    const Epetra_Comm& comm_;
    //! PID
    int myrank_;
    //! no. of processors
    int numproc_;
    //! sending information
    std::vector<std::set<int>> sendplan_;
#if 0
  //! receiving information
  std::vector<std::set<int> >   recvplan_;
#endif

    /// Internal helper class for Exporter that encapsulates packing and unpacking
    /*!
      The communication algorithm we use to export a map of objects is the same
      independent of the actual type of objects we have. However, different
      object types require different packing and unpacking routines. So we put
      the type specific stuff in a helper class and get away with one clean
      communication algorithm. Nice.
     */
    class ExporterHelper
    {
     public:
      /// have a virtual destructor
      virtual ~ExporterHelper() {}

      /// validations performed before the communication
      virtual void PreExportTest(Exporter* exporter) = 0;

      /// Pack one object
      /*!
        Get the object by gid, pack it and append it to the sendblock. We only
        pack it if we know about it.
       */
      virtual bool PackObject(int gid, DRT::PackBuffer& sendblock) = 0;

      /// Unpack one object
      /*!
        After receiving we know the gid of the object and have its packed data
        at position index in recvblock. index must be incremented by the objects
        size.
       */
      virtual void UnpackObject(
          int gid, std::vector<char>::size_type& index, const std::vector<char>& recvblock) = 0;

      /// after communication remove all objects that are not in the target map
      virtual void PostExportCleanup(Exporter* exporter) = 0;
    };


    /// Concrete helper class that handles Teuchos::RCPs to ParObjects
    template <typename T>
    class ParObjectExporterHelper : public ExporterHelper
    {
     public:
      explicit ParObjectExporterHelper(std::map<int, Teuchos::RCP<T>>& parobjects)
          : parobjects_(parobjects)
      {
      }

      virtual void PreExportTest(Exporter* exporter)
      {
        // test whether type T implements ParObject
        typename std::map<int, Teuchos::RCP<T>>::iterator curr = parobjects_.begin();
        if (curr != parobjects_.end())
        {
          T* ptr = curr->second.get();
          ParObject* tester = dynamic_cast<ParObject*>(ptr);
          if (!tester)
            dserror(
                "typename T in template does not implement class ParObject (dynamic_cast failed)");
        }
      }

      virtual bool PackObject(int gid, DRT::PackBuffer& sendblock)
      {
        typename std::map<int, Teuchos::RCP<T>>::const_iterator curr = parobjects_.find(gid);
        if (curr != parobjects_.end())
        {
          ParObject::AddtoPack(sendblock, *curr->second);
          return true;
        }
        return false;
      }

      virtual void UnpackObject(
          int gid, std::vector<char>::size_type& index, const std::vector<char>& recvblock)
      {
        std::vector<char> data;
        ParObject::ExtractfromPack(index, recvblock, data);

        DRT::ParObject* o = DRT::UTILS::Factory(data);
        T* ptr = dynamic_cast<T*>(o);
        if (!ptr)
          dserror("typename T in template does not implement ParObject (dynamic_cast failed)");
        Teuchos::RCP<T> refptr = Teuchos::rcp(ptr);
        // add object to my map
        parobjects_[gid] = refptr;
      }

      virtual void PostExportCleanup(Exporter* exporter)
      {
        // loop map and kick out everything that's not in TargetMap()
        std::map<int, Teuchos::RCP<T>> newmap;
        typename std::map<int, Teuchos::RCP<T>>::const_iterator fool;
        for (fool = parobjects_.begin(); fool != parobjects_.end(); ++fool)
          if (exporter->TargetMap().MyGID(fool->first)) newmap[fool->first] = fool->second;
        swap(newmap, parobjects_);
      }

     private:
      std::map<int, Teuchos::RCP<T>>& parobjects_;
    };


    /// Concrete helper class that handles Teuchos::RCPs to any object
    /**!
       The objects considered here must have a default constructor and must be
       supported by ParObject::AddtoPack and ParObject::ExtractfromPack
       functions.

       Ideally one would manage ParObject with this helper as well, however, the
       cast to the concrete ParObject type prevents that.
     */
    template <typename T>
    class AnyObjectExporterHelper : public ExporterHelper
    {
     public:
      explicit AnyObjectExporterHelper(std::map<int, Teuchos::RCP<T>>& objects) : objects_(objects)
      {
      }

      virtual void PreExportTest(Exporter* exporter) {}

      virtual bool PackObject(int gid, DRT::PackBuffer& sendblock)
      {
        typename std::map<int, Teuchos::RCP<T>>::const_iterator curr = objects_.find(gid);
        if (curr != objects_.end())
        {
          ParObject::AddtoPack(sendblock, *curr->second);
          return true;
        }
        return false;
      }

      virtual void UnpackObject(
          int gid, std::vector<char>::size_type& index, const std::vector<char>& recvblock)
      {
        Teuchos::RCP<T> obj = Teuchos::rcp(new T);
        ParObject::ExtractfromPack(index, recvblock, *obj);

        // add object to my map
        objects_[gid] = obj;
      }

      virtual void PostExportCleanup(Exporter* exporter)
      {
        // loop map and kick out everything that's not in TargetMap()
        std::map<int, Teuchos::RCP<T>> newmap;
        typename std::map<int, Teuchos::RCP<T>>::const_iterator fool;
        for (fool = objects_.begin(); fool != objects_.end(); ++fool)
          if (exporter->TargetMap().MyGID(fool->first)) newmap[fool->first] = fool->second;
        swap(newmap, objects_);
      }

     private:
      std::map<int, Teuchos::RCP<T>>& objects_;
    };



    /// Concrete helper class that handles plain old data (POD) objects
    template <typename T>
    class PODExporterHelper : public ExporterHelper
    {
     public:
      explicit PODExporterHelper(std::map<int, T>& objects) : objects_(objects) {}

      virtual void PreExportTest(Exporter* exporter)
      {
        // Nothing to do. We do not check for T to be POD.
      }

      virtual bool PackObject(int gid, DRT::PackBuffer& sendblock)
      {
        typename std::map<int, T>::const_iterator curr = objects_.find(gid);
        if (curr != objects_.end())
        {
          ParObject::AddtoPack(sendblock, curr->second);
          return true;
        }
        return false;
      }

      virtual void UnpackObject(
          int gid, std::vector<char>::size_type& index, const std::vector<char>& recvblock)
      {
        memcpy(&objects_[gid], &recvblock[index], sizeof(T));
        index += sizeof(T);
      }

      virtual void PostExportCleanup(Exporter* exporter)
      {
        // loop map and kick out everything that's not in TargetMap()
        std::map<int, T> newmap;
        typename std::map<int, T>::const_iterator fool;
        for (fool = objects_.begin(); fool != objects_.end(); ++fool)
          if (exporter->TargetMap().MyGID(fool->first)) newmap[fool->first] = fool->second;
        swap(newmap, objects_);
      }

     private:
      std::map<int, T>& objects_;
    };


    /// Concrete helper class that handles vectors of plain old data (POD) objects
    template <typename T>
    class PODVectorExporterHelper : public ExporterHelper
    {
     public:
      explicit PODVectorExporterHelper(std::map<int, std::vector<T>>& objects) : objects_(objects)
      {
      }

      virtual void PreExportTest(Exporter* exporter)
      {
        // Nothing to do. We do not check for T to be POD.
      }

      virtual bool PackObject(int gid, DRT::PackBuffer& sendblock)
      {
        typename std::map<int, std::vector<T>>::const_iterator curr = objects_.find(gid);
        if (curr != objects_.end())
        {
          ParObject::AddtoPack(sendblock, curr->second);
          return true;
        }
        return false;
      }

      virtual void UnpackObject(
          int gid, std::vector<char>::size_type& index, const std::vector<char>& recvblock)
      {
        ParObject::ExtractfromPack(index, recvblock, objects_[gid]);
      }

      virtual void PostExportCleanup(Exporter* exporter)
      {
        // loop map and kick out everything that's not in TargetMap()
        std::map<int, std::vector<T>> newmap;
        typename std::map<int, std::vector<T>>::iterator fool;
        for (fool = objects_.begin(); fool != objects_.end(); ++fool)
          if (exporter->TargetMap().MyGID(fool->first)) swap(newmap[fool->first], fool->second);
        swap(newmap, objects_);
      }

     private:
      std::map<int, std::vector<T>>& objects_;
    };

    /// Concrete helper class that handles sets of plain old data (POD) objects
    template <typename T>
    class PODSetExporterHelper : public ExporterHelper
    {
     public:
      explicit PODSetExporterHelper(std::map<int, std::set<T>>& objects) : objects_(objects) {}

      virtual void PreExportTest(Exporter* exporter)
      {
        // Nothing to do. We do not check for T to be POD.
      }

      virtual bool PackObject(int gid, DRT::PackBuffer& sendblock)
      {
        typename std::map<int, std::set<T>>::const_iterator curr = objects_.find(gid);
        if (curr != objects_.end())
        {
          ParObject::AddtoPack(sendblock, curr->second);
          return true;
        }
        return false;
      }

      virtual void UnpackObject(
          int gid, std::vector<char>::size_type& index, const std::vector<char>& recvblock)
      {
        ParObject::ExtractfromPack(index, recvblock, objects_[gid]);
      }

      virtual void PostExportCleanup(Exporter* exporter)
      {
        // loop map and kick out everything that's not in TargetMap()
        std::map<int, std::set<T>> newmap;
        typename std::map<int, std::set<T>>::iterator fool;
        for (fool = objects_.begin(); fool != objects_.end(); ++fool)
          if (exporter->TargetMap().MyGID(fool->first)) swap(newmap[fool->first], fool->second);
        swap(newmap, objects_);
      }

     private:
      std::map<int, std::set<T>>& objects_;
    };

    /// Concrete helper class that handles maps of plain old data (POD) objects
    ///                                                          nagler 07/2012
    template <typename T, typename U>
    class PODMapExporterHelper : public ExporterHelper
    {
     public:
      explicit PODMapExporterHelper(std::map<int, std::map<T, U>>& objects) : objects_(objects) {}

      virtual void PreExportTest(Exporter* exporter)
      {
        // Nothing to do. We do not check for T to be POD.
      }

      virtual bool PackObject(int gid, DRT::PackBuffer& sendblock)
      {
        typename std::map<int, std::map<T, U>>::const_iterator curr = objects_.find(gid);
        if (curr != objects_.end())
        {
          ParObject::AddtoPack(sendblock, curr->second);
          return true;
        }
        return false;
      }

      virtual void UnpackObject(
          int gid, std::vector<char>::size_type& index, const std::vector<char>& recvblock)
      {
        ParObject::ExtractfromPack(index, recvblock, objects_[gid]);
      }

      virtual void PostExportCleanup(Exporter* exporter)
      {
        // loop map and kick out everything that's not in TargetMap()
        std::map<int, std::map<T, U>> newmap;
        typename std::map<int, std::map<T, U>>::iterator fool;
        for (fool = objects_.begin(); fool != objects_.end(); ++fool)
          if (exporter->TargetMap().MyGID(fool->first)) swap(newmap[fool->first], fool->second);
        swap(newmap, objects_);
      }

     private:
      std::map<int, std::map<T, U>>& objects_;
    };

  };  // class Exporter
}  // namespace DRT

/*----------------------------------------------------------------------*
 |  communicate objects (public)                             mwgee 11/06|
 *----------------------------------------------------------------------*/
template <typename T>
void DRT::Exporter::Export(std::map<int, Teuchos::RCP<T>>& parobjects)
{
  ParObjectExporterHelper<T> helper(parobjects);
  GenericExport(helper);
}

/*----------------------------------------------------------------------*
 |  communicate objects (public)                             mwgee 11/06|
 *----------------------------------------------------------------------*/
template <typename T>
void DRT::Exporter::Export(std::map<int, std::vector<T>>& data)
{
  PODVectorExporterHelper<T> helper(data);
  GenericExport(helper);
}

/*----------------------------------------------------------------------*
 |  communicate objects (public)                                   07/11|
 *----------------------------------------------------------------------*/
template <typename T>
void DRT::Exporter::Export(std::map<int, std::set<T>>& data)
{
  PODSetExporterHelper<T> helper(data);
  GenericExport(helper);
}

/*----------------------------------------------------------------------*
 |  communicate objects (public)                           nagler  07/12|
 *----------------------------------------------------------------------*/
template <typename T, typename U>
void DRT::Exporter::Export(std::map<int, std::map<T, U>>& data)
{
  PODMapExporterHelper<T, U> helper(data);
  GenericExport(helper);
}

#endif  // #ifndef EXPORTER_H
