/*----------------------------------------------------------------------------*/
/*! \file

\brief collect data in parallel

\level 1


*/
/*----------------------------------------------------------------------------*/

#ifndef SRC_DRT_LIB_DATA_COLLECTOR_H_
#define SRC_DRT_LIB_DATA_COLLECTOR_H_

#include "drt_exporter.H"

namespace DRT
{
  /** Collect data over multiple processors
   *
   *  \author hiermeier \date 11/17 */
  template <typename T>
  class DataCollector
  {
    DataCollector() = delete;
    ~DataCollector() = delete;
    DataCollector(const DataCollector<T>& collector) = delete;

   public:
    static void CollectData(const Epetra_Comm& comm, const T& my_data, T& collected_data)
    {
      DRT::PackBuffer pack_data;

      CountMyData(my_data, pack_data);
      PackMyData(my_data, pack_data);

      RoundRobinLoop(comm, pack_data, collected_data);
    }

   private:
    static void RoundRobinLoop(
        const Epetra_Comm& comm, DRT::PackBuffer& pack_data, T& collected_data)
    {
      // collect the information over all procs
      std::vector<char> mydata;

      // swap into std::vector
      std::swap(mydata, pack_data());

      std::vector<int> mysize(1, mydata.size());

      std::vector<int> receivedsize;
      std::vector<char> receiveddata;

      // create an exporter for point to point communication
      DRT::Exporter exporter(comm);
      const int numprocs = comm.NumProc();

      for (int p = 0; p < numprocs; ++p)
      {
        switch (p)
        {
          case 0:
          {
            std::swap(receivedsize, mysize);
            std::swap(receiveddata, mydata);
            break;
          }
          default:
          {
            SendToNextProc(p, exporter, mysize, mydata, receivedsize, receiveddata);

            break;
          }
        }

        // unpack received block
        UnpackReceivedBlock(receiveddata, collected_data);

        // the received data will be sent to the next proc
        std::swap(receivedsize, mysize);
        std::swap(receiveddata, mydata);

        // we need a new receive buffer
        receivedsize.clear();
        receiveddata.clear();
      }

      //      for (int p=0; p<numprocs; ++p)
      //      {
      //        if (comm.MyPID()==p)
      //        {
      //          std::cout << "--------------- Processor " << p << " ---------------------" <<
      //          std::endl; std::cout << "---------------  size = " << collected_data.size() <<
      //              " ----------------------" << std::endl;
      //          for ( auto& cd : collected_data )
      //          {
      //            std::cout << "first: " << cd.first << " | second: " << cd.second
      //                << "\n";
      //          }
      //          std::cout << "\n\n\n\n";
      //        }
      //        comm.Barrier();
      //      }
    }

    static void UnpackReceivedBlock(const std::vector<char>& receiveddata, T& collected_data)
    {
      std::vector<char>::size_type index = 0;
      int j = 0;
      while (index < receiveddata.size())
      {
        // the set gets cleared at the beginning of the ExtractfromPack routine!
        T rs;
        DRT::ParObject::ExtractfromPack(index, receiveddata, rs);
        collected_data.insert(rs.begin(), rs.end());
        ++j;
      }
      // sanity check
      if (index > receiveddata.size())
        dserror(
            "Something is messed up in the received data block! Expected "
            "size = %d <--> received size = %d",
            receiveddata.size(), index);
    }

    static void SendToNextProc(const int p, DRT::Exporter& exporter, const std::vector<int>& mysize,
        const std::vector<char>& mydata, std::vector<int>& receivedsize,
        std::vector<char>& receiveddata)
    {
      const Epetra_Comm& comm = exporter.Comm();

      const int numprocs = comm.NumProc();
      const int myrank = comm.MyPID();
      int tag = myrank;

      int frompid = myrank;
      const int topid = (myrank + 1) % numprocs;

      // send size
      MPI_Request sizerequest;
      exporter.ISend(frompid, topid, mysize.data(), mysize.size(), tag, sizerequest);

      // send data
      MPI_Request datarequest;
      exporter.ISend(frompid, topid, mydata.data(), mydata.size(), tag * 10, datarequest);

      // make sure that you do not think you received something if
      // you didn't
      if (not receiveddata.empty() or not receivedsize.empty())
        dserror("Received data objects are not empty!");

      // receive from predecessor
      frompid = (myrank + numprocs - 1) % numprocs;

      // receive size information
      int length = 0;
      exporter.ReceiveAny(frompid, tag, receivedsize, length);
      if (length != static_cast<int>(mysize.size()) or tag != frompid)
        dserror(
            "Size information got mixed up!\n"
            "Received length = %d, Expected length = %d \n"
            "Received tag    = %d, Expected tag    = %d",
            length, mysize.size(), tag, frompid);

      exporter.Wait(sizerequest);

      // receive the gids
      exporter.ReceiveAny(frompid, tag, receiveddata, length);
      if (length != receivedsize[0] or tag != frompid * 10)
        dserror(
            "Data information got mixed up! \n"
            "Received length = %d, Expected length = %d \n"
            "Received tag    = %d, Expected tag    = %d",
            length, receivedsize[0], tag, frompid * 10);

      exporter.Wait(datarequest);

      return;
    }

    static void CountMyData(const T& my_data, DRT::PackBuffer& pack_data)
    {
      DRT::ParObject::AddtoPack(pack_data, my_data);
    }

    static void PackMyData(const T& my_data, DRT::PackBuffer& pack_data)
    {
      pack_data.StartPacking();
      DRT::ParObject::AddtoPack(pack_data, my_data);
    }
  };

  template <typename T>
  void CollectData(const Epetra_Comm& comm, const T& my_data, T& collected_data)
  {
    DataCollector<T>::CollectData(comm, my_data, collected_data);
  }

}  // namespace DRT

#endif /* SRC_DRT_LIB_DATA_COLLECTOR_H_ */
