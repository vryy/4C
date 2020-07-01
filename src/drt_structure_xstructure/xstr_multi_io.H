/*----------------------------------------------------------------------------*/
/** \file

\brief Input/output routines for the multi discretization wrapper


\level 3

*/
/*----------------------------------------------------------------------------*/

#ifndef SRC_DRT_STRUCTURE_XSTRUCTURE_XSTR_MULTI_IO_H_
#define SRC_DRT_STRUCTURE_XSTRUCTURE_XSTR_MULTI_IO_H_

#include "../drt_io/io.H"
#include "../drt_xfem/xfem_multi_field_mapextractor.H"

namespace DRT
{
  class DiscretizationInterface;
}  // namespace DRT
namespace XSTR
{
  class MultiDiscretizationWrapper;
  namespace IO
  {
    class DiscretizationReader;
    class DiscretizationWriter;

    /*-------------------------------------------------------------------------*/
    /** \brief create a temporal partial system multi-vector pointer
     *
     *  \param vt       (in): desired vector type
     *  \param num_vecs (in): number of vectors of the multi vector object
     *  \param dis      (in): discretization */
    Epetra_MultiVector* CreatePartialMultiVector(const enum ::IO::VectorType& vt,
        const int& num_vecs, const DRT::DiscretizationInterface& dis);

    /** \class MultiDiscretizationReader
     *
     *  \brief Wrapper class to manage the reading process of more than one
     *  discretization at once.
     *
     *  \author hiermeier \date 10/16 */
    class MultiDiscretizationReader : public ::IO::DiscretizationReader
    {
      typedef std::map<enum XFEM::FieldName, Teuchos::RCP<XSTR::IO::DiscretizationReader>>
          XReaderMap;

     public:
      /** \brief constructor
       *
       *   reader for a given discretization to read a particular time step */
      MultiDiscretizationReader(const Teuchos::RCP<XSTR::MultiDiscretizationWrapper>& dis_wrapper,
          const Teuchos::RCP<::IO::InputControl>& input, int step);

      /** \brief constructor
       *
       *  construct reader for a given discretization to read a particular time step */
      MultiDiscretizationReader(
          const Teuchos::RCP<XSTR::MultiDiscretizationWrapper>& dis_wrapper, int step);

      /// destructor
      virtual ~MultiDiscretizationReader(){};

      /// read into given vector  [derived]
      virtual void ReadVector(Teuchos::RCP<Epetra_Vector> vec, std::string name);

      /// read into given vector  [derived]
      virtual void ReadMultiVector(Teuchos::RCP<Epetra_MultiVector> vec, std::string name);

      /// read into given std::map<int, Teuchos::RCP<Epetra_SerialDenseMatrix> >  [derived]
      virtual void ReadSerialDenseMatrix(
          Teuchos::RCP<std::map<int, Teuchos::RCP<Epetra_SerialDenseMatrix>>> mapdata,
          std::string name);

      /// read an integer value from the control file  [derived]
      virtual int ReadInt(std::string name);

      /// read a double value from the control file  [derived]
      virtual double ReadDouble(std::string name);

      /// read into the discretization given in the constructor  [derived]
      virtual void ReadMesh(int step);

      /// read nodes into the discretization given in the constructor  [derived]
      virtual void ReadNodesOnly(int step);

      /// Read the history data of elements and nodes from restart files  [derived]
      virtual void ReadHistoryData(int step);

      /// read a non discretisation based vector of doubles  [derived]
      virtual void ReadRedundantDoubleVector(
          Teuchos::RCP<std::vector<double>>& doublevec, const std::string name);

      /// read a non discretisation based vector of integers  [derived]
      virtual void ReadRedundantIntVector(
          Teuchos::RCP<std::vector<int>>& intvec, const std::string name);

      /// return number of procs which were used for restart output (read from control file)
      /// [derived]
      virtual int GetNumOutputProc(int step);

     protected:
      /// access the communicator (read-only)
      virtual const Epetra_Comm& Comm() const;

     private:
      /// pointer to the multi discretization wrapper
      Teuchos::RCP<XSTR::MultiDiscretizationWrapper> dis_wrapper_;

      XReaderMap readers_;
    };  // class MultiDiscretizationReader

    /*-------------------------------------------------------------------------*/
    /** \class DiscretizationReader
     *
     *  \brief Derived DiscretizationReader which implements additional functionality
     *  necessary for the MultiDiscretizationReader wrapper class.
     *
     *  \author hiermeier \date 10/16 */
    class DiscretizationReader : public ::IO::DiscretizationReader
    {
     public:
      /// construct reader for a given discretization to read a particular time step [derived]
      DiscretizationReader(Teuchos::RCP<DRT::DiscretizationInterface> dis,
          Teuchos::RCP<::IO::InputControl> input, int step);

      /// construct reader for a given discretization to read a particular time step [derived]
      DiscretizationReader(Teuchos::RCP<DRT::DiscretizationInterface> dis, int step);

      /// destructor
      virtual ~DiscretizationReader(){};

      /** \brief get the vector type of a given vector name
       *  \param name (in): name of the vector */
      enum ::IO::VectorType VectorType(const std::string& name);
    };  // class DiscretizationReader

    /*-------------------------------------------------------------------------*/
    /** \class MultiDiscretizationWriter
     *
     *  \brief Wrapper class to manage the writing process of more than one
     *  discretization at once.
     *
     *  \author hiermeier \date 10/16 */
    class MultiDiscretizationWriter : public ::IO::DiscretizationWriter
    {
      typedef std::map<enum XFEM::FieldName, Teuchos::RCP<XSTR::IO::DiscretizationWriter>>
          XWriterMap;

     public:
      /// constructor
      MultiDiscretizationWriter(const Teuchos::RCP<XSTR::MultiDiscretizationWrapper>& dis_wrapper);

      /// destructor
      virtual ~MultiDiscretizationWriter(){};

      /// @name Output methods
      /// @{
      /// write result header to control file [derived]
      virtual void NewStep(const int step, const double time);

      /// write a result double to control file [derived]
      virtual void WriteDouble(const std::string name, const double value);

      /// write a result integer to control file [derived]
      virtual void WriteInt(const std::string name, const int value);

      /// write a result vector [derived]
      virtual void WriteVector(const std::string name, Teuchos::RCP<const Epetra_MultiVector> vec,
          ::IO::VectorType vt = ::IO::dofvector);

      /// write a result vector [derived]
      virtual void WriteVector(const std::string name, const std::vector<char>& vec,
          const Epetra_Map& elemap, ::IO::VectorType vt = ::IO::dofvector);

      /// write new mesh and result file next time it is possible
      virtual void CreateNewResultAndMeshFile();

      /// write new "field" group to control file including node and element chunks
      virtual void WriteMesh(const int step, const double time);

      /** for MLMC purposes do not write new meshfile but write name
       *  of base mesh file to controlfile [derived] */
      virtual void WriteMesh(const int step, const double time, std::string name_base_file);

      /// write element data to file [derived]
      virtual void WriteElementData(bool writeowner);

      /// write node data to file [derived]
      virtual void WriteNodeData(bool writeowner);

      /// write a non discretisation based vector of doubles [derived]
      virtual void WriteRedundantDoubleVector(
          const std::string name, Teuchos::RCP<std::vector<double>> doublevec);

      /// write a non discretisation based vector of integers [derived]
      virtual void WriteRedundantIntVector(
          const std::string name, Teuchos::RCP<std::vector<int>> intvec);

      /// overwrite result files for the inverse analysis [derived]
      virtual void OverwriteResultFile();

      /// creating new result files for the inverse analysis [derived]
      virtual void NewResultFile(int numb_run);

      /// creating new result files for the mlmc [derived]
      virtual void NewResultFile(std::string name_appendix, int numb_run);

      /// @}

      /// @name Data management
      /// @{

      /// clear all stored map data [derived]
      virtual void ClearMapCache();

      /// @}

      /// get output control [derived]
      virtual Teuchos::RCP<::IO::OutputControl> Output() const;

      /// set output control [derived]
      virtual void SetOutput(Teuchos::RCP<::IO::OutputControl> output);

     protected:
      /// access the Epetra_Comm object [derived]
      virtual const Epetra_Comm& Comm() const;

      /// pointer to the multi discretization wrapper
      Teuchos::RCP<XSTR::MultiDiscretizationWrapper> dis_wrapper_;

      XWriterMap writers_;

      int curr_step_;
    };  // MultiDiscretizationWriter

    /*-------------------------------------------------------------------------*/
    /** \class DiscretizationWriter
     *
     *  \brief Derived DiscretizationWriter which implements additional functionality
     *  necessary for the MultiDiscretizationWriter wrapper class.
     *
     *  \author hiermeier \date 11/16 */
    class DiscretizationWriter : public ::IO::DiscretizationWriter
    {
     public:
      /** Writer constructor always needs a discretization
       *  Use default control file from global problem. */
      DiscretizationWriter(Teuchos::RCP<DRT::DiscretizationInterface> dis);

      /// destructor
      virtual ~DiscretizationWriter(){};

      /** Set the output stream at the right position in the control file  */
      unsigned AugmentControlFile(int step, std::vector<char>& buffer);

      void AugmentControlFile(const char* buffer, unsigned bufferlength);

    };  // class DiscretizationReader

  }  // namespace IO
}  // namespace XSTR


#endif /* SRC_DRT_STRUCTURE_XSTRUCTURE_XSTR_MULTI_IO_H_ */
