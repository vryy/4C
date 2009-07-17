/*----------------------------------------------------------------------*/
/*!
\file post_drt_monitor.cpp

\brief monitoring filter for one data

<pre>
Maintainer: Christiane Förster
            foerster@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/foerster
            089 - 289-15262
</pre>

*/
/*----------------------------------------------------------------------*/

/*!
\addtogroup Monitoring
*//*! @{ (documentation module open)*/
#ifdef CCADISCRET

#ifndef POST_DRT_MONITOR_H_
#define POST_DRT_MONITOR_H_

#include <string>
#include <Teuchos_CommandLineProcessor.hpp>

#include "../post_drt_common/post_drt_common.H"
#include "../drt_lib/drt_discret.H"


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
/*!
 * \brief pure virtual class to do monitoring
 */
class MonWriter
{
public:
  
  //! constructor
  MonWriter(
    PostProblem& problem,
    string& infieldtype,
    int node
    );

  //! destructor
  virtual ~MonWriter()
  {}

  //! write something
  virtual void WriteMonFile(
    PostProblem& problem,
    string& infieldtype,
    int node
    );

  //! write something : stress a point
  void WriteMonStressFile(
    PostProblem& problem,
    string& infieldtype,
    string stresstype,
    int node
    );

  //! write something : strain a point
  void WriteMonStrainFile(
    PostProblem& problem,
    string& infieldtype,
    string straintype,
    int node
    );

protected:

  virtual PostField* GetFieldPtr(PostProblem& problem) = 0;

  virtual void CheckInfieldType(string& infieldtype) = 0;

  virtual void FieldError(int node) = 0;

  virtual void WriteHeader(ofstream& outfile) = 0;

  virtual void WriteTableHead(ofstream& outfile, int dim) = 0;

  virtual void WriteResult(
    ofstream& outfile,
    PostResult& result,
    std::vector<int>& gdof,
    int dim
    ) = 0;

  void WriteMonStrFile(
    const string& filename,
    PostProblem& problem,
    string& infieldtype,
    const string strname,
    const string strtype,
    vector<string> groupnames,
    int node
    );

  virtual void WriteStrTableHead(
    ofstream& outfile,
    const string strname,
    const string strtype,
    const int dim
    )
  {
    dserror("Not impl.");
  }

  virtual void WriteStrResults(
    ofstream& outfile,
    PostProblem& problem,
    PostResult& result,
    std::vector<int>& gdof,
    int dim,
    string strtype,
    string groupname,
    const int node)
  {
    dserror("Not impl.");
  }

  const int myrank_; //! local processor id
  bool nodeowner_;   //! only true if proc owns the node

private:
  // undesired copy constructor
  MonWriter(const MonWriter& old);
  // undesired = operator
  MonWriter& operator= (const MonWriter& old);

}; // end of class MonWriter



/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
class FieldMonWriter : public MonWriter
{
public:
  //! constructor
  FieldMonWriter(
    PostProblem& problem,
    string& infieldtype,
    int node
    )
    : MonWriter(problem,infieldtype,node)
  {}

  //! destructor
  virtual ~FieldMonWriter()
  {}

protected:

  virtual PostField* GetFieldPtr(PostProblem& problem);

private:
}; // end of class FieldMonWriter



/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
class FluidMonWriter : public FieldMonWriter
{
public:
  //! constructor
  FluidMonWriter(
    PostProblem& problem,
    string& infieldtype,
    int node)
    : FieldMonWriter(problem,infieldtype,node)
  {}

  //! destructor
  virtual ~FluidMonWriter()
  {}

protected:

  virtual void CheckInfieldType(string& infieldtype);

  void FieldError(int node);

  virtual void WriteHeader(ofstream& outfile);

  void WriteTableHead(ofstream& outfile, int dim);

  virtual void WriteResult(
    ofstream& outfile,
    PostResult& result,
    std::vector<int>& gdof,
    int dim
    );

private:
}; // end of class FluidMonWriter



/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
class StructMonWriter : public FieldMonWriter
{
public:
  //! constructor
  StructMonWriter(
    PostProblem& problem,
    string& infieldtype,
    int node
    )
    : FieldMonWriter(problem,infieldtype,node)
  {}

  //! destructor
  virtual ~StructMonWriter()
  {}

protected:

  virtual void CheckInfieldType(string& infieldtype);

  void FieldError(int node);

  virtual void WriteHeader(ofstream& outfile);

  void WriteTableHead(ofstream& outfile, int dim);

  void WriteResult(
    ofstream& outfile,
    PostResult& result,
    std::vector<int>& gdof,
    int dim
    );

  virtual void WriteStrTableHead(
    ofstream& outfile,
    const string strname,
    const string strtype,
    const int dim
    );

  virtual void WriteStrResults(
    ofstream& outfile,
    PostProblem& problem,
    PostResult& result,
    std::vector<int>& gdof,
    int dim,
    string strtype,
    string groupname,
    const int node
    );

  void WriteStrResult(
    std::ofstream& file,
    PostField*& field,
    PostResult& result,
    const std::string groupname,
    const std::string name,
    const int numdf,
    const int node
    ) const;

private:
}; // end of class StructMonWriter



/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
class AleMonWriter : public FieldMonWriter
{
public:
  //! constructor
  AleMonWriter(
    PostProblem& problem,
    string& infieldtype,
    int node
    )
    : FieldMonWriter(problem,infieldtype,node)
  {}

  //! destructor
  virtual ~AleMonWriter()
  {}

protected:

  virtual void CheckInfieldType(string& infieldtype);

  void FieldError(int node);

  void WriteHeader(ofstream& outfile);

  void WriteTableHead(ofstream& outfile, int dim);

  void WriteResult(
    ofstream& outfile,
    PostResult& result,
    std::vector<int>& gdof,
    int dim
    );

private:
}; // end of class AleMonWriter



/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
class FsiFluidMonWriter : public FluidMonWriter
{
public:
  //! constructor
  FsiFluidMonWriter(
    PostProblem& problem,
    string& infieldtype,
    int node
    )
    : FluidMonWriter(problem,infieldtype,node)
  {}

  //! destructor
  virtual ~FsiFluidMonWriter()
  {}

protected:

  void CheckInfieldType(string& infieldtype){};

  PostField* GetFieldPtr(PostProblem& problem);

  void WriteHeader(ofstream& outfile);

  void WriteTableHead(ofstream& outfile, int dim);

  void WriteResult(
    ofstream& outfile,
    PostResult& result,
    std::vector<int>& gdof,
    int dim
    );

private:
}; // end of class FsiFluidMonWriter



/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
class FsiStructMonWriter : public StructMonWriter
{
public:
  //! constructor
  FsiStructMonWriter(
    PostProblem& problem,
    string& infieldtype,
    int node
    )
    : StructMonWriter(problem,infieldtype,node)
  {}

  //! destructor
  virtual ~FsiStructMonWriter()
  {}

protected:

  void CheckInfieldType(string& infieldtype){};

  PostField* GetFieldPtr(PostProblem& problem);

  void WriteHeader(ofstream& outfile);

private:
}; // end of class FsiStructMonWriter



/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
class FsiAleMonWriter : public AleMonWriter
{
public:
  //! constructor
  FsiAleMonWriter(
    PostProblem& problem,
    string& infieldtype,
    int node
    )
    : AleMonWriter(problem,infieldtype,node)
  {}

  //! destructor
  virtual ~FsiAleMonWriter()
  {}

protected:

  void CheckInfieldType(string& infieldtype){};

  PostField* GetFieldPtr(PostProblem& problem);

  void WriteHeader(ofstream& outfile);

private:
}; // end of class FsiAleMonWriter



#endif  /*POST_DRT_MONITOR_H_*/
#endif  /*CCADISCRET*/
/*! @} (documentation module close)*/
