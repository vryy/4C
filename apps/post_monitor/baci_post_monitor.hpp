/*----------------------------------------------------------------------*/
/*! \file

\brief monitoring filter for one data


\level 2

*/
/*----------------------------------------------------------------------*/

/*!
\addtogroup Monitoring
*//*! @{ (documentation module open)*/

#ifndef FOUR_C_POST_MONITOR_HPP
#define FOUR_C_POST_MONITOR_HPP

#include "baci_config.hpp"

#include "baci_lib_discret.hpp"
#include "baci_post_common.hpp"

#include <Teuchos_CommandLineProcessor.hpp>

#include <string>

BACI_NAMESPACE_OPEN

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
  MonWriter(PostProblem& problem, std::string& infieldtype, int node);

  //! destructor
  virtual ~MonWriter() = default;
  //! write something
  virtual void WriteMonFile(PostProblem& problem, std::string& infieldtype, int node);

  //! write something : stress a point
  void WriteMonStressFile(
      PostProblem& problem, std::string& infieldtype, std::string stresstype, int node);

  //! write something : strain a point
  void WriteMonStrainFile(
      PostProblem& problem, std::string& infieldtype, std::string straintype, int node);

  //! write something : strain a point
  void WriteMonPlStrainFile(
      PostProblem& problem, std::string& infieldtype, std::string straintype, int node);

  //! write something : heatflux a point
  void WriteMonHeatfluxFile(
      PostProblem& problem, std::string& infieldtype, std::string heatfluxtype, int node);

  //! write something : temperature gradient a point
  void WriteMonTempgradFile(
      PostProblem& problem, std::string& infieldtype, std::string tempgradtype, int node);

 protected:
  virtual PostField* GetFieldPtr(PostProblem& problem) = 0;

  virtual void CheckInfieldType(std::string& infieldtype) = 0;

  virtual void FieldError(int node) = 0;

  virtual void WriteHeader(std::ofstream& outfile) = 0;

  virtual void WriteTableHead(std::ofstream& outfile, int dim) = 0;

  virtual void WriteResult(
      std::ofstream& outfile, PostResult& result, std::vector<int>& gdof, int dim) = 0;

  void WriteMonStrFile(const std::string& filename, PostProblem& problem, std::string& infieldtype,
      const std::string strname, const std::string strtype, std::vector<std::string> groupnames,
      int node);

  virtual void WriteStrTableHead(
      std::ofstream& outfile, const std::string strname, const std::string strtype, const int dim)
  {
    dserror("Not impl.");
  }

  virtual void WriteStrResults(std::ofstream& outfile, PostProblem& problem, PostResult& result,
      std::vector<int>& gdof, int dim, std::string strtype, std::string groupname, const int node)
  {
    dserror("Not impl.");
  }

  void WriteMonThrFile(const std::string& filename, PostProblem& problem, std::string& infieldtype,
      const std::string thrname, const std::string thrtype, std::vector<std::string> groupnames,
      int node);

  virtual void WriteThrTableHead(
      std::ofstream& outfile, const std::string thrname, const std::string thrtype, const int dim)
  {
    dserror("Not impl.");
  }

  virtual void WriteThrResults(std::ofstream& outfile, PostProblem& problem, PostResult& result,
      std::vector<int>& gdof, int dim, std::string thrtype, std::string groupname, const int node)
  {
    dserror("Not impl.");
  }

  const int myrank_;  //! local processor id
  bool nodeowner_;    //! only true if proc owns the node

 private:
  // undesired copy constructor
  MonWriter(const MonWriter& old);
  // undesired = operator
  MonWriter& operator=(const MonWriter& old);

};  // end of class MonWriter



/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
class FieldMonWriter : public MonWriter
{
 public:
  //! constructor
  FieldMonWriter(PostProblem& problem, std::string& infieldtype, int node)
      : MonWriter(problem, infieldtype, node)
  {
  }

 protected:
  PostField* GetFieldPtr(PostProblem& problem) override;

 private:
};  // end of class FieldMonWriter



/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
class FluidMonWriter : public FieldMonWriter
{
 public:
  //! constructor
  FluidMonWriter(PostProblem& problem, std::string& infieldtype, int node)
      : FieldMonWriter(problem, infieldtype, node)
  {
  }

 protected:
  void CheckInfieldType(std::string& infieldtype) override;

  void FieldError(int node) override;

  void WriteHeader(std::ofstream& outfile) override;

  void WriteTableHead(std::ofstream& outfile, int dim) override;

  void WriteResult(
      std::ofstream& outfile, PostResult& result, std::vector<int>& gdof, int dim) override;

 private:
};  // end of class FluidMonWriter


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
class RedAirwayMonWriter : public FieldMonWriter
{
 public:
  //! constructor
  RedAirwayMonWriter(PostProblem& problem, std::string& infieldtype, int node)
      : FieldMonWriter(problem, infieldtype, node)
  {
  }

 protected:
  void CheckInfieldType(std::string& infieldtype) override;

  void FieldError(int node) override;

  void WriteHeader(std::ofstream& outfile) override;

  void WriteTableHead(std::ofstream& outfile, int dim) override;

  void WriteResult(
      std::ofstream& outfile, PostResult& result, std::vector<int>& gdof, int dim) override;

 private:
};  // end of class RedAirwayMonWriter


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
class StructMonWriter : public FieldMonWriter
{
 public:
  //! constructor
  StructMonWriter(PostProblem& problem, std::string& infieldtype, int node)
      : FieldMonWriter(problem, infieldtype, node)
  {
  }

 protected:
  void CheckInfieldType(std::string& infieldtype) override;

  void FieldError(int node) override;

  void WriteHeader(std::ofstream& outfile) override;

  void WriteTableHead(std::ofstream& outfile, int dim) override;

  void WriteResult(
      std::ofstream& outfile, PostResult& result, std::vector<int>& gdof, int dim) override;

  void WriteStrTableHead(std::ofstream& outfile, const std::string strname,
      const std::string strtype, const int dim) override;

  void WriteStrResults(std::ofstream& outfile, PostProblem& problem, PostResult& result,
      std::vector<int>& gdof, int dim, std::string strtype, std::string groupname,
      const int node) override;

  void WriteStrResult(std::ofstream& file, PostField*& field, PostResult& result,
      const std::string groupname, const std::string name, const int numdf, const int node) const;

 private:
};  // end of class StructMonWriter



/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
class AleMonWriter : public FieldMonWriter
{
 public:
  //! constructor
  AleMonWriter(PostProblem& problem, std::string& infieldtype, int node)
      : FieldMonWriter(problem, infieldtype, node)
  {
  }

 protected:
  void CheckInfieldType(std::string& infieldtype) override;

  void FieldError(int node) override;

  void WriteHeader(std::ofstream& outfile) override;

  void WriteTableHead(std::ofstream& outfile, int dim) override;

  void WriteResult(
      std::ofstream& outfile, PostResult& result, std::vector<int>& gdof, int dim) override;

 private:
};  // end of class AleMonWriter


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
class ScatraMonWriter : public FieldMonWriter
{
 public:
  //! constructor
  ScatraMonWriter(PostProblem& problem, std::string& infieldtype, int node)
      : FieldMonWriter(problem, infieldtype, node)
  {
  }

 protected:
  void CheckInfieldType(std::string& infieldtype) override;

  void FieldError(int node) override;

  void WriteHeader(std::ofstream& outfile) override;

  void WriteTableHead(std::ofstream& outfile, int dim) override;

  PostField* GetFieldPtr(PostProblem& problem) override;

  void WriteResult(
      std::ofstream& outfile, PostResult& result, std::vector<int>& gdof, int dim) override;

 private:
};  // end of class ScatraMonWriter


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
class ThermoMonWriter : public FieldMonWriter
{
 public:
  //! constructor
  ThermoMonWriter(PostProblem& problem, std::string& infieldtype, int node)
      : FieldMonWriter(problem, infieldtype, node)
  {
  }

 protected:
  void CheckInfieldType(std::string& infieldtype) override;

  void FieldError(int node) override;

  void WriteHeader(std::ofstream& outfile) override;

  void WriteTableHead(std::ofstream& outfile, int dim) override;

  void WriteResult(
      std::ofstream& outfile, PostResult& result, std::vector<int>& gdof, int dim) override;

  void WriteThrTableHead(std::ofstream& outfile, const std::string thrname,
      const std::string thrtype, const int dim) override;

  void WriteThrResults(std::ofstream& outfile, PostProblem& problem, PostResult& result,
      std::vector<int>& gdof, int dim, std::string thrtype, std::string groupname,
      const int node) override;

  void WriteThrResult(std::ofstream& file, PostField*& field, PostResult& result,
      const std::string groupname, const std::string name, const int dim, const int node) const;

 private:
};  // end of class ThermoMonWriter



/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
class FsiFluidMonWriter : public FluidMonWriter
{
 public:
  //! constructor
  FsiFluidMonWriter(PostProblem& problem, std::string& infieldtype, int node)
      : FluidMonWriter(problem, infieldtype, node)
  {
  }

 protected:
  void CheckInfieldType(std::string& infieldtype) override{};

  PostField* GetFieldPtr(PostProblem& problem) override;

  void WriteHeader(std::ofstream& outfile) override;

  void WriteTableHead(std::ofstream& outfile, int dim) override;

  void WriteResult(
      std::ofstream& outfile, PostResult& result, std::vector<int>& gdof, int dim) override;

 private:
};  // end of class FsiFluidMonWriter



/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
class FsiStructMonWriter : public StructMonWriter
{
 public:
  //! constructor
  FsiStructMonWriter(PostProblem& problem, std::string& infieldtype, int node)
      : StructMonWriter(problem, infieldtype, node)
  {
  }

 protected:
  void CheckInfieldType(std::string& infieldtype) override{};

  PostField* GetFieldPtr(PostProblem& problem) override;

  void WriteHeader(std::ofstream& outfile) override;

  void WriteTableHead(std::ofstream& outfile, int dim) override;

  void WriteResult(
      std::ofstream& outfile, PostResult& result, std::vector<int>& gdof, int dim) override;

 private:
};  // end of class FsiStructMonWriter



/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
class FsiAleMonWriter : public AleMonWriter
{
 public:
  //! constructor
  FsiAleMonWriter(PostProblem& problem, std::string& infieldtype, int node)
      : AleMonWriter(problem, infieldtype, node)
  {
  }

 protected:
  void CheckInfieldType(std::string& infieldtype) override{};

  PostField* GetFieldPtr(PostProblem& problem) override;

  void WriteHeader(std::ofstream& outfile) override;

 private:
};  // end of class FsiAleMonWriter



/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
class TsiStructMonWriter : public StructMonWriter
{
 public:
  //! constructor
  TsiStructMonWriter(PostProblem& problem, std::string& infieldtype, int node)
      : StructMonWriter(problem, infieldtype, node)
  {
  }

 protected:
  void CheckInfieldType(std::string& infieldtype) override{};

  PostField* GetFieldPtr(PostProblem& problem) override;

  void WriteHeader(std::ofstream& outfile) override;

 private:
};  // end of class TsiStructMonWriter


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
class TsiThermoMonWriter : public ThermoMonWriter
{
 public:
  //! constructor
  TsiThermoMonWriter(PostProblem& problem, std::string& infieldtype, int node)
      : ThermoMonWriter(problem, infieldtype, node)
  {
  }

 protected:
  void CheckInfieldType(std::string& infieldtype) override{};

  PostField* GetFieldPtr(PostProblem& problem) override;

  void WriteHeader(std::ofstream& outfile) override;

 private:
};  // end of class TsiThermoMonWriter

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
/**
 * \brief monitoring for porofluidmultiphase-problems
 */
class PoroFluidMultiMonWriter : public FieldMonWriter
{
 public:
  //! constructor
  PoroFluidMultiMonWriter(PostProblem& problem, std::string& infieldtype, int node)
      : FieldMonWriter(problem, infieldtype, node), poro_dof_(node)
  {
    output_porosity_ = false;
    if (problem.optquantitytype() == "porosity") output_porosity_ = true;
  }

 protected:
  //! check if infieldtype is correct
  void CheckInfieldType(std::string& infieldtype) override;

  //! print out error
  void FieldError(int node) override;

  //! write header into output file
  void WriteHeader(std::ofstream& outfile) override;

  //! write table head into output file
  void WriteTableHead(std::ofstream& outfile, int dim) override;

  //! write result into output file
  void WriteResult(
      std::ofstream& outfile, PostResult& result, std::vector<int>& gdof, int dim) override;

 private:
  //! flag for output of porosity
  bool output_porosity_;
  //! dof for output of porosity ( = node number, since defined as nodal quantity)
  int poro_dof_;
};  // end of class FluidMonWriter

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
/**
 * \brief monitoring for porofluidmultiphase-quantities as part of a
 *        Poromultiphase(Scatra)-problem
 */
class PoroMultiElastScatraFluidMonWriter : public PoroFluidMultiMonWriter
{
 public:
  //! constructor
  PoroMultiElastScatraFluidMonWriter(PostProblem& problem, std::string& infieldtype, int node)
      : PoroFluidMultiMonWriter(problem, infieldtype, node)
  {
  }

 protected:
  //! get pointer to field
  PostField* GetFieldPtr(PostProblem& problem) override;

  //! write header into output file
  void WriteHeader(std::ofstream& outfile) override;

 private:
};  // end of class PoroMultiElastFluidMonWriter

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
/**
 * \brief monitoring for scatra-quantities as part of a
 *        Poromultiphase(Scatra)-problem
 */
class PoroMultiElastScatraScatraMonWriter : public ScatraMonWriter
{
 public:
  //! constructor
  PoroMultiElastScatraScatraMonWriter(PostProblem& problem, std::string& infieldtype, int node)
      : ScatraMonWriter(problem, infieldtype, node)
  {
  }

 protected:
  //! get pointer to field
  PostField* GetFieldPtr(PostProblem& problem) override;

  //! write header into output file
  void WriteHeader(std::ofstream& outfile) override;

 private:
};  // end of class PoroMultiElastScatraScatraMonWriter

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
/**
 * \brief monitoring for artery_scatra-quantities (1D discretization) as part of a
 *        Poromultiphase(Scatra)-problem
 */
class PoroMultiElastScatraArteryScatraMonWriter : public ScatraMonWriter
{
 public:
  //! constructor
  PoroMultiElastScatraArteryScatraMonWriter(
      PostProblem& problem, std::string& infieldtype, int node)
      : ScatraMonWriter(problem, infieldtype, node)
  {
  }

 protected:
  //! get pointer to field
  PostField* GetFieldPtr(PostProblem& problem) override;

  //! write header into output file
  void WriteHeader(std::ofstream& outfile) override;

 private:
};  // end of class PoroMultiElastScatraArteryScatraMonWriter


#endif  // POST_MONITOR_H
/*! @} (documentation module close)*/

BACI_NAMESPACE_CLOSE
