/*----------------------------------------------------------------------*/
/*! \file

\brief Managing of space- and/or time-dependent functions

\level 0

*/
/*----------------------------------------------------------------------*/

#ifndef DRT_FUNCTION_MANAGER_H
#define DRT_FUNCTION_MANAGER_H


#include <iostream>
#include <vector>
#include "drt_function.H"
#include "drt_dserror.H"
#include "drt_parser.H"

namespace DRT
{
  // class Condition;
  // class Problem;

  namespace INPUT
  {
    class Lines;
    class DatFileReader;

  }  // namespace INPUT

  namespace UTILS
  {
    // forward declaration
    class Function;

    /// class that manages all defined spatial functions
    class FunctionManager
    {
     public:
      Teuchos::RCP<DRT::INPUT::Lines> ValidFunctionLines();

      /// read the baci input file and set up all functions
      void ReadInput(DRT::INPUT::DatFileReader& reader);

      /// return a specified function
      DRT::UTILS::Function& Funct(int num);

     private:
      /// all functions defined. The functions may be vector valued!
      std::vector<Teuchos::RCP<Function>> functions_;
    };

    /// add valid function lines
    void AddValidFunctionFunctionLines(Teuchos::RCP<DRT::INPUT::Lines> lines);
  }  // namespace UTILS

}  // namespace DRT

void PrintFunctionDatHeader();

#endif