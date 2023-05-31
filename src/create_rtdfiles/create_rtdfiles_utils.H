/*---------------------------------------------------------------------*/
/*! \file

\brief A collection of helper functions and classes for the readthedocs parser

\level 0

*/
/*---------------------------------------------------------------------*/

#ifndef CREATE_RTDFILES_UTILS_H
#define CREATE_RTDFILES_UTILS_H
#include <string>
#include <iostream>
#include <vector>
#include <map>
#include "lib_materialdefinition.H"
#include "lib_conditiondefinition.H"
#include "inpar_validparameters.H"
#include "contact_constitutivelaw_constitutivelaw_definition.H"
#include <Teuchos_RCP.hpp>

namespace DRT
{
  namespace RTD
  {
    /*!
     * /brief The Table class is used to create a table in a readthedocs document (written in
     * restructuredText).
     *
     * The table for readthedocs is a so-called list-table.
     * One may enter directives (like "align", "header-rows"), and also enter the relative widths of
     * columns, which is then converted into a directive as well.
     * Then, one may enter table content row by row.
     *
     * The construction parameter defines the number of columns in the table.
     * Of course, if you add a row with a different number of cells than predefined, an error is
     * thrown.

    */
    class Table
    {
     public:
      /*!
       * The table constructor doesn't do much besides defining the number of columns of the table.
       *
       * @param size Number of columns
       */
      Table(const unsigned& size);
      /*!
       * The content of the table is entered row by row.
       *
       * @param row This vector of strings fills the row. Of course, the size of this vector must be
       * equal to the size of the table (construction parameter). Otherwise an error is thrown.
       */
      void AddRow(const std::vector<std::string>& row);
      /*!
       * The relative widths of the table columns are set here. The numbers should be given as
       * characters after which a new line is created, that is, each cell gets an enforced line
       * break at a space closest to the i'th width character position. Can be set to 0 for any
       * cell. For this case the cell is not split! If the size of the vector does not correspond to
       * the size of the table, an error is thrown.
       *
       * @param widths vector of unsigned ints, defining the characters that should fit into each
       * cell.
       */
      void SetWidths(const std::vector<unsigned>& widths);
      /*!
       * The directives are entered as a key-value pair, which is also the way, how it is written in
       * a list-table. If the type of table is changed some day, the key-value pairs can probably
       * still be used to provide the controls needed in any other table style (csv-table, ...)
       *
       * @param key The control parameter of the table (e.g., width, widths, header-rows, align,...)
       * @param value The value of the control parameter, given as a string no matter what it
       *              actually is.
       */
      void AddDirective(const std::string& key, const std::string& value);
      /*!
       * A getter method for the current number of rows in the table
       *
       * @return Number of rows of the table (that is, the size of the private vector tablerows_
       */
      unsigned GetRows() const;
      /*!
       * Here the table is printed as a restructuredText list-table
       *
       * @param stream The output file stream of the restructuredText file.
       */
      void Print(std::ostream& stream) const;

     private:
      /*!
       * The construction parameter: Number of columns of the table
       */
      const unsigned tablewidth_;
      std::vector<unsigned> widths_;
      /*!
       * tablerows_ is a vector, that contains the content rows of table as each entry.
       *            Each row is again a vector of strings that contains the cells.
       */
      std::vector<std::vector<std::string>> tablerows_;
      /*!
       * directives_ is a mapping containing all table controls which are given for the table.
       */
      std::map<std::string, std::string> directives_;
    };

    /*!
      \brief Write a link target to the restructuredText stream

      \param[in] stream: stream for the restructuredText file
      \param[in] line: The link target to be printed

    */
    void WriteLinktarget(std::ostream& stream, const std::string& line);

    /*----------------------------------------------------------------------*/
    /*!
      \brief Write a header of a specific level to the restructuredText stream

      \param[in] stream: stream for the restructuredText file
      \param[in] level: header level
      \param[in] line: The link target to be printed

    */
    void WriteHeader(std::ostream& stream, unsigned level, const std::string& line);

    /*----------------------------------------------------------------------*/
    /*!
      \brief Write a paragraph to the restructuredText stream, with an optional additional
      indentation

      \param[in] stream: stream for the restructuredText file
      \param[in] line: The link target to be printed
      \param[in] indent: Optional indentation of the paragraph (default 0)

    */
    void WriteParagraph(std::ostream& stream, std::string& line, size_t indent = 0);
    /*----------------------------------------------------------------------*/
    /*!
      \brief Write a code block to the restructuredText stream

      \param[in] stream: stream for the restructuredText file
      \param[in] lines: A vector of strings that make up the code

    */
    void WriteCode(std::ostream& stream, const std::vector<std::string>& lines);
    /*----------------------------------------------------------------------*/
    /*!
      \brief Write a note block to the restructuredText stream

      \param[in] stream: stream for the restructuredText file
      \param[in] line: The link target to be printed

    */
    void WriteNote(std::ostream& stream, const std::string& paragraph);

    /*!
     *  \brief write all known material sections for readthedocs
     *
     *  \param[in] stream: stream for the restructuredText file
     *  \param[in] matlist: vector containing all material definitions
     */
    void WriteMaterialReference(std::ostream& stream,
        const std::vector<Teuchos::RCP<DRT::INPUT::MaterialDefinition>>& matlist);

    /*!
     *  \brief write the information of a single material to readthedocs
     *
     *  \param[in] stream: stream for the restructuredText file
     *  \param[in] material: single material definition
     */
    void WriteSingleMaterialReadTheDocs(
        std::ostream& stream, const Teuchos::RCP<INPUT::MaterialDefinition> material);
    /*!
     *  \brief write all parameters of the header sections for readthedocs
     *
     *  \param[in] stream: stream for the restructuredText file
     *  \param[in] list: vector containing all paramters in the current section
     *  \param[in] parentname: name of the parent section (initially empty string)
     */
    /// print flag sections of dat file with given list
    void WriteHeaderReference(
        std::ostream& stream, const Teuchos::ParameterList& list, std::string parentname = "");

    /*!
     *  write all known condition sections including explanations (if available) to a .rst file for
     * ReadTheDocs
     *
     *  @param[in] stream restructuredText file for prescribed conditions.
     *  @param[in] condlist List of prescribed conditions to be written to that file
     */
    void WriteConditionsReference(std::ostream& stream,
        const std::vector<Teuchos::RCP<DRT::INPUT::ConditionDefinition>>& condlist);

    /*!
     *  write a single condition including explanations (if available) to a .rst file for
     * ReadTheDocs
     *
     *  @param[in] stream restructuredText file for prescribed conditions.
     *  @param[in] condition Single prescribed condition to be written to that file
     */
    void WriteSingleConditionReadTheDocs(
        std::ostream& stream, const Teuchos::RCP<INPUT::ConditionDefinition> condition);

    /*!
     *  write a single contact law including explanations (if available) to a .rst file for
     * ReadTheDocs
     *
     *  @param[in] stream restructuredText file for prescribed contact law.
     *  @param[in] coconstlaw Single contact law to be written to that file
     */
    void WriteSingleContactLawReadTheDocs(std::ostream& stream,
        const Teuchos::RCP<CONTACT::CONSTITUTIVELAW::LawDefinition> coconstlaw);

    /*!
     *  write all known contact laws including explanations (if available) to a .rst file for
     * ReadTheDocs
     *
     *  @param[in] stream restructuredText file for prescribed contact law.
     */
    void WriteContactLawReference(std::ostream& stream,
        const std::vector<Teuchos::RCP<CONTACT::CONSTITUTIVELAW::LawDefinition>>& coconstlawlist);
    /*!
     *  write various other parameters including explanations (if available) to a .rst file for
     * ReadTheDocs
     *
     *  @param[in] stream restructuredText file for functions.
     */
    void WriteVariousReference(std::ostream& stream);
  }  // namespace RTD
}  // namespace DRT

std::ostream& operator<<(std::ostream& os, const DRT::RTD::Table& table);


#endif  // CREATE_RTDFILES_UTILS_H
