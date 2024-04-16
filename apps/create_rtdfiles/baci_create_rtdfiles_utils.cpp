/*---------------------------------------------------------------------*/
/*! \file

\brief A collection of helper methods for the readthedocs parser

\level 0

*/
/*---------------------------------------------------------------------*/

#include "baci_create_rtdfiles_utils.hpp"

#include "baci_discretization_fem_general_cell_type_traits.hpp"
#include "baci_io_linedefinition.hpp"
#include "baci_io_utils_reader.hpp"
#include "baci_lib_utils_createdis.hpp"
#include "baci_utils_exceptions.hpp"

#include <boost/algorithm/string.hpp>
#include <boost/format.hpp>
#include <unistd.h>

BACI_NAMESPACE_OPEN


namespace RTD
{
  /*----------------------------------------------------------------------*/
  /*----------------------------------------------------------------------*/
  Table::Table(const unsigned &size) : tablewidth_(size)
  {
    for (unsigned i = 0; i < tablewidth_; ++i)
    {
      widths_.push_back(0);
    }
  }
  /*----------------------------------------------------------------------*/
  /*----------------------------------------------------------------------*/
  void Table::AddRow(const std::vector<std::string> &row)
  {
    if (row.size() != tablewidth_)
    {
      dserror("Trying to add %i row elements into a table with %i rows", row.size(), tablewidth_);
    }
    tablerows_.push_back(row);
  }
  /*----------------------------------------------------------------------*/
  /*----------------------------------------------------------------------*/
  void Table::SetWidths(const std::vector<unsigned> &widths)
  {
    if (widths.size() != tablewidth_)
    {
      dserror(
          "Number of given cell widths (%i) "
          "does not correspond to the number of rows (%i)",
          widths.size(), tablewidth_);
    }
    widths_ = widths;
  }
  /*----------------------------------------------------------------------*/
  /*----------------------------------------------------------------------*/
  void Table::AddDirective(const std::string &key, const std::string &value)
  {
    directives_[key] = value;
  }
  /*----------------------------------------------------------------------*/
  /*----------------------------------------------------------------------*/
  unsigned Table::GetRows() const { return tablerows_.size(); }
  /*----------------------------------------------------------------------*/
  /*----------------------------------------------------------------------*/
  void Table::Print(std::ostream &stream) const
  {
    // if the widths are not set, they are currently set to 100/tablewidth_
    unsigned defaultcolsize = 100 / tablewidth_;
    bool isWidthDirectiveGiven = false;
    stream << ".. list-table::\n";
    // write directives
    for (const auto &directive : directives_)
    {
      stream << "   :" << directive.first << ": " << directive.second << "\n";
      if (directive.first.substr(0, 5) == "width") isWidthDirectiveGiven = true;
    }
    // add the cell widths to the directives if not already given
    if (!isWidthDirectiveGiven)
    {
      stream << "   :widths: ";
      unsigned wd;
      for (auto itw = widths_.begin(); itw != widths_.end(); ++itw)
      {
        if (itw != widths_.begin()) stream << ",";
        wd = (*itw == 0) ? defaultcolsize : *itw;
        stream << wd;
      }
      stream << "\n";
    }
    stream << "\n";
    //
    // now write table content (split if necessary, i.e., more characters than given in widths_)
    for (const auto &tablerow : tablerows_)
    {
      for (unsigned i = 0; i < tablewidth_; ++i)
      {
        std::string cellstring = (i == 0) ? "   * - " : "     - ";
        if ((widths_[i] != 0) and (tablerow[i].length() > widths_[i]))
        {
          std::string cellstringPart = tablerow[i];
          std::size_t spacepos = cellstringPart.rfind(" ", widths_[i]);
          if (spacepos < cellstringPart.npos)
          {
            cellstring += cellstringPart.substr(0, spacepos) + " |break| \n";
            cellstringPart = cellstringPart.substr(spacepos + 1);
            // print the rest of the description with two empty columns before
            while (cellstringPart.length() > widths_[i])
            {
              spacepos = cellstringPart.rfind(" ", widths_[i]);
              if (spacepos == cellstringPart.npos) break;
              cellstring += "       " + cellstringPart.substr(0, spacepos) + " |break| \n";
              cellstringPart = cellstringPart.substr(spacepos + 1);
            }
            cellstring += "       ";
          }
          cellstring += cellstringPart;
        }
        else
        {
          cellstring += tablerow[i];
        }
        stream << cellstring << "\n";
      }
    }
    stream << "\n";
  }
  /*----------------------------------------------------------------------*/
  /*----------------------------------------------------------------------*/
  void WriteLinktarget(std::ostream &stream, const std::string &line)
  {
    stream << ".. _" << line << ":\n\n";
  }
  /*----------------------------------------------------------------------*/
  /*----------------------------------------------------------------------*/
  void WriteHeader(std::ostream &stream, unsigned level, const std::string &line)
  {
    const std::vector<char> headerchar{'=', '-', '~', '^'};
    unsigned headerlength = line.length();
    stream << line << "\n";
    if (level > headerchar.size())
    {
      dserror("Header level for ReadTheDocs output must be [0,3], but is %i", level);
    }
    stream << std::string(headerlength, headerchar[level]);
    stream << "\n\n";
  }
  /*----------------------------------------------------------------------*/
  /*----------------------------------------------------------------------*/
  void WriteParagraph(std::ostream &stream, std::string &paragraph, size_t indent)
  {
    size_t mathstartpos = paragraph.find("\f$");
    size_t mathendpos = 0;
    while (mathstartpos != paragraph.npos)
    {
      mathendpos = paragraph.find("\f$", mathstartpos + 1);
      if (mathendpos == paragraph.npos)
      {
        dserror(
            "Math tags in a ReadTheDocs paragraph must occur pairwise. "
            "Error found in: \n" +
            paragraph);
      }
      paragraph.replace(mathendpos, 2, "`");
      paragraph.replace(mathstartpos, 2, ":math:`");
      mathstartpos = paragraph.find("\f$");
    }
    stream << std::string(" ", indent) << paragraph << "\n\n";
  }
  /*----------------------------------------------------------------------*/
  /*----------------------------------------------------------------------*/
  void WriteCode(std::ostream &stream, const std::vector<std::string> &lines)
  {
    stream << "::\n\n";
    for (const auto &line : lines)
    {
      stream << "   " << line << "\n";
    }
    stream << "\n";
  }
  /*----------------------------------------------------------------------*/
  /*----------------------------------------------------------------------*/
  void WriteNote(std::ostream &stream, const std::string &paragraph)
  {
    stream << ".. note::\n\n";
    stream << "   " << paragraph << "\n\n";
  }

  /*----------------------------------------------------------------------*
   *----------------------------------------------------------------------*/
  std::ostream &operator<<(std::ostream &stream, const Table &table)
  {
    table.Print(stream);
    return stream;
  }

  /*----------------------------------------------------------------------*/
  /*----------------------------------------------------------------------*/
  void WriteCelltypeReference(std::ostream &stream)
  {
    WriteLinktarget(stream, "celltypes");
    WriteHeader(stream, 1, "Cell types");

    // We run the loop over the cell types four times to sort the cell types after their dimension
    for (unsigned outputdim = 0; outputdim < 4; ++outputdim)
    {
      WriteLinktarget(stream, boost::str(boost::format("%1dD_cell_types") % outputdim));
      WriteHeader(stream, 2, boost::str(boost::format("%1dD cell types") % outputdim));

      for (auto celltype : CORE::FE::celltype_array<CORE::FE::all_physical_celltypes>)
      {
        std::string celltypename = CORE::FE::CellTypeToString(celltype);
        // Skip the cell type if it has not the desired dimension
        const unsigned celldimension = CORE::FE::getDimension(celltype);
        if (celldimension != outputdim) continue;

        std::string celltypelinkname = boost::algorithm::to_lower_copy(celltypename);
        WriteLinktarget(stream, celltypelinkname);
        WriteHeader(stream, 3, celltypename);

        std::stringstream celltypeinfostream;
        celltypeinfostream << "- Nodes: " << CORE::FE::getNumberOfElementNodes(celltype)
                           << std::endl;
        celltypeinfostream << "- Dimension: " << celldimension << std::endl;
        if (CORE::FE::getOrder(celltype, -1) >= 0)
        {
          celltypeinfostream << "- Shape function order (element): "
                             << CORE::FE::getDegree(celltype) << std::endl;
          celltypeinfostream << "- Shape function order (edges): " << CORE::FE::getOrder(celltype)
                             << std::endl;
        }
        std::string celltypeinformation = celltypeinfostream.str();
        WriteParagraph(stream, celltypeinformation);

        if (celldimension >= 2)
        {
          const std::string figurename("reference_images/" + celltypename + ".png");
          std::string captionstring = "**" + celltypename + ":** ";
          if (celldimension == 2)
            captionstring += "Line and node numbering";
          else
            captionstring += "Left: Line and node numbering, right: Face numbering";
          std::string figureincludestring = ".. figure:: " + figurename + "\n";
          figureincludestring += "    :alt: Figure not available for " + celltypename + "\n";
          figureincludestring += "    :width: ";
          figureincludestring += (outputdim == 3) ? "100%" : "50%";
          figureincludestring += "\n\n";
          figureincludestring += "    " + captionstring;
          WriteParagraph(stream, figureincludestring);
        }
      }
    }
  }

  /*----------------------------------------------------------------------*/
  /*----------------------------------------------------------------------*/
  void WriteMaterialReference(
      std::ostream &stream, const std::vector<Teuchos::RCP<INPUT::MaterialDefinition>> &matlist)
  {
    WriteLinktarget(stream, "materialsreference");
    WriteHeader(stream, 0, "Material reference");

    std::vector<std::string> materialsectionstring{std::string(58, '-') + "MATERIALS"};
    WriteCode(stream, materialsectionstring);

    for (auto &material : matlist)
    {
      WriteSingleMaterialReadTheDocs(stream, material);
    }
    //
    // adding the section for the CLONING MATERIAL MAP
    WriteLinktarget(stream, "cloningmaterialsreference");
    WriteHeader(stream, 0, "Cloning material reference");
    const INPUT::Lines lines = DRT::UTILS::ValidCloningMaterialMapLines();
    std::stringstream cloningMatStream;
    lines.Print(cloningMatStream);
    const std::vector<std::string> cloningMatList = DRT::UTILS::Split(cloningMatStream.str(), "\n");

    WriteCode(stream, cloningMatList);
  }


  void WriteSingleMaterialReadTheDocs(
      std::ostream &stream, const Teuchos::RCP<INPUT::MaterialDefinition> material)
  {
    /* Each entry consists of a number of fields:
    - header
    - description
    - code line
    - parameter description */

    // the Material title
    WriteLinktarget(stream, material->Name());
    WriteHeader(stream, 1, material->Name());

    // the description of the material
    std::string materialDescription = material->Description();
    WriteParagraph(stream, materialDescription);

    // the material line as it occurs in the dat file
    std::string parameter = "   MAT <matID>  " + material->Name();
    std::vector<std::string> materialcode;
    //
    // Also: create the table from the parameter descriptions (based on the so-called
    // separatorComponents) table header
    const unsigned tablesize = 3;
    Table parametertable(tablesize);
    std::vector<std::string> tablerow(tablesize);
    tablerow = {"Parameter", "optional", "Description"};
    parametertable.AddRow(tablerow);

    for (auto &parameterterm : material->Inputline())
    {
      if (auto *separator = dynamic_cast<INPUT::SeparatorComponent *>(parameterterm.get()))
      {
        parametertable.AddRow(separator->WriteReadTheDocsTableRow());

        if (parameter.length() > 60)
        {
          parameter += " \\";
          materialcode.push_back(parameter);
          parameter = "   ";
        }
      }
      std::ostringstream parameterstream;
      parameterterm->DefaultLine(parameterstream);
      parameter += " " + parameterstream.str();
    }
    materialcode.push_back(parameter);
    WriteCode(stream, materialcode);
    //
    // Now printing the parameter table
    parametertable.SetWidths({10, 10, 50});
    parametertable.AddDirective("header-rows", "1");

    if (parametertable.GetRows() == 1)
    {
      tablerow = {"no parameters", "", ""};
      parametertable.AddRow(tablerow);
    }
    parametertable.Print(stream);

    return;
  }
  /*----------------------------------------------------------------------*/
  /*----------------------------------------------------------------------*/
  void WriteHeaderReference(
      std::ostream &stream, const Teuchos::ParameterList &list, std::string parentname)
  {
    // prevent invalid ordering of parameters caused by alphabetical output:
    // in the first run, print out all list elements that are not a sublist
    // in the second run, do the recursive call for all the sublists in the list

    for (int j = 0; j < 2; ++j)
    {
      // bool loop_isEntry = (j == 0);
      bool loop_isList = (j == 1);
      for (Teuchos::ParameterList::ConstIterator it = list.begin(); it != list.end(); ++it)
      {
        const Teuchos::ParameterEntry &entry = list.entry(it);
        if (entry.isList() != loop_isList) continue;
        const std::string &name = list.name(it);
        Teuchos::RCP<const Teuchos::ParameterEntryValidator> validator = entry.validator();

        std::string doc = (entry.docString() == "") ? "no description yet" : entry.docString();

        std::string fullname = parentname;
        bool issubsection = false;
        if (fullname != "")
        {
          fullname += "/";
          issubsection = true;
        }
        fullname += name;
        std::string linktarget = boost::algorithm::replace_all_copy(fullname, "/", "_");
        linktarget = Teuchos::StrUtils::removeAllSpaces(DRT::UTILS::ToLower(linktarget));

        if (entry.isList())  // it is a section header
        {
          unsigned l = fullname.length();
          // write link:
          WriteLinktarget(stream, "SEC" + linktarget);
          // write section header
          unsigned level = (issubsection) ? 2 : 1;
          WriteHeader(stream, level, fullname);

          WriteParagraph(stream, doc);

          std::vector<std::string> codelines;
          codelines.push_back("--" + std::string(std::max<int>(65 - l, 0), '-') + fullname);
          WriteCode(stream, codelines);

          if (INPUT::NeedToPrintEqualSign(list.sublist(name)))
          {
            WriteNote(stream,
                "   The parameters in this section need an equal sign (=) "
                "between the parameter name and its value!");
          }

          WriteHeaderReference(stream, list.sublist(name), fullname);
        }
        else  // it is a parameter entry
        {
          WriteLinktarget(stream, linktarget);

          const Teuchos::any &v = entry.getAny(false);
          boost::format parstring = boost::format{"**%s** | *default:* %s |break| %s"} % name %
                                    Teuchos::toString(v) % doc;
          std::string s = parstring.str();
          WriteParagraph(stream, s);
          if (validator != Teuchos::null)  // it can only take specific values
          {
            Teuchos::RCP<const Teuchos::Array<std::string>> values = validator->validStringValues();
            if (values != Teuchos::null)
            {
              stream << "   **Possible values:**\n\n";
              for (auto val : *values) stream << "   - " << val << "\n";
              stream << "\n";
            }
          }
        }
      }
    }
  }

  /*----------------------------------------------------------------------*/
  /*----------------------------------------------------------------------*/
  void WriteConditionsReference(
      std::ostream &stream, const std::vector<Teuchos::RCP<INPUT::ConditionDefinition>> &condlist)
  {
    WriteLinktarget(stream, "prescribedconditionreference");
    WriteHeader(stream, 0, "Prescribed Condition Reference");

    for (auto &condition : condlist)
    {
      WriteSingleConditionReadTheDocs(stream, condition);
    }
  }


  void WriteSingleConditionReadTheDocs(
      std::ostream &stream, const Teuchos::RCP<INPUT::ConditionDefinition> condition)
  {
    /* Each entry consists of a number of fields:
    - Part 1: link target and header
    - Part 2: description
    - Part 3: code lines
    - Part 4: table for description and admissible values of string
      and conditionComponentBundleSelector parameters
    - Part 5: finally additional code lines for model specific parameters
      or the complex ConditionComponentBundleSelectors
    */

    std::string sectionname = condition->SectionName();
    const std::string sectionlinktarget =
        Teuchos::StrUtils::removeAllSpaces(DRT::UTILS::ToLower(sectionname));
    //
    // boundary condition header
    //
    /*------ PART 1 --------------------------
     * Boundary condition header (incl. link target)
     */
    // link target line
    WriteLinktarget(stream, sectionlinktarget);
    // condition name as section header
    WriteHeader(stream, 1, sectionname);

    /*------ PART 2 -------------------------
     * boundary condition description string
     */
    std::string descriptionline =
        (condition->Description() == "") ? "no description yet" : condition->Description();
    WriteParagraph(stream, descriptionline);

    /*------ PART 3 -------------------------
     * boundary condition input lines
     * In this section, the table for parameter description is filled as well.
     */
    // First line: condition name as a section
    unsigned l = sectionname.length();
    std::vector<std::string> conditioncode{
        "--" + std::string(std::max<int>(65 - l, 0), '-') + sectionname};
    // second line: geometry type
    std::string name;
    switch (condition->GeometryType())
    {
      case DRT::Condition::Point:
        conditioncode.push_back("DPOINT  0");
        break;
      case DRT::Condition::Line:
        conditioncode.push_back("DLINE  0");
        break;
      case DRT::Condition::Surface:
        conditioncode.push_back("DSURF  0");
        break;
      case DRT::Condition::Volume:
        conditioncode.push_back("DVOL  0");
        break;
      default:
        dserror("geometry type unspecified");
        break;
    }
    // Collecting information for the final code line (conditioncodeline)
    // Also:
    // store admissible values for string parameters (vector<string> parametertable) -> Part 4
    // store options for the CondCompBundles (vector<string> condCompStrings) -> Part 5
    std::string conditioncodeline = "E <setnumber> -";
    bool isNewlinePossible = false;
    //
    // also: Generate the table rows for admissible values of string parameters
    const unsigned tablesize = 3;
    Table parametertable(tablesize);
    std::vector<std::string> tablerow = {"Parameter", "Default", "Admissible values"};
    std::vector<std::string> condCompStrings;
    std::string condCompName("");
    parametertable.AddRow(tablerow);
    for (auto &condparameter : condition->Inputline())
    {
      // newline after some 60 characters, but no newline after a separator condition
      if (isNewlinePossible)
      {
        conditioncode.push_back(conditioncodeline + " \\ ");
        conditioncodeline = "   ";  // start a new line
      }
      conditioncodeline += " " + condparameter->WriteReadTheDocs();
      isNewlinePossible = (conditioncodeline.length() > 60);
      if (auto *previousparameter = dynamic_cast<INPUT::SeparatorComponent *>(condparameter.get()))
      {
        previousparameter->GetOptions();  // just needed to prevent an unusedVariable warning
        isNewlinePossible = false;
      }
      // If the component is a string component, store the admissible parameters in the table:
      if (auto *stringComponent = dynamic_cast<INPUT::SelectionComponent *>(condparameter.get()))
      {
        tablerow[0] = stringComponent->Name();
        std::ostringstream parametercell;
        stringComponent->DefaultLine(parametercell);
        tablerow[1] = parametercell.str();
        Teuchos::Array<std::string> datfilevalues = stringComponent->GetOptions();
        tablerow[2] = boost::algorithm::join(datfilevalues, ", ");
        parametertable.AddRow(tablerow);
      }
      // if the component is a bundleselector (bundle of variables following a string keyword):
      if (auto *compBundleSelector = dynamic_cast<INPUT::SwitchComponent *>(condparameter.get()))
      {
        condCompName = compBundleSelector->Name();
        std::vector<std::string> bundle = compBundleSelector->WriteReadTheDocsLines();
        condCompStrings.insert(condCompStrings.end(), bundle.begin(), bundle.end());
        tablerow[0] = condCompName;
        Teuchos::Array<std::string> datfilevalues = compBundleSelector->GetOptions();
        tablerow[1] = datfilevalues[0];
        tablerow[2] = boost::algorithm::join(datfilevalues, ", ");
        parametertable.AddRow(tablerow);
      }
    }
    // Now write the complete code of this condition to the readthedocs file
    conditioncode.push_back(conditioncodeline);
    WriteCode(stream, conditioncode);

    /*------ PART 4 -------------------------
     * Now write a table for the options of the string variables, if any have been stored above
     */
    if (parametertable.GetRows() > 1)
    {
      std::string optionheaderstring("**String options:**");
      WriteParagraph(stream, optionheaderstring);
      // table header for the options in string parameters
      parametertable.SetWidths({0, 0, 50});
      parametertable.AddDirective("header-rows", "1");

      parametertable.Print(stream);
    }

    /*------ PART 5 -------------------------
     * Finally add the model specific parameters for the complex ConditionComponentBundleSelectors
     */
    if (condCompStrings.size() > 0)
    {
      std::string optionheaderstring =
          "The following parameter sets are possible for `<" + condCompName + ">`:";
      WriteParagraph(stream, optionheaderstring);
      conditioncode.clear();
      for (auto &condCompString : condCompStrings)
      {
        conditioncode.push_back(condCompString);
      }
      WriteCode(stream, conditioncode);
    }
    return;
  }

  /*----------------------------------------------------------------------*/
  /*----------------------------------------------------------------------*/
  void WriteContactLawReference(std::ostream &stream,
      const std::vector<Teuchos::RCP<CONTACT::CONSTITUTIVELAW::LawDefinition>> &coconstlawlist)
  {
    WriteLinktarget(stream, "contactconstitutivelawreference");
    WriteHeader(stream, 0, "Contact Constitutive Law Reference");


    std::vector<std::string> contactlawsectionstring{
        std::string(43, '-') + "CONTACT CONSTITUTIVE LAW"};
    WriteCode(stream, contactlawsectionstring);

    for (auto &contactlaw : coconstlawlist)
    {
      WriteSingleContactLawReadTheDocs(stream, contactlaw);
    }
  }


  /*----------------------------------------------------------------------*/
  /*----------------------------------------------------------------------*/
  void WriteSingleContactLawReadTheDocs(
      std::ostream &stream, const Teuchos::RCP<CONTACT::CONSTITUTIVELAW::LawDefinition> contactlaw)
  {
    /* Each entry consists of a number of fields:
    - header
    - description
    - code line
    - parameter description */

    // the Law title
    WriteLinktarget(stream, contactlaw->Name());
    WriteHeader(stream, 1, contactlaw->Name());

    // the description of the contact law
    std::string contactDescription = contactlaw->Description();
    WriteParagraph(stream, contactDescription);

    // the material line as it occurs in the dat file
    std::string parameter = "LAW <lawID>   " + contactlaw->Name();
    std::vector<std::string> contactlawCode;
    //
    // Also: create the table from the parameter descriptions (based on the so-called
    // separatorComponents) table header
    const unsigned tablesize = 3;
    Table parametertable(tablesize);
    std::vector<std::string> tablerow(tablesize);
    tablerow = {"Parameter", "optional", "Description"};
    parametertable.AddRow(tablerow);

    for (auto &parameterterm : contactlaw->Inputline())
    {
      if (auto *separator = dynamic_cast<INPUT::SeparatorComponent *>(parameterterm.get()))
      {
        parametertable.AddRow(separator->WriteReadTheDocsTableRow());

        if (parameter.length() > 60)
        {
          parameter += " \\";
          contactlawCode.push_back(parameter);
          parameter = "   ";
        }
      }
      std::ostringstream parameterstream;
      parameterterm->DefaultLine(parameterstream);
      parameter += " " + parameterstream.str();
    }
    contactlawCode.push_back(parameter);
    WriteCode(stream, contactlawCode);
    //
    // Now printing the parameter table
    parametertable.SetWidths({10, 10, 50});
    parametertable.AddDirective("header-rows", "1");

    if (parametertable.GetRows() == 1)
    {
      tablerow = {"no parameters", "", ""};
      parametertable.AddRow(tablerow);
    }
    parametertable.Print(stream);

    return;
  }

  /*----------------------------------------------------------------------*/
  /*----------------------------------------------------------------------*/
  void WriteVariousReference(std::ostream &stream)
  {
    //
    // adding the sections for the RESULT DESCRIPTION
    {
      WriteLinktarget(stream, "restultdescriptionreference");
      WriteHeader(stream, 0, "Result description reference");
      DRT::ResultTestManager resulttestmanager;
      const INPUT::Lines lines = resulttestmanager.ValidResultLines();
      std::string sectionDescription = lines.Description();
      WriteParagraph(stream, sectionDescription);
      std::stringstream resultDescriptionStream;
      lines.Print(resultDescriptionStream);
      const std::vector<std::string> resultDescriptionList =
          DRT::UTILS::Split(resultDescriptionStream.str(), "\n");
      WriteCode(stream, resultDescriptionList);
    }
    //
    // adding the sections for the FUNCTION
    {
      WriteLinktarget(stream, "functionreference");
      WriteHeader(stream, 0, "Functions reference");
      const auto lines = CORE::UTILS::FunctionManager().ValidFunctionLines();
      std::string sectionDescription = lines.Description();
      WriteParagraph(stream, sectionDescription);
      std::stringstream functionStream;
      lines.Print(functionStream);
      const std::vector<std::string> functionList = DRT::UTILS::Split(functionStream.str(), "\n");
      WriteCode(stream, functionList);
    }
  }
  void WriteYamlCellTypeInformation(std::ostream &yamlfile)
  {
    for (auto celltype : CORE::FE::celltype_array<CORE::FE::all_physical_celltypes>)
    {
      std::string celltypename = CORE::FE::CellTypeToString(celltype);
      std::string yamlcelltypestring = celltypename + ":\n";
      // 0. information: dimension of the element
      yamlcelltypestring +=
          "  dimension: " + std::to_string(CORE::FE::getDimension(celltype)) + "\n";
      // 1. information: nodal coordinates
      CORE::LINALG::SerialDenseMatrix coordmap;
      try
      {
        coordmap = CORE::FE::getEleNodeNumbering_nodes_paramspace(celltype);
      }
      catch (...)
      {
        std::cout << "could not read coords\n";
        continue;
      }
      const unsigned num_nodes = coordmap.numCols();
      yamlcelltypestring += "  nodes:\n";
      for (unsigned node = 0; node < num_nodes; ++node)
      {
        yamlcelltypestring += "    - [";
        for (int indx = 0; indx < coordmap.numRows(); ++indx)
        {
          if (indx > 0) yamlcelltypestring += ",";
          yamlcelltypestring += boost::str(boost::format("%6.2f") % coordmap[node][indx]);
        }
        yamlcelltypestring += "]\n";
      }
      // 2. information: line vectors of internal node numbers
      bool nodes_exist = true;
      std::vector<std::vector<int>> linevector;
      try
      {
        linevector = CORE::FE::getEleNodeNumberingLines(celltype);
      }
      catch (...)
      {
        std::cout << "could not read lines\n";
        continue;
      }
      yamlcelltypestring += "  lines:\n";
      for (auto line : linevector)
      {
        yamlcelltypestring += "    - [";
        for (size_t indx = 0; indx < line.size(); ++indx)
        {
          if (indx > 0) yamlcelltypestring += ",";
          yamlcelltypestring += boost::str(boost::format("%3d") % line[indx]);
          if ((unsigned int)line[indx] >= num_nodes) nodes_exist = false;
        }
        yamlcelltypestring += "]\n";
      }
      if (not nodes_exist)
      {
        std::cout << "line nodes are not contained\n";
        continue;
      }
      // 3. information: surface vectors of internal node numbers (for 3D elements)
      if (CORE::FE::getDimension(celltype) == 3)
      {
        std::vector<std::vector<int>> surfacevector;
        try
        {
          surfacevector = CORE::FE::getEleNodeNumberingSurfaces(celltype);
        }
        catch (...)
        {
          std::cout << "could not read surfaces\n";
          continue;
        }
        yamlcelltypestring += "  surfaces:\n";
        for (auto surface : surfacevector)
        {
          yamlcelltypestring += "    - [";
          for (size_t indx = 0; indx < surface.size(); ++indx)
          {
            if (indx > 0) yamlcelltypestring += ",";
            yamlcelltypestring += boost::str(boost::format("%3d") % surface[indx]);
            if ((unsigned)surface[indx] >= num_nodes) nodes_exist = false;
          }
          yamlcelltypestring += "]\n";
        }
        if (not nodes_exist)
        {
          std::cout << "surface nodes are not contained\n";
          continue;
        }
        // 4. information: vector of number of nodes for all surfaces
        std::vector<int> surfacecorners;
        try
        {
          surfacecorners = CORE::FE::getNumberOfFaceElementCornerNodes(celltype);
        }
        catch (...)
        {
          std::cout << "could not read surface corners\n";
          continue;
        }
        yamlcelltypestring += "  surfacecorners: [";
        for (size_t indx = 0; indx < surfacecorners.size(); ++indx)
        {
          if (indx > 0) yamlcelltypestring += ",";
          yamlcelltypestring += boost::str(boost::format("%3d") % surfacecorners[indx]);
        }
        yamlcelltypestring += "]\n";
      }
      std::cout << "Writing information on cell type " << celltypename << " to yaml file\n";
      yamlfile << yamlcelltypestring;
    }
  }
}  // namespace RTD
BACI_NAMESPACE_CLOSE
