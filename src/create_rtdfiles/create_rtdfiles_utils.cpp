/*---------------------------------------------------------------------*/
/*! \file

\brief A collection of helper methods for the readthedocs parser

\level 0

*/
/*---------------------------------------------------------------------*/

#include "create_rtdfiles_utils.H"
#include "utils_exceptions.H"
#include "lib_utils_reader.H"
#include <boost/format.hpp>

namespace DRT
{
  namespace RTD
  {
    /*----------------------------------------------------------------------*/
    /*----------------------------------------------------------------------*/
    Table::Table(const unsigned &t) : tablewidth_(t)
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
    unsigned Table::GetRows() { return tablerows_.size(); }
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
    void WriteHeader(std::ostream &stream, const unsigned level, const std::string &line)
    {
      const std::vector<char> headerchar{'=', '-', '~', '^'};
      unsigned headerlength = line.length();
      stream << line << "\n";
      if (level < 0 or level > headerchar.size())
      {
        dserror("Header level for ReadTheDocs output must be [0,3], but is %i", level);
      }
      stream << std::string(headerlength, headerchar[level]);
      stream << "\n\n";
    }
    /*----------------------------------------------------------------------*/
    /*----------------------------------------------------------------------*/
    void WriteParagraph(std::ostream &stream, std::string &paragraph, const size_t indent)
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
  }  // namespace RTD
  /*----------------------------------------------------------------------*/
  /*----------------------------------------------------------------------*/
  void RTD::WriteMaterialReference(
      std::ostream &stream, std::vector<Teuchos::RCP<INPUT::MaterialDefinition>> &matlist)
  {
    RTD::WriteLinktarget(stream, "materialsreference");
    RTD::WriteHeader(stream, 0, "Material reference");

    std::vector<std::string> materialsectionstring{std::string(58, '-') + "MATERIALS"};
    RTD::WriteCode(stream, materialsectionstring);

    for (auto &material : matlist)
    {
      material->WriteReadTheDocs(stream);
    }
  }



  /*----------------------------------------------------------------------*/
  /*----------------------------------------------------------------------*/
  void RTD::WriteHeaderReference(
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
        if (name == INPUT::PrintEqualSign()) continue;
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
        linktarget = Teuchos::StrUtils::removeAllSpaces(UTILS::ToLower(linktarget));

        if (entry.isList())  // it is a section header
        {
          unsigned l = fullname.length();
          // write link:
          DRT::RTD::WriteLinktarget(stream, "SEC" + linktarget);
          // write section header
          unsigned level = (issubsection) ? 2 : 1;
          DRT::RTD::WriteHeader(stream, level, fullname);

          DRT::RTD::WriteParagraph(stream, doc);

          std::vector<std::string> codelines;
          codelines.push_back("--" + std::string(std::max<int>(65 - l, 0), '-') + fullname);
          DRT::RTD::WriteCode(stream, codelines);

          if (INPUT::NeedToPrintEqualSign(list))
          {
            DRT::RTD::WriteNote(stream,
                "   The parameters in this section need an equal sign (=) "
                "between the parameter name and its value!");
          }

          WriteHeaderReference(stream, list.sublist(name), fullname);
        }
        else  // it is a parameter entry
        {
          DRT::RTD::WriteLinktarget(stream, linktarget);

          const Teuchos::any &v = entry.getAny(false);
          boost::format parstring = boost::format{"**%s** | *default:* %s |break| %s"} % name %
                                    Teuchos::toString(v) % doc;
          std::string s = parstring.str();
          DRT::RTD::WriteParagraph(stream, s);
          if (validator != Teuchos::null)  // it can only take specific values
          {
            Teuchos::RCP<const Teuchos::Array<std::string>> values = validator->validStringValues();
            if (values != Teuchos::null)
            {
              stream << "   **Possible values:**\n\n";
              for (const auto val : *values) stream << "   - " << val << "\n";
              stream << "\n";
            }
          }
        }
      }
    }
  }

  /*----------------------------------------------------------------------*/
  /*----------------------------------------------------------------------*/
  void RTD::WriteConditionsReference(
      std::ostream &stream, std::vector<Teuchos::RCP<INPUT::ConditionDefinition>> &condlist)
  {
    RTD::WriteLinktarget(stream, "prescribedconditionreference");
    RTD::WriteHeader(stream, 0, "Prescribed Condition Reference");

    for (unsigned i = 0; i < condlist.size(); ++i)
    {
      condlist[i]->WriteReadTheDocs(stream);
    }
  }

}  // namespace DRT