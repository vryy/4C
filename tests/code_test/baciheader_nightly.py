# HEADER
import json
import os
import sys

class Header(object):
  def __init__(self, filetext):
    self.header     = False
    self.start      = ""
    self.brief      = ""
    self.maintainer = ""
    self.compliant_maintainer = ""
    self.level      = -1
    # parse the file header for the required tags
    blk = Header._extract_first_doxy_block(filetext)
    if len(blk) > 0:
      self.header = True
   # check for correct header start
      for line in blk:
        if (("/*!" in line or "/**" in line)  and ("*/" not in line)):
          start = line
          if len(start) >= 1:
            self.start=start
          break
     # check for brief tag
      for line in blk:
        if "\\brief " in line:
          brief = line.split("\\brief ",1)[1].strip()
          if len(brief) >= 5:
            self.brief = brief
          break
      # check for maintainer tag
      for line in blk:
        if "\\maintainer " in line:
          maint = line.split("\\maintainer ",1)[1].strip()
          if len(maint) >= 1:
            self.maintainer = maint
          break
      # check for level tag
      for line in blk:
        if "\\level " in line:
          lvl = line.split("\\level ",1)[1].strip()
          try:
            self.level = int(lvl)
          except ValueError:
            pass
          break
      # check for compliant developers as maintainers
    with open(os.path.join(sys.path[0],'../../utilities/git_hooks/baci_developers.json'),'r') as f:
        developers = json.load(f)
    #maintainer_list = self.maintainer.split(",")
    #maintainer_list = [item.strip() for item in maintainer_list]
    maintainer_clean = [self.maintainer.strip()]
    if self.maintainer:
      flag = False
      compliant_names = [item['name'] for item in developers]
      flag = set(maintainer_clean).issubset(compliant_names)
      if (flag): self.compliant_maintainer = maintainer_clean
    else:
      self.compliant_maintainer='dummy_string'


  def has_header(self):
    return self.header

  def get_start(self):
    return self.start

  def get_brief(self):
    return self.brief

  def get_maintainer(self):
    return self.maintainer

  def get_compliant_maintainer(self):
   return self.compliant_maintainer

  def get_level(self):
    return self.level

  @staticmethod
  def get_example():
    return ["","An appropriate header block looks as follows:",
            "/*----------------------------------------------------------------------------*/",
            "/*!",
            "\\brief This header provides the interface for all FE simulations",
            "",
            "\\level 3",
            "",
            "\\maintainer Max Mustermann",
            "*/",
            "/*----------------------------------------------------------------------------*/"]

  @staticmethod
  def _extract_first_doxy_block(cont):
    blk = []
    marker = False
    for line in cont:
      if (not marker) and ("/*!" in line or "/**" in line or "/*" in line or "/!*") and ("*/" not in line):
        marker = True
      if marker:
        blk.append(line)
        if "*/" in line:
          return blk
    return []
