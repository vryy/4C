# HEADER

class Header(object):
  def __init__(self, filetext):
    self.header     = False
    self.start      = ""
    self.file       = ""
    self.brief      = ""
    self.maintainer = ""
    self.level      = -1
    # parse the file header for the required tags
    blk = Header._extract_first_doxy_block(filetext)
    if len(blk) > 0:
      self.header = True
      # check for correct header start
      for line in blk:
        if (("/*" in line or "/!*" in line) and ("*/" not in line) and ("/*!" not in line)):
          start = line
          if len(start)>=1:
            self.start=start
          break
      # check for file tag
      for line in blk:
        if "\\file " in line:
          file = line.split("\\file ",1)[1].strip()
          if len(file) >= 1:
            self.file = file
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
          if len(maint) >= 5 and " " in maint:
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

  def has_header(self):
    return self.header

  def get_file(self):
    return self.file

  def get_start(self):
    return self.start

  def get_brief(self):
    return self.brief

  def get_maintainer(self):
    return self.maintainer

  def get_level(self):
    return self.level

  @staticmethod
  def get_example():
    return ["","An appropriate header block looks as follows:",
            "/*----------------------------------------------------------------------------*/",
            "/*!",
            "\\file fe_simulation.H",
            "",
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
    for line in cont.splitlines():
      if (not marker) and ("/*!" in line or "/**" in line or "/*" in line or "/!*" in line) and ("*/" not in line):
        marker = True
      if marker:
        blk.append(line)
        if "*/" in line:
          return blk
    return []

