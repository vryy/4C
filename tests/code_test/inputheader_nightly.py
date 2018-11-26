# HEADER

class Header(object):
  def __init__(self, filetext):
    self.header     = False
    self.maintainer = ""
    # parse the file header for the required tags
    blk = Header._extract_header(filetext)
    if len(blk) > 0:
      self.header = True
      # check for maintainer tag
      for line in blk:
        if line.startswith("// Maintainer: "):
          maint = line[len("// Maintainer: "):].strip()
          if len(maint) >= 5 and " " in maint:
            self.maintainer = maint
          break

  def has_header(self):
    return self.header

  def get_maintainer(self):
    return self.maintainer

  @staticmethod
  def get_example():
    return ["","An appropriate header must be before any section and contain:",
            "// Maintainer: Max Mustermann"]

  @staticmethod
  def _extract_header(cont):
    blk = []
    marker = False
    for line in cont:
      if line.startswith("--"):
        break
      blk.append(line)
    return blk

