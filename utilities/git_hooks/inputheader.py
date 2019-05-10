# HEADER
import json
import os
import sys

class Header(object):
  def __init__(self, filetext):
    self.header     = False
    self.maintainer = ""
    self.compliant_maintainer=""
    # parse the file header for the required tags
    blk = Header._extract_header(filetext)
    if len(blk) > 0:
      self.header = True
      # check for maintainer tag
      for line in blk:
        if line.startswith("// Maintainer: "):
          maint = line[len("// Maintainer: "):].strip()
          if len(maint) >= 1:
            self.maintainer = maint
          break
    #Check for compliant maintainers
    with open(os.path.join(sys.path[0],'baci_developers.json'),'r') as f:
        developers = json.load(f)
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

  def get_maintainer(self):
    return self.maintainer

  def get_compliant_maintainer(self):
    return self.compliant_maintainer

  @staticmethod
  def get_example():
    return ["","An appropriate header must be before any section and contain:",
            "// Maintainer: Max Mustermann"]

  @staticmethod
  def _extract_header(cont):
    blk = []
    marker = False
    for line in cont.splitlines():
      if line.startswith("--"):
        break
      blk.append(line)
    return blk

