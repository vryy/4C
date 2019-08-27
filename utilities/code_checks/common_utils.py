import subprocess, os, sys

#UTILS

def command_output(cmd):
  " Capture a command's standard output. "
  message = subprocess.Popen(cmd.split(), stdout=subprocess.PIPE).communicate()[0]
  return message

def path_split_all(path):
  folders = []
  while 1:
    path, file = os.path.split(path)
    if file == "" and (path == "" or path == os.path.sep):
      break
    if file != "":
      folders.append(file)
  folders.reverse()
  while folders[0] == '.':
    folders.pop(0)
  return folders

def is_unittest_file(fname):
  parts = path_split_all(fname)

  if len(parts) < 3:
    return False
  if parts[0] != 'Unittests':
    return False
  if parts[1] != 'src':
    return False
    
  return os.path.splitext(parts[-1])[1] in ".h .H .hpp".split()

def is_unittest_cmakefile(fname):
  parts = path_split_all(fname)

  if parts[0] != 'Unittests':
    return False

  return parts[-1] == 'CMakeLists.txt'

def is_source_file(fname):
  return os.path.splitext(fname)[1] in ".c .cpp .cxx .h .H .hpp".split()

def is_support_file(fname):
  return os.path.splitext(fname)[1] in ".dat .cmake .config .md".split()

def path_contains(test, path):
  head, tail = os.path.split(path)
  if test == tail or len(head) <= 0:
    return (test == tail)
  return path_contains(test, head)

def is_input_file(fname):
  return path_contains("Input", fname) and os.path.splitext(fname)[1] == ".dat"

def is_checked_file(fname):
  return is_source_file(fname) or is_support_file(fname)

def files_changed(look_cmd):
  """ List the files added or updated by this transaction.
  """
  def filename(line):
    return line
  looked_files  = command_output(look_cmd).split("\n")
  changed_files = [filename(line) for line in looked_files]
  return changed_files

def file_contents(filename):
  " Return a file's contents for this transaction. "
  with open (filename) as myfile:
    output=myfile.readlines()
  return output

def pretty_print_error_stderr(allerrors):
  _common_write_error(allerrors, sys.stderr, ['Your commit was rejected due to the following reason(s):'])

def pretty_print_error_file(allerrors, errfile):
  with open(errfile, 'w') as file:
    _common_write_error(allerrors, file)

def _common_write_error(allerrors, deststream, header = None):
  max_width = max([56]+[len(line) for line in allerrors])

  print_line = lambda line : deststream.write('E '+line+' '*(max_width-len(line))+' E\n')

  # print header
  deststream.write('\n'+'E'*(max_width+4)+'\nE'+' '*(max_width+2)+'E\n')
  if header is not None:
    for line in header:
      print_line(line)
    print_line('')
  
  # print body
  for line in allerrors:
    print_line(line)
  
  # footer
  deststream.write('E'+' '*(max_width+2)+'E\n'+'E'*(max_width+4)+'\n')

def is_alphabetically(items, err_entries = None):
  if len(items) <= 1:
    return True

  for i in range(1, len(items)):
    if items[i - 1] > items[i]:
      if err_entries is not None:
        err_entries.append(items[i - 1])
        err_entries.append(items[i])
      return False

  return True
  
