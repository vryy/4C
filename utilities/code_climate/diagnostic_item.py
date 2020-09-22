import os
import hashlib
import json


class DiagnosticItem:
    def __init__(self, src):
        self.name = None
        self.message = None
        self.file = None
        self.offset = None
        self.srcdir = src

        self.line_number = None
        self.previous_line = None
        self.line = None
        self.line_offset = None

    def hash(self):
        hash_source = "{0}{1}{2}{3}".format(
            self.name, self.message, self.get_relative_path(), self.offset
        )
        return str(hashlib.md5(hash_source.encode()).hexdigest())

    def global_hash(self, hash_set):
        hash_source = "{0}/{1}/{2}/{3}/{4}".format(
            self.name,
            self.message,
            self.get_relative_path(),
            self.line,
            self.previous_line,
        )
        hash_str = str(hashlib.md5(hash_source.encode()).hexdigest()[:9])
        for i in range(1000):
            my_hash = "{0}-{1}".format(hash_str, i)

            if my_hash not in hash_set:
                hash_set.add(my_hash)
                return my_hash

        return RuntimeError("Should never end here")

    def get_relative_path(self):
        return os.path.relpath(self.file, self.srcdir)

    def is_of_interest(self):
        return self.get_relative_path().startswith("src/") and self.file != ""

    def write(self, fout, hash_set):
        item = {
            "description": "{0} ({1})".format(self.message, self.name),
            "fingerprint": self.global_hash(hash_set),
            "location": {
                "path": self.get_relative_path(),
                "lines": {"begin": self.line_number},
            },
        }
        fout.write(json.dumps(item))

    @staticmethod
    def parse(object, src):
        diagnosticItem = DiagnosticItem(src)
        diagnosticItem.name = object["DiagnosticName"]
        diagnosticItem.message = object["DiagnosticMessage"]["Message"]

        diagnosticItem.file = object["DiagnosticMessage"]["FilePath"]
        diagnosticItem.offset = object["DiagnosticMessage"]["FileOffset"]

        return diagnosticItem
