import re
import os
import snakewrap.exception as exception

class TemplateFile():
    def __init__(self, name):
        self.name = name
        self.raw_name = None

    def get_wildcard(self, key):
        return self.wildcards[key]

    def assign(self, name):
        self.filename = name
        self.matched = True
    
    def infer_raw_name(self):
        if not self.matched:
            raise exception.FileUnassignedException('Please run `assign` method before inferring raw name.')

        if self.raw_name is None:
            return self.filename
        else:
            return self.raw_name.format(**self.wildcards)

    def rename_command(self):
        if not self.matched:
            raise exception.FileUnassignedException('Please run `assign` method before inferring raw name.')
        
        if self.raw_name is not None:
            return 'mv %s %s' % (self.infer_raw_name(), self.filename)
        else:
            # We don't need renaming when raw_name is not specified.
            return None

class SimpleTemplateFile(TemplateFile):
    def __init__(self, name):
        super(SimpleTemplateFile, self).__init__(name)

class RenamedTemplateFile(TemplateFile):
    def __init__(self, name, regex, raw_name):
        super(RenamedTemplateFile, self).__init__(name)

        self.regex, self.compiled_regex = regex, re.compile(regex)
        self.raw_name = raw_name
        self.matched = False
        self.wildcards = None

    def assign(self, name):
        found = re.match(self.regex, name)
        if found and not self.matched:
            self.matched = True
            self.filename, self.wildcards = name, found.groupdict()
        
        if not self.matched:
            raise exception.FileNameMismatchException()

class SimpleTemplateDirectory(TemplateFile):
    def __init__(self, name):
        super(SimpleTemplateDirectory, self).__init__(name)

    def check_subdirectory(self, dir):
        subdirs = os.listdir(self.filename)
        if dir in subdirs and os.path.isdir(os.path.join(self.filename, dir)):
            return True
        return False
        
