import re
import os
import snakewrap.exception as exception

class File():
    def __init__(self, name, regex, raw_name=None):
        self.name = name
        self.regex, self.compiled_regex = regex, re.compile(regex)
        self.raw_name = raw_name
        self.matched = False
        self.wildcards = None
    
    def get_wildcard(self, key):
        return self.wildcards[key]

    def match(self, names):
        for name in names:
            found = re.match(self.regex, name)
            if found and not self.matched:
                self.matched = True
                self.filename, self.wildcards = name, found.groupdict()
            elif found:
                raise exception.AmbiguousFileNameException()
        
        if not self.matched:
            raise exception.NoMatchingFileNameException()
    
    def infer_raw_name(self):
        if not self.matched:
            raise exception.FileUnmatchedException('Please run `match` method before inferring raw name.')

        if self.raw_name is None:
            return self.filename
        else:
            return self.raw_name.format(**self.wildcards)

    def rename_command(self):
        if not self.matched:
            raise exception.FileUnmatchedException('Please run `match` method before inferring raw name.')
        
        if self.raw_name is not None:
            return 'mv %s %s' % (self.infer_raw_name(), self.filename)
        else:
            # We don't need renaming when raw_name is not specified.
            return None

class Directory(File):
    def __init__(self, name, regex, raw_name=None):
        super(Directory, self).__init__(name, regex, raw_name)

    def check_subdirectory(self, dir):
        subdirs = os.listdir(self.filename)
        if dir in subdirs and os.path.isdir(os.path.join(self.filename, dir)):
            return True
        return False
        
