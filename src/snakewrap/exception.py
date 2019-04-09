class RuleInputException(Exception):
    pass

class RuleOutputException(Exception):
    pass

class RuleParameterException(Exception):
    pass

class RuleNotAttachedException(Exception):
    pass

class NoMatchingFileNameException(Exception):
    pass

class AmbiguousFileNameException(Exception):
    pass

class FileUnmatchedException(Exception):
    pass