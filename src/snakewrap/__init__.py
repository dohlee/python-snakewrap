__version__ = '0.0.0'

from snakewrap.files import SimpleTemplateFile, RenamedTemplateFile, SimpleTemplateDirectory
from snakewrap.ruleinput import SimpleRuleInput, MultiRuleInput, ListInput

__all__ = [
    'files',
    'ruleinput',
    'ruleoutput',
    'rulethreads',
    'ruleparams',
    'wrapper',
]