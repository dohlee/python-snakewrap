__version__ = '0.0.0'

from snakewrap.files import SimpleTemplateFile, RenamedTemplateFile, SimpleTemplateDirectory
from snakewrap.ruleinput import SimpleRuleInput, MultiRuleInput, ListInput
from snakewrap.ruleoutput import SimpleRuleOutput
from snakewrap.ruleparams import SimpleRuleParams
from snakewrap.rulethreads import SimpleRuleThreads, ScaledRuleThreads
from snakewrap.wrapper import Wrapper

__all__ = [
    'files',
    'ruleinput',
    'ruleoutput',
    'rulethreads',
    'ruleparams',
    'wrapper',
]