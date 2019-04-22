import snakewrap as sw
from os.path import join
import os

class NamedList():
    def __init__(self, d):
        self.d = d
        for k, v in d.items():
            setattr(self, k, v)
    
    def items(self):
        return self.d.items()

class InputFiles(NamedList):
    pass

class OutputFiles(NamedList):
    pass

class Params(NamedList):
    pass

class Wildcards(NamedList):
    pass

class MockSnakemake():
    __slots__ = ['input', 'output', 'params', 'wildcards', 'log', 'threads']

    def __init__(self, rule_dict):
        setattr(self, 'input', InputFiles(rule_dict.get('input', {})))
        setattr(self, 'output', OutputFiles(rule_dict.get('output', {})))
        setattr(self, 'params', Params(rule_dict.get('params', {'extra': ''})))
        setattr(self, 'wildcards', Wildcards(rule_dict.get('wildcards', {})))
        setattr(self, 'log', rule_dict.get('log', ''))
        setattr(self, 'threads', rule_dict.get('threads', 1))

    def log_fmt_shell(self, stdout=True, stderr=True, append=False):
        """Taken from https://bitbucket.org/snakemake/snakemake."""
        if not self.log:
            return ""
        lookup = {
            (True, True, True): " >> {0} 2>&1",
            (True, False, True): " >> {0}",
            (False, True, True): " 2>> {0}",
            (True, True, False): " > {0} 2>&1",
            (True, False, False): " > {0}",
            (False, True, False): " 2> {0}",
        }
        return lookup[(stdout, stderr, append)].format(self.log)

def test_samtools_sort():
    snakemake = MockSnakemake({
        'input': {'bam': 'test.bam'},
        'output': {'sorted_bam': 'test.sorted.bam'},
        'log': 'logs/samtools_sort/test.log',
    })

    bam = sw.SimpleTemplateFile('bam')
    sorted_bam = sw.SimpleTemplateFile('sorted_bam')

    # Create RuleInput object.
    input = sw.SimpleRuleInput({
        'bam': (bam, None),
    })

    # Create RuleOutput object.
    output = sw.SimpleRuleOutput({
        'sorted_bam': [(sorted_bam, '-o')],
    })

    # Create RuleParams object.
    params = sw.SimpleRuleParams(
        extra=True,
        prefix=(lambda input, output: os.path.splitext(output.sorted_bam)[0], '-T')
    )

    # Create RuleThreads object.
    threads = sw.SimpleRuleThreads(command_key='-@')

    wrapper = sw.Wrapper(
        snakemake,
        command='samtools sort',
        input=input,
        output=output,
        params=params,
        threads=threads,
    )

    expected_shell_command = '( samtools sort test.bam -o test.sorted.bam -T test.sorted -@ 1 ) 2> logs/samtools_sort/test.log'
    assert wrapper.shell_command() == expected_shell_command
    wrapper.run()