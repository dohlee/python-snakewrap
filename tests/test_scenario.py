import snakewrap as sw
from snakemake.workflow import Workflow
from snakemake.script import Snakemake
from os.path import join
import os


DATA_PATH = os.path.realpath(join(os.path.dirname(__file__), 'data'))

def get_snakemake_object(snakefile, rule_name):
    w = Workflow(os.path.join(DATA_PATH, snakefile))
    w.include(os.path.join(DATA_PATH, snakefile))
    w.check()
    r = w._rules[rule_name]
    snakemake = Snakemake(
        input=r._input,
        output=r._output,
        params=r._params,
        wildcards=r._wildcard_names,
        threads=r.resources['_cores'],
        resources=r.resources,
        log=r.log,
        config=None,
        rulename=r.name,
        bench_iteration=1,
    )
    print(r.wildcard_constraints)

    return snakemake

def test_samtools_sort():
    snakemake = get_snakemake_object('samtools_sort_snakefile', 'samtools_sort')

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
        prefix=(lambda input, output: os.path.splitext(output['sorted_bam'])[0], '-T')
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