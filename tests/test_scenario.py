import snakewrap as sw
from snakewrap.parameter import Parameter
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
        setattr(self, 'input', InputFiles(rule_dict.get('input', dict())))
        setattr(self, 'output', OutputFiles(rule_dict.get('output', dict())))
        setattr(self, 'params', Params(rule_dict.get('params', {'extra': ''})))
        setattr(self, 'wildcards', Wildcards(rule_dict.get('wildcards', dict())))
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
        'bam': Parameter(f=bam, option=None),
    })

    # Create RuleOutput object.
    output = sw.SimpleRuleOutput({
        'sorted_bam': [Parameter(f=sorted_bam, option='-o')],
    })

    # Create RuleParams object.
    params = sw.SimpleRuleParams(
        extra=True,
        prefix=Parameter(f=lambda sn: os.path.splitext(sn.output.sorted_bam)[0], option='-T')
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


def test_samtools_sort_with_extra():
    snakemake = MockSnakemake({
        'input': {'bam': 'test.bam'},
        'output': {'sorted_bam': 'test.sorted.bam'},
        'params': {'extra': '-l 9'},
        'log': 'logs/samtools_sort/test.log',
    })

    bam = sw.SimpleTemplateFile('bam')
    sorted_bam = sw.SimpleTemplateFile('sorted_bam')

    # Create RuleInput object.
    input = sw.SimpleRuleInput({
        'bam': Parameter(f=bam, option=None),
    })

    # Create RuleOutput object.
    output = sw.SimpleRuleOutput({
        'sorted_bam': [Parameter(f=sorted_bam, option='-o')],
    })

    # Create RuleParams object.
    params = sw.SimpleRuleParams(
        extra=True,
        prefix=Parameter(f=lambda sn: os.path.splitext(sn.output.sorted_bam)[0], option='-T')
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

    expected_shell_command = '( samtools sort test.bam -o test.sorted.bam -l 9 -T test.sorted -@ 1 ) 2> logs/samtools_sort/test.log'
    assert wrapper.shell_command() == expected_shell_command
    wrapper.run()


def test_bismark_single_with_unzipped_fastq():
    snakemake = MockSnakemake({
        'input': {
            'fastq': 'data/test.fastq',
            'reference_dir': 'reference/hg38_bismark',
            'bisulfite_genome_dir': 'reference/hg38_bismark/Bisulfite_Genome',
        },
        'output': {
            'bam': 'result/test.bismark.bam',
            'report': 'result/test.bismark_report.txt',
        },
        'wildcards': {
            'sample': 'test',
        },
        'threads': 6,
        'params': {'extra': ''},
        'log': 'logs/bismark/test.log',
    })

    # Define files.
    fastq = sw.SimpleTemplateFile('fastq')
    reference_dir = sw.SimpleTemplateDirectory('reference_dir')
    bisulfite_genome_dir = sw.SimpleTemplateDirectory('bisulfite_genome_dir')
    bam = sw.RenamedTemplateFile(
        name='bam',
        regex=r'(?P<prefix>.+?).bismark.bam',
        raw_name='{prefix}_bismark_bt2.bam',
    )
    report = sw.RenamedTemplateFile(
        name='report',
        regex=r'(?P<prefix>.+?).bismark_report.txt',
        raw_name='{prefix}_bismark_bt2_SE_report.txt'
    )

    # Define input, output, parameters, threads.
    input = sw.SimpleRuleInput({
        'reference_dir': Parameter(f=reference_dir, option=None),
        'fastq': Parameter(f=fastq, option=None),
        'bisulfite_genome_dir': Parameter(f=bisulfite_genome_dir, option=None, used=False)
    })

    output = sw.SimpleRuleOutput({
        'bam': Parameter(f=bam, option=None, used=False),
        'report': Parameter(f=report, option=None, used=False),
    })

    params = sw.SimpleRuleParams(
        extra=True,
        outdir=Parameter(f=lambda sn: os.path.dirname(sn.output.bam), option='-o', used=True)
    )

    threads = sw.ScaledRuleThreads(
        command_key='--multicore',
        scale=lambda sn: 1/2 if '--bowtie1' in sn.params.extra else 1/3
    )

    wrapper = sw.Wrapper(
        snakemake,
        command='bismark',
        input=input,
        output=output,
        params=params,
        threads=threads,
    )

    assert bam.rename_command() == 'mv result/test_bismark_bt2.bam result/test.bismark.bam'
    assert report.rename_command() == 'mv result/test_bismark_bt2_SE_report.txt result/test.bismark_report.txt'

    assert wrapper.shell_command() == \
        '( bismark reference/hg38_bismark data/test.fastq '\
        '-o result ' \
        '--multicore 2 ' \
        '&& mv result/test_bismark_bt2.bam result/test.bismark.bam ' \
        '&& mv result/test_bismark_bt2_SE_report.txt result/test.bismark_report.txt ' \
        ') 2> logs/bismark/test.log' \



def test_bismark_single_with_unzipped_fastq_bowtie1():
    snakemake = MockSnakemake({
        'input': {
            'fastq': 'data/test.fastq',
            'reference_dir': 'reference/hg38_bismark',
            'bisulfite_genome_dir': 'reference/hg38_bismark/Bisulfite_Genome',
        },
        'output': {
            'bam': 'result/test.bismark.bam',
            'report': 'result/test.bismark_report.txt',
        },
        'wildcards': {
            'sample': 'test',
        },
        'threads': 6,
        'params': {'extra': '--bowtie1'},
        'log': 'logs/bismark/test.log',
    })

    # Define files.
    fastq = sw.SimpleTemplateFile('fastq')
    reference_dir = sw.SimpleTemplateDirectory('reference_dir')
    bisulfite_genome_dir = sw.SimpleTemplateDirectory('bisulfite_genome_dir')
    bam = sw.RenamedTemplateFile(
        name='bam',
        regex=r'(?P<prefix>.+?).bismark.bam',
        raw_name='{prefix}_bismark_bt2.bam',
    )
    report = sw.RenamedTemplateFile(
        name='report',
        regex=r'(?P<prefix>.+?).bismark_report.txt',
        raw_name='{prefix}_bismark_bt2_SE_report.txt'
    )

    # Define input, output, parameters, threads.
    input = sw.SimpleRuleInput({
        'reference_dir': Parameter(f=reference_dir, option=None),
        'fastq': Parameter(f=fastq, option=None),
        'bisulfite_genome_dir': Parameter(f=bisulfite_genome_dir, option=None, used=False),
    })

    output = sw.SimpleRuleOutput({
        'bam': Parameter(f=bam, option=None, used=False),
        'report': Parameter(f=report, option=None, used=False),
    })

    params = sw.SimpleRuleParams(
        extra=True,
        outdir=Parameter(lambda sn: os.path.dirname(sn.output.bam), '-o', True)
    )

    threads = sw.ScaledRuleThreads(
        command_key='--multicore',
        scale=lambda sn: 1/2 if '--bowtie1' in sn.params.extra else 1/3
    )

    wrapper = sw.Wrapper(
        snakemake,
        command='bismark',
        input=input,
        output=output,
        params=params,
        threads=threads,
    )

    assert bam.rename_command() == 'mv result/test_bismark_bt2.bam result/test.bismark.bam'
    assert report.rename_command() == 'mv result/test_bismark_bt2_SE_report.txt result/test.bismark_report.txt'

    assert wrapper.shell_command() == \
        '( bismark reference/hg38_bismark data/test.fastq '\
        '--bowtie1 ' \
        '-o result ' \
        '--multicore 3 ' \
        '&& mv result/test_bismark_bt2.bam result/test.bismark.bam ' \
        '&& mv result/test_bismark_bt2_SE_report.txt result/test.bismark_report.txt ' \
        ') 2> logs/bismark/test.log' \



def test_bismark_paired_with_unzipped_fastq():
    snakemake = MockSnakemake({
        'input': {
            'fastq': ['data/test.read1.fastq', 'data/test.read2.fastq'],
            'reference_dir': 'reference/hg38_bismark',
            'bisulfite_genome_dir': 'reference/hg38_bismark/Bisulfite_Genome',
        },
        'output': {
            'bam': 'result/test.bismark.bam',
            'report': 'result/test.bismark_report.txt',
        },
        'wildcards': {
            'sample': 'test',
        },
        'threads': 6,
        'params': {'extra': '--bowtie1'},
        'log': 'logs/bismark/test.log',
    })

    # Define files.
    read1 = sw.SimpleTemplateFile('fastq')
    read2 = sw.SimpleTemplateFile('fastq')
    reference_dir = sw.SimpleTemplateDirectory('reference_dir')
    bisulfite_genome_dir = sw.SimpleTemplateDirectory('bisulfite_genome_dir')
    bam = sw.RenamedTemplateFile(
        name='bam',
        regex=r'(?P<prefix>.+?).bismark.bam',
        raw_name='{prefix}.read1_bismark_bt2_pe.bam',
    )
    report = sw.RenamedTemplateFile(
        name='report',
        regex=r'(?P<prefix>.+?).bismark_report.txt',
        raw_name='{prefix}.read1_bismark_bt2_PE_report.txt'
    )

    # Define input, output, parameters, threads.
    input = sw.SimpleRuleInput({
        'reference_dir': Parameter(f=reference_dir, option=None),
        'fastq': [Parameter(f=read1, option='-1'), Parameter(f=read2, option='-2')],
        'bisulfite_genome_dir': Parameter(f=bisulfite_genome_dir, option=None, used=False),
    })

    output = sw.SimpleRuleOutput({
        'bam': Parameter(f=bam, option=None, used=False),
        'report': Parameter(f=report, option=None, used=False),
    })

    params = sw.SimpleRuleParams(
        extra=True,
        outdir=Parameter(f=lambda sn: os.path.dirname(sn.output.bam), option='-o')
    )

    threads = sw.ScaledRuleThreads(
        command_key='--multicore',
        scale=lambda sn: 1/2 if '--bowtie1' in sn.params.extra else 1/3
    )

    wrapper = sw.Wrapper(
        snakemake,
        command='bismark',
        input=input,
        output=output,
        params=params,
        threads=threads,
    )

    assert wrapper.shell_command() == \
        '( bismark reference/hg38_bismark ' \
        '-1 data/test.read1.fastq -2 data/test.read2.fastq '\
        '--bowtie1 ' \
        '-o result ' \
        '--multicore 3 ' \
        '&& mv result/test.read1_bismark_bt2_pe.bam result/test.bismark.bam ' \
        '&& mv result/test.read1_bismark_bt2_PE_report.txt result/test.bismark_report.txt ' \
        ') 2> logs/bismark/test.log' \


def test_bismark_paired_with_zipped_fastq():
    snakemake = MockSnakemake({
        'input': {
            'fastq': ['data/test.read1.fastq.gz', 'data/test.read2.fastq.gz'],
            'reference_dir': 'reference/hg38_bismark',
            'bisulfite_genome_dir': 'reference/hg38_bismark/Bisulfite_Genome',
        },
        'output': {
            'bam': 'result/test.bismark.bam',
            'report': 'result/test.bismark_report.txt',
        },
        'wildcards': {
            'sample': 'test',
        },
        'threads': 6,
        'params': {'extra': '--bowtie1'},
        'log': 'logs/bismark/test.log',
    })

    # Define files.
    read1 = sw.SimpleTemplateFile('fastq')
    read2 = sw.SimpleTemplateFile('fastq')
    reference_dir = sw.SimpleTemplateDirectory('reference_dir')
    bisulfite_genome_dir = sw.SimpleTemplateDirectory('bisulfite_genome_dir')
    bam = sw.RenamedTemplateFile(
        name='bam',
        regex=r'(?P<prefix>.+?).bismark.bam',
        raw_name='{prefix}.read1_bismark_bt2_pe.bam',
    )
    report = sw.RenamedTemplateFile(
        name='report',
        regex=r'(?P<prefix>.+?).bismark_report.txt',
        raw_name='{prefix}.read1_bismark_bt2_PE_report.txt'
    )

    input = sw.SimpleRuleInput({
        'reference_dir': Parameter(reference_dir, None),
        'fastq': [Parameter(read1, '-1'), Parameter(read2, '-2')],
        'bisulfite_genome_dir': Parameter(bisulfite_genome_dir, None, False)
    })

    output = sw.SimpleRuleOutput({
        'bam': Parameter(bam, None, False),
        'report': Parameter(report, None, False),
    })

    params = sw.SimpleRuleParams(
        extra=True,
        outdir=Parameter(lambda sn: os.path.dirname(sn.output.bam), '-o', True)
    )

    threads = sw.ScaledRuleThreads(
        command_key='--multicore',
        scale=lambda sn: 1/2 if '--bowtie1' in sn.params.extra else 1/3
    )

    wrapper = sw.Wrapper(
        snakemake,
        command='bismark',
        input=input,
        output=output,
        params=params,
        threads=threads,
    )

    assert wrapper.shell_command() == \
        '( bismark reference/hg38_bismark ' \
        '-1 data/test.read1.fastq.gz -2 data/test.read2.fastq.gz '\
        '--bowtie1 ' \
        '-o result ' \
        '--multicore 3 ' \
        '&& mv result/test.read1_bismark_bt2_pe.bam result/test.bismark.bam ' \
        '&& mv result/test.read1_bismark_bt2_PE_report.txt result/test.bismark_report.txt ' \
        ') 2> logs/bismark/test.log' \
    

def test_bismark_genome_preparation():
    # Define mock snakemake object.
    snakemake = MockSnakemake({
        'input': {
            'reference_dir': 'reference/hg38',
        },
        'output': {
            'bisulfite_genome_dir': 'reference/hg38/Bisulfite_Genome'
        },
        'params': {
            'extra': '',
        },
        'log': 'logs/bismark_genome_preparation/log.log'
    })

    # Define files.
    reference_dir = sw.SimpleTemplateDirectory('reference_dir')
    bisulfite_genome_dir = sw.SimpleTemplateDirectory('bisulfite_genome_dir')

    # Define input, output, and parameters.
    input = sw.SimpleRuleInput({
        'reference_dir': Parameter(reference_dir, None),
    })

    output = sw.SimpleRuleOutput({
        'bisulfite_genome_dir': Parameter(bisulfite_genome_dir, None, False),
    })

    params = sw.SimpleRuleParams(extra=True)

    wrapper = sw.Wrapper(
        snakemake=snakemake,
        command='bismark_genome_preparation',
        input=input,
        output=output,
        params=params,
    )
    assert wrapper.shell_command() == '( bismark_genome_preparation reference/hg38 ) 2> logs/bismark_genome_preparation/log.log'

def test_bwa_mem_single():
    # Define mock snakemake object.
    snakemake = MockSnakemake({
        'input': {
            'reference': 'reference/hg38/hg38.fasta',
            'reads': ['test.fastq.gz'],
        },
        'output': {
            'bam': 'test.bam',
        },
        'params': {
            'extra': '',
        },
        'log': 'logs/bwa_mem/test.log',
        'threads': 12,
    })

    expected_command = '( bwa mem reference/hg38/hg38 test.fastq.gz -t 12 | samtools view -Sb - > test.bam ) 2> logs/bwa_mem/test.log'

    # Define files.
    reference = sw.SimpleTemplateFile('reference')
    read = sw.SimpleTemplateFile('reads')
    bam = sw.SimpleTemplateFile('bam')

    # Define input, output, parameters and threads.
    input = sw.SimpleRuleInput({
        'reference': Parameter(f=reference, option=None, used=False),
        'reads': [Parameter(f=read, option=None)],
    })

    output = sw.SimpleRuleOutput({
        'bam': Parameter(f=bam, option=None, used=True, redirected=True),
    })

    threads = sw.SimpleRuleThreads(command_key='-t')

    params = sw.SimpleRuleParams(
        extra=True,
        db_prefix=Parameter(f=lambda sn: os.path.splitext(sn.input.reference)[0], option=None, priority=1),
        pipe_command=Parameter(f=lambda sn: '| samtools view -Sb -' if sn.output.bam.endswith('.bam') else '', option=None, priority=100)
    )

    wrapper = sw.Wrapper(
        snakemake=snakemake,
        command='bwa mem',
        input=input,
        output=output,
        threads=threads,
        params=params,
    )

    assert wrapper.shell_command() == expected_command