<h1 align="center">snakewrap</h1>
<p align="center">Write snakemake wrappers in more organized way.</p>
<p align="center">

</p>

<h2 align="center">Quickstart</h2>

## Before...

```python
# wrapper.py for samtools sort.
from os import path
from snakemake.shell import shell

# Extract log.
log = snakemake.log_fmt_shell(stdout=False, stderr=True)

# Extract parameters.
extra = snakemake.params.get('extra', '')

# Extract required arguments.
bam = snakemake.input.bam
sorted_bam = snakemake.output.sorted_bam
prefix = path.splitext(sorted_bam)[0]

# Execute shell command.
shell(
    "("
    "samtools sort "
    "{extra} "
    "-@ {snakemake.threads} "
    "-o {sorted_bam} "
    "-T {prefix} "
    "{bam}"
    ") "
    "{log}"
)
```

## After...

```python
import snakewrap as sw

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
wrapper.run()
```
