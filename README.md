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
bam_output = sw.SimpleTemplateFile('sorted_bam')

# Create RuleInput object.
input = sw.SimpleRuleInput(
    rule_key='bam',
)
# Set template to the object.
input.set_template({
    'bam': [(bam, None)],
})

# Create RuleOutput object.
output = sw.SimpleRuleOutput(
    rule_key='sorted_bam',
)
# Set template to the object.
output.set_template({
    'sorted_bam': [(sorted_bam, '-o')],
})

# Create RuleParams object.
params = sw.SimpleRuleParams(
    extra=True,
)

# Create RuleThreads object.
threads = sw.SimpleRuleThreads()

wrapper = Wrapper(snakemake)
wrapper.run()

```
