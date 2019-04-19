import snakewrap as sw

def test_samtools():
    snakemake = {
        'input': {
            'bam': 'test.bam',
        },
        'output': {
            'sorted_bam': 'test.sorted.bam',
        },
        'threads': 1,
        'params': {
            'extra': '',
        },
        'log': 'logs/log.txt'
    }

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

    expected_shell_command = '(samtools sort test.bam -o test.sorted.bam -T test.sorted -@ 1) 2> logs/log.txt'
    assert wrapper.shell_command() == expected_shell_command
    wrapper.run()