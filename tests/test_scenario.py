import snakewrap as sw

def test_samtools():
    bam = sw.SimpleTemplateFile('bam')
    bam_output = sw.SimpleTemplateFile('sorted_bam')

    # Create RuleInput object.
    input = sw.SimpleRuleInput(
        rule_keys='bam',
    )
    # Set template to the object.
    input.set_template({
        'bam': (bam, None),
    })

    # Create RuleOutput object.
    output = sw.SimpleRuleOutput(
        rule_keys='sorted_bam',
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

    wrapper = Wrapper(snakemake, command='samtools sort')
    wrapper.run()