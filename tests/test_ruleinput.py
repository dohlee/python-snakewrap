import pytest
from snakewrap import ruleinput, exception, files

def test_base_ruleinput_class_init():
    read_input = ruleinput.RuleInput('reads', '-i')
    assert read_input.rule_key == 'reads'
    assert read_input.command_keys == ['-i']

def test_base_ruleinput_class_init_with_name_and_desc():
    read_input = ruleinput.RuleInput('reads', '-i', name='read_input', desc='Input read.')
    assert read_input.rule_key == 'reads'
    assert read_input.command_keys == ['-i']
    assert read_input.name == 'read_input'
    assert read_input.desc == 'Input read.'

def test_base_ruleinput_class_init_with_multiple_command_keys():
    pairedend_read_input = ruleinput.RuleInput('reads', '-1', '-2')
    assert pairedend_read_input.rule_key == 'reads'
    assert pairedend_read_input.command_keys == ['-1', '-2']

    pairedend_read_input = ruleinput.RuleInput('reads', '-a', '-b')
    assert pairedend_read_input.rule_key == 'reads'
    assert pairedend_read_input.command_keys == ['-a', '-b']

def test_simpleruleinput_class_init():
    simple_read_input = ruleinput.SimpleRuleInput('read', '-1')
    assert simple_read_input.rule_key == 'read'
    assert simple_read_input.command_keys == ['-1']

def test_simpleruleinput_class_init_failed():
    with pytest.raises(exception.RuleInputException):
        ruleinput.SimpleRuleInput('reads', '-1', '-2')

def test_simpleruleinput_attach():
    simple_read_input = ruleinput.SimpleRuleInput('reads', '-1')
    
    read = files.File('read.fastq', r'(?P<sample>.+?).fastq$', raw_name=None)
    mock_snakemake_input = {'reads': [read]}
    simple_read_input.attach(mock_snakemake_input)

    assert simple_read_input.names == ['read.fastq']

def test_simpleruleinput_attach_failed_because_of_unknown_key():
    simple_read_input = ruleinput.SimpleRuleInput('reads', '-i')

    read = files.File('read.fastq', r'(?P<sample>.+?).fastq$', raw_name=None)
    mock_snakemake_input = {'read': [read]}
    
    with pytest.raises(exception.RuleInputException):
        simple_read_input.attach(mock_snakemake_input)

def test_simpleruleinput_attach_failed_because_of_more_than_one_input():
    simple_read_input = ruleinput.SimpleRuleInput('reads', '-i')

    read1 = files.File('read.read1.fastq', r'(?P<sample>.+?).read1.fastq$', raw_name=None)
    read2 = files.File('read.read2.fastq', r'(?P<sample>.+?).read2.fastq$', raw_name=None)
    mock_snakemake_input = {'reads': [read1, read2]}

    with pytest.raises(exception.RuleInputException):
        simple_read_input.attach(mock_snakemake_input)

def test_simplerule_input_string():
    simple_read_input = ruleinput.SimpleRuleInput('reads', '-i')

    read = files.File('read.fastq', r'(?P<sample>.+?).fastq$', raw_name=None)
    mock_snakemake_input = {'reads': [read]}
    simple_read_input.attach(mock_snakemake_input)

    assert str(simple_read_input) == '-i read.fastq'

def test_multiruleinput_class_init():
    paired_end_input = ruleinput.MultiRuleInput('reads', '-1', '-2')
    assert paired_end_input.rule_key == 'reads'
    assert paired_end_input.command_keys == ['-1', '-2']

def test_multiruleinput_attach():
    paired_end_input = ruleinput.MultiRuleInput('reads', '-1', '-2')

    read1 = files.File('sample.read1.fastq', r'(?P<sample>.+?).read1.fastq$', raw_name=None)
    read2 = files.File('sample.read2.fastq', r'(?P<sample>.+?).read2.fastq$', raw_name=None)
    mock_snakemake_input = {'reads': [read1, read2]}
    paired_end_input.attach(mock_snakemake_input)

    assert paired_end_input.names == ['sample.read1.fastq', 'sample.read2.fastq']