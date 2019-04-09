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
    
    read = files.File('read', r'(?P<sample>.+?).fastq$', raw_name=None)
    mock_snakemake_input = {'reads': [read]}
    simple_read_input.attach(mock_snakemake_input)

    assert simple_read_input.names == ['read']

def test_simpleruleinput_attach_failed_because_of_unknown_key():
    simple_read_input = ruleinput.SimpleRuleInput('reads', '-i')

    read = files.File('read', r'(?P<sample>.+?).fastq$', raw_name=None)
    mock_snakemake_input = {'read': [read]}
    
    with pytest.raises(exception.RuleInputException):
        simple_read_input.attach(mock_snakemake_input)

def test_simpleruleinput_attach_failed_because_of_more_than_one_input():
    simple_read_input = ruleinput.SimpleRuleInput('reads', '-i')

    read1 = files.File('read1', r'(?P<sample>.+?).read1.fastq$', raw_name=None)
    read2 = files.File('read2', r'(?P<sample>.+?).read2.fastq$', raw_name=None)
    mock_snakemake_input = {'reads': [read1, read2]}

    with pytest.raises(exception.RuleInputException):
        simple_read_input.attach(mock_snakemake_input)

def test_multiruleinput_class_init():
    paired_end_input = ruleinput.MultiRuleInput('reads', '-1', '-2')
    assert paired_end_input.rule_key == 'reads'
    assert paired_end_input.command_keys == ['-1', '-2']

def test_multiruleinput_attach():
    paired_end_input = ruleinput.MultiRuleInput('reads', '-1', '-2')

    read1 = files.File('read1', r'(?P<sample>.+?).read1.fastq$', raw_name=None)
    read2 = files.File('read2', r'(?P<sample>.+?).read2.fastq$', raw_name=None)
    mock_snakemake_input = {'reads': [read1, read2]}
    paired_end_input.attach(mock_snakemake_input)

    assert paired_end_input.names == ['read1', 'read2']

def test_simpleruleinput_access_file():
    si = ruleinput.SimpleRuleInput('reads', '-i')

    read = files.File('read', r'(?P<sample>.+?).fastq$', raw_name=None)
    mock_snakemake_input = {'reads': [read]}
    si.attach(mock_snakemake_input)

    assert si.files[0].name == 'read'

def test_simpleruleinput_match():
    si = ruleinput.SimpleRuleInput('read', '-i')

    read = files.File('read', r'(?P<sample>.+?).fastq', raw_name=None)
    si.attach({'read': [read]})

    fs = {'read': ['mysample.fastq']}
    si.match(fs)

    assert si.files[0].wildcards['sample'] == 'mysample'

def test_simpleruleinput_match_failed_because_of_not_attached_ruleinput():
    si = ruleinput.SimpleRuleInput('read', '-i')

    read = files.File('read', r'(?P<sample>.+?).fastq', raw_name=None)
    fs = {'read': ['mysample.fastq']}

    with pytest.raises(exception.RuleNotAttachedException):
        si.match(fs)

def test_simpleruleinput_string():
    si = ruleinput.SimpleRuleInput('read', '-i')

    read = files.File('read', r'(?P<sample>.+?).fastq', raw_name=None)

    mock_templates = {
        'read': [read],
    }
    si.attach(mock_templates)

    mock_input = {
        'read': ['mysample.fastq'],
    }
    si.match(mock_input)

    assert str(si) == '-i mysample.fastq'

def test_multiruleinput_workflow():
    si = ruleinput.MultiRuleInput('reads', '-1', '-2')

    read1 = files.File('read1', r'(?P<sample>.+?).read1.fastq', raw_name=None)
    read2 = files.File('read2', r'(?P<sample>.+?).read2.fastq', raw_name=None)

    mock_templates = {
        'reads': [read1, read2],
    }
    si.attach(mock_templates)

    mock_input = {
        'reads': ['mysample.read1.fastq', 'mysample.read2.fastq'],
    }
    si.match(mock_input)

    assert str(si) == '-1 mysample.read1.fastq -2 mysample.read2.fastq'