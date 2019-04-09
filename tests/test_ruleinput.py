import pytest
from snakewrap import ruleinput, exception
from snakewrap.files import SimpleTemplateFile, RenamedTemplateFile, SimpleTemplateDirectory

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

def test_simpleruleinput_set_template():
    simple_read_input = ruleinput.SimpleRuleInput('reads', '-1')
    
    read = SimpleTemplateFile('read')
    mock_snakemake_input = {'reads': [read]}
    simple_read_input.set_template(mock_snakemake_input)

    assert simple_read_input.names == ['read']

def test_simpleruleinput_set_template_failed_because_of_unknown_key():
    simple_read_input = ruleinput.SimpleRuleInput('reads', '-i')

    read = SimpleTemplateFile('read')
    mock_snakemake_input = {'read': [read]}
    
    with pytest.raises(exception.RuleInputException):
        simple_read_input.set_template(mock_snakemake_input)

def test_simpleruleinput_set_template_failed_because_of_more_than_one_input():
    simple_read_input = ruleinput.SimpleRuleInput('reads', '-i')

    read1 = SimpleTemplateFile('read1')
    read2 = SimpleTemplateFile('read2')
    mock_snakemake_input = {'reads': [read1, read2]}

    with pytest.raises(exception.RuleInputException):
        simple_read_input.set_template(mock_snakemake_input)

def test_multiruleinput_class_init():
    paired_end_input = ruleinput.MultiRuleInput('reads', '-1', '-2')
    assert paired_end_input.rule_key == 'reads'
    assert paired_end_input.command_keys == ['-1', '-2']

def test_multiruleinput_set_template():
    paired_end_input = ruleinput.MultiRuleInput('reads', '-1', '-2')

    read1 = SimpleTemplateFile('read1')
    read2 = SimpleTemplateFile('read2')
    mock_snakemake_input = {'reads': [read1, read2]}
    paired_end_input.set_template(mock_snakemake_input)

    assert paired_end_input.names == ['read1', 'read2']

def test_simpleruleinput_access_file():
    si = ruleinput.SimpleRuleInput('reads', '-i')

    read = SimpleTemplateFile('read')
    mock_snakemake_input = {'reads': [read]}
    si.set_template(mock_snakemake_input)

    assert si.files[0].name == 'read'

def test_simpleruleinput_assign():
    si = ruleinput.SimpleRuleInput('read', '-i')

    read = SimpleTemplateFile('read')
    si.set_template({'read': [read]})

    fs = {'read': ['mysample.fastq']}
    si.assign(fs)

    assert si.files[0].filename == 'mysample.fastq'

def test_simpleruleinput_assign_failed_because_of_not_set_template():
    si = ruleinput.SimpleRuleInput('read', '-i')

    read = SimpleTemplateFile('read')
    fs = {'read': ['mysample.fastq']}

    with pytest.raises(exception.TemplateNotSetException):
        si.assign(fs)

def test_simpleruleinput_string():
    si = ruleinput.SimpleRuleInput('read', '-i')

    read = SimpleTemplateFile('read')

    mock_templates = {
        'read': [read],
    }
    si.set_template(mock_templates)

    mock_input = {
        'read': ['mysample.fastq'],
    }
    si.assign(mock_input)

    assert str(si) == '-i mysample.fastq'

def test_multiruleinput_simpletemplatefile_workflow():
    si = ruleinput.MultiRuleInput('reads', '-1', '-2')

    read1 = SimpleTemplateFile('read1')
    read2 = SimpleTemplateFile('read2')

    mock_templates = {
        'reads': [read1, read2],
    }
    si.set_template(mock_templates)

    mock_input = {
        'reads': ['mysample.read1.fastq', 'mysample.read2.fastq'],
    }
    si.assign(mock_input, sequential=True)

    assert str(si) == '-1 mysample.read1.fastq -2 mysample.read2.fastq'

def test_listinput_init():
    li = ruleinput.ListInput('reads', '-i')

def test_listinput_init_failed_because_of_multiple_command_keys():
    with pytest.raises(AttributeError):
        li = ruleinput.ListInput('reads', '-1', '-2')

def test_listinput_set_template():
    li = ruleinput.ListInput('reads', '-i')
    read1 = SimpleTemplateFile('read1')
    read2 = SimpleTemplateFile('read2')

    mock_templates = {
        'reads': [read1, read2],
    }
    li.set_template(mock_templates)

def test_listinput_assign():
    li = ruleinput.ListInput('reads', '-i')
    read1 = SimpleTemplateFile('read1')
    read2 = SimpleTemplateFile('read2')

    mock_templates = {
        'reads': [read1, read2],
    }
    li.set_template(mock_templates)

    mock_input = {
        'reads': ['mysample.read1.fastq', 'mysample.read2.fastq'],
    }
    li.assign(mock_input, sequential=True)

    assert li.files[0].filename == 'mysample.read1.fastq'
    assert li.files[1].filename == 'mysample.read2.fastq'

def test_listinput_string():
    li = ruleinput.ListInput('reads', '-i')
    read1 = SimpleTemplateFile('read1')
    read2 = SimpleTemplateFile('read2')

    mock_templates = {
        'reads': [read1, read2],
    }
    li.set_template(mock_templates)

    mock_input = {
        'reads': ['mysample.read1.fastq', 'mysample.read2.fastq'],
    }
    li.assign(mock_input, sequential=True)

    assert str(li) == '-i mysample.read1.fastq mysample.read2.fastq'