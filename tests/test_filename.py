from snakewrap import filename

def test_get_common_prefixes_with_no_extension():
    filenames = [
        'prefix.txt',
        'prefix.bam',
        'prefix.vcf',
    ]

    prefix = filename.get_common_prefix(filenames, include_ext=False)
    assert prefix == 'prefix'

def test_get_common_prefixes_with_extension():
    filenames = [
        'prefix.txt',
        'prefix.tat',
        'prefix.tbt',
        'prefix.tct',
    ]

    prefix = filename.get_common_prefix(filenames, include_ext=True)
    assert prefix == 'prefix.t'

def test_get_common_prefixes_with_full_path():
    filenames = [
        'a/b/c/d/prefix.txt',
        'a/b/c/d/prefix.bam',
        'a/b/c/d/prefix.vcf',
    ]

    prefix = filename.get_common_prefix(filenames, full_path=True)
    assert prefix == 'a/b/c/d/prefix'

def test_get_core_name():
    name = 'test.txt'
    core_name = filename.get_core_name(name)
    assert core_name == 'test'

def test_get_core_name_with_full_path():
    name = '/a/b/c/d/test.txt'
    core_name = filename.get_core_name(name)
    assert core_name == 'test'

def test_get_core_name_with_multiple_dots():
    name = 't.e.s.t.txt'
    core_name = filename.get_core_name(name)
    assert core_name == 't.e.s.t'

def test_extract_with_regex():
    pattern = 'test_(.+?).txt'
    name = 'test_test1.txt'
    assert filename.extract_with_regex(name, pattern) == 'test1'

def test_extract_with_regex_sra_run():
    pattern = '(SRR[0-9]+)'
    name = 'SRR1020421_trimmed.fastq'
    assert filename.extract_with_regex(name, pattern) == 'SRR1020421'