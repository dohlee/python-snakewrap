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

