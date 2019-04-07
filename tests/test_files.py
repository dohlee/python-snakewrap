import os
import pytest
from snakewrap import files, exception

def test_file_init():
    f = files.File(
        name='bam_output',
        regex=r'(?P<basename>.+?)_output.bam',
    )

    assert f.name == 'bam_output'
    assert f.regex == r'(?P<basename>.+?)_output.bam'

def test_file_match():
    f = files.File(
        name='bam_output',
        regex=r'(?P<basename>.+?)_output.bam',
    )

    names = ['test_output_report.txt', 'test_output.bam']
    f.match(names)

    assert f.filename == 'test_output.bam'
    assert f.get_wildcard('basename') == 'test'

def test_file_match_failed_because_of_multiple_match():
    f = files.File(
        name='bam_output',
        regex=r'(?P<basename>.+?)_output.bam'
    )

    names = ['test1_output.bam', 'test2_output.bam']

    with pytest.raises(exception.AmbiguousFileNameException):
        f.match(names)

def test_file_match_failed_because_of_no_match():
    f = files.File(
        name='bam_output',
        regex=r'(?P<basename>.+?)_output.bam'
    )

    names = ['test1_output.bam', 'test2_output.bam']

    with pytest.raises(exception.AmbiguousFileNameException):
        f.match(names)

def test_file_init_with_raw_name():
    f = files.File(
        name='bam_output',
        regex=r'(?P<basename>.+?)_output.bam',
        raw_name='{basename}.bam'
    )

    assert f.raw_name == '{basename}.bam'

def test_file_infer_raw_name():
    f = files.File(
        name='bam_output',
        regex=r'(?P<basename>.+?)_output.bam',
        raw_name='{basename}.bam'
    )

    names = ['test1_output.bam', 'test1_report.bam']
    f.match(names)

    assert f.infer_raw_name() == 'test1.bam'

def test_file_infer_raw_name_failed_because_of_unmatched_file():
    f = files.File(
        name='bam_output',
        regex=r'(?P<basename>.+?)_output.bam',
        raw_name='{basename}.bam'
    )

    with pytest.raises(exception.FileUnmatchedException):
        assert f.infer_raw_name() == 'test1.bam'

def test_file_rename_command():
    f = files.File(
        name='bam_output',
        regex=r'(?P<basename>.+?)_output.bam',
        raw_name='{basename}.bam'
    )

    names = ['test1_output.bam', 'test1_report.bam']
    f.match(names)

    assert f.rename_command() == 'mv test1.bam test1_output.bam'

def test_file_with_multiple_group_regex():
    f = files.File(
        name='chip_output',
        regex=r'(?P<cell>.+?)_(?P<target>.+?).bam',
        raw_name='{cell}_{target}_output.bam',
    )

    names = ['HEK293T_H3K27ac.bam', 'HEK293T_H3K27ac_report.txt']
    f.match(names)

    assert f.filename == 'HEK293T_H3K27ac.bam'
    assert f.get_wildcard('cell') == 'HEK293T'
    assert f.get_wildcard('target') == 'H3K27ac'
    assert f.infer_raw_name() == 'HEK293T_H3K27ac_output.bam'
    assert f.rename_command() == 'mv HEK293T_H3K27ac_output.bam HEK293T_H3K27ac.bam'

def test_file_with_directory():
    f = files.File(
        name='chip_output',
        regex=r'/project/test/mapping_result/(?P<cell>.+?)_(?P<target>.+?).bam',
        raw_name='/project/test/mapping_result/{cell}_{target}_output.bam',
    )

    bam_name = '/project/test/mapping_result/HEK293T_H3K27ac.bam'
    raw_bam_name = '/project/test/mapping_result/HEK293T_H3K27ac_output.bam'
    report_name = '/project/test/mapping_result/HEK293T_H3K27ac_report.txt'

    names = [bam_name, report_name]
    f.match(names)

    assert f.filename == bam_name
    assert f.get_wildcard('cell') == 'HEK293T'
    assert f.get_wildcard('target') == 'H3K27ac'
    assert f.infer_raw_name() == raw_bam_name
    assert f.rename_command() == 'mv %s %s' % (raw_bam_name, bam_name)

def test_directory_init():
    dir = files.Directory(
        name='output_dir',
        regex=r'(?P<basename>.+?)_output$',
    )

    names = ['test_output.bam', 'test_output_report.txt', 'test_output']
    dir.match(names)

    assert dir.filename == 'test_output'
    assert dir.get_wildcard('basename') == 'test'

def test_directory_check_subdirectory():
    dir = files.Directory(
        name='output_dir',
        regex=r'(?P<basename>.+?)_output$',
    )

    if not os.path.exists('test_output'):
        os.makedirs('test_output/test1')
    names = ['test_output.bam', 'test_output_report.txt', 'test_output']
    dir.match(names)
    
    assert dir.check_subdirectory('test1') == True
    assert dir.check_subdirectory('test2') == False

    os.removedirs('test_output/test1')


