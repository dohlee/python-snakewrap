import os
import pytest
from snakewrap.files import SimpleTemplateFile, RenamedTemplateFile, SimpleTemplateDirectory
from snakewrap import files, exception

def test_simpletemplatefile_init():
    f = SimpleTemplateFile(
        name='bam_output',
    )

    assert f.name == 'bam_output'

def test_simpletemplatefile_assign():
    f = SimpleTemplateFile(
        name='bam_output',
    )

    f.assign('mysample.bam')
    assert f.filename == 'mysample.bam'

def test_renamedtemplatefile_init():
    f = RenamedTemplateFile(
        name='bam_output',
        regex=r'(?P<basename>.+?)_output.bam',
        raw_name='{basename}.bam',
    )

def test_renamedtemplatefile_assign():
    f = RenamedTemplateFile(
        name='bam_output',
        regex=r'(?P<basename>.+?)_output.bam',
        raw_name='{basename}.bam',
    )

    name = 'test_output.bam'
    f.assign(name)

    assert f.filename == 'test_output.bam'
    assert f.get_wildcard('basename') == 'test'

def test_renamedtemplatefile_assign_failed_because_of_no_match():
    f = RenamedTemplateFile(
        name='bam_output',
        regex=r'(?P<basename>.+?)_output.bam',
        raw_name='{basename}.bam',
    )

    name = 'test1_report.bam'
    with pytest.raises(exception.FileNameMismatchException):
        f.assign(name)

def test_renamedtemplatefile_infer_raw_name():
    f = RenamedTemplateFile(
        name='bam_output',
        regex=r'(?P<basename>.+?)_output.bam',
        raw_name='{basename}.bam'
    )

    name = 'test1_output.bam'
    f.assign(name)

    assert f.infer_raw_name() == 'test1.bam'

def test_renamedtemplatefile_infer_raw_name_failed_because_of_unassigned_file():
    f = RenamedTemplateFile(
        name='bam_output',
        regex=r'(?P<basename>.+?)_output.bam',
        raw_name='{basename}.bam'
    )

    with pytest.raises(exception.FileUnassignedException):
        assert f.infer_raw_name() == 'test1.bam'

def test_renamedtemplatefile_rename_command():
    f = RenamedTemplateFile(
        name='bam_output',
        regex=r'(?P<basename>.+?)_output.bam',
        raw_name='{basename}.bam'
    )

    name = 'test1_output.bam'
    f.assign(name)

    assert f.rename_command() == 'mv test1.bam test1_output.bam'

def test_renamedtemplatefile_with_multiple_group_regex():
    f = RenamedTemplateFile(
        name='chip_output',
        regex=r'(?P<cell>.+?)_(?P<target>.+?).bam',
        raw_name='{cell}_{target}_output.bam',
    )

    names = 'HEK293T_H3K27ac.bam'
    f.assign(names)

    assert f.filename == 'HEK293T_H3K27ac.bam'
    assert f.get_wildcard('cell') == 'HEK293T'
    assert f.get_wildcard('target') == 'H3K27ac'
    assert f.infer_raw_name() == 'HEK293T_H3K27ac_output.bam'
    assert f.rename_command() == 'mv HEK293T_H3K27ac_output.bam HEK293T_H3K27ac.bam'

def test_file_with_directory():
    f = RenamedTemplateFile(
        name='chip_output',
        regex=r'/project/test/mapping_result/(?P<cell>.+?)_(?P<target>.+?).bam',
        raw_name='/project/test/mapping_result/{cell}_{target}_output.bam',
    )

    bam_name = '/project/test/mapping_result/HEK293T_H3K27ac.bam'
    raw_bam_name = '/project/test/mapping_result/HEK293T_H3K27ac_output.bam'

    name = bam_name
    f.assign(name)

    assert f.filename == bam_name
    assert f.get_wildcard('cell') == 'HEK293T'
    assert f.get_wildcard('target') == 'H3K27ac'
    assert f.infer_raw_name() == raw_bam_name
    assert f.rename_command() == 'mv %s %s' % (raw_bam_name, bam_name)

def test_simpletemplatedirectory_init():
    dir = SimpleTemplateDirectory(
        name='output_dir',
    )

def test_simpletemplatedirectory_assign():
    dir = SimpleTemplateDirectory(
        name='output_dir',
    )

    name = 'test_output'
    dir.assign(name)

    assert dir.filename == 'test_output'

def test_simpletemplatedirectory_check_subdirectory():
    dir = SimpleTemplateDirectory(
        name='output_dir',
    )

    if not os.path.exists('test_output'):
        os.makedirs('test_output/test1')

    name = 'test_output'
    dir.assign(name)
    
    assert dir.check_subdirectory('test1') == True
    assert dir.check_subdirectory('test2') == False

    os.removedirs('test_output/test1')


