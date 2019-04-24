import pytest

from snakewrap import rulethreads

class MockSnakemake:
    def __init__(self, input=None, output=None, params=None, threads=None):
        self.input = input
        self.output = output
        self.params = params
        self.threads = threads

def test_simplerulethreads_init():
    rt = rulethreads.SimpleRuleThreads(command_key='-t')

def test_simplerulethreads_assign():
    sn = MockSnakemake(threads=3)
    rt = rulethreads.SimpleRuleThreads('-t')
    rt.assign(sn)

def test_rulethreads_attach_failed_because_of_non_int_assignment():
    rt = rulethreads.SimpleRuleThreads('-t')

    sn = MockSnakemake(threads='a')
    with pytest.raises(TypeError):
        rt.assign(sn)

    sn = MockSnakemake(threads=[1, 2])
    with pytest.raises(TypeError):
        rt.assign(sn)
    
    sn = MockSnakemake(threads={'thread': 4})
    with pytest.raises(TypeError):
        rt.assign(sn)

def test_rulethreads_string():
    sn = MockSnakemake(threads=4)
    rt = rulethreads.SimpleRuleThreads('-t')
    rt.assign(sn)

    assert str(rt) == '-t 4'

def test_rulethreads_assign_success_given_integer_convertible_string():
    sn = MockSnakemake(threads='2')
    rt = rulethreads.SimpleRuleThreads('-t')
    rt.assign(sn)

    assert str(rt) == '-t 2'

