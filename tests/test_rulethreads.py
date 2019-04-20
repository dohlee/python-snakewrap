import pytest

from snakewrap import rulethreads

def test_simplerulethreads_init():
    rt = rulethreads.SimpleRuleThreads('-t')

def test_simplerulethreads_assign():
    rt = rulethreads.SimpleRuleThreads('-t')
    rt.assign(4)

def test_rulethreads_attach_failed_because_of_non_int_assignment():
    rt = rulethreads.SimpleRuleThreads('-t')

    with pytest.raises(TypeError):
        rt.assign('a')

    with pytest.raises(TypeError):
        rt.assign([1, 2])
    
    with pytest.raises(TypeError):
        rt.assign({'thread': 4})

def test_rulethreads_string():
    rt = rulethreads.SimpleRuleThreads('-t')
    rt.assign(4)

    assert str(rt) == '-t 4'

def test_rulethreads_assign_success_given_integer_convertible_string():
    rt = rulethreads.SimpleRuleThreads('-t')
    rt.assign('2')

    assert str(rt) == '-t 2'

