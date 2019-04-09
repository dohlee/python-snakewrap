import pytest

from snakewrap import rulethreads

def test_rulethreads_init():
    rt = rulethreads.RuleThreads('-t')

def test_rulethreads_attach():
    rt = rulethreads.RuleThreads('-t')
    rt.attach(4)

def test_rulethreads_attach_failed_because_of_non_int_attachment():
    rt = rulethreads.RuleThreads('-t')

    with pytest.raises(TypeError):
        rt.attach('a')

    with pytest.raises(TypeError):
        rt.attach([1, 2])
    
    with pytest.raises(TypeError):
        rt.attach({'thread': 4})

def test_rulethreads_string():
    rt = rulethreads.RuleThreads('-t')
    rt.attach(4)

    assert str(rt) == '-t 4'

def test_rulethreads_attach_success_given_integer_convertible_string():
    rt = rulethreads.RuleThreads('-t')
    rt.attach('2')

    assert str(rt) == '-t 2'

