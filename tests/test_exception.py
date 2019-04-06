import pytest
from snakewrap import exception

def test_rule_input_exception_can_be_raised():
    with pytest.raises(exception.RuleInputException):
        raise exception.RuleInputException


def test_rule_output_exception_can_be_raised():
    with pytest.raises(exception.RuleOutputException):
        raise exception.RuleOutputException

def test_rule_parameter_exception_can_be_raised():
    with pytest.raises(exception.RuleParameterException):
        raise exception.RuleParameterException