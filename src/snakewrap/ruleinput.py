from snakewrap import exception

class RuleInput():
    def __init__(self, rule_key, *command_keys, name=None, desc=None):
        self.rule_key = rule_key
        self.command_keys = list(command_keys)
        self.name = name
        self.desc = desc

    def describe(self):
        if self.name:
            return 'Input %s (for %s)' % (self.rule_key, self.name)
        else:
            return 'Input %s' % self.rule_key

    def __str__(self):
        raise NotImplementedError

    def attach(self, snakemake_input):
        raise NotImplementedError

class SimpleRuleInput(RuleInput):
    def __init__(self, rule_key, *command_keys, name=None, desc=None):
        super(SimpleRuleInput, self).__init__(rule_key, *command_keys, name=name, desc=desc)

        # Sanity check for the numbers of command keys.
        if len(command_keys) != 1:
            raise exception.RuleInputException('%s requires only one command key.' % self.describe())

    def attach(self, snakemake_input):
        try:
            self.names = snakemake_input[self.rule_key]
        except KeyError:
            raise exception.RuleInputException('Missing required input: %s' % self.rule_key)

        # Sanity check for the numbers of input files.
        if len(self.names) != 1:
            raise exception.RuleInputException('%s requires only one input.' % self.describe())

    def __str__(self):
        return '-i %s' % self.names[0]

class MultiRuleInput(RuleInput):
    def __init__(self, rule_key, *command_keys, name=None, desc=None):
        super(MultiRuleInput, self).__init__(rule_key, *command_keys, name=name, desc=desc)

    def attach(self, snakemake_input):
        try:
            self.names = snakemake_input[self.rule_key]
        except KeyError:
            raise exception.RuleInputException('Missing required input: %s' % self.rule_key)
        
        if len(self.names) < 2:
            raise exception.RuleInputException('%s requires more than one inputs.' % self.describe())
