from snakewrap import exception

class RuleInput():
    def __init__(self, rule_key, *command_keys, name=None, desc=None):
        self.rule_key = rule_key
        self.command_keys = list(command_keys)
        self.files, self.names = None, None
        self.name = name
        self.desc = desc

    def describe(self):
        if self.name:
            return 'Input %s (for %s)' % (self.rule_key, self.name)
        else:
            return 'Input %s' % self.rule_key

    def __str__(self):
        raise NotImplementedError

    def set_template(self, snakemake_input):
        raise NotImplementedError

    def assign(self, snakemake_input, sequential=True):
        if self.files is None:
            raise exception.TemplateNotSetException('Set template files before assigning.')

        if sequential: 
            for f, input_f in zip(self.files, snakemake_input[self.rule_key]):
                f.assign(input_f)

class SimpleRuleInput(RuleInput):
    def __init__(self, rule_key, *command_keys, name=None, desc=None):
        super(SimpleRuleInput, self).__init__(rule_key, *command_keys, name=name, desc=desc)

        # Sanity check for the numbers of command keys.
        if len(command_keys) != 1:
            raise exception.RuleInputException('%s requires only one command key.' % self.describe())

    def set_template(self, snakemake_input):
        try:
            self.files = snakemake_input[self.rule_key]
            self.names = [f.name for f in self.files]
        except KeyError:
            raise exception.RuleInputException('Missing required input: %s' % self.rule_key)

        # Sanity check for the numbers of input files.
        if len(self.names) != 1:
            raise exception.RuleInputException('%s requires only one input.' % self.describe())

    def __str__(self):
        return '-i %s' % (self.files[0].infer_raw_name())

class MultiRuleInput(RuleInput):
    def __init__(self, rule_key, *command_keys, name=None, desc=None):
        super(MultiRuleInput, self).__init__(rule_key, *command_keys, name=name, desc=desc)

    def set_template(self, snakemake_input):
        try:
            self.files = snakemake_input[self.rule_key]
            self.names = [f.name for f in self.files]
        except KeyError:
            raise exception.RuleInputException('Missing required input: %s' % self.rule_key)
        
        if len(self.names) < 2:
            raise exception.RuleInputException('%s requires more than one inputs.' % self.describe())

    def __str__(self):
        return ' '.join([' '.join([k, f.infer_raw_name()]) for k, f in zip(self.command_keys, self.files)])