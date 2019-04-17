from snakewrap import exception

class RuleOutput():
    def __init__(self, template, name=None, desc=None):
        self.template = template
        self.rule_keys = list(template.keys())
        self.n_anonymous_keys = 0
        self.name = name
        self.desc = desc

    def describe(self):
        if self.name:
            return 'Output %s (for %s)' % (self.rule_keys, self.name)
        else:
            return 'Output %s' % self.rule_keys

    def __str__(self):
        raise NotImplementedError

    def assign(self, snakemake_input, sequential=True):
        if self.files is None:
            raise exception.TemplateNotSetException('Set template files before assigning.')

        if sequential: 
            for f, input_f in zip(self.files, snakemake_input[self.rule_key]):
                f.assign(input_f)

class SimpleRuleOutput(RuleOutput):
    def __init__(self, template, name=None, desc=None):
        super(SimpleRuleOutput, self).__init__(template, name, desc)

        # Sanity check for the numbers of command keys.
        if len(self.rule_keys) != 1:
            raise exception.RuleInputException('%s requires only one command key.' % self.describe())
    
    def set_template(self, template):
        try:
            self.files = template[self.rule_key]
            self.names = [f.name for f in self.files]
        except KeyError:
            raise exception.RuleInputException('Missing required output: %s' % self.rule_key)

        # Sanity check for the numbers of input files.
        if len(self.names) != 1:
            raise exception.RuleInputException('%s requires only one output.' % self.describe())

    def __str__(self):
        return '%s %s' % (self.command_keys[0], self.files[0].infer_raw_name())