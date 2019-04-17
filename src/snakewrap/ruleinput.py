from snakewrap import exception, util

class RuleInput():
    def __init__(self, template, name=None, desc=None):
        self.template = template
        self.rule_keys = list(template.keys())
        self.n_anonymous_keys = 0
        self.name = name
        self.desc = desc

    def describe(self):
        if self.name:
            return 'Input %s (for %s)' % (self.rule_keys, self.name)
        else:
            return 'Input %s' % self.rule_keys

    def __str__(self):
        raise NotImplementedError

    def assign(self, snakemake_input, sequential=True):
        if self.files is None:
            raise exception.TemplateNotSetException('Set template files before assigning.')

        if sequential: 
            for f, input_f in zip(self.files, snakemake_input[self.rule_key]):
                f.assign(input_f)

class SimpleRuleInput(RuleInput):
    def __init__(self, template, name=None, desc=None):
        super(SimpleRuleInput, self).__init__(template, name=name, desc=desc)

        # Sanity check for the numbers of command keys.
        if len(self.rule_keys) > 1:
            raise exception.RuleInputException('%s does not accept more than one command key.' % self.describe())
    
    def _template_reader(self, template):
        for template_files, command_keys in template:
            template_files = util.alwayslist(template_files)

            if command_key is None:
                self.n_anonymous_keys += 1
                command_key = '_anonymous_%d' % self.n_anonymous_keys
            
            yield template_files, command_keys

    def __str__(self):
        if len(self.command_keys) == 0:
            return '%s' % (self.files[0].infer_raw_name())
        else:
            return '%s %s' % (self.command_keys[0], self.files[0].infer_raw_name())

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

class ListInput(RuleInput):
    def __init__(self, rule_key, *command_keys, name=None, desc=None):
        super(ListInput, self).__init__(rule_key, *command_keys, name=name, desc=desc)

        if len(command_keys) != 1:
            raise AttributeError('The number of commandline option for ListInput should be one.')
    
    def set_template(self, snakemake_input):
        try:
            self.files = snakemake_input[self.rule_key]
            self.names = [f.name for f in self.files]
        except KeyError:
            raise exception.RuleInputException('Missing required input: %s' % self.rule_key)
    
    def __str__(self):
        return '%s %s' % (self.command_keys[0], ' '.join([f.filename for f in self.files]))