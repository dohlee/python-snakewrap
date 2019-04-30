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

    def get_options(self):
        raise NotImplementedError

    def _input_reader(self, snakemake_input):
        for key, filenames in snakemake_input.items():
            yield key, util.alwayslist(filenames)
    
    def assign(self, snakemake_input):
        for key, filenames in self._input_reader(snakemake_input):
            if key not in self.rule_keys:
                raise exception.RuleInputException('Unexpected input rule key: %s' % key)

            this_template = util.alwayslist(self.template[key])
            for filename, parameter in zip(filenames, this_template):
                parameter.f.assign(filename)

class SimpleRuleInput(RuleInput):
    def __init__(self, template, name=None, desc=None):
        super(SimpleRuleInput, self).__init__(template, name=name, desc=desc)

    def get_options(self):
        options = []
        for _, this_template in self.template.items():
            for parameter in util.alwayslist(this_template):
                if not parameter.used:
                    continue

                fname, priority = parameter.f.infer_raw_name(), parameter.priority
                if parameter.option is not None:
                    option = ' '.join([parameter.option, fname])
                    options.append((option, priority))
                else:
                    options.append((fname, priority))

        return options

class MultiRuleInput(RuleInput):
    def __init__(self, rule_key, *command_keys, name=None, desc=None):
        super(MultiRuleInput, self).__init__(rule_key, *command_keys, name=name, desc=desc)

    def set_template(self, snakemake_input):
        try:
            self.files = snakemake_input[self.rule_keys]
            self.names = [f.name for f in self.files]
        except KeyError:
            raise exception.RuleInputException('Missing required input: %s' % self.rule_keys)
        
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