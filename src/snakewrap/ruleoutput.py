from snakewrap import exception, util

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

    def _input_reader(self, snakemake_input):
        for key, filenames in snakemake_input.items():
            yield key, util.alwayslist(filenames)

    def assign(self, snakemake_input):
        for key, filenames in self._input_reader(snakemake_input):
            if key not in self.rule_keys:
                raise exception.RuleInputException('Unexpected input rule key: %s' % key)

            template_files_and_command_keys = util.alwayslist(self.template[key])
            for filename, (template_file, _) in zip(filenames, template_files_and_command_keys):
                template_file.assign(filename)

class SimpleRuleOutput(RuleOutput):
    def __init__(self, template, name=None, desc=None):
        super(SimpleRuleOutput, self).__init__(template, name, desc)

        # Sanity check for the numbers of command keys.
        if len(self.rule_keys) != 1:
            raise exception.RuleInputException('%s requires only one command key.' % self.describe())

    def __str__(self):
        tmp = []
        for _, template_files_and_command_keys in self.template.items():
            for template_file, command_key in util.alwayslist(template_files_and_command_keys):
                if command_key is not None:
                    tmp.append('%s %s' % (command_key, template_file.infer_raw_name()))
                else:
                    tmp.append(template_file.infer_raw_name())
        
        return ' '.join(tmp)
