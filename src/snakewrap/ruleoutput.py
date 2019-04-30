from snakewrap import exception, util
from snakewrap.files import RenamedTemplateFile

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

class SimpleRuleOutput(RuleOutput):
    def __init__(self, template, name=None, desc=None):
        super(SimpleRuleOutput, self).__init__(template, name, desc)

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
    
    def rename_command(self):
        tmp = []
        for _, this_template in self.template.items():
            for parameter in util.alwayslist(this_template):
                if isinstance(parameter.f, RenamedTemplateFile):
                    tmp.append(parameter.f.rename_command())
        
        return ' && '.join(tmp)
