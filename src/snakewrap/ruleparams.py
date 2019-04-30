from snakewrap.parameter import Parameter

class RuleParams():
    def __init__(self, extra=True, **kwargs):
        """Base class for rule parameters.
        Keyword arguments can get predefined values for parameters,
        or functions that accepts input and output rules(snakemake.input, snakemake.output) as parameter.
        All the keyword argument parameters should be given with corresponding
        parameter option.

        e.g. (lambda input: os.path.splitext(input)[0], '-T',)
        """
        self.include_extra = extra
        self.extra = None
        self.params = dict()

        for k, v in kwargs.items():
            self.params[k] = v
    
    def assign(self, snakemake):
        if self.include_extra:
            self.extra = snakemake.params.extra
        
        for key, parameter in self.params.items():
            if callable(parameter.f):
                self.params[key] = Parameter(
                    self._infer_parameter_with_function(parameter.f, snakemake), 
                    parameter.option,
                    parameter.used,
                    parameter.redirected,
                    parameter.priority
                )

    def _infer_parameter_with_function(self, func, snakemake):
        return func(snakemake)

    def get_options(self):
        options = [] if not self.extra else [(self.extra, None)]
        for _, parameter in self.params.items():
            if not parameter.used:
                continue
            
            if parameter.option is None:
                options.append((parameter.f, parameter.priority))
            else:
                options.append((' '.join([parameter.option, parameter.f]), parameter.priority))
        
        return options

class SimpleRuleParams(RuleParams):
    def __init__(self, extra=True, **kwargs):
        super(SimpleRuleParams, self).__init__(extra, **kwargs)
