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
        
        for k, v in self.params.items():
            value, opt = v
            if callable(value):
                self.params[k] = (self._infer_parameter_with_function(value, snakemake.input, snakemake.output), opt)

    def _infer_parameter_with_function(self, func, snakemake_input, snakemake_output):
        return func(snakemake_input, snakemake_output)

    def __str__(self):
        tmp = [] if not self.extra else [self.extra]
        for _, v in self.params.items():
            value, opt = v
            if opt is None:
                tmp.append(value)
            else:
                tmp.append('%s %s' % (opt, value))

        return ' '.join(tmp)

class SimpleRuleParams(RuleParams):
    def __init__(self, extra=True, **kwargs):
        super(SimpleRuleParams, self).__init__(extra, **kwargs)

