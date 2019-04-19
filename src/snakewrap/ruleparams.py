class RuleParams():
    def __init__(self, extra=True):
        self.include_extra = extra
        self.extra = None
        pass
    
    def assign(self, snakemake_params):
        if self.include_extra:
            self.extra = snakemake_params['extra']


class SimpleRuleParams(RuleParams):
    def __init__(self, extra=True):
        super(SimpleRuleParams, self).__init__(extra)