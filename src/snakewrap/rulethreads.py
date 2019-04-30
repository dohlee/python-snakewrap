from snakewrap import exception

class RuleThreads():
    def __init__(self, command_key, name=None, desc=None):
        self.command_key = command_key
        self.num_threads = None
        self.name, self.desc = name, desc
            
    def assign(self, snakemake_threads):
        raise NotImplementedError

    def get_options(self):
        return [('%s %d' % (self.command_key, self.num_threads), None)]

    def __str__(self):
        return '%s %d' % (self.command_key, self.num_threads)

class SimpleRuleThreads(RuleThreads):
    def __init__(self, command_key, name=None, desc=None):
        super(SimpleRuleThreads, self).__init__(command_key, name=name, desc=desc)
    
    def assign(self, snakemake):
        try:
            self.num_threads = int(snakemake.threads)
        except ValueError:
            raise TypeError('Please give integer for the number of threads.')

class ScaledRuleThreads(RuleThreads):
    def __init__(self, command_key, scale, name=None, desc=None):
        super(ScaledRuleThreads, self).__init__(command_key, name=name, desc=desc)
        self.scale = scale
    
    def assign(self, snakemake):
        # Compute scaler value.
        scaler = self.scale(snakemake) if callable(self.scale) else self.scale

        try:
            self.num_threads = int(snakemake.threads)
        except ValueError:
            raise TypeError('Please give integer for the number of threads.')

        self.num_threads = max(1, int(self.num_threads * scaler))


        
