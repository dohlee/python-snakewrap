from snakewrap import exception

class RuleThreads():
    def __init__(self, command_key, name=None, desc=None):
        self.command_key = command_key
        self.num_threads = None
        self.name, self.desc = name, desc
            
    def assign(self, snakemake_threads):
        raise NotImplementedError

    def __str__(self):
        return '%s %d' % (self.command_key, self.num_threads)

class SimpleRuleThreads(RuleThreads):
    def __init__(self, command_key, name=None, desc=None):
        super(SimpleRuleThreads, self).__init__(command_key, name=None, desc=None)
    
    def assign(self, snakemake_threads):
        try:
            self.num_threads = int(snakemake_threads)
        except ValueError:
            raise TypeError('Please give integer for the number of threads.')