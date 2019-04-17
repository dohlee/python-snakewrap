
class RuleThreads():
    def __init__(self, command_key, name=None, desc=None):
        self.command_key = command_key
        self.num_thread = None
        self.name, self.desc = name, desc
    
    def attach(self, num_thread):
        try:
            self.num_thread = int(num_thread)
        except ValueError:
            raise TypeError('Please give integer to threads')
            

    def __str__(self):
        return '%s %d' % (self.command_key, self.num_thread)

class SimpleRuleThreads(RuleThreads):
    def __init__(self, command_key, name=None, desc=None):
        super(SimpleRuleThreads, self).__init__(command_key, name=None, desc=None)