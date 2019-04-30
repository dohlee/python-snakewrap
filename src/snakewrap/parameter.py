class Parameter:
    def __init__(self, f, option, used=True, redirected=False, priority=None):
        self.f = f
        self.option = option
        self.used = used
        self.redirected = redirected
        self.priority = priority
    
    