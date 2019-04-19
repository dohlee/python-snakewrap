class Wrapper:
    def __init__(self, snakemake, command, input, output, params=None, threads=None):
        self.snakemake = snakemake
        self.input = input
        self.output = output
        self.params = params
        self.threads = threads
        self.command = command

        self.input.assign(snakemake['input'])
        self.output.assign(snakemake['output'])
        self.params.assign(snakemake['params'])
        self.threads.assign(snakemake['threads'])
    
    def run(self):
        pass
    
    def shell_command(self):
        command = [self.command]
        command.append(str(self.input))

        return ' '.join(command)