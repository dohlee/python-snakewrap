class Wrapper:
    def __init__(self, snakemake, command, input, output, params=None, threads=None):
        self.snakemake = snakemake
        self.input = input
        self.output = output
        self.params = params
        self.threads = threads
        self.command = command

        self.input.assign(snakemake.input)
        self.output.assign(snakemake.output)
        self.params.assign(snakemake)
        self.threads.assign(snakemake)
    
    def run(self):
        pass
    
    def shell_command(self):
        command = ['(', self.command] if self.snakemake.log else [self.command]
        command.append(str(self.input))
        command.append(str(self.output))
        command.append(str(self.params))
        command.append(str(self.threads))

        rename_command = self.output.rename_command()
        if rename_command:
            command.append('&& ' + rename_command)
        if self.snakemake.log:
            command.append(')' + self.snakemake.log_fmt_shell(stdout=False, stderr=True))
        return ' '.join(filter(lambda x: x, command))  # Discard empty strings, then join.