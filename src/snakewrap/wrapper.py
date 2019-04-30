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
        if self.params:
            self.params.assign(snakemake)
        if self.threads: 
            self.threads.assign(snakemake)
    
    def run(self):
        pass
    
    def shell_command(self):
        # -1 denotes absolute priority. Those commands SHOULD come at the first of the command.
        options = [('(', -1), (self.command, -1)] if self.snakemake.log else [(self.command, -1)]
        if self.input:
            options.extend(self.input.get_options())
        if self.output:
            options.extend(self.output.get_options())
        if self.params:
            options.extend(self.params.get_options())
        if self.threads:
            options.extend(self.threads.get_options())
        
        # Sort commands according to the priorities.
        options_with_high_priority = list(sorted([opt for opt, pr in options if pr is not None and pr < 100]))
        options_with_low_priority = list(sorted([opt for opt, pr in options if pr is not None and pr >= 100]))
        options_without_priority = [opt for opt, pr in options if pr is None]
        options = options_with_high_priority + options_without_priority + options_with_low_priority

        options.append(self.output.redirect_command())
        rename_command = self.output.rename_command()
        if rename_command:
            options.append('&& ' + rename_command)
        if self.snakemake.log:
            options.append(')' + self.snakemake.log_fmt_shell(stdout=False, stderr=True))
        return ' '.join(filter(lambda x: x, options))  # Discard empty strings, then join.