"""
Entrypoint for phage_genome_assembly

Check out the wiki for a detailed look at customising this file:
https://github.com/beardymcjohnface/Snaketool/wiki/Customising-your-Snaketool
"""

import os
import click

from .util import snake_base, print_version, copy_config, run_snakemake, OrderedCommands, print_citation


def common_options(func):
    """Common command line args
    Define common command line args here, and include them with the @common_options decorator below.
    """
    options = [
        click.option('--output', help='Output directory', type=click.Path(),
                     default='phage_genome_assembly.out', show_default=True),
        click.option('--configfile', default='config.yaml', help='Custom config file', show_default=True),
        click.option('--threads', help='Number of threads to use', default=1, show_default=True),
        click.option('--use-conda/--no-use-conda', default=True, help='Use conda for Snakemake rules',
                     show_default=True),
        click.option('--conda-prefix', default=snake_base(os.path.join('workflow', 'conda')),
                     help='Custom conda env directory', type=click.Path(), show_default=False),
        click.option('--snake-default', multiple=True,
                     default=['--rerun-incomplete', '--printshellcmds', '--nolock', '--show-failed-logs'],
                     help="Customise Snakemake runtime args", show_default=True),
        click.argument('snake_args', nargs=-1)
    ]
    for option in reversed(options):
        func = option(func)
    return func


@click.group(cls=OrderedCommands)
def cli():
    """For more options, run:
    phage_genome_assembly command --help"""
    pass


help_msg_extra = """
\b
INSTALLING DATABASES REQUIRED
This command downloads the databases to the directory 'database' 
\b
phage_genome_assembly install 
\b
\b
CLUSTER EXECUTION:
phage_genome_assembly run ... --profile [profile]
For information on Snakemake profiles see:
https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles
\b
RUN EXAMPLES:
Required:           phage_genome_assembly run --input [file]
Specify threads:    phage_genome_assembly run ... --threads [threads]
Disable conda:      phage_genome_assembly run ... --no-use-conda 
Change defaults:    phage_genome_assembly run ... --snake-default="-k --nolock"
Add Snakemake args: phage_genome_assembly run ... --dry-run --keep-going --touch
Specify targets:    phage_genome_assembly run ... all print_targets
Available targets:
    all             Run everything (default)
\b
\b
PHAGE CONTIG QUALITY CHECK
Step 2 of the workflow checking the quality of the assembled contigs
\b
phage_genome_assembly contig ... --profile [profile]
For information on Snakemake profiles see:
https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles
\b
RUN EXAMPLES:
Required:           phage_genome_assembly contig --input [file] --contigs [file]
Specify threads:    phage_genome_assembly contig ... --threads [threads]
Disable conda:      phage_genome_assembly contig ... --no-use-conda 
Change defaults:    phage_genome_assembly contig ... --snake-default="-k --nolock"
Add Snakemake args: phage_genome_assembly contig ... --dry-run --keep-going --touch
Specify targets:    phage_genome_assembly contig ... all genomes
Available targets:
    all             Run everything (default)
\b
\b
PHAGE TAXA ASSIGNMENT
Step 3 of the workflow checking the quality of the assembled contigs
\b
phage_genome_assembly taxa ... --profile [profile]
For information on Snakemake profiles see:
https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles
\b
RUN EXAMPLES:
Required:           phage_genome_assembly taxa --phage [file]
Specify threads:    phage_genome_assembly taxa ... --threads [threads]
Disable conda:      phage_genome_assembly taxa ... --no-use-conda 
Change defaults:    phage_genome_assembly taxa ... --snake-default="-k --nolock"
Add Snakemake args: phage_genome_assembly taxa ... --dry-run --keep-going --touch
Specify targets:    phage_genome_assembly taxa ... all 
Available targets:
    all             Run everything (default)
\b
\b
PHAGE ANNOTATION
Step 4 of the workflow checking the quality of the assembled contigs
\b
phage_genome_assembly annotate ... --profile [profile]
For information on Snakemake profiles see:
https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles
\b
RUN EXAMPLES:
Required:           phage_genome_assembly annotate --phage [file]
Specify threads:    phage_genome_assembly annotate ... --threads [threads]
Disable conda:      phage_genome_assembly annotate ... --no-use-conda 
Change defaults:    phage_genome_assembly annotate ... --snake-default="-k --nolock"
Add Snakemake args: phage_genome_assembly annotate ... --dry-run --keep-going --touch
Specify targets:    phage_genome_assembly annotate ... all 
Available targets:
    all             Run everything (default)
\b
\b
"""

@click.command(epilog=help_msg_extra, context_settings=dict(help_option_names=["-h", "--help"], ignore_unknown_options=True))

@common_options
def install(configfile, threads, use_conda, conda_prefix, snake_default, **kwargs ):
    """
    The install function for databases
    """
    print("Checking and installing databases for atavide to directory")

    # run!
    run_snakemake(
        snakefile_path=snake_base(os.path.join('workflow', 'install.smk')),   # Full path to Snakefile
        configfile=snake_base(os.path.join('config', 'databases.yaml')),
        snake_default_args=snake_default,
        use_conda=use_conda,
    )

@click.command(epilog=help_msg_extra, context_settings=dict(help_option_names=["-h", "--help"], ignore_unknown_options=True))
@click.option('--input', '_input', help='Input file/directory', type=str, required=True)
@click.option('--preprocess', help="sequencing method", default='paired', show_default=True,
                     type=click.Choice(['paired', 'longread']))

@common_options
def run(_input, preprocess, configfile, output, threads, use_conda, conda_prefix, snake_default,
        snake_args, **kwargs):
    """Run phage_genome_assembly"""

    # copy default config file if missing
    copy_config(configfile, system_config=snake_base(os.path.join('config', 'config.yaml')))

    # Config to add or update in configfile
    merge_config = {
        'input': _input,
        'output': output,
        'sequencing': preprocess,
        }

    # run!
    run_snakemake(
        snakefile_path=snake_base(os.path.join('workflow', 'run.smk')),   # Full path to Snakefile
        configfile=configfile,
        merge_config=merge_config,
        threads=threads,
        use_conda=use_conda,
        conda_prefix=conda_prefix,
        snake_default_args=snake_default,
        snake_extra=snake_args,
    )

@click.command(epilog=help_msg_extra, context_settings=dict(help_option_names=["-h", "--help"], ignore_unknown_options=True))
@click.option('--input', '_input', help='Input file/directory', type=str, required=True)
@click.option('--preprocess', help="sequencing method", default='paired', show_default=True,
                     type=click.Choice(['paired', 'longread']))
@click.option('--phage-contigs', '_contigs', help="phage contigs picked from assemblies", required=True)

@common_options
def contig(_input, _contigs, preprocess, configfile, output, threads, use_conda, conda_prefix, snake_default,
        snake_args, **kwargs):
    """Run phage_contig_quality_check"""

    # copy default config file if missing
    copy_config(configfile, system_config=snake_base(os.path.join('config', 'config.yaml')))

    # Config to add or update in configfile
    merge_config = {
        'input': _input,
        'output': output,
        'sequencing': preprocess,
        'contigs': _contigs,
        }

    # run!
    run_snakemake(
        snakefile_path=snake_base(os.path.join('workflow', 'contig.smk')),   # Full path to Snakefile
        configfile=configfile,
        merge_config=merge_config,
        threads=threads,
        use_conda=use_conda,
        conda_prefix=conda_prefix,
        snake_default_args=snake_default,
        snake_extra=snake_args,
    )

@click.command(epilog=help_msg_extra, context_settings=dict(help_option_names=["-h", "--help"], ignore_unknown_options=True))
@click.option('--phage', '_phage', help='Input file/directory', type=str, required=True)

@common_options
def annotate(_phage, configfile, output, threads, use_conda, conda_prefix, snake_default,
        snake_args, **kwargs):
    """Run phage_contig_quality_check"""

    # copy default config file if missing
    copy_config(configfile, system_config=snake_base(os.path.join('config', 'config.yaml')))

    # Config to add or update in configfile
    merge_config = {
        'phage': _phage,
        'output': output,
        }

    # run!
    run_snakemake(
        snakefile_path=snake_base(os.path.join('workflow', 'annotate.smk')),   # Full path to Snakefile
        configfile=configfile,
        merge_config=merge_config,
        threads=threads,
        use_conda=use_conda,
        conda_prefix=conda_prefix,
        snake_default_args=snake_default,
        snake_extra=snake_args,
    )

@click.command(epilog=help_msg_extra, context_settings=dict(help_option_names=["-h", "--help"], ignore_unknown_options=True))
@click.option('--phage', '_phage', help='Input file/directory', type=str, required=True)

@common_options
def taxa(_phage, configfile, output, threads, use_conda, conda_prefix, snake_default,
        snake_args, **kwargs):
    """Run phage_contig_quality_check"""

    # copy default config file if missing
    copy_config(configfile, system_config=snake_base(os.path.join('config', 'config.yaml')))

    # Config to add or update in configfile
    merge_config = {
        'phage': _phage,
        'output': output,
        }

    # run!
    run_snakemake(
        snakefile_path=snake_base(os.path.join('workflow', 'taxa.smk')),   # Full path to Snakefile
        configfile=configfile,
        merge_config=merge_config,
        threads=threads,
        use_conda=use_conda,
        conda_prefix=conda_prefix,
        snake_default_args=snake_default,
        snake_extra=snake_args,
    )


@click.command()
@click.option('--configfile', default='config.yaml', help='Copy template config to file', show_default=True)
def config(configfile, **kwargs):
    """Copy the system default config file"""
    copy_config(configfile, system_config=snake_base(os.path.join('config', 'config.yaml')))


@click.command()
def citation(**kwargs):
    """Print the citation(s) for this tool"""
    print_citation()


cli.add_command(run)
cli.add_command(install)
cli.add_command(contig)
cli.add_command(taxa)
cli.add_command(annotate)
cli.add_command(config)
cli.add_command(citation)


def main():
    print_version()
    cli()


if __name__ == '__main__':
    main()
