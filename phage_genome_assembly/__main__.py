"""
Entrypoint for spae

Check out the wiki for a detailed look at customising this file:
https://github.com/beardymcjohnface/Snaketool/wiki/Customising-your-Snaketool
"""

import os
import click

from snaketool_utils.cli_utils import OrderedCommands, run_snakemake, copy_config, echo_click


def snake_base(rel_path):
    """Get the filepath to a Snaketool system file (relative to __main__.py)"""
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), rel_path)


def get_version():
    with open(snake_base("phage_genome_assembly.VERSION"), "r") as f:
        version = f.readline()
    return version


def print_citation():
    with open(snake_base("phage_genome_assembly.CITATION"), "r") as f:
        for line in f:
            echo_click(line)


def default_to_output(ctx, param, value):
    """Callback for click options; places value in output directory unless specified"""
    if param.default == value:
        return os.path.join(ctx.params["output"], value)
    return value


def common_options(func):
    """Common command line args
    Define common command line args here, and include them with the @common_options decorator below.
    """
    options = [
        click.option('--output', help='Output directory', type=click.Path(),
                     default='spae.out', show_default=True),
        click.option("--configfile", default="spae.config.yaml", show_default=False, callback=default_to_output,
                     help="Custom config file [default: (outputDir)/spae.config.yaml]",),
        click.option('--threads', help='Number of threads to use', default=1, show_default=True),
        click.option('--use-conda/--no-use-conda', default=True, help='Use conda for Snakemake rules',
                     show_default=True),
        click.option('--conda-prefix', default=snake_base(os.path.join('workflow', 'conda')),
                     help='Custom conda env directory', type=click.Path(), show_default=False),
        click.option('--snake-default', multiple=True,
                     default=['--rerun-incomplete', '--printshellcmds', '--nolock', '--show-failed-logs'],
                     help="Customise Snakemake runtime args", show_default=True),
        click.option("--log", default="spae.log", callback=default_to_output, hidden=True,),
        click.argument('snake_args', nargs=-1)
    ]
    for option in reversed(options):
        func = option(func)
    return func


@click.group(cls=OrderedCommands, context_settings=dict(help_option_names=["-h", "--help"]))
@click.version_option(get_version(), "-v", "--version", is_flag=True)
def cli():
    """Assembling pure culture phages from both Illumina and Nanopore sequencing technology
    \b
    For more options, run:
    spae command --help"""
    pass


help_msg_extra = """
\b
INSTALLING DATABASES REQUIRED
This command downloads the databases to the directory 'database' 
\b
spae install 
\b
\b
CLUSTER EXECUTION:
spae run ... --profile [profile]
For information on Snakemake profiles see:
https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles
\b
RUN EXAMPLES:
Required:           spae run --input [file]
Specify threads:    spae run ... --threads [threads]
Disable conda:      spae run ... --no-use-conda 
Change defaults:    spae run ... --snake-default="-k --nolock"
Add Snakemake args: spae run ... --dry-run --keep-going --touch
Specify targets:    spae run ... all print_targets
Available targets:
    all             Run everything (default)
\b
\b
PHAGE CONTIG QUALITY CHECK
Step 2 of the workflow checking the quality of the assembled contigs
\b
spae contig ... --profile [profile]
For information on Snakemake profiles see:
https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles
\b
RUN EXAMPLES:
Required:           spae contig --input [file] --contigs [file]
Specify threads:    spae contig ... --threads [threads]
Disable conda:      spae contig ... --no-use-conda 
Change defaults:    spae contig ... --snake-default="-k --nolock"
Add Snakemake args: spae contig ... --dry-run --keep-going --touch
Specify targets:    spae contig ... all genomes
Available targets:
    all             Run everything (default)
\b
\b
PHAGE TAXA ASSIGNMENT
Step 3 of the workflow checking the quality of the assembled contigs
\b
spae taxa ... --profile [profile]
For information on Snakemake profiles see:
https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles
\b
RUN EXAMPLES:
Required:           spae taxa --phage [file]
Specify threads:    spae taxa ... --threads [threads]
Disable conda:      spae taxa ... --no-use-conda 
Change defaults:    spae taxa ... --snake-default="-k --nolock"
Add Snakemake args: spae taxa ... --dry-run --keep-going --touch
Specify targets:    spae taxa ... all 
Available targets:
    all             Run everything (default)
\b
\b
PHAGE ANNOTATION
Step 4 of the workflow checking the quality of the assembled contigs
\b
spae annotate ... --profile [profile]
For information on Snakemake profiles see:
https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles
\b
RUN EXAMPLES:
Required:           spae annotate --phage [file]
Specify threads:    spae annotate ... --threads [threads]
Disable conda:      spae annotate ... --no-use-conda 
Change defaults:    spae annotate ... --snake-default="-k --nolock"
Add Snakemake args: spae annotate ... --dry-run --keep-going --touch
Specify targets:    spae annotate ... all 
Available targets:
    all             Run everything (default)
\b
\b
"""


@click.command(epilog=help_msg_extra, context_settings=dict(help_option_names=["-h", "--help"], ignore_unknown_options=True))
@common_options
def install(**kwargs ):
    """The install function for databases"""

    # run!
    run_snakemake(
        snakefile_path=snake_base(os.path.join('workflow', 'install.smk')),   # Full path to Snakefile
        system_config=snake_base(os.path.join('config', 'databases.yaml')),
        **kwargs
    )


@click.command(epilog=help_msg_extra, context_settings=dict(help_option_names=["-h", "--help"], ignore_unknown_options=True))
@click.option('--input', '_input', help='Input file/directory', type=str, required=True)
@click.option('--preprocess', help="sequencing method", default='paired', show_default=True,
                     type=click.Choice(['paired', 'longread']))
@common_options
def run(**kwargs):
    """Run spae"""

    # Config to add or update in configfile
    merge_config = {
        'input': kwargs["_input"],
        'output': kwargs["output"],
        'sequencing': kwargs["preprocess"],
        }

    # run!
    run_snakemake(
        snakefile_path=snake_base(os.path.join('workflow', 'run.smk')),   # Full path to Snakefile
        system_config=snake_base(os.path.join('config', 'config.yaml')),
        merge_config=merge_config,
        **kwargs
    )


@click.command(epilog=help_msg_extra, context_settings=dict(help_option_names=["-h", "--help"], ignore_unknown_options=True))
@click.option('--input', '_input', help='Input file/directory', type=str, required=True)
@click.option('--preprocess', help="sequencing method", default='paired', show_default=True,
                     type=click.Choice(['paired', 'longread']))
@click.option('--phage-contigs', '_contigs', help="phage contigs picked from assemblies", required=True)
@common_options
def contig(**kwargs):
    """Run phage_contig_quality_check"""

    # Config to add or update in configfile
    merge_config = {
        'input': kwargs["_input"],
        'output': kwargs["output"],
        'sequencing': kwargs["preprocess"],
        'contigs': kwargs["_contigs"],
        }

    # run!
    run_snakemake(
        snakefile_path=snake_base(os.path.join('workflow', 'contig.smk')),   # Full path to Snakefile
        system_config=snake_base(os.path.join('config', 'config.yaml')),
        merge_config=merge_config,
        **kwargs
    )


@click.command(epilog=help_msg_extra, context_settings=dict(help_option_names=["-h", "--help"], ignore_unknown_options=True))
@click.option('--phage', '_phage', help='Input file/directory', type=str, required=True)
@common_options
def annotate(**kwargs):
    """Run phage_contig_quality_check"""

    # Config to add or update in configfile
    merge_config = {
        'phage': kwargs["_phage"],
        'output': kwarags["output"],
        }

    # run!
    run_snakemake(
        snakefile_path=snake_base(os.path.join('workflow', 'annotate.smk')),   # Full path to Snakefile
        system_config=snake_base(os.path.join('config', 'config.yaml')),
        merge_config=merge_config,
        **kwargs
    )


@click.command(epilog=help_msg_extra, context_settings=dict(help_option_names=["-h", "--help"], ignore_unknown_options=True))
@click.option('--phage', '_phage', help='Input file/directory', type=str, required=True)
@common_options
def taxa(**kwargs):
    """Run phage_contig_quality_check"""

    # Config to add or update in configfile
    merge_config = {
        'phage': kwargs["_phage"],
        'output': kwargs["output"],
        }

    # run!
    run_snakemake(
        snakefile_path=snake_base(os.path.join('workflow', 'taxa.smk')),   # Full path to Snakefile
        system_config=snake_base(os.path.join('config', 'config.yaml')),
        merge_config=merge_config,
        **kwargs
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
    cli()


if __name__ == '__main__':
    main()
