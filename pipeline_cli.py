import sys
import click
import pandas as pd
import time

from database_pipeline import (
    Pipeline,
    Walker,
    Concatenation,
    CDRComputation,
    AntigenComputation,
    FlattenDuplicates,
    Write,
    RmPurificationTags,
    AssignIDs,
    ComputeRelationships,
    WorkWithDatabase,
)

from function_dump import (
    inspect_summary,
    inspect_verbose,
)

@click.group(context_settings={"help_option_names": ["--help"]})
@click.version_option("1.5")
def cli():
    """Database Pipeline CLI."""
    pass

RECIPES = {
    "update": [
        Walker,
        WorkWithDatabase
    ],
    "complete": [
        Walker,
        Concatenation,
        FlattenDuplicates,
        RmPurificationTags,
        CDRComputation,
        AntigenComputation,
        AssignIDs,
        ComputeRelationships,
        Write,
    ],
        "half": [
        Walker,
        Concatenation,
        FlattenDuplicates,
        RmPurificationTags,
        CDRComputation,
        AntigenComputation,
        Write,
    ],
}

def build_pipeline(recipe: str, path: str) -> Pipeline:
    steps = [step(path) if step is Walker else step() for step in RECIPES[recipe]]
    return Pipeline(steps)

@cli.command(name="run", help="Run a configured pipeline recipe on all CSVs under a directory.")
@click.option(
    "--complete",
    "recipe",
    flag_value="complete",
    help='Pefrom parsing, filtration, concatenation, chemical characteristic calculations, purification tag removal, novel ID assignment and print to new antibody/antigen/relationships CSV files that ready for SQL'
)
@click.option(
    "--half",
    "recipe",
    flag_value="half",
    help='Pefrom parsing, filtration, concatenation, chemical characteristic calculations and print to new antigen/antibody CSV files'
)
@click.argument("path", type=click.Path(exists=True, file_okay=False))
def run(recipe: str, path: str):
    click.echo(f"▶ Using recipe '{recipe}' on '{path}'")
    start = time.time()
    pipeline = build_pipeline(recipe, path)
    pipeline.run()
    elapsed = time.time() - start
    click.echo(f"✔ Done in {elapsed:.2f} seconds.")

@cli.command(name="inspect", help="Inspect a single CSV file without running the pipeline.")
@click.argument("csv_path", type=click.Path(exists=True, dir_okay=False))
@click.option("-v", "--verbose", is_flag=True, help="Print every field for each row of the CSV file")
def inspect_cmd(csv_path: str, verbose: bool):
    df = pd.read_csv(csv_path)
    if verbose:
        click.echo("🛈 Verbose inspection:\n")
        inspect_verbose(df)
    else:
        click.echo("🛈 Summary inspection:\n")
        inspect_summary(df)

@cli.command(name="database", help="Update your SQL database automatically given DDL instructions within the sql_files/ directory, with data computed from this pipeline.")
@click.option("--update", "recipe", flag_value="update", help='Perform full database update with computed data'
)
@click.argument("path", type=click.Path(exists=True, file_okay=False))
def database_work(recipe: str, path: str):
    click.echo(f"▶ Using recipe '{recipe}' on '{path}'")
    start = time.time()
    pipeline = build_pipeline(recipe, path)
    pipeline.run()
    elapsed = time.time() - start
    click.echo(f"✔ Done in {elapsed:.2f} seconds.")


# Add a 'help' subcommand to show help for all commands or a specific command
@cli.command(name="help", help="Show help for all commands or a specific command.")
@click.argument("command", required=False)
def help_cmd(command):
    """
    Show the top-level help or help for a specific subcommand.
    """
    ctx = click.get_current_context()
    # If no command specified, show top-level help
    if not command:
        click.echo(ctx.parent.get_help())
    else:
        cmd = cli.get_command(ctx, command)
        if cmd is None:
            click.echo(f"No such command: {command}", err=True)
        else:
            # Create a new context for the subcommand to retrieve its help
            cmd_ctx = click.Context(cmd, info_name=command)
            click.echo(cmd.get_help(cmd_ctx))

if __name__ == "__main__":
    cli()