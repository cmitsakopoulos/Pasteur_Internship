import click
import pandas as pd
import time
import os

from app_components.database_pipeline import (
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
    CleanUp,
    PreWalker,
    LevenshteinDistance,
)

from app_components.function_dump import (
    inspect_summary,
    inspect_verbose,
)

@click.group(context_settings={"help_option_names": ["--help"]})
@click.version_option("2.0")
def cli():
    """Database Pipeline CLI."""
    pass

RECIPES = {
    "update": [
        Walker,
        WorkWithDatabase
    ],
    "normal": [
        CleanUp,
        PreWalker,
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
    "rerun": [
        PreWalker,
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
    "distance": [PreWalker, Walker, Concatenation, LevenshteinDistance],
}
def build_pipeline(recipe: str, path: str) -> Pipeline:
    steps = [step(path) if step is PreWalker else step() for step in RECIPES[recipe]]
    return Pipeline(steps)

@cli.command(name="run", help="Run a configured pipeline recipe on all CSVs under a directory.")
@click.option(
    "--normal",
    "recipe",
    flag_value="normal",
    default="normal",
    help='Perform all advertised functions apart from SQL injection, account for prior application runs to prevent recomputations.'
)
@click.option(
    "--rerun",
    "recipe",
    flag_value="rerun",
    help='Perform all advertised functions apart from SQL injection, WIPE all application component files AND results if you have made radical chages to code.'
)
@click.argument("path", required=False)
def run(recipe: str, path: str):
    default_internal = os.path.join(os.path.dirname(__file__), "Internal_Files")
    if path is None:
        if recipe == "rerun":
            path = default_internal
        else:
            raise click.BadParameter("PATH argument is required when --normal is selected")
    else:
        if not os.path.isdir(path):
            raise click.BadParameter(f"Path '{path}' does not exist or is not a directory")
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

@cli.command(name="distance", help="Calculate Levenshtein distance. Uses internal files directory if no path is provided.")
@click.argument("path", required=False, type=click.Path(exists=True, dir_okay=False))
def distance(path: str):
    recipe = "distance"
    if path is None:
        path = os.path.join(os.path.dirname(__file__), "Internal_Files")
        if not os.path.isdir(path):
            raise click.BadParameter(f"The default directory '{path}' does not exist. Please create it or provide a path.")
    click.echo(f"▶ Using recipe '{recipe}' on '{path}'")
    start = time.time()
    pipeline = build_pipeline(recipe, path)
    pipeline.run()
    elapsed = time.time() - start
    click.echo(f"✔ Done in {elapsed:.2f} seconds.")

@cli.command(name="help", help="Show help for all commands or a specific command.")
@click.argument("command", required=False)
def help_cmd(command):
    """
    Show the top-level help or help for a specific subcommand.
    """
    ctx = click.get_current_context()
    if not command:
        click.echo(ctx.parent.get_help())
    else:
        cmd = cli.get_command(ctx, command)
        if cmd is None:
            click.echo(f"No such command: {command}", err=True)
        else:
            cmd_ctx = click.Context(cmd, info_name=command)
            click.echo(cmd.get_help(cmd_ctx))

if __name__ == "__main__":
    cli()