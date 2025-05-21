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
)

from function_dump import (
    inspect_summary,
    inspect_verbose,
)

@click.group(context_settings={"help_option_names": ["--help"]})
@click.version_option("1.0")
def cli():
    """Database Pipeline CLI."""
    pass

RECIPES = {
    "test": [
        Walker,
        Concatenation,
        Write,
    ],
    "complete": [
        Walker,
        Concatenation,
        FlattenDuplicates,
        RmPurificationTags,
        CDRComputation,
        AntigenComputation,
        AssignIDs,
        Write,
    ],
}

def build_pipeline(recipe: str, path: str) -> Pipeline:
    steps = [step(path) if step is Walker else step() for step in RECIPES[recipe]]
    return Pipeline(steps)

@cli.command(name="run", help="Run a configured pipeline recipe on all CSVs under a directory.")
@click.option(
    "--test",
    "recipe",
    flag_value="test",
    default=True,
    help='Parametre for testing new functions added to the pipeline'
)
@click.option(
    "--complete",
    "recipe",
    flag_value="complete",
    help='Pefrom parsing, filtration, concatenation, chemical characteristic calculations, purification tag removal, novel ID assignment and print to new csv ready for SQL'
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
@click.option("-v", "--verbose", is_flag=True, help="Print every field for each row")
def inspect_cmd(csv_path: str, verbose: bool):
    """Quickly inspect a single CSV."""
    df = pd.read_csv(csv_path)
    if verbose:
        click.echo("🛈 Verbose inspection:\n")
        inspect_verbose(df)
    else:
        click.echo("🛈 Summary inspection:\n")
        inspect_summary(df)

if __name__ == "__main__":
    cli()