import sys
import click
import pandas as pd

from database_pipeline import (
    Pipeline,
    Walker,
    Parser,
    Concatenation,
    CDRComputation,
    AntigenComputation,
    FlattenDuplicates,
    InspectionStep,
)

from function_dump import (
    inspect_summary,
    inspect_verbose,
)

RECIPES = {
    "extract": [
        Walker,
        Parser,
    ],
    "chem": [
        Walker,
        Parser,
        Concatenation,
        CDRComputation,
        AntigenComputation,
    ],
    "full": [
        Walker,
        Parser,
        Concatenation,
        CDRComputation,
        AntigenComputation,
        FlattenDuplicates,
    ],
    # An inspect recipe—we’ll inject the InspectionStep in build_pipeline()
    "inspect": [
        Walker,
        Parser,
        Concatenation,
        InspectionStep,
    ],
}

def build_pipeline(recipe: str, path: str, inspect_mode: str = "summary") -> Pipeline:
    """
    Instantiate and return a Pipeline configured according to the chosen recipe.
    """
    steps = []
    for step in RECIPES[recipe]:
        if step is Walker:
            steps.append(step(path))
        elif step is InspectionStep: #INGNORE FOR NOW< STILL PLAECHOLDER
            steps.append(step(mode=inspect_mode))
        else:
            steps.append(step())
    return Pipeline(steps)

@click.group()
@click.version_option()
def cli():
    """Database Pipeline CLI."""
    pass

@cli.command()
@click.argument("recipe", type=click.Choice(list(RECIPES.keys())))
@click.argument("path",   type=click.Path(exists=True, file_okay=False))
@click.option(
    "--inspect", "-i",
    type=click.Choice(["summary", "verbose"]),
    default="summary",
    show_default=True,
    help="Read CSVs, no matter the format.",
)

def run(recipe: str, path: str, inspect_mode: str):
    """
    Run a configured pipeline RECIPE on all CSVs under PATH.
    """
    click.echo(f"▶ Running recipe '{recipe}' on '{path}'")
    pipeline = build_pipeline(recipe, path, inspect_mode)
    data = pipeline.run()
    product = getattr(pipeline, "product", None)
    if isinstance(product, pd.DataFrame):
        click.echo(f"✔ Done: final product has {len(product)} rows")
    else:
        click.echo("✔ Done.")

@cli.command()
@click.argument("csv_path", type=click.Path(exists=True, dir_okay=False))
@click.option("--verbose", "-v", is_flag=True, help="Print every field for each row")
def inspect(csv_path: str, verbose: bool):
    """
    Quickly inspect a single CSV file (bypassing the full pipeline).
    """
    df = pd.read_csv(csv_path)
    if verbose:
        click.echo("🛈 Verbose inspection:\n")
        inspect_verbose(df)
    else:
        click.echo("🛈 Summary inspection:\n")
        inspect_summary(df)

if __name__ == "__main__":
    cli()