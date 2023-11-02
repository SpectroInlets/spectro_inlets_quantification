"""Definition of invoke tasks for task and tool running automation"""

import platform
from pathlib import Path
from shutil import rmtree

try:
    from rich import print as rprint
except ImportError:
    rprint = print
from invoke import task


THIS_DIR = Path(__file__).parent
SOURCE_DIR = THIS_DIR / "src" / "spectro_inlets_quantification"


@task(aliases=("ruff",))
def lint(context, fix=False):
    """Lint the source code with ruff"""
    rprint("\n[bold]Linting...")
    with context.cd(THIS_DIR):
        arg_string = ""
        if fix:
            arg_string += " --fix"
        result = context.run(f"ruff {SOURCE_DIR} {arg_string}")
        if result.return_code == 0:
            rprint("[bold green]Files linted. No errors.")
    return result.return_code


@task(aliases=("mypy", "tc"))
def type_check(context):
    """Run the static type checker on the source code"""
    rprint("\n[bold]Checking types...")
    with context.cd(THIS_DIR):
        result = context.run(f"mypy src/")
        if result.return_code == 0:
            rprint("[bold green]Files type checked. No errors.")
    return result.return_code


@task(aliases=("fc", "black"))
def format_code(context):
    """Format all of zilien_qt with black"""
    context.run(f"black {SOURCE_DIR}")


@task(aliases=("check_black", "cb", "ccf"))
def check_code_format(context):
    """Check that the code is black formatted"""
    rprint("\n[bold]Checking code style...")
    with context.cd(THIS_DIR):
        result = context.run(f"black --check {SOURCE_DIR}")
        if result.return_code == 0:
            rprint("[bold green]Code format checked. No issues.")
    return result.return_code


@task(
    aliases=["tests"],
    help={
        "color": "Whether to display pytest output in color, 'yes' or 'no'",
        "verbose": "Makes the pytest output verbose",
        "s_no_capture": (
            "Prevents pytest from capturing output (making it possible to see prints etc.)"
        ),
        "k_only_run": (
            "Only run tests that matches the expression in STRING. See the help for pytest's "
            "`-k` option to read more about the options for expression"
        ),
        "x_exit_on_first_error": "Make pytest exit on first error",
        "also_slow": "Also run slow tests, disabled by default",
    },
)
def test(
    context,
    color="yes",
    verbose=False,
    s_no_capture=False,
    k_only_run=None,
    x_exit_on_first_error=False,
    also_slow=False,
):
    """Run tests with pytest"""
    if platform.system() == "Windows":
        color = "no"
    args = []
    if verbose:
        args.append("--verbose")
    if s_no_capture:
        args.append("-s")
    if k_only_run:
        args.append(f"-k {k_only_run}")
    if x_exit_on_first_error:
        args.append("-x")
    if not also_slow:
        args += ["-m", r'"not slow"']
    print("### Testing ...")
    result = context.run(f'pytest --color "{color}" {" ".join(args)} {THIS_DIR}/tests')
    return result.return_code


@task(aliases=["check", "c"])
def checks(context):
    """Check code with black, flake8, mypy and run tests"""
    combined_return_code = check_code_format(context)
    combined_return_code += lint(context)
    # combined_return_code += type_check(context)
    combined_return_code += test(context)
    if combined_return_code == 0:
        rprint()
        rprint(r"+----------+")
        rprint(r"| All good |")
        rprint(r"+----------+")
    else:
        rprint()
        rprint(r"+---------------------+")
        rprint(r"| [bold red]Some checks [blink]FAILED![/blink][/bold red] |")
        rprint(r"| [bold]Check output above[/bold]  |")
        rprint(r"+---------------------+")


@task(aliases=("bd",))
def build_docs(context):
    with context.cd(THIS_DIR / "docs"):
        context.run("sphinx-build -M html source build")


@task(aliases=("pip", "deps"))
def dependencies(context):
    """Install all requirements and development requirements"""
    with context.cd(THIS_DIR):
        context.run("python -m pip install --upgrade pip")
        context.run("python -m pip install .")  # Also installs core dependencies
        context.run("python -m pip install -r requirements-dev.txt")
        context.run("python -m pip install -r requirements-doc.txt")


CLEAN_PATTERNS = ("__pycache__", "*.pyc", "*.pyo", ".mypy_cache", "build", "dist")


@task
def clean(context, dryrun=False):
    """Clean the repository"""
    if dryrun:
        print("CLEANING DRYRUN")
    for clean_pattern in CLEAN_PATTERNS:
        for cleanpath in THIS_DIR.glob("**/" + clean_pattern):
            if cleanpath.is_dir():
                print("DELETE DIR :", cleanpath)
                if not dryrun:
                    rmtree(cleanpath)
            else:
                print("DELETE FILE:", cleanpath)
                if not dryrun:
                    cleanpath.unlink()
                    pass


@task
def build(context):
    """Build the package"""
    clean(context)
    with context.cd(THIS_DIR):
        context.run("python -m build")
    print("Upload with: python3 -m twine upload --repository testpypi dist/*")
