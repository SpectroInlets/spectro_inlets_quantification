"""Definition of invoke tasks for task and tool running automation"""

import platform
from pathlib import Path
from shutil import rmtree

from invoke import task


THIS_DIR = Path(__file__).parent



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


@task(aliases=("bd",))
def build_docs(context):
    with context.cd(THIS_DIR / "docs"):
        context.run("sphinx-build -M html source build")


@task(aliases=("pip", "deps", "requirements"))
def dependencies(context):
    """Install all requirements and development requirements"""
    context.run("python -m pip install --upgrade pip")
    context.run("python -m pip install --upgrade -r requirements.txt")
    context.run("python -m pip install --upgrade -r requirements-dev.txt")



CLEAN_PATTERNS = ("__pycache__", "*.pyc", "*.pyo", ".mypy_cache", "build")


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
