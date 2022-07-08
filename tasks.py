"""Definition of invoke tasks for task and tool running automation"""

from pathlib import Path

from invoke import task


THIS_DIR = Path(__file__).parent


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
