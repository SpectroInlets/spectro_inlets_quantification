"""This module contains shared fixtures"""

import pytest

from spectro_inlets_quantification.config import Config
from spectro_inlets_quantification.tools import Singleton

SINGLETONS = (Config,)


@pytest.fixture(scope="function")
def reset_singletons():
    for singleton in SINGLETONS:
        if singleton in Singleton._instances:
            del Singleton._instances[singleton]
    yield
    for singleton in SINGLETONS:
        if singleton in Singleton._instances:
            del Singleton._instances[singleton]
