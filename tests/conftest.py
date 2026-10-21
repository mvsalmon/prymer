import pytest

@pytest.fixture
def sample_template():
    """Fixed 60bp test sequence template."""
    return "ATGC" * 15