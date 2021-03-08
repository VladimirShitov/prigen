import pytest

from prigen.generators import PrimersGenerator


def test_primers_generator_with_wrong_length():
    with pytest.raises(ValueError):
        PrimersGenerator(length=-1, gc_percentage=.5)
