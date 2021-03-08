import numpy as np
import pytest

from prigen.generators import PrimersGenerator


def test_primers_generator_with_wrong_length():
    with pytest.raises(ValueError):
        PrimersGenerator(length=-1, gc_percentage=.5)

    with pytest.raises(ValueError):
        PrimersGenerator(length=0, gc_percentage=.5)


def test_primers_generator_with_wrong_number():
    with pytest.raises(ValueError):
        PrimersGenerator(length=1, gc_percentage=.5, number_of_primers=0)

    with pytest.raises(ValueError):
        PrimersGenerator(length=1, gc_percentage=.5, number_of_primers=-999)


def test_primers_generator_with_wrong_temperature_bounds():
    with pytest.raises(ValueError):
        PrimersGenerator(length=1, gc_percentage=.5, min_temperature=5, max_temperature=3)

    with pytest.raises(ValueError):
        PrimersGenerator(length=1, gc_percentage=.5, min_temperature=5, max_temperature=5)


def test_primers_generator_with_wrong_gc_percentage():
    for percent in np.linspace(-2, -0.001, 30):
        with pytest.raises(ValueError):
            PrimersGenerator(length=1, gc_percentage=percent)

    for percent in np.linspace(1.001, 10, 30):
        with pytest.raises(ValueError):
            PrimersGenerator(length=1, gc_percentage=percent)


def test_primers_generator_with_any_temperature():
    for n_primers in range(1, 30):
        for length in range(4, 30):
            generator = PrimersGenerator(
                length=length, number_of_primers=n_primers, gc_percentage=.5
            )
            primers = generator.generate_primers()

            assert len(primers) == n_primers
            assert all(len(primer) == length for primer in primers)
