import numpy as np
import pytest

from prigen.generators import PrimersGenerator
from prigen.utils import gc_content


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
    for n_primers in range(1, 20):
        for length in range(4, 20):
            for gc_percentage in np.linspace(0, 1, 20):
                generator = PrimersGenerator(
                    length=length,
                    number_of_primers=n_primers,
                    gc_percentage=gc_percentage
                )
                primers = generator.generate_primers()

                assert len(primers) <= n_primers
                assert all(len(primer) == length for primer in primers)
                assert all(
                    (gc_content(primer) - gc_percentage < 1/length) for primer in primers
                )


def test_primers_generator_with_temperature_bounds():
    for gc_percentage in np.linspace(0, 1, 20):
        for min_temperature in np.linspace(5, 60):
            max_temperature = min_temperature + 5

            generator = PrimersGenerator(
                length=20,
                number_of_primers=20,
                gc_percentage=gc_percentage,
                min_temperature=min_temperature,
                max_temperature=max_temperature
            )
            primers = generator.generate_primers()

            for temperature in primers.values():
                assert min_temperature <= temperature <= max_temperature
