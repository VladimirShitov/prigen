import numpy as np
import pytest

from prigen.generators import PrimersGenerator
from prigen.utils import gc_content, parse_blast_result, remove_keys_from_dict


@pytest.fixture
def blast_result():
    return "gene1\tchr7\t100.000\t96\t0\t0\t1\t96\t5529033\t5528938\t4.22e-43\t178\n" + \
           "gene2\tchr17\t100.000\t111\t0\t0\t1\t111\t7686764\t7686654\t2.36e-51\t206\n"


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


def test_number_of_generated_primers():
    for n_primers in range(1, 20):
        generator = PrimersGenerator(
            length=20,
            number_of_primers=n_primers,
            gc_percentage=.5
        )
        primers = generator.generate_primers()

        assert len(primers) <= n_primers


def test_length_of_generated_primers():
    for length in range(4, 20):
        generator = PrimersGenerator(
            length=length,
            number_of_primers=20,
            gc_percentage=.5
        )
        primers = generator.generate_primers()

        assert all(len(primer) == length for primer in primers)


def test_primers_generator_with_any_temperature():
    for length in range(2, 50):
        for gc_percentage in np.linspace(0, 1, 30):
            generator = PrimersGenerator(
                length=length,
                number_of_primers=20,
                gc_percentage=gc_percentage
            )
            primers = generator.generate_primers()

            for primer in primers:
                assert abs(gc_content(primer) - gc_percentage) < 1/length


def test_primers_generator_with_temperature_bounds():
    min_temperature, max_temperature = 30, 80

    for gc_percentage in np.linspace(0, 1, 20):
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


def test_gc_content():
    assert gc_content("A") == 0

    length = 20

    for n_g in range(10):
        for n_c in range(10):
            primer = "G" * n_g + "C" * n_c + "A" * (length - n_c - n_g)

            assert gc_content(primer) == (n_c + n_g) / length

    primer = "GCGCGCATGCAGTACGA"

    assert gc_content(primer) == gc_content(primer.lower())


def test_parse_blast_result(blast_result):
    assert parse_blast_result(blast_result) == {"gene1", "gene2"}


def test_parse_empty_blast_result():
    assert parse_blast_result("") == set()


def test_remove_keys_from_dict():
    d = {"a": 1, "b": 2, "c": 3}

    assert remove_keys_from_dict(d, ["b"]) == {"a": 1, "c": 3}
    assert remove_keys_from_dict(d, ["not_existing"]) == d
    assert remove_keys_from_dict(d, ["a", "not_existing"]) == {"b": 2, "c": 3}
    assert remove_keys_from_dict(d, d.keys()) == {}
