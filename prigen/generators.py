import numpy as np

from prigen.validators import (
    check_length,
    check_gc_percentage,
    check_number_of_primers,
    check_temperature_bounds
)


class PrimersGenerator:
    def __init__(
            self,
            length: int,
            gc_percentage: float,
            number_of_primers: int = 20,
            temperature_from: float = -999,
            temperature_to: float = 999
    ):
        self.length = length
        self.gc_percentage = gc_percentage
        self.number_of_primers = number_of_primers
        self.temperature_from = temperature_from
        self.temperature_to = temperature_to

        self._check_params()

    def _check_params(self):
        check_length(self.length)
        check_gc_percentage(self.gc_percentage)
        check_number_of_primers(self.number_of_primers)
        check_temperature_bounds(self.temperature_from, self.temperature_to)

    @staticmethod
    def generate_primer(length: int, gc_percentage: float) -> str:
        """Generate a nucleotide sequence with given GC-contend and length

        Parameters
        ----------
        length: int
            Length of the nucleotide sequence
        gc_percentage: float
            GC-content of the nucleotide sequence. Must be a number between 0 and 1
        """
        indexes = np.arange(length)
        np.random.shuffle(indexes)

        desired_gc = int(length * gc_percentage)
        desired_at = length - desired_gc

        if desired_at:
            desired_a = np.random.randint(desired_at + 1)
        else:
            desired_a = 0

        if desired_gc:
            desired_g = np.random.randint(desired_gc + 1)
        else:
            desired_g = 0

        g_indexes = indexes[0: desired_g]
        c_indexes = indexes[desired_g: desired_gc]
        a_indexes = indexes[desired_gc: desired_a + desired_gc]

        sequence = []
        for i in range(length):
            if i in g_indexes:
                sequence.append("G")
            elif i in c_indexes:
                sequence.append("C")
            elif i in a_indexes:
                sequence.append("A")
            else:
                sequence.append("T")

        return "".join(sequence)

    def generate(self):
        pass
