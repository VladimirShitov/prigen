from Bio.SeqUtils import MeltingTemp as MT
import numpy as np

from prigen.exceptions import NoPrimersGeneratedError
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

    @staticmethod
    def generate_primers(
            number_of_primers: int,
            length: int,
            gc_percentage: float,
            min_temperature: float,
            max_temperature: float,
            max_iterations: int = 10000
    ) -> dict[str, float]:
        """Generate nucleotide sequences with given GC-content and melting temperature

        Parameters
        ----------
        number_of_primers : int
            Number of nucleotide sequences to be generated
        length : int
            Length of the sequences
        gc_percentage : float
            Number between 0 and 1, a percent of G and C letters in primer
        min_temperature : float
            Minimal acceptable melting temperature
        max_temperature : float
            Maximal acceptable melting temperature
        max_iterations : int
            Maximum number of attempts to generate primer with given properties

        Returns
        -------
        primers : dict[str, float] -
            Dictionary, where keys are nucleotide sequences with desired properties,
            and values are their melting temperatures

        Raises
        ------
        ValueError
            If `min_temperature` > `max_temperature`, `gc_percentage` is not between 0 and 1, or
            if `length` < 0
        NoPrimersGeneratedError
            If `max_iterations` attempts were made, but no primers with desired properties
            were generated
        """
        # Give at least 100 attempts to generate each primer
        max_iterations = max(max_iterations, number_of_primers * 100)

        iteration = 0

        check_temperature_bounds(min_temperature, max_temperature)
        check_gc_percentage(gc_percentage)
        check_length(length)

        primers: dict[str, float] = {}
        temperatures = np.zeros(max_iterations, dtype="float")

        while len(primers) < number_of_primers and iteration < max_iterations:
            primer = PrimersGenerator.generate_primer(length, gc_percentage)

            if primer not in primers:
                melting_temperature = MT.Tm_NN(primer)
                temperatures[iteration] = melting_temperature

                if min_temperature <= melting_temperature <= max_temperature:
                    primers[primer] = melting_temperature

            iteration += 1

        if iteration == max_iterations and not primers:
            avg_temperature = round(np.mean(temperatures), 2)
            temperature_std = round(np.std(temperatures), 2)

            error_message = f"No primers were generated. Try to change parameters. "\
                            f"Average melting temperature of generated primers: {avg_temperature} "\
                            f"with standard deviation: {temperature_std}"

            raise NoPrimersGeneratedError(error_message)

        return primers

    def generate(self):
        pass
