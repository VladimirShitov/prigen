from prigen.validators import check_length, check_gc_percentage


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

        if self.number_of_primers < 0:
            raise ValueError("Number of primers must be greater than 0")

        if self.temperature_from >= self.temperature_to:
            raise ValueError(
                "The lower bound of melting temperature must be less than the upper bound"
            )

    def generate(self):
        pass
