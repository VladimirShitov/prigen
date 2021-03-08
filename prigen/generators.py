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
        if self.length < 0:
            raise ValueError("Length must be greater than 0")

        if self.gc_percentage < 0 or self.gc_percentage > 1:
            raise ValueError("GC-percentage must be a number between 0 and 1")

        if self.number_of_primers < 0:
            raise ValueError("Number of primers must be greater than 0")

        if self.temperature_from >= self.temperature_to:
            raise ValueError(
                "The lower bound of melting temperature must be less than the upper bound"
            )

    def generate(self):
        pass
