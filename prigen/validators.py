def check_length(length: int):
    if length <= 0:
        raise ValueError("Length must be greater than 0")


def check_gc_percentage(gc_percentage):
    if gc_percentage < 0 or gc_percentage > 1:
        raise ValueError("GC-percentage must be a number between 0 and 1")

