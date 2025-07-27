"""Helper functions for the library"""

from time import perf_counter
import logging

logger = logging.getLogger(__name__)


class SimpleTimer:
    """
    A simple context manager to time a block of code.

    This context manager measures the execution time of the code within a
    `with` block. After the block finishes, the total duration in seconds is
    available in the `.duration` attribute.

    Attributes
    ----------
    duration : float
        The total execution time of the block in seconds.

    Example
    -------
    >>> timer = SimpleTimer()
    >>> with timer:
    ...     time.sleep(1)
    >>> print(f"Operation took {timer.duration:.2f} seconds.")
    Operation took 1.00 seconds.
    """
    def __init__(self):
        self.start: float = 0.0
        self.end: float = 0.0
        self.duration: float = 0.0

    def __enter__(self):
        self.start = perf_counter()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.end = perf_counter()
        self.duration = self.end - self.start
