"""Internal module: logging configuration for VarSim.

Provides a consistent logging setup across all package modules.
Import ``get_logger(__name__)`` in each module to obtain a logger.
"""

import logging
import time
from functools import wraps
from contextlib import contextmanager

# Module-level default: WARNING — users must opt in to see INFO/DEBUG
_DEFAULT_LEVEL = logging.WARNING

_handler: logging.Handler | None = None


def _ensure_handler() -> logging.Handler:
    """Ensure there is exactly one StreamHandler for the 'varsim' logger."""
    global _handler
    root = logging.getLogger("varsim")
    if _handler is None:
        _handler = logging.StreamHandler()
        _handler.setFormatter(logging.Formatter(
            "[%(asctime)s] %(levelname)-7s %(name)s | %(message)s",
            datefmt="%H:%M:%S",
        ))
        _handler.setLevel(logging.DEBUG)  # handler passes everything; logger level controls
        root.addHandler(_handler)
        root.setLevel(_DEFAULT_LEVEL)
        root.propagate = False  # don't bubble to root logger
    return _handler


def set_log_level(level: int | str) -> None:
    """Set the logging level for all VarSim loggers.

    Parameters
    ----------
    level : int or str
        One of ``logging.DEBUG``, ``logging.INFO``, ``logging.WARNING``,
        ``logging.ERROR``, or the equivalent string like ``"DEBUG"``.

    Examples
    --------
    >>> from varsim._logging import set_log_level
    >>> import logging
    >>> set_log_level(logging.DEBUG)   # show everything
    >>> set_log_level("INFO")          # show info and above
    >>> set_log_level(logging.WARNING) # default — warnings and errors only
    """
    _ensure_handler()
    if isinstance(level, str):
        level = getattr(logging, level.upper(), logging.WARNING)
    logging.getLogger("varsim").setLevel(level)


def get_logger(name: str) -> logging.Logger:
    """Get a logger for a VarSim sub-module.

    Parameters
    ----------
    name : str
        Typically ``__name__`` from the calling module.

    Returns
    -------
    logging.Logger
        A child logger of ``"varsim"`` with a shared handler and formatter.
    """
    _ensure_handler()
    return logging.getLogger(name)


def log_duration(logger: logging.Logger, message: str, level: int = logging.DEBUG):
    """Context manager that logs the duration of a code block.

    Usage::

        with log_duration(logger, "Fetching G6PD"):
            result = nm("G6PD")

    Logs ``message ... done (1.23s)`` at *level* on exit.
    """
    @contextmanager
    def _ctx():
        t0 = time.perf_counter()
        yield
        dt = time.perf_counter() - t0
        logger.log(level, "%s ... done (%.2fs)", message, dt)
    return _ctx()
