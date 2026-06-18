"""Shared logging configuration for IFCB command-line workflows."""

from __future__ import annotations

from datetime import datetime
import logging
from pathlib import Path
import sys
from typing import Mapping

SECRET_TERMS = ("password", "token", "secret", "api_key", "apikey")


def _format_config_value(value: object) -> str:
    """Format one workflow setting for a compact, readable log entry."""
    if isinstance(value, (list, tuple)):
        return ", ".join(str(item) for item in value)
    return str(value)


def redact_command_line(argv: list[str]) -> str:
    """Return a command string with values for secret-like options redacted."""
    redacted: list[str] = []
    redact_next = False
    for argument in argv:
        if redact_next:
            redacted.append("<redacted>")
            redact_next = False
            continue

        option, separator, value = argument.partition("=")
        is_secret = any(term in option.lower() for term in SECRET_TERMS)
        if is_secret and separator:
            redacted.append(f"{option}=<redacted>")
        else:
            redacted.append(argument)
            redact_next = is_secret
    return " ".join(redacted)


def log_run_configuration(logger: logging.Logger, settings: Mapping[str, object]) -> None:
    """Log workflow inputs while redacting values with secret-like names."""
    logger.info("Run configuration:")
    for key, value in settings.items():
        safe_value = "<redacted>" if any(term in key.lower() for term in SECRET_TERMS) else value
        logger.info("  %s: %s", key, _format_config_value(safe_value))


def setup_logging(
    log_dir: str | Path,
    name: str,
    level: int = logging.INFO,
    console: bool = True,
    file: bool = True,
) -> logging.Logger:
    """Configure timestamped file logging and optional console output."""
    log_dir = Path(log_dir)
    log_dir.mkdir(parents=True, exist_ok=True)

    logger = logging.getLogger()
    logger.setLevel(level)

    # Replace handlers created by an earlier IFCB workflow while preserving
    # handlers owned by embedding applications such as notebooks.
    for handler in list(logger.handlers):
        if getattr(handler, "_ifcb_handler", False):
            logger.removeHandler(handler)
            handler.close()

    formatter = logging.Formatter("%(asctime)s | %(levelname)s | %(name)s | %(message)s")
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

    if file:
        out_path = log_dir / f"{name}_{timestamp}.out.log"
        err_path = log_dir / f"{name}_{timestamp}.err.log"

        out_handler = logging.FileHandler(out_path)
        out_handler._ifcb_handler = True
        out_handler.setLevel(level)
        out_handler.setFormatter(formatter)

        err_handler = logging.FileHandler(err_path)
        err_handler._ifcb_handler = True
        err_handler.setLevel(logging.ERROR)
        err_handler.setFormatter(formatter)

        logger.addHandler(out_handler)
        logger.addHandler(err_handler)

    if console:
        console_handler = logging.StreamHandler()
        console_handler._ifcb_handler = True
        console_handler.setLevel(level)
        console_handler.setFormatter(formatter)
        logger.addHandler(console_handler)

    if file:
        logger.info("Logging to: %s", out_path)
        logger.info("Errors to: %s", err_path)

    # Command-line workflows may fail outside an explicit try/except block.
    # Route those fatal exceptions through logging so the error file remains a
    # complete diagnostic record while retaining normal KeyboardInterrupt behavior.
    previous_hook = sys.excepthook

    def log_uncaught_exception(exc_type, exc_value, exc_traceback):
        if issubclass(exc_type, KeyboardInterrupt):
            previous_hook(exc_type, exc_value, exc_traceback)
            return
        logging.getLogger(name).critical(
            "Unhandled exception",
            exc_info=(exc_type, exc_value, exc_traceback),
        )

    sys.excepthook = log_uncaught_exception

    return logger
