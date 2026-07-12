# logging_utils.py
from __future__ import annotations

import logging
import sys
from datetime import datetime
from pathlib import Path


def setup_logging(out_dir: Path, *, name: str = "optimization") -> logging.Logger:
    """
    Create a logger that writes:
      - <out_dir>/logs/<name>_<timestamp>.log
      - stdout
    """
    out_dir = Path(out_dir)
    log_dir = out_dir / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_file = log_dir / f"{name}_{timestamp}.log"

    logger = logging.getLogger(name)
    logger.setLevel(logging.INFO)
    logger.propagate = False  # prevent double-printing via root logger

    # If this logger already has handlers (re-runs), don’t add duplicates.
    if not logger.handlers:
        fmt = logging.Formatter("%(asctime)s | %(message)s")

        fh = logging.FileHandler(log_file)
        fh.setLevel(logging.INFO)
        fh.setFormatter(fmt)

        sh = logging.StreamHandler(sys.stdout)
        sh.setLevel(logging.INFO)
        sh.setFormatter(fmt)

        logger.addHandler(fh)
        logger.addHandler(sh)

    logger.info("🚀 Combo Optimization started")
    logger.info(f"Log file: {log_file}")
    return logger
