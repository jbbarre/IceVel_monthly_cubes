"""
logger.py

Author: Jean Baptiste Barré
Created: 2025-04-18
Description:
    logger utility module.
    Provides a centralized logger configuration for the whole project.
"""

import logging
import sys



def get_logger(name: str = "IceVel_monthlycubes") -> logging.Logger:
    """Create and configure a logger instance.

    Args:
        name (str): Name of the logger. Default is "IceVel_monthlycubes".
    Returns:
        logging.logger: Configured logger instance.
    """
    logger = logging.getLogger(name)

    if not logger.hasHandlers():
        logger.setLevel(logging.DEBUG)
        # Create a console handler
        console_handler = logging.StreamHandler(sys.stdout)
        console_handler.setLevel(logging.INFO)

        # create a formatter
        formatter = logging.Formatter(
            "%(asctime)s - %(name)s - %(levelname)s - %(message)s",
            datefmt="%Y-%m-%d %H:%M:%S",
        )

        console_handler.setFormatter(formatter)
        # Add the handler to the logger
        logger.addHandler(console_handler)

        # File handler (DEBUG level)
        file_handler = logging.FileHandler(f"{name.lower()}.log")
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)

    return logger
