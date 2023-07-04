import logging

try:
    import bayesflow
except ImportError as err:
    logger = logging.getLogger()
    logger.setLevel(logging.WARNING)
    logger.warn(str(err))
    logger.warn("Some dependencies failed to load.")

# from .exg2 import prior, simulator

__version__ = "0.0.1"
__author__ = "Leonhard Volz"
__email__ = "l.volz@uva.nl"
