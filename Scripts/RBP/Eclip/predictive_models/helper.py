import logging


def get_logger(name, logger_path=None):
    """Setup logger
    """
    logger = logging.getLogger(name)

    formatter = logging.Formatter('%(levelname)s:%(asctime)s:%(name)s] %(message)s')
    # file logger
    if logger_path is None:
        fh = logging.FileHandler(logger_path)
        fh.setLevel(logging.DEBUG)
        fh.setFormatter(formatter)
        logger.addHandler(fh)

    # console logger
    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)
    ch.setFormatter(formatter)
    logger.addHandler(ch)

    return logger
