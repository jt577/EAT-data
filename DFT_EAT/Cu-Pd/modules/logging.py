import os

def log(message, log_file):
    """
    Logs a message to the log file.

    Args:
        message (str): The message to log.
        log_file (str): The path to the log file.
    """
    with open(log_file, 'a') as f:
        f.write(message + '\n')