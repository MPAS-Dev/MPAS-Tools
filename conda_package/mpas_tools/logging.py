import sys
import logging
import subprocess


def check_call(args, logger, log_command=True, **kwargs):
    """
    Call the given subprocess with logging to the given logger.

    Parameters
    ----------
    args : list or str
        A list or string of argument to the subprocess.  If ``args`` is a
        string, you must pass ``shell=True`` as one of the ``kwargs``.

    logger : logging.Logger
        The logger to write output to

    log_command : bool, optional
        Whether to write the command that is running ot the logger

    **kwargs : dict
        Keyword arguments to pass to subprocess.Popen

    Raises
    ------
    subprocess.CalledProcessError
        If the given subprocess exists with nonzero status

    """

    if isinstance(args, str):
        print_args = args
    else:
        print_args = ' '.join(args)
    if log_command:
        logger.info(f'Running: {print_args}')

    process = subprocess.Popen(args, stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE, **kwargs)
    stdout, stderr = process.communicate()

    if stdout:
        stdout = stdout.decode('utf-8')
        for line in stdout.split('\n'):
            logger.info(line)
    if stderr:
        stderr = stderr.decode('utf-8')
        for line in stderr.split('\n'):
            logger.error(line)

    if process.returncode != 0:
        raise subprocess.CalledProcessError(process.returncode,
                                            print_args)


class LoggingContext(object):

    """
    A context manager for creating a logger or using an existing logger

    Attributes
    ----------
    logger : logging.Logger
        A logger that sends output to a log file or stdout/stderr
    """

    def __init__(self, name, logger=None, log_filename=None):
        """
        If ``logger`` is ``None``, create a new logger either to a log file
        or stdout/stderr.  If ``logger`` is anything else, just set the logger
        attribute

        Parameters
        ----------
        name : str
            A unique name for the logger (e.g. ``__name__`` of the calling
            module)

        logger : logging.Logger, optional
           An existing logger that sends output to a log file or stdout/stderr
           to be used in this context

        log_filename : str, optional
            The name of a file where output should be written.  If none is
            supplied, output goes to stdout/stderr
        """
        self.logger = logger
        self.name = name
        self.log_filename = log_filename
        self.handler = None
        self.old_stdout = None
        self.old_stderr = None
        self.existing_logger = logger is not None

    def __enter__(self):
        if not self.existing_logger:
            if self.log_filename is not None:
                # redirect output to a log file
                logger = logging.getLogger(self.name)
                handler = logging.FileHandler(self.log_filename)
            else:
                logger = logging.getLogger(self.name)
                handler = logging.StreamHandler(sys.stdout)

            formatter = MpasFormatter()
            handler.setFormatter(formatter)
            logger.addHandler(handler)
            logger.setLevel(logging.INFO)
            logger.propagate = False
            self.logger = logger
            self.handler = handler

            if self.log_filename is not None:
                self.old_stdout = sys.stdout
                self.old_stderr = sys.stderr
                sys.stdout = StreamToLogger(logger, logging.INFO)
                sys.stderr = StreamToLogger(logger, logging.ERROR)
        return self.logger

    def __exit__(self, exc_type, exc_val, exc_tb):
        if not self.existing_logger:
            if self.old_stdout is not None:
                self.handler.close()
                # restore stdout and stderr
                sys.stdout = self.old_stdout
                sys.stderr = self.old_stderr

            # remove the handlers from the logger (probably only necessary if
            # writeLogFile==False)
            self.logger.handlers = []

        self.stdout = self.original_stdout = sys.stdout
        self.stderr = self.original_stderr = sys.stderr


class MpasFormatter(logging.Formatter):
    """
    A custom formatter for logging
    Modified from:
    https://stackoverflow.com/a/8349076/7728169
    https://stackoverflow.com/a/14859558/7728169
    """
    # Authors
    # -------
    # Xylar Asay-Davis

    # printing error messages without a prefix because they are sometimes
    # errors and sometimes only warnings sent to stderr
    dbg_fmt = "DEBUG: %(module)s: %(lineno)d: %(msg)s"
    info_fmt = "%(msg)s"
    err_fmt = info_fmt

    def __init__(self, fmt=info_fmt):
        logging.Formatter.__init__(self, fmt)

    def format(self, record):

        # Save the original format configured by the user
        # when the logger formatter was instantiated
        format_orig = self._fmt

        # Replace the original format with one customized by logging level
        if record.levelno == logging.DEBUG:
            self._fmt = MpasFormatter.dbg_fmt

        elif record.levelno == logging.INFO:
            self._fmt = MpasFormatter.info_fmt

        elif record.levelno == logging.ERROR:
            self._fmt = MpasFormatter.err_fmt

        # Call the original formatter class to do the grunt work
        result = logging.Formatter.format(self, record)

        # Restore the original format configured by the user
        self._fmt = format_orig

        return result


class StreamToLogger(object):
    """
    Modified based on code by:
    https://www.electricmonk.nl/log/2011/08/14/redirect-stdout-and-stderr-to-a-logger-in-python/
    Copyright (C) 2011 Ferry Boender
    License: GPL, see https://www.electricmonk.nl/log/posting-license/
    Fake file-like stream object that redirects writes to a logger instance.
    """

    def __init__(self, logger, log_level=logging.INFO):
        self.logger = logger
        self.log_level = log_level
        self.linebuf = ''

    def write(self, buf):
        for line in buf.rstrip().splitlines():
            self.logger.log(self.log_level, str(line.rstrip()))

    def flush(self):
        pass
