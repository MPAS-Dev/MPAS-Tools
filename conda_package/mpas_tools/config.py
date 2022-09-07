from configparser import RawConfigParser, ConfigParser, ExtendedInterpolation
import os
from importlib import resources
import inspect
import sys
import numpy as np
import ast
from io import StringIO


class MpasConfigParser:
    """
    A "meta" config parser that keeps a dictionary of config parsers and their
    sources to combine when needed.  The custom config parser allows provenance
    of the source of different config options and allows the "user" config
    options to always take precedence over other config options (even if they
    are added later).

    Attributes
    ----------
    combined : {None, configparser.ConfigParser}
        The combined config options

    combined_comments : {None, dict}
        The combined comments associated with sections and options

    sources : {None, dict}
        The source of each section or option
    """

    _np_allowed = dict(linspace=np.linspace, xrange=range,
                       range=range, array=np.array, arange=np.arange,
                       pi=np.pi, Pi=np.pi, int=int, __builtins__=None)

    def __init__(self):
        """
        Make a new (empty) config parser
        """

        self._configs = dict()
        self._user_config = dict()
        self._comments = dict()
        self.combined = None
        self.combined_comments = None
        self.sources = None

    def add_user_config(self, filename):
        """
        Add a the contents of a user config file to the parser.  These options
        take precedence over all other options.

        Parameters
        ----------
        filename : str
            The relative or absolute path to the config file
        """
        self._add(filename, user=True)

    def add_from_file(self, filename):
        """
        Add the contents of a config file to the parser.

        Parameters
        ----------
        filename : str
            The relative or absolute path to the config file
        """
        self._add(filename, user=False)

    def add_from_package(self, package, config_filename, exception=True):
        """
        Add the contents of a config file to the parser.

        Parameters
        ----------
        package : str or Package
            The package where ``config_filename`` is found

        config_filename : str
            The name of the config file to add

        exception : bool, optional
            Whether to raise an exception if the config file isn't found
        """
        try:
            with resources.path(package, config_filename) as path:
                self._add(path, user=False)
        except (ModuleNotFoundError, FileNotFoundError, TypeError):
            if exception:
                raise

    def get(self, section, option):
        """
        Get an option value for a given section.

        Parameters
        ----------
        section : str
            The name of the config section

        option : str
            The name of the config option

        Returns
        -------
        value : str
            The value of the config option
        """
        if self.combined is None:
            self.combine()
        return self.combined.get(section, option)

    def getint(self, section, option):
        """
        Get an option integer value for a given section.

        Parameters
        ----------
        section : str
            The name of the config section

        option : str
            The name of the config option

        Returns
        -------
        value : int
            The value of the config option
        """
        if self.combined is None:
            self.combine()
        return self.combined.getint(section, option)

    def getfloat(self, section, option):
        """
        Get an option float value for a given section.

        Parameters
        ----------
        section : str
            The name of the config section

        option : str
            The name of the config option

        Returns
        -------
        value : float
            The value of the config option
        """
        if self.combined is None:
            self.combine()
        return self.combined.getfloat(section, option)

    def getboolean(self, section, option):
        """
        Get an option boolean value for a given section.

        Parameters
        ----------
        section : str
            The name of the config section

        option : str
            The name of the config option

        Returns
        -------
        value : bool
            The value of the config option
        """
        if self.combined is None:
            self.combine()
        return self.combined.getboolean(section, option)

    def getlist(self, section, option, dtype=str):
        """
        Get an option value as a list for a given section.

        Parameters
        ----------
        section : str
            The name of the config section

        option : str
            The name of the config option

        dtype : {Type[str], Type[int], Type[float]}
            The type of the elements in the list

        Returns
        -------
        value : list
            The value of the config option parsed into a list
        """
        values = self.get(section, option)
        values = [dtype(value) for value in values.replace(',', ' ').split()]
        return values

    def getexpression(self, section, option, dtype=None, use_numpyfunc=False):
        """
        Get an option as an expression (typically a list, though tuples and
        dicts are also available).  The expression is required to have valid
        python syntax, so that string entries are required to be in single or
        double quotes.

        Parameters
        ----------
        section : str
            The section in the config file

        option : str
            The option in the config file

        dtype : {Type[bool], Type[int], Type[float], Type[list], Type[tuple], Type[str]}, optional
            If supplied, each element in a list or tuple, or
            each value in a dictionary are cast to this type.  This is likely
            most useful for ensuring that all elements in a list of numbers are
            of type float, rather than int, when the distinction is important.

        use_numpyfunc : bool, optional
            If ``True``, the expression is evaluated including functionality
            from the numpy package (which can be referenced either as ``numpy``
            or ``np``).
        """

        expression_string = self.get(section, option)
        if use_numpyfunc:
            assert '__' not in expression_string, \
                f'"__" is not allowed in {expression_string} ' \
                f'for use_numpyfunc=True'
            sanitized_str = expression_string.replace('np.', '') \
                .replace('numpy.', '')
            result = eval(sanitized_str, MpasConfigParser._np_allowed)
        else:
            result = ast.literal_eval(expression_string)

        if dtype is not None:
            if isinstance(result, list):
                result = [dtype(element) for element in result]
            elif isinstance(result, tuple):
                result = (dtype(element) for element in result)
            elif isinstance(result, dict):
                for key in result:
                    result[key] = dtype(result[key])

        return result

    def has_section(self, section):
        """
        Whether the given section is part of the config

        Parameters
        ----------
        section : str
            The name of the config section

        Returns
        -------
        found : bool
            Whether the option was found in the section
        """
        if self.combined is None:
            self.combine()
        return self.combined.has_section(section)

    def has_option(self, section, option):
        """
        Whether the given section has the given option

        Parameters
        ----------
        section : str
            The name of the config section

        option : str
            The name of the config option

        Returns
        -------
        found : bool
            Whether the option was found in the section
        """
        if self.combined is None:
            self.combine()
        return self.combined.has_option(section, option)

    def set(self, section, option, value=None, comment=None, user=False):
        """
        Set the value of the given option in the given section.  The file from
         which this function was called is also retained for provenance.

        Parameters
        ----------
        section : str
            The name of the config section

        option : str
            The name of the config option

        value : str, optional
            The value to set the option to

        comment : str, optional
            A comment to include with the config option when it is written
            to a file

        user : bool, optional
            Whether this config option was supplied by the user (e.g. through
            a command-line flag) and should take priority over other sources
        """
        option = option.lower()
        calling_frame = inspect.stack(context=2)[1]
        filename = os.path.abspath(calling_frame.filename)

        if user:
            config_dict = self._user_config
        else:
            config_dict = self._configs
        if filename not in config_dict:
            config_dict[filename] = RawConfigParser()
        config = config_dict[filename]
        if not config.has_section(section):
            config.add_section(section)
        config.set(section, option, value)
        self.combined = None
        self.combined_comments = None
        self.sources = None
        if filename not in self._comments:
            self._comments[filename] = dict()
        if comment is None:
            comment = ''
        else:
            comment = ''.join([f'# {line}\n' for line in comment.split('\n')])
        self._comments[filename][(section, option)] = comment

    def write(self, fp, include_sources=True, include_comments=True):
        """
        Write the config options to the given file pointer.

        Parameters
        ----------
        fp : typing.TestIO
            The file pointer to write to.

        include_sources : bool, optional
            Whether to include a comment above each option indicating the
            source file where it was defined

        include_comments : bool, optional
            Whether to include the original comments associated with each
            section or option
        """
        if self.combined is None:
            self.combine()
        for section in self.combined.sections():
            section_items = self.combined.items(section=section)
            if include_comments and section in self.combined_comments:
                fp.write(self.combined_comments[section])
            fp.write(f'[{section}]\n\n')
            for option, value in section_items:
                if include_comments:
                    fp.write(self.combined_comments[(section, option)])
                if include_sources:
                    source = self.sources[(section, option)]
                    fp.write(f'# source: {source}\n')
                value = str(value).replace('\n', '\n\t').replace('$', '$$')
                fp.write(f'{option} = {value}\n\n')
            fp.write('\n')

    def list_files(self):
        """
        Get a list of files contributing to the combined config options

        Returns
        -------
        filenames : list of str
            A list of file paths

        """
        filenames = list(self._configs.keys()) + list(self._user_config.keys())
        return filenames

    def copy(self):
        """
        Get a deep copy of the config parser

        Returns
        -------
        config_copy : mpas_tools.config.MpasConfigParser
            The deep copy
        """
        config_copy = MpasConfigParser()
        for filename, config in self._configs.items():
            config_copy._configs[filename] = MpasConfigParser._deepcopy(config)

        for filename, config in self._user_config.items():
            config_copy._user_config[filename] = \
                MpasConfigParser._deepcopy(config)

        config_copy._comments = dict(self._comments)
        return config_copy

    def __getitem__(self, section):
        """
        Get get the config options for a given section.

        Parameters
        ----------
        section : str
            The name of the section to retrieve.

        Returns
        -------
        section_proxy : configparser.SectionProxy
            The config options for the given section.
        """
        if self.combined is None:
            self.combine()
        return self.combined[section]

    def _add(self, filename, user):
        filename = os.path.abspath(filename)
        config = RawConfigParser()
        if not os.path.exists(filename):
            raise FileNotFoundError(f'Config file does not exist: {filename}')
        config.read(filenames=filename)
        with open(filename) as fp:
            comments = self._parse_comments(fp, filename, comments_before=True)

        if user:
            self._user_config[filename] = config
        else:
            self._configs[filename] = config
        self._comments[filename] = comments
        self.combined = None
        self.combined_comments = None
        self.sources = None

    def combine(self):
        """
        Combine the config files into one.  This is normally handled
        automatically.
        """
        self.combined = ConfigParser(interpolation=ExtendedInterpolation())
        self.sources = dict()
        self.combined_comments = dict()
        for configs in [self._configs, self._user_config]:
            for source, config in configs.items():
                for section in config.sections():
                    if section in self._comments[source]:
                        self.combined_comments[section] = \
                            self._comments[source][section]
                    if not self.combined.has_section(section):
                        self.combined.add_section(section)
                    for option, value in config.items(section):
                        self.sources[(section, option)] = source
                        self.combined.set(section, option, value)
                        self.combined_comments[(section, option)] = \
                            self._comments[source][(section, option)]

    @staticmethod
    def _parse_comments(fp, filename, comments_before=True):
        """ Parse the comments in a config file into a dictionary """
        comments = dict()
        current_comment = ''
        section_name = None
        option_name = None
        indent_level = 0
        for line_number, line in enumerate(fp, start=1):
            value = line.strip()
            is_comment = value.startswith('#')
            if is_comment:
                current_comment = current_comment + line
            if len(value) == 0 or is_comment:
                # end of value
                indent_level = sys.maxsize
                continue

            cur_indent_level = len(line) - len(line.lstrip())
            is_continuation = cur_indent_level > indent_level
            # a section header or option header?
            if section_name is None or option_name is None or \
                    not is_continuation:
                indent_level = cur_indent_level
                # is it a section header?
                is_section = value.startswith('[') and value.endswith(']')
                if is_section:
                    if not comments_before:
                        if option_name is None:
                            comments[section_name] = current_comment
                        else:
                            comments[(section_name, option_name)] = \
                                current_comment
                    section_name = value[1:-1].strip()
                    option_name = None

                    if comments_before:
                        comments[section_name] = current_comment
                    current_comment = ''
                # an option line?
                else:
                    delimiter_index = value.find('=')
                    if delimiter_index == -1:
                        raise ValueError(f'Expected to find "=" on line '
                                         f'{line_number} of {filename}')

                    if not comments_before:
                        if option_name is None:
                            comments[section_name] = current_comment
                        else:
                            comments[(section_name, option_name)] = \
                                current_comment

                    option_name = value[:delimiter_index].strip().lower()

                    if comments_before:
                        comments[(section_name, option_name)] = current_comment
                    current_comment = ''

        return comments

    @staticmethod
    def _deepcopy(config):
        """ Make a deep copy of the ConfigParser object """
        config_string = StringIO()
        config.write(config_string)
        # We must reset the buffer to make it ready for reading.
        config_string.seek(0)
        new_config = ConfigParser()
        new_config.read_file(config_string)
        return new_config
