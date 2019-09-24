import logging
import os
import platform
from contextlib import contextmanager
from datetime import date
from datetime import datetime as dt
from pathlib import Path
from types import GeneratorType
from typing import Iterable
from typing import List
from typing import Sequence
from typing import Union

KiB = int(pow(2, 10))
MiB = int(pow(2, 20))
GiB = int(pow(2, 30))


def creation_date(path_to_file: Union[Path, str]) -> float or date:
    """
    Try to get the date that a file was created, falling back to when it was
    last modified if that isn't possible.
    See http://stackoverflow.com/a/39501288/1709587 for explanation.
    """
    if platform.system() == "Windows":
        return Path(path_to_file).stat().st_ctime
    else:
        stat = Path(path_to_file).stat()
        try:
            return date.fromtimestamp(stat.st_ctime)
        except AttributeError:
            # We're probably on Linux. No easy way to get creation dates here,
            # so we'll settle for when its content was last modified.
            return date.fromtimestamp(stat.st_mtime)


def read_privileges(location: Union[Path, str]) -> bool:
    if os.access(location, os.R_OK):
        logging.info(
            "{}: {} Read OK!".format(dt.now().strftime("%Y-%m-%d %X"), location)
        )
    else:
        msg = "Ensure read privileges. Exiting."
        logging.error(msg)
        raise Exception(msg)
    return True


@contextmanager
def working_directory(directory):
    """
    This function momentarily changes the working directory within the
     context and reverts to the file working directory when the code block
     it is acting upon exits
    """
    owd = os.getcwd()
    try:
        os.chdir(directory)
        yield directory
    finally:
        os.chdir(owd)
    return


def find_filepaths(
    source: Union[Path, str, GeneratorType, List[Union[Path, str]]],
    recursive: bool = True,
    file_suffixes: List[str] = None,
    **_
) -> List[Path]:

    if file_suffixes is None:
        file_suffixes = list().append(["*", ".*"])

    file_list = list()
    if isinstance(source, (GeneratorType, List)):
        file_list = [Path(f).expanduser() for f in source]
    else:
        for pattern in file_suffixes:
            if recursive:
                found = [f for f in Path(source).expanduser().rglob(pattern)]
            elif not recursive:
                found = [f for f in Path(source).expanduser().glob(pattern)]
            else:
                raise ValueError("Recursive: {}".format(recursive))
            file_list.append(found)
    return file_list


def single_item_list(iterable: Iterable) -> bool:
    """
    See: https://stackoverflow.com/a/16801605/7322852
    """
    iterator = iter(iterable)
    has_true = any(iterator)  # consume from "i" until first true or it's exhausted
    has_another_true = any(
        iterator
    )  # carry on consuming until another true value / exhausted
    return has_true and not has_another_true  # True if exactly one true found


def make_local_dirs(pathway):
    """Create directories recursively, unless they already exist.
    Parameters
    ----------
    pathway : Union[Path, str]
      Path of folders to create.
    """

    pathway = Path(pathway)

    if not pathway.exists():
        try:
            pathway.mkdir(parents=True)
        except OSError:
            raise


def set_comparions(set1: Sequence, set2: Sequence) -> bool:
    """Compare two sequences of non hashable objects as if they were sets.

    Parameters
    ----------
    set1 : Sequence
      First sequence of objects.
    set2 : Sequence
      Second sequence of objects.
    Returns
    -------
    out : bool
      True if two sets are identical, that is, they contain the same elements.
    """

    for item1 in set1:
        if item1 not in set2:
            return False
    for item2 in set2:
        if item2 not in set1:
            return False
    return True


########################################################################################


def yesno_prompt(query):
    """Prompt user for a yes/no answer.
    Parameters
    ----------
    query : str
        the yes/no question to ask the user.

    Returns
    -------
    out : bool
        True (yes) or False (otherwise).
    """

    user_input = input("{} (y/n) ".format(query))
    if user_input.lower() == "y":
        return True
    elif user_input.lower() == "n":
        return False
    else:
        raise ValueError("{} not in (y, n)".format(user_input))


def verbose_fn(message, verbose=True):
    """Trigger verbose mode.
    Parameters
    ----------
    message : str
    verbose : bool
        flag for whether of not to output the message (default: True).
    """

    if verbose:
        print(message)


def list_paths_with_elements(base_paths: List[str], elements: List[str]):
    """List a given path structure.
    Parameters
    ----------
    base_paths : List[str]
        list of paths from which to start the search.
    elements : List[str]
        ordered list of the expected elements.
    Returns
    -------
    out : list of dictionaries
        the keys are 'path' and each of the members of the given elements,
        the path is the absolute path.
    Notes
    -----
    Suppose you have the following structure:
    /base_path/{color}/{shape}
    The resulting list would look like:
    [{'path':/base_path/red/square, 'color':'red', 'shape':'square'},
     {'path':/base_path/red/circle, 'color':'red', 'shape':'circle'},
     {'path':/base_path/blue/triangle, 'color':'blue', 'shape':'triangle'},
     ...
    ]
    Obviously, 'path' should not be in the input list of elements.
    """

    # Make sure the base_paths input is a list of absolute path
    if not hasattr(base_paths, "__iter__"):
        base_paths = [base_paths]
    base_paths = map(os.path.abspath, base_paths)
    # If elements list is empty, return empty list (end of recursion).
    if not elements:
        return []
    #
    paths_elements = []
    for base_path in base_paths:
        try:
            path_content = os.listdir(base_path)
        except NotADirectoryError:
            continue
        path_content.sort()
        next_base_paths = []
        for path_item in path_content:
            next_base_paths.append(os.path.join(base_path, path_item))
        next_pe = list_paths_with_elements(next_base_paths, elements[1:])
        if next_pe:
            for i, one_pe in enumerate(next_pe):
                relative_path = next_pe[i]["path"].replace(base_path, "", 1)
                new_element = relative_path.split("/")[1]
                next_pe[i][elements[0]] = new_element
            paths_elements.extend(next_pe)
        elif len(elements) == 1:
            for my_path, my_item in zip(next_base_paths, path_content):
                paths_elements.append({"path": my_path, elements[0]: my_item})
    return paths_elements


def get_info_var(variable_name):
    """fonction qui retourne differentes informations en fonction de la variable voulue"""

    if variable_name == "tas":
        var_code = 78
        unites = "degC"
    elif variable_name == "hourly_rainfall":
        var_code = 123
        unites = "mm"
    elif variable_name == "precipitation_amount":
        var_code = 262
        unites = "mm"
    else:
        raise RuntimeError
    fact_mlt = 0.1
    fact_add = 0.0
    missing_flags = ["M"]
    least_sig_digit = None

    return dict(
        code_var=var_code,
        unites=unites,
        fact_mlt=fact_mlt,
        fact_add=fact_add,
        flag_manquants=missing_flags,
        least_significant_digit=least_sig_digit,
    )
