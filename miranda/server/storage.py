"""
=====================
Disk space management
=====================

Classes:

 * DiskSpaceError - the exception raised on failure.
 * :class:`FileMeta` - file and its size.
 * :class:`StorageState` - storage capacity and availability of a medium.

Functions:

 * :func:`total_size` - get total size of a list of files.
 * :func:`size_division` - divide files based on number and size restrictions.

"""
import os


class DiskSpaceError(Exception):
    pass


class FileMeta:
    """File path and size."""

    django = {
        "path": ["CharField", "max_length=512"],
        "size": ["IntegerField", "null=True", "blank=True"],
    }

    def __init__(self, path, size=-1):
        """Initialize file meta.

        Parameters
        ----------
        path : str
            full path of the file.
        size : int
            size of file in bytes (default: will obtain from os.path.getsize
            if file exists, set to 0 otherwise).

        """

        # Make sure we have the full path of the file
        self.path = os.path.abspath(path)

        # Get size of file if it is not specified
        if (-1 == size) and os.path.isfile(self.path):
            try:
                self.size = os.path.getsize(self.path)
            except OSError:
                raise DiskSpaceError("Cannot get size of " + self.path + ".")
        elif -1 == size:
            self.size = 0
        else:
            self.size = size

    def __eq__(self, other):
        if self.path == other.path:
            return True
        else:
            return False


#
class StorageState:
    """Information regarding the storage capacity of a disk."""

    def __init__(self, base_path, capacity=-1, used_space=-1, free_space=-1):
        """Initialize storage state.

        Parameters
        ----------
        base_path : str
            base path of the storage medium.
        capacity : int
            capacity of medium in bytes (default: will obtain from system
            call to 'df').
        used_space : int
            space currently used on the medium (default: will obtain from
            system call to 'df').
        free_space : int
            space available on the medium (default: will obtain from system
            call to 'df').

        """

        # Make sure we have the full base path
        self.base_path = os.path.abspath(base_path)

        # Get attributes from 'df' function if they are not specified
        if (-1 == capacity) or (-1 == used_space) or (-1 == free_space):
            if not os.path.isdir(base_path):
                raise DiskSpaceError("Cannot analyse " + base_path + ".")
            elif not os.path.isfile("/bin/df"):
                raise DiskSpaceError("/bin/df does not exist.")
            df_output = os.popen("/bin/df -P " + base_path + " | tail -1").read()
            if not df_output:
                raise DiskSpaceError("/bin/df call failed.")
            df_output_split = df_output.split()
        if -1 == capacity:
            try:
                self.capacity = int(df_output_split[1]) * 1000
            except:
                raise DiskSpaceError("/bin/df output not as expected.")
        else:
            self.capacity = capacity
        if -1 == used_space:
            try:
                self.used_space = int(df_output_split[2]) * 1000
            except:
                raise DiskSpaceError("/bin/df output not as expected.")
        else:
            self.used_space = used_space
        if -1 == free_space:
            try:
                self.free_space = int(df_output_split[3]) * 1000
            except:
                raise DiskSpaceError("/bin/df output not as expected.")
        else:
            self.free_space = free_space


#
def total_size(file_list):
    """Total size of files.

    Parameters
    ----------
    file_list : list of strings or FileMeta objects

    Returns
    -------
    out : int
        total size of files in bytes.

    """

    if file_list:
        size = 0
        for file_to_add in file_list:
            # If file paths are given, convert to FileMeta objects first
            if not isinstance(file_to_add, FileMeta):
                try:
                    file_to_add = FileMeta(file_to_add)
                except DiskSpaceError:
                    raise
            size = size + file_to_add.size
        return size
    else:
        return 0


#
def size_division(
    files_to_divide,
    size_limit=0,
    file_limit=0,
    check_name_repetition=False,
    preserve_order=False,
):
    """Divide files according to size and number limits.

    Parameters
    ----------
    files_to_divide : list of strings or FileMeta objects
    size_limit : int
        size limit of divisions in bytes (default: no limit).
    file_limit : int
        number of files limit of divisions (default: no limit).
    check_name_repetition : bool
        flag to prevent file name repetitions (default: off).
    preserve_order : bool
        flag to force files to be restored in the order they are given
        (default: off).

    Returns
    -------
    out - 2d nested lists
        list of divisions (each division is a list of FileMeta objects).

    """

    divisions = []
    for file_divide in files_to_divide:
        # If file paths are given, convert to FileMeta objects first
        if not isinstance(file_divide, FileMeta):
            try:
                file_divide = FileMeta(file_divide)
            except DiskSpaceError:
                raise
        # Loop through divisions and try to add file according to limitations
        for i, division in enumerate(divisions):
            size = file_divide.size
            file_count = 1
            flag_skip = 0
            for file_divided in division:
                if check_name_repetition and (
                    os.path.basename(file_divided.path)
                    == os.path.basename(file_divide.path)
                ):
                    flag_skip = 1
                size = size + file_divided.size
                file_count = file_count + 1
            if (
                (size > size_limit != 0)
                or (file_count > file_limit != 0)
                or flag_skip == 1
            ):
                continue
            elif preserve_order and (i != len(divisions) - 1):
                continue
            else:
                divisions[i].append(file_divide)
                break
        else:
            divisions.append([file_divide])
    return divisions


#
