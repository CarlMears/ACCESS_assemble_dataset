from os import PathLike
from pathlib import Path
from time import sleep
from typing import Literal, Union

from netCDF4 import Dataset
from typing_extensions import TypeAlias


StrPath: TypeAlias = Union[str, PathLike[str]]


class LockedDataset:
    """Simple wrapper to lock a netCDF file using an external lock file.

    >>> with LockedDataset("file.nc", "w") as f:
    ...     f.setncattr("test", 1)
    """

    def __init__(
        self, filename: StrPath, mode: Literal["r", "w", "a"] = "a", timeout: float = 30
    ) -> None:
        """Wrap the netCDF file with a lock.

        filename: the netCDF file to read/write

        mode: the mode string passed to netCDF

        timeout: the number of seconds to wait before retries. If it's a
        negative value then the lock will fail immediately instead of retrying.
        """

        self.filename = Path(filename)
        self.mode = mode
        self.timeout = timeout
        self.lock = self.filename.with_suffix(".nc.lock")

    def __enter__(self) -> Dataset:
        while True:
            try:
                self.lock.touch(exist_ok=False)
            except FileExistsError as e:
                if self.timeout >= 0:
                    print(f'collision! - waiting {self.timeout} seconds')
                    sleep(self.timeout)
                else:
                    raise e
            else:
                break

        self._dataset = Dataset(self.filename, self.mode)
        return self._dataset

    def __exit__(self, _exc_type, _exc_value, _traceback) -> None:
        self._dataset.close()
        self.lock.unlink()
