# Run with, for example: "pytest -v"

from pathlib import Path

import pytest

from locked_dataset import LockedDataset


def test_creation(tmp_path: Path) -> None:
    fname = tmp_path / "test.nc"
    lock = tmp_path / "test.nc.lock"

    # Neither the file nor the lock file should be present
    assert not fname.exists()
    assert not lock.exists()

    att_val = 1
    with LockedDataset(fname, "w") as f:
        f.setncattr("test", att_val)
        assert lock.exists()

    # At this point the lock file should be gone
    assert fname.exists()
    assert not lock.exists()

    with LockedDataset(fname) as f:
        assert f.getncattr("test") == att_val


def test_concurrent(tmp_path: Path) -> None:
    fname = tmp_path / "test.nc"
    assert not fname.exists()

    # Create the file
    with LockedDataset(fname, "w") as f:
        f.setncattr("test", 1)

        # While the file is opened/locked, check that a second case can't open
        # it
        with pytest.raises(FileExistsError):
            with LockedDataset(fname, "a", -1) as f2:
                # This should never be reached
                f2.setncattr("test", 2)

    with LockedDataset(fname) as f:
        assert f.getncattr("test", 1)
