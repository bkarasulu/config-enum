#!/usr/bin/env python
"""config_enum command-line script (formerly config_enum-v2.py).

Equivalent to the ``config_enum`` command installed with ``pip install .``; all command-line options are unchanged.

This file has the same name as the ``config_enum`` package. Run from this folder, Python would find this file first for
``import config_enum`` (e.g. ``python -m pytest`` or ``python -c "import config_enum"``); it then loads the real package
in its place, so both uses work.
"""
import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))


def _load_package():
    """Import the config_enum package (installed, or src/ next to this script), skipping this file."""
    saved = list(sys.path)
    sys.path[:] = [p for p in sys.path if os.path.abspath(p or os.curdir) != _HERE]
    shadow = sys.modules.pop("config_enum", None)
    try:
        try:
            import config_enum as pkg
        except ImportError:  # not installed: use the source tree next to this script
            sys.path.insert(0, os.path.join(_HERE, "src"))
            import config_enum as pkg
    except Exception:
        if shadow is not None:
            sys.modules["config_enum"] = shadow
        raise
    finally:
        sys.path[:] = saved
    sys.modules["config_enum"] = pkg
    return pkg


if __name__ == "__main__":
    _load_package()
    from config_enum.cli import main
    sys.exit(main())
else:  # imported as "config_enum" by mistake: replace this module by the package
    _load_package()
