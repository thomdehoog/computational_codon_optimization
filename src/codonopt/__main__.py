"""Allow ``python -m codonopt``."""

import sys

from .cli import main

sys.exit(main())
