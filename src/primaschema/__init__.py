"""Infrastructure for tiled amplicon PCR primer scheme definitions"""

import logging
import logging.config
import os

from pathlib import Path


__version__ = "1.0.0a0"

pkg_dir = Path(
    os.environ.get(
        "PRIMASCHEMA_ROOT_PATH", Path(__file__).absolute().parent.parent.parent
    )
)
schema_dir = pkg_dir / "src" / "primaschema" / "schema"
manifest_schema_path = schema_dir / "manifest.json"
header_path = schema_dir / "manifest-header.yml"

logging_config = {
    "version": 1,
    "disable_existing_loggers": False,
    "formatters": {
        "default": {
            "format": "%(name)s %(levelname)s: %(message)s",
        },
    },
    "handlers": {
        "console": {
            "class": "logging.StreamHandler",
            "formatter": "default",
            "level": "INFO",
            "stream": "ext://sys.stderr",
        },
    },
    "root": {
        "handlers": ["console"],
        "level": "INFO",
    },
    "loggers": {
        "primaschema": {
            "handlers": ["console"],
            "level": "INFO",
            "propagate": False,
        },
    },
}

logging.config.dictConfig(logging_config)
logger = logging.getLogger("primaschema")

logger.debug(f"{pkg_dir=} {schema_dir=}")
