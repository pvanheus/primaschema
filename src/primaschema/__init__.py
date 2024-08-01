"""A toolkit for tiling primer scheme definitions"""

import logging
import logging.config
import os

from pathlib import Path


__version__ = "1.0.0a0"

pkg_dir = Path(
    os.environ.get(
        "PRIMER_SCHEMES_PATH", Path(__file__).absolute().parent.parent.parent
    )
)
schema_dir = pkg_dir / "test" / "data" / "primer-schemes" / "schema"
info_schema_path = schema_dir / "info-schema.yml"
manifest_schema_path = schema_dir / "manifest-schema.json"
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
