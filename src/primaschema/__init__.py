"""A toolkit for tiling primer scheme definitions"""

import logging
import logging.config

from pathlib import Path


__version__ = "0.2.0"

pkg_dir = Path(__file__).resolve().parent.parent.parent
schema_dir = pkg_dir / "schema"
primer_scheme_schema_path = schema_dir / "primer-scheme_schema.yml"
manifest_schema_path = schema_dir / "manifest_schema.json"
organisms_path = schema_dir / "organisms.yml"
# print(f"{pkg_dir=} {schema_dir=} {schema_path=}")

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
            "level": "INFO",  # Set level to INFO
            "stream": "ext://sys.stderr",
        },
    },
    "root": {
        "handlers": ["console"],
        "level": "INFO",  # Set level to INFO
    },
    "loggers": {
        "primaschema": {
            "handlers": ["console"],
            "level": "INFO",  # Set level to INFO
            "propagate": False,
        },
    },
}

logging.config.dictConfig(logging_config)
