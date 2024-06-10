# Auto generated from primer_scheme.yml by pythongen.py version: 0.0.1
# Generation date: 2024-06-08T10:58:04
# Schema: primer-scheme
#
# id: https://github.com/pha4ge/primer-schemes/schemas/primer-scheme
# description: Data model for tiling PCR primer scheme definitions
# license: https://creativecommons.org/publicdomain/zero/1.0/

import dataclasses
import re
from jsonasobj2 import as_dict
from typing import Optional, List, Union, Dict, ClassVar, Any
from dataclasses import dataclass
from linkml_runtime.linkml_model.meta import EnumDefinition, PermissibleValue

from linkml_runtime.utils.slot import Slot
from linkml_runtime.utils.metamodelcore import empty_list
from linkml_runtime.utils.yamlutils import YAMLRoot, extended_str
from linkml_runtime.utils.dataclass_extensions_376 import (
    dataclasses_init_fn_with_kwargs,
)
from linkml_runtime.utils.enumerations import EnumDefinitionImpl
from rdflib import URIRef
from linkml_runtime.utils.curienamespace import CurieNamespace

metamodel_version = "1.7.0"
version = "0.9.0"

# Overwrite dataclasses _init_fn to add **kwargs in __init__
dataclasses._init_fn = dataclasses_init_fn_with_kwargs

# Namespaces
GENEPIO = CurieNamespace("GENEPIO", "http://purl.obolibrary.org/obo/GENEPIO_")
IAO = CurieNamespace("IAO", "https://bioregistry.io/reference/iao:")
ORCID = CurieNamespace("ORCID", "http://identifiers.org/orcid/")
LINKML = CurieNamespace("linkml", "https://w3id.org/linkml/")
SCHEMA = CurieNamespace("schema", "http://schema.org/")
DEFAULT_ = CurieNamespace(
    "", "https://github.com/pha4ge/primer-schemes/schemas/primer-scheme/"
)


# Types


# Class references
class PrimerSchemeName(extended_str):
    pass


@dataclass
class PrimerScheme(YAMLRoot):
    """
    A description of an amplicon primer scheme
    """

    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = URIRef(
        "https://github.com/pha4ge/primer-schemes/schemas/primer-scheme/PrimerScheme"
    )
    class_class_curie: ClassVar[str] = None
    class_name: ClassVar[str] = "PrimerScheme"
    class_model_uri: ClassVar[URIRef] = URIRef(
        "https://github.com/pha4ge/primer-schemes/schemas/primer-scheme/PrimerScheme"
    )

    name: Union[str, PrimerSchemeName] = None
    schema_version: str = None
    organism: str = None
    developers: Union[Union[dict, "Entity"], List[Union[dict, "Entity"]]] = None
    repository_url: str = None
    organism_aliases: Optional[Union[str, List[str]]] = empty_list()
    display_name: Optional[str] = None
    primer_scheme_status: Optional[Union[str, "PrimerSchemeStatus"]] = "PUBLISHED"
    aliases: Optional[Union[str, List[str]]] = empty_list()
    derived_from: Optional[str] = None
    vendors: Optional[Union[Union[dict, "Vendor"], List[Union[dict, "Vendor"]]]] = (
        empty_list()
    )
    amplicon_size: Optional[int] = None
    notes: Optional[Union[str, List[str]]] = empty_list()
    citations: Optional[Union[str, List[str]]] = empty_list()
    primer_checksum: Optional[str] = None
    reference_checksum: Optional[str] = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self._is_empty(self.name):
            self.MissingRequiredField("name")
        if not isinstance(self.name, PrimerSchemeName):
            self.name = PrimerSchemeName(self.name)

        if self._is_empty(self.schema_version):
            self.MissingRequiredField("schema_version")
        if not isinstance(self.schema_version, str):
            self.schema_version = str(self.schema_version)

        if self._is_empty(self.organism):
            self.MissingRequiredField("organism")
        if not isinstance(self.organism, str):
            self.organism = str(self.organism)

        if self._is_empty(self.developers):
            self.MissingRequiredField("developers")
        if not isinstance(self.developers, list):
            self.developers = [self.developers] if self.developers is not None else []
        self.developers = [
            v if isinstance(v, Entity) else Entity(**as_dict(v))
            for v in self.developers
        ]

        if self._is_empty(self.repository_url):
            self.MissingRequiredField("repository_url")
        if not isinstance(self.repository_url, str):
            self.repository_url = str(self.repository_url)

        if not isinstance(self.organism_aliases, list):
            self.organism_aliases = (
                [self.organism_aliases] if self.organism_aliases is not None else []
            )
        self.organism_aliases = [
            v if isinstance(v, str) else str(v) for v in self.organism_aliases
        ]

        if self.display_name is not None and not isinstance(self.display_name, str):
            self.display_name = str(self.display_name)

        if self.primer_scheme_status is not None and not isinstance(
            self.primer_scheme_status, PrimerSchemeStatus
        ):
            self.primer_scheme_status = PrimerSchemeStatus(self.primer_scheme_status)

        if not isinstance(self.aliases, list):
            self.aliases = [self.aliases] if self.aliases is not None else []
        self.aliases = [v if isinstance(v, str) else str(v) for v in self.aliases]

        if self.derived_from is not None and not isinstance(self.derived_from, str):
            self.derived_from = str(self.derived_from)

        if not isinstance(self.vendors, list):
            self.vendors = [self.vendors] if self.vendors is not None else []
        self.vendors = [
            v if isinstance(v, Vendor) else Vendor(**as_dict(v)) for v in self.vendors
        ]

        if self.amplicon_size is not None and not isinstance(self.amplicon_size, int):
            self.amplicon_size = int(self.amplicon_size)

        if not isinstance(self.notes, list):
            self.notes = [self.notes] if self.notes is not None else []
        self.notes = [v if isinstance(v, str) else str(v) for v in self.notes]

        if not isinstance(self.citations, list):
            self.citations = [self.citations] if self.citations is not None else []
        self.citations = [v if isinstance(v, str) else str(v) for v in self.citations]

        if self.primer_checksum is not None and not isinstance(
            self.primer_checksum, str
        ):
            self.primer_checksum = str(self.primer_checksum)

        if self.reference_checksum is not None and not isinstance(
            self.reference_checksum, str
        ):
            self.reference_checksum = str(self.reference_checksum)

        super().__post_init__(**kwargs)


@dataclass
class Vendor(YAMLRoot):
    """
    Vendor of the primers described in the amplicon scheme or a kit containing these primers
    """

    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = GENEPIO["0100674"]
    class_class_curie: ClassVar[str] = "GENEPIO:0100674"
    class_name: ClassVar[str] = "Vendor"
    class_model_uri: ClassVar[URIRef] = URIRef(
        "https://github.com/pha4ge/primer-schemes/schemas/primer-scheme/Vendor"
    )

    organisation_name: str = None
    kit_name: Optional[str] = None
    home_page: Optional[str] = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self._is_empty(self.organisation_name):
            self.MissingRequiredField("organisation_name")
        if not isinstance(self.organisation_name, str):
            self.organisation_name = str(self.organisation_name)

        if self.kit_name is not None and not isinstance(self.kit_name, str):
            self.kit_name = str(self.kit_name)

        if self.home_page is not None and not isinstance(self.home_page, str):
            self.home_page = str(self.home_page)

        super().__post_init__(**kwargs)


class Entity(YAMLRoot):
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = URIRef(
        "https://github.com/pha4ge/primer-schemes/schemas/primer-scheme/Entity"
    )
    class_class_curie: ClassVar[str] = None
    class_name: ClassVar[str] = "Entity"
    class_model_uri: ClassVar[URIRef] = URIRef(
        "https://github.com/pha4ge/primer-schemes/schemas/primer-scheme/Entity"
    )


@dataclass
class Person(Entity):
    """
    A natural person
    """

    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = URIRef(
        "https://github.com/pha4ge/primer-schemes/schemas/primer-scheme/Person"
    )
    class_class_curie: ClassVar[str] = None
    class_name: ClassVar[str] = "Person"
    class_model_uri: ClassVar[URIRef] = URIRef(
        "https://github.com/pha4ge/primer-schemes/schemas/primer-scheme/Person"
    )

    person_name: str = None
    orcid: Optional[str] = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self._is_empty(self.person_name):
            self.MissingRequiredField("person_name")
        if not isinstance(self.person_name, str):
            self.person_name = str(self.person_name)

        if self.orcid is not None and not isinstance(self.orcid, str):
            self.orcid = str(self.orcid)

        super().__post_init__(**kwargs)


@dataclass
class Organisation(Entity):
    """
    An organisation
    """

    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = URIRef(
        "https://github.com/pha4ge/primer-schemes/schemas/primer-scheme/Organisation"
    )
    class_class_curie: ClassVar[str] = None
    class_name: ClassVar[str] = "Organisation"
    class_model_uri: ClassVar[URIRef] = URIRef(
        "https://github.com/pha4ge/primer-schemes/schemas/primer-scheme/Organisation"
    )

    organisation_name: str = None
    home_page: Optional[str] = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self._is_empty(self.organisation_name):
            self.MissingRequiredField("organisation_name")
        if not isinstance(self.organisation_name, str):
            self.organisation_name = str(self.organisation_name)

        if self.home_page is not None and not isinstance(self.home_page, str):
            self.home_page = str(self.home_page)

        super().__post_init__(**kwargs)


# Enumerations
class PrimerSchemeStatus(EnumDefinitionImpl):
    """
    Status of this amplicon primer scheme
    """

    PUBLISHED = PermissibleValue(text="PUBLISHED")
    DEPRECATED = PermissibleValue(text="DEPRECATED")
    DRAFT = PermissibleValue(text="DRAFT")

    _defn = EnumDefinition(
        name="PrimerSchemeStatus",
        description="Status of this amplicon primer scheme",
    )


# Slots
class slots:
    pass


slots.schema_version = Slot(
    uri=DEFAULT_.schema_version,
    name="schema_version",
    curie=DEFAULT_.curie("schema_version"),
    model_uri=DEFAULT_.schema_version,
    domain=None,
    range=str,
)

slots.name = Slot(
    uri=GENEPIO["0001456"],
    name="name",
    curie=GENEPIO.curie("0001456"),
    model_uri=DEFAULT_.name,
    domain=None,
    range=URIRef,
    pattern=re.compile(r"^[\da-z0-9_.-]+$"),
)

slots.display_name = Slot(
    uri=DEFAULT_.display_name,
    name="display_name",
    curie=DEFAULT_.curie("display_name"),
    model_uri=DEFAULT_.display_name,
    domain=None,
    range=Optional[str],
)

slots.organism = Slot(
    uri=GENEPIO["0100682"],
    name="organism",
    curie=GENEPIO.curie("0100682"),
    model_uri=DEFAULT_.organism,
    domain=None,
    range=str,
)

slots.organism_aliases = Slot(
    uri=DEFAULT_.organism_aliases,
    name="organism_aliases",
    curie=DEFAULT_.curie("organism_aliases"),
    model_uri=DEFAULT_.organism_aliases,
    domain=None,
    range=Optional[Union[str, List[str]]],
)

slots.aliases = Slot(
    uri=GENEPIO["0100670"],
    name="aliases",
    curie=GENEPIO.curie("0100670"),
    model_uri=DEFAULT_.aliases,
    domain=None,
    range=Optional[Union[str, List[str]]],
)

slots.developers = Slot(
    uri=GENEPIO["0100673"],
    name="developers",
    curie=GENEPIO.curie("0100673"),
    model_uri=DEFAULT_.developers,
    domain=None,
    range=Union[Union[dict, Entity], List[Union[dict, Entity]]],
)

slots.vendors = Slot(
    uri=DEFAULT_.vendors,
    name="vendors",
    curie=DEFAULT_.curie("vendors"),
    model_uri=DEFAULT_.vendors,
    domain=None,
    range=Optional[Union[Union[dict, Vendor], List[Union[dict, Vendor]]]],
)

slots.amplicon_size = Slot(
    uri=GENEPIO["0001449"],
    name="amplicon_size",
    curie=GENEPIO.curie("0001449"),
    model_uri=DEFAULT_.amplicon_size,
    domain=None,
    range=Optional[int],
)

slots.repository_url = Slot(
    uri=GENEPIO["0100683"],
    name="repository_url",
    curie=GENEPIO.curie("0100683"),
    model_uri=DEFAULT_.repository_url,
    domain=None,
    range=str,
)

slots.notes = Slot(
    uri=GENEPIO["0100672"],
    name="notes",
    curie=GENEPIO.curie("0100672"),
    model_uri=DEFAULT_.notes,
    domain=None,
    range=Optional[Union[str, List[str]]],
)

slots.primer_scheme_status = Slot(
    uri=GENEPIO["0100681"],
    name="primer_scheme_status",
    curie=GENEPIO.curie("0100681"),
    model_uri=DEFAULT_.primer_scheme_status,
    domain=None,
    range=Optional[Union[str, "PrimerSchemeStatus"]],
)

slots.citations = Slot(
    uri=IAO["0000301"],
    name="citations",
    curie=IAO.curie("0000301"),
    model_uri=DEFAULT_.citations,
    domain=None,
    range=Optional[Union[str, List[str]]],
)

slots.primer_checksum = Slot(
    uri=GENEPIO["0100675"],
    name="primer_checksum",
    curie=GENEPIO.curie("0100675"),
    model_uri=DEFAULT_.primer_checksum,
    domain=None,
    range=Optional[str],
)

slots.reference_checksum = Slot(
    uri=DEFAULT_.reference_checksum,
    name="reference_checksum",
    curie=DEFAULT_.curie("reference_checksum"),
    model_uri=DEFAULT_.reference_checksum,
    domain=None,
    range=Optional[str],
)

slots.derived_from = Slot(
    uri=GENEPIO["0100671"],
    name="derived_from",
    curie=GENEPIO.curie("0100671"),
    model_uri=DEFAULT_.derived_from,
    domain=None,
    range=Optional[str],
)

slots.person_name = Slot(
    uri=SCHEMA.name,
    name="person_name",
    curie=SCHEMA.curie("name"),
    model_uri=DEFAULT_.person_name,
    domain=None,
    range=str,
)

slots.orcid = Slot(
    uri=DEFAULT_.orcid,
    name="orcid",
    curie=DEFAULT_.curie("orcid"),
    model_uri=DEFAULT_.orcid,
    domain=None,
    range=Optional[str],
)

slots.organisation_name = Slot(
    uri=DEFAULT_.organisation_name,
    name="organisation_name",
    curie=DEFAULT_.curie("organisation_name"),
    model_uri=DEFAULT_.organisation_name,
    domain=None,
    range=str,
)

slots.home_page = Slot(
    uri=DEFAULT_.home_page,
    name="home_page",
    curie=DEFAULT_.curie("home_page"),
    model_uri=DEFAULT_.home_page,
    domain=None,
    range=Optional[str],
)

slots.vendor__kit_name = Slot(
    uri=GENEPIO["0100693"],
    name="vendor__kit_name",
    curie=GENEPIO.curie("0100693"),
    model_uri=DEFAULT_.vendor__kit_name,
    domain=None,
    range=Optional[str],
)
