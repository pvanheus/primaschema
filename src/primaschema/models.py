from pydantic import BaseModel, ConfigDict, computed_field, model_validator
from typing import Literal
from typing_extensions import Self


from . import logger


class PrimerModel(BaseModel):
    """A primer as represented by a BED record in Primal Scheme v3 format"""

    model_config = ConfigDict(str_strip_whitespace=True)
    chrom: str
    chrom_start: int
    chrom_end: int
    name: str
    pool_name: int
    strand: Literal["+", "-"]
    sequence: str

    @computed_field
    @property
    def name_parts(self) -> list[str]:
        parts = self.name.split("_")
        if len(parts) != 4:
            raise ValueError(
                "Name must be in the format '{name}_{amplicon-number}_{LEFT|RIGHT}_{primer-number}'"
            )
        return parts

    @computed_field
    @property
    def number(self) -> int:
        return int(self.name.split("_")[3])


class AmpliconModel(BaseModel):
    """An amplicon as represented by two or more primer records"""

    primers: list[PrimerModel]

    @computed_field
    @property
    def number(self) -> int:
        return int(self.primers[0].name_parts[1])

    @computed_field
    @property
    def min_start(self) -> int:
        return min(p.chrom_start for p in self.primers)

    @computed_field
    @property
    def max_end(self) -> int:
        return max(p.chrom_end for p in self.primers)


class BedModel(BaseModel):
    """A BED file as represented by a collection of amplicons each comprising primer records"""

    amplicons: dict[str, list[AmpliconModel]]

    @model_validator(mode="after")
    def check_tiling_bed(self) -> Self:
        chroms_amplicon_boundaries = {}
        for chrom, amplicons in self.amplicons.items():
            chrom = chrom.partition(" ")[0]
            chroms_amplicon_boundaries[chrom] = [
                (a.min_start, a.max_end) for a in amplicons
            ]

        for chrom, amplicon_boundaries in chroms_amplicon_boundaries.items():
            BedModel.check_tiling(amplicon_boundaries)
        return self

    @staticmethod
    def check_overlap(interval1: tuple[int, int], interval2: tuple[int, int]) -> bool:
        """Check if two intervals overlap"""
        return max(interval1[0], interval2[0]) <= min(interval1[1], interval2[1])

    @staticmethod
    def check_tiling(intervals):
        """Verify if each interval overlaps exactly once with the interval before it and after it,
        and does not overlap with any other interval"""
        logger.debug(f"{intervals=}")
        n = len(intervals)
        if n < 2:
            raise ValueError("Fewer than two amplicons detected")

        # Check first interval
        if not BedModel.check_overlap(intervals[0], intervals[1]):
            raise ValueError("First and second amplicons do not overlap")
        for j in range(2, n):
            if BedModel.check_overlap(intervals[0], intervals[j]):
                raise ValueError(
                    f"First amplicon overlaps with more than one amplicon ({j})"
                )

        # Check last interval
        if not BedModel.check_overlap(intervals[-1], intervals[-2]):
            raise ValueError("Penultimate and last amplicons do not overlap")
        for j in range(n - 2):
            if BedModel.check_overlap(intervals[-1], intervals[j]):
                raise ValueError(
                    f"Last amplicon overlaps with more than one amplicon ({j})"
                )

        # Check other intervals
        for i in range(1, n - 1):
            if not (
                BedModel.check_overlap(intervals[i], intervals[i - 1])
                and BedModel.check_overlap(intervals[i], intervals[i + 1])
            ):
                raise ValueError(f"Amplicons {i-1} and {i} do not overlap")

            # Ensure interval overlaps only with previous and next
            for j in range(n):
                if j != i - 1 and j != i + 1 and j != i:
                    if BedModel.check_overlap(intervals[i], intervals[j]):
                        raise ValueError(
                            f"Amplicon {i} overlaps with more than two other amplicons ({j})"
                        )
