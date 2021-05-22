import argparse
import csv
import sqlite3
from collections import defaultdict
from dataclasses import dataclass
from itertools import combinations
from pathlib import Path
from typing import TextIO

import numpy as np


def diff_l2(w1: dict[str, float], w2: dict[str, float]) -> float:
    keys = set(w1.keys()) | set(w2.keys())
    w = {k: w1.get(k, 0.0) - w2.get(k, 0.0) for k in keys}
    return np.linalg.norm(list(w.values()), ord=2)


@dataclass
class PretrainedModel:
    method: str
    tissue: str
    gene_weights: dict[str, dict[str, float]]

    @classmethod
    def from_db(cls, db_path: Path):
        gene_weights = defaultdict(dict)
        with sqlite3.connect(db_path) as con:
            cursor = con.cursor()
            for rsid, gene, weight in cursor.execute(
                "SELECT rsid, gene, weight FROM weights"
            ):
                gene_weights[gene][rsid] = weight
        method, tissue = db_path.stem.split("_", maxsplit=1)
        return cls(method, tissue, gene_weights)

    @property
    def genes(self):
        return set(self.gene_weights.keys())


@dataclass
class WeightsSummaryWriter:
    file: TextIO
    delimiter: str = ","

    def __post_init__(self):
        self._writer = csv.writer(self.file, delimiter=self.delimiter)
        self._writer.writerow(("method1", "method2", "tissue", "gene", "diff_l2"))

    def compare_and_write(
        self, m1: PretrainedModel, m2: PretrainedModel, geneid: str
    ) -> None:
        assert m1.tissue == m2.tissue
        w1 = m1.gene_weights.get(geneid, {})
        w2 = m2.gene_weights.get(geneid, {})
        if len(w1) > 0 or len(w2) > 0:
            row = (m1.method, m2.method, m1.tissue, geneid, diff_l2(w1, w2))
            self._writer.writerow(row)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path)
    parser.add_argument("dbs", nargs="*", type=Path)
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()

    # Load models and use union of genes as geneset
    models = [PretrainedModel.from_db(db) for db in args.dbs]
    geneset = set.union(*tuple(m.genes for m in models))

    with open(args.out, mode="w") as f:
        writer = WeightsSummaryWriter(f)
        for geneid in geneset:
            for m1, m2 in combinations(models, 2):
                writer.compare_and_write(m1, m2, geneid)
