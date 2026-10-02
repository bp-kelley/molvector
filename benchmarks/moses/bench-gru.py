"""The MOSES benchmark with a GRU decoder in place of the transformer.

Tests whether molvector's advantage in 2019 came from the decoders of the
time: a recurrent model has to carry ring-closure digits and open branches
across the whole SMILES string, while molvector bonds are local offsets.

Same commands and flags as bench.py. Differences:
  --work-dir         defaults to runs_gru, so transformer runs are untouched
  --prepared-dir     tokenized data is linked from here (default runs) when
                     runs_gru/<rep> has none, so there is no need to prepare
                     again
  --d-model/--layers default to 768 and 3, about the size of the 6x384
                     transformer (~10.6M parameters)

  python bench-gru.py train  --rep all
  python bench-gru.py sample --rep all --fcd-device cpu
  python bench-gru.py report
"""
import argparse
import os
import sys

import torch.nn as nn

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import bench  # noqa: E402
from representations import REPRESENTATIONS  # noqa: E402

PREPARED_FILES = ("meta.json", "train.npz", "heldout.npz")


class GRULM(nn.Module):
    """Drop-in replacement for bench.GPT; `heads` is ignored."""

    def __init__(self, vocab, max_len, d=768, layers=3, heads=None, dropout=0.1):
        super().__init__()
        self.max_len = max_len
        self.tok = nn.Embedding(vocab, d)
        self.gru = nn.GRU(d, d, num_layers=layers, dropout=dropout if layers > 1 else 0.0, batch_first=True)
        self.drop = nn.Dropout(dropout)
        self.head = nn.Linear(d, vocab, bias=False)
        self.head.weight = self.tok.weight
        nn.init.normal_(self.tok.weight, std=0.02)

    def forward(self, x):
        h, _ = self.gru(self.tok(x))
        return self.head(self.drop(h))


def link_prepared(prepared_dir, work_dir):
    for rep in REPRESENTATIONS:
        src, dst = os.path.join(prepared_dir, rep), os.path.join(work_dir, rep)
        if os.path.exists(os.path.join(dst, "meta.json")) or not os.path.exists(os.path.join(src, "meta.json")):
            continue
        os.makedirs(dst, exist_ok=True)
        for f in PREPARED_FILES:
            os.symlink(os.path.abspath(os.path.join(src, f)), os.path.join(dst, f))


def main():
    pre = argparse.ArgumentParser(add_help=False)
    pre.add_argument("--prepared-dir", default="runs")
    pre.add_argument("--work-dir", default="runs_gru")
    known, rest = pre.parse_known_args()
    if rest and rest[0] != "prepare":
        link_prepared(known.prepared_dir, known.work_dir)
    # defaults first so flags given on the command line override them
    sys.argv = [sys.argv[0]] + rest[:1] + ["--work-dir", known.work_dir, "--d-model", "768", "--layers", "3"] + rest[1:]
    bench.GPT = GRULM
    bench.main()


if __name__ == "__main__":
    main()
