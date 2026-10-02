"""Chunk crossover: does splicing a run of molvector blocks transfer chemistry?

For each pair of MOSES molecules, one contiguous chunk of parent A replaces
one contiguous chunk of parent B, and the child is scored on three things:

  valid      sanitizes to a single-fragment molecule
  changed    differs from both parents
  transfer   carries substructures unique to A *and* unique to B, so the
             chunk really moved rather than the child collapsing back to
             one parent

Methods:
  mv_splice         runs of molvector atom blocks, decoded strictly
  mv_splice_repair  the same splice through the valence-repairing decoder
  smiles_splice     runs of SMILES tokens (string crossover, GA-of-2019 style)
  selfies_splice    runs of SELFIES symbols (STONED style)
  brics             BRICS fragment recombination, the chemistry-aware baseline

  python crossover_chunks.py --data-dir data
"""
import argparse
import json
import random
import statistics

from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import BRICS, rdFingerprintGenerator

import selfies as sf

import bench  # noqa: F401  (also puts the repo root on sys.path)
import molvector as mv
from representations import BLOCK, decode_with_repair, molvector_of, tokenize

RDLogger.DisableLog("rdApp.*")
FPGEN = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
METHODS = ("mv_splice", "mv_splice_repair", "smiles_splice", "selfies_splice", "brics")


def chunk_bounds(n, rng, max_frac=0.5):
    """A contiguous [start, end) of length 1..max_frac*n."""
    size = max(1, min(n, int(rng.random() * max_frac * n) + 1))
    start = rng.randrange(n - size + 1)
    return start, start + size


def splice(a, b, rng):
    """Replace a chunk of b with a chunk of a."""
    i, j = chunk_bounds(len(a), rng)
    k, l = chunk_bounds(len(b), rng)
    return b[:k] + a[i:j] + b[l:]


def mv_child(ma, mb, rng, repair):
    va, vb = molvector_of(ma, False), molvector_of(mb, False)
    blocks_a = [va[i:i + BLOCK] for i in range(0, len(va), BLOCK)]
    blocks_b = [vb[i:i + BLOCK] for i in range(0, len(vb), BLOCK)]
    v = [x for block in splice(blocks_a, blocks_b, rng) for x in block]
    if repair:
        return decode_with_repair(v)
    try:
        mol = mv.decode(v)
        Chem.SanitizeMol(mol)
        return Chem.MolToSmiles(mol)
    except Exception:
        return None


def token_child(ma, mb, rng, rep):
    ta, tb = tokenize(ma, rep), tokenize(mb, rep)
    toks = splice(ta, tb, rng)
    if rep == "selfies":
        try:
            smi = sf.decoder("".join(toks))
        except Exception:
            return None
    else:
        smi = "".join(toks)
    mol = Chem.MolFromSmiles(smi)
    return Chem.MolToSmiles(mol) if mol and mol.GetNumAtoms() else None


def brics_child(ma, mb, rng):
    """One BRICS recombination of fragments drawn from both parents."""
    frags = list(BRICS.BRICSDecompose(ma)) + list(BRICS.BRICSDecompose(mb))
    mols = [Chem.MolFromSmiles(f) for f in frags]
    mols = [m for m in mols if m is not None]
    if len(mols) < 2:
        return None
    rng.shuffle(mols)
    try:
        for built in BRICS.BRICSBuild(mols, scrambleReagents=True, maxDepth=1):
            built.UpdatePropertyCache(strict=False)
            Chem.SanitizeMol(built)
            return Chem.MolToSmiles(built)
    except Exception:
        return None
    return None


def bits(smi):
    mol = Chem.MolFromSmiles(smi)
    return set(FPGEN.GetSparseFingerprint(mol).GetOnBits()) if mol else set()


def similarity(x, y):
    return DataStructs.TanimotoSimilarity(FPGEN.GetFingerprint(Chem.MolFromSmiles(x)),
                                          FPGEN.GetFingerprint(Chem.MolFromSmiles(y)))


def summarize(method, records, min_unique=3):
    n = len(records)
    valid = [(a, b, c) for a, b, c in records if c and "." not in c]
    changed = [(a, b, c) for a, b, c in valid if c != a and c != b]
    transferred, sims = [], []
    for a, b, c in changed:
        ba, bb, bc = bits(a), bits(b), bits(c)
        if len(bc & (ba - bb)) >= min_unique and len(bc & (bb - ba)) >= min_unique:
            transferred.append((a, b, c))
        sims.append(max(similarity(c, a), similarity(c, b)))
    sanitizable = [r for r in records if r[2]]
    return dict(method=method, n=n,
                sanitizable=len(sanitizable) / n,
                valid=len(valid) / n,
                changed=len(changed) / n,
                transfer=len(transferred) / n,
                transfer_of_changed=len(transferred) / max(1, len(changed)),
                median_max_sim=statistics.median(sims) if sims else float("nan"))


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--data-dir", default="data")
    p.add_argument("--work-dir", default="runs")
    p.add_argument("--pairs", type=int, default=2000)
    p.add_argument("--seed", type=int, default=0)
    a = p.parse_args()

    smiles = random.Random(1).sample(bench.moses_smiles(a.data_dir, "test"), 2 * a.pairs)
    mols = [Chem.MolFromSmiles(s) for s in smiles]
    pairs = [(mols[2 * i], mols[2 * i + 1]) for i in range(a.pairs)]

    rows = []
    for method in METHODS:
        rng = random.Random(a.seed)
        records = []
        for ma, mb in pairs:
            pa, pb = Chem.MolToSmiles(ma), Chem.MolToSmiles(mb)
            if method == "mv_splice":
                child = mv_child(ma, mb, rng, repair=False)
            elif method == "mv_splice_repair":
                child = mv_child(ma, mb, rng, repair=True)
            elif method == "smiles_splice":
                child = token_child(ma, mb, rng, "smiles")
            elif method == "selfies_splice":
                child = token_child(ma, mb, rng, "selfies")
            else:
                child = brics_child(ma, mb, rng)
            records.append((pa, pb, child))
        rows.append(summarize(method, records))
        print(json.dumps(rows[-1]))

    cols = ["method", "sanitizable", "valid", "changed", "transfer", "transfer_of_changed", "median_max_sim"]
    lines = ["| " + " | ".join(cols) + " |", "|" + "---|" * len(cols)]
    for r in rows:
        lines.append("| " + " | ".join(r[c] if isinstance(r[c], str) else "%.3f" % r[c] for c in cols) + " |")
    text = "\n".join(lines)
    with open(a.work_dir + "/crossover_chunks.md", "w") as fh:
        fh.write(text + "\n")
    print(text)


if __name__ == "__main__":
    main()
