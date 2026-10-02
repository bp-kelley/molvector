"""Single-edit mutation test: validity and locality per representation.

For held-out MOSES molecules, apply one random edit to the token sequence
(replace, insert or delete a token drawn uniformly from the training
vocabulary) and record whether the result decodes to a valid single-fragment
molecule, whether it differs from the parent, and the Morgan Tanimoto
similarity between parent and child.

`mv_int` edits one integer of the raw molvector instead of a token: an atom
field, a bond type or an offset, replaced with a value seen in that field.
Inserting or deleting a whole atom block are its insert and delete edits.

  python mutation_locality.py --work-dir runs --data-dir data
"""
import argparse
import collections
import json
import os
import random
import statistics

from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import rdFingerprintGenerator

import bench
import molvector as mv
from representations import BLOCK, REPRESENTATIONS, detokenize, molvector_of, tokenize

RDLogger.DisableLog("rdApp.*")
FPGEN = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)


def fingerprint(smi):
    return FPGEN.GetFingerprint(Chem.MolFromSmiles(smi))


def mutate_tokens(tokens, alphabet, rng):
    t = list(tokens)
    op = rng.choice(("replace", "insert", "delete"))
    i = rng.randrange(len(t) + (op == "insert"))
    if op == "replace":
        t[i] = rng.choice(alphabet)
    elif op == "insert":
        t.insert(i, rng.choice(alphabet))
    elif len(t) > 1:
        del t[i]
    return t


def mv_int_values(mols):
    """Values observed at each position of a molvector block."""
    seen = [collections.Counter() for _ in range(BLOCK)]
    for m in mols:
        v = molvector_of(m, False)
        for k, x in enumerate(v):
            seen[k % BLOCK][x] += 1
    return [sorted(c) for c in seen]


def mutate_mv_int(v, values, rng):
    v = list(v)
    n = len(v) // BLOCK
    op = rng.choice(("replace", "insert", "delete"))
    if op == "replace":
        k = rng.randrange(len(v))
        v[k] = rng.choice(values[k % BLOCK])
    elif op == "insert":
        i = rng.randrange(n + 1) * BLOCK
        block = [rng.choice(values[k]) for k in range(BLOCK)]
        v[i:i] = block
    elif n > 1:
        i = rng.randrange(n) * BLOCK
        del v[i:i + BLOCK]
    return v


def decode_mv_int(v):
    try:
        mol = mv.decode(v)
        Chem.SanitizeMol(mol)
        mol = Chem.MolFromSmiles(Chem.MolToSmiles(mol))
        return Chem.MolToSmiles(mol) if mol and mol.GetNumAtoms() else None
    except Exception:
        return None


def summarize(rep, parents, children):
    n = len(children)
    valid = [(p, c) for p, c in zip(parents, children) if c and "." not in c]
    changed = [(p, c) for p, c in valid if p != c]
    sims = [DataStructs.TanimotoSimilarity(fingerprint(p), fingerprint(c)) for p, c in changed]
    return dict(rep=rep, n=n,
                valid=len(valid) / n,
                valid_and_changed=len(changed) / n,
                median_sim=statistics.median(sims) if sims else float("nan"),
                frac_sim_ge_0_6=sum(s >= 0.6 for s in sims) / max(1, len(sims)),
                useful=sum(s >= 0.6 for s in sims) / n)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--work-dir", default="runs")
    p.add_argument("--data-dir", default="data")
    p.add_argument("--n-mols", type=int, default=2000)
    p.add_argument("--edits", type=int, default=5, help="independent single edits per molecule")
    p.add_argument("--seed", type=int, default=0)
    a = p.parse_args()

    rng = random.Random(a.seed)
    smiles = random.Random(1).sample(bench.moses_smiles(a.data_dir, "test"), a.n_mols)
    mols = [Chem.MolFromSmiles(s) for s in smiles]
    canon = [Chem.MolToSmiles(m) for m in mols]
    rows = []
    for rep in REPRESENTATIONS:
        meta_path = os.path.join(a.work_dir, rep, "meta.json")
        if not os.path.exists(meta_path):
            print("%s: not prepared, skipping" % rep)
            continue
        alphabet = [t for t in json.load(open(meta_path))["vocab"] if not t.startswith("<")]
        parents, children = [], []
        for m, c in zip(mols, canon):
            toks = tokenize(m, rep)
            for _ in range(a.edits):
                parents.append(c)
                children.append(detokenize(mutate_tokens(toks, alphabet, rng), rep))
        rows.append(summarize(rep, parents, children))
        print(json.dumps(rows[-1]))
    values = mv_int_values(mols)
    parents, children = [], []
    for m, c in zip(mols, canon):
        v = molvector_of(m, False)
        for _ in range(a.edits):
            parents.append(c)
            children.append(decode_mv_int(mutate_mv_int(v, values, rng)))
    rows.append(summarize("mv_int", parents, children))
    print(json.dumps(rows[-1]))

    cols = ["rep", "valid", "valid_and_changed", "median_sim", "frac_sim_ge_0_6", "useful"]
    lines = ["| " + " | ".join(cols) + " |", "|" + "---|" * len(cols)]
    for r in rows:
        lines.append("| " + " | ".join(r[c] if isinstance(r[c], str) else "%.3f" % r[c] for c in cols) + " |")
    text = "\n".join(lines)
    with open(os.path.join(a.work_dir, "mutation_locality.md"), "w") as fh:
        fh.write(text + "\n")
    print(text)


if __name__ == "__main__":
    main()
