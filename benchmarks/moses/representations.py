"""Tokenizers for the MOSES benchmark.

Each representation turns an RDKit molecule into a list of string tokens
and turns a list of tokens back into a SMILES string (or None when the
tokens do not decode to a sanitizable molecule).

  smiles    regex-tokenized SMILES (Schwaller et al. tokenizer)
  selfies   SELFIES symbols
  mv_atom   molvector, one composite token per atom block
  mv_field  molvector, one token for the atom fields plus one token
            per non-empty (bond type, offset) record
  mv_back   like mv_field, but each bond is written once, from its later
            atom (negative offsets only), so the two ends cannot disagree
  mv_lean   like mv_back, on the Kekule form with no hydrogen count, so the
            model never has to get aromaticity or H counts consistent
"""
import os
import re
import sys

import selfies as sf
from rdkit import Chem

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
import molvector as mv  # noqa: E402

REPRESENTATIONS = ("smiles", "selfies", "mv_atom", "mv_field", "mv_back", "mv_lean")

SMILES_REGEX = re.compile(
    r"(\[[^\]]+]|Br?|Cl?|N|O|S|P|F|I|b|c|n|o|s|p|\(|\)|\.|=|#|-|\+|\\|\/|:|~|@|\?|>|\*|\$|%[0-9]{2}|[0-9])")

BLOCK = mv.atom_size + mv.bond_chunk_size


def atom_order(mol, randomize):
    """Atom output order of a canonical or random SMILES traversal."""
    Chem.MolToSmiles(mol, canonical=not randomize, doRandom=randomize)
    return tuple(int(i) for i in mol.GetProp("_smilesAtomOutputOrder").strip("[]").split(",") if i.strip())


def molvector_of(mol, randomize):
    order = atom_order(mol, randomize)
    return mv.encode(mol, lambda m: iter([order]))[0]


def tokenize(mol, rep, randomize=False):
    if rep == "smiles":
        smi = Chem.MolToSmiles(mol, canonical=not randomize, doRandom=randomize)
        return SMILES_REGEX.findall(smi)
    if rep == "selfies":
        smi = Chem.MolToSmiles(mol, canonical=not randomize, doRandom=randomize)
        return list(sf.split_selfies(sf.encoder(smi)))
    if rep == "mv_lean":
        mol = Chem.Mol(mol)
        Chem.Kekulize(mol, clearAromaticFlags=True)
    v = molvector_of(mol, randomize)
    blocks = [v[i:i + BLOCK] for i in range(0, len(v), BLOCK)]
    if rep == "mv_atom":
        return [",".join(map(str, b)) for b in blocks]
    if rep in ("mv_field", "mv_back", "mv_lean"):
        toks = []
        for b in blocks:
            if rep == "mv_lean":
                toks.append("A%d,%d" % tuple(b[:2]))
            else:
                toks.append("A%d,%d,%d" % tuple(b[:mv.atom_size]))
            for k in range(mv.max_num_bonds):
                btype, off = b[mv.atom_size + 2 * k], b[mv.atom_size + 2 * k + 1]
                if off and (rep == "mv_field" or off < 0):
                    toks.append("B%d,%d" % (btype, off))
        return toks
    raise ValueError(rep)


def _blocks_from_fields(tokens):
    blocks = []
    for t in tokens:
        if t.startswith("A"):
            blocks.append([int(x) for x in t[1:].split(",")])
        elif t.startswith("B") and blocks:
            # bond records beyond max_num_bonds are dropped
            if len(blocks[-1]) < BLOCK:
                blocks[-1].extend(int(x) for x in t[1:].split(","))
    return blocks


def _mirror_bonds(blocks):
    """Add the forward half of every backward bond record.

    molvector.decode only forms a bond when both atoms record it, while
    mv_back writes each bond once. Records past max_num_bonds are dropped.
    """
    n = len(blocks)
    out = [b[:mv.atom_size] for b in blocks]
    for i, b in enumerate(blocks):
        for btype, off in zip(b[mv.atom_size::2], b[mv.atom_size + 1::2]):
            if not off:
                continue
            j = (i + off) % n
            for src, dst in ((i, j), (j, i)):
                if len(out[src]) < BLOCK:
                    out[src].extend([btype, dst - src])
    return out


def detokenize(tokens, rep):
    """Return a canonical SMILES string or None."""
    try:
        if rep == "smiles":
            mol = Chem.MolFromSmiles("".join(tokens))
        elif rep == "selfies":
            mol = Chem.MolFromSmiles(sf.decoder("".join(tokens)))
        else:
            if rep == "mv_atom":
                blocks = [[int(x) for x in t.split(",")] for t in tokens]
            else:
                blocks = _blocks_from_fields(tokens)
            if rep == "mv_lean":
                # hydrogen count 0: RDKit fills in implicit hydrogens
                blocks = [b[:2] + [0] + b[2:] for b in blocks]
            if rep in ("mv_back", "mv_lean"):
                blocks = _mirror_bonds(blocks)
            v = []
            for b in blocks:
                v.extend(b + [0] * (BLOCK - len(b)))
            if not v:
                return None
            mol = mv.decode(v)
            Chem.SanitizeMol(mol)
            mol = Chem.MolFromSmiles(Chem.MolToSmiles(mol))
        if mol is None or mol.GetNumAtoms() == 0:
            return None
        return Chem.MolToSmiles(mol)
    except Exception:
        return None
