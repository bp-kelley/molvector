"""Tokenizers for the MOSES benchmark.

Each representation turns an RDKit molecule into a list of string tokens
and turns a list of tokens back into a SMILES string (or None when the
tokens do not decode to a sanitizable molecule).

  smiles    regex-tokenized SMILES (Schwaller et al. tokenizer)
  selfies   SELFIES symbols
  mv_atom   molvector, one composite token per atom block
  mv_field  molvector, one token for the atom fields plus one token
            per non-empty (bond type, offset) record
"""
import os
import re
import sys

import selfies as sf
from rdkit import Chem

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
import molvector as mv  # noqa: E402

REPRESENTATIONS = ("smiles", "selfies", "mv_atom", "mv_field")

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
    v = molvector_of(mol, randomize)
    blocks = [v[i:i + BLOCK] for i in range(0, len(v), BLOCK)]
    if rep == "mv_atom":
        return [",".join(map(str, b)) for b in blocks]
    if rep == "mv_field":
        toks = []
        for b in blocks:
            toks.append("A%d,%d,%d" % tuple(b[:mv.atom_size]))
            for k in range(mv.max_num_bonds):
                btype, off = b[mv.atom_size + 2 * k], b[mv.atom_size + 2 * k + 1]
                if off:
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
