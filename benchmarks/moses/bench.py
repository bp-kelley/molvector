"""MOSES benchmark: the same GPT trained on SMILES, SELFIES and molvector.

  python bench.py prepare --rep all
  python bench.py train   --rep mv_atom
  python bench.py sample  --rep mv_atom
  python bench.py report

See README.md in this directory for the full protocol.
"""
import argparse
import csv
import functools
import json
import math
import os
import random
import time
import urllib.request
from multiprocessing import Pool

import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
from rdkit import Chem, RDLogger

from representations import REPRESENTATIONS, detokenize, tokenize

RDLogger.DisableLog("rdApp.*")

MOSES_URL = "https://media.githubusercontent.com/media/molecularsets/moses/master/data/%s.csv"
PAD, BOS, EOS, UNK = "<pad>", "<bos>", "<eos>", "<unk>"


# ----------------------------------------------------------------- data

def moses_smiles(data_dir, split):
    os.makedirs(data_dir, exist_ok=True)
    path = os.path.join(data_dir, split + ".csv")
    if not os.path.exists(path):
        print("downloading", split)
        urllib.request.urlretrieve(MOSES_URL % split, path)
    with open(path) as fh:
        return [r["SMILES"] for r in csv.DictReader(fh)]


def _tokenize_one(smi, rep, n_random):
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return []
    try:
        out = [tokenize(mol, rep, randomize=False)]
        seen = {" ".join(out[0])}
        for _ in range(n_random * 3):
            if len(out) > n_random:
                break
            t = tokenize(mol, rep, randomize=True)
            if " ".join(t) not in seen:
                seen.add(" ".join(t))
                out.append(t)
        return out
    except Exception:
        return []


def encode_split(smiles, rep, n_random, workers):
    fn = functools.partial(_tokenize_one, rep=rep, n_random=n_random)
    with Pool(workers) as pool:
        return pool.map(fn, smiles, chunksize=512)


def pack(seqs, stoi):
    unk = stoi[UNK]
    ids = [[stoi.get(t, unk) for t in s] for s in seqs]
    lengths = np.array([len(s) for s in ids], dtype=np.int32)
    flat = np.array([i for s in ids for i in s], dtype=np.int32)
    return flat, lengths


def cmd_prepare(a):
    train = moses_smiles(a.data_dir, "train")
    test = moses_smiles(a.data_dir, "test")
    moses_smiles(a.data_dir, "test_scaffolds")
    if a.n_train:
        random.Random(0).shuffle(train)
        train = train[:a.n_train]
    heldout = random.Random(1).sample(test, min(a.n_heldout, len(test)))
    reps = REPRESENTATIONS if a.rep == "all" else [a.rep]
    for rep in reps:
        t0 = time.time()
        out = os.path.join(a.work_dir, rep)
        os.makedirs(out, exist_ok=True)
        tr = encode_split(train, rep, a.augment, a.workers)
        n_fail = sum(1 for x in tr if not x)
        tr_seqs = [s for x in tr for s in x]
        ho_seqs = [x[0] for x in encode_split(heldout, rep, 0, a.workers) if x]
        counts = {}
        for s in tr_seqs:
            for t in s:
                counts[t] = counts.get(t, 0) + 1
        vocab = [PAD, BOS, EOS, UNK] + sorted(counts, key=lambda t: -counts[t])
        stoi = {t: i for i, t in enumerate(vocab)}
        for name, seqs in (("train", tr_seqs), ("heldout", ho_seqs)):
            flat, lengths = pack(seqs, stoi)
            np.savez(os.path.join(out, name + ".npz"), flat=flat, lengths=lengths)
        meta = dict(rep=rep, vocab=vocab, n_train_mols=len(train), n_encode_fail=n_fail,
                    n_train_seqs=len(tr_seqs), augment=a.augment,
                    mean_len=float(np.mean([len(s) for s in tr_seqs])),
                    max_len=max(len(s) for s in tr_seqs))
        with open(os.path.join(out, "meta.json"), "w") as fh:
            json.dump(meta, fh)
        print("%-9s vocab=%d seqs=%d mean_len=%.1f max_len=%d encode_fail=%d (%.0fs)" % (
            rep, len(vocab), len(tr_seqs), meta["mean_len"], meta["max_len"], n_fail, time.time() - t0))


class Seqs:
    def __init__(self, path):
        z = np.load(path)
        self.flat = torch.from_numpy(z["flat"].astype(np.int64))
        self.lengths = z["lengths"]
        self.starts = np.concatenate([[0], np.cumsum(self.lengths)[:-1]])

    def __len__(self):
        return len(self.lengths)

    def batch(self, idx, max_len):
        """BOS + tokens + EOS, right padded with PAD (=0)."""
        L = min(int(self.lengths[idx].max()) + 2, max_len)
        x = torch.zeros(len(idx), L, dtype=torch.long)
        for row, i in enumerate(idx):
            s, n = self.starts[i], self.lengths[i]
            seq = torch.cat([torch.tensor([1]), self.flat[s:s + n], torch.tensor([2])])[:L]
            x[row, :len(seq)] = seq
        return x


# ----------------------------------------------------------------- model

class Block(nn.Module):
    def __init__(self, d, heads, dropout):
        super().__init__()
        self.ln1, self.ln2 = nn.LayerNorm(d), nn.LayerNorm(d)
        self.attn = nn.MultiheadAttention(d, heads, dropout=dropout, batch_first=True)
        self.mlp = nn.Sequential(nn.Linear(d, 4 * d), nn.GELU(), nn.Linear(4 * d, d), nn.Dropout(dropout))

    def forward(self, x, mask):
        h = self.ln1(x)
        x = x + self.attn(h, h, h, attn_mask=mask, need_weights=False, is_causal=True)[0]
        return x + self.mlp(self.ln2(x))


class GPT(nn.Module):
    def __init__(self, vocab, max_len, d=512, layers=8, heads=8, dropout=0.1):
        super().__init__()
        self.max_len = max_len
        self.tok = nn.Embedding(vocab, d)
        self.pos = nn.Embedding(max_len, d)
        self.blocks = nn.ModuleList(Block(d, heads, dropout) for _ in range(layers))
        self.ln = nn.LayerNorm(d)
        self.head = nn.Linear(d, vocab, bias=False)
        self.head.weight = self.tok.weight
        self.apply(self._init)

    @staticmethod
    def _init(m):
        if isinstance(m, (nn.Linear, nn.Embedding)):
            nn.init.normal_(m.weight, std=0.02)
            if getattr(m, "bias", None) is not None:
                nn.init.zeros_(m.bias)

    def forward(self, x):
        T = x.shape[1]
        mask = torch.triu(torch.full((T, T), float("-inf"), device=x.device), 1)
        h = self.tok(x) + self.pos(torch.arange(T, device=x.device))
        for b in self.blocks:
            h = b(h, mask)
        return self.head(self.ln(h))


def device_of(name):
    if name != "auto":
        return torch.device(name)
    if torch.cuda.is_available():
        return torch.device("cuda")
    if torch.backends.mps.is_available():
        return torch.device("mps")
    return torch.device("cpu")


@torch.no_grad()
def heldout_nll(model, data, dev, max_len, bs=500):
    """Mean NLL per molecule in nats (sum over tokens incl. EOS)."""
    model.eval()
    total = 0.0
    for s in range(0, len(data), bs):
        x = data.batch(np.arange(s, min(s + bs, len(data))), max_len).to(dev)
        logits = model(x[:, :-1])
        nll = F.cross_entropy(logits.transpose(1, 2), x[:, 1:], ignore_index=0, reduction="sum")
        total += nll.item()
    model.train()
    return total / len(data)


@torch.no_grad()
def sample(model, n, dev, bs=1000, temperature=1.0):
    model.eval()
    out = []
    while len(out) < n:
        b = min(bs, n - len(out))
        x = torch.ones(b, 1, dtype=torch.long, device=dev)
        done = torch.zeros(b, dtype=torch.bool, device=dev)
        for _ in range(model.max_len - 1):
            logits = model(x)[:, -1] / temperature
            logits[:, 0] = -float("inf")  # never PAD
            logits[:, 1] = -float("inf")  # never BOS
            nxt = torch.multinomial(F.softmax(logits, -1), 1).squeeze(1)
            nxt[done] = 0
            x = torch.cat([x, nxt[:, None]], 1)
            done |= nxt == 2
            if done.all():
                break
        out.extend(x[:, 1:].tolist())
    model.train()
    return out


def ids_to_tokens(ids, vocab):
    toks = []
    for i in ids:
        if i in (0, 2):
            break
        toks.append(vocab[i])
    return toks


def quick_validity(seqs, vocab, rep):
    smis = [detokenize(ids_to_tokens(s, vocab), rep) for s in seqs]
    valid = [s for s in smis if s]
    return len(valid) / len(smis), len(set(valid)) / max(1, len(valid))


def cmd_train(a):
    d = os.path.join(a.work_dir, a.rep)
    meta = json.load(open(os.path.join(d, "meta.json")))
    vocab = meta["vocab"]
    train, heldout = Seqs(os.path.join(d, "train.npz")), Seqs(os.path.join(d, "heldout.npz"))
    max_len = meta["max_len"] + 2 + a.len_margin
    dev = device_of(a.device)
    torch.manual_seed(a.seed)
    model = GPT(len(vocab), max_len, a.d_model, a.layers, a.heads, a.dropout).to(dev)
    n_params = sum(p.numel() for p in model.parameters())
    opt = torch.optim.AdamW(model.parameters(), lr=a.lr, weight_decay=0.01, betas=(0.9, 0.95))
    sched = torch.optim.lr_scheduler.LambdaLR(opt, lambda s: min(1, (s + 1) / a.warmup) * 0.5 * (
        1 + math.cos(math.pi * min(1.0, s / a.steps))))
    amp = dev.type == "cuda"
    log_path = os.path.join(d, "train_log.jsonl")
    log = open(log_path, "w")
    print("%s: %d params, vocab %d, max_len %d, device %s" % (a.rep, n_params, len(vocab), max_len, dev))
    rng = np.random.default_rng(a.seed)
    t0 = time.time()
    tokens_seen = 0
    for step in range(a.steps + 1):
        if step % a.eval_every == 0 or step == a.steps:
            rec = dict(step=step, time=time.time() - t0, tokens=tokens_seen,
                       heldout_nll=heldout_nll(model, heldout, dev, max_len), params=n_params)
            if a.eval_samples:
                rec["valid"], rec["unique"] = quick_validity(sample(model, a.eval_samples, dev), vocab, a.rep)
            log.write(json.dumps(rec) + "\n")
            log.flush()
            print(json.dumps(rec))
            if step == a.steps:
                break
        x = train.batch(rng.integers(0, len(train), a.batch_size), max_len).to(dev)
        tokens_seen += int((x[:, 1:] != 0).sum())
        with torch.autocast(dev.type, dtype=torch.bfloat16, enabled=amp):
            logits = model(x[:, :-1])
            loss = F.cross_entropy(logits.float().transpose(1, 2), x[:, 1:], ignore_index=0)
        opt.zero_grad(set_to_none=True)
        loss.backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
        opt.step()
        sched.step()
    torch.save(dict(model=model.state_dict(), args=vars(a), max_len=max_len, vocab_size=len(vocab)),
               os.path.join(d, "model.pt"))


def cmd_sample(a):
    d = os.path.join(a.work_dir, a.rep)
    meta = json.load(open(os.path.join(d, "meta.json")))
    ck = torch.load(os.path.join(d, "model.pt"), map_location="cpu")
    ta = ck["args"]
    dev = device_of(a.device)
    model = GPT(ck["vocab_size"], ck["max_len"], ta["d_model"], ta["layers"], ta["heads"], 0.0)
    model.load_state_dict(ck["model"])
    model.to(dev)
    torch.manual_seed(a.seed)
    t0 = time.time()
    seqs = sample(model, a.n, dev, a.sample_batch, a.temperature)
    t_sample = time.time() - t0
    smis = [detokenize(ids_to_tokens(s, meta["vocab"]), a.rep) for s in seqs]
    with open(os.path.join(d, "samples.txt"), "w") as fh:
        fh.write("\n".join(s or "" for s in smis))
    res = score(smis, a.data_dir, a.fcd_device or str(dev))
    res.update(rep=a.rep, sample_seconds=t_sample, n_params=sum(p.numel() for p in model.parameters()))
    json.dump(res, open(os.path.join(d, "results.json"), "w"), indent=1)
    print(json.dumps(res, indent=1))


def cmd_score(a):
    """Re-score an existing samples.txt without sampling again."""
    d = os.path.join(a.work_dir, a.rep)
    with open(os.path.join(d, "samples.txt")) as fh:
        smis = [s or None for s in fh.read().split("\n")]
    ck = torch.load(os.path.join(d, "model.pt"), map_location="cpu")
    res = score(smis, a.data_dir, a.fcd_device or "cpu")
    res.update(rep=a.rep, n_params=sum(t.numel() for k, t in ck["model"].items() if k != "head.weight"))
    json.dump(res, open(os.path.join(d, "results.json"), "w"), indent=1)
    print(json.dumps(res, indent=1))


def frechet_distance(p, q, eps=1e-6):
    """FCD from ChemNet statistics, as fcd_torch computes it.

    Done here because fcd_torch passes sqrtm(disp=False), which newer SciPy
    no longer accepts.
    """
    from scipy import linalg
    if not p or not q:
        return float("nan")
    diff = p["mu"] - q["mu"]
    covmean = linalg.sqrtm(p["sigma"].dot(q["sigma"]))
    if not np.isfinite(covmean).all():
        offset = np.eye(p["sigma"].shape[0]) * eps
        covmean = linalg.sqrtm((p["sigma"] + offset).dot(q["sigma"] + offset))
    covmean = np.real(covmean)
    return float(diff.dot(diff) + np.trace(p["sigma"]) + np.trace(q["sigma"]) - 2 * np.trace(covmean))


def score(smis, data_dir, fcd_device):
    valid = [s for s in smis if s]
    uniq = set(valid)
    train = canonical_train(data_dir)
    res = dict(n=len(smis),
               valid=len(valid) / len(smis),
               valid_single_fragment=sum(1 for s in valid if "." not in s) / len(smis),
               unique_at_1k=len(set(valid[:1000])) / max(1, len(valid[:1000])),
               unique_at_10k=len(set(valid[:10000])) / max(1, len(valid[:10000])),
               novelty=len(uniq - train) / max(1, len(uniq)))
    if not hasattr(np, "row_stack"):
        np.row_stack = np.vstack  # removed in NumPy 2.4, still used by fcd_torch
    try:
        from fcd_torch import FCD
        fcd = FCD(device=fcd_device, n_jobs=1)
        gen = list(uniq)[:10000] if len(uniq) >= 10000 else list(uniq)
        for split in ("test", "test_scaffolds"):
            ref = random.Random(0).sample(moses_smiles(data_dir, split), 10000)
            res["fcd_" + split] = frechet_distance(fcd.precalc(ref), fcd.precalc(gen))
    except ImportError:
        print("fcd_torch not installed; skipping FCD")
    return res


def canonical_train(data_dir):
    """MOSES train SMILES canonicalized with the installed RDKit (cached)."""
    path = os.path.join(data_dir, "train_canonical.txt")
    if not os.path.exists(path):
        with Pool(os.cpu_count()) as pool:
            can = pool.map(_canonical, moses_smiles(data_dir, "train"), chunksize=2048)
        with open(path, "w") as fh:
            fh.write("\n".join(c for c in can if c))
    return set(open(path).read().split("\n"))


def _canonical(smi):
    mol = Chem.MolFromSmiles(smi)
    return Chem.MolToSmiles(mol) if mol else None


def cmd_report(a):
    rows = []
    for rep in REPRESENTATIONS:
        d = os.path.join(a.work_dir, rep)
        if not os.path.exists(os.path.join(d, "results.json")):
            continue
        r = json.load(open(os.path.join(d, "results.json")))
        meta = json.load(open(os.path.join(d, "meta.json")))
        logs = [json.loads(l) for l in open(os.path.join(d, "train_log.jsonl"))]
        r.update(vocab=len(meta["vocab"]), mean_len=meta["mean_len"], final_nll=logs[-1]["heldout_nll"],
                 train_seconds=logs[-1]["time"])
        thr = [l["step"] for l in logs if l.get("valid", 0) >= 0.9]
        r["step_valid_90"] = thr[0] if thr else None
        rows.append(r)
    cols = ["rep", "vocab", "mean_len", "n_params", "final_nll", "step_valid_90", "valid",
            "valid_single_fragment", "unique_at_10k", "novelty", "fcd_test", "fcd_test_scaffolds", "train_seconds"]
    lines = ["| " + " | ".join(cols) + " |", "|" + "---|" * len(cols)]
    for r in rows:
        lines.append("| " + " | ".join(
            ("%.3f" % r[c] if isinstance(r.get(c), float) else str(r.get(c, ""))) for c in cols) + " |")
    text = "\n".join(lines)
    open(os.path.join(a.work_dir, "report.md"), "w").write(text + "\n")
    print(text)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("cmd", choices=["prepare", "train", "sample", "score", "report"])
    p.add_argument("--rep", default="all", choices=("all",) + REPRESENTATIONS)
    p.add_argument("--data-dir", default="data")
    p.add_argument("--work-dir", default="runs")
    p.add_argument("--device", default="auto")
    p.add_argument("--seed", type=int, default=0)
    # prepare
    p.add_argument("--n-train", type=int, default=0, help="subsample MOSES train (0 = all 1.58M)")
    p.add_argument("--n-heldout", type=int, default=10000)
    p.add_argument("--augment", type=int, default=0, help="random atom orders per molecule, on top of canonical")
    p.add_argument("--workers", type=int, default=os.cpu_count())
    # train
    p.add_argument("--steps", type=int, default=50000)
    p.add_argument("--batch-size", type=int, default=256)
    p.add_argument("--lr", type=float, default=3e-4)
    p.add_argument("--warmup", type=int, default=1000)
    p.add_argument("--d-model", type=int, default=512)
    p.add_argument("--layers", type=int, default=8)
    p.add_argument("--heads", type=int, default=8)
    p.add_argument("--dropout", type=float, default=0.1)
    p.add_argument("--len-margin", type=int, default=8)
    p.add_argument("--eval-every", type=int, default=1000)
    p.add_argument("--eval-samples", type=int, default=1000)
    # sample
    p.add_argument("--n", type=int, default=30000)
    p.add_argument("--sample-batch", type=int, default=1000)
    p.add_argument("--temperature", type=float, default=1.0)
    p.add_argument("--fcd-device", default=None)
    a = p.parse_args()
    if a.cmd != "prepare" and a.cmd != "report" and a.rep == "all":
        for rep in REPRESENTATIONS:
            a.rep = rep
            globals()["cmd_" + a.cmd](a)
        return
    globals()["cmd_" + a.cmd](a)


if __name__ == "__main__":
    main()
