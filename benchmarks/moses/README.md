# MOSES benchmark: SMILES vs SELFIES vs molvector

Trains the same decoder-only transformer on MOSES with four tokenizations
and scores the samples with the standard MOSES metrics. The question it
answers is whether molvector is worth continuing: does a model learn it
faster, or produce better molecules, than it learns SMILES or SELFIES?

| rep        | one token is                                         | mean tokens / mol |
|------------|------------------------------------------------------|-------------------|
| `smiles`   | a SMILES regex token                                 | ~35               |
| `selfies`  | a SELFIES symbol                                     | ~35               |
| `mv_atom`  | a whole molvector atom block (13 ints)               | ~22               |
| `mv_field` | the atom fields, or one `(bond type, offset)` record | ~68               |
| `mv_back`  | as `mv_field`, but each bond written once (to an earlier atom) | ~45     |
| `mv_lean`  | as `mv_back`, on the Kekulé form with no hydrogen-count field  | ~45     |

All four use the canonical RDKit atom order by default. `--augment K` adds
K random-order variants per molecule (randomized SMILES/SELFIES, or random
molvector atom orders), which is how molvector was meant to be trained.

## Run

```bash
pip install torch rdkit selfies fcd_torch numpy
cd benchmarks/moses

python bench.py prepare --rep all            # downloads MOSES into ./data, writes ./runs/<rep>/
python bench.py train   --rep all            # 50k steps, batch 256, 8x512 GPT, one rep after another
python bench.py sample  --rep all            # 30k samples each, MOSES metrics + FCD
python bench.py report                        # writes runs/report.md
```

Each rep is independent, so on several GPUs run them in parallel:

```bash
python bench.py prepare --rep all
i=0; for r in smiles selfies mv_atom mv_field mv_back mv_lean; do
  (export CUDA_VISIBLE_DEVICES=$i; python bench.py train --rep $r && python bench.py sample --rep $r) &
  i=$((i+1))
done; wait; python bench.py report
```

`python bench.py --help` lists every knob (model size, steps, learning rate,
augmentation, eval cadence).

## What gets measured

During training, every `--eval-every` steps, `runs/<rep>/train_log.jsonl` records:

- `heldout_nll`: mean negative log-likelihood per molecule (nats, summed over
  tokens) on held-out MOSES test molecules. Because each rep encodes a molecule
  deterministically in canonical order, this is comparable across reps and is
  the cleanest convergence curve.
- `valid` / `unique` on `--eval-samples` quick samples, so you can see how
  many steps each rep needs to reach a given validity (`step_valid_90` in the
  report).

After sampling, `runs/<rep>/results.json` holds MOSES validity, validity of
single-fragment molecules, uniqueness@1k/@10k, novelty against the train set,
and FCD against MOSES Test and TestSF.

## Things to know when reading the results

- `valid_single_fragment` matters for molvector. Its decoder always yields a
  graph, but that graph can be disconnected, and MOSES `valid` counts
  disconnected molecules as valid.
- `mv_atom` has a large vocabulary (about 6.6k blocks for 200k canonical
  molecules, more with augmentation), so its embedding table adds parameters.
  The report lists parameter counts.
- molvector writes each bond twice, once from each end. A model can emit
  contradictory pairs, and in a small CPU run this caused 80 to 95% of invalid
  molvector samples. `mv_back` removes the duplicate to test that directly.
- molvector currently drops stereochemistry, and so does this benchmark for
  every rep (MOSES is mostly stereo-free).
