#!/usr/bin/env bash
# Stage 1 Step 3 — AlphaFold2 on the BINDER ALONE (foldability / self-consistency filter).
# For every MPNN sequence: thread it onto its dimer-arm binder backbone (chain B), emit the
# binder alone (as chain A), batch-predict with af2_initial_guess in monomer mode, and score by
# plddt_binder (fold confidence) + binder_aligned_rmsd (does it fold back to the RFd design?).
# This does NOT test binding — the interface energy check (Rosetta, Step 4) does. It prunes the
# 1600 MPNN sequences down to a foldable shortlist to carry forward.
#
# Default run splits the work across BOTH GPUs (backbones sharded by design-number parity),
# then merges the two scorefiles and ranks.
#
# Input : <campaign>/02_mpnn/seqs/design_*.fa  +  <campaign>/01_rfd/design_*.pdb
# Output: <campaign>/03_af2/{in,out}_s{0,1}/ , scores_s{0,1}.sc , ranked_binderonly.csv
#         where <campaign> = outputs/C3_symmetric_Binder_2026_07_02/ (one campaign = one folder)
# Usage : bash scripts/22_af2_dimer_arm.sh            # two-GPU orchestrator (GPU0+GPU1)
#         SHARD=0 NSHARD=2 GPU=0 bash scripts/22_af2_dimer_arm.sh _worker   # one shard (internal)
#         bash scripts/22_af2_dimer_arm.sh rank       # merge shard scores + rank only
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
source "$HERE/lib/common.sh"

BCHAIN="${BCHAIN:-B}"                       # binder chain in the dimer-arm backbones
CAMPAIGN="$PROJECT/outputs/C3_symmetric_Binder_2026_07_02"   # one campaign = one folder
SEQDIR="${SEQDIR:-$CAMPAIGN/02_mpnn/seqs}"   # override for round 2 (06_mpnn_buns/seqs)
BBDIR="${BBDIR:-$CAMPAIGN/01_rfd}"
OUT="${OUT:-$CAMPAIGN/03_af2}"               # override for round 2 (07_af2)
export SEQDIR BBDIR OUT                       # propagate overrides to _worker child shells
LOGDIR="$CAMPAIGN/logs"; mkdir -p "$OUT" "$LOGDIR"
MODE="${1:-orchestrate}"

# ---------------------------------------------------------------- merge + rank
if [ "$MODE" = "rank" ]; then
  echo "→ Merging shard scorefiles + ranking by plddt/rmsd..."
  "$IG_PY" - "$OUT/ranked_binderonly.csv" "$OUT"/scores_s*.sc <<'PY'
import sys, csv, re, glob
outcsv = sys.argv[1]; scfiles = sys.argv[2:]
rows=[]
for sc in scfiles:
    hdr=None
    for line in open(sc):
        if not line.startswith("SCORE:"): continue
        t=line.split()[1:]
        if t[0]=="binder_aligned_rmsd": hdr=t; continue
        if hdr and len(t)==len(hdr):
            d=dict(zip(hdr,t)); tag=re.sub(r"_af2pred$","",d["description"])
            m=re.match(r"(design_\d+)_s\d+",tag); bb=m.group(1) if m else tag
            rows.append({"tag":tag,"backbone":bb,
                "plddt_binder":round(float(d["plddt_binder"]),2),
                "binder_rmsd":round(float(d["binder_aligned_rmsd"]),2),
                "foldable":bool(float(d["plddt_binder"])>80 and float(d["binder_aligned_rmsd"])<2.0)})
rows.sort(key=lambda r:(not r["foldable"], -r["plddt_binder"], r["binder_rmsd"]))
keys=["rank","tag","backbone","plddt_binder","binder_rmsd","foldable"]
with open(outcsv,"w",newline="") as fh:
    w=csv.DictWriter(fh,fieldnames=keys,lineterminator="\n"); w.writeheader()
    for i,r in enumerate(rows,1): w.writerow({"rank":i,**r})
nf=sum(r["foldable"] for r in rows)
nbb=len({r["backbone"] for r in rows if r["foldable"]})
print(f"  {nf}/{len(rows)} sequences foldable (plddt>80 AND rmsd<2A), spanning {nbb} backbones -> {outcsv}")
for r in rows[:10]:
    print(f"  {r['tag']}: plddt={r['plddt_binder']} rmsd={r['binder_rmsd']} foldable={r['foldable']}")
PY
  echo "→ Done. Ranked: $OUT/ranked_binderonly.csv"
  exit 0
fi

# ---------------------------------------------------------------- one shard worker
if [ "$MODE" = "_worker" ]; then
  SHARD="${SHARD:-0}"; NSHARD="${NSHARD:-2}"; GPU="${GPU:-0}"
  INDIR="$OUT/in_s${SHARD}"; PREDOUT="$OUT/out_s${SHARD}"; SC="$OUT/scores_s${SHARD}.sc"
  mkdir -p "$INDIR" "$PREDOUT"
  echo "[shard $SHARD/$NSHARD, GPU$GPU] threading backbones (design_num % $NSHARD == $SHARD)..."
  "$IG_PY" - "$SEQDIR" "$BBDIR" "$INDIR" "$BCHAIN" "$SHARD" "$NSHARD" <<'PY'
import sys, glob, os, re
from pathlib import Path
ONE_TO_THREE = {"A":"ALA","R":"ARG","N":"ASN","D":"ASP","C":"CYS","Q":"GLN","E":"GLU",
"G":"GLY","H":"HIS","I":"ILE","L":"LEU","K":"LYS","M":"MET","F":"PHE","P":"PRO",
"S":"SER","T":"THR","W":"TRP","Y":"TYR","V":"VAL"}
seqdir, bbdir, indir, bchain, shard, nshard = sys.argv[1:7]
shard=int(shard); nshard=int(nshard)
def read_seqs(fa):  # record 0 is the native poly-G placeholder; return the designed seqs in order
    recs=[]; cur=None
    for line in open(fa):
        if line.startswith(">"): cur=[]; recs.append(cur)
        elif cur is not None: cur.append(line.strip())
    seqs=["".join(r) for r in recs]
    return seqs[1:] if len(seqs)>1 else []
n=0
for fa in sorted(glob.glob(os.path.join(seqdir,"design_*.fa"))):
    name=os.path.basename(fa)[:-3]
    num=int(re.match(r"design_(\d+)",name).group(1))
    if num % nshard != shard: continue
    bb=os.path.join(bbdir,name+".pdb")
    if not os.path.exists(bb): continue
    # Only chain B is designed, so each designed FASTA record IS the 120-aa binder sequence.
    atoms=[l for l in open(bb) if l.startswith("ATOM") and l[21]==bchain]
    order=[]; seen=set()
    for l in atoms:
        if l[22:26] not in seen: seen.add(l[22:26]); order.append(l[22:26])
    ridx={r:i for i,r in enumerate(order)}
    for i,seq in enumerate(read_seqs(fa)):
        if len(seq)!=len(order): continue
        out=[]; ser=0
        for l in atoms:
            ri=ridx[l[22:26]]; ser+=1
            out.append(f"{l[:6]}{ser:5d}{l[11:17]}{ONE_TO_THREE[seq[ri]]:>3s} A{ri+1:4d}{l[26:]}")
        out+=["TER","END"]
        Path(os.path.join(indir,f"{name}_s{i}.pdb")).write_text("\n".join(out)+"\n")
        n+=1
print(f"  [shard {shard}] threaded {n} binder-only structures")
PY
  NIN=$(find "$INDIR" -name '*.pdb' | wc -l)
  echo "[shard $SHARD, GPU$GPU] af2_initial_guess (monomer) over $NIN structures..."
  cd "$IG_DIR"
  CUDA_VISIBLE_DEVICES=$GPU "$IG_PY" predict.py \
    -pdbdir "$INDIR" -outpdbdir "$PREDOUT" -scorefilename "$SC" \
    -checkpoint_name "$OUT/check_s${SHARD}.point"
  echo "[shard $SHARD, GPU$GPU] done."
  exit 0
fi

# ---------------------------------------------------------------- orchestrate (default): both GPUs
echo "→ Two-GPU run: shard 0 → GPU0, shard 1 → GPU1."
SHARD=0 NSHARD=2 GPU=0 bash "$0" _worker > "$LOGDIR/run_af2_s0.log" 2>&1 &
P0=$!
SHARD=1 NSHARD=2 GPU=1 bash "$0" _worker > "$LOGDIR/run_af2_s1.log" 2>&1 &
P1=$!
echo "  shard0 PID $P0 (GPU0), shard1 PID $P1 (GPU1). Waiting..."
wait $P0; R0=$?
wait $P1; R1=$?
echo "  shard0 exit=$R0, shard1 exit=$R1"
if [ "$R0" -ne 0 ] || [ "$R1" -ne 0 ]; then
  echo "!! a shard failed — see run_af2_s0.log / run_af2_s1.log"; exit 1
fi
bash "$0" rank
