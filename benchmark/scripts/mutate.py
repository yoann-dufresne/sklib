#!/usr/bin/env python3
"""Apply random single-base substitutions to a FASTA, for controlled-overlap set-op pairs.

    mutate.py <in.fa> <rate> [seed] > out.fa

Each ACGT base is substituted to a different base with probability <rate> (0..1), deterministically
given <seed>. Headers and non-ACGT characters are preserved; line wrapping is normalized to 70 cols.
At substitution rate p and k-mer length k, a k-mer survives unmutated with prob ~ (1-p)^k, so the
two sets share ~ (1-p)^k of their k-mers — a knob for the intersection ratio.
"""
import random
import sys

A = "ACGT"
_other = {c: [d for d in A if d != c] for c in A}


def main() -> int:
    if len(sys.argv) < 3:
        sys.stderr.write("usage: mutate.py <in.fa> <rate> [seed]\n")
        return 1
    path, rate = sys.argv[1], float(sys.argv[2])
    rng = random.Random(int(sys.argv[3]) if len(sys.argv) > 3 else 0)
    out = sys.stdout

    def flush(seq_chars):
        for i in range(0, len(seq_chars), 70):
            out.write("".join(seq_chars[i:i + 70]))
            out.write("\n")

    seq = []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if seq:
                    flush(seq)
                    seq = []
                out.write(line if line.endswith("\n") else line + "\n")
                continue
            for c in line.strip():
                if c in _other and rng.random() < rate:
                    seq.append(rng.choice(_other[c]))
                else:
                    seq.append(c)
    if seq:
        flush(seq)
    return 0


if __name__ == "__main__":
    sys.exit(main())
