#!/usr/bin/env python3
"""Helpers for scripts/large_scale_e2e.sh.

Subcommands:
  sanitize <fasta>             Uppercase + split each record at runs of non-ACGT
                               characters into pure-ACGT sub-records (stdout).
                               Makes sklib and KMC agree (sklib silently encodes
                               N as G, KMC drops N-containing k-mers).
  randkmers <k> <n> <seed>     Print n random ACGT k-mers (deterministic seed).
  canon <kmers_file>           Print canonical form min(fwd, revcomp) per line,
                               preserving order (matches KMC's canonical k-mers).
  truth <present_set> <qcanon> Print 0/1 per query line: 1 iff its canonical form
                               is in <present_set>, in order.
  prefix <fasta> <max_bp>      Emit the first <max_bp> ACGT bases of the FASTA as
                               records, never crossing a record boundary (so every
                               emitted k-mer is a real genome k-mer). Used to take a
                               bounded but structure-preserving self-query sample.
  bincount <list.sskm>         Print the #skmers stored in the binary header.
  binheadersize <list.sskm>    Print the header byte size (32 for V2; 40+16*n_buckets for V3).
  dropshort <fasta> <k>        Emit only FASTA records with sequence length >= k (feeds CBL).

Benchmark query-workload helpers (scripts/bench/bench.sh):
  sample_positive <fa> <k> <n> <seed>
                               Print n k-mers (one per line) drawn uniformly at random
                               from the positions of <fa> -- every emitted k-mer is
                               guaranteed present (positive query workload, scattered).
  simreads <fa> <len> <n> <err> <seed>
                               Emit n FASTA reads of <len> bp sampled from <fa> with a
                               per-base substitution error rate <err> in [0,1]
                               (realistic mixed-status query workload).
  shuffle <file> <seed>        Print the lines of <file> in a deterministic random order
                               (used to break query locality: streaming vs shuffled).
  count_kmers <fa> <k>         Print the total number of k-mers over all records
                               (throughput denominator).

Note: sample_positive/simreads load the FASTA into memory (~genome size); fine up to a
human genome on this box. For >RAM inputs, use a real read set instead.
"""
import sys
import random
import re
import struct

_COMP = str.maketrans("ACGT", "TGCA")


def _canon(s):
    r = s.translate(_COMP)[::-1]
    return s if s <= r else r


def _write_fasta(w, header, seq, width=60):
    # Wrap to fixed-width lines: multi-megabyte single lines choke grep/awk and
    # some FASTA readers; 60 nt/line is the conventional, safe width.
    w(header + "\n")
    for j in range(0, len(seq), width):
        w(seq[j:j + width])
        w("\n")


def sanitize(path):
    w = sys.stdout.write
    name = None
    parts = []

    def flush():
        if name is None:
            return
        seq = "".join(parts).upper()
        i = 0
        for frag in re.split(r"[^ACGT]+", seq):
            if frag:
                i += 1
                _write_fasta(w, ">%s_%d" % (name, i), frag)

    with open(path) as f:
        for line in f:
            if line[:1] == ">":
                flush()
                hdr = line[1:].strip()
                name = hdr.split()[0] if hdr else "seq"
                parts = []
            else:
                parts.append(line.strip())
    flush()


def randkmers(k, n, seed):
    random.seed(seed)
    w = sys.stdout.write
    buf = []
    for _ in range(n):
        buf.append("".join(random.choice("ACGT") for _ in range(k)))
        if len(buf) >= 65536:
            w("\n".join(buf) + "\n")
            buf = []
    if buf:
        w("\n".join(buf) + "\n")


def canon(path):
    w = sys.stdout.write
    with open(path) as f:
        for line in f:
            s = line.strip()
            if s:
                w(_canon(s) + "\n")


def truth(present_path, qcanon_path):
    present = set()
    with open(present_path) as f:
        for line in f:
            x = line.strip()
            if x:
                present.add(x)
    w = sys.stdout.write
    with open(qcanon_path) as f:
        for line in f:
            x = line.strip()
            if x:
                w(("1" if x in present else "0") + "\n")


def prefix(path, max_bp):
    max_bp = int(max_bp)
    w = sys.stdout.write
    emitted = 0
    name = None
    buf = []
    buflen = 0

    def flush():
        nonlocal emitted
        if name is None or emitted >= max_bp:
            return
        seq = "".join(buf)
        take = min(len(seq), max_bp - emitted)
        if take > 0:
            _write_fasta(w, ">%s" % name, seq[:take])
            emitted += take

    with open(path) as f:
        for line in f:
            if emitted >= max_bp:
                break
            if line[:1] == ">":
                flush()
                hdr = line[1:].strip()
                name = hdr.split()[0] if hdr else "seq"
                buf = []
                buflen = 0
            else:
                need = max_bp - emitted
                if buflen < need:
                    s = line.strip()
                    buf.append(s)
                    buflen += len(s)
    flush()


def bincount(path):
    with open(path, "rb") as f:
        hdr = f.read(32)
    if len(hdr) < 32:
        sys.exit("file too short to hold a VirtualSkmer header: " + path)
    _magic, _k, _m, count = struct.unpack("<4Q", hdr)
    print(count)


def binheadersize(path):
    # Byte size of the on-disk header (everything before the skmer payload), so callers can
    # compute the true payload size. VSKMER_2: magic+k+m+count = 32 B. VSKMER_3 adds
    # n_buckets(8) + a per-bucket directory of 16 B each.
    V2 = 0x56534B4D45525F32  # "VSKMER_2"
    V3 = 0x56534B4D45525F33  # "VSKMER_3"
    with open(path, "rb") as f:
        hdr = f.read(40)
    if len(hdr) < 32:
        sys.exit("file too short to hold a VirtualSkmer header: " + path)
    magic = struct.unpack("<Q", hdr[:8])[0]
    if magic == V2:
        print(32)
    elif magic == V3:
        if len(hdr) < 40:
            sys.exit("file too short to hold a VSKMER_3 header: " + path)
        n_buckets = struct.unpack("<Q", hdr[32:40])[0]
        print(40 + 16 * n_buckets)
    else:
        sys.exit("unknown VirtualSkmer magic in: " + path)


def _iter_records(path):
    """Yield the concatenated sequence of each FASTA record (uppercased)."""
    name = None
    parts = []
    with open(path) as f:
        for line in f:
            if line[:1] == ">":
                if name is not None:
                    yield "".join(parts).upper()
                name = line[1:].strip()
                parts = []
            else:
                parts.append(line.strip())
    if name is not None:
        yield "".join(parts).upper()


def dropshort(path, k):
    # Emit only records whose sequence length >= k. Records shorter than k hold no
    # k-mers, so dropping them is k-mer-set-neutral. Used to feed CBL, which aborts on
    # any sequence shorter than K (e.g. 1-bp fragments left between N-runs by sanitize),
    # unlike the other tools that simply skip short records.
    i = 0
    for seq in _iter_records(path):
        if len(seq) >= k:
            sys.stdout.write(">%d\n%s\n" % (i, seq))
            i += 1


def _load_seqs(path, min_len):
    """Load records >= min_len into memory; return (seqs, cumulative_positions).

    cumulative[i] is the number of valid start positions in seqs[0..i-1], so the
    total over all records is cumulative[-1] and a global index g in
    [0, cumulative[-1]) maps to (record, offset) by bisecting cumulative.
    """
    seqs = []
    cumulative = [0]
    total = 0
    for seq in _iter_records(path):
        nstarts = len(seq) - min_len + 1
        if nstarts > 0:
            seqs.append(seq)
            total += nstarts
            cumulative.append(total)
    return seqs, cumulative


def _global_to_local(g, cumulative):
    """Map global start index g to (record_index, offset) via the cumulative table."""
    import bisect
    rec = bisect.bisect_right(cumulative, g) - 1
    return rec, g - cumulative[rec]


def sample_positive(path, k, n, seed):
    seqs, cumulative = _load_seqs(path, k)
    total = cumulative[-1]
    if total == 0:
        sys.exit("no %d-mer in %s" % (k, path))
    rng = random.Random(seed)
    w = sys.stdout.write
    buf = []
    for _ in range(n):
        rec, off = _global_to_local(rng.randrange(total), cumulative)
        buf.append(seqs[rec][off:off + k])
        if len(buf) >= 65536:
            w("\n".join(buf) + "\n")
            buf = []
    if buf:
        w("\n".join(buf) + "\n")


def simreads(path, read_len, n, err, seed):
    seqs, cumulative = _load_seqs(path, read_len)
    total = cumulative[-1]
    if total == 0:
        sys.exit("no read of length %d fits in %s" % (read_len, path))
    rng = random.Random(seed)
    others = {"A": "CGT", "C": "AGT", "G": "ACT", "T": "ACG"}
    w = sys.stdout.write
    for i in range(n):
        rec, off = _global_to_local(rng.randrange(total), cumulative)
        read = list(seqs[rec][off:off + read_len])
        if err > 0:
            for j in range(read_len):
                if rng.random() < err:
                    read[j] = rng.choice(others[read[j]])
        _write_fasta(w, ">read_%d" % i, "".join(read))


def shuffle(path, seed):
    with open(path) as f:
        lines = f.read().splitlines()
    random.Random(seed).shuffle(lines)
    sys.stdout.write("\n".join(lines) + ("\n" if lines else ""))


def count_kmers(path, k):
    total = 0
    for seq in _iter_records(path):
        nk = len(seq) - k + 1
        if nk > 0:
            total += nk
    print(total)


def allkmers(path, k):
    """Emit every k-mer of every record, one per line (streaming, low memory)."""
    w = sys.stdout.write
    buf = []
    for seq in _iter_records(path):
        for j in range(len(seq) - k + 1):
            buf.append(seq[j:j + k])
            if len(buf) >= 65536:
                w("\n".join(buf) + "\n")
                buf = []
    if buf:
        w("\n".join(buf) + "\n")


def main():
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    cmd = sys.argv[1]
    if cmd == "sanitize":
        sanitize(sys.argv[2])
    elif cmd == "randkmers":
        randkmers(int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]))
    elif cmd == "canon":
        canon(sys.argv[2])
    elif cmd == "truth":
        truth(sys.argv[2], sys.argv[3])
    elif cmd == "prefix":
        prefix(sys.argv[2], sys.argv[3])
    elif cmd == "bincount":
        bincount(sys.argv[2])
    elif cmd == "binheadersize":
        binheadersize(sys.argv[2])
    elif cmd == "dropshort":
        dropshort(sys.argv[2], int(sys.argv[3]))
    elif cmd == "sample_positive":
        sample_positive(sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5]))
    elif cmd == "simreads":
        simreads(sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), float(sys.argv[5]), int(sys.argv[6]))
    elif cmd == "shuffle":
        shuffle(sys.argv[2], int(sys.argv[3]))
    elif cmd == "count_kmers":
        count_kmers(sys.argv[2], int(sys.argv[3]))
    elif cmd == "allkmers":
        allkmers(sys.argv[2], int(sys.argv[3]))
    else:
        sys.exit("unknown subcommand: %s" % cmd)


if __name__ == "__main__":
    main()
