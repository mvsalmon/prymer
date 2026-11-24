# TODO take file with coordinates as input for multi job runs
# TODO handle a design job where the input is a single range, not just a single nucleotide
# TODO add PRIMER_EXPLAIN info to output file?

from prymer.cli import main
from prymer.sequence import reverse_complement, build_breakpoint

if __name__ == "__main__":
    raise SystemExit(main())
