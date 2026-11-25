"""DNA string manipulation functions."""

from primer3 import p3helpers


# p3helpers only included in v1.2.0+
def reverse_complement(sequence: str) -> str:
    """Wrapper around primer3.p3helpers.reverse_complement()."""
    return p3helpers.reverse_complement(sequence)

def build_breakpoint(start_seq: str, end_seq: str, rev_comp: str | None = None) -> str:
    """Join breakpoint sequences for fusion primers.
    Parameters:
        start_seq: String. Sequence of start (first) breakpoint.
        end_seq: String. Sequence of the end (second) breakpoint.
        rev_comp: String. Controls reverse complementation. Can be 'start', 'end', 'both', or None/empty.

    Returns:
        String of the full breakpoint sequence for primer design.
    """
    if rev_comp in ("start", "both"):
        start_seq = reverse_complement(start_seq)
    if rev_comp in ("end", "both"):
        end_seq = reverse_complement(end_seq)
    return start_seq + end_seq
