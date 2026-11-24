"""UCSC API access functions"""

import json
import requests

from prymer.exceptions import UCSCError


def build_ucsc_url(genome: str, chrom: str, start: int, end: int) -> str:
    """Construct URL for UCSC API request."""
    return(f"https://api.genome.ucsc.edu/getData/sequence?"
           f"genome={genome};chrom={chrom};start={start};end={end}")

def fetch_sequence(url: str, timeout: int = 10) -> str:
    """Fetch sequence from UCSC API. Returns string containing target sequence."""
    try:
        response = requests.get(url, timeout=timeout)
        response.raise_for_status()
    except requests.RequestException as err:
        raise UCSCError(f"UCSC request failed: {err}\n {json.loads(err.response.text)["error"]}") from err
    return response.json()['dna'].upper()

# SNP checking
def build_ucsc_track_url(genome: str, track: str, chrom: str, start: int, end: int) -> str:
    """Construct URL for UCSC track request."""
    return f"https://api.genome.ucsc.edu/getData/track?genome={genome};track={track};chrom={chrom};start={start};end={end}"

def fetch_snps(url: str, track: str, timeout: int = 10) -> list[dict]:
    """Fetch variant records from UCSC track for given region."""
    try:
        response = requests.get(url, timeout=timeout)
        response.raise_for_status()
    except requests.RequestException as err:
        raise UCSCError(f"UCSC request failed: {err}") from err
    return response.json().get(track, [])
