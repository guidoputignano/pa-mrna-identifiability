#!/usr/bin/env python3
"""
Build paper/refs.bib from paper/refs.sources.

Design rule: nothing is invented. Each entry is fetched from the
registrar, either CrossRef (via doi.org content negotiation) or NCBI
(for PubMed IDs). Anything that fails to resolve is written to
paper/refs.todo with the reason, and is NOT written to refs.bib.

Hand-written entries that have no DOI (GeneReviews, standards, theses)
live in paper/refs.manual.bib and are concatenated unchanged.

Usage:
    python tools/fetch_refs.py                 # normal run
    python tools/fetch_refs.py --check         # exit 1 if anything unresolved
    python tools/fetch_refs.py --force         # ignore the cache
"""

import argparse
import json
import re
import sys
import time
from pathlib import Path

import requests

ROOT = Path(__file__).resolve().parents[1]
SOURCES = ROOT / "paper" / "refs.sources"
OUT = ROOT / "paper" / "refs.bib"
TODO = ROOT / "paper" / "refs.todo"
MANUAL = ROOT / "paper" / "refs.manual.bib"
CACHE = ROOT / ".refcache.json"

#: the line that separates refs.bib's fetched half from its hand-written one.
#: It is a marker, not decoration: the offline path below splits on it.
MANUAL_MARK = "\n% ---- hand-written entries, from paper/refs.manual.bib ----\n"

MAILTO = "you@example.org"        # CrossRef asks for a contact; be polite
UA = f"pa-identifiability-refs/1.0 (mailto:{MAILTO})"
PAUSE = 0.4                        # seconds between requests


# ------------------------------------------------------------------ input
def read_sources(path):
    """Parse 'citekey  doi:... | pmid:...' lines, ignoring comments."""
    entries = []
    for lineno, raw in enumerate(path.read_text().splitlines(), 1):
        line = raw.split("#", 1)[0].strip()
        if not line:
            continue
        parts = line.split()
        if len(parts) != 2 or ":" not in parts[1]:
            print(f"  ! {path.name}:{lineno}: cannot parse -> {raw.strip()!r}")
            continue
        key, ident = parts
        kind, _, value = ident.partition(":")
        if kind not in ("doi", "pmid"):
            print(f"  ! {path.name}:{lineno}: unknown id type {kind!r}")
            continue
        entries.append((key, kind, value))
    return entries


# ------------------------------------------------------------------ fetch
def fetch_doi(doi):
    """CrossRef content negotiation returns BibTeX directly."""
    r = requests.get(
        f"https://doi.org/{doi}",
        headers={"Accept": "application/x-bibtex", "User-Agent": UA},
        timeout=30,
        allow_redirects=True,
    )
    if r.status_code != 200:
        raise RuntimeError(f"HTTP {r.status_code} from doi.org")
    text = r.content.decode("utf-8", errors="replace").strip()
    if not text.startswith("@"):
        raise RuntimeError("response was not BibTeX")
    return text


def bibtex_name(name):
    """
    Put an NCBI author name into a form BibTeX parses correctly.

    esummary returns names as "Surname II", e.g. "Al-Hamed MH". BibTeX reads
    an unpunctuated name as "First Last", so that string becomes given name
    "Al-Hamed", family name "MH", and renders as "A.-H. MH". Rewriting it as
    "Al-Hamed, M. H." fixes the parse.

    This reformats the value NCBI returned; it does not add information. A
    name whose final token is not an initials group is braced instead, so
    that BibTeX treats it as one literal unit rather than guessing wrongly.
    """
    head, _, last = name.rpartition(" ")
    if head and re.fullmatch(r"[A-Z]{1,3}", last):
        return f"{head}, " + " ".join(f"{c}." for c in last)
    return f"{{{name}}}"


def fetch_pmid(pmid):
    """
    NCBI esummary gives structured fields; assemble an @article from them.
    Every field written here comes from the response, never from memory.
    """
    r = requests.get(
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi",
        params={"db": "pubmed", "id": pmid, "retmode": "json"},
        headers={"User-Agent": UA},
        timeout=30,
    )
    if r.status_code != 200:
        raise RuntimeError(f"HTTP {r.status_code} from NCBI")
    rec = r.json()["result"][str(pmid)]
    if "error" in rec:
        raise RuntimeError(f"NCBI: {rec['error']}")

    authors = " and ".join(bibtex_name(a["name"])
                           for a in rec.get("authors", [])
                           if a.get("authtype") == "Author")
    year = (rec.get("pubdate", "") or "").split()[0]
    doi = next((i["value"] for i in rec.get("articleids", [])
                if i.get("idtype") == "doi"), "")

    fields = {
        "author": authors,
        "title": rec.get("title", "").rstrip("."),
        "journal": rec.get("fulljournalname") or rec.get("source", ""),
        "year": year,
        "volume": rec.get("volume", ""),
        "number": rec.get("issue", ""),
        "pages": expand_pages(rec.get("pages", "")),
        "doi": doi,
        "pmid": str(pmid),
    }
    body = ",\n".join(f"  {k} = {{{v}}}" for k, v in fields.items() if v)
    return "@article{PLACEHOLDER,\n" + body + "\n}"


def expand_pages(pages):
    """
    NCBI abbreviates the closing page against the opening one: "60-5" means
    60-65 and "27719-27" means 27719-27727. biblatex has no way to know that,
    so it printed "pp. 60-5" in the reference list, which is wrong rather than
    merely terse. Expand the closing number by borrowing the missing leading
    digits from the opening one, and leave anything unrecognised alone.
    """
    m = re.fullmatch(r"\s*(\d+)\s*-\s*(\d+)\s*", pages or "")
    if not m:
        return (pages or "").strip()
    start, end = m.group(1), m.group(2)
    if len(end) < len(start):
        end = start[: len(start) - len(end)] + end
        # a borrowed prefix that runs backwards means the guess was wrong
        if int(end) < int(start):
            return f"{start}--{m.group(2)}"
    return f"{start}--{end}"


def rekey(bibtex, citekey):
    """Replace whatever key the registrar chose with ours."""
    return re.sub(r"^(@\w+\s*\{)[^,]*,", rf"\g<1>{citekey},", bibtex, count=1)


# CrossRef returns real Unicode. pdflatex chokes on Greek letters in
# titles (e.g. the alpha6beta6 structure paper), so map them to math mode.
SANITIZE = {
    "\u03b1": r"$\alpha$", "\u03b2": r"$\beta$", "\u03b3": r"$\gamma$",
    "\u03ba": r"$\kappa$", "\u03bc": r"$\mu$", "\u2013": "--",
    "\u2014": "---", "\u2019": "'", "\u2018": "'",
    "\u201c": "``", "\u201d": "''", "\u00a0": " ",
    # Publishers hyphenate double-barrelled surnames with U+2010 HYPHEN
    # rather than ASCII '-' (Dionisi-Vici, Scholl-Buergi in the MMA/PA
    # guidelines). Recent TeX Live accepts it; older ones do not, and it
    # is a plain hyphen either way.
    "\u2010": "-", "\u2011": "-", "\u00ad": "",
    "\u2212": "$-$", "\u2032": "$'$",
}

# Letters in the Latin-1 supplement (accented names) are handled by
# inputenc + T1 and are left alone. Anything else surviving sanitise is
# a codepoint we have no rule for, and is reported rather than shipped
# silently into a build that may or may not accept it.
def unmapped(text):
    return sorted({c for c in text
                   if ord(c) > 127 and not (0x00C0 <= ord(c) <= 0x00FF)})


# BibTeX's month macros are three letters. CrossRef emits some months in
# full ("month=July"), which is an undefined macro: biber warns and drops
# the field. Map the value onto the macro standing for the same month.
# This re-spells registrar data, it does not supply any.
MONTHS = ("jan", "feb", "mar", "apr", "may", "jun",
          "jul", "aug", "sep", "oct", "nov", "dec")


def normalise_months(bibtex):
    def sub(m):
        macro = m.group(2).strip().lower()[:3]
        return m.group(1) + macro if macro in MONTHS else m.group(0)
    # Only unbraced values, which is the form the macro warning comes from.
    return re.sub(r"(month\s*=\s*)([A-Za-z]+)", sub, bibtex)


def sanitize(bibtex):
    for bad, good in SANITIZE.items():
        bibtex = bibtex.replace(bad, good)
    return normalise_months(bibtex)


# ------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true",
                    help="exit non-zero if any source is unresolved")
    ap.add_argument("--force", action="store_true", help="bypass the cache")
    args = ap.parse_args()

    if not SOURCES.exists():
        sys.exit(f"missing {SOURCES}")

    cache = {} if args.force or not CACHE.exists() else json.loads(CACHE.read_text())
    entries = read_sources(SOURCES)
    print(f"{len(entries)} sources listed\n")

    out, todo = [], []
    for key, kind, value in entries:
        cache_id = f"{kind}:{value}"
        if cache_id in cache:
            out.append(rekey(cache[cache_id], key))
            print(f"  = {key:24s} (cached)")
            continue
        try:
            bib = fetch_doi(value) if kind == "doi" else fetch_pmid(value)
            cache[cache_id] = bib
            out.append(rekey(bib, key))
            print(f"  + {key:24s} {cache_id}")
        except Exception as exc:                       # noqa: BLE001
            todo.append(f"{key}\t{cache_id}\t{exc}")
            print(f"  ! {key:24s} {cache_id}  FAILED: {exc}")
        time.sleep(PAUSE)

    def citekeys(s):
        return set(re.findall(r"^@\w+\{([^,]+),", s, re.M))

    header = (
        "% Generated by tools/fetch_refs.py. DO NOT EDIT.\n"
        "% Edit paper/refs.sources and re-run.\n"
        "% Entries fetched from CrossRef (doi.org) and NCBI PubMed.\n\n"
    )
    manual = ""
    if MANUAL.exists():
        body = MANUAL.read_text().strip()
        #  A PROVISIONAL entry is one whose key ALSO has a live line in
        #  refs.sources. It exists so that a source added from an environment
        #  that cannot reach a registrar still produces a document that builds
        #  and a claim register that checks; the moment a run resolves that
        #  key for real, the fetched entry is authoritative and the
        #  hand-written one must go, or the two collide on the same key.
        #  Dropping it here rather than asking a human to remember is the
        #  whole point: the version that survives is the fetched one.
        live = {k for k, _, _ in entries}
        dropped_prov = []
        for key in sorted(live):
            m = re.search(r"^@\w+\{" + re.escape(key) + r",.*?^\}\n?",
                          body, re.M | re.S)
            if m and key in citekeys("\n\n".join(out)):
                body = body[:m.start()] + body[m.end():]
                dropped_prov.append(key)
        if dropped_prov:
            print(f"\n  provisional entries superseded by a real fetch and "
                  f"removed from the output: {', '.join(dropped_prov)}")
            print(f"    delete them from {MANUAL.relative_to(ROOT)} as well")
        manual = (MANUAL_MARK + "\n" + body.strip() + "\n")

    text = sanitize(header + "\n\n".join(out) + "\n" + manual)

    # refs.bib is tracked so the paper builds without network access. An
    # unconditional overwrite defeats exactly that: one unreachable registrar
    # and every entry it could not re-fetch vanishes from the bibliography,
    # and the next build quietly produces a paper full of [?] citations. So
    # never write a refs.bib carrying fewer citekeys than the one on disk,
    # and never write at all under --check, which is a check and not a build.
    existing = OUT.read_text() if OUT.exists() else ""
    dropped = citekeys(existing) - citekeys(text)

    if args.check:
        wrote = False
        print(f"\n--check: {OUT.relative_to(ROOT)} left untouched")
    elif dropped:
        wrote = False
        print(f"\n  ! NOT writing {OUT.relative_to(ROOT)}: this run resolved "
              f"{len(citekeys(text))} entries, the file on disk has "
              f"{len(citekeys(existing))}, and writing would lose "
              f"{', '.join(sorted(dropped))}.")
        print("    The committed bibliography is kept, so the build still "
              "succeeds. Re-run with network access to the registrars.")
        # The hand-written half of refs.bib is not fetched from anywhere, so
        # refusing the whole write also froze it, and an entry added to
        # refs.manual.bib could never reach the bibliography from an
        # environment the registrars are unreachable from -- which is the only
        # kind of environment that file exists for. Splice that half in on its
        # own. The fetched half is copied across untouched, and the splice is
        # abandoned rather than allowed to drop a key.
        head = existing.split(MANUAL_MARK)[0].rstrip("\n") + "\n"
        spliced = sanitize(head + manual) if manual else existing
        if spliced != existing and not citekeys(existing) - citekeys(spliced):
            OUT.write_text(spliced)
            print(f"    refreshed the hand-written half from "
                  f"{MANUAL.relative_to(ROOT)}; the fetched half is unchanged.")
    else:
        OUT.write_text(text)
        wrote = True

    CACHE.write_text(json.dumps(cache, indent=1))

    leftover = unmapped(text)
    if leftover:
        print("\n  ! unmapped non-ASCII in refs.bib, add a SANITIZE rule: "
              + " ".join(f"{c!r} (U+{ord(c):04X})" for c in leftover))

    if todo:
        TODO.write_text("citekey\tidentifier\treason\n" + "\n".join(todo) + "\n")
        print(f"\n{len(todo)} unresolved -> {TODO.relative_to(ROOT)}")
    else:
        TODO.unlink(missing_ok=True)

    if wrote:
        print(f"\nwrote {OUT.relative_to(ROOT)}  ({len(out)} entries"
              f"{', plus refs.manual.bib' if manual else ''})")

    if args.check and todo:
        sys.exit(1)


if __name__ == "__main__":
    main()
