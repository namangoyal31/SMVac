"""Check relative links and file references in the active Markdown documentation.

Scans README.md, CHANGELOG.md, CONTRIBUTING.md, docs/*.md, figures/README.md, legacy/README.md and audit/README.md.
Reports every Markdown link target, image and back-quoted repository path (e.g. `validation/results/x.csv`) that does
not exist. Exit code 1 if anything is missing.
Usage: python validation/check_links.py
"""
from __future__ import annotations

import fnmatch
import glob
import os
import re
import subprocess
import sys

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
FILES = ["README.md", "CHANGELOG.md", "CONTRIBUTING.md", "figures/README.md", "legacy/README.md", "audit/README.md"]
FILES += sorted(os.path.relpath(p, ROOT) for p in glob.glob(os.path.join(ROOT, "docs", "*.md")))
LINK = re.compile(r"!?\[[^\]]*\]\(([^)\s]+)\)")
CODE = re.compile(r"`([A-Za-z0-9_./\-*]+\.(?:md|py|cpp|hpp|csv|json|toml|txt|png|cff|yml|gz))`")
# The review documents of the original implementation describe the audited revision a7c070a. A path they cite is
# accepted if it exists now, if it was moved to legacy/ by the clean-up, or if it existed in that revision (files
# removed later remain in the git history).
AUDITED_REVISION = "a7c070a"
RELOCATED = {"apps/": "legacy/apps/", "figures/": "legacy/figures/", "data/": "legacy/data/",
             "scripts/": "legacy/scripts/", "SMVacuumDecay/": "legacy/SMVacuumDecay/"}
HISTORICAL = {"REPOSITORY_AUDIT.md", "SCIENTIFIC_METHODOLOGY.md", "NUMERICAL_METHODS_REVIEW.md", "STAGE_C_PLAN.md",
              "README.md"}   # README.md here means legacy/README.md and audit/README.md (paths of the old layout)
PREFIXES = ("docs/", "validation/", "results/", "figures/", "python/", "src/", "include/", "apps/", "tests/",
            "configs/", "legacy/", "audit/", ".github/")


def audited_paths():
    """Paths tracked in the audited revision (empty if git or the revision is unavailable)."""
    try:
        out = subprocess.run(["git", "ls-tree", "-r", "--name-only", AUDITED_REVISION], cwd=ROOT,
                             capture_output=True, text=True, check=True).stdout
        return set(out.split())
    except (OSError, subprocess.CalledProcessError):
        return set()


AUDITED = audited_paths()


def in_audited_revision(target):
    return any(fnmatch.fnmatch(p, target) for p in AUDITED)


def exists(base, target):
    target = target.split("#")[0]
    if not target or "://" in target or target.startswith("mailto:"):
        return True
    for cand in (os.path.join(base, target), os.path.join(ROOT, target)):
        if glob.glob(cand) or glob.glob(cand + ".gz") or os.path.exists(cand):
            return True
    return False


def main():
    missing = []
    for f in FILES:
        p = os.path.join(ROOT, f)
        if not os.path.exists(p):
            missing.append((f, "(file itself)"))
            continue
        base = os.path.dirname(p)
        text = open(p, encoding="utf-8").read()
        for t in LINK.findall(text):
            if not exists(base, t):
                missing.append((f, t))
        for t in CODE.findall(text):
            if not t.startswith(PREFIXES) or exists(base, t):
                continue
            if os.path.basename(f) in HISTORICAL:
                moved = [t.replace(a, b, 1) for a, b in RELOCATED.items() if t.startswith(a)]
                if any(exists(ROOT, m) for m in moved) or in_audited_revision(t):
                    continue
            missing.append((f, t))
    for f, t in missing:
        print(f"MISSING  {f}: {t}")
    print(f"checked {len(FILES)} files; {len(missing)} missing references")
    sys.exit(1 if missing else 0)


if __name__ == "__main__":
    main()
