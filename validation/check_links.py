"""Check relative links and file references in the active Markdown documentation.

Scans README.md, CHANGELOG.md, CONTRIBUTING.md, docs/*.md, figures/README.md, legacy/README.md and audit/README.md
(docs/archive/ is a frozen historical record and is not checked). Reports every Markdown link target, image and
back-quoted repository path (e.g. `validation/results/x.csv`) that does not exist. Exit code 1 if anything is missing.
Usage: python validation/check_links.py
"""
from __future__ import annotations

import glob
import os
import re
import sys

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
FILES = ["README.md", "CHANGELOG.md", "CONTRIBUTING.md", "figures/README.md", "legacy/README.md", "audit/README.md"]
FILES += sorted(os.path.relpath(p, ROOT) for p in glob.glob(os.path.join(ROOT, "docs", "*.md")))
LINK = re.compile(r"!?\[[^\]]*\]\(([^)\s]+)\)")
CODE = re.compile(r"`([A-Za-z0-9_./\-*]+\.(?:md|py|cpp|hpp|csv|json|toml|txt|png|cff|yml|gz))`")
# The Stage A/B review documents describe the audited revision a7c070a; files moved by the clean-up are resolved
# through this map, and files deliberately removed (kept in git history) are listed as such, not as missing.
RELOCATED = {"apps/": "legacy/apps/", "figures/": "legacy/figures/", "docs/": "docs/archive/", "data/": "legacy/data/",
             "scripts/": "legacy/scripts/", "SMVacuumDecay/": "legacy/SMVacuumDecay/"}
HISTORICAL = {"REPOSITORY_AUDIT.md", "SCIENTIFIC_METHODOLOGY.md", "NUMERICAL_METHODS_REVIEW.md", "STAGE_C_PLAN.md",
              "README.md"}   # README.md here means legacy/README.md and audit/README.md (paths of the old layout)
REMOVED = {"docs/AGENTS.md", "docs/prompt_for_claude.md", "docs/critique_prompt_for_claude.md",
           "docs/reformulate_prompt_for_claude.md", "docs/README"}
PREFIXES = ("docs/", "validation/", "results/", "figures/", "python/", "src/", "include/", "apps/", "tests/",
            "configs/", "legacy/", "audit/", ".github/")


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
            if not t.startswith(PREFIXES) or exists(base, t) or t in REMOVED:
                continue
            moved = [t.replace(a, b, 1) for a, b in RELOCATED.items() if t.startswith(a)] \
                if os.path.basename(f) in HISTORICAL else []
            if not any(exists(ROOT, m) for m in moved):
                missing.append((f, t))
    for f, t in missing:
        print(f"MISSING  {f}: {t}")
    print(f"checked {len(FILES)} files; {len(missing)} missing references")
    sys.exit(1 if missing else 0)


if __name__ == "__main__":
    main()
