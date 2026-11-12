# -*- coding: utf-8 -*-
"""349b — negative-control for the docs/148 §2 remediation gates.

Core lesson (2026-11-12): a green gate does NOT prove it covers the edit. This
script injects a targeted regression into ONE carrier (docs/08, the true source),
asserts the corresponding gate FAILS, then byte-exactly restores the original and
asserts the gate returns to green. No deliverable is left modified; safe & idempotent.

Cases
-----
1. 310 section (G) "Table S33 索引行存在" — scope is md+mirror only (docx strips
   the supplementary index). Tamper the S33 index line in docs/08 -> 310 must FAIL.
2. 182 P4 "measurement artefact + cardiovascular risk" — revert the strontium
   clause to the drifted "as an artefact" form -> 182 must FAIL.

Run: python3 code/349b_negative_control.py   (never writes a final deliverable)
"""
import subprocess, hashlib, re, sys

F = "docs/08_投稿稿_OA_BMD_fracture.md"
PY = sys.executable


def run(script):
    r = subprocess.run([PY, script], capture_output=True, text=True)
    return r.stdout + r.stderr


def sha(b):
    return hashlib.sha256(b).hexdigest()


def green(out):
    return "FAIL: 0" in out or "PASSED" in out or "SAFE" in out


G310 = "code/310_assert_no_regression.py"
G182 = "code/182_final_qc.py"

CASES = [
    # label, gate, inject_old, inject_new, needle expected in FAIL output
    ("310 §(G) S33 index scope (md/mirror-only)",
     G310,
     "**Table S33.** Section-2 review remediation audit",
     "**Table S33.** (negative-control tamper)",
     "Table S33 索引行存在"),
    ("182 P4 strontium 'measurement artefact'",
     G182,
     "as a measurement artefact and is restricted for cardiovascular risk",
     "as an artefact and is restricted for cardiovascular risk",
     "P4"),
]

orig = open(F, "rb").read()
h0 = sha(orig)

# baseline: every gate must be green before injection
for _, script, _, _, _ in CASES:
    base = run(script)
    if not green(base):
        print(f"[ABORT] baseline not green for {script}\n{base[-400:]}")
        sys.exit(2)

results = []
try:
    for label, script, old, new, needle in CASES:
        ob, nb = old.encode(), new.encode()
        assert orig.count(ob) == 1, f"inject anchor count != 1 for {label} (got {orig.count(ob)})"
        open(F, "wb").write(orig.replace(ob, nb, 1))
        out = run(script)
        m = re.search(r"FAIL:\s*(\d+)", out)
        nfail = int(m.group(1)) if m else (-1 if "[FAIL]" in out else 0)
        hit = needle in out
        results.append((label, nfail != 0, nfail, hit, needle))
        open(F, "wb").write(orig)  # restore before next case
finally:
    open(F, "wb").write(orig)  # guarantee restore even on exception

h1 = sha(open(F, "rb").read())
assert h0 == h1, f"[FATAL] SHA not restored: {h0} -> {h1}"

# post-restore: gates must be green again
restored = {}
for _, script, *_ in CASES:
    out = run(script)
    restored[script] = green(out)
    assert green(out), f"[FATAL] gate not restored green: {script}\n{out[-400:]}"

print("=" * 74)
print("349b · §2 remediation gate negative-control")
print("=" * 74)
for label, ok, nfail, hit, needle in results:
    tag = "PASS" if (ok and hit) else "FAIL"
    print(f"[{tag}] {label}")
    print(f"        injected -> gate FAIL={nfail}; needle '{needle}' in output={hit}")
print(f"[PASS] byte-exact restore: {h0[:16]}… == {h1[:16]}…  ({h0 == h1})")
for script, g in restored.items():
    print(f"[{'PASS' if g else 'FAIL'}] {script} green after restore")
if all(ok and hit for _, ok, _, hit, _ in results) and h0 == h1 and all(restored.values()):
    print("VERDICT: NEGATIVE CONTROL VALID — both §2 assertions are non-no-op")
else:
    print("VERDICT: NEGATIVE CONTROL FAILED")
    sys.exit(1)
