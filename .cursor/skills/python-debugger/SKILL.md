---
name: python-debugger
description: Systematically debug Python errors using reproduction, minimal tests, logging, and environment checks. Use when errors occur (tracebacks), unexpected outputs, crashes, or performance regressions in Python/ML code.
---

# Python Debugger

## Core loop

1. Reproduce reliably:
   - capture command, input files, seed, and environment
2. Reduce:
   - isolate to smallest script/function that still fails
3. Diagnose:
   - read full traceback; identify first user-code frame
   - validate assumptions (types, shapes, ranges, paths)
4. Fix:
   - make the smallest change that addresses the root cause
5. Verify:
   - rerun reproduction and a small regression check

## Fast checks (ML-friendly)
- Print/assert: shapes, dtypes, device (`cpu/cuda`), `model.training`.
- Seed: `random`, `numpy`, `torch` (and cudnn determinism if needed).
- Data: confirm normalization ranges and label alignment.
- Performance: watch accidental `.cpu().numpy()` in hot loops, and dataloader workers.

## Output format
- **Root cause** (one sentence)
- **Fix** (what to change)
- **Proof** (how to verify)
