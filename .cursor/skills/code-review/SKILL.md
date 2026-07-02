---
name: code-review
description: Review code changes for correctness, security, performance, maintainability, and ML/DS pitfalls. Use when the user asks for a code review, PR review, or feedback on diffs/patches.
---

# Code Review (Cursor)

## What to do

1. Identify the change scope (files, entrypoints, behavior change).
2. Review in this order:
   - Correctness & edge cases
   - Security & secrets
   - Performance (time/memory/I/O)
   - Maintainability (API, naming, structure)
   - Tests & reproducibility
3. Give feedback grouped by severity:
   - **Must fix**: likely bug, security, data loss, incorrect results
   - **Should fix**: maintainability/perf issues, flaky behavior
   - **Nice to have**: style, minor refactors

## Checklists

### Python/ML specific
- Validate shapes/dtypes, device placement, random seeds, train/eval mode.
- Avoid data leakage; verify split logic and preprocessing fit/transform order.
- Watch for silent broadcasting, float precision, NaNs/Infs, overflow.
- Ensure metrics are computed on correct targets and averaged correctly.

### Security & hygiene
- Don’t commit secrets (`.env`, tokens, keys, credentials).
- Prefer parameterized queries; avoid shell injection patterns.

### Tests
- If behavior changes, propose minimal tests (unit or smoke).
- For scripts, add a quick “dry-run” path or tiny fixture input when feasible.

## Output format

Provide:
- **Summary** (1–3 bullets)
- **Must fix**
- **Should fix**
- **Nice to have**
- **Test plan** (explicit commands / steps)
