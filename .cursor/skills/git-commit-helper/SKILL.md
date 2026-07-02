---
name: git-commit-helper
description: Generate high-quality git commit messages by analyzing git diff and intent. Use when the user asks to write commit messages, prepare commits, or summarize changes for version control.
---

# Git Commit Helper

## Workflow

1. Inspect repository state:
   - `git status`
   - `git diff` (unstaged)
   - `git diff --staged` (staged)
   - `git log -5 --oneline` (style)
2. Identify:
   - What changed? (surface area)
   - Why changed? (intent)
   - Risk? (breaking change, migration, data format)
3. Draft message following repo style; default to Conventional Commits if unclear:
   - `feat:`, `fix:`, `refactor:`, `docs:`, `test:`, `chore:`

## Message template (default)

Title (<= 72 chars): `<type>(<scope>): <intent>`

Body:
- What/why, not a line-by-line summary
- Mention breaking changes explicitly
- Mention key trade-offs or follow-ups

## Quality bar
- Avoid “update”/“changes” titles.
- Prefer intent verbs: “prevent”, “support”, “reduce”, “ensure”, “correct”.
- If multiple unrelated changes, recommend splitting commits.
