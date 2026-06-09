# AGENTS.md

## Project

This repository is an R package: powRICLPM.

## Branch rules

- `main` is a clean copy of upstream. Do not edit `main` directly.
- `playground` is for experiments, exploratory changes, and Codex-assisted trials.
- `integration` is for clean PR-ready changes.
- The usual flow is `playground` -> `integration` -> final PR-ready branch/change set.
- `upstream` is the package owner's repository and should be treated as read-only. The push URL may intentionally be disabled.
- `origin` is the user's fork on GitHub. Push `playground` changes to `origin/playground` only when the user asks.

## Working rules

- Keep changes small and reviewable.
- Do not make broad refactors unless explicitly asked.
- Do not add new package dependencies without explaining why.
- Do not change public behavior without updating tests and documentation.
- Prefer preserving the existing package style.
- Before editing, explain the intended change when the task is ambiguous.
- Do not commit or push unless the user asks for it.
- When committing, keep the commit scoped to the requested work and leave unrelated local changes alone.

## Project memory and Obsidian logging

The project memory lives in the user's existing Obsidian vault, not in a separate vault:

`C:\Users\Admin\Documents\Obsidian Vault\UU\powRICLPM`

At the start of a new Codex chat, read this file first:

`C:\Users\Admin\Documents\Obsidian Vault\UU\powRICLPM\powRICLPM.md`

That note is the top-level runner. It contains compact timelines for ideas, decisions, work logs, and verifications. Use it to understand what has recently happened and which deeper notes are relevant before editing code.

### Obsidian folders

- `Decisions\`: design and API decisions. These notes explain what changed and why.
- `Work logs\`: implementation records. These notes summarize what changed in each touched file, why it changed, and how it was done.
- `Verifications\`: checks that were run. These notes record what was checked, how it was checked, and the result.
- `Ideas\`: user-owned thinking space. Do not create or edit notes in `Ideas\` unless the user explicitly asks.

### Obsidian naming conventions

- Shared Codex-created notes use `YYYY-MM-DD-N title.md`.
- `N` is the same-day sequence number. Inspect the runner and folders before choosing the next number.
- Related decision, work-log, and verification notes for the same change should use the same date and number.
- Do not use timestamps in note titles.
- Keep runner-table descriptions short, ideally one sentence and one line.

Examples:

- `Decisions\2026-06-05-1 clarify renamed argument conflict errors.md`
- `Work logs\2026-06-05-1 renamed argument conflict message implementation.md`
- `Verifications\2026-06-05-1 renamed argument conflict message verification.md`

### Obsidian note content

Decision notes should usually include:

- `Short description`
- `What changed`
- `Rationale`
- A sentence link to the related work log and verification note where useful.

Work logs should usually include:

- `Summary`
- `Files changed`
- `Implementation notes`
- File-by-file prose that is detailed enough to understand the change later without reading the diff.

Verification notes should usually include:

- `What was checked`
- `How it was checked`
- `Result`
- Commands or relevant snippets when they help future interpretation.

### Updating the runner

After creating shared notes, update:

`C:\Users\Admin\Documents\Obsidian Vault\UU\powRICLPM\powRICLPM.md`

Add the new notes to the appropriate timeline table. Keep the runner readable and compact. Link notes with Obsidian wiki links like `[[2026-06-05-1 clarify renamed argument conflict errors]]`.

## R package workflow

When changing exported functions:
- update roxygen documentation if needed
- update tests if behavior changes
- update examples if relevant

Useful checks:
- devtools::document()
- devtools::test()
- devtools::check()

If a check cannot be run, explain why.

## Review guidelines

When reviewing changes, check for:
- broken tests
- missing documentation
- accidental API changes
- unnecessary dependencies
- R CMD check issues
