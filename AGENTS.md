# AGENTS.md

## Project

This repository is an R package: powRICLPM.

## Branch rules

- main is a clean copy of upstream. Do not edit main directly.
- playground is for experiments, exploratory changes, and Codex-assisted trials.
- integration is for clean PR-ready changes.

## Working rules

- Keep changes small and reviewable.
- Do not make broad refactors unless explicitly asked.
- Do not add new package dependencies without explaining why.
- Do not change public behavior without updating tests and documentation.
- Prefer preserving the existing package style.
- Before editing, explain the intended change when the task is ambiguous.

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