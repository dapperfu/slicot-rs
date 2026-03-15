# Git workflow (from .cursor/rules/git/)

Follow these steps for every commit. See `.cursor/rules/git/*.mdc` for full rules.

## 1. Before committing

**Upstream sync** (`upstream-sync.mdc`):

```bash
git fetch upstream && git fetch origin
branch=$(git symbolic-ref refs/remotes/upstream/HEAD 2>/dev/null | sed 's@^refs/remotes/upstream/@@') || branch="main"
git merge --no-edit "upstream/$branch"
```

**User config** (`user-config.mdc`):

```bash
git config user.name "$(whoami) | Cursor.sh | Auto"
git config user.email "$(whoami)@$(hostname).local"
```

## 2. Commit message format (`commit-format.mdc`)

- Line 1: Brief one-liner summary
- Blank line
- List of changes (each line prefixed with `- `)
- Blank line
- Separator: exactly 5 dashes `-----`
- Technical attribution:
  - Prompt: {{prompt}}
  - Context: {{brief_description}}
  - Justification: {{justification}}
  - Technical details: Model, IDE, Generation method, Code style, Dependencies, Tokens

## 3. Atomicity (`commit-atomicity.mdc`)

- Prefer **one commit per file**.
- Exceptions: header+impl, test+fix, tightly coupled files, requirements+code, TODO checkoff, atomic feature.

## 4. After commit (`push-requirement.mdc`)

```bash
git push origin main   # if remote exists
```

## Reference

- `commit-requirement.mdc` — when to commit (after every file / after each prompt with changes)
- `commit-format.mdc` — full message format
- `commit-atomicity.mdc` — single-file vs grouped commits
- `push-requirement.mdc` — push after commit
- `upstream-sync.mdc` — sync before commit
- `user-config.mdc` — user.name / user.email for AI commits
