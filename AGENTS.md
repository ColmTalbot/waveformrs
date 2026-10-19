# Notes for agents working on this repository

- **Keep README.md up to date.** When you add a feature or a model, or remove
  one, update the lists in `README.md` in the same change. The README says
  what is implemented, in plain language: no technical details and no timing
  or accuracy numbers. Those belong in `CHANGELOG.md` and the report
  (`report/`).
- **Record ported code in NOTICE.md.** The crate is GPL-3.0-or-later because
  much of it is translated from LALSuite (GPL-2.0-or-later) and GSL
  (GPL-3.0-or-later). Before porting code from anywhere else, check that its
  license is compatible with GPL-3.0-or-later. When you port or follow code
  from another file or project, add it to `NOTICE.md` with its copyright
  holders and license. New source files start with the same three-line
  SPDX/copyright header as the existing ones.
