# Documentation has moved

The project documentation is the Quarto site in [`docs/`](../docs/) (`quarto preview docs`).

A few error messages in pipeline code still point here; where to look instead:

| Message mentions | See |
|---|---|
| "Adding a new method" | `docs/methods.qmd`, "Adding a method" |
| the config reference / "Preprocessing script contract" | `docs/data.qmd` |
| the "sPCA" section | `docs/methods.qmd`, "sPCA: choosing the sparsity" |

(Some of those messages are inside code the pipeline's targets depend on, and rewording them would
rerun targets -- in sPCA's case every sPCA fit. See `TODOS.md` for what each would rerun.)
