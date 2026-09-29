# Conventions

- Project dir: `/u/project/cluo/terencew/claude/project_ideas/latent_genos/`. This repo (`github/demux_latent_genos`) is a hand-maintained mirror, and the `md/` docs live only here.
- Simulated pools, and the scripts that run on them, stay in `project_ideas/pool_design/`:
  `scripts/latent_genos/qsub/A03*` → copied here under `scripts/pool_design/qsub/`.
  Outputs go next to the demuxlet outputs: `<tree>/<pool>/demux/souporcell/<mod>/<sample>/`.
- Task list rows: `<tree> <pool> <k> <gex|atac>` (`pool_design/txt/latent_genos_souporcell_tasks.txt`).
- Envs: `souporcell` (conda), `monopogen` (conda, exists), `allcools` for Python analysis. `claude_env` is Node.js.
- Script naming: pool_design already uses A00–A02 for the simulation → demuxlet steps, so this project starts at A03 (souporcell); imputation will be A04 and Monopogen A05.
- Array jobs: `# ID=1` / `ID=$SGE_TASK_ID`, `-tc`, skip-if-exists, `time`, and no email flags.
- Git: `module load git/2.40.0`. Run `git fetch` before editing docs, because the user sometimes edits on GitHub.
