REMOTE_BASE="htcf:/scratch/mblab/chasem/latent_factor_eda"

for d in ./slurm_bundles; do
    d="${d%/}"
    name="$(basename "$d")"
    rclone sync "$d" "$REMOTE_BASE/$name" --progress
done
