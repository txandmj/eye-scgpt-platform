# Mock Pipeline Mode

Enable a fast, fake pipeline for development and UI/DB testing. All API routes and responses remain unchanged. The server will produce tiny demo artifacts quickly instead of running the real model.

## Enable

Set the environment variable before starting the backend:

```bash
# macOS/Linux
export PIPELINE_MODE=mock
python start.py
```

Or put it into your shell/PM2/systemd environment. Omit or set to `real` to run the actual pipeline.

## Behavior in mock mode
- Upload: same as real; creates an `inferences` row with `pending` and a `raw_h5ad` artifact.
- Start processing: status becomes `processing` briefly, then within ~1–2 seconds:
  - Creates `results/{job_id}/predictions.csv` with tiny demo content.
  - Creates `results/{job_id}/umap/umap_0.png` and a small `{base}_enriched.h5ad`.
  - Inserts corresponding `artifacts` rows.
  - Marks the inference `completed`.

Front-end behavior and all API responses are identical to real mode.

## Notes
- Files are stored with relative paths in the DB and resolved at read time.
- Logs will include a clear `[mock]` prefix when mock mode is active.



