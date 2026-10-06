"""PIE (Arc, arXiv/bioRxiv 2026.10.02.756297) zero-shot inference on Modal.

The replogle_xdataset checkpoint trains on Tahoe, Jiang, ARC VCC 25 and X-Atlas/Orion only, so Jurkat (Nadig)
responses are unseen: a fair test against our Jurkat local eval.
  modal volume put vcc-pie ~/vcc/pie_work/vccjurkat /vccjurkat
  modal volume put vcc-pie ~/vcc/pie_work/query_jurkat.json /query_jurkat.json
  modal run modal_pie.py::infer --query query_jurkat.json --dirs vccjurkat --out jurkat
  modal volume get vcc-pie /out/jurkat.parquet ~/vcc/pie_work/
"""
import modal

app = modal.App("vcc-pie")
vol = modal.Volume.from_name("vcc-pie", create_if_missing=True)
img = modal.Image.debian_slim(python_version="3.12").pip_install("arc-pie==1.0.0", "huggingface_hub")
ENV = {"PIE_DATA_ROOT": "/vol/data", "PIE_RUNS_ROOT": "/vol/runs", "PIE_CACHE_DIR": "/vol/cache",
       "HF_HOME": "/vol/hf", "WANDB_ENTITY": "none", "WANDB_PROJECT": "none", "WANDB_MODE": "disabled"}


@app.function(image=img, volumes={"/vol": vol}, gpu="L40S", cpu=8, memory=98304, timeout=4 * 3600,
              env=ENV)
def infer(query: str, dirs: str, out: str, experiment: str = "replogle_xdataset"):
    import os
    import subprocess

    from huggingface_hub import hf_hub_download
    run = f"/vol/runs/{experiment}"
    os.makedirs(run, exist_ok=True)
    for f in ("best_auprc.ckpt", "config.yaml", "data_stats.json"):
        if not os.path.exists(f"{run}/{f}"):
            hf_hub_download(f"arcinstitute/PIE_{experiment}", f, local_dir=run)
    vol.commit()
    os.makedirs("/vol/out", exist_ok=True)
    pre = "[" + ",".join(f"/vol/{d}" for d in dirs.split(",")) + "]"
    subprocess.run(["pie", "infer", f"experiment_name={experiment}", "rows_kind=query",
                    f"rows_path=/vol/{query}", f"preprocessed_dirs={pre}", "device=cuda",
                    f"output_path=/vol/out/{out}.parquet", "overwrite=true"], check=True)
    vol.commit()
