# Legion GPU feasibility and offline preflight proposal

Prepared for review on 15 October 2026. **The script has not been executed
on the ThinkPad or the Legion.** No CUDA/PyTorch/Boltz installation, model
download, inference or GPU hosting has been performed for this proposal.

The target is the user's **Legion Pro 7 Gen 10, RTX 5070 Ti Laptop GPU,
12 GB GDDR7 VRAM, 64 GB system RAM, Windows 11**, currently powered off
until its charger arrives. These are user-provided specifications.
ThinkPad hardware observations are excluded from the Legion assessment.
The Legion's driver, available VRAM and installed software remain unknown.

The separate Phase 1 CI gate is complete: PR #1 is at
`cba754944b8ce80028fb9e0c2602fcbc7c35ff0c`, and both Windows 2022 and
Ubuntu 24.04 passed all 45 tests, including the eight CSV/summary byte
comparisons and subsequent regeneration/drift checks.
[CI run](https://github.com/pkaysantana/BRD4-structure-guided-discovery/actions/runs/38784555615).
Phase 1 is ready for user merge approval; the Legion is not a dependency
of that decision. No merge or Phase 2 push was performed.
Three historical review threads from `f1a9966` remain open in GitHub.
Their report/CSV issues are addressed by `5c06db2` and `cba7549`; two
threads are marked outdated, and native CI now confirms the CSV fix.
Discussion state was not changed and no review comment was posted.

## Revised feasibility

A 12 GB GPU is a plausible candidate for the small 4MEN structure-only
pilot: 127 protein residues, 20 ligand heavy atoms, one seed, one diffusion
sample and one parallel sample. **A 24 GB remote GPU is not presumed
necessary.** Actual peak memory with the full MSA has not been measured.
64 GB host RAM helps input/checkpoint loading but does not replace VRAM.
The reported power allowance of up to 140 W affects throughput and thermal
behaviour, not the 12 GB capacity; this script changes no clocks or power limits.

| Option | Assessment |
|---|---|
| Local 12 GB Legion | Reasonable first candidate after driver/runtime validation. Limited VRAM headroom and Windows/WSL setup require checking; no claim yet that Boltz fits. |
| Remote 24 GB Linux GPU | More VRAM headroom and potentially simpler runtime setup. Requires its own compatible software, timing, cost approval and separate run directory. A fallback option rather than an automatic requirement. |

The existing config's `recommended_vram_GB: 24` is a provisioning
recommendation, not a hard VRAM gate in the implementation. The revised
assessment supports evaluating the 12 GB local option. The scientific
config and primary MSA policy have not been rewritten for this proposal.

The primary condition remains a **real cached full MSA**, with the existing
8192-sequence limit. The no-MSA condition remains a smoke test only.
If a full-MSA run later exhausts memory, stop and record the failure.
Either retain that condition on a larger GPU or explicitly review a
separate MSA-depth ablation. The current validator deliberately rejects
limits other than 8192: do not bypass it or silently change the command.

## Verified build choice and limits of verification

[NVIDIA's compute-capability table](https://developer.nvidia.com/cuda/gpus)
places the RTX 5070 Ti family at **compute capability 12.0 (`sm_120`)**.
The script will check the actual laptop device rather than assume this
from the product name alone.

The proposed reproducible stack is:

| Component | Proposed pin |
|---|---|
| Python | 64-bit Python 3.12 |
| PyTorch | **`2.7.1+cu128`**, from the official CUDA 12.8 wheel index |
| NumPy | `1.26.4` |
| Boltz, later Linux environment stage | `2.2.1`, without the optional `[cuda]` extras |
| RDKit, when Boltz is installed | `2026.3.6`, matching the registered atom-naming implementation |

[PyTorch 2.7 introduced Blackwell support with CUDA 12.8 wheels](https://pytorch.org/blog/pytorch-2-7/).
The official index was checked for the actual Python 3.12
`torch-2.7.1+cu128` wheels for both `win_amd64` and
`manylinux_2_28_x86_64`. Only index metadata was read; neither wheel was
downloaded. See the [official version instructions](https://pytorch.org/get-started/previous-versions/).

Boltz 2.2.1 declares Torch ≥2.2 and NumPy ≥1.26,<2.0.
Its pinned Lightning 2.5.0 declares Torch ≥2.1.0. The proposed versions
satisfy those metadata constraints. This is not an end-to-end Boltz
runtime validation or a complete resolved environment lock.

The inspected Boltz source uses `bf16-mixed` for Boltz-2.
`--no_kernels` passes `use_kernels=False`; it avoids the optional optimized
kernel path, but does not make an unsupported PyTorch CUDA build compatible.
The preflight therefore tests both FP32 arithmetic and a tiny BF16 matrix
multiplication on CUDA. It never imports Boltz or loads a checkpoint.

Use a current NVIDIA/Lenovo Windows driver supporting the actual laptop.
The CUDA 12.8 GA toolkit baseline is Windows driver 570.65 / Linux 570.26;
these are baseline checks, **not a claim that those early drivers support
every RTX 5070 Ti Laptop configuration**. Device recognition and actual
CUDA operations must also pass. A newer compatible driver is acceptable.
`nvidia-smi`'s CUDA version describes driver capability, whereas
`torch.version.cuda` identifies the wheel runtime.

## Files and isolation

The standalone script is
[`scripts/boltz_gpu_preflight.py`](../../scripts/boltz_gpu_preflight.py).
It can be copied by itself to the Legion; it does not require this
repository or the Phase 1 environment.

Eight CPU-only tests of the parser and decision logic passed using
synthetic observations. Subprocess calls and Torch imports were blocked
in those tests. They did not call the hardware collector, run this
preflight entry point or measure either computer's GPU.

It refuses an arbitrary interpreter or the Phase 1 `.venv` before querying
hardware. Use a fresh virtual environment named **`brd4-boltz221-cu128`**.
The script has no network, installation, MSA, checkpoint or inference code.
Importing it performs no diagnostics.

The commands below are **for approval and later execution on the Legion**.
They have not been executed here. Framework installation downloads packages
and can be large; it does not download Boltz model weights.

## Windows / PowerShell: GPU-only probe

After the Legion is powered and appropriately connected to its charger:
copy `boltz_gpu_preflight.py` into a folder and open PowerShell there.
Use 64-bit Python 3.12; if `py -3.12` is unavailable, install that Python
version separately before continuing.

```powershell
py -3.12 --version
if ($LASTEXITCODE -ne 0) { throw "Python 3.12 is required." }
$Brd4GpuEnv = Join-Path $env:USERPROFILE ".venvs\brd4-boltz221-cu128"
if (Test-Path $Brd4GpuEnv) {
    throw "Environment already exists. Stop and review it; do not overwrite Phase 1 or an earlier GPU environment."
}
py -3.12 -m venv "$Brd4GpuEnv"
if ($LASTEXITCODE -ne 0) { throw "Virtual environment creation failed." }
$Brd4GpuPython = Join-Path $Brd4GpuEnv "Scripts\python.exe"
& $Brd4GpuPython -m pip install "torch==2.7.1+cu128" --index-url https://download.pytorch.org/whl/cu128
if ($LASTEXITCODE -ne 0) { throw "Torch installation failed; stop." }
& $Brd4GpuPython -m pip install "numpy==1.26.4"
if ($LASTEXITCODE -ne 0) { throw "NumPy installation failed; stop." }
& $Brd4GpuPython -m pip check
if ($LASTEXITCODE -ne 0) { throw "Dependency check failed; stop." }
& $Brd4GpuPython .\boltz_gpu_preflight.py --output .\legion-windows-preflight.json
if ($LASTEXITCODE -ne 0) { throw "Preflight failed or is incomplete. Preserve its JSON and stop." }
& $Brd4GpuPython -m pip freeze | Out-File -Encoding utf8 .\legion-windows-packages.txt
```

Run the commands in order and stop if any guard throws. A report return code of 2 means failed or
incomplete preflight, not permission to proceed.
Do not use `--require-boltz` for this Windows GPU-only stage.

A native Windows CUDA pass establishes a useful hardware/runtime check.
The existing **Boltz pilot runtime gate still requires Linux**. Passing
this probe does not silently change that protocol to native Windows.

## Optional WSL2 / Linux route

Do not assume WSL is installed on the Legion. Check there, not on the
ThinkPad:

```powershell
wsl --list --verbose
```

If WSL/Ubuntu is absent, the following is a **separate approved OS setup
step**, run in administrator PowerShell, followed by any requested reboot
and initial Ubuntu account setup:

```powershell
wsl --install -d Ubuntu-24.04
```

Confirm the distribution is WSL **version 2** before proceeding. If it is
version 1 or GPU access is unavailable, stop and review the WSL setup.
Per [NVIDIA's WSL guide](https://docs.nvidia.com/cuda/wsl-user-guide/),
the Windows NVIDIA driver provides GPU access: do not install a Linux
display driver inside WSL. The prebuilt Torch route does not require a
separate CUDA toolkit or a source build.

Copy the standalone script to a directory accessible from Ubuntu, and
open a Linux shell in that directory. On Ubuntu 24.04, install venv support
only if missing:

```bash
# Optional prerequisite setup, after approval:
sudo apt-get update
sudo apt-get install -y python3.12-venv
```

Then create the separate Linux environment. These commands also work on a
native Linux host, but the script intentionally expects the Legion GPU:

```bash
set -euo pipefail
BRD4_GPU_ENV="$HOME/.venvs/brd4-boltz221-cu128"
test ! -e "$BRD4_GPU_ENV"
python3.12 -m venv "$BRD4_GPU_ENV"
BRD4_GPU_PY="$BRD4_GPU_ENV/bin/python"
"$BRD4_GPU_PY" -m pip install "torch==2.7.1+cu128" --index-url https://download.pytorch.org/whl/cu128
"$BRD4_GPU_PY" -m pip install "numpy==1.26.4"
"$BRD4_GPU_PY" -m pip check
"$BRD4_GPU_PY" ./boltz_gpu_preflight.py --output ./legion-wsl-preflight.json
"$BRD4_GPU_PY" -m pip freeze > ./legion-wsl-packages.txt
```

If the GPU probe passes and installing the Boltz **package** is approved,
the optional Linux package-metadata stage is:

```bash
"$BRD4_GPU_PY" -m pip install \
  "boltz==2.2.1" "rdkit==2026.3.6" "numpy==1.26.4" "torch==2.7.1+cu128" \
  --extra-index-url https://download.pytorch.org/whl/cu128
"$BRD4_GPU_PY" -m pip check
"$BRD4_GPU_PY" ./boltz_gpu_preflight.py \
  --require-boltz --output ./legion-wsl-boltz-metadata-preflight.json
"$BRD4_GPU_PY" -m pip freeze > ./legion-wsl-boltz-packages.txt
```

This checks package pins and the Linux requirement, not a Boltz prediction.
Do not invoke `boltz predict` or `capture-runtime` to obtain weights.
The latter requires already verified assets; upstream prediction may
automatically download them even for a structure-only request.

## Pass/fail criteria and what to send back

The JSON report records driver information, PyTorch CUDA build and native
architectures, actual GPU compute capability, device-wide memory, this
process's allocated/reserved memory, and the tiny probe's allocation peaks.
It does not collect process lists, credentials or general environment variables.
Existing report files are refused; use a new filename for a later attempt.

| Check | Passing condition |
|---|---|
| Environment | Dedicated named venv, Python 3.12; Phase 1 environment untouched |
| GPU | One visible RTX 5070 Ti **Laptop** GPU, compute capability `[12, 0]` |
| Driver/runtime | Driver baseline met, CUDA available, Torch `2.7.1+cu128`, CUDA build `12.8`, native `sm_120` listed |
| Model arithmetic prerequisite | BF16 supported; both tiny FP32 and BF16 operations correct |
| VRAM class | At least 11 GiB reported total, allowing driver reservations on a nominal 12 GB card |
| Starting headroom policy | At least 9 GiB free, using the lower of the snapshots before/after the tiny probe, both after CUDA initialization |
| Optional Boltz metadata stage | Linux, Boltz 2.2.1, RDKit 2026.3.6, Biopython 1.84; pip dependency check clean |

The **9 GiB threshold is an explicit conservative policy**, not a measured
Boltz requirement. A headroom failure with otherwise passing CUDA checks
does not prove the model cannot fit; review competing GPU usage and the
policy before any model run. Do not silently lower the threshold.
The script will not attempt even its tiny operations with less than 1 GiB
free after initialization.

Send the JSON report and the `pip check` result. The package listing is
useful if there is a dependency mismatch. A preflight pass does not
guarantee that every Boltz operation or the full-MSA input fits in memory.
The script's allocation peaks describe only this tiny probe.

## Troubleshooting

| Result | Next step |
|---|---|
| Wrong interpreter / environment | Use the explicit dedicated Python path above; do not install into `.venv`. |
| `nvidia-smi` absent or wrong GPU | Verify this is the powered-on Legion and its NVIDIA driver is installed. Stop; do not reinterpret ThinkPad results. |
| CPU Torch, CUDA ≤12.6, or no `sm_120` | Recheck the official `cu128` wheel and exact pin in the isolated environment. `--no_kernels` cannot fix an incompatible Torch build. |
| CUDA unavailable in WSL | Check WSL2 and Windows driver GPU support. Do not install a Linux display driver into WSL. |
| Tiny operation error / “no kernel image” | Preserve the JSON/error. Review driver, wheel and architecture before proceeding; no model trial or automatic retry. |
| Low free VRAM | Close unnecessary GPU applications and recheck later with a new report name. Host RAM does not enlarge device VRAM. |
| Power telemetry is N/A | This can occur on laptops/WSL; it is recorded but is not a CUDA-compatibility failure. |
| Package installation or `pip check` fails | Stop and share the dependency error; do not silently upgrade Torch, NumPy or RDKit. |
| Windows probe passes but Linux requirement fails | Hardware/runtime probe and the approved Linux Boltz environment are separate gates. Review WSL setup or a remote option. |

## Monitoring and explicit OOM stop for a later approved model run

Keep the smallest planned structure-only settings unchanged: one seed
`20261019`, one diffusion/parallel sample, 3 recycles, 200 sampling steps,
step scale 1.5, `--no_kernels`, no affinity, templates, constraints or
steering, and the pinned full MSA. No model run is authorised by this document.

For a later approved trial, record idle device memory, then sample global
VRAM once per second in a separate terminal:

```powershell
nvidia-smi --query-gpu=timestamp,memory.total,memory.used,memory.free,utilization.gpu --format=csv --loop-ms=1000 |
    Tee-Object -FilePath .\vram.csv
```

```bash
nvidia-smi --query-gpu=timestamp,memory.total,memory.used,memory.free,utilization.gpu \
  --format=csv --loop-ms=1000 > vram.csv
```

Stop that monitor with Ctrl+C after the trial. Save model stdout/stderr,
exit code, start/end times and the monitor log in that experiment's
directory. These samples include other GPU users and may miss brief peaks;
do not call them exact PyTorch allocator peaks.

**Stop on the first CUDA out-of-memory, driver/device error, nonzero model
exit, or the existing 30-minute ceiling. Preserve logs and do not retry
automatically.** Review sustained near-exhaustion/paging rather than
pressing on. Do not reduce MSA depth, sampling, precision or sequence
coverage in place. Any reduced-depth condition needs a reviewed validator/
config change, a distinct experiment ID and explicit approval. A 24 GB
remote run also needs approval and its own directory; none has been booked.
