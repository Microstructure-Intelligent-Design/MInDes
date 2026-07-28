# MInDes
\- **M**icrostructure **In**telligent **Des**ign Core

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free 
Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.

Cooperation Contact:

Dr. Qi Huang  - qihuang0908@163.com (e-mail) - hq5088028 (WeChat)

## Publications
**If you use **MInDes** solver in your research, please cite the following publications:**
1. **Huang, Qi**, ..., Yong Du*, et al. Phase-field simulation for voltage profile of LixSn nanoparticle during lithiation/delithiation, *Computational Materials Science*, 220 (2023): 112047.

2. **Huang, Qi***, et al. Multiphase transformation and mechanical analysis of polycrystalline CuxLiySn nanoparticle during lithiation via phase diagram-guided phase-field approach, *Electrochimica Acta*, 495 (2024): 144471.

3. Zhixuan Zhang, ..., **Huang, Qi***, et al. Multiscale phase-field simulation framework for spinodal decomposition behavior in composite carbide ceramics, *Journal of Materials Science & Technology*, 246 (2025): 58-75.

4. Xi Liu, **Huang, Qi***, et al. Phase-field simulation of grain coarsening and Mo-induced inhibition in NbC-Ni cermets during liquid phase sintering, *Journal of the American Ceramic Society*, 108 (2025): e70075.

5. Yiqi Guan, ..., **Qi Huang***, Weibin Zhang*, Yong Du*. Fracture mode and toughening mechanism induced by microstructure in binderless WC cemented carbides: A phase-field simulation integrating energy dissipation analysis, *Acta Materialia*, 2026, 305: 121834

## Building and Running MInDes

MInDes is officially supported on Windows x86/x64 and Linux x86_64. Release builds require a configured official License V2 public key in `MInDes/MInDes/src/modules/base/license_public_key.h`.

### Windows

Requirements:

- Windows 10 or later.
- Visual Studio 2022 with the **Desktop development with C++** workload and the MSVC v143 toolset.
- The FFTW headers, import libraries, and DLLs included under `MInDes/MInDes/lib/x64`.

To build with Visual Studio, open `MInDes/MInDes.sln`, select `x64` and either `Debug` or `Release`, and then build the solution.

To build from a Visual Studio Developer PowerShell opened at the repository root:

```powershell
msbuild MInDes/MInDes.sln /m /p:Configuration=Release /p:Platform=x64
```

The executable is written to `MInDes/x64/Release/MInDes.exe` for a Release build or `MInDes/x64/Debug/MInDes.exe` for a Debug build. Copy the FFTW runtime DLLs beside the executable before running it:

```powershell
Copy-Item MInDes/MInDes/lib/x64/*.dll MInDes/x64/Release/
```

For a Debug build, replace `Release` with `Debug` in both commands.

### Linux

The Linux build script supports x86_64 systems and uses GCC by default. On Ubuntu or Debian, install the required compiler, OpenSSL, and package-discovery tools:

```bash
sudo apt update
sudo apt install build-essential pkg-config libssl-dev
```

TPM 2.0 support is optional at build time. Install TPM2-TSS development files to enable it when available:

```bash
sudo apt install libtss2-dev
```

The bundled FFTW header and static library under `MInDes/MInDes/lib/linux` are used automatically.

Run the script from any directory. A Release build is the default:

```bash
bash MInDes/MInDes/scripts/compile_exec.sh
```

Common build variants are:

```bash
BUILD_TYPE=Debug bash MInDes/MInDes/scripts/compile_exec.sh
TPM_MODE=required bash MInDes/MInDes/scripts/compile_exec.sh
TPM_MODE=off bash MInDes/MInDes/scripts/compile_exec.sh
NPROC=8 CXX=clang++ CC=clang bash MInDes/MInDes/scripts/compile_exec.sh
```

`TPM_MODE=auto` is the default. It enables TPM support when `pkg-config` finds `tss2-fapi`; otherwise the script prints a warning and builds an SMBIOS-only client. `TPM_MODE=required` fails instead of falling back at build time, while `TPM_MODE=off` explicitly disables TPM support.

The resulting executable is written to:

```text
MInDes/MInDes/.build/linux/Release/MInDes
MInDes/MInDes/.build/linux/Debug/MInDes
```

### License Activation and Use

Keep the executable in a writable portable directory. The activation state, activation request, and signed license are always resolved relative to the executable directory, not the current working directory.

Display the command-line help:

```powershell
MInDes.exe -H
```

```bash
./MInDes -H
```

Open the information and license-management menu:

```powershell
MInDes.exe -I
```

```bash
./MInDes -I
```

Use the activation menu to register the current operating-system user and generate `mindes.request`. Send that request to the developer, then place the returned officially signed `mindes.license` beside the executable. The same directory will contain:

```text
mindes.activate
mindes.request
mindes.license
```

The activation file alone never grants solver permission; both the local activation check and the official license-signature check must pass. When deploying or moving MInDes, copy the complete portable directory, including the executable, runtime libraries, and license files.

On Linux, finding TPM2-TSS during compilation does not guarantee that TPM enrollment will succeed at runtime. The current user must also have access to the TPM device and an initialized TPM2-TSS FAPI environment. If initial TPM enrollment fails, License V2 records the reason and attempts SMBIOS binding; an existing TPM-bound license never falls back to SMBIOS during runtime validation.
