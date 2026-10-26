# nmr-endgroup-mn

**JEOL 一维氢谱预处理、逐谱积分复核与聚合物端基法分子量估算。**

把可重复的数值处理交给 Python，把溶剂核对、峰归属、积分边界与重叠判断交给 Codex skill。每份样品都经过 `survey → inspect → quantify`，不将上一份样品的区间直接套用到下一份。

目前支持 **JEOL JDF 一维复数时域 FID**；内置分子量公式仅适用于 **PDLLA–PEG–PDLLA 线性大分子二醇**。不支持直接处理已变换频域 JDF、二维谱或扩链聚氨酯的分子量。

## 功能

- 读取溶剂元数据并统一别名；缺失、未知、混合溶剂需要进一步核对。
- 保留原始横坐标绘制全谱和参照候选图；确认参照身份、原始峰位和依据后才允许校准。
- FID 起始点处理、指数加窗、补零、FFT；按实际谱线选择相位锚点和无峰基线区。
- 输出候选峰、累计积分、左右边界独立移动的敏感性图；不自动认定峰归属。
- 标记已指定溶剂/水/杂质区与基线或积分区的冲突；不自动删除信号。
- 有符号实部积分、参数敏感性检查、公式计算、图表与 JSON/CSV 输出。
- 源文件与配置哈希绑定复核记录，避免误用另一份数据或旧参数的审查。

仓库只包含通用代码、示例配置与合成测试，不包含实验原始数据、个人文件路径或样品结果。

## 安装

需要 Python 3.11 或更新版本。

```bash
git clone https://github.com/klu6163-droid/nmr-endgroup-mn.git
cd nmr-endgroup-mn
python -m venv .venv
```

Windows PowerShell 激活：

```powershell
.\.venv\Scripts\Activate.ps1
python -m pip install -r scripts/requirements.txt
```

macOS/Linux 激活：

```bash
source .venv/bin/activate
python -m pip install -r scripts/requirements.txt
```

如需作为 Codex skill 使用，将此目录放到 `$CODEX_HOME/skills/nmr-endgroup-mn`（未设置 CODEX_HOME 时使用 `~/.codex/skills/nmr-endgroup-mn`）。然后可请求：

> 使用 $nmr-endgroup-mn 处理这份核磁数据，先核对溶剂和积分边界，再按适用结构估算分子量。

## 三阶段工作流

### 1. Survey：先看溶剂与原始谱线

创建本地 `sample.json`。第一阶段只需以下处理字段：

```json
{
  "input_jdf": "/absolute/path/to/sample.jdf",
  "fid_drop_points": "auto",
  "line_broadening_hz": 0.3,
  "zero_fill_factor": 2,
  "reference": null
}
```

```bash
python scripts/workflow.py survey --config sample.json --out survey
```

读取 `survey/solvent_survey.json`，查看 `survey/solvent_survey.png`。幅度图仅用于定位信号，不用于定量。溶剂元数据属于采集设置，仍需核对实际残余峰。

### 2. Inspect：处理并复核积分区

参考 [完整示例字段](references/sample_config.example.json)，为本次样品补充相位、基线、宽搜索窗、待审查积分窗和结构参数。**示例中的数值只是字段演示，不能作为新样品的默认处理参数。**

无可靠参照时保持 `reference: null`。确认参照后，使用原坐标中的峰群中心填写：

```json
{
  "mode": "reviewed_peak",
  "kind": "solvent",
  "solvent": "DMSO-d6",
  "observed_ppm": 2.46,
  "target_ppm": 2.50,
  "source_sha256": "COPY_FROM_THIS_SAMPLE_SURVEY",
  "basis": "写明本谱参照峰身份、峰群中心和重叠检查依据"
}
```

以上数字仅为示例，不代表待测样品。该对象填入 `sample.json` 的 `reference` 字段。

```bash
python scripts/workflow.py inspect --config sample.json --out preview
```

实际查看 `preview/interval_review.png`、候选峰和边界扫描。调整区间后重新运行。完整峰群覆盖、基线与邻峰重叠需要判读；稳定面积不能证明归属正确。

### 3. Quantify：有依据地计算

将 `preview/review_draft.json` 另存为 `review.json`，填写真实观察依据。支持条件性计算时设置 `decision` 为 `quantify_provisionally`。

```bash
python scripts/workflow.py quantify --config sample.json --review review.json --out results
```

复核记录是分析者/agent 的技术记录，不是用户审批流程，也不能机械填字代替看图。若峰群不可分离或公式不适用，应保留 `needs_review` 并说明原因。

输出包括 `report.md`、`spectrum_review.png/.pdf`、`spectrum.csv`、`processing_arrays.npz`、`result.json`、配置与复核记录。原始 JDF 仅读取。重复使用同一输出目录会更新对应结果；保留不同分析版本时使用不同目录。

## 分子量公式与边界

对每条链含一个 PEG 核、两个 PDLLA 臂及两个乳酸羟基末端的线性二醇：

```text
Mn = (6 × I_main / I_end + 2) × 72 + PEG_Mn
```

`I_main` 为非末端乳酸 CH 面积；`I_end` 假定完整涵盖两个末端乳酸 CH（2H）和两侧 PEG 酯相邻 CH2（4H）。`+2` 补回末端乳酸单元。PEG_Mn 和重复单元摩尔质量可配置。

端基氢数、峰归属、纯度和链模型必须与实际样品一致。环状/支化物、单端接枝、游离小分子、重叠峰或扩链聚氨酯不能直接套用该式。无 T1 或其他定量验证时，报告的是条件性估算。

FID 起始截断为经验估计，不是已确认的 JEOL 群延迟常数；其敏感性会被检查。参数敏感性范围不是统计置信区间。算法可给出诊断，不能自动证明化学身份或绝对分子量。

详细字段、判读规则和来源见 [SKILL.md](SKILL.md) 与 [工作流说明](references/workflow.md)。

## 测试

```bash
python -m unittest discover -s scripts -p "test_*.py" -v
```

15 项测试覆盖已知相位回收、FFT 轴方向、有符号积分、基线后面积比、峰位漂移、复核失效、溶剂别名、参照冲突与重叠标记。测试不依赖真实样品。

## 依赖与来源

第三方依赖通过 pip 安装，仓库不复制第三方源码。

- [nmrglue JEOL reader](https://github.com/jjhelmus/nmrglue/blob/master/nmrglue/fileio/jeol.py)
- [Sigma 溶剂与杂质位移表](https://b2b.sigmaaldrich.com/BR/en/technical-documents/technical-article/analytical-chemistry/nuclear-magnetic-resonance/1h-nmr-and-13c-nmr-chemical-shifts-of-impurities-chart)
- [JEOL qNMR 与弛豫条件](https://www.jeol.com/solutions/applications/details/qNMR-issue8.php)
