# 配置和判读记录

## 配置

流程为 survey → inspect → quantify。survey 只需 input_jdf、fid_drop_points、line_broadening_hz、zero_fill_factor，可把 reference 设为 null；它读取溶剂并绘制未校准幅度谱，不做积分。幅度图仅定位，真正积分仍用校正后的吸收实部。

溶剂别名统一后保存原值、规范名及来源。暂支持 DMSO-d6、CDCl3、acetone-d6、CD3CN、CD3OD、C6D6 的常见残余峰候选；D2O 不给固定 HOD 参照。其余/混合溶剂保留原值并提示澄清，不能静默回退到 DMSO。参照表来自 [Sigma 溶剂杂质位移表](https://b2b.sigmaaldrich.com/BR/en/technical-documents/technical-article/analytical-chemistry/nuclear-magnetic-resonance/1h-nmr-and-13c-nmr-chemical-shifts-of-impurities-chart)，只是定位起点，环境条件和参照约定可能影响实际位置。

确认峰身份后设置 reference：`mode: reviewed_peak`、`kind: solvent`（或有依据的 `standard`）、`solvent`、`observed_ppm`（未平移坐标中的多重峰中心）、`target_ppm`、`basis`、`source_sha256`。默认允许平移不超过0.25 ppm；超出时需检查原因，不能只为通过检查而提高阈值。校准使用明确记录的位置，不重新选择窗口中的最高峰。旧版仅含 search_ppm 的自动参照配置将被拒绝。

若文件溶剂字段填错，可以设置 `solvent_override: {name: ..., basis: ...}`，依据应来自用户或实验记录，保留原值供追踪。没有可靠参照可保持 reference=null。所有相位锚点、基线区和积分区都应使用最终选定的坐标系。

`observed_contaminant_regions` 为实际观察到的水/杂质范围列表，各项含 label 和 bounds_ppm；不要仅按溶剂名自动生成固定水峰位置。脚本标记这些范围及参照溶剂峰附近与基线/积分区的交叠。处理时重新选无峰基线区，定量区有重叠则说明可否分离；不自动扣除整个区间。存在交叠但仍可作条件估算时，在 review 的 solvent_overlap_resolution 里记录具体依据。

以 `sample_config.example.json` 的字段结构创建一份样品配置，替换输入路径和数值。该例中的数字来自一份 DMSO-d6 样品，不能无检查迁移到其他文件。

`review_search_windows` 是宽搜索窗；`main_integral_ppm` 与 `end_integral_ppm` 是待复核积分窗，两者不要混淆。搜索窗应涵盖完整候选峰群、峰脚与两侧基线。参考峰、相位锚点及基线区同样是样品相关参数。

inspect 生成的峰候选依赖噪声和突出度阈值，可能漏掉低信噪宽峰；不能用候选峰列表替代原始谱线。累计积分应在峰群前后出现合理平缓区域，但基线偏移、其他峰重叠或截尾可能造成误导。结合相位后实部、虚部、基线以及未经平滑的谱线共同检查。

边界扫描在当前下边界和上边界各 ±0.05 ppm 上独立移动（11×11），对照有符号积分。它检验的是局部边界敏感性，不证明信号归属，也不覆盖所有可能漂移。根据峰宽、邻峰和实际漂移幅度扩大诊断区，必要时调整脚本扫描范围。不要跨入相邻信号后把变化简单称为测量误差。

## 复核 JSON

inspect 的 `review_draft.json` 中哈希已计算。查看图后另存 `review.json`，写明：

- `sample_identity_basis`：样品身份来自用户/实验记录，哪些只是推断。
- `solvent_and_reference_assessment`：原始溶剂字段、规范名、残余峰核对、是否采用参照及其依据、可疑水峰/重叠；不能将元数据视为实验组成证明。
- `phase_assessment`：峰脚、负尾及相位拟合残差；弱峰区是否同时改善。
- `baseline_assessment`：锚点是否无峰，峰两侧是否回到基线。
- `main_assignment_and_bounds`：主链归属、最终上下限、是否覆盖完整多重峰及排除其他信号的依据。
- `end_assignment_and_bounds`：端基峰群归属、6H 假设、最终边界依据。
- `overlap_assessment`：可见重叠、杂质或低信噪问题；无法排除的说明不确定性，不能写成已证明无杂质。
- `boundary_sensitivity_assessment`：独立左右移动、平移及宽度改变后的面积变化和选择理由。
- `formula_applicability`：所用链模型、PEG Mn、每条链的参考氢数、末端补偿及弛豫限制。

支持条件性计算时使用 `decision=quantify_provisionally`；无法可靠选区时保持 `needs_review` 并报告原因。任何影响处理的配置变化都会使旧哈希失效。更换输入文件也需要重新判读；更改输出目录不影响哈希。

## PDLLA–PEG–PDLLA 公式

令 A 为非末端乳酸 CH 面积，E 为两个末端乳酸 CH（2H）与两侧酯相邻 PEG CH2（4H）的合计面积：

`Mn = (6*A/E + 2)*72 + PEG_Mn`

必须涵盖这组完整的 6H，不能只根据窗口内出现若干峰就确认氢数。不适用于已经扩链的聚氨酯。72 与 PEG400 是用户所给图片的计算约定；若改用其他摩尔质量，报告实际数值。

原图片存在两处不一致：端基范围写作 4.23–4.46 及 4.23–4.36 ppm；17.08/1.00 按公式得到 7922.56 g/mol，而文中写 7700。它们是来源问题，不是校正本谱的目标。

## 实现与来源

engine 使用 nmrglue JEOL reader，读取符号约定与 NumPy FFT 轴方向已通过已知频率信号验证。相位可用复数线形锚点加线性回归或候选区熵/负峰惩罚；后者曾在本样品残留色散尾，因此必须看图。基线采用无峰区中位数插值；初始点截断是经验估计，要做敏感性检查。

- [nmrglue JEOL reader](https://github.com/jjhelmus/nmrglue/blob/master/nmrglue/fileio/jeol.py)
- [DMSO 残余峰参照](https://chem.vt.edu/analytical-services/nmr.html)
- [JEOL 定量与弛豫条件](https://www.jeol.com/solutions/applications/details/qNMR-issue8.php)
