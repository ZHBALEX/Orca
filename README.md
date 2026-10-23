# Orca — 代码与研究报告

本仓库保存助手为 Orca 编写或修改的代码及对应报告，供远程查看。仿真数据、个人文档和论文全文保留本地。

## 从这里开始

- **[最新研究进度、创新点与新增分析](reports/2026-10-23-evidence-audit/README.md)**
- **[23 份去重历史报告目录](reports/history/README.md)**：包含历史草稿与片段，已标明不能代替最新结论。
- **[英文论文结果与讨论草稿](reports/2026-10-23-evidence-audit/MANUSCRIPT_DRAFT.md)**
- **[下一步研究任务与验收标准](reports/2026-10-23-evidence-audit/NEXT.md)**

本轮完成：平均/相位压力闭合复核、RV⁺ 同点释放检查、固定材料邻域的压力与受力检验，以及小/大邻域动量项比较。最重要的进展是将涡结构和局部表面压力建立空间联系，同时发现不能把全上表面 Lamb 主导的均值结论直接套到最靠近涡核的小邻域。

![局部压力和材料面片分析](reports/2026-10-23-evidence-audit/03_local_pressure_bridge.png)

## 代码

- `code/run__evidenceAudit.py`：读取保存结果，复核平均量、相位闭合及轨迹，生成图 1–2。
- `code/run__localPressureBridge.py`：读取本地 profile 数据，生成固定材料邻域及动量预算图 3–4。
- `code/vortex_atlas/run__vortexAtlas.py`：既有 v10 三幅涡结构图的维护版代码快照。它使用本地缓存；本轮未重新运行该图谱流程。

先在 runner 的 `main()` 设置本机数据路径，再运行。新分析依赖 NumPy 和 Matplotlib；图谱还需要 SciPy。输入文件未随仓库上传；报告中的输入哈希和方法说明用于核对来源，不表示任意机器可无数据复算。

图内科学内容以栅格层嵌入 PDF/SVG，文字和坐标保留可编辑元素；图像是报告的一部分。GitHub 可直接阅读 Markdown 与 PNG；HTML 报告可下载后打开。

项目规则见 [PROJECT.md](docs/agent-loop/PROJECT.md)。研究记录与英文草稿均非投稿定稿。
